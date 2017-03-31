#!/usr/bin/env python

import argparse
import sys
import os
import re
from mltreemap import os_type, is_exe, which


class ReferenceSequence():
    def __init__(self):
        self.accession = ""
        self.description = ""
        self.definition = ""
        self.short_id = ""
        self.sequence = ""
        self.locus = ""

    def get_info(self):
        sys.stdout.write("accession = " + self.accession + "\t")
        sys.stdout.write("locus = " + self.locus + "\t")
        sys.stdout.write("description = " + self.description + "\t")
        sys.stdout.write("mltree_id = " + self.short_id + "\n")
        sys.stdout.flush()


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fasta_file",
                        help="FASTA file that will be used to create reference data for TreeSAPP", required=True)
    parser.add_argument("-c", "--code_name",
                        help="Unique name to be used by TreeSAPP internally.\n"
                             "Refer to the first column of 'cog_list.txt' under the '#functional cogs' section)",
                        required=True)
    parser.add_argument("-m", "--min_length",
                        help="The minimum length of a protein to be used in building reference data. [ DEFAULT = 0 ]",
                        required=False,
                        default=0)
    parser.add_argument("-T", "--num_threads",
                        help="The number of threads for RAxML to use [ DEFAULT = 4 ]",
                        required=False,
                        default=4)

    args = parser.parse_args()
    args.mltreemap = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + os.sep

    return args


def find_executables(args):
    """
    Finds the executables in a user's path to alleviate the requirement of a sub_binaries directory
    :param args: command-line arguments objects
    :return: exec_paths beings the absolute path to each executable
    """
    exec_paths = dict()
    dependencies = ["raxmlHPC", "makeblastdb", "muscle", "hmmbuild"]

    if os_type() == "linux":
        args.executables = args.mltreemap + "sub_binaries" + os.sep + "ubuntu"
    if os_type() == "mac":
        args.executables = args.mltreemap + "sub_binaries" + os.sep + "mac"
    elif os_type() == "win" or os_type() is None:
        sys.exit("ERROR: Unsupported OS")

    for dep in dependencies:
        if is_exe(args.executables + os.sep + dep):
            exec_paths[dep] = str(args.executables + os.sep + dep)
        # For rpkm and potentially other executables that are compiled ad hoc
        elif is_exe(args.mltreemap + "sub_binaries" + os.sep + dep):
            exec_paths[dep] = str(args.mltreemap + "sub_binaries" + os.sep + dep)
        elif which(dep):
            exec_paths[dep] = which(dep)
        else:
            sys.stderr.write("Could not find a valid executable for " + dep + ". ")
            sys.exit("Bailing out.")

    args.executables = exec_paths
    return args


def format_read_fasta(fasta_file, args):
    """
    Splits the input file into multiple files, each containing a maximum number of sequences as specified by the user.
    Ensures each sequence and sequence name is valid.
    :param args: Command-line argument object from get_arguments
    :return A list of the files produced from the input file.
    """
    fasta = open(fasta_file, 'r')
    formatted_fasta_dict = dict()
    header = ""
    sequence = ""
    convert_upper = lambda pat: pat.group(1).upper()
    reg_aa_ambiguity = re.compile(r'[bjxzBJXZ-]')
    reg_amino = re.compile(r'[acdefghiklmnpqrstuvwyACDEFGHIKLMNPQRSTUVWY*]')
    substitutions = ""
    count_total = 0
    seq_count = 0
    count_xn = 0
    num_headers = 0
    for line in fasta:
        # If the line is a header...
        if line[0] == '>':
            if len(header) > 0:
                if len(sequence) > args.min_length:
                    formatted_fasta_dict[header] = sequence
                else:
                    num_headers -= 1

            sequence = ""

            header = line.strip()
            num_headers += 1

        # Else, if the line is a sequence...
        else:
            characters = line.strip()
            if len(characters) == 0:
                continue
            # Remove all non-characters from the sequence
            re.sub(r'[^a-zA-Z]', '', characters)

            # Count the number of {atcg} and {xn} in all the sequences
            count_total += len(characters)

            characters = re.sub(r'([a-z])', convert_upper, characters)
            aminos = len(reg_amino.findall(characters))
            ambiguity = len(reg_aa_ambiguity.findall(characters))
            if (aminos + ambiguity) != len(characters):
                sys.stderr.write("ERROR: " + header.strip() + " contains unknown characters: ")
                unknown = ""
                for c in characters:
                    if c not in "abcdefghiklmnpqrstvwxyzABCDEFGHIKLMNPQRSTVWXYZ":
                        unknown += c
                sys.stderr.write(unknown + "\n")
                sys.exit()
            else:
                seq_count += aminos
                count_xn += ambiguity
            sequence += characters
    formatted_fasta_dict[header] = sequence

    # Check for duplicate headers
    if len(formatted_fasta_dict.keys()) != num_headers:
        sys.stderr.write("ERROR: duplicate header names were detected in " + fasta_file + "!\n")
        sys.stderr.flush()
        sys.exit(2)

    if count_total == 0:
        sys.exit('ERROR: Your input file appears to be corrupted. No sequences were found!\n')

    # Exit the program if all sequences are composed only of X or N
    elif count_xn == count_total:
        sys.exit('ERROR: Your sequence(s) contain only X or N!\n')

    if len(substitutions) > 0:
        sys.stderr.write("WARNING: " + str(len(substitutions)) + " ambiguous character substitutions were made!\n")
        sys.stderr.flush()

    # Close the files
    fasta.close()

    return formatted_fasta_dict


def create_new_fasta(code_name, fasta_dict, out_fasta, dictionary, dashes=True):
    """
    Writes a new FASTA file with the headers
    :param code_name:
    :param fasta_dict:
    :param out_fasta:
    :param dictionary:
    :param dashes:
    :return:
    """
    out_fasta_handle = open(out_fasta, "w")

    for header in fasta_dict.keys():
        header_type = get_header_format(header, code_name)
        short_id = ""

        if header_type == "ncbi":
            header_match = re.match(">gi\|(\d+)\|(\w+)\|(\S+(\.\d+)*)\|(.*)$", header)
            accession = header_match.group(3)
        elif header_type == "fungene":
            header_match = re.match("^>(\d+)\s+coded_by=(.+),organism=(.+),definition=(.+)$", header)
            accession = header_match.group(1)
            coded_by = header_match.group(2)
            for mltree_id in dictionary:
                if dictionary[mltree_id].accession == accession and dictionary[mltree_id].locus == coded_by:
                    short_id = dictionary[mltree_id].short_id
                    break
        elif header_type == "mltree":
            header_match = re.match("^>(\d+)_" + re.escape(code_name), header)
            mltree_id = header_match.group(1)
            accession = mltree_id
            short_id = dictionary[mltree_id].short_id
        else:
            print "ERROR: Incorrect regex matching of header!"
            sys.exit(3)

        if short_id:
            out_fasta_handle.write(">%s\n" % short_id)
            if dashes:
                out_fasta_handle.write(fasta_dict[header] + "\n")
            else:
                purified = re.sub('-', '', fasta_dict[header])
                out_fasta_handle.write(purified + "\n")
        else:
            print accession, "not in dictionary"

    out_fasta_handle.close()
    return


def remove_dashes_from_msa(fasta_in, fasta_out):
    dashed_fasta = open(fasta_in, 'r')
    fasta = open(fasta_out, 'w')
    sequence = ""

    line = dashed_fasta.readline()
    while line:
        if line[0] == '>':
            if sequence:
                fasta.write(sequence + "\n")
                sequence = ""
            fasta.write(line)
        else:
            sequence += re.sub('-', '', line.strip())
        line = dashed_fasta.readline()
    fasta.write(sequence + "\n")
    dashed_fasta.close()
    fasta.close()
    return


def get_header_format(header, code_name):
    """
    Used to decipher which formatting style was used: NCBI, FunGenes, or other
    :param header: A sequences header from a FASTA file
    :param code_name:
    :return:
    """
    ncbi_re = re.compile(">gi\|(\d+)\|(\w+)\|(\S+(\.\d+)*)\|(.*)$")
    fungene_re = re.compile("^>(\d+)\s+coded_by=(.+),organism=(.+),definition=(.+)$")
    mltree_re = re.compile("^>(\d+)_" + re.escape(code_name))
    if ncbi_re.match(header):
        return "ncbi"
    if fungene_re.search(header):
        return "fungene"
    if mltree_re.match(header):
        return "mltree"
    else:
        return None


def get_sequence_info(code_name, fasta_dict):
    """
    This function is used to find the accession ID and description of each sequence from the FASTA file
    :param code_name: code_name from the command-line parameters
    :param fasta_dict: a dictionary with headers as keys and sequences as values (returned by format_read_fasta)
    :return:
    """

    fasta_mltree_repl_dict = {}
    mltree_id_accumulator = 1

    print "******************** Generating information for formatting purposes ********************"

    for header in fasta_dict.keys():
        mltree_id = str(mltree_id_accumulator)
        mltree_id_accumulator += 1
        ref_seq = ReferenceSequence()
        ref_seq.sequence = fasta_dict[header]
        header_format = get_header_format(header, code_name)
        if header_format == "fungene":
            accession, info = header[1:].split("  ")
            ref_seq.accession = accession
            fungene_info = re.match("^coded_by=(.+),organism=(.+),definition=(.+)$", info)
            if fungene_info:
                ref_seq.locus = fungene_info.group(1)
                # organism = fungene_info.group(2)
                # definition = fungene_info.group(3)
                ref_seq.description = fungene_info.group(2)
            else:
                print "Fail."
                sys.exit()
            short_id = mltree_id + '_' + code_name
            ref_seq.short_id = short_id
        else:
            print "Unable to handle header:", header
            sys.exit()

        fasta_mltree_repl_dict[mltree_id] = ref_seq

    return fasta_mltree_repl_dict


def write_tax_ids(code_name, fasta_mltree_repl_dict):
    tree_taxa_list = "tax_ids_%s.txt" % code_name

    tree_tax_list_handle = open(tree_taxa_list, "w")

    for mltree_id_key in sorted(fasta_mltree_repl_dict.keys(), key=int):
        tree_tax_list_handle.write("%s\t%s\n" % (mltree_id_key, fasta_mltree_repl_dict[mltree_id_key].description))
    tree_tax_list_handle.close()
    return tree_taxa_list


def swap_tree_names(tree_to_swap, final_mltree, code_name):
    original_tree = open(tree_to_swap, 'r')
    raxml_tree = open(final_mltree, 'w')

    tree = original_tree.readlines()
    original_tree.close()
    if len(tree) > 1:
        sys.stderr.write("ERROR: >1 line contained in RAxML tree " + tree_to_swap)

    new_tree = re.sub('_' + re.escape(code_name), '', str(tree[0]))
    raxml_tree.write(new_tree)

    raxml_tree.close()
    return


def main():
    args = get_arguments()
    args = find_executables(args)

    input_fasta = args.fasta_file
    fasta_dict = format_read_fasta(input_fasta, args)
    code_name = args.code_name

    fasta_mltree_repl_dict = get_sequence_info(code_name, fasta_dict)

    tree_taxa_list = write_tax_ids(code_name, fasta_mltree_repl_dict)
    
    print "******************** %s generated ********************\n" % tree_taxa_list

    fasta_replaced = code_name + ".fc.repl.fasta"
    
    create_new_fasta(code_name, fasta_dict, fasta_replaced, fasta_mltree_repl_dict)

    print "************************** FASTA file, %s generated *************************\n" % fasta_replaced
    
    print "******************** Aligning the sequences using MUSCLE ********************\n"
    
    fasta_replaced_align = code_name + ".fc.repl.aligned.fasta"

    muscle_align_command = "%s -in %s -out %s" %\
                           (args.executables["muscle"], fasta_replaced, fasta_replaced_align)
    
    print muscle_align_command, "\n"
    os.system(muscle_align_command)
    
    fasta_mltree = code_name + ".fa"

    remove_dashes_from_msa(fasta_replaced_align, fasta_mltree)
    
    print "******************** FASTA file, %s generated ********************\n" % fasta_mltree
    
    makeblastdb_command = "%s -in %s -dbtype prot -input_type fasta -out %s" %\
                          (args.executables["makeblastdb"], fasta_mltree, fasta_mltree)
    os.system(makeblastdb_command)
    
    print "******************** BLAST DB for %s generated ********************\n" % code_name
    
    hmm_build_command = "%s -s %s.hmm %s" %\
                        (args.executables["hmmbuild"], code_name, fasta_replaced_align)
    os.system(hmm_build_command)
    
    print "******************** HMM file for %s generated ********************\n" % code_name
    
    phylip_command = "java -cp %s/sub_binaries/readseq.jar run -a -f=12 %s" % (args.mltreemap, fasta_replaced_align)
    os.system(phylip_command)
    
    phylip_file = code_name + ".phy"
    os.system('mv %s.phylip %s' % (fasta_replaced_align, phylip_file))
    
    raxml_out = "%s_phy_files" % code_name

    if not os.path.exists(raxml_out):
        os.system("mkdir %s" % raxml_out)
    
    raxml_command = "%s -f a -p 12345 -x 12345 -# 100 -m PROTGAMMAWAG -s %s -n %s -w %s -T %s" %\
                    (args.executables["raxmlHPC"], phylip_file, code_name, args.mltreemap + raxml_out, args.num_threads)
    os.system(raxml_command)

    tree_to_swap = "%s/RAxML_bestTree.%s" % (raxml_out, code_name)
    final_mltree = "%s_tree.txt" % code_name

    # TODO: Update the data/tree_data/cog_list.txt file with the new marker gene
    swap_tree_names(tree_to_swap, final_mltree, code_name)

    final_output_folder = "MLTreeMap_files_%s" % code_name
    os.system("mkdir %s" % final_output_folder)

    os.system("mv %s.fa %s.fa.p* %s" % (code_name, code_name, final_output_folder))
    os.system("mv %s.hmm %s %s %s" % (code_name, tree_taxa_list, final_mltree, final_output_folder))

main()
