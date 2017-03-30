#!/usr/bin/env python

import argparse
import string
import sys
import os
import re


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fasta_file",
                        help="FASTA file that will be used to create reference data for TreeSAPP", required=True)
    # parser.add_argument("-t", "--table",
    #                     help="Tab-separated columns listing sequences to be used in reference data")
    parser.add_argument("-c", "--code_name",
                        help="Unique name to be used by TreeSAPP internally.\n"
                             "Refer to the first column of 'cog_list.txt' under the '#functional cogs' section)",
                        required=True)
    parser.add_argument("-m", "--min_length",
                        help="The minimum length of a protein to be used in building reference data. [ DEFAULT = 0 ]",
                        required=False,
                        default=0)

    args = parser.parse_args()
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
    reg_aa_ambiguity = re.compile(r'[bjxzBJXZ]')
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
        sys.stderr.write("ERROR: duplicate header names were detected in " + args.input + "!\n")
        sys.stderr.flush()
        sys.exit(2)

    if count_total == 0:
        sys.exit('ERROR: Your input file appears to be corrupted. No sequences were found!\n')

    # Exit the program if all sequences are composed only of X or N
    elif count_xn == count_total:
        sys.exit('ERROR: Your sequence(s) contain only X or N!\n')

    # Exit the program if less than half of the characters are nucleotides
    elif float(seq_count / float(count_total)) < 0.5:
        sys.exit('ERROR: Your sequence(s) most likely contain no AA!\n')

    if len(substitutions) > 0:
        sys.stderr.write("WARNING: " + str(len(substitutions)) + " ambiguous character substitutions were made!\n")
        sys.stderr.flush()

    # Close the files
    fasta.close()

    return formatted_fasta_dict


def create_new_fasta(fasta_dict, out_fasta, dictionary):
    """
    Writes a new FASTA file with the headers
    :param fasta_dict:
    :param out_fasta:
    :param dictionary:
    :return:
    """
    out_fasta_handle = open(out_fasta, "w")

    for header in fasta_dict.keys():
        header_type = get_header_format(header)

        if header_type == "ncbi":
            header_match = re.match(">gi\|(\d+)\|(\w+)\|(\S+(\.\d+)*)\|(.*)$", header)
            id = header_match.group(3)
        elif header_type == "fungene":
            header_match = re.match("^>(\d+)  coded_by=(.+),organism=(.+),definition=(.+)$", header)
            id = header_match.group(1)
        elif header_type == "mltree":
            header_match = re.match(">(r_\d+)$", header)
            id = header_match.group(1)
        else:
            print "ERROR: Incorrect regex matching of header!"
            sys.exit(3)

        if id in dictionary:
            out_fasta_handle.write(">%s\n" % (dictionary[id]))
        else:
            print id, "not in dictionary"
            out_fasta_handle.write(fasta_dict[header])

    out_fasta_handle.close()
    return


def get_header_format(header):
    """
    Used to decipher which formatting style was used: NCBI, FunGenes, or other
    :param header: A sequences header from a FASTA file
    :return:
    """
    ncbi_re = re.compile(">gi\|(\d+)\|(\w+)\|(\S+(\.\d+)*)\|(.*)$")
    fungene_re = re.compile("^>(\d+)  coded_by=(.+),organism=(.+),definition=(.+)$")
    if ncbi_re.match(header):
        return "ncbi"
    if fungene_re.search(header):
        return "fungene"
    else:
        return None


def get_sequence_info(code_name, fasta_dict):
    """
    This function is used to find the accession ID and description of each sequence from the FASTA file
    :param code_name: code_name from the command-line parameters
    :param fasta_dict: a dictionary with headers as keys and sequences as values (returned by format_read_fasta)
    :return:
    """

    fasta_repl_dict = {}
    fasta_mltree_repl_dict = {}
    tree_repl_dict = {}
    mltree_dict = {}
    tree_name_dict = {}
    mltree_id_accumulator = 1

    print "******************** Generating information for formatting purposes ********************"

    for header in fasta_dict.keys():
        mltree_id = str(mltree_id_accumulator)
        mltree_id_accumulator += 1
        header_format = get_header_format(header)
        if header_format == "fungene":
            accession, info = header[1:].split("  ")
            fungene_info = re.match("^coded_by=(.+),organism=(.+),definition=(.+)$", info)
            if fungene_info:
                # coded = fungene_info.group(1)
                # organism = fungene_info.group(2)
                description = fungene_info.group(2)
                # definition = fungene_info.group(3)
            else:
                print "Fail."
                sys.exit()
            short_id = code_name + '_' + mltree_id
            tree_id = mltree_id
        else:
            sys.exit()
        print "accession =", accession
        print "short_id =", short_id
        print "mltree_id =", mltree_id
        print "description =", description
        # sys.exit()

        fasta_repl_dict[accession] = short_id
        fasta_mltree_repl_dict[short_id] = mltree_id
        tree_repl_dict[short_id] = description
        mltree_dict[short_id] = tree_id
        tree_name_dict[tree_id] = description

    return tree_name_dict, mltree_dict, fasta_repl_dict, fasta_mltree_repl_dict


def write_tax_ids(code_name, tree_name_dict):
    tree_taxa_list = "tax_ids_%s.txt" % code_name

    tree_tax_list_handle = open(tree_taxa_list, "w")

    for mltree_id_key in sorted(tree_name_dict.keys(), key=int):
        tree_tax_list_handle.write("%s\t%s\n" % (mltree_id_key, tree_name_dict[mltree_id_key]))
    tree_tax_list_handle.close()
    return tree_taxa_list


def main():
    args = get_arguments()

    input_fasta = args.fasta_file
    fasta_dict = format_read_fasta(input_fasta, args)
    code_name = args.code_name

    tree_name_dict, mltree_dict, fasta_repl_dict, fasta_mltree_repl_dict = get_sequence_info(code_name, fasta_dict)

    tree_taxa_list = write_tax_ids(code_name, tree_name_dict)

    fasta_replace_names = code_name + "_fasta_replace.names"
    fasta_mltree_names = code_name + "_fasta_mltree.names"
    tree_replace_names = code_name + "_tree_replace.names"
    
    print "******************** %s generated ********************\n" % tree_taxa_list
    
    tree_names_list = "%s_tree_replace.names" % code_name
    
    tree_names_list_handle = open(tree_names_list, "w")
    
    for short_id_key in sorted(mltree_dict.keys()):
        tree_names_list_handle.write("%s:\t%s:\n" % (short_id_key, mltree_dict[short_id_key]))
    
    tree_names_list_handle.close()
    
    fasta_replaced = code_name + ".fc.repl.fasta"
    
    create_new_fasta(fasta_dict, fasta_replaced, fasta_repl_dict)
    
    print "******************** %s generated ********************\n" % tree_names_list

    print "******************** FASTA file, %s generated ********************\n" % fasta_replaced
    
    print "******************** Aligning the sequences using MUSCLE ********************\n"
    
    fasta_replaced_align = code_name + ".fc.repl.aligned.fasta"

    muscle_align_command = "muscle -in %s -out %s" % (fasta_replaced, fasta_replaced_align)
    
    print muscle_align_command, "\n"
    os.system(muscle_align_command)
    fasta_replaced_align_dict = format_read_fasta(code_name + ".fc.repl.aligned.fasta", args)
    
    fasta_mltree = code_name + ".fa"
    
    create_new_fasta(fasta_replaced_align_dict, fasta_mltree, fasta_mltree_repl_dict)
    
    print "******************** FASTA file, %s generated ********************\n" % fasta_mltree
    
    makeblastdb_command = "makeblastdb -in %s -dbtype prot -input_type fasta -out %s" % (fasta_mltree, fasta_mltree)
    os.system(makeblastdb_command)
    
    print "******************** BLAST DB for %s generated ********************\n" % code_name
    
    hmm_build_command = "hmmbuild -s %s.hmm %s" % (code_name, fasta_mltree)
    os.system(hmm_build_command)
    
    print "******************** HMM file for %s generated ********************\n" % code_name
    
    phylip_command = "java -cp readseq.jar run -a -f=12 %s" % fasta_replaced_align
    os.system(phylip_command)
    
    phylip_file = code_name + ".phy"
    os.system('mv %s.phylip %s' % (fasta_replaced_align, phylip_file))
    
    raxml_out = "%s_phy_files" % code_name
    os.system("mkdir %s" % raxml_out)
    
    raxml_command = "raxmlHPC-PTHREADS -f a -x 12345 -# 100 -m PROTGAMMAWAG -s %s -n %s -w %s -T 4" % (phylip_file, code_name, raxml_out)
    os.system(raxml_command)

    tree_to_swap = "%s/RAxML_bestTree.%s" % (raxml_out, code_name)
    final_mltree = "%s_tree.txt" % code_name
    
    swapTree_command = "swapTreeNames.pl -t %s -l %s -o %s" % (tree_to_swap, tree_names_list, final_mltree)
    os.system(swapTree_command)

    final_output_folder = "MLTreeMap_files_%s" % (code_name)
    os.system("mkdir %s" % (final_output_folder))
    
    os.system("mv %s.fa %s.fa.p* %s" % (code_name, code_name, final_output_folder))
    os.system("mv %s.hmm %s %s %s" % (code_name, tree_taxa_list, final_mltree, final_output_folder))

main()
