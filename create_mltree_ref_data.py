#!/usr/bin/env python

import argparse
import sys
import os
import re
import subprocess
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
                        help="FASTA file that will be used to create reference data for TreeSAPP",
                        required=True)
    parser.add_argument("-u", "--uc",
                        help="The USEARCH cluster format file produced from clustering reference sequences",
                        required=False,
                        default=None)
    parser.add_argument("-c", "--code_name",
                        help="Unique name to be used by TreeSAPP internally. NOTE: Must be >=5 characters.\n"
                             "Refer to the first column of 'cog_list.txt' under the '#functional cogs' section)",
                        required=True)
    parser.add_argument("-m", "--min_length",
                        help="The minimum length of a protein to be used in building reference data. [ DEFAULT = 0 ]",
                        required=False,
                        default=0,
                        type=int)
    parser.add_argument("-b", "--bootstraps",
                        help="The number of bootstrap replicates RAxML should perform [ DEFAULT = 100 ]",
                        required=False,
                        default=100,
                        type=int)
    parser.add_argument("-T", "--num_threads",
                        help="The number of threads for RAxML to use [ DEFAULT = 4 ]",
                        required=False,
                        default=4)

    args = parser.parse_args()
    args.mltreemap = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + os.sep

    # if len(args.code_name) < 5:
    #     sys.stderr.write("ERROR: code_name must be >= 5 characters!\n")
    #     sys.stderr.flush()
    #     sys.exit(-1)

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


def format_read_fasta(fasta_file, args, swappers = None):
    """
    Splits the input file into multiple files, each containing a maximum number of sequences as specified by the user.
    Ensures each sequence and sequence name is valid.
    :param fasta_file: 
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
                    if swappers and header in swappers.keys():
                            formatted_fasta_dict[swappers[header]] = sequence
                    else:
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
            accession = header_match.group(1)
            for mltree_id in dictionary:
                if dictionary[mltree_id].accession == accession:
                    short_id = dictionary[mltree_id].short_id
                    break
        elif header_type == "fungene":
            header_match = re.match("^>([A-Z0-9.]+)\s+coded_by=(.+),organism=(.+),definition=(.+)$", header)
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


def read_uc_present_options(uc_file):
    """
    Function to read a USEARCH cluster (.uc) file and present the user with optional headers from identical sequences
    :param uc_file: Path to a .uc file produced by USEARCH
    :return: swappers, a dictionary where the key is the representative cluster header and the value is the one to use
    """
    cluster_dict = dict()
    swappers = dict()
    candidates = dict()
    try:
        uc = open(uc_file, 'r')
    except:
        raise IOError("Unable to open " + uc_file + " for reading! Exiting...")

    line = uc.readline()
    # Find all clusters with multiple identical sequences
    while line:
        cluster_type, _, length, identity, _, _, _, cigar, header, representative = line.strip().split("\t")
        if cluster_type != "C":
            try:
                identity = float(identity)
            except ValueError:
                identity = "*"
            if cluster_type == "S":
                cluster_dict[header] = list()
            if cluster_type == "H" and identity == 100.0 and cigar == '=':
                cluster_dict[representative].append(header)
        line = uc.readline()

    # Present the headers of identical sequences to user
    for rep in cluster_dict:
        candidates.clear()
        subs = cluster_dict[rep]
        if len(subs) >= 1:
            sys.stderr.write("Found multiple identical sequences in cluster file:\n")
            candidates[str(1)] = rep
            acc = 2
            for candidate in subs:
                candidates[str(acc)] = candidate
                acc += 1
            for num in sorted(candidates.keys(), key=int):
                sys.stderr.write(num + ". " + candidates[num] + "\n")
            sys.stderr.flush()
            best = raw_input("Number of the best representative? ")
            while best not in candidates.keys():
                best = raw_input("Invalid number. Number of the best representative? ")
            if best != str(1):
                swappers[rep] = candidates[best]

    return swappers


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
    fungene_re = re.compile("^>([A-Z0-9.]+)\s+coded_by=(.+),organism=(.+),definition=(.+)$")
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
                sys.stderr.write("ERROR: Failed to parse suspected header from FunGenes repository!.\n")
                sys.stderr.write("Problem header: " + header + "\n")
                sys.stderr.flush()
                sys.exit()
            short_id = mltree_id + '_' + code_name
            ref_seq.short_id = short_id
        elif header_format == "ncbi":
            if re.match(">gi\|(\d+)\|(\w+)\|(\S+(\.\d+)*)\| (.*) \[(.*)\]$", header):
                ncbi_info = re.match(">gi\|(\d+)\|(\w+)\|(\S+(\.\d+)*)\| (.*) \[(.*)\]$", header)
                ref_seq.accession = ncbi_info.group(1)
                ref_seq.description = ncbi_info.group(6)
                # organism = ncbi_info.group(6)
            elif re.match(">gi\|(\d+)\|(\w+)\|(\S+(\.\d+)*)\|.*: (.*;) .*", header):
                ncbi_info = re.match(">gi\|(\d+)\|(\w+)\|(\S+(\.\d+)*)\|.*: (.*;) .*$", header)
                ref_seq.accession = ncbi_info.group(1)
                # definition = ncbi_info.group(5)
                ref_seq.description = "Unclassified"
            elif re.match(">gi\|(\d+)\|(\w+)\|(\S+(\.\d+)*)\|.*: Full=(.*)", header):
                ncbi_info = re.match(">gi\|(\d+)\|(\w+)\|(\S+(\.\d+)*)\|.*: Full=(.*)$", header)
                ref_seq.accession = ncbi_info.group(1)
                # definition = ncbi_info.group(5)
                ref_seq.description = "Unclassified"
            else:
                sys.stderr.write("ERROR: Failed to parse suspected header from Genbank!.\n")
                sys.stderr.write("Problem header: " + header + "\n")
                sys.stderr.flush()
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
        tree_tax_list_handle.write("%s\t%s | %s\n" % (mltree_id_key,
                                                    fasta_mltree_repl_dict[mltree_id_key].description,
                                                    fasta_mltree_repl_dict[mltree_id_key].accession))
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

    if args.uc:
        swappers = read_uc_present_options(args.uc)
    else:
        swappers = None
    input_fasta = args.fasta_file
    fasta_dict = format_read_fasta(input_fasta, args, swappers)
    code_name = args.code_name

    fasta_mltree_repl_dict = get_sequence_info(code_name, fasta_dict)

    tree_taxa_list = write_tax_ids(code_name, fasta_mltree_repl_dict)

    print "******************** %s generated ********************\n" % tree_taxa_list

    fasta_replaced = code_name + ".fc.repl.fasta"

    create_new_fasta(code_name, fasta_dict, fasta_replaced, fasta_mltree_repl_dict)

    print "************************** FASTA file, %s generated *************************\n" % fasta_replaced

    print "******************** Aligning the sequences using MUSCLE ********************\n"

    fasta_replaced_align = code_name + ".fc.repl.aligned.fasta"

    muscle_align_command = [args.executables["muscle"]]
    muscle_align_command += ["-in", fasta_replaced]
    muscle_align_command += ["-out", fasta_replaced_align]

    # print muscle_align_command, "\n"

    muscle_pro = subprocess.Popen(' '.join(muscle_align_command), shell=True, preexec_fn=os.setsid)
    muscle_pro.wait()
    if muscle_pro.returncode != 0:
        sys.stderr.write("ERROR: Multiple sequence alignment using " + args.executables["muscle"] +
                         " did not complete successfully! Command used:\n" + ' '.join(muscle_align_command) + "\n")
        sys.exit()

    fasta_mltree = code_name + ".fa"

    remove_dashes_from_msa(fasta_replaced_align, fasta_mltree)

    print "******************** FASTA file, %s generated ********************\n" % fasta_mltree

    makeblastdb_command = [args.executables["makeblastdb"]]
    makeblastdb_command += ["-in", fasta_mltree]
    makeblastdb_command += ["-out", fasta_mltree]
    makeblastdb_command += ["-dbtype", "prot"]
    makeblastdb_command += ["-input_type", "fasta"]

    makeblastdb_pro = subprocess.Popen(' '.join(makeblastdb_command), shell=True, preexec_fn=os.setsid)
    makeblastdb_pro.wait()
    if makeblastdb_pro.returncode != 0:
        sys.stderr.write("ERROR: BLAST database was unable to be made using " + args.executables["makeblastdb"] +
                         "! Command used:\n" + ' '.join(makeblastdb_command) + "\n")
        sys.exit()

    print "******************** BLAST DB for %s generated ********************\n" % code_name

    os.system("mv %s %s" % (fasta_replaced_align, fasta_mltree))

    hmm_build_command = [args.executables["hmmbuild"]]
    hmm_build_command += ["-s", code_name + ".hmm"]
    hmm_build_command.append(fasta_mltree)
    # hmm_build_command = "%s -s %s.hmm %s" %\
    #                     (args.executables["hmmbuild"], code_name, fasta_mltree)
    hmmbuild_pro = subprocess.Popen(' '.join(hmm_build_command), shell=True, preexec_fn=os.setsid)
    hmmbuild_pro.wait()

    if hmmbuild_pro.returncode != 0:
        sys.stderr.write("ERROR: hmmbuild did not complete successfully for:\n")
        sys.stderr.write(' '.join(hmm_build_command) + "\n")
        sys.exit()

    print "******************** HMM file for %s generated ********************\n" % code_name

    phylip_command = "java -cp %s/sub_binaries/readseq.jar run -a -f=12 %s" % (args.mltreemap, fasta_mltree)
    os.system(phylip_command)

    phylip_file = code_name + ".phy"
    os.system('mv %s.phylip %s' % (fasta_mltree, phylip_file))

    raxml_out = "%s_phy_files" % code_name

    if not os.path.exists(raxml_out):
        os.system("mkdir %s" % raxml_out)

    raxml_command = "%s -f a -p 12345 -x 12345 -# %s -m PROTGAMMAWAG -s %s -n %s -w %s -T %s" %\
                    (args.executables["raxmlHPC"], str(args.bootstraps), phylip_file, code_name,
                     args.mltreemap + raxml_out, args.num_threads)
    os.system(raxml_command)

    tree_to_swap = "%s/RAxML_bestTree.%s" % (raxml_out, code_name)
    final_mltree = "%s_tree.txt" % code_name
    os.system("mv %s %s" % (phylip_file, raxml_out))
    os.system("rm %s" % fasta_replaced)

    # TODO: Update the data/tree_data/cog_list.txt file with the new marker gene
    swap_tree_names(tree_to_swap, final_mltree, code_name)

    final_output_folder = "MLTreeMap_files_%s" % code_name
    os.system("mkdir %s" % final_output_folder)

    os.system("mv %s.fa %s.fa.p* %s" % (code_name, code_name, final_output_folder))
    os.system("mv %s.hmm %s %s %s" % (code_name, tree_taxa_list, final_mltree, final_output_folder))

    sys.stdout.write("Data for " + code_name + " has been generated succesfully.\n\n")
    sys.stdout.write("To integrate these data for use in TreeSAPP, the following steps must be performed:\n")
    sys.stdout.write("1. Modify data/tree_data/cog_list.tsv to include a properly formatted 'denominator' code\n")
    sys.stdout.write("2. $ cp " + final_output_folder + os.sep + "tax_ids_%s.txt" % code_name + " data/tree_data/\n")
    sys.stdout.write("3. $ cp " + final_output_folder + os.sep + code_name + "_tree.txt data/tree_data/\n")
    sys.stdout.write("4. $ cp " + final_output_folder + os.sep + code_name + ".hmm data/hmm_data/\n")
    sys.stdout.write("5. $ cp " + final_output_folder + os.sep + code_name + ".fa* data/alignment_data/\n")
    sys.stdout.write("6. $ cp " + final_output_folder + os.sep + code_name +
                     "_tree.txt MLTreeMap_imagemaker_2_061/tree_data/\n")
    sys.stdout.write("7. $ cp " + final_output_folder + os.sep + "tax_ids_%s.txt" % code_name +
                     " MLTreeMap_imagemaker_2_061/tree_data/\n")
    sys.stdout.write("8. Create a file called MLTreeMap_imagemaker_2_061/tree_data/domain_and_color_descriptions_" +
                     code_name + ".txt to add colours to clades in the new reference tree.\n")
    sys.stdout.write("9. Modify MLTreeMap_imagemaker_2_061/tree_data/drawing_info.txt following the obvious format\n")
    sys.stdout.flush()

main()
