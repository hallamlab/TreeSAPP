#!/usr/bin/env python

import argparse
import sys
import os
import re
import subprocess
from mltreemap import os_type, is_exe, which


class ReferenceSequence:
    def __init__(self):
        self.accession = ""
        self.description = ""
        self.organism = ""
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


def format_read_fasta(fasta_file, args):
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


def create_new_fasta(out_fasta, ref_seq_dict, dashes=True):
    """
    Writes a new FASTA file using a dictionary of ReferenceSequence class objects
    :param out_fasta: Name of the FASTA file to write to
    :param ref_seq_dict: Dictionary containing ReferenceSequence objects, numbers are keys
    :param dashes: Flag indicating whether hyphens should be removed from sequences
    :return:
    """
    out_fasta_handle = open(out_fasta, "w")
    num_seqs_written = 0

    for mltree_id in sorted(ref_seq_dict, key=int):
        ref_seq = ref_seq_dict[mltree_id]
        if dashes:
            sequence = re.sub('-', '', ref_seq.sequence)
        else:
            sequence = ref_seq.sequence
        out_fasta_handle.write(">" + ref_seq.short_id + "\n")
        out_fasta_handle.write(sequence + "\n")
        num_seqs_written += 1

    out_fasta_handle.close()

    if num_seqs_written == 0:
        sys.stderr.write("ERROR: No sequences written to " + out_fasta + ".\n")
        sys.stderr.write("The headers in your input file are probably not accommodated in the regex patterns used. "
                         "Function responsible: get_header_format. Please make an issue on the GitHub page.\n")
        sys.stderr.flush()
        sys.exit(5)

    return


def read_uc(uc_file):
    """
    Function to read a USEARCH cluster (.uc) file
    :param uc_file: Path to a .uc file produced by USEARCH
    :return: Dictionary where keys are representative cluster headers and the values are headers of identical sequences
    """
    cluster_dict = dict()
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
                cluster_dict['>' + header] = list()
            if cluster_type == "H" and identity == 100.0 and cigar == '=':
                cluster_dict['>' + representative].append('>' + header)
        line = uc.readline()
    return cluster_dict


def present_cluster_rep_options(cluster_dict):
    """
    Present the headers of identical sequences to user for them to decide on representative header
    :param cluster_dict: dictionary from read_uc(uc_file)
    :return:
    """
    swappers = dict()
    candidates = dict()

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
            # Useful for testing - no need to pick which sequence name is best!
            # best = str(1)
            while best not in candidates.keys():
                best = raw_input("Invalid number. Number of the best representative? ")
            if best != str(1):
                swappers['>' + rep] = '>' + candidates[best]

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
    # The regular expressions with the accession and organism name grouped
    gi_re = re.compile(">gi\|(\d+)\|(\w+)\|(\S+(\.\d+)*)\|(.*)$")
    dbj_re = re.compile(">dbj\|(.*)\|.*\[(.*)\]")
    emb_re = re.compile(">emb\|(.*)\|.*\[(.*)\]")
    gb_re = re.compile(">gb\|(.*)\|.*\[(.*)\]")
    ref_re = re.compile(">ref\|(.*)\|.*\[(.*)\]")
    pdb_re = re.compile(">pdb\|(.*)\|(.*)$")
    pir_re = re.compile(">pir\|\|(\w+).* - (.*)$")
    sp_re = re.compile(">sp\|(.*)\|.*Full=(.*); AltName:.*$")
    fungene_re = re.compile("^>([A-Z0-9.]+)\s+coded_by=(.+),organism=(.+),definition=(.+)$")
    # TODO: Find the description field for the mltree_re
    mltree_re = re.compile("^>(\d+)_" + re.escape(code_name))

    header_format_regexi = [gi_re, dbj_re, emb_re, gb_re, pdb_re, pir_re, ref_re, sp_re, fungene_re, mltree_re]
    for regex in header_format_regexi:
        if regex.match(header):
            return regex

    return None


def get_sequence_info(code_name, fasta_dict, fasta_replace_dict, swappers=None):
    """
    This function is used to find the accession ID and description of each sequence from the FASTA file
    :param code_name: code_name from the command-line parameters
    :param fasta_dict: a dictionary with headers as keys and sequences as values (returned by format_read_fasta)
    :param fasta_replace_dict:
    :param swappers:
    :return: fasta_replace_dict with a complete ReferenceSequence() value for every mltree_id key
    """

    print "******************** Generating information for formatting purposes ********************"
    mltree_id_accumulator = 1
    if len(fasta_replace_dict.keys()) > 0:
        for mltree_id in fasta_replace_dict:
            ref_seq = fasta_replace_dict[mltree_id]
            ref_seq.short_id = mltree_id + '_' + code_name
            tmp_ref_def = re.sub('\)|\(', '', ref_seq.description) # Remove parentheses for comparisons
            # This `swappers` is actually cluster_dict
            # keys are rep. headers, values are list of identical sequence names
            for header in swappers.keys():
                tmp_header = re.sub('\)|\(', '', header) # Remove parentheses for comparisons
                # Need to check both keys and values since it is unknown whether the rep was selected or not
                if re.search(ref_seq.accession, header):
                    if re.search(tmp_ref_def, tmp_header):
                        ref_seq.sequence = fasta_dict[header]
                if len(swappers[header]) > 0:
                    for constituent in swappers[header]:
                        if re.search(tmp_ref_def, constituent) and re.search(ref_seq.accession, constituent):
                            ref_seq.sequence = fasta_dict[header]

            if not ref_seq.sequence:
                sys.exit("Unable to find header for " + ref_seq.accession)

    else:  # if fasta_replace_dict needs to be populated, this is a new run
        for header in fasta_dict.keys():
            mltree_id = str(mltree_id_accumulator)
            ref_seq = ReferenceSequence()
            ref_seq.sequence = fasta_dict[header]
            if swappers and header in swappers.keys():
                header = swappers[header]
            header_format_re = get_header_format(header, code_name)
            if header_format_re is None:
                raise AssertionError("Unable to parse header: " + header)
            sequence_info = header_format_re.match(header)
            if sequence_info:
                if len(sequence_info.groups()) == 2:
                    ref_seq.accession = sequence_info.group(1)
                    ref_seq.description = sequence_info.group(2)
                elif len(sequence_info.groups()) == 4:
                    # From FunGenes
                    ref_seq.accession = sequence_info.group(1)
                    ref_seq.locus = sequence_info.group(2)
                    # ref_seq.organism = sequence_info.group(3)
                    ref_seq.description = sequence_info.group(3)
            else:
                print "Unable to handle header: ", header
                sys.exit()

            ref_seq.short_id = mltree_id + '_' + code_name
            fasta_replace_dict[mltree_id] = ref_seq

            mltree_id_accumulator += 1

    return fasta_replace_dict


def write_tax_ids(fasta_replace_dict, tree_taxa_list):
    """
    Write the number, organism and accession ID, if possible
    :param fasta_replace_dict:
    :param tree_taxa_list: The name of the output file
    :return:
    """
    tree_tax_list_handle = open(tree_taxa_list, "w")
    for mltree_id_key in sorted(fasta_replace_dict.keys(), key=int):
        tree_tax_list_handle.write("%s\t%s | %s\n" % (mltree_id_key,
                                                      fasta_replace_dict[mltree_id_key].description,
                                                      fasta_replace_dict[mltree_id_key].accession))

    tree_tax_list_handle.close()
    return


def read_tax_ids(tree_taxa_list):
    """
    Reads the taxonomy and accession ID affiliated with each sequence number.
    This information is used to avoid horrible manual work if the pipeline is ran multiple times
    :param tree_taxa_list: The name of the tax_ids file to read
    :return:
    """
    try:
        tree_tax_list_handle = open(tree_taxa_list, 'r')
    except:
        raise IOError("Unable to open " + tree_taxa_list + " for reading! Exiting.")
    fasta_replace_dict = dict()
    line = tree_tax_list_handle.readline()
    while line:
        mltree_id_key, seq_info = line.strip().split("\t")
        ref_seq = ReferenceSequence()
        ref_seq.description = seq_info.split(" | ")[0]
        ref_seq.accession = seq_info.split(" | ")[1]
        fasta_replace_dict[mltree_id_key] = ref_seq
        line = tree_tax_list_handle.readline()
    tree_tax_list_handle.close()

    return fasta_replace_dict


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
    tree_taxa_list = "tax_ids_%s.txt" % code_name

    if args.uc and os.path.exists(tree_taxa_list):
        use_previous_names = raw_input(tree_taxa_list + " found from a previous attempt. "
                                                        "Should it be used for this run? [y|n] ")
        while use_previous_names != "y" and use_previous_names != "n":
            use_previous_names = raw_input("Incorrect response. Please input either 'y' or 'n'. ")
    else:
        use_previous_names = 'n'

    fasta_replace_dict = dict()
    if args.uc:
        cluster_dict = read_uc(args.uc)
        if use_previous_names == 'n':
            swappers = present_cluster_rep_options(cluster_dict)
            fasta_replace_dict = get_sequence_info(code_name, fasta_dict, fasta_replace_dict, swappers)
            write_tax_ids(fasta_replace_dict, tree_taxa_list)
        if use_previous_names == 'y':
            fasta_replace_dict = read_tax_ids(tree_taxa_list)
            if len(fasta_replace_dict.keys()) != len(fasta_dict.keys()):
                raise AssertionError("Number of sequences in new FASTA input and " + tree_taxa_list + " are not equal!")
            fasta_replace_dict = get_sequence_info(code_name, fasta_dict, fasta_replace_dict, cluster_dict)
    else:
        # args.uc is None and use_previous_names == 'n'
        fasta_replace_dict = get_sequence_info(code_name, fasta_dict, fasta_replace_dict)
        write_tax_ids(fasta_replace_dict, tree_taxa_list)

    print "******************** %s generated ********************\n" % tree_taxa_list
    
    fasta_replaced_file = code_name + ".fc.repl.fasta"

    create_new_fasta(fasta_replaced_file, fasta_replace_dict)

    print "************************** FASTA file, %s generated *************************\n" % fasta_replaced_file

    print "******************** Aligning the sequences using MUSCLE ********************\n"

    fasta_replaced_align = code_name + ".fc.repl.aligned.fasta"

    muscle_align_command = [args.executables["muscle"]]
    muscle_align_command += ["-in", fasta_replaced_file]
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
    os.system("rm %s" % fasta_replaced_file)

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
