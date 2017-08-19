#!/usr/bin/env python

import argparse
import sys
import os
import re
import subprocess
from time import gmtime, strftime
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
                        help="Unique name to be used by TreeSAPP internally. NOTE: Must be <=6 characters.\n"
                             "Refer to the first column of 'cog_list.txt' under the '#functional cogs' section)",
                        required=True)
    parser.add_argument("-m", "--min_length",
                        help="The minimum length of a protein to be used in building reference data. [ DEFAULT = 10 ]",
                        required=False,
                        default=10,
                        type=int)
    parser.add_argument("-i", "--identity",
                        help="The percent identity which the input sequences were clustered",
                        required=True,
                        type=str)
    parser.add_argument("-b", "--bootstraps",
                        help="The number of bootstrap replicates RAxML should perform [ DEFAULT = autoMR ]",
                        required=False,
                        default="autoMR")
    parser.add_argument("-T", "--num_threads",
                        help="The number of threads for RAxML to use [ DEFAULT = 4 ]",
                        required=False,
                        default=4)

    args = parser.parse_args()
    args.mltreemap = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + os.sep

    if len(args.code_name) > 6:
        sys.stderr.write("ERROR: code_name must be <= 6 characters!\n")
        sys.stderr.flush()
        sys.exit(-1)

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


def launch_write_command(cmd_list, collect_all=True):
    """
    Wrapper function for opening subprocesses through subprocess.Popen()
    :param cmd_list: A list of strings forming a complete command call
    :param collect_all: A flag determining whether stdout and stderr are returned 
    via stdout or just stderr is returned leaving stdout to be written to the screen
    :return: A string with stdout and/or stderr text and the returncode of the executable
    """
    stdout = ""
    if collect_all:
        proc = subprocess.Popen(' '.join(cmd_list),
                                shell=True,
                                preexec_fn=os.setsid,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT)
        stdout = proc.communicate()[0].decode("utf-8")
    else:
        proc = subprocess.Popen(' '.join(cmd_list),
                                shell=True,
                                preexec_fn=os.setsid)
        proc.wait()
    return stdout, proc.returncode


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
    fungene_gi_bad = re.compile("^>[0-9]+\s+coded_by=.+,organism=.+,definition=.+$")
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
            if fungene_gi_bad.match(header):
                sys.stderr.write("\nWARNING: Input sequences use 'GIs' which are obsolete and may be non-unique. "
                                 "For everyone's sanity, please download sequences with the `accno` instead.\n")
                sys.exit()

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
    # The regular expressions with the accession and organism name grouped

    gi_re = re.compile(">gi\|(\d+)\|[a-z]+\|\w.+\|(.*)$")
    gi_prepend_proper_re = re.compile(">gi\|([0-9]+)\|[a-z]+\|[_A-Z0-9.]+\|.*\[(.*)\]$")
    gi_prepend_mess_re = re.compile(">gi\|([0-9]+)\|pir\|\|(.*)$")
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

    header_format_regexi = [dbj_re, emb_re, gb_re, pdb_re, pir_re, ref_re, sp_re,
                            fungene_re, mltree_re, gi_prepend_proper_re, gi_prepend_mess_re, gi_re]
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
            tmp_ref_def = re.sub('[)(\[\]]', '', ref_seq.description)  # Remove parentheses for comparisons
            # This `swappers` is actually cluster_dict
            # keys are rep. headers, values are list of identical sequence names
            for header in swappers.keys():
                tmp_header = re.sub('[)(\[\]]', '', header)  # Remove parentheses for comparisons
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


def annotate_partition_tree(code_name, fasta_replace_dict, bipart_tree):
    try:
        tree_txt = open(bipart_tree, 'r')
    except IOError:
        raise IOError("Unable to open " + bipart_tree + " for reading.")

    tree = tree_txt.readline()
    tree_txt.close()
    for mltree_id_key in fasta_replace_dict.keys():
        tree = re.sub('\(' + mltree_id_key + "_" + code_name, '(' + fasta_replace_dict[mltree_id_key].description, tree)
        tree = re.sub(',' + mltree_id_key + "_" + code_name, ',' + fasta_replace_dict[mltree_id_key].description, tree)

    raxml_out = os.path.dirname(bipart_tree)
    annotated_tree_name = raxml_out + os.sep + "RAxML_bipartitions_annotated." + code_name
    try:
        annotated_tree = open(annotated_tree_name, 'w')
    except IOError:
        raise IOError("Unable to open " + annotated_tree_name + " for writing!")

    annotated_tree.write(tree)
    annotated_tree.close()

    return


def find_model_used(raxml_info_file):
    model_statement_re = re.compile(r".* model: ([A-Z]+) likelihood.*")
    model = ""
    with open(raxml_info_file) as raxml_info:
        for line in raxml_info:
            if model_statement_re.search(line):
                model = model_statement_re.search(line).group(1)
                break
            else:
                pass

    if model == "":
        sys.stderr.write("WARNING: Unable to parse model used from " + raxml_info_file + "!")
        sys.stderr.flush()
    return model


def update_build_parameters(args, code_name, aa_model):
    """
    Function to update the data/tree_data/ref_build_parameters.tsv file with information on this new reference sequence
    Format of file is "code_name       denominator     aa_model        cluster_identity        last_updated"
    :param args: 
    :param code_name: 
    :param aa_model: 
    :return: 
    """
    param_file = args.mltreemap + "data" + os.sep + "tree_data" + os.sep + "ref_build_parameters.tsv"
    try:
        params = open(param_file, 'a')
    except IOError:
        raise IOError("Unable to open " + param_file + "for appending!")

    date = strftime("%d_%b_%Y", gmtime())

    build_list = [code_name, "Z1111", "PROTGAMMA" + aa_model, args.identity, date]
    params.write("\t".join(build_list))

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

    final_output_folder = "MLTreeMap_files_%s" % code_name
    if not os.path.exists(final_output_folder):
        os.makedirs(final_output_folder)
    else:
        sys.stderr.write("WARNING: Output directory already exists. Previous outputs will be overwritten.\n")
        sys.stderr.flush()

    log = open("create_" + code_name + "_treesapp_data_log.txt", 'w')
    log.write("Command used:\n" + ' '.join(sys.argv) + "\n\n")

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

    stdout, muscle_pro_returncode = launch_write_command(muscle_align_command, False)

    if muscle_pro_returncode != 0:
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

    stdout, makeblastdb_pro_returncode = launch_write_command(makeblastdb_command)

    if makeblastdb_pro_returncode != 0:
        sys.stderr.write("ERROR: BLAST database was unable to be made using " + args.executables["makeblastdb"] +
                         "! Command used:\n" + ' '.join(makeblastdb_command) + "\n")
        sys.exit()

    log.write("\n### MAKEBLASTDB ###" + stdout)

    print "******************** BLAST DB for %s generated ********************\n" % code_name

    os.system("mv %s %s" % (fasta_replaced_align, fasta_mltree))

    hmm_build_command = [args.executables["hmmbuild"]]
    hmm_build_command += ["-s", code_name + ".hmm"]
    hmm_build_command.append(fasta_mltree)

    stdout, hmmbuild_pro_returncode = launch_write_command(hmm_build_command)

    if hmmbuild_pro_returncode != 0:
        sys.stderr.write("ERROR: hmmbuild did not complete successfully for:\n")
        sys.stderr.write(' '.join(hmm_build_command) + "\n")
        sys.exit()

    log.write("\n### HMMBUILD ###\n\n" + stdout)
    log.close()

    os.rename(code_name + ".hmm", final_output_folder + os.sep + code_name + ".hmm")

    print "******************** HMM file for %s generated ********************\n" % code_name

    phylip_command = "java -cp %s/sub_binaries/readseq.jar run -a -f=12 %s" % (args.mltreemap, fasta_mltree)
    os.system(phylip_command)

    phylip_file = code_name + ".phy"
    os.rename(fasta_mltree + ".phylip", phylip_file)

    raxml_out = "%s_phy_files" % code_name

    if not os.path.exists(raxml_out):
        os.system("mkdir %s" % raxml_out)
    else:
        sys.stderr.write("ERROR: " + raxml_out + " already exists from a previous run! "
                                                 "Please delete or rename it and try again.\n")
        sys.exit()

    raxml_command = [args.executables["raxmlHPC"]]
    raxml_command += ["-f", "a"]
    raxml_command += ["-p", "12345"]
    raxml_command += ["-x", "12345"]
    raxml_command += ["-#", args.bootstraps]
    raxml_command += ["-m", "PROTGAMMAAUTO"]
    raxml_command += ["-s", phylip_file]
    raxml_command += ["-n", code_name]
    raxml_command += ["-w", args.mltreemap + raxml_out]
    raxml_command += ["-T", args.num_threads]

    stdout, raxml_returncode = launch_write_command(raxml_command, False)

    if raxml_returncode != 0:
        sys.stderr.write("ERROR: RAxML did not complete successfully! "
                         "Look in " + args.mltreemap + raxml_out + "RAxML_info." + code_name + " for error message.")
        sys.exit(3)

    tree_to_swap = "%s/RAxML_bestTree.%s" % (raxml_out, code_name)
    final_mltree = "%s_tree.txt" % code_name
    os.system("mv %s %s" % (phylip_file, raxml_out))
    os.system("rm %s %s" % (fasta_replaced_file, phylip_file + ".reduced"))

    swap_tree_names(tree_to_swap, final_mltree, code_name)

    os.system("mv %s.fa %s.fa.p* %s" % (code_name, code_name, final_output_folder))
    os.system("mv %s %s %s" % (tree_taxa_list, final_mltree, final_output_folder))

    annotate_partition_tree(code_name, fasta_replace_dict, raxml_out + os.sep + "RAxML_bipartitions." + code_name)
    aa_model = find_model_used(raxml_out + os.sep + "RAxML_info." + code_name)
    update_build_parameters(args, code_name, aa_model)

    sys.stdout.write("Data for " + code_name + " has been generated succesfully.\n\n")
    sys.stdout.write("To integrate these data for use in TreeSAPP, the following steps must be performed:\n")
    sys.stdout.write("1. Include properly formatted 'denominator' codes "
                     "in data/tree_data/cog_list.tsv and data/tree_data/ref_build_parameters.tsv\n")
    sys.stdout.write("2. $ cp " + final_output_folder + os.sep + "tax_ids_%s.txt" % code_name + " data/tree_data/\n")
    sys.stdout.write("3. $ cp " + final_output_folder + os.sep + code_name + "_tree.txt data/tree_data/\n")
    sys.stdout.write("4. $ cp " + final_output_folder + os.sep + code_name + ".hmm data/hmm_data/\n")
    sys.stdout.write("5. $ cp " + final_output_folder + os.sep + code_name + ".fa* data/alignment_data/\n")
    sys.stdout.write("6. $ cp " + final_output_folder + os.sep + code_name +
                     "_tree.txt imagemaker_2_061/tree_data/\n")
    sys.stdout.write("7. $ cp " + final_output_folder + os.sep + "tax_ids_%s.txt" % code_name +
                     " imagemaker_2_061/tree_data/\n")
    sys.stdout.write("8. Create a file called imagemaker_2_061/tree_data/domain_and_color_descriptions_" +
                     code_name + ".txt to add colours to clades in the new reference tree.\n")
    sys.stdout.write("9. Modify imagemaker_2_061/tree_data/drawing_info.txt following the obvious format\n")
    sys.stdout.flush()

main()