#!/usr/bin/env python3


__author__ = "Connor Morgan-Lang and Kishori Konwar"
__maintainer__ = "Connor Morgan-Lang"
__license__ = "GPL-3.0"
__version__ = "1.1.0"

try:
    import profile
    import argparse
    import sys
    import os
    import shutil
    import re
    import glob
    import time
    import traceback
    import subprocess
    import logging
    from ete3 import Tree
    from multiprocessing import Pool, Process, Lock, Queue, JoinableQueue
    from os import path
    from os import listdir
    from os.path import isfile, join
    from time import gmtime, strftime

    from utilities import Autovivify, os_type, which, find_executables, generate_blast_database, clean_lineage_string,\
        reformat_string, available_cpu_count, write_phy_file, reformat_fasta_to_phy, profile_aligner, run_papara, \
        launch_evolutionary_placement_queries
    from classy import CommandLineWorker, CommandLineFarmer, ItolJplace, NodeRetrieverWorker,\
        TreeLeafReference, TreeProtein, ReferenceSequence, prep_logging
    from fasta import format_read_fasta, get_headers, write_new_fasta, trim_multiple_alignment, read_fasta_to_dict
    from entish import create_tree_info_hash, deconvolute_assignments, read_and_understand_the_reference_tree,\
        get_node, annotate_partition_tree, find_cluster
    from external_command_interface import launch_write_command, setup_progress_bar
    from lca_calculations import *
    from jplace_utils import *
    from file_parsers import *
    from phylo_dist import *
    from update_refpkg import CreateFuncTreeUtility

    import _tree_parser
    import _fasta_reader
except ImportWarning:
    sys.stderr.write("Could not load some user defined module functions")
    sys.stderr.write(traceback.print_exc(10))
    sys.exit(3)


def get_options():
    """
    Returns the parser to interpret user options.
    """
    parser = argparse.ArgumentParser(description='Phylogenetically informed insertion of sequence into a reference tree'
                                                 ' using a Maximum Likelihood algorithm.')
    parser.add_argument('-i', '--fasta_input', required=True,
                        help='Your sequence input file in FASTA format')
    parser.add_argument('-o', '--output', default='./output/', required=False,
                        help='output directory [DEFAULT = ./output/]')
    parser.add_argument('-c', '--composition', default="meta", choices=["meta", "single"],
                        help="Sample composition being either a single organism or a metagenome.")
    parser.add_argument("--trim_align", default=False, action="store_true",
                        help="Flag to turn on position masking of the multiple sequence alignmnet [DEFAULT = False]")
    parser.add_argument('-g', '--min_seq_length', default=30, type=int,
                        help='minimal sequence length after alignment trimming [DEFAULT = 30]')
    parser.add_argument('-R', '--reftree', default='p', type=str,
                        help='Reference tree (p = MLTreeMap reference phylogenetic tree [DEFAULT])'
                             ' Change to code to map query sequences to specific phylogenetic tree.')
    parser.add_argument('-t', '--targets', default='ALL', type=str,
                        help='A comma-separated list specifying which marker genes to query in input by'
                             ' the "denominator" column in data/tree_data/cog_list.tsv'
                             ' - e.g., M0701,D0601 for mcrA and nosZ\n[DEFAULT = ALL]')
    parser.add_argument('-m', '--molecule', default='dna', choices=['prot', 'dna'],
                        help='the type of input sequences (prot = Protein; dna = Nucleotide [DEFAULT])')
    parser.add_argument("-l", "--min_likelihood", default=0.2, type=float,
                        help="The minimum likelihood weight ratio required for a RAxML placement. "
                             "[DEFAULT = 0.2]")
    parser.add_argument("-P", "--placement_parser", default="best", type=str, choices=["best", "lca"],
                        help="Algorithm used for parsing each sequence's potential RAxML placements. "
                             "[DEFAULT = 'best']")

    rpkm_opts = parser.add_argument_group('RPKM options')
    rpkm_opts.add_argument("--rpkm", action="store_true", default=False,
                           help="Flag indicating RPKM values should be calculated for the gene sequences detected")
    rpkm_opts.add_argument("-r", "--reads", required=False,
                           help="FASTQ file containing to be aligned to predicted genes using BWA MEM")
    rpkm_opts.add_argument("-2", "--reverse", required=False,
                           help="FASTQ file containing to reverse mate-pair reads to be aligned using BWA MEM")
    rpkm_opts.add_argument("-p", "--pairing", required=False, default='pe', choices=['pe', 'se'],
                           help="Indicating whether the reads are paired-end (pe) or single-end (se)")

    update_tree = parser.add_argument_group('Update-tree options')
    # output will by treesapp_output/update_tree
    update_tree.add_argument("--update_tree", action="store_true", default=False,
                             help="Flag indicating the reference tree specified by `--reftree` "
                                  "is to be updated using the sequences found in TreeSAPP output")
    update_tree.add_argument("--uclust", required=False, default=False, action="store_true",
                             help="Cluster sequences that mapped to the reference tree prior to updating")
    # update_tree.add_argument("-u", "--uclust_identity", required=False, default=0.97, type=float,
    #                          help="Sequence identity value to be used in uclust [DEFAULT = 0.97]")
    update_tree.add_argument("-a", "--alignment_mode", required=False, default='d', type=str, choices=['d', 'p'],
                             help="Alignment mode: 'd' for default and 'p' for profile-profile alignment")

    miscellaneous_opts = parser.add_argument_group("Miscellaneous options")
    miscellaneous_opts.add_argument("--reclassify", action="store_true", default=False,
                                    help="Flag indicating current outputs should be used to generate "
                                         "all outputs downstream of phylogenetic placement.")
    miscellaneous_opts.add_argument('--overwrite', action='store_true', default=False,
                                    help='overwrites previously processed output folders')
    miscellaneous_opts.add_argument('-v', '--verbose', action='store_true',  default=False,
                                    help='Prints a more verbose runtime log')
    miscellaneous_opts.add_argument("--check_trees", action="store_true", default=False,
                                    help="Quality-check the reference trees before running TreeSAPP")
    miscellaneous_opts.add_argument('-T', '--num_threads', default=2, type=int,
                                    help='specifies the number of CPU threads to use in RAxML and BLAST '
                                         'and processes throughout the pipeline [DEFAULT = 2]')
    miscellaneous_opts.add_argument('-d', '--delete', default=False, action="store_true",
                                    help='Delete intermediate file to save disk space\n'
                                         'Recommended for large metagenomes!')

    args = parser.parse_args()

    return args


def check_parser_arguments(args):
    """
    Ensures the command-line arguments returned by argparse are sensible
    :param args: object with parameters returned by argparse.parse_args()
    :return: 'args', a summary of TreeSAPP settings.
    """

    # Add (or replace a trailing (back)slash with) the os.sep to the end of the output directory
    while re.search(r'/\Z', args.output) or re.search(r'\\\Z', args.output):
        args.output = args.output[:-1]
    args.output += os.sep
    if not os.path.isdir(args.output):
        os.makedirs(args.output)

    # Setup the global logger and main log file
    log_file_name = args.output + os.sep + "TreeSAPP_log.txt"
    prep_logging(log_file_name, args.verbose)

    # Ensure files contain more than 0 sequences
    args.treesapp = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + os.sep

    # Set the reference data file prefix and the reference tree name
    args.reference_tree = args.reftree

    args.targets = args.targets.split(',')
    if args.targets != ['ALL']:
        for marker in args.targets:
            if not re.match('[A-Z][0-9]{4}', marker):
                logging.error("Incorrect format for target: " + str(marker) +
                              "\nRefer to column 'Denominator' in " + args.treesapp + "data/tree_data/" +
                              "cog_list.tsv for identifiers that can be used.\n")
                sys.exit()

    args = find_executables(args)

    if sys.version_info > (2, 9):
        args.py_version = 3
    else:
        args.py_version = 2

    args.output_dir_var = args.output + 'various_outputs' + os.sep
    args.output_dir_raxml = args.output + 'final_RAxML_outputs' + os.sep
    args.output_dir_final = args.output + 'final_outputs' + os.sep

    if args.num_threads > available_cpu_count():
        logging.warning("Number of threads specified is greater than those available! "
                        "Using maximum threads available (" + str(available_cpu_count()) + ")\n")
        args.num_threads = available_cpu_count()

    if args.rpkm:
        if not args.reads:
            logging.error("At least one FASTQ file must be provided if -rpkm flag is active!")
            sys.exit()
        if args.reverse and not args.reads:
            logging.error("File containing reverse reads provided but forward mates file missing!")
            sys.exit()

    if args.molecule == "prot" and args.rpkm:
        logging.error("Unable to calculate RPKM values for protein sequences.\n")
        sys.exit()

    # Parameterizing the hmmsearch output parsing:
    args.min_acc = 0.7
    args.min_e = 0.0001
    args.perc_aligned = 15

    return args


def check_previous_output(args):
    """
    Prompts the user to determine how to deal with a pre-existing output directory.
    :rtype: Namespace object
    :param args: Command-line argument object from get_options and check_parser_arguments
    :return An updated version of 'args', a summary of TreeSAPP settings.
    """

    main_output_dirs = [args.output_dir_var, args.output_dir_raxml, args.output_dir_final]
    args.skip = 'n'

    if os.path.isdir(args.output_dir_final) and args.overwrite:
        for output_dir in main_output_dirs:
            shutil.rmtree(output_dir)
            os.mkdir(output_dir)
    elif os.path.isdir(args.output_dir_final) and not args.overwrite:
        workflows = list()
        if args.update_tree:
            if os.path.isfile(args.output_dir_final + os.sep + "marker_contig_map.tsv"):
                args.skip = 'y'
                workflows.append("updating")
            else:
                logging.warning("Update-tree impossible as " + args.output + " is missing input files.\n")
        elif args.rpkm:
            if os.path.isfile(args.reads) and os.path.isfile(args.output_dir_final + os.sep + "marker_contig_map.tsv"):
                args.skip = 'y'
                workflows.append("calculating RPKM")
            else:
                logging.warning("RPKM impossible as " + args.output + " is missing input files.\n")
        elif args.reclassify:
            if os.path.isdir(args.output_dir_var):
                jplace_files = glob.glob(args.output_dir_var + os.sep + "*jplace")
                if len(jplace_files) >= 1:
                    args.skip = 'y'
                    workflows.append("reclassifying")
                else:
                    logging.warning("reclassify impossible as " + args.output + " is missing input files.\n")
        else:
            # Warn user then remove all main output directories, leaving log in output
            logging.warning("Removing previous outputs in '" + args.output + "'. " +
                            "You have 10 seconds to hit Ctrl-C before this proceeds.\n")
            time.sleep(10)
            for output_dir in main_output_dirs:
                shutil.rmtree(output_dir)
                os.mkdir(output_dir)

        if len(workflows) >= 1:
            sys.stdout.write(','.join(workflows) + ".\n")
            sys.stdout.flush()
    else:
        # Create the output directories
        for output_dir in main_output_dirs:
            os.mkdir(output_dir)

    return args


def get_hmm_length(hmm_file):
    """
    Function to open the ref_tree's hmm file and determine its length
    :param hmm_file: The HMM file produced by hmmbuild to parse for the HMM length
    :return: The length (int value) of the HMM
    """
    try:
        hmm = open(hmm_file, 'r')
    except IOError:
        raise IOError("Unable to open " + hmm_file + " for reading! Exiting.")

    line = hmm.readline()
    length = 0
    while line:
        # LENG XXX
        if re.match(r"^LENG\s+([0-9]+)", line):
            length = int(line.split()[1])
        line = hmm.readline()
    if length > 0:
        return length
    else:
        raise AssertionError("Unable to parse the HMM length from " + hmm_file + ". Exiting.")


def align_ref_queries(args, new_ref_queries, update_tree):
    """
    Function queries the candidate set of proteins to be used for updating the tree against the reference set
    The output feeds into find_novel_refs. Necessary to determine whether there are interesting new proteins or
    just more of the same
    :param args: Command-line argument object from get_options and check_parser_arguments
    :param new_ref_queries:
    :param update_tree:
    :return: Path to tabular alignment file
    """
    alignments = update_tree.Output + "candidate_alignments.tsv"

    ref_fasta = os.sep.join([args.treesapp, "data",  "alignment_data",  update_tree.COG + ".fa"])
    db_prefix = update_tree.Output + os.sep + update_tree.COG
    # Make a temporary BLAST database to see what is novel
    # Needs a path to write the temporary unaligned FASTA file
    generate_blast_database(args, ref_fasta, "prot", db_prefix)

    logging.info("Aligning the candidate sequences to the current reference sequences using blastp... ")

    align_cmd = [args.executables["blastp"]]
    align_cmd += ["-query", new_ref_queries]
    align_cmd += ["-db", db_prefix + ".fa"]
    align_cmd += ["-outfmt", str(6)]
    align_cmd += ["-out", alignments]
    align_cmd += ["-num_alignments", str(1)]

    launch_write_command(align_cmd)

    # Remove the temporary BLAST database
    db_suffixes = ['', ".phr", ".pin", ".psq"]
    for db_file in db_suffixes:
        if os.path.isfile(db_prefix + ".fa" + db_file):
            os.remove(db_prefix + ".fa" + db_file)

    logging.info("done.\n")

    return alignments


def find_novel_refs(ref_candidate_alignments, aa_dictionary, create_func_tree):
    new_refs = dict()
    try:
        alignments = open(ref_candidate_alignments, 'r')
    except IOError:
        raise IOError("Unable to open " + ref_candidate_alignments + " for reading! Exiting.")

    line = alignments.readline()
    while line:
        fields = line.split("\t")
        if float(fields[2]) <= create_func_tree.cluster_id:
            query = '>' + fields[0]
            new_refs[query] = aa_dictionary[query]
        else:
            pass
        line = alignments.readline()

    alignments.close()
    return new_refs


def build_hmm(args, msa_file, hmm_output):

    sys.stdout.write("Building HMM... ")
    sys.stdout.flush()

    if os.path.isfile(hmm_output):
        os.remove(hmm_output)
    command = [args.executables["hmmbuild"], hmm_output, msa_file]

    launch_write_command(command)

    sys.stdout.write("done.\n")

    return


def validate_inputs(args, marker_build_dict):
    """
    This function filters the files in data/alignment_data/ for sequences that are entirely ambiguity characters
    or if there are any sequences in the MSA that are not the consistent length
    :param args: the command-line and default options
    :param marker_build_dict: A dictionary (indexed by marker 5-character 'denominator's) mapping MarkerBuild objects
    :return: list of files that were edited
    """
    logging.info("Testing validity of reference trees... ")
    ref_trees = glob.glob(args.treesapp + os.sep + "data/tree_data/*_tree.txt")
    ref_tree_dict = dict()
    for tree_file in ref_trees:
        marker = os.path.basename(tree_file).strip("_tree.txt").strip("_")
        for denominator in marker_build_dict:
            if marker_build_dict[denominator].cog == marker:
                ref_tree_dict[denominator] = tree_file
    status = pparse_ref_trees(denominator_ref_tree_dict=ref_tree_dict, args=args)
    logging.info("done.\n")
    if status is None:
        logging.error("Reference trees do not appear to be formatted correctly!\n")
        sys.exit(3)
    else:
        logging.info("Reference trees appear to be formatted correctly. Continuing with TreeSAPP.\n")
    return


def run_prodigal(args, fasta_file, output_file, nucleotide_orfs=None):
    prodigal_command = [args.executables["prodigal"]]
    prodigal_command += ["-i", fasta_file]
    prodigal_command += ["-p", "meta"]
    prodigal_command += ["-a", output_file]
    if nucleotide_orfs:
        prodigal_command += ["-d", nucleotide_orfs]
    stdout, proc_code = launch_write_command(prodigal_command)

    if proc_code != 0:
        logging.error("Prodigal did not complete successfully!\n" +
                      "Command used:\n" + ' '.join(prodigal_command), "err", "\n")
        sys.exit(3)
    return


def predict_orfs(args):
    """
    Predict ORFs from the input FASTA file using Prodigal

    :param args: Command-line argument object from get_options and check_parser_arguments
    :return:
    """

    logging.info("Predicting open-reading frames in the genomes using Prodigal... ")

    start_time = time.time()

    sample_prefix = '.'.join(os.path.basename(args.fasta_input).split('.')[:-1])
    if args.num_threads > 1 and args.composition == "meta":
        # Split the input FASTA into num_threads files to run Prodigal in parallel
        # TODO: low-priority - split into multiple files based on total sequence length rather than number of sequences
        input_fasta_dict = format_read_fasta(args.fasta_input, args.molecule, args.output)
        n_seqs = len(input_fasta_dict.keys())
        chunk_size = int(n_seqs / args.num_threads) + (n_seqs % args.num_threads)
        split_files = write_new_fasta(input_fasta_dict,
                                      args.output_dir_var + sample_prefix,
                                      chunk_size)
    else:
        split_files = [args.fasta_input]

    task_list = list()
    for fasta_chunk in split_files:
        chunk_prefix = args.output_dir_final + '.'.join(os.path.basename(fasta_chunk).split('.')[:-1])
        prodigal_command = [args.executables["prodigal"]]
        prodigal_command += ["-i", fasta_chunk]
        prodigal_command += ["-p", args.composition]
        prodigal_command += ["-a", chunk_prefix + "_ORFs.faa"]
        prodigal_command += ["-d", chunk_prefix + "_ORFs.fna"]
        prodigal_command += ["1>/dev/null", "2>/dev/null"]
        task_list.append(prodigal_command)

    num_tasks = len(task_list)
    if num_tasks > 0:
        cl_farmer = CommandLineFarmer("Prodigal -p " + args.composition, args.num_threads)
        cl_farmer.add_tasks_to_queue(task_list)

        cl_farmer.task_queue.close()
        cl_farmer.task_queue.join()

    # Concatenate outputs
    aa_orfs_file = args.output_dir_final + sample_prefix + "_ORFs.faa"
    nuc_orfs_file = args.output_dir_final + sample_prefix + "_ORFs.fna"
    if not os.path.isfile(aa_orfs_file) and not os.path.isfile(nuc_orfs_file):
        tmp_prodigal_aa_orfs = glob.glob(args.output_dir_final + sample_prefix + "*_ORFs.faa")
        tmp_prodigal_nuc_orfs = glob.glob(args.output_dir_final + sample_prefix + "*_ORFs.fna")
        os.system("cat " + ' '.join(tmp_prodigal_aa_orfs) + " > " + aa_orfs_file)
        os.system("cat " + ' '.join(tmp_prodigal_nuc_orfs) + " > " + nuc_orfs_file)
        intermediate_files = list(tmp_prodigal_aa_orfs + tmp_prodigal_nuc_orfs + split_files)
        for tmp_file in intermediate_files:
            os.remove(tmp_file)

    logging.info("done.\n")

    args.fasta_input = aa_orfs_file
    args.nucleotide_orfs = nuc_orfs_file

    end_time = time.time()
    hours, remainder = divmod(end_time - start_time, 3600)
    minutes, seconds = divmod(remainder, 60)
    logging.debug("\tProdigal time required: " +
                  ':'.join([str(hours), str(minutes), str(round(seconds, 2))]) + "\n")

    return args


def hmmsearch_orfs(args, marker_build_dict):
    hmm_domtbl_files = list()
    nucl_target_hmm_files = list()
    prot_target_hmm_files = list()

    # Find all of the available HMM profile files
    hmm_dir = args.treesapp + os.sep + 'data' + os.sep + "hmm_data" + os.sep
    try:
        os.path.isdir(hmm_dir)
    except IOError:
        logging.error(hmm_dir + "does not exist!")
        sys.exit(3)
    hmm_files = glob.glob(hmm_dir + "*.hmm")

    if len(hmm_files) == 0:
        logging.error(hmm_dir + "does not contain any files with '.hmm' extension... so no HMMs.\n")
        sys.exit(3)

    # Filter the HMM files to only the target markers
    for marker_code in marker_build_dict:
        ref_marker = marker_build_dict[marker_code]
        hmm_profile = hmm_dir + ref_marker.cog + ".hmm"
        if hmm_profile not in hmm_files:
            logging.error("Unable to locate HMM-profile for " + ref_marker.cog + "(" + marker_code + ").\n")
        else:
            if ref_marker.molecule == "prot":
                prot_target_hmm_files.append(hmm_profile)
            else:
                nucl_target_hmm_files.append(hmm_profile)

    acc = 0.0
    logging.info("Searching for marker proteins in ORFs using hmmsearch.\n")
    step_proportion = setup_progress_bar(len(prot_target_hmm_files) + len(nucl_target_hmm_files))

    # Create and launch the hmmsearch commands iteratively.
    # Since its already rippin' fast, don't need to run in parallel
    hmmsearch_command_base = [args.executables["hmmsearch"]]
    hmmsearch_command_base += ["--cpu", str(args.num_threads)]
    hmmsearch_command_base.append("--noali")
    for hmm_file in prot_target_hmm_files:
        rp_marker = re.sub(".hmm", '', os.path.basename(hmm_file))
        domtbl = args.output_dir_var + rp_marker + "_to_ORFs_domtbl.txt"
        hmm_domtbl_files.append(domtbl)
        final_hmmsearch_command = hmmsearch_command_base + ["--domtblout", domtbl]
        final_hmmsearch_command += [hmm_file, args.formatted_input_file]
        stdout, ret_code = launch_write_command(final_hmmsearch_command)
        if ret_code != 0:
            logging.error("hmmsearch did not complete successfully! Output:\n" + stdout + "\n" +
                          "Command used:\n" + ' '.join(final_hmmsearch_command) + "\n")
            sys.exit(3)

        # Update the progress bar
        acc += 1.0
        if acc >= step_proportion:
            acc -= step_proportion
            time.sleep(0.1)
            sys.stdout.write("-")
            sys.stdout.flush()

    sys.stdout.write("-]\n")
    return hmm_domtbl_files


def extract_hmm_matches(args, hmm_matches: dict, fasta_dict: dict):
    """
    Function writes the sequences identified by the HMMs to output files in FASTA format.
    Full-length query sequences with homologous regions are put into two FASTA files:
        1. Numeric index headers, used for downstream phylogenetic placement
        2. Contains the original headers along with marker, start, and end positions included
    The negative integers (or numeric indexes) are stored in a dictionary and returned
    Sequences are grouped based on the location on the HMM profile they mapped to

    :param args: Command-line argument object from get_options and check_parser_arguments
    :param hmm_matches: Contains lists of HmmMatch objects mapped to the marker they matched
    :param fasta_dict: Stores either the original or ORF-predicted input FASTA. Headers are keys, sequences are values
    :return: List of files that go on to placement stage, dictionary mapping marker-specific numbers to contig names
    """
    logging.info("Extracting the quality-controlled protein sequences... ")
    hmmalign_input_fastas = list()
    marker_gene_dict = dict()
    numeric_contig_index = dict()
    trimmed_query_bins = dict()
    bins = dict()

    for marker in hmm_matches:
        if len(hmm_matches[marker]) == 0:
            continue
        if marker not in numeric_contig_index.keys():
            numeric_contig_index[marker] = dict()
        numeric_decrementor = -1
        if marker not in marker_gene_dict:
            marker_gene_dict[marker] = dict()

        # Algorithm for binning sequences:
        # 1. Sort HmmMatches by the proportion of the HMM profile they covered in increasing order (full-length last)
        # 2. For HmmMatch in sorted matches, determine overlap between HmmMatch and each bin's representative HmmMatch
        # 3. If overlap exceeds 80% of representative's aligned length add it to the bin, else continue
        # 4. When bins are exhausted create new bin with HmmMatch
        for hmm_match in sorted(hmm_matches[marker], key=lambda x: x.end - x.start):
            if hmm_match.desc != '-':
                contig_name = hmm_match.orf + '_' + hmm_match.desc
            else:
                contig_name = hmm_match.orf
            # Add the query sequence to the index map
            orf_coordinates = str(hmm_match.start) + '_' + str(hmm_match.end)
            numeric_contig_index[marker][numeric_decrementor] = contig_name + '_' + orf_coordinates
            # Add the FASTA record of the trimmed sequence - this one moves on for placement
            full_sequence = fasta_dict[reformat_string('>' + contig_name)]
            binned = False
            for bin_num in sorted(bins):
                bin_rep = bins[bin_num][0]
                overlap = min(hmm_match.pend, bin_rep.pend) - max(hmm_match.pstart, bin_rep.pstart)
                if (100*overlap)/(bin_rep.pend - bin_rep.pstart) > 80:
                    bins[bin_num].append(hmm_match)
                    trimmed_query_bins[bin_num] += '>' + str(numeric_decrementor) + "\n" + \
                                                   full_sequence[hmm_match.start - 1:hmm_match.end] + "\n"
                    binned = True
                    break
            if not binned:
                bin_num = len(bins)
                bins[bin_num] = list()
                bins[bin_num].append(hmm_match)
                trimmed_query_bins[bin_num] = '>' + str(numeric_decrementor) + "\n" + \
                                              full_sequence[hmm_match.start - 1:hmm_match.end] + "\n"

            # Now for the header format to be used in the bulk FASTA:
            # >contig_name|marker_gene|start_end
            bulk_header = '>' + contig_name + '|' +\
                          hmm_match.target_hmm + '|' +\
                          orf_coordinates
            marker_gene_dict[marker][bulk_header] = full_sequence[hmm_match.start-1:hmm_match.end]
            numeric_decrementor -= 1

        # Write all the homologs to the FASTA file
        for group in trimmed_query_bins:
            if trimmed_query_bins[group]:
                marker_query_fa = args.output_dir_var + marker + "_hmm_purified_group" + str(group) + ".faa"
                try:
                    homolog_seq_fasta = open(marker_query_fa, 'w')
                except IOError:
                    logging.error("Unable to open " + marker_query_fa + " for writing.\n")
                    sys.exit(3)
                hmmalign_input_fastas.append(marker_query_fa)
                homolog_seq_fasta.write(trimmed_query_bins[group])
                homolog_seq_fasta.close()
        trimmed_query_bins.clear()
        bins.clear()
    logging.info("done.\n")

    # Now write a single FASTA file with all identified markers
    for marker in marker_gene_dict:
        trimmed_hits_fasta = args.output_dir_final + marker + "_hmm_purified.faa"
        logging.debug("\tWriting " + marker + " sequences to " + trimmed_hits_fasta + "\n")
        write_new_fasta(marker_gene_dict[marker], trimmed_hits_fasta)
    return hmmalign_input_fastas, numeric_contig_index

 
def collect_blast_outputs(args):
    """
    Deletes empty BLAST results files.

    :param args: Command-line argument object from get_options and check_parser_arguments
    Returns a list of non-empty BLAST results files.
    """
    cog_blast_result = args.output_dir_var + path.basename(args.fasta_input) + "_formatted.BLAST_results_raw.txt"
    rrna_blast_result = args.output_dir_var + path.basename(args.fasta_input) + "_formatted.rRNA_BLAST_results_raw.txt"

    blast_tables = [cog_blast_result, rrna_blast_result]
    total_size = 0
    for blast_out in blast_tables:
        if os.path.exists(blast_out):
            total_size += path.getsize(blast_out)
            if path.getsize(blast_out) <= 0:
                os.remove(blast_out)
                blast_tables.pop()
        else:
            blast_tables.pop()

    if total_size <= 0:
        logging.error("No marker genes detected in input! Exiting...\n")
        sys.exit(3)

    return blast_tables


def subsequence(fasta_dictionary, contig_name, start, end):
    """
    Extracts a sub-sequence from `start` to `end` of `contig_name` in `fasta_dictionary`
     with headers for keys and sequences as values. `contig_name` does not contain the '>' character
    :param fasta_dictionary:
    :param contig_name:
    :param start:
    :param end:
    :return: A string representing the sub-sequence of interest
    """
    subseq = fasta_dictionary['>' + contig_name][start:end]
    return subseq


def write_nuc_sequences(args, gene_coordinates, formatted_fasta_dict):
    """
    Function to write the nucleotide sequences representing the BLAST alignment region for each hit in the fasta
    :param args: Command-line argument object from get_options and check_parser_arguments
    :param gene_coordinates:
    :param formatted_fasta_dict:
    :return: None
    """
    # Header format:
    # >contig_name|marker_gene|start_end
    # input_multi_fasta = re.match(r'\A.*\/(.*)', args.fasta_input).group(1)
    input_multi_fasta = path.basename(args.fasta_input)
    orf_nuc_fasta = args.output_dir_var + '.'.join(input_multi_fasta.split('.')[:-1]) + "_genes.fna"
    try:
        fna_output = open(orf_nuc_fasta, 'w')
    except IOError:
        raise IOError("Unable to open " + orf_nuc_fasta + " for writing!")

    output_fasta_string = ""

    for contig_name in gene_coordinates:
        start = 0
        end = 0
        for coords_start in sorted(gene_coordinates[contig_name].keys()):
            start = coords_start
            for coords_end in gene_coordinates[contig_name][coords_start].keys():
                end = coords_end
                cog = gene_coordinates[contig_name][coords_start][coords_end]
                output_fasta_string += '>' + contig_name + '|' + cog + '|' + str(start) + '_' + str(end) + "\n"
                output_fasta_string += subsequence(formatted_fasta_dict, contig_name, start, end) + "\n"

    fna_output.write(output_fasta_string)
    fna_output.close()

    return


def fprintf(opened_file, fmt, *args):
    """
    A helper function used to print to a specified file.
    :param opened_file: A file object that has already been opened using open()
    """
    opened_file.write(fmt % args)


def parse_infernal_table(infernal_table):
    """
    Function to parse the `--tblout` file from Infernal, generating a list of sequence names
    that meet pre-defined thresholds and their start and end positions on their contig
    :return: A dictionary containing contig names, and the start and end positions of an rRNA sequence (tuple)
    """

    rrna_seqs = dict()

    return rrna_seqs


def extract_rrna_sequences(rrna_seqs, rrna_marker, fasta_dictionary):
    """
    TEMPLATE
    :param rrna_seqs:
    :param rrna_marker:
    :return:
    """
    rrna_fasta_dict = dict()
    for contig_name in rrna_seqs:
        for coordinates in sorted(rrna_seqs[contig_name].keys()):
            start, end = coordinates
            rrna_fasta_dict['>' + contig_name + '|' + str(start) + '_' + str(end)] = subsequence(fasta_dictionary,
                                                                                                 contig_name,
                                                                                                 start,
                                                                                                 end)
    return rrna_fasta_dict


def detect_ribrna_sequences(args, cog_list, formatted_fasta_dict):
    """

    :param args: Command-line argument object from get_options and check_parser_arguments
    :param cog_list:
    :param formatted_fasta_dict:
    :return:
    """
    logging.info("Retrieving rRNA hits with Infernal... ")

    function_start_time = time.time()

    num_rrna = 0

    for rrna_marker in cog_list["phylogenetic_cogs"]:
        cov_model = ""
        infernal_table = ""

        cmsearch_command = [args.executables["cmsearch"]]
        cmsearch_command += ["--tblout", infernal_table]
        cmsearch_command += ["--noali", "--cpu", str(args.num_threads)]
        cmsearch_command += [cov_model, args.input]

        rrna_seqs = parse_infernal_table(infernal_table)
        num_rrna += len(rrna_seqs.values())
        extract_rrna_sequences(rrna_seqs, rrna_marker, formatted_fasta_dict)

    logging.info("done.\n")

    function_end_time = time.time()
    hours, remainder = divmod(function_end_time - function_start_time, 3600)
    minutes, seconds = divmod(remainder, 60)
    logging.debug("\trRNA-identification time required: " +
                  ':'.join([str(hours), str(minutes), str(round(seconds, 2))]) + "\n")
    logging.debug("\t" + str(num_rrna) + " rRNA sequences found.\n\n")

    return


def get_ribrna_hit_sequences(args, blast_hits_purified, genewise_summary_files, cog_list):
    if args.py_version == 3:
        izip = zip
    else:
        from itertools import izip
    # TODO: It doesn't need to return anything; OR this function could return the contig_rrna_coordinates instead of writing files
    """
    rRNA does not get translated into protein. Regardless, we want to take the
    rRNA and summarize it in a way that is parallel to the Genewise summary files.
    This function does that using the provided Autovivification of purified BLAST
    hits, list of COGs, and Autovivification of Genewise summary files.

    Returns an Autovivification summarizing the coordinates for each rRNA hit.
    Returns a list of the rRNA summary files.
    """

    logging.info("Retrieving rRNA hits... ")

    contig_rrna_coordinates = Autovivify()
    rRNA_hit_files = {}
    rrna_seqs = 0

    function_start_time = time.time()
    
    for contig in sorted(blast_hits_purified.keys()):
        # note: We skipped the Genewise step (we are dealing with rRNA) but we bring the rRNA files in the
        # same structure as the Genewise summary files and bring them back into the ordinary pipeline.
        for identifier in sorted(blast_hits_purified[contig].keys()):
            # if not re.search("rRNA", blast_hits_purified[contig][identifier]['cog']):
            if blast_hits_purified[contig][identifier]['cog'] not in cog_list["phylogenetic_rRNA_cogs"]:
                continue

            start = blast_hits_purified[contig][identifier]["start"]
            end = blast_hits_purified[contig][identifier]["end"]
            cog = blast_hits_purified[contig][identifier]["cog"]
            direction = blast_hits_purified[contig][identifier]["direction"]
            contig_rrna_coordinates[contig][identifier]["start"] = start
            contig_rrna_coordinates[contig][identifier]["end"] = end
            contig_rrna_coordinates[contig][identifier]["cog"] = cog
            contig_rrna_coordinates[contig][identifier]["direction"] = direction
            outfile_name = args.output_dir_var + contig + '_rRNA_result_summary.txt'
            contig_rrna_coordinates[contig][identifier]["outfile"] = outfile_name
            genewise_summary_files[contig][outfile_name] = 1
            try:
                outfile = open(outfile_name, 'w')
                outfile.close()
            except IOError:
                logging.error("Unable to create create " + outfile_name + ' for writing!\n')
                sys.exit(3)

    # This overwrites the original log file
    fasta_list = _fasta_reader._read_format_fasta(args.fasta_input, args.min_seq_length, args.output, 'dna', 110)
    if not fasta_list:
        sys.exit()
    tmp_iterable = iter(fasta_list)
    formatted_fasta_dict = dict(izip(tmp_iterable, tmp_iterable))

    for contig_name in formatted_fasta_dict:
        sequence = formatted_fasta_dict[contig_name]
        contig_name = contig_name[1:]
        if contig_name in contig_rrna_coordinates:
            # start searching for the information to shorten the file.
            for identifier in sorted(contig_rrna_coordinates[contig_name].keys()):
                start = contig_rrna_coordinates[contig_name][identifier]["start"]
                end = contig_rrna_coordinates[contig_name][identifier]["end"]
                cog = contig_rrna_coordinates[contig_name][identifier]["cog"]
                direction = contig_rrna_coordinates[contig_name][identifier]["direction"]
                outfile = contig_rrna_coordinates[contig_name][identifier]['outfile']
                count = -1
                shortened_sequence = ""
                for nucleotide in sequence:
                    count += 1
                    if not (start <= count <= end):
                        continue
                    shortened_sequence += nucleotide

                if direction == 'reverse':
                    # ok, our hit has been on the opposite strand of the reference.
                    # to get a proper alignment, we thus have to produce a negative strand version of the input
                    nucleotides2 = ''.join(reversed(shortened_sequence))
                    shortened_sequence = ""
                    nucleotides2 = nucleotides2.lower()
                    for nucleotide in nucleotides2:
                        if nucleotide == 't':
                            nucleotide = 'a'
                        elif nucleotide == 'a':
                            nucleotide = 't'
                        elif nucleotide == 'c':
                            nucleotide = 'g'
                        elif nucleotide == 'g':
                            nucleotide = 'c'

                        shortened_sequence += nucleotide
                rrna_seqs += 1
                try:
                    out = open(outfile, 'a')
                    fprintf(out, '%s\t%s\t%s\t%s\t%s\n', cog, start, end, 'n/a', shortened_sequence)
                    out.close()
                except IOError:
                    sys.stderr.write("ERROR: Can't create " + outfile + '!\n')
                    sys.exit(0)

            try:
                output_file = open(args.output_dir_var + contig_name + '_sequence.txt', 'w')
                fprintf(output_file, '>%s\n%s', contig_name, sequence)
                output_file.close()
            except IOError:
                sys.stderr.write("ERROR: Can't create " + args.output_dir_var + contig_name + '_sequence.txt!\n')
                sys.exit(12)

    logging.info("done.\n")

    function_end_time = time.time()
    hours, remainder = divmod(function_end_time - function_start_time, 3600)
    minutes, seconds = divmod(remainder, 60)
    logging.debug("\trRNA-identification time required: " +
                  ':'.join([str(hours), str(minutes), str(round(seconds, 2))]) + "\n")
    logging.debug("\t" + str(rrna_seqs) + " rRNA sequences found.\n\n")

    return contig_rrna_coordinates, rRNA_hit_files


def get_sequence_counts(concatenated_mfa_files, ref_alignment_dimensions, verbosity, file_type):
    alignment_length_dict = dict()
    for denominator in concatenated_mfa_files:
        if denominator not in ref_alignment_dimensions:
            logging.error("Unrecognized code '" + denominator + "'.")
            sys.exit(3)

        ref_n_seqs, ref_seq_length = ref_alignment_dimensions[denominator]
        for msa_file in concatenated_mfa_files[denominator]:
            if file_type == "Fasta":
                seq_dict = read_fasta_to_dict(msa_file)
            elif file_type == "Phylip":
                seq_dict = read_phylip_to_dict(msa_file)
            elif file_type == "Stockholm":
                seq_dict = read_stockholm_to_dict(msa_file)
            else:
                logging.error("File type '" + file_type + "' is not recognized.")
                sys.exit(3)
            num_seqs, sequence_length = validate_multiple_alignment_utility(seq_dict, msa_file)
            alignment_length_dict[msa_file] = sequence_length

            # Warn user if the multiple sequence alignment has grown significantly
            if verbosity and ref_seq_length*1.5 < sequence_length:
                logging.warning("Multiple alignment of " + denominator +
                                "caused >150% increase in the number of columns (" +
                                str(ref_seq_length) + '->' + str(sequence_length) + ").\n")
    return alignment_length_dict


def validate_multiple_alignment_utility(seq_dict, mfa_file):
    """
    Checks to ensure all sequences are the same length and returns a tuple of (nrow, ncolumn)

    :param seq_dict: A dictionary containing headers as keys and sequences as values
    :param mfa_file: The name of the multiple alignment FASTA file being validated
    :return: tuple = (nrow, ncolumn)
    """
    sequence_length = 0
    for seq_name in seq_dict:
        sequence = seq_dict[seq_name]
        if sequence_length == 0:
            sequence_length = len(sequence)
        elif sequence_length != len(sequence) and sequence_length > 0:
            logging.error("Number of aligned columns is inconsistent in " + mfa_file + "!\n")
            sys.exit(3)
        else:
            pass
            # Sequence is the right length, carrying on
    return len(seq_dict), sequence_length


def get_alignment_dims(args, marker_build_dict):
    alignment_dimensions_dict = dict()
    alignment_data_dir = os.sep.join([args.treesapp, 'data', 'alignment_data' + os.sep])
    try:
        os.path.isdir(alignment_data_dir)
    except IOError:
        logging.error(alignment_data_dir + "does not exist.")
        sys.exit(3)

    all_markers = [marker_build_dict[marker_code].cog for marker_code in marker_build_dict]

    for fasta in glob.glob(alignment_data_dir + "*fa"):
        cog = os.path.basename(fasta).split('.')[0]
        if cog in all_markers:
            for marker_code in marker_build_dict:
                if marker_build_dict[marker_code].cog == cog:
                    seq_dict = read_fasta_to_dict(fasta)
                    alignment_dimensions_dict[marker_code] = (validate_multiple_alignment_utility(seq_dict, fasta))
    return alignment_dimensions_dict


def multiple_alignments(args, single_query_sequence_files, marker_build_dict, tool="hmmalign"):
    """
    Wrapper function for the multiple alignment functions - only purpose is to make an easy decision at this point...

    :param args: Command-line argument object from get_options and check_parser_arguments
    :param single_query_sequence_files:
    :param marker_build_dict:
    :param tool: Tool to use for aligning query sequences to a reference multiple alignment [hmmalign|papara]
    :return: list of multiple sequence alignment files, generated by `tool`
    """
    if tool == "papara":
        singlehit_files = prepare_and_run_papara(args, single_query_sequence_files, marker_build_dict)
    elif tool == "hmmalign":
        singlehit_files = prepare_and_run_hmmalign(args, single_query_sequence_files, marker_build_dict)
    else:
        logging.error("Unrecognized tool '" + str(tool) + "' for multiple sequence alignment.\n")
        sys.exit(3)
    return singlehit_files


def create_ref_phy_files(args, single_query_fasta_files, marker_build_dict, ref_alignment_dimensions):
    """
    Creates a phy file for every reference marker that was matched by a query sequence
    :param args:
    :param single_query_fasta_files:
    :param marker_build_dict:
    :param ref_alignment_dimensions:
    :return:
    """
    treesapp_resources = args.treesapp + os.sep + 'data' + os.sep

    # Convert the reference sequence alignments to .phy files for every marker identified
    for query_fasta in single_query_fasta_files:
        marker = re.match("(.*)_hmm_purified.*", os.path.basename(query_fasta)).group(1)
        denominator = None
        for denominator in marker_build_dict.keys():
            if marker_build_dict[denominator].cog == marker:
                break

        ref_alignment_phy = args.output_dir_var + marker + ".phy"
        if os.path.isfile(ref_alignment_phy):
            continue
        aligned_fasta = treesapp_resources + "alignment_data" + os.sep + marker + ".fa"

        num_ref_seqs, ref_align_len = ref_alignment_dimensions[denominator]
        aligned_fasta_dict = read_fasta_to_dict(aligned_fasta)
        dict_for_phy = dict()
        for seq_name in aligned_fasta_dict:
            dict_for_phy[seq_name.split('_')[0]] = aligned_fasta_dict[seq_name]
        phy_dict = reformat_fasta_to_phy(dict_for_phy)

        write_phy_file(ref_alignment_phy, phy_dict, (num_ref_seqs, ref_align_len))
    return


def prepare_and_run_papara(args, single_query_fasta_files, marker_build_dict):
    """
    Uses the Parsimony-based Phylogeny-aware short Read Alignment (PaPaRa) tool

    :param args:
    :param single_query_fasta_files:
    :param marker_build_dict:
    :return:
    """
    treesapp_resources = args.treesapp + os.sep + 'data' + os.sep
    query_alignment_files = dict()
    logging.info("Running PaPaRa... ")
    # TODO: Parallelize; PaPaRa's posix threading is not guaranteed with binary
    start_time = time.time()

    # Convert the reference sequence alignments to .phy files for every marker identified
    for query_fasta in sorted(single_query_fasta_files):
        file_name_info = re.match("(.*)_hmm_purified.*\.(f.*)$", os.path.basename(query_fasta))
        if file_name_info:
            marker, extension = file_name_info.groups()
        else:
            logging.error("Unable to parse information from file name:" + "\n" + str(query_fasta) + "\n")
            sys.exit(3)

        ref_marker = None
        for denominator in marker_build_dict:
            if marker == marker_build_dict[denominator].cog:
                ref_marker = marker_build_dict[denominator]
                break
        query_multiple_alignment = re.sub('.' + re.escape(extension) + r"$", ".phy", query_fasta)
        tree_file = treesapp_resources + "tree_data" + os.sep + marker + "_tree.txt"
        ref_alignment_phy = args.output_dir_var + marker + ".phy"
        if not os.path.isfile(ref_alignment_phy):
            logging.error("Phylip file '" + ref_alignment_phy + "' not found.\n")
            sys.exit(3)

        run_papara(args.executables["papara"], tree_file, ref_alignment_phy, query_fasta, ref_marker.molecule)
        shutil.copy("papara_alignment.default", query_multiple_alignment)
        os.remove("papara_alignment.default")
        if ref_marker.denominator not in query_alignment_files:
            query_alignment_files[ref_marker.denominator] = []
        query_alignment_files[ref_marker.denominator].append(query_multiple_alignment)

    logging.info("done.\n")
    os.remove("papara_log.default")
    os.remove("papara_quality.default")

    end_time = time.time()
    hours, remainder = divmod(end_time - start_time, 3600)
    minutes, seconds = divmod(remainder, 60)
    logging.debug("\tPaPaRa time required: " +
                  ':'.join([str(hours), str(minutes), str(round(seconds, 2))]) + "\n")

    return query_alignment_files


def prepare_and_run_hmmalign(args, single_query_fasta_files, marker_build_dict):
    """
    Runs `hmmalign` to add the query sequences into the reference FASTA multiple alignments

    :param args:
    :param single_query_fasta_files:
    :param marker_build_dict:
    :return: list of multiple sequence alignment files (in Clustal format) generated by hmmalign.
    """

    hmmalign_singlehit_files = dict()
    logging.info("Running hmmalign... ")

    start_time = time.time()
    # task_list = list()

    # Run hmmalign on each fasta file
    for query_fasta in sorted(single_query_fasta_files):
        file_name_info = re.match("(.*)_hmm_purified.*\.(f.*)$", os.path.basename(query_fasta))
        if file_name_info:
            marker, extension = file_name_info.groups()
        else:
            logging.error("Unable to parse information from file name:" + "\n" + str(query_fasta) + "\n")
            sys.exit(3)

        query_multiple_alignment = re.sub('.' + re.escape(extension) + r"$", ".sto", query_fasta)

        ref_marker = None
        for denominator in marker_build_dict:
            if marker == marker_build_dict[denominator].cog:
                ref_marker = marker_build_dict[denominator]
                break
        if not ref_marker:
            logging.error("Unable to match marker '" + marker + "' to a code.\n")
            sys.exit(3)

        # Get the paths to either the HMM or CM profile files
        ref_alignment = os.sep.join([args.treesapp, 'data', "alignment_data", ref_marker.cog + ".fa"])
        ref_profile = os.sep.join([args.treesapp, 'data', "hmm_data", ref_marker.cog])
        if ref_marker.kind == "phylogenetic_rRNA":
            ref_profile += ".cm"
        else:
            ref_profile += ".hmm"
        profile_aligner(args.executables, ref_alignment, ref_profile,
                        query_fasta, query_multiple_alignment, ref_marker.kind)

        if ref_marker.denominator not in hmmalign_singlehit_files:
            hmmalign_singlehit_files[ref_marker.denominator] = []
        mfa_file = re.sub("\.sto$", ".mfa", query_multiple_alignment)
        tmp_dict = read_stockholm_to_dict(query_multiple_alignment)
        seq_dict = dict()
        for seq_name in tmp_dict:
            seq_dict[seq_name.split('_')[0]] = tmp_dict[seq_name]
        write_new_fasta(seq_dict, mfa_file)
        hmmalign_singlehit_files[ref_marker.denominator].append(mfa_file)

    # num_tasks = len(task_list)
    # if num_tasks > 0:
    #     cl_farmer = CommandLineFarmer("cmalign/hmmalign --mapali", args.num_threads)
    #     cl_farmer.add_tasks_to_queue(task_list)
    #
    #     cl_farmer.task_queue.close()
    #     cl_farmer.task_queue.join()

    logging.info("done.\n")

    end_time = time.time()
    hours, remainder = divmod(end_time - start_time, 3600)
    minutes, seconds = divmod(remainder, 60)
    logging.debug("\thmmalign time required: " +
                  ':'.join([str(hours), str(minutes), str(round(seconds, 2))]) + "\n")

    return hmmalign_singlehit_files


def cat_hmmalign_singlehit_files(args, hmmalign_singlehit_files):
    """
    The hmmalign command write a Stockholm file (because an MFA is not an option) and therefore we need to convert it
    to a FASTA file. To do so, TreeSAPP concatenates the individual lines (~80 characters long) for each sequence.
    :param args: Command-line argument object from get_options and check_parser_arguments
    :param hmmalign_singlehit_files:
    Returns a list of the files containing the concatenated hmmalign results.
    Returns a list of the number of sequences found in each file.
    """

    # For each type of gene...
    concatenated_mfa_files = {}
    nrs_of_sequences = {}

    logging.info("Reformatting hmmalign output files to FASTA... ")

    for clustal_mfa_file in sorted(hmmalign_singlehit_files):
        # Determine what type of gene is currently represented, or raise an error
        file_name_info = re.match("^([A-Z][0-9]{4}_.*)_\d+_\d+_([A-Za-z0-9_]+).cl$",
                                  os.path.basename(clustal_mfa_file))
        if file_name_info:
            f_contig, cog = file_name_info.groups()
        else:
            logging.error("Regular expression unable to pull contig and marker information from file name!\n" +
                          "Offending file:\n\t" + clustal_mfa_file + "\n")
            sys.exit(3)
        sequences = dict()
        query_sequence = ""
        parsing_order = dict()
        cog_rep_sequences = dict()
        acc = 0
        if cog not in cog_rep_sequences.keys():
            acc += 1
        cog_rep_sequences[cog] = set()

        if f_contig not in concatenated_mfa_files:
            concatenated_mfa_files[f_contig] = list()
        cog_len = 0

        # Begin reading the file
        try:
            hmmalign_msa = open(clustal_mfa_file, 'r')
        except IOError:
            logging.error("Can't open " + clustal_mfa_file + " for reading!\n")
            sys.exit(3)

        # Get sequence from file
        for _line in hmmalign_msa:
            line = _line.strip()
            reached_data_part = re.match(r'\A(.+) (\S+)\Z', line)
            if not reached_data_part:
                continue
            search_result = re.search(r'\A(.+) (\S+)\Z', line)
            if search_result:
                name_long = search_result.group(1)
                sequence_part = search_result.group(2)
                sequence_name = ''
                if re.search(r'query', name_long):
                    query_sequence += sequence_part
                    cog_len += len(sequence_part)

                elif re.search(r'(\d+)_', name_long):
                    sequence_name = re.search(r'(\d+)_', name_long).group(1)
                    cog_rep_sequences[cog].add(sequence_name)
                    if sequence_name not in sequences.keys():
                        sequences[sequence_name] = dict()
                    if cog not in sequences[sequence_name].keys():
                        sequences[sequence_name][cog] = ""
                    sequences[sequence_name][cog] += sequence_part

        parsing_order[acc] = cog, cog_len
        hmmalign_msa.close()

        concatenated_mfa_files[f_contig].append(args.output_dir_var + f_contig + '.mfa')
        # Write to the output file
        try:
            output = open(args.output_dir_var + f_contig + '.mfa', 'w')
        except IOError:
            logging.error("Can't create " + args.output_dir_var + f_contig + '.mfa\n')
            sys.exit(3)
        output.write('>query\n' + query_sequence + '\n')
        nrs_of_sequences[f_contig] = 1
        qlen = len(query_sequence)

        for sequence_name in sequences.keys():
            nrs_of_sequences[f_contig] += 1
            sequence = ""
            for p_order in sorted(parsing_order.keys(), key=int):
                cog, cog_len = parsing_order[p_order]
                # print f_contig, sequence_name, p_order, cog
                if sequence_name not in cog_rep_sequences[cog]:
                    sequence += "." * cog_len
                else:
                    sequence += sequences[sequence_name][cog]
            output.write('>' + sequence_name + '\n' + sequence + '\n')
            if len(sequence) != qlen:
                output.close()
                logging.error("Inconsistent sequence lengths between query and concatenated HMM alignments!\n" +
                              "Check " + args.output_dir_var + f_contig + ".mfa for bad sequence " + sequence_name)
                sys.exit(3)

        output.close()

    logging.info("done.\n")

    return concatenated_mfa_files, nrs_of_sequences


def filter_multiple_alignments(args, concatenated_mfa_files, marker_build_dict, tool="BMGE"):
    """
    Runs BMGE using the provided lists of the concatenated hmmalign files, and the number of sequences in each file.

    :param args:
    :param concatenated_mfa_files: A dictionary containing f_contig keys mapping to a FASTA or Phylip sequential file
    :param marker_build_dict:
    :param tool:
    :return: A list of files resulting from BMGE multiple sequence alignment masking.
    """
    # TODO: Parallelize with multiprocessing

    logging.info("Running " + tool + "... ")

    start_time = time.time()

    trimmed_output_files = {}

    for denominator in sorted(concatenated_mfa_files.keys()):
        if denominator not in trimmed_output_files:
            trimmed_output_files[denominator] = []
        mfa_files = concatenated_mfa_files[denominator]
        for concatenated_mfa_file in mfa_files:
            trimmed_msa_file = trim_multiple_alignment(args.executables["BMGE.jar"], concatenated_mfa_file,
                                                       marker_build_dict[denominator].molecule, tool)
            trimmed_output_files[denominator].append(trimmed_msa_file)

    logging.info("done.\n")

    end_time = time.time()
    hours, remainder = divmod(end_time - start_time, 3600)
    minutes, seconds = divmod(remainder, 60)
    logging.debug("\t" + tool + " time required: " +
                  ':'.join([str(hours), str(minutes), str(round(seconds, 2))]) + "\n")
    return trimmed_output_files


def check_for_removed_sequences(args, mfa_files: dict, marker_build_dict: dict):
    """
    Reads the multiple alignment files (either Phylip or FASTA formatted) and looks for both reference and query
    sequences that have been removed. Multiple alignment files are removed from `mfa_files` if:
        1. all query sequences were removed; a DEBUG message is issued
        2. at least one reference sequence was removed
    This quality-control function is necessary for placing short query sequences onto reference trees.

    :param args:
    :param mfa_files:
    :param marker_build_dict:
    :return: A dictionary of denominators, with multiple alignment dictionaries as values. Example:
        {M0702: { "McrB_hmm_purified.phy-BMGE.fasta": {'1': seq1, '2': seq2}}}
    """
    qc_ma_dict = dict()
    num_successful_alignments = 0
    discarded_seqs_string = ""
    logging.debug("Validating trimmed multiple sequence alignment files... ")

    for denominator in sorted(mfa_files.keys()):
        marker = marker_build_dict[denominator].cog
        # Create a set of the reference sequence names
        ref_headers = get_headers(os.sep.join([args.treesapp, "data", "alignment_data", marker + ".fa"]))
        unique_refs = set([re.sub('_' + re.escape(marker), '', x)[1:] for x in ref_headers])
        msa_passed, summary_str = validate_alignment_trimming(mfa_files[denominator], unique_refs, True, args.min_seq_length)
        num_successful_alignments += len(msa_passed)
        qc_ma_dict[denominator] = msa_passed
        discarded_seqs_string += summary_str

    logging.debug("done.\n")
    logging.debug("\tSequences <" + str(args.min_seq_length) + " characters removed:" + discarded_seqs_string + "\n")

    if num_successful_alignments == 0:
        logging.error("No quality alignment files to analyze after trimming. Exiting now.\n")
        sys.exit(0)  # Should be 3, but this allows Clade_exclusion_analyzer to continue after exit

    return qc_ma_dict


def evaluate_trimming_performance(qc_ma_dict, alignment_length_dict, concatenated_msa_files, tool):
    """

    :param qc_ma_dict: A dictionary mapping denominators to files to multiple alignment dictionaries
    :param alignment_length_dict:
    :param concatenated_msa_files: Dictionary with markers indexing original (untrimmed) multiple alignment files
    :param tool: The name of the tool that was appended to the original, untrimmed or unmasked alignment files
    :return: None
    """
    trimmed_length_dict = dict()
    for denominator in sorted(qc_ma_dict.keys()):
        if len(concatenated_msa_files[denominator]) >= 1:
            of_ext = concatenated_msa_files[denominator][0].split('.')[-1]
        else:
            continue
        if denominator not in trimmed_length_dict:
            trimmed_length_dict[denominator] = list()
        for multi_align_file in qc_ma_dict[denominator]:
            file_type = multi_align_file.split('.')[-1]
            multi_align = qc_ma_dict[denominator][multi_align_file]
            num_seqs, trimmed_seq_length = validate_multiple_alignment_utility(multi_align, multi_align_file)

            original_multi_align = re.sub('-' + tool + '.' + file_type, '.' + of_ext, multi_align_file)
            raw_align_len = alignment_length_dict[original_multi_align]
            diff = raw_align_len - trimmed_seq_length
            if diff < 0:
                logging.warning("MSA length increased after " + tool + " processing for " + multi_align_file + "\n")
            else:
                trimmed_length_dict[denominator].append(diff)

    trimming_performance_string = "\tAverage columns removed:\n"
    for denominator in trimmed_length_dict:
        trimming_performance_string += "\t\t" + denominator + "\t"
        trimming_performance_string += str(round(sum(trimmed_length_dict[denominator])/len(trimmed_length_dict[denominator]), 1)) + "\n"

    logging.debug(trimming_performance_string + "\n")
    return


def produce_phy_files(args, qc_ma_dict):
    """
    Produces phy files from the provided list of alignment files

    :param args:
    :param qc_ma_dict:
    :return: Dictionary containing the names of the produced phy files mapped to its f_contig
    """

    phy_files = dict()
    sequence_lengths = dict()

    logging.debug("Writing filtered multiple alignment files to Phylip... ")

    # Open each alignment file
    for denominator in sorted(qc_ma_dict.keys()):
        sequence_lengths[denominator] = set()
        # Prepare the phy file for writing
        if denominator not in phy_files.keys():
            phy_files[denominator] = list()

        for multi_align_file in qc_ma_dict[denominator]:
            multi_align_dict = qc_ma_dict[denominator][multi_align_file]
            final_phy_file_name = re.sub(".fasta$|.phy$", "-qcd.phy", multi_align_file)
            sequences_for_phy = dict()

            for name in sorted(multi_align_dict.keys()):
                seq_name = name.strip()
                seq_name = seq_name.split('_')[0]

                sequence = multi_align_dict[name]
                sequence = re.sub(r' ', '', sequence)
                sequence_lengths[denominator].add(len(sequence))
                # Ensure the sequences contain only valid characters for RAxML
                sequence = re.sub(r'\.', 'X', sequence)
                sequence = re.sub(r'\*', 'X', sequence)
                sequence = re.sub('-', 'X', sequence)

                if re.search(r'\AX+\Z', sequence):
                    sequence = re.sub('X', 'V', sequence, 1)

                if args.molecule != "prot":
                    sequence = re.sub('U', 'T', sequence)  # Got error from RAxML when encountering Uracil

                sequences_for_phy[seq_name] = sequence

            # Write the sequences to the phy file
            phy_dict = reformat_fasta_to_phy(sequences_for_phy)
            if len(sequence_lengths[denominator]) != 1:
                logging.error("Sequence lengths varied in " + multi_align_file + "\n")
            phy_string = ' ' + str(len(multi_align_dict.keys())) + '  ' + str(sequence_lengths[denominator].pop()) + '\n'
            for count in sorted(phy_dict.keys(), key=int):
                for seq_name in sorted(phy_dict[count].keys()):
                    sequence_part = phy_dict[count][seq_name]
                    if count == 0:
                        phy_string += str(seq_name)
                        length = len(str(seq_name))
                        c = length
                        while c < 10:
                            phy_string += ' '
                            c += 1
                    phy_string += sequence_part + '\n'

                phy_string += '\n'

            with open(final_phy_file_name, 'w') as phy_output:
                phy_output.write(phy_string)
            phy_files[denominator].append(final_phy_file_name)
    logging.debug("done.\n")

    return phy_files


def pparse_ref_trees(denominator_ref_tree_dict, args):
    ref_trees_dict = dict()

    pool = Pool(processes=int(args.num_threads))

    def log_tree(result):
        marker, terminal_children_of_reference = result
        if terminal_children_of_reference is None:
            logging.warning("Letting threads finish before exiting... ")
        ref_trees_dict[marker] = terminal_children_of_reference

    for denominator in denominator_ref_tree_dict:
        reference_tree_file = denominator_ref_tree_dict[denominator]
        pool.apply_async(func=read_and_understand_the_reference_tree,
                         args=(reference_tree_file, denominator, ),
                         callback=log_tree)
    pool.close()
    pool.join()
    for marker in ref_trees_dict:
        if ref_trees_dict[marker] is None:
            logging.info("done.\n")
            return None
        else:
            pass
    return ref_trees_dict


def pparse_raxml_out_trees(labelled_trees, args):
    """
    The wrapper command for parsing all trees of a gene family (denominator) in parallel
    :param labelled_trees: Dictionary containing labelled tree files for each f_contig
    :param args: args object (for num_threads)
    :return: Dictionary containing all parsed trees for each contig
    """
    raxml_tree_dict = dict()

    pool = Pool(processes=int(args.num_threads))

    def log_tree(result):
        f_contig, rooted_labelled_trees, insertion_point_node_hash = result
        if rooted_labelled_trees is None:
            pool.terminate()
            sys.exit()
        raxml_tree_dict[f_contig] = [rooted_labelled_trees, insertion_point_node_hash]

    def no_tree_handler(error):
        logging.error(error + "-->{}\n<--".format(error.__cause__))
        pool.terminate()
        sys.exit(5)

    for f_contig in labelled_trees:
        tree_file = labelled_trees[f_contig]
        if args.py_version == 3:
            pool.apply_async(func=read_understand_and_reroot_the_labelled_tree,
                             args=(tree_file, f_contig, ),
                             callback=log_tree,
                             error_callback=no_tree_handler)
        else:
            pool.apply_async(func=read_understand_and_reroot_the_labelled_tree,
                             args=(tree_file, f_contig, ),
                             callback=log_tree)

    pool.close()
    pool.join()
    return raxml_tree_dict


def read_understand_and_reroot_the_labelled_tree(labelled_tree_file, f_contig):
    labelled_tree_elements, insertion_point_node_hash = read_the_raxml_out_tree(labelled_tree_file)
    # # Old and slow:
    # labelled_tree_info = create_tree_info_hash()
    # labelled_tree_info = get_node_subtrees(labelled_tree_elements, labelled_tree_info)
    # labelled_tree_info = assign_parents_and_children(labelled_tree_info, f_contig)
    # if labelled_tree_info is None:
    #     return [f_contig, None, insertion_point_node_hash]
    # labelled_tree_info = build_tree_info_quartets(labelled_tree_info)
    # rooted_labelled_trees = build_newly_rooted_trees(labelled_tree_info)
    # return [f_contig, rooted_labelled_trees, insertion_point_node_hash]

    # Using the C++ _tree_parser extension:
    labelled_tree_assignments = _tree_parser._get_parents_and_children(labelled_tree_elements)
    if labelled_tree_assignments == "$":
        sys.stderr.write("Poison pill received from " + f_contig + "\n")
        sys.stderr.flush()
        return [f_contig, None, insertion_point_node_hash]
    else:
        labelled_tree_info, terminal_children_of_labelled_tree = deconvolute_assignments(labelled_tree_assignments)
        labelled_tree_info['subtree_of_node'] = terminal_children_of_labelled_tree
        labelled_tree_info = build_tree_info_quartets(labelled_tree_info)
        rooted_labelled_trees = build_newly_rooted_trees(labelled_tree_info)
        return [f_contig, rooted_labelled_trees, insertion_point_node_hash]


def identify_the_correct_terminal_children_of_each_assignment(terminal_children_of_reference,
                                                              rooted_labelled_trees,
                                                              insertion_point_node_hash,
                                                              assignments, num_threads, parse_log):
    terminal_children_of_assignments = build_terminal_children_strings_of_assignments(rooted_labelled_trees,
                                                                                      insertion_point_node_hash,
                                                                                      assignments,
                                                                                      num_threads,
                                                                                      parse_log)
    real_terminal_children_of_assignments = compare_terminal_children_strings(terminal_children_of_assignments,
                                                                              terminal_children_of_reference,
                                                                              parse_log)
    return real_terminal_children_of_assignments


def read_the_raxml_out_tree(labelled_tree_file):
    """
    Reads and reformats the labelled_tree_file for downstream interpretation
    :param labelled_tree_file: RAxML output f_contig.originalRAxML_labelledTree.txt file in various_outputs directory
    :return: An easily interpretable labelled tree and a collection of
    """

    insertion_point_node_hash = Autovivify()
    try:
        raxml_tree = open(labelled_tree_file, 'r')
    except IOError:
        logging.error("Could not open " + labelled_tree_file + " for reading!\n")
        sys.exit(5)
    tree_string = ''

    for line in raxml_tree:
        line = line.strip()
        tree_string += line

    raxml_tree.close()
    tree_symbols_raw_1 = list(tree_string)
    bracket_diff = 0
    tree_string_neu = '('
    comma_count = 0

    for tree_symbol_raw_1 in tree_symbols_raw_1:
        if comma_count < 2:
            if tree_symbol_raw_1 == '(':
                bracket_diff += 1
            if tree_symbol_raw_1 == ')':
                bracket_diff -= 1
            if tree_symbol_raw_1 == ',' and bracket_diff == 1:
                comma_count += 1
            if comma_count == 2:
                tree_string_neu += '):1.0[I666999666]'
        tree_string_neu += tree_symbol_raw_1

    tree_string = tree_string_neu
    tree_string = re.sub('\(', 'L', tree_string)
    tree_string = re.sub('\)', 'R', tree_string)
    tree_string = re.sub('\[', 'Q', tree_string)

    # Remove the branch lengths
    try:
        tree_string = re.sub(":[.0-9]+Q", 'Q', tree_string)
    except AssertionError:
        raise AssertionError("Unable to remove branch lengths!")

    while re.search(r'((\D(\d+))QI(\d+)])', tree_string):
        to_be_replaced = re.search(r'((\D(\d+))QI(\d+)])', tree_string).group(1)
        replacement = re.search(r'((\D(\d+))QI(\d+)])', tree_string).group(2)
        terminal_leaf = re.search(r'((\D(\d+))QI(\d+)])', tree_string).group(3)
        insertion_point = re.search(r'((\D(\d+))QI(\d+)])', tree_string).group(4)
        if int(terminal_leaf) <= 0:
            sys.stderr.write("ERROR: Your tree has terminal leaves with numbers <= 0. "
                             "Please change them to positive values!\n")
            sys.stderr.flush()
            sys.exit(-1)
        insertion_point_node_hash[insertion_point] = terminal_leaf
        tree_string = re.sub(to_be_replaced, replacement, tree_string)
    count = -2
    while re.search(r'QI(\d+)]', tree_string):
        insertion_point_node_hash[re.search(r'QI(\d+)]', tree_string).group(1)] = count
        tree_string = re.sub(r'QI(\d+)]', str(count), tree_string, 1)
        count += -1

    tree_string = re.sub('L', '(', tree_string)
    tree_string = re.sub('R', ')', tree_string)
    tree_string = re.sub('Q', '[', tree_string)
    # Remove these lines when using the C++ extension:
    # tree_elements = split_tree_string(tree_string)
    # return tree_elements, insertion_point_node_hash
    return tree_string, insertion_point_node_hash


def split_tree_string(tree_string):
    tree_symbols_raw = list(str(tree_string))
    count = -1
    previous_symbol = ''
    tree_elements = Autovivify()
    
    for tree_symbol_raw in tree_symbols_raw:
        if re.search(r'\d', tree_symbol_raw) and (re.search(r'\d', previous_symbol) or previous_symbol == '-'):
            tree_elements[count] += tree_symbol_raw
        else:
            count += 1
            tree_elements[count] = tree_symbol_raw
        previous_symbol = tree_symbol_raw
    return tree_elements


def build_tree_info_quartets(tree_info):
    for node in sorted(tree_info['parent_of_node'].keys(), key=int):
        parent = tree_info['parent_of_node'][node]
        if int(parent) == -1:
            for roots_child in sorted(tree_info['children_of_node']['-1'].keys(), key=int):
                if roots_child == node:
                    continue
                parent = roots_child

        tree_info['quartets'][node][parent] = 1
        if node in tree_info['children_of_node']:
            for child in sorted(tree_info['children_of_node'][node].keys(), key=int):
                tree_info['quartets'][node][child] = 1

    return tree_info


def build_newly_rooted_trees(tree_info):
    """
    Builds a new tree that is re-rooted on every node in the tree
    :param tree_info:
    :return:
    """

    tree_number = 0
    list_of_already_used_attachments = Autovivify()
    rooted_trees = Autovivify()
    
    for node in sorted(tree_info['quartets'].keys(), key=int):
        if node in list_of_already_used_attachments:
            continue
        for attachment in sorted(tree_info['quartets'][node].keys(), key=int):
            list_of_already_used_attachments[attachment] = 1
            tree_string = ''
            node_infos = Autovivify()
            node_infos['previous_node'] = ''
            node_infos['node'] = ';'
            node_infos['open_attachments'][node] = 1
            node_infos['open_attachments'][attachment] = 1
            new_tree = recursive_tree_builder(tree_info, node_infos, tree_string)
            rooted_trees[tree_number] = new_tree
            tree_number += 1
    return rooted_trees


def recursive_tree_builder(tree_info, node_infos, tree_string):
    node = node_infos['node']
    count = 0

    for attachment in sorted(node_infos['open_attachments'].keys(), key=int):
        count += 1
        if count == 1:
            tree_string += '('
        node_infos2 = Autovivify()
        node_infos2['previous_node'] = node
        node_infos2['node'] = attachment
        count2 = 0

        for attachment_of_used_attachment in sorted(tree_info['quartets'][attachment].keys()):
            if attachment_of_used_attachment in node_infos['open_attachments']:
                continue
            if attachment_of_used_attachment == node:
                continue
            count2 += 1
            node_infos2['open_attachments'][attachment_of_used_attachment] = 1

        if count2 > 0:
            tree_string = recursive_tree_builder(tree_info, node_infos2, tree_string)
        else:
            tree_string += str(attachment)
        if count == 1:
            tree_string += ','
        if count == 2:
            tree_string += ')' + str(node)

    return tree_string


def parallel_subtree_node_retriever(rooted_trees, num_threads, parse_log):
    """
    Run `get_node_subtrees` in parallel for each of the elements in rooted_trees
    :param rooted_trees: Dictionary of rooted trees
    :param num_threads: Number of threads to use
    :param parse_log: The file object to write parsing information to
    :return: rooted_tree_nodes - a list of results from get_node_subtrees(), one for each rooted_tree
    """
    job_queue = JoinableQueue()
    result_queue = Queue()
    rooted_tree_nodes = list()

    worker_group = [NodeRetrieverWorker(job_queue, result_queue) for i in range(int(num_threads))]
    for worker in worker_group:
        worker.start()

    # tasks = [split_tree_string(rooted_trees[rooted_tree]) for rooted_tree in rooted_trees.keys()]
    # tasks = [rooted_trees[rooted_tree] for rooted_tree in rooted_trees.keys()]
    # for task in tasks:
    #     print "Input: " + task
    tasks = rooted_trees.values()
    parse_log.write("Number of subtrees = " + str(len(tasks)) + "\n")
    parse_log.flush()
    for rooted_tree_elements in tasks:
        if job_queue.full():
            sys.exit("ERROR: multiprocessing.Queue full in parallel_subtree_node_retriever!")
        job_queue.put(rooted_tree_elements)

    for i in range(int(num_threads)):
        job_queue.put(None)

    for i in range(len(rooted_trees.keys())):
        rooted_tree_nodes.append(result_queue.get())

    job_queue.close()
    result_queue.close()
    result_queue.join_thread()
    return rooted_tree_nodes


def build_terminal_children_strings_of_assignments(rooted_trees, insertion_point_node_hash,
                                                   assignments, num_threads, parse_log):
    """
    Performed for each gene (f_contig) identified
    :param rooted_trees: All possible rooted trees for a given tree (with sequence inserted)
    :param insertion_point_node_hash:
    :param assignments: The node that is inserted into the RAxML tree - found in *RAxML_classification.txt for f_contig
    :param num_threads: Number of threads to use for parsing the subtrees of each node in parallel
    :param parse_log: Name of the RAxML_output parse log file to write to
    :return:
    """
    terminal_children_strings_of_assignments = Autovivify()

    for assignment in sorted(assignments.keys()):
        internal_node_of_assignment = insertion_point_node_hash[assignment]
        # parse_log.write("Starting to retrieve all subtrees at " + time.ctime())
        rooted_tree_nodes = parallel_subtree_node_retriever(rooted_trees, num_threads, parse_log)
        # parse_log.write("Finished retrieving subtrees at " + time.ctime() + "\n")
        for rooted_tree_info in rooted_tree_nodes:
            assignment_subtree = str(rooted_tree_info['subtree_of_node'][str(internal_node_of_assignment)])
            terminal_children = Autovivify()

            if re.search(r'\A(\d+)\Z', assignment_subtree):
                terminal_children[re.search(r'\A(\d+)\Z', assignment_subtree).group(1)] = 1
            else:
                for each_hit in re.findall(r'(\D)(\d+)', assignment_subtree):
                    if each_hit[0] == '-':
                        continue
                    terminal_children[each_hit[1]] = 1

            terminal_children_string_of_assignment = ''

            # terminal_children_string_of_assignment = ' '.join(sorted(terminal_children.keys(), key=int))
            for terminal_child_of_assignment in sorted(terminal_children.keys(), key=int):
                terminal_children_string_of_assignment += str(terminal_child_of_assignment) + ' '

            terminal_children_strings_of_assignments[assignment][terminal_children_string_of_assignment] = 1

    return terminal_children_strings_of_assignments


def build_terminal_children_strings_of_reference_nodes(reference_tree_info):
    terminal_children_strings_of_reference = Autovivify()

    for node in sorted(reference_tree_info['subtree_of_node'].keys()):
        reference_subtree = reference_tree_info['subtree_of_node'][node]
        terminal_children = Autovivify()
        if re.search(r'\A(\d+)\Z', str(reference_subtree)):
            terminal_children[re.search(r'\A(\d+)\Z', str(reference_subtree)).group(1)] = 1
        else:

            for each_hit in re.findall(r'(.)(\d+)', str(reference_subtree)):
                if each_hit[0] == '-':
                    continue
                terminal_children[each_hit[1]] = 1

        terminal_children_string_of_reference = ''
        # terminal_children_string_of_reference = ' '.join(sorted(terminal_children.keys(), key=int))
        for terminal_child_of_reference in sorted(terminal_children.keys(), key=int):
            terminal_children_string_of_reference += str(terminal_child_of_reference) + ' '

        terminal_children_strings_of_reference[terminal_children_string_of_reference] = 1

    return terminal_children_strings_of_reference


def compare_terminal_children_strings(terminal_children_of_assignments, terminal_children_of_reference, parse_log):
    real_terminal_children_of_assignments = Autovivify()
    there_was_a_hit = 0
    parse_log.write("compare_terminal_children_strings\tstart: ")
    parse_log.write(time.ctime())
    for assignment in sorted(terminal_children_of_assignments.keys()):
        real_terminal_children_string = ''

        for terminal_children_string_of_assignment in sorted(terminal_children_of_assignments[assignment].keys()):
            if terminal_children_string_of_assignment in terminal_children_of_reference:
                real_terminal_children_string = terminal_children_string_of_assignment
                real_terminal_children_of_assignments[assignment] = real_terminal_children_string
                there_was_a_hit = 1
                break

        if str(real_terminal_children_string) == '' and not str(assignment) == 'mp_root':
            sys.exit('ERROR: The RAxML output tree could not be rooted correctly!!!\n')

    if there_was_a_hit <= 0:
        sys.exit('ERROR: The RAxML output tree could not be rooted correctly!!!\n')

    parse_log.write("\tstop: " + time.ctime() + "\n")
    parse_log.flush()
    return real_terminal_children_of_assignments


def delete_files(args, section):
    files_to_be_deleted = []
    if args.delete:
        if section == 1:
            files_to_be_deleted += glob.glob(args.output_dir_var + '*BLAST_results*')
            files_to_be_deleted += glob.glob(args.output_dir_var + '*blast_result_purified.txt')
        if section == 2:
            files_to_be_deleted += glob.glob(args.output_dir_var + '*_sequence.txt')
            files_to_be_deleted += glob.glob(args.output_dir_var + '*sequence_shortened.txt')
        if section == 3:
            files_to_be_deleted += glob.glob(args.output_dir_var + '*genewise.txt')
            files_to_be_deleted += glob.glob(args.output_dir_var + '*genewise_result_summary.txt')
            files_to_be_deleted += glob.glob(args.output_dir_var + '*rRNA_result_summary.txt')
        if section == 4:
            files_to_be_deleted += glob.glob(args.output_dir_var + '*.mfa')
            files_to_be_deleted += glob.glob(args.output_dir_var + '*.mfa-gb')
            files_to_be_deleted += glob.glob(args.output_dir_var + '*.mfa-gb.txt')
            if args.rpkm:
                files_to_be_deleted += glob.glob(args.output + "RPKM_outputs" + os.sep + "*.sam")
        if section == 5:
            files_to_be_deleted += glob.glob(args.output_dir_var + '*_exit_after_trimal.txt')
            files_to_be_deleted += glob.glob(args.output_dir_var + '*.mfa-trimal')
            files_to_be_deleted += glob.glob(args.output_dir_var + '*.cl')
            files_to_be_deleted += glob.glob(args.output_dir_var + '*_RAxML.txt')
            files_to_be_deleted += glob.glob(args.output_dir_var + 'RAxML_entropy.*')
            files_to_be_deleted += glob.glob(args.output_dir_var + '*RAxML_info.txt')
            files_to_be_deleted += glob.glob(args.output_dir_var + '*RAxML_labelledTree.txt')
            files_to_be_deleted += glob.glob(args.output_dir_var + 'RAxML_classificationLikelihoodWeights*')
            files_to_be_deleted += glob.glob(args.output_dir_var + '*.phy')
            files_to_be_deleted += glob.glob(args.output_dir_var + '*.phy.reduced')
            files_to_be_deleted += glob.glob(args.output_dir_var + '*RAxML_classification.txt')
            # Need this for annotate_extra_treesapp.py
            # files_to_be_deleted += glob.glob(args.output_dir_var + '*.jplace')

    for useless_file in files_to_be_deleted:
        if path.exists(useless_file):
            os.remove(useless_file)


def single_family_msa(args, cog_list, formatted_fasta_dict):
    """
    A wrapper function for hmmalign -- to generate a multiple-sequence alignment with the reference sequences
    of the gene family being updated
    :param args: Command-line argument object returned by argparse
    :param cog_list: The reference gene family to be updated
    :param formatted_fasta_dict: keys are fasta headers; values are fasta sequences. Returned by format_read_fasta
    :return: An Autovivification mapping the summary files to each contig
    """
    hmmalign_singlehit_files = Autovivify()
    logging.info("Running hmmalign... ")

    cog = list(cog_list["all_cogs"].keys())[0]
    denominator = cog_list["all_cogs"][cog]

    start = 0

    # Imitate the Genewise / blastp_summary_files output
    oddly_long = False
    for contig in formatted_fasta_dict.keys():
        header = contig[1:]
        sequence = formatted_fasta_dict[contig]
        end = len(sequence)

        if end > 3000:
            oddly_long = True

        f_contig = denominator + "_" + header
        genewise_singlehit_file = args.output_dir_var + os.sep + \
                                  f_contig + '_' + cog + "_" + str(start) + "_" + str(end)
        hmmalign_singlehit_files[f_contig][genewise_singlehit_file + ".mfa"] = True
        genewise_singlehit_file_fa = genewise_singlehit_file + ".fa"
        try:
            outfile = open(genewise_singlehit_file_fa, 'w')
            fprintf(outfile, '>query\n%s\n', sequence)
            outfile.close()
        except IOError:
            sys.stderr.write('Can\'t create ' + genewise_singlehit_file_fa + '\n')
            sys.exit(0)
        treesapp_resources = args.treesapp + os.sep + 'data' + os.sep
        hmmalign_command = [args.executables["hmmalign"], '-m', '--mapali',
                            treesapp_resources + 'alignment_data' + os.sep + cog + '.fa',
                            '--outformat', 'Clustal',
                            treesapp_resources + 'hmm_data' + os.sep + cog + '.hmm',
                            genewise_singlehit_file_fa, '>', genewise_singlehit_file + '.mfa']
        os.system(' '.join(hmmalign_command))

    logging.info("done.\n")

    if oddly_long:
        sys.stderr.write("WARNING: These sequences look awfully long for a gene... "
                         "are you sure you want to be running in this mode?\n")
        sys.stderr.flush()

    return hmmalign_singlehit_files


def num_sequences_fasta(fasta):
    fasta_file_handle = open(fasta, "r")
    fasta_lines = fasta_file_handle.readlines()
    fasta_file_handle.close()

    num_seqs = 0
    for fasta_line in fasta_lines:
        if re.search("^>", fasta_line):
            num_seqs += 1

    return num_seqs


def read_marker_classification_table(assignment_file):
    """
    Function for reading the tabular assignments file (currently marker_contig_map.tsv)
    Assumes column 2 is the TreeSAPP assignment and column 3 is the sequence header
    (leaving 1 for marker name and 4 for numerical abundance)

    :param assignment_file: Path to the file containing sequence phylogenetic origin and assignment
    :return: dictionary whose keys are phylogenetic origin and values are lists of TreeSAPP assignments
    """
    assignments = dict()
    n_classified = 0
    assignments_handle = open(assignment_file, 'r')
    header = "Sample\tQuery\tMarker\tLength\tTaxonomy\tConfident_Taxonomy\tAbundance\tiNode\tLWR\tEvoDist\tDistances\n"
    # This is the header line
    if assignments_handle.readline() != header:
        logging.error("Header of classification file is unexpected!\n")
        raise AssertionError

    # First line in the table containing data
    line = assignments_handle.readline()
    while line:
        fields = line.strip().split('\t')
        try:
            _, header, marker, lineage, _, _, _, _, _, _ = fields
            header = '>' + header
            if marker and lineage:
                n_classified += 1
                if marker not in assignments:
                    assignments[marker] = dict()
                if header in assignments[marker]:
                    raise AssertionError("ERROR: multiple " + marker + " annotations for " + header[1:] + '.')
                assignments[marker][header] = lineage
        except ValueError:
            sys.stderr.write("ERROR: Unable to parse line:\n")
            sys.stderr.write(str(line))
            sys.exit(1)
        line = assignments_handle.readline()

    assignments_handle.close()
    return assignments, n_classified


def get_new_ref_sequences(args, update_tree):
    """
    Function for retrieving the protein sequences from the TreeSAPP hmm-purified outputs
    :param args: Command-line argument object from get_options and check_parser_arguments
    :param update_tree: An instance of CreateFuncTreeUtility class
    :return: aa_dictionary is a dictionary of fasta sequences with headers as keys and protein sequences as values
    """
    logging.info("Retrieving candidate reference sequences... ")

    candidate_fa = update_tree.InputData + os.sep + "final_outputs" + os.sep + update_tree.COG + "_hmm_purified.fasta"
    aa_dictionary = dict()
    if os.path.isfile(candidate_fa):
        candidate_fa_handler = open(candidate_fa, 'r')
        sequence = ""
        header = ""
        line = candidate_fa_handler.readline().strip()
        while line:
            if line[0] == '>':
                if header and sequence:
                    aa_dictionary[header] = sequence
                header = line.split('|')[0]
                sequence = ""
            else:
                sequence += line.strip()
            line = candidate_fa_handler.readline().strip()
        candidate_fa_handler.close()
        aa_dictionary[header] = sequence
    else:
        logging.error("File necessary for updating the reference tree (" + candidate_fa + ") is missing!\n")
        sys.exit(3)

    logging.info("done.\n")
    logging.debug("\t" + str(len(aa_dictionary)) + " candidate " + update_tree.COG + " reference sequences.\n")
    sys.stdout.flush()

    return aa_dictionary


def cluster_new_reference_sequences(update_tree, args, new_ref_seqs_fasta):
    logging.info("Clustering sequences at %s percent identity with USEARCH... " % str(update_tree.cluster_id))

    usearch_command = [args.executables["usearch"]]
    usearch_command += ["-sortbylength", new_ref_seqs_fasta]
    usearch_command += ["-fastaout", update_tree.Output + "usearch_sorted.fasta"]
    usearch_command += ["--log", update_tree.Output + os.sep + "usearch_sort.log"]
    # usearch_command += ["1>", "/dev/null", "2>", "/dev/null"]

    launch_write_command(usearch_command)

    uclust_id = "0." + str(int(update_tree.cluster_id))
    try:
        float(uclust_id)
    except ValueError:
        logging.error("Weird formatting of cluster_id: " + uclust_id + "\n")

    uclust_command = [args.executables["usearch"]]
    uclust_command += ["-cluster_fast", update_tree.Output + "usearch_sorted.fasta"]
    uclust_command += ["--id", uclust_id]
    uclust_command += ["--centroids", update_tree.Output + "uclust_" + update_tree.COG + ".fasta"]
    uclust_command += ["--uc", update_tree.Output + "uclust_" + update_tree.COG + ".uc"]
    uclust_command += ["--log", update_tree.Output + os.sep + "usearch_cluster.log"]
    # uclust_command += ["1>", "/dev/null", "2>", "/dev/null"]

    launch_write_command(uclust_command)

    logging.info("done.\n")

    return


def filter_short_sequences(args, aa_dictionary, length_threshold):
    """
    Removes all sequences shorter than length_threshold from a dictionary
    :param args: Command-line argument object from get_options and check_parser_arguments
    :param aa_dictionary: Dictionary containing all candidate reference sequences from a TreeSAPP analysis
    :param length_threshold: Minimum number of AA a sequence must contain to be included in further analyses
    :return: dictionary with sequences only longer than length_threshold
    """
    long_queries = dict()
    short_seqs = 0
    logging.info("Removing all sequences shorter than " + str(length_threshold) + " amino acids... ")

    for seq in aa_dictionary:
        if len(aa_dictionary[seq]) >= length_threshold:
            long_queries[seq] = aa_dictionary[seq]
        else:
            short_seqs += 1

    logging.info("done.\n")
    logging.info("\t" + str(short_seqs) + " were removed.\n")
    if len(long_queries.keys()) == 0:
        logging.warning("No sequences passed the minimum length threshold! Skipping updating.\n")
        return

    return long_queries


def write_classified_nuc_sequences(tree_saps, nuc_orfs_formatted_dict, orf_nuc_fasta):
    """
    Function to write the nucleotide sequences representing the full-length ORF for each classified sequence
    :param tree_saps: A dictionary of gene_codes as keys and TreeSap objects as values
    :param nuc_orfs_formatted_dict:
    :param orf_nuc_fasta:
    :return: None
    """
    # Header format:
    # >contig_name|marker_gene

    output_fasta_string = ""

    for denominator in tree_saps:
        for placed_sequence in tree_saps[denominator]:
            if placed_sequence.classified:
                try:
                    output_fasta_string += '>' + placed_sequence.contig_name + '|' + placed_sequence.name + "\n"
                    output_fasta_string += nuc_orfs_formatted_dict['>' + placed_sequence.contig_name] + "\n"
                except KeyError:
                    logging.error("Unable to find '>" + placed_sequence.contig_name + "' in predicted ORFs file!\n" +
                                  "Example headers in the predicted ORFs file:\n\t" +
                                  '\n\t'.join(list(nuc_orfs_formatted_dict.keys())[:6]) + "\n")
                    sys.exit(3)

    if output_fasta_string:
        try:
            fna_output = open(orf_nuc_fasta, 'w')
        except IOError:
            logging.error("Unable to open " + orf_nuc_fasta + " for writing!")
            sys.exit(3)

        fna_output.write(output_fasta_string)
        fna_output.close()

    return


def align_reads_to_nucs(args, reference_fasta):
    """
    Align the predicted ORFs to the reads using BWA MEM
    :param args: Command-line argument object from get_options and check_parser_arguments
    :param reference_fasta: A FASTA file containing the sequences to be aligned to
    :return: Path to the SAM file
    """
    rpkm_output_dir = args.output + "RPKM_outputs" + os.sep
    if not os.path.exists(rpkm_output_dir):
        try:
            os.makedirs(rpkm_output_dir)
        except OSError:
            if os.path.exists(rpkm_output_dir):
                logging.warning("Overwriting files in " + rpkm_output_dir + ".\n")
            else:
                raise OSError("Unable to make " + rpkm_output_dir + "!\n")

    logging.info("Aligning reads to ORFs with BWA MEM... ")

    sam_file = rpkm_output_dir + '.'.join(os.path.basename(reference_fasta).split('.')[0:-1]) + ".sam"
    index_command = [args.executables["bwa"], "index"]
    index_command += [reference_fasta]
    index_command += ["1>", "/dev/null", "2>", args.output + "treesapp_bwa_index.stderr"]

    launch_write_command(index_command)

    bwa_command = [args.executables["bwa"], "mem"]
    bwa_command += ["-t", str(args.num_threads)]
    if args.pairing == "pe" and not args.reverse:
        bwa_command.append("-p")
        logging.debug("FASTQ file containing reverse mates was not provided - assuming the reads are interleaved!\n")
    elif args.pairing == "se":
        bwa_command += ["-S", "-P"]

    bwa_command.append(reference_fasta)
    bwa_command.append(args.reads)
    if args.pairing == "pe" and args.reverse:
        bwa_command.append(args.reverse)
    bwa_command += ["1>", sam_file, "2>", args.output + "treesapp_bwa_mem.stderr"]

    p_bwa = subprocess.Popen(' '.join(bwa_command), shell=True, preexec_fn=os.setsid)
    p_bwa.wait()
    if p_bwa.returncode != 0:
        logging.error("bwa mem did not complete successfully for:\n" +
                      str(' '.join(bwa_command)) + "\n")

    logging.info("done.\n")

    return sam_file


def run_rpkm(args, sam_file, orf_nuc_fasta):
    """
    Calculate RPKM values using the rpkm executable
    :param args: Command-line argument object from get_options and check_parser_arguments
    :param sam_file:
    :param orf_nuc_fasta:
    :return: Path to the RPKM output csv file
    """
    logging.info("Calculating RPKM values for each ORF... ")

    rpkm_output_file = '.'.join(sam_file.split('.')[0:-1]) + ".csv"
    rpkm_output_dir = args.output + "RPKM_outputs" + os.sep

    rpkm_command = [args.executables["rpkm"]]
    rpkm_command += ["-c", orf_nuc_fasta]
    rpkm_command += ["-a", sam_file]
    rpkm_command += ["-o", rpkm_output_file]
    rpkm_command += ["1>", rpkm_output_dir + "rpkm_stdout.txt", "2>", rpkm_output_dir + "rpkm_stderr.txt"]

    p_rpkm = subprocess.Popen(' '.join(rpkm_command), shell=True, preexec_fn=os.setsid)
    p_rpkm.wait()
    if p_rpkm.returncode != 0:
        logging.error("RPKM calculation did not complete successfully for:\n" +
                      str(' '.join(rpkm_command)) + "\n")
        sys.exit(3)
    logging.info("done.\n")

    return rpkm_output_file


def summarize_placements_rpkm(args, rpkm_output_file, marker_build_dict):
    """
    Recalculates the percentages for each marker gene final output based on the RPKM values
    Recapitulates MLTreeMap standard out summary
    :param args: Command-line argument object from get_options and check_parser_arguments
    :param rpkm_output_file: CSV file containing contig names and RPKM values
    :type marker_build_dict: dict
    :param marker_build_dict:
    :return:
    """
    contig_rpkm_map = dict()
    marker_contig_map = dict()
    contig_placement_map = dict()
    placement_rpkm_map = dict()
    marker_rpkm_map = dict()

    # Pull the RPKM values for each marker predicted; seq_name format is contig|marker
    rpkm_values = read_rpkm(rpkm_output_file)
    for seq_name in rpkm_values:
        contig, marker = seq_name.split('|')
        contig_rpkm_map[contig] = rpkm_values[seq_name]
        if marker not in marker_contig_map:
            marker_contig_map[marker] = list()
        marker_contig_map[marker].append(contig)

    # Identify the internal node each sequence was placed, used for
    final_raxml_outputs = os.listdir(args.output_dir_raxml)
    for raxml_contig_file in final_raxml_outputs:
        contig_name = '_'.join(re.sub("_RAxML_parsed.txt", '', raxml_contig_file).split('_')[1:])
        try:
            contig_placement = open(args.output_dir_raxml + raxml_contig_file, 'r')
        except IOError:
            raise IOError("Unable to open " + args.output_dir_raxml + raxml_contig_file + " for reading!")
        line = contig_placement.readline()
        while not line.startswith("Placement"):
            line = contig_placement.readline().strip()

        placement = re.sub("^.*: ", '', line)
        contig_placement_map[contig_name] = placement
        if placement not in placement_rpkm_map:
            placement_rpkm_map[placement] = 0
        contig_placement.close()

    for marker in marker_contig_map:
        marker_rpkm_total = 0
        marker_rpkm_map[marker] = dict()
        for contig in marker_contig_map[marker]:
            if contig in contig_placement_map:
                placement = contig_placement_map[contig]
                placement_rpkm_map[placement] += float(contig_rpkm_map[contig])
                marker_rpkm_total += float(contig_rpkm_map[contig])
                marker_rpkm_map[marker][placement] = 0
        for placement in marker_rpkm_map[marker]:
            try:
                percentage = (placement_rpkm_map[placement]*100)/marker_rpkm_total
            except ZeroDivisionError:
                percentage = 0
            marker_rpkm_map[marker][placement] = percentage

    for marker in marker_rpkm_map:
        ref_marker = None
        for marker_code in marker_build_dict:
            if marker_build_dict[marker_code].cog == marker:
                ref_marker = marker_build_dict[marker_code]
                break
        if not ref_marker:
            raise AssertionError("Unable to map marker '" + marker + "' to a MarkerBuild instance.")

        final_output_file = args.output_dir_final + str(ref_marker.denominator) + "_concatenated_RAxML_outputs.txt"
        # Not all of the genes predicted will have made it to the RAxML stage
        if os.path.isfile(final_output_file):
            shutil.move(final_output_file, args.output_dir_final + ref_marker.denominator + "_concatenated_counts.txt")
            try:
                cat_output = open(final_output_file, 'w')
            except IOError:
                raise IOError("Unable to open " + final_output_file + " for writing!")

            description_text = '# ' + str(ref_marker.kind) + '\n\n'
            cat_output.write(description_text)

            for placement in sorted(marker_rpkm_map[marker].keys(), reverse=True):
                relative_weight = marker_rpkm_map[marker][placement]
                if relative_weight > 0:
                    cat_output.write('Placement weight ')
                    cat_output.write('%.2f' % relative_weight + "%: ")
                    cat_output.write(placement + "\n")

            cat_output.close()
            
    return


def get_reference_sequence_dict(args, update_tree):
    # Determine the name of the aligned FASTA with reference sequences
    ref_alignment_fasta = "data" + os.sep + "alignment_data" + os.sep + update_tree.COG + ".fa"
    # Read the FASTA to get headers and sequences
    ref_fasta_dict = format_read_fasta(ref_alignment_fasta, update_tree.marker_molecule, args.output)
    # Strip the '-'s since these will be re-aligned again
    unaligned_ref_seqs = {header: re.sub('-', '', ref_fasta_dict[header]) for header in ref_fasta_dict}

    return unaligned_ref_seqs


def update_func_tree_workflow(args, ref_marker: MarkerBuild):

    # Load information essential to updating the reference data into a CreateFuncTreeUtility class object
    update_tree = CreateFuncTreeUtility(args.output, ref_marker.denominator)

    # Get HMM, sequence, reference build, and taxonomic information for the original sequences
    ref_hmm_file = args.treesapp + os.sep + 'data' + os.sep + "hmm_data" + os.sep + update_tree.COG + ".hmm"
    hmm_length = get_hmm_length(ref_hmm_file)
    unaligned_ref_seqs = get_reference_sequence_dict(args, update_tree)
    # read_species_translation_files expects the entire marker_build_dict, so we're making a mock one
    ref_organism_lineage_info = read_species_translation_files(args, {ref_marker.denominator: ref_marker})

    # Set up the output directories
    time_of_run = strftime("%d_%b_%Y_%H_%M", gmtime())
    project_folder = update_tree.Output + str(time_of_run) + os.sep
    raxml_destination_folder = project_folder + "phy_files_%s" % update_tree.COG
    final_tree_dir = project_folder + "tree_data" + os.sep
    alignment_files_dir = project_folder + "alignment_data" + os.sep
    hmm_files_dir = project_folder + "hmm_data" + os.sep
    classification_table = update_tree.InputData + os.sep + "final_outputs" + os.sep + "marker_contig_map.tsv"
    os.makedirs(project_folder)
    os.makedirs(final_tree_dir)
    os.makedirs(alignment_files_dir)
    os.makedirs(hmm_files_dir)

    # Begin finding and filtering the new candidate reference sequences
    aa_dictionary = get_new_ref_sequences(args, update_tree)
    assignments, n_classified = read_marker_classification_table(classification_table)
    if len(aa_dictionary) == 0:
        sys.stderr.write("WARNING: No new " + update_tree.COG + " sequences. Skipping update.\n")
        return
    if n_classified == 0 or update_tree.COG not in assignments.keys():
        sys.stderr.write("WARNING: No " + update_tree.COG + " sequences were classified. Skipping update.\n")
        return
    aa_dictionary = filter_short_sequences(args, aa_dictionary, 0.5*hmm_length)
    if not aa_dictionary:
        return
    new_ref_seqs_fasta = update_tree.Output + os.path.basename(update_tree.InputData) + \
                         "_" + update_tree.COG + "_unaligned.fasta"
    # Write only the sequences that have been properly classified
    write_new_fasta(aa_dictionary, new_ref_seqs_fasta, None, list(assignments[update_tree.COG].keys()))
    # Make sure the tree is updated only if there are novel sequences (i.e. <97% similar to ref sequences)
    ref_candidate_alignments = align_ref_queries(args, new_ref_seqs_fasta, update_tree)
    # Get the sequences that pass the similarity threshold
    new_refs = find_novel_refs(ref_candidate_alignments, aa_dictionary, update_tree)
    write_new_fasta(new_refs, new_ref_seqs_fasta)
    if args.uclust and len(new_refs.keys()) > 1:
        cluster_new_reference_sequences(update_tree, args, new_ref_seqs_fasta)
        centroids_fasta = update_tree.Output + "uclust_" + update_tree.COG + ".fasta"
    else:
        if len(aa_dictionary) == 1 and args.uclust:
            sys.stderr.write("WARNING: Not clustering new " + update_tree.COG + " since there is 1 sequence\n")
            sys.stderr.flush()
        centroids_fasta = new_ref_seqs_fasta

    # The candidate set has been finalized. Begin rebuilding!
    update_tree.load_new_refs_fasta(args, centroids_fasta, ref_organism_lineage_info)
    aligned_fasta = update_tree.align_multiple_sequences(unaligned_ref_seqs, args)
    trimal_file = trim_multiple_alignment(args.executables["BMGE.jar"], aligned_fasta, update_tree.marker_molecule)

    shutil.move(trimal_file, alignment_files_dir + update_tree.COG + ".fa")
    aligned_fasta = alignment_files_dir + update_tree.COG + ".fa"
    update_tree.update_tax_ids(args, ref_organism_lineage_info, assignments)

    new_hmm_file = update_tree.Output + os.sep + update_tree.COG + ".hmm"
    build_hmm(args, alignment_files_dir + update_tree.COG + ".fa", new_hmm_file)
    new_hmm_length = get_hmm_length(new_hmm_file)
    logging.debug("\tOld HMM length = " + str(hmm_length) + "\n" +
                  "\tNew HMM length = " + str(new_hmm_length) + "\n")

    os.system('java -cp sub_binaries/readseq.jar run -a -f=12 %s' % aligned_fasta)

    phylip_file = update_tree.Output + "%s.phy" % update_tree.COG
    os.system('mv %s.phylip %s' % (aligned_fasta, phylip_file))

    update_tree.execute_raxml(phylip_file, raxml_destination_folder, args)

    # Organize outputs
    shutil.move(new_hmm_file, hmm_files_dir)
    shutil.move(update_tree.Output + "tax_ids_" + update_tree.COG + ".txt", final_tree_dir)

    best_tree = raxml_destination_folder + "/RAxML_bestTree." + update_tree.COG
    bootstrap_tree = raxml_destination_folder + "/RAxML_bipartitionsBranchLabels." + update_tree.COG
    best_tree_nameswap = final_tree_dir + update_tree.COG + "_tree.txt"
    bootstrap_nameswap = final_tree_dir + update_tree.COG + "_bipartitions.txt"
    update_tree.swap_tree_names(best_tree, best_tree_nameswap)
    update_tree.swap_tree_names(bootstrap_tree, bootstrap_nameswap)
    annotate_partition_tree(update_tree.COG,
                            update_tree.master_reference_index,
                            raxml_destination_folder + os.sep + "RAxML_bipartitions." + update_tree.COG)

    prefix = update_tree.Output + update_tree.COG
    os.system('mv %s* %s' % (prefix, project_folder))

    if args.uclust:
        uclust_output_dir = prefix + "_uclust"
        os.system('mkdir %s' % uclust_output_dir)

        os.system('mv %suclust_* %s' % (update_tree.Output, uclust_output_dir))
        os.system('mv %susearch_* %s' % (update_tree.Output, uclust_output_dir))

    intermediate_files = [project_folder + update_tree.COG + ".phy",
                          project_folder + update_tree.COG + "_gap_removed.fa",
                          project_folder + update_tree.COG + "_d_aligned.fasta"]
    for useless_file in intermediate_files:
        try:
            os.remove(useless_file)
        except OSError:
            sys.stderr.write("WARNING: unable to remove intermediate file " + useless_file + "\n")


def create_itol_labels(args, marker):
    """
    
    :param args: Command-line argument object from get_options and check_parser_arguments
    :param marker:
    :return: None
    """
    itol_base_dir = args.output + 'iTOL_output' + os.sep
    itol_label_file = itol_base_dir + os.sep + marker + os.sep + marker + "_labels.txt"
    tax_ids_file = os.sep.join([args.treesapp, "data", "tree_data", "tax_ids_" + marker + ".txt"])

    if os.path.exists(itol_label_file):
        return

    try:
        label_f = open(itol_label_file, 'w')
    except IOError:
        raise IOError("Unable to open " + itol_label_file + " for writing! Exiting now.")

    try:
        tax_ids = open(tax_ids_file, 'r')
    except IOError:
        raise IOError("Unable to open " + tax_ids_file + " for reading!")

    label_f.write("LABELS\nSEPARATOR COMMA\nDATA\n#NODE_ID,LABEL\n")
    for line in tax_ids:
        line = line.strip()
        try:
            fields = line.split("\t")
        except ValueError:
            sys.stderr.write('ValueError: .split(\'\\t\') on ' + str(line) +
                             " generated " + str(len(line.split("\t"))) + " fields.")
            sys.exit(9)
        if len(fields) == 2:
            number, translation = fields
        elif len(fields) == 3:
            number, translation, lineage = fields
        else:
            sys.stderr.write("ValueError: Unexpected number of fields in " + tax_ids_file +
                             ".\nInvoked .split(\'\\t\') on line " + str(line))
            raise ValueError
        label_f.write(number + ',' + translation + "\n")

    tax_ids.close()
    label_f.close()

    return


def generate_simplebar(target_marker, tree_protein_list, itol_bar_file, all_rpkm_values=None):
    """
    From the basic RPKM output csv file, generate an iTOL-compatible simple bar-graph file for each leaf

    :param all_rpkm_values: A dictionary mapping seq_names to RPKM floats
    :param target_marker:
    :param tree_protein_list: A list of TreeProtein objects, for single sequences
    :param itol_bar_file: The name of the file to write the simple-bar data for iTOL
    :return:
    """
    rpkm_values = dict()
    leaf_rpkm_sums = dict()

    if all_rpkm_values:
        for seq_name in all_rpkm_values:
            contig, marker = seq_name.split('|')
            # Filter out RPKMs for contigs not associated with the target marker
            if target_marker == marker:
                rpkm_values[contig] = all_rpkm_values[seq_name]
    else:
        for tree_sap in tree_protein_list:
            rpkm_values[tree_sap.contig_name] = 1.0

    for tree_sap in tree_protein_list:
        if tree_sap.classified:
            if tree_sap.contig_name in rpkm_values:
                tree_sap.abundance = rpkm_values[tree_sap.contig_name]
            else:
                tree_sap.abundance = 0
            leaf_rpkm_sums = tree_sap.sum_rpkms_per_node(leaf_rpkm_sums)

    # Only make the file if there is something to write
    if len(leaf_rpkm_sums.keys()) > 0:
        try:
            itol_rpkm_out = open(itol_bar_file, 'w')
        except IOError:
            logging.error("Unable to open " + itol_bar_file + " for writing.\n")
            sys.exit(3)

        # Write the header
        header = "DATASET_SIMPLEBAR\nSEPARATOR COMMA\nDATASET_LABEL,RPKM\nCOLOR,#ff0000\n"
        itol_rpkm_out.write(header)
        # Write the RPKM sums for each leaf
        itol_rpkm_out.write("DATA\n")
        data_lines = [','.join([str(k), str(v)]) for k, v in leaf_rpkm_sums.items()]
        itol_rpkm_out.write("\n".join(data_lines))

        itol_rpkm_out.close()

    return tree_protein_list


def lowest_confident_taxonomy(lct, depth):
    """
    Truncates the initial taxonomic assignment to rank of depth.

    :param lct: String for the taxonomic lineage ('; ' separated)
    :param depth: The recommended depth to truncate the taxonomy
    :return: String representing 'confident' taxonomic assignment for the sequence
    """

    purified_lineage_list = clean_lineage_string(lct).split("; ")
    confident_assignment = "; ".join(purified_lineage_list[:depth])
    # For debugging
    # rank_depth = {1: "Kingdom", 2: "Phylum", 3: "Class", 4: "Order", 5: "Family", 6: "Genus", 7: "Species", 8: "Strain"}
    # if clean_lineage_string(lct) == confident_assignment:
    #     print("Unchanged: (" + rank_depth[depth] + ')', confident_assignment)
    # else:
    #     print("Adjusted: (" + rank_depth[depth] + ')', confident_assignment)

    return confident_assignment


def enumerate_taxonomic_lineages(lineage_list):
    rank = 1
    taxonomic_counts = dict()
    while rank < 9:
        for lineage in lineage_list:
            if len(lineage) < rank:
                continue
            taxonomy = "; ".join(lineage[:rank])
            if taxonomy not in taxonomic_counts:
                taxonomic_counts[taxonomy] = 0
            taxonomic_counts[taxonomy] += 1
        rank += 1
    return taxonomic_counts


def filter_placements(args, tree_saps, marker_build_dict, unclassified_counts):
    """
    Determines the total distance of each placement from its branch point on the tree
    and removes the placement if the distance is deemed too great

    :param args: Command-line argument object from get_options and check_parser_arguments
    :param tree_saps: A dictionary containing TreeProtein objects
    :param marker_build_dict: A dictionary of MarkerBuild objects (used here for lowest_confident_rank)
    :param unclassified_counts: A dictionary tracking the number of putative markers that were not classified
    :return:
    """

    logging.info("Filtering low-quality placements... ")
    distant_seqs = dict()

    for denominator in tree_saps:
        marker = marker_build_dict[denominator].cog
        distant_seqs[marker] = list()
        tree = Tree(os.sep.join([args.treesapp, "data", "tree_data", marker + "_tree.txt"]))
        for tree_sap in tree_saps[denominator]:
            # max_dist_threshold equals the maximum path length from root to tip in its clade
            max_dist_threshold = tree.get_farthest_leaf()[1]  # Too permissive of a threshold, but good for first pass
            if tree_sap.name not in unclassified_counts.keys():
                unclassified_counts[tree_sap.name] = 0
            if not tree_sap.placements:
                unclassified_counts[tree_sap.name] += 1
            elif len(tree_sap.placements) > 1:
                logging.warning("More than one placement for a single contig:\n" +
                                tree_sap.summarize())
                tree_sap.classified = False
                continue
            elif tree_sap.placements[0] == '{}':
                unclassified_counts[tree_sap.name] += 1
                tree_sap.classified = False
                continue

            distal_length = float(tree_sap.get_jplace_element("distal_length"))
            pendant_length = float(tree_sap.get_jplace_element("pendant_length"))
            # Find the length of the edge this sequence was placed onto
            leaf_children = tree_sap.node_map[int(tree_sap.inode)]
            # Find the distance away from this edge's bifurcation (if internal) or tip (if leaf)
            if len(leaf_children) > 1:
                # We need to find the LCA in the Tree instance to find the distances to tips for ete3
                parent = tree.get_common_ancestor(leaf_children)
                tip_distances = parent_to_tip_distances(parent, leaf_children)
            else:
                tree_leaf = tree.get_leaves_by_name(leaf_children[0])[0]
                sister = tree_leaf.get_sisters()[0]
                parent = tree_leaf.get_common_ancestor(sister)
                tip_distances = [0.0]

            tree_sap.avg_evo_dist = round(distal_length + pendant_length + (sum(tip_distances) / len(tip_distances)), 4)
            tree_sap.distances = str(distal_length) + ',' +\
                                 str(pendant_length) + ',' +\
                                 str(sum(tip_distances) / len(tip_distances))
            # Discard this placement as a false positive if the avg_evo_dist exceeds max_dist_threshold
            if pendant_length > max_dist_threshold:
                # print("Global", tree_sap.summarize())
                unclassified_counts[tree_sap.name] += 1
                distant_seqs[marker].append(tree_sap.contig_name)
                tree_sap.classified = False
                continue

            # Estimate the branch lengths of the clade to factor heterogeneous substitution rates
            ancestor, clade_tip_distances = find_cluster(parent)
            if not clade_tip_distances:
                # Either the parent or ancestor is the root, so there are many children. This doesn't scale well.
                # TODO: Decrease the number of leaves sampled for this distance
                leaf, dist = parent.get_farthest_leaf()
                clade_tip_distances.append(dist)
            # If the longest root-to-tip distance from the ancestral node (one-up from LCA) is exceeded, discard
            if pendant_length > max(clade_tip_distances) * 1.2 and \
                    rank_recommender(pendant_length, marker_build_dict[denominator].pfit) < 0:
                unclassified_counts[tree_sap.name] += 1
                distant_seqs[marker].append(tree_sap.contig_name)
                tree_sap.classified = False
                # print("Local", tree_sap.summarize())
    logging.info("done.\n")

    for marker in distant_seqs:
        # unclassified_counts[marker] will always be >= distant_seqs[marker]
        if unclassified_counts[marker] > 0:
            logging.debug("\t" + str(unclassified_counts[marker]) + " " + marker +
                          " sequence(s) detected but not classified.\n" +
                          marker + " queries with extremely long placement distances:\n\t" +
                          "\n\t".join(distant_seqs[marker]) + "\n")

    return tree_saps


def write_tabular_output(args, tree_saps, tree_numbers_translation, marker_build_dict):
    """


    :param args: Command-line argument object from get_options and check_parser_arguments
    :param tree_saps: A dictionary containing TreeProtein objects
    :param tree_numbers_translation: Dictionary containing taxonomic information for each leaf in the reference tree
    :param marker_build_dict: A dictionary of MarkerBuild objects (used here for lowest_confident_rank)
    :return:
    """
    leaf_taxa_map = dict()
    mapping_output = args.output_dir_final + os.sep + "marker_contig_map.tsv"
    sample_name = os.path.basename(args.output)
    if not sample_name:
        sample_name = args.output.split(os.sep)[-2]
    tab_out_string = "Sample\tQuery\tMarker\tLength\tTaxonomy\tConfident_Taxonomy\tAbundance\tiNode\tLWR\tEvoDist\tDistances\n"
    try:
        tab_out = open(mapping_output, 'w')
    except IOError:
        logging.error("Unable to open " + mapping_output + " for writing!\n")
        sys.exit(3)

    for denominator in tree_saps:
        # All the leaves for that tree [number, translation, lineage]
        leaves = tree_numbers_translation[denominator]
        lineage_complete = False
        lineage_list = list()
        # Test if the reference set have lineage information
        for leaf in leaves:
            lineage_list.append(clean_lineage_string(leaf.lineage).split('; '))
            if leaf.complete:
                lineage_complete = True
            leaf_taxa_map[leaf.number] = leaf.lineage
        taxonomic_counts = enumerate_taxonomic_lineages(lineage_list)

        for tree_sap in tree_saps[denominator]:
            if not tree_sap.classified:
                continue

            tree_sap.lineage_list = children_lineage(leaves, tree_sap.placements[0], tree_sap.node_map)

            # Based on the calculated distance from the leaves, what rank is most appropriate?
            recommended_rank = rank_recommender(tree_sap.avg_evo_dist,
                                                marker_build_dict[denominator].pfit)
            # Sequence likely isn't a FP but is highly divergent from reference set so set to Kingdom classification
            if recommended_rank < 1:
                recommended_rank = 1

            if len(tree_sap.lineage_list) == 0:
                logging.error("Unable to find lineage information for marker " +
                              denominator + ", contig " + tree_sap.contig_name + "!\n")
                sys.exit(3)
            if len(tree_sap.lineage_list) == 1:
                tree_sap.lct = tree_sap.lineage_list[0]
                tree_sap.wtd = 0.0
            if len(tree_sap.lineage_list) > 1:
                if lineage_complete:
                    lca = tree_sap.megan_lca()
                    # algorithm options are "MEGAN", "LCAp", and "LCA*" (default)
                    tree_sap.lct = lowest_common_taxonomy(tree_sap.lineage_list, lca, taxonomic_counts, "LCA*")
                    tree_sap.wtd, status = compute_taxonomic_distance(tree_sap.lineage_list, tree_sap.lct)
                    if status > 0:
                        tree_sap.summarize()
                else:
                    tree_sap.lct = "Lowest common ancestor of: "
                    tree_sap.lct += ', '.join(tree_sap.lineage_list)
                    tree_sap.wtd = 1

            # tree_sap.summarize()
            tab_out_string += '\t'.join([sample_name,
                                         tree_sap.contig_name,
                                         tree_sap.name,
                                         str(tree_sap.seq_len),
                                         clean_lineage_string(tree_sap.lct),
                                         lowest_confident_taxonomy(tree_sap.lct, recommended_rank),
                                         str(tree_sap.abundance),
                                         str(tree_sap.inode),
                                         str(tree_sap.lwr),
                                         str(tree_sap.avg_evo_dist),
                                         tree_sap.distances]) + "\n"
    tab_out.write(tab_out_string)
    tab_out.close()

    return


def parse_raxml_output(args, marker_build_dict):
    """

    :param args: Command-line argument object from get_options and check_parser_arguments
    :param marker_build_dict: Dictionary of MarkerBuild instances indexed by denominator (refpkg code e.g. M0701)
    :return:
        1. Dictionary of TreeProtein instances indexed by denominator (refpkg code e.g. M0701)
        2. Dictionary of an ItolJplace instance (values) mapped to marker name
        3. Dictionary mapping the number of sequences that were unclassified (value) to a marker (key)
    """

    logging.info('Parsing the RAxML outputs... ')

    function_start_time = time.time()

    jplace_files = glob.glob(args.output_dir_var + '*.jplace')
    jplace_collection = organize_jplace_files(jplace_files)
    itol_data = dict()  # contains all pqueries, indexed by marker name (e.g. McrA, nosZ, 16srRNA)
    tree_saps = dict()  # contains individual pquery information for each mapped protein (N==1), indexed by denominator
    unclassified_counts = dict()  # A dictionary tracking the number of putative markers that were not classified
    # Use the jplace files to guide which markers iTOL outputs should be created for
    classified_seqs = 0
    for denominator in jplace_collection:
        marker = marker_build_dict[denominator].cog
        if denominator not in tree_saps:
            tree_saps[denominator] = list()
        if marker not in unclassified_counts.keys():
            unclassified_counts[marker] = 0
        for filename in jplace_collection[denominator]:
            # Load the JSON placement (jplace) file containing >= 1 pquery into ItolJplace object
            jplace_data = jplace_parser(filename)
            # Demultiplex all pqueries in jplace_data into individual TreeProtein objects
            tree_placement_queries = demultiplex_pqueries(jplace_data)
            # Filter the placements, determine the likelihood associated with the harmonized placement
            for pquery in tree_placement_queries:
                pquery.name = marker
                pquery.filter_min_weight_threshold(args.min_likelihood)
                if not pquery.classified:
                    unclassified_counts[marker] += 1
                    logging.debug("A putative " + marker +
                                  " sequence has been unclassified due to low placement likelihood weights. " +
                                  "More info:\n" +
                                  pquery.summarize())
                    continue
                if re.match(".*_(\d+)_(\d+)$", pquery.contig_name):
                    start, end = re.match(".*_(\d+)_(\d+)$", pquery.contig_name).groups()
                    pquery.seq_len = int(end) - int(start)
                    pquery.contig_name = re.sub(r"_(\d+)_(\d+)$", '', pquery.contig_name)
                pquery.create_jplace_node_map()
                if args.placement_parser == "best":
                    pquery.filter_max_weight_placement()
                else:
                    pquery.harmonize_placements(args.treesapp)
                if pquery.classified and len(pquery.placements) != 1:
                    logging.error("Number of JPlace pqueries is " + str(len(pquery.placements)) +
                                  " when only 1 is expected at this point.\n" +
                                  pquery.summarize())
                    sys.exit(3)
                pquery.inode = str(pquery.get_jplace_element("edge_num"))
                pquery.lwr = float(pquery.get_jplace_element("like_weight_ratio"))
                pquery.likelihood = float(pquery.get_jplace_element("likelihood"))
                tree_saps[denominator].append(pquery)
                classified_seqs += 1

            if marker not in itol_data:
                itol_data[marker] = jplace_data
                itol_data[marker].name = marker
            else:
                # If a JPlace file for that tree has already been parsed, just append the placements
                itol_data[marker].placements = itol_data[marker].placements + jplace_data.placements

            # I have decided to not remove the original JPlace files since some may find these useful
            # os.remove(filename)

    logging.info("done.\n")

    function_end_time = time.time()
    hours, remainder = divmod(function_end_time - function_start_time, 3600)
    minutes, seconds = divmod(remainder, 60)
    logging.debug("\tTree parsing time required: " +
                  ':'.join([str(hours), str(minutes), str(round(seconds, 2))]) + "\n")
    logging.debug("\t" + str(len(jplace_files)) + " RAxML output files.\n" +
                  "\t" + str(classified_seqs) + " sequences classified by TreeSAPP.\n\n")

    return tree_saps, itol_data, unclassified_counts


def produce_itol_inputs(args, tree_saps, marker_build_dict, itol_data, rpkm_output_file=None):
    """
    Function to create outputs for the interactive tree of life (iTOL) webservice.
    There is a directory for each of the marker genes detected to allow the user to "drag-and-drop" all files easily

    :param args:
    :param tree_saps:
    :param marker_build_dict:
    :param itol_data:
    :param rpkm_output_file:
    :return: None
    """
    logging.debug("Generating inputs for iTOL... ")
    
    itol_base_dir = args.output + 'iTOL_output' + os.sep
    if not os.path.exists(itol_base_dir):
        os.mkdir(itol_base_dir)  # drwxr-xr-x
    # Now that all the JPlace files have been loaded, generate the abundance stats for each marker
    
    strip_missing = []
    style_missing = []
    for denominator in tree_saps:
        if len(tree_saps[denominator]) == 0:
            # No sequences that were mapped met the minimum likelihood weight ration threshold. Skipping!
            continue
        marker = marker_build_dict[denominator].cog
        if not os.path.exists(itol_base_dir + marker):
            os.mkdir(itol_base_dir + marker)

        bipartition_file = os.sep.join([args.treesapp, "data", "tree_data", marker + "_bipartitions.txt"])
        if os.path.isfile(bipartition_file):
            itol_data[marker] = add_bipartitions(itol_data[marker], bipartition_file)

        # Make a master jplace file from the set of placements in all jplace files for each marker
        master_jplace = itol_base_dir + marker + os.sep + marker + "_complete_profile.jplace"
        itol_data[marker] = filter_jplace_data(itol_data[marker], tree_saps[denominator])
        write_jplace(itol_data[marker], master_jplace)
        itol_data[marker].clear_object()
        # Create a labels file from the tax_ids_marker.txt
        create_itol_labels(args, marker)

        annotation_style_files = glob.glob(os.sep.join([args.treesapp, "data", "iTOL_datasets", marker + "*"]))
        # Copy the respective colours and styles files for each marker found to the itol_output directories
        colours_styles = os.sep.join([args.treesapp, "data", "iTOL_datasets", marker + "_colours_style.txt"])
        colour_strip = os.sep.join([args.treesapp, "data", "iTOL_datasets", marker + "_colour_strip.txt"])
        if colours_styles not in annotation_style_files:
            style_missing.append(marker)
        if colour_strip not in annotation_style_files:
            strip_missing.append(marker)

        for annotation_file in annotation_style_files:
            shutil.copy(annotation_file, itol_base_dir + marker)
        itol_bar_file = os.sep.join([args.output, "iTOL_output", marker, marker + "_abundance_simplebar.txt"])
        if args.rpkm:
            all_rpkm_values = read_rpkm(rpkm_output_file)
            generate_simplebar(marker, tree_saps[denominator], itol_bar_file, all_rpkm_values)
        else:
            generate_simplebar(marker, tree_saps[denominator], itol_bar_file)

    logging.debug("done.\n")
    if style_missing:
        logging.debug("A colours_style.txt file does not yet exist for markers:\n\t" +
                      "\n\t".join(style_missing) + "\n")
    if strip_missing:
        logging.debug("A colours_strip.txt file does not yet exist for markers:\n\t" +
                      "\n\t".join(strip_missing) + "\n")

    return


def main(argv):
    sys.stdout.write("\n##\t\t\t\tTreeSAPP\t\t\t\t##\n\n")
    sys.stdout.flush()
    # STAGE 1: Prompt the user and prepare files and lists for the pipeline
    args = get_options()
    args = check_parser_arguments(args)
    args = check_previous_output(args)

    marker_build_dict = parse_ref_build_params(args)
    marker_build_dict = parse_cog_list(args, marker_build_dict)
    tree_numbers_translation = read_species_translation_files(args, marker_build_dict)
    if args.check_trees:
        validate_inputs(args, marker_build_dict)
    if args.skip == 'n':
        # STAGE 2: Predict open reading frames (ORFs) if the input is an assembly, read, format and write the FASTA
        if args.molecule == "dna":
            # args.fasta_input is set to the predicted ORF protein sequences
            args = predict_orfs(args)
        logging.info("Formatting " + args.fasta_input + " for pipeline... ")
        formatted_fasta_dict = format_read_fasta(args.fasta_input, "prot", args.output)
        logging.info("done.\n")

        logging.info("\tTreeSAPP will analyze the " + str(len(formatted_fasta_dict)) + " sequences found in input.\n")
        if re.match(r'\A.*\/(.*)', args.fasta_input):
            input_multi_fasta = os.path.basename(args.fasta_input)
        else:
            input_multi_fasta = args.fasta_input
        args.formatted_input_file = args.output_dir_var + input_multi_fasta + "_formatted.fasta"
        formatted_fasta_files = write_new_fasta(formatted_fasta_dict, args.formatted_input_file)
        ref_alignment_dimensions = get_alignment_dims(args, marker_build_dict)

        # STAGE 3: Run hmmsearch on the query sequences to search for marker homologs
        hmm_domtbl_files = hmmsearch_orfs(args, marker_build_dict)
        hmm_matches = parse_domain_tables(args, hmm_domtbl_files)
        homolog_seq_files, numeric_contig_index = extract_hmm_matches(args, hmm_matches, formatted_fasta_dict)

        # STAGE 4: Run hmmalign or PaPaRa, and optionally BMGE, to produce the MSAs required to for the ML estimations
        create_ref_phy_files(args, homolog_seq_files, marker_build_dict, ref_alignment_dimensions)
        concatenated_msa_files = multiple_alignments(args, homolog_seq_files, marker_build_dict)
        file_types = set()
        for mc in concatenated_msa_files:
            sample_msa_file = concatenated_msa_files[mc][0]
            f_ext = sample_msa_file.split('.')[-1]
            if re.match("phy|phylip", f_ext):
                file_types.add("Phylip")
            elif re.match("sto|stockholm", f_ext):
                file_types.add("Stockholm")
            elif re.match("mfa|fa|fasta", f_ext):
                file_types.add("Fasta")
            else:
                logging.error("Unrecognized file extension: '" + f_ext + "'")
                sys.exit(3)
        if len(file_types) > 1:
            logging.error("Multiple file types detected in multiple alignment files:\n" + ','.join(file_types) + "\n")
            sys.exit(3)
        elif len(file_types) == 0:
            logging.error("No alignment files were generated!\n")
            sys.exit(3)
        else:
            file_type = file_types.pop()
        alignment_length_dict = get_sequence_counts(concatenated_msa_files, ref_alignment_dimensions,
                                                    args.verbose, file_type)

        if args.trim_align:
            tool = "BMGE"
            mfa_files = filter_multiple_alignments(args, concatenated_msa_files, marker_build_dict, tool)
            qc_ma_dict = check_for_removed_sequences(args, mfa_files, marker_build_dict)
            evaluate_trimming_performance(qc_ma_dict, alignment_length_dict, concatenated_msa_files, tool)
            phy_files = produce_phy_files(args, qc_ma_dict)
        else:
            phy_files = concatenated_msa_files
        delete_files(args, 3)

        # STAGE 5: Run RAxML to compute the ML estimations
        launch_evolutionary_placement_queries(args, phy_files, marker_build_dict)
        sub_indices_for_seq_names_jplace(args, numeric_contig_index, marker_build_dict)
    tree_saps, itol_data, unclassified_counts = parse_raxml_output(args, marker_build_dict)
    tree_saps = filter_placements(args, tree_saps, marker_build_dict, unclassified_counts)

    abundance_file = None
    if args.molecule == "dna":
        sample_name = '.'.join(os.path.basename(re.sub("_ORFs", '', args.fasta_input)).split('.')[:-1])
        orf_nuc_fasta = args.output_dir_final + sample_name + "_classified_seqs.fna"
        if not os.path.isfile(orf_nuc_fasta):
            logging.info("Creating nucleotide FASTA file of classified sequences '" + orf_nuc_fasta + "'... ")
            genome_nuc_genes_file = args.output_dir_final + sample_name + "_ORFs.fna"
            if os.path.isfile(genome_nuc_genes_file):
                nuc_orfs_formatted_dict = format_read_fasta(genome_nuc_genes_file, 'dna', args.output)
                write_classified_nuc_sequences(tree_saps, nuc_orfs_formatted_dict, orf_nuc_fasta)
                logging.info("done.\n")
            else:
                logging.info("failed.\nWARNING: Unable to read '" + genome_nuc_genes_file + "'.\n" +
                             "Cannot create the nucleotide FASTA file of classified sequences!\n")
        if args.rpkm:
            sam_file = align_reads_to_nucs(args, orf_nuc_fasta)
            abundance_file = run_rpkm(args, sam_file, orf_nuc_fasta)
            summarize_placements_rpkm(args, abundance_file, marker_build_dict)
    else:
        pass

    write_tabular_output(args, tree_saps, tree_numbers_translation, marker_build_dict)
    produce_itol_inputs(args, tree_saps, marker_build_dict, itol_data, abundance_file)
    delete_files(args, 4)

    # STAGE 6: Optionally update the reference tree
    if args.update_tree:
        for marker_code in args.targets:
            update_func_tree_workflow(args, marker_build_dict[marker_code])

    delete_files(args, 5)
    logging.info("TreeSAPP has finished successfully.\n")


if __name__ == "__main__":
    main(sys.argv[1:])
