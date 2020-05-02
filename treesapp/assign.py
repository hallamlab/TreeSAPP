#!/usr/bin/env python3

__author__ = "Connor Morgan-Lang and Kishori Konwar"
__maintainer__ = "Connor Morgan-Lang"
__license__ = "GPL-3.0"

try:
    import profile
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

    from .classy import CommandLineWorker, CommandLineFarmer, ItolJplace, NodeRetrieverWorker,\
        TreeLeafReference, TreeProtein, MarkerBuild
    from .fasta import get_headers, write_new_fasta, read_fasta_to_dict, FASTA
    from .entish import create_tree_info_hash, deconvolute_assignments, read_and_understand_the_reference_tree,\
        get_node, annotate_partition_tree, find_cluster, tree_leaf_distances, index_tree_edges
    from .external_command_interface import launch_write_command
    from .lca_calculations import lowest_common_taxonomy, weighted_taxonomic_distance
    from .jplace_utils import children_lineage, jplace_parser, demultiplex_pqueries, write_jplace, filter_jplace_data,\
        organize_jplace_files, add_bipartitions
    from .file_parsers import read_stockholm_to_dict, read_phylip_to_dict,\
        multiple_alignment_dimensions, validate_alignment_trimming
    from .phylo_dist import parent_to_tip_distances, rank_recommender
    from . import utilities
    from . import wrapper

    import _tree_parser
    import _fasta_reader
except ImportWarning:
    sys.stderr.write("Could not load some user defined module functions")
    sys.stderr.write(traceback.print_exc(10))
    sys.exit(3)


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


def replace_contig_names(numeric_contig_index: dict, fasta: FASTA):
    for marker in numeric_contig_index:
        assign_re = re.compile(r"(.*)\|{0}\|(\d+_\d+)$".format(marker))
        for neg_num_id in numeric_contig_index[marker]:
            assign_name = numeric_contig_index[marker][neg_num_id]
            seq_name, coords = assign_re.match(assign_name).groups()
            try:
                original_name = fasta.header_registry[seq_name].original
            except KeyError:
                logging.error("Unable to find TreeSAPP numerical ID '" + seq_name + "' in header registry.\n")
                sys.exit(3)
            numeric_contig_index[marker][neg_num_id] = original_name + '|' + marker + '|' + coords
    return numeric_contig_index


def extract_hmm_matches(hmm_matches: dict, fasta_dict: dict):
    """
    Function writes the sequences identified by the HMMs to output files in FASTA format.
    Full-length query sequences with homologous regions are put into two FASTA files:
        1. Numeric index headers, used for downstream phylogenetic placement
        2. Contains the original headers along with marker, start, and end positions included
    The negative integers (or numeric indexes) are stored in a dictionary and returned
    Sequences are grouped based on the location on the HMM profile they mapped to

    :param hmm_matches: Contains lists of HmmMatch objects mapped to the marker they matched
    :param fasta_dict: Stores either the original or ORF-predicted input FASTA. Headers are keys, sequences are values
    :return: List of files that go on to placement stage, dictionary mapping marker-specific numbers to contig names
    """
    logging.info("Extracting and grouping the quality-controlled sequences... ")
    extracted_seq_dict = dict()  # Keys are markers -> bin_num -> negative integers -> extracted sequences
    numeric_contig_index = dict()  # Keys are markers -> negative integers -> headers
    bins = dict()

    for marker in hmm_matches:
        if len(hmm_matches[marker]) == 0:
            continue
        if marker not in numeric_contig_index.keys():
            numeric_contig_index[marker] = dict()
        numeric_decrementor = -1
        if marker not in extracted_seq_dict:
            extracted_seq_dict[marker] = dict()

        # Algorithm for binning sequences:
        # 1. Sort HmmMatches by the proportion of the HMM profile they covered in increasing order (full-length last)
        # 2. For HmmMatch in sorted matches, determine overlap between HmmMatch and each bin's representative HmmMatch
        # 3. If overlap exceeds 80% of representative's aligned length add it to the bin, else continue
        # 4. When bins are exhausted create new bin with HmmMatch
        for hmm_match in sorted(hmm_matches[marker], key=lambda x: x.end - x.start):
            if hmm_match.desc != '-':
                contig_name = hmm_match.orf + ' ' + hmm_match.desc
            else:
                contig_name = hmm_match.orf
            # Add the query sequence to the index map
            orf_coordinates = str(hmm_match.start) + '_' + str(hmm_match.end)
            numeric_contig_index[marker][numeric_decrementor] = contig_name + '|' + marker + '|' + orf_coordinates
            # Add the FASTA record of the trimmed sequence - this one moves on for placement
            full_sequence = fasta_dict[contig_name]
            binned = False
            for bin_num in sorted(bins):
                bin_rep = bins[bin_num][0]
                overlap = min(hmm_match.pend, bin_rep.pend) - max(hmm_match.pstart, bin_rep.pstart)
                if (100*overlap)/(bin_rep.pend - bin_rep.pstart) > 80:  # 80 refers to overlap proportion with seed
                    bins[bin_num].append(hmm_match)
                    extracted_seq_dict[marker][bin_num][numeric_decrementor] = full_sequence[
                                                                               hmm_match.start - 1:hmm_match.end]
                    binned = True
                    break
            if not binned:
                bin_num = len(bins)
                bins[bin_num] = list()
                extracted_seq_dict[marker][bin_num] = dict()
                bins[bin_num].append(hmm_match)
                extracted_seq_dict[marker][bin_num][numeric_decrementor] = full_sequence[
                                                                           hmm_match.start - 1:hmm_match.end]
            numeric_decrementor -= 1

        bins.clear()
    logging.info("done.\n")

    return extracted_seq_dict, numeric_contig_index


def write_grouped_fastas(extracted_seq_dict: dict, numeric_contig_index: dict, marker_build_dict, output_dir):
    hmmalign_input_fastas = list()
    bulk_marker_fasta = dict()
    bin_fasta = dict()

    group_size_string = "Number of query sequences in each marker's group:\n"
    for marker in extracted_seq_dict:
        for group in sorted(extracted_seq_dict[marker]):
            if extracted_seq_dict[marker][group]:
                group_size_string += "\t".join([marker, str(group), str(len(extracted_seq_dict[marker][group]))]) + "\n"
    logging.debug(group_size_string + "\n")

    logging.info("Writing the grouped sequences to FASTA files... ")

    for marker in extracted_seq_dict:
        # Find the reference package instance
        ref_marker = utilities.fish_refpkg_from_build_params(marker, marker_build_dict)

        f_acc = 0  # For counting the number of files for a marker. Will exceed groups if len(queries) > len(references)
        for group in sorted(extracted_seq_dict[marker]):
            if extracted_seq_dict[marker][group]:
                group_sequences = extracted_seq_dict[marker][group]
                for num in group_sequences:
                    # Add the query sequence to the master marker FASTA with the full sequence name
                    bulk_marker_fasta[numeric_contig_index[marker][num]] = group_sequences[num]

                    # Add the query sequence to this bin's FASTA file
                    bin_fasta[str(num)] = group_sequences[num]
                    # Ensuring the number of query sequences doesn't exceed the number of reference sequences
                    if len(bin_fasta) >= ref_marker.num_reps:
                        write_new_fasta(bin_fasta, output_dir + marker + "_hmm_purified_group" + str(f_acc) + ".faa")
                        hmmalign_input_fastas.append(output_dir + marker + "_hmm_purified_group" + str(f_acc) + ".faa")
                        f_acc += 1
                        bin_fasta.clear()
                if len(bin_fasta) >= 1:
                    write_new_fasta(bin_fasta, output_dir + marker + "_hmm_purified_group" + str(f_acc) + ".faa")
                    hmmalign_input_fastas.append(output_dir + marker + "_hmm_purified_group" + str(f_acc) + ".faa")
            f_acc += 1
            bin_fasta.clear()

        # Now write a single FASTA file with all identified markers
        if len(bulk_marker_fasta) >= 1:
            trimmed_hits_fasta = output_dir + marker + "_hmm_purified.faa"
            write_new_fasta(bulk_marker_fasta, trimmed_hits_fasta)
        bulk_marker_fasta.clear()
    logging.info("done.\n")
    return hmmalign_input_fastas


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
            num_seqs, sequence_length = multiple_alignment_dimensions(seq_dict, msa_file)
            alignment_length_dict[msa_file] = sequence_length

            # Warn user if the multiple sequence alignment has grown significantly
            if verbosity and ref_seq_length*1.5 < sequence_length:
                logging.warning("Multiple alignment of " + denominator +
                                "caused >150% increase in the number of columns (" +
                                str(ref_seq_length) + '->' + str(sequence_length) + ").\n")
    return alignment_length_dict


def get_alignment_dims(treesapp_dir, marker_build_dict):
    alignment_dimensions_dict = dict()
    alignment_data_dir = os.sep.join([treesapp_dir, "data", 'alignment_data' + os.sep])
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
                    alignment_dimensions_dict[marker_code] = (multiple_alignment_dimensions(seq_dict, fasta))
    return alignment_dimensions_dict


def multiple_alignments(executables, treesapp_data_dir, output_dir, single_query_sequence_files, marker_build_dict,
                        tool="hmmalign", num_proc=4):
    """
    Wrapper function for the multiple alignment functions - only purpose is to make an easy decision at this point...

    :param executables: Dictionary mapping software names to their executables
    :param treesapp_data_dir:
    :param output_dir:
    :param single_query_sequence_files:
    :param marker_build_dict:
    :param tool: Tool to use for aligning query sequences to a reference multiple alignment [hmmalign|papara]
    :param num_proc: The number of alignment jobs to run in parallel
    :return: list of multiple sequence alignment files, generated by `tool`
    """
    if tool == "papara":
        concatenated_msa_files = prepare_and_run_papara(executables, treesapp_data_dir, single_query_sequence_files,
                                                        marker_build_dict, output_dir)
    elif tool == "hmmalign":
        hmm_dir = treesapp_data_dir + os.sep + "hmm_data"
        alignment_dir = treesapp_data_dir + os.sep + "alignment_data"
        concatenated_msa_files = prepare_and_run_hmmalign(executables, hmm_dir, alignment_dir,
                                                          single_query_sequence_files, marker_build_dict, num_proc)
    else:
        logging.error("Unrecognized tool '" + str(tool) + "' for multiple sequence alignment.\n")
        sys.exit(3)
    return concatenated_msa_files


def create_ref_phy_files(alignment_dir, output_dir, single_query_fasta_files, marker_build_dict, ref_aln_dimensions):
    """
    Creates a phy file for every reference marker that was matched by a query sequence
    :param alignment_dir:
    :param output_dir:
    :param single_query_fasta_files:
    :param marker_build_dict:
    :param ref_aln_dimensions:
    :return:
    """

    # Convert the reference sequence alignments to .phy files for every marker identified
    for query_fasta in single_query_fasta_files:
        marker = re.match("(.*)_hmm_purified.*", os.path.basename(query_fasta)).group(1)
        denominator = None
        for denominator in marker_build_dict.keys():
            if marker_build_dict[denominator].cog == marker:
                break

        ref_alignment_phy = output_dir + marker + ".phy"
        if os.path.isfile(ref_alignment_phy):
            continue
        aligned_fasta = alignment_dir + os.sep + marker + ".fa"

        num_ref_seqs, ref_align_len = ref_aln_dimensions[denominator]
        aligned_fasta_dict = read_fasta_to_dict(aligned_fasta)
        dict_for_phy = dict()
        for seq_name in aligned_fasta_dict:
            dict_for_phy[seq_name.split('_')[0]] = aligned_fasta_dict[seq_name]
        phy_dict = utilities.reformat_fasta_to_phy(dict_for_phy)

        utilities.write_phy_file(ref_alignment_phy, phy_dict, (num_ref_seqs, ref_align_len))
    return


def prepare_and_run_papara(executables, treesapp_resources, single_query_fasta_files, marker_build_dict, output_dir):
    """
    Uses the Parsimony-based Phylogeny-aware short Read Alignment (PaPaRa) tool

    :param executables:
    :param treesapp_resources:
    :param single_query_fasta_files:
    :param marker_build_dict:
    :param output_dir:
    :return:
    """
    query_alignment_files = dict()
    logging.info("Running PaPaRa... ")
    # TODO: Parallelize; PaPaRa's posix threading is not guaranteed with binary
    start_time = time.time()

    # Convert the reference sequence alignments to .phy files for every marker identified
    for query_fasta in sorted(single_query_fasta_files):
        file_name_info = re.match(r"(.*)_hmm_purified.*\.(f.*)$", os.path.basename(query_fasta))
        if file_name_info:
            marker, extension = file_name_info.groups()
        else:
            logging.error("Unable to parse information from file name:" + "\n" + str(query_fasta) + "\n")
            sys.exit(3)

        ref_marker = utilities.fish_refpkg_from_build_params(marker, marker_build_dict)
        query_multiple_alignment = re.sub('.' + re.escape(extension) + r"$", ".phy", query_fasta)
        tree_file = treesapp_resources + "tree_data" + os.sep + marker + "_tree.txt"
        ref_alignment_phy = output_dir + marker + ".phy"
        if not os.path.isfile(ref_alignment_phy):
            logging.error("Phylip file '" + ref_alignment_phy + "' not found.\n")
            sys.exit(3)

        wrapper.run_papara(executables["papara"], tree_file, ref_alignment_phy, query_fasta, ref_marker.molecule)
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


def prepare_and_run_hmmalign(execs, hmm_dir, alignment_dir, single_query_fasta_files, marker_build_dict, n_proc=2):
    """
    Runs `hmmalign` to add the query sequences into the reference FASTA multiple alignments

    :param execs:
    :param hmm_dir:
    :param alignment_dir:
    :param single_query_fasta_files:
    :param marker_build_dict:
    :param n_proc: The number of alignment jobs to run in parallel
    :return: list of multiple sequence alignment files (in Clustal format) generated by hmmalign.
    """

    hmmalign_singlehit_files = dict()
    mfa_out_dict = dict()
    logging.info("Running hmmalign... ")

    start_time = time.time()
    task_list = list()

    # Run hmmalign on each fasta file
    for query_fa_in in sorted(single_query_fasta_files):
        file_name_info = re.match(r"(.*)_hmm_purified.*\.(f.*)$", os.path.basename(query_fa_in))
        if file_name_info:
            marker, extension = file_name_info.groups()
        else:
            logging.error("Unable to parse information from file name:" + "\n" + str(query_fa_in) + "\n")
            sys.exit(3)

        query_mfa_out = re.sub('.' + re.escape(extension) + r"$", ".sto", query_fa_in)

        ref_marker = utilities.fish_refpkg_from_build_params(marker, marker_build_dict)
        if ref_marker.denominator not in hmmalign_singlehit_files:
            hmmalign_singlehit_files[ref_marker.denominator] = []
        try:
            mfa_out_dict[ref_marker.denominator].append(query_mfa_out)
        except KeyError:
            mfa_out_dict[ref_marker.denominator] = [query_mfa_out]

        # Get the paths to either the HMM or CM profile files
        ref_alignment = alignment_dir + os.sep + ref_marker.cog + ".fa"
        ref_profile = hmm_dir + os.sep + ref_marker.cog
        if ref_marker.kind == "phylogenetic_rRNA":
            ref_profile += ".cm"
        else:
            ref_profile += ".hmm"

        if ref_marker.kind == "phylogenetic_rRNA":
            task_list.append(wrapper.hmmalign_command(execs["cmalign"],
                                                      ref_alignment, ref_profile, query_fa_in, query_mfa_out))
        else:
            task_list.append(wrapper.hmmalign_command(execs["hmmalign"],
                                                      ref_alignment, ref_profile, query_fa_in, query_mfa_out))

    if len(task_list) > 0:
        cl_farmer = CommandLineFarmer("cmalign/hmmalign --mapali", n_proc)
        cl_farmer.add_tasks_to_queue(task_list)

        cl_farmer.task_queue.close()
        cl_farmer.task_queue.join()

    logging.info("done.\n")

    for refpkg_code in mfa_out_dict:
        for query_mfa_out in mfa_out_dict[refpkg_code]:
            mfa_file = re.sub(r"\.sto$", ".mfa", query_mfa_out)
            tmp_dict = read_stockholm_to_dict(query_mfa_out)
            seq_dict = dict()
            for seq_name in tmp_dict:
                seq_dict[seq_name.split('_')[0]] = tmp_dict[seq_name]
            write_new_fasta(seq_dict, mfa_file)
            hmmalign_singlehit_files[refpkg_code].append(mfa_file)

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


def check_for_removed_sequences(aln_dir, trimmed_msa_files: dict, msa_files: dict, marker_build_dict: dict, min_len=10):
    """
    Reads the multiple alignment files (either Phylip or FASTA formatted) and looks for both reference and query
    sequences that have been removed. Multiple alignment files are removed from `mfa_files` if:
        1. all query sequences were removed; a DEBUG message is issued
        2. at least one reference sequence was removed
    This quality-control function is necessary for placing short query sequences onto reference trees.

    :param aln_dir: Path to the directory containing TreeSAPP reference package multiple sequence alignments
    :param trimmed_msa_files:
    :param msa_files: A dictionary containing the untrimmed MSA files indexed by reference package code (denominator)
    :param marker_build_dict: A dictionary of MarkerBuild objects indexed by their refpkg codes/denominators
    :param min_len: The minimum allowable sequence length after trimming (not including gap characters)
    :return: A dictionary of denominators, with multiple alignment dictionaries as values. Example:
        {M0702: { "McrB_hmm_purified.phy-BMGE.fasta": {'1': seq1, '2': seq2}}}
    """
    qc_ma_dict = dict()
    num_successful_alignments = 0
    discarded_seqs_string = ""
    trimmed_away_seqs = dict()
    untrimmed_msa_failed = []
    logging.debug("Validating trimmed multiple sequence alignment files... ")

    for denominator in sorted(trimmed_msa_files.keys()):
        marker = marker_build_dict[denominator].cog
        trimmed_away_seqs[marker] = 0
        # Create a set of the reference sequence names
        ref_headers = get_headers(aln_dir + os.sep + marker + ".fa")
        unique_refs = set([re.sub('_' + re.escape(marker), '', x)[1:] for x in ref_headers])
        msa_passed, msa_failed, summary_str = validate_alignment_trimming(trimmed_msa_files[denominator], unique_refs,
                                                                          True, min_len)

        # Report the number of sequences that are removed by BMGE
        for trimmed_msa_file in trimmed_msa_files[denominator]:
            try:
                prefix, tool = re.search('(' + re.escape(marker) + r"_.*_group\d+)-(BMGE|trimAl).fasta$",
                                         os.path.basename(trimmed_msa_file)).groups()
            except TypeError:
                logging.error("Unexpected file name format for a trimmed MSA.\n")
                sys.exit(3)
            pair = ""
            for msa_file in msa_files[denominator]:
                if re.search(re.escape(prefix) + '\.', msa_file):
                    pair = msa_file
                    break
            if pair:
                if trimmed_msa_file in msa_failed:
                    untrimmed_msa_failed.append(pair)
                trimmed_away_seqs[marker] += len(set(get_headers(pair)).difference(set(get_headers(trimmed_msa_file))))
            else:
                logging.error("Unable to map trimmed MSA file '" + trimmed_msa_file + "' to its original MSA.\n")
                sys.exit(5)

        if len(msa_failed) > 0:
            if len(untrimmed_msa_failed) != len(msa_failed):
                logging.error("Not all of the failed (" + str(len(msa_failed)) + '/' +
                              str(len(trimmed_msa_files[denominator])) +
                              "), trimmed MSA files were mapped to their original MSAs.\n")
                sys.exit(3)
            untrimmed_msa_passed, _, _ = validate_alignment_trimming(untrimmed_msa_failed, unique_refs,
                                                                     True, min_len)
            msa_passed.update(untrimmed_msa_passed)
        num_successful_alignments += len(msa_passed)
        qc_ma_dict[denominator] = msa_passed
        discarded_seqs_string += summary_str
        untrimmed_msa_failed.clear()

    logging.debug("done.\n")
    logging.debug("\tSequences removed during trimming:\n\t\t" +
                  '\n\t\t'.join([k + ": " + str(trimmed_away_seqs[k]) for k in trimmed_away_seqs.keys()]) + "\n")

    logging.debug("\tSequences <" + str(min_len) + " characters removed after trimming:" +
                  discarded_seqs_string + "\n")

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
            num_seqs, trimmed_seq_length = multiple_alignment_dimensions(multi_align, multi_align_file)

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
        n_trimmed_files = len(trimmed_length_dict[denominator])
        if n_trimmed_files > 0:
            trimming_performance_string += str(round(sum(trimmed_length_dict[denominator])/n_trimmed_files, 1)) + "\n"
        else:
            trimming_performance_string += str(0.0) + "\n"

    logging.debug(trimming_performance_string + "\n")
    return


def produce_phy_files(qc_ma_dict, molecule_type="prot"):
    """
    Produces phy files from the provided list of alignment files

    :param qc_ma_dict:
    :param molecule_type:
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

                if molecule_type != "prot":
                    sequence = re.sub('U', 'T', sequence)  # Got error from RAxML when encountering Uracil

                sequences_for_phy[seq_name] = sequence

            # Write the sequences to the phy file
            phy_dict = utilities.reformat_fasta_to_phy(sequences_for_phy)
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

    insertion_point_node_hash = utilities.Autovivify()
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
    tree_elements = utilities.Autovivify()
    
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
    list_of_already_used_attachments = utilities.Autovivify()
    rooted_trees = utilities.Autovivify()
    
    for node in sorted(tree_info['quartets'].keys(), key=int):
        if node in list_of_already_used_attachments:
            continue
        for attachment in sorted(tree_info['quartets'][node].keys(), key=int):
            list_of_already_used_attachments[attachment] = 1
            tree_string = ''
            node_infos = utilities.Autovivify()
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
        node_infos2 = utilities.Autovivify()
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
    terminal_children_strings_of_assignments = utilities.Autovivify()

    for assignment in sorted(assignments.keys()):
        internal_node_of_assignment = insertion_point_node_hash[assignment]
        # parse_log.write("Starting to retrieve all subtrees at " + time.ctime())
        rooted_tree_nodes = parallel_subtree_node_retriever(rooted_trees, num_threads, parse_log)
        # parse_log.write("Finished retrieving subtrees at " + time.ctime() + "\n")
        for rooted_tree_info in rooted_tree_nodes:
            assignment_subtree = str(rooted_tree_info['subtree_of_node'][str(internal_node_of_assignment)])
            terminal_children = utilities.Autovivify()

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
    terminal_children_strings_of_reference = utilities.Autovivify()

    for node in sorted(reference_tree_info['subtree_of_node'].keys()):
        reference_subtree = reference_tree_info['subtree_of_node'][node]
        terminal_children = utilities.Autovivify()
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
    real_terminal_children_of_assignments = utilities.Autovivify()
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


def delete_files(clean_up, output_dir_var, section):
    files_to_be_deleted = []
    if clean_up:
        if section == 1:
            files_to_be_deleted += glob.glob(output_dir_var + '*BLAST_results*')
            files_to_be_deleted += glob.glob(output_dir_var + '*blast_result_purified.txt')
        if section == 2:
            files_to_be_deleted += glob.glob(output_dir_var + '*_sequence.txt')
            files_to_be_deleted += glob.glob(output_dir_var + '*sequence_shortened.txt')
        if section == 3:
            files_to_be_deleted += glob.glob(output_dir_var + '*genewise.txt')
            files_to_be_deleted += glob.glob(output_dir_var + '*genewise_result_summary.txt')
            files_to_be_deleted += glob.glob(output_dir_var + '*rRNA_result_summary.txt')
        if section == 4:
            files_to_be_deleted += glob.glob(output_dir_var + '*.mfa')
            files_to_be_deleted += glob.glob(output_dir_var + '*.mfa-gb')
            files_to_be_deleted += glob.glob(output_dir_var + '*.mfa-gb.txt')
            files_to_be_deleted += glob.glob(output_dir_var + "*.sam")
        if section == 5:
            files_to_be_deleted += glob.glob(output_dir_var + '*_exit_after_trimal.txt')
            files_to_be_deleted += glob.glob(output_dir_var + '*.mfa-trimal')
            files_to_be_deleted += glob.glob(output_dir_var + '*.cl')
            files_to_be_deleted += glob.glob(output_dir_var + '*_RAxML.txt')
            files_to_be_deleted += glob.glob(output_dir_var + 'RAxML_entropy.*')
            files_to_be_deleted += glob.glob(output_dir_var + '*RAxML_info.txt')
            files_to_be_deleted += glob.glob(output_dir_var + '*RAxML_labelledTree.txt')
            files_to_be_deleted += glob.glob(output_dir_var + 'RAxML_classificationLikelihoodWeights*')
            files_to_be_deleted += glob.glob(output_dir_var + '*.phy')
            files_to_be_deleted += glob.glob(output_dir_var + '*.phy.reduced')
            files_to_be_deleted += glob.glob(output_dir_var + '*RAxML_classification.txt')
            # Need this for annotate_extra_treesapp.py
            # files_to_be_deleted += glob.glob(output_dir_var + '*.jplace')

    for useless_file in files_to_be_deleted:
        if path.exists(useless_file):
            os.remove(useless_file)


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


def filter_short_sequences(aa_dictionary, length_threshold):
    """
    Removes all sequences shorter than length_threshold from a dictionary
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


def align_reads_to_nucs(bwa_exe: str, reference_fasta: str, aln_output_dir: str,
                        reads: str, pairing: str, reverse=None, num_threads=2) -> str:
    """
    Align the predicted ORFs to the reads using BWA MEM

    :param bwa_exe: Path to the BWA executable
    :param reference_fasta: A FASTA file containing the sequences to be aligned to
    :param aln_output_dir: Path to the directory to write the index and SAM files
    :param reads: FASTQ file containing reads to be aligned to the reference FASTA file
    :param pairing: Either 'se' or 'pe' indicating the reads are single-end or paired-end, respectively
    :param reverse: Path to reverse-orientation mate pair reads [OPTIONAL]
    :param num_threads: Number of threads for BWA MEM to use
    :return: Path to the SAM file
    """
    if not os.path.exists(aln_output_dir):
        try:
            os.makedirs(aln_output_dir)
        except OSError:
            if os.path.exists(aln_output_dir):
                logging.warning("Overwriting files in " + aln_output_dir + ".\n")
            else:
                raise OSError("Unable to make " + aln_output_dir + "!\n")

    logging.info("Aligning reads to ORFs with BWA MEM... ")

    sam_file = aln_output_dir + '.'.join(os.path.basename(reference_fasta).split('.')[0:-1]) + ".sam"
    if os.path.isfile(sam_file):
        logging.info("output found.\n")
        return sam_file
    index_command = [bwa_exe, "index"]
    index_command += [reference_fasta]
    index_command += ["1>", "/dev/null", "2>", aln_output_dir + "treesapp_bwa_index.stderr"]

    launch_write_command(index_command)

    bwa_command = [bwa_exe, "mem"]
    bwa_command += ["-t", str(num_threads)]
    if pairing == "pe" and not reverse:
        bwa_command.append("-p")
        logging.debug("FASTQ file containing reverse mates was not provided - assuming the reads are interleaved!\n")
    elif pairing == "se":
        bwa_command += ["-S", "-P"]

    bwa_command.append(reference_fasta)
    bwa_command.append(reads)
    if pairing == "pe" and reverse:
        bwa_command.append(reverse)
    bwa_command += ["1>", sam_file, "2>", aln_output_dir + "treesapp_bwa_mem.stderr"]

    p_bwa = subprocess.Popen(' '.join(bwa_command), shell=True, preexec_fn=os.setsid)
    p_bwa.wait()
    if p_bwa.returncode != 0:
        logging.error("bwa mem did not complete successfully for:\n" +
                      str(' '.join(bwa_command)) + "\n")

    logging.info("done.\n")

    return sam_file


def summarize_placements_rpkm(tree_saps: dict, abundance_dict: dict, marker_build_dict: dict, final_output_dir: str):
    """
    Recalculates the percentages for each marker gene final output based on the RPKM values
    The abundance_dict contains RPKM values of contigs whereas tree_saps may be fragments of contigs,
    and if multiple fragments are classified this could "inflate" the RPKM values. Currently, this is not handled.

    :param tree_saps: A dictionary of ItolJplace instances, indexed by their respective RefPkg codes (denominators)
    :param abundance_dict: A dictionary mapping predicted (not necessarily classified) seq_names to abundance values
    :param marker_build_dict:
    :param final_output_dir:
    :return: None
    """
    placement_rpkm_map = dict()  # Used to map the internal nodes to the total RPKM of all descendent nodes
    marker_rpkm_map = dict()  # Used to hold the RPKM sums for each marker
    orf_rpkms = dict()  # Essentially a duplicate of abundance_dict after the first for-loop

    # Identify the internal node each sequence was placed, used for iTOL outputs
    for denominator in tree_saps:
        for placed_sequence in tree_saps[denominator]:  # type ItolJplace
            if not placed_sequence.classified:
                continue
            seq_name = re.sub(r"\|{0}\|\d+_\d+$".format(placed_sequence.name), '', placed_sequence.contig_name)
            try:
                placed_sequence.abundance = abundance_dict.pop(seq_name)
                orf_rpkms[seq_name] = placed_sequence.abundance
            except KeyError:
                if seq_name in orf_rpkms:
                    placed_sequence.abundance = orf_rpkms[seq_name]
                else:
                    logging.error("Unable to find sequence '" + seq_name +
                                  "' in RPKM abundance dictionary keys. Examples:\n" +
                                  "\n".join("'" + f + "'" for f in list(abundance_dict.keys())[0:6]) + "\n")
                    sys.exit(3)
            if placed_sequence.inode not in placement_rpkm_map:
                placement_rpkm_map[placed_sequence.inode] = 0

    if len(abundance_dict) > 0:
        logging.warning(str(len(abundance_dict)) + " sequence names remain in the RPKM abundance dictionary.\n")
        logging.debug("Leftover sequences in abundance dict:\n" + "\n".join(abundance_dict.keys()) + "\n")
    orf_rpkms.clear()

    # Calculate the percentage contribution of each placed sequence
    for denominator in tree_saps:
        marker = marker_build_dict[denominator].cog
        marker_rpkm_total = 0
        marker_rpkm_map[marker] = dict()
        for placed_sequence in tree_saps[denominator]:
            if not placed_sequence.classified:
                continue
            placement_rpkm_map[placed_sequence.inode] += float(placed_sequence.abundance)
            marker_rpkm_total += float(placed_sequence.abundance)
            marker_rpkm_map[marker][placed_sequence.inode] = 0
        for placement in marker_rpkm_map[marker]:
            try:
                percentage = (placement_rpkm_map[placement]*100)/marker_rpkm_total
            except ZeroDivisionError:
                percentage = 0
            marker_rpkm_map[marker][placement] = percentage

    # TODO: currently doesn't work as required files are missing. Fix or abandon?
    for marker in marker_rpkm_map:
        ref_marker = utilities.fish_refpkg_from_build_params(marker, marker_build_dict)

        final_output_file = final_output_dir + str(ref_marker.denominator) + "_concatenated_RAxML_outputs.txt"
        # Not all of the genes predicted will have made it to the RAxML stage
        if os.path.isfile(final_output_file):
            shutil.move(final_output_file, final_output_dir + ref_marker.denominator + "_concatenated_counts.txt")
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


def abundify_tree_saps(tree_saps: dict, abundance_dict: dict):
    """
    Add abundance (RPKM or presence count) values to the TreeProtein instances (abundance variable)

    :param tree_saps: Dictionary mapping refpkg codes to all TreeProtein instances for classified sequences
    :param abundance_dict: Dictionary mapping sequence names to floats
    :return: None
    """
    abundance_mapped_acc = 0
    for refpkg_code in tree_saps:
        for placed_seq in tree_saps[refpkg_code]:  # type: TreeProtein
            if not placed_seq.abundance:
                # Filter out RPKMs for contigs not associated with the target marker
                try:
                    placed_seq.abundance = abundance_dict[placed_seq.contig_name]
                    abundance_mapped_acc += 1
                except KeyError:
                    placed_seq.abundance = 0.0
            else:
                abundance_mapped_acc += 1

    if abundance_mapped_acc == 0:
        logging.warning("No placed sequences with abundances identified.\n")

    return


def create_itol_labels(marker, itol_base_dir, tree_data_dir):
    """
    Create the marker_labels.txt file for each marker gene that was used for classification

    :param marker: Name of the reference package
    :param itol_base_dir: Path to the directory for where these outputs should be written to
    :param tree_data_dir: Path to the directory containing reference package trees for TreeSAPP
    :return: None
    """
    itol_label_file = itol_base_dir + os.sep + marker + os.sep + marker + "_labels.txt"
    tax_ids_file = tree_data_dir + os.sep + "tax_ids_" + marker + ".txt"

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
            logging.error("split(\'\\t\') on " + str(line) + " generated " + str(len(line.split("\t"))) + " fields.")
            sys.exit(9)
        if len(fields) == 2:
            number, translation = fields
        elif len(fields) == 3:
            number, translation, lineage = fields
        else:
            logging.error("Unexpected number of fields in " + tax_ids_file +
                          ".\nInvoked .split(\'\\t\') on line " + str(line))
            sys.exit(3)
        label_f.write(number + ',' + translation + "\n")

    tax_ids.close()
    label_f.close()

    return


def generate_simplebar(target_marker, tree_protein_list, itol_bar_file):
    """
    From the basic RPKM output csv file, generate an iTOL-compatible simple bar-graph file for each leaf

    :param target_marker:
    :param tree_protein_list: A list of TreeProtein objects, for single sequences
    :param itol_bar_file: The name of the file to write the simple-bar data for iTOL
    :return:
    """
    leaf_rpkm_sums = dict()

    for tree_sap in tree_protein_list:
        if tree_sap.name == target_marker and tree_sap.classified:
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


def filter_placements(tree_saps, marker_build_dict, tree_data_dir: str, min_likelihood: float):
    """
    Determines the total distance of each placement from its branch point on the tree
    and removes the placement if the distance is deemed too great

    :param tree_saps: A dictionary containing TreeProtein objects
    :param marker_build_dict: A dictionary of MarkerBuild objects (used here for lowest_confident_rank)
    :param tree_data_dir: Directory containing reference package tree files (Newick)
    :param min_likelihood: Likelihood-weight-ratio (LWR) threshold for filtering pqueries
    :return:
    """

    logging.info("Filtering low-quality placements... ")
    unclassified_seqs = dict()  # A dictionary tracking the seqs unclassified for each marker

    for denominator in tree_saps:
        marker = marker_build_dict[denominator].cog
        unclassified_seqs[marker] = dict()
        unclassified_seqs[marker]["low_lwr"] = list()
        unclassified_seqs[marker]["np"] = list()
        unclassified_seqs[marker]["beyond"] = list()
        unclassified_seqs[marker]["far_beyond"] = list()

        tree = Tree(tree_data_dir + os.sep + marker + "_tree.txt")
        max_dist, leaf_ds = tree_leaf_distances(tree)
        # Find the maximum distance and standard deviation of distances from the root to all leaves
        max_dist_threshold = max_dist
        mean_dist_threshold = utilities.mean(leaf_ds)
        logging.debug(denominator + " maximum pendant length distance threshold: " + str(max_dist_threshold) + "\n")

        for tree_sap in tree_saps[denominator]:
            tree_sap.filter_min_weight_threshold(min_likelihood)
            if not tree_sap.classified:
                unclassified_seqs[marker]["low_lwr"].append(tree_sap)
                continue
            if not tree_sap.placements:
                unclassified_seqs[tree_sap.name]["np"].append(tree_sap)
                continue
            elif len(tree_sap.placements) > 1:
                logging.warning("More than one placement for a single contig:\n" +
                                tree_sap.summarize())
                tree_sap.classified = False
                continue
            elif tree_sap.placements[0] == '{}':
                unclassified_seqs[marker]["np"].append(tree_sap)
                tree_sap.classified = False
                continue

            pendant_length = float(tree_sap.get_jplace_element("pendant_length"))
            # Discard this placement as a false positive if the pendant_length exceeds max_dist_threshold
            if pendant_length > max_dist_threshold:
                unclassified_seqs[tree_sap.name]["far_beyond"].append(tree_sap)
                tree_sap.classified = False
                continue
            tree_sap.lwr = float(tree_sap.get_jplace_element("like_weight_ratio"))
            if pendant_length > mean_dist_threshold and tree_sap.lwr < 0.66:
                unclassified_seqs[marker]["beyond"].append(tree_sap)
                tree_sap.classified = False
                continue

            # Find the length of the edge this sequence was placed onto
            leaf_children = tree_sap.node_map[int(tree_sap.inode)]
            distal_length = float(tree_sap.get_jplace_element("distal_length"))
            # Find the distance away from this edge's bifurcation (if internal) or tip (if leaf)

            if len(leaf_children) > 1:
                # We need to find the LCA in the Tree instance to find the distances to tips for ete3
                parent = tree.get_common_ancestor(leaf_children)
                tip_distances = parent_to_tip_distances(parent, leaf_children)
            else:
                tip_distances = [0.0]

            tree_sap.avg_evo_dist = round(distal_length + pendant_length + (sum(tip_distances) / len(tip_distances)), 4)
            tree_sap.distances = str(distal_length) + ',' +\
                                 str(pendant_length) + ',' +\
                                 str(sum(tip_distances) / len(tip_distances))

    logging.info("done.\n")

    declass_summary = ""
    for marker in unclassified_seqs:
        # unclassified_counts[marker] will always be >= distant_seqs[marker]
        for declass in unclassified_seqs[marker]:
            declass_summary += marker + '\t' + declass + '\t' + str(len(unclassified_seqs[marker][declass])) + "\n"

    logging.debug(declass_summary)

    return tree_saps


def write_tabular_output(tree_saps, tree_numbers_translation, marker_build_dict, sample_name, output_file):
    """
    Write the final classification table

    :param tree_saps: A dictionary containing TreeProtein objects
    :param tree_numbers_translation: Dictionary containing taxonomic information for each leaf in the reference tree
    :param marker_build_dict: A dictionary of MarkerBuild objects (used here for lowest_confident_rank)
    :param sample_name: String representing the name of the sample (i.e. Assign.sample_prefix)
    :param output_file: Path to write the classification table
    :return: None
    """
    leaf_taxa_map = dict()
    # TODO: Add the start and stop positions of the extracted sequence to the classification table
    tab_out_string = "Sample\tQuery\tMarker\tLength\tTaxonomy\tConfident_Taxonomy\tAbundance\tiNode\tLWR\tEvoDist\tDistances\n"

    for denominator in tree_saps:
        # All the leaves for that tree [number, translation, lineage]
        leaves = tree_numbers_translation[denominator]
        lineage_list = list()
        # Test if the reference set have lineage information
        for leaf in leaves:
            lineage_list.append(utilities.clean_lineage_string(leaf.lineage).split('; '))
            leaf_taxa_map[leaf.number] = leaf.lineage
        taxonomic_counts = enumerate_taxonomic_lineages(lineage_list)

        for tree_sap in tree_saps[denominator]:  # type: TreeProtein
            if not tree_sap.classified:
                continue

            if tree_sap.placements and not tree_sap.recommended_lineage:
                tree_sap.lineage_list = children_lineage(leaf_taxa_map, tree_sap.placements[0], tree_sap.node_map)

                if len(tree_sap.lineage_list) == 0:
                    logging.error("Unable to find lineage information for marker " +
                                  denominator + ", contig " + tree_sap.contig_name + "!\n")
                    sys.exit(3)
                elif len(tree_sap.lineage_list) == 1:
                    tree_sap.lct = tree_sap.lineage_list[0]
                    tree_sap.wtd = 0.0
                else:
                    lca = tree_sap.megan_lca()
                    # algorithm options are "MEGAN", "LCAp", and "LCA*" (default)
                    tree_sap.lct = lowest_common_taxonomy(tree_sap.lineage_list, lca, taxonomic_counts, "LCA*")
                    tree_sap.wtd, status = weighted_taxonomic_distance(tree_sap.lineage_list, tree_sap.lct)
                    if status > 0:
                        tree_sap.summarize()

                # Based on the calculated distance from the leaves, what rank is most appropriate?
                recommended_rank = rank_recommender(tree_sap.avg_evo_dist,
                                                    marker_build_dict[denominator].pfit)
                if tree_sap.lct.split("; ")[0] != "Root":
                    tree_sap.lct = "Root; " + tree_sap.lct
                    recommended_rank += 1
                tree_sap.recommended_lineage = tree_sap.lowest_confident_taxonomy(recommended_rank)
                # tree_sap.summarize()

            tab_out_string += '\t'.join([sample_name,
                                         re.sub(r"\|{0}\|\d+_\d+$".format(tree_sap.name), '', tree_sap.contig_name),
                                         tree_sap.name,
                                         str(tree_sap.seq_len),
                                         utilities.clean_lineage_string(tree_sap.lct),
                                         tree_sap.recommended_lineage,
                                         str(tree_sap.abundance),
                                         str(tree_sap.inode),
                                         str(tree_sap.lwr),
                                         str(tree_sap.avg_evo_dist),
                                         tree_sap.distances]) + "\n"
    try:
        tab_out = open(output_file, 'w')
    except IOError:
        logging.error("Unable to open " + output_file + " for writing!\n")
        sys.exit(3)

    tab_out.write(tab_out_string)
    tab_out.close()

    return


def parse_raxml_output(epa_output_dir, tree_data_dir, marker_build_dict, parsing_method="best"):
    """

    :param epa_output_dir: Directory where EPA wrote the JPlace files
    :param tree_data_dir: Directory containing reference package tree files (Newick)
    :param marker_build_dict: Dictionary of MarkerBuild instances indexed by denominator (refpkg code e.g. M0701)
    :param parsing_method: Method for parsing the final placement edge position [best | lca]
    :return:
        1. Dictionary of TreeProtein instances indexed by denominator (refpkg code e.g. M0701)
        2. Dictionary of an ItolJplace instance (values) mapped to marker name
        3. Dictionary mapping the number of sequences that were unclassified (value) to a marker (key)
    """

    logging.info('Parsing the RAxML outputs... ')

    function_start_time = time.time()

    jplace_files = glob.glob(epa_output_dir + '*.jplace')
    jplace_collection = organize_jplace_files(jplace_files)
    itol_data = dict()  # contains all pqueries, indexed by marker name (e.g. McrA, nosZ, 16srRNA)
    tree_saps = dict()  # contains individual pquery information for each mapped protein (N==1), indexed by denominator
    # Use the jplace files to guide which markers iTOL outputs should be created for
    classified_seqs = 0
    for denominator in jplace_collection:
        marker = marker_build_dict[denominator].cog
        if denominator not in tree_saps:
            tree_saps[denominator] = list()
        for filename in jplace_collection[denominator]:
            # Load the JSON placement (jplace) file containing >= 1 pquery into ItolJplace object
            jplace_data = jplace_parser(filename)
            tree_index = index_tree_edges(jplace_data.tree)
            # Demultiplex all pqueries in jplace_data into individual TreeProtein objects
            tree_placement_queries = demultiplex_pqueries(jplace_data)
            # Filter the placements, determine the likelihood associated with the harmonized placement
            for pquery in tree_placement_queries:
                pquery.name = marker
                seq_info = re.match(r"(.*)\|" + re.escape(marker) + r"\|(\\d+)_(\\d+)$", pquery.contig_name)
                if seq_info:
                    pquery.contig_name = seq_info.group(1)
                    start, end = seq_info.groups()[1:]
                    pquery.seq_len = int(end) - int(start)
                pquery.create_jplace_node_map()
                pquery.check_jplace(tree_index)
                if parsing_method == "best":
                    pquery.filter_max_weight_placement()
                else:
                    pquery.harmonize_placements(tree_data_dir)
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
                  "\t" + str(classified_seqs) + " sequences placed into trees by RAxML.\n\n")

    return tree_saps, itol_data


def produce_itol_inputs(tree_saps, marker_build_dict, itol_data, output_dir: str, treesapp_data_dir: str):
    """
    Function to create outputs for the interactive tree of life (iTOL) webservice.
    There is a directory for each of the marker genes detected to allow the user to "drag-and-drop" all files easily

    :param tree_saps:
    :param marker_build_dict:
    :param itol_data:
    :param output_dir:
    :param treesapp_data_dir:
    :return: None
    """
    logging.info("Generating inputs for iTOL... ")

    itol_base_dir = output_dir + 'iTOL_output' + os.sep
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
        marker_placements = itol_data[marker]
        marker_treesaps = tree_saps[denominator]

        bipartition_file = os.sep.join([treesapp_data_dir, "tree_data", marker + "_bipartitions.txt"])
        if os.path.isfile(bipartition_file):
            marker_placements = add_bipartitions(marker_placements, bipartition_file)

        # Make a master jplace file from the set of placements in all jplace files for each marker
        master_jplace = itol_base_dir + marker + os.sep + marker + "_complete_profile.jplace"
        marker_placements = filter_jplace_data(marker_placements, marker_treesaps)
        # TODO: validate no distal lengths exceed their corresponding edge lengths

        write_jplace(marker_placements, master_jplace)
        itol_data[marker].clear_object()
        marker_placements.clear_object()
        # Create a labels file from the tax_ids_marker.txt
        create_itol_labels(marker, itol_base_dir, treesapp_data_dir + os.sep + "tree_data")

        annotation_style_files = glob.glob(os.sep.join([treesapp_data_dir, "iTOL_data", marker + "*"]))
        # Copy the respective colours and styles files for each marker found to the itol_output directories
        colours_styles = os.sep.join([treesapp_data_dir, "iTOL_data", marker + "_colours_style.txt"])
        colour_strip = os.sep.join([treesapp_data_dir, "iTOL_data", marker + "_colour_strip.txt"])
        if colours_styles not in annotation_style_files:
            style_missing.append(marker)
        if colour_strip not in annotation_style_files:
            strip_missing.append(marker)

        for annotation_file in annotation_style_files:
            shutil.copy(annotation_file, itol_base_dir + marker)
        itol_bar_file = os.sep.join([output_dir, "iTOL_output", marker, marker + "_abundance_simplebar.txt"])
        generate_simplebar(marker, marker_treesaps, itol_bar_file)

    logging.info("done.\n")
    if style_missing:
        logging.debug("A colours_style.txt file does not yet exist for markers:\n\t" +
                      "\n\t".join(style_missing) + "\n")
    if strip_missing:
        logging.debug("A colours_strip.txt file does not yet exist for markers:\n\t" +
                      "\n\t".join(strip_missing) + "\n")

    return

