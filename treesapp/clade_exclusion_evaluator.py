#!/usr/bin/env python3

__author__ = 'Connor Morgan-Lang'

import argparse
import sys
import os
import inspect
import shutil
from ete3 import Tree
from glob import glob
from random import randint

cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile(inspect.currentframe()))[0]))
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)
sys.path.insert(0, cmd_folder + os.sep + ".." + os.sep)
from fasta import format_read_fasta, write_new_fasta, get_headers, read_fasta_to_dict
from create_treesapp_refpkg import get_header_format, finalize_ref_seq_lineages
from utilities import return_sequence_info_groups, find_executables
from external_command_interface import launch_write_command
from file_parsers import parse_ref_build_params, tax_ids_file_to_leaves,\
    read_graftm_classifications, read_marker_classification_table, parse_assignments
from classy import prep_logging, get_header_info, register_headers, MarkerTest
from entrez_utils import *
from phylo_dist import trim_lineages_to_rank
from lca_calculations import all_possible_assignments, grab_graftm_taxa, identify_excluded_clade

# TODO: Ensure this dictionary works for every taxonomic hierarchy scheme
_RANK_DEPTH_MAP = {0: "Cellular organisms", 1: "Kingdom",
                   2: "Phylum", 3: "Class", 4: "Order",
                   5: "Family", 6: "Genus", 7: "Species", 8: "Strain"}


def get_arguments():
    parser = argparse.ArgumentParser(add_help=False)
    required_args = parser.add_argument_group("Required arguments")
    required_args.add_argument('-o', "--output",
                               required=True,
                               help='Path to a directory for writing output files')
    required_args.add_argument('-i', '--fasta_input',
                               help='Your sequence input file (for TreeSAPP) in FASTA format',
                               required=True)
    required_args.add_argument("-r", "--reference_marker",
                               help="Short-form name of the marker gene to be tested (e.g. mcrA, pmoA, nosZ)",
                               required=True)

    optopt = parser.add_argument_group("Optional options")
    optopt.add_argument("--fresh", default=False, required=False, action="store_true",
                        help="Recalculate a fresh phylogenetic tree with the target clades removed instead of"
                             " removing the leaves corresponding to targets from the reference tree.")
    optopt.add_argument("--tool", default="treesapp", required=False,
                        choices=["treesapp", "graftm", "diamond"],
                        help="Classify using one of the tools: treesapp [DEFAULT], graftm, or diamond.")
    optopt.add_argument("-t", "--taxon_rank",
                        help="Comma-separated list of the taxonomic ranks to test " +
                        "choices = [Phylum, Class, Order, Family, Genus, Species] " +
                        "(DEFAULT = Species)",
                        default="Species",
                        required=False)
    optopt.add_argument('-m', '--molecule',
                        help='the type of input sequences (prot = Protein [DEFAULT]; dna = Nucleotide; rrna = rRNA)',
                        default='prot',
                        choices=['prot', 'dna', 'rrna'])
    optopt.add_argument("-l", "--length",
                        required=False, type=int, default=0,
                        help="Arbitrarily slice the input sequences to this length. "
                             "Useful for testing classification accuracy for fragments.")

    miscellaneous_opts = parser.add_argument_group("Miscellaneous options")
    miscellaneous_opts.add_argument('--overwrite', action='store_true', default=False,
                                    help='overwrites previously processed output folders')
    miscellaneous_opts.add_argument("-T", "--threads", required=False, default=4, type=int,
                                    help="The number of threads to be used for classification.")
    miscellaneous_opts.add_argument('-v', '--verbose', action='store_true', default=False,
                                    help='Prints a more verbose runtime log')
    miscellaneous_opts.add_argument("-h", "--help",
                                    action="help",
                                    help="Show this help message and exit")

    args = parser.parse_args()

    if args.output[-1] != os.sep:
        args.output += os.sep
    args.treesapp = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + os.sep + ".." + os.sep

    if args.overwrite:
        if os.path.exists(args.output):
            shutil.rmtree(args.output)

    args.min_seq_length = 1

    if sys.version_info > (2, 9):
        args.py_version = 3
    else:
        args.py_version = 2

    return args


def parse_distances(classification_lines):
    distances = dict()
    distances["distal"] = list()
    distances["pendant"] = list()
    distances["tip"] = list()
    for fields in classification_lines:
        dist_fields = fields[-1].split(',')
        distances["distal"].append(float(dist_fields[0]))
        distances["pendant"].append(float(dist_fields[1]))
        distances["tip"].append(float(dist_fields[2]))
    return distances


def read_intermediate_assignments(inter_class_file):
    logging.debug("Reading " + inter_class_file + " for saved intermediate data... ")
    assignments = dict()
    n_saved = 0
    try:
        file_handler = open(inter_class_file, 'r')
    except OSError:
        logging.error("Unable to open " + inter_class_file + " for reading.\n")
        sys.exit(21)
    line = file_handler.readline()
    while line:
        line = line.strip()
        try:
            marker, ref, queries = line.split('\t')
        except ValueError:
            logging.warning("\tWARNING: " + inter_class_file + " is incorrectly formatted so regenerating instead.\n" +
                            "Violating line with " + str(len(line.split('\t'))) + " elements:\n" + line + "\n")
            return assignments, n_saved
        queries_list = queries.split(',')
        n_saved += len(queries_list)

        if marker not in assignments:
            assignments[marker] = dict()
        assignments[marker][ref] = queries_list
        line = file_handler.readline()

    file_handler.close()
    logging.debug("done.\n")
    return assignments, n_saved


def write_intermediate_assignments(inter_class_file, assignments):
    logging.debug("Saving intermediate data... ")
    try:
        file_handler = open(inter_class_file, 'w')
    except OSError:
        logging.error("Unable to open " + inter_class_file + " for writing!\n")
        sys.exit(21)
    assignments_string = ""
    for marker in assignments.keys():
        for ref in assignments[marker].keys():
            queries = assignments[marker][ref]
            if len(queries) == 0:
                logging.warning("No lineage information was downloaded for queries assigned to " +
                                ref + "\n")
            else:
                cleaned_ref = clean_lineage_string(ref)
                if not cleaned_ref:
                    logging.debug("Reference assignment is empty after cleaning the lineage:\n" +
                                  ref + "\n")
                else:
                    assignments_string += marker + "\t" + cleaned_ref + "\t"
                    assignments_string += ','.join(queries) + "\n"

    file_handler.write(assignments_string)
    file_handler.close()
    logging.debug("done.\n")
    return


def determine_offset(classified, optimal):
    # Figure out which taxonomic lineage is longer
    offset = 0
    while classified != optimal and offset < 7:
        offset += 1
        classified_ranks = classified.split("; ")
        optimal_ranks = optimal.split("; ")
        if len(classified_ranks) > len(optimal_ranks):
            classified = "; ".join(classified_ranks[:-1])
        elif len(classified_ranks) < len(optimal_ranks):
            optimal = "; ".join(optimal_ranks[:-1])
        else:
            optimal = "; ".join(optimal_ranks[:-1])
            classified = "; ".join(classified_ranks[:-1])
    return offset


def summarize_taxonomic_diversity(marker_eval_instance):
    """
    Function for summarizing the taxonomic diversity of a reference dataset by rank

    :param marker_eval_instance: An instance of the MarkerTest class.
        The following instance variables need to be populated:
        taxa_tests: list of TaxonTest objects, resulting from a GraftM or TreeSAPP analysis
    :return:
    """
    depth = 1  # Accumulator for parsing _RANK_DEPTH_MAP; not really interested in Cellular Organisms or Strains.
    info_str = ""
    while depth < 8:
        rank = _RANK_DEPTH_MAP[depth]
        unique_taxa = marker_eval_instance.get_unique_taxa_tested(rank)
        if unique_taxa:
            buffer = " "
            while len(rank) + len(str(len(unique_taxa))) + len(buffer) < 12:
                buffer += ' '
            info_str += "\t" + rank + buffer + str(len(unique_taxa)) + "\n"
        else:
            pass
        depth += 1
    logging.info("Number of unique lineages tested:\n" + info_str)
    return


def get_classification_performance(marker_eval_instance):
    """
    Correct if: optimal_assignment == query_lineage

    :param marker_eval_instance: An instance of the MarkerTest class.
        The following instance variables need to be populated:
        marker: Name of the marker currently being evaluated (e.g., nifHc1, mcrA)
        classifications: The rank-wise classification dictionary
    :return:
    """
    clade_exclusion_tabular_string = ""
    std_out_report_string = ""
    clade_exclusion_strings = list()
    rank_assigned_dict = marker_eval_instance.classifications

    sys.stdout.write("Rank-level performance of " + marker_eval_instance.target_marker + ":\n")
    sys.stdout.write("\tRank\tQueries\tClassified\tCorrect\tD=1\tD=2\tD=3\tD=4\tD=5\tD=6\tD=7\n")

    for depth in sorted(_RANK_DEPTH_MAP):
        rank = _RANK_DEPTH_MAP[depth]
        if rank == "Cellular organisms":
            continue
        correct = 0
        incorrect = 0
        taxonomic_distance = dict()
        n_queries, n_classified, sensitivity = marker_eval_instance.get_sensitivity(rank)
        for dist in range(0, 8):
            taxonomic_distance[dist] = 0
        std_out_report_string += "\t" + rank + "\t"
        if rank not in rank_assigned_dict or len(rank_assigned_dict[rank]) == 0:
            std_out_report_string += "0\t0\t\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n"
        else:
            acc = 0
            for assignments in rank_assigned_dict[rank]:
                for classified in assignments:
                    acc += 1
                    if classified.split("; ")[0] == "Cellular organisms":
                        logging.error("Lineage string cleaning has gone awry somewhere. "
                                      "The root rank should be a Kingdom (e.g. Bacteria or Archaea) but nope.\n")
                        sys.exit(21)
                    optimal, query = assignments[classified]
                    if optimal == classified:
                        offset = 0
                        correct += 1
                    else:
                        offset = determine_offset(classified, optimal)
                        incorrect += 1
                    if offset > 7:
                        # This shouldn't be possible since there are no more than 7 taxonomic ranks
                        logging.error("Offset found to be greater than what is possible (" + str(offset) + ").\n" +
                                      "Classified: " + classified + "\n" +
                                      "Optimal: " + optimal + "\n" +
                                      "Query: " + query + "\n")
                    taxonomic_distance[offset] += 1
            std_out_report_string += str(n_queries) + "\t" + str(n_classified) + "\t\t"

            dist_sum = 0
            for dist in taxonomic_distance:
                dist_sum += taxonomic_distance[dist]
                if taxonomic_distance[dist] > 0:
                    if n_classified == 0:
                        logging.error("No sequences were classified at rank '" + rank +
                                      "' but optimal placements were pointed here. " +
                                      "This is a bug - please alert the developers!\n")
                        sys.exit(21)
                    else:
                        taxonomic_distance[dist] = round(float((taxonomic_distance[dist]*100)/n_classified), 1)
                else:
                    taxonomic_distance[dist] = 0.0
                clade_exclusion_tabular_string += marker_eval_instance.target_marker + "\t" + rank + "\t"
                clade_exclusion_tabular_string += str(n_queries) + "\t" + str(n_classified) + "\t"
                clade_exclusion_tabular_string += str(dist) + "\t" + str(taxonomic_distance[dist])
                clade_exclusion_strings.append(clade_exclusion_tabular_string)
                clade_exclusion_tabular_string = ""
            if dist_sum != n_classified:
                logging.error("Discrepancy between classified sequences at each distance (" + str(dist_sum) +
                              ") and total (" + str(n_classified) + ").\n")
                sys.exit(15)

            std_out_report_string += '\t'.join([str(val) for val in taxonomic_distance.values()]) + "\n"
            if sum(taxonomic_distance.values()) > 101.0:
                logging.error("Sum of proportional assignments at all distances is greater than 100.\n" +
                              "\n".join(["Rank = " + rank,
                                         "Queries = " + str(n_queries),
                                         "Classified = " + str(n_classified),
                                         "Classifications = " + str(len(rank_assigned_dict[rank]))]) + "\n")
                sys.exit(21)

    sys.stdout.write(std_out_report_string)

    return clade_exclusion_strings


def determine_containment(marker_eval_inst: MarkerTest):
    """
    Determines the accuracy of sequence classifications of all sequences contained at different taxonomic ranks

    :param marker_eval_inst: An instance of the MarkerTest class.
        The following instance variables need to be populated:
        marker: Name of the marker currently being evaluated (e.g., nifHc1, mcrA)
        classifications: The rank-wise classification dictionary with
        tuples of (optimal assignment, actual assignment) as values.
        E.g. {"Phylum": {"Proteobacteria": ("Proteobacteria", "Proteobacteria; Alphaproteobacteria")}}
    """
    rank_assigned_dict = marker_eval_inst.classifications
    marker = marker_eval_inst.target_marker
    # Set up collection for this analysis
    lowest_rank = ""
    for depth in sorted(_RANK_DEPTH_MAP):
        if _RANK_DEPTH_MAP[depth] in marker_eval_inst.ranks:
            if marker_eval_inst.get_sensitivity(_RANK_DEPTH_MAP[depth])[1] > 0:
                lowest_rank = _RANK_DEPTH_MAP[depth]
    containment_strings = list()
    n_queries, n_classified, _ = marker_eval_inst.get_sensitivity(lowest_rank)

    sys.stdout.write("Accuracy of " + str(n_classified) + ' ' + marker + " classifications at " + lowest_rank + ":\n")
    sys.stdout.write("\tRank\tCorrect\tToo Shallow (%)\n")

    # Begin parsing through the depths
    if lowest_rank == "Cellular organisms":
        return
    for depth in sorted(_RANK_DEPTH_MAP):
        incorrect_assignments = dict()
        rank = _RANK_DEPTH_MAP[depth]
        if rank in ["Cellular organisms", "Species", "Strain"]:
            continue
        elif rank == lowest_rank:
            break
        correct = 0
        incorrect = 0
        too_shallow = 0
        parse_depth = depth-1
        for assignments in rank_assigned_dict[lowest_rank]:
            for classified in assignments:
                status = 0  # 0 == incorrect, 1 == correct
                optimal, query = assignments[classified]
                classified_lineage = classified.split("; ")
                query_lineage = query.split("; ")
                # print("\nReference:\t", classified)
                # print("Query:\t\t", query)
                # print("Optimal:\t", optimal)
                if len(classified_lineage) > parse_depth:
                    try:
                        if classified_lineage[parse_depth] == query_lineage[parse_depth]:
                            status += 1
                        else:
                            if query not in incorrect_assignments.keys():
                                incorrect_assignments[query] = 0
                            incorrect_assignments[query] += 1
                            incorrect += 1
                    except IndexError:
                        too_shallow += 1
                        # This indicates TreeSAPP placed the sequence too deep in the tree according to its NCBI lineage
                elif len(optimal.split("; ")) > parse_depth:
                    if query not in incorrect_assignments.keys():
                        incorrect_assignments[query] = 0
                    incorrect_assignments[query] += 1
                    incorrect += 1
                else:
                    too_shallow += 1
                correct += status

        percentage_too_shallow = round(float((too_shallow * 100) / n_classified), 1)
        if percentage_too_shallow == 100.0:
            containment_strings.append("\t" + rank + "\t" + str(0.0) + "\t" + str(100.0))
        else:
            containment_strings.append("\t" + rank + "\t" +
                                       str(round((correct * 100) / (correct + incorrect), 1)) + "\t" +
                                       str(percentage_too_shallow))

    sys.stdout.write("\n".join(containment_strings) + "\n")

    return containment_strings


def clean_classification_names(assignments):
    """
    Removes certain ranks, taxonomic super-groups, etc. so strings are more comparable

    :param assignments:
    :return:
    """
    cleaned_assignments = dict()
    for marker in assignments.keys():
        cleaned_assignments[marker] = dict()
        for ref, queries in assignments[marker].items():
            ref = clean_lineage_string(ref)
            cleaned_assignments[marker][ref] = list()
            for query in queries:
                try:
                    query = re.sub('_', "; ", query)
                except TypeError:
                    logging.error("Unexpected type of lineage string (" + str(type(query)) + "):\n" +
                                  str(query) + "\n")
                    sys.exit(21)
                query = clean_lineage_string(query)
                cleaned_assignments[marker][ref].append(query)
    return cleaned_assignments


def map_full_headers(fasta_headers, header_map, assignments, molecule_type):
    """
    Since the headers used throughout the TreeSAPP pipeline are truncated,
    we read the FASTA file and use those headers instead of their short version
    in case valuable information was discarded
    :param fasta_headers: full-length headers that will replace those in assignments
    :param header_map:
    :param assignments: A dictionary of reference (lineage) and query names
    :param molecule_type: prot, nucl, or rrna? Parsed from command-line arguments
    :return: assignments with a full lineage for both reference (keys) and queries (values)
    """
    genome_position_re = re.compile("^([A-Z0-9._]+):.[0-9]+-[0-9]+[_]+.*")
    marker_assignments = dict()
    entrez_query_list = list()
    num_queries = 0
    for robust_classification in assignments.keys():
        marker_assignments[robust_classification] = list()
        for query in assignments[robust_classification]:
            try:
                original_header = header_map['>' + query]
            except KeyError:
                logging.error("Unable to find " + query +
                              " in header_map (constructed from either the input FASTA or .uc file).\n" +
                              "This is probably an error stemming from `reformat_string()`.\n")
                sys.exit(21)
            # print(query, original_header)
            database = molecule_type
            q_accession = ""
            if not re.search(r"Bacteria|Archaea", query):
                # Check for genomic coordinates, indicating these genes will be in the 'nucleotide' database
                if genome_position_re.match(query):
                    q_accession = genome_position_re.match(query).group(1)
                    database = "dna"
                else:
                    header_format_re, header_db, header_molecule = get_header_format(original_header)
                    sequence_info = header_format_re.match(original_header)
                    q_accession, _, _, _, _ = return_sequence_info_groups(sequence_info, header_db, original_header)
                # print(q_accession)
                entrez_query_list.append((q_accession, database))
                marker_assignments[robust_classification].append(q_accession)
                num_queries += 1
            else:
                n_match = 0
                for header in fasta_headers:
                    f_accession = header[1:].split('_')[0]
                    # Useful for headers containing the full lineage
                    if q_accession == f_accession:
                        marker_assignments[robust_classification].append('_'.join(header.split('_')[1:]))
                        n_match += 1
                    else:
                        pass
                if n_match == 0:
                    logging.error("Unable to find matching header for " + query + " in fasta!\n")
                    sys.exit(21)
                elif n_match > 1:
                    logging.error("Headers with identical accessions were identified in fasta!\n" +
                                  "Offending accession: " + q_accession + "\n")
                sys.exit(21)

    return marker_assignments, entrez_query_list


def assign_lineages(complete_ref_seqs, assignments):
    """

    :param assignments:
    :param complete_ref_seqs:
    :return: lineages_list = full_assignments[marker][ref_classification]
    """
    for marker in assignments:
        for robust_classification in assignments[marker]:
            tmp_lineage_list = []
            for accession in assignments[marker][robust_classification]:
                for ref_seq in complete_ref_seqs:
                    if ref_seq.accession == accession:
                        tmp_lineage_list.append(ref_seq.lineage)
            assignments[marker][robust_classification] = tmp_lineage_list
    return assignments


def map_headers_to_lineage(assignments, ref_sequences):
    """
    The alternative function to map_full_headers. Using ReferenceSequence objects,
    :param assignments:
    :param ref_sequences:
    :return:
    """
    lineage_assignments = dict()
    for marker in assignments:
        lineage_assignments[marker] = dict()
        for assigned_lineage in assignments[marker].keys():
            classified_headers = assignments[marker][assigned_lineage]
            c_lineage = clean_lineage_string(assigned_lineage)
            lineage_assignments[marker][c_lineage] = list()
            for query in classified_headers:
                for ref_seq in ref_sequences:
                    if ref_seq.accession == query:
                        lineage_assignments[marker][c_lineage].append(clean_lineage_string(ref_seq.lineage))
                        break
            if len(lineage_assignments[marker][c_lineage]) > len(classified_headers):
                logging.error(str(len(classified_headers)) + " accessions mapped to " +
                              str(len(lineage_assignments[marker][c_lineage])) + " lineages.\n")
                sys.exit(21)
            elif len(lineage_assignments[marker][c_lineage]) < len(classified_headers):
                logging.debug(str(len(classified_headers)) + " accessions mapped to " +
                              str(len(lineage_assignments[marker][c_lineage])) + " lineages.\n")
    return lineage_assignments


def get_unclassified_rank(pos, split_lineage):
    """
    Recursive function to retrieve the first rank at which the
    :param pos:
    :param split_lineage:
    :return:
    """
    if re.search("nclassified", split_lineage[pos]):
        return pos
    else:
        pos = get_unclassified_rank(pos+1, split_lineage)
    return pos


def pick_taxonomic_representatives(ref_seqs_list, taxonomic_filter_stats, max_cluster_size=5):
    """
    Removes queries with duplicate taxa - to prevent the taxonomic composition of the input
    from biasing the accuracy to over- or under-perform by classifying many sequences from very few groups.
    Also removes taxonomies with "*nclassified" in their lineage or are derived from environmental samples

    :param ref_seqs_list: A dictionary mapping accessions to lineages that need to be filtered
    :param taxonomic_filter_stats: A dictionary for tracking the number sequences filtered, those unique, etc.
    :param max_cluster_size: The maximum number of sequences representing a taxonomic cluster
    :return: dereplicated_lineages dict with lineages mapping to a (short) list of accessions
    """
    good_classified_lineages = dict()
    dereplicated_lineages = dict()
    num_rep_seqs = 0
    for ref_seq in ref_seqs_list:
        query_taxonomy = clean_lineage_string(ref_seq.lineage)
        if len(clean_lineage_string(query_taxonomy).split("; ")) < 7:
            continue
        if query_taxonomy not in good_classified_lineages:
            good_classified_lineages[query_taxonomy] = list()
        if re.search("unclassified|environmental sample", query_taxonomy, re.IGNORECASE):
            # # Remove taxonomic lineages that are unclassified at the Phylum level or higher
            # unclassified_depth = get_unclassified_rank(0, query_taxonomy.split("; "))
            # if unclassified_depth > 4:
            #     good_classified_lineages[query_taxonomy].append(ref_seq.accession)
            # else:
            taxonomic_filter_stats["Unclassified"] += 1
        else:
            taxonomic_filter_stats["Classified"] += 1
            good_classified_lineages[query_taxonomy].append(ref_seq.accession)

    if taxonomic_filter_stats["Unclassified"] == len(ref_seqs_list):
        logging.error("All sequences provided are derived from uncultured, unclassified organisms.\n")
        sys.exit(21)

    # In order to maintain consistency among multiple runs with the same input
    # a separate loop is required for sorting
    taxonomic_filter_stats["Max"] = 0
    taxonomic_filter_stats["Min"] = 1000
    taxonomic_filter_stats["Mean"] = 0
    for query_taxonomy in good_classified_lineages:
        dereplicated_lineages[query_taxonomy] = list()
        i = 0
        lineage_candidates = sorted(good_classified_lineages[query_taxonomy])
        while i < len(lineage_candidates) and i < max_cluster_size:
            dereplicated_lineages[query_taxonomy].append(lineage_candidates[i])
            i += 1

        # Generate stats on the taxonomic clusters
        cluster_size = len(dereplicated_lineages[query_taxonomy])
        num_rep_seqs += cluster_size
        if cluster_size == 0:
            # This is a taxonomic lineage removed by the "Unclassified" filter
            continue
        if cluster_size > taxonomic_filter_stats["Max"]:
            taxonomic_filter_stats["Max"] = cluster_size
        if cluster_size < taxonomic_filter_stats["Min"]:
            taxonomic_filter_stats["Min"] = cluster_size
        taxonomic_filter_stats["Mean"] = round(float(num_rep_seqs / len(dereplicated_lineages)), 2)

    taxonomic_filter_stats["Unique_taxa"] += len(dereplicated_lineages)

    logging.info("\t" + str(num_rep_seqs) + " representative sequences will be used for TreeSAPP analysis.\n")

    logging.debug("Representative sequence stats:\n\t" +
                  "Maximum representative sequences for a taxon " + str(taxonomic_filter_stats["Max"]) + "\n\t" +
                  "Minimum representative sequences for a taxon " + str(taxonomic_filter_stats["Min"]) + "\n\t" +
                  "Mean representative sequences for a taxon " + str(taxonomic_filter_stats["Mean"]) + "\n")

    return dereplicated_lineages, taxonomic_filter_stats


def select_rep_seqs(deduplicated_assignments: dict, test_sequences: list, taxon=None):
    """
    Function for creating a fasta-formatted dict from the accessions representing unique taxa in the test sequences

    :param deduplicated_assignments: dict of lineages mapped to a list of accessions
    :param test_sequences: list of ReferenceSequence objects
    :param taxon: A taxonomic lineage to filter the lineages (and sequences) in deduplicated assignments
    :return: Dictionary containing accessions as keys and sequences as values
    """
    if taxon:
        filtered_assignments = dict()
        for lineage in sorted(deduplicated_assignments):
            if re.search(taxon, lineage):
                filtered_assignments[lineage] = deduplicated_assignments[lineage]
    else:
        filtered_assignments = deduplicated_assignments

    deduplicated_fasta_dict = dict()
    for lineage in sorted(filtered_assignments):
        for accession in filtered_assignments[lineage]:
            matched = False
            for ref_seq in test_sequences:
                if ref_seq.accession == accession:
                    deduplicated_fasta_dict[accession] = ref_seq.sequence
                    matched = True
            if not matched:
                logging.error("Unable to find accession (" + accession + ") in accession-lineage map\n")
                sys.exit(21)
    return deduplicated_fasta_dict


def map_seqs_to_lineages(accession_lineage_map, deduplicated_fasta_dict):
    seq_taxa_map = dict()
    for seq_name in deduplicated_fasta_dict:
        if seq_name not in accession_lineage_map:
            logging.error("Unable to find matching key for '" + seq_name + "' in accession_lineage_map.\n")
            sys.exit(21)
        else:
            seq_taxa_map[seq_name] = accession_lineage_map[seq_name]
    return seq_taxa_map


def filter_queries_by_taxonomy(taxonomic_lineages):
    """
    Removes queries with duplicate taxa - to prevent the taxonomic composition of the input
    from biasing the accuracy to over- or under-perform by classifying many sequences from very few groups.
    Also removes taxonomies with "*nclassified" in their lineage

    :param taxonomic_lineages: A list of lineages that need to be filtered
    :return: A list that contains no more than 3 of each query taxonomy (arbitrarily normalized) and counts
    """
    unclassifieds = 0
    classified = 0
    unique_query_taxonomies = 0
    lineage_enumerator = dict()
    normalized_lineages = list()
    for query_taxonomy in sorted(taxonomic_lineages):
        can_classify = True
        if re.search("unclassified", query_taxonomy, re.IGNORECASE):
            # Remove taxonomic lineages that are unclassified at the Phylum level or higher
            unclassified_depth = get_unclassified_rank(0, query_taxonomy.split("; "))
            if unclassified_depth <= 3:
                can_classify = False

        if can_classify:
            if query_taxonomy not in lineage_enumerator:
                lineage_enumerator[query_taxonomy] = 0
                classified += 1
            if lineage_enumerator[query_taxonomy] < 3:
                normalized_lineages.append(query_taxonomy)
            lineage_enumerator[query_taxonomy] += 1
        else:
            unclassifieds += 1
    unique_query_taxonomies += len(lineage_enumerator)

    return normalized_lineages, unclassifieds, classified, unique_query_taxonomies


def write_containment_table(args, containment_table, containment_strings):
    try:
        output_handler = open(containment_table, 'w')
    except IOError:
        logging.error("Unable to open " + containment_table + " for writing.\n")
        sys.exit(21)

    output_handler.write("# Input file for testing: " + args.fasta_input + "\n")
    output_name = os.path.dirname(args.output)
    for line in containment_strings:
        # Line has a "\t" prefix already
        line = output_name + "\t" + args.reference_marker + "\t" + args.tool + line + "\n"
        output_handler.write(line)

    output_handler.close()
    return


def write_performance_table(args, performance_table, clade_exclusion_strings):
    try:
        output_handler = open(performance_table, 'w')
    except IOError:
        logging.error("Unable to open " + performance_table + " for writing.\n")
        sys.exit(21)

    output_handler.write("# Input file for testing: " + args.fasta_input + "\n")
    output_name = os.path.dirname(args.output)
    for line in clade_exclusion_strings:
        line = output_name + "\t" + args.tool + "\t" + line + "\n"
        output_handler.write(line)

    output_handler.close()
    return


def correct_accession(description):
    header_format_re, header_db, header_molecule = get_header_format(description)
    sequence_info = header_format_re.match(description)
    accession, _, _, _, _ = return_sequence_info_groups(sequence_info, header_db, description)
    return accession


def load_ref_seqs(fasta_dict, header_registry, ref_seq_dict):
    """
    Function for adding sequences from a fasta-formatted dictionary into dictionary of ReferenceSequence objects

    :param fasta_dict:
    :param header_registry: An optional dictionary of Header objects
    :param ref_seq_dict: A dictionary indexed by arbitrary integers mapping to ReferenceSequence instances
    :return:
    """
    missing = list()
    if len(header_registry) != len(fasta_dict):
        logging.warning("Number of records in FASTA collection and header list differ.\n" +
                        "Chances are these were short sequences that didn't pass the filter. Carrying on.\n")

    for num_id in ref_seq_dict.keys():
        ref_seq = ref_seq_dict[num_id]
        formatted_header = header_registry[num_id].formatted
        try:
            ref_seq.sequence = fasta_dict[formatted_header]
        except KeyError:
            if len(header_registry) == len(fasta_dict):
                logging.error(formatted_header + " not found in FASTA records due to format incompatibilities.\n")
                sys.exit(21)
            missing.append(str(header_registry[num_id].original))
    if len(missing) > 0:
        logging.debug("The following sequences have been removed from further analyses:\n\t" +
                      "\n\t".join(missing) + "\n")
    return ref_seq_dict


def prep_graftm_ref_files(treesapp_dir, intermediate_dir, target_clade, marker, depth):
    # Move the original FASTA, tree and tax_ids files to a temporary location
    marker_fa = os.sep.join([treesapp_dir, "data", "alignment_data", marker + ".fa"])
    marker_tax_ids = os.sep.join([treesapp_dir, "data", "tree_data", "tax_ids_" + marker + ".txt"])
    off_target_ref_leaves = dict()
    # tax_ids file
    ref_tree_leaves = tax_ids_file_to_leaves(marker_tax_ids)
    with open(intermediate_dir + "tax_ids_" + marker + ".txt", 'w') as tax_ids_handle:
        tax_ids_strings = list()
        for ref_leaf in ref_tree_leaves:
            c_lineage = clean_lineage_string(ref_leaf.lineage)
            if re.search(target_clade, c_lineage):
                continue
            sc_lineage = c_lineage.split("; ")
            if len(sc_lineage) < depth:
                continue
            if re.search("unclassified|environmental sample", c_lineage, re.IGNORECASE):
                i = 0
                while i <= depth:
                    if re.search("unclassified|environmental sample", sc_lineage[i], re.IGNORECASE):
                        i -= 1
                        break
                    i += 1
                if i < depth:
                    continue
            organism, accession = ref_leaf.description.split(" | ")
            off_target_ref_leaves[ref_leaf.number] = accession
            tax_ids_strings.append(accession + "\t" + clean_lineage_string(ref_leaf.lineage))
        tax_ids_handle.write("\n".join(tax_ids_strings) + "\n")

    # fasta
    ref_fasta_dict = read_fasta_to_dict(marker_fa)
    accession_fasta_dict = dict()
    for key_id in ref_fasta_dict:
        num_key = key_id.split('_')[0]
        if num_key in off_target_ref_leaves.keys():
            accession_fasta_dict[off_target_ref_leaves[num_key]] = ref_fasta_dict[key_id]
        else:
            pass
    write_new_fasta(accession_fasta_dict, intermediate_dir + marker + ".mfa")
    for acc in accession_fasta_dict:
        accession_fasta_dict[acc] = re.sub('-', '', accession_fasta_dict[acc])
    write_new_fasta(accession_fasta_dict, intermediate_dir + marker + ".fa")
    return


def exclude_clade_from_ref_files(treesapp_dir, marker, intermediate_dir, target_clade, depth,
                                 executables, fresh, molecule):
    # Move the original FASTA, tree and tax_ids files to a temporary location
    marker_fa = os.sep.join([treesapp_dir, "data", "alignment_data", marker + ".fa"])
    marker_hmm = os.sep.join([treesapp_dir, "data", "hmm_data", marker + ".hmm"])
    marker_tree = os.sep.join([treesapp_dir, "data", "tree_data", marker + "_tree.txt"])
    marker_bipart_tree = os.sep.join([treesapp_dir, "data", "tree_data", marker + "_bipartitions.txt"])
    marker_tax_ids = os.sep.join([treesapp_dir, "data", "tree_data", "tax_ids_" + marker + ".txt"])
    intermediate_prefix = intermediate_dir + "ORIGINAL"

    shutil.copy(marker_fa, intermediate_prefix + ".fa")
    shutil.copy(marker_hmm, intermediate_prefix + ".hmm")
    shutil.copy(marker_tree, intermediate_prefix + "_tree.txt")
    if os.path.isfile(marker_bipart_tree):
        shutil.copy(marker_bipart_tree, intermediate_prefix + "_bipartitions.txt")
        os.remove(marker_bipart_tree)
    shutil.copy(marker_tax_ids, intermediate_prefix + "_tax_ids.txt")

    off_target_ref_leaves = list()
    n_match = 0
    n_shallow = 0
    n_unclassified = 0
    # tax_ids file
    ref_tree_leaves = tax_ids_file_to_leaves(marker_tax_ids)
    with open(marker_tax_ids, 'w') as tax_ids_handle:
        for ref_leaf in ref_tree_leaves:
            tax_ids_string = ""
            c_lineage = clean_lineage_string(ref_leaf.lineage)
            if re.search(target_clade, c_lineage):
                n_match += 1
                continue
            sc_lineage = c_lineage.split("; ")
            if len(sc_lineage) < depth:
                n_shallow += 1
                continue
            if re.search("unclassified|environmental sample", c_lineage, re.IGNORECASE):
                i = 0
                while i <= depth:
                    if re.search("unclassified|environmental sample", sc_lineage[i], re.IGNORECASE):
                        i -= 1
                        break
                    i += 1
                if i < depth:
                    n_unclassified += 1
                    continue
            off_target_ref_leaves.append(ref_leaf.number)
            tax_ids_string += "\t".join([ref_leaf.number, ref_leaf.description, ref_leaf.lineage])
            tax_ids_handle.write(tax_ids_string + "\n")

    # print("Matching =", n_match, "\tUnclassified =", n_unclassified, "\tShallow =", n_shallow)

    # fasta
    ref_fasta_dict = read_fasta_to_dict(marker_fa)
    off_target_headers = [num_id + '_' + marker for num_id in off_target_ref_leaves]
    split_files = write_new_fasta(ref_fasta_dict, marker_fa, len(off_target_ref_leaves)+1, off_target_headers)
    if len(split_files) > 1:
        logging.error("Only one FASTA file should have been written.\n")
        sys.exit(21)
    else:
        shutil.copy(split_files[0], marker_fa)

    # HMM profile
    hmm_build_command = [executables["hmmbuild"],
                         marker_hmm,
                         marker_fa]
    launch_write_command(hmm_build_command)

    # Trees
    if fresh:
        tree_build_cmd = [executables["FastTree"]]
        if molecule == "rrna" or molecule == "dna":
            tree_build_cmd += ["-nt", "-gtr"]
        else:
            tree_build_cmd += ["-lg", "-wag"]
        tree_build_cmd += ["-out", marker_tree]
        tree_build_cmd.append(marker_fa)
        logging.info("Building Approximately-Maximum-Likelihood tree with FastTree... ")
        stdout, returncode = launch_write_command(tree_build_cmd, True)
        with open(intermediate_dir + os.sep + "FastTree_info." + marker, 'w') as fast_info:
            fast_info.write(stdout + "\n")
        logging.info("done.\n")
    else:
        ref_tree = Tree(marker_tree)
        ref_tree.prune(off_target_ref_leaves)
        logging.debug("\t" + str(len(ref_tree.get_leaves())) + " leaves in pruned tree.\n")
        ref_tree.write(outfile=marker_tree, format=5)

    return intermediate_prefix


def validate_ref_package_files(treesapp_dir, marker, intermediate_dir):
    # Move the original FASTA, tree and tax_ids files to a temporary location
    marker_fa = os.sep.join([treesapp_dir, "data", "alignment_data", marker + ".fa"])
    marker_tree = os.sep.join([treesapp_dir, "data", "tree_data", marker + "_tree.txt"])
    marker_bipart_tree = os.sep.join([treesapp_dir, "data", "tree_data", marker + "_bipartitions.txt"])
    marker_tax_ids = os.sep.join([treesapp_dir, "data", "tree_data", "tax_ids_" + marker + ".txt"])
    intermediate_prefix = intermediate_dir + "ORIGINAL"

    # This prevents users from quitting mid-analysis then using bad reference package files with clades removed
    remnants = glob(intermediate_prefix + "*")
    if len(remnants) > 0:
        logging.warning("Intermediate files from incomplete analysis exist:\n" + "\n".join(remnants) + "\n")
        logging.info("Attempting recovery... ")
        # Easiest way to recover is to move the originals back to proper location and proceed
        try:
            shutil.copy(intermediate_prefix + ".fa", marker_fa)
            shutil.copy(intermediate_prefix + "_tree.txt", marker_tree)
            if os.path.isfile(intermediate_prefix + "_bipartitions.txt"):
                shutil.copy(intermediate_prefix + "_bipartitions.txt", marker_bipart_tree)
            shutil.copy(intermediate_prefix + "_tax_ids.txt", marker_tax_ids)
            logging.info("succeeded.\n")
        except FileNotFoundError:
            logging.info("failed. Redownload data files and start over.\n")
            sys.exit(21)
    return


def remove_clade_exclusion_files(intermediate_dir):
    remnants = glob(intermediate_dir + "ORIGINAL*")
    for tmp_file in remnants:
        os.remove(tmp_file)


def classify_excluded_taxon(args, prefix, output_dir, marker, min_seq_length, test_rep_taxa_fasta):
    """
      Prepares TreeSAPP tree, alignment and taxonomic identification map (tax_ids) files for clade exclusion analysis,
    and performs classification with TreeSAPP

    :return: The paths to the classification table and the taxon-excluded tax_ids file
    """

    # Classify representative sequences using TreeSAPP
    classify_command = [args.treesapp + "/treesapp.py", "-i", test_rep_taxa_fasta,
                        "-o", output_dir,
                        "-m", args.molecule,
                        "-T", str(args.threads),
                        "--min_seq_length", min_seq_length,
                        "--trim_align",
                        "--overwrite",
                        "--delete"]
    logging.debug("Command used:\n" + ' '.join(classify_command) + "\n")
    launch_write_command(classify_command, False)
    # Move the original FASTA, tree and tax_ids files back to the proper directories
    shutil.copy(prefix + "_tree.txt", os.sep.join([args.treesapp, "data", "tree_data", marker + "_tree.txt"]))
    if os.path.isfile(prefix + "_bipartitions.txt"):
        shutil.copy(prefix + "_bipartitions.txt",
                    os.sep.join([args.treesapp, "data", "tree_data", marker + "_bipartitions.txt"]))
    # The edited tax_ids file with clade excluded is required for performance analysis
    shutil.copy(os.sep.join([args.treesapp, "data", "tree_data", "tax_ids_" + marker + ".txt"]), output_dir)
    shutil.copy(prefix + "_tax_ids.txt", os.sep.join([args.treesapp, "data", "tree_data", "tax_ids_" + marker + ".txt"]))
    shutil.copy(prefix + ".fa", os.sep.join([args.treesapp, "data", "alignment_data", marker + ".fa"]))
    shutil.copy(prefix + ".hmm", os.sep.join([args.treesapp, "data", "hmm_data", marker + ".hmm"]))

    return


def build_graftm_package(target_marker, output_dir, mfa_file, fa_file, threads):
    create_command = ["graftM", "create"]
    create_command += ["--threads", str(threads)]
    create_command += ["--alignment", mfa_file]
    create_command += ["--sequences", fa_file]
    create_command += ["--taxonomy", output_dir + os.sep + "tax_ids_" + target_marker + ".txt"]
    create_command += ["--output", output_dir + target_marker + ".gpkg"]
    create_command.append("--force")

    logging.debug("Command used:\n" + ' '.join(create_command) + "\n")
    launch_write_command(create_command, False)


def graftm_classify(test_rep_taxa_fasta, gpkg_path, output_dir, threads, tool):
    classify_command = ["graftM", "graft"]
    classify_command += ["--forward", test_rep_taxa_fasta]
    classify_command += ["--graftm_package", gpkg_path]
    classify_command += ["--threads", str(threads)]
    if tool == "graftm":
        classify_command += ["--assignment_method", "pplacer"]
        classify_command += ["--search_method", "hmmsearch"]
    elif tool == "diamond":
        classify_command += ["--assignment_method", "diamond"]
        classify_command += ["--search_method", "diamond"]
    classify_command += ["--output_directory", output_dir]
    classify_command += ["--input_sequence_type", "aminoacid"]
    classify_command.append("--force")

    logging.debug("Command used:\n" + ' '.join(classify_command) + "\n")
    launch_write_command(classify_command, False)

    return


def main():
    """
    Method for running this script:
        Provide it a FASTA file for which it will determine the taxonomic lineage for each sequence
         and run all taxonomic representative sequences with TreeSAPP then analyze via clade exclusion

    :return:
    """
    sys.stdout.write("\n##\t\t\tBeginning clade exclusion analysis\t\t\t##\n")

    args = get_arguments()
    args.targets = ["ALL"]
    args = find_executables(args)

    os.makedirs(args.output, exist_ok=True)
    log_file_name = args.output + os.sep + "Clade_exclusion_analyzer_log.txt"
    prep_logging(log_file_name, args.verbose)
    logging.debug("Command used:\n" + ' '.join(sys.argv) + "\n")

    marker_build_dict = parse_ref_build_params(args)
    ref_leaves = tax_ids_file_to_leaves(os.sep.join([args.treesapp, 'data',  'tree_data',
                                                     "tax_ids_" + args.reference_marker + ".txt"]))
    ref_lineages = dict()
    for leaf in ref_leaves:
        ref_lineages[leaf.number] = leaf.lineage
    marker_eval_inst = MarkerTest(args.reference_marker)
    taxa_choices = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]
    args.taxon_rank = args.taxon_rank.split(',')
    for rank in args.taxon_rank:
        if rank not in taxa_choices:
            logging.error(rank + " not an available option for `--taxon_rank`.\n")
            sys.exit(21)
        else:
            marker_eval_inst.ranks.append(rank)

    # Alert the user if the denominator format was (incorrectly) provided to Clade_exclusion_analyzer
    if re.match("^[A-Z][0-9]{4}$", args.reference_marker):
        code_name = args.reference_marker
        try:
            marker = marker_build_dict[code_name].cog
        except KeyError:
            logging.error("Unable to find '" + code_name + "' in ref_build_parameters collection!\n" +
                          "Has it been added to data/tree_data/ref_build_parameters.tsv and cog_list.tsv?\n")
            sys.exit(21)
    elif len(args.reference_marker) <= 6:
        # Assuming the provided args.reference_marker is a short gene name (e.g. mcrA, nifH)
        marker = args.reference_marker
        code_name = ""
        for denominator in marker_build_dict:
            if marker_build_dict[denominator].cog == marker:
                code_name = denominator
                break
        if not code_name:
            logging.error("Unable to identify the gene name from the code name '" + args.reference_marker + "'.")
            sys.exit(21)
    else:
        logging.error("Wrong format for the reference code_name provided: " + args.reference_marker + "\n")
        sys.exit(21)

    ##
    # Define locations of files TreeSAPP outputs
    ##
    inter_class_file = args.output + os.sep + "tmp_clade_exclusion_assignments.tsv"
    accession_map_file = args.output + os.sep + "accession_id_lineage_map.tsv"
    test_rep_taxa_fasta = args.output + os.sep + "representative_taxa_sequences.fasta"
    performance_table = args.output + os.sep + "clade_exclusion_performance.tsv"
    containment_table = args.output + os.sep + "accuracy.tsv"
    # Determine the analysis stage and user's intentions with four booleans
    accessions_downloaded = False  # Whether the accessions have been downloaded from Entrez
    classified = False  # Has TreeSAPP been completed for these sequences?
    treesapp_output_dir = args.output + "TreeSAPP_output" + os.sep
    # Working with a pre-existing TreeSAPP output directory
    if os.path.exists(treesapp_output_dir):
        extant = True
    else:
        extant = False
    if extant:
        if os.path.isfile(accession_map_file):
            accessions_downloaded = True

    if extant:
        logging.debug("The output directory has been detected.\n")

    if not os.path.isdir(treesapp_output_dir):
        os.makedirs(treesapp_output_dir)

    if args.tool in ["diamond", "graftm"]:
        graftm_files = glob(treesapp_output_dir + os.sep + "*" + os.sep + "*_read_tax.tsv")
        if len(graftm_files) == 1:
            classification_table = glob(treesapp_output_dir + os.sep + "*" + os.sep + "*_read_tax.tsv")[0]
        else:
            classification_table = ''
    else:
        classification_table = treesapp_output_dir + os.sep + "final_outputs" + os.sep + "marker_contig_map.tsv"
    if os.path.isfile(classification_table):
        classified = True

    # Load FASTA data
    fasta_dict = format_read_fasta(args.fasta_input, args.molecule, args.output, 110)
    if args.length:
        for seq_id in fasta_dict:
            if len(fasta_dict[seq_id]) < args.length:
                logging.warning(seq_id + " sequence is shorter than " + str(args.length) + "\n")
            else:
                max_stop = len(fasta_dict[seq_id]) - args.length
                random_start = randint(0, max_stop)
                fasta_dict[seq_id] = fasta_dict[seq_id][random_start:random_start+args.length]
    header_registry = register_headers(get_headers(args.fasta_input))
    # Load the query test sequences as ReferenceSequence objects
    complete_ref_sequences = get_header_info(header_registry)
    complete_ref_sequences = load_ref_seqs(fasta_dict, header_registry, complete_ref_sequences)

    logging.debug("\tNumber of input sequences =\t" + str(len(complete_ref_sequences)) + "\n")

    # Checkpoint one: does anything exist?
    # If not, begin by downloading lineage information for each sequence accession
    args.output = treesapp_output_dir
    if not extant and not accessions_downloaded and not classified:
        # User hasn't analyzed anything, sequences will be taxonomically screened then analyzed
        # NOTE: This is only supported for analyzing a single marker gene
        extant = True
        entrez_query_list, num_lineages_provided = build_entrez_queries(complete_ref_sequences)
        if len(entrez_query_list) > 0:
            entrez_records = get_multiple_lineages(entrez_query_list, args.molecule)
            accession_lineage_map = entrez_records_to_accession_lineage_map(entrez_records)
            all_accessions = entrez_records_to_accessions(entrez_records, entrez_query_list)

        elif num_lineages_provided > 0:
            logging.info("No Entrez queries are necessary - all sequences have provided lineage information.\n")
            accession_lineage_map = dict()
            all_accessions = []
        else:
            logging.error("No accessions were parsed from FASTA records in " + args.fasta_input)
            sys.exit(21)
        # Download lineages separately for those accessions that failed,
        # map proper accession to lineage from the tuple keys (accession, accession.version)
        #  in accession_lineage_map returned by get_multiple_lineages.
        complete_ref_sequences, accession_lineage_map = verify_lineage_information(accession_lineage_map,
                                                                                   all_accessions,
                                                                                   complete_ref_sequences,
                                                                                   num_lineages_provided)
        write_accession_lineage_map(accession_map_file, accession_lineage_map)
        complete_ref_sequences = finalize_ref_seq_lineages(complete_ref_sequences, accession_lineage_map)
        accessions_downloaded = True
    elif extant and accessions_downloaded:
        # File being read should contain accessions mapped to their lineages for all sequences in input FASTA
        accession_lineage_map = read_accession_taxa_map(accession_map_file)
        complete_ref_sequences = finalize_ref_seq_lineages(complete_ref_sequences, accession_lineage_map)
    else:
        logging.error("Logic error.")
        sys.exit(21)
    fasta_record_objects = complete_ref_sequences.values()

    logging.info("Selecting representative sequences for each taxon.\n")

    # Filter the sequences from redundant taxonomic lineages, picking up to 5 representative sequences
    representative_seqs, marker_eval_inst.taxa_filter = pick_taxonomic_representatives(fasta_record_objects,
                                                                                       marker_eval_inst.taxa_filter)
    deduplicated_fasta_dict = select_rep_seqs(representative_seqs, fasta_record_objects)
    write_new_fasta(deduplicated_fasta_dict, test_rep_taxa_fasta)
    rep_accession_lineage_map = map_seqs_to_lineages(accession_lineage_map, deduplicated_fasta_dict)

    # Checkpoint three: We have accessions linked to taxa, and sequences to analyze with TreeSAPP, but not classified

    if extant and accessions_downloaded and not classified:
        # Run TreeSAPP against the provided tax_ids file and the unique taxa FASTA file
        if args.length:
            min_seq_length = str(min(args.length - 10, 30))
        else:
            min_seq_length = str(30)

        validate_ref_package_files(args.treesapp, marker, args.output)

        ranks = {"Kingdom": 0, "Phylum": 1, "Class": 2, "Order": 3, "Family": 4, "Genus": 5, "Species": 6}
        for rank in args.taxon_rank:
            leaf_trimmed_taxa_map = trim_lineages_to_rank(ref_lineages, rank)
            unique_ref_lineages = sorted(set(leaf_trimmed_taxa_map.values()))
            unique_query_lineages = sorted(set(trim_lineages_to_rank(rep_accession_lineage_map, rank).values()))
            depth = ranks[rank]
            for lineage in unique_query_lineages:
                taxon = re.sub(r"([ /])", '_', lineage.split("; ")[-1])
                rank_tax = rank[0] + '_' + taxon
                treesapp_output = args.output + os.sep + rank_tax + os.sep

                optimal_lca_taxonomy = "; ".join(lineage.split("; ")[:-1])
                if optimal_lca_taxonomy not in ["; ".join(tl.split("; ")[:-1]) for tl in unique_ref_lineages
                                                if tl != lineage]:
                    logging.debug("Optimal placement target '" + optimal_lca_taxonomy + "' not in pruned tree.\n")
                    continue

                # Select representative sequences belonging to the taxon being tested
                taxon_rep_seqs = select_rep_seqs(representative_seqs, fasta_record_objects, lineage)
                # Decide whether to continue analyzing taxon based on number of query sequences
                if len(taxon_rep_seqs.keys()) == 0:
                    logging.debug("No query sequences for " + lineage + ".\n")
                    continue

                logging.info("Classifications for the " + rank + " '" + taxon + "' put " + treesapp_output + "\n")
                test_obj = marker_eval_inst.new_taxa_test(rank, lineage)
                test_obj.queries = taxon_rep_seqs.keys()
                test_rep_taxa_fasta = args.output + rank_tax + ".fa"

                if args.tool in ["graftm", "diamond"]:
                    classification_table = os.sep.join([treesapp_output, rank_tax, rank_tax + "_read_tax.tsv"])

                    if not os.path.isfile(classification_table):
                        tax_ids_file = os.sep.join([args.output,
                                                    marker + ".gpkg",
                                                    marker + ".gpkg.refpkg",
                                                    marker + "_taxonomy.csv"])
                        # Copy reference files, then exclude all clades belonging to the taxon being tested
                        prep_graftm_ref_files(args.treesapp, args.output, lineage, marker_eval_inst.target_marker, depth)
                        build_graftm_package(marker_eval_inst.target_marker,
                                             args.output,
                                             mfa_file=args.output + marker + ".mfa",
                                             fa_file=args.output + marker + ".fa",
                                             threads=args.threads)
                        # Write the query sequences
                        write_new_fasta(taxon_rep_seqs, test_rep_taxa_fasta)

                        graftm_classify(test_rep_taxa_fasta, args.output + os.sep + marker + ".gpkg",
                                        treesapp_output, args.threads, args.tool)

                        if not os.path.isfile(classification_table):
                            # The TaxonTest object is maintained for record-keeping (to track # queries & classifieds)
                            logging.warning("GraftM did not generate output for " + lineage + ". Skipping.\n")
                            # shutil.rmtree(treesapp_output)
                            continue

                        shutil.copy(tax_ids_file, treesapp_output + os.sep + rank_tax + os.sep)

                    tax_ids_file = os.sep.join([treesapp_output, rank_tax, marker + "_taxonomy.csv"])
                    test_obj.taxonomic_tree = grab_graftm_taxa(tax_ids_file)
                    graftm_assignments = read_graftm_classifications(classification_table)
                    test_obj.assignments = {marker: graftm_assignments}
                    test_obj.filter_assignments(marker_eval_inst.target_marker)
                else:
                    tax_ids_file = treesapp_output + "tax_ids_" + marker + ".txt"
                    classification_table = treesapp_output + "final_outputs" + os.sep + "marker_contig_map.tsv"

                    if not os.path.isfile(classification_table) or not os.path.isfile(tax_ids_file):
                        # Copy reference files, then exclude all clades belonging to the taxon being tested
                        prefix = exclude_clade_from_ref_files(args.treesapp, marker, args.output, lineage, depth,
                                                              args.executables, args.fresh, args.molecule)
                        # Write the query sequences
                        write_new_fasta(taxon_rep_seqs, test_rep_taxa_fasta)
                        classify_excluded_taxon(args, prefix, treesapp_output, marker, min_seq_length, test_rep_taxa_fasta)
                        if not os.path.isfile(classification_table):
                            # The TaxonTest object is maintained for record-keeping (to track # queries & classifieds)
                            logging.warning("TreeSAPP did not generate output for " + lineage + ". Skipping.\n")
                            shutil.rmtree(treesapp_output)
                            continue
                    else:
                        # Valid number of queries and these sequences have already been classified
                        pass

                    test_obj.taxonomic_tree = all_possible_assignments(tax_ids_file)
                    if os.path.isfile(classification_table):
                        assigned_lines = read_marker_classification_table(classification_table)
                        test_obj.assignments = parse_assignments(assigned_lines)
                        test_obj.filter_assignments(marker_eval_inst.target_marker)
                        test_obj.distances = parse_distances(assigned_lines)
                    else:
                        logging.error("marker_contig_map.tsv is missing from output directory '" +
                                      os.path.basename(classification_table) + "'\n" +
                                      "Please remove this directory and re-run.\n")
                        sys.exit(21)
        classified = True
    remove_clade_exclusion_files(args.output)

    # Checkpoint four: everything has been prepared, only need to parse the classifications and map lineages
    logging.info("Finishing up the mapping of classified, filtered taxonomic sequences.\n")
    for rank in sorted(marker_eval_inst.taxa_tests):
        for test_obj in marker_eval_inst.taxa_tests[rank]:
            if test_obj.assignments:
                marker_assignments = map_headers_to_lineage(test_obj.assignments, fasta_record_objects)
                # Return the number of correct, classified, and total sequences of that taxon at the current rank
                # Identify the excluded rank for each query sequence
                if len(marker_assignments) == 0:
                    logging.debug("No sequences were classified for " + test_obj.taxon + "\n")
                    continue

                for marker in marker_assignments:
                    marker_eval_inst.markers.add(marker)

                rank_assignments = identify_excluded_clade(marker_assignments,
                                                           test_obj.taxonomic_tree,
                                                           marker_eval_inst.target_marker)
                for a_rank in rank_assignments:
                    if a_rank != rank and len(rank_assignments[a_rank]) > 0:
                        logging.warning(rank + "-level clade excluded but classifications were found to be " + a_rank +
                                        "-level.\nAssignments were: " + str(rank_assignments[a_rank]) + "\n")
                        continue
                    if a_rank not in marker_eval_inst.classifications:
                        marker_eval_inst.classifications[a_rank] = list()
                    if len(rank_assignments[a_rank]) > 0:
                        marker_eval_inst.classifications[a_rank] += rank_assignments[a_rank]

    # TODO: In the case of a prior TreeSAPP analysis without taxonomic sequence filtering (external of Clade_exclusion)
    # if extant and classified:
    #     logging.info("Outputs from a previous unsupervised TreeSAPP analysis found.\n")
    #     # Read the classification table
    #     logging.error("This functionality is now deprecated.\n" +
    #                   "Software no longer analyzes TreeSAPP outputs created without Clade_exclusion_analyzer.py\n")
    #     sys.exit(21)

    ##
    # On to the standard clade-exclusion analysis...
    ##
    if marker_eval_inst.taxa_filter["Classified"] != marker_eval_inst.taxa_filter["Unique_taxa"]:
        logging.debug("\n\t" + str(marker_eval_inst.taxa_filter["Classified"] -
                                   marker_eval_inst.taxa_filter["Unique_taxa"]) +
                      " duplicate query taxonomies removed.\n")

    if marker_eval_inst.taxa_filter["Unclassified"] > 0:
        logging.debug("\t" + str(marker_eval_inst.taxa_filter["Unclassified"]) +
                      " query sequences with unclassified taxonomies were removed.\n" +
                      "This is not a problem, its just they have 'unclassified' somewhere in their lineages\n" +
                      "(e.g. Unclassified Bacteria) and this is not good for assessing placement accuracy.\n\n")

    # Write the intermediate classifications to a file
    # write_intermediate_assignments(inter_class_file, marker_assignments)

    if code_name not in marker_eval_inst.markers and marker not in marker_eval_inst.markers:
        logging.error("No sequences were classified as " + marker + ".\n")
        sys.exit(21)

    # For debugging:
    # for rank in marker_eval_inst.ranks:
    #     distal, pendant, tip = marker_eval_inst.summarize_rankwise_distances(rank)

    # Determine the specificity for each rank
    clade_exclusion_strings = get_classification_performance(marker_eval_inst)
    write_performance_table(args, performance_table, clade_exclusion_strings)
    summarize_taxonomic_diversity(marker_eval_inst)
    containment_strings = determine_containment(marker_eval_inst)
    write_containment_table(args, containment_table, containment_strings)


if __name__ == "__main__":
    main()
