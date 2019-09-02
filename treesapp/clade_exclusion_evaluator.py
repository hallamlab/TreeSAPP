#!/usr/bin/env python3

__author__ = 'Connor Morgan-Lang'

import os
import shutil
from ete3 import Tree
from glob import glob

from .fasta import write_new_fasta, read_fasta_to_dict
from .utilities import return_sequence_info_groups
from .external_command_interface import launch_write_command
from .file_parsers import tax_ids_file_to_leaves
from .classy import get_header_format, Evaluator, MarkerBuild
from .entrez_utils import *

_RANK_DEPTH_MAP = {0: "Cellular organisms", 1: "Kingdom",
                   2: "Phylum", 3: "Class", 4: "Order",
                   5: "Family", 6: "Genus", 7: "Species", 8: "Strain"}


def load_rank_depth_map(evaluator: Evaluator):
    evaluator.rank_depth_map = _RANK_DEPTH_MAP


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


def determine_containment(marker_eval_inst):
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
    # Set up collection for this analysis
    lowest_rank = ""
    for depth in sorted(_RANK_DEPTH_MAP):
        if _RANK_DEPTH_MAP[depth] in marker_eval_inst.ranks:
            if marker_eval_inst.get_sensitivity(_RANK_DEPTH_MAP[depth])[1] > 0:
                lowest_rank = _RANK_DEPTH_MAP[depth]
    containment_strings = list()
    n_queries, n_classified, _ = marker_eval_inst.get_sensitivity(lowest_rank)

    sys.stdout.write("Accuracy of " + str(n_classified) + ' ' + marker_eval_inst.target_marker.cog +
                     " classifications at " + lowest_rank + ":\n" +
                     "\tRank\tCorrect\tToo Shallow (%)\n")

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
                    q_accession = return_sequence_info_groups(sequence_info, header_db, original_header).accession
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
                mapped = False
                for treesapp_id, ref_seq in ref_sequences.items():
                    if ref_seq.accession == query:
                        lineage_assignments[marker][c_lineage].append(clean_lineage_string(ref_seq.lineage))
                        mapped = True
                        break
                if not mapped:
                    logging.error("Unable to map classified sequence '" + query + "' to a lineage.\n")
                    sys.exit(3)
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
    if re.search("unclassified", split_lineage[pos], re.IGNORECASE):
        return pos
    else:
        pos = get_unclassified_rank(pos+1, split_lineage)
    return pos


def pick_taxonomic_representatives(ref_seqs: dict, taxonomic_filter_stats: dict, max_cluster_size=5):
    """
    Removes queries with duplicate taxa - to prevent the taxonomic composition of the input
    from biasing the accuracy to over- or under-perform by classifying many sequences from very few groups.
    Also removes taxonomies with "*nclassified" in their lineage or are derived from environmental samples

    :param ref_seqs: A dictionary mapping accessions to lineages that need to be filtered
    :param taxonomic_filter_stats: A dictionary for tracking the number sequences filtered, those unique, etc.
    :param max_cluster_size: The maximum number of sequences representing a taxonomic cluster
    :return: dereplicated_lineages dict with lineages mapping to a (short) list of accessions
    """
    good_classified_lineages = dict()
    dereplicated_lineages = dict()
    num_rep_seqs = 0
    for tree_id, ref_seq in ref_seqs.items():
        query_taxonomy = clean_lineage_string(ref_seq.lineage)
        if len(query_taxonomy.split("; ")) < 5:
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

    if taxonomic_filter_stats["Unclassified"] == len(ref_seqs):
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

    logging.info("\t" + str(num_rep_seqs) + " representative sequences will be used for TreeSAPP evaluate analysis.\n")

    logging.debug("Representative sequence stats:\n\t" +
                  "Maximum representative sequences for a taxon " + str(taxonomic_filter_stats["Max"]) + "\n\t" +
                  "Minimum representative sequences for a taxon " + str(taxonomic_filter_stats["Min"]) + "\n\t" +
                  "Mean representative sequences for a taxon " + str(taxonomic_filter_stats["Mean"]) + "\n")

    return dereplicated_lineages, taxonomic_filter_stats


def str_list_index(lst: list, query: str):
    i = 0
    while i < len(lst):
        if query == lst[i]:
            break
        i += 1
    return i


def same_lineage(target_lineage: str, candidate_lineage: str) -> bool:
    """
    Identify which candidate lineages are belong to target lineage and are at least as specific
    1. "Root; Bacteria", "Root"                           False
    2. "Root; Bacteria", "Root; Archaea"                  False
    3. "Root; Bacteria", "Root; Bacteria"                 True
    4. "Root; Bacteria", "Root; Bacteria; Proteobacteria" True
    5. "Bacteria", "Root; Bacteria; Proteobacteria"       True
    Runs in O(n) complexity
    :param target_lineage:
    :param candidate_lineage:
    :return: Boolean
    """
    target_lineage = clean_lineage_string(target_lineage).split("; ")
    candidate_lineage = clean_lineage_string(candidate_lineage).split("; ")
    tlen = len(target_lineage)
    clen = len(candidate_lineage)
    k = 0  # For tracking the overlapping positions
    i = str_list_index(target_lineage, candidate_lineage[0])  # Lineage offset for the target
    j = str_list_index(candidate_lineage, target_lineage[0])  # Lineage offset for the candidate

    if i == tlen and j == clen:
        if set(target_lineage).intersection(set(candidate_lineage)):
            logging.debug("Target '%s' and candidate '%s' lineages were not overlapped " %
                          (target_lineage, candidate_lineage))
        return False

    # There is a common string
    while i+k < tlen and j+k < clen and target_lineage[i+k] == candidate_lineage[j+k]:
        k += 1

    # The verdict...
    if i+k == len(target_lineage):
        return True
    else:
        return False


def select_rep_seqs(deduplicated_assignments: dict, test_sequences: dict, target_lineage=None):
    """
    Function for creating a fasta-formatted dict from the accessions representing unique taxa in the test sequences

    :param deduplicated_assignments: dict of lineages mapped to a list of accessions
    :param test_sequences: list of ReferenceSequence objects
    :param target_lineage: A taxonomic lineage to filter the lineages (and sequences) in deduplicated assignments
    :return: Dictionary containing accessions as keys and sequences as values
    """
    if target_lineage:
        filtered_assignments = dict()
        for candidate_lineage in sorted(deduplicated_assignments):
            if same_lineage(target_lineage, candidate_lineage):
                filtered_assignments[candidate_lineage] = deduplicated_assignments[candidate_lineage]
    else:
        filtered_assignments = deduplicated_assignments

    deduplicated_fasta_dict = dict()
    for lineage in sorted(filtered_assignments):
        for accession in filtered_assignments[lineage]:
            matched = False
            for treesapp_id, ref_seq in test_sequences.items():
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


def correct_accession(description):
    header_format_re, header_db, header_molecule = get_header_format(description)
    sequence_info = header_format_re.match(description)
    return return_sequence_info_groups(sequence_info, header_db, description).accession


def prep_graftm_ref_files(treesapp_dir: str, intermediate_dir: str, target_taxon: str, marker: MarkerBuild, depth: int):
    """
    From the original TreeSAPP reference package files, the necessary GraftM create input files are generated
    with all reference sequences related to the target_taxon removed from the multiple sequence alignment,
    unaligned reference FASTA file and the tax_ids file.
    :param treesapp_dir: Path to the TreeSAPP reference package directory
    :param intermediate_dir:  Path to write the intermediate files with target references removed
    :param target_taxon: Name of the taxon that is being tested in the current clade exclusion iteration
    :param marker: MarkerBuild instance for the reference package being tested
    :param depth: Depth of the current taxonomic rank in hierarchy (e.g. Phylum = 2, Class = 3, etc.)
    :return: None
    """
    # Move the original FASTA, tree and tax_ids files to a temporary location
    marker_fa = os.sep.join([treesapp_dir, "data", "alignment_data", marker.cog + ".fa"])
    marker_tax_ids = os.sep.join([treesapp_dir, "data", "tree_data", "tax_ids_" + marker.cog + ".txt"])
    off_target_ref_leaves = dict()
    # tax_ids file
    ref_tree_leaves = tax_ids_file_to_leaves(marker_tax_ids)
    with open(intermediate_dir + "tax_ids_" + marker.cog + ".txt", 'w') as tax_ids_handle:
        tax_ids_strings = list()
        for ref_leaf in ref_tree_leaves:
            c_lineage = clean_lineage_string(ref_leaf.lineage)
            if re.search(target_taxon, c_lineage):
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
    write_new_fasta(accession_fasta_dict, intermediate_dir + marker.cog + ".mfa")
    for acc in accession_fasta_dict:
        accession_fasta_dict[acc] = re.sub('-', '', accession_fasta_dict[acc])
    write_new_fasta(accession_fasta_dict, intermediate_dir + marker.cog + ".fa")
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
    target_clade = clean_lineage_string(target_clade, ["Root; "])
    with open(marker_tax_ids, 'w') as tax_ids_handle:
        for ref_leaf in ref_tree_leaves:
            tax_ids_string = ""
            c_lineage = clean_lineage_string(ref_leaf.lineage, ["Root; "])
            sc_lineage = c_lineage.split("; ")
            if len(sc_lineage) < depth:
                n_shallow += 1
                continue
            if target_clade == '; '.join(sc_lineage[:depth+1]):
                n_match += 1
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

    logging.debug("Reference sequence filtering stats for " + target_clade + "\n" +
                  "\n".join(["Match taxon\t" + str(n_match),
                             "Unclassified\t" + str(n_unclassified),
                             "Too shallow\t" + str(n_shallow),
                             "Remaining\t" + str(len(off_target_ref_leaves))]) + "\n")

    # fasta
    ref_fasta_dict = read_fasta_to_dict(marker_fa)
    off_target_headers = [num_id + '_' + marker for num_id in off_target_ref_leaves]
    if len(off_target_headers) == 0:
        logging.error("No reference sequences were retained for building testing " + target_clade + "\n")
        sys.exit(19)
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
    remnants = glob(intermediate_dir + "ORIGINAL*", recursive=True)
    if not remnants:
        logging.warning("No remaining original reference package files found for destruction.\n")
    for tmp_file in remnants:
        os.remove(tmp_file)
    return


def restore_reference_package(treesapp_dir, prefix, output_dir, marker):
    """
      Prepares TreeSAPP tree, alignment and taxonomic identification map (tax_ids) files for clade exclusion analysis,
    and performs classification with TreeSAPP

    :return: The paths to the classification table and the taxon-excluded tax_ids file
    """
    # The edited tax_ids file with clade excluded is required for performance analysis
    # Copy the edited, clade-excluded tax_ids file to the output directory
    shutil.copy(os.sep.join([treesapp_dir, "data", "tree_data", "tax_ids_" + marker + ".txt"]), output_dir)

    # Move the original FASTA, tree and tax_ids files back to the proper directories
    shutil.copy(prefix + "_tree.txt", os.sep.join([treesapp_dir, "data", "tree_data", marker + "_tree.txt"]))
    if os.path.isfile(prefix + "_bipartitions.txt"):
        shutil.copy(prefix + "_bipartitions.txt",
                    os.sep.join([treesapp_dir, "data", "tree_data", marker + "_bipartitions.txt"]))
    shutil.copy(prefix + "_tax_ids.txt", os.sep.join([treesapp_dir, "data", "tree_data", "tax_ids_" + marker + ".txt"]))
    shutil.copy(prefix + ".fa", os.sep.join([treesapp_dir, "data", "alignment_data", marker + ".fa"]))
    shutil.copy(prefix + ".hmm", os.sep.join([treesapp_dir, "data", "hmm_data", marker + ".hmm"]))

    return


def build_graftm_package(gpkg_path: str, tax_file: str, mfa_file: str, fa_file: str, threads: int):
    create_command = ["graftM", "create"]
    create_command += ["--threads", str(threads)]
    create_command += ["--alignment", mfa_file]
    create_command += ["--sequences", fa_file]
    create_command += ["--taxonomy", tax_file]
    create_command += ["--output", gpkg_path]
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
