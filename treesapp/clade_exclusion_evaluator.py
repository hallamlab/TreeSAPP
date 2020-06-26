#!/usr/bin/env python3

__author__ = 'Connor Morgan-Lang'

import os
import shutil
import logging
import sys
import re
from glob import glob

from treesapp.fasta import write_new_fasta, read_fasta_to_dict, load_fasta_header_regexes, sequence_info_groups
from treesapp.external_command_interface import launch_write_command
from treesapp.classy import get_header_format, Evaluator
from treesapp.refpkg import ReferencePackage
from treesapp.phylo_dist import trim_lineages_to_rank
from treesapp.entrez_utils import EntrezRecord

_RANK_DEPTH_MAP = {0: "cellular organisms", 1: "domain",
                   2: "phylum", 3: "class", 4: "order",
                   5: "family", 6: "genus", 7: "species", 8: "strain"}


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
                assignments_string += marker + "\t" + ref + "\t"
                assignments_string += ','.join(queries) + "\n"

    file_handler.write(assignments_string)
    file_handler.close()
    logging.debug("done.\n")
    return


def determine_containment(marker_eval_inst: Evaluator):
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
    for depth in sorted(marker_eval_inst.rank_depth_map):
        if marker_eval_inst.rank_depth_map[depth] in marker_eval_inst.ranks:
            if marker_eval_inst.get_sensitivity(marker_eval_inst.rank_depth_map[depth])[1] > 0:
                lowest_rank = marker_eval_inst.rank_depth_map[depth]
    containment_strings = list()
    n_queries, n_classified, _ = marker_eval_inst.get_sensitivity(lowest_rank)

    sys.stdout.write("Accuracy of " + str(n_classified) + ' ' + marker_eval_inst.ref_pkg.prefix +
                     " classifications at " + lowest_rank + ":\n" +
                     "\tRank\tCorrect\tToo Shallow (%)\n")

    # Begin parsing through the depths
    if lowest_rank == "root":
        return
    for depth in sorted(marker_eval_inst.rank_depth_map):
        incorrect_assignments = dict()
        rank = marker_eval_inst.rank_depth_map[depth]
        if rank in ["root", "species", "strain"]:
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
            cleaned_assignments[marker][ref] = list()
            for query in queries:
                try:
                    query = re.sub('_', "; ", query)
                except TypeError:
                    logging.error("Unexpected type of lineage string (" + str(type(query)) + "):\n" +
                                  str(query) + "\n")
                    sys.exit(21)
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
    header_regex_dict = load_fasta_header_regexes()
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
                    header_format_re, header_db, header_molecule = get_header_format(original_header, header_regex_dict)
                    sequence_info = header_format_re.match(original_header)
                    q_accession = sequence_info_groups(sequence_info, header_db, original_header).accession
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


def map_headers_to_lineage(assignments: dict, ref_sequences: dict) -> dict:
    """
    The alternative function to map_full_headers. Using ReferenceSequence objects,

    :param assignments: A dictionary of ReferencePackage.prefix keys mapped to a dictionary of assigned lineage keys
    mapped to the headers of sequences assigned as their respective lineages.
    :param ref_sequences:
    :return:  A dictionary mapping ReferencePackage.prefix's to a dictionary of assigned lineage keys to a list of
    reference sequences (values). An example is
    """
    lineage_assignments = dict()
    for refpkg_name in assignments:
        lineage_assignments[refpkg_name] = dict()
        for assigned_lineage in assignments[refpkg_name].keys():
            classified_headers = assignments[refpkg_name][assigned_lineage]
            lineage_assignments[refpkg_name][assigned_lineage] = list()
            for query in classified_headers:
                mapped = False
                for treesapp_id, ref_seq in ref_sequences.items():  # type: (int, EntrezRecord)
                    if ref_seq.accession == query:
                        lineage_assignments[refpkg_name][assigned_lineage].append(ref_seq.lineage)
                        mapped = True
                        break
                if not mapped:
                    logging.error("Unable to map classified sequence '" + query + "' to a lineage.\n")
                    sys.exit(3)
            if len(lineage_assignments[refpkg_name][assigned_lineage]) > len(classified_headers):
                logging.error(str(len(classified_headers)) + " accessions mapped to " +
                              str(len(lineage_assignments[refpkg_name][assigned_lineage])) + " lineages.\n")
                sys.exit(21)
            elif len(lineage_assignments[refpkg_name][assigned_lineage]) < len(classified_headers):
                logging.debug(str(len(classified_headers)) + " accessions mapped to " +
                              str(len(lineage_assignments[refpkg_name][assigned_lineage])) + " lineages.\n")
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


def get_testable_lineages_for_rank(ref_lineage_map: dict, query_lineage_map: dict, rank: str) -> list:
    """


    :param ref_lineage_map: A dictionary mapping all reference sequences to their respective taxonomic lineages
    :param query_lineage_map: A dictionary mapping representative query sequences to taxonomic lineages
    :param rank: Name of a taxonomic rank to test and is used to guide taxonomic lineage trimming with
     {"Kingdom": 1, "Phylum": 2, "Class": 3, "Order": 4, "Family": 5, "Genus": 6, "Species": 7}
    :return: List of lineages where the optimal rank exists in the reference tree after clade exclusion
    """
    lineages = []
    leaf_trimmed_taxa_map = trim_lineages_to_rank(ref_lineage_map, rank)
    unique_ref_lineages = sorted(set(leaf_trimmed_taxa_map.values()))
    unique_query_lineages = sorted(set(trim_lineages_to_rank(query_lineage_map, rank).values()))
    for lineage in unique_query_lineages:
        # Is the optimal placement in the pruned reference tree?
        optimal_lca_taxonomy = "; ".join(lineage.split("; ")[:-1])
        if optimal_lca_taxonomy not in ["; ".join(tl.split("; ")[:-1]) for tl in unique_ref_lineages
                                        if tl != lineage]:
            logging.debug("Optimal placement target '" + optimal_lca_taxonomy + "' not in pruned tree.\n")
        else:
            lineages.append(lineage)
    return lineages


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
        if len(ref_seq.lineage.split("; ")) < 5:
            continue
        if ref_seq.lineage not in good_classified_lineages:
            good_classified_lineages[ref_seq.lineage] = list()
        if re.search("unclassified|environmental sample", ref_seq.lineage, re.IGNORECASE):
            # # Remove taxonomic lineages that are unclassified at the Phylum level or higher
            # unclassified_depth = get_unclassified_rank(0, query_taxonomy.split("; "))
            # if unclassified_depth > 4:
            #     good_classified_lineages[query_taxonomy].append(ref_seq.accession)
            # else:
            taxonomic_filter_stats["Unclassified"] += 1
        else:
            taxonomic_filter_stats["Classified"] += 1
            good_classified_lineages[ref_seq.lineage].append(ref_seq.accession)

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
    target_lineage = target_lineage.split("; ")
    candidate_lineage = candidate_lineage.split("; ")
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
            for treesapp_id, ref_seq in test_sequences.items():  # type: (str, EntrezRecord)
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


def prep_graftm_ref_files(intermediate_dir: str, target_taxon: str, refpkg: ReferencePackage, depth: int):
    """
    From the original TreeSAPP reference package files, the necessary GraftM create input files are generated
    with all reference sequences related to the target_taxon removed from the multiple sequence alignment,
    unaligned reference FASTA file and the tax_ids file.

    :param intermediate_dir:  Path to write the intermediate files with target references removed
    :param target_taxon: Name of the taxon that is being tested in the current clade exclusion iteration
    :param refpkg: A ReferencePackage instance for the reference package being tested
    :param depth: Depth of the current taxonomic rank in hierarchy (e.g. Phylum = 2, Class = 3, etc.)
    :return: None
    """
    # Move the original FASTA, tree and tax_ids files to a temporary location
    off_target_ref_leaves = dict()
    # tax_ids file
    ref_tree_leaves = refpkg.generate_tree_leaf_references_from_refpkg()
    with open(intermediate_dir + "tax_ids_" + refpkg.prefix + ".txt", 'w') as tax_ids_handle:
        tax_ids_strings = list()
        for ref_leaf in ref_tree_leaves:
            if re.search(target_taxon, ref_leaf.lineage):
                continue
            sc_lineage = ref_leaf.lineage.split("; ")
            if len(sc_lineage) < depth:
                continue
            if re.search("unclassified|environmental sample", ref_leaf.lineage, re.IGNORECASE):
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
            tax_ids_strings.append(accession + "\t" + ref_leaf.lineage)
        tax_ids_handle.write("\n".join(tax_ids_strings) + "\n")

    # fasta
    ref_fasta_dict = read_fasta_to_dict(refpkg.msa)
    accession_fasta_dict = dict()
    for key_id in ref_fasta_dict:
        num_key = key_id.split('_')[0]
        if num_key in off_target_ref_leaves.keys():
            accession_fasta_dict[off_target_ref_leaves[num_key]] = ref_fasta_dict[key_id]
        else:
            pass
    write_new_fasta(accession_fasta_dict, intermediate_dir + refpkg.prefix + ".mfa")
    for acc in accession_fasta_dict:
        accession_fasta_dict[acc] = re.sub('-', '', accession_fasta_dict[acc])
    write_new_fasta(accession_fasta_dict, intermediate_dir + refpkg.prefix + ".fa")
    return


def validate_ref_package_files(refpkg: ReferencePackage, intermediate_dir):
    intermediate_prefix = intermediate_dir + "ORIGINAL"

    # This prevents users from quitting mid-analysis then using bad reference package files with clades removed
    remnants = glob(intermediate_prefix + "*")
    if len(remnants) > 0:
        logging.warning("Intermediate files from incomplete analysis exist:\n" + "\n".join(remnants) + "\n")
        logging.info("Attempting recovery... ")
        # Easiest way to recover is to move the originals back to proper location and proceed
        try:
            shutil.copy(intermediate_prefix + ".fa", refpkg.msa)
            shutil.copy(intermediate_prefix + "_tree.txt", refpkg.tree)
            if os.path.isfile(intermediate_prefix + "_bipartitions.txt"):
                shutil.copy(intermediate_prefix + "_bipartitions.txt", refpkg.boot_tree)
            shutil.copy(intermediate_prefix + "_tax_ids.txt", refpkg.lineage_ids)
            logging.info("succeeded.\n")
        except FileNotFoundError:
            logging.info("failed. Redownload data files and start over.\n")
            sys.exit(21)
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
