#!/usr/bin/env python3

import os
import sys
import argparse
import logging
import re
from ete3 import Tree
import numpy as np
from glob import glob

from fasta import read_fasta_to_dict, write_new_fasta, deduplicate_fasta_sequences, trim_multiple_alignment,\
    get_headers, summarize_fasta_sequences, format_read_fasta
from file_parsers import tax_ids_file_to_leaves, read_stockholm_to_dict, validate_alignment_trimming,\
    parse_ref_build_params, multiple_alignment_dimensions, parse_domain_tables
import utilities
from entrez_utils import read_accession_taxa_map, get_multiple_lineages, build_entrez_queries, \
    write_accession_lineage_map, verify_lineage_information, \
    entrez_records_to_accessions, entrez_records_to_accession_lineage_map
from phylo_dist import trim_lineages_to_rank, cull_outliers, parent_to_tip_distances, regress_ranks
from external_command_interface import setup_progress_bar
from jplace_utils import jplace_parser
from classy import prep_logging, register_headers, get_header_info, ReferencePackage
from entish import map_internal_nodes_leaves

__author__ = 'Connor Morgan-Lang'


def get_options():
    parser = argparse.ArgumentParser(description="Workflow for estimating calibrating the edge distances corresponding"
                                                 " to taxonomic ranks by iterative leave-one-out validation")
    required_args = parser.add_argument_group("Required arguments")
    seqop_args = parser.add_argument_group("Sequence operation arguments")
    taxa_args = parser.add_argument_group("Taxonomic-lineage arguments")
    miscellaneous_opts = parser.add_argument_group("Miscellaneous arguments")

    required_args.add_argument("-f", "--fasta_input", required=True,
                               help='The raw, unclustered and unfiltered FASTA file to train the reference package.')
    required_args.add_argument("-n", "--name", required=True,
                               help="Prefix name of the reference package "
                                    "(i.e. McrA for McrA.fa, McrA.hmm, McrA_tree.txt)")
    required_args.add_argument("-p", "--pkg_path", required=True,
                               help="The path to the TreeSAPP-formatted reference package.")
    seqop_args.add_argument("-d", "--domain",
                            help="An HMM profile representing a specific domain.\n"
                                 "Domains will be excised from input sequences based on hmmsearch alignments.",
                            required=False, default=None)
    taxa_args.add_argument("-l", "--lineages",
                           help="The accession lineage map downloaded during reference package generation.",
                           required=False)
    miscellaneous_opts.add_argument('-m', '--molecule', default='prot', choices=['prot', 'dna', 'rrna'],
                                    help='the type of input sequences (prot = Protein [DEFAULT]; dna = Nucleotide )')
    miscellaneous_opts.add_argument("-T", "--num_threads", required=False, default=4, type=int,
                                    help="The number of threads to be used by RAxML.")
    miscellaneous_opts.add_argument("-o", "--output_dir", required=False, default='.',
                                    help="Path to directory for writing outputs.")
    miscellaneous_opts.add_argument("-v", "--verbose", action='store_true',  default=False,
                                    help='Prints a more verbose runtime log')
    miscellaneous_opts.add_argument("-O", "--overwrite", default=False, action="store_true",
                                    help="Force recalculation of placement distances for query sequences.")
    args = parser.parse_args()
    args.targets = ["ALL"]
    return args


class PQuery:
    def __init__(self, lineage_str, rank_str):
        # Inferred from JPlace file
        self.pendant = 0.0
        self.mean_tip = 0.0
        self.distal = 0.0
        self.likelihood = 0.0
        self.lwr = 0.0
        self.inode = ""
        self.parent_node = ""
        self.name = ""

        # Known from outer scope
        self.lineage = lineage_str
        self.rank = rank_str

    def total_distance(self):
        return round(sum([self.pendant, self.mean_tip, self.distal]), 5)

    def summarize_placement(self):
        summary_string = "Placement of " + self.name + " at rank " + self.rank + \
                         ":\nLineage = " + self.lineage + \
                         "\nInternal node = " + str(self.inode) + \
                         "\nDistances:" + \
                         "\n\tDistal = " + str(self.distal) +\
                         "\n\tPendant = " + str(self.pendant) +\
                         "\n\tTip = " + str(self.mean_tip) +\
                         "\nLikelihood = " + str(self.likelihood) +\
                         "\nL.W.R. = " + str(self.lwr) +\
                         "\n"
        return summary_string


def write_placement_table(pqueries, placement_table_file, marker):
    header = ["Marker", "Rank", "Lineage", "Query.Name", "Internal.Node", "Placement.LWR", "Tree.Likelihood",
              "Dist.Distal", "Dist.Pendant", "Dist.MeanTip", "Dist.Total"]
    placement_info_strs = list()
    for pquery in pqueries:
        if pquery:
            placement_info_strs.append("\t".join(
                [marker, str(pquery.rank), str(pquery.lineage), str(pquery.name), str(pquery.inode),
                 str(pquery.lwr), str(pquery.likelihood),
                 str(pquery.distal), str(pquery.pendant), str(pquery.mean_tip), str(pquery.total_distance())])
            )

    with open(placement_table_file, 'w') as file_handler:
        file_handler.write('#' + "\t".join(header) + "\n")
        file_handler.write("\n".join(placement_info_strs) + "\n")
    return


def rarefy_rank_distances(rank_distances):
    rarefied_dists = dict()
    min_samples = min([len(rank_distances[rank]) for rank in rank_distances])
    for rank in rank_distances:
        slist = sorted(rank_distances[rank])
        if len(slist) == min_samples:
            rarefied_dists[rank] = slist
        else:
            rarefied_dists[rank] = list()
            i = 0
            while i < min_samples:
                rarefied_dists[rank].append(slist.pop(np.random.randint(0, len(slist))))
                i += 1
    return rarefied_dists


def read_placement_summary(placement_summary_file):
    """
    Reads a specially-formatted file and returns the rank-wise clade-exclusion placement distances
    :param placement_summary_file:
    :return:
    """
    taxonomic_placement_distances = dict()
    with open(placement_summary_file, 'r') as place_summary:
        rank = ""
        line = place_summary.readline()
        while line:
            line = line.strip()
            if line:
                if line[0] == '#':
                    rank = line.split(' ')[1]
                elif line[0] == '[':
                    dist_strings = re.sub('\[|\]', '', line).split(", ")
                    dists = [float(dist) for dist in dist_strings]
                    if len(dists) > 1:
                        taxonomic_placement_distances[rank] = dists
            line = place_summary.readline()
    return taxonomic_placement_distances


def complete_regression(taxonomic_placement_distances, taxonomic_ranks=None):
    """
    Wrapper for performing outlier removal, normalization via rarefaction, and regression

    :param taxonomic_placement_distances:
    :param taxonomic_ranks: A dictionary mapping rank names (e.g. Phylum)
    to rank depth values where Kingdom is 0, Phylum is 1, etc.
    :return:
    """
    if not taxonomic_placement_distances:
        return []

    if not taxonomic_ranks:
        taxonomic_ranks = {"Phylum": 1, "Class": 2, "Order": 3, "Family": 4, "Genus": 5, "Species": 6, "Strain": 7}

    filtered_pds = dict()
    for rank in taxonomic_placement_distances:
        init_s = len(list(taxonomic_placement_distances[rank]))
        if init_s <= 3:
            logging.warning("Insufficient placement distance samples (" + str(init_s) + ") for " + rank + ".\n")
            return []
        # print(rank, "raw", np.median(list(taxonomic_placement_distances[rank])))
        filtered_pds[rank] = cull_outliers(list(taxonomic_placement_distances[rank]))
        # print(rank, "filtered", np.median(list(filtered_pds[rank])))
        if len(filtered_pds[rank]) == 0:
            logging.warning("Ranks have 0 samples after filtering outliers.\n")
            return []

    # Rarefy the placement distances to the rank with the fewest samples
    rarefied_pds = rarefy_rank_distances(filtered_pds)
    for rank in rarefied_pds:
        if len(rarefied_pds[rank]) == 0:
            logging.warning("Ranks have 0 samples after rarefaction.\n")
            return []
    # for rank in rarefied_pds:
    #     print(rank, "rarefied", np.median(list(rarefied_pds[rank])))
    return regress_ranks(rarefied_pds, taxonomic_ranks)


def prepare_training_data(fasta_input: str, output_dir: str, executables: dict,
                          leaf_taxa_map: dict, accession_lineage_map: dict, taxonomic_ranks):
    """
    Function for creating a non-redundant inventory of sequences to be used for training the rank-placement distance
    linear model. Removes sequences that share an identical accession, are more than 95% similar and limits the
    number of taxonomically-identical sequences to 30.

    :param fasta_input: Path to a FASTA-formatted file
    :param output_dir: Path to write intermediate output files (such as UCLUST outputs)
    :param executables: A dictionary mapping software to a path of their respective executable
    :param leaf_taxa_map: A dictionary mapping TreeSAPP numeric sequence identifiers to taxonomic lineages
    :param accession_lineage_map: A dictionary mapping NCBI accession IDs to full NCBI taxonomic lineages
    :param taxonomic_ranks: A dictionary mapping rank names (e.g. Phylum)
     to rank depth values where Kingdom is 0, Phylum is 1, etc.
    :return: A dictionary storing the sequence names being used to test each taxon within each rank
    """
    rank_training_seqs = dict()
    optimal_placement_missing = list()
    taxon_training_queries = list()
    similarity = 0.97  # The proportional similarity to cluster the training sequences
    max_reps = 30  # The maximum number of representative sequences from a specific taxon for training
    uclust_prefix = output_dir + os.sep + "uclust" + str(similarity)

    # Cluster the training sequences to mitigate harmful redundancy
    utilities.cluster_sequences(executables["usearch"], fasta_input, uclust_prefix, similarity)
    uclust_fasta_dict = read_fasta_to_dict(uclust_prefix + ".fa")
    logging.debug("\t" + str(len(uclust_fasta_dict.keys())) + " sequence clusters\n")

    logging.info("Preparing deduplicated sequence set for training... ")
    # Remove sequences with duplicate accessions
    nr_fasta_dict = deduplicate_fasta_sequences(uclust_fasta_dict)
    uclust_fasta_dict.clear()
    for seq_name in nr_fasta_dict.keys():
        uclust_fasta_dict[seq_name.split(" ")[0]] = nr_fasta_dict[seq_name]
    nr_fasta_dict.clear()

    # Determine the set of reference sequences to use at each rank
    for rank in taxonomic_ranks:
        rank_training_seqs[rank] = dict()
        leaf_trimmed_taxa_map = trim_lineages_to_rank(leaf_taxa_map, rank)
        unique_taxonomic_lineages = sorted(set(leaf_trimmed_taxa_map.values()))

        # Remove all sequences belonging to a taxonomic rank from tree and reference alignment
        for taxonomy in unique_taxonomic_lineages:
            optimal_lca_taxonomy = "; ".join(taxonomy.split("; ")[:-1])
            if optimal_lca_taxonomy not in ["; ".join(tl.split("; ")[:-1]) for tl in unique_taxonomic_lineages if
                                            tl != taxonomy]:
                optimal_placement_missing.append(optimal_lca_taxonomy)
            else:
                for seq_name in sorted(accession_lineage_map):
                    # Not all keys in accession_lineage_map are in fasta_dict (duplicate sequences were removed)
                    if re.search(taxonomy,
                                 utilities.clean_lineage_string(accession_lineage_map[seq_name])) and seq_name in uclust_fasta_dict:
                        taxon_training_queries.append(seq_name)
                    if len(taxon_training_queries) == max_reps:
                        break
                if len(taxon_training_queries) > 0:
                    rank_training_seqs[rank][taxonomy] = list(taxon_training_queries)
                    taxon_training_queries.clear()
    logging.info("done.\n")

    logging.debug("Optimal placement target is not found in the pruned tree for following taxa:\n\t" +
                  "\n\t".join(optimal_placement_missing) + "\n")

    return rank_training_seqs, uclust_fasta_dict


def train_placement_distances(rank_training_seqs: dict, taxonomic_ranks: dict,
                              ref_fasta_dict: dict, dedup_fasta_dict: dict,
                              ref_pkg: ReferencePackage, leaf_taxa_map: dict,
                              molecule: str, executables: dict, raxml_threads=4):
    """
    Function for iteratively performing leave-one-out analysis for every taxonomic lineage represented in the tree,
    yielding an estimate of placement distances corresponding to taxonomic ranks.

    :param rank_training_seqs: A dictionary storing the sequence names being used to test each taxon within each rank
    :param taxonomic_ranks: A dictionary mapping rank names (e.g. Phylum)
     to rank depth values where Kingdom is 0, Phylum is 1, etc.
    :param ref_fasta_dict: A dictionary with headers as keys and sequences as values containing only reference sequences
    :param dedup_fasta_dict: Dictionary with headers as keys and sequences as values for deduplicated training sequences
    :param ref_pkg: A ReferencePackage instance
    :param leaf_taxa_map: A dictionary mapping TreeSAPP numeric sequence identifiers to taxonomic lineages
    :param executables: A dictionary mapping software to a path of their respective executable
    :param molecule: Molecule type [prot | dna | rrna]
    :param raxml_threads: Number of threads to be used by RAxML for parallel computation

    :return:
    """

    logging.info("\nEstimating branch-length placement distances for taxonomic ranks. Progress:\n")
    taxonomic_placement_distances = dict()
    taxonomy_filtered_query_seqs = dict()
    pruned_ref_fasta_dict = dict()
    query_seq_name_map = dict()
    seq_dict = dict()
    pqueries = list()
    intermediate_files = list()
    aligner = "hmmalign"

    temp_tree_file = "tmp_tree.txt"
    temp_ref_aln_prefix = "taxonomy_filtered_ref_seqs"
    temp_query_fasta_file = "queries.fasta"
    query_multiple_alignment = aligner + "_queries_aligned.phy"

    # Read the tree as ete3 Tree instance
    ref_tree = Tree(ref_pkg.tree)

    bmge_file = executables["BMGE.jar"]
    if not os.path.exists(bmge_file):
        raise FileNotFoundError("Cannot find " + bmge_file)

    num_training_queries = 0
    for rank in rank_training_seqs:
        num_rank_training_seqs = 0
        for taxonomy in rank_training_seqs[rank]:
            num_rank_training_seqs += len(rank_training_seqs[rank][taxonomy])
        if len(rank_training_seqs[rank]) == 0:
            logging.error("No sequences available for estimating " + rank + "-level placement distances.\n")
            return taxonomic_placement_distances, pqueries
        else:
            logging.debug(str(num_rank_training_seqs) + " sequences to train " + rank + "-level placement distances\n")
        num_training_queries += num_rank_training_seqs

    if num_training_queries < 30:
        logging.error("Too few (" + str(num_training_queries) + ") sequences for training placement distance model.\n")
        sys.exit(19)
    if num_training_queries < 50:
        logging.warning("Only " + str(num_training_queries) + " sequences for training placement distance model.\n")
    step_proportion = setup_progress_bar(num_training_queries)
    acc = 0.0

    # For each rank from Class to Species (Kingdom & Phylum-level classifications to be inferred by LCA):
    for rank in sorted(rank_training_seqs, reverse=True):
        if rank not in taxonomic_ranks:
            logging.error("Rank '" + rank + "' not found in ranks being used for training.\n")
            sys.exit(33)
        taxonomic_placement_distances[rank] = list()
        leaf_trimmed_taxa_map = trim_lineages_to_rank(leaf_taxa_map, rank)
        
        # Add the lineages to the Tree instance
        for leaf in ref_tree:
            leaf.add_features(lineage=leaf_trimmed_taxa_map.get(leaf.name, "none"))

        # Remove all sequences belonging to a taxonomic rank from tree and reference alignment
        for taxonomy in sorted(rank_training_seqs[rank]):
            logging.debug("Testing placements for " + taxonomy + ":\n")
            query_name = re.sub(r"([ /])", '_', taxonomy.split("; ")[-1])
            leaves_excluded = 0

            # Write query FASTA containing sequences belonging to `taxonomy`
            query_seq_decrementor = -1
            for seq_name in rank_training_seqs[rank][taxonomy]:
                query_seq_name_map[query_seq_decrementor] = seq_name
                taxonomy_filtered_query_seqs[str(query_seq_decrementor)] = dedup_fasta_dict[seq_name]
                query_seq_decrementor -= 1
            logging.debug("\t" + str(len(taxonomy_filtered_query_seqs.keys())) + " query sequences.\n")
            acc += len(taxonomy_filtered_query_seqs.keys())
            write_new_fasta(taxonomy_filtered_query_seqs, fasta_name=temp_query_fasta_file)
            intermediate_files.append(temp_query_fasta_file)

            for key in ref_fasta_dict.keys():
                node = key.split('_')[0]
                # Node with truncated and/or unclassified lineages are not in `leaf_trimmed_taxa_map`
                if node in leaf_trimmed_taxa_map and not re.match(taxonomy, leaf_trimmed_taxa_map[node]):
                    pruned_ref_fasta_dict[node] = ref_fasta_dict[key]
                else:
                    leaves_excluded += 1

            unique_ref_headers = set([re.sub('_' + re.escape(ref_pkg.prefix), '', x) for x in pruned_ref_fasta_dict.keys()])
            logging.debug("\t" + str(leaves_excluded) + " sequences pruned from tree.\n")

            # Copy the tree since we are removing leaves of `taxonomy` and don't want this to be permanent
            tmp_tree = ref_tree.copy(method="deepcopy")
            # iteratively detaching the monophyletic clades generates a bad tree, so do it all at once
            tmp_tree.prune(pruned_ref_fasta_dict.keys(), preserve_branch_length=True)
            # Resolve any multifurcations
            tmp_tree.resolve_polytomy()
            logging.debug("\t" + str(len(tmp_tree.get_leaves())) + " leaves in pruned tree.\n")

            # Write the new reference tree with sequences from `taxonomy` removed
            tmp_tree.write(outfile=temp_tree_file, format=5)
            intermediate_files.append(temp_tree_file)

            ##
            # Run hmmalign, BMGE and RAxML to map sequences from the taxonomic rank onto the tree
            ##
            if aligner == "papara":
                temp_ref_phylip_file = temp_ref_aln_prefix + ".phy"
                # Write the reference MSA with sequences of `taxonomy` removed
                phy_dict = utilities.reformat_fasta_to_phy(pruned_ref_fasta_dict)
                utilities.write_phy_file(temp_ref_phylip_file, phy_dict)
                aln_stdout = utilities.run_papara(executables["papara"],
                                                  temp_tree_file, temp_ref_phylip_file, temp_query_fasta_file,
                                        "prot")
                intermediate_files.append(temp_ref_phylip_file)
                os.rename("papara_alignment.default", query_multiple_alignment)
            elif aligner == "hmmalign":
                temp_ref_fasta_file = temp_ref_aln_prefix + ".fasta"
                temp_ref_profile = temp_ref_aln_prefix + ".hmm"
                sto_file = re.sub("\.phy$", ".sto", query_multiple_alignment)
                # Write the pruned reference FASTA file
                write_new_fasta(pruned_ref_fasta_dict, temp_ref_fasta_file)
                # Build the HMM profile that doesn't include pruned reference sequences
                utilities.build_hmm_profile(executables["hmmbuild"], temp_ref_fasta_file, temp_ref_profile)
                # Currently not supporting rRNA references (phylogenetic_rRNA)
                aln_stdout = utilities.profile_aligner(executables, temp_ref_fasta_file, temp_ref_profile,
                                                       temp_query_fasta_file, sto_file)
                # Reformat the Stockholm format created by cmalign or hmmalign to Phylip
                sto_dict = read_stockholm_to_dict(sto_file)
                for seq_name in sto_dict:
                    try:
                        int(seq_name.split('_')[0])
                        seq_dict[seq_name.split('_')[0]] = sto_dict[seq_name]
                    except ValueError:
                        seq_dict[seq_name] = sto_dict[seq_name]
                write_new_fasta(seq_dict, query_multiple_alignment)
                intermediate_files += [temp_ref_fasta_file, temp_ref_profile, sto_file, query_multiple_alignment]
            else:
                logging.error("Unrecognised alignment tool '" + aligner + "'. Exiting now.\n")
                sys.exit(33)
            logging.debug(str(aln_stdout) + "\n")

            query_filtered_multiple_alignment = trim_multiple_alignment(bmge_file, query_multiple_alignment,
                                                                        molecule, "BMGE")
            intermediate_files += glob(query_filtered_multiple_alignment + "*")

            # Ensure reference sequences haven't been removed
            msa_dict, summary_str = validate_alignment_trimming([query_filtered_multiple_alignment], unique_ref_headers, True)
            nrow, ncolumn = multiple_alignment_dimensions(seq_dict=read_fasta_to_dict(query_filtered_multiple_alignment),
                                                          mfa_file=query_filtered_multiple_alignment)
            logging.debug("Columns = " + str(ncolumn) + "\n")
            if query_filtered_multiple_alignment not in msa_dict.keys():
                logging.debug("Placements for '" + taxonomy + "' are being skipped after failing MSA validation.\n")
                for old_file in intermediate_files:
                    os.remove(old_file)
                    intermediate_files.clear()
                continue
            logging.debug("Number of sequences discarded: " + summary_str + "\n")

            # Run RAxML with the parameters specified
            raxml_files = utilities.raxml_evolutionary_placement(executables["raxmlHPC"], temp_tree_file,
                                                                 query_filtered_multiple_alignment, ref_pkg.sub_model, "./",
                                                                 query_name, raxml_threads)

            # Parse the JPlace file to pull distal_length+pendant_length for each placement
            jplace_data = jplace_parser(raxml_files["jplace"])
            placement_tree = jplace_data.tree
            node_map = map_internal_nodes_leaves(placement_tree)
            for pquery in jplace_data.placements:
                top_lwr = 0.1
                top_placement = PQuery(taxonomy, rank)
                for name, info in pquery.items():
                    if name == 'p':
                        for placement in info:
                            # Only record the best placement's distance
                            lwr = float(placement[2])
                            if lwr > top_lwr:
                                top_lwr = lwr
                                top_placement.inode = placement[0]
                                top_placement.likelihood = placement[1]
                                top_placement.lwr = lwr
                                top_placement.distal = round(float(placement[3]), 6)
                                top_placement.pendant = round(float(placement[4]), 6)
                                leaf_children = node_map[int(top_placement.inode)]
                                if len(leaf_children) > 1:
                                    # Reference tree with clade excluded
                                    parent = tmp_tree.get_common_ancestor(leaf_children)
                                    tip_distances = parent_to_tip_distances(parent, leaf_children)
                                    top_placement.mean_tip = round(float(sum(tip_distances)/len(tip_distances)), 6)
                    elif name == 'n':
                        top_placement.name = query_seq_name_map[int(info.pop())]
                    else:
                        logging.error("Unexpected variable in pquery keys: '" + name + "'\n")
                        sys.exit(33)

                if top_placement.lwr >= 0.5:  # The minimum likelihood weight ration a placement requires to be included
                    pqueries.append(top_placement)
                    taxonomic_placement_distances[rank].append(top_placement.total_distance())

            # Remove intermediate files from the analysis of this taxon
            intermediate_files += list(raxml_files.values())
            for old_file in intermediate_files:
                os.remove(old_file)
            # Clear collections
            taxonomy_filtered_query_seqs.clear()
            intermediate_files.clear()
            pruned_ref_fasta_dict.clear()
            seq_dict.clear()
            query_seq_name_map.clear()

            while acc > step_proportion:
                acc -= step_proportion
                sys.stdout.write('-')
                sys.stdout.flush()

        if len(taxonomic_placement_distances[rank]) == 0:
            logging.debug("No samples available for " + rank + ".\n")
        else:
            stats_string = "RANK: " + rank + "\n"
            stats_string += "\tSamples = " + str(len(taxonomic_placement_distances[rank])) + "\n"
            stats_string += "\tMedian = " + str(round(utilities.median(taxonomic_placement_distances[rank]), 4)) + "\n"
            stats_string += "\tMean = " + str(round(float(sum(taxonomic_placement_distances[rank])) /
                                                    len(taxonomic_placement_distances[rank]), 4)) + "\n"
            logging.debug(stats_string)
    sys.stdout.write("-]\n")
    return taxonomic_placement_distances, pqueries


def regress_rank_distance(args, ref_pkg: ReferencePackage, accession_lineage_map, ref_fasta_dict, training_ranks=None):
    """

    :param args:
    :param ref_pkg: A ReferencePackage instance
    :param accession_lineage_map:
    :param ref_fasta_dict:
    :param training_ranks:
    :return:
    """
    if not training_ranks:
        training_ranks = {"Class": 2, "Species": 6}
    # Read the taxonomic map; the final sequences used to build the tree are inferred from this
    leaf_taxa_map = dict()
    ref_taxa_map = tax_ids_file_to_leaves(ref_pkg.lineage_ids)
    for ref_seq in ref_taxa_map:
        leaf_taxa_map[ref_seq.number] = ref_seq.lineage
    # Find non-redundant set of diverse sequences to train
    rank_training_seqs, dedup_fasta_dict = prepare_training_data(args.fasta_input, args.output_dir, args.executables,
                                                                 leaf_taxa_map, accession_lineage_map, training_ranks)
    # Perform the rank-wise clade exclusion analysis for estimating placement distances
    taxonomic_placement_distances, pqueries = train_placement_distances(rank_training_seqs, training_ranks,
                                                                        ref_fasta_dict, dedup_fasta_dict,
                                                                        ref_pkg, leaf_taxa_map,
                                                                        args.molecule, args.executables,
                                                                        args.num_threads)
    # Finish up
    pfit_array = complete_regression(taxonomic_placement_distances, training_ranks)
    if pfit_array:
        logging.info("Placement distance regression model complete.\n")
    else:
        logging.info("Unable to complete phylogenetic distance and rank correlation.\n")

    return pfit_array, taxonomic_placement_distances, pqueries


def main():
    args = get_options()

    if not os.path.isdir(args.output_dir):
        os.makedirs(args.output_dir)
    args.treesapp = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + os.sep
    args = utilities.find_executables(args)

    # Limit this to just Class, Family, and Species - other ranks are inferred through regression
    training_ranks = {"Class": 2, "Species": 6}

    sys.stdout.write("\n##\t\t\tEstimate " + args.name + " taxonomic rank placement distances\t\t\t##\n")
    prep_logging(args.output_dir + os.sep + "placement_trainer_log.txt", args.verbose)

    ref_pkg = ReferencePackage()
    ref_pkg.gather_package_files(args.name, args.pkg_path, "flat")
    ref_pkg.validate()
    # TODO: wrap this into a function along with dependency version logging
    logging.debug("ANALYSIS SPECIFICATIONS:\n" +
                  "\tQuery FASTA: " + args.fasta_input + "\n" +
                  "\tTaxonomy map: " + ref_pkg.lineage_ids + "\n" +
                  "\tReference tree: " + ref_pkg.tree + "\n" +
                  "\tReference FASTA: " + ref_pkg.msa + "\n" +
                  "\tLineage map: " + str(args.lineages) + "\n" +
                  "\tRanks tested: " + ','.join(training_ranks.keys()) + "\n")

    # Get the model to be used for phylogenetic placement
    marker_build_dict = parse_ref_build_params(args)
    for denominator in marker_build_dict:
        marker_build = marker_build_dict[denominator]
        if marker_build.cog == args.name and args.molecule == marker_build.molecule:
            ref_pkg.sub_model = marker_build_dict[denominator].model
            break
    if not ref_pkg.sub_model:
        logging.error("Unable to find the substitution model used for " + args.name + ".\n")
        sys.exit(33)

    if args.domain:
        hmm_purified_fasta = args.output_dir + args.name + "_hmm_purified.fasta"
        logging.info("Searching for domain sequences... ")
        hmm_domtbl_files = utilities.hmmsearch_input_references(args, args.fasta_input)
        logging.info("done.\n")
        hmm_matches = parse_domain_tables(args, hmm_domtbl_files)
        fasta_dict = format_read_fasta(args.fasta_input, args.molecule, args.output_dir)
        header_registry = register_headers(get_headers(args.fasta_input))
        marker_gene_dict = utilities.extract_hmm_matches(hmm_matches, fasta_dict, header_registry)
        write_new_fasta(marker_gene_dict, hmm_purified_fasta)
        summarize_fasta_sequences(hmm_purified_fasta)
        utilities.hmm_pile(hmm_matches)
        # Point all future operations to the HMM purified FASTA file as the original input
        args.fasta_input = hmm_purified_fasta

    # Get the lineage information for the training/query sequences
    if args.lineages:
        accession_lineage_map = read_accession_taxa_map(args.lineages)
    else:
        header_registry = register_headers(get_headers(args.fasta_input))
        fasta_record_objects = get_header_info(header_registry)
        query_accession_list, num_lineages_provided = build_entrez_queries(fasta_record_objects)
        entrez_records = get_multiple_lineages(query_accession_list, args.molecule)
        accession_lineage_map = entrez_records_to_accession_lineage_map(entrez_records)
        all_accessions = entrez_records_to_accessions(entrez_records, query_accession_list)
        fasta_record_objects, accession_lineage_map = verify_lineage_information(accession_lineage_map,
                                                                                 all_accessions,
                                                                                 fasta_record_objects,
                                                                                 num_lineages_provided)
        write_accession_lineage_map(args.output_dir + os.sep + "placement_trainer_accession_lineage_map.tsv",
                                    accession_lineage_map)

    # Read in the reference fasta file
    ref_fasta_dict = read_fasta_to_dict(ref_pkg.msa)

    placement_table_file = args.output_dir + os.sep + "placement_info.tsv"
    placement_summary_file = args.output_dir + os.sep + "placement_trainer_results.txt"
    taxonomic_placement_distances = dict()

    # Goal is to use the distances already calculated but re-print
    if os.path.isfile(placement_summary_file) and not args.overwrite:
        # Read the summary file and pull the phylogenetic distances for each rank
        taxonomic_placement_distances = read_placement_summary(placement_summary_file)
    if len(taxonomic_placement_distances) == 0:
        pfit_array, taxonomic_placement_distances, pqueries = regress_rank_distance(args,
                                                                                    ref_pkg,
                                                                                    accession_lineage_map,
                                                                                    ref_fasta_dict,
                                                                                    training_ranks)
        # Write the tab-delimited file with metadata included for each placement
        write_placement_table(pqueries, placement_table_file, args.name)
    else:
        pfit_array = complete_regression(taxonomic_placement_distances)
        if pfit_array:
            logging.info("Placement distance regression model complete.\n")
        else:
            logging.info("Unable to complete phylogenetic distance and rank correlation.\n")

    # Write the text file containing distances used in the regression analysis
    with open(placement_summary_file, 'w') as out_handler:
        trained_string = "Regression parameters = " + re.sub(' ', '', str(pfit_array)) + "\n"
        ranks = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]
        for rank in ranks:
            trained_string += "# " + rank + "\n"
            if rank in taxonomic_placement_distances:
                trained_string += str(sorted(taxonomic_placement_distances[rank], key=float)) + "\n"
            trained_string += "\n"
        out_handler.write(trained_string)


if __name__ == "__main__":
    main()
