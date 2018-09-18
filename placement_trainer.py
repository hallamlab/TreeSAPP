#!/usr/bin/env python3

import os
import sys
import argparse
import logging
import re
from ete3 import Tree
import numpy as np

from fasta import read_fasta_to_dict, write_new_fasta, deduplicate_fasta_sequences, trim_multiple_alignment
from file_parsers import tax_ids_file_to_leaves
from utilities import reformat_fasta_to_phy, write_phy_file, median, clean_lineage_string
from entrez_utils import read_accession_taxa_map, get_multiple_lineages, build_entrez_queries, \
    write_accession_lineage_map, verify_lineage_information
from phylo_dist import trim_lineages_to_rank, cull_outliers, parent_to_tip_distances, regress_ranks
from external_command_interface import launch_write_command, setup_progress_bar
from jplace_utils import jplace_parser
from treesapp import run_papara
from classy import prep_logging, register_headers, get_header_info, get_headers
from entish import map_internal_nodes_leaves

__author__ = 'Connor Morgan-Lang'


def get_options():
    parser = argparse.ArgumentParser(description="Workflow for estimating calibrating the edge distances corresponding"
                                                 " to taxonomic ranks by iterative leave-one-out validation")
    parser.add_argument("-f", "--fasta", required=True,
                        help='The raw, unclusterd and unfiltered FASTA file used to generated the reference package.')
    parser.add_argument("-m", "--taxa_map", required=True,
                        help="The tax_ids file generated from the marker's reference package.")
    parser.add_argument("-t", "--tree", required=True,
                        help="The Newick formatted tree generated from the marker's reference package.")
    parser.add_argument("-r", "--ref_seqs", required=True,
                        help="A FASTA file containing the aligned reference sequences.")
    parser.add_argument("-l", "--lineages",
                        help="The accession lineage map downloaded during reference package generation.",
                        required=False)
    parser.add_argument("-T", "--threads", required=False, default=4, type=int,
                        help="The number of threads to be used by RAxML.")
    parser.add_argument("-o", "--output_dir", required=False, default='.',
                        help="Path to directory for writing outputs.")
    args = parser.parse_args()
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
        return sum([self.pendant, self.mean_tip, self.distal])

    def summarize_placement(self):
        summary_string = "Placement of " + self.lineage + " at rank " + self.rank + ":\n" + \
                         "Distances:" + \
                         "\n\tDistal = " + str(self.distal) +\
                         "\n\tPendant = " + str(self.pendant) +\
                         "\n\tTip = " + str(self.mean_tip) +\
                         "\nLikelihood = " + str(self.likelihood) +\
                         "\nL.W.R. = " + str(self.lwr) +\
                         "\n"
        return summary_string


def write_placement_table(pqueries, placement_table_file, marker):

    placement_info_strs = list()
    for pquery in pqueries:
        if pquery:
            placement_info_strs.append("\t".join(
                [marker, str(pquery.rank), str(pquery.lineage), str(pquery.name), str(pquery.inode),
                 str(pquery.lwr), str(pquery.likelihood),
                 str(pquery.distal), str(pquery.pendant), str(pquery.mean_tip), str(pquery.total_distance())])
            )

    with open(placement_table_file, 'w') as file_handler:
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


def train_placement_distances(fasta_dict: dict, ref_fasta_dict: dict, ref_tree_file: str, tax_ids_file: str,
                              accession_lineage_map: dict, molecule: str, raxml_threads=4):
    """
    Function for iteratively performing leave-one-out analysis for every taxonomic lineage represented in the tree,
    yielding an estimate of placement distances corresponding to taxonomic ranks.

    :param fasta_dict: A dictionary with headers as keys and sequences as values for all potential reference sequences
    :param ref_fasta_dict: A dictionary with headers as keys and sequences as values containing only reference sequences
    :param ref_tree_file: A Newick-formatted phylogenetic tree with branch length distances (no internal nodes)
    :param tax_ids_file: A tabular file created by create_treesapp_ref_data.py
    :param accession_lineage_map: A dictionary mapping NCBI accession IDs to full NCBI taxonomic lineages
    :param molecule: Molecule type [prot | dna | rrna]
    :param raxml_threads: Number of threads to be used by RAxML for parallel computation

    :return:
    """

    logging.info("\nEstimating branch-length placement distances for taxonomic ranks. Progress:\n")

    leaf_taxa_map = dict()
    taxonomic_placement_distances = dict()
    taxonomy_filtered_query_seqs = dict()
    dict_for_phy = dict()
    rank_distance_ranges = dict()
    pqueries = list()
    # Limit this to just Class, Family, and Species - other ranks are inferred through regression
    taxonomic_ranks = {"Class": 2, "Species": 6}

    temp_tree_file = "tmp_tree.txt"
    temp_ref_phylip_file = "taxonomy_filtered_ref_seqs.phy"
    temp_query_fasta_file = "queries.fasta"
    query_multiple_alignment = "papara_queries_aligned.phy"

    # Read the tree as ete3 Tree instance
    ref_tree = Tree(ref_tree_file)
    # Read the taxonomic map; the final sequences used to build the tree are inferred from this
    ref_taxa_map = tax_ids_file_to_leaves(tax_ids_file)

    bmge_file = "sub_binaries" + os.sep + "BMGE.jar"
    if not os.path.exists(bmge_file):
        raise FileNotFoundError("Cannot find " + bmge_file)

    for ref_seq in ref_taxa_map:
        leaf_taxa_map[ref_seq.number] = ref_seq.lineage

    # Remove duplicate sequences to prevent biasing the distance estimates
    nr_fasta_dict = deduplicate_fasta_sequences(fasta_dict)
    fasta_dict.clear()
    for seq_name in nr_fasta_dict.keys():
        fasta_dict[seq_name.split(" ")[0]] = nr_fasta_dict[seq_name]
    nr_fasta_dict.clear()

    setup_progress_bar(len(taxonomic_ranks))
    # For each rank from Class to Species (Kingdom & Phylum-level classifications to be inferred by LCA):
    for rank in taxonomic_ranks:
        taxonomic_placement_distances[rank] = list()
        leaf_trimmed_taxa_map = trim_lineages_to_rank(leaf_taxa_map, rank)
        unique_taxonomic_lineages = sorted(set(leaf_trimmed_taxa_map.values()))

        # Add the lineages to the Tree instance
        for leaf in ref_tree:
            leaf.add_features(lineage=leaf_trimmed_taxa_map.get(leaf.name, "none"))

        # Remove all sequences belonging to a taxonomic rank from tree and reference alignment
        for taxonomy in unique_taxonomic_lineages:
            logging.debug("Testing placements for " + taxonomy + ":\n")
            # Clear collections
            taxonomy_filtered_query_seqs.clear()
            dict_for_phy.clear()
            leaves_excluded = 0

            optimal_lca_taxonomy = "; ".join(taxonomy.split("; ")[:-1])
            if optimal_lca_taxonomy not in ["; ".join(tl.split("; ")[:-1]) for tl in unique_taxonomic_lineages if
                                            tl != taxonomy]:
                logging.debug("Optimal placement target '" + optimal_lca_taxonomy + "' not found in pruned tree.\n")
                continue

            # Write query FASTA containing sequences belonging to `taxonomy`
            for seq_name in sorted(accession_lineage_map):
                # Not all keys in accession_lineage_map are in fasta_dict (duplicate sequences were removed)
                if re.search(taxonomy, clean_lineage_string(accession_lineage_map[seq_name])) and seq_name in fasta_dict:
                    taxonomy_filtered_query_seqs[seq_name] = fasta_dict[seq_name]
            logging.debug("\t" + str(len(taxonomy_filtered_query_seqs.keys())) + " query sequences.\n")
            if len(taxonomy_filtered_query_seqs) == 0:
                continue
            write_new_fasta(taxonomy_filtered_query_seqs, fasta_name=temp_query_fasta_file)

            for key in ref_fasta_dict.keys():
                node = key.split('_')[0]
                # Node with truncated and/or unclassified lineages are not in `leaf_trimmed_taxa_map`
                if node in leaf_trimmed_taxa_map and not re.match(taxonomy, leaf_trimmed_taxa_map[node]):
                    dict_for_phy[node] = ref_fasta_dict[key]
                else:
                    leaves_excluded += 1
            logging.debug("\t" + str(leaves_excluded) + " sequences pruned from tree.\n")

            # Write the reference MSA with sequences of `taxonomy` removed
            phy_dict = reformat_fasta_to_phy(dict_for_phy)
            write_phy_file(temp_ref_phylip_file, phy_dict)

            # Copy the tree since we are removing leaves of `taxonomy` and don't want this to be permanent
            tmp_tree = ref_tree.copy(method="deepcopy")
            tmp_tree.prune(dict_for_phy.keys())  # iteratively detaching the monophyletic clades generates a bad tree
            logging.debug("\t" + str(len(tmp_tree.get_leaves())) + " leaves in pruned tree.\n")
            # Write the new reference tree with sequences from `taxonomy` removed
            tmp_tree.write(outfile=temp_tree_file, format=5)

            # Run PaPaRa, BMGE and RAxML to map sequences from the taxonomic rank onto the tree
            papara_stdout = run_papara("papara",
                                       temp_tree_file, temp_ref_phylip_file, temp_query_fasta_file,
                                       "prot")
            os.rename("papara_alignment.default", query_multiple_alignment)
            logging.debug(str(papara_stdout) + "\n")

            query_filtered_multiple_alignment = trim_multiple_alignment(bmge_file, query_multiple_alignment,
                                                                        molecule, "BMGE")
            query_name = re.sub(' ', '_', taxonomy.split("; ")[-1])
            raxml_command = ["raxmlHPC",
                             "-m", "PROTGAMMALG",
                             "-p", str(12345),
                             '-T', str(raxml_threads),
                             '-s', query_filtered_multiple_alignment,
                             '-t', temp_tree_file,
                             '-G', str(0.2),
                             '-f', 'v',
                             '-n', query_name,
                             '>', 'RAxML.txt']
            logging.debug("RAxML placement command:\n" + ' '.join(raxml_command) + "\n")
            launch_write_command(raxml_command)
            # Parse the JPlace file to pull distal_length+pendant_length for each placement
            jplace_file = "RAxML_portableTree." + query_name + ".jplace"
            jplace_data = jplace_parser(jplace_file)
            placement_tree = jplace_data.tree
            node_map = map_internal_nodes_leaves(placement_tree)
            for pquery in jplace_data.placements:
                top_lwr = 0.5
                distance = 100
                top_placement = None
                seq_name = ''
                for name, info in pquery.items():
                    if name == 'p':
                        for placement in info:
                            # Only record the best placement's distance
                            lwr = float(placement[2])
                            if lwr > top_lwr:
                                top_lwr = lwr
                                top_placement = PQuery(taxonomy, rank)
                                top_placement.inode = placement[0]
                                top_placement.likelihood = placement[1]
                                top_placement.lwr = lwr
                                top_placement.distal = float(placement[3])
                                top_placement.pendant = float(placement[4])
                                leaf_children = node_map[int(top_placement.inode)]
                                if len(leaf_children) > 1:
                                    # Reference tree with clade excluded
                                    parent = tmp_tree.get_common_ancestor(leaf_children)
                                    tip_distances = parent_to_tip_distances(parent, leaf_children)
                                    top_placement.mean_tip = float(sum(tip_distances)/len(tip_distances))
                                distance = top_placement.total_distance()
                    elif name == 'n':
                        seq_name = info.pop()
                    else:
                        logging.error("Unexpected variable in pquery keys: '" + name + "'\n")
                        sys.exit(33)

                    if top_placement:
                        top_placement.name = seq_name
                        pqueries.append(top_placement)
                if distance < 100:
                    taxonomic_placement_distances[rank].append(distance)
            os.system("rm papara_* RAxML*")
        if len(taxonomic_placement_distances[rank]) == 0:
            logging.debug("No samples available for " + rank + ".\n")
        else:
            stats_string = "RANK: " + rank + "\n"
            stats_string += "\tSamples = " + str(len(taxonomic_placement_distances[rank])) + "\n"
            stats_string += "\tMedian = " + str(round(median(taxonomic_placement_distances[rank]), 4)) + "\n"
            stats_string += "\tMean = " + str(round(float(sum(taxonomic_placement_distances[rank])) /
                                                    len(taxonomic_placement_distances[rank]), 4)) + "\n"
            logging.debug(stats_string)
            rank_distance_ranges[rank] = cull_outliers(list(taxonomic_placement_distances[rank]))
        sys.stdout.write('-')
        sys.stdout.flush()
    sys.stdout.write("-]\n")
    os.system("rm taxonomy_filtered_ref_seqs.phy queries.fasta tmp_tree.txt")

    # Rarefy the placement distances to the rank with the fewest samples
    rank_distance_ranges = rarefy_rank_distances(rank_distance_ranges)
    pfit_array = regress_ranks(rank_distance_ranges, taxonomic_ranks)

    return pfit_array, taxonomic_placement_distances, pqueries


def main():
    args = get_options()

    if not os.path.isdir(args.output_dir):
        os.makedirs(args.output_dir)

    seq_m = os.path.basename(args.ref_seqs).split('.')[0]
    tree_m = re.sub("_tree", '', os.path.basename(args.tree)).split('.')[0]
    tax_m = re.sub("tax_ids_", '', os.path.basename(args.taxa_map)).split('.')[0]

    if len({seq_m, tree_m, tax_m}) != 1:
        logging.error("Marker gene names are inconsistent between reference package inputs.\n")
        sys.exit(33)

    molecule = "prot"
    sys.stdout.write("\n##\t\t\tEstimate taxonomic rank placement distances\t\t\t##\n")
    prep_logging("placement_trainer_log.txt", False)
    if args.lineages:
        accession_lineage_map = read_accession_taxa_map(args.lineages)
    else:
        header_registry = register_headers(get_headers(args.fasta))
        fasta_record_objects = get_header_info(header_registry)
        query_accession_list, num_lineages_provided = build_entrez_queries(fasta_record_objects)
        accession_lineage_map, all_accessions = get_multiple_lineages(query_accession_list, molecule)
        fasta_record_objects, accession_lineage_map = verify_lineage_information(accession_lineage_map,
                                                                                 all_accessions,
                                                                                 fasta_record_objects,
                                                                                 num_lineages_provided,
                                                                                 molecule)
        write_accession_lineage_map(args.output_dir + os.sep + "placement_trainer_accession_lineage_map.tsv",
                                    accession_lineage_map)

    # Read in the original fasta file
    fasta_dict = read_fasta_to_dict(args.fasta)
    ref_fasta_dict = read_fasta_to_dict(args.ref_seqs)

    pfit_array, taxonomic_placement_distances, pqueries = train_placement_distances(fasta_dict,
                                                                                    ref_fasta_dict,
                                                                                    args.tree,
                                                                                    args.taxa_map,
                                                                                    accession_lineage_map,
                                                                                    molecule,
                                                                                    args.threads)
    with open(args.output_dir + os.sep + "placement_trainer_results.txt", 'w') as out_handler:
        trained_string = "Polynomial params = " + str(pfit_array) + "\n"
        ranks = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]
        for rank in ranks:
            trained_string += "# " + rank + "\n"
            if rank in taxonomic_placement_distances:
                trained_string += str(sorted(taxonomic_placement_distances[rank], key=float)) + "\n"
            trained_string += "\n"
        out_handler.write(trained_string)

    placement_table_file = args.output_dir + os.sep + "placement_info.tsv"
    write_placement_table(pqueries, placement_table_file, seq_m)


if __name__ == "__main__":
    main()
