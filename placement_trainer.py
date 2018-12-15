#!/usr/bin/env python3

import os
import sys
import argparse
import logging
import re
from ete3 import Tree
import numpy as np
from glob import glob

from fasta import read_fasta_to_dict, write_new_fasta, deduplicate_fasta_sequences, trim_multiple_alignment
from file_parsers import tax_ids_file_to_leaves, read_stockholm_to_dict
from utilities import reformat_fasta_to_phy, write_phy_file, median, clean_lineage_string,\
    find_executables, cluster_sequences, profile_aligner, run_papara
from entrez_utils import read_accession_taxa_map, get_multiple_lineages, build_entrez_queries, \
    write_accession_lineage_map, verify_lineage_information
from phylo_dist import trim_lineages_to_rank, cull_outliers, parent_to_tip_distances, regress_ranks
from external_command_interface import launch_write_command, setup_progress_bar
from jplace_utils import jplace_parser
from classy import prep_logging, register_headers, get_header_info, get_headers, ReferencePackage
from entish import map_internal_nodes_leaves

__author__ = 'Connor Morgan-Lang'


def get_options():
    parser = argparse.ArgumentParser(description="Workflow for estimating calibrating the edge distances corresponding"
                                                 " to taxonomic ranks by iterative leave-one-out validation")
    parser.add_argument("-f", "--fasta_input", required=True,
                        help='The raw, unclustered and unfiltered FASTA file to train the reference package.')
    parser.add_argument("-n", "--name", required=True,
                        help="Prefix name of the reference package (i.e. McrA for McrA.fa, McrA.hmm, McrA_tree.txt)")
    parser.add_argument("-p", "--pkg_path", required=True,
                        help="The path to the TreeSAPP-formatted reference package.")
    parser.add_argument("-l", "--lineages",
                        help="The accession lineage map downloaded during reference package generation.",
                        required=False)
    parser.add_argument('-m', '--molecule', default='prot', choices=['prot', 'dna', 'rrna'],
                        help='the type of input sequences (prot = Protein [DEFAULT]; dna = Nucleotide )')
    parser.add_argument("-T", "--num_threads", required=False, default=4, type=int,
                        help="The number of threads to be used by RAxML.")
    parser.add_argument("-o", "--output_dir", required=False, default='.',
                        help="Path to directory for writing outputs.")
    parser.add_argument("-O", "--overwrite", default=False, action="store_true",
                        help="Force recalculation of placement distances for query sequences.")
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
    if not taxonomic_ranks:
        taxonomic_ranks = {"Phylum": 1, "Class": 2, "Order": 3, "Family": 4, "Genus": 5, "Species": 6, "Strain": 7}

    filtered_pds = dict()
    for rank in taxonomic_placement_distances:
        # print(rank, "raw", np.median(list(taxonomic_placement_distances[rank])))
        filtered_pds[rank] = cull_outliers(list(taxonomic_placement_distances[rank]))
        # print(rank, "filtered", np.median(list(filtered_pds[rank])))

    # Rarefy the placement distances to the rank with the fewest samples
    rarefied_pds = rarefy_rank_distances(filtered_pds)
    # for rank in rarefied_pds:
    #     print(rank, "rarefied", np.median(list(rarefied_pds[rank])))
    pfit_array = regress_ranks(rarefied_pds, taxonomic_ranks)
    return pfit_array


def prepare_training_data(fasta_input: str, output_dir: str, executables: dict,
                          leaf_taxa_map: dict, accession_lineage_map: dict, taxonomic_ranks):
    """
    Function for creating a non-redundant inventory of sequences to be used for training the rank-placement distance
    linear model. Removes sequences that share an identical accession, are more than 95% similar and limits the
    number of taxonomically-identical sequences to 30.

    :param fasta_input: 
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
    cluster_sequences(executables["usearch"], fasta_input, uclust_prefix, similarity)
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
                                 clean_lineage_string(accession_lineage_map[seq_name])) and seq_name in uclust_fasta_dict:
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
    dict_for_phy = dict()
    seq_dict = dict()
    pqueries = list()
    intermediate_files = list()
    aligner = "hmmalign"

    temp_tree_file = "tmp_tree.txt"
    temp_ref_phylip_file = "taxonomy_filtered_ref_seqs.phy"
    temp_query_fasta_file = "queries.fasta"
    query_multiple_alignment = "papara_queries_aligned.phy"

    # Read the tree as ete3 Tree instance
    ref_tree = Tree(ref_pkg.tree)

    bmge_file = executables["BMGE.jar"]
    if not os.path.exists(bmge_file):
        raise FileNotFoundError("Cannot find " + bmge_file)

    num_training_queries = 0
    for rank in rank_training_seqs:
        for taxonomy in rank_training_seqs[rank]:
            num_training_queries += len(rank_training_seqs[rank][taxonomy])

    if num_training_queries < 30:
        logging.error("Too few (" + str(num_training_queries) + ") sequences for training placement distance model.\n")
        sys.exit(19)
    if num_training_queries < 50:
        logging.warning("Only " + str(num_training_queries) + " sequences for training placement distance model.\n")
    step_proportion = setup_progress_bar(num_training_queries)
    acc = 0.0

    # For each rank from Class to Species (Kingdom & Phylum-level classifications to be inferred by LCA):
    for rank in rank_training_seqs:
        if rank not in taxonomic_ranks:
            logging.error("Rank '" + rank + "' not found in ranks being used for training.\n")
            sys.exit(33)
        taxonomic_placement_distances[rank] = list()
        leaf_trimmed_taxa_map = trim_lineages_to_rank(leaf_taxa_map, rank)
        
        # Add the lineages to the Tree instance
        for leaf in ref_tree:
            leaf.add_features(lineage=leaf_trimmed_taxa_map.get(leaf.name, "none"))

        # Remove all sequences belonging to a taxonomic rank from tree and reference alignment
        for taxonomy in rank_training_seqs[rank]:
            logging.debug("Testing placements for " + taxonomy + ":\n")
            leaves_excluded = 0

            # Write query FASTA containing sequences belonging to `taxonomy`
            for seq_name in rank_training_seqs[rank][taxonomy]:
                taxonomy_filtered_query_seqs[seq_name] = dedup_fasta_dict[seq_name]
            logging.debug("\t" + str(len(taxonomy_filtered_query_seqs.keys())) + " query sequences.\n")
            acc += len(taxonomy_filtered_query_seqs.keys())
            write_new_fasta(taxonomy_filtered_query_seqs, fasta_name=temp_query_fasta_file)

            for key in ref_fasta_dict.keys():
                node = key.split('_')[0]
                # Node with truncated and/or unclassified lineages are not in `leaf_trimmed_taxa_map`
                if node in leaf_trimmed_taxa_map and not re.match(taxonomy, leaf_trimmed_taxa_map[node]):
                    dict_for_phy[node] = ref_fasta_dict[key]
                else:
                    leaves_excluded += 1
            logging.debug("\t" + str(leaves_excluded) + " sequences pruned from tree.\n")

            # Copy the tree since we are removing leaves of `taxonomy` and don't want this to be permanent
            tmp_tree = ref_tree.copy(method="deepcopy")
            tmp_tree.prune(dict_for_phy.keys())  # iteratively detaching the monophyletic clades generates a bad tree
            logging.debug("\t" + str(len(tmp_tree.get_leaves())) + " leaves in pruned tree.\n")
            # Write the new reference tree with sequences from `taxonomy` removed
            tmp_tree.write(outfile=temp_tree_file, format=5)

            ##
            # Run hmmalign, BMGE and RAxML to map sequences from the taxonomic rank onto the tree
            ##
            if aligner == "papara":
                # Write the reference MSA with sequences of `taxonomy` removed
                phy_dict = reformat_fasta_to_phy(dict_for_phy)
                write_phy_file(temp_ref_phylip_file, phy_dict)
                aln_stdout = run_papara(executables["papara"],
                                        temp_tree_file, temp_ref_phylip_file, temp_query_fasta_file,
                                        "prot")
                intermediate_files.append(temp_ref_phylip_file)
                os.rename("papara_alignment.default", query_multiple_alignment)
            elif aligner == "hmmalign":
                sto_file = re.sub("\.phy$", ".sto", query_multiple_alignment)
                intermediate_files.append(sto_file)
                # Currently not supporting rRNA references (phylogenetic_rRNA)
                aln_stdout = profile_aligner(executables, ref_pkg.msa, ref_pkg.profile,
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
            else:
                logging.error("Unrecognised alignment tool '" + aligner + "'. Exiting now.\n")
                sys.exit(33)
            logging.debug(str(aln_stdout) + "\n")

            query_filtered_multiple_alignment = trim_multiple_alignment(bmge_file, query_multiple_alignment,
                                                                        molecule, "BMGE")
            query_name = re.sub(r"([ /])", '_', taxonomy.split("; ")[-1])
            raxml_command = [executables["raxmlHPC"],
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
            intermediate_files += glob("RAxML_*" + query_name + "*")
            intermediate_files.append("RAxML.txt")
            intermediate_files += glob(query_filtered_multiple_alignment + "*")
            intermediate_files += [query_multiple_alignment, temp_tree_file, temp_query_fasta_file]
            for old_file in intermediate_files:
                os.remove(old_file)
            # Clear collections
            taxonomy_filtered_query_seqs.clear()
            intermediate_files.clear()
            dict_for_phy.clear()
            seq_dict.clear()

            while acc > step_proportion:
                acc -= step_proportion
                sys.stdout.write('-')
                sys.stdout.flush()

        if len(taxonomic_placement_distances[rank]) == 0:
            logging.debug("No samples available for " + rank + ".\n")
        else:
            stats_string = "RANK: " + rank + "\n"
            stats_string += "\tSamples = " + str(len(taxonomic_placement_distances[rank])) + "\n"
            stats_string += "\tMedian = " + str(round(median(taxonomic_placement_distances[rank]), 4)) + "\n"
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

    logging.info("Placement distance regression model complete.\n")

    return pfit_array, taxonomic_placement_distances, pqueries


def main():
    args = get_options()

    if not os.path.isdir(args.output_dir):
        os.makedirs(args.output_dir)
    args.treesapp = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + os.sep
    args = find_executables(args)

    # Limit this to just Class, Family, and Species - other ranks are inferred through regression
    training_ranks = {"Class": 2, "Species": 6}

    sys.stdout.write("\n##\t\t\tEstimate taxonomic rank placement distances\t\t\t##\n")
    prep_logging(args.output_dir + os.sep + "placement_trainer_log.txt", False)

    ref_pkg = ReferencePackage()
    ref_pkg.gather_package_files(args.name, args.pkg_path)
    ref_pkg.validate()
    logging.debug("ANALYSIS SPECIFICATIONS:\n" +
                  "\tQuery FASTA: " + args.fasta_input + "\n" +
                  "\tTaxonomy map: " + ref_pkg.lineage_ids + "\n" +
                  "\tReference tree: " + ref_pkg.tree + "\n" +
                  "\tReference FASTA: " + ref_pkg.msa + "\n" +
                  "\tLineage map: " + str(args.lineages) + "\n" +
                  "\tRanks tested: " + ','.join(training_ranks.keys()) + "\n")
    if args.lineages:
        accession_lineage_map = read_accession_taxa_map(args.lineages)
    else:
        header_registry = register_headers(get_headers(args.fasta_input))
        fasta_record_objects = get_header_info(header_registry)
        query_accession_list, num_lineages_provided = build_entrez_queries(fasta_record_objects)
        accession_lineage_map, all_accessions = get_multiple_lineages(query_accession_list, args.molecule)
        fasta_record_objects, accession_lineage_map = verify_lineage_information(accession_lineage_map,
                                                                                 all_accessions,
                                                                                 fasta_record_objects,
                                                                                 num_lineages_provided,
                                                                                 args.molecule)
        write_accession_lineage_map(args.output_dir + os.sep + "placement_trainer_accession_lineage_map.tsv",
                                    accession_lineage_map)

    # Read in the reference fasta file
    ref_fasta_dict = read_fasta_to_dict(ref_pkg.msa)

    placement_table_file = args.output_dir + os.sep + "placement_info.tsv"
    placement_summary_file = args.output_dir + os.sep + "placement_trainer_results.txt"
    taxonomic_placement_distances = dict()
    pfit_array = list()

    # Goal is to use the distances already calculated but re-print
    if os.path.isfile(placement_summary_file) and not args.overwrite:
        # Read the summary file and pull the phylogenetic distances for each rank
        taxonomic_placement_distances = read_placement_summary(placement_summary_file)
        pfit_array = complete_regression(taxonomic_placement_distances)

    if len(taxonomic_placement_distances) == 0:
        pfit_array, taxonomic_placement_distances, pqueries = regress_rank_distance(args,
                                                                                    ref_pkg,
                                                                                    accession_lineage_map,
                                                                                    ref_fasta_dict,
                                                                                    training_ranks)

        # Write the tab-delimited file with metadata included for each placement
        write_placement_table(pqueries, placement_table_file, args.name)

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
