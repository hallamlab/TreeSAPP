#!/usr/bin/env python3

import os
import sys
import logging
import re
import random

from tqdm import tqdm

from treesapp import file_parsers
from treesapp import utilities
from treesapp import wrapper
from treesapp import fasta
from treesapp.phylo_seq import PQuery
from treesapp.phylo_dist import cull_outliers, regress_ranks
from treesapp.taxonomic_hierarchy import TaxonomicHierarchy
from treesapp.refpkg import ReferencePackage
from treesapp.training_utils import rarefy_rank_distances, generate_pquery_data_for_trainer


def fail_training(msg) -> None:
    boilerplate = "Clade-exclusion analysis could not be performed for training the reference package models:\n"
    logging.error(boilerplate + msg + "\n")
    return


def write_placement_table(pqueries: dict, placement_table_file: str, marker: str) -> None:
    header = ["Marker", "Rank", "Lineage", "Query.Name", "Internal.Node", "Placement.LWR", "Tree.Likelihood",
              "Dist.Distal", "Dist.Pendant", "Dist.MeanTip", "Dist.Total"]
    placement_info_strs = list()
    for rank, taxa in pqueries.items():
        for taxon in taxa:
            for pquery in taxa[taxon]:  # type: PQuery
                if pquery:
                    placement_info_strs.append("\t".join(
                        [marker, str(pquery.rank), str(pquery.lineage),
                         str(pquery.ref_name), str(pquery.consensus_placement.edge_num),
                         str(pquery.consensus_placement.like_weight_ratio), str(pquery.consensus_placement.likelihood),
                         str(pquery.consensus_placement.distal_length), str(pquery.consensus_placement.pendant_length),
                         str(pquery.consensus_placement.mean_tip_length),
                         str(pquery.consensus_placement.total_distance())])
                    )

    with open(placement_table_file, 'w') as file_handler:
        file_handler.write('#' + "\t".join(header) + "\n")
        file_handler.write("\n".join(placement_info_strs) + "\n")
    return


def flatten_pquery_dict(pqueries: dict, refpkg_prefix: str) -> dict:
    """
    Takes a dictionary storing PQuery instances in values and extracts them to a new
    dictionary where the key is refpkg_prefix and the value is a list of PQuery instances.

    :param pqueries: A dictionary mapping keys (function is agnostic though some examples are taxonomic rank or lineage)
    to an iterable - either a list containing ITolJPlace instances (or a subclass) or a dictionary whose values are
    a list of PQuery instances.
    :param refpkg_prefix: A ReferencePackage.prefix to which the PQuery instances were classified by
    :return: A dictionary indexed by refpkg_prefix mapped to a list of PQuery instances
    """
    refpkg_pqueries = {refpkg_prefix: []}
    for _, taxa in pqueries.items():  # type: (str, dict)
        for taxon in taxa:  # type: str
            try:
                for pquery in taxa[taxon]:  # type: PQuery
                    refpkg_pqueries[refpkg_prefix].append(pquery)
            except TypeError:
                if isinstance(taxon, PQuery):
                    refpkg_pqueries[refpkg_prefix].append(taxon)
                else:
                    logging.error("An instance of type PQuery was expected, found '{}' instead.\n".format(type(taxon)))
                    raise TypeError
    return refpkg_pqueries


def read_placement_summary(placement_summary_file: str) -> dict:
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
                    dist_strings = re.sub(r'[\[\]]', '', line).split(", ")
                    try:
                        dists = [float(dist) for dist in dist_strings]
                    except ValueError:
                        logging.error("Looks like treesapp train did not complete successfully.\n"
                                      "Please re-run with the flag '--overwrite'.\n")
                        sys.exit(5)
                    if len(dists) > 1:
                        taxonomic_placement_distances[rank] = dists
            line = place_summary.readline()
    return taxonomic_placement_distances


def complete_regression(taxonomic_placement_distances, taxonomic_ranks=None) -> (float, float):
    """
    Wrapper for performing outlier removal, normalization via rarefaction, and regression

    :param taxonomic_placement_distances:
    :param taxonomic_ranks: A dictionary mapping rank names to rank depth values where domain is 0, phylum is 1, etc.
    :return: Tuple of floats representing the slope and intercept estimated from linear regression
    """
    if not taxonomic_placement_distances:
        return []

    if not taxonomic_ranks:
        taxonomic_ranks = {"phylum": 2, "class": 3, "order": 4, "family": 5, "genus": 6, "species": 7, "strain": 8}

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

    return regress_ranks(rarefied_pds, taxonomic_ranks)


def prepare_training_data(test_seqs: fasta.FASTA, output_dir: str, executables: dict, leaf_taxa_map: dict,
                          t_hierarchy: TaxonomicHierarchy, accession_lineage_map: dict, taxonomic_ranks=None) -> dict:
    """
    Function for creating a non-redundant inventory of sequences to be used for training the rank-placement distance
    linear model. Removes sequences that share an identical accession, are more than 97% similar and limits the
    number of taxonomically-identical sequences to 30.

    :param test_seqs: A FASTA object. All headers in FASTA.header_registry must have their accession attribute filled
    :param output_dir: Path to write intermediate output files (such as UCLUST outputs)
    :param executables: A dictionary mapping software to a path of their respective executable
    :param t_hierarchy: A populated TaxonomicHierarchy instance for the reference package
    :param leaf_taxa_map: A dictionary mapping TreeSAPP numeric identifiers of reference sequences to taxonomic lineages
    :param accession_lineage_map: A dictionary mapping header accession IDs to full NCBI taxonomic lineages
    :param taxonomic_ranks: A set of rank names (e.g. Phylum) the NCBI taxonomic hierarchy
     to that could be mapped to rank depth values where Kingdom is 0, Phylum is 1, etc.
    :return: A dictionary storing the sequence accession names being used to test each taxon within each rank,
     so the structure is {'rank': {'taxon': [accession_1, accession_2]}}
    """
    similarity = 0.99  # The proportional similarity to cluster the training sequences
    max_reps = 30  # The maximum number of representative sequences from a specific taxon for training
    warning_threshold = 10  # Emit a warning if the number of taxa representing a rank drops below this proportion
    min_refpkg_size = 10  # Minimum number of sequences in a clade-excluded reference package for it to be used
    min_refpkg_proportion = int(len(leaf_taxa_map)*0.5)
    clustering_prefix = output_dir + os.sep + "cluster" + str(similarity)
    clustering_input = output_dir + os.sep + "clustering_input.fasta"
    lin_sep = t_hierarchy.lin_sep

    rank_training_seqs = dict()
    optimal_assignment_missing = set()
    too_short = list()
    taxon_training_queries = list()
    taxon_contributions = dict()
    unrelated_queries = list()
    related_queries = list()
    optimal_lineages_present = 0
    represented_taxa = 0
    test_seq_found = 0

    if not taxonomic_ranks:
        taxonomic_ranks = set([rank for rank, depth in t_hierarchy.accepted_ranks_depths.items() if depth > 1])

    # Cluster the training sequences to mitigate harmful redundancy
    # Remove fasta records with duplicate accessions
    test_seqs.dedup_by_accession()
    # Remove fasta records with duplicate sequences
    test_seqs.dedup_by_sequences()
    test_seqs.change_dict_keys("accession")

    # Remove sequences that are not related at the rank of Domain
    ref_domains = t_hierarchy.rank_representatives("domain", True)
    for seq_name in sorted(accession_lineage_map):  # type: str
        query_domain = t_hierarchy.get_rank_from_lineage(accession_lineage_map[seq_name].split(lin_sep)[-1], "domain")
        if not query_domain:
            continue
        elif query_domain.prefix_taxon() not in ref_domains:
            unrelated_queries.append(seq_name)
        else:
            related_queries.append(seq_name)
    if not related_queries:
        logging.error("No sequences were retained after filtering reference sequences by domains '%s'\n" %
                      str(', '.join(ref_domains)))
        sys.exit(5)

    related_queries_in_fasta = list(set(test_seqs.fasta_dict.keys()).intersection(set(related_queries)))
    test_seqs.keep_only(related_queries_in_fasta)
    test_seqs.change_dict_keys("accession")
    [accession_lineage_map.pop(seq_name) for seq_name in unrelated_queries]

    # Calculate the number of sequences that cannot be used in clade exclusion analysis due to no coverage in the input
    test_taxa_summary = []
    for rank in taxonomic_ranks:
        trimmed_ref_lineages = t_hierarchy.trim_lineages_to_rank(leaf_taxa_map, rank)
        if not trimmed_ref_lineages:
            logging.warning("No reference sequences are resolved to the rank '{}' in reference package.\n".format(rank))
            continue
        test_taxa_summary.append("Sequences available for training %s-level placement distances:" % rank)
        unique_ref_lineages = sorted(set(trimmed_ref_lineages.values()))

        # Remove all sequences belonging to a taxonomic rank from tree and reference alignment
        for taxonomy in unique_ref_lineages:
            unrelated_ref_names = ReferencePackage.get_unrelated_taxa(leaf_taxa_map, taxonomy)
            optimal_lca_taxonomy = lin_sep.join(taxonomy.split(lin_sep)[:-1])
            unrelated_refs = set([lin_sep.join(trimmed_ref_lineages[ref_name].split(lin_sep)[:-1])
                                  for ref_name in unrelated_ref_names.intersection(set(trimmed_ref_lineages))])
            # Ensure a representative of the optimal taxonomic assignment is present in the truncated refpkg
            if optimal_lca_taxonomy not in unrelated_refs:
                optimal_assignment_missing.add(optimal_lca_taxonomy)
            # Ensure there are a sufficient number of reference sequences in the truncated refpkg
            elif len(unrelated_ref_names) < max([min_refpkg_size, min_refpkg_proportion]):
                too_short.append(taxonomy)
            else:
                for seq_name in sorted(accession_lineage_map, key=lambda x: accession_lineage_map[x]):
                    # Not all keys in accession_lineage_map are in fasta_dict (duplicate sequences were removed)
                    if re.search(taxonomy, accession_lineage_map[seq_name]) and seq_name in test_seqs.fasta_dict:
                        taxon_training_queries.append(seq_name)
                        test_seq_found = 1
                represented_taxa += test_seq_found
                optimal_lineages_present += 1
            test_seq_found = 0

            test_taxa_summary.append("\t{}\t{}".format(len(taxon_training_queries), taxonomy))
            taxon_training_queries.clear()

        taxonomic_coverage = round(float(represented_taxa*100/len(unique_ref_lineages)), 2)

        if taxonomic_coverage < warning_threshold:
            logging.warning("Only {0}% of unique {1}-level taxa can be used represent"
                            " {1} phylogenetic placements.\n".format(taxonomic_coverage, rank))

        test_taxa_summary.append("%d/%d unique %s-level taxa have training sequences.\n" % (represented_taxa,
                                                                                            len(unique_ref_lineages),
                                                                                            rank))
        logging.debug("%.1f%% of optimal %s lineages are present in the pruned trees.\n" %
                      (round(float(optimal_lineages_present*100/len(unique_ref_lineages)), 1), rank))
        optimal_lineages_present = 0
        represented_taxa = 0

    logging.debug("Optimal placement target was not found in the pruned tree for following taxa:\n\t" +
                  "\n\t".join(optimal_assignment_missing) + "\n")
    
    logging.debug("Unable to generate placement data for the following taxa since the refpkg would be too small:\n\t" +
                  "\n\t".join(too_short) + "\n")

    logging.debug("\n".join(test_taxa_summary) + "\n")

    test_seqs.change_dict_keys("num")
    fasta.write_new_fasta(test_seqs.fasta_dict, clustering_input)
    wrapper.cluster_sequences(executables["mmseqs"], clustering_input, clustering_prefix, similarity)
    cluster_dict = file_parsers.create_mmseqs_clusters(clusters_tbl=clustering_prefix + "_cluster.tsv",
                                                       aln_tbl=clustering_prefix + "_cluster_aln.tsv")
    test_seqs.keep_only([cluster_dict[clust_id].representative for clust_id in cluster_dict.keys()])
    logging.debug("\t" + str(len(test_seqs.fasta_dict.keys())) + " sequence clusters\n")

    logging.info("Preparing deduplicated sequence set for training... ")
    test_seqs.change_dict_keys("accession")

    # Determine the set of reference sequences to use at each rank
    for rank in taxonomic_ranks:
        rank_training_seqs[rank] = dict()
        leaf_trimmed_taxa_map = t_hierarchy.trim_lineages_to_rank(leaf_taxa_map, rank)
        unique_taxonomic_lineages = sorted(set(leaf_trimmed_taxa_map.values()))

        # Remove all sequences belonging to a taxonomic rank from tree and reference alignment
        for taxonomy in unique_taxonomic_lineages:
            optimal_assignment = lin_sep.join(taxonomy.split(lin_sep)[:-1])
            if optimal_assignment not in optimal_assignment_missing and taxonomy not in too_short:
                for seq_name in sorted(accession_lineage_map):
                    # Not all keys in accession_lineage_map are in fasta_dict (duplicate sequences were removed)
                    seq_lineage = accession_lineage_map[seq_name]
                    if re.search(taxonomy, seq_lineage) and seq_name in test_seqs.fasta_dict:
                        try:
                            contrib = taxon_contributions[seq_lineage]
                        except KeyError:
                            contrib = 0
                            taxon_contributions[seq_lineage] = contrib
                        if contrib < int(0.3*max_reps):
                            taxon_training_queries.append(seq_name)
                            taxon_contributions[seq_lineage] += 1
                    if len(taxon_training_queries) == max_reps:
                        break
                if len(taxon_training_queries) > 0:
                    rank_training_seqs[rank][taxonomy] = list(taxon_training_queries)
                    taxon_training_queries.clear()
        taxon_contributions.clear()
    logging.info("done.\n")

    return rank_training_seqs


def clade_exclusion_phylo_placement(rank_training_seqs: dict,
                                    test_fasta: fasta.FASTA, ref_pkg: ReferencePackage,
                                    executables: dict, output_dir="./", raxml_threads=4, min_seqs=30) -> dict:
    """
    Function for iteratively performing leave-one-out analysis for every taxonomic lineage represented in the tree,
    yielding an estimate of placement distances corresponding to taxonomic ranks.

    :param rank_training_seqs: A dictionary storing the sequence names being used to test each taxon within each rank
     to rank depth values where Kingdom is 0, Phylum is 1, etc.
    :param test_fasta: Dictionary with headers as keys and sequences as values for deduplicated training sequences
    :param ref_pkg: A ReferencePackage instance
    :param executables: A dictionary mapping software to a path of their respective executable
    :param output_dir: Path to directory where all intermediate files should be written
    :param raxml_threads: Number of threads to be used by RAxML for parallel computation
    :param min_seqs: The minimum number of sequences required for training all taxonomic ranks of a reference package

    :return: Dictionary of ranks indexing a dictionary of taxa of that rank indexing a list of PQuery instances
    """
    pqueries = dict()

    if output_dir[-1] != os.sep:
        output_dir += os.sep

    logging.debug("Calculating the total number of queries to be used for training... ")
    num_training_queries = 0
    for rank in rank_training_seqs:
        num_rank_training_seqs = 0
        for taxonomy in rank_training_seqs[rank]:
            num_rank_training_seqs += len(rank_training_seqs[rank][taxonomy])
        if len(rank_training_seqs[rank]) == 0:
            logging.debug("No sequences available for estimating {}-level placement distances.\n".format(rank))
            continue
        else:
            logging.debug("{} sequences to train {}-level placement distances\n".format(num_rank_training_seqs, rank))
        num_training_queries += num_rank_training_seqs

    if num_training_queries < min_seqs:
        fail_training("Too few ({}) sequences for training placement distance model.\n".format(num_training_queries))
        return pqueries
    if num_training_queries < 50:
        logging.warning("Only {} sequences for training placement distance model.\n".format(num_training_queries))
    logging.debug("done.\n")

    logging.info("Estimating branch-length placement distances for taxonomic ranks\n")
    pbar = tqdm(total=num_training_queries, ncols=100)

    for rank in sorted(rank_training_seqs, reverse=True):
        pbar.set_description("Processing %s" % rank)
        pqueries[rank] = {}
        for taxonomy in sorted(rank_training_seqs[rank]):
            logging.debug("Testing placements for {}:\n".format(taxonomy))
            pqueries[rank][taxonomy] = generate_pquery_data_for_trainer(ref_pkg, taxonomy, test_fasta,
                                                                        rank_training_seqs[rank][taxonomy], rank,
                                                                        executables, output_dir, pbar, raxml_threads)

        if len(pqueries[rank]) == 0:
            logging.debug("No samples available for " + rank + ".\n")

    pbar.close()

    return pqueries


def evo_dists_from_pqueries(pqueries: dict, training_ranks=None) -> dict:
    """
    Pull the evolutionary distance values from each PQuery instance for each rank.

    :param pqueries: A dictionary containing PQuery instances, indexed by both the rank their lineage was excluded
    from the reference package and their taxonomy. Example:
    {"class": {"Methanosarcinales": [PQuery.inst, PQuery.inst]}}
    :param training_ranks: A dictionary mapping the name of a taxonomic rank to its depth in the hierarchy
    :return: A dictionary of taxonomic ranks as keys and a list of evolutionary distances for each PQuery as values
    """
    taxonomic_placement_distances = {}
    if not training_ranks:
        training_ranks = {"class": 3, "species": 7}

    # Populate dictionary of evolutionary distances indexed by rank and taxon
    for rank in training_ranks:  # type: str
        if rank not in pqueries:
            fail_training(msg="No phylogenetic placements were performed for taxonomic rank '{}'.\n"
                              "Either find more reference sequences or increase percent similarity to "
                              "increase the chances of having multiple representatives of the same rank".format(rank))
            return {}

        taxonomic_placement_distances[rank] = []
        if rank in pqueries and len(pqueries[rank]) > 0:
            for taxon in pqueries[rank]:
                taxonomic_placement_distances[rank] += [pquery.consensus_placement.total_distance() for pquery
                                                        in pqueries[rank][taxon]]

        if len(taxonomic_placement_distances[rank]) == 0:
            fail_training(msg="No phylogenetic placements succeeded for taxonomic rank '{}'.\n".format(rank))
            return {}

        median_dist = round(utilities.median(taxonomic_placement_distances[rank]), 4)
        mean_dist = round(utilities.mean(taxonomic_placement_distances[rank]), 4)

        logging.debug("RANK: {}\n"
                      "\tSamples = {}\n"
                      "\tMedian = {}\n"
                      "\tMean = {}\n"
                      "".format(rank, len(taxonomic_placement_distances[rank]), median_dist, mean_dist))
    return taxonomic_placement_distances


def reduce_examples(candidate_seqs: dict, max_examples: int) -> dict:
    """
    When building reference packages for large, diverse protein families users may encounter many thousands of sequences
    that could be used for training the various models TreeSAPP uses to improve classifications.

    This function is used to equitably sample sequence names (i.e. examples for training) across
    taxa and taxonomic ranks.

    :param candidate_seqs: A dictionary containing the names of sequences to be used for training.
     The taxonomic rank indexes a sub-dictionary with taxa names as keys and a list of sequence names for values.
    :param max_examples: The maximum number of sequence names (to sample from candidate_seqs)
    :return: A dictionary following the same format as candidate_seqs but containing the number of examples desired
    """
    random.seed(987)
    rank_training_seqs = {}
    r_i = 0  # The index for iterating over the ranks
    n = 0  # The number of candidate sequences added to rank_training_seqs

    # Evenly sample the sequences from the rank_training_seqs dict across ranks and taxa
    ranks = list(candidate_seqs.keys())
    while n < max_examples:
        # Loop back to the beginning if the rank index has exceeded the number of ranks
        if r_i >= len(ranks):
            r_i = 0
        # Quite iterating if there are no more ranks
        if len(ranks) == 0:
            break
        rank = ranks[r_i]  # type: str
        if rank not in rank_training_seqs:
            rank_training_seqs[rank] = {}
        if len(candidate_seqs[rank]) == 0:
            candidate_seqs.pop(rank)
            ranks = list(candidate_seqs.keys())
            continue
        taxa_seqs = list(candidate_seqs[rank].keys())
        try:
            taxon = taxa_seqs[random.randint(a=0, b=len(taxa_seqs)-1)]  # type: str
        except ValueError:
            taxon = taxa_seqs[0]
        try:
            rank_training_seqs[rank][taxon].append(candidate_seqs[rank][taxon].pop())
        except KeyError:
            rank_training_seqs[rank][taxon] = [candidate_seqs[rank][taxon].pop()]

        # Clean up the empty lists and dictionaries following popping
        if len(candidate_seqs[rank][taxon]) == 0:
            candidate_seqs[rank].pop(taxon)
        if len(candidate_seqs[rank]) == 0:
            candidate_seqs.pop(rank)
            ranks = list(candidate_seqs.keys())
        r_i += 1
        n += 1

    return rank_training_seqs


def gen_cladex_data(fasta_input: str, executables: dict, ref_pkg: ReferencePackage,
                    accession_lineage_map: dict, output_dir: str,
                    max_examples=0, num_threads=2) -> dict:
    """
    Generate pquery instances resulting from phylogenetic placement using clade-excluded reference packages.
    These instances represent phylogenetic placements onto phylogenies lacking close relatives,
    as would be expected when classifying sequences from an environmental metagenome.

    :param fasta_input: Path to a FASTA file containing query sequences to be used during training
    :param executables: A dictionary containing executable names mapped to absolute paths of the executables
    :param ref_pkg: A ReferencePackage instance
    :param accession_lineage_map: Path to a file mapping query sequence accessions to their taxonomic lineage
    :param output_dir: Path the a directory to write the temporary files
    :param max_examples: The maximum number of examples to use for training the reference package's models
    :param num_threads: The number of threads to be used by the various dependencies during phylogenetic placement
    :return:
    """
    # Read the taxonomic map; the final sequences used to build the tree are inferred from this
    leaf_taxa_map = {}
    ref_pkg.load_taxonomic_hierarchy()
    for ref_seq in ref_pkg.generate_tree_leaf_references_from_refpkg():
        leaf_taxa_map[ref_seq.number] = ref_seq.lineage

    # Fix the taxonomic lineages of the query sequences in accession_lineage_map
    for seq_id, lineage in accession_lineage_map.items():
        taxon_name = lineage.split(ref_pkg.taxa_trie.lin_sep)[-1]
        if ref_pkg.taxa_trie.emit(taxon_name):
            accession_lineage_map[seq_id] = ref_pkg.taxa_trie.emit(taxon_name, with_prefix=True)
        else:
            pass

    # Load the query FASTA and unalign the sequences, in case the fasta is a MSA
    test_seqs = fasta.FASTA(fasta_input)
    test_seqs.load_fasta()
    test_seqs.unalign()
    test_seqs.add_accession_to_headers(ref_pkg.prefix)
    # Find non-redundant set of diverse sequences to train for all taxonomic ranks
    rank_training_seqs = prepare_training_data(test_seqs, output_dir, executables, leaf_taxa_map,
                                               ref_pkg.taxa_trie, accession_lineage_map)
    if max_examples > 0:
        rank_training_seqs = reduce_examples(rank_training_seqs, max_examples)

    if len(rank_training_seqs) == 0:
        return {}

    # Perform the rank-wise clade exclusion analysis for estimating placement distances
    pqueries = clade_exclusion_phylo_placement(rank_training_seqs, test_seqs, ref_pkg,
                                               executables, output_dir, num_threads)

    return pqueries
