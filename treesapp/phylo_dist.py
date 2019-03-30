#!/usr/bin/env python3

import sys
import re
import logging
from ete3 import Tree
from .file_parsers import tax_ids_file_to_leaves
from .utilities import clean_lineage_string, median

import numpy as np
import scipy.optimize as so
np.random.seed(0)

__author__ = 'Connor Morgan-Lang'


def cull_outliers(data: list, dev=3):
    """
    Returns the Interquartile Range (IQR) of a list after filtering outliers
    based on log transformed data where outliers are farther than 1 std-dev from the median and
    an un-transformed distribution where outliers are farther than `dev` standard deviations from the median.

    :param data: A list of floats
    :param dev: Number of acceptable deviations from the median; beyond this, values are outliers and removed
    :return: A smaller list of floats
    """
    # Reject outliers from ln-transformed distribution
    ln_a = np.log10(1.0 * np.array(data))
    noo_a = np.power(10, ln_a[abs(ln_a - np.median(ln_a)) < 2 * np.std(ln_a)])

    # Reject outliers from untransformed distribution
    d = np.abs(noo_a - np.median(noo_a))
    mdev = np.median(d)
    s = d / mdev if mdev else 0
    noo_a = noo_a[s < dev]

    return list(noo_a)


def regress_ranks(rank_distance_ranges, taxonomic_ranks):

    all_distances = list()
    for rank in rank_distance_ranges:
        all_distances += list(rank_distance_ranges[rank])

    # Prep arrays for regression
    rank_depth_list = list()
    dist_list = list()
    depth_dist_dict = dict()
    for rank in rank_distance_ranges:
        depth = taxonomic_ranks[rank]
        rank_distances = rank_distance_ranges[rank]
        n_samples = len(rank_distances)
        if n_samples > 3:
            depth_dist_dict[depth] = rank_distances
            dist_list += rank_distances

    if len(depth_dist_dict.keys()) <= 1:
        logging.error("Only " + str(len(depth_dist_dict.keys())) + " ranks available for modelling.\n")
        sys.exit(33)

    dist_list.clear()
    for depth in sorted(depth_dist_dict, key=int):
        rank_depth_list += [depth] * len(depth_dist_dict[depth])
        dist_list += list(depth_dist_dict[depth])

    # For TreeSAPP predictions
    # opt_slope, intercept = [round(float(x), 4) for x in list(np.polyfit(dist_list, rank_depth_list, 1))]
    intercept = 7.0
    opt_slope = round(float(so.fmin(lambda m, x, y: ((m * x - y + intercept) ** 2).sum(),
                                    x0=-6.0,
                                    args=(dist_list, rank_depth_list))), 4)

    return opt_slope, intercept


def rank_recommender(phylo_dist: float, taxonomic_rank_pfit: list):
    """
    Determines the rank depth (for example Class == 2) a taxonomic lineage should be truncated to
     based on which rank distance range (in taxonomic_rank_intervals) phylo_dist falls into

    :param phylo_dist: Float representing the branch distance from the nearest node
    :param taxonomic_rank_pfit: Dictionary with rank keys (e.g. Class) and distance ranges (min, max) as values
    :return: int
    """

    if not taxonomic_rank_pfit:
        return 7

    # For a polynomial
    # polyreg = np.poly1d(taxonomic_rank_pfit)
    # depth = int(round(polyreg(phylo_dist)))

    slope, intercept = taxonomic_rank_pfit
    depth = int(round(phylo_dist*slope + intercept))

    return depth


def trim_lineages_to_rank(leaf_taxa_map: dict, rank: str):
    """
    Iterates a dictionary and trims the lineage string values to a specified rank.
     Also removes lineages that are unclassified at the desired rank or higher (closer to root)

    :param leaf_taxa_map: Maps indices to lineages
    :param rank: The taxonomic rank lineages need to be trimmed to
    :return: Dictionary of keys mapped to trimmed lineages
    """
    trimmed_lineage_map = dict()
    # ranks is offset by 1 (e.g. Kingdom is the first index and therefore should be 1) for the final trimming step
    ranks = {"Kingdom": 1, "Phylum": 2, "Class": 3, "Order": 4, "Family": 5, "Genus": 6, "Species": 7}
    unknowns_re = re.compile("unclassified|environmental sample", re.IGNORECASE)
    depth = ranks[rank]
    truncated = 0
    unclassified = 0
    for node_name in sorted(leaf_taxa_map):
        c_lineage = clean_lineage_string(leaf_taxa_map[node_name])
        c_lineage_s = c_lineage.split("; ")

        # Remove lineage from testing if the rank doesn't exist (unclassified at a high rank)
        if len(c_lineage_s) < depth:
            truncated += 1
            continue

        try:
            unknowns_re.search(c_lineage)
        except TypeError:
            logging.error("Unexpected type (" + str(type(c_lineage)) + ") for '" + str(c_lineage) + "'\n")
            sys.exit(33)

        if unknowns_re.search(c_lineage):
            i = 0
            while i < depth:
                try:
                    taxon = c_lineage_s[i]
                except IndexError:
                    logging.error(rank + " position (" + str(depth) + ") unavailable in " + c_lineage + " ")
                    break
                if unknowns_re.search(taxon):
                    i -= 1
                    break
                i += 1
            if i < depth:
                unclassified += 1
                continue
        trimmed_lineage_map[node_name] = "; ".join(c_lineage_s[:depth])

    logging.debug(str(truncated) + " lineages truncated before " + rank + " were removed during lineage trimming.\n" +
                  str(unclassified) + " lineages unclassified at or before " + rank + " also removed.\n")
    return trimmed_lineage_map


def prune_branches(tree, leaf_taxa_map: dict, rank="Genus"):
    """
    Function for removing leaves of unclassified and polyphyletic lineages

    :type tree: Tree()
    :param tree: An Environment for Tree Exploration (ETE) Tree object
    :param leaf_taxa_map: A dictionary mapping tree leaf number keys to NCBI lineage strings
    :param rank: A taxonomic rank to test for monophyly
    :return: pruned_nodes dict() of Tree() nodes
    """
    pruned_nodes = dict()
    if not isinstance(tree, Tree):
        logging.error("Tree is not ete tree object.\n")
        raise AssertionError()
    # Check to see if the two collections are comparable
    for leaf in tree:
        if leaf.name not in leaf_taxa_map.keys():
            logging.error(str(leaf.name) + " not found in leaf_taxa_map.\n")
            raise AssertionError("Leaves in tree and tax_ids file are disparate sets.\n")
    # Raw lineages are too specific to test for monophyly, so try at a deeper rank by trimming the lineages
    leaf_taxa_map = trim_lineages_to_rank(leaf_taxa_map, rank)

    # Add the lineages to the Tree instance
    for leaf in tree:
        leaf.add_features(lineage=leaf_taxa_map.get(leaf.name, "none"))

    # Print the tree for debugging:
    # print(tree.get_ascii(attributes=["name", "lineage"], show_internal=False))

    # Add all the monophyletic leaf node numbers to pruned_nodes
    unique_lineages = sorted(list(set(leaf_taxa_map.values())))
    for lineage in unique_lineages:
        pruned_nodes[lineage] = dict()
        acc = 1  # In case lineages are scattered (e.g. paralogs in the tree) these need to be indexed
        for node in tree.get_monophyletic(values=[lineage], target_attr="lineage"):
            leaf_descendents = node.get_leaves()
            if len(leaf_descendents) == 1:
                pass
            elif len(leaf_descendents) > 1:
                # This is an internal node
                pruned_nodes[lineage][acc] = leaf_descendents
                acc += 1
            else:
                raise AssertionError("Expected at least one leaf leading from node " + str(node.name))
    return pruned_nodes


def parent_to_tip_distances(parent: Tree, children: Tree, estimate=False):
    """
    Function utilizing ete3's tree object for calculating distances between a reference node (parent)
     and query nodes (children).
    The `estimate` flag will cause the parent's edge length to be included in the distance calculation.

    :param parent: A reference node Tree instance
    :param children: A list of query nodes, also Tree instances
    :param estimate: Boolean indicating whether these distances are to be used for estimating the edge length ranges
    :return: list() of all branch distances between the parent node and the tips
    """
    branch_distances = list()
    # Calculate distance between parent and all descendants
    for child_node in children:
        if isinstance(child_node, Tree):
            distal_length = parent.get_distance(child_node.name)
        elif isinstance(child_node, str):
            distal_length = parent.get_distance(child_node)
        elif isinstance(child_node, int):
            distal_length = parent.get_distance(str(child_node))
        else:
            logging.error("Cannot handle type '" + type(child_node) + "' for child.")
            raise AssertionError()
        if estimate:
            distal_length += parent.dist
        branch_distances.append(distal_length)
    return branch_distances


def bound_taxonomic_branch_distances(tree, leaf_taxa_map):
    # Seed the final dictionary to be used for bounding
    taxonomic_rank_distances = dict()
    lca_nodes = dict()
    taxonomic_rank_intervals = dict()
    ranks = {1: "Phylum", 2: "Class", 3: "Order", 4: "Family", 5: "Genus", 6: "Species"}
    for depth in sorted(ranks, key=int, reverse=True):
        rank = ranks[depth]
        lca_nodes[depth] = list()
        taxonomic_rank_distances[rank] = list()
        # Remove the leaves belonging to "Unclassified" lineages and/or polyphyletic clades
        monophyletic_clades = prune_branches(tree, leaf_taxa_map, rank)
        for lineage in monophyletic_clades:
            for para_group in monophyletic_clades[lineage]:
                edge_lengths = []
                clade = monophyletic_clades[lineage][para_group]
                if len(clade) == 1:
                    # Unable to calculate distance with a single leaf
                    continue
                lca = clade[0].get_common_ancestor(clade[1:])
                lca_nodes[depth].append(lca)
                if depth + 1 in lca_nodes:
                    # Ensure this parent node wasn't already used to estimate the previous rank's bounds
                    if lca not in lca_nodes[depth + 1]:
                        edge_lengths = parent_to_tip_distances(lca, clade, True)
                else:
                    edge_lengths += parent_to_tip_distances(lca, clade, True)
                taxonomic_rank_distances[rank] += edge_lengths
        if len(taxonomic_rank_distances[rank]) >= 5:
            noo_taxa_rank_dists = cull_outliers(taxonomic_rank_distances[rank])
            taxonomic_rank_intervals[rank] = (round(np.percentile(noo_taxa_rank_dists, q=2.5), 4),
                                              round(np.percentile(noo_taxa_rank_dists, q=97.5), 4))
        else:
            taxonomic_rank_intervals[rank] = ()
    return taxonomic_rank_intervals


def validate_rank_intervals(taxonomic_rank_intervals):
    """
    Placeholder, for now.
        This idea has been abandoned but a similar function will be required soon.
    :param taxonomic_rank_intervals:
    :return:
    """
    validated_intervals = dict()
    ranks = {1: "Phylum", 2: "Class", 3: "Order", 4: "Family", 5: "Genus", 6: "Species"}
    ci_ranges = [max_ci-min_ci for min_ci, max_ci in taxonomic_rank_intervals.values()]
    ci_maxes = [max_ci for min_ci, max_ci in taxonomic_rank_intervals.values()]
    p_min = p_max = 1
    for depth in sorted(ranks):
        rank = ranks[depth]
        if rank not in taxonomic_rank_intervals:
            if depth == 1:
                min_ci = max(ci_maxes)
                max_ci = min_ci + median(ci_ranges)
            else:
                max_ci = p_min
                min_ci = max_ci - median(ci_ranges)
        else:
            min_ci, max_ci = taxonomic_rank_intervals[rank]
            # Fix the ranges with CIs that exceed the previous rank's CIs
            if max_ci > p_min:
                max_ci = p_min
                min_ci = max_ci - median(ci_ranges)
        p_min, p_max = min_ci, max_ci
        validated_intervals[rank] = (min_ci, max_ci)
        # Fill in those that are missing; estimate the length corresponding to a rank's range
    return validated_intervals


def main(args: list):
    newick_tree_file, tax_ids_file = args[:]
    taxa = tax_ids_file_to_leaves(tax_ids_file)
    leaf_taxa_map = dict()
    for taxa_leaf in taxa:
        leaf_taxa_map[taxa_leaf.number] = taxa_leaf.lineage
    # Read the NEWICK-formatted phylogenetic tree
    tree = Tree(newick_tree_file)
    taxonomic_rank_intervals = bound_taxonomic_branch_distances(tree, leaf_taxa_map)
    # taxonomic_rank_intervals = validate_rank_intervals(taxonomic_rank_intervals)
    sys.stdout.write("Rank\tMin.   - Max.\n")
    ranks = {1: "Phylum", 2: "Class", 3: "Order", 4: "Family", 5: "Genus", 6: "Species"}
    for depth in sorted(ranks, key=int, reverse=True):
        rank = ranks[depth]
        sys.stdout.write(rank + "\t" + ' - '.join([str(x) for x in taxonomic_rank_intervals[rank]]) + "\n")
    sys.stdout.flush()


if __name__ == "__main__":
    main(sys.argv[1:])
