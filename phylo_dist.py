#!/usr/bin/env python3

import sys
import re
from ete3 import Tree
from file_parsers import tax_ids_file_to_leaves
from utilities import clean_lineage_string, median

from scipy import stats
import scipy as sp
import numpy as np

# TODO: replace confidence_interval to remove scipy and numpy dependency

__author__ = 'Connor Morgan-Lang'


def confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), stats.sem(a)
    h = se * sp.stats.t._ppf((1 + confidence) / 2., n - 1)
    return round(float(m - h), 4), round(float(m + h), 4)


def rank_recommender(phylo_dist: float, taxonomic_rank_intervals: dict):
    """
    Determines the rank depth (for example Class == 2) a taxonomic lineage should be truncated to
     based on which rank distance range (in taxonomic_rank_intervals) phylo_dist falls into
    
    :param phylo_dist: Float < 1.0 representing the branch distance from the nearest node
    :param taxonomic_rank_intervals: Dictionary with rank keys (e.g. Class) and distance ranges (min, max) as values
    :return: int
    """    
    ranks = ["Kingdom", "Phylum", "Class", "Order",
             "Family", "Genus", "Species", "Strain"]

    if not taxonomic_rank_intervals:
        return 7

    depth = 2  # Start at Class
    while depth < 7:
        rank = ranks[depth]
        if taxonomic_rank_intervals[rank]:
            min_dist, max_dist = taxonomic_rank_intervals[rank]
            if min_dist < phylo_dist < max_dist:
                depth += 1
                break
            elif phylo_dist > max_dist:
                break
        depth += 1
    return depth


def trim_lineages_to_rank(leaf_taxa_map, rank):
    trimmed_lineage_map = dict()
    ranks = {"Kingdom": 1, "Phylum": 2, "Class": 3, "Order": 4, "Family": 5, "Genus": 6, "Species": 7}
    depth = ranks[rank]
    for node_name in leaf_taxa_map:
        lineage = clean_lineage_string(leaf_taxa_map[node_name])
        if re.search("unclassified|environmental sample", lineage, re.IGNORECASE):
            continue
        lineage = lineage.split("; ")
        if len(lineage) < depth:
            continue
        trimmed_lineage_map[node_name] = "; ".join(lineage[:depth])
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
    assert isinstance(tree, Tree)
    # Check to see if the two collections are comparable
    for leaf in tree:
        if leaf.name not in leaf_taxa_map.keys():
            sys.stderr.write(str(leaf.name) + " not found in leaf_taxa_map.\n")
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
        if type(child_node) is Tree:
            distal_length = parent.get_distance(child_node.name)
        elif type(child_node) is str:
            distal_length = parent.get_distance(child_node)
        elif type(child_node) is int:
            distal_length = parent.get_distance(str(child_node))
        else:
            raise AssertionError("Cannot handle type '" + type(child_node) + "' for child.")
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
            min_dist, max_dist = confidence_interval(taxonomic_rank_distances[rank])
            taxonomic_rank_intervals[rank] = (round(min_dist, 4), round(max_dist, 4))
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
    # # As inclusive as possible and alignment trimming, 95% placement distance range
    # validated_intervals = {"Class": (0.4383, 0.509),
    #                        "Order": (0.3835, 0.4673),
    #                        "Family": (0.1854, 0.2183),
    #                        "Genus": (0.1383, 0.1637),
    #                        "Species": (0.0652, 0.0972)}
    # # Filtering those queries where a taxonomic ancestor is not present, 95% placement distance range
    # validated_intervals = {"Class": (0.4391, 0.5094),
    #                        "Order": (0.4759, 0.6276),
    #                        "Family": (0.1792, 0.2056),
    #                        "Genus": (0.1187, 0.1403),
    #                        "Species": (0.0447, 0.0935)}
    # # Optimized above - now need to automate via smoothing function?
    # validated_intervals = {"Class": (0.4391, 0.5094),
    #                        "Order": (0.206, 0.4390),
    #                        "Family": (0.13, 0.2056),
    #                        "Genus": (0.1187, 0.1403),
    #                        "Species": (0.0447, 0.0935)}
    # return validated_intervals

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
