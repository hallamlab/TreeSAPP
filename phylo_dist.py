#!/usr/bin/env python3

import sys
import re
from ete3 import Tree
from file_parsers import tax_ids_file_to_leaves
from utilities import clean_lineage_string

from scipy import stats
import scipy as sp
import numpy as np

# TODO: replace mean_confidence_interval to remove scipy and numpy dependency

__author__ = 'Connor Morgan-Lang'


def confidence_interval(data, confidence=0.80):
    a = 1.0*np.array(data)
    n = len(a)
    m, se = np.mean(a), stats.sem(a)
    h = se * sp.stats.t._ppf((1+confidence)/2., n-1)
    return m-h, m+h


def rank_recommender(phylogenetic_distance: float, taxonomic_rank_intervals: dict):
    ranks = ["Kingdom", "Phylum", "Class", "Order",
             "Family", "Genus", "Species", "Strain"]
    depth = 1  # Start at Phylum
    while depth < 7:
        rank = ranks[depth]
        if taxonomic_rank_intervals[rank]:
            min_dist, max_dist = taxonomic_rank_intervals[rank]
            if min_dist < phylogenetic_distance < max_dist:
                depth += 1
                break
            elif phylogenetic_distance > max_dist:
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


def prune_branches(tree: Tree, leaf_taxa_map: dict, rank="Genus"):
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
                clade = monophyletic_clades[lineage][para_group]
                if len(clade) == 1:
                    # Unable to calculate distance with a single leaf
                    continue
                lca = clade[0].get_common_ancestor(clade[1:])
                lca_nodes[depth].append(lca)
                if depth+1 in lca_nodes:
                    if lca not in lca_nodes[depth+1]:
                        taxonomic_rank_distances[rank] += parent_to_tip_distances(lca, clade, True)
                else:
                    taxonomic_rank_distances[rank] += parent_to_tip_distances(lca, clade, True)
        if len(taxonomic_rank_distances[rank]) >= 5:
            min_dist, max_dist = confidence_interval(taxonomic_rank_distances[rank])
            taxonomic_rank_intervals[rank] = (round(min_dist, 4), round(max_dist, 4))
        else:
            taxonomic_rank_intervals[rank] = ()
    return taxonomic_rank_intervals


def main(args: list):
    newick_tree_file, tax_ids_file = args[:]
    taxa = tax_ids_file_to_leaves(tax_ids_file)
    leaf_taxa_map = dict()
    for taxa_leaf in taxa:
        leaf_taxa_map[taxa_leaf.number] = taxa_leaf.lineage
    # Read the NEWICK-formatted phylogenetic tree
    tree = Tree(newick_tree_file)
    taxonomic_rank_intervals = bound_taxonomic_branch_distances(tree, leaf_taxa_map)
    sys.stdout.write("Rank\tMin.   - Max.\n")
    ranks = {1: "Phylum", 2: "Class", 3: "Order", 4: "Family", 5: "Genus", 6: "Species"}
    for depth in sorted(ranks, key=int, reverse=True):
        rank = ranks[depth]
        sys.stdout.write(rank + "\t" + ' - '.join([str(x) for x in taxonomic_rank_intervals[rank]]) + "\n")
    sys.stdout.flush()


if __name__ == "__main__":
    main(sys.argv[1:])
