#!/usr/bin/env python3

import sys
import logging
from ete3 import Tree

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


def regress_ranks(rank_distance_ranges: dict, taxonomic_ranks: dict) -> (float, float):
    """
    Uses linear regression to correlate phylogenetic distance with taxonomic rank

    :param rank_distance_ranges: A dictionary with taxonomic ranks as keys and a list of distances (floats) as values
    :param taxonomic_ranks: A dictionary mapping taxonomic ranks (keys) to their depth (int) in the taxonomic hierarchy
    :return: A tuple with the slope and intercept estimated by a linear regression
    """

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
        logging.error("Only {} ranks available for modelling.\n".format(len(depth_dist_dict.keys())))
        sys.exit(33)

    dist_list.clear()
    for depth in sorted(depth_dist_dict, key=int):
        rank_depth_list += [depth] * len(depth_dist_dict[depth])
        dist_list += list(depth_dist_dict[depth])

    # For TreeSAPP predictions
    # opt_slope, intercept = [round(float(x), 4) for x in list(np.polyfit(dist_list, rank_depth_list, 1))]
    intercept = 8.0
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

    slope, intercept = [float(i) for i in taxonomic_rank_pfit]
    depth = int(round(phylo_dist*slope + intercept))

    return depth


def parent_to_tip_distances(parent: Tree, children: Tree, estimate=False):
    """
    Function utilizing ete3's tree object for calculating distances between a reference node (parent)
     and query nodes (children).
    The `estimate` flag will cause the parent's edge length to be included in the distance calculation.

    :param parent: A reference node Tree instance
    :param children: A list of query nodes, also Tree instances
    :param estimate: Boolean indicating whether these distances are to be used for estimating the edge length ranges
    :return: List of all branch distances between the parent node and the tips
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
