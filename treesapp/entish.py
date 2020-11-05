__author__ = 'Connor Morgan-Lang'

import sys
import re
import os
import logging

from ete3 import Tree, TreeNode


def get_node(tree: str, pos: int) -> (int, int):
    """
    Retrieves an internal node name from a Newick tree string

    :param tree: A string of an already read Newick tree file. Tree text exists on a single line.
    :param pos: Position in the string to start parsing from
    :return: Tuple of an integer representing the node and an integer for the new position in the string
    """
    node = ""
    pos += 1
    c = tree[pos]
    while c != '}':
        node += c
        pos += 1
        c = tree[pos]
    return int(node), pos


def label_internal_nodes_ete(ete_tree: Tree) -> None:
    i = 0
    if len(ete_tree.children) > 2:
        ete_tree.resolve_polytomy(recursive=True)
    for n in ete_tree.traverse(strategy="postorder"):  # type: Tree
        # Name the edge by it's unique internal node number
        if not n.name:
            n.name = str(i)
        i += 1
    return


def edge_from_node_name(ete_tree: Tree, node_name) -> int:
    edge_name = 0
    if len(ete_tree.children) > 2:
        ete_tree.resolve_polytomy(recursive=False)
    ete_tree = ete_tree.get_tree_root()
    for node in ete_tree.traverse(strategy="postorder"):  # type: Tree
        if len(node.children) > 2:
            node.resolve_polytomy(recursive=False)
        if str(node_name) == node.name:
            return edge_name
        edge_name += 1
    return -1


def get_ete_edge(ete_tree: Tree, edge_name) -> (TreeNode, TreeNode):
    """
    Traverses an ETE3 Tree structure in post-order, looking to match the desired edge_num to the current edge number,
    which is equal to node.number. Edge numbers are zero-indexed.

    :param ete_tree: An ETE3 Tree instance
    :param edge_name: An integer representing the desired edge number.
    :return: A tuple of the two immediately adjacent TreeNode instances for the corresponding branch/edge.
    """

    edge_n = 0
    if len(ete_tree.children) > 2:
        ete_tree.resolve_polytomy(recursive=False)
    ete_tree = ete_tree.get_tree_root()
    for node in ete_tree.traverse(strategy="postorder"):  # type: Tree
        if len(node.children) > 2:
            node.resolve_polytomy(recursive=False)
        if edge_n == int(edge_name):
            return node.up, node
        edge_n += 1
    return


def find_mean_pairwise_distances(children):
    pairwise_dists = list()
    for rleaf in children:
        for qleaf in children:
            if rleaf.name != qleaf.name:
                pairwise_dists.append(rleaf.get_distance(qleaf))
    return sum(pairwise_dists) / len(pairwise_dists)


def get_tip_distances(parent_node):
    children = parent_node.get_leaves()
    distances = [parent_node.get_distance(child) for child in children]
    return distances


# def find_cluster(lost_node: Tree, intra_distances=None):
#     """
#     Recursively calculates whether a node in a tree is the root of a cluster,
#     where a cluster is defined as a sub-tree (clade) whose members satisfy the condition:
#      the distance to the parent of the current subtree's root multiplied by the logarithm base 2 of the number of
#      cousins is greater than the mean(intra-cluster root-to-tip distances).
#     Adding a large number of members to the clade is penalized, as well as large distances from the existing subtree
#     to a new parent.
#
#     :param lost_node: A node within a tree, for which we want to orient
#     :param intra_distances: A list with the current set of leaf-tip distances
#     :return: Tree node, a list of float distances
#     """
#     if intra_distances is None:
#         intra_distances = []
#     parent = lost_node.up
#     if lost_node.is_root() or parent.is_root():
#         return lost_node, intra_distances
#
#     if not intra_distances:
#         # If this is the initial attempt at finding lost_node's cluster, find the intra-cluster leaf distances
#         intra_distances = get_tip_distances(lost_node)
#
#     # Penalty for increasing the size of the clade is log-base 2
#     cousins = lost_node.get_sisters()[0].get_leaf_names()
#     parent_dist = parent.get_distance(lost_node)
#     cost = parent_dist * log2(len(cousins) + 1)
#
#     if mean(intra_distances) > cost:
#         return lost_node, intra_distances
#
#     # Add the distance from the parent to the
#     intra_distances = [dist+parent_dist for dist in intra_distances]
#     for cousin in cousins:
#         intra_distances.append(parent.get_distance(cousin))
#     return find_cluster(parent, intra_distances)


def load_ete3_tree(newick_tree) -> Tree:
    if isinstance(newick_tree, str):
        return Tree(re.sub(r"{\d+}", '', newick_tree))
    else:
        return newick_tree


def map_internal_nodes_leaves(tree: str) -> dict:
    """
    Loads a Newick-formatted tree into a dictionary of all internal nodes (keys) and a list of child leaves (values).
    NOTE: the Newick tree must already contain internal nodes in braces. These trees are returned by EPA-NG.

    :param tree: A string of an already read Newick tree file. Tree text exists on a single line.
    :return: Dictionary of all internal nodes (keys) and a list of child leaves (values)
    """
    no_length_tree = re.sub(r":[0-9.]+(\[\d+\])?{", ":{", tree)
    node_map = dict()
    node_stack = list()
    leaf_stack = list()
    x = 0
    num_buffer = ""
    while x < len(no_length_tree):
        c = no_length_tree[x]
        if re.match(r"\d", c):
            while c != ':':
                num_buffer += c
                x += 1
                c = no_length_tree[x]
            node_stack.append([str(num_buffer)])
            num_buffer = ""
            x -= 1
        elif c == ':':
            # Append the most recent leaf
            current_node, x = get_node(no_length_tree, x + 1)
            if current_node in node_map:
                logging.error("Key '" + str(current_node) + "' being overwritten in internal-node map\n")
                sys.exit(11)
            node_map[current_node] = node_stack.pop()
            leaf_stack.append(current_node)
        elif c == ')':
            # Set the child leaves to the leaves of the current node's two children
            while c == ')' and x < len(no_length_tree):
                if no_length_tree[x + 1] == ';':
                    break
                while c != '{':
                    x += 1
                    c = no_length_tree[x]
                current_node, x = get_node(no_length_tree, x)
                if current_node in node_map:
                    logging.error("Key '" + str(current_node) + "' being overwritten in internal-node map\n")
                    sys.exit(11)
                node_map[current_node] = node_map[leaf_stack.pop()] + node_map[leaf_stack.pop()]
                leaf_stack.append(current_node)
                x += 1
                c = no_length_tree[x]
        x += 1

    if node_stack:
        logging.error("Node stack not empty by end of loading internal-node map:\n" + str(node_stack) + "\n")
        sys.exit(11)

    # validate_internal_node_map(node_map)
    return node_map


def annotate_partition_tree(refpkg_name: str, leaf_nodes: list, bipart_tree: str):
    try:
        tree_file = open(bipart_tree, 'r')
    except (FileNotFoundError, IOError):
        raise IOError("Unable to open RAxML bipartition tree " + bipart_tree + " for reading.")

    tree = tree_file.readline()
    tree_file.close()
    for leaf_node in leaf_nodes:
        if not re.search(r"[,(]{0}_{1}".format(leaf_node.number, refpkg_name), tree):
            logging.warning("Unable to find '{}' in {}.\n"
                            "The bipartition tree will not be annotated"
                            " (no effect on reference package).\n".format(leaf_node.number + '_' + refpkg_name,
                                                                          bipart_tree))
            break
        tree = re.sub(r"[,(]{0}_{1}".format(leaf_node.number, refpkg_name),
                      '(' + leaf_node.description,
                      tree)

    tree_output_dir = os.path.dirname(bipart_tree)
    annotated_tree_name = tree_output_dir + os.sep + "RAxML_bipartitions_annotated." + refpkg_name
    try:
        annotated_tree = open(annotated_tree_name, 'w')
    except IOError:
        raise IOError("Unable to open the annotated RAxML tree " + annotated_tree_name + " for writing.")

    annotated_tree.write(tree)
    annotated_tree.close()

    return


def tree_leaf_distances(tree: Tree) -> (float, list):
    """
    Calculates the maximum branch length distance observed in the tree
    Compiles a list of all distances from the root to the leaves

    :param tree: An ete3 Tree instance
    :return:
    """
    leaf_distances = []
    for leaf in tree.get_leaves():
        leaf_distances.append(tree.get_distance(leaf))
    max_dist = max([n.dist for n in tree.traverse("postorder")])
    return max_dist, leaf_distances


def index_tree_edges(tree: str):
    edge_index = dict()
    dist = ""
    edge = ""
    i = 0
    n = len(tree)
    while i < n:
        if tree[i] in [':', '{']:
            i += 1
            if dist:
                while re.match(r"[0-9]", tree[i]):
                    edge += tree[i]
                    i += 1
                edge_index[edge] = float(dist)
                dist = ""
                edge = ""
            else:
                while re.match(r"[0-9.]", tree[i]):
                    dist += tree[i]
                    i += 1
        else:
            i += 1

    return edge_index


def verify_bifurcations(newick_tree: str) -> str:
    ete_tree = load_ete3_tree(newick_tree)
    ete_tree.resolve_polytomy(recursive=True)
    return ete_tree.write(format=1)
