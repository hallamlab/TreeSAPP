__author__ = 'Connor Morgan-Lang'

import sys
import re
from classy import ItolJplace, TreeProtein
from json import load, loads
from utilities import clean_lineage_string


def children_lineage(leaves, pquery, node_map):
    """
    From the jplace placement field ()
    :param leaves:
    :param pquery:
    :param node_map:
    :return:
    """
    children = list()
    placement = loads(pquery, encoding="utf-8")
    loci = placement['p']
    for locus in loci:
        jplace_node = locus[0]
        tree_leaves = node_map[jplace_node]
        for tree_leaf in tree_leaves:
            # ref_leaf is a TreeLeafReference object
            for ref_leaf in leaves:
                if ref_leaf.number == tree_leaf:
                    if ref_leaf.complete:
                        children.append(clean_lineage_string(ref_leaf.lineage))
                    else:
                        children.append(clean_lineage_string(ref_leaf.description))
    return children


def pquery_likelihood_weight_ratio(pquery, position):
    """
    Determines the likelihood weight ratio (LWR) for a single placement. There may be multiple placements
    (or 'pquery's) in a single .jplace file. Therefore, this function is usually looped over.
    :param pquery:
    :param position: The position of "like_weight_ration" in the pquery fields
    :return: The float(LWR) of a single placement
    """
    lwr = 0.0
    placement = loads(str(pquery), encoding="utf-8")
    for k, v in placement.items():
        if k == 'p':
            acc = 0
            while acc < len(v):
                pquery_fields = v[acc]
                lwr = float(pquery_fields[position])
                acc += 1
    return lwr


def jplace_parser(filename):
    """
    Parses the jplace file using the load function from the JSON library
    :param filename: jplace file output by RAxML
    :return: ItolJplace object
    """
    itol_datum = ItolJplace()
    with open(filename) as jplace:
        jplace_dat = load(jplace, encoding="utf-8")
        itol_datum.tree = jplace_dat["tree"]
        # A list of strings
        if sys.version_info > (2, 9):
            itol_datum.fields = jplace_dat["fields"]
        else:
            itol_datum.fields = [x.decode("utf-8") for x in jplace_dat["fields"]]
        itol_datum.version = jplace_dat["version"]
        itol_datum.metadata = jplace_dat["metadata"]
        # A list of dictionaries of where the key is a string and the value is a list of lists
        itol_datum.placements = jplace_dat["placements"]

    jplace_dat.clear()

    return itol_datum


def demultiplex_pqueries(jplace_data):
    tree_placement_queries = list()
    for pquery in jplace_data.placements:
        pquery_obj = TreeProtein()
        pquery_obj.transfer(jplace_data)
        pquery_obj.placements = [pquery]
        pquery_obj.correct_decoding()
        tree_placement_queries.append(pquery_obj)

    return tree_placement_queries


def write_jplace(itol_datum, jplace_file):
    """
    A hacky function for writing jplace files with concatenated placements
     which are also compatible with iTOL's jplace parser
    :param itol_datum: A ItolJplace class object
    :param jplace_file:
    :return:
    """
    try:
        jplace_out = open(jplace_file, 'w')
    except IOError:
        raise IOError("Unable to open " + jplace_file + " for writing! Exiting now.")

    itol_datum.correct_decoding()
    # itol_datum.filter_min_weight_threshold(0.3)
    itol_datum.filter_max_weight_placement()

    # Begin writing elements to the jplace file
    jplace_out.write('{\n\t"tree": "')
    jplace_out.write(itol_datum.tree + "\", \n")
    jplace_out.write("\t\"placements\": [\n\t")
    jplace_out.write(", ".join(itol_datum.placements))
    jplace_out.write("\n\t],\n")
    jplace_out.write("\t\"version\": " + str(itol_datum.version) + ",\n")
    jplace_out.write("\t\"fields\": [\n\t")
    jplace_out.write(", ".join(itol_datum.fields) + "\n\t]\n}")

    jplace_out.close()
    return


def add_bipartitions(itol_datum, bipartition_file):
    """
    Adds bootstrap values read from a NEWICK tree-file where they are indicated by values is square-brackets
    and inserts them into a JPlace-formatted NEWICK tree (Internal nodes in curly braces are required)
    :param itol_datum: An ItolJplace object
    :param bipartition_file: Path to file containing the bootstrapped tree
    :return: Updated ItolJPlace object
    """
    with open(bipartition_file) as bootstrap_tree_handler:
        bootstrap_tree = bootstrap_tree_handler.readline().strip()
        no_length_tree = re.sub(":[0-9.]+", '', bootstrap_tree)
        node_stack = list()
        internal_node_bipart_map = dict()
        x = 0
        i_node = 0
        num_buffer = ""
        while x < len(no_length_tree):
            c = no_length_tree[x]
            if c == '[':
                x += 1
                c = no_length_tree[x]
                while re.search(r"[0-9]", c):
                    num_buffer += c
                    x += 1
                    c = no_length_tree[x]
                bootstrap = num_buffer
                num_buffer = ""
                x -= 1
                internal_node_bipart_map[node_stack.pop()] = bootstrap
            elif re.search(r"[0-9]", c):
                while re.search(r"[0-9]", c):
                    x += 1
                    c = no_length_tree[x]
                node_stack.append(str(i_node))
                i_node += 1
                x -= 1
            elif c == ')':
                node_stack.append(str(i_node))
                i_node += 1
            x += 1
    x = 0
    bootstrapped_jplace_tree = ''
    num_buffer = ''
    while x < len(itol_datum.tree):
        c = itol_datum.tree[x]
        if c == '{':
            x += 1
            c = itol_datum.tree[x]
            while re.search(r"[0-9]", c):
                num_buffer += c
                x += 1
                c = itol_datum.tree[x]
            if num_buffer in internal_node_bipart_map.keys():
                bootstrapped_jplace_tree += '[' + internal_node_bipart_map[num_buffer] + ']'
            bootstrapped_jplace_tree += '{' + num_buffer
            num_buffer = ''
            x -= 1
        else:
            bootstrapped_jplace_tree += c
        x += 1
    itol_datum.tree = bootstrapped_jplace_tree
    return itol_datum
