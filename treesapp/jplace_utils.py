__author__ = 'Connor Morgan-Lang'

import sys
import re
import glob
import os
import logging
from json import load, loads, dumps

from treesapp.phylo_seq import JPlace, PQuery


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


def jplace_parser(filename: str) -> JPlace:
    """
    Parses the jplace file using the load function from the JSON library

    :param filename: jplace file output by RAxML
    :return: JPlace object
    """
    itol_datum = JPlace()
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


def demultiplex_pqueries(jplace_data: JPlace, pquery_map=None) -> list:
    """
    Demultiplexes each placed query sequence (PQuery) into its own PQuery instance,
    copying over all JPlace information and the set of possible placements.

    :param jplace_data: A JPlace instance
    :param pquery_map: A dictionary mapping placed query sequence names to their respective PQuery instances
    :return: List of PQuery instances
    """
    tree_placement_queries = list()
    for pquery in jplace_data.placements:
        pquery_obj = PQuery()
        # Copy the essential information to the PQuery instance
        pquery_obj.placements = [pquery]
        pquery_obj.name_placed_sequence()

        if pquery_map:
            pquery_obj = pquery_map[pquery_obj.place_name]
            pquery_obj.placements = [pquery]
            # Placed sequence has already been named

        pquery_obj.transfer_metadata(jplace_data)
        pquery_obj.correct_decoding()

        tree_placement_queries.append(pquery_obj)

    return tree_placement_queries


def filter_jplace_data(jplace_data: JPlace, tree_saps: list):
    """
    Removes unclassified pqueries from the placements element in jplace_data

    :param jplace_data:
    :param tree_saps: List of TreeProtein objects
    :return:
    """
    jplace_data.correct_decoding()
    jplace_data.filter_max_weight_placement()
    new_placement_collection = list()
    sapling_map = {sapling.place_name: sapling for sapling in tree_saps}

    for d_place in jplace_data.placements:
        dict_strings = list()
        classified = False
        for k, v in loads(d_place).items():
            if k == 'n':
                # Find the TreeProtein that matches the placement (same contig name)
                try:
                    sapling = sapling_map[v[0]]
                except KeyError:
                    logging.error("Unable to find sequence '" + str(v[0]) + "' in sapling-map keys.\n")
                    sys.exit(5)
                # If the TreeProtein is classified, flag to append
                classified = sapling.classified
                dict_strings.append(dumps(k) + ":" + dumps(v))
            elif k == 'p':
                dict_strings.append(dumps(k) + ":" + dumps(v))
            else:
                logging.error("Unrecognized key '" + str(k) + "' in Jplace \"placements\".")
                sys.exit(9)
        if classified:
            new_placement_collection.append('{' + ', '.join(dict_strings) + '}')
    jplace_data.placements = new_placement_collection
    return jplace_data


def write_jplace(itol_datum: JPlace, jplace_file: str):
    """
    A hacky function for writing jplace files with concatenated placements
     which are also compatible with iTOL's jplace parser

    :param itol_datum: A JPlace class object
    :param jplace_file: A JPlace file path to write to
    :return:
    """

    if len(itol_datum.placements) == 0:
        return

    try:
        jplace_out = open(jplace_file, 'w')
    except IOError:
        logging.error("Unable to open " + jplace_file + " for writing.\n")
        sys.exit(9)

    itol_datum.correct_decoding()

    # Begin writing elements to the jplace file
    jplace_out.write('{\n\t"tree": "')
    jplace_out.write(itol_datum.tree + "\", \n")
    jplace_out.write("\t\"placements\": [\n\t")
    # The [] is lost from 'n': ["query"] during loads
    new_placement_collection = list()
    for d_place in itol_datum.placements:
        dict_strings = list()
        for k, v in loads(d_place).items():
            if k == 'n':
                dict_strings.append(dumps(k) + ":[" + dumps(v) + "]")
            elif k == 'p':
                dict_strings.append(dumps(k) + ":" + dumps(v))
            else:
                logging.error("Unrecognized key '" + str(k) + "' in Jplace \"placements\".")
                sys.exit(9)
        new_placement_collection.append('{' + ', '.join(dict_strings) + '}')
    jplace_out.write(",\n\t".join(new_placement_collection))
    jplace_out.write("\n\t],\n")
    jplace_out.write("\t\"metadata\": " + re.sub('\'', '"', str(itol_datum.metadata)) + ",\n")
    jplace_out.write("\t\"version\": " + str(itol_datum.version) + ",\n")
    jplace_out.write("\t\"fields\": [\n\t")
    jplace_out.write(", ".join(itol_datum.fields) + "\n\t]\n}\n")

    jplace_out.close()
    return


def organize_jplace_files(jplace_files: list) -> dict:
    """
    Given a list of jplace files, construct a dictionary of jplace files indexed by their refpkg prefix (denominator).

    :param jplace_files: list object containing full paths to jplace files generated by RAxML
    :return: Dictionary indexed by the refpkg prefix (marker code)
    """
    jplace_marker_re = re.compile(r".*epa_result.(.*)_hmm_purified_.*.jplace$")
    jplace_collection = dict()
    for filename in jplace_files:
        file_name_info = jplace_marker_re.match(filename)
        try:
            refpkg_name = file_name_info.group(1)
        except AttributeError:
            logging.error("Regex parsing marker information from jplace files was unsuccessful!\n"
                          "The offending file name: " + filename + "\n")
            sys.exit(7)

        if refpkg_name not in jplace_collection:
            jplace_collection[refpkg_name] = list()
        jplace_collection[refpkg_name].append(filename)
    return jplace_collection


def sub_indices_for_seq_names_jplace(jplace_dir, numeric_contig_index, refpkg_dict):
    """
    Ugly script for running re.sub on a set of jplace files

    :param jplace_dir:
    :param numeric_contig_index:
    :param refpkg_dict: A dictionary of ReferencePackage instances indexed by their prefix values
    :return:
    """
    jplace_files = glob.glob(jplace_dir + '*.jplace')
    jplace_collection = organize_jplace_files(jplace_files)
    for refpkg_name, jplace_files in jplace_collection.items():
        try:
            refpkg = refpkg_dict[refpkg_name]
        except KeyError:
            logging.warning("Intermediate files found from a previous run will be skipped:\n\t" +
                            "\n\t".join(jplace_files) + "\n")
            continue
        for jplace_path in jplace_files:
            jplace_data = jplace_parser(jplace_path)
            for pquery in jplace_data.placements:
                pquery["n"] = numeric_contig_index[refpkg.prefix][int(pquery["n"][0])]
            write_jplace(jplace_data, jplace_dir + os.sep + "tmp.jplace")
            os.rename(jplace_dir + os.sep + "tmp.jplace", jplace_path)
    return


def add_bipartitions(itol_datum, bipartition_file):
    """
    Adds bootstrap values read from a NEWICK tree-file where they are indicated by values is square-brackets
    and inserts them into a JPlace-formatted NEWICK tree (Internal nodes in curly braces are required)

    :param itol_datum: An JPlace object
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


def find_edge_length(jplace_tree: str, node: str):
    """
    Parses a Newick-formatted tree with depth-first-search internal nodes enveloped by curly braces
    :param jplace_tree: Newick tree
    :param node: The number of an internal node
    :return: float
    """
    edge_component = re.search(r':([0-9.]+)\{' + re.escape(node) + '}', jplace_tree)
    if edge_component:
        edge_length = float(edge_component.group(1))
    else:
        raise AssertionError("Unable to find node '" + node + "' in JPlace-formatted tree.")
    return edge_length
