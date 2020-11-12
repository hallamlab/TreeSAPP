__author__ = 'Connor Morgan-Lang'

import sys
import re
import glob
import os
import logging
from json import load, dumps

from treesapp.phylo_seq import PQuery, PhyloPlace, split_placements
from treesapp.entish import load_ete3_tree
from treesapp.phylo_dist import parent_to_tip_distances


# def pquery_likelihood_weight_ratio(pquery, position):
#     """
#     Determines the likelihood weight ratio (LWR) for a single placement. There may be multiple placements
#     (or 'pquery's) in a single .jplace file. Therefore, this function is usually looped over.
#
#     :param pquery:
#     :param position: The position of "like_weight_ration" in the pquery fields
#     :return: The float(LWR) of a single placement
#     """
#     lwr = 0.0
#     placement = loads(str(pquery), encoding="utf-8")
#     for k, v in placement.items():
#         if k == 'p':
#             acc = 0
#             while acc < len(v):
#                 pquery_fields = v[acc]
#                 lwr = float(pquery_fields[position])
#                 acc += 1
#     return lwr

class JPlace:
    """
    A class to hold all data relevant to a jplace file to be viewed in iTOL
    """
    ##
    # Information derived from Jplace pqueries:
    ##
    fields = list()
    pqueries = list()
    ref_name = ""
    file_name = ""
    tree = ""  # NEWICK tree
    metadata = ""
    version = ""  # Jplace version

    def __init__(self):
        return

    def summarize(self):
        """
        Prints a summary of the JPlace object (equivalent to a single marker) to stderr
        Summary include the number of marks found, the tree used, and the tree-placement of each sequence identified
        Written solely for testing purposes

        :return:
        """
        summary_string = "\nInformation for JPlace file '{}'\n" \
                         "{} sequence(s) grafted onto the {} tree.\n" \
                         "JPlace fields:\n\t{}\n" \
                         "Placement information:\n" \
                         "".format(self.file_name, len(self.pqueries), self.ref_name, self.fields)
        if not self.pqueries:
            summary_string += "\tNone.\n"
        elif self.pqueries[0] == '{}':
            summary_string += "\tNone.\n"
        else:
            for placement in self.pqueries:  # type: PhyloPlace
                summary_string += placement.summary()
        summary_string += "\n"
        return summary_string

    def clear_object(self):
        self.pqueries.clear()
        self.fields.clear()
        self.tree = ""
        self.metadata = ""
        self.version = ""

    def write_jplace(self, jplace_file: str) -> None:
        """
        Writes A JPlace file with concatenated placements that is compatible with iTOL's JPlace viewer

        :param jplace_file: A JPlace file path to write to
        :return: None
        """

        if len(self.pqueries) == 0:
            return

        jplace_str = ""
        # Begin writing elements to the jplace file
        jplace_str += '{\n\t"tree": "'
        jplace_str += self.tree + "\", \n"

        # Format the placements as JSON
        jplace_str += "\t\"placements\": [\n\t"
        new_placement_collection = []
        for pquery in self.pqueries:  # type: PQuery
            if pquery.classified:
                new_placement_collection.append(dumps(PhyloPlace.format_pplace_to_jplace(pquery.placements)))
        if len(new_placement_collection) == 0:
            return
        jplace_str += ",\n\t".join(new_placement_collection)
        jplace_str += "\n\t],\n"

        jplace_str += "\t\"metadata\": " + re.sub('\'', '"', str(self.metadata)) + ",\n"
        jplace_str += "\t\"version\": " + str(self.version) + ",\n"
        jplace_str += "\t\"fields\": [\n\t"
        jplace_str += ", ".join(['"' + x + '"' for x in self.fields]) + "\n\t]\n}\n"

        try:
            jplace_out = open(jplace_file, 'w')
        except IOError:
            logging.error("Unable to open " + jplace_file + " for writing.\n")
            sys.exit(9)

        jplace_out.write(jplace_str)
        jplace_out.close()
        return


def jplace_parser(filename: str) -> JPlace:
    """
    Parses the jplace file using the load function from the JSON library

    :param filename: jplace file output by RAxML
    :return: JPlace object
    """
    jplace_data = JPlace()
    with open(filename) as jplace:
        jplace_json = load(jplace, encoding="utf-8")
        jplace_data.tree = jplace_json["tree"]
        # A list of strings
        if sys.version_info > (2, 9):
            jplace_data.fields = jplace_json["fields"]
        else:
            jplace_data.fields = [x.decode("utf-8") for x in jplace_json["fields"]]
        jplace_data.version = jplace_json["version"]
        jplace_data.metadata = jplace_json["metadata"]
        # A list of dictionaries of where the key is a string and the value is a list of lists
        jplace_data.pqueries = jplace_json["placements"]

    jplace_json.clear()

    return jplace_data


def demultiplex_pqueries(jplace_data: JPlace, pquery_map=None) -> list:
    """
    Demultiplexes each placed query sequence (PQuery) into its own PQuery instance,
     copying over all JPlace information and the set of possible placements.
    The format of the PQuery's placements is modified from a list of dictionaries to a list of PhyloPlace instances.

    :param jplace_data: A JPlace instance
    :param pquery_map: A dictionary mapping placed query sequence names to their respective PQuery instances
    :return: List of PQuery instances
    """
    tree_placement_queries = list()
    for pquery in jplace_data.pqueries:
        pquery_obj = PQuery()
        # Copy the essential information to the PQuery instance
        pquery_obj.placements = split_placements(pquery)
        pquery_obj.name_placed_sequence()

        if pquery_map:
            # Placed sequence has already been named
            mapped_pquery = pquery_map[pquery_obj.place_name]
            # Prevent rerunning split_placements
            mapped_pquery.placements = pquery_obj.placements
            pquery_obj = mapped_pquery

        # pquery_obj.transfer_metadata(jplace_data)

        tree_placement_queries.append(pquery_obj)

    return tree_placement_queries


def calc_pquery_mean_tip_distances(jplace_data: JPlace, internal_node_leaf_map: dict) -> None:
    """
    One of the attributes of a PhyloPlace instance that is not precalculated by EPA is the mean-tip distance from a
    query sequence's placement position on an edge. This must be calculated by TreeSAPP.
    This function calls PhyloPlace.calc_mean_tip_length() on each PhyloPlace object in all PQuery.placements collection
    for all PQuery's in a JPlace's pqueries collection. It's just a convenience function.
    """
    jplace_tree = load_ete3_tree(jplace_data.tree)
    parent_leaf_memoizer = dict()

    for pquery in jplace_data.pqueries:  # type: PQuery
        for pplace in pquery.placements:  # type: PhyloPlace
            # Populate the memoization table to speed up length calculations for large JPlace files
            leaf_children = internal_node_leaf_map[pplace.edge_num]
            if len(leaf_children) > 1:
                if int(pplace.edge_num) not in parent_leaf_memoizer:
                    parent_leaf_memoizer[int(pplace.edge_num)] = parent_to_tip_distances(
                        jplace_tree.get_common_ancestor(leaf_children), leaf_children
                    )
                # Find the distance away from this edge's bifurcation (if internal) or tip (if leaf)
                pplace.calc_mean_tip_length(internal_leaf_node_map=internal_node_leaf_map, ref_tree=jplace_tree,
                                            memoization_map=parent_leaf_memoizer)
            else:
                pplace.mean_tip_length = 0.0
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


def sub_indices_for_seq_names_jplace(jplace_dir, numeric_contig_index, refpkg_dict) -> None:
    """
    Replaces the numerical identifier for each query sequence with their true names for all JPlace files in a directory

    :param jplace_dir: A directory containing JPlace files (extension is .jplace)
    :param numeric_contig_index: A dictionary mapping numerical identifiers to sequence name strings
    :param refpkg_dict: A dictionary of ReferencePackage instances indexed by their prefix values
    :return: None
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
            for pquery in jplace_data.pqueries:
                pquery["n"] = numeric_contig_index[refpkg.prefix][int(pquery["n"][0])]
            jplace_data.pqueries = demultiplex_pqueries(jplace_data)
            jplace_data.write_jplace(jplace_dir + os.sep + "tmp.jplace")
            os.rename(jplace_dir + os.sep + "tmp.jplace", jplace_path)
    return


def add_bipartitions(jplace_data, bipartition_file):
    """
    Adds bootstrap values read from a NEWICK tree-file where they are indicated by values is square-brackets
    and inserts them into a JPlace-formatted NEWICK tree (Internal nodes in curly braces are required)

    :param jplace_data: An JPlace object
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
    while x < len(jplace_data.tree):
        c = jplace_data.tree[x]
        if c == '{':
            x += 1
            c = jplace_data.tree[x]
            while re.search(r"[0-9]", c):
                num_buffer += c
                x += 1
                c = jplace_data.tree[x]
            if num_buffer in internal_node_bipart_map.keys():
                bootstrapped_jplace_tree += '[' + internal_node_bipart_map[num_buffer] + ']'
            bootstrapped_jplace_tree += '{' + num_buffer
            num_buffer = ''
            x -= 1
        else:
            bootstrapped_jplace_tree += c
        x += 1
    jplace_data.tree = bootstrapped_jplace_tree
    return jplace_data


def find_edge_length(jplace_tree: str, node: str):
    """
    Parses a Newick-formatted tree with depth-first-search internal nodes enveloped by curly braces

    :param jplace_tree: Newick tree
    :param node: The number of an internal node
    :return: float
    """
    edge_component = re.search(r':([0-9.]+){' + re.escape(node) + '}', jplace_tree)
    if edge_component:
        edge_length = float(edge_component.group(1))
    else:
        raise AssertionError("Unable to find node '" + node + "' in JPlace-formatted tree.")
    return edge_length
