#!/usr/bin/env python3

__author__ = 'Connor Morgan-Lang'

import sys
import os
import re
import argparse
import logging
import ete3
from seaborn import color_palette
from treesapp.classy import TreeLeafReference, prep_logging
from treesapp.utilities import clean_lineage_string
from treesapp.entish import map_internal_nodes_leaves

rank_depth_map = {0: "Cellular organisms", 1: "Kingdom",
                  2: "Phylum", 3: "Class", 4: "Order",
                  5: "Family", 6: "Genus", 7: "Species",
                  8: "Strain"}


class Clade:
    def __init__(self):
        self.taxon = ""
        self.colour = ""
        self.leaves = []  # Leaf names (not internal nodes)
        self.i_nodes = []  # List of all internal nodes that
        return


class TaxaColours:
    def __init__(self, tax_ids_table):
        self.marker = re.match("tax_ids_(.*).txt$", os.path.basename(tax_ids_table)).group(1)
        self.output = ""
        self.tax_ids_file = tax_ids_table
        self.num_seqs = 0
        self.num_taxa = 0
        self.tree_leaves = list()  # List of 'TreeLeafReference's
        self.taxon_leaf_map = dict()  # taxon -> list(TreeLeafReference.number)
        self.unique_clades = dict()  # Unique numeric ID -> taxon
        self.internal_node_map = dict()  # Internal node -> child leaf names

    def get_clades(self, target_depth):
        clades = dict()
        self.num_taxa = 0

        for leaf in sorted(self.tree_leaves, key=lambda x: x.lineage):  # type: TreeLeafReference
            if not leaf.lineage:
                continue
            lineage_path = clean_lineage_string(leaf.lineage).split('; ')
            if len(lineage_path) > target_depth:
                taxon = lineage_path[target_depth]
                try:
                    self.taxon_leaf_map[taxon].append(leaf.number)
                except KeyError:
                    self.taxon_leaf_map[taxon] = [leaf.number]
                if len(clades) == 0 and self.num_taxa == 0:
                    clades[self.num_taxa] = taxon
                elif clades[self.num_taxa] != taxon:
                    self.num_taxa += 1
                    clades[self.num_taxa] = taxon
            else:
                pass
        self.num_taxa += 1
        self.unique_clades = clades

        logging.debug("\t" + self.marker + ": " + str(self.num_taxa) + " unique clades.\n")

        return

    def remove_taxon_from_colours(self, taxon):
        for n in self.unique_clades:
            if self.unique_clades[n] == taxon:
                self.unique_clades.pop(n)
                break
        self.taxon_leaf_map.pop(taxon)
        return

    def filter_rare_groups(self, min_p: float):
        filter_str = "Taxa filtered due to proportion < " + str(min_p) + "\n"
        filtered_taxa = []
        for taxon in self.taxon_leaf_map:
            leaf_nodes = set(self.taxon_leaf_map[taxon])
            if len(leaf_nodes)/self.num_seqs < min_p:
                filter_str += "\t" + taxon + "= " + str(round(len(leaf_nodes)/self.num_seqs, 5)) + "\n"
                filtered_taxa.append(taxon)
        for taxon in filtered_taxa:
            self.num_taxa -= 1
            self.num_seqs -= len(self.taxon_leaf_map[taxon])
            self.remove_taxon_from_colours(taxon)
        filter_str += "Remaining sequences for colouring\t" + str(self.num_seqs) + "\n"
        filter_str += "Remaining taxa for colouring\t" + str(self.num_taxa) + "\n"
        logging.info(filter_str)
        return

    def filter_polyphyletic_groups(self):
        for taxon in self.taxon_leaf_map:
            ancestral_node = None
            leaf_nodes = set(self.taxon_leaf_map[taxon])

            # Sort the internal nodes by the number of children (moving from shallow to deep)
            for i_node in sorted(self.internal_node_map, key=lambda n: len(self.internal_node_map[n])):
                if len(self.internal_node_map[i_node]) < len(leaf_nodes):
                    continue
                elif len(set(self.internal_node_map[i_node]).intersection(leaf_nodes)) == len(leaf_nodes):
                    ancestral_node = i_node
                    break
                else:
                    pass
            if not ancestral_node:
                logging.warning("Unable to find internal node ancestral to all children of " + taxon + ".\n")
                self.remove_taxon_from_colours(taxon)
                continue
            if len(self.internal_node_map[ancestral_node]) > len(leaf_nodes):
                logging.warning(taxon + " clade was found to be polyphyletic.\n")
                self.remove_taxon_from_colours(taxon)

        return


def linearize_tree_leaves(tree_string):
    """
    Parses a tree by depth-first search and returns the list
    :return:
    """
    tree_root = ete3.Tree(tree_string)
    leaf_order = []
    for node in tree_root.traverse():
        if node.name:
            leaf_order.append(node.name)
    return leaf_order


def get_arguments():
    parser = argparse.ArgumentParser(add_help=False)
    required_args = parser.add_argument_group("Required arguments")
    required_args.add_argument("-t", "--tree", required=True,
                               help="The tree for the reference package in NEWICK format.")
    required_args.add_argument("-l", "--tax_ids_table", required=True,
                               help="tax_ids table for the TreeSAPP marker under investigation. "
                                    "Found in data/tree_data/tax_ids*.txt.")

    optopt = parser.add_argument_group("Optional options")
    optopt.add_argument('-r', "--rank", default="Order", required=False,
                        help="The rank to generate unique colours for [ DEFAULT = 'Order' ]")
    optopt.add_argument("-o", "--output", default=None, required=False,
                        help="Name of the output file. [ DEFAULT = ./${RefPkg}_colours_style.txt ]")
    optopt.add_argument('-p', "--palette", default="BrBG", required=False,
                        help="The Seaborn colour palette to use [ DEFAULT = BrBG ]")
    optopt.add_argument('-m', '--min_proportion', dest="min_prop", default=0.0, required=False, type=float,
                        help="Minimum proportion of sequences a group contains to assign colour [ DEFAULT = 0 ]")
    optopt.add_argument("--no_polyphyletic", dest="no_poly", default=False, action="store_true", required=False,
                        help="Flag forcing the omission of all polyphyletic taxa from the colours file.")

    miscellaneous_opts = parser.add_argument_group("Miscellaneous options")
    miscellaneous_opts.add_argument('-v', '--verbose', action='store_true', default=False,
                                    help='Prints a more verbose runtime log')
    miscellaneous_opts.add_argument("-h", "--help",
                                    action="help",
                                    help="Show this help message and exit")

    args = parser.parse_args()

    return args


def validate_command(args):
    if sys.version_info < (2, 9):
        logging.error("Python version '" + sys.version_info + "' not supported.\n")
        sys.exit(3)

    if args.rank not in rank_depth_map.values():
        logging.error("Rank '" + args.rank + "' not accepted! Please choose one of the following:\n" +
                      ", ".join([rank_depth_map[depth] for depth in rank_depth_map]) + "\n")
        sys.exit(1)
    return


def read_tax_ids_file(taxa_colours):
    """

    :param taxa_colours: A TaxaColours instance
    :return: None
    """
    try:
        cog_tax_ids = open(taxa_colours.tax_ids_file, 'r', encoding='utf-8')
    except IOError:
        logging.error('Can\'t open ' + str(taxa_colours.tax_ids_file) + '!\n')
        sys.exit(3)
    leaves = list()
    for line in cog_tax_ids:
        line = line.strip()
        try:
            fields = line.split("\t")
        except ValueError:
            logging.error('.split(\'\\t\') on ' + str(line) +
                          " generated " + str(len(line.split("\t"))) + " fields.")
            sys.exit(9)
        if len(fields) == 2:
            number, translation = fields
            lineage = ""
        elif len(fields) == 3:
            number, translation, lineage = fields
        else:
            logging.error("Unexpected number of fields in " + taxa_colours.tax_ids_file +
                          ".\nInvoked .split(\'\\t\') on line " + str(line))
            raise ValueError
        leaf = TreeLeafReference(number, translation)
        if lineage:
            leaf.lineage = lineage
            leaf.complete = True
        leaves.append(leaf)
        taxa_colours.num_seqs += 1

    taxa_colours.tree_leaves = leaves
    cog_tax_ids.close()
            
    return


def get_colours(args, taxa_colours, palette):
    clades = set()
    for num in taxa_colours.unique_clades:
        clades.add(taxa_colours.unique_clades[num])

    if len(clades) <= 1:
        logging.error(str(len(clades)) + " unique clade(s) identified at rank '" + args.rank + "'\n" +
                      "Consider changing rank to something more specific.\n"
                      "It may also be worth glancing at " + args.tax_ids_table + " to see if its okay.\n")
        sys.exit(4)
    else:
        logging.debug("Identified " + str(len(clades)) + " unique clades at " + args.rank + " rank.\n")

    n_clade = len(clades)
    try:
        colours = color_palette(palette, n_clade).as_hex()
        logging.info("Using palette '" + args.palette + "' for styling.\n")
    except ValueError as e:
        logging.error(str(e) + "\n")
        sys.exit(4)

    if len(colours) != n_clade:
        logging.error("Bad colour parsing! len(colours) != number of clades.\n")
        sys.exit(8)
    return colours


def map_colours_to_taxa(taxa_order, colours):
    """
    Function for mapping
    :param taxa_order: Dictionary indexed by numerical order of the taxa
    :param colours:
    :return:
    """
    palette_taxa_map = dict()
    alpha_sort_i = 0

    if len(taxa_order) != len(colours):
        logging.error("Number of taxa (%d) and colours (%d) are not equal!\n" % (len(taxa_order), len(colours)))
        sys.exit(7)

    for taxon_num in sorted(taxa_order, key=int):
        taxon = taxa_order[taxon_num]
        palette_taxa_map[taxon] = colours[alpha_sort_i]
        alpha_sort_i += 1

    return palette_taxa_map


def write_colours_styles(taxa_colours: TaxaColours, palette_taxa_map: dict) -> None:
    colours_style_header = "TREE_COLORS\nSEPARATOR TAB\nDATA\n"

    colourless = set()

    colours_style_string = ""
    for taxon in palette_taxa_map:
        col = palette_taxa_map[taxon]
        for leaf_num in taxa_colours.taxon_leaf_map[taxon]:
            colours_style_line = [str(leaf_num), "range", col, taxon]
            if len(colours_style_line) > 0:
                colours_style_string += "\t".join(colours_style_line) + "\n"
    # TODO: Find new way to track the uncoloured leaves

    # Open the output file
    try:
        cs_handler = open(taxa_colours.output, 'w')
    except IOError:
        logging.error("Unable to open " + taxa_colours.output + " for writing!\n")
        raise IOError

    logging.info("Output is available in " + taxa_colours.output + ".\n")
    logging.debug(str(len(colourless)) + " lineages were not assigned to a colour:\n\t" +
                  "\n\t".join(colourless) + "\n")
    # Write the iTOL-formatted header and node, colour and taxon fields
    cs_handler.write(colours_style_header)
    cs_handler.write(colours_style_string)
    cs_handler.close()

    return


def init_taxa_colours(args):
    if re.match(".*tax_ids_(.*).txt$", args.tax_ids_table):
        taxa_colours = TaxaColours(args.tax_ids_table)
        logging.debug("\tGenerating colour palette for " + taxa_colours.marker + ".\n")
    else:
        logging.error("Name of tax_ids file is formatted oddly.\n")
        sys.exit(2)

    if args.output:
        taxa_colours.output = args.output
    else:
        taxa_colours.output = taxa_colours.marker + "_colours_style.txt"

    # Create the internal-node to leaf-node map
    if not os.path.isfile(args.tree):
        logging.error("Unable to find tree file '" + args.tree + "'\n")
    tree_string = ""
    with open(args.tree) as tree_handler:
        for line in tree_handler:
            tree_string += line.strip()

    # add_internal_node function:
    annotated_tree = ""
    x = 0
    i_node = 0
    while x < len(tree_string):
        c = tree_string[x]
        if c in [',', ')', ';']:
            annotated_tree += '{' + str(i_node) + '}'
            i_node += 1
        annotated_tree += c
        x += 1
    annotated_tree = re.sub(r"\)[0-9.]+:", "):", annotated_tree)
    taxa_colours.internal_node_map = map_internal_nodes_leaves(annotated_tree)

    return taxa_colours


def find_rank_depth(args):
    target_depth = 0
    for depth in rank_depth_map:
        if rank_depth_map[depth] == args.rank:
            target_depth = depth - 1
            break
    if target_depth < 0:
        logging.error("Rank '" + args.rank + "' not accepted!\n")
        sys.exit(3)
    return target_depth


def order_taxa(taxon_leaf_map: dict, leaf_order: list):
    leftmost_taxon = dict()
    taxon_order = dict()
    indices = dict()  # Could just use a list but may as well keep track of the left-most node with dict
    i = 0
    for taxon in taxon_leaf_map:
        for leaf_name in taxon_leaf_map[taxon]:
            indices[leaf_order.index(leaf_name)] = leaf_name
        leftmost_leaf = min(indices.keys())
        leftmost_taxon[leftmost_leaf] = taxon
        indices.clear()
    for index in sorted(leftmost_taxon):
        taxon_order[i] = leftmost_taxon[index]
        i += 1
    return taxon_order


def main():
    args = get_arguments()

    log_file_name = os.path.abspath("./TreeSAPP_auto-colour_log.txt")
    prep_logging(log_file_name, args.verbose)
    logging.info("\n##\t\t\tGenerating colour-style file for iTOL\t\t\t##\n")
    validate_command(args)

    taxa_colours = init_taxa_colours(args)  # type: TaxaColours
    read_tax_ids_file(taxa_colours)
    target_depth = find_rank_depth(args)
    taxa_colours.get_clades(target_depth)
    if args.no_poly:
        # Optionally not colour polyphyletic clades based on args.no_poly
        taxa_colours.filter_polyphyletic_groups()
    if args.min_prop:
        taxa_colours.filter_rare_groups(args.min_prop)
    leaf_order = linearize_tree_leaves(args.tree)
    colours = get_colours(args, taxa_colours, args.palette)
    # Sort the nodes by their internal node order
    taxon_order = order_taxa(taxa_colours.taxon_leaf_map, leaf_order)
    palette_taxa_map = map_colours_to_taxa(taxon_order, colours)
    write_colours_styles(taxa_colours, palette_taxa_map)


if __name__ == "__main__":
    main()
