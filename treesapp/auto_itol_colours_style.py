#!/usr/bin/env python3

__author__ = 'Connor Morgan-Lang'

import sys
import os
import re
import argparse
import logging
import ete3
from collections import namedtuple
from seaborn import color_palette
from treesapp.classy import prep_logging
from treesapp.phylo_seq import TreeLeafReference
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
        self.style_output = ""
        self.strip_output = ""
        self.tax_ids_file = tax_ids_table
        self.num_seqs = 0
        self.num_taxa = 0
        self.tree_leaves = list()  # List of 'TreeLeafReference's
        self.taxon_leaf_map = dict()  # taxon -> list(TreeLeafReference.number)
        self.unique_clades = dict()  # Unique numeric ID -> taxon
        self.internal_node_map = dict()  # Internal node -> child leaf names

    def get_clades(self, target_depth: int) -> None:
        """
        Loads unique taxa found in the tree into a dictionary where the values are all leaves of a taxon.
        num_taxa is incremented here and represents the number of unique taxa.
        The number of keys in taxon_leaf_map represents the number of

        :param target_depth: A numeric value representing the taxonomic rank threshold, i.e., 0 is Kingdom, 7 is species
        :return: None
        """
        clades = dict()
        truncated_lineages = dict()
        self.num_taxa = 0
        potential_domains = ["Archaea", "Bacteria", "Eukaryota"]
        potential_luca = ["Root", "cellular organisms"]

        for leaf in sorted(self.tree_leaves, key=lambda x: x.lineage):  # type: TreeLeafReference
            if not leaf.lineage:
                continue
            lineage_path = leaf.lineage.split('; ')
            if lineage_path[0] in potential_domains:
                depth = target_depth - 1
            elif lineage_path[0] in potential_luca:
                depth = target_depth
            else:
                logging.debug("Unable to find root or domain name for '%s'. Skipping.\n" % leaf.lineage)
                continue
            if len(lineage_path) > depth:
                taxon = lineage_path[depth]
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
                try:
                    truncated_lineages[leaf.lineage] += 1
                except KeyError:
                    truncated_lineages[leaf.lineage] = 1
        self.num_taxa += 1
        self.unique_clades = clades
        logging.info(self.marker + ": " + str(self.num_taxa) + " unique clades.\n" +
                     str(len(truncated_lineages)) + " truncated lineages (" +
                     str(sum(truncated_lineages.values())) + " sequences) that will not be coloured.\n")
        return

    def remove_taxon_from_colours(self, taxon):
        for n in self.unique_clades:
            if self.unique_clades[n] == taxon:
                self.unique_clades.pop(n)
                break
        self.taxon_leaf_map.pop(taxon)
        return

    def filter_unwanted_taxa(self, taxa_filter, target_depth):
        """
        A lineage-based string filter for excluding taxa from the output colour palette.
        Even though the palette is for taxa, and not complete lineages,
        we may want to filter out a very large group of organisms, such as Eukaryota.
        This facilitates that operation so not all taxa from that larger group need to be individually named
        :param taxa_filter:
        :param target_depth:
        :return:
        """
        filter_str = str(taxa_filter) + " taxa that won't be coloured:\n\t"
        filtered_taxa = set()
        filtered_seqs = 0
        filter_terms = taxa_filter.split(',')
        for leaf in self.tree_leaves:  # type: TreeLeafReference
            if len(filter_terms) > 0:
                for term in filter_terms:
                    if re.search(term, leaf.lineage):
                        lineage_path = leaf.lineage.split('; ')
                        taxon = lineage_path[target_depth]
                        filtered_taxa.add(taxon)
        for taxon in filtered_taxa:
            self.num_taxa -= 1
            filtered_seqs += len(self.taxon_leaf_map[taxon])
            self.num_seqs -= len(self.taxon_leaf_map[taxon])
            self.remove_taxon_from_colours(taxon)
        filter_str += ", ".join(filtered_taxa) + "\n"
        filter_str += "Sequences from unwanted taxa removed\t" + str(filtered_seqs) + "\n"
        filter_str += "Unwanted taxa removed\t" + str(len(filtered_taxa)) + "\n"
        filter_str += "Remaining sequences for colouring\t" + str(self.num_seqs) + "\n"
        filter_str += "Remaining taxa for colouring\t" + str(self.num_taxa) + "\n"
        logging.info(filter_str)
        return

    def filter_rare_groups(self, min_p: float):
        filter_str = "Taxa filtered due to proportion < " + str(min_p) + ":\n"
        filtered_taxa = []
        filtered_taxa_strings = []
        for taxon in self.taxon_leaf_map:
            leaf_nodes = set(self.taxon_leaf_map[taxon])
            if len(leaf_nodes)/self.num_seqs < min_p:
                filtered_taxa_strings.append("\t" + taxon + " = " + str(round(len(leaf_nodes)/self.num_seqs, 5)))
                filtered_taxa.append(taxon)
        if not filtered_taxa_strings:
            filtered_taxa_strings.append("\tNone")
        for taxon in filtered_taxa:
            self.num_taxa -= 1
            self.num_seqs -= len(self.taxon_leaf_map[taxon])
            self.remove_taxon_from_colours(taxon)
        filter_str += "\n".join(filtered_taxa_strings) + "\n"
        filter_str += "Low-abundance taxa removed\t" + str(len(filtered_taxa)) + "\n"
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

    def find_mono_clades(self) -> dict:
        taxa_internal_node_map = dict()
        failed_taxa = []
        for taxon in self.taxon_leaf_map:
            taxa_internal_node_map[taxon] = []
            taxon_leaves = set(self.taxon_leaf_map[taxon])
            # Sorting in descending order of number of child nodes
            for i_node in sorted(self.internal_node_map, key=lambda n: len(self.internal_node_map[n]), reverse=True):
                i_node_leaves = set(self.internal_node_map[i_node])
                # if "39" in i_node_leaves:
                #     print(i_node_leaves)
                if len(i_node_leaves) > len(taxon_leaves):
                    continue
                # taxon_leaves can either be a superset or match
                # Detect whether the internal node contains a set of only the leaves belonging to the taxon
                if taxon_leaves.issuperset(i_node_leaves):
                    taxa_internal_node_map[taxon].append(i_node_leaves)
                    # Pop from the list of all taxon leaves
                    for leaf in i_node_leaves:
                        taxon_leaves.remove(leaf)
                else:
                    pass
                if len(taxon_leaves) == 0:
                    break
            if taxon_leaves:
                failed_taxa.append(taxon)
        if failed_taxa:
            logging.warning("Unable to find internal nodes for all children of:\n\t" +
                            ", ".join(failed_taxa) + "\n")
        return taxa_internal_node_map


def linearize_tree_leaves(tree_string):
    """
    Parses a tree by depth-first search and returns the list
    :return:
    """
    tree_root = ete3.Tree(tree_string)
    leaf_order = []
    for node in tree_root.traverse(strategy="postorder"):
        if node.name:
            leaf_order.append(node.name)
    return leaf_order


def convert_clades_to_ranges(taxa_clades: dict, leaf_order: list) -> dict:
    taxa_ranges = dict()
    for taxon in taxa_clades:
        taxa_ranges[taxon] = []
        for clade_leaves in taxa_clades[taxon]:
            if len(clade_leaves) > 1:
                positions = clade_leaf_sides(clade_leaves, leaf_order)
                taxa_ranges[taxon].append([positions.left, positions.right])
            else:
                taxa_ranges[taxon].append(clade_leaves)
    return taxa_ranges


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
                        help="Path to the output directory to write the output files. [ DEFAULT = ./ ]")
    optopt.add_argument('-p', "--palette", default="BrBG", required=False,
                        help="The Seaborn colour palette to use [ DEFAULT = BrBG ]")
    optopt.add_argument('-m', '--min_proportion', dest="min_prop", default=0.0, required=False, type=float,
                        help="Minimum proportion of sequences a group contains to assign colour [ DEFAULT = 0 ]")
    optopt.add_argument("--no_polyphyletic", dest="no_poly", default=False, action="store_true", required=False,
                        help="Flag forcing the omission of all polyphyletic taxa from the colours file.")
    optopt.add_argument("-f", "--filter", dest="taxa_filter", default="", required=False,
                        help="Keywords for excluding specific taxa from the colour palette.\n"
                             "[ DEFAULT is no filter ]")

    miscellaneous_opts = parser.add_argument_group("Miscellaneous options")
    miscellaneous_opts.add_argument('-v', '--verbose', action='store_true', default=False,
                                    help='Prints a more verbose runtime log')
    miscellaneous_opts.add_argument("-h", "--help",
                                    action="help",
                                    help="Show this help message and exit")

    args = parser.parse_args()

    if not args.output:
        args.output = os.path.abspath(os.path.curdir)
    if args.output[-1] != os.sep:
        args.output += os.sep

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


def create_write_file(file_name: str, text: str) -> None:
    # Open the output file
    try:
        file_handler = open(file_name, 'w')
    except IOError:
        logging.error("Unable to open file '" + file_name + "' for writing!\n")
        raise IOError

    # Write the iTOL-formatted header and node, colour and taxon fields
    file_handler.write(text)
    file_handler.close()
    return


def read_tax_ids_file(taxa_colours: TaxaColours) -> None:
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
        leaf_node = number + "_" + taxa_colours.marker
        leaf = TreeLeafReference(leaf_node, translation)
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
    colourless = set()

    colours_style_string = "TREE_COLORS\nSEPARATOR TAB\nDATA\n"
    for taxon in palette_taxa_map:
        col = palette_taxa_map[taxon]
        for leaf_num in taxa_colours.taxon_leaf_map[taxon]:
            colours_style_line = [str(leaf_num), "range", col, taxon]
            if len(colours_style_line) > 0:
                colours_style_string += "\t".join(colours_style_line) + "\n"
    # TODO: Find new way to track the uncoloured leaves

    logging.info("Output is available in " + taxa_colours.style_output + ".\n")
    logging.debug(str(len(colourless)) + " lineages were not assigned to a colour:\n\t" +
                  "\n\t".join(colourless) + "\n")

    create_write_file(taxa_colours.style_output, colours_style_string)

    return


def write_colour_strip(taxa_nodes: dict, palette_taxa_map: dict, colour_strip_file: str) -> None:
    colourless = set()

    colour_strip_text = "DATASET_COLORSTRIP\nSEPARATOR SPACE\nDATASET_LABEL clade_colours\n" +\
                        "STRIP_WIDTH 75\nMARGIN 0\nSHOW_INTERNAL 1\nDATA\n"
    for taxon in palette_taxa_map:
        col = palette_taxa_map[taxon]
        for leaf_range in taxa_nodes[taxon]:
            colours_style_line = ['|'.join(leaf_range), col, taxon]
            if len(colours_style_line) > 0:
                colour_strip_text += " ".join(colours_style_line) + "\n"

    logging.info("Output is available in " + colour_strip_file + ".\n")
    logging.debug(str(len(colourless)) + " lineages were not assigned to a colour:\n\t" +
                  "\n\t".join(colourless) + "\n")

    create_write_file(colour_strip_file, colour_strip_text)
    return


def init_taxa_colours(args):
    if re.match(".*tax_ids_(.*).txt$", args.tax_ids_table):
        taxa_colours = TaxaColours(args.tax_ids_table)
        logging.debug("\tGenerating colour palette for " + taxa_colours.marker + ".\n")
    else:
        logging.error("Name of tax_ids file is formatted oddly.\n")
        sys.exit(2)

    if not os.path.isdir(args.output):
        os.makedirs(args.output)

    taxa_colours.style_output = args.output + taxa_colours.marker + "_colours_style.txt"
    taxa_colours.strip_output = args.output + taxa_colours.marker + "_colour_strip.txt"

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
            target_depth = depth
            break
    if target_depth == 0:
        logging.error("Rank '" + args.rank + "' not accepted!\n")
        sys.exit(3)
    return target_depth


def clade_leaf_sides(clade_leaves, leaf_order: list) -> namedtuple:
    positions = namedtuple("positions",
                           ["left", "right"])
    indices = dict()
    for leaf in clade_leaves:
        indices[leaf_order.index(leaf)] = leaf
    p = positions
    p.left = indices[min(indices.keys())]
    p.right = indices[max(indices.keys())]
    return p


def order_taxa(taxon_leaf_map: dict, leaf_order: list):
    leftmost_taxon = dict()
    taxon_order = dict()
    indices = dict()  # Could just use a list but may as well keep track of the left-most node with dict
    for taxon in taxon_leaf_map:
        for leaf_name in taxon_leaf_map[taxon]:
            indices[int(leaf_order.index(leaf_name))] = leaf_name
        leftmost_leaf = min(indices.keys())
        leftmost_taxon[leftmost_leaf] = taxon
        indices.clear()
    i = 0
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
    if args.taxa_filter:
        taxa_colours.filter_unwanted_taxa(args.taxa_filter, target_depth)
    if args.no_poly:
        # Optionally not colour polyphyletic clades based on args.no_poly
        taxa_colours.filter_polyphyletic_groups()
    if args.min_prop:
        taxa_colours.filter_rare_groups(args.min_prop)
    leaf_order = linearize_tree_leaves(args.tree)
    colours = get_colours(args, taxa_colours, args.palette)
    # Sort the nodes by their internal node order
    taxa_order = order_taxa(taxa_colours.taxon_leaf_map, leaf_order)
    palette_taxa_map = map_colours_to_taxa(taxa_order, colours)
    write_colours_styles(taxa_colours, palette_taxa_map)
    # Find the minimum set of monophyletic internal nodes for each taxon
    taxa_clades = taxa_colours.find_mono_clades()
    taxa_ranges = convert_clades_to_ranges(taxa_clades, leaf_order)
    write_colour_strip(taxa_ranges, palette_taxa_map, taxa_colours.strip_output)


if __name__ == "__main__":
    main()
