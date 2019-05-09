#!/usr/bin/env python3

__author__ = 'Connor Morgan-Lang'

import sys
import os
import re
import argparse
import logging
# import seaborn as sns
from treesapp.classy import TreeLeafReference, prep_logging
from treesapp.utilities import clean_lineage_string
from treesapp.entish import map_internal_nodes_leaves

rank_depth_map = {0: "Cellular organisms", 1: "Kingdom",
                  2: "Phylum", 3: "Class", 4: "Order",
                  5: "Family", 6: "Genus", 7: "Species",
                  8: "Strain"}


class TaxaColours:
    def __init__(self, tax_ids_table):
        self.marker = re.match("tax_ids_(.*).txt$", os.path.basename(tax_ids_table)).group(1)
        self.output = self.marker + "_colours_style.txt"
        self.tax_ids_file = tax_ids_table
        self.num_seqs = 0
        self.num_taxa = 0
        self.tree_leaves = list()
        self.taxon_leaf_map = dict()
        self.unique_clades = dict()
        self.internal_node_map = dict()

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
        return

    def tree_order(self):
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
    optopt.add_argument('-p', "--palette", default="BrBG",
                        help="The ColorBrewer palette to use. Currently supported:"
                             " 'BrBG', 'PRGn', 'Set3', 'Greys', 'Paired'. [ DEFAULT = BrBG ]")
    optopt.add_argument("--no_polyphyletic", dest="no_poly", default=False, action="store_true", required=False,
                        help="Flag forcing the omission of all polyphyletic taxa from the colours file.")
    optopt.add_argument("-h", "--help",
                        action="help",
                        help="Show this help message and exit")

    miscellaneous_opts = parser.add_argument_group("Miscellaneous options")
    miscellaneous_opts.add_argument('-v', '--verbose', action='store_true', default=False,
                                    help='Prints a more verbose runtime log')

    args = parser.parse_args()

    return args


def validate_command(args):
    if sys.version_info < (2, 9):
        logging.error("Python version '" + sys.version_info + "' not supported.\n")
        sys.exit(3)

    palette_options = ['BrBG', 'PRGn', 'Set3', "Greys", "Paired"]
    if args.palette not in palette_options:
        logging.error("Palette '" + args.palette + "' is not supported. Please select one of the following options:\n" +
                      ', '.join(palette_options) + "\n")
        sys.exit(4)
    else:
        logging.info("Using palette '" + args.palette + "' for styling.\n")

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

    taxa_colours.tree_leaves = leaves
    cog_tax_ids.close()
            
    return


def get_colours(args, taxa_colours, palette):
    # TODO: Use seaborn so we are able to generate a palette on the fly

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

#    brbg = ["#543005", "#643906", "#744207", "#844C09", "#93570E", "#A16518", "#B07323",
#            "#BF812C", "#C89343", "#D1A65A", "#DAB871", "#E2C786", "#E8D29A", "#EFDDAE",
#            "#F5E7C2", "#F5EBD1", "#F5EFDF", "#F5F3ED", "#EEF3F2", "#E1F0EE", "#D4EDE9",
#            "#C7EAE5", "#B2E1DA", "#9ED9D0", "#8AD1C6", "#75C5B9", "#5FB5AB", "#4AA69D",
#            "#35978F", "#268981", "#177B73", "#086D65", "#006057", "#00544A", "#00483D",
#            "#003C30"]
    brbg = ["#543005", "#5D3505", "#663A06", "#6F3F07", "#784508", "#814A09", "#8A4F09",
            "#92570E", "#9A5E14", "#A36619", "#AB6E1F", "#B37625", "#BB7D2A", "#C28734",
            "#C79141", "#CC9C4E", "#D1A65B", "#D6B168", "#DBBB75", "#E0C481", "#E4CA8C",
            "#E7D098", "#EBD6A3", "#EFDCAE", "#F3E3B9", "#F5E8C4", "#F5EACC", "#F5ECD4",
            "#F5EEDC", "#F5F0E4", "#F5F2EC", "#F5F5F5", "#EDF3F2", "#E6F1EF", "#DEEFED",
            "#D7EDEA", "#CFECE8", "#C8EAE5", "#BDE6E0", "#B2E1DA", "#A6DCD4", "#9BD8CE",
            "#90D3C9", "#84CEC3", "#78C7BC", "#6CBFB4", "#60B6AC", "#54ADA3", "#48A49B",
            "#3C9C93", "#31938B", "#298B83", "#20847C", "#187C74", "#10746C", "#076C64",
            "#00645C", "#005D55", "#00574D", "#005046", "#00493E", "#004237", "#003C30"]
    set3 = ["#8DD3C7", "#A0DAC3", "#B3E1C0", "#C6E9BC", "#DAF0B9", "#EDF8B6", "#FDFDB3",
            "#F2F2BA", "#E7E6C0", "#DCDAC7", "#D1CFCE", "#C6C3D4", "#C0B8D6", "#CAAEC4",
            "#D4A4B3", "#DF9AA1", "#E9908F", "#F3867E", "#F48276", "#DF8A87", "#CB9397",
            "#B69BA8", "#A1A3B8", "#8CACC9", "#88B1CB", "#9DB1B8", "#B2B2A5", "#C8B291",
            "#DDB37E", "#F2B36B", "#F6B762", "#EABE63", "#DDC564", "#D1CC66", "#C4D467",
            "#B8DB68", "#BADC75", "#C6D98A", "#D3D69F", "#DFD3B4", "#EBD0C9", "#F8CDDE",
            "#F7CEE3", "#F1D0E1", "#EBD2DF", "#E6D4DD", "#E0D6DB", "#DAD8D9", "#D5CCD5",
            "#D0BDD0", "#CBAECB", "#C69FC6", "#C190C2", "#BC81BD", "#BE90BE", "#C1A2BF",
            "#C3B4C0", "#C6C6C2", "#C9D8C3", "#CCEBC5"]
    prgn = ["#40004B", "#490755", "#540F5F", "#5E176A", "#681F74", "#71267E", "#793187",
            "#803E8E", "#864B96", "#8D589D", "#9365A5", "#9A71AC", "#A27BB3", "#A985B9",
            "#B18FC0", "#B899C7", "#C0A3CD", "#C7ABD2", "#CEB4D7", "#D5BDDB", "#DCC6E0",
            "#E2CEE5", "#E8D6E9", "#EBDDEB", "#EEE3EE", "#F1EAF1", "#F4F0F4", "#F7F7F7",
            "#F1F5F0", "#EBF4E9", "#E6F3E3", "#E0F1DC", "#DBF0D5", "#D3EDCD", "#C9E9C3",
            "#C0E5BA", "#B7E2B1", "#ADDEA7", "#A3D99D", "#95D192", "#87C886", "#78C07A",
            "#6AB86F", "#5CAF63", "#50A65A", "#459C53", "#39924B", "#2D8843", "#227E3B",
            "#197434", "#146A2F", "#0F602A", "#0A5725", "#054D20", "#00441B"]
    greys = ["#FFFFFF", "#F6F6F6", "#ECECEC", "#DFDFDF", "#D1D1D1", "#C1C1C1", "#ACACAC",
             "#969696", "#828282", "#6E6E6E", "#5B5B5B", "#454545", "#2B2B2B", "#151515"]
    paired = ["#A6CEE3", "#69A7CD", "#2C80B8", "#529CA5", "#94CA92", "#92CF72", "#59B248",
              "#519F3C", "#AB9C6D", "#F99392", "#EF595A", "#E42022", "#ED5C3D", "#F9A662",
              "#FDAB4D", "#FE8E1B", "#F4892A", "#DCA08B", "#C0A6CF", "#9571B4", "#6A3D9A"]

    color_brewer_palettes = {"BrBG": brbg, "Set3": set3, "PRGn": prgn, "Greys": greys, "Paired": paired}

    brewed = color_brewer_palettes[palette]

    n_clade = len(clades)
    n_brewed = len(brewed)
    if n_clade > n_brewed:
        logging.error("Number of clades (" + str(n_clade) + ") exceeds the number of colours in " +
                      palette + " (" + str(n_brewed) + ")!\n" +
                      "\tUnfortunately we have hard-coded palettes and they do not scale.\n"
                      "\tWe are hoping to integrate a ColorBrewer API in the near future.\n"
                      "\tIn the meantime, would you please select a higher clade? Thank you!\n")
        sys.exit(7)
    else:
        colours_pre = list()
        colours_suf = list()
        colours = list()
        i = 0
        while i < n_brewed and len(colours) < n_clade:
            if i % 2:
                colours_pre.append(brewed[i])
            else:
                colours_suf.append(brewed[n_brewed-i-1])
            i += 1
            colours = colours_pre + colours_suf

    if len(colours) != n_clade:
        logging.error("Bad colour parsing! len(colours) != number of clades.\n")
        sys.exit(8)

    return colours


def map_colours_to_taxa(clades, lineages, colours):
    """
    Function for mapping
    :param clades:
    :param lineages:
    :param colours:
    :return:
    """
    palette_taxa_map = dict()
    alpha_sort_i = 0
    num_clades = len(clades)
    for lineage in sorted(lineages):
        i = 0
        while i < len(clades):
            clade = clades[i]
            if re.search(clean_lineage_string(clade), clean_lineage_string(lineage)):
                palette_taxa_map[clade] = colours[alpha_sort_i]
                clades.pop(i)
                alpha_sort_i += 1
                i = len(clades)
            else:
                pass
            i += 1

    if len(palette_taxa_map) != num_clades:
        logging.error("Some clades were not identified in the lineages provided!\n" +
                      "Clades = " + str(num_clades) + "\n" +
                      "Clades mapped = " + str(len(palette_taxa_map)) + "\n" +
                      str(clades) + "\n")
        raise AssertionError

    return palette_taxa_map


def write_colours_styles(colours, taxa_colours: TaxaColours):

    clades = set()
    lineages = set()
    for leaf in taxa_colours.tree_leaves:
        lineages.add(leaf.lineage)
    for num in taxa_colours.unique_clades:
        clades.add(taxa_colours.unique_clades[num])

    palette_taxa_map = map_colours_to_taxa(list(clades), list(lineages), colours)

    colours_style_header = "TREE_COLORS\nSEPARATOR TAB\nDATA\n"

    colourless = set()

    colours_style_string = ""
    for leaf in taxa_colours.tree_leaves:
        coloured = False
        colours_style_line = list()
        for taxon in palette_taxa_map:
            if re.search(taxon, leaf.lineage):
                col = palette_taxa_map[taxon]
                coloured = True
                colours_style_line = [str(leaf.number), "range", col, taxon]
                break
        if len(colours_style_line) > 0:
            colours_style_string += "\t".join(colours_style_line) + "\n"
        if coloured is False:
            colourless.add(leaf.lineage)
    try:
        cs_handler = open(taxa_colours.output, 'w')
    except IOError:
        logging.error("Unable to open " + taxa_colours.output + " for writing!\n")
        raise IOError

    logging.info("Output is available in " + taxa_colours.output + ".\n")
    logging.debug(str(len(colourless)) + " lineages were not assigned to a colour:\n\t" +
                  "\n\t".join(colourless) + "\n")

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
        if c in [',', ')']:
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
    # TODO: Sort the nodes by their internal node order
    colours = get_colours(args, taxa_colours, args.palette)
    write_colours_styles(colours, taxa_colours)


if __name__ == "__main__":
    main()
