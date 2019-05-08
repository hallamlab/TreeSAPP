#!/usr/bin/env python3

__author__ = 'Connor Morgan-Lang'

import sys
import os
import re
import argparse
import logging
from treesapp.classy import TreeLeafReference, prep_logging
from treesapp.utilities import clean_lineage_string

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
        self.unique_clades = dict()

    def get_clades(self, target_depth):
        clades = dict()
        lineages = list()
        
        for leaf in self.tree_leaves:
            if leaf.lineage:
                lineages.append(clean_lineage_string(leaf.lineage))

        self.num_taxa = 0
        for lineage in sorted(lineages):
            lineage_path = lineage.split('; ')
            if len(lineage_path) > target_depth:
                if len(clades) == 0 and self.num_taxa == 0:
                    clades[self.num_taxa] = lineage_path[target_depth]
                elif clades[self.num_taxa] != lineage_path[target_depth]:
                    self.num_taxa += 1
                    clades[self.num_taxa] = lineage_path[target_depth]
            else:
                pass
        self.num_taxa += 1
        self.unique_clades = clades

        logging.debug("\t" + self.marker + ": " + str(self.num_taxa) + " unique clades.\n")

        return


def get_arguments():
    parser = argparse.ArgumentParser(add_help=False)
    required_args = parser.add_argument_group("Required arguments")
    required_args.add_argument("-t", "--tax_ids_table", nargs='+',
                               help="tax_ids table for the TreeSAPP marker under investigation. "
                                    "Found in data/tree_data/tax_ids*.txt.",
                               required=True)
    # TODO: Optionally not colour polyphyletic clades
    optopt = parser.add_argument_group("Optional options")
    optopt.add_argument('-r', "--rank", default="Order", required=False,
                        help="The rank to generate unique colours for [ DEFAULT = 'Order' ]")
    optopt.add_argument('-p', "--palette", default="BrBG",
                        help="The ColorBrewer palette to use. Currently supported:"
                             " 'BrBG', 'PRGn', 'Set3', 'Greys', 'Paired'. [ DEFAULT = BrBG ]")
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

    :param taxa_colours: A dictionary of TaxaColours instances
    :return: 
    """
    for marker in taxa_colours:
        ref = taxa_colours[marker]
        try:
            cog_tax_ids = open(ref.tax_ids_file, 'r', encoding='utf-8')
        except IOError:
            logging.error('Can\'t open ' + str(ref.tax_ids_file) + '!\n')
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
                logging.error("Unexpected number of fields in " + ref.tax_ids_file +
                              ".\nInvoked .split(\'\\t\') on line " + str(line))
                raise ValueError
            leaf = TreeLeafReference(number, translation)
            if lineage:
                leaf.lineage = lineage
                leaf.complete = True
            leaves.append(leaf)
        
        ref.tree_leaves = leaves
        cog_tax_ids.close()
            
    return taxa_colours


def get_colours(args, taxa_colours, palette):
    # TODO: find a decent API for colorbrewer so we are able to generate a palette on the fly

    clades = set()
    for marker in taxa_colours:
        ref = taxa_colours[marker]
        for num in ref.unique_clades:
            clades.add(ref.unique_clades[num])

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


def write_colours_styles(colours, taxa_colours):

    clades = set()
    lineages = set()
    for marker in taxa_colours:
        ref = taxa_colours[marker]
        for leaf in ref.tree_leaves:
            lineages.add(leaf.lineage)
        for num in ref.unique_clades:
            clades.add(ref.unique_clades[num])

    palette_taxa_map = map_colours_to_taxa(list(clades), list(lineages), colours)

    colours_style_header = "TREE_COLORS\nSEPARATOR TAB\nDATA\n"

    colourless = set()

    for marker in taxa_colours:
        colours_style_string = ""
        ref = taxa_colours[marker]
        for leaf in ref.tree_leaves:
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
            cs_handler = open(ref.output, 'w')
        except IOError:
            logging.error("Unable to open " + ref.output + " for writing!\n")
            raise IOError

        logging.info("Output is available in " + ref.output + ".\n")
        logging.debug(str(len(colourless)) + " lineages were not assigned to a colour:\n\t" +
                      "\n\t".join(colourless) + "\n")

        cs_handler.write(colours_style_header)
        cs_handler.write(colours_style_string)
        cs_handler.close()

    return


def init_taxa_colours(args):
    taxa_colours = dict()
    for tax_ids_table in args.tax_ids_table:
        if re.match(".*tax_ids_(.*).txt$", tax_ids_table):
            ref = TaxaColours(tax_ids_table)
            taxa_colours[ref.marker] = ref
            logging.debug("\tGenerating colour palette for " + ref.marker + ".\n")
        else:
            logging.error("Name of tax_ids file is formatted oddly.\n")
            sys.exit(2)

    return taxa_colours


def find_rank_depth(args):
    target_depth = 0
    for depth in rank_depth_map:
        if rank_depth_map[depth] == args.rank:
            target_depth = depth - 1
            break
    if target_depth == 0:
        logging.error("Rank '" + args.rank + "' not accepted!\n")
        sys.exit(3)
    return target_depth


def main():
    args = get_arguments()

    log_file_name = os.path.abspath("./TreeSAPP_auto-colour_log.txt")
    prep_logging(log_file_name, args.verbose)
    logging.info("\n##\t\t\tGenerating colour-style file for iTOL\t\t\t##\n")
    validate_command(args)

    taxa_colours = init_taxa_colours(args)
    taxa_colours = read_tax_ids_file(taxa_colours)
    target_depth = find_rank_depth(args)
    for marker in sorted(taxa_colours):
        # TODO: Sort the nodes by their internal node order
        ref = taxa_colours[marker]
        ref.get_clades(target_depth)
    colours = get_colours(args, taxa_colours, args.palette)
    write_colours_styles(colours, taxa_colours)


if __name__ == "__main__":
    main()
