#!/usr/bin/env python3

__author__ = 'Connor Morgan-Lang'

import sys
import os
import re
import argparse
import inspect

cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile(inspect.currentframe()))[0]))
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)
sys.path.insert(0, cmd_folder + os.sep + ".." + os.sep)
from classy import TreeLeafReference
from utilities import clean_lineage_string

rank_depth_map = {0: "Cellular organisms", 1: "Kingdom",
                  2: "Phylum", 3: "Class", 4: "Order",
                  5: "Family", 6: "Genus", 7: "Species",
                  8: "Strain"}


def get_arguments():
    parser = argparse.ArgumentParser(add_help=False)
    required_args = parser.add_argument_group("Required arguments")
    required_args.add_argument("-t", "--tax_ids_table",
                               help="tax_ids table for the TreeSAPP marker under investigation. "
                                    "Found in data/tree_data/tax_ids*.txt.",
                               required=True)

    optopt = parser.add_argument_group("Optional options")
    optopt.add_argument('-r', "--rank", default="Class", required=False,
                        help="The rank to generate unique colours for [ DEFAULT = 'Class' ]")
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
    args.treesapp = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + os.sep

    if re.match(".*tax_ids_(.*).txt$", args.tax_ids_table):
        args.marker = re.match(".*tax_ids_(.*).txt$", args.tax_ids_table).group(1)
        args.output = args.marker + "_colours_style.txt"
        if args.verbose:
            sys.stdout.write("\tGenerating colour palette for " + args.marker + ".\n")
    else:
        sys.stderr.write("ERROR: name of tax_ids file is formatted oddly.\n")
        sys.exit(2)

    if sys.version_info > (2, 9):
        args.py_version = 3
    else:
        args.py_version = 2

    palette_options = ['BrBG', 'PRGn', 'Set3', "Greys", "Paired"]
    if args.palette not in palette_options:
        sys.stderr.write("ERROR: " + args.palette + " is not supported. Please select one of the following options:\n")
        sys.stderr.write(', '.join(palette_options) + "\n")
        sys.exit(4)
    else:
        sys.stdout.write("\tUsing palette '" + args.palette + "' for styling.\n")
        sys.stdout.flush()

    if args.rank not in rank_depth_map.values():
        sys.stderr.write("ERROR: rank '" + args.rank + "' not accepted! Please choose one of the following:\n")
        i = 1
        while i < 8:
            sys.stderr.write(rank_depth_map[i] + ',')
            i += 1
        sys.stderr.write("\n")
        sys.exit(1)

    return args


def read_tax_ids_file(args, filename):
    tree_leaves = list()
    try:
        if args.py_version == 3:
            cog_tax_ids = open(filename, 'r', encoding='utf-8')
        else:
            cog_tax_ids = open(filename, 'r')
    except IOError:
        sys.exit('ERROR: Can\'t open ' + str(filename) + '!\n')

    for line in cog_tax_ids:
        line = line.strip()
        try:
            fields = line.split("\t")
        except ValueError:
            sys.stderr.write('ValueError: .split(\'\\t\') on ' + str(line) +
                             " generated " + str(len(line.split("\t"))) + " fields.")
            sys.exit(9)
        if len(fields) == 2:
            number, translation = fields
            lineage = ""
        elif len(fields) == 3:
            number, translation, lineage = fields
        else:
            sys.stderr.write("ValueError: Unexpected number of fields in " + filename +
                             ".\nInvoked .split(\'\\t\') on line " + str(line))
            raise ValueError
        leaf = TreeLeafReference(number, translation)
        if lineage:
            leaf.lineage = lineage
            leaf.complete = True
        tree_leaves.append(leaf)

    cog_tax_ids.close()
    return tree_leaves


def get_clades(args, tree_number_translations):
    clades = dict()
    lineages = list()
    target_depth = 0
    for depth in rank_depth_map:
        if rank_depth_map[depth] == args.rank:
            target_depth = depth
            break
    if target_depth == 0:
        sys.stderr.write("ERROR: rank '" + args.rank + "' not accepted!\n")
        sys.exit(3)

    for leaf in tree_number_translations:
        if leaf.lineage:
            lineages.append(clean_lineage_string(leaf.lineage))

    n_unique_clades = 0
    for lineage in sorted(lineages):
        lineage_path = lineage.split('; ')
        if len(lineage_path) > target_depth:
            if len(clades) == 0 and n_unique_clades == 0:
                clades[n_unique_clades] = lineage_path[target_depth]
            elif clades[n_unique_clades] != lineage_path[target_depth]:
                n_unique_clades += 1
                clades[n_unique_clades] = lineage_path[target_depth]
        else:
            pass
    n_unique_clades += 1

    if len(clades) <= 1:
        sys.stderr.write("ERROR: " + str(n_unique_clades) + " unique clade(s) identified at rank '" + args.rank + "'\n")
        sys.stderr.write("Consider changing rank to something more specific.\n"
                         "It may also be worth glancing at " + args.tax_ids_table + " to see if its okay.\n")
        sys.exit(4)
    else:
        if args.verbose:
            sys.stdout.write("\tIdentified " + str(n_unique_clades) + " unique clades at " + args.rank + " rank.\n")
            sys.stdout.flush()

    return clades


def get_colours(clades, palette):
    # TODO: find a decent API for colorbrewer so we are able to generate a palette on the fly

    brbg = ["#543005", "#643906", "#744207", "#844C09", "#93570E", "#A16518", "#B07323",
            "#BF812C", "#C89343", "#D1A65A", "#DAB871", "#E2C786", "#E8D29A", "#EFDDAE",
            "#F5E7C2", "#F5EBD1", "#F5EFDF", "#F5F3ED", "#EEF3F2", "#E1F0EE", "#D4EDE9",
            "#C7EAE5", "#B2E1DA", "#9ED9D0", "#8AD1C6", "#75C5B9", "#5FB5AB", "#4AA69D",
            "#35978F", "#268981", "#177B73", "#086D65", "#006057", "#00544A", "#00483D",
            "#003C30"]
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
        sys.stderr.write("ERROR: Number of clades (" + str(n_clade) + ") exceeds the number of colours in " +
                         palette + "(" + str(n_brewed) + ")!\n")
        sys.stderr.write("\tUnfortunately we have hard-coded palettes and they do not scale.\n"
                         "\tWe are hoping to integrate a ColorBrewer API in the near future.\n"
                         "\tIn the meantime, would you please select a higher clade? Thank you!\n")
        sys.exit(7)
    else:
        colours = list()
        n = int(n_brewed / n_clade)
        i = n
        p = 0
        while p < n_brewed and len(colours) < n_clade:
            if i == n:
                colours.append(brewed[p])
                i = 0
            i += 1
            p += 1

    if len(colours) != n_clade:
        sys.stderr.write("ERROR: Bad colour parsing! len(colours) != number of clades.\n")
        sys.exit(8)

    return colours


def map_colours_to_taxa(clades, colours):
    palette_taxa_map = dict()
    for level in sorted(clades, key=int):
        taxon = clades[level]
        palette_taxa_map[taxon] = colours[level]
    return palette_taxa_map


def write_colours_styles(args, colours, clades, tree_number_translations):
    palette_taxa_map = map_colours_to_taxa(clades, colours)

    colours_style_header = "TREE_COLORS\nSEPARATOR TAB\nDATA\n"
    colours_style_string = ""
    for leaf in tree_number_translations:
        colours_style_line = list()
        for taxon in palette_taxa_map:
            if re.search(taxon, leaf.lineage):
                col = palette_taxa_map[taxon]
                # colours_style_line = [str(leaf.number), "label", col, "bold", str(2)]
                colours_style_line = [str(leaf.number), "range", col, taxon]
                break
        if len(colours_style_line) > 0:
            colours_style_string += "\t".join(colours_style_line) + "\n"

    try:
        cs_handler = open(args.output, 'w')
    except IOError:
        sys.stderr.write("ERROR: Unable to open " + args.output + " for writing!\n")
        raise IOError

    if args.verbose:
        sys.stdout.write("\tOutput is available in " + args.output + ".\n")

    cs_handler.write(colours_style_header)
    cs_handler.write(colours_style_string)
    cs_handler.close()

    return


def main():
    args = get_arguments()
    tree_number_translations = read_tax_ids_file(args, args.tax_ids_table)
    clades = get_clades(args, tree_number_translations)
    colours = get_colours(clades, args.palette)
    write_colours_styles(args, colours, clades, tree_number_translations)


if __name__ == "__main__":
    main()
