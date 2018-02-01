__author__ = 'Connor Morgan-Lang'


import sys
import os
import re
from classy import TreeLeafReference
import argparse

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
    optopt.add_argument('-p', "--palette", default="BrBG")
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
        print(args.marker)
    else:
        sys.stderr.write("ERROR: name of tax_ids file is formatted oddly.\n")
        sys.exit(2)

    if sys.version_info > (2, 9):
        args.py_version = 3
    else:
        args.py_version = 2

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


def get_clades(rank, tree_number_translations):
    clades = set()
    target_depth = 0
    for depth in rank_depth_map:
        if rank_depth_map[depth] == rank:
            target_depth = depth
            break
    if target_depth == 0:
        sys.stderr.write("ERROR: rank '" + rank + "' not accepted!\n")
        sys.exit(3)

    for leaf in tree_number_translations:
        if leaf.lineage:
            lineage_path = leaf.lineage.split('; ')
            if len(lineage_path) > target_depth:
                clades.add(lineage_path[target_depth])
            else:
                pass

    if len(clades) <= 1:
        sys.stderr.write("ERROR: " + str(len(clades)) + " unique clade(s) identified using rank '" + rank + "'\n")
        sys.stderr.write("Consider changing rank to something more specific.\n")
        sys.exit(4)

    return clades


def get_colours(clades, palette):
    colours = list()
    return colours


def write_colours_styles(output, colours):
    return


def main():
    args = get_arguments()
    tree_number_translations = read_tax_ids_file(args, args.tax_ids_table)
    clades = get_clades(args.rank, tree_number_translations)
    colours = get_colours(clades, args.palette)
    write_colours_styles(args.output, colours)


if __name__ == "__main__":
    main()
