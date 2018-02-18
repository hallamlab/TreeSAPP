#!/usr/bin/env python3

__author__ = 'Connor Morgan-Lang'


import sys
import argparse
import os
import re
import inspect
cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile(inspect.currentframe()))[0]))
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)
sys.path.insert(0, cmd_folder + os.sep + ".." + os.sep)
from entish import read_and_understand_the_reference_tree
from fasta import format_read_fasta
from classy import TreeLeafReference


def get_arguments():
    parser = argparse.ArgumentParser(add_help=False,
                                     description="Script for extracting all sequences in tree clades "
                                                 "based on the first and last leaves for each clade. "
                                                 "File with coordinates is an exported iTOL color.txt")
    required_args = parser.add_argument_group("Required arguments")
    required_args.add_argument("--colours_style",
                               help="The colours_style file exported from iTOL with the coordinates for each cluster",
                               required=True)
    required_args.add_argument("--prefix",
                               help="Prefix of the FASTA files to be exported",
                               required=True)
    required_args.add_argument("--tree",
                               help="Tree file in Newick format",
                               required=True)
    required_args.add_argument("--tax_ids",
                               help="tax_ids table for the TreeSAPP marker under investigation",
                               required=True)
    required_args.add_argument("--fasta",
                               help="FASTA file with the sequences used to create the tree. "
                                    "Sequences will be pulled from here and written into cluster-specific FASTAs.",
                               required=True)

    miscellaneous_opts = parser.add_argument_group("Miscellaneous options")
    miscellaneous_opts.add_argument("-m", "--molecule",
                                    help='the type of input sequences '
                                         '(prot = Protein [DEFAULT]; dna = Nucleotide; rrna = rRNA)',
                                    default='prot',
                                    choices=['prot', 'dna', 'rrna'])
    miscellaneous_opts.add_argument("-o", "--output",
                                    help="The output directory. [ DEFAULT = ./output/ ]",
                                    default="./output/")
    miscellaneous_opts.add_argument('-v', '--verbose',
                                    action='store_true',
                                    default=False,
                                    help='Prints a more verbose runtime log')
    miscellaneous_opts.add_argument("-h", "--help",
                                    action="help",
                                    help="Show this help message and exit")

    args = parser.parse_args()
    args.treesapp = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + os.sep

    # Adding a dummy value for format_read_fasta args namespace
    args.min_seq_length = 10

    if sys.version_info > (2, 9):
        args.py_version = 3
    else:
        args.py_version = 2

    return args


def read_colours_file(args):
    """

    :param args:
    :return: A dictionary of lists where each list is populated by tuples with start and end leaves
    """
    # Determine the number of clusters
    clusters = dict()
    try:
        style_handler = open(args.colours_style, 'r')
    except IOError:
        sys.stderr.write("ERROR: Unable to open " + args.colours_style + " for reading!\n")
        sys.exit()

    data_line = re.compile("^([0-9]+)\|([0-9]+)	.*	.*	(.*)$")
    line = style_handler.readline()
    # Skip the header
    while line and not data_line.match(line):
        line = style_handler.readline()
    # Begin parsing the data from 4 columns
    while line:
        style_data = data_line.match(line)
        if style_data:
            description = re.sub(' ', '_', style_data.group(3))
            if description not in clusters.keys():
                clusters[description] = list()
            clusters[description].append((style_data.group(1), style_data.group(2)))
        line = style_handler.readline()

    if args.verbose:
        sys.stdout.write("\tParsed " + str(len(clusters)) + " clusters from " + args.colours_style + "\n")

    style_handler.close()
    return clusters


def parse_clades_from_tree(args, clusters):
    clade_members = dict()
    terminal_children_clades = list()
    leaves_in_clusters = 0

    denominator, terminal_children = read_and_understand_the_reference_tree(args.tree, "D0601")

    for clade in terminal_children.keys():
        terminal_children_clades.append([x for x in clade.strip().split(' ')])
    for cluster in clusters.keys():
        clade_members[cluster] = list()
        for frond_tips in clusters[cluster]:
            start, end = frond_tips
            # Find the minimum set that includes both start and end
            warm_front = dict()
            for clade in terminal_children_clades:
                if start in clade:
                    warm_front[len(clade)] = clade
            for size in sorted(warm_front, key=int):
                if end in warm_front[size]:
                    leaves_in_clusters += len(warm_front[size])
                    clade_members[cluster] += warm_front[size]
                    break

    if args.verbose:
        sys.stdout.write("\tCaptured " + str(leaves_in_clusters) + " leaves in clusters.\n")

    return clade_members


def read_tax_ids_table(args):
    tree_numbers_translation = list()
    try:
        if args.py_version == 3:
            cog_tax_ids = open(args.tax_ids, 'r', encoding='utf-8')
        else:
            cog_tax_ids = open(args.tax_ids, 'r')
    except IOError:
        sys.exit('ERROR: Can\'t open ' + str(args.tax_ids) + '!\n')

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
            sys.stderr.write("ValueError: Unexpected number of fields in " + args.tax_ids +
                             ".\nInvoked .split(\'\\t\') on line " + str(line))
            raise ValueError
        leaf = TreeLeafReference(number, translation)
        if lineage:
            leaf.lineage = lineage
            leaf.complete = True
        tree_numbers_translation.append(leaf)

    cog_tax_ids.close()

    return tree_numbers_translation


def strip_header_prefix(fasta_dict):
    new_sequence_dict = dict()
    for header in fasta_dict:
        treesapp_header = re.match(">([0-9]+)_.*$", header)
        if not treesapp_header:
            sys.stderr.write("ERROR: Unexpected header format encountered: " + header + "\n")
            sys.exit()
        new_sequence_dict[treesapp_header.group(1)] = fasta_dict[header]

    if len(new_sequence_dict) != len(fasta_dict):
        sys.stderr.write("ERROR: Not all sequences were converted from `fasta_dict`!\n")
        sys.exit()

    return new_sequence_dict


def remove_dashes_from_msa_dict(fasta_dict):
    """

    :param fasta_dict:
    :return: stripped_fasta_dict, a dictionary with headers as keys and unaligned fasta sequences for values
    """
    stripped_fasta_dict = dict()
    for header in fasta_dict:
        stripped_fasta_dict[header] = re.sub('[-.]', '', fasta_dict[header])
    return stripped_fasta_dict


def write_cluster_fasta(args, clade_members, fasta_dict, tax_ids):
    if not os.path.isdir(args.output):
        os.makedirs(args.output)

    for cluster in clade_members:
        output_fasta = args.output + os.sep + args.prefix + str(cluster) + ".fasta"
        try:
            fasta_handler = open(output_fasta, 'w')
        except IOError:
            sys.stderr.write("ERROR: Unable to open " + output_fasta + " for writing!\n")
            raise IOError
        for leaf in clade_members[cluster]:
            header_string = ""
            for ref_tree_leaf in tax_ids:
                if ref_tree_leaf.number == leaf:
                    accession = ref_tree_leaf.description.split(" | ")[1]
                    header_string = '>' + accession + " [" + ref_tree_leaf.description.split(" | ")[0] + ']'
                    break
            if header_string:
                fasta_handler.write(header_string + "\n" + fasta_dict[leaf] + "\n")
            else:
                sys.stderr.write("ERROR: Unable to find a matching tree leaf for " +
                                 str(leaf) + " in " + args.tax_ids + "\n")
                sys.exit(5)

        fasta_handler.close()

    return


def main():
    args = get_arguments()
    clusters = read_colours_file(args)
    clade_members = parse_clades_from_tree(args, clusters)
    fasta_dict = format_read_fasta(args.fasta, args.molecule, args, 1000)
    fasta_dict = strip_header_prefix(fasta_dict)
    fasta_dict = remove_dashes_from_msa_dict(fasta_dict)
    tax_ids = read_tax_ids_table(args)
    write_cluster_fasta(args, clade_members, fasta_dict, tax_ids)


main()
