#!/usr/bin/env python3

__author__ = 'Connor Morgan-Lang'


import sys
import argparse
import os
import re
from .entish import read_and_map_internal_nodes_from_newick_tree
from .fasta import format_read_fasta
from phylo_seq import TreeLeafReference
from .file_parsers import read_colours_file, parse_ref_build_params


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
    required_args.add_argument("--name",
                               help="Name of the marker, either code (e.g. D0601) or short form (e.g. nifH)",
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

    if not os.path.isdir(args.output):
        os.makedirs(args.output)

    return args


def parse_clades_from_tree(args, clusters, tree_file, marker):
    clade_members = dict()
    leaves_in_clusters = 0

    internal_node_map = read_and_map_internal_nodes_from_newick_tree(tree_file, marker)

    for cluster in clusters.keys():
        clade_members[cluster] = list()
        for frond_tips in clusters[cluster]:
            start, end = frond_tips
            clade_members[cluster] = internal_node_map[int(start)]
            leaves_in_clusters += len(internal_node_map[int(start)])
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
    for cluster in clade_members:
        output_fasta = args.output + os.sep + args.prefix + '_' + str(re.sub(' ', '_', cluster)) + ".fasta"
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
                    header_string = '>' + accession + ' ' + cluster + " [" +\
                                    ref_tree_leaf.description.split(" | ")[0] + ']'
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
    marker_build_dict = parse_ref_build_params(args.treesapp)
    # Find the marker name and the code name for the name inputted
    if re.match("[A-Z][0-9]{4}", args.name):
        args.code_name = args.name
        args.marker = marker_build_dict[args.code_name].cog
    elif len(args.name) <= 6:
        # Assuming the provided name is a short gene name (e.g. mcrA, nifH)
        args.marker = args.name
        args.code_name = ""
        for denominator in marker_build_dict:
            if marker_build_dict[denominator].cog == args.marker:
                args.code_name = denominator
                break
        if not args.code_name:
            sys.stderr.write("ERROR: Unable to identify the gene name from the code name '" + args.code_name + "'.\n")
            sys.stderr.flush()
            sys.exit()
    else:
        sys.stderr.write("ERROR: Wrong format for the reference code_name provided: " + args.name + "\n")
        sys.exit(9)

    args.fasta = os.sep.join([args.treesapp, "data", "alignment_data", args.marker + ".fa"])
    args.tax_ids =os.sep.join([args.treesapp, "data", "tree_data", "tax_ids_" + args.marker + ".txt"])

    clusters, internal_nodes = read_colours_file(args.colours_style)
    if internal_nodes:
        # The tree file needs to be read with internal nodes as keys
        pass

    tree_file = os.sep.join([args.treesapp, "data", "tree_data", args.marker + "_tree.txt"])
    clade_members = parse_clades_from_tree(args, clusters, tree_file, args.marker)

    fasta_dict = format_read_fasta(args.fasta, args.molecule)
    fasta_dict = strip_header_prefix(fasta_dict)
    fasta_dict = remove_dashes_from_msa_dict(fasta_dict)
    tax_ids = read_tax_ids_table(args)
    write_cluster_fasta(args, clade_members, fasta_dict, tax_ids)


main()
