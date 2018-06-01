#!/usr/bin/env python3

__author__ = 'Connor Morgan-Lang'


import sys
import argparse
import os
import re
import inspect
import glob
cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile(inspect.currentframe()))[0]))
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)
sys.path.insert(0, cmd_folder + os.sep + ".." + os.sep)
from entish import get_node
from treesapp import parse_ref_build_params, jplace_parser, tax_ids_file_to_leaves
from classy import TreeProtein
from utilities import read_colours_file, annotate_internal_nodes, convert_outer_to_inner_nodes


def get_arguments():
    parser = argparse.ArgumentParser(add_help=False,
                                     description="This script is generally used for layering extra annotations "
                                                 "beyond taxonomy (such as Subgroup or Metabolism) to TreeSAPP outputs."
                                                 " This is accomplished by adding an extra column (to all rows) of an "
                                                 "existing marker_contig_map.tsv and annotating the relevant sequences")
    required_args = parser.add_argument_group("Required arguments")
    required_args.add_argument("-c", "--colours_style",
                               help="The colours_style file exported from iTOL with the annotation information. "
                                    "For the variable name to be automatically inferred (rather than through `names`). "
                                    "Format of the file should be `marker`_`var`.txt. For example: mcrA_Metabolism.txt "
                                    "would create a new column in marker_contig_map.tsv named 'Metabolism'.",
                               required=True,
                               nargs='+')
    required_args.add_argument("-o", "--output",
                               help="The TreeSAPP output directory.",
                               required=True)

    # optopt = parser.add_argument_group("Optional options")
    # optopt.add_argument('-n', '--names',
    #                     help='The names corresponding to each of the colours_style files.'
    #                          ' Provide a comma-separated list if multiple colours_style files.',
    #                     required=False,
    #                     default=None,
    #                     nargs='+')

    miscellaneous_opts = parser.add_argument_group("Miscellaneous options")
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


def check_arguments(args):
    """
    Check that the required files (e.g. jplace, marker_contig_map, annotation files) exist
    :param args:
    :return:
    """
    various_dir = args.output + os.sep + "various_outputs" + os.sep
    jplace_files = glob.glob(various_dir + "*jplace")
    if len(jplace_files) == 0:
        sys.stderr.write("ERROR: could not find .jplace files in " + various_dir + "\n")
        sys.exit()
    if not os.path.isfile(os.sep.join([args.output, "final_outputs", "marker_contig_map.tsv"])):
        sys.stderr.write("ERROR: could not find a classification file in " + args.output + os.sep + "final_outputs\n")
        sys.exit()
    for annot_f in args.colours_style:
        if not os.path.isfile(annot_f):
            sys.stderr.write("ERROR: " + annot_f + " does not exist!\n")
            sys.exit()
    return jplace_files


def parse_marker_classification_table(args):
    """
    Function to read marker_contig_map.tsv and gather the relevant information for adding extra annotations
    This function is different from Clade_exclusion_analyzer::read_marker_classification_table(assignment_file)
    as we are interested in all fields in this function.
    :param args:
    :return:
    """
    master_dat = dict()
    field_order = dict()
    marker_classification_file = os.sep.join([args.output, "final_outputs", "marker_contig_map.tsv"])
    try:
        classifications = open(marker_classification_file, 'r')
    except IOError:
        sys.stderr.write("ERROR: Unable to open " + marker_classification_file + " for reading!\n")
        sys.exit()

    header = classifications.readline()
    header_fields = header.strip().split("\t")
    x = 0
    for field in header_fields:
        field_order[x] = field
        x += 1
    line = classifications.readline()
    while line:
        fields = line.strip().split("\t")
        if len(fields) != len(header_fields):
            sys.stderr.write("ERROR: Inconsistent number of columns in table! Offending line:\n" + line + "\n")
            sys.exit()
        master_dat[fields[0]] = dict()
        y = 1
        while y < len(fields):
            master_dat[fields[0]][header_fields[y]] = fields[y]
            y += 1
        line = classifications.readline()

    classifications.close()
    # for query in master_dat:
    #     print(master_dat[query])

    if args.verbose:
        sys.stdout.write("\tClassifications read:\t" + str(len(master_dat.keys())) + "\n")

    return master_dat, field_order


def names_for_nodes(clusters, node_map, taxa_map):
    """
    This function is used to convert from a string name of a leaf (e.g. Methylocapsa_acidiphila_|_CAJ01617)
    to an internal node number when all other nodes are internal nodes. Because consistent parsing is preferred!
    :param clusters:
    :param node_map:
    :param taxa_map:
    :return:
    """
    node_only_clusters = dict()
    for annotation in clusters:
        node_only_clusters[annotation] = list()
        for inodes in clusters[annotation]:
            node_1, node_2 = inodes
            try:
                int(node_1)
            except ValueError:
                # print(node_1)
                for leaf in taxa_map:
                    if re.sub(' ', '_', leaf.description) == node_1:
                        leaf_node = leaf.number
                        for inode_key, clade_value in node_map.items():
                            if clade_value[0] == leaf_node:
                                node_1 = inode_key
                                break
                        break
                    else:
                        pass
                # print(node_1)
            node_only_clusters[annotation].append((node_1, node_2))
    return node_only_clusters


def create_node_map(jplace_tree_string):
    """
    Loads a mapping between all nodes (internal and leaves) and all leaves
    :return:
    """
    no_length_tree = re.sub(":[0-9.]+{", ":{", jplace_tree_string)
    node_map = dict()
    node_stack = list()
    leaf_stack = list()
    x = 0
    num_buffer = ""
    while x < len(no_length_tree):
        c = no_length_tree[x]
        if re.search(r"[0-9]", c):
            while re.search(r"[0-9]", c):
                num_buffer += c
                x += 1
                c = no_length_tree[x]
            node_stack.append([str(num_buffer)])
            num_buffer = ""
            x -= 1
        elif c == ':':
            # Append the most recent leaf
            current_node, x = get_node(no_length_tree, x + 1)
            node_map[current_node] = node_stack.pop()
            leaf_stack.append(current_node)
        elif c == ')':
            # Set the child leaves to the leaves of the current node's two children
            while c == ')' and x < len(no_length_tree):
                if no_length_tree[x + 1] == ';':
                    break
                current_node, x = get_node(no_length_tree, x + 2)
                node_map[current_node] = node_map[leaf_stack.pop()] + node_map[leaf_stack.pop()]
                leaf_stack.append(current_node)
                x += 1
                c = no_length_tree[x]
        x += 1
    return node_map


def map_queries_to_annotations(marker_tree_info, marker_build_dict, jplace_files_to_parse, master_dat):
    num_unclassified = 0
    for jplace in jplace_files_to_parse:
        file_name = os.path.basename(jplace)
        jplace_info = re.match("RAxML_portableTree.([A-Z][0-9]{4})_(.*).jplace", file_name)
        gene_code, contig_name = jplace_info.groups()
        if contig_name not in master_dat.keys():
            num_unclassified += 1
            continue
        marker = marker_build_dict[gene_code].cog
        for data_type in marker_tree_info:
            if marker in marker_tree_info[data_type]:
                metadata_placement = set()
                query_obj = TreeProtein()
                itol_datum = jplace_parser(jplace)
                query_obj.transfer(itol_datum)
                query_obj.correct_decoding()
                query_obj.filter_max_weight_placement()
                for node in query_obj.list_placements():
                    for group in marker_tree_info[data_type][marker]:
                        if int(node) in marker_tree_info[data_type][marker][group]:
                            metadata_placement.add(group)

                if len(metadata_placement) == 0:
                    metadata_placement.add("Unknown")
                master_dat[contig_name][data_type] = ';'.join(sorted(metadata_placement))
    if num_unclassified > 0:
        sys.stdout.write("Number of placed sequences that were unclassified: " + str(num_unclassified) + "\n")
    return master_dat


def write_classification_table(args, field_order, master_dat):
    """
    Writes data in master_dat to a new tabular file with original and extra annotation information
    :param args:
    :param field_order:
    :param master_dat:
    :return:
    """
    output_file = os.sep.join([args.output, "final_outputs", "extra_annotated_marker_contig_map.tsv"])
    try:
        table_handler = open(output_file, 'w')
    except IOError:
        sys.stderr.write("ERROR: Unable to open " + output_file + " for writing!\n")
        sys.exit()

    fields = list()
    # Prepare the new header and write it to the new classification table
    for order in sorted(field_order.keys()):
        fields.append(field_order[order])
    table_handler.write("\t".join(fields) + "\n")

    # Now parse the classification data in master_dat
    for query_name in master_dat:
        query_classification_list = [query_name]
        field_acc = 1
        while field_acc < len(fields):
            query_classification_list.append(str(master_dat[query_name][field_order[field_acc]]))
            field_acc += 1
        table_handler.write("\t".join(query_classification_list) + "\n")

    table_handler.close()

    return


def main():
    ##
    # Worklow:
    #   1. Read data/tree_data/ref_build_parameters.tsv to get marker codes, denominators, and more (oh my!)
    #   2. Read the marker_contig_map.tsv file from the output directory to create the master data structure
    #   3. For each of the colours_styles files provided (potentially multiple for the same marker):
    #       3.1) Add the annotation variable to master_dat for every sequence (instantiate with "NA")
    #       3.2) Read the .jplace file for every sequence classified as marker
    #       3.3) Add the annotation information to every sequence classified as marker in master_dat
    #   4. Write the new classification file called "extra_annotated_marker_contig_map.tsv"
    ##
    args = get_arguments()
    jplace_files = check_arguments(args)
    var_dir = args.output + os.sep + "various_outputs" + os.sep

    marker_subgroups = dict()
    unique_markers_annotated = set()
    marker_tree_info = dict()
    jplace_files_to_parse = list()
    jplace_tree_strings = dict()
    internal_nodes = dict()
    marker_build_dict = parse_ref_build_params(args)
    master_dat, field_order = parse_marker_classification_table(args)
    # structure of master dat:
    # {"Sequence_1": {"Field1": x, "Field2": y, "Extra": n},
    #  "Sequence_2": {"Field1": i, "Field2": j, "Extra": n}}
    for annot_f in args.colours_style:
        # Determine the marker being annotated
        marker = data_type = ""
        for gene_code in marker_build_dict:
            marker = marker_build_dict[gene_code].cog
            if re.search("^" + re.escape(marker) + "_\w+.txt$", os.path.basename(annot_f)):
                data_type = re.search("^" + re.escape(marker) + "_(\w+).txt$", os.path.basename(annot_f)).group(1)
                unique_markers_annotated.add(gene_code)
                break
            else:
                marker = data_type = ""
        if marker and data_type:
            if data_type not in marker_subgroups:
                marker_subgroups[data_type] = dict()
                internal_nodes[data_type] = dict()
            marker_subgroups[data_type][marker], internal_nodes[data_type][marker] = read_colours_file(args, annot_f)
        else:
            sys.stderr.write("ERROR: Unable to parse the marker and/or annotation type from " + annot_f + "\n")
            sys.exit()
    # Instantiate every query sequence in marker_contig_map with an empty string for each data_type
    for data_type in marker_subgroups:
        for query_seq in master_dat:
            master_dat[query_seq][data_type] = "NA"
    # Update the field_order dictionary with new fields
    field_acc = len(field_order)
    for new_datum in sorted(marker_subgroups.keys()):
        field_order[field_acc] = new_datum
        field_acc += 1

    # Load the query sequence annotations
    for data_type in marker_subgroups:
        if data_type not in marker_tree_info:
            marker_tree_info[data_type] = dict()
        for gene_code in unique_markers_annotated:
            for jplace in jplace_files:
                if re.match(re.escape(var_dir) + "RAxML_portableTree." + re.escape(gene_code) + "_.*.jplace", jplace):
                    jplace_files_to_parse.append(jplace)
                    if gene_code not in jplace_tree_strings:
                        jplace_tree_strings[gene_code] = jplace_parser(jplace).tree

            marker = marker_build_dict[gene_code].cog
            if marker in marker_subgroups[data_type]:
                # Create the dictionary mapping an internal node to all child nodes
                internal_node_map = create_node_map(jplace_tree_strings[gene_code])

                # Routine for exchanging any organism designations for their respective node number
                tax_ids_file = os.sep.join([args.treesapp, "data", "tree_data", "tax_ids_" + marker + ".txt"])
                taxa_map = tax_ids_file_to_leaves(args, tax_ids_file)
                clusters = names_for_nodes(marker_subgroups[data_type][marker], internal_node_map, taxa_map)

                if not internal_nodes[data_type][marker]:
                    # Convert the leaf node ranges to internal nodes for consistency
                    clusters = convert_outer_to_inner_nodes(clusters, internal_node_map)

                marker_tree_info[data_type][marker], leaves_in_clusters = annotate_internal_nodes(args,
                                                                                                  internal_node_map,
                                                                                                  clusters)
                diff = len(taxa_map) - len(leaves_in_clusters)
                if diff != 0:
                    unannotated = set()
                    sys.stderr.write("WARNING: the following leaf nodes were not mapped to annotation groups:\n")
                    for inode in internal_node_map:
                        for leaf in internal_node_map[inode]:
                            if leaf not in leaves_in_clusters:
                                unannotated.add(str(leaf))
                    sys.stderr.write("\t" + ', '.join(sorted(unannotated, key=int)) + "\n")
                    sys.stderr.flush()
            else:
                pass
    marker_subgroups.clear()
    master_dat = map_queries_to_annotations(marker_tree_info, marker_build_dict, jplace_files_to_parse, master_dat)
    write_classification_table(args, field_order, master_dat)


main()
