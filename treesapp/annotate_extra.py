#!/usr/bin/env python3

__author__ = 'Connor Morgan-Lang'


import sys
import os
import re
import glob

import logging

from treesapp.classy import Layerer


def check_arguments(layerer: Layerer, args):
    """
    Check that the required files (e.g. jplace, marker_contig_map, annotation files) exist

    :param layerer:
    :param args:
    :return:
    """
    layerer.treesapp_output = args.output
    if layerer.treesapp_output[-1] != os.sep:
        layerer.treesapp_output += os.sep
    layerer.var_output_dir = layerer.treesapp_output + "intermediates" + os.sep
    layerer.final_output_dir = layerer.treesapp_output + "final_outputs" + os.sep
    if not os.path.isfile(layerer.final_output_dir + "marker_contig_map.tsv"):
        logging.error("Could not find a classification file in " + layerer.final_output_dir + "\n")
        sys.exit(3)
    if args.colours_style:
        for annot_f in args.colours_style:
            if not os.path.isfile(annot_f):
                logging.error(annot_f + " does not exist!\n")
                sys.exit(3)
            layerer.annot_files.append(annot_f)
    else:
        # If a directory containing annotation files isn't given, set it to the default data/iTOL_data directory
        if args.annot_dir is None:
            args.annot_dir = layerer.itol_dir
        annotation_files = glob.glob(args.annot_dir + '*')
        # Add all files in the annot_dir to the colours_style list
        for af in annotation_files:
            if not layerer.c_strip_re.match(af) and not layerer.c_style_re.match(af):
                layerer.annot_files.append(af)
    return


def identify_field_position(field_name: str, header_fields: list):
    x = 0
    for field in header_fields:
        if field == field_name:
            return x
        x += 1
    logging.error("Unable to find field name '" + field_name + "' in marker_contig_map.tsv header!\n")
    sys.exit()


class ClassifiedSequence:
    def __init__(self, assigned_refpkg):
        self.refpkg = assigned_refpkg
        self.query_name = ""
        self.i_node = ""
        self.assignment_fields = None
        self.expected_header = ['Sample', 'Query', 'Marker', 'Start_pos', 'End_pos', 'Taxonomy',
                                'Abundance', 'iNode', 'E-value', 'LWR', 'EvoDist', 'Distances']
        self.layers = dict()
        return

    def load_assignment_line(self, fields, header_fields, query_pos, node_pos):
        if header_fields != self.expected_header:
            logging.error("Header in marker_contig_map.tsv is unexpected!\n")
            sys.exit(7)
        self.query_name = fields[query_pos]
        self.i_node = fields[node_pos]
        self.assignment_fields = fields
        return


def parse_marker_classification_table(marker_classification_file):
    """
    Function to read marker_contig_map.tsv and gather the relevant information for adding extra annotations
    This function is different from Clade_exclusion_analyzer::read_classification_table(assignment_file)
    as we are interested in all fields in this function.
    :param marker_classification_file:
    :return:
    """
    master_dat = dict()
    field_order = dict()
    try:
        classifications = open(marker_classification_file, 'r')
    except IOError:
        logging.error("Unable to open " + marker_classification_file + " for reading!\n")
        sys.exit(3)

    header_fields = classifications.readline().strip().split("\t")
    x = 0
    for field in header_fields:
        field_order[x] = field
        x += 1
    marker_pos = identify_field_position("Marker", header_fields)
    node_pos = identify_field_position("iNode", header_fields)
    query_pos = identify_field_position("Query", header_fields)

    line = classifications.readline()
    while line:
        fields = line.strip().split("\t")
        if len(fields) != len(header_fields):
            logging.error("Inconsistent number of columns in table! Offending line:\n" + line + "\n")
            sys.exit(3)
        if fields[marker_pos] not in master_dat:
            master_dat[fields[marker_pos]] = list()
        jplace_seq = ClassifiedSequence(fields[marker_pos])
        jplace_seq.load_assignment_line(fields, header_fields, query_pos, node_pos)
        master_dat[fields[marker_pos]].append(jplace_seq)
        line = classifications.readline()

    classifications.close()

    return master_dat, field_order


def names_for_nodes(clusters: dict, node_map: dict, taxa_map: list) -> dict:
    """
    This function is used to convert from a string name of a leaf (e.g. Methylocapsa_acidiphila_|_CAJ01617)
    to an internal node number when all other nodes are internal nodes. Because consistent parsing is preferred!

    :param clusters: A dictionary of lists where each list is populated by tuples with start and end leaves
    :param node_map: Dictionary of all internal nodes (keys) and a list of child leaves (values)
    :param taxa_map: List of TreeLeafReference instances parsed from tax_ids file
    :return:
    """
    node_only_clusters = dict()
    for annotation in clusters:
        node_only_clusters[annotation] = list()
        for inodes in clusters[annotation]:
            node_1, node_2 = inodes
            try:
                int(node_1)  # This is an internal node
            except ValueError:
                for leaf in taxa_map:
                    if re.sub(' ', '_', leaf.description) == node_1:
                        for inode_key, clade_value in node_map.items():
                            if len(clade_value) == 1:
                                leaf_node = clade_value[0]
                                if int(leaf_node.split('_')[0]) == int(leaf.number):
                                    node_1, node_2 = inode_key, inode_key
                                    break
                        continue
                    else:
                        pass
            node_only_clusters[annotation].append((node_1, node_2))
    return node_only_clusters


def map_queries_to_annotations(marker_tree_info, master_dat):
    """

    :param marker_tree_info:
    :param master_dat:
    :return:
    """
    num_unclassified = 0
    metadata_placement = set()
    for data_type in marker_tree_info:
        for marker in master_dat:
            if marker in marker_tree_info[data_type]:
                for query_obj in master_dat[marker]:  # type: ClassifiedSequence
                    for group in marker_tree_info[data_type][marker]:
                        if int(query_obj.i_node) in marker_tree_info[data_type][marker][group]:
                            metadata_placement.add(group)

                    if len(metadata_placement) == 0:
                        metadata_placement.add("Unknown")
                        num_unclassified += 1
                    query_obj.layers[data_type] = ';'.join(sorted(metadata_placement))
                    metadata_placement.clear()
            else:
                num_unclassified += len(master_dat[marker])
                continue
    if num_unclassified > 0:
        logging.debug("Number of placed sequences that were unclassified: " + str(num_unclassified) + "\n")
    return master_dat


def annotate_internal_nodes(internal_node_map: dict, clusters: dict) -> (dict, set):
    """
    A function for mapping the clusters to all internal nodes of the tree.
    It also adds overlapping functional annotations for deep internal nodes and ensures all the leaves are annotated.

    :param internal_node_map: A dictionary mapping the internal nodes (keys) to the leaf nodes (values)
    :param clusters: Dictionary with the cluster names for keys and a list of internal nodes as values
    :return: A dictionary of the annotation (AKA group) as keys and internal nodes as values
    """
    annotated_clade_members = dict()
    leaf_group_members = dict()
    leaves_in_clusters = set()

    # Create a dictionary to map the cluster name (e.g. Function, Activity, Class, etc) to all the leaf nodes
    for annotation in clusters:
        if annotation not in annotated_clade_members:
            annotated_clade_members[annotation] = set()
        if annotation not in leaf_group_members:
            leaf_group_members[annotation] = set()
        for i_node in clusters[annotation]:
            try:
                for leaf in internal_node_map[int(i_node)]:
                    leaf_group_members[annotation].add(leaf)
                    leaves_in_clusters.add(leaf)
            except ValueError:
                # TODO: Convert headers to internal nodes where an annotation cluster is a single leaf
                logging.warning("Unable to assign '{}' to an internal node ID.\n".format(i_node))
            except KeyError:
                logging.error("Unable to find internal node '{}' in internal node map.\n".format(i_node))
                sys.exit(7)
        # Find the set of internal nodes that are children of this annotated clade
        for i_node in internal_node_map:
            if leaf_group_members[annotation].issuperset(internal_node_map[i_node]):
                annotated_clade_members[annotation].add(i_node)

    logging.debug("\tCaptured {} nodes in clusters.\n".format(len(leaves_in_clusters)))

    return annotated_clade_members, leaves_in_clusters


def write_classification_table(output_dir, field_order, master_dat):
    """
    Writes data in master_dat to a new tabular file with original and extra annotation information

    :param output_dir:
    :param field_order:
    :param master_dat:
    :return:
    """
    fields = list()
    # Prepare the new header and write it to the new classification table
    for order in sorted(field_order.keys()):
        fields.append(field_order[order])
    new_classification_lines = ["\t".join(fields)]

    # Now parse the classification data in master_dat
    for refpkg_code in master_dat:
        for assignment in master_dat[refpkg_code]:  # type: ClassifiedSequence
            for order in sorted(field_order.keys()):
                field = field_order[order]
                if field in assignment.layers:
                    assignment.assignment_fields.append(assignment.layers[field])
            new_classification_lines.append("\t".join(assignment.assignment_fields))

    output_file = output_dir + "extra_annotated_marker_contig_map.tsv"
    try:
        table_handler = open(output_file, 'w')
    except IOError:
        logging.error("Unable to open " + output_file + " for writing!\n")
        sys.exit(3)
    table_handler.write("\n".join(new_classification_lines) + "\n")
    table_handler.close()

    return
