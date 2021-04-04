#!/usr/bin/env python3

__author__ = 'Connor Morgan-Lang'


import sys
import os

import logging

from treesapp.classy import TreeSAPP
from treesapp.clade_annotation import CladeAnnotation


class Layerer(TreeSAPP):
    def __init__(self):
        super(Layerer, self).__init__("layer")
        self.stages = {}
        self.target_refpkgs = list()
        self.treesapp_output = ""
        self.layered_table = ""

    def check_arguments(self, args) -> None:
        """Check that the required files (e.g. jplace, classifications, annotation files) exist."""
        self.treesapp_output = args.output
        if self.treesapp_output[-1] != os.sep:
            self.treesapp_output += os.sep
        self.var_output_dir = self.treesapp_output + "intermediates" + os.sep
        self.final_output_dir = self.treesapp_output + "final_outputs" + os.sep

        if not os.path.isfile(self.final_output_dir + self.classification_tbl_name):
            logging.error("Could not find a classification file in " + self.final_output_dir + "\n")
            sys.exit(3)
        if args.refpkg_dir:
            self.refpkg_dir = args.refpkg_dir

        self.layered_table = self.final_output_dir + "layered_" + self.classification_tbl_name
        return


def identify_field_position(field_name: str, header_fields: list):
    x = 0
    for field in header_fields:
        if field == field_name:
            return x
        x += 1
    logging.error("Unable to find field name '" + field_name + "' in classifications.tsv header!\n")
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
            logging.error("Header in classifications.tsv is unexpected!\n")
            sys.exit(7)
        self.query_name = fields[query_pos]
        self.i_node = fields[node_pos]
        self.assignment_fields = fields
        return


def parse_marker_classification_table(marker_classification_file):
    """
    Function to read classifications.tsv and gather the relevant information for adding extra annotations
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


def map_queries_to_annotations(marker_tree_info: dict, master_dat: dict):
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


def annotate_internal_nodes(internal_node_map: dict, clade_annotations: list) -> (dict, set):
    """
    A function for mapping the clusters to all internal nodes of the tree.
    It also adds overlapping functional annotations for deep internal nodes and ensures all the leaves are annotated.

    :param internal_node_map: A dictionary mapping the internal nodes (keys) to the leaf nodes (values)
    :param clade_annotations: A list of CladeAnnotation instances from a single feature annotation type i.e.
    all of their 'feature' attributes should be the same.
    :return: A dictionary of the annotation (AKA group) as keys and internal nodes as values
    """
    annotated_clade_members = dict()
    annotation_clusters = dict()
    leaf_group_members = dict()
    leaves_in_clusters = set()

    if len(clade_annotations) == 0:
        logging.error("No clade annotations provided for layering.\n")
        raise AssertionError(17)

    for clade_annot in clade_annotations:  # type: CladeAnnotation
        annotation_internal_nodes = clade_annot.get_internal_nodes(internal_node_map)
        if len(annotation_internal_nodes) == 0:
            logging.error("Unable to match leaf node names to internal nodes for the clade annotation:\n"
                          "{}.".format(str(clade_annot.feature)))
            sys.exit(17)
        annotation_clusters.update({clade_annot.name: clade_annot.get_internal_nodes(internal_node_map)})

    # Create a dictionary to map the cluster name (e.g. Function, Activity, Class, etc) to all the leaf nodes
    for annotation in annotation_clusters:
        if annotation not in annotated_clade_members:
            annotated_clade_members[annotation] = set()
        if annotation not in leaf_group_members:
            leaf_group_members[annotation] = set()
        for i_node in annotation_clusters[annotation]:
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


def write_classification_table(table_name: str, field_order: dict, master_dat: dict) -> None:
    """
    Writes data in master_dat to a new tabular file with original and extra annotation information

    :param table_name: Path to the layered classification table to write to
    :param field_order: A dictionary mapping the order of the table's fields to the field name
    :param master_dat: A dictionary mapping reference package prefixes to the list of all ClassifiedSequence instances
    :return: None
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

    try:
        table_handler = open(table_name, 'w')
    except IOError:
        logging.error("Unable to open " + table_name + " for writing!\n")
        sys.exit(3)
    table_handler.write("\n".join(new_classification_lines) + "\n")
    table_handler.close()

    return
