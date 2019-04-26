#!/usr/bin/env python3

__author__ = 'Connor Morgan-Lang'


import sys
import logging
import os
import re
from .jplace_utils import jplace_parser
from .classy import TreeProtein, Layerer


def check_arguments(layerer: Layerer, args):
    """
    Check that the required files (e.g. jplace, marker_contig_map, annotation files) exist
    :param args:
    :return:
    """
    layerer.treesapp_output = args.output
    layerer.var_output_dir = layerer.treesapp_output + "intermediates" + os.sep
    layerer.final_output_dir = layerer.treesapp_output + "final_outputs" + os.sep
    if not os.path.isfile(layerer.final_output_dir + "marker_contig_map.tsv"):
        logging.error("Could not find a classification file in " + layerer.final_output_dir + "\n")
        sys.exit(3)
    for annot_f in args.colours_style:
        if not os.path.isfile(annot_f):
            logging.error(annot_f + " does not exist!\n")
            sys.exit(3)
    return


def parse_marker_classification_table(marker_classification_file):
    """
    Function to read marker_contig_map.tsv and gather the relevant information for adding extra annotations
    This function is different from Clade_exclusion_analyzer::read_marker_classification_table(assignment_file)
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

    header = classifications.readline()
    header_fields = header.strip().split("\t")
    x = 0
    query_index = None
    for field in header_fields:
        field_order[x] = field
        if field == "Query":
            query_index = x
        x += 1
    line = classifications.readline()
    while line:
        fields = line.strip().split("\t")
        if len(fields) != len(header_fields):
            logging.error("Inconsistent number of columns in table! Offending line:\n" + line + "\n")
            sys.exit(3)
        master_dat[fields[query_index]] = dict()
        y = 1
        while y < len(fields):
            # TODO: format the master_dat as a dictionary indexed by the refpkg containing a list of TreeProteins
            master_dat[fields[query_index]][header_fields[y]] = fields[y]
            y += 1
        line = classifications.readline()

    classifications.close()
    # for query in master_dat:
    #     print(master_dat[query])

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


def map_queries_to_annotations(marker_tree_info, marker_build_dict, master_dat):
    num_unclassified = 0
    for jplace in jplace_files_to_parse:
        file_name = os.path.basename(jplace)
        jplace_info = re.match(r"RAxML_portableTree.([A-Z][0-9]{4})_.*.jplace", file_name)
        refpkg_code = jplace_info.group(1)
        if refpkg_code not in master_dat.keys():
            num_unclassified += 1
            continue
        marker = marker_build_dict[refpkg_code].cog
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
                master_dat[refpkg_code][data_type] = ';'.join(sorted(metadata_placement))
    if num_unclassified > 0:
        logging.warning("Number of placed sequences that were unclassified: " + str(num_unclassified) + "\n")
    return master_dat


def annotate_internal_nodes(internal_node_map, clusters):
    """
    A function for mapping the clusters to all internal nodes of the tree.
    It also adds overlapping functional annotations for deep internal nodes and ensures all the leaves are annotated.
    :param internal_node_map: A dictionary mapping the internal nodes (keys) to the leaf nodes (values)
    :param clusters: Dictionary with the cluster names for keys and a tuple containing leaf boundaries as values
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
        for i_node_range in clusters[annotation]:
            for i_node in i_node_range:
                try:
                    for leaf in internal_node_map[int(i_node)]:
                        leaf_group_members[annotation].add(leaf)
                        leaves_in_clusters.add(leaf)
                except ValueError:
                    # TODO: Convert headers to internal nodes where an annotation cluster is a single leaf
                    pass
                except KeyError:
                    logging.error("Unable to find internal node " + i_node + " in internal node map.\n")
                    sys.exit(7)
        # Find the set of internal nodes that are children of this annotated clade
        for i_node in internal_node_map:
            if leaf_group_members[annotation].issuperset(internal_node_map[i_node]):
                annotated_clade_members[annotation].add(i_node)

    logging.debug("\tCaptured " + str(len(leaves_in_clusters)) + " nodes in clusters.\n")

    return annotated_clade_members, leaves_in_clusters


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
        logging.error("Unable to open " + output_file + " for writing!\n")
        sys.exit(3)

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
