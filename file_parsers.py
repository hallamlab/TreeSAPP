#!/usr/bin/env python3

import sys
import os
import re
from classy import TreeLeafReference, MarkerBuild
from HMMER_domainTblParser import DomainTableParser, format_split_alignments, filter_incomplete_hits, filter_poor_hits

__author__ = 'Connor Morgan-Lang'


def parse_ref_build_params(args):
    """
    Returns a dictionary of MarkerBuild objects storing information pertaining to the build parameters of each marker.
    :param args: Command-line argument object returned by get_options and check_parser_arguments
    """
    ref_build_parameters = args.treesapp + 'data' + os.sep + 'tree_data' + os.sep + 'ref_build_parameters.tsv'
    try:
        param_handler = open(ref_build_parameters, 'r')
    except IOError:
        sys.exit('ERROR: Can\'t open ' + ref_build_parameters + '!\n')

    header_re = re.compile("code_name\tdenominator\taa_model\tcluster_identity\tlowest_confident_rank\tlast_updated$")
    if not header_re.match(param_handler.readline().strip()):
        raise AssertionError("ERROR: Header of '" + ref_build_parameters + "' is unexpected!")

    marker_build_dict = dict()
    for line in param_handler:
        if header_re.match(line):
            continue
        if line[0] == '#':
            continue
        line = line.strip()
        marker_build = MarkerBuild(line)
        marker_build.check_rank()
        if args.targets != ["ALL"] and marker_build.denominator not in args.targets:
            continue
        marker_build_dict[marker_build.denominator] = marker_build

    param_handler.close()

    return marker_build_dict


def parse_cog_list(args, marker_build_dict):
    """
    Loads the TreeSAPP COG list file into marker_build_dict and check that the args.reftree exists
    :param args: The command-line and default arguments object
    :param marker_build_dict: A dictionary (indexed by marker 5-character 'denominator's) mapping MarkerBuild objects
    :return: marker_build_dict with updated information
    """
    cog_list_file = args.treesapp + os.sep + 'data' + os.sep + 'tree_data' + os.sep + 'cog_list.tsv'
    cog_input_list = open(cog_list_file, 'r', encoding="latin1")
    if args.reftree not in ['i', 'p', 'g']:
        alignment_set = ''
    else:
        alignment_set = args.reftree
    # Load lines from the COG list file and close
    cog_list_lines = [x.strip() for x in cog_input_list.readlines()]
    # Close the COG list file
    cog_input_list.close()

    for marker_input in cog_list_lines:
        if re.match(r'\A#', marker_input):
            continue

        if not re.match(r'\w+\t[A-Z][0-9]{4}\t\w+', marker_input):
            sys.stderr.write("ERROR: entry in cog_list.tsv is incorrectly formatted! Violating line:\n")
            sys.stderr.write(str(marker_input) + "\n")
            sys.stderr.flush()
            sys.exit()

        marker, denominator, description = marker_input.split("\t")
        if args.targets != ["ALL"] and denominator not in args.targets:
            continue

        for ref_code in marker_build_dict:
            if denominator == ref_code:
                marker_build_dict[denominator].description = description
                if description == "phylogenetic_cogs":
                    marker_build_dict[denominator].kind = "phylogenetic_cogs"
                elif description == "rRNA_marker":
                    marker_build_dict[denominator].kind = "phylogenetic_rRNA_cogs"
                else:
                    marker_build_dict[denominator].kind = "functional_cogs"

    if args.reftree not in ['i', 'g', 'p'] and args.reftree not in marker_build_dict.keys():
        sys.stderr.write("ERROR: " + args.reftree + " not found in " + cog_list_file +
                         "! Please use a valid reference tree ID!\n")
        sys.stderr.flush()
        sys.exit()
    return marker_build_dict


def best_match(matches):
    """
    Function for finding the best alignment in a list of HmmMatch() objects
    The best match is based off of the full sequence score
    :param matches: A list of HmmMatch() objects
    :return: The best HmmMatch
    """
    # TODO: Incorporate the alignment intervals to allow for proteins with multiple different functional domains
    # Code currently only permits multi-domains of the same gene
    best_target_hmm = ""
    best_alignment = None
    top_score = 0
    for match in matches:
        # match.print_info()
        if match.full_score > top_score:
            best_alignment = match
            best_target_hmm = match.target_hmm
            top_score = match.full_score
    return best_target_hmm, best_alignment


def parse_domain_tables(args, hmm_domtbl_files, log=None):
    # Check if the HMM filtering thresholds have been set
    if not hasattr(args, "min_e"):
        args.min_e = 0.01
        args.min_acc = 0.6
        args.perc_aligned = 80
    # Print some stuff to inform the user what they're running and what thresholds are being used.
    if args.verbose:
        sys.stdout.write("Filtering HMM alignments using the following thresholds:\n")
        sys.stdout.write("\tMinimum E-value = " + str(args.min_e) + "\n")
        sys.stdout.write("\tMinimum acc = " + str(args.min_acc) + "\n")
        sys.stdout.write("\tMinimum percentage of the HMM covered = " + str(args.perc_aligned) + "%\n")
    sys.stdout.write("Parsing domain tables generated by HMM searches for high-quality matches... ")

    raw_alignments = 0
    seqs_identified = 0
    dropped = 0
    fragmented = 0
    glued = 0
    multi_alignments = 0  # matches of the same query to a different HMM (>1 lines)
    hmm_matches = dict()
    orf_gene_map = dict()

    # TODO: Capture multimatches across multiple domain table files
    for domtbl_file in hmm_domtbl_files:
        rp_marker, reference = re.sub("_domtbl.txt", '', os.path.basename(domtbl_file)).split("_to_")
        domain_table = DomainTableParser(domtbl_file)
        domain_table.read_domtbl_lines()
        distinct_matches, fragmented, glued, multi_alignments, raw_alignments = format_split_alignments(domain_table,
                                                                                                        fragmented,
                                                                                                        glued,
                                                                                                        multi_alignments,
                                                                                                        raw_alignments)
        purified_matches, dropped = filter_poor_hits(args, distinct_matches, dropped)
        complete_gene_hits, dropped = filter_incomplete_hits(args, purified_matches, dropped)

        for match in complete_gene_hits:
            match.genome = reference
            if match.orf not in orf_gene_map:
                orf_gene_map[match.orf] = dict()
            orf_gene_map[match.orf][match.target_hmm] = match
            if match.target_hmm not in hmm_matches.keys():
                hmm_matches[match.target_hmm] = list()

    for orf in orf_gene_map:
        if len(orf_gene_map[orf]) == 1:
            target_hmm = list(orf_gene_map[orf].keys())[0]
            hmm_matches[target_hmm].append(orf_gene_map[orf][target_hmm])
        else:
            optional_matches = [orf_gene_map[orf][target_hmm] for target_hmm in orf_gene_map[orf]]
            target_hmm, match = best_match(optional_matches)
            hmm_matches[target_hmm].append(match)
            multi_alignments += 1
            dropped += (len(optional_matches) - 1)

            if log:
                dropped_annotations = list()
                for optional in optional_matches:
                    if optional.target_hmm != target_hmm:
                        dropped_annotations.append(optional.target_hmm)
                log.write("HMM search annotations for " + orf + ":\n")
                log.write("\tRetained\t" + target_hmm + "\n")
                log.write("\tDropped\t" + ','.join(dropped_annotations) + "\n")

        seqs_identified += 1

    sys.stdout.write("done.\n")

    if seqs_identified == 0 and dropped == 0:
        sys.stderr.write("\tWARNING: No alignments found!\n")
        sys.stderr.write("TreeSAPP is exiting now.\n")
        sys.exit(11)
    if seqs_identified == 0 and dropped > 0:
        sys.stderr.write("\tWARNING: No alignments met the quality cut-offs!\n")
        sys.stderr.write("TreeSAPP is exiting now.\n")
        sys.exit(13)

    sys.stdout.write("\tNumber of markers identified:\n")
    for marker in sorted(hmm_matches):
        sys.stdout.write("\t\t" + marker + "\t" + str(len(hmm_matches[marker])) + "\n")
        # For debugging:
        # for match in hmm_matches[marker]:
        #     match.print_info()
    if args.verbose:
        sys.stdout.write("\tInitial alignments:\t" + str(raw_alignments) + "\n")
        sys.stdout.write("\tAlignments discarded:\t" + str(dropped) + "\n")
        sys.stdout.write("\tFragmented alignments:\t" + str(fragmented) + "\n")
        sys.stdout.write("\tAlignments scaffolded:\t" + str(glued) + "\n")
        sys.stdout.write("\tMulti-alignments:\t" + str(multi_alignments) + "\n")
        sys.stdout.write("\tSequences identified:\t" + str(seqs_identified) + "\n")

    sys.stdout.flush()
    return hmm_matches


def read_colours_file(args, annotation_file):
    """
    Read annotation data from 'annotation_file' and store it in marker_subgroups under the appropriate
    marker and data_type.
    :param args:
    :param annotation_file:
    :return: A dictionary of lists where each list is populated by tuples with start and end leaves
    """
    try:
        style_handler = open(annotation_file, 'r')
    except IOError:
        sys.stderr.write("ERROR: Unable to open " + annotation_file + " for reading!\n")
        sys.exit()

    clusters = dict()
    field_sep = ''
    internal_nodes = True

    line = style_handler.readline()
    # Skip the header
    while line.strip() != "DATA":
        header_fields = line.strip().split(' ')
        if header_fields[0] == "SEPARATOR":
            if header_fields[1] == "SPACE":
                field_sep = ' '
            elif header_fields[1] == "TAB":
                field_sep = '\t'
            else:
                sys.stderr.write("ERROR: Unknown separator used in " + annotation_file + ": " + header_fields[1] + "\n")
                sys.stderr.flush()
                sys.exit()
        line = style_handler.readline()
    # For RGB
    range_line_rgb = re.compile("^(\d+)\|(\d+)" + re.escape(field_sep) +
                                "range" + re.escape(field_sep) +
                                ".*\)" + re.escape(field_sep) +
                                "(.*)$")
    single_node_rgb = re.compile("^(\d+)" + re.escape(field_sep) +
                                 "range" + re.escape(field_sep) +
                                 ".*\)" + re.escape(field_sep) +
                                 "(.*)$")
    lone_node_rgb = re.compile("^(.*)" + re.escape(field_sep) +
                               "range" + re.escape(field_sep) +
                               ".*\)" + re.escape(field_sep) +
                               "(.*)$")

    # For hexadecimal
    range_line = re.compile("^(\d+)\|(\d+)" + re.escape(field_sep) +
                            "range" + re.escape(field_sep) +
                            "#[0-9A-Za-z]{6}" + re.escape(field_sep) +
                            "(.*)$")
    single_node = re.compile("^(\d+)" + re.escape(field_sep) +
                             "range" + re.escape(field_sep) +
                             "#[0-9A-Za-z]{6}" + re.escape(field_sep) +
                             "(.*)$")
    lone_node = re.compile("^(.*)" + re.escape(field_sep) +
                           "range" + re.escape(field_sep) +
                           "#[0-9A-Za-z]{6}" + re.escape(field_sep) +
                           "(.*)$")

    # Begin parsing the data from 4 columns
    line = style_handler.readline().strip()
    while line:
        if range_line.match(line):
            style_data = range_line.match(line)
            start, end, description = style_data.groups()
            internal_nodes = False
        elif range_line_rgb.match(line):
            style_data = range_line_rgb.match(line)
            start, end, description = style_data.groups()
            internal_nodes = False
        elif single_node.match(line):
            style_data = single_node.match(line)
            start, end, description = style_data.group(1), style_data.group(1), style_data.group(2)
        elif single_node_rgb.match(line):
            style_data = single_node_rgb.match(line)
            start, end, description = style_data.group(1), style_data.group(1), style_data.group(2)
        elif lone_node.match(line):
            style_data = lone_node.match(line)
            start, end, description = style_data.group(1), style_data.group(1), style_data.group(2)
        elif lone_node_rgb.match(line):
            style_data = lone_node_rgb.match(line)
            start, end, description = style_data.group(1), style_data.group(1), style_data.group(2)
        else:
            sys.stderr.write("ERROR: Unrecognized line formatting in " + annotation_file + ":\n")
            sys.stderr.write(line + "\n")
            sys.exit()

        description = style_data.groups()[-1]
        if description not in clusters.keys():
            clusters[description] = list()
        clusters[description].append((start, end))

        line = style_handler.readline().strip()

    style_handler.close()

    if args.verbose:
        sys.stdout.write("\tParsed " + str(len(clusters)) +
                         " clades from " + annotation_file + "\n")

    return clusters, internal_nodes


def tax_ids_file_to_leaves(tax_ids_file):
    tree_leaves = list()
    try:
        if sys.version_info > (2, 9):
            cog_tax_ids = open(tax_ids_file, 'r', encoding='utf-8')
        else:
            cog_tax_ids = open(tax_ids_file, 'r')
    except IOError:
        sys.exit('ERROR: Can\'t open ' + str(tax_ids_file) + '!\n')

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
            sys.stderr.write("ValueError: Unexpected number of fields in " + tax_ids_file +
                             ".\nInvoked .split(\'\\t\') on line " + str(line))
            raise ValueError
        leaf = TreeLeafReference(number, translation)
        if lineage:
            leaf.lineage = lineage
            leaf.complete = True
        tree_leaves.append(leaf)

    cog_tax_ids.close()
    return tree_leaves


def read_species_translation_files(args, marker_build_dict):
    """
    :param args:
    :param marker_build_dict: A dictionary (indexed by marker 5-character 'denominator's) mapping MarkerBuild objects
    :return: The taxonomic identifiers for each of the organisms in a tree for all trees
    """

    tree_numbers_translation = dict()
    translation_files = dict()
    tree_resources_dir = os.sep.join([args.treesapp, "data", "tree_data"]) + os.sep

    for denominator in sorted(marker_build_dict.keys()):
        marker_build_obj = marker_build_dict[denominator]
        filename = 'tax_ids_' + str(marker_build_obj.cog) + '.txt'
        translation_files[denominator] = tree_resources_dir + filename

    for denominator in sorted(translation_files.keys()):
        filename = translation_files[denominator]
        tree_numbers_translation[denominator] = tax_ids_file_to_leaves(filename)

    return tree_numbers_translation


def xml_parser(xml_record, term):
    """
    Recursive function for parsing individual xml records
    :param xml_record:
    :param term:
    :return:
    """
    # TODO: Finish this off - would be great for consistently extracting data from xml
    value = None
    if type(xml_record) == str:
        return value
    if term not in xml_record.keys():
        for record in xml_record:
            value = xml_parser(record, term)
            if value:
                return value
            else:
                continue
    else:
        return xml_record[term]
    return value
