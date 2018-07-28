#!/usr/bin/env python3

import sys
import os
import re
import logging
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
        logging.error("\tUnable to open " + ref_build_parameters + " for reading.\n")
        sys.exit(5)

    header_re = re.compile("^name\tcode\tmolecule\tsub_model\tcluster_identity\t"
                           "class_dist\torder_dist\tfamily_dist\tgenus_dist\tspecies_dist\t"
                           "lowest_confident_rank\tdate$")
    if not header_re.match(param_handler.readline().strip()):
        logging.error("Header of '" + ref_build_parameters + "' is unexpected!")
        sys.exit(5)

    logging.info("Reading build parameters of reference markers... ")
    skipped_lines = []
    missing_info = []
    marker_build_dict = dict()
    for line in param_handler:
        if header_re.match(line):
            continue
        line = line.strip()
        if line[0] == '#':
            skipped_lines.append(line)
            continue
        marker_build = MarkerBuild(line)
        if args.targets != ["ALL"] and marker_build.denominator not in args.targets:
            skipped_lines.append(line)
        else:
            marker_build_dict[marker_build.denominator] = marker_build
            if marker_build.load_rank_distances(line) > 0:
                missing_info.append(marker_build)
            marker_build.check_rank()
    param_handler.close()

    logging.info("done.\n")

    logging.debug("Rank distance information missing for:\n\t" +
                  "\n\t".join([mb.cog + '-' + mb.denominator for mb in missing_info]) + "\n")
    logging.debug("Skipped the following lines:\n\t" +
                  "\n\t".join(skipped_lines) + "\n")
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
    # Load lines from the COG list file and close
    cog_list_lines = [x.strip() for x in cog_input_list.readlines()]
    # Close the COG list file
    cog_input_list.close()

    for marker_input in cog_list_lines:
        if re.match(r'\A#', marker_input):
            continue

        if not re.match(r'\w+\t[A-Z][0-9]{4}\t\w+', marker_input):
            message = "Entry in cog_lit.tsv is incorrectly formatted! Violating line:\n" + str(marker_input)
            logging.error(message)
            sys.exit(5)

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
        logging.error(args.reftree + " not found in " + cog_list_file + "! Please use a valid reference tree ID!\n")
        sys.exit(5)
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


def parse_domain_tables(args, hmm_domtbl_files):
    # Check if the HMM filtering thresholds have been set
    if not hasattr(args, "min_e"):
        args.min_e = 0.01
        args.min_acc = 0.6
        args.perc_aligned = 80
    # Print some stuff to inform the user what they're running and what thresholds are being used.
    info_string = "Filtering HMM alignments using the following thresholds:\n"
    info_string += "\tMinimum E-value = " + str(args.min_e) + "\n"
    info_string += "\tMinimum acc = " + str(args.min_acc) + "\n"
    info_string += "\tMinimum percentage of the HMM covered = " + str(args.perc_aligned) + "%\n"
    logging.info(info_string)

    logging.info("Parsing HMMER domain tables for high-quality matches... ")

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

            dropped_annotations = list()
            for optional in optional_matches:
                if optional.target_hmm != target_hmm:
                    dropped_annotations.append(optional.target_hmm)
            logging.debug("HMM search annotations for " + orf +
                          ":\n\tRetained\t" + target_hmm +
                          "\n\tDropped\t" + ','.join(dropped_annotations) + "\n")

        seqs_identified += 1

    logging.info("done.\n")

    if seqs_identified == 0 and dropped == 0:
        logging.warning("No alignments found! TreeSAPP is exiting now.\n")
        sys.exit(5)
    if seqs_identified == 0 and dropped > 0:
        logging.warning("No alignments met the quality cut-offs! TreeSAPP is exiting now.\n")
        sys.exit(5)

    alignment_stat_string = "\tNumber of markers identified:\n"
    for marker in sorted(hmm_matches):
        alignment_stat_string += "\t\t" + marker + "\t" + str(len(hmm_matches[marker])) + "\n"
        # For debugging:
        # for match in hmm_matches[marker]:
        #     match.print_info()

    alignment_stat_string += "\tInitial alignments:\t" + str(raw_alignments) + "\n"
    alignment_stat_string += "\tAlignments discarded:\t" + str(dropped) + "\n"
    alignment_stat_string += "\tFragmented alignments:\t" + str(fragmented) + "\n"
    alignment_stat_string += "\tAlignments scaffolded:\t" + str(glued) + "\n"
    alignment_stat_string += "\tMulti-alignments:\t" + str(multi_alignments) + "\n"
    alignment_stat_string += "\tSequences identified:\t" + str(seqs_identified) + "\n"

    logging.debug(alignment_stat_string)
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
        logging.error("Unable to open " + annotation_file + " for reading!\n")
        sys.exit(5)

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
                logging.error("ERROR: Unknown separator used in " + annotation_file + ": " + header_fields[1] + "\n")
                sys.exit(5)
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
            logging.error("ERROR: Unrecognized line formatting in " + annotation_file + ":\n" + line + "\n")
            sys.exit(5)

        description = style_data.groups()[-1]
        if description not in clusters.keys():
            clusters[description] = list()
        clusters[description].append((start, end))

        line = style_handler.readline().strip()

    style_handler.close()

    logging.debug("\tParsed " + str(len(clusters)) + " clades from " + annotation_file + "\n")

    return clusters, internal_nodes


def tax_ids_file_to_leaves(tax_ids_file):
    tree_leaves = list()
    try:
        if sys.version_info > (2, 9):
            cog_tax_ids = open(tax_ids_file, 'r', encoding='utf-8')
        else:
            cog_tax_ids = open(tax_ids_file, 'r')
    except IOError:
        logging.error("Unable to open " + tax_ids_file + "\n")
        sys.exit(5)

    for line in cog_tax_ids:
        line = line.strip()
        try:
            fields = line.split("\t")
        except ValueError:
            logging.error('ValueError: .split(\'\\t\') on ' + str(line) +
                          " generated " + str(len(line.split("\t"))) + " fields.\n")
            sys.exit(5)
        if len(fields) == 2:
            number, translation = fields
            lineage = ""
        elif len(fields) == 3:
            number, translation, lineage = fields
        else:
            logging.error("ValueError: Unexpected number of fields in " + tax_ids_file +
                          ".\nInvoked .split(\'\\t\') on line " + str(line) + "\n")
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


def read_phylip(phylip_input):
    header_dict = dict()
    alignment_dict = dict()
    x = 0

    try:
        phylip = open(phylip_input, 'r')
    except IOError:
        logging.error("Unable to open the Phylip file (" + phylip_input + ") provided for reading!\n")
        sys.exit(5)

    line = phylip.readline()
    try:
        num_sequences, aln_length = line.strip().split(' ')
        num_sequences = int(num_sequences)
        aln_length = int(aln_length)
    except ValueError:
        logging.error("Phylip file is not formatted correctly!\n" +
                      "Header must contain 2 space-separated fields (number of sequences and alignment length).\n")
        sys.exit(5)

    line = phylip.readline()
    while line:
        line = line.strip()
        if len(line.split()) == 2:
            # This is the introduction set: header, sequence
            header, sequence = line.split()
            header_dict[x] = header
            alignment_dict[x] = sequence
            x += 1
        elif 60 >= len(line) >= 1:
            alignment_dict[x] += line
            x += 1
        elif line == "":
            # Reset accumulator on blank lines
            x = 0
        else:
            logging.error("Unexpected line in Phylip file:\n" + line + "\n")
            sys.exit(5)
        line = phylip.readline()

        if x > num_sequences:
            logging.error("Accumulator has exceeded the number of sequences in the file (according to header)!\n")
            sys.exit(5)

    # Check that the alignment length matches that in the header line
    x = 0
    while x < num_sequences-1:
        if len(alignment_dict[x]) != aln_length:
            logging.error(header_dict[x] +
                          " sequence length exceeds the stated multiple alignment length (according to header)!\n" +
                          "sequence length = " + str(len(alignment_dict[x])) +
                          ", alignment length = " + str(aln_length) + "\n")
            sys.exit(5)
        else:
            pass
        x += 1

    phylip.close()
    return header_dict, alignment_dict
