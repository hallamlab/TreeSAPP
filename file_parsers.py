#!/usr/bin/env python3

import sys
import os
import re
import logging
from classy import TreeLeafReference, MarkerBuild, Cluster
from utilities import Autovivify, calculate_overlap
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
                           "polynomial_params\tlowest_confident_rank\tdate$")
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
            if marker_build.load_pfit_params(line):
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
    logging.debug(info_string)

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
        sys.exit(0)
    if seqs_identified == 0 and dropped > 0:
        logging.warning("No alignments met the quality cut-offs! TreeSAPP is exiting now.\n")
        sys.exit(0)

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
    unknown = 0
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
        else:
            unknown += 1
        tree_leaves.append(leaf)

    if len(tree_leaves) == unknown:
        logging.error("Lineage information was not properly loaded for " + tax_ids_file + "\n")
        sys.exit(5)

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


def read_phylip_to_dict(phylip_input):
    header_dict = dict()
    tmp_seq_dict = dict()
    seq_dict = dict()
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
            tmp_seq_dict[x] = sequence
            x += 1
        elif 60 >= len(line) >= 1:
            tmp_seq_dict[x] += line
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
    if num_sequences != len(tmp_seq_dict):
        logging.error("Number of lines declared in Phylip header (" + str(num_sequences) +
                      ") does not match number of sequences parsed (" + str(len(tmp_seq_dict)) + ")!\n")
        sys.exit(5)

    x = 0
    while x < num_sequences-1:
        if len(tmp_seq_dict[x]) != aln_length:
            logging.error(header_dict[x] +
                          " sequence length exceeds the stated multiple alignment length (according to header)!\n" +
                          "sequence length = " + str(len(tmp_seq_dict[x])) +
                          ", alignment length = " + str(aln_length) + "\n")
            sys.exit(5)
        else:
            pass
        x += 1

    phylip.close()
    for x in header_dict:
        seq_dict[header_dict[x]] = tmp_seq_dict[x]
    return seq_dict


def read_stockholm_to_dict(sto_file):
    """

    :param sto_file: A Stockholm-formatted multiple alignment file
    :return: A dictionary with sequence headers as keys and sequences as values
    """
    seq_dict = dict()

    try:
        sto_handler = open(sto_file, 'r')
    except IOError:
        logging.error("Unable to open " + sto_file + " for reading!\n")
        sys.exit(3)

    line = sto_handler.readline()
    while line:
        line = line.strip()
        if re.match("^[#|/].*", line):
            # Skip the header (first line) as well as secondary structure lines
            pass
        elif not line:
            pass
        else:
            try:
                seq_name, sequence = line.split()
            except ValueError:
                logging.error("Unexpected line format in " + sto_file + ":\n" + line + "\n")
                sys.exit(3)

            if seq_name not in seq_dict:
                seq_dict[seq_name] = ""
            seq_dict[seq_name] += re.sub('\.', '-', sequence.upper())
        line = sto_handler.readline()

    return seq_dict


def parse_blast_results(args, blast_tables, cog_list):
    """
    Returns an Autovivification of purified (eg. non-redundant) BLAST hits.
    :param args: Command-line argument object from get_options and check_parser_arguments
    :param blast_tables: file produced by BLAST alignment
    :param cog_list: list of COGs included in analysis pipeline
    """

    logging.info("Parsing BLAST results... ")

    # reg_cog_id = re.compile(r'.*(.{5})\Z')
    counter = 0
    purified_blast_hits = Autovivify()
    contigs = {}
    hit_logger = dict()
    alignment_count = 0

    for blast_table in blast_tables:
        try:
            blast_results = open(blast_table, 'r')
        except IOError:
            logging.error("Cannot open BLAST output file " + blast_table + "\n")
            sys.exit(3)

        identifier = 0
        for line in blast_results:
            # Clear variables referencing the contig, COG, qstart, qend, reference start, reference end, and bitscore
            # Interpret the BLAST hit, and assign the details accordingly
            alignment_count += 1
            temp_contig, temp_detailed_cog, _, _, _, _, temp_query_start, temp_query_end, temp_ref_start, temp_ref_end, _, temp_bitscore = line.split('\t')
            temp_ref_end = int(temp_ref_end)
            temp_ref_start = int(temp_ref_start)
            temp_query_end = int(temp_query_end)
            temp_query_start = int(temp_query_start)
            temp_bitscore = float(temp_bitscore)

            # Skip to next BLAST hit if bit score is less than user-defined minimum
            if temp_bitscore <= args.bitscore:
                continue

            # Determine the direction of the hit relative to the reference
            direction = 'forward'
            if temp_ref_start > temp_ref_end:
                temp = temp_ref_start
                temp_ref_start = temp_ref_end
                temp_ref_end = temp
                direction = 'reverse'
            if temp_query_start > temp_query_end:
                temp = temp_query_start
                temp_query_start = temp_query_end
                temp_query_end = temp
                if direction == 'reverse':
                    logging.error("Confusing BLAST result!\n" +
                                  "Please notify the authors about " +
                                  temp_contig + ' at ' +
                                  temp_detailed_cog +
                                  " q(" + str(temp_query_end) + '..' + str(temp_query_start) + ")," +
                                  " r(" + str(temp_ref_end) + '..' + str(temp_ref_start) + ")")
                    sys.exit(3)
                direction = 'reverse'

            # This limitation is so-far not necessary
            # result = reg_cog_id.match(temp_detailed_cog)
            # if result:
            #     tempCOG = result.group(1)
            result = '_'.join(temp_detailed_cog.split('_')[1:])
            if result:
                tempCOG = result
            else:
                sys.exit('ERROR: Could not detect the COG of sequence ' + temp_detailed_cog)

            # Save contig details to the list
            if temp_contig not in contigs:
                contigs[temp_contig] = {}

            if identifier not in contigs[temp_contig]:
                contigs[temp_contig][identifier] = {}

            contigs[temp_contig][identifier]['bitscore'] = temp_bitscore
            contigs[temp_contig][identifier]['cog'] = tempCOG
            contigs[temp_contig][identifier]['seq_start'] = temp_query_start
            contigs[temp_contig][identifier]['seq_end'] = temp_query_end
            contigs[temp_contig][identifier]['direction'] = direction
            contigs[temp_contig][identifier]['validity'] = True
            identifier += 1

        # Close the file
        blast_results.close()

    # Purify the BLAST hits
    # For each contig sorted by their string-wise comparison...
    for contig in sorted(contigs.keys()):
        identifier = 0

        # create tuple array to sort
        IDs = []
        for raw_identifier in sorted(contigs[contig].keys()):
            base_start = contigs[contig][raw_identifier]['seq_start']
            IDs.append((raw_identifier, base_start))
        _IDs = sorted(IDs, key=lambda x: x[1])
        IDs = [x[0] for x in _IDs]

        base_blast_result_raw_identifier = IDs.pop()
        contigs[contig][base_blast_result_raw_identifier]['validity'] = True
        base_bitscore = contigs[contig][base_blast_result_raw_identifier]['bitscore']
        base_cog = contigs[contig][base_blast_result_raw_identifier]['cog']
        base_start = contigs[contig][base_blast_result_raw_identifier]['seq_start']
        base_end = contigs[contig][base_blast_result_raw_identifier]['seq_end']
        direction = contigs[contig][base_blast_result_raw_identifier]['direction']
        base_length = base_end - base_start

        # Compare the BLAST hit (base) against all others
        # There may be several opinions about how to do this. This way is based on the original MLTreeMap
        # ----A----  --C--
        #        ---B---
        # A kills B, B kills C. (Another approach would be to let C live,
        # but the original MLTreeMap authors don't expect C to be useful)
        for check_blast_result_raw_identifier in IDs:
            check_bitscore = contigs[contig][check_blast_result_raw_identifier]['bitscore']
            check_cog = contigs[contig][check_blast_result_raw_identifier]['cog']
            check_start = contigs[contig][check_blast_result_raw_identifier]['seq_start']
            check_end = contigs[contig][check_blast_result_raw_identifier]['seq_end']
            check_length = check_end - check_start

            # Compare the base and check BLAST hits
            info = Autovivify()
            info['base']['start'] = base_start
            info['base']['end'] = base_end
            info['check']['start'] = check_start
            info['check']['end'] = check_end
            overlap = calculate_overlap(info)
            counter += 1

            # Check for validity for hits with overlap
            if overlap == 0:
                base_blast_result_raw_identifier = check_blast_result_raw_identifier
                base_bitscore = check_bitscore
                base_cog = check_cog
                base_start = check_start
                base_end = check_end
                base_length = check_length
                contigs[contig][base_blast_result_raw_identifier]['validity'] = True
            else:
                if overlap > 0.5*base_length and base_bitscore < check_bitscore:
                    contigs[contig][base_blast_result_raw_identifier]['validity'] = False
                    base_blast_result_raw_identifier = check_blast_result_raw_identifier
                    base_bitscore = check_bitscore
                    base_cog = check_cog
                    base_start = check_start
                    base_end = check_end
                    base_length = check_length
                    contigs[contig][base_blast_result_raw_identifier]['validity'] = True
                elif overlap > 0.5*check_length and check_bitscore < base_bitscore:
                    contigs[contig][check_blast_result_raw_identifier]['validity'] = False
                elif base_start == check_start and base_end == check_end:
                    # If both are the same, keep only the one with the smaller identifier
                    if check_blast_result_raw_identifier > base_blast_result_raw_identifier:
                        contigs[contig][check_blast_result_raw_identifier]['validity'] = False

        # Set validity to 0 if COG is not in list of TreeSAPP COGs
        if base_cog not in cog_list['all_cogs']:
            contigs[contig][base_blast_result_raw_identifier]['validity'] = False
            logging.warning("WARNING: " + base_cog + " not in list of TreeSAPP markers\n")

        # Save purified hits for valid base hits
        for base_blast_result_raw_identifier in IDs:
            base_bitscore = contigs[contig][base_blast_result_raw_identifier]['bitscore']
            base_cog = contigs[contig][base_blast_result_raw_identifier]['cog']
            base_start = contigs[contig][base_blast_result_raw_identifier]['seq_start']
            base_end = contigs[contig][base_blast_result_raw_identifier]['seq_end']
            direction = contigs[contig][base_blast_result_raw_identifier]['direction']
            if contigs[contig][base_blast_result_raw_identifier]['validity']:
                purified_blast_hits[contig][identifier]['bitscore'] = base_bitscore
                purified_blast_hits[contig][identifier]['cog'] = base_cog
                purified_blast_hits[contig][identifier]['start'] = base_start
                purified_blast_hits[contig][identifier]['end'] = base_end
                purified_blast_hits[contig][identifier]['direction'] = direction
                purified_blast_hits[contig][identifier]['is_already_placed'] = False
                identifier += 1

    # Print the BLAST results for each contig
    for contig in sorted(purified_blast_hits.keys()):
        outfile = args.output_dir_var + contig + '_blast_result_purified.txt'
        out = open(outfile, 'w')
        sorting_hash = {}

        # Identify the first instance of each bitscore
        for identifier in sorted(purified_blast_hits[contig].keys()):
            if not purified_blast_hits[contig][identifier]['bitscore'] in sorting_hash:
                sorting_hash[purified_blast_hits[contig][identifier]['bitscore']] = {}
            sorting_hash[purified_blast_hits[contig][identifier]['bitscore']][identifier] = 1

        # Print the (potentially reduced set of) BLAST results ordered by decreasing bitscore
        for bitscore in sorted(sorting_hash.keys(), reverse=True):
            for identifier in sorted(sorting_hash[bitscore]):
                marker = purified_blast_hits[contig][identifier]['cog']
                if marker not in hit_logger:
                    hit_logger[marker] = 0
                hit_logger[marker] += 1
                out.write(contig + '\t' + str(purified_blast_hits[contig][identifier]['start']) + '\t' +
                          str(purified_blast_hits[contig][identifier]['end']) + '\t' +
                          str(purified_blast_hits[contig][identifier]['direction']) + '\t' +
                          purified_blast_hits[contig][identifier]['cog'] + '\t' + str(bitscore) + '\n')

        out.close()
    logging.info("done.\n")

    logging.debug("\t" + str(alignment_count) + " intial BLAST alignments found.\n")
    total = 0
    for n in hit_logger.values():
        total += n
    logging.debug("\t" + str(total) + " purified BLAST alignments:\n" +
                  "\n".join(["\t\t" + str(hit_logger[marker]) + " " + marker for marker in hit_logger]))

    return purified_blast_hits


def blastp_parser(args, blast_hits_purified):
    """
    For each contig, produces a file similar to the Genewise output file
    (this is in cases where Genewise is unnecessary because it is already an AA sequence.
    :param args: Command-line argument object from get_options and check_parser_arguments
    :param blast_hits_purified: Parsed blastp outputs
    :return blastp_summary_files: Autovivification of the output file for each contig.
    """

    blastp_summary_files = Autovivify()

    reg_header = re.compile(r'\A>')

    for contig in sorted(blast_hits_purified.keys()):
        output_file = args.output_dir_var + contig + '_blast_result_summary.txt'
        try:
            output = open(output_file, 'w')
        except IOError:
            sys.exit('ERROR: Unable to open ' + output_file + '!\n')
        blastp_summary_files[contig][output_file] = 1
        shortened_sequence_file = args.output_dir_var + contig + '_sequence_shortened.txt'
        try:
            sequence_file = open(shortened_sequence_file, 'r')
        except IOError:
            sys.exit('ERROR: Could not open ' + shortened_sequence_file + '!\n')
        flag_seq = 0
        sequence = ''

        # Get the sequence from the shortened sequence file
        for line in sequence_file:
            if reg_header.search(line):
                if flag_seq == 1:
                    sys.exit('ERROR: Unexpected multiple shortened sequences found!\n')
                flag_seq = 1
                continue
            else:
                line.strip()
                sequence += line

        # Write the output file to imitate the Genewise results
        for count in sorted(blast_hits_purified[contig].keys()):
            output.write(str(blast_hits_purified[contig][count]['cog']) + '\t')
            output.write(str(blast_hits_purified[contig][count]['start']) + '\t')
            output.write(str(blast_hits_purified[contig][count]['end']) + '\t')
            output.write(str(blast_hits_purified[contig][count]['direction']) + '\t')
            output.write(str(sequence) + '\n')
        sequence_file.close()
        output.close()

    return blastp_summary_files


def read_uc(uc_file):
    """
    Function to read a USEARCH cluster (.uc) file

    :param uc_file: Path to a .uc file produced by USEARCH
    :return: Dictionary where keys are representative cluster headers and the values are headers of identical sequences
    """
    cluster_dict = dict()
    rep_len_map = dict()
    try:
        uc = open(uc_file, 'r')
    except IOError:
        logging.error("Unable to open USEARCH cluster file " + uc_file + " for reading!\n")
        sys.exit(13)

    line = uc.readline()
    # Find all clusters with multiple identical sequences
    while line:
        cluster_type, num_id, length, identity, _, _, _, cigar, header, representative = line.strip().split("\t")
        if cluster_type == "S":
            cluster_dict[num_id] = Cluster('>' + header)
            rep_len_map['>' + header] = length
        elif cluster_type == "H":
            cluster_dict[num_id].members.append(['>' + header, identity])
        elif cluster_type == "C":
            pass
        else:
            logging.error("Unexpected cluster type '" + str(cluster_type) + "' in " + uc_file + "\n")
            sys.exit(13)
        line = uc.readline()
    return cluster_dict
