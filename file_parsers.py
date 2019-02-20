#!/usr/bin/env python3

import sys
import os
import re
import logging
from classy import TreeLeafReference, MarkerBuild, Cluster
from HMMER_domainTblParser import DomainTableParser, format_split_alignments, filter_incomplete_hits, filter_poor_hits
from fasta import read_fasta_to_dict

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

    header_re = re.compile("\t".join(["name", "code",
                                      "molecule", "sub_model", "marker_info", "cluster_identity", "ref_sequences",
                                      "tree_tool", "poly-params", "lowest_reliable_rank",
                                      "last_updated", "description"]))
    if not header_re.match(param_handler.readline().strip()):
        logging.error("Header of '" + ref_build_parameters + "' is unexpected!")
        sys.exit(5)

    logging.debug("Reading build parameters of reference markers... ")
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
        marker_build = MarkerBuild()
        marker_build.load_build_params(line)
        if args.targets != ["ALL"] and marker_build.denominator not in args.targets:
            skipped_lines.append(line)
        else:
            if marker_build.denominator in marker_build_dict:
                logging.debug("Multiple '" + marker_build.denominator + "' codes in " + ref_build_parameters +
                              ". Previous entry in marker_build_dict being overwritten...\n")
            marker_build_dict[marker_build.denominator] = marker_build
            if marker_build.load_pfit_params(line):
                missing_info.append(marker_build)
            marker_build.check_rank()
    param_handler.close()

    logging.debug("done.\n")

    if missing_info:
        logging.debug("Rank distance information missing for:\n\t" +
                      "\n\t".join([mb.cog + '-' + mb.denominator for mb in missing_info]) + "\n")
    if skipped_lines:
        logging.debug("Skipped the following lines:\n\t" +
                      "\n\t".join(skipped_lines) + "\n")

    if len(marker_build_dict) == 0:
        logging.error("No reference package information was parsed.\n" +
                      "Is your target '" + ','.join(args.targets) + "' in " + ref_build_parameters + "?\n")
        sys.exit(3)
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
                    marker_build_dict[denominator].kind = "phylogenetic"
                elif description == "rRNA_marker":
                    marker_build_dict[denominator].kind = "phylogenetic_rRNA"
                else:
                    marker_build_dict[denominator].kind = "functional"

    if args.reftree not in ['i', 'g', 'p'] and args.reftree not in marker_build_dict.keys():
        logging.error(args.reftree + " not found in " + cog_list_file + "! Please use a valid reference tree ID!\n")
        sys.exit(5)
    return marker_build_dict


def read_graftm_classifications(assignment_file):
    """
    Function for reading the _read_tax.tsv file generated by graftM.
    Sequences that have either multiple genes and/or subunits encoded or have homologous regions separated by
     a divergent sequence are recorded as "_split_" and only the first split is recorded and analyzed.

    :param assignment_file: Path to the _read_tax.tsv file
    :return: Dictionary indexed by taxonomic lineage whose values are headers of classified sequences
    """
    assignments = dict()
    assignments_handle = open(assignment_file, 'r')
    tax_lines = assignments_handle.readlines()
    assignments_handle.close()

    for line in tax_lines:
        fields = line.strip().split('\t')
        try:
            header, classified = fields
            if re.search("_split_", header):
                split = int(header.split('_')[-1])
                if split > 1:
                    continue
                else:
                    header = re.sub("_split_.*", '', header)
            classified = '; '.join([re.sub('e\d+$', '', taxon) for taxon in classified.split('; ')])
            if header and classified:
                if classified not in assignments:
                    assignments[classified] = list()
                assignments[classified].append(header)
        except ValueError:
            logging.error("Unable to parse line:" + str(line))
            sys.exit(21)

    return assignments


def parse_assignments(classified_lines: list):
    """
    Parses the marker_contig_map.tsv lines loaded to retrieve lineage assignment and marker information

    :param classified_lines: A list of classification lines returned by read_marker_classification_table
    :return: A dictionary of lineage information for each assignment, indexed by the marker gene it was classified as
    """
    assignments = dict()
    for fields in classified_lines:
        _, header, marker, _, raw_tax, rob_class, _, _, _, _, _ = fields
        if marker and rob_class:
            if marker not in assignments:
                assignments[marker] = dict()
            if rob_class not in assignments[marker]:
                assignments[marker][rob_class] = list()
            assignments[marker][rob_class].append(header)
    return assignments


def read_marker_classification_table(assignment_file):
    """
    Function for reading the tabular assignments file (currently marker_contig_map.tsv)
    Assumes column 2 is the TreeSAPP assignment and column 3 is the sequence header
    (leaving 1 for marker name and 4 for numerical abundance)

    :param assignment_file: Path to the file containing sequence phylogenetic origin and assignment
    :return: A list of lines that have been split by tabs into lists themselves
    """
    classified_lines = list()
    header = "Sample\tQuery\tMarker\tLength\tTaxonomy\tConfident_Taxonomy\tAbundance\tiNode\tLWR\tEvoDist\tDistances\n"

    assignments_handle = open(assignment_file, 'r')
    # This is the header line
    if assignments_handle.readline() != header:
        logging.error("Header of assignments file is unexpected!\n")
        sys.exit(21)

    # First line in the table containing data
    line = assignments_handle.readline()
    n_fields = len(header.split("\t"))
    while line:
        fields = line.strip().split('\t')
        if len(fields) == n_fields:
            classified_lines.append(fields)
        else:
            logging.error("Unable to parse line:\n" + str(line))
            sys.exit(21)
        line = assignments_handle.readline()
    assignments_handle.close()

    return classified_lines


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


def parse_domain_tables(args, hmm_domtbl_files, single=True):
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
            if single is False:
                for match in optional_matches:
                    hmm_matches[match.target_hmm].append(match)
                seqs_identified += len(optional_matches) - 1
            else:
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
                              "\n\tDropped\t\t" + ','.join(dropped_annotations) + "\n")

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
    if isinstance(xml_record, str):
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


def read_uc(uc_file):
    """
    Function to read a USEARCH cluster (.uc) file

    :param uc_file: Path to a .uc file produced by USEARCH
    :return: Dictionary where keys are numerical identifiers and values are Cluster objects
        The Cluster object
    """
    cluster_dict = dict()
    rep_len_map = dict()
    try:
        uc = open(uc_file, 'r')
    except IOError:
        logging.error("Unable to open USEARCH cluster file " + uc_file + " for reading!\n")
        sys.exit(13)

    logging.debug("Reading usearch cluster file... ")

    # Find all clusters with multiple identical sequences
    for line in uc:
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

    uc.close()
    logging.debug("done.\n")
    return cluster_dict


def read_rpkm(rpkm_output_file):
    """
    Read the CSV file written by rpkm. A header and line with unmapped reads is expected and are skipped.
    Each line is expected to have 4 elements: Sample ID, sequence name, number of reads recruited, RPKM
    :param rpkm_output_file: A file path
    :return: Dictionary mapping contig names to floats
    """
    rpkm_values = dict()

    try:
        rpkm_stats = open(rpkm_output_file)
    except IOError:
        logging.error("Unable to open " + rpkm_output_file + " for reading!\n")
        sys.exit(13)

    # Skip the header
    next(rpkm_stats)
    # Skip the line with unaligned reads
    next(rpkm_stats)
    for line in rpkm_stats:
        # Line format is Sample ID (output file name), sequence name, number of reads recruited, RPKM
        try:
            _, seq_name, _, rpkm = line.strip().split(',')
        except ValueError:
            n_values = str(len(line.split(',')))
            logging.error("Unexpected line format in RPKM file - should contain 4 elements, "
                          "" + n_values + " encountered. Offending line:\n" + line + "\n")
            sys.exit(13)
        rpkm_values[seq_name] = float(rpkm)
    rpkm_stats.close()
    return rpkm_values


def validate_alignment_trimming(msa_files: list, unique_ref_headers: set, queries_mapped=False, min_seq_length=30):
    """
    Parse a list of multiple sequence alignment (MSA) files and determine whether the multiple alignment:
        1. is shorter then the min_seq_length (30 by default)
        2. is missing any reference sequences
    The number of query sequences discarded - these may have been added by hmmalign or PaPaRa - is returned via a string
    NOTE: Initially design for sequences records with numeric names (e.g. >4889) but accomodates other TreeSAPP formats
    :param msa_files: A list of either Phylip- or FASTA-formatted MSA files
    :param unique_ref_headers: A set of all headers that were in the untrimmed MSA
    :param queries_mapped: Boolean indicating whether sequences should be present in addition to reference sequences.
        While query sequences _could_ be identified as any that are not in unique_ref_headers,
        queries have names that are negative integers for more rapid and scalable identification
    :param min_seq_length: Optional minimum unaligned (no '-'s) length a sequence must exceed to be retained
    :return: 1. Dictionary indexed by MSA file name mapping to FASTA-dictionaries and
    2. A string mapping the number of query sequences removed from each MSA file
    """
    discarded_seqs_string = ""
    successful_multiple_alignments = dict()
    for multi_align_file in msa_files:
        filtered_multi_align = dict()
        discarded_seqs = list()
        num_queries_retained = 0
        f_ext = multi_align_file.split('.')[-1]

        # Read the multiple alignment file
        if re.search("phy", f_ext):  # File is in Phylip format
            seq_dict = read_phylip_to_dict(multi_align_file)
        elif re.match("^f", f_ext):  # This is meant to match all fasta extensions
            seq_dict = read_fasta_to_dict(multi_align_file)
        else:
            logging.error("Unable to detect file format of " + multi_align_file + ".\n")
            sys.exit(13)

        # Parse the MSA dict and ensure headers are integer-compatible
        multi_align = dict()
        for seq_name in seq_dict:
            seq = seq_dict[seq_name]
            try:
                int(seq_name)
            except ValueError:
                if re.match("^_\d+", seq_name):
                    seq_name = re.sub("^_", '-', seq_name)
                # The section of regular expresion after '_' needs to match denominator and refpkg names
                elif re.match("^\d+_\w{3,7}$", seq_name):
                    seq_name = seq_name.split('_')[0]
                else:
                    logging.error("Unexpected sequence name " + seq_name +
                                  " detected in " + multi_align_file + ".\n")
                    sys.exit(13)
            multi_align[seq_name] = seq
        if len(multi_align) == 0:
            logging.error("No sequences were read from " + multi_align_file + ".\n")
            sys.exit(3)
        # The numeric identifiers make it easy to maintain order in the Phylip file by a numerical sort
        for seq_name in sorted(multi_align, key=int):
            seq_dummy = re.sub('-', '', multi_align[seq_name])
            if len(seq_dummy) < min_seq_length:
                discarded_seqs.append(seq_name)
            else:
                filtered_multi_align[seq_name] = multi_align[seq_name]
                # The negative integers indicate this is a query sequence
                if seq_name[0] == '-':
                    num_queries_retained += 1

        discarded_seqs_string += "\n\t\t" + multi_align_file + " = " + str(len(discarded_seqs))
        multi_align_seq_names = set(multi_align.keys())
        filtered_multi_align_seq_names = set(filtered_multi_align.keys())
        if len(discarded_seqs) == len(multi_align.keys()):
            # Throw an error if the final trimmed alignment is shorter than min_seq_length, and therefore empty
            logging.warning("Multiple sequence alignment in " + multi_align_file +
                            " is shorter than minimum sequence length threshold (" + str(min_seq_length) +
                            ").\nThese sequences will not be analyzed.\n")
        elif not unique_ref_headers.issubset(multi_align_seq_names):
            # Testing whether there were more sequences in the untrimmed alignment than the trimmed one
            logging.error("Reference sequences in " + multi_align_file + " were removed during alignment trimming.\n" +
                          "This suggests either truncated sequences or the initial reference alignment was terrible.\n")
            sys.exit(3)
        # Calculate the number of reference sequences removed
        elif not unique_ref_headers.issubset(filtered_multi_align_seq_names):
            logging.warning("Reference sequences shorter than the minimum character length (" +
                            str(min_seq_length) + ") in " + multi_align_file +
                            " were removed after alignment trimming.\n" +
                            "These sequences will not be analyzed.\n")
        # Ensure that there is at least 1 query sequence retained after trimming the multiple alignment
        elif queries_mapped and num_queries_retained == 0:
            logging.debug("No query sequences in " + multi_align_file + " were retained after trimming.\n")
        else:
            successful_multiple_alignments[multi_align_file] = filtered_multi_align

        if multi_align_file in successful_multiple_alignments:
            discarded_seqs_string += " (retained)"
        else:
            discarded_seqs_string += " (removed)"

    return successful_multiple_alignments, discarded_seqs_string


def multiple_alignment_dimensions(seq_dict, mfa_file):
    """
    Checks to ensure all sequences are the same length and returns a tuple of (nrow, ncolumn)

    :param seq_dict: A dictionary containing headers as keys and sequences as values
    :param mfa_file: The name of the multiple alignment FASTA file being validated
    :return: tuple = (nrow, ncolumn)
    """
    sequence_length = 0
    for seq_name in seq_dict:
        sequence = seq_dict[seq_name]
        if sequence_length == 0:
            sequence_length = len(sequence)
        elif sequence_length != len(sequence) and sequence_length > 0:
            logging.error("Number of aligned columns is inconsistent in " + mfa_file + "!\n")
            sys.exit(3)
        else:
            pass
            # Sequence is the right length, carrying on
    return len(seq_dict), sequence_length
