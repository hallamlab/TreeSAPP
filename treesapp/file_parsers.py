import sys
import os
import re
import logging
from glob import glob

from collections import namedtuple

from treesapp import fasta
from treesapp import hmmer_tbl_parser
from treesapp import phylo_seq
from treesapp import utilities
from treesapp import logger

LOGGER = logging.getLogger(logger.logger_name())


def write_classification_table(tree_saps: dict, sample_name: str, output_file: str, append=False) -> None:
    """
    Write the final classification table

    :param tree_saps: A dictionary containing PQuery objects
    :param sample_name: String representing the name of the sample (i.e. Assign.sample_prefix)
    :param output_file: Path to write the classification table
    :param append: Boolean controlling whether a file is overwritten or appended to
    :return: None
    """
    try:
        if append:
            tab_out = open(output_file, 'a')
            tab_out_string = ""
        else:
            tab_out = open(output_file, 'w')
            tab_out_string = "Sample\tQuery\tMarker\tStart_pos\tEnd_pos\tTaxonomy\tAbundance\t" \
                             "iNode\tE-value\tLWR\tEvoDist\tDistances\n"
    except IOError:
        LOGGER.error("Unable to open " + output_file + " for writing!\n")
        sys.exit(3)

    for refpkg_name in tree_saps:
        for tree_sap in tree_saps[refpkg_name]:  # type: phylo_seq.PQuery
            if not tree_sap.classified:
                continue
            pplace = tree_sap.consensus_placement
            tab_out_string += '\t'.join([sample_name,
                                         re.sub(r"\|{0}\|\d+_\d+$".format(tree_sap.ref_name), '', tree_sap.place_name),
                                         tree_sap.ref_name,
                                         str(tree_sap.start),
                                         str(tree_sap.end),
                                         tree_sap.recommended_lineage,
                                         str(tree_sap.abundance),
                                         str(pplace.edge_num),
                                         str(tree_sap.evalue),
                                         str(pplace.like_weight_ratio),
                                         str(tree_sap.avg_evo_dist),
                                         tree_sap.distances]) + "\n"

    tab_out.write(tab_out_string)
    tab_out.close()

    return


def parse_assignments(classified_lines: list) -> dict:
    """
    Parses the classifications.tsv lines loaded to retrieve lineage assignment and marker information

     Now also looks for fragments of identical parent sequences that were individually classified.
     If these are found, the longest fragment is selected for classification.
     Removing redundant fragments is performed to ease downstream classifications as
     the number of query sequences is calculated separately and we don't want to exceed 100% classifications!
      Alternatively, the number of query sequences could be calculated from the classification tables
      but we don't think this is the best route as unclassified seqs would wreak havoc.

    :param classified_lines: A list of classification lines returned by read_classification_table
    :return: A dictionary of lineage information for each assignment, indexed by the marker gene it was classified as
    """
    classified = namedtuple("classified", ["refpkg", "taxon", "length"])
    assignments = dict()
    unique_headers = dict()  # Temporary storage for classified sequences prior to filtering
    dups = set()  # For storing the names of all query sequences that were split and classified separately
    for fields in classified_lines:
        _, header, marker, start_pos, end_pos, rec_tax, _, _, _, _, _, _ = fields
        length = int(end_pos) - int(start_pos)
        if marker and rec_tax:
            if marker not in assignments:
                assignments[marker] = dict()
            if rec_tax not in assignments[marker]:
                assignments[marker][rec_tax] = list()
            if header not in unique_headers:
                unique_headers[header] = None
            # If fragments from the same parent query had identical lengths these would be overwritten anyway
            unique_headers[header] = {int(length): classified(refpkg=marker, taxon=rec_tax, length=int(length))}
        else:
            LOGGER.error("Bad line in classification table - no robust taxonomic classification:\n" +
                         '\t'.join(fields) + "\n")
            sys.exit(21)
    for header in unique_headers:
        if len(unique_headers[header]) > 1:
            dups.add(header)
        max_len = max(unique_headers[header])
        best_dat = unique_headers[header][max_len]
        assignments[best_dat.refpkg][best_dat.taxon].append(header)

    if dups:
        LOGGER.debug(str(len(dups)) + " fragments from identical parent sequences were removed post-classification.\n")
    return assignments


def read_classification_table(assignment_file) -> list:
    """
    Function for reading the tabular assignments file (currently classifications.tsv)
    Assumes column 2 is the TreeSAPP assignment and column 3 is the sequence header
    (leaving 1 for marker name and 4 for numerical abundance)

    :param assignment_file: Path to the file containing sequence phylogenetic origin and assignment
    :return: A list of lines that have been split by tabs into lists themselves
    """
    classified_lines = list()
    header = "Sample\tQuery\tMarker\tStart_pos\tEnd_pos\tTaxonomy\tAbundance\tiNode\tE-value\tLWR\tEvoDist\tDistances\n"

    try:
        assignments_handle = open(assignment_file, 'r')
    except IOError:
        LOGGER.error("Unable to open classification file '" + assignment_file + "' for reading.\n")
        sys.exit(21)

    header_line = assignments_handle.readline()
    if not header_line:
        LOGGER.error("Classification file '{}' is empty!\n".format(assignment_file))
        sys.exit(21)

    # This is the header line
    if not header_line.startswith(header.strip()):
        LOGGER.error("Header of assignments file is unexpected!\n")
        sys.exit(21)

    # First line in the table containing data
    line = assignments_handle.readline()
    n_fields = len(header_line.split("\t"))
    while line:
        fields = line.strip().split('\t')
        if len(fields) == n_fields:
            classified_lines.append(fields)
        else:
            LOGGER.error("Unable to parse line:\n" + str(line))
            sys.exit(21)
        line = assignments_handle.readline()
    assignments_handle.close()

    return classified_lines


def load_classified_sequences_from_assign_output(assign_output_dir: str, assigner_cls, refpkg_name=None) -> dict:
    """
    Reads the classified sequence names from Assigner.classification_tbl_name (e.g. classifications.tsv).
    These are converted into phylo_seq.PQuery instances and their 'seq' attributes are populated by the FASTA file
    holding sequence information of the classified queries.

    :param assign_output_dir: Path to a treesapp assign output directory
    :param assigner_cls: A treesapp.assign.Assigner
    :param refpkg_name: The prefix attribute of a ReferencePackage
    :return: A dictionary of PQuery instances indexed by their respective reference package
    """
    assigner_cls.final_output_dir = os.path.join(assign_output_dir, "final_outputs")
    classification_tbl = os.path.join(assigner_cls.final_output_dir, assigner_cls.classification_tbl_name)
    try:
        classified_seqs_fa = glob(os.path.join(assigner_cls.final_output_dir, "*_classified.faa")).pop()
    except IndexError:
        LOGGER.error("Classified sequences FASTA file was not found in output '{}'.\n".format(assign_output_dir))
        sys.exit(7)
    assigner_cls.sample_prefix = re.sub("_classified.faa", '', os.path.basename(classified_seqs_fa))

    # Read the lines from the classification table and convert assignments to PQuery instances
    pqueries = phylo_seq.assignments_to_pqueries(read_classification_table(classification_tbl))
    if refpkg_name:
        # Remove sequences that were not classified as refpkg_name
        if refpkg_name not in pqueries:
            LOGGER.warning("No queries were classified as '{}' in sample {}.\n"
                           "".format(refpkg_name, assigner_cls.sample_prefix))
            return {}
        else:
            for key_name in set(pqueries.keys()).difference({refpkg_name}):
                pqueries.pop(key_name)

    # Read the sequences into the 'seq' attribute of the pqueries
    tmp_dict = {}
    for pquery_list in pqueries.values():
        tmp_dict.update({pq.place_name: pq for pq in pquery_list})
    pquery_fasta = fasta.read_fasta_to_dict(classified_seqs_fa)
    if len(pquery_fasta) == 0:
        LOGGER.warning("Unable to read classified query sequences from FASTA file '{}'.\n".format(classified_seqs_fa))
    for pq_name, pq in tmp_dict.items():  # type: (str, phylo_seq.PQuery)
        pq.seq = pquery_fasta[pq_name]
    return pqueries


def best_discrete_matches(matches: list) -> list:
    """
    Function for finding the best alignment in a list of HmmMatch() objects
    The best match is based off of the full sequence score

    :param matches: A list of HmmMatch() objects
    :return: List of the best HmmMatch's
    """
    # Code currently only permits multi-domains of the same gene
    dropped_annotations = list()
    len_sorted_matches = sorted(matches, key=lambda x: x.end - x.start)
    i = 0
    orf = len_sorted_matches[0].orf
    while i + 1 < len(len_sorted_matches):
        j = i + 1
        a_match = len_sorted_matches[i]  # type HmmMatch
        while j < len(len_sorted_matches):
            b_match = len_sorted_matches[j]  # type HmmMatch
            if a_match.target_hmm != b_match.target_hmm:
                if hmmer_tbl_parser.detect_orientation(a_match.start, a_match.end,
                                                       b_match.start, b_match.end) != "satellite":
                    if a_match.full_score > b_match.full_score:
                        dropped_annotations.append(len_sorted_matches.pop(j))
                        j -= 1
                    else:
                        dropped_annotations.append(len_sorted_matches.pop(i))
                        j = len(len_sorted_matches)
                        i -= 1
            j += 1
        i += 1

    if len(len_sorted_matches) == 0:
        LOGGER.error("All alignments were discarded while deciding the best discrete HMM-match.\n")
        sys.exit(3)

    LOGGER.debug("HMM search annotations for " + orf +
                 ":\n\tRetained\t" +
                 ', '.join([match.target_hmm +
                            " (%d-%d)" % (match.start, match.end) for match in len_sorted_matches]) +
                 "\n\tDropped\t\t" +
                 ', '.join([match.target_hmm +
                            " (%d-%d)" % (match.start, match.end) for match in dropped_annotations]) + "\n")
    return len_sorted_matches


def parse_domain_tables(thresholds, hmm_domtbl_files: dict) -> dict:
    """
    Parses HMMER domain tables using predetermined thresholds

    :param thresholds: A namedtuple instance: namedtuple("thresholds", "max_e max_ie min_acc min_score perc_aligned")
    :param hmm_domtbl_files: A list of domain table files written by hmmsearch
    :return: Dictionary of HmmMatch objects indexed by their reference package and/or HMM name
    """
    # Check if the HMM filtering thresholds have been set
    LOGGER.info("Parsing HMMER domain tables for high-quality matches... ")

    search_stats = hmmer_tbl_parser.HmmSearchStats()
    hmm_matches = dict()
    orf_gene_map = dict()
    optional_matches = list()

    # TODO: Capture multimatches across multiple domain table files
    for r_q, domtbl_file in hmm_domtbl_files.items():
        _prefix, reference = r_q
        domain_table = hmmer_tbl_parser.DomainTableParser(domtbl_file)
        domain_table.read_domtbl_lines()
        distinct_hits = hmmer_tbl_parser.format_split_alignments(domain_table, search_stats)
        purified_hits = hmmer_tbl_parser.filter_poor_hits(thresholds, distinct_hits, search_stats)
        complete_hits = hmmer_tbl_parser.filter_incomplete_hits(thresholds, purified_hits, search_stats)
        hmmer_tbl_parser.renumber_multi_matches(complete_hits)

        for match in complete_hits:
            match.genome = reference
            if match.orf not in orf_gene_map:
                orf_gene_map[match.orf] = dict()
            try:
                orf_gene_map[match.orf][match.target_hmm].append(match)
            except KeyError:
                orf_gene_map[match.orf][match.target_hmm] = [match]
            if match.target_hmm not in hmm_matches.keys():
                hmm_matches[match.target_hmm] = list()
    search_stats.num_dropped()
    for orf in orf_gene_map:
        if len(orf_gene_map[orf]) == 1:
            for target_hmm in orf_gene_map[orf]:
                for match in orf_gene_map[orf][target_hmm]:
                    hmm_matches[target_hmm].append(match)
                    search_stats.seqs_identified += 1
        else:
            search_stats.multi_alignments += 1
            # Remove all the overlapping domains - there can only be one highlander
            for target_hmm in orf_gene_map[orf]:
                optional_matches += orf_gene_map[orf][target_hmm]
            retained = 0
            for discrete_match in best_discrete_matches(optional_matches):
                hmm_matches[discrete_match.target_hmm].append(discrete_match)
                retained += 1
            search_stats.dropped += (len(optional_matches) - retained)
            search_stats.seqs_identified += retained
            optional_matches.clear()

    LOGGER.info("done.\n")

    alignment_stat_string = search_stats.summarize()

    if search_stats.seqs_identified == 0 and search_stats.dropped == 0:
        LOGGER.warning("No alignments found! TreeSAPP is exiting now.\n")
        sys.exit(0)
    if search_stats.seqs_identified == 0 and search_stats.dropped > 0:
        LOGGER.warning("No alignments (" + str(search_stats.seqs_identified) + '/' + str(search_stats.dropped) +
                       ") met the quality cut-offs! TreeSAPP is exiting now.\n")
        alignment_stat_string += "\tPoor quality alignments:\t" + str(search_stats.bad) + "\n"
        alignment_stat_string += "\tShort alignments:\t" + str(search_stats.short) + "\n"
        LOGGER.debug(alignment_stat_string)
        sys.exit(0)

    alignment_stat_string += "\n\tNumber of markers identified:\n"
    for marker in sorted(hmm_matches):
        alignment_stat_string += "\t\t" + marker + "\t" + str(len(hmm_matches[marker])) + "\n"
        # For debugging:
        # for match in hmm_matches[marker]:
        #     match.print_info()

    LOGGER.debug(alignment_stat_string)
    return hmm_matches


def read_colours_file(annotation_file: str, refpkg_name: str) -> dict:
    """
    Read annotation data from 'annotation_file' and store it in marker_subgroups under the appropriate
    marker and data_type.

    :param annotation_file: Path to an iTOL-compatible annotation file (e.g. colours_styles_file.txt)
    :param refpkg_name: Name of the reference package
    :return: A dictionary of lists where each list is populated by tuples with start and end leaves
    """
    try:
        style_handler = open(annotation_file, 'r')
    except IOError:
        LOGGER.error("Unable to open " + annotation_file + " for reading!\n")
        sys.exit(5)

    clusters = dict()
    field_sep = ''

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
                LOGGER.error("Unknown separator used in " + annotation_file + ": " + header_fields[1] + "\n")
                sys.exit(5)
        line = style_handler.readline()

    # For RGB
    range_line_rgb = re.compile(r"^(\d+)_" + re.escape(refpkg_name) + r"\|(\d+)_" + re.escape(refpkg_name) +
                                re.escape(field_sep) + "range" + re.escape(field_sep) + r".*\)" + re.escape(field_sep) +
                                "(.*)$")
    single_node_rgb = re.compile(r"^(\d+)_" + re.escape(refpkg_name) + re.escape(field_sep) +
                                 "range" + re.escape(field_sep) +
                                 r".*\)" + re.escape(field_sep) +
                                 "(.*)$")
    internal_node_rgb = re.compile(r"^(\d+)" + re.escape(field_sep) +
                                   "range" + re.escape(field_sep) +
                                   r".*\)" + re.escape(field_sep) +
                                   "(.*)$")

    # For hexadecimal
    range_line = re.compile(r"^(\d+)_" + re.escape(refpkg_name) + r"\|(\d+)_" +
                            re.escape(refpkg_name) + re.escape(field_sep) +
                            "range" + re.escape(field_sep) +
                            "#[0-9A-Za-z]{6}" + re.escape(field_sep) +
                            "(.*)$")
    single_node = re.compile(r"^(\d+)_" + re.escape(refpkg_name) + re.escape(field_sep) +
                             "range" + re.escape(field_sep) +
                             "#[0-9A-Za-z]{6}" + re.escape(field_sep) +
                             "(.*)$")
    internal_node = re.compile(r"^(\d+)" + re.escape(field_sep) +
                               "range" + re.escape(field_sep) +
                               "#[0-9A-Za-z]{6}" + re.escape(field_sep) +
                               "(.*)$")

    single_name = re.compile(r"^([\w\d|.]+)" + re.escape(field_sep) +
                             "range" + re.escape(field_sep) + r"[#(), \w]+" +
                             re.escape(field_sep) + "(.*)$")

    # Begin parsing the data from 4 columns
    line = style_handler.readline().strip()
    while line:
        if range_line.match(line):
            style_data = range_line.match(line)
            start, end, description = style_data.groups()
            label_node = True
        elif range_line_rgb.match(line):
            style_data = range_line_rgb.match(line)
            start, end, description = style_data.groups()
            label_node = True
        elif single_node.match(line):
            style_data = single_node.match(line)
            start, end, description = style_data.group(1), style_data.group(1), style_data.group(2)
            label_node = True
        elif single_node_rgb.match(line):
            style_data = single_node_rgb.match(line)
            start, end, description = style_data.group(1), style_data.group(1), style_data.group(2)
            label_node = True
        elif internal_node.match(line):
            style_data = internal_node.match(line)
            start, end, description = style_data.group(1), style_data.group(1), style_data.group(2)
            label_node = False
        elif internal_node_rgb.match(line):
            style_data = internal_node_rgb.match(line)
            start, end, description = style_data.group(1), style_data.group(1), style_data.group(2)
            label_node = False
        elif single_name.match(line):
            style_data = single_name.match(line)
            start, end, description = style_data.group(1), style_data.group(1), style_data.group(2)
            label_node = False
        else:
            LOGGER.error("Unrecognized line formatting in '{}':\n{}\n".format(annotation_file, line))
            sys.exit(5)

        description = style_data.groups()[-1]
        if description not in clusters.keys():
            clusters[description] = list()

        if label_node:
            clusters[description].append((start + "_" + refpkg_name, end + "_" + refpkg_name))
        else:
            clusters[description].append((start, end))

        line = style_handler.readline().strip()

    style_handler.close()

    LOGGER.debug("\tParsed {} clades from '{}'\n".format(len(clusters), annotation_file))

    return clusters


def read_phylip_to_dict(phylip_input: str) -> dict:
    header_dict = dict()
    tmp_seq_dict = dict()
    seq_dict = dict()
    x = 0

    try:
        phylip = open(phylip_input, 'r')
    except IOError:
        LOGGER.error("Unable to open the Phylip file (" + phylip_input + ") provided for reading!\n")
        sys.exit(5)

    line = phylip.readline()
    try:
        num_sequences, aln_length = line.strip().split(' ')
        num_sequences = int(num_sequences)
        aln_length = int(aln_length)
    except ValueError:
        LOGGER.error("Phylip file is not formatted correctly!\n" +
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
            LOGGER.error("Unexpected line in Phylip file:\n" + line + "\n")
            sys.exit(5)
        line = phylip.readline()

        if x > num_sequences:
            LOGGER.error("Accumulator has exceeded the number of sequences in the file (according to header)!\n")
            sys.exit(5)

    # Check that the alignment length matches that in the header line
    if num_sequences != len(tmp_seq_dict):
        LOGGER.error("Number of lines declared in Phylip header ({}) does not match number of sequences parsed ({})!\n"
                     "".format(num_sequences, len(tmp_seq_dict)))
        sys.exit(5)

    x = 0
    while x < num_sequences - 1:
        if len(tmp_seq_dict[x]) != aln_length:
            LOGGER.error("{} sequence length exceeds the stated multiple alignment length (according to header)!\n"
                         "sequence length = {}, alignment length = {}\n"
                         "".format(header_dict[x], len(tmp_seq_dict[x]), aln_length))
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
        LOGGER.error("Unable to open " + sto_file + " for reading!\n")
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
                LOGGER.error("Unexpected line format in " + sto_file + ":\n" + line + "\n")
                sys.exit(3)

            if seq_name not in seq_dict:
                seq_dict[seq_name] = ""
            seq_dict[seq_name] += re.sub(r'\.', '-', sequence.upper())
        line = sto_handler.readline()

    return seq_dict


def check_seq_name_integer_compatibility(seq_dict: dict) -> (dict, int):
    # Parse the MSA dict and ensure headers are integer-compatible
    multi_align = {}
    n_msa_refs = 0
    for seq_name, seq in seq_dict.items():
        try:
            if int(seq_name) > 0:
                n_msa_refs += 1
        except ValueError:
            if re.match(r"^_\d+", seq_name):
                leaf_num = re.sub("^_", '-', seq_name)
            # The section of regular expresion after '_' needs to match denominator and refpkg names
            elif re.match(r"^\d+_\w{2,10}$", seq_name):
                leaf_num = seq_name.split('_')[0]
            else:
                return {seq_name: ""}, -1
            if int(leaf_num) > 0:
                n_msa_refs += 1
        multi_align[seq_name] = seq
    return multi_align, n_msa_refs


def validate_alignment_trimming(msa_files: list, unique_ref_headers: set,
                                queries_mapped=False, min_seq_length=30) -> (dict, list, str):
    """
    Parse a list of multiple sequence alignment (MSA) files and determine whether the multiple alignment:
        1. is shorter than the min_seq_length (30 by default)
        2. is missing any reference sequences
    The number of query sequences discarded - these may have been added by hmmalign or PaPaRa - is returned via a string

    NOTE: Initially designed for sequence records with numeric names (e.g. >488) but accommodates other TreeSAPP formats

    :param msa_files: A list of either Phylip or FASTA formatted MSA files
    :param unique_ref_headers: A set of all headers that were in the untrimmed MSA
    :param queries_mapped: Boolean indicating whether sequences should be present in addition to reference sequences.
           While query sequences _could_ be identified as any that are not in unique_ref_headers,
           queries have names that are negative integers for more rapid and scalable identification
    :param min_seq_length: Optional minimum unaligned (no '-'s) length a sequence must exceed to be retained
    :return: 1. Dictionary indexed by MSA file name mapping to FASTA-dictionaries
             2. A string mapping the number of query sequences removed from each MSA file
             3. A string describing the number of sequences discarded
    """
    discarded_seqs_string = ""
    successful_multiple_alignments = dict()
    failed_multiple_alignments = list()
    n_refs = len(unique_ref_headers)
    for multi_align_file in msa_files:
        filtered_multi_align = dict()
        discarded_seqs = list()
        num_queries_retained = 0
        n_retained_refs = 0
        f_ext = multi_align_file.split('.')[-1]

        # Read the multiple alignment file
        if re.search("phy", f_ext):  # File is in Phylip format
            seq_dict = read_phylip_to_dict(multi_align_file)
        elif re.match("^f", f_ext):  # This is meant to match all fasta extensions
            seq_dict = fasta.read_fasta_to_dict(multi_align_file)
        elif f_ext == "mfa":  # This is meant to match a multiple alignment in FASTA format
            seq_dict = fasta.read_fasta_to_dict(multi_align_file)
        else:
            LOGGER.error("Unable to detect file format of " + multi_align_file + ".\n")
            sys.exit(13)

        multi_align, n_msa_refs = check_seq_name_integer_compatibility(seq_dict)
        if n_msa_refs < 0:
            LOGGER.error("Unexpected sequence name ('{}') detected in {}.\n"
                         "".format(multi_align.popitem()[0], multi_align_file))
            sys.exit(13)
        if len(multi_align) == 0:
            LOGGER.warning("No sequences were read from {}. "
                           "The untrimmed alignment will be used instead.\n".format(multi_align_file))
            failed_multiple_alignments.append(multi_align_file)
            continue
        # The numeric identifiers make it easy to maintain order in the Phylip file by a numerical sort
        for seq_name in sorted(multi_align, key=lambda x: int(x.split('_')[0])):
            seq_dummy = re.sub('-', '', multi_align[seq_name])
            if len(seq_dummy) < min_seq_length:
                discarded_seqs.append(seq_name)
            else:
                filtered_multi_align[seq_name] = multi_align[seq_name]
                # The negative integers indicate this is a query sequence
                if seq_name[0] == '-':
                    num_queries_retained += 1
                else:
                    n_retained_refs += 1
        discarded_seqs_string += "\n\t\t" + multi_align_file + " = " + str(len(discarded_seqs))
        if len(discarded_seqs) == len(multi_align.keys()):
            # Throw an error if the final trimmed alignment is shorter than min_seq_length, and therefore empty
            LOGGER.warning("Multiple sequence alignment in {} is shorter than minimum sequence length threshold ({})."
                           "\nThe untrimmed MSA will be used instead.\n".format(multi_align_file, min_seq_length))
            failed_multiple_alignments.append(multi_align_file)
        elif n_refs > n_msa_refs:
            # Testing whether there were more sequences in the untrimmed alignment than the trimmed one
            LOGGER.warning("Reference sequences in " + multi_align_file + " were removed during alignment trimming " +
                           "suggesting either truncated sequences or the initial reference alignment was terrible.\n" +
                           "The untrimmed alignment will be used instead.\n")
            failed_multiple_alignments.append(multi_align_file)
        elif n_refs > n_retained_refs:
            LOGGER.warning("Reference sequences shorter than the minimum character length ({})"
                           " in {} were removed after alignment trimming.\n".format(min_seq_length, multi_align_file) +
                           "The untrimmed alignment will be used instead.\n")
            failed_multiple_alignments.append(multi_align_file)
        # Ensure that there is at least 1 query sequence retained after trimming the multiple alignment
        elif queries_mapped and num_queries_retained == 0:
            LOGGER.warning("No query sequences in " + multi_align_file + " were retained after trimming.\n")
        else:
            successful_multiple_alignments[multi_align_file] = filtered_multi_align

        if multi_align_file in successful_multiple_alignments:
            discarded_seqs_string += " (retained)"
        else:
            discarded_seqs_string += " (removed)"

    return successful_multiple_alignments, failed_multiple_alignments, discarded_seqs_string


def read_annotation_mapping_file(annot_map_file: str) -> dict:
    """
    Used for reading a file mapping the reference package name to all true positive orthologs in the query input
    The first column is the ReferencePackage.prefix.
    The second column is the ortholog name used by the database.
    The third column is the name of a true positive.

    :param annot_map_file: Path to a tab-delimited file
    :return: A dictionary containing database sequence names mapped to their respective reference package names
    """
    annot_map = dict()
    try:
        annot_map_handler = open(annot_map_file)
    except IOError:
        LOGGER.error("Unable to open annotation file '{}' for reading!\n".format(annot_map_file))
        sys.exit(3)

    # Assuming the first column is the reference package name and the second is the database annotation name
    n = 0
    for line in annot_map_handler:
        n += 1
        if line[0] == '#':
            continue
        elif not line:
            continue
        else:
            try:
                refpkg_name, og, query_name = line.strip().split("\t")
            except ValueError:
                LOGGER.error("Unexpected number of fields on line {} in {}!\n".format(n, annot_map_file) +
                             "File must have the reference package name and the database name in"
                             " the first two columns, respectively. Any number of columns can follow.\n")
                sys.exit(9)
            if query_name not in annot_map:
                annot_map[query_name] = set()
            annot_map[query_name].add(refpkg_name)

    annot_map_handler.close()
    return annot_map


def read_phenotypes(phenotypes_file: str, comment_char='#') -> dict:
    taxa_phenotype_map = {}
    try:
        file_handler = open(phenotypes_file, 'r')
    except IOError:
        LOGGER.error("Unable to open taxa-phenotype table '{}' for reading.\n".format(phenotypes_file))
        sys.exit(7)

    for line in file_handler:
        if not line or line[0] == comment_char:
            continue
        if line.find(comment_char) >= 0:
            line = line[:line.find(comment_char)]
        try:
            taxon_name, phenotype = line.rstrip().split("\t")
        except ValueError:
            LOGGER.error("Unable to parse line in {}:\n{}\n".format(phenotypes_file, line))
            sys.exit(9)
        if taxon_name in taxa_phenotype_map:
            LOGGER.warning("Taxon '{}' found in {} multiple times and will be overwritten.\n"
                           "".format(taxon_name, phenotypes_file))
        taxa_phenotype_map[taxon_name.strip()] = phenotype.strip()

    file_handler.close()

    return taxa_phenotype_map


def read_lineage_map(lineage_table: str) -> dict:
    if not os.path.isfile(lineage_table):
        LOGGER.error("lineage mapping table '{}' does not exist.\n".format(lineage_table))
        sys.exit(5)

    sep = utilities.get_field_delimiter(lineage_table)
    lines = [line.strip().split(sep) for line in utilities.get_file_lines(lineage_table)]
    if re.search(r"^#", lines[0][0]):
        lines.pop(0)
    lineage_map = {lin[0]: lin[1] for lin in lines}

    return lineage_map


def read_lineage_ids(lineage_table: str) -> dict:
    lineage_map = dict()
    if not os.path.isfile(lineage_table):
        LOGGER.error("lineage_ids table '{}' does not exist.\n".format(lineage_table))
        sys.exit(5)

    lin_id_lines = utilities.get_file_lines(lineage_table)
    counter = 1
    for line in lin_id_lines:
        try:
            num, desc, lineage = line.strip().split("\t")
        except (IndexError, ValueError):
            LOGGER.error("Line {} in lineage table '{}' is misformatted:\n{}\n".format(counter, lineage_table, line))
            sys.exit(15)
        lineage_map[num] = desc + "\t" + lineage
        counter += 1

    return lineage_map
