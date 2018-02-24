#!/usr/bin/env python3

__author__ = 'Connor Morgan-Lang'

import argparse
import sys
import pygtrie
import re
import os
import inspect

cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile(inspect.currentframe()))[0]))
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)
sys.path.insert(0, cmd_folder + os.sep + ".." + os.sep)
from treesapp import format_read_fasta
from create_treesapp_ref_data import get_lineage, get_header_format, \
    return_sequence_info_groups, map_good_headers_to_ugly, get_headers
from utilities import clean_lineage_string
from external_command_interface import setup_progress_bar

rank_depth_map = {0: "Cellular organisms", 1: "Kingdom",
                  2: "Phylum", 3: "Class", 4: "Order",
                  5: "Family", 6: "Genus", 7: "Species",
                  8: "Strain"}


def get_arguments_():
    parser = argparse.ArgumentParser(add_help=False)
    required_args = parser.add_argument_group("Required arguments")
    required_args.add_argument("-a", "--assignment_table",
                               help="Tabular file containing each sequence's phylogenetic origin"
                                    " and its TreeSAPP assignment.",
                               required=True)
    required_args.add_argument("-t", "--tax_ids_table",
                               help="tax_ids table for the TreeSAPP marker under investigation. "
                                    "Found in data/tree_data/tax_ids*.txt.",
                               required=True,
                               nargs='+')
    required_args.add_argument('-i', '--fasta_input',
                               help='Your sequence input file (for TreeSAPP) in FASTA format',
                               required=True)

    optopt = parser.add_argument_group("Optional options")
    optopt.add_argument('-o', '--output', default='./output/', required=False,
                        help='output directory [DEFAULT = ./output/]')
    optopt.add_argument('-m', '--molecule',
                        help='the type of input sequences (prot = Protein [DEFAULT]; dna = Nucleotide; rrna = rRNA)',
                        default='prot',
                        choices=['prot', 'dna', 'rrna'])
    optopt.add_argument("-h", "--help",
                        action="help",
                        help="Show this help message and exit")

    miscellaneous_opts = parser.add_argument_group("Miscellaneous options")
    miscellaneous_opts.add_argument('--overwrite', action='store_true', default=False,
                                    help='overwrites previously processed output folders')
    miscellaneous_opts.add_argument('-v', '--verbose', action='store_true', default=False,
                                    help='Prints a more verbose runtime log')

    args = parser.parse_args()

    if sys.version_info > (2, 9):
        args.py_version = 3
    else:
        args.py_version = 2

    return args


def read_table(assignment_file):
    """
    Function for reading the tabular assignments file (currently marker_contig_map.tsv)
    Assumes column 2 is the TreeSAPP assignment and column 3 is the sequence header
    (leaving 1 for marker name and 4 for numerical abundance)
    :param assignment_file: Path to the file containing sequence phylogenetic origin and assignment
    :return: dictionary whose keys are phylogenetic origin and values are lists of TreeSAPP assignments
    """
    assignments = dict()
    n_classified = 0
    assignments_handle = open(assignment_file, 'r')
    # This is the header line
    if not re.match("^Query\tMarker\tTaxonomy\tConfident_Taxonomy\tAbundance$", assignments_handle.readline()):
        sys.stderr.write("ERROR: header of assignments file is unexpected!\n")
        raise AssertionError

    # First line in the table containing data
    line = assignments_handle.readline()
    while line:
        fields = line.strip().split('\t')
        try:
            header, marker, classified, rob_class, abundance = fields
            if marker and rob_class:
                n_classified += 1
                if marker not in assignments:
                    assignments[marker] = dict()
                if rob_class not in assignments[marker]:
                    assignments[marker][rob_class] = list()
                assignments[marker][rob_class].append(header)
        except ValueError:
            sys.stderr.write("ERROR: Unable to parse line:")
            sys.stderr.write(str(line))
            sys.exit(1)
        line = assignments_handle.readline()

    assignments_handle.close()
    return assignments, n_classified


def read_intermediate_assignments(args, tmp_file):
    if args.verbose:
        sys.stderr.write("NOTIFICATION: Reading " + tmp_file + " for saved intermediate data... ")
    assignments = dict()
    n_saved = 0
    try:
        file_handler = open(tmp_file, 'r')
    except OSError:
        sys.stderr.write("ERROR: Unable to open " + tmp_file + " for reading!\n")
        raise OSError
    line = file_handler.readline()
    while line:
        line = line.strip()
        try:
            marker, ref, queries = line.split('\t')
        except ValueError:
            sys.stderr.write("\tWARNING: " + tmp_file + " is incorrectly formatted so regenerating instead.\n")
            sys.stderr.write("Violating line with " + str(len(line.split('\t'))) + " elements:\n" + line + "\n")
            return assignments, n_saved
        queries_list = queries.split(',')
        n_saved += len(queries_list)

        if marker not in assignments:
            assignments[marker] = dict()
        assignments[marker][ref] = queries_list
        line = file_handler.readline()

    file_handler.close()
    if args.verbose:
        sys.stderr.write("done.\n")
        sys.stderr.flush()
    return assignments, n_saved


def write_intermediate_assignments(args, assignments):
    tmp_file = args.output + os.sep + "tmp_clade_exclusion_assignments.tsv"
    if args.verbose:
        sys.stderr.write("Saving intermediate data... ")
    try:
        file_handler = open(tmp_file, 'w')
    except OSError:
        sys.stderr.write("ERROR: Unable to open " + tmp_file + " for writing!\n")
        raise OSError
    assignments_string = ""
    for marker in assignments.keys():
        for ref in assignments[marker].keys():
            queries = assignments[marker][ref]
            if len(queries) == 0:
                sys.stderr.write("\tWARNING: No lineage information was downloaded for queries assigned to " +
                                 ref + "\n")
            else:
                cleaned_ref = clean_lineage_string(ref)
                if not cleaned_ref:
                    if args.verbose:
                        sys.stderr.write("\tWARNING: reference assignment is empty after cleaning the lineage:\n" +
                                         ref + "\n")
                else:
                    assignments_string += marker + "\t" + cleaned_ref + "\t"
                    assignments_string += ','.join(queries) + "\n"

    file_handler.write(assignments_string)
    file_handler.close()
    if args.verbose:
        sys.stderr.write("done.\n")
        sys.stderr.flush()
    return


def all_possible_assignments(args, tax_ids_file):
    taxonomic_tree = pygtrie.StringTrie(separator='; ')
    if os.path.exists(tax_ids_file):
        file_name = os.path.basename(tax_ids_file)
        if re.match("^tax_ids_(.*).txt", file_name):
            marker = re.match("^tax_ids_(.*).txt", file_name).group(1)
        else:
            sys.stderr.write("ERROR: Format of tax_ids file (" + tax_ids_file +
                             ") is unexpected. Unable to parse marker name! Exiting...\n")
            sys.exit(7)
    else:
        raise IOError("File doesn't exist: " + tax_ids_file + "\n")
    try:
        if args.py_version == 3:
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
        if len(fields) == 3:
            number, translation, lineage = fields
            lineage = clean_lineage_string(lineage)
        else:
            sys.stderr.write("ValueError: Unexpected number of fields in " + tax_ids_file +
                             ".\nInvoked .split(\'\\t\') on line " + str(line))
            raise ValueError

        i = 0
        ranks = len(lineage)
        while i < len(lineage):
            taxonomic_tree["; ".join(lineage.split("; ")[:ranks - i])] = True
            i += 1

    cog_tax_ids.close()
    return taxonomic_tree, marker


def identify_excluded_clade(args, assignment_dict, trie, marker):
    """
    Using the taxonomic information from the sequence headers and the lineages of the reference sequence,
    this function determines the rank at which each sequence's clade is excluded.
    These data are returned and sorted in the form of a dictionary.
    :param: assignment_dict:
    :param: trie: A pygtrie.StringTrie object containing all reference sequence lineages
    :return: rank_assigned_dict; key is rank, values are dictionaries with assigned (reference) lineage as key and
    tuples of (optimal assignment, actual assignment) as values.
    E.g. {"Phylum": {"Proteobacteria": ("Proteobacteria", "Proteobacteria; Alphaproteobacteria")}}
    """
    rank_assigned_dict = dict()
    for depth in rank_depth_map:
        rank_assigned_dict[rank_depth_map[depth]] = list()
    if args.verbose:
        sys.stdout.write("Unique taxonomies sequences were assigned = " + str(len(assignment_dict[marker].keys())) + "\n")
    for ref_lineage in assignment_dict[marker]:
        if args.verbose:
            sys.stdout.write("Reference lineage: " + ref_lineage + "\n")
        for query_lineage in assignment_dict[marker][ref_lineage]:
            if query_lineage == ref_lineage:
                if args.verbose:
                    sys.stdout.write("\tQuery lineage: " + query_lineage + ", ")
                    sys.stdout.write("Optimal lineage: " + ref_lineage + "\n")
                    sys.stdout.flush()
            # While the query_lineage contains clades which are not in the reference trie,
            # remove the taxonomic rank and traverse again. Complexity: O(ln(n))
            contained_taxonomy = query_lineage
            while not trie.__contains__(contained_taxonomy) and len(contained_taxonomy.split('; ')) > 1:
                contained_taxonomy = "; ".join(contained_taxonomy.split('; ')[:-1])
            if args.verbose and contained_taxonomy != ref_lineage:
                sys.stdout.write("\tQuery lineage: " + query_lineage + ", ")
                sys.stdout.write("Optimal lineage: " + contained_taxonomy + "\n")
            if len(contained_taxonomy.split("; ")) <= 7:
                rank_assigned_dict[rank_depth_map[len(contained_taxonomy.split("; ")) + 1]].append(
                    {ref_lineage: (contained_taxonomy, query_lineage)})
            else:
                # TODO: Fix the handling of this. Currently printed for strains when it shouldn't be
                sys.stderr.write("\tWARNING: number of ranks in lineage '" + contained_taxonomy + "' is ridiculous.\n")
                sys.stderr.write("\tThis sequence will be removed from clade exclusion calculations\n")
    return rank_assigned_dict


def determine_offset(classified, optimal):
    # Figure out which taxonomic lineage is longer
    offset = 0
    while classified != optimal and offset < 7:
        offset += 1
        classified_ranks = classified.split("; ")
        optimal_ranks = optimal.split("; ")
        if len(classified_ranks) > len(optimal_ranks):
            classified = "; ".join(classified_ranks[:-1])
        elif len(classified_ranks) < len(optimal_ranks):
            optimal = "; ".join(optimal_ranks[:-1])
        else:
            optimal = "; ".join(optimal_ranks[:-1])
            classified = "; ".join(classified_ranks[:-1])
    return offset


def summarize_taxonomic_diversity(query_lineages):
    """
    Function for summarizing the taxonomic diversity of a reference dataset by rank
    :param query_lineages: A list of lineages represented in the test dataset sequences
    :return:
    """
    # Not really interested in Cellular Organisms or Strains.
    rank_depth_map = {0: "Kingdoms", 1: "Phyla", 2: "Classes", 3: "Orders", 4: "Families", 5: "Genera", 6: "Species"}
    taxa_counts = dict()
    for depth in rank_depth_map:
        name = rank_depth_map[depth]
        taxa_counts[name] = set()
    for lineage in sorted(query_lineages):
        position = 0
        taxa = lineage.split('; ')
        while position < len(taxa) and position < 7:
            taxa_counts[rank_depth_map[position]].add(taxa[position])
            position += 1
    sys.stdout.write("Number of unique lineages:\n")
    for depth in rank_depth_map:
        rank = rank_depth_map[depth]
        buffer = " "
        while len(rank) + len(str(len(taxa_counts[rank]))) + len(buffer) < 12:
            buffer += ' '
        sys.stdout.write("\t" + rank + buffer + str(len(taxa_counts[rank])) + "\n")
    sys.stdout.flush()

    return


def determine_specificity(rank_assigned_dict, marker, clade_exclusion_strings):
    """
    Correct if: optimal_assignment == query_lineage
    Correct if:
    :param rank_assigned_dict:
    :param marker: Name of the marker currently being evaluated (e.g., nifHc1, mcrA)
    :param clade_exclusion_strings: A list of strings that is appended to. Finally used for writing the tabular output.
    :return:
    """
    sys.stdout.write("Clade-level specificities for " + marker + ":\n")
    sys.stdout.write("\tRank\tTotal Evaluated\tCorrect\tD=1\tD=2\tD=3\tD=4\tD=5\tD=6\tD=7\n")
    clade_exclusion_tabular_string = ""
    clades_tested = list()
    for depth in sorted(rank_depth_map):
        rank = rank_depth_map[depth]
        if rank == "Cellular organisms":
            continue
        correct = 0
        incorrect = 0
        taxonomic_distance = dict()
        for i in range(0, 8):
            taxonomic_distance[i] = 0
        sys.stdout.write("\t" + rank + "\t")
        if len(rank_assigned_dict[rank]) == 0:
            sys.stdout.write("0\t\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n")
        else:
            for assignments in rank_assigned_dict[rank]:
                # if rank == "Family":
                #     print(assignments)
                for classified in assignments:
                    if classified.split("; ")[0] == "Cellular organisms":
                        sys.stderr.write("ERROR: lineage string cleaning has gone awry somewhere. "
                                         "The root rank should be a Kingdom (e.g. Bacteria or Archaea) but nope.\n")
                        sys.exit(9)
                    optimal, query = assignments[classified]
                    clades_tested.append(query)
                    if optimal == classified:
                        offset = 0
                        correct += 1
                    else:
                        offset = determine_offset(classified, optimal)
                        incorrect += 1
                    taxonomic_distance[offset] += 1
            rank_total = correct + incorrect
            sys.stdout.write(str(rank_total) + "\t\t")

            for dist in taxonomic_distance:
                if taxonomic_distance[dist] > 0:
                    taxonomic_distance[dist] = str(round(float((taxonomic_distance[dist]*100)/rank_total), 1))
                else:
                    taxonomic_distance[dist] = str(0.0)
                clade_exclusion_tabular_string += marker + "\t" + rank + "\t"
                clade_exclusion_tabular_string += str(rank_total) + "\t"
                clade_exclusion_tabular_string += str(dist) + "\t" + str(taxonomic_distance[dist]) + "\t"
                clade_exclusion_strings.append(clade_exclusion_tabular_string)
                clade_exclusion_tabular_string = ""
            sys.stdout.write('\t'.join(taxonomic_distance.values()) + "\n")
    summarize_taxonomic_diversity(clades_tested)

    return clade_exclusion_strings


def determine_containment(args, marker, rank_assigned_dict):
    """
    Determines the accuracy of sequence classifications of all sequences contained at different taxonomic ranks
    :param args:
    :param marker: Name of the marker currently being evaluated (e.g., nifHc1, mcrA)
    :param rank_assigned_dict: key is rank, values are dictionaries with assigned (reference) lineage as key and
    tuples of (optimal assignment, actual assignment) as values.
    E.g. {"Phylum": {"Proteobacteria": ("Proteobacteria", "Proteobacteria; Alphaproteobacteria")}}
    """
    sys.stdout.write("Clade-level containments for " + marker + ":\n")
    sys.stdout.write("\tRank\t\tTotal Evaluated\tCorrect (%)\tToo Shallow (%)\n")
    # Set up collection for this analysis
    all_assignments = list()
    for rank in rank_assigned_dict:
        if len(rank_assigned_dict[rank]) > 0:
            all_assignments += rank_assigned_dict[rank]
    # Begin parsing through the depths
    for depth in sorted(rank_depth_map):
        rank = rank_depth_map[depth]
        depth = depth - 1
        if rank in ["Cellular organisms", "Species", "Strain"]:
            continue
        correct = 0
        incorrect = 0
        too_shallow = 0
        rank_total = 0
        incorrect_assignments = dict()
        sys.stdout.write("\t" + rank + "\t\t")
        for assignments in all_assignments:
            for classified in assignments:
                rank_total += 1
                optimal, query = assignments[classified]
                classified_lineage = classified.split("; ")
                query_lineage = query.split("; ")
                if len(classified_lineage) > depth:
                    # print("\nReference:\t", classified)
                    # print("Query:\t\t", query)
                    # print("Optimal:\t", optimal)
                    try:
                        if classified_lineage[depth] == query_lineage[depth]:
                            correct += 1
                            # print("Correct")
                        else:
                            # print("Incorrect - 1")
                            if query not in incorrect_assignments.keys():
                                incorrect_assignments[query] = 0
                            incorrect_assignments[query] += 1
                            incorrect += 1
                    except IndexError:
                        too_shallow += 1
                        # This indicates TreeSAPP placed the sequence too deep in the tree according to its NCBI lineage
                        # sys.stderr.write(str(depth) + " " + str(classified_lineage) + " " + str(query_lineage) + "\n")
                elif len(optimal.split("; ")) > depth:
                    # print("\nReference:\t", classified)
                    # print("Query:\t\t", query)
                    # print("Optimal:\t", optimal)
                    # print("Incorrect - 2")
                    if query not in incorrect_assignments.keys():
                        incorrect_assignments[query] = 0
                    incorrect_assignments[query] += 1
                    incorrect += 1
                else:
                    too_shallow += 1

        # if incorrect == correct + incorrect:
        #     print(incorrect_assignments)

        sys.stdout.write(str(rank_total) + "\t\t")
        sys.stdout.write(str(round((correct * 100) / (correct + incorrect), 0)) + "\t\t")
        percentage_too_shallow = round(float((too_shallow * 100)/rank_total), 1)
        sys.stdout.write(str(percentage_too_shallow) + "\n")
    sys.stdout.write("\n")

    return


def clean_classification_names(assignments):
    """
    Removes certain ranks, taxonomic super-groups, etc. so strings are more comparable
    :param assignments:
    :return:
    """
    cleaned_assignments = dict()
    for marker in assignments.keys():
        cleaned_assignments[marker] = dict()
        for ref, queries in assignments[marker].items():
            ref = clean_lineage_string(ref)
            cleaned_assignments[marker][ref] = list()
            for query in queries:
                try:
                    query = re.sub('_', "; ", query)
                except TypeError:
                    sys.stderr.write("ERROR: Unexpected type of lineage string (" + str(type(query)) + "):\n")
                    sys.stderr.write(str(query) + "\n")
                    raise TypeError
                query = clean_lineage_string(query)
                cleaned_assignments[marker][ref].append(query)
    return cleaned_assignments


def map_full_headers(fasta_headers, header_map, assignments, molecule_type):
    """
    Since the headers used throughout the TreeSAPP pipeline are truncated,
    we read the FASTA file and use those headers instead of their short version
    in case valuable information was discarded
    :param fasta_headers: full-length headers that will replace those in assignments
    :param header_map:
    :param assignments: A dictionary of reference (lineage) and query names
    :param molecule_type: prot, nucl, or rrna? Parsed from command-line arguments
    :return: assignments with a full lineage for both reference (keys) and queries (values)
    """
    genome_position_re = re.compile("^([A-Z0-9._]+):.[0-9]+-[0-9]+[_]+.*")
    full_assignments = dict()
    for marker in assignments.keys():
        full_assignments[marker] = dict()
        entrez_query_list = dict()
        num_queries = 0
        for robust_classification in assignments[marker].keys():
            full_assignments[marker][robust_classification] = list()
            entrez_query_list[robust_classification] = list()
            for query in assignments[marker][robust_classification]:
                try:
                    original_header = header_map['>' + query]
                except KeyError:
                    sys.stderr.write("ERROR: unable to find " + query +
                                     " in header_map (constructed from either the input FASTA or .uc file).\n")
                    sys.stderr.write("This is probably an error stemming from `reformat_string()`.\n")
                    sys.exit(5)
                # print(query, original_header)
                database = molecule_type
                query_fields = query.split('_')
                q_accession = query_fields[0]
                if not re.search("Bacteria\|Archaea", query):
                    # Check for genomic coordinates, indicating these genes will be in the 'nucleotide' database
                    if genome_position_re.match(query):
                        q_accession = genome_position_re.match(query).group(1)
                        database = "dna"
                    else:
                        header_format_re, header_db, header_molecule = get_header_format(original_header, marker)
                        sequence_info = header_format_re.match(original_header)
                        q_accession, _, _, _ = return_sequence_info_groups(sequence_info, header_db, original_header)
                    # print(q_accession)
                    entrez_query_list[robust_classification].append((q_accession, database))
                    num_queries += 1
                else:
                    n_match = 0
                    for header in fasta_headers:
                        f_accession = header[1:].split('_')[0]
                        # Useful for headers containing the full lineage
                        if q_accession == f_accession:
                            full_assignments[marker][robust_classification].append('_'.join(header.split('_')[1:]))
                            n_match += 1
                        else:
                            pass
                    if n_match == 0:
                        sys.stderr.write("ERROR: Unable to find matching header for " + query + " in fasta!\n")
                        sys.exit(2)
                    elif n_match > 1:
                        sys.stderr.write("ERROR: headers with identical accessions were identified in fasta!\n")
                        sys.stderr.write("Offending accession: " + q_accession + "\n")
                    sys.exit(3)

        if num_queries > 1:
            sys.stdout.write("Downloading lineage information for " + marker + ":\n")
            step_proportion = setup_progress_bar(num_queries)
            acc = 0.0
            for ref in entrez_query_list:
                for entrez_query in entrez_query_list[ref]:
                    q_accession, database = entrez_query
                    lineage = get_lineage(q_accession, database)
                    if lineage:
                        full_assignments[marker][ref].append(lineage)

                    # Update the progress bar
                    acc += 1.0
                    if acc >= step_proportion:
                        acc -= step_proportion
                        sys.stdout.write('-')
                        sys.stdout.flush()
            if acc > 1:
                sys.stdout.write('-')
            sys.stdout.write("]\n")
        else:
            sys.stderr.write("ERROR: No accessions could be parsed from the FASTA headers.\n")
            sys.exit(11)

    return full_assignments


def filter_queries_by_taxonomy(assignments):
    """
    Removes queries with duplicate taxa - to prevent the taxonomic composition of the input
    from biasing the accuracy to over- or under-perform by classifying many sequences from very few groups.
    Also removes taxonomies with "*nclassified" in their lineage
    :param assignments:
    :return:
    """
    deduplicated_assignments = dict()
    unclassifieds = 0
    for marker in assignments:
        deduplicated_assignments[marker] = dict()
        for ref in assignments[marker]:
            deduplicated_assignments[marker][ref] = set()
            for query_taxonomy in assignments[marker][ref]:
                if re.search("nclassified", query_taxonomy):
                    unclassifieds += 1
                else:
                    deduplicated_assignments[marker][ref].add(query_taxonomy)
    if unclassifieds > 0:
        sys.stdout.write("\n\t" + str(unclassifieds) + " query sequences with unclassified taxonomies were removed.\n")
        sys.stdout.write("This is not a problem, its just they have unclassified somewhere in their lineages\n"
                         "(e.g. Unclassified Bacteria) and this is not good for assessing placement accuracy.\n\n")

    return deduplicated_assignments


def write_performance_table(args, clade_exclusion_strings, sensitivity):
    output = args.output + os.sep + "clade_exclusion_performance.tsv"
    try:
        output_handler = open(output, 'w')
    except IOError:
        sys.stderr.write("ERROR: Unable to open " + output + " for writing!")
        raise IOError

    output_handler.write("# Input file for testing: " + args.fasta_input + "\n")
    output_name = os.path.dirname(args.output)
    for line in clade_exclusion_strings:
        line += sensitivity + "\n"
        line = output_name + "\t" + line
        output_handler.write(line)

    output_handler.close()
    return


def main():
    args = get_arguments_()
    args.min_seq_length = 1
    # Load all the assignments from TreeSAPP output
    tmp_file = args.output + os.sep + "tmp_clade_exclusion_assignments.tsv"
    if os.path.exists(tmp_file):
        pre_assignments, n_saved = read_intermediate_assignments(args, tmp_file)
    else:
        n_saved = 0
        pre_assignments = dict()
    assignments, n_classified = read_table(args.assignment_table)
    fasta_headers = format_read_fasta(args.fasta_input, args.molecule, args, 10000).keys()
    # Used for sensitivity
    total_sequences = len(fasta_headers)
    sensitivity = str(round(float(n_classified / total_sequences), 3))
    sys.stdout.write("Sensitivity = " + sensitivity + "\n")
    if n_classified == n_saved:
        sys.stdout.write("Using taxonomic lineage information previously downloaded.\n")
        assignments = pre_assignments
    else:
        original_headers = get_headers(args.fasta_input)
        header_map = map_good_headers_to_ugly(original_headers)
        # Retrieve taxonomic lineage for each sequence from Entrez API or parsing the header
        assignments = map_full_headers(fasta_headers, header_map, assignments, args.molecule)
    write_intermediate_assignments(args, assignments)
    # This function could be replaced by a set in map_full_headers but I'd rather perform this task explicitly,
    # to make it more obvious this is being performed rather than hiding it with sets. Also easier to turn off :)
    assignments = filter_queries_by_taxonomy(assignments)
    # Get rid of some names, replace underscores with semi-colons
    assignments = clean_classification_names(assignments)
    # Load the reference lineages into a trie (prefix trie)
    clade_exclusion_strings = list()
    for table in args.tax_ids_table:
        taxonomic_tree, marker = all_possible_assignments(args, table)
        if marker not in assignments.keys():
            sys.stderr.write("WARNING: You provided a tax_ids file for " + marker +
                             " but no sequences were classified for this marker. "
                             "Everything else is okay, don't worry.\n")
            sys.stderr.flush()
        else:
            # Identify at which rank each sequence's clade was excluded: [K,P,C,O,F,G,S]
            rank_assigned_dict = identify_excluded_clade(args, assignments, taxonomic_tree, marker)
            # Determine the specificity for each rank
            clade_exclusion_strings = determine_specificity(rank_assigned_dict, marker, clade_exclusion_strings)
            determine_containment(args, marker, rank_assigned_dict)
    write_performance_table(args, clade_exclusion_strings, sensitivity)
    # Remove the intermediate file since this run completed successfully
    # if os.path.exists(tmp_file):
    #     os.remove(tmp_file)


if __name__ == "__main__":
    main()
