#!/usr/bin/env python3

__author__ = 'Connor Morgan-Lang'

import argparse
import sys
import pygtrie
import re
import os
import inspect
import shutil
from glob import glob
from random import randint

cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile(inspect.currentframe()))[0]))
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)
sys.path.insert(0, cmd_folder + os.sep + ".." + os.sep)
from fasta import format_read_fasta, write_new_fasta, get_headers
from create_treesapp_ref_data import get_header_format, register_headers
from utilities import clean_lineage_string, return_sequence_info_groups
from external_command_interface import setup_progress_bar, launch_write_command
from classy import ReferenceSequence
from file_parsers import parse_ref_build_params
from entrez_utils import get_lineage, read_accession_taxa_map, write_accession_lineage_map

rank_depth_map = {0: "Cellular organisms", 1: "Kingdom",
                  2: "Phylum", 3: "Class", 4: "Order",
                  5: "Family", 6: "Genus", 7: "Species",
                  8: "Strain"}


def get_arguments_():
    parser = argparse.ArgumentParser(add_help=False)
    required_args = parser.add_argument_group("Required arguments")
    required_args.add_argument('-o', '--output',
                               required=True,
                               help='TreeSAPP output directory from a previous analysis or clade exclusion analysis')
    required_args.add_argument('-i', '--fasta_input',
                               help='Your sequence input file (for TreeSAPP) in FASTA format',
                               required=True)
    required_args.add_argument("-r", "--reference_markers",
                               help="Short-form name of the marker gene to be tested (e.g. mcrA, pmoA, nosZ)",
                               required=True,
                               nargs='+')

    optopt = parser.add_argument_group("Optional options")
    optopt.add_argument('-m', '--molecule',
                        help='the type of input sequences (prot = Protein [DEFAULT]; dna = Nucleotide; rrna = rRNA)',
                        default='prot',
                        choices=['prot', 'dna', 'rrna'])
    optopt.add_argument("-l", "--length",
                        required=False, type=int, default=0,
                        help="Arbitrarily slice the input sequences to this length. "
                             "Useful for testing classification accuracy for fragments.")

    miscellaneous_opts = parser.add_argument_group("Miscellaneous options")
    miscellaneous_opts.add_argument("--graftm", default=None, required=False,
                                    help="The path to a GraftM package. Classifications performed using GraftM.")
    miscellaneous_opts.add_argument('--overwrite', action='store_true', default=False,
                                    help='overwrites previously processed output folders')
    miscellaneous_opts.add_argument('-v', '--verbose', action='store_true', default=False,
                                    help='Prints a more verbose runtime log')
    miscellaneous_opts.add_argument("-h", "--help",
                                    action="help",
                                    help="Show this help message and exit")

    args = parser.parse_args()

    if args.output[-1] != os.sep:
        args.output += os.sep
    args.treesapp = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + os.sep + ".." + os.sep

    if args.overwrite:
        if os.path.exists(args.output):
            shutil.rmtree(args.output)

    args.min_seq_length = 1

    if sys.version_info > (2, 9):
        args.py_version = 3
    else:
        args.py_version = 2

    return args


def read_graftm_classifications(assignment_file):
    assignments = dict()
    n_classified = 0
    assignments_handle = open(assignment_file, 'r')
    line = assignments_handle.readline()
    while line:
        fields = line.strip().split('\t')
        try:
            header, classified = fields
            if header and classified:
                n_classified += 1
                if classified not in assignments:
                    assignments[classified] = list()
                assignments[classified].append(header)
        except ValueError:
            sys.stderr.write("ERROR: Unable to parse line:")
            sys.stderr.write(str(line))
            sys.exit(1)
        line = assignments_handle.readline()

    assignments_handle.close()
    return assignments, n_classified


def read_marker_classification_table(assignment_file, marker=None):
    """
    Function for reading the tabular assignments file (currently marker_contig_map.tsv)
    Assumes column 2 is the TreeSAPP assignment and column 3 is the sequence header
    (leaving 1 for marker name and 4 for numerical abundance)
    :param assignment_file: Path to the file containing sequence phylogenetic origin and assignment
    :param marker: Optionally, the marker name can be explicitly provided. Necessary when analyzing GraftM outputs
    :return: dictionary whose keys are phylogenetic origin and values are lists of TreeSAPP assignments
    """
    assignments = dict()
    n_classified = 0
    if re.match(r".*_read_tax.tsv$", assignment_file):
        # This is a GraftM classification table
        marker_assignments, n_classified = read_graftm_classifications(assignment_file)
        assignments[marker] = marker_assignments
    else:
        assignments_handle = open(assignment_file, 'r')
        # This is the header line
        if not re.match("^Sample\tQuery\tMarker\tTaxonomy\tConfident_Taxonomy\tAbundance\tInternal_node\tLikelihood\tLWR\tWTD$",
                        assignments_handle.readline()):
            sys.stderr.write("ERROR: header of assignments file is unexpected!\n")
            raise AssertionError

        # First line in the table containing data
        line = assignments_handle.readline()
        while line:
            fields = line.strip().split('\t')
            try:
                _, header, marker, _, rob_class, _, _, _, _, _ = fields
                if marker and rob_class:
                    n_classified += 1
                    if marker not in assignments:
                        assignments[marker] = dict()
                    if rob_class not in assignments[marker]:
                        assignments[marker][rob_class] = list()
                    assignments[marker][rob_class].append(header)
            except ValueError:
                sys.stderr.write("ERROR: Unable to parse line:\n")
                sys.stderr.write(str(line))
                sys.exit(1)
            line = assignments_handle.readline()

        assignments_handle.close()
    return assignments, n_classified


def read_intermediate_assignments(args, inter_class_file):
    if args.verbose:
        sys.stderr.write("Reading " + inter_class_file + " for saved intermediate data... ")
    assignments = dict()
    n_saved = 0
    try:
        file_handler = open(inter_class_file, 'r')
    except OSError:
        sys.stderr.write("ERROR: Unable to open " + inter_class_file + " for reading!\n")
        raise OSError
    line = file_handler.readline()
    while line:
        line = line.strip()
        try:
            marker, ref, queries = line.split('\t')
        except ValueError:
            sys.stderr.write("\tWARNING: " + inter_class_file + " is incorrectly formatted so regenerating instead.\n")
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


def write_intermediate_assignments(args, inter_class_file, assignments):
    if args.verbose:
        sys.stderr.write("Saving intermediate data... ")
    try:
        file_handler = open(inter_class_file, 'w')
    except OSError:
        sys.stderr.write("ERROR: Unable to open " + inter_class_file + " for writing!\n")
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


def grab_graftm_taxa(tax_ids_file):
    taxonomic_tree = pygtrie.StringTrie(separator='; ')
    with open(tax_ids_file) as tax_ids:
        header = tax_ids.readline().strip()
        if header != "tax_id,parent_id,rank,tax_name,root,kingdom,phylum,class,order,family,genus,species":
            raise AssertionError("ERROR: Unable to handle format of " + tax_ids_file + "!")
        line = tax_ids.readline().strip()
        while line:
            try:
                _, _, _, _, _, k_, p_, c_, o_, f_, g_, s_ = line.split(',')
            except IndexError:
                raise AssertionError("ERROR: Unexpected format of line in " + tax_ids_file + ":\n" + line)
            ranks = [k_, p_, c_, o_, f_, g_, s_]
            lineage_list = []
            for rank in ranks:
                if rank:
                    lineage_list.append(rank)
            lineage = clean_lineage_string('; '.join(lineage_list))
            i = 0
            ranks = len(lineage)
            while i < len(lineage):
                taxonomic_tree["; ".join(lineage.split("; ")[:ranks - i])] = True
                i += 1

            line = tax_ids.readline().strip()

    return taxonomic_tree


def all_possible_assignments(args, tax_ids_file):
    taxonomic_tree = pygtrie.StringTrie(separator='; ')
    # if os.path.exists(tax_ids_file):
    #     file_name = os.path.basename(tax_ids_file)
    #     if re.match("^tax_ids_(.*).txt", file_name):
    #         marker = re.match("^tax_ids_(.*).txt", file_name).group(1)
    #     else:
    #         sys.stderr.write("ERROR: Format of tax_ids file (" + tax_ids_file +
    #                          ") is unexpected. Unable to parse marker name! Exiting...\n")
    #         sys.exit(7)
    # else:
    #     raise IOError("File doesn't exist: " + tax_ids_file + "\n")
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
    return taxonomic_tree


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
        sys.stdout.write("Number of unique taxonomies that sequences were assigned to = " +
                         str(len(assignment_dict[marker].keys())) + "\n")
    for ref_lineage in assignment_dict[marker]:
        if args.verbose:
            sys.stdout.write("Reference lineage: " + ref_lineage + "\n")
        for query_lineage in assignment_dict[marker][ref_lineage]:
            # if query_lineage == ref_lineage:
            #     if args.verbose:
            #         sys.stdout.write("\tQuery lineage: " + query_lineage + ", ")
            #         sys.stdout.write("Optimal lineage: " + ref_lineage + "\n")
            #         sys.stdout.flush()
            # While the query_lineage contains clades which are not in the reference trie,
            # remove the taxonomic rank and traverse again. Complexity: O(ln(n))
            contained_taxonomy = query_lineage
            while not trie.__contains__(contained_taxonomy) and len(contained_taxonomy.split('; ')) > 1:
                contained_taxonomy = "; ".join(contained_taxonomy.split('; ')[:-1])
            if len(contained_taxonomy.split("; ")) <= 7:
                rank_excluded = rank_depth_map[len(contained_taxonomy.split("; ")) + 1]
                if args.verbose and contained_taxonomy != ref_lineage:
                    sys.stdout.write("\tRank excluded: " + rank_excluded + "\n")
                    sys.stdout.write("\t\tQuery lineage:   " + query_lineage + "\n")
                    sys.stdout.write("\t\tOptimal lineage: " + contained_taxonomy + "\n")
                rank_assigned_dict[rank_excluded].append({ref_lineage: (contained_taxonomy, query_lineage)})
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


def download_taxonomic_lineage_for_queries(entrez_query_list):

    sys.stdout.write("Downloading lineage information for each sequence accession from Entrez\n")

    num_queries = len(entrez_query_list)
    download_accumulator = 0
    full_assignments = dict()
    if num_queries > 1:
        step_proportion = setup_progress_bar(num_queries)
        acc = 0.0
        # for ref in entrez_query_list:
        #     if ref not in full_assignments:
        #         full_assignments[ref] = list()
        #     for entrez_query in entrez_query_list[ref]:
        for entrez_query in entrez_query_list:
            q_accession, database = entrez_query
            lineage = get_lineage(q_accession, database)
            if lineage:
                full_assignments[q_accession] = lineage
                download_accumulator += 1

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

    sys.stdout.write("\t" + str(download_accumulator) + " taxonomic lineages downloaded\n")

    return full_assignments


def finalize_ref_seq_lineages(test_ref_sequences, accession_lineages):
    for ref_seq in test_ref_sequences:
        if not ref_seq.lineage:
            try:
                ref_seq.lineage = accession_lineages[ref_seq.accession]
            except KeyError:
                sys.stderr.write("ERROR: Lineage information was not retrieved for " + ref_seq.accession + "!\n")
        else:
            pass
    return test_ref_sequences


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
    marker_assignments = dict()
    entrez_query_list = list()
    num_queries = 0
    for robust_classification in assignments.keys():
        marker_assignments[robust_classification] = list()
        for query in assignments[robust_classification]:
            try:
                original_header = header_map['>' + query]
            except KeyError:
                sys.stderr.write("ERROR: unable to find " + query +
                                 " in header_map (constructed from either the input FASTA or .uc file).\n")
                sys.stderr.write("This is probably an error stemming from `reformat_string()`.\n")
                sys.exit(5)
            # print(query, original_header)
            database = molecule_type
            q_accession = ""
            if not re.search("Bacteria\|Archaea", query):
                # Check for genomic coordinates, indicating these genes will be in the 'nucleotide' database
                if genome_position_re.match(query):
                    q_accession = genome_position_re.match(query).group(1)
                    database = "dna"
                else:
                    header_format_re, header_db, header_molecule = get_header_format(original_header)
                    sequence_info = header_format_re.match(original_header)
                    q_accession, _, _, _, _ = return_sequence_info_groups(sequence_info, header_db, original_header)
                # print(q_accession)
                entrez_query_list.append((q_accession, database))
                marker_assignments[robust_classification].append(q_accession)
                num_queries += 1
            else:
                n_match = 0
                for header in fasta_headers:
                    f_accession = header[1:].split('_')[0]
                    # Useful for headers containing the full lineage
                    if q_accession == f_accession:
                        marker_assignments[robust_classification].append('_'.join(header.split('_')[1:]))
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

    return marker_assignments, entrez_query_list


def assign_lineages(complete_ref_seqs, assignments):
    """

    :param assignments:
    :param complete_ref_seqs:
    :return: lineages_list = full_assignments[marker][ref_classification]
    """
    for marker in assignments:
        for robust_classification in assignments[marker]:
            tmp_lineage_list = []
            for accession in assignments[marker][robust_classification]:
                for ref_seq in complete_ref_seqs:
                    if ref_seq.accession == accession:
                        tmp_lineage_list.append(ref_seq.lineage)
            assignments[marker][robust_classification] = tmp_lineage_list
    return assignments


def map_headers_to_lineage(assignments, ref_sequences):
    """
    The alternative function to map_full_headers. Using ReferenceSequence objects,
    :param assignments:
    :param ref_sequences:
    :return:
    """
    lineage_assignments = dict()
    for marker in assignments:
        lineage_assignments[marker] = dict()
        for assigned_lineage in assignments[marker].keys():
            classified_headers = assignments[marker][assigned_lineage]
            lineage_assignments[marker][assigned_lineage] = list()
            for query in classified_headers:
                for ref_seq in ref_sequences:
                    if ref_seq.accession == query:
                        lineage_assignments[marker][assigned_lineage].append(ref_seq.lineage)
    return lineage_assignments


def get_unclassified_rank(pos, split_lineage):
    """
    Recursive function to retrieve the first rank at which the
    :param pos:
    :param split_lineage:
    :return:
    """
    if re.search("nclassified", split_lineage[pos]):
        return pos
    else:
        pos = get_unclassified_rank(pos+1, split_lineage)
    return pos


def pick_taxonomic_representatives(ref_seqs_list, taxonomic_filter_stats, max_cluster_size=5):
    """
    Removes queries with duplicate taxa - to prevent the taxonomic composition of the input
    from biasing the accuracy to over- or under-perform by classifying many sequences from very few groups.
    Also removes taxonomies with "*nclassified" in their lineage
    :param ref_seqs_list: A dictionary mapping accessions to lineages that need to be filtered
    :param taxonomic_filter_stats: A dictionary for tracking the number sequences filtered, those unique, etc.
    :param max_cluster_size: The maximum number of sequences representing a taxonomic cluster
    :return: dereplicated_lineages dict with lineages mapping to a single accession
    """
    good_classified_lineages = dict()
    dereplicated_lineages = dict()
    num_rep_seqs = 0
    for ref_seq in ref_seqs_list:
        query_taxonomy = ref_seq.lineage
        if query_taxonomy not in good_classified_lineages:
            good_classified_lineages[query_taxonomy] = list()
        if re.search("unclassified|environmental sample", query_taxonomy, re.IGNORECASE):
            # # Remove taxonomic lineages that are unclassified at the Phylum level or higher
            # unclassified_depth = get_unclassified_rank(0, query_taxonomy.split("; "))
            # if unclassified_depth > 4:
            #     good_classified_lineages[query_taxonomy].append(ref_seq.accession)
            # else:
            taxonomic_filter_stats["Unclassified"] += 1
        else:
            taxonomic_filter_stats["Classified"] += 1
            good_classified_lineages[query_taxonomy].append(ref_seq.accession)

    # In order to maintain consistency among multiple runs with the same input
    # a separate loop is required for sorting
    taxonomic_filter_stats["Max"] = 0
    taxonomic_filter_stats["Min"] = 1000
    taxonomic_filter_stats["Mean"] = 0
    for query_taxonomy in good_classified_lineages:
        dereplicated_lineages[query_taxonomy] = list()
        i = 0
        lineage_candidates = sorted(good_classified_lineages[query_taxonomy])
        while i < len(lineage_candidates) and i < max_cluster_size:
            dereplicated_lineages[query_taxonomy].append(lineage_candidates[i])
            i += 1

        # Generate stats on the taxonomic clusters
        cluster_size = len(dereplicated_lineages[query_taxonomy])
        num_rep_seqs += cluster_size
        if cluster_size == 0:
            # This is a taxonomic lineage removed by the "Unclassified" filter
            continue
        if cluster_size > taxonomic_filter_stats["Max"]:
            taxonomic_filter_stats["Max"] = cluster_size
        if cluster_size < taxonomic_filter_stats["Min"]:
            taxonomic_filter_stats["Min"] = cluster_size
        taxonomic_filter_stats["Mean"] = round(float(num_rep_seqs / len(dereplicated_lineages)), 2)

    taxonomic_filter_stats["Unique_taxa"] += len(dereplicated_lineages)

    sys.stdout.write("\t" + str(num_rep_seqs) + " representative sequences will be used for TreeSAPP analysis.\n")
    sys.stdout.flush()

    sys.stdout.write("Representative sequence stats:\n\t")
    sys.stdout.write("Maximum representative sequences for a taxon " + str(taxonomic_filter_stats["Max"]) + "\n\t")
    sys.stdout.write("Minimum representative sequences for a taxon " + str(taxonomic_filter_stats["Min"]) + "\n\t")
    sys.stdout.write("Mean representative sequences for a taxon " + str(taxonomic_filter_stats["Mean"]) + "\n")

    return dereplicated_lineages, taxonomic_filter_stats


def select_rep_seqs(deduplicated_assignments, test_sequences):
    """
    Function for creating a fasta-formatted dict from the accessions representing unique taxa in the test sequences
    :param deduplicated_assignments:
    :param test_sequences:
    :return:
    """
    deduplicated_fasta_dict = dict()
    for lineage in sorted(deduplicated_assignments):
        for accession in deduplicated_assignments[lineage]:
            matched = False
            for ref_seq in test_sequences:
                if ref_seq.accession == accession:
                    deduplicated_fasta_dict[accession] = ref_seq.sequence
                    matched = True
            if not matched:
                sys.stderr.write("ERROR: Unable to find accession (" + accession + ") in accession-lineage map\n")
    return deduplicated_fasta_dict


def filter_queries_by_taxonomy(taxonomic_lineages):
    """
    Removes queries with duplicate taxa - to prevent the taxonomic composition of the input
    from biasing the accuracy to over- or under-perform by classifying many sequences from very few groups.
    Also removes taxonomies with "*nclassified" in their lineage
    :param taxonomic_lineages: A list of lineages that need to be filtered
    :return: A list that contains no more than 3 of each query taxonomy (arbitrarily normalized) and counts
    """
    unclassifieds = 0
    classified = 0
    unique_query_taxonomies = 0
    lineage_enumerator = dict()
    normalized_lineages = list()
    for query_taxonomy in sorted(taxonomic_lineages):
        can_classify = True
        if re.search("unclassified", query_taxonomy, re.IGNORECASE):
            # Remove taxonomic lineages that are unclassified at the Phylum level or higher
            unclassified_depth = get_unclassified_rank(0, query_taxonomy.split("; "))
            if unclassified_depth <= 3:
                can_classify = False

        if can_classify:
            if query_taxonomy not in lineage_enumerator:
                lineage_enumerator[query_taxonomy] = 0
                classified += 1
            if lineage_enumerator[query_taxonomy] < 3:
                normalized_lineages.append(query_taxonomy)
            lineage_enumerator[query_taxonomy] += 1
        else:
            unclassifieds += 1
    unique_query_taxonomies += len(lineage_enumerator)

    return normalized_lineages, unclassifieds, classified, unique_query_taxonomies


def write_performance_table(args, performance_table, clade_exclusion_strings, sensitivity):
    try:
        output_handler = open(performance_table, 'w')
    except IOError:
        sys.stderr.write("ERROR: Unable to open " + performance_table + " for writing!")
        raise IOError

    output_handler.write("# Input file for testing: " + args.fasta_input + "\n")
    output_name = os.path.dirname(args.output)
    for line in clade_exclusion_strings:
        line += sensitivity + "\n"
        line = output_name + "\t" + line
        output_handler.write(line)

    output_handler.close()
    return


def correct_accession(description):
    header_format_re, header_db, header_molecule = get_header_format(description)
    sequence_info = header_format_re.match(description)
    accession, _, _, _, _ = return_sequence_info_groups(sequence_info, header_db, description)
    return accession


def build_entrez_queries(test_sequences, molecule_type):
    """
    Function to create data collections to fulfill entrez query searches
    :param test_sequences: A list of ReferenceSequence objects - lineage information to be filled
    :param molecule_type:
    :return: dictionary with
    """
    genome_position_re = re.compile("^([A-Z0-9._]+):.[0-9]+-[0-9]+[ ]+.*")
    entrez_query_list = list()
    for ref_seq in test_sequences:
        database = molecule_type
        if not ref_seq.lineage:
            # Check for genomic coordinates, indicating these genes will be in the 'nucleotide' database
            # TODO: test this with the unformatted headers
            if genome_position_re.match(ref_seq.description):
                q_accession = genome_position_re.match(ref_seq.description).group(1)
                database = "dna"
            else:
                q_accession = ref_seq.accession
            entrez_query_list.append((q_accession, database))

    return entrez_query_list


def load_ref_seqs(fasta_dict, header_registry=None):
    """
    Function for loading a fasta-formatted dictionary into ReferenceSequence objects
    :param fasta_dict:
    :param header_registry: An optional dictionary of Header objects
    :return:
    """
    ref_seq_list = list()
    for header in fasta_dict.keys():
        accession = header[1:].split('_')
        ref_seq = ReferenceSequence()
        ref_seq.accession = accession  # Could be incorrect, still needs to be edited in build_entrez_queries
        ref_seq.sequence = fasta_dict[header]
        if header_registry:
            # format of header_map: header_map[reformat_string(original_header)] = original_header
            # This is used to recover the real accession ID later on
            for num_id in header_registry:
                if header_registry[num_id].formatted == header:
                    ref_seq.description = header_registry[num_id].original
            if ref_seq.description == "":
                sys.stderr.write("ERROR: unable to find " + header + " in header_map!\n")
                sys.exit(13)

            if re.search("Bacteria\|Archaea", header):
                # TODO: test this. Assuming taxonomic lineage immediately follows sequence accession
                ref_seq.lineage = ' '.join(header.split(' ')[1:])

        # Lineage, organism, and description data are still potentially missing at this point
        ref_seq_list.append(ref_seq)
    return ref_seq_list


def main():
    """
    Different ways to run this script:
        1. Provide it a FASTA file for which it will determine the taxonomic lineage for each sequence
         and run all taxonomic representative sequences with TreeSAPP then analyze via clade exclusion
        2. Run TreeSAPP on a set of sequences then provide Clade_exclusion_analyzer.py with the output directory
         and the input FASTA file analyzed. No taxonomic filtering will be performed.
    :return:
    """
    args = get_arguments_()
    if args.verbose:
        sys.stdout.write("\nBeginning clade exclusion analysis\n")
        sys.stdout.flush()

    ##
    # Define locations of files TreeSAPP outputs
    ##
    inter_class_file = args.output + os.sep + "tmp_clade_exclusion_assignments.tsv"
    accession_map_file = args.output + os.sep + "accession_id_lineage_map.tsv"
    test_rep_taxa_fasta = args.output + os.sep + "representative_taxa_sequences.fasta"
    performance_table = args.output + os.sep + "clade_exclusion_performance.tsv"
    # Determine the analysis stage and user's intentions with four booleans
    # Working with a pre-existing TreeSAPP output directory
    if os.path.exists(args.output):
        extant = True
        sys.stdout.write("The output directory has been detected.\n")
    else:
        extant = False
    taxa_filter = False  # Whether the test sequences were screened for redundant taxa
    accessions_downloaded = False  # Whether the accessions have been downloaded from Entrez
    classified = False  # Has TreeSAPP been completed for these sequences?
    treesapp_output_dir = args.output + "TreeSAPP_output" + os.sep
    if extant:
        if os.path.isfile(accession_map_file):
            accessions_downloaded = True
        if os.path.isfile(test_rep_taxa_fasta):
            taxa_filter = True
            args.output = treesapp_output_dir
    else:
        os.makedirs(args.output)
        args.output = treesapp_output_dir

    if not os.path.isdir(treesapp_output_dir):
        os.makedirs(treesapp_output_dir)

    if args.graftm:
        graftm_files = glob(args.output + os.sep + "*" + os.sep + "*_read_tax.tsv")
        if len(graftm_files) == 1:
            classification_table = glob(args.output + os.sep + "*" + os.sep + "*_read_tax.tsv")[0]
        else:
            classification_table = ''
    else:
        classification_table = args.output + os.sep + "final_outputs" + os.sep + "marker_contig_map.tsv"
    if os.path.isfile(classification_table):
        classified = True

    # Load FASTA data
    fasta_dict = format_read_fasta(args.fasta_input, args.molecule, args, 110)
    if args.length:
        for seq_id in fasta_dict:
            if len(fasta_dict[seq_id]) < args.length:
                sys.stderr.write("WARNING: " + seq_id + " sequence is shorter than " + str(args.length))
            else:
                max_stop = len(fasta_dict[seq_id]) - args.length
                random_start = randint(0, max_stop)
                fasta_dict[seq_id] = fasta_dict[seq_id][random_start:random_start+args.length]
    header_registry = register_headers(args, get_headers(args.fasta_input))
    # Load the query test sequences as ReferenceSequence objects
    test_ref_sequences = load_ref_seqs(fasta_dict, header_registry)
    complete_ref_sequences = list()
    for ref_seq in test_ref_sequences:
        ref_seq.accession = correct_accession(ref_seq.description)
        complete_ref_sequences.append(ref_seq)
    test_ref_sequences.clear()

    if args.verbose:
        sys.stdout.write("\tNumber of input sequences =\t" + str(len(complete_ref_sequences)) + "\n")
        sys.stdout.flush()

    n_classified = 0
    marker_assignments = {}
    # Set up a dictionary for tracking the taxonomic filtering
    taxonomic_filter_stats = dict()
    taxonomic_filter_stats["Unclassified"] = 0
    taxonomic_filter_stats["Classified"] = 0
    taxonomic_filter_stats["Unique_taxa"] = 0

    if args.graftm:
        gpkg_re_match = re.match("[0-9.]+.*\.(\w+)\.gpkg", os.path.basename(args.graftm))
        marker = gpkg_re_match.group(1)
        # try:
        #     marker = gpkg_re_match.group(1)
        # except:
        #     raise Exception("ERROR: Unable to parse marker from GraftM package name: " + args.graftm + "\n")
    else:
        marker = None

    # TODO: fix bug where marker_contig_map.tsv is missing when accession_lineage_map present
    # Checkpoint one: does anything exist?
    # If not, begin by downloading lineage information for each sequence accession
    if not extant and not accessions_downloaded and not taxa_filter and not classified:
        # User hasn't analyzed anything, sequences will be taxonomically screened then analyzed
        # NOTE: This is only supported for analyzing a single marker gene
        extant = True
        entrez_query_list = build_entrez_queries(complete_ref_sequences, args.molecule)
        if len(entrez_query_list) > 0:
            accession_lineage_map = download_taxonomic_lineage_for_queries(entrez_query_list)
        else:
            sys.stdout.write("No sequences with lineage information that needs to be downloaded from Entrez.\n")
            accession_lineage_map = dict()
        complete_ref_sequences = finalize_ref_seq_lineages(complete_ref_sequences, accession_lineage_map)
        write_accession_lineage_map(accession_map_file, accession_lineage_map)
        accessions_downloaded = True
    elif extant and accessions_downloaded:
        # File being read should contain accessions mapped to their lineages for all sequences in input FASTA
        read_accession_lineage_dict = read_accession_taxa_map(accession_map_file)
        complete_ref_sequences = finalize_ref_seq_lineages(complete_ref_sequences, read_accession_lineage_dict)

    # Checkpoint two: Do we have accessions? Are the sequences filtered by taxonomy or are the sequences classified?
    if extant and accessions_downloaded and not taxa_filter and not classified:
        sys.stdout.write("Selecting representative sequences for each taxon from downloaded lineage information.\n")

        # Filter the sequences from redundant taxonomic lineages, picking up to 5 representative sequences
        representative_seqs, taxonomic_filter_stats = pick_taxonomic_representatives(complete_ref_sequences,
                                                                                     taxonomic_filter_stats)
        deduplicated_fasta_dict = select_rep_seqs(representative_seqs, complete_ref_sequences)
        total_sequences = len(deduplicated_fasta_dict.keys())
        write_new_fasta(deduplicated_fasta_dict, test_rep_taxa_fasta)
        taxa_filter = True
    elif extant and accessions_downloaded and taxa_filter:
        deduplicated_fasta_dict = format_read_fasta(test_rep_taxa_fasta, args.molecule, args)
        total_sequences = len(deduplicated_fasta_dict.keys())
    else:

        total_sequences = len(fasta_dict)

    # Checkpoint three: We have accessions linked to taxa, and sequences to analyze with TreeSAPP, but not classified
    if extant and accessions_downloaded and taxa_filter and not classified:
        sys.stdout.write("Analyzing the " + str(total_sequences) + " representative sequences with TreeSAPP.\n")
        # Run TreeSAPP against the provided tax_ids file and the unique taxa FASTA file
        if args.length:
            min_seq_length = str(min(args.length - 10, 50))
        else:
            min_seq_length = str(50)

        if args.graftm:
            classify_command = ["graftM", "graft"]
            classify_command += ["--forward", test_rep_taxa_fasta]
            classify_command += ["--graftm_package", args.graftm]
            classify_command += ["--threads", str(4)]
            classify_command += ["--assignment_method", "pplacer"]
            # classify_command += ["--assignment_method", "kraken"]
            classify_command += ["--search_method", "hmmsearch"]
            # classify_command += ["--search_method", "kraken"]
            classify_command += ["--output_directory", treesapp_output_dir]
            classify_command += ["--input_sequence_type", "aminoacid"]
            classify_command.append("--force")

            graftm_prefix = '.'.join(os.path.basename(test_rep_taxa_fasta).split('.')[:-1])
            classification_table = os.sep.join([args.output, graftm_prefix, graftm_prefix + "_read_tax.tsv"])

        else:
            classify_command = ["./treesapp.py", "-i", test_rep_taxa_fasta,
                                "-o", treesapp_output_dir,
                                "-m", args.molecule,
                                "-T", str(4),
                                "--filter_align",
                                "--min_likelihood", str(0.1),
                                "--placement_parser", "lca",
                                "--min_seq_length", min_seq_length,
                                "--verbose",
                                "--overwrite",
                                "--delete"]
        sys.stdout.write("Command used:\n" + ' '.join(classify_command) + "\n")
        launch_write_command(classify_command, False)
        classified = True

    # Checkpoint four: everything has been prepared, only need to parse the classifications and map lineages
    if extant and accessions_downloaded and taxa_filter and classified:
        sys.stdout.write("Finishing up the mapping of classified, filtered taxonomic sequences.\n")
        if os.path.isfile(classification_table):
            assignments, n_classified = read_marker_classification_table(classification_table, marker)
        else:
            sys.stderr.write("\nERROR: marker_contig_map.tsv is missing from output directory '" +
                             os.path.basename(classification_table) + "'\n")
            sys.stderr.write("Please remove this directory and re-run.\n")
            sys.exit(3)
        marker_assignments = map_headers_to_lineage(assignments, complete_ref_sequences)

    # In the case of a prior TreeSAPP analysis without taxonomic sequence filtering (external of Clade_exclusion)
    # User is just interested in seeing the clade exclusion analysis results again
    if extant and not taxa_filter and classified:
        sys.stdout.write("Outputs from a previous unsupervised TreeSAPP analysis found.\n")
        # Read the classification table
        assignments, n_classified = read_marker_classification_table(classification_table)
        redownload = True
        sys.exit("Untested.\n")

        # if os.path.isfile(inter_class_file):
        #     pre_assignments, n_saved = read_intermediate_assignments(args, inter_class_file)

        if accessions_downloaded and len(complete_ref_sequences) == len(fasta_dict):
            # Accessions have already been loaded into complete_ref_sequences above
            sys.stdout.write("Using accessions that have been previously downloaded.\n")
            for marker in assignments:
                for rob_class in assignments[marker]:
                    bad_headers = assignments[marker][rob_class]
                    assignments[marker][rob_class] = [correct_accession('>' + header) for header in bad_headers]
            assignments = assign_lineages(complete_ref_sequences, assignments)
            redownload = False

        if redownload:
            # Retrieve taxonomic lineage for each sequence from Entrez API or parsing the header
            # Only necessary if the outputs are from an externally-generated TreeSAPP analysis
            full_assignments = dict()
            for marker in assignments.keys():
                full_assignments[marker], entrez_query_list = map_full_headers(fasta_dict.keys(), header_map,
                                                                               assignments[marker], args.molecule)
                if len(entrez_query_list) > 0:
                    accession_lineage_map = download_taxonomic_lineage_for_queries(entrez_query_list)
                    complete_ref_sequences = finalize_ref_seq_lineages(complete_ref_sequences, accession_lineage_map)

            assignments = assign_lineages(complete_ref_sequences, full_assignments)
            write_accession_lineage_map(accession_map_file, complete_ref_sequences)

        # This function could be replaced by a set in map_full_headers but I'd rather perform this task explicitly,
        # to make it obvious this is being performed rather than hiding it with sets. Also easier to turn off :)
        marker_assignments = dict()
        for marker in assignments:
            marker_assignments[marker] = dict()
            for ref in assignments[marker]:
                lineages_list = assignments[marker][ref]
                marker_assignments[marker][ref], i, j, k = filter_queries_by_taxonomy(lineages_list)
                taxonomic_filter_stats["Unclassified"] += i
                taxonomic_filter_stats["Classified"] += j
                taxonomic_filter_stats["Unique_taxa"] += k

    ##
    # On to the standard clade-exclusion analysis...
    ##
    if taxonomic_filter_stats["Classified"] != taxonomic_filter_stats["Unique_taxa"]:
        sys.stdout.write(
            "\n\t" + str(taxonomic_filter_stats["Classified"] - taxonomic_filter_stats["Unique_taxa"]) +
            " duplicate query taxonomies removed.\n")

    if taxonomic_filter_stats["Unclassified"] > 0:
        sys.stdout.write(
            "\t" + str(taxonomic_filter_stats["Unclassified"]) +
            " query sequences with unclassified taxonomies were removed.\n")
        sys.stdout.write("This is not a problem, its just they have 'unclassified' somewhere in their lineages\n"
                         "(e.g. Unclassified Bacteria) and this is not good for assessing placement accuracy.\n\n")

    if n_classified == 0 and marker_assignments == {}:
        sys.stderr.write("\nERROR: Outputs from previous TreeSAPP or "
                         "Clade-exclusion analysis were incorrectly inferred to be present. Oopsie!\n")
        sys.stderr.write("The developer(s) should be made aware of this - "
                         "please post an issue to the GitHub page with your command used. Thanks!\n")
        sys.exit()

    sensitivity = str(round(float(n_classified / total_sequences), 3))
    sys.stdout.write("Sensitivity = " + sensitivity + "\n")

    # Write the intermediate classifications to a file
    write_intermediate_assignments(args, inter_class_file, marker_assignments)

    # Get rid of some names, replace underscores with semi-colons
    assignments = clean_classification_names(marker_assignments)

    # Load the reference lineages into a trie (prefix trie)
    clade_exclusion_strings = list()
    # tax_ids_tables = list()
    args.targets = ["ALL"]
    marker_build_dict = parse_ref_build_params(args)

    for name in args.reference_markers:
        # Alert the user if the denominator format was (incorrectly) provided to Clade_exclusion_analyzer
        if re.match("[A-Z][0-9]{4}", name):
            code_name = name
            marker = marker_build_dict[code_name].cog
        elif len(name) <= 6:
            # Assuming the provided name is a short gene name (e.g. mcrA, nifH)
            marker = name
            code_name = ""
            for denominator in marker_build_dict:
                if marker_build_dict[denominator].cog == marker:
                    code_name = denominator
                    break
            if not code_name:
                raise AssertionError("Unable to identify the gene name from the code name '" + name + "'.")
        else:
            sys.stderr.write("ERROR: Wrong format for the reference code_name provided: " + name + "\n")
            sys.exit(9)

        if args.graftm:
            tax_ids_file = glob(os.sep.join([args.graftm, "*refpkg", "*_taxonomy.csv"]))[0]
        else:
            tax_ids_file = os.sep.join([args.treesapp, "data", "tree_data", "tax_ids_" + marker + ".txt"])
        if os.path.exists(tax_ids_file):
            if args.graftm:
                taxonomic_tree = grab_graftm_taxa(tax_ids_file)
            else:
                taxonomic_tree = all_possible_assignments(args, tax_ids_file)
        else:
            sys.stderr.write("ERROR: Unable to find taxonomy IDs table: " + tax_ids_file + "\n")
            sys.exit(3)
        if code_name not in assignments.keys() and marker not in assignments.keys():
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
    write_performance_table(args, performance_table, clade_exclusion_strings, sensitivity)


if __name__ == "__main__":
    main()
