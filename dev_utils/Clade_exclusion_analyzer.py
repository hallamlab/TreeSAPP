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
from create_treesapp_ref_data import get_lineage

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
                               required=True)
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
    if not re.match("^Marker\tTaxonomy\tQuery\tAbundance$", assignments_handle.readline()):
        sys.stderr.write("ERROR: header of assignments file is unexpected!\n")
        raise AssertionError
    # First line in the table containing data
    line = assignments_handle.readline()
    while line:
        fields = line.strip().split('\t')
        try:
            marker, classified, header, abundance = fields
            if classified:
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


def read_intermediate_assignments(args, tmp_file):
    if args.verbose:
        sys.stderr.write("NOTIFICATION: Reading " + tmp_file + " for saved intermediate data... ")
    assignments = dict()
    try:
        file_handler = open(tmp_file, 'r')
    except OSError:
        sys.stderr.write("ERROR: Unable to open " + tmp_file + " for reading!\n")
        raise OSError
    line = file_handler.readline()
    while line:
        line = line.strip()
        ref, queries = line.split('\t')
        queries_list = queries.split(',')
        assignments[ref] = queries_list
        line = file_handler.readline()

    file_handler.close()
    if args.verbose:
        sys.stderr.write("done.\n")
        sys.stderr.flush()
    return assignments


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
    for ref, queries in assignments.items():
        cleaned_ref = clean_lineage_string(ref)
        if not cleaned_ref:
            if args.verbose:
                sys.stderr.write("WARNING: reference assignment is empty after cleaning the lineage:\n" + ref + "\n")
        else:
            assignments_string += cleaned_ref + "\t"
            assignments_string += ','.join(queries) + "\n"

    file_handler.write(assignments_string)
    file_handler.close()
    if args.verbose:
        sys.stderr.write("done.\n")
        sys.stderr.flush()
    return


def all_possible_assignments(args, tax_ids_file):
    taxonomic_tree = pygtrie.StringTrie(separator='; ')
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


def identify_excluded_clade(args, assignment_dict, trie):
    """
    Using the taxonomic information from the sequence headers and the lineages of the reference sequence,
    this function determines the rank at which each sequence's clade is excluded.
    These data are returned and sorted in the form of a dictionary.
    :param: assignment_dict:
    :param: trie: A pygtrie.StringTrie object containing all reference sequence lineages
    :return: rank_assigned_dict; key is rank, values are tuples of optimal assignment and actual assignment
    """
    rank_assigned_dict = dict()
    for depth in rank_depth_map:
        rank_assigned_dict[rank_depth_map[depth]] = list()
    if args.verbose:
        sys.stdout.write("Unique taxonomies sequences were assigned = " + str(len(assignment_dict.keys())) + "\n")
    for ref_lineage in assignment_dict:
        if args.verbose:
            sys.stdout.write("Reference lineage: " + ref_lineage + "\n")
        for query_lineage in assignment_dict[ref_lineage]:
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
            rank_assigned_dict[rank_depth_map[len(contained_taxonomy.split("; ")) + 1]].append(
                {ref_lineage: (contained_taxonomy, query_lineage)})
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


def determine_specificity(rank_assigned_dict):
    """
    Correct if: optimal_assignment == query_lineage
    Correct if:
    :param rank_assigned_dict:
    :return:
    """
    sys.stdout.write("Clade-level specificities:\n")
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
            sys.stdout.write("NA\n")
        else:
            for assignments in rank_assigned_dict[rank]:
                # if rank == "Family":
                #     print(assignments)
                for classified in assignments:
                    optimal, query = assignments[classified]
                    if optimal == classified:
                        offset = 0
                        correct += 1
                    else:
                        offset = determine_offset(classified, optimal)
                        incorrect += 1
                    taxonomic_distance[offset] += 1
            rank_total = correct + incorrect
            sys.stdout.write(str(rank_total) + "\t")
            for dist in taxonomic_distance:
                if taxonomic_distance[dist] > 0:
                    taxonomic_distance[dist] = str(round(float(taxonomic_distance[dist]/rank_total), 3))
                else:
                    taxonomic_distance[dist] = str(0.0)
            sys.stdout.write('\t'.join(taxonomic_distance.values()) + "\n")
    return


def clean_lineage_string(lineage):
    bad_strings = ["cellular organisms; ", "TACK group; ", "DPANN group; "]
    for bs in bad_strings:
        lineage = re.sub(bs, '', lineage)
    return lineage


def clean_classification_names(assignments):
    """
    Removes certain ranks, taxonomic super-groups, etc. so strings are more comparable
    :param assignments:
    :return:
    """
    cleaned_assignments = dict()
    for ref, queries in assignments.items():
        ref = clean_lineage_string(ref)
        cleaned_assignments[ref] = list()
        for query in queries:
            try:
                query = re.sub('_', "; ", query)
            except TypeError:
                sys.stderr.write("ERROR: Unexpected type of lineage string (" + str(type(query)) + "):\n")
                sys.stderr.write(str(query) + "\n")
                raise TypeError
            query = clean_lineage_string(query)
            cleaned_assignments[ref].append(query)
    return cleaned_assignments


def map_full_headers(fasta_headers, assignments, molecule_type):
    """
    Since the headers used throughout the TreeSAPP pipeline are truncated,
    we read the FASTA file and use those headers instead of their short version
    in case valuable information was discarded
    :param fasta_headers: full-length headers that will replace those in assignments
    :param assignments: A dictionary of reference (lineage) and query names
    :param molecule_type: prot, nucl, or rrna? Parsed from command-line arguments
    :return: assignments with a full lineage for both reference (keys) and queries (values)
    """
    genome_position_re = re.compile("^([A-Z0-9._]+):.[0-9]+-[0-9]+[_]+.*")
    full_assignments = dict()
    for ref, queries in assignments.items():
        full_assignments[ref] = list()
        for query in queries:
            database = molecule_type
            query_fields = query.split('_')
            q_accession = query_fields[0]
            if not re.search("Bacteria\|Archaea", query):
                # Check for genomic coordinates, indicating these genes will be in the 'nucleotide' database
                if genome_position_re.match(query):
                    q_accession = genome_position_re.match(query).group(1)
                    database = "dna"
                else:
                    q_accession = q_accession.split('.')[0]
                lineage = get_lineage(q_accession, database)
                full_assignments[ref].append(lineage)
            else:
                n_match = 0
                for header in fasta_headers:
                    f_accession = header[1:].split('_')[0]
                    # Useful for headers containing the full lineage
                    if q_accession == f_accession:
                        full_assignments[ref].append('_'.join(header.split('_')[1:]))
                        n_match += 1
                    else:
                        pass
                if n_match == 0:
                    sys.stderr.write("ERROR: Unable to find matching header for " + query + " in fasta_input!\n")
                    sys.exit(2)
                elif n_match > 1:
                    sys.stderr.write("ERROR: headers with identical accessions were identified in fasta_input!\n")
                    sys.stderr.write("Offending accession: " + q_accession + "\n")
                    sys.exit(3)
    return full_assignments


def main():
    args = get_arguments_()
    args.min_seq_length = 1
    # Load all the assignments from TreeSAPP output
    pre_assignments = dict()
    assignments, n_classified = read_table(args.assignment_table)
    tmp_file = args.output + os.sep + "tmp_clade_exclusion_assignments.tsv"
    if os.path.exists(tmp_file):
        pre_assignments = read_intermediate_assignments(args, tmp_file)
    fasta_headers = format_read_fasta(args.fasta_input, args.molecule, args, 1000).keys()
    # Used for sensitivity
    total_sequences = len(fasta_headers)
    sys.stdout.write("Sensitivity = " + str(round(float(n_classified / total_sequences), 3)) + "\n")
    if len(pre_assignments) == len(assignments):
        assignments = pre_assignments
    else:
        assignments = map_full_headers(fasta_headers, assignments, args.molecule)
    write_intermediate_assignments(args, assignments)
    # Get rid of some names, replace underscores with semi-colons
    assignments = clean_classification_names(assignments)
    # Load the reference lineages into a trie (prefix trie)
    taxonomic_tree = all_possible_assignments(args, args.tax_ids_table)
    # Identify at which rank each sequence's clade was excluded: [K,P,C,O,F,G,S]
    rank_assigned_dict = identify_excluded_clade(args, assignments, taxonomic_tree)
    # Determine the specificity for each rank
    determine_specificity(rank_assigned_dict)
    # Remove the intermediate file since this run completed successfully
    # if os.path.exists(tmp_file):
    #     os.remove(tmp_file)


if __name__ == "__main__":
    main()
