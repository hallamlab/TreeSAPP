#!/usr/bin/env python3

__author__ = 'Connor Morgan-Lang'

import argparse
import sys
import os
import inspect
import shutil
import glob

cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile(inspect.currentframe()))[0]))
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)
sys.path.insert(0, cmd_folder + os.sep + ".." + os.sep)
from fasta import format_read_fasta, write_new_fasta, get_headers, read_fasta_to_dict
from create_treesapp_refpkg import get_header_format, finalize_ref_seq_lineages
from utilities import return_sequence_info_groups, find_executables
from external_command_interface import launch_write_command
import file_parsers
from classy import prep_logging, get_header_info, register_headers, MarkerTest
from entrez_utils import *


class ClassifiedSequence:
    def __init__(self, header):
        self.name = header
        self.ref = ""
        self.ncbi_tax = ""
        self.assigned_lineage = ""
        self.true_lineage = ""
        self.tax_dist = 0


class ConfusionTest:
    def __init__(self, gene_list):
        self.header_regex = None
        self.test_markers = gene_list
        self.fn = {key: [] for key in gene_list}
        self.fp = {key: [] for key in gene_list}
        self.tp = {key: [] for key in gene_list}  # This will be a list of ClassifiedSequence instances
        self.tax_lineage_map = dict()
        self.dist_wise_tp = dict()
        self.tn = 0

    def get_info(self):
        info_string = ""
        return info_string

    def marker_classifcation_summary(self, g_name):
        """
        Provide a classification summary for a specific marker gene, g_name
        :param g_name:
        :return:
        """
        summary_string = ""
        return summary_string

    def retrieve_lineages(self, group=1):
        if not self.header_regex:
            logging.error("Unable to parse taxonomic identifiers from header without a regular expression.\n")
            sys.exit(19)

        # Gather the unique taxonomy IDs and store in EntrezRecord instances
        entrez_records = list()
        for marker in self.tp:
            for header in self.tp[marker]:
                try:
                    tax_id = self.header_regex.search(header).group(group)
                except (KeyError, TypeError):
                    logging.error("Header '" + header + "' in test FASTA file does not match the supported format.\n")
                    sys.exit(5)

                # Only make Entrez records for new NCBI taxonomy IDs
                if tax_id not in self.tax_lineage_map:
                    e_record = EntrezRecord(header, "")
                    e_record.ncbi_tax = tax_id
                    entrez_records.append(e_record)

        # Query the Entrez database for these unique taxonomy IDs
        fetch_lineages_from_taxids(entrez_records)

        for e_record in entrez_records:  # type: EntrezRecord
            self.tax_lineage_map[e_record.ncbi_tax] = e_record.lineage
        return

    def bin_true_positives_by_taxdist(self):
        for marker in self.tp:
            for tp_inst in self.tp[marker]:  # type: ClassifiedSequence
                tp_inst.tax_dist = taxonomic_distance(tp_inst.assigned_lineage, self.tax_lineage_map[tp_inst.ncbi_tax])
                try:
                    self.dist_wise_tp[tp_inst.tax_dist] += 1
                except KeyError:
                    self.dist_wise_tp[tp_inst.tax_dist] = 1
        return

    def get_true_positives_at_dist(self, max_distance=7, g_name=None):
        total_tp = 0
        for tax_dist in sorted(self.dist_wise_tp, key=int):
            if tax_dist <= max_distance:
                total_tp += self.dist_wise_tp[tax_dist]
        return total_tp

    def get_true_negatives(self, all_queries: int):
        acc = 0
        for marker in self.test_markers:
            acc += len(self.fn[marker])
            acc += len(self.fp[marker])
            acc += len(self.tp[marker])
        self.tn = all_queries - acc

    def get_false_positives(self, g_name=None):
        if g_name:
            return len(self.fp[g_name])
        else:
            return sum([len(self.fp[marker]) for marker in self.test_markers])

    def get_false_negatives(self, g_name=None):
        if g_name:
            return len(self.fn[g_name])
        else:
            return sum([len(self.fn[marker]) for marker in self.test_markers])


def get_arguments():
    parser = argparse.ArgumentParser(add_help=False,
                                     description="A wrapper script for calculating Matthews correlation coefficient for"
                                                 "TreeSAPP or GraftM. Currently only supports testing proteins.")
    required_args = parser.add_argument_group("Required arguments")
    required_args.add_argument("--fasta_input", required=True,
                               help="FASTA-formatted file used for testing the classifiers")
    # required_args.add_argument("--reference_marker", required=True,
    #                            help="Short-form name of the marker gene to be tested (e.g. mcrA, pmoA, nosZ)")
    required_args.add_argument("--annot_map", required=True,
                               help="Path to a tabular file mapping markers being tested to their database annotations."
                                    " First column is the ")

    optopt = parser.add_argument_group("Optional options")
    optopt.add_argument("--tool", default="treesapp", required=False,
                        choices=["treesapp", "graftm", "diamond"],
                        help="Classify using one of the tools: treesapp [DEFAULT], graftm, or diamond.")
    optopt.add_argument("--gpkg_dir", required=False, default=None,
                        help="Path to a directory containing all GraftM reference packages to test."
                             " Files must follow 'name.gpkg' scheme and 'name' is in the first column of --annot_map")
    optopt.add_argument("--output", required=False, default="./MCC_output/",
                        help="Path to a directory for writing output files")
    # optopt.add_argument('-m', '--molecule',
    #                     help='the type of input sequences (prot = Protein [DEFAULT]; dna = Nucleotide; rrna = rRNA)',
    #                     default='prot',
    #                     choices=['prot', 'dna', 'rrna'])
    # optopt.add_argument("-l", "--length",
    #                     required=False, type=int, default=0,
    #                     help="Arbitrarily slice the input sequences to this length. "
    #                          "Useful for testing classification accuracy for fragments.")

    miscellaneous_opts = parser.add_argument_group("Miscellaneous options")
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

    if sys.version_info > (2, 9):
        args.py_version = 3
    else:
        args.py_version = 2

    if args.gpkg_dir and args.gpkg_dir[-1] != os.sep:
        args.gpkf_dir += os.sep
    return args


def validate_command(args):
    if args.overwrite:
        if os.path.exists(args.output):
            shutil.rmtree(args.output)

    if args.tool in ["diamond", "graftm"]:
        if not args.gpkg_dir:
            logging.error(args.tool + " specified but not GraftM reference package directory specified.\n")
            sys.exit(17)
        elif not os.path.isdir(args.gpkg_dir):
            logging.error(args.gpkg_dir + " GraftM reference package directory does not exist!\n")
            sys.exit(17)
        elif len(glob.glob(args.gpkg_dir + "*gpkg")) == 0:
            logging.error("Not GraftM reference packages found in " + args.gpkg_dir + ".\n")
            sys.exit(17)
    elif args.tool == "treesapp" and args.gpkg_dir:
        logging.warning("--gpkg specific but tool selected is 'treesapp'... it will be ignored.\n")

    return


def read_annotation_mapping_file(annot_map_file):
    annot_map = dict()
    try:
        annot_map_handler = open(annot_map_file)
    except IOError:
        logging.error("Unable to open " + annot_map_file + " for reading!\n")
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
            fields = line.split("\t")
            try:
                annot_map[fields[0]] = fields[1].split(',')
            except KeyError:
                logging.error("Insufficient number of fields on line " + str(n) + " in " + annot_map_file + "!\n" +
                              "File must have the reference package name and the database name in"
                              " the first two columns, respectively. Any number of columns can follow.\n")
                sys.exit(9)

    annot_map_handler.close()
    return annot_map


def bin_headers(test_seq_names, assignments, annot_map, confusion_obj: ConfusionTest):
    """
    Function for sorting/binning the classified sequences at T/F positives/negatives based on the
    :param test_seq_names: List of all headers in the input FASTA file
    :param assignments: Dictionary mapping taxonomic lineages to a list of headers that were classified to this lineage
    :param annot_map: Dictionary mapping reference package (gene) name keys to database names values
    :return:
    """
    # False positives: those that do not belong to the annotation matching a reference package name
    # False negatives: those whose annotations match a reference package name and were not classified
    # True positives: those with a matching annotation and reference package name and were classified
    # True negatives: those that were not classified and should not have been

    mapping_dict = dict()
    for refpkg in annot_map:
        orthos = annot_map[refpkg].split(',')
        for gene in orthos:
            mapping_dict[gene] = refpkg

    # test the format of the header in test_seq_names
    eggnog_re = re.compile(r"^(COG[A-Z0-9]+|ENOG[A-Z0-9]+)_(\d+)\..*")
    positive_queries = dict()
    for header in test_seq_names:
        try:
            ref_g, tax_id = eggnog_re.match(header).groups()
        except (TypeError, KeyError):
            logging.error("Header '" + header + "' in test FASTA file does not match the supported format.\n")
            sys.exit(5)
        # Keep the name in
        if ref_g in mapping_dict:
            marker = mapping_dict[ref_g]
            if marker not in positive_queries:
                positive_queries[marker] = []
            positive_queries[marker].append(header)

    for marker in assignments:
        positives = set(positive_queries[marker])
        true_positives = set()
        for tax_lin in assignments:
            classified_seqs = assignments[marker][tax_lin]
            for seq_name in classified_seqs:
                try:
                    ref_g, tax_id = eggnog_re.match(seq_name).groups()
                except (TypeError, KeyError):
                    logging.error("Classified sequence name '" + seq_name + "' does not match the supported format.\n")
                    sys.exit(5)
                if annot_map[marker] == ref_g:
                    # Populate the relevant information for the classified sequence
                    tp_inst = ClassifiedSequence(seq_name)
                    tp_inst.ncbi_tax = tax_id
                    tp_inst.ref = marker
                    tp_inst.assigned_lineage = tax_lin

                    # Add the True Positive to the relevant collections
                    confusion_obj.tp[marker].append(tp_inst)
                    true_positives.add(seq_name)
                else:
                    confusion_obj.fp[marker].append(seq_name)

        # Identify the False Negatives using set difference - those that were not classified but should have been
        confusion_obj.fn[marker] = list(positives.difference(true_positives))
    return


def main():
    args = get_arguments()
    log_name = args.output + os.sep + "MCC_log.txt"
    prep_logging(log_name, args.verbose)
    validate_command(args)

    ##
    # Read the file mapping reference package name to the database annotations
    ##
    pkg_name_dict = read_annotation_mapping_file(args.annot_map)
    test_obj = ConfusionTest(pkg_name_dict.keys())

    ##
    # TODO: Run the specified taxonomic analysis tool and collect the classifications
    ##
    assignments = {}
    test_fa_prefix = '.'.join(args.fasta_input.split('.')[:-1])
    if args.tool == "treesapp":
        ref_pkgs = ','.join(pkg_name_dict.keys())
        # TODO: Replace with a call to the treesapp main function
        classify_call = [args.treesapp + "treesapp.py", "-i", args.fasta_input,
                         "-t", ref_pkgs, "-T", str(4),
                         "-m", "prot", "--output", args.output + "TreeSAPP_output" + os.sep,
                         "--trim_align", "--overwrite"]
        launch_write_command(classify_call, False)
        classification_table = os.sep.join([args.output, "TreeSAPP_output", "final_outputs", "marker_contig_map.tsv"])
        classification_lines = file_parsers.read_marker_classification_table(classification_table)
        assignments = file_parsers.parse_assignments(classification_lines)
    else:
        # Since you are only able to analyze a single reference package at a time with GraftM, this is ran iteratively
        for gpkg in glob.glob(args.gpkg_dir + "*gpkg"):
            pkg_name = os.path.basename(gpkg).split('.')[0]
            if pkg_name not in pkg_name_dict:
                logging.warning(gpkg + " not in " + args.annot_map + " and will be skipped...\n")
                continue
            output_dir = os.sep.join([args.output, "GraftM_output", pkg_name]) + os.sep
            classify_call = ["graftM", "graft",
                             "--forward", args.fasta_input,
                             "--graftm_package", gpkg,
                             "--input_sequence_type", "aminoacid",
                             "--threads", str(8),
                             "--output_directory", output_dir,
                             "--force"]
            launch_write_command(classify_call, False)
            classification_table = output_dir + test_fa_prefix + os.sep + test_fa_prefix + "_read_tax.tsv"
            assignments[pkg_name] = file_parsers.read_graftm_classifications(classification_table)

    test_seq_names = get_headers(args.fasta_input)

    ##
    # Bin the test sequence names into their respective confusion categories (TP, TN, FP, FN)
    ##
    bin_headers(test_seq_names, assignments, pkg_name_dict, test_obj)
    test_seq_names.clear()

    ##
    # Parse the taxonomic IDs from EggNOG headers and map taxonomic lineage information to classified sequences
    ##
    test_obj.retrieve_lineages(2)
    test_obj.bin_true_positives_by_taxdist()

    ##
    # TODO: Report the MCC score across different taxonomic distances - should increase with greater allowed distance
    ##
    # d = 0
    # while d < 7:
    #     num_tp = test_obj.get_true_positives_at_dist(d)
    #     mcc = calculate_matthews_correlation_coefficient(num_tp,
    #                                                      test_obj.get_true_negatives(),
    #                                                      test_obj.get_false_negatives()
    #                                                      test_obj.get_false_positives())
    #     d += 1


main()

if __name__ == "__main__":
    main()
