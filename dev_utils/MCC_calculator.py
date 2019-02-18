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
import treesapp
from phylo_dist import trim_lineages_to_rank


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


def main():
    args = get_arguments()
    log_name = args.output + os.sep + "MCC_log.txt"
    prep_logging(log_name, args.verbose)
    validate_command(args)

    ##
    # Read the file mapping reference package name to the database annotations
    ##
    pkg_name_dict = read_annotation_mapping_file(args.annot_map)

    ##
    # TODO: Run the specified taxonomic analysis tool and collect the classifications
    ##
    assignments = []
    if args.tool == "treesapp":
        ref_pkgs = ','.join(pkg_name_dict.keys())
        # TODO: Replace with a call to the treesapp main function
        classify_call = [args.treesapp + "treesapp.py", "-i", args.fasta_input,
                         "-t", ref_pkgs, "-T", str(4),
                         "-m", "prot", "--output", args.output + "TreeSAPP_output" + os.sep,
                         "--trim_align", "--overwrite"]
        launch_write_command(classify_call, False)
        classification_table = os.sep.join([args.output, "TreeSAPP_output", "final_outputs", "marker_contig_map.tsv"])
        assignments = file_parsers.read_marker_classification_table(classification_table)
    else:
        # Since you are only able to analyze a single reference package at a time with GraftM, this is ran iteratively
        for gpkg in glob.glob(args.gpkg_dir + "*gpkg"):
            pkg_name = os.path.basename(gpkg).split('.')[0]
            if pkg_name not in pkg_name_dict:
                logging.warning(gpkg + " not in " + args.annot_map + " and will be skipped...\n")
        assignments += file_parsers.read_graftm_classifications()

    sys.exit()
    test_seq_names = get_headers(args.fasta_input)

    ##
    # TODO: Bin the test sequence names into their respective confusion categories (TP, TN, FP, FN)
    ##
    # False positives: those that do not belong to the annotation matching a reference package name
    # False negatives: those whose annotations match a reference package name and were not classified
    # True positives: those with a matching annotation and reference package name and were classified
    # True negatives: those that were not classified and should not have been
    test_seq_names.clear()

    ##
    # TODO: Parse the taxonomic IDs from EggNOG headers and map taxonomic lineage information to classified sequences
    ##

    ##
    # TODO: Report the MCC score across different taxonomic distances - should increase with greater allowed distance
    ##


main()

if __name__ == "__main__":
    main()
