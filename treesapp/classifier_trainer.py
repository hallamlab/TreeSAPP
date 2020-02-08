#!/usr/bin/env python3

__author__ = 'Connor Morgan-Lang'

import os
import sys
import logging
import re
import argparse
import numpy
import joblib
from sklearn import model_selection, svm, metrics
from treesapp import file_parsers
from treesapp.classy import prep_logging, ReferencePackage
from treesapp.fasta import get_headers
from treesapp.utilities import get_hmm_length
from treesapp import MCC_calculator


def get_arguments():
    parser = argparse.ArgumentParser(add_help=False,
                                     description="A tool for training a support vector machine from TreeSAPP output")
    required_args = parser.add_argument_group("Required arguments")
    required_args.add_argument("--fasta_input", required=True, dest="input",
                               help="FASTA-formatted file used for testing the classifiers")
    required_args.add_argument("--annot_map", required=True,
                               help="Path to a tabular file mapping markers being tested to their database annotations."
                                    " First column is the ")
    required_args.add_argument("--treesapp_output", required=False, dest="output",
                               help="Path to the `treesapp assign` output directory")

    optopt = parser.add_argument_group("Optional options")
    optopt.add_argument("-m", "--model_file", required=False, default="./treesapp_svm.pkl",
                        help="Path to a directory for writing output files")
    optopt.add_argument("-p", "--pkg_path", required=False, default=None,
                        help="The path to the TreeSAPP-formatted reference package(s) [ DEFAULT = TreeSAPP/data/ ].")

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
    return args


def validate_command(args, sys_args):
    logging.debug("Command used:\n" + ' '.join(sys_args) + "\n")

    args.treesapp = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + os.sep
    args.pkg_path = args.treesapp + "data" + os.sep

    return


def vectorize_placement_data(condition_names: dict, classifieds: list,
                             internal_nodes: dict, hmm_lengths: dict, refpkg_map: dict) -> numpy.array:
    """

    :param condition_names:
    :param classifieds:
    :param internal_nodes:
    :param hmm_lengths:
    :param refpkg_map:
    :return: Summary of the distances and internal nodes for each of the false positives
    """
    #
    features = []
    for fields in classifieds:
        _, header, refpkg, length, _, _, _, i_node, lwr, _, dists = fields
        if header in condition_names[refpkg_map[refpkg]]:
            # Calculate the features
            distal, pendant, avg = [round(float(x), 3) for x in dists.split(',')]
            hmm_perc = round((int(length)*100)/hmm_lengths[refpkg], 1)
            descendents = len(internal_nodes[refpkg][i_node])
            lwr_bin = round(float(lwr), 2)

            features.append(numpy.array([hmm_perc, descendents, lwr_bin, distal, pendant, avg]))

    return numpy.array(features)


def main():
    args = get_arguments()
    log_name = args.output + os.sep + "TreeSAPP_classifier_trainer_log.txt"
    prep_logging(log_name, args.verbose)
    logging.info("\n##\t\t\tBeginning SVM classifier build\t\t\t##\n")
    validate_command(args, sys.argv)

    pkg_name_dict = MCC_calculator.read_annotation_mapping_file(args.annot_map)
    marker_build_dict = file_parsers.parse_ref_build_params(args.treesapp, [])
    test_obj = MCC_calculator.ConfusionTest(pkg_name_dict.keys())
    test_obj.map_data(output_dir=args.output, tool="treesapp")

    ##
    # Create internal node to leaf, HMM model length and refpkg name maps for each refpkg
    ##
    internal_nodes_dict = dict()
    refpkg_map = dict()
    hmm_lengths = dict()
    for pkg_name in test_obj.ref_packages:
        refpkg = test_obj.ref_packages[pkg_name]  # type: ReferencePackage
        refpkg.prefix = marker_build_dict[pkg_name].cog
        refpkg.gather_package_files(args.pkg_path)
        refpkg_map[refpkg.prefix] = pkg_name
        internal_nodes_dict[refpkg.prefix] = MCC_calculator.internal_node_leaf_map(refpkg.tree)
        hmm_lengths[refpkg.prefix] = get_hmm_length(refpkg.profile)

    ##
    # Read the classification lines from the output
    ##
    classification_table = os.sep.join([args.output, "TreeSAPP_output", "final_outputs", "marker_contig_map.tsv"])
    if not os.path.isfile(classification_table):
        logging.error("Path to classification table '%s' doesn't exist.\n" % classification_table)
        sys.exit(3)
    classification_lines = file_parsers.read_marker_classification_table(classification_table)
    assignments = file_parsers.parse_assignments(classification_lines)

    logging.info("Reading headers in " + args.input + "... ")
    test_seq_names = [seq_name[1:] if seq_name[0] == '>' else seq_name for seq_name in get_headers(args.input)]
    logging.info("done.\n")
    test_obj.num_total_queries = len(test_seq_names)
    eggnog_re = re.compile(r"^>?(COG[A-Z0-9]+|ENOG[A-Z0-9]+)_(\d+)\..*")
    test_obj.header_regex = eggnog_re

    ##
    # Bin the test sequence names into their respective confusion categories (TP, TN, FP, FN)
    ##
    test_obj.bin_headers(test_seq_names, assignments, pkg_name_dict, marker_build_dict)
    test_seq_names.clear()
    test_obj.validate_false_positives()
    test_obj.validate_false_negatives(pkg_name_dict)

    ##
    # Convert the true positive dictionary to the same format, flattening the ClassifiedSequence instances
    ##
    flattened_tp = dict()
    for refpkg in test_obj.tp:
        flattened_tp[refpkg] = set()
        for classified_seq in test_obj.tp[refpkg]:  # type: MCC_calculator.ClassifiedSequence
            flattened_tp[refpkg].add(classified_seq.name)
    test_obj.tp = flattened_tp

    logging.info("Extracting features from TreeSAPP classifications... ")
    fp = vectorize_placement_data(condition_names=test_obj.fp, classifieds=classification_lines,
                                  hmm_lengths=hmm_lengths, internal_nodes=internal_nodes_dict, refpkg_map=refpkg_map)
    tp = vectorize_placement_data(condition_names=test_obj.tp, classifieds=classification_lines,
                                  hmm_lengths=hmm_lengths, internal_nodes=internal_nodes_dict, refpkg_map=refpkg_map)
    classified_data = numpy.append(fp, tp, axis=0)
    conditions = numpy.append(numpy.array([0]*len(fp)), numpy.array([1]*len(tp)))
    if len(conditions) != len(classified_data):
        logging.error("Inconsistent array lengths between data points (%d) and targets (%d).\n"
                      % (len(classified_data), len(conditions)))
        sys.exit(5)
    logging.info("done.\n")

    logging.info("Training the SVM with a linear kernel... ")
    # Split dataset into training set and test set - 40% training and 60% test
    x_train, x_test, y_train, y_test = model_selection.train_test_split(classified_data, conditions,
                                                                        test_size=0.6, random_state=12345)
    # Create a SVM Classifier
    clf = svm.SVC(kernel='linear')  # Linear Kernel
    # Train the model using the training sets
    clf.fit(x_train, y_train)
    logging.info("done.\n")

    logging.info("Classifying test data... ")
    # Predict the response for test dataset
    y_pred = clf.predict(x_test)
    logging.info("done.\n")
    # Model Accuracy: how often is the classifier correct?
    logging.info("Accuracy\t" + str(round(metrics.accuracy_score(y_test, y_pred), 3)) + "\n")
    # Model Precision: what percentage of positive tuples are labeled as such?
    logging.info("Precision\t" + str(round(metrics.precision_score(y_test, y_pred), 3)) + "\n")
    # Model Recall: what percentage of positive tuples are labelled as such?
    logging.info("Recall\t" + str(round(metrics.recall_score(y_test, y_pred), 3)) + "\n")

    logging.info("Pickling model... ")
    try:
        pkl_handler = open(args.model_file, 'wb')
    except IOError:
        logging.error("Unable to open model file '%s' for writing.\n" % args.model_file)
        sys.exit(3)
    joblib.dump(value=clf, filename=pkl_handler)
    pkl_handler.close()
    logging.info("done.\n")

    return


main()
