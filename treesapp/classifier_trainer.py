#!/usr/bin/env python3

__author__ = 'Connor Morgan-Lang'

import os
import sys
import logging
import re
import argparse
import joblib
import seaborn
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import model_selection, svm, metrics, preprocessing, manifold
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
                                    " Columns are:  1. RefPkg code name 2. OG 3. PFam Accession.Version"
                                    " 4. PFam Accession 5. RefPkg Name 6. Description 7. HMM Length")
    required_args.add_argument("--treesapp_output", required=False, dest="output",
                               help="Path to the `treesapp assign` output directory")

    optopt = parser.add_argument_group("Optional options")
    optopt.add_argument("-n", "--num_procs", required=False, default=2, dest="procs", type=int,
                        help="The number of parallel processes to use during grid search [ DEFAULT = 2 ]")
    optopt.add_argument("-m", "--model_file", required=False, default="./treesapp_svm.pkl",
                        help="Path to a directory for writing output files")
    optopt.add_argument("-p", "--pkg_path", required=False, default=None,
                        help="The path to the TreeSAPP-formatted reference package(s) [ DEFAULT = None ].")
    optopt.add_argument("-k", "--svm_kernel", required=False, default="lin",
                        choices=["lin", "rbf", "poly"], dest="kernel",
                        help="Specifies the kernel type to be used in the SVM algorithm."
                             "It must be either 'lin' 'poly' or 'rbf'. [ DEFAULT = lin ]")
    optopt.add_argument("--grid_search", default=False, required=False, action="store_true",
                        help="Perform a grid search across hyperparameters.")
    optopt.add_argument("--tsne", default=False, required=False, action="store_true",
                        help="Generate a tSNE plot. Output will be in the same directory as the model file.")

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
    if not args.pkg_path:
        args.pkg_path = args.treesapp + "data" + os.sep

    return


def vectorize_placement_data(condition_names: dict, classifieds: list,
                             internal_nodes: dict, hmm_lengths: dict, refpkg_map: dict) -> np.array:
    """
    Parses out the relevant fields from the *treesapp assign* classification table (marker_contig_map.tsv) to create
a numpy array from each query's data:

1. the percentage of HMM profile covered
2. number of nodes in reference phylogeny
3. number of leaf node descendents of the query's placement edge
4. the likelihood weight ratio and
5. 6, 7, distal, pendant and average distance from placement position on an edge to all leaf tips.

    :param condition_names: A dictionary of header names indexed by refpkg names. Used to determine whether the query
     belongs to this class condition (e.g. true positive, false positive)
    :param classifieds: A list of fields from the *treesapp assign* classification table
    :param internal_nodes: Dictionary mapping the internal tree nodes to descendent leaf nodes indexed by refpkg
    :param hmm_lengths: Dictionary mapping refpkg names to HMM profile lengths
    :param refpkg_map: Dictionary mapping refpkg names to
    :return: Summary of the distances and internal nodes for each of the false positives
    """
    #
    features = []
    tree_size_dict = {}
    for fields in classifieds:
        _, header, refpkg, length, _, _, _, i_node, lwr, _, dists = fields
        try:
            num_nodes = tree_size_dict[refpkg]
        except KeyError:
            num_nodes = len(internal_nodes[refpkg].values())
            tree_size_dict[refpkg] = num_nodes
        if header in condition_names[refpkg_map[refpkg]]:
            # Calculate the features
            distal, pendant, avg = [round(float(x), 3) for x in dists.split(',')]
            hmm_perc = round((int(length)*100)/hmm_lengths[refpkg], 1)
            try:
                descendents = len(internal_nodes[refpkg][i_node])
            except KeyError:
                logging.error("Unable to find internal node '%d' in the %s node-leaf map indicating a discrepancy "
                              "between reference package versions used by treesapp assign and those used here.\n"
                              "Was the correct output directory provided?" %
                              (i_node, refpkg))
                sys.exit(5)
            lwr_bin = round(float(lwr), 2)
            features.append(np.array([descendents, lwr_bin, distal, pendant, avg]))

    return preprocessing.normalize(np.array(features), norm='l1')


def generate_tsne(x, y, tsne_file):
    feat_cols = ['pixel' + str(i) for i in range(x.shape[1])]
    df = pd.DataFrame(x, columns=feat_cols)
    df['y'] = y
    df['label'] = df['y'].apply(lambda i: str(i))

    time_start = time.time()
    tsne = manifold.TSNE(n_components=2, verbose=1, perplexity=30, n_iter=300)
    tsne_results = tsne.fit_transform(df)
    print('t-SNE done! Time elapsed: {} seconds'.format(time.time() - time_start))

    df['tsne-2d-one'] = tsne_results[:, 0]
    df['tsne-2d-two'] = tsne_results[:, 1]
    plt.figure(figsize=(16, 10))
    seaborn.scatterplot(
        x="tsne-2d-one", y="tsne-2d-two",
        hue="y",
        palette=seaborn.color_palette("hls", 2),
        data=df,
        legend="full",
        alpha=0.3
    )
    plt.savefig(tsne_file, format="png")


def evaluate_grid_scores(kernel, x_train, x_test, y_train, y_test, score="f1", jobs=4):
    print("# Tuning hyper-parameters for '%s'" % kernel)
    print()

    if kernel == "lin":
        grid_params = [{'C': [0.1, 1, 10, 100, 1000]}]
        svc = svm.LinearSVC(dual=False)
    elif kernel == "rbf":
        grid_params = [{'C': [0.1, 1, 10, 100, 1000], 'gamma': [0.01, 0.001, 0.001, 0.0001]}]
        svc = svm.SVC(kernel="rbf")
    elif kernel == "poly":
        grid_params = [
            {'C': [0.1, 1, 10, 100, 1000], 'gamma': [0.01, 0.001, 0.001, 0.0001], 'degree': [1, 3, 5, 7]}]
        svc = svm.SVC(kernel="poly")
    else:
        logging.error("Unsupported kernel %s" % kernel)
        sys.exit()

    print("# Parameters: %s" % svc.get_params())

    clf = model_selection.GridSearchCV(estimator=svc, param_grid=grid_params, scoring=score, n_jobs=jobs)

    clf.fit(x_train, y_train)

    print("Best parameters set found on development set:")
    print()
    print(clf.best_params_)
    print()
    print("Grid scores on development set:")
    print()
    means = clf.cv_results_['mean_test_score']
    stds = clf.cv_results_['std_test_score']
    for mean, std, params in zip(means, stds, clf.cv_results_['params']):
        print("%0.3f (+/-%0.03f) for %r"
              % (mean, std * 2, params))
    print()

    print("Detailed classification report:")
    print()
    print("The model is trained on the full development set.")
    print("The scores are computed on the full evaluation set.")
    print()
    y_true, y_pred = y_test, clf.predict(x_test)
    print(metrics.classification_report(y_true, y_pred))
    print()


def main():
    args = get_arguments()
    log_name = args.output + os.sep + "TreeSAPP_classifier_trainer_log.txt"
    pickled_model = args.model_file
    tsne_file = '.'.join(pickled_model.split('.')[:-1]) + "_t-SNE.png"
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
    classification_table = os.sep.join([args.output, "final_outputs", "marker_contig_map.tsv"])
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
    classified_data = np.append(fp, tp, axis=0)
    conditions = np.append(np.array([0]*len(fp)), np.array([1]*len(tp)))
    if len(conditions) != len(classified_data):
        logging.error("Inconsistent array lengths between data points (%d) and targets (%d).\n"
                      % (len(classified_data), len(conditions)))
        sys.exit(5)
    logging.info("done.\n")

    logging.info("Using %d true positives and %d false positives to train and test.\n" % (len(tp), len(fp)))

    # Split dataset into the two training and testing sets - 60% training and 40% testing
    x_train, x_test, y_train, y_test = model_selection.train_test_split(classified_data, conditions,
                                                                        test_size=0.2, random_state=12345)

    if args.tsne:
        generate_tsne(classified_data, conditions, tsne_file)

    if args.grid_search:
        kernels = ["lin", "rbf", "poly"]
        for k in kernels:
            evaluate_grid_scores(k, x_train=x_train, x_test=x_test, y_train=y_train, y_test=y_test, jobs=args.procs)
        return

    # Create a SVM Classifier
    if args.kernel == "lin":
        clf = svm.LinearSVC(random_state=12345, max_iter=1E7, tol=1E-5, dual=False, C=10)  # Linear Kernel
        k_name = "linear"
    elif args.kernel == "rbf":
        clf = svm.SVC(kernel="rbf", tol=1E-5, gamma="auto", C=100)
        k_name = "Radial Basis Function (RBF)"
    elif args.kernel == "poly":
        clf = svm.SVC(kernel="poly", tol=1E-5, gamma="auto", C=100)
        k_name = "polynomial"
    else:
        logging.error("Unknown SVM kernel '%s'.\n" % args.kernel)
        # poly_clf = svm.SVC(kernel="poly", tol=1E-3, max_iter=1E6, degree=6, gamma="auto")
        sys.exit(3)

    logging.info("Training the SVM with a %s kernel... " % k_name)
    # Train the model using the training sets
    clf.fit(x_train, y_train)
    logging.info("done.\n")

    logging.info("Classifying test data... ")
    # Predict the response for test dataset
    y_pred = clf.predict(x_test)
    logging.info("done.\n")

    scores = model_selection.cross_val_score(clf, classified_data, conditions, cv=10, scoring='f1')
    logging.info("F1-score\t%0.2f (+/- %0.2f)\n" % (scores.mean(), scores.std() * 2))
    # Model Accuracy: how often is the classifier correct?
    logging.info("Accuracy\t" + str(round(metrics.accuracy_score(y_test, y_pred), 2)) + "\n")
    # Model Precision: what percentage of positive tuples are labeled as such?
    logging.info("Precision\t" + str(round(metrics.precision_score(y_test, y_pred), 2)) + "\n")
    # Model Recall: what percentage of positive tuples are labelled as such?
    logging.info("Recall\t\t" + str(round(metrics.recall_score(y_test, y_pred), 2)) + "\n")

    logging.info("Pickling model... ")
    try:
        pkl_handler = open(pickled_model, 'wb')
    except IOError:
        logging.error("Unable to open model file '%s' for writing.\n" % pickled_model)
        sys.exit(3)
    joblib.dump(value=clf, filename=pkl_handler)
    pkl_handler.close()
    logging.info("done.\n")

    return


main()
