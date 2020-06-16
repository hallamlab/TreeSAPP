#!/usr/bin/env python3

__author__ = 'Connor Morgan-Lang'

import os
import sys
import logging
import joblib
import time
import glob

import numpy as np
import seaborn
import pandas as pd
import matplotlib
from sklearn import model_selection, svm, metrics, preprocessing, manifold

from treesapp.classy import prep_logging, PhyTrainer
from treesapp.refpkg import ReferencePackage
from treesapp.fasta import FASTA, write_new_fasta
from treesapp.commands import assign, evaluate
from treesapp.phylo_seq import assignments_to_treesaps, ItolJplace
from treesapp import MCC_calculator as ts_MCC
from treesapp import file_parsers as ts_fp
from treesapp import treesapp_args


def add_classifier_arguments(parser: treesapp_args.TreeSAPPArgumentParser) -> None:
    """
    Adds command-line arguments that are specific to *treesapp train*

    :return: None
    """
    parser.add_refpkg_opt()
    parser.add_refpkg_targets()
    treesapp_args.add_trainer_arguments(parser)

    parser.optopt.add_argument("-k", "--svm_kernel", required=False, default="lin",
                               choices=["lin", "rbf", "poly"], dest="kernel",
                               help="Specifies the kernel type to be used in the SVM algorithm."
                                    "It must be either 'lin' 'poly' or 'rbf'. [ DEFAULT = lin ]")
    parser.optopt.add_argument("--grid_search", default=False, required=False, action="store_true",
                               help="Perform a grid search across hyperparameters.")
    parser.optopt.add_argument("--tsne", default=False, required=False, action="store_true",
                               help="Generate a tSNE plot. Output will be in the same directory as the model file.")
    parser.optopt.add_argument("--classifier", required=False, choices=["occ", "bin"], default="bin",
                               help="Specify the kind of classifier to be trained: one-class classifier (OCC) or "
                                    "a binary classifier (bin).")
    parser.optopt.add_argument("--annot_map", required=False, default=None,
                               help="Path to a tabular file mapping markers being tested to their database annotations."
                                    " Columns are:  1. RefPkg name 2. Database name."
                                    " All other columns will be ignored.")
    parser.optopt.add_argument("--treesapp_output", dest="ts_out", required=False, default=None,
                               help="Path to the directory containing TreeSAPP outputs to be used for training.")

    return


def validate_command(args, trainer: PhyTrainer, sys_args) -> None:
    logging.debug("Command used:\n" + ' '.join(sys_args) + "\n")

    args.treesapp = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + os.sep
    if not args.refpkg_dir:
        args.refpkg_dir = args.treesapp + "data" + os.sep

    # Load the reference package directory if it is specified, otherwise load the package's path
    if args.pkg_path:
        target_refpkg = ReferencePackage()
        target_refpkg.f__json = args.pkg_path
        target_refpkg.slurp()
        trainer.refpkg_dir = os.path.dirname(target_refpkg.f__json) + os.sep
        trainer.target_refpkgs = [target_refpkg.prefix]
    else:
        trainer.refpkg_dir = args.refpkg_dir

    if args.targets:
        trainer.target_refpkgs = args.targets.split(',')

    if args.acc_to_lin:
        trainer.acc_to_lin = args.acc_to_lin
    else:
        trainer.acc_to_lin = trainer.var_output_dir + os.sep + "accession_id_lineage_map.tsv"
    return


def generate_training_data(ts_trainer: PhyTrainer, refpkg_dict: dict, accession_lineage_map: str,
                           trim_align=False, num_procs=2) -> dict:
    """
    Gathers the classification data for each reference package and returns a dictionary from assignments_to_treesaps.
    The dictionary has ReferencePackage.prefix keys and a list of ITolJplace for values.

    :param ts_trainer: A PhyTrainer instance
    :param refpkg_dict: A dictionary containing ReferencePackage instances indexed by their respective prefix attributes
    :param accession_lineage_map: Path to a file mapping sequence accessions to taxonomic lineages. This is to be used
    by treesapp evaluate during clade exclusion analysis.
    :param trim_align: Boolean indicating whether the multiple alignments are trimmed prior to phylogenetic placement
    :param num_procs: Number of processors/threads available to TreeSAPP and its dependencies
    :return: A dictionary returned by assignments_to_treesaps
    """
    classification_lines = []
    file_name, suffix1 = os.path.splitext(os.path.basename(ts_trainer.input_sequences))
    ts_trainer.sample_prefix = file_name
    # Run a vanilla treesapp assign on the inputs to get distance values without clade exclusion
    for name, refpkg in refpkg_dict.items():  # type: str, ReferencePackage
        # Set up the output paths
        assign_prefix = os.path.join(ts_trainer.var_output_dir, refpkg.prefix + "_assign")
        clade_exclusion_prefix = os.path.join(ts_trainer.var_output_dir, refpkg.prefix + "_ce")

        # Parameterize the commands for both treesapp assign and treesapp evaluate
        assign_params = ["-i", ts_trainer.input_sequences,
                         "-o", assign_prefix,
                         "--num_procs", str(num_procs),
                         "--refpkg_dir", os.path.dirname(refpkg.f__json),
                         "--targets", refpkg.prefix,
                         "--molecule", refpkg.molecule,
                         "--delete", "--no_svm"]
        # TODO: Allow the ranks to be controlled via command-line
        ce_params = ["-i", clade_exclusion_prefix + ".faa",
                     "-o", clade_exclusion_prefix,
                     "--molecule", refpkg.molecule,
                     "--accession2lin", accession_lineage_map,
                     "--refpkg_path", refpkg.f__json,
                     "--taxon_rank", "class",  # "order", "family", "genus", "species",
                     "--num_procs", str(num_procs),
                     "--delete"]
        if trim_align:
            assign_params.append("--trim_align")
            ce_params.append("--trim_align")

        # Run the commands if either of the output directories do not exist
        if not os.path.isdir(assign_prefix):
            assign(assign_params)
        if not os.path.isdir(clade_exclusion_prefix):
            # Create a new file for clade exclusion using just the classified sequences in case the input was unspecific
            ce_fasta_in = FASTA(
                os.path.join(assign_prefix, "final_outputs", ts_trainer.sample_prefix + "_classified.faa"))
            ce_fasta_in.load_fasta()
            ce_fasta_in.change_dict_keys("first_split")
            write_new_fasta(ce_fasta_in.fasta_dict, clade_exclusion_prefix + ".faa")

            evaluate(ce_params)

        assign_table = os.path.join(assign_prefix, "final_outputs", "marker_contig_map.tsv")
        clade_exclusion_outputs = glob.glob(os.path.join(clade_exclusion_prefix, "intermediates", '*',
                                                         "treesapp_output", "final_outputs", "marker_contig_map.tsv"))
        # Read the classification tables to gather the training data
        classification_tables = clade_exclusion_outputs
        classification_tables.append(assign_table)
        for table in classification_tables:
            classification_lines += ts_fp.read_marker_classification_table(table)

    assignments = assignments_to_treesaps(classification_lines, refpkg_dict)

    return assignments


def vectorize_placement_data(condition_names: dict, classifieds: dict, refpkg_map: dict, annot_map=None) -> np.array:
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
    :param classifieds: A dictionary of ItolJplace instances, indexed by their respective RefPkg codes (denominators)
    :param refpkg_map: Dictionary mapping refpkg names to
    :return: Summary of the distances and internal nodes for each of the false positives
    """
    features = []

    for refpkg_name, pqueries in classifieds.items():
        refpkg = refpkg_map[refpkg_name]  # type: ReferencePackage
        if not refpkg.profile_length:
            refpkg.hmm_length()
        internal_nodes = refpkg.get_internal_node_leaf_map()

        for pquery in pqueries:  # type: ItolJplace
            if annot_map:
                ref_name = annot_map[pquery.name]
            else:
                ref_name = pquery.name
            if pquery.contig_name in condition_names[ref_name]:
                # Calculate the features
                distal, pendant, avg = [round(float(x), 3) for x in pquery.distances.split(',')]
                lwr_bin = round(float(pquery.lwr), 2)
                hmm_perc = round((int(pquery.seq_len)*100)/refpkg.profile_length, 1)

                try:
                    leaf_children = len(internal_nodes[int(pquery.inode)])
                except KeyError:
                    logging.error("Unable to find internal node '%d' in the %s node-leaf map indicating a discrepancy "
                                  "between reference package versions used by treesapp assign and those used here.\n"
                                  "Was the correct output directory provided?".format(pquery.inode, pquery.name))
                    sys.exit(5)

                features.append(np.array([leaf_children, lwr_bin, distal, pendant, avg]))

    if len(features) == 0:
        return np.array(features)
    else:
        return preprocessing.normalize(np.array(features), norm='l1')


def generate_tsne(x, y, tsne_file):
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
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


def instantiate_classifier(kernel_name, occ=False):
    """
    Instantiate a particular support vector machine (SVM) classifier with a particular kernel based on the
    command-line parameter 'kernel'.

    :param kernel_name: Name of the kernel
    :param occ: Boolean flag indicating whether the desired classifier is Binary (False) or OneClass (True)
    :return: A svm.SVC instance
    """
    if kernel_name == "lin":
        if occ:
            clf = svm.OneClassSVM(kernel="linear", max_iter=1E7, tol=1E-5)
        else:
            clf = svm.LinearSVC(random_state=12345, max_iter=1E7, tol=1E-5, dual=False, C=10)  # Linear Kernel
        k_name = "linear"
    elif kernel_name == "rbf":
        if occ:
            clf = svm.OneClassSVM(kernel="rbf", tol=1E-5, gamma="auto")
        else:
            clf = svm.SVC(kernel="rbf", tol=1E-5, gamma="auto", C=100)
        k_name = "Radial Basis Function (RBF)"
    elif kernel_name == "poly":
        if occ:
            clf = svm.OneClassSVM(kernel="poly", tol=1E-5, gamma="auto")
        else:
            clf = svm.SVC(kernel="poly", tol=1E-5, gamma="auto", C=100)
        k_name = "polynomial"
    else:
        logging.error("Unknown SVM kernel '%s'.\n" % kernel_name)
        # poly_clf = svm.SVC(kernel="poly", tol=1E-3, max_iter=1E6, degree=6, gamma="auto")
        sys.exit(3)

    logging.info("Using a '%s' kernel for the SVM classifier ".format(k_name))

    return clf


def evaluate_binary_classifier(clf: svm.SVC, x_test: np.array, y_test: np.array,
                               classified_data: np.array, conditions: np.array) -> None:
    """
    Performs 10-fold cross-validation, computing F1-scores, as well as accuracy, precision and recall.
    All stats are written to the logging info stream

    :param clf: A scikit-learn classifier, either a svm.SVC or svm.LinearSVC
    :param x_test: A Numpy array containing feature vectors of the true positives
    :param y_test: A Numpy array containing feature vectors of the false positives
    :param classified_data: A Numpy array containing feature vectors for all classified data (TP and FP)
    :param conditions: A Numpy array containing the class condition of each feature vector in classified_data
    :return: None
    """
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
    return


def package_classifier(clf: svm.SVC, pickle_file: str) -> None:
    logging.info("Pickling model... ")
    try:
        pkl_handler = open(pickle_file, 'wb')
    except IOError:
        logging.error("Unable to open model file '%s' for writing.\n" % pickle_file)
        sys.exit(3)
    joblib.dump(value=clf, filename=pkl_handler)
    pkl_handler.close()
    logging.info("done.\n")
    return


def summarize_training_rank_coverage(mcc_test: ts_MCC.ConfusionTest) -> None:
    """
    Enumerate the optimal placement ranks to determine the relative distance of reference packages to test dataset

    :param mcc_test: A ConfusionTest instance
    :return: None
    """
    rank_representation = mcc_test.enumerate_optimal_rank_distances()
    logging.info("Rank coverage of clade exclusion:\n" + mcc_test.summarise_rank_coverage(rank_representation))
    return


def train_binary_classifier(tp: np.array, fp: np.array, kernel: str, tsne: str, grid_search: bool, num_procs=2):
    classified_data = np.append(fp, tp, axis=0)
    conditions = np.append(np.array([0] * len(fp)), np.array([1] * len(tp)))
    if len(conditions) != len(classified_data):
        logging.error("Inconsistent array lengths between data points (%d) and targets (%d).\n"
                      % (len(classified_data), len(conditions)))
        sys.exit(5)

    logging.info("Using %d true positives and %d false positives to train and test.\n" % (len(tp), len(fp)))

    # Split dataset into the two training and testing sets - 60% training and 40% testing
    x_train, x_test, y_train, y_test = model_selection.train_test_split(classified_data, conditions,
                                                                        test_size=0.2, random_state=12345)

    if tsne:
        generate_tsne(classified_data, conditions, tsne)

    if grid_search:
        # Test all the available kernels and return/exit
        kernels = ["lin", "rbf", "poly"]
        for k in kernels:
            evaluate_grid_scores(k, x_train=x_train, x_test=x_test, y_train=y_train, y_test=y_test, jobs=num_procs)
        return

    # Create a SVM Classifier
    clf = instantiate_classifier(kernel)  # type: svm.SVC

    logging.info("Training the classifier... ")
    # Train the model using the training sets
    clf.fit(x_train, y_train)
    logging.info("done.\n")

    evaluate_binary_classifier(clf, x_test, y_test, classified_data, conditions)
    return clf


def train_oc_classifier(tp: np.array, kernel: str):
    x_train, x_test, _, _ = model_selection.train_test_split(tp, np.array([0] * len(tp)),
                                                             test_size=0.4, random_state=12345)
    clf = instantiate_classifier(kernel, occ=True)  # type: svm.OneClassSVM
    logging.info("Training the classifier... ")
    clf.fit(x_train)
    logging.info("done.\n")

    clf.predict(x_test)
    return clf


def train_classification_filter(assignments: dict, rp_true_positives: dict, false_positives: dict, refpkg_map: dict,
                                kernel="lin", tsne="", grid_search=False, num_procs=2) -> dict:
    """
    Trains a sklearn classifier (e.g. Support Vector Machine or SVC) using the TreeSAPP assignments and writes
    the trained model to a pickle file.
    The classifier is tested using 10-fold cross-validation and reports precision, recall, accuracy and F1-scores.

    Optionally, a grid search can be performed that will automatically test mutliple different classifiers and kernels,
    not just the one specified, and return (classifier isn't written).

    :param assignments: A dictionary of ItolJplace instances, indexed by their respective RefPkg codes (denominators)
    :param rp_true_positives: A dictionary of true positive query names indexed by refpkg names
    :param false_positives: A dictionary of false positive query names indexed by refpkg names
    :param refpkg_map: A dictionary of ReferencePackage instances indexed by their respective prefix values
    :param kernel: Specifies the kernel type to be used in the SVM algorithm. Choices are 'lin' 'poly' or 'rbf'.
    [ DEFAULT = lin ]
    :param tsne: Path to a file for writing the TSNE plot. Also used to decide whether to create the TSNE plot
    :param grid_search: Flag indicating whether to perform a grid search across available classifier kernels
    :param num_procs: The number of threads to run the grid search
    :return: None
    """
    ##
    # Create internal node to leaf, HMM model length and refpkg name maps for each refpkg
    ##

    # Convert the true positive dictionary to the same format, flattening the ClassifiedSequence instances
    flattened_tp = dict()
    classifiers = dict()
    for refpkg in rp_true_positives:
        flattened_tp[refpkg] = set()
        for classified_seq in rp_true_positives[refpkg]:  # type: ts_MCC.ClassifiedSequence
            flattened_tp[refpkg].add(classified_seq.name)

    logging.info("Extracting features from TreeSAPP classifications... ")
    fp = vectorize_placement_data(condition_names=false_positives, classifieds=assignments, refpkg_map=refpkg_map)
    tp = vectorize_placement_data(condition_names=flattened_tp, classifieds=assignments, refpkg_map=refpkg_map)
    logging.info("done.\n")

    if len(fp) > 0:
        clf = train_binary_classifier(tp, fp, tsne, kernel, grid_search, num_procs)
        for refpkg in refpkg_map:  # type: str
            if refpkg in flattened_tp:
                classifiers[refpkg] = clf
    else:
        for refpkg_prefix in flattened_tp:
            refpkg = refpkg_map[refpkg_prefix]  # type: ReferencePackage
            clf = train_oc_classifier(tp, kernel)
            classifiers[refpkg.prefix] = clf

    return classifiers


def classifier_trainer(sys_args):
    """
    The classifier trainer is meant for training a classifier model within SciKit-learn on evolutionary placement data
    to improve the precision while reducing the harm of recall

    :param sys_args:
    :return:
    """
    parser = treesapp_args.TreeSAPPArgumentParser(description='Build a classifier for evolutionary placement data')
    add_classifier_arguments(parser)
    args = parser.parse_args(sys_args)

    ts_trainer = PhyTrainer()
    ts_trainer.furnish_with_arguments(args)
    ts_trainer.check_previous_output(args.overwrite)

    log_name = args.output + os.sep + "TreeSAPP_classifier_trainer_log.txt"
    if args.tsne:
        tsne_file = args.output + os.sep + "t-SNE.png"
    else:
        tsne_file = ""

    logging.info("Reading '{}'... ".format(ts_trainer.input_sequences))
    # Read all sequences that are related to the reference packages from the test data (e.g. EggNOG) into FASTA
    test_fasta = FASTA(ts_trainer.input_sequences)
    test_fasta.load_fasta()
    logging.info("done.\n")

    prep_logging(log_name, args.verbose)
    logging.info("\n##\t\t\tBeginning SVM classifier build\t\t\t##\n")
    validate_command(args, ts_trainer, sys.argv)
    ts_trainer.validate_continue(args)

    if args.annot_map:
        pkg_dbname_dict = ts_MCC.read_annotation_mapping_file(args.annot_map)
        if len(ts_trainer.target_refpkgs) == 0:
            ts_trainer.target_refpkgs = list(pkg_dbname_dict.keys())
    else:
        pkg_dbname_dict = None

    refpkg_dict = ts_fp.gather_ref_packages(ts_trainer.refpkg_dir, ts_trainer.target_refpkgs)
    if not ts_trainer.target_refpkgs:
        ts_trainer.target_refpkgs = list(refpkg_dict.keys())
    test_obj = ts_MCC.ConfusionTest(ts_trainer.target_refpkgs)
    test_obj.ref_packages = refpkg_dict
    if args.ts_out:
        test_obj.map_data(output_dir=args.ts_out, tool="treesapp")

    entrez_record_dict = ts_trainer.fetch_entrez_lineages(ref_seqs=test_fasta, molecule=ts_trainer.molecule_type,
                                                          acc_to_taxid=args.acc_to_taxid)

    assignments = generate_training_data(ts_trainer, test_obj.ref_packages, ts_trainer.acc_to_lin,
                                         args.trim_align, args.num_threads)

    test_seq_names = [seq_name[1:] if seq_name[0] == '>' else seq_name for seq_name in test_fasta.get_seq_names()]
    test_obj.num_total_queries = len(test_seq_names)
    test_obj.bin_headers(test_seq_names, assignments, pkg_dbname_dict)
    test_seq_names.clear()

    ##
    # Bin the test sequence names into their respective confusion categories (TP, TN, FP, FN)
    ##
    test_obj.populate_tax_lineage_map(entrez_record_dict)
    test_obj.map_lineages()
    test_obj.validate_false_positives()
    if args.annot_map:
        test_obj.validate_false_negatives(args.annot_map)

    summarize_training_rank_coverage(test_obj)

    classifiers = train_classification_filter(assignments, test_obj.tp, test_obj.fp, test_obj.ref_packages,
                                              args.kernel, tsne_file, args.grid_search, args.num_threads)

    for refpkg_name in classifiers:
        try:
            refpkg = refpkg_dict[refpkg_name]  # type: ReferencePackage
        except KeyError:
            logging.warning("A classifier wasn't trained for the ReferencePackage '{}'"
                            " likely because no true positives were identified.\n".format(refpkg_name))
            continue
        refpkg.svc = classifiers[refpkg_name]
        refpkg.f__json = ts_trainer.final_output_dir + os.path.basename(refpkg.f__json)
        refpkg.write_json()

    return


if __name__ == '__main__':
    classifier_trainer(sys.argv)
