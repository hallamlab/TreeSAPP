
import os
import logging
import sys
import time
import itertools
from glob import glob

from ete3 import Tree
from tqdm import tqdm
import numpy as np
import seaborn
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from sklearn import model_selection, svm, metrics, preprocessing, manifold

from treesapp import file_parsers
from treesapp import wrapper
from treesapp import fasta
from treesapp.phylo_seq import PQuery
from treesapp.external_command_interface import launch_write_command, create_dir_from_taxon_name
from treesapp.jplace_utils import jplace_parser
from treesapp.entish import map_internal_nodes_leaves
from treesapp.phylo_dist import parent_to_tip_distances
from treesapp.refpkg import ReferencePackage


def rarefy_rank_distances(rank_distances: dict) -> dict:
    """
    The number of observations (phylogenetic distances) for each key (taxonomic rank) are rarefied to
    number of observations found for the rank with the fewest observations.
    First, the minimum number is identified by finding the smallest list in the rank_distances.values().
    Then observations from the input dictionary are randomly copied into a new dictionary for each rank.

    :param rank_distances: A dictionary of floats indexed by taxonomic rank
    :return: Dictionary of floats indexed by taxonomic rank
    """
    rarefied_dists = dict()
    min_samples = min([len(rank_distances[rank]) for rank in rank_distances])
    for rank in rank_distances:
        slist = sorted(rank_distances[rank])
        if len(slist) == min_samples:
            rarefied_dists[rank] = slist
        else:
            rarefied_dists[rank] = list()
            i = 0
            while i < min_samples:
                rarefied_dists[rank].append(slist.pop(np.random.randint(0, len(slist))))
                i += 1
    return rarefied_dists


def generate_pquery_data_for_trainer(ref_pkg: ReferencePackage, taxon: str,
                                     test_fasta: fasta.FASTA, training_seqs: list, rank: str,
                                     executables: dict, output_dir: str, pbar: tqdm, num_threads=2) -> list:
    """
    An all-in-one function for performing phylogenetic placement using query sequences from an excluded clade.
    The original ReferencePackage instance is cloned (so it is not modified) and query sequences belonging to taxon
    are removed from it.

    :param ref_pkg: A ReferencePackage instance that has not been modified
    :param taxon: A taxonomic lineage of the query sequences that will be removed and placed
    :param test_fasta: A FASTA instance containing the query sequences to be placed
    :param training_seqs:  A list of the query sequence names in test_fasta to be placed.
    :param rank: The taxonomic rank the taxon represents
    :param executables: A dictionary of executable names mapped to their respective absolute paths
    :param output_dir: Path to a directory to write the output files
    :param pbar: A tqdm.tqdm progress bar instance
    :param num_threads: The number of threads to use during phylogenetic placement
    :return: A list of PQuery instances representing the best phylogenetic placement for each query sequence
    """
    intermediate_files = list()
    pqueries = list()
    taxonomy_filtered_query_seqs = dict()
    query_seq_name_map = dict()
    # Clean up the query taxon's name
    taxon_test_dir = create_dir_from_taxon_name(taxon, output_dir)
    query_name = os.path.split(taxon_test_dir[:-1])[1]

    # Create the cloned ReferencePackage to be used for this taxon's trials
    clade_exclusion_json = taxon_test_dir + ref_pkg.prefix + ref_pkg.refpkg_suffix
    ce_refpkg = ref_pkg.clone(clade_exclusion_json)
    # Remove all sequences belonging to a taxonomic rank from the reference package
    ce_refpkg.exclude_clade_from_ref_files(tmp_dir=taxon_test_dir, target_clade=taxon, executables=executables)

    ce_fasta = fasta.FASTA(ce_refpkg.f__msa)
    ce_fasta.load_fasta()
    ce_tree = Tree(ce_refpkg.f__tree)

    # Paths to temporary files
    query_fasta_file = taxon_test_dir + "queries.fa"
    query_sto_file = os.path.splitext(query_fasta_file)[0] + ".sto"
    all_msa = os.path.splitext(query_fasta_file)[0] + ".mfa"
    query_msa = taxon_test_dir + "queries.mfa"
    ref_msa = taxon_test_dir + "references.mfa"
    intermediate_files += [query_fasta_file, query_sto_file, all_msa, query_msa, ref_msa]

    # Write query FASTA containing sequences belonging to `taxon`
    query_seq_decrementor = -1
    for seq_name in training_seqs:
        query_seq_name_map[query_seq_decrementor] = seq_name
        taxonomy_filtered_query_seqs[str(query_seq_decrementor)] = test_fasta.fasta_dict[seq_name]
        query_seq_decrementor -= 1
    logging.debug("\t{} query sequences.\n".format(len(taxonomy_filtered_query_seqs.keys())))

    fasta.write_new_fasta(taxonomy_filtered_query_seqs, fasta_name=query_fasta_file)

    ##
    # Run hmmalign, BMGE and EPA-NG to map sequences from the taxonomic rank onto the tree
    ##
    aln_stdout = wrapper.profile_aligner(executables, ce_refpkg.f__msa, ce_refpkg.f__profile,
                                         query_fasta_file, query_sto_file)
    # Reformat the Stockholm format created by cmalign or hmmalign to FASTA
    sto_dict = file_parsers.read_stockholm_to_dict(query_sto_file)
    fasta.write_new_fasta(sto_dict, all_msa)

    logging.debug(str(aln_stdout) + "\n")

    trim_command, combined_msa = wrapper.get_msa_trim_command(executables, all_msa, ce_refpkg.molecule)
    launch_write_command(trim_command)
    intermediate_files += glob(combined_msa + "*")

    # Ensure reference sequences haven't been removed during MSA trimming
    msa_dict, failed_msa_files, summary_str = file_parsers.validate_alignment_trimming([combined_msa],
                                                                                       set(ce_fasta.fasta_dict),
                                                                                       True)
    nrow, ncolumn = fasta.multiple_alignment_dimensions(mfa_file=combined_msa,
                                                        seq_dict=fasta.read_fasta_to_dict(combined_msa))
    logging.debug("Columns = " + str(ncolumn) + "\n")
    if combined_msa not in msa_dict.keys():
        logging.debug("Placements for '{}' are being skipped after failing MSA validation.\n".format(taxon))
        for old_file in intermediate_files:
            os.remove(old_file)
            intermediate_files.clear()
        return pqueries
    logging.debug("Number of sequences discarded: " + summary_str + "\n")

    # Create the query-only FASTA file required by EPA-ng
    fasta.split_combined_ref_query_fasta(combined_msa, query_msa, ref_msa)

    raxml_files = wrapper.raxml_evolutionary_placement(epa_exe=executables["epa-ng"],
                                                       refpkg_tree=ce_refpkg.f__tree,
                                                       refpkg_msa=ref_msa,
                                                       refpkg_model=ce_refpkg.f__model_info,
                                                       query_msa=query_msa, query_name=query_name,
                                                       output_dir=output_dir, num_threads=num_threads)

    # Parse the JPlace file to pull distal_length+pendant_length for each placement
    jplace_data = jplace_parser(raxml_files["jplace"])
    placement_tree = jplace_data.tree
    node_map = map_internal_nodes_leaves(placement_tree)
    for pquery in jplace_data.placements:
        top_lwr = 0.1
        top_placement = PQuery(taxon, rank)
        top_placement.ref_name = ref_pkg.prefix
        for name, info in pquery.items():
            if name == 'p':
                for placement in info:
                    # Only record the best placement's distance
                    lwr = float(placement[2])
                    if lwr > top_lwr:
                        top_lwr = lwr
                        top_placement.inode = placement[0]
                        top_placement.likelihood = placement[1]
                        top_placement.lwr = lwr
                        top_placement.distal = round(float(placement[3]), 6)
                        top_placement.pendant = round(float(placement[4]), 6)
                        leaf_children = node_map[int(top_placement.inode)]
                        if len(leaf_children) > 1:
                            # Reference tree with clade excluded
                            parent = ce_tree.get_common_ancestor(leaf_children)
                            tip_distances = parent_to_tip_distances(parent, leaf_children)
                            top_placement.mean_tip = round(float(sum(tip_distances) / len(tip_distances)), 6)
            elif name == 'n':
                top_placement.place_name = query_seq_name_map[int(info.pop())]
            else:
                logging.error("Unexpected variable in pquery keys: '" + name + "'\n")
                sys.exit(33)

        if top_placement.lwr >= 0.5:  # The minimum likelihood weight ration a placement requires to be included
            pqueries.append(top_placement)

    # Remove intermediate files from the analysis of this taxon
    intermediate_files += list(raxml_files.values())
    for old_file in intermediate_files:
        if os.path.isfile(old_file):
            os.remove(old_file)

    pbar.update(len(taxonomy_filtered_query_seqs.keys()))

    return pqueries


def augment_training_set(row: np.array, n_reps=3, feature_scale=0.2):
    n_features = len(row)
    # make copies of row
    training_set = np.array([])
    for _ in range(n_reps):
        # create vector of random gaussians
        gauss = np.random.normal(loc=0.0, scale=feature_scale, size=len(row))
        # add to test case
        new_row = row + gauss
        # store in test set
        training_set = np.append(training_set, [new_row.clip(min=0)])
    return training_set.reshape(n_reps, n_features)


def vectorize_placement_data(condition_names: dict, classifieds: dict, refpkg_map: dict,
                             annot_map=None) -> (np.array, list):
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
    :param classifieds: A dictionary of JPlace instances, indexed by their respective RefPkg codes (denominators)
    :param refpkg_map: A dictionary of ReferencePackage instances indexed by their respective prefix values
    :param annot_map: A dictionary mapping ReferencePackage.prefix values to the names of database sequence names
    :return: Summary of the distances and internal nodes for each of the false positives
    """
    rank_feature_vectors = {}
    feature_vectors = np.array([])
    pqueries_used = []

    # Generate feature vectors for each rank and index them in the rank_feature_vectors dictionary
    for refpkg_name, pqueries in classifieds.items():
        refpkg = refpkg_map[refpkg_name]  # type: ReferencePackage
        if not refpkg.profile_length:
            refpkg.hmm_length()
        internal_nodes = refpkg.get_internal_node_leaf_map()

        for pquery in pqueries:  # type: PQuery
            if annot_map:
                ref_name = annot_map[pquery.ref_name]
            else:
                ref_name = pquery.ref_name
            if pquery.place_name in condition_names[ref_name]:
                # Calculate the features
                try:
                    distal, pendant, avg = pquery.distal, pquery.pendant, pquery.mean_tip
                except AttributeError:
                    distal, pendant, avg = [round(float(x), 3) for x in pquery.distances.split(',')]

                lwr_bin = round(float(pquery.lwr), 2)
                # hmm_perc = round((int(pquery.seq_len)*100)/refpkg.profile_length, 1)

                try:
                    leaf_children = len(internal_nodes[int(pquery.inode)])
                except KeyError:
                    logging.error("Unable to find internal node '%d' in the %s node-leaf map indicating a discrepancy "
                                  "between reference package versions used by treesapp assign and those used here.\n"
                                  "Was the correct output directory provided?".format(pquery.inode, pquery.ref_name))
                    sys.exit(5)

                raw_array = np.array([leaf_children, lwr_bin, distal, pendant, avg])
                try:
                    rank_feature_vectors[pquery.rank].append(raw_array)
                except KeyError:
                    rank_feature_vectors[pquery.rank] = [raw_array]
                pqueries_used.append(pquery)

    # Match the normalized feature vectors with their respective PQuery instance and set PQuery.feature_vec
    if len(rank_feature_vectors) == 0:
        return np.array(feature_vectors), pqueries_used
    else:
        for rank in rank_feature_vectors:
            if len(rank_feature_vectors[rank]) > 0:
                rank_feature_vectors[rank] = preprocessing.normalize(rank_feature_vectors[rank], norm='l1')
                i = 0
                end = len(rank_feature_vectors[rank])
                while i < end:
                    pqueries_used[i].feature_vec = rank_feature_vectors[rank][i]
                    i += 1

    # Augment the training data based on the number of samples available at each rank
    limit = max(len(rank_feature_vectors[r]) for r in rank_feature_vectors)
    for rank in rank_feature_vectors:
        for sample in rank_feature_vectors[rank]:
            if len(rank_feature_vectors[rank]) < limit:
                synthetic_samples = augment_training_set(sample)
                rank_feature_vectors[rank] = np.append(rank_feature_vectors[rank], synthetic_samples, axis=0)
            else:
                break
        # TODO: replace this, seems hacky
        try:
            feature_vectors = np.append(feature_vectors, rank_feature_vectors[rank], axis=0)
        except ValueError:
            feature_vectors = rank_feature_vectors[rank]

    return feature_vectors, pqueries_used


def generate_tsne(x, y, tsne_file):
    matplotlib.use('Agg')
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
        logging.error("Unknown SVM kernel '{}'.\n".format(kernel_name))
        # poly_clf = svm.SVC(kernel="poly", tol=1E-3, max_iter=1E6, degree=6, gamma="auto")
        sys.exit(3)

    logging.info("Using a '{}' kernel for the SVM classifier\n".format(k_name))

    return clf


def evaluate_classifier(clf, x_test: np.array, y_test: np.array,
                        classified_data: np.array, conditions: np.array) -> np.array:
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

    return y_pred


def summarize_pquery_ranks(pqueries: list) -> None:
    rank_counts = {}
    n_skipped = 0
    for pq in pqueries:  # type: PQuery
        try:
            rank_counts[pq.rank] += 1
        except KeyError:
            rank_counts[pq.rank] = 1
        except AttributeError:
            n_skipped += 1

    logging.info("\n".join(["{}\t{}".format(rank, rank_counts[rank]) for
                            rank in sorted(rank_counts, key=lambda x: rank_counts[x])]) + "\n")

    logging.debug("{} PQuery instances not included in rank summary.\n".format(n_skipped))

    return


def characterize_predictions(preds: np.array, pquery_test_map: dict) -> None:
    # Summarise the number of samples for each rank
    logging.info("Rank representation across clade-exclusion samples tested:\n")
    summarize_pquery_ranks(list(itertools.chain.from_iterable(pquery_test_map.values())))

    # Summarise the number of false negatives for each rank
    fn_pqueries = []
    for i in pquery_test_map:
        if preds[i] == -1:
            fn_pqueries += pquery_test_map[i]
    logging.info("Rank representation across clade-exclusion false negatives:\n")
    summarize_pquery_ranks(fn_pqueries)
    return


def map_test_indices_to_pqueries(feature_vectors, pqueries) -> dict:
    """

    :param feature_vectors:
    :param pqueries:
    :return:
    """
    index_pquery_map = {}
    i = 0
    while i < len(feature_vectors):
        vec = feature_vectors[i]  # type: np.ndarray
        index_pquery_map[i] = []
        for pquery in pqueries:  # type: PQuery
            if np.array_equal(vec, pquery.feature_vec):
                index_pquery_map[i].append(pquery)
        i += 1
    return index_pquery_map


def generate_train_test_data(true_ps: np.array, true_pqueries: list, test_pr=0.4, false_ps=None):
    """
    Used for splitting the total dataset into training ana testing arrays and mapping indices
    in the positive test fraction to PQueries. This map is to be later used for characterizing predictions by
    rank and lineage.

    Conditions [0|1] are set for the dataset based on whether they are false positives (0) or true positives (1).
    If no false positives are present, a numpy.array of length true_ps is used.


    :param true_ps:
    :param true_pqueries:
    :param test_pr:
    :param false_ps:
    :return:
    """
    if false_ps:
        classified_data = np.append(false_ps, true_ps, axis=0)
        conditions = np.append(np.array([0] * len(false_ps)),
                               np.array([1] * len(true_ps)))
    else:
        classified_data = true_ps
        conditions = np.array([1] * len(true_ps))

    if len(conditions) != len(classified_data):
        logging.error("Inconsistent array lengths between data points (%d) and targets (%d).\n"
                      % (len(classified_data), len(conditions)))
        sys.exit(5)

    # Split dataset into the two training and testing sets - 60% training and 40% testing
    x_train, x_test, y_train, y_test = model_selection.train_test_split(classified_data, conditions,
                                                                        test_size=test_pr, random_state=12345)

    index_pquery_map = map_test_indices_to_pqueries(x_test, true_pqueries)

    return classified_data, conditions, x_train, x_test, y_train, y_test, index_pquery_map


def train_binary_classifier(x_train: np.array, y_train: np.array, kernel: str):
    # Create a SVM Classifier
    clf = instantiate_classifier(kernel)  # type: svm.SVC

    logging.info("Training the classifier... ")
    # Train the model using the training sets
    clf.fit(x_train, y_train)
    logging.info("done.\n")

    return clf


def train_oc_classifier(x_train: np.array, kernel: str):
    clf = instantiate_classifier(kernel, occ=True)  # type: svm.OneClassSVM

    logging.info("Training the classifier... ")
    clf.fit(x_train)
    logging.info("done.\n")

    return clf


def train_classification_filter(assignments: dict, true_pos: dict, refpkg_map: dict, kernel: str,
                                tsne="", grid_search=False, false_pos=None, num_procs=2) -> (dict, dict):
    """
    Trains a sklearn classifier (e.g. Support Vector Machine) using the TreeSAPP assignments and writes
    the trained model to a pickle file.
    The classifier is tested using 10-fold cross-validation and reports precision, recall, accuracy and F1-scores.

    Optionally, a grid search can be performed that will automatically test multiple different classifiers and kernels,
    not just the one specified, and return (classifier isn't written).

    :param assignments: A dictionary of JPlace instances, indexed by their respective RefPkg codes (denominators)
    :param true_pos: A dictionary of true positive query names indexed by refpkg names
    :param false_pos: A dictionary of false positive query names indexed by refpkg names
    :param refpkg_map: A dictionary of ReferencePackage instances indexed by their respective prefix values
    :param kernel: Specifies the kernel type to be used in the SVM algorithm. Choices are 'lin' 'poly' or 'rbf'.
    [ DEFAULT = lin ]
    :param tsne: Path to a file for writing the TSNE plot. Also used to decide whether to create the TSNE plot
    :param grid_search: Flag indicating whether to perform a grid search across available classifier kernels
    :param num_procs: The number of threads to run the grid search
    :return: None
    """
    classifiers = dict()

    logging.info("Extracting features from TreeSAPP classifications... ")
    tp, t_pqs = vectorize_placement_data(condition_names=true_pos, classifieds=assignments, refpkg_map=refpkg_map)
    if false_pos:
        fp, f_pqs = vectorize_placement_data(condition_names=false_pos, classifieds=assignments, refpkg_map=refpkg_map)
    else:
        fp, f_pqs = np.array([]), []
    logging.info("done.\n")

    if len(tp) == 0:
        logging.error("No true positives are available for training.\n")
        sys.exit(17)
    elif len(fp) > 0:
        classified_data, conditions, x_train, x_test, y_train, y_test, pquery_test_map = generate_train_test_data(tp, t_pqs, false_ps=fp)
        if tsne:
            generate_tsne(classified_data, conditions, tsne)

        if grid_search:
            # Test all the available kernels and return/exit
            kernels = ["lin", "rbf", "poly"]
            for k in kernels:
                evaluate_grid_scores(k, x_train=x_train, x_test=x_test, y_train=y_train, y_test=y_test, jobs=num_procs)
            return

        clf = train_binary_classifier(tp, fp, kernel)

        for refpkg_prefix in refpkg_map:  # type: str
            if refpkg_prefix in true_pos:
                classifiers[refpkg_prefix] = clf
                preds = evaluate_classifier(clf, x_test, y_test, classified_data, conditions)
                characterize_predictions(preds, pquery_test_map)
    else:
        for refpkg_prefix in true_pos:
            classified_data, conditions, x_train, x_test, y_train, y_test, pquery_test_map = generate_train_test_data(tp, t_pqs)
            clf = train_oc_classifier(x_train, kernel)
            classifiers[refpkg_prefix] = clf
            preds = evaluate_classifier(clf, x_test, y_test, classified_data, conditions)
            characterize_predictions(preds, pquery_test_map)

    return classifiers
