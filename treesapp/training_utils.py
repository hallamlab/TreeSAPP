
import os
import logging
import sys
import time
import itertools
from glob import glob

from tqdm import tqdm
import numpy as np
import seaborn
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from sklearn import model_selection, svm, metrics, manifold
from collections import namedtuple

from treesapp import file_parsers
from treesapp import wrapper
from treesapp import fasta
from treesapp import classy
from treesapp import phylo_seq
from treesapp.external_command_interface import launch_write_command, create_dir_from_taxon_name
from treesapp.jplace_utils import jplace_parser, demultiplex_pqueries, calc_pquery_mean_tip_distances
from treesapp.entish import map_internal_nodes_leaves
from treesapp.refpkg import ReferencePackage


_named_vec = namedtuple("_named_vec", ["refpkg", "tp", "tax_rank", "evalue", "hmm_cov",
                                       "leaves", "lwr", "distal", "pendant", "avg_tip_dist"])


class QuerySequence(phylo_seq.PQuery):
    def __init__(self, header: str, lineage_str="", rank_str=""):
        """
        Initialization function for QuerySequence objects. Sets self.ref_name attribute to header if provided

        :param header: Name of the classified sequence, used to set self.ref_name
        """
        super(QuerySequence, self).__init__(lineage_str, rank_str)
        self.place_name = header
        self.ncbi_tax = ""
        self.assigned_lineage = ""
        self.true_lineage = ""
        self.optimal_lineage = ""
        self.tax_dist = 0
        return


class PhyTrainer(classy.TreeSAPP):
    def __init__(self):
        super(PhyTrainer, self).__init__("train")
        self.hmm_purified_seqs = ""  # If an HMM profile of the gene is provided its a path to FASTA with homologs
        self.placement_table = ""
        self.placement_summary = ""
        self.clade_ex_pquery_pkl = ""
        self.plain_pquery_pkl = ""
        self.feature_vector_file = ""
        self.conditions_file = ""
        self.tsne_plot = ""
        self.pkg_dbname_dict = dict()
        self.target_refpkgs = list()
        self.training_ranks = {}
        self.pqueries = {}

        # Stage names only holds the required stages; auxiliary stages (e.g. RPKM, update) are added elsewhere
        self.stages = {0: classy.ModuleFunction("clean", 0),
                       1: classy.ModuleFunction("search", 1),
                       2: classy.ModuleFunction("lineages", 2),
                       3: classy.ModuleFunction("place", 3),
                       4: classy.ModuleFunction("train", 4),
                       5: classy.ModuleFunction("update", 5)}

    def decide_stage(self, args):
        if not args.profile:
            self.change_stage_status("search", False)

        if args.acc_to_lin:
            self.acc_to_lin = args.acc_to_lin
            if os.path.isfile(self.acc_to_lin):
                self.change_stage_status("lineages", False)
            else:
                logging.error("Unable to find accession-lineage mapping file '{}'\n".format(self.acc_to_lin))
                sys.exit(3)

        if os.path.isfile(self.clade_ex_pquery_pkl) and os.path.isfile(self.plain_pquery_pkl):
            self.change_stage_status("place", False)

        self.validate_continue(args)
        return

    def check_trainer_arguments(self, args):
        self.find_sequence_molecule_type()
        self.ref_pkg.f__pkl = args.pkg_path
        self.ref_pkg.slurp()
        self.ref_pkg.validate()

        # Make the directory for storing intermediate outputs
        if not os.path.isdir(self.var_output_dir):
            os.makedirs(self.var_output_dir)

        for rank in args.taxon_rank:
            self.training_ranks[rank] = self.ref_pkg.taxa_trie.accepted_ranks_depths[rank]

        # Check whether the parameters for the classifier make sense
        if args.classifier == "bin":
            if not args.annot_map:
                logging.error("An annotation mapping file is required when building a binary classifier.\n")
                sys.exit(3)
            else:
                self.annot_map = args.annot_map
        elif args.classifier == "occ" and args.annot_map:
            logging.warning("Annotation mapping file is ignored when building a One-Class Classifier (OCC).\n")

        if args.tsne:
            self.tsne_plot = os.path.join(self.var_output_dir, "train", "tSNE.png")

        return

    def set_file_paths(self) -> None:
        """
        Define the file path locations of treesapp train outputs

        :return: None
        """
        self.placement_table = os.path.join(self.final_output_dir, "placement_info.tsv")
        self.placement_summary = os.path.join(self.final_output_dir, "placement_trainer_results.txt")
        self.clade_ex_pquery_pkl = os.path.join(self.final_output_dir, "clade_exclusion_pqueries.pkl")
        self.plain_pquery_pkl = os.path.join(self.final_output_dir, "raw_refpkg_pqueries.pkl")
        self.formatted_input = os.path.join(self.var_output_dir, "clean", self.ref_pkg.prefix + "_formatted.fa")
        self.hmm_purified_seqs = os.path.join(self.var_output_dir, "search", self.ref_pkg.prefix + "_hmm_purified.fa")
        self.acc_to_lin = os.path.join(self.var_output_dir, "lineages", "accession_id_lineage_map.tsv")
        self.conditions_file = os.path.join(self.var_output_dir, "train", "conditions.npy")
        self.feature_vector_file = os.path.join(self.var_output_dir, "train", "examples.npy")
        return

    def get_info(self):
        info_string = "PhyTrainer instance summary:\n"
        info_string += super(PhyTrainer, self).get_info() + "\n\t"
        info_string += "\n\t".join(["Target reference packages = " + str(self.ref_pkg.prefix),
                                    "Taxonomy map: " + self.ref_pkg.lineage_ids,
                                    "Reference tree: " + self.ref_pkg.tree,
                                    "Reference FASTA: " + self.ref_pkg.msa,
                                    "Lineage map: " + str(self.acc_to_lin),
                                    "Ranks tested: " + ','.join(self.training_ranks.keys())]) + "\n"

        return info_string

def summarize_query_classes(positives: set, query_seq_names: set) -> None:
    logging.info("Enumeration of potential query sequence classes:\n")
    false_pos = set()
    false_neg = set()
    true_pos = set()
    total = positives.union(query_seq_names)
    if not positives:
        positives = query_seq_names

    for seq_name in positives:
        if seq_name in query_seq_names:
            true_pos.add(seq_name)
        else:
            false_neg.add(seq_name)
    for seq_name in query_seq_names.difference(true_pos):
        if seq_name not in false_neg:
            false_pos.add(seq_name)

    logging.info("\tTrue positives:          {:>9}\n"
                 "\tFalse positives:         {:>9}\n"
                 "\tFalse negatives:         {:>9}\n"
                 "\tSequence names provided: {:>9}\n"
                 "\n".format(len(true_pos), len(false_pos), len(false_neg), len(total)))
    return


def bin_headers(assignments: dict, annot_map: dict, entrez_query_dict: dict, refpkgs: dict) -> (dict, dict, dict):
    """
    Function for sorting/binning the classified sequences at T/F positives/negatives based on the annot_map file
    that specifies which queries belong to which orthologous group/protein family and reference package.

    :param assignments: Dictionary mapping the ReferencePackage.prefix to TreeProtein instances
    :param annot_map: Dictionary mapping reference package (gene) name keys to database names values
    :param entrez_query_dict: Dictionary mapping sequence names to their respective EntrezRecord instances
    :param refpkgs: A dictionary mapping reference package prefixes to its ReferencePackage instance
    :return: None
    """
    # False positives: those that do not belong to the annotation matching a reference package name
    # False negatives: those whose annotations match a reference package name and were not classified
    # True positives: those with a matching annotation and reference package name and were classified
    # True negatives: those that were not classified and should not have been
    positive_queries = {}
    tp = {}
    fp = {}
    fn = {}
    logging.info("Assigning test sequences to the four class conditions... ")
    missing = set()
    # Create a dictionary for rapid look-ups of queries that are homologs of sequences in the refpkg
    for qname in annot_map:
        for rname in annot_map[qname]:
            try:
                positive_queries[rname].add(qname)
            except KeyError:
                positive_queries[rname] = {qname}

    for refpkg_name, pqueries in assignments.items():  # type: (str, list)
        try:
            positives = positive_queries[refpkg_name]
        except KeyError:
            refpkg = refpkgs[refpkg_name]  # type: ReferencePackage
            if refpkg.refpkg_code in positive_queries:
                positives = positive_queries[refpkg.refpkg_code]
            else:
                logging.error("Unable to find '{}' in the set of positive queries:\n".format(refpkg_name) +
                              ", ".join([str(n) for n in positive_queries.keys()]) + "\n")
                sys.exit(5)
        true_positives = set()
        tp[refpkg_name] = list()
        fp[refpkg_name] = set()
        fn[refpkg_name] = set()
        for pquery in sorted(pqueries, key=lambda x: x.seq_name):  # type: phylo_seq.PQuery
            # Populate the relevant information for the classified sequence
            tp_inst = QuerySequence(pquery.seq_name)
            tp_inst.ref = refpkg_name
            tp_inst.assigned_lineage = pquery.recommended_lineage
            tp_inst.rank = pquery.rank
            # Bin it if the mapping_dict is present, otherwise classify it as a TP
            if tp_inst.place_name in positives:
                try:
                    e_record = entrez_query_dict[pquery.seq_name]
                except KeyError:
                    missing.add(pquery.seq_name)
                    continue
                tp_inst.ncbi_tax = e_record.ncbi_tax
                if not e_record.lineage:
                    true_positives.add(tp_inst.place_name)  # So these queries are not classified as false negatives
                else:
                    # Add the True Positive to the relevant collections
                    tp[refpkg_name].append(tp_inst)
                    true_positives.add(tp_inst.place_name)
            else:
                fp[refpkg_name].add(tp_inst)

        # Identify the False Negatives using set difference - those that were not classified but should have been
        for seq_name in list(positives.difference(true_positives)):
            qseq = QuerySequence(seq_name)
            try:
                e_record = entrez_query_dict[seq_name]
            except KeyError:
                fn[refpkg_name].add(qseq)  # These are false negatives that were not aligned to by the profile HMM
                continue
            if not e_record.lineage:
                continue
            qseq.ncbi_tax = e_record.ncbi_tax
            fn[refpkg_name].add(qseq)
    logging.info("done.\n")

    if missing:
        logging.warning("Unable to find {}/{} sequence accessions in the Entrez records.\n".format(len(missing),
                                                                                                   len(entrez_query_dict)))
        logging.debug("Unable to find the following sequence accessions in the Entrez records:\n"
                      "{}\n".format(', '.join(missing)))

    return tp, fp, fn


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
    ce_tree = ce_refpkg.taxonomically_label_tree()

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
    jplace_data.pqueries = demultiplex_pqueries(jplace_data=jplace_data)
    calc_pquery_mean_tip_distances(jplace_data, internal_node_leaf_map=node_map)
    for pquery in jplace_data.pqueries:  # type: phylo_seq.PQuery
        pquery.ref_name = ref_pkg.prefix
        pquery.rank = rank
        pquery.lineage = taxon
        pquery.process_max_weight_placement(ce_tree)

        if pquery.consensus_placement.like_weight_ratio >= 0.5:
            pqueries.append(pquery)

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
    if not n_reps:
        n_reps = 1
    for _i in range(n_reps):
        # create vector of random gaussians
        gauss = np.random.normal(loc=0.0, scale=feature_scale, size=len(row))
        # add to test case
        new_row = row + gauss
        # store in test set
        training_set = np.append(training_set, [new_row.clip(min=0)])
    return training_set.reshape(n_reps, n_features)


def pquery_to_vector(pquery: phylo_seq.PQuery, ref_pkg: ReferencePackage,
                     positives=None, internal_nodes=None) -> namedtuple:
    if not internal_nodes:
        internal_nodes = ref_pkg.get_internal_node_leaf_map()
    if not positives:
        tp = True
    else:
        tp = pquery.seq_name in positives[ref_pkg.prefix]
    pquery_vec = _named_vec(refpkg=ref_pkg.prefix, tp=tp,
                            tax_rank=pquery.rank, evalue=pquery.evalue,
                            hmm_cov=round((int(pquery.seq_len)*100)/ref_pkg.profile_length, 1),
                            leaves=len(internal_nodes[int(pquery.consensus_placement.edge_num)]),
                            lwr=round(float(pquery.consensus_placement.like_weight_ratio), 3),
                            distal=pquery.consensus_placement.distal_length,
                            pendant=pquery.consensus_placement.pendant_length,
                            avg_tip_dist=pquery.consensus_placement.mean_tip_length)
    return pquery_vec


def load_training_data_frame(pqueries: dict, refpkg_map: dict, refpkg_positive_annots: dict) -> pd.DataFrame:
    """
    Uses attributes from PQuery and ReferencePackage instances to generate data that can be used for training a
    classifier, such as sklearn.svm().

    :param pqueries: A dictionary of ReferencePackage.prefix keys indexing a list of all PQuery instances
    classified as that reference package.
    :param refpkg_map: A dictionary of ReferencePackage.prefix keys indexing that reference package.
    :param refpkg_positive_annots: A dictionary listing the PQuery.seq_name strings indexed by ReferencePackage.prefix
    :return: A pandas.DataFrame with the relevant data for training a classifier.
    """
    training_dat = dict(refpkg=[], tp=[], tax_rank=[],
                        evalue=[], hmm_cov=[],
                        leaves=[], lwr=[], distal=[], pendant=[], avg_tip_dist=[])

    if len(refpkg_positive_annots) == 0:
        return pd.DataFrame(training_dat)

    for refpkg_name, pqueries in pqueries.items():
        ref_pkg = refpkg_map[refpkg_name]  # type: ReferencePackage
        if not ref_pkg.profile_length:
            ref_pkg.hmm_length()
        internal_nodes = ref_pkg.get_internal_node_leaf_map()

        for pquery in pqueries:  # type: phylo_seq.PQuery
            if not pquery.classified:
                continue

            if int(pquery.consensus_placement.edge_num) not in internal_nodes:
                logging.error("Unable to find internal node '{}' in the {} node-leaf map indicating a discrepancy "
                              "between reference package versions used by treesapp assign and those used here.\n"
                              "Was the correct output directory provided?"
                              "".format(pquery.consensus_placement.edge_num, pquery.ref_name))
                sys.exit(5)

            classification_vec = pquery_to_vector(pquery, ref_pkg, refpkg_positive_annots, internal_nodes)

            training_dat["refpkg"].append(classification_vec.refpkg)
            training_dat["tp"].append(classification_vec.tp)
            training_dat["tax_rank"].append(classification_vec.tax_rank)
            training_dat["evalue"].append(classification_vec.evalue)
            training_dat["hmm_cov"].append(classification_vec.hmm_cov)
            training_dat["leaves"].append(classification_vec.leaves)
            training_dat["lwr"].append(classification_vec.lwr)
            training_dat["distal"].append(classification_vec.distal)
            training_dat["pendant"].append(classification_vec.pendant)
            training_dat["avg_tip_dist"].append(classification_vec.avg_tip_dist)

    training_df = pd.DataFrame(training_dat)
    return training_df


def split_placement_data_to_class_vectors(placement_training_df: pd.DataFrame) -> (np.array, np.array):
    """
    Creates a numpy array for true positives and false positives from a DataFrame.

    :param placement_training_df: A pandas.DataFrame with the relevant data for training a classifier.
    :return: Two numpy.array() instances holding data for true positives and false positives
    """
    tp_feature_vectors = np.array([])
    fp_feature_vectors = np.array([])

    if len(placement_training_df) == 0:
        return np.array(tp_feature_vectors), np.array(fp_feature_vectors)

    tp_feature_vectors = vectorize_placement_data_by_rank(placement_training_df[placement_training_df["tp"] == True])
    fp_feature_vectors = vectorize_placement_data_by_rank(placement_training_df[placement_training_df["tp"] == False])

    return tp_feature_vectors, fp_feature_vectors


def vectorize_placement_data_by_rank(placement_training_df: pd.DataFrame) -> np.array:
    """
    Creates an array of numpy array for placement data and indexed by taxonomic rank.

    :param placement_training_df: A pandas.DataFrame with the relevant data for training a classifier.
    :return: Two numpy.array() instances holding data for true positives and false positives
    """
    rank_feature_vectors = {}
    feature_vectors = np.array([])
    columns = ["evalue", "hmm_cov", "leaves", "lwr", "distal", "pendant", "avg_tip_dist"]

    if len(placement_training_df) == 0:
        return np.array(feature_vectors)

    # Generate feature vectors for each rank and index them in the rank_feature_vectors dictionary
    for rank in placement_training_df.tax_rank.unique():
        if not rank:
            continue
        rank_df = placement_training_df[placement_training_df["tax_rank"] == rank].filter(columns)
        rank_feature_vectors[rank] = rank_df.to_numpy()

    # Match the normalized feature vectors with their respective PQuery instance and set PQuery.feature_vec
    if len(rank_feature_vectors) == 0:
        return np.array(feature_vectors)

    # Augment the training data based on the number of examples available at each rank
    limit = max(len(rank_feature_vectors[r]) for r in rank_feature_vectors)
    for rank, examples in rank_feature_vectors.items():
        diff = limit - len(examples)
        for row in examples:
            if len(examples) < limit:
                synthetic_samples = augment_training_set(row, n_reps=round(diff/len(rank_feature_vectors[rank])))
                examples = np.append(examples, synthetic_samples, axis=0)
            else:
                break
        # TODO: replace this, seems hacky
        try:
            feature_vectors = np.append(feature_vectors, examples, axis=0)
        except ValueError:
            feature_vectors = np.array(examples)

    return feature_vectors


def generate_tsne(x, y, tsne_file):
    matplotlib.use('Agg')
    feat_cols = ['pixel' + str(i) for i in range(x.shape[1])]
    df = pd.DataFrame(x, columns=feat_cols)
    df['y'] = y
    df['label'] = df['y'].apply(lambda i: str(i))

    time_start = time.time()
    tsne = manifold.TSNE(n_components=2, verbose=1, perplexity=30, n_iter=300)
    tsne_results = tsne.fit_transform(df)
    logging.debug('t-SNE done! Time elapsed: {} seconds\n'.format(time.time() - time_start))

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
    logging.info("tSNE plot '{}' created.\n".format(tsne_file))
    return


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

    logging.debug("Using a '{}' kernel for the SVM classifier\n".format(k_name))

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
    logging.info("Evaluating model predictions with test data... ")
    # Predict the response for test dataset
    y_pred = clf.predict(x_test)
    scores = model_selection.cross_val_score(clf, classified_data, conditions, cv=10, scoring='f1')
    logging.info("done.\n")

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
    for pq in pqueries:  # type: phylo_seq.PQuery
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
    :param pqueries: A list of PQuery instances
    :return:
    """
    index_pquery_map = {}
    i = 0
    while i < len(feature_vectors):
        vec = feature_vectors[i]  # type: np.ndarray
        index_pquery_map[i] = []
        for pquery in pqueries:  # type: phylo_seq.PQuery
            if np.array_equal(vec, pquery.feature_vec):
                index_pquery_map[i].append(pquery)
        i += 1
    return index_pquery_map


def save_np_array(a: np.array, file_name: str):
    np.save(file_name, a)
    return


def generate_train_test_data(true_ps: np.array, false_ps: np.array, test_pr=0.4):
    """
    Used for splitting the total dataset into training ana testing arrays and mapping indices
    in the positive test fraction to PQueries. This map is to be later used for characterizing predictions by
    rank and lineage.

    Conditions [0|1] are set for the dataset based on whether they are false positives (0) or true positives (1).
    If no false positives are present, a numpy.array of length true_ps is used.

    :param true_ps: A numpy.array instance holding data for true positives
    :param false_ps: A numpy.array instance holding data for false positives
    :param test_pr: Proportion of the training data to use for testing (default = 0.4)
    :return:
    """
    if len(false_ps) > 0:
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

    return classified_data.astype(np.float64), conditions.astype(np.float64),\
           x_train.astype(np.float64), x_test.astype(np.float64), y_train.astype(np.float64), y_test.astype(np.float64)


def train_binary_classifier(x_train: np.array, y_train: np.array, kernel: str):
    # Create a SVM Classifier
    clf = instantiate_classifier(kernel)  # type: svm.SVC

    # logging.info("Using a '{}' kernel for the SVM classifier\n".format(k_name))
    logging.info("Training the '{}' classifier... ".format(kernel))
    # Train the model using the training sets
    clf.fit(x_train, y_train)
    logging.info("done.\n")

    return clf


def train_oc_classifier(x_train: np.array, kernel: str):
    clf = instantiate_classifier(kernel, occ=True)  # type: svm.OneClassSVM

    logging.info("Training the '{}' classifier... ".format(kernel))
    clf.fit(x_train)
    logging.info("done.\n")

    return clf


def train_classification_filter(classified_data: np.array, conditions: np.array,
                                x_train: np.array, x_test: np.array,
                                y_train: np.array, y_test: np.array,
                                kernel: str, index_names: set,
                                grid_search=False, num_procs=2) -> dict:
    """
    Trains a sklearn classifier (e.g. Support Vector Machine) using the TreeSAPP assignments and writes
    the trained model to a pickle file.
    The classifier is tested using 10-fold cross-validation and reports precision, recall, accuracy and F1-scores.

    Optionally, a grid search can be performed that will automatically test multiple different classifiers and kernels,
    not just the one specified, and return (classifier isn't written).

    :param index_names: A dictionary of true positive query names indexed by refpkg names
    :param classified_data: A numpy array of feature vectors for the true positive examples
    :param conditions:
    :param x_train:
    :param x_test:
    :param y_train:
    :param y_test:
    :param kernel: Specifies the kernel type to be used in the SVM algorithm. Choices are 'lin' 'poly' or 'rbf'.
    [ DEFAULT = lin ]
    :param grid_search: Flag indicating whether to perform a grid search across available classifier kernels
    :param num_procs: The number of threads to run the grid search
    :return: None
    """
    classifiers = dict()
    if len(x_train) == 0:
        logging.error("No positive examples are available for training.\n")
        sys.exit(17)
    elif 0 in y_train:  # See if there are any negative examples for binary classification
        if grid_search:
            # Test all the available kernels and return/exit
            kernels = ["lin", "rbf", "poly"]
            for k in kernels:
                evaluate_grid_scores(k, x_train=x_train, x_test=x_test, y_train=y_train, y_test=y_test, jobs=num_procs)
            return classifiers

        # pquery_test_map = map_test_indices_to_pqueries(x_test, true_pqueries)
        clf = train_binary_classifier(x_train, y_train, kernel)

        for refpkg_prefix in index_names:
            classifiers[refpkg_prefix] = clf
            preds = evaluate_classifier(clf, x_test, y_test, classified_data, conditions)
            # characterize_predictions(preds, pquery_test_map)
    else:
        for refpkg_prefix in index_names:
            # pquery_test_map = map_test_indices_to_pqueries(x_test, true_pqueries)

            clf = train_oc_classifier(x_train, kernel)

            classifiers[refpkg_prefix] = clf
            preds = evaluate_classifier(clf, x_test, y_test, classified_data, conditions)
            # characterize_predictions(preds, pquery_test_map)

    return classifiers


def train_classifier_from_dataframe(training_df: pd.DataFrame, kernel: str, num_threads: int,
                                    tsne=False, grid_search=False, trainer=None, refpkg_names=None) -> dict:
    """

    :param training_df:
    :param kernel: Specifies the kernel type to be used in the SVM algorithm. Choices are 'lin' 'poly' or 'rbf'.
    :param num_threads: The number of threads to run the grid search
    :param tsne: Path to a file for writing the tSNE plot. Also used to decide whether to create the tSNE plot
    :param grid_search: Flag indicating whether to perform a grid search across available classifier kernels
    :param trainer: A PhyTrainer instance
    :param refpkg_names: An iterable of reference package names that the classifier will be applied to
    :return: A dictionary of reference package names indexing scikit classifiers
    """
    tp, fp = split_placement_data_to_class_vectors(training_df)
    if refpkg_names is None:
        refpkg_names = set(training_df["refpkg"].unique())
    # Split the training and testing data from the classifications
    classified_data, conditions, x_train, x_test, y_train, y_test = generate_train_test_data(tp, fp)
    # Train the classifier
    classifiers = train_classification_filter(classified_data, conditions,
                                              x_train, x_test, y_train, y_test,
                                              kernel=kernel, index_names=set(refpkg_names),
                                              grid_search=grid_search,
                                              num_procs=num_threads)
    if trainer:
        # Save the feature vectors and condition vectors
        save_np_array(classified_data, trainer.feature_vector_file)
        save_np_array(conditions, trainer.conditions_file)
        if tsne and 0 in y_train:
            generate_tsne(classified_data, conditions, trainer.tsne_plot)
    return classifiers
