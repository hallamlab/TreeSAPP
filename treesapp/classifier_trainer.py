#!/usr/bin/env python3

__author__ = 'Connor Morgan-Lang'

import os
import sys
import logging
from glob import glob

from sklearn import svm
import joblib

from treesapp.classy import prep_logging, PhyTrainer
from treesapp.refpkg import ReferencePackage
from treesapp.fasta import FASTA, write_new_fasta
from treesapp.commands import assign, evaluate
from treesapp.phylo_seq import assignments_to_treesaps
from treesapp.training_utils import train_classification_filter
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
    parser.add_classifier_model_params()
    treesapp_args.add_trainer_arguments(parser)

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
        clade_exclusion_outputs = glob(os.path.join(clade_exclusion_prefix, "intermediates", '*',
                                                    "treesapp_output", "final_outputs", "marker_contig_map.tsv"))
        # Read the classification tables to gather the training data
        classification_tables = clade_exclusion_outputs
        classification_tables.append(assign_table)
        for table in classification_tables:
            classification_lines += ts_fp.read_classification_table(table)

    assignments = assignments_to_treesaps(classification_lines, refpkg_dict)

    return assignments


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
        pkg_dbname_dict = ts_fp.read_annotation_mapping_file(args.annot_map)
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
    test_obj.bin_headers(assignments, pkg_dbname_dict)
    test_seq_names.clear()

    ##
    # Bin the test sequence names into their respective confusion categories (TP, TN, FP, FN)
    ##
    test_obj.populate_tax_lineage_map(entrez_record_dict)
    test_obj.map_true_lineages()
    test_obj.validate_false_positives()
    if args.annot_map:
        test_obj.validate_false_negatives(args.annot_map)

    summarize_training_rank_coverage(test_obj)

    # Convert the true positive dictionary to the same format, flattening the QuerySequence instances
    flattened_tp = dict()
    for refpkg in test_obj.tp:
        flattened_tp[refpkg] = set()
        for classified_seq in test_obj.tp[refpkg]:  # type: ts_MCC.QuerySequence
            flattened_tp[refpkg].add(classified_seq.name)

    classifiers = train_classification_filter(assignments, flattened_tp, test_obj.ref_packages,
                                              args.kernel, tsne_file, args.grid_search, test_obj.fp, args.num_threads)

    for refpkg_name in classifiers:
        try:
            refpkg = refpkg_dict[refpkg_name]  # type: ReferencePackage
        except KeyError:
            logging.warning("A classifier wasn't trained for the ReferencePackage '{}'"
                            " likely because no true positives were identified.\n".format(refpkg_name))
            continue
        refpkg.svc = classifiers[refpkg_name]
        refpkg.f__json = ts_trainer.final_output_dir + os.path.basename(refpkg.f__json)
        refpkg.pickle_package()

    return


if __name__ == '__main__':
    classifier_trainer(sys.argv)
