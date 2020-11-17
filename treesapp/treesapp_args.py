import argparse
import os
import sys
import re
import logging
from glob import glob
from datetime import datetime as dt

from treesapp.classy import Assigner, Evaluator, Creator, PhyTrainer, Updater, Purity
from treesapp.utilities import available_cpu_count


class TreeSAPPArgumentParser(argparse.ArgumentParser):
    """
    A base argparse ArgumentParser for TreeSAPP with functions to furnish with common arguments.
    This standardizes the interface for a unified aesthetic across all functions: create, classify, evaluate, etc.
    """
    def __init__(self, **kwargs):
        """
        Instantiate the argparse argument-parser and create three broad argument groups:
            reqs - for the required parameters
            optopt - for the optional parameters
            miscellany - for the miscellaneous parameters that are module agnostic,
            for example verbose, help, num_threads
        :param kwargs:
        """
        _prog = " ".join(os.path.basename(x) for x in sys.argv[0:2])
        super(TreeSAPPArgumentParser, self).__init__(add_help=False, prog=_prog, **kwargs)
        self.reqs = self.add_argument_group("Required parameters")
        self.seqops = self.add_argument_group("Sequence operation arguments")
        self.pplace_args = self.add_argument_group("Phylogenetic placement arguments")
        self.rpkm_opts = self.add_argument_group("RPKM options")
        self.io = self.add_argument_group("Inputs and Outputs")
        self.aes = self.add_argument_group("Aesthetic options")
        self.optopt = self.add_argument_group("Optional options")
        self.taxa_args = self.add_argument_group("Taxonomic-lineage arguments")
        self.miscellany = self.add_argument_group("Miscellaneous options")

        self.miscellany.add_argument("-v", "--verbose", action="store_true", default=False,
                                     help="Prints a more verbose runtime log")
        self.miscellany.add_argument("-h", "--help",
                                     action="help",
                                     help="Show this help message and exit")

    def parse_args(self, args=None, namespace=None):
        args = super(TreeSAPPArgumentParser, self).parse_args(args=args, namespace=namespace)

        return args

    # The following are building-block functions
    def add_delete(self):
        self.miscellany.add_argument('--delete', default=False, action="store_true",
                                     help='Delete all intermediate files to save disk space.')

    def add_io(self):
        self.reqs.add_argument('-i', '--fastx_input', required=True, dest="input",
                               help='An input file containing DNA or protein sequences in either FASTA or FASTQ format')
        self.optopt.add_argument('-o', '--output', default='./output/', required=False,
                                 help='Path to an output directory [DEFAULT = ./output/]')
        self.add_delete()
        return

    def add_refpkg_opt(self):
        self.optopt.add_argument("--refpkg_dir", dest="refpkg_dir", default=None,
                                 help="Path to the directory containing reference package pickle (.pkl) files. "
                                      "[ DEFAULT = treesapp/data/ ]")

    def add_refpkg_targets(self):
        self.optopt.add_argument('-t', '--targets', default='', type=str, dest="targets",
                                 help="A comma-separated list specifying which reference packages to use. "
                                      "They are to be referenced by their 'prefix' attribute. "
                                      "Use `treesapp info -v` to get the available list [ DEFAULT = ALL ]")

    def add_annot_map(self, required=False):
        self.optopt.add_argument("--annot_map", required=required, default=None, dest="annot_map",
                                 help="Path to a tabular file mapping reference (refpkg) package names being tested to "
                                      "database corresponding sequence names, indicating a true positive relationship."
                                      " First column is the refpkg name, second is the orthologous group name and third"
                                      " is the query sequence name.")

    def add_refpkg_file_param(self):
        self.reqs.add_argument("-r", "--refpkg_path", dest="pkg_path", required=True, nargs='+',
                               help="Path to the reference package pickle (.pkl) file.\n")

    def add_seq_params(self):
        self.optopt.add_argument("--trim_align", default=False, action="store_true",
                                 help="Flag to turn on position masking of the multiple sequence alignment"
                                      " [DEFAULT = False]")
        self.optopt.add_argument('-w', '--min_seq_length', default=30, type=int,
                                 help='minimal sequence length after alignment trimming [DEFAULT = 30]')
        self.optopt.add_argument('-m', '--molecule', default='dna', choices=['prot', 'dna', 'rrna'],
                                 help="Type of input sequences "
                                      "(prot = protein; dna = nucleotide [DEFAULT]; rrna = rRNA)")

    def add_rpkm_params(self):
        self.rpkm_opts.add_argument("-r", "--reads", required=False,
                                    help="FASTQ file containing to be aligned to predicted genes using BWA MEM")
        self.rpkm_opts.add_argument("-2", "--reverse", required=False,
                                    help="FASTQ file containing to reverse mate-pair reads to be aligned using BWA MEM")
        self.rpkm_opts.add_argument("-p", "--pairing", required=False, default='pe', choices=['pe', 'se'],
                                    help="Indicating whether the reads are paired-end (pe) or single-end (se)")

    def add_search_params(self):
        self.optopt.add_argument("-s", "--stringency", choices=["relaxed", "strict"], default="relaxed", required=False,
                                 help="HMM-threshold mode affects the number of query sequences that advance "
                                      "[DEFAULT = relaxed]")

    def add_phylogeny_params(self):
        self.optopt.add_argument("-b", "--bootstraps",
                                 help="The maximum number of bootstrap replicates RAxML-NG should perform using "
                                      "the autoMRE algorithm.\n[ DEFAULT = 0 ]",
                                 required=False, default=0, type=int)
        self.optopt.add_argument("-e", "--raxml_model",
                                 help="The evolutionary model for RAxML-NG to use\n"
                                      "[ Proteins = LG+G4 | Nucleotides =  GTR+G ]",
                                 required=False, default=None)
        self.optopt.add_argument("--fast",
                                 help="A flag indicating the tree should be built rapidly, using FastTree.",
                                 required=False, default=False, action="store_true")
        # Doesn't _really_ fit in here but good enough. Needs to be used by create and update.
        self.optopt.add_argument("--outdet_align", default=False, action="store_true", dest="od_seq",
                                 help="Flag to activate outlier detection and removal from multiple sequence alignments"
                                      " using OD-seq. [DEFAULT = False]")

    def add_pplace_params(self):
        self.pplace_args.add_argument("-l", "--min_like_weight_ratio", default=0.1, type=float, dest="min_lwr",
                                      help="The minimum likelihood weight ratio required for an EPA placement. "
                                           "[DEFAULT = 0.1]")
        self.pplace_args.add_argument("--placement_summary", default="max_lwr", choices=["aelw", "max_lwr"],
                                      dest="p_sum",
                                      help="Controls the algorithm for consolidating multiple phylogenetic placements."
                                           "Max LWR will take use the phylogenetic placement with greatest LWR."
                                           "aELW uses the taxon with greatest accumulated LWR across placements.")

    def add_compute_miscellany(self):
        self.miscellany.add_argument('--overwrite', action='store_true', default=False,
                                     help='overwrites previously processed output folders')
        self.miscellany.add_argument('-n', '--num_procs', dest="num_threads", default=2, type=int,
                                     help='The number of CPU threads or parallel processes '
                                          'to use in various pipeline steps [DEFAULT = 2]')

    def add_accession_params(self):
        self.taxa_args.add_argument("--accession2taxid", dest="acc_to_taxid", required=False, default=None,
                                    help="Path to an NCBI accession2taxid file "
                                         "for more rapid accession-to-lineage mapping.\n")
        self.taxa_args.add_argument("-a", "--accession2lin", dest="acc_to_lin", required=False, default=None,
                                    help="Path to a file that maps sequence accessions to taxonomic lineages, "
                                         "possibly made by `treesapp create`...")

    def add_lineage_table_param(self):
        self.taxa_args.add_argument("--seqs2lineage", dest="seq_names_to_taxa", required=False, default=None,
                                    help="Path to a file mapping sequence names to taxonomic lineages.\n")

    def add_cluster_args(self):
        self.seqops.add_argument("--cluster",
                                 help="Cluster input sequences at the proportional similarity indicated by identity",
                                 action="store_true",
                                 required=False, default=False)
        self.seqops.add_argument("-p", "--similarity",
                                 help="Proportional similarity (between 0.50 and 1.0) to cluster sequences.",
                                 required=False, default=1.0, type=float)

    def add_taxa_args(self):
        self.taxa_args.add_argument("-s", "--screen",
                                    help="Keywords for including specific taxa in the tree.\n"
                                         "To only include Bacteria and Archaea use `--screen Bacteria,Archaea`\n"
                                         "[ DEFAULT is no screen ]",
                                    default="", required=False)
        self.taxa_args.add_argument("-f", "--filter",
                                    help="Keywords for removing specific taxa; the opposite of `--screen`.\n"
                                         "[ DEFAULT is no filter ]",
                                    default="", required=False)
        self.taxa_args.add_argument("--min_taxonomic_rank",
                                    required=False, default='k', choices=['k', 'p', 'c', 'o', 'f', 'g', 's'],
                                    help="The minimum taxonomic resolution for reference sequences [ DEFAULT = k ].\n")
        self.taxa_args.add_argument("--taxa_lca",
                                    help="Set taxonomy of representative sequences to LCA of cluster member's taxa.\n"
                                         "[ --cluster or --uc REQUIRED ]",
                                    default=False, required=False, action="store_true")
        self.taxa_args.add_argument("--taxa_norm",
                                    help="[ IN DEVELOPMENT ] Subsample leaves by taxonomic lineage.\n"
                                         "A comma-separated argument with the Rank (e.g. Phylum) and\n"
                                         "number of representatives is required.\n")

    def add_taxa_ranks_param(self):
        self.optopt.add_argument("--taxonomic_ranks", dest="taxon_rank", required=False, nargs='+',
                                 default=["class", "species"],
                                 help="A list of the taxonomic ranks (space-separated) to test."
                                      " [ DEFAULT = class species ]",
                                 choices=["domain", "phylum", "class", "order", "family", "genus", "species"])

    def add_classifier_model_params(self):
        self.optopt.add_argument("--max_examples", required=False, default=1E3, type=int,
                                 help="Limits the number of examples used for training models. [ DEFAULT = 1E3 ]")
        self.optopt.add_argument("-k", "--svm_kernel", required=False, default="lin",
                                 choices=["lin", "rbf", "poly"], dest="kernel",
                                 help="Specifies the kernel type to be used in the SVM algorithm. "
                                      "It must be either 'lin' 'poly' or 'rbf'. [ DEFAULT = lin ]")
        self.optopt.add_argument("--grid_search", default=False, required=False, action="store_true",
                                 help="Perform a grid search across hyperparameters. Binary classifier only.")
        self.optopt.add_argument("--tsne", default=False, required=False, action="store_true",
                                 help="Generate a tSNE plot. Output will be in the same directory as the model file. "
                                      "Binary classifier only.")
        self.optopt.add_argument("--classifier", required=False, choices=["occ", "bin"], default="occ",
                                 help="Specify the kind of classifier to be trained: one-class classifier (OCC) or "
                                      "a binary classifier (bin). [ DEFAULT = occ ]")


def add_info_arguments(parser: TreeSAPPArgumentParser):
    parser.add_refpkg_opt()
    return


def add_package_arguments(pkg_parser: TreeSAPPArgumentParser, attributes: list):
    pkg_parser.add_refpkg_file_param()

    pkg_parser.reqs.add_argument("attributes", nargs="+",
                                 help="One or more reference package attributes to view. "
                                      "Note: edit will only modify a single attribute at a time. "
                                      "Choices include: {}\n".format(', '.join(attributes)))
    pkg_parser.optopt.add_argument('-o', '--output', default="./", required=False,
                                   help='Path to an output directory. '
                                        'Default is the current working directory.')
    pkg_parser.optopt.add_argument("--overwrite", default=False, required=False, action="store_true",
                                   help="When editing a reference package, should the current file be overwritten?")
    return


def add_colour_arguments(colour_parser: TreeSAPPArgumentParser) -> None:
    colour_parser.add_refpkg_file_param()

    colour_parser.io.add_argument("-n", "--name", required=False, default=None,
                                  help="The prefix name to use when creating iTOL-compatible output files. "
                                       "By default the taxonomic rank is used or file name of taxa_map if provided.")
    colour_parser.io.add_argument("-o", "--output_dir", dest="output", default="./", required=False,
                                  help="Path to the output directory to write the output files. [ DEFAULT = ./ ]")
    colour_parser.io.add_argument("-t", "--taxa_map", dest="phenotypes", required=False, default=None,
                                  help="A file mapping unique taxonomic labels to non-unique features "
                                       "(e.g. activity, pathway, or other phenotype)")

    colour_parser.aes.add_argument('-l', "--rank_level", dest="rank", default="order", required=False,
                                   help="The rank to generate unique colours for [ DEFAULT = 'order' ]")
    colour_parser.aes.add_argument('-p', "--palette", default="BrBG", required=False,
                                   help="The Seaborn colour palette to use [ DEFAULT = BrBG ]")

    colour_parser.optopt.add_argument('-m', '--min_proportion', dest="min_prop",
                                      default=0.0, required=False, type=float,
                                      help="Minimum proportion of sequences a group contains to assign colour"
                                           " [ DEFAULT = 0 ]")
    colour_parser.optopt.add_argument("-f", "--filter", dest="taxa_filter", default="", required=False,
                                      help="Keywords for excluding specific taxa from the colour palette.\n"
                                           "[ DEFAULT is no filter ]")
    colour_parser.optopt.add_argument("-s", "--taxa_set_operation", dest="set_op", required=False,
                                      choices=['u', 'i'], default='u',
                                      help="When multiple reference packages are provided, should the union (u) or"
                                           " intersection (i) of all labelled taxa (post-filtering) be coloured?"
                                           " [ DEFAULT = 'u' ]")
    colour_parser.optopt.add_argument("--no_polyphyletic", dest="no_poly",
                                      default=False, action="store_true", required=False,
                                      help="Flag forcing the omission of all polyphyletic taxa from the colours file.")
    return


def add_layer_arguments(parser: TreeSAPPArgumentParser):
    parser.add_refpkg_opt()
    parser.reqs.add_argument("-o", "--treesapp_output", dest="output", required=True,
                             help="The TreeSAPP output directory.")
    parser.optopt.add_argument("-c", "--colours_style", required=False, nargs='+',
                               help="The colours_style file exported from iTOL with the annotation information. "
                                     "To automatically infer the variable name (rather than through `names`). "
                                     "File name format should be `marker`_`var`.txt. For example: McrA_Metabolism.txt "
                                     "would create a new column in marker_contig_map.tsv named 'Metabolism'.")
    parser.optopt.add_argument("-d", "--annot_dir", required=False, default=None,
                               help="Path to a directory containing iTOL annotation files for layering.")
    return


def add_classify_arguments(assign_parser: TreeSAPPArgumentParser) -> None:
    """
    Adds command-line arguments that are specific to *treesapp assign*

    :return: None
    """
    assign_parser.add_io()
    assign_parser.add_refpkg_opt()
    assign_parser.add_refpkg_targets()
    assign_parser.add_rpkm_params()
    assign_parser.add_seq_params()
    assign_parser.add_search_params()
    assign_parser.add_pplace_params()
    assign_parser.add_compute_miscellany()
    # The required parameters... for which there are currently none. But they would go here!

    assign_parser.optopt.add_argument("--svm", default=False, required=False, action="store_true",
                                      help="Uses the support vector machine (SVM) classification filter. "
                                           "WARNING: Unless you *really* know your refpkg, you probably don't want this.")

    # The optionals
    assign_parser.optopt.add_argument('-c', '--composition', default="meta", choices=["meta", "single"],
                                      help="Sample composition being either a single organism or a metagenome.")
    assign_parser.optopt.add_argument("--stage", default="continue", required=False,
                                      choices=["continue", "orf-call", "search", "align", "place", "classify"],
                                      help="The stage(s) for TreeSAPP to execute [DEFAULT = continue]")
    assign_parser.rpkm_opts.add_argument("--rpkm", action="store_true", default=False,
                                         help="Flag indicating RPKM values should be calculated for the sequences detected")

    # The miscellany
    assign_parser.miscellany.add_argument('-R', '--reftree', required=False, default="", type=str,
                                          help="[IN PROGRESS] Reference package that all queries should be immediately and "
                                               "directly classified as (i.e. homology search step is skipped).")
    return


def add_abundance_arguments(parser: TreeSAPPArgumentParser):
    parser.add_refpkg_opt()
    parser.add_rpkm_params()
    parser.add_compute_miscellany()
    parser.add_delete()
    parser.reqs.add_argument("--treesapp_output", dest="output", required=True,
                             help="Path to the directory containing TreeSAPP outputs, "
                                  "including sequences to be used for the update.")
    # TODO: Include an option to append new values to the classification table
    parser.optopt.add_argument("--report", choices=["update", "nothing"], required=False, default="update",
                               help="What should be done with the abundance values? The TreeSAPP classification table "
                                    "can be overwritten (update) or left unchanged. "
                                    "[ DEFAULT = update ]")


def add_create_arguments(parser: TreeSAPPArgumentParser) -> None:
    """
    Adds command-line arguments that are specific to *treesapp create*

    :return: None
    """
    parser.add_io()
    parser.add_seq_params()
    parser.add_taxa_args()
    parser.add_cluster_args()
    parser.add_lineage_table_param()
    parser.add_phylogeny_params()
    parser.add_accession_params()
    parser.add_compute_miscellany()

    parser.reqs.add_argument("-c", "--refpkg_name",
                             help="Unique name to be used by TreeSAPP internally. NOTE: Must be <=6 characters.\n"
                                  "Examples are 'McrA', 'DsrAB', and 'p_amoA'.",
                             required=True)

    parser.seqops.add_argument("--multiple_alignment",
                               help='The FASTA input is also the multiple alignment file to be used.\n'
                                    'In this workflow, alignment with MAFFT is skipped and this file is used instead.',
                               action="store_true",
                               default=False)
    parser.seqops.add_argument("-d", "--profile", dest="profile",
                               help="An HMM profile representing a specific domain.\n"
                                    "Domains will be excised from input sequences based on hmmsearch alignments.",
                               required=False, default=None)
    parser.seqops.add_argument("-g", "--guarantee",
                               help="A FASTA file containing sequences that need to be included \n"
                                    "in the tree after all clustering and filtering",
                               default=None,
                               required=False)

    parser.optopt.add_argument("--kind", default="functional", choices=["functional", "taxonomic"], required=False,
                               help="The broad classification of marker gene type, either "
                                    "functional or taxonomic. [ DEFAULT = functional ]")
    parser.optopt.add_argument("--stage", default="continue", required=False,
                               choices=["continue", "search", "lineages", "clean", "cluster", "build",
                                        "evaluate", "support", "train", "update"],
                               help="The stage(s) for TreeSAPP to execute [DEFAULT = continue]")

    parser.miscellany.add_argument('--pc', action='store_true', default=False,
                                   help='Prints the final commands to complete\n'
                                        'installation for a provided `code_name`')
    parser.miscellany.add_argument("--headless", action="store_true", default=False,
                                   help="Do not require any user input during runtime.")


def add_purity_arguments(parser: TreeSAPPArgumentParser) -> None:
    """
    Adds command-line arguments that are specific to *treesapp purity*

    :return: None
    """
    parser.add_io()
    parser.add_refpkg_file_param()
    parser.add_seq_params()
    parser.add_compute_miscellany()
    parser.optopt.add_argument("-x", "--extra_info", required=False, default=None,
                               help="File mapping header prefixes to description information.")
    parser.optopt.add_argument("--stage", default="continue", required=False,
                               choices=["continue", "lineages", "classify", "calculate"],
                               help="The stage(s) for TreeSAPP to execute [DEFAULT = continue]")
    # TODO: Remove --trim_align from command-line options in parser.add_seq_params()
    return


def add_evaluate_arguments(parser: TreeSAPPArgumentParser) -> None:
    """
    Adds command-line arguments that are specific to *treesapp evaluate*

    :return: None
    """
    parser.add_io()
    parser.add_seq_params()
    parser.add_refpkg_file_param()
    parser.add_accession_params()
    parser.add_compute_miscellany()
    parser.add_taxa_ranks_param()

    parser.optopt.add_argument("--fresh", default=False, required=False, action="store_true",
                               help="Recalculate a fresh phylogenetic tree with the target clades removed instead of"
                                    " removing the leaves corresponding to targets from the reference tree.")
    parser.optopt.add_argument("--tool", default="treesapp", required=False,
                               choices=["treesapp", "graftm", "diamond"],
                               help="Classify using one of the tools: treesapp [DEFAULT], graftm, or diamond.")
    parser.optopt.add_argument("-l", "--length",
                               required=False, type=int, default=0,
                               help="Arbitrarily slice the input sequences to this length. "
                                    "Useful for testing classification accuracy for fragments.")
    parser.optopt.add_argument("--stage", default="continue", required=False,
                               choices=["continue", "lineages", "classify", "calculate"],
                               help="The stage(s) for TreeSAPP to execute [DEFAULT = continue]")
    return


def add_update_arguments(parser: TreeSAPPArgumentParser) -> None:
    """
    Adds command-line arguments that are specific to *treesapp update*

    :return: None
    """
    parser.add_io()  # i, o
    parser.add_seq_params()  # w, m
    parser.add_cluster_args()  # p
    parser.add_taxa_args()  # s, f, t
    parser.add_refpkg_file_param()  # r
    parser.add_lineage_table_param()
    parser.add_phylogeny_params()  # b, e
    parser.add_compute_miscellany()  # n
    parser.reqs.add_argument("--treesapp_output", dest="ts_out", required=True,
                             help="Path to the directory containing TreeSAPP outputs, "
                                  "including sequences to be used for the update.")
    parser.optopt.add_argument("-l", "--min_lwr", dest="min_lwr", required=False, default=0.0, type=float,
                               help="The minimum likelihood weight ratio for a sequence to be included in update.")
    parser.optopt.add_argument("--skip_assign", default=False, required=False, action="store_true",
                               help="The assigned sequences are from a database and their database lineages "
                                    "should be used instead of the TreeSAPP-assigned lineages.")
    parser.optopt.add_argument("--resolve", default=False, required=False, action="store_true",
                               help="Flag indicating candidate references with better resolved lineages and"
                                    " comparable sequence lengths can replace old references."
                                    " Useful when updating with sequences from isolates, SAGs and maybe quality MAGs.")
    parser.optopt.add_argument("--stage", default="continue", required=False,
                               choices=["continue", "lineages", "rebuild"],
                               help="The stage(s) for TreeSAPP to execute [DEFAULT = continue]")
    parser.miscellany.add_argument("--headless", action="store_true", default=False,
                                   help="Do not require any user input during runtime.")


def add_trainer_arguments(parser: TreeSAPPArgumentParser) -> None:
    """
    Adds command-line arguments that are specific to *treesapp train*

    :return: None
    """
    parser.add_io()
    parser.add_refpkg_file_param()
    parser.add_seq_params()
    parser.add_accession_params()
    parser.add_taxa_ranks_param()
    parser.add_compute_miscellany()
    parser.add_classifier_model_params()
    parser.add_annot_map()

    parser.seqops.add_argument("-d", "--profile", required=False, default=False, action="store_true",
                               help="Flag indicating input sequences need to be purified using an HMM profile.")
    parser.optopt.add_argument("--stage", default="continue", required=False,
                               choices=["continue", "clean", "search", "lineages", "place", "train", "update"],
                               help="The stage(s) for TreeSAPP to execute [DEFAULT = continue]")
    return


def check_parser_arguments(args, sys_args):
    """
    Function for checking arguments that are found in args.namespace()
    This is the only parser validation function used by clade exclusion evaluator

    :param args: Parsed command-line arguments
    :param sys_args: Unparsed command-line arguments passed to the current TreeSAPP module
    :return: None
    """
    ##
    # Remove the output directory if it exists and overwrite permission granted.
    ##
    logging.info("Arguments used:\n" + ' '.join(sys_args) + "\n")
    if re.match(r"^/$", args.output):
        logging.error("Output directory specified as root. Bailing out to prevent future catastrophe!\n")
        sys.exit(1)
    # Add (or replace a trailing (back)slash with) the os.sep to the end of the output directory
    while re.search(r'.*/$', args.output) or re.search(r'.*\\$', args.output):
        args.output = args.output[:-1]
    args.output += os.sep
    if not re.match(r'^/.*', args.output):
        args.output = os.getcwd() + os.sep + args.output  # args.output is now the absolute path

    if "min_seq_length" not in vars(args):
        args.min_seq_length = 1

    if sys.version_info <= (2, 9):
        logging.error("Python 2 is not supported by TreeSAPP.\n")
        sys.exit(3)

    if "num_threads" in vars(args) and args.num_threads > available_cpu_count():
        logging.warning("Number of threads specified is greater than those available! "
                        "Using maximum threads available (" + str(available_cpu_count()) + ")\n")
        args.num_threads = available_cpu_count()

    return


def check_purity_arguments(purity_instance: Purity, args):
    purity_instance.ref_pkg.f__json = args.pkg_path
    purity_instance.ref_pkg.slurp()
    purity_instance.refpkg_dir = os.path.dirname(purity_instance.ref_pkg.f__json)

    ##
    # Define locations of files TreeSAPP outputs
    ##
    purity_instance.assign_dir = purity_instance.var_output_dir + "assign" + os.sep
    purity_instance.summarize_dir = purity_instance.var_output_dir + "summarize" + os.sep
    purity_instance.classifications = purity_instance.assign_dir + "final_outputs" + os.sep + "marker_contig_map.tsv"
    purity_instance.metadata_file = args.extra_info

    if not os.path.isdir(purity_instance.var_output_dir):
        os.makedirs(purity_instance.var_output_dir)

    return


def check_evaluate_arguments(evaluator_instance: Evaluator, args):
    for rank in args.taxon_rank:
        evaluator_instance.ranks.append(rank)

    evaluator_instance.ref_pkg.f__json = args.pkg_path
    evaluator_instance.ref_pkg.slurp()

    if args.acc_to_lin:
        evaluator_instance.acc_to_lin = args.acc_to_lin
        if os.path.isfile(evaluator_instance.acc_to_lin):
            evaluator_instance.change_stage_status("lineages", False)
        else:
            logging.error("Unable to find accession-lineage mapping file '{}'\n".format(evaluator_instance.acc_to_lin))
            sys.exit(3)
    else:
        evaluator_instance.acc_to_lin = evaluator_instance.var_output_dir + os.sep + "accession_id_lineage_map.tsv"

    ##
    # Define locations of files TreeSAPP outputs
    ##
    evaluator_instance.test_rep_taxa_fasta = evaluator_instance.final_output_dir + "representative_taxa_sequences.fasta"
    evaluator_instance.performance_table = evaluator_instance.final_output_dir + "clade_exclusion_performance.tsv"
    evaluator_instance.recall_table = evaluator_instance.final_output_dir + "taxonomic_recall.tsv"
    evaluator_instance.containment_table = evaluator_instance.final_output_dir + "accuracy.tsv"
    evaluator_instance.var_output_dir = args.output + "intermediates" + os.sep

    if not os.path.isdir(evaluator_instance.var_output_dir):
        os.makedirs(evaluator_instance.var_output_dir)

    return


def check_trainer_arguments(phy_trainer: PhyTrainer, args):
    phy_trainer.ref_pkg.f__json = args.pkg_path
    phy_trainer.ref_pkg.slurp()
    phy_trainer.ref_pkg.validate()

    # Make the directory for storing intermediate outputs
    if not os.path.isdir(phy_trainer.var_output_dir):
        os.makedirs(phy_trainer.var_output_dir)

    for rank in args.taxon_rank:
        phy_trainer.training_ranks[rank] = phy_trainer.ref_pkg.taxa_trie.accepted_ranks_depths[rank]

    # Check whether the parameters for the classifier make sense
    if args.classifier == "bin":
        if not args.annot_map:
            logging.error("An annotation mapping file is required when building a binary classifier.\n")
            sys.exit(3)
        else:
            phy_trainer.annot_map = args.annot_map
    elif args.classifier == "occ" and args.annot_map:
        logging.warning("Annotation mapping file is ignored when building a One-Class Classifier (OCC).\n")

    if args.tsne:
        phy_trainer.tsne_plot = os.path.join(phy_trainer.var_output_dir, "train", "tSNE.png")

    return


def check_classify_arguments(assigner: Assigner, args):
    """
    Ensures the command-line arguments returned by argparse are sensible.

    :param assigner: An instantiated Assigner object
    :param args: object with parameters returned by argparse.parse_args()
    :return: 'args', a summary of TreeSAPP settings.
    """
    assigner.aa_orfs_file = assigner.final_output_dir + assigner.sample_prefix + "_ORFs.faa"
    assigner.nuc_orfs_file = assigner.final_output_dir + assigner.sample_prefix + "_ORFs.fna"
    assigner.classified_aa_seqs = assigner.final_output_dir + assigner.sample_prefix + "_classified.faa"
    assigner.classified_nuc_seqs = assigner.final_output_dir + assigner.sample_prefix + "_classified.fna"
    
    if args.targets:
        assigner.target_refpkgs = args.targets.split(',')
    else:
        assigner.target_refpkgs = []

    assigner.validate_refpkg_dir(args.refpkg_dir)

    if args.molecule == "prot":
        assigner.change_stage_status("orf-call", False)
        if args.rpkm:
            logging.error("Unable to calculate RPKM values for protein sequences.\n")
            sys.exit(3)

    if args.svm:
        assigner.svc_filter = True

    # TODO: transfer all of this HMM-parsing stuff to the assigner_instance
    # Parameterizing the hmmsearch output parsing:
    args.perc_aligned = 10
    args.min_acc = 0.7
    if args.stringency == "relaxed":
        args.max_e = 1E-3
        args.max_ie = 1E-1
        args.min_score = 15
    elif args.stringency == "strict":
        args.max_e = 1E-5
        args.max_ie = 1E-3
        args.min_score = 30
    else:
        logging.error("Unknown HMM-parsing stringency argument '" + args.stringency + "'.\n")
        sys.exit(3)

    return args


def check_create_arguments(creator: Creator, args) -> None:
    # Populate ReferencePackage attributes from command-line arguments
    if args.fast:
        creator.ref_pkg.tree_tool = "FastTree"
    else:
        creator.ref_pkg.tree_tool = "RAxML-NG"
    creator.ref_pkg.prefix = args.refpkg_name
    creator.ref_pkg.pid = args.similarity
    creator.ref_pkg.molecule = args.molecule
    creator.ref_pkg.kind = args.kind
    creator.ref_pkg.sub_model = args.raxml_model
    creator.ref_pkg.date = dt.now().strftime("%Y-%m-%d")
    creator.ref_pkg.f__json = creator.final_output_dir + creator.ref_pkg.prefix + creator.ref_pkg.refpkg_suffix
    # TODO: Create placement trainer output directory and make it an attribute
    if not args.output:
        args.output = os.getcwd() + os.sep + creator.ref_pkg.prefix + "_treesapp_refpkg" + os.sep

    if len(creator.ref_pkg.prefix) > 10:
        logging.error("Name should be <= 10 characters.\n")
        sys.exit(13)

    # TODO: Check the substitution model for compatibility with RAxML-NG

    if args.cluster:
        if args.multiple_alignment:
            logging.error("--cluster and --multiple_alignment are mutually exclusive!\n")
            sys.exit(13)
        if not 0.5 <= float(args.similarity) <= 1.0:
            if 0.5 < float(args.similarity)/100 < 1.0:
                args.similarity = str(float(args.similarity)/100)
                logging.warning("--similarity  set to {} for compatibility with VSEARCH.\n".format(args.similarity))
            else:
                logging.error("--similarity {} is not between the supported range [0.5-1.0].\n".format(args.similarity))
                sys.exit(13)

    if args.taxa_lca and not args.cluster:
        logging.error("Unable to perform LCA for representatives without clustering information: " +
                      "either with a provided VSEARCH file or by clustering within the pipeline.\n")
        sys.exit(13)

    if args.guarantee and not args.cluster:
        logging.error("--guarantee used but without clustering there is no reason for it.\n" +
                      "Include all sequences in " + args.guarantee +
                      " in " + creator.input_sequences + " and re-run without --guarantee\n")
        sys.exit(13)

    if args.profile:
        if not os.path.isfile(args.profile):
            logging.error("Unable to find HMM profile at '" + args.profile + "'.\n")
            sys.exit(3)
        creator.hmm_profile = args.profile

    # Names of files and directories to be created
    creator.phy_dir = os.path.abspath(creator.var_output_dir) + os.sep + "phylogeny_files" + os.sep
    creator.training_dir = os.path.abspath(creator.var_output_dir) + os.sep + "placement_trainer" + os.sep
    creator.hmm_purified_seqs = creator.var_output_dir + creator.ref_pkg.prefix + "_hmm_purified.fasta"
    creator.filtered_fasta = creator.var_output_dir + creator.sample_prefix + "_filtered.fa"
    creator.cluster_input = creator.var_output_dir + creator.sample_prefix + "_cluster_input.fasta"
    creator.clusters_prefix = creator.var_output_dir + creator.sample_prefix + "_cluster" + str(creator.ref_pkg.pid)
    creator.unaln_ref_fasta = creator.var_output_dir + creator.ref_pkg.prefix + "_ref.fa"
    creator.phylip_file = creator.var_output_dir + creator.ref_pkg.prefix + ".phy"

    # Ensure the phylogenetic tree output directory from a previous run isn't going to be over-written
    if not os.path.exists(creator.phy_dir):
        os.mkdir(creator.phy_dir)
    else:
        logging.error(creator.phy_dir + " already exists from a previous run! " +
                      "Please delete or rename it and try again.\n")
        sys.exit(13)

    if not os.path.isdir(creator.training_dir):
        os.mkdir(creator.training_dir)

    return


def check_updater_arguments(updater: Updater, args):
    updater.ref_pkg.f__json = args.pkg_path
    updater.ref_pkg.slurp()
    updater.updated_refpkg_path = os.path.join(updater.output_dir, "final_outputs",
                                               os.path.basename(updater.ref_pkg.f__json))
    updater.ref_pkg.disband(os.path.join(updater.output_dir, "intermediates"))
    updater.seq_names_to_taxa = args.seq_names_to_taxa
    # updater.rank_depth_map = {'k': 1, 'p': 2, 'c': 3, 'o': 4, 'f': 5, 'g': 6, 's': 7}

    if args.similarity == 1.0:
        updater.prop_sim = updater.ref_pkg.pid
    else:
        updater.prop_sim = args.similarity

    if args.cluster:
        if not 0.5 <= float(args.similarity) <= 1.0:
            if 0.5 < float(args.similarity) / 100 < 1.0:
                args.similarity = str(float(args.similarity) / 100)
                logging.warning("--similarity set to {} for compatibility.\n".format(args.similarity))
            else:
                logging.error("--similarity {} is not between the supported range [0.5-1.0].\n".format(args.similarity))
                sys.exit(13)

    if updater.seq_names_to_taxa and not os.path.isfile(updater.seq_names_to_taxa):
        logging.error("Unable to find file mapping sequence names to taxonomic lineages '" +
                      updater.seq_names_to_taxa + "'.\n")

    # TODO: Write a TreeSAPP function for validating outputs

    # Reset these values to reflect the paths for the TreeSAPP output that will be parsed
    updater.treesapp_output = args.ts_out
    if updater.treesapp_output[-1] != os.sep:
        updater.treesapp_output += os.sep
    updater.final_output_dir = updater.treesapp_output + "final_outputs" + os.sep
    # updater.var_output_dir = updater.treesapp_output + "intermediates" + os.sep
    updater.old_ref_fasta = updater.var_output_dir + "original_refs.fasta"
    updater.combined_fasta = updater.var_output_dir + "all_refs.fasta"
    updater.lineage_map_file = updater.var_output_dir + "accession_id_lineage_map.tsv"
    updater.assignment_table = updater.final_output_dir + "marker_contig_map.tsv"
    updater.cluster_input = updater.var_output_dir + updater.sample_prefix + "_cluster_input.fasta"
    updater.clusters_prefix = updater.var_output_dir + updater.sample_prefix + "_cluster" + str(updater.prop_sim)
    classified_seqs = glob(updater.final_output_dir + "*_classified.faa")

    if len(classified_seqs) == 1:
        updater.query_sequences = classified_seqs.pop()
    elif len(classified_seqs) == 0:
        logging.error("No classified sequence files found in {}.\n".format(updater.final_output_dir))
    else:
        logging.error("Multiple classified sequence files in '{}'"
                      " where only one expected.\n".format(updater.final_output_dir))
        sys.exit(5)

    return
