import argparse
import os
import sys
import re
import logging
from glob import glob

from treesapp.classy import Evaluator
from treesapp.update_refpkg import Updater
from treesapp.utilities import available_cpu_count
from treesapp import logger

LOGGER = logging.getLogger(logger.logger_name())


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
        self.hmmer_args = self.add_argument_group("Homology search arguments")
        self.pplace_args = self.add_argument_group("Phylogenetic placement arguments")
        self.fpkm_opts = self.add_argument_group("Abundance options")
        self.io = self.add_argument_group("Inputs and Outputs")
        self.aes = self.add_argument_group("Aesthetic options")
        self.editors = self.add_argument_group("Package edit options")
        self.optopt = self.add_argument_group("Optional arguments")
        self.taxa_args = self.add_argument_group("Taxonomic-lineage arguments")
        self.svc_opts = self.add_argument_group("Classifier arguments")
        self.miscellany = self.add_argument_group("Miscellaneous options")

        self.miscellany.add_argument('--overwrite', action='store_true', default=False,
                                     help='Overwrites previously written output files and directories')
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

    def add_input_fastx(self):
        self.reqs.add_argument('-i', '--fastx_input', required=True, dest="input", nargs='+',
                               help='An input file containing DNA or protein sequences in either FASTA or FASTQ format')

    def add_output_dir(self):
        self.optopt.add_argument('-o', '--output', default='./output/', required=False,
                                 help='Path to an output directory [DEFAULT = ./output/]')

    def add_io(self):
        self.add_input_fastx()
        self.add_output_dir()
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
        self.optopt.add_argument('-w', '--min_seq_length', default=0, type=int,
                                 help='minimal sequence length after alignment trimming [DEFAULT = 0]')
        self.optopt.add_argument('-m', '--molecule', choices=['prot', 'dna', 'rrna'], required=False, default="",
                                 help="Type of input sequences (prot = protein; dna = nucleotide; rrna = rRNA). "
                                      "TreeSAPP will guess by default but this may be required if ambiguous.")

    def add_search_params(self):
        self.hmmer_args.add_argument("-s", "--stringency",
                                     choices=["relaxed", "strict"], default="relaxed", required=False,
                                     help="HMM-threshold mode affects the number of query sequences that advance "
                                          "[DEFAULT = relaxed]")
        self.hmmer_args.add_argument("-P", "--hmm_coverage", type=int, default=80, required=False,
                                     help="Minimum percentage of a profile HMM that a query alignment must cover "
                                          "for it to be considered. [ DEFAULT = 80 ]")
        self.hmmer_args.add_argument("-Q", "--query_coverage", type=int, default=80, required=False,
                                     help="Minimum percentage of a query sequence that an alignment must cover "
                                          "to be retained. [ DEFAULT = 80 ]")

    def add_abundance_params(self):
        self.fpkm_opts.add_argument("--metric", required=False, default="tpm", choices=["fpkm", "tpm"],
                                    help="Selects which normalization metric to use, FPKM or TPM. [ DEFAULT = tpm ]")
        self.fpkm_opts.add_argument("-r", "--reads", required=False, nargs='+',
                                    help="FASTQ file containing to be aligned to predicted genes using BWA MEM")
        self.fpkm_opts.add_argument("-2", "--reverse", required=False, nargs='+', default=[],
                                    help="FASTQ file containing to reverse mate-pair reads to be aligned using BWA MEM")
        self.fpkm_opts.add_argument("-p", "--pairing", required=False, default='pe', choices=['pe', 'se'],
                                    help="Indicating whether the reads are paired-end (pe) or single-end (se). "
                                         "[ DEFAULT = pe ]")

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
                                      " using OD-seq. [ DEFAULT = False ]")

    def add_pplace_filter_params(self):
        self.pplace_args.add_argument("--min_like_weight_ratio", default=0.1, type=float, dest="min_lwr",
                                      help="The minimum likelihood weight ratio required for an EPA placement. "
                                           "[ DEFAULT = 0.1 ]")
        self.pplace_args.add_argument("--max_pendant_length", default=None, type=float, dest="max_pd",
                                      help="The maximum pendant length distance threshold, "
                                           "beyond which EPA placements are unclassified. [ DEFAULT = Inf ]")
        self.pplace_args.add_argument("--max_evol_distance", default=None, type=float, dest="max_evo",
                                      help="The maximum total evolutionary distance between a query and reference(s), "
                                           "beyond which EPA placements are unclassified. [ DEFAULT = Inf ]")

    def add_pplace_params(self):
        self.add_pplace_filter_params()
        self.pplace_args.add_argument("--placement_summary", default="max_lwr", choices=["aelw", "max_lwr", "lca"],
                                      dest="p_sum",
                                      help="Controls the algorithm for consolidating multiple phylogenetic placements. "
                                           "Max LWR will take use the phylogenetic placement with greatest LWR. "
                                           "aELW uses the taxon with greatest accumulated LWR across placements.")

    def add_compute_miscellany(self):
        self.miscellany.add_argument('-n', '--num_procs', dest="num_threads", default=2, type=int,
                                     help='The number of CPU threads or parallel processes '
                                          'to use in various pipeline steps [ DEFAULT = 2 ]')

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
                                    required=False, default='r', choices=['r', 'd', 'p', 'c', 'o', 'f', 'g', 's'],
                                    help="The minimum taxonomic resolution for reference sequences [ DEFAULT = r ].\n")
        self.taxa_args.add_argument("--taxa_lca",
                                    help="Set taxonomy of representative sequences to LCA of cluster member's taxa.\n"
                                         "[ --cluster REQUIRED ]",
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

    def add_classifier_kernel_param(self):
        self.svc_opts.add_argument("-k", "--svm_kernel", required=False, default="lin",
                                   choices=["lin", "rbf", "poly"], dest="kernel",
                                   help="Specifies the kernel type to be used in the SVM algorithm. "
                                        "It must be either 'lin' 'poly' or 'rbf'. [ DEFAULT = lin ]")
        return

    def add_basic_classifier_model_params(self):
        self.add_classifier_kernel_param()
        self.svc_opts.add_argument("--max_examples", required=False, default=1000, type=int,
                                   help="Limits the number of examples used for training models. [ DEFAULT = 1000 ]")
        return

    def add_advanced_classifier_params(self):
        self.svc_opts.add_argument("--grid_search", default=False, required=False, action="store_true",
                                   help="Perform a grid search across hyperparameters. Binary classifier only.")
        self.svc_opts.add_argument("--tsne", default=False, required=False, action="store_true",
                                   help="Generate a tSNE plot. Output will be in the same directory as the model file. "
                                        "Binary classifier only.")
        self.svc_opts.add_argument("--classifier", required=False, choices=["occ", "bin"], default="occ",
                                   help="Specify the kind of classifier to be trained: one-class classifier (OCC) or "
                                        "a binary classifier (bin). [ DEFAULT = occ ]")
        return


def add_info_arguments(parser: TreeSAPPArgumentParser):
    parser.add_refpkg_opt()
    return


def add_package_arguments(pkg_parser: TreeSAPPArgumentParser, attributes: list):
    pkg_parser.add_refpkg_file_param()

    pkg_parser.reqs.add_argument("attributes", nargs="+",
                                 help="One or more reference package attributes to view or edit. "
                                      "Note: edit will only modify a single attribute at a time. "
                                      "Choices include: {}\n".format(', '.join(attributes)))
    pkg_parser.editors.add_argument("-t", "--taxa_map", dest="phenotypes", required=False, default=None,
                                    help="A file mapping unique taxonomic labels to non-unique features "
                                         "(e.g. activity, pathway, or other phenotype)")
    pkg_parser.editors.add_argument("--reset", default=False, action="store_true",
                                    help="Flag to reset the reference package attribute when editing.")
    pkg_parser.editors.add_argument("--join", default=False, action="store_true",
                                    help="Flag indicating the input needs to be joined to the existing attribute. "
                                         "Applicable to: 'lineage_ids'. [ DEFAULT = False ]")
    pkg_parser.optopt.add_argument('-o', '--output', default=None, required=False,
                                   help='Path to an output directory. '
                                        'Default is the current working directory.')
    return


def add_colour_arguments(colour_parser: TreeSAPPArgumentParser) -> None:
    colour_parser.add_refpkg_file_param()

    colour_parser.io.add_argument("-n", "--attribute", required=False, default="taxonomy",
                                  help="The reference package attribute to colour by. "
                                       "Either 'taxonomy' or a reference package's layering annotation name. "
                                       "[ DEFAULT = 'taxonomy' ]")
    colour_parser.io.add_argument("-o", "--output_dir", dest="output", default="./", required=False,
                                  help="Path to the output directory to write the output files. [ DEFAULT = ./ ]")

    colour_parser.aes.add_argument('-l', "--rank_level", dest="rank", default="order", required=False,
                                   choices=["domain", "phylum", "class", "order", "family", "genus", "species"],
                                   help="The rank to generate unique colours for [ DEFAULT = 'order' ]")
    colour_parser.aes.add_argument('-p', "--palette", default="BrBG", required=False,
                                   help="The Seaborn colour palette to use [ DEFAULT = BrBG ]")
    colour_parser.aes.add_argument("--unknown_colour", default="", required=False, dest="un_col",
                                   help="Colour of the 'Unknown' category. [ DEFAULT = None ]")

    colour_parser.optopt.add_argument('-m', '--min_proportion', dest="min_prop",
                                      default=0.0, required=False, type=float,
                                      help="Minimum proportion of sequences a group contains to be coloured"
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
    return


def add_classify_arguments(assign_parser: TreeSAPPArgumentParser) -> None:
    """
    Adds command-line arguments that are specific to *treesapp assign*

    :return: None
    """
    assign_parser.add_io()
    assign_parser.add_refpkg_opt()
    assign_parser.add_refpkg_targets()
    assign_parser.add_abundance_params()
    assign_parser.add_seq_params()
    assign_parser.add_search_params()
    assign_parser.add_pplace_params()
    assign_parser.add_compute_miscellany()
    assign_parser.add_classifier_kernel_param()
    # The required parameters... for which there are currently none. But they would go here!

    assign_parser.optopt.add_argument("--svm", default=False, required=False, action="store_true",
                                      help="Uses the support vector machine (SVM) classification filter. "
                                           "WARNING: Unless you *really* know your refpkg, you don't want this.")

    # The optionals
    assign_parser.optopt.add_argument('-c', '--composition', default="meta", choices=["meta", "single"],
                                      help="Sample composition being either a single organism or a metagenome.")
    assign_parser.optopt.add_argument("--stage", default="continue", required=False,
                                      choices=["continue", "orf-call", "search", "align", "place", "classify"],
                                      help="The stage(s) for TreeSAPP to execute [DEFAULT = continue]")
    assign_parser.fpkm_opts.add_argument("--rel_abund", action="store_true", default=False,
                                         help="Flag indicating relative abundance values should be"
                                              " calculated for the sequences detected")

    # The miscellany
    assign_parser.miscellany.add_argument('-R', '--reftree', required=False, default="", type=str,
                                          help="[IN PROGRESS] Reference package that all queries should be immediately"
                                               " and directly classified as (i.e. homology search step is skipped).")
    assign_parser.miscellany.add_argument("--silent", action="store_true", default=False,
                                          help="treesapp assign will not log anything to the console if used.")
    return


def add_abundance_arguments(abundance_parser: TreeSAPPArgumentParser):
    abundance_parser.add_abundance_params()
    abundance_parser.add_refpkg_opt()
    abundance_parser.add_compute_miscellany()
    abundance_parser.add_delete()
    abundance_parser.reqs.add_argument("--treesapp_output", dest="output", required=True,
                                       help="Path to the directory containing TreeSAPP outputs, "
                                            "including sequences to be used for the update.")
    abundance_parser.optopt.add_argument("--report", choices=["update", "nothing", "append"],
                                         required=False, default="append",
                                         help="What should be done with the abundance values? "
                                              "The TreeSAPP classification table can be overwritten (update), "
                                              "appended or left unchanged. [ DEFAULT = append ]")
    abundance_parser.optopt.add_argument("--stage", default="continue", required=False,
                                         choices=["continue", "align_map", "sam_sum", "summarise"],
                                         help="The stage(s) for TreeSAPP to execute [DEFAULT = continue]")
    return


def add_create_arguments(crt_parser: TreeSAPPArgumentParser) -> None:
    """
    Adds command-line arguments that are specific to *treesapp create*

    :return: None
    """
    crt_parser.add_io()
    crt_parser.add_seq_params()
    crt_parser.add_taxa_args()
    crt_parser.add_cluster_args()
    crt_parser.add_lineage_table_param()
    crt_parser.add_phylogeny_params()
    crt_parser.add_accession_params()
    crt_parser.add_compute_miscellany()
    crt_parser.add_basic_classifier_model_params()

    crt_parser.reqs.add_argument("-c", "--refpkg_name",
                                 help="Unique name to be used by TreeSAPP internally. NOTE: Must be <=6 characters.\n"
                                      "Examples are 'McrA', 'DsrAB', and 'p_amoA'.",
                                 required=True)

    crt_parser.seqops.add_argument("--deduplicate", default=False, dest="dedup", action="store_true",
                                   help="Deduplicate the input sequences at 99.9 percent similarity. "
                                        "This is a pre-processing step to require fewer Entrez queries - "
                                        "clustering at lower resolution with '--cluster' is still suggested.")
    crt_parser.seqops.add_argument("--multiple_alignment",
                                   help='The FASTA input is also the multiple alignment file to be used.\n',
                                   action="store_true",
                                   default=False)
    crt_parser.seqops.add_argument("-d", "--profile", dest="profile",
                                   help="An HMM profile representing a specific domain.\n"
                                        "Domains will be excised from input sequences based on hmmsearch alignments.",
                                   required=False, default=None)
    crt_parser.seqops.add_argument("-g", "--guarantee",
                                   help="A FASTA file containing sequences that need to be included \n"
                                        "in the tree after all clustering and filtering",
                                   default=None,
                                   required=False)

    crt_parser.optopt.add_argument("--kind", default="functional", choices=["functional", "taxonomic"], required=False,
                                   help="The broad classification of marker gene type, either "
                                        "functional or taxonomic. [ DEFAULT = functional ]")
    crt_parser.optopt.add_argument("--stage", default="continue", required=False,
                                   choices=["continue", "search", "lineages", "clean", "cluster", "build",
                                            "evaluate", "support", "train", "update"],
                                   help="The stage(s) for TreeSAPP to execute [DEFAULT = continue]")

    crt_parser.miscellany.add_argument("--headless", action="store_true", default=False,
                                       help="Do not require any user input during runtime.")
    return


def add_purity_arguments(pur_parser: TreeSAPPArgumentParser) -> None:
    """
    Adds command-line arguments that are specific to *treesapp purity*

    :return: None
    """
    pur_parser.add_io()
    pur_parser.add_refpkg_file_param()
    pur_parser.add_seq_params()
    pur_parser.add_compute_miscellany()
    pur_parser.optopt.add_argument("-x", "--extra_info", required=False, default=None,
                                   help="File mapping header prefixes to description information.")
    pur_parser.optopt.add_argument("--stage", default="continue", required=False,
                                   choices=["continue", "lineages", "classify", "calculate"],
                                   help="The stage(s) for TreeSAPP to execute [DEFAULT = continue]")
    # TODO: Remove --trim_align from command-line options in pur_parser.add_seq_params()
    return


def add_evaluate_arguments(eval_parser: TreeSAPPArgumentParser) -> None:
    """
    Adds command-line arguments that are specific to *treesapp evaluate*

    :return: None
    """
    eval_parser.add_io()
    eval_parser.add_seq_params()
    eval_parser.add_refpkg_file_param()
    eval_parser.add_accession_params()
    eval_parser.add_compute_miscellany()
    eval_parser.add_pplace_params()
    eval_parser.add_taxa_ranks_param()

    eval_parser.optopt.add_argument("--fresh", default=False, required=False, action="store_true",
                                    help="Recalculate a fresh phylogenetic tree with the target clades removed instead"
                                         " of removing the leaves corresponding to targets from the reference tree.")
    eval_parser.optopt.add_argument("--tool", default="treesapp", required=False,
                                    choices=["treesapp", "graftm", "diamond"],
                                    help="Classify using one of the tools: treesapp [DEFAULT], graftm, or diamond.")
    eval_parser.optopt.add_argument("-l", "--length",
                                    required=False, type=int, default=0,
                                    help="Arbitrarily slice the input sequences to this length. "
                                         "Useful for testing classification accuracy for fragments.")
    eval_parser.optopt.add_argument("--stage", default="continue", required=False,
                                    choices=["continue", "lineages", "classify", "calculate"],
                                    help="The stage(s) for TreeSAPP to execute [DEFAULT = continue]")
    return


def add_update_arguments(up_parser: TreeSAPPArgumentParser) -> None:
    """
    Adds command-line arguments that are specific to *treesapp update*

    :return: None
    """
    up_parser.add_output_dir()
    up_parser.add_delete()
    up_parser.add_seq_params()  # w, m
    up_parser.add_cluster_args()  # p
    up_parser.add_taxa_args()  # s, f, t
    up_parser.add_refpkg_file_param()  # r
    up_parser.add_pplace_filter_params()
    up_parser.add_lineage_table_param()
    up_parser.add_phylogeny_params()  # b, e
    up_parser.add_compute_miscellany()  # n
    up_parser.add_basic_classifier_model_params()
    up_parser.reqs.add_argument("--treesapp_output", dest="ts_out", required=True,
                                help="Path to the directory containing TreeSAPP outputs, "
                                     "including sequences to be used for the update.")
    up_parser.optopt.add_argument('-i', '--fastx_input', required=False, dest="input", default=[''], nargs='+',
                                  help='An input file containing candidate reference sequences in either FASTA format. '
                                       'Will trigger re-training the reference package if provided.')
    up_parser.optopt.add_argument("--skip_assign", default=False, required=False, action="store_true",
                                  help="The assigned sequences are from a database and their database lineages "
                                       "should be used instead of the TreeSAPP-assigned lineages.")
    up_parser.optopt.add_argument("--resolve", default=False, required=False, action="store_true",
                                  help="Flag indicating candidate references with better resolved lineages and"
                                       " comparable sequence lengths can replace old references."
                                       " Useful when updating with sequences from isolates, SAGs and quality MAGs.")
    up_parser.optopt.add_argument("--stage", default="continue", required=False,
                                  choices=["continue", "lineages", "rebuild"],
                                  help="The stage(s) for TreeSAPP to execute [DEFAULT = continue]")
    up_parser.miscellany.add_argument("--headless", action="store_true", default=False,
                                      help="Do not require any user input during runtime.")
    return


def add_trainer_arguments(parser: TreeSAPPArgumentParser) -> None:
    """
    Adds command-line arguments that are specific to *treesapp train*

    :return: None
    """
    parser.add_io()
    parser.add_refpkg_file_param()
    parser.add_seq_params()
    parser.add_accession_params()
    parser.add_lineage_table_param()
    parser.add_taxa_ranks_param()
    parser.add_compute_miscellany()
    parser.add_basic_classifier_model_params()
    parser.add_advanced_classifier_params()
    parser.add_annot_map()

    parser.seqops.add_argument("-d", "--profile", required=False, default=False, action="store_true",
                               help="Input sequences will be purified with the reference package's profile HMM.")
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
    LOGGER.info("Arguments used:\n" + ' '.join(sys_args) + "\n")
    if re.match(r"^/$", args.output):
        LOGGER.error("Output directory specified as root. Bailing out to prevent future catastrophe!\n")
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
        LOGGER.error("Python 2 is not supported by TreeSAPP.\n")
        sys.exit(3)

    if "num_threads" in vars(args) and args.num_threads > available_cpu_count():
        LOGGER.warning("Number of threads specified is greater than those available! "
                        "Using maximum threads available (" + str(available_cpu_count()) + ")\n")
        args.num_threads = available_cpu_count()

    return


def check_evaluate_arguments(evaluator_instance: Evaluator, args) -> None:
    for rank in args.taxon_rank:
        evaluator_instance.ranks.append(rank)

    evaluator_instance.find_sequence_molecule_type()
    evaluator_instance.ref_pkg.f__pkl = args.pkg_path
    evaluator_instance.ref_pkg.slurp()

    if args.length:
        evaluator_instance.min_seq_length = str(min(args.length - 10, 30))
    else:
        evaluator_instance.min_seq_length = str(30)

    ##
    # Define locations of files TreeSAPP outputs
    ##
    evaluator_instance.test_rep_taxa_fasta = evaluator_instance.final_output_dir + "representative_taxa_sequences.fasta"
    evaluator_instance.performance_table = evaluator_instance.final_output_dir + "clade_exclusion_performance.tsv"
    evaluator_instance.recall_table = evaluator_instance.final_output_dir + "taxonomic_recall.tsv"
    evaluator_instance.containment_table = evaluator_instance.final_output_dir + "accuracy.tsv"

    if not os.path.isdir(evaluator_instance.var_output_dir):
        os.makedirs(evaluator_instance.var_output_dir)

    return


def check_updater_arguments(updater: Updater, args):
    updater.ref_pkg.f__pkl = args.pkg_path
    updater.ref_pkg.slurp()
    updater.updated_refpkg_path = os.path.join(updater.output_dir, "final_outputs",
                                               os.path.basename(updater.ref_pkg.f__pkl))
    updater.ref_pkg.disband(os.path.join(updater.output_dir, "intermediates"))
    updater.seq_names_to_taxa = args.seq_names_to_taxa
    # updater.rank_depth_map = {'k': 1, 'p': 2, 'c': 3, 'o': 4, 'f': 5, 'g': 6, 's': 7}

    if updater.input_sequences:
        updater.find_sequence_molecule_type()
    else:
        updater.molecule_type = updater.ref_pkg.molecule

    if args.similarity == 1.0:
        updater.prop_sim = updater.ref_pkg.pid
    else:
        updater.prop_sim = args.similarity

    if args.cluster:
        if not 0.5 <= float(args.similarity) <= 1.0:
            if 0.5 < float(args.similarity) / 100 < 1.0:
                args.similarity = str(float(args.similarity) / 100)
                LOGGER.warning("--similarity set to {} for compatibility.\n".format(args.similarity))
            else:
                LOGGER.error("--similarity {} is not between the supported range [0.5-1.0].\n".format(args.similarity))
                sys.exit(13)

    if updater.seq_names_to_taxa and not os.path.isfile(updater.seq_names_to_taxa):
        LOGGER.error("Unable to find file mapping sequence names to taxonomic lineages '" +
                      updater.seq_names_to_taxa + "'.\n")

    # TODO: Write a TreeSAPP function for validating outputs

    # Reset these values to reflect the paths for the TreeSAPP output that will be parsed
    updater.treesapp_output = args.ts_out
    if updater.treesapp_output[-1] != os.sep:
        updater.treesapp_output += os.sep
    updater.final_output_dir = updater.treesapp_output + "final_outputs" + os.sep
    # updater.var_output_dir = updater.treesapp_output + "intermediates" + os.sep
    updater.training_dir = updater.var_output_dir + "train"
    updater.old_ref_fasta = updater.var_output_dir + "original_refs.fasta"
    updater.combined_fasta = updater.var_output_dir + "all_refs.fasta"
    updater.lineage_map_file = updater.var_output_dir + "accession_id_lineage_map.tsv"
    updater.assignment_table = updater.final_output_dir + updater.classification_tbl_name
    updater.cluster_input = updater.var_output_dir + updater.sample_prefix + "_cluster_input.fasta"
    updater.clusters_prefix = updater.var_output_dir + updater.sample_prefix + "_cluster" + str(updater.prop_sim)
    classified_seqs = glob(updater.final_output_dir + "*_classified.faa")

    if len(classified_seqs) == 1:
        updater.query_sequences = classified_seqs.pop()
    elif len(classified_seqs) == 0:
        LOGGER.error("No classified sequence files found in {}.\n".format(updater.final_output_dir))
    else:
        LOGGER.error("Multiple classified sequence files in '{}'"
                      " where only one expected.\n".format(updater.final_output_dir))
        sys.exit(5)

    return
