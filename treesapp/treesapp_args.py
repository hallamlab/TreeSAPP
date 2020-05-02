import argparse
import os
import sys
import re
import logging
from glob import glob
from .classy import Assigner, Evaluator, Creator, PhyTrainer, Updater, Purity
from .utilities import available_cpu_count, get_refpkg_build


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
        super(TreeSAPPArgumentParser, self).__init__(add_help=False, **kwargs)
        self.reqs = self.add_argument_group("Required parameters")
        self.seqops = self.add_argument_group("Sequence operation arguments")
        self.rpkm_opts = self.add_argument_group("RPKM options")
        self.optopt = self.add_argument_group("Optional options")
        self.miscellany = self.add_argument_group("Miscellaneous options")
        self.taxa_args = self.add_argument_group("Taxonomic-lineage arguments")

        self.miscellany.add_argument("-v", "--verbose", action="store_true", default=False,
                                     help="Prints a more verbose runtime log")
        self.miscellany.add_argument("-h", "--help",
                                     action="help",
                                     help="Show this help message and exit")

    def parse_args(self, args=None, namespace=None):
        args = super(TreeSAPPArgumentParser, self).parse_args(args=args, namespace=namespace)

        return args

    # The following are building-block functions
    def add_io(self):
        self.reqs.add_argument('-i', '--fastx_input', required=True, dest="input",
                               help='An input file containing DNA or protein sequences in either FASTA or FASTQ format')
        self.optopt.add_argument('-o', '--output', default='./output/', required=False,
                                 help='Path to an output directory [DEFAULT = ./output/]')
        return

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

    def add_compute_miscellany(self):
        self.miscellany.add_argument('--overwrite', action='store_true', default=False,
                                     help='overwrites previously processed output folders')
        self.miscellany.add_argument('-n', '--num_procs', dest="num_threads", default=2, type=int,
                                     help='The number of CPU threads or parallel processes '
                                          'to use in various pipeline steps [DEFAULT = 2]')

    def add_accession_params(self):
        self.optopt.add_argument("--accession2taxid", dest="acc_to_taxid", required=False, default=None,
                                 help="Path to an NCBI accession2taxid file "
                                      "for more rapid accession-to-lineage mapping.\n")
        self.optopt.add_argument("-a", "--accession2lin", dest="acc_to_lin", required=False, default=None,
                                 help="Path to a file that maps sequence accessions to taxonomic lineages, "
                                      "possibly made by `treesapp create`...")

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
        self.taxa_args.add_argument("-t", "--min_taxonomic_rank",
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


def add_layer_arguments(parser: TreeSAPPArgumentParser):
    parser.reqs.add_argument("-o", "--treesapp_output", dest="output", required=True,
                             help="The TreeSAPP output directory.")
    parser.optopt.add_argument("-c", "--colours_style", required=False, nargs='+',
                               help="The colours_style file exported from iTOL with the annotation information. "
                                     "For the variable name to be automatically inferred (rather than through `names`). "
                                     "Format of the file should be `marker`_`var`.txt. For example: mcrA_Metabolism.txt "
                                     "would create a new column in marker_contig_map.tsv named 'Metabolism'.")
    parser.optopt.add_argument("-d", "--annot_dir", required=False, default=None,
                               help="Path to a directory containing iTOL annotation files for layering.")
    return


def add_classify_arguments(parser: TreeSAPPArgumentParser):
    """
    Returns the parser to interpret user options.
    """
    parser.add_io()
    parser.add_rpkm_params()
    parser.add_seq_params()
    parser.add_search_params()
    parser.add_compute_miscellany()
    # The required parameters... for which there are currently none. But they would go here!

    # The optionals
    parser.optopt.add_argument('-c', '--composition', default="meta", choices=["meta", "single"],
                               help="Sample composition being either a single organism or a metagenome.")
    parser.optopt.add_argument("-l", "--min_likelihood", default=0.1, type=float,
                               help="The minimum likelihood weight ratio required for a RAxML placement. "
                               "[DEFAULT = 0.1]")
    parser.optopt.add_argument("-P", "--placement_parser", default="best", type=str, choices=["best", "lca"],
                               help="Algorithm used for parsing each sequence's potential RAxML placements. "
                               "[DEFAULT = 'best']")
    parser.optopt.add_argument('-t', '--targets', default='', type=str,
                               help='A comma-separated list specifying which marker genes to query in input by'
                               ' the "denominator" column in data/tree_data/cog_list.tsv'
                               ' - e.g., M0701,D0601 for mcrA and nosZ\n[DEFAULT = ALL]')
    parser.optopt.add_argument("--stage", default="continue", required=False,
                               choices=["continue", "orf-call", "search", "align", "place", "classify"],
                               help="The stage(s) for TreeSAPP to execute [DEFAULT = continue]")
    parser.rpkm_opts.add_argument("--rpkm", action="store_true", default=False,
                                  help="Flag indicating RPKM values should be calculated for the sequences detected")

    # The miscellany
    parser.miscellany.add_argument('-R', '--reftree', default='p', type=str,
                                   help='Reference tree (p = MLTreeMap reference phylogenetic tree [DEFAULT])'
                                   ' Change to code to map query sequences to specific phylogenetic tree.')
    parser.miscellany.add_argument("--check_trees", action="store_true", default=False,
                                   help="Quality-check the reference trees before running TreeSAPP")
    parser.miscellany.add_argument('-d', '--delete', default=False, action="store_true",
                                   help='Delete intermediate file to save disk space. '
                                        'Recommended for large metagenomes!')

    return


def add_abundance_arguments(parser: TreeSAPPArgumentParser):
    parser.add_rpkm_params()
    parser.add_compute_miscellany()
    parser.reqs.add_argument("--treesapp_output", dest="output", required=True,
                             help="Path to the directory containing TreeSAPP outputs, "
                                  "including sequences to be used for the update.")
    # TODO: Include an option to append new values to the classification table
    parser.optopt.add_argument("--report", choices=["update", "nothing"], required=False, default="nothing",
                               help="What should be done with the abundance values? The TreeSAPP classification table "
                                    "can overwritten (update) or left unchanged. "
                                    "[ DEFAULT = nothing ]")


def add_create_arguments(parser: TreeSAPPArgumentParser):
    parser.add_io()
    parser.add_seq_params()
    parser.add_taxa_args()
    parser.add_accession_params()
    parser.add_compute_miscellany()
    # The required parameters
    parser.reqs.add_argument("-c", "--refpkg_name",
                             help="Unique name to be used by TreeSAPP internally. NOTE: Must be <=6 characters.\n"
                                  "Examples are 'McrA', 'DsrAB', and 'p_amoA'.",
                             required=True)
    parser.reqs.add_argument("-p", "--identity",
                             help="Fractional identity value (between 0.50 and 1.0)\n"
                                  "the input sequences were clustered at.",
                             required=True,
                             type=str)

    parser.seqops.add_argument("--cluster",
                               help="Flag indicating usearch should be used to cluster sequences\n"
                                    "at the fractional similarity indicated by identity (`-p`)",
                               action="store_true",
                               default=False)
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
    parser.seqops.add_argument('-r', "--rfam_cm",
                               help="The covariance model of the RNA family being packaged.\n"
                                    "REQUIRED if molecule is rRNA!",
                               default=None,
                               required=False)

    parser.optopt.add_argument("-b", "--bootstraps",
                               help="The number of bootstrap replicates RAxML should perform\n"
                                    "[ DEFAULT = autoMR ]",
                               required=False, default="autoMR")
    parser.optopt.add_argument("-e", "--raxml_model",
                               help="The evolutionary model for RAxML to use\n"
                                    "[ Proteins = PROTGAMMAAUTO | Nucleotides =  GTRGAMMA ]",
                               required=False, default=None)
    parser.optopt.add_argument("-u", "--uc",
                               help="The USEARCH cluster format file produced from clustering reference sequences.\n"
                                    "This can be used for selecting representative headers from identical sequences.",
                               required=False, default=None)
    parser.optopt.add_argument("--fast",
                               help="A flag indicating the tree should be built rapidly, using FastTree.",
                               default=False, required=False,
                               action="store_true")
    parser.optopt.add_argument("--kind", default="functional", choices=["functional", "taxonomic"], required=False,
                               help="The broad classification of marker gene type, either "
                                    "functional or taxonomic. [ DEFAULT = functional ]")
    parser.optopt.add_argument("--stage", default="continue", required=False,
                               choices=["continue", "lineages", "clean", "cluster", "build", "train", "cc"],
                               help="The stage(s) for TreeSAPP to execute [DEFAULT = continue]")

    parser.miscellany.add_argument('--pc', action='store_true', default=False,
                                   help='Prints the final commands to complete\n'
                                        'installation for a provided `code_name`')
    parser.miscellany.add_argument("--headless", action="store_true", default=False,
                                   help="Do not require any user input during runtime.")


def add_purity_arguments(parser: TreeSAPPArgumentParser):
    parser.add_io()
    parser.add_seq_params()
    parser.add_compute_miscellany()
    parser.reqs.add_argument("-r", "--reference_marker", dest="refpkg", required=True,
                             help="Short-form name of the marker gene to be tested (e.g. mcrA, pmoA, nosZ)")
    parser.reqs.add_argument("-p", "--pkg_path", dest="pkg_path", required=True,
                             help="Path to the reference package.\n")
    parser.optopt.add_argument("-x", "--extra_info", required=False, default=None,
                               help="File mapping header prefixes to description information.")
    parser.optopt.add_argument("--stage", default="continue", required=False,
                               choices=["continue", "lineages", "classify", "calculate"],
                               help="The stage(s) for TreeSAPP to execute [DEFAULT = continue]")
    # TODO: Remove --trim_align from command-line options in parser.add_seq_params()
    return


def add_evaluate_arguments(parser: TreeSAPPArgumentParser):
    parser.add_io()
    parser.add_seq_params()
    parser.add_accession_params()
    parser.add_compute_miscellany()
    parser.reqs.add_argument("-r", "--reference_marker",
                             help="Short-form name of the marker gene to be tested (e.g. mcrA, pmoA, nosZ)",
                             required=True)

    parser.optopt.add_argument("--fresh", default=False, required=False, action="store_true",
                               help="Recalculate a fresh phylogenetic tree with the target clades removed instead of"
                                    " removing the leaves corresponding to targets from the reference tree.")
    parser.optopt.add_argument("--tool", default="treesapp", required=False,
                               choices=["treesapp", "graftm", "diamond"],
                               help="Classify using one of the tools: treesapp [DEFAULT], graftm, or diamond.")
    parser.optopt.add_argument("-t", "--taxon_rank",
                               help="Comma-separated list of the taxonomic ranks to test " +
                                    "choices = [Phylum, Class, Order, Family, Genus, Species] (DEFAULT = Species)",
                               default="Species",
                               required=False)
    parser.optopt.add_argument("-l", "--length",
                               required=False, type=int, default=0,
                               help="Arbitrarily slice the input sequences to this length. "
                                    "Useful for testing classification accuracy for fragments.")
    parser.optopt.add_argument("--stage", default="continue", required=False,
                               choices=["continue", "lineages", "classify", "calculate"],
                               help="The stage(s) for TreeSAPP to execute [DEFAULT = continue]")
    return


def add_update_arguments(parser: TreeSAPPArgumentParser):
    parser.add_io()
    parser.add_seq_params()
    parser.add_taxa_args()
    parser.add_compute_miscellany()
    parser.reqs.add_argument("-c", "--refpkg_name", dest="name", required=True,
                             help="Unique name to be used by TreeSAPP internally. NOTE: Must be <=6 characters.\n"
                                  "Examples are 'McrA', 'DsrAB', and 'p_amoA'.")
    parser.reqs.add_argument("--treesapp_output", dest="ts_out", required=True,
                             help="Path to the directory containing TreeSAPP outputs, "
                                  "including sequences to be used for the update.")
    parser.optopt.add_argument("-l", "--min_lwr", dest="min_lwr", required=False, default=0.0, type=float,
                               help="The minimum likelihood weight ratio for a sequence to be included in update.")
    parser.optopt.add_argument("-a", "--seqs2taxa", dest="seq_names_to_taxa", required=False, default=None,
                               help="Path to a file mapping sequence names (i.e. contig headers) to taxonomic lineages")
    parser.optopt.add_argument("--fast", default=False, required=False, action="store_true",
                               help="A flag indicating the tree should be built rapidly, using FastTree.")
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
    parser.seqops.add_argument("--cluster", required=False, default=False, action="store_true",
                               help="Cluster sequences that mapped to the reference tree prior to updating")
    parser.seqops.add_argument("-p", "--identity", required=False, type=float,
                               help="Fractional similarity (between 0.50 and 1.0) to cluster sequences.")
    parser.miscellany.add_argument("--headless", action="store_true", default=False,
                                   help="Do not require any user input during runtime.")


def add_trainer_arguments(parser: TreeSAPPArgumentParser):
    parser.add_io()
    parser.add_seq_params()
    parser.add_accession_params()
    parser.add_compute_miscellany()
    parser.reqs.add_argument("-c", "--refpkg_name", dest="name", required=True,
                             help="Unique name to be used by TreeSAPP internally. NOTE: Must be <=6 characters.\n"
                                  "Examples are 'McrA', 'DsrAB', and 'p_amoA'.")
    parser.reqs.add_argument("-p", "--pkg_path", required=True,
                             help="Path to the reference package.\n")
    parser.seqops.add_argument("-d", "--profile", required=False, default=False, action="store_true",
                               help="Flag indicating input sequences need to be purified using an HMM profile.")
    parser.optopt.add_argument("--stage", default="continue", required=False,
                               choices=["continue", "search", "lineages", "place", "regress"],
                               help="The stage(s) for TreeSAPP to execute [DEFAULT = continue]")


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

    if args.num_threads > available_cpu_count():
        logging.warning("Number of threads specified is greater than those available! "
                        "Using maximum threads available (" + str(available_cpu_count()) + ")\n")
        args.num_threads = available_cpu_count()

    return


def check_purity_arguments(purity_instance: Purity, args, marker_build_dict: dict):
    purity_instance.ref_pkg.prefix = args.refpkg
    purity_instance.refpkg_build = get_refpkg_build(purity_instance.ref_pkg.prefix,
                                                    marker_build_dict,
                                                    purity_instance.refpkg_code_re)
    purity_instance.pkg_path = args.pkg_path
    purity_instance.ref_pkg.gather_package_files(purity_instance.pkg_path, purity_instance.molecule_type)

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


def check_evaluate_arguments(evaluator_instance: Evaluator, args, marker_build_dict):
    taxa_choices = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]
    args.taxon_rank = args.taxon_rank.split(',')
    for rank in args.taxon_rank:
        if rank not in taxa_choices:
            logging.error(rank + " not an available option for `--taxon_rank`.\n")
            sys.exit(21)
        else:
            evaluator_instance.ranks.append(rank)
    evaluator_instance.target_marker = get_refpkg_build(args.reference_marker,
                                                        marker_build_dict,
                                                        evaluator_instance.refpkg_code_re)
    evaluator_instance.targets.append(evaluator_instance.target_marker.denominator)

    if args.acc_to_lin:
        evaluator_instance.acc_to_lin = args.acc_to_lin
        if os.path.isfile(evaluator_instance.acc_to_lin):
            evaluator_instance.change_stage_status("lineages", False)
        else:
            logging.error("Unable to find accession-lineage mapping file '" + evaluator_instance.acc_to_lin + "'\n")
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


def check_trainer_arguments(trainer_instance: PhyTrainer, args, marker_build_dict):
    target_marker_build = get_refpkg_build(args.name,
                                           marker_build_dict,
                                           trainer_instance.refpkg_code_re)
    trainer_instance.ref_pkg.prefix = target_marker_build.cog
    trainer_instance.ref_pkg.refpkg_code = target_marker_build.denominator

    ##
    # Define locations of files TreeSAPP outputs
    ##
    trainer_instance.placement_table = trainer_instance.output_dir + "placement_info.tsv"
    trainer_instance.placement_summary = trainer_instance.output_dir + "placement_trainer_results.txt"
    trainer_instance.hmm_purified_seqs = trainer_instance.output_dir + trainer_instance.ref_pkg.prefix + "_hmm_purified.fasta"

    if not os.path.isdir(trainer_instance.var_output_dir):
        os.makedirs(trainer_instance.var_output_dir)

    return


def check_classify_arguments(assigner: Assigner, args):
    """
    Ensures the command-line arguments returned by argparse are sensible
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
        for marker in assigner.target_refpkgs:
            if not assigner.refpkg_code_re.match(marker):
                logging.error("Incorrect format for target: " + str(marker) +
                              "\nRefer to column 'Denominator' in " + assigner.treesapp_dir +
                              "data/ref_build_parameters.tsv for identifiers that can be used.\n")
                sys.exit(3)
    else:
        assigner.target_refpkgs = []

    if args.molecule == "prot":
        assigner.change_stage_status("orf-call", False)
        if args.rpkm:
            logging.error("Unable to calculate RPKM values for protein sequences.\n")
            sys.exit(3)

    # TODO: transfer all of this HMM-parsing stuff to the assigner_instance
    # Parameterizing the hmmsearch output parsing:
    args.perc_aligned = 10
    args.min_acc = 0.7
    if args.stringency == "relaxed":
        args.min_e = 1E-3
        args.min_ie = 1E-1
        args.min_score = 15
    elif args.stringency == "strict":
        args.min_e = 1E-7
        args.min_ie = 1E-5
        args.min_score = 30
    else:
        logging.error("Unknown HMM-parsing stringency argument '" + args.stringency + "'.\n")
        sys.exit(3)

    return args


def check_create_arguments(creator: Creator, args):
    creator.ref_pkg.prefix = args.refpkg_name
    if not args.output:
        args.output = os.getcwd() + os.sep + creator.ref_pkg.prefix + "_treesapp_refpkg" + os.sep

    if len(creator.ref_pkg.prefix) > 6:
        logging.error("Name must be <= 6 characters!\n")
        sys.exit(13)

    if args.rfam_cm is None and args.molecule == "rrna":
        logging.error("Covariance model file must be provided for rRNA data!\n")
        sys.exit(13)

    # Check the RAxML model
    raxml_models = ["PROTGAMMAWAG", "PROTGAMMAAUTO", "PROTGAMMALG", "GTRCAT", "GTRCATI ", "GTRCATX", "GTRGAMMA",
                    "ASC_GTRGAMMA", "ASC_GTRCAT", "BINGAMMA", "PROTGAMMAILGX", "PROTGTRGAMMA"]
    if args.raxml_model:
        valid = False
        for model in raxml_models:
            if re.search(model, args.raxml_model):
                valid = True
                break
        if not valid:
            logging.error("Phylogenetic substitution model '" + args.raxml_model + "' is not valid!\n" +
                          "If this model is valid (not a typo), add it to `raxml_models` list and re-run.\n")
            sys.exit(13)
        else:
            creator.ref_pkg.sub_model = args.raxml_model

    if args.cluster:
        if args.multiple_alignment:
            logging.error("--cluster and --multiple_alignment are mutually exclusive!\n")
            sys.exit(13)
        if args.uc:
            logging.error("--cluster and --uc are mutually exclusive!\n")
            sys.exit(13)
        if not 0.5 <= float(args.identity) <= 1.0:
            if 0.5 < float(args.identity)/100 < 1.0:
                args.identity = str(float(args.identity)/100)
                logging.warning("--identity  set to " + args.identity + " for compatibility with USEARCH \n")
            else:
                logging.error("--identity " + args.identity + " is not between the supported range [0.5-1.0]\n")
                sys.exit(13)
        creator.prop_sim = args.identity

    if args.taxa_lca:
        if not args.cluster and not args.uc:
            logging.error("Unable to perform LCA for representatives without clustering information: " +
                          "either with a provided UCLUST file or by clustering within the pipeline.\n")
            sys.exit(13)

    if args.guarantee:
        if not args.uc and not args.cluster:
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
    creator.hmm_purified_seqs = creator.var_output_dir + creator.ref_pkg.prefix + "_hmm_purified.fasta"
    creator.filtered_fasta = creator.var_output_dir + creator.sample_prefix + "_filtered.fa"
    creator.cluster_input = creator.var_output_dir + creator.sample_prefix + "_uclust_input.fasta"
    creator.uclust_prefix = creator.var_output_dir + creator.sample_prefix + "_uclust" + str(creator.prop_sim)
    creator.unaln_ref_fasta = creator.var_output_dir + creator.ref_pkg.prefix + "_ref.fa"
    creator.phylip_file = creator.var_output_dir + creator.ref_pkg.prefix + ".phy"

    return


def check_updater_arguments(updater: Updater, args, marker_build_dict):
    updater.ref_pkg.prefix = args.name
    updater.seq_names_to_taxa = args.seq_names_to_taxa
    updater.rank_depth_map = {'k': 1, 'p': 2, 'c': 3, 'o': 4, 'f': 5, 'g': 6, 's': 7}
    updater.target_marker = get_refpkg_build(updater.ref_pkg.prefix, marker_build_dict, updater.refpkg_code_re)
    if not args.identity:
        updater.prop_sim = updater.target_marker.pid
    else:
        updater.prop_sim = args.identity

    if args.cluster:
        if not 0.5 <= float(updater.prop_sim) <= 1.0:
            if 0.5 < float(updater.prop_sim)/100 < 1.0:
                updater.prop_sim = str(float(updater.prop_sim) / 100)
                logging.warning("--identity  set to " + updater.prop_sim + " for compatibility with USEARCH \n")
            else:
                logging.error("--identity " + updater.prop_sim + " is not between the supported range [0.5-1.0]\n")
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
    updater.var_output_dir = updater.treesapp_output + "intermediates" + os.sep
    updater.old_ref_fasta = updater.output_dir + "original_refs.fasta"
    updater.combined_fasta = updater.output_dir + "all_refs.fasta"
    updater.lineage_map_file = updater.output_dir + "accession_id_lineage_map.tsv"
    updater.assignment_table = updater.final_output_dir + "marker_contig_map.tsv"
    updater.cluster_input = updater.var_output_dir + updater.sample_prefix + "_uclust_input.fasta"
    updater.uclust_prefix = updater.var_output_dir + updater.sample_prefix + "_uclust" + str(updater.prop_sim)
    classified_seqs = glob(updater.final_output_dir + "*_classified.faa")
    if len(classified_seqs) == 1:
        updater.query_sequences = classified_seqs.pop()
    else:
        logging.error("Multiple classified sequence files in '" + updater.final_output_dir + "' but only 1 expected.\n")
        sys.exit(5)

    return
