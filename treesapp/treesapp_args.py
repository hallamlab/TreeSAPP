import argparse
import os
import sys
import re
import logging
from .classy import Assigner, Evaluator, Creator
from .utilities import available_cpu_count, check_previous_output, get_refpkg_build


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
        self.reqs.add_argument('-i', '--fasta_input', required=True, dest="input",
                               help='An input file containing DNA or protein sequences in FASTA format')
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
        self.rpkm_opts.add_argument("--rpkm", action="store_true", default=False,
                                    help="Flag indicating RPKM values should be calculated for the sequences detected")
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
        self.miscellany.add_argument('-d', '--delete', default=False, action="store_true",
                                     help='Delete intermediate file to save disk space\n'
                                          'Recommended for large metagenomes!')
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
    parser.optopt.add_argument("--update_tree", action="store_true", default=False,
                               help="Flag indicating the reference tree specified by `--reftree` "
                                    "is to be updated using the sequences found in TreeSAPP output")
    parser.optopt.add_argument("--stage", default="continue", required=False,
                               choices=["continue", "orf-call", "search", "align", "place", "classify"],
                               help="The stage(s) for TreeSAPP to execute [DEFAULT = continue]")

    # The miscellany
    parser.miscellany.add_argument('-R', '--reftree', default='p', type=str,
                                   help='Reference tree (p = MLTreeMap reference phylogenetic tree [DEFAULT])'
                                   ' Change to code to map query sequences to specific phylogenetic tree.')
    parser.miscellany.add_argument("--check_trees", action="store_true", default=False,
                                   help="Quality-check the reference trees before running TreeSAPP")

    return


def add_create_arguments(parser: TreeSAPPArgumentParser):
    parser.add_io()
    parser.add_seq_params()
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
    parser.seqops.add_argument("-d", "--domain",
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

    taxa_args = parser.add_argument_group("Taxonomic-lineage arguments")
    taxa_args.add_argument("-s", "--screen",
                           help="Keywords for including specific taxa in the tree.\n"
                                "Example: to only include Bacteria and Archaea do `--screen Bacteria,Archaea`\n"
                                "[ DEFAULT is no screen ]",
                           default="", required=False)
    taxa_args.add_argument("-f", "--filter",
                           help="Keywords for removing specific taxa; the opposite of `--screen`.\n"
                                "[ DEFAULT is no filter ]",
                           default="", required=False)
    taxa_args.add_argument("-t", "--min_taxonomic_rank",
                           required=False, default='k', choices=['k', 'p', 'c', 'o', 'f', 'g', 's'],
                           help="The minimum taxonomic lineage resolution for reference sequences [ DEFAULT = k ].\n")
    taxa_args.add_argument("--taxa_lca",
                           help="Set taxonomy of representative sequences to LCA of cluster member's taxa.\n"
                                "[ --cluster or --uc REQUIRED ]",
                           default=False, required=False, action="store_true")
    taxa_args.add_argument("--taxa_norm",
                           help="[ IN DEVELOPMENT ] Perform taxonomic normalization on the provided sequences.\n"
                                "A comma-separated argument with the Rank (e.g. Phylum) and\n"
                                "number of representatives is required.\n")

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

    parser.miscellany.add_argument('--pc', action='store_true', default=False,
                                   help='Prints the final commands to complete\n'
                                        'installation for a provided `code_name`')
    parser.miscellany.add_argument("--headless", action="store_true", default=False,
                                   help="Do not require any user input during runtime.")


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
    parser.optopt.add_argument("--cluster", required=False, default=False, action="store_true",
                               help="Cluster sequences that mapped to the reference tree prior to updating")


def check_parser_arguments(args):
    """
    Function for checking arguments that are found in args.namespace()
    This is the only parser validation function used by clade exclusion evaluator
    :param args:
    :return:
    """
    ##
    # Remove the output directory if it exists and overwrite permission granted.
    ##
    if re.match(r"^/$", args.output):
        logging.error("Output directory specified as root. Bailing out to prevent future catastrophe!\n")
        sys.exit(1)
    # Add (or replace a trailing (back)slash with) the os.sep to the end of the output directory
    while re.search(r'.*/$', args.output) or re.search(r'.*\\$', args.output):
        args.output = args.output[:-1]
    args.output += os.sep
    if not re.match(r'^/.*', args.output):
        args.output = os.getcwd() + os.sep + args.output  # args.output is now the absolute path

    args.var_output_dir = args.output + 'intermediates' + os.sep
    args.final_output_dir = args.output + 'final_outputs' + os.sep

    # Determine whether the output directory should be removed, creates output directory if it doesn't exist
    check_previous_output(args)

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

    if args.acc_to_lin and os.path.isfile(args.acc_to_lin):
        evaluator_instance.acc_to_lin = args.acc_to_lin
        evaluator_instance.change_stage_status("lineages", False)
    else:
        evaluator_instance.acc_to_lin = args.output + os.sep + "accession_id_lineage_map.tsv"

    ##
    # Define locations of files TreeSAPP outputs
    ##
    evaluator_instance.test_rep_taxa_fasta = args.output + os.sep + "representative_taxa_sequences.fasta"
    evaluator_instance.performance_table = args.output + os.sep + "clade_exclusion_performance.tsv"
    evaluator_instance.containment_table = args.output + os.sep + "accuracy.tsv"
    evaluator_instance.var_output_dir = args.output + "TreeSAPP_output" + os.sep

    if not os.path.isdir(evaluator_instance.var_output_dir):
        os.makedirs(evaluator_instance.var_output_dir)

    return


def check_classify_arguments(assigner_instance: Assigner, args):
    """
    Ensures the command-line arguments returned by argparse are sensible
    :param assigner_instance: An instantiated Assigner object
    :param args: object with parameters returned by argparse.parse_args()
    :return: 'args', a summary of TreeSAPP settings.
    """
    if args.targets:
        assigner_instance.target_refpkgs = args.targets.split(',')
        for marker in assigner_instance.target_refpkgs:
            if not assigner_instance.refpkg_code_re.match(marker):
                logging.error("Incorrect format for target: " + str(marker) +
                              "\nRefer to column 'Denominator' in " + assigner_instance.treesapp_dir +
                              "data/ref_build_parameters.tsv for identifiers that can be used.\n")
                sys.exit(3)
    else:
        assigner_instance.target_refpkgs = []

    if args.molecule == "prot":
        assigner_instance.change_stage_status("orf-call", False)
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


def check_create_arguments(create_instance: Creator, args):
    if not args.output:
        args.output = os.getcwd() + os.sep + args.code_name + "_treesapp_refpkg"

    # Names of files and directories to be created
    create_instance.refpkg_output = create_instance.final_output_dir + "TreeSAPP_files_%s" % args.code_name + os.sep
    create_instance.phy_dir = create_instance.output_dir + "phylogeny_files" + os.sep
    create_instance.acc_to_lin = create_instance.output_dir + os.sep + "accession_id_lineage_map.tsv"
    create_instance.hmm_purified_seqs = create_instance.output_dir + create_instance.refpkg_name + "_hmm_purified.fasta"
    create_instance.filtered_fasta = create_instance.output_dir + create_instance.sample_prefix + "_filtered.fa"
    create_instance.uclust_prefix = create_instance.output_dir + create_instance.sample_prefix + "_uclust" + create_instance.prop_sim
    create_instance.unaln_ref_fasta = create_instance.output_dir + create_instance.refpkg_name + "_ref.fa"
    create_instance.phylip_file = create_instance.output_dir + create_instance.refpkg_name + ".phy"

    if len(args.code_name) > 6:
        logging.error("Name must be <= 6 characters!\n")
        sys.exit(13)

    if args.rfam_cm is None and args.molecule == "rrna":
        logging.error("Covariance model file must be provided for rRNA data!\n")
        sys.exit(13)

    # Check the RAxML model
    raxml_models = ["PROTGAMMAWAG", "PROTGAMMAAUTO", "PROTGAMMALGF", "GTRCAT", "GTRCATI ", "GTRCATX", "GTRGAMMA",
                    "ASC_GTRGAMMA", "ASC_GTRCAT", "BINGAMMA", "PROTGAMMAILGX", "PROTGTRGAMMA"]
    if args.raxml_model:
        if args.raxml_model not in raxml_models:
            logging.error("Phylogenetic substitution model '" + args.raxml_model + "' is not valid!\n" +
                          "If this model is valid (not a typo), add it to `raxml_models` list and re-run.\n")
            sys.exit(13)
        else:
            create_instance.ref_pkg.sub_model = args.raxml_model

    if args.cluster:
        if args.multiple_alignment:
            logging.error("--cluster and --multiple_alignment are mutually exclusive!\n")
            sys.exit(13)
        if args.uc:
            logging.error("--cluster and --uc are mutually exclusive!\n")
            sys.exit(13)
        if not 0.5 < float(args.identity) < 1.0:
            if 0.5 < float(args.identity)/100 < 1.0:
                args.identity = str(float(args.identity)/100)
                logging.warning("--identity  set to " + args.identity + " for compatibility with USEARCH \n")
            else:
                logging.error("--identity " + args.identity + " is not between the supported range [0.5-1.0]\n")
                sys.exit(13)

    if args.taxa_lca:
        if not args.cluster and not args.uc:
            logging.error("Unable to perform LCA for representatives without clustering information: " +
                          "either with a provided UCLUST file or by clustering within the pipeline.\n")
            sys.exit(13)

    if args.guarantee:
        if not args.uc and not args.cluster:
            logging.error("--guarantee used but without clustering there is no reason for it.\n" +
                          "Include all sequences in " + args.guarantee +
                          " in " + args.fasta_input + " and re-run without --guarantee\n")
            sys.exit(13)

    return

