import argparse
import os
import sys
import re
import logging
import shutil
from .classy import prep_logging
from .utilities import find_executables, executable_dependency_versions, available_cpu_count, check_previous_output


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
        self.reqs.add_argument('-i', '--fasta_input', required=True,
                               help='An input file containing DNA or protein sequences in FASTA format')
        self.optopt.add_argument('-o', '--output', default='./output/', required=False,
                                 help='Path to an output directory [DEFAULT = ./output/]')
        return

    def add_seq_params(self):
        self.optopt.add_argument("--trim_align", default=False, action="store_true",
                                 help="Flag to turn on position masking of the multiple sequence alignment"
                                      " [DEFAULT = False]")
        self.optopt.add_argument('-g', '--min_seq_length', default=30, type=int,
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
        self.add_argument("-s", "--stringency", choices=["relaxed", "strict"], default="relaxed", required=False,
                          help="HMM-threshold mode affects the number of query sequences that advance.")

    def add_compute_miscellany(self):
        self.miscellany.add_argument('--overwrite', action='store_true', default=False,
                                     help='overwrites previously processed output folders')
        self.miscellany.add_argument('-T', '--num_threads', default=2, type=int,
                                     help='specifies the number of CPU threads to use in RAxML and BLAST '
                                          'and processes throughout the pipeline [DEFAULT = 2]')
        self.miscellany.add_argument('-d', '--delete', default=False, action="store_true",
                                     help='Delete intermediate file to save disk space\n'
                                          'Recommended for large metagenomes!')


def add_classify_arguments(parser: TreeSAPPArgumentParser):
    """
    Returns the parser to interpret user options.
    """
    parser.add_io()
    parser.add_rpkm_params()
    parser.add_seq_params()
    parser.add_search_params()
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
    parser.optopt.add_argument('-t', '--targets', default='ALL', type=str,
                               help='A comma-separated list specifying which marker genes to query in input by'
                               ' the "denominator" column in data/tree_data/cog_list.tsv'
                               ' - e.g., M0701,D0601 for mcrA and nosZ\n[DEFAULT = ALL]')
    parser.optopt.add_argument("--update_tree", action="store_true", default=False,
                               help="Flag indicating the reference tree specified by `--reftree` "
                                    "is to be updated using the sequences found in TreeSAPP output")
    parser.optopt.add_argument("--reclassify", action="store_true", default=False,
                               help="Flag indicating current outputs should be used to generate "
                                    "all outputs downstream of phylogenetic placement.")

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
    parser.reqs.add_argument("-c", "--code_name",
                             help="Unique name to be used by TreeSAPP internally. NOTE: Must be <=6 characters.\n"
                                  "(Refer to first column of 'cog_list.txt' under the '#functional cogs' section)",
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
                           help="Keywords (taxonomic regular expressions) for including specific taxa in the tree.\n"
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
    taxa_args.add_argument("--accession2taxid", required=False, default=None,
                           help="Path to an NCBI accession2taxid file for more rapid accession-to-lineage mapping.\n")

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

    # TODO: fix this... perhaps import treesapp, treesapp.__path__?
    args.treesapp = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + os.sep

    args.min_seq_length = 1

    if sys.version_info > (2, 9):
        args.py_version = 3
    else:
        logging.error("Python 2 is not supported by TreeSAPP.\n")
        sys.exit(3)

    return


def check_classify_arguments(args):
    """
    Ensures the command-line arguments returned by argparse are sensible
    :param args: object with parameters returned by argparse.parse_args()
    :return: 'args', a summary of TreeSAPP settings.
    """
    # Setup the global logger and main log file
    log_file_name = args.output + os.sep + "TreeSAPP_log.txt"
    prep_logging(log_file_name, args.verbose)
    logging.debug("Command used:\n" + ' '.join(sys.argv) + "\n")

    # Set the reference data file prefix and the reference tree name
    args.reference_tree = args.reftree

    args.targets = args.targets.split(',')
    if args.targets != ['ALL']:
        for marker in args.targets:
            if not re.match('[A-Z][0-9]{4}', marker):
                logging.error("Incorrect format for target: " + str(marker) +
                              "\nRefer to column 'Denominator' in " + args.treesapp + "data/tree_data/" +
                              "cog_list.tsv for identifiers that can be used.\n")
                sys.exit()

    args = find_executables(args)
    logging.debug(executable_dependency_versions(args.executables))

    if args.num_threads > available_cpu_count():
        logging.warning("Number of threads specified is greater than those available! "
                        "Using maximum threads available (" + str(available_cpu_count()) + ")\n")
        args.num_threads = available_cpu_count()

    if args.rpkm:
        if not args.reads:
            logging.error("At least one FASTQ file must be provided if -rpkm flag is active!")
            sys.exit()
        if args.reverse and not args.reads:
            logging.error("File containing reverse reads provided but forward mates file missing!")
            sys.exit()

    if args.molecule == "prot" and args.rpkm:
        logging.error("Unable to calculate RPKM values for protein sequences.\n")
        sys.exit()

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


def check_create_arguments(args):
    args.treesapp = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + os.sep
    if not args.output:
        args.output = os.getcwd() + os.sep + args.code_name + "_treesapp_refpkg"
    args.final_output_dir = args.output + "TreeSAPP_files_%s" % args.code_name + os.sep

    if len(args.code_name) > 6:
        logging.error("code_name must be <= 6 characters!\n")
        sys.exit(13)

    if args.rfam_cm is None and args.molecule == "rrna":
        logging.error("Covariance model file must be provided for rRNA data!\n")
        sys.exit(13)

    # Check the RAxML model
    raxml_models = ["PROTGAMMAWAG", "PROTGAMMAAUTO", "PROTGAMMALGF", "GTRCAT", "GTRCATI ", "GTRCATX", "GTRGAMMA",
                    "ASC_GTRGAMMA", "ASC_GTRCAT", "BINGAMMA", "PROTGAMMAILGX", "PROTGTRGAMMA"]
    if args.raxml_model and args.raxml_model not in raxml_models:
        sys.stderr.write("ERROR: --raxml_model (" + args.raxml_model + ") not valid!\n")
        sys.stderr.write("If this model is valid (not a typo), add if to `raxml_models` list and re-run.\n")
        sys.exit(13)

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
                sys.stderr.write("WARNING: --identity  set to " + args.identity + " for compatibility with USEARCH \n")
            else:
                sys.exit("ERROR: --identity " + args.identity + " is not between the supported range [0.5-1.0]\n")

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

    return args

