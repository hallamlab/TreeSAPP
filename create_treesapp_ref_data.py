#!/usr/bin/env python3

__author__ = "Connor Morgan-Lang"
__maintainer__ = "Connor Morgan-Lang"
__license__ = "GPL"
__version__ = "0.1.0"

try:
    import argparse
    import logging
    import sys
    import os
    import shutil
    import re
    import traceback

    from time import gmtime, strftime, sleep

    from utilities import os_type, which, find_executables, reformat_string, return_sequence_info_groups,\
        reformat_fasta_to_phy, write_phy_file, cluster_sequences
    from fasta import format_read_fasta, get_headers, get_header_format, write_new_fasta, summarize_fasta_sequences,\
        trim_multiple_alignment, read_fasta_to_dict
    from classy import ReferenceSequence, ReferencePackage, Header, Cluster, MarkerBuild,\
        prep_logging, register_headers, get_header_info
    from external_command_interface import launch_write_command
    from entish import annotate_partition_tree
    from lca_calculations import megan_lca, lowest_common_taxonomy, clean_lineage_list
    from entrez_utils import get_multiple_lineages, get_lineage_robust, verify_lineage_information,\
        read_accession_taxa_map, write_accession_lineage_map, build_entrez_queries
    from file_parsers import parse_domain_tables, read_phylip_to_dict, read_uc
    from placement_trainer import regress_rank_distance

except ImportError:
    sys.stderr.write("Could not load some user defined module functions:\n")
    sys.stderr.write(str(traceback.print_exc(10)))
    sys.exit(13)


def get_arguments():
    parser = argparse.ArgumentParser(add_help=False, formatter_class=argparse.RawTextHelpFormatter)
    required_args = parser.add_argument_group("Required arguments")
    required_args.add_argument("-i", "--fasta_input",
                               help="FASTA file that will be used to create reference data for TreeSAPP",
                               required=True)
    required_args.add_argument("-c", "--code_name",
                               help="Unique name to be used by TreeSAPP internally. NOTE: Must be <=6 characters.\n"
                                    "(Refer to first column of 'cog_list.txt' under the '#functional cogs' section)",
                               required=True)
    required_args.add_argument("-p", "--identity",
                               help="Fractional identity value (between 0.50 and 1.0)\n"
                                    "the input sequences were clustered at.",
                               required=True,
                               type=str)

    seqop_args = parser.add_argument_group("Sequence operation arguments")
    seqop_args.add_argument("--cluster",
                            help="Flag indicating usearch should be used to cluster sequences\n"
                                 "at the fractional similarity indicated by identity (`-p`)",
                            action="store_true",
                            default=False)
    seqop_args.add_argument("--multiple_alignment",
                            help='The FASTA input is also the multiple alignment file to be used.\n'
                                 'In this workflow, alignment with MAFFT is skipped and this file is used instead.',
                            action="store_true",
                            default=False)
    seqop_args.add_argument("--trim_align",
                            help="Flag indicating BMGE should be used to trim the non-conserved\n"
                                 "alignment positions in the multiple sequence alignment.\n",
                            action="store_true",
                            default=False)
    seqop_args.add_argument("-d", "--domain",
                            help="An HMM profile representing a specific domain.\n"
                                 "Domains will be excised from input sequences based on hmmsearch alignments.",
                            required=False, default=None)
    seqop_args.add_argument('-l', '--min_seq_length',
                            help='Minimal sequence length [DEFAULT = 100]',
                            required=False,
                            default=100,
                            type=int)
    seqop_args.add_argument('-m', '--molecule',
                            help='The type of input sequences:\n'
                                 'prot = Protein [DEFAULT]; dna = Nucleotide; rrna = rRNA',
                            default='prot',
                            choices=['prot', 'dna', 'rrna'])
    seqop_args.add_argument('-r', "--rfam_cm",
                            help="The covariance model of the RNA family being packaged.\n"
                                 "REQUIRED if molecule is rRNA!",
                            default=None)
    seqop_args.add_argument("-s", "--screen",
                            help="Keywords (taxonomic regular expressions) for including specific taxa in the tree.\n"
                                 "Example: to only include Bacteria and Archaea do `--screen Bacteria,Archaea`\n"
                                 "[ DEFAULT is no screen ]",
                            default="",
                            required=False)
    seqop_args.add_argument("-f", "--filter",
                            help="Keywords for removing specific taxa; the opposite of `--screen`.\n"
                                 "[ DEFAULT is no filter ]",
                            default="",
                            required=False)
    seqop_args.add_argument("-g", "--guarantee",
                            help="A FASTA file containing sequences that need to be included \n"
                                 "in the tree after all clustering and filtering",
                            default=None,
                            required=False)

    optopt = parser.add_argument_group("Optional options")
    optopt.add_argument("-b", "--bootstraps",
                        help="The number of bootstrap replicates RAxML should perform\n"
                             "[ DEFAULT = autoMR ]",
                        required=False, default="autoMR")
    optopt.add_argument("-e", "--raxml_model",
                        help="The evolutionary model for RAxML to use\n"
                             "[ Proteins = PROTGAMMAAUTO | Nucleotides =  GTRGAMMA ]",
                        required=False, default=None)
    optopt.add_argument("-o", "--output_dir",
                        help="Path to a directory for all outputs [ DEFAULT = ./ ]",
                        default="./", required=False)
    optopt.add_argument("-u", "--uc",
                        help="The USEARCH cluster format file produced from clustering reference sequences.\n"
                             "This can be used for selecting representative headers from identical sequences.",
                        required=False, default=None)
    optopt.add_argument("--taxa_lca",
                        help="Set taxonomy of representative sequences to LCA of cluster member's taxa.\n"
                             "[ --cluster or --uc REQUIRED ]",
                        default=False, required=False,
                        action="store_true")
    optopt.add_argument("--taxa_norm",
                        help="BETA: Perform taxonomic normalization on the provided sequences.\n"
                             "A comma-separated argument with the Rank (e.g. Phylum) and\n"
                             "number of representatives is required.\n")
    optopt.add_argument("--fast",
                        help="A flag indicating the tree should be built rapidly, using FastTree.",
                        default=False, required=False,
                        action="store_true")
    optopt.add_argument("--kind",
                        help="The broad classification of marker gene type, either "
                             "functional, phylogenetic, or phylogenetic_rRNA. [ DEFAULT = functional ]",
                        default="functional", required=False)

    miscellaneous_opts = parser.add_argument_group("Miscellaneous options")
    miscellaneous_opts.add_argument("-h", "--help",
                                    action="help",
                                    help="Show this help message and exit")
    miscellaneous_opts.add_argument('--overwrite', action='store_true', default=False,
                                    help='Overwrites previously processed output folders')
    miscellaneous_opts.add_argument('--pc', action='store_true', default=False,
                                    help='Prints the final commands to complete\n'
                                         'installation for a provided `code_name`')
    miscellaneous_opts.add_argument('--add_lineage', action='store_true', default=False,
                                    help='If the tax_ids file exists for the code_name,\n'
                                         'the third (lineage) column is appended then exits,\n'
                                         'leaving all other files.')
    miscellaneous_opts.add_argument("--headless", action="store_true", default=False,
                                    help="Do not require any user input.")
    miscellaneous_opts.add_argument("-T", "--num_threads",
                                    help="The number of threads for RAxML to use [ DEFAULT = 4 ]",
                                    required=False,
                                    default=str(4),
                                    type=str)
    miscellaneous_opts.add_argument('-v', '--verbose', action='store_true', default=False,
                                    help='Prints a more verbose runtime log')

    args = parser.parse_args()
    args.treesapp = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + os.sep
    if args.output_dir[0] != os.sep:
        # The user didn't provide a full path
        args.output_dir = os.getcwd() + os.sep + args.output_dir
    if args.output_dir[-1] != os.sep:
        args.output_dir += os.sep
    args.final_output_dir = args.output_dir + "TreeSAPP_files_%s" % args.code_name + os.sep

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


def generate_cm_data(args, unaligned_fasta):
    """
    Using the input unaligned FASTA file:
     1. align the sequences using cmalign against a reference Rfam covariance model to generate a Stockholm file
     2. use the Stockholm file (with secondary structure annotated) to build a covariance model
     3. align the sequences using cmalign against a reference Rfam covariance model to generate an aligned fasta (AFA)
    :param args: command-line arguments objects
    :param unaligned_fasta:
    :return:
    """
    logging.info("Running cmalign to build Stockholm file with secondary structure annotations... ")

    cmalign_base = [args.executables["cmalign"],
                    "--mxsize", str(3084),
                    "--informat", "FASTA",
                    "--cpu", str(args.num_threads)]
    # First, generate the stockholm file
    cmalign_sto = cmalign_base + ["-o", args.code_name + ".sto"]
    cmalign_sto += [args.rfam_cm, unaligned_fasta]

    stdout, cmalign_pro_returncode = launch_write_command(cmalign_sto)

    if cmalign_pro_returncode != 0:
        logging.error("cmalign did not complete successfully for:\n" + ' '.join(cmalign_sto) + "\n")
        sys.exit(13)

    logging.info("done.\n")
    logging.info("Running cmbuild... ")

    # Build the CM
    cmbuild_command = [args.executables["cmbuild"]]
    cmbuild_command += ["-n", args.code_name]
    cmbuild_command += [args.code_name + ".cm", args.code_name + ".sto"]

    stdout, cmbuild_pro_returncode = launch_write_command(cmbuild_command)

    if cmbuild_pro_returncode != 0:
        logging.error("cmbuild did not complete successfully for:\n" +
                      ' '.join(cmbuild_command) + "\n")
        sys.exit(13)
    os.rename(args.code_name + ".cm", args.final_output_dir + os.sep + args.code_name + ".cm")
    if os.path.isfile(args.final_output_dir + os.sep + args.code_name + ".sto"):
        logging.warning("Overwriting " + args.final_output_dir + os.sep + args.code_name + ".sto\n")
        os.remove(args.final_output_dir + os.sep + args.code_name + ".sto")
    shutil.move(args.code_name + ".sto", args.final_output_dir)

    logging.info("done.\n")
    logging.info("Running cmalign to build MSA... ")

    # Generate the aligned FASTA file which will be used to build the BLAST database and tree with RAxML
    aligned_fasta = args.code_name + ".fa"
    cmalign_afa = cmalign_base + ["--outformat", "Phylip"]
    cmalign_afa += ["-o", args.code_name + ".phy"]
    cmalign_afa += [args.rfam_cm, unaligned_fasta]

    stdout, cmalign_pro_returncode = launch_write_command(cmalign_afa)

    if cmalign_pro_returncode != 0:
        logging.error("cmalign did not complete successfully for:\n" + ' '.join(cmalign_afa) + "\n")
        sys.exit(13)

    # Convert the Phylip file to an aligned FASTA file for downstream use
    seq_dict = read_phylip_to_dict(args.code_name + ".phy")
    write_new_fasta(seq_dict, aligned_fasta)

    logging.info("done.\n")

    return aligned_fasta


def run_mafft(mafft_exe: str, fasta_in: str, fasta_out: str, num_threads):
    """
    Wrapper function for the MAFFT multiple sequence alignment tool.
    Runs MAFFT using `--auto` and checks if the output is empty.
    :param mafft_exe: Path to the executable for mafft
    :param fasta_in: An unaligned FASTA file
    :param fasta_out: The path to a file MAFFT will write aligned sequences to
    :param num_threads: Integer (or string) for the number of threads MAFFT can use
    :return:
    """
    mafft_align_command = [mafft_exe]
    mafft_align_command += ["--maxiterate", str(1000)]
    mafft_align_command += ["--thread", str(num_threads)]
    mafft_align_command.append("--auto")
    mafft_align_command += [fasta_in, '1>' + fasta_out]
    mafft_align_command += ["2>", "/dev/null"]

    stdout, mafft_proc_returncode = launch_write_command(mafft_align_command, False)

    if mafft_proc_returncode != 0:
        logging.error("Multiple sequence alignment using " + mafft_exe +
                      " did not complete successfully! Command used:\n" + ' '.join(mafft_align_command) + "\n")
        sys.exit(7)
    else:
        mfa = read_fasta_to_dict(fasta_out)
        if len(mfa.keys()) < 1:
            logging.error("MAFFT did not generate a proper FASTA file. " +
                          "Check the output by running:\n" + ' '.join(mafft_align_command) + "\n")
            sys.exit(7)

    return


def run_odseq(odseq_exe, fasta_in, outliers_fa, num_threads):
    odseq_command = [odseq_exe]
    odseq_command += ["-i", fasta_in]
    odseq_command += ["-f", "fasta"]
    odseq_command += ["-o", outliers_fa]
    odseq_command += ["-m", "linear"]
    odseq_command += ["--boot-rep", str(1000)]
    odseq_command += ["--threads", str(num_threads)]
    odseq_command += ["--score", str(5)]
    odseq_command.append("--full")

    stdout, odseq_proc_returncode = launch_write_command(odseq_command)

    if odseq_proc_returncode != 0:
        logging.error("Outlier detection using " + odseq_exe +
                      " did not complete successfully! Command used:\n" + ' '.join(odseq_command) + "\n")
        sys.exit(7)

    return


def create_new_ref_fasta(out_fasta, ref_seq_dict, dashes=False):
    """
    Writes a new FASTA file using a dictionary of ReferenceSequence class objects

    :param out_fasta: Name of the FASTA file to write to
    :param ref_seq_dict: Dictionary containing ReferenceSequence objects, numbers are keys
    :param dashes: Flag indicating whether hyphens should be retained in sequences
    :return:
    """
    out_fasta_handle = open(out_fasta, "w")
    num_seqs_written = 0

    for mltree_id in sorted(ref_seq_dict, key=int):
        ref_seq = ref_seq_dict[mltree_id]
        if dashes is False:
            sequence = re.sub('[-.]', '', ref_seq.sequence)
        else:
            # sequence = re.sub('\.', '', ref_seq.sequence)
            sequence = ref_seq.sequence
        out_fasta_handle.write(">" + ref_seq.short_id + "\n" + sequence + "\n")
        num_seqs_written += 1

    out_fasta_handle.close()

    if num_seqs_written == 0:
        logging.error("No sequences written to " + out_fasta + ".\n" +
                      "The headers in your input file are probably not accommodated in the regex patterns used. " +
                      "Function responsible: get_header_format. Please make an issue on the GitHub page.\n")
        sys.exit(5)

    return


def hmmsearch_input_references(args, fasta_replaced_file):
    """
    Function for searching a fasta file with an hmm profile
    :param args:
    :param fasta_replaced_file:
    :return:
    """
    # Find the name of the HMM. Use it to name the output file
    rp_marker = re.sub(".hmm", '', os.path.basename(args.domain))
    domtbl = args.output_dir + rp_marker + "_to_ORFs_domtbl.txt"

    # Basic hmmsearch command
    hmmsearch_command_base = [args.executables["hmmsearch"]]
    hmmsearch_command_base += ["--cpu", str(args.num_threads)]
    hmmsearch_command_base.append("--noali")
    # Customize the command for this input and HMM
    final_hmmsearch_command = hmmsearch_command_base + ["--domtblout", domtbl]
    final_hmmsearch_command += [args.domain, fasta_replaced_file]
    stdout, ret_code = launch_write_command(final_hmmsearch_command)

    # Check to ensure the job finished properly
    if ret_code != 0:
        logging.error("hmmsearch did not complete successfully!\n" + stdout + "\n" +
                      "Command used:\n" + ' '.join(final_hmmsearch_command) + "\n")
        sys.exit(13)

    return [domtbl]


def extract_hmm_matches(hmm_matches, fasta_dict, header_registry):
    """
    Function for slicing sequences guided by alignment co-ordinates.
    :param hmm_matches: Dictionary containing a list HmmMatch() objects as values for each 'marker' key
    :param fasta_dict: A dictionary with headers as keys and sequences as values
    :param header_registry: A list of Header() objects, each used to map various header formats to each other
    :return:
    """
    sys.stdout.write("Extracting the quality-controlled protein sequences... ")
    sys.stdout.flush()
    marker_gene_dict = dict()
    for marker in hmm_matches:
        if len(hmm_matches.keys()) > 1:
            logging.error("Number of markers found from HMM alignments is >1\n" +
                          "Does your HMM file contain more than 1 profile? TreeSAPP is unprepared for this.\n")
            sys.exit(13)
        if marker not in marker_gene_dict:
            marker_gene_dict[marker] = dict()

        for hmm_match in hmm_matches[marker]:
            # Now for the header format to be used in the bulk FASTA:
            # >contig_name|marker_gene|start_end
            sequence = ""
            bulk_header = ""
            for num in header_registry:
                if hmm_match.orf == header_registry[num].first_split[1:]:
                    query_names = header_registry[num]
                    sequence = fasta_dict[query_names.formatted]
                    if hmm_match.of > 1:
                        query_names.post_align = \
                            ' '.join([query_names.first_split, str(hmm_match.num) + '.' + str(hmm_match.of), re.sub(re.escape(query_names.first_split), '', query_names.original)])
                    else:
                        query_names.post_align = query_names.original
                    bulk_header = query_names.post_align
                    break
            if sequence:
                if bulk_header in marker_gene_dict[marker]:
                    logging.warning(bulk_header + " being overwritten by an alternative alignment!\n" +
                                    hmm_match.get_info())
                marker_gene_dict[marker][bulk_header] = sequence[hmm_match.start-1:hmm_match.end]
            else:
                logging.error("Unable to map " + hmm_match.orf + " to a sequence in the input FASTA.\n")
                sys.exit(13)
    sys.stdout.write("done.\n")
    return marker_gene_dict[marker]


def hmm_pile(hmm_matches):
    """
    Function to inspect the placement of query sequences on the reference HMM
    :param hmm_matches:
    :return:
    """
    hmm_bins = dict()
    window_size = 2

    for marker in hmm_matches:
        for hmm_match in hmm_matches[marker]:
            # Initialize the hmm_bins using the HMM profile length
            if not hmm_bins:
                i = 1
                hmm_length = int(hmm_match.hmm_len)
                while i < hmm_length:
                    if i+window_size-1 > hmm_length:
                        hmm_bins[(i, hmm_length)] = 0
                    else:
                        hmm_bins[(i, i+window_size-1)] = 0
                    i += window_size
            # Skip ahead to the HMM profile position where the query sequence began aligning
            for bin_start, bin_end in hmm_bins:
                if hmm_match.pstart <= bin_start and hmm_match.pend >= bin_end:
                    hmm_bins[(bin_start, bin_end)] += 1
                else:
                    pass
        low_coverage_start = 0
        low_coverage_stop = 0
        maximum_coverage = 0
        sys.stdout.write("Low coverage HMM windows (start-stop):\n")
        for window in sorted(hmm_bins.keys()):
            height = hmm_bins[window]
            if height > maximum_coverage:
                maximum_coverage = height
            if height < len(hmm_matches[marker])/2:
                begin, end = window
                if low_coverage_start == low_coverage_stop:
                    low_coverage_start = begin
                    low_coverage_stop = end
                else:
                    low_coverage_stop = end
            elif height > len(hmm_matches[marker])/2 and low_coverage_stop != 0:
                sys.stdout.write("\t" + str(low_coverage_start) + '-' + str(low_coverage_stop) + "\n")
                low_coverage_start = 0
                low_coverage_stop = 0
            else:
                pass
        if low_coverage_stop != low_coverage_start:
            sys.stdout.write("\t" + str(low_coverage_start) + "-end\n")
        sys.stdout.write("Maximum coverage = " + str(maximum_coverage) + " sequences\n")
    return


def regenerate_cluster_rep_swaps(args, cluster_dict, fasta_replace_dict):
    """
    Function to regenerate the swappers dictionary with the original headers as keys and
    the new header (swapped in the previous attempt based on USEARCH's uc file) as a value
    :param args: command-line arguments objects
    :param cluster_dict: Dictionary where keys are centroid headers and values are headers of identical sequences
    :param fasta_replace_dict: Immature (lacking sequences) dictionary with header information parsed from tax_ids file
    :return:
    """
    swappers = dict()
    if args.verbose:
        sys.stderr.write("Centroids with identical sequences in the unclustered input file:\n")
    for rep in sorted(cluster_dict):
        matched = False
        subs = cluster_dict[rep]
        # If its entry in cluster_dict == 0 then there were no identical
        # sequences and the header could not have been swapped
        if len(subs) >= 1:
            # If there is the possibility the header could have been swapped,
            # check if the header is in fasta_replace_dict
            for mltree_id in fasta_replace_dict:
                if matched:
                    break
                ref_seq = fasta_replace_dict[mltree_id]
                # If the accession from the tax_ids file is the same as the representative
                # this one has not been swapped for an identical sequence's header since it is in use
                if re.search(ref_seq.accession, rep):
                    if args.verbose:
                        sys.stderr.write("\tUnchanged: " + rep + "\n")
                        matched = True
                    break
                # The original representative is no longer in the reference sequences
                # so it was replaced, with this sequence...
                for candidate in subs:
                    if rep in swappers or matched:
                        break

                    # parse the accession from the header
                    header_format_re, header_db, header_molecule = get_header_format(candidate, args.code_name)
                    sequence_info = header_format_re.match(candidate)
                    if sequence_info:
                        candidate_acc = sequence_info.group(1)
                    else:
                        logging.error("Unable to handle header: " + candidate + "\n")
                        sys.exit(13)

                    # Now compare...
                    if candidate_acc == ref_seq.accession:
                        if args.verbose:
                            sys.stderr.write("\tChanged: " + candidate + "\n")
                        swappers[rep] = candidate
                        matched = True
                        break
            sys.stderr.flush()
    return swappers


def finalize_cluster_reps(cluster_dict: dict, refseq_objects, header_registry):
    """
        Transfer information from the cluster data (representative sequence, identity and cluster taxonomic LCA) to the
    dictionary of ReferenceSequence objects. The sequences not representing a cluster will have their `cluster_rep`
    flags remain *False* so as to not be analyzed further.

    :param cluster_dict:
    :param refseq_objects:
    :param header_registry: A list of Header() objects, each used to map various header formats to each other
    :return: Dictionary of ReferenceSequence objects with complete clustering information
    """
    for cluster_id in sorted(cluster_dict, key=int):
        cluster_info = cluster_dict[cluster_id]
        for treesapp_id in sorted(refseq_objects, key=int):
            if header_registry[treesapp_id].formatted == cluster_info.representative:
                refseq_objects[treesapp_id].cluster_rep_similarity = '*'
                refseq_objects[treesapp_id].cluster_rep = True
                refseq_objects[treesapp_id].cluster_lca = cluster_info.lca
    return refseq_objects


def present_cluster_rep_options(cluster_dict, refseq_objects, header_registry, important_seqs=None):
    """
    Present the headers of identical sequences to user for them to decide on representative header

    :param cluster_dict: dictionary from read_uc(uc_file)
    :param refseq_objects:
    :param header_registry: A list of Header() objects, each used to map various header formats to each other
    :param important_seqs: If --guarantee is provided, a dictionary mapping headers to seqs from format_read_fasta()
    :return:
    """
    if not important_seqs:
        important_seqs = dict()
    candidates = dict()
    for cluster_id in sorted(cluster_dict, key=int):
        cluster_info = cluster_dict[cluster_id]
        acc = 1
        candidates.clear()
        for num_id in sorted(refseq_objects, key=int):
            if header_registry[num_id].formatted == cluster_info.representative:
                refseq_objects[num_id].cluster_rep_similarity = '*'
                candidates[str(acc)] = refseq_objects[num_id]
                acc += 1
                break
        if acc != 2:
            raise AssertionError("Unable to find " + cluster_info.representative + " in ReferenceSequence objects!")

        if len(cluster_info.members) >= 1 and cluster_info.representative not in important_seqs.keys():
            for cluster_member_info in cluster_info.members:
                for treesapp_id in sorted(refseq_objects, key=int):
                    formatted_header = header_registry[treesapp_id].formatted
                    if formatted_header == cluster_member_info[0]:
                        refseq_objects[treesapp_id].cluster_rep_similarity = cluster_member_info[1]
                        candidates[str(acc)] = refseq_objects[treesapp_id]
                        acc += 1
                        break

            sys.stderr.write("Sequences in '" + cluster_info.lca + "' cluster:\n")
            for num in sorted(candidates.keys(), key=int):
                sys.stderr.write("\t" + num + ". ")
                sys.stderr.write('\t'.join([candidates[num].organism + " | " + candidates[num].accession + "\t",
                                            str(len(candidates[num].sequence)) + "bp",
                                            str(candidates[num].cluster_rep_similarity)]) + "\n")
            sys.stderr.flush()

            best = input("Number of the best representative? ")
            # Useful for testing - no need to pick which sequence name is best!
            # best = str(1)
            while best not in candidates.keys():
                best = input("Invalid number. Number of the best representative? ")
            candidates[best].cluster_rep = True
            candidates[best].cluster_lca = cluster_info.lca
        else:
            refseq_objects[num_id].cluster_rep = True
            refseq_objects[num_id].cluster_lca = cluster_info.lca

    return refseq_objects


def reformat_headers(header_dict):
    """
    Imitate format_read_fasta header name reformatting
    :param header_dict: Dictionary of old header : new header key : value pairs
    :return:
    """
    swappers = dict()

    for old, new in header_dict.items():
        swappers[reformat_string(old)] = reformat_string(new)
    return swappers


def get_sequence_info(code_name, fasta_dict, fasta_replace_dict, header_registry, swappers=None):
    """
    This function is used to find the accession ID and description of each sequence from the FASTA file

    :param code_name: code_name from the command-line parameters
    :param fasta_dict: a dictionary with headers as keys and sequences as values (returned by format_read_fasta)
    :param fasta_replace_dict:
    :param header_registry:
    :param swappers: A dictionary containing representative clusters (keys) and their constituents (values)
    :return: fasta_replace_dict with a complete ReferenceSequence() value for every mltree_id key
    """

    logging.info("Extracting information from headers for formatting purposes... ")
    fungene_gi_bad = re.compile("^>[0-9]+\s+coded_by=.+,organism=.+,definition=.+$")
    swapped_headers = []
    if len(fasta_replace_dict.keys()) > 0:
        for mltree_id in sorted(fasta_replace_dict):
            ref_seq = fasta_replace_dict[mltree_id]
            ref_seq.short_id = mltree_id + '_' + code_name
            for header in fasta_dict:
                # Find the matching header in the header_registry
                original_header = header_registry[mltree_id].original
                header_format_re, header_db, header_molecule = get_header_format(original_header, code_name)
                sequence_info = header_format_re.match(original_header)
                _, fasta_header_organism, _, _, _ = return_sequence_info_groups(sequence_info, header_db, header)
                if re.search(ref_seq.accession, header):
                    if re.search(reformat_string(ref_seq.organism), reformat_string(fasta_header_organism)):
                        ref_seq.sequence = fasta_dict[header]
                    else:
                        logging.warning("Accession '" + ref_seq.accession + "' matches, organism differs:\n" +
                                        "'" + ref_seq.organism + "' versus '" + fasta_header_organism + "'\n")
            if not ref_seq.sequence:
                # TODO: test this case (uc file provided, both fresh attempt and re-attempt)
                # Ensure the header isn't a value within the swappers dictionary
                for swapped in swappers.keys():
                    header = swappers[swapped]
                    original_header = ""
                    # Find the original header of the swapped header
                    for num in header_registry:
                        if header_registry[num].first_split[1:] == header:
                            original_header = header_registry[num].original
                        elif re.search(header_registry[num].first_split[1:], header):
                            original_header = header_registry[num].original
                        else:
                            pass
                    if not original_header:
                        logging.error("Unable to find the original header for " + header + "\n")
                        sys.exit(13)
                    if re.search(ref_seq.accession, header) and re.search(ref_seq.organism, original_header):
                        # It is and therefore the header was swapped last run
                        ref_seq.sequence = fasta_dict[swapped]
                        break
                if not ref_seq.sequence:
                    # Unable to find sequence in swappers too
                    logging.error("Unable to find header for " + ref_seq.accession)
                    sys.exit(13)

    else:  # if fasta_replace_dict needs to be populated, this is a new run
        for header in sorted(fasta_dict.keys()):
            if fungene_gi_bad.match(header):
                logging.warning("Input sequences use 'GIs' which are obsolete and may be non-unique. " +
                                "For everyone's sanity, please download sequences with the `accno` instead.\n")

            # Try to find the original header in header_registry
            original_header = ""
            mltree_id = ""
            for num in header_registry:
                if header == header_registry[num].formatted:
                    original_header = header_registry[num].original
                    mltree_id = str(num)
                    break
            ref_seq = ReferenceSequence()
            ref_seq.sequence = fasta_dict[header]

            if swappers and header in swappers.keys():
                header = swappers[header]
                swapped_headers.append(header)
            if original_header and mltree_id:
                pass
            else:
                logging.error("Unable to find the header:\n\t" + header +
                              "\nin header_map (constructed from either the input FASTA or .uc file).\n" +
                              "There is a chance this is due to the FASTA file and .uc being generated separately.\n")
                sys.exit(13)
            header_format_re, header_db, header_molecule = get_header_format(original_header, code_name)
            sequence_info = header_format_re.match(original_header)
            ref_seq.accession,\
                ref_seq.organism,\
                ref_seq.locus,\
                ref_seq.description,\
                ref_seq.lineage = return_sequence_info_groups(sequence_info, header_db, original_header)

            ref_seq.short_id = mltree_id + '_' + code_name
            fasta_replace_dict[mltree_id] = ref_seq

        if swappers and len(swapped_headers) != len(swappers):
            logging.error("Some headers that were meant to be replaced could not be compared!\n")
            for header in swappers.keys():
                if header not in swapped_headers:
                    sys.stdout.write(header + "\n")
            sys.exit(13)

    logging.info("done.\n")

    return fasta_replace_dict


def screen_filter_taxa(args, fasta_replace_dict):
    if args.screen == "" and args.filter == "":
        return fasta_replace_dict
    else:
        if args.screen:
            screen_terms = args.screen.split(',')
        else:
            screen_terms = ''
        if args.filter:
            filter_terms = args.filter.split(',')
        else:
            filter_terms = ''

    purified_fasta_dict = dict()
    num_filtered = 0
    num_screened = 0

    for mltree_id in fasta_replace_dict:
        screen_pass = False
        filter_pass = True
        ref_seq = fasta_replace_dict[mltree_id]
        # Screen
        if len(screen_terms) > 0:
            for term in screen_terms:
                # If any term is found in the lineage, it will pass... unless it fails the filter
                if re.search(term, ref_seq.lineage):
                    screen_pass = True
                    break
        else:
            screen_pass = True
        # Filter
        if len(filter_terms) > 0:
            for term in filter_terms:
                if re.search(term, ref_seq.lineage):
                    filter_pass = False

        if filter_pass and screen_pass:
            purified_fasta_dict[mltree_id] = ref_seq
        else:
            if screen_pass is False:
                num_screened += 1
            if filter_pass is False:
                num_filtered += 1

    logging.debug('\t' + str(num_screened) + " sequences removed after failing screen.\n" +
                  '\t' + str(num_filtered) + " sequences removed after failing filter.\n" +
                  '\t' + str(len(purified_fasta_dict.keys())) + " sequences retained for building tree.\n")

    return purified_fasta_dict


def remove_duplicate_records(fasta_record_objects):
    nr_record_dict = dict()
    accessions = dict()
    dups = False
    for treesapp_id in fasta_record_objects:
        ref_seq = fasta_record_objects[treesapp_id]
        if ref_seq.accession not in accessions:
            accessions[ref_seq.accession] = 0
            nr_record_dict[treesapp_id] = ref_seq
        else:
            dups = True
        accessions[ref_seq.accession] += 1
    if dups:
        logging.warning("Redundant accessions have been detected in your input FASTA.\n" +
                        "The duplicates have been removed leaving a single copy for further analysis.\n" +
                        "Please view the log file for the list of redundant accessions and their copy numbers.\n")
        msg = "Redundant accessions found and copies:\n"
        for acc in accessions:
            if accessions[acc] > 1:
                msg += "\t" + acc + "\t" + str(accessions[acc]) + "\n"
        logging.debug(msg)
    return nr_record_dict


def order_dict_by_lineage(fasta_object_dict):
    """
    Re-order the fasta_record_objects by their lineages (not phylogenetic, just alphabetical sort)
    Remove the cluster members since they will no longer be used

    :param fasta_object_dict: A dictionary mapping `treesapp_id`s (integers) to ReferenceSequence objects
    :return: An ordered, filtered version of the input dictionary
    """
    # Create a new dictionary with lineages as keys
    lineage_dict = dict()
    sorted_lineage_dict = dict()
    accessions = list()
    for treesapp_id in fasta_object_dict:
        ref_seq = fasta_object_dict[treesapp_id]
        if ref_seq.accession in accessions:
            logging.error("Uh oh... duplicate accession identifiers '" + ref_seq.accession + "' found!\n" +
                          "TreeSAPP should have removed these by now. " +
                          "Please alert the developers so they can cobble a fix together.\n")
            sys.exit(13)
        else:
            accessions.append(ref_seq.accession)
        # Skip the redundant sequences that are not cluster representatives
        if not ref_seq.cluster_rep:
            continue
        if ref_seq.lineage not in lineage_dict.keys():
            # Values of the new dictionary are lists of ReferenceSequence instances
            lineage_dict[ref_seq.lineage] = list()
        lineage_dict[ref_seq.lineage].append(ref_seq)

    # Now re-write the fasta_object_dict, but the numeric keys are now sorted by lineage
    #  AND it doesn't contain redundant fasta objects
    num_key = 1
    for lineage in sorted(lineage_dict.keys(), key=str):
        for ref_seq in lineage_dict[lineage]:
            if ref_seq.cluster_rep:
                # Replace the treesapp_id object
                code = '_'.join(ref_seq.short_id.split('_')[1:])
                ref_seq.short_id = str(num_key) + '_' + code
                sorted_lineage_dict[str(num_key)] = ref_seq
                num_key += 1

    return sorted_lineage_dict


def threshold(lst, confidence="low"):
    """

    :param lst:
    :param confidence:
    :return:
    """
    if confidence == "low":
        # Majority calculation
        index = round(len(lst)*0.51)-1
    elif confidence == "medium":
        # >=75% of the list is reported
        index = round(len(lst)*0.75)-1
    else:
        # confidence is "high" and >=90% of the list is reported
        index = round(len(lst)*0.9)-1
    return sorted(lst, reverse=True)[index]


def estimate_taxonomic_redundancy(args, reference_dict):
    """

    :param args:
    :param reference_dict:
    :return:
    """
    # TODO: Factor proximity of leaves in the tree into this measure
    # For instance, if the two or so species of the same genus are in the tree,
    # are they also beside each other in the same clade or are they located in different clusters?
    lowest_reliable_rank = "Strain"
    rank_depth_map = {1: "Kingdoms", 2: "Phyla", 3: "Classes", 4: "Orders", 5: "Families", 6: "Genera", 7: "Species"}
    taxa_counts = dict()
    for depth in rank_depth_map:
        name = rank_depth_map[depth]
        taxa_counts[name] = dict()
    for mltree_id_key in sorted(reference_dict.keys(), key=int):
        lineage = reference_dict[mltree_id_key].lineage
        position = 1
        taxa = lineage.split('; ')
        while position < len(taxa) and position < 8:
            if taxa[position] not in taxa_counts[rank_depth_map[position]]:
                taxa_counts[rank_depth_map[position]][taxa[position]] = 0
            taxa_counts[rank_depth_map[position]][taxa[position]] += 1
            position += 1
    for depth in rank_depth_map:
        rank = rank_depth_map[depth]
        redundancy = list()
        for taxon in taxa_counts[rank]:
            redundancy.append(taxa_counts[rank][taxon])
        if threshold(redundancy, "medium") == 1:
            lowest_reliable_rank = rank
            break

    logging.info("Lowest reliable rank for taxonomic classification is: " + lowest_reliable_rank + "\n")

    return lowest_reliable_rank


def summarize_reference_taxa(args, reference_dict):
    """
        Function for enumerating the representation of each taxonomic rank within the finalized reference sequences
    :param args:
    :param reference_dict:
    :return: A formatted, human-readable string stating the number of unique taxa at each rank
    """
    taxonomic_summary_string = ""
    # Not really interested in Cellular Organisms or Strains.
    rank_depth_map = {1: "Kingdoms", 2: "Phyla", 3: "Classes", 4: "Orders", 5: "Families", 6: "Genera", 7: "Species"}
    taxa_counts = dict()
    unclassifieds = 0

    for depth in rank_depth_map:
        name = rank_depth_map[depth]
        taxa_counts[name] = set()
    for num_id in sorted(reference_dict.keys(), key=int):
        if args.taxa_lca and reference_dict[num_id].cluster_lca:
            lineage = reference_dict[num_id].cluster_lca
        else:
            lineage = reference_dict[num_id].lineage

        if re.search("unclassified", lineage, re.IGNORECASE):
            unclassifieds += 1

        position = 1
        taxa = lineage.split('; ')
        while position < len(taxa) and position < 8:
            taxa_counts[rank_depth_map[position]].add(taxa[position])
            position += 1

    taxonomic_summary_string += "Number of unique lineages:\n"
    for depth in rank_depth_map:
        rank = rank_depth_map[depth]
        buffer = " "
        while len(rank) + len(str(len(taxa_counts[rank]))) + len(buffer) < 12:
            buffer += ' '
        taxonomic_summary_string += "\t" + rank + buffer + str(len(taxa_counts[rank])) + "\n"
    # Report number of "Unclassified" lineages
    taxonomic_summary_string += "Unclassified lineages account for " +\
                                str(unclassifieds) + '/' + str(len(reference_dict.keys())) + ' (' +\
                                str(round(float(unclassifieds*100)/len(reference_dict.keys()), 1)) + "%) references.\n"

    return taxonomic_summary_string


def write_tax_ids(args, fasta_replace_dict, tax_ids_file):
    """
    Write the number, organism and accession ID, if possible
    :param args: command-line arguments object
    :param fasta_replace_dict: Dictionary mapping numbers (internal treesapp identifiers) to ReferenceSequence objects
    :param tax_ids_file: The name of the output file
    :return: Nothing
    """

    tree_taxa_string = ""
    warning_string = ""
    no_lineage = list()

    for mltree_id_key in sorted(fasta_replace_dict.keys(), key=int):
        # Definitely will not uphold phylogenetic relationships but at least sequences
        # will be in the right neighbourhood rather than ordered by their position in the FASTA file
        reference_sequence = fasta_replace_dict[mltree_id_key]
        if args.taxa_lca:
            lineage = reference_sequence.cluster_lca
        else:
            lineage = reference_sequence.lineage
        if not lineage:
            no_lineage.append(reference_sequence.accession)
            lineage = ''

        tree_taxa_string += "\t".join([str(mltree_id_key),
                                      reference_sequence.organism + " | " + reference_sequence.accession,
                                       lineage]) + "\n"

    # Write the tree_taxa_string to the tax_ids file
    tree_tax_list_handle = open(tax_ids_file, "w")
    tree_tax_list_handle.write(tree_taxa_string)
    tree_tax_list_handle.close()

    if len(no_lineage) > 0:
        warning_string += str(len(no_lineage)) + " reference sequences did not have an associated lineage!\n\t"
        warning_string += "\n\t".join(no_lineage)

    return warning_string


def read_tax_ids(tree_taxa_list):
    """
    Reads the taxonomy and accession ID affiliated with each sequence number.
    This information is used to avoid horrible manual work if the pipeline is ran multiple times
    :param tree_taxa_list: The name of the tax_ids file to read
    :return:
    """
    try:
        tree_tax_list_handle = open(tree_taxa_list, 'r')
    except IOError:
        logging.error("Unable to open taxa list file '" + tree_taxa_list + "' for reading!\n")
        sys.exit(13)
    fasta_replace_dict = dict()
    line = tree_tax_list_handle.readline()
    while line:
        fields = line.strip().split("\t")
        if len(fields) == 3:
            mltree_id_key, seq_info, lineage = fields
        else:
            mltree_id_key, seq_info = fields
            lineage = ""
        ref_seq = ReferenceSequence()
        try:
            ref_seq.organism = seq_info.split(" | ")[0]
            ref_seq.accession = seq_info.split(" | ")[1]
            ref_seq.lineage = lineage
        except IndexError:
            ref_seq.organism = seq_info
        fasta_replace_dict[mltree_id_key] = ref_seq
        line = tree_tax_list_handle.readline()
    tree_tax_list_handle.close()

    return fasta_replace_dict


def swap_tree_names(tree_to_swap, final_mltree, code_name):
    original_tree = open(tree_to_swap, 'r')
    raxml_tree = open(final_mltree, 'w')

    tree = original_tree.readlines()
    original_tree.close()
    if len(tree) > 1:
        logging.error(">1 line contained in RAxML tree " + tree_to_swap + "\n")
        sys.exit(13)

    new_tree = re.sub('_' + re.escape(code_name), '', str(tree[0]))
    raxml_tree.write(new_tree)

    raxml_tree.close()
    return


def find_model_used(raxml_info_file):
    model_statement_re = re.compile(r".* model: ([A-Z]+) likelihood.*")
    model = ""
    command_line = ""
    with open(raxml_info_file) as raxml_info:
        for line in raxml_info:
            if model_statement_re.search(line):
                model = model_statement_re.search(line).group(1)
                break
            elif re.match('^.*/raxml.*-m ([A-Z]+)$', line):
                command_line = line
            else:
                pass
    if model == "":
        if command_line == "":
            logging.warning("Unable to parse model used from " + raxml_info_file + "!\n")
        else:
            model = re.match('^.*/raxml.*-m ([A-Z]+)$', command_line).group(1)
    return model


def update_build_parameters(param_file, marker_package: MarkerBuild):
    """
    Function to update the data/tree_data/ref_build_parameters.tsv file with information on this new reference sequence
    Format of file is:
     "\t".join(["name","code","molecule","sub_model","cluster_identity","ref_sequences","tree-tool","poly-params",
     "lowest_reliable_rank","last_updated","description"])

    :param param_file: Path to the ref_build_parameters.tsv file used by TreeSAPP for storing refpkg metadata
    :param marker_package: A MarkerBuild instance
    :return: 
    """
    try:
        params = open(param_file, 'a')
    except IOError:
        logging.error("Unable to open " + param_file + "for appending.\n")
        sys.exit(13)

    marker_package.update = strftime("%d_%b_%Y", gmtime())

    if marker_package.molecule == "prot":
        marker_package.model = "PROTGAMMA" + marker_package.model

    build_list = [marker_package.cog, marker_package.denominator, marker_package.molecule, marker_package.model,
                  marker_package.kind, str(marker_package.pid), str(marker_package.num_reps), marker_package.tree_tool,
                  ','.join([str(param) for param in marker_package.pfit]),
                  marker_package.lowest_confident_rank, marker_package.update, marker_package.description]
    params.write("\t".join(build_list) + "\n")

    return


def terminal_commands(final_output_folder, code_name):
    logging.info("\nTo integrate these data for use in TreeSAPP, the following steps must be performed:\n" +
                 "1. Include properly formatted 'denominator' codes " +
                 "in data/tree_data/cog_list.tsv and data/tree_data/ref_build_parameters.tsv\n" +
                 "2. $ cp " + final_output_folder + os.sep + "tax_ids_%s.txt" % code_name + " data/tree_data/\n" +
                 "3. $ cp " + final_output_folder + os.sep + code_name + "_tree.txt data/tree_data/\n" +
                 "4. $ cp " + final_output_folder + os.sep + code_name + ".hmm data/hmm_data/\n" +
                 "5. $ cp " + final_output_folder + os.sep + code_name + ".fa data/alignment_data/\n")
    return


def construct_tree(args, multiple_alignment_file, tree_output_dir):
    """
    Wrapper script for generating phylogenetic trees with either RAxML or FastTree from a multiple alignment

    :param args:
    :param multiple_alignment_file:
    :param tree_output_dir:
    :return:
    """

    # Decide on the command to build the tree, make some directories and files when necessary
    final_mltree = args.final_output_dir + args.code_name + "_tree.txt"
    if args.fast:
        tree_build_cmd = [args.executables["FastTree"]]
        if args.molecule == "rrna" or args.molecule == "dna":
            tree_build_cmd += ["-nt", "-gtr"]
        else:
            tree_build_cmd += ["-lg", "-wag"]
        tree_build_cmd += ["-out", final_mltree]
        tree_build_cmd.append(multiple_alignment_file)
        tree_builder = "FastTree"
    else:
        tree_build_cmd = [args.executables["raxmlHPC"]]
        tree_build_cmd += ["-f", "a"]
        tree_build_cmd += ["-p", "12345"]
        tree_build_cmd += ["-x", "12345"]
        tree_build_cmd += ["-#", args.bootstraps]
        tree_build_cmd += ["-s", multiple_alignment_file]
        tree_build_cmd += ["-n", args.code_name]
        tree_build_cmd += ["-w", tree_output_dir]
        tree_build_cmd += ["-T", args.num_threads]

        if args.raxml_model:
            tree_build_cmd += ["-m", args.raxml_model]
        elif args.molecule == "prot":
            tree_build_cmd += ["-m", "PROTGAMMAAUTO"]
        elif args.molecule == "rrna" or args.molecule == "dna":
            tree_build_cmd += ["-m", "GTRGAMMA"]
        else:
            logging.error("A substitution model could not be specified with the 'molecule' argument: " + args.molecule)
            sys.exit(13)
        tree_builder = "RAxML"

    # Ensure the tree from a previous run isn't going to be over-written
    if not os.path.exists(tree_output_dir):
        os.makedirs(tree_output_dir)
    else:
        logging.error(tree_output_dir + " already exists from a previous run! " +
                      "Please delete or rename it and try again.\n")
        sys.exit(13)

    if args.fast:
        logging.info("Building Approximately-Maximum-Likelihood tree with FastTree... ")
        stdout, returncode = launch_write_command(tree_build_cmd, True)
        with open(tree_output_dir + os.sep + "FastTree_info." + args.code_name, 'w') as fast_info:
            fast_info.write(stdout + "\n")
    else:
        logging.info("Building phylogenetic tree with RAxML... ")
        stdout, returncode = launch_write_command(tree_build_cmd, False)
    logging.info("done.\n")

    if returncode != 0:
        logging.error(tree_builder + " did not complete successfully! " +
                      "Look in " + tree_output_dir + os.sep +
                      tree_builder + "_info." + args.code_name + " for an error message.\n" +
                      tree_builder + " command used:\n" + ' '.join(tree_build_cmd) + "\n")
        sys.exit(13)

    if not args.fast:
        raw_newick_tree = "%s/RAxML_bestTree.%s" % (tree_output_dir, args.code_name)
        bootstrap_tree = tree_output_dir + os.sep + "RAxML_bipartitionsBranchLabels." + args.code_name
        bootstrap_nameswap = args.final_output_dir + args.code_name + "_bipartitions.txt"
        shutil.copy(multiple_alignment_file, tree_output_dir)
        os.remove(multiple_alignment_file)
        swap_tree_names(raw_newick_tree, final_mltree, args.code_name)
        swap_tree_names(bootstrap_tree, bootstrap_nameswap, args.code_name)

    return tree_builder


def update_tax_ids_with_lineage(args, tree_taxa_list):
    tax_ids_file = args.treesapp + os.sep + "data" + os.sep + "tree_data" + os.sep + "tax_ids_%s.txt" % args.code_name
    if not os.path.exists(tax_ids_file):
        logging.error("Unable to find " + tax_ids_file + "!\n")
        raise FileNotFoundError
    else:
        fasta_replace_dict = read_tax_ids(tax_ids_file)
        # Determine how many sequences already have lineage information:
        lineage_info_complete = 0
        for mltree_id_key in fasta_replace_dict:
            ref_seq = fasta_replace_dict[mltree_id_key]
            if ref_seq.lineage:
                lineage_info_complete += 1
        # There are some that are already complete. Should they be over-written?
        if lineage_info_complete >= 1:
            if sys.version_info > (2, 9):
                overwrite_lineages = input(tree_taxa_list + " contains some sequences with complete lineages. "
                                                            "Should they be over-written? [y|n] ")
                while overwrite_lineages != "y" and overwrite_lineages != "n":
                    overwrite_lineages = input("Incorrect response. Please input either 'y' or 'n'. ")
            else:
                overwrite_lineages = raw_input(tree_taxa_list + " contains some sequences with complete lineages."
                                                                "Should they be over-written? [y|n] ")
                while overwrite_lineages != "y" and overwrite_lineages != "n":
                    overwrite_lineages = raw_input("Incorrect response. Please input either 'y' or 'n'. ")
            if overwrite_lineages == 'y':
                ref_seq_dict = dict()
                for mltree_id_key in fasta_replace_dict:
                    ref_seq = fasta_replace_dict[mltree_id_key]
                    if ref_seq.lineage:
                        ref_seq.lineage = ""
                    ref_seq_dict[mltree_id_key] = ref_seq
        # write_tax_ids(args, fasta_replace_dict, tax_ids_file, args.molecule)
    return


def finalize_ref_seq_lineages(fasta_record_objects, accession_lineages):
    """
    Adds lineage information from accession_lineages to fasta_record_objects

    :param fasta_record_objects: dict() indexed by TreeSAPP numeric identifiers mapped to ReferenceSequence instances
    :param accession_lineages: a dictionary mapping {accession: lineage}
    :return:
    """
    for treesapp_id in fasta_record_objects:
        ref_seq = fasta_record_objects[treesapp_id]
        if not ref_seq.lineage:
            try:
                lineage = accession_lineages[ref_seq.accession]
            except KeyError:
                logging.error("Lineage information was not retrieved for " + ref_seq.accession + "!\n" +
                              "Please remove the output directory and restart.\n")
                sys.exit(13)
            # Add the species designation since it is often not included in the sequence record's lineage
            ref_seq.lineage = lineage
        if not ref_seq.organism and ref_seq.lineage:
            ref_seq.organism = ref_seq.lineage.split("; ")[-1]
        else:
            pass
    return fasta_record_objects


def main():
    # TODO: Record each external software command and version in log
    ##
    # STAGE 0: PARAMETERIZE - retrieve command line arguments, query user about settings if necessary
    ##
    args = get_arguments()

    if os.path.isdir(args.output_dir) and args.overwrite:
        shutil.rmtree(args.output_dir)
        os.mkdir(args.output_dir)
    os.makedirs(args.output_dir, exist_ok=True)

    log_file_name = args.output_dir + os.sep + "create_" + args.code_name + "_TreeSAPP_log.txt"
    prep_logging(log_file_name, args.verbose)
    logging.info("\n##\t\t\tCreating TreeSAPP reference package for '" + args.code_name + "' \t\t\t##\n")
    logging.info("Command used:\n" + ' '.join(sys.argv) + "\n")

    args = find_executables(args)

    if args.pc:
        terminal_commands(args.final_output_dir, args.code_name)
        sys.exit(0)

    # Names of files to be created
    tree_output_dir = args.output_dir + args.code_name + "_phy_files" + os.sep
    accession_map_file = args.output_dir + os.sep + "accession_id_lineage_map.tsv"
    hmm_purified_fasta = args.output_dir + args.code_name + "_hmm_purified.fasta"
    filtered_fasta_name = args.output_dir + '.'.join(os.path.basename(args.fasta_input).split('.')[:-1]) + "_filtered.fa"
    uclust_prefix = args.output_dir + '.'.join(os.path.basename(filtered_fasta_name).split('.')[:-1]) + "_uclust" + args.identity
    clustered_fasta = uclust_prefix + ".fa"
    clustered_uc = uclust_prefix + ".uc"
    od_input = args.output_dir + "od_input.fasta"
    unaln_ref_fasta = args.output_dir + args.code_name + "_ref.fa"  # FASTA file of unaligned reference sequences
    phylip_file = args.output_dir + args.code_name + ".phy"  # Used for building the phylogenetic tree with RAxML

    # Gather all the final TreeSAPP reference files
    ref_pkg = ReferencePackage()
    ref_pkg.gather_package_files(args.code_name, args.final_output_dir, "flat")

    # Create a new MarkerBuild instance to hold all relevant information for recording in ref_build_parameters.tsv
    marker_package = MarkerBuild()
    marker_package.pid = args.identity
    marker_package.cog = args.code_name
    marker_package.molecule = args.molecule
    marker_package.kind = args.kind
    marker_package.denominator = "Z1111"

    # TODO: Restore this functionality
    # if args.add_lineage:
    #     update_tax_ids_with_lineage(args, tree_taxa_list)
    #     terminal_commands(args.final_output_dir, code_name)
    #     sys.exit(0)

    if not os.path.exists(args.final_output_dir):
        try:
            os.makedirs(args.final_output_dir, exist_ok=False)
        except OSError:
            logging.warning("Making all directories in path " + args.final_output_dir + "\n")
            os.makedirs(args.final_output_dir, exist_ok=True)
    else:
        logging.warning("Output directory already exists. " +
                        "You have 10 seconds to hit Ctrl-C before previous outputs will be overwritten.\n")
        sleep(10)
        # Comment out to skip tree-building
        if os.path.exists(tree_output_dir):
            shutil.rmtree(tree_output_dir)

    ##
    # STAGE 2: FILTER - begin filtering sequences by homology and taxonomy
    ##
    if args.domain:
        if os.path.isfile(hmm_purified_fasta):
            logging.info("Using " + hmm_purified_fasta + " from a previous attempt.\n")
        else:
            logging.info("Searching for domain sequences... ")
            hmm_domtbl_files = hmmsearch_input_references(args, args.fasta_input)
            logging.info("done.\n")
            hmm_matches = parse_domain_tables(args, hmm_domtbl_files)
            # If we're screening a massive fasta file, we don't want to read every sequence - just those with hits
            # TODO: Implement a screening procedure in _fasta_reader._read_format_fasta()
            fasta_dict = format_read_fasta(args.fasta_input, args.molecule, args.output_dir)
            header_registry = register_headers(get_headers(args.fasta_input))
            marker_gene_dict = extract_hmm_matches(hmm_matches, fasta_dict, header_registry)
            write_new_fasta(marker_gene_dict, hmm_purified_fasta)
            summarize_fasta_sequences(hmm_purified_fasta)
            hmm_pile(hmm_matches)

        fasta_dict = format_read_fasta(hmm_purified_fasta, args.molecule, args.output_dir)
        header_registry = register_headers(get_headers(hmm_purified_fasta))
        # Point all future operations to the HMM purified FASTA file as the original input
        args.fasta_input = hmm_purified_fasta
    else:
        fasta_dict = format_read_fasta(args.fasta_input, args.molecule, args.output_dir)
        header_registry = register_headers(get_headers(args.fasta_input))
    unprocessed_fasta_dict = read_fasta_to_dict(args.fasta_input)

    ##
    # Synchronize records between fasta_dict and header_registry (e.g. short ones may be removed by format_read_fasta())
    ##
    if len(fasta_dict.keys()) != len(header_registry):
        sync_header_registry = dict()
        excluded_headers = list()
        for num_id in header_registry:
            if header_registry[num_id].formatted not in fasta_dict:
                excluded_headers.append(header_registry[num_id].original)
            else:
                sync_header_registry[num_id] = header_registry[num_id]
        logging.warning("The following sequences have been excluded from downstream analyses:\n" +
                        "\n".join(excluded_headers))
        header_registry = sync_header_registry

    ##
    # If there are sequences that needs to be guaranteed to be included,
    #  add them now as its easier to work with more sequences than repeat everything
    ##
    if args.guarantee:
        if not os.path.isfile(args.guarantee):
            logging.error("File '" + args.guarantee + "' does not exist!\n")
            sys.exit(13)
        important_seqs = format_read_fasta(args.guarantee, args.molecule, args.output_dir)
        important_headers = register_headers(get_headers(args.guarantee))
        fasta_dict.update(important_seqs)
        acc = max([int(x) for x in header_registry.keys()])
        for num_id in sorted(important_headers, key=int):
            acc += 1
            header_registry[str(acc)] = important_headers[num_id]
    else:
        important_seqs = None

    ##
    # Determine the format of each sequence name (header) and pull important info (e.g. accession, lineage)
    ##
    fasta_record_objects = get_header_info(header_registry, args.code_name)

    ##
    # Read lineages corresponding to accessions for each sequence if available, otherwise download them
    ##
    query_accession_list, num_lineages_provided = build_entrez_queries(fasta_record_objects)
    if os.path.isfile(accession_map_file):
        logging.info("Reading cached lineages in '" + accession_map_file + "'... ")
        accession_lineage_map = read_accession_taxa_map(accession_map_file)
        logging.info("done.\n")
    else:
        accession_lineage_map, all_accessions = get_multiple_lineages(query_accession_list,
                                                                      args.molecule)
        # Download lineages separately for those accessions that failed
        # Map proper accession to lineage from the tuple keys (accession, accession.version)
        #  in accession_lineage_map returned by get_multiple_lineages.
        fasta_record_objects, accession_lineage_map = verify_lineage_information(accession_lineage_map, all_accessions,
                                                                                 fasta_record_objects, num_lineages_provided,
                                                                                 args.molecule)
        write_accession_lineage_map(accession_map_file, accession_lineage_map)
    # Add lineage information to the ReferenceSequence() objects in fasta_record_objects if not contained
    finalize_ref_seq_lineages(fasta_record_objects, accession_lineage_map)

    ##
    # Remove the sequences failing 'filter' and/or only retain the sequences in 'screen'
    ##
    if args.add_lineage:
        if args.screen or args.filter:
            logging.warning("Skipping taxonomic filtering and screening in `--add_lineage` mode.\n")
    else:
        fasta_record_objects = screen_filter_taxa(args, fasta_record_objects)

    if len(fasta_record_objects.keys()) < 2:
        logging.error(str(len(fasta_record_objects)) + " sequences post-homology + taxonomy filtering\n")
        sys.exit(11)

    fasta_record_objects = remove_duplicate_records(fasta_record_objects)

    # Add the respective protein or nucleotide sequence string to each ReferenceSequence object
    for num_id in fasta_record_objects:
        refseq_object = fasta_record_objects[num_id]
        treesapp_id = refseq_object.short_id[1:].split('_')[0]
        refseq_object.sequence = fasta_dict[header_registry[treesapp_id].formatted]
    # Write a new FASTA file containing the sequences that passed the homology and taxonomy filters
    filtered_fasta_dict = dict()
    for num_id in fasta_record_objects:
        refseq_object = fasta_record_objects[num_id]
        formatted_header = header_registry[num_id].formatted
        filtered_fasta_dict[formatted_header] = refseq_object.sequence
    write_new_fasta(filtered_fasta_dict, filtered_fasta_name)

    ##
    # Optionally cluster the input sequences using USEARCH at the specified identity
    ##
    if args.cluster:
        cluster_sequences(args.executables["usearch"], filtered_fasta_name, uclust_prefix, args.identity)
        args.uc = clustered_uc

    ##
    # Read the uc file if present
    ##
    if args.uc:
        cluster_dict = read_uc(args.uc)

        # Ensure the headers in cluster_dict have been reformatted if UC file was not generated internally
        if not args.cluster:
            members = list()
            for num_id in cluster_dict:
                cluster = cluster_dict[num_id]
                cluster.representative = reformat_string(cluster.representative)
                for member in cluster.members:
                    header, identity = member
                    members.append([reformat_string(header), identity])
                cluster.members = members
                members.clear()
        logging.debug("\t" + str(len(cluster_dict.keys())) + " sequence clusters\n")
        ##
        # Calculate LCA of each cluster to represent the taxonomy of the representative sequence
        ##
        lineages = list()
        for cluster_id in sorted(cluster_dict, key=int):
            members = [cluster_dict[cluster_id].representative]
            # format of member list is: [header, identity, member_seq_length/representative_seq_length]
            members += [member[0] for member in cluster_dict[cluster_id].members]
            # Create a lineage list for all sequences in the cluster
            for num_id in fasta_record_objects:
                if header_registry[num_id].formatted in members:
                    lineages.append(fasta_record_objects[num_id].lineage)
            cleaned_lineages = clean_lineage_list(lineages)
            cluster_dict[cluster_id].lca = megan_lca(cleaned_lineages)
            # For debugging
            # if len(lineages) != len(cleaned_lineages) and len(lineages) > 1:
            #     print("Before:")
            #     for l in lineages:
            #         print(l)
            #     print("After:")
            #     for l in cleaned_lineages:
            #         print(l)
            #     print("LCA:", cluster_dict[cluster_id].lca)
            lineages.clear()
    else:
        cluster_dict = None

    ##
    # Swap sequences in 'guarantee' for the representatives, creating new clusters
    ##
    if args.guarantee and args.uc:
        # We don't want to make the tree redundant so instead of simply adding the sequences in guarantee,
        #  we will swap them for their respective representative sequences.
        # All important sequences become representative, even if multiple are in the same cluster
        num_swaps = 0
        nonredundant_guarantee_cluster_dict = dict()  # Will be used to replace cluster_dict
        expanded_cluster_id = 0
        for cluster_id in sorted(cluster_dict, key=int):
            if len(cluster_dict[cluster_id].members) == 0:
                nonredundant_guarantee_cluster_dict[cluster_id] = cluster_dict[cluster_id]
            else:
                contains_important_seq = False
                # The case where a member of a cluster is a guaranteed sequence, but not the representative
                representative = cluster_dict[cluster_id].representative
                for member in cluster_dict[cluster_id].members:
                    if member[0] in important_seqs.keys():
                        nonredundant_guarantee_cluster_dict[expanded_cluster_id] = Cluster(member[0])
                        nonredundant_guarantee_cluster_dict[expanded_cluster_id].members = []
                        nonredundant_guarantee_cluster_dict[expanded_cluster_id].lca = cluster_dict[cluster_id].lca
                        expanded_cluster_id += 1
                        contains_important_seq = True
                if contains_important_seq and representative not in important_seqs.keys():
                    num_swaps += 1
                elif contains_important_seq and representative in important_seqs.keys():
                    # So there is no opportunity for the important representative sequence to be swapped, clear members
                    cluster_dict[cluster_id].members = []
                    nonredundant_guarantee_cluster_dict[cluster_id] = cluster_dict[cluster_id]
                else:
                    nonredundant_guarantee_cluster_dict[cluster_id] = cluster_dict[cluster_id]
            expanded_cluster_id += 1
        cluster_dict = nonredundant_guarantee_cluster_dict

    # TODO: Taxonomic normalization

    ##
    # Set the cluster-specific values for ReferenceSequence objects
    ##
    if args.uc and not args.headless:
        # Allow user to select the representative sequence based on organism name, sequence length and similarity
        fasta_record_objects = present_cluster_rep_options(cluster_dict,
                                                           fasta_record_objects,
                                                           header_registry,
                                                           important_seqs)
    elif args.uc and args.headless:
        finalize_cluster_reps(cluster_dict, fasta_record_objects, header_registry)
    else:
        for num_id in fasta_record_objects:
            fasta_record_objects[num_id].cluster_rep = True
            # fasta_record_objects[num_id].cluster_lca is left empty

    logging.info("Detecting outlier reference sequences... ")
    outlier_test_fasta_dict = order_dict_by_lineage(fasta_record_objects)
    create_new_ref_fasta(od_input, outlier_test_fasta_dict)
    od_input_m = '.'.join(od_input.split('.')[:-1]) + ".mfa"
    od_output = args.output_dir + "outliers.fasta"
    # Perform MSA with MAFFT
    run_mafft(args.executables["mafft"], od_input, od_input_m, args.num_threads)
    # Run OD-seq on MSA to identify outliers
    run_odseq(args.executables["OD-seq"], od_input_m, od_output, args.num_threads)
    # Remove outliers from fasta_record_objects collection
    outlier_seqs = read_fasta_to_dict(od_output)
    outlier_names = list()
    for seq_name in outlier_seqs:
        for seq_num_id in fasta_record_objects:
            ref_seq = fasta_record_objects[seq_num_id]
            if ref_seq.short_id == seq_name:
                ref_seq.cluster_rep = False
                outlier_names.append(ref_seq.accession)
    logging.info("done.\n")
    logging.debug(str(len(outlier_seqs)) + " outlier sequences detected and discarded.\n\t" +
                  "\n\t".join([outseq for outseq in outlier_names]) + "\n")

    ##
    # Re-order the fasta_record_objects by their lineages (not phylogenetic, just alphabetical sort)
    # Remove the cluster members since they will no longer be used
    ##
    fasta_replace_dict = order_dict_by_lineage(fasta_record_objects)

    # For debugging. This is the finalized set of reference sequences:
    # for num_id in sorted(fasta_replace_dict, key=int):
    #     fasta_replace_dict[num_id].get_info()

    warnings = write_tax_ids(args, fasta_replace_dict, ref_pkg.lineage_ids)
    if warnings:
        logging.warning(warnings + "\n")

    logging.info("Generated the taxonomic lineage map " + ref_pkg.lineage_ids + "\n")
    taxonomic_summary = summarize_reference_taxa(args, fasta_replace_dict)
    logging.info(taxonomic_summary)
    marker_package.lowest_confident_rank = estimate_taxonomic_redundancy(args, fasta_replace_dict)

    ##
    # Perform multiple sequence alignment
    ##
    if args.multiple_alignment:
        create_new_ref_fasta(unaln_ref_fasta, fasta_replace_dict, True)
    else:
        create_new_ref_fasta(unaln_ref_fasta, fasta_replace_dict)

    if args.molecule == 'rrna':
        generate_cm_data(args, unaln_ref_fasta)
        args.multiple_alignment = True
    elif args.multiple_alignment is False:
        logging.info("Aligning the sequences using MAFFT... ")
        run_mafft(args.executables["mafft"], unaln_ref_fasta, ref_pkg.msa, args.num_threads)
        logging.info("done.\n")
    else:
        pass
    aligned_fasta_dict = read_fasta_to_dict(ref_pkg.msa)
    marker_package.num_reps = len(aligned_fasta_dict.keys())

    ##
    # Build the HMM profile from the aligned reference FASTA file
    ##
    if args.molecule == "rrna":
        # A .cm file has already been generated, no need for HMM
        pass
    else:
        logging.info("Building HMM profile... ")
        hmm_build_command = [args.executables["hmmbuild"],
                             args.final_output_dir + args.code_name + ".hmm",
                             ref_pkg.msa]
        stdout, hmmbuild_pro_returncode = launch_write_command(hmm_build_command)
        logging.info("done.\n")
        logging.debug("\n### HMMBUILD ###\n\n" + stdout)

        if hmmbuild_pro_returncode != 0:
            logging.error("hmmbuild did not complete successfully for:\n" +
                          ' '.join(hmm_build_command) + "\n")
            sys.exit(7)
    ##
    # Optionally trim with BMGE and create the Phylip multiple alignment file
    ##
    dict_for_phy = dict()
    if args.trim_align:
        logging.info("Running BMGE... ")
        trimmed_msa_file = trim_multiple_alignment(args.executables["BMGE.jar"], ref_pkg.msa, args.molecule)
        logging.info("done.\n")
        trimmed_aligned_fasta_dict = read_fasta_to_dict(trimmed_msa_file)
        if len(trimmed_aligned_fasta_dict) == 0:
            logging.warning("Trimming removed all your sequences. " +
                            "This could mean you have many non-homologous sequences or they are very dissimilar.\n" +
                            "Proceeding with the untrimmed multiple alignment instead.\n")
            for seq_name in aligned_fasta_dict:
                dict_for_phy[seq_name.split('_')[0]] = aligned_fasta_dict[seq_name]
        else:
            for seq_name in aligned_fasta_dict:
                dict_for_phy[seq_name.split('_')[0]] = trimmed_aligned_fasta_dict[seq_name]
        os.remove(trimmed_msa_file)
    else:
        for seq_name in aligned_fasta_dict:
            dict_for_phy[seq_name.split('_')[0]] = aligned_fasta_dict[seq_name]
    phy_dict = reformat_fasta_to_phy(dict_for_phy)
    write_phy_file(phylip_file, phy_dict)

    ##
    # Build the tree using either RAxML or FastTree
    ##
    marker_package.tree_tool = construct_tree(args, phylip_file, tree_output_dir)

    if os.path.exists(unaln_ref_fasta):
        os.remove(unaln_ref_fasta)
    if os.path.exists(phylip_file + ".reduced"):
        os.remove(phylip_file + ".reduced")
    if os.path.exists(args.final_output_dir + "fasta_reader_log.txt"):
        os.remove(args.final_output_dir + "fasta_reader_log.txt")

    if args.fast:
        if args.molecule == "prot":
            marker_package.model = "LG"
        else:
            marker_package.model = "GTRGAMMA"
    else:
        annotate_partition_tree(args.code_name,
                                fasta_replace_dict,
                                tree_output_dir + os.sep + "RAxML_bipartitions." + args.code_name)
        marker_package.model = find_model_used(tree_output_dir + os.sep + "RAxML_info." + args.code_name)

    # Build the regression model of placement distances to taxonomic ranks
    marker_package.pfit, _, _ = regress_rank_distance(args, ref_pkg, accession_lineage_map, aligned_fasta_dict)

    ##
    # Finish validating the file and append the reference package build parameters to the master table
    ##
    ref_pkg.validate(marker_package.num_reps)
    param_file = args.treesapp + "data" + os.sep + "tree_data" + os.sep + "ref_build_parameters.tsv"
    update_build_parameters(param_file, marker_package)

    logging.info("Data for " + args.code_name + " has been generated successfully.\n")
    terminal_commands(args.final_output_dir, args.code_name)


if __name__ == "__main__":
    main()
