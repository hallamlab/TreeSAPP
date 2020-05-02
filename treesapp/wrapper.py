import sys
import os
import time
import re
import glob
import logging
from shutil import copy

from treesapp.external_command_interface import launch_write_command, setup_progress_bar
from treesapp.classy import CommandLineFarmer
from .fasta import read_fasta_to_dict
from .utilities import remove_dashes_from_msa


def construct_tree(executables: dict, molecule: str, multiple_alignment_file: str,
                   tree_output_dir, tree_file, tree_prefix, args):
    """
    Wrapper script for generating phylogenetic trees with either RAxML or FastTree from a multiple alignment

    :param executables: Dictionary containing paths to executables, crucially FastTree and RAxML
    :param molecule: Molecule type of the sequences being used to infer the phylogeny
    :param multiple_alignment_file: Path to the multiple sequence alignment file
    :param tree_output_dir: Path to the directory where output files should be written to
    :param tree_file: Path to write the inferred phylogenetic tree
    :param tree_prefix: Prefix to be used for the outputs
    :param args: Command-line arguments parsed using ArgParse
    :return: Stylized name of the tree-building software used
    """

    # Decide on the command to build the tree, make some directories and files when necessary
    if args.fast:
        tree_build_cmd = [executables["FastTree"]]
        if molecule == "rrna" or molecule == "dna":
            tree_build_cmd += ["-nt", "-gtr"]
        else:
            tree_build_cmd += ["-lg", "-wag"]
        tree_build_cmd += ["-out", tree_file]
        tree_build_cmd.append(multiple_alignment_file)
        tree_builder = "FastTree"
    else:
        tree_build_cmd = [executables["raxmlHPC"]]
        tree_build_cmd += ["-f", "a"]
        tree_build_cmd += ["-p", "12345"]
        tree_build_cmd += ["-x", "12345"]
        tree_build_cmd += ["-#", str(args.bootstraps)]
        tree_build_cmd += ["-s", multiple_alignment_file]
        tree_build_cmd += ["-n", tree_prefix]
        tree_build_cmd += ["-w", tree_output_dir]
        tree_build_cmd += ["-T", str(args.num_threads)]

        if args.raxml_model:
            tree_build_cmd += ["-m", args.raxml_model]
        elif args.molecule == "prot":
            tree_build_cmd += ["-m", "PROTGAMMAAUTO"]
        elif args.molecule == "rrna" or molecule == "dna":
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

    logging.info("Building phylogenetic tree with " + tree_builder + "... ")
    if args.fast:
        stdout, returncode = launch_write_command(tree_build_cmd, True)
        with open(tree_output_dir + os.sep + "FastTree_info." + tree_prefix, 'w') as fast_info:
            fast_info.write(stdout + "\n")
    else:
        stdout, returncode = launch_write_command(tree_build_cmd, False)
    logging.info("done.\n")

    if returncode != 0:
        logging.error(tree_builder + " did not complete successfully! " +
                      "Look in " + tree_output_dir + os.sep +
                      tree_builder + "_info." + tree_prefix + " for an error message.\n" +
                      tree_builder + " command used:\n" + ' '.join(tree_build_cmd) + "\n")
        sys.exit(13)

    return tree_builder


def launch_evolutionary_placement_queries(executables, tree_dir, phy_files, marker_build_dict, output_dir, num_threads):
    """
    Run EPA through RAxML using Phylip files containing the reference and query sequences, and the reference trees

    :param executables:
    :param tree_dir:
    :param phy_files:
    :param marker_build_dict:
    :param output_dir:
    :param num_threads:
    :return: None
    """
    logging.info("Running RAxML... coffee?\n")

    start_time = time.time()

    raxml_calls = 0
    # Maximum-likelihood sequence placement analyses
    denominator_reference_tree_dict = dict()
    for denominator in sorted(phy_files.keys()):
        if not isinstance(denominator, str):
            logging.error(str(denominator) + " is not string but " + str(type(denominator)) + "\n")
            raise AssertionError()
        # Establish the reference tree file to be used for this contig
        ref_marker = marker_build_dict[denominator]
        reference_tree_file = tree_dir + os.sep + ref_marker.cog + '_tree.txt'
        if denominator not in denominator_reference_tree_dict.keys():
            denominator_reference_tree_dict[denominator] = reference_tree_file
        for phy_file in phy_files[denominator]:
            query_name = re.sub("_hmm_purified.phy.*$", '', os.path.basename(phy_file))
            query_name = re.sub(marker_build_dict[denominator].cog, denominator, query_name)
            raxml_evolutionary_placement(executables["raxmlHPC"], reference_tree_file, phy_file,
                                         ref_marker.model, output_dir, query_name, num_threads)
            raxml_calls += 1

    end_time = time.time()
    hours, remainder = divmod(end_time - start_time, 3600)
    minutes, seconds = divmod(remainder, 60)
    logging.debug("\tRAxML time required: " +
                  ':'.join([str(hours), str(minutes), str(round(seconds, 2))]) + "\n")
    logging.debug("\tRAxML was called " + str(raxml_calls) + " times.\n")

    return


def raxml_evolutionary_placement(raxml_exe: str, reference_tree_file: str, multiple_alignment: str, model: str,
                                 output_dir: str, query_name: str, num_threads=2):
    """
    A wrapper for RAxML's evolutionary placement algorithm (EPA)
        1. checks to ensure the output files do not already exist, and removes them if they do
        2. ensures the output directory is an absolute path, satisfying RAxML
        3. Runs RAxML with the provided parameters
        4. Renames the files for consistency in TreeSAPP
    :param raxml_exe: Path to the RAxML executable to be used
    :param reference_tree_file: The reference tree for evolutionary placement to operate on
    :param multiple_alignment: Path to a multiple alignment file containing reference and query sequences
    :param model: The substitution model to be used by RAxML e.g. PROTGAMMALG, GTRCAT
    :param output_dir: Path to write the EPA outputs
    :param query_name: Prefix name for all of the output files
    :param num_threads: Number of threads EPA should use (default = 2)
    :return: A dictionary of files written by RAxML's EPA that are used by TreeSAPP. For example epa_files["jplace"]
    """
    epa_files = dict()
    ##
    # Start with some housekeeping - are the inputs looking alright?
    # Do the outputs already exist?
    # Is the output directory an absolute path?
    ##
    if not os.path.isabs(output_dir):
        output_dir = os.getcwd() + os.sep + output_dir
    if output_dir[-1] != os.sep:
        output_dir += os.sep

    if model is None:
        logging.error("No substitution model provided for evolutionary placement of " + query_name + ".\n")
        raise AssertionError()

    # Determine the output file names, and remove any pre-existing output files
    if not isinstance(reference_tree_file, str):
        logging.error(str(reference_tree_file) + " is not string but " + str(type(reference_tree_file)) + "\n")
        raise AssertionError()

    if len(reference_tree_file) == 0:
        logging.error("Could not find reference tree for " + query_name + " to be used by EPA.\n")
        raise AssertionError()

    # This is the final set of files that will be written by RAxML's EPA algorithm
    epa_files["stdout"] = output_dir + query_name + '_RAxML.txt'
    epa_info = output_dir + 'RAxML_info.' + query_name
    epa_files["info"] = output_dir + query_name + '.RAxML_info.txt'
    epa_labelled_tree = output_dir + 'RAxML_labelledTree.' + query_name
    epa_tree = output_dir + 'RAxML_originalLabelledTree.' + query_name
    epa_files["tree"] = output_dir + query_name + '.originalRAxML_labelledTree.txt'
    epa_classification = output_dir + 'RAxML_classification.' + query_name
    epa_files["classification"] = output_dir + query_name + '.RAxML_classification.txt'
    epa_files["jplace"] = output_dir + "RAxML_portableTree." + query_name + ".jplace"
    epa_entropy = output_dir + "RAxML_entropy." + query_name
    epa_weights = output_dir + "RAxML_classificationLikelihoodWeights." + query_name

    for raxml_file in [epa_info, epa_labelled_tree, epa_tree, epa_classification, epa_entropy, epa_weights]:
        try:
            os.remove(raxml_file)
        except OSError:
            pass

    # Set up the command to run RAxML
    raxml_command = [raxml_exe,
                     '-m', model,
                     '-T', str(int(num_threads)),
                     '-s', multiple_alignment,
                     '-t', reference_tree_file,
                     '-G', str(0.2),
                     "--epa-prob-threshold=" + str(0.10),
                     '-f', 'v',
                     '-n', query_name,
                     '-w', output_dir,
                     '>', epa_files["stdout"]]
    launch_write_command(raxml_command)

    # Rename the RAxML output files
    if os.path.exists(epa_info):
        copy(epa_info, epa_files["info"])
        os.remove(epa_info)
    if os.path.exists(epa_classification):
        copy(epa_classification, epa_files["classification"])
        os.remove(epa_classification)
    if os.path.exists(epa_tree):
        copy(epa_tree, epa_files["tree"])
        os.remove(epa_tree)
    else:
        logging.error("Some files were not successfully created for " + query_name + "\n" +
                      "Check " + epa_files["stdout"] + " for an error!\n")
        sys.exit(3)
    # Remove useless files
    if os.path.exists(epa_labelled_tree):
        os.remove(epa_labelled_tree)
        os.remove(epa_weights)
        os.remove(epa_entropy)

    return epa_files


def trimal_command(executable, mfa_file, trimmed_msa_file):
    trim_command = [executable,
                    '-in', mfa_file,
                    '-out', trimmed_msa_file,
                    '-gt', str(0.02)]
    return trim_command


def bmge_command(executable, mfa_file, trimmed_msa_file, molecule):
    if molecule == "prot":
        bmge_settings = ["-t", "AA", "-m", "BLOSUM30"]
    else:
        bmge_settings = ["-t", "DNA", "-m", "DNAPAM100:2"]
    trim_command = ["java", "-Xmx512m", "-jar", executable]
    trim_command += bmge_settings
    trim_command += ["-g", "0.99:0.33"]  # Specifying the gap rate per_sequence:per_character
    trim_command += ['-i', mfa_file,
                     '-of', trimmed_msa_file]
    return trim_command


def hmmalign_command(executable, ref_aln, ref_profile, input_fasta, output_multiple_alignment):
    malign_command = [executable,
                      '--mapali', ref_aln,
                      '--outformat', 'Stockholm',
                      ref_profile, input_fasta,
                      '>', output_multiple_alignment]

    return malign_command


def profile_aligner(executables, ref_aln, ref_profile, input_fasta, output_sto, kind="functional"):
    """
    A wrapper for both cmalign and hmmalign for performing profile-based multiple sequence alignment
    :param executables: A dictionary containing keys "cmalign" and "hmmalign"
    :param ref_aln: Path to a FASTA or Stockholm file with the multiple alignment pattern
    :param ref_profile: Path to the HMM or CM profile for the reference gene
    :param input_fasta: Path to the FASTA containing query sequences
    :param output_sto: Name of the output Stockholm formatted file
    :param kind: The type of marker gene being analyzed [functional (default), phylogenetic, phylogenetic_rRNA]
    :return:
    """

    if kind == "phylogenetic_rRNA":
        malign_command = hmmalign_command(executables["cmalign"], ref_aln, ref_profile, input_fasta, output_sto)
    else:
        malign_command = hmmalign_command(executables["hmmalign"], ref_aln, ref_profile, input_fasta, output_sto)

    stdout, returncode = launch_write_command(malign_command)
    if returncode != 0:
        logging.error("Multiple alignment failed for " + input_fasta + ". Command used:\n" +
                      ' '.join(malign_command) + " output:\n" + stdout + "\n")
        sys.exit(3)
    return stdout


def run_papara(executable, tree_file, ref_alignment_phy, query_fasta, molecule):
    papara_command = [executable]
    papara_command += ["-t", tree_file]
    papara_command += ["-s", ref_alignment_phy]
    papara_command += ["-q", query_fasta]
    if molecule == "prot":
        papara_command.append("-a")

    stdout, ret_code = launch_write_command(papara_command)
    if ret_code != 0:
        logging.error("PaPaRa did not complete successfully!\n" +
                      "Command used:\n" + ' '.join(papara_command) + "\n")
        sys.exit(3)
    return stdout


def cluster_new_reference_sequences(update_tree, args, new_ref_seqs_fasta):
    logging.info("Clustering sequences at %s percent identity with USEARCH... " % str(update_tree.cluster_id))

    usearch_command = [args.executables["usearch"]]
    usearch_command += ["-sortbylength", new_ref_seqs_fasta]
    usearch_command += ["-fastaout", update_tree.Output + "usearch_sorted.fasta"]
    usearch_command += ["--log", update_tree.Output + os.sep + "usearch_sort.log"]
    # usearch_command += ["1>", "/dev/null", "2>", "/dev/null"]

    launch_write_command(usearch_command)

    uclust_id = "0." + str(int(update_tree.cluster_id))
    try:
        float(uclust_id)
    except ValueError:
        logging.error("Weird formatting of cluster_id: " + uclust_id + "\n")

    uclust_command = [args.executables["usearch"]]
    uclust_command += ["-cluster_fast", update_tree.Output + "usearch_sorted.fasta"]
    uclust_command += ["--id", uclust_id]
    uclust_command += ["--centroids", update_tree.Output + "uclust_" + update_tree.COG + ".fasta"]
    uclust_command += ["--uc", update_tree.Output + "uclust_" + update_tree.COG + ".uc"]
    uclust_command += ["--log", update_tree.Output + os.sep + "usearch_cluster.log"]
    # uclust_command += ["1>", "/dev/null", "2>", "/dev/null"]

    launch_write_command(uclust_command)

    logging.info("done.\n")

    return


def cluster_sequences(uclust_exe, fasta_input, uclust_prefix, similarity=0.60):
    """
    Wrapper function for clustering a FASTA file at some similarity using usearch's cluster_fast algorithm

    :param uclust_exe: Path to the usearch executable
    :param fasta_input: FASTA file for which contained sequences will be clustered
    :param uclust_prefix: Prefix for the output files
    :param similarity: The proportional similarity to cluster input sequences
    :return: None
    """
    logging.info("Clustering sequences with UCLUST... ")
    uclust_cmd = [uclust_exe]
    uclust_cmd += ["-cluster_fast", fasta_input]
    uclust_cmd += ["-id", str(similarity)]
    uclust_cmd += ["-sort", "length"]
    uclust_cmd += ["-centroids", uclust_prefix + ".fa"]
    uclust_cmd += ["--uc", uclust_prefix + ".uc"]
    logging.info("done.\n")
    stdout, returncode = launch_write_command(uclust_cmd)

    if returncode != 0:
        logging.error("UCLUST did not complete successfully! Command used:\n" +
                      ' '.join(uclust_cmd) + "\n")
        sys.exit(13)

    logging.debug(stdout)
    return


def build_hmm_profile(hmmbuild_exe, msa_in, output_hmm):
    logging.debug("Building HMM profile... ")
    hmm_build_command = [hmmbuild_exe, output_hmm, msa_in]
    stdout, hmmbuild_pro_returncode = launch_write_command(hmm_build_command)
    logging.debug("done.\n")

    if hmmbuild_pro_returncode != 0:
        logging.error("hmmbuild did not complete successfully for:\n" +
                      ' '.join(hmm_build_command) + "\n")
        sys.exit(7)


def run_prodigal(args, fasta_file, output_file, nucleotide_orfs=None):
    prodigal_command = [args.executables["prodigal"]]
    prodigal_command += ["-i", fasta_file]
    prodigal_command += ["-p", "meta"]
    prodigal_command += ["-a", output_file]
    if nucleotide_orfs:
        prodigal_command += ["-d", nucleotide_orfs]
    stdout, proc_code = launch_write_command(prodigal_command)

    if proc_code != 0:
        logging.error("Prodigal did not complete successfully!\n" +
                      "Command used:\n" + ' '.join(prodigal_command), "err", "\n")
        sys.exit(3)
    return


def run_hmmsearch(hmmsearch_exe: str, hmm_profile: str, query_fasta: str, output_dir: str, num_threads=2):
    """
    Function for searching a fasta file with a profile HMM

    :param hmmsearch_exe: Path to the executable for hmmsearch
    :param hmm_profile: Path to the HMM profile file
    :param query_fasta: Path to the FASTA file to be queried by the profile
    :param output_dir: Path to the directory for writing the outputs
    :param num_threads: Number of threads to be used by hmmsearch
    :return:
    """
    # Find the name of the HMM. Use it to name the output file
    rp_marker = re.sub(".hmm", '', os.path.basename(hmm_profile))
    domtbl = output_dir + rp_marker + "_to_ORFs_domtbl.txt"

    # Basic hmmsearch command
    hmmsearch_command_base = [hmmsearch_exe]
    hmmsearch_command_base += ["--cpu", str(num_threads)]
    hmmsearch_command_base.append("--noali")
    # Customize the command for this input and HMM
    final_hmmsearch_command = hmmsearch_command_base + ["--domtblout", domtbl]
    final_hmmsearch_command += [hmm_profile, query_fasta]
    stdout, ret_code = launch_write_command(final_hmmsearch_command)

    # Check to ensure the job finished properly
    if ret_code != 0:
        logging.error("hmmsearch did not complete successfully! Output:\n" + stdout + "\n" +
                      "Command used:\n" + ' '.join(final_hmmsearch_command) + "\n")
        sys.exit(13)

    return [domtbl]


def hmmsearch_orfs(hmmsearch_exe, hmm_dir, marker_build_dict, fasta_file, output_dir, num_threads=2):
    hmm_domtbl_files = list()
    nucl_target_hmm_files = list()
    prot_target_hmm_files = list()

    # Find all of the available HMM profile files
    try:
        os.path.isdir(hmm_dir)
    except IOError:
        logging.error(hmm_dir + "does not exist!")
        sys.exit(3)
    hmm_files = glob.glob(hmm_dir + "*.hmm")

    if len(hmm_files) == 0:
        logging.error(hmm_dir + "does not contain any files with '.hmm' extension... so no HMMs.\n")
        sys.exit(3)

    # Filter the HMM files to only the target markers
    for marker_code in marker_build_dict:
        ref_marker = marker_build_dict[marker_code]
        hmm_profile = hmm_dir + ref_marker.cog + ".hmm"
        if hmm_profile not in hmm_files:
            logging.error("Unable to locate HMM-profile for " + ref_marker.cog + "(" + marker_code + ").\n")
        else:
            if ref_marker.molecule == "prot":
                prot_target_hmm_files.append(hmm_profile)
            else:
                nucl_target_hmm_files.append(hmm_profile)

    acc = 0.0
    logging.info("Searching for marker proteins in ORFs using hmmsearch.\n")
    step_proportion = setup_progress_bar(len(prot_target_hmm_files) + len(nucl_target_hmm_files))

    # Create and launch the hmmsearch commands iteratively.
    for hmm_file in prot_target_hmm_files:
        # TODO: Parallelize this by allocating no more than 2 threads per process
        hmm_domtbl_files += run_hmmsearch(hmmsearch_exe, hmm_file, fasta_file, output_dir, num_threads)

        # Update the progress bar
        acc += 1.0
        if acc >= step_proportion:
            acc -= step_proportion
            time.sleep(0.1)
            sys.stdout.write("-")
            sys.stdout.flush()
    sys.stdout.write("-]\n")

    return hmm_domtbl_files


def generate_blast_database(args, fasta, molecule, prefix, multiple=True):
    """

    :param args:
    :param fasta: File to make a BLAST database for
    :param molecule: 'prot' or 'nucl' - necessary argument for makeblastdb
    :param prefix: prefix string for the output BLAST database
    :param multiple: Flag indicating the input `fasta` is a MSA. Alignment information is removed prior to makeblastdb
    :return:
    """

    # Remove the multiple alignment information from fasta_replaced_file and write to fasta_mltree
    blastdb_out = prefix + ".fa"
    if multiple:
        if blastdb_out == fasta:
            logging.error("prefix.fa is the same as " + fasta + " and would be overwritten!\n")
            sys.exit(13)
        remove_dashes_from_msa(fasta, blastdb_out)
        blastdb_in = blastdb_out
    else:
        blastdb_in = fasta

    logging.info("Making the BLAST database for " + blastdb_in + "... ")

    # Format the `makeblastdb` command
    makeblastdb_command = [args.executables["makeblastdb"]]
    makeblastdb_command += ["-in", blastdb_in]
    makeblastdb_command += ["-out", blastdb_out]
    makeblastdb_command += ["-input_type", "fasta"]
    makeblastdb_command += ["-dbtype", molecule]

    # Launch the command
    stdout, makeblastdb_pro_returncode = launch_write_command(makeblastdb_command)

    logging.info("done\n")

    return stdout, blastdb_out


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
    mafft_align_command += ["--randomseed", str(12345)]
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


def get_msa_trim_command(executables, mfa_file, molecule, tool="BMGE"):
    """
    Trims/masks/filters the multiple sequence alignment using either BMGE or trimAl

    :param executables: A dictionary mapping software to a path of their respective executable
    :param mfa_file: Name of a MSA file
    :param molecule: prot | dna
    :param tool: Name of the software to use for trimming [BMGE|trimAl]
    Returns file name of the trimmed multiple alignment file in FASTA format
    """
    f_ext = mfa_file.split('.')[-1]
    if not re.match("mfa|fasta|phy|fa", f_ext):
        logging.error("Unsupported file format: '" + f_ext + "'\n")
        sys.exit(5)

    trimmed_msa_file = re.sub('.' + re.escape(f_ext), '-' + re.escape(tool) + ".fasta", mfa_file)
    if tool == "trimAl":
        trim_command = trimal_command(executables["trimal"], mfa_file, trimmed_msa_file)
    elif tool == "BMGE":
        trim_command = bmge_command(executables["BMGE.jar"], mfa_file, trimmed_msa_file, molecule)
    else:
        logging.error("Unsupported trimming software requested: '" + tool + "'")
        sys.exit(5)

    return trim_command, trimmed_msa_file


def filter_multiple_alignments(executables, concatenated_mfa_files, marker_build_dict, n_proc=1, tool="BMGE"):
    """
    Runs BMGE using the provided lists of the concatenated hmmalign files, and the number of sequences in each file.

    :param executables: A dictionary mapping software to a path of their respective executable
    :param concatenated_mfa_files: A dictionary containing f_contig keys mapping to a FASTA or Phylip sequential file
    :param marker_build_dict:
    :param n_proc: The number of parallel processes to be launched for alignment trimming
    :param tool: The software to use for alignment trimming
    :return: A list of files resulting from BMGE multiple sequence alignment masking.
    """
    logging.info("Running " + tool + "... ")

    start_time = time.time()
    task_list = list()
    trimmed_output_files = {}

    for denominator in sorted(concatenated_mfa_files.keys()):
        if denominator not in trimmed_output_files:
            trimmed_output_files[denominator] = []
        mfa_files = concatenated_mfa_files[denominator]
        for concatenated_mfa_file in mfa_files:
            trim_command, trimmed_msa_file = get_msa_trim_command(executables, concatenated_mfa_file,
                                                                  marker_build_dict[denominator].molecule, tool)
            trimmed_output_files[denominator].append(trimmed_msa_file)
            task_list.append(trim_command)

    if len(task_list) > 0:
        cl_farmer = CommandLineFarmer("Multiple alignment trimming with " + tool, n_proc)
        cl_farmer.add_tasks_to_queue(task_list)

        cl_farmer.task_queue.close()
        cl_farmer.task_queue.join()

    logging.info("done.\n")

    end_time = time.time()
    hours, remainder = divmod(end_time - start_time, 3600)
    minutes, seconds = divmod(remainder, 60)
    logging.debug("\t" + tool + " time required: " +
                  ':'.join([str(hours), str(minutes), str(round(seconds, 2))]) + "\n")
    return trimmed_output_files
