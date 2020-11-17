import sys
import os
import time
import re
import glob
import logging
from shutil import copy, rmtree

from tqdm import tqdm

from treesapp.external_command_interface import launch_write_command, CommandLineFarmer
from treesapp.fasta import read_fasta_to_dict


def estimate_ml_model(modeltest_exe: str, msa: str, output_prefix: str, molecule: str, threads=1) -> str:
    """
    Eventually this function will be a wrapper for ModelTest-ng or IQTree's ModelFinder.

    :param modeltest_exe: Path to the executable for ModelTest-NG
    :param msa: Path to a Multiple Sequence Alignment file
    :param output_prefix: The prefix for the output files (includes directory paths)
    :param molecule: A string indicating the molecule-type of the reference package: 'rrna', 'prot' or 'dna'
    :param threads: Number of threads for ModelTest-NG to use
    :return: A RAxML-ng and EPA-ng compatible string representing the substitution model
    """
    if molecule == "prot":
        dt = "aa"
        model_candidates = "JTT,LG,WAG"
    else:
        dt = "nt"
        model_candidates = "GTR"

    model_find_cmd = [modeltest_exe]
    model_find_cmd += ["--input", msa]
    model_find_cmd += ["--output", output_prefix]
    model_find_cmd += ["--processes", threads]
    model_find_cmd += ["--topology", "fixed-mp"]
    model_find_cmd += ["--datatype", dt]
    model_find_cmd += ["--template", "raxml"]
    model_find_cmd += ["--frequencies", "f"]
    model_find_cmd += ["--models", model_candidates]

    stdout, returncode = launch_write_command(model_find_cmd)
    if returncode != 0:
        logging.error("{} could not determine a substitution model for your sequences in '{}'.\n".format(modeltest_exe,
                                                                                                         molecule))
        sys.exit(13)

    # TODO: Parse the optimal model determined by ModelTest-NG and return it
    evo_model = ""

    return evo_model


def model_parameters(raxml_exe: str, ref_msa: str, tree_file: str, output_prefix: str, model: str, threads=2) -> str:
    """
    Wrapper function for RAxML-ng's `evaluate` sub-command that generates a file to be used by EPA-ng.
    This file, in reality, is part of a reference package as it is reference tree and MSA dependent
    and will need to be provided for phylogenetic placement.

    :param raxml_exe: Path to the raxml-ng executable
    :param ref_msa: Path to the reference package's multiple sequence alignment (FASTA)
    :param tree_file: Path to the reference package's tree (Newick)
    :param output_prefix: Prefix path for the output files
    :param model: Substitution model (and other parameters e.g. Gamma model of rate heterogeneity) used to build tree
    :param threads: The number of threads that should be used by RAxML-NG
    :return: Path to the bestModel file that can be used by epa-ng for phylogenetic placement
    """
    output_prefix += "_evaluate"
    model_params_file = output_prefix + ".raxml.bestModel"
    model_eval_cmd = [raxml_exe, "--evaluate"]
    model_eval_cmd += ["--msa", ref_msa]
    model_eval_cmd += ["--tree", tree_file]
    model_eval_cmd += ["--prefix", output_prefix]
    model_eval_cmd += ["--model", model]
    model_eval_cmd += ["--threads", "auto{{{}}}".format(threads)]
    model_eval_cmd += ["--seed", str(12345)]
    model_eval_cmd += ["--workers", "auto{{{}}}".format(threads)]
    model_eval_cmd.append("--force")

    logging.debug("Evaluating phylogenetic tree with RAxML-NG... ")
    stdout, returncode = launch_write_command(model_eval_cmd)
    logging.debug("done.\n")

    if returncode != 0:
        logging.error("{} did not complete successfully! Look in {}_info.txt for an error message.\n"
                      "RAxML-NG command used:\n{}\n".format(raxml_exe, output_prefix, ' '.join(model_eval_cmd)))
        sys.exit(13)

    return model_params_file


def bootstrap_tree_raxml(raxml_exe: str, multiple_alignment: str, model: str, tree_prefix: str,
                         mre=True, bootstraps=1000, num_threads=2) -> str:
    bootstrap_cmd = [raxml_exe, "--bootstrap"]
    bootstrap_cmd += ["--msa", multiple_alignment]
    bootstrap_cmd += ["--model", model]
    if mre:
        bootstrap_cmd += ["--bs-trees", "autoMRE{%d}" % bootstraps]
    else:
        bootstrap_cmd += ["--bs-trees", str(bootstraps)]
    bootstrap_cmd += ["--prefix", tree_prefix]
    bootstrap_cmd += ["--seed", str(12345)]
    bootstrap_cmd += ["--threads", str(num_threads)]

    logging.info("Bootstrapping reference tree with RAxML-NG... ")
    launch_write_command(bootstrap_cmd)
    logging.info("done.\n")

    bootstraps_file = tree_prefix + ".raxml.bootstraps"
    if not os.path.isfile(bootstraps_file):
        logging.error("Unable to find bootstrap file '%s'.\n" % bootstraps_file)
        sys.exit(17)

    return bootstraps_file


def support_tree_raxml(raxml_exe: str, ref_tree: str, ref_msa: str, model: str, tree_prefix: str,
                       mre=False, n_bootstraps=1000, num_threads=2) -> str:
    bootstraps = bootstrap_tree_raxml(raxml_exe, ref_msa, model, tree_prefix, mre, n_bootstraps, num_threads)

    support_cmd = [raxml_exe, "--support"]
    support_cmd += ["--tree", ref_tree]
    support_cmd += ["--bs-trees", bootstraps]
    support_cmd += ["--prefix", tree_prefix]
    support_cmd += ["--threads", str(num_threads)]

    logging.info("Calculating bootstrap support for node in reference tree with RAxML-NG... ")
    launch_write_command(support_cmd)
    logging.info("done.\n")

    support_file = tree_prefix + ".raxml.support"
    if not os.path.isfile(support_file):
        logging.error("Unable to find support file '%s'.\n" % support_file)
        sys.exit(17)

    return support_file


def construct_tree(tree_builder: str, executables: dict, evo_model: str, multiple_alignment_file: str,
                   tree_output_dir: str, tree_prefix: str, num_threads=2) -> str:
    """
    Wrapper script for generating phylogenetic trees with either RAxML or FastTree from a multiple alignment

    :param tree_builder: String indicating which phylogeny inference software is to be used. Current options are
    FastTree and RAxML-NG
    :param executables: Dictionary containing paths to executables, crucially FastTree and RAxML
    :param evo_model: The substitution model (and possible gamma rate heterogeneity) string (e.g. GTR+G4)
    :param multiple_alignment_file: Path to the multiple sequence alignment file
    :param tree_output_dir: Path to the directory where output files should be written to
    :param tree_prefix: Prefix to be used for the outputs
    :param num_threads: Number of threads to use (for RAxML-NG only)
    :return: Stylized name of the tree-building software used
    """

    # Decide on the command to build the tree, make some directories and files when necessary
    logging.info("Building phylogenetic tree with " + tree_builder + "... ")
    if tree_builder == "FastTree":
        best_tree = tree_output_dir + tree_prefix + ".FastTree.nwk"
        tree_build_cmd = [executables["FastTree"]]
        if re.search(r"GTR", evo_model):
            tree_build_cmd += ["-nt", "-gtr"]
        else:
            tree_build_cmd += ["-lg", "-gamma", "-cat", str(4)]
        tree_build_cmd += ["-out", best_tree]
        tree_build_cmd.append(multiple_alignment_file)

        stdout, returncode = launch_write_command(tree_build_cmd)
        with open(tree_output_dir + tree_prefix + ".FastTree.log", 'w') as fast_info:
            fast_info.write(stdout + "\n")
    elif tree_builder == "RAxML-NG":
        best_tree = tree_output_dir + tree_prefix + ".raxml.bestTree"
        tree_build_cmd = [executables["raxml-ng"], "--search"]
        tree_build_cmd += ["--prefix", tree_output_dir + tree_prefix]
        tree_build_cmd += ["--msa", multiple_alignment_file]
        tree_build_cmd += ["--model", evo_model]
        # tree_build_cmd += ["--msa-format", "PHYLIP"]  # File isn't read properly with this parameter, use auto-detect
        tree_build_cmd += ["--seed", str(12345)]
        tree_build_cmd += ["--threads", str(num_threads)]
        # tree_build_cmd += ["--tree", "rand{1},pars{1}"]  # For debugging, alternatively could use '--search1'

        stdout, returncode = launch_write_command(tree_build_cmd)
    else:
        logging.error("Unrecognized software '{}'.\n".format(tree_builder))
        sys.exit(5)

    logging.info("done.\n")
    logging.debug(stdout + "\n")

    if returncode != 0:
        logging.error("{0} did not complete successfully! Look in {1} for an error message.\n"
                      "{0} command used:\n{2}\n".format(tree_builder,
                                                        tree_output_dir + '.'.join([tree_prefix, tree_builder, "log"]),
                                                        ' '.join(tree_build_cmd)))
        sys.exit(13)

    return best_tree


def launch_evolutionary_placement_queries(executables: dict, split_msa_files: dict,
                                          refpkg_dict: dict, output_dir: str,
                                          num_threads: int) -> None:
    """
    Run EPA-ng using FASTA files containing the reference and query sequences, and the reference trees

    :param executables: Dictionary of executables where executable name strings are keys and paths are values
    :param split_msa_files: Dictionary of TreeSAPP refpkg code (denominator) keys indexing a list of
     namedtuple instances called MSAs. Each instance has 'ref' and 'query' variables referring to the
     MSA files (FASTA formatted) with aligned reference sequences and aligned query sequences, respectively.
    :param refpkg_dict: Dictionary of ReferencePackage instances indexed by their TreeSAPP refpkg code (denominator)
    :param output_dir: Path to write the EPA-ng outputs
    :param num_threads: Number of threads to use during placement
    :return: None
    """
    logging.info("Running EPA... ")

    start_time = time.time()

    epa_calls = 0
    # Maximum-likelihood sequence placement analyses
    for refpkg_name in sorted(split_msa_files.keys()):
        if not isinstance(refpkg_name, str):
            logging.error("RefPkg name key '{}' is not string but {}\n.".format(refpkg_name, type(refpkg_name)))
            raise AssertionError
        ref_pkg = refpkg_dict[refpkg_name]
        for split_msa in split_msa_files[refpkg_name]:
            query_name = re.sub("_queries.mfa", '', os.path.basename(split_msa.query))
            query_name = re.sub(ref_pkg.prefix, refpkg_name, query_name)
            # Find the query names
            raxml_evolutionary_placement(executables["epa-ng"], ref_pkg.f__tree, split_msa.ref, ref_pkg.f__model_info,
                                         split_msa.query, query_name, output_dir, num_threads)
            epa_calls += 1

    end_time = time.time()
    hours, remainder = divmod(end_time - start_time, 3600)
    minutes, seconds = divmod(remainder, 60)
    logging.info("done.\n")

    logging.debug("\tEPA-ng time required: " +
                  ':'.join([str(hours), str(minutes), str(round(seconds, 2))]) + "\n")
    logging.debug("\tEPA-ng was called " + str(epa_calls) + " times.\n")

    return


def raxml_evolutionary_placement(epa_exe: str, refpkg_tree: str, refpkg_msa: str, refpkg_model: str,
                                 query_msa: str, query_name: str, output_dir: str, num_threads=2):
    """
    A wrapper for evolutionary placement algorithm (EPA) next-generation
        1. checks to ensure the output files do not already exist, and removes them if they do
        2. ensures the output directory is an absolute path, satisfying EPA
        3. Runs EPA with the provided parameters
        4. Renames the files for consistency in TreeSAPP

    :param epa_exe: Path to the EPA-ng executable to be used
    :param refpkg_tree: The reference tree for evolutionary placement to operate on
    :param refpkg_msa: The reference multiple sequence alignment for the reference package (FASTA)
    :param refpkg_model: The substitution model to be used by EPA e.g. PROTGAMMALG, GTRCAT
    :param query_msa: Path to a multiple alignment file containing aligned query sequences (FASTA)
    :param query_name: Prefix name for all of the output files
    :param output_dir: Path to write the EPA outputs
    :param num_threads: Number of threads EPA should use (default = 2)
    :return: A dictionary of files written by EPA-ng that are used by TreeSAPP. For example epa_files["jplace"]
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

    if refpkg_model is None:
        logging.error("No substitution model provided for evolutionary placement of " + query_name + ".\n")
        raise AssertionError()

    # Determine the output file names, and remove any pre-existing output files
    if not isinstance(refpkg_tree, str):
        logging.error(str(refpkg_tree) + " is not string but " + str(type(refpkg_tree)) + "\n")
        raise AssertionError()

    if len(refpkg_tree) == 0:
        logging.error("Could not find reference tree for " + query_name + " to be used by EPA-ng.\n")
        raise AssertionError()

    # This is the final set of files that will be written by EPA-ng
    epa_files["stdout"] = output_dir + query_name + '_EPA.txt'
    epa_info = output_dir + 'epa_info.log'
    epa_files["info"] = output_dir + query_name + '.EPA_info.txt'
    epa_jplace = output_dir + "epa_result.jplace"
    epa_files["jplace"] = output_dir + "epa_result." + query_name + ".jplace"

    for raxml_file in [epa_info, epa_jplace]:
        try:
            os.remove(raxml_file)
        except OSError:
            pass

    # Set up the command to run EPA-ng
    epa_command = [epa_exe,
                   '-s', refpkg_msa,
                   '-t', refpkg_tree,
                   '-q', query_msa,
                   "--model", refpkg_model,
                   "--no-pre-mask",
                   "--dyn-heur", str(0.9),
                   # "--fix-heur", str(0.2),
                   "--preserve-rooting", "on",
                   "--filter-min-lwr", str(0.01),
                   "--outdir", output_dir,
                   '-T', str(num_threads),
                   '>', epa_files["stdout"]]
    launch_write_command(epa_command)

    # Rename the RAxML output files
    if os.path.exists(epa_info):
        copy(epa_info, epa_files["info"])
        os.remove(epa_info)
    if os.path.exists(epa_jplace):
        copy(epa_jplace, epa_files["jplace"])
        os.remove(epa_jplace)
    else:
        logging.error("Some files were not successfully created for " + query_name + "\n" +
                      "Check " + epa_files["stdout"] + " for an error!\n")
        sys.exit(3)

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


def generate_mmseqs_cluster_alignments(mmseqs_exe: str, db_name: str,
                                       clusters_file: str, align_db: str, align_tab: str) -> None:
    """
    Wrapper for generating an alignment file from MMSeqs clusters.

    :param mmseqs_exe: Path to the mmseqs executable
    :param db_name: Root path to the MMSeqs database files
    :param clusters_file: Path to the file defining MMSeqs clusters
    :param align_db: Name of the output alignment database
    :param align_tab: Name of the output alignment table
    :return: None
    """
    # Align the member sequences to each representative for a cluster
    launch_write_command([mmseqs_exe, "align", "--alignment-mode", str(3), db_name, db_name, clusters_file, align_db])

    # Convert the MMSeqs alignments to a twelve-column BLAST-like table
    stdout, mmseqs_retcode = launch_write_command([mmseqs_exe, "convertalis",
                                                   db_name, db_name, align_db, align_tab])
    if mmseqs_retcode != 0:
        logging.error("{} sequence cluster alignments failed. MMSeqs output:\n"
                      "{}\n".format(mmseqs_exe, stdout))
        sys.exit(7)
    return


def generate_mmseqs_cluster_fasta(mmseqs_exe: str, db_name: str,
                                  clusters_file: str, reps: str, reps_fasta: str) -> None:
    """
    Wrapper for generating a FASTA file containing just the representative sequences from an MMSeqs clustering analysis.

    :param mmseqs_exe: Path to the mmseqs executable
    :param db_name: Root path to the MMSeqs database files
    :param clusters_file: Path to the file defining MMSeqs clusters
    :param reps: Prefix for the file defining the cluster representatives
    :param reps_fasta: Path to the FASTA file storing the sequences of the cluster representatives
    :return: None
    """
    # Subset the database to just include the sequence for representatives
    launch_write_command([mmseqs_exe, "createsubdb", clusters_file, db_name, reps])

    # Convert the subsetted database to a FASTA file
    stdout, mmseqs_retcode = launch_write_command([mmseqs_exe, "convert2fasta", reps, reps_fasta])
    if mmseqs_retcode != 0:
        logging.error("{} failed to generate a cluster fasta file. MMSeqs output:\n"
                      "{}\n".format(mmseqs_exe, stdout))
        sys.exit(7)
    return


def generate_mmseqs_cluster_table(mmseqs_exe, db_name, int_clusters, cluster_tbl):
    stdout, mmseqs_retcode = launch_write_command([mmseqs_exe, "createtsv",
                                                   db_name, db_name,
                                                   int_clusters, cluster_tbl])
    if mmseqs_retcode != 0:
        logging.error("{} failed to generate a cluster fasta file. MMSeqs output:\n"
                      "{}\n".format(mmseqs_exe, stdout))
        sys.exit(7)
    return


def run_linclust(mmseqs_exe: str, fa_in: list, output_prefix: str, prop_sim: float,
                 aln_cov=0.95, num_threads=2, db_type="aa", tmp_dir="") -> None:
    """
    Wrapper for MMSeqs2's LinClust algorithm for clustering amino acid and nucleotide sequences in linear time.
    In addition to just generating the clusters file, this wrapper generates:
1. a FASTA file with the cluster representatives
2. a table including the sequence alignment similarities between the representatives and members

    :param mmseqs_exe: Path to the mmseqs executable
    :param fa_in: Path to the input FASTA to be clustered
    :param output_prefix: Prefix to use for naming the output files
    :param prop_sim: List matches above this sequence identity (for clustering) (range 0.0-1.0)
    :param aln_cov: List matches above this fraction of aligned (covered) residues
    :param num_threads: The number of computational threads to use when clustering
    :param db_type: Molecule type for the database, default being amino acids (aa). These are translated for MMSeqs.
    :param tmp_dir: The path to a temporary output directory. It does not need to already exist; MMSeqs2 will create it.
    :return: None
    """
    if not 0.0 < aln_cov < 1.0:
        logging.error("Provided alignment sequence similarity {} not in expected range, 0.0-1.0.\n".format(prop_sim))
        sys.exit(5)

    if not tmp_dir:
        tmp_dir = os.path.join(os.path.dirname(output_prefix), "tmp")

    # Set up MMSeqs parameters
    if db_type == "aa":
        db_type = 1
    elif db_type == "nuc":
        db_type = 2
    else:
        db_type = 0
    # Intermediate files
    db_name = output_prefix + "_mmseqs_db"
    clusters_file = db_name + "_clu"
    reps = clusters_file + "_reps"
    align_db = clusters_file + "_aln"
    # Final outputs
    align_tab = output_prefix + "_cluster_aln.tsv"
    clust_tab = output_prefix + "_cluster.tsv"
    reps_fasta = output_prefix + "_rep_seq.fasta"

    # Create the MMSeqs database
    mmseqs_db_cmd = [mmseqs_exe, "createdb",
                     "--dbtype", str(db_type),
                     ' '.join(fa_in), db_name]

    stdout, mmseqs_retcode = launch_write_command(mmseqs_db_cmd)
    if mmseqs_retcode != 0:
        logging.error("MMSeqs database creation with {} failed. Command used:\n"
                      "{}\n"
                      "MMSeqs output:\n{}\n".format(mmseqs_exe, ' '.join(mmseqs_db_cmd), stdout))
        sys.exit(7)

    # Cluster the MMSeqs database
    mmseqs_cluster_cmd = [mmseqs_exe, "linclust",
                          "--min-seq-id", str(prop_sim),
                          "-c", str(aln_cov),
                          "--threads", str(num_threads),
                          "--cluster-mode", str(2),
                          "--cov-mode", str(1),
                          db_name, clusters_file, tmp_dir]

    stdout, mmseqs_retcode = launch_write_command(mmseqs_cluster_cmd)
    if mmseqs_retcode != 0:
        logging.error("Linclust sequence clustering with {} failed."
                      " Command used:\n{}\n".format(mmseqs_exe, ' '.join(mmseqs_cluster_cmd)))
        sys.exit(7)

    # Format the desired outputs
    generate_mmseqs_cluster_fasta(mmseqs_exe, db_name, clusters_file, reps, reps_fasta)
    generate_mmseqs_cluster_alignments(mmseqs_exe, db_name, clusters_file, align_db, align_tab)
    generate_mmseqs_cluster_table(mmseqs_exe, db_name, clusters_file, clust_tab)

    # Remove the intermediate files
    rmtree(tmp_dir)
    for f_path in [reps, align_db, db_name]:
        if glob.glob(f_path + "*"):
            for f_name in glob.glob(f_path + "*"):
                if os.path.isfile(f_name):
                    os.remove(f_name)

    return


def run_vsearch_clustering(uclust_exe, fasta_input, uclust_prefix, similarity=0.60):
    """
    Wrapper function for clustering a FASTA file at some similarity using vsearch's cluster_fast algorithm

    :param uclust_exe: Path to the vsearch executable
    :param fasta_input: FASTA file for which contained sequences will be clustered
    :param uclust_prefix: Prefix for the output files
    :param similarity: The proportional similarity to cluster input sequences
    :return: None
    """
    uclust_cmd = [uclust_exe]
    uclust_cmd += ["--cluster_fast", fasta_input]
    uclust_cmd += ["--id", str(similarity)]
    uclust_cmd += ["--centroids", uclust_prefix + ".fa"]
    uclust_cmd += ["--uc", uclust_prefix + ".uc"]

    stdout, returncode = launch_write_command(uclust_cmd)

    if returncode != 0:
        logging.error("VSEARCH did not complete successfully! Command used:\n" +
                      ' '.join(uclust_cmd) + "\n")
        sys.exit(13)

    logging.debug(stdout)
    return


def cluster_sequences(software_path: str, fasta_input: str, output_prefix: str, similarity=0.60, num_threads=2) -> None:
    """
    Wrapper function for clustering a FASTA file at some similarity using some tool.

    :param software_path: Path to the software's executable
    :param fasta_input: FASTA file for which contained sequences will be clustered
    :param output_prefix: Prefix for the output files
    :param similarity: The proportional similarity to cluster input sequences
    :param num_threads: The number of threads for the clustering software to use
    :return: None
    """
    if "mmseqs" in software_path:
        logging.info("Clustering sequences with MMSeqs' Linclust... ")
        run_linclust(mmseqs_exe=software_path, fa_in=[fasta_input], output_prefix=output_prefix, prop_sim=similarity,
                     num_threads=num_threads)
    elif "vsearch" in software_path:
        logging.info("Clustering sequences with VSEARCH's cluster_fast algorithm... ")
        run_vsearch_clustering(software_path, fasta_input, output_prefix, similarity)
    else:
        logging.error("Unsupported software for clustering sequences: '{}'.\n".format(software_path))
        sys.exit(19)

    logging.info("done.\n")

    return


def build_hmm_profile(hmmbuild_exe: str, msa_in: str, output_hmm: str, name=None, graceful=False) -> None:
    """
    Wrapper for hmmbuild, used to build profile hidden Markov models (HMMs).

    :param hmmbuild_exe: Path to the hmmbuild executable
    :param msa_in: Path to the multiple sequence alignment to be used for inferring the profile HMM
    :param output_hmm: Path to the output HMM
    :param name: Optional name for the HMM. If not the name is taken from the MSA file name
    :param graceful: Optionally, TreeSAPP will not exit if hmmbuild failed and graceful is True; it will return instead.
    :return: None
    """
    logging.debug("Building HMM profile... ")
    hmm_build_command = [hmmbuild_exe]
    if name:
        hmm_build_command += ["-n", str(name)]
    hmm_build_command += [output_hmm, msa_in]
    stdout, hmmbuild_pro_returncode = launch_write_command(hmm_build_command)
    logging.debug("done.\n")

    if hmmbuild_pro_returncode != 0 and not graceful:
        logging.error("hmmbuild did not complete successfully for:\n" +
                      ' '.join(hmm_build_command) + "\n")
        sys.exit(7)
    return


def run_hmmsearch(hmmsearch_exe: str, hmm_profile: str, query_fasta: str, output_dir: str,
                  num_threads=2, e_value=1) -> list:
    """
    Function for searching a fasta file with a profile HMM

    :param hmmsearch_exe: Path to the executable for hmmsearch
    :param hmm_profile: Path to the HMM profile file
    :param query_fasta: Path to the FASTA file to be queried by the profile
    :param output_dir: Path to the directory for writing the outputs
    :param num_threads: Number of threads to be used by hmmsearch
    :param e_value: report sequences <= this E-value threshold in output
    :return: A list of domain tables created
    """
    # Find the name of the HMM. Use it to name the output file
    rp_marker = re.sub(r".hmm", '', os.path.basename(hmm_profile), flags=re.IGNORECASE)
    domtbl = output_dir + rp_marker + "_to_ORFs_domtbl.txt"

    # Basic hmmsearch command
    hmmsearch_command_base = [hmmsearch_exe]
    hmmsearch_command_base += ["--cpu", str(num_threads)]
    hmmsearch_command_base += ["-E", str(e_value)]
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


def hmmsearch_orfs(hmmsearch_exe: str, refpkg_dict: dict, fasta_file: str, output_dir: str,
                   num_threads=2, e_value=1) -> list:
    hmm_domtbl_files = list()
    nucl_target_hmm_files = list()
    prot_target_hmm_files = list()

    # Filter the HMM files to only the target markers
    for refpkg_name in refpkg_dict:
        refpkg = refpkg_dict[refpkg_name]  # type: refpkg.ReferencePackage
        if not os.path.exists(refpkg.f__search_profile):
            logging.error("Unable to locate HMM-profile for '{}'.\n".format(refpkg.prefix))
            sys.exit(3)
        else:
            if refpkg.molecule == "prot":
                prot_target_hmm_files.append(refpkg.f__search_profile)
            else:
                nucl_target_hmm_files.append(refpkg.f__search_profile)

    logging.info("Searching for marker proteins in ORFs using hmmsearch.\n")
    if logging.getLogger().disabled:
        pbar = None
    else:
        pbar = tqdm(total=len(prot_target_hmm_files) + len(nucl_target_hmm_files), ncols=120)

    # Create and launch the hmmsearch commands iteratively.
    for hmm_file in prot_target_hmm_files:
        if pbar:
            pbar.set_description("Processing {}".format(os.path.basename(hmm_file)))

        # TODO: Parallelize this by allocating no more than 2 threads per process
        hmm_domtbl_files += run_hmmsearch(hmmsearch_exe, hmm_file, fasta_file, output_dir, num_threads, e_value)

        if pbar:
            pbar.update()

    if pbar:
        pbar.close()

    return hmm_domtbl_files


def run_mafft(mafft_exe: str, fasta_in: str, fasta_out: str, num_threads) -> None:
    """
    Wrapper function for the MAFFT multiple sequence alignment tool.
    Runs MAFFT using `--auto` and checks if the output is empty.

    :param mafft_exe: Path to the executable for mafft
    :param fasta_in: An unaligned FASTA file
    :param fasta_out: The path to a file MAFFT will write aligned sequences to
    :param num_threads: Integer (or string) for the number of threads MAFFT can use
    :return:
    """
    mafft_align_command = [mafft_exe, "--auto", "--anysymbol"]
    mafft_align_command += ["--maxiterate", str(1000)]
    mafft_align_command += ["--thread", str(num_threads)]
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


def run_odseq(odseq_exe: str, fasta_in: str, outliers_fa: str, num_threads: int) -> None:
    """
    Wrapper for OD-Seq, software for detecting outliers in a multiple sequence alignment.
    1000 bootstrap replicates are used, along with a linear metric for selecting the outliers.
    These parameters were found to perform the best, empirically.

    :param odseq_exe: Path to the OD-Seq executable
    :param fasta_in: Path to the multiple sequence alignment file in FASTA format.
    :param outliers_fa: Path to the fasta file containing the outliers detected
    :param num_threads: Number of threads for OD-Seq to use
    :return: None
    """
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

    trimmed_msa_file = '.'.join(mfa_file.split('.')[:-1]) + '-' + re.escape(tool) + ".fasta"
    if tool == "trimAl":
        trim_command = trimal_command(executables["trimal"], mfa_file, trimmed_msa_file)
    elif tool == "BMGE":
        trim_command = bmge_command(executables["BMGE.jar"], mfa_file, trimmed_msa_file, molecule)
    else:
        logging.error("Unsupported trimming software requested: '" + tool + "'")
        sys.exit(5)

    return trim_command, trimmed_msa_file


def filter_multiple_alignments(executables, concatenated_mfa_files, refpkg_dict, n_proc=1, tool="BMGE"):
    """
    Runs BMGE using the provided lists of the concatenated hmmalign files, and the number of sequences in each file.

    :param executables: A dictionary mapping software to a path of their respective executable
    :param concatenated_mfa_files: A dictionary containing f_contig keys mapping to a FASTA or Phylip sequential file
    :param refpkg_dict: A dictionary of ReferencePackage instances indexed by their respective denominators
    :param n_proc: The number of parallel processes to be launched for alignment trimming
    :param tool: The software to use for alignment trimming
    :return: A list of files resulting from BMGE multiple sequence alignment masking.
    """
    logging.info("Running " + tool + "... ")

    start_time = time.time()
    task_list = list()
    trimmed_output_files = {}

    for refpkg_code in sorted(concatenated_mfa_files.keys()):
        if refpkg_code not in trimmed_output_files:
            trimmed_output_files[refpkg_code] = []
        mfa_files = concatenated_mfa_files[refpkg_code]
        for concatenated_mfa_file in mfa_files:
            trim_command, trimmed_msa_file = get_msa_trim_command(executables, concatenated_mfa_file,
                                                                  refpkg_dict[refpkg_code].molecule, tool)
            trimmed_output_files[refpkg_code].append(trimmed_msa_file)
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
