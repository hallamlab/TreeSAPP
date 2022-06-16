import sys
import os.path
import time
import logging

from treesapp import logger
from treesapp import refpkg
from treesapp import external_command_interface as eci
from treesapp import clipkit_helper as ckh

LOGGER = logging.getLogger(logger.logger_name())


def trim_multiple_alignment_clipkit(msa_file: str, ref_pkg: refpkg.ReferencePackage,
                                    min_seq_length: int) -> ckh.ClipKitHelper:
    # Modes can be one of 'smart-gap', 'kpi', 'kpic', 'gappy', 'kpi-smart-gap', 'kpi-gappy'
    trimmer = ckh.ClipKitHelper(fasta_in=msa_file,
                                output_dir=os.path.dirname(msa_file),
                                mode="gappy",
                                gap_prop=0.25)
    trimmer.refpkg_name = ref_pkg.prefix
    trimmer.run()
    trimmer.compare_original_and_trimmed_multiple_alignments(min_seq_length, ref_pkg)
    trimmer.write_qc_trimmed_multiple_alignment()
    return trimmer


def summarise_trimming(msa_trimmers: list) -> None:
    """Summarises various outcomes of trimming MSAs."""
    refpkg_trimming_stats = {trimmer.refpkg_name: {
        "msa_files": 0,
        "cols_removed": [],
        "seqs_removed": [],
        "successes": 0,
    }
        for trimmer in msa_trimmers}
    num_successful_alignments = 0

    LOGGER.debug("Validating trimmed multiple sequence alignment files... ")
    for trimmer in msa_trimmers:  # type: ckh.ClipKitHelper
        # Gather all useful stats for each trimmer instance
        refpkg_trimming_stats[trimmer.refpkg_name]["msa_files"] += 1
        if trimmer.success:
            refpkg_trimming_stats[trimmer.refpkg_name]["successes"] += 1
            num_successful_alignments += 1
        else:
            continue
        refpkg_trimming_stats[trimmer.refpkg_name]["cols_removed"].append(trimmer.num_msa_cols - trimmer.num_trim_cols)
        refpkg_trimming_stats[trimmer.refpkg_name]["seqs_removed"].append(trimmer.num_msa_seqs - trimmer.num_trim_seqs)

    # Summarise trimming by reference package
    for refpkg_name, stats in refpkg_trimming_stats.items():
        trim_summary = "\t\t{} trimming stats:\n".format(refpkg_name)
        if stats["msa_files"] == 0:
            continue
        trim_summary += "Multiple alignment files   = {}\n".format(stats["msa_files"])
        trim_summary += "Files successfully trimmed = {}\n".format(stats["successes"])
        trim_summary += "Average columns removed    = {}\n".format(round(sum(stats["cols_removed"]) /
                                                                         len(stats["cols_removed"])))
        trim_summary += "Average sequences removed  = {}\n".format(round(sum(stats["seqs_removed"]) /
                                                                         len(stats["seqs_removed"])))

        LOGGER.debug(trim_summary + "\n")

    LOGGER.debug("done.\n")

    if num_successful_alignments == 0:
        LOGGER.error("No quality alignment files to analyze after trimming. Exiting now.\n")
        sys.exit(0)  # Should be 3, but this allows Clade_exclusion_analyzer to continue after exit
    return


def gather_multiple_alignments(msa_trimmers: list) -> dict:
    """
    Creates a dictionary of MSA files indexed by reference package names.
    These files are trimmed outputs if trimming was successful, or the original if not.
    """
    trimmed_output_files = {}
    for trimmer in msa_trimmers:  # type: ckh.ClipKitHelper
        try:
            trimmed_output_files[trimmer.refpkg_name].append(trimmer.get_qc_output())
        except KeyError:
            trimmed_output_files[trimmer.refpkg_name] = [trimmer.get_qc_output()]
    return trimmed_output_files


def trim_multiple_alignment_farmer(concatenated_mfa_files: dict, min_seq_length: int, ref_pkgs: dict,
                                   n_proc=1, silent=False) -> dict:
    """
    Runs ClipKit using the provided lists of the concatenated hmmalign files, and the number of sequences in each file.

    :param concatenated_mfa_files: A dictionary containing f_contig keys mapping to a FASTA or Phylip sequential file
    :param min_seq_length: Minimum length for a sequence to be retained in the MSA
    :param ref_pkgs: A dictionary of reference package names mapped to ReferencePackage instances
    :param n_proc: The number of parallel processes to be launched for alignment trimming
    :param silent: A boolean indicating whether the
    :return: A list of files resulting from multiple sequence alignment masking.
    """
    start_time = time.time()
    task_list = list()

    for refpkg_code, mfa_files in sorted(concatenated_mfa_files.items()):
        for msa in mfa_files:
            task_list.append({"msa_file": msa,
                              "ref_pkg": ref_pkgs[refpkg_code],
                              "min_seq_length": min_seq_length})

    msa_trimmers = eci.run_apply_async_multiprocessing(func=trim_multiple_alignment_clipkit,
                                                       arguments_list=task_list,
                                                       num_processes=n_proc,
                                                       pbar_desc="Multiple alignment trimming",
                                                       disable=silent)

    end_time = time.time()
    hours, remainder = divmod(end_time - start_time, 3600)
    minutes, seconds = divmod(remainder, 60)
    LOGGER.debug("\tMultiple alignment trimming time required: " +
                 ':'.join([str(hours), str(minutes), str(round(seconds, 2))]) + "\n")

    summarise_trimming(msa_trimmers)
    # Collect the trimmed (or untrimmed if reference sequences were removed) output files
    trimmed_output_files = gather_multiple_alignments(msa_trimmers)

    return trimmed_output_files
