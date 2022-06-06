import os.path
import time
import logging

from treesapp import logger
from treesapp import refpkg
from treesapp import external_command_interface as eci
from treesapp import clipkit_helper as ckh


LOGGER = logging.getLogger(logger.logger_name())


def trim_multiple_alignment_clipkit(msa_file: str, ref_pkg: refpkg.ReferencePackage, min_seq_length: int) -> ckh.ClipKitHelper:
    trimmer = ckh.ClipKitHelper(fasta_in=msa_file,
                                output_dir=os.path.dirname(msa_file))
    trimmer.refpkg_name = ref_pkg.prefix
    trimmer.run()
    trimmer.compare_original_and_trimmed_multiple_alignments(min_seq_length, ref_pkg)
    trimmer.write_qc_trimmed_multiple_alignment()
    return trimmer


def summarise_trimming(msa_trimmers: list) -> None:
    """Summarises various outcomes of trimming MSAs."""
    num_successful_alignments = 0
    discarded_seqs_string = ""
    trimmed_away_seqs = dict()
    untrimmed_msa_failed = []
    LOGGER.debug("Validating trimmed multiple sequence alignment files... ")
    for trimmer in msa_trimmers:  # type: ckh.ClipKitHelper
        # TODO: Gather all useful stats for each trimmer instance
        if trimmer.success:
            num_successful_alignments += 1

    # TODO: Summarise trimming by reference package
    trimming_performance_string = "\tAverage columns removed:\n"
    for refpkg_name in trimmed_length_dict:
        trimming_performance_string += "\t\t" + refpkg_name + "\t"
        n_trimmed_files = len(trimmed_length_dict[denominator])
        if n_trimmed_files > 0:
            trimming_performance_string += str(
                round(sum(trimmed_length_dict[denominator]) / n_trimmed_files, 1)) + "\n"
        else:
            trimming_performance_string += str(0.0) + "\n"

    LOGGER.debug(trimming_performance_string + "\n")

    discarded_seqs_string += "\n\t\t" + self.mfa_out + " = " + str(len(discarded_seqs))
    num_successful_alignments += len(msa_passed)
    qc_ma_dict[ref_pkg.prefix] = msa_passed
    discarded_seqs_string += summary_str
    untrimmed_msa_failed.clear()

    LOGGER.debug("done.\n")
    LOGGER.debug("\tSequences removed during trimming:\n\t\t" +
                 '\n\t\t'.join([k + ": " + str(trimmed_away_seqs[k]) for k in trimmed_away_seqs.keys()]) + "\n")

    LOGGER.debug("\tSequences <" + str(min_len) + " characters removed after trimming:" +
                 discarded_seqs_string + "\n")

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
