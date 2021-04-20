__author__ = "Connor Morgan-Lang"
__maintainer__ = "Connor Morgan-Lang"
__license__ = "GPL-3.0"

import sys
import os
import shutil
import re
import glob
import time
import logging
from collections import namedtuple

import pyfastx
import pandas as pd
from numpy import array as np_array

from treesapp import abundance
from treesapp import classy
from treesapp import phylo_seq
from treesapp import refpkg
from treesapp import treesapp_args
from treesapp.entish import index_tree_edges, map_internal_nodes_leaves
from treesapp import lca_calculations as ts_lca
from treesapp import jplace_utils
from treesapp import file_parsers
from treesapp import phylo_dist
from treesapp import utilities
from treesapp import wrapper
from treesapp import fasta
from treesapp import training_utils
from treesapp.hmmer_tbl_parser import HmmMatch


class Assigner(classy.TreeSAPP):
    def __init__(self):
        """Instantiate a class instance used by `treesapp assign`."""
        super(Assigner, self).__init__("assign")
        self.reference_tree = None
        self.svc_filter = False
        self.fasta_full_name = True
        self.aa_orfs_file = ""
        self.nuc_orfs_file = ""
        self.classified_aa_seqs = ""
        self.classified_nuc_seqs = ""
        self.classification_table = ""
        self.itol_out = ""
        self.composition = ""
        self.target_refpkgs = list()

        self.stages = {0: classy.ModuleFunction("orf-call", 0, self.predict_orfs),
                       1: classy.ModuleFunction("clean", 1, self.clean),
                       2: classy.ModuleFunction("search", 2, self.search),
                       3: classy.ModuleFunction("align", 3, self.align),
                       4: classy.ModuleFunction("place", 4, self.place),
                       5: classy.ModuleFunction("classify", 5, self.classify),
                       6: classy.ModuleFunction("abundance", 6)}
        self.current_stage = self.stages[0]
        return

    def check_classify_arguments(self, args):
        """
        Ensures the command-line arguments returned by argparse are sensible.

        :param args: object with parameters returned by argparse.parse_args()
        :return: 'args', a summary of TreeSAPP settings.
        """
        if not os.path.isfile(self.input_sequences):
            logging.error("FASTX input file '{}' doesn't exist.\n".format(self.input_sequences))
            sys.exit(5)

        self.find_sequence_molecule_type()
        self.classification_table = self.final_output_dir + os.sep + self.classification_tbl_name
        self.itol_out = self.output_dir + 'iTOL_output' + os.sep
        self.classified_aa_seqs = self.final_output_dir + self.sample_prefix + "_classified.faa"
        self.classified_nuc_seqs = self.final_output_dir + self.sample_prefix + "_classified.fna"

        if args.targets:
            self.target_refpkgs = args.targets.split(',')
        else:
            self.target_refpkgs = []

        self.validate_refpkg_dir(args.refpkg_dir)

        if self.molecule_type == "prot":
            self.query_sequences = self.input_sequences
            self.change_stage_status("orf-call", False)
            if args.rel_abund:
                logging.error("Unable to calculate abundance values for protein sequences.\n")
                sys.exit(3)

        self.formatted_input = self.stage_lookup("clean").dir_path + self.sample_prefix + "_formatted.fasta"

        if args.svm:
            self.svc_filter = True

        return args

    def load_refpkg_classifiers(self, refpkg_dict: dict, kernel: str, threads: int, combine=False) -> None:
        if not self.svc_filter:
            return

        untrained_refpkgs = []
        classifiers = dict()
        if combine:
            training_frames = []
            refpkg_names = set()
            for name, ref_pkg in refpkg_dict.items():  # type: refpkg.ReferencePackage
                refpkg_names.add(name)
                if len(ref_pkg.training_df) == 0:
                    continue
                training_frames.append(ref_pkg.training_df)
            if len(training_frames) == 0:
                logging.error("All reference package training data frames are empty.\n"
                              "Unable to train a classifier from combined training data.\n")
                sys.exit(5)

            classifiers = training_utils.train_classifier_from_dataframe(training_df=pd.concat(training_frames),
                                                                         kernel=kernel,
                                                                         num_threads=threads,
                                                                         refpkg_names=refpkg_names)
        else:
            for _name, ref_pkg in refpkg_dict.items():  # type: refpkg.ReferencePackage
                if len(ref_pkg.training_df) == 0:
                    continue
                classifiers.update(training_utils.train_classifier_from_dataframe(training_df=ref_pkg.training_df,
                                                                                  kernel=kernel,
                                                                                  num_threads=threads))

        # Set the svc attribute of each ReferencePackage to its respective classifier
        for _name, ref_pkg in refpkg_dict.items():
            try:
                ref_pkg.svc = classifiers[ref_pkg.prefix]
            except KeyError:
                if not ref_pkg.svc:
                    untrained_refpkgs.append(ref_pkg.prefix)

        if untrained_refpkgs:
            logging.warning("Unable to train classifiers for {} reference packages.\n".format(len(untrained_refpkgs)))
            logging.debug("Reference packages that will not use SVC for filtering placements:\n\t"
                          "{}\n".format("\n\t".join(untrained_refpkgs)))

        return


    @staticmethod
    def define_hmm_domtbl_thresholds(args):
        thresholds_nt = namedtuple("thresholds", ["perc_aligned", "min_acc", "max_e", "max_ie", "min_score"])

        # Parameterizing the hmmsearch output parsing:
        if args.stringency == "relaxed":
            domtbl_thresholds = thresholds_nt(perc_aligned=args.hmm_coverage,
                                              min_acc=0.7, max_e=1E-3, max_ie=1E-1, min_score=15)
        elif args.stringency == "strict":
            domtbl_thresholds = thresholds_nt(perc_aligned=args.hmm_coverage,
                                              min_acc=0.7, max_e=1E-5, max_ie=1E-3, min_score=30)
        else:
            logging.error("Unknown HMM-parsing stringency argument '" + args.stringency + "'.\n")
            sys.exit(3)
        return domtbl_thresholds

    def decide_stage(self, args) -> None:
        """
        Bases the stage(s) to run on args.stage which is broadly set to either 'continue' or any other valid stage

        This function ensures all the required inputs are present for beginning at the desired first stage,
        otherwise, the pipeline begins at the first possible stage to continue and ends once the desired stage is done.

        :return: None
        """
        self.set_stage_dir()
        if self.stage_status("orf-call"):
            self.fasta_full_name = False
            self.aa_orfs_file = self.stage_output_dir + self.sample_prefix + "_ORFs.faa"
            self.nuc_orfs_file = self.stage_output_dir + self.sample_prefix + "_ORFs.fna"
            if os.path.isfile(self.aa_orfs_file) and os.path.isfile(self.nuc_orfs_file):
                self.change_stage_status("orf-call", False)
                self.query_sequences = self.aa_orfs_file

        if args.rel_abund and self.molecule_type == "dna":
            if not args.reads:
                if args.reverse:
                    logging.error("File containing reverse reads provided but forward mates file missing!\n")
                    sys.exit(3)
                else:
                    logging.error("At least one FASTQ file must be provided if --rel_abund flag is active!\n")
                    sys.exit(3)
            for fp in args.reads:
                if not os.path.isfile(fp):
                    logging.error("Path to forward reads ('%s') doesn't exist.\n" % fp)
                    sys.exit(3)
            for fp in args.reverse:
                if not os.path.isfile(fp):
                    logging.error("Path to reverse reads ('%s') doesn't exist.\n" % fp)
                    sys.exit(3)
            self.change_stage_status("abundance", True)
        else:
            self.change_stage_status("abundance", False)

        self.validate_continue(args)
        return

    def get_info(self):
        info_string = "Assigner instance summary:\n"
        info_string += super(Assigner, self).get_info() + "\n\t"
        info_string += "\n\t".join(["ORF protein sequences = " + self.aa_orfs_file,
                                    "Target reference packages = " + str(self.target_refpkgs),
                                    "Composition of input = " + self.composition]) + "\n"

        return info_string

    def predict_orfs(self, composition: str, num_threads: int) -> None:
        """
        Predict ORFs from the input FASTA file using Prodigal

        :param composition: Sample composition being either a single organism or a metagenome [single | meta]
        :param num_threads: The number of CPU threads to use
        :return: None
        """

        logging.info("Predicting open-reading frames using Prodigal... ")

        start_time = time.time()

        if num_threads > 1 and composition == "meta":
            # Split the input FASTA into num_threads files to run Prodigal in parallel
            split_files = fasta.fastx_split(self.input_sequences, self.stage_output_dir, num_threads)
        else:
            split_files = [self.input_sequences]

        task_list = list()
        for fasta_chunk in split_files:
            chunk_prefix = self.stage_output_dir + '.'.join(os.path.basename(fasta_chunk).split('.')[:-1])
            prodigal_command = [self.executables["prodigal"], "-q"]
            prodigal_command += ["-i", fasta_chunk]
            prodigal_command += ["-p", composition]
            prodigal_command += ["-a", chunk_prefix + "_ORFs.faa"]
            prodigal_command += ["-d", chunk_prefix + "_ORFs.fna"]
            task_list.append(prodigal_command)

        num_tasks = len(task_list)
        if num_tasks > 0:
            cl_farmer = wrapper.CommandLineFarmer("Prodigal -p " + composition, num_threads)
            cl_farmer.add_tasks_to_queue(task_list)

            cl_farmer.task_queue.close()
            cl_farmer.task_queue.join()

        tmp_prodigal_aa_orfs = glob.glob(self.stage_output_dir + self.sample_prefix + "*_ORFs.faa")
        tmp_prodigal_nuc_orfs = glob.glob(self.stage_output_dir + self.sample_prefix + "*_ORFs.fna")
        if not tmp_prodigal_aa_orfs or not tmp_prodigal_nuc_orfs:
            logging.error("Prodigal outputs were not generated:\n"
                          "Amino acid ORFs: " + ", ".join(tmp_prodigal_aa_orfs) + "\n" +
                          "Nucleotide ORFs: " + ", ".join(tmp_prodigal_nuc_orfs) + "\n")
            sys.exit(5)

        # Concatenate outputs
        if not os.path.isfile(self.aa_orfs_file) and not os.path.isfile(self.nuc_orfs_file):
            aa_cat_tmp = self.aa_orfs_file + ".tmp"
            nuc_cat_tmp = self.nuc_orfs_file + ".tmp"
            tmp_cat_files = [aa_cat_tmp, nuc_cat_tmp]
            utilities.concatenate_files(tmp_prodigal_aa_orfs, aa_cat_tmp)
            utilities.concatenate_files(tmp_prodigal_nuc_orfs, nuc_cat_tmp)

            # Remove Prodigal header tags from the ORF FASTA files
            fasta.format_fasta(fasta_input=aa_cat_tmp, output_fasta=self.aa_orfs_file,
                               molecule="prot", true_name=True, full_name=self.fasta_full_name)
            fasta.format_fasta(fasta_input=nuc_cat_tmp, output_fasta=self.nuc_orfs_file,
                               molecule="dna", true_name=True, full_name=self.fasta_full_name)

            # Remove intermediate files
            intermediate_files = list(tmp_prodigal_aa_orfs + tmp_prodigal_nuc_orfs + split_files + tmp_cat_files)
            for tmp_file in intermediate_files:
                if tmp_file != self.input_sequences:
                    os.remove(tmp_file)

        logging.info("done.\n")

        end_time = time.time()
        hours, remainder = divmod(end_time - start_time, 3600)
        minutes, seconds = divmod(remainder, 60)
        logging.debug("\tProdigal time required: " +
                      ':'.join([str(hours), str(minutes), str(round(seconds, 2))]) + "\n")

        self.query_sequences = self.aa_orfs_file
        # self.change_stage_status("clean", False)
        self.increment_stage_dir()
        return

    def clean(self):
        return

    def search(self):
        return

    def align(self):
        return

    def place(self):
        return

    def classify(self):
        return

    def write_classified_orfs(self, pqueries: dict, extracted_seqs: dict) -> None:
        # Write the amino acid sequences
        write_classified_sequences(pqueries, extracted_seqs, self.classified_aa_seqs)
        # Write the nucleotide sequences
        if self.molecule_type == "dna":
            # Create the set of classified sequence names to subset the fasta - prevent reading in entire file
            classified_seq_names = set()
            for refpkg_name in pqueries:
                classified_seq_names.update({pq.seq_name for pq in pqueries[refpkg_name]})
            if os.path.isfile(self.nuc_orfs_file):
                nuc_orfs = fasta.FASTA(self.nuc_orfs_file)
                nuc_orfs.fasta_dict = fasta.format_read_fasta(self.nuc_orfs_file,
                                                              molecule="dna",
                                                              subset=classified_seq_names,
                                                              full_name=self.fasta_full_name)
                nuc_orfs.header_registry = fasta.register_headers(nuc_orfs.fasta_dict.keys())
                nuc_orfs.change_dict_keys()
                if not os.path.isfile(self.classified_nuc_seqs):
                    logging.info("Creating nucleotide FASTA file of classified sequences '{}'... "
                                 "".format(self.classified_nuc_seqs))
                    write_classified_sequences(pqueries, nuc_orfs.fasta_dict, self.classified_nuc_seqs)
                    logging.info("done.\n")
            else:
                logging.warning("Unable to read '" + self.nuc_orfs_file + "'.\n" +
                                "Cannot create the nucleotide FASTA file of classified sequences!\n")
        return

    def fetch_hmmsearch_outputs(self, target_refpkg_prefixes: set) -> list:
        if self.current_stage.name != "search":
            logging.error("Unable to fetch hmmsearch outputs as the current stage ({}) is incorrect.\n"
                          "".format(self.current_stage.name))
            sys.exit(3)

        # Collect all files from output directory
        hmm_domtbls = glob.glob(self.stage_output_dir + "*_search_to_ORFs_domtbl.txt")

        # Ensure all of the domain tables are present compared to the reference packages
        searched = set([os.path.basename(f_path).split('_')[0] for f_path in hmm_domtbls])
        if not searched or target_refpkg_prefixes.difference(searched):
            return []
        else:
            return hmm_domtbls


def replace_contig_names(numeric_contig_index: dict, fasta_obj: fasta.FASTA):
    for marker in numeric_contig_index:
        assign_re = re.compile(r"(.*)\|{0}\|(\d+_\d+)$".format(marker))
        for neg_num_id in numeric_contig_index[marker]:
            assign_name = numeric_contig_index[marker][neg_num_id]
            seq_name, coords = assign_re.match(assign_name).groups()
            try:
                original_name = fasta_obj.header_registry[seq_name].original
            except KeyError:
                logging.error("Unable to find TreeSAPP numerical ID '" + seq_name + "' in header registry.\n")
                sys.exit(3)
            numeric_contig_index[marker][neg_num_id] = original_name + '|' + marker + '|' + coords
    return numeric_contig_index


def load_pqueries(hmm_matches: dict, query_seq_fasta: fasta.FASTA) -> list:
    logging.debug("Instantiating the PQuery instances... ")

    pqueries = []
    query_seq_fasta.change_dict_keys("num")
    for refpkg_name, refpkg_matches in hmm_matches.items():  # type: (str, list)
        for hmm_match in refpkg_matches:  # type: HmmMatch
            if hmm_match.desc != '-':
                seq_name = hmm_match.orf + ' ' + hmm_match.desc
            else:
                seq_name = hmm_match.orf
            # Load the homology search data
            qseq = phylo_seq.PQuery()
            qseq.ref_name = refpkg_name
            qseq.evalue, qseq.start, qseq.end = hmm_match.eval, hmm_match.start, hmm_match.end
            pqueries.append(qseq)
            # Load the query's sequence
            qseq.seq = query_seq_fasta.fasta_dict[seq_name]
            header = query_seq_fasta.header_registry[seq_name]  # type: fasta.Header
            qseq.seq_name = header.original
            qseq.place_name = "{}|{}|{}_{}".format(qseq.seq_name, qseq.ref_name, qseq.start, qseq.end)

    logging.debug("done.\n")

    return pqueries


def load_homologs(hmm_matches: dict, hmmsearch_query_fasta: str, query_seq_fasta: fasta.FASTA) -> None:
    """
    Loads the fasta_dict attribute in query_seq_fasta guided by the homologous sequences.
    This reduces the RAM usage when the query FASTA contains many non-homologous sequences, which would not be loaded.
    Homologous sequences are the names of sequences in HmmMatch instances stored by the dictionary hmm_matches.

    :param hmm_matches: A dictionary of lists of HmmMatch instances, indexed by reference package names
    :param hmmsearch_query_fasta: Path to a FASTA file that was used as a query against the profile HMM(s),
    reflecting those alignments stored in hmm_matches.
    :param query_seq_fasta: A FASTA instance, typically with an empty fasta_dict attribute
    :return: None
    """
    logging.info("Loading homologous sequences identified... ")
    # Create a set of sequence names that matched a profile HMM
    matched_query_names = set()
    query_seq_fasta.fasta_dict.clear()
    for refpkg_name, refpkg_matches in hmm_matches.items():  # type: (str, list)
        for hmm_match in refpkg_matches:  # type: HmmMatch
            if hmm_match.desc != '-':
                seq_name = hmm_match.orf + ' ' + hmm_match.desc
            else:
                seq_name = hmm_match.orf
            matched_query_names.add(seq_name)

    # Load just homologous sequences into the FASTA.fasta_dict
    for name, seq in pyfastx.Fasta(hmmsearch_query_fasta, build_index=False):  # type: (str, str)
        if name in matched_query_names:
            query_seq_fasta.fasta_dict[name] = seq

    # Keep only the homologous sequences in FASTA.header_registry
    query_seq_fasta.synchronize_seqs_n_headers()

    logging.info("done.\n")
    return


def bin_hmm_matches(hmm_matches: dict, fasta_dict: dict) -> (dict, dict):
    """
    Used for extracting query sequences that mapped to reference package HMM profiles. These are binned into groups
    based on the location on the HMM profile they mapped to such that MSAs downstream will have more conserved positions

    The first nested dictionary returned "extracted_seq_dict" contains marker (i.e. ref_pkg names) strings mapped
    to bin numbers mapped to query sequence negative integer code names mapped to their extracted, or sliced, sequence.

    The second dictionary returned "numeric_contig_index" is used for mapping query sequence negative integer code names
    mapped to their original header names with the alignment coordinates appended at the end for each marker.

    :param hmm_matches: Contains lists of HmmMatch objects mapped to the marker they matched
    :param fasta_dict: Stores either the original or ORF-predicted input FASTA. Headers are keys, sequences are values
    :return: List of files that go on to placement stage, dictionary mapping marker-specific numbers to contig names
    """
    logging.info("Extracting and grouping the quality-controlled sequences... ")
    extracted_seq_dict = dict()  # Keys are markers -> bin_num -> negative integers -> extracted sequences
    numeric_contig_index = dict()  # Keys are markers -> negative integers -> headers
    bins = dict()

    for marker in hmm_matches:
        if len(hmm_matches[marker]) == 0:
            continue
        if marker not in numeric_contig_index.keys():
            numeric_contig_index[marker] = dict()
        numeric_decrementor = -1
        if marker not in extracted_seq_dict:
            extracted_seq_dict[marker] = dict()

        # Algorithm for binning sequences:
        # 1. Sort HmmMatches by the proportion of the HMM profile they covered in increasing order (full-length last)
        # 2. For HmmMatch in sorted matches, determine overlap between HmmMatch and each bin's representative HmmMatch
        # 3. If overlap exceeds 80% of representative's aligned length add it to the bin, else continue
        # 4. When bins are exhausted create new bin with HmmMatch
        for hmm_match in sorted(hmm_matches[marker], key=lambda x: x.end - x.start):  # type: HmmMatch
            if hmm_match.desc != '-':
                contig_name = hmm_match.orf + ' ' + hmm_match.desc
            else:
                contig_name = hmm_match.orf
            # Add the query sequence to the index map
            orf_coordinates = str(hmm_match.start) + '_' + str(hmm_match.end)
            numeric_contig_index[marker][numeric_decrementor] = contig_name + '|' + marker + '|' + orf_coordinates
            # Add the FASTA record of the trimmed sequence - this one moves on for placement
            full_sequence = fasta_dict[contig_name]
            binned = False
            for bin_num in sorted(bins):
                bin_rep = bins[bin_num][0]
                overlap = min(hmm_match.pend, bin_rep.pend) - max(hmm_match.pstart, bin_rep.pstart)
                if (100*overlap)/(bin_rep.pend - bin_rep.pstart) > 80:  # 80 refers to overlap proportion with seed
                    bins[bin_num].append(hmm_match)
                    extracted_seq_dict[marker][bin_num][numeric_decrementor] = full_sequence[
                                                                               hmm_match.start - 1:hmm_match.end]
                    binned = True
                    break
            if not binned:
                bin_num = len(bins)
                bins[bin_num] = list()
                extracted_seq_dict[marker][bin_num] = dict()
                bins[bin_num].append(hmm_match)
                extracted_seq_dict[marker][bin_num][numeric_decrementor] = full_sequence[
                                                                           hmm_match.start - 1:hmm_match.end]
            numeric_decrementor -= 1

        bins.clear()
    logging.info("done.\n")

    return extracted_seq_dict, numeric_contig_index


def write_grouped_fastas(extracted_seq_dict: dict, numeric_contig_index: dict, refpkg_dict: dict, output_dir: str):
    hmmalign_input_fastas = list()
    bulk_marker_fasta = dict()
    bin_fasta = dict()

    group_size_string = "Number of query sequences in each marker's group:\n"
    for marker in extracted_seq_dict:
        for group in sorted(extracted_seq_dict[marker]):
            if extracted_seq_dict[marker][group]:
                group_size_string += "\t".join([marker, str(group), str(len(extracted_seq_dict[marker][group]))]) + "\n"
    logging.debug(group_size_string + "\n")

    logging.info("Writing the grouped sequences to FASTA files... ")

    for marker in extracted_seq_dict:
        ref_pkg = refpkg_dict[marker]  # type: refpkg.ReferencePackage
        f_acc = 0  # For counting the number of files for a marker. Will exceed groups if len(queries) > len(references)
        for group in sorted(extracted_seq_dict[marker]):
            if extracted_seq_dict[marker][group]:
                group_sequences = extracted_seq_dict[marker][group]
                for num in group_sequences:
                    # Add the query sequence to the master marker FASTA with the full sequence name
                    bulk_marker_fasta[numeric_contig_index[marker][num]] = group_sequences[num]

                    # Add the query sequence to this bin's FASTA file
                    bin_fasta[str(num)] = group_sequences[num]
                    # Ensuring the number of query sequences doesn't exceed the number of reference sequences
                    if len(bin_fasta) >= ref_pkg.num_seqs:
                        fasta.write_new_fasta(bin_fasta, output_dir + marker + "_hmm_purified_group" + str(f_acc) + ".faa")
                        hmmalign_input_fastas.append(output_dir + marker + "_hmm_purified_group" + str(f_acc) + ".faa")
                        f_acc += 1
                        bin_fasta.clear()
                if len(bin_fasta) >= 1:
                    fasta.write_new_fasta(bin_fasta, output_dir + marker + "_hmm_purified_group" + str(f_acc) + ".faa")
                    hmmalign_input_fastas.append(output_dir + marker + "_hmm_purified_group" + str(f_acc) + ".faa")
            f_acc += 1
            bin_fasta.clear()

        # Now write a single FASTA file with all identified markers
        if len(bulk_marker_fasta) >= 1:
            trimmed_hits_fasta = output_dir + marker + "_hmm_purified.faa"
            fasta.write_new_fasta(bulk_marker_fasta, trimmed_hits_fasta)
        bulk_marker_fasta.clear()
    logging.info("done.\n")
    return hmmalign_input_fastas


def subsequence(fasta_dictionary, contig_name, start, end):
    """
    Extracts a sub-sequence from `start` to `end` of `contig_name` in `fasta_dictionary`
     with headers for keys and sequences as values. `contig_name` does not contain the '>' character

    :param fasta_dictionary:
    :param contig_name:
    :param start:
    :param end:
    :return: A string representing the sub-sequence of interest
    """
    subseq = fasta_dictionary['>' + contig_name][start:end]
    return subseq


def get_sequence_counts(concatenated_mfa_files: dict, ref_alignment_dimensions: dict, verbosity: bool, file_type: str):
    alignment_length_dict = dict()
    for refpkg_name in concatenated_mfa_files:
        if refpkg_name not in ref_alignment_dimensions:
            logging.error("Unrecognized code '" + refpkg_name + "'.\n")
            sys.exit(3)

        ref_n_seqs, ref_seq_length = ref_alignment_dimensions[refpkg_name]
        for msa_file in concatenated_mfa_files[refpkg_name]:
            if file_type == "Fasta":
                seq_dict = fasta.read_fasta_to_dict(msa_file)
            elif file_type == "Phylip":
                seq_dict = file_parsers.read_phylip_to_dict(msa_file)
            elif file_type == "Stockholm":
                seq_dict = file_parsers.read_stockholm_to_dict(msa_file)
            else:
                logging.error("File type '" + file_type + "' is not recognized.")
                sys.exit(3)
            num_seqs, sequence_length = fasta.multiple_alignment_dimensions(msa_file, seq_dict)
            alignment_length_dict[msa_file] = sequence_length

            # Warn user if the multiple sequence alignment has grown significantly
            if verbosity and ref_seq_length*1.5 < sequence_length:
                logging.warning("Multiple alignment of '{}' caused >150% increase in the number of columns"
                                " ({} -> {}).\n".format(refpkg_name, ref_seq_length, sequence_length))
    return alignment_length_dict


def prep_reference_packages_for_assign(refpkg_dict: dict, output_dir: str) -> None:
    """
    Write the individual reference package files for each ReferencePackage instance in refpkg_dict

    :param refpkg_dict: A dictionary of ReferencePackage instances indexed by their prefix values
    :param output_dir: A main directory to write the files for each ReferencePackage
    :return: None
    """
    for refpkg_name in refpkg_dict:
        ref_pkg = refpkg_dict[refpkg_name]  # type: refpkg.ReferencePackage
        prefix_dir = os.path.join(output_dir, ref_pkg.prefix + "_RefPkg")
        ref_pkg.disband(prefix_dir)
        ref_pkg.pickle_package(new_output_dir=prefix_dir)
    return


def get_alignment_dims(refpkg_dict: dict):
    """
    Used for collecting the dimensions of multiple sequence alignment (MSA) files for ReferencePackage instances

    :param refpkg_dict: A dictionary of ReferencePackage.prefix keys mapped to ReferencePackage instances
    :return: A dictionary of ReferencePackage.prefix keys and (nrow, ncolumn) values
    """
    alignment_dimensions_dict = dict()
    for refpkg_name in refpkg_dict:  # type: str
        ref_pkg = refpkg_dict[refpkg_name]  # type: refpkg.ReferencePackage
        alignment_dimensions_dict[ref_pkg.prefix] = ref_pkg.alignment_dims()

    return alignment_dimensions_dict


def multiple_alignments(executables: dict, single_query_sequence_files: list, refpkg_dict: dict,
                        tool="hmmalign", output_dir="", num_proc=4) -> dict:
    """
    Wrapper function for the multiple alignment functions - only purpose is to make an easy decision at this point...

    :param executables: Dictionary mapping software names to their executables
    :param single_query_sequence_files: List of unaligned query sequences in FASTA format
    :param refpkg_dict: A dictionary of ReferencePackage instances indexed by their respective prefix
    :param tool: Tool to use for aligning query sequences to a reference multiple alignment [hmmalign|papara]
    :param output_dir: Path to write the new MSA files
    :param num_proc: The number of alignment jobs to run in parallel
    :return: Dictionary of multiple sequence alignment (FASTA) files indexed by denominator
    """
    if tool == "hmmalign":
        concatenated_msa_files = prepare_and_run_hmmalign(executables, single_query_sequence_files, refpkg_dict,
                                                          output_dir, num_proc)
    else:
        logging.error("Unrecognized tool '" + str(tool) + "' for multiple sequence alignment.\n")
        sys.exit(3)
    return concatenated_msa_files


def create_ref_phy_files(refpkgs: dict, output_dir: str, query_fasta_files: list, ref_aln_dimensions: dict) -> None:
    """
    Creates a phy file for every reference marker that was matched by a query sequence

    :param refpkgs: A dictionary of ReferencePackage instances indexed by their prefix values
    :param output_dir: Path to a directory for writing the Phylip files
    :param query_fasta_files: A list containing paths to FASTA files that are to be converted to Phylip format
    :param ref_aln_dimensions: A dictionary of ref_pkg.prefix keys mapping to a tuple of the nrow, ncol for a MSA
    :return: None
    """

    # Convert the reference sequence alignments to .phy files for every marker identified
    for query_fasta in query_fasta_files:
        marker = re.match("(.*)_hmm_purified.*", os.path.basename(query_fasta)).group(1)
        ref_pkg = refpkgs[marker]  # type: refpkg.ReferencePackage

        ref_alignment_phy = output_dir + marker + ".phy"
        if os.path.isfile(ref_alignment_phy):
            continue

        num_ref_seqs, ref_align_len = ref_aln_dimensions[ref_pkg.prefix]
        aligned_fasta_dict = fasta.read_fasta_to_dict(ref_pkg.f__msa)
        phy_dict = utilities.reformat_fasta_to_phy(aligned_fasta_dict)

        utilities.write_phy_file(ref_alignment_phy, phy_dict, (num_ref_seqs, ref_align_len))
    return


def prepare_and_run_hmmalign(execs: dict, single_query_fasta_files: list, refpkg_dict: dict,
                             output_dir="", n_proc=2) -> dict:
    """
    Runs `hmmalign` to add the query sequences into the reference FASTA multiple alignments

    :param execs: Dictionary of executable file paths indexed by the software names
    :param single_query_fasta_files: List of unaligned query sequences in FASTA format
    :param refpkg_dict: A dictionary of ReferencePackage instances indexed by their respective prefix attributes
    :param output_dir: Where to write the multiple alignment files containing reference and query sequences
    :param n_proc: The number of alignment jobs to run in parallel
    :return: Dictionary of multiple sequence alignment (FASTA) files generated by hmmalign indexed by denominator
    """

    hmmalign_singlehit_files = dict()
    mfa_out_dict = dict()
    logging.info("Running hmmalign... ")

    start_time = time.time()
    task_list = list()

    # Run hmmalign on each fasta file
    for query_fa_in in sorted(single_query_fasta_files):
        file_name_info = re.match(r"(.*)_hmm_purified.*\.(f.*)$", os.path.basename(query_fa_in))
        if file_name_info:
            refpkg_name, extension = file_name_info.groups()
        else:
            logging.error("Unable to parse information from file name:" + "\n" + str(query_fa_in) + "\n")
            sys.exit(3)

        query_mfa_out = os.path.join(output_dir,
                                     re.sub('.' + re.escape(extension) + r"$", ".sto", os.path.basename(query_fa_in)))

        ref_pkg = refpkg_dict[refpkg_name]  # type: refpkg.ReferencePackage
        if ref_pkg.prefix not in hmmalign_singlehit_files:
            hmmalign_singlehit_files[ref_pkg.prefix] = []
        try:
            mfa_out_dict[ref_pkg.prefix].append(query_mfa_out)
        except KeyError:
            mfa_out_dict[ref_pkg.prefix] = [query_mfa_out]

        # Get the paths to either the HMM or CM profile files
        if ref_pkg.kind == "phylogenetic_rRNA":
            task_list.append(wrapper.hmmalign_command(execs["cmalign"], ref_pkg.f__msa, ref_pkg.f__profile,
                                                      query_fa_in, query_mfa_out))
        else:
            task_list.append(wrapper.hmmalign_command(execs["hmmalign"], ref_pkg.f__msa, ref_pkg.f__profile,
                                                      query_fa_in, query_mfa_out))

    if len(task_list) > 0:
        cl_farmer = wrapper.CommandLineFarmer("cmalign/hmmalign --mapali", n_proc)
        cl_farmer.add_tasks_to_queue(task_list)

        cl_farmer.task_queue.close()
        cl_farmer.task_queue.join()

    logging.info("done.\n")

    for prefix in mfa_out_dict:
        for query_mfa_out in mfa_out_dict[prefix]:
            mfa_file = re.sub(r"\.sto$", ".mfa", query_mfa_out)
            seq_dict = file_parsers.read_stockholm_to_dict(query_mfa_out)
            fasta.write_new_fasta(seq_dict, mfa_file)
            hmmalign_singlehit_files[prefix].append(mfa_file)

    end_time = time.time()
    hours, remainder = divmod(end_time - start_time, 3600)
    minutes, seconds = divmod(remainder, 60)
    logging.debug("\thmmalign time required: " +
                  ':'.join([str(hours), str(minutes), str(round(seconds, 2))]) + "\n")

    return hmmalign_singlehit_files


def gather_split_msa(refpkg_names: list, align_dir: str) -> dict:
    """
    Collects the multiple sequence alignment files (FASTA format) in a directory and orders them into a namedtuple.
    The resulting dictionary is used by launch_evolutionary_placement_queries().

    :param refpkg_names: A list of reference package prefix attributes to look for in the align_dir
    :param align_dir: Dictionary path containing the
    :return: A dictionary containing paired reference and query multiple alignments:
        d = {ref_pkg.prefix: MSA(ref.mfa, query.mfa)}
    """
    split_msa_map = {}
    MSAs = namedtuple("MSAs", "ref query")
    for refpkg_name in refpkg_names:
        if not glob.glob(align_dir + refpkg_name + "*"):
            continue
        split_msa_map[refpkg_name] = []
        for ref_msa, query_msa in dict(zip(glob.glob(align_dir + re.escape(refpkg_name) + "*_references.mfa"),
                                           glob.glob(align_dir + re.escape(refpkg_name) + "*_queries.mfa"))).items():
            split_msa = MSAs(ref_msa, query_msa)
            split_msa_map[refpkg_name].append(split_msa)
    return split_msa_map


def check_for_removed_sequences(trimmed_msa_files: dict, msa_files: dict, refpkg_dict: dict, min_len=10):
    """
    Reads the multiple alignment files (either Phylip or FASTA formatted) and looks for both reference and query
    sequences that have been removed. Multiple alignment files are removed from `mfa_files` if:
        1. all query sequences were removed; a DEBUG message is issued
        2. at least one reference sequence was removed
    This quality-control function is necessary for placing short query sequences onto reference trees.

    :param trimmed_msa_files:
    :param msa_files: A dictionary containing the untrimmed MSA files indexed by reference package code (denominator)
    :param refpkg_dict: A dictionary of ReferencePackage objects indexed by their ref_pkg names
    :param min_len: The minimum allowable sequence length after trimming (not including gap characters)
    :return: A dictionary of denominators, with multiple alignment dictionaries as values. Example:
        {M0702: { "McrB_hmm_purified.phy-BMGE.fasta": {'1': seq1, '2': seq2}}}
    """
    qc_ma_dict = dict()
    num_successful_alignments = 0
    discarded_seqs_string = ""
    trimmed_away_seqs = dict()
    untrimmed_msa_failed = []
    logging.debug("Validating trimmed multiple sequence alignment files... ")

    for refpkg_name in sorted(trimmed_msa_files.keys()):
        ref_pkg = refpkg_dict[refpkg_name]  # type: refpkg.ReferencePackage
        trimmed_away_seqs[ref_pkg.prefix] = 0
        # Create a set of the reference sequence names
        ref_headers = fasta.get_headers(ref_pkg.f__msa)
        unique_refs = set([re.sub('_' + re.escape(ref_pkg.prefix), '', x)[1:] for x in ref_headers])
        msa_passed, msa_failed, summary_str = file_parsers.validate_alignment_trimming(trimmed_msa_files[ref_pkg.prefix],
                                                                                       unique_refs, True, min_len)

        # Report the number of sequences that are removed by BMGE
        for trimmed_msa_file in trimmed_msa_files[ref_pkg.prefix]:
            try:
                prefix = re.search('(' + re.escape(ref_pkg.prefix) + r"_.*_group\d+)-(BMGE|trimAl).fasta$",
                                   os.path.basename(trimmed_msa_file)).group(1)
            except TypeError:
                logging.error("Unexpected file name format for a trimmed MSA.\n")
                sys.exit(3)
            # Find the untrimmed query sequence MSA file - the trimmed MSA file's 'pair'
            pair = ""
            for msa_file in msa_files[ref_pkg.prefix]:
                if re.search(re.escape(prefix) + r'\.', msa_file):
                    pair = msa_file
                    break
            if pair:
                if trimmed_msa_file in msa_failed:
                    untrimmed_msa_failed.append(pair)
                trimmed_away_seqs[ref_pkg.prefix] += len(set(fasta.get_headers(pair)).difference(set(fasta.get_headers(trimmed_msa_file))))
            else:
                logging.error("Unable to map trimmed MSA file '" + trimmed_msa_file + "' to its original MSA.\n")
                sys.exit(5)

        if len(msa_failed) > 0:
            if len(untrimmed_msa_failed) != len(msa_failed):
                logging.error("Not all of the failed ({}/{}),"
                              " trimmed MSA files were mapped to their original MSAs."
                              "\n".format(len(msa_failed), len(trimmed_msa_files[ref_pkg.prefix])))
                sys.exit(3)
            untrimmed_msa_passed, _, _ = file_parsers.validate_alignment_trimming(untrimmed_msa_failed, unique_refs,
                                                                                  True, min_len)
            msa_passed.update(untrimmed_msa_passed)
        num_successful_alignments += len(msa_passed)
        qc_ma_dict[ref_pkg.prefix] = msa_passed
        discarded_seqs_string += summary_str
        untrimmed_msa_failed.clear()

    logging.debug("done.\n")
    logging.debug("\tSequences removed during trimming:\n\t\t" +
                  '\n\t\t'.join([k + ": " + str(trimmed_away_seqs[k]) for k in trimmed_away_seqs.keys()]) + "\n")

    logging.debug("\tSequences <" + str(min_len) + " characters removed after trimming:" +
                  discarded_seqs_string + "\n")

    if num_successful_alignments == 0:
        logging.error("No quality alignment files to analyze after trimming. Exiting now.\n")
        sys.exit(0)  # Should be 3, but this allows Clade_exclusion_analyzer to continue after exit

    return qc_ma_dict


def evaluate_trimming_performance(qc_ma_dict, alignment_length_dict, concatenated_msa_files, tool):
    """

    :param qc_ma_dict: A dictionary mapping denominators to files to multiple alignment dictionaries
    :param alignment_length_dict:
    :param concatenated_msa_files: Dictionary with markers indexing original (untrimmed) multiple alignment files
    :param tool: The name of the tool that was appended to the original, untrimmed or unmasked alignment files
    :return: None
    """
    trimmed_length_dict = dict()
    for denominator in sorted(qc_ma_dict.keys()):
        if len(concatenated_msa_files[denominator]) >= 1:
            of_ext = concatenated_msa_files[denominator][0].split('.')[-1]
        else:
            continue
        if denominator not in trimmed_length_dict:
            trimmed_length_dict[denominator] = list()
        for multi_align_file in qc_ma_dict[denominator]:
            file_type = multi_align_file.split('.')[-1]
            multi_align = qc_ma_dict[denominator][multi_align_file]
            num_seqs, trimmed_seq_length = fasta.multiple_alignment_dimensions(multi_align_file, multi_align)

            original_multi_align = re.sub('-' + tool + '.' + file_type, '.' + of_ext, multi_align_file)
            raw_align_len = alignment_length_dict[original_multi_align]
            diff = raw_align_len - trimmed_seq_length
            if diff < 0:
                logging.warning("MSA length increased after {} processing for {}\n".format(tool, multi_align_file))
            else:
                trimmed_length_dict[denominator].append(diff)

    trimming_performance_string = "\tAverage columns removed:\n"
    for denominator in trimmed_length_dict:
        trimming_performance_string += "\t\t" + denominator + "\t"
        n_trimmed_files = len(trimmed_length_dict[denominator])
        if n_trimmed_files > 0:
            trimming_performance_string += str(round(sum(trimmed_length_dict[denominator])/n_trimmed_files, 1)) + "\n"
        else:
            trimming_performance_string += str(0.0) + "\n"

    logging.debug(trimming_performance_string + "\n")
    return


def delete_files(clean_up: bool, root_dir: str, section: int) -> None:
    files_to_be_deleted = []
    if clean_up:
        if section == 1:  # search
            files_to_be_deleted += glob.glob(root_dir + '*_search_to_ORFs_domtbl.txt')
        if section == 2:  # search
            files_to_be_deleted += glob.glob(root_dir + '*_hmm_purified*.faa')
            files_to_be_deleted += glob.glob(root_dir + '*_hmm_purified*.fna')
        if section == 3:  # align
            files_to_be_deleted += glob.glob(root_dir + '*.mfa')
            files_to_be_deleted += glob.glob(root_dir + '*.sto')
            files_to_be_deleted += glob.glob(root_dir + '*_hmm_purified*.fasta')
            files_to_be_deleted += glob.glob(root_dir + '*.phy')
        if section == 4:  # place
            files_to_be_deleted += glob.glob(root_dir + '*_EPA.txt')
            files_to_be_deleted += glob.glob(root_dir + '*EPA_info.txt')
        if section == 5:  # intermediates
            files_to_be_deleted += glob.glob(root_dir + "**/*_formatted.fasta", recursive=True)
            files_to_be_deleted += glob.glob(root_dir + "**/*.sam", recursive=True)
            # Need this for annotate_extra_treesapp.py
            # files_to_be_deleted += glob.glob(root_dir + '*.jplace')

    for useless_file in files_to_be_deleted:
        if os.path.exists(useless_file):
            os.remove(useless_file)
    return


def filter_placements(tree_saps: dict, refpkg_dict: dict, svc: bool, min_lwr: float, max_pendant: float) -> None:
    """
    Determines the total distance of each placement from its branch point on the tree
    and removes the placement if the distance is deemed too great

    :param tree_saps: A dictionary containing PQuery objects
    :param refpkg_dict: A dictionary of ReferencePackage instances indexed by their prefix values
    :param svc: A boolean indicating whether placements should be filtered using ReferencePackage.svc
    :param min_lwr: Likelihood-weight-ratio (LWR) threshold for filtering pqueries
    :param max_pendant: The maximum pendant length distance threshold.
    PQueries with pendant length > max_pendant will be unclassified.
    :return: None
    """
    # The following list must match that of training_utils.vectorize_placement_data_by_rank()
    features = ["evalue", "hmm_cov", "leaves", "lwr", "distal", "pendant", "avg_tip_dist"]

    logging.info("Filtering low-quality placements... ")
    unclassified_seqs = dict()  # A dictionary tracking the seqs unclassified for each marker

    for refpkg_name, pqueries in tree_saps.items():  # type: (str, list)
        ref_pkg = refpkg_dict[refpkg_name]  # type: refpkg.ReferencePackage
        unclassified_seqs[ref_pkg.prefix] = dict()
        unclassified_seqs[ref_pkg.prefix]["low_lwr"] = list()
        unclassified_seqs[ref_pkg.prefix]["big_pendant"] = list()
        unclassified_seqs[ref_pkg.prefix]["svm"] = list()
        svc_attempt = False

        for tree_sap in sorted(pqueries, key=lambda x: x.seq_name):  # type: phylo_seq.PQuery
            tree_sap.filter_min_weight_threshold(min_lwr)
            if not tree_sap.classified:
                unclassified_seqs[ref_pkg.prefix]["low_lwr"].append(tree_sap)
                continue
            elif not tree_sap.placements:
                continue

            tree_sap.avg_evo_dist = tree_sap.consensus_placement.total_distance()
            tree_sap.string_distances()

            if tree_sap.consensus_placement.pendant_length > max_pendant:
                unclassified_seqs[ref_pkg.prefix]["big_pendant"].append(tree_sap)
                tree_sap.classified = False

            if svc:
                if ref_pkg.svc is None or len(features) != ref_pkg.svc.n_features_in_:
                    svc_attempt = True
                    call = 1
                else:
                    classification_tup = training_utils.pquery_to_vector(tree_sap, ref_pkg)
                    call = ref_pkg.svc.predict(np_array([getattr(classification_tup, feat) for feat in features]).reshape(1, -1))
                # Discard this placement as a false positive if classifier calls this a 0
                if call == 0:
                    unclassified_seqs[tree_sap.ref_name]["svm"].append(tree_sap)
                    tree_sap.classified = False

        if svc_attempt:
            logging.warning("SVM classifier unavailable for reference package '{}'\n".format(ref_pkg.prefix))

    logging.info("done.\n")

    declass_summary = ""
    for marker in unclassified_seqs:
        # unclassified_counts[marker] will always be >= distant_seqs[marker]
        for declass in unclassified_seqs[marker]:
            declass_summary += marker + '\t' + declass + '\t' + str(len(unclassified_seqs[marker][declass])) + "\n"

    logging.debug(declass_summary)

    return


def select_query_placements(pquery_dict: dict, refpkg_dict: dict, mode="max_lwr"):
    """


    :return: Dictionary of PQuery instances indexed by denominator (ref_pkg code e.g. M0701)
    """

    logging.info('Selecting the optimal query placements... ')

    function_start_time = time.time()
    classified_seqs = 0

    for refpkg_code in pquery_dict:  # type: str
        ref_pkg = refpkg_dict[refpkg_code]  # type: refpkg.ReferencePackage
        taxa_tree = ref_pkg.taxonomically_label_tree()
        for pquery in pquery_dict[refpkg_code]:  # type: phylo_seq.PQuery
            if mode == "max_lwr":
                pquery.process_max_weight_placement(taxa_tree)
            elif mode == "aelw":
                pquery.calculate_consensus_placement(taxa_tree)
            else:
                logging.error("Unknown PQuery consensus algorithm provided: '{}'.\n".format(mode))
                raise ValueError

            classified_seqs += 1
            pquery.placements = [pquery.consensus_placement]

            # I have decided to not remove the original JPlace files since some may find these useful
            # os.remove(filename)

    logging.info("done.\n")

    function_end_time = time.time()
    hours, remainder = divmod(function_end_time - function_start_time, 3600)
    minutes, seconds = divmod(remainder, 60)
    logging.debug("\tPQuery parsing time required: " +
                  ':'.join([str(hours), str(minutes), str(round(seconds, 2))]) + "\n")
    logging.debug("\t" + str(classified_seqs) + " sequences placed into trees by EPA-NG.\n")

    return pquery_dict


def parse_raxml_output(epa_output_dir: str, refpkg_dict: dict, pqueries=None):
    """
    For every JPlace file found in the directory **epa_output_dir**, all placed query sequences in the JPlace
    are demultiplexed, each becoming a PQuery instance.
    The JPlace data are validated by ensuring the distal placement lengths reported by EPA are less than or equal to
    the corresponding edge length in the JPlace tree.

    :param epa_output_dir: Directory where EPA wrote the JPlace files
    :param refpkg_dict: A dictionary of ReferencePackage instances indexed by their prefix values
    :param pqueries: A list of instantiated PQuery instances
    :return:
        1. Dictionary of PQuery instances indexed by denominator (ref_pkg code e.g. M0701)
        2. Dictionary of an JPlace instance (values) mapped to a ref_pkg prefix.
         These instances store all PQuery instances from all JPlace files
    """

    logging.info('Parsing the EPA-NG outputs... ')

    function_start_time = time.time()

    jplace_files = glob.glob(epa_output_dir + '*.jplace')
    itol_data = dict()  # contains all pqueries, indexed by marker ref_name (e.g. McrA, nosZ, 16srRNA)
    tree_saps = dict()  # contains individual pquery information for each mapped protein (N==1), indexed by denominator
    # Use the jplace files to guide which markers iTOL outputs should be created for
    if pqueries:
        pquery_map = {pq.place_name: pq for pq in pqueries}
    else:
        pquery_map = None

    for refpkg_name, jplace_list in jplace_utils.organize_jplace_files(jplace_files).items():
        ref_pkg = refpkg_dict[refpkg_name]
        if ref_pkg.prefix not in tree_saps:
            tree_saps[ref_pkg.prefix] = list()
        for filename in jplace_list:
            # Load the JSON placement (jplace) file containing >= 1 pquery into JPlace object
            jplace_data = jplace_utils.jplace_parser(filename)
            edge_dist_index = index_tree_edges(jplace_data.tree)
            internal_node_leaf_map = map_internal_nodes_leaves(jplace_data.tree)
            # Demultiplex all pqueries in jplace_data into individual PQuery objects
            jplace_data.pqueries = jplace_utils.demultiplex_pqueries(jplace_data, pquery_map)
            jplace_utils.calc_pquery_mean_tip_distances(jplace_data, internal_node_leaf_map)
            for pquery in jplace_data.pqueries:  # type: phylo_seq.PQuery
                # Flesh out the internal-leaf node map
                pquery.ref_name = ref_pkg.prefix
                if not pquery.seq_name:
                    seq_info = re.match(r"(.*)\|" + re.escape(pquery.ref_name) + r"\|(\d+)_(\d+)$", pquery.place_name)
                    pquery.seq_name, pquery.start, pquery.end = seq_info.groups()
                pquery.seq_len = int(pquery.end) - int(pquery.start)
                pquery.node_map = internal_node_leaf_map
                pquery.check_jplace_edge_lengths(edge_dist_index)
                tree_saps[ref_pkg.prefix].append(pquery)

            if ref_pkg.prefix not in itol_data:
                itol_data[ref_pkg.prefix] = jplace_data
                itol_data[ref_pkg.prefix].ref_name = ref_pkg.prefix
            else:
                # If a JPlace file for that tree has already been parsed, just append the placements
                itol_data[ref_pkg.prefix].pqueries = itol_data[ref_pkg.prefix].pqueries + jplace_data.pqueries

            # I have decided to not remove the original JPlace files since some may find these useful
            # os.remove(filename)

    logging.info("done.\n")

    function_end_time = time.time()
    hours, remainder = divmod(function_end_time - function_start_time, 3600)
    minutes, seconds = divmod(remainder, 60)
    logging.debug("\tJPlace parsing time required: " +
                  ':'.join([str(hours), str(minutes), str(round(seconds, 2))]) + "\n")
    logging.debug("\t" + str(len(jplace_files)) + " JPlace files.\n")

    return tree_saps, itol_data


def write_classified_sequences(tree_saps: dict, formatted_fasta_dict: dict, fasta_file: str) -> None:
    """
    Function to write the nucleotide sequences representing the full-length ORF for each classified sequence
    Sequence names are from JPlace.contig_name values so output format is:

     >contig_name|RefPkg|StartCoord_StopCoord

    :param tree_saps: A dictionary of gene_codes as keys and TreeSap objects as values
    :param formatted_fasta_dict: A dictionary with headers/sequence names as keys and sequences as values
    :param fasta_file: Path to a file to write the sequences to in FASTA format
    :return: None
    """
    output_fasta_dict = dict()
    len_parsing_problem = False
    prefix = ''  # For adding a '>' if the formatted_fasta_dict sequences have them
    for seq_name in formatted_fasta_dict:
        if seq_name[0] == '>':
            prefix = '>'
        break

    for denominator in tree_saps:
        for placed_sequence in tree_saps[denominator]:  # type: phylo_seq.PQuery
            if placed_sequence.classified:
                output_fasta_dict[placed_sequence.place_name] = ""
                try:
                    output_fasta_dict[placed_sequence.place_name] = formatted_fasta_dict[prefix +
                                                                                         placed_sequence.place_name]
                except KeyError:
                    seq_name = re.sub(r"\|{0}\|\d+_\d+.*".format(placed_sequence.ref_name),
                                      '',
                                      placed_sequence.place_name)
                    try:
                        output_fasta_dict[placed_sequence.place_name] = formatted_fasta_dict[prefix + seq_name]
                    except KeyError:
                        logging.error("Unable to find '" + prefix + placed_sequence.place_name +
                                      "' in predicted ORFs file!\nExample headers in the predicted ORFs file:\n\t" +
                                      '\n\t'.join(list(formatted_fasta_dict.keys())[:6]) + "\n")
                        sys.exit(3)

                if not placed_sequence.seq_len:
                    placed_sequence.seq_len = len(output_fasta_dict[placed_sequence.place_name])
                    len_parsing_problem = True

    if output_fasta_dict:
        fasta.write_new_fasta(output_fasta_dict, fasta_file)

    if len_parsing_problem:
        logging.warning("Problem parsing homologous subsequence lengths from headers of classified sequences.\n")

    return


def determine_confident_lineage(tree_saps: dict, refpkg_dict: dict, mode="max_lwr") -> None:
    """
    Determines the best taxonomic lineage for classified sequences based on their

1. placement in the phylogeny
2. the lowest common ancestor of all children to the placement edge
3. the optimal rank recommended by the linear model

    :param tree_saps: A dictionary containing PQuery objects
    :param refpkg_dict: A dictionary of ReferencePackage instances indexed by their prefix values
    :param mode: The phylogenetic placement and taxonomic assignment setting, either 'max_lwr' or 'aelw'
    :return: None
    """
    leaf_taxa_map = dict()
    for refpkg_name in tree_saps:
        # All the leaves for that tree [number, translation, lineage]
        ref_pkg = refpkg_dict[refpkg_name]  # type: refpkg.ReferencePackage
        for leaf in ref_pkg.generate_tree_leaf_references_from_refpkg():
            leaf_taxa_map[leaf.number] = leaf.lineage

        for pquery in tree_saps[refpkg_name]:  # type: phylo_seq.PQuery
            if not pquery.classified:
                continue

            lineage_list = pquery.children_lineage(leaf_taxa_map)
            # algorithm options are "MEGAN", "LCAp", and "LCA*" (default)
            # pquery.lct = lowest_common_taxonomy(lineage_list, lca, taxonomic_counts, "LCA*")
            pquery.wtd, status = ts_lca.weighted_taxonomic_distance(lineage_list, pquery.lct)
            if status > 0:
                pquery.summarize()

            recommended_rank = phylo_dist.rank_recommender(pquery.avg_evo_dist, refpkg_dict[refpkg_name].pfit)
            if pquery.lct.split(ref_pkg.taxa_trie.lin_sep)[0] != "r__Root":
                pquery.lct = "r__Root; " + pquery.lct
                recommended_rank += 1

            if mode in ["max_lwr"]:
                # Based on the calculated distance from the leaves, what rank is most appropriate?
                pquery.recommended_lineage = pquery.lowest_confident_taxonomy(recommended_rank)
            else:
                pquery.recommended_lineage = pquery.lct
        leaf_taxa_map.clear()
    return


def produce_itol_inputs(pqueries: dict, refpkg_dict: dict, jplaces: dict,
                        itol_base_dir: str, treesapp_data_dir: str) -> None:
    """
    Function to create outputs for the interactive tree of life (iTOL) webservice.
    There is a directory for each of the marker genes detected to allow the user to "drag-and-drop" all files easily

    :param pqueries: Dictionary of PQuery instances indexed by denominator (ref_pkg code e.g. M0701)
    :param refpkg_dict: A dictionary of ReferencePackage instances indexed by their prefix values
    :param jplaces: Dictionary of an JPlace instance (values) mapped to marker name
    :param itol_base_dir: Output directory to write the iTOL_output files, with the outputs stored within a directory
     named after their respective reference package prefix (e.g. McrA) 
    :param treesapp_data_dir: Path to the directory containing reference packages and a iTOL_data directory
    :return: None
    """
    logging.info("Generating inputs for iTOL... ")

    if not os.path.exists(itol_base_dir):
        os.mkdir(itol_base_dir)  # drwxr-xr-x
    # Now that all the JPlace files have been loaded, generate the abundance stats for each marker

    strip_missing = []
    style_missing = []
    for refpkg_name in pqueries:
        if len(pqueries[refpkg_name]) == 0:
            # No sequences that were mapped met the minimum likelihood weight ration threshold. Skipping!
            continue
        ref_pkg = refpkg_dict[refpkg_name]  # type: refpkg.ReferencePackage
        if not os.path.exists(itol_base_dir + ref_pkg.prefix):
            os.mkdir(itol_base_dir + ref_pkg.prefix)
        jplace_data = jplaces[ref_pkg.prefix]
        refpkg_pqueries = pqueries[refpkg_name]

        if os.path.isfile(ref_pkg.f__boot_tree):
            jplace_utils.add_bipartitions(jplace_data, ref_pkg.f__boot_tree)

        # Make a master jplace file from the set of placements in all jplace files for each marker
        master_jplace = os.path.join(itol_base_dir, ref_pkg.prefix, ref_pkg.prefix + "_complete_profile.jplace")
        # TODO: validate no distal lengths exceed their corresponding edge lengths

        jplace_data.write_jplace(master_jplace)
        jplaces[ref_pkg.prefix].clear_object()
        jplace_data.clear_object()
        # Create a labels file from the tax_ids_marker.txt
        ref_pkg.create_itol_labels(itol_base_dir)

        annotation_style_files = glob.glob(os.sep.join([treesapp_data_dir, "iTOL_data", ref_pkg.prefix + "*"]))
        # Copy the respective colours and styles files for each marker found to the itol_output directories
        colours_styles = os.sep.join([treesapp_data_dir, "iTOL_data", ref_pkg.prefix + "_colours_style.txt"])
        colour_strip = os.sep.join([treesapp_data_dir, "iTOL_data", ref_pkg.prefix + "_colour_strip.txt"])
        if colours_styles not in annotation_style_files:
            style_missing.append(ref_pkg.prefix)
        if colour_strip not in annotation_style_files:
            strip_missing.append(ref_pkg.prefix)

        for annotation_file in annotation_style_files:
            shutil.copy(annotation_file, itol_base_dir + ref_pkg.prefix)
        itol_bar_file = os.path.join(itol_base_dir, ref_pkg.prefix, ref_pkg.prefix + "_abundance_simplebar.txt")
        abundance.generate_simplebar(ref_pkg.prefix, refpkg_pqueries, itol_bar_file)

    logging.info("done.\n")
    if style_missing:
        logging.debug("A colours_style.txt file does not yet exist for markers:\n\t" +
                      "\n\t".join(style_missing) + "\n")
    if strip_missing:
        logging.debug("A colours_strip.txt file does not yet exist for markers:\n\t" +
                      "\n\t".join(strip_missing) + "\n")

    return


def alert_for_refpkg_feature_annotations(pqueries: dict, refpkg_dict: dict) -> None:
    feature_positive = []
    for refpkg_name in pqueries:
        if len(pqueries[refpkg_name]) == 0:
            continue
        ref_pkg = refpkg_dict[refpkg_name]  # type: refpkg.ReferencePackage
        if len(ref_pkg.feature_annotations) > 0:
            feature_positive.append(ref_pkg.prefix)

    if len(feature_positive) > 0:
        logging.info("Alert: {} reference package(s) have feature annotations. "
                     "Consider running treesapp layer.\n".format(len(feature_positive)))
        logging.debug("Reference packages with clade_annotations attribute filled:\n\t{}\n"
                      "".format("\n\t".join(feature_positive)))

    return


def assign(sys_args):
    # STAGE 1: Prompt the user and prepare files and lists for the pipeline
    parser = treesapp_args.TreeSAPPArgumentParser(description='Classify sequences through evolutionary placement.')
    treesapp_args.add_classify_arguments(parser)
    args = parser.parse_args(sys_args)

    ts_assign = Assigner()
    ts_assign.furnish_with_arguments(args)
    ts_assign.check_previous_output(args.overwrite)

    log_file_name = args.output + os.sep + "TreeSAPP_classify_log.txt"
    classy.prep_logging(log_file_name, args.verbose)
    logging.info("\n##\t\t\t\tAssigning sequences with TreeSAPP\t\t\t\t##\n\n")

    treesapp_args.check_parser_arguments(args, sys_args)
    ts_assign.check_classify_arguments(args)
    hmm_parsing_thresholds = ts_assign.define_hmm_domtbl_thresholds(args)
    ts_assign.decide_stage(args)
    n_proc = args.num_threads

    refpkg_dict = refpkg.gather_ref_packages(ts_assign.refpkg_dir, ts_assign.target_refpkgs)
    prep_reference_packages_for_assign(refpkg_dict, ts_assign.var_output_dir)
    ref_alignment_dimensions = get_alignment_dims(refpkg_dict)
    ts_assign.load_refpkg_classifiers(refpkg_dict, kernel=args.kernel, threads=n_proc)

    ##
    # STAGE 2: Predict open reading frames (ORFs) if the input is an assembly, read, format and write the FASTA
    ##
    if ts_assign.stage_status("orf-call"):
        ts_assign.predict_orfs(args.composition, n_proc)

    query_seqs = fasta.FASTA(ts_assign.query_sequences)
    # Read the query sequences provided and (by default) write a new FASTA file with formatted headers
    if ts_assign.stage_status("clean"):
        logging.info("Reading and formatting {}... ".format(ts_assign.query_sequences))
        query_seqs.header_registry = fasta.format_fasta(fasta_input=ts_assign.query_sequences, molecule="prot",
                                                        output_fasta=ts_assign.formatted_input,
                                                        full_name=ts_assign.fasta_full_name)
        logging.info("done.\n")
    else:
        query_seqs.file = ts_assign.formatted_input
        query_seqs.header_registry = fasta.register_headers(fasta.get_headers(query_seqs.file))
    logging.info("\tTreeSAPP will analyze the " + str(len(query_seqs.header_registry)) + " sequences found in input.\n")
    ts_assign.increment_stage_dir()

    ##
    # STAGE 3: Run hmmsearch on the query sequences to search for marker homologs
    ##
    if ts_assign.stage_status("search"):
        hmm_domtbl_files = wrapper.hmmsearch_orfs(ts_assign.executables["hmmsearch"],
                                                  refpkg_dict, ts_assign.formatted_input,
                                                  ts_assign.stage_output_dir, n_proc, hmm_parsing_thresholds.max_e)
    else:
        hmm_domtbl_files = ts_assign.fetch_hmmsearch_outputs(set(refpkg_dict.keys()))

    # Load alignment information
    hmm_matches = file_parsers.parse_domain_tables(hmm_parsing_thresholds, hmm_domtbl_files)
    load_homologs(hmm_matches, ts_assign.formatted_input, query_seqs)
    pqueries = load_pqueries(hmm_matches, query_seqs)
    query_seqs.change_dict_keys("num")
    extracted_seq_dict, numeric_contig_index = bin_hmm_matches(hmm_matches, query_seqs.fasta_dict)
    numeric_contig_index = replace_contig_names(numeric_contig_index, query_seqs)
    homolog_seq_files = write_grouped_fastas(extracted_seq_dict, numeric_contig_index,
                                             refpkg_dict, ts_assign.stage_output_dir)
    # TODO: Replace this merge_fasta_dicts_by_index with FASTA - only necessary for writing the classified sequences
    extracted_seq_dict = fasta.merge_fasta_dicts_by_index(extracted_seq_dict, numeric_contig_index)
    delete_files(args.delete, ts_assign.stage_output_dir, 1)
    ts_assign.increment_stage_dir()

    ##
    # STAGE 4: Run hmmalign, and optionally BMGE, to produce the MSAs for phylogenetic placement
    ##
    MSAs = namedtuple("MSAs", "ref query")
    if ts_assign.stage_status("align"):
        combined_msa_files = {}
        split_msa_files = {}
        create_ref_phy_files(refpkg_dict, ts_assign.stage_output_dir,
                             homolog_seq_files, ref_alignment_dimensions)
        concatenated_msa_files = multiple_alignments(ts_assign.executables, homolog_seq_files,
                                                     refpkg_dict, "hmmalign",
                                                     ts_assign.stage_output_dir, num_proc=n_proc)
        file_type = utilities.find_msa_type(concatenated_msa_files)
        alignment_length_dict = get_sequence_counts(concatenated_msa_files, ref_alignment_dimensions,
                                                    args.verbose, file_type)

        if args.trim_align:
            tool = "BMGE"
            trimmed_mfa_files = wrapper.filter_multiple_alignments(ts_assign.executables, concatenated_msa_files,
                                                                   refpkg_dict, n_proc, tool)
            qc_ma_dict = check_for_removed_sequences(trimmed_mfa_files, concatenated_msa_files,
                                                     refpkg_dict, args.min_seq_length)
            evaluate_trimming_performance(qc_ma_dict, alignment_length_dict, concatenated_msa_files, tool)
            combined_msa_files.update(qc_ma_dict)
        else:
            combined_msa_files.update(concatenated_msa_files)

        # Subset the multiple alignment of reference sequences and queries to just contain query sequences
        for refpkg_name in combined_msa_files:
            split_msa_files[refpkg_name] = []
            for combined_msa in combined_msa_files[refpkg_name]:
                split_msa = MSAs(os.path.dirname(combined_msa) + os.sep +
                                 os.path.basename('.'.join(combined_msa.split('.')[:-1])) + "_references.mfa",
                                 os.path.dirname(combined_msa) + os.sep +
                                 os.path.basename('.'.join(combined_msa.split('.')[:-1])) + "_queries.mfa")
                fasta.split_combined_ref_query_fasta(combined_msa, split_msa.query, split_msa.ref)
                split_msa_files[refpkg_name].append(split_msa)
        combined_msa_files.clear()
        delete_files(args.delete, ts_assign.stage_lookup("search").dir_path, 2)
    else:
        split_msa_files = gather_split_msa(list(refpkg_dict.keys()), ts_assign.stage_output_dir)
    ts_assign.increment_stage_dir()

    ##
    # STAGE 5: Run EPA-ng to compute the ML estimations
    ##
    if ts_assign.stage_status("place"):
        wrapper.launch_evolutionary_placement_queries(ts_assign.executables, split_msa_files, refpkg_dict,
                                                      ts_assign.stage_output_dir, n_proc)
        jplace_utils.sub_indices_for_seq_names_jplace(ts_assign.stage_output_dir, numeric_contig_index, refpkg_dict)
        delete_files(args.delete, ts_assign.stage_lookup("align").dir_path, 3)

    if ts_assign.stage_status("classify"):
        tree_saps, itol_data = parse_raxml_output(ts_assign.stage_output_dir, refpkg_dict, pqueries)
        ts_assign.increment_stage_dir()
        # Set PQuery.consensus_placement attributes
        select_query_placements(tree_saps, refpkg_dict, mode=args.p_sum)
        filter_placements(tree_saps, refpkg_dict, ts_assign.svc_filter, args.min_lwr, args.max_pd)
        determine_confident_lineage(tree_saps, refpkg_dict)

        ts_assign.write_classified_orfs(tree_saps, extracted_seq_dict)
        abundance_dict = dict()
        for refpkg_code in tree_saps:
            for placed_seq in tree_saps[refpkg_code]:  # type: phylo_seq.PQuery
                abundance_dict[placed_seq.place_name] = 1.0

        # Run the abundance subcommand on the classified sequences
        if ts_assign.stage_status("abundance"):
            abundance_args = ["--treesapp_output", ts_assign.output_dir,
                              "--reads", ' '.join(args.reads),
                              "--pairing", args.pairing,
                              "--num_procs", str(n_proc),
                              "--report", "nothing",
                              "--metric", args.metric]
            if args.reverse:
                abundance_args += ["--reverse", ' '.join(args.reverse)]
            sample_id, abundance_dict = abundance.abundance(abundance_args).popitem()
        phylo_seq.abundify_tree_saps(tree_saps, abundance_dict)

        file_parsers.write_classification_table(tree_saps, ts_assign.sample_prefix,
                                                output_file=ts_assign.classification_table)

        produce_itol_inputs(tree_saps, refpkg_dict, itol_data, ts_assign.itol_out, ts_assign.refpkg_dir)
        delete_files(args.delete, ts_assign.stage_lookup("place").dir_path, 4)

        alert_for_refpkg_feature_annotations(tree_saps, refpkg_dict)

    # Clear out the rest of the intermediates
    delete_files(args.delete, ts_assign.var_output_dir, 5)

    return
