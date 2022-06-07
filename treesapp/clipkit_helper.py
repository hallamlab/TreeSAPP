import sys
import logging
import os.path

from clipkit import clipkit as ck
from clipkit import modes as ck_modes

from treesapp import logger
from treesapp import fasta
from treesapp import refpkg
from treesapp import file_parsers


class ClipKitHelper:
    CLIPKIT_MODES = {"smart-gap"}

    def __init__(self, fasta_in: str, output_dir: str, mode="smart-gap", gap_prop=0.9):
        self.input = fasta_in
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)

        prefix, ext = os.path.splitext(os.path.basename(fasta_in))
        self.mfa_out = os.path.join(output_dir, prefix + ".trim" + ext)
        self.qc_mfa_out = os.path.join(output_dir, prefix + ".trim.qc" + ext)

        self.logger = logging.getLogger(logger.logger_name())
        self.mode = ck_modes.TrimmingMode(mode)
        self.gap_prop = gap_prop

        self.ff_in = "fasta"
        self.ff_out = "fasta"
        self.refpkg_name = ''
        self.min_unaligned_seq_length = 1

        # Attributes used in evaluating trimming performance
        self.success = True
        self.num_msa_seqs = 0
        self.num_msa_cols = 0
        self.num_trim_seqs = 0
        self.num_trim_cols = 0
        self.trim_qc_seqs = []  # These sequences passed the min_unaligned_seq_length filter

        # Specific to MSAs for phylogenetic placement
        self.num_queries_failed_trimming = 0
        self.num_refs_failed_trimming = 0
        self.num_queries_failed_qc = 0
        self.num_refs_failed_qc = 0
        self.num_queries_retained = 0
        self.num_refs_retained = 0
        return

    def run(self):
        ck.execute(input_file=self.input,
                   input_file_format=self.ff_in,
                   output_file=self.mfa_out,
                   output_file_format=self.ff_out,
                   gaps=self.gap_prop,
                   complement=False,
                   mode=self.mode,
                   use_log=False)
        return

    def summarise_trimming(self):
        if self.num_trim_seqs == 0:
            self.logger.warning("No sequences were read from {}.\n".format(self.mfa_out))

        if self.num_trim_cols < self.min_unaligned_seq_length:
            # Throw an error if the final trimmed alignment is shorter than min_seq_length, and therefore empty
            self.logger.warning(
                "Multiple sequence alignment in {} is shorter than minimum sequence length threshold ({}).\n"
                "".format(self.mfa_out, self.min_unaligned_seq_length))
        elif self.num_refs_failed_trimming:
            # Testing whether there were more sequences in the untrimmed alignment than the trimmed one
            self.logger.warning(
                "{} reference sequences in {} were removed during alignment trimming " +
                "suggesting either truncated sequences or the initial reference alignment was terrible.\n"
                "".format(self.num_refs_failed_trimming, self.mfa_out))
        elif self.num_refs_failed_qc:
            self.logger.warning("{} reference sequences in {} were shorter than the minimum character length ({})"
                                " and removed after alignment trimming.\n"
                                "".format(self.num_refs_failed_qc, self.mfa_out, self.min_unaligned_seq_length))

        # Ensure that there is at least 1 query sequence retained after trimming the multiple alignment
        elif self.num_queries_retained == 0:
            self.logger.warning("No query sequences in {} were retained after trimming.\n".format(self.mfa_out))

        if self.success is False:
            self.logger.debug("The untrimmed MSA will be used instead.\n")
        return

    def quantify_refs_and_pqueries(self, unique_ref_headers: set, msa_fasta: fasta.FASTA = None):
        if not unique_ref_headers:
            return

        if not msa_fasta:
            msa_fasta = self.read_trimmed_msa()

        for seq_name in msa_fasta.fasta_dict:
            if seq_name[0] == '-':  # The negative integers indicate this is a query sequence
                if seq_name in self.trim_qc_seqs:
                    self.num_queries_retained += 1
                else:
                    self.num_queries_failed_qc += 1
            elif seq_name in unique_ref_headers:
                if seq_name in self.trim_qc_seqs:
                    self.num_refs_retained += 1
                else:
                    self.num_refs_failed_qc += 1
            else:
                raise RuntimeError("Unsure what to do with sequence '{}'.\n".format(seq_name))
        return

    def read_trimmed_msa(self) -> fasta.FASTA:
        msa_records = fasta.FASTA(file_name=self.mfa_out)
        if self.ff_out == "phylip":
            msa_records.fasta_dict = file_parsers.read_phylip_to_dict(self.mfa_out)
        elif self.ff_out == "fasta":
            msa_records.fasta_dict = fasta.read_fasta_to_dict(self.mfa_out)
        else:
            self.logger.error("Unsupported file format ('{}') of {}.\n".format(self.ff_out, self.mfa_out))
            sys.exit(1)

        msa_records.header_registry = fasta.register_headers(list(msa_records.fasta_dict.keys()),
                                                             drop=True)
        return msa_records
    
    def validate_alignment_trimming(self):
        msa_fasta = self.read_trimmed_msa()
        if self.num_trim_seqs == 0:
            self.success = False

        if self.num_trim_cols < self.min_unaligned_seq_length:
            self.success = False

        if self.num_trim_cols > self.num_msa_cols:
            self.logger.warning("MSA length increased after trimming {}\n".format(self.input))
            self.success = False

        msa_fasta.unalign()
        for seq_name, seq in msa_fasta.fasta_dict.items():
            if len(seq) >= self.min_unaligned_seq_length:
                self.trim_qc_seqs.append(seq_name)

        return

    def compare_original_and_trimmed_multiple_alignments(self, min_len: int, ref_pkg=None):
        """Summarises the number of character positions trimmed and new dimensions between the input and output MSA."""

        self.num_trim_seqs, self.num_trim_cols = fasta.multiple_alignment_dimensions(self.mfa_out)
        self.num_msa_seqs, self.num_msa_cols = fasta.multiple_alignment_dimensions(self.input)

        self.min_unaligned_seq_length = min_len
        self.validate_alignment_trimming()

        if ref_pkg is not None:  # type: refpkg.ReferencePackage
            # Create a set of the reference sequence names
            unique_ref_headers = set(ref_pkg.get_fasta().get_seq_names())
            self.quantify_refs_and_pqueries(unique_ref_headers)

        return

    def get_qc_output(self) -> str:
        if self.success:
            return self.qc_mfa_out
        else:
            return self.input

    def get_qc_trimmed_fasta(self) -> fasta.FASTA:
        if not self.success:
            return

        msa_fasta = self.read_trimmed_msa()
        msa_fasta.keep_only(header_subset=self.trim_qc_seqs)
        return msa_fasta

    def write_qc_trimmed_multiple_alignment(self) -> None:
        msa_fasta = self.get_qc_trimmed_fasta()
        fasta.write_new_fasta(fasta_dict=msa_fasta.fasta_dict,
                              fasta_name=self.qc_mfa_out)
        return
