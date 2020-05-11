#!/usr/bin/env python3

import logging
import os
import argparse

from . import fasta
from treesapp.refpkg import ReferencePackage
from .wrapper import build_hmm_profile, run_mafft
from .phylo_dist import trim_lineages_to_rank
from .utilities import base_file_prefix


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--refpkg",
                        required=True,
                        help="Name of the reference package")
    parser.add_argument("-o", "--output_hmm",
                        required=True,
                        help="HMM file to build")
    parser.add_argument("-d", "--refpkg_dir",
                        required=True,
                        help="The directory containing the reference package files")
    parser.add_argument("-r", "--rank")
    parser.add_argument("--he")
    parser.add_argument("--me")
    parser.add_argument("-n", "--num_threads",
                        required=False, default=2)
    parser.add_argument("-x", "--tmp_dir")

    return parser.parse_args()


def make_dereplicated_hmm(refpkg_name: str, package_path: str, dereplication_rank: str,
                          hmmbuild: str, mafft: str, hmm_file: str, n_threads=2, intermediates_dir=None) -> None:
    """
    Function to create a taxonomically-dereplicated hidden Markov model (HMM) profile. This reduces the bias from
    potentially over-represented clades, increasing the weight of their conserved positions.

    :param refpkg_name: Short-form gene/protein name of the reference package (e.g. McrA, DsrAB)
    :param package_path: Path to directory containing reference package files
    :param dereplication_rank: The taxonomic rank to dereplicate the reference sequences at
    :param hmmbuild: Path to an hmmbuild executable
    :param mafft: Path to a MAFFT executable
    :param hmm_file: Name of the HMM file to be written
    :param n_threads: Number of threads MAFFT should use
    :param intermediates_dir: A path to a directory to write intermediate files
    :return: None
    """

    logging.info("Creating taxonomically-dereplicated HMM... ")

    refpkg = ReferencePackage(refpkg_name)
    try:
        # This is the most frequent case, where make_dereplicated_hmm is called in treesapp create
        refpkg.gather_package_files(package_path, layout="flat")
    except AssertionError:
        refpkg.gather_package_files(package_path)

    if not intermediates_dir:
        intermediates_dir = os.path.dirname(refpkg.msa)
    if intermediates_dir[-1] != os.sep:
        intermediates_dir += os.sep
    derep_fa = intermediates_dir + base_file_prefix(refpkg.msa) + "_derep.fa"
    derep_aln = intermediates_dir + base_file_prefix(refpkg.msa) + "_derep.mfa"
    intermediates = [derep_aln, derep_fa]

    lineage_reps = []
    t = {}

    mfa = fasta.FASTA(refpkg.msa)
    mfa.load_fasta()

    leaf_nodes = refpkg.tax_ids_file_to_leaves()
    # Trim the taxonomic lineages to the dereplication level
    leaf_taxa_map = {leaf.number + "_" + refpkg_name: leaf.lineage for leaf in leaf_nodes}
    trimmed_lineages = trim_lineages_to_rank(leaf_taxa_map, dereplication_rank)
    # Add back the sequences with truncated lineages
    for leaf_name in leaf_taxa_map:
        if leaf_name not in trimmed_lineages:
            trimmed_lineages[leaf_name] = leaf_taxa_map[leaf_name]

    # Find the longest sequence that each lineage
    for seq_name in mfa.get_seq_names():
        taxon = trimmed_lineages[seq_name]
        if taxon in t:
            curr_rep = t[taxon]
            if len(mfa.fasta_dict[curr_rep]) > len(mfa.fasta_dict[seq_name]):
                continue
        t[taxon] = seq_name
    for taxon in t:
        lineage_reps.append(t[taxon])

    logging.debug("%i %s-dereplicated sequences retained for building HMM profile.\n" %
                  (len(lineage_reps), dereplication_rank))

    # Remove all sequences from the FASTA instance that are not representatives
    mfa.keep_only(lineage_reps)
    mfa.unalign()

    # Write the dereplicated FASTA file
    fasta.write_new_fasta(fasta_dict=mfa.fasta_dict, fasta_name=derep_fa)

    # Re-align the sequences
    run_mafft(mafft_exe=mafft, fasta_in=derep_fa, fasta_out=derep_aln, num_threads=n_threads)

    # Build the new HMM profile
    build_hmm_profile(hmmbuild_exe=hmmbuild, msa_in=derep_aln, output_hmm=hmm_file, name=refpkg_name)

    # Clean up intermediates
    for f_path in intermediates:
        if os.path.isfile(f_path):
            os.remove(f_path)

    logging.info("done.\n")

    return


if __name__ == "__main__":
    args = get_arguments()
    make_dereplicated_hmm(args.refpkg, args.refpkg_dir, args.rank, args.he, args.me, args.output_hmm,
                          args.num_threads, args.tmp_dir)
