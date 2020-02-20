#!/usr/bin/env python3

import logging
import os
import argparse
import re
import sys

from treesapp import fasta
from treesapp.file_parsers import tax_ids_file_to_leaves
from treesapp.wrapper import build_hmm_profile, run_mafft
from treesapp.phylo_dist import trim_lineages_to_rank


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_fa",
                        required=True,
                        help="Multiple alignment FASTA file.")
    parser.add_argument("-o", "--output_hmm",
                        required=True,
                        help="HMM file to build")
    parser.add_argument("-t", "--tax_ids",
                        required=True,
                        help="The tax_ids file created by treesapp create")
    parser.add_argument("-r", "--rank")
    parser.add_argument("--he")
    parser.add_argument("--me")
    parser.add_argument("-n", "--num_threads",
                        required=False, default=2)
    parser.add_argument("-x", "--tmp_dir")

    return parser.parse_args()


def make_dereplicated_hmm(aln_file: str, taxonomic_ids: str, dereplication_rank: str,
                          hmmbuild: str, mafft: str, hmm_file: str, n_threads=2, intermediates_dir=None) -> None:
    """
    Function to create a taxonomically-dereplicated hidden Markov model (HMM) profile. This reduces the bias from
    potentially over-represented clades, increasing the weight of their conserved positions.

    :param aln_file: Path to a reference multiple sequence alignment file made by *treesapp create* as part of a refpkg
    :param taxonomic_ids: Path to a tax_ids file created by *treesapp create* as part of a refpkg
    :param dereplication_rank: The taxonomic rank to dereplicate the reference sequences at
    :param hmmbuild: Path to an hmmbuild executable
    :param mafft: Path to a MAFFT executable
    :param hmm_file: Name of the HMM file to be written
    :param n_threads: Number of threads MAFFT should use
    :param intermediates_dir: A path to a directory to write intermediate files
    :return: None
    """

    logging.info("Creating taxonomically-dereplicated HMM... ")

    aln_pattern_match = re.match(r"(\w+).fa", os.path.basename(aln_file))
    tax_ids_pattern_match = re.match(r"^tax_ids_(\w+).txt", os.path.basename(taxonomic_ids))

    if aln_pattern_match.group(1) == tax_ids_pattern_match.group(1):
        refpkg_name = aln_pattern_match.group(1)
    elif not aln_pattern_match:
        logging.error("File name format is unexpected for '%s'.\n" % os.path.basename(aln_file))
        sys.exit(3)
    elif not tax_ids_pattern_match:
        logging.error("File name format is unexpected for '%s'.\n" % os.path.basename(taxonomic_ids))
        sys.exit(3)
    else:
        logging.error("File names suggest '%s' and '%s' represent different reference packages.\n" %
                      (os.path.basename(aln_file), os.path.basename(taxonomic_ids)))
        sys.exit()

    if not intermediates_dir:
        intermediates_dir = os.path.dirname(aln_file)
    if intermediates_dir[-1] != os.sep:
        intermediates_dir += os.sep
    derep_fa = intermediates_dir + '.'.join(os.path.basename(aln_file).split('.')[:-1]) + "_derep.fa"
    derep_aln = intermediates_dir + '.'.join(os.path.basename(aln_file).split('.')[:-1]) + "_derep.mfa"
    intermediates = [derep_aln, derep_fa]

    lineage_reps = []
    t = {}

    mfa = fasta.FASTA(aln_file)
    mfa.load_fasta()

    leaf_nodes = tax_ids_file_to_leaves(taxonomic_ids)
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
    make_dereplicated_hmm(args.input_fa, args.tax_ids, args.rank, args.he, args.me, args.output_hmm,
                          args.num_threads, args.tmp_dir)
