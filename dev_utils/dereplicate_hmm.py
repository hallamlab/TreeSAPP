#!/usr/bin/env python3

import logging
import os
import argparse

from treesapp import fasta
from treesapp.refpkg import ReferencePackage
from treesapp.utilities import base_file_prefix, which


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--refpkg", required=True,
                        help="Name of the reference package")
    parser.add_argument("-o", "--output_hmm", required=True,
                        help="HMM file to build")
    parser.add_argument("-d", "--refpkg_dir", required=True,
                        help="The directory containing the reference package files")
    parser.add_argument("-r", "--rank", required=False,
                        default="genus",
                        help="Taxonomic rank to dereplicate the reference sequences to. [ DEFAULT = genus ]")
    parser.add_argument("--he", required=False,
                        help="Path to hmmbuild executable.")
    parser.add_argument("--me", required=False,
                        help="Path to MAFFT executable.")
    parser.add_argument("-n", "--num_threads", required=False,
                        default=2,
                        help="Number of threads for MAFFT to use. [ DEFAULT = 2 ]")
    parser.add_argument("-x", "--tmp_dir")

    return parser.parse_args()


def make_dereplicated_hmm(refpkg_name: str, package_path: str, dereplication_rank: str, hmm_file: str,
                          hmmbuild=None, mafft=None, n_threads=2, intermediates_dir=None) -> None:
    """
    Function to create a taxonomically-dereplicated hidden Markov model (HMM) profile. This reduces the bias from
    potentially over-represented clades, increasing the weight of their conserved positions.

    :param refpkg_name: Short-form gene/protein name of the reference package (e.g. McrA, DsrAB)
    :param package_path: Path to directory containing the reference package JSON
    :param dereplication_rank: The taxonomic rank to dereplicate the reference sequences at
    :param hmmbuild: Path to an hmmbuild executable
    :param mafft: Path to a MAFFT executable
    :param hmm_file: Name of the HMM file to be written
    :param n_threads: Number of threads MAFFT should use
    :param intermediates_dir: A path to a directory to write intermediate files
    :return: None
    """
    # TODO: Change command-line parameters to include fasta file

    refpkg = ReferencePackage(refpkg_name)
    refpkg.change_file_paths(package_path)
    refpkg.slurp()

    if not hmmbuild:
        hmmbuild = which("hmmbuild")
    if not mafft:
        mafft = which("mafft")

    refpkg.dereplicate_hmm(dereplication_rank, hmmbuild, mafft, n_threads, intermediates_dir)

    return


if __name__ == "__main__":
    args = get_arguments()
    make_dereplicated_hmm(args.refpkg, args.refpkg_dir, args.rank, args.he, args.me, args.output_hmm,
                          args.num_threads, args.tmp_dir)
