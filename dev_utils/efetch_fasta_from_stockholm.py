#!/usr/bin/env python3

import time
import re
import sys
import argparse
from urllib import error
from Bio import Entrez

__author__ = 'Connor Morgan-Lang'


def get_options():
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--stockholm",
                        required=True,
                        help="Stockholm file containing NCBI accession IDs to be downloaded")
    parser.add_argument("-m", "--molecule",
                        required=True,
                        help="The molecule type of the sequences to be downloaded. "
                             "This effects which database is queried, "
                             "in turn effecting the success of this run...",
                        choices=["protein", "nucleotide"])
    parser.add_argument("-o", "--output",
                        required=False,
                        default="stockholm_downloads.fasta",
                        help="FASTA file to write the sequences to. [DEFAULT=stockholm_downloads.fasta]")
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        default=False,
                        help="Print what is happening at every stage.")

    args = parser.parse_args()

    return args


def read_stockholm(args):
    accessions = list()
    try:
        stklm_handler = open(args.stockholm, 'r')
    except IOError:
        sys.stderr.write("ERROR: Unable to open " + args.stockholm + " for reading!\n")
        sys.exit(1)

    if args.verbose:
        sys.stdout.write("\tReading stockholm file... ")
        sys.stdout.flush()

    stockholm_header_re = re.compile("^#=GS (.*)/([0-9])+-([0-9])+\w+.*")
    # We are able to include the start and end coordinates for downstream sequence slicing
    line = stklm_handler.readline()
    while not stockholm_header_re.match(line):
        line = stklm_handler.readline()
    while stockholm_header_re.match(line):
        accessions.append(stockholm_header_re.match(line).group(1))
        line = stklm_handler.readline()

    stklm_handler.close()

    if args.verbose:
        sys.stdout.write("done.\n")
        sys.stdout.write(str(len(accessions)) + " accessions parsed from " + args.stockholm + "\n")

    return accessions


def setup_progress_bar(num_items):
    if num_items > 100:
        progress_bar_width = 100
        step_proportion = float(num_items) / progress_bar_width
    else:
        progress_bar_width = num_items
        step_proportion = 1

    sys.stdout.write("[%s ]" % (" " * progress_bar_width))
    sys.stdout.write("\b" * (progress_bar_width + 3))
    sys.stdout.flush()

    return step_proportion


def fetch_sequences(args, accessions):
    fasta_string = ""
    Entrez.email = "A.N.Other@example.com"

    if args.verbose:
        sys.stdout.write("\tTesting Entrez efetch and your internet connection... ")
        sys.stdout.flush()
    try:
        Entrez.efetch(db="Taxonomy", id="158330", retmode="xml")
    except error.URLError:
        raise AssertionError("ERROR: Unable to serve Entrez query. Are you connected to the internet?")

    if args.verbose:
        sys.stdout.write("done.\n")
        sys.stdout.write("\tFetching sequences from Entrez. Here is a sweet progress bar:\n")
        sys.stdout.flush()

    step_proportion = setup_progress_bar(len(accessions))
    acc = 0.0

    for accession in accessions:
        if not accession:
            sys.stderr.write("WARNING: Blank accession ID found for unknown sequence.\n")
        handle = Entrez.efetch(db=args.molecule, id=str(accession), rettype="fasta", retmode="text")
        fasta_seq = handle.read()
        if fasta_seq[0] == '>':
            fasta_string += fasta_seq.strip() + "\n"
        else:
            sys.stderr.write("Something bad happened when fetching " + accession + " from " + args.molecule + ".\n")
            sys.exit(4)

        # Update the progress bar
        acc += 1.0
        if acc >= step_proportion:
            acc -= step_proportion
            time.sleep(0.1)
            sys.stdout.write("-")
            sys.stdout.flush()

    sys.stdout.write("-]\n")

    return fasta_string


def write_fasta(args, fasta_dict):
    try:
        fa_out_handler = open(args.output, 'w')
    except IOError:
        sys.stderr.write("ERROR: Unable to open " + args.output + " for writing!\n")
        sys.exit(5)

    fa_out_handler.write(fasta_dict)
    fa_out_handler.close()
    return


def main():
    args = get_options()
    accessions = read_stockholm(args)
    fasta_dict = fetch_sequences(args, accessions)
    write_fasta(args, fasta_dict)


main()
