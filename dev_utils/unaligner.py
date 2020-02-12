#!/usr/bin/env python3

import os
import argparse
import re
import logging

from treesapp.classy import prep_logging
from treesapp.fasta import write_new_fasta, FASTA

__author__ = 'Connor Morgan-Lang'


def get_options():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_fa",
                        required=True,
                        help="Multiple alignment FASTA file.")
    parser.add_argument("-o", "--output_fa",
                        required=False,
                        help="FASTA file to write the sequences to. [DEFAULT=input_unaligned.fasta]")
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        default=False,
                        help="Print what is happening at every stage.")

    args = parser.parse_args()
    if not args.output_fa:
        input_prefix = '.'.join(args.input_fa.split('.')[:-1])
        ext = args.input_fa.split('.')[-1]
        args.output_fa = input_prefix + "_unaligned." + ext

    return args


def strip_alignment(fasta_dict):
    unaligned_dict = {}
    for header in fasta_dict:
        seq = fasta_dict[header]
        unaligned_dict[header] = re.sub("[-.]", '', seq)
    return unaligned_dict


def main():
    args = get_options()
    prep_logging(os.path.dirname(args.output_fa) + os.sep + "unaligner_log.txt", args.verbose)
    logging.info("Reading FASTA file '" + args.input_fa + "'... ")
    msa_fa = FASTA(args.input_fa)
    msa_fa.load_fasta()
    logging.info("done.\n")

    logging.info("Stripping alignment characters from sequences... ")
    msa_fa.unalign()
    logging.info("done.\n")
    msa_fa.change_dict_keys()

    logging.info("Writing unaligned sequences to " + args.output_fa + "... ")
    write_new_fasta(msa_fa.fasta_dict, args.output_fa)
    logging.info("done.\n")

    logging.info("The unaligner has finished successfully!\n")


main()
