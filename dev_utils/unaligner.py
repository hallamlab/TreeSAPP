#!/usr/bin/env python3

import sys
import os
import argparse
import inspect
import re
import logging

cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile(inspect.currentframe()))[0]))
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)
sys.path.insert(0, cmd_folder + os.sep + ".." + os.sep)
from classy import prep_logging
from fasta import read_fasta_to_dict, write_new_fasta

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
    fasta_dict = read_fasta_to_dict(args.input_fa)
    logging.info("done.\n")
    logging.info("Stripping alignment characters from sequences... ")
    unaligned_dict = strip_alignment(fasta_dict)
    logging.info("done.\n")
    logging.info("Writing unaligned sequences to " + args.output_fa + "... ")
    write_new_fasta(unaligned_dict, args.output_fa)
    logging.info("done.\n")
    logging.info("The unaligner has finished successfully!\n")


main()
