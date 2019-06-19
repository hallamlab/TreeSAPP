#!/usr/bin/env python

import argparse
import sys
import logging
from treesapp.fasta import read_fasta_to_dict, write_new_fasta
from treesapp.classy import prep_logging


def get_options():
    parser = argparse.ArgumentParser(description="Used to replace truncated headers (likely just accession identifiers)"
                                                 " with their respective complete headers in a FASTA file.")
    parser.add_argument("-f", "--old_fa", dest="old_fasta", required=True,
                        help="Name of the FASTA file with the bad/truncated headers.")
    parser.add_argument("-l", "--names", dest="name_map", required=True,
                        help="Path to a csv file with old sequence names in the first field and new ones in the second")
    parser.add_argument("-o", "--new_fa", dest="output", required=True,
                        help="The output FASTA file")
    parser.add_argument("-v", "--verbose",
                        action="store_true", required=False, default=True)
    args = parser.parse_args()
    return args


def find_replace_headers(stunted: dict, header_map: dict) -> dict:
    new_fa = dict()
    for short_head in stunted.keys():
        try:
            new_name = header_map.pop(short_head)
        except KeyError:
            logging.error("Sequence name '" + short_head + "' was not found in the header map.\n")
            sys.exit(5)
        new_fa[new_name] = short_head[stunted]

    logging.debug("Matched " + str(len(new_fa.keys())) + "/" + str(len(stunted.keys())) + " headers.\n")
    logging.debug(str(len(header_map.keys())) + " keys remain unpopped in header map.\n")

    return new_fa


def read_name_map(name_map_file: str) -> dict:
    seq_name_map = dict()
    with open(name_map_file) as name_map_handler:
        for line in name_map_handler:
            try:
                old, new = line.strip().split(',')
            except ValueError:
                logging.error("Line in " + name_map_file + " is not formatted properly:\n" + line)
                sys.exit(3)
            seq_name_map[old] = new
    return seq_name_map


def main():
    args = get_options()
    prep_logging()
    stunted_fa = read_fasta_to_dict(args.old_fasta)
    header_map = read_name_map(args.name_map)
    new_fa = find_replace_headers(stunted_fa, header_map)
    write_new_fasta(new_fa, args.output)


main()
