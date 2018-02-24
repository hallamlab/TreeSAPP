#!/usr/bin/env python

import argparse
import re
import sys
from Bio import SeqIO


def get_options():
    parser = argparse.ArgumentParser(description="Used to replace truncated headers (likely just accession identifiers)"
                                                 "with their respective complete headers in a FASTA file.")
    parser.add_argument("--stunted",
                        help="Name of the FASTA file with the stunted headers.",
                        required=True)
    parser.add_argument("--full",
                        help="Name of the FASTA file with the full-length headers.",
                        required=True)
    parser.add_argument("-o", "--output",
                        help="The output FASTA file",
                        required=True)
    parser.add_argument("-l", "--min_length",
                        help="The minimum length of a sequence to be included [ DEFAULT = 1 ]",
                        required=False, type=int, default=1)
    parser.add_argument("-v", "--verbose",
                        action="store_true", required=False, default=True)
    args = parser.parse_args()
    return args


def format_read_fasta(args, fasta_file):
    """
    Reads a FASTA file. Does not check for proper formatting apart from headers beginning with '>'
    :param args: Command-line argument object from get_options and check_parser_arguments
    :param fasta_file: Name of the fasta file to parse
    :return A dictionary with headers as keys and sequences as values
    """
    formatted_fasta_dict = dict()
    too_short = 0

    try:
        fasta_handle = open(fasta_file, 'rU')
    except IOError:
        sys.stderr.write("ERROR: Unable to open FASTA " + fasta_file + " for reading!")
        sys.exit(3)

    for record in SeqIO.parse(fasta_handle, "fasta"):
        if record.__len__() >= args.min_length:
            formatted_fasta_dict[str(record.description)] = str(record.seq)
        else:
            too_short += 1

    fasta_handle.close()

    if too_short > 0:
        sys.stdout.write(str(too_short) + " sequences failed to pass length threshold.\n")

    if args.verbose:
        sys.stdout.write("Parsed " + str(len(formatted_fasta_dict)) + " sequences from " + fasta_file + "\n")
        sys.stdout.flush()

    return formatted_fasta_dict


def find_replace_headers(args, stunted, full):
    new_fa = dict()
    for short_head in stunted.keys():
        found = False
        for header in full.keys():
            if re.search(short_head, header):
                if stunted[short_head] == full[header]:
                    found = True
                    new_fa[header] = stunted[short_head]
                    break
        if not found:
            sys.stderr.write("\tUnable to find match for " + short_head + "\n")

    if args.verbose:
        sys.stdout.write("Matched " + str(len(new_fa.keys())) + "/" + str(len(stunted.keys())) + " headers.\n")
        sys.stdout.flush()

    return new_fa


def write_new_fasta(fasta_dict, fasta_name):
    """
    Function for writing sequences stored in dictionary to file in FASTA format; optional filtering with headers list
    :param fasta_dict: A dictionary containing headers as keys and sequences as values
    :param fasta_name: Name of the FASTA file to write to
    :return:
    """
    try:
        fa_out = open(fasta_name, 'w')
    except IOError:
        raise IOError("Unable to open " + fasta_name + " for writing!")

    for name in fasta_dict.keys():
        seq = fasta_dict[name]
        if name[0] == '>':
            header = name
        else:
            header = '>' + name

        fa_out.write(header + "\n" + seq + "\n")

    fa_out.close()
    return


def main():
    args = get_options()
    stunted_fa = format_read_fasta(args, args.stunted)
    full_fa = format_read_fasta(args, args.full)
    new_fa = find_replace_headers(args, stunted_fa, full_fa)
    write_new_fasta(new_fa, args.output)


main()
