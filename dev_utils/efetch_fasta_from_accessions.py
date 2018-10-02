#!/usr/bin/env python3

import time
import re
import sys
import os
import argparse
import inspect
import logging
from urllib import error
from Bio import Entrez

cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile(inspect.currentframe()))[0]))
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)
sys.path.insert(0, cmd_folder + os.sep + ".." + os.sep)
from external_command_interface import setup_progress_bar
from classy import prep_logging
from fasta import get_headers

__author__ = 'Connor Morgan-Lang'


def get_options():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--format",
                        required=False, default="list",
                        help="Format of the input file",
                        choices=["list", "stockholm", "fasta"])
    parser.add_argument("-i", "--input",
                        required=True,
                        help="File containing NCBI accession IDs to be downloaded")
    parser.add_argument("-m", "--molecule_in",
                        required=True,
                        help="The molecule type of the accessions. "
                             "This effects which database is queried, "
                             "in turn effecting the success of this run...",
                        choices=["protein", "nucleotide"])
    parser.add_argument("-s", "--seq_out",
                        required=False, default="protein", choices=["protein", "nucleotide"],
                        help="The molecule type of the sequences to be written to the FASTA file")
    parser.add_argument("-o", "--output",
                        required=False,
                        default="entrez_downloads.fasta",
                        help="FASTA file to write the sequences to. [DEFAULT=entrez_downloads.fasta]")
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        default=False,
                        help="Print what is happening at every stage.")

    args = parser.parse_args()

    return args


def read_accession_list(args):
    accessions = list()
    try:
        list_handler = open(args.input, 'r')
    except IOError:
        sys.stderr.write("\tERROR: Unable to open " + args.input + " for reading!\n")
        sys.exit(1)

    if args.verbose:
        sys.stdout.write("Reading accession list file... ")
        sys.stdout.flush()
    line = list_handler.readline()
    while line:
        if line[0] == '>':
            line = line.strip()[1:]
        else:
            line = line.strip()
        accessions.append(line)
        line = list_handler.readline()
        
    list_handler.close()

    if args.verbose:
        sys.stdout.write("done.\n")
        sys.stdout.write(str(len(accessions)) + " accessions parsed from " + args.input + "\n")

    return accessions


def read_stockholm(args):
    accessions = list()
    try:
        stklm_handler = open(args.input, 'r')
    except IOError:
        sys.stderr.write("\tERROR: Unable to open " + args.input + " for reading!\n")
        sys.exit(1)

    if args.verbose:
        sys.stdout.write("Reading stockholm file... ")
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
        sys.stdout.write(str(len(accessions)) + " accessions parsed from " + args.input + "\n")

    return accessions


def parse_entrez_xml(args, xml_string):
    sequence = ""
    organism = ""
    accession = ""
    # Header format will be: > + accession + '_' + GBSeq_organism
    record = Entrez.read(xml_string)
    if args.molecule_in == "protein":
        sys.stderr.write("\tThis is untested :) Here is the information you need:\n")
        print(record)
        sys.exit()
    else:
        if not record:
            return ""
        for sub_info in record:
            if 'GBSeq_feature-table' in sub_info.keys():
                for feature in sub_info['GBSeq_feature-table']:
                    for feature_key in feature:
                        if feature_key == "GBFeature_key":
                            if feature[feature_key] == "CDS":
                                for element in feature["GBFeature_quals"]:
                                    if element['GBQualifier_name'] == "translation":
                                        sequence = element['GBQualifier_value']
                                    else:
                                        pass
                            elif feature[feature_key] == "source":
                                for element in feature["GBFeature_quals"]:
                                    if element['GBQualifier_name'] == "organism":
                                        organism = element['GBQualifier_value']
                                    else:
                                        pass
                        elif feature_key == 'GBFeature_intervals' and accession == "":
                            for element in feature['GBFeature_intervals']:
                                accession = element['GBInterval_accession']
                if sequence == "" or accession == "":
                    return [accession, organism, sequence, record]
                else:
                    return '>' + accession + ' [' + organism + "]\n" + sequence
            else:
                pass
    return ""


def fetch_sequences(args, accessions):
    """

    :param args:
    :param accessions:
    :return:
    """
    fasta_string = ""
    alternative_molecule = ""
    failures = list()
    blanks = list()
    Entrez.email = "A.N.Other@example.com"

    logging.debug("Testing Entrez efetch and your internet connection... ")

    try:
        Entrez.efetch(db="Taxonomy", id="158330", retmode="xml")
    except error.URLError:
        logging.error("Unable to serve Entrez query. Are you connected to the internet?\n")
        sys.exit(11)

    logging.debug("done.\n")
    logging.info("Fetching sequences from Entrez.\n")

    for mol in ["protein", "nucleotide"]:
        if mol != args.molecule_in:
            alternative_molecule = mol
            break

    step_proportion = setup_progress_bar(len(accessions))
    acc = 0.0

    for accession in accessions:
        if not accession:
            logging.warning("Blank accession ID found for unknown sequence.\n")

        if args.molecule_in == args.seq_out:
            try:
                handle = Entrez.efetch(db=args.molecule_in,
                                       id=str(accession),
                                       rettype="fasta",
                                       retmode="text")
            except error.HTTPError:
                # Accession ID is for another database, so try and find the correct accession
                try:
                    genbank_handle = Entrez.efetch(db=alternative_molecule,
                                                   id=str(accession),
                                                   rettype="gb",
                                                   retmode="text")
                    gb_record = genbank_handle.read()
                except error.HTTPError:
                    failures.append(accession)
                    acc += 1.0
                    continue
                if args.molecule_in == "protein":
                    accession = re.search("protein_id=\"(.*)\"", gb_record).group(1)
                else:
                    print("New territory. Figure this out!")
                    print(gb_record)
                handle = Entrez.efetch(db=args.molecule_in, id=str(accession), rettype="fasta", retmode="text")
            fasta_seq = handle.read()
            if fasta_seq == "":
                blanks.append(accession)
                acc += 1.0
                continue
        else:
            handle = (Entrez.efetch(db=args.molecule_in, id=str(accession), retmode="xml"))
            fasta_seq = parse_entrez_xml(args, handle)
            if type(fasta_seq) is list:
                logging.warning("Unable to find the converted accession or sequence for " + accession + "\n")
                logging.debug("\tAccession: " + fasta_seq[0] + "\n" +
                              "\tOrganism: " + fasta_seq[1] + "\n" +
                              "\tSequence: " + str(fasta_seq[2]) + "\n" +
                              "\tRecord: " + str(fasta_seq[3]) + "\n")
            elif fasta_seq == "":
                blanks.append(accession)
                acc += 1.0
                continue

        if fasta_seq and fasta_seq[0] == '>':
            fasta_string += fasta_seq.strip() + "\n"
        else:
            acc += 1.0
            failures.append(accession)

        # Update the progress bar
        acc += 1.0
        if acc >= step_proportion:
            acc -= step_proportion
            time.sleep(0.1)
            sys.stdout.write("-")
            sys.stdout.flush()

    sys.stdout.write("-]\n")

    logging.debug("Unable to fetch information from NCBI:\n" +
                  '\n'.join(failures) + "\n")
    logging.debug("Entrez server returned empty fasta sequences:\n" +
                  "\n".join(blanks) + "\n")

    return fasta_string


def write_fasta(args, fasta_dict):
    try:
        fa_out_handler = open(args.output, 'w')
    except IOError:
        logging.error("Unable to open " + args.output + " for writing!\n")
        sys.exit(11)

    fa_out_handler.write(fasta_dict)
    fa_out_handler.close()
    return


def main():
    args = get_options()
    log_file_name = "Efetch_fasta_log.txt"
    prep_logging(log_file_name, args.verbose)
    if args.format == "stockholm":
        accessions = read_stockholm(args)
    elif args.format == "fasta":
        # Need to remove the '>'
        accessions = [header[1:] for header in get_headers(args.input)]
    elif args.format == "list":
        accessions = read_accession_list(args)
    else:
        logging.error("Unrecognized file format '" + args.format + "'.\n")
        sys.exit(11)

    if len(accessions) == 0:
        logging.error("No accessions were read from '" + args.input + "'.\n")
        sys.exit(11)
    fasta_dict = fetch_sequences(args, accessions)
    write_fasta(args, fasta_dict)


main()
