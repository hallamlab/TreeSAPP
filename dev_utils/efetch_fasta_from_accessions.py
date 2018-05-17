#!/usr/bin/env python3

import time
import re
import sys
import os
import argparse
import inspect
from urllib import error
from Bio import Entrez

cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile(inspect.currentframe()))[0]))
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)
sys.path.insert(0, cmd_folder + os.sep + ".." + os.sep)
from external_command_interface import setup_progress_bar

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
    Entrez.email = "A.N.Other@example.com"

    if args.verbose:
        sys.stdout.write("Testing Entrez efetch and your internet connection... ")
        sys.stdout.flush()
    try:
        Entrez.efetch(db="Taxonomy", id="158330", retmode="xml")
    except error.URLError:
        raise AssertionError("ERROR: Unable to serve Entrez query. Are you connected to the internet?")

    if args.verbose:
        sys.stdout.write("done.\n")
        sys.stdout.write("Fetching sequences from Entrez. Here is a sweet progress bar:\n")
        sys.stdout.flush()

    for mol in ["protein", "nucleotide"]:
        if mol != args.molecule_in:
            alternative_molecule = mol
            break

    step_proportion = setup_progress_bar(len(accessions))
    acc = 0.0
    skipped_seqs = 0

    for accession in accessions:
        if not accession:
            sys.stderr.write("\tWARNING: Blank accession ID found for unknown sequence.\n")

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
                    sys.exit("\nERROR: Unable to fetch information from NCBi for accession: " + str(accession) + "\n")
                if args.molecule_in == "protein":
                    accession = re.search("protein_id=\"(.*)\"", gb_record).group(1)
                else:
                    print("New territory. Figure this out!")
                    print(gb_record)
                handle = Entrez.efetch(db=args.molecule_in, id=str(accession), rettype="fasta", retmode="text")
            fasta_seq = handle.read()
            if fasta_seq == "":
                sys.stderr.write("\tEntrez server returned empty fasta sequence for "
                                 + accession + " from " + args.molecule_in + " database.\n")
                skipped_seqs += 1
        else:
            handle = (Entrez.efetch(db=args.molecule_in, id=str(accession), retmode="xml"))
            fasta_seq = parse_entrez_xml(args, handle)
            if type(fasta_seq) is list:
                sys.stderr.write("\tWARNING: Unable to find the converted accession or sequence for " + accession + "\n")
                if args.verbose:
                    sys.stderr.write("\tAccession: " + fasta_seq[0] + "\n")
                    sys.stderr.write("\tOrganism: " + fasta_seq[1] + "\n")
                    sys.stderr.write("\tSequence: " + str(fasta_seq[2]) + "\n")
                    sys.stderr.write("\tRecord: " + str(fasta_seq[3]) + "\n")
            elif fasta_seq == "":
                sys.stderr.write("Entrez record for " + accession + " is empty! Skipping...\n")

        if fasta_seq and fasta_seq[0] == '>':
            fasta_string += fasta_seq.strip() + "\n"
        else:
            sys.stderr.write("\tSomething bad happened when fetching " +
                             accession + " from " + args.molecule_in + ". Skipping...\n")
            skipped_seqs += 1

        # Update the progress bar
        acc += 1.0
        if acc >= step_proportion:
            acc -= step_proportion
            time.sleep(0.1)
            sys.stdout.write("-")
            sys.stdout.flush()

    sys.stdout.write("-]\n")

    sys.stdout.write("Number of sequences skipped = " + str(skipped_seqs) + "\n")

    return fasta_string


def write_fasta(args, fasta_dict):
    try:
        fa_out_handler = open(args.output, 'w')
    except IOError:
        sys.stderr.write("\tERROR: Unable to open " + args.output + " for writing!\n")
        sys.exit(5)

    fa_out_handler.write(fasta_dict)
    fa_out_handler.close()
    return


def main():
    args = get_options()
    if args.format == "stockholm":
        accessions = read_stockholm(args)
    elif args.format == "list":
        accessions = read_accession_list(args)
    else:
        accessions = list()
    fasta_dict = fetch_sequences(args, accessions)
    write_fasta(args, fasta_dict)


main()
