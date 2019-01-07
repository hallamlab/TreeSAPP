#!/usr/bin/env python3

import re
import sys
import os
import argparse
import inspect
import logging
from urllib import error
from Bio import Entrez
from tqdm import tqdm

cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile(inspect.currentframe()))[0]))
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)
sys.path.insert(0, cmd_folder + os.sep + ".." + os.sep)
from classy import prep_logging
from fasta import get_headers, get_header_format
from utilities import return_sequence_info_groups, complement_nucs

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
        logging.error("Unable to open " + args.input + " for reading!\n")
        sys.exit(1)

    logging.debug("Reading accession list file... ")
    line = list_handler.readline()
    while line:
        if line[0] == '>':
            line = line.strip()[1:]
        else:
            line = line.strip()
        accessions.append(line)
        line = list_handler.readline()
        
    list_handler.close()

    logging.debug("done.\n")
    logging.info(str(len(accessions)) + " accessions parsed from " + args.input + "\n")

    return accessions


def read_stockholm(args):
    accessions = list()
    try:
        stklm_handler = open(args.input, 'r')
    except IOError:
        logging.error("Unable to open " + args.input + " for reading!\n")
        sys.exit(1)

    logging.debug("Reading stockholm file... ")

    stockholm_header_re = re.compile("^#=GS (.*)/([0-9])+-([0-9])+\w+.*")
    # We are able to include the start and end coordinates for downstream sequence slicing
    line = stklm_handler.readline()
    while not stockholm_header_re.match(line):
        line = stklm_handler.readline()
    while stockholm_header_re.match(line):
        accessions.append(stockholm_header_re.match(line).group(1))
        line = stklm_handler.readline()

    stklm_handler.close()

    logging.debug("done.\n")
    logging.debug(str(len(accessions)) + " accessions parsed from " + args.input + "\n")

    return accessions


def grab_coding_sequence_entrez(cds):
    try:
        seq_name, positions = cds.split(':')
        start, end = positions.split("..")
    except ValueError:
        logging.error("Unexpected formatting of GBQualifier_value:\n\t" + cds)
        sys.exit(11)
    try:
        handle = Entrez.efetch(db="nucleotide", id=str(seq_name), rettype="fasta")
        scf = handle.read()
    except error.HTTPError:
        return ""

    long_scf = ''.join(scf.split("\n")[1:])
    try:
        return long_scf[(int(start) - 1):int(end)]
    except ValueError:
        # logging.warning("ValueError for CDS: '" + cds + "'\n\t" +
        #                 "start position = " + start +
        #                 ", end position = " + end + "\n")
        start = int(re.sub("[<>]", '', start))
        end = int(re.sub("[<>]", '', end))
        return long_scf[(start - 1):end]


def parse_entrez_xml(args, xml_string):
    """
    This function is only ever used if the input molecule type is not identical to the desired output molecule type.

    Parse the xml string returned by Entrez and return either a list or empty string if the query was unsuccessful
    or a string with the format:

     >accession [organism]
     sequence
    :param args:
    :param xml_string:
    :return:
    """
    if not xml_string:
        return ""
    sequence = ""
    organism = ""
    accession = ""

    if 'GBSeq_feature-table' in xml_string.keys():
        for feature in xml_string['GBSeq_feature-table']:
            for feature_key in feature:
                if feature_key == "GBFeature_key":
                    if feature[feature_key] == "CDS":
                        for element in feature["GBFeature_quals"]:
                            if args.molecule_in == "protein" and element['GBQualifier_name'] == "coded_by":
                                cds = element['GBQualifier_value']
                                if re.search("complement|join", cds):
                                    tmp_seq = ""
                                    instructions = cds.split('(')[:-1]
                                    loci = cds.strip(')').split('(')[-1]
                                    for locus in loci.split(','):
                                        tmp_seq += grab_coding_sequence_entrez(locus)
                                    if "complement" in instructions:
                                        sequence = complement_nucs(tmp_seq)
                                else:
                                    sequence = grab_coding_sequence_entrez(cds)
                            elif element['GBQualifier_name'] == "translation":
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
            return [accession, organism, sequence, xml_string]
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
    chunk_size = 1000
    num_queries = len(accessions)
    fasta_string = ""
    alternative_molecule = ""
    query_acc_list = list()
    failures = list()
    blanks = list()
    Entrez.email = "c.morganlang@gmail.com"

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

    # Make the lists of query accessions
    for i in range(0, num_queries, chunk_size):
        query_acc_list.append(", ".join(accessions[i:i+chunk_size]))

    for acc_list_chunk in tqdm(query_acc_list):
        if not acc_list_chunk:
            logging.warning("Blank accession chunk encountered.\n")
            continue

        if args.molecule_in == args.seq_out:
            # TODO: test this
            try:
                handle = Entrez.efetch(db=args.molecule_in,
                                       id=str(acc_list_chunk),
                                       rettype="fasta",
                                       retmode="text")
            except error.HTTPError:
                # Accession ID is for another database, so try and find the correct accession
                try:
                    genbank_handle = Entrez.efetch(db=alternative_molecule,
                                                   id=str(acc_list_chunk),
                                                   rettype="gb",
                                                   retmode="text")
                    gb_record = genbank_handle.read()
                except error.HTTPError:
                    failures.append(acc_list_chunk)
                    continue
                if args.molecule_in == "protein":
                    acc_list_chunk = re.search("protein_id=\"(.*)\"", gb_record).group(1)
                else:
                    print("New territory. Figure this out!")
                    print(gb_record)
                handle = Entrez.efetch(db=args.molecule_in, id=str(acc_list_chunk), rettype="fasta", retmode="text")
            fasta_seq = handle.read()
            if fasta_seq == "":
                blanks.append(acc_list_chunk)
                continue
        else:
            try:
                handle = (Entrez.efetch(db=args.molecule_in, id=str(acc_list_chunk), retmode="xml"))
                records = Entrez.read(handle)
            except error.HTTPError:
                logging.error("HTTPError! Here is the id used for this failed query:\n" + acc_list_chunk + "\n")
                sys.exit(3)
            for xml_string in records:
                # TODO: group the records into a single query (involving sig. changes to parse_entrez_xml) for speed
                fasta_seq = parse_entrez_xml(args, xml_string)
                if type(fasta_seq) is list:
                    logging.debug("Unable to find the converted acc_list_chunk or sequence for " + fasta_seq[0] + "\n" +
                                  "\tAccession: " + fasta_seq[0] + "\n" +
                                  "\tOrganism: " + fasta_seq[1] + "\n" +
                                  "\tSequence: " + str(fasta_seq[2]) + "\n" +
                                  "\tRecord: " + str(fasta_seq[3]) + "\n")
                elif fasta_seq == "":
                    blanks.append(fasta_seq[0])

                if fasta_seq and fasta_seq[0] == '>':
                    fasta_string += fasta_seq.strip() + "\n"
                else:
                    failures.append(fasta_seq[0])

    if failures:
        logging.debug("Unable to fetch information from NCBI for " +
                      str(len(failures)) + '/' + str(num_queries) + ":\n" +
                      '\n'.join(failures) + "\n")
    if blanks:
        logging.debug("Entrez server returned empty fasta sequences for " +
                      str(len(blanks)) + '/' + str(num_queries) + ":\n" +
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
        accessions = []
        for header in get_headers(args.input):
            header_format_re, header_db, header_molecule = get_header_format(header)
            sequence_info = header_format_re.match(header)
            accession, _, _, _, _ = return_sequence_info_groups(sequence_info, header_db, header)
            accessions.append(accession)
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
