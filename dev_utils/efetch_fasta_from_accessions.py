#!/usr/bin/env python3

import re
import sys
import argparse
import logging
from time import sleep, time
from urllib import error
from Bio import Entrez
from tqdm import tqdm
from treesapp.classy import prep_logging
from treesapp.fasta import get_headers, get_header_format, load_fasta_header_regexes, sequence_info_groups
from treesapp.utilities import complement_nucs

__author__ = 'Connor Morgan-Lang'


def get_options():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input",
                        required=True,
                        help="File containing NCBI accession IDs to be downloaded")
    parser.add_argument("-f", "--format",
                        required=False, default="list",
                        help="Format of the input file",
                        choices=["list", "stockholm", "fasta"])
    parser.add_argument("-m", "--molecule_in",
                        required=False, default="protein", choices=["protein", "nucleotide"],
                        help="The molecule type of the accessions. "
                             "This effects which database is queried,  in turn effecting the success of this run... "
                             "[DEFAULT=protein]")
    parser.add_argument("-s", "--seq_out",
                        required=False, default="protein", choices=["protein", "nucleotide"],
                        help="The molecule type of the sequences to be written to the FASTA file. "
                             "[DEFAULT=protein]")
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


class SeqRecord:
    def __init__(self):
        self.prot_acc = ""
        self.nuc_acc = ""
        self.organism = ""
        self.sequence = ""
        self.positions = []
        self.instructions = []

    def finish_sequence(self):
        if not self.sequence:
            logging.warning("No sequence found for:\n" + self.get_info())
            return

        if self.positions:
            seq_chunks = []
            for coords in self.positions:
                if not re.search('\.\.', coords):
                    logging.warning("Invalid range '" + coords + "' encountered in record:" + "'\n" + self.get_info())
                    continue
                try:
                    start, end = coords.split("..")
                except ValueError:
                    logging.warning("Bad coords '" + coords + "' slipped through:\n" + self.get_info() + "\n")
                    sys.exit()
                try:
                    seq_chunks.append(self.sequence[(int(start) - 1):int(end)])
                except ValueError:
                    # logging.warning("ValueError for CDS: '" + cds + "'\n\t" +
                    #                 "start position = " + start +
                    #                 ", end position = " + end + "\n")
                    start = int(re.sub("[<>]", '', start))
                    end = int(re.sub("[<>]", '', end))
                    seq_chunks.append(self.sequence[(start - 1):end])
            self.sequence = ''.join(seq_chunks)

        if self.instructions:
            for operation in self.instructions:
                if operation == "complement":
                    self.sequence = complement_nucs(self.sequence)
                elif operation == "join":
                    pass  # This is done first
                else:
                    logging.warning("Unknown instruction '" + str(operation) + "'\n")
        return

    def get_desired_accession(self, seq_out):
        if seq_out == "protein":
            return self.prot_acc
        elif seq_out == "nucleotide":
            return self.nuc_acc
        else:
            logging.error("Unknown seq_out parameter '" + seq_out + "'")
            sys.exit(3)

    def get_info(self):
        return "\tProtein accession: " + self.prot_acc + "\n" + \
               "\tNucleotide accession: " + self.nuc_acc + "\n" + \
               "\tOrganism: " + self.organism + "\n" + \
               "\tSequence positions: " + ", ".join(self.positions) + "\n" + \
               "\tSequence instructions: " + ", ".join(self.instructions) + "\n" + \
               "\tSequence is " + str(self.sequence.isalnum()) + " and length is " + str(len(self.sequence)) + "\n"

    def fastafy(self, seq_out):
        header = self.get_desired_accession(seq_out)
        if not header:
            return None
        if header[0] != '>':
            header = '>' + header
        if self.organism:
            header += " [" + self.organism + "]"
        if self.sequence:
            return {header: self.sequence}
        else:
            logging.warning("No sequence available for:" + self.get_info() + "\n")
            return None


def read_accession_list(args):
    accessions = list()
    try:
        list_handler = open(args.input, 'r')
    except IOError:
        logging.error("Unable to open " + args.input + " for reading!\n")
        sys.exit(1)

    logging.debug("Reading accession list file... ")
    for line in list_handler:
        if line[0] == '>':
            line = line.strip()[1:]
        else:
            line = line.strip()
        accessions.append(line)

    list_handler.close()

    logging.debug("done.\n")
    logging.info(str(len(accessions)) + " accessions parsed from " + args.input + "\n")
    uniq_accs = set(accessions)
    logging.info(str(len(uniq_accs)) + " unique accessions " + args.input + "\n")

    return uniq_accs


def read_stockholm(args):
    accessions = list()
    try:
        stklm_handler = open(args.input, 'r')
    except IOError:
        logging.error("Unable to open " + args.input + " for reading!\n")
        sys.exit(1)

    logging.debug("Reading stockholm file... ")

    stockholm_header_re = re.compile(r"^#=GS (.*)/([0-9])+-([0-9])+\w+.*")
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
    uniq_accs = set(accessions)
    logging.info(str(len(uniq_accs)) + " unique accessions " + args.input + "\n")

    return uniq_accs


def parse_entrez_xml(molecule_in: str, xml_string):
    """
    This function is only ever used if the input molecule type is not identical to the desired output molecule type.

    Parse the xml string returned by Entrez and return a list containing the accession of the desired sequence,
    the organism, instructions for generating the final sequence, and potentially the sequence itself (for nuc -> prot)

    :param molecule_in: Either 'protein' or 'nucleotide'
    :param xml_string: A single Entrez xml record
    :return: [prot_accession, nuc_accession, organism, positions, instructions, sequence]
    """
    if not xml_string:
        return ""
    organism = ""
    prot_acc = ""
    nuc_acc = ""
    sequence = ""
    positions = []
    instructions = []

    if 'GBSeq_feature-table' in xml_string.keys():
        for feature in xml_string['GBSeq_feature-table']:
            for feature_key in feature:
                if feature_key == "GBFeature_key":
                    if feature[feature_key] == "CDS":
                        for element in feature["GBFeature_quals"]:
                            if molecule_in == "protein" and element['GBQualifier_name'] == "coded_by":
                                cds = element['GBQualifier_value']
                                if re.search("(|)", cds):
                                    instructions = cds.split('(')[:-1]
                                cds = re.sub(r".*\(", '', re.sub("\)", '', cds))
                                for piece in cds.split(','):
                                    nuc_acc, piece_pos = piece.split(':')
                                    positions.append(piece_pos)
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
                elif feature_key == 'GBFeature_intervals':
                    for element in feature['GBFeature_intervals']:
                        prot_acc = element['GBInterval_accession']
        return [prot_acc, nuc_acc, organism, positions, instructions, sequence]
    else:
        pass
    return ""


def fetch_sequences(args, accessions: set):
    """

    :param args:
    :param accessions:
    :return:
    """
    alternative_molecule = ""
    primary_accs = set()
    secondary_accs = set()
    blanks = 0
    seq_record_list = dict()
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
    accessions = sorted(list(accessions))
    num_queries = len(accessions)
    chunk_size = int(round(num_queries/100)) + 1
    for i in range(0, num_queries, chunk_size):
        primary_accs.add(", ".join(accessions[i:i+chunk_size]))

    for acc_list_chunk in tqdm(sorted(primary_accs)):
        t_start = time()
        if not acc_list_chunk:
            logging.warning("Blank accession chunk encountered.\n")
            continue

        if args.molecule_in == args.seq_out:
            try:
                handle = Entrez.efetch(db=args.molecule_in,
                                       id=str(acc_list_chunk),
                                       rettype="fasta",
                                       retmode="text")
                # Load the current batch of SequenceRecords in acc_list_chunk into seq_record_list
                logging.debug("Reading... ")
                fasta_seq = handle.read()
                handle.close()
                logging.debug("done.\n")
            except error.HTTPError:
                # Accession ID is for another database, so try and find the correct accession
                fasta_seq = ""
                for chunk_acc in acc_list_chunk.split(", "):
                    # There appears to be a bad egg in here - query them individually
                    try:
                        gb_handle = Entrez.efetch(db=args.molecule_in,
                                                  id=str(chunk_acc),
                                                  rettype="gb",
                                                  retmode="text")
                        gb_record = gb_handle.read()
                    except error.HTTPError:
                        # Alrighty, maybe we're looking in the wrong database?
                        try:
                            gb_handle = Entrez.efetch(db=alternative_molecule,
                                                      id=str(chunk_acc),
                                                      rettype="gb",
                                                      retmode="text")
                            gb_record = gb_handle.read()
                        except error.HTTPError:
                            blanks += 1
                            continue
                    if args.molecule_in == "protein":
                        chunk_acc = re.search("protein_id=\"(.*)\"", gb_record).group(1)
                    else:
                        print("New territory. Figure this out!")
                        print(gb_record)
                    handle = Entrez.efetch(db=args.molecule_in, id=str(chunk_acc),
                                           rettype="fasta", retmode="text")
                    fasta_seq += handle.read()

            fasta_dict = fasta_string_to_dict(fasta_seq, divide=False)
            for acc_out in sorted(fasta_dict.keys()):
                sr = SeqRecord()
                if args.seq_out == "protein":
                    sr.prot_acc = acc_out
                elif args.seq_out == "nucleotide":
                    sr.nuc_acc = acc_out
                sr.sequence = fasta_dict[acc_out]
                seq_record_list[acc_out] = sr
        else:
            # This method requires two separate Entrez queries and so is a bit different
            logging.debug("Downloading primary... ")
            try:
                handle = (Entrez.efetch(db=args.molecule_in, id=str(acc_list_chunk), retmode="xml"))
                records = Entrez.read(handle)
            except error.HTTPError:
                # This is for debugging - the query failed for one of the accessions, but which one?
                logging.error("HTTPError! Here are the failed queries:\n" + acc_list_chunk + "\n")
                for failure in acc_list_chunk.split(", "):
                    sys.stderr.write("Now testing download of: " + failure + "\n")
                    Entrez.efetch(db=args.seq_out, id=failure, retmode="xml")
                sys.exit(3)
            logging.debug("done.\n")

            # Parse the individual Entrez xml-formatted records generated from the input-molecule
            for xml_string in records:
                sr = SeqRecord()
                distilled_record = parse_entrez_xml(args.molecule_in, xml_string)
                if distilled_record:
                    sr.prot_acc, sr.nuc_acc, sr.organism, sr.positions, sr.instructions, sr.sequence = distilled_record
                else:
                    blanks += 1
                    continue
                if not sr.get_desired_accession(args.seq_out):
                    logging.debug("Unable to find the converted accession for:\n" + sr.get_info())
                    blanks += 1
                    continue
                elif not sr.sequence:
                    secondary_accs.add(sr.get_desired_accession(args.seq_out))
                seq_record_list[sr.get_desired_accession(args.seq_out)] = sr

            # Download records for the output-molecule
            if secondary_accs:
                logging.debug("Downloading secondary... ")
                try:
                    handle = Entrez.efetch(db=args.seq_out, id=", ".join(secondary_accs),
                                           rettype="fasta", retmode="text")
                except error.HTTPError:
                    logging.error("HTTPError! Here are the failed queries:\n" + str(secondary_accs) + "\n")
                    # for failure in secondary_accs:
                    #     print(failure)
                    #     handle = Entrez.efetch(db=args.seq_out, id=failure, rettype="fasta", retmode="text")
                    sys.exit(3)
                logging.debug("done.\n")
                secondary_accs.clear()
            else:
                continue
            logging.debug("Reading... ")
            fasta_seq = handle.read()
            handle.close()
            logging.debug("done.\n")
            fasta_dict = fasta_string_to_dict(fasta_seq)
            for acc_out in sorted(fasta_dict.keys()):
                try:
                    sr = seq_record_list[acc_out]
                except KeyError:
                    if acc_out.find('.'):
                        acc_out_trunc, ver = acc_out.split('.')
                        sr = seq_record_list[acc_out_trunc]
                    else:
                        logging.error("ERROR: '" + acc_out + "' not found in SequenceRecords dictionary!\n")
                        sys.exit(1)
                if not sr.sequence:
                    sr.sequence = fasta_dict[acc_out]
                    sr.finish_sequence()
        blanks += len(acc_list_chunk.split(", ")) - len(fasta_dict.keys())
        fasta_dict.clear()

        t_end = time()
        duration = t_end - t_start

        # Sleep for a while to ensure NCBI doesn't slap us with a 'HTTP Error 429: Too Many Requests'
        if duration < 1:
            sleep(1-duration)

    if blanks:
        logging.info("Entrez server returned empty records for " +
                     str(blanks) + '/' + str(num_queries) + "\n")

    return seq_record_list


def write_fasta(args, fasta_dict: dict):
    logging.info("Writing " + str(len(fasta_dict.keys())) + " sequence records to " + args.output + "\n")
    fasta_string = ""
    for header in fasta_dict:
        seq = fasta_dict[header]
        if header[0] != '>':
            header = '>' + header
        fasta_string += header + "\n" + seq + "\n"

    try:
        fa_out_handler = open(args.output, 'w')
    except IOError:
        logging.error("Unable to open " + args.output + " for writing!\n")
        sys.exit(11)

    fa_out_handler.write(fasta_string)
    fa_out_handler.close()
    return


def fasta_string_to_dict(fasta_record, divide=True):
    fasta_dict = dict()
    header = ""
    seq = ""
    logging.debug("Converting FASTA string to a dictionary... ")
    for line in fasta_record.split("\n"):
        if not line:
            continue
        elif line[0] == '>':
            if seq:
                fasta_dict[header] = seq
            seq = ""
            # Just get the accession ID, the first element split on whitespace
            header = line.strip()[1:]
            if divide:
                header = header.split()[0]
        else:
            seq += line
    fasta_dict[header] = seq
    logging.debug("done.\n")
    return fasta_dict


def main():
    args = get_options()
    log_file_name = "Efetch_fasta_log.txt"
    prep_logging(log_file_name, args.verbose)
    header_regexes = load_fasta_header_regexes()
    if args.format == "stockholm":
        accessions = read_stockholm(args)
    elif args.format == "fasta":
        accessions = set()
        input_headers = get_headers(args.input)
        for header in input_headers:
            header_format_re, header_db, header_molecule = get_header_format(header, header_regexes)
            sequence_info = header_format_re.match(header)
            seq_info_tuple = sequence_info_groups(sequence_info, header_db, header, header_regexes)
            accessions.add(seq_info_tuple.accession)
        logging.debug(str(len(input_headers) - len(accessions)) + " duplicate accessions found in " + args.input + "\n")
    elif args.format == "list":
        accessions = read_accession_list(args)
    else:
        logging.error("Unrecognized file format '" + args.format + "'.\n")
        sys.exit(11)

    if len(accessions) == 0:
        logging.error("No accessions were read from '" + args.input + "'.\n")
        sys.exit(11)

    seq_record_list = fetch_sequences(args, accessions)

    # Generate a fasta dictionary from each of the seq_record objects
    fasta_dict = dict()
    failures = list()
    for seq_record in seq_record_list.values():
        seq_dict = seq_record.fastafy(args.seq_out)
        if seq_dict:
            fasta_dict.update(seq_dict)
        else:
            failures.append(seq_record.get_desired_accession(args.molecule_in))

    if failures:
        logging.info("Unable to fetch information from NCBI for " +
                     str(len(failures)) + '/' + str(len(accessions)) + ":\n" +
                     '\n'.join(failures) + "\n")

    write_fasta(args, fasta_dict)


main()
