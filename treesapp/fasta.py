__author__ = 'Connor Morgan-Lang'

import sys
import re
import os
import logging
from time import sleep, time

from math import ceil
import pyfastx
from pyfastxcli import fastx_format_check

from .utilities import median, reformat_string, rekey_dict, return_sequence_info_groups


# No bioinformatic software would be complete without a contribution from Heng Li.
# Adapted from his readfq generator
def generate_fasta(fasta_handler):  # this is a generator function
    last = None  # this is a buffer keeping the last unprocessed line
    while True:  # mimic closure; is it a bad idea?
        if not last:
            for line in fasta_handler:  # search for the start of the next record
                if line[0] == '>':  # fasta header line
                    last = line[:-1]  # save this line
                    break
        if not last:
            break
        name, seqs, last = last[1:], [], None
        for line in fasta_handler:  # read the sequence
            if line[0] == '>':
                last = line[:-1]
                break
            seqs.append(line[:-1])
        if not last or last[0] != '+':  # this is a fasta record
            yield name, ''.join(seqs)  # yield a fasta record
            if not last:
                break
        else:
            seq, seqs = ''.join(seqs), []
            for line in fasta_handler:  # read the quality
                seqs.append(line[:-1])
            if last:  # reach EOF before reading enough quality
                yield name, seq  # yield a fasta record instead
                break


def fastx_split(fastx: str, outdir: str, file_num=1) -> list:
    fastx_type = fastx_format_check(fastx)

    if fastx_type == 'fasta':
        split_files = split_fa(fastx, outdir, file_num)
    elif fastx_type == 'fastq':
        split_files = fq2fa(fastx, outdir, file_num)
    else:
        logging.error("Unknown fastx type: '{}'\n".format(fastx_type))
        sys.exit(3)
    return split_files


def split_seq_writer_helper(fasta_string, fh, seq_write, seqs_num):
    fh.write(fasta_string)
    fasta_string = ""

    # Reduce the number of checks
    if seq_write > seqs_num:
        fh.write(fasta_string)
        fh.close()
        seq_write = 0
    return fasta_string, seq_write


def spawn_new_fasta(file_num, outdir, file_name, digit):
    file_num += 1

    subfile = "{}.{}.{}".format(file_name, str(file_num).zfill(digit), "fasta")
    if outdir is not None:
        subfile = os.path.join(outdir, subfile)

    fh = open(subfile, 'w')
    logging.debug("Writing split FASTA to file '{}'.\n".format(subfile))
    return fh, subfile, file_num


def split_fa(fastx: str, outdir: str, file_num=1, seq_count=0):
    """
    Credit: Lianming Du of pyfastx (v0.6.10)

    Differences are:
1. Format of the output files is always uncompressed FASTA, compatible with ORF prediction tools e.g. Prodigal
2. Batch writing of reads instead of writing a single read at a time to reduce number of  I/O operations

    :param fastx: Path to a FASTQ file to be read and converted to FASTA
    :param outdir: Path to write the output files
    :param file_num: Number of files to split the input FASTQ file into
    :param seq_count: Optionally, instead of splitting into N number of files FASTQ can be split by number of sequences
    :return: List of fasta files generated
    """
    start = time()
    # Fix path with tilde
    fastx = os.path.expanduser(fastx)
    outdir = os.path.expanduser(outdir)
    outputs = []
    fa = pyfastx.Fasta(file_name=fastx, build_index=False)

    # Determine the number of reads in the FASTQ file - faster than building an index
    count = 0
    for _ in fa:
        count += 1

    if file_num > 1:
        seqs_num = ceil(count / file_num)
        parts_num = file_num
    else:
        seqs_num = seq_count
        parts_num = ceil(count / seqs_num)

    file_name, suffix1 = os.path.splitext(os.path.basename(fastx))

    if fa.is_gzip:
        file_name, suffix2 = os.path.splitext(file_name)

    digit = len(str(parts_num))

    seq_write = 0
    max_buffer_size = 1E4
    fasta_string = ""
    fh = None
    file_num = 0

    for name, seq in fa:
        if seq_write == 0:
            fh, subfile, file_num = spawn_new_fasta(file_num, outdir, file_name, digit)
            outputs.append(subfile)

        fasta_string += ">%s\n%s\n" % (name, seq)
        seq_write += 1
        # Batch write
        if seq_write % max_buffer_size == 0:
            fasta_string, seq_write = split_seq_writer_helper(fasta_string, fh, seq_write, seqs_num)

    fh.write(fasta_string)
    fh.close()
    end = time()

    logging.debug("{} completed split_fa in {}s.\n".format(fastx, end - start))
    return outputs


def fq2fa(fastx: str, outdir: str, file_num=1, seq_count=0) -> list:
    """
    Credit: Lianming Du of pyfastx (v0.6.10)
    A modified version of the function fastq_split in the pyfastx python package: https://github.com/lmdu/pyfastx

    Differences are:
1. Format of the output files is always uncompressed FASTA, compatible with ORF prediction tools e.g. Prodigal
2. Batch writing of reads instead of writing a single read at a time to reduce number of  I/O operations

    :param fastx: Path to a FASTQ file to be read and converted to FASTA
    :param outdir: Path to write the output files
    :param file_num: Number of files to split the input FASTQ file into
    :param seq_count: Optionally, instead of splitting into N number of files FASTQ can be split by number of sequences
    :return: List of fasta files generated
    """
    start = time()
    # Fix path with tilde
    fastx = os.path.expanduser(fastx)
    outdir = os.path.expanduser(outdir)
    outputs = []
    fq = pyfastx.Fastq(file_name=fastx, build_index=False)

    # Determine the number of reads in the FASTQ file - faster than building an index
    read_count = 0
    for _ in fq:
        read_count += 1

    if file_num > 1:
        seqs_num = ceil(read_count / file_num)
        parts_num = file_num
    else:
        seqs_num = seq_count
        parts_num = ceil(read_count / seqs_num)

    name, suffix1 = os.path.splitext(os.path.basename(fastx))

    if fq.is_gzip:
        name, suffix2 = os.path.splitext(name)

    digit = len(str(parts_num))

    seq_write = 0
    max_buffer_size = 1E4
    fasta_string = ""
    fh = None
    file_num = 0

    for read_name, seq, _ in fq:
        if seq_write == 0:
            fh, subfile, file_num = spawn_new_fasta(file_num, outdir, name, digit)
            outputs.append(subfile)

        fasta_string += ">%s\n%s\n" % (read_name, seq)
        seq_write += 1
        # Batch write
        if seq_write % max_buffer_size == 0:
            fasta_string, seq_write = split_seq_writer_helper(fasta_string, fh, seq_write, seqs_num)

    fh.write(fasta_string)
    fh.close()
    end = time()

    logging.debug("{} completed fq2fa in {}s.\n".format(fastx, end-start))

    return outputs


def read_fasta_to_dict(fasta_file):
    """
    Reads any fasta file using a generator function (generate_fasta) into a dictionary collection

    :param fasta_file: Path to a FASTA file to be read into a dict
    :return: Dict where headers/record names are keys and sequences are the values
    """
    fasta_dict = dict()
    try:
        fasta_handler = open(fasta_file, 'r')
    except IOError:
        logging.error("Unable to open " + fasta_file + " for reading!\n")
        sys.exit(5)
    for record in generate_fasta(fasta_handler):
        name, sequence = record
        fasta_dict[reformat_string(name)] = sequence.upper()
    return fasta_dict


class Header:
    def __init__(self, header):
        self.original = header
        self.treesapp_num_id = 0
        self.formatted = ""
        self.post_align = ""
        self.first_split = ""
        self.accession = ""

    def get_info(self):
        info_string = "TreeSAPP ID = '" + str(self.treesapp_num_id) + "'\tPrefix = '" + self.first_split + "'\n"
        info_string += "Original =  " + self.original + "\nFormatted = " + self.formatted
        if self.accession:
            info_string += "Accession = " + self.accession + "\n"
        return info_string

    def find_accession(self, refpkg_name=""):
        header_format_re, header_db, header_molecule = get_header_format(self.original, refpkg_name)
        sequence_info = header_format_re.match(self.original)
        self.accession = return_sequence_info_groups(sequence_info, header_db, self.original).accession


def register_headers(header_list, drop=True):
    acc = 1
    header_registry = dict()
    dup_checker = set()
    dups = []
    for header in header_list:
        if drop and header[0] == '>':
            header = header[1:]
        if header in dup_checker:
            dups.append(header)
            continue
        else:
            dup_checker.add(header)
        new_header = Header(header)
        new_header.formatted = reformat_string(header)
        new_header.first_split = header.split()[0]
        new_header.treesapp_num_id = str(acc)
        header_registry[str(acc)] = new_header
        acc += 1

    if len(dups) > 0:
        logging.warning(str(len(dups)) + " duplicate sequence headers found:\n" +
                        ", ".join(dups) + "\n" +
                        "TreeSAPP will proceed as if these duplicate sequences were never seen in 5 seconds or"
                        " you can stop it here with Ctrl-c... the choice is yours.\n")
        sleep(5)

    return header_registry


class FASTA:
    def __init__(self, file_name):
        self.file = file_name
        self.fasta_dict = dict()
        self.header_registry = dict()
        self.amendments = set()  # Set of the TreeSAPP numerical identifiers for all guaranteed sequences
        self.index_form = None

    def load_fasta(self):
        self.fasta_dict = read_fasta_to_dict(self.file)
        self.header_registry = register_headers(get_headers(self.file), True)
        if len(self.fasta_dict) == 0 and len(self.header_registry) == 0:
            logging.error("FASTA file '" + str(self.file) + "' is empty or corrupted - no sequences were found!\n")
            sys.exit(3)

    def add_accession_to_headers(self, refpkg_name=""):
        for acc in self.header_registry:
            header = self.header_registry[acc]  # type: Header
            header.find_accession(refpkg_name)

    def mapping_error(self, bad_headers):
        logging.error("No classified sequences were mapped to '" + self.file + "' FASTA dictionary.\n" +
                      "Here are some example names from the mapping list:\n\t" +
                      "\n\t".join(sorted(bad_headers)[0:6]) + "\n" +
                      "And example names from FASTA dict:\n\t" +
                      "\n\t".join(list(sorted(self.fasta_dict.keys()))[0:6]) + "\n")
        sys.exit(3)

    def n_seqs(self):
        if len(self.header_registry) != len(self.fasta_dict):
            logging.warning("FASTA header registry and fasta dictionary are out of sync!\n")
        return len(self.fasta_dict.keys())

    def original_header_map(self):
        return {self.header_registry[index].original: [self.header_registry[index]] for index in self.header_registry}

    def formatted_header_map(self):
        header_map = dict()
        for index in self.header_registry:
            header = self.header_registry[index]
            f_h = header.formatted
            if f_h in header_map:
                header_map[f_h].append(header)
            else:
                header_map[f_h] = [header]
        return header_map

    def first_split_header_map(self):
        header_map = dict()
        for index in self.header_registry:
            header = self.header_registry[index]
            fs_h = header.first_split
            if fs_h in header_map:
                header_map[fs_h].append(header)
            else:
                header_map[fs_h] = [header]
        return header_map

    def unversion_first_split_header_map(self) -> dict:
        first_splits = self.first_split_header_map()
        accession_header_map = dict()
        for accession in first_splits:
            accession_header_map[re.sub(r'\.\d$', '', accession)] = first_splits[accession]
        return accession_header_map

    def get_seq_names(self, name_format="original") -> list:
        if name_format == "original":
            return [self.header_registry[index].original for index in sorted(self.header_registry, key=int)]
        elif name_format == "first_split":
            return [self.header_registry[index].first_split for index in sorted(self.header_registry, key=int)]
        elif name_format == "formatted":
            return [self.header_registry[index].formatted for index in sorted(self.header_registry, key=int)]
        elif name_format == "num":
            return [index for index in sorted(self.header_registry, key=int)]
        else:
            logging.error("Unrecognized Header format '" + name_format + "'." +
                          " Options are 'original', 'formatted', 'first_split' and 'num'.\n")
            sys.exit(5)

    def keep_only(self, header_subset: list, superset=False):
        """
        Removes all entries from self.fasta_dict and self.header_registry that are not in header_subset.
        :param header_subset: The list of headers found in self.fasta_dict that are to be kept
        :param superset: The header_subset list is a superset so not all headers will be found - do not emit a warning
        :return: None
        """

        if not header_subset:
            logging.debug("List of headers to retain was empty. fasta_dict cleared.\n")
            self.fasta_dict.clear()
            self.synchronize_seqs_n_headers()
            return

        pruned_fasta_dict = dict()
        unmapped = list()
        for seq_name in header_subset:
            try:
                pruned_fasta_dict[seq_name] = self.fasta_dict[seq_name]
            except KeyError:
                unmapped.append(seq_name)
        if len(pruned_fasta_dict) == 0:
            self.mapping_error(header_subset)

        if superset:
            unmapped.clear()

        if unmapped:
            logging.warning(str(len(unmapped)) + " sequences were not mapped to FASTA dictionary.\n")
            logging.debug("Headers that were not mapped to FASTA dictionary:\n\t" +
                          "\n\t".join(unmapped) + "\n")

        self.fasta_dict = pruned_fasta_dict
        self.synchronize_seqs_n_headers()

        return

    def remove_shorter_than(self, min_len: int):
        long_seqs = list()
        dropped = 0
        for seq_name in self.fasta_dict:
            if len(self.fasta_dict[seq_name]) > min_len:
                long_seqs.append(seq_name)
            else:
                dropped += 1
        if dropped >= 1:
            logging.debug(str(dropped) + " sequences were found to be shorter than " + str(min_len) + " and removed.\n")
        self.keep_only(long_seqs)
        return

    def get_header_mapping_dict(self) -> dict:
        """
        Creates a dictionary mapping all forms of the header formats (original, first_split, num, and formatted)
        to the unique, numeric TreeSAPP IDs for rapid look-ups
        :return: Dictionary for look-ups
        """
        mapping_dict = dict()
        mapping_dict.update(self.original_header_map())
        mapping_dict.update(self.formatted_header_map())
        mapping_dict.update(self.first_split_header_map())
        return mapping_dict

    def change_dict_keys(self, index_replace="original"):
        # TODO: Include a value to track the fasta dict key-type (e.g. num, original)
        # TODO: Use utilities.rekey_dict
        repl_fasta_dict = dict()
        for acc in sorted(self.header_registry, key=int):
            header = self.header_registry[acc]  # type: Header
            # Get the new header
            if index_replace == "original":
                new_header = header.original
            elif index_replace == "first_split":
                new_header = header.first_split
            elif index_replace == "formatted":
                new_header = header.formatted
            elif index_replace == "num":
                new_header = acc
            elif index_replace == "accession":
                new_header = header.accession
            else:
                logging.error("Unknown replacement type.\n")
                sys.exit(3)
            # Find the old header to be replaced
            # NOTE: Using .first_split may be lossy from duplicates. This will _not_ propagate under current scheme.
            if header.original in self.fasta_dict:
                repl_fasta_dict[new_header] = self.fasta_dict[header.original]
            elif header.formatted in self.fasta_dict:
                repl_fasta_dict[new_header] = self.fasta_dict[header.formatted]
            elif acc in self.fasta_dict:
                repl_fasta_dict[new_header] = self.fasta_dict[acc]
            elif header.first_split in self.fasta_dict:
                repl_fasta_dict[new_header] = self.fasta_dict[header.first_split]
            elif header.accession in self.fasta_dict:
                repl_fasta_dict[new_header] = self.fasta_dict[header.accession]
            else:
                pass

        if not repl_fasta_dict:
            logging.error("Unable to change dictionary keys as no headers in '" + self.file + "' were found in dict.\n")
            sys.exit(3)
        self.fasta_dict = repl_fasta_dict
        self.index_form = index_replace
        return

    def synchronize_seqs_n_headers(self):
        """
        If the header_registry and fasta_dict objects are of different sizes,
        the header registry is remade, excluding sequences that are not found in fasta_dict.
        The num_id is static during the synchronization so sequences can be mapped from before-and-after.
        :return:
        """
        excluded_headers = list()
        self.change_dict_keys("num")
        if len(self.fasta_dict.keys()) != len(self.header_registry):
            sync_header_registry = dict()
            sync_fasta_dict = dict()
            header_num_set = set(self.header_registry.keys())
            fasta_num_set = set(self.fasta_dict.keys())
            for num_id in header_num_set.intersection(fasta_num_set):
                sync_header_registry[num_id] = self.header_registry[num_id]
                sync_fasta_dict[num_id] = self.fasta_dict[num_id]
            for num_id in header_num_set.difference(fasta_num_set):
                try:
                    excluded_headers.append(self.header_registry[num_id].original)
                except KeyError:
                    logging.error("Unable to find TreeSAPP ID '%s' in header_registry.\n" % num_id)
                    sys.exit()
            self.header_registry = sync_header_registry
            self.fasta_dict = sync_fasta_dict
        if len(excluded_headers) >= 1:
            logging.debug("The following sequences were excluded after synchronizing FASTA:\n\t" +
                          "\n\t".join(excluded_headers) + "\n")
        if len(self.header_registry) == 0:
            logging.error("All sequences were discarded during header_registry and fasta_dict synchronization.\n")
            sys.exit()
        self.change_dict_keys()
        return

    def swap_headers(self, header_map):
        self.fasta_dict = rekey_dict(self.fasta_dict, header_map)
        self.header_registry = register_headers(list(self.fasta_dict.keys()), True)
        return

    def custom_lineage_headers(self, header_lineage_map: dict):
        """
        Converts a header to the TreeSAPP custom header format (below) using a dictionary

        Custom fasta header with taxonomy:
         First group = contig/sequence name, second = full taxonomic lineage, third = description for tree
         There are no character restrictions on the first and third groups
         The lineage must be formatted like:
          cellular organisms; Bacteria; Proteobacteria; Gammaproteobacteria
        :param header_lineage_map: A dictionary mapping sequence names/headers to NCBI-formatted taxonomic lineages
        :return: None
        """
        swap_header_map = dict()
        skipped = []
        header_re = re.compile(r"^>?([A-Za-z0-9.\-_]+\.?[0-9]?)\s+(\[.*\])$")
        for seq_name in header_lineage_map:
            try:
                accession, organism = header_re.match(seq_name).groups()
                swap_header_map[seq_name] = accession + " lineage=" + header_lineage_map[seq_name] + ' ' + organism
            except TypeError:
                skipped.append(seq_name)
                swap_header_map[seq_name] = seq_name
        self.swap_headers(swap_header_map)
        if skipped:
            logging.warning("Skipped converting sequence names to custom lineage format for the following:\n" +
                            ", ".join(skipped) + "\n")
        return

    def summarize_fasta_sequences(self):
        num_headers = 0
        longest = 0
        shortest = 0
        sequence_lengths = []
        for name in self.fasta_dict:
            sequence = self.fasta_dict[name]
            # Calculate stats
            num_headers += 1
            if longest == 0 and shortest == 0:
                longest = len(sequence)
                shortest = len(sequence)
            elif len(sequence) > longest:
                longest = len(sequence)
            elif len(sequence) < shortest:
                shortest = len(sequence)
            else:
                pass
            sequence_lengths.append(len(sequence))

        stats_string = "\tNumber of sequences: " + str(num_headers) + "\n"
        stats_string += "\tLongest sequence length: " + str(longest) + "\n"
        stats_string += "\tShortest sequence length: " + str(shortest) + "\n"
        stats_string += "\tMean sequence length: " + str(round(sum(sequence_lengths) / num_headers, 1)) + "\n"
        stats_string += "\tMedian sequence length: " + str(median(sequence_lengths)) + "\n"

        return stats_string

    def dedup_by_sequences(self):
        """
        Removes exact duplicate sequences from a FASTA-formatted dictionary of sequence records
        """
        dedup_dict = dict()
        duplicates = list()
        logging.debug("Checking for redundant FASTA records with duplicate sequences... ")
        if len(set(self.fasta_dict.values())) == len(self.fasta_dict):
            return
        else:
            for header in self.fasta_dict.keys():
                if self.fasta_dict[header] not in dedup_dict.values():
                    dedup_dict[header] = self.fasta_dict[header]
                else:
                    duplicates.append(header)
            self.fasta_dict = dedup_dict
        logging.debug("done.\n")

        if duplicates:
            logging.debug("Removed " + str(len(duplicates)) + " sequences with duplicate sequences.\n")
        self.synchronize_seqs_n_headers()

        return

    def dedup_by_accession(self):
        count_dict = dict()
        dedup_header_dict = dict()
        duplicates = list()
        logging.debug("Checking for redundant sequences with duplicate accessions.\n")
        for acc in self.header_registry:
            header = self.header_registry[acc]  # type: Header
            if not header.accession:
                continue
            try:
                count_dict[header.accession].append(acc)
            except KeyError:
                count_dict[header.accession] = [acc]
        for accession in count_dict:
            if len(count_dict[accession]) > 1:
                while len(count_dict[accession]) > 1:
                    # Discard the duplicates
                    duplicates.append(count_dict[accession].pop())
            rep_acc = count_dict[accession].pop()
            dedup_header_dict[rep_acc] = self.header_registry[rep_acc]
        if duplicates:
            logging.debug("Removed " + str(len(duplicates)) + " sequences with duplicate accessions:\n\t" +
                          "\n\t".join(duplicates) + "\n")
        self.header_registry = dedup_header_dict
        self.synchronize_seqs_n_headers()
        return

    def update(self, fasta, file=True) -> None:
        """
        Used to append sequences to the current FASTA object. All new sequence records are assigned a new numerical ID,
        continuing from the largest ID found in the original FASTA instance.

        :param fasta: Either a fasta-formatted dictionary (keys are sequence names (i.e. headers) and values are seqs)
        or the name of a file, which will be read to generate a fasta-formatted dictionary
        :param file: Flag indicating whether the update is using a file or a fasta-formatted dictionary
        :return: None
        """
        # Format the inputs
        if file:
            if not os.path.isfile(fasta):
                logging.error("File '" + fasta + "' does not exist!\n")
                sys.exit(13)
            new_fasta = FASTA(fasta)
            new_fasta.load_fasta()
        else:
            new_fasta = FASTA("dummy_name")
            new_fasta.fasta_dict = fasta
            new_fasta.header_registry = register_headers(fasta.keys())

        # Guarantee the index type is original for both self.fasta_dict and new_fasta
        self.change_dict_keys()
        header_map = self.get_header_mapping_dict()
        new_fasta.change_dict_keys()

        # TODO: Potentially refactor the following code into a function called 'join' used to merge two FASTA objects
        # Load the new fasta and headers
        acc = max([int(x) for x in self.header_registry.keys()]) + 1
        for num_id in sorted(new_fasta.header_registry, key=int):
            header = new_fasta.header_registry[num_id]  # type: Header
            if header.original in header_map:
                for h_i in header_map[header.original]:  # type: Header
                    self.amendments.add(str(h_i.treesapp_num_id))
            elif header.first_split in header_map:
                for h_i in header_map[header.first_split]:  # type: Header
                    self.amendments.add(str(h_i.treesapp_num_id))
            else:
                self.fasta_dict[header.original] = new_fasta.fasta_dict[header.original]
                ts_id = acc
                self.header_registry[str(ts_id)] = new_fasta.header_registry[num_id]
                self.header_registry[str(ts_id)].treesapp_num_id = ts_id
                self.amendments.add(str(ts_id))
                acc += 1
        self.synchronize_seqs_n_headers()
        return

    def unalign(self):
        unaligned_dict = {}
        for header in self.fasta_dict:
            seq = self.fasta_dict[header]
            unaligned_dict[header] = re.sub("[-.]", '', seq)
        self.fasta_dict = unaligned_dict
        return


def merge_fasta_dicts_by_index(extracted_seq_dict, numeric_contig_index):
    merged_extracted_seq_dict = dict()
    for marker in extracted_seq_dict:
        for bin_num in extracted_seq_dict[marker]:
            for i in extracted_seq_dict[marker][bin_num]:
                merged_extracted_seq_dict[numeric_contig_index[marker][i]] = extracted_seq_dict[marker][bin_num][i]
    return merged_extracted_seq_dict


def write_classified_sequences(tree_saps: dict, formatted_fasta_dict: dict, fasta_file: str):
    """
    Function to write the nucleotide sequences representing the full-length ORF for each classified sequence
    Sequence names are from ItolJplace.contig_name values so output format is:

     >contig_name|RefPkg|StartCoord_StopCoord
    :param tree_saps: A dictionary of gene_codes as keys and TreeSap objects as values
    :param formatted_fasta_dict: A dictionary with headers/sequence names as keys and sequences as values
    :param fasta_file: Path to a file to write the sequences to in FASTA format
    :return: None
    """
    output_fasta_dict = dict()
    prefix = ''  # For adding a '>' if the formatted_fasta_dict sequences have them
    for seq_name in formatted_fasta_dict:
        if seq_name[0] == '>':
            prefix = '>'
        break

    for denominator in tree_saps:
        for placed_sequence in tree_saps[denominator]:  # type ItolJplace
            if placed_sequence.classified:
                output_fasta_dict[placed_sequence.contig_name] = ""
                try:
                    output_fasta_dict[placed_sequence.contig_name] = formatted_fasta_dict[prefix +
                                                                                          placed_sequence.contig_name]
                except KeyError:
                    seq_name = re.sub(r"\|{0}\|\d+_\d+.*".format(placed_sequence.name), '', placed_sequence.contig_name)
                    try:
                        output_fasta_dict[placed_sequence.contig_name] = formatted_fasta_dict[prefix + seq_name]
                    except KeyError:
                        logging.error("Unable to find '" + prefix + placed_sequence.contig_name +
                                      "' in predicted ORFs file!\nExample headers in the predicted ORFs file:\n\t" +
                                      '\n\t'.join(list(formatted_fasta_dict.keys())[:6]) + "\n")
                        sys.exit(3)

                if not placed_sequence.seq_len:
                    placed_sequence.seq_len = len(output_fasta_dict[placed_sequence.contig_name])

    if output_fasta_dict:
        write_new_fasta(output_fasta_dict, fasta_file)

    return


def format_read_fasta(fasta_input: str, molecule: str, subset=None, min_seq_length=10):
    """
    Reads a FASTA file, ensuring each sequence and sequence name is valid.

    :param fasta_input: Absolute path of the FASTA file to be read
    :param molecule: Molecule type of the sequences ['prot', 'dna', 'rrna']
    :param subset: A set for filtering sequences. Only sequences with names in subset will be included in the dictionary
    :type subset: set
    :param min_seq_length: All sequences shorter than this will not be included in the returned list.
    :return: A Python dictionary with headers as keys and sequences as values
    """
    start = time()

    if molecule == "prot":
        bad_chars = re.compile(r"[OU]")
    else:
        bad_chars = re.compile(r"[EFIJLOPQZ]")
    bad_seqs = set()

    if subset and type(subset) is not set:
        logging.error("Unexpected type of subset: {} instead of set.\n".format(type(subset)))
        sys.exit(13)

    formatted_fasta_dict = {}
    for name, seq in pyfastx.Fasta(fasta_input, build_index=False):  # type: (str, str)
        if len(seq) < min_seq_length:
            continue
        if subset:
            if name in subset:
                formatted_fasta_dict[name] = seq
        else:
            if bad_chars.search(seq):
                bad_seqs.add(name)
                continue
            formatted_fasta_dict[name] = seq
    end = time()
    logging.debug("{} read by pyfastx in {} seconds.\n".format(fasta_input, end-start))

    if len(formatted_fasta_dict) == 0:
        logging.error("No sequences in FASTA {0} were saved.\n"
                      "Either the molecule type specified ({1}) or minimum sequence length ({2}) may be unsuitable.\n"
                      "Consider changing these before rerunning.\n".format(fasta_input, molecule, min_seq_length))
        sys.exit(13)

    if len(bad_seqs) > 0:
        logging.debug("The following sequences were removed due to bad characters:\n" +
                      "\n".join(bad_seqs) + "\n")

    return formatted_fasta_dict


def get_headers(fasta_file):
    """
    Reads a FASTA file and returns a list of all headers it found in the file. No reformatting or filtering performed.
    :param fasta_file: Path to the FASTA file to be read.
    :return:
    """
    original_headers = list()
    try:
        fasta = open(fasta_file, 'r')
    except IOError:
        logging.error("Unable to open the FASTA file '" + fasta_file + "' for reading!\n")
        sys.exit(5)

    n_headers = 0
    for name, _ in generate_fasta(fasta):
        n_headers += 1
        original_headers.append('>' + str(name))

    fasta.close()
    if len(original_headers) == 0:
        logging.error("No sequence headers read from FASTA file " + fasta_file + "\n")
        sys.exit(3)
    logging.debug("Read " + str(n_headers) + " headers from " + fasta_file + ".\n")

    return original_headers


def write_new_fasta(fasta_dict, fasta_name, max_seqs=None, headers=None):
    """
    Function for writing sequences stored in dictionary to file in FASTA format; optional filtering with headers list

    :param fasta_dict: A dictionary containing headers as keys and sequences as values
    :param fasta_name: Name of the FASTA file to write to
    :param max_seqs: If not None, the maximum number of sequences to write to a single FASTA file
    :param headers: Optional list of sequence headers. Only fasta_dict keys in headers will be written
    :return:
    """
    split_files = list()
    file_counter = 1
    sequence_accumulator = 0
    n_char = 0
    fasta_string = ""

    # Check for '>' leading sequence names. Strip them if present.
    if headers:
        side_chevy = headers[0][0] == '>'
        for header in headers:
            state = header[0] == '>'
            if state is not side_chevy:
                logging.error("Inconsistent header names in headers list object\n")
                sys.exit(5)
        if side_chevy:
            headers = [header[1:] for header in headers]

    if max_seqs is not None:
        fasta_name = fasta_name + '_' + str(file_counter) + ".fasta"

    try:
        fa_out = open(fasta_name, 'w')
    except IOError:
        logging.error("Unable to open " + fasta_name + " for writing!\n")
        sys.exit(5)

    for name in sorted(fasta_dict.keys()):
        seq = fasta_dict[name]
        if name[0] != '>':
            name = '>' + name
        sequence_accumulator += 1
        n_char += len(seq) + len(name)
        if max_seqs and sequence_accumulator > max_seqs:
            # If input is to be split and number of sequences per file has been exceeded begin writing to new file
            fa_out.write(fasta_string)
            fa_out.close()
            fasta_string = ""
            split_files.append(fasta_name)
            file_counter += 1
            sequence_accumulator = 1
            fasta_name = re.sub(r"_\d+.fasta$", '_' + str(file_counter) + ".fasta", fasta_name)
            fa_out = open(fasta_name, 'w')

        # Only write the records specified in `headers`, or all records if `headers` isn't provided
        if headers is None:
            fasta_string += name + "\n" + seq + "\n"
        elif name[1:] in headers:
            fasta_string += name + "\n" + seq + "\n"
        else:
            sequence_accumulator -= 1
        # Write the string buffer if it hold more than 100 million characters
        if n_char > 1E7:
            fa_out.write(fasta_string)
            fasta_string = ""
            n_char = 0

    fa_out.write(fasta_string)
    fa_out.close()
    split_files.append(fasta_name)

    return split_files


def get_header_format(header, code_name=""):
    """
    Used to decipher which formatting style was used and parse information, ideally reliably
    HOW TO ADD A NEW REGULAR EXPRESSION:
        1. create a new compiled regex pattern, like below
        2. add the name of the compiled regex pattern to the header_regexes dictionary
        3. if the regex groups are new and complicated (parsing more than the accession and organism info),
        alter return_sequence_info_groups in create_treesapp_ref_data to add another case

    :param header: A sequences header from a FASTA file
    :param code_name:
    :return:
    """
    # The regular expressions with the accession and organism name grouped
    # Protein databases:
    gi_re = re.compile(r">?gi\|(\d+)\|[a-z]+\|[_A-Z0-9.]+\|.* RecName: Full=([A-Za-z1-9 _\-]+);?.*$")  # a
    gi_prepend_proper_re = re.compile(r">?gi\|\d+\|[a-z]{2,4}\|([_A-Z0-9.]+)\| (.*) \[(.*)\]$")  # a, d, o
    gi_prepend_mess_re = re.compile(r">?gi\|(\d+)\|[a-z]{2,4}\|.*\|([\w\s.,\-()]+)$")  # a
    dbj_re = re.compile(r">?dbj\|(.*)\|.*\[(.*)\]")  # a, o
    emb_re = re.compile(r">?emb\|(.*)\|.*\[(.*)\]")
    gb_re = re.compile(r">?gb\|(.*)\|.*\[(.*)\]")
    ref_re = re.compile(r">?ref\|(.*)\|.*\[(.*)\]")
    pdb_re = re.compile(r">?pdb\|(.*)\|.+$")  # a
    pir_re = re.compile(r">?pir\|.*\|(\w+).* - (.*)$")  # a, o
    presf_re = re.compile(r">?prf\|.*\|([A-Z0-9]+)\s+.*$")  # a
    sp_re = re.compile(r">?sp\|(.*)\|.*$")  # a
    tr_re = re.compile(r">?tr\|(\w+)\|\w+_\w+ .* OS=(.*) GN=.*$")  # a, o
    fungene_gi_bad = re.compile(r"^>?[0-9]+\s+coded_by=.+,organism=.+,definition=.+$")
    treesapp_re = re.compile(r"^>?(\d+)_" + re.escape(code_name) + "$")
    pfam_re = re.compile(r"^>?([A-Za-z0-9_|]+)/[0-9]+-[0-9]+$")  # a
    eggnog_re = re.compile(r"^>?(\d+)\.([A-Za-z][-A-Za-z0-9_]+)(\.\d)?(\s\[.*\])?$")  # t, o
    eggnot_re = re.compile(r">?eggnog\|(\d+)\|(.*)")  # a
    # Nucleotide databases:
    # silva_arb_re = re.compile("^>([A-Z0-9]+)\.([0-9]+)\.([0-9]+)_(.*)$")
    # refseq_nuc_re = re.compile("^>([A-Z]+_[0-9]+\.[0-9])_.+$")  # a
    # nr_re = re.compile("^>([A-Z0-9]+\.[0-9])_.*$")  # a

    # Ambiguous:
    # genbank_exact_genome = re.compile("^>([A-Z]{1,2}[0-9]{5,6}\.?[0-9]?) .* \[(.*)\]$")  # a, o
    accession_only = re.compile(r"^>?([A-Z]+_?[0-9]+\.?[0-9]?)$")  # a
    ncbi_ambiguous = re.compile(r"^>?([A-Za-z0-9.\-_]+)\s+.*(?<!])$")  # a
    ncbi_org = re.compile(r"^>?([A-Z][A-Za-z0-9.\-_]+\.?[0-9]?)\s+(?!lineage=).*\[.*\]$")  # a
    assign_re = re.compile(r"^>?(.*)\|({0})\|(\d+_\d+)$".format(re.escape(code_name)))  # a, d, l

    # Custom fasta header with taxonomy:
    # First group = contig/sequence name, second = full taxonomic lineage, third = description for tree
    # There are no character restrictions on the first and third groups
    # The lineage must be formatted like:
    #   cellular organisms; Bacteria; Proteobacteria; Gammaproteobacteria
    custom_tax = re.compile(r"^>?(.*) lineage=([A-Za-z ]+.*) \[(.*)\]$")  # a, l, o

    header_regexes = {"prot": {dbj_re: "dbj",
                               emb_re: "emb",
                               gb_re: "gb",
                               pdb_re: "pdb",
                               pir_re: "pir",
                               ref_re: "ref",
                               sp_re: "sp",
                               tr_re: "tr",
                               gi_re: "gi_re",
                               gi_prepend_proper_re: "gi_proper",
                               gi_prepend_mess_re: "gi_mess",
                               pfam_re: "pfam",
                               presf_re: "prf",
                               eggnog_re: "eggnog",
                               eggnot_re: "eggnot"},
                      "dna": {treesapp_re: "treesapp"},
                      "ambig": {accession_only: "bare",
                                ncbi_ambiguous: "ncbi_ambig",
                                ncbi_org: "ncbi_org",
                                custom_tax: "custom",
                                assign_re: "ts_assign"}
                      }

    if fungene_gi_bad.match(header):
        logging.warning(header + " uses GI numbers which are now unsupported by the NCBI! " +
                        "Consider switching to Accession.Version identifiers instead.\n")

    header_format_re = None
    header_db = None
    header_molecule = None
    format_matches = list()
    for molecule in header_regexes:
        for regex in header_regexes[molecule]:
            if regex.match(header):
                header_format_re = regex
                header_db = header_regexes[molecule][regex]
                header_molecule = molecule
                format_matches.append(header_db)
            else:
                pass
    if len(format_matches) > 1:
        if "ts_assign" in format_matches:
            return assign_re, "ts_assign", "ambig"
        logging.error("Header '" + header + "' matches multiple potential formats:\n\t" +
                      ", ".join(format_matches) + "\n" +
                      "TreeSAPP is unable to parse necessary information properly.\n")
        sys.exit(5)

    if header_format_re is None:
        logging.error("Unable to parse header '" + header + "'\n")
        sys.exit(5)

    return header_format_re, header_db, header_molecule


def summarize_fasta_sequences(fasta_file):
    try:
        fasta_handler = open(fasta_file, 'r')
    except IOError:
        logging.error("Unable to open " + fasta_file + " for reading!\n")
        sys.exit(5)

    num_headers = 0
    longest = 0
    shortest = 0
    sequence = None
    sequence_lengths = []
    line = fasta_handler.readline()
    while line:
        if line[0] == '>':
            num_headers += 1
            if sequence is not None:
                if longest == 0 and shortest == 0:
                    longest = len(sequence)
                    shortest = len(sequence)
                if len(sequence) > longest:
                    longest = len(sequence)
                if len(sequence) < shortest:
                    shortest = len(sequence)
                else:
                    pass
                sequence_lengths.append(len(sequence))
            sequence = ""
        else:
            sequence += line.strip()
        line = fasta_handler.readline()
    # Log the last sequence in the file
    if longest == 0 and shortest == 0:
        longest = len(sequence)
        shortest = len(sequence)
    if len(sequence) > longest:
        longest = len(sequence)
    if len(sequence) < shortest:
        shortest = len(sequence)
    sequence_lengths.append(len(sequence))

    stats_string = "\tNumber of sequences: " + str(num_headers) + "\n"
    stats_string += "\tLongest sequence length: " + str(longest) + "\n"
    stats_string += "\tShortest sequence length: " + str(shortest) + "\n"
    stats_string += "\tMean sequence length: " + str(round(sum(sequence_lengths)/num_headers, 1)) + "\n"
    stats_string += "\tMedian sequence length: " + str(median(sequence_lengths)) + "\n"

    logging.info(stats_string)
    return


def rename_cluster_headers(cluster_dict, header_registry):
    """
    Map the numerical TreeSAPP IDs to each sequence's original header
    cluster.representative and header are both 'treesapp_id's
    :param cluster_dict:
    :param header_registry:
    :return:
    """
    for num_id in cluster_dict:
        members = list()
        cluster = cluster_dict[num_id]
        try:
            cluster.representative = header_registry[cluster.representative].original
        except KeyError:
            logging.error("Unable to find '" + cluster.representative + "' in formatted header-registry names.\n")
            sys.exit(7)
        for member in cluster.members:
            header, identity = member
            members.append([header_registry[header].original, identity])
        cluster.members = members
    return
