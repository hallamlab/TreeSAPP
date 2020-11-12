__author__ = 'Connor Morgan-Lang'

import sys
import re
import os
import logging
from time import sleep, time
from math import ceil

from pyfastx import Fasta, Fastq
from pyfastxcli import fastx_format_check
from collections import namedtuple

from treesapp.utilities import median, reformat_string, rekey_dict


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


def split_seq_writer_helper(fasta_string: str, fh, max_file_size: int):
    fh.write(fasta_string)

    # Reduce the number of checks
    fh.seek(0, 2)
    if fh.tell() > max_file_size:
        fh.close()
        fh = None
    return fh


def spawn_new_file(file_num, outdir, file_name, ext="fasta"):
    """
    Creates a new file name based on the file_num and prefix.

    :param file_num: Previous file number component of the file name in the series. Will be incremented.
    :param outdir: Path to the directory to write files to
    :param file_name: Prefix name of the new file
    :param ext: File extension to use for the new file
    :return: tuple(file_handler, name of the new file, number of file in the series)
    """
    file_num += 1

    subfile = "{}.{}.{}".format(file_name, str(file_num), ext)
    if outdir is not None:
        subfile = os.path.join(outdir, subfile)

    fh = open(subfile, 'w')
    logging.debug("Writing split {} to file '{}'.\n".format(ext, subfile))

    return fh, file_num


def parameterize_sub_file_size(file_name: str, num_files: int, max_seqs: int):
    # Determine the size of the FASTA file instead of building an index
    f_size = os.path.getsize(file_name)

    if num_files > 1:
        max_file_size = ceil(f_size / num_files)
        max_seq_count = 0
    elif max_seqs > 1:
        max_file_size = 0
        max_seq_count = max_seqs
    else:
        logging.error("Unable to split FASTA into {} files and {} sequences.\n".format(num_files, max_seqs))
        sys.exit(3)

    return max_file_size, max_seq_count


def split_fa(fastx: str, outdir: str, file_num=1, max_seq_count=0):
    """
    Credit: Lianming Du of pyfastx (v0.6.10)

    Differences are:
1. Format of the output files is always uncompressed FASTA, compatible with ORF prediction tools e.g. Prodigal
2. Batch writing of reads instead of writing a single read at a time to reduce number of  I/O operations

    :param fastx: Path to a FASTQ file to be read and converted to FASTA
    :param outdir: Path to write the output files
    :param file_num: Number of files to split the input FASTQ file into
    :param max_seq_count: Instead of splitting into N number of files FASTQ can be split by number of sequences
    :return: List of fasta files generated
    """
    start = time()
    # Fix path with tilde
    fastx = os.path.expanduser(fastx)
    outdir = os.path.expanduser(outdir)
    outputs = []
    fa = Fasta(file_name=fastx, build_index=False, full_name=True)

    file_name, suffix1 = os.path.splitext(os.path.basename(fastx))
    if fa.is_gzip:
        file_name, suffix2 = os.path.splitext(file_name)

    max_file_size, max_seq_count = parameterize_sub_file_size(fastx, file_num, max_seq_count)

    seq_write = 0
    max_buffer_size = min(1E5, max_file_size)  # Controls how often the strings are written to the file
    fasta_string = ""  # Stores the intermediate fasta strings
    fh = None  # File handler
    file_num = 0  # The current file number

    for name, seq in fa:  # type: (str, str)
        fasta_string += ">%s\n%s\n" % (name, seq)
        seq_write += 1
        # Batch write
        if len(fasta_string) >= max_buffer_size and seq_write >= max_seq_count:
            if not fh:
                fh, file_num = spawn_new_file(file_num, outdir, file_name)
                outputs.append(fh.name)
            fh = split_seq_writer_helper(fasta_string, fh, max_file_size)
            fasta_string = ""
            seq_write = 0

    if fasta_string:
        try:
            fh.write(fasta_string)
            fh.close()
        except AttributeError:
            fh, file_num = spawn_new_file(file_num, outdir, file_name)
            outputs.append(fh.name)
            fh.write(fasta_string)
            fh.close()
    end = time()

    logging.debug("{} completed split_fa in {}s.\n".format(fastx, end - start))
    return outputs


def fq2fa(fastx: str, outdir: str, file_num=1, max_seq_count=0) -> list:
    """
    Credit: Lianming Du of pyfastx (v0.6.10)
    A modified version of the function fastq_split in the pyfastx python package: https://github.com/lmdu/pyfastx

    Differences are:
1. Format of the output files is always uncompressed FASTA, compatible with ORF prediction tools e.g. Prodigal
2. Batch writing of reads instead of writing a single read at a time to reduce number of  I/O operations

    :param fastx: Path to a FASTQ file to be read and converted to FASTA
    :param outdir: Path to write the output files
    :param file_num: Number of files to split the input FASTQ file into
    :param max_seq_count: Instead of splitting into N number of files FASTQ can be split by number of sequences
    :return: List of fasta files generated
    """
    start = time()
    # Fix path with tilde
    fastx = os.path.expanduser(fastx)
    outdir = os.path.expanduser(outdir)
    outputs = []
    fq = Fastq(file_name=fastx, build_index=False)

    file_name, suffix1 = os.path.splitext(os.path.basename(fastx))
    if fq.is_gzip:
        file_name, suffix2 = os.path.splitext(file_name)

    max_file_size, max_seq_count = parameterize_sub_file_size(fastx, file_num, max_seq_count)

    seq_write = 0
    max_buffer_size = min(1E5, max_file_size)  # Controls how often the strings are written to the file
    fasta_string = ""  # Stores the intermediate fasta strings
    fh = None  # File handler
    file_num = 0  # The current file number

    for read_name, seq, _ in fq:  # type: (str, str)
        fasta_string += ">%s\n%s\n" % (read_name, seq)
        seq_write += 1
        # Batch write
        if len(fasta_string) >= max_buffer_size and seq_write >= max_seq_count:
            if not fh:
                fh, file_num = spawn_new_file(file_num, outdir, file_name)
                outputs.append(fh.name)
            fh = split_seq_writer_helper(fasta_string, fh, max_file_size)
            fasta_string = ""
            seq_write = 0

    if fasta_string:
        try:
            fh.write(fasta_string)
            fh.close()
        except AttributeError:
            fh, file_num = spawn_new_file(file_num, outdir, file_name)
            outputs.append(fh.name)
            fh.write(fasta_string)
            fh.close()
    end = time()

    logging.debug("{} completed fq2fa in {}s.\n".format(fastx, end-start))

    return outputs


def read_fasta_to_dict(fasta_file: str) -> dict:
    """
    Reads any fasta file using the pyfastx library

    :param fasta_file: Path to a FASTA file to be read into a dict
    :return: Dict where headers/record names are keys and sequences are the values
    """
    fasta_dict = dict()

    if not os.path.exists(fasta_file):
        logging.debug("'{}' fasta file doesn't exist.\n".format(fasta_file))
        return fasta_dict

    try:
        py_fa = Fasta(fasta_file, build_index=False, full_name=True)
    except RuntimeError as error:
        logging.debug(str(error)+"\n")
        return fasta_dict

    for name, seq in py_fa:  # type: (str, str)
        fasta_dict[name] = seq.upper()

    return fasta_dict


class Header:
    def __init__(self, header):
        self.original = header
        self.treesapp_num_id = 0
        self.formatted = ""
        self.post_align = ""
        self.first_split = ""
        self.accession = ""
        self.version = ""

    def get_info(self):
        info_string = "TreeSAPP ID = '%s'\tPrefix = '%s'\n" % (str(self.treesapp_num_id), self.first_split)
        info_string += "Original =  %s\nFormatted = %s\n" % (self.original, self.formatted)
        if self.accession:
            info_string += "Accession = " + self.accession + "\n"
        return info_string

    def find_accession(self, refpkg_name="") -> None:
        header_regexes = load_fasta_header_regexes(refpkg_name)
        header_format_re, header_db, header_molecule = get_header_format(self.original, header_regexes)
        sequence_info = header_format_re.match(self.original)
        self.accession = sequence_info_groups(sequence_info, header_db, self.original, header_regexes).accession
        return


def register_headers(header_list: list, drop=True) -> dict:
    """
    Instantiates a Header instance for each sequence name (i.e. header) in header_list.
    The attributes 'formatted', 'first_split', 'original' and 'treesapp_num_id' are populated

    :param header_list: A list of headers to parse
    :param drop: A flag indicating whether the '>' character should be dropped from the sequence names
    :return: A dictionary of Header instances indexed by a numerical identifier
    """
    acc = 1
    header_registry = dict()
    dup_checker = set()
    dups = []
    for header in header_list:
        header = header.strip()
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
        self.header_registry = dict()  # A dictionary of Header instances indexed by a unique numerical identifier
        self.amendments = set()  # Set of the TreeSAPP numerical identifiers for all guaranteed sequences
        self.index_form = None

    def clone(self, fasta) -> None:
        """
        Used for cloning a one FASTA instance into the current one.

        :param fasta: A FASTA instance with attributes to copy into the current instance
        :return: None
        """
        self.file = fasta.file
        self.fasta_dict.update(fasta.fasta_dict)
        self.header_registry.update(fasta.header_registry)
        self.amendments.update(fasta.amendments)
        self.index_form = fasta.index_form
        return

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
        logging.error("No classified sequences were mapped to '{}' FASTA dictionary.\n"
                      "Here are some example names from the mapping list:\n\t{}\n"
                      "And example names from FASTA dict:\n\t{}\n".format(self.file,
                                                                          "\n\t".join(sorted(bad_headers)[0:6]),
                                                                          "\n\t".join(list(sorted(
                                                                              self.fasta_dict.keys()))[0:6])))
        raise AssertionError

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

    def get_accession_header_map(self) -> dict:
        accession_header_map = dict()
        for index, header in self.header_registry.items():  # type: (str, Header)
            if len(header.accession) == 0:
                header.find_accession()
            try:
                accession_header_map[header.accession].append(header)
            except KeyError:
                accession_header_map[header.accession] = [header]
            except TypeError:
                if not header.accession:
                    logging.error("Attempting to create an accession:header dictionary but"
                                  " accession could not be set for header '{}'.\n".format(header.original))
                    sys.exit(17)
                else:
                    raise TypeError
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

    def create_header_mapping_table(self):
        table_str = "\t".join(["TreeSAPP number", "Accession", "Original", "Formatted"]) + "\n"
        for acc in sorted(self.header_registry):
            header = self.header_registry[acc]  # type: Header
            table_str += "\t".join([header.treesapp_num_id, header.accession, header.original, header.formatted]) + "\n"
        return table_str

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
            logging.warning("{} sequences were not mapped to FASTA dictionary.\n".format(str(len(unmapped))))
            logging.debug("Headers that were not mapped to FASTA dictionary:\n\t{}\n".format("\n\t".join(unmapped)))

        self.fasta_dict = pruned_fasta_dict
        self.synchronize_seqs_n_headers()

        return

    def replace_ambiguity_chars(self, molecule, replace_char='X'):
        if molecule == "prot":
            invalid = {'U', 'O'}
        else:
            logging.debug("FASTA.replace_ambiguity_chars is not equipped to handle molecule type '%s'.\n" % molecule)
            return
        invalid_re = re.compile('|'.join(invalid))
        bad_seqs = 0
        for seq_name in self.fasta_dict:
            if invalid_re.search(self.fasta_dict[seq_name]):
                bad_seqs += 1
            self.fasta_dict[seq_name] = invalid_re.sub(replace_char, self.fasta_dict[seq_name])

        logging.debug("Identified and replaced invalid ambiguity characters in %d sequences.\n" % bad_seqs)

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
            logging.debug("{} sequences were found to be shorter than {} and removed.\n".format(dropped, min_len))
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
        if len(self.header_registry) == 0:
            logging.error("FASTA.header_registry is empty. Unable to change dictionary keys.\n")
            raise AssertionError

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
            # NOTE: Using .first_split may be lossy if duplicates exist. This will _not_ propagate under current scheme.
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

    def synchronize_seqs_n_headers(self) -> None:
        """
        If the header_registry and fasta_dict objects are of different sizes,
        the header registry is remade, excluding sequences that are not found in fasta_dict.
        The num_id is static during the synchronization so sequences can be mapped from before-and-after.

        :return: None
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
            sys.exit(-1)
        if len(self.fasta_dict) == 0:
            logging.error("No fasta sequence names were mapped to the header registry!\n")
            sys.exit(-1)
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
            logging.debug("Removed {} sequences with duplicate sequences.\n".format(len(duplicates)))
        self.synchronize_seqs_n_headers()

        return

    def dedup_by_accession(self) -> None:
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
            logging.debug("Removed {} sequences with duplicate accessions:\n"
                          "\t{}\n".format(len(duplicates), "\n\t".join(duplicates)))
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

    def unalign(self) -> None:
        """
        Removes common multiple sequence alignments characters ('-', '.') from all sequences in self.fasta_dict

        :return: None
        """
        for header, seq in self.fasta_dict.items():
            if seq.find('-') >= 0:
                self.fasta_dict[header] = re.sub("[-.]", '', seq)
            else:
                self.fasta_dict[header] = seq
        return


def merge_fasta_dicts_by_index(extracted_seq_dict, numeric_contig_index):
    merged_extracted_seq_dict = dict()
    for marker in extracted_seq_dict:
        for bin_num in extracted_seq_dict[marker]:
            for i in extracted_seq_dict[marker][bin_num]:
                merged_extracted_seq_dict[numeric_contig_index[marker][i]] = extracted_seq_dict[marker][bin_num][i]
    return merged_extracted_seq_dict


def format_fasta(fasta_input: str, molecule: str, output_fasta: str, min_seq_length=10) -> dict:
    """
    Reads a FASTA file, ensuring each sequence and sequence name is valid, and writes the valid sequence to a new FASTA.
    Only headers are read into memory and the formatted FASTA is saved to a buffer before written to a file and cleared.

    :param fasta_input: Absolute path of the FASTA file to be read
    :param molecule: Molecule type of the sequences ['prot', 'dna', 'rrna']
    :param output_fasta: Path to the formatted FASTA file to write
    :param min_seq_length: All sequences shorter than this will not be included in the returned list.
    :return: A dictionary of Header instances indexed by a numerical identifier
    """
    start = time()

    # Select the alphabet to use when determining whether there are any bad characters
    if molecule == "prot":
        bad_chars = re.compile(r"[OUou\d]")
    else:
        bad_chars = re.compile(r"[EFIJLOPQZefijlopqz\d]")
    bad_seqs = set()

    # Open the output FASTA for writing
    try:
        fa_out_handle = open(output_fasta, 'w')
    except IOError:
        logging.error("Unable to open '{}' for writing.\n")
        sys.exit(15)

    headers = []
    max_buffer_size = 1E4
    seq_acc = 0
    fasta_string = ""
    for name, seq in Fasta(fasta_input, build_index=False, full_name=True):  # type: (str, str)
        if len(seq) < min_seq_length:
            continue
        if bad_chars.search(seq):
            bad_seqs.add(name)
            continue

        seq_acc += 1
        headers.append(name)
        fasta_string += ">{}\n{}\n".format(seq_acc, seq)

        # Write the fasta_string to the output fasta if the size exceeds the max_buffer_size
        if len(fasta_string) > max_buffer_size:
            fa_out_handle.write(fasta_string)
            fasta_string = ""

    # Write the final chunk in the FASTA file
    fa_out_handle.write(fasta_string)

    end = time()
    logging.debug("{} read by pyfastx in {} seconds.\n".format(fasta_input, round(end-start, 2)))

    if len(headers) == 0:
        logging.error("No sequences in FASTA {0} were saved.\n"
                      "Either the molecule type specified ({1}) or minimum sequence length ({2}) may be unsuitable.\n"
                      "Consider changing these before rerunning.\n".format(fasta_input, molecule, min_seq_length))
        sys.exit(13)

    if len(bad_seqs) > 0:
        logging.debug("The following sequences were removed due to bad characters:\n" +
                      "\n".join(bad_seqs) + "\n")

    header_registry = register_headers(headers, True)
    # if len(header_registry) != seq_acc:
    #     logging.error("The number of sequences read ({}) does not equal"
    #                   " the number of sequence names registered ({}).\n".format(seq_acc, len(header_registry)))
    #     sys.exit(13)

    return header_registry


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
        bad_chars = re.compile(r"[OUou\d]")
    else:
        bad_chars = re.compile(r"[EFIJLOPQZefijlopqz\d]")
    bad_seqs = set()

    if subset and type(subset) is not set:
        logging.error("Unexpected type of subset: {} instead of set.\n".format(type(subset)))
        sys.exit(13)

    formatted_fasta_dict = {}
    for name, seq in Fasta(fasta_input, build_index=False, full_name=True):  # type: (str, str)
        if len(seq) < min_seq_length:
            continue
        if subset:
            if name in subset:
                formatted_fasta_dict[name] = seq
        else:
            if bad_chars.search(seq):
                bad_seqs.add(name)
                continue
            formatted_fasta_dict[name.rstrip()] = seq
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


def get_headers(fasta_file: str) -> list:
    """
    Reads a FASTA file and returns a list of all headers it found in the file. No reformatting or filtering performed.

    :param fasta_file: Path to the FASTA file to be read.
    :return:
    """
    original_headers = list()
    if not os.path.exists(fasta_file):
        logging.error("'{}' fasta file doesn't exist.\n".format(fasta_file))

    n_headers = 0
    try:
        fa = Fasta(file_name=fasta_file, build_index=False, full_name=True)
    except RuntimeError:
        logging.warning("Pyfastx is unable to open '{}' to read headers. "
                        "There is a chance this is an empty file and will be skipped.\n".format(fasta_file))
        return original_headers

    for name, seq in fa:  # type: (str, str)
        n_headers += 1
        original_headers.append('>' + str(name))

    if len(original_headers) == 0:
        # Not a good idea to exit right from here, handle it case-by-case
        logging.warning("No sequence headers read from FASTA file '{}'\n".format(fasta_file))
    else:
        logging.debug("Read {} headers from FASTA file '{}'.\n".format(n_headers, fasta_file))

    return original_headers


def write_new_fasta(fasta_dict: dict, fasta_name: str, max_seqs=None, headers=None) -> list:
    """
    Function for writing sequences stored in dictionary to file in FASTA format; optional filtering with headers list

    :param fasta_dict: A dictionary containing headers as keys and sequences as values
    :param fasta_name: Name of the FASTA file to write to
    :param max_seqs: If not None, the maximum number of sequences to write to a single FASTA file
    :param headers: Optional list of sequence headers. Only fasta_dict keys in headers will be written
    :return: List of FASTA files written to
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


def load_fasta_header_regexes(code_name="") -> dict:
    """
    Create the dictionary of all currently known regular expressions that match database-specific fasta headers
    HOW TO ADD A NEW REGULAR EXPRESSION:
        1. create a new compiled regex pattern, like below
        2. add the name of the compiled regex pattern to the header_regexes dictionary
        3. if the regex groups are new and complicated (parsing more than the accession and organism info),
        alter sequence_info_groups in create_treesapp_ref_data to add another case

    :param code_name: Reference package name/prefix (e.g. DsrAB, p_amoA)
    :return: Dictionary of regular expressions that match database-specific fasta headers indexed by molecule type (str)
    """
    # The regular expressions with the accession and organism name grouped
    # Protein databases:
    gi_re = re.compile(r">?gi\|(\d+)\|[a-z]+\|[_A-Z0-9.]+\|.* RecName: Full=([A-Za-z1-9 _\-]+);?.*$")  # a
    gi_prepend_proper_re = re.compile(r">?gi\|\d+\|[a-z]{2,4}\|([_A-Z0-9.]+)\| (.*) \[(.*)\]$")  # a, d, o
    gi_prepend_mess_re = re.compile(r">?gi\|(\d+)[|/]?.*$")  # a
    dbj_re = re.compile(r">?dbj\|(.*)\|.*\[(.*)\]")  # a, o
    emb_re = re.compile(r">?emb\|(.*)\|.*\[(.*)\]")
    gb_re = re.compile(r">?gb\|(.*)\|.*\[(.*)\]")
    ref_re = re.compile(r">?ref\|(.*)\|.*\[(.*)\]")
    pdb_re = re.compile(r">?pdb\|(.*)\|.+$")  # a
    pir_re = re.compile(r">?pir\|.*\|(\w+).* - (.*)$")  # a, o
    presf_re = re.compile(r">?prf\|.*\|([A-Z0-9]+)\s+.*$")  # a
    sp_re = re.compile(r">?sp\|(.*)\|.*$")  # a
    tr_re = re.compile(r">?tr\|(\w+)\|\w+_\w+ .* OS=(.*) GN=.*$")  # a, o
    treesapp_re = re.compile(r"^>?(\d+)_" + re.escape(code_name) + "$")
    pfam_re = re.compile(r">?([A-Z]+\|)?([A-Z0-9]+)(\.\d)?(_[A-Z0-9]+)?/\d+-\d+$")  # a
    eggnog_re = re.compile(r"^>?(\d+)\.([-\w]+)((\.[\w-]+){0,2})?(?!\s\[.*\])$")  # t, o
    eggnot_re = re.compile(r">?eggnog\|(\d+)\|(.*)")  # a
    # Nucleotide databases:
    # silva_arb_re = re.compile("^>([A-Z0-9]+)\.([0-9]+)\.([0-9]+)_(.*)$")
    # refseq_nuc_re = re.compile("^>([A-Z]+_[0-9]+\.[0-9])_.+$")  # a
    # nr_re = re.compile("^>([A-Z0-9]+\.[0-9])_.*$")  # a

    # Ambiguous:
    # genbank_exact_genome = re.compile("^>([A-Z]{1,2}[0-9]{5,6}\.?[0-9]?) .* \[(.*)\]$")  # a, o
    # accession_only = re.compile(r"^>?([A-Za-z_0-9.]+\.?[0-9]?)$")  # a
    ncbi_ambiguous = re.compile(r"^>?([A-Za-z]{1,6}_[0-9.\-]{5,11})\s+.*(?<!])$")  # a
    ncbi_org = re.compile(r"^>?([A-Z][A-Za-z0-9.\-_]+\.?[0-9]?)\s+(?!lineage=).*\[.*\]$")  # a
    assign_re = re.compile(r"^>?(\w+)?(.*)\|({0})\|(\d+_\d+)$".format(re.escape(code_name)))  # a, d, l

    # Custom fasta header with taxonomy:
    # First group = contig/sequence name, second = full taxonomic lineage, third = description for tree
    # There are no character restrictions on the first and third groups
    # The lineage must be formatted like:
    #   cellular organisms; Bacteria; Proteobacteria; Gammaproteobacteria
    custom_tax = re.compile(r"^>?(.*) lineage=([A-Za-z ]+.*) \[(.*)\]$")  # a, l, o
    unformat_re = re.compile(r"^>?(\w+)?(.*)")

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
                      "ambig": {ncbi_ambiguous: "ncbi_ambig",
                                # accession_only: "bare",
                                ncbi_org: "ncbi_org",
                                custom_tax: "custom",
                                assign_re: "ts_assign",
                                unformat_re: "unformatted"}
                      }
    return header_regexes


def sequence_info_groups(regex_match_groups, header_db: str, header: str, header_regexes=None):
    """
    Depending on the header formats, returns a namedtuple with certain fields filled

    :param regex_match_groups: regular expression (re) match groups
    :param header_db: The name of the assumed database/source of the sequence
    :param header: Header i.e. sequence name that was analyzed
    :param header_regexes: Dictionary of regular expressions matching database-specific headers indexed by molecule type
    :return: namedtuple called seq_info with "description", "locus", "organism", "lineage" and "taxid" fields
    """
    seq_info = namedtuple(typename="seq_info",
                          field_names=["accession", "version", "description", "locus", "organism", "lineage", "taxid"])
    accession = ""
    locus = ""
    organism = ""
    lineage = ""
    taxid = ""
    version = ""

    if regex_match_groups:
        if header_db == "custom":
            lineage = regex_match_groups.group(2)
            organism = regex_match_groups.group(3)
        elif header_db in ["eggnog", "eggnot"]:
            taxid = regex_match_groups.group(1)
            accession = regex_match_groups.group(1) + '.' + regex_match_groups.group(2)
        elif header_db == "ts_assign":
            stripped_header = '|'.join(header.split('|')[:-2])
            header_format_re, stripped_header_db, _ = get_header_format(stripped_header, header_regexes)
            stripped_info = sequence_info_groups(header_format_re.match(stripped_header),
                                                 stripped_header_db, stripped_header)
            accession, version = stripped_info.accession, stripped_info.version
            locus = regex_match_groups.group(3)
        elif header_db == "unformatted":
            if regex_match_groups.group(1):
                accession = regex_match_groups.group(1)
            else:
                accession = re.sub(r"^>", '', header)
        elif header_db == "silva":
            locus = str(regex_match_groups.group(2)) + '-' + str(regex_match_groups.group(3))
            lineage = regex_match_groups.group(4)
        elif header_db == "pfam":
            accession = str(regex_match_groups.group(2))
        elif len(regex_match_groups.groups()) == 3:
            organism = regex_match_groups.group(3)
        elif len(regex_match_groups.groups()) == 2:
            organism = regex_match_groups.group(2)
        if not accession:
            accession = regex_match_groups.group(1)
        if accession.find('.') >= 0:
            pieces = accession.split('.')
            version_match = re.match(r"^(\d{1,2})( (.*)?)?", pieces[1])
            if version_match:
                accession = pieces[0]
                version = '.'.join([accession, version_match.group(1)])

    else:
        logging.error("Unable to handle header: '" + header + "'\n")
        sys.exit(13)

    if not (accession or organism or lineage or taxid):
        logging.error("Insufficient information was loaded for header:\n" +
                      header + "\n" + "regex_match: " + header_db + '\n')
        sys.exit(13)

    if not version:
        version = re.sub(r"^>", '', header.split()[0])

    seq_info = seq_info(accession, version, header, locus, organism, lineage, taxid)

    return seq_info


def get_header_format(header: str, header_regexes: dict) -> (re.compile, str, str):
    """
    Used to decipher which formatting style was used and parse information, ideally reliably

    :param header: A sequences header from a FASTA file
    :param header_regexes: Dictionary of regular expressions matching database-specific headers indexed by molecule type
    :return: Tuple containing the compiled regular expression, matched database name and assumed molecule type
    """
    if re.match(r"^>?[0-9]+\s+coded_by=.+,organism=.+,definition=.+$", header):
        logging.warning(header + " uses GI numbers which are now unsupported by the NCBI! " +
                        "Consider switching to Accession.Version identifiers instead.\n")

    # Gather all possible matches
    format_matches = dict()
    for molecule in header_regexes:
        for regex in header_regexes[molecule]:
            if regex.match(header):
                header_db = header_regexes[molecule][regex]
                format_matches[header_db] = (molecule, regex)
            else:
                pass

    # Exit if there were no matches
    if len(format_matches) == 0:
        logging.error("Unable to parse header '{}'. Unknown format.\n".format(header))
        sys.exit(5)

    # Sort through the matches to find the most specific
    if len(format_matches) == 2 and "unformatted" in format_matches:
        format_matches.pop("unformatted")
    if len(format_matches) > 1:
        if "ts_assign" in format_matches:
            format_matches = {"ts_assign": format_matches["ts_assign"]}  # ts_assign over-rules all others
        else:
            logging.error("Header '{}' matches multiple potential formats:\n\t{}\n"
                          "TreeSAPP is unable to parse necessary information.\n".format(header,
                                                                                        ", ".join(format_matches)))
            sys.exit(5)

    header_db, info = format_matches.popitem()
    header_molecule, header_format_re = info
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


def split_combined_ref_query_fasta(combined_msa, query_msa_file, ref_msa_file) -> None:
    combined_fasta = FASTA(combined_msa)
    combined_fasta.load_fasta()
    seq_names = combined_fasta.get_seq_names()
    write_new_fasta(combined_fasta.fasta_dict, query_msa_file, None,
                    [seq_name for seq_name in seq_names if int(seq_name.split('_')[0]) < 0])
    write_new_fasta(combined_fasta.fasta_dict, ref_msa_file, None,
                    [seq_name for seq_name in seq_names if int(seq_name.split('_')[0]) > 0])
    return


def multiple_alignment_dimensions(mfa_file: str, seq_dict=None) -> (int, int):
    """
    Checks to ensure all sequences are the same length and returns a tuple of (nrow, ncolumn)

    :param seq_dict: A dictionary containing headers as keys and sequences as values
    :param mfa_file: The name of the multiple alignment FASTA file being validated
    :return: tuple = (nrow, ncolumn)
    """
    if not seq_dict:
        seq_dict = read_fasta_to_dict(mfa_file)
    sequence_length = 0
    for seq_name in seq_dict:
        sequence = seq_dict[seq_name]
        if sequence_length == 0:
            sequence_length = len(sequence)
        elif sequence_length != len(sequence) and sequence_length > 0:
            logging.error("Number of aligned columns is inconsistent in " + mfa_file + "!\n")
            sys.exit(3)
        else:
            pass
            # Sequence is the right length, carrying on
    return len(seq_dict), sequence_length
