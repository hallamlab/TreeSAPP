__author__ = 'Connor Morgan-Lang'

import sys
import re

import _fasta_reader


def format_read_fasta(fasta_input, molecule, args, max_header_length=110):
    """
    Reads a FASTA file, ensuring each sequence and sequence name is valid.
    :param fasta_input: Absolute path of the FASTA file to be read
    :param molecule: Molecule type of the sequences ['prot', 'dna', 'rrna']
    :param args: Command-line argument object from get_options and check_parser_arguments
    :param max_header_length: The length of the header string before all characters after this length are removed
    :return A Python dictionary with headers as keys and sequences as values
    """

    if sys.version_info > (2, 9):
        py_version = 3
    else:
        py_version = 2
        from itertools import izip

    fasta_list = _fasta_reader._read_format_fasta(fasta_input,
                                                  args.min_seq_length,
                                                  args.output,
                                                  molecule,
                                                  max_header_length)
    if not fasta_list:
        sys.exit()
    tmp_iterable = iter(fasta_list)
    if py_version == 2:
        formatted_fasta_dict = dict(izip(tmp_iterable, tmp_iterable))
    elif py_version == 3:
        formatted_fasta_dict = dict(zip(tmp_iterable, tmp_iterable))
    else:
        raise AssertionError("Unexpected Python version detected")

    for header in formatted_fasta_dict.keys():
        if len(header) > max_header_length:
            sys.stderr.write(header + " is too long!\nThere is a bug in _read_format_fasta - please report!\n")
            sys.exit()

    return formatted_fasta_dict


def get_headers(fasta_file):
    original_headers = list()
    try:
        fasta = open(fasta_file, 'r')
    except IOError:
        raise IOError("ERROR: Unable to open the FASTA file (" + fasta_file + ") provided for reading!")
    line = fasta.readline()
    while line:
        line = line.strip()
        if line and line[0] == '>':
            original_headers.append(line)
        else:
            pass
        line = fasta.readline()

    fasta.close()
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
    file_counter = 0
    sequence_accumulator = 0

    if max_seqs is not None:
        fasta_name = fasta_name + '_' + str(max_seqs)

    try:
        fa_out = open(fasta_name, 'w')
    except IOError:
        raise IOError("Unable to open " + fasta_name + " for writing!")

    for name in fasta_dict.keys():
        seq = fasta_dict[name]
        if name[0] != '>':
            name = '>' + name
        sequence_accumulator += 1
        if max_seqs and sequence_accumulator > max_seqs:
            # If input is to be split and number of sequences per file has been exceeded begin writing to new file
            fa_out.close()
            split_files.append(fasta_name)
            file_counter += 1
            sequence_accumulator = 1
            fasta_name = re.sub(r'_d+$', '_' + str(file_counter), fasta_name)
            fa_out = open(fasta_name, 'w')

        if headers is None:
            fa_out.write(name + "\n")
            fa_out.write(seq + "\n")
        elif name[1:] in headers:
            fa_out.write(name + "\n")
            fa_out.write(seq + "\n")

    fa_out.close()
    split_files.append(fasta_name)
    file_counter += 1
    return split_files


def get_header_format(header, code_name=""):
    """
    Used to decipher which formatting style was used: NCBI, FunGenes, or other
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
    gi_re = re.compile(">gi\|(\d+)\|[a-z]+\|?\w.+\|(.*)$")
    gi_prepend_proper_re = re.compile(">gi\|([0-9]+)\|[a-z]+\|[_A-Z0-9.]+\|.*\[(.*)\]$")
    gi_prepend_mess_re = re.compile(">gi\|([0-9]+)\|pir\|\|(.*)$")
    dbj_re = re.compile(">dbj\|(.*)\|.*\[(.*)\]")
    emb_re = re.compile(">emb\|(.*)\|.*\[(.*)\]")
    gb_re = re.compile(">gb\|(.*)\|.*\[(.*)\]")
    ref_re = re.compile(">ref\|(.*)\|.*\[(.*)\]")
    pdb_re = re.compile(">pdb\|(.*)\|(.*)$")
    pir_re = re.compile(">pir\|\|(\w+).* - (.*)$")
    sp_re = re.compile(">sp\|(.*)\|.*Full=(.*); AltName:.*$")
    fungene_re = re.compile("^>([A-Z0-9.]+)[ ]+coded_by=(.+)[,]+organism=(.+)[,]+definition=(.+)$")
    fungene_trunc_re = re.compile("^>([A-Z0-9.]+)[ ]+organism=(.+)[,]+definition=(.+)$")
    mltree_re = re.compile("^>(\d+)_" + re.escape(code_name))
    treesapp_re = re.compile("^>([A-Z0-9.]+) .* \[(.*)\]$")
    refseq_prot_re = re.compile("^>([A-Z]{2}_[0-9]+\.[0-9]) (.*) \[(.*)\]$")

    # Nucleotide databases:
    genbank_exact_genome = re.compile("^>([A-Z0-9]+\.[0-9]) (.*) \[(.*)\]$")
    ncbi_ambiguous = re.compile("^>([A-Z0-9]+\.[0-9]) (.*)$")
    silva_arb_re = re.compile("^>([A-Z0-9]+)\.([0-9]+)\.([0-9]+)_(.*)$")
    refseq_nuc_re = re.compile("^>([A-Z]+_[0-9]+\.[0-9])_(.*)$")
    nr_re = re.compile("^>([A-Z0-9]+\.[0-9])_(.*)$")

    header_regexes = {"prot": {dbj_re: "dbj",
                               emb_re: "emb",
                               gb_re: "gb",
                               pdb_re: "pdb",
                               pir_re: "pir",
                               ref_re: "ref",
                               sp_re: "sp",
                               gi_re: "gi_re",
                               gi_prepend_proper_re: "gi_proper",
                               gi_prepend_mess_re: "gi_mess",
                               refseq_prot_re: "refseq_prot"},
                      "dna": {fungene_re: "fungene",
                              fungene_trunc_re: "fungene_truncated",
                              mltree_re: "mltree",
                              silva_arb_re: "silva",
                              refseq_nuc_re: "refseq_nuc",
                              nr_re: "nr"},
                      "ambig": {ncbi_ambiguous: "ncbi_ambig",
                                genbank_exact_genome: "gen_genome",
                                treesapp_re: "treesapp"}
                      }

    header_format_re = None
    header_db = None
    header_molecule = None
    for molecule in header_regexes:
        for regex in header_regexes[molecule]:
            if regex.match(header):
                header_format_re = regex
                header_db = header_regexes[molecule][regex]
                header_molecule = molecule
                break
            else:
                pass

    if header_format_re is None:
        raise AssertionError("Unable to parse header: " + header)

    return header_format_re, header_db, header_molecule


