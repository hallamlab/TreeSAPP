#!/usr/bin/env python

import argparse
import re
import sys


def get_options():
    parser = argparse.ArgumentParser(description="Used to replace shortened headers with their respective complete headers in a FASTA file.")
    parser.add_argument("--stunted")
    parser.add_argument("--full")
    parser.add_argument("-o", "--output", help="The output FASTA file", required=True)
    parser.add_argument("-m", "--molecule", required=True)
    parser.add_argument("-l", "--min_length", required=False, default=0)
    parser.add_argument("-v", "--verbose", required=False, default=True)
    args = parser.parse_args()
    return args


def format_read_fasta(args, duplicates=False):
    """
    Splits the input file into multiple files, each containing a maximum number of sequences as specified by the user.
    Ensures each sequence and sequence name is valid.
    :param args: Command-line argument object from get_options and check_parser_arguments
    :param duplicates: A flag indicating the function should be duplicate-aware
    :return A list of the files produced from the input file.
    """
    if args.verbose:
        sys.stdout.write("Formatting " + args.input + " for pipeline... ")
        sys.stdout.flush()

    fasta = open(args.input, 'r')
    formatted_fasta_dict = dict()
    header = ""
    sequence = ""
    convert_upper = lambda pat: pat.group(1).upper()
    reg_nuc = re.compile(r'[acgtACGT]')
    reg_nuc_ambiguity = re.compile(r'[xnXN]')
    reg_aa_ambiguity = re.compile(r'[bjxzBJXZ]')
    reg_amino = re.compile(r'[acdefghiklmnpqrstuvwyACDEFGHIKLMNPQRSTUVWY*]')
    iupac_map = {'R': 'A', 'Y': 'C', 'S': 'G', 'W': 'A', 'K': 'G', 'M': 'A', 'B': 'C', 'D': 'A', 'H': 'A', 'V': 'A'}
    substitutions = ""
    count_total = 0
    seq_count = 0
    count_xn = 0
    header_clash = False
    num_headers = 0
    if duplicates:
        duplicate_headers = dict()
    for line in fasta:
        # If the line is a header...
        if line[0] == '>':
            if len(header) > 0:
                if len(sequence) > args.min_length:
                    formatted_fasta_dict[header] = sequence
                else:
                    num_headers -= 1

            sequence = ""
            # Replace all non a-z, A-Z, 0-9, or . characters with a _
            # Then replace the initial _ with a >
            # line = re.sub(r'[^a-zA-Z0-9.\r\n]', '_', line)
            # line = re.sub(r'\A_', '>', line)

            # Because RAxML can only work with file names having length <= 125,
            # Ensure that the sequence name length is <= 100
            # if line.__len__() > 125:
            #     line = line[0:120]

            header = line.strip()
            if duplicates:
                if header in formatted_fasta_dict.keys():
                    if header not in duplicate_headers.keys():
                        duplicate_headers[header] = 1
                    duplicate_headers[header] += 1
                    header = header + "_" + str(duplicate_headers[header])
            num_headers += 1

        # Else, if the line is a sequence...
        else:
            characters = line.strip()
            if len(characters) == 0:
                continue
            # Remove all non-characters from the sequence
            re.sub(r'[^a-zA-Z]', '', characters)

            # Count the number of {atcg} and {xn} in all the sequences
            count_total += len(characters)

            if args.molecule == 'n':
                nucleotides = len(reg_nuc.findall(characters))
                ambiguity = len(reg_nuc_ambiguity.findall(characters))
                if (nucleotides + ambiguity) != len(characters):
                    substituted_chars = ""
                    for char in characters:
                        if not reg_nuc.search(char):
                            if char not in iupac_map.keys():
                                sys.stderr.write("ERROR: " + header.strip() + " contains unknown characters!\n")
                                sys.exit()
                            else:
                                substituted_chars += iupac_map[char]
                                substitutions += char
                        else:
                            substituted_chars += char
                    characters = substituted_chars
                else:
                    seq_count += nucleotides
                    count_xn += ambiguity
            elif args.molecule == 'a':
                characters = re.sub(r'([a-z])', convert_upper, characters)
                aminos = len(reg_amino.findall(characters))
                ambiguity = len(reg_aa_ambiguity.findall(characters))
                if (aminos + ambiguity) != len(characters):
                    sys.stderr.write("ERROR: " + header.strip() + " contains unknown characters: ")
                    unknown = ""
                    for c in characters:
                        if c not in "abcdefghiklmnpqrstvwxyzABCDEFGHIKLMNPQRSTVWXYZ":
                            unknown += c
                    sys.stderr.write(unknown + "\n")
                    sys.exit()
                else:
                    seq_count += aminos
                    count_xn += ambiguity
            sequence += characters
    formatted_fasta_dict[header] = sequence

    # Check for duplicate headers and change them to make unique headers
    if len(formatted_fasta_dict.keys()) != num_headers:
        header_clash = True

    if count_total == 0:
        sys.exit('ERROR: Your input file appears to be corrupted. No sequences were found!\n')

    # Exit the program if all sequences are composed only of X or N
    elif count_xn == count_total:
        sys.exit('ERROR: Your sequence(s) contain only X or N!\n')

    # Exit the program if less than half of the characters are nucleotides
    elif float(seq_count / float(count_total)) < 0.5:
        if args.molecule == 'n':
            sys.exit('ERROR: Your sequence(s) most likely contain no DNA!\n')
        elif args.molecule == 'a':
            sys.exit('ERROR: Your sequence(s) most likely contain no AA!\n')

    if len(substitutions) > 0:
        sys.stderr.write("WARNING: " + str(len(substitutions)) + " ambiguous character substitutions were made!\n")
        sys.stderr.flush()

    # Close the files
    fasta.close()
    if args.verbose:
        sys.stdout.write("done.\n")

    if header_clash:
        sys.stderr.write("WARNING: duplicate header names were detected in " + args.input + "!\n")
        sys.stderr.write("Formatting input FASTA with duplicates... ")
        sys.stderr.flush()
        formatted_fasta_dict = format_read_fasta(args, True)
        sys.stderr.write("done.\n")
        sys.stderr.flush()

    return formatted_fasta_dict


def find_replace_headers(stunted, full):
    new_fa = dict()
    for short_head in stunted.keys():
        for header in full.keys():
            if short_head.split('_')[0] == header.split(' ')[0]:
                if stunted[short_head] == full[header]:
                    new_fa[header] = stunted[short_head]
    return new_fa


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
    except:
        raise IOError("Unable to open " + fasta_name + " for writing!")

    for name in fasta_dict.keys():
        seq = fasta_dict[name]
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


def main():
    args = get_options()
    args.input = args.stunted
    stunted_fa = format_read_fasta(args)
    print "Stunted contains", len(stunted_fa), "sequences"
    args.input = args.full
    full_fa = format_read_fasta(args)
    print "Full contains", len(full_fa), "sequences"
    new_fa = find_replace_headers(stunted_fa, full_fa)
    write_new_fasta(new_fa, args.output)

main()
