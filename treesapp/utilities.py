__author__ = 'Connor Morgan-Lang'

import os
import re
import sys
import subprocess
import logging
import time
from collections import namedtuple
from shutil import rmtree
from .external_command_interface import launch_write_command
from pygtrie import StringTrie


def load_taxonomic_trie(lineages: list) -> StringTrie:
    taxonomic_tree = StringTrie(separator='; ')

    for lineage in lineages:
        i = 0
        ranks = len(lineage)
        while i < len(lineage):
            taxonomic_tree["; ".join(lineage.split("; ")[:ranks - i])] = True
            i += 1

    return taxonomic_tree


def prepend_deep_rank(seq_lineage_dict):
    cellular_organisms = ["Bacteria", "Archaea", "Eukaryota"]
    for seq_name in seq_lineage_dict:
        lineage = seq_lineage_dict[seq_name]
        if lineage.split("; ")[0] in cellular_organisms:
            seq_lineage_dict[seq_name] = "cellular organisms; " + lineage
    return


def rekey_dict(og_dict: dict, key_map: dict) -> dict:
    """
    Creates a new dictionary with new keys, indicated by a map, mapped to the original values.
    Logs a warning if not all of the original keys are popped.
    :param og_dict: The original dictionary with keys found in key_map. Values are retained.
    :param key_map: A dictionary mapping old keys to new keys
    :return: Dictionary with new keys, same values
    """
    updated_dict = dict()
    unmapped = list()

    if len(og_dict) != len(key_map):
        logging.error("Key map (" + str(len(key_map)) + ") and original dictionary (" + str(len(og_dict)) +
                      ") are different sizes. Unable to re-key.\n")
        sys.exit(5)

    og_keys = sorted(list(og_dict.keys()))
    for old_key in og_keys:
        try:
            new_key = key_map[old_key]
        except KeyError:
            unmapped.append(old_key)
            continue
        updated_dict[new_key] = og_dict.pop(old_key)

    if len(unmapped) > 0:
        logging.warning("Dictionary rekey incomplete - " + str(len(unmapped)) + " keys were not found in dictionary.\n")

    if len(og_dict) > 0:
        logging.warning(str(len(og_dict)) + " keys in original dictionary were not popped during re-keying process.\n")

    return updated_dict


def reluctant_remove_replace(dir_path):
    # DO NOT USE LOGGING - this function can (and does) appear before the logger is instantiated
    # Warn user then remove all main output directories, leaving log in output
    sys.stderr.write("WARNING:\nRemoving previous outputs in '" + dir_path + "'. " +
                     "You have 10 seconds to hit Ctrl-C before this proceeds.\n")
    time.sleep(10)
    rmtree(dir_path)
    os.mkdir(dir_path)
    return


def get_refpkg_build(name: str, marker_build_dict: dict, refpkg_code_re):
    """
    Find and return the MarkerBuild instance with a matching name
    :param name:
    :param marker_build_dict: A dictionary of MarkerBuild objects indexed by their refpkg codes/denominators
    :param refpkg_code_re: A compiled regular expression (re) for matching refpkg code names
    :return: MarkerBuild
    """
    if refpkg_code_re.match(name):
        try:
            return marker_build_dict[name]
        except KeyError:
            # Alert the user if the denominator format was (incorrectly) provided to Clade_exclusion_analyzer
            logging.error("Unable to find '" + name + "' in ref_build_parameters collection!\n" +
                          "Has it been added to data/tree_data/ref_build_parameters.tsv?\n")
            sys.exit(21)
    else:
        for denominator in marker_build_dict:
            refpkg_build = marker_build_dict[denominator]
            if refpkg_build.cog == name:
                return refpkg_build
            elif name == denominator:
                return refpkg_build
        logging.error("Wrong format for the reference code_name provided: " + name + "\n")
        sys.exit(21)


def is_exe(fpath):
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)


def which(program):
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path_element in os.environ["PATH"].split(os.pathsep):
            path_element = path_element.strip('"')
            exe_file = os.path.join(path_element, program)
            if is_exe(exe_file):
                return exe_file
    return None


def os_type():
    """Return the operating system of the user."""
    x = sys.platform
    if x:

        hits = re.search(r'darwin', x, re.I)
        if hits:
            return 'mac'

        hits = re.search(r'win', x, re.I)
        if hits:
            return 'win'

        hits = re.search(r'linux', x, re.I)
        if hits:
            return 'linux'


def available_cpu_count():
    """ Number of available virtual or physical CPUs on this system, i.e.
    user/real as output by time(1) when called with an optimally scaling
    userspace-only program"""

    # cpuset
    # cpuset may restrict the number of *available* processors
    try:
        m = re.search(r'(?m)^Cpus_allowed:\s*(.*)$',
                      open('/proc/self/status').read())
        if m:
            res = bin(int(m.group(1).replace(',', ''), 16)).count('1')
            if res > 0:
                return res
    except IOError:
        pass

    # Python 2.6+
    try:
        import multiprocessing
        return multiprocessing.cpu_count()
    except (ImportError, NotImplementedError):
        pass

    # http://code.google.com/p/psutil/
    try:
        import psutil
        return psutil.NUM_CPUS
    except (ImportError, AttributeError):
        pass

    # POSIX
    try:
        res = int(os.sysconf('SC_NPROCESSORS_ONLN'))

        if res > 0:
            return res
    except (AttributeError, ValueError):
        pass

    # Windows
    try:
        res = int(os.environ['NUMBER_OF_PROCESSORS'])

        if res > 0:
            return res
    except (KeyError, ValueError):
        pass

    # BSD
    try:
        sysctl = subprocess.Popen(['sysctl', '-n', 'hw.ncpu'],
                                  stdout=subprocess.PIPE)
        scStdout = sysctl.communicate()[0]
        res = int(scStdout)

        if res > 0:
            return res
    except (OSError, ValueError):
        pass

    # Linux
    try:
        res = open('/proc/cpuinfo').read().count('processor\t:')

        if res > 0:
            return res
    except IOError:
        pass

    # Solaris
    try:
        pseudoDevices = os.listdir('/devices/pseudo/')
        res = 0
        for pd in pseudoDevices:
            if re.match(r'^cpuid@[0-9]+$', pd):
                res += 1

        if res > 0:
            return res
    except OSError:
        pass

    # Other UNIXes (heuristic)
    try:
        try:
            dmesg = open('/var/run/dmesg.boot').read()
        except IOError:
            dmesgProcess = subprocess.Popen(['dmesg'], stdout=subprocess.PIPE)
            dmesg = dmesgProcess.communicate()[0]

        res = 0
        while '\ncpu' + str(res) + ':' in dmesg:
            res += 1

        if res > 0:
            return res
    except OSError:
        pass

    logging.error('Can not determine number of CPUs on this system')


def executable_dependency_versions(exe_dict):
    """
    Function for retrieving the version numbers for each executable in exe_dict
    :param exe_dict: A dictionary mapping names of software to the path to their executable
    :return: A formatted string with the executable name and its respective version found
    """
    versions_dict = dict()
    versions_string = "Software versions used:\n"

    simple_v = ["prodigal", "raxmlHPC"]
    version_param = ["trimal", "mafft"]
    no_params = ["usearch", "papara"]
    help_param = ["hmmbuild", "hmmalign", "hmmsearch", "OD-seq"]
    version_re = re.compile(r"[Vv]\d+.\d|version \d+.\d|\d\.\d\.\d|HMMER")

    for exe in exe_dict:
        ##
        # Get the help/version statement for the software
        ##
        versions_dict[exe] = ""
        if exe in simple_v:
            stdout, returncode = launch_write_command([exe_dict[exe], "-v"])
        elif exe in version_param:
            stdout, returncode = launch_write_command([exe_dict[exe], "--version"])
        elif exe in help_param:
            stdout, returncode = launch_write_command([exe_dict[exe], "-h"])
        elif exe in no_params:
            stdout, returncode = launch_write_command([exe_dict[exe]])
        elif exe == "FastTree":
            stdout, returncode = launch_write_command([exe_dict[exe], "-expert"])
        elif exe == "BMGE.jar":
            stdout, returncode = launch_write_command(["java", "-Xmx10m", "-jar", exe_dict[exe], "-?"])
        else:
            logging.warning("Unknown version command for " + exe + ".\n")
            continue
        ##
        # Identify the line with the version number (since often more than a single line is returned)
        ##
        for line in stdout.split("\n"):
            if version_re.search(line):
                # If a line was identified, try to get just the string with the version number
                for word in line.split(" "):
                    if re.search(r"\d\.\d", word):
                        versions_dict[exe] = re.sub(r"[,:()[\]]", '', word)
                        break
                break
            else:
                pass
        if not versions_dict[exe]:
            logging.debug("Unable to find version for " + exe + ".\n")

    ##
    # Format the string with the versions of all software
    ##
    for exe in sorted(versions_dict):
        n_spaces = 12-len(exe)
        versions_string += "\t" + exe + ' '*n_spaces + versions_dict[exe] + "\n"

    return versions_string


def reformat_string(string):
    if string and string[0] == '>':
        header = True
    else:
        header = False
    string = re.sub(r"\[|\]|\(|\)|\/|\\\\|'|<|>", '', string)
    if header:
        string = '>' + string
    string = re.sub(r"\s|;|,|\|", '_', string)
    if len(string) > 110:
        string = string[0:109]
    while string and string[-1] == '.':
        string = string[:-1]
    return string


class Autovivify(dict):
    """In cases of Autovivify objects, enable the referencing of variables (and sub-variables)
    without explicitly declaring those variables beforehand."""

    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value


def return_sequence_info_groups(regex_match_groups, header_db: str, header: str):
    """
    Depending on the header formats, returns a namedtuple with certain fields filled
    :param regex_match_groups: regular expression (re) match groups
    :param header_db: The name of the assumed database/source of the sequence
    :param header: Header i.e. sequence name that was analyzed
    :return: namedtuple called seq_info with "description", "locus", "organism", "lineage" and "taxid" fields
    """
    seq_info = namedtuple(typename="seq_info",
                          field_names=["accession", "version", "description", "locus", "organism", "lineage", "taxid"])
    description = ""
    locus = ""
    organism = ""
    lineage = ""
    taxid = ""
    version = ""

    if regex_match_groups:
        accession = regex_match_groups.group(1)
        if accession.find('.'):
            pieces = accession.split('.')
            if len(pieces) == 2 and re.match(r"\d+", pieces[1]):
                version = accession
                accession = pieces[0]
        if header_db == "custom":
            lineage = regex_match_groups.group(2)
            organism = regex_match_groups.group(3)
            description = regex_match_groups.group(3)
        elif header_db in ["eggnog", "eggnot"]:
            taxid = regex_match_groups.group(1)
            accession = regex_match_groups.group(1) + '.' + regex_match_groups.group(2)
        elif header_db == "ts_assign":
            accession = '|'.join(regex_match_groups.groups())
            description = regex_match_groups.group(1)
            locus = regex_match_groups.group(3)
        elif header_db == "silva":
            locus = str(regex_match_groups.group(2)) + '-' + str(regex_match_groups.group(3))
            lineage = regex_match_groups.group(4)
            description = regex_match_groups.group(4)
        elif len(regex_match_groups.groups()) == 3:
            description = regex_match_groups.group(2)
            organism = regex_match_groups.group(3)
        elif len(regex_match_groups.groups()) == 2:
            organism = regex_match_groups.group(2)
    else:
        logging.error("Unable to handle header: " + header + "\n")
        sys.exit(13)

    if not (accession or organism or lineage or taxid):
        logging.error("Insufficient information was loaded for header:\n" +
                      header + "\n" + "regex_match: " + header_db + '\n')
        sys.exit(13)
    if not accession:
        accession = re.sub(r"^>", '', header)
    seq_info = seq_info(accession, version, description, locus, organism, lineage, taxid)

    return seq_info


def remove_dashes_from_msa(fasta_in, fasta_out):
    """
    fasta_out is the new FASTA file written with no dashes (unaligned)
    There are no line breaks in this file, whereas there may have been in fasta_in
    :param fasta_in: Multiply-aligned FASTA file
    :param fasta_out: FASTA file to write
    :return:
    """
    dashed_fasta = open(fasta_in, 'r')
    fasta = open(fasta_out, 'w')
    sequence = ""

    line = dashed_fasta.readline()
    while line:
        if line[0] == '>':
            if sequence:
                fasta.write(sequence + "\n")
                sequence = ""
            fasta.write(line)
        else:
            sequence += re.sub('[-.]', '', line.strip())
        line = dashed_fasta.readline()
    fasta.write(sequence + "\n")
    dashed_fasta.close()
    fasta.close()
    return


def clean_lineage_string(lineage: str, also=None) -> str:
    """
    Removes superfluous taxonomic ranks and characters that make lineage comparisons difficult

    :param lineage: A taxonomic lineage string where each rank is separated by a semi-colon
    :param also: An optional list containing other strings to filter out e.g. ["Root; "]
    :return: String with the purified taxonomic lineage adhering to the NCBI hierarchy
    """
    non_standard_names_re = re.compile(" group| cluster| complex", re.IGNORECASE)
    bad_strings = ["cellular organisms; ", "delta/epsilon subdivisions; ", "\(miscellaneous\)", "[a-p]__"]
    if also:
        bad_strings += also
    for bs in bad_strings:
        lineage = re.sub(bs, '', lineage)
    # filter 'group' and 'cluster'
    if non_standard_names_re.search(lineage):
        reconstructed_lineage = ""
        ranks = lineage.split("; ")
        for rank in ranks:
            if not non_standard_names_re.search(rank):
                reconstructed_lineage = reconstructed_lineage + str(rank) + '; '
        reconstructed_lineage = re.sub('; $', '', reconstructed_lineage)
        lineage = reconstructed_lineage
    return re.sub(r"^; ", '', lineage)


def remove_elongated_lineages(fasta_records: dict, guarantees=None):
    """
    Meant to address the issues caused by unreasonable long lineages -- those that are expanded beyond the canonical 7
    (R; K; P; C; O; F; G; S). Since taxon labels are not linked or validated with a database we simply infer their
    rank from their position in a lineage.
    If a lineage has sub-ranks that haven't been filtered by `clean_lineage_string()` this can cause a problem in
    `treesapp train` as the phylogenetic distance data wouldn't reflect the rank its supposed to, and `treesapp evaluate`
    where the sequences being evaluated wouldn't reflect the performance of the target rank.
    :param fasta_records: Dictionary mapping numerical treesapp_id values to EntrezRecord instances
    :param guarantees: Set object containing treesapp_id values corresponding to sequences that need to be retained
    :return:
    """
    num_removed = 0
    fasta_replace_dict = dict()

    for treesapp_id in fasta_records:
        ref_seq = fasta_records[treesapp_id]
        if guarantees and treesapp_id in guarantees:
            fasta_replace_dict[treesapp_id] = ref_seq
            continue
        # Not much can be done here without searching for each taxon's rank
        if len(clean_lineage_string(ref_seq.lineage, ["Root; "]).split("; ")) >= 7:
            lp = clean_lineage_string(ref_seq.lineage, ["Root; ", "Candidatus ", "sp. "]).split("; ")
            # Check if the Species-position matches a typical species RE
            if not re.match(r"^[A-z][a-z]+ [.a-z0-9]+$", lp[6]):
                num_removed += 1
            else:
                fasta_replace_dict[treesapp_id] = ref_seq
        else:
            fasta_replace_dict[treesapp_id] = ref_seq

    logging.debug('\t' + str(num_removed) + " sequences removed with elongated taxonomic lineages.\n" +
                  '\t' + str(len(fasta_replace_dict) - num_removed) + " sequences retained for building tree.\n")

    return fasta_replace_dict


def median(num_list: list):
    n = len(num_list)
    if n < 1:
            return None
    if n % 2 == 1:
            return sorted(num_list)[n//2]
    else:
            return sum(sorted(num_list)[n//2-1:n//2+1])/2.0


def mean(num_list: list):
    """
    Simple function for a returning a floating-point integer for the mean of a list of numbers
    :param num_list: List of numbers
    :return: Float
    """
    return float(sum(num_list) / len(num_list))


def convert_outer_to_inner_nodes(clusters: dict, internal_node_map: dict):
    """
    Find the lowest common ancestor (internal node) for all leaves in the range.
    This is only necessary if the original nodes parsed from the colours_style.txt file were leaves.
    :param clusters: A dictionary mapping start and end leaves of a clade for a single marker's colours_style.txt layer
    :param internal_node_map: A dictionary mapping each internal node to a list of all of its descendent leaves
    :return:
    """
    leaf_annotation_map = dict()
    for cluster in clusters.keys():
        if cluster not in leaf_annotation_map:
            leaf_annotation_map[cluster] = list()
        for frond_tips in clusters[cluster]:
            start, end = frond_tips
            # Find the minimum set that includes both start and end
            warm_front = dict()
            # Add all the potential internal nodes
            for inode in internal_node_map:
                clade = internal_node_map[inode]
                if start in clade:
                    warm_front[inode] = clade
            for inode in sorted(warm_front, key=lambda x: len(warm_front[x])):
                if end in warm_front[inode]:
                    leaf_annotation_map[cluster].append((inode, inode))
                    break
    return leaf_annotation_map


def reformat_fasta_to_phy(fasta_dict):
    phy_dict = dict()
    for seq_name in fasta_dict:
        sequence = fasta_dict[seq_name]
        sub_sequences = re.findall(r'.{1,50}', sequence)
        count = 0
        for sub_sequence in sub_sequences:
            if count not in phy_dict:
                phy_dict[count] = dict()
            phy_dict[count][seq_name] = sub_sequence
            count += 1
    return phy_dict


def write_phy_file(phy_output_file: str, phy_dict: dict, alignment_dims=None):
    """
    Writes a Phylip-formatted alignment file
    PaPaRa is EXTREMELY particular about the input of its Phylip file. Don't mess.
    :param phy_output_file: File path to write the Phylip file
    :param phy_dict: Dictionary of sequences to write
    :param alignment_dims: Tuple containing (num_seqs, alignment_len)
    :return:
    """
    if not alignment_dims:
        seq_chunks = [len(aligned_seq) for aligned_seq in phy_dict.values()]
        if min(seq_chunks) != max(seq_chunks):
            logging.error("Inconsistent number of sequences in Phylip dictionary keys.")
        num_seqs = seq_chunks[0]
        aligned_seqs = dict()
        for seq_chunk in phy_dict.keys():
            for seq_name in phy_dict[seq_chunk].keys():
                if seq_name not in aligned_seqs:
                    aligned_seqs[seq_name] = ""
                aligned_seqs[seq_name] += phy_dict[seq_chunk][seq_name]
        aligned_seq_lengths = [len(aligned_seqs[aligned_seq]) for aligned_seq in aligned_seqs.keys()]
        if min(aligned_seq_lengths) != max(aligned_seq_lengths):
            logging.error("Lengths of aligned sequences are heterogeneous.")
        alignment_len = aligned_seq_lengths[0]
        aligned_seqs.clear()
        seq_chunks.clear()
    else:
        num_seqs, alignment_len = alignment_dims

    with open(phy_output_file, 'w') as phy_output:
        phy_string = ' ' + str(num_seqs) + ' ' + str(alignment_len) + '\n'
        for count in sorted(phy_dict.keys(), key=int):
            for seq_name in sorted(phy_dict[count].keys()):
                sequence_part = re.sub('X', '-', phy_dict[count][seq_name])
                if count == 0:
                    phy_string += str(seq_name)
                    length = len(str(seq_name))
                    c = length
                    while c < 11:
                        phy_string += ' '
                        c += 1
                else:
                    phy_string += 11*' '
                phy_string += ' '.join(re.findall(r'.{1,10}', sequence_part)) + '\n'
            phy_string += "\n"

        phy_output.write(phy_string)
    return


def extract_hmm_matches(hmm_matches, fasta_dict, header_registry):
    """
    Function for slicing sequences guided by alignment co-ordinates.

    :param hmm_matches: Dictionary containing a list HmmMatch() objects as values for each 'marker' key
    :param fasta_dict: A dictionary with headers as keys and sequences as values
    :param header_registry: A list of Header() objects, each used to map various header formats to each other
    :return:
    """

    if len(hmm_matches.keys()) > 1:
        logging.error("Number of markers found from HMM alignments is >1\n" +
                      "Does your HMM file contain more than 1 profile? TreeSAPP is unprepared for this.\n")
        sys.exit(13)

    marker_gene_dict = dict()
    header_matching_dict = dict()

    logging.debug("Creating a temporary dictionary for rapid sequence name look-ups... ")
    for num in header_registry:
        if header_registry[num].first_split[0] == '>':
            header_matching_dict[header_registry[num].first_split[1:]] = header_registry[num]
        else:
            header_matching_dict[header_registry[num].first_split] = header_registry[num]
    logging.debug("done.\n")

    logging.info("Extracting the quality-controlled protein sequences... ")

    for marker in hmm_matches:
        if marker not in marker_gene_dict:
            marker_gene_dict[marker] = dict()

        for hmm_match in hmm_matches[marker]:
            # Now for the header format to be used in the bulk FASTA:
            # >contig_name|marker_gene|start_end
            query_names = header_matching_dict[hmm_match.orf]
            try:
                sequence = fasta_dict[query_names.first_split]
            except KeyError:
                logging.debug("Unable to map " + hmm_match.orf + " to a sequence in the input FASTA.\n")
                continue
            if hmm_match.of > 1:
                query_names.post_align = ' '.join([query_names.first_split,
                                                   str(hmm_match.num) + '.' + str(hmm_match.of),
                                                   re.sub(re.escape(query_names.first_split), '', query_names.original)])
            else:
                query_names.post_align = query_names.original
            bulk_header = query_names.post_align

            if bulk_header in marker_gene_dict[marker]:
                logging.warning(bulk_header + " being overwritten by an alternative alignment!\n" + hmm_match.get_info())
            marker_gene_dict[marker][bulk_header] = sequence[hmm_match.start-1:hmm_match.end]

    logging.info("done.\n")
    return marker_gene_dict[marker]


def hmm_pile(hmm_matches):
    """
    Function to inspect the placement of query sequences on the reference HMM
    :param hmm_matches:
    :return:
    """
    hmm_bins = dict()
    window_size = 2

    for marker in hmm_matches:
        for hmm_match in hmm_matches[marker]:
            # Initialize the hmm_bins using the HMM profile length
            if not hmm_bins:
                i = 1
                hmm_length = int(hmm_match.hmm_len)
                while i < hmm_length:
                    if i+window_size-1 > hmm_length:
                        hmm_bins[(i, hmm_length)] = 0
                    else:
                        hmm_bins[(i, i+window_size-1)] = 0
                    i += window_size
            # Skip ahead to the HMM profile position where the query sequence began aligning
            for bin_start, bin_end in hmm_bins:
                if hmm_match.pstart <= bin_start and hmm_match.pend >= bin_end:
                    hmm_bins[(bin_start, bin_end)] += 1
                else:
                    pass
        low_coverage_start = 0
        low_coverage_stop = 0
        maximum_coverage = 0
        low_cov_summary = ""
        for window in sorted(hmm_bins.keys()):
            height = hmm_bins[window]
            if height > maximum_coverage:
                maximum_coverage = height
            if height < len(hmm_matches[marker])/2:
                begin, end = window
                if low_coverage_start == low_coverage_stop:
                    low_coverage_start = begin
                    low_coverage_stop = end
                else:
                    low_coverage_stop = end
            elif height > len(hmm_matches[marker])/2 and low_coverage_stop != 0:
                low_cov_summary += "\t" + str(low_coverage_start) + '-' + str(low_coverage_stop) + "\n"
                low_coverage_start = 0
                low_coverage_stop = 0
            else:
                pass
        if low_coverage_stop != low_coverage_start:
            low_cov_summary += "\t" + str(low_coverage_start) + "-end\n"

        if low_cov_summary:
            logging.info("Low coverage HMM windows (start-stop):\n" + low_cov_summary)
        logging.info("Maximum coverage = " + str(maximum_coverage) + " sequences\n")
    return


def complement_nucs(nuc_str: str):
    replacements = []
    comp_str = ""
    trans_map = {"A": "T", "T": "A", "C": "G", "G": "C", "U": "A", 'N': 'N'}
    for c in nuc_str.upper():
        try:
            comp_str += trans_map[c]
        except KeyError:
            if c == '.' or c == '-':
                comp_str += c
            else:
                replacements.append(c)
                comp_str += 'N'

    if replacements:
        logging.warning(str(len(replacements)) +
                        " ambiguity character(s) (" + ', '.join(sorted(set(replacements))) +
                        ") replaced by 'N' while complementing\n")
    return comp_str


def reverse_complement(nuc_sequence: str):
    return complement_nucs(nuc_sequence[::-1])


def fish_refpkg_from_build_params(bait: str, marker_build_dict: dict):
    """
    Using a marker gene name (first column in ref_build_parameters.tsv, e.g. RecA, McrA) as bait,
    find and return the reference package for that marker gene name by matching the 'cog' elements.

    :param bait: A marker gene name
    :param marker_build_dict: A dictionary of refpkg name keys mapping to MarkerBuild instances
    :return: MarkerBuild object
    """
    refpkg = None
    for denominator in marker_build_dict:
        if bait == marker_build_dict[denominator].cog:
            refpkg = marker_build_dict[denominator]
            break
    if refpkg is None:
        logging.error("Unable to find '" + bait + "' in marker_build_dict!\n")
        sys.exit(13)
    else:
        return refpkg


def find_msa_type(msa_files):
    file_types = set()
    for mc in msa_files:
        sample_msa_file = msa_files[mc][0]
        f_ext = sample_msa_file.split('.')[-1]
        if re.match("phy|phylip", f_ext):
            file_types.add("Phylip")
        elif re.match("sto|stockholm", f_ext):
            file_types.add("Stockholm")
        elif re.match("mfa|fa|fasta", f_ext):
            file_types.add("Fasta")
        else:
            logging.error("Unrecognized file extension: '" + f_ext + "'")
            sys.exit(3)
    if len(file_types) > 1:
        logging.error(
            "Multiple file types detected in multiple alignment files:\n" + ','.join(file_types) + "\n")
        sys.exit(3)
    elif len(file_types) == 0:
        logging.error("No alignment files were generated!\n")
        sys.exit(3)
    else:
        return file_types.pop()


def swap_tree_names(tree_to_swap, final_mltree, code_name):
    original_tree = open(tree_to_swap, 'r')
    raxml_tree = open(final_mltree, 'w')

    tree = original_tree.readlines()
    original_tree.close()
    if len(tree) > 1:
        logging.error(">1 line contained in RAxML tree " + tree_to_swap + "\n")
        sys.exit(13)

    new_tree = re.sub('_' + re.escape(code_name), '', str(tree[0]))
    raxml_tree.write(new_tree)

    raxml_tree.close()
    return


def match_target_marker(refpkg_name: str, headers: list) -> list:
    """
    Returns the list of sequences with TreeSAPP classification tags matching the refpkg_name
    :param refpkg_name: The refpkg name (e.g. McrA) not code (e.g. M0701) for desired classified sequences
    :param headers: List of classified sequences
    :return: List of headers that match the refpkg_name
    """
    matches = list()
    classified_target = re.compile(r'\|{0}\|\d+_\d+$'.format(re.escape(refpkg_name)))
    for seq_name in headers:
        if classified_target.search(seq_name):
            matches.append(seq_name)
    return matches


def get_hmm_length(hmm_file):
    """
    Function to open the ref_tree's hmm file and determine its length
    :param hmm_file: The HMM file produced by hmmbuild to parse for the HMM length
    :return: The length (int value) of the HMM
    """
    try:
        hmm = open(hmm_file, 'r')
    except IOError:
        raise IOError("Unable to open " + hmm_file + " for reading! Exiting.")

    line = hmm.readline()
    length = 0
    while line:
        # LENG XXX
        if re.match(r"^LENG\s+([0-9]+)", line):
            length = int(line.split()[1])
        line = hmm.readline()
    if length > 0:
        return length
    else:
        raise AssertionError("Unable to parse the HMM length from " + hmm_file + ". Exiting.")


def write_dict_to_table(data_dict: dict, output_file: str, sep="\t") -> None:
    """
    Function for writing a dictionary of key: value pairs separated to a file.
    :param data_dict: Dictionary containing keys (e.g. accessions) and values (e.g. lineages)
    :param output_file: Path to a file to write to
    :param sep: Separator to use between the keys and values. Default is tab.
    :return: None
    """
    data_strings = []
    for key in data_dict:
        values = data_dict[key]
        if type(values) is str:
            data_strings.append(sep.join([key, values]))
        elif type(values) is list:
            data_strings.append(sep.join([key, sep.join(values)]))
        else:
            logging.error("Unable to tabularize values of type '" + str(type(values)) + "'\n")
            sys.exit(5)
    try:
        handler = open(output_file, 'w')
    except IOError:
        logging.error("Unable to open file '" + output_file + "' for writing.\n")
        sys.exit(3)
    handler.write("\n".join(data_strings) + "\n")
    handler.close()

    return
