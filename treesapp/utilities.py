import os
import re
import sys
import logging
import shutil
from glob import glob
from csv import Sniffer

from pygtrie import StringTrie
import multiprocessing

from treesapp import external_command_interface as eci
from treesapp import logger

LOGGER = logging.getLogger(logger.logger_name())


def base_file_prefix(file_path: str) -> str:
    return os.path.splitext(os.path.basename(file_path))[0]


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
        LOGGER.error("Key map (" + str(len(key_map)) + ") and original dictionary (" + str(len(og_dict)) +
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
        LOGGER.warning("Dictionary rekey incomplete - " + str(len(unmapped)) + " keys were not found in dictionary.\n")

    if len(og_dict) > 0:
        LOGGER.warning(str(len(og_dict)) + " keys in original dictionary were not popped during re-keying process.\n")

    return updated_dict


def is_exe(fpath):
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)


def which(program: str):
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


def globify_path(path_name: str) -> str:
    """Ensures the path name is compatible with glob.glob()."""
    # Properly escape square brackets for glob.glob()
    i = 0
    while i < 10 and re.search(r'\w\[\w', path_name):
        re_m = re.search(r'\w(\[)\w', path_name)
        s, e = re_m.regs[1]
        path_name = path_name[0:s] + '[[]' + path_name[e:]
        i += 1
    return path_name


def match_file(glob_pattern: str) -> str:
    """
    Using a valid glob pattern, glob.glob is used to find a single file and
    return the path to the file that matches the pattern.

    :param glob_pattern: A string representing a glob pattern. Used to search for files.
    :return: Path to the single file matching the glob pattern
    """
    glob_pattern = globify_path(glob_pattern)

    file_matches = glob(glob_pattern)

    if len(file_matches) > 1:
        LOGGER.error("Multiple files match glob pattern '{}':\n{}".format(glob_pattern, ", ".join(file_matches)))
        sys.exit(17)
    elif len(file_matches) == 0:
        LOGGER.error("Unable to find file matching glob pattern '{}'.\n".format(glob_pattern))
        sys.exit(19)
    else:
        return file_matches.pop()


def available_cpu_count():
    """ Number of available virtual or physical CPUs on this system, i.e.
    user/real as output by time(1) when called with an optimally scaling
    userspace-only program"""

    try:
        res = multiprocessing.cpu_count()
        return res
    except (ImportError, NotImplementedError):
        pass

    try:
        m = re.search(r'(?m)^Cpus_allowed:\s*(.*)$',
                      open('/proc/self/status').read())
        if m:
            res = bin(int(m.group(1).replace(',', ''), 16)).count('1')
            if res > 0:
                return res
    except IOError:
        pass

    try:
        return os.cpu_count()
    except (ImportError, NotImplementedError):
        pass

    LOGGER.error('Can not determine number of CPUs on this system')
    return


def executable_dependency_versions(exe_dict: dict) -> str:
    """
    Function for retrieving the version numbers for each executable in exe_dict

    :param exe_dict: A dictionary mapping names of software to the path to their executable
    :return: A formatted string with the executable name and its respective version found
    """
    versions_dict = dict()
    versions_string = "Software versions used:\n"

    simple_v = ["prodigal", "raxmlHPC", "epa-ng"]
    version_param = ["trimal", "mafft", "raxml-ng"]
    no_params = ["vsearch", "papara"]
    help_param = ["hmmbuild", "hmmalign", "hmmsearch", "OD-seq"]
    version_re = re.compile(r"[Vv]\d+.\d|version \d+.\d|\d\.\d\.\d|HMMER")

    for exe in exe_dict:
        ##
        # Get the help/version statement for the software
        ##
        versions_dict[exe] = ""
        if exe in simple_v:
            stdout, returncode = eci.launch_write_command([exe_dict[exe], "-v"])
        elif exe in version_param:
            stdout, returncode = eci.launch_write_command([exe_dict[exe], "--version"])
        elif exe in help_param:
            stdout, returncode = eci.launch_write_command([exe_dict[exe], "-h"])
        elif exe in no_params:
            stdout, returncode = eci.launch_write_command([exe_dict[exe]])
        elif exe == "mmseqs":
            stdout, returncode = eci.launch_write_command([exe_dict[exe], "version"])
            versions_dict[exe] = stdout.strip()
        elif exe == "FastTree":
            stdout, returncode = eci.launch_write_command([exe_dict[exe], "-expert"])
        elif exe == "BMGE.jar":
            stdout, returncode = eci.launch_write_command(["java", "-Xmx10m", "-jar", exe_dict[exe], "-?"])
        else:
            LOGGER.warning("Unknown version command for " + exe + ".\n")
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
            LOGGER.debug("Unable to find version for " + exe + ".\n")

    ##
    # Format the string with the versions of all software
    ##
    for exe in sorted(versions_dict):
        n_spaces = 12-len(exe)
        versions_string += "\t" + exe + ' '*n_spaces + versions_dict[exe] + "\n"

    return versions_string


def median(num_list: list):
    n = len(num_list)
    if n < 1:
        return None
    elif n == 1:
        return num_list[0]
    elif n % 2 == 1:
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


def reformat_fasta_to_phy(fasta_dict: dict) -> dict:
    """
    The fasta_dict input is a dictionary of sequence names (seq_name) keys indexing their respective sequences.
    Each sequence in the dictionary is split into subsequences of 50 characters and indexed first by the line count
    as well as the sequence name. These are stored in a dictionary of dictionaries and fairly trivial to iterate the
    sorted keys.

    :param fasta_dict: A dictionary of sequence names (seq_name) keys indexing their respective sequences.
    :return: A dictionary of line counts indexing sequence names indexing that sequence's subsequence for that line
    """
    phy_dict = dict()
    for seq_name in fasta_dict:
        sequence = fasta_dict[seq_name]
        sub_sequences = re.findall(r'.{1,50}', sequence)
        count = 0
        for sub_sequence in sub_sequences:
            try:
                phy_dict[count][seq_name] = sub_sequence
            except KeyError:
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
            LOGGER.error("Inconsistent number of sequences in Phylip dictionary keys.")
        num_seqs = seq_chunks[0]
        aligned_seqs = dict()
        for seq_chunk in phy_dict.keys():
            for seq_name in phy_dict[seq_chunk].keys():
                if seq_name not in aligned_seqs:
                    aligned_seqs[seq_name] = ""
                aligned_seqs[seq_name] += phy_dict[seq_chunk][seq_name]
        aligned_seq_lengths = [len(aligned_seqs[aligned_seq]) for aligned_seq in aligned_seqs.keys()]
        if min(aligned_seq_lengths) != max(aligned_seq_lengths):
            LOGGER.error("Lengths of aligned sequences are heterogeneous.")
        alignment_len = aligned_seq_lengths[0]
        aligned_seqs.clear()
        seq_chunks.clear()
    else:
        num_seqs, alignment_len = alignment_dims

    longest_seq_name = max([len(seq_name) for seq_name in phy_dict[0]])
    with open(phy_output_file, 'w') as phy_output:
        phy_string = ' ' + str(num_seqs) + ' ' + str(alignment_len) + '\n'
        for count in sorted(phy_dict.keys(), key=int):
            for seq_name in sorted(phy_dict[count].keys()):
                sequence_part = re.sub('X', '-', phy_dict[count][seq_name])
                if count == 0:
                    phy_string += str(seq_name)
                    length = len(str(seq_name))
                    c = length
                    while c < longest_seq_name + 1:
                        phy_string += ' '
                        c += 1
                else:
                    phy_string += (longest_seq_name+1)*' '
                phy_string += ' '.join(re.findall(r'.{1,10}', sequence_part)) + '\n'
            phy_string += "\n"

        phy_output.write(phy_string)
    return


def extract_hmm_matches(hmm_matches: dict, fasta_dict: dict, header_registry: dict) -> dict:
    """
    Function for slicing sequences guided by alignment co-ordinates.

    :param hmm_matches: A dictionary containing a list of HmmMatch() objects indexed by reference package name
    :param fasta_dict: A dictionary with headers as keys and sequences as values
    :param header_registry: A list of Header() objects, each used to map various header formats to each other
    :return: Dictionary containing a modified query name mapped to its corresponding (sub)sequence that was aligned
     to the HMM, indexed by the HMM profile name (i.e. marker or reference package name)
    """

    marker_gene_dict = dict()
    header_matching_dict = dict()

    LOGGER.debug("Creating a temporary dictionary for rapid sequence name look-ups... ")
    for num in header_registry:
        if header_registry[num].first_split[0] == '>':
            header_matching_dict[header_registry[num].first_split[1:]] = header_registry[num]
        else:
            header_matching_dict[header_registry[num].first_split] = header_registry[num]
    LOGGER.debug("done.\n")

    LOGGER.info("Extracting the quality-controlled protein sequences... ")

    for refpkg_name in hmm_matches:
        extracted_loci = dict()
        for hmm_match in hmm_matches[refpkg_name]:
            # Now for the header format to be used in the bulk FASTA:
            # >contig_name|marker_gene|start_end
            q_header = header_matching_dict[hmm_match.orf]
            if hmm_match.of > 1:
                q_header.post_align = ' '.join([q_header.first_split,
                                                str(hmm_match.num) + '.' + str(hmm_match.of),
                                                re.sub(re.escape(q_header.first_split), '', q_header.original)]).strip()
            else:
                q_header.post_align = q_header.original

            if q_header.post_align in extracted_loci:
                LOGGER.warning("Query '{}' being overwritten by an alternative alignment:\n"
                               "{}\n".format(q_header.post_align, hmm_match.get_info()))
            try:
                extracted_loci[q_header.post_align] = fasta_dict[q_header.original][hmm_match.start-1:hmm_match.end]
            except KeyError:
                LOGGER.debug("Unable to map '{}' to a sequence in the input FASTA.\n".format(hmm_match.orf))

        marker_gene_dict[refpkg_name] = extracted_loci

    LOGGER.info("done.\n")
    return marker_gene_dict


def hmm_pile(hmm_matches: dict) -> None:
    """
    Function to inspect the placement of query sequences on the reference HMM

    :param hmm_matches: A dictionary of HmmMatch instances indexed by the HMM profile name they mapped to
    :return: None
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
            LOGGER.info("Low coverage " + marker + " profile windows (start-stop):\n" + low_cov_summary)
        LOGGER.info("Maximum coverage for " + marker + " = " + str(maximum_coverage) + " sequences\n")
    return


def concatenate_files(input_files: list, output_path: str) -> None:
    """A Pythonic replacement for UNIX's cat utility."""
    try:
        output_handler = open(output_path, 'wb')
    except IOError:
        LOGGER.error("Unable to open output file '{}' for writing.\n".format(output_path))
        sys.exit(17)

    if isinstance(input_files, str):
        raise AssertionError("File list is str, when list was expected.")
    for input_f in input_files:
        with open(input_f, 'rb') as fh:
            shutil.copyfileobj(fh, output_handler)

    output_handler.close()
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
        LOGGER.warning("{} ambiguity character(s) ({}) replaced by 'N' while complementing\n"
                       "".format(len(replacements), ', '.join(sorted(set(replacements)))))
    return comp_str


def reverse_complement(nuc_sequence: str):
    return complement_nucs(nuc_sequence[::-1])


def find_msa_type(msa_files: dict):
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
            LOGGER.error("Unrecognized file extension: '" + f_ext + "'")
            sys.exit(3)
    if len(file_types) > 1:
        LOGGER.error(
            "Multiple file types detected in multiple alignment files:\n" + ','.join(file_types) + "\n")
        sys.exit(3)
    elif len(file_types) == 0:
        LOGGER.error("No alignment files were generated!\n")
        sys.exit(3)
    else:
        return file_types.pop()


def get_file_lines(file_path: str, re_pattern=None, max_matches=0) -> list:
    """
    Opens a file and returns it's contents as a list using readlines().
    Optionally a pattern can be provided and used for filtering the lines of the file.
     Only lines matching this pattern will be returned

    :param file_path: Path to a file
    :param re_pattern: An optional compiled re pattern to filter the lines by
    :param max_matches: An optional limit to the number of lines that match the re_pattern
    :return: The lines of a file as a list
    """
    try:
        file_handler = open(file_path, 'r')
    except IOError:
        raise IOError("Unable to open file '{}' for reading! Exiting.".format(file_path))

    lines = []
    n_matches = 0
    if not re_pattern:
        lines += file_handler.readlines()
    else:
        line = file_handler.readline()
        while line:
            if re_pattern.match(line):
                lines.append(line)
                n_matches += 1
            if 0 < max_matches < n_matches:
                break
            line = file_handler.readline()
    file_handler.close()

    return lines


def get_hmm_value(profile_hmm, attribute_name: str) -> str:
    """
    Function to parse an attribute from a profile HMM file made by HMMER

    :param profile_hmm: Either a list of strings from the profile HMM file or the path to the file produced by hmmbuild
    :param attribute_name: The name of the attribute to parse from the profile HMM file.
     Examples are 'length', 'name', 'alphabet', and 'num_seq'.
    :return: The profile HMM's value for the requested attribute
    """
    attribute_name_map = {"length": "LENG", "name": "NAME", "alphabet": "ALPH", "num_seqs": "NSEQ"}
    try:
        hmm_attr_name = attribute_name_map.get(attribute_name)
    except KeyError:
        LOGGER.error("Unsupported profile HMM feature '{}'".format(attribute_name))
        raise KeyError()

    # Compiled regular expression for the requested attribute
    attr_re = re.compile(r"^" + re.escape(hmm_attr_name) + r"\s+(\w+)")

    # Gather the lines that match the attribute name
    attribute_lines = []
    if isinstance(profile_hmm, str):
        attribute_lines = get_file_lines(profile_hmm, re_pattern=attr_re, max_matches=1)
    else:
        for line in profile_hmm:
            if attr_re.match(line):
                attribute_lines.append(line)

    if not attribute_lines:
        LOGGER.error("Unable to parse '{}' from the profile HMM.\n".format(attribute_name))
        raise AssertionError()
    elif len(attribute_lines) > 1:
        LOGGER.error("More than one profile HMM included in file.\n")
        raise AssertionError()
    else:
        attr_line = attribute_lines.pop()
        return ' '.join(attr_line.split()[1:])


def get_hmm_length(hmm_file: str) -> int:
    return int(get_hmm_value(hmm_file, "length"))


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
        if isinstance(values, str):
            data_strings.append(sep.join([key, values]))
        elif isinstance(values, list):
            data_strings.append(sep.join([key, sep.join(values)]))
        else:
            LOGGER.error("Unable to tabularize values of type '{}'\n".format(type(values)))
            sys.exit(5)
    try:
        handler = open(output_file, 'w')
    except IOError:
        LOGGER.error("Unable to open file '" + output_file + "' for writing.\n")
        sys.exit(3)
    handler.write("\n".join(data_strings) + "\n")
    handler.close()

    return


def dict_diff(snap_one: dict, snap_two: dict) -> list:
    """
    Takes two dictionaries and based on the difference between the two keys,
    a list is returned of values of keys that are only present in the first dictionary.

    :param snap_one: A dictionary that is a superset of snap_two
    :param snap_two: A dictionary that is a subset of snap_one
    :return: A list of values that are only present in the first dictionary
    """
    one_keys = set(snap_one.keys())
    two_keys = set(snap_two.keys())
    diff_indexes = one_keys.difference(two_keys)

    return [snap_one[index] for index in diff_indexes]


def get_list_positions(ref_list: list, queries: list, all_queries=False) -> dict:
    """
    Used for finding the position of strings (queries) in a list (ref_list)

    :param ref_list: A list of strings
    :param queries: A list of strings
    :param all_queries: If True, all queries will be returned in the dictionary even if they are not in ref_list
    :return: Dictionary of strings in queries as keys and their position in ref_list as values
    """
    # Set the default value if all queries are to be returned
    if all_queries:
        positions = {query: None for query in queries}
    else:
        positions = {}

    # Ensure the references and queries are both lowercase and can be compared
    ref_list = [i.strip().lower() for i in ref_list]
    for query in queries:  # type: str
        try:
            positions[query] = ref_list.index(query.lower())
        except ValueError:
            continue

    return positions


def get_field_delimiter(file_path: str, sniff_size=50) -> str:
    with open(file_path, newline='') as csvfile:
        sniffer = Sniffer()
        sniffer.preferred = ["\t", ","]

        # Ensure there are at least sniff_size lines in the file
        x = 1
        line = csvfile.readline()
        sample = ""
        while line and x <= sniff_size:
            sample += line
            x += 1
            line = csvfile.readline()

        # Return to the beginning
        csvfile.seek(0)
        if x < sniff_size:
            LOGGER.info("{} contains {} lines which were used to determine field delimiter.\n".format(file_path, x))

        # Read the file to determine its delimiter
        dialect = sniffer.sniff(sample, delimiters=',;\t')
    return str(dialect.delimiter)


def validate_new_dir(output_dir: str) -> str:
    output_dir = os.path.abspath(output_dir)
    # Check whether the output path exists
    up, down = os.path.split(output_dir.rstrip(os.sep))
    if output_dir[-1] != os.sep:
        output_dir += os.sep
    if not os.path.isdir(up):
        LOGGER.error("The directory above output ({}) does not exist.\n"
                      "Please make these as TreeSAPP only creates a single new directory.".format(up))
        sys.exit(3)
    return output_dir


def fetch_executable_path(exe_name, treesapp_dir):
    if is_exe(os.path.join(treesapp_dir, "sub_binaries", exe_name)):
        return str(os.path.join(treesapp_dir, "sub_binaries", exe_name))
    elif which(exe_name):
        return which(exe_name)
    else:
        LOGGER.error("Could not find a valid executable for '{}'.\n".format(exe_name))
        sys.exit(13)


def elegant_pair(x: int, y: int, sort=False) -> int:
    if sort:
        x, y = sorted([x, y])
    if y > x:
        return y*y + x
    else:
        return x*x + x + y
