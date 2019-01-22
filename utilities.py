__author__ = 'Connor Morgan-Lang'

import os
import re
import sys
import subprocess
import logging
from external_command_interface import launch_write_command


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


def find_executables(args):
    """
    Finds the executables in a user's path to alleviate the requirement of a sub_binaries directory
    :param args: command-line arguments objects
    :return: exec_paths beings the absolute path to each executable
    """
    exec_paths = dict()
    dependencies = ["prodigal", "hmmbuild", "hmmalign", "hmmsearch", "raxmlHPC", "usearch", "trimal", "BMGE.jar", "papara"]

    # Extra executables necessary for certain modes of TreeSAPP
    if hasattr(args, "rpkm") and args.rpkm:
        dependencies += ["bwa", "rpkm"]

    if hasattr(args, "update_tree"):
        if args.update_tree:
            dependencies += ["usearch", "blastn", "blastp", "makeblastdb", "mafft"]

    if hasattr(args, "cluster") or hasattr(args, "multiple_alignment") or hasattr(args, "fast"):
        dependencies += ["mafft", "OD-seq"]
        if args.cluster:
            dependencies.append("usearch")
        if args.fast:
            dependencies.append("FastTree")

    if args.molecule == "rrna":
        dependencies += ["cmalign", "cmsearch", "cmbuild"]

    if os_type() == "linux":
        args.executables = args.treesapp + "sub_binaries" + os.sep + "ubuntu"
    if os_type() == "mac":
        args.executables = args.treesapp + "sub_binaries" + os.sep + "mac"
    elif os_type() == "win" or os_type() is None:
        logging.error("Unsupported OS")
        sys.exit(13)

    for dep in dependencies:
        if is_exe(args.executables + os.sep + dep):
            exec_paths[dep] = str(args.executables + os.sep + dep)
        # For rpkm and potentially other executables that are compiled ad hoc
        elif is_exe(args.treesapp + "sub_binaries" + os.sep + dep):
            exec_paths[dep] = str(args.treesapp + "sub_binaries" + os.sep + dep)
        elif which(dep):
            exec_paths[dep] = which(dep)
        else:
            logging.error("Could not find a valid executable for " + dep + ".\n")
            sys.exit(13)

    args.executables = exec_paths
    return args


def reformat_string(string):
    if string and string[0] == '>':
        header = True
    else:
        header = False
    string = re.sub("\[|\]|\(|\)|\/|\\\\|'|<|>", '', string)
    if header:
        string = '>' + string
    string = re.sub("\s|;|,|\|", '_', string)
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


def return_sequence_info_groups(regex_match_groups, header_db, header):
    description = ""
    locus = ""
    organism = ""
    lineage = ""
    if regex_match_groups:
        accession = regex_match_groups.group(1)
        if header_db == "custom":
            lineage = regex_match_groups.group(2)
            organism = regex_match_groups.group(3)
            description = regex_match_groups.group(3)
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

    if not accession and not organism:
        logging.error("Insufficient information was loaded for header:\n" +
                      header + "\n" + "regex_match: " + header_db + '\n')
        sys.exit(13)

    return accession, organism, locus, description, lineage


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


def generate_blast_database(args, fasta, molecule, prefix, multiple=True):
    """

    :param args:
    :param fasta: File to make a BLAST database for
    :param molecule: 'prot' or 'nucl' - necessary argument for makeblastdb
    :param prefix: prefix string for the output BLAST database
    :param multiple: Flag indicating the input `fasta` is a MSA. Alignment information is removed prior to makeblastdb
    :return:
    """

    # Remove the multiple alignment information from fasta_replaced_file and write to fasta_mltree
    blastdb_out = prefix + ".fa"
    if multiple:
        if blastdb_out == fasta:
            logging.error("prefix.fa is the same as " + fasta + " and would be overwritten!\n")
            sys.exit(13)
        remove_dashes_from_msa(fasta, blastdb_out)
        blastdb_in = blastdb_out
    else:
        blastdb_in = fasta

    logging.info("Making the BLAST database for " + blastdb_in + "... ")

    # Format the `makeblastdb` command
    makeblastdb_command = [args.executables["makeblastdb"]]
    makeblastdb_command += ["-in", blastdb_in]
    makeblastdb_command += ["-out", blastdb_out]
    makeblastdb_command += ["-input_type", "fasta"]
    makeblastdb_command += ["-dbtype", molecule]

    # Launch the command
    stdout, makeblastdb_pro_returncode = launch_write_command(makeblastdb_command)

    logging.info("done\n")

    return stdout, blastdb_out


def clean_lineage_string(lineage):
    non_standard_names_re = re.compile(" group| cluster| complex", re.IGNORECASE)
    bad_strings = ["cellular organisms; ", "delta/epsilon subdivisions; ", "\(miscellaneous\)", "Root; ", "[a-p]__"]
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
    return lineage


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


def convert_outer_to_inner_nodes(clusters, internal_node_map):
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
                    leaf_annotation_map[cluster].append(inode)
                    break
    return leaf_annotation_map


def annotate_internal_nodes(args, internal_node_map, clusters):
    """
    A function for mapping the clusters to all internal nodes of the tree.
    It also adds overlapping functional annotations for deep internal nodes and ensures all the leaves are annotated.
    :param args:
    :param internal_node_map: A dictionary mapping the internal nodes (keys) to the leaf nodes (values)
    :param clusters: Dictionary with the cluster names for keys and a tuple containing leaf boundaries as values
    :return: A dictionary of the annotation (AKA group) as keys and internal nodes as values
    """
    annotated_clade_members = dict()
    leaf_group_members = dict()
    leaves_in_clusters = set()

    # Create a dictionary to map the cluster name (e.g. Function, Activity, Class, etc) to all the leaf nodes
    for annotation in clusters.keys():
        if annotation not in annotated_clade_members:
            annotated_clade_members[annotation] = set()
        if annotation not in leaf_group_members:
            leaf_group_members[annotation] = set()
        for i_node in internal_node_map:
            if i_node in clusters[annotation]:
                for leaf in internal_node_map[i_node]:
                    leaf_group_members[annotation].add(leaf)
                    leaves_in_clusters.add(leaf)
        # Find the set of internal nodes that are children of this annotated clade
        for i_node in internal_node_map:
            if leaf_group_members[annotation].issuperset(internal_node_map[i_node]):
                annotated_clade_members[annotation].add(i_node)

    logging.debug("\tCaptured " + str(len(leaves_in_clusters)) + " nodes in clusters.\n")

    return annotated_clade_members, leaves_in_clusters


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


def calculate_overlap(info):
    """
    Returns the overlap length of the base and the check sequences.
    :param info: Autovivify() object holding start and end sequence coordinates for overlapping sequences
    :return overlap: The number of overlapping bases between the sequences
    """

    overlap = 0
    base_start = info['base']['start']
    base_end = info['base']['end']
    check_start = info['check']['start']
    check_end = info['check']['end']

    # Calculate the overlap based on the relative positioning of the base and check sequences
    assert isinstance(base_end, (int, int, float, complex))
    if base_start <= check_start:
        if check_end >= base_end >= check_start:
            # Base     ----
            # Check      -------
            overlap = base_end - check_start
        elif check_end <= base_end:
            # Base     --------
            # Check        --
            overlap = check_end - check_start
    elif check_start <= base_start:
        if base_start <= check_end <= base_end:
            # Base         -----
            # Check    -----
            overlap = check_end - base_start
        elif base_end <= check_end:
            # Base       --
            # Check    --------
            overlap = base_end - base_start

    return overlap


def cluster_sequences(uclust_exe, fasta_input, uclust_prefix, similarity=0.60):
    """
    Wrapper function for clustering a FASTA file at some similarity using usearch's cluster_fast algorithm

    :param uclust_exe: Path to the usearch executable
    :param fasta_input: FASTA file for which contained sequences will be clustered
    :param uclust_prefix: Prefix for the output files
    :param similarity: The proportional similarity to cluster input sequences
    :return: None
    """
    logging.info("Clustering sequences with UCLUST... ")
    uclust_cmd = [uclust_exe]
    uclust_cmd += ["-cluster_fast", fasta_input]
    uclust_cmd += ["-id", str(similarity)]
    uclust_cmd += ["-sort", "length"]
    uclust_cmd += ["-centroids", uclust_prefix + ".fa"]
    uclust_cmd += ["--uc", uclust_prefix + ".uc"]
    logging.info("done.\n")
    stdout, returncode = launch_write_command(uclust_cmd)

    if returncode != 0:
        logging.error("UCLUST did not complete successfully! Command used:\n" +
                      ' '.join(uclust_cmd) + "\n")
        sys.exit(13)
    return


def profile_aligner(executables, ref_aln, ref_profile, input_fasta, output_multiple_alignment, kind="functional"):
    """
    A wrapper for both cmalign and hmmalign for performing profile-based multiple sequence alignment
    :param executables: A dictionary containing keys "cmalign" and "hmmalign"
    :param ref_aln: Path to a FASTA or Stockholm file with the multiple alignment pattern
    :param ref_profile: Path to the HMM or CM profile for the reference gene
    :param input_fasta: Path to the FASTA containing query sequences
    :param output_multiple_alignment: Name of the output Stockholm formatted file
    :param kind: The type of marker gene being analyzed [functional (default), phylogenetic, phylogenetic_rRNA]
    :return:
    """

    if kind == "phylogenetic_rRNA":
        malign_command = [executables["cmalign"]]
    else:
        malign_command = [executables["hmmalign"]]

    malign_command += ['--mapali', ref_aln,
                       '--outformat', 'Stockholm',
                       ref_profile, input_fasta,
                       '>', output_multiple_alignment]

    stdout, returncode = launch_write_command(malign_command)
    if returncode != 0:
        logging.error("Multiple alignment failed for " + input_fasta + ". Command used:\n" +
                      ' '.join(malign_command) + " output:\n" + stdout + "\n")
        sys.exit(3)
    return stdout


def run_papara(executable, tree_file, ref_alignment_phy, query_fasta, molecule):
    papara_command = [executable]
    papara_command += ["-t", tree_file]
    papara_command += ["-s", ref_alignment_phy]
    papara_command += ["-q", query_fasta]
    if molecule == "prot":
        papara_command.append("-a")

    stdout, ret_code = launch_write_command(papara_command)
    if ret_code != 0:
        logging.error("PaPaRa did not complete successfully!\n" +
                      "Command used:\n" + ' '.join(papara_command) + "\n")
        sys.exit(3)
    return stdout


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
