__author__ = 'Connor Morgan-Lang'

import os
import re
import sys
import subprocess

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

    raise Exception('Can not determine number of CPUs on this system')


def find_executables(args):
    """
    Finds the executables in a user's path to alleviate the requirement of a sub_binaries directory
    :param args: command-line arguments objects
    :return: exec_paths beings the absolute path to each executable
    """
    exec_paths = dict()
    dependencies = ["prodigal", "hmmbuild", "hmmalign", "hmmsearch", "raxmlHPC", "trimal", "BMGE.jar", "papara"]
    # old_dependencies = ["blastn", "blastx", "blastp", "genewise", "Gblocks", "makeblastdb", "muscle"]

    # Extra executables necessary for certain modes of TreeSAPP
    if hasattr(args, "rpkm") and args.rpkm:
        dependencies += ["bwa", "rpkm"]

    if hasattr(args, "update_tree"):
        if args.update_tree:
            dependencies += ["usearch", "blastn", "blastp", "makeblastdb", "mafft"]

    if hasattr(args, "cluster") or hasattr(args, "multiple_alignment") or hasattr(args, "fast"):
        if args.cluster:
            dependencies.append("usearch")
        dependencies.append("mafft")
        if args.fast:
            dependencies.append("FastTree")

    if args.molecule == "rrna":
        dependencies += ["cmalign", "cmsearch", "cmbuild"]

    if os_type() == "linux":
        args.executables = args.treesapp + "sub_binaries" + os.sep + "ubuntu"
    if os_type() == "mac":
        args.executables = args.treesapp + "sub_binaries" + os.sep + "mac"
    elif os_type() == "win" or os_type() is None:
        sys.exit("ERROR: Unsupported OS")

    for dep in dependencies:
        if is_exe(args.executables + os.sep + dep):
            exec_paths[dep] = str(args.executables + os.sep + dep)
        # For rpkm and potentially other executables that are compiled ad hoc
        elif is_exe(args.treesapp + "sub_binaries" + os.sep + dep):
            exec_paths[dep] = str(args.treesapp + "sub_binaries" + os.sep + dep)
        elif which(dep):
            exec_paths[dep] = which(dep)
        else:
            sys.stderr.write("Could not find a valid executable for " + dep + ". ")
            sys.exit("Bailing out.")

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
    accession = ""
    description = ""
    locus = ""
    organism = ""
    lineage = ""
    if regex_match_groups:
        if len(regex_match_groups.groups()) == 2:
            accession = regex_match_groups.group(1)
            organism = regex_match_groups.group(2)
            description = regex_match_groups.group(2)
        elif header_db in ["ncbi_ambig", "refseq_prot", "gen_genome"]:
            accession = regex_match_groups.group(1)
            description = regex_match_groups.group(2)
            organism = regex_match_groups.group(3)
        elif header_db == "silva":
            accession = regex_match_groups.group(1)
            locus = str(regex_match_groups.group(2)) + '-' + str(regex_match_groups.group(3))
            lineage = regex_match_groups.group(4)
            description = regex_match_groups.group(4)
        elif header_db == "fungene":
            accession = regex_match_groups.group(1)
            locus = regex_match_groups.group(2)
            organism = regex_match_groups.group(3)
            description = regex_match_groups.group(3)
        elif header_db == "fungene_truncated":
            accession = regex_match_groups.group(1)
            organism = regex_match_groups.group(2)
            description = regex_match_groups.group(3)
        elif header_db == "custom":
            description = regex_match_groups.group(1)
            lineage = regex_match_groups.group(2)
            organism = regex_match_groups.group(3)
    else:
        sys.stderr.write("Unable to handle header: " + header + "\n")
        sys.exit()

    if not accession and not organism:
        sys.stderr.write("ERROR: Insufficient information was loaded for header:\n" + header + "\n")
        sys.stderr.write("regex_match: " + header_db + '\n')
        sys.exit(33)

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
            sys.stderr.write("ERROR: prefix.fa is the same as " + fasta + " and would be overwritten!\n")
            sys.exit(11)
        remove_dashes_from_msa(fasta, blastdb_out)
        blastdb_in = blastdb_out
    else:
        blastdb_in = fasta

    sys.stdout.write("Making the BLAST database for " + blastdb_in + "... ")

    # Format the `makeblastdb` command
    makeblastdb_command = [args.executables["makeblastdb"]]
    makeblastdb_command += ["-in", blastdb_in]
    makeblastdb_command += ["-out", blastdb_out]
    makeblastdb_command += ["-input_type", "fasta"]
    makeblastdb_command += ["-dbtype", molecule]

    # Launch the command
    stdout, makeblastdb_pro_returncode = launch_write_command(makeblastdb_command)

    sys.stdout.write("done\n")
    sys.stdout.flush()

    return stdout, blastdb_out


def clean_lineage_string(lineage):
    non_standard_names_re = re.compile("group; | cluster; ")
    bad_strings = ["cellular organisms; ", "delta/epsilon subdivisions; ", "\(miscellaneous\)", "Root; ", "[a-p]__"]
    for bs in bad_strings:
        lineage = re.sub(bs, '', lineage)
    # filter 'group' and 'cluster'
    if non_standard_names_re.search(lineage):
        reconstructed_lineage = ""
        ranks = lineage.split("; ")
        for rank in ranks:
            if not (re.search("group$", rank) or re.search("cluster$", rank)):
                reconstructed_lineage = reconstructed_lineage + str(rank) + '; '
        reconstructed_lineage = re.sub('; $', '', reconstructed_lineage)
        lineage = reconstructed_lineage
    return lineage


def median(lst):
    n = len(lst)
    if n < 1:
            return None
    if n % 2 == 1:
            return sorted(lst)[n//2]
    else:
            return sum(sorted(lst)[n//2-1:n//2+1])/2.0


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

    if args.verbose:
        sys.stdout.write("\tCaptured " + str(len(leaves_in_clusters)) + " nodes in clusters.\n")

    return annotated_clade_members, leaves_in_clusters
