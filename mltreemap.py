#!/usr/bin/python

__author__ = "Connor Morgan-Lang, Kishori Konwar and Young Song"
__maintainer__ = "Connor Morgan-Lang"
__license__ = "GPL"
__version__ = "1.1.0"

try:
    import argparse
    import sys
    import os
    from os import path
    from os import listdir
    from os.path import isfile, join
    import errno
    import shutil
    import re
    import glob
    import subprocess
    import signal
    import time
    import traceback
    from multiprocessing import Pool, Process, Lock, Queue, JoinableQueue
    import string
    import random
    from time import gmtime, strftime
except ImportWarning:
    sys.stderr.write("Could not load some user defined module functions")
    sys.stderr.write(traceback.print_exc(10))
    sys.exit(3)


# Classes begin:

class Autovivify(dict):
    """In cases of Autovivify objects, enable the referencing of variables (and sub-variables)
    without explicitly declaring those variables beforehand."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value


class GenewiseWorker(Process):
    """
    A worker that will launch genewise processes in its queue
    """

    def __init__(self, task_queue):
        Process.__init__(self)
        self.task_queue = task_queue

    def run(self):
        while True:
            next_task = self.task_queue.get()
            if next_task is None:
                # Poison pill means shutdown
                self.task_queue.task_done()
                break
            p_genewise = subprocess.Popen(' '.join(next_task), shell=True, preexec_fn=os.setsid)
            p_genewise.wait()
            if p_genewise.returncode != 0:
                sys.stderr.write("ERROR: Genewise did not complete successfully for:\n")
                sys.stderr.write(str(' '.join(next_task)))
                sys.stderr.flush()
            self.task_queue.task_done()
        return


class NodeRetrieverWorker(Process):
    """
    Doug Hellman's Consumer class for handling processes via queues
    """

    def __init__(self, task_queue, result_queue):
        Process.__init__(self)
        self.task_queue = task_queue
        self.result_queue = result_queue

    def run(self):
        while True:
            next_task = self.task_queue.get()
            if next_task is None:
                # Poison pill means shutdown
                self.task_queue.task_done()
                break
            result = get_node_subtrees(next_task, create_tree_info_hash())
            self.task_queue.task_done()
            self.result_queue.put(result)
        return


class CreateFuncTreeUtility:
    """
    Output is the directory to write the outputs for the updated tree
    InputData is the path to the MLTReeMap output folder containing, various_outputs/ and final_RAxML_outputs/
    RefTree is the second column in cog_list.txt for the gene to update
    Cluster is a flag indicating whether the protein sequences for the RefTree in InputData is to be clustered at 97%
    """
    def __init__(self, input_data, ref_tree):
        if os.path.isabs(input_data):
            self.InputData = input_data
        else:
            self.InputData = os.getcwd() + os.sep + input_data

        if self.InputData[-1] == '/':
            self.InputData = self.InputData[:-1]

        self.Output = self.InputData + os.sep + "updated_func_tree" + os.sep
        self.Denominator = ref_tree
        self.COG = ""
        self.ContigDict = dict()
        self.names = list()
        self.ref_names = list()
        # Automatically remove the last attempt at updating the reference tree
        if os.path.isdir(self.Output):
            shutil.rmtree(self.Output)
        try:
            os.makedirs(self.Output)
        except:
            raise IOError("Unable to make the directory " + str(self.Output))

    def get_contigs_for_ref(self):
        """
        Uses self.InputData to find all the RAxML_outputs for each protein sequence for self.RefTree
        :return: list of file names with a protein sequence of self.RefTree
        """
        raxml_outputs = os.listdir(self.InputData + os.sep + "final_RAxML_outputs" + os.sep)
        for placement in raxml_outputs:
            ref_tree = os.path.basename(placement).split('_')[0]
            if ref_tree == self.Denominator:
                suffix = re.sub("%s_" % ref_tree, '', placement)
                predicted_orf = re.sub("_RAxML_parsed.txt", '', suffix)
                self.names.append(predicted_orf)
        return

    def find_cog_name(self, cog_list):
        for cog in cog_list["all_cogs"]:
            denominator = cog_list["all_cogs"][cog]
            if denominator == self.Denominator:
                self.COG = cog
                break
        return

    def write_reference_names(self):
        """
        Generate the mapping between reference taxa IDs and descriptions for the gene being updated
        :return: Dictionary containing alignment_data names as keys and tax_ids descriptions as values
        """
        # TODO: Potentially replace this function with two class variables ref_names and new_names
        ref_tax_id_map = {}

        # Check to see whether we need to use the COG alignment files from GEBA
        geba_ref_match = re.match("g_COG(\d+)", self.COG)

        if geba_ref_match:
            cog_number = geba_ref_match.group(1)
            cog_id = "COG" + cog_number
            ref_alignment_fasta = "data" + os.sep + "geba_alignment_data" + os.sep + cog_id + ".fa"
            ref_alignment_handle = open(ref_alignment_fasta)
        else:
            cog_id = self.COG
            ref_alignment_fasta = "data" + os.sep + "alignment_data" + os.sep + cog_id + ".fa"
            ref_alignment_handle = open(ref_alignment_fasta)
        ref_align_lines = ref_alignment_handle.readlines()

        for line in ref_align_lines:
            line = line.strip()
            header_match = re.match("^>(\d+)_%s" % cog_id, line)

            if header_match:
                header_trimmed = re.sub("^>", "", line)
                ref_tax_id_map[header_trimmed] = ""

        ref_alignment_handle.close()

        # Handle tax ids for COG here #
        cog_input_match = re.match("COG\d+", self.COG)
        geba_cog_match = re.match("g_COG\d+", self.COG)

        if cog_input_match:
            ref_tax_ids_handle = open("data/tree_data/tax_ids_nr.txt", "rb")
        elif geba_cog_match:
            ref_tax_ids_handle = open("data/tree_data/tax_ids_geba_tree.txt", "rb")
        else:
            ref_tax_ids_handle = open("data/tree_data/tax_ids_%s.txt" % self.COG, "rb")

        ref_tax_ids_lines = ref_tax_ids_handle.readlines()

        for each_ref_tax_ids_line in ref_tax_ids_lines:
            each_ref_tax_ids_line = each_ref_tax_ids_line.strip()

            ids_desc_match = re.match("(\S+)\t(\S+(\s+\S+)*)", each_ref_tax_ids_line)

            if ids_desc_match:
                num = ids_desc_match.group(1)
                ref_id = num + "_" + cog_id
                description = ids_desc_match.group(2)

                if ref_id in ref_tax_id_map.keys():
                    ref_tax_id_map[ref_id] = description
                else:
                    AssertionError("Unknown reference number " + str(ref_id) + " in " + ref_alignment_fasta)

        return ref_tax_id_map

    def align_sequences(self, alignment_mode, ref_align, query_fasta, args):
        """
        Call MUSCLE to perform a multiple sequence alignment of the reference sequences and the
        gene sequences identified by MLTreeMap
        :param alignment_mode:
        :param query_fasta: Name of the FASTA file containing the MLTreeMap-identified genes
        :param ref_align: FASTA file containing
        :return: Name of the FASTA file containing the MSA
        """
        if args.verbose:
            sys.stdout.write("Aligning the reference and identified " + self.COG + " sequences using MUSCLE... ")
            sys.stdout.flush()

        # Default alignment #
        if alignment_mode == "d":
            ref_align_gap_removed = self.write_unaligned_ref_fasta(ref_align)

            self.scan_unaligned_ref_fasta(ref_align_gap_removed)
            ref_align_gap_rm_scan = self.Output + self.COG + "_gap_rm_scan.fa"

            concat_fasta = self.Output + self.COG + "_concat.fasta"
            os.system('cat %s %s > %s' % (query_fasta, ref_align_gap_rm_scan, concat_fasta))

            aligned_fasta = self.Output + self.COG + "_d_aligned.fasta"
            muscle_align_command = "muscle -in %s -out %s 1>/dev/null 2>/dev/null" % (concat_fasta, aligned_fasta)

        # Profile-Profile alignment #
        elif alignment_mode == "p":

            query_align = self.Output + self.COG + "_query_aligned.fasta"

            muscle_align_command = "muscle -in %s -out %s 1>/dev/null 2>/dev/null" % (query_fasta, query_align)

            os.system(muscle_align_command)

            aligned_fasta = self.Output + self.COG + "_p_aligned.fasta"
            muscle_align_command = "muscle -profile -in1 %s -in2 %s -out %s 1>/dev/null 2>/dev/null" % \
                                   (query_align, ref_align, aligned_fasta)

        else:
            sys.exit("ERROR: --alignment_mode was not properly assigned!")

        os.system(muscle_align_command)

        if args.verbose:
            sys.stdout.write("done.\n")
            sys.stdout.flush()

        return aligned_fasta

    def write_unaligned_ref_fasta(self, ref_align):
        ref_align_gap_removed = self.Output + self.COG + "_gap_removed.fa"

        ref_align_handle = open(ref_align, "rb")
        ref_align_gap_rm_handle = open(ref_align_gap_removed, "w")

        first_fas_line = ref_align_handle.readline()
        first_fas_line = first_fas_line.strip()

        first_header_match = re.match("^>", first_fas_line)

        if first_header_match:
            ref_align_gap_rm_handle.write(first_fas_line + "\n")

        fasta_in_lines = ref_align_handle.readlines()

        alignment_gap_removed = ""

        for each_fas_line in fasta_in_lines:
            each_fas_line = each_fas_line.strip()

            fasta_header_match = re.match("^>", each_fas_line)

            if fasta_header_match:
                ref_align_gap_rm_handle.write(alignment_gap_removed + "\n")
                ref_align_gap_rm_handle.write(each_fas_line + "\n")

                alignment_gap_removed = ""
            else:

                alignment_gap_removed += each_fas_line

                if re.search("[\-]+", alignment_gap_removed):
                    alignment_gap_removed = re.sub("-", "", alignment_gap_removed)

        ref_align_gap_rm_handle.write(alignment_gap_removed + "\n")

        ref_align_gap_rm_handle.close()
        ref_align_handle.close()
        return ref_align_gap_removed

    def scan_unaligned_ref_fasta(self, ref_align_gap_removed):
        """

        :param ref_align_gap_removed:
        :return:
        """
        ref_align_handle = open(ref_align_gap_removed, "rb")
        ref_align_gap_rm_scan = self.Output + self.COG + "_gap_rm_scan.fa"

        ref_align_scan_handle = open(ref_align_gap_rm_scan, "w")

        fasta_map = {}

        fasta_in_lines = ref_align_handle.readlines()

        sequence_id = ""
        sequence = ""

        for each_fas_line in fasta_in_lines:
            each_fas_line = each_fas_line.strip()

            fasta_header_match = re.match("^>(\S+)", each_fas_line)

            if fasta_header_match:
                sequence_id = each_fas_line
            else:
                sequence = each_fas_line

            fasta_map[sequence_id] = sequence

        for each_sequence_id in fasta_map.keys():
            each_sequence = fasta_map[each_sequence_id]

            line_of_x = re.match("^X((X)+)*$", each_sequence)

            if not line_of_x:
                ref_align_scan_handle.write(each_sequence_id + "\n")
                ref_align_scan_handle.write(each_sequence + "\n")

        ref_align_handle.close()
        ref_align_scan_handle.close()

    def randomize_fasta_id(self, fasta):
        """
        Create a random hash for every reference and new name
        :param fasta: A FASTA file
        :return: Name of fasta_random - the FASTA file with random identifiers
        """
        original_random_dict = {}
        rfive_list = list()

        fasta_handle = open(fasta, "rb")
        fasta_lines = fasta_handle.readlines()

        fasta_random = self.Output + self.COG + "_concat_rfive.fasta"
        fasta_random_handle = open(fasta_random, "w")

        for line in fasta_lines:
            line = line.strip()
            if line[0] == '>':
                original_id = line[1:]
                if original_id not in self.names:
                    self.ref_names.append(original_id)
                rfive_header = "ID"
                rfive = ''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(5))
                while rfive in rfive_list:
                    rfive = ''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(5))
                rfive_list.append(rfive)
                rfive_header += rfive
                original_random_dict[original_id] = rfive_header

                fasta_random_handle.write('>' + rfive_header + "\n")
            else:
                fasta_random_handle.write(line + "\n")

        assert len(set(original_random_dict.values())) == len(original_random_dict.values())

        fasta_random_handle.close()
        fasta_handle.close()
        
        return fasta_random, original_random_dict

    def create_random_names(self, random_map, ref_tax_id_map):
        """
        Write the _concat_rand.names file which contains
        :param random_map: A dictionary with keys from the MSA containing reference and query sequences
        :param ref_tax_id_map: Dictionary with contig names and descriptions for reference sequences
        :return:
        """
        concat_random_names = self.Output + self.COG + "_concat_rand.names"
        concat_rand_names_handle = open(concat_random_names, "w")

        for original_id in random_map.keys():
            names_line = random_map[original_id] + "\t" + ref_tax_id_map[original_id] + "\n"
            concat_rand_names_handle.write(names_line)

        concat_rand_names_handle.close()

    def execute_raxml(self, phylip_file, raxml_destination_folder, args, bootstraps=100):
        os.makedirs(raxml_destination_folder)

        if self.Denominator == "a":
            model_to_be_used = "GTRGAMMA"
        else:
            model_to_be_used = "PROTGAMMAWAG"

        raxml_command = [args.executables["raxmlHPC"], '-m', model_to_be_used]
        # Run RAxML using multiple threads, if CPUs available
        if args.num_threads:
            raxml_command += ['-T', str(int(args.num_threads))]
        else:
            raxml_command += ['-T', '2']
        raxml_command += ['-s', phylip_file,
                          '-f', 'a',
                          '-x', str(12345),
                          '-#', str(bootstraps),
                          '-n', self.COG,
                          '-w', raxml_destination_folder,
                          '-p', str(8),
                          '>', raxml_destination_folder + os.sep + 'RAxML_log.txt']

        raxml_pro = subprocess.Popen(' '.join(raxml_command), shell=True, preexec_fn=os.setsid)
        raxml_pro.wait()

        return
# Classes end


def retrieve_data_size(aligned_fasta):
    # TODO: Replace if possible (just counts number of sequences) and destroy
    num_seqs = 0

    fasta_file_handle = open(aligned_fasta, "rb")

    fasta_lines = fasta_file_handle.readlines()

    for each_fa_line in fasta_lines:
        if re.search(">", each_fa_line):
            num_seqs += 1

    return num_seqs


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


def get_options(): 
    """
    Returns the parser to interpret user options.
    """
    parser = argparse.ArgumentParser(description='Phylogenetically informed insertion of sequence into a reference tree'
                                                 ' using a Maximum Likelihood algorithm.')
    parser.add_argument('-i', '--input', required=True,
                        help='Your sequence input file in FASTA format')
    parser.add_argument('-o', '--output', default='./output/', required=False,
                        help='output directory [DEFAULT = ./output/]')
    parser.add_argument('-b', '--bootstraps', default=0, type=int,
                        help='the number of Bootstrap replicates [DEFAULT = 0]')
    parser.add_argument('-f', '--phylogeny', default='v', choices=['v', 'p'],
                        help='RAxML algorithm (v = Maximum Likelihood [DEFAULT]; p = Maximum Parsimony)')
    parser.add_argument('-g', '--gblocks', default=50, type=int,
                        help='minimal sequence length after Gblocks [DEFAULT = 50]')
    parser.add_argument('-s', '--bitscore', default=60, type=int,
                        help='minimum bitscore for the blast hits [DEFAULT = 60]')
    parser.add_argument('-t', '--reftree', default='p', type=str,
                        help='reference tree (p = MLTreeMap reference tree [DEFAULT]; '
                             'g = GEBA reference tree; i = fungi tree; '
                             'other gene family in data/tree_data/cog_list.txt - e.g., y for mcrA]')
    parser.add_argument('-m', '--molecule', default='n', choices=['a', 'n'],
                        help='the type of input sequences (a = Amino Acid; n = Nucleotide [DEFAULT])')

    rpkm_opts = parser.add_argument_group('RPKM options')
    rpkm_opts.add_argument("--rpkm", action="store_true", default=False,
                           help="Flag indicating RPKM values should be calculated for the gene sequences detected")
    rpkm_opts.add_argument("-r", "--reads", required=False,
                           help="FASTQ file containing to be aligned to predicted genes using BWA MEM")
    rpkm_opts.add_argument("-2", "--reverse", required=False,
                           help="FASTQ file containing to reverse mate-pair reads to be aligned using BWA MEM")
    rpkm_opts.add_argument("-p", "--pairing", required=False, default='pe', choices=['pe', 'se'],
                           help="Indicating whether the reads are paired-end (pe) or single-end (se)")

    update_tree = parser.add_argument_group('Update-tree options')
    # mltreemap_output uses the output argument
    # output will by mltreemap_output/update_tree
    update_tree.add_argument("--update_tree", action="store_true", default=False,
                             help="Flag indicating the reference tree specified by `--reftree` "
                                  "is to be updated using the sequences found in MLTreeMap output")
    update_tree.add_argument("--uclust", required=False, default=False, action="store_true",
                             help="Cluster sequences that mapped to the reference tree prior to updating")
    update_tree.add_argument("--gap_removal", required=False, default=False, action="store_true",
                             help="Remove minor gaps using Gblocks?")
    update_tree.add_argument("-u", "--uclust_identity", required=False, default=0.97, type=float,
                             help="Sequence identity value to be used in uclust [DEFAULT = 0.97]")
    update_tree.add_argument("-a", "--alignment_mode", required=False, default='d', type=str, choices=['d', 'p'],
                             help="Alignment mode: 'd' for default and 'p' for profile-profile alignment")

    miscellaneous_opts = parser.add_argument_group("Miscellaneous options")
    miscellaneous_opts.add_argument('--overwrite', action='store_true', default=False,
                                    help='overwrites previously processed output folders')
    miscellaneous_opts.add_argument('-v', '--verbose', action='store_true',  default=False,
                                    help='Prints a more verbose runtime log')
    miscellaneous_opts.add_argument("--check_trees", action="store_true", default=False,
                                    help="Quality-check the reference trees before running MLTreeMap")
    miscellaneous_opts.add_argument('-T', '--num_threads', default=2, type=int,
                                    help='specifies the number of CPU threads to use in RAxML and BLAST [DEFAULT = 2]')
    miscellaneous_opts.add_argument('-d', '--delete', default=None,
                                    help='the sections of files to be deleted, as separated by colons '
                                         '(1 = Sequence Files; '
                                         ' 2 = BLAST Results; '
                                         ' 3 = Genewise Results; '
                                         ' 4 = hmmalign and Gblocks Results; '
                                         ' 5 = Unparsed RAxML Results)')

    return parser


def find_executables(args):
    """
    Finds the executables in a user's path to alleviate the requirement of a sub_binaries directory
    :param args: command-line arguments objects
    :return: exec_paths beings the absolute path to each executable
    """
    exec_paths = dict()
    dependencies = ["blastn", "blastx", "blastp", "genewise", "Gblocks", "raxmlHPC", "hmmalign", "hmmbuild"]

    if args.rpkm:
        dependencies += ["bwa", "rpkm"]

    if args.update_tree:
        dependencies.append("usearch")

    if os_type() == "linux":
        args.executables = args.mltreemap + "sub_binaries" + os.sep + "ubuntu"
    if os_type() == "mac":
        args.executables = args.mltreemap + "sub_binaries" + os.sep + "mac"
    elif os_type() == "win" or os_type() is None:
        sys.exit("ERROR: Unsupported OS")

    for dep in dependencies:
        if is_exe(args.executables + os.sep + dep):
            exec_paths[dep] = str(args.executables + os.sep + dep)
        elif which(dep):
            exec_paths[dep] = which(dep)
        else:
            sys.stderr.write("Could not find a valid executable for " + dep + ". ")
            sys.exit("Bailing out.")

    args.executables = exec_paths
    return args


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


def check_parser_arguments(parser):
    """
    Ensures the command-line arguments returned by argparse are sensical
    :param parser: object with parameters returned by argparse
    :return 'args', a summary of MLTreeMap settings.
    """

    # Ensure files contain more than 0 sequences
    args = parser.parse_args()
    args.mltreemap = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + os.sep

    # Set the reference data file prefix and the reference tree name
    if args.reftree == 'g':
        args.reference_data_prefix = 'geba_'
        args.reference_tree = 'geba.tree'
    elif args.reftree == 'i':
        args.reference_data_prefix = 'fungi_'
        args.reference_tree = 'fungitr_tree.txt'
    elif args.reftree == 'p':
        args.reference_data_prefix = ''
        args.reference_tree = 'MLTreeMap_reference.tree'
    else:
        args.reference_data_prefix = ''
        args.reference_tree = args.reftree + "ref_tree.txt"

    # Notify the user that bootstraps cannot be used with the Maximum Parsimony settings of RAxML.
    if args.bootstraps > 1 and args.phylogeny == 'p':
        sys.stderr.write('WARNING: You intended to do ' + str(args.bootstraps) +
                         ' bootstrap replicates. Bootstrapping is disabled in the parsimony mode of MLTreeMap.' +
                         ' The pipeline will continue without bootstrapping.\n')
        sys.stderr.flush()
        args.bootstraps = 1

    args = find_executables(args)

    # Add (or replace a trailing (back)slash with) the os.sep to the end of the output directory
    while re.search(r'/\Z', args.output) or re.search(r'\\\Z', args.output):
        args.output = args.output[:-1]
    args.output += os.sep

    if sys.version_info > (2, 9):
        args.py_version = 3
    else:
        args.py_version = 2

    args.output_dir_var = args.output + 'various_outputs' + os.sep
    args.output_dir_raxml = args.output + 'final_RAxML_outputs' + os.sep
    args.output_dir_final = args.output + 'final_outputs' + os.sep

    mltreemap_dir = args.mltreemap + os.sep + 'data' + os.sep
    genewise_support = mltreemap_dir + os.sep + 'genewise_support_files' + os.sep

    if args.num_threads >= available_cpu_count():
        sys.stderr.write("WARNING: Number of threads specified is greater than those available! "
                         "Using maximum threads available (" + str(available_cpu_count()) + ")\n")
        sys.stderr.flush()
        args.num_threads = available_cpu_count()

    # TODO: make this soultion a bit better
    if os.getenv("WISECONFIGDIR") is None:
        sys.stderr.write("ERROR: WISECONFIGDIR not set!\n")
        sys.exit("export WISECONFIGDIR=" + genewise_support + os.sep + "wisecfg")

    if args.rpkm:
        if not args.reads:
            sys.stderr.write("ERROR: At least one FASTQ file must be provided if -rpkm flag is active!")
            sys.exit()
        if args.reverse and not args.reads:
            sys.stderr.write("ERROR: File containing reverse reads provided but forward mates file missing!")
            sys.exit()

    return args


def get_response(py_version, response_string=""):
    if py_version == 3:
        return input(response_string)
    if py_version == 2:
        return raw_input(response_string)


def check_previous_output(args):
    """
    Prompts the user to determine how to deal with a pre-existing output directory.
    :param args: Command-line argument object from get_options and check_parser_arguments
    :return An updated version of 'args', a summary of MLTreeMap settings.
    """

    # delete previous output folders by force
    if args.overwrite:
        if path.exists(args.output):
            shutil.rmtree(args.output)

    args.skip = 'n'
    if path.exists(args.output):
        sys.stdout.write("MLTreeMap output directory " + args.output + " already exists.\n")
        sys.stdout.flush()
        if args.update_tree:
            args.skip = get_response(args.py_version, "Should this be used for updating? [y|n]")
            while not args.skip == 'y' and not args.skip == 'n':
                args.skip = get_response(args.py_version, "Invalid response. Should this be used for updating? [y|n]")
        elif args.rpkm:
            args.skip = get_response(args.py_version, "Should this be used for RPKM calculation? [y|n]")
            while not args.skip == 'y' and not args.skip == 'n':
                args.skip = get_response(args.py_version, "Invalid response. Should this be used for updating? [y|n]")
        else:
            # Prompt the user to deal with the pre-existing output directory
            while os.path.isdir(args.output):
                sys.stdout.write('Overwrite [1], quit [2], or change directory [3]?\n')
                answer = int(get_response(args.py_version))

                while not answer == 1 and not answer == 2 and not answer == 3:
                    answer = int(get_response(args.py_version, 'Invalid input. Please choose 1, 2, or 3.\n'))
                if answer == 1:
                    sys.stdout.write('Do you really want to overwrite the old output directory?\n')
                    sys.stdout.write('All data in it will be lost!\n')
                    answer2 = get_response(args.py_version, 'Yes [y] or no [n]?\n')
                    while not answer2 == 'y' and not answer2 == 'n':
                        answer2 = get_response(args.py_version, 'Invalid input. Please choose y or n.\n')
                    if answer2 == 'y':
                        shutil.rmtree(args.output)
                    else:
                        sys.exit('Exit MLTreeMap\n')
                elif answer == 2:
                    sys.exit('Exit MLTreeMap\n')
                else:
                    args.output = get_response(args.py_version, 'Please enter the path to the new directory.\n')
    
    # Create the output directories
    if not os.path.isdir(args.output):
        os.makedirs(args.output)
        os.mkdir(args.output_dir_var)
        os.mkdir(args.output_dir_raxml)
        os.mkdir(args.output_dir_final)

    return args


def create_cog_list(args):
    """
    Loads the MLTreeMap COG list file and check that the args.reftree exists
    :param args: The command-line and default arguments object
    :return: An autovivification of the COGs in cog_list.txt. This also includes their short-form name (termed
    denominator e.g. mcrA, mcrB, a, b, cdhA, etc.) and a list of output text precursors based on the analysis type.
    The denominator is equal to the command-line reference tree specifier argument (p, g, or i) if phylogenetic COGs
    """
    
    cog_list = Autovivify()
    text_of_analysis_type = Autovivify()
    cog_list_file = args.mltreemap + os.sep + 'data' + os.sep + 'tree_data' + os.sep + 'cog_list.txt'
    cog_input_list = open(cog_list_file, 'r')
    # TODO: alter function so the GEBA and fungi trees can be used instead of just the MLTreeMap reference
    if args.reftree not in ['i', 'p', 'g']:
        alignment_set = ''
    else:
        alignment_set = args.reftree
    kind_of_cog = ''

    # For each line in the COG list file...

    cog_list_lines = [x.strip() for x in cog_input_list.readlines()]
    # Close the COG list file
    cog_input_list.close()

    for cogInput in cog_list_lines:
        # Get the kind of COG if cogInput is a header line        
        if re.match(r'\A#(.*)', cogInput):
            kind_of_cog = re.match(r'\A#(.*)', cogInput).group(1)
            continue

        if re.match(r'\A%(.*)', cogInput):
            continue

        # Add data to COG list based on the kind of COG it is
        if kind_of_cog == 'phylogenetic_cogs':
            cog_list[kind_of_cog][cogInput] = alignment_set
            cog_list['all_cogs'][cogInput] = alignment_set
            text_inset = ''
            if alignment_set == 'g':
                text_inset = ' based on the GEBA reference'
            if alignment_set == 'i':
                text_inset = ' focusing only on fungi'
            text_of_analysis_type[alignment_set] = 'Phylogenetic analysis' + text_inset + ':'
        elif kind_of_cog == 'phylogenetic_rRNA_cogs':
            cog, denominator, text = cogInput.split('\t')
            cog_list[kind_of_cog][cog] = denominator
            cog_list['all_cogs'][cog] = denominator
            text_of_analysis_type[denominator] = 'Phylogenetic analysis, ' + text + ':'
        elif kind_of_cog == 'functional_cogs':
            cog, denominator, text = cogInput.split('\t')
            cog_list[kind_of_cog][cog] = denominator
            cog_list['all_cogs'][cog] = denominator
            text_of_analysis_type[denominator] = 'Functional analysis, ' + text + ':'

    if args.reftree not in ['i', 'g', 'p']:
        if args.reftree not in cog_list['all_cogs'].values():
            sys.stderr.write("ERROR: " + args.reftree +
                             " not found in " + cog_list_file + "! Please use a valid reference tree ID!")
            sys.stderr.flush()
            sys.exit()
    elif args.reftree in ['i', 'g', 'p'] and args.update_tree:
        sys.stderr.write("ERROR: Unable to update all reference trees at the same time! "
                         "Please specify the reference tree to update with the '--reftree' argument and retry.\n")
        sys.stderr.flush()
        sys.exit()

    return cog_list, text_of_analysis_type


def single_cog_list(reftree, cog_list, text_of_analysis_type):
    """
    Copies the relevant information from the cog_list and text_of_analysis_type to new dictionaries
    :param reftree: The reference gene family
    :param cog_list: The vivification of cog names, reference names, and their respective type
    :param text_of_analysis_type: Mapping of cogs and analysis type based on their type
    :return: Pared down versions of cog_list and text_of_analysis_type containing only information of gene
    """
    new_list = Autovivify()
    single_text_analysis = Autovivify()

    # Parse the cog_list
    for cog_type in cog_list:
        for cog, denominator in cog_list[cog_type].iteritems():
            if denominator == reftree:
                new_list[cog_type][cog] = denominator
                break

    # Parse the text_of_analysis_type
    for denominator, analysis in text_of_analysis_type.iteritems():
        if denominator == reftree:
            single_text_analysis[denominator] = analysis
            break

    return new_list, single_text_analysis


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
    assert isinstance(base_end, (int, long, float, complex))
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
    reg_nuc = re.compile(r'[acgtACGT]')
    reg_nuc_ambiguity = re.compile(r'[xnXN]')
    reg_aa_ambiguity = re.compile(r'[bxzBXZ]')
    reg_amino = re.compile(r'[acdefghiklmnpqrstvwyACDEFGHIKLMNPQRSTVWY*]')
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
                if len(sequence) > args.gblocks:
                    formatted_fasta_dict[header] = sequence
                else:
                    num_headers -= 1

            sequence = ""
            # Replace all non a-z, A-Z, 0-9, or . characters with a _
            # Then replace the initial _ with a >
            line = re.sub(r'[^a-zA-Z0-9.\r\n]', '_', line)
            line = re.sub(r'\A_', '>', line)
    
            # Because RAxML can only work with file names having length <= 125,
            # Ensure that the sequence name length is <= 100
            if line.__len__() > 100:
                line = line[0:100]
    
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
        sys.stderr.write("ERROR: duplicate header names were detected in " + args.input + "!\n")
        sys.stderr.flush()
        response = raw_input("Do you want the headers to be renamed "
                             "(The headers in the original input and formatted fasta will not match) [Y/N]?\n")
        if response == 'Y':
            sys.stderr.write("Formatting input FASTA with duplicates... ")
            sys.stderr.flush()
            formatted_fasta_dict = format_read_fasta(args, True)
            sys.stderr.write("done.\n")
            sys.stderr.flush()
        else:
            sys.stderr.write("Bailing out. ")
            sys.stderr.flush()

    return formatted_fasta_dict


def build_hmm(msa_file, args):
    gene_family = ".".join(msa_file.split("/")[-1].split('.')[0:-1])
    sys.stdout.write("realigning sequences for " + gene_family)
    hmm_output = args.mltreemap + "/data/hmm_data/" + gene_family + ".hmm"
    if os.path.isfile(hmm_output):
        os.remove(hmm_output)
    command = [args.executables["hmmbuild"], "-s", "--verbose", hmm_output, msa_file, ">> /dev/null"]
    hmmbuild_process = subprocess.Popen(' '.join(command), shell=True, preexec_fn=os.setsid)
    hmmbuild_process.wait()
    return


def validate_inputs(args, cog_list):
    """
    This function filters the files in data/alignment_data/ for sequences that are entirely ambiguity characters
    or if there are any sequences in the MSA that are not the consistent length
    :param args: the command-line and default options
    :param cog_list: Dictionary containing cog identifiers sorted into their classes
    :return: list of files that were edited
    """
    sys.stdout.write("Testing validity of reference trees... ")
    sys.stdout.flush()
    ref_trees = glob.glob(args.mltreemap + os.sep + "data/tree_data/*_tree.txt")
    ref_tree_dict = dict()
    f_cogs = [cog.strip("_") for cog in cog_list["functional_cogs"].keys()]
    for tree_file in ref_trees:
        denominator = os.path.basename(tree_file).strip("_tree.txt")
        denominator = denominator.strip("_")
        if denominator in f_cogs:
            ref_tree_dict[denominator] = tree_file
    status = pparse_ref_trees(denominator_ref_tree_dict=ref_tree_dict, args=args)
    if status is None:
        sys.exit()
    else:
        sys.stdout.write("Reference trees appear to be formatted correctly. Continuing with MLTreeMap.\n")
        sys.stdout.flush()
    return


def run_blast(args, split_files, cog_list):
    """
    Runs the BLAST algorithm on each of the split input files.
    :param args: Command-line argument object from get_options and check_parser_arguments
    :param cog_list: Dictionary containing cog identifiers sorted into their classes
    :param split_files: List of all files that need to be individually used for BLAST calls
    """

    sys.stdout.write("Running BLAST... ")
    sys.stdout.flush()

    excluded_cogs = list()

    # For each file containing a maximum of the specified number of sequences...
    alignment_data_dir = args.mltreemap + os.sep + \
        'data' + os.sep + \
        args.reference_data_prefix + 'alignment_data' + os.sep + \
        "*.fa"
    try:
        os.path.isdir(alignment_data_dir)
    except IOError:
        sys.stderr.write("ERROR: " + alignment_data_dir + "does not exist!")
        sys.stderr.flush()
        sys.exit()

    db_nt = '-db "'
    db_aa = '-db "'

    for fasta in glob.glob(alignment_data_dir):
        cog = os.path.basename(fasta).split('.')[0]
        if cog in cog_list["all_cogs"].keys():
            if re.match(r'.*rRNA\.fa\Z', fasta):
                db_nt += fasta + ' '
            else:
                db_aa += fasta + ' '
        else:
            excluded_cogs.append(cog)

    db_nt += '"'
    db_aa += '"'

    if len(excluded_cogs) > 0:
        with open(args.output+"mltreemap_BLAST_log.txt", 'w') as blast_log:
            blast_log.write("WARNING:\nThe following markers were excluded from the analysis since they were " +
                            "found in " + alignment_data_dir + " but not in " +
                            args.mltreemap + "/data/tree_data/cog_list.txt:\n")
            for ec in excluded_cogs:
                blast_log.write(ec + "\n")
                
    for split_fasta in sorted(split_files):

        # Ensure split_fasta is a .fasta file; save file name if so, die otherwise
        
        if not re.match(r'\A.+/(.+)\.fasta\Z', split_fasta):
            sys.exit('ERROR: Something is wrong with the directory of the BLAST input file!\n')
        else:
            blast_input_file_name = re.match(r'\A.+/(.+)\.fasta\Z', split_fasta).group(1)

        # Run the appropriate BLAST command(s) based on the input sequence type
        if args.molecule == 'n':
            command = args.executables["blastx"] + " " + \
                '-query ' + split_fasta + ' ' + db_aa + ' ' + \
                '-evalue 0.01 -max_target_seqs 20000 ' + \
                '-dbsize 1000000 -outfmt 6 '
            if args.num_threads:
                command += '-num_threads ' + str(int(args.num_threads)) + ' '
            command += '>> ' + args.output_dir_var + blast_input_file_name + '.BLAST_results_raw.txt'
            command += " 2>/dev/null"
            os.system(command)
            command = args.executables["blastn"] + " " + \
                '-query ' + split_fasta + ' ' + db_nt + ' ' + \
                '-evalue 0.01 -max_target_seqs 20000 ' + \
                '-dbsize 1000000 -outfmt 6 '
            if args.num_threads:
                command += '-num_threads ' + str(int(args.num_threads)) + ' '
            command += '>> ' + args.output_dir_var + blast_input_file_name + '.rRNA_BLAST_results_raw.txt'
            command += " 2>/dev/null"
            os.system(command)

        elif args.molecule == 'a':
            command = args.executables["blastp"] + " " + \
                      '-query ' + split_fasta + ' ' + db_aa + ' ' + \
                      '-evalue 0.01 -max_target_seqs 20000 ' + \
                      '-dbsize 1000000 -outfmt 6 '
            if args.num_threads:
                command += '-num_threads ' + str(int(args.num_threads)) + ' '
            command += '>> ' + args.output_dir_var + blast_input_file_name + '.BLAST_results_raw.txt'
            command += " 2> /dev/null"
            os.system(command)

    sys.stdout.write("done.\n")

 
def collect_blast_outputs(args):
    """
    Deletes empty BLAST results files.
    :param args: Command-line argument object from get_options and check_parser_arguments
    Returns a list of non-empty BLAST results files.
    """
    raw_blast_result_files = []

    for blast_result in glob.glob(args.output_dir_var + '*BLAST_results_raw.txt'):
        blast_result.rstrip('\r\n')
        if path.getsize(blast_result) <= 0:
            os.remove(blast_result)
        else:
            raw_blast_result_files.append(blast_result)
    
    return raw_blast_result_files


def parse_blast_results(args, raw_blast_results, cog_list):
    """
    Returns an Autovivification of purified (eg. non-redundant) BLAST hits.
    :param args: Command-line argument object from get_options and check_parser_arguments
    :param raw_blast_results: list of files produced from BLAST alignment
    :param cog_list: list of COGs included in analysis pipeline
    """

    if args.verbose:
        sys.stdout.write("Parsing BLAST results... ")
        sys.stdout.flush()

    reg_cog_id = re.compile(r'.*(.{7})\Z')
    counter = 0
    purified_blast_hits = Autovivify()

    for blast_table in sorted(raw_blast_results):
        try:     
            blast_results = open(blast_table, 'r')
        except IOError:
            sys.stdout.write("ERROR: Cannot open BLAST outputfile " + blast_table)
            continue

        contigs = {}
        identifier = 0
        for line in blast_results:
            # Clear variables referencing the contig, COG, qstart, qend, reference start, reference end, and bitscore
            # Interpret the BLAST hit, and assign the details accordingly
            temp_contig, tempDetailedCOG, _, _, _, _, tempQStart, tempQEnd, tempRStart, tempREnd, _, tempBitScore = line.split('\t')
            tempREnd = int(tempREnd)
            tempRStart = int(tempRStart)
            tempQEnd = int(tempQEnd)
            tempQStart = int(tempQStart)
            tempBitScore = float(tempBitScore)

            # Skip to next BLAST hit if bit score is less than user-defined minimum
            if tempBitScore <= args.bitscore:
                continue

            # Determine the direction of the hit relative to the reference
            direction = 'forward'
            if tempRStart > tempREnd:
                temp = tempRStart
                tempRStart = tempREnd
                tempREnd = temp
                direction = 'reverse'
            if tempQStart > tempQEnd:
                temp = tempQStart
                tempQStart = tempQEnd
                tempQEnd = temp
                if direction == 'reverse':
                    sys.stderr.write("ERROR: Confusing BLAST result!\n")
                    sys.stderr.write("Please notify the authors about " +
                                     temp_contig + ' at ' +
                                     tempDetailedCOG +
                                     " q(" + str(tempQEnd) + '..' + str(tempQStart) + ")," +
                                     " r(" + str(tempREnd) + '..' + str(tempRStart) + ")")
                    sys.stderr.flush()
                    sys.exit()
                direction = 'reverse'

            # Trim COG name to last 7 characters of detailed COG name
            # TK - This will be important to note in the user's manual,
            # especially if we enable people to add their own COGs later
            result = reg_cog_id.match(tempDetailedCOG)
            if result:
                tempCOG = result.group(1)
            else:
                sys.exit('ERROR: Could not detect the COG of sequence ' + tempDetailedCOG)

            # Save contig details to the list
            if temp_contig not in contigs:
                contigs[temp_contig] = {}

            if identifier not in contigs[temp_contig]:
                contigs[temp_contig][identifier] = {}

            contigs[temp_contig][identifier]['bitscore'] = tempBitScore
            contigs[temp_contig][identifier]['cog'] = tempCOG
            contigs[temp_contig][identifier]['seq_start'] = tempQStart
            contigs[temp_contig][identifier]['seq_end'] = tempQEnd
            contigs[temp_contig][identifier]['direction'] = direction
            contigs[temp_contig][identifier]['validity'] = True
            identifier += 1

        # Close the file
        blast_results.close()

        # Purify the BLAST hits
        # For each contig sorted by their string-wise comparison...
        for contig in sorted(contigs.keys()):
            identifier = 0

            # create tuple array to sort
            IDs = []
            for raw_identifier in sorted(contigs[contig].keys()):
                base_start = contigs[contig][raw_identifier]['seq_start']
                IDs.append((raw_identifier, base_start))
            _IDs = sorted(IDs, key=lambda x: x[1])
            IDs = [x[0] for x in _IDs]

            base_blast_result_raw_identifier = IDs.pop()
            contigs[contig][base_blast_result_raw_identifier]['validity'] = True
            base_bitscore = contigs[contig][base_blast_result_raw_identifier]['bitscore']
            base_cog = contigs[contig][base_blast_result_raw_identifier]['cog']
            base_start = contigs[contig][base_blast_result_raw_identifier]['seq_start']
            base_end = contigs[contig][base_blast_result_raw_identifier]['seq_end']
            direction = contigs[contig][base_blast_result_raw_identifier]['direction']
            base_length = base_end - base_start

            # Compare the BLAST hit (base) against all others
            # There may be several opinions about how to do this. This way is based on the original MLTreeMap
            # ----A----  --C--
            #        ---B---
            # A kills B, B kills C. (Another approach would be to let C live,
            # but the original MLTreeMap authors don't expect C to be useful)
            for check_blast_result_raw_identifier in IDs:
                check_bitscore = contigs[contig][check_blast_result_raw_identifier]['bitscore']
                check_cog = contigs[contig][check_blast_result_raw_identifier]['cog']
                check_start = contigs[contig][check_blast_result_raw_identifier]['seq_start']
                check_end = contigs[contig][check_blast_result_raw_identifier]['seq_end']
                check_length = check_end - check_start

                # Compare the base and check BLAST hits
                info = Autovivify()
                info['base']['start'] = base_start
                info['base']['end'] = base_end
                info['check']['start'] = check_start
                info['check']['end'] = check_end
                overlap = calculate_overlap(info)
                counter += 1

                # Check for validity for hits with overlap
                if overlap == 0:
                    base_blast_result_raw_identifier = check_blast_result_raw_identifier
                    base_bitscore = check_bitscore
                    base_cog = check_cog
                    base_start = check_start
                    base_end = check_end
                    base_length = check_length
                    contigs[contig][base_blast_result_raw_identifier]['validity'] = True
                else:
                    if overlap > 0.5*base_length and base_bitscore < check_bitscore:
                        contigs[contig][base_blast_result_raw_identifier]['validity'] = False
                        base_blast_result_raw_identifier = check_blast_result_raw_identifier
                        base_bitscore = check_bitscore
                        base_cog = check_cog
                        base_start = check_start
                        base_end = check_end
                        base_length = check_length
                        contigs[contig][base_blast_result_raw_identifier]['validity'] = True
                    elif overlap > 0.5*check_length and check_bitscore < base_bitscore:
                        contigs[contig][check_blast_result_raw_identifier]['validity'] = False
                    elif base_start == check_start and base_end == check_end:
                        # If both are the same, keep only the one with the smaller identifier
                        if check_blast_result_raw_identifier > base_blast_result_raw_identifier:
                            contigs[contig][check_blast_result_raw_identifier]['validity'] = False

            # Set validity to 0 if COG is not in list of MLTreeMap COGs
            if base_cog not in cog_list['all_cogs']:
                contigs[contig][base_blast_result_raw_identifier]['validity'] = False

            # Save purified hits for valid base hits
            for base_blast_result_raw_identifier in IDs:
                base_bitscore = contigs[contig][base_blast_result_raw_identifier]['bitscore']
                base_cog = contigs[contig][base_blast_result_raw_identifier]['cog']
                base_start = contigs[contig][base_blast_result_raw_identifier]['seq_start']
                base_end = contigs[contig][base_blast_result_raw_identifier]['seq_end']
                direction = contigs[contig][base_blast_result_raw_identifier]['direction']
                if contigs[contig][base_blast_result_raw_identifier]['validity']:
                    purified_blast_hits[contig][identifier]['bitscore'] = base_bitscore
                    purified_blast_hits[contig][identifier]['cog'] = base_cog
                    purified_blast_hits[contig][identifier]['start'] = base_start
                    purified_blast_hits[contig][identifier]['end'] = base_end
                    purified_blast_hits[contig][identifier]['direction'] = direction
                    purified_blast_hits[contig][identifier]['is_already_placed'] = False
                    identifier += 1

    # Print the BLAST results for each contig
    for contig in sorted(purified_blast_hits.keys()):
        outfile = args.output_dir_var + contig + '_blast_result_purified.txt'
        out = open(outfile, 'w')
        sorting_hash = {}

        # Identify the first instance of each bitscore
        for identifier in sorted(purified_blast_hits[contig].keys()):
            if not purified_blast_hits[contig][identifier]['bitscore'] in sorting_hash:
                sorting_hash[purified_blast_hits[contig][identifier]['bitscore']] = {}
            sorting_hash[purified_blast_hits[contig][identifier]['bitscore']][identifier] = 1

        # Print the (potentially reduced set of) BLAST results ordered by decreasing bitscore
        for bitscore in sorted(sorting_hash.keys(), reverse=True):
            for identifier in sorted(sorting_hash[bitscore]):
                out.write(contig + '\t' + str(purified_blast_hits[contig][identifier]['start']) + '\t' +
                          str(purified_blast_hits[contig][identifier]['end']) + '\t' +
                          str(purified_blast_hits[contig][identifier]['direction']) + '\t' +
                          purified_blast_hits[contig][identifier]['cog'] + '\t' + str(bitscore) + '\n')

        out.close()
    if args.verbose:
        sys.stdout.write("done.\n")
    return purified_blast_hits


def blastpParser(args, blast_hits_purified):
    """
    For each contig, produces a file similar to the Genewise output file
    (this is in cases where Genewise is unnecessary because it is already an AA sequence.
    :param args: Command-line argument object from get_options and check_parser_arguments
    :return blastpSummaryFiles: Autovivification of the output file for each contig.
    """

    blastpSummaryFiles = Autovivify()

    regHEADER = re.compile(r'\A>')

    for contig in sorted(blast_hits_purified.keys()):
        output_file = args.output_dir_var + contig + '_blast_result_summary.txt'
        try:
            output = open(output_file, 'w')
        except IOError:
            sys.exit('ERROR: Unable to open ' + output_file + '!\n')
        blastpSummaryFiles[contig][output_file] = 1
        shortened_sequence_file = args.output_dir_var + contig + '_sequence_shortened.txt'
        try:
            sequence_file = open(shortened_sequence_file, 'r')
        except IOError:
            sys.exit('ERROR: Could not open ' + shortened_sequence_file + '!\n')
        flagSeq = 0
        sequence = ''

        # Get the sequence from the shortened sequence file
        for line in sequence_file:
            if regHEADER.search(line):
                if flagSeq == 1:
                    sys.exit('ERROR: Unexpected multiple shortened sequences found!\n')
                flagSeq = 1
                continue
            else:
                line.strip()
                sequence += line

        # Write the output file to imitate the Genewise results
        for count in sorted(blast_hits_purified[contig].keys()):
            output.write(str(blast_hits_purified[contig][count]['cog']) + '\t')
            output.write(str(blast_hits_purified[contig][count]['start']) + '\t')
            output.write(str(blast_hits_purified[contig][count]['end']) + '\t')
            output.write(str(blast_hits_purified[contig][count]['direction']) + '\t')
            output.write(str(sequence) + '\n')
        sequence_file.close()
        output.close()

    return blastpSummaryFiles


def make_genewise_inputs(args, blast_hits_purified, formatted_fasta_dict):
    """
    Takes an Autovivification of purified BLAST hits and uses these to produce the input files needed for Genewise.

    Returns an Autovivification mapping the contig to its sequence's start and end positions for Genewise.
    Returns a list of files to be run through Genewise.
    """
    if args.verbose:
        sys.stdout.write("Producing Genewise input files... ")
        sys.stdout.flush()

    flanking_length = 1000  # Recommended: 1000
    prae_contig_coordinates = Autovivify()
    contig_coordinates = Autovivify()
    gene_coordinates = Autovivify()
    shortened_sequence_files = {}

    for contig in sorted(blast_hits_purified.keys()):
        nr_of_blast_hits = len(blast_hits_purified[contig].keys())
        for base_identifier in sorted(blast_hits_purified[contig].keys()):
            # Skip rRNA hits for now (we work with them later)
            # if re.search("rRNA", blast_hits_purified[contig][base_identifier]['cog']):
            #     continue

            # Skip hits which have already been placed; otherwise, mark them as placed
            if blast_hits_purified[contig][base_identifier]['is_already_placed']:
                continue

            blast_hits_purified[contig][base_identifier]['is_already_placed'] = True
            base_start = blast_hits_purified[contig][base_identifier]['start'] - flanking_length
            base_end = blast_hits_purified[contig][base_identifier]['end'] + flanking_length
            check_identifier = 0
            while check_identifier < nr_of_blast_hits:
                # Skip rRNA hits for now (we work with them later)
                # if re.search(r'rRNA', blast_hits_purified[contig][check_identifier]['cog']):
                #     check_identifier += 1
                #     continue

                # Skip hits which have already been placed; otherwise, mark them as placed
                if blast_hits_purified[contig][check_identifier]['is_already_placed']:
                    check_identifier += 1
                    continue

                check_start = blast_hits_purified[contig][check_identifier]['start'] - flanking_length
                check_end = blast_hits_purified[contig][check_identifier]['end'] + flanking_length

                # Check for overlap
                if base_start <= check_start and check_start <= base_end and base_end <= check_end:
                    # Base  --------
                    # Check     --------
                    base_end = check_end
                    blast_hits_purified[contig][check_identifier]['is_already_placed'] = True
                    check_identifier = 0
                    continue
                elif base_start <= check_start and check_end <= base_end:
                    # Base  --------
                    # Check   ----
                    blast_hits_purified[contig][check_identifier]['is_already_placed'] = True
                    check_identifier = 0
                    continue
                elif check_start <= base_start and base_start <= check_end and check_end <= base_end:
                    # Base      --------
                    # Check --------
                    base_start = check_start
                    blast_hits_purified[contig][check_identifier]['is_already_placed'] = True
                    check_identifier = 0
                    continue
                elif check_start <= base_start and base_end <= check_end:
                    # Base    ----
                    # Check --------
                    base_start = check_start
                    base_end = check_end
                    blast_hits_purified[contig][check_identifier]['is_already_placed'] = True
                    check_identifier = 0
                    continue
                check_identifier += 1

            prae_contig_coordinates[contig][base_start][base_end] = blast_hits_purified[contig][base_identifier]['cog']

    # Produce the input files for Genewise
    for contig_name in prae_contig_coordinates:
        sequence = formatted_fasta_dict[">" + contig_name]
        sequence_length = len(sequence)
        shortened_sequence = ""
        # Start searching for the information to shorten the file.
        # Creates a chimera of the sequence if there are multiple hits.
        # shortened_sequence = sequence + sequence[start_blast_hit2:] + ... + sequence[start_blast_hitn:]
        for start_blast in sorted(prae_contig_coordinates[contig_name].keys()):
            for end_blast in sorted(prae_contig_coordinates[contig_name][start_blast].keys()):
                marker_gene = prae_contig_coordinates[contig_name][start_blast][end_blast]

                # Ok, now we have all information about the hit. Correct start and end if needed:
                if start_blast < 0:
                    start_blast = 0

                if end_blast >= sequence_length:
                    end_blast = sequence_length - 1

                gene_coordinates[contig_name][start_blast][end_blast] = marker_gene
                # Skip rRNA hits for now (we work with them later)
                if re.search("rRNA", marker_gene):
                    continue

                # Note: Genewise (gw) positions start with 1, blast positions with 0 ->
                # thus differentiate between start_blast and start_gw
                shortened_start_gw = len(shortened_sequence) + 1
                count = -1
                for nucleotide in sequence:
                    count += 1
                    if not count >= start_blast and count <= end_blast:
                        continue
                    shortened_sequence += nucleotide

                shortened_end_gw = len(shortened_sequence)
                addition_factor = (start_blast + 1) - shortened_start_gw  # $start_B + 1 == $start_GW
                contig_coordinates[contig_name][shortened_start_gw][shortened_end_gw] = addition_factor

        # Skip rRNA hits for now (we work with them later)
        if re.search("rRNA", marker_gene):
            continue

        try:
            with open(args.output_dir_var + contig_name + "_sequence.txt", 'w') as f:
                fprintf(f, "%s\n", ">" + contig_name + "\n" + sequence)
            f.close()
        except:
            raise IOError("Can't create " + args.output_dir_var + contig_name + "_sequence.txt!")

        try:
            with open(args.output_dir_var + contig_name + "_sequence_shortened.txt", 'w') as f:
                fprintf(f, "%s\n", ">" + contig_name + "\n" + shortened_sequence)
            f.close()
            shortened_sequence_files[args.output_dir_var + contig_name + "_sequence_shortened.txt"] = contig_name
        except:
            raise IOError("Can't create " + args.output_dir_var + contig_name + "_sequence_shortened.txt!")

    if args.verbose:
        sys.stdout.write("done.\n")
    return contig_coordinates, shortened_sequence_files, gene_coordinates


def write_nuc_sequences(args, gene_coordinates, formatted_fasta_dict):
    """
    Function to write the nucleotide sequences representing the BLAST alignment region for each hit in the fasta
    :param args:
    :param gene_coordinates:
    :param formatted_fasta_dict:
    :return: nothing
    """
    # Header format:
    # >contig_name|marker_gene|start_end
    input_multi_fasta = re.match(r'\A.*\/(.*)', args.input).group(1)
    orf_nuc_fasta = args.output_dir_var + '.'.join(input_multi_fasta.split('.')[:-1]) + "_genes.fna"
    try:
        fna_output = open(orf_nuc_fasta, 'w')
    except:
        raise IOError("Unable to open " + orf_nuc_fasta + " for writing!")

    for contig_name in gene_coordinates:
        start = 0
        end = 0
        for coords_start in sorted(gene_coordinates[contig_name].keys()):
            start = coords_start
            for coords_end in gene_coordinates[contig_name][coords_start].keys():
                end = coords_end
                cog = gene_coordinates[contig_name][coords_start][coords_end]
                fna_output.write('>' + contig_name + '|' + cog + '|' + str(start) + '_' + str(end) + "\n")
                fna_output.write(formatted_fasta_dict['>' + contig_name][start:end] + "\n")

    return


def fprintf(opened_file, fmt, *args):
    """
    A helper function used to print to a specified file.
    :param opened_file: A file object that has already been opened using open()
    """
    opened_file.write(fmt % args)


def add_tasks_to_queue(task_list, task_queue, num_threads):
    """
    Function for adding genewise commands from task_list to task_queue while ensuring space in the JoinableQueue
    :param task_list: List of genewise commands
    :param task_queue: JoinableQueue object with a maximum size of 32767
    :param num_threads: Number of threads to be used
    :return: Nothing
    """
    num_tasks = len(task_list)
    if num_tasks == 0:
        raise AssertionError("Task list is empty - nothing to do. Exiting.")

    task = task_list.pop()
    while task:
        if not task_queue.full():
            task_queue.put(task)
            if num_tasks > 1:
                task = task_list.pop()
                num_tasks -= 1
            else:
                task = None

    i = int(num_threads)
    while i:
        if not task_queue.full():
            task_queue.put(None)
            i -= 1

    return


def start_genewise(args, shortened_sequence_files, blast_hits_purified):
    """
    Runs Genewise on the provided list of sequence files.
    (The provided Autovivification of purified BLAST hits is used for file naming purposes).

    Returns an Autovivification mapping the Genewise output files to each contig.
    """

    max_size = 32767  # The actual size limit of a JoinableQueue
    task_queue = JoinableQueue(max_size)
    task_list = list()
    genewise_process_queues = [GenewiseWorker(task_queue) for i in range(int(args.num_threads))]
    for process in genewise_process_queues:
        process.start()

    if args.verbose:
        sys.stdout.write("Running Genewise... ")
        sys.stdout.flush()

    mltreemap_dir = args.mltreemap + os.sep + 'data' + os.sep
    genewise_support = mltreemap_dir + os.sep + 'genewise_support_files' + os.sep
    hmm_dir = mltreemap_dir + "hmm_data" + os.sep

    genewise_outputfiles = Autovivify()

    # This is not working on linux machines. Having to revert to command-line calls
    # if os.getenv("WISECONFIGDIR") is None:
    #     sys.stderr.write("genewise exception\n")
    #     sys.stderr.flush()
    #     os.environ["WISECONFIGDIR"] = genewise_support + os.sep + "wisecfg"

    hmm_dir_files = [f for f in os.listdir(hmm_dir) if os.path.isfile(join(hmm_dir, f))]

    cog_hmms = ['.'.join(hmmF.split('.')[:-1]) for hmmF in hmm_dir_files]

    # For each file which has been shortened by make_genewise_inputs...
    for shortened_sequence_file in sorted(shortened_sequence_files.keys()):
        contig = shortened_sequence_files[shortened_sequence_file]
    
        # For each identifier associated with this contig in the output of parse_blast_results
        for identifier in sorted(blast_hits_purified[contig].keys()):
            cog = blast_hits_purified[contig][identifier]['cog']
            if cog not in cog_hmms:
                continue

            # Prepare the output file name, and store it
            genewise_outputfile = args.output_dir_var + contig + '_' + cog + '_genewise.txt'
            genewise_outputfiles[contig][genewise_outputfile] = 1

            # Prepare the Genewise command and run it
            genewise_command = [args.executables["genewise"], 
                                hmm_dir + cog + ".hmm"]
            genewise_command += [shortened_sequence_file, "-init", "local", "-quiet"]
            genewise_command += ["-gene", genewise_support + 'human.gf']
            genewise_command += ["-matrix", genewise_support + "blosum62.bla"]
            genewise_command += ["-codon", genewise_support + "codon.table"]
            genewise_command.append("-hmmer")
            genewise_command += ["-subs", str(0.01)]
            genewise_command += ["-indel", str(0.01)]
            genewise_command += ["-gap", str(11)]
            genewise_command += ["-ext", str(1)]
            genewise_command += ["-both", "-pep", "-sum", ">", genewise_outputfile]

            task_list.append(genewise_command)

    add_tasks_to_queue(task_list, task_queue, args.num_threads)

    task_queue.close()
    task_queue.join()

    # Return the list of output files for each contig
    if args.verbose:
        sys.stdout.write("done.\n")
    return genewise_outputfiles


def parse_genewise_results(args, genewise_outputfiles, contig_coordinates):
    """
    Uses the provided Autovivification of Genewise output files and the provided
    Autovivification mapping the contig to its Genewise sequence's start and end
    points to produce files summarizing the purified Genewise results.

    Returns an Autovivification mapping the summary files to each contig.
    """

    if args.verbose:
        sys.stdout.write("Parsing Genewise outputs... ")
        sys.stdout.flush()

    genewise_summary_files = Autovivify()

    # For each contig analyzed by Genewise...
    for contig in sorted(genewise_outputfiles.keys()):
        genewise_results_raw = Autovivify()
        genewise_results = Autovivify()
        at_least_one_hit = 0
        count = 0

        # Parse each output file of that contig
        for genewise_outputfile in sorted(genewise_outputfiles[contig].keys()):
            try:     
                genewise_file = open(genewise_outputfile, 'r')
            except IOError:
                sys.stdout.write("ERROR: Cannot open Genewise output file " + genewise_outputfile + "\n")
                sys.exit()

            header_count = 0
            sequence_count = -1

            for line in genewise_file:
                line.strip()

                # If the line starts with a digit, parse it
                if re.match(r'\A\d', line):

                    # Split the results based on one or more spaces between the desired data
                    bitscore, query, _, _, _, start, end, _, _ = re.split(' +', line)
                    bitscore = float(bitscore)
                    start = int(start)
                    end = int(end)

                    # If there is at least one query, take note for future use
                    if query is not None:
                        at_least_one_hit = 1

                    # Determine the direction of the predicted amino acid sequence
                    direction = 'forward'
                    if start > end:
                        temp = start
                        start = end
                        end = temp
                        direction = 'reverse'

                    # Genewise is run on a shortened sequence, so the true positions must be calculated
                    for coords_start in sorted(contig_coordinates[contig].keys()):
                        if start >= coords_start:
                            for coords_end in sorted(contig_coordinates[contig][coords_start].keys()):
                                if end <= coords_end:
                                    addition_factor = contig_coordinates[contig][coords_start][coords_end]
                                    start += addition_factor
                                    end += addition_factor
                                    break

                    genewise_results_raw[contig][genewise_outputfile][header_count]['start'] = start
                    genewise_results_raw[contig][genewise_outputfile][header_count]['end'] = end
                    genewise_results_raw[contig][genewise_outputfile][header_count]['cog'] = query
                    genewise_results_raw[contig][genewise_outputfile][header_count]['bitscore'] = bitscore
                    genewise_results_raw[contig][genewise_outputfile][header_count]['direction'] = direction
                    header_count += 1

                # Otherwise, if the line starts with a '>', prepare to intake the sequence
                elif re.match(r'\A>', line):
                    sequence_count += 1
                    genewise_results_raw[contig][genewise_outputfile][sequence_count]['sequence'] = ''

                # If the line begins with any word character, and contains neither 'Bits' (a title line)
                # nor 'Making' (a Genewise comment about the treatment of introns)
                elif re.match(r'\A\w', line) and not re.match(r'\ABits', line) and not re.match(r'\AMaking', line):
                    genewise_results_raw[contig][genewise_outputfile][sequence_count]['sequence'] += line.strip()

            genewise_file.close()

        # Skip to next contig if there isn't at least 1 hit
        if at_least_one_hit != 1:
            continue

        for base_genewise_output_file in sorted(genewise_results_raw[contig].keys()):

            # For each count of the genewise_outputfile...
            for base_count in sorted(genewise_results_raw[contig][base_genewise_output_file].keys()):
                base_start = genewise_results_raw[contig][base_genewise_output_file][base_count]['start']
                base_end = genewise_results_raw[contig][base_genewise_output_file][base_count]['end']
                base_cog = genewise_results_raw[contig][base_genewise_output_file][base_count]['cog']
                base_bitscore = genewise_results_raw[contig][base_genewise_output_file][base_count]['bitscore']
                base_direction = genewise_results_raw[contig][base_genewise_output_file][base_count]['direction']
                base_sequence = genewise_results_raw[contig][base_genewise_output_file][base_count]['sequence']

                # Ensure that the base_cog, base_start, and base_end are defined
                if base_cog is None or base_start is None or base_end is None:
                    error_string = 'ERROR: The file "' + base_genewise_output_file + '" cannot be parsed!\n' +\
                                   'Please contact the authors about it. As a quick solution to the problem, ' +\
                                   'try to remove the sequence which produced this hit from your input file.\n'
                    sys.exit(error_string)
                is_valid = 1

                base_length = base_end - base_start

                # Check against all other genewise_outputfiles for that contig
                # For each genewise_outputfile for the contig...
                for check_genewise_outputfile in sorted(genewise_results_raw[contig].keys()):

                    # For each count of the genewise_outputfile...
                    for check_count in sorted(genewise_results_raw[contig][check_genewise_outputfile].keys()):

                        # Skip to next iteration if comparing the base to itself
                        if base_count == check_count:
                            continue
                        check_start = genewise_results_raw[contig][check_genewise_outputfile][check_count]['start']
                        check_end = genewise_results_raw[contig][check_genewise_outputfile][check_count]['end']
                        check_cog = genewise_results_raw[contig][check_genewise_outputfile][check_count]['cog']

                        # Ensure that the check_cog, check_start, and check_end are defined
                        if check_cog is None or check_start is None or check_end is None:
                            error_string = 'ERROR: The file "' + check_genewise_outputfile + '" cannot be parsed!\n' +\
                                           'Please contact the authors about it. As a quick solution to the problem, ' +\
                                           'try to remove the sequence which produced this hit from your input file.\n'
                            sys.exit(error_string)

                        check_length = check_end - check_start
                        info = Autovivify()
                        info['base']['start'] = base_start
                        info['base']['end'] = base_end
                        info['check']['start'] = check_start
                        info['check']['end'] = check_end
                        overlap = calculate_overlap(info)

                        # Purify the results
                        # If the size of the overlap is more than half the size of the hit...
                        if float(overlap) / float(base_length) > 0.5:

                            # And if the hit and check are the same COG...
                            if base_cog == check_cog:

                                # Keep only the longer hit, since the major difference between the hits is the length 
                                if base_length < check_length:
                                    is_valid = 0

                            # The COGs are different,
                            # so only skip the hit if it is less than half the length of the check
                            elif base_length < check_length / 2:
                                is_valid = 0

                        # But if the overlap is not more than half the size of the hit, and the hit remains valid...
                        if is_valid and base_cog == check_cog:

                            # Invalidate the hit if it is a side hit of the same COG
                            if base_length < check_length * 0.7:
                                is_valid = 0

                # If the hit is valid, save it
                if is_valid == 1:
                    genewise_results[contig][count]['start'] = base_start
                    genewise_results[contig][count]['end'] = base_end
                    genewise_results[contig][count]['cog'] = base_cog
                    genewise_results[contig][count]['direction'] = base_direction
                    genewise_results[contig][count]['sequence'] = base_sequence
                    count += 1

        # Skip to next hit if there are no valid hits
        if count <= 0:
            sys.stdout.write("Number of valid hits for " + contig + " = "  + str(count))
            continue

        # Write the summary file
        genewise_summary_file = args.output_dir_var + contig + '_genewise_result_summary.txt'
        try: 
            output = open(genewise_summary_file, 'w')
        except IOError:
            sys.stdout.write('ERROR: Cannot open Genewise summary file ' + genewise_summary_file + ' for writing')
            sys.exit(0)

        genewise_summary_files[contig][genewise_summary_file] = 1
        for count in sorted(genewise_results[contig].keys()):
            output.write(genewise_results[contig][count]['cog'] + '\t' +
                         str(genewise_results[contig][count]['start']) + '\t' +
                         str(genewise_results[contig][count]['end']) + '\t' +
                         genewise_results[contig][count]['direction'] + '\t' +
                         genewise_results[contig][count]['sequence'] + '\n')

        output.close()

    if args.verbose:
        sys.stdout.write("done.\n")
        sys.stdout.flush()
    return genewise_summary_files


def get_rRNA_hit_sequences(args, blast_hits_purified, cog_list, genewise_summary_files):
    """
    rRNA does not get translated into protein. Regardless, we want to take the
    rRNA and summarize it in a way that is parallel to the Genewise summary files.
    This function does that using the provided Autovivification of purified BLAST
    hits, list of COGs, and Autovivification of Genewise summary files.

    Returns an Autovivification summarizing the coordinates for each rRNA hit.
    Returns a list of the rRNA summary files.
    """

# TK: ...the list of rRNA hit files is empty.
    sys.stdout.write("Retrieving rRNA hits... ")
    sys.stdout.flush()
    contig_rRNA_coordinates = Autovivify()
    rRNA_hit_files = {}
    
    for contig in sorted(blast_hits_purified.keys()):
        # note: We skipped the Genewise step (we are dealing with rRNA) but we bring the rRNA files in the
        # same structure as the Genewise summary files and bring them back into the ordinary pipeline.
        for identifier in sorted(blast_hits_purified[contig].keys()):
            if not re.search("rRNA", blast_hits_purified[contig][identifier]['cog']):
                continue

            start = blast_hits_purified[contig][identifier]["start"]
            end = blast_hits_purified[contig][identifier]["end"]
            cog = blast_hits_purified[contig][identifier]["cog"]
            direction = blast_hits_purified[contig][identifier]["direction"]
            contig_rRNA_coordinates[contig][identifier]["start"] = start
            contig_rRNA_coordinates[contig][identifier]["end"] = end
            contig_rRNA_coordinates[contig][identifier]["cog"] = cog
            contig_rRNA_coordinates[contig][identifier]["direction"] = direction
            outfile_name = args.output_dir_var + contig + '_rRNA_result_summary.txt'
            contig_rRNA_coordinates[contig][identifier]["outfile"] = outfile_name
            genewise_summary_files[contig][outfile_name] = 1     
            try:
                outfile = open(outfile_name, 'w')
                outfile.close()
            except IOError:
                sys.stderr.write("ERROR: Can't create " + outfile_name + '!\n')
                sys.exit(0)

    try:
        input = open(args.input, 'r')
    except IOError:
        sys.stdout.write("ERROR: Can't create " + args.input + '!\n')
        sys.exit(0)
    contig_name = ''
    sequence = ''

    line = 'x'
    while line:
        line = input.readline()
        line = line.strip()
        line = re.sub(r'\s', '_', line)
        searchmatch = re.search(r'\A>(.+)', line)

        if searchmatch or not line:
            if not line:
                sequence += line

            if contig_name in contig_rRNA_coordinates:
                sequence_length = len(sequence)
                # start searching for the information to shorten the file.
                for identifier in sorted(contig_rRNA_coordinates[contig_name].keys()):
                    start = contig_rRNA_coordinates[contig_name][identifier]["start"]
                    end = contig_rRNA_coordinates[contig_name][identifier]["end"]
                    cog = contig_rRNA_coordinates[contig_name][identifier]["cog"]
                    direction = contig_rRNA_coordinates[contig_name][identifier]["direction"]
                    outfile = contig_rRNA_coordinates[contig_name][identifier]['outfile']
                    denominator = cog_list['all_cogs'][cog]
                    count = -1
                    shortened_sequence = ""
                    for nucleotide in sequence: 
                        count += 1
                        if not (count >= start and count <= end):
                            continue
                        shortened_sequence += nucleotide

                    if direction == 'reverse':
                        # ok, our hit has been on the opposite strand of the reference.
                        # to get a proper alignment, we thus have to produce a negative strand version of the input
                        nucleotides2 = ''.join(reversed(shortened_sequence))
                        shortened_sequence = ""
                        nucleotides2 = nucleotides2.lower()
                        for nucleotide in nucleotides2:
                            if nucleotide == 't':
                                nucleotide = 'a'
                            elif nucleotide == 'a':
                                nucleotide = 't'
                            elif nucleotide == 'c':
                                nucleotide = 'g'
                            elif nucleotide == 'g':
                                nucleotide = 'c'

                            shortened_sequence += nucleotide
                    try:
                        out = open(outfile, 'a')
                        fprintf(out, '%s\t%s\t%s\t%s\t%s\n', cog, start, end, 'n/a', shortened_sequence)
                        out.close()
                    except IOError:
                        sys.stderr.write("ERROR: Can't create " + outfile + '!\n')
                        sys.exit(0)

                try:
                    output_file = open(args.output_dir_var + contig_name + '_sequence.txt', 'w')
                    fprintf(output_file, '>%s\n%s', contig_name, sequence)
                    output_file.close()
                except IOError:
                    sys.stderr.write("ERROR: Can't create " + args.output_dir_var + contig_name + '_sequence.txt!\n')
                    sys.exit(0)

            if searchmatch:
                contig_name = searchmatch.group(1)
                sequence = ""
        else:
            sequence += line
    input.close()
    sys.stdout.write("done.\n")
    return contig_rRNA_coordinates, rRNA_hit_files


def prepare_and_run_hmmalign(args, genewise_summary_files, cog_list):
    """
    Runs hmmalign using the provided COG list and summary of Genewise files.

    Returns an Autovivification of the resulting files from hmmalign.
    """

    reference_data_prefix = args.reference_data_prefix
    hmmalign_singlehit_files = Autovivify()
    if args.verbose:
        sys.stdout.write("Running hmmalign... ")
        sys.stdout.flush()

    # Run hmmalign on each Genewise summary file
    for contig in sorted(genewise_summary_files.keys()):

        for genewise_summary_file in sorted(genewise_summary_files[contig].keys()):
            try:
                input = open(genewise_summary_file, 'r')
            except IOError:
                sys.stderr.write("ERROR: Can't open " + genewise_summary_file + "!\n")
                sys.exit(0)

            line = input.readline()
            line = line.strip()

            while line:
                cog, start, end, _, sequence = line.split('\t')
                denominator = cog_list["all_cogs"][cog]
                f_contig = denominator + "_" + contig
                genewise_singlehit_file = args.output_dir_var + os.sep + \
                                          f_contig + '_' + cog + "_" + str(start) + "_" + str(end)
                hmmalign_singlehit_files[f_contig][genewise_singlehit_file + ".mfa"] = True 
                genewise_singlehit_file_fa = genewise_singlehit_file + ".fa" 
                try:
                    outfile = open(genewise_singlehit_file_fa, 'w')
                    fprintf(outfile, '>query\n%s', sequence)
                    outfile.close()
                except IOError:
                    sys.stderr.write('Can\'t create ' + genewise_singlehit_file_fa + '\n')
                    sys.exit(0)
                mltreemap_resources = args.mltreemap + os.sep + 'data' + os.sep
                hmmalign_command = [args.executables["hmmalign"], '--mapali',
                                    mltreemap_resources + reference_data_prefix + 'alignment_data' +
                                    os.sep + cog + '.fa',
                                    '--outformat', 'Clustal',
                                    mltreemap_resources + reference_data_prefix + 'hmm_data' + os.sep + cog + '.hmm',
                                    genewise_singlehit_file_fa, '>', genewise_singlehit_file + '.mfa']
                os.system(' '.join(hmmalign_command))

                line = input.readline()
                line = line.strip()

            input.close()
    if args.verbose:
        sys.stdout.write("done.\n")
        sys.stdout.flush()
    return hmmalign_singlehit_files
                   

def get_non_wag_cogs(args):
    """
    Returns an Autovivification listing the COGs which don't follow the WAG evolutionary model.
    :param args: Command-line argument object returned by get_options and check_parser_arguments
    """
    denominator = ""
    non_wag_cog_list = Autovivify()
    non_wag_cogs_file = args.mltreemap + os.sep + 'data' + os.sep + 'tree_data' + os.sep + 'non_wag_cogs.txt'
    try:
        cogin = open(non_wag_cogs_file, 'r')
    except IOError:
        sys.exit('ERROR: Can\'t open ' + non_wag_cogs_file + '!\n')

    for line in cogin:
        line = line.strip()
        if re.search(r'\A#(.+)', line):
            denominator = re.search(r'\A#(.+)', line).group(1)
        else:
            cog, model = line.split('\t')
            non_wag_cog_list[denominator][cog] = model

    cogin.close()
    return non_wag_cog_list


def concatenate_hmmalign_singlehits_files(args, hmmalign_singlehit_files, non_wag_cog_list):
    """
    Concatenates the hmmalign files using the provided Autovivifications of hmmalign files and non-WAG COGs.
    :param args: Command-line argument object from get_options and check_parser_arguments
    Returns a list of the files containing the concatenated hmmalign results.
    Returns a list of the model used for each file.
    Returns a list of the number of sequences found in each file.
    """

    # For each type of gene...
    concatenated_mfa_files = {}
    models_to_be_used = {}
    nrs_of_sequences = {}

    if args.verbose:
        sys.stdout.write("Concatenating hmmalign files... ")
        sys.stdout.flush()

    for f_contig in sorted(hmmalign_singlehit_files.keys()):
        # Determine what type of gene is currently represented, or raise an error
        sequences = dict()
        model_to_be_used = ""
        query_sequence = ""
        parsing_order = dict()
        cog_rep_sequences = dict()
        acc = 0

        if re.search(r'\A(.)', f_contig):
            # An issue if there were denominators with underscores
            denominator = f_contig.split('_')[0]
        else:
            sys.exit('ERROR: The analysis type could not be parsed from ' + f_contig + '!\n')

        for hmmalign_singlehit_file in sorted(hmmalign_singlehit_files[f_contig].keys()):
            cog_len = 0
            try:
                hmmalign_msa = open(hmmalign_singlehit_file, 'r')
            except IOError:
                sys.exit('Can\'t open ' + hmmalign_singlehit_file + '!\n')
            reached_data_part = False
            # Determine the best AA model
            if re.search(r'\A.+_(.{7})_\d+_\d+\.mfa\Z', hmmalign_singlehit_file):
                cog = re.search(r'\A.+_(.{7})_\d+_\d+\.mfa\Z', hmmalign_singlehit_file).group(1)
                if cog not in cog_rep_sequences.keys():
                    acc += 1
                cog_rep_sequences[cog] = set()
            else:
                sys.exit('ERROR: The COG could not be parsed from ' + hmmalign_singlehit_file + '!\n')

            if non_wag_cog_list[denominator][cog] and model_to_be_used != 'PROTGAMMAWAG':
                model_to_be_used = non_wag_cog_list[denominator][cog]
            else:
                model_to_be_used = 'PROTGAMMAWAG'
            # Get sequence from file
            for _line in hmmalign_msa:
                line = _line.strip()
                if re.search(r'query', line):
                    reached_data_part = True
                if not reached_data_part:
                    continue
                searchResult = re.search(r'\A(.+) (\S+)\Z', line)
                if searchResult:
                    name_long = searchResult.group(1)
                    sequence_part = searchResult.group(2)
                    sequence_name = ''
                    if re.search(r'query', name_long):
                        query_sequence += sequence_part
                        cog_len += len(sequence_part)

                    elif re.search(r'(\d+)_', name_long):
                        sequence_name = re.search(r'(\d+)_', name_long).group(1)
                        cog_rep_sequences[cog].add(sequence_name)
                        if sequence_name not in sequences.keys():
                            sequences[sequence_name] = dict()
                        if cog not in sequences[sequence_name].keys():
                            sequences[sequence_name][cog] = ""
                        sequences[sequence_name][cog] += sequence_part

            parsing_order[acc] = cog, cog_len
            hmmalign_msa.close()

        models_to_be_used[f_contig] = model_to_be_used
        concatenated_mfa_files[f_contig] = args.output_dir_var + f_contig + '.mfa'
        # Write to the output file
        try:
            output = open(args.output_dir_var + f_contig + '.mfa', 'w')
        except IOError:
            sys.exit('ERROR: Can\'t create ' + args.output_dir_var + f_contig + '.mfa\n')
        output.write('>query\n' + query_sequence + '\n')
        nrs_of_sequences[f_contig] = 1
        qlen = len(query_sequence)

        for sequence_name in sequences.keys():
            nrs_of_sequences[f_contig] += 1
            sequence = ""
            for p_order in sorted(parsing_order.keys(), key=int):
                cog, cog_len = parsing_order[p_order]
                # print f_contig, sequence_name, p_order, cog
                if sequence_name not in cog_rep_sequences[cog]:
                    sequence += "." * cog_len
                else:
                    sequence += sequences[sequence_name][cog]
            output.write('>' + sequence_name + '\n' + sequence + '\n')
            if len(sequence) != qlen:
                output.close()
                sys.stderr.write("ERROR: inconsistent sequence lengths between query and concatenated HMM alignments!\n")
                sys.exit("Check " + args.output_dir_var + f_contig + ".mfa for offending sequence " + sequence_name)

        output.close()

    if args.verbose:
        sys.stdout.write("done.\n")

    return concatenated_mfa_files, nrs_of_sequences, models_to_be_used


def start_gblocks(args, concatenated_mfa_files, nrs_of_sequences):
    """
    Runs Gblocks using the provided lists of the concatenated hmmalign files, and the number of sequences in each file.

    Returns a list of files resulting from Gblocks.
    """

    gblocks_files = {}
    if args.verbose:
        sys.stdout.write("Running Gblocks... ")
        sys.stdout.flush()
    
    for f_contig in sorted(concatenated_mfa_files.keys()):
        concatenated_mfa_file = concatenated_mfa_files[f_contig]
        nr_of_sequences = nrs_of_sequences[f_contig]
        min_flank_pos = int(nr_of_sequences * 0.55)
        gblocks_file = concatenated_mfa_file + "-gb"
        log = args.output + os.sep + "mltreemap.gblocks_log.txt"
        gblocks_files[f_contig] = gblocks_file
        gblocks_command = [args.executables["Gblocks"], concatenated_mfa_file]
        gblocks_command += ['-t=p', '-s=y', '-u=n', '-p=t', '-b3=15',
                            '-b4=3', '-b5=h', '-b2=' + str(min_flank_pos),
                            '-e=-gb', '>', log]
        os.system(' '.join(gblocks_command))
        if not os.path.isfile(gblocks_file):
            sys.exit("ERROR: " + gblocks_file + " was not successfully created! Check " + log)
    if args.verbose:
        sys.stdout.write("done.\n")
        sys.stdout.flush()
    return gblocks_files


def produce_phy_file(args, gblocks_files, nrs_of_sequences):
    """
    Produces phy files from the provided list of Gblocks result files, and the number of sequences in each file.

    Returns an Autovivification containing the names of the produced phy files.
    """

    phy_files = Autovivify()
    sequence_lengths = Autovivify()

    # Open each Gblocks result file
    for f_contig in sorted(gblocks_files.keys()):
        sequences_for_phy = Autovivify()
        do_not_continue = 0
        sequences_raw = Autovivify()
        gblocks_file = gblocks_files[f_contig]

        try:
            input = open(gblocks_file, 'r')
        except IOError:
            sys.exit('ERROR: Can\'t open ' + gblocks_file + '!\n')

        for line in input:
            line = line.strip()
            seq_name_search = re.search(r'\A>(.+)', line)
            if seq_name_search:
                seq_name = seq_name_search.group(1)
                # Flag the user if the reference alignment contains the number -666, which is needed later in the code
                if seq_name == -666:
                    sys.exit('ERROR: Your reference alignment contains element with the number -666. ' +\
                             'Please change it, because this number is needed for internal purposes.\n')
                if seq_name == 'query':
                    seq_name = -666
            else:
                line = re.sub(r' ', '', line)
                if seq_name == "":
                    sys.stderr.write("ERROR: The Gblocks output " + gblocks_file + "is not in the required format.")
                    sys.stderr.write("Please make sure that your versions of hmmalign and gblocks are compatible with MLTreeMap.")
                    sys.exit()
                if seq_name in sequences_raw:
                    sequences_raw[seq_name] += line
                else:
                    sequences_raw[seq_name] = line

        input.close()

        # Ensure the sequences contain only valid characters for RAxML
        # for seq_name in sorted(sequences_raw.keys()):
        for seq_name in sequences_raw.keys():
            if do_not_continue == 1:
                continue
            sequence = sequences_raw[seq_name]
            count = 0
            sequence_lengths[f_contig] = len(sequence)
            sequence = re.sub(r'\.', 'X', sequence)
            sequence = re.sub(r'\*', 'X', sequence)
            sequence = re.sub('-', 'X', sequence)

            if re.search(r'\AX+\Z', sequence):
                sequence = re.sub('X', 'V', sequence, 1)
            if seq_name == -666:
                seq_dummy = re.sub('X', '', sequence)
                if len(seq_dummy) < args.gblocks:
                    do_not_continue = 1
                    exit_file_name = args.output_dir_var + f_contig + '_exit_after_Gblocks.txt'
                    try:
                        output = open(exit_file_name, 'w')
                    except IOError:
                        sys.exit('ERROR: Can\'t open ' + exit_file_name + '!\n')
                    output.write('final alignment after gblocks is too short (<' + str(args.gblocks) + 'AAs) ' +
                                 '-  insufficient number of marker gene residues in query sequence.\n')
                    output.close()
                    continue
            #
            # if sequence.count('X') > (0.99*len(sequence)) and seq_name != -666:
            #     print "WARNING: More than 99% of", seq_name, "is unknown sequence!"
            #     print "Removing it from further processing to prevent errors with RAxML."
            #     do_not_continue = 1
            #     exit_file_name = args.output_dir_var + f_contig + '_exit_after_Gblocks.txt'
            #     try:
            #         output = open(exit_file_name, 'w')
            #     except IOError:
            #         sys.exit('ERROR: Can\'t open ' + exit_file_name + '!\n')
            #     output.write(seq_name + 'contained an insufficient number of marker gene residues in alignment ' +
            #                  '- this would cause an error in Gblocks and RAxML.\n')
            #     output.close()
            #     continue

            sub_sequences = re.findall(r'.{1,50}', sequence)

            for sub_sequence in sub_sequences:
                sub_sequence = re.sub('U', 'T', sub_sequence)  # Got error from RAxML when encountering Uracil
                sequences_for_phy[f_contig][count][int(seq_name)] = sub_sequence
                count += 1

        if do_not_continue == 1:
            continue

        # Write the sequences to the phy file
        phy_file_name = args.output_dir_var + f_contig + '.phy'
        phy_files[f_contig] = phy_file_name
        try:
            output = open(phy_file_name, 'w')
        except IOError:
            sys.exit('ERROR: Can\'t open ' + phy_file_name + '!\n')
        nr_of_sequences = nrs_of_sequences[f_contig]
        output.write(' ' + str(nr_of_sequences) + '  ' + str(sequence_lengths[f_contig]) + '\n')

        for count in sorted(sequences_for_phy[f_contig].keys()):
            for seq_name in sorted(sequences_for_phy[f_contig][count].keys()):
                sequence_part = sequences_for_phy[f_contig][count][seq_name]
                if count == 0:
                    print_seqname = seq_name
                    if seq_name == -666:
                        print_seqname = 'query'
                    output.write(str(print_seqname))
                    length = len(str(print_seqname))
                    c = length
                    while c < 10:
                        output.write(' ')
                        c += 1
                output.write(sequence_part + '\n')

            output.write('\n')
        output.close()

    return phy_files


def start_RAxML(args, phy_files, cog_list, models_to_be_used):
    """
    Run RAxML using the provided Autovivifications of phy files and COGs, as well as the list of models used for each COG.

    Returns an Autovivification listing the output files of RAxML.
    Returns an Autovivification containing the reference tree file associated with each functional or rRNA COG.
    """
    sys.stdout.write("Running RAxML... coffee?\n")
    sys.stdout.flush()

    raxml_outfiles = Autovivify()

    bootstrap_replicates = args.bootstraps
    denominator_reference_tree_dict = dict()
    mltree_resources = args.mltreemap + os.sep + 'data' + os.sep
    if os.path.isabs(args.output_dir_var):
        output_dir = args.output_dir_var
    else:
        output_dir = os.getcwd() + os.sep + args.output_dir_var
    for f_contig in sorted(phy_files.keys()):
        # Establish the reference tree file to be used for this contig
        reference_tree_file = mltree_resources + 'tree_data' + os.sep + args.reference_tree
        phy_file = phy_files[f_contig]
        if re.search(r'\A(.)', f_contig):
            denominator = f_contig.split('_')[0]
        if not denominator == 'p' and not denominator == 'g' and not denominator == 'i':
            for cog in sorted(cog_list['all_cogs'].keys()):
                if not cog_list['all_cogs'][cog] == denominator:
                    continue
                reference_tree_file = mltree_resources + 'tree_data' + os.sep + cog + '_tree.txt'
                break

        # Determine the output file names, and remove any pre-existing output files
        if type(denominator) is not str:
            sys.exit("ERROR: " + str(denominator) + " is not string but " + str(type(denominator)))
        if type(reference_tree_file) is not str:
            sys.exit("ERROR: " + str(reference_tree_file) + " is not string but " + str(type(reference_tree_file)))

        if len(reference_tree_file) == 0:
            sys.exit("ERROR: could not find reference tree for " + denominator)
        if denominator not in denominator_reference_tree_dict.keys():
            denominator_reference_tree_dict[denominator] = reference_tree_file
        raxml_files = [output_dir + 'RAxML_info.' + f_contig,
                       output_dir + 'RAxML_labelledTree.' + f_contig,
                       output_dir + 'RAxML_classification.' + f_contig]

        for raxml_file in raxml_files:
            try:
                shutil.rmtree(raxml_file) 
            except OSError:
                pass

        raxml_option = args.phylogeny
        model_to_be_used = models_to_be_used[f_contig]
        if model_to_be_used is None:
            sys.exit('ERROR: No best AA model could be detected for the ML step!\n')
        # Set up the command to run RAxML
        raxml_command = [args.executables["raxmlHPC"], '-m', model_to_be_used]
        if bootstrap_replicates > 1:
            raxml_command += ["-p 12345 -b 12345 -#", str(bootstrap_replicates)]
        # Run RAxML using multiple threads, if CPUs available
        if args.num_threads:
            raxml_command += ['-T', str(int(args.num_threads))]
        else:
            raxml_command += ['-T', '2']
        raxml_command += ['-s', phy_file,
                          '-t', reference_tree_file,
                          '-f', str(raxml_option),
                          '-n', str(f_contig),
                          '-w', str(output_dir),
                          '>', str(output_dir) + str(f_contig) + '_RAxML.txt']
        raxml_pro = subprocess.Popen(' '.join(raxml_command), shell=True, preexec_fn=os.setsid)
        raxml_pro.wait()

    # Rename the RAxML output files

    for f_contig in sorted(phy_files.keys()):
        denominator = ''
        if re.match(r'\A(.)', f_contig):
            denominator = f_contig.split('_')[0]
        move_command = ['mv', str(output_dir) + 'RAxML_info.' + str(f_contig),
                        str(output_dir) + str(f_contig) + '.RAxML_info.txt']
        if os.path.exists(str(output_dir) + 'RAxML_info.' + str(f_contig)):
            os.system(' '.join(move_command))
        if raxml_option == 'v':
            raxml_outfiles[denominator][f_contig]['classification'] = str(output_dir) + \
                                                                      str(f_contig) + \
                                                                      '.RAxML_classification.txt'
            raxml_outfiles[denominator][f_contig]['labelled_tree'] = str(output_dir) + \
                                                                     str(f_contig) + \
                                                                     '.originalRAxML_labelledTree.txt'
            move_command1 = ['mv', str(output_dir) + 'RAxML_classification.' + str(f_contig),
                             str(raxml_outfiles[denominator][f_contig]['classification'])]
            move_command2 = ['mv', str(output_dir) + 'RAxML_originalLabelledTree.' + str(f_contig),
                             str(raxml_outfiles[denominator][f_contig]['labelled_tree'])]
            remove_command = ['rm', str(output_dir) + 'RAxML_labelledTree.' + str(f_contig)]
            if os.path.exists(str(output_dir) + 'RAxML_classification.' + str(f_contig)):
                os.system(' '.join(move_command1))
            if os.path.exists(str(output_dir) + 'RAxML_originalLabelledTree.' + str(f_contig)):
                os.system(' '.join(move_command2))
            if os.path.exists(str(output_dir) + 'RAxML_labelledTree.' + str(f_contig)):
                os.system(' '.join(remove_command))
            else:
                sys.stderr.write("Some files were not successfully created for " + str(f_contig) + "\n")
                sys.stderr.write("Check " + str(output_dir) + str(f_contig) + "_RAxML.txt for an error!\n")
                sys.exit("Bailing out!")
        elif raxml_option == 'p':
            raxml_outfiles[denominator][f_contig] = str(output_dir) + str(f_contig) + '.RAxML_parsimonyTree.txt'
            move_command1 = ['mv', str(output_dir) + 'RAxML_parsimonyTree.' + str(f_contig),
                             str(raxml_outfiles[denominator][f_contig])]
            os.system(' '.join(move_command1))
        else:
            sys.exit('ERROR: The chosen RAxML mode is invalid. This should have been noticed earlier by MLTreeMap.' +
                     'Please notify the authors\n')

    return raxml_outfiles, denominator_reference_tree_dict, len(phy_files.keys())


def pparse_ref_trees(denominator_ref_tree_dict, args):
    pool = Pool(processes=int(args.num_threads))
    ref_trees_dict = dict()

    def log_tree(result):
        marker, terminal_children_of_reference = result
        if terminal_children_of_reference is None:
            sys.stdout.write("Letting threads finish before exiting... ")
            sys.stdout.flush()
        ref_trees_dict[marker] = terminal_children_of_reference

    for denominator in denominator_ref_tree_dict:
        reference_tree_file = denominator_ref_tree_dict[denominator]
        pool.apply_async(func=read_and_understand_the_reference_tree,
                         args=(reference_tree_file, denominator, ),
                         callback=log_tree)
    pool.close()
    pool.join()
    for marker in ref_trees_dict:
        if ref_trees_dict[marker] is None:
            sys.stdout.write("done.\n")
            return None
        else:
            pass
    return ref_trees_dict


def pparse_raxml_out_trees(labelled_trees, args):
    """
    The wrapper command for parsing all trees of a gene family (denominator) in parallel
    :param labelled_trees: Dictionary containing labelled tree files for each f_contig
    :param args: args object (for num_threads)
    :return: Dictionary containing all parsed trees for each contig
    """
    pool = Pool(processes=int(args.num_threads))
    raxml_tree_dict = dict()

    def log_tree(result):
        f_contig, rooted_labelled_trees, insertion_point_node_hash = result
        if rooted_labelled_trees is None:
            pool.terminate()
            sys.exit()
        raxml_tree_dict[f_contig] = [rooted_labelled_trees, insertion_point_node_hash]

    for f_contig in labelled_trees:
        tree_file = labelled_trees[f_contig]
        pool.apply_async(func=read_understand_and_reroot_the_labelled_tree,
                         args=(tree_file, f_contig, ),
                         callback=log_tree)
    pool.close()
    pool.join()
    return raxml_tree_dict


def parse_RAxML_output(args, denominator_reference_tree_dict, tree_numbers_translation,
                       raxml_outfiles, text_of_analysis_type, num_raxml_outputs):
    """
    Parse the RAxML output files.
    :param args: Command-line argument object from get_options and check_parser_arguments
    :param denominator_reference_tree_dict:
    :param tree_numbers_translation:
    :param raxml_outfiles:
    :param text_of_analysis_type:
    :return: An Autovivification of the final RAxML output files.
    """

    raxml_option = args.phylogeny

    sys.stdout.write('Parsing the RAxML outputs...\n')
    sys.stdout.flush()

    final_raxml_output_files = Autovivify()

    if num_raxml_outputs > 50:
        progress_bar_width = 50
        step_proportion = float(num_raxml_outputs) / progress_bar_width
    else:
        progress_bar_width = num_raxml_outputs
        step_proportion = 1

    sys.stdout.write("[%s ]" % (" " * progress_bar_width))
    sys.stdout.write("%")
    sys.stdout.write("\b" * (progress_bar_width + 3))
    sys.stdout.flush()

    acc = 0.0

    try:
        parse_log = open(args.output + os.sep + "mltreemap_parse_RAxML_log.txt", 'w')
    except IOError:
        sys.stderr.write("WARNING: Unable to open " + args.output + os.sep + "mltreemap_parse_RAxML_log.txt!")
        sys.stderr.flush()
        parse_log = sys.stdout

    parse_log.write("Parsing each gene reference tree file found in the input sequences in parallel...")
    parse_log.flush()
    terminal_children_strings_of_ref_denominators = pparse_ref_trees(denominator_reference_tree_dict, args)
    parse_log.write(" done.\n")
    if terminal_children_strings_of_ref_denominators is None:
        sys.exit()
    parse_log.write(time.ctime() + "\n")

    if sorted(denominator_reference_tree_dict.keys()) != sorted(terminal_children_strings_of_ref_denominators.keys()):
        sys.stderr.write("input: " + str(denominator_reference_tree_dict.keys()) + "\n")
        sys.stderr.write("output: " + str(terminal_children_strings_of_ref_denominators.keys()) + "\n")
        sys.exit("ERROR: Not all of the reference trees were parsed!")

    for denominator in sorted(raxml_outfiles.keys()):
        description_text = '# ' + str(text_of_analysis_type[denominator]) + '\n'

        # Retrieve the parsed reference tree from the dictionary of parsed reference trees
        if args.verbose:
            parse_log.write("Retrieving the reference tree for " + denominator + "... ")
        terminal_children_strings_of_reference = terminal_children_strings_of_ref_denominators[denominator]
        if args.verbose:
            parse_log.write("done.\n")

        content_of_previous_labelled_tree_file = ''
        previous_f_contig = ""
        rooted_labelled_trees = ''
        insertion_point_node_hash = ''
        final_assignment_target_strings = Autovivify()

        # Parse all labelled tree files for denominator in parallel
        labelled_tree_files = dict()
        for f_contig in raxml_outfiles[denominator].keys():
            if not os.path.isfile(raxml_outfiles[denominator][f_contig]['labelled_tree']):
                parse_log.write("WARNING: " + str(raxml_outfiles[denominator][f_contig]['labelled_tree']) +
                                "was included in RAxML output files but is not a file. Continuing...\n")
            else:
                labelled_tree_files[f_contig] = raxml_outfiles[denominator][f_contig]['labelled_tree']

        parse_log.write("Parsing the " + str(len(labelled_tree_files.keys())) +
                        " trees for " + denominator + " in parallel... ")
        start_time = time.time()
        parse_log.flush()
        raxml_tree_dict = pparse_raxml_out_trees(labelled_tree_files, args)
        end_time = time.time()
        parse_log.write("done.\n")
        hours, remainder = divmod(end_time - start_time, 3600)
        minutes, seconds = divmod(remainder, 60)
        parse_log.write("Parsing required " + str(hours) + ":" + str(minutes) + ":" + str(seconds) + "\n")
        parse_log.flush()

        for f_contig in sorted(raxml_outfiles[denominator].keys()):
            # Update the progress bar
            acc += 1.0
            if acc >= step_proportion:
                acc -= step_proportion
                time.sleep(0.1)
                sys.stdout.write("-")
                sys.stdout.flush()

            denominator = ''
            if re.search(r'\A(.)', f_contig):
                denominator = f_contig.split('_')[0]
            content_of_labelled_tree_file = ''
            assignments = Autovivify()

            if raxml_option == 'v':
                # Maximum-likelihood analysis
                classification_file = raxml_outfiles[denominator][f_contig]['classification']
                labelled_tree_file = raxml_outfiles[denominator][f_contig]['labelled_tree']
                try:
                    RAxML_labelled_tree = open(labelled_tree_file, 'r')
                except IOError:
                    sys.exit('ERROR: Can\'t open ' + str(labelled_tree_file) + '!\n')

                for line in RAxML_labelled_tree:
                    line = line.strip()
                    content_of_labelled_tree_file += str(line)

                RAxML_labelled_tree.close()
                if not content_of_labelled_tree_file == content_of_previous_labelled_tree_file:
                    parse_log.write("Retrieving the labelled tree " + labelled_tree_file + "... ")
                    parse_log.flush()
                    if f_contig not in raxml_tree_dict.keys():
                        sys.exit("ERROR: " + f_contig + " is not found in not in raxml_tree_dict.keys():"
                                                        " \n" + str(raxml_tree_dict.keys()))
                    rooted_labelled_trees, insertion_point_node_hash = raxml_tree_dict[f_contig]

                    parse_log.write("done.\n")
                    parse_log.flush()
                    final_assignment_target_strings = Autovivify()
                    nr_of_assignments = 0  # This does not exist in the original MLTreeMap perl code
                else:
                    if args.verbose:
                        parse_log.write("Identical RAxML classifications between" + str(f_contig) +
                                        "and" + previous_f_contig + "!")

                new_assignments = Autovivify()
                at_least_one_new_assignment = 0
                try:
                    RAxML_classification = open(classification_file, 'r')
                except IOError:
                    sys.exit('ERROR: Can\'t open ' + str(classification_file) + '!\n')

                for line in RAxML_classification:
                    line = line.strip()
                    query, insertion_point_l, weight = line.split(' ')[0:3]
                    assignment = ''
                    if re.search(r'I(\d+)', insertion_point_l):
                        assignment = re.search(r'I(\d+)', insertion_point_l).group(1)
                    assignments[assignment] = weight
                    if assignment not in final_assignment_target_strings.keys():
                        new_assignments[assignment] = 1
                        at_least_one_new_assignment = 1
                        final_assignment_target_strings[assignment] = ""
                        nr_of_assignments += 1

                RAxML_classification.close()
                if at_least_one_new_assignment > 0:
                    parse_log.write("identifying the terminal children of each assignment for "+f_contig+"... ")
                    parse_log.write(time.ctime() + "\n")
                    parse_log.flush()
                    prae_assignment_target_strings = identify_the_correct_terminal_children_of_each_assignment(
                        terminal_children_strings_of_reference, rooted_labelled_trees, insertion_point_node_hash,
                        new_assignments, args.num_threads, parse_log)
                    parse_log.write("done.\n")

                    for assignment in sorted(prae_assignment_target_strings.keys()):
                        assignment_target_string = prae_assignment_target_strings[assignment]
                        final_assignment_target_strings[assignment] = assignment_target_string

                parse_log.write("Finished parsing " + f_contig + "'s RAxML output at " + time.ctime() + "\n")

            elif raxml_option == 'p':
                # Maximum parsimony analysis
                mp_tree_file = raxml_outfiles[denominator][f_contig]
                assignment = 'mp_root'
                assignments[assignment] = 1
                nr_of_assignments = 1
                prae_assignment_target_strings = get_correct_mp_assignment(terminal_children_strings_of_reference,
                                                                           mp_tree_file, assignments)
                assignment_target_string = prae_assignment_target_strings[assignment]
                final_assignment_target_strings[assignment] = assignment_target_string

            final_RAxML_filename = str(args.output_dir_raxml) + str(f_contig) + '_RAxML_parsed.txt'
            final_raxml_output_files[denominator][final_RAxML_filename] = 1
            
            try:
                output = open(final_RAxML_filename, 'w')
            except IOError:
                sys.exit('ERROR: Can\'t create ' + str(final_RAxML_filename) + '!\n')
            output.write(str(description_text) + '\n')
            
            for assignment in sorted(assignments.keys()):
                assignment_target_string = final_assignment_target_strings[assignment]
                weight = float(assignments[assignment])
                relative_weight = float(weight * 100.0 / float(nr_of_assignments))
                assignment_terminal_targets = assignment_target_string.split(' ')
                nr_of_terminal_targets = len(assignment_terminal_targets) - 1
                output.write('Placement weight ' + '%.2f' % relative_weight + '%: Assignment of query to ')
                if not nr_of_terminal_targets == 1:
                    output.write('the lowest common ancestor of ')
                count = 1

                while count <= nr_of_terminal_targets:
                    assignment_terminal_target = assignment_terminal_targets[count - 1]
                    name_of_terminal_target = tree_numbers_translation[denominator][assignment_terminal_target]
                    try:
                        name_of_terminal_target
                    except NameError:
                        sys.exit('ERROR: ' + str(assignment_terminal_target) +
                                 ' could not be located in the tree with the denominator ' +
                                 str(denominator) + '!\n')
                    output.write(str(name_of_terminal_target) + ' (' + str(assignment_terminal_target) + ')')
                    if count < nr_of_terminal_targets - 1:
                        output.write(', ')
                    if count == nr_of_terminal_targets - 1:
                        output.write(' and ')
                    if count == nr_of_terminal_targets:
                        output.write('.\n')
                    count += 1

            output.close()
            content_of_previous_labelled_tree_file = content_of_labelled_tree_file
            previous_f_contig = f_contig

    sys.stdout.write("-]%\n")
    parse_log.close()
    return final_raxml_output_files


def read_and_understand_the_reference_tree(reference_tree_file, denominator):
    reference_tree_elements = read_the_reference_tree(reference_tree_file)
    reference_tree_info = create_tree_info_hash()
    reference_tree_info = get_node_subtrees(reference_tree_elements, reference_tree_info)
    reference_tree_info = assign_parents_and_children(reference_tree_info, denominator)
    if reference_tree_info is None:
        return denominator, None
    terminal_children_of_reference = build_terminal_children_strings_of_reference_nodes(reference_tree_info)
    return denominator, terminal_children_of_reference


def read_understand_and_reroot_the_labelled_tree(labelled_tree_file, f_contig):
    labelled_tree_elements, insertion_point_node_hash = read_the_raxml_out_tree(labelled_tree_file)
    labelled_tree_info = create_tree_info_hash()
    labelled_tree_info = get_node_subtrees(labelled_tree_elements, labelled_tree_info)
    labelled_tree_info = assign_parents_and_children(labelled_tree_info, f_contig)
    if labelled_tree_info is None:
        return [f_contig, None, insertion_point_node_hash]
    labelled_tree_info = build_tree_info_quartets(labelled_tree_info)
    rooted_labelled_trees = build_newly_rooted_trees(labelled_tree_info)
    return [f_contig, rooted_labelled_trees, insertion_point_node_hash]


def identify_the_correct_terminal_children_of_each_assignment(terminal_children_of_reference,
                                                              rooted_labelled_trees,
                                                              insertion_point_node_hash,
                                                              assignments, num_threads, parse_log):
    terminal_children_of_assignments = build_terminal_children_strings_of_assignments(rooted_labelled_trees,
                                                                                      insertion_point_node_hash,
                                                                                      assignments, num_threads, parse_log)
    real_terminal_children_of_assignments = compare_terminal_children_strings(terminal_children_of_assignments,
                                                                                      terminal_children_of_reference)
    return real_terminal_children_of_assignments


def get_correct_mp_assignment(terminal_children_of_reference, mp_tree_file, assignments):
    potential_terminal_children_strings = read_the_raxml_mp_out_tree(mp_tree_file, assignments)
    real_terminal_children_strings_of_assignments = compare_terminal_children_strings(potential_terminal_children_strings, terminal_children_of_reference)
    return real_terminal_children_strings_of_assignments


def read_the_reference_tree(reference_tree_file):
    try:
        reference_tree = open(reference_tree_file, 'r')
    except IOError:
        sys.exit('ERROR: Could not open ' + reference_tree_file + '!\n')
    tree_string = ''

    for line in reference_tree:
        line = line.strip()
        tree_string += line

    reference_tree.close()

    tree_string = re.sub('\(', 'L', tree_string)
    tree_string = re.sub('\)', 'R', tree_string)
    tree_string = re.sub(r':\d+\.\d+', '', tree_string)
    count = -2

    while re.search('R', tree_string):
        tree_string = re.sub('R', 'Q' + str(count), tree_string, 1)
        count += -1

    tree_string = re.sub(r'Q-\d+;', 'Q;', tree_string)
    tree_string = re.sub('L', '(', tree_string)
    tree_string = re.sub('Q', ')', tree_string)
    reference_tree_elements = split_tree_string(tree_string)
    return reference_tree_elements


def read_the_raxml_out_tree(labelled_tree_file):
    """
    Reads the labelled_tree_file and reformats it for downstream interpretation
    :param labelled_tree_file: RAxML output f_contig.originalRAxML_labelledTree.txt file in various_outputs directory
    :return: An easily interpretable labelled tree
    """

    insertion_point_node_hash = Autovivify()
    try:
        raxml_tree = open(labelled_tree_file, 'r')
    except IOError:
        sys.exit('ERROR: Could not open ' + labelled_tree_file + '!\n')
    tree_string = ''

    for line in raxml_tree:
        line = line.strip()
        tree_string += line

    raxml_tree.close()
    tree_symbols_raw_1 = list(tree_string)
    bracket_diff = 0
    tree_string_neu = '('
    comma_count = 0

    for tree_symbol_raw_1 in tree_symbols_raw_1:
        if comma_count < 2:
            if tree_symbol_raw_1 == '(':
                bracket_diff += 1
            if tree_symbol_raw_1 == ')':
                bracket_diff -= 1
            if tree_symbol_raw_1 == ',' and bracket_diff == 1:
                comma_count += 1
            if comma_count == 2:
                tree_string_neu += '):1.0[I666999666]'
        tree_string_neu += tree_symbol_raw_1

    tree_string = tree_string_neu
    tree_string = re.sub('\(', 'L', tree_string)
    tree_string = re.sub('\)', 'R', tree_string)
    tree_string = re.sub('\[', 'Q', tree_string)

    # Remove the branch lengths
    tree_string = re.sub(":[.0-9]+Q", 'Q', tree_string)

    while re.search(r'((\D(\d+))QI(\d+)])', tree_string):
        to_be_replaced = re.search(r'((\D(\d+))QI(\d+)])', tree_string).group(1)
        replacement = re.search(r'((\D(\d+))QI(\d+)])', tree_string).group(2)
        terminal_leaf = re.search(r'((\D(\d+))QI(\d+)])', tree_string).group(3)
        insertion_point = re.search(r'((\D(\d+))QI(\d+)])', tree_string).group(4)
        if terminal_leaf <= 0:
            sys.exit('ERROR: Your tree has terminal leaves with numbers <= 0. Please change them to positive values!\n')
        insertion_point_node_hash[insertion_point] = terminal_leaf
        tree_string = re.sub(to_be_replaced, replacement, tree_string)
    count = -2

    while re.search(r'QI(\d+)]', tree_string):
        insertion_point_node_hash[re.search(r'QI(\d+)]', tree_string).group(1)] = count
        tree_string = re.sub(r'QI(\d+)]', str(count), tree_string, 1)
        count += -1
    
    tree_string = re.sub('L', '(', tree_string)
    tree_string = re.sub('R', ')', tree_string)
    tree_string = re.sub('Q', '[', tree_string)
    tree_elements = split_tree_string(tree_string)
    return tree_elements, insertion_point_node_hash


def read_the_raxml_mp_out_tree(mp_tree_file, assignments):
    """
    A function for specifically reading the maximum-parsimony tree from RAxML
    :param mp_tree_file: The tree file built by RAxML using the maximum-parsimony based algorithm
    :param assignments: A dictionary for holding nodes -- currently just the root
    :return: Autovivification of all potential terminal children
    """
    potential_terminal_children_strings = Autovivify()
    assignment = ''

    for assig in sorted(assignments.keys()):
        assignment = assig
        sys.stdout.write(assig + "\n")
        sys.stdout.flush()
        break

    try:
        mp_tree = open(mp_tree_file, 'r')
    except IOError:
        sys.exit('ERROR: Can\'t open ' + str(mp_tree_file) + '\n')
    tree_string = ''

    for line in mp_tree:
        line = line.strip()
        tree_string += line

    mp_tree.close()
    tree_string = re.sub('\(', 'L', tree_string)
    tree_string = re.sub('\)', 'R', tree_string)
    if not re.search(r',queryR;\Z', tree_string):
        sys.exit('ERROR: The query is not at the root of ' + str(mp_tree_file) + '!\n')
    else:
        tree_string = re.sub(r',queryR;\Z', 'R;', tree_string)
    tree_string = re.sub(r':\d+\.\d+', '', tree_string)
    count = -2

    while re.search('R', tree_string):
        tree_string = re.sub('R', 'Q' + str(count), tree_string, 1)
        count += -1

    tree_string = re.sub(r'Q-\d+;', 'Q;', tree_string)
    tree_string = re.sub('L', '(', tree_string)
    tree_string = re.sub('Q', ')', tree_string)
    tree_symbols = list(tree_string)
    bracket_diff = 0
    comma_count = 0
    substrings = ['', ',']

    for tree_symbol in tree_symbols:
        if comma_count < 1:
            if tree_symbol == '(':
                bracket_diff += 1
            if tree_symbol == ')':
                bracket_diff += -1
            if tree_symbol == ',' and bracket_diff == 1:
                comma_count += 1
            substrings[0] += tree_symbol
        else:
            substrings[1] += tree_symbol

    for substring in substrings:
        terminal_children = Autovivify()

        for eachGroup in re.findall(r'(\D)(\d+)', str(substring)):
            if eachGroup[0] == '-':
                continue
            terminal_children[eachGroup[1]] = 1

        potential_terminal_children_string = ''

        for potential_terminal_child in sorted(terminal_children.keys(), key=int):
            potential_terminal_children_string += str(potential_terminal_child) + ' '

        potential_terminal_children_strings[assignment][potential_terminal_children_string] = 1

    return potential_terminal_children_strings


def split_tree_string(tree_string):
    tree_symbols_raw = list(str(tree_string))
    count = -1
    previous_symbol = ''
    tree_elements = Autovivify()
    
    for tree_symbol_raw in tree_symbols_raw:
        if re.search(r'\d', tree_symbol_raw) and (re.search(r'\d', previous_symbol) or previous_symbol == '-'):
            tree_elements[count] += tree_symbol_raw
        else:
            count += 1
            tree_elements[count] = tree_symbol_raw
        previous_symbol = tree_symbol_raw
    return tree_elements


def create_tree_info_hash():
    tree_info = Autovivify()
    return tree_info


def get_node_subtrees(tree_elements, tree_info):
    bracket_l_count = 0
    bracket_r_count = 0
    parents_of_node = Autovivify()
    tree_element_nr = -1

    for tree_element in tree_elements.values():
        tree_element_nr += 1
        
        if str(tree_element) == '(':
            bracket_l_count = 1
            bracket_r_count = 0
            tree_sub_element_nr = tree_element_nr
            subtree_string = '('
            
            while True:
                tree_sub_element_nr += 1
                tree_sub_element = tree_elements[tree_sub_element_nr]
                if str(tree_sub_element) == '(':
                    bracket_l_count += 1
                if str(tree_sub_element) == ')':
                    bracket_r_count += 1
                if bracket_l_count == bracket_r_count:
                    node_name = tree_elements[tree_sub_element_nr + 1]
                    if str(node_name) == ';':
                        node_name = -1
                    subtree_string += ')' + str(node_name)
                    tree_info['subtree_of_node'][node_name] = subtree_string
                    break
                else:
                    subtree_string += str(tree_sub_element)
    
    for tree_element in tree_elements.values():
        if not re.search(r'\d+', str(tree_element)):
            continue
        if tree_element in tree_info['subtree_of_node'].keys():
            continue
        tree_info['subtree_of_node'][tree_element] = tree_element
    return tree_info


def assign_parents_and_children(tree_info, source):
    """
    :param tree_info: Autovivification of a tree from get_node_subtrees
    :return: tree info with parent and child relationships included
    """

    # parse_log.write("assigning_parents_and_children... \nStart:\t" + time.ctime() + "\n")
    tree_nodes = sorted(list(tree_info['subtree_of_node'].keys()))
    for node in tree_nodes:
        if node == -1:
            continue
        subtree = str(tree_info['subtree_of_node'][node])
        parent = None
        for potential_parent in tree_nodes:
            if node == potential_parent:
                continue
            potential_parent_subtree = str(tree_info['subtree_of_node'][potential_parent])
            subtree = re.sub('\(', 'L', subtree)
            subtree = re.sub('\)', '#', subtree)
            potential_parent_subtree = re.sub('\(', 'L', potential_parent_subtree)
            potential_parent_subtree = re.sub('\)', '#', potential_parent_subtree)
            potential_parent = str(potential_parent)
            if re.search(r'\AL'+re.escape(subtree)+r',.+#'+re.escape(potential_parent)+r'\Z', potential_parent_subtree) or \
               re.search(r'\AL.+,'+re.escape(subtree)+r'#'+re.escape(potential_parent)+r'\Z', potential_parent_subtree):
                parent = potential_parent
                break
        if parent is None:
            sys.stderr.write("ERROR: No parent assigned for " + node + " for " + source + "\n")
            sys.stderr.write("This is due to either an incompatibility with your RAxML version or ")
            sys.stderr.write("a formatting error in the reference tree.\n")
            sys.stderr.write("Please post an issue on the github page!\n")
            sys.stderr.flush()
            return None
            # TODO: handle this better when dealing with multiple processes
        else:
            tree_info['parent_of_node'][node] = parent
            tree_info['children_of_node'][parent][node] = 1

    # parse_log.write("End:  \t" + time.ctime() + "\n")
    return tree_info


def build_tree_info_quartets(tree_info):
    for node in sorted(tree_info['parent_of_node'].keys(), key=int):
        parent = tree_info['parent_of_node'][node]
        if int(parent) == -1:
            for roots_child in sorted(tree_info['children_of_node']['-1'].keys(), key=int):
                if roots_child == node:
                    continue
                parent = roots_child

        tree_info['quartets'][node][parent] = 1
        if node in tree_info['children_of_node']:
            for child in sorted(tree_info['children_of_node'][node].keys(), key=int):
                tree_info['quartets'][node][child] = 1

    return tree_info


def build_newly_rooted_trees(tree_info):
    tree_number = 0
    list_of_already_used_attachments = Autovivify()
    rooted_trees = Autovivify()
    
    for node in sorted(tree_info['quartets'].keys(), key=int):
        if node in list_of_already_used_attachments:
            continue
        for attachment in sorted(tree_info['quartets'][node].keys(), key=int):
            list_of_already_used_attachments[attachment] = 1
            tree_string = ''
            node_infos = Autovivify()
            node_infos['previous_node'] = ''
            node_infos['node'] = ';'
            node_infos['open_attachments'][node] = 1
            node_infos['open_attachments'][attachment] = 1
            new_tree = recursive_tree_builder(tree_info, node_infos, tree_string)
            rooted_trees[tree_number] = new_tree
            tree_number += 1
    
    return rooted_trees


def recursive_tree_builder(tree_info, node_infos, tree_string):
    node = node_infos['node']
    count = 0

    for attachment in sorted(node_infos['open_attachments'].keys(), key=int):
        count += 1
        if count == 1:
            tree_string += '('
        node_infos2 = Autovivify()
        node_infos2['previous_node'] = node
        node_infos2['node'] = attachment
        count2 = 0

        for attachment_of_used_attachment in sorted(tree_info['quartets'][attachment].keys()):
            if attachment_of_used_attachment in node_infos['open_attachments']:
                continue
            if attachment_of_used_attachment == node:
                continue
            count2 += 1
            node_infos2['open_attachments'][attachment_of_used_attachment] = 1

        if count2 > 0:
            tree_string = recursive_tree_builder(tree_info, node_infos2, tree_string)
        else:
            tree_string += str(attachment)
        if count == 1:
            tree_string += ','
        if count == 2:
            tree_string += ')' + str(node)

    return tree_string


def parallel_subtree_node_retriever(rooted_trees, num_threads, parse_log):
    """
    Run `get_node_subtrees` in parallel for each of the elements in rooted_trees
    :param rooted_trees: Dictionary of rooted trees
    :param num_threads: Number of threads to use
    :param parse_log: The file object to write parsing information to
    :return: rooted_tree_nodes - a list of results from get_node_subtrees()
    """
    job_queue = JoinableQueue()
    result_queue = Queue()
    rooted_tree_nodes = list()

    worker_group = [NodeRetrieverWorker(job_queue, result_queue) for i in range(int(num_threads))]
    for worker in worker_group:
        worker.start()

    tasks = [split_tree_string(rooted_trees[rooted_tree]) for rooted_tree in rooted_trees.keys()]
    parse_log.write("Number of subtrees = " + str(len(tasks)) + "\n")
    parse_log.flush()
    for rooted_tree_elements in tasks:
        if job_queue.full():
            sys.exit("ERROR: multiprocessing.Queue full in parallel_subtree_node_retriever!")
        job_queue.put(rooted_tree_elements)

    for i in range(int(num_threads)):
        job_queue.put(None)

    for i in range(len(rooted_trees.keys())):
        rooted_tree_nodes.append(result_queue.get())

    job_queue.close()
    result_queue.close()
    result_queue.join_thread()

    return rooted_tree_nodes


def build_terminal_children_strings_of_assignments(rooted_trees, insertion_point_node_hash,
                                                   assignments, num_threads, parse_log):
    terminal_children_strings_of_assignments = Autovivify()

    for assignment in sorted(assignments.keys()):
        internal_node_of_assignment = insertion_point_node_hash[assignment]

        rooted_tree_nodes = parallel_subtree_node_retriever(rooted_trees, num_threads, parse_log)

        for rooted_tree_info in rooted_tree_nodes:
            assignment_subtree = str(rooted_tree_info['subtree_of_node'][str(internal_node_of_assignment)])
            terminal_children = Autovivify()
            if re.search(r'\A(\d+)\Z', assignment_subtree):
                terminal_children[re.search(r'\A(\d+)\Z', assignment_subtree).group(1)] = 1
            else:

                for each_hit in re.findall(r'(\D)(\d+)', assignment_subtree):
                    if each_hit[0] == '-':
                        continue
                    terminal_children[each_hit[1]] = 1

            terminal_children_string_of_assignment = ''

            for terminal_child_of_assignment in sorted(terminal_children.keys(), key=int):
                terminal_children_string_of_assignment += str(terminal_child_of_assignment) + ' '

            terminal_children_strings_of_assignments[assignment][terminal_children_string_of_assignment] = 1

    return terminal_children_strings_of_assignments


def build_terminal_children_strings_of_reference_nodes(reference_tree_info):
    terminal_children_strings_of_reference = Autovivify()

    for node in sorted(reference_tree_info['subtree_of_node'].keys()):
        reference_subtree = reference_tree_info['subtree_of_node'][node]
        terminal_children = Autovivify()
        if re.search(r'\A(\d+)\Z', str(reference_subtree)):
            terminal_children[re.search(r'\A(\d+)\Z', str(reference_subtree)).group(1)] = 1
        else:

            for each_hit in re.findall(r'(.)(\d+)', str(reference_subtree)):
                if each_hit[0] == '-':
                    continue
                terminal_children[each_hit[1]] = 1

        terminal_children_string_of_reference = ''

        for terminal_child_of_reference in sorted(terminal_children.keys(), key=int):
            terminal_children_string_of_reference += str(terminal_child_of_reference) + ' '

        terminal_children_strings_of_reference[terminal_children_string_of_reference] = 1

    return terminal_children_strings_of_reference


def compare_terminal_children_strings(terminal_children_strings_of_assignments, terminal_children_strings_of_reference):
    real_terminal_children_strings_of_assignments = Autovivify()
    there_was_a_hit = 0

    for assignment in sorted(terminal_children_strings_of_assignments.keys()):
        real_terminal_children_string = ''

        for terminal_children_string_of_assignment in sorted(terminal_children_strings_of_assignments[assignment].keys()):
            if terminal_children_string_of_assignment in terminal_children_strings_of_reference:
                real_terminal_children_string = terminal_children_string_of_assignment
                real_terminal_children_strings_of_assignments[assignment] = real_terminal_children_string
                there_was_a_hit = 1
                break

        if str(real_terminal_children_string) == '' and not str(assignment) == 'mp_root':
            sys.exit('ERROR: The RAxML output tree could not be rooted correctly!!!\n')

    if there_was_a_hit <= 0:
        sys.exit('ERROR: The RAxML output tree could not be rooted correctly!!!\n')
    return real_terminal_children_strings_of_assignments


def concatenate_RAxML_output_files(args, final_raxml_output_files, text_of_analysis_type):
    if args.verbose:
        sys.stdout.write("Concatenating the RAxML outputs for each marker gene class...\n")
    output_directory_final = args.output_dir_final
    
    for denominator in sorted(final_raxml_output_files.keys()):
        nr_of_files = 0
        assignments = Autovivify()
        description_text = '# ' + str(text_of_analysis_type[denominator]) + '\n'
        final_output_file_name = str(output_directory_final) + str(denominator) + '_concatenated_RAxML_outputs.txt'
        
        for final_RAxML_output_file in sorted(final_raxml_output_files[denominator].keys()):
            nr_of_files += 1
            try:
                final_raxml_output = open(final_RAxML_output_file, 'r')
            except IOError:
                sys.exit('ERROR: Can\'t open ' + str(final_RAxML_output_file) + '!\n')
            
            for line in final_raxml_output:
                line = line.strip()
                if re.search(r'Placement weight (\d+\.\d+)%: (.+)\Z', line):
                    weight = float(re.search(r'Placement weight (\d+\.\d+)%: (.+)\Z', line).group(1))
                    assignment = re.search(r'Placement weight (\d+\.\d+)%: (.+)\Z', line).group(2)
                    if assignment in assignments.keys():
                        assignments[assignment] += weight
                    else:
                        assignments[assignment] = weight
                else:
                    continue

            final_raxml_output.close()

        assignments_with_relative_weights = Autovivify()

        for assignment in sorted(assignments.keys(), reverse=True):
            weight = assignments[assignment]
            relative_weight = weight / float(nr_of_files)
            assignments_with_relative_weights[relative_weight][assignment] = 1

        try:
            output = open(final_output_file_name, 'w')
        except IOError:
            sys.exit('ERROR: Can\'t create ' + str(final_output_file_name) + '!\n')
        if args.verbose:
            sys.stdout.write(str(denominator) + '_ results concatenated:\n')
        output.write(str(description_text) + '\n')
        sum_of_relative_weights = 0

        for relative_weight in sorted(assignments_with_relative_weights.keys(), reverse=True):

            for assignment in sorted(assignments_with_relative_weights[relative_weight].keys(), reverse=True):
                sum_of_relative_weights += relative_weight
                sys.stdout.write('Placement weight ')
                sys.stdout.write('%.2f' % relative_weight + "%: ")
                sys.stdout.write(assignment + "\n")
                output.write('Placement weight ' + str(relative_weight) + '%: ' + str(assignment) + '\n')

        output.close()
        sys.stdout.write('_' + str(denominator) + '_ sum of placement weights (should be 100): ')
        sys.stdout.write(str(int(sum_of_relative_weights + 0.5)) + "\n")
        sys.stdout.flush()


def read_species_translation_files(args, cog_list):
    """
    :param cog_list: The list on COGs used for tre insertion
    :return: The taxonomic identifiers for each of the organisms in a tree for all trees
    """

    tree_numbers_translation = Autovivify()
    translation_files = Autovivify()
    phylogenetic_denominator = args.reftree
    if phylogenetic_denominator == 'g':
        translation_files[phylogenetic_denominator] = args.mltreemap + os.sep + \
                                                      'data' + os.sep + 'tree_data' + \
                                                      os.sep + 'tax_ids_geba_tree.txt'
    elif phylogenetic_denominator == 'i':
        translation_files[phylogenetic_denominator] = args.mltreemap + os.sep + \
                                                      'data' + os.sep + 'tree_data' + \
                                                      os.sep + 'tax_ids_fungitr.txt'
    elif phylogenetic_denominator == 'p':
        translation_files[phylogenetic_denominator] = args.mltreemap + os.sep + \
                                                      'data' + os.sep + 'tree_data' + \
                                                      os.sep + 'tax_ids_nr.txt'

    for functional_cog in sorted(cog_list['functional_cogs'].keys()):
        denominator = cog_list['functional_cogs'][functional_cog]
        filename = 'tax_ids_' + str(functional_cog) + '.txt'
        translation_files[denominator] = args.mltreemap + os.sep + \
                                         'data' + os.sep + 'tree_data' + \
                                         os.sep + filename

    for phylogenetic_rRNA_cog in sorted(cog_list['phylogenetic_rRNA_cogs'].keys()):
        denominator = cog_list['phylogenetic_rRNA_cogs'][phylogenetic_rRNA_cog]
        filename = 'tax_ids_' + str(phylogenetic_rRNA_cog) + '.txt'
        translation_files[denominator] = args.mltreemap + os.sep + \
                                         'data' + os.sep + 'tree_data' + \
                                         os.sep + filename

    for denominator in sorted(translation_files.keys()):
        filename = translation_files[denominator]
        try:
            if args.py_version == 3:
                cog_tax_ids = open(filename, 'r', encoding='utf-8')
            else:
                cog_tax_ids = open(filename, 'r')
        except IOError:
            sys.exit('ERROR: Can\'t open ' + str(filename) + '!\n')

        for line in cog_tax_ids:
            line = line.strip()
            try:
                number, translation = line.split('\t')
            except ValueError:
                sys.exit('ValueError: .split(\'\\t\') on ' + str(line))
            tree_numbers_translation[denominator][number] = translation

        cog_tax_ids.close()

    return tree_numbers_translation


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

    # jython
    try:
        from java.lang import Runtime
        runtime = Runtime.getRuntime()
        res = runtime.availableProcessors()
        if res > 0:
            return res
    except ImportError:
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


def delete_files(args):
    sys.stdout.write('Deleting files as requested\n')
    sections_to_be_deleted = []
    if args.delete:
        sections_to_be_deleted = args.delete.split(':')

    files_to_be_deleted = []

    for section in sections_to_be_deleted:
        if section == '1':
            files_to_be_deleted += glob.glob(args.output_dir_var + '*.fa')
            files_to_be_deleted += glob.glob(args.output_dir_var + '*sequence.txt')
            files_to_be_deleted += glob.glob(args.output_dir_var + '*sequence_shortened.txt')
            files_to_be_deleted += glob.glob(args.output_dir_var + '*.fasta_formatted.txt')
        if section == '2':
            files_to_be_deleted += glob.glob(args.output_dir_var + '*BLAST_results*')
            files_to_be_deleted += glob.glob(args.output_dir_var + '*blast_result_purified.txt')
            files_to_be_deleted += glob.glob(args.output_dir_var + '*rRNA_result_summary.txt')
        if section == '3':
            files_to_be_deleted += glob.glob(args.output_dir_var + '*genewise.txt')
            files_to_be_deleted += glob.glob(args.output_dir_var + '*genewise_result_summary.txt')
        if section == '4':
            files_to_be_deleted += glob.glob(args.output_dir_var + '*.mfa')
            files_to_be_deleted += glob.glob(args.output_dir_var + '*.mfa-gb')
            files_to_be_deleted += glob.glob(args.output_dir_var + '*.mfa-gb.txt')
        if section == '5':
            files_to_be_deleted += glob.glob(args.output_dir_var + '*_exit_after_Gblocks.txt')
            files_to_be_deleted += glob.glob(args.output_dir_var + '*_RAxML.txt')
            files_to_be_deleted += glob.glob(args.output_dir_var + '*RAxML_classification.txt')
            files_to_be_deleted += glob.glob(args.output_dir_var + '*RAxML_info.txt')
            files_to_be_deleted += glob.glob(args.output_dir_var + '*RAxML_labelledTree.txt')
            files_to_be_deleted += glob.glob(args.output_dir_var + '*.phy')
            files_to_be_deleted += glob.glob(args.output_dir_var + '*.phy.reduced')

    for file in files_to_be_deleted:
        if path.exists(file):
            os.remove(file)


def single_family_msa(args, cog_list, formatted_fasta_dict):
    """
    A wrapper function for hmmalign -- to generate a multiple-sequence alignment with the reference sequences
    of the gene family being updated
    :param args: Command-line argument object returned by argparse
    :param cog_list: The reference gene family to be updated
    :param formatted_fasta_dict: keys are fasta headers; values are fasta sequences. Returned by format_read_fasta
    :return: An Autovivification mapping the summary files to each contig
    """
    reference_data_prefix = args.reference_data_prefix
    hmmalign_singlehit_files = Autovivify()
    if args.verbose:
        sys.stdout.write("Running hmmalign... ")

    cog = cog_list["all_cogs"].keys()[0]
    denominator = cog_list["all_cogs"][cog]

    start = 0

    # Imitate the Genewise / blastpSummaryFiles output
    for contig in formatted_fasta_dict.keys():
        header = contig[1:]
        sequence = formatted_fasta_dict[contig]
        end = len(sequence)

        f_contig = denominator + "_" + header
        genewise_singlehit_file = args.output_dir_var + os.sep + \
                                  f_contig + '_' + cog + "_" + str(start) + "_" + str(end)
        hmmalign_singlehit_files[f_contig][genewise_singlehit_file + ".mfa"] = True
        genewise_singlehit_file_fa = genewise_singlehit_file + ".fa"
        try:
            outfile = open(genewise_singlehit_file_fa, 'w')
            fprintf(outfile, '>query\n%s\n', sequence)
            outfile.close()
        except IOError:
            sys.stderr.write('Can\'t create ' + genewise_singlehit_file_fa + '\n')
            sys.exit(0)
        mltreemap_resources = args.mltreemap + os.sep + 'data' + os.sep
        hmmalign_command = [args.executables["hmmalign"], '-m', '--mapali',
                            mltreemap_resources + reference_data_prefix + 'alignment_data' +
                            os.sep + cog + '.fa',
                            '--outformat', 'Clustal',
                            mltreemap_resources + reference_data_prefix + 'hmm_data' + os.sep + cog + '.hmm',
                            genewise_singlehit_file_fa, '>', genewise_singlehit_file + '.mfa']
        os.system(' '.join(hmmalign_command))

    if args.verbose:
        sys.stdout.write("done.\n")
    return hmmalign_singlehit_files


def execute_gblocks(args, aligned_fasta):
    data_size = retrieve_data_size(aligned_fasta)
    min_flank_pos = str(0.55 * data_size)
    gblock_command = args.executables["Gblocks"] + \
        "%s -t=p -s=y -u=n -p=t -b3=15 -b4=3 -b5=h -b2=%s" % (aligned_fasta, min_flank_pos)

    print gblock_command, "\n"

    os.system(gblock_command)


def get_new_ref_sequences(update_tree):
    """
    Function for retrieving the protein sequences from the MLTreeMap various_outputs
    :param update_tree: An instance of CreateFuncTreeUtility class
    :return: aa_dictionary is a dictionary of fasta sequences with headers as keys and protein sequences as values
    """
    aa_dictionary = dict()
    various_files = os.listdir(update_tree.InputData + os.sep + "various_outputs" + os.sep)

    for var_file in various_files:
        var_file_path = update_tree.InputData + os.sep + "various_outputs" + os.sep + var_file
        prefix = var_file.split('_')[0]
        if prefix == update_tree.Denominator:
            # Pull out the contig ID between the denominator and the COG ID for those that are in final_RAxML_outputs/
            orf_fa = re.match("%s_(\S+)_%s(_\d+)*.fa" % (update_tree.Denominator, update_tree.COG), var_file)
            if orf_fa:
                suffix = re.sub("%s_" % update_tree.Denominator, '', var_file)
                seq_name = re.sub("_%s(_\d+)*.fa" % update_tree.COG, '', suffix)
                if seq_name in update_tree.names:
                    line_counter = 0
                    aa_dictionary['>' + seq_name] = ""
                    try:
                        fasta = open(var_file_path, 'r')
                    except:
                        raise IOError("Unable to open " + var_file_path + " for reading!")

                    line = fasta.readline()
                    if line[0] != '>':
                        sys.stderr.write("ERROR: first line in " + var_file_path + " is not a proper FASTA file!")
                        sys.stderr.flush()
                        sys.exit()
                    # TODO: Include line to quality-check the sequences prior to saving in aa_dictionary
                    while line:
                        line_counter += 1
                        if line[0] == '>':
                            pass
                        else:
                            aa_dictionary['>' + seq_name] = line.strip()
                        line = fasta.readline()
                    fasta.close()

                    if line_counter > 2:
                        sys.stderr.write("ERROR: " + var_file_path + " contains more than 1 sequence when 1 is expected!")
                        sys.stderr.flush()
                        sys.exit()

    return aa_dictionary


def cluster_new_reference_sequences(update_tree, args, new_ref_seqs_fasta):
    if args.verbose:
        sys.stdout.write("Running usearch to cluster sequences at %s percent identity... " % str(args.uclust_identity))
        sys.stdout.flush()

    usearch_command = [args.executables["usearch"]]
    usearch_command += ["-sortbylength", new_ref_seqs_fasta]
    usearch_command += ["--output", update_tree.Output + "usearch_sorted.fasta"]
    usearch_command += ["--log", update_tree.Output + os.sep + "usearch_sort.log"]
    usearch_command += ["1>", "/dev/null", "2>", "/dev/null"]

    p_usort = subprocess.Popen(' '.join(usearch_command), shell=True, preexec_fn=os.setsid)
    p_usort.wait()
    if p_usort.returncode != 0:
        sys.stderr.write("ERROR: usearch did not complete successfully for:\n")
        sys.stderr.write(str(' '.join(usearch_command)))
        sys.stderr.flush()

    uclust_command = [args.executables["usearch"]]
    uclust_command += ["-cluster_fast", update_tree.Output + "usearch_sorted.fasta"]
    uclust_command += ["--id", str(args.uclust_identity)]
    uclust_command += ["--centroids", update_tree.Output + "uclust_" + update_tree.COG + ".fasta"]
    uclust_command += ["--uc", update_tree.Output + "uclust_" + update_tree.COG + ".uc"]
    uclust_command += ["--log", update_tree.Output + os.sep + "usearch_cluster.log"]
    uclust_command += ["1>", "/dev/null", "2>", "/dev/null"]

    p_uclust = subprocess.Popen(' '.join(uclust_command), shell=True, preexec_fn=os.setsid)
    p_uclust.wait()
    if p_uclust.returncode != 0:
        sys.stderr.write("ERROR: usearch did not complete successfully for:\n")
        sys.stderr.write(str(' '.join(uclust_command)))
        sys.stderr.flush()

    if args.verbose:
        sys.stdout.write("done.\n")
        sys.stdout.flush()

    return


def swap_tree_names(tree, tree_swap_name, random_map, ref_tax_id_map):
    """
    Function used for replacing unique identifiers in a NEWICK tree file
    :param tree:
    :param tree_swap_name:
    :param random_map:
    :param ref_tax_id_map:
    :return:
    """
    try:
        old_tree = open(tree, 'r')
    except:
        raise IOError("Unable to open " + tree + " for reading!")
    try:
        new_tree = open(tree_swap_name, 'w')
    except:
        raise IOError("Unable to open " + tree_swap_name + " for writing!")

    newick_tree = old_tree.readlines()
    old_tree.close()

    if len(newick_tree) > 1:
        raise AssertionError("ERROR: " + tree + " should only contain one line of text to be a NEWICK tree!")
    else:
        newick_tree = str(newick_tree[0])

    for node_id in random_map:
        newick_tree = re.sub(random_map[node_id], ref_tax_id_map[node_id], newick_tree)

    new_tree.write(newick_tree + "\n")
    new_tree.close()

    return


def filter_short_sequences(query_fasta, length_threshold=110):
    """
    Writes all sequences longer than length_threshold to a new FASTA file before performing MSA
    :param query_fasta: FASTA file which may contain sequences shorter than length_threshold
    :param length_threshold: Minimum number of AA a sequence must contain to be included in further analyses
    :return: name of a new FASTA file with only the long sequences
    """
    try:
        fasta_in = open(query_fasta, 'r')
    except:
        raise IOError("Unable to open " + query_fasta + " for reading!")

    new_faa = '.'.join(query_fasta.split('.')[0:-1]) + "_" + str(length_threshold) + ".faa"
    try:
        fasta_out = open(new_faa, 'w')
    except:
        raise IOError("Unable to open " + query_fasta + " for reading!")

    line = fasta_in.readline()
    header = ""
    sequence = ""
    while line:
        if line[0] == '>':
            if len(sequence) > length_threshold:
                fasta_out.write(header + "\n")
                fasta_out.write(sequence + "\n")
            header = line.strip()
            sequence = ""
        else:
            sequence += line.strip()
        line = fasta_in.readline()

    fasta_in.close()
    fasta_out.close()

    return new_faa


def align_reads_to_nucs(args):
    """
    Align the BLAST-predicted ORFs to the reads using BWA MEM
    :param args:
    :return: Path to the SAM file
    """
    input_multi_fasta = re.match(r'\A.*\/(.*)', args.input).group(1)
    orf_nuc_fasta = args.output_dir_var + '.'.join(input_multi_fasta.split('.')[:-1]) + "_genes.fna"
    rpkm_output_dir = args.output + "RPKM_outputs" + os.sep
    if not os.path.exists(rpkm_output_dir):
        try:
            os.makedirs(rpkm_output_dir)
        except:
            raise IOError("Unable to make " + rpkm_output_dir)

    if args.verbose:
        sys.stdout.write("Aligning reads to ORFs with BWA MEM... ")
        sys.stdout.flush()

    sam_file = rpkm_output_dir + '.'.join(os.path.basename(orf_nuc_fasta).split('.')[0:-1]) + ".sam"
    index_command = [args.executables["bwa"], "index"]
    index_command += [orf_nuc_fasta]
    index_command += ["1>", "/dev/null", "2>", args.output + "mltreemap_bwa_index.stderr"]

    p_index = subprocess.Popen(' '.join(index_command), shell=True, preexec_fn=os.setsid)
    p_index.wait()
    if p_index.returncode != 0:
        sys.stderr.write("ERROR: bwa index did not complete successfully for:\n")
        sys.stderr.write(str(' '.join(index_command)) + "\n")
        sys.stderr.flush()

    bwa_command = [args.executables["bwa"], "mem"]
    bwa_command += ["-t", str(args.num_threads)]
    if args.pairing == "pe" and not args.reverse:
        bwa_command.append("-p")
        sys.stderr.write("FASTQ file containing reverse mates was not provided - assuming the reads are interleaved!\n")
        sys.stderr.flush()
    elif args.pairing == "se":
        bwa_command += ["-S", "-P"]

    bwa_command.append(orf_nuc_fasta)
    bwa_command.append(args.reads)
    if args.pairing == "pe" and args.reverse:
        bwa_command.append(args.reverse)
    bwa_command += ["1>", sam_file, "2>", args.output + "mltreemap_bwa_mem.stderr"]

    p_bwa = subprocess.Popen(' '.join(bwa_command), shell=True, preexec_fn=os.setsid)
    p_bwa.wait()
    if p_bwa.returncode != 0:
        sys.stderr.write("ERROR: bwa mem did not complete successfully for:\n")
        sys.stderr.write(str(' '.join(bwa_command)) + "\n")
        sys.stderr.flush()

    if args.verbose:
        sys.stdout.write("done.\n")
        sys.stdout.flush()

    return sam_file, orf_nuc_fasta


def run_rpkm(args, sam_file, orf_nuc_fasta):
    """
    Calculate RPKM values using the rpkm executable
    :param args:
    :param sam_file:
    :param orf_nuc_fasta:
    :return: Path to the RPKM output csv file
    """
    if args.verbose:
        sys.stdout.write("Calculating RPKM values for each ORF... ")
        sys.stdout.flush()

    rpkm_output_file = '.'.join(sam_file.split('.')[0:-1]) + ".csv"
    rpkm_output_dir = args.output + "RPKM_outputs" + os.sep

    rpkm_command = [args.executables["rpkm"]]
    rpkm_command += ["-c", orf_nuc_fasta]
    rpkm_command += ["-a", sam_file]
    rpkm_command += ["-o", rpkm_output_file]
    rpkm_command += ["1>", rpkm_output_dir + "rpkm_stdout.txt", "2>", rpkm_output_dir + "rpkm_stderr.txt"]

    p_rpkm = subprocess.Popen(' '.join(rpkm_command), shell=True, preexec_fn=os.setsid)
    p_rpkm.wait()
    if p_rpkm.returncode != 0:
        sys.stderr.write("ERROR: RPKM calculation did not complete successfully for:\n")
        sys.stderr.write(str(' '.join(rpkm_command)) + "\n")
        sys.stderr.flush()

    if args.verbose:
        sys.stdout.write("done.\n")
        sys.stdout.flush()

    return rpkm_output_file


def main(argv):
    # STAGE 1: Prompt the user and prepare files and lists for the pipeline
    parser = get_options()
    args = check_parser_arguments(parser)
    args = check_previous_output(args)
    cog_list, text_of_analysis_type = create_cog_list(args)
    non_wag_cog_list = get_non_wag_cogs(args)
    if args.check_trees:
        validate_inputs(args, cog_list)

    if args.skip == 'n':
        # Read and format the input FASTA
        formatted_fasta_dict = format_read_fasta(args)
        if re.match(r'\A.*\/(.*)', args.input):
            input_multi_fasta = os.path.basename(args.input)
        else:
            input_multi_fasta = args.input
        args.formatted_input_file = args.output_dir_var + input_multi_fasta + "_formatted.fasta"
        formatted_fasta_files = write_new_fasta(formatted_fasta_dict, args.formatted_input_file)
        # UPDATE GENE FAMILY TREE MODE:
        if args.reftree not in ['i', 'g', 'p']:
            cog_list, text_of_analysis_type = single_cog_list(args.reftree, cog_list, text_of_analysis_type)
            hmmalign_singlehit_files = single_family_msa(args, cog_list, formatted_fasta_dict)
        else:
            # STAGE 2: Run BLAST to determine which COGs are present in the input sequence(s)
            run_blast(args, formatted_fasta_files, cog_list)
            raw_blast_results = collect_blast_outputs(args)
            blast_hits_purified = parse_blast_results(args, raw_blast_results, cog_list)

            # STAGE 3: Produce amino acid sequences based on the COGs found in the input sequence(s)
            contig_coordinates, shortened_sequence_files, gene_coordinates = make_genewise_inputs(args,
                                                                                                  blast_hits_purified,
                                                                                                  formatted_fasta_dict)
            write_nuc_sequences(args, gene_coordinates, formatted_fasta_dict)
            formatted_fasta_dict.clear()
            genewise_summary_files = Autovivify()
            if args.molecule == 'n':
                genewise_outputfiles = start_genewise(args, shortened_sequence_files, blast_hits_purified)
                genewise_summary_files = parse_genewise_results(args, genewise_outputfiles, contig_coordinates)
                get_rRNA_hit_sequences(args, blast_hits_purified, cog_list, genewise_summary_files)
            elif args.molecule == 'a':
                genewise_summary_files = blastpParser(args, blast_hits_purified)

            # STAGE 4: Run hmmalign and Gblocks to produce the MSAs required to perform the subsequent ML/MP estimations
            hmmalign_singlehit_files = prepare_and_run_hmmalign(args, genewise_summary_files, cog_list)
        concatenated_mfa_files, nrs_of_sequences, models_to_be_used = concatenate_hmmalign_singlehits_files(args,
                                                                                                            hmmalign_singlehit_files,
                                                                                                            non_wag_cog_list)
        gblocks_files = start_gblocks(args, concatenated_mfa_files, nrs_of_sequences)
        phy_files = produce_phy_file(args, gblocks_files, nrs_of_sequences)

        # STAGE 5: Run RAxML to compute the ML/MP estimations
        raxml_outfiles, denominator_reference_tree_dict, num_raxml_outputs = start_RAxML(args, phy_files,
                                                                                         cog_list, models_to_be_used)
        tree_numbers_translation = read_species_translation_files(args, cog_list)
        final_raxml_output_files = parse_RAxML_output(args, denominator_reference_tree_dict, tree_numbers_translation,
                                                      raxml_outfiles, text_of_analysis_type, num_raxml_outputs)
        concatenate_RAxML_output_files(args, final_raxml_output_files, text_of_analysis_type)

    if args.rpkm:
        sam_file, orf_nuc_fasta = align_reads_to_nucs(args)
        run_rpkm(args, sam_file, orf_nuc_fasta)
        # normalize_rpkm_values()

    # TODO: Provide stats file with proportion of sequences detected to have marker genes, N50, map contigs to genes
    # STAGE 6: Optionally update the reference tree
    if args.update_tree:
        # TODO: Include a minimum length threshold for a sequence to be considered a new reference sequence
        update_tree = CreateFuncTreeUtility(args.output, args.reftree)
        update_tree.find_cog_name(cog_list)
        update_tree.get_contigs_for_ref()
        aa_dictionary = get_new_ref_sequences(update_tree)
        new_ref_seqs_fasta = update_tree.Output + os.path.basename(update_tree.InputData) +\
                             "_" + update_tree.COG + "_unaligned.fasta"
        write_new_fasta(aa_dictionary, new_ref_seqs_fasta)
        ref_tax_id_map = update_tree.write_reference_names()
        for seq_name in update_tree.names:
            ref_tax_id_map[seq_name] = seq_name
        if args.uclust and len(aa_dictionary) > 1:
            cluster_new_reference_sequences(update_tree, args, new_ref_seqs_fasta)
            query_fasta = update_tree.Output + "uclust_" + update_tree.COG + ".fasta"
            # TODO: Make sure the tree is updated only if there are novel sequences (i.e. <97% similar to ref sequences)
        else:
            if len(aa_dictionary) == 1 and args.uclust:
                sys.stderr.write("WARNING: Not clustering new " + update_tree.COG + " since there is 1 sequence\n")
                sys.stderr.flush()
            query_fasta = new_ref_seqs_fasta

        # Remove short sequences
        # TODO: Make the length_threshold for filter_short_sequences adaptable to reference input
        query_fasta = filter_short_sequences(query_fasta, 150)

        ref_align = "data/alignment_data/" + update_tree.COG + ".fa"
        aligned_fasta = update_tree.align_sequences(args.alignment_mode, ref_align, query_fasta, args)
        fasta_random, original_random_dict = update_tree.randomize_fasta_id(aligned_fasta)
        update_tree.create_random_names(original_random_dict, ref_tax_id_map)

        if args.gap_removal == "y":
            if args.verbose:
                sys.stdout.write("Executing Gblocks... ")
                sys.stdout.flush()
            execute_gblocks(args, fasta_random)
            if args.verbose:
                sys.stdout.write("done.\n")
                sys.stdout.flush()
            os.system('cp %s-gb %s' % (fasta_random, fasta_random))

        os.system('java -cp sub_binaries/readseq.jar run -a -f=12 %s' % fasta_random)

        phylip_file = update_tree.Output + "%s.phy" % update_tree.COG
        os.system('mv %s.phylip %s' % (fasta_random, phylip_file))

        time_of_run = strftime("%d_%b_%Y_%H_%M", gmtime())
        project_folder = update_tree.Output + str(time_of_run) + os.sep
        os.makedirs(project_folder)
        raxml_destination_folder = project_folder + "phy_files_%s" % update_tree.COG
        final_tree_dir = project_folder + "final_tree_files" + os.sep
        alignment_files_dir = project_folder + "alignment_files" + os.sep

        if args.verbose:
            sys.stdout.write("Executing RAxML... ")
            sys.stdout.flush()
        update_tree.execute_raxml(phylip_file, raxml_destination_folder, args)
        if args.verbose:
            sys.stdout.write("done.\n")
            sys.stdout.flush()

        # Organize Output Files #

        os.makedirs(final_tree_dir)
        os.makedirs(alignment_files_dir)

        shutil.move(aligned_fasta, alignment_files_dir)
        shutil.move(phylip_file, alignment_files_dir)

        best_tree = raxml_destination_folder + "/RAxML_bestTree." + update_tree.COG
        bootstrap_tree = raxml_destination_folder + "/RAxML_bipartitions." + update_tree.COG
        best_tree_nameswap = final_tree_dir + update_tree.COG + "_best.tree"
        bootstrap_nameswap = final_tree_dir + update_tree.COG + "_bootstrap.tree"
        swap_tree_names(best_tree, best_tree_nameswap, original_random_dict, ref_tax_id_map)
        swap_tree_names(bootstrap_tree, bootstrap_nameswap, original_random_dict, ref_tax_id_map)

        prefix = update_tree.Output + update_tree.COG
        os.system('mv %s* %s' % (prefix, project_folder))

        if args.uclust == "y":
            os.system('mkdir %s_uclust' % update_tree.Output)

            os.system('mv uclust_%s %s_uclust' % (update_tree.Output, update_tree.Output))
            os.system('mv usort_%s %s_uclust' % (update_tree.Output, update_tree.Output))

            os.system('mv %s_uclust %s' % (update_tree.Output, project_folder))

    # STAGE 7: Delete files as determined by the user
    delete_files(args)

    sys.stdout.write("MLTreeMap has finished successfully.\n")
    sys.stdout.flush()


if __name__ == "__main__":
    main(sys.argv[1:])
