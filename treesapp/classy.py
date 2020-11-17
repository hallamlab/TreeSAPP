
__author__ = 'Connor Morgan-Lang'

import sys
import os
import re
import logging
import time
from datetime import datetime as dt
from shutil import rmtree, copy
from glob import glob
from collections import namedtuple
from numpy import var

from treesapp.phylo_seq import convert_entrez_to_tree_leaf_references, PQuery
from treesapp.refpkg import ReferencePackage
from treesapp.fasta import fastx_split, get_header_format, FASTA, load_fasta_header_regexes, sequence_info_groups
from treesapp.utilities import median, write_dict_to_table, validate_new_dir, fetch_executable_path
from treesapp.lca_calculations import determine_offset, optimal_taxonomic_assignment
from treesapp import entrez_utils
from treesapp.wrapper import CommandLineFarmer, estimate_ml_model


class ModuleFunction:
    def __init__(self, name: str, order: int, func=None):
        """
        Create a new instance of the ModuleFunction class

        :param name: Name of the function's module (e.g. clean, lineages, assign, build)
        :param order: The order in which the module is supposed to be ran
        :param func: An optional function to call when this module is reached
        """
        self.order = order
        self.name = name
        self.function = func
        self.run = True

    def get_info(self):
        info_string = "Information for '" + self.name + "':\n"
        info_string += "\tOrder: " + str(self.order) + "\n"
        info_string += "\tRun: " + str(self.run) + "\n"
        return info_string


# class NodeRetrieverWorker(Process):
#     """
#     Doug Hellman's Consumer class for handling processes via queues
#     """
#
#     def __init__(self, task_queue, result_queue):
#         Process.__init__(self)
#         self.task_queue = task_queue
#         self.result_queue = result_queue
#
#     def run(self):
#         while True:
#             next_task = self.task_queue.get()
#             if next_task is None:
#                 # Poison pill means shutdown
#                 self.task_queue.task_done()
#                 break
#             result = _tree_parser._build_subtrees_newick(next_task)
#             subtrees = subtrees_to_dictionary(result, create_tree_info_hash())
#             self.task_queue.task_done()
#             self.result_queue.put(subtrees)
#         return


def get_header_info(header_registry: dict, code_name=''):
    """

    :param header_registry: A dictionary of Header instances, indexed by numerical treesapp_id
    :param code_name: [OPTIONAL] The code_name of the reference package (marker gene/domain/family/protein)
    :return: Dictionary where keys are numerical treesapp_ids and values are EntrezRecord instances
    """
    logging.info("Extracting information from headers... ")
    ref_records = dict()
    header_regexes = load_fasta_header_regexes(code_name)
    # TODO: Fix parsing of combined EggNOG and custom headers such that the taxid is parsed from the "accession"
    for treesapp_id in sorted(header_registry.keys(), key=int):  # type: str
        original_header = header_registry[treesapp_id].original
        header_format_re, header_db, header_molecule = get_header_format(original_header, header_regexes)
        sequence_info = header_format_re.match(original_header)
        seq_info_tuple = sequence_info_groups(sequence_info, header_db, original_header, header_regexes)

        # Load the parsed sequences info into the EntrezRecord objects
        ref_seq = entrez_utils.EntrezRecord(seq_info_tuple.accession, seq_info_tuple.version)
        ref_seq.organism = seq_info_tuple.organism
        ref_seq.lineage = seq_info_tuple.lineage
        ref_seq.ncbi_tax = seq_info_tuple.taxid
        ref_seq.description = seq_info_tuple.description
        ref_seq.locus = seq_info_tuple.locus
        ref_seq.short_id = treesapp_id + '_' + code_name
        ref_seq.tracking_stamp()
        ref_records[treesapp_id] = ref_seq
    logging.info("done.\n")

    return ref_records


def dedup_records(ref_seqs: FASTA, ref_seq_records: dict) -> dict:
    """
    Entrez records with identical versioned accessions are grouped together and those grouped sequence records
    are compared further. The Entrez records in ref_seq_records are deduplicated based on:
 1. Their full-length original headers
 2. Their EntrezRecord.bitflag values where the records with the most information attributed to them (e.g. NCBI tax_id,
    lineage) at this stage are retained.

    :param ref_seqs: A FASTA instance with populated fasta_dict and header_registry attributes
    :param ref_seq_records: A dictionary mapping TreeSAPP numerical identifiers to their EntrezRecord instances
    :return: ref_seq_records dictionary with duplicate records removed
    """
    duplicate_treesapp_ids = set()
    all_accessions = dict()

    if ref_seqs.index_form != "num":
        ref_seqs.change_dict_keys("num")

    # Create a dictionary mapping versioned accessions to EntrezRecords for identifying duplicate accessions
    for treesapp_id, record in ref_seq_records.items():  # type: (str, entrez_utils.EntrezRecord)
        if record.versioned not in all_accessions:
            all_accessions[record.versioned] = []
        all_accessions[record.versioned].append(treesapp_id)

    if len(all_accessions) != ref_seqs.n_seqs():
        logging.debug("{}/{} unique versioned accessions were loaded for deduplication.\n"
                      "".format(len(all_accessions), ref_seqs.n_seqs()))

    # If the sequences are identical across the records -> take record with greater bitflag
    # If the sequences are different -> keep
    for accession, dup_list in all_accessions.items():  # type: (str, list)
        # TODO: Deduplicate based on whether sequences or headers are substrings of others, using a Trie
        if len(dup_list) > 1:
            ##
            # Remove records with redundant accessions by removing those with the lower bitflag
            ##
            max_bitflag = max([ref_seq_records[treesapp_id].bitflag for treesapp_id in dup_list])
            x = 0
            while x < len(dup_list):
                treesapp_id = dup_list[x]
                ref_seq = ref_seq_records[treesapp_id]  # type: entrez_utils.EntrezRecord
                if ref_seq.bitflag < max_bitflag:
                    duplicate_treesapp_ids.add(treesapp_id)
                    dup_list.pop(x)
                else:
                    x += 1
            # Check for records with identical headers, keeping only one
            og_headers = dict()
            for ts_id in dup_list:  # type: str
                try:
                    og_headers[ref_seqs.header_registry[ts_id].original].append(ts_id)
                except KeyError:
                    og_headers[ref_seqs.header_registry[ts_id].original] = [ts_id]
            for seq_name, ts_ids in og_headers.items():
                if len(ts_ids) > 1:
                    for treesapp_id in ts_ids[1:]:
                        duplicate_treesapp_ids.add(treesapp_id)
        else:
            pass

    # Actually remove the duplicates from the ref_seq_records dictionary
    if duplicate_treesapp_ids:
        deduped_accessions = []
        for ts_id in duplicate_treesapp_ids:
            ref_seq_records.pop(ts_id)
            deduped_accessions.append(ref_seqs.header_registry[ts_id].original)

        logging.warning("The following sequences were removed during deduplication of Entrez records:\n\t" +
                        "\n\t".join(deduped_accessions) + "\n")

    return ref_seq_records


class Cluster:
    def __init__(self, rep_name):
        self.representative = rep_name
        self.members = list()
        self.lca = ''

    def get_info(self):
        info_string = "Representative: " + str(self.representative) + "\n" + \
                      "LCA: " + self.lca + "\n" + \
                      "Members:\n\t" + "\n\t".join([', '.join(member) for member in self.members]) + "\n"
        return info_string


class BlastAln:
    def __init__(self):
        self.subject = ""
        self.query = ""
        self.start = 0
        self.end = 0
        self.alnlen = 0
        self.gapopen = 0
        self.mismatch = 0
        self.qstart = 0
        self.qend = 0
        self.tstart = 0
        self.tend = 0
        self.evalue = 0.0
        self.pident = 0.0
        self.bits = 0.0

    def load_blast_tab(self, line: str, sep="\t") -> None:
        try:
            fields = line.strip().split(sep)
        except ValueError:
            logging.error("Unable to parse line in from alignment table:\n{}\n".format(line))
            sys.exit(7)

        try:
            self.subject, self.query = fields[0], fields[1]
            self.pident, self.alnlen, self.mismatch, self.gapopen = fields[2:6]
            self.qstart, self.qend, self.tstart, self.tend = fields[6:10]
            self.evalue, self.bits = fields[10], fields[11]
        except ValueError:
            logging.error("Incorrect format for line in alignment file. Twelve were expected, found {}.\n{}\n"
                          "".format(len(fields), line))
            sys.exit(7)

        return


class MyFormatter(logging.Formatter):

    error_fmt = "%(levelname)s - %(module)s, line %(lineno)d:\n%(message)s"
    warning_fmt = "%(levelname)s:\n%(message)s"
    debug_fmt = "%(asctime)s\n%(message)s"
    info_fmt = "%(message)s"

    def __init__(self):
        super().__init__(fmt="%(levelname)s: %(message)s",
                         datefmt="%d/%m %H:%M:%S")

    def format(self, record):

        # Save the original format configured by the user
        # when the logger formatter was instantiated
        format_orig = self._style._fmt

        # Replace the original format with one customized by logging level
        if record.levelno == logging.DEBUG:
            self._style._fmt = MyFormatter.debug_fmt

        elif record.levelno == logging.INFO:
            self._style._fmt = MyFormatter.info_fmt

        elif record.levelno == logging.ERROR:
            self._style._fmt = MyFormatter.error_fmt

        elif record.levelno == logging.WARNING:
            self._style._fmt = MyFormatter.warning_fmt

        # Call the original formatter class to do the grunt work
        result = logging.Formatter.format(self, record)

        # Restore the original format configured by the user
        self._style._fmt = format_orig

        return result


def prep_logging(log_file=None, verbosity=False, stream=sys.stderr) -> None:
    """
    Allows for multiple file handlers to be added to the root logger, but only a single stream handler.
    The new file handlers must be removed outside of this function explicitly

    :param log_file: Path to a file to write the TreeSAPP log
    :param verbosity: Whether debug-level information should be written (True) or not (False)
    :param stream: Which stream, sys.stdout or sys.stderr, should the console logger write to?
    :return: None
    """
    if verbosity:
        logging_level = logging.DEBUG
    else:
        logging_level = logging.INFO

    # Detect whether handlers are already present and return if true
    logger = logging.getLogger()
    if len(logger.handlers):
        return

    formatter = MyFormatter()
    # Set the console handler normally writing to stdout/stderr
    console_handler = logging.StreamHandler(stream=stream)
    console_handler.setLevel(logging_level)
    console_handler.terminator = ''
    console_handler.setFormatter(formatter)

    if log_file:
        if not os.path.isabs(log_file):
            log_file = os.path.join(os.getcwd(), os.path.dirname(log_file), os.path.basename(log_file))
        output_dir = os.path.dirname(log_file)
        try:
            if output_dir and not os.path.isdir(output_dir):
                os.mkdir(output_dir)
        except (IOError, OSError):
            sys.stderr.write("ERROR: Unable to make directory '" + output_dir + "'.\n")
            sys.exit(3)
        logging.basicConfig(level=logging.DEBUG,
                            filename=log_file,
                            filemode='w',
                            datefmt="%d/%m %H:%M:%S",
                            format="%(asctime)s %(levelname)s:\n%(message)s")
        logging.getLogger('').addHandler(console_handler)
        logging.getLogger('').propagate = False
    else:
        logging.basicConfig(level=logging_level,
                            datefmt="%d/%m %H:%M:%S",
                            format="%(asctime)s %(levelname)s:\n%(message)s")
    return


class TreeSAPP:
    """
    Abstract class for each of the different analysis types - create, evaluate, assign, update and train
    """
    def __init__(self, cmd):
        # Static values
        self.command = cmd
        # self.refpkg_code_re = re.compile(r'[A-Z][0-9]{4,5}')
        self.treesapp_dir = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + os.sep
        self.refpkg_dir = self.treesapp_dir + 'data' + os.sep
        self.itol_dir = self.treesapp_dir + 'data' + os.sep + "iTOL_data" + os.sep
        # Necessary for Evaluator, Creator and PhyTrainer:
        self.seq_lineage_map = dict()  # Dictionary holding the accession-lineage mapping information
        self.acc_to_lin = ""  # Path to an accession-lineage mapping file
        self.ref_pkg = ReferencePackage()

        # Values derived from the command-line arguments
        self.input_sequences = ""
        self.query_sequences = ""  # This could be any of the input_sequences aa_orfs or nuc_orfs files
        self.sample_prefix = ""
        self.formatted_input = ""
        self.molecule_type = ""
        self.output_dir = ""
        self.final_output_dir = ""
        self.var_output_dir = ""
        self.stage_output_dir = ""
        self.executables = dict()
        # Values that need to be entered later, in the command-specific class
        self.stages = dict()  # Used to track what progress stages need to be completed
        self.stage_file = ""  # The file to write progress updates to
        self.current_stage = None

    def get_info(self):
        info_string = "Executables:\n\t" + "\n\t".join([k + ": " + v for k, v in self.executables.items()]) + "\n"

        for module_step in self.stages:
            info_string += self.stages[module_step].get_info()
        info_string += "\n\t".join(["TreeSAPP directory = " + self.treesapp_dir,
                                    "FASTA input sequences = " + self.input_sequences,
                                    "Formatted input = " + self.formatted_input,
                                    "Output directory = " + self.output_dir,
                                    "Input molecule type = " + self.molecule_type])
        return info_string

    def furnish_with_arguments(self, args) -> None:
        """
        Carries over the basic TreeSAPP arguments to the respective TreeSAPP-subclass.
        All auxiliary arguments are pushed to the TreeSAPP classes in check_module_arguments

        :param args: arguments from argparse.ParseArgs() with output, input and molecule attributes
        :return: None
        """
        if self.command != "info":
            self.output_dir = validate_new_dir(args.output)
            self.final_output_dir = self.output_dir + "final_outputs" + os.sep
            self.var_output_dir = self.output_dir + "intermediates" + os.sep
            if set(vars(args)).issuperset({"molecule", "input"}):
                self.input_sequences = args.input
                self.molecule_type = args.molecule
                file_name, suffix1 = os.path.splitext(os.path.basename(self.input_sequences))
                if suffix1 == ".gz":
                    file_name, suffix2 = os.path.splitext(file_name)
                self.sample_prefix = file_name
                self.formatted_input = self.var_output_dir + self.sample_prefix + "_formatted.fasta"

        if self.command != "colour" and "pkg_path" in vars(args):
            if len(args.pkg_path) > 1:
                logging.warning("Multiple reference packages cannot be used by treesapp {}.\n"
                                "Only the first one ({}) will be used.\n".format(self.command, args.pkg_path))
            args.pkg_path = args.pkg_path.pop(0)

        self.executables = self.find_executables(args)
        return

    def validate_refpkg_dir(self, refpkg_dir: str):
        if refpkg_dir:
            if not os.path.isdir(refpkg_dir):
                logging.error("Directory containing reference packages ({}) does not exist.\n".format(refpkg_dir))
                sys.exit(5)
            self.refpkg_dir = refpkg_dir
        return

    def check_previous_output(self, overwrite=True) -> None:
        """
        Prompts the user to determine how to deal with a pre-existing output directory.
        By the end of this function, all directories should exist and be in the correct state for a new analysis run

        :param overwrite: Boolean flag controlling whether output directories are removed or not
        :return None
        """

        main_output_dirs = [self.var_output_dir, self.final_output_dir]

        # Identify all the various reasons someone may not want to have their previous results overwritten
        if not os.path.isdir(self.output_dir):
            os.mkdir(self.output_dir)
        elif not overwrite and glob(self.final_output_dir + "*"):
            # reluctant_remove_replace(self.output_dir)
            pass
        elif not overwrite and not glob(self.final_output_dir + "*"):
            # Last run didn't complete so use the intermediates if possible
            pass
        elif overwrite:
            if os.path.isdir(self.output_dir):
                rmtree(self.output_dir)
            os.mkdir(self.output_dir)
        # Create the output directories
        for output_dir in main_output_dirs:
            if not os.path.isdir(output_dir):
                os.mkdir(output_dir)
        return

    def log_progress(self):
        """
        Write the current stage's information (e.g. number of output files) to a JSON file

        :return: None
        """
        # TODO: Finish this to enable continuing part-way through an analysis
        if os.path.isdir(self.var_output_dir):
            jplace_files = glob(self.var_output_dir + os.sep + "*jplace")
            if len(jplace_files) >= 1:
                self.stages["classify"] = True
            else:
                logging.warning("Reclassify impossible as " + self.output_dir + " is missing input files.\n")
        return

    def stage_lookup(self, name: str, tolerant=False) -> ModuleFunction:
        """
        Used for looking up a stage in self.stages by its stage.name

        :param name: Name of a stage
        :param tolerant: Boolean controlling whether the function exits if look-up failed (defualt) or returns None
        :return: ModuleFunction instance that matches the name, or None if failed and tolerant, exit otherwise
        """
        found = False
        stage_order = 0
        for module_step in sorted(self.stages, key=int):  # type: int
            stage = self.stages[module_step]  # type: ModuleFunction
            stage_order = stage.order
            if stage.name == name:
                return stage

        if not found and not tolerant:
            logging.error("Unable to find '{}' in {} stages.\n".format(name, self.command))
            sys.exit(3)
        else:
            logging.warning("Unable to find '{}' stage. Returning a new one instead.\n".format(name))
            return ModuleFunction(name=name, order=stage_order+1)

    def first_stage(self):
        for x in sorted(self.stages, key=int):  # type: int
            stage = self.stages[x]  # type: ModuleFunction
            if stage.run:
                return stage.order
        logging.error("No stages are set to run!\n")
        sys.exit(3)

    # def read_progress_log(self):
    #     """
    #     Read the object's 'stage_file' and determine the completed stage
    #
    #     :return: An integer corresponding to the last completed stage's rank in 'stage_order'
    #     """
    #     completed_int = self.first_stage()
    #     try:
    #         progress_handler = open(self.stage_file)
    #     except IOError:
    #         logging.debug("Unable to open stage file '" + self.stage_file +
    #                       "' for reading. Defaulting to stage " + str(completed_int) + ".\n")
    #         return completed_int
    #     # TODO: Finish this to enable continuing part-way through an analysis
    #     progress_handler.close()
    #     return completed_int

    def find_stage_dirs(self) -> int:
        """
        Selects the earliest checkpoint that the workflow can be started from.
        Stages (ModuleFunction instances) are skipped by setting their 'run' attribute to False.
        This is the stage preceding the first stage whose directory does not exist.

        :return: None
        """
        for i, module in sorted(self.stages.items()):  # type: (int, ModuleFunction)
            if not module.run:
                continue
            else:
                return module.order
        return 0

    def set_stage_dir(self) -> None:
        self.stage_output_dir = os.path.join(self.var_output_dir, self.current_stage.name) + os.sep
        if not os.path.isdir(self.stage_output_dir):
            os.mkdir(self.stage_output_dir)
        return

    def increment_stage_dir(self) -> None:
        """
        Updates self.stage_output_dir with the directory path of the next ModuleFunction.

        :return: None
        """
        # Find the next stage
        self.current_stage.run = False
        while not self.current_stage.run and self.current_stage.order < max(self.stages.keys()):
            self.current_stage = self.stages[self.current_stage.order+1]

        # Update the output directory for this stage
        self.set_stage_dir()
        return

    def stage_status(self, name):
        return self.stage_lookup(name).run

    def change_stage_status(self, name: str, new_status: bool):
        stage = self.stage_lookup(name)
        stage.run = new_status
        return

    def edit_stages(self, start, end=None):
        for x in sorted(self.stages, key=int):  # type: int
            stage = self.stages[x]  # type: ModuleFunction
            if end is not None:
                if x < start or x > end:
                    stage.run = False
            elif x < start:
                stage.run = False
            else:
                # These stages need to be ran
                pass

        return

    def validate_continue(self, args) -> None:
        """
        Bases the stage(s) to run on args.stage which is broadly set to either 'continue' or any other valid stage

        This function ensures all the required inputs are present for beginning at the desired first stage,
        otherwise, the pipeline begins at the first possible stage to continue and ends once the desired stage is done.

        If a directory exists with the name of a stage, that step is skipped.
        :return: None
        """
        # TODO: Summarise the steps to be taken and write to log
        for i, module in sorted(self.stages.items()):  # type: (int, ModuleFunction)
            if os.path.isdir(os.path.join(self.var_output_dir, module.name)):
                module.run = False

        if args.overwrite:
            last_valid_stage = self.first_stage()
        else:
            last_valid_stage = self.find_stage_dirs()

        self.current_stage = self.stages[last_valid_stage]
        self.set_stage_dir()
        if args.stage == "continue":
            logging.debug("Continuing with stage '{}'\n".format(self.current_stage.name))
            self.edit_stages(last_valid_stage)
            return

        # Update the stage status
        desired_stage = self.stage_lookup(args.stage).order
        if desired_stage > last_valid_stage:
            logging.warning("Unable to run '{}' as it is ahead of the last completed stage.\n"
                            "Continuing with stage '{}'\n".format(args.stage, self.current_stage.name))
            self.edit_stages(last_valid_stage, desired_stage)
        elif desired_stage < last_valid_stage:
            logging.info("Wow - its your lucky day:\n"
                         "All stages up to and including '{}' have already been completed!\n".format(args.stage))
            self.edit_stages(desired_stage, desired_stage)
        else:
            # Proceed with running the desired stage
            logging.debug("Proceeding with '{}'\n".format(args.stage))
            self.edit_stages(desired_stage, desired_stage)

        return

    def find_executables(self, args):
        """
        Finds the executables in a user's path to alleviate the requirement of a sub_binaries directory

        :param args: The parsed command-line arguments
        :return: exec_paths beings the absolute path to each executable
        """
        exec_paths = dict()
        dependencies = ["prodigal", "hmmbuild", "hmmalign", "hmmsearch", "epa-ng", "raxml-ng", "BMGE.jar"]

        # Extra executables necessary for certain modes of TreeSAPP
        if self.command == "abundance":
            dependencies += ["bwa"]

        if self.command in ["create", "update", "train", "evaluate"]:
            dependencies += ["mmseqs", "mafft"]
            if hasattr(args, "fast") and args.fast:
                dependencies.append("FastTree")

            if hasattr(args, "od_seq") and args.od_seq:
                dependencies.append("OD-seq")

        if self.molecule_type == "rrna":
            dependencies += ["cmalign", "cmsearch", "cmbuild"]

        for dep in dependencies:
            exec_paths[dep] = fetch_executable_path(dep, self.treesapp_dir)

        return exec_paths

    def fetch_entrez_lineages(self, ref_seqs: FASTA, molecule: str, acc_to_taxid=None, seqs_to_lineage=None) -> dict:
        """
        The root function orchestrating download of taxonomic lineage information using BioPython's Entrez API.
        In addition to this, the function supports a couple different tables containing taxonomic lineage information
        such as organism names, NCBI 'taxid's, and complete lineages associated with their respective sequence names.

        :param ref_seqs: A populated FASTA instance. From this instance, the header_registry is parsed and the format
        (i.e. database-specific header format) is assumed using a host of regular expressions to pull the relevant
        information that can be used for downloading the whole taxonomic lineage using the Entrez API.
        :param molecule: The molecule type of the sequences in ref_seqs. Either 'prot', 'dna', 'rrna', or 'ambig'.
        :param acc_to_taxid: Path to a table mapping NCBI accessions to NCBI taxids
        :param seqs_to_lineage: Path to a table mapping sequence names to their respective taxonomic lineage
        :return: A dictionary mapping unique numerical TreeSAPP identifiers to EntrezRecord instances
        """
        # Get the lineage information for the training/query sequences
        entrez_record_dict = get_header_info(ref_seqs.header_registry, self.ref_pkg.prefix)
        entrez_record_dict = dedup_records(ref_seqs, entrez_record_dict)
        ref_seqs.change_dict_keys("formatted")
        entrez_utils.load_ref_seqs(ref_seqs.fasta_dict, ref_seqs.header_registry, entrez_record_dict)
        logging.debug("\tNumber of input sequences = {}\n".format(len(entrez_record_dict)))

        # Seed the seq_lineage_map with any lineages that were parsed from the FASTA file
        for ts_id in entrez_record_dict:
            e_record = entrez_record_dict[ts_id]  # type: entrez_utils.EntrezRecord
            if e_record.lineage:
                self.seq_lineage_map[e_record.accession] = e_record.lineage
                e_record.lineage = ""

        if seqs_to_lineage:
            lineage_map, refs_mapped = entrez_utils.map_orf_lineages(seqs_to_lineage, ref_seqs.header_registry)
            # Add lineage information to entrez records for each reference sequence
            entrez_utils.fill_ref_seq_lineages(entrez_record_dict, lineage_map,
                                               complete=(len(refs_mapped) == ref_seqs.n_seqs()))
            ref_leaf_nodes = convert_entrez_to_tree_leaf_references(entrez_record_dict)
            self.ref_pkg.taxa_trie.feed_leaf_nodes(ref_leaf_nodes)
            entrez_utils.sync_record_and_hierarchy_lineages(ref_leaf_nodes, entrez_record_dict)
            self.ref_pkg.taxa_trie.validate_rank_prefixes()
            self.ref_pkg.taxa_trie.build_multifurcating_trie()
        if self.stage_status("lineages"):
            entrez_query_list, num_lineages_provided = entrez_utils.build_entrez_queries(entrez_record_dict)
            logging.debug("\tNumber of queries =\t" + str(len(entrez_query_list)) + "\n")
            if len(entrez_query_list) == 0:
                return entrez_record_dict
            entrez_utils.map_accessions_to_lineages(entrez_query_list, self.ref_pkg.taxa_trie, molecule, acc_to_taxid)
            # Repair entrez_record instances either lacking lineages or whose lineages do not contain rank-prefixes
            entrez_utils.repair_lineages(entrez_record_dict, self.ref_pkg.taxa_trie)
            self.seq_lineage_map = entrez_utils.entrez_records_to_accession_lineage_map(entrez_query_list)
            # Map proper accession to lineage from the tuple keys (accession, accession.version)
            #  in accession_lineage_map returned by entrez_utils.get_multiple_lineages.
            entrez_utils.verify_lineage_information(self.seq_lineage_map, entrez_record_dict,
                                                    self.ref_pkg.taxa_trie, num_lineages_provided)

            self.seq_lineage_map = entrez_utils.accession_lineage_map_from_entrez_records(entrez_record_dict)
            # Ensure the accession IDs are stripped of '>'s
            for accession in sorted(self.seq_lineage_map):
                if accession[0] == '>':
                    self.seq_lineage_map[accession[1:]] = self.seq_lineage_map.pop(accession)
            write_dict_to_table(self.seq_lineage_map, self.acc_to_lin)
            self.increment_stage_dir()
        else:
            logging.info("Reading cached lineages in '{}'... ".format(self.acc_to_lin))
            self.seq_lineage_map.update(entrez_utils.read_accession_taxa_map(self.acc_to_lin))
            logging.info("done.\n")

        ref_seqs.change_dict_keys()
        return entrez_record_dict


class Updater(TreeSAPP):
    def __init__(self):
        super(Updater, self).__init__("update")
        self.seq_names_to_taxa = ""  # Optional user-provided file mapping query sequence contigs to lineages
        self.lineage_map_file = ""  # File that is passed to create() containing lineage info for all sequences
        self.treesapp_output = ""  # Path to the TreeSAPP output directory - modified by args
        self.assignment_table = ""  # Path to the marker_contig_map.tsv file written by treesapp assign
        self.combined_fasta = ""  # Holds the newly identified candidate reference sequences and the original ref seqs
        self.old_ref_fasta = ""  # Contains only the original reference sequences
        self.cluster_input = ""  # Used only if resolve is True
        self.clusters_prefix = ""  # Used only if resolve is True
        self.updated_refpkg_path = ""
        # self.rank_depth_map = None
        self.prop_sim = 1.0
        self.min_length = 0  # The minimum sequence length for a classified sequence to be included in the refpkg
        self.updated_refpkg = ReferencePackage()

        # Stage names only holds the required stages; auxiliary stages (e.g. RPKM, update) are added elsewhere
        self.stages = {0: ModuleFunction("lineages", 0),
                       1: ModuleFunction("rebuild", 1)}

    def decide_stage(self, args):
        """
        Bases the stage(s) to run on args.stage which is broadly set to either 'continue' or any other valid stage

        This function ensures all the required inputs are present for beginning at the desired first stage,
        otherwise, the pipeline begins at the first possible stage to continue and ends once the desired stage is done.

        If a
        :return: None
        """
        self.acc_to_lin = self.var_output_dir + os.sep + "accession_id_lineage_map.tsv"
        self.validate_continue(args)
        return

    def get_info(self):
        info_string = "Updater instance summary:\n"
        info_string += super(Updater, self).get_info() + "\n\t"
        info_string += "\n\t".join(["Target reference packages = " + str(self.ref_pkg.prefix),
                                    "Taxonomy map: " + self.ref_pkg.lineage_ids,
                                    "Reference tree: " + self.ref_pkg.tree,
                                    "Reference FASTA: " + self.ref_pkg.msa,
                                    "Lineage map: " + str(self.seq_names_to_taxa)]) + "\n"

        return info_string

    def update_refpkg_fields(self) -> None:
        """
        Using the original ReferencePackage as a template modify the following updated ReferencePackage attributes:
1. original creation date
2. update date
3. code name
4. description
        :return: None
        """
        # Change the creation and update dates, code name and description
        self.updated_refpkg.date = self.ref_pkg.date
        self.updated_refpkg.update = dt.now().strftime("%Y-%m-%d")
        self.updated_refpkg.refpkg_code = self.ref_pkg.refpkg_code
        self.updated_refpkg.description = self.ref_pkg.description
        self.updated_refpkg.pickle_package()

        logging.info("Summary of the updated reference package:\n" + self.updated_refpkg.get_info() + "\n")

        logging.debug("\tNew sequences  = " + str(self.updated_refpkg.num_seqs - self.ref_pkg.num_seqs) + "\n" +
                      "\tOld HMM length = " + str(self.ref_pkg.hmm_length()) + "\n" +
                      "\tNew HMM length = " + str(self.updated_refpkg.hmm_length()) + "\n")

        return


class Creator(TreeSAPP):
    def __init__(self):
        super(Creator, self).__init__("create")
        self.candidates = dict()  # Dictionary tracking all candidate ReferenceSequences
        self.phy_dir = ""  # Directory for intermediate or unnecessary files created during phylogeny inference
        self.hmm_purified_seqs = ""  # If an HMM profile of the gene is provided its a path to FASTA with homologs
        self.filtered_fasta = ""
        self.hmm_profile = ""  # HMM profile used for screening the input sequences
        self.clusters_prefix = ""  # FASTA file prefix for cluster centroids
        self.clusters_table = ""  # Output file defining the clusters, members, representatives, and similarity
        self.cluster_input = ""  # Name of the file to be used for clustering
        self.unaln_ref_fasta = ""  # FASTA file of unaligned reference sequences
        self.phylip_file = ""  # Used for building the phylogenetic tree with RAxML
        self.min_tax_rank = "Kingdom"  # Minimum taxonomic rank
        self.metadata_file = ""
        self.training_dir = ""

        # Stage names only holds the required stages; auxiliary stages (e.g. RPKM, update) are added elsewhere
        self.stages = {0: ModuleFunction("search", 0),
                       1: ModuleFunction("lineages", 1),
                       2: ModuleFunction("clean", 2),
                       3: ModuleFunction("cluster", 3),
                       4: ModuleFunction("build", 4),
                       5: ModuleFunction("evaluate", 5),
                       6: ModuleFunction("support", 6),
                       7: ModuleFunction("train", 7),
                       8: ModuleFunction("update", 8)}

    def decide_stage(self, args):
        """
        Bases the stage(s) to run on args.stage which is broadly set to either 'continue' or any other valid stage

        This function ensures all the required inputs are present for beginning at the desired first stage,
        otherwise, the pipeline begins at the first possible stage to continue and ends once the desired stage is done.

        If a
        :return: None
        """
        if not args.profile:
            self.change_stage_status("search", False)
        if args.pc:
            self.edit_stages(self.stage_lookup("update").order)
        # TODO: Allow users to provide sequence-lineage maps for a subset of the query sequences
        if args.acc_to_lin:
            self.acc_to_lin = args.acc_to_lin
            if os.path.isfile(self.acc_to_lin):
                self.change_stage_status("lineages", False)
            else:
                logging.error("Unable to find accession-lineage mapping file '{}'\n".format(self.acc_to_lin))
                sys.exit(3)
        else:
            self.acc_to_lin = self.var_output_dir + os.sep + "accession_id_lineage_map.tsv"

        if args.bootstraps == 0:
            self.change_stage_status("support", False)
        if not args.fast:
            self.change_stage_status("evaluate", False)
        self.validate_continue(args)
        return

    def get_info(self):
        info_string = "Creator instance summary:\n"
        info_string += super(Creator, self).get_info() + "\n\t"
        return info_string

    def remove_intermediates(self, detonate=False) -> None:
        if os.path.exists(self.unaln_ref_fasta):
            os.remove(self.unaln_ref_fasta)
        if os.path.exists(self.phylip_file + ".reduced"):
            os.remove(self.phylip_file + ".reduced")
        if os.path.exists(self.final_output_dir + "fasta_reader_log.txt"):
            os.remove(self.final_output_dir + "fasta_reader_log.txt")
        if os.path.exists(self.phylip_file):
            copy(self.phylip_file, self.phy_dir)
            os.remove(self.phylip_file)

        if detonate and os.path.isdir(self.var_output_dir):
            rmtree(self.var_output_dir)
        return

    def determine_model(self, refpkg: ReferencePackage, estimate=False) -> None:
        if refpkg.tree_tool == "FastTree":
            if refpkg.molecule == "prot":
                evo_model = "LG+G4"
            elif refpkg.molecule == "rrna" or refpkg.molecule == "dna":
                evo_model = "GTR+G"
            else:
                logging.error("Unrecognized reference package molecule type: '{}'.\n".format(refpkg.molecule))
                sys.exit(3)
            if refpkg.sub_model:
                logging.warning("Model provided '{}' will be ignored when FastTree is used to infer phylogeny.\n"
                                "".format(refpkg.sub_model))
        elif refpkg.tree_tool == "RAxML-NG":
            if estimate:
                evo_model = estimate_ml_model(modeltest_exe=self.executables["ModelTest-NG"],
                                              msa=refpkg.f__msa, output_prefix=self.phy_dir, molecule=refpkg.molecule)
            elif not refpkg.sub_model:
                if refpkg.molecule == "prot":
                    evo_model = "LG+G4"
                else:
                    evo_model = "GTR+G"
            else:
                logging.debug("Using specified RAxML-NG-compatible model: '{}'.\n".format(refpkg.sub_model))
                evo_model = refpkg.sub_model
        else:
            logging.error("Unexpected phylogenetic inference tool: '{}'.\n".format(refpkg.tree_tool))
            sys.exit(3)

        refpkg.sub_model = evo_model
        return

    def print_terminal_commands(self):
        logging.info("\nTo integrate this package for use in TreeSAPP the following steps must be performed:\n"
                     "1. Replace the current refpkg_code 'Z1111' with:\n"
                     "`treesapp package edit refpkg_code $code --overwrite --refpkg_path {0}`"
                     " where $code is a unique identifier.\n"
                     "2. Copy {0} to a directory containing other reference packages you want to analyse. "
                     "This may be in {1}/data/ or elsewhere\n"
                     "".format(self.ref_pkg.f__json, self.treesapp_dir))
        return


class Purity(TreeSAPP):
    def __init__(self):
        super(Purity, self).__init__("purity")
        self.assign_dir = ""
        self.classifications = ""
        self.summarize_dir = ""
        self.metadata_file = ""
        self.assignments = None
        self.stages = {0: ModuleFunction("assign", 0),
                       1: ModuleFunction("summarize", 1)}

    def decide_stage(self, args):
        """
        Bases the stage(s) to run on args.stage which is broadly set to either 'continue' or any other valid stage

        This function ensures all the required inputs are present for beginning at the desired first stage,
        otherwise, the pipeline begins at the first possible stage to continue and ends once the desired stage is done.

        If a
        :return: None
        """
        self.validate_continue(args)
        return

    def summarize_groups_assigned(self, ortholog_map: dict, metadata=None):
        unique_orthologs = dict()
        tree_leaves = self.ref_pkg.generate_tree_leaf_references_from_refpkg()
        for refpkg, info in self.assignments.items():
            for lineage in info:
                for seq_name in info[lineage]:
                    bits_pieces = seq_name.split('_')
                    if bits_pieces[0] not in unique_orthologs:
                        unique_orthologs[bits_pieces[0]] = []
                    unique_orthologs[bits_pieces[0]].append('_'.join(bits_pieces[1:]))

        # Summarize the counts in the log file
        summary_str = "Ortholog\tHits\tLeaves\tTree-coverage\tDescription\n" + '-'*80 + "\n"
        for og_name in sorted(ortholog_map, key=lambda x: len(ortholog_map[x])):
            n_leaves = len(ortholog_map[og_name])
            perc_coverage = round(float((n_leaves*100)/len(tree_leaves)), 1)
            try:
                desc = metadata[og_name].de
            except (TypeError, KeyError):
                desc = "NA"
            summary_str += "\t".join([og_name,
                                      str(len(unique_orthologs[og_name])), str(n_leaves), str(perc_coverage),
                                      desc]) + "\n"
        logging.info(summary_str + "\n")

        return

    def load_metadata(self) -> dict:
        metadat_dict = dict()
        xtra_dat = namedtuple("xtra_dat", ["db_id", "ac", "de"])
        if not os.path.isfile(self.metadata_file):
            logging.error("Extra information file '" + self.metadata_file + "' doesn't exist!\n")
            sys.exit(3)
        try:
            metadata_handler = open(self.metadata_file, 'r')
        except IOError:
            logging.error("Unable to open extra information file '" + self.metadata_file + "' for reading!\n")
            sys.exit(3)
        for line in metadata_handler:
            try:
                db_id, accession, desc = line.strip().split("\t")
            except ValueError:
                logging.error("Bad format for '" + self.metadata_file + "'. Three tab-separated fields expected.\n" +
                              "Example line:\n" +
                              str(line) + "\n")
                sys.exit(7)
            metadat_dict[accession] = xtra_dat(db_id, accession, desc)
        metadata_handler.close()
        return metadat_dict

    def assign_leaves_to_orthologs(self, p_queries: list, internal_node_map: dict) -> dict:
        ortholog_map = dict()
        leaf_map = dict()
        tree_leaves = self.ref_pkg.generate_tree_leaf_references_from_refpkg()
        for leaf in tree_leaves:
            leaf_map[leaf.number + "_" + self.ref_pkg.prefix] = leaf.description
        for p_query in p_queries:  # type: PQuery
            p_query.ref_name = self.ref_pkg.prefix
            if type(p_query.place_name) is list and len(p_query.place_name):
                p_query.place_name = p_query.place_name[0]
            seq_info = re.match(r"(.*)\|" + re.escape(p_query.ref_name) + r"\|(\\d+)_(\\d+)$", p_query.place_name)
            if seq_info:
                p_query.place_name = seq_info.group(1)
            leaves = internal_node_map[p_query.consensus_placement.edge_num]
            ortholog_name = p_query.place_name.split('_')[0]
            if ortholog_name not in ortholog_map:
                ortholog_map[ortholog_name] = []
            for descendent in leaves:
                ortholog_map[ortholog_name].append(leaf_map[descendent])

        return ortholog_map


class TaxonTest:
    def __init__(self, name):
        self.lineage = name
        self.taxon = name.split('; ')[-1]
        self.queries = list()
        self.classifieds = list()
        self.distances = dict()
        self.assignments = dict()
        self.taxonomic_tree = None
        self.intermediates_dir = ""
        self.temp_files_prefix = ""
        self.test_query_fasta = ""
        self.test_tax_ids_file = ""
        self.classifier_output = ""

    def get_optimal_assignment(self):
        if self.lineage.split('; ')[0] != "r__Root":
            self.lineage = "; ".join(["r__Root"] + self.lineage.split("; "))
        return optimal_taxonomic_assignment(self.taxonomic_tree, self.lineage)

    def summarise_taxon_test(self):
        summary_string = "Test for taxonomic lineage '" + self.lineage + "':\n" + \
                         "\tNumber of query sequences = " + str(len(self.queries)) + "\n" + \
                         "\tNumber of classified queries = " + str(len(self.classifieds)) + "\n"
        if self.assignments:
            for marker in self.assignments:
                summary_string += "Sequences classified as marker '" + marker + "':\n"
                for lineage in self.assignments[marker]:
                    summary_string += str(len(self.assignments[marker][lineage])) + "\t'" + lineage + "'\n"
        if self.taxonomic_tree:
            summary_string += "Optimal taxonomic assignment: '" + self.get_optimal_assignment() + "'\n"
        return summary_string

    def filter_assignments(self, target_marker: str):
        """
        Filters the assignments from TreeSAPP for the target marker.
        Off-target classifications are accounted for and reported.
        TaxonTest.classifieds only include the headers of the correctly annotated sequences

        :param target_marker:
        :return: None
        """
        off_targets = dict()
        num_classified = 0
        for marker in self.assignments:
            for lineage in self.assignments[marker]:
                classifieds = self.assignments[marker][lineage]
                num_classified += len(classifieds)
                if marker == target_marker:
                    self.classifieds += classifieds
                else:
                    if marker not in off_targets:
                        off_targets[marker] = list()
                    off_targets[marker] += classifieds
        if off_targets:
            for marker in off_targets:
                logging.warning(str(len(off_targets[marker])) + '/' + str(num_classified) +
                                " sequences were classified as " + marker + ":\n" +
                                "\t\n".join(off_targets[marker]) + "\n")
        return


class Evaluator(TreeSAPP):
    def __init__(self):
        super(Evaluator, self).__init__("evaluate")
        self.rank_depth_map = None
        self.ranks = list()
        self.markers = set()
        self.taxa_filter = dict()
        self.taxa_filter["Unclassified"] = 0
        self.taxa_filter["Classified"] = 0
        self.taxa_filter["Unique_taxa"] = 0
        self.taxa_tests = dict()
        self.classifications = dict()

        self.test_rep_taxa_fasta = ""
        self.performance_table = ""
        self.containment_table = ""
        self.recall_table = ""

        # Stage names only holds the required stages; auxiliary stages (e.g. RPKM, update) are added elsewhere
        self.stages = {0: ModuleFunction("lineages", 0),
                       1: ModuleFunction("classify", 1),
                       2: ModuleFunction("calculate", 2)}

    def decide_stage(self, args):
        """
        Bases the stage(s) to run on args.stage which is broadly set to either 'continue' or any other valid stage

        This function ensures all the required inputs are present for beginning at the desired first stage,
        otherwise, the pipeline begins at the first possible stage to continue and ends once the desired stage is done.

        If a
        :return: None
        """
        self.validate_continue(args)
        return

    def get_info(self):
        info_string = "Evaluator instance summary:\n"
        info_string += super(Evaluator, self).get_info() + "\n\t"
        info_string += "\n\t".join(["Accession-to-lineage map = " + self.acc_to_lin,
                                    "Clade-exclusion table = " + self.performance_table,
                                    "Target marker = " + str(self.ref_pkg.prefix)]) + "\n"

        return info_string

    def new_taxa_test(self, lineage) -> TaxonTest:
        """
        Creates a new TaxonTest instance for clade exclusion analysis

        :param lineage: The lineage being tested
        :return: A TaxonTest instance specific to lineage
        """
        # Determine the rank
        rank = self.ref_pkg.taxa_trie.resolved_to(lineage)
        if not rank:
            logging.error("Unable to find the rank the '{}' was resolved to.\n".format(lineage))
            sys.exit(5)

        if rank not in self.taxa_tests:
            self.taxa_tests[rank] = list()
        taxa_test_inst = TaxonTest(lineage)
        self.taxa_tests[rank].append(taxa_test_inst)
        return taxa_test_inst

    def delete_test(self, rank, lineage):
        i = 0
        for tti in self.taxa_tests[rank]:
            if re.match(tti.lineage, lineage):
                self.taxa_tests[rank].pop(i)
                break
            i += 1
        return

    def get_sensitivity(self, rank):
        total_queries = 0
        total_classified = 0
        if rank in self.taxa_tests:
            for tt in self.taxa_tests[rank]:  # type: TaxonTest
                total_queries += len(tt.queries)
                total_classified += len(tt.classifieds)
            return total_queries, total_classified, float(total_classified/total_queries)
        else:
            return 0, 0, 0.0

    def get_unique_taxa_tested(self, rank):
        taxa = set()
        if rank in self.taxa_tests:
            for tt in self.taxa_tests[rank]:
                taxa.add(tt.taxon)
            return taxa
        else:
            return None

    def summarize_rankwise_distances(self, rank):
        distals = list()
        pendants = list()
        tips = list()
        totals = list()
        n_dists = 0
        if rank in self.taxa_tests:
            for tt in self.taxa_tests[rank]:
                distals += tt.distances["distal"]
                pendants += tt.distances["pendant"]
                tips += tt.distances["tip"]
                n_dists += 1
            n_dists = len(distals)
            x = 0
            while x < n_dists:
                totals.append(sum([distals[x], pendants[x], tips[x]]))
                x += 1
            if len(pendants) != n_dists or len(tips) != n_dists:
                logging.error("Unequal number of values found between distal-, pendant- and tip-distances.\n")
                sys.exit(17)
            distance_summary = ["Rank\tType\tMean\tMedian\tVariance",
                                "\t".join([rank, "Distal",
                                           str(round(sum(distals) / float(n_dists), 4)),
                                           str(round(median(distals), 4)),
                                           str(round(float(var(distals)), 4))]),
                                "\t".join([rank, "Pendant",
                                           str(round(sum(pendants) / float(n_dists), 4)),
                                           str(round(median(pendants), 4)),
                                           str(round(float(var(pendants)), 4))]),
                                "\t".join([rank, "Tip",
                                           str(round(sum(tips) / float(n_dists), 4)),
                                           str(round(median(tips), 4)),
                                           str(round(float(var(tips)), 4))]),
                                "\t".join([rank, "Total",
                                           str(round(sum(totals) / float(n_dists), 4)),
                                           str(round(median(totals), 4)),
                                           str(round(float(var(totals)), 4))])]
            sys.stdout.write("\n".join(distance_summary) + "\n")
            return distals, pendants, tips
        else:
            return None, None, None

    def summarize_taxonomic_diversity(self):
        """
        Function for summarizing the taxonomic diversity of a reference dataset by rank

        :return: None
        """
        depth = 1  # Accumulator for parsing _RANK_DEPTH_MAP; not really interested in Cellular Organisms or Strains.
        info_str = ""
        while depth < 8:
            rank = self.rank_depth_map[depth]
            unique_taxa = self.get_unique_taxa_tested(rank)
            if unique_taxa:
                buffer = " "
                while len(rank) + len(str(len(unique_taxa))) + len(buffer) < 12:
                    buffer += ' '
                info_str += "\t" + rank + buffer + str(len(unique_taxa)) + "\n"
            else:
                pass
            depth += 1
        logging.info("Number of unique lineages tested:\n" + info_str)
        return

    def taxonomic_recall_tree(self):
        tree_summary_str = "Recall taxonomic tree:\n\n"
        taxa_tests_list = []
        observed = set()
        for rank in self.taxa_tests:
            taxa_tests_list += self.taxa_tests[rank]
        for tt in sorted(taxa_tests_list, key=lambda t: t.lineage):  # type: TaxonTest
            if len(tt.queries) == 0:
                continue
            lineage = tt.lineage.split("; ")
            x = 0
            while x < len(lineage)-1:
                if lineage[x] not in observed:
                    tree_summary_str += "\t"*x + lineage[x] + " = NA\n"
                    observed.add(lineage[x])
                x += 1
            tree_summary_str += "\t"*x
            tree_summary_str += "%s = %.1f\n" % (lineage[-1], float(100*len(tt.classifieds)/len(tt.queries)))
            observed.add(lineage[-1])

        logging.debug(tree_summary_str + "\n")
        logging.info("An alphabetically-sorted tree displaying recall of all taxa evaluated is in the log file.\n")
        return

    def taxonomic_recall_table(self):
        taxa_recall_str = "Rank\tTaxon\tRecall\n"
        for depth in sorted(self.rank_depth_map):
            rank = self.rank_depth_map[depth]
            if rank == "Cellular organisms":
                continue
            if rank not in self.classifications or len(self.classifications[rank]) == 0:
                continue
            for tt in sorted(self.taxa_tests[rank], key=lambda x: x.lineage):
                if len(tt.queries) == 0:
                    continue
                taxa_recall_str += "\t".join([rank, tt.taxon,
                                              str(round(float(100*len(tt.classifieds))/len(tt.queries), 2))]) + "\n"
        with open(self.recall_table, 'w') as recall_tbl_handler:
            recall_tbl_handler.write(taxa_recall_str)

        logging.info("Wrote taxon-wise recall to %s.\n" % self.recall_table)
        return

    def get_classification_performance(self):
        """
        Correct if: optimal_assignment == query_lineage

        :return: List of strings to be written to Evaluator.performance_table
        """
        std_out_report_string = ""
        clade_exclusion_strings = list()
        rank_assigned_dict = self.classifications

        sys.stdout.write("Rank-level performance of " + self.ref_pkg.prefix + ":\n")
        sys.stdout.write("\tRank\tQueries\tClassified\tCorrect\tD=1\tD=2\tD=3\tD=4\tD=5\tD=6\tD=7\n")

        for depth in sorted(self.rank_depth_map):
            rank = self.rank_depth_map[depth]
            if rank == "root":
                continue
            taxonomic_distance = dict()
            n_queries, n_classified, sensitivity = self.get_sensitivity(rank)
            for dist in range(0, 8):
                taxonomic_distance[dist] = EvaluateStats(self.ref_pkg.prefix, rank, dist)
                taxonomic_distance[dist].n_qs = n_queries
            std_out_report_string += "\t" + rank + "\t"
            if rank not in rank_assigned_dict or len(rank_assigned_dict[rank]) == 0:
                std_out_report_string += "0\t0\t\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n"
            else:
                acc = 0
                for assignments in rank_assigned_dict[rank]:
                    for classified in assignments:
                        acc += 1
                        optimal, query = assignments[classified]
                        if optimal == classified:
                            offset = 0
                        else:
                            offset = determine_offset(classified, optimal)
                        if offset > 7:
                            # This shouldn't be possible since there are no more than 7 taxonomic ranks
                            logging.error("Offset found to be greater than what is possible (" + str(offset) + ").\n" +
                                          "Classified: " + classified + "\n" +
                                          "Optimal: " + optimal + "\n" +
                                          "Query: " + query + "\n")
                            continue
                        for dist in taxonomic_distance:
                            eval_stats = taxonomic_distance[dist]  # type: EvaluateStats
                            if offset > dist:
                                if len(classified.split("; ")) > len(optimal.split("; ")):
                                    eval_stats.over_p += 1
                                else:
                                    eval_stats.under_p += 1
                            elif dist == offset:
                                eval_stats.correct += 1

                std_out_report_string += str(n_queries) + "\t" + str(n_classified) + "\t\t"

                classified_acc = 0
                for dist in sorted(taxonomic_distance.keys()):  # type: int
                    eval_stats = taxonomic_distance[dist]  # type: EvaluateStats
                    eval_stats.cumulative = eval_stats.correct + classified_acc
                    classified_acc += eval_stats.correct
                    if eval_stats.correct > 0:
                        if n_classified == 0:
                            logging.error("No sequences were classified at rank '" + rank +
                                          "' but optimal placements were pointed here. " +
                                          "This is a bug - please alert the developers!\n")
                            sys.exit(21)
                        else:
                            eval_stats.proportion = float(eval_stats.correct / n_classified)
                    else:
                        eval_stats.proportion = 0.0
                    clade_exclusion_strings.append(eval_stats.linear_stats())
                if classified_acc != n_classified:
                    logging.error("Discrepancy between classified sequences at each distance (" + str(classified_acc) +
                                  ") and total (" + str(n_classified) + ").\n")
                    sys.exit(15)

                std_out_report_string += '\t'.join([str(round(val.proportion*100, 1)) for val in
                                                    taxonomic_distance.values()]) + "\n"
                if (classified_acc*100)/n_queries > 101.0:
                    logging.error("Sum of proportional assignments at all distances is greater than 100.\n" +
                                  "\n".join(["Rank = " + rank,
                                             "Queries = " + str(n_queries),
                                             "Classified = " + str(n_classified),
                                             "Classifications = " + str(len(rank_assigned_dict[rank]))]) + "\n")
                    sys.exit(21)

        sys.stdout.write(std_out_report_string)

        return clade_exclusion_strings

    # TODO: Merge the two table-writing functions below
    def write_containment_table(self, containment_strings, tool):
        try:
            output_handler = open(self.containment_table, 'w')
        except IOError:
            logging.error("Unable to open " + self.containment_table + " for writing.\n")
            sys.exit(21)

        output_handler.write("# Input file for testing: " + self.input_sequences + "\n")
        trial_name = os.path.basename(self.output_dir)
        for line in containment_strings:
            # Line has a "\t" prefix already
            line = trial_name + "\t" + self.ref_pkg.prefix + "\t" + tool + line + "\n"
            output_handler.write(line)

        output_handler.close()
        return

    def write_performance_table(self, clade_exclusion_strings: list, tool):
        try:
            output_handler = open(self.performance_table, 'w')
        except IOError:
            logging.error("Unable to open " + self.performance_table + " for writing.\n")
            sys.exit(21)

        header = ["Trial", "Tool", "RefPkg", "Rank", "TaxDist", "Queries", "Correct", "Cumulative", "Over", "Under"]
        output_handler.write("# Input file for testing: " + self.input_sequences + "\n")
        output_handler.write("\t".join(header) + "\n")
        if self.output_dir[-1] == os.sep:
            trial_name = self.output_dir.split(os.sep)[-2]
        else:
            trial_name = self.output_dir.split(os.sep)[-1]
        for line in clade_exclusion_strings:
            line = trial_name + "\t" + tool + "\t" + line + "\n"
            output_handler.write(line)

        output_handler.close()
        return

    def clade_exclusion_outputs(self, lineage, rank, refpkg_name) -> TaxonTest:
        """
        Creates a TaxonTest instance that stores file paths and settings relevant to a clade exclusion analysis

        :param lineage: Taxonomic lineage which is going to be used in this clade exclusion analysis
        :param rank: The taxonomic rank (to which `lineage` belongs to) that is being excluded
        :param refpkg_name: Name (not code) of the ReferencePackage object that is to be used for clade exclusion
        :return: TaxonTest instance
        """
        # Refpkg input files in ts_evaluate.var_output_dir/refpkg_name/rank_tax/
        # Refpkg built in ts_evaluate.var_output_dir/refpkg_name/rank_tax/{refpkg_name}_{rank_tax}.gpkg/
        taxon = re.sub(r"([ /])", '_', lineage.split("; ")[-1])
        rank_tax = rank[0] + '_' + taxon

        test_obj = self.new_taxa_test(lineage)
        test_obj.intermediates_dir = self.var_output_dir + refpkg_name + os.sep + rank_tax + os.sep
        if not os.path.isdir(test_obj.intermediates_dir):
            os.makedirs(test_obj.intermediates_dir)

        logging.info("Classifications for the " + rank + " '" + taxon + "' put " + test_obj.intermediates_dir + "\n")
        test_obj.test_query_fasta = test_obj.intermediates_dir + rank_tax + ".fa"
        test_obj.test_tax_ids_file = test_obj.intermediates_dir + "tax_ids_" + refpkg_name + ".txt"
        test_obj.classifier_output = test_obj.intermediates_dir + "TreeSAPP_output" + os.sep
        return test_obj


class Layerer(TreeSAPP):
    def __init__(self):
        super(Layerer, self).__init__("layer")
        self.c_style_re = re.compile(".*_style.txt$")
        self.c_strip_re = re.compile(".*_strip.txt$")
        self.stages = {}
        self.annot_files = list()
        self.target_refpkgs = list()
        self.treesapp_output = ""
        self.colours_file = ""


class Abundance(TreeSAPP):
    def __init__(self):
        super(Abundance, self).__init__("abundance")
        self.stages = {}
        self.target_refpkgs = list()
        self.classified_nuc_seqs = ""
        self.classifications = ""
        self.fq_suffix_re = re.compile(r"([._-])+(pe|fq|fastq|fwd|R1)$")

    def check_arguments(self, args):
        ##
        # Define locations of files TreeSAPP outputs
        ##
        if len(glob(self.final_output_dir + "*_classified.fna")) < 1:
            logging.error("Unable to find classified ORF nucleotide sequences in '{}'.\n".format(self.final_output_dir))
            sys.exit(5)
        self.classified_nuc_seqs = glob(self.final_output_dir + "*_classified.fna")[0]
        if not os.path.isfile(self.classified_nuc_seqs):
            logging.error("Unable to find classified sequences FASTA file in %s.\n" % self.final_output_dir)
        self.classifications = self.output_dir + "final_outputs" + os.sep + "marker_contig_map.tsv"

        if not os.path.isdir(self.var_output_dir):
            os.makedirs(self.var_output_dir)
        self.validate_refpkg_dir(args.refpkg_dir)

        return


class Assigner(TreeSAPP):
    def __init__(self):
        """

        """
        super(Assigner, self).__init__("assign")
        self.reference_tree = None
        self.svc_filter = False
        self.aa_orfs_file = ""
        self.nuc_orfs_file = ""
        self.classified_aa_seqs = ""
        self.classified_nuc_seqs = ""
        self.composition = ""
        self.target_refpkgs = list()

        # Stage names only holds the required stages; auxiliary stages (e.g. RPKM, update) are added elsewhere
        self.stages = {0: ModuleFunction("orf-call", 0, self.predict_orfs),
                       1: ModuleFunction("clean", 1, self.clean),
                       2: ModuleFunction("search", 2, self.search),
                       3: ModuleFunction("align", 3, self.align),
                       4: ModuleFunction("place", 4, self.place),
                       5: ModuleFunction("classify", 5, self.classify)}

    def decide_stage(self, args):
        """
        Bases the stage(s) to run on args.stage which is broadly set to either 'continue' or any other valid stage

        This function ensures all the required inputs are present for beginning at the desired first stage,
        otherwise, the pipeline begins at the first possible stage to continue and ends once the desired stage is done.

        If a
        :return: None
        """
        if args.rpkm:
            if not args.reads:
                if args.reverse:
                    logging.error("File containing reverse reads provided but forward mates file missing!\n")
                    sys.exit(3)
                else:
                    logging.error("At least one FASTQ file must be provided if -rpkm flag is active!\n")
                    sys.exit(3)
            elif args.reads and not os.path.isfile(args.reads):
                logging.error("Path to forward reads ('%s') doesn't exist.\n" % args.reads)
                sys.exit(3)
            elif args.reverse and not os.path.isfile(args.reverse):
                logging.error("Path to reverse reads ('%s') doesn't exist.\n" % args.reverse)
                sys.exit(3)
            else:
                self.stages[len(self.stages)] = ModuleFunction("rpkm", len(self.stages))

        self.validate_continue(args)
        return

    def get_info(self):
        info_string = "Assigner instance summary:\n"
        info_string += super(Assigner, self).get_info() + "\n\t"
        info_string += "\n\t".join(["ORF protein sequences = " + self.aa_orfs_file,
                                    "Target reference packages = " + str(self.target_refpkgs),
                                    "Composition of input = " + self.composition]) + "\n"

        return info_string

    def predict_orfs(self, composition: str, num_threads: int) -> None:
        """
        Predict ORFs from the input FASTA file using Prodigal

        :param composition: Sample composition being either a single organism or a metagenome [single | meta]
        :param num_threads: The number of CPU threads to use
        :return: None
        """

        logging.info("Predicting open-reading frames using Prodigal... ")

        start_time = time.time()

        if num_threads > 1 and composition == "meta":
            # Split the input FASTA into num_threads files to run Prodigal in parallel
            split_files = fastx_split(self.input_sequences, self.output_dir, num_threads)
        else:
            split_files = [self.input_sequences]

        task_list = list()
        for fasta_chunk in split_files:
            chunk_prefix = self.var_output_dir + '.'.join(os.path.basename(fasta_chunk).split('.')[:-1])
            prodigal_command = [self.executables["prodigal"]]
            prodigal_command += ["-i", fasta_chunk]
            prodigal_command += ["-p", composition]
            prodigal_command += ["-a", chunk_prefix + "_ORFs.faa"]
            prodigal_command += ["-d", chunk_prefix + "_ORFs.fna"]
            prodigal_command += ["1>/dev/null", "2>/dev/null"]
            task_list.append(prodigal_command)

        num_tasks = len(task_list)
        if num_tasks > 0:
            cl_farmer = CommandLineFarmer("Prodigal -p " + composition, num_threads)
            cl_farmer.add_tasks_to_queue(task_list)

            cl_farmer.task_queue.close()
            cl_farmer.task_queue.join()

        tmp_prodigal_aa_orfs = glob(self.var_output_dir + self.sample_prefix + "*_ORFs.faa")
        tmp_prodigal_nuc_orfs = glob(self.var_output_dir + self.sample_prefix + "*_ORFs.fna")
        if not tmp_prodigal_aa_orfs or not tmp_prodigal_nuc_orfs:
            logging.error("Prodigal outputs were not generated:\n"
                          "Amino acid ORFs: " + ", ".join(tmp_prodigal_aa_orfs) + "\n" +
                          "Nucleotide ORFs: " + ", ".join(tmp_prodigal_nuc_orfs) + "\n")
            sys.exit(5)

        # Concatenate outputs
        if not os.path.isfile(self.aa_orfs_file) and not os.path.isfile(self.nuc_orfs_file):
            os.system("cat " + ' '.join(tmp_prodigal_aa_orfs) + " > " + self.aa_orfs_file)
            os.system("cat " + ' '.join(tmp_prodigal_nuc_orfs) + " > " + self.nuc_orfs_file)
            intermediate_files = list(tmp_prodigal_aa_orfs + tmp_prodigal_nuc_orfs + split_files)
            for tmp_file in intermediate_files:
                if tmp_file != self.input_sequences:
                    os.remove(tmp_file)

        logging.info("done.\n")

        end_time = time.time()
        hours, remainder = divmod(end_time - start_time, 3600)
        minutes, seconds = divmod(remainder, 60)
        logging.debug("\tProdigal time required: " +
                      ':'.join([str(hours), str(minutes), str(round(seconds, 2))]) + "\n")

        return

    def clean(self):
        return

    def search(self):
        return

    def align(self):
        return

    def place(self):
        return

    def classify(self):
        return


class PhyTrainer(TreeSAPP):
    def __init__(self):
        super(PhyTrainer, self).__init__("train")
        self.hmm_purified_seqs = ""  # If an HMM profile of the gene is provided its a path to FASTA with homologs
        self.placement_table = ""
        self.placement_summary = ""
        self.clade_ex_pquery_pkl = ""
        self.plain_pquery_pkl = ""
        self.feature_vector_file = ""
        self.conditions_file = ""
        self.tsne_plot = ""
        self.pkg_dbname_dict = dict()
        self.target_refpkgs = list()
        self.training_ranks = {}
        self.pqueries = {}

        # Stage names only holds the required stages; auxiliary stages (e.g. RPKM, update) are added elsewhere
        self.stages = {0: ModuleFunction("clean", 0),
                       1: ModuleFunction("search", 1),
                       2: ModuleFunction("lineages", 2),
                       3: ModuleFunction("place", 3),
                       4: ModuleFunction("train", 4),
                       5: ModuleFunction("update", 5)}

    def decide_stage(self, args):
        if not args.profile:
            self.change_stage_status("search", False)

        if args.acc_to_lin:
            self.acc_to_lin = args.acc_to_lin
            if os.path.isfile(self.acc_to_lin):
                self.change_stage_status("lineages", False)
            else:
                logging.error("Unable to find accession-lineage mapping file '{}'\n".format(self.acc_to_lin))
                sys.exit(3)

        if os.path.isfile(self.clade_ex_pquery_pkl) and os.path.isfile(self.plain_pquery_pkl):
            self.change_stage_status("place", False)

        self.validate_continue(args)
        return
    
    def set_file_paths(self) -> None:
        """
        Define the file path locations of treesapp train outputs

        :return: None
        """
        self.placement_table = os.path.join(self.final_output_dir, "placement_info.tsv")
        self.placement_summary = os.path.join(self.final_output_dir, "placement_trainer_results.txt")
        self.clade_ex_pquery_pkl = os.path.join(self.final_output_dir, "clade_exclusion_pqueries.pkl")
        self.plain_pquery_pkl = os.path.join(self.final_output_dir, "raw_refpkg_pqueries.pkl")
        self.formatted_input = os.path.join(self.var_output_dir, "clean", self.ref_pkg.prefix + "_formatted.fa")
        self.hmm_purified_seqs = os.path.join(self.var_output_dir, "search", self.ref_pkg.prefix + "_hmm_purified.fa")
        self.acc_to_lin = os.path.join(self.var_output_dir, "lineages", "accession_id_lineage_map.tsv")
        self.conditions_file = os.path.join(self.var_output_dir, "train", "conditions.npy")
        self.feature_vector_file = os.path.join(self.var_output_dir, "train", "examples.npy")
        return

    def get_info(self):
        info_string = "PhyTrainer instance summary:\n"
        info_string += super(PhyTrainer, self).get_info() + "\n\t"
        info_string += "\n\t".join(["Target reference packages = " + str(self.ref_pkg.prefix),
                                    "Taxonomy map: " + self.ref_pkg.lineage_ids,
                                    "Reference tree: " + self.ref_pkg.tree,
                                    "Reference FASTA: " + self.ref_pkg.msa,
                                    "Lineage map: " + str(self.acc_to_lin),
                                    "Ranks tested: " + ','.join(self.training_ranks.keys())]) + "\n"

        return info_string


class EvaluateStats:
    def __init__(self, refpkg: str, rank: str, dist: int):
        self.refpkg = refpkg
        self.rank = rank
        self.offset = dist
        self.n_qs = 0  # The number of query sequences used to evaluate this Rank
        self.correct = 0  # The number of query sequences that were correctly assigned at this taxonomic distance
        self.cumulative = 0  # Number of query sequences assigned at this taxonomic distance or less
        self.over_p = 0  # For the over-predictions: the assignment is wrong and more ranks are included than optimal
        self.under_p = 0  # For the under-predictions: the assignment is wrong and deeper than optimal
        self.proportion = 0.0

    def get_info(self):
        info_str = "# Evaluation stats for '" + self.refpkg + "' at rank: " + self.rank + "\n"
        info_str += "Taxonomic rank distance = " + str(self.offset) + "\n"
        info_str += "Queries\tCorrect\tCumulative\tOver-Pred\tUnder-Pred\n"
        info_str += "\t".join([str(val) for val in
                               [self.n_qs, self.correct, self.cumulative, self.over_p, self.under_p]]) + "\n"
        return info_str

    def linear_stats(self):
        stat_fields = [self.refpkg, self.rank, self.offset,
                       self.n_qs, self.correct, self.cumulative, self.over_p, self.under_p]
        return "\t".join([str(val) for val in stat_fields])

    def recall(self):
        return self.correct/self.n_qs

    def precision(self):
        return self.correct/sum([self.correct, self.over_p, self.under_p])
