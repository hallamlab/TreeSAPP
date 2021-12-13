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

from treesapp import phylo_seq
from treesapp import refpkg
from treesapp import fasta
from treesapp import utilities as ts_utils
from treesapp import lca_calculations
from treesapp import entrez_utils
from treesapp import logger
from treesapp.wrapper import estimate_ml_model

LOGGER = logging.getLogger(logger.logger_name())


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
        self.dir_path = ""
        self.inputs = []
        self.outputs = []
        self.function = func
        self.run = True

    def __str__(self) -> str:
        return "Stage '{}' is order number {} with run set to {}.".format(self.name, self.order, self.run)


def get_header_info(header_registry: dict, code_name=''):
    """

    :param header_registry: A dictionary of Header instances, indexed by numerical treesapp_id
    :param code_name: [OPTIONAL] The code_name of the reference package (marker gene/domain/family/protein)
    :return: Dictionary where keys are numerical treesapp_ids and values are EntrezRecord instances
    """
    LOGGER.info("Extracting information from headers... ")
    ref_records = dict()
    header_regexes = fasta.load_fasta_header_regexes(code_name)
    for treesapp_id in sorted(header_registry.keys(), key=int):  # type: str
        original_header = header_registry[treesapp_id].original
        header_format_re, header_db, header_molecule = fasta.get_header_format(original_header, header_regexes)
        sequence_info = header_format_re.match(original_header)
        seq_info_tuple = fasta.sequence_info_groups(sequence_info, header_db, original_header, header_regexes)

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
    LOGGER.info("done.\n")

    return ref_records


def dedup_records(ref_seqs: fasta.FASTA, ref_seq_records: dict) -> dict:
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
        LOGGER.debug("{}/{} unique versioned accessions were loaded for deduplication.\n"
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

        LOGGER.warning("The following sequences were removed during deduplication of Entrez records:\n\t" +
                       "\n\t".join(deduped_accessions) + "\n")

    return ref_seq_records


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
        self.ref_pkg = refpkg.ReferencePackage()
        self.classification_tbl_name = "classifications.tsv"

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
        self.ts_logger = logging.getLogger(logger.logger_name())

    def get_info(self):
        info_string = "Executables:\n\t" + "\n\t".join([k + ": " + v for k, v in self.executables.items()]) + "\n"

        for module_step in self.stages:
            info_string += str(self.stages[module_step])
        info_string += "\n\t".join(["TreeSAPP directory = " + self.treesapp_dir,
                                    "FASTA input sequences = " + self.input_sequences,
                                    "Formatted input = " + self.formatted_input,
                                    "Output directory = " + self.output_dir,
                                    "Input molecule type = " + self.molecule_type])
        return info_string

    def set_output_dirs(self) -> None:
        self.output_dir = ts_utils.validate_new_dir(self.output_dir)
        self.final_output_dir = self.output_dir + "final_outputs" + os.sep
        self.var_output_dir = self.output_dir + "intermediates" + os.sep
        return

    def set_sample_prefix(self, files: list) -> None:
        file_prefixes = []
        for f_name in files:
            file_name, _suffix1 = os.path.splitext(os.path.basename(f_name))
            if _suffix1 == ".gz":
                file_name, _suffix2 = os.path.splitext(file_name)
            file_prefixes.append(file_name)
        self.sample_prefix = '_'.join(file_prefixes)
        return

    def furnish_with_arguments(self, args) -> None:
        """
        Carries over the basic TreeSAPP arguments to the respective TreeSAPP-subclass.
        All auxiliary arguments are pushed to the TreeSAPP classes in check_module_arguments

        :param args: arguments from argparse.ParseArgs() with output, input and molecule attributes
        :return: None
        """
        if self.command != "info":
            self.output_dir = args.output
            self.set_output_dirs()
            if set(vars(args)).issuperset({"molecule", "input"}):
                self.set_sample_prefix(args.input)
                # Handle the fastx input file list
                if self.command == "create":
                    self.input_sequences = args.input
                else:
                    if len(args.input) > 1:
                        self.ts_logger.error("treesapp {} is unable to handle more than one fastx-input file.\n"
                                             "".format(self.command))
                        sys.exit(7)
                    self.input_sequences = args.input.pop(0)

                if args.molecule:
                    self.molecule_type = args.molecule

        if self.command != "colour" and "pkg_path" in vars(args):
            if len(args.pkg_path) > 1:
                self.ts_logger.warning("Multiple reference packages cannot be used by treesapp {}.\n"
                                       "Only the first one ({}) will be used.\n".format(self.command, args.pkg_path))
            args.pkg_path = args.pkg_path.pop(0)

        self.executables = self.find_executables(args)
        return

    def validate_refpkg_dir(self, refpkg_dir: str):
        if refpkg_dir:
            if not os.path.isdir(refpkg_dir):
                self.ts_logger.error(
                    "Directory containing reference packages ({}) does not exist.\n".format(refpkg_dir))
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

        for stage_order, stage in self.stages.items():  # type: (int, ModuleFunction)
            stage.dir_path = self.var_output_dir + stage.name + os.sep
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
                self.ts_logger.warning("Reclassify impossible as " + self.output_dir + " is missing input files.\n")
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
            self.ts_logger.error("Unable to find '{}' in {} stages.\n".format(name, self.command))
            sys.exit(3)
        else:
            self.ts_logger.warning("Unable to find '{}' stage. Returning a new one instead.\n".format(name))
            return ModuleFunction(name=name, order=stage_order + 1)

    def get_first_stage(self, optional_start=False) -> ModuleFunction:
        """
        Selects the earliest checkpoint that the workflow can be started from.
        Stages (ModuleFunction instances) are skipped by setting their 'run' attribute to False.
        This is the stage preceding the first stage whose directory does not exist.

        :return: None
        """
        for i, module in sorted(self.stages.items()):  # type: (int, ModuleFunction)
            if module.run:
                return module
        if optional_start:
            return self.stages[0]
        else:
            self.ts_logger.error("No stages are set to run!\n")
            sys.exit(3)

    def past_last_stage(self, stage_name=None):
        if stage_name:
            x = self.stage_lookup(name=stage_name).order
        else:
            x = self.current_stage.order

        while x < len(self.stages):
            if self.stages[x].run is True:
                return False
            x += 1
        return True

    def set_stage_dir(self) -> None:
        self.stage_output_dir = self.current_stage.dir_path
        if not os.path.isdir(self.stage_output_dir):
            os.mkdir(self.stage_output_dir)
        return

    def increment_stage_dir(self, checkpoint=None) -> None:
        """
        Updates self.stage_output_dir with the directory path of the next ModuleFunction.

        :param checkpoint: A string representing a ModuleFunction name.
        If provided, the 'current_stage' and 'stage_output_dir' are incremented if the order of the current_stage
        is greater than the checkpoint. Prevents run-away increments.
        :return: None
        """
        curr_order = self.current_stage.order
        if checkpoint:
            chk_stage = self.stage_lookup(name=checkpoint)
            if chk_stage.order > curr_order:
                return

        # Find the next stage
        while curr_order < max(self.stages.keys()):
            next_stage = self.stages[curr_order + 1]  # type: ModuleFunction
            if next_stage.run:
                self.current_stage = next_stage
                break
            curr_order += 1

        # Update the output directory for this stage
        self.set_stage_dir()
        return

    def stage_status(self, name):
        return self.stage_lookup(name).run

    def change_stage_status(self, name: str, new_status: bool):
        stage = self.stage_lookup(name)
        stage.run = new_status
        return

    def edit_stage_run_range(self, start: int, end=-1) -> None:
        """
        Alters ModuleFunction's run attribute within the self.stages dictionary.
        All stages between start and end have their 'run' attribute set to True.

        :param start: The order of the stage to start processing from
        :param end: The order of the stage to end processing
        :return: None
        """
        for x in sorted(self.stages, key=int):  # type: int
            stage = self.stages[x]  # type: ModuleFunction
            if end >= 0:
                if x < start or x > end:
                    stage.run = False
                else:
                    stage.run = True
            else:
                if x < start:
                    stage.run = False
                else:
                    stage.run = True

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
        # Set ModuleFunction run attribute according to the contents of the module's output directory
        for i, module in sorted(self.stages.items()):  # type: (int, ModuleFunction)
            if os.path.isdir(module.dir_path):
                if len(os.listdir(module.dir_path)) > 0:
                    module.run = False

        if args.overwrite:
            self.current_stage = self.get_first_stage(optional_start=True)
        else:
            self.current_stage = self.get_first_stage(optional_start=False)

        self.set_stage_dir()
        if "stage" not in vars(args) or args.stage == "continue":
            self.ts_logger.debug("Continuing with stage '{}'\n".format(self.current_stage.name))
            self.edit_stage_run_range(self.current_stage.order)
            return

        # Update the stage status
        desired_stage = self.stage_lookup(args.stage)
        desired_stage.run = True
        if desired_stage.order > self.current_stage.order:
            self.ts_logger.warning("Unable to run '{}' as it is ahead of the last completed stage.\n"
                                   "Continuing with stage '{}'\n".format(args.stage, self.current_stage.name))
            self.edit_stage_run_range(self.current_stage.order, desired_stage.order)
        elif desired_stage.order < self.current_stage.order:
            self.ts_logger.debug("Proceeding with '{}'\n".format(args.stage))
            self.edit_stage_run_range(desired_stage.order)
        else:
            # Proceed with running the desired stage
            self.ts_logger.debug("Proceeding with '{}'\n".format(args.stage))
            self.edit_stage_run_range(desired_stage.order, desired_stage.order)

        self.current_stage = self.get_first_stage(optional_start=True)
        return

    def find_sequence_molecule_type(self):
        if not self.molecule_type:
            self.molecule_type = fasta.guess_sequence_type(fastx_file=self.input_sequences)
            if not self.molecule_type:
                self.ts_logger.error("Unable to automatically detect the molecule type of '{}'.\n"
                                     "Please rerun with the argument '--molecule'.\n".format(self.input_sequences))
                sys.exit(7)
        return

    def find_executables(self, args) -> dict:
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

        if self.command in ["create", "update", "train", "evaluate", "info", "phylotu"]:
            dependencies += ["mmseqs", "mafft", "FastTree"]

            if hasattr(args, "od_seq") and args.od_seq:
                dependencies.append("OD-seq")

        if self.molecule_type == "rrna":
            dependencies += ["cmalign", "cmsearch", "cmbuild"]

        for dep in dependencies:
            exec_paths[dep] = ts_utils.fetch_executable_path(dep, self.treesapp_dir)

        return exec_paths

    def fetch_entrez_lineages(self, ref_seqs: fasta.FASTA, molecule: str, acc_to_taxid=None,
                              seqs_to_lineage=None) -> dict:
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
        ref_seqs.change_dict_keys("original")
        entrez_utils.load_ref_seqs(ref_seqs.fasta_dict, ref_seqs.header_registry, entrez_record_dict)
        self.ts_logger.debug("\tNumber of input sequences = {}\n".format(len(entrez_record_dict)))

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
            ref_leaf_nodes = phylo_seq.convert_entrez_to_tree_leaf_references(entrez_record_dict)
            self.ref_pkg.taxa_trie.feed_leaf_nodes(ref_leaf_nodes)
            entrez_utils.sync_record_and_hierarchy_lineages(ref_leaf_nodes, entrez_record_dict)
            self.ref_pkg.taxa_trie.validate_rank_prefixes()
            self.ref_pkg.taxa_trie.build_multifurcating_trie()

        if self.stage_status("lineages"):
            entrez_query_list, num_lineages_provided = entrez_utils.build_entrez_queries(entrez_record_dict)
            self.ts_logger.debug("\tNumber of queries =\t" + str(len(entrez_query_list)) + "\n")

            if len(entrez_query_list) >= 1:
                entrez_utils.map_accessions_to_lineages(entrez_query_list, self.ref_pkg.taxa_trie,
                                                        molecule, acc_to_taxid)
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

            # Write the accession-lineage mapping file - essential for training too
            ts_utils.write_dict_to_table(self.seq_lineage_map, self.acc_to_lin)
            self.increment_stage_dir()
        elif self.stage_status("lineages") is False and os.path.isfile(self.acc_to_lin):
            self.ts_logger.info("Reading cached lineages in '{}'... ".format(self.acc_to_lin))
            self.seq_lineage_map.update(entrez_utils.read_accession_taxa_map(self.acc_to_lin))
            self.ts_logger.info("done.\n")

        ref_seqs.change_dict_keys()
        return entrez_record_dict


class Updater(TreeSAPP):
    def __init__(self):
        super(Updater, self).__init__("update")
        self.seq_names_to_taxa = ""  # Optional user-provided file mapping query sequence contigs to lineages
        self.lineage_map_file = ""  # File that is passed to create() containing lineage info for all sequences
        self.treesapp_output = ""  # Path to the TreeSAPP output directory - modified by args
        self.assignment_table = ""  # Path to the classifications.tsv file written by treesapp assign
        self.combined_fasta = ""  # Holds the newly identified candidate reference sequences and the original ref seqs
        self.old_ref_fasta = ""  # Contains only the original reference sequences
        self.cluster_input = ""  # Used only if resolve is True
        self.clusters_prefix = ""  # Used only if resolve is True
        self.updated_refpkg_path = ""
        self.training_dir = ""
        # self.rank_depth_map = None
        self.prop_sim = 1.0
        self.min_length = 0  # The minimum sequence length for a classified sequence to be included in the refpkg
        self.updated_refpkg = refpkg.ReferencePackage()

        # Stage names only holds the required stages; auxiliary stages (e.g. RPKM, update) are added elsewhere
        self.stages = {0: ModuleFunction("lineages", 0),
                       1: ModuleFunction("rebuild", 1),
                       2: ModuleFunction("train", 2)}

    def decide_stage(self, args):
        """
        Bases the stage(s) to run on args.stage which is broadly set to either 'continue' or any other valid stage

        This function ensures all the required inputs are present for beginning at the desired first stage,
        otherwise, the pipeline begins at the first possible stage to continue and ends once the desired stage is done.

        :return: None
        """
        self.validate_continue(args)
        self.acc_to_lin = self.var_output_dir + os.sep + "accession_id_lineage_map.tsv"

        if not self.input_sequences:
            self.change_stage_status("train", False)

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

    def update_refpkg_fields(self, output_dir=None) -> None:
        """
        Using the original ReferencePackage as a template modify the following updated ReferencePackage attributes:
1. original creation date
2. update date
3. code name
4. description
        :return: None
        """
        if not output_dir:
            output_dir = self.final_output_dir
        # Change the creation and update dates, code name and description
        self.updated_refpkg.change_file_paths(output_dir)
        self.updated_refpkg.date = self.ref_pkg.date
        self.updated_refpkg.update = dt.now().strftime("%Y-%m-%d")
        self.updated_refpkg.refpkg_code = self.ref_pkg.refpkg_code
        self.updated_refpkg.description = self.ref_pkg.description
        self.updated_refpkg.pickle_package()

        self.ts_logger.info("Summary of the updated reference package:\n" + self.updated_refpkg.get_info() + "\n")

        self.ts_logger.debug("\tNew sequences  = " + str(self.updated_refpkg.num_seqs - self.ref_pkg.num_seqs) + "\n" +
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
        self.stages = {0: ModuleFunction("deduplicate", 0),
                       1: ModuleFunction("search", 1),
                       2: ModuleFunction("lineages", 2),
                       3: ModuleFunction("clean", 3),
                       4: ModuleFunction("cluster", 4),
                       5: ModuleFunction("build", 5),
                       6: ModuleFunction("evaluate", 6),
                       7: ModuleFunction("support", 7),
                       8: ModuleFunction("train", 8),
                       9: ModuleFunction("update", 9)}

    def decide_stage(self, args):
        """
        Bases the stage(s) to run on args.stage which is broadly set to either 'continue' or any other valid stage

        This function ensures all the required inputs are present for beginning at the desired first stage,
        otherwise, the pipeline begins at the first possible stage to continue and ends once the desired stage is done.

        :return: None
        """
        self.validate_continue(args)
        if not args.dedup:
            self.change_stage_status("deduplicate", False)

        if not args.profile:
            self.change_stage_status("search", False)
        if args.pc:
            self.edit_stage_run_range(self.stage_lookup("update").order)
        # TODO: Allow users to provide sequence-lineage maps for a subset of the query sequences
        if args.acc_to_lin:
            self.acc_to_lin = args.acc_to_lin
            if os.path.isfile(self.acc_to_lin):
                self.change_stage_status("lineages", False)
            else:
                self.ts_logger.error("Unable to find accession-lineage mapping file '{}'\n".format(self.acc_to_lin))
                sys.exit(3)
        else:
            self.acc_to_lin = self.var_output_dir + os.sep + "accession_id_lineage_map.tsv"

        if args.bootstraps == 0:
            self.change_stage_status("support", False)
        if not args.fast:
            self.change_stage_status("evaluate", False)
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

    def combine_input_files(self) -> None:
        cat_fastx_file = os.path.join(self.var_output_dir, self.sample_prefix + ".fa")
        ts_utils.concatenate_files(input_files=self.input_sequences,
                                   output_path=cat_fastx_file)
        self.input_sequences = cat_fastx_file
        return

    def determine_model(self, ref_pkg: refpkg.ReferencePackage, estimate=False) -> None:
        if ref_pkg.tree_tool == "FastTree":
            if ref_pkg.molecule == "prot":
                evo_model = "LG+G4"
            elif ref_pkg.molecule == "rrna" or ref_pkg.molecule == "dna":
                evo_model = "GTR+G"
            else:
                self.ts_logger.error("Unrecognized reference package molecule type: '{}'.\n".format(ref_pkg.molecule))
                sys.exit(3)
            if ref_pkg.sub_model:
                self.ts_logger.warning("Model provided '{}' will be ignored when FastTree is used to infer phylogeny.\n"
                                       "".format(ref_pkg.sub_model))
        elif ref_pkg.tree_tool == "RAxML-NG":
            if estimate:
                evo_model = estimate_ml_model(modeltest_exe=self.executables["ModelTest-NG"],
                                              msa=ref_pkg.f__msa, output_prefix=self.phy_dir, molecule=ref_pkg.molecule)
            elif not ref_pkg.sub_model:
                if ref_pkg.molecule == "prot":
                    evo_model = "LG+G4"
                else:
                    evo_model = "GTR+G"
            else:
                self.ts_logger.debug("Using specified RAxML-NG-compatible model: '{}'.\n".format(ref_pkg.sub_model))
                evo_model = ref_pkg.sub_model
        else:
            self.ts_logger.error("Unexpected phylogenetic inference tool: '{}'.\n".format(ref_pkg.tree_tool))
            sys.exit(3)

        ref_pkg.sub_model = evo_model
        return

    def print_terminal_commands(self):
        self.ts_logger.info(
            "\nTo integrate this package for use in TreeSAPP you must copy {0} to a directory containing other"
            " reference packages you want to analyse. This may be in {1}/data/ or elsewhere\n"
            "".format(self.ref_pkg.f__pkl, self.treesapp_dir))
        return

    def overcluster_warning(self, pre_count: int, post_count: int) -> None:
        low_count = 100
        proportion = 0.5
        if pre_count > low_count > post_count and post_count / pre_count < proportion:
            self.ts_logger.warning("Clustering at {} similarity removed >{}% of input sequences -"
                                   " consider increasing this threshold.\n"
                                   "You have five seconds to cancel and restart.\n".format(self.ref_pkg.pid,
                                                                                           100 * (1 - proportion)))
            time.sleep(5)
        return


class Purity(TreeSAPP):
    def __init__(self):
        """Class instance for validating a reference package's purity with a reference sequence database"""
        super(Purity, self).__init__("purity")
        self.classifications = ""
        self.metadata_file = ""
        self.assignments = None
        self.assign_jplace_file = ""
        self.stages = {0: ModuleFunction("clean", 0),
                       1: ModuleFunction("assign", 1),
                       2: ModuleFunction("summarize", 2)}

    def check_purity_arguments(self, args):
        self.find_sequence_molecule_type()
        self.ref_pkg.f__pkl = args.pkg_path
        self.ref_pkg.slurp()
        self.refpkg_dir = os.path.dirname(os.path.realpath(self.ref_pkg.f__pkl))

        self.formatted_input = self.stage_lookup("clean").dir_path + self.sample_prefix + "_formatted.fasta"

        ##
        # Define locations of files TreeSAPP outputs
        ##
        self.classifications = os.path.join(self.stage_lookup("assign").dir_path,
                                            "final_outputs",
                                            self.classification_tbl_name)
        self.assign_jplace_file = os.path.join(self.stage_lookup("assign").dir_path,
                                               "iTOL_output",
                                               self.ref_pkg.prefix,
                                               self.ref_pkg.prefix + "_complete_profile.jplace")
        self.metadata_file = args.extra_info

        if not os.path.isdir(self.var_output_dir):
            os.makedirs(self.var_output_dir)

        return

    def decide_stage(self, args):
        """
        Bases the stage(s) to run on args.stage which is broadly set to either 'continue' or any other valid stage

        This function ensures all the required inputs are present for beginning at the desired first stage,
        otherwise, the pipeline begins at the first possible stage to continue and ends once the desired stage is done.

        :return: None
        """
        self.validate_continue(args)
        return

    def summarize_groups_assigned(self, ortholog_map: dict, metadata=None) -> None:
        unique_orthologs = dict()
        tree_leaves = self.ref_pkg.generate_tree_leaf_references_from_refpkg()
        for _ref_pkg, info in self.assignments.items():
            for lineage in info:
                for seq_name in info[lineage]:
                    bits_pieces = seq_name.split('_')
                    if bits_pieces[0] not in unique_orthologs:
                        unique_orthologs[bits_pieces[0]] = []
                    unique_orthologs[bits_pieces[0]].append('_'.join(bits_pieces[1:]))

        # Summarize the counts in the log file
        summary_str = "Ortholog\tHits\tLeaves\tTree-coverage\tDescription\n" + '-' * 80 + "\n"
        for og_name in sorted(ortholog_map, key=lambda x: len(ortholog_map[x])):
            n_leaves = len(set(ortholog_map[og_name]))
            perc_coverage = round(float((n_leaves * 100) / len(tree_leaves)), 1)
            try:
                desc = metadata[og_name].de
            except (TypeError, KeyError):
                desc = "NA"
            summary_str += "\t".join([og_name,
                                      str(len(unique_orthologs[og_name])), str(n_leaves), str(perc_coverage),
                                      desc]) + "\n"
        self.ts_logger.info(summary_str + "\n")

        return

    def load_metadata(self) -> dict:
        metadat_dict = dict()
        xtra_dat = namedtuple("xtra_dat", ["db_id", "ac", "de"])
        if not os.path.isfile(self.metadata_file):
            self.ts_logger.error("Extra information file '" + self.metadata_file + "' doesn't exist!\n")
            sys.exit(3)
        try:
            metadata_handler = open(self.metadata_file, 'r')
        except IOError:
            self.ts_logger.error("Unable to open extra information file '" + self.metadata_file + "' for reading!\n")
            sys.exit(3)
        for line in metadata_handler:
            try:
                db_id, accession, desc = line.strip().split("\t")
            except ValueError:
                self.ts_logger.error(
                    "Bad format for '" + self.metadata_file + "'. Three tab-separated fields expected.\n" +
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
        for p_query in p_queries:  # type: phylo_seq.PQuery
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
    def __init__(self, name: str):
        self.lineage = name
        self.taxon = name.split('; ')[-1]
        self.classifier = ""
        self.refpkg_path = ""
        self.logger = logging.getLogger(logger.logger_name())

        # Collections
        self.queries = list()
        self.classifieds = list()
        self.distances = dict()
        self.assignments = dict()
        self.taxonomic_tree = None

        # Directories
        self.intermediates_dir = ""
        self.classifications_root = ""

        # Output files
        self.test_query_fasta = ""
        self.classification_table = ""
        return

    def get_optimal_assignment(self):
        if self.lineage.split('; ')[0] != "r__Root":
            self.lineage = "; ".join(["r__Root"] + self.lineage.split("; "))
        return lca_calculations.optimal_taxonomic_assignment(self.taxonomic_tree, self.lineage)

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
                self.logger.warning(str(len(off_targets[marker])) + '/' + str(num_classified) +
                                    " sequences were classified as " + marker + ":\n" +
                                    "\t\n".join(off_targets[marker]) + "\n")
        return

    def clade_exclusion_outputs(self, ref_pkg: refpkg.ReferencePackage, output_dir: str, tool: str) -> None:
        """
        Creates a TaxonTest instance that stores file paths and settings relevant to a clade exclusion analysis

        :param ref_pkg: The ReferencePackage object that is to be used for clade exclusion
        :param output_dir: Root directory for the various outputs for this TaxonTest
        :param tool: Name of the tool used for classifying query sequences: 'graft', 'diamond' or 'treesapp'
        :return: TaxonTest instance
        """
        taxon_path = re.sub(r"([ /])", '_', self.taxon)
        taxon_path = re.sub(r"([()'\[\]])", '', taxon_path)

        self.intermediates_dir = os.path.join(output_dir, ref_pkg.prefix, taxon_path) + os.sep
        if not os.path.isdir(self.intermediates_dir):
            os.makedirs(self.intermediates_dir)

        self.test_query_fasta = self.intermediates_dir + taxon_path + ".fa"
        self.classifications_root = self.intermediates_dir + tool + "_output" + os.sep

        if tool in ["graftm", "diamond"]:
            self.classification_table = self.classifications_root + taxon_path + os.sep + taxon_path + "_read_tax.tsv"
            self.refpkg_path = self.intermediates_dir + ref_pkg.prefix + '_' + taxon_path + ".gpkg"
        else:
            self.classification_table = self.classifications_root + "final_outputs" + os.sep + "classifications.tsv"
            self.refpkg_path = os.path.join(self.intermediates_dir, ref_pkg.prefix + ref_pkg.refpkg_suffix)

        return


class Evaluator(TreeSAPP):
    def __init__(self):
        super(Evaluator, self).__init__("evaluate")
        self.rank_depth_map = None
        self.min_seq_length = 0
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

        :return: None
        """
        self.acc_to_lin = self.stage_lookup("lineages").dir_path + "accession_id_lineage_map.tsv"
        if args.acc_to_lin:
            if os.path.isfile(args.acc_to_lin):
                self.acc_to_lin = args.acc_to_lin
                self.change_stage_status("lineages", False)
            else:
                self.ts_logger.error(
                    "Unable to find accession-lineage mapping file '{}'\n".format(args.acc_to_lin))
                sys.exit(3)
        elif os.path.isfile(self.acc_to_lin):
            self.ts_logger.info("An accession-lineage mapping file from a previous run ('{}') was found"
                                " and will attempt to be used.\n".format(self.acc_to_lin))
            self.change_stage_status("lineages", False)
        self.validate_continue(args)
        return

    def get_info(self):
        info_string = "Evaluator instance summary:\n"
        info_string += super(Evaluator, self).get_info() + "\n\t"
        info_string += "\n\t".join(["Accession-to-lineage map = " + self.acc_to_lin,
                                    "Clade-exclusion table = " + self.performance_table,
                                    "Target marker = " + str(self.ref_pkg.prefix)]) + "\n"

        return info_string

    def new_taxa_test(self, lineage: str, tool: str) -> TaxonTest:
        """
        Creates a new TaxonTest instance for clade exclusion analysis

        :param lineage: The lineage being tested
        :param tool: Name of the tool used for classifying query sequences: 'graft', 'diamond' or 'treesapp'
        :return: A TaxonTest instance specific to lineage
        """
        # Determine the rank
        rank = self.ref_pkg.taxa_trie.resolved_to(lineage)
        if not rank:
            self.ts_logger.error("Unable to find the rank the '{}' was resolved to.\n".format(lineage))
            sys.exit(5)

        if rank not in self.taxa_tests:
            self.taxa_tests[rank] = list()
        taxa_test_inst = TaxonTest(lineage)
        self.taxa_tests[rank].append(taxa_test_inst)

        taxa_test_inst.clade_exclusion_outputs(output_dir=self.var_output_dir, ref_pkg=self.ref_pkg, tool=tool)

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
            return total_queries, total_classified, float(total_classified / total_queries)
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
                self.ts_logger.error("Unequal number of values found between distal-, pendant- and tip-distances.\n")
                sys.exit(17)
            distance_summary = ["Rank\tType\tMean\tMedian\tVariance",
                                "\t".join([rank, "Distal",
                                           str(round(sum(distals) / float(n_dists), 4)),
                                           str(round(ts_utils.median(distals), 4)),
                                           str(round(float(var(distals)), 4))]),
                                "\t".join([rank, "Pendant",
                                           str(round(sum(pendants) / float(n_dists), 4)),
                                           str(round(ts_utils.median(pendants), 4)),
                                           str(round(float(var(pendants)), 4))]),
                                "\t".join([rank, "Tip",
                                           str(round(sum(tips) / float(n_dists), 4)),
                                           str(round(ts_utils.median(tips), 4)),
                                           str(round(float(var(tips)), 4))]),
                                "\t".join([rank, "Total",
                                           str(round(sum(totals) / float(n_dists), 4)),
                                           str(round(ts_utils.median(totals), 4)),
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
        self.ts_logger.info("Number of unique lineages tested:\n" + info_str)
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
            while x < len(lineage) - 1:
                if lineage[x] not in observed:
                    tree_summary_str += "\t" * x + lineage[x] + " = NA\n"
                    observed.add(lineage[x])
                x += 1
            tree_summary_str += "\t" * x
            tree_summary_str += "%s = %.1f\n" % (lineage[-1], float(100 * len(tt.classifieds) / len(tt.queries)))
            observed.add(lineage[-1])

        self.ts_logger.debug(tree_summary_str + "\n")
        self.ts_logger.info(
            "An alphabetically-sorted tree displaying recall of all taxa evaluated is in the log file.\n")
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
                                              str(round(float(100 * len(tt.classifieds)) / len(tt.queries), 2))]) + "\n"
        with open(self.recall_table, 'w') as recall_tbl_handler:
            recall_tbl_handler.write(taxa_recall_str)

        self.ts_logger.info("Wrote taxon-wise recall to %s.\n" % self.recall_table)
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
                            offset = lca_calculations.determine_offset(classified, optimal)
                        if offset > 7:
                            # This shouldn't be possible since there are no more than 7 taxonomic ranks
                            self.ts_logger.error(
                                "Offset found to be greater than what is possible (" + str(offset) + ").\n" +
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
                            self.ts_logger.error("No sequences were classified at rank '" + rank +
                                                 "' but optimal placements were pointed here. " +
                                                 "This is a bug - please alert the developers!\n")
                            sys.exit(21)
                        else:
                            eval_stats.proportion = float(eval_stats.correct / n_classified)
                    else:
                        eval_stats.proportion = 0.0
                    clade_exclusion_strings.append(eval_stats.linear_stats())
                if classified_acc != n_classified:
                    self.ts_logger.error(
                        "Discrepancy between classified sequences at each distance (" + str(classified_acc) +
                        ") and total (" + str(n_classified) + ").\n")
                    sys.exit(15)

                std_out_report_string += '\t'.join([str(round(val.proportion * 100, 1)) for val in
                                                    taxonomic_distance.values()]) + "\n"
                if (classified_acc * 100) / n_queries > 101.0:
                    self.ts_logger.error("Sum of proportional assignments at all distances is greater than 100.\n" +
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
            self.ts_logger.error("Unable to open " + self.containment_table + " for writing.\n")
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
            self.ts_logger.error("Unable to open " + self.performance_table + " for writing.\n")
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


class Abundance(TreeSAPP):
    def __init__(self):
        super(Abundance, self).__init__("abundance")
        self.target_refpkgs = dict()
        self.ref_nuc_seqs = ""
        self.classifications = ""
        self.aln_file = ""
        self.append_abundance = True
        self.fq_suffix_re = re.compile(r"([._-])+(pe|fq|fastq|fwd|R1|1)$")
        self.idx_extensions = ["amb", "ann", "bwt", "pac", "sa"]
        self.stages = {0: ModuleFunction("align_map", 0),
                       1: ModuleFunction("sam_sum", 1),
                       2: ModuleFunction("summarise", 2)}
        return

    def check_arguments(self, args):
        ##
        # Define locations of files TreeSAPP outputs
        ##
        if len(glob(self.final_output_dir + "*_classified.fna")) < 1:
            self.ts_logger.error(
                "Unable to find classified ORF nucleotide sequences in '{}'.\n".format(self.final_output_dir))
            sys.exit(5)
        self.ref_nuc_seqs = ts_utils.match_file(os.path.join(self.var_output_dir, "orf-call", "*_ORFs.fna"))

        self.classifications = self.final_output_dir + self.classification_tbl_name

        if not os.path.isdir(self.var_output_dir):
            os.makedirs(self.var_output_dir)

        # Set the directory paths for each stage. Usually done in TreeSAPP.check_previous_output() but can't use here.
        for stage_order, stage in self.stages.items():  # type: (int, ModuleFunction)
            stage.dir_path = self.var_output_dir + stage.name + os.sep

        # Remove the directory containing SAM and BWA index files
        if os.path.isdir(self.stage_lookup("align_map").dir_path) and args.overwrite:
            self.ts_logger.debug(
                "Removing directory with BWA outputs '{}'.\n".format(self.stage_lookup("align_map").dir_path))
            rmtree(self.stage_lookup("align_map").dir_path)

        self.aln_file = self.stage_lookup("align_map").dir_path + \
                        '.'.join(os.path.basename(self.ref_nuc_seqs).split('.')[0:-1]) + ".sam"

        # Ensure all the FASTQ file paths are valid for the forward reads
        for reads_file in args.reads:
            if not os.path.isfile(reads_file):
                self.ts_logger.error("Unable to find forward reads file '{}'.\n".format(reads_file))
                sys.exit(1)

        if len(args.reverse) > 0:
            if len(args.reads) != len(args.reverse):
                self.ts_logger.error("Number of fastq files differs between reads ({}) and reverse ({}) arguments!.\n"
                                     "".format(len(args.reads), len(args.reverse)))
                sys.exit(3)
            # Ensure all the FASTQ file paths are valid for the reverse reads
            for reads_file in args.reverse:
                if not os.path.isfile(reads_file):
                    self.ts_logger.error("Unable to find reverse reads files '{}'.\n".format(reads_file))
                    sys.exit(5)

        if args.report == "append":
            self.append_abundance = True
        else:
            self.append_abundance = False

        return

    def decide_stage(self, args):
        """
        Bases the stage(s) to run on args.stage which is broadly set to either 'continue' or any other valid stage

        This function ensures all the required inputs are present for beginning at the desired first stage,
        otherwise, the pipeline begins at the first possible stage to continue and ends once the desired stage is done.

        :return: None
        """
        self.validate_continue(args)

        # Not necessary for either 'update'
        for reads_file in args.reads:
            self.strip_file_to_sample_name(reads_file)
            if not os.path.isfile(self.stage_lookup("align_map").dir_path + self.sample_prefix + ".sam"):
                self.change_stage_status("align_map", True)
                break

        return

    def delete_intermediates(self, clean_up: bool) -> None:
        if not clean_up:
            return

        files_to_be_deleted = glob(self.stage_lookup("align_map").dir_path + "*.sam")
        files_to_be_deleted += glob(self.stage_lookup("align_map").dir_path + "*.stderr")
        for ext in self.idx_extensions:
            files_to_be_deleted += glob(self.stage_lookup("align_map").dir_path + "*." + ext)

        if self.aln_file not in files_to_be_deleted:
            files_to_be_deleted.append(self.aln_file)

        for file_path in set(files_to_be_deleted):
            if os.path.exists(file_path):
                os.remove(file_path)
        return

    def strip_file_to_sample_name(self, file_path) -> None:
        file_name, suffix = os.path.splitext(os.path.basename(file_path))
        if suffix == ".gz":
            file_name, suffix = os.path.splitext(file_name)
        file_prefix = '.'.join(file_name.split('.'))
        self.sample_prefix = self.fq_suffix_re.sub('', file_prefix)
        return

    def fetch_refpkgs_used(self, refpkg_dir=None) -> None:
        if refpkg_dir:
            self.target_refpkgs = refpkg.gather_ref_packages(refpkg_data_dir=refpkg_dir)
            self.refpkg_dir = refpkg_dir
        else:
            # Load the reference packages in intermediates/, or the reference packages provided through refpkg_dir
            self.target_refpkgs = refpkg.load_refpkgs_from_assign_output(self.var_output_dir)
            if not self.target_refpkgs:
                self.target_refpkgs = refpkg.gather_ref_packages(self.refpkg_dir)
            else:
                self.refpkg_dir = self.var_output_dir
        if os.sep != self.refpkg_dir[-1]:
            self.refpkg_dir += os.sep
        return


class EvaluateStats:
    def __init__(self, ref_pkg: str, rank: str, dist: int):
        self.ref_pkg = ref_pkg
        self.rank = rank
        self.offset = dist
        self.n_qs = 0  # The number of query sequences used to evaluate this Rank
        self.correct = 0  # The number of query sequences that were correctly assigned at this taxonomic distance
        self.cumulative = 0  # Number of query sequences assigned at this taxonomic distance or less
        self.over_p = 0  # For the over-predictions: the assignment is wrong and more ranks are included than optimal
        self.under_p = 0  # For the under-predictions: the assignment is wrong and deeper than optimal
        self.proportion = 0.0

    def get_info(self):
        info_str = "# Evaluation stats for '" + self.ref_pkg + "' at rank: " + self.rank + "\n"
        info_str += "Taxonomic rank distance = " + str(self.offset) + "\n"
        info_str += "Queries\tCorrect\tCumulative\tOver-Pred\tUnder-Pred\n"
        info_str += "\t".join([str(val) for val in
                               [self.n_qs, self.correct, self.cumulative, self.over_p, self.under_p]]) + "\n"
        return info_str

    def linear_stats(self):
        stat_fields = [self.ref_pkg, self.rank, self.offset,
                       self.n_qs, self.correct, self.cumulative, self.over_p, self.under_p]
        return "\t".join([str(val) for val in stat_fields])

    def recall(self):
        return self.correct / self.n_qs

    def precision(self):
        return self.correct / sum([self.correct, self.over_p, self.under_p])
