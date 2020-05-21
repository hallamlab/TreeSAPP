import logging
import os
import re
import sys
from glob import glob
from shutil import copy
import json

from ete3 import Tree

from treesapp.phylo_seq import TreeLeafReference
from treesapp.entish import annotate_partition_tree
from treesapp.external_command_interface import launch_write_command
from treesapp.fasta import get_headers, read_fasta_to_dict, write_new_fasta
from treesapp.taxonomic_hierarchy import TaxonomicHierarchy
from treesapp.utilities import swap_tree_names, get_hmm_length
from treesapp.wrapper import model_parameters
from treesapp import __version__ as ts_version


class ReferencePackage:
    def __init__(self, refpkg_name=""):
        self.prefix = refpkg_name
        self.refpkg_code = "Z1111"  # AKA denominator

        # These are files (with '_f' suffix) and their respective data (read with file.readlines())
        self.f__json = self.prefix + "_build.json"  # Path to the JSON reference package file
        self.msa = ""  # Reference MSA FASTA
        self.f__msa = self.prefix + ".fa"
        self.profile = ""
        self.f__profile = self.prefix + ".hmm"  # HMM file
        self.search_profile = ""
        self.f__search_profile = self.prefix + '_' + "search.hmm"  # profile HMM that has been dereplicated
        self.tree = ""  # Reference tree
        self.f__tree = self.prefix + ".nwk"
        self.boot_tree = ""  # Reference tree with support values
        self.f__boot_tree = self.prefix + "_bipart.nwk"
        self.model_info = ""
        self.f__model_info = self.prefix + "_epa.model"  # RAxML-NG --evaluate model file
        self.lineage_ids = dict()  # Reference sequence lineage map

        # These are metadata values
        self.ts_version = ts_version
        self.sub_model = ""  # EPA-NG compatible substitution model
        self.date = ""  # Date the reference package was created
        self.update = ""  # Date the reference package was last updated
        self.num_seqs = 0  # Number of reference sequences in the MSA, phylogeny
        self.profile_length = 0  # LENG of the profile HMM (not search profile)
        self.molecule = ""  # nucleotide or amino acid sequence
        self.kind = ""  # Taxonomic or functional anchor?
        self.tree_tool = ""  # Software used for inferring the reference phylogeny
        self.description = ""
        self.pid = 1.0  # Proportional sequence similarity inputs were clustered at
        self.pfit = []  # Parameters for the polynomial regression function
        self.cmd = ""  # The command used for building the reference package

        # These are attributes only used during runtime
        self.core_ref_files = [self.f__msa, self.f__profile, self.f__search_profile,
                               self.f__tree, self.f__boot_tree, self.f__model_info]
        self.taxa_trie = TaxonomicHierarchy()

    def __iter__(self):
        for attr, value in self.__dict__.items():
            yield attr, value

    def band(self):
        """
        Writes a JSON file containing all of the reference package files and metadata

        :return:
        """
        refpkg_dict = {}
        if len(self.f__json) == 0:
            logging.error("ReferencePackage.f__json not set. ReferencePackage band cannot be completed.\n")
            sys.exit(11)

        # Read the all of the individual reference package files that are available (e.g. MSA, HMM, phylogeny)
        for a, v in self.__iter__():
            if a.startswith('f__') and os.path.isfile(v):
                dest = a.lstrip("f__")
                if dest in self.__dict__:
                    with open(v) as fh:
                        self.__dict__[dest] = fh.readlines()

        try:
            refpkg_handler = open(self.f__json, 'w')
        except IOError:
            logging.error("Unable to open reference package JSON file '{}' for writing.\n".format(self.f__json))
            sys.exit(11)

        non_primitives = ["core_ref_files", "taxa_trie"]
        for a, v in self.__iter__():
            if a not in non_primitives:
                refpkg_dict[a] = v

        json.dump(obj=refpkg_dict, fp=refpkg_handler)

        refpkg_handler.close()

        return

    def write_refpkg_component(self, file_name: str, text) -> None:
        """
        Writes text to a file.

        :param file_name: Name of the file, presumably an individual reference package file (e.g. MSA, HMM, tree)
        :param text: Either a list (created by file_handler.readlines()) or a string to be written to file_name
        :return: None
        """
        try:
            file_h = open(file_name, 'w')
            if type(text) is list:
                text = ''.join(text)
            file_h.write(text)
            file_h.close()
        except IOError:
            logging.warning("Unable to write reference package component to '{}' for '{}'.\n".format(file_name,
                                                                                                     self.prefix))
        return

    # def copy_refpkg_file_to_dest(self, destination_dir, prefix=None) -> None:
    #     if prefix:
    #         intermediate_prefix = destination_dir + os.sep + prefix
    #     else:
    #         intermediate_prefix = destination_dir + os.sep + self.prefix
    #     copy(self.msa, intermediate_prefix + ".fa")
    #     copy(self.profile, intermediate_prefix + ".hmm")
    #     copy(self.search_profile, intermediate_prefix + "_search.hmm")
    #     copy(self.tree, intermediate_prefix + "_tree.txt")
    #     copy(self.model_info, intermediate_prefix + "_bestModel.txt")
    #     if os.path.isfile(self.boot_tree):
    #         copy(self.boot_tree, intermediate_prefix + "_bipartitions.txt")
    #         os.remove(self.boot_tree)
    #     return

    def change_file_paths(self, new_dir: str, move=False) -> None:
        """
        Used for changing all of the reference package file paths to another directory.
        These include: f__msa, f__profile, f__search_profile, f__tree, f__boot_tree, f__model_info

        :param new_dir: Path to another directory where the file does (or should) exist
        :param move: Boolean indicating whether the file should be moved
        :return: None
        """
        for a, v in self.__iter__():
            if a.startswith('f__'):
                new_path = new_dir + os.path.basename(v)
                if move:
                    copy(v, new_path)
                    os.remove(v)
                self.__dict__[a] = new_path
        return

    def update_file_names(self):
        for a, v in self.__iter__():
            if a.startswith('f__'):
                path, name = os.path.split(v)
                self.__dict__[a] = os.path.join(path, self.prefix + re.sub(self.prefix, '', name))
        return

    def disband(self, output_dir: str):
        """
        From a ReferencePackage's JSON file, the individual file components (e.g. profile HMM, MSA, phylogeny) are
        written to their separate files in a new directory, a sub-directory of where the JSON is located.
        The directory name follows the format: self.prefix + '_' self.refpkg_code + '_' + self.date

        :param output_dir: Output directory to write the individual reference package files
        :return:
        """
        try:
            if not os.path.isdir(output_dir):
                os.mkdir(output_dir)
        except (FileNotFoundError, IOError):
            logging.error("Unable to create directory '{0}' for ReferencePackage '{1}'.\n"
                          "TreeSAPP will only create a single directory at a time.\n".format(output_dir, self.prefix))
            sys.exit(3)

        output_prefix = os.path.join(output_dir, '_'.join([self.prefix, self.refpkg_code, self.date])) + os.sep

        if not os.path.isdir(output_prefix):
            os.mkdir(output_prefix)

        # Add the output directory prefix to each file name
        self.change_file_paths(output_prefix)

        self.write_refpkg_component(self.f__msa, self.msa)
        self.write_refpkg_component(self.f__profile, self.profile)
        self.write_refpkg_component(self.f__search_profile, self.search_profile)
        self.write_refpkg_component(self.f__model_info, self.model_info)
        self.write_refpkg_component(self.f__tree, self.tree)
        if len(self.boot_tree) > 0:
            self.write_refpkg_component(output_prefix + self.f__boot_tree, self.boot_tree)

        return

    def slurp(self) -> None:
        """
        Reads the reference package's JSON-formatted file and stores the elements in their respective variables

        :return: None
        """
        if len(self.f__json) == 0:
            logging.error("ReferencePackage.f__json was not set.\n")
            sys.exit(11)

        if not os.path.isfile(self.f__json):
            logging.error("ReferencePackage JSON file '{}' doesn't exist.\n".format(self.f__json))
            sys.exit(7)

        try:
            refpkg_handler = open(self.f__json, 'r')
        except IOError:
            logging.error("Unable to open reference package JSON file '{}' for reading.\n".format(self.f__json))
            sys.exit(3)
        refpkg_data = json.load(refpkg_handler)
        refpkg_handler.close()

        for a, v in refpkg_data.items():
            self.__dict__[a] = v

        self.load_taxonomic_hierarchy()

        return

    def validate(self, check_files=False):
        """
        Function that ensures the number of sequences is equal across all files and that in the ref_build_parameters.tsv

        :return: Boolean
        """
        # Check to ensure all files exist
        if check_files:
            for ref_file in self.core_ref_files:
                if not os.path.isfile(ref_file):
                    logging.error("File '" + ref_file + "' does not exist for reference package: " + self.prefix + "\n")
                    sys.exit(17)
        # TODO: Compare the number of sequences in the multiple sequence alignment
        # self.num_seqs = len(get_headers(self.f__msa))
        # TODO: Compare the number of sequences in the Hidden-Markov model
        # TODO: Compare the number of sequences in the Tree files
        # TODO: Compare the number of sequences in the tax_ids file
        # num_taxa = len(self.tax_ids_file_to_leaves())
        # if num_taxa != num_ref_seqs != self.num_seqs:
        #     logging.error("Number of reference sequences for reference package '" +
        #                   self.prefix + "' is inconsistent!\n")
        #     sys.exit(3)
        if not self.sub_model:
            logging.error("Unable to find the substitution model used for ReferencePackage '{}'.\n".format(self.prefix))
            sys.exit(33)
        return True

    def generate_tree_leaf_references_from_refpkg(self) -> list:
        """
        From the dictionary containing lineage and organism information of reference sequences (self.lineage_ids)
        this function creates a list of TreeLeafReference instances for every reference sequence

        :return: List of TreeLeafReference instances
        """
        ref_leaf_nodes = []
        for treesapp_id in sorted(self.lineage_ids, key=int):
            seq_name, lineage = self.lineage_ids[treesapp_id].split("\t")
            ref = TreeLeafReference(treesapp_id, seq_name)
            ref.lineage = lineage
            ref_leaf_nodes.append(ref)
        return ref_leaf_nodes

    def load_taxonomic_hierarchy(self) -> None:
        """
        Loads the ReferencePackage's taxonomic lineages into it's TaxonomicHierarchy (self.taxa_trie) using
        TaxonomicHierarchy.feed_leaf_nodes.

        :return: None
        """
        ref_leaf_nodes = self.generate_tree_leaf_references_from_refpkg()
        self.taxa_trie.feed_leaf_nodes(ref_leaf_nodes)
        self.taxa_trie.validate_rank_prefixes()
        return

    # def tax_ids_file_to_leaves(self):
    #     tree_leaves = list()
    #     unknown = 0
    #     try:
    #         tax_ids_handler = open(self.lineage_ids, 'r', encoding='utf-8')
    #     except IOError:
    #         logging.error("Unable to open " + self.lineage_ids + "\n")
    #         sys.exit(5)
    #
    #     for line in tax_ids_handler:
    #         line = line.strip()
    #
    #         try:
    #             number, seq_name, lineage = line.split("\t")
    #         except (ValueError, IndexError):
    #             logging.error("Unexpected number of fields in " + self.lineage_ids +
    #                           ".\nInvoked .split(\'\\t\') on line " + str(line) + "\n")
    #             sys.exit(5)
    #         leaf = TreeLeafReference(number, seq_name)
    #         if lineage:
    #             leaf.lineage = lineage
    #             leaf.complete = True
    #         else:
    #             unknown += 1
    #         tree_leaves.append(leaf)
    #
    #     if len(tree_leaves) == unknown:
    #         logging.error("Lineage information was not properly loaded for " + self.lineage_ids + "\n")
    #         sys.exit(5)
    #
    #     tax_ids_handler.close()
    #     return tree_leaves

    def hmm_length(self):
        self.profile_length = get_hmm_length(self.f__profile)

    def remove_taxon_from_lineage_ids(self, target_taxon) -> list:
        """
        Removes all sequences/leaves from the reference package that match the target taxon. Leaves with that have a
        taxonomic resolution lower than the target are also removed as their taxonomic provenance is uncertain.

        :param target_taxon: A '; '-separated taxonomic lineage for which all matches and descendents are removed
        :return: A list of LeafNode objects that don't match the target_taxon and have sufficient lineage depth
        """
        off_target_ref_leaves = list()
        depth = len(target_taxon.split("; "))
        n_match = 0
        n_shallow = 0
        n_unclassified = 0
        for ref_leaf in self.tax_ids_file_to_leaves():
            sc_lineage = ref_leaf.lineage.split("; ")
            if len(sc_lineage) < depth:
                n_shallow += 1
                continue
            if target_taxon == '; '.join(sc_lineage[:depth + 1]):
                n_match += 1
                continue
            if re.search("unclassified|environmental sample", ref_leaf.lineage, re.IGNORECASE):
                i = 0
                while i <= depth:
                    if re.search("unclassified|environmental sample", sc_lineage[i], re.IGNORECASE):
                        i -= 1
                        break
                    i += 1
                if i < depth:
                    n_unclassified += 1
                    continue
            off_target_ref_leaves.append(ref_leaf)

        logging.debug("Reference sequence filtering stats for " + target_taxon + "\n" +
                      "\n".join(["Match taxon\t" + str(n_match),
                                 "Unclassified\t" + str(n_unclassified),
                                 "Too shallow\t" + str(n_shallow),
                                 "Remaining\t" + str(len(off_target_ref_leaves))]) + "\n")
        return off_target_ref_leaves

    def clean_up_raxmlng_outputs(self, phylogeny_dir: str, fasta_replace_dict: dict) -> None:
        output_prefix = phylogeny_dir + self.prefix
        # Gather the best tree file
        try:
            raw_newick_tree = glob(output_prefix + ".*.bestTree")[0]
        except IndexError:
            logging.error("Unable to find " + output_prefix + ".*.bestTree generated by either FastTree or RAxML-NG.\n")
            sys.exit(17)
        # Gather the best model file
        try:
            model_info = glob(output_prefix + ".*.bestModel")[0]
        except IndexError:
            logging.error("Unable to find " + output_prefix + ".*.bestModel generated RAxML-NG.\n")
            sys.exit(17)

        copy(model_info, self.f__model_info)
        swap_tree_names(raw_newick_tree, self.f__tree)
        bootstrap_tree = output_prefix + ".raxml.support"
        if os.path.isfile(bootstrap_tree):
            annotate_partition_tree(self.prefix, fasta_replace_dict, bootstrap_tree)
            swap_tree_names(bootstrap_tree, self.f__boot_tree)

        intermediates = [raw_newick_tree, model_info,
                         output_prefix + ".raxml.log",
                         output_prefix + ".raxml.rba",
                         output_prefix + ".raxml.reduced.phy",
                         output_prefix + ".raxml.startTree"]
        for f in intermediates:
            try:
                os.remove(f)
            except OSError:
                logging.debug("Unable to remove %s as it doesn't exist.\n" % f)

        return

    def exclude_clade_from_ref_files(self, treesapp_refpkg_dir: str, molecule: str,
                                     original_storage_dir: str, target_clade: str, executables: dict,
                                     fresh=False) -> str:
        """
        Removes all reference sequences/leaf nodes from a reference package that are descendents of a target clade.
        All reference package files are regenerated without these sequences and written to the original reference
        package directory.

        Original reference package files (so, still containing descendent sequences) are saved to a specified location.

        :param treesapp_refpkg_dir: Path to the treesapp/data (or reference package) directory for RAxML-NG to write
         output files from its `--evaluate` routine.
        :param molecule: Any of 'prot', 'rrna' or 'dna'
        :param original_storage_dir: Path to directory where original reference package files need to be moved to
        :param target_clade: Taxonomic lineage of the clade that is being excluded from the reference package.
        :param executables: Dictionary of paths to dependency executables indexed by their names. Must include:
         'hmmbuild', 'FastTree' and 'raxml-ng'.
        :param fresh: Boolean indicating whether the reference package's tree should be built from scratch (True) or
         if the clades that are descendents of 'target_clade' should just be pruned (False) by ETE3
        :return: Path to all the original reference package files that have been copied to prevent overwriting
        """
        intermediate_prefix = original_storage_dir + os.sep + "ORIGINAL"
        self.copy_refpkg_file_to_dest(original_storage_dir, "ORIGINAL")

        # tax_ids
        off_target_ref_leaves = self.remove_taxon_from_lineage_ids(target_clade)
        with open(self.lineage_ids, 'w') as tax_ids_handle:
            tax_ids_string = ""
            for ref_leaf in off_target_ref_leaves:
                tax_ids_string += "\t".join([ref_leaf.number, ref_leaf.description, ref_leaf.lineage]) + "\n"
            tax_ids_handle.write(tax_ids_string)

        # fasta
        ref_fasta_dict = read_fasta_to_dict(self.msa)
        off_target_ref_headers = [ref_leaf.number + '_' + self.prefix for ref_leaf in off_target_ref_leaves]
        if len(off_target_ref_headers) == 0:
            logging.error("No reference sequences were retained for building testing " + target_clade + "\n")
            sys.exit(19)
        split_files = write_new_fasta(ref_fasta_dict, self.msa, headers=off_target_ref_headers)
        if len(split_files) > 1:
            logging.error("Only one FASTA file should have been written.\n")
            sys.exit(21)

        # HMM profile
        hmm_build_command = [executables["hmmbuild"], self.profile, self.msa]
        launch_write_command(hmm_build_command)

        # Trees
        if fresh:
            tree_build_cmd = [executables["FastTree"]]
            if molecule == "rrna" or molecule == "dna":
                tree_build_cmd += ["-nt", "-gtr"]
            else:
                tree_build_cmd += ["-lg", "-wag"]
            tree_build_cmd += ["-out", self.tree]
            tree_build_cmd.append(self.msa)
            logging.info("Building Approximately-Maximum-Likelihood tree with FastTree... ")
            stdout, returncode = launch_write_command(tree_build_cmd, True)
            with open(original_storage_dir + os.sep + "FastTree_info." + self.prefix, 'w') as fast_info:
                fast_info.write(stdout + "\n")
            logging.info("done.\n")
        else:
            ref_tree = Tree(self.tree)
            ref_tree.prune(off_target_ref_headers)
            logging.debug("\t" + str(len(ref_tree.get_leaves())) + " leaves in pruned tree.\n")
            ref_tree.write(outfile=self.tree, format=5)
        # Model parameters
        model_parameters(executables["raxml-ng"], self.msa, self.tree,
                         treesapp_refpkg_dir + os.sep + "tree_data" + os.sep + self.prefix, self.sub_model)
        self.clean_up_raxmlng_outputs(treesapp_refpkg_dir + os.sep + "tree_data" + os.sep, {})
        return intermediate_prefix

    # def restore_reference_package(self, prefix: str, output_dir: str) -> None:
    #     """
    #
    #     :param prefix: Prefix (path and basename) of the stored temporary files
    #     :param output_dir: Path to the output directory for any temporary files that should be stored
    #     :return: None
    #     """
    #     # The edited tax_ids file with clade excluded is required for performance analysis
    #     # Copy the edited, clade-excluded tax_ids file to the output directory
    #     copy(self.lineage_ids, output_dir)
    #
    #     # Move the original reference package files back to the proper directories
    #     copy(prefix + "_tree.txt", self.tree)
    #     copy(prefix + "_bestModel.txt", self.model_info)
    #     if os.path.isfile(prefix + "_bipartitions.txt"):
    #         copy(prefix + "_bipartitions.txt", self.boot_tree)
    #     copy(prefix + "_tax_ids.txt", self.lineage_ids)
    #     copy(prefix + ".fa", self.msa)
    #     copy(prefix + ".hmm", self.profile)
    #     copy(prefix + "_search.hmm", self.search_profile)
    #
    #     return

    def load_pfit_params(self, build_param_line):
        build_param_fields = build_param_line.split("\t")
        if build_param_fields[8]:
            self.pfit = [float(x) for x in build_param_fields[8].split(',')]
        return

    def get_info(self):
        return "\n\t".join(["MarkerBuild instance of %s (%s):" % (self.prefix, self.refpkg_code),
                            "Molecule type:                                      " + self.molecule,
                            "Substitution model used for phylogenetic inference: " + self.sub_model,
                            "Number of reference sequences (leaf nodes):         " + str(self.num_seqs),
                            "Software used to infer phylogeny:                   " + self.tree_tool,
                            "Date of last update:                                " + self.update,
                            "Description:                                        '%s'" % self.description]) + "\n"
