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


class ReferencePackage:
    def __init__(self, refpkg_name=""):
        self.prefix = refpkg_name
        self.refpkg_code = ""  # AKA denominator

        # These are files (with '_f' suffix) and their respective data (read with file.readlines())
        self.json_path = ""  # Path to the JSON reference package file
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
        self.lineage_ids = dict()  # Reference sequence lineage map (tax_ids)

        # These are metadata values
        self.ts_version = ""
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

    def to_dict(self) -> dict:
        # Remove all of the non-primitives from the dictionary
        refpkg_dict = self.__dict__
        non_primitives = ["core_ref_files", "taxa_trie"]
        for bad_attr in non_primitives:
            refpkg_dict.pop(bad_attr)
        return refpkg_dict

    def band(self):
        """
        Writes a JSON file containing all of the reference package files and metadata

        :return:
        """
        if len(self.json_path) == 0:
            logging.error("ReferencePackage JSON file not set. ReferencePackage band cannot be completed.\n")
            sys.exit(11)

        # Read the MSA, profile HMMs, and phylogenies
        with open(self.f__msa) as fh:
            self.msa = fh.readlines()
        with open(self.f__profile) as fh:
            self.profile = fh.readlines()
        with open(self.f__search_profile) as fh:
            self.search_profile = fh.readlines()
        with open(self.f__tree) as fh:
            self.tree = fh.readlines()
        with open(self.f__model_info) as fh:
            self.model_info = fh.readlines()
        if os.path.isfile(self.f__boot_tree):  # This file isn't guaranteed to exist in all cases
            with open(self.f__boot_tree) as fh:
                self.boot_tree = fh.readlines()

        try:
            refpkg_handler = open(self.json_path, 'w')
        except IOError:
            logging.error("Unable to open reference package JSON file '{}' for writing.\n".format(self.json_path))
            sys.exit(11)

        json.dump(obj=self.to_dict(), fp=refpkg_handler)

        refpkg_handler.close()

        return

    def write_refpkg_component(self, file_name, text):
        try:
            file_h = open(file_name, 'w')
            file_h.write(text)
            file_h.close()
        except IOError:
            logging.warning("Unable to write reference package component to '{}' for '{}'.\n".format(file_name,
                                                                                                     self.prefix))
        return

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

    def disband(self, output_dir):
        """
        From a ReferencePackage's JSON file, the individual file components (e.g. profile HMM, MSA, phylogeny) are
        written to their separate files in a new directory, a sub-directory of where the JSON is located.
        The directory name follows the format: self.prefix + '_' self.refpkg_code + '_' + self.date


        :return:
        """
        output_prefix = os.path.join(output_dir, '_'.join([self.prefix, self.refpkg_code, self.date])) + os.sep

        if not os.path.isdir(output_prefix):
            os.mkdir(output_prefix)

        # Add the output directory prefix to each file name
        self.change_file_paths(output_prefix)
        print(self.f__msa)
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
        try:
            refpkg_handler = open(self.json_path, 'r')
        except IOError:
            logging.error("Unable to open reference package JSON file '{}' for reading.\n".format(self.json_path))
            sys.exit(11)
        refpkg_data = json.load(refpkg_handler)
        refpkg_handler.close()

        for a, v in refpkg_data.items():
            self.__dict__[a] = v

        return

    def validate(self, num_ref_seqs=0):
        """
        Function that ensures the number of sequences is equal across all files and that in the ref_build_parameters.tsv

        :return: Boolean
        """
        # Check to ensure all files exist
        for ref_file in self.core_ref_files:
            if not os.path.isfile(ref_file):
                logging.error("File '" + ref_file + "' does not exist for reference package: " + self.prefix + "\n")
                sys.exit(17)
        # TODO: Compare the number of sequences in the multiple sequence alignment
        self.num_seqs = len(get_headers(self.msa))
        if not num_ref_seqs:
            num_ref_seqs = self.num_seqs
        # TODO: Compare the number of sequences in the Hidden-Markov model
        # TODO: Compare the number of sequences in the Tree files
        # TODO: Compare the number of sequences in the tax_ids file
        num_taxa = len(self.tax_ids_file_to_leaves())
        if num_taxa != num_ref_seqs != self.num_seqs:
            logging.error("Number of reference sequences for reference package '" +
                          self.prefix + "' is inconsistent!\n")
            sys.exit(3)
        return True

    def tax_ids_file_to_leaves(self):
        tree_leaves = list()
        unknown = 0
        try:
            tax_ids_handler = open(self.lineage_ids, 'r', encoding='utf-8')
        except IOError:
            logging.error("Unable to open " + self.lineage_ids + "\n")
            sys.exit(5)

        for line in tax_ids_handler:
            line = line.strip()

            try:
                number, seq_name, lineage = line.split("\t")
            except (ValueError, IndexError):
                logging.error("Unexpected number of fields in " + self.lineage_ids +
                              ".\nInvoked .split(\'\\t\') on line " + str(line) + "\n")
                sys.exit(5)
            leaf = TreeLeafReference(number, seq_name)
            if lineage:
                leaf.lineage = lineage
                leaf.complete = True
            else:
                unknown += 1
            tree_leaves.append(leaf)

        if len(tree_leaves) == unknown:
            logging.error("Lineage information was not properly loaded for " + self.lineage_ids + "\n")
            sys.exit(5)

        tax_ids_handler.close()
        return tree_leaves

    def gather_package_files(self, pkg_path: str, molecule="prot", layout=None) -> None:
        """
        Populates a ReferencePackage instances fields with files based on 'pkg_format' where hierarchical indicates
        files are sorted into 'alignment_data', 'hmm_data' and 'tree_data' directories and flat indicates they are all
        in the same directory.

        :param pkg_path: Path to the reference package
        :param molecule: A string indicating the molecule type of the reference package. If 'rRNA' profile is CM.
        :param layout: An optional string indicating what layout to use (flat | hierarchical)
        :return: None
        """
        if not self.prefix:
            logging.error("ReferencePackage.prefix not set - unable to gather package files.\n")
            sys.exit(3)
        if molecule == "rRNA":
            profile_ext = ".cm"
        else:
            profile_ext = ".hmm"
        if pkg_path[-1] != os.sep:
            pkg_path += os.sep

        flat = {"msa": pkg_path + self.prefix + ".fa",
                "tree": pkg_path + self.prefix + "_tree.txt",
                "profile": pkg_path + self.prefix + profile_ext,
                "search_profile": pkg_path + self.prefix + "_search" + profile_ext,
                "taxid": pkg_path + "tax_ids_" + self.prefix + ".txt",
                "bestModel": pkg_path + self.prefix + "_bestModel.txt"}
        hierarchical = {"msa": pkg_path + "alignment_data" + os.sep + self.prefix + ".fa",
                        "profile": pkg_path + "hmm_data" + os.sep + self.prefix + profile_ext,
                        "search_profile": pkg_path + "hmm_data" + os.sep + self.prefix + "_search" + profile_ext,
                        "taxid": pkg_path + "tree_data" + os.sep + "tax_ids_" + self.prefix + ".txt",
                        "tree": pkg_path + "tree_data" + os.sep + self.prefix + "_tree.txt",
                        "bestModel": pkg_path + "tree_data" + os.sep + self.prefix + "_bestModel.txt"}

        if layout == "flat":
            layout = flat
        elif layout == "hierarchical":
            layout = hierarchical
        else:
            # Exhaustively test whether all predicted files exist for each layout
            layout = None
            for option in [flat, hierarchical]:
                acc = 0
                for f_type in option:
                    if os.path.exists(option[f_type]):
                        acc += 1
                    else:
                        break
                if acc == len(option):
                    layout = option
                    break

        if layout:
            self.f__msa = layout["msa"]
            self.f__tree = layout["tree"]
            self.f__profile = layout["profile"]
            self.f__search_profile = layout["search_profile"]
            self.f__model_info = layout["bestModel"]
            self.f__boot_tree = os.path.dirname(layout["tree"]) + os.sep + self.prefix + "_bipartitions.txt"
        else:
            logging.error("Unable to gather reference package files for " + self.prefix + " from '" + pkg_path + "'\n")
            raise AssertionError()

        return

    def hmm_length(self):
        self.profile_length = get_hmm_length(self.profile)

    def copy_refpkg_file_to_dest(self, destination_dir, prefix=None) -> None:
        if prefix:
            intermediate_prefix = destination_dir + os.sep + prefix
        else:
            intermediate_prefix = destination_dir + os.sep + self.prefix
        copy(self.msa, intermediate_prefix + ".fa")
        copy(self.profile, intermediate_prefix + ".hmm")
        copy(self.search_profile, intermediate_prefix + "_search.hmm")
        copy(self.tree, intermediate_prefix + "_tree.txt")
        copy(self.model_info, intermediate_prefix + "_bestModel.txt")
        if os.path.isfile(self.boot_tree):
            copy(self.boot_tree, intermediate_prefix + "_bipartitions.txt")
            os.remove(self.boot_tree)
        copy(self.lineage_ids, intermediate_prefix + "_tax_ids.txt")
        return

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

        copy(model_info, self.model_info)
        swap_tree_names(raw_newick_tree, self.tree)
        bootstrap_tree = output_prefix + ".raxml.support"
        if os.path.isfile(bootstrap_tree):
            annotate_partition_tree(self.prefix, fasta_replace_dict, bootstrap_tree)
            swap_tree_names(bootstrap_tree, self.boot_tree)

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

    def restore_reference_package(self, prefix: str, output_dir: str) -> None:
        """

        :param prefix: Prefix (path and basename) of the stored temporary files
        :param output_dir: Path to the output directory for any temporary files that should be stored
        :return: None
        """
        # The edited tax_ids file with clade excluded is required for performance analysis
        # Copy the edited, clade-excluded tax_ids file to the output directory
        copy(self.lineage_ids, output_dir)

        # Move the original reference package files back to the proper directories
        copy(prefix + "_tree.txt", self.tree)
        copy(prefix + "_bestModel.txt", self.model_info)
        if os.path.isfile(prefix + "_bipartitions.txt"):
            copy(prefix + "_bipartitions.txt", self.boot_tree)
        copy(prefix + "_tax_ids.txt", self.lineage_ids)
        copy(prefix + ".fa", self.msa)
        copy(prefix + ".hmm", self.profile)
        copy(prefix + "_search.hmm", self.search_profile)

        return


class MarkerBuild:
    def __init__(self):
        self.cog = ""
        self.denominator = ""
        self.molecule = ""
        self.model = ""
        self.lowest_confident_rank = ""
        self.update = ""
        self.kind = ""
        self.tree_tool = ""
        self.description = ""
        self.pid = 1.0
        self.num_reps = 0
        self.pfit = []

    def load_build_params(self, build_param_line, n_fields):
        build_param_fields = build_param_line.split('\t')
        if len(build_param_fields) != n_fields:
            logging.error("Incorrect number of values (" + str(len(build_param_fields)) +
                          ") in ref_build_parameters.tsv. Line:\n" + build_param_line)
            sys.exit(17)

        self.cog = build_param_fields[0]
        self.denominator = build_param_fields[1]
        self.molecule = build_param_fields[2]
        self.model = build_param_fields[3]
        self.kind = build_param_fields[4]
        self.pid = float(build_param_fields[5])
        self.num_reps = int(build_param_fields[6])
        self.tree_tool = build_param_fields[7]
        self.lowest_confident_rank = build_param_fields[9]
        self.update = build_param_fields[10]
        self.description = build_param_fields[-1].strip()

    def load_pfit_params(self, build_param_line):
        build_param_fields = build_param_line.split("\t")
        if build_param_fields[8]:
            self.pfit = [float(x) for x in build_param_fields[8].split(',')]
        return

    def check_rank(self):
        taxonomies = ["NA", "Kingdoms", "Phyla", "Classes", "Orders", "Families", "Genera", "Species"]

        if self.lowest_confident_rank not in list(taxonomies):
            logging.error("Unable to find '" + self.lowest_confident_rank + "' in taxonomic map!\n")
            sys.exit(17)

        return

    def get_info(self):
        return "\n\t".join(["MarkerBuild instance of %s (%s):" % (self.cog, self.denominator),
                            "Molecule type:                                      " + self.molecule,
                            "Substitution model used for phylogenetic inference: " + self.model,
                            "Number of reference sequences (leaf nodes):         " + str(self.num_reps),
                            "Software used to infer phylogeny:                   " + self.tree_tool,
                            "Date of last update:                                " + self.update,
                            "Description:                                        '%s'" % self.description]) + "\n"

    def attributes_to_dict(self) -> dict:
        metadata = {"RefPkg":
                        {"name": self.cog,
                         "code": self.denominator,
                         "molecule": self.molecule,
                         "sequences": self.num_reps,
                         "model": self.model,
                         "lowest_confident_rank": self.lowest_confident_rank,
                         "kind": self.kind,
                         "tree-tool": self.tree_tool,
                         "regression-parameters": self.pfit,
                         "cluster-similarity": self.pid,
                         "description": self.description}
                    }
        return metadata

    def dict_to_attributes(self, build_params: dict) -> None:
        refpkg_params = build_params["RefPkg"]
        self.cog = refpkg_params["name"]
        self.denominator = refpkg_params["code"]
        self.description = refpkg_params["description"]
        self.kind = refpkg_params["kind"]
        self.lowest_confident_rank = refpkg_params["lowest_confident_rank"]
        self.model = refpkg_params["model"]
        self.molecule = refpkg_params["molecule"]
        self.num_reps = refpkg_params["sequences"]
        self.pfit = refpkg_params["regression-parameters"]
        self.pid = refpkg_params["cluster-similarity"]
        self.tree_tool = refpkg_params["tree-tool"]

        return
