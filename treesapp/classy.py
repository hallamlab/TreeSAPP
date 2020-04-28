
__author__ = 'Connor Morgan-Lang'

import sys
import os
import re
import logging
import time
from shutil import rmtree, copy
from copy import deepcopy
from multiprocessing import Process
from glob import glob
from json import loads, dumps
from collections import namedtuple
from numpy import var

from ete3 import Tree

from .fasta import format_read_fasta, write_new_fasta, get_header_format, FASTA, get_headers, load_fasta_header_regexes, read_fasta_to_dict
from .utilities import median, which, is_exe, return_sequence_info_groups, write_dict_to_table,\
    load_pickle, swap_tree_names, fish_refpkg_from_build_params
from .entish import create_tree_info_hash, subtrees_to_dictionary, annotate_partition_tree
from .lca_calculations import determine_offset, optimal_taxonomic_assignment
from .external_command_interface import launch_write_command
from . import entrez_utils
from .wrapper import model_parameters, CommandLineFarmer
from .taxonomic_hierarchy import TaxonomicHierarchy

import _tree_parser


class ModuleFunction:
    def __init__(self, name, order, func=None):
        self.order = order
        self.name = name
        self.function = func
        self.run = True

    def get_info(self):
        info_string = "Information for '" + self.name + "':\n"
        info_string += "\tOrder: " + str(self.order) + "\n"
        info_string += "\tRun: " + str(self.run) + "\n"
        return info_string


class ReferencePackage:
    def __init__(self, refpkg_name=""):
        self.prefix = refpkg_name
        self.refpkg_code = ""  # AKA denominator
        self.msa = ""  # Reference MSA FASTA
        self.profile = ""  # HMM file
        self.search_profile = ""  # HMM file
        self.tree = ""  # Reference tree
        self.boot_tree = ""  # Reference tree with support values
        self.lineage_ids = ""  # Reference sequence lineage map (tax_ids)
        self.taxa_trie = TaxonomicHierarchy()
        self.sub_model = ""  # EPA-NG compatible substitution model
        self.model_info = ""  # RAxML-NG --evaluate model file
        self.core_ref_files = list()
        self.num_seqs = 0
        self.profile_length = 0  # LENG of the HMM profile

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
                fields = line.split("\t")
            except ValueError:
                logging.error('ValueError: .split(\'\\t\') on ' + str(line) +
                              " generated " + str(len(line.split("\t"))) + " fields.\n")
                sys.exit(5)
            try:
                number, seq_name, lineage = fields
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
            self.msa = layout["msa"]
            self.tree = layout["tree"]
            self.profile = layout["profile"]
            self.search_profile = layout["search_profile"]
            self.lineage_ids = layout["taxid"]
            self.model_info = layout["bestModel"]
            self.boot_tree = os.path.dirname(layout["tree"]) + os.sep + self.prefix + "_bipartitions.txt"
        else:
            logging.error("Unable to gather reference package files for " + self.prefix + " from '" + pkg_path + "'\n")
            raise AssertionError()

        self.core_ref_files += [self.msa, self.profile, self.tree, self.lineage_ids, self.model_info]

        return

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


class ItolJplace:
    """
    A class to hold all data relevant to a jplace file to be viewed in iTOL
    """
    fields = list()

    def __init__(self):
        self.contig_name = ""  # Sequence name (from FASTA header)
        self.name = ""  # Code name of the tree it mapped to (e.g. mcrA)
        self.abundance = None  # Either the number of occurences, or the FPKM of that sequence
        self.node_map = dict()  # A dictionary mapping internal nodes (Jplace) to all leaf nodes
        self.seq_len = 0
        ##
        # Taxonomic information:
        ##
        self.lineage_list = list()  # List containing each child's lineage
        self.wtd = 0
        self.lct = ""  # The LCA taxonomy derived from lineage_list
        self.recommended_lineage = ""
        ##
        # Information derived from Jplace pqueries:
        ##
        self.placements = list()
        self.lwr = 0  # Likelihood weight ratio of an individual placement
        self.likelihood = 0
        self.avg_evo_dist = 0.0
        self.distances = ""
        self.classified = True
        self.inode = ""
        self.tree = ""  # NEWICK tree
        self.metadata = ""
        self.version = ""  # Jplace version

    def summarize(self):
        """
        Prints a summary of the ItolJplace object (equivalent to a single marker) to stderr
        Summary include the number of marks found, the tree used, and the tree-placement of each sequence identified
        Written solely for testing purposes

        :return:
        """
        summary_string = ""
        summary_string += "\nInformation for query sequence '" + str(self.contig_name) + "'\n"
        summary_string += str(len(self.placements)) + " sequence(s) grafted onto the " + self.name + " tree.\n"
        # summary_string += "Reference tree:\n")
        # summary_string += self.tree + "\n")
        summary_string += "JPlace fields:\n\t" + str(self.fields) + "\n"
        summary_string += "Placement information:\n"
        if not self.placements:
            summary_string += "\tNone.\n"
        elif self.placements[0] == '{}':
            summary_string += "\tNone.\n"
        else:
            if self.likelihood and self.lwr and self.inode:
                summary_string += "\tInternal node\t" + str(self.inode) + "\n"
                summary_string += "\tLikelihood\t" + str(self.likelihood) + "\n"
                summary_string += "\tL.W.R\t\t" + str(self.lwr) + "\n"
            else:
                for pquery in self.placements:
                    placement = loads(str(pquery), encoding="utf-8")
                    for k, v in placement.items():
                        if k == 'p':
                            summary_string += '\t' + str(v) + "\n"
        summary_string += "Non-redundant lineages of child nodes:\n"
        if len(self.lineage_list) > 0:
            for lineage in sorted(set(self.lineage_list)):
                summary_string += '\t' + str(lineage) + "\n"
        else:
            summary_string += "\tNone.\n"
        summary_string += "Lowest common taxonomy:\n"
        if self.lct:
            summary_string += "\t" + str(self.lct) + "\n"
        else:
            summary_string += "\tNone.\n"
        if self.abundance:
            summary_string += "Abundance:\n\t" + str(self.abundance) + "\n"
        if self.distances:
            summary_string += "Distal, pendant and tip distances:\n\t" + self.distances + "\n"
        summary_string += "\n"
        return summary_string

    def list_placements(self):
        """
        Returns a list of all the nodes contained in placements
        :return:
        """
        nodes = list()
        for d_place in self.placements:
            if isinstance(d_place, str):
                for k, v in loads(d_place).items():
                    if k == 'p':
                        for pquery in v:
                            nodes.append(str(pquery[0]))
            else:
                logging.error("Unable to handle type " + type(d_place) + "\n")
                sys.exit(17)
        return nodes

    def correct_decoding(self) -> None:
        """
        Since the JSON decoding is unable to decode recursively, this needs to be fixed for each placement
        Formatting and string conversion are also performed here

        :return: None
        """
        new_placement_collection = []  # a list of dictionary-like strings
        placement_string = ""  # e.g. {"p":[[226, -31067.028237, 0.999987, 0.012003, 2e-06]], "n":["query"]}
        for d_place in self.placements:
            if not isinstance(d_place, str):
                dict_strings = list()  # e.g. "n":["query"]
                for k, v in d_place.items():
                    dict_strings.append(dumps(k) + ':' + dumps(v))
                    placement_string = ', '.join(dict_strings)
                new_placement_collection.append('{' + placement_string + '}')
            else:
                new_placement_collection.append(d_place)
        self.placements = new_placement_collection

        decoded_fields = list()
        for field in self.fields:
            if not re.match('".*"', field):
                decoded_fields.append(dumps(field))
            else:
                decoded_fields.append(field)
        self.fields = decoded_fields
        return

    def name_placed_sequence(self):
        for d_place in self.placements:
            for key, value in d_place.items():
                if key == 'n':
                    self.contig_name = value[0]
        return

    def get_field_position_from_jplace_fields(self, field_name) -> int:
        """
        Find the position of a specific field in self.fields

        :return: The integer position of a field name (string)
        """
        x = 0
        # Find the position of field_name in the placements from fields descriptor
        quoted_field = '"' + field_name + '"'
        for field in self.fields:
            if str(field) == quoted_field:
                break
            else:
                x += 1
        if x == len(self.fields):
            logging.warning("Unable to find '" + field_name + "' in the jplace \"field\" string!\n")
            return None
        return x

    def get_jplace_element(self, element_name) -> str:
        """
        Determines the element value (e.g. likelihood, edge_num) for a single placement.
        There may be multiple placements (or 'pquery's) in a single .jplace file, therefore, this function is usually looped over.

        :param element_name:
        :return:
        """
        position = self.get_field_position_from_jplace_fields(element_name)
        placement = loads(self.placements[0], encoding="utf-8")
        element_value = None
        for k, v in placement.items():
            if k == 'p':
                acc = 0
                while acc < len(v):
                    pquery_fields = v[acc]
                    element_value = pquery_fields[position]
                    acc += 1
        return element_value

    def filter_min_weight_threshold(self, threshold=0.1) -> None:
        """
        Sets the instance's *classified* attribute to False if the likelihood weight ratio (LWR)
        threshold is not met or exceeded.

        :param threshold: The threshold which all placements with LWRs less than this are removed
        :return: None
        """
        if len(self.placements) != 1:
            logging.error("Only one placement is expected here, but %d were found.\n%s" %
                          (len(self.placements), self.summarize()))
            sys.exit(5)

        if self.lwr < threshold:
            self.classified = False
        return

    def sum_rpkms_per_node(self, leaf_rpkm_sums):
        """
        Function that adds the RPKM value of a contig to the node it was placed.
        For contigs mapping to internal nodes: the proportional RPKM assigned is summed for all children.

        :param leaf_rpkm_sums: A dictionary mapping tree leaf numbers to abundances (RPKM sums)
        :return: dict()
        """
        for pquery in self.placements:
            placement = loads(pquery, encoding="utf-8")
            for k, v in placement.items():
                if k == 'p':
                    for locus in v:
                        jplace_node = locus[0]
                        tree_leaves = self.node_map[jplace_node]
                        try:
                            normalized_abundance = float(self.abundance/len(tree_leaves))
                        except TypeError:
                            logging.warning("Unable to find abundance for " + self.contig_name + "... setting to 0.\n")
                            normalized_abundance = 0.0
                        for tree_leaf in tree_leaves:
                            if tree_leaf not in leaf_rpkm_sums.keys():
                                leaf_rpkm_sums[tree_leaf] = 0.0
                            leaf_rpkm_sums[tree_leaf] += normalized_abundance
        return leaf_rpkm_sums

    def filter_max_weight_placement(self) -> None:
        """
        Removes all secondary placements of each pquery, leaving only the placement with the maximum like_weight_ratio

        :return: None
        """
        # Find the position of like_weight_ratio in the placements from fields descriptor
        x = self.get_field_position_from_jplace_fields("like_weight_ratio")
        if not x:
            return

        # Filter the placements
        new_placement_collection = list()
        placement_string = ""
        for pquery in self.placements:
            placement = loads(pquery, encoding="utf-8")
            if placement:
                dict_strings = list()
                max_lwr = 0
                if len(placement["p"]) > 1:
                    for k, v in placement.items():
                        if k == 'p':
                            acc = 0
                            tmp_placements = deepcopy(v)
                            while acc < len(tmp_placements):
                                candidate = tmp_placements[acc]
                                if float(candidate[x]) > max_lwr:
                                    v = [tmp_placements.pop(acc)]
                                    max_lwr = candidate[x]
                                else:
                                    acc += 1
                        dict_strings.append(dumps(k) + ':' + dumps(v))
                        placement_string = ', '.join(dict_strings)
                    # Add the filtered placements back to the object.placements
                    new_placement_collection.append('{' + placement_string + '}')
                else:
                    new_placement_collection.append(pquery)
        self.placements = new_placement_collection
        return

    def check_jplace(self, tree_index) -> None:
        """
        Currently validates a pquery's JPlace distal length, ensuring it is less than or equal to the edge length
        This is necessary to handle a case found in RAxML v8.2.12 (and possibly older versions) where the distal length
        of a placement is greater than the corresponding branch length in some rare cases.

        :return: None
        """
        distal_pos = self.get_field_position_from_jplace_fields("distal_length")
        edge_pos = self.get_field_position_from_jplace_fields("edge_num")
        for pquery in self.placements:
            placement = loads(pquery, encoding="utf-8")
            if placement:
                if len(placement["p"]) > 1:
                    for k, v in placement.items():
                        if k == 'p':
                            for edge_placement in v:
                                place_len = float(edge_placement[distal_pos])
                                edge = edge_placement[edge_pos]
                                tree_len = tree_index[str(edge)]
                                if place_len > tree_len:
                                    logging.debug("Distal length adjusted to fit JPlace " +
                                                  self.name + " tree for " + self.contig_name + ".\n")
                                    edge_placement[distal_pos] = tree_len
                else:
                    pass

        return

    def clear_object(self):
        self.placements.clear()
        self.fields.clear()
        self.node_map.clear()
        self.contig_name = ""
        self.name = ""
        self.tree = ""
        self.metadata = ""
        self.version = ""
        self.lineage_list = list()
        self.lct = ""
        self.abundance = None

    def lowest_confident_taxonomy(self, depth):
        """
        Truncates the initial taxonomic assignment to rank of depth.
        Uses self.lct - a string for the taxonomic lineage ('; ' separated)

        :param depth: The recommended depth to truncate the taxonomy
        :return: String representing 'confident' taxonomic assignment for the sequence
        """
        # Sequence likely isn't a FP but is highly divergent from reference set
        confident_assignment = "Root"
        if depth < 1:
            return confident_assignment

        confident_assignment = "; ".join(self.lct.split("; ")[:depth])

        return confident_assignment


class TreeProtein(ItolJplace):
    """
    A class for sequences that were properly mapped to its gene tree.
    While it mostly contains RAxML outputs, functions are used to make 'biological' sense out of these outputs.
    """
    def transfer(self, itol_jplace_object):
        self.placements = itol_jplace_object.placements
        self.tree = itol_jplace_object.tree
        self.fields = itol_jplace_object.fields
        self.version = itol_jplace_object.version
        self.metadata = itol_jplace_object.metadata

    def megan_lca(self):
        """
        Using the lineages of all leaves to which this sequence was mapped (n >= 1),
        A lowest common ancestor is found at the point which these lineages converge.
        This emulates the LCA algorithm employed by the MEtaGenome ANalyzer (MEGAN).
        :return:
        """
        # If there is only one child, return the joined string
        if len(self.lineage_list) == 1:
            return "; ".join(self.lineage_list[0])

        listed_lineages = [lineage.strip().split("; ") for lineage in self.lineage_list]
        max_depth = max([len(lineage) for lineage in listed_lineages])
        lca_set = set()
        lca_lineage_strings = list()
        i = 0
        while i < max_depth:
            contributors = 0
            for lineage in sorted(listed_lineages):
                try:
                    lca_set.add(lineage[i])
                    contributors += 1
                except IndexError:
                    pass

            if len(lca_set) == 1 and contributors == len(listed_lineages):
                lca_lineage_strings.append(list(lca_set)[0])
                i += 1
                lca_set.clear()
            else:
                i = max_depth

        return "; ".join(lca_lineage_strings)

    def children_lineage(self, leaves_taxa_map: dict):
        """
        From the jplace placement field ()

        :param leaves_taxa_map: Dictionary mapping tree leaf nodes to taxonomic lineages
        :return:
        """
        children = list()
        pquery = self.placements[0]
        placement = loads(pquery, encoding="utf-8")
        loci = placement['p']
        for locus in loci:
            jplace_node = locus[0]
            tree_leaves = self.node_map[jplace_node]
            for leaf_node in tree_leaves:
                try:
                    leaf_num = leaf_node.split('_')[0]
                except TypeError:
                    logging.error("Unexpected format of leaf node: '" + str(leaf_node) + "'.\n")
                    sys.exit(3)
                try:
                    ref_lineage = leaves_taxa_map[leaf_num]
                except KeyError:
                    logging.error("Unable to find '" + leaf_num + "' in leaf-lineage map.\n")
                    sys.exit(3)
                if ref_lineage:
                    children.append(ref_lineage)
                else:
                    logging.warning("No lineage information available for " + leaf_node + ".\n")

        return children


class TreeLeafReference:
    """
    Objects for each leaf in a tree
    """
    def __init__(self, number, description):
        self.number = number
        self.description = description
        self.lineage = ""
        self.accession = ""
        self.complete = False

    def summarize_tree_leaf(self):
        summary_string = "Leaf ID:\n\t" + str(self.number) + "\n" +\
                         "Description:\n\t" + str(self.description) + "\n"
        summary_string += "Accession:\n\t'" + self.accession + "'\n"
        if self.complete:
            summary_string += "Lineage:\n\t" + str(self.lineage) + "\n"
        return summary_string

    class MarkerInfo:
        """
        Class serves to store information pertaining to each COG in data/tree_data/cog_list.tsv
        """

        def __init__(self, marker, denominator, description):
            self.marker = marker
            self.denominator = denominator  # alphanumeric unique ID, R0016 for example
            self.marker_class = ""  # phylogenetic rRNA
            self.num_reference_seqs = 0
            self.description = description
            self.analysis_type = ""


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
            result = _tree_parser._build_subtrees_newick(next_task)
            subtrees = subtrees_to_dictionary(result, create_tree_info_hash())
            self.task_queue.task_done()
            self.result_queue.put(subtrees)
        return


def get_header_info(header_registry, code_name=''):
    """

    :param header_registry: A dictionary of Header instances, indexed by numerical treesapp_id
    :param code_name: [OPTIONAL] The code_name of the reference package (marker gene/domain/family/protein)
    :return: Dictionary where keys are numerical treesapp_ids and values are ReferenceSequence instances
    """
    logging.info("Extracting information from headers... ")
    ref_records = dict()
    header_regexes = load_fasta_header_regexes(code_name)
    # TODO: Fix parsing of combined EggNOG and custom headers such that the taxid is parsed from the "accession"
    for treesapp_id in sorted(header_registry.keys(), key=int):
        original_header = header_registry[treesapp_id].original
        header_format_re, header_db, header_molecule = get_header_format(original_header, header_regexes)
        sequence_info = header_format_re.match(original_header)
        seq_info_tuple = return_sequence_info_groups(sequence_info, header_db, original_header)

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


def dedup_records(ref_seqs: FASTA, ref_seq_records: dict):
    if ref_seqs.index_form != "num":
        ref_seqs.change_dict_keys("num")

    deduped_accessions = list()
    all_accessions = dict()
    for treesapp_id in ref_seq_records:
        record = ref_seq_records[treesapp_id]  # type: entrez_utils.EntrezRecord
        if record.accession not in all_accessions:
            all_accessions[record.accession] = []
        all_accessions[record.accession].append(treesapp_id)

    # If the sequences are identical across the records -> take record with greater bitflag
    # If the sequences are different -> keep
    for accession in all_accessions:
        dup_list = all_accessions[accession]
        if len(dup_list) > 1:
            dup_seqs_dict = dict()
            # Load the occurrence of each duplicate record's sequence
            for treesapp_id in dup_list:
                # Could just pop from dup_list if the sequence isn't in dup_seqs_dict but bitflag filter is smarter
                record_seq = ref_seqs.fasta_dict[treesapp_id]
                try:
                    dup_seqs_dict[record_seq] += 1
                except KeyError:
                    dup_seqs_dict[record_seq] = 1
            x = 0
            # Remove the sequences from dup_list if their sequence was found once
            while x < len(dup_list):
                treesapp_id = dup_list[x]
                record_seq = ref_seqs.fasta_dict[treesapp_id]
                if dup_seqs_dict[record_seq] == 1:
                    dup_list.pop(x)
                else:
                    x += 1
            # Move on to the next accession if dup_list is empty (because all duplicates have unique sequences)
            if not dup_list:
                continue

            # Remove records with redundant accessions by removing those with the lower bitflag
            max_bitflag = max([ref_seq_records[treesapp_id].bitflag for treesapp_id in dup_list])
            x = 0
            while x < len(dup_list):
                treesapp_id = dup_list[x]
                ref_seq = ref_seq_records[treesapp_id]
                if ref_seq.bitflag < max_bitflag:
                    ref_seq_records.pop(treesapp_id)
                    dup_list.pop(x)
                    deduped_accessions.append(ref_seqs.header_registry[treesapp_id].original)
                else:
                    x += 1
            x = 1
            # Check for records with identical headers, keeping only one
            while len(dup_list) > 1:
                treesapp_id = dup_list[x]
                ref_seq_records.pop(treesapp_id)
                dup_list.pop(x)  # Don't increment
                deduped_accessions.append(ref_seqs.header_registry[treesapp_id].original)
        else:
            pass

    if deduped_accessions:
        logging.warning("The following sequences were removed during deduplication of accession IDs:\n\t" +
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


def prep_logging(log_file_name=None, verbosity=False) -> None:
    """
    Allows for multiple file handlers to be added to the root logger, but only a single stream handler.
    The new file handlers must be removed outside of this function explicitly

    :param log_file_name:
    :param verbosity:
    :return: None
    """
    if verbosity:
        logging_level = logging.DEBUG
    else:
        logging_level = logging.INFO

    # Detect whether a handlers are already present and return if true
    logger = logging.getLogger()
    if len(logger.handlers):
        return

    formatter = MyFormatter()
    # Set the console handler normally writing to stdout/stderr
    ch = logging.StreamHandler()
    ch.setLevel(logging_level)
    ch.terminator = ''
    ch.setFormatter(formatter)

    if log_file_name:
        output_dir = os.path.dirname(log_file_name)
        try:
            if output_dir and not os.path.isdir(output_dir):
                os.makedirs(output_dir)
        except (IOError, OSError):
            sys.stderr.write("ERROR: Unable to make directory '" + output_dir + "'.\n")
            sys.exit(3)
        logging.basicConfig(level=logging.DEBUG,
                            filename=log_file_name,
                            filemode='w',
                            datefmt="%d/%m %H:%M:%S",
                            format="%(asctime)s %(levelname)s:\n%(message)s")
        logging.getLogger('').addHandler(ch)
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
        self.refpkg_code_re = re.compile(r'[A-Z][0-9]{4,5}')
        self.treesapp_dir = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + os.sep
        self.refpkg_dir = self.treesapp_dir + 'data' + os.sep
        self.tree_dir = self.treesapp_dir + 'data' + os.sep + "tree_data" + os.sep
        self.hmm_dir = self.treesapp_dir + 'data' + os.sep + "hmm_data" + os.sep
        self.aln_dir = self.treesapp_dir + 'data' + os.sep + "alignment_data" + os.sep
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
        self.executables = dict()
        # Values that need to be entered later, in the command-specific class
        self.stages = dict()  # Used to track what progress stages need to be completed
        self.stage_file = ""  # The file to write progress updates to

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

        :param args: arguments from argarse.ParseArgs() with output, input and molecule attributes
        :return: None
        """
        if self.command != "info":
            self.output_dir = args.output
            if self.output_dir[-1] != os.sep:
                self.output_dir += os.sep
            self.final_output_dir = self.output_dir + "final_outputs" + os.sep
            self.var_output_dir = self.output_dir + "intermediates" + os.sep
            if set(vars(args)).issuperset({"molecule", "input"}):
                self.input_sequences = args.input
                self.molecule_type = args.molecule
                self.sample_prefix = '.'.join(os.path.basename(self.input_sequences).split('.')[:-1])
                self.formatted_input = self.var_output_dir + self.sample_prefix + "_formatted.fasta"
        self.executables = self.find_executables(args)
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

    def stage_lookup(self, name: str, tolerant=False):
        """
        Used for looking up a stage in self.stages by its stage.name

        :param name: Name of a stage
        :param tolerant: Boolean controlling whether the function exits if look-up failed (defualt) or returns None
        :return: ModuleFunction instance that matches the name, or None if failed and tolerant, exit otherwise
        """
        fingerprint = False
        for module_step in sorted(self.stages, key=int):
            stage = self.stages[module_step]  # type: ModuleFunction
            if stage.name == name:
                return stage

        if not fingerprint and not tolerant:
            logging.error("Unable to find '" + name + "' in " + self.command + " stages.\n")
            sys.exit(3)
        return

    def first_stage(self):
        for x in sorted(self.stages, key=int):  # type: int
            stage = self.stages[x]  # type: ModuleFunction
            if stage.run:
                return stage.order
        logging.error("No stages are set to run!\n")
        sys.exit(3)

    def read_progress_log(self):
        """
        Read the object's 'stage_file' and determine the completed stage

        :return: An integer corresponding to the last completed stage's rank in 'stage_order'
        """
        completed_int = self.first_stage()
        try:
            progress_handler = open(self.stage_file)
        except IOError:
            logging.debug("Unable to open stage file '" + self.stage_file +
                          "' for reading. Defaulting to stage " + str(completed_int) + ".\n")
            return completed_int
        # TODO: Finish this to enable continuing part-way through an analysis
        progress_handler.close()
        return completed_int

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

    def validate_continue(self, args):
        """
        Bases the stage(s) to run on args.stage which is broadly set to either 'continue' or any other valid stage

        This function ensures all the required inputs are present for beginning at the desired first stage,
        otherwise, the pipeline begins at the first possible stage to continue and ends once the desired stage is done.

        If a
        :return: None
        """
        if self.command == "assign":
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
        elif self.command == "create":
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
                    logging.error("Unable to find accession-lineage mapping file '" + self.acc_to_lin + "'\n")
                    sys.exit(3)
            else:
                self.acc_to_lin = self.var_output_dir + os.sep + "accession_id_lineage_map.tsv"

        elif self.command == "evaluate":
            pass
        elif self.command == "train":
            if not args.profile:
                self.change_stage_status("search", False)
            # TODO: Remove duplicate code once the log-file parsing is implemented
            if args.acc_to_lin:
                self.acc_to_lin = args.acc_to_lin
                if os.path.isfile(self.acc_to_lin):
                    self.change_stage_status("lineages", False)
                else:
                    logging.error("Unable to find accession-lineage mapping file '" + self.acc_to_lin + "'\n")
                    sys.exit(3)
            else:
                self.acc_to_lin = self.var_output_dir + os.sep + "accession_id_lineage_map.tsv"
        elif self.command == "update":
            self.acc_to_lin = self.var_output_dir + os.sep + "accession_id_lineage_map.tsv"
        elif self.command == "purity":
            pass
        else:
            logging.error("Unknown sub-command: " + str(self.command) + "\n")
            sys.exit(3)

        # TODO: Summarise the steps to be taken and write to log
        if args.overwrite:
            last_valid_stage = self.first_stage()
        else:
            last_valid_stage = self.read_progress_log()

        if args.stage == "continue":
            logging.debug("Continuing with stage '" + self.stages[last_valid_stage].name + "'\n")
            self.edit_stages(last_valid_stage)
            return

        # Update the stage status
        desired_stage = self.stage_lookup(args.stage).order
        if desired_stage > last_valid_stage:
            logging.warning("Unable to run '" + args.stage + "' as it is ahead of the last completed stage.\n" +
                            "Continuing with stage '" + self.stages[last_valid_stage].name + "'\n")
            self.edit_stages(last_valid_stage, desired_stage)
        elif desired_stage < last_valid_stage:
            logging.info("Wow - its your lucky day:\n"
                         "All stages up to and including '" + args.stage + "' have already been completed!\n")
            self.edit_stages(desired_stage, desired_stage)
        else:
            # Proceed with running the desired stage
            logging.debug("Proceeding with '" + args.stage + "'\n")
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

        if self.command == "update":
            dependencies += ["usearch", "mafft", "OD-seq"]

        if self.command == "create":
            dependencies += ["usearch", "mafft", "OD-seq"]
            if args.fast:
                dependencies.append("FastTree")

        if self.molecule_type == "rrna":
            dependencies += ["cmalign", "cmsearch", "cmbuild"]

        for dep in dependencies:
            # For rpkm and potentially other executables that are compiled ad hoc
            if is_exe(self.treesapp_dir + "sub_binaries" + os.sep + dep):
                exec_paths[dep] = str(self.treesapp_dir + "sub_binaries" + os.sep + dep)
            elif which(dep):
                exec_paths[dep] = which(dep)
            else:
                logging.error("Could not find a valid executable for " + dep + ".\n")
                sys.exit(13)

        return exec_paths

    def fetch_entrez_lineages(self, ref_seqs: FASTA, molecule, acc_to_taxid=None):
        # Get the lineage information for the training/query sequences
        ref_seq_records = get_header_info(ref_seqs.header_registry, self.ref_pkg.prefix)
        ref_seq_records = dedup_records(ref_seqs, ref_seq_records)
        ref_seqs.change_dict_keys("formatted")
        entrez_utils.load_ref_seqs(ref_seqs.fasta_dict, ref_seqs.header_registry, ref_seq_records)
        logging.debug("\tNumber of input sequences =\t" + str(len(ref_seq_records)) + "\n")

        if self.stage_status("lineages"):
            entrez_query_list, num_lineages_provided = entrez_utils.build_entrez_queries(ref_seq_records)
            logging.debug("\tNumber of queries =\t" + str(len(entrez_query_list)) + "\n")
            entrez_records = entrez_utils.map_accessions_to_lineages(entrez_query_list, self.ref_pkg.taxa_trie,
                                                                     molecule, acc_to_taxid)
            self.seq_lineage_map = entrez_utils.entrez_records_to_accession_lineage_map(entrez_records)
            # Download lineages separately for those accessions that failed
            # Map proper accession to lineage from the tuple keys (accession, accession.version)
            #  in accession_lineage_map returned by entrez_utils.get_multiple_lineages.
            entrez_utils.repair_lineages(ref_seq_records, self.ref_pkg.taxa_trie)
            entrez_utils.verify_lineage_information(self.seq_lineage_map, ref_seq_records,
                                                    self.ref_pkg.taxa_trie, num_lineages_provided)
            self.seq_lineage_map = entrez_utils.accession_lineage_map_from_entrez_records(ref_seq_records)

            # Ensure the accession IDs are stripped of '>'s
            for accession in sorted(self.seq_lineage_map):
                if accession[0] == '>':
                    self.seq_lineage_map[accession[1:]] = self.seq_lineage_map.pop(accession)
            write_dict_to_table(self.seq_lineage_map, self.acc_to_lin)
            # Add lineage information to the ReferenceSequence() objects in ref_seq_records if not contained
        else:
            logging.info("Reading cached lineages in '" + self.acc_to_lin + "'... ")
            self.seq_lineage_map.update(entrez_utils.read_accession_taxa_map(self.acc_to_lin))
            logging.info("done.\n")
        ref_seqs.change_dict_keys()
        return ref_seq_records


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
        self.uclust_prefix = ""  # Used only if resolve is True
        self.target_marker = None
        self.rank_depth_map = None
        self.prop_sim = 1.0
        self.min_length = 0  # The minimum sequence length for a classified sequence to be included in the refpkg

        # Stage names only holds the required stages; auxiliary stages (e.g. RPKM, update) are added elsewhere
        self.stages = {0: ModuleFunction("lineages", 0),
                       1: ModuleFunction("rebuild", 1)}

    def get_info(self):
        info_string = "Updater instance summary:\n"
        info_string += super(Updater, self).get_info() + "\n\t"
        info_string += "\n\t".join(["Target reference packages = " + str(self.ref_pkg.prefix),
                                    "Taxonomy map: " + self.ref_pkg.lineage_ids,
                                    "Reference tree: " + self.ref_pkg.tree,
                                    "Reference FASTA: " + self.ref_pkg.msa,
                                    "Lineage map: " + str(self.seq_names_to_taxa)]) + "\n"

        return info_string

    def map_orf_lineages(self, seq_lineage_map: dict, header_registry: dict) -> (dict, list):
        """
        The classified sequences have a signature at the end (|RefPkg_name|start_stop) that needs to be removed
        Iterates over the dictionary of sequence names and attempts to match those with headers in the registry.
        If a match is found the header is assigned the corresponding lineage in seq_lineage_map.

        :param seq_lineage_map: A dictionary mapping contig sequence names to taxonomic lineages
        :param header_registry: A dictionary mapping numerical TreeSAPP identifiers to Header instances
        :return: A dictionary mapping each classified sequence to a lineage and list of TreeSAPP IDs that were mapped
        """
        logging.info("Mapping assigned sequences to provided taxonomic lineages... ")
        classified_seq_lineage_map = dict()
        treesapp_nums = list(header_registry.keys())
        mapped_treesapp_nums = []
        for seq_name in seq_lineage_map:
            # Its slow to perform so many re.search's but without having a guaranteed ORF pattern
            # we can't use hash-based data structures to bring it to O(N)
            parent_re = re.compile(seq_name)
            x = 0
            while x < len(treesapp_nums):
                header = header_registry[treesapp_nums[x]].original
                assigned_seq_name = re.sub(r"\|{0}\|\d+_\d+.*".format(self.ref_pkg.prefix), '', header)
                if parent_re.search(assigned_seq_name):
                    classified_seq_lineage_map[header] = seq_lineage_map[seq_name]
                    mapped_treesapp_nums.append(treesapp_nums.pop(x))
                else:
                    x += 1
            if len(treesapp_nums) == 0:
                logging.info("done.\n")
                return classified_seq_lineage_map, mapped_treesapp_nums
        logging.debug("Unable to find parent for " + str(len(treesapp_nums)) + " ORFs in sequence-lineage map:\n" +
                      "\n".join([header_registry[n].original for n in treesapp_nums]) + "\n")
        return classified_seq_lineage_map, mapped_treesapp_nums


class Creator(TreeSAPP):
    def __init__(self):
        super(Creator, self).__init__("create")
        self.prop_sim = 1.0
        self.candidates = dict()  # Dictionary tracking all candidate ReferenceSequences
        self.phy_dir = ""  # Directory for intermediate or unnecessary files created during phylogeny inference
        self.hmm_purified_seqs = ""  # If an HMM profile of the gene is provided its a path to FASTA with homologs
        self.filtered_fasta = ""
        self.hmm_profile = ""  # HMM profile used for screening the input sequences
        self.uclust_prefix = ""  # FASTA file prefix for cluster centroids
        self.uc = ""  # UCLUST-specific output file defining the clusters, members, representatives, etc.
        self.cluster_input = ""  # Name of the file to be used for clustering
        self.unaln_ref_fasta = ""  # FASTA file of unaligned reference sequences
        self.phylip_file = ""  # Used for building the phylogenetic tree with RAxML
        self.min_tax_rank = "Kingdom"  # Minimum taxonomic rank
        self.metadata_file = ""

        # Stage names only holds the required stages; auxiliary stages (e.g. RPKM, update) are added elsewhere
        self.stages = {0: ModuleFunction("search", 0),
                       1: ModuleFunction("lineages", 1),
                       2: ModuleFunction("clean", 2),
                       3: ModuleFunction("cluster", 3),
                       4: ModuleFunction("build", 4),
                       5: ModuleFunction("train", 5),
                       6: ModuleFunction("update", 6)}

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

    def determine_model(self, fast):
        model = ""
        if fast:
            if self.molecule_type == "prot":
                model = "LG"
            else:
                model = "GTRGAMMA"
        else:
            raxml_info_file = self.phy_dir + "RAxML_info." + self.ref_pkg.prefix
            model_statement_re = re.compile(r".* model: ([A-Z]+) likelihood.*")
            command_line = ""
            with open(raxml_info_file) as raxml_info:
                for line in raxml_info:
                    if model_statement_re.search(line):
                        model = model_statement_re.search(line).group(1)
                        break
                    elif re.match('^.*/raxml.*-m ([A-Z]+)$', line):
                        command_line = line
                    else:
                        pass
            if model == "":
                if command_line == "":
                    logging.warning("Unable to parse model used from " + raxml_info_file + "!\n")
                else:
                    model = re.match('^.*/raxml.*-m ([A-Z]+)$', command_line).group(1)
        if self.molecule_type == "prot":
            model = "PROTGAMMA" + model
        self.ref_pkg.sub_model = model
        return model

    def print_terminal_commands(self):
        param_file = self.treesapp_dir + "data" + os.sep + "ref_build_parameters.tsv"
        logging.info("\nTo integrate this package for use in TreeSAPP the following steps must be performed:\n" +
                     "1. Write a properly formatted reference package 'code' in " + param_file + "\n" +
                     "2. $ cp " + ' '.join([self.ref_pkg.lineage_ids, self.ref_pkg.tree, self.ref_pkg.model_info]) +
                     ' ' + self.tree_dir + "\n" +
                     "3. $ cp " + self.ref_pkg.profile + self.ref_pkg.search_profile + ' ' + self.hmm_dir + "\n" +
                     "4. $ cp " + self.ref_pkg.msa + ' ' + self.aln_dir + "\n")
        return


class Purity(TreeSAPP):
    def __init__(self):
        super(Purity, self).__init__("purity")
        self.assign_dir = ""
        self.classifications = ""
        self.summarize_dir = ""
        self.metadata_file = ""
        self.assignments = None
        self.refpkg_build = MarkerBuild()
        self.stages = {0: ModuleFunction("assign", 0),
                       1: ModuleFunction("summarize", 1)}

    def summarize_groups_assigned(self, ortholog_map: dict, metadata=None):
        unique_orthologs = dict()
        tree_leaves = self.ref_pkg.tax_ids_file_to_leaves()
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
        xtra_dat = namedtuple("xtra_dat", ["id", "ac", "de"])
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
                id, accession, desc = line.strip().split("\t")
            except ValueError:
                logging.error("Bad format for '" + self.metadata_file + "'. Three tab-separated fields expected.\n" +
                              "Example line:\n" +
                              str(line) + "\n")
                sys.exit(7)
            metadat_dict[accession] = xtra_dat(id, accession, desc)
        metadata_handler.close()
        return metadat_dict

    def assign_leaves_to_orthologs(self, p_queries: list, internal_node_map: dict) -> dict:
        ortholog_map = dict()
        leaf_map = dict()
        tree_leaves = self.ref_pkg.tax_ids_file_to_leaves()
        for leaf in tree_leaves:
            leaf_map[leaf.number + "_" + self.ref_pkg.prefix] = leaf.description
        for p_query in p_queries:  # type: TreeProtein
            p_query.name = self.ref_pkg.prefix
            if type(p_query.contig_name) is list and len(p_query.contig_name):
                p_query.contig_name = p_query.contig_name[0]
            seq_info = re.match(r"(.*)\|" + re.escape(p_query.name) + r"\|(\\d+)_(\\d+)$", p_query.contig_name)
            if seq_info:
                p_query.contig_name = seq_info.group(1)
            p_query.inode = int(p_query.get_jplace_element("edge_num"))
            leaves = internal_node_map[p_query.inode]
            ortholog_name = p_query.contig_name.split('_')[0]
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
        if self.lineage.split('; ')[0] != "Root":
            self.lineage = "; ".join(["Root"] + self.lineage.split("; "))
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
        :return:
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
        self.targets = []  # Left empty to appease parse_ref_build_parameters()
        self.target_marker = MarkerBuild()
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

    def get_info(self):
        info_string = "Evaluator instance summary:\n"
        info_string += super(Evaluator, self).get_info() + "\n\t"
        info_string += "\n\t".join(["Accession-to-lineage map = " + self.acc_to_lin,
                                    "Clade-exclusion table = " + self.performance_table,
                                    "Target marker(s) = " + str(self.targets)]) + "\n"

        return info_string

    def new_taxa_test(self, rank, lineage) -> TaxonTest:
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

        sys.stdout.write("Rank-level performance of " + self.target_marker.cog + ":\n")
        sys.stdout.write("\tRank\tQueries\tClassified\tCorrect\tD=1\tD=2\tD=3\tD=4\tD=5\tD=6\tD=7\n")

        for depth in sorted(self.rank_depth_map):
            rank = self.rank_depth_map[depth]
            if rank == "Cellular organisms":
                continue
            taxonomic_distance = dict()
            n_queries, n_classified, sensitivity = self.get_sensitivity(rank)
            for dist in range(0, 8):
                taxonomic_distance[dist] = EvaluateStats(self.target_marker.cog, rank, dist)
                taxonomic_distance[dist].n_qs = n_queries
            std_out_report_string += "\t" + rank + "\t"
            if rank not in rank_assigned_dict or len(rank_assigned_dict[rank]) == 0:
                std_out_report_string += "0\t0\t\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n"
            else:
                acc = 0
                for assignments in rank_assigned_dict[rank]:
                    for classified in assignments:
                        acc += 1
                        if classified.split("; ")[0] == "Cellular organisms":
                            logging.error("Lineage string cleaning has gone awry somewhere. "
                                          "The root rank should be a Kingdom (e.g. Bacteria or Archaea) but nope.\n")
                            sys.exit(21)
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
            line = trial_name + "\t" + self.target_marker.cog + "\t" + tool + line + "\n"
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

        test_obj = self.new_taxa_test(rank, lineage)
        test_obj.intermediates_dir = self.var_output_dir + refpkg_name + os.sep + rank_tax + os.sep
        if not os.path.isdir(test_obj.intermediates_dir):
            os.makedirs(test_obj.intermediates_dir)

        logging.info("Classifications for the " + rank + " '" + taxon + "' put " + test_obj.intermediates_dir + "\n")
        test_obj.test_query_fasta = test_obj.intermediates_dir + rank_tax + ".fa"
        test_obj.test_tax_ids_file = test_obj.intermediates_dir + "tax_ids_" + refpkg_name + ".txt"
        test_obj.classifier_output = test_obj.intermediates_dir + "TreeSAPP_output" + os.sep
        return test_obj

    def prep_for_clade_exclusion(self, refpkg: ReferencePackage, taxa_test: TaxonTest, lineage: str, lineage_seqs: dict,
                                 trim_align=False, no_svm=False, min_seq_len=0, targeted=False, num_threads=2) -> list:
        """
        Calls ReferencePackage.exlude_clade_from_ref_files() to remove all reference sequences/leaf nodes from the
        reference package that are descendents of `lineage`.

        Creates the TreeSAPP assign arguments list, which is to be called outside of this function.

        :param refpkg: The ReferencePackage object that is to be used for clade exclusion, of which reference sequences
         that are descendents of `lineage` will be removed.
        :param lineage: Taxonomic lineage which is going to be used in this clade exclusion analysis
        :param lineage_seqs: A fasta-like dictionary where headers are sequence names (accessions) and values are
         their respective nucleotide or amino acid sequences.
        :param taxa_test:
        :param trim_align: Flag determining whether TreeSAPP should use BMGE to trim the multiple sequence alignments
         prior to phylogenetic placement
        :param no_svm: Flag controlling whether SVM-filtering is applied to placements
        :param min_seq_len: The minimum sequence length argument for TreeSAPP
        :param num_threads: The maximum number or threads and processes TreeSAPP and its dependencies can use
        :param targeted: Boolean controlling whether homology search against the query sequences just uses this
         reference package's HMM (True) or all HMMs available in TreeSAPP (False)
        :return: A list of arguments to be called by TreeSAPP's assign function and a TaxonTest object
        """
        # Copy reference files, then exclude all clades belonging to the taxon being tested
        taxa_test.temp_files_prefix = refpkg.exclude_clade_from_ref_files(self.refpkg_dir, self.molecule_type,
                                                                          self.var_output_dir + refpkg.prefix + os.sep,
                                                                          lineage, self.executables)
        # Write the query sequences
        taxa_test.queries = lineage_seqs.keys()
        write_new_fasta(lineage_seqs, taxa_test.test_query_fasta)

        assign_args = ["-i", taxa_test.test_query_fasta, "-o", taxa_test.classifier_output,
                       "-m", self.molecule_type, "-n", str(num_threads),
                       "--overwrite", "--delete"]
        if trim_align:
            assign_args.append("--trim_align")
        if no_svm:
            assign_args.append("--no_svm")
        if min_seq_len:
            assign_args += ["--min_seq_length", str(min_seq_len)]
        if targeted:
            assign_args += ["--targets", refpkg.refpkg_code]

        return assign_args


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
        self.fq_suffix_re = re.compile(r"([._-])+[pe|fq|fastq|fwd|R1]+$")

    def check_arguments(self, args):
        ##
        # Define locations of files TreeSAPP outputs
        ##
        self.classified_nuc_seqs = glob(self.final_output_dir + "*_classified.fna")[0]
        if not os.path.isfile(self.classified_nuc_seqs):
            logging.error("Unable to find classified sequences FASTA file in %s.\n" % self.final_output_dir)
        self.classifications = self.output_dir + "final_outputs" + os.sep + "marker_contig_map.tsv"

        if not os.path.isdir(self.var_output_dir):
            os.makedirs(self.var_output_dir)

        return

    def assignments_to_treesaps(self, classified_lines: list, marker_build_dict: dict) -> dict:
        """

        :param classified_lines:
        :param marker_build_dict:
        :return: A dictionary of ItolJplace instances, indexed by their respective RefPkg codes (denominators)
        """
        pqueries = dict()
        for fields in classified_lines:
            tree_sap = ItolJplace()
            try:
                _, tree_sap.contig_name, tree_sap.name, tree_sap.seq_len, tree_sap.lct, tree_sap.recommended_lineage,\
                _, tree_sap.inode, tree_sap.lwr, tree_sap.avg_evo_dist, tree_sap.distances = fields
            except ValueError:
                logging.error("Bad line in classification table:\n" + '\t'.join(fields) + "\n")
                sys.exit(21)
            refpkg = fish_refpkg_from_build_params(tree_sap.name, marker_build_dict).denominator
            try:
                pqueries[refpkg].append(tree_sap)
            except KeyError:
                pqueries[refpkg] = [tree_sap]
        return pqueries


class Assigner(TreeSAPP):
    def __init__(self):
        """

        """
        super(Assigner, self).__init__("assign")
        self.reference_tree = None
        self.aa_orfs_file = ""
        self.nuc_orfs_file = ""
        self.classified_aa_seqs = ""
        self.classified_nuc_seqs = ""
        self.composition = ""
        self.target_refpkgs = list()
        self.clf = load_pickle(self.refpkg_dir + "treesapp_svm.pkl")

        # Stage names only holds the required stages; auxiliary stages (e.g. RPKM, update) are added elsewhere
        self.stages = {0: ModuleFunction("orf-call", 0, self.predict_orfs),
                       1: ModuleFunction("clean", 1, self.clean),
                       2: ModuleFunction("search", 2, self.search),
                       3: ModuleFunction("align", 3, self.align),
                       4: ModuleFunction("place", 4, self.place),
                       5: ModuleFunction("classify", 5, self.classify)}

    def get_info(self):
        info_string = "Assigner instance summary:\n"
        info_string += super(Assigner, self).get_info() + "\n\t"
        info_string += "\n\t".join(["ORF protein sequences = " + self.aa_orfs_file,
                                    "Target reference packages = " + str(self.target_refpkgs),
                                    "Composition of input = " + self.composition]) + "\n"

        return info_string

    def predict_orfs(self, composition, num_threads):
        """
        Predict ORFs from the input FASTA file using Prodigal

        :param composition: Sample composition being either a single organism or a metagenome [single | meta]
        :param num_threads: The number of CPU threads to use
        :return: None
        """

        logging.info("Predicting open-reading frames in the genomes using Prodigal... ")

        start_time = time.time()

        if num_threads > 1 and composition == "meta":
            # Split the input FASTA into num_threads files to run Prodigal in parallel
            input_fasta_dict = format_read_fasta(self.input_sequences, self.molecule_type, self.output_dir)
            n_seqs = len(input_fasta_dict.keys())
            chunk_size = int(n_seqs / num_threads) + (n_seqs % num_threads)
            split_files = write_new_fasta(input_fasta_dict,
                                          self.var_output_dir + self.sample_prefix,
                                          chunk_size)
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

        # Concatenate outputs
        if not os.path.isfile(self.aa_orfs_file) and not os.path.isfile(self.nuc_orfs_file):
            tmp_prodigal_aa_orfs = glob(self.var_output_dir + self.sample_prefix + "*_ORFs.faa")
            tmp_prodigal_nuc_orfs = glob(self.var_output_dir + self.sample_prefix + "*_ORFs.fna")
            os.system("cat " + ' '.join(tmp_prodigal_aa_orfs) + " > " + self.aa_orfs_file)
            os.system("cat " + ' '.join(tmp_prodigal_nuc_orfs) + " > " + self.nuc_orfs_file)
            intermediate_files = list(tmp_prodigal_aa_orfs + tmp_prodigal_nuc_orfs + split_files)
            for tmp_file in intermediate_files:
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

        # Limit this to just Class, Family, and Species - other ranks are inferred through regression
        self.training_ranks = {"Class": 3, "Species": 7}

        # Stage names only holds the required stages; auxiliary stages (e.g. RPKM, update) are added elsewhere
        self.stages = {0: ModuleFunction("search", 0),
                       1: ModuleFunction("lineages", 1),
                       2: ModuleFunction("place", 2),
                       3: ModuleFunction("regress", 3)}

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
