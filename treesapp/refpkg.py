import logging
import os
import re
import sys
import json
import inspect
from shutil import copy

from packaging import version
from ete3 import Tree
import joblib

from treesapp.phylo_seq import TreeLeafReference
from treesapp.entish import annotate_partition_tree, label_internal_nodes_ete, verify_bifurcations
from treesapp.external_command_interface import launch_write_command
from treesapp.fasta import read_fasta_to_dict, write_new_fasta, multiple_alignment_dimensions, FASTA, register_headers
from treesapp.taxonomic_hierarchy import TaxonomicHierarchy, Taxon
from treesapp.utilities import base_file_prefix, load_taxonomic_trie, match_file, get_hmm_value
from treesapp import wrapper
from treesapp import __version__ as ts_version


_COMPATIBLE_VERSION = "0.8.3"


class ReferencePackage:
    """
    The ReferencePackage class is a collection of all attributes required by a TreeSAPP reference package as well as
    the suite of functions for building, manipulating, and accessing different attributes.
    """
    def __init__(self, refpkg_name=""):
        self.prefix = refpkg_name
        self.refpkg_code = "Z1111"  # AKA denominator

        # These are files (with '_f' suffix) and their respective data (read with file.readlines())
        # TODO: Rename to f__pkl
        self.refpkg_suffix = "_build.pkl"
        self.f__json = self.prefix + self.refpkg_suffix  # Path to the pickled reference package file
        self.msa = []  # Reference MSA FASTA
        self.f__msa = self.prefix + ".fa"
        self.profile = []
        self.f__profile = self.prefix + ".hmm"  # HMM file
        self.search_profile = []
        self.f__search_profile = self.prefix + '_' + "search.hmm"  # profile HMM that has been dereplicated
        self.tree = []  # Reference tree
        self.f__tree = self.prefix + ".nwk"
        self.boot_tree = []  # Reference tree with support values
        self.f__boot_tree = self.prefix + "_bipart.nwk"
        self.model_info = []
        self.f__model_info = self.prefix + "_epa.model"  # RAxML-NG --evaluate model file
        self.svc = None
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
        self.__core_ref_files = [self.f__msa, self.f__profile, self.f__search_profile,
                                 self.f__tree, self.f__boot_tree, self.f__model_info]
        self.taxa_trie = TaxonomicHierarchy()

    def __iter__(self):
        """
        Used for iterative through an instance's dictionary of attributes, returning both the key and value pairs

        :return: A tuple of a attribute name (key) and its corresponding value
        """
        for attr, value in self.__dict__.items():
            yield attr, value

    def get_public_attributes(self) -> list:
        return [attr[0] for attr in inspect.getmembers(self) if
                type(attr[1]) in [int, float, str, list, dict] and
                not attr[0].startswith('_')]

    def get_info(self):
        return "\n\t".join(["ReferencePackage instance of {} ({}):".format(self.prefix, self.refpkg_code),
                            "Molecule type:                                      '{}'".format(self.molecule),
                            "TreeSAPP version:                                   '{}'".format(self.ts_version),
                            "Profile HMM length:                                 '{}'".format(self.profile_length),
                            "Substitution model used for phylogenetic inference: '{}'".format(self.sub_model),
                            "Number of reference sequences (leaf nodes):          {}".format(self.num_seqs),
                            "Software used to infer phylogeny:                   '{}'".format(self.tree_tool),
                            "Date of last update:                                '{}'".format(self.update),
                            "Description:                                        '{}'\n".format(self.description)
                            ])

    def bail(self, msg=""):
        """
        Function to consistently bail on processing from within ReferencePackage operations.

        :return: None
        """
        logging.error(msg + "\n" + self.get_info())
        return

    def clone(self, clone_path: str):
        refpkg_clone = ReferencePackage()
        if not os.path.isfile(self.f__json):
            self.pickle_package()
        refpkg_clone.f__json = self.f__json
        refpkg_clone.slurp()
        refpkg_clone.f__json = clone_path
        return refpkg_clone

    def pickle_package(self) -> None:
        """
        Dumps a new joblib pickled file of the ReferencePackage instance dictionary to self.f__json.

        :return: None.
        """
        if len(self.f__json) == 0:
            self.bail("ReferencePackage.f__json not set. ReferencePackage band() cannot be completed.\n")
            raise AttributeError

        if not os.path.isdir(os.path.dirname(self.f__json)):
            os.mkdir(os.path.dirname(self.f__json))

        try:
            refpkg_handler = open(self.f__json, 'wb')
        except IOError:
            self.bail("Unable to open reference package pickled file '{}' for writing.\n".format(self.f__json))
            raise IOError

        refpkg_dict = {}
        non_primitives = ["__core_ref_files", "taxa_trie"]
        for a, v in self.__iter__():
            if a not in non_primitives:
                refpkg_dict[a] = v

        joblib.dump(value=refpkg_dict, filename=refpkg_handler)

        refpkg_handler.close()
        return

    def band(self) -> None:
        """
        Reads each of the individual reference package component files (e.g. 'f__msa', 'f__tree', 'f__model_info') and
        writes a pickle file containing all of the reference package files and metadata.

        :return: None
        """
        # Read the all of the individual reference package files that are available (e.g. MSA, HMM, phylogeny)
        for a, v in self.__iter__():
            if a.startswith('f__') and os.path.isfile(v):
                dest = a.lstrip("f__")
                if dest in self.__dict__:
                    with open(v) as fh:
                        self.__dict__[dest] = fh.readlines()
        self.pickle_package()

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

    def update_file_names(self) -> None:
        for a, v in self.__iter__():
            if a.startswith('f__'):
                path, name = os.path.split(v)
                self.__dict__[a] = os.path.join(path, self.prefix + re.sub(self.prefix, '', name))
        return

    def disband(self, output_dir: str) -> None:
        """
        From a ReferencePackage's pickled file, the individual file components (e.g. profile HMM, MSA, phylogeny) are
        written to their separate files in a new directory, a sub-directory of where the pickle is located.
        The directory name follows the format: self.prefix + '_' self.refpkg_code + '_' + self.date

        Their file paths (e.g. self.f__msa, self.f__tree) are updated with the output_dir as the directory path.

        :param output_dir: Output directory to write the individual reference package files
        :return: None
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
        Reads the reference package's pickled-formatted file and stores the elements in their respective variables

        :return: None
        """
        new_path = self.f__json
        if len(self.f__json) == 0:
            logging.error("ReferencePackage.f__json was not set.\n")
            sys.exit(11)

        if not os.path.isfile(self.f__json):
            logging.error("ReferencePackage pickle file '{}' doesn't exist.\n".format(self.f__json))
            sys.exit(7)

        try:
            refpkg_data = joblib.load(self.f__json)
        except KeyError:
            refpkg_handler = open(self.f__json, 'r')
            refpkg_data = json.load(refpkg_handler)
            refpkg_handler.close()
        except EOFError:
            logging.error("Joblib was unable to load reference package pickle '{}'.\n".format(self.f__json))
            sys.exit(17)

        for a, v in refpkg_data.items():
            self.__dict__[a] = v

        # Fix the pickle path
        self.f__json = new_path

        if type(self.tree) is list:
            if len(self.tree) > 0:
                self.tree = self.tree[0]
            else:
                logging.error("Unable to load the phylogenetic tree for reference package '{}' ({})\n."
                              "".format(self.prefix, self.f__json))
                return

        self.load_taxonomic_hierarchy()

        return

    def validate(self, check_files=False):
        """
        Function that ensures the number of sequences is equal across all files and
        the version of TreeSAPP used to create this reference package is compatible with the current version.

        :return: Boolean
        """
        # Check to ensure all files exist
        if check_files:
            for ref_file in self.__core_ref_files:
                if not os.path.isfile(ref_file):
                    self.bail("File '{}' does not exist for ReferencePackage {}\n".format(ref_file, self.prefix))
                    return False

        if version.parse(self.ts_version) < version.parse(_COMPATIBLE_VERSION):
            self.bail("'{}' reference package (created with version '{}')"
                      " is not compatible with this version of TreeSAPP ('{}').\n".format(self.prefix,
                                                                                          self.ts_version, ts_version))
            return False

        # Compare the number of sequences in the multiple sequence alignment
        refpkg_fa = self.get_fasta()  # type: FASTA
        if refpkg_fa.n_seqs() == 0:
            return False
        if self.num_seqs != refpkg_fa.n_seqs():
            self.bail("Number of sequences in {} "
                      "ReferencePackage.num_seqs ({}) and MSA ({}) differ.\n".format(self.prefix,
                                                                                     self.num_seqs,
                                                                                     refpkg_fa.n_seqs()))
            return False

        # Compare the number of sequences in the Hidden-Markov model
        hmm_nseq = self.num_hmm_seqs()
        if hmm_nseq != self.num_seqs:
            self.bail("Number of sequences in {} "
                      "ReferencePackage.num_seqs ({}) and HMM ({}) differ.\n".format(self.prefix,
                                                                                     self.num_seqs,
                                                                                     hmm_nseq))
            return False

        # Compare the number of sequences in the Tree files
        rt = self.get_ete_tree()
        if len(rt) != self.num_seqs:
            self.bail("Number of sequences in {} "
                      "ReferencePackage.num_seqs ({}) and tree ({}) differ.\n".format(self.prefix,
                                                                                      self.num_seqs,
                                                                                      len(rt)))
            return False
        self.tree = verify_bifurcations(self.tree)

        # Compare the number of sequences in lineage IDs
        n_leaf_nodes = len(self.generate_tree_leaf_references_from_refpkg())
        if n_leaf_nodes != self.num_seqs:
            self.bail("Number of sequences in {} "
                      "ReferencePackage.num_seqs ({}) and leaf nodes ({}) differ.\n".format(self.prefix,
                                                                                            self.num_seqs,
                                                                                            n_leaf_nodes))
            return False

        if not self.sub_model:
            self.bail("Unable to find the substitution model used for ReferencePackage '{}'.\n".format(self.prefix))
            return False

        return True

    def get_fasta(self) -> FASTA:
        """
        Used for reading the reference package's self.msa lines to instantiate a FASTA instance

        :return: A FASTA instance populated by the reference package's msa attribute
        """
        refpkg_fa = FASTA(self.f__msa)
        if not self.msa:
            logging.debug("ReferencePackage MSA hasn't been slurped, unable to read FASTA.\n")
            return refpkg_fa

        name, seq = "", ""
        for line in self.msa:
            line = line.rstrip()  # Strip the newline character
            if not line:
                continue
            if line[0] == '>':
                if name:  # This is the first sequence
                    refpkg_fa.fasta_dict[name] = seq
                name, seq = "", ""
                name = line[1:]
            else:
                seq += line
        refpkg_fa.fasta_dict[name] = seq
        refpkg_fa.header_registry = register_headers(list(refpkg_fa.fasta_dict.keys()), True)

        if len(refpkg_fa.fasta_dict) == 0 and len(refpkg_fa.fasta_dict) == 0:
            logging.error("ReferencePackage.msa is empty or corrupted - no sequences were found!\n")
            sys.exit(3)

        return refpkg_fa

    def get_ete_tree(self) -> Tree:
        if not self.tree:
            logging.error("Unable to load tree - '{}' reference package hasn't been slurped yet.\n".format(self.prefix))
            sys.exit(5)
        if type(self.tree) is list:
            if len(self.tree) > 0:
                self.tree = self.tree.pop(0)

        rt = Tree(self.tree)
        label_internal_nodes_ete(rt)
        return rt

    def get_internal_node_leaf_map(self) -> dict:
        node_map = dict()
        leaf_stack = list()
        rt = self.get_ete_tree()

        i = 0
        for inode in rt.traverse(strategy="postorder"):
            if inode.children:
                node_map[i] = leaf_stack.pop() + leaf_stack.pop()
                leaf_stack.append(node_map[i])
            else:
                node_map[i] = [inode.name]
                leaf_stack.append(node_map[i])
            i += 1
        node_map[i] = leaf_stack.pop()
        if leaf_stack:
            logging.error("Some leaves remain unpopped from the stack - internal node to leaf map is incomplete.\n")
            sys.exit(17)

        return node_map

    def leaf_node_order(self) -> list:
        """
        Parses a tree by depth-first search and returns the list

        :return: A list of internal node identifiers sorted by post-order traversal
        """
        tree_root = self.get_ete_tree()
        leaf_order = []
        for node in tree_root.traverse(strategy="postorder"):
            if node.name:
                leaf_order.append(node.name)
        return leaf_order

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

    def create_itol_labels(self, output_dir) -> None:
        """
        Create the marker_labels.txt file for each marker gene that was used for classification

        :param output_dir: Path to the directory for where these outputs should be written to
        :return: None
        """
        itol_label_file = os.path.join(output_dir, self.prefix, self.prefix + "_labels.txt")

        if os.path.exists(itol_label_file):
            return

        try:
            label_f = open(itol_label_file, 'w')
        except IOError:
            raise IOError("Unable to open " + itol_label_file + " for writing! Exiting now.")

        leaf_node_lineages = self.generate_tree_leaf_references_from_refpkg()
        label_f.write("LABELS\nSEPARATOR COMMA\nDATA\n#NODE_ID,LABEL\n")
        for leaf_node in leaf_node_lineages:  # type: TreeLeafReference
            label_f.write(leaf_node.number + '_' + self.prefix + ',' + leaf_node.description + "\n")

        label_f.close()

        return

    def hmm_length(self) -> int:
        if not self.profile:
            self.bail("ReferencePackage attribute 'profile' was not set when trying to find profile HMM length.\n")
            raise ValueError
        self.profile_length = int(get_hmm_value(self.profile, "length"))
        return self.profile_length

    def num_hmm_seqs(self) -> int:
        return int(get_hmm_value(self.profile, "num_seqs"))

    def alignment_dims(self):
        return multiple_alignment_dimensions(self.f__msa)

    def screen_refs_by_lineage(self, target_lineage: str) -> list:
        """
        Returns the TreeLeafReference instances from the reference package that match or are children of target_lineage

        :param target_lineage: A string representing the lineage to fish for reference leaf nodes (e.g. d__Archaea)
        :return: A list of all TreeLeafReference instances that belong to the target_lineage
        """
        children = []
        lineage_re = re.compile(target_lineage)
        for leaf_node in self.generate_tree_leaf_references_from_refpkg():  # type: TreeLeafReference
            if lineage_re.search(leaf_node.lineage):
                children.append(leaf_node)
        return children

    def get_hierarchy_taxon_for_leaf(self, leaf: TreeLeafReference):
        if not leaf.lineage:
            return None
        else:
            return self.taxa_trie.get_taxon(leaf.lineage.split(self.taxa_trie.lin_sep)[-1])

    def get_leaf_nodes_for_taxon(self, taxon_name: str) -> list:
        """
        Find all of the leaf nodes (TreeLeafReference objects) in a reference package that represent a taxon

        :param taxon_name: A rank-prefixed taxon name, such as 'd__Bacteria'
        :return: A list of TreeLeafReference instances with taxon_name in their lineages
        """
        taxon_leaf_nodes = []
        tree_leaves = sorted(self.generate_tree_leaf_references_from_refpkg(), key=lambda x: x.lineage)
        for leaf in tree_leaves:  # type: TreeLeafReference
            if not leaf.lineage:
                continue
            if taxon_name in leaf.lineage.split(self.taxa_trie.lin_sep):
                taxon_leaf_nodes.append(leaf)
        return taxon_leaf_nodes

    def get_leaf_node_by_name(self, leaf_name: str) -> list:
        leaves = []
        for leaf_node in self.generate_tree_leaf_references_from_refpkg():  # type: TreeLeafReference
            if leaf_name == leaf_node.description:
                leaves.append(leaf_node)
            desc, acc = leaf_node.description.split(" | ")
            if leaf_name == desc or leaf_name == acc:
                leaves.append(leaf_node)

        if len(leaves) == 0:
            logging.warning("Unable to find leaf '{}' in {} reference package leaves.\n".format(leaf_name, self.prefix))

        return leaves

    def map_rank_representatives_to_leaves(self, rank_name: str) -> (dict, dict):
        """
        Loads unique taxa found in the tree into a dictionary where the values are all leaves of a taxon.

        :param rank_name: Name of the rank to search for taxa
        :return:
          1. A dictionary mapping unique taxa of rank 'rank_name' mapped to lists of leaf node names
           (format = leaf.number_self.prefix)
          2. A dictionary mapping leaf.numbers to Taxon instances
        """
        # Check rank name to see if its present in the hierarchy
        if rank_name not in self.taxa_trie.accepted_ranks_depths:
            self.taxa_trie.so_long_and_thanks_for_all_the_fish("Rank '{}' is not in hierarchy's accepted names.\n"
                                                               "".format(rank_name))

        taxon_leaf_map = {}
        unique_taxa = {}

        rank_representatives = self.taxa_trie.rank_representatives(rank_name, with_prefix=True)

        for taxon_name in rank_representatives:
            taxon_leaf_nodes = self.get_leaf_nodes_for_taxon(taxon_name)
            taxon_leaf_map[taxon_name] = [leaf.number + '_' + self.prefix for leaf in taxon_leaf_nodes]

            for leaf in taxon_leaf_nodes:  # type: TreeLeafReference
                unique_taxa[leaf.number] = self.taxa_trie.hierarchy[taxon_name]

        logging.info("{}: {} unique taxa.\n".format(self.prefix, len(taxon_leaf_map)))

        return taxon_leaf_map, unique_taxa

    @staticmethod
    def get_unrelated_taxa(ref_leaf_nodes: dict, target_lineage: str, lin_sep="; ") -> set:
        """
        Returns a set of all ref_leaf_nodes keys whose corresponding lineages (value) are unrelated to target_lineage
        If a lineage (value) is truncated and is related to the target_lineage up until its last rank it is not
        included in the returned set (off_target_ref_leaves).
        This is because the resolution is too poor to know whether the ref_leaf is related to the target_lineage or not.

        :param ref_leaf_nodes: A dictionary of str keys mapped to taxonomic lineages
        :param target_lineage: A taxonomic lineage whose relatives are meant to be excluded
        :param lin_sep: An optional string that is used to separate the taxonomic ranks in the lineage
        :return: A set of all key names whose corresponding values are unrelated to the target_lineage
        """
        n_match = 0
        n_shallow = 0
        n_unclassified = 0
        off_target_ref_leaves = set()

        s_target = target_lineage.split(lin_sep)

        for ref_leaf_num, lineage in ref_leaf_nodes.items():  # type: (str, str)
            s_query = lineage.split(lin_sep)
            x = 0
            while x < min([len(s_query), len(s_target)]):
                if re.search("unclassified|environmental sample", s_query[x], re.IGNORECASE):
                    n_unclassified += 1
                    break
                if s_query[x] != s_target[x]:
                    off_target_ref_leaves.add(ref_leaf_num)
                    break
                x += 1
            if ref_leaf_num not in off_target_ref_leaves:
                if len(s_query) < len(s_target):
                    n_shallow += 1
                else:
                    n_match += 1

        logging.debug("Reference sequence filtering stats for taxon '{}'\n".format(s_target[-1]) +
                      "\n".join(["Match taxon\t" + str(n_match),
                                 "Unclassified\t" + str(n_unclassified),
                                 "Too shallow\t" + str(n_shallow),
                                 "Remaining\t" + str(len(off_target_ref_leaves))]) + "\n")
        return off_target_ref_leaves

    def get_monophyletic_clades(self, taxon_name: str, leaf_nodes=None) -> list:
        taxa_internal_node_map = []
        if not leaf_nodes:
            # Find all leaf nodes that represent taxon_name (i.e. p__Proteobacteria)
            leaf_nodes = set([leaf.number + '_' + self.prefix for leaf in self.get_leaf_nodes_for_taxon(taxon_name)])

        # Determine the deepest internal node that contains a monophyletic set of leaf_nodes
        # Sorting in descending order of number of child nodes
        internal_node_map = self.get_internal_node_leaf_map()
        for i_node in sorted(internal_node_map, key=lambda n: len(internal_node_map[n]), reverse=True):
            i_node_leaves = set(internal_node_map[i_node])
            if len(i_node_leaves) > len(leaf_nodes):
                continue
            # taxon_leaves can either be a superset or match
            # Detect whether the internal node contains a set of only the leaves belonging to the taxon
            if leaf_nodes.issuperset(i_node_leaves):
                taxa_internal_node_map.append(i_node_leaves)
                # Pop from the list of all taxon leaves
                for leaf in i_node_leaves:
                    leaf_nodes.remove(leaf)
            else:
                pass
            if len(leaf_nodes) == 0:
                break

        if leaf_nodes:
            logging.error("Unable to find an internal node containing the leaf nodes {} in reference package '{}'.\n"
                          "".format(', '.join(leaf_nodes), self.prefix))
            sys.exit(13)

        return taxa_internal_node_map

    def remove_taxon_from_lineage_ids(self, taxon: str) -> None:
        """
        Removes all sequences/leaves from the reference package that match the target taxon. Leaves with that have a
        taxonomic resolution lower than the target are also removed as their taxonomic provenance is uncertain.

        :param taxon: A '; '-separated taxonomic lineage for which all matches and descendents are removed
        :return: Nothing
        """
        leaf_lineage_map = {}
        off_target_lineage_ids = {}

        # Make the leaf name: lineage map for get_unrelated_taxa
        ref_leaf_nodes = self.generate_tree_leaf_references_from_refpkg()
        for ref_leaf in ref_leaf_nodes:  # type: TreeLeafReference
            leaf_lineage_map[ref_leaf.number] = ref_leaf.lineage

        # Get the ref_leaf.number for all TreeLeafReference instances that are unrelated to taxon
        off_target_ref_leaves = self.get_unrelated_taxa(ref_leaf_nodes=leaf_lineage_map,
                                                        target_lineage=taxon,
                                                        lin_sep=self.taxa_trie.lin_sep)

        # Format the lineage_ids strings
        for ref_leaf in ref_leaf_nodes:  # type: TreeLeafReference
            if ref_leaf.number in off_target_ref_leaves:
                off_target_lineage_ids[ref_leaf.number] = "{0} | {1}\t{2}".format(ref_leaf.description,
                                                                                  ref_leaf.accession,
                                                                                  ref_leaf.lineage)

        # Update self.lineage_ids with the remaining reference leaf nodes
        self.lineage_ids = off_target_lineage_ids
        self.num_seqs = len(self.lineage_ids)

        return

    def infer_phylogeny(self, input_msa: str, executables: dict, phylogeny_dir: str, num_threads=2) -> str:
        """
        Selects the substitution model and parameters for the phylogenetic inference tools then
        builds a phylogeny using the software specified by the ReferencePackage's tree_tool

        :param input_msa: A multiple sequence alignment (MSA) compatible with any of the supported tree building tools
        :param executables: A dictionary of executable names mapped to the
        :param phylogeny_dir: A directory for writing the outputs of this workflow
        :param num_threads: Number of threads to use while building the tree and bootstrapping
        :return: Path to the final phylogeny with greatest likelihood
        """
        best_tree = wrapper.construct_tree(self.tree_tool, executables, self.sub_model, input_msa,
                                           phylogeny_dir, self.prefix, num_threads)

        if self.tree_tool == "FastTree":
            etree = Tree(best_tree)
            farthest_node, farthest_dist = etree.get_farthest_node()
            etree.set_outgroup(farthest_node)
            etree.resolve_polytomy(recursive=True)  # Remove any polytomies (where a node divides into >two edges)
            etree.write(outfile=best_tree, format=5)

        return best_tree

    def recover_raxmlng_model_outputs(self, phylogeny_dir: str) -> None:
        # Find the best model file
        model_info = match_file(phylogeny_dir + "*.bestModel")

        # Import the tree and model info files into the reference package
        copy(model_info, self.f__model_info)

        return

    def recover_raxmlng_tree_outputs(self, phylogeny_dir: str) -> None:
        # Find the best tree file
        raw_newick_tree = match_file(phylogeny_dir + "*.bestTree")

        # Import the tree files into the reference package
        copy(raw_newick_tree, self.f__tree)
        return

    def recover_raxmlng_supports(self, phylogeny_dir: str) -> None:
        # Annotate the bootstrapped phylogeny
        bootstrap_tree = match_file(phylogeny_dir + "*.raxml.support")
        annotate_partition_tree(self.prefix, self.generate_tree_leaf_references_from_refpkg(), bootstrap_tree)
        copy(bootstrap_tree, self.f__boot_tree)
        return

    def exclude_clade_from_ref_files(self, tmp_dir: str, target_clade: str, executables: dict,
                                     fresh=False) -> None:
        """
        Removes all reference sequences/leaf nodes from a reference package that are descendents of a target clade.
        All reference package files are regenerated without these sequences and written to the original reference
        package directory.

        Original reference package files (so, still containing descendent sequences) are saved to a specified location.

        :param tmp_dir: Path to the treesapp/data (or reference package) directory for RAxML-NG to write
         output files from its `--evaluate` routine.
        :param target_clade: Taxonomic lineage of the clade that is being excluded from the reference package.
        :param executables: Dictionary of paths to dependency executables indexed by their names. Must include:
         'hmmbuild', 'FastTree' and 'raxml-ng'.
        :param fresh: Boolean indicating whether the reference package's tree should be built from scratch (True) or
         if the clades that are descendents of 'target_clade' should just be pruned (False) by ETE3
        :return: None
        """
        if tmp_dir[-1] != os.sep:
            tmp_dir = tmp_dir + os.sep
        self.disband(tmp_dir)

        split_lineage = target_clade.split(self.taxa_trie.lin_sep)
        base_taxon = self.taxa_trie.get_taxon(split_lineage[0])
        if not base_taxon:
            logging.error("Unable to find taxon '{}' in taxonomic hierarchy."
                          " Unable to exclude from reference package '{}'.\n".format(split_lineage[0], self.prefix))
            sys.exit(11)

        # tax_ids
        self.remove_taxon_from_lineage_ids(target_clade)

        # fasta
        ref_fasta_dict = read_fasta_to_dict(self.f__msa)
        if len(ref_fasta_dict) == 0:
            logging.error("No sequences were read from Reference Package {}'s FASTA '{}'.\n"
                          "".format(self.prefix, self.f__msa))
            sys.exit(13)
        off_target_ref_headers = [ref_num + '_' + self.prefix for ref_num in self.lineage_ids]
        if len(off_target_ref_headers) == 0:
            logging.error("No reference sequences were retained for building testing " + target_clade + "\n")
            sys.exit(19)
        split_files = write_new_fasta(ref_fasta_dict, self.f__msa, headers=off_target_ref_headers)
        if len(split_files) > 1:
            logging.error("Only one FASTA file should have been written.\n")
            sys.exit(21)

        # profile HMM
        wrapper.build_hmm_profile(executables["hmmbuild"],
                                  output_hmm=self.f__profile, msa_in=self.f__msa, graceful=True)

        # profile HMM used for homology search
        self.dereplicate_hmm(dereplication_rank="genus",
                             hmmbuild_exe=executables["hmmbuild"], mafft_exe=executables["mafft"],
                             realign=False)

        # Trees
        if fresh:
            tree_build_cmd = [executables["FastTree"]]
            if self.molecule == "rrna" or self.molecule == "dna":
                tree_build_cmd += ["-nt", "-gtr"]
            else:
                tree_build_cmd += ["-lg", "-gamma"]
            tree_build_cmd += ["-out", self.f__tree]
            tree_build_cmd.append(self.f__msa)
            logging.info("Building Approximately-Maximum-Likelihood tree with FastTree... ")
            stdout, returncode = launch_write_command(tree_build_cmd, True)
            with open(tmp_dir + os.sep + "FastTree_info." + self.prefix, 'w') as fast_info:
                fast_info.write(stdout + "\n")
            logging.info("done.\n")
        else:
            ref_tree = self.get_ete_tree()
            ref_tree.prune(off_target_ref_headers, preserve_branch_length=True)
            logging.debug("\t" + str(len(ref_tree.get_leaves())) + " leaves in pruned tree.\n")
            ref_tree.write(outfile=self.f__tree, format=5)

        # Model parameters
        model_output_prefix = os.path.join(tmp_dir, "tree_data")
        wrapper.model_parameters(executables["raxml-ng"],
                                 self.f__msa, self.f__tree, model_output_prefix, self.sub_model)
        self.recover_raxmlng_model_outputs(model_output_prefix)

        self.band()

        return

    def dereplicate_hmm(self, dereplication_rank: str, hmmbuild_exe, mafft_exe,
                        n_threads=2, intermediates_dir=None, realign=True) -> None:
        """
        Function to create a taxonomically-dereplicated hidden Markov model (HMM) profile. This reduces the bias from
        potentially over-represented clades, increasing the weight of their conserved positions.

        :param dereplication_rank: The taxonomic rank to dereplicate the reference sequences at
        :param hmmbuild_exe: Path to an hmmbuild executable
        :param mafft_exe: Path to a MAFFT executable
        :param n_threads: Number of threads MAFFT should use
        :param intermediates_dir: A path to a directory to write intermediate files
        :param realign: Boolean indicating whether the sequences with a lineage excluded should be realigned using MAFFT
        :return: None
        """

        logging.debug("Creating taxonomically-dereplicated HMM... ")

        if not intermediates_dir:
            intermediates_dir = os.path.dirname(self.f__msa)
        if intermediates_dir[-1] != os.sep:
            intermediates_dir += os.sep
        derep_fa = intermediates_dir + base_file_prefix(self.f__msa) + "_derep.fa"
        derep_aln = intermediates_dir + base_file_prefix(self.f__msa) + "_derep.mfa"
        intermediates = [derep_aln, derep_fa]

        lineage_reps = []
        t = {}

        mfa = FASTA(self.f__msa)
        mfa.load_fasta()

        # Trim the taxonomic lineages to the dereplication level
        leaf_taxa_map = {leaf.number + "_" + self.prefix: leaf.lineage for leaf in
                         self.generate_tree_leaf_references_from_refpkg()}
        trimmed_lineages = self.taxa_trie.trim_lineages_to_rank(leaf_taxa_map, dereplication_rank)
        # Add back the sequences with truncated lineages
        for leaf_name in leaf_taxa_map:
            if leaf_name not in trimmed_lineages:
                trimmed_lineages[leaf_name] = leaf_taxa_map[leaf_name]

        # Find the longest sequence that each lineage
        for seq_name in mfa.get_seq_names():
            taxon = trimmed_lineages[seq_name]
            if taxon in t:
                curr_rep = t[taxon]
                if len(mfa.fasta_dict[curr_rep]) > len(mfa.fasta_dict[seq_name]):
                    continue
            t[taxon] = seq_name
        for taxon in t:
            lineage_reps.append(t[taxon])

        # Remove all sequences from the FASTA instance that are not representatives
        mfa.keep_only(lineage_reps)
        if realign:
            mfa.unalign()
            # Write the dereplicated FASTA file
            write_new_fasta(fasta_dict=mfa.fasta_dict, fasta_name=derep_fa)
            # Re-align the sequences
            wrapper.run_mafft(mafft_exe=mafft_exe, fasta_in=derep_fa, fasta_out=derep_aln, num_threads=n_threads)
        else:
            write_new_fasta(fasta_dict=mfa.fasta_dict, fasta_name=derep_aln)

        # Build the new HMM profile
        wrapper.build_hmm_profile(hmmbuild_exe=hmmbuild_exe, msa_in=derep_aln,
                                  output_hmm=self.f__search_profile, name=self.prefix)

        # Clean up intermediates
        for f_path in intermediates:
            if os.path.isfile(f_path):
                os.remove(f_path)

        logging.debug("done.\n")

        logging.debug("%i %s-dereplicated sequences retained for building HMM profile.\n" %
                      (len(lineage_reps), dereplication_rank))
        return

    def load_pfit_params(self, build_param_line):
        """
        Loads the polynomial classifier parameters that are stored in the ReferencePackage.pfit attribute
        """
        build_param_fields = build_param_line.split("\t")
        if build_param_fields[8]:
            self.pfit = [float(x) for x in build_param_fields[8].split(',')]
        return

    def enumerate_taxonomic_lineages(self) -> dict:
        """
        At times, it is handy to know how many of each taxonomic label are present in the reference package.
        enumerate_taxonomic_lineages() iterates over the ranks (as the numeric representation) and sums the number of
        taxa that are descendents of each broader taxonomic rank.
        This leads to multiple-counts of each taxon by design.

        :return: A dictionary mapping unique taxa in the ReferencePackage to their respective count
        """
        rank = 1
        taxonomic_counts = dict()
        while rank < 9:
            for leaf in self.generate_tree_leaf_references_from_refpkg():
                lineage = leaf.lineage.split(self.taxa_trie.lin_sep)
                if len(lineage) < rank:
                    continue
                taxonomy = "; ".join(lineage[:rank])
                if taxonomy not in taxonomic_counts:
                    taxonomic_counts[taxonomy] = 0
                taxonomic_counts[taxonomy] += 1
            rank += 1
        return taxonomic_counts

    def all_possible_assignments(self):
        """
        Returns a StringTrie (prefix-tree) representation of the taxonomic lineages of the ReferencePackage
        """
        if len(self.lineage_ids) == 0:
            logging.error("ReferencePackage.lineage_ids is empty - information hasn't been slurped up yet.\n")
            sys.exit(17)

        lineage_list = list()
        for ref_leaf_node in self.generate_tree_leaf_references_from_refpkg():  # type: TreeLeafReference
            lineage = ref_leaf_node.lineage
            if not re.match(r"^r__Root.*", lineage):
                lineage = "r__Root" + self.taxa_trie.lin_sep + lineage
            lineage_list.append(lineage)

        return load_taxonomic_trie(lineage_list)

    def taxonomically_label_tree(self) -> Tree:
        """
        When deciding what the taxonomic label should be assigned to a query sequence, algorithms may require a tree
        (ete3.Tree instance) for which each node contains a 'taxon' instance, itself being a taxonomic_hierarchy.Taxon
        instance.
        This function assures that the phylogeny is strictly bifurcating, by using ete3's resolve_polytomy() function.
        To label the nodes of a phylogeny,
        """
        # Generate an instance of the ETE3 Master Tree class
        rt = self.get_ete_tree()

        # Ensure that every taxonomic lineage is rooted by setting the parent of each 'domain' taxon to "r__Root"
        root_taxon = self.taxa_trie.find_root_taxon()
        self.taxa_trie.root_domains(root_taxon)

        # Propagate a 'taxon' feature - None by default - to all TreeNodes for holding Taxon instances
        for n in rt.traverse(strategy="postorder"):  # type: Tree
            n.add_feature(pr_name="taxon", pr_value=None)

        # Add the taxonomic lineages as features to each of the leaf nodes
        for leaf_node in self.generate_tree_leaf_references_from_refpkg():  # type: TreeLeafReference
            taxon = self.get_hierarchy_taxon_for_leaf(leaf_node)
            tree_node = rt.get_leaves_by_name(name="{}_{}".format(leaf_node.number, self.prefix)).pop()
            tree_node.add_feature(pr_name="taxon", pr_value=taxon)

        # Label the internal nodes based on the LCA of the leaf nodes
        for n in rt.traverse(strategy="postorder"):
            if len(n.children) == 0:  # this is a leaf node
                continue
            elif len(n.children) == 1:
                self.bail("Unexpected number of children (1) for Tree.Node {}.\n".format(n.name))
                raise AssertionError
            elif len(n.children) > 2:
                n.resolve_polytomy(recursive=False)
                for c in n.children:
                    if not hasattr(c, "taxon"):
                        c.add_feature(pr_name="taxon", pr_value=None)
                        c.taxon = Taxon.lca(c.children[0].taxon, c.children[1].taxon)
            l_node, r_node = n.get_children()
            lca = Taxon.lca(l_node.taxon, r_node.taxon)
            n.taxon = lca

        return rt


def view(refpkg: ReferencePackage, attributes: list) -> None:
    view_dict = {}
    for attr in attributes:
        try:
            view_dict[attr] = refpkg.__dict__[attr]
        except KeyError:
            logging.error("Attribute '{}' doesn't exist in ReferencePackage.\n".format(attr))
            sys.exit(1)

    for k, v in view_dict.items():
        if type(v) is list and k not in ["pfit"]:
            v = ''.join(v)
        if type(v) is dict:
            v = "\n" + "\n".join([sk + "\t" + sv for sk, sv in v.items()])
        logging.info("{}\t{}\n".format(k, v))

    # TODO: optionally use ReferencePackage.write_refpkg_component
    return


def edit(refpkg: ReferencePackage, attributes: list, output_dir, overwrite: bool) -> None:
    if len(attributes) > 2:
        logging.error("`treesapp package edit` only edits a single attribute at a time.\n")
        sys.exit(3)
    elif len(attributes) == 1:
        logging.error("`treesapp package edit` requires a value to change.\n")
        sys.exit(3)
    else:
        k, v = attributes  # type: (str, str)

    try:
        current_v = refpkg.__dict__[k]
    except KeyError:
        logging.error("Attribute '{}' doesn't exist in ReferencePackage.\n".format(k))
        sys.exit(1)

    logging.info("Replacing attribute '{}' (currently '{}')\n".format(k, current_v))

    if k.startswith("f__"):  # This is a file that is being modified and the contents need to be read
        k_content = re.sub(r"^f__", '', k)
        if not os.path.isfile(v):
            logging.error("Unable to find path to new file '{}'. "
                          "Exiting and the reference package '{}' will remain unchanged.\n".format(v, refpkg.f__json))
            sys.exit(5)
        with open(v) as content_handler:
            v_content = content_handler.readlines()
        try:
            refpkg.__dict__[k_content] = v_content
        except KeyError:
            logging.error("Unable to find ReferencePackage attribute '{}' in '{}'.\n".format(k_content, refpkg.f__json))
            sys.exit(7)
    refpkg.__dict__[k] = v
    if not overwrite:
        refpkg.f__json = os.path.join(output_dir, os.path.basename(refpkg.f__json))
        if os.path.isfile(refpkg.f__json):
            logging.warning("RefPkg file '{}' already exists.\n".format(refpkg.f__json))
            return
    refpkg.pickle_package()

    return
