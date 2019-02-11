__author__ = 'Connor Morgan-Lang'

import sys
import os
import re
import copy
import subprocess
import logging
from multiprocessing import Process, JoinableQueue
from json import loads, dumps
from fasta import get_header_format
from utilities import reformat_string, return_sequence_info_groups, median
from entish import get_node, create_tree_info_hash, subtrees_to_dictionary

import _tree_parser


class ReferencePackage:
    def __init__(self):
        self.prefix = ""
        self.msa = ""
        self.profile = ""
        self.tree = ""
        self.boot_tree = ""
        self.lineage_ids = ""
        self.sub_model = ""
        self.core_ref_files = list()
        self.num_seqs = 0

    def validate(self, num_ref_seqs=None):
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
        # TODO: Compare the number of sequences in the Hidden-Markov model
        # TODO: Compare the number of sequences in the Tree files
        # TODO: Compare the number of sequences in the tax_ids file
        return True

    def gather_package_files(self, ref_name: str, pkg_path: str, pkg_format="hierarchical", molecule="prot"):
        """
        Populates a ReferencePackage instances fields with files based on 'pkg_format' where hierarchical indicates
         files are sorted into 'alignment_data', 'hmm_data' and 'tree_data' directories and flat indicates they are all
         in the same directory.
        :param ref_name: Prefix name of all the files of a reference package
        :param pkg_path: Path to the reference package
        :param pkg_format: The format of the files within the pkg_path directory
        :param molecule: A string indicating the molecule type of the reference package. If 'rRNA' profile is CM.
        :return:
        """
        self.prefix = ref_name
        if pkg_format == "flat":
            self.msa = pkg_path + os.sep + ref_name + ".fa"
            self.profile = pkg_path + os.sep + ref_name
            self.tree = pkg_path + os.sep + ref_name + "_tree.txt"
            self.boot_tree = pkg_path + os.sep + ref_name + "_bipartitions.txt"
            self.lineage_ids = pkg_path + os.sep + "tax_ids_" + ref_name + ".txt"
        elif pkg_format == "hierarchical":
            self.msa = pkg_path + os.sep + "alignment_data" + os.sep + ref_name + ".fa"
            self.profile = pkg_path + os.sep + "hmm_data" + os.sep + ref_name
            self.tree = pkg_path + os.sep + "tree_data" + os.sep + ref_name + "_tree.txt"
            self.boot_tree = pkg_path + os.sep + "tree_data" + os.sep + ref_name + "_bipartitions.txt"
            self.lineage_ids = pkg_path + os.sep + "tree_data" + os.sep + "tax_ids_" + ref_name + ".txt"
        else:
            logging.error("Unrecognised reference package format '" + pkg_format + "'\n")
            sys.exit(17)

        if molecule == "rRNA":
            self.profile += ".cm"
        else:
            self.profile += ".hmm"

        self.core_ref_files += [self.msa, self.profile, self.tree, self.lineage_ids]

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

    def load_build_params(self, build_param_line):
        build_param_fields = build_param_line.split('\t')
        if len(build_param_fields) != 11:
            logging.error("Incorrect number of values (" + str(len(build_param_fields)) +
                          ") in ref_build_parameters.tsv. Line:\n" + build_param_line)
            sys.exit(17)

        self.cog = build_param_fields[0]
        self.denominator = build_param_fields[1]
        self.molecule = build_param_fields[2]
        self.model = build_param_fields[3]
        self.kind = build_param_fields[4]
        self.pid = build_param_fields[5]
        self.num_reps = build_param_fields[6]
        self.tree_tool = build_param_fields[7]
        self.lowest_confident_rank = build_param_fields[9]
        self.update = build_param_fields[10]
        self.description = build_param_fields[-1]

    def load_pfit_params(self, build_param_line):
        build_param_fields = build_param_line.split("\t")
        if build_param_fields[8]:
            self.pfit = [float(x) for x in build_param_fields[8].split(',')]
        return

    def check_rank(self):
        taxonomies = ["NA", "Kingdoms", "Phyla", "Classes", "Orders", "Families", "Genera", "Species"]

        if self.lowest_confident_rank not in list(taxonomies):
            logging.error("Unable to find " + self.lowest_confident_rank + " in taxonomic map!")
            sys.exit(17)

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
        summary_string += "\nInformation for query sequence '" + self.contig_name + "'\n"
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
            summary_string += "Distances:\n\t" + self.distances + "\n"
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

    def correct_decoding(self):
        """
        Since the JSON decoding is unable to decode recursively, this needs to be fixed for each placement
        Formatting and string conversion are also performed here

        :return:
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

    def rename_placed_sequence(self, seq_name):
        new_placement_collection = dict()
        for d_place in self.placements:
            for key, value in d_place.items():
                if key == 'n':
                    new_placement_collection['n'] = [seq_name]
                else:
                    new_placement_collection[key] = value
        self.placements = [new_placement_collection]
        return

    def name_placed_sequence(self):
        for d_place in self.placements:
            for key, value in d_place.items():
                if key == 'n':
                    self.contig_name = value[0]
        return

    def get_field_position_from_jplace_fields(self, field_name):
        """
        Find the position in self.fields of 'like_weight_ratio'
        :return: position in self.fields of 'like_weight_ratio'
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

    def get_jplace_element(self, element_name):
        """
        Determines the element value (e.g. likelihood, edge_num) for a single placement.
        There may be multiple placements (or 'pquery's) in a single .jplace file.
        Therefore, this function is usually looped over.
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

    def filter_min_weight_threshold(self, threshold=0.1):
        """
        Remove all placements with likelihood weight ratios less than threshold
        :param threshold: The threshold which all placements with LWRs less than this are removed
        :return:
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
            dict_strings = list()
            if len(placement["p"]) > 1:
                for k, v in placement.items():
                    if k == 'p':
                        # For debugging:
                        # sys.stderr.write(str(v) + "\nRemoved:\n")
                        acc = 0
                        tmp_placements = copy.deepcopy(v)
                        while acc < len(tmp_placements):
                            candidate = tmp_placements[acc]
                            if float(candidate[x]) < threshold:
                                removed = tmp_placements.pop(acc)
                                # For debugging:
                                # sys.stderr.write("\t".join([self.name, str(removed[0]),
                                #                             str(float(removed[x]))]) + "\n")
                            else:
                                acc += 1
                            # sys.stderr.flush()
                        # If no placements met the likelihood filter then the sequence cannot be classified
                        # Alternatively: first two will be returned and used for LCA - can test...
                        if len(tmp_placements) > 0:
                            v = tmp_placements
                            dict_strings.append(dumps(k) + ':' + dumps(v))
                            placement_string = ', '.join(dict_strings)
                        else:
                            self.classified = False
                # Add the filtered placements back to the object.placements
                new_placement_collection.append('{' + placement_string + '}')
            else:
                # If there is only one placement, the LWR is 1.0 so no filtering required!
                new_placement_collection.append(pquery)
        if self.classified:
            self.placements = new_placement_collection
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
                        normalized_abundance = float(self.abundance/len(tree_leaves))
                        for tree_leaf in tree_leaves:
                            if tree_leaf not in leaf_rpkm_sums.keys():
                                leaf_rpkm_sums[tree_leaf] = 0.0
                            leaf_rpkm_sums[tree_leaf] += normalized_abundance
        return leaf_rpkm_sums

    def filter_max_weight_placement(self):
        """
        Removes all secondary placements of each pquery,
        leaving only the placement with maximum like_weight_ratio
        :return:
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
                            tmp_placements = copy.deepcopy(v)
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

    def create_jplace_node_map(self):
        """
        Loads a mapping between all nodes (internal and leaves) and all leaves
        :return:
        """
        no_length_tree = re.sub(":[0-9.]+{", ":{", self.tree)
        self.node_map.clear()
        node_stack = list()
        leaf_stack = list()
        x = 0
        num_buffer = ""
        while x < len(no_length_tree):
            c = no_length_tree[x]
            if re.search(r"[0-9]", c):
                while re.search(r"[0-9]", c):
                    num_buffer += c
                    x += 1
                    c = no_length_tree[x]
                node_stack.append([str(num_buffer)])
                num_buffer = ""
                x -= 1
            elif c == ':':
                # Append the most recent leaf
                current_node, x = get_node(no_length_tree, x + 1)
                self.node_map[current_node] = node_stack.pop()
                leaf_stack.append(current_node)
            elif c == ')':
                # Set the child leaves to the leaves of the current node's two children
                while c == ')' and x < len(no_length_tree):
                    if no_length_tree[x + 1] == ';':
                        break
                    current_node, x = get_node(no_length_tree, x + 2)
                    self.node_map[current_node] = self.node_map[leaf_stack.pop()] + self.node_map[leaf_stack.pop()]
                    leaf_stack.append(current_node)
                    x += 1
                    c = no_length_tree[x]
            x += 1
        return

    def harmonize_placements(self, treesapp_dir):
        """
        Often times, the placements field in a jplace file contains multiple possible tree locations.
        In order to consolidate these into a single tree location, the LCA algorithm is utilized.
        Since all placements are valid, there is no need to be uncertain about including all nodes during LCA compute
        :return: The single internal node which is the parent node of all possible placements is returned.
        """
        if self.name == "nr":
            self.name = "COGrRNA"
        reference_tree_file = os.sep.join([treesapp_dir, "data", "tree_data"]) + os.sep + self.name + "_tree.txt"
        reference_tree_elements = _tree_parser._read_the_reference_tree(reference_tree_file)
        lwr_pos = self.get_field_position_from_jplace_fields("like_weight_ratio")
        if not lwr_pos:
            return
        singular_placements = list()
        for pquery in self.placements:
            placement = loads(pquery, encoding="utf-8")
            dict_strings = list()
            for k, v in placement.items():
                if len(v) > 1:
                    lwr_sum = 0
                    loci = list()
                    for locus in v:
                        lwr_sum += float(locus[lwr_pos])
                        loci.append(str(self.node_map[locus[0]][0]))
                    ancestral_node = _tree_parser._lowest_common_ancestor(reference_tree_elements, ','.join(loci))
                    # Create a placement from the ancestor, and the first locus in loci fields
                    v = [[ancestral_node, v[0][1], round(lwr_sum, 2), 0, 0]]
                dict_strings.append(dumps(k) + ':' + dumps(v))
            singular_placements.append('{' + ','.join(dict_strings) + '}')

        self.placements = singular_placements
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


class TreeLeafReference:
    """
    Objects for each leaf in a tree
    """
    def __init__(self, number, description):
        # self.tree = ""
        self.number = number
        self.description = description
        self.lineage = ""
        self.complete = False

    def summarize_tree_leaf(self):
        summary_string = "Leaf ID:\n\t" + str(self.number) + "\n" +\
                         "Description:\n\t" + str(self.description) + "\n"

        # summary_string += "Tree:\n\t" + str(self.tree) + "\n"
        if self.complete:
            summary_string += "Lineage:\n\t" + str(self.lineage) + "\n"
        logging.debug(summary_string)

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


class ReferenceSequence:
    def __init__(self):
        self.accession = ""
        self.description = ""
        self.organism = ""
        self.lineage = ""
        self.short_id = ""
        self.sequence = ""
        self.locus = ""
        self.cluster_rep = False
        self.cluster_rep_similarity = 0
        self.cluster_lca = None

    def get_info(self):
        """
        Returns a string with the ReferenceSequence instance's current fields

        :return: str
        """
        info_string = ""
        info_string += "accession = " + self.accession + ", " + "mltree_id = " + self.short_id + "\n"
        info_string += "description = " + self.description + ", " + "locus = " + self.locus + "\n"
        info_string += "organism = " + self.organism + "\n"
        info_string += "lineage = " + self.lineage + "\n"
        return info_string


class EntrezRecord(ReferenceSequence):
    def __init__(self, acc, ver):
        super().__init__()
        self.accession = acc
        self.versioned = ver
        self.ncbi_tax = ""
        self.bitflag = 0  # For monitoring progress during download stage


class CommandLineWorker(Process):
    def __init__(self, task_queue, commander):
        Process.__init__(self)
        self.task_queue = task_queue
        self.master = commander

    def run(self):
        while True:
            next_task = self.task_queue.get()
            if next_task is None:
                # Poison pill means shutdown
                self.task_queue.task_done()
                break
            p_instance = subprocess.Popen(' '.join(next_task), shell=True, preexec_fn=os.setsid)
            p_instance.wait()
            if p_instance.returncode != 0:
                logging.error(self.master + " did not complete successfully for:\n" +
                              str(' '.join(next_task)) + "\n")
            self.task_queue.task_done()
        return


class CommandLineFarmer:
    """
    A worker that will launch command-line jobs using multiple processes in its queue
    """

    def __init__(self, command, num_threads):
        """
        Instantiate a CommandLineFarmer object to oversee multiprocessing of command-line jobs
        :param command:
        :param num_threads:
        """
        self.max_size = 32767  # The actual size limit of a JoinableQueue
        self.task_queue = JoinableQueue(self.max_size)
        self.num_threads = int(num_threads)

        genewise_process_queues = [CommandLineWorker(self.task_queue, command) for i in range(int(self.num_threads))]
        for process in genewise_process_queues:
            process.start()

    def add_tasks_to_queue(self, task_list):
        """
        Function for adding commands from task_list to task_queue while ensuring space in the JoinableQueue
        :param task_list: List of commands
        :return: Nothing
        """
        num_tasks = len(task_list)

        task = task_list.pop()
        while task:
            if not self.task_queue.full():
                self.task_queue.put(task)
                if num_tasks > 1:
                    task = task_list.pop()
                    num_tasks -= 1
                else:
                    task = None

        i = self.num_threads
        while i:
            if not self.task_queue.full():
                self.task_queue.put(None)
                i -= 1

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
            result = _tree_parser._build_subtrees_newick(next_task)
            subtrees = subtrees_to_dictionary(result, create_tree_info_hash())
            self.task_queue.task_done()
            self.result_queue.put(subtrees)
        return


class Header:
    def __init__(self, header):
        self.original = header
        self.formatted = ""
        self.treesapp_name = ""
        self.post_align = ""
        self.first_split = ""

    def get_info(self):
        info_string = "TreeSAPP ID = '" + self.treesapp_name + "'\tPrefix = '" + self.first_split + "'\n"
        info_string += "Original =  " + self.original + "\nFormatted = " + self.formatted
        return info_string


def register_headers(header_list):
    header_registry = dict()
    acc = 1
    for header in header_list:
        new_header = Header(header)
        new_header.formatted = reformat_string(header)
        new_header.first_split = header.split()[0]
        header_registry[str(acc)] = new_header
        acc += 1
    return header_registry


def get_header_info(header_registry, code_name=''):
    """

    :param header_registry: A dictionary of Header instances, indexed by numerical treesapp_id
    :param code_name: [OPTIONAL] The code_name of the reference package (marker gene/domain/family/protein)
    :return: Dictionary where keys are numerical treesapp_ids and values are ReferenceSequence instances
    """
    logging.info("Extracting information from headers... ")
    fasta_record_objects = dict()
    for treesapp_id in sorted(header_registry.keys(), key=int):
        original_header = header_registry[treesapp_id].original
        formatted_header = header_registry[treesapp_id].formatted
        header_format_re, header_db, header_molecule = get_header_format(original_header, code_name)
        sequence_info = header_format_re.match(original_header)
        accession, organism, locus, description, lineage = return_sequence_info_groups(sequence_info,
                                                                                       header_db,
                                                                                       formatted_header)
        ref_seq = ReferenceSequence()
        ref_seq.organism = organism
        ref_seq.accession = accession
        ref_seq.lineage = lineage
        ref_seq.description = description
        ref_seq.locus = locus
        ref_seq.short_id = '>' + treesapp_id + '_' + code_name
        fasta_record_objects[treesapp_id] = ref_seq

    logging.info("done.\n")

    return fasta_record_objects


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


def prep_logging(log_file_name, verbosity):
    logging.basicConfig(level=logging.DEBUG,
                        filename=log_file_name,
                        filemode='w',
                        datefmt="%d/%m %H:%M:%S",
                        format="%(asctime)s %(levelname)s:\n%(message)s")
    if verbosity:
        logging_level = logging.DEBUG
    else:
        logging_level = logging.INFO

    # Set the console handler normally writing to stdout/stderr
    ch = logging.StreamHandler()
    ch.setLevel(logging_level)
    ch.terminator = ''

    formatter = MyFormatter()
    ch.setFormatter(formatter)
    logging.getLogger('').addHandler(ch)

    return


class MarkerTest:
    def __init__(self, marker_name):
        self.target_marker = marker_name
        self.ranks = list()
        self.markers = set()
        self.taxa_filter = dict()
        self.taxa_filter["Unclassified"] = 0
        self.taxa_filter["Classified"] = 0
        self.taxa_filter["Unique_taxa"] = 0
        self.taxa_tests = dict()
        self.classifications = dict()

    def new_taxa_test(self, rank, lineage):
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
            for tt in self.taxa_tests[rank]:
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
                                           str(median(distals))]),
                                "\t".join([rank, "Pendant",
                                           str(round(sum(pendants) / float(n_dists), 4)),
                                           str(median(pendants))]),
                                "\t".join([rank, "Tip",
                                           str(round(sum(tips) / float(n_dists), 4)),
                                           str(median(tips))]),
                                "\t".join([rank, "Total",
                                           str(round(sum(totals) / float(n_dists), 4)),
                                           str(median(totals))], )]
            sys.stdout.write("\n".join(distance_summary) + "\n")
            return distals, pendants, tips
        else:
            return None, None, None


class TaxonTest:
    def __init__(self, name):
        self.lineage = name
        self.taxon = name.split('; ')[-1]
        self.queries = list()
        self.classifieds = list()
        self.distances = dict()
        self.assignments = dict()
        self.taxonomic_tree = None

    def summarise_taxon_test(self):
        summary_string = "Test for taxonomic lineage '" + self.lineage + "':\n" + \
                         "\tNumber of query sequences = " + str(len(self.queries)) + "\n" + \
                         "\tNumber of classified queries = " + str(len(self.classifieds)) + "\n"
        if self.assignments:
            for marker in self.assignments:
                summary_string += "Sequences classified as marker '" + marker + "':\n"
                for lineage in self.assignments[marker]:
                    summary_string += str(len(self.assignments[marker][lineage])) + "\t'" + lineage + "'\n"
        return summary_string

    def filter_assignments(self, target_marker):
        """
        Filters the assignments from TreeSAPP for the target marker.
        Off-target classifications are accounted for and reported.
        TaxonTest.classifieds only include the headers of the correctly annotated sequences
        :param target_marker:
        :return:
        """
        off_targets = dict()
        for marker in self.assignments:
            for lineage in self.assignments[marker]:
                classifieds = self.assignments[marker][lineage]
                if marker == target_marker:
                    self.classifieds += classifieds
                else:
                    if marker not in off_targets:
                        off_targets[marker] = list()
                    off_targets[marker] += classifieds
        if off_targets:
            for marker in off_targets:
                logging.warning(str(len(off_targets)) + " sequences were classified as " + marker + ":\n" +
                                "\t\n".join(off_targets[marker]) + "\n")
        return
