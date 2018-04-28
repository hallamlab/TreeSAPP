__author__ = 'Connor Morgan-Lang'

import sys
import os
import shutil
import re
import random
import copy
import subprocess
from multiprocessing import Process, JoinableQueue
from json import loads, dumps

from fasta import format_read_fasta, get_headers, write_new_fasta, get_header_format
from utilities import reformat_string, get_lineage
from entish import get_node, create_tree_info_hash, subtrees_to_dictionary

import _tree_parser


class CreateFuncTreeUtility:
    """
    Output is the directory to write the outputs for the updated tree
    InputData is the path to the TreeSAPP output folder containing, various_outputs/ and final_RAxML_outputs/
    RefTree is the second column in cog_list.tsv for the gene to update
    Cluster is a flag indicating whether the protein sequences for the RefTree in InputData is to be clustered at 97%
    """
    def __init__(self, input_data, ref_tree):
        if os.path.isabs(input_data):
            self.InputData = input_data
        else:
            self.InputData = os.getcwd() + os.sep + input_data

        if self.InputData[-1] == '/':
            self.InputData = self.InputData[:-1]

        self.Output = self.InputData + os.sep + "updated_" + ref_tree + "_tree" + os.sep
        self.Denominator = ref_tree
        self.marker_molecule = ""  # prot, dna, or rrna
        self.COG = ""
        self.raxml_model = ""  # Rather than determining the best model again, use the one that was previously used
        self.ContigDict = dict()  # Used for storing the final candidate reference sequences
        self.names = list()
        self.header_id_map = dict()  # Used for mapping the original header of the sequence to internal numeric header
        self.cluster_id = 0  # The percent similarity the original reference sequences were clustered at
        # Automatically remove the last attempt at updating the reference tree
        if os.path.isdir(self.Output):
            shutil.rmtree(self.Output)
        try:
            os.makedirs(self.Output)
        except IOError:
            raise IOError("Unable to make the directory " + str(self.Output) + "\n")

    def get_raxml_files_for_ref(self):
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
        if not self.COG:
            sys.stderr.write("ERROR: Invalid marker code provided! Unable to find matching gene in cog_list.txt!\n")
            sys.exit(100)
        return

    def find_marker_type(self, cog_list):
        for marker_type in cog_list:
            if marker_type == "all_cogs":
                continue
            if self.COG in cog_list[marker_type]:
                if marker_type == "phylogenetic_rRNA_cogs":
                    self.marker_molecule = "rrna"
                elif marker_type == "functional_cogs":
                    self.marker_molecule = "prot"
                else:
                    sys.stderr.write("ERROR: Unrecognized marker class: " + marker_type)
                    sys.exit(101)
                break

        return

    def load_new_refs_fasta(self, args, centroids_fasta, ref_organism_lineage_info):
        # Determine whether the numerical IDs are a series and sorted
        acc = 0
        ref_numeric_ids = list()
        for leaf in ref_organism_lineage_info[self.Denominator]:
            ref_numeric_ids.append(int(leaf.number))
            acc += 1

        # Read the FASTA to get headers and sequences
        centroids_fasta_dict = format_read_fasta(centroids_fasta, self.marker_molecule, args)

        # Create the final contig dictionary of new internal TreeSAPP headers (keys) and sequences (values)
        additional = acc
        for new_ref_seq in centroids_fasta_dict:
            if acc == sorted(ref_numeric_ids, key=int)[-1]:
                additional += 1
                # These sequences form a series, therefore continue the series for new reference sequences
                internal_id = '>' + str(additional) + '_' + self.COG
            else:
                rfive = ''.join(str(x) for x in random.sample(range(10), 5))
                while rfive in ref_numeric_ids:
                    rfive = ''.join(str(x) for x in random.sample(range(10), 5))
                internal_id = '>' + rfive + '_' + self.COG
            # Map the new reference headers to their numerical IDs
            self.header_id_map[internal_id] = new_ref_seq
            self.ContigDict[internal_id] = centroids_fasta_dict[new_ref_seq]

        if additional == acc:
            # The reference sequence identifiers are random
            sys.stderr.write("WARNING: numerical TreeSAPP identifiers for " + self.Denominator +
                             "are not in the format of a sequential series!\n")
            sys.stderr.write("Generating random numerical unique identifiers for the new sequence(s).\n")
            sys.stderr.flush()

        if args.verbose:
            sys.stdout.write("\t" + str(len(self.header_id_map)) + " new " + self.COG + " reference sequences.\n")

        return

    def update_tax_ids(self, args, ref_organism_lineage_info):
        """
        Write the number, organism and accession ID, if possible
        :param args:
        :param ref_organism_lineage_info:
        :return:
        """
        sys.stdout.write("Writing updated tax_ids file... ")
        sys.stdout.flush()

        tree_taxa_string = ""
        original_to_formatted_header_map = dict()

        if self.Denominator not in ref_organism_lineage_info.keys():
            raise ValueError(self.Denominator + " not included in data from tax_ids files!\n")

        # Load the original reference sequences first since this shouldn't change
        for leaf in ref_organism_lineage_info[self.Denominator]:
            if leaf.lineage:
                tree_taxa_string += '\t'.join([str(leaf.number), leaf.description, leaf.lineage]) + "\n"
            else:
                sys.stderr.write("Unable to retrieve lineage information for sequence " + str(leaf.number) + "\n")

        # Build a map of original headers to formatted ones
        original_candidate_headers = get_headers(args.fasta_input)
        for original in original_candidate_headers:
            original_to_formatted_header_map[reformat_string(original)] = original

        # Now figure out what to do with the new reference sequences
        for reference in sorted(self.ContigDict.keys()):
            num_id = reference[1:].split('_')[0]
            reformatted_header = self.header_id_map[reference]
            header = original_to_formatted_header_map[reformatted_header]
            # Check to see if this header contains an accession value or if it is a contig name with no useful info
            header_format_re, header_db, header_molecule = get_header_format(header, self.Denominator)
            sequence_info = header_format_re.match(header)
            if sequence_info:
                candidate_acc = sequence_info.group(1)
                lineage = get_lineage(candidate_acc, header_molecule)
                try:
                    taxonomic_lineage = lineage.split('; ')
                    # Try to get the species name
                    if len(taxonomic_lineage) >= 8:
                        description = ' '.join(taxonomic_lineage[6:8])
                    # We'll settle for the Phylum
                    elif len(taxonomic_lineage) > 3:
                        description = taxonomic_lineage[3]
                    else:
                        description = taxonomic_lineage[-1]
                except ValueError:
                    sys.stderr.write("WARNING: Attempt to parse species from lineage failed for:\n" + lineage + "\n")
                    description = sequence_info.group(2)
                description = description + " | " + candidate_acc
            else:
                # This sequence is probably an uninformative contig name
                lineage = ";"
                description = self.InputData.split('/')[-1]

            tree_taxa_string += num_id + "\t" + description + "\t" + lineage + "\n"

        # Write the new TreeSAPP numerical IDs, descriptions and lineages
        tree_taxa_list = self.Output + "tax_ids_" + self.COG + ".txt"
        tree_tax_list_handle = open(tree_taxa_list, "w")
        tree_tax_list_handle.write(tree_taxa_string)
        tree_tax_list_handle.close()

        sys.stdout.write("done.\n")
        sys.stdout.flush()

        return

    def swap_tree_names(self, tree, tree_swap_name):
        """
        Function used for replacing unique identifiers in a NEWICK tree file
        :param tree: The tree with leaf names that need to be replaced
        :param tree_swap_name: Name of the output tree file
        :return:
        """
        try:
            old_tree = open(tree, 'r')
        except IOError:
            raise IOError("Unable to open " + tree + " for reading!")
        try:
            new_tree = open(tree_swap_name, 'w')
        except IOError:
            raise IOError("Unable to open " + tree_swap_name + " for writing!")

        newick_tree = old_tree.readlines()
        old_tree.close()

        if len(newick_tree) > 1:
            raise AssertionError("ERROR: " + tree + " should only contain one line of text to be a NEWICK tree!")
        else:
            newick_tree = str(newick_tree[0])

        new_tree.write(re.sub('_' + self.COG, '', newick_tree) + "\n")
        new_tree.close()

        return

    def align_sequences(self, alignment_mode, ref_align, unaligned_ref_seqs, args):
        """
        Call MUSCLE to perform a multiple sequence alignment of the reference sequences and the
        gene sequences identified by TreeSAPP
        :param args: Command-line argument object from get_options and check_parser_arguments
        :param unaligned_ref_seqs:
        :param alignment_mode: d (default; re-do the MSA) or p (profile; use the reference MSA)
        :param ref_align: FASTA file containing
        :return: Name of the FASTA file containing the MSA
        """
        if len(self.ContigDict) <= 1 and args.alignment_mode == "p":
            alignment_mode = "d"
            sys.stderr.write("WARNING: Default multiple alignment since only a single new reference sequence.\n")

        if args.verbose:
            sys.stdout.write("Aligning the reference and identified " + self.COG + " sequences using MUSCLE... ")
            sys.stdout.flush()

        # Default alignment #
        if alignment_mode == "d":
            # Combine the reference and candidate sequence dictionaries
            unaligned_ref_seqs.update(self.ContigDict)
            ref_unaligned = self.Output + self.COG + "_gap_removed.fa"
            write_new_fasta(unaligned_ref_seqs, ref_unaligned)

            aligned_fasta = self.Output + self.COG + "_d_aligned.fasta"
            muscle_align_command = "muscle -in %s -out %s 1>/dev/null 2>/dev/null" % (ref_unaligned, aligned_fasta)

        # Profile-Profile alignment #
        elif alignment_mode == "p":
            query_fasta = self.Output + self.COG + "_query_unaligned.fasta"
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

    def execute_raxml(self, phylip_file, raxml_destination_folder, args):
        os.makedirs(raxml_destination_folder)
        # No difference between this command and that in create_treesapp_ref_data
        raxml_command = [args.executables["raxmlHPC"], '-m', self.raxml_model]
        # Run RAxML using multiple threads, if CPUs available
        raxml_command += ['-T', str(int(args.num_threads))]
        if args.bootstraps == 0:
            nboot = "autoMR"
        else:
            nboot = str(args.bootstraps)
        raxml_command += ['-s', phylip_file,
                          '-f', 'a',
                          '-x', str(12345),
                          '-#', nboot,
                          '-n', self.COG,
                          '-w', raxml_destination_folder,
                          '-p', str(12345)] #,
                          # '>', raxml_destination_folder + os.sep + 'RAxML_log.txt']

        if args.verbose:
            sys.stdout.write("RAxML command:\n\t" + ' '.join(raxml_command) + "\n")
            sys.stdout.write("Inferring Maximum-Likelihood tree with RAxML... ")
            sys.stdout.flush()

        raxml_pro = subprocess.Popen(' '.join(raxml_command), shell=True, preexec_fn=os.setsid)
        raxml_pro.wait()

        if args.verbose:
            sys.stdout.write("done.\n")
            sys.stdout.flush()

        return


class ItolJplace:
    """
    A class to hold all data relevant to a jplace file to be viewed in iTOL
    """
    fields = list()

    def __init__(self):
        # Sequence name (from FASTA header)
        self.contig_name = ""
        # Code name of the tree it mapped to (e.g. mcrA)
        self.name = ""
        # NEWICK tree
        self.tree = ""
        self.metadata = ""
        self.version = ""
        # List of lineages for each child identified by RAxML.
        self.lineage_list = list()
        self.node_map = dict()
        self.placements = list()
        # The LCA taxonomy derived from lineage_list
        self.lct = ""
        self.lwr = 0
        self.wtd = 0
        # Either the number of times that sequence was observed, or the FPKM of that sequence
        self.abundance = None

    def summarize(self):
        """
        Prints a summary of the ItolJplace object (equivalent to a single marker) to stderr
        Summary include the number of marks found, the tree used, and the tree-placement of each sequence identified
        Written solely for testing purposes
        :return:
        """
        sys.stderr.write("\nInformation for query sequence '" + self.contig_name + "'\n")
        sys.stderr.write(str(len(self.placements)) + " sequence(s) grafted onto the " + self.name + " tree.\n")
        # sys.stderr.write("Reference tree:\n")
        # sys.stderr.write(self.tree + "\n")
        sys.stderr.write("JPlace fields:\n\t" + str(self.fields) + "\n")
        sys.stderr.write("Placement information:\n")
        if not self.placements:
            sys.stderr.write("\tNone.\n")
        elif self.placements[0] == '{}':
            sys.stderr.write("\tNone.\n")
        else:
            for pquery in self.placements:
                placement = loads(pquery, encoding="utf-8")
                for k, v in placement.items():
                    if k == 'p':
                        sys.stderr.write('\t' + str(v) + "\n")
        sys.stderr.write("Non-redundant lineages of child nodes:\n")
        if len(self.lineage_list) > 0:
            for lineage in sorted(set(self.lineage_list)):
                sys.stderr.write('\t' + str(lineage) + "\n")
        else:
            sys.stderr.write("\tNone.\n")
        sys.stderr.write("Lowest common taxonomy:\n")
        if self.lct:
            sys.stderr.write("\t" + str(self.lct) + "\n")
        else:
            sys.stderr.write("\tNone.\n")
        if self.abundance:
            sys.stderr.write("Abundance:\n\t" + str(self.abundance) + "\n")
        sys.stderr.write("\n")
        sys.stderr.flush()
        return

    def list_placements(self):
        """
        Returns a list of all the nodes contained in placements
        :return:
        """
        nodes = list()
        for d_place in self.placements:
            if type(d_place) == str:
                for k, v in loads(d_place).items():
                    if k == 'p':
                        for pquery in v:
                            nodes.append(str(pquery[0]))
            else:
                raise AssertionError("Unable to handle type " + type(d_place) + "\n")
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
            if type(d_place) != str:
                dict_strings = list()  # e.g. "n":["query"]
                for k, v in d_place.items():
                    dict_strings.append(dumps(k) + ':' + dumps(v))
                    placement_string = ', '.join(dict_strings)
                new_placement_collection.append('{' + placement_string + '}')
            else:
                new_placement_collection.append(d_place)
        self.placements = new_placement_collection

        self.fields = [dumps(x) for x in self.fields]
        return

    def get_lwr_position_from_jplace_fields(self):
        """
        Find the position in self.fields of 'like_weight_ratio'
        :return: position in self.fields of 'like_weight_ratio'
        """
        x = 0
        # Find the position of like_weight_ratio in the placements from fields descriptor
        for field in self.fields:
            if field == '"like_weight_ratio"':
                break
            else:
                x += 1
        if x == len(self.fields):
            sys.stderr.write("Unable to find \"like_weight_ratio\" in the jplace string!\n")
            sys.stderr.write("WARNING: Skipping filtering with `filter_min_weight_threshold`\n")
            return None
        return x

    def filter_min_weight_threshold(self, threshold=0.1):
        """
        Remove all placements with likelihood weight ratios less than threshold
        :param threshold: The threshold which all placements with LWRs less than this are removed
        :return:
        """
        unclassified = 0
        # Find the position of like_weight_ratio in the placements from fields descriptor
        x = self.get_lwr_position_from_jplace_fields()
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
                                # sys.stderr.write("\t".join([self.name, str(removed[0]), str(float(removed[x]))]) + "\n")
                            else:
                                acc += 1
                            sys.stderr.flush()
                        # If no placements met the likelihood filter then the sequence cannot be classified
                        # Alternatively: first two will be returned and used for LCA - can test...
                        if len(tmp_placements) > 0:
                            v = tmp_placements
                            dict_strings.append(dumps(k) + ':' + dumps(v))
                            placement_string = ', '.join(dict_strings)
                        else:
                            unclassified += 1
                # Add the filtered placements back to the object.placements
                new_placement_collection.append('{' + placement_string + '}')
            else:
                new_placement_collection.append(pquery)
        self.placements = new_placement_collection
        return unclassified

    def sum_rpkms_per_node(self, leaf_rpkm_sums):
        for pquery in self.placements:
            placement = loads(pquery, encoding="utf-8")
            for k, v in placement.items():
                if k == 'p':
                    for locus in v:
                        jplace_node = locus[0]
                        tree_leaves = self.node_map[jplace_node]
                        for tree_leaf in tree_leaves:
                            if tree_leaf not in leaf_rpkm_sums.keys():
                                leaf_rpkm_sums[tree_leaf] = 0.0
                            leaf_rpkm_sums[tree_leaf] += self.abundance
        return leaf_rpkm_sums

    def filter_max_weight_placement(self):
        """
        Removes all secondary placements of each pquery,
        leaving only the placement with maximum likelihood_weight_ratio
        :return:
        """
        # Find the position of like_weight_ratio in the placements from fields descriptor
        x = self.get_lwr_position_from_jplace_fields()
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
        lwr_pos = self.get_lwr_position_from_jplace_fields()
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
        placements = list()
        fields = list()
        node_map = dict()
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
    While it mostly contains RAxML outputs,
    several functions are used to make 'biological' sense out of these outputs.
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
        # cellular organisms; Kingdom; Phylum; Class; Order; Family; Genus; Species

        return


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
        # sys.stderr.write("Tree:\n\t" + str(self.tree) + "\n")
        sys.stderr.write("Leaf ID:\n\t" + str(self.number) + "\n")
        sys.stderr.write("Description:\n\t" + str(self.description) + "\n")
        if self.complete:
            sys.stderr.write("Lineage:\n\t" + str(self.lineage) + "\n")

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

    def get_info(self):
        info_string = ""
        info_string += "accession = " + self.accession + ", " + "mltree_id = " + self.short_id + "\n"
        info_string += "description = " + self.description + ", " + "locus = " + self.locus + "\n"
        info_string += "organism = " + self.organism + "\n"
        info_string += "lineage = " + self.lineage + "\n"
        sys.stdout.write(info_string)
        sys.stdout.flush()


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
            p_genewise = subprocess.Popen(' '.join(next_task), shell=True, preexec_fn=os.setsid)
            p_genewise.wait()
            if p_genewise.returncode != 0:
                sys.stderr.write("ERROR: " + self.master + " did not complete successfully for:\n")
                sys.stderr.write(str(' '.join(next_task)) + "\n")
                sys.stderr.flush()
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


class MarkerBuild:
    def __init__(self, build_param_line):
        build_param_fields = build_param_line.split('\t')
        if len(build_param_fields) != 6:
            sys.stderr.write("ERROR: Incorrect number of values in ref_build_parameters.tsv line:\n" + build_param_line)
            raise ValueError
        self.cog = build_param_fields[0]
        self.denominator = build_param_fields[1]
        self.model = build_param_fields[2]
        self.pid = build_param_fields[3]
        self.lowest_confident_rank = build_param_fields[4]
        self.update = build_param_fields[5]

    def check_rank(self):
        taxonomies = ["NA", "Kingdoms", "Phyla", "Classes", "Orders", "Families", "Genera", "Species"]

        if self.lowest_confident_rank not in list(taxonomies):
            raise AssertionError("Unable to find " + self.lowest_confident_rank + " in taxonomic map!")

        return


class Header:
    def __init__(self, header):
        self.original = header
        self.formatted = ""
        self.treesapp_name = ""
        self.post_align = ""
        self.first_split = ""
