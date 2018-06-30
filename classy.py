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
from utilities import reformat_string, return_sequence_info_groups
from entish import get_node, create_tree_info_hash, subtrees_to_dictionary
from external_command_interface import launch_write_command
from entrez_utils import get_lineage

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
        self.master_reference_index = dict()  # Used for storing ReferenceSequence objects, indexed by numeric IDs
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
            ref_seq = ReferenceSequence()
            ref_seq.lineage = leaf.lineage
            ref_seq.organism = leaf.description.split(' | ')[0]
            ref_seq.accession = leaf.description.split(' | ')[1]
            ref_seq.short_id = '>' + leaf.number + '_' + self.COG
            self.master_reference_index[leaf.number] = ref_seq
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
            ref_seq = ReferenceSequence()
            ref_seq.short_id = internal_id
            # Map the new reference headers to their numerical IDs
            self.header_id_map[internal_id] = new_ref_seq
            self.ContigDict[internal_id] = centroids_fasta_dict[new_ref_seq]
            ref_seq.sequence = centroids_fasta_dict[new_ref_seq]
            self.master_reference_index[str(additional)] = ref_seq

        if additional == acc:
            # The reference sequence identifiers are random
            sys.stderr.write("WARNING: numerical TreeSAPP identifiers for " + self.Denominator +
                             "are not in the format of a sequential series!\n")
            sys.stderr.write("Generating random numerical unique identifiers for the new sequence(s).\n")
            sys.stderr.flush()

        if args.verbose:
            sys.stdout.write("\t" + str(len(self.header_id_map)) + " new " + self.COG + " reference sequences.\n")

        return

    def update_tax_ids(self, args, ref_organism_lineage_info, assignments):
        """
        Write the number, organism and accession ID, if possible
        :param args:
        :param ref_organism_lineage_info: A dictionary mapping the
        :param assignments: A dictionary containing marker name as keys, and header: [lineages] mappings as values
        :return:
        """
        sys.stdout.write("Writing updated tax_ids file... ")
        sys.stdout.flush()

        tree_taxa_string = ""
        original_to_formatted_header_map = dict()
        unclassified_seqs = list()

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
            header_format_re, header_db, header_molecule = '', '', ''
            description, lineage, accession, organism = '', '', '', ''
            # Check to see if this header contains an accession value or if it is a contig name with no useful info
            try:
                header_format_re, header_db, header_molecule = get_header_format(header, self.Denominator)
            except AssertionError:
                # Just a contig name, so it will need to be formatted with explicit lineage
                for candidate_header_name in assignments[self.COG]:
                    if header == candidate_header_name:
                        assigned_lineage = assignments[self.COG][candidate_header_name]
                        # Determine if this is a cellular organism or not
                        if assigned_lineage.split(';')[0] in ["Bacteria", "Archaea", "Eukaryota"]:
                            assigned_lineage = "cellular organisms; " + assigned_lineage
                        annotated_header = ' '.join([header,
                                                     "lineage=" + assigned_lineage,
                                                     "[" + re.sub('>', '', header),
                                                     self.InputData.split(os.sep)[-1],
                                                     assigned_lineage.split(';')[-1] + "]"])
                        header_format_re, header_db, header_molecule = get_header_format(annotated_header, self.Denominator)
                        header = annotated_header
                        break
                if not header_format_re:
                    unclassified_seqs.append(header)
            if header_format_re and header_db and header_molecule:
                sequence_info = header_format_re.match(header)
                if sequence_info:
                    accession, organism, locus, description, lineage = return_sequence_info_groups(sequence_info, header_db, header)
                    if lineage:
                        pass
                    elif accession:
                        lineage = get_lineage(accession, header_molecule)
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
                    else:
                        sys.stderr.write("ERROR: TreeSAPP is unsure what to do with header '" + header + "'\n")
                        raise AssertionError()
            else:
                # This sequence is probably an uninformative contig name
                lineage = "unclassified sequences; "
                description = self.InputData.split('/')[-1]

            if not description or not lineage:
                sys.stderr.write("WARNING: description is unavailable for sequence '" + header + "'\n")
            tree_taxa_string += num_id + "\t" + description + " | " + accession + "\t" + lineage + "\n"
            self.master_reference_index[num_id].organism = organism
            self.master_reference_index[num_id].description = description
            self.master_reference_index[num_id].accession = accession
            self.master_reference_index[num_id].lineage = lineage

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

    def align_multiple_sequences(self, unaligned_ref_seqs, args):
        """
        Call MAFFT to perform a multiple sequence alignment of the reference sequences and the
        gene sequences identified by TreeSAPP
        :param args: Command-line argument object from get_options and check_parser_arguments
        :param unaligned_ref_seqs:
        :param ref_align: FASTA file containing
        :return: Name of the FASTA file containing the MSA
        """

        if args.verbose:
            sys.stdout.write("Aligning the reference and identified " + self.COG + " sequences using MAFFT... ")
            sys.stdout.flush()

        # Combine the reference and candidate sequence dictionaries
        unaligned_ref_seqs.update(self.ContigDict)
        ref_unaligned = self.Output + self.COG + "_gap_removed.fa"
        write_new_fasta(unaligned_ref_seqs, ref_unaligned)

        aligned_fasta = self.Output + self.COG + "_d_aligned.fasta"

        mafft_align_command = [args.executables["mafft"]]
        mafft_align_command += ["--maxiterate", str(1000)]
        mafft_align_command += ["--thread", str(args.num_threads)]
        if len(self.ContigDict) > 700:
            mafft_align_command.append("--auto")
        else:
            mafft_align_command.append("--localpair")
        mafft_align_command += [ref_unaligned, '1>' + aligned_fasta]
        mafft_align_command += ["2>", "/dev/null"]

        stdout, mafft_proc_returncode = launch_write_command(mafft_align_command, False)

        if mafft_proc_returncode != 0:
            sys.stderr.write("ERROR: Multiple sequence alignment using " + args.executables["mafft"] +
                             " did not complete successfully! Command used:\n" + ' '.join(
                mafft_align_command) + "\n")
            sys.exit()

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
        summary_string += "\n"
        return summary_string

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
                    self.contig_name = value
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
            sys.stderr.write("Unable to find '" + field_name + "' in the jplace string!\n")
            sys.stderr.write("WARNING: Skipping filtering with `filter_min_weight_threshold`\n")
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

    def get_inode(self):
        self.inode = str(self.get_jplace_element("edge_num"))
        return

    def get_placement_lwr(self):
        self.lwr = float(self.get_jplace_element("like_weight_ratio"))
        return

    def get_placement_likelihood(self):
        self.likelihood = float(self.get_jplace_element("likelihood"))
        return

    def filter_min_weight_threshold(self, threshold=0.1):
        """
        Remove all placements with likelihood weight ratios less than threshold
        :param threshold: The threshold which all placements with LWRs less than this are removed
        :return:
        """
        unclassified = 0
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
                # If there is only one placement, the LWR is 1.0 so no filtering required!
                new_placement_collection.append(pquery)
        if unclassified == 0:
            self.placements = new_placement_collection
        return unclassified

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
        self.cluster_rep = False
        self.cluster_rep_similarity = 0
        self.cluster_lca = None

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
        self.description = ""
        self.kind = ""
        self.molecule = "prot"

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

    def get_info(self):
        info_string = "TreeSAPP ID = '" + self.treesapp_name + "'\tPrefix = '" + self.first_split + "'\n"
        info_string += "Original =  " + self.original + "\nFormatted = " + self.formatted
        return info_string


class Cluster:
    def __init__(self, rep_name):
        self.representative = rep_name
        self.members = list()
        self.lca = ''
