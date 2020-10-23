import logging
import sys
from copy import deepcopy
from json import loads, dumps


class PhyloPlace:
    def __init__(self, placement: dict, field_positions=None):
        self.name = ""
        self.edge_num = 0
        self.like_weight_ratio = 0
        self.likelihood = 0
        self.distal_length = 0
        self.pendant_length = 0
        self.mean_tip_length = 0
        n_key = 'n'
        p_key = 'p'
        # Load the placed sequence name
        try:
            if len(placement[n_key]) != 1:
                logging.error("Multiple sequence names found in PQuery name:\n{}".format("\t{}\n".join(placement['n'])))
                sys.exit(11)
            self.name = placement[n_key][0]
        except KeyError:
            logging.error("Unable to find name key ('{}') in placement: {}\n".format(n_key, placement))
            sys.exit(13)

        # Load the placement's PQuery information
        try:
            fields = placement[p_key][0]
        except KeyError:
            logging.error("Unable to find placement info key ('{}') in placement: {}\n".format(p_key, placement))
            sys.exit(15)

        # Load the placement information fields
        if not field_positions:
            field_positions = ['edge_num', 'likelihood', 'like_weight_ratio', 'distal_length', 'pendant_length']
        x = 0
        while x < len(field_positions):
            try:
                self.__dict__[field_positions[x]] = fields[x]
            except KeyError:
                logging.error("Field '{}' doesn't exist in PhyloPlace class attributes.\n".format(field_positions[x]))
                sys.exit(17)
            x += 1

        return

    def total_distance(self) -> int:
        return round(sum([self.pendant_length, self.mean_tip_length, self.distal_length]), 5)

    def summary(self) -> str:
        # TODO: Finish implementing this
        summary_string = "Summary for placement '{}' on edge {}\n".format(self.name, self.edge_num)
        if self.likelihood and self.like_weight_ratio and self.edge_num:
            summary_string += "\tLikelihood\t" + str(self.likelihood) + "\n"
            summary_string += "\tL.W.R.\t\t" + str(self.like_weight_ratio) + "\n"
            summary_string += ""
        return summary_string


class JPlace:
    """
    A class to hold all data relevant to a jplace file to be viewed in iTOL
    """
    fields = list()

    def __init__(self):
        self.place_name = ""  # A unique name for the query sequence placed (in case multiple subsequences are placed)
        self.ref_name = ""  # Code name of the tree it mapped to (e.g. mcrA)
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
        self.distances = ""
        self.classified = True
        self.inode = ""
        self.tree = ""  # NEWICK tree
        self.metadata = ""
        self.version = ""  # Jplace version

    def summarize(self):
        """
        Prints a summary of the JPlace object (equivalent to a single marker) to stderr
        Summary include the number of marks found, the tree used, and the tree-placement of each sequence identified
        Written solely for testing purposes

        :return:
        """
        summary_string = ""
        summary_string += "\nInformation for query sequence '" + str(self.place_name) + "'\n"
        summary_string += str(len(self.placements)) + " sequence(s) grafted onto the " + self.ref_name + " tree.\n"
        # summary_string += "Reference tree:\n")
        # summary_string += self.tree + "\n")
        summary_string += "JPlace fields:\n\t" + str(self.fields) + "\n"
        summary_string += "Placement information:\n"
        if not self.placements:
            summary_string += "\tNone.\n"
        elif self.placements[0] == '{}':
            summary_string += "\tNone.\n"
        else:
            for placement in self.placements:  # type: PhyloPlace
                summary_string += placement.summary()
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

    # def correct_decoding(self) -> None:
    #     """
    #     Since the JSON decoding is unable to decode recursively, this needs to be fixed for each placement
    #     Formatting and string conversion are also performed here
    #
    #     :return: None
    #     """
    #     new_placement_collection = []  # a list of dictionary-like strings
    #     placement_string = ""  # e.g. {"p":[[226, -31067.028237, 0.999987, 0.012003, 2e-06]], "n":["query"]}
    #     for d_place in self.placements:
    #         if not isinstance(d_place, str):
    #             dict_strings = list()  # e.g. "n":["query"]
    #             for k, v in d_place.items():
    #                 dict_strings.append(dumps(k) + ':' + dumps(v))
    #                 placement_string = ', '.join(dict_strings)
    #             new_placement_collection.append('{' + placement_string + '}')
    #         else:
    #             new_placement_collection.append(d_place)
    #     self.placements = new_placement_collection
    #
    #     decoded_fields = list()
    #     for field in self.fields:
    #         if not re.match('".*"', field):
    #             decoded_fields.append(dumps(field))
    #         else:
    #             decoded_fields.append(field)
    #     self.fields = decoded_fields
    #     return

    def name_placed_sequence(self) -> None:
        names = set()
        for pplace in self.placements:  # type: PhyloPlace
            names.add(pplace.name)

        if len(names) > 1:
            logging.error("Multiple names encountered for a single PQuery:\n{}\n".format(','.join(names)))
            raise AssertionError
        self.place_name = names.pop()
        return

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
                            logging.warning("Unable to find abundance for " + self.place_name + "... setting to 0.\n")
                            normalized_abundance = 0.0
                        for tree_leaf in tree_leaves:
                            if tree_leaf not in leaf_rpkm_sums.keys():
                                leaf_rpkm_sums[tree_leaf] = 0.0
                            leaf_rpkm_sums[tree_leaf] += normalized_abundance
        return leaf_rpkm_sums

    def check_jplace_edge_lengths(self, tree_index) -> None:
        """
        Currently validates a pquery's JPlace distal length, ensuring it is less than or equal to the edge length
        This is necessary to handle a case found in RAxML v8.2.12 (and possibly older versions) where the distal length
        of a placement is greater than the corresponding branch length in some rare cases.

        :return: None
        """
        for pquery in self.placements:  # type: PhyloPlace
            place_len = float(pquery.distal_length)
            tree_len = tree_index[str(pquery.edge_num)]
            if place_len > tree_len:
                logging.debug("Distal length adjusted to fit JPlace {} tree for {}.\n"
                              "".format(self.ref_name, self.place_name))
                pquery.distal_length = tree_len

        return

    def clear_object(self):
        self.placements.clear()
        self.fields.clear()
        self.node_map.clear()
        self.place_name = ""
        self.ref_name = ""
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


class PQuery(JPlace):
    """
    A class for sequences that were properly mapped to its gene tree.
    While it mostly contains EPA outputs, functions are used to make 'biological' sense out of these outputs.
    """
    def __init__(self, lineage_str="", rank_str=""):
        super(PQuery, self).__init__()
        self.seq_name = ""  # Full sequence name (from FASTA header)
        # Sourced from phylogenetic placement (JPlace file)
        self.consensus_placement = None
        self.parent_node = ""

        # Known from outer scope
        self.lineage = lineage_str
        self.rank = rank_str
        self.feature_vec = None

        # Features from homology search
        self.seq = ""
        self.evalue = 0.0
        self.start = 0
        self.end = 0

    def transfer_metadata(self, itol_jplace_object):
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

    def calculate_consensus_placement(self) -> None:
        for pplace in self.placements:  # type: PhyloPlace
            pass
        self.
        return

    def filter_max_weight_placement(self) -> None:
        """
        Removes all secondary placements of each pquery, leaving only the placement with the maximum like_weight_ratio

        :return: None
        """
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
        self.consensus_placement = self.placements[0]
        return


def assignments_to_treesaps(classified_lines: list, refpkg_dict: dict) -> dict:
    """
    Used for converting the TreeSAPP-assignment information of classified sequences (found in self.classifications)
    into JPlace instances such that these can be reproducibly modified and written again, if needed.

    :param classified_lines: A list of lists. Each sub-list represents a line from self.classifications
    :param refpkg_dict: A dictionary of MarkerBuild instances indexed by their RefPkg codes
    :return: A dictionary of JPlace instances, indexed by their respective names (ReferencePackage.prefix)
    """
    pqueries = dict()
    # "Sample\tQuery\tMarker\tStart_pos\tEnd_pos\tTaxonomy\tAbundance\tiNode\tE-value\tLWR\tEvoDist\tDistances\n"
    for fields in classified_lines:
        pquery = PQuery()
        try:
            _, pquery.place_name, pquery.ref_name, pquery.start, pquery.end, pquery.recommended_lineage, \
            _, pquery.inode, pquery.evalue, pquery.lwr, pquery.avg_evo_dist, pquery.distances = fields
        except ValueError:
            logging.error("Bad line in classification table:\n" +
                          '\t'.join(fields) + "\n")
            sys.exit(21)
        pquery.end = int(pquery.end)
        pquery.start = int(pquery.start)
        pquery.seq_len = pquery.end - pquery.start
        pquery.lct = pquery.recommended_lineage
        pquery.distal, pquery.pendant, pquery.mean_tip = [float(d) for d in pquery.distances.split(',')]
        refpkg = refpkg_dict[pquery.ref_name]  # type: refpkg.ReferencePackage
        try:
            pqueries[refpkg.prefix].append(pquery)
        except KeyError:
            pqueries[refpkg.prefix] = [pquery]
    return pqueries


class TreeLeafReference:
    """
    Objects for each leaf in a tree
    """
    def __init__(self, number: str, description: str):
        self.number = number
        self.description = description
        self.lineage = ""
        self.accession = ""
        self.complete = False

    def summarize_tree_leaf(self):
        summary_string = "Leaf ID:\n\t{}\n".format(str(self.number)) +\
                         "Description:\n\t'{}'\n".format(str(self.description))
        summary_string += "Accession:\n\t'{}'\n".format(self.accession)
        summary_string += "Lineage:\n\t'{}'\n".format(str(self.lineage))
        return summary_string


def convert_entrez_to_tree_leaf_references(entrez_records: dict) -> list:
    """
    From the dictionary containing lineage and organism information of reference sequences (self.lineage_ids)
    this function creates a list of TreeLeafReference instances for every reference sequence

    :return: List of TreeLeafReference instances
    """
    ref_leaf_nodes = []
    for treesapp_id in sorted(entrez_records, key=int):
        record = entrez_records[treesapp_id]
        ref = TreeLeafReference(treesapp_id, record.description)
        ref.accession = record.accession
        ref.lineage = record.lineage
        ref_leaf_nodes.append(ref)
    return ref_leaf_nodes
