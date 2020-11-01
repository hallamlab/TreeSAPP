import logging
import sys

from ete3 import Tree

from treesapp.phylo_dist import parent_to_tip_distances


class PhyloPlace:
    n_key = 'n'
    p_key = 'p'

    def __init__(self, placement=None, field_positions=None):
        self.name = ""
        self.edge_num = -1
        self.like_weight_ratio = 0.0
        self.likelihood = 0.0
        self.distal_length = 0.0
        self.pendant_length = 0.0
        self.mean_tip_length = 0.0
        # Load the placed sequence name
        if placement:
            try:
                # Ensure the name value type is always a string
                if isinstance(placement[self.n_key], list):
                    if len(placement[self.n_key]) == 1:
                        name_val = placement[self.n_key][0]
                    else:
                        logging.error(
                            "{} sequence names found in PQuery name value when one was expected :\n{}"
                            "".format(len(placement[self.n_key]),
                                      "\t{}\n".join(placement['n'])))
                        sys.exit(11)
                else:
                    name_val = placement[self.n_key]
                self.name = name_val
            except KeyError:
                logging.error("Unable to find name key ('{}') in placement: {}\n".format(self.n_key, placement))
                sys.exit(13)

            # Load the placement's PQuery information
            try:
                fields = placement[self.p_key][0]
            except KeyError:
                logging.error("Unable to find placement key ('{}') in placement: {}\n".format(self.p_key, placement))
                sys.exit(15)

            # Load the placement information fields
            if not field_positions:
                field_positions = ['edge_num', 'likelihood', 'like_weight_ratio', 'distal_length', 'pendant_length']
            x = 0
            while x < len(field_positions):
                try:
                    self.__dict__[field_positions[x]] = fields[x]
                except KeyError:
                    logging.error("Field '{}' not found in PhyloPlace class attributes.\n".format(field_positions[x]))
                    sys.exit(17)
                x += 1

        return

    def calc_mean_tip_length(self, internal_leaf_node_map: dict, ref_tree: str) -> None:
        if isinstance(ref_tree, str):
            ref_tree = Tree(ref_tree)
        leaf_children = internal_leaf_node_map[self.edge_num]
        if len(leaf_children) > 1:
            # Reference tree with clade excluded
            parent = ref_tree.get_common_ancestor(leaf_children)
            tip_distances = parent_to_tip_distances(parent, leaf_children)
            self.mean_tip_length = round(float(sum(tip_distances) / len(tip_distances)), 6)
        return

    def total_distance(self) -> float:
        return round(sum([self.pendant_length, self.mean_tip_length, self.distal_length]), 5)

    def summary(self) -> str:
        summary_string = "Summary for placement '{}' on edge {}\n".format(self.name, self.edge_num)
        if self.likelihood and self.like_weight_ratio and self.edge_num:
            summary_string += "\tLikelihood\t" + str(self.likelihood) + "\n"
            summary_string += "\tL.W.R.\t\t" + str(self.like_weight_ratio) + "\n"
            summary_string += "Distal, pendant and mean-tip lengths are {}, {} and {}\n" \
                              "".format(self.distal_length, self.pendant_length, self.mean_tip_length)
        return summary_string

    @staticmethod
    def format_pplace_to_jplace(pplaces: list, field_positions=None) -> dict:
        jplace_dict = {}
        placements_list = []
        name = set()

        if not field_positions:
            field_positions = ['edge_num', 'likelihood', 'like_weight_ratio', 'distal_length', 'pendant_length']

        for pplace in pplaces:  # type: PhyloPlace
            x = 0
            placement = []
            name.add(pplace.name)
            while x < len(field_positions):
                placement.append(pplace.__dict__[field_positions[x]])
                x += 1
            placements_list.append(placement)

        if len(name) > 1:
            logging.error("More than one query name provided when rebuilding JPlace file:\n{}\n".format(','.join(name)))
            sys.exit(13)
        else:
            jplace_dict['n'] = [name.pop()]
        jplace_dict['p'] = placements_list

        return jplace_dict


class PQuery:
    """
    A class for sequences that were properly mapped to its gene tree.
    While it mostly contains EPA outputs, functions are used to make 'biological' sense out of these outputs.
    """
    def __init__(self, lineage_str="", rank_str=""):
        self.seq_name = ""  # Full sequence name (from FASTA header)
        self.place_name = ""  # A unique name for the query sequence placed (in case multiple subsequences are placed)
        self.ref_name = ""  # Code name of the tree it mapped to (e.g. McrA)
        self.abundance = None  # Either the number of occurrences, or the FPKM of that sequence
        self.node_map = dict()  # A dictionary mapping internal nodes (Jplace) to all leaf nodes
        ##
        # Taxonomic information:
        ##
        self.lineage_list = list()  # List containing each child's lineage
        self.wtd = 0
        self.lct = ""  # The LCA taxonomy derived from lineage_list
        self.recommended_lineage = ""
        ##
        # Information derived from JPlace 'PQuery's:
        ##
        self.placements = list()
        self.classified = True
        # Sourced from phylogenetic placement (JPlace file)
        self.consensus_placement = None
        self.parent_node = ""
        self.avg_evo_dist = 0.0
        self.distances = ""

        # Known from outer scope
        self.lineage = lineage_str
        self.rank = rank_str
        self.feature_vec = None

        # Features from homology search
        self.seq = ""
        self.evalue = 0.0
        self.start = 0
        self.end = 0
        self.seq_len = 0

    def clear(self):
        self.node_map.clear()
        self.placements.clear()
        self.lineage_list.clear()
        return

    # def transfer_metadata(self, jplace_inst):
    #     self.tree = jplace_inst.tree
    #     self.fields = jplace_inst.fields
    #     self.version = jplace_inst.version
    #     self.metadata = jplace_inst.metadata

    def summarize(self) -> str:
        summary_str = "Summarizing placement of query '{}' onto the {} reference tree:\n" \
                      "".format(self.seq_name, self.ref_name)
        summary_str += "\tNumber of placement locations = {}\n".format(len(self.placements))
        if self.consensus_placement:
            summary_str += "\tConsensus placement identified = True\n"
            summary_str += "\tConsensus placement information:\n{}\n".format(self.consensus_placement.summary())
        else:
            summary_str += "\tConsensus placement identified = False\n"
        summary_str += "\tQuery sequence abundance approx. = {}\n".format(self.abundance)
        summary_str += "\tSequence length = {}\n"
        summary_str += "\thmmsearch E-value = {}\n"
        summary_str += "\tLCA taxonomic lineage is {}\n"
        summary_str += "\tTaxonomic rank resolved to is {}\n"
        return summary_str

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
        i = 0
        while i < len(self.placements):
            pplace = self.placements[i]  # type: PhyloPlace
            if pplace.like_weight_ratio < threshold:
                self.placements.pop(i)
            else:
                i += 1
        if self.consensus_placement.like_weight_ratio < threshold:
            self.classified = False
        return

    def sum_rpkms_per_node(self, leaf_rpkm_sums: dict) -> dict:
        """
        Function that adds the RPKM value of a contig to the node it was placed.
        For contigs mapping to internal nodes: the proportional RPKM assigned is summed for all children.

        :param leaf_rpkm_sums: A dictionary mapping tree leaf numbers to abundances (RPKM sums)
        :return: A dictionary mapping numerical leaf node identifiers to normalized abundance values
        """
        jplace_node = self.consensus_placement.edge_num
        tree_leaves = self.node_map[jplace_node]
        try:
            normalized_abundance = float(self.abundance/len(tree_leaves))
        except TypeError:
            logging.warning("Unable to find abundance for " + self.place_name + "... setting to 0.\n")
            normalized_abundance = 0.0

        for tree_leaf in tree_leaves:  # type: str
            if tree_leaf not in leaf_rpkm_sums.keys():
                leaf_rpkm_sums[tree_leaf] = 0.0
            leaf_rpkm_sums[tree_leaf] += normalized_abundance
        return leaf_rpkm_sums

    def check_jplace_edge_lengths(self, tree_index: dict) -> None:
        """
        Currently validates a pquery's JPlace distal length, ensuring it is less than or equal to the edge length
        This is necessary to handle a case found in RAxML v8.2.12 (and possibly older versions) where the distal length
        of a placement is greater than the corresponding branch length in some rare cases.

        :param tree_index: A dictionary mapping a reference package's tree edge numbers to their respective lengths
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

    def children_lineage(self, leaves_taxa_map: dict) -> list:
        """
        From the jplace placement field ()

        :param leaves_taxa_map: Dictionary mapping tree leaf nodes to taxonomic lineages
        :return:
        """
        children = list()
        pplace = self.consensus_placement  # type: PhyloPlace
        try:
            tree_leaves = self.node_map[pplace.edge_num]
        except KeyError:
            logging.error("Unable to find placement edge '{}' in the PhyloPlace's node_map.\n".format(pplace.edge_num))
            sys.exit(13)

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

    def lowest_confident_taxonomy(self, depth: int) -> str:
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

    def calculate_consensus_placement(self) -> None:
        """
        Often times, a single PQuery (query sequence mapped onto a phylogeny) may be inserted into the phylogeny
        at multiple edges with similar likelihood.

        The PQuery.consensus_placement attribute is either an existing placement (PhyloPlace instance) or the...

        :return: None
        """
        # TODO: finish implementing this.
        for pplace in self.placements:  # type: PhyloPlace
            self.consensus_placement = pplace

        return

    def filter_max_weight_placement(self) -> PhyloPlace:
        """
        Often times, a single PQuery (query sequence mapped onto a phylogeny) may be inserted into the phylogeny
        at multiple edges with similar likelihood. This function aims to select the single best placement based on
        it's respective Likelihood Weight Ratio (LWR) or PQuery like_weight_ratio attribute.

        The PQuery.consensus_placement attribute is set to this placement with maximum LWR. The PQuery.placements
        attribute is unmodified.

        :return: None
        """
        # Filter the placements
        max_lwr = 0
        for pplace in self.placements:  # type: PhyloPlace
            if pplace.name:
                if pplace.like_weight_ratio > max_lwr:
                    max_lwr = pplace.like_weight_ratio
                    self.consensus_placement = pplace
            else:
                logging.error("Unexpected state of PhyloPlace instance!\n{}\n".format(pplace.summary()))
                sys.exit(3)
        return self.consensus_placement


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
        con_place = PhyloPlace()
        try:
            _, pquery.place_name, pquery.ref_name, pquery.start, pquery.end, pquery.recommended_lineage, \
            _, con_place.edge_num, pquery.evalue, con_place.like_weight_ratio, pquery.avg_evo_dist, pquery.distances = fields
        except ValueError:
            logging.error("Bad line in classification table:\n" +
                          '\t'.join(fields) + "\n")
            sys.exit(21)
        pquery.end = int(pquery.end)
        pquery.start = int(pquery.start)
        pquery.seq_len = pquery.end - pquery.start
        pquery.lct = pquery.recommended_lineage
        con_place.distal_length, con_place.pendant_length, con_place.mean_tip_length = [float(d) for d in
                                                                                        pquery.distances.split(',')]
        pquery.consensus_placement = con_place
        refpkg = refpkg_dict[pquery.ref_name]  # type: refpkg.ReferencePackage
        try:
            pqueries[refpkg.prefix].append(pquery)
        except KeyError:
            pqueries[refpkg.prefix] = [pquery]
    return pqueries


def split_placements(placements: dict) -> list:
    """
    When the JPlace data is loaded by json.load(), it is returned as a list.
    This accepts the dictionary from that list collection, for which there is only a single element, and generates a
    PhyloPlace instance for each placement in the 'p' value.

    :param placements: A dictionary with two keys: 'n' (corresponding to the query sequence name) and
     'p' (corresponding to the placement data e.g. likelihood weight ratio, distal length).
    :return: A list of PhyloPlace instances, where the attributes have been populated from each of the placement
     locations on the phylogeny.
    """
    # Reformat the dictionary so each one is a unique placement
    phylo_places = []
    for pplace in placements['p']:  # type: list
        phylo_places.append(PhyloPlace({'p': [pplace], 'n': placements['n']}))
    return phylo_places


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
