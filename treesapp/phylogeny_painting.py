import sys
import re
import os
import logging

from collections import namedtuple
from matplotlib import colors as mpl_colours
import seaborn

from treesapp.phylo_seq import TreeLeafReference
from treesapp.classy import TreeSAPP
from treesapp.refpkg import ReferencePackage
from treesapp.taxonomic_hierarchy import Taxon
from treesapp import logger

LOGGER = logging.getLogger(logger.logger_name())


class PhyPainter(TreeSAPP):
    def __init__(self):
        super(PhyPainter, self).__init__("colour")
        self.refpkg_dict = {}
        self.refpkg_leaf_nodes_to_colour = {}
        self.taxa_to_colour = set()
        self.rank = ""
        self.palette = ""
        self.unknown_label = None
        self.unknown_col = ""
        self.set_op = ""
        self.feature_name = ""
        self.num_taxa = 0
        self.num_seqs = 0
        self.rank_depth = -1
        self.order_method = "phylo"

    def primer(self, args) -> None:
        rank_map = {}
        self.ref_pkg = None
        for pkl_file in args.pkg_path:
            ref_pkg = ReferencePackage()
            ref_pkg.f__pkl = pkl_file
            ref_pkg.slurp()
            if ref_pkg.validate():
                self.refpkg_dict[ref_pkg.prefix] = ref_pkg
            if not rank_map:
                rank_map = ref_pkg.taxa_trie.accepted_ranks_depths

        self.ts_logger.debug("\tGenerating colour palette for reference package(s): {}.\n"
                      "".format(', '.join(self.refpkg_dict.keys())))

        self.set_op = args.set_op

        self.output_dir = args.output
        if not os.path.isdir(self.output_dir):
            try:
                os.makedirs(self.output_dir)
            except IOError:
                self.ts_logger.error("Unable to create output directory '{}'.\n".format(self.output_dir))
                sys.exit(1)

        self.rank = args.rank
        self.palette = args.palette

        # Determine the prefix for the colour files
        if args.attribute == "taxonomy":
            self.feature_name = self.rank
        else:
            self.feature_name = args.attribute
            self.rank = None
            self.order_method = "alpha"

        if args.un_col:
            self.unknown_label = "Unknown"
            self.set_unknown_colour_from_mpl(args.un_col)

        return

    def set_unknown_colour_from_mpl(self, requested_colour: str):
        try:
            self.unknown_col = mpl_colours.get_named_colors_mapping()[requested_colour]
        except KeyError:
            self.ts_logger.error("Colour '{}' is not available in matplotlib.colors.get_named_colors_mapping().\n"
                                 "Unable to find hexcode for requested unknown colour.\n")
            sys.exit(11)
        return

    def find_rank_depth(self, ref_pkg: ReferencePackage, depth: int) -> None:
        """
        Depending on the lineages parsed by `treesapp create` and users, there may be some discrepancies between a rank
        and its depth (distance from the root) in a TaxonomicHierarhcy (ref_pkg.taxa_trie.accepted_ranks_depths)

        This function is meant to check whether all of the lineages begin at the same rank and therefore all lineages
        can be sliced consistently, e.g. the third position in a lineage will always represent the taxon's class

        :param depth: The depth specified to represent a taxonomic rank
        :param ref_pkg: A ReferencePackage instance
        :return: None
        """
        offsets = set()
        domain_reps = ref_pkg.taxa_trie.rank_representatives("domain", with_prefix=True)
        luca_reps = ref_pkg.taxa_trie.rank_representatives("root", with_prefix=True)

        ref_leaf_nodes = sorted(ref_pkg.generate_tree_leaf_references_from_refpkg(), key=lambda x: x.lineage)
        for leaf in ref_leaf_nodes:  # type: TreeLeafReference
            if not leaf.lineage:
                continue
            lineage_path = leaf.lineage.split(ref_pkg.taxa_trie.lin_sep)
            if len(lineage_path) < depth-1:
                continue
            if lineage_path[0] in domain_reps:
                offsets.add(depth - 1)
            elif lineage_path[0] in luca_reps:
                offsets.add(depth)
            else:
                self.ts_logger.debug("Unable to find root or domain name for '%s'. Skipping.\n" % leaf.lineage)
                continue

        if len(offsets) > 1:
            self.ts_logger.error("Inconsistent rank positions encountered across lineages.\n")
            sys.exit(13)
        self.rank_depth = offsets.pop()

        return

    def remove_taxa_from_colours(self, taxon_leaf_map: dict, unique_taxa: dict, taxa: set) -> None:
        """
        At various stages during the `treesapp colour` workflow, some taxa are removed/filtered from collections
        with names of taxa that should be annotated on the phylogeny. They may be removed for multiple reasons,
        such as too few members or being a polyphyletic group.

        This function ensures that these taxa are removed from _all_ collections that list taxa so they are
        guaranteed to not be included in any annotation (i.e. iTOL-compatible colour) files.

        :param taxa:
        :param unique_taxa:
        :param taxon_leaf_map:
        :return: None
        """
        for taxon in taxa:
            self.num_taxa -= 1
            self.num_seqs -= len(taxon_leaf_map[taxon])
            for n, name in unique_taxa.items():  # type: (int, Taxon)
                if name == taxon:
                    unique_taxa.pop(n)
                    break
            taxon_leaf_map.pop(taxon)

        self.ts_logger.info("Remaining sequences for colouring\t{}\n"
                     "Remaining taxa for colouring\t{}\n"
                     .format(self.num_seqs, self.num_taxa))
        return

    @staticmethod
    def filter_unwanted_taxa(taxon_leaf_map: dict, unique_taxa: dict, taxa_filter: str) -> set:
        """
        A lineage-based string filter for excluding taxa from the output colour palette.
        Even though the palette is for taxa, and not complete lineages,
        we may want to filter out a very large group of organisms, such as Eukaryota.
        This facilitates that operation so not all taxa from that larger group need to be individually named

        :param taxon_leaf_map:
        :param unique_taxa:
        :param taxa_filter: A comma-separated list of taxa to remove from those that should be coloured
        :return: A set containing taxa with rank-prefixes
        """
        filter_str = str(taxa_filter) + " taxa that won't be coloured:\n\t"
        filtered_taxa = set()
        filtered_seqs = 0
        filter_terms = taxa_filter.split(',')
        for leaf_number, taxon in unique_taxa.items():  # type: (int, Taxon)
            if len(filter_terms) > 0:
                for term in filter_terms:
                    if re.search(term, ''.join([t.prefix_taxon() for t in taxon.lineage()])):
                        filtered_taxa.add(taxon.prefix_taxon())

        for taxon in filtered_taxa:
            filtered_seqs += len(taxon_leaf_map[taxon])

        filter_str += ", ".join(filtered_taxa) + "\n"
        filter_str += "Sequences from unwanted taxa removed\t" + str(filtered_seqs) + "\n"
        filter_str += "Unwanted taxa removed\t" + str(len(filtered_taxa)) + "\n"
        LOGGER.info(filter_str)

        return filtered_taxa

    @staticmethod
    def filter_rare_groups(taxon_leaf_map: dict, num_refpkg_seqs: int, min_p: float) -> set:
        filtered_taxa = set()
        filtered_taxa_strings = []
        for taxon in taxon_leaf_map:  # type: str
            leaf_nodes = set(taxon_leaf_map[taxon])
            if len(leaf_nodes)/num_refpkg_seqs < min_p:
                filtered_taxa_strings.append("\t" + taxon + " = " + str(round(len(leaf_nodes)/num_refpkg_seqs, 5)))
                filtered_taxa.add(taxon)
        if not filtered_taxa_strings:
            filtered_taxa_strings.append("\tNone")

        LOGGER.info("Taxa filtered due to proportion < {}:\n"
                     "\t{}\nLow-abundance taxa removed\t{}\n"
                     "".format(min_p, "\n\t".join(filtered_taxa_strings), len(filtered_taxa)))

        return filtered_taxa

    @staticmethod
    def filter_polyphyletic_groups(taxon_leaf_map: dict, internal_node_map: dict) -> set:
        bad_taxa = set()
        for taxon in taxon_leaf_map:  # type: str
            ancestral_node = None
            leaf_nodes = set(taxon_leaf_map[taxon])

            # Sort the internal nodes by the number of children (moving from shallow to deep)
            for i_node in sorted(internal_node_map, key=lambda n: len(internal_node_map[n])):
                if len(internal_node_map[i_node]) < len(leaf_nodes):
                    continue
                elif len(set(internal_node_map[i_node]).intersection(leaf_nodes)) == len(leaf_nodes):
                    ancestral_node = i_node
                    break
                else:
                    pass
            if not ancestral_node:
                LOGGER.warning("Unable to find internal node ancestral to all children of " + taxon + ".\n")
                bad_taxa.add(taxon)
                continue
            if len(internal_node_map[ancestral_node]) > len(leaf_nodes):
                LOGGER.warning(taxon + " clade was found to be polyphyletic.\n")
                bad_taxa.add(taxon)

        return bad_taxa

    def add_unknowns_to_feature_leaf_map(self, feature_leaf_map: dict, ref_pkg: ReferencePackage) -> None:
        # Assign 'Unknown' feature to all leaves missing from the taxa_leaf_nodes_map if set
        if self.unknown_col:
            present = set(sum(feature_leaf_map.values(), []))
            absent = set(ref_pkg.ref_names()).difference(present)
            try:
                feature_leaf_map[self.unknown_label] += list(absent)
            except KeyError:
                feature_leaf_map[self.unknown_label] = list(absent)
        return

    @staticmethod
    def find_mono_clades(taxon_leaf_map: dict, ref_pkg: ReferencePackage) -> dict:
        """

        :param taxon_leaf_map: A dictionary with rank-prefixed taxon names as keys mapping to all leaf nodes
         (e.g. 111_PuhA) whose lineage contains that taxon (i.e. the leaf nodes are sequences from the taxon)
        :param ref_pkg: A ReferencePackage instance
        :return: A dictionary mapping taxon names in taxon_leaf_map to the minimal set of internal nodes in the
         ReferencePackage's tree that cover all leaf nodes in the taxon_leaf_map[taxon].
        """
        taxa_leaf_nodes_map = {}

        for taxon, taxon_leaves in taxon_leaf_map.items():
            taxa_leaf_nodes_map[taxon] = ref_pkg.get_monophyletic_clades(taxon_name=taxon, leaf_nodes=set(taxon_leaves))
        return taxa_leaf_nodes_map

    def harmonize_taxa_colours(self, taxon_leaf_map: dict, set_operation: str) -> None:
        """

        :param taxon_leaf_map:
        :param set_operation:
        :return:
        """
        if not self.taxa_to_colour:
            self.taxa_to_colour.update(set(taxon_leaf_map.keys()))
            return

        if set_operation == 'u':
            self.taxa_to_colour.union(set(taxon_leaf_map.keys()))
        elif set_operation == 'i':
            self.taxa_to_colour.intersection(set(taxon_leaf_map.keys()))
        else:
            self.ts_logger.error("Unsupported set operation specified: '{}'\n".format(set_operation))
            sys.exit(17)

        return

    def get_colours(self) -> list:
        n_clade = len(self.taxa_to_colour)
        if n_clade <= 1 and self.rank:
            self.ts_logger.warning("{} unique clade(s) identified at rank '{}'."
                                   " Consider using a more resolved rank or adjusting other filtering parameters.\n"
                                   "".format(n_clade, self.feature_name))
        else:
            self.ts_logger.debug("Identified {} unique clades for '{}'.\n".format(n_clade, self.feature_name))

        try:
            if self.unknown_label in self.taxa_to_colour:
                colours = seaborn.color_palette(self.palette, n_clade-1).as_hex()
                colours.append(self.unknown_col)
            else:
                colours = seaborn.color_palette(self.palette, n_clade).as_hex()
            self.ts_logger.info("Using palette '{}' for styling.\n".format(self.palette))
        except ValueError as e:
            self.ts_logger.error(str(e) + "\n")
            sys.exit(4)

        if len(colours) != n_clade:
            self.ts_logger.error("Bad colour parsing! len(colours) != number of clades.\n")
            raise ValueError
        return colours

    def map_colours_to_taxa(self, taxa_order: dict):
        """
        Function for mapping

        :param taxa_order: Dictionary indexed by numerical order of the taxa
        :return:
        """

        colours = self.get_colours()
        palette_taxa_map = dict()
        alpha_sort_i = 0

        if len(taxa_order) != len(colours):
            self.ts_logger.error("Number of taxa (%d) and colours (%d) are not equal!\n" % (len(taxa_order), len(colours)))
            sys.exit(7)

        for taxon_num in sorted(taxa_order, key=int):
            taxon = taxa_order[taxon_num]
            palette_taxa_map[taxon] = colours[alpha_sort_i]
            alpha_sort_i += 1

        if self.unknown_label and self.unknown_col:
            palette_taxa_map[self.unknown_label] = self.unknown_col

        return palette_taxa_map


def convert_clades_to_ranges(taxa_clades: dict, leaf_order: list) -> dict:
    taxa_ranges = dict()
    for taxon in taxa_clades:
        taxa_ranges[taxon] = []
        for clade_leaves in taxa_clades[taxon]:
            if len(clade_leaves) > 1:
                positions = clade_leaf_sides(clade_leaves, leaf_order)
                taxa_ranges[taxon].append([positions.left, positions.right])
            else:
                taxa_ranges[taxon].append(clade_leaves)
    return taxa_ranges


def create_write_file(file_name: str, text: str) -> None:
    # Open the output file
    try:
        file_handler = open(file_name, 'w')
    except IOError:
        LOGGER.error("Unable to open file '{}' for writing!\n".format(file_name))
        raise IOError

    # Write the iTOL-formatted header and node, colour and taxon fields
    file_handler.write(text)
    file_handler.close()
    return


def write_colours_styles(taxon_leaf_map: dict, palette_taxa_map: dict, style_output: str) -> None:
    colourless = set()

    colours_style_string = "TREE_COLORS\n" \
                           "SEPARATOR TAB\n" \
                           "DATA\n"
    for taxon, col in palette_taxa_map.items():
        try:
            leaves = taxon_leaf_map[taxon]
        except KeyError:
            continue
        for leaf_num in leaves:
            colours_style_line = [str(leaf_num), "range", col, taxon]
            if len(colours_style_line) > 0:
                colours_style_string += "\t".join(colours_style_line) + "\n"

    LOGGER.info("iTOL colours style input written to " + style_output + ".\n")
    LOGGER.debug("{} lineages were not assigned to a colour:\n\t{}\n".format(len(colourless), "\n\t".join(colourless)))

    create_write_file(style_output, colours_style_string)

    return


def write_colour_strip(taxa_nodes: dict, palette_taxa_map: dict, colour_strip_file: str, data_label: str) -> None:
    colour_strip_text = "DATASET_COLORSTRIP\n" \
                        "SEPARATOR SPACE\n" \
                        "DATASET_LABEL {}\n".format(data_label) +\
                        "STRIP_WIDTH 75\n" \
                        "BORDER_WIDTH 1" \
                        "MARGIN 10\n" \
                        "SHOW_INTERNAL 1\n" \
                        "DATA\n"
    for taxon, col in palette_taxa_map.items():
        try:
            leaf_ranges = taxa_nodes[taxon]
        except KeyError:
            continue
        for leaves in leaf_ranges:
            colours_style_line = ['|'.join(leaves), col, taxon]
            if len(colours_style_line) > 0:
                colour_strip_text += " ".join(colours_style_line) + "\n"

    LOGGER.info("iTOL colour strip input written to " + colour_strip_file + ".\n")

    create_write_file(colour_strip_file, colour_strip_text)
    return


def clade_leaf_sides(clade_leaves, leaf_order: list) -> namedtuple:
    positions = namedtuple("positions",
                           ["left", "right"])
    indices = dict()
    for leaf in clade_leaves:
        indices[leaf_order.index(leaf)] = leaf
    p = positions
    p.left = indices[min(indices.keys())]
    p.right = indices[max(indices.keys())]
    return p


def order_taxa(taxa_to_colour: set, taxon_leaf_map: dict, leaf_order: list, method="phylo"):
    # TODO: improve this function to look at mean/median distances from the left-most taxon
    leftmost_taxon = dict()
    taxon_order = dict()
    indices = dict()  # Could just use a list but may as well keep track of the left-most node with dict
    if method == "phylo":
        for taxon, leaves in taxon_leaf_map.items():
            if taxon not in taxa_to_colour:
                continue
            for leaf_name in leaves:
                indices[int(leaf_order.index(leaf_name))] = leaf_name
            leftmost_leaf = min(indices.keys())
            leftmost_taxon[leftmost_leaf] = taxon
            indices.clear()
        # Collect the order of the taxa in the phylogeny by their left-most positions
        i = 0
        for index in sorted(leftmost_taxon):
            taxon_order[i] = leftmost_taxon[index]
            i += 1
        # Ensure all of the taxa that need to be coloured are in taxon_order
        for taxon in taxa_to_colour.difference(set(leftmost_taxon.values())):
            taxon_order[i] = taxon
            i += 1
    elif method == "alpha":
        i = 0
        for feature in sorted(taxa_to_colour):
            taxon_order[i] = feature
            i += 1
    return taxon_order


def convert_taxa_to_phenotypes(taxon_leaf_map: dict, phenotypes_map: dict) -> dict:
    """Creates a dictionary mapping feature annotation names to a list of leaf node names."""
    phenotype_leaf_map = {}
    while taxon_leaf_map:
        taxon, leaves = taxon_leaf_map.popitem()  # type: (str, list)
        phenotype = phenotypes_map[taxon]  # type: str
        try:
            phenotype_leaf_map[phenotype].extend(leaves)
        except KeyError:
            phenotype_leaf_map[phenotype] = leaves
    return phenotype_leaf_map
