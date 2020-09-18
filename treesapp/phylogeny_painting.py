#!/usr/bin/env python3

__author__ = 'Connor Morgan-Lang'

import sys
import re
import os
import logging

import ete3
from collections import namedtuple
from seaborn import color_palette

from treesapp.phylo_seq import TreeLeafReference
from treesapp.classy import TreeSAPP


class Clade:
    def __init__(self):
        self.taxon = ""
        self.colour = ""
        self.leaves = []  # Leaf names (not internal nodes)
        self.i_nodes = []  # List of all internal nodes that
        return


class PhyPainter(TreeSAPP):
    def __init__(self):
        super(PhyPainter, self).__init__("colour")
        self.style_output = ""
        self.strip_output = ""
        self.num_seqs = 0
        self.num_taxa = 0
        self.tree_leaves = list()  # List of 'TreeLeafReference's
        self.taxon_leaf_map = dict()  # taxon -> list(TreeLeafReference.number)
        self.unique_clades = dict()  # Unique numeric ID -> taxon
        self.internal_node_map = dict()  # Internal node -> child leaf names
        self.rank = ""
        self.palette = ""
        self.rank_depth = 0

    def primer(self, args) -> None:
        self.ref_pkg.f__json = args.pkg_path
        self.ref_pkg.slurp()
        self.ref_pkg.validate()

        logging.debug("\tGenerating colour palette for " + self.ref_pkg.prefix + ".\n")

        self.output_dir = args.output
        if not os.path.isdir(self.output_dir):
            try:
                os.makedirs(self.output_dir)
            except IOError:
                logging.error("Unable to create output directory '{}'.\n".format(self.output_dir))
                sys.exit(1)

        self.style_output = self.output_dir + self.ref_pkg.prefix + "_colours_style.txt"
        self.strip_output = self.output_dir + self.ref_pkg.prefix + "_colour_strip.txt"

        rank_map = self.ref_pkg.taxa_trie.accepted_ranks_depths
        if args.rank not in rank_map:
            logging.error("Rank '{}' not accepted! Please choose one of the following:\n"
                          "{}\n".format(args.rank, ", ".join(rank_map.keys())))
            sys.exit(1)
        else:
            self.rank = args.rank
            self.rank_depth = rank_map[self.rank]
        self.palette = args.palette

        self.internal_node_map = self.ref_pkg.get_internal_node_leaf_map()
        self.tree_leaves = self.ref_pkg.generate_tree_leaf_references_from_refpkg()

        return

    def check_rank_depth(self) -> None:
        """
        Depending on the lineages parsed by `treesapp create` and users, there may be some discrepancies between a rank
        and its depth (distance from the root) in a TaxonomicHierarhcy (ref_pkg.taxa_trie.accepted_ranks_depths)

        This function is meant to check whether all of the lineages begin at the same rank and therefore all lineages
        can be sliced consistently, e.g. the third position in a lineage will always represent the taxon's class

        :return: None
        """
        offsets = set()
        domain_reps = self.ref_pkg.taxa_trie.rank_representatives("domain", with_prefix=True)
        luca_reps = self.ref_pkg.taxa_trie.rank_representatives("root", with_prefix=True)
        for leaf in sorted(self.tree_leaves, key=lambda x: x.lineage):  # type: TreeLeafReference
            if not leaf.lineage:
                continue
            lineage_path = leaf.lineage.split(self.ref_pkg.taxa_trie.lin_sep)
            if len(lineage_path) < self.rank_depth-1:
                continue
            if lineage_path[0] in domain_reps:
                offsets.add(self.rank_depth - 1)
            elif lineage_path[0] in luca_reps:
                offsets.add(self.rank_depth)
            else:
                logging.debug("Unable to find root or domain name for '%s'. Skipping.\n" % leaf.lineage)
                continue

        if len(offsets) > 1:
            logging.error("Inconsistent rank positions encountered across lineages.\n")
            sys.exit(13)
        self.rank_depth = offsets.pop()

        return

    def get_clades(self) -> None:
        """
        Loads unique taxa found in the tree into a dictionary where the values are all leaves of a taxon.
        num_taxa is incremented here and represents the number of unique taxa.
        The number of keys in taxon_leaf_map represents the number of

        :return: None
        """
        clades = dict()
        truncated_lineages = dict()
        self.num_taxa = 0

        for leaf in sorted(self.tree_leaves, key=lambda x: x.lineage):  # type: TreeLeafReference
            if not leaf.lineage:
                continue
            lineage_path = leaf.lineage.split(self.ref_pkg.taxa_trie.lin_sep)
            if len(lineage_path) > self.rank_depth:
                taxon = lineage_path[self.rank_depth]
                try:
                    self.taxon_leaf_map[taxon].append(leaf.number + '_' + self.ref_pkg.prefix)
                except KeyError:
                    self.taxon_leaf_map[taxon] = [leaf.number + '_' + self.ref_pkg.prefix]
                if len(clades) == 0 and self.num_taxa == 0:
                    clades[self.num_taxa] = taxon
                elif clades[self.num_taxa] != taxon:
                    self.num_taxa += 1
                    clades[self.num_taxa] = taxon
            else:
                try:
                    truncated_lineages[leaf.lineage] += 1
                except KeyError:
                    truncated_lineages[leaf.lineage] = 1
        self.num_taxa += 1
        self.unique_clades = clades
        logging.info("{}: {} unique clades.\n"
                     "{} truncated lineages ({} sequences) will not be coloured.\n"
                     "".format(self.ref_pkg.prefix, self.num_taxa,
                               len(truncated_lineages), sum(truncated_lineages.values())))
        return

    def remove_taxon_from_colours(self, taxon: str) -> None:
        """
        At various stages during the `treesapp colour` workflow, some taxa are removed/filtered from collections
        with names of taxa that should be annotated on the phylogeny. They may be removed for multiple reasons,
        such as too few members or being a polyphyletic group.

        This function ensures that these taxa are removed from _all_ collections that list taxa so they are
        guaranteed to not be included in any annotation (i.e. iTOL-compatible colour) files.

        :param taxon: Name of a taxon to remove from the various collections
        :return: None
        """
        for n, name in self.unique_clades.items():  # type: (int, str)
            if name == taxon:
                self.unique_clades.pop(n)
                break
        self.taxon_leaf_map.pop(taxon)
        return

    def filter_unwanted_taxa(self, taxa_filter: str) -> None:
        """
        A lineage-based string filter for excluding taxa from the output colour palette.
        Even though the palette is for taxa, and not complete lineages,
        we may want to filter out a very large group of organisms, such as Eukaryota.
        This facilitates that operation so not all taxa from that larger group need to be individually named

        :param taxa_filter: A comma-separated list of taxa to remove from those that should be coloured
        :return: None
        """
        filter_str = str(taxa_filter) + " taxa that won't be coloured:\n\t"
        filtered_taxa = set()
        filtered_seqs = 0
        filter_terms = taxa_filter.split(',')
        for leaf in self.tree_leaves:  # type: TreeLeafReference
            if not leaf.lineage:
                continue
            if len(filter_terms) > 0:
                for term in filter_terms:
                    if re.search(term, leaf.lineage):
                        lineage_path = leaf.lineage.split(self.ref_pkg.taxa_trie.lin_sep)
                        taxon = lineage_path[self.rank_depth]
                        filtered_taxa.add(taxon)

        for taxon in filtered_taxa:
            self.num_taxa -= 1
            filtered_seqs += len(self.taxon_leaf_map[taxon])
            self.num_seqs -= len(self.taxon_leaf_map[taxon])
            self.remove_taxon_from_colours(taxon)

        filter_str += ", ".join(filtered_taxa) + "\n"
        filter_str += "Sequences from unwanted taxa removed\t" + str(filtered_seqs) + "\n"
        filter_str += "Unwanted taxa removed\t" + str(len(filtered_taxa)) + "\n"
        filter_str += "Remaining sequences for colouring\t" + str(self.num_seqs) + "\n"
        filter_str += "Remaining taxa for colouring\t" + str(self.num_taxa) + "\n"
        logging.info(filter_str)
        return

    def filter_rare_groups(self, min_p: float):
        filter_str = "Taxa filtered due to proportion < " + str(min_p) + ":\n"
        filtered_taxa = []
        filtered_taxa_strings = []
        for taxon in self.taxon_leaf_map:
            leaf_nodes = set(self.taxon_leaf_map[taxon])
            if len(leaf_nodes)/self.num_seqs < min_p:
                filtered_taxa_strings.append("\t" + taxon + " = " + str(round(len(leaf_nodes)/self.num_seqs, 5)))
                filtered_taxa.append(taxon)
        if not filtered_taxa_strings:
            filtered_taxa_strings.append("\tNone")
        for taxon in filtered_taxa:
            self.num_taxa -= 1
            self.num_seqs -= len(self.taxon_leaf_map[taxon])
            self.remove_taxon_from_colours(taxon)
        filter_str += "\n".join(filtered_taxa_strings) + "\n"
        filter_str += "Low-abundance taxa removed\t" + str(len(filtered_taxa)) + "\n"
        filter_str += "Remaining sequences for colouring\t" + str(self.num_seqs) + "\n"
        filter_str += "Remaining taxa for colouring\t" + str(self.num_taxa) + "\n"
        logging.info(filter_str)
        return

    def filter_polyphyletic_groups(self) -> None:
        bad_taxa = []
        for taxon in self.taxon_leaf_map:
            ancestral_node = None
            leaf_nodes = set(self.taxon_leaf_map[taxon])

            # Sort the internal nodes by the number of children (moving from shallow to deep)
            for i_node in sorted(self.internal_node_map, key=lambda n: len(self.internal_node_map[n])):
                if len(self.internal_node_map[i_node]) < len(leaf_nodes):
                    continue
                elif len(set(self.internal_node_map[i_node]).intersection(leaf_nodes)) == len(leaf_nodes):
                    ancestral_node = i_node
                    break
                else:
                    pass
            if not ancestral_node:
                logging.warning("Unable to find internal node ancestral to all children of " + taxon + ".\n")
                bad_taxa.append(taxon)
                continue
            if len(self.internal_node_map[ancestral_node]) > len(leaf_nodes):
                logging.warning(taxon + " clade was found to be polyphyletic.\n")
                bad_taxa.append(taxon)

        for taxon in bad_taxa:
            self.remove_taxon_from_colours(taxon)

        return

    def find_mono_clades(self) -> dict:
        taxa_internal_node_map = dict()
        failed_taxa = []
        for taxon in self.taxon_leaf_map:
            taxa_internal_node_map[taxon] = []
            taxon_leaves = set(self.taxon_leaf_map[taxon])
            # Sorting in descending order of number of child nodes
            for i_node in sorted(self.internal_node_map, key=lambda n: len(self.internal_node_map[n]), reverse=True):
                i_node_leaves = set(self.internal_node_map[i_node])
                if len(i_node_leaves) > len(taxon_leaves):
                    continue
                # taxon_leaves can either be a superset or match
                # Detect whether the internal node contains a set of only the leaves belonging to the taxon
                if taxon_leaves.issuperset(i_node_leaves):
                    taxa_internal_node_map[taxon].append(i_node_leaves)
                    # Pop from the list of all taxon leaves
                    for leaf in i_node_leaves:
                        taxon_leaves.remove(leaf)
                else:
                    pass
                if len(taxon_leaves) == 0:
                    break
            if taxon_leaves:
                failed_taxa.append(taxon)
        if failed_taxa:
            logging.warning("Unable to find internal nodes for all children of:\n\t" +
                            ", ".join(failed_taxa) + "\n")
        return taxa_internal_node_map


def linearize_tree_leaves(tree_string) -> list:
    """
    Parses a tree by depth-first search and returns the list

    :return: A list of internal node identifiers sorted by post-order traversal
    """
    tree_root = ete3.Tree(tree_string)
    leaf_order = []
    for node in tree_root.traverse(strategy="postorder"):
        if node.name:
            leaf_order.append(node.name)
    return leaf_order


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
        logging.error("Unable to open file '{}' for writing!\n".format(file_name))
        raise IOError

    # Write the iTOL-formatted header and node, colour and taxon fields
    file_handler.write(text)
    file_handler.close()
    return


def get_colours(painter: PhyPainter):
    clades = set()
    for num in painter.unique_clades:
        clades.add(painter.unique_clades[num])

    if len(clades) <= 1:
        logging.error("{} unique clade(s) identified at rank '{}'."
                      " Consider changing rank to something more specific.\n".format(len(clades), painter.rank))
        sys.exit(4)
    else:
        logging.debug("Identified {} unique clades at {} rank.\n".format(len(clades), painter.rank))

    n_clade = len(clades)
    try:
        colours = color_palette(painter.palette, n_clade).as_hex()
        logging.info("Using palette '{}' for styling.\n".format(painter.palette))
    except ValueError as e:
        logging.error(str(e) + "\n")
        sys.exit(4)

    if len(colours) != n_clade:
        logging.error("Bad colour parsing! len(colours) != number of clades.\n")
        sys.exit(8)
    return colours


def map_colours_to_taxa(taxa_order, colours):
    """
    Function for mapping

    :param taxa_order: Dictionary indexed by numerical order of the taxa
    :param colours:
    :return:
    """
    palette_taxa_map = dict()
    alpha_sort_i = 0

    if len(taxa_order) != len(colours):
        logging.error("Number of taxa (%d) and colours (%d) are not equal!\n" % (len(taxa_order), len(colours)))
        sys.exit(7)

    for taxon_num in sorted(taxa_order, key=int):
        taxon = taxa_order[taxon_num]
        palette_taxa_map[taxon] = colours[alpha_sort_i]
        alpha_sort_i += 1

    return palette_taxa_map


def write_colours_styles(taxa_colours: PhyPainter, palette_taxa_map: dict) -> None:
    colourless = set()

    colours_style_string = "TREE_COLORS\nSEPARATOR TAB\nDATA\n"
    for taxon in palette_taxa_map:
        col = palette_taxa_map[taxon]
        for leaf_num in taxa_colours.taxon_leaf_map[taxon]:
            colours_style_line = [str(leaf_num), "range", col, taxon]
            if len(colours_style_line) > 0:
                colours_style_string += "\t".join(colours_style_line) + "\n"
    # TODO: Find new way to track the uncoloured leaves

    logging.info("Output is available in " + taxa_colours.style_output + ".\n")
    logging.debug("{} lineages were not assigned to a colour:\n\t{}\n".format(len(colourless), "\n\t".join(colourless)))

    create_write_file(taxa_colours.style_output, colours_style_string)

    return


def write_colour_strip(taxa_nodes: dict, palette_taxa_map: dict, colour_strip_file: str) -> None:
    colour_strip_text = "DATASET_COLORSTRIP\nSEPARATOR SPACE\nDATASET_LABEL clade_colours\n" +\
                        "STRIP_WIDTH 75\nMARGIN 0\nSHOW_INTERNAL 1\nDATA\n"
    for taxon in palette_taxa_map:
        col = palette_taxa_map[taxon]
        for leaf_range in taxa_nodes[taxon]:
            colours_style_line = ['|'.join(leaf_range), col, taxon]
            if len(colours_style_line) > 0:
                colour_strip_text += " ".join(colours_style_line) + "\n"

    logging.info("Output is available in " + colour_strip_file + ".\n")

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


def order_taxa(taxon_leaf_map: dict, leaf_order: list):
    leftmost_taxon = dict()
    taxon_order = dict()
    indices = dict()  # Could just use a list but may as well keep track of the left-most node with dict
    for taxon in taxon_leaf_map:
        for leaf_name in taxon_leaf_map[taxon]:
            indices[int(leaf_order.index(leaf_name))] = leaf_name
        leftmost_leaf = min(indices.keys())
        leftmost_taxon[leftmost_leaf] = taxon
        indices.clear()
    i = 0
    for index in sorted(leftmost_taxon):
        taxon_order[i] = leftmost_taxon[index]
        i += 1
    return taxon_order
