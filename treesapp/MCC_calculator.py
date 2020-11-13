#!/usr/bin/env python3

__author__ = 'Connor Morgan-Lang'

import os
import shutil
import sys
import logging
from time import sleep
from glob import glob
from numpy import sqrt

from ete3 import Tree

from treesapp import file_parsers
from treesapp.phylo_seq import assignments_to_treesaps
from treesapp.refpkg import ReferencePackage
from treesapp.commands import assign
from treesapp.fasta import get_headers, register_headers
from treesapp.external_command_interface import launch_write_command
from treesapp.classy import prep_logging, get_header_info
from treesapp.entrez_utils import EntrezRecord, get_multiple_lineages
from treesapp.lca_calculations import compute_taxonomic_distance, optimal_taxonomic_assignment, grab_graftm_taxa
from treesapp.treesapp_args import TreeSAPPArgumentParser
from treesapp.taxonomic_hierarchy import TaxonomicHierarchy
from treesapp.training_utils import bin_headers, QuerySequence


class ConfusionTest:
    def __init__(self, gene_list):
        self._MAX_TAX_DIST = -1
        self.ref_packages = {key: ReferencePackage() for key in gene_list}
        self.fn = {key: [] for key in gene_list}
        self.fp = {key: set() for key in gene_list}
        self.tp = {key: [] for key in gene_list}  # This will be a list of QuerySequence instances
        self.redundant_queries = set()
        self.tax_lineage_map = {}
        self.tp_lineage_map = {}
        self.dist_wise_tp = {}
        self.header_registry = {}
        self.entrez_query_dict = {}
        self.all_queries = set()
        self.rank_depth_map = {"Domain": 1, "Phylum": 2, "Class": 3, "Order": 4, "Family": 5, "Genus": 6, "Species": 7}
        self.classification_table = ""
        self.treesapp_dir = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + os.sep
        self.refpkg_dir = self.treesapp_dir + 'data' + os.sep
        self.data_dir = ""

    def get_info(self, verbose=False):
        info_string = "\nReference packages being tested:\n"
        info_string += "\t" + ", ".join(list(self.ref_packages.keys())) + "\n"

        self.check_dist()
        info_string += "Stats based on taxonomic distance < " + str(self._MAX_TAX_DIST) + "\n"

        if len(self.all_queries) > 0:
            info_string += "{} query sequences being used for testing.\n".format(len(self.all_queries))

        if verbose:
            for refpkg in sorted(self.ref_packages):
                info_string += self.marker_classification_summary(refpkg)

        return info_string

    def map_data(self, output_dir, tool: str):
        if tool == "treesapp":
            self.data_dir = output_dir + os.sep + "TreeSAPP_output" + os.sep
            self.classification_table = self.data_dir + "final_outputs" + os.sep + "marker_contig_map.tsv"
        elif tool == "graftm":
            self.data_dir = output_dir + os.sep + "GraftM_output" + os.sep
        elif tool == "diamond":
            self.data_dir = output_dir + os.sep + "DIAMOND_output" + os.sep
        else:
            logging.error("Unrecognized tool: " + tool + "\n")
            sys.exit(1)
        return

    def summarise_rank_coverage(self, rank_coverage: dict) -> str:
        summary_string = "Rank\tClassifications\n"
        for rank in sorted(self.rank_depth_map, key=lambda x: self.rank_depth_map[x]):
            try:
                summary_string += rank + "\t" + str(rank_coverage[self.rank_depth_map[rank]]) + "\n"
            except KeyError:
                continue
        return summary_string

    def marker_classification_summary(self, refpkg_name):
        """
        Provide a classification summary for a specific marker gene, refpkg_name

        :param refpkg_name:
        :return: A string summarizing the classification performance of a single marker/refpkg_name
        """
        self.check_dist()
        self.check_refpkg_name(refpkg_name)
        tp, remain = self.get_true_positives_at_dist(refpkg_name)
        summary_string = "\nSummary for reference package '" + str(refpkg_name) + "':\n"
        summary_string += "\tTrue positives\t\t" + str(len(self.tp[refpkg_name])) + "\n"
        summary_string += "Stats based on taxonomic distance <" + str(self._MAX_TAX_DIST) + ":\n"
        summary_string += "\tTrue positives\t\t" + str(len(tp)) + "\n"
        summary_string += "\tFalse positives\t\t" + str(len(self.get_false_positives(refpkg_name)) + len(remain)) + "\n"
        summary_string += "\tFalse negatives\t\t" + str(len(self.get_false_negatives(refpkg_name))) + "\n"
        summary_string += "\tTrue negatives\t\t" + str(self.get_true_negatives(refpkg_name)) + "\n"
        return summary_string

    def populate_tax_lineage_map(self, entrez_record_dict: dict) -> None:
        """
        Used for pulling the NCBI taxid from a sequence name (header). This is most effective for handling
        EggNOG-formatted headers where the NCBI taxid is the first part, followed by a period, then the organism name

        :param entrez_record_dict: A dictionary mapping unique numerical TreeSAPP identifiers to EntrezRecord instances
        :return: None
        """
        for ts_id in entrez_record_dict:  # type: str
            e_record = entrez_record_dict[ts_id]  # type: EntrezRecord
            if not e_record.description:
                logging.error("No description for the following EntrezRecord:\n{}\n".format(e_record.get_info()))
                sys.exit(3)

            # Only make Entrez records for new NCBI taxonomy IDs
            if e_record.description not in self.tax_lineage_map:
                if not e_record.lineage:
                    logging.warning("Lineage information unavailable for accession '{}'\n".format(e_record.accession))
                    self.tax_lineage_map[e_record.description] = ''
                else:
                    self.tax_lineage_map[e_record.description] = e_record.lineage

        return

    def populate_true_positive_lineage_map(self) -> None:
        for refpkg_code in self.ref_packages:
            self.tp_lineage_map[refpkg_code] = {tp_inst.ref_name: tp_inst.true_lineage
                                                for tp_inst in self.tp[refpkg_code]}
        return

    def generate_entrez_queries(self) -> None:
        entrez_record_dict = get_header_info(self.header_registry)
        for index, e_record in entrez_record_dict.items():  # type: EntrezRecord
            self.entrez_query_dict[e_record.description] = e_record
        return

    def retrieve_lineages(self) -> None:
        """
        Determines the format of the query sequence header to extract the accession and/or NCBI taxid then
        queries the Entrez database using these information to collect the taxonomic lineage for each unique NCBI taxid
        NCBI taxid and lineage information are stored in self.tax_lineage_map

        :return: None
        """
        # Gather the unique taxonomy IDs and store in EntrezRecord instances
        t_hierarchy = TaxonomicHierarchy()
        entrez_records = [self.entrez_query_dict[index] for index in self.entrez_query_dict]
        # Query the Entrez database for these unique taxonomy IDs
        get_multiple_lineages(entrez_records, t_hierarchy, "prot")

        for e_record in entrez_records:  # type: EntrezRecord
            if not e_record.lineage:
                logging.warning("Lineage information unavailable for taxonomy ID '" + e_record.ncbi_tax + "'\n")
                self.tax_lineage_map[e_record.ncbi_tax] = ''
            else:
                self.tax_lineage_map[e_record.ncbi_tax] = "r__Root; " + e_record.lineage
        return

    def map_true_lineages(self) -> None:
        for refpkg_name in self.ref_packages:
            if refpkg_name in self.tp:
                map_lineages(self.tp[refpkg_name], self.tax_lineage_map)
            if refpkg_name in self.fn:
                map_lineages(self.fn[refpkg_name], self.tax_lineage_map)
        return

    def summarise_reference_taxa(self, taxa_file: str, classification_file: str, rank="Phylum"):
        lineage_list = []
        info_string = "RefPkg\tName\tTaxDist\tClassified\tTrueLineage\tAssignedLineage\tOptimalLineage\n"
        for marker in self.tp:
            for tp_inst in self.tp[marker]:  # type: QuerySequence
                lineage_list.append(tp_inst.true_lineage)
                info_string += "\t".join([marker, tp_inst.place_name, str(tp_inst.tax_dist), "True",
                                          tp_inst.true_lineage, tp_inst.assigned_lineage,
                                          tp_inst.optimal_lineage]) + "\n"
        for marker in self.fn:
            for tp_inst in self.fn[marker]:  # type: QuerySequence
                # TODO: Support headers from databases other than EggNOG
                tax_id = tp_inst.ncbi_tax
                try:
                    true_lineage = self.tax_lineage_map[tax_id]
                except KeyError:
                    continue
                lineage_list.append(true_lineage)
                info_string += "\t".join([marker, tp_inst.place_name, "NA", "False", true_lineage, "NA", "NA"]) + "\n"

        with open(classification_file, 'w') as info_handler:
            info_handler.write(info_string)
        logging.debug("Taxonomic lineage distribution of " + str(len(lineage_list)) + " 'true' reference sequences.\n")
        input_taxa_dist = summarize_taxonomy(taxa_list=lineage_list,
                                             rank=rank,
                                             rank_depth_map=self.rank_depth_map)
        write_dict_to_table(input_taxa_dist, taxa_file, "\t")
        summary_string = ""
        for summary_lineage in sorted(input_taxa_dist, key=lambda x: input_taxa_dist[x]):
            summary_string += "\t'" + summary_lineage + "'\t" + str(input_taxa_dist[summary_lineage]) + "\n"
        logging.debug(summary_string)
        return

    def bin_true_positives_by_taxdist(self):
        """
        Defines the number of true positives at each taxonomic distance x where 0 <= x <= 7,
        since there are 7 ranks in the NCBI taxonomic hierarchy.
        All sequences correctly classified (at the gene level) are assigned a taxonomic distance,
        so the sum of dist_wise_tp[x] for all x will equal the number of all true positives.

        :return: None
        """
        for marker in self.tp:
            refpkg = self.ref_packages[marker]  # type: ReferencePackage
            refpkg.taxa_trie.trie_check()
            if not refpkg.taxa_trie.rooted:
                refpkg.taxa_trie.root_domains(root=refpkg.taxa_trie.find_root_taxon())
            self.dist_wise_tp[marker] = dict()
            for tp_inst in self.tp[marker]:  # type: QuerySequence
                # Find the optimal taxonomic assignment
                tp_inst.true_lineage = refpkg.taxa_trie.clean_lineage_string(tp_inst.true_lineage)
                optimal_taxon = optimal_taxonomic_assignment(refpkg.taxa_trie.trie, tp_inst.true_lineage)
                if not optimal_taxon:
                    logging.debug("Optimal taxonomic assignment '{}' for {}"
                                  " not found in reference hierarchy.\n".format(tp_inst.true_lineage, tp_inst.place_name))
                    continue
                tp_inst.optimal_lineage = optimal_taxon
                tp_inst.tax_dist, status = compute_taxonomic_distance(tp_inst.assigned_lineage, optimal_taxon)
                if status > 0:
                    logging.debug("Lineages didn't converge between:\n"
                                  "'{}' and '{}' (taxid: {})\n".format(tp_inst.assigned_lineage,
                                                                       optimal_taxon, tp_inst.ncbi_tax))
                try:
                    self.dist_wise_tp[marker][tp_inst.tax_dist].append(tp_inst.place_name)
                except KeyError:
                    self.dist_wise_tp[marker][tp_inst.tax_dist] = [tp_inst.place_name]
        return

    def enumerate_optimal_rank_distances(self) -> dict:
        """
        Calculates the number of query sequences representing an excluded rank.

        For each query sequence in self.tp (so just the true positives), the optimal taxonomy is identified based on
        the query's true taxonomic lineage and all the taxonomic lineages in the ReferencePackage.taxa_trie.
        The optimal taxonomic assignment's depth (i.e. numerical description of a rank where domain = 1,
         phylum = 2, etc.) is found by determining the length of the list after splitting the lineage string by '; '.

        :return: Dictionary mapping the numerical description of a taxonomic rank (i.e. depth) to the number of queries
        that represent that rank.
        """
        rank_representation = dict()
        for marker in self.tp:
            refpkg = self.ref_packages[marker]  # type: ReferencePackage
            refpkg.taxa_trie.trie_check()

            for tp_inst in self.tp[marker]:  # type: QuerySequence
                # Find the optimal taxonomic assignment
                optimal_taxon = optimal_taxonomic_assignment(refpkg.taxa_trie.trie, tp_inst.true_lineage)
                if not optimal_taxon:
                    logging.debug("Optimal taxonomic assignment '{}' for {} "
                                  "not found in reference hierarchy.\n".format(tp_inst.true_lineage, tp_inst.place_name))
                    continue
                try:
                    rank_representation[len(optimal_taxon.split("; "))] += 1
                except KeyError:
                    rank_representation[len(optimal_taxon.split("; "))] = 1
        return rank_representation

    def check_dist(self):
        if self._MAX_TAX_DIST < 0:
            logging.error("ConfusionTest's _MAX_TAX_DIST has yet to be set.\n")
            sys.exit(5)
        return

    def check_refpkg_name(self, refpkg_name):
        if refpkg_name not in self.ref_packages:
            logging.error(refpkg_name + " is not found in the names of markers to be tested.\n")
            sys.exit(9)
        return

    def get_true_positives_at_dist(self, refpkg_name=None):
        """
        Calculates the sum of all true positives at a specified maximum taxonomic distance and less.
        Sequences classified at a distance greater than self._MAX_TAX_DIST are counted as false negatives,
        since it is as if the sequences were not classified at all.

        :return: The sum of true positives at taxonomic distance <= max_distance
        """
        self.check_dist()
        all_tp_headers = set()
        remainder_headers = set()
        if refpkg_name:
            marker_set = [refpkg_name]
        else:
            marker_set = self.dist_wise_tp
        for ref_name in marker_set:
            for tax_dist in sorted(self.dist_wise_tp[ref_name], key=int):  # type: int
                if tax_dist <= self._MAX_TAX_DIST:
                    all_tp_headers.update(set(self.dist_wise_tp[ref_name][tax_dist]))
                else:
                    remainder_headers.update(set(self.dist_wise_tp[ref_name][tax_dist]))
        return all_tp_headers, remainder_headers

    def get_true_negatives(self, refpkg_name=None):
        if refpkg_name:
            marker_set = [refpkg_name]
        else:
            marker_set = self.ref_packages
        unique_fn = set()
        unique_fp = set()
        unique_tp = set()
        for marker in marker_set:
            unique_fn.update(qseq.place_name for qseq in self.fn[marker])
            unique_fp.update(qseq.place_name for qseq in self.fp[marker])
            unique_tp.update(qseq.place_name for qseq in self.tp[marker])
        return len(self.all_queries) - sum([len(unique_fn), len(unique_fp), len(unique_tp)])

    def get_false_positives(self, refpkg_name=None):
        if refpkg_name:
            return self.fp[refpkg_name]
        else:
            unique_fp = set()
            for marker in self.ref_packages:
                unique_fp.update(qseq.place_name for qseq in self.fp[marker])
            return unique_fp

    def get_false_negatives(self, refpkg_name=None):
        if refpkg_name:
            return self.fn[refpkg_name]
        else:
            unique_fn = set()
            for marker in self.ref_packages:
                # Remove sequences that were classified by a homologous marker
                unique_fn.update(qseq.place_name for qseq in self.fn[marker])
            return unique_fn

    def bin_headers(self, assignments: dict, annot_map: dict) -> None:
        binned_tp, binned_fp, binned_fn = bin_headers(assignments, annot_map, self.entrez_query_dict, self.ref_packages)

        self.tp.update(binned_tp)
        self.fp.update(binned_fp)
        self.fn.update(binned_fn)

        return

    def validate_false_positives(self):
        # Get all the original Orthologous Group (OG) headers for sequences classified as TP or FN
        tp_names = set()
        for marker in list(self.ref_packages.keys()):
            tp_names = tp_names.union(set([pquery.place_name for pquery in self.tp[marker]]))
            tp_names = tp_names.union(set([pquery.place_name for pquery in self.fn[marker]]))

        for marker in self.fp:
            validated_fp = set()
            for qseq in self.fp[marker]:  # type: QuerySequence
                seq_name = qseq.place_name
                if seq_name[0] == '>':
                    seq_name = seq_name[1:]
                if seq_name not in tp_names:
                    validated_fp.add(qseq)
            self.fp[marker] = self.all_queries.intersection(validated_fp)
        return

    def validate_false_negatives(self, refpkg_dbname_dict: dict):
        # Invert the dictionary
        refpkg_og_map = dict()
        for refpkg_name in refpkg_dbname_dict:
            ogs = refpkg_dbname_dict[refpkg_name]
            for og in ogs:
                try:
                    refpkg_og_map[og].append(refpkg_name)
                except KeyError:
                    refpkg_og_map[og] = [refpkg_name]

        for marker in self.fn:
            if marker in refpkg_dbname_dict:
                homologous_tp_names = set()
                for og in refpkg_dbname_dict[marker]:
                    for homolgous_marker in refpkg_og_map[og]:
                        if homolgous_marker != marker:
                            homologous_tp_names.update(set(qs.place_name for qs in self.tp[homolgous_marker]))
                validated_fns = set()
                for qseq in self.fn[marker]:  # type: QuerySequence
                    if qseq.place_name not in homologous_tp_names:
                        validated_fns.add(qseq)
                # Remove all sequences from this marker's false negatives that are found in a homologous TP set
                self.fn[marker] = self.all_queries.intersection(validated_fns)
        return

    def summarise_type_one_placements(self, classification_lines) -> str:
        """
        Its nice to understand the placement distances and where on the phylogeny false positives were inserted.
        First, figure out which headers are false positives, then

        :param classification_lines:
        :return: Summary of the distances and internal nodes for each of the false positives
        """
        # Read internal node maps for each refpkg
        internal_nodes_dict = dict()
        summary_dict = {"leaves": {}, "LWR": {}, "distances": {}}
        refpkg_map = dict()
        hmm_lengths = dict()
        for name in self.ref_packages:
            refpkg = self.ref_packages[name]  # type: ReferencePackage
            refpkg_map[refpkg.prefix] = name
            try:
                internal_nodes_dict[refpkg.prefix] = refpkg.get_internal_node_leaf_map()
            except IndexError:
                logging.error("Unable to read tree for reference package %s from '%s'.\n" % (name, refpkg.tree))
                sys.exit(3)
            hmm_lengths[refpkg.prefix] = refpkg.profile_length
        #
        for fields in classification_lines:
            _, header, refpkg_name, start_pos, end_pos, _, _, i_node, e_val, lwr, evo_dist, dists = fields
            if header in self.fp[refpkg_map[refpkg_name]]:
                distal, pendant, avg = [round(float(x), 3) for x in dists.split(',')]
                # hmm_perc = round((int(length)*100)/hmm_lengths[refpkg], 0)
                descendents = len(internal_nodes_dict[refpkg_name][i_node])
                if descendents not in summary_dict["leaves"]:
                    summary_dict["leaves"][descendents] = 0
                summary_dict["leaves"][descendents] += 1

                lwr_bin = round(float(lwr), 2)
                if lwr_bin not in summary_dict["LWR"]:
                    summary_dict["LWR"][lwr_bin] = 0
                summary_dict["LWR"][lwr_bin] += 1

                dist_bin = round(float(pendant), 1)
                if dist_bin not in summary_dict["distances"]:
                    summary_dict["distances"][dist_bin] = 0
                summary_dict["distances"][dist_bin] += 1

        # Convert the dictionary into a human-readable string
        summary_str = ""
        for measure in summary_dict:
            summary_str += "Summary of false positive %s:\n" % measure
            for n in sorted(summary_dict[measure], key=float):
                summary_str += "\t" + str(n) + "\t" + str(summary_dict[measure][n]) + "\n"
        return summary_str

    def summarize_type_two_taxa(self, rank="Phylum"):
        lineage_list = []
        summary_string = ""
        for marker in self.fn:
            if len(self.fn[marker]) == 0:
                continue
            summary_string += "Number of false negatives for '" + marker + "':\n"
            for qseq in self.fn[marker]:  # type: QuerySequence
                try:
                    lineage_list.append(qseq.true_lineage)
                except KeyError:
                    continue
            lineage_count_dict = summarize_taxonomy(lineage_list, rank, self.rank_depth_map)
            for summary_lineage in sorted(lineage_count_dict):
                summary_string += "\t'" + summary_lineage + "'\t" + str(lineage_count_dict[summary_lineage]) + "\n"
            lineage_count_dict.clear()
            lineage_list.clear()

        return summary_string

    def true_positive_taxonomic_summary(self, rank="Phylum", aggregate=False):
        lineage_list = []
        summary_string = ""
        for marker in self.tp:
            if len(self.tp[marker]) == 0:
                continue
            lineage_list += [tp_inst.true_lineage for tp_inst in self.tp[marker]]

            if aggregate is False:
                lineage_count_dict = summarize_taxonomy(lineage_list, rank, self.rank_depth_map)
                summary_string += "Taxonomic distribution of true positives for '" + marker + "':\n"
                for summary_lineage in sorted(lineage_count_dict, key=lambda x: lineage_count_dict[x]):
                    summary_string += "\t'" + summary_lineage + "'\t" + str(lineage_count_dict[summary_lineage]) + "\n"
                lineage_count_dict.clear()
                lineage_list.clear()
        if aggregate:
            lineage_count_dict = summarize_taxonomy(lineage_list, rank, self.rank_depth_map)
            summary_string += "Taxonomic distribution of true positives:\n"
            for summary_lineage in sorted(lineage_count_dict, key=lambda x: lineage_count_dict[x]):
                summary_string += "\t'" + summary_lineage + "'\t" + str(lineage_count_dict[summary_lineage]) + "\n"
            lineage_count_dict.clear()

        return summary_string


def map_lineages(qseq_collection: set, tax_lineage_map: dict) -> None:
    for qseq in qseq_collection:  # type: QuerySequence
        try:
            qseq.true_lineage = tax_lineage_map[qseq.ncbi_tax]
        except KeyError:
            continue
    return


def summarize_taxonomy(taxa_list, rank, rank_depth_map=None):
    """
    Given a list of taxonomic lineages and a taxonomic rank for which to summarise at
    it will count the number of instances for each taxon at the desired rank.
    E.g. {''}

    :return: A dictionary mapping a taxon to the number of representatives seen
    """
    if not rank_depth_map:
        rank_depth_map = {"Domain": 1, "Phylum": 2, "Class": 3, "Order": 4, "Family": 5, "Genus": 6, "Species": 7}

    try:
        depth = rank_depth_map[rank] + 1
    except KeyError:
        logging.error("Rank '{}' not present in rank-depth map:\n"
                      "{}\n".format(rank, rank_depth_map))
        sys.exit(3)

    taxa_census = dict()
    empty = 0
    acc = 0
    for taxon in taxa_list:
        if not taxon:
            empty += 1
            continue
        taxon_path = taxon.split("; ")
        if depth > len(taxon_path):
            summary_lineage = taxon
        else:
            summary_lineage = "; ".join(taxon_path[:depth])
        try:
            taxa_census[summary_lineage] += 1
        except KeyError:
            taxa_census[summary_lineage] = 1
        acc += 1
    logging.debug(str(empty) + " empty lineages encountered.\n")

    return taxa_census


def write_dict_to_table(data_dict, file_name, sep="\t"):
    """
    Basic function for writing the key, value pairs into a file with a specified separator

    :param data_dict: A dictionary object
    :param file_name: Path to a file for the dictionary to be written to
    :param sep: The separator to use. Tabs by default
    :return: None
    """
    data_string = ""
    for key in sorted(data_dict):
        data_string += sep.join([str(key), str(data_dict[key])]) + "\n"
    with open(file_name, 'w') as data_out:
        data_out.write(data_string)

    return


def get_arguments(sys_args):
    parser = TreeSAPPArgumentParser(description="Calculate Matthews' correlation coefficient from classifications")
    ts = parser.add_argument_group("TreeSAPP options")

    parser.add_io()
    parser.add_refpkg_opt()
    parser.add_refpkg_targets()
    parser.add_seq_params()
    parser.add_compute_miscellany()
    parser.add_pplace_params()
    parser.add_annot_map(required=True)

    parser.optopt.add_argument("--tool", default="treesapp", required=False,
                               choices=["treesapp", "graftm", "diamond"],
                               help="Classify using one of the tools: treesapp [DEFAULT], graftm, or diamond.")

    ts.add_argument("--svm", default=False, required=False, action="store_true",
                    help="Uses the support vector machine (SVM) classification filter. "
                         "WARNING: Unless you *really* know your refpkg, you probably don't want this.")
    ts.add_argument("-s", "--stringency", choices=["relaxed", "strict"], default="relaxed", required=False,
                    help="HMM-threshold mode affects the number of query sequences that advance [DEFAULT = relaxed]")

    args = parser.parse_args(sys_args)

    if args.output[-1] != os.sep:
        args.output += os.sep
    return args


def internal_node_leaf_map(tree: str) -> dict:
    """
    Loads a Newick-formatted tree with internal nodes into a dictionary of
    all internal nodes (keys) and a list of child leaves (values).

    :param tree: Path to a Newick tree file.
    :return: Dictionary of all internal nodes (keys) and a list of child leaves (values)
    """
    node_map = dict()
    x = 0
    ete_tree = Tree(tree)
    for node in ete_tree.traverse(strategy="postorder"):
        node_map[str(x)] = node.get_leaf_names()
        x += 1

    return node_map


def validate_command(args, sys_args):
    logging.debug("Command used:\n" + ' '.join(sys_args) + "\n")

    if not os.path.isdir(args.output):
        os.mkdir(args.output)

    args.treesapp = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + os.sep
    if args.tool == "treesapp" and not args.refpkg_dir:
        args.refpkg_dir = args.treesapp + "data" + os.sep

    if sys.version_info < (2, 9):
        logging.error("Python version '" + str(sys.version_info) + "' not supported.\n")
        sys.exit(3)

    if args.refpkg_dir and args.refpkg_dir[-1] != os.sep:
        args.refpkg_dir += os.sep
    if args.tool in ["diamond", "graftm"]:
        if not args.refpkg_dir:
            logging.error(args.tool + " specified but a GraftM reference package directory was not provided.\n")
            sys.exit(17)
        elif not os.path.isdir(args.refpkg_dir):
            logging.error(args.refpkg_dir + " GraftM reference package directory does not exist!\n")
            sys.exit(17)
        elif len(glob(args.refpkg_dir + "*gpkg")) == 0:
            logging.error("No GraftM reference packages found in " + args.refpkg_dir + ".\n")
            sys.exit(17)

    if args.targets:
        args.targets = args.targets.split(',')
    else:
        args.targets = []

    return


def calculate_matthews_correlation_coefficient(tp: int, fp: int, fn: int, tn: int):
    numerator = float((tp * tn) - (fp * fn))
    denominator = float((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    if numerator == 0 or denominator == 0:
        return 0
    else:
        return round(numerator/sqrt(denominator), 3)


def check_previous_output(output_dir: str, files: list, overwrite=False) -> None:
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    elif overwrite:
        if os.path.isdir(output_dir):
            logging.warning("Overwriting directory '{}' in 5 seconds. Press Ctrl-C to cancel.\n".format(output_dir))
            sleep(5)
            shutil.rmtree(output_dir)
        os.mkdir(output_dir)
        for f in files:
            if os.path.isfile(f):
                os.remove(f)
    return


def filter_redundant_og(query_og_map: dict) -> set:
    redundants = set()
    ortho_counts = {}
    ortho_filters = {}
    for query_name, ortho_groups in query_og_map.items():
        for og in ortho_groups:
            try:
                ortho_counts[og] += 1
            except KeyError:
                ortho_counts[og] = 1
        if len(ortho_groups) > 1:
            for og in ortho_groups:
                try:
                    ortho_filters[og] += 1
                except KeyError:
                    ortho_filters[og] = 1
            redundants.add(query_name)

    if redundants:
        logging.warning("{}/{} query sequences were annotated as multiple reference packages and have been removed.\n"
                        "".format(len(redundants), len(query_og_map)))
        for og in ortho_filters:
            logging.info("{}% of annotated sequences were filtered for '{}'.\n"
                         "".format(round((ortho_filters[og]*100)/ortho_counts[og], 1), og))

    for query_name in redundants:
        query_og_map.pop(query_name)

    return redundants


def create_refpkg_query_map(query_name_map: dict, refpkg_dict: dict) -> dict:
    refpkg_query_map = {}
    tmp_map = {}
    for rp_name, rp in refpkg_dict.items():  # type: (str, ReferencePackage)
        try:
            tmp_map[rp.refpkg_code].add(rp_name)
        except KeyError:
            tmp_map[rp.refpkg_code] = {rp_name}
    for qname, ref_codes in query_name_map.items():  # type: (str, str)
        for rp_code in ref_codes:
            for rp_name in tmp_map[rp_code]:
                try:
                    refpkg_query_map[rp_name].add(qname)
                except KeyError:
                    refpkg_query_map[rp_name] = {qname}
    return refpkg_query_map


def match_queries_to_refpkgs(query_name_map: dict, refpkg_dict: dict) -> (dict, dict):
    """
    Used for removing any query sequence names from the query name map if the reference package they should be assigned
    is not present in the refpkg_dict.

    :param query_name_map: A dictionary mapping query sequence names to either ReferencePackage
    'prefix' or 'refpkg_code' values.
    :param refpkg_dict: A dictionary mapping reference package prefixes to their respective ReferencePackage instance
    :return: None
   """
    # Create the dictionaries for rapid look-up
    prefix_to_code = {refpkg.prefix: refpkg.refpkg_code for name, refpkg in refpkg_dict.items()}
    code_to_prefix = {refpkg.refpkg_code: set() for name, refpkg in refpkg_dict.items()}
    for rp_code in code_to_prefix:
        for prefix, refpkg in refpkg_dict.items():  # type: (str, ReferencePackage)
            if rp_code == refpkg.refpkg_code:
                code_to_prefix[rp_code].add(prefix)
    query_to_prefix = {}
    query_to_code = {}
    unmapped_queries = []
    missing_matches = []
    for query, matches in query_name_map.items():  # type: (str, set)
        query_to_prefix[query] = set()
        query_to_code[query] = set()
        for match in matches:  # type: str
            if match in prefix_to_code:  # Match must be a ReferencePackage.prefix
                query_to_prefix[query].add(match)
                query_to_code[query].add(prefix_to_code[match])
            elif match in code_to_prefix:  # Match must be a ReferencePackage.refpkg_code
                query_to_code[query].add(match)
                query_to_prefix[query].update(code_to_prefix[match])
            else:
                unmapped_queries.append(query)
                missing_matches.append(match)
    # Remove the queries that don't have a valid reference package name
    for qname in unmapped_queries:
        query_to_prefix.pop(qname)
        query_to_code.pop(qname)

    if missing_matches:
        logging.info("{} unique reference packages ({} total) referred to in the annotation file were not found"
                     " in the reference package list.\n".format(len(set(missing_matches)), len(missing_matches)))
        logging.debug("Unique annotation names that couldn't be mapped to reference packages:\n{}\n"
                      "".format(set(missing_matches)))

    return query_to_prefix, query_to_code


def mcc_calculator(sys_args):
    args = get_arguments(sys_args)

    summary_rank = "Phylum"
    output_prefix = args.output + os.path.splitext(os.path.basename(args.input))[0]
    log_name = output_prefix + "_MCC_log.txt"
    mcc_file = output_prefix + "_MCC_table.tsv"
    taxa_dist_output = output_prefix + '_' + summary_rank + "_dist.tsv"
    classification_info_output = output_prefix + "_classifications.tsv"

    # Instantiate the logging instance and write the log
    prep_logging(log_name, args.verbose)

    logging.info("\n##\t\t\tBeginning Matthews Correlation Coefficient analysis\t\t\t##\n")
    validate_command(args, sys.argv)

    ##
    # Read the file mapping reference package name to the database annotations
    ##
    query_name_dict = file_parsers.read_annotation_mapping_file(args.annot_map)
    refpkg_dict = file_parsers.gather_ref_packages(args.refpkg_dir, args.targets)
    test_obj = ConfusionTest(refpkg_dict.keys())
    test_obj.map_data(output_dir=args.output, tool=args.tool)

    # Remove outputs from past runs if overwriting
    check_previous_output(output_dir=test_obj.data_dir, files=[mcc_file, taxa_dist_output, classification_info_output],
                          overwrite=args.overwrite)

    # Remove the query names that are mapped to multiple OGs or reference packages
    test_obj.redundant_queries = filter_redundant_og(query_name_dict)
    query_to_prefix, query_to_code = match_queries_to_refpkgs(query_name_dict, refpkg_dict)
    if len(query_to_code) == 0 or len(query_to_prefix) == 0:
        logging.error("Matching reference package annotations to queries failed.\n")
        sys.exit(13)

    refpkg_name_query_map = create_refpkg_query_map(query_to_code, refpkg_dict)

    # Creates the dictionary of Header instances
    test_obj.header_registry = register_headers(list(query_to_code.keys()))
    test_obj.generate_entrez_queries()
    # Downloads lineage information for the header instances
    test_obj.retrieve_lineages()

    ##
    # Load the taxonomic trie for each reference package
    ##
    if args.tool == "treesapp":
        test_obj.ref_packages = refpkg_dict
    elif args.tool in ["graftm", "diamond"]:
        for gpkg in glob(args.refpkg_dir + "*gpkg"):
            gpkg_name = str(os.path.basename(gpkg).split('.')[0])
            if gpkg_name in test_obj.ref_packages:
                try:
                    tax_ids_file = glob(os.path.join(gpkg, gpkg_name + ".gpkg.refpkg",
                                                     gpkg_name + "*taxonomy.csv")).pop()
                    test_obj.ref_packages[gpkg_name].taxa_trie = grab_graftm_taxa(tax_ids_file)
                except IndexError:
                    logging.warning("No GraftM taxonomy file found for {}. Is this gpkg complete?\n".format(gpkg_name))
    else:
        logging.error("Unrecognized tool name '{}'\n".format(args.tool))
        sys.exit(3)

    ##
    # Run the specified taxonomic analysis tool and collect the classifications
    ##
    assignments = {}
    test_fa_prefix = '.'.join(os.path.basename(args.input).split('.')[:-1])
    classification_lines = []
    if args.tool == "treesapp":
        ref_pkgs = ','.join(test_obj.ref_packages.keys())
        classification_table = os.path.join(args.output, "TreeSAPP_output", "final_outputs", "marker_contig_map.tsv")
        if not os.path.isfile(classification_table):
            classify_args = ["-i", args.input,
                             "-t", ref_pkgs,
                             "-n", str(args.num_threads),
                             "-m", "prot",
                             "--output", test_obj.data_dir,
                             "--stringency", args.stringency,
                             "--placement_summary", args.p_sum,
                             "--min_like_weight_ratio", str(args.min_lwr),
                             "--overwrite", "--delete"]
            if args.trim_align:
                classify_args.append("--trim_align")
            if args.svm:
                classify_args.append("--svm")
            assign(classify_args)
        classification_lines = file_parsers.read_classification_table(classification_table)
        assignments = assignments_to_treesaps(classification_lines, test_obj.ref_packages)
    else:
        # Since you are only able to analyze a single reference package at a time with GraftM, this is ran iteratively
        for gpkg in glob(args.refpkg_dir + "*gpkg"):
            pkg_name = str(os.path.basename(gpkg).split('.')[0])
            if not pkg_name:
                logging.error("Unable to parse marker name from gpkg '{}'\n".format(gpkg))
                sys.exit(5)
            if pkg_name not in test_obj.ref_packages:
                logging.warning("'{}' not in {} and will be skipped...\n".format(pkg_name, args.annot_map))
                continue
            output_dir = test_obj.data_dir + pkg_name + os.sep
            if not os.path.isdir(output_dir):
                os.makedirs(output_dir)
            classification_table = output_dir + test_fa_prefix + os.sep + test_fa_prefix + "_read_tax.tsv"
            if not os.path.isfile(classification_table):
                classify_call = ["graftM", "graft",
                                 "--forward", args.input,
                                 "--graftm_package", gpkg,
                                 "--input_sequence_type", "aminoacid",
                                 "--threads", str(args.num_threads),
                                 "--output_directory", output_dir,
                                 "--force"]
                if args.tool == "diamond":
                    classify_call += ["--assignment_method", "diamond"]
                    classify_call += ["--search_method", "diamond"]
                launch_write_command(classify_call, False)

            # TODO: Figure out how to convert GraftM classifications into JPlace objects
            assignments[pkg_name] = file_parsers.read_graftm_classifications(classification_table)

    if len(assignments) == 0:
        logging.error("No sequences were classified by " + args.tool + ".\n")
        sys.exit(3)

    logging.info("Reading headers in " + args.input + "... ")
    test_obj.all_queries = set([seq_name[1:] if seq_name[0] == '>' else seq_name for
                                seq_name in get_headers(args.input)])
    logging.info("done.\n")

    ##
    # Bin the test sequence names into their respective confusion categories (TP, TN, FP, FN)
    ##
    test_obj.bin_headers(assignments, query_to_code)
    test_obj.map_true_lineages()

    test_obj.bin_true_positives_by_taxdist()
    test_obj.validate_false_positives()
    test_obj.validate_false_negatives(refpkg_name_query_map)

    test_obj.summarise_reference_taxa(taxa_dist_output, classification_info_output, summary_rank)
    if args.tool == "treesapp" and classification_lines:
        logging.debug(test_obj.summarise_type_one_placements(classification_lines))
    logging.debug(test_obj.summarize_type_two_taxa(summary_rank))
    logging.debug(test_obj.true_positive_taxonomic_summary(summary_rank, True))

    ##
    # Report the MCC score across different taxonomic distances - should increase with greater allowed distance
    ##
    test_obj._MAX_TAX_DIST = 6
    logging.debug(test_obj.get_info(True))
    d = 0
    mcc_string = "Tax.dist\tMCC\tTrue.Pos\tTrue.Neg\tFalse.Pos\tFalse.Neg\n"
    while d < 8:
        test_obj._MAX_TAX_DIST = d
        tp, remainder = test_obj.get_true_positives_at_dist()
        num_tp = len(tp)
        num_fp = len(test_obj.get_false_positives()) + len(remainder)
        num_fn = len(test_obj.get_false_negatives())
        num_tn = test_obj.get_true_negatives()
        mcc = calculate_matthews_correlation_coefficient(num_tp, num_fp, num_fn, num_tn)
        mcc_string += "\t".join([str(x) for x in [d, mcc, num_tp, num_tn, num_fp, num_fn]]) + "\n"
        d += 1
    logging.info(mcc_string)
    with open(mcc_file, 'w') as mcc_handler:
        mcc_handler.write(mcc_string)
    return


if __name__ == '__main__':
    mcc_calculator(sys.argv[1:])
