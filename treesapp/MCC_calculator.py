#!/usr/bin/env python3

__author__ = 'Connor Morgan-Lang'

import argparse
import os
import sys
import re
import logging
import shutil
from glob import glob
from numpy import sqrt
from treesapp.commands import assign
from treesapp.fasta import get_headers
from treesapp.external_command_interface import launch_write_command
from treesapp import file_parsers
from treesapp.classy import prep_logging, ReferencePackage
from treesapp.entrez_utils import EntrezRecord, fetch_lineages_from_taxids
from treesapp.lca_calculations import compute_taxonomic_distance, all_possible_assignments, \
    optimal_taxonomic_assignment, grab_graftm_taxa
from treesapp.utilities import fish_refpkg_from_build_params, clean_lineage_string


class ClassifiedSequence:
    def __init__(self, header):
        self.name = header
        self.ref = ""
        self.ncbi_tax = ""
        self.assigned_lineage = ""
        self.optimal_lineage = ""
        self.true_lineage = ""
        self.optimal_lineage = ""
        self.tax_dist = 0


class ConfusionTest:
    def __init__(self, gene_list):
        self._MAX_TAX_DIST = -1
        self.header_regex = None
        self.ref_packages = {key: ReferencePackage() for key in gene_list}
        self.fn = {key: [] for key in gene_list}
        self.fp = {key: [] for key in gene_list}
        self.tp = {key: [] for key in gene_list}  # This will be a list of ClassifiedSequence instances
        self.tax_lineage_map = dict()
        self.dist_wise_tp = dict()
        self.num_total_queries = 0
        self.rank_depth_map = {"Kingdom": 1, "Phylum": 2, "Class": 3, "Order": 4, "Family": 5, "Genus": 6, "Species": 7}
        self.classification_table = ""
        self.data_dir = ""

    def get_info(self, verbose=False):
        info_string = "\nReference packages being tested:\n"
        info_string += "\t" + ", ".join(list(self.ref_packages.keys())) + "\n"

        self.check_dist()
        info_string += "Stats based on taxonomic distance < " + str(self._MAX_TAX_DIST) + "\n"

        if self.num_total_queries > 0:
            info_string += str(self.num_total_queries) + " query sequences being used for testing.\n"

        fp_ogs = set()
        for marker in self.fp:
            if len(self.fp[marker]) > 0:
                info_string += "False positive OGs classified as " + marker + ":\n"
                for query_name in self.fp[marker]:
                    fp_ogs.add(self.header_regex.match(query_name).group(1))
                info_string += "\t" + ', '.join(fp_ogs) + "\n"
            fp_ogs.clear()

        if verbose:
            for refpkg in sorted(self.ref_packages):
                info_string += self.marker_classification_summary(refpkg)

        return info_string

    def map_data(self, output_dir, tool):
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

    def populate_tax_lineage_map(self, classified_seq_list, records_list, group=1):
        for seq_name in classified_seq_list:  # type: str
            try:
                tax_id = self.header_regex.search(seq_name).group(group)
            except (KeyError, TypeError):
                logging.error("Header '" + str(seq_name) + "' in test FASTA doesn't match the supported format.\n")
                sys.exit(5)

            # Only make Entrez records for new NCBI taxonomy IDs
            if tax_id not in self.tax_lineage_map:
                e_record = EntrezRecord(seq_name, "")
                e_record.ncbi_tax = tax_id
                e_record.bitflag = 3
                records_list.append(e_record)
                self.tax_lineage_map[e_record.ncbi_tax] = "Root; "
        return records_list

    def retrieve_lineages(self, group=1):
        if not self.header_regex:
            logging.error("Unable to parse taxonomic identifiers from header without a regular expression.\n")
            sys.exit(19)

        # Gather the unique taxonomy IDs and store in EntrezRecord instances
        entrez_records = list()
        for marker in self.tp:
            entrez_records = self.populate_tax_lineage_map([tp_inst.name for tp_inst in self.tp[marker]],
                                                           entrez_records, group)
        for marker in self.fn:
            entrez_records = self.populate_tax_lineage_map(self.fn[marker], entrez_records, group)

        # Query the Entrez database for these unique taxonomy IDs
        fetch_lineages_from_taxids(entrez_records)

        for e_record in entrez_records:  # type: EntrezRecord
            if not e_record.lineage:
                logging.warning("Lineage information unavailable for taxonomy ID '" + e_record.ncbi_tax + "'\n")
                self.tax_lineage_map[e_record.ncbi_tax] = ''
            else:
                self.tax_lineage_map[e_record.ncbi_tax] += clean_lineage_string(e_record.lineage)
        return

    def map_lineages(self):
        for marker in self.tp:
            for tp_inst in self.tp[marker]:  # type: ClassifiedSequence
                try:
                    tp_inst.true_lineage = self.tax_lineage_map[tp_inst.ncbi_tax]
                except KeyError:
                    continue
        return

    def summarise_reference_taxa(self, taxa_file, classification_file, rank="Phylum"):
        lineage_list = []
        info_string = "RefPkg\tName\tTaxDist\tClassified\tTrueLineage\tAssignedLineage\tOptimalLineage\n"
        for marker in self.tp:
            for tp_inst in self.tp[marker]:  # type: ClassifiedSequence
                lineage_list.append(tp_inst.true_lineage)
                info_string += "\t".join([marker, tp_inst.name, str(tp_inst.tax_dist), "True",
                                          tp_inst.true_lineage, tp_inst.assigned_lineage,
                                          tp_inst.optimal_lineage]) + "\n"
        for marker in self.fn:
            for seq_name in self.fn[marker]:
                ref_g, tax_id = self.header_regex.match(seq_name).groups()
                try:
                    true_lineage = self.tax_lineage_map[tax_id]
                except KeyError:
                    continue
                lineage_list.append(true_lineage)
                info_string += "\t".join([marker, seq_name, "NA", "False", true_lineage, "NA", "NA"]) + "\n"

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
            self.dist_wise_tp[marker] = dict()
            for tp_inst in self.tp[marker]:  # type: ClassifiedSequence
                # Find the optimal taxonomic assignment
                optimal_taxon = optimal_taxonomic_assignment(self.ref_packages[marker].taxa_trie, self.tax_lineage_map[tp_inst.ncbi_tax])
                if not optimal_taxon:
                    logging.debug("Optimal taxonomic assignment '" + tp_inst.true_lineage + "' for " +
                                  tp_inst.name + " not found in reference hierarchy.\n")
                    continue
                tp_inst.optimal_lineage = optimal_taxon
                tp_inst.tax_dist, status = compute_taxonomic_distance(tp_inst.assigned_lineage, optimal_taxon)
                if status > 0:
                    logging.debug("Lineages didn't converge between:\n'" +
                                  tp_inst.assigned_lineage + "' and '" + optimal_taxon +
                                  "' (taxid: " + tp_inst.ncbi_tax + ")\n")
                try:
                    self.dist_wise_tp[marker][tp_inst.tax_dist].append(tp_inst.name)
                except KeyError:
                    self.dist_wise_tp[marker][tp_inst.tax_dist] = [tp_inst.name]
        return

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
            unique_fn.update(self.fn[marker])
            unique_fp.update(self.fp[marker])
            unique_tp.update(self.tp[marker])
        return self.num_total_queries - sum([len(unique_fn), len(unique_fp), len(unique_tp)])

    def get_false_positives(self, refpkg_name=None):
        if refpkg_name:
            return self.fp[refpkg_name]
        else:
            unique_fp = set()
            for marker in self.ref_packages:
                unique_fp.update(self.fp[marker])
            return unique_fp

    def get_false_negatives(self, refpkg_name=None):
        if refpkg_name:
            return self.fn[refpkg_name]
        else:
            unique_fn = set()
            for marker in self.ref_packages:
                # Remove sequences that were classified by a homologous marker
                unique_fn.update(self.fn[marker])
            return unique_fn

    def bin_headers(self, test_seq_names, assignments, annot_map, marker_build_dict):
        """
        Function for sorting/binning the classified sequences at T/F positives/negatives based on the
        :param test_seq_names: List of all headers in the input FASTA file
        :param assignments: Dictionary mapping taxonomic lineages to a list of headers that were classified as lineage
        :param annot_map: Dictionary mapping reference package (gene) name keys to database names values
        :param marker_build_dict:
        :return: None
        """
        # False positives: those that do not belong to the annotation matching a reference package name
        # False negatives: those whose annotations match a reference package name and were not classified
        # True positives: those with a matching annotation and reference package name and were classified
        # True negatives: those that were not classified and should not have been
        mapping_dict = dict()
        positive_queries = dict()

        for refpkg in annot_map:
            marker = marker_build_dict[refpkg].cog
            orthos = annot_map[refpkg]  # List of all orthologous genes corresponding to a reference package
            for gene in orthos:
                if gene not in mapping_dict:
                    mapping_dict[gene] = []
                mapping_dict[gene].append(marker)
        og_names_mapped = list(mapping_dict.keys())

        logging.info("Labelling true test sequences... ")
        for header in test_seq_names:
            try:
                ref_g, tax_id = self.header_regex.match(header).groups()
            except (TypeError, KeyError):
                logging.error("Header '" + header + "' in test FASTA file does not match the supported format.\n")
                sys.exit(5)
            if ref_g in mapping_dict:
                markers = mapping_dict[ref_g]
                ##
                # This leads to double-counting and is therefore deduplicated later
                ##
                for marker in markers:
                    if not marker:
                        logging.error("Bad marker name in " + str(mapping_dict.keys()) + "\n")
                        sys.exit(5)
                    if marker not in positive_queries:
                        positive_queries[marker] = []
                        if ref_g in og_names_mapped:
                            i = 0
                            while i < len(og_names_mapped):
                                if og_names_mapped[i] == ref_g:
                                    og_names_mapped.pop(i)
                                    i = len(og_names_mapped)
                                i += 1
                    positive_queries[marker].append(header)
        logging.info("done.\n")

        # Ensure all reference genes in mapping_dict have been used
        if len(og_names_mapped) > 0:
            logging.warning("Some orthologous groups in the annotation mapping file were not found in the FASTA file." +
                            " Perhaps a mistake was made when making this file? The following OGs will be skipped:\n" +
                            '\n'.join([str(og) + ": " + str(mapping_dict[og]) for og in og_names_mapped]) + "\n")

        logging.info("Assigning test sequences to the four class conditions... ")
        for marker in assignments:
            try:
                positives = set(positive_queries[marker])
            except KeyError:
                logging.error("Unable to find '" + marker + "' in the set of positive queries:\n" +
                              ", ".join([str(n) for n in positive_queries.keys()]) + "\n")
                sys.exit(5)
            true_positives = set()
            refpkg = fish_refpkg_from_build_params(marker, marker_build_dict).denominator
            for tax_lin in assignments[marker]:
                classified_seqs = assignments[marker][tax_lin]
                for seq_name in classified_seqs:
                    try:
                        ref_g, tax_id = self.header_regex.match(seq_name).groups()
                    except (TypeError, KeyError, AttributeError):
                        logging.error(
                            "Classified sequence name '" + seq_name + "' does not match the supported format.\n")
                        sys.exit(5)
                    if ref_g in mapping_dict and marker in mapping_dict[ref_g]:
                        # Populate the relevant information for the classified sequence
                        tp_inst = ClassifiedSequence(seq_name)
                        tp_inst.ncbi_tax = tax_id
                        tp_inst.ref = marker
                        tp_inst.assigned_lineage = tax_lin

                        # Add the True Positive to the relevant collections
                        self.tp[refpkg].append(tp_inst)
                        true_positives.add(seq_name)
                    else:
                        self.fp[refpkg].append(seq_name)

            # Identify the False Negatives using set difference - those that were not classified but should have been
            self.fn[refpkg] = list(positives.difference(true_positives))
        logging.info("done.\n")
        return

    def og_names(self, seq_names):
        og_names = []
        for query in seq_names:
            try:
                seq_name = query.name
            except AttributeError:
                seq_name = query
            og = self.header_regex.match(seq_name).group(1)
            original_name = re.sub(og, '', seq_name)
            if original_name[0] == '>':
                original_name = original_name[1:]
            if original_name[0] == '_':
                original_name = original_name[1:]
            og_names.append(original_name)
        return og_names

    def validate_false_positives(self):
        # Get all the original Orthologous Group (OG) headers for sequences classified as TP or FN
        tp_names = []
        for marker in list(self.ref_packages.keys()):
            tp_names += self.og_names(self.tp[marker])
            tp_names += self.og_names(self.fn[marker])

        for marker in self.fp:
            validated_fp = []
            for seq_name in self.fp[marker]:
                og_name = self.header_regex.match(seq_name).group(1)
                original_name = re.sub(og_name + '_', '', seq_name)
                if original_name[0] == '>':
                    original_name = original_name[1:]
                if original_name not in tp_names:
                    validated_fp.append(seq_name)
            self.fp[marker] = validated_fp
        return

    def validate_false_negatives(self, pkg_name_dict):
        # Invert the dictionary
        refpkg_og_map = dict()
        for refpkg_name in pkg_name_dict:
            ogs = pkg_name_dict[refpkg_name]
            for og in ogs:
                try:
                    refpkg_og_map[og].append(refpkg_name)
                except KeyError:
                    refpkg_og_map[og] = [refpkg_name]

        for marker in self.fn:
            homologous_tps = set()
            for og in pkg_name_dict[marker]:
                for homolgous_marker in refpkg_og_map[og]:
                    if homolgous_marker != marker:
                        homologous_tps.update(set(cs.name for cs in self.tp[homolgous_marker]))
            # Remove all sequences from this marker's false negatives that are found in a homologous TP set
            self.fn[marker] = set(self.fn[marker]).difference(homologous_tps)
        return

    def summarize_type_two_taxa(self, rank="Phylum"):
        lineage_list = []
        summary_string = ""
        for marker in self.fn:
            if len(self.fn[marker]) == 0:
                continue
            summary_string += "Number of false negatives for '" + marker + "':\n"
            for seq_name in self.fn[marker]:
                ref_g, tax_id = self.header_regex.match(seq_name).groups()
                try:
                    lineage_list.append(self.tax_lineage_map[tax_id])
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


def summarize_taxonomy(taxa_list, rank, rank_depth_map=None):
    """
    Given a list of taxonomic lineages and a taxonomic rank for which to summarise at
    it will count the number of instances for each taxon at the desired rank.
    E.g. {''}
    :return: A dictionary mapping a taxon to the number of representatives seen
    """
    if not rank_depth_map:
        rank_depth_map = {"Kingdom": 1, "Phylum": 2, "Class": 3, "Order": 4, "Family": 5, "Genus": 6, "Species": 7}

    try:
        depth = rank_depth_map[rank] + 1
    except KeyError:
        logging.error("Rank '" + rank + "' not present in rank-depth map:\n" + str(rank_depth_map) + "\n")
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


def get_arguments():
    parser = argparse.ArgumentParser(add_help=False,
                                     description="A wrapper script for calculating Matthews correlation coefficient for"
                                                 "TreeSAPP or GraftM. Currently only supports testing proteins.")
    required_args = parser.add_argument_group("Required arguments")
    required_args.add_argument("--fasta_input", required=True, dest="input",
                               help="FASTA-formatted file used for testing the classifiers")
    required_args.add_argument("--annot_map", required=True,
                               help="Path to a tabular file mapping markers being tested to their database annotations."
                                    " First column is the ")

    optopt = parser.add_argument_group("Optional options")
    optopt.add_argument("--tool", default="treesapp", required=False,
                        choices=["treesapp", "graftm", "diamond"],
                        help="Classify using one of the tools: treesapp [DEFAULT], graftm, or diamond.")
    optopt.add_argument("--gpkg_dir", required=False, default=None,
                        help="Path to a directory containing all GraftM reference packages to test."
                             " Files must follow 'name.gpkg' scheme and 'name' is in the first column of --annot_map")
    optopt.add_argument("--output", required=False, default="./MCC_output/",
                        help="Path to a directory for writing output files")
    optopt.add_argument("-p", "--pkg_path", required=False, default=None,
                        help="The path to the TreeSAPP-formatted reference package(s) [ DEFAULT = TreeSAPP/data/ ].")

    miscellaneous_opts = parser.add_argument_group("Miscellaneous options")
    miscellaneous_opts.add_argument("-n", "--num_threads", default=4, required=False,
                                    help="The number of threads to use when running either TreeSAPP or GraftM")
    miscellaneous_opts.add_argument('--overwrite', action='store_true', default=False,
                                    help='overwrites previously processed output folders')
    miscellaneous_opts.add_argument('-v', '--verbose', action='store_true', default=False,
                                    help='Prints a more verbose runtime log')
    miscellaneous_opts.add_argument("-h", "--help",
                                    action="help",
                                    help="Show this help message and exit")

    args = parser.parse_args()

    if args.output[-1] != os.sep:
        args.output += os.sep
    return args


def validate_command(args, sys_args):
    logging.debug("Command used:\n" + ' '.join(sys_args) + "\n")

    if args.overwrite:
        if os.path.exists(args.output):
            shutil.rmtree(args.output)
        os.mkdir(args.output)

    args.treesapp = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + os.sep
    if args.tool == "treesapp" and not args.pkg_path:
        args.pkg_path = args.treesapp + "data" + os.sep

    if sys.version_info < (2, 9):
        logging.error("Python version '" + sys.version_info + "' not supported.\n")
        sys.exit(3)

    if args.gpkg_dir and args.gpkg_dir[-1] != os.sep:
        args.gpkg_dir += os.sep
    if args.tool in ["diamond", "graftm"]:
        if not args.gpkg_dir:
            logging.error(args.tool + " specified but a GraftM reference package directory was not provided.\n")
            sys.exit(17)
        elif not os.path.isdir(args.gpkg_dir):
            logging.error(args.gpkg_dir + " GraftM reference package directory does not exist!\n")
            sys.exit(17)
        elif len(glob(args.gpkg_dir + "*gpkg")) == 0:
            logging.error("Not GraftM reference packages found in " + args.gpkg_dir + ".\n")
            sys.exit(17)
    elif args.tool == "treesapp" and args.gpkg_dir:
        logging.warning("--gpkg specified but tool selected is 'treesapp'... it will be ignored.\n")

    return


def read_annotation_mapping_file(annot_map_file):
    annot_map = dict()
    try:
        annot_map_handler = open(annot_map_file)
    except IOError:
        logging.error("Unable to open " + annot_map_file + " for reading!\n")
        sys.exit(3)

    # Assuming the first column is the reference package name and the second is the database annotation name
    n = 0
    for line in annot_map_handler:
        n += 1
        if line[0] == '#':
            continue
        elif not line:
            continue
        else:
            fields = line.split("\t")
            try:
                annot_map[fields[0]] = fields[1].split(',')
            except KeyError:
                logging.error("Insufficient number of fields on line " + str(n) + " in " + annot_map_file + "!\n" +
                              "File must have the reference package name and the database name in"
                              " the first two columns, respectively. Any number of columns can follow.\n")
                sys.exit(9)

    annot_map_handler.close()
    return annot_map


def calculate_matthews_correlation_coefficient(tp: int, fp: int, fn: int, tn: int):
    numerator = float((tp * tn) - (fp * fn))
    denominator = float((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    if numerator == 0 or denominator == 0:
        return 0
    else:
        return round(numerator/sqrt(denominator), 3)


def mcc_calculator():
    args = get_arguments()
    log_name = args.output + os.sep + "MCC_log.txt"
    mcc_file = args.output + os.sep + "MCC_table.tsv"
    summary_rank = "Phylum"
    taxa_dist_output = args.output + '.'.join(os.path.basename(args.input).split('.')[:-1]) + '_' + summary_rank + "_dist.tsv"
    classification_info_output = args.output + '.'.join(os.path.basename(args.input).split('.')[:-1]) + "_classifications.tsv"
    prep_logging(log_name, args.verbose)
    logging.info("\n##\t\t\tBeginning Matthews Correlation Coefficient analysis\t\t\t##\n")
    validate_command(args, sys.argv)

    ##
    # Read the file mapping reference package name to the database annotations
    ##
    pkg_name_dict = read_annotation_mapping_file(args.annot_map)
    marker_build_dict = file_parsers.parse_ref_build_params(args.treesapp, [])
    test_obj = ConfusionTest(pkg_name_dict.keys())
    test_obj.map_data(output_dir=args.output, tool=args.tool)
    if args.overwrite and os.path.exists(test_obj.data_dir):
            shutil.rmtree(test_obj.data_dir)

    ##
    # Load the taxonomic trie for each reference package
    ##
    if args.tool == "treesapp":
        for pkg_name in test_obj.ref_packages:
            refpkg = test_obj.ref_packages[pkg_name]
            marker = marker_build_dict[pkg_name].cog
            refpkg.prefix = marker
            refpkg.gather_package_files(args.pkg_path)
            test_obj.ref_packages[pkg_name].taxa_trie = all_possible_assignments(test_obj.ref_packages[pkg_name].lineage_ids)
    else:
        for gpkg in glob(args.gpkg_dir + "*gpkg"):
            marker = str(os.path.basename(gpkg).split('.')[0])
            pkg_name = fish_refpkg_from_build_params(marker, marker_build_dict).denominator
            if pkg_name in pkg_name_dict:
                try:
                    tax_ids_file = glob(os.sep.join([gpkg, marker + ".gpkg.refpkg", marker + "*_taxonomy.csv"])).pop()
                    test_obj.ref_packages[pkg_name].taxa_trie = grab_graftm_taxa(tax_ids_file)
                except IndexError:
                    logging.warning("No GraftM taxonomy file found for " + marker + ". Is this refpkg incomplete?\n")

    ##
    # Run the specified taxonomic analysis tool and collect the classifications
    ##
    assignments = {}
    test_fa_prefix = '.'.join(os.path.basename(args.input).split('.')[:-1])
    if args.tool == "treesapp":
        ref_pkgs = ','.join(pkg_name_dict.keys())
        classification_table = os.sep.join([args.output, "TreeSAPP_output", "final_outputs", "marker_contig_map.tsv"])
        if not os.path.isfile(classification_table):
            classify_args = ["-i", args.input,
                             "-t", ref_pkgs, "-n", str(args.num_threads),
                             "-m", "prot",
                             "--output", test_obj.data_dir,
                             "--trim_align", "--overwrite"]
            assign(classify_args)
        classification_lines = file_parsers.read_marker_classification_table(classification_table)
        assignments = file_parsers.parse_assignments(classification_lines)
    else:
        # Since you are only able to analyze a single reference package at a time with GraftM, this is ran iteratively
        for gpkg in glob(args.gpkg_dir + "*gpkg"):
            marker = str(os.path.basename(gpkg).split('.')[0])
            if not marker:
                logging.error("Unable to parse marker name from gpkg '" + gpkg + "'\n")
                sys.exit(5)
            pkg_name = fish_refpkg_from_build_params(marker, marker_build_dict).denominator
            if pkg_name not in pkg_name_dict:
                logging.warning("'" + pkg_name + "' not in " + args.annot_map + " and will be skipped...\n")
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

            assignments[marker] = file_parsers.read_graftm_classifications(classification_table)

    if len(assignments) == 0:
        logging.error("No sequences were classified by " + args.tool + ".\n")
        sys.exit(3)

    logging.info("Reading headers in " + args.input + "... ")
    test_seq_names = [seq_name[1:] if seq_name[0] == '>' else seq_name for seq_name in get_headers(args.input)]
    logging.info("done.\n")
    test_obj.num_total_queries = len(test_seq_names)
    eggnog_re = re.compile(r"^>?(COG[A-Z0-9]+|ENOG[A-Z0-9]+)_(\d+)\..*")
    test_obj.header_regex = eggnog_re

    ##
    # Bin the test sequence names into their respective confusion categories (TP, TN, FP, FN)
    ##
    test_obj.bin_headers(test_seq_names, assignments, pkg_name_dict, marker_build_dict)
    test_seq_names.clear()

    ##
    # Parse the taxonomic IDs from EggNOG headers and map taxonomic lineage information to classified sequences
    ##
    _TAXID_GROUP = 2
    test_obj.retrieve_lineages(_TAXID_GROUP)
    test_obj.map_lineages()

    test_obj.bin_true_positives_by_taxdist()
    test_obj.validate_false_positives()
    test_obj.validate_false_negatives(pkg_name_dict)

    test_obj.summarise_reference_taxa(taxa_dist_output, classification_info_output, summary_rank)
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
    mcc_calculator()
