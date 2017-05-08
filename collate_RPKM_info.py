#!/usr/bin/env python

import sys
import re
import os
import argparse


def get_arguments():
    parser = argparse.ArgumentParser(description="Used to collect RPKM values, gene names, contig name, and phylogeny into a single file from MLTreeMap outputs")
    parser.add_argument('-t', '--reftree', default='p', type=str,
                        help='reference tree (p = MLTreeMap reference tree [DEFAULT]; '
                             'g = GEBA reference tree; i = fungi tree; '
                             'other gene family in data/tree_data/cog_list.txt - e.g., y for mcrA]')
    parser.add_argument("--update_tree", action="store_true", default=False,
                        help="Flag indicating the reference tree specified by `--reftree` "
                             "is to be updated using the sequences found in MLTreeMap output")
    parser.add_argument('-r', "--rpkm_output",
                        help="RPKM output file (.csv)", required=True)
    parser.add_argument('-o', "--output_file",
                        help="Name of the rpkm output file", required=True)
    parser.add_argument("--output_dir_raxml", required=True)

    args = parser.parse_args()
    args.mltreemap = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + os.sep
    return args


class Autovivify(dict):
    """In cases of Autovivify objects, enable the referencing of variables (and sub-variables)
    without explicitly declaring those variables beforehand."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value


def create_cog_list(args):
    """
    Loads the MLTreeMap COG list file and check that the args.reftree exists
    :param args: The command-line and default arguments object
    :return: An autovivification of the COGs in cog_list.txt. This also includes their short-form name (termed
    denominator e.g. mcrA, mcrB, a, b, cdhA, etc.) and a list of output text precursors based on the analysis type.
    The denominator is equal to the command-line reference tree specifier argument (p, g, or i) if phylogenetic COGs
    """

    cog_list = Autovivify()
    text_of_analysis_type = Autovivify()
    cog_list_file = args.mltreemap + os.sep + 'data' + os.sep + 'tree_data' + os.sep + 'cog_list.txt'
    cog_input_list = open(cog_list_file, 'r')
    # TODO: alter function so the GEBA and fungi trees can be used instead of just the MLTreeMap reference
    if args.reftree not in ['i', 'p', 'g']:
        alignment_set = ''
    else:
        alignment_set = args.reftree
    kind_of_cog = ''

    # For each line in the COG list file...

    cog_list_lines = [x.strip() for x in cog_input_list.readlines()]
    # Close the COG list file
    cog_input_list.close()

    for cogInput in cog_list_lines:
        # Get the kind of COG if cogInput is a header line
        if re.match(r'\A#(.*)', cogInput):
            kind_of_cog = re.match(r'\A#(.*)', cogInput).group(1)
            continue

        if re.match(r'\A%(.*)', cogInput):
            continue

        # Add data to COG list based on the kind of COG it is
        if kind_of_cog == 'phylogenetic_cogs':
            cog_list[kind_of_cog][cogInput] = alignment_set
            cog_list['all_cogs'][cogInput] = alignment_set
            text_inset = ''
            if alignment_set == 'g':
                text_inset = ' based on the GEBA reference'
            if alignment_set == 'i':
                text_inset = ' focusing only on fungi'
            text_of_analysis_type[alignment_set] = 'Phylogenetic analysis' + text_inset + ':'
        elif kind_of_cog == 'phylogenetic_rRNA_cogs':
            cog, denominator, text = cogInput.split('\t')
            cog_list[kind_of_cog][cog] = denominator
            cog_list['all_cogs'][cog] = denominator
            text_of_analysis_type[denominator] = 'Phylogenetic analysis, ' + text + ':'
        elif kind_of_cog == 'functional_cogs':
            cog, denominator, text = cogInput.split('\t')
            cog_list[kind_of_cog][cog] = denominator
            cog_list['all_cogs'][cog] = denominator
            text_of_analysis_type[denominator] = 'Functional analysis, ' + text + ':'

    if args.reftree not in ['i', 'g', 'p']:
        if args.reftree not in cog_list['all_cogs'].values():
            sys.stderr.write("ERROR: " + args.reftree +
                             " not found in " + cog_list_file + "! Please use a valid reference tree ID!")
            sys.stderr.flush()
            sys.exit()
    elif args.reftree in ['i', 'g', 'p'] and args.update_tree:
        sys.stderr.write("ERROR: Unable to update all reference trees at the same time! "
                         "Please specify the reference tree to update with the '--reftree' argument and retry.\n")
        sys.stderr.flush()
        sys.exit()

    return cog_list, text_of_analysis_type


def collate_rpkm_info(rpkm_output_file, args, output, cog_list):
    contig_rpkm_map = dict()
    marker_map = dict()
    for marker in cog_list["all_cogs"]:
        denominator = cog_list["all_cogs"][marker]
        marker_map[denominator] = marker

    try:
        rpkm_values = open(rpkm_output_file, 'r')
    except:
        raise IOError("Unable to open " + rpkm_output_file + " for reading!")

    try:
        rpkm_out = open(output, 'w')
    except:
        raise IOError("Unable to open " + output + " for writing!")
    rpkm_out.write("Contig\tGene\tPhylogeny\tRPKM\n")

    for line in rpkm_values:
        contig, rpkm = line.strip().split(',')
        name, marker, start_end = contig.split('|')
        marker = re.sub("ref", '', marker)
        if name not in contig_rpkm_map:
            contig_rpkm_map[name] = dict()
        contig_rpkm_map[name][marker] = rpkm

    rpkm_values.close()

    final_raxml_outputs = os.listdir(args.output_dir_raxml)
    for raxml_contig_file in final_raxml_outputs:
        contig_name = '_'.join(re.sub("_RAxML_parsed.txt", '', raxml_contig_file).split('_')[1:])
        denominator = raxml_contig_file.split("_")[0]
        gene = re.sub("ref", '', marker_map[denominator])

        # By the end of the pipeline the gene was classified differently than the BLAST hit recommended
        if gene not in contig_rpkm_map[contig_name]:
            print contig_name, gene
            rpkm = 0
        else:
            rpkm = contig_rpkm_map[contig_name][gene]

        try:
            contig_placement = open(args.output_dir_raxml + raxml_contig_file, 'r')
        except:
            raise IOError("Unable to open " + args.output_dir_raxml + raxml_contig_file + " for reading!")
        line = contig_placement.readline()
        while not line.startswith("Placement"):
            line = contig_placement.readline().strip()

        placement = re.sub("^.*: Assignment of query to ", '', line)
        if contig_name not in contig_rpkm_map:
            raise AssertionError(contig_name + " not found in " + rpkm_output_file)
        rpkm_out.write('\t'.join([contig_name, gene, placement, str(rpkm)]) + "\n")

        contig_placement.close()

    return


def main():
    args = get_arguments()
    cog_list, text_of_analysis_type = create_cog_list(args)
    rpkm_file = args.rpkm_output
    output = args.output_file
    collate_rpkm_info(rpkm_file, args, output, cog_list)

main()
