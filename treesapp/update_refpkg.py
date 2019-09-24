import os
import sys
import re
import logging
from pygtrie import StringTrie
from .external_command_interface import launch_write_command
from treesapp.utilities import clean_lineage_string, load_taxonomic_trie
from . import wrapper


def reformat_ref_seq_descriptions(original_header_map):
    reformatted_header_map = dict()
    for treesapp_id in original_header_map:
        try:
            organism_info, accession = original_header_map[treesapp_id].split(" | ")
            if organism_info and accession:
                reformatted_header_map[treesapp_id] = accession + " [" + organism_info + "]"
            else:
                reformatted_header_map[treesapp_id] = original_header_map[treesapp_id]
        except IndexError:
            reformatted_header_map[treesapp_id] = original_header_map[treesapp_id]
        # Remove the side-chevron character
        if reformatted_header_map[treesapp_id][0] == '>':
            reformatted_header_map[treesapp_id] = reformatted_header_map[treesapp_id][1:]
    return reformatted_header_map


def map_classified_seqs(ref_pkg_name: str, assignments: dict, unmapped_seqs: list) -> dict:
    """
    An algorithmically slow, O(n^2), process to assign sequences to their respective taxonomic lineages.
    This function is necessary as the sequence names in assignments and unmapped_seqs are not identical!

    :param ref_pkg_name: Name of the reference package (e.g. McrA). Necessary for matching regex of classified sequence
    :param assignments: Dictionary mapping queries assigned to taxonomic lineages, indexed by refpkg name
    :param unmapped_seqs: List of sequence names
    :return: Dictionary mapping classified query sequences to their respective taxonomic lineage
    """
    classified_seq_lineage_map = dict()
    for lineage in assignments[ref_pkg_name]:
        # Map the classified sequence to the header in FASTA
        x = 0
        while x < len(unmapped_seqs):
            seq_name = unmapped_seqs[x]
            original_name = re.sub(r"\|{0}\|\d+_\d+$".format(ref_pkg_name), '', seq_name)
            if original_name in assignments[ref_pkg_name][lineage]:
                classified_seq_lineage_map[seq_name] = lineage
                unmapped_seqs.pop(x)
            else:
                x += 1
    # Ensure all the classified sequences were mapped to lineages
    if unmapped_seqs:
        logging.error("Unable to map all classified sequences in to a lineage. These are missing:\n" +
                      "\n".join(unmapped_seqs) + "\n")
        sys.exit(5)
    return classified_seq_lineage_map


def validate_mixed_lineages(mixed_seq_lineage_map: dict) -> None:
    """
    Function to ensure all lineages begin at the same rank, typically either 'Root' or 'cellular organisms' if appropriate
    :param mixed_seq_lineage_map: A dictionary mapping sequence names (keys) to taxonomic lineages (values)
    :return: None
    """
    superfluous_prefixes = set()
    
    lineage_list = list(mixed_seq_lineage_map.values())
    taxa_trie = load_taxonomic_trie(lineage_list)  # type: StringTrie
    for taxon in taxa_trie:
        # Find all the prefixes that are inconsistent across lineages,
        # by seeing if the entire subtrees are present in the trie
        i = 0
        lineage_ranks = taxon.split('; ')
        while i < len(lineage_ranks) and taxa_trie.has_subtrie('; '.join(lineage_ranks[i+1:])):
            superfluous_prefixes.add(lineage_ranks[i])
            i += 1

    # Don't want to use `clean_lineage_string` here to maintain auxiliary ranks - just remove prefixes
    prefix_re = re.compile("|".join([prefix + "; " for prefix in superfluous_prefixes]))
    for seq_name in sorted(mixed_seq_lineage_map, key=lambda x: mixed_seq_lineage_map[x]):
        mixed_seq_lineage_map[seq_name] = prefix_re.sub('', mixed_seq_lineage_map[seq_name])

    return


def strip_assigment_pattern(seq_names: list, refpkg_name: str):
    """
    Strips the |RefPkg|start_stop pattern from the end of sequence names
    :param seq_names: A list of sequence names (headers) that were assigned using TreeSAPP
    :param refpkg_name: Name of the reference package that sequences were assigned to (e.g. McrA, nosZ)
    :return: Dictionary mapping the original headers to the new headers
    """
    return {seq_name: re.sub(r"\|{0}\|\d+_\d+$".format(refpkg_name), '', seq_name) for seq_name in seq_names}


def filter_by_lwr(classified_lines: list, min_lwr: float) -> set:
    high_lwr_placements = set()
    num_filtered = 0
    target_field = 8
    for classification in classified_lines:
        try:
            lwr = float(classification[target_field])
        except TypeError:
            logging.error("Unable to convert all classifications from column " + str(target_field) + " to a float.\n")
            sys.exit(17)

        assert 0.0 < lwr <= 1.0
        if lwr >= min_lwr:
            high_lwr_placements.add(classification[1])
        else:
            num_filtered += 1

    logging.debug(str(num_filtered) + " classified sequences did not meet minimum LWR of "
                  + str(min_lwr) + " for updating\n")

    return high_lwr_placements


def filter_by_lineage_depth(classified_lines: list, min_lineage_depth: int) -> set:
    resolved_placements = set()
    num_filtered = 0
    target_field = 5
    for classification in classified_lines:
        # Since we're dealing with classified sequences, the 'Root' prefix will also need to be removed
        lineage = clean_lineage_string(str(classification[target_field]), ["Root"])

        if lineage and len(lineage.split("; ")) >= min_lineage_depth:
            resolved_placements.add(classification[1])
        else:
            num_filtered += 1

    logging.debug(str(num_filtered) + " classified sequences did not meet lineage depth threshold of "
                  + str(min_lineage_depth) + " for updating\n")

    return resolved_placements


def intersect_incomparable_lists(superset, subset, name_map: dict) -> list:
    """

    :param superset: A list or set of strings
    :param subset: A list or set of strings
    :param name_map: A dictionary whose keys are in superset and values are in subset
    :return: A list of strings that are in both superset and subset
    """
    intersection = list()
    for seq_name in superset:  # type: str
        alt_name = name_map[seq_name]  # type: str
        if alt_name in subset:
            intersection.append(seq_name)
    return intersection
