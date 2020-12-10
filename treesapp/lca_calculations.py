__author__ = 'Connor Morgan-Lang'


import re
import logging

from pygtrie import StringTrie

from .utilities import median


def optimal_taxonomic_assignment(trie: StringTrie, query_taxon: str):
    while not trie.__contains__(query_taxon) and len(query_taxon.split('; ')) > 1:
        query_taxon = "; ".join(query_taxon.split('; ')[:-1])
    if not trie.__contains__(query_taxon):
        query_taxon = "r__Root"
    return query_taxon


def identify_excluded_clade(assignment_dict: dict, trie: StringTrie, marker: str) -> dict:
    """
    Using the taxonomic information from the sequence headers and the lineages of the reference sequence,
    this function determines the rank at which each sequence's clade is excluded.
    These data are returned and sorted in the form of a dictionary.

    :param assignment_dict:
    :param trie: A pygtrie.StringTrie object containing all reference sequence lineages
    :param marker: Name of the marker gene being tested

    :return: rank_assigned_dict; key is rank, values are dictionaries with assigned (reference) lineage as key and
      tuples of (optimal assignment, actual assignment) as values.
      E.g. {"Phylum": {"Proteobacteria": ("Proteobacteria", "Proteobacteria; Alphaproteobacteria")}}
    """
    rank_assigned_dict = dict()
    _RANK_DEPTH_MAP = {0: "root", 1: "domain",
                       2: "phylum", 3: "class", 4: "order",
                       5: "family", 6: "genus", 7: "species", 8: "strain"}

    if marker not in assignment_dict:
        logging.debug("No sequences assigned as " + marker + "\n")
        return rank_assigned_dict

    for depth in _RANK_DEPTH_MAP:
        rank_assigned_dict[_RANK_DEPTH_MAP[depth]] = list()
    log_stats = "Number of unique taxonomies that sequences were assigned to = " + \
                str(len(assignment_dict[marker].keys())) + "\n"

    for ref_lineage in assignment_dict[marker]:
        log_stats += "Assigned reference lineage: " + ref_lineage + "\n"
        for query_lineage in assignment_dict[marker][ref_lineage]:  # type: str
            if query_lineage.split('; ')[0] != "r__Root":
                query_lineage = "; ".join(["r__Root"] + query_lineage.split("; "))
            # if query_lineage == ref_lineage:
            #     logging.debug("\tQuery lineage: " + query_lineage + ", " +
            #                   "Optimal lineage: " + ref_lineage + "\n")
            # While the query_lineage contains clades which are not in the reference trie,
            # remove the taxonomic rank and traverse again. Complexity: O(ln(n))
            contained_taxonomy = optimal_taxonomic_assignment(trie, query_lineage)
            if len(contained_taxonomy.split("; ")) <= 7:
                rank_excluded = _RANK_DEPTH_MAP[len(contained_taxonomy.split("; "))]
                if contained_taxonomy != ref_lineage:
                    log_stats += "\tRank excluded: " + rank_excluded + "\n"
                    log_stats += "\t\tQuery lineage:   " + query_lineage + "\n"
                    log_stats += "\t\tOptimal lineage: " + contained_taxonomy + "\n"
                rank_assigned_dict[rank_excluded].append({ref_lineage: (contained_taxonomy, query_lineage)})
            else:
                logging.warning("Number of ranks in lineage '{}' is ridiculous.\n"
                                "This will not be used in clade exclusion calculations\n".format(contained_taxonomy))
    logging.debug(log_stats + "\n")
    return rank_assigned_dict


def megan_lca(lineage_list: list):
    """
    Using the lineages of all leaves to which this sequence was mapped (n >= 1),
    A lowest common ancestor is found at the point which these lineages converge.
    This emulates the LCA algorithm employed by the MEtaGenome ANalyzer (MEGAN).

    :param lineage_list: List of '; '-separated lineage strings
    :return:
    """
    # If there is only one child, return the joined string
    if len(lineage_list) == 1:
        return "; ".join(lineage_list[0])

    listed_lineages = [lineage.strip().split("; ") for lineage in lineage_list]
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
    if len(lca_lineage_strings) == 0:
        logging.debug("Empty LCA from lineages:\n\t" + "\n\t".join(lineage_list) + "\n")
        lca_lineage_strings.append("Unclassified")

    return "; ".join(lca_lineage_strings)


def weighted_taxonomic_distance(lineage_list, common_ancestor):
    """
    Input is a list >= 2, potentially either a leaf string or NCBI lineage

    :param lineage_list: Lineages of all leaves for this sequence
    :param common_ancestor: The common ancestor for the elements in lineage_list
    :return:
    """
    numerator = 0
    status = 0
    _MAX_DIST = 7
    for lineage in lineage_list:
        distance, status = compute_taxonomic_distance(lineage, common_ancestor)
        if status:
            logging.debug("Taxonomic lineages didn't converge between " + common_ancestor + " and " + lineage + ".\n")
        status += 1

        numerator += 2**distance

    wtd = round(float(numerator/(len(lineage_list) * 2**_MAX_DIST)), 5)
    return wtd, status


def compute_taxonomic_distance(ref_lineage: str, query_lineage: str):
    """
    Calculates the number of taxonomic ranks need to be climbed in the taxonomic hierarchy before a common ancestor
    is identified between the two taxa.
    If no common ancestor is reached, the distance is returned along with a status of 1 to indicate non-convergence.

    :param ref_lineage: A taxonomic lineage string, where each rank is separated by semi-colons (;)
    :param query_lineage: Another taxonomic lineage string, where each rank is separated by semi-colons (;)
    :return: Tuple of (distance, status)
    """
    distance = 0
    l1 = ref_lineage.split("; ")
    l2 = query_lineage.split("; ")
    # Compare the last elements of each list to see if the lineage is equal
    try:
        while l1 != l2:
            if len(l1) > len(l2):
                l1 = l1[:-1]
            elif len(l2) > len(l1):
                l2 = l2[:-1]
            else:
                # They are the same length, but disagree
                l2 = l2[:-1]
                l1 = l1[:-1]
            distance += 1
    except IndexError:
        return distance, 1

    return distance, 0


def determine_offset(classified: str, optimal: str) -> int:
    # Figure out which taxonomic lineage is longer
    offset = 0
    while classified != optimal and offset < 7:
        offset += 1
        classified_ranks = classified.split("; ")
        optimal_ranks = optimal.split("; ")
        if len(classified_ranks) > len(optimal_ranks):
            classified = "; ".join(classified_ranks[:-1])
        elif len(classified_ranks) < len(optimal_ranks):
            optimal = "; ".join(optimal_ranks[:-1])
        else:
            optimal = "; ".join(optimal_ranks[:-1])
            classified = "; ".join(classified_ranks[:-1])
    return offset


def clean_lineage_list(lineage_list):
    """
    Removes deeply unclassified sequences:
        If first rank == "unclassified sequences"
        If unclassified depth is < median unclassified depth
        cellular organisms; *; environmental samples

    :param lineage_list:
    :return: A list of lineages with mostly or entirely classified sequences, as long as they comprise the majority
    """
    if len(lineage_list) <= 1:
        return [lineage_list]

    classified_lineages = list()
    unclassified_depths = list()

    # Determine the median depth of "Unclassification" - maybe the entire clade is evolutionarily divergent
    for lineage_string in lineage_list:
        depth = 0
        lineage_split = lineage_string.split("; ")
        while depth < len(lineage_split):
            if re.search("unclassified", lineage_split[depth], re.IGNORECASE):
                break
            else:
                depth += 1
        unclassified_depths.append(depth)
    un_median = median(unclassified_depths)

    # Begin filtering
    for lineage_string in lineage_list:
        # Add lineage string if its all classified
        if not re.search("unclassified", lineage_string, re.IGNORECASE):
            classified_lineages.append(lineage_string)
            continue

        # If the depth of the unclassified rank is deep (close to the root), then skip
        lineage_split = lineage_string.split("; ")
        i = 0
        while i < len(lineage_split):
            if lineage_split[i] in ["unclassified sequences", "environmental samples"]:
                break
            elif re.search("unclassified", lineage_split[i], re.IGNORECASE) and i < un_median:
                break
            i += 1
        if i > un_median:
            classified_lineages.append(lineage_string)
            continue

    # Decide whether to return the cleaned list or not by determining whether the majority of lineages are included
    if float(len(classified_lineages)/len(lineage_list)) > 0.5:
        return classified_lineages
    else:
        return lineage_list
