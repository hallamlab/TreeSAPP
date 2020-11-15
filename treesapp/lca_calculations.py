__author__ = 'Connor Morgan-Lang'


import sys
import re
import logging

from pygtrie import StringTrie

from .utilities import median


def grab_graftm_taxa(tax_ids_file):
    taxonomic_tree = StringTrie(separator='; ')
    with open(tax_ids_file) as tax_ids:
        header = tax_ids.readline().strip()
        last_rank = int(header[-1])
        final_index = 6 - last_rank
        if not re.search("parent_id,rank,tax_name,root,rank_0,rank_1,rank_2,rank_3,rank_4,rank_5,rank_6", header):
            logging.error("Unable to handle format of " + tax_ids_file + "!")
            sys.exit(21)
        line = tax_ids.readline().strip()
        while line:
            fields = line.split(',')
            if final_index < 0:
                fields = line.split(',')[:final_index]
            try:
                _, _, _, _, _, k_, p_, c_, o_, f_, g_, s_, = fields
            except (IndexError, ValueError):
                logging.error("Unexpected format of line with %d fields in " % len(line.split(',')) +
                              tax_ids_file + ":\n" + line)
                sys.exit(21)
            ranks = ["Root", k_, p_, c_, o_, f_, g_, s_]
            lineage_list = []
            # In case there are missing ranks... which is likely
            for rank in ranks:
                if rank:
                    # GraftM seems to append an 'e1' to taxa that are duplicated in the taxonomic lineage.
                    # For example: Bacteria; Aquificae; Aquificaee1; Aquificales
                    lineage_list.append(re.sub(r'_graftm_\d+$', '', rank))
                    # lineage_list.append(rank)
            lineage = re.sub('_', ' ', '; '.join(lineage_list))
            i = 0
            ranks = len(lineage)
            while i < len(lineage):
                taxonomic_tree["; ".join(lineage.split("; ")[:ranks - i])] = True
                i += 1

            line = tax_ids.readline().strip()
    return taxonomic_tree


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


def disseminate_vote_weights(base_lca, taxonomic_counts, lineages_list):
    lca_depth = 0
    max_depth = max([len(lineage) for lineage in lineages_list])
    nucleus = float(len(lineages_list)/taxonomic_counts[base_lca])
    vote_weights = dict()
    taxonomic_tree = dict()
    subtree_taxonomic_counts = dict()
    vote_weights[base_lca] = nucleus
    taxonomic_tree[base_lca] = set()

    # Find the depth of the LCA
    while "; ".join(lineages_list[0][:lca_depth]) != base_lca:
        lca_depth += 1

    # Load the taxonomic subtree into a dictionary
    while lca_depth < max_depth:
        for lineage in lineages_list:
            taxonomic_lineage = "; ".join(lineage[0:lca_depth+1])
            parent = "; ".join(lineage[0:lca_depth])
            if parent == taxonomic_lineage:
                continue
            if taxonomic_lineage not in taxonomic_tree:
                taxonomic_tree[taxonomic_lineage] = set()
            if taxonomic_lineage not in subtree_taxonomic_counts:
                subtree_taxonomic_counts[taxonomic_lineage] = 0
            subtree_taxonomic_counts[taxonomic_lineage] += 1
            taxonomic_tree[parent].add(taxonomic_lineage)
        lca_depth += 1

    # Calculate the vote weights for each taxonomic lineage based on it's parent's and polyphyly
    for node in sorted(taxonomic_tree, key=lambda x: x.count(';')):
        children = taxonomic_tree[node]
        n_taxa_split = len(children)
        for child in children:
            taxa_portion = float(subtree_taxonomic_counts[child]/taxonomic_counts[child])
            vote_weights[child] = vote_weights[node] * float(taxa_portion/n_taxa_split)

    return vote_weights, taxonomic_tree, nucleus/2


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


def lowest_common_taxonomy(children, base_lca, taxonomic_counts, algorithm="LCA*"):
    """
    Input is a list >= 2, potentially either a leaf string or NCBI lineage

    :param children: Lineages of all leaves for this sequence
    :param base_lca:
    :param taxonomic_counts:
    :param algorithm: A string indicating what lowest common ancestor algorithm should be used [ MEGAN | LCA* | LCAp ]
    :return: string - represent the consensus lineage for that node
    """
    lineages_considered = list()
    # Check that children have lineage information and discard those that don't have a known lineage
    # (e.g. unclassified sequences; metagenomes; ecological metagenomes)
    max_ranks = 0
    for child in children:
        ranks = child.split("; ")
        num_ranks = len(ranks)
        if num_ranks > 3:
            lineages_considered.append(ranks)
        if num_ranks > max_ranks:
            max_ranks = num_ranks
    # print("All children:", children)
    # print("Lineages considered:", lineages_considered)
    if len(lineages_considered) == 0:
        lineages_considered = [child.split("; ") for child in children]

    if algorithm == "MEGAN":
        # Already calculated by tree_sap.megan_lca()
        lineage_string = base_lca
    elif algorithm == "LCA*":
        # approximate LCA* (no entropy calculations):
        hits = dict()
        consensus = list()
        i = 0  # The accumulator to guarantee the lineages are parsed from Kingdom -> Strain
        while i < max_ranks:
            hits.clear()
            lineages_used = 0
            elected = False
            for ranks in lineages_considered:
                # If the ranks of this hit has not been exhausted (i.e. if it has a depth of 4 and max_ranks >= 5)
                if len(ranks) > i:
                    taxonomy = ranks[i]
                    if taxonomy not in hits.keys():
                        hits[taxonomy] = 0
                    hits[taxonomy] += 1
                    lineages_used += 1
            for taxonomy in hits.keys():
                if hits[taxonomy] > float(len(lineages_considered)/2):
                    consensus.append(str(taxonomy))
                    elected = True

            # If there is no longer a consensus, break the loop
            if not elected:
                i = max_ranks
            i += 1
        lineage_string = "; ".join(consensus)
    elif algorithm == "LCAp":
        lineage_string = base_lca
        vote_weights, t_tree, majority = disseminate_vote_weights(base_lca, taxonomic_counts, lineages_considered)
        for parent in sorted(t_tree, key=lambda x: x.count(';')):
            if len(t_tree[parent]) == 0:
                break
            children = t_tree[parent]
            for child in children:
                if vote_weights[child] > majority:
                    lineage_string = child
                else:
                    pass
    else:
        logging.error("Common ancestor algorithm '" + algorithm + "' is not known.\n")
        sys.exit(33)

    return lineage_string


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
