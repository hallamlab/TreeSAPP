__author__ = 'Connor Morgan-Lang'


import sys
import re
import logging
from utilities import median


def disseminate_vote_weights(megan_lca, taxonomic_counts, lineages_list):
    lca_depth = 0
    max_depth = max([len(lineage) for lineage in lineages_list])
    nucleus = float(len(lineages_list)/taxonomic_counts[megan_lca])
    vote_weights = dict()
    taxonomic_tree = dict()
    subtree_taxonomic_counts = dict()
    vote_weights[megan_lca] = nucleus
    taxonomic_tree[megan_lca] = set()

    # Find the depth of the LCA
    while "; ".join(lineages_list[0][:lca_depth]) != megan_lca:
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


def lowest_common_taxonomy(children, megan_lca, taxonomic_counts, algorithm="LCA*"):
    """
    Input is a list >= 2, potentially either a leaf string or NCBI lineage

    :param children: Lineages of all leaves for this sequence
    :param megan_lca:
    :param taxonomic_counts:
    :param algorithm: A string indicating what lowest common ancestor algorithm should be used [ MEGAN | LCA* | LCAp ]
    :return: string - represent the consensus lineage for that node
    """
    lineages_considered = list()
    # Check that children have lineage information and discard those that don't have a known lineage
    # (e.g. unclassified sequences; metagenomes; ecological metagenomes)
    max_ranks = 0
    for child in children:
        # child = clean_lineage_string(child)
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
        lineage_string = megan_lca
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
        lineage_string = megan_lca
        vote_weights, t_tree, majority = disseminate_vote_weights(megan_lca, taxonomic_counts, lineages_considered)
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


def compute_taxonomic_distance(lineage_list, common_ancestor):
    """
    Input is a list >= 2, potentially either a leaf string or NCBI lineage

    :param lineage_list: Lineages of all leaves for this sequence
    :param common_ancestor: The common ancestor for the elements in lineage_list
    :return:
    """
    numerator = 0
    status = 0
    max_dist = 7
    lca_lineage = common_ancestor.split("; ")
    for lineage in lineage_list:
        # Compute the distance to common_ancestor
        lineage_path = lineage.split("; ")
        ref = lca_lineage
        query = lineage_path
        distance = 0
        # Compare the last elements of each list to see if the lineage is equal
        try:
            while ref[-1] != query[-1]:
                if len(ref) > len(query):
                    ref = ref[:-1]
                elif len(query) > len(ref):
                    query = query[:-1]
                else:
                    # They are the same length, but disagree
                    query = query[:-1]
                    ref = ref[:-1]
                distance += 1
        except IndexError:
            logging.debug("Taxonomic lineages didn't converge between " +
                          common_ancestor + " and " + lineage +
                          ". Taxonomic distance set to length of LCA lineage (" + str(len(lca_lineage)) + ").\n")
            distance = len(lca_lineage)
            status += 1

        numerator += 2**distance

    wtd = round(float(numerator/(len(lineage_list) * 2**max_dist)), 5)
    return wtd, status


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
