import sys
import re
import logging

from pygtrie import StringTrie

from treesapp.utilities import load_taxonomic_trie
from treesapp.classy import Cluster, Updater
from treesapp.entrez_utils import EntrezRecord
from treesapp.fasta import FASTA
from treesapp.taxonomic_hierarchy import TaxonomicHierarchy
from treesapp import refpkg as ts_ref_pkg


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
    Function to ensure all lineages begin at the same rank, typically either 'Root' or 'cellular organisms'

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

    prefix_re = re.compile("|".join([prefix + "; " for prefix in superfluous_prefixes]))
    for seq_name in sorted(mixed_seq_lineage_map, key=lambda x: mixed_seq_lineage_map[x]):
        mixed_seq_lineage_map[seq_name] = prefix_re.sub('', mixed_seq_lineage_map[seq_name])

    return


def strip_assigment_pattern(seq_names: list, refpkg_name: str) -> dict:
    """
    Strips the |RefPkg|start_stop pattern from the end of sequence names
    :param seq_names: A list of sequence names (headers) that were assigned using TreeSAPP
    :param refpkg_name: Name of the reference package that sequences were assigned to (e.g. McrA, nosZ)
    :return: Dictionary mapping the original headers to the new headers
    """
    return {seq_name: re.sub(r"\|{0}\|\d+_\d+$".format(refpkg_name), '', seq_name) for seq_name in seq_names}


def filter_by_lwr(classified_lines: list, min_lwr: float) -> set:
    """

    :param classified_lines:
    :param min_lwr: A float representing the minimum acceptable Likelihood Weight Ratio as calculated by EPA-NG
    :return: A set of sequence names with corresponding placement likelihoods greater than min_lwr
    """
    high_lwr_placements = set()
    num_filtered = 0  # Number of placements with LWR less than min_lwr
    target_field = 9  # The field index (from the classification table) storing LWR
    for classification in classified_lines:
        try:
            lwr = float(classification[target_field])
        except TypeError:
            logging.error("Unable to convert all classifications from column {} to a float.\n".format(target_field))
            sys.exit(17)

        assert 0.0 < lwr <= 1.0
        if lwr >= min_lwr:
            high_lwr_placements.add(classification[1])
        else:
            num_filtered += 1

    logging.debug("{} classified sequences did not meet minimum LWR of {} for updating\n".format(num_filtered, min_lwr))

    return high_lwr_placements


def decide_length_filter(refpkg: ts_ref_pkg.ReferencePackage, proposed_min_length=30, min_hmm_proportion=0.66) -> int:
    # Ensure the HMM fraction provided is a proportion
    if not 0 < min_hmm_proportion < 1:
        logging.warning("Minimum HMM fraction provided ({}) isn't a proportion. Converting to {}.\n"
                        "".format(min_hmm_proportion, min_hmm_proportion/100))
        min_hmm_proportion = min_hmm_proportion/100

    hmm_length = refpkg.hmm_length()
    # Use the smaller of the minimum sequence length or some proportion of the profile HMM to remove sequence fragments
    if proposed_min_length < min_hmm_proportion*hmm_length:
        min_length = int(round(min_hmm_proportion*hmm_length))
        logging.debug("New minimum sequence length threshold set to {}% of HMM length ({}) instead of {}\n"
                      "".format(min_hmm_proportion*100, min_length, proposed_min_length))
    else:
        min_length = proposed_min_length

    return min_length


def intersect_incomparable_lists(superset, subset, name_map: dict) -> list:
    """
    Function for identifying the intersection of two lists by
    using a proxy identifier (alt_name) for elements in superset.

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


def drop_queries_by_accession(query_seqs: list, ref_seq_leaves: list):
    ref_seq_accessions = {leaf.accession for leaf in ref_seq_leaves}
    i = 0
    while i < len(query_seqs):
        if query_seqs[i].split()[0] in ref_seq_accessions:
            query_seqs.pop(i)
        else:
            i += 1
    return


def guided_header_lineage_map(header_registry: dict, entrez_records: dict) -> dict:
    """Generate a dictionary mapping sequence names to lineages guided by a header registry."""
    seq_lineage_map = {}
    missing = []
    # Generate a dictionary between the EntrezRecord sequence names and lineages, for sequences in the header_registry
    for treesapp_id in sorted(header_registry.keys(), key=int):
        try:
            record = entrez_records[treesapp_id]  # type: EntrezRecord
        except KeyError:
            # Log those that were not found
            missing.append(treesapp_id)
            continue
        seq_lineage_map[record.versioned] = record.lineage
    if missing:
        logging.warning(str(len(missing)) + " sequences were not assigned a taxonomic lineage.\n" +
                        "This should match the number of accessions deduplicated while fetching lineages.\n")
        for treesapp_id in missing:
            logging.debug("Unable to find '" + treesapp_id + "' in fasta records. More info:\n" +
                          header_registry[treesapp_id].original + "\n")
            header_registry.pop(treesapp_id)
        missing.clear()
    return seq_lineage_map


def simulate_entrez_records(fasta_records: FASTA, seq_lineage_map: dict) -> dict:
    """
    Creates new EntrezRecord instances for each sequence: lineage pair in seq_lineage_map.
    This function circumvents downloading lineage information for accessions, organism names or NCBI taxonomic IDs
    if the lineage is provided via other means, and enables compatibility with downstream functions.

    :param fasta_records: FASTA object containing sequences that can be mapped to
    :param seq_lineage_map: A dictionary mapping parsed sequence accessions to taxonomic lineages
    :return: A dictionary of EntrezRecord instances indexed by their respective TreeSAPP numerical IDs
    """
    entrez_records = dict()
    header_map = fasta_records.get_acc_ver_header_map()
    for seq_accession in sorted(seq_lineage_map):
        for header in header_map[seq_accession]:
            er = EntrezRecord(seq_accession, "")
            er.lineage = seq_lineage_map[seq_accession]
            er.organism = er.lineage.split("; ")[-1]
            er.description = " ".join(header.original.split(" ")[1:])
            er.versioned = header.original.split(" ")[0]
            er.sequence = fasta_records.fasta_dict[str(header.treesapp_num_id)]
            entrez_records[str(header.treesapp_num_id)] = er
    return entrez_records


def resolve_cluster_lineages(cluster_dict: dict, entrez_records: dict, taxa_trie: TaxonomicHierarchy) -> None:
    """
    Sets the 'cluster_rep' attribute to True for the EntrezRecord with the most resolved lineage out of the
    cluster members and cluster representative. If the cluster representative has a less resolved lineage its
    'cluster_rep' attribute is set to False and it is moved into the members list.

    :param entrez_records: Dictionary mapping unique TreeSAPP numerical IDs to Cluster instances
    :param cluster_dict: Dictionary mapping unique cluster IDs to Cluster instances
    :param taxa_trie: A TaxonomicHierarchy instance of the ReferencePackage being updated
    :return: None
    """
    # A temporary dictionary for rapid mapping of sequence names to lineages
    er_lookup = {er.rebuild_header(): er for (num_id, er) in entrez_records.items()}

    for cluster_id in cluster_dict:
        cluster = cluster_dict[cluster_id]  # type: Cluster
        ref_er = er_lookup[cluster.representative]  # type: EntrezRecord
        ref_depth = taxa_trie.accepted_ranks_depths[taxa_trie.resolved_to(ref_er.lineage)]
        if len(cluster.members) >= 1:
            validated_cluster_members = []
            for member in cluster.members:
                seq_name, seq_similarity = member
                member_er = er_lookup[seq_name]  # type: EntrezRecord
                member_depth = taxa_trie.accepted_ranks_depths[taxa_trie.resolved_to(member_er.lineage)]
                if member_depth <= ref_depth:
                    member_er.cluster_rep = False
                    validated_cluster_members.append(member)
                else:
                    validated_cluster_members.append((cluster.representative, seq_similarity))
                    ref_er.cluster_rep = False
                    cluster.representative = seq_name
                    ref_er = member_er
                    ref_depth = member_depth
            cluster.members = validated_cluster_members
    return


def prefilter_clusters(cluster_dict: dict, lineage_lookup: dict, priority: list, lineage_collapse=True) -> None:
    """
    Switches the representative sequence of a Cluster instance based on a priority list.

    Optionally, with the lineage_collapse flag, Cluster.members can be emptied if all members
    (including the representative) have identical taxonomic lineages.

    :param cluster_dict: Dictionary mapping unique cluster IDs to Cluster instances
    :param lineage_lookup: Dictionary mapping original sequence names to taxonomic lineages
    :param priority: List of sequences that should be centroids, if not already
    :param lineage_collapse: Flag indicating whether clusters whose members have identical lineages are removed
    :return: Sequence names in `priority` that were members of a cluster represented by another priority sequence.
    These can be used to identify which clusters should be broken such that all 'priority' sequences will be centroids
    """
    # cluster_ids list is used for iterating through dictionary keys and allowing dict to change size with 'pop's
    cluster_ids = sorted(list(cluster_dict.keys()), key=int)
    # Track the number of priority sequences that remained members of clusters
    guaranteed_redundant = []
    cluster_num = len(cluster_dict)

    for cluster_id in cluster_ids:
        cluster = cluster_dict[cluster_id]  # type: Cluster
        if len(cluster.members) == 0:
            continue
        # Ensure the centroids/representatives are the original reference sequences
        if cluster.representative in priority:
            rep_found = True
        else:
            rep_found = False
        i = 0
        while i < len(cluster.members):
            seq_name, seq_similarity = cluster.members[i]
            if seq_name in priority:
                if rep_found:
                    # Save the reference sequence from being absorbed into another reference sequence's cluster
                    cluster_break = Cluster(seq_name)
                    while str(cluster_num) in cluster_dict:
                        cluster_num += 1
                    guaranteed_redundant.append(cluster_break)
                    cluster_dict[str(cluster_num)] = cluster_break
                    cluster.members.pop(i)
                    i -= 1
                else:
                    cluster.members[i] = [cluster.representative, seq_similarity]
                    cluster.representative = seq_name
                    rep_found = True
            i += 1
        # Remove the cluster members from the dictionary if the lineages are identical
        if lineage_collapse:
            identical = True
            for member_seq in cluster.members:
                seq_name, seq_similarity = member_seq
                if lineage_lookup[seq_name] != cluster.lca:
                    identical = False
            if identical:
                cluster.members = []

    if guaranteed_redundant:
        logging.debug("{} original reference sequences saved from clustering:\n\t"
                      "{}\n".format(len(guaranteed_redundant),
                                    "\n\t".join(clust.representative for clust in guaranteed_redundant)))

    return


def formulate_create_command(ts_updater: Updater, args) -> list:
    create_cmd = ["-i", ts_updater.combined_fasta,
                  "-c", ts_updater.ref_pkg.prefix,
                  "-p", str(ts_updater.prop_sim),
                  "-m", ts_updater.molecule_type,
                  "--guarantee", ts_updater.old_ref_fasta,
                  "-o", ts_updater.output_dir,
                  "--accession2lin", ts_updater.lineage_map_file,
                  "--num_procs", str(args.num_threads),
                  "--bootstraps", str(args.bootstraps)]
    if args.trim_align:
        create_cmd.append("--trim_align")
    if args.od_seq:
        create_cmd.append("--outdet_align")
    if args.raxml_model:
        create_cmd += ["--raxml_model", args.raxml_model]
    if args.fast:
        create_cmd.append("--fast")
    if args.taxa_lca:
        create_cmd.append("--taxa_lca")
    if args.cluster:
        create_cmd.append("--cluster")
    if args.headless:
        create_cmd.append("--headless")
    if args.screen:
        create_cmd += ["--screen", args.screen]
    if args.filter:
        create_cmd += ["--filter", args.filter]
    if args.min_taxonomic_rank:
        create_cmd += ["--min_taxonomic_rank", args.min_taxonomic_rank]
    return create_cmd
