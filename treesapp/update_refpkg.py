import os
import re
import logging

from datetime import datetime as dt
from pygtrie import StringTrie

from treesapp.utilities import load_taxonomic_trie
from treesapp.classy import TreeSAPP, ModuleFunction
from treesapp.seq_clustering import Cluster
from treesapp.entrez_utils import EntrezRecord
from treesapp.fasta import FASTA
from treesapp.taxonomic_hierarchy import TaxonomicHierarchy
from treesapp import refpkg as ts_ref_pkg
from treesapp import clade_annotation
from treesapp import logger

LOGGER = logging.getLogger(logger.logger_name())


class Updater(TreeSAPP):
    def __init__(self):
        super(Updater, self).__init__("update")
        self.seq_names_to_taxa = ""  # Optional user-provided file mapping query sequence contigs to lineages
        self.lineage_map_file = ""  # File that is passed to create() containing lineage info for all sequences
        self.treesapp_output = ""  # Path to the TreeSAPP output directory - modified by args
        self.assignment_table = ""  # Path to the classifications.tsv file written by treesapp assign
        self.combined_fasta = ""  # Holds the newly identified candidate reference sequences and the original ref seqs
        self.old_ref_fasta = ""  # Contains only the original reference sequences
        self.cluster_input = ""  # Used only if resolve is True
        self.clusters_prefix = ""  # Used only if resolve is True
        self.updated_refpkg_path = ""
        self.training_dir = ""
        # self.rank_depth_map = None
        self.prop_sim = 1.0
        self.min_length = 0  # The minimum sequence length for a classified sequence to be included in the ref_pkg
        self.updated_refpkg = ts_ref_pkg.ReferencePackage()

        # Stage names only holds the required stages; auxiliary stages (e.g. RPKM, update) are added elsewhere
        self.stages = {0: ModuleFunction("lineages", 0),
                       1: ModuleFunction("rebuild", 1),
                       2: ModuleFunction("train", 2)}

    def decide_stage(self, args):
        """
        Bases the stage(s) to run on args.stage which is broadly set to either 'continue' or any other valid stage

        This function ensures all the required inputs are present for beginning at the desired first stage,
        otherwise, the pipeline begins at the first possible stage to continue and ends once the desired stage is done.

        :return: None
        """
        self.validate_continue(args)
        self.acc_to_lin = self.var_output_dir + os.sep + "accession_id_lineage_map.tsv"

        if not self.input_sequences:
            self.change_stage_status("train", False)

        return

    def get_info(self):
        info_string = "Updater instance summary:\n"
        info_string += super(Updater, self).get_info() + "\n\t"
        info_string += "\n\t".join(["Target reference packages = " + str(self.ref_pkg.prefix),
                                    "Taxonomy map: " + self.ref_pkg.lineage_ids,
                                    "Reference tree: " + self.ref_pkg.tree,
                                    "Reference FASTA: " + self.ref_pkg.msa,
                                    "Lineage map: " + str(self.seq_names_to_taxa)]) + "\n"

        return info_string

    def update_refpkg_fields(self, output_dir=None) -> None:
        """
        Using the original ReferencePackage as a template modify the following updated ReferencePackage attributes:
1. original creation date
2. update date
3. code name
4. description
        :return: None
        """
        if not output_dir:
            output_dir = self.final_output_dir
        # Change the creation and update dates, code name and description
        self.updated_refpkg.change_file_paths(output_dir)
        self.updated_refpkg.date = self.ref_pkg.date
        self.updated_refpkg.update = dt.now().strftime("%Y-%m-%d")
        self.updated_refpkg.refpkg_code = self.ref_pkg.refpkg_code
        self.updated_refpkg.description = self.ref_pkg.description
        self.updated_refpkg.pickle_package()

        self.ts_logger.info("Summary of the updated reference package:\n" + self.updated_refpkg.get_info() + "\n")

        self.ts_logger.debug("\tNew sequences  = " + str(self.updated_refpkg.num_seqs - self.ref_pkg.num_seqs) + "\n" +
                             "\tOld HMM length = " + str(self.ref_pkg.hmm_length()) + "\n" +
                             "\tNew HMM length = " + str(self.updated_refpkg.hmm_length()) + "\n")

        return


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


def filter_by_placement_thresholds(pqueries: dict, min_lwr: float, max_pendant=None, max_evo_distance=None) -> list:
    """

    :param pqueries: A dictionary of PQuery instances indexed by their respective ReferencePackage's
    :param min_lwr: Minimum acceptable Likelihood Weight Ratio as calculated by EPA-NG (float)
    :param max_pendant: Maximum pendant length for a PQuery as calculated by EPA-NG (float)
    :param max_evo_distance: Maximum total evolutionary distance length for a PQuery (float)
    :return: A set of PQuery instances that passed all phylogenetic placement thresholds
    """
    good_placements = list()
    num_filtered = 0  # Number of placements removed
    for refpkg_pqueries in pqueries.values():
        for pquery in refpkg_pqueries:
            if pquery.consensus_placement.like_weight_ratio < min_lwr:
                num_filtered += 1
            elif max_pendant and pquery.consensus_placement.pendant_length > max_pendant:
                num_filtered += 1
            elif max_evo_distance and pquery.avg_evo_dist > max_evo_distance:
                num_filtered += 1
            else:
                good_placements.append(pquery)

    LOGGER.debug("{} classified sequences did not meet minimum LWR of {} for updating\n".format(num_filtered, min_lwr))

    return good_placements


def decide_length_filter(ref_pkg: ts_ref_pkg.ReferencePackage, proposed_min_length=0, min_hmm_proportion=0.66) -> int:
    # Ensure the HMM fraction provided is a proportion
    if not 0 < min_hmm_proportion < 1:
        LOGGER.warning("Minimum HMM fraction provided ({}) isn't a proportion. Converting to {}.\n"
                       "".format(min_hmm_proportion, min_hmm_proportion/100))
        min_hmm_proportion = min_hmm_proportion/100

    hmm_length = ref_pkg.hmm_length()
    if proposed_min_length > 0:
        min_length = proposed_min_length
        if min_length > hmm_length:
            LOGGER.warning("Minimum sequence length () is greater than profile HMM length {}. "
                           "This may cause an inordinate number of sequences to be filtered out."
                           "\n".format(min_length, hmm_length))
    # Use the smallest of the minimum sequence length or some proportion of the profile HMM to remove sequence fragments
    else:
        min_length = int(round(min_hmm_proportion*hmm_length))

    LOGGER.debug("Minimum sequence length threshold set to {} ({}% of HMM length {}).\n"
                 "".format(min_length, 100*min_length/hmm_length, hmm_length))

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
        LOGGER.warning(str(len(missing)) + " sequences were not assigned a taxonomic lineage.\n" +
                        "This should match the number of accessions deduplicated while fetching lineages.\n")
        for treesapp_id in missing:
            LOGGER.debug("Unable to find '" + treesapp_id + "' in fasta records. More info:\n" +
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
            er.sequence = fasta_records.fasta_dict[str(header.num_id)]
            entrez_records[str(header.num_id)] = er
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
        LOGGER.debug("{} original reference sequences saved from clustering:\n\t"
                      "{}\n".format(len(guaranteed_redundant),
                                    "\n\t".join(clust.representative for clust in guaranteed_redundant)))

    return


def formulate_create_command(ts_updater: Updater, args, final_stage) -> list:
    create_cmd = ["-i", ts_updater.combined_fasta,
                  "-c", ts_updater.ref_pkg.prefix,
                  "-p", str(ts_updater.prop_sim),
                  "-m", ts_updater.molecule_type,
                  "--guarantee", ts_updater.old_ref_fasta,
                  "-o", ts_updater.output_dir,
                  "--accession2lin", ts_updater.lineage_map_file,
                  "--num_procs", str(args.num_threads),
                  "--bootstraps", str(args.bootstraps),
                  "--stage", final_stage]
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


def update_features(old_refpkg: ts_ref_pkg.ReferencePackage, new_refpkg: ts_ref_pkg.ReferencePackage) -> None:
    # Match leaf node identifiers to each other
    leaf_node_name_map = {}
    old_leaf_desc_map = {leaf.description: leaf for leaf in old_refpkg._ref_leaves.values()}
    # Dictionary is indexed by the new leaf node names
    for desc, leaf in {ln.description: ln for ln in new_refpkg._ref_leaves.values()}.items():
        if desc in old_leaf_desc_map:
            leaf_node_name_map[old_leaf_desc_map[desc].number + '_' + old_refpkg.prefix] = leaf.number + '_' + new_refpkg.prefix

    for feature in old_refpkg.feature_annotations:
        new_refpkg.feature_annotations[feature] = []
        for clade_annot in old_refpkg.feature_annotations[feature]:  # type: clade_annotation.CladeAnnotation
            new_annotation = clade_annotation.CladeAnnotation(name=clade_annot.name, key=feature)
            # try:
            for leaf_node, rank in clade_annot.members.items():
                new_annotation.members[leaf_node_name_map[leaf_node]] = rank
            # except AttributeError:
            #     for leaf_node in clade_annot.members:
            #         new_annotation.members[leaf_node_name_map[leaf_node]] = 7
            new_annotation.taxa = clade_annot.taxa
            new_annotation.colour = clade_annot.colour
            new_refpkg.feature_annotations[feature].append(new_annotation)
    return
