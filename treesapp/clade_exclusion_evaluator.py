import os
import logging
import sys
import re
import shutil

from treesapp import refpkg
from treesapp import entrez_utils as ez_utils
from treesapp import fasta
from treesapp import classy
from treesapp import taxonomic_hierarchy
from treesapp import assign
from treesapp import file_parsers
from treesapp import logger
from treesapp import graftm_utils

LOGGER = logging.getLogger(logger.logger_name())


def parse_distances(classification_lines):
    distances = dict()
    distances["distal"] = list()
    distances["pendant"] = list()
    distances["tip"] = list()
    for fields in classification_lines:
        dist_fields = fields[-1].split(',')
        distances["distal"].append(float(dist_fields[0]))
        distances["pendant"].append(float(dist_fields[1]))
        distances["tip"].append(float(dist_fields[2]))
    return distances


def determine_containment(marker_eval_inst: classy.Evaluator):
    """
    Determines the accuracy of sequence classifications of all sequences contained at different taxonomic ranks

    :param marker_eval_inst: An instance of the MarkerTest class.
        The following instance variables need to be populated:
        marker: Name of the marker currently being evaluated (e.g., nifHc1, mcrA)
        classifications: The rank-wise classification dictionary with
        tuples of (optimal assignment, actual assignment) as values.
        E.g. {"Phylum": {"Proteobacteria": ("Proteobacteria", "Proteobacteria; Alphaproteobacteria")}}
    """
    rank_assigned_dict = marker_eval_inst.classifications
    # Set up collection for this analysis
    lowest_rank = ""
    for depth in sorted(marker_eval_inst.rank_depth_map):
        if marker_eval_inst.rank_depth_map[depth] in marker_eval_inst.ranks:
            if marker_eval_inst.get_sensitivity(marker_eval_inst.rank_depth_map[depth])[1] > 0:
                lowest_rank = marker_eval_inst.rank_depth_map[depth]
    containment_strings = list()
    n_queries, n_classified, _ = marker_eval_inst.get_sensitivity(lowest_rank)

    sys.stdout.write("Accuracy of " + str(n_classified) + ' ' + marker_eval_inst.ref_pkg.prefix +
                     " classifications at " + lowest_rank + ":\n" +
                     "\tRank\tCorrect\tToo Shallow (%)\n")

    # Begin parsing through the depths
    if lowest_rank == "root":
        return
    for depth in sorted(marker_eval_inst.rank_depth_map):
        incorrect_assignments = dict()
        rank = marker_eval_inst.rank_depth_map[depth]
        if rank in ["root", "species", "strain"]:
            continue
        elif rank == lowest_rank:
            break
        correct = 0
        incorrect = 0
        too_shallow = 0
        parse_depth = depth - 1
        for assignments in rank_assigned_dict[lowest_rank]:
            for classified in assignments:
                status = 0  # 0 == incorrect, 1 == correct
                optimal, query = assignments[classified]
                classified_lineage = classified.split("; ")
                query_lineage = query.split("; ")
                # print("\nReference:\t", classified)
                # print("Query:\t\t", query)
                # print("Optimal:\t", optimal)
                if len(classified_lineage) > parse_depth:
                    try:
                        if classified_lineage[parse_depth] == query_lineage[parse_depth]:
                            status += 1
                        else:
                            if query not in incorrect_assignments.keys():
                                incorrect_assignments[query] = 0
                            incorrect_assignments[query] += 1
                            incorrect += 1
                    except IndexError:
                        too_shallow += 1
                        # This indicates TreeSAPP placed the sequence too deep in the tree according to its NCBI lineage
                elif len(optimal.split("; ")) > parse_depth:
                    if query not in incorrect_assignments.keys():
                        incorrect_assignments[query] = 0
                    incorrect_assignments[query] += 1
                    incorrect += 1
                else:
                    too_shallow += 1
                correct += status

        percentage_too_shallow = round(float((too_shallow * 100) / n_classified), 1)
        if percentage_too_shallow == 100.0:
            containment_strings.append("\t" + rank + "\t" + str(0.0) + "\t" + str(100.0))
        else:
            containment_strings.append("\t" + rank + "\t" +
                                       str(round((correct * 100) / (correct + incorrect), 1)) + "\t" +
                                       str(percentage_too_shallow))

    sys.stdout.write("\n".join(containment_strings) + "\n")

    return containment_strings


def map_headers_to_lineage(assignments: dict, ref_sequences: dict) -> dict:
    """
    The alternative function to map_full_headers. Using ReferenceSequence objects,

    :param assignments: A dictionary of ReferencePackage.prefix keys mapped to a dictionary of assigned lineage keys
    mapped to the headers of sequences assigned as their respective lineages.
    :param ref_sequences:
    :return:  A dictionary mapping ReferencePackage.prefix's to a dictionary of assigned lineage keys to a list of
    reference sequences (values). An example is
    """
    lineage_assignments = dict()
    for refpkg_name in assignments:
        lineage_assignments[refpkg_name] = dict()
        for assigned_lineage in assignments[refpkg_name].keys():
            classified_headers = assignments[refpkg_name][assigned_lineage]
            lineage_assignments[refpkg_name][assigned_lineage] = list()
            for query in classified_headers:
                mapped = False
                for _, ref_seq in ref_sequences.items():  # type: (int, ez_utils.EntrezRecord)
                    if ref_seq.accession == query:
                        lineage_assignments[refpkg_name][assigned_lineage].append(ref_seq.lineage)
                        mapped = True
                        break
                if not mapped:
                    LOGGER.error("Unable to map classified sequence '" + query + "' to a lineage.\n")
                    sys.exit(3)
            if len(lineage_assignments[refpkg_name][assigned_lineage]) > len(classified_headers):
                LOGGER.error(str(len(classified_headers)) + " accessions mapped to " +
                             str(len(lineage_assignments[refpkg_name][assigned_lineage])) + " lineages.\n")
                sys.exit(21)
            elif len(lineage_assignments[refpkg_name][assigned_lineage]) < len(classified_headers):
                LOGGER.debug(str(len(classified_headers)) + " accessions mapped to " +
                             str(len(lineage_assignments[refpkg_name][assigned_lineage])) + " lineages.\n")
    return lineage_assignments


def check_lineage_compatibility(lineage_map: dict, taxonomy: taxonomic_hierarchy.TaxonomicHierarchy) -> None:
    """Ensures that all the lineages are rooted."""
    unrooted = False
    for num_id, lineage in lineage_map.items():
        # Ensure the query lineages are rooted and emit a warning if this is not the case
        query_split = lineage.split(taxonomy.lin_sep)
        if query_split[0] != taxonomy.root_taxon:
            unrooted = True
            lineage_map[num_id] = taxonomy.root_taxon + taxonomy.lin_sep + lineage
    if unrooted:
        LOGGER.warning("Unrooted lineages have been rooted.\n")
    return


def get_testable_lineages_for_rank(ref_lineage_map: dict, query_lineage_map: dict, rank: str,
                                   taxonomy: taxonomic_hierarchy.TaxonomicHierarchy) -> list:
    """
    Clade exclusion analysis requires a set of query sequences that can be compared to the reference sequences at some
    taxonomic rank. Critically, the query sequences need to be represented by the reference sequences at the rank.
    This function determines which query lineages are represented and returns them as a list.

    :param ref_lineage_map: A dictionary mapping all reference sequences to their respective taxonomic lineages
    :param query_lineage_map: A dictionary mapping representative query sequences to taxonomic lineages
    :param rank: Name of a taxonomic rank to test and is used to guide taxonomic lineage trimming with
     {"Kingdom": 1, "Phylum": 2, "Class": 3, "Order": 4, "Family": 5, "Genus": 6, "Species": 7}
    :param taxonomy: A TaxonomicHierarchy instance with all reference and query lineages loaded
    :return: List of lineages where the optimal rank exists in the reference tree after clade exclusion
    """
    lineages = []
    leaf_trimmed_taxa_map = taxonomy.trim_lineages_to_rank(ref_lineage_map, rank)
    unique_ref_lineages = sorted(set(leaf_trimmed_taxa_map.values()))
    unique_query_lineages = sorted(set(taxonomy.trim_lineages_to_rank(query_lineage_map, rank).values()))
    for lineage in unique_query_lineages:
        # Is the optimal placement in the pruned reference tree?
        optimal_lca_taxonomy = "; ".join(lineage.split("; ")[:-1])
        if optimal_lca_taxonomy not in ["; ".join(tl.split("; ")[:-1]) for tl in unique_ref_lineages
                                        if tl != lineage]:
            LOGGER.debug("Optimal placement target '" + optimal_lca_taxonomy + "' not in pruned tree.\n")
        else:
            lineages.append(lineage)
    return lineages


def pick_taxonomic_representatives(ref_seqs: dict, taxonomic_filter_stats: dict, max_cluster_size=5):
    """
    Removes queries with duplicate taxa - to prevent the taxonomic composition of the input
    from biasing the accuracy to over- or under-perform by classifying many sequences from very few groups.
    Also removes taxonomies with "*nclassified" in their lineage or are derived from environmental samples

    :param ref_seqs: A dictionary mapping accessions to lineages that need to be filtered
    :param taxonomic_filter_stats: A dictionary for tracking the number sequences filtered, those unique, etc.
    :param max_cluster_size: The maximum number of sequences representing a taxonomic cluster
    :return: dereplicated_lineages dict with lineages mapping to a (short) list of accessions
    """
    good_classified_lineages = dict()
    dereplicated_lineages = dict()
    num_rep_seqs = 0
    for tree_id, ref_seq in ref_seqs.items():
        if len(ref_seq.lineage.split("; ")) < 5:
            continue
        if ref_seq.lineage not in good_classified_lineages:
            good_classified_lineages[ref_seq.lineage] = list()
        if re.search("unclassified|environmental sample", ref_seq.lineage, re.IGNORECASE):
            # # Remove taxonomic lineages that are unclassified at the Phylum level or higher
            # unclassified_depth = get_unclassified_rank(0, query_taxonomy.split("; "))
            # if unclassified_depth > 4:
            #     good_classified_lineages[query_taxonomy].append(ref_seq.accession)
            # else:
            taxonomic_filter_stats["Unclassified"] += 1
        else:
            taxonomic_filter_stats["Classified"] += 1
            good_classified_lineages[ref_seq.lineage].append(ref_seq.accession)

    if taxonomic_filter_stats["Unclassified"] == len(ref_seqs):
        LOGGER.error("All sequences provided are derived from uncultured, unclassified organisms.\n")
        sys.exit(21)

    # In order to maintain consistency among multiple runs with the same input
    # a separate loop is required for sorting
    taxonomic_filter_stats["Max"] = 0
    taxonomic_filter_stats["Min"] = 1000
    taxonomic_filter_stats["Mean"] = 0
    for query_taxonomy in good_classified_lineages:
        dereplicated_lineages[query_taxonomy] = list()
        i = 0
        lineage_candidates = sorted(good_classified_lineages[query_taxonomy])
        while i < len(lineage_candidates) and i < max_cluster_size:
            dereplicated_lineages[query_taxonomy].append(lineage_candidates[i])
            i += 1

        # Generate stats on the taxonomic clusters
        cluster_size = len(dereplicated_lineages[query_taxonomy])
        num_rep_seqs += cluster_size
        if cluster_size == 0:
            # This is a taxonomic lineage removed by the "Unclassified" filter
            continue
        if cluster_size > taxonomic_filter_stats["Max"]:
            taxonomic_filter_stats["Max"] = cluster_size
        if cluster_size < taxonomic_filter_stats["Min"]:
            taxonomic_filter_stats["Min"] = cluster_size
        taxonomic_filter_stats["Mean"] = round(float(num_rep_seqs / len(dereplicated_lineages)), 2)

    taxonomic_filter_stats["Unique_taxa"] += len(dereplicated_lineages)

    LOGGER.debug("Representative sequence stats:\n\t" +
                 "Maximum representative sequences for a taxon " + str(taxonomic_filter_stats["Max"]) + "\n\t" +
                 "Minimum representative sequences for a taxon " + str(taxonomic_filter_stats["Min"]) + "\n\t" +
                 "Mean representative sequences for a taxon " + str(taxonomic_filter_stats["Mean"]) + "\n")

    return dereplicated_lineages, taxonomic_filter_stats


def select_rep_seqs(deduplicated_assignments: dict, test_sequences: dict,
                    taxon_hierarchy: taxonomic_hierarchy.TaxonomicHierarchy, target_lineage=None) -> dict:
    """
    Function for creating a fasta-formatted dict from the accessions representing unique taxa in the test sequences

    :param deduplicated_assignments: dict of lineages mapped to a list of accessions
    :param test_sequences: Dictionary of ReferenceSequence objects indexed by their numerical identifiers
    :param taxon_hierarchy: A TaxonomicHierarchy instance
    :param target_lineage: A taxonomic lineage to filter the lineages (and sequences) in deduplicated assignments
    :return: Dictionary containing accessions as keys and sequences as values
    """
    if target_lineage:
        descendents = set()
        filtered_assignments = dict()
        # Gather the set of all descendent lineages
        target_taxon = target_lineage.split(taxon_hierarchy.lin_sep)[-1]
        for line in [taxon.lineage() for taxon in taxon_hierarchy.get_taxon_descendents(target_taxon)]:  # type: list
            descendents.add(taxon_hierarchy.lin_sep.join([taxon.prefix_taxon() for taxon in line]))
        # Query the candidate lineages for whether they're in the target lineage's list of descendents
        for candidate_lineage in sorted(deduplicated_assignments):
            if candidate_lineage in descendents:
                filtered_assignments[candidate_lineage] = deduplicated_assignments[candidate_lineage]
    else:
        filtered_assignments = deduplicated_assignments

    test_seq_lookup = {ref_seq.accession: ref_seq.sequence for ref_seq in test_sequences.values()}

    deduplicated_fasta_dict = dict()
    for lineage in sorted(filtered_assignments):
        for accession in filtered_assignments[lineage]:
            try:
                deduplicated_fasta_dict[accession] = test_seq_lookup[accession]
            except KeyError:
                LOGGER.error("Unable to find accession (" + accession + ") in accession-lineage map\n")
                sys.exit(21)
    return deduplicated_fasta_dict


def map_seqs_to_lineages(accession_lineage_map, deduplicated_fasta_dict):
    seq_taxa_map = dict()
    for seq_name in deduplicated_fasta_dict:
        if seq_name not in accession_lineage_map:
            LOGGER.error("Unable to find matching key for '" + seq_name + "' in accession_lineage_map.\n")
            sys.exit(21)
        else:
            seq_taxa_map[seq_name] = accession_lineage_map[seq_name]
    return seq_taxa_map


def run_clade_exclusion_treesapp(tt_obj: classy.TaxonTest, taxon_rep_seqs: dict, ref_pkg: refpkg.ReferencePackage,
                                 cl_assign_args, executables, fresh=False, min_seq_length=30) -> None:
    if not os.path.isfile(tt_obj.classification_table):
        # Copy reference files, then exclude all clades belonging to the taxon being tested
        ce_refpkg = ref_pkg.clone(tt_obj.refpkg_path)
        ce_refpkg.exclude_clade_from_ref_files(tt_obj.intermediates_dir, tt_obj.lineage,
                                               executables, fresh)
        # Write the query sequences
        fasta.write_new_fasta(taxon_rep_seqs, tt_obj.test_query_fasta)
        assign_args = ["-i", tt_obj.test_query_fasta, "-o", tt_obj.classifications_root,
                       "--refpkg_dir", os.path.dirname(ce_refpkg.f__pkl),
                       "-n", str(cl_assign_args.num_threads),
                       "--min_seq_length", str(min_seq_length),
                       "--placement_summary", cl_assign_args.p_sum,
                       "--overwrite", "--delete", "--silent"]
        if cl_assign_args.trim_align:
            assign_args.append("--trim_align")
        if cl_assign_args.max_evo:
            assign_args += ["--max_evol_distance", cl_assign_args.max_evo]
        if cl_assign_args.max_pd:
            assign_args += ["--max_pendant_length", cl_assign_args.max_pd]
        if cl_assign_args.min_lwr:
            assign_args += ["--min_like_weight_ratio", str(cl_assign_args.min_lwr)]

        try:
            assign.assign(assign_args)
        except:  # Just in case treesapp assign fails, just continue
            pass
        if not os.path.isfile(tt_obj.classification_table):
            # The TaxonTest object is maintained for record-keeping (to track # queries & classifieds)
            LOGGER.warning("TreeSAPP did not generate output for '{}'. Skipping.\n".format(tt_obj.lineage))
            shutil.rmtree(tt_obj.classifications_root)
            return
    else:
        # Valid number of queries and these sequences have already been classified
        ce_refpkg = refpkg.ReferencePackage()
        ce_refpkg.f__pkl = os.path.join(tt_obj.intermediates_dir,
                                        ref_pkg.refpkg_code,
                                        ref_pkg.prefix + ref_pkg.refpkg_suffix)
        ce_refpkg.slurp()

    tt_obj.taxonomic_tree = ce_refpkg.all_possible_assignments()
    if os.path.isfile(tt_obj.classification_table):
        assigned_lines = file_parsers.read_classification_table(tt_obj.classification_table)
        tt_obj.assignments = file_parsers.parse_assignments(assigned_lines)
        tt_obj.filter_assignments(ref_pkg.prefix)
        tt_obj.distances = parse_distances(assigned_lines)
    else:
        LOGGER.error("{} is missing from output directory '{}'\n"
                     "Please remove this directory and re-run.\n"
                     "".format(os.path.basename(tt_obj.classification_table),
                               os.path.dirname(tt_obj.classification_table)))
        sys.exit(21)
    return


def run_clade_exclusion_graftm(tt_obj: classy.TaxonTest, taxon_rep_seqs, ref_pkg: refpkg.ReferencePackage,
                               graftm_classifier: str, num_threads: int, executables: dict) -> None:
    # GraftM refpkg output files:
    gpkg_refpkg_path = os.path.join(tt_obj.refpkg_path,
                                    ref_pkg.prefix +
                                    '_' + re.sub(' ', '_', tt_obj.taxon) + ".gpkg.refpkg") + os.sep
    gpkg_tax_ids_file = gpkg_refpkg_path + ref_pkg.prefix + "_taxonomy.csv"
    if not os.path.isfile(tt_obj.classification_table):
        # Copy reference files, then exclude all clades belonging to the taxon being tested
        output_paths = graftm_utils.prep_graftm_ref_files(tmp_dir=tt_obj.intermediates_dir,
                                                          target_clade=tt_obj.lineage,
                                                          ref_pkg=ref_pkg,
                                                          executables=executables)

        if not os.path.isdir(tt_obj.refpkg_path):
            graftm_utils.build_graftm_package(gpkg_path=tt_obj.refpkg_path,
                                              tax_file=output_paths["filtered_tax_ids"],
                                              mfa_file=output_paths["filtered_mfa"],
                                              fa_file=output_paths["filtered_fasta"],
                                              threads=num_threads)
        # Write the query sequences
        fasta.write_new_fasta(taxon_rep_seqs, tt_obj.test_query_fasta)

        graftm_utils.run_graftm_graft(input_path=tt_obj.test_query_fasta, gpkg_path=tt_obj.refpkg_path,
                                      output_dir=tt_obj.classifications_root,
                                      num_threads=num_threads, classifier=graftm_classifier)

        if not os.path.isfile(tt_obj.classification_table):
            # The TaxonTest object is maintained for record-keeping (to track # queries & classifieds)
            LOGGER.warning("GraftM did not generate output for " + tt_obj.lineage + ". Skipping.\n")
            shutil.rmtree(tt_obj.intermediates_dir)

    tt_obj.taxonomic_tree = graftm_utils.grab_graftm_taxa(gpkg_tax_ids_file)
    graftm_assignments = graftm_utils.read_graftm_classifications(tt_obj.classification_table)
    tt_obj.assignments = {ref_pkg.prefix: graftm_assignments}
    tt_obj.filter_assignments(ref_pkg.prefix)

    return
