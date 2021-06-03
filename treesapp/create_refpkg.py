import logging
import sys
import re
import os

from treesapp import seq_clustering
from treesapp.wrapper import run_odseq, run_mafft
from treesapp.lca_calculations import megan_lca, clean_lineage_list
from treesapp.taxonomic_hierarchy import TaxonomicHierarchy
from treesapp import entrez_utils
from treesapp import fasta
from treesapp import classy
from treesapp import logger

LOGGER = logging.getLogger(logger.logger_name())


def check_create_arguments(creator: classy.Creator, args) -> None:
    creator.combine_input_files()
    creator.find_sequence_molecule_type()
    # Populate ReferencePackage attributes from command-line arguments
    if args.fast:
        creator.ref_pkg.tree_tool = "FastTree"
    else:
        creator.ref_pkg.tree_tool = "RAxML-NG"
    creator.ref_pkg.prefix = args.refpkg_name
    creator.ref_pkg.pid = args.similarity
    creator.ref_pkg.molecule = creator.molecule_type
    creator.ref_pkg.kind = args.kind
    creator.ref_pkg.sub_model = args.raxml_model
    creator.ref_pkg.f__pkl = creator.final_output_dir + creator.ref_pkg.prefix + creator.ref_pkg.refpkg_suffix
    # TODO: Create placement trainer output directory and make it an attribute
    if not args.output:
        args.output = os.getcwd() + os.sep + creator.ref_pkg.prefix + "_treesapp_refpkg" + os.sep

    if len(creator.ref_pkg.prefix) > 10:
        LOGGER.error("Name should be <= 10 characters.\n")
        sys.exit(13)

    # TODO: Check the substitution model for compatibility with RAxML-NG

    if args.cluster:
        if args.multiple_alignment:
            LOGGER.error("--cluster and --multiple_alignment are mutually exclusive!\n")
            sys.exit(13)
        if not 0.5 <= float(args.similarity) <= 1.0:
            if 0.5 < float(args.similarity) / 100 < 1.0:
                args.similarity = str(float(args.similarity) / 100)
                LOGGER.warning("--similarity  set to {} for compatibility with VSEARCH.\n".format(args.similarity))
            else:
                LOGGER.error("--similarity {} is not between the supported range [0.5-1.0].\n".format(args.similarity))
                sys.exit(13)

    if args.taxa_lca and not args.cluster:
        LOGGER.error("Unable to perform LCA for representatives without clustering information: " +
                      "either with a provided VSEARCH file or by clustering within the pipeline.\n")
        sys.exit(13)

    if args.guarantee and not args.cluster:
        LOGGER.error("--guarantee used but without clustering there is no reason for it.\n" +
                      "Include all sequences in " + args.guarantee +
                      " in " + creator.input_sequences + " and re-run without --guarantee\n")
        sys.exit(13)

    if args.profile:
        if not os.path.isfile(args.profile):
            LOGGER.error("Unable to find HMM profile at '" + args.profile + "'.\n")
            sys.exit(3)
        creator.hmm_profile = args.profile

    # Names of files and directories to be created
    creator.phy_dir = os.path.abspath(creator.var_output_dir) + os.sep + "phylogeny_files" + os.sep
    creator.training_dir = os.path.abspath(creator.var_output_dir) + os.sep + "placement_trainer" + os.sep
    creator.hmm_purified_seqs = creator.var_output_dir + creator.ref_pkg.prefix + "_hmm_purified.fasta"
    creator.filtered_fasta = creator.var_output_dir + creator.sample_prefix + "_filtered.fa"
    creator.cluster_input = creator.var_output_dir + creator.sample_prefix + "_cluster_input.fasta"
    creator.clusters_prefix = creator.var_output_dir + creator.sample_prefix + "_cluster" + str(creator.ref_pkg.pid)
    creator.unaln_ref_fasta = creator.var_output_dir + creator.ref_pkg.prefix + "_ref.fa"
    creator.phylip_file = creator.var_output_dir + creator.ref_pkg.prefix + ".phy"

    # Ensure the phylogenetic tree output directory from a previous run isn't going to be over-written
    if not os.path.exists(creator.phy_dir):
        os.mkdir(creator.phy_dir)
    else:
        LOGGER.error(creator.phy_dir + " already exists from a previous run! " +
                      "Please delete or rename it and try again.\n")
        sys.exit(13)

    if not os.path.isdir(creator.training_dir):
        os.mkdir(creator.training_dir)

    return


def create_new_ref_fasta(out_fasta, ref_seq_dict, dashes=False):
    """
    Writes a new FASTA file using a dictionary of ReferenceSequence class objects

    :param out_fasta: Name of the FASTA file to write to
    :param ref_seq_dict: Dictionary containing ReferenceSequence objects, numbers are keys
    :param dashes: Flag indicating whether hyphens should be retained in sequences
    :return:
    """
    out_fasta_handle = open(out_fasta, "w")
    num_seqs_written = 0

    for treesapp_id in sorted(ref_seq_dict, key=int):
        ref_seq = ref_seq_dict[treesapp_id]
        if dashes is False:
            sequence = re.sub('[-.]', '', ref_seq.sequence)
        else:
            # sequence = re.sub('\.', '', ref_seq.sequence)
            sequence = ref_seq.sequence
        out_fasta_handle.write(">" + ref_seq.short_id + "\n" + sequence + "\n")
        num_seqs_written += 1

    out_fasta_handle.close()

    if num_seqs_written == 0:
        LOGGER.error("No sequences written to " + out_fasta + ".\n" +
                      "The headers in your input file are probably not accommodated in the regex patterns used. " +
                      "Function responsible: get_header_format. Please make an issue on the GitHub page.\n")
        sys.exit(5)

    return


def finalize_cluster_reps(cluster_dict: dict, refseq_objects: dict, header_registry: dict) -> None:
    """
    Transfer information from the cluster data (representative sequence, identity and cluster taxonomic LCA) to the
    dictionary of ReferenceSequence objects. The sequences not representing a cluster will have their `cluster_rep`
    flags set to *False* so as to not be analyzed further.

    :param cluster_dict: A dictionary of unique cluster IDs mapped to Cluster instances
    :param refseq_objects: A dictionary of numerical TreeSAPP IDs mapped to EntrezRecord instances
    :param header_registry: A dictionary of Header() objects indexed by TreeSAPP IDs,
     each used to map various header formats to each other
    :return: Dictionary of ReferenceSequence objects with complete clustering information
    """
    LOGGER.debug("Finalizing representative sequence clusters... ")

    cluster_reps = dict()  # A set of unique sequence names (original headers) to rapidly query
    for cluster_id in cluster_dict:
        cluster_reps[cluster_dict[cluster_id].representative] = cluster_dict[cluster_id]
    for treesapp_id in header_registry:
        header = header_registry[treesapp_id]  # type: fasta.Header
        try:
            ref_seq = refseq_objects[treesapp_id]  # type: entrez_utils.EntrezRecord
        except KeyError:
            continue  # Sequence was likely removed
        if header.original not in cluster_reps:
            ref_seq.cluster_rep = False
        else:
            ref_seq.cluster_lca = cluster_reps[header.original].lca

    LOGGER.debug("done.\n")
    return


def present_cluster_rep_options(cluster_dict: dict, refseq_objects: dict, header_registry: dict,
                                default=0, important_seqs=None, each_lineage=False) -> None:
    """
    Present the headers of identical sequences to user for them to decide on representative header

    :param cluster_dict: dictionary from read_uc(uc_file)
    :param refseq_objects: A dictionary of numerical TreeSAPP IDs mapped to EntrezRecord instances
    :param header_registry: A list of Header() objects, each used to map various header formats to each other
    :param default: The default cluster member number to use as the cluster representative.
        Running default set to 1 is equivalent to using the --headless flag.
    :param important_seqs: If --guarantee is provided, a dictionary mapping headers to seqs from format_read_fasta()
    :param each_lineage: If set to True, each candidate's lineage is shown as well as the cluster's LCA
    :return: None
    """
    if not important_seqs:
        important_seqs = set()
    candidates = dict()
    seq_names = list(refseq_objects.keys())
    for cluster_id in sorted(cluster_dict, key=int):
        cluster_info = cluster_dict[cluster_id]
        acc = 1
        candidates.clear()
        for num_id in sorted(refseq_objects, key=int):
            if header_registry[num_id].original == cluster_info.representative:
                refseq_objects[num_id].cluster_rep_similarity = '*'
                refseq_objects[num_id].cluster_lca = cluster_info.lca
                candidates[str(acc)] = refseq_objects[num_id]
                acc += 1
                break
        if acc != 2:
            raise AssertionError("Unable to find " + cluster_info.representative + " in ReferenceSequence objects!")

        if len(cluster_info.members) >= 1 and cluster_info.representative not in important_seqs:
            if len(cluster_info.members) == 1:
                member_id, perc_id = cluster_info.members[0]
                if cluster_info.representative == member_id:
                    continue
            # Find the EntrezRecords corresponding to each member so they can be displayed
            for cluster_member_info in cluster_info.members:
                if cluster_member_info[0] == cluster_info.representative:
                    continue
                for treesapp_id in sorted(seq_names, key=int):
                    if header_registry[treesapp_id].original == cluster_member_info[0]:
                        refseq_objects[treesapp_id].cluster_rep_similarity = cluster_member_info[1]
                        candidates[str(acc)] = refseq_objects[treesapp_id]
                        acc += 1
                        seq_names.pop(seq_names.index(treesapp_id))
                        break
            sys.stderr.write("Sequences in '" + cluster_info.lca + "' cluster:\n")
            for num in sorted(candidates.keys(), key=int):
                sys.stderr.write("\t{}. ".format(num) +
                                 '\t'.join(["{} | {}\t".format(candidates[num].organism, candidates[num].accession),
                                            str(len(candidates[num].sequence)) + " bp or aa",
                                            str(candidates[num].cluster_rep_similarity)]))
                if each_lineage:
                    sys.stderr.write("\t" + "(lineage = " + candidates[num].lineage + ")")
                sys.stderr.write("\n")
            sys.stderr.flush()

            if default:
                best = str(default)
            else:
                best = input("Number of the best representative? ")
            while best not in candidates.keys():
                best = input("Invalid number. Number of the best representative? ")
            for num_id in candidates:
                if num_id != best:
                    candidates[num_id].cluster_rep = False

    return


def screen_filter_taxa(fasta_records: dict, screen_strs="", filter_strs="", guarantees=None) -> dict:
    """
    Searches the ReferenceSequence lineages in fasta_records for taxa that are either
    not meant to be retained (controlled by screen_strs) or supposed to be filtered out (controlled by filter_strs).

    :param fasta_records: A dictionary mapping `treesapp_id`s (integers) to ReferenceSequence objects
    :param screen_strs: The comma-separated string of taxa that should be RETAINED
    :param filter_strs: The comma-separated string of taxa that should be REMOVED
    :param guarantees: Optional set of numerical treesapp IDs that should not be removed from fasta_records
    :return: A new fasta_records dictionary with only the target taxa
    """
    fasta_replace_dict = dict()
    if not screen_strs and not filter_strs:
        return fasta_records

    if screen_strs:
        screen_terms = screen_strs.split(',')
    else:
        screen_terms = []
    if filter_strs:
        filter_terms = filter_strs.split(',')
    else:
        filter_terms = []

    if len(set(screen_terms).intersection(filter_terms)) > 0:
        LOGGER.error("Taxon name(s) {} present in both search and filter terms. This is confusing, please fix.\n"
                      "".format(', '.join(set(screen_terms).intersection(filter_terms))))
        sys.exit(13)

    num_filtered = 0
    num_screened = 0
    saved = set()
    for treesapp_id in fasta_records:
        screen_pass = False
        filter_pass = True
        ref_seq = fasta_records[treesapp_id]  # type: entrez_utils.EntrezRecord
        # Screen
        if len(screen_terms) > 0:
            for term in screen_terms:
                # If any term is found in the lineage, it will pass... unless it fails the filter
                if re.search(term, ref_seq.lineage):
                    screen_pass = True
                    break
        else:
            screen_pass = True
        # Filter
        if len(filter_terms) > 0:
            for term in filter_terms:
                if re.search(term, ref_seq.lineage):
                    filter_pass = False

        if filter_pass and screen_pass:
            fasta_replace_dict[treesapp_id] = ref_seq
        elif guarantees and treesapp_id in guarantees:
            saved.add(treesapp_id)
            fasta_replace_dict[treesapp_id] = ref_seq
        else:
            if screen_pass is False:
                num_screened += 1
            if filter_pass is False:
                num_filtered += 1

    LOGGER.debug('\t' + str(num_screened) + " sequences removed after failing screen.\n" +
                  '\t' + str(num_filtered) + " sequences removed after failing filter.\n" +
                  '\t' + str(len(fasta_replace_dict)) + " sequences retained.\n")
    if saved:
        LOGGER.debug('\t' + str(len(saved)) + " guaranteed sequences saved from taxonomic filtering.\n")

    return fasta_replace_dict


def strip_rank_prefix_from_organisms(entrez_record_dict: dict, taxa_trie: TaxonomicHierarchy) -> None:
    """
    Used for removing the rank-prefix (e.g. n__, d__) from EntrezRecord.organism attributes.
    This is purely for aesthetic reasons as the organism names only show up in the visualized phylogeny.

    :param entrez_record_dict: A dictionary of EntrezRecord values indexed by numerical identifiers
    :param taxa_trie: A TaxonomicHierarchy instance with the canonical_prefix attribute of a compiled re object
    :return: None
    """
    for num_id in entrez_record_dict:  # type: entrez_utils.EntrezRecord
        e_record = entrez_record_dict[num_id]
        if taxa_trie.canonical_prefix.search(e_record.organism):
            e_record.organism = taxa_trie.canonical_prefix.sub('', e_record.organism)
    return


def remove_by_truncated_lineages(fasta_records: dict, min_taxonomic_rank: str, taxa_hierarchy: TaxonomicHierarchy,
                                 guarantees=None) -> dict:
    if min_taxonomic_rank == 'r':
        return fasta_records

    num_removed = 0
    fasta_replace_dict = dict()

    rank_name = taxa_hierarchy.rank_prefix_map[min_taxonomic_rank]

    for treesapp_id in fasta_records:
        ref_seq = fasta_records[treesapp_id]
        # Keep all sequences that are guaranteed to be in the final reference package
        if guarantees and treesapp_id in guarantees:
            fasta_replace_dict[treesapp_id] = ref_seq
            continue
        # Check whether the reference sequence is resolved to at least the rank
        if not taxa_hierarchy.resolved_as(ref_seq.lineage, rank_name):
            num_removed += 1
        else:
            fasta_replace_dict[treesapp_id] = ref_seq

    LOGGER.debug('\t' + str(num_removed) + " sequences removed with truncated taxonomic lineages.\n" +
                  '\t' + str(len(fasta_replace_dict) - num_removed) + " sequences retained for building tree.\n")

    return fasta_replace_dict


def order_dict_by_lineage(ref_seqs: dict) -> dict:
    """
    Re-order the fasta_record_objects by their lineages (not phylogenetic, just alphabetical sort)
    Remove the cluster members since they will no longer be used

    :param ref_seqs: A dictionary mapping `treesapp_id`s (integers) to ReferenceSequence objects
    :return: An ordered, filtered version of the input dictionary
    """
    # Create a new dictionary with lineages as keys
    LOGGER.debug("Re-enumerating the reference sequences in taxonomic order... ")
    lineage_dict = dict()
    sorted_lineage_dict = dict()
    for treesapp_id, e_record in ref_seqs.items():  # type: (str, entrez_utils.EntrezRecord)
        # Skip the redundant sequences that are not cluster representatives
        if not e_record.cluster_rep:
            continue
        try:
            lineage_dict[e_record.lineage].append(e_record)
        except KeyError:
            lineage_dict[e_record.lineage] = [e_record]

    # Now re-write the ref_seq_dict, but the numeric keys are now sorted by lineage
    #  AND it doesn't contain redundant fasta objects
    num_key = 1
    for lineage in sorted(lineage_dict.keys(), key=str):
        for ref_seq in lineage_dict[lineage]:  # type: entrez_utils.EntrezRecord
            if ref_seq.cluster_rep:
                # Replace the treesapp_id object
                code = '_'.join(ref_seq.short_id.split('_')[1:])
                ref_seq.short_id = str(num_key) + '_' + code
                sorted_lineage_dict[str(num_key)] = ref_seq
                num_key += 1

    LOGGER.debug("done.\n")
    return sorted_lineage_dict


def summarize_reference_taxa(reference_dict: dict, t_hierarchy: TaxonomicHierarchy, cluster_lca=False):
    """
    Function for enumerating the representation of each taxonomic rank within the finalized reference sequences

    :param reference_dict: A dictionary holding ReferenceSequence objects indexed by their unique numerical identifier
    :param t_hierarchy: A TaxonomicHierarchy instance
    :param cluster_lca: Boolean specifying whether a cluster's LCA should be used for calculation or not
    :return: A formatted, human-readable string stating the number of unique taxa at each rank
    """
    unclassifieds = 0

    taxonomic_summary_string = t_hierarchy.summarize_taxa()

    for num_id in sorted(reference_dict.keys(), key=int):
        if cluster_lca and reference_dict[num_id].cluster_lca:
            lineage = reference_dict[num_id].cluster_lca
        else:
            lineage = reference_dict[num_id].lineage

        if re.search("unclassified", lineage, re.IGNORECASE) or not t_hierarchy.resolved_as(lineage, "species"):
            unclassifieds += 1

    # Report number of "Unclassified" lineages
    taxonomic_summary_string += "Unclassified and incomplete lineages account for " +\
                                str(unclassifieds) + '/' + str(len(reference_dict.keys())) + ' (' +\
                                str(round(float(unclassifieds*100)/len(reference_dict.keys()), 1)) + "%) references.\n"

    return taxonomic_summary_string


def lineages_to_dict(fasta_entrez_records: dict, taxa_lca=False) -> dict:
    """
    Populates the organism, accession ID and lineage information contained in ReferencePackage.lineage_ids

    :param fasta_entrez_records: Dictionary mapping numbers (internal treesapp identifiers) to ReferenceSequence objects
    :param taxa_lca: Flag indicating whether a cluster's lineage is just the representatives or the LCA of all members
    :return: A dictionary mapping TreeSAPP reference node IDs to a string with their organism, accession and lineage
    """
    no_lineage = list()
    ref_lineage_map = {}
    for treesapp_id in sorted(fasta_entrez_records.keys(), key=int):  # type: str
        # Definitely will not uphold phylogenetic relationships but at least sequences
        # will be in the right neighbourhood rather than ordered by their position in the FASTA file
        reference_sequence = fasta_entrez_records[treesapp_id]  # type: entrez_utils.EntrezRecord

        if not reference_sequence.versioned:
            LOGGER.error("'versioned' attribute has not been set for the following EntrezRecord:\n"
                          "{}\n".format(reference_sequence.get_info()))
            raise AssertionError("Unable to create the sequence lineage map.")

        if taxa_lca:
            lineage = reference_sequence.cluster_lca
        else:
            lineage = reference_sequence.lineage
        if not lineage:
            no_lineage.append(reference_sequence.versioned)
            lineage = ''

        ref_lineage_map[treesapp_id] = "{0} | {1}\t{2}".format(reference_sequence.organism,
                                                               reference_sequence.versioned,
                                                               lineage)

    if len(no_lineage) > 0:
        LOGGER.warning("{0} reference sequences did not have a lineage:\n\t{1}\n".format(len(no_lineage),
                                                                                          "\n\t".join(no_lineage)))

    return ref_lineage_map


def remove_outlier_sequences(fasta_record_objects: dict, od_seq_exe: str, mafft_exe: str,
                             output_dir="./outliers", num_threads=2) -> None:
    od_input = output_dir + "od_input.fasta"
    od_output = output_dir + "outliers.fasta"
    outlier_names = list()
    tmp_dict = dict()

    outlier_test_fasta_dict = order_dict_by_lineage(fasta_record_objects)

    LOGGER.info("Detecting outlier reference sequences... ")
    create_new_ref_fasta(od_input, outlier_test_fasta_dict)
    od_input_m = '.'.join(od_input.split('.')[:-1]) + ".mfa"
    # Perform MSA with MAFFT
    run_mafft(mafft_exe, od_input, od_input_m, num_threads)
    # Run OD-seq on MSA to identify outliers
    run_odseq(od_seq_exe, od_input_m, od_output, num_threads)
    # Remove outliers from fasta_record_objects collection
    outlier_seqs = fasta.read_fasta_to_dict(od_output)
    for seq_num_id in fasta_record_objects:
        ref_seq = fasta_record_objects[seq_num_id]
        tmp_dict[ref_seq.short_id] = ref_seq

    for seq_name in outlier_seqs:
        ref_seq = tmp_dict[seq_name]  # type: entrez_utils.EntrezRecord
        ref_seq.cluster_rep = False
        outlier_names.append(ref_seq.accession)

    LOGGER.info("done.\n")
    LOGGER.debug(str(len(outlier_seqs)) + " outlier sequences detected and discarded.\n\t" +
                  "\n\t".join([outseq for outseq in outlier_names]) + "\n")

    return


def guarantee_ref_seqs(cluster_dict: dict, important_seqs: set) -> dict:
    """
    Ensures all "guaranteed sequences" are representative sequences, swapping non-guaranteed sequences for the
    Cluster.representative where necessary. Also makes sure all guaranteed sequences are accounted for.

    :param cluster_dict:
    :param important_seqs:
    :return:
    """
    num_swaps = 0
    important_finds = set()
    nonredundant_guarantee_cluster_dict = dict()  # Will be used to replace cluster_dict
    expanded_cluster_id = 0
    for cluster_id in sorted(cluster_dict, key=int):
        cluster_inst = cluster_dict[cluster_id]  # type: seq_clustering.Cluster
        representative = cluster_inst.representative
        if len(cluster_inst.members) == 0:
            if representative in important_seqs:
                important_finds.add(representative)
            nonredundant_guarantee_cluster_dict[expanded_cluster_id] = cluster_inst
        else:
            contains_important_seq = False
            # The case where a member of a cluster is a guaranteed sequence, but not the representative
            x = 0
            while x < len(cluster_inst.members):
                member = cluster_inst.members[x]
                if member[0] in important_seqs:
                    nonredundant_guarantee_cluster_dict[expanded_cluster_id] = seq_clustering.Cluster(member[0])
                    nonredundant_guarantee_cluster_dict[expanded_cluster_id].members = []
                    nonredundant_guarantee_cluster_dict[expanded_cluster_id].lca = cluster_inst.lca
                    expanded_cluster_id += 1
                    contains_important_seq = True
                    important_finds.add(member[0])
                    cluster_inst.members.pop(x)
                else:
                    x += 1
            if representative in important_seqs:
                # So there is no opportunity for the important representative sequence to be swapped, clear members
                nonredundant_guarantee_cluster_dict[expanded_cluster_id] = cluster_inst
                important_finds.add(representative)
            elif contains_important_seq and representative not in important_seqs:
                num_swaps += 1
            else:
                nonredundant_guarantee_cluster_dict[expanded_cluster_id] = cluster_inst
        expanded_cluster_id += 1

    # Some final accounting - in case the header formats are altered!
    if important_seqs.difference(important_finds):
        LOGGER.error(str(len(important_finds)) + '/' + str(len(important_seqs)) +
                      " guaranteed sequences found in cluster output file. The following are missing:\n" +
                      ", ".join([vis for vis in list(important_seqs.difference(important_finds))]) + "\n")
        sys.exit(7)
    LOGGER.debug(str(num_swaps) + " former representative sequences were succeeded by 'guaranteed-sequences'.\n")

    return nonredundant_guarantee_cluster_dict


def find_cluster_lca(cluster_dict: dict, fasta_record_objects: dict, header_registry: dict) -> None:
    """
    Populates the cluster_lca attribute for Cluster instances by calculating the lowest common ancestor (LCA)
    across all lineages of sequences in a cluster (inferred using a sequence clustering tool such as VSEARCH)

    :param cluster_dict: A dictionary of Cluster instanced indexed by their numerical cluster IDs
    :param fasta_record_objects: A dictionary of EntrezRecord instances indexed by TreeSAPP numerical IDs
    :param header_registry: A dictionary of Header instances indexed by TreeSAPP numerical IDs used for mapping
    EntrezRecord headers to their original
    :return: None
    """
    # Create a temporary dictionary for faster mapping
    formatted_to_num_map = dict()
    for num_id in fasta_record_objects:
        formatted_to_num_map[header_registry[num_id].original] = num_id

    lineages = list()
    for cluster_id in sorted(cluster_dict, key=int):
        cluster_inst = cluster_dict[cluster_id]  # type: seq_clustering.Cluster
        members = [cluster_inst.representative]
        # format of member list is: [header, identity, member_seq_length/representative_seq_length]
        members += [member[0] for member in cluster_inst.members]
        # Create a lineage list for all sequences in the cluster
        for member in members:
            try:
                num_id = formatted_to_num_map[member]
                lineages.append(fasta_record_objects[num_id].lineage)
            except KeyError:
                LOGGER.warning("Unable to map '{}' to a TreeSAPP numeric ID. "
                                "It will not be used in determining the cluster LCA.\n".format(member))

        cleaned_lineages = clean_lineage_list(lineages)
        cluster_inst.lca = megan_lca(cleaned_lineages)

        lineages.clear()
    formatted_to_num_map.clear()
    return


def formulate_train_command(input_seqs: str, ref_pkg, output_dir: str, args, acc_to_lin=None, seqs_to_lin=None) -> list:
    trainer_cmd = ["-i", input_seqs,
                   "-r", ref_pkg.f__pkl,
                   "-o", output_dir,
                   "-m", ref_pkg.molecule,
                   "--num_procs", str(args.num_threads),
                   "--max_examples", str(args.max_examples),
                   "--svm_kernel", args.kernel]
    if args.trim_align:
        trainer_cmd.append("--trim_align")
    if acc_to_lin:
        trainer_cmd += ["-a", acc_to_lin]
    if seqs_to_lin:
        trainer_cmd += ["--seqs2lineage", seqs_to_lin]

    return trainer_cmd
