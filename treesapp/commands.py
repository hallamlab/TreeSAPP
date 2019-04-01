
import logging
import sys
import re
import os
import shutil
from glob import glob
from random import randint
from . import file_parsers
from .fasta import format_read_fasta, trim_multiple_alignment, write_new_fasta, summarize_fasta_sequences, get_headers,\
    read_fasta_to_dict
from .treesapp_args import TreeSAPPArgumentParser, add_classify_arguments, add_create_arguments,\
    add_evaluate_arguments, add_update_arguments, check_parser_arguments, check_classify_arguments
from . import utilities
from . import entrez_utils
from . import entish
from . import lca_calculations
from . import placement_trainer
from .phylo_dist import trim_lineages_to_rank
from .classy import TreeProtein, MarkerBuild, MarkerTest, ReferencePackage, EntrezRecord, TreeSAPP, Assigner,\
    register_headers, get_header_info, prep_logging
from . import create_refpkg
from .assign import abundify_tree_saps, delete_files, validate_inputs, predict_orfs,\
    get_alignment_dims, hmmsearch_orfs, extract_hmm_matches, write_grouped_fastas, create_ref_phy_files,\
    multiple_alignments, get_sequence_counts, filter_multiple_alignments, check_for_removed_sequences,\
    evaluate_trimming_performance, produce_phy_files, parse_raxml_output, filter_placements, align_reads_to_nucs,\
    write_classified_nuc_sequences, summarize_placements_rpkm, run_rpkm, write_tabular_output, produce_itol_inputs,\
    update_func_tree_workflow
from .external_command_interface import launch_write_command
from .jplace_utils import sub_indices_for_seq_names_jplace
from .clade_exclusion_evaluator import load_ref_seqs, pick_taxonomic_representatives, select_rep_seqs,\
    map_seqs_to_lineages, prep_graftm_ref_files, build_graftm_package, get_classification_performance,\
    write_containment_table, map_headers_to_lineage, graftm_classify, validate_ref_package_files, \
    classify_excluded_taxon, exclude_clade_from_ref_files, determine_containment, parse_distances,\
    write_performance_table, summarize_taxonomic_diversity, remove_clade_exclusion_files


def info(args):
    """

    """
    parser = TreeSAPPArgumentParser(description="Return package, executable and refpkg information.")
    args = parser.parse_args(args)
    prep_logging()
    ts = TreeSAPP("info")

    import treesapp
    import Bio
    import numpy
    import scipy
    import ete3
    logging.info("TreeSAPP version " + treesapp.version + ".\n")

    # Write the version of all python deps
    py_deps = {"BioPython": Bio.__version__,
               "ETE3": ete3.__version__,
               "numpy": numpy.__version__,
               "scipy": scipy.__version__}

    logging.info("Python package dependency versions:\n" +
                 "\n\t".join([k + ": " + v for k, v in py_deps.items()]) + "\n")

    # Write the version of executable deps
    ts.furnish_with_arguments(args)
    logging.info(utilities.executable_dependency_versions(ts.executables))

    if args.verbose:
        # TODO: Write relevant reference package information (e.g. codes, gene names, descriptions)
        pass

    return


def create(args):
    # TODO: Record each external software command and version in log
    ##
    # STAGE 0: PARAMETERIZE - retrieve command line arguments, query user about settings if necessary
    ##
    parser = TreeSAPPArgumentParser(description='Create a reference package for TreeSAPP.')
    add_create_arguments(parser)
    args = parser.parse_args(args)

    log_file_name = args.output + os.sep + "create_" + args.code_name + "_TreeSAPP_log.txt"
    prep_logging(log_file_name, args.verbose)
    logging.info("\n##\t\t\tCreating TreeSAPP reference package for '" + args.code_name + "' \t\t\t##\n")
    logging.info("Command used:\n" + ' '.join(sys.argv) + "\n")

    args = utilities.find_executables(args)

    if args.pc:
        create_refpkg.terminal_commands(args.final_output_dir, args.code_name)
        sys.exit(0)

    # Names of files to be created
    tree_output_dir = args.output + args.code_name + "_phy_files" + os.sep
    accession_map_file = args.output + os.sep + "accession_id_lineage_map.tsv"
    hmm_purified_fasta = args.output + args.code_name + "_hmm_purified.fasta"
    filtered_fasta_name = args.output + '.'.join(
        os.path.basename(args.fasta_input).split('.')[:-1]) + "_filtered.fa"
    uclust_prefix = args.output + '.'.join(
        os.path.basename(filtered_fasta_name).split('.')[:-1]) + "_uclust" + args.identity
    unaln_ref_fasta = args.output + args.code_name + "_ref.fa"  # FASTA file of unaligned reference sequences
    phylip_file = args.output + args.code_name + ".phy"  # Used for building the phylogenetic tree with RAxML

    # Gather all the final TreeSAPP reference files
    ref_pkg = ReferencePackage()
    ref_pkg.gather_package_files(args.code_name, args.final_output_dir, "flat")

    # Create a new MarkerBuild instance to hold all relevant information for recording in ref_build_parameters.tsv
    marker_package = MarkerBuild()
    marker_package.pid = args.identity
    marker_package.cog = args.code_name
    marker_package.molecule = args.molecule
    marker_package.kind = args.kind
    marker_package.denominator = "Z1111"

    ##
    # STAGE 2: FILTER - begin filtering sequences by homology and taxonomy
    ##
    if args.domain:
        if os.path.isfile(hmm_purified_fasta):
            logging.info("Using " + hmm_purified_fasta + " from a previous attempt.\n")
        else:
            logging.info("Searching for domain sequences... ")
            hmm_domtbl_files = utilities.hmmsearch_input_references(args, args.fasta_input)
            logging.info("done.\n")
            hmm_matches = file_parsers.parse_domain_tables(args, hmm_domtbl_files)
            # If we're screening a massive fasta file, we don't want to read every sequence - just those with hits
            # TODO: Implement a screening procedure in _fasta_reader._read_format_fasta()
            fasta_dict = format_read_fasta(args.fasta_input, args.molecule, args.output)
            header_registry = register_headers(get_headers(args.fasta_input))
            marker_gene_dict = utilities.extract_hmm_matches(hmm_matches, fasta_dict, header_registry)
            write_new_fasta(marker_gene_dict, hmm_purified_fasta)
            summarize_fasta_sequences(hmm_purified_fasta)
            utilities.hmm_pile(hmm_matches)

        fasta_dict = format_read_fasta(hmm_purified_fasta, args.molecule, args.output, 110, args.min_seq_length)
        header_registry = register_headers(get_headers(hmm_purified_fasta))
        # Point all future operations to the HMM purified FASTA file as the original input
        args.fasta_input = hmm_purified_fasta
    else:
        fasta_dict = format_read_fasta(args.fasta_input, args.molecule, args.output, 110, args.min_seq_length)
        header_registry = register_headers(get_headers(args.fasta_input))

    ##
    # Synchronize records between fasta_dict and header_registry (e.g. short ones may be removed by format_read_fasta())
    ##
    if len(fasta_dict.keys()) != len(header_registry):
        sync_header_registry = dict()
        excluded_headers = list()
        for num_id in header_registry:
            if header_registry[num_id].formatted not in fasta_dict:
                excluded_headers.append(header_registry[num_id].original)
            else:
                sync_header_registry[num_id] = header_registry[num_id]
        logging.debug("The following " + str(len(excluded_headers)) + " sequences have been excluded:\n" +
                      "\n".join(excluded_headers) + "\n")
        header_registry = sync_header_registry

    ##
    # If there are sequences that needs to be guaranteed to be included,
    #  add them now as its easier to work with more sequences than repeat everything
    ##
    if args.guarantee:
        if not os.path.isfile(args.guarantee):
            logging.error("File '" + args.guarantee + "' does not exist!\n")
            sys.exit(13)
        important_seqs = format_read_fasta(args.guarantee, args.molecule, args.output)
        important_headers = register_headers(get_headers(args.guarantee))
        fasta_dict.update(important_seqs)
        acc = max([int(x) for x in header_registry.keys()])
        for num_id in sorted(important_headers, key=int):
            acc += 1
            header_registry[str(acc)] = important_headers[num_id]
    else:
        important_seqs = None

    ##
    # Determine the format of each sequence name (header) and pull important info (e.g. accession, lineage)
    ##
    fasta_record_objects = get_header_info(header_registry, args.code_name)

    ##
    # Read lineages corresponding to accessions for each sequence if available, otherwise download them
    ##
    query_accession_list, num_lineages_provided = entrez_utils.build_entrez_queries(fasta_record_objects)
    if os.path.isfile(accession_map_file):
        logging.info("Reading cached lineages in '" + accession_map_file + "'... ")
        accession_lineage_map = entrez_utils.read_accession_taxa_map(accession_map_file)
        logging.info("done.\n")
    else:
        if args.accession2taxid:
            # Determine find the query accessions that are located in the provided accession2taxid file
            entrez_record_dict = entrez_utils.map_accession2taxid(query_accession_list, args.accession2taxid)
            # Map lineages to taxids for successfully-mapped query sequences
            entrez_utils.fetch_lineages_from_taxids(entrez_record_dict.values())
            # Use the normal querying functions to obtain lineage information for the unmapped queries
            unmapped_queries = entrez_utils.pull_unmapped_entrez_records(entrez_record_dict.values())
            if len(unmapped_queries) > 0:
                # This tends to be a minority so shouldn't be too taxing
                for e_record in entrez_utils.get_multiple_lineages(unmapped_queries, args.molecule):  # type: EntrezRecord
                    try:
                        entrez_record_dict[e_record.accession] = e_record
                    except KeyError:
                        logging.warning(e_record.accession + " not found in original query list.\n")
                        continue
            entrez_records = list(entrez_record_dict.values())
            entrez_record_dict.clear()
            unmapped_queries.clear()
        else:
            entrez_records = entrez_utils.get_multiple_lineages(query_accession_list, args.molecule)
        accession_lineage_map = entrez_utils.entrez_records_to_accession_lineage_map(entrez_records)
        all_accessions = entrez_utils.entrez_records_to_accessions(entrez_records, query_accession_list)

        # Download lineages separately for those accessions that failed
        # Map proper accession to lineage from the tuple keys (accession, accession.version)
        #  in accession_lineage_map returned by entrez_utils.get_multiple_lineages.
        fasta_record_objects, accession_lineage_map = entrez_utils.verify_lineage_information(accession_lineage_map,
                                                                                              all_accessions,
                                                                                              fasta_record_objects,
                                                                                              num_lineages_provided)
        entrez_utils.write_accession_lineage_map(accession_map_file, accession_lineage_map)
    # Add lineage information to the ReferenceSequence() objects in fasta_record_objects if not contained
    create_refpkg.finalize_ref_seq_lineages(fasta_record_objects, accession_lineage_map)

    ##
    # Perform taxonomic lineage-based filtering and screening based on command-line arguments
    ##
    if args.add_lineage:
        if args.screen or args.filter:
            logging.warning("Skipping taxonomic filtering and screening in `--add_lineage` mode.\n")
    else:
        # Remove the sequences failing 'filter' and/or only retain the sequences in 'screen'
        fasta_record_objects = create_refpkg.screen_filter_taxa(args, fasta_record_objects)
        # Remove the sequence records with low resolution lineages, according to args.min_taxonomic_rank
        fasta_record_objects = create_refpkg.remove_by_truncated_lineages(args.min_taxonomic_rank, fasta_record_objects)

    if len(fasta_record_objects.keys()) < 2:
        logging.error(str(len(fasta_record_objects)) + " sequences post-homology + taxonomy filtering\n")
        sys.exit(11)

    fasta_record_objects = create_refpkg.remove_duplicate_records(fasta_record_objects)

    # Add the respective protein or nucleotide sequence string to each ReferenceSequence object
    for num_id in fasta_record_objects:
        refseq_object = fasta_record_objects[num_id]
        treesapp_id = refseq_object.short_id[1:].split('_')[0]
        refseq_object.sequence = fasta_dict[header_registry[treesapp_id].formatted]
    # Write a new FASTA file containing the sequences that passed the homology and taxonomy filters
    filtered_fasta_dict = dict()
    for num_id in fasta_record_objects:
        refseq_object = fasta_record_objects[num_id]
        formatted_header = header_registry[num_id].formatted
        filtered_fasta_dict[formatted_header] = refseq_object.sequence
    write_new_fasta(filtered_fasta_dict, filtered_fasta_name)

    ##
    # Optionally cluster the input sequences using USEARCH at the specified identity
    ##
    if args.cluster:
        utilities.cluster_sequences(args.executables["usearch"], filtered_fasta_name, uclust_prefix, args.identity)
        args.uc = uclust_prefix + ".uc"

    ##
    # Read the uc file if present
    ##
    if args.uc:
        cluster_dict = file_parsers.read_uc(args.uc)

        # Ensure the headers in cluster_dict have been reformatted if UC file was not generated internally
        if not args.cluster:
            members = list()
            for num_id in cluster_dict:
                cluster = cluster_dict[num_id]
                cluster.representative = utilities.reformat_string(cluster.representative)
                for member in cluster.members:
                    header, identity = member
                    members.append([utilities.reformat_string(header), identity])
                cluster.members = members
                members.clear()
        logging.debug("\t" + str(len(cluster_dict.keys())) + " sequence clusters\n")
        ##
        # Calculate LCA of each cluster to represent the taxonomy of the representative sequence
        ##
        # Create a temporary dictionary for faster mapping
        formatted_to_num_map = dict()
        for num_id in fasta_record_objects:
            formatted_to_num_map[header_registry[num_id].formatted] = num_id

        lineages = list()
        for cluster_id in sorted(cluster_dict, key=int):
            members = [cluster_dict[cluster_id].representative]
            # format of member list is: [header, identity, member_seq_length/representative_seq_length]
            members += [member[0] for member in cluster_dict[cluster_id].members]
            # Create a lineage list for all sequences in the cluster
            for member in members:
                try:
                    num_id = formatted_to_num_map[member]
                    lineages.append(fasta_record_objects[num_id].lineage)
                except KeyError:
                    logging.warning("Unable to map " + str(member) + " to a TreeSAPP numeric ID.\n")
            cleaned_lineages = lca_calculations.clean_lineage_list(lineages)
            cluster_dict[cluster_id].lca = lca_calculations.megan_lca(cleaned_lineages)
            # For debugging
            # if len(lineages) != len(cleaned_lineages) and len(lineages) > 1:
            #     print("Before:")
            #     for l in lineages:
            #         print(l)
            #     print("After:")
            #     for l in cleaned_lineages:
            #         print(l)
            #     print("LCA:", cluster_dict[cluster_id].lca)
            lineages.clear()
        formatted_to_num_map.clear()
    else:
        cluster_dict = None

    ##
    # Swap sequences in 'guarantee' for the representatives, creating new clusters
    ##
    if args.guarantee and args.uc:
        # We don't want to make the tree redundant so instead of simply adding the sequences in guarantee,
        #  we will swap them for their respective representative sequences.
        # All important sequences become representative, even if multiple are in the same cluster
        cluster_dict = create_refpkg.guarantee_ref_seqs(cluster_dict, important_seqs)

    # TODO: Taxonomic normalization

    ##
    # Set the cluster-specific values for ReferenceSequence objects
    ##
    if args.uc and not args.headless:
        # Allow user to select the representative sequence based on organism name, sequence length and similarity
        fasta_record_objects = create_refpkg.present_cluster_rep_options(cluster_dict,
                                                                         fasta_record_objects,
                                                                         header_registry,
                                                                         important_seqs)
    elif args.uc and args.headless:
        create_refpkg.finalize_cluster_reps(cluster_dict, fasta_record_objects, header_registry)
    else:
        for num_id in fasta_record_objects:
            fasta_record_objects[num_id].cluster_rep = True
            # fasta_record_objects[num_id].cluster_lca is left empty

    fasta_record_objects = create_refpkg.remove_outlier_sequences(fasta_record_objects,
                                                                  args.executables["OD-seq"], args.executables["mafft"],
                                                                  args.output, args.num_threads)

    ##
    # Re-order the fasta_record_objects by their lineages (not phylogenetic, just alphabetical sort)
    # Remove the cluster members since they will no longer be used
    ##
    fasta_replace_dict = create_refpkg.order_dict_by_lineage(fasta_record_objects)

    # For debugging. This is the finalized set of reference sequences:
    # for num_id in sorted(fasta_replace_dict, key=int):
    #     fasta_replace_dict[num_id].get_info()

    warnings = create_refpkg.write_tax_ids(args, fasta_replace_dict, ref_pkg.lineage_ids)
    if warnings:
        logging.warning(warnings + "\n")

    logging.info("Generated the taxonomic lineage map " + ref_pkg.lineage_ids + "\n")
    taxonomic_summary = create_refpkg.summarize_reference_taxa(fasta_replace_dict, args.taxa_lca)
    logging.info(taxonomic_summary)
    marker_package.lowest_confident_rank = create_refpkg.estimate_taxonomic_redundancy(args, fasta_replace_dict)

    ##
    # Perform multiple sequence alignment
    ##
    if args.multiple_alignment:
        create_refpkg.create_new_ref_fasta(unaln_ref_fasta, fasta_replace_dict, True)
    else:
        create_refpkg.create_new_ref_fasta(unaln_ref_fasta, fasta_replace_dict)

    if args.molecule == 'rrna':
        create_refpkg.generate_cm_data(args, unaln_ref_fasta)
        args.multiple_alignment = True
    elif args.multiple_alignment is False:
        logging.info("Aligning the sequences using MAFFT... ")
        create_refpkg.run_mafft(args.executables["mafft"], unaln_ref_fasta, ref_pkg.msa, args.num_threads)
        logging.info("done.\n")
    else:
        pass
    ref_aligned_fasta_dict = read_fasta_to_dict(ref_pkg.msa)
    marker_package.num_reps = len(ref_aligned_fasta_dict.keys())
    n_rows, n_cols = file_parsers.multiple_alignment_dimensions(seq_dict=ref_aligned_fasta_dict,
                                                                mfa_file=ref_pkg.msa)
    logging.debug("Reference alignment contains " +
                  str(n_rows) + " sequences with " +
                  str(n_cols) + " character positions.\n")

    ##
    # Build the HMM profile from the aligned reference FASTA file
    ##
    if args.molecule == "rrna":
        # A .cm file has already been generated, no need for HMM
        pass
    else:
        logging.info("Building HMM profile... ")
        hmm_build_command = [args.executables["hmmbuild"],
                             args.final_output_dir + args.code_name + ".hmm",
                             ref_pkg.msa]
        stdout, hmmbuild_pro_returncode = launch_write_command(hmm_build_command)
        logging.info("done.\n")
        logging.debug("\n### HMMBUILD ###\n\n" + stdout)

        if hmmbuild_pro_returncode != 0:
            logging.error("hmmbuild did not complete successfully for:\n" +
                          ' '.join(hmm_build_command) + "\n")
            sys.exit(7)
    ##
    # Optionally trim with BMGE and create the Phylip multiple alignment file
    ##
    dict_for_phy = dict()
    if args.trim_align:
        logging.info("Running BMGE... ")
        trimmed_msa_file = trim_multiple_alignment(args.executables["BMGE.jar"], ref_pkg.msa, args.molecule)
        logging.info("done.\n")

        unique_ref_headers = set(
            [re.sub('_' + re.escape(ref_pkg.prefix), '', x) for x in ref_aligned_fasta_dict.keys()])
        msa_dict, failed_trimmed_msa, summary_str = file_parsers.validate_alignment_trimming([trimmed_msa_file], unique_ref_headers)
        logging.debug("Number of sequences discarded: " + summary_str + "\n")
        if trimmed_msa_file not in msa_dict.keys():
            # At least one of the reference sequences were discarded and therefore this package is invalid.
            logging.error("Trimming removed reference sequences. This indicates you have non-homologous sequences.\n" +
                          "Please improve sequence quality-control and/or re-run without the '--trim_align' flag.\n")
            sys.exit(13)
        aligned_fasta_dict = msa_dict[trimmed_msa_file]
        os.remove(trimmed_msa_file)
    else:
        aligned_fasta_dict = ref_aligned_fasta_dict

    for seq_name in aligned_fasta_dict:
        dict_for_phy[seq_name.split('_')[0]] = aligned_fasta_dict[seq_name]
    phy_dict = utilities.reformat_fasta_to_phy(dict_for_phy)
    utilities.write_phy_file(phylip_file, phy_dict)

    ##
    # Build the tree using either RAxML or FastTree
    ##
    marker_package.tree_tool = create_refpkg.construct_tree(args, phylip_file, tree_output_dir)

    if os.path.exists(unaln_ref_fasta):
        os.remove(unaln_ref_fasta)
    if os.path.exists(phylip_file + ".reduced"):
        os.remove(phylip_file + ".reduced")
    if os.path.exists(args.final_output_dir + "fasta_reader_log.txt"):
        os.remove(args.final_output_dir + "fasta_reader_log.txt")

    if args.fast:
        if args.molecule == "prot":
            marker_package.model = "LG"
        else:
            marker_package.model = "GTRGAMMA"
    else:
        entish.annotate_partition_tree(args.code_name,
                                       fasta_replace_dict,
                                       tree_output_dir + os.sep + "RAxML_bipartitions." + args.code_name)
        marker_package.model = create_refpkg.find_model_used(tree_output_dir + os.sep + "RAxML_info." + args.code_name)
    if marker_package.molecule == "prot":
        marker_package.model = "PROTGAMMA" + marker_package.model
    ref_pkg.sub_model = marker_package.model

    # Build the regression model of placement distances to taxonomic ranks
    marker_package.pfit, _, _ = placement_trainer.regress_rank_distance(args, ref_pkg,
                                                                        accession_lineage_map, ref_aligned_fasta_dict)

    ##
    # Finish validating the file and append the reference package build parameters to the master table
    ##
    ref_pkg.validate(marker_package.num_reps)
    param_file = args.treesapp + "data" + os.sep + "tree_data" + os.sep + "ref_build_parameters.tsv"
    create_refpkg.update_build_parameters(param_file, marker_package)

    logging.info("Data for " + args.code_name + " has been generated successfully.\n")
    create_refpkg.terminal_commands(args.final_output_dir, args.code_name)

    return


def evaluate(args):
    """
    Method for running this script:
        Provide it a FASTA file for which it will determine the taxonomic lineage for each sequence
         and run all taxonomic representative sequences with TreeSAPP then analyze via clade exclusion

    :return:
    """
    parser = TreeSAPPArgumentParser(description='Evaluate classification performance using clade-exclusion analysis.')
    add_evaluate_arguments(parser)
    args = parser.parse_args(args)
    sys.stdout.write("\n##\t\t\tBeginning clade exclusion analysis\t\t\t##\n")

    args = check_parser_arguments(args)
    args.targets = ["ALL"]
    exe_dict = utilities.find_executables(args)

    os.makedirs(args.output, exist_ok=True)
    log_file_name = args.output + os.sep + "Clade_exclusion_analyzer_log.txt"
    prep_logging(log_file_name, args.verbose)
    logging.debug("Command used:\n" + ' '.join(sys.argv) + "\n")

    marker_build_dict = file_parsers.parse_ref_build_params(args)
    ref_leaves = file_parsers.tax_ids_file_to_leaves(os.sep.join([args.treesapp, 'data', 'tree_data',
                                                     "tax_ids_" + args.reference_marker + ".txt"]))
    ref_lineages = dict()
    for leaf in ref_leaves:
        ref_lineages[leaf.number] = leaf.lineage
    marker_eval_inst = MarkerTest(args.reference_marker)
    taxa_choices = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]
    args.taxon_rank = args.taxon_rank.split(',')
    for rank in args.taxon_rank:
        if rank not in taxa_choices:
            logging.error(rank + " not an available option for `--taxon_rank`.\n")
            sys.exit(21)
        else:
            marker_eval_inst.ranks.append(rank)

    # Alert the user if the denominator format was (incorrectly) provided to Clade_exclusion_analyzer
    if re.match("^[A-Z][0-9]{4}$", args.reference_marker):
        code_name = args.reference_marker
        try:
            marker = marker_build_dict[code_name].cog
        except KeyError:
            logging.error("Unable to find '" + code_name + "' in ref_build_parameters collection!\n" +
                          "Has it been added to data/tree_data/ref_build_parameters.tsv and cog_list.tsv?\n")
            sys.exit(21)
    elif len(args.reference_marker) <= 6:
        # Assuming the provided args.reference_marker is a short gene name (e.g. mcrA, nifH)
        marker = args.reference_marker
        code_name = ""
        for denominator in marker_build_dict:
            if marker_build_dict[denominator].cog == marker:
                code_name = denominator
                break
        if not code_name:
            logging.error("Unable to identify the gene name from the code name '" + args.reference_marker + "'.")
            sys.exit(21)
    else:
        logging.error("Wrong format for the reference code_name provided: " + args.reference_marker + "\n")
        sys.exit(21)

    ##
    # Define locations of files TreeSAPP outputs
    ##
    inter_class_file = args.output + os.sep + "tmp_clade_exclusion_assignments.tsv"
    accession_map_file = args.output + os.sep + "accession_id_lineage_map.tsv"
    test_rep_taxa_fasta = args.output + os.sep + "representative_taxa_sequences.fasta"
    performance_table = args.output + os.sep + "clade_exclusion_performance.tsv"
    containment_table = args.output + os.sep + "accuracy.tsv"
    # Determine the analysis stage and user's intentions with four booleans
    accessions_downloaded = False  # Whether the accessions have been downloaded from Entrez
    classified = False  # Has TreeSAPP been completed for these sequences?
    treesapp_output_dir = args.output + "TreeSAPP_output" + os.sep
    # Working with a pre-existing TreeSAPP output directory
    if os.path.exists(treesapp_output_dir):
        extant = True
    else:
        extant = False
    if extant:
        if os.path.isfile(accession_map_file):
            accessions_downloaded = True

    if extant:
        logging.debug("The output directory has been detected.\n")

    if not os.path.isdir(treesapp_output_dir):
        os.makedirs(treesapp_output_dir)

    if args.tool in ["diamond", "graftm"]:
        graftm_files = glob(treesapp_output_dir + os.sep + "*" + os.sep + "*_read_tax.tsv")
        if len(graftm_files) == 1:
            classification_table = glob(treesapp_output_dir + os.sep + "*" + os.sep + "*_read_tax.tsv")[0]
        else:
            classification_table = ''
    else:
        classification_table = treesapp_output_dir + os.sep + "final_outputs" + os.sep + "marker_contig_map.tsv"
    if os.path.isfile(classification_table):
        classified = True

    # Load FASTA data
    fasta_dict = format_read_fasta(args.fasta_input, args.molecule, args.output, 110)
    if args.length:
        for seq_id in fasta_dict:
            if len(fasta_dict[seq_id]) < args.length:
                logging.warning(seq_id + " sequence is shorter than " + str(args.length) + "\n")
            else:
                max_stop = len(fasta_dict[seq_id]) - args.length
                random_start = randint(0, max_stop)
                fasta_dict[seq_id] = fasta_dict[seq_id][random_start:random_start + args.length]
    header_registry = register_headers(get_headers(args.fasta_input))
    # Load the query test sequences as ReferenceSequence objects
    complete_ref_sequences = get_header_info(header_registry)
    complete_ref_sequences = load_ref_seqs(fasta_dict, header_registry, complete_ref_sequences)

    logging.debug("\tNumber of input sequences =\t" + str(len(complete_ref_sequences)) + "\n")

    # Checkpoint one: does anything exist?
    # If not, begin by downloading lineage information for each sequence accession
    args.output = treesapp_output_dir
    if not extant and not accessions_downloaded and not classified:
        # User hasn't analyzed anything, sequences will be taxonomically screened then analyzed
        # NOTE: This is only supported for analyzing a single marker gene
        extant = True
        entrez_query_list, num_lineages_provided = entrez_utils.build_entrez_queries(complete_ref_sequences)
        if len(entrez_query_list) > 0:
            entrez_records = entrez_utils.get_multiple_lineages(entrez_query_list, args.molecule)
            accession_lineage_map = entrez_utils.entrez_records_to_accession_lineage_map(entrez_records)
            all_accessions = entrez_utils.entrez_records_to_accessions(entrez_records, entrez_query_list)

        elif num_lineages_provided > 0:
            logging.info("No Entrez queries are necessary - all sequences have provided lineage information.\n")
            accession_lineage_map = dict()
            all_accessions = []
        else:
            logging.error("No accessions were parsed from FASTA records in " + args.fasta_input)
            sys.exit(21)
        # Download lineages separately for those accessions that failed,
        # map proper accession to lineage from the tuple keys (accession, accession.version)
        #  in accession_lineage_map returned by get_multiple_lineages.
        complete_ref_sequences, accession_lineage_map = entrez_utils.verify_lineage_information(accession_lineage_map,
                                                                                                all_accessions,
                                                                                                complete_ref_sequences,
                                                                                                num_lineages_provided)
        create_refpkg.write_accession_lineage_map(accession_map_file, accession_lineage_map)
        complete_ref_sequences = create_refpkg.finalize_ref_seq_lineages(complete_ref_sequences, accession_lineage_map)
        accessions_downloaded = True
    elif extant and accessions_downloaded:
        # File being read should contain accessions mapped to their lineages for all sequences in input FASTA
        accession_lineage_map = entrez_utils.read_accession_taxa_map(accession_map_file)
        complete_ref_sequences = create_refpkg.finalize_ref_seq_lineages(complete_ref_sequences, accession_lineage_map)
    else:
        logging.error("Logic error.")
        sys.exit(21)
    fasta_record_objects = complete_ref_sequences.values()

    logging.info("Selecting representative sequences for each taxon.\n")

    # Filter the sequences from redundant taxonomic lineages, picking up to 5 representative sequences
    representative_seqs, marker_eval_inst.taxa_filter = pick_taxonomic_representatives(fasta_record_objects,
                                                                                       marker_eval_inst.taxa_filter)
    deduplicated_fasta_dict = select_rep_seqs(representative_seqs, fasta_record_objects)
    write_new_fasta(deduplicated_fasta_dict, test_rep_taxa_fasta)
    rep_accession_lineage_map = map_seqs_to_lineages(accession_lineage_map, deduplicated_fasta_dict)

    # Checkpoint three: We have accessions linked to taxa, and sequences to analyze with TreeSAPP, but not classified

    if extant and accessions_downloaded and not classified:
        # Run TreeSAPP against the provided tax_ids file and the unique taxa FASTA file
        if args.length:
            min_seq_length = str(min(args.length - 10, 30))
        else:
            min_seq_length = str(30)

        validate_ref_package_files(args.treesapp, marker, args.output)

        ranks = {"Kingdom": 0, "Phylum": 1, "Class": 2, "Order": 3, "Family": 4, "Genus": 5, "Species": 6}
        for rank in args.taxon_rank:
            leaf_trimmed_taxa_map = trim_lineages_to_rank(ref_lineages, rank)
            unique_ref_lineages = sorted(set(leaf_trimmed_taxa_map.values()))
            unique_query_lineages = sorted(set(trim_lineages_to_rank(rep_accession_lineage_map, rank).values()))
            depth = ranks[rank]
            for lineage in unique_query_lineages:
                taxon = re.sub(r"([ /])", '_', lineage.split("; ")[-1])
                rank_tax = rank[0] + '_' + taxon
                treesapp_output = args.output + os.sep + rank_tax + os.sep

                optimal_lca_taxonomy = "; ".join(lineage.split("; ")[:-1])
                if optimal_lca_taxonomy not in ["; ".join(tl.split("; ")[:-1]) for tl in unique_ref_lineages
                                                if tl != lineage]:
                    logging.debug("Optimal placement target '" + optimal_lca_taxonomy + "' not in pruned tree.\n")
                    continue

                # Select representative sequences belonging to the taxon being tested
                taxon_rep_seqs = select_rep_seqs(representative_seqs, fasta_record_objects, lineage)
                # Decide whether to continue analyzing taxon based on number of query sequences
                if len(taxon_rep_seqs.keys()) == 0:
                    logging.debug("No query sequences for " + lineage + ".\n")
                    continue

                logging.info("Classifications for the " + rank + " '" + taxon + "' put " + treesapp_output + "\n")
                test_obj = marker_eval_inst.new_taxa_test(rank, lineage)
                test_obj.queries = taxon_rep_seqs.keys()
                test_rep_taxa_fasta = args.output + rank_tax + ".fa"

                if args.tool in ["graftm", "diamond"]:
                    classification_table = os.sep.join([treesapp_output, rank_tax, rank_tax + "_read_tax.tsv"])

                    if not os.path.isfile(classification_table):
                        tax_ids_file = os.sep.join([args.output,
                                                    marker + ".gpkg",
                                                    marker + ".gpkg.refpkg",
                                                    marker + "_taxonomy.csv"])
                        # Copy reference files, then exclude all clades belonging to the taxon being tested
                        prep_graftm_ref_files(args.treesapp, args.output, lineage, marker_eval_inst.target_marker,
                                              depth)
                        build_graftm_package(marker_eval_inst.target_marker,
                                             args.output,
                                             mfa_file=args.output + marker + ".mfa",
                                             fa_file=args.output + marker + ".fa",
                                             threads=args.threads)
                        # Write the query sequences
                        write_new_fasta(taxon_rep_seqs, test_rep_taxa_fasta)

                        graftm_classify(test_rep_taxa_fasta, args.output + os.sep + marker + ".gpkg",
                                        treesapp_output, args.threads, args.tool)

                        if not os.path.isfile(classification_table):
                            # The TaxonTest object is maintained for record-keeping (to track # queries & classifieds)
                            logging.warning("GraftM did not generate output for " + lineage + ". Skipping.\n")
                            # shutil.rmtree(treesapp_output)
                            continue

                        shutil.copy(tax_ids_file, treesapp_output + os.sep + rank_tax + os.sep)

                    tax_ids_file = os.sep.join([treesapp_output, rank_tax, marker + "_taxonomy.csv"])
                    test_obj.taxonomic_tree = lca_calculations.grab_graftm_taxa(tax_ids_file)
                    graftm_assignments = file_parsers.read_graftm_classifications(classification_table)
                    test_obj.assignments = {marker: graftm_assignments}
                    test_obj.filter_assignments(marker_eval_inst.target_marker)
                else:
                    tax_ids_file = treesapp_output + "tax_ids_" + marker + ".txt"
                    classification_table = treesapp_output + "final_outputs" + os.sep + "marker_contig_map.tsv"

                    if not os.path.isfile(classification_table) or not os.path.isfile(tax_ids_file):
                        # Copy reference files, then exclude all clades belonging to the taxon being tested
                        prefix = exclude_clade_from_ref_files(args.treesapp, marker, args.output, lineage, depth,
                                                              args.executables, args.fresh, args.molecule)
                        # Write the query sequences
                        write_new_fasta(taxon_rep_seqs, test_rep_taxa_fasta)
                        classify_excluded_taxon(args, prefix, treesapp_output, marker, min_seq_length,
                                                test_rep_taxa_fasta)
                        if not os.path.isfile(classification_table):
                            # The TaxonTest object is maintained for record-keeping (to track # queries & classifieds)
                            logging.warning("TreeSAPP did not generate output for " + lineage + ". Skipping.\n")
                            shutil.rmtree(treesapp_output)
                            continue
                    else:
                        # Valid number of queries and these sequences have already been classified
                        pass

                    test_obj.taxonomic_tree = lca_calculations.all_possible_assignments(tax_ids_file)
                    if os.path.isfile(classification_table):
                        assigned_lines = file_parsers.read_marker_classification_table(classification_table)
                        test_obj.assignments = file_parsers.parse_assignments(assigned_lines)
                        test_obj.filter_assignments(marker_eval_inst.target_marker)
                        test_obj.distances = parse_distances(assigned_lines)
                    else:
                        logging.error("marker_contig_map.tsv is missing from output directory '" +
                                      os.path.basename(classification_table) + "'\n" +
                                      "Please remove this directory and re-run.\n")
                        sys.exit(21)
        classified = True
    remove_clade_exclusion_files(args.output)

    # Checkpoint four: everything has been prepared, only need to parse the classifications and map lineages
    logging.info("Finishing up the mapping of classified, filtered taxonomic sequences.\n")
    for rank in sorted(marker_eval_inst.taxa_tests):
        for test_obj in marker_eval_inst.taxa_tests[rank]:
            if test_obj.assignments:
                marker_assignments = map_headers_to_lineage(test_obj.assignments, fasta_record_objects)
                # Return the number of correct, classified, and total sequences of that taxon at the current rank
                # Identify the excluded rank for each query sequence
                if len(marker_assignments) == 0:
                    logging.debug("No sequences were classified for " + test_obj.taxon + "\n")
                    continue

                for marker in marker_assignments:
                    marker_eval_inst.markers.add(marker)

                rank_assignments = lca_calculations.identify_excluded_clade(marker_assignments,
                                                           test_obj.taxonomic_tree,
                                                           marker_eval_inst.target_marker)
                for a_rank in rank_assignments:
                    if a_rank != rank and len(rank_assignments[a_rank]) > 0:
                        logging.warning(
                            rank + "-level clade excluded but classifications were found to be " + a_rank +
                            "-level.\nAssignments were: " + str(rank_assignments[a_rank]) + "\n")
                        continue
                    if a_rank not in marker_eval_inst.classifications:
                        marker_eval_inst.classifications[a_rank] = list()
                    if len(rank_assignments[a_rank]) > 0:
                        marker_eval_inst.classifications[a_rank] += rank_assignments[a_rank]

    ##
    # On to the standard clade-exclusion analysis...
    ##
    if marker_eval_inst.taxa_filter["Classified"] != marker_eval_inst.taxa_filter["Unique_taxa"]:
        logging.debug("\n\t" + str(marker_eval_inst.taxa_filter["Classified"] -
                                   marker_eval_inst.taxa_filter["Unique_taxa"]) +
                      " duplicate query taxonomies removed.\n")

    if marker_eval_inst.taxa_filter["Unclassified"] > 0:
        logging.debug("\t" + str(marker_eval_inst.taxa_filter["Unclassified"]) +
                      " query sequences with unclassified taxonomies were removed.\n" +
                      "This is not a problem, its just they have 'unclassified' somewhere in their lineages\n" +
                      "(e.g. Unclassified Bacteria) and this is not good for assessing placement accuracy.\n\n")

    if code_name not in marker_eval_inst.markers and marker not in marker_eval_inst.markers:
        logging.error("No sequences were classified as " + marker + ".\n")
        sys.exit(21)

    # For debugging:
    # for rank in marker_eval_inst.ranks:
    #     distal, pendant, tip = marker_eval_inst.summarize_rankwise_distances(rank)

    # Determine the specificity for each rank
    clade_exclusion_strings = get_classification_performance(marker_eval_inst)
    write_performance_table(args, performance_table, clade_exclusion_strings)
    summarize_taxonomic_diversity(marker_eval_inst)
    containment_strings = determine_containment(marker_eval_inst)
    write_containment_table(args, containment_table, containment_strings)

    return


def assign(args):
    # STAGE 1: Prompt the user and prepare files and lists for the pipeline
    parser = TreeSAPPArgumentParser(description='Taxonomically classify sequences through evolutionary placement.')
    add_classify_arguments(parser)
    args = parser.parse_args(args)

    log_file_name = args.output + os.sep + "TreeSAPP_classify_log.txt"
    prep_logging(log_file_name, args.verbose)

    check_parser_arguments(args)
    assign_this = Assigner(args)
    check_classify_arguments(assign_this, args)

    marker_build_dict = file_parsers.parse_ref_build_params(args)
    tree_numbers_translation = file_parsers.read_species_translation_files(args, marker_build_dict)
    if args.check_trees:
        validate_inputs(args, marker_build_dict)
    if args.skip == 'n':
        # STAGE 2: Predict open reading frames (ORFs) if the input is an assembly, read, format and write the FASTA
        if args.molecule == "dna":
            # args.fasta_input is set to the predicted ORF protein sequences
            args = predict_orfs(args)
        logging.info("Formatting " + args.fasta_input + " for pipeline... ")
        formatted_fasta_dict = format_read_fasta(args.fasta_input, "prot", args.output)
        logging.info("done.\n")

        logging.info(
            "\tTreeSAPP will analyze the " + str(len(formatted_fasta_dict)) + " sequences found in input.\n")
        if re.match(r'\A.*\/(.*)', args.fasta_input):
            input_multi_fasta = os.path.basename(args.fasta_input)
        else:
            input_multi_fasta = args.fasta_input
        args.formatted_input_file = args.var_output_dir + input_multi_fasta + "_formatted.fasta"
        logging.debug("Writing formatted FASTA file to " + args.formatted_input_file + "... ")
        formatted_fasta_files = write_new_fasta(formatted_fasta_dict, args.formatted_input_file)
        logging.debug("done.\n")
        ref_alignment_dimensions = get_alignment_dims(args, marker_build_dict)

        # STAGE 3: Run hmmsearch on the query sequences to search for marker homologs
        hmm_domtbl_files = hmmsearch_orfs(args, marker_build_dict)
        hmm_matches = file_parsers.parse_domain_tables(args, hmm_domtbl_files)
        extracted_seq_dict, numeric_contig_index = extract_hmm_matches(hmm_matches, formatted_fasta_dict)
        homolog_seq_files = write_grouped_fastas(extracted_seq_dict, numeric_contig_index,
                                                 marker_build_dict, args.var_output_dir)

        # STAGE 4: Run hmmalign or PaPaRa, and optionally BMGE, to produce the MSAs required to for the ML estimations
        create_ref_phy_files(args, homolog_seq_files, marker_build_dict, ref_alignment_dimensions)
        concatenated_msa_files = multiple_alignments(args, homolog_seq_files, marker_build_dict)
        # TODO: Wrap this into a function
        file_types = set()
        for mc in concatenated_msa_files:
            sample_msa_file = concatenated_msa_files[mc][0]
            f_ext = sample_msa_file.split('.')[-1]
            if re.match("phy|phylip", f_ext):
                file_types.add("Phylip")
            elif re.match("sto|stockholm", f_ext):
                file_types.add("Stockholm")
            elif re.match("mfa|fa|fasta", f_ext):
                file_types.add("Fasta")
            else:
                logging.error("Unrecognized file extension: '" + f_ext + "'")
                sys.exit(3)
        if len(file_types) > 1:
            logging.error(
                "Multiple file types detected in multiple alignment files:\n" + ','.join(file_types) + "\n")
            sys.exit(3)
        elif len(file_types) == 0:
            logging.error("No alignment files were generated!\n")
            sys.exit(3)
        else:
            file_type = file_types.pop()
        alignment_length_dict = get_sequence_counts(concatenated_msa_files, ref_alignment_dimensions,
                                                    args.verbose, file_type)

        if args.trim_align:
            tool = "BMGE"
            trimmed_mfa_files = filter_multiple_alignments(args, concatenated_msa_files, marker_build_dict, tool)
            qc_ma_dict = check_for_removed_sequences(args, trimmed_mfa_files, concatenated_msa_files,
                                                     marker_build_dict)
            evaluate_trimming_performance(qc_ma_dict, alignment_length_dict, concatenated_msa_files, tool)
            phy_files = produce_phy_files(args, qc_ma_dict)
        else:
            phy_files = concatenated_msa_files
        delete_files(args, 3)

        # STAGE 5: Run RAxML to compute the ML estimations
        utilities.launch_evolutionary_placement_queries(args, phy_files, marker_build_dict)
        sub_indices_for_seq_names_jplace(args, numeric_contig_index, marker_build_dict)
    tree_saps, itol_data = parse_raxml_output(args, marker_build_dict)
    tree_saps = filter_placements(args, tree_saps, marker_build_dict)

    abundance_dict = dict()
    if args.molecule == "dna":
        sample_name = '.'.join(os.path.basename(re.sub("_ORFs", '', args.fasta_input)).split('.')[:-1])
        orf_nuc_fasta = args.final_output_dir + sample_name + "_classified_seqs.fna"
        if not os.path.isfile(orf_nuc_fasta):
            logging.info("Creating nucleotide FASTA file of classified sequences '" + orf_nuc_fasta + "'... ")
            genome_nuc_genes_file = args.final_output_dir + sample_name + "_ORFs.fna"
            if os.path.isfile(genome_nuc_genes_file):
                nuc_orfs_formatted_dict = format_read_fasta(genome_nuc_genes_file, 'dna', args.output)
                write_classified_nuc_sequences(tree_saps, nuc_orfs_formatted_dict, orf_nuc_fasta)
                logging.info("done.\n")
            else:
                logging.info("failed.\nWARNING: Unable to read '" + genome_nuc_genes_file + "'.\n" +
                             "Cannot create the nucleotide FASTA file of classified sequences!\n")
        if args.rpkm:
            sam_file = align_reads_to_nucs(args, orf_nuc_fasta)
            rpkm_output_file = run_rpkm(args, sam_file, orf_nuc_fasta)
            abundance_dict = file_parsers.read_rpkm(rpkm_output_file)
            summarize_placements_rpkm(args, abundance_dict, marker_build_dict)
    else:
        for refpkg_code in tree_saps:
            for placed_seq in tree_saps[refpkg_code]:  # type: TreeProtein
                abundance_dict[placed_seq.contig_name + '|' + placed_seq.name] = 1.0

    abundify_tree_saps(tree_saps, abundance_dict)
    write_tabular_output(args, tree_saps, tree_numbers_translation, marker_build_dict)
    produce_itol_inputs(args, tree_saps, marker_build_dict, itol_data)
    delete_files(args, 4)

    # STAGE 6: Optionally update the reference tree
    if args.update_tree:
        for marker_code in args.targets:
            update_func_tree_workflow(args, marker_build_dict[marker_code])

    delete_files(args, 5)

    return


def update(args):
    parser = TreeSAPPArgumentParser(description='Update a TreeSAPP reference package with newly identified sequences.')
    add_update_arguments(parser)
    args = parser.parse_args(args)
    return


def train(args):
    parser = TreeSAPPArgumentParser(description='Model phylogenetic distances across taxonomic ranks.')
    args = parser.parse_args(args)
    return
