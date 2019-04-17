
import logging
import sys
import re
import os
import shutil
from random import randint
from . import file_parsers
from .fasta import format_read_fasta, trim_multiple_alignment, write_new_fasta, get_headers,\
    read_fasta_to_dict, register_headers, write_classified_sequences, FASTA
from .treesapp_args import TreeSAPPArgumentParser, add_classify_arguments, add_create_arguments,\
    add_evaluate_arguments, add_update_arguments, check_parser_arguments, check_evaluate_arguments,\
    check_classify_arguments, check_create_arguments, add_trainer_arguments, check_trainer_arguments, check_updater_arguments
from . import utilities
from . import wrapper
from . import entrez_utils
from . import entish
from . import lca_calculations
from . import placement_trainer
from .phylo_dist import trim_lineages_to_rank
from .classy import TreeProtein, MarkerBuild, TreeSAPP, Assigner, Evaluator, Creator, PhyTrainer, Updater,\
    get_header_info, prep_logging
from . import create_refpkg
from .assign import abundify_tree_saps, delete_files, validate_inputs,\
    get_alignment_dims, extract_hmm_matches, write_grouped_fastas, create_ref_phy_files,\
    multiple_alignments, get_sequence_counts, filter_multiple_alignments, check_for_removed_sequences,\
    evaluate_trimming_performance, produce_phy_files, parse_raxml_output, filter_placements, align_reads_to_nucs,\
    summarize_placements_rpkm, run_rpkm, write_tabular_output, produce_itol_inputs
from .jplace_utils import sub_indices_for_seq_names_jplace
from .clade_exclusion_evaluator import load_ref_seqs, pick_taxonomic_representatives, select_rep_seqs,\
    map_seqs_to_lineages, prep_graftm_ref_files, build_graftm_package, map_headers_to_lineage, graftm_classify,\
    validate_ref_package_files, restore_reference_package, exclude_clade_from_ref_files, determine_containment,\
    parse_distances, remove_clade_exclusion_files, load_rank_depth_map


def info(sys_args):
    """

    """
    parser = TreeSAPPArgumentParser(description="Return package, executable and refpkg information.")
    args = parser.parse_args(sys_args)
    prep_logging()
    ts = TreeSAPP("info")

    import treesapp
    import Bio
    import numpy
    import scipy
    import ete3
    logging.info("TreeSAPP version " + treesapp.version + ".\n")

    # Write the version of all python deps
    py_deps = {"biopython": Bio.__version__,
               "ete3": ete3.__version__,
               "numpy": numpy.__version__,
               "scipy": scipy.__version__}

    logging.info("Python package dependency versions:\n\t" +
                 "\n\t".join([k + ": " + v for k, v in py_deps.items()]) + "\n")

    # Write the version of executable deps
    ts.furnish_with_arguments(args)
    logging.info(utilities.executable_dependency_versions(ts.executables))

    if args.verbose:
        # TODO: Write relevant reference package information (e.g. codes, gene names, descriptions)
        pass

    return


def train(sys_args):
    parser = TreeSAPPArgumentParser(description='Model phylogenetic distances across taxonomic ranks.')
    add_trainer_arguments(parser)
    args = parser.parse_args(sys_args)
    # TODO: Prevent hmmalign_queries_aligned-BMGE.fasta.reduced file from being written to cwd
    log_file_name = args.output + os.sep + "TreeSAPP_trainer_log.txt"
    prep_logging(log_file_name, args.verbose)
    logging.info("\n##\t\t\tTrain taxonomic rank-placement distance model\t\t\t##\n")

    check_parser_arguments(args, sys_args)
    ts_trainer = PhyTrainer()
    ts_trainer.furnish_with_arguments(args)
    marker_build_dict = file_parsers.parse_ref_build_params(ts_trainer.treesapp_dir, [])
    check_trainer_arguments(ts_trainer, args, marker_build_dict)
    ts_trainer.ref_pkg.gather_package_files(args.pkg_path)
    ts_trainer.ref_pkg.validate()

    ref_seqs = FASTA(args.input)

    # Get the model to be used for phylogenetic placement
    ts_trainer.validate_continue(args)
    for denominator in marker_build_dict:
        marker_build = marker_build_dict[denominator]
        if marker_build.cog == ts_trainer.ref_pkg.prefix and args.molecule == marker_build.molecule:
            ts_trainer.ref_pkg.sub_model = marker_build_dict[denominator].model
            break
    if not ts_trainer.ref_pkg.sub_model:
        logging.error("Unable to find the substitution model used for " + ts_trainer.ref_pkg.prefix + ".\n")
        sys.exit(33)

    if ts_trainer.stage_status("search"):
        # Read the FASTA into a dictionary - homologous sequences will be extracted from this
        ref_seqs.fasta_dict = format_read_fasta(ts_trainer.input_sequences, ts_trainer.molecule_type, ts_trainer.output_dir)
        ref_seqs.header_registry = register_headers(get_headers(ts_trainer.input_sequences))

        logging.info("Searching for domain sequences... ")
        hmm_domtbl_files = wrapper.run_hmmsearch(ts_trainer.executables["hmmsearch"],
                                                 ts_trainer.ref_pkg.profile,
                                                 ts_trainer.input_sequences,
                                                 ts_trainer.var_output_dir)
        logging.info("done.\n")
        hmm_matches = file_parsers.parse_domain_tables(args, hmm_domtbl_files)
        marker_gene_dict = utilities.extract_hmm_matches(hmm_matches, ref_seqs.fasta_dict, ref_seqs.header_registry)
        ref_seqs.summarize_fasta_sequences()
        write_new_fasta(marker_gene_dict, ts_trainer.hmm_purified_seqs)
        utilities.hmm_pile(hmm_matches)
    else:
        ref_seqs.load_fasta()
        ref_seqs.change_dict_keys("formatted")
        ts_trainer.hmm_purified_seqs = ts_trainer.input_sequences

    # Get the lineage information for the training/query sequences
    fasta_record_objects = get_header_info(ref_seqs.header_registry)
    fasta_record_objects = load_ref_seqs(ref_seqs.fasta_dict, ref_seqs.header_registry, fasta_record_objects)
    entrez_query_list, num_lineages_provided = entrez_utils.build_entrez_queries(fasta_record_objects)

    if ts_trainer.stage_status("lineages"):
        entrez_records = entrez_utils.map_accessions_to_lineages(entrez_query_list, args.molecule, args.acc_to_taxid)
        accession_lineage_map = entrez_utils.entrez_records_to_accession_lineage_map(entrez_records)
        all_accessions = entrez_utils.entrez_records_to_accessions(entrez_records, entrez_query_list)

        # Download lineages separately for those accessions that failed
        # Map proper accession to lineage from the tuple keys (accession, accession.version)
        #  in accession_lineage_map returned by entrez_utils.get_multiple_lineages.
        fasta_record_objects, accession_lineage_map = entrez_utils.verify_lineage_information(accession_lineage_map,
                                                                                              all_accessions,
                                                                                              fasta_record_objects,
                                                                                              num_lineages_provided)
        entrez_utils.write_accession_lineage_map(ts_trainer.acc_to_lin, accession_lineage_map)
        # Add lineage information to the ReferenceSequence() objects in fasta_record_objects if not contained
    else:
        logging.info("Reading cached lineages in '" + ts_trainer.acc_to_lin + "'... ")
        accession_lineage_map = entrez_utils.read_accession_taxa_map(ts_trainer.acc_to_lin)
        logging.info("done.\n")

    # Read in the reference fasta file
    ref_fasta_dict = read_fasta_to_dict(ts_trainer.ref_pkg.msa)

    taxa_evo_dists = dict()

    # Goal is to use the distances already calculated but re-print
    if os.path.isfile(ts_trainer.placement_summary) and not args.overwrite:
        # Read the summary file and pull the phylogenetic distances for each rank
        taxa_evo_dists = placement_trainer.read_placement_summary(ts_trainer.placement_summary)
        # Remove any ranks that are not to be used in this estimation
        estimated_ranks = set(taxa_evo_dists.keys())
        for rank_key in estimated_ranks.difference(set(ts_trainer.training_ranks.keys())):
            taxa_evo_dists.pop(rank_key)

    if len(set(ts_trainer.training_ranks.keys()).difference(set(taxa_evo_dists.keys()))) > 0:
        pfit_array, taxa_evo_dists, pqueries = placement_trainer.regress_rank_distance(ts_trainer.hmm_purified_seqs,
                                                                                       ts_trainer.executables,
                                                                                       ts_trainer.ref_pkg,
                                                                                       accession_lineage_map,
                                                                                       ref_fasta_dict,
                                                                                       ts_trainer.var_output_dir,
                                                                                       ts_trainer.molecule_type,
                                                                                       ts_trainer.training_ranks,
                                                                                       args.num_threads)
        # Write the tab-delimited file with metadata included for each placement
        placement_trainer.write_placement_table(pqueries, ts_trainer.placement_table, args.name)
    else:
        pfit_array = placement_trainer.complete_regression(taxa_evo_dists, ts_trainer.training_ranks)
        if pfit_array:
            logging.info("Placement distance regression model complete.\n")
        else:
            logging.info("Unable to complete phylogenetic distance and rank correlation.\n")

    # Write the text file containing distances used in the regression analysis
    with open(ts_trainer.placement_summary, 'w') as out_handler:
        trained_string = "Regression parameters = " + re.sub(' ', '', str(pfit_array)) + "\n"
        ranks = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]
        for rank in ranks:
            trained_string += "# " + rank + "\n"
            if rank in taxa_evo_dists:
                trained_string += str(sorted(taxa_evo_dists[rank], key=float)) + "\n"
            trained_string += "\n"
        out_handler.write(trained_string)

    return


def create(sys_args):
    parser = TreeSAPPArgumentParser(description="Create a reference package for TreeSAPP.")
    add_create_arguments(parser)
    args = parser.parse_args(sys_args)

    log_file_name = args.output + os.sep + "TreeSAPP_create_" + args.refpkg_name + "_log.txt"
    prep_logging(log_file_name, args.verbose)
    logging.info("\n##\t\t\tCreating TreeSAPP reference package\t\t\t##\n")
    # TODO: prevent the log file being removed by overwrite

    check_parser_arguments(args, sys_args)
    ts_create = Creator()
    ts_create.furnish_with_arguments(args)
    check_create_arguments(ts_create, args)
    ts_create.validate_continue(args)

    # Gather all the final TreeSAPP reference files
    ts_create.ref_pkg.gather_package_files(ts_create.final_output_dir, ts_create.molecule_type, "flat")

    # Create a new MarkerBuild instance to hold all relevant information for recording in ref_build_parameters.tsv
    # TODO: Merge the MarkerBuild and ReferencePackage classes
    marker_package = MarkerBuild()
    marker_package.pid = args.identity
    marker_package.cog = ts_create.ref_pkg.prefix
    marker_package.molecule = args.molecule
    marker_package.kind = args.kind
    marker_package.denominator = "Z1111"

    ref_seqs = FASTA(args.input)

    if ts_create.stage_status("search"):
        # Read the FASTA into a dictionary - homologous sequences will be extracted from this
        ref_seqs.fasta_dict = format_read_fasta(args.input, ts_create.molecule_type, ts_create.output_dir)
        ref_seqs.header_registry = register_headers(get_headers(args.input))

        logging.info("Searching for domain sequences... ")
        hmm_domtbl_files = wrapper.run_hmmsearch(ts_create.executables["hmmsearch"], args.domain, args.input,
                                                 ts_create.var_output_dir)
        logging.info("done.\n")
        hmm_matches = file_parsers.parse_domain_tables(args, hmm_domtbl_files)
        marker_gene_dict = utilities.extract_hmm_matches(hmm_matches, ref_seqs.fasta_dict, ref_seqs.header_registry)
        ref_seqs.summarize_fasta_sequences()
        write_new_fasta(marker_gene_dict, ts_create.hmm_purified_seqs)
        utilities.hmm_pile(hmm_matches)
    else:
        ts_create.hmm_purified_seqs = ts_create.input_sequences

    ##
    # Synchronize records between fasta_dict and header_registry (e.g. short ones may be removed by format_read_fasta())
    ##
    ref_seqs.file = ts_create.hmm_purified_seqs
    ref_seqs.fasta_dict = format_read_fasta(ref_seqs.file, marker_package.molecule, ts_create.output_dir,
                                            110, args.min_seq_length)
    ref_seqs.header_registry = register_headers(get_headers(ref_seqs.file))
    ref_seqs.synchronize_seqs_n_headers()

    ##
    # If there are sequences that needs to be guaranteed to be included,
    #  add them now as its easier to work with more sequences than repeat everything
    ##
    if args.guarantee:
        ref_seqs.update(args.guarantee)

    fasta_record_objects = get_header_info(ref_seqs.header_registry, ts_create.ref_pkg.prefix)
    fasta_record_objects = load_ref_seqs(ref_seqs.fasta_dict, ref_seqs.header_registry, fasta_record_objects)
    entrez_query_list, num_lineages_provided = entrez_utils.build_entrez_queries(fasta_record_objects)

    logging.debug("\tNumber of input sequences =\t" + str(len(fasta_record_objects)) + "\n")

    if ts_create.stage_status("lineages"):
        entrez_records = entrez_utils.map_accessions_to_lineages(entrez_query_list, args.molecule, args.acc_to_taxid)
        accession_lineage_map = entrez_utils.entrez_records_to_accession_lineage_map(entrez_records)
        all_accessions = entrez_utils.entrez_records_to_accessions(entrez_records, entrez_query_list)

        # Download lineages separately for those accessions that failed
        # Map proper accession to lineage from the tuple keys (accession, accession.version)
        #  in accession_lineage_map returned by entrez_utils.get_multiple_lineages.
        fasta_record_objects, accession_lineage_map = entrez_utils.verify_lineage_information(accession_lineage_map,
                                                                                              all_accessions,
                                                                                              fasta_record_objects,
                                                                                              num_lineages_provided)
        entrez_utils.write_accession_lineage_map(ts_create.acc_to_lin, accession_lineage_map)
        # Add lineage information to the ReferenceSequence() objects in fasta_record_objects if not contained
    else:
        logging.info("Reading cached lineages in '" + ts_create.acc_to_lin + "'... ")
        accession_lineage_map = entrez_utils.read_accession_taxa_map(ts_create.acc_to_lin)
        logging.info("done.\n")
    create_refpkg.finalize_ref_seq_lineages(fasta_record_objects, accession_lineage_map)

    if ts_create.stage_status("clean"):
        # Remove the sequences failing 'filter' and/or only retain the sequences in 'screen'
        fasta_record_objects = create_refpkg.screen_filter_taxa(args, fasta_record_objects)
        # Remove the sequence records with low resolution lineages, according to args.min_taxonomic_rank
        fasta_record_objects = create_refpkg.remove_by_truncated_lineages(args.min_taxonomic_rank, fasta_record_objects)

        fasta_record_objects = create_refpkg.remove_duplicate_records(fasta_record_objects)

        if len(fasta_record_objects.keys()) < 2:
            logging.error(str(len(fasta_record_objects)) + " sequences post-homology + taxonomy filtering\n")
            sys.exit(11)

        # Write a new FASTA file containing the sequences that passed the homology and taxonomy filters
        filtered_fasta_dict = dict()
        for num_id in fasta_record_objects:
            refseq_object = fasta_record_objects[num_id]
            formatted_header = ref_seqs.header_registry[num_id].formatted
            filtered_fasta_dict[formatted_header] = refseq_object.sequence
        write_new_fasta(filtered_fasta_dict, ts_create.filtered_fasta)

    ##
    # Optionally cluster the input sequences using USEARCH at the specified identity
    ##
    if ts_create.stage_status("cluster"):
        if args.cluster:
            wrapper.cluster_sequences(ts_create.executables["usearch"], ts_create.filtered_fasta,
                                      ts_create.uclust_prefix, ts_create.prop_sim)
            ts_create.uc = ts_create.uclust_prefix + ".uc"
        # Read the uc file if present
        if ts_create.uc:
            cluster_dict = file_parsers.read_uc(ts_create.uc)

            # Ensure the headers in cluster_dict have been reformatted if UC file was not generated internally
            if not args.cluster:
                create_refpkg.rename_cluster_headers(cluster_dict)
            logging.debug("\t" + str(len(cluster_dict.keys())) + " sequence clusters\n")
            ##
            # Calculate LCA of each cluster to represent the taxonomy of the representative sequence
            ##
            create_refpkg.cluster_lca(cluster_dict, fasta_record_objects, ref_seqs.header_registry)
        else:
            cluster_dict = None

        ##
        # Swap sequences in 'guarantee' for the representatives, creating new clusters
        ##
        if args.guarantee and ts_create.uc:
            # We don't want to make the tree redundant so instead of simply adding the sequences in guarantee,
            #  we will swap them for their respective representative sequences.
            # All important sequences become representative, even if multiple are in the same cluster
            cluster_dict = create_refpkg.guarantee_ref_seqs(cluster_dict, ref_seqs.amendments)

        ##
        # Set the cluster-specific values for ReferenceSequence objects
        ##
        if ts_create.uc and not args.headless:
            # Allow user to select the representative sequence based on organism name, sequence length and similarity
            fasta_record_objects = create_refpkg.present_cluster_rep_options(cluster_dict,
                                                                             fasta_record_objects,
                                                                             ref_seqs.header_registry,
                                                                             ref_seqs.amendments)
        elif ts_create.uc and args.headless:
            create_refpkg.finalize_cluster_reps(cluster_dict, fasta_record_objects, ref_seqs.header_registry)
        else:
            for num_id in fasta_record_objects:
                fasta_record_objects[num_id].cluster_rep = True
                # fasta_record_objects[num_id].cluster_lca is left empty

    if ts_create.stage_status("build"):
        fasta_record_objects = create_refpkg.remove_outlier_sequences(fasta_record_objects,
                                                                      ts_create.executables["OD-seq"],
                                                                      ts_create.executables["mafft"],
                                                                      ts_create.var_output_dir, args.num_threads)

        ##
        # Re-order the fasta_record_objects by their lineages (not phylogenetic, just alphabetical sort)
        # Remove the cluster members since they will no longer be used
        ##
        fasta_replace_dict = create_refpkg.order_dict_by_lineage(fasta_record_objects)

        # For debugging. This is the finalized set of reference sequences:
        # for num_id in sorted(fasta_replace_dict, key=int):
        #     fasta_replace_dict[num_id].get_info()

        warnings = create_refpkg.write_tax_ids(fasta_replace_dict, ts_create.ref_pkg.lineage_ids, args.taxa_lca)
        if warnings:
            logging.warning(warnings + "\n")

        logging.info("Generated the taxonomic lineage map " + ts_create.ref_pkg.lineage_ids + "\n")
        taxonomic_summary = create_refpkg.summarize_reference_taxa(fasta_replace_dict, args.taxa_lca)
        logging.info(taxonomic_summary)
        marker_package.lowest_confident_rank = create_refpkg.estimate_taxonomic_redundancy(fasta_replace_dict)

        ##
        # Perform multiple sequence alignment
        ##
        if args.multiple_alignment:
            create_refpkg.create_new_ref_fasta(ts_create.unaln_ref_fasta, fasta_replace_dict, True)
        else:
            create_refpkg.create_new_ref_fasta(ts_create.unaln_ref_fasta, fasta_replace_dict)

        if ts_create.molecule_type == 'rrna':
            create_refpkg.generate_cm_data(args, ts_create.unaln_ref_fasta)
            args.multiple_alignment = True
        elif args.multiple_alignment is False:
            logging.info("Aligning the sequences using MAFFT... ")
            create_refpkg.run_mafft(ts_create.executables["mafft"],
                                    ts_create.unaln_ref_fasta, ts_create.ref_pkg.msa, args.num_threads)
            logging.info("done.\n")
        else:
            pass
        ref_aligned_fasta_dict = read_fasta_to_dict(ts_create.ref_pkg.msa)
        marker_package.num_reps = len(ref_aligned_fasta_dict.keys())
        n_rows, n_cols = file_parsers.multiple_alignment_dimensions(seq_dict=ref_aligned_fasta_dict,
                                                                    mfa_file=ts_create.ref_pkg.msa)
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
            wrapper.build_hmm_profile(ts_create.executables["hmmbuild"],
                                      ts_create.ref_pkg.msa,
                                      ts_create.ref_pkg.profile)
        ##
        # Optionally trim with BMGE and create the Phylip multiple alignment file
        ##
        dict_for_phy = dict()
        if args.trim_align:
            logging.info("Running BMGE... ")
            trimmed_msa_file = trim_multiple_alignment(ts_create.executables["BMGE.jar"],
                                                       ts_create.ref_pkg.msa,
                                                       ts_create.molecule_type)
            logging.info("done.\n")

            unique_ref_headers = set(
                [re.sub('_' + re.escape(ts_create.ref_pkg.prefix), '', x) for x in ref_aligned_fasta_dict.keys()])
            msa_dict, failed_trimmed_msa, summary_str = file_parsers.validate_alignment_trimming([trimmed_msa_file],
                                                                                                 unique_ref_headers)
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
        utilities.write_phy_file(ts_create.phylip_file, phy_dict)

        ##
        # Build the tree using either RAxML or FastTree
        ##
        marker_package.tree_tool = wrapper.construct_tree(ts_create.executables, ts_create.molecule_type,
                                                          ts_create.phylip_file, ts_create.phy_dir,
                                                          ts_create.ref_pkg.tree, ts_create.ref_pkg.prefix, args)
        marker_package.model = ts_create.determine_model(args.fast)
        if not args.fast:
            entish.annotate_partition_tree(ts_create.ref_pkg.prefix,
                                           fasta_replace_dict,
                                           ts_create.final_output_dir + os.sep + "RAxML_bipartitions." + ts_create.ref_pkg.prefix)

    if ts_create.stage_status("train"):
        # Build the regression model of placement distances to taxonomic ranks
        trainer_cmd = ["-i", ts_create.input_sequences,
                       "-c", ts_create.ref_pkg.prefix,
                       "-p", ts_create.final_output_dir,
                       "-o", ts_create.var_output_dir + "placement_trainer" + os.sep,
                       "-m", ts_create.molecule_type,
                       "-a", ts_create.acc_to_lin,
                       "-n", str(args.num_threads)]
        if args.trim_align:
            trainer_cmd.append("--trim_align")
        train(trainer_cmd)

    ts_create.remove_intermediates()
    ##
    # Finish validating the file and append the reference package build parameters to the master table
    ##
    ts_create.ref_pkg.validate(marker_package.num_reps)
    param_file = ts_create.treesapp_dir + "data" + os.sep + "ref_build_parameters.tsv"
    create_refpkg.update_build_parameters(param_file, marker_package)

    logging.info("Data for " + ts_create.ref_pkg.prefix + " has been generated successfully.\n")
    if ts_create.stage_status("cc"):
        ts_create.print_terminal_commands()

    return


def update(sys_args):
    parser = TreeSAPPArgumentParser(description='Update a TreeSAPP reference package with newly identified sequences.')
    add_update_arguments(parser)
    args = parser.parse_args(sys_args)

    log_file_name = args.output + os.sep + "TreeSAPP_trainer_log.txt"
    prep_logging(log_file_name, args.verbose)
    logging.info("\n##\t\t\tTrain taxonomic rank-placement distance model\t\t\t##\n")

    check_parser_arguments(args, sys_args)
    ts_updater = Updater()
    ts_updater.furnish_with_arguments(args)
    marker_build_dict = file_parsers.parse_ref_build_params(ts_updater.treesapp_dir, [])
    check_updater_arguments(ts_updater, args, marker_build_dict)
    ts_updater.ref_pkg.gather_package_files()
    ts_updater.ref_pkg.validate()

    ##
    # Pull out sequences from TreeSAPP output
    ##

    ##
    # Add lineages - use taxa if provided with a table mapping contigs to taxa, TreeSAPP-assigned taxonomy otherwise
    ##

    ##
    # Call create to create a new, updated reference package where the new sequences are guaranteed
    ##
    #
    # def update_func_tree_workflow(args, ref_marker: MarkerBuild):
    #
    #     # Load information essential to updating the reference data into a CreateFuncTreeUtility class object
    #     update_tree = CreateFuncTreeUtility(args.output, ref_marker.denominator)
    #
    #     # Get HMM, sequence, reference build, and taxonomic information for the original sequences
    #     ref_hmm_file = args.treesapp + os.sep + 'data' + os.sep + "hmm_data" + os.sep + update_tree.COG + ".hmm"
    #     hmm_length = get_hmm_length(ref_hmm_file)
    #     unaligned_ref_seqs = get_reference_sequence_dict(args, update_tree)
    #     # read_species_translation_files expects the entire marker_build_dict, so we're making a mock one
    #     ref_organism_lineage_info = read_species_translation_files(args, {ref_marker.denominator: ref_marker})
    #
    #     # Set up the output directories
    #     time_of_run = strftime("%d_%b_%Y_%H_%M", gmtime())
    #     project_folder = update_tree.Output + str(time_of_run) + os.sep
    #     raxml_destination_folder = project_folder + "phy_files_%s" % update_tree.COG
    #     final_tree_dir = project_folder + "tree_data" + os.sep
    #     alignment_files_dir = project_folder + "alignment_data" + os.sep
    #     hmm_files_dir = project_folder + "hmm_data" + os.sep
    #     classification_table = update_tree.InputData + os.sep + "final_outputs" + os.sep + "marker_contig_map.tsv"
    #     os.makedirs(project_folder)
    #     os.makedirs(final_tree_dir)
    #     os.makedirs(alignment_files_dir)
    #     os.makedirs(hmm_files_dir)
    #
    #     # Begin finding and filtering the new candidate reference sequences
    #     aa_dictionary = get_new_ref_sequences(args, update_tree)
    #     assignments, n_classified = read_marker_classification_table(classification_table)
    #     if len(aa_dictionary) == 0:
    #         sys.stderr.write("WARNING: No new " + update_tree.COG + " sequences. Skipping update.\n")
    #         return
    #     if n_classified == 0 or update_tree.COG not in assignments.keys():
    #         sys.stderr.write("WARNING: No " + update_tree.COG + " sequences were classified. Skipping update.\n")
    #         return
    #     aa_dictionary = filter_short_sequences(args, aa_dictionary, 0.5 * hmm_length)
    #     if not aa_dictionary:
    #         return
    #     new_ref_seqs_fasta = update_tree.Output + os.path.basename(update_tree.InputData) + \
    #                          "_" + update_tree.COG + "_unaligned.fasta"
    #     # Write only the sequences that have been properly classified
    #     write_new_fasta(aa_dictionary, new_ref_seqs_fasta, None, list(assignments[update_tree.COG].keys()))
    #     # Make sure the tree is updated only if there are novel sequences (i.e. <97% similar to ref sequences)
    #     ref_candidate_alignments = align_ref_queries(args, new_ref_seqs_fasta, update_tree)
    #     # Get the sequences that pass the similarity threshold
    #     new_refs = find_novel_refs(ref_candidate_alignments, aa_dictionary, update_tree)
    #     write_new_fasta(new_refs, new_ref_seqs_fasta)
    #     if args.uclust and len(new_refs.keys()) > 1:
    #         cluster_new_reference_sequences(update_tree, args, new_ref_seqs_fasta)
    #         centroids_fasta = update_tree.Output + "uclust_" + update_tree.COG + ".fasta"
    #     else:
    #         if len(aa_dictionary) == 1 and args.uclust:
    #             sys.stderr.write("WARNING: Not clustering new " + update_tree.COG + " since there is 1 sequence\n")
    #             sys.stderr.flush()
    #         centroids_fasta = new_ref_seqs_fasta
    #
    #     # The candidate set has been finalized. Begin rebuilding!
    #     update_tree.load_new_refs_fasta(args, centroids_fasta, ref_organism_lineage_info)
    #     aligned_fasta = update_tree.align_multiple_sequences(unaligned_ref_seqs, args)
    #     trimal_file = trim_multiple_alignment(args.executables["BMGE.jar"], aligned_fasta, update_tree.marker_molecule)
    #
    #     shutil.move(trimal_file, alignment_files_dir + update_tree.COG + ".fa")
    #     aligned_fasta = alignment_files_dir + update_tree.COG + ".fa"
    #     update_tree.update_tax_ids(args, ref_organism_lineage_info, assignments)
    #
    #     new_hmm_file = update_tree.Output + os.sep + update_tree.COG + ".hmm"
    #     build_hmm(args, alignment_files_dir + update_tree.COG + ".fa", new_hmm_file)
    #     new_hmm_length = get_hmm_length(new_hmm_file)
    #     logging.debug("\tOld HMM length = " + str(hmm_length) + "\n" +
    #                   "\tNew HMM length = " + str(new_hmm_length) + "\n")
    #
    #     os.system('java -cp sub_binaries/readseq.jar run -a -f=12 %s' % aligned_fasta)
    #
    #     phylip_file = update_tree.Output + "%s.phy" % update_tree.COG
    #     os.system('mv %s.phylip %s' % (aligned_fasta, phylip_file))
    #
    #     update_tree.execute_raxml(phylip_file, raxml_destination_folder, args)
    #
    #     # Organize outputs
    #     shutil.move(new_hmm_file, hmm_files_dir)
    #     shutil.move(update_tree.Output + "tax_ids_" + update_tree.COG + ".txt", final_tree_dir)
    #
    #     best_tree = raxml_destination_folder + "/RAxML_bestTree." + update_tree.COG
    #     bootstrap_tree = raxml_destination_folder + "/RAxML_bipartitionsBranchLabels." + update_tree.COG
    #     best_tree_nameswap = final_tree_dir + update_tree.COG + "_tree.txt"
    #     bootstrap_nameswap = final_tree_dir + update_tree.COG + "_bipartitions.txt"
    #     update_tree.swap_tree_names(best_tree, best_tree_nameswap)
    #     update_tree.swap_tree_names(bootstrap_tree, bootstrap_nameswap)
    #     annotate_partition_tree(update_tree.COG,
    #                             update_tree.master_reference_index,
    #                             raxml_destination_folder + os.sep + "RAxML_bipartitions." + update_tree.COG)
    #
    #     prefix = update_tree.Output + update_tree.COG
    #     os.system('mv %s* %s' % (prefix, project_folder))
    #
    #     if args.uclust:
    #         uclust_output_dir = prefix + "_uclust"
    #         os.system('mkdir %s' % uclust_output_dir)
    #
    #         os.system('mv %suclust_* %s' % (update_tree.Output, uclust_output_dir))
    #         os.system('mv %susearch_* %s' % (update_tree.Output, uclust_output_dir))
    #
    #     intermediate_files = [project_folder + update_tree.COG + ".phy",
    #                           project_folder + update_tree.COG + "_gap_removed.fa",
    #                           project_folder + update_tree.COG + "_d_aligned.fasta"]
    #     for useless_file in intermediate_files:
    #         try:
    #             os.remove(useless_file)
    #         except OSError:
    #             sys.stderr.write("WARNING: unable to remove intermediate file " + useless_file + "\n")

    return


def assign(sys_args):
    # STAGE 1: Prompt the user and prepare files and lists for the pipeline
    parser = TreeSAPPArgumentParser(description='Taxonomically classify sequences through evolutionary placement.')
    add_classify_arguments(parser)
    args = parser.parse_args(sys_args)

    log_file_name = args.output + os.sep + "TreeSAPP_classify_log.txt"
    prep_logging(log_file_name, args.verbose)
    logging.info("\n##\t\t\t\tAssigning sequences with TreeSAPP\t\t\t\t##\n\n")

    check_parser_arguments(args, sys_args)
    ts_assign = Assigner()
    ts_assign.furnish_with_arguments(args)
    check_classify_arguments(ts_assign, args)
    ts_assign.validate_continue(args)

    marker_build_dict = file_parsers.parse_ref_build_params(ts_assign.treesapp_dir,
                                                            ts_assign.target_refpkgs)
    ref_alignment_dimensions = get_alignment_dims(ts_assign.treesapp_dir, marker_build_dict)
    tree_numbers_translation = file_parsers.read_species_translation_files(ts_assign.treesapp_dir, marker_build_dict)
    if args.check_trees:
        validate_inputs(args, marker_build_dict)

    ##
    # STAGE 2: Predict open reading frames (ORFs) if the input is an assembly, read, format and write the FASTA
    ##
    if ts_assign.stage_status("orf-call"):
        ts_assign.predict_orfs(args.composition, args.num_threads)
    else:
        ts_assign.orf_file = ts_assign.input_sequences

    if ts_assign.stage_status("clean"):
        logging.info("Formatting " + ts_assign.input_sequences + " for pipeline... ")
        formatted_fasta_dict = format_read_fasta(ts_assign.input_sequences, "prot", ts_assign.output_dir)
        logging.info("done.\n")
        logging.info("\tTreeSAPP will analyze the " + str(len(formatted_fasta_dict)) + " sequences found in input.\n")
        logging.info("Writing formatted FASTA file to " + ts_assign.formatted_input + "... ")
        write_new_fasta(formatted_fasta_dict, ts_assign.formatted_input)
        logging.info("done.\n")

    ##
    # STAGE 3: Run hmmsearch on the query sequences to search for marker homologs
    ##
    if ts_assign.stage_status("search"):
        hmm_domtbl_files = wrapper.hmmsearch_orfs(ts_assign.executables["hmmsearch"], ts_assign.hmm_dir,
                                                  marker_build_dict, ts_assign.formatted_input, ts_assign.var_output_dir,
                                                  args.num_threads)
        hmm_matches = file_parsers.parse_domain_tables(args, hmm_domtbl_files)
        extracted_seq_dict, numeric_contig_index = extract_hmm_matches(hmm_matches, formatted_fasta_dict)
        homolog_seq_files = write_grouped_fastas(extracted_seq_dict, numeric_contig_index,
                                                 marker_build_dict, args.var_output_dir)

    ##
    # STAGE 4: Run hmmalign or PaPaRa, and optionally BMGE, to produce the MSAs required to for the ML estimations
    ##
    if ts_assign.stage_status("align"):
        create_ref_phy_files(ts_assign.aln_dir, ts_assign.var_output_dir,
                             homolog_seq_files, marker_build_dict, ref_alignment_dimensions)
        concatenated_msa_files = multiple_alignments(ts_assign.executables, ts_assign.refpkg_dir, ts_assign.var_output_dir,
                                                     homolog_seq_files, marker_build_dict)
        file_type = utilities.find_msa_type(concatenated_msa_files)
        alignment_length_dict = get_sequence_counts(concatenated_msa_files, ref_alignment_dimensions,
                                                    args.verbose, file_type)

        if args.trim_align:
            tool = "BMGE"
            trimmed_mfa_files = filter_multiple_alignments(ts_assign.executables, concatenated_msa_files,
                                                           marker_build_dict, tool)
            qc_ma_dict = check_for_removed_sequences(ts_assign.aln_dir, trimmed_mfa_files, concatenated_msa_files,
                                                     marker_build_dict, args.min_seq_length)
            evaluate_trimming_performance(qc_ma_dict, alignment_length_dict, concatenated_msa_files, tool)
            phy_files = produce_phy_files(qc_ma_dict)
        else:
            phy_files = concatenated_msa_files
        delete_files(args.delete, ts_assign.var_output_dir, 3)

    ##
    # STAGE 5: Run RAxML to compute the ML estimations
    ##
    if ts_assign.stage_status("place"):
        wrapper.launch_evolutionary_placement_queries(ts_assign.executables, ts_assign.tree_dir,
                                                      phy_files, marker_build_dict,
                                                      ts_assign.var_output_dir, args.num_threads)
        sub_indices_for_seq_names_jplace(ts_assign.var_output_dir, numeric_contig_index, marker_build_dict)

    if ts_assign.stage_status("classify"):
        tree_saps, itol_data = parse_raxml_output(ts_assign.var_output_dir, ts_assign.tree_dir, marker_build_dict)
        tree_saps = filter_placements(tree_saps, marker_build_dict, ts_assign.tree_dir, args.min_likelihood)
        # TODO: Write a FASTA file containing the classified sequences
        # write_classified_sequences(tree_saps, nuc_orfs_formatted_dict, ts_assign.classified_aa_seqs)
        abundance_dict = dict()
        if args.molecule == "dna":
            if not os.path.isfile(ts_assign.classified_nuc_seqs):
                logging.info("Creating nucleotide FASTA file of classified sequences '" +
                             ts_assign.classified_nuc_seqs + "'... ")
                if os.path.isfile(ts_assign.nuc_orfs_file):
                    nuc_orfs_formatted_dict = format_read_fasta(ts_assign.nuc_orfs_file, 'dna', args.output)
                    write_classified_sequences(tree_saps, nuc_orfs_formatted_dict, ts_assign.classified_nuc_seqs)
                    logging.info("done.\n")
                else:
                    logging.info("failed.\nWARNING: Unable to read '" + ts_assign.nuc_orfs_file + "'.\n" +
                                 "Cannot create the nucleotide FASTA file of classified sequences!\n")
            if args.rpkm:
                sam_file = align_reads_to_nucs(args, ts_assign.classified_nuc_seqs)
                rpkm_output_file = run_rpkm(args, sam_file, ts_assign.classified_nuc_seqs)
                abundance_dict = file_parsers.read_rpkm(rpkm_output_file)
                summarize_placements_rpkm(args, abundance_dict, marker_build_dict)
        else:
            for refpkg_code in tree_saps:
                for placed_seq in tree_saps[refpkg_code]:  # type: TreeProtein
                    abundance_dict[placed_seq.contig_name + '|' + placed_seq.name] = 1.0

        abundify_tree_saps(tree_saps, abundance_dict)
        assign_out = ts_assign.final_output_dir + os.sep + "marker_contig_map.tsv"
        write_tabular_output(tree_saps, tree_numbers_translation, marker_build_dict, ts_assign.sample_prefix, assign_out)
        produce_itol_inputs(tree_saps, marker_build_dict, itol_data, ts_assign.output_dir, ts_assign.refpkg_dir)
        delete_files(args.delete, ts_assign.var_output_dir, 4, args.rpkm)

    ##
    # STAGE 6: Optionally update the reference tree
    ##
    if "update" in ts_assign.stages and ts_assign.stage_status("update"):
        for marker_code in args.targets:
            update_args = ["-i", ts_assign.input_sequences,
                           "-c", marker_code,
                           "-t", ts_assign.output_dir,
                           "-o", ts_assign.output_dir + marker_code + "_update" + os.sep,
                           "-m", ts_assign.molecule_type,
                           "-n", str(args.num_threads)]
            if args.cluster:
                update_args.append("--cluster")
            if args.trim_align:
                update_args.append("--trim_align")
            update(update_args)

    delete_files(args.delete, ts_assign.var_output_dir, 5)

    return


def evaluate(sys_args):
    """
    Method for running this script:
        Provide it a FASTA file for which it will determine the taxonomic lineage for each sequence
         and run all taxonomic representative sequences with TreeSAPP then analyze via clade exclusion

    :return:
    """
    parser = TreeSAPPArgumentParser(description='Evaluate classification performance using clade-exclusion analysis.')
    add_evaluate_arguments(parser)
    args = parser.parse_args(sys_args)

    log_file_name = args.output + os.sep + "TreeSAPP_evaluation_log.txt"
    prep_logging(log_file_name, args.verbose)
    logging.info("\n##\t\t\tBeginning clade exclusion analysis\t\t\t##\n")

    check_parser_arguments(args, sys_args)
    ts_evaluate = Evaluator()
    ts_evaluate.furnish_with_arguments(args)

    marker_build_dict = file_parsers.parse_ref_build_params(ts_evaluate.treesapp_dir, ts_evaluate.targets)
    check_evaluate_arguments(ts_evaluate, args, marker_build_dict)
    ts_evaluate.validate_continue(args)
    load_rank_depth_map(ts_evaluate)

    ref_leaves = file_parsers.tax_ids_file_to_leaves(ts_evaluate.tree_dir + "tax_ids_" + args.reference_marker + ".txt")
    ref_lineages = dict()
    for leaf in ref_leaves:
        ref_lineages[leaf.number] = leaf.lineage

    # Load FASTA data
    fasta_dict = format_read_fasta(args.input, args.molecule, ts_evaluate.output_dir, 110)
    if args.length:
        for seq_id in fasta_dict:
            if len(fasta_dict[seq_id]) < args.length:
                logging.warning(seq_id + " sequence is shorter than " + str(args.length) + "\n")
            else:
                max_stop = len(fasta_dict[seq_id]) - args.length
                random_start = randint(0, max_stop)
                fasta_dict[seq_id] = fasta_dict[seq_id][random_start:random_start + args.length]
    header_registry = register_headers(get_headers(args.input))
    # Load the query test sequences as ReferenceSequence objects
    complete_ref_sequences = get_header_info(header_registry)
    complete_ref_sequences = load_ref_seqs(fasta_dict, header_registry, complete_ref_sequences)
    entrez_query_list, num_lineages_provided = entrez_utils.build_entrez_queries(complete_ref_sequences)

    logging.debug("\tNumber of input sequences =\t" + str(len(complete_ref_sequences)) + "\n")

    if ts_evaluate.stage_status("lineages"):
        entrez_records = entrez_utils.map_accessions_to_lineages(entrez_query_list, args.molecule, args.acc_to_taxid)
        accession_lineage_map = entrez_utils.entrez_records_to_accession_lineage_map(entrez_records)
        all_accessions = entrez_utils.entrez_records_to_accessions(entrez_records, entrez_query_list)

        # Download lineages separately for those accessions that failed
        # Map proper accession to lineage from the tuple keys (accession, accession.version)
        #  in accession_lineage_map returned by entrez_utils.get_multiple_lineages.
        fasta_record_objects, accession_lineage_map = entrez_utils.verify_lineage_information(accession_lineage_map,
                                                                                              all_accessions,
                                                                                              complete_ref_sequences,
                                                                                              num_lineages_provided)
        entrez_utils.write_accession_lineage_map(ts_evaluate.acc_to_lin, accession_lineage_map)
        # Add lineage information to the ReferenceSequence() objects in fasta_record_objects if not contained

    else:
        logging.info("Reading cached lineages in '" + ts_evaluate.acc_to_lin + "'... ")
        accession_lineage_map = entrez_utils.read_accession_taxa_map(ts_evaluate.acc_to_lin)
        logging.info("done.\n")
        create_refpkg.finalize_ref_seq_lineages(complete_ref_sequences, accession_lineage_map)

    fasta_record_objects = complete_ref_sequences.values()

    logging.info("Selecting representative sequences for each taxon.\n")

    # Filter the sequences from redundant taxonomic lineages, picking up to 5 representative sequences
    representative_seqs, ts_evaluate.taxa_filter = pick_taxonomic_representatives(fasta_record_objects,
                                                                                  ts_evaluate.taxa_filter)
    deduplicated_fasta_dict = select_rep_seqs(representative_seqs, fasta_record_objects)
    write_new_fasta(deduplicated_fasta_dict, ts_evaluate.test_rep_taxa_fasta)
    rep_accession_lineage_map = map_seqs_to_lineages(accession_lineage_map, deduplicated_fasta_dict)

    # Checkpoint three: We have accessions linked to taxa, and sequences to analyze with TreeSAPP, but not classified
    if ts_evaluate.stage_status("classify"):
        # Run TreeSAPP against the provided tax_ids file and the unique taxa FASTA file
        if args.length:
            min_seq_length = str(min(args.length - 10, 30))
        else:
            min_seq_length = str(30)

        validate_ref_package_files(ts_evaluate.treesapp_dir, ts_evaluate.target_marker.cog, ts_evaluate.output_dir)

        ranks = {"Kingdom": 0, "Phylum": 1, "Class": 2, "Order": 3, "Family": 4, "Genus": 5, "Species": 6}
        for rank in args.taxon_rank:
            leaf_trimmed_taxa_map = trim_lineages_to_rank(ref_lineages, rank)
            unique_ref_lineages = sorted(set(leaf_trimmed_taxa_map.values()))
            unique_query_lineages = sorted(set(trim_lineages_to_rank(rep_accession_lineage_map, rank).values()))
            depth = ranks[rank]
            for lineage in unique_query_lineages:
                taxon = re.sub(r"([ /])", '_', lineage.split("; ")[-1])
                rank_tax = rank[0] + '_' + taxon
                treesapp_output = ts_evaluate.var_output_dir + os.sep + rank_tax + os.sep

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
                test_obj = ts_evaluate.new_taxa_test(rank, lineage)
                test_obj.queries = taxon_rep_seqs.keys()
                test_rep_taxa_fasta = ts_evaluate.var_output_dir + rank_tax + ".fa"

                if args.tool in ["graftm", "diamond"]:
                    classification_table = os.sep.join([treesapp_output, rank_tax, rank_tax + "_read_tax.tsv"])

                    if not os.path.isfile(classification_table):
                        tax_ids_file = os.sep.join([ts_evaluate.var_output_dir,
                                                    ts_evaluate.target_marker.cog + ".gpkg",
                                                    ts_evaluate.target_marker.cog + ".gpkg.refpkg",
                                                    ts_evaluate.target_marker.cog + "_taxonomy.csv"])
                        # Copy reference files, then exclude all clades belonging to the taxon being tested
                        prep_graftm_ref_files(ts_evaluate.treesapp_dir, ts_evaluate.var_output_dir, lineage, ts_evaluate.target_marker,
                                              depth)
                        build_graftm_package(ts_evaluate.target_marker,
                                             ts_evaluate.var_output_dir,
                                             mfa_file=ts_evaluate.var_output_dir + ts_evaluate.target_marker.cog + ".mfa",
                                             fa_file=ts_evaluate.var_output_dir + ts_evaluate.target_marker.cog + ".fa",
                                             threads=args.num_threads)
                        # Write the query sequences
                        write_new_fasta(taxon_rep_seqs, test_rep_taxa_fasta)

                        graftm_classify(test_rep_taxa_fasta, ts_evaluate.var_output_dir + os.sep + ts_evaluate.target_marker.cog + ".gpkg",
                                        treesapp_output, args.num_threads, args.tool)

                        if not os.path.isfile(classification_table):
                            # The TaxonTest object is maintained for record-keeping (to track # queries & classifieds)
                            logging.warning("GraftM did not generate output for " + lineage + ". Skipping.\n")
                            # shutil.rmtree(treesapp_output)
                            continue

                        shutil.copy(tax_ids_file, treesapp_output + os.sep + rank_tax + os.sep)

                    tax_ids_file = os.sep.join([treesapp_output, rank_tax, ts_evaluate.target_marker.cog + "_taxonomy.csv"])
                    test_obj.taxonomic_tree = lca_calculations.grab_graftm_taxa(tax_ids_file)
                    graftm_assignments = file_parsers.read_graftm_classifications(classification_table)
                    test_obj.assignments = {ts_evaluate.target_marker.cog: graftm_assignments}
                    test_obj.filter_assignments(ts_evaluate.target_marker)
                else:
                    tax_ids_file = treesapp_output + "tax_ids_" + ts_evaluate.target_marker.cog + ".txt"
                    classification_table = treesapp_output + "final_outputs" + os.sep + "marker_contig_map.tsv"

                    if not os.path.isfile(classification_table) or not os.path.isfile(tax_ids_file):
                        # Copy reference files, then exclude all clades belonging to the taxon being tested
                        prefix = exclude_clade_from_ref_files(ts_evaluate.treesapp_dir, ts_evaluate.target_marker.cog,
                                                              ts_evaluate.var_output_dir, lineage, depth,
                                                              ts_evaluate.executables, args.fresh, args.molecule)
                        # Write the query sequences
                        write_new_fasta(taxon_rep_seqs, test_rep_taxa_fasta)
                        assign_args = ["-i", test_rep_taxa_fasta, "-o", treesapp_output,
                                       "-m", ts_evaluate.molecule_type, "-n", str(args.num_threads),
                                       "--min_seq_length", str(min_seq_length), "--overwrite", "--delete"]
                        if args.trim_align:
                            assign_args.append("--trim_align")
                        assign(assign_args)
                        restore_reference_package(ts_evaluate.treesapp_dir, prefix,
                                                  treesapp_output, ts_evaluate.target_marker.cog)
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
                        test_obj.filter_assignments(ts_evaluate.target_marker.cog)
                        test_obj.distances = parse_distances(assigned_lines)
                    else:
                        logging.error("marker_contig_map.tsv is missing from output directory '" +
                                      os.path.basename(classification_table) + "'\n" +
                                      "Please remove this directory and re-run.\n")
                        sys.exit(21)
        remove_clade_exclusion_files(ts_evaluate.var_output_dir)

    if ts_evaluate.stage_status("calculate"):
        # everything has been prepared, only need to parse the classifications and map lineages
        logging.info("Finishing up the mapping of classified, filtered taxonomic sequences.\n")
        for rank in sorted(ts_evaluate.taxa_tests):
            for test_obj in ts_evaluate.taxa_tests[rank]:
                if test_obj.assignments:
                    marker_assignments = map_headers_to_lineage(test_obj.assignments, fasta_record_objects)
                    # Return the number of correct, classified, and total sequences of that taxon at the current rank
                    # Identify the excluded rank for each query sequence
                    if len(marker_assignments) == 0:
                        logging.debug("No sequences were classified for " + test_obj.taxon + "\n")
                        continue

                    for marker in marker_assignments:
                        ts_evaluate.markers.add(marker)

                    rank_assignments = lca_calculations.identify_excluded_clade(marker_assignments,
                                                                                test_obj.taxonomic_tree,
                                                                                ts_evaluate.target_marker.cog)
                    for a_rank in rank_assignments:
                        if a_rank != rank and len(rank_assignments[a_rank]) > 0:
                            logging.warning(
                                rank + "-level clade excluded but classifications were found to be " + a_rank +
                                "-level.\nAssignments were: " + str(rank_assignments[a_rank]) + "\n")
                            continue
                        if a_rank not in ts_evaluate.classifications:
                            ts_evaluate.classifications[a_rank] = list()
                        if len(rank_assignments[a_rank]) > 0:
                            ts_evaluate.classifications[a_rank] += rank_assignments[a_rank]

        ##
        # On to the standard clade-exclusion analysis...
        ##
        if ts_evaluate.taxa_filter["Classified"] != ts_evaluate.taxa_filter["Unique_taxa"]:
            logging.debug("\n\t" + str(ts_evaluate.taxa_filter["Classified"] -
                                       ts_evaluate.taxa_filter["Unique_taxa"]) +
                          " duplicate query taxonomies removed.\n")

        if ts_evaluate.taxa_filter["Unclassified"] > 0:
            logging.debug("\t" + str(ts_evaluate.taxa_filter["Unclassified"]) +
                          " query sequences with unclassified taxonomies were removed.\n" +
                          "This is not a problem, its just they have 'unclassified' somewhere in their lineages\n" +
                          "(e.g. Unclassified Bacteria) and this is not good for assessing placement accuracy.\n\n")

        if ts_evaluate.target_marker.denominator not in ts_evaluate.markers and ts_evaluate.target_marker.cog not in ts_evaluate.markers:
            logging.error("No sequences were classified as " + ts_evaluate.target_marker.cog + ".\n")
            sys.exit(21)

        # For debugging:
        # for rank in ts_evaluate.ranks:
        #     distal, pendant, tip = ts_evaluate.summarize_rankwise_distances(rank)

        # Determine the specificity for each rank
        clade_exclusion_strings = ts_evaluate.get_classification_performance()
        ts_evaluate.write_performance_table(clade_exclusion_strings, args.tool)
        ts_evaluate.summarize_taxonomic_diversity()
        containment_strings = determine_containment(ts_evaluate)
        ts_evaluate.write_containment_table(containment_strings, args.tool)

    return
