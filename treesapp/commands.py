
import logging
import sys
import re
import os
import shutil
from random import randint

from samsum.commands import ref_sequence_abundances

from . import entrez_utils
from . import file_parsers
from . import fasta
from .treesapp_args import TreeSAPPArgumentParser, add_classify_arguments, add_create_arguments, add_layer_arguments,\
    add_evaluate_arguments, add_update_arguments, check_parser_arguments, check_evaluate_arguments,\
    check_classify_arguments, check_create_arguments, add_trainer_arguments, check_trainer_arguments,\
    check_updater_arguments, check_purity_arguments, add_purity_arguments, add_abundance_arguments
from . import utilities
from . import wrapper
from . import entish
from . import lca_calculations
from . import placement_trainer
from . import update_refpkg
from . import annotate_extra
from .phylo_dist import trim_lineages_to_rank
from .classy import TreeProtein, MarkerBuild, TreeSAPP, Assigner, Evaluator, Creator, PhyTrainer, Updater, Layerer,\
    prep_logging, dedup_records, TaxonTest, Purity, Abundance
from . import create_refpkg
from .assign import abundify_tree_saps, delete_files, validate_inputs,\
    get_alignment_dims, extract_hmm_matches, write_grouped_fastas, create_ref_phy_files,\
    multiple_alignments, get_sequence_counts, check_for_removed_sequences,\
    evaluate_trimming_performance, produce_phy_files, parse_raxml_output, filter_placements, align_reads_to_nucs,\
    summarize_placements_rpkm, write_tabular_output, produce_itol_inputs, replace_contig_names
from .jplace_utils import sub_indices_for_seq_names_jplace, jplace_parser, demultiplex_pqueries
from .clade_exclusion_evaluator import pick_taxonomic_representatives, select_rep_seqs,\
    map_seqs_to_lineages, prep_graftm_ref_files, build_graftm_package, map_headers_to_lineage, graftm_classify,\
    validate_ref_package_files, restore_reference_package, exclude_clade_from_ref_files, determine_containment,\
    parse_distances, remove_clade_exclusion_files, load_rank_depth_map


def info(sys_args):
    """

    """
    parser = TreeSAPPArgumentParser(description="Return package, executable and refpkg information.")
    args = parser.parse_args(sys_args)
    prep_logging()
    ts_info = TreeSAPP("info")

    import treesapp
    import Bio
    import numpy
    import scipy
    import ete3
    import samsum
    logging.info("TreeSAPP version " + treesapp.__version__ + ".\n")

    # Write the version of all python deps
    py_deps = {"biopython": Bio.__version__,
               "ete3": ete3.__version__,
               "numpy": numpy.__version__,
               "scipy": scipy.__version__,
               "samsum": samsum.__version__}

    logging.info("Python package dependency versions:\n\t" +
                 "\n\t".join([k + ": " + v for k, v in py_deps.items()]) + "\n")

    # Write the version of executable deps
    ts_info.furnish_with_arguments(args)
    logging.info(utilities.executable_dependency_versions(ts_info.executables))

    if args.verbose:
        marker_build_dict = file_parsers.parse_ref_build_params(ts_info.treesapp_dir, [])
        refpkg_summary_str = "\t".join(["Name", "Code-name", "Molecule", "RefPkg-type", "Description", "Last-updated"])
        refpkg_summary_str += "\n"
        for refpkg_code in sorted(marker_build_dict, key= lambda x: marker_build_dict[x].cog):
            refpkg = marker_build_dict[refpkg_code]  # type: MarkerBuild
            refpkg_summary_str += refpkg_code + " -> " + ", ".join(
                [refpkg.cog, refpkg.molecule, refpkg.kind, refpkg.description, refpkg.update]
            ) + "\n"
        logging.info(refpkg_summary_str)

    return


def train(sys_args):
    parser = TreeSAPPArgumentParser(description='Model phylogenetic distances across taxonomic ranks.')
    add_trainer_arguments(parser)
    args = parser.parse_args(sys_args)

    ts_trainer = PhyTrainer()
    ts_trainer.furnish_with_arguments(args)
    ts_trainer.check_previous_output(args)

    # TODO: Prevent hmmalign_queries_aligned-BMGE.fasta.reduced file from being written to cwd
    log_file_name = args.output + os.sep + "TreeSAPP_trainer_log.txt"
    prep_logging(log_file_name, args.verbose)
    logging.info("\n##\t\t\tTrain taxonomic rank-placement distance model\t\t\t##\n")

    check_parser_arguments(args, sys_args)
    marker_build_dict = file_parsers.parse_ref_build_params(ts_trainer.treesapp_dir, [])
    check_trainer_arguments(ts_trainer, args, marker_build_dict)
    ts_trainer.ref_pkg.gather_package_files(args.pkg_path)
    ts_trainer.ref_pkg.validate()

    ref_seqs = fasta.FASTA(ts_trainer.input_sequences)

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
        ref_seqs.fasta_dict = fasta.format_read_fasta(ts_trainer.input_sequences, ts_trainer.molecule_type)
        ref_seqs.header_registry = fasta.register_headers(fasta.get_headers(ts_trainer.input_sequences))

        logging.info("Searching for domain sequences... ")
        hmm_domtbl_files = wrapper.run_hmmsearch(ts_trainer.executables["hmmsearch"],
                                                 ts_trainer.ref_pkg.profile,
                                                 ts_trainer.input_sequences,
                                                 ts_trainer.var_output_dir)
        logging.info("done.\n")
        hmm_matches = file_parsers.parse_domain_tables(args, hmm_domtbl_files)
        marker_gene_dict = utilities.extract_hmm_matches(hmm_matches, ref_seqs.fasta_dict, ref_seqs.header_registry)
        logging.info(ref_seqs.summarize_fasta_sequences())
        fasta.write_new_fasta(marker_gene_dict, ts_trainer.hmm_purified_seqs)
        utilities.hmm_pile(hmm_matches)
    else:
        ref_seqs.load_fasta()
        ref_seqs.change_dict_keys("formatted")
        ts_trainer.hmm_purified_seqs = ts_trainer.input_sequences

    ts_trainer.fetch_entrez_lineages(ref_seqs, args.molecule, args.acc_to_taxid)

    # Read in the reference fasta file
    ref_fasta_dict = fasta.read_fasta_to_dict(ts_trainer.ref_pkg.msa)

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
                                                                                       ts_trainer.seq_lineage_map,
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
            logging.info("Placement distance model complete.\n")
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

    ts_create = Creator()
    ts_create.furnish_with_arguments(args)
    ts_create.check_previous_output(args)

    log_file_name = args.output + os.sep + "TreeSAPP_create_" + args.refpkg_name + "_log.txt"
    prep_logging(log_file_name, args.verbose)
    logging.info("\n##\t\t\tCreating TreeSAPP reference package\t\t\t##\n")

    check_parser_arguments(args, sys_args)
    check_create_arguments(ts_create, args)
    ts_create.validate_continue(args)

    # Gather all the final TreeSAPP reference files
    ts_create.ref_pkg.gather_package_files(ts_create.final_output_dir, ts_create.molecule_type, "flat")

    # Create a new MarkerBuild instance to hold all relevant information for recording in ref_build_parameters.tsv
    # TODO: Merge the MarkerBuild and ReferencePackage classes
    marker_package = MarkerBuild()
    marker_package.pid = ts_create.prop_sim
    marker_package.cog = ts_create.ref_pkg.prefix
    marker_package.molecule = args.molecule
    marker_package.kind = args.kind
    marker_package.denominator = "Z1111"

    ref_seqs = fasta.FASTA(args.input)

    if ts_create.stage_status("search"):
        # Read the FASTA into a dictionary - homologous sequences will be extracted from this
        ref_seqs.fasta_dict = fasta.format_read_fasta(args.input, ts_create.molecule_type)
        ref_seqs.header_registry = fasta.register_headers(fasta.get_headers(args.input))
        logging.debug("Raw, unfiltered sequence summary:\n" + ref_seqs.summarize_fasta_sequences())

        logging.info("Searching for domain sequences... ")
        hmm_domtbl_files = wrapper.run_hmmsearch(ts_create.executables["hmmsearch"], ts_create.hmm_profile,
                                                 args.input, ts_create.var_output_dir)
        logging.info("done.\n")
        hmm_matches = file_parsers.parse_domain_tables(args, hmm_domtbl_files)
        marker_gene_dict = utilities.extract_hmm_matches(hmm_matches, ref_seqs.fasta_dict, ref_seqs.header_registry)
        fasta.write_new_fasta(marker_gene_dict, ts_create.hmm_purified_seqs)
        utilities.hmm_pile(hmm_matches)
    else:
        ts_create.hmm_purified_seqs = ts_create.input_sequences

    ##
    # Synchronize records between fasta_dict and header_registry (e.g. short ones may be removed by format_read_fasta())
    ##
    ref_seqs.file = ts_create.hmm_purified_seqs
    ref_seqs.fasta_dict = fasta.format_read_fasta(fasta_input=ref_seqs.file, molecule=marker_package.molecule,
                                                  min_seq_length=args.min_seq_length)
    ref_seqs.header_registry = fasta.register_headers(fasta.get_headers(ref_seqs.file))
    ref_seqs.synchronize_seqs_n_headers()
    logging.info("Sequence summary:\n" + ref_seqs.summarize_fasta_sequences())

    ##
    # If there are sequences that needs to be guaranteed to be included,
    #  add them now as its easier to work with more sequences than repeat everything
    ##
    if args.guarantee:
        ref_seqs.update(args.guarantee)
        ref_seqs.change_dict_keys("formatted")

    ##
    # Save all sequence names in the header registry as EntrezRecord instances
    # Using the accession-lineage-map (if available) map the sequence names to their respective lineages
    # Proceed with creating the Entrez-queries for sequences lacking lineage information
    ##
    fasta_records = ts_create.fetch_entrez_lineages(ref_seqs, args.molecule, args.acc_to_taxid)
    entrez_utils.fill_ref_seq_lineages(fasta_records, ts_create.seq_lineage_map)

    if ts_create.stage_status("clean"):
        # Remove the sequences failing 'filter' and/or only retain the sequences in 'screen'
        fasta_records = create_refpkg.screen_filter_taxa(fasta_records, args.screen, args.filter, ref_seqs.amendments)
        # Remove the sequence records with low resolution lineages, according to args.min_taxonomic_rank
        fasta_records = create_refpkg.remove_by_truncated_lineages(fasta_records, args.min_taxonomic_rank, ref_seqs.amendments)
        # Ensure there are no records with redundant headers and sequences
        fasta_records = dedup_records(ref_seqs, fasta_records)

        if len(fasta_records.keys()) < 2:
            logging.error(str(len(fasta_records)) + " sequences post-homology + taxonomy filtering\n")
            sys.exit(11)
        # Write a new FASTA file containing the sequences that passed the homology and taxonomy filters

    ref_seqs.file = ts_create.filtered_fasta
    # NOTE: original header must be used as this is being passed to train
    ref_seqs.change_dict_keys("original")
    filtered_headers = [ref_seqs.header_registry[num_id].original for num_id in fasta_records]
    # ref_seqs.keep_only(filtered_headers)  # Currently avoiding this as it causes a KeyError for guaranteed seqs
    fasta.write_new_fasta(fasta_dict=ref_seqs.fasta_dict, fasta_name=ref_seqs.file, headers=filtered_headers)

    ##
    # Optionally cluster the input sequences using USEARCH at the specified identity
    ##
    if ts_create.stage_status("cluster"):
        ref_seqs.change_dict_keys("num")
        # Write a FASTA for clustering containing the formatted headers since
        # not all clustering tools + versions keep whole header - spaces are replaced with underscores
        fasta.write_new_fasta(fasta_dict=ref_seqs.fasta_dict,
                              fasta_name=ts_create.cluster_input,
                              headers=list(fasta_records.keys()))
        if args.cluster:
            wrapper.cluster_sequences(ts_create.executables["usearch"], ts_create.cluster_input,
                                      ts_create.uclust_prefix, ts_create.prop_sim)
            ts_create.uc = ts_create.uclust_prefix + ".uc"
        # Read the uc file if present
        if ts_create.uc:
            cluster_dict = file_parsers.read_uc(ts_create.uc)

            # Revert headers in cluster_dict from 'formatted' back to 'original'
            fasta.rename_cluster_headers(cluster_dict, ref_seqs.header_registry)
            logging.debug("\t" + str(len(cluster_dict.keys())) + " sequence clusters\n")
            ##
            # Calculate LCA of each cluster to represent the taxonomy of the representative sequence
            ##
            create_refpkg.cluster_lca(cluster_dict, fasta_records, ref_seqs.header_registry)
        else:
            cluster_dict = None

        ##
        # Swap sequences in 'guarantee' for the representatives, creating new clusters
        ##
        if args.guarantee and ts_create.uc:
            # We don't want to make the tree redundant so instead of simply adding the sequences in guarantee,
            #  we will swap them for their respective representative sequences.
            # All important sequences become representative, even if multiple are in the same cluster
            very_important_seqs = set([ref_seqs.header_registry[num].original for num in ref_seqs.amendments])
            cluster_dict = create_refpkg.guarantee_ref_seqs(cluster_dict, very_important_seqs)

        ##
        # Set the cluster-specific values for ReferenceSequence objects
        ##
        if ts_create.uc and not args.headless:
            # Allow user to select the representative sequence based on organism name, sequence length and similarity
            fasta_records = create_refpkg.present_cluster_rep_options(cluster_dict, fasta_records,
                                                                      ref_seqs.header_registry, ref_seqs.amendments)
        elif ts_create.uc and args.headless:
            fasta_records = create_refpkg.finalize_cluster_reps(cluster_dict, fasta_records, ref_seqs.header_registry)
        else:
            for num_id in fasta_records:
                fasta_records[num_id].cluster_rep = True
                # fasta_records[num_id].cluster_lca is left empty

    if ts_create.stage_status("build"):
        # TODO: Have a command-line flag to toggle this on (DEFAULT) and off
        fasta_records = create_refpkg.remove_outlier_sequences(fasta_records,
                                                               ts_create.executables["OD-seq"],
                                                               ts_create.executables["mafft"],
                                                               ts_create.var_output_dir, args.num_threads)

        # This precautionary measure is for `create` called from `update` and reference seqs have the assign signature
        accession_ids = [fasta_records[num_id].accession for num_id in fasta_records]
        name_map = update_refpkg.strip_assigment_pattern(accession_ids, ts_create.ref_pkg.prefix)
        for num_id in fasta_records:
            record = fasta_records[num_id]
            record.accession = name_map[record.accession]
        ##
        # Re-order the fasta_records by their lineages (not phylogenetic, just alphabetical sort)
        # Remove the cluster members since they will no longer be used
        ##
        fasta_replace_dict = create_refpkg.order_dict_by_lineage(fasta_records)

        # For debugging. This is the finalized set of reference sequences:
        # for num_id in sorted(fasta_replace_dict, key=int):
        #     fasta_replace_dict[num_id].get_info()

        warnings = create_refpkg.write_tax_ids(fasta_replace_dict, ts_create.ref_pkg.lineage_ids, args.taxa_lca)
        if warnings:
            logging.warning(warnings + "\n")

        logging.info("Generated the taxonomic lineage map " + ts_create.ref_pkg.lineage_ids + "\n")
        taxonomic_summary = create_refpkg.summarize_reference_taxa(fasta_replace_dict, args.taxa_lca)
        logging.info(taxonomic_summary)

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
        ref_aligned_fasta_dict = fasta.read_fasta_to_dict(ts_create.ref_pkg.msa)
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
            trimmed_mfa_files = wrapper.filter_multiple_alignments(ts_create.executables,
                                                                   {marker_package.denominator:
                                                                    [ts_create.ref_pkg.msa]},
                                                                   {marker_package.denominator:
                                                                    marker_package})
            trimmed_mfa_file = trimmed_mfa_files[marker_package.denominator]
            unique_ref_headers = set(
                [re.sub(r"_{0}".format(ts_create.ref_pkg.prefix), '', x) for x in ref_aligned_fasta_dict.keys()])
            qc_ma_dict, failed_trimmed_msa, summary_str = file_parsers.validate_alignment_trimming(trimmed_mfa_file,
                                                                                                   unique_ref_headers)
            logging.debug("Number of sequences discarded: " + summary_str + "\n")
            if len(qc_ma_dict.keys()) == 0:
                # At least one of the reference sequences were discarded and therefore this package is invalid.
                logging.error("Trimming removed reference sequences. This indicates non-homologous sequences.\n" +
                              "Please improve sequence quality-control and/or re-run without the '--trim_align' flag.\n")
                sys.exit(13)
            elif len(qc_ma_dict.keys()) > 1:
                logging.error("Multiple trimmed alignment files are found when only one is expected:\n" +
                              "\n".join([str(k) + ": " + str(qc_ma_dict[k]) for k in qc_ma_dict]))
                sys.exit(13)
            # There is only a single trimmed-MSA file in the dictionary
            for trimmed_msa_file in qc_ma_dict:
                dict_for_phy = qc_ma_dict[trimmed_msa_file]
                os.remove(trimmed_msa_file)
        else:
            for seq_name in ref_aligned_fasta_dict:
                dict_for_phy[seq_name.split('_')[0]] = ref_aligned_fasta_dict[seq_name]

        phy_dict = utilities.reformat_fasta_to_phy(dict_for_phy)
        utilities.write_phy_file(ts_create.phylip_file, phy_dict)

        ##
        # Build the tree using either RAxML or FastTree
        ##
        wrapper.construct_tree(ts_create.executables, ts_create.molecule_type, ts_create.phylip_file,
                               ts_create.phy_dir, ts_create.ref_pkg.tree, ts_create.ref_pkg.prefix, args)
        if not args.fast:
            raw_newick_tree = ts_create.phy_dir + "RAxML_bestTree." + ts_create.ref_pkg.prefix
            bootstrap_tree = ts_create.phy_dir + "RAxML_bipartitionsBranchLabels." + ts_create.ref_pkg.prefix
            entish.annotate_partition_tree(ts_create.ref_pkg.prefix, fasta_replace_dict, bootstrap_tree)
            bootstrap_nameswap = ts_create.final_output_dir + ts_create.ref_pkg.prefix + "_bipartitions.txt"
            utilities.swap_tree_names(raw_newick_tree, ts_create.ref_pkg.tree, ts_create.ref_pkg.prefix)
            utilities.swap_tree_names(bootstrap_tree, bootstrap_nameswap, ts_create.ref_pkg.prefix)

    if args.raxml_model:
        marker_package.model = args.raxml_model
    else:
        marker_package.model = ts_create.determine_model(args.fast)
    param_file = ts_create.treesapp_dir + "data" + os.sep + "ref_build_parameters.tsv"
    refpkg_lineages = [ref.lineage for ref in ts_create.ref_pkg.tax_ids_file_to_leaves()]
    marker_package.lowest_confident_rank = create_refpkg.estimate_taxonomic_redundancy(refpkg_lineages)
    create_refpkg.update_build_parameters(param_file, marker_package)

    # Build the regression model of placement distances to taxonomic ranks
    trainer_cmd = ["-i", ts_create.filtered_fasta,
                   "-c", ts_create.ref_pkg.prefix,
                   "-p", ts_create.final_output_dir,
                   "-o", ts_create.var_output_dir + "placement_trainer" + os.sep,
                   "-m", ts_create.molecule_type,
                   "-a", ts_create.acc_to_lin,
                   "-n", str(args.num_threads)]
    if args.trim_align:
        trainer_cmd.append("--trim_align")
    if ts_create.stage_status("train"):
        train(trainer_cmd)
    else:
        logging.info("Skipping training:\n$ treesapp train" + ' '.join(trainer_cmd))

    ##
    # Finish validating the file and append the reference package build parameters to the master table
    ##
    if ts_create.stage_status("update"):
        if args.fast:
            marker_package.tree_tool = "FastTree"
        else:
            marker_package.tree_tool = "RAxML"
        marker_package.pfit = create_refpkg.parse_model_parameters(ts_create.var_output_dir + "placement_trainer" +
                                                                   os.sep + "placement_trainer_results.txt")
        ts_create.ref_pkg.validate(marker_package.num_reps)
        if not marker_package.num_reps:
            marker_package.num_reps = ts_create.ref_pkg.num_seqs
        if not marker_package.pfit:
            logging.warning("Linear regression parameters could not be estimated. " +
                            "Taxonomic ranks will not be distance-adjusted during classification for this package.\n")
            marker_package.pfit = [0.0, 7.0]
        create_refpkg.update_build_parameters(param_file, marker_package)
        ts_create.print_terminal_commands()

    logging.info("Data for " + ts_create.ref_pkg.prefix + " has been generated successfully.\n")
    ts_create.remove_intermediates()

    return


def update(sys_args):
    parser = TreeSAPPArgumentParser(description='Update a TreeSAPP reference package with newly identified sequences.')
    add_update_arguments(parser)
    args = parser.parse_args(sys_args)

    ts_updater = Updater()
    ts_updater.furnish_with_arguments(args)
    ts_updater.check_previous_output(args)

    log_file_name = args.output + os.sep + "TreeSAPP_updater_log.txt"
    prep_logging(log_file_name, args.verbose)
    logging.info("\n##\t\t\tUpdating TreeSAPP reference package\t\t\t##\n")

    check_parser_arguments(args, sys_args)
    marker_build_dict = file_parsers.parse_ref_build_params(ts_updater.treesapp_dir, [])
    check_updater_arguments(ts_updater, args, marker_build_dict)
    ts_updater.validate_continue(args)
    ts_updater.ref_pkg.gather_package_files(ts_updater.refpkg_dir, ts_updater.molecule_type, "hierarchical")
    ts_updater.ref_pkg.validate()
    ref_seq_lineage_info = file_parsers.tax_ids_file_to_leaves(ts_updater.ref_pkg.lineage_ids)

    ##
    # Pull out sequences from TreeSAPP output
    ##
    classified_fasta = fasta.FASTA(ts_updater.query_sequences)  # These are the classified sequences
    classified_fasta.load_fasta()
    classified_lines = file_parsers.read_marker_classification_table(ts_updater.assignment_table)
    high_likelihood_seqs = update_refpkg.filter_by_lwr(classified_lines, args.min_lwr)
    resolved_seqs = update_refpkg.filter_by_lineage_depth(classified_lines,
                                                          ts_updater.rank_depth_map[args.min_taxonomic_rank])
    candidate_update_seqs = high_likelihood_seqs.intersection(resolved_seqs)
    classified_targets = utilities.match_target_marker(ts_updater.ref_pkg.prefix, classified_fasta.get_seq_names())
    name_map = update_refpkg.strip_assigment_pattern(classified_fasta.get_seq_names(), ts_updater.ref_pkg.prefix)
    classified_targets = update_refpkg.intersect_incomparable_lists(classified_targets, candidate_update_seqs, name_map)
    # Remove classified sequences that are already in the reference package
    update_refpkg.drop_queries_by_accession(classified_targets, ref_seq_lineage_info)

    if len(classified_targets) == 0:
        logging.error("No new candidate reference sequences. Skipping update.\n")
        return
    classified_fasta.change_dict_keys("original")

    ##
    # Filter out sequences that shouldn't be used in the update: different refpkg, too short, low LWR, etc.
    ##
    classified_fasta.keep_only(classified_targets)
    logging.info(classified_fasta.summarize_fasta_sequences())
    hmm_length = utilities.get_hmm_length(ts_updater.ref_pkg.profile)
    # Use the smaller of the minimum sequence length or 2/3 HMM profile to remove sequence fragments
    if args.min_seq_length < 0.66*hmm_length:
        ts_updater.min_length = int(round(0.66*hmm_length))
        logging.debug("New minimum sequence length threshold set to 2/3 of HMM length (" +
                      str(ts_updater.min_length) + ") instead of " + str(args.min_seq_length) + "\n")
    else:
        ts_updater.min_length = args.min_seq_length
    classified_fasta.remove_shorter_than(ts_updater.min_length)
    if classified_fasta.n_seqs() == 0:
        logging.error("No classified sequences exceed minimum length threshold of " + str(ts_updater.min_length) + ".\n")
        return

    ##
    # Add lineages - use taxa if provided with a table mapping contigs to taxa, TreeSAPP-assigned taxonomy otherwise
    ##
    classified_seq_lineage_map = dict()
    # need_lineage_list = set(classified_fasta.header_registry.keys())  # TreeSAPP IDs that still need lineages
    querying_classified_fasta = classified_fasta
    if ts_updater.seq_names_to_taxa:
        seq_lineage_map = file_parsers.read_seq_taxa_table(ts_updater.seq_names_to_taxa)
        lineage_map, mapped_treesapp_ids = ts_updater.map_orf_lineages(seq_lineage_map,
                                                                       querying_classified_fasta.header_registry)
        classified_seq_lineage_map.update(lineage_map)
        for ts_num in mapped_treesapp_ids:
            querying_classified_fasta.header_registry.pop(ts_num)
    if args.skip_assign:
        name_map = update_refpkg.strip_assigment_pattern(querying_classified_fasta.get_seq_names(),
                                                         ts_updater.ref_pkg.prefix)
        querying_classified_fasta.synchronize_seqs_n_headers()
        querying_classified_fasta.swap_headers(name_map)
        fasta_records = ts_updater.fetch_entrez_lineages(querying_classified_fasta, args.molecule)
        entrez_utils.fill_ref_seq_lineages(fasta_records, classified_seq_lineage_map)
        deduped = []
        for treesapp_id in sorted(querying_classified_fasta.header_registry.keys(), key=int):
            try:
                record = fasta_records[treesapp_id]  # type: entrez_utils.EntrezRecord
            except KeyError:
                deduped.append(treesapp_id)
                continue
            classified_seq_lineage_map[record.accession] = record.lineage
        if deduped:
            logging.warning(str(len(deduped)) + " sequences were not assigned a taxonomic lineage.\n" +
                            "This should match the number of accessions deduplicated while fetching lineage information.\n")
            for treesapp_id in deduped:
                logging.debug("Unable to find '" + treesapp_id + "' in fasta records. More info:\n" +
                              querying_classified_fasta.header_registry[treesapp_id].original + "\n")
                querying_classified_fasta.header_registry.pop(treesapp_id)
            deduped.clear()
            querying_classified_fasta.synchronize_seqs_n_headers()
    else:
        # Map candidate reference sequence names to their TreeSAPP-assigned taxonomies
        assignments = file_parsers.parse_assignments(classified_lines)
        classified_seq_lineage_map.update(update_refpkg.map_classified_seqs(ts_updater.ref_pkg.prefix,
                                                                            assignments,
                                                                            querying_classified_fasta.get_seq_names()))
    classified_seq_indices = classified_fasta.get_seq_names("num")

    ref_header_map = {leaf.number + '_' + ts_updater.ref_pkg.prefix: leaf.description for leaf in ref_seq_lineage_info}
    ref_header_map = update_refpkg.reformat_ref_seq_descriptions(ref_header_map)
    ref_seq_lineage_map = {ref_header_map[leaf.number + '_' + ts_updater.ref_pkg.prefix]:
                           leaf.lineage for leaf in ref_seq_lineage_info}
    num_assigned_candidates = len(classified_seq_lineage_map)
    num_ref_seqs = len(ref_seq_lineage_map)
    classified_seq_lineage_map.update({ref_header_map[leaf.number + '_' + ts_updater.ref_pkg.prefix].split(' ')[0]:
                                       leaf.lineage for leaf in ref_seq_lineage_info})
    diff = num_ref_seqs + num_assigned_candidates - len(classified_seq_lineage_map)
    if diff > 0:
        logging.warning(str(diff) + "candidate sequences are already in the reference package.\n")
    elif diff < 0:
        logging.error("Something's not adding up between the reference (%d), candidate (%d) and complete (%d) "
                      "sequence collections. Reference and candidate should sum to equal complete.\n" %
                      (num_ref_seqs, num_assigned_candidates, len(classified_seq_lineage_map)))
        sys.exit()

    update_refpkg.validate_mixed_lineages(classified_seq_lineage_map)
    utilities.prepend_deep_rank(classified_seq_lineage_map)

    utilities.write_dict_to_table(classified_seq_lineage_map, ts_updater.lineage_map_file)

    ref_fasta = fasta.FASTA(ts_updater.ref_pkg.msa)
    ref_fasta.load_fasta()
    # Update the original reference headers using info from the tax_ids file
    ref_fasta.swap_headers(ref_header_map)
    ref_fasta.custom_lineage_headers(ref_seq_lineage_map)

    classified_fasta.update(ref_fasta.fasta_dict, False)
    classified_fasta.unalign()
    
    if args.resolve:
        classified_fasta.change_dict_keys("num")
        # Write a FASTA for clustering containing the formatted headers since
        # not all clustering tools + versions keep whole header - spaces are replaced with underscores
        fasta.write_new_fasta(fasta_dict=classified_fasta.fasta_dict,
                              fasta_name=ts_updater.cluster_input)
        wrapper.cluster_sequences(ts_updater.executables["usearch"], ts_updater.cluster_input,
                                  ts_updater.uclust_prefix, ts_updater.prop_sim)
        ts_updater.uc = ts_updater.uclust_prefix + ".uc"

        cluster_dict = file_parsers.read_uc(ts_updater.uc)

        # Revert headers in cluster_dict from 'formatted' back to 'original'
        fasta.rename_cluster_headers(cluster_dict, classified_fasta.header_registry)
        logging.debug("\t" + str(len(cluster_dict.keys())) + " sequence clusters\n")

        # Calculate LCA of each cluster to represent the taxonomy of the representative sequence
        entrez_records = update_refpkg.simulate_entrez_records(classified_fasta, classified_seq_lineage_map)
        create_refpkg.cluster_lca(cluster_dict, entrez_records, classified_fasta.header_registry)

        # Ensure centroids are the original reference sequences and skip clusters with identical lineages
        collapsed = update_refpkg.prefilter_clusters(cluster_dict, entrez_records,
                                                     list(ref_fasta.original_header_map().keys()))
        if collapsed:
            logging.warning(str(len(collapsed)) + " original reference sequences removed while resolving:\n\t" +
                            "\n\t".join(collapsed) + "\n")
        # Allow user to select the representative sequence based on organism name, sequence length and similarity
        entrez_records = create_refpkg.present_cluster_rep_options(cluster_dict, entrez_records,
                                                                   classified_fasta.header_registry,
                                                                   classified_fasta.amendments, True)
        # Set all the newly classified candidate reference sequences to cluster_reps to make sure they're still included
        update_refpkg.break_clusters(entrez_records, classified_seq_indices)

        # Remove sequences that were replaced by resolve from ts_updater.old_ref_fasta
        still_repping = []
        for num_id in entrez_records:
            ref_seq = entrez_records[num_id]  # type: entrez_utils.EntrezRecord
            if ref_seq.cluster_rep:
                still_repping.append(ref_seq.versioned + ' ' + ref_seq.description)
        ref_fasta.keep_only(still_repping, True)  # This removes the original reference sequences to be replaced
        classified_fasta.change_dict_keys("original")
        classified_fasta.keep_only(still_repping)

    # Write only the sequences that have been properly classified
    classified_fasta.change_dict_keys("original")
    fasta.write_new_fasta(classified_fasta.fasta_dict, ts_updater.combined_fasta)
    fasta.write_new_fasta(ref_fasta.fasta_dict, ts_updater.old_ref_fasta)

    ##
    # Call create to create a new, updated reference package where the new sequences are guaranteed
    ##
    create_cmd = ["-i", ts_updater.combined_fasta,
                  "-c", ts_updater.ref_pkg.prefix,
                  "-p", str(ts_updater.prop_sim),
                  "-m", ts_updater.molecule_type,
                  "--guarantee", ts_updater.old_ref_fasta,
                  "-o", ts_updater.output_dir,
                  "--accession2lin", ts_updater.lineage_map_file,
                  "--num_procs", str(args.num_threads)]
    if args.trim_align:
        create_cmd.append("--trim_align")
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
    create(create_cmd)

    ##
    # Summarize some key parts of the new reference package, compared to the old one
    ##
    new_hmm_length = utilities.get_hmm_length(ts_updater.output_dir + "final_outputs" + os.sep +
                                              ts_updater.ref_pkg.prefix + ".hmm")
    logging.debug("\tOld HMM length = " + str(hmm_length) + "\n" +
                  "\tNew HMM length = " + str(new_hmm_length) + "\n")

    return


def layer(sys_args):
    # STAGE 1: Prompt the user and prepare files and lists for the pipeline
    parser = TreeSAPPArgumentParser(description="This script is generally used for layering extra annotations "
                                                "beyond taxonomy (such as Subgroup or Metabolism) to TreeSAPP outputs."
                                                " This is accomplished by adding an extra column (to all rows) of an "
                                                "existing marker_contig_map.tsv and annotating the relevant sequences")
    add_layer_arguments(parser)
    args = parser.parse_args(sys_args)

    ts_layer = Layerer()

    log_file_name = args.output + os.sep + "TreeSAPP_layer_log.txt"
    prep_logging(log_file_name, args.verbose)
    logging.info("\n##\t\t\t\tLayering extra annotations on TreeSAPP classifications\t\t\t\t##\n\n")

    annotate_extra.check_arguments(ts_layer, args)

    ##
    # Worklow:
    #   1. Read data/tree_data/ref_build_parameters.tsv to get marker codes, denominators, and more (oh my!)
    #   2. Read the marker_contig_map.tsv file from the output directory to create the master data structure
    #   3. For each of the colours_styles files provided (potentially multiple for the same marker):
    #       3.1) Add the annotation variable to master_dat for every sequence (instantiate with "NA")
    #       3.2) Read the .jplace file for every sequence classified as marker
    #       3.3) Add the annotation information to every sequence classified as marker in master_dat
    #   4. Write the new classification file called "extra_annotated_marker_contig_map.tsv"
    ##
    marker_subgroups = dict()
    unique_markers_annotated = set()
    marker_tree_info = dict()
    internal_nodes = dict()
    marker_build_dict = file_parsers.parse_ref_build_params(ts_layer.treesapp_dir, [])
    master_dat, field_order = annotate_extra.parse_marker_classification_table(ts_layer.final_output_dir +
                                                                               "marker_contig_map.tsv")
    # structure of master dat:
    # {"Sequence_1": {"Field1": x, "Field2": y, "Extra": n},
    #  "Sequence_2": {"Field1": i, "Field2": j, "Extra": n}}
    for annot_f in ts_layer.annot_files:
        # Determine the marker being annotated
        marker = data_type = refpkg = ""
        for code in marker_build_dict:
            marker = marker_build_dict[code].cog
            annot_marker_re = re.compile(r"^{0}_(\w+).txt$".format(marker))
            if annot_marker_re.match(os.path.basename(annot_f)):
                data_type = annot_marker_re.match(os.path.basename(annot_f)).group(1)
                refpkg = code
                break
            else:
                marker = data_type = refpkg = ""
        if marker not in master_dat.keys():
            continue
        if marker and data_type:
            unique_markers_annotated.add(refpkg)
            if data_type not in marker_subgroups:
                marker_subgroups[data_type] = dict()
                internal_nodes[data_type] = dict()
            marker_subgroups[data_type][marker], internal_nodes[data_type][marker] = file_parsers.read_colours_file(annot_f)
        else:
            logging.warning("Unable to parse the marker and/or annotation type from " + annot_f + ".\n" +
                            "Is it possible this reference package is not in " +
                            ts_layer.treesapp_dir + os.sep + "data" + os.sep + "ref_build_parameters.tsv?\n")
    # Instantiate every query sequence in marker_contig_map with an empty string for each data_type
    for data_type in marker_subgroups:
        for marker in master_dat:
            for assignment in master_dat[marker]:  # type: annotate_extra.ClassifiedSequence
                assignment.layers[data_type] = "NA"
    # Update the field_order dictionary with new fields
    field_acc = len(field_order)
    for new_datum in sorted(marker_subgroups.keys()):
        field_order[field_acc] = new_datum
        field_acc += 1

    # Load the query sequence annotations
    for data_type in marker_subgroups:
        if data_type not in marker_tree_info:
            marker_tree_info[data_type] = dict()
        for refpkg_code in unique_markers_annotated:
            marker = marker_build_dict[refpkg_code].cog
            jplace = os.sep.join([ts_layer.treesapp_output, "iTOL_output", marker, marker + "_complete_profile.jplace"])

            if marker in marker_subgroups[data_type]:
                # Create the dictionary mapping an internal node to all leaves
                internal_node_map = entish.map_internal_nodes_leaves(jplace_parser(jplace).tree)
                # Routine for exchanging any organism designations for their respective node number
                tax_ids_file = ts_layer.tree_dir + "tax_ids_" + marker + ".txt"
                taxa_map = file_parsers.tax_ids_file_to_leaves(tax_ids_file)
                clusters = annotate_extra.names_for_nodes(marker_subgroups[data_type][marker],
                                                          internal_node_map,
                                                          taxa_map)

                if not internal_nodes[data_type][marker]:
                    # Convert the leaf node ranges to internal nodes for consistency
                    clusters = utilities.convert_outer_to_inner_nodes(clusters, internal_node_map)

                marker_tree_info[data_type][marker], leaves_in_clusters = annotate_extra.annotate_internal_nodes(internal_node_map,
                                                                                                                 clusters)
                diff = len(taxa_map) - len(leaves_in_clusters)
                if diff != 0:
                    unannotated = set()
                    for inode in internal_node_map:
                        for leaf in internal_node_map[inode]:
                            if leaf not in leaves_in_clusters:
                                unannotated.add(str(leaf))
                    logging.warning("The following leaf nodes were not mapped to annotation groups:\n" +
                                    "\t" + ', '.join(sorted(unannotated, key=int)) + "\n")
            else:
                pass
    marker_subgroups.clear()
    master_dat = annotate_extra.map_queries_to_annotations(marker_tree_info, master_dat)
    annotate_extra.write_classification_table(ts_layer.final_output_dir, field_order, master_dat)

    return


def assign(sys_args):
    # STAGE 1: Prompt the user and prepare files and lists for the pipeline
    parser = TreeSAPPArgumentParser(description='Taxonomically classify sequences through evolutionary placement.')
    add_classify_arguments(parser)
    args = parser.parse_args(sys_args)

    ts_assign = Assigner()
    ts_assign.furnish_with_arguments(args)
    ts_assign.check_previous_output(args)

    log_file_name = args.output + os.sep + "TreeSAPP_classify_log.txt"
    prep_logging(log_file_name, args.verbose)
    logging.info("\n##\t\t\t\tAssigning sequences with TreeSAPP\t\t\t\t##\n\n")

    check_parser_arguments(args, sys_args)
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
        ts_assign.query_sequences = ts_assign.aa_orfs_file
    else:
        ts_assign.query_sequences = ts_assign.input_sequences

    query_seqs = fasta.FASTA(ts_assign.query_sequences)
    # Read the query sequences provided and (by default) write a new FASTA file with formatted headers
    if ts_assign.stage_status("clean"):
        logging.info("Reading and formatting " + ts_assign.query_sequences + "... ")
        query_seqs.fasta_dict = fasta.format_read_fasta(ts_assign.query_sequences, "prot")
        query_seqs.header_registry = fasta.register_headers(fasta.get_headers(ts_assign.query_sequences), True)
        query_seqs.change_dict_keys("num")
        logging.info("done.\n")
        logging.info("Writing formatted FASTA file to " + ts_assign.formatted_input + "... ")
        fasta.write_new_fasta(query_seqs.fasta_dict, ts_assign.formatted_input)
        logging.info("done.\n")
    else:
        ts_assign.formatted_input = ts_assign.query_sequences
        query_seqs.load_fasta()
        query_seqs.change_dict_keys("num")  # Swap the formatted headers for the numerical IDs for quick look-ups
    logging.info("\tTreeSAPP will analyze the " + str(len(query_seqs.fasta_dict)) + " sequences found in input.\n")

    ##
    # STAGE 3: Run hmmsearch on the query sequences to search for marker homologs
    ##
    if ts_assign.stage_status("search"):
        hmm_domtbl_files = wrapper.hmmsearch_orfs(ts_assign.executables["hmmsearch"], ts_assign.hmm_dir,
                                                  marker_build_dict, ts_assign.formatted_input,
                                                  ts_assign.var_output_dir, args.num_threads)
        hmm_matches = file_parsers.parse_domain_tables(args, hmm_domtbl_files)
        extracted_seq_dict, numeric_contig_index = extract_hmm_matches(hmm_matches, query_seqs.fasta_dict)
        numeric_contig_index = replace_contig_names(numeric_contig_index, query_seqs)
        homolog_seq_files = write_grouped_fastas(extracted_seq_dict, numeric_contig_index,
                                                 marker_build_dict, ts_assign.var_output_dir)

    ##
    # STAGE 4: Run hmmalign or PaPaRa, and optionally BMGE, to produce the MSAs required to for the ML estimations
    ##
    if ts_assign.stage_status("align"):
        create_ref_phy_files(ts_assign.aln_dir, ts_assign.var_output_dir,
                             homolog_seq_files, marker_build_dict, ref_alignment_dimensions)
        concatenated_msa_files = multiple_alignments(ts_assign.executables, ts_assign.refpkg_dir,
                                                     ts_assign.var_output_dir, homolog_seq_files, marker_build_dict,
                                                     "hmmalign", args.num_threads)
        file_type = utilities.find_msa_type(concatenated_msa_files)
        alignment_length_dict = get_sequence_counts(concatenated_msa_files, ref_alignment_dimensions,
                                                    args.verbose, file_type)

        if args.trim_align:
            tool = "BMGE"
            trimmed_mfa_files = wrapper.filter_multiple_alignments(ts_assign.executables, concatenated_msa_files,
                                                                   marker_build_dict, args.num_threads, tool)
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
        # TODO: Replace this merge_fasta_dicts_by_index with FASTA - only necessary for writing the classified sequences
        extracted_seq_dict = fasta.merge_fasta_dicts_by_index(extracted_seq_dict, numeric_contig_index)
        fasta.write_classified_sequences(tree_saps, extracted_seq_dict, ts_assign.classified_aa_seqs)
        abundance_dict = dict()
        for refpkg_code in tree_saps:
            for placed_seq in tree_saps[refpkg_code]:  # type: TreeProtein
                abundance_dict[placed_seq.contig_name] = 1.0
        if args.molecule == "dna":
            if os.path.isfile(ts_assign.nuc_orfs_file):
                nuc_orfs = fasta.FASTA(ts_assign.nuc_orfs_file)
                nuc_orfs.load_fasta()
                nuc_orfs.change_dict_keys()
                if not os.path.isfile(ts_assign.classified_nuc_seqs):
                    logging.info("Creating nucleotide FASTA file of classified sequences '" +
                                 ts_assign.classified_nuc_seqs + "'... ")
                    fasta.write_classified_sequences(tree_saps, nuc_orfs.fasta_dict, ts_assign.classified_nuc_seqs)
                    logging.info("done.\n")
            else:
                logging.warning("Unable to read '" + ts_assign.nuc_orfs_file + "'.\n" +
                                "Cannot create the nucleotide FASTA file of classified sequences!\n")
            if args.rpkm:
                abundance_args = ["--treesapp_output", ts_assign.output_dir,
                                  "--reads", args.reads,
                                  "--pairing", args.pairing,
                                  "--num_procs", str(args.num_threads)]
                if args.reverse:
                    abundance_args += ["--reverse", args.reverse]
                abundance_dict = abundance(abundance_args)
                summarize_placements_rpkm(tree_saps, abundance_dict, marker_build_dict, ts_assign.final_output_dir)

        abundify_tree_saps(tree_saps, abundance_dict)
        assign_out = ts_assign.final_output_dir + os.sep + "marker_contig_map.tsv"
        write_tabular_output(tree_saps, tree_numbers_translation, marker_build_dict, ts_assign.sample_prefix, assign_out)
        produce_itol_inputs(tree_saps, marker_build_dict, itol_data, ts_assign.output_dir, ts_assign.refpkg_dir)
        delete_files(args.delete, ts_assign.var_output_dir, 4)

    delete_files(args.delete, ts_assign.var_output_dir, 5)

    return


def abundance(sys_args):
    """
    TreeSAPP subcommand that is used to add read-inferred abundance information (e.g. FPKM, TPM) to classified sequences
    Command requires:

1. Path to TreeSAPP output directory that contains classified sequences (FASTA-format) in the final_outputs/
2. Path to read file(s) in FASTQ format
3. Parameters indicating whether a) the reads are paired-end or single-end and b) the FASTQ is interleaved

    With these arguments and the option `--report update` TreeSAPP would run BWA MEM and samsum to
    calculate the desired abundance values (FPKM by default) and update the classification table with the
    Sample (first column) of the classification table matching the prefix of the FASTQ.

    Optionally, `treesapp abundance` can be called with `--report nothing` and a dictionary containing the abundance
    values would be returned.

    :param sys_args: treesapp abundance arguments with the treesapp subcommand removed
    :return: A dictionary containing the abundance values indexed by the reference sequence (e.g. ORF, contig) names
    """
    parser = TreeSAPPArgumentParser(description="Validate the functional purity of a reference package.")
    add_abundance_arguments(parser)
    args = parser.parse_args(sys_args)

    ts_abund = Abundance()
    ts_abund.furnish_with_arguments(args)
    abundance_dict = {}

    log_file_name = args.output + os.sep + "TreeSAPP_purity_log.txt"
    prep_logging(log_file_name, args.verbose)
    logging.info("\n##\t\t\tCalculating abundance of classified sequences\t\t\t##\n")

    check_parser_arguments(args, sys_args)
    marker_build_dict = file_parsers.parse_ref_build_params(ts_abund.treesapp_dir)
    tree_numbers_translation = file_parsers.read_species_translation_files(ts_abund.treesapp_dir, marker_build_dict)
    ts_abund.check_arguments(args)
    # TODO: Implement check-pointing for abundance
    # ts_abund.validate_continue(args)

    sam_file = align_reads_to_nucs(ts_abund.executables["bwa"], ts_abund.classified_nuc_seqs,
                                   ts_abund.var_output_dir, args.reads, args.pairing, args.reverse, args.num_threads)
    ts_abund.sample_prefix = ts_abund.fq_suffix_re.sub('', '.'.join(os.path.basename(args.reads).split('.')[:-1]))

    if os.path.isfile(sam_file):
        ref_seq_abunds = ref_sequence_abundances(aln_file=sam_file, seq_file=ts_abund.classified_nuc_seqs,
                                                 min_aln=10, p_cov=50, map_qual=1, multireads=False)
        for ref_name, ref_seq in ref_seq_abunds.items():
            abundance_dict[re.sub(r"\|(.*){2,10}\|\d+_\d+$", '', ref_seq.name)] = ref_seq.fpkm
        ref_seq_abunds.clear()
    else:
        logging.warning("SAM file '%s' was not generated.\n" % sam_file)
        return abundance_dict

    # TODO: Add delete argument to io, use args.delete instead
    delete_files(True, ts_abund.var_output_dir, 4)

    # TODO: Index each TreeProtein's abundance by the dataset name, write a new row for each dataset's abundance
    if args.report != "nothing" and os.path.isfile(ts_abund.classifications):
        assignments = file_parsers.read_marker_classification_table(ts_abund.classifications)
        # Convert assignments to TreeProtein instances
        tree_saps = ts_abund.assignments_to_treesaps(assignments, marker_build_dict)
        summarize_placements_rpkm(tree_saps, abundance_dict, marker_build_dict, ts_abund.final_output_dir)
        write_tabular_output(tree_saps, tree_numbers_translation, marker_build_dict, ts_abund.sample_prefix,
                             ts_abund.classifications)

    return abundance_dict


def purity(sys_args):
    """

    :param sys_args:
    :return:
    """
    parser = TreeSAPPArgumentParser(description="Validate the functional purity of a reference package.")
    add_purity_arguments(parser)
    args = parser.parse_args(sys_args)

    ts_purity = Purity()
    ts_purity.furnish_with_arguments(args)
    ts_purity.check_previous_output(args)

    log_file_name = args.output + os.sep + "TreeSAPP_purity_log.txt"
    prep_logging(log_file_name, args.verbose)
    logging.info("\n##\t\t\tBeginning purity analysis\t\t\t##\n")

    check_parser_arguments(args, sys_args)
    marker_build_dict = file_parsers.parse_ref_build_params(ts_purity.treesapp_dir)
    check_purity_arguments(ts_purity, args, marker_build_dict)
    ts_purity.validate_continue(args)

    # Load FASTA data
    ref_seqs = fasta.FASTA(args.input)
    ref_seqs.load_fasta()

    if ts_purity.stage_status("assign"):
        assign_args = ["-i", ts_purity.input_sequences, "-o", ts_purity.assign_dir,
                       "-m", ts_purity.molecule_type, "-n", str(args.num_threads),
                       "-t", ts_purity.refpkg_build.denominator,
                       "--overwrite", "--delete"]
        try:
            assign(assign_args)
        except:  # Just in case treesapp assign fails, just continue
            logging.error("TreeSAPP failed.\n")

    if ts_purity.stage_status("summarize"):
        metadat_dict = dict()
        # Parse classification table and identify the groups that were assigned
        if os.path.isfile(ts_purity.classifications):
            assigned_lines = file_parsers.read_marker_classification_table(ts_purity.classifications)
            ts_purity.assignments = file_parsers.parse_assignments(assigned_lines)
        else:
            logging.error("marker_contig_map.tsv is missing from output directory '" +
                          os.path.dirname(ts_purity.classifications) + "'\n" +
                          "Please remove this directory and re-run.\n")
            sys.exit(5)

        logging.info("\nSummarizing assignments for reference package " + ts_purity.refpkg_build.cog + "\n")
        # If an information table was provided, map the metadata to classified markers
        if ts_purity.metadata_file:
            metadat_dict.update(ts_purity.load_metadata())
        # Identify the number of sequences that are descendents of each orthologous group
        jplace_file = os.sep.join([ts_purity.assign_dir, "iTOL_output", ts_purity.refpkg_build.cog,
                                   ts_purity.refpkg_build.cog + "_complete_profile.jplace"])
        jplace_data = jplace_parser(jplace_file)
        tree_placement_queries = demultiplex_pqueries(jplace_data)
        placement_tree = jplace_data.tree
        node_map = entish.map_internal_nodes_leaves(placement_tree)
        ortholog_map = ts_purity.assign_leaves_to_orthologs(tree_placement_queries, node_map)
        ts_purity.summarize_groups_assigned(ortholog_map, metadat_dict)

        # Write each sequence name that can be assigned to an ortholog to the log
        summary_str = ""
        for ortholog_name in sorted(ortholog_map, key=lambda x: len(ortholog_map[x])):
            summary_str += ortholog_name + ":\n\t"
            summary_str += "\n\t".join(ortholog_map[ortholog_name]) + "\n"
        logging.debug(summary_str)

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

    ts_evaluate = Evaluator()
    ts_evaluate.furnish_with_arguments(args)
    ts_evaluate.check_previous_output(args)

    log_file_name = args.output + os.sep + "TreeSAPP_evaluation_log.txt"
    prep_logging(log_file_name, args.verbose)
    logging.info("\n##\t\t\tBeginning clade exclusion analysis\t\t\t##\n")

    check_parser_arguments(args, sys_args)
    marker_build_dict = file_parsers.parse_ref_build_params(ts_evaluate.treesapp_dir, ts_evaluate.targets)
    check_evaluate_arguments(ts_evaluate, args, marker_build_dict)
    ts_evaluate.validate_continue(args)
    load_rank_depth_map(ts_evaluate)
    refpkg_name = ts_evaluate.target_marker.cog

    ref_leaves = file_parsers.tax_ids_file_to_leaves(ts_evaluate.tree_dir + "tax_ids_" + args.reference_marker + ".txt")
    ref_lineages = dict()
    for leaf in ref_leaves:
        ref_lineages[leaf.number] = leaf.lineage

    # Load FASTA data
    ref_seqs = fasta.FASTA(args.input)
    ref_seqs.fasta_dict = fasta.format_read_fasta(ref_seqs.file, args.molecule)
    if args.length:
        for seq_id in ref_seqs.fasta_dict:
            if len(ref_seqs.fasta_dict[seq_id]) < args.length:
                logging.warning(seq_id + " sequence is shorter than " + str(args.length) + "\n")
            else:
                max_stop = len(ref_seqs.fasta_dict[seq_id]) - args.length
                random_start = randint(0, max_stop)
                ref_seqs.fasta_dict[seq_id] = ref_seqs.fasta_dict[seq_id][random_start:random_start + args.length]
    ref_seqs.header_registry = fasta.register_headers(fasta.get_headers(ref_seqs.file))

    fasta_records = ts_evaluate.fetch_entrez_lineages(ref_seqs, args.molecule, args.acc_to_taxid)
    entrez_utils.fill_ref_seq_lineages(fasta_records, ts_evaluate.seq_lineage_map)

    fasta_records = utilities.remove_elongated_lineages(fasta_records)

    logging.info("Selecting representative sequences for each taxon.\n")

    # Filter the sequences from redundant taxonomic lineages, picking up to 5 representative sequences
    representative_seqs, ts_evaluate.taxa_filter = pick_taxonomic_representatives(fasta_records,
                                                                                  ts_evaluate.taxa_filter)
    deduplicated_fasta_dict = select_rep_seqs(representative_seqs, fasta_records)
    fasta.write_new_fasta(deduplicated_fasta_dict, ts_evaluate.test_rep_taxa_fasta)
    rep_accession_lineage_map = map_seqs_to_lineages(ts_evaluate.seq_lineage_map, deduplicated_fasta_dict)

    # Checkpoint three: We have accessions linked to taxa, and sequences to analyze with TreeSAPP, but not classified
    if ts_evaluate.stage_status("classify"):
        # Run TreeSAPP against the provided tax_ids file and the unique taxa FASTA file
        if args.length:
            min_seq_length = str(min(args.length - 10, 30))
        else:
            min_seq_length = str(30)

        validate_ref_package_files(ts_evaluate.treesapp_dir, refpkg_name, ts_evaluate.output_dir)

        ranks = {"Kingdom": 0, "Phylum": 1, "Class": 2, "Order": 3, "Family": 4, "Genus": 5, "Species": 6}
        for rank in args.taxon_rank:
            leaf_trimmed_taxa_map = trim_lineages_to_rank(ref_lineages, rank)
            unique_ref_lineages = sorted(set(leaf_trimmed_taxa_map.values()))
            unique_query_lineages = sorted(set(trim_lineages_to_rank(rep_accession_lineage_map, rank).values()))
            depth = ranks[rank]
            for lineage in unique_query_lineages:
                # Check number one: Is the optimal placement in the pruned reference tree?
                optimal_lca_taxonomy = "; ".join(lineage.split("; ")[:-1])
                if optimal_lca_taxonomy not in ["; ".join(tl.split("; ")[:-1]) for tl in unique_ref_lineages
                                                if tl != lineage]:
                    logging.debug("Optimal placement target '" + optimal_lca_taxonomy + "' not in pruned tree.\n")
                    continue

                # Select representative sequences belonging to the taxon being tested
                taxon_rep_seqs = select_rep_seqs(representative_seqs, fasta_records, lineage)
                # Check number 2: Decide whether to continue analyzing taxon based on number of query sequences
                if len(taxon_rep_seqs.keys()) == 0:
                    logging.debug("No query sequences for " + lineage + ".\n")
                    continue

                # Continuing with classification
                # Refpkg input files in ts_evaluate.var_output_dir/refpkg_name/rank_tax/
                # Refpkg built in ts_evaluate.var_output_dir/refpkg_name/rank_tax/{refpkg_name}_{rank_tax}.gpkg/
                taxon = re.sub(r"([ /])", '_', lineage.split("; ")[-1])
                rank_tax = rank[0] + '_' + taxon
                intermediates_path = ts_evaluate.var_output_dir + refpkg_name + os.sep + rank_tax + os.sep
                if not os.path.isdir(intermediates_path):
                    os.makedirs(intermediates_path)

                logging.info("Classifications for the " + rank + " '" + taxon + "' put " + intermediates_path + "\n")
                test_obj = ts_evaluate.new_taxa_test(rank, lineage)
                test_obj.queries = taxon_rep_seqs.keys()
                test_rep_taxa_fasta = intermediates_path + rank_tax + ".fa"
                test_refpkg_prefix = refpkg_name + '_' + rank_tax
                classifier_output = intermediates_path + args.tool + "_output" + os.sep

                if args.tool in ["graftm", "diamond"]:
                    tax_ids_file = intermediates_path + os.sep + refpkg_name + "_taxonomy.csv"
                    classification_table = classifier_output + rank_tax + os.sep + rank_tax + "_read_tax.tsv"
                    gpkg_path = intermediates_path + test_refpkg_prefix + ".gpkg"

                    if not os.path.isfile(classification_table):
                        # GraftM refpkg input paths:
                        filtered_gpkg_tax_ids = intermediates_path + "tax_ids_" + refpkg_name + ".txt"
                        filtered_mfa = intermediates_path + refpkg_name + ".mfa"
                        filtered_fasta = intermediates_path + refpkg_name + ".fa"
                        # GraftM refpkg output files:
                        gpkg_refpkg_path = gpkg_path + os.sep + test_refpkg_prefix + ".gpkg.refpkg" + os.sep
                        gpkg_tax_ids_file = gpkg_refpkg_path + refpkg_name + "_taxonomy.csv"

                        # Copy reference files, then exclude all clades belonging to the taxon being tested
                        prep_graftm_ref_files(treesapp_dir=ts_evaluate.treesapp_dir,
                                              intermediate_dir=intermediates_path,
                                              target_taxon=lineage,
                                              marker=ts_evaluate.target_marker,
                                              depth=depth)

                        if not os.path.isdir(gpkg_path):
                            build_graftm_package(gpkg_path=gpkg_path,
                                                 tax_file=filtered_gpkg_tax_ids,
                                                 mfa_file=filtered_mfa,
                                                 fa_file=filtered_fasta,
                                                 threads=args.num_threads)
                        shutil.copy(gpkg_tax_ids_file, tax_ids_file)
                        # Write the query sequences
                        fasta.write_new_fasta(taxon_rep_seqs, test_rep_taxa_fasta)

                        graftm_classify(test_rep_taxa_fasta,
                                        gpkg_path,
                                        classifier_output,
                                        args.num_threads, args.tool)

                        if not os.path.isfile(classification_table):
                            # The TaxonTest object is maintained for record-keeping (to track # queries & classifieds)
                            logging.warning("GraftM did not generate output for " + lineage + ". Skipping.\n")
                            shutil.rmtree(intermediates_path)
                            continue

                    test_obj.taxonomic_tree = lca_calculations.grab_graftm_taxa(tax_ids_file)
                    graftm_assignments = file_parsers.read_graftm_classifications(classification_table)
                    test_obj.assignments = {refpkg_name: graftm_assignments}
                    test_obj.filter_assignments(refpkg_name)
                else:
                    tax_ids_file = intermediates_path + "tax_ids_" + refpkg_name + ".txt"
                    classification_table = classifier_output + "final_outputs" + os.sep + "marker_contig_map.tsv"

                    if not os.path.isfile(classification_table) or not os.path.isfile(tax_ids_file):
                        # Copy reference files, then exclude all clades belonging to the taxon being tested
                        prefix = exclude_clade_from_ref_files(ts_evaluate.treesapp_dir, refpkg_name,
                                                              ts_evaluate.var_output_dir + refpkg_name + os.sep,
                                                              lineage, depth,
                                                              ts_evaluate.executables, args.fresh, args.molecule)
                        # Write the query sequences
                        fasta.write_new_fasta(taxon_rep_seqs, test_rep_taxa_fasta)
                        assign_args = ["-i", test_rep_taxa_fasta, "-o", classifier_output,
                                       "-m", ts_evaluate.molecule_type, "-n", str(args.num_threads),
                                       "--min_seq_length", str(min_seq_length), "--overwrite", "--delete"]
                        if args.trim_align:
                            assign_args.append("--trim_align")
                        try:
                            assign(assign_args)
                        except:  # Just in case treesapp assign fails, just continue
                            pass
                        restore_reference_package(ts_evaluate.treesapp_dir, prefix,
                                                  intermediates_path, refpkg_name)
                        if not os.path.isfile(classification_table):
                            # The TaxonTest object is maintained for record-keeping (to track # queries & classifieds)
                            logging.warning("TreeSAPP did not generate output for " + lineage + ". Skipping.\n")
                            shutil.rmtree(classifier_output)
                            continue
                    else:
                        # Valid number of queries and these sequences have already been classified
                        pass

                    test_obj.taxonomic_tree = lca_calculations.all_possible_assignments(tax_ids_file)
                    if os.path.isfile(classification_table):
                        assigned_lines = file_parsers.read_marker_classification_table(classification_table)
                        test_obj.assignments = file_parsers.parse_assignments(assigned_lines)
                        test_obj.filter_assignments(refpkg_name)
                        test_obj.distances = parse_distances(assigned_lines)
                    else:
                        logging.error("marker_contig_map.tsv is missing from output directory '" +
                                      os.path.dirname(classification_table) + "'\n" +
                                      "Please remove this directory and re-run.\n")
                        sys.exit(21)
        # TODO: Currently emits warning for GraftM and DIAMOND - only needed when running TreeSAPP
        remove_clade_exclusion_files(ts_evaluate.var_output_dir + refpkg_name + os.sep)

    if ts_evaluate.stage_status("calculate"):
        # everything has been prepared, only need to parse the classifications and map lineages
        logging.info("Finishing up the mapping of classified, filtered taxonomic sequences.\n")
        for rank in sorted(ts_evaluate.taxa_tests):
            for test_obj in ts_evaluate.taxa_tests[rank]:  # type: TaxonTest
                if test_obj.assignments:
                    marker_assignments = map_headers_to_lineage(test_obj.assignments, fasta_records)
                    # Return the number of correct, classified, and total sequences of that taxon at the current rank
                    # Identify the excluded rank for each query sequence
                    if len(marker_assignments) == 0:
                        logging.debug("No sequences were classified for " + test_obj.taxon + "\n")
                        continue

                    for marker in marker_assignments:
                        ts_evaluate.markers.add(marker)
                    rank_assignments = lca_calculations.identify_excluded_clade(marker_assignments,
                                                                                test_obj.taxonomic_tree,
                                                                                refpkg_name)
                    for a_rank in rank_assignments:
                        if a_rank != rank and len(rank_assignments[a_rank]) > 0:
                            logging.warning(
                                rank + "-level clade excluded but optimal classification was found to be " + a_rank +
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

        if ts_evaluate.target_marker.denominator not in ts_evaluate.markers and refpkg_name not in ts_evaluate.markers:
            logging.error("No sequences were classified as " + refpkg_name + ".\n")
            sys.exit(21)

        # For debugging:
        # for rank in ts_evaluate.ranks:
        #     distal, pendant, tip = ts_evaluate.summarize_rankwise_distances(rank)

        # Determine the specificity for each rank
        clade_exclusion_strings = ts_evaluate.get_classification_performance()
        ts_evaluate.taxonomic_recall_table()
        ts_evaluate.taxonomic_recall_tree()
        ts_evaluate.write_performance_table(clade_exclusion_strings, args.tool)
        ts_evaluate.summarize_taxonomic_diversity()
        containment_strings = determine_containment(ts_evaluate)
        ts_evaluate.write_containment_table(containment_strings, args.tool)

    return
