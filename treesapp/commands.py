
import logging
import sys
import re
import os
import shutil
from random import randint

from collections import namedtuple
from samsum.commands import ref_sequence_abundances

from treesapp import entrez_utils
from treesapp import file_parsers
from treesapp import fasta
from treesapp import create_refpkg
from treesapp import utilities
from treesapp import wrapper
from treesapp import entish
from treesapp import lca_calculations
from treesapp import placement_trainer
from treesapp import update_refpkg
from treesapp import annotate_extra
from treesapp import treesapp_args
from treesapp import classy
from treesapp.phylo_seq import TreeProtein, assignments_to_treesaps
from treesapp.refpkg import ReferencePackage, view, edit
from treesapp.training_utils import train_classification_filter
from treesapp.assign import abundify_tree_saps, delete_files, prep_reference_packages_for_assign,\
    get_alignment_dims, bin_hmm_matches, write_grouped_fastas, create_ref_phy_files,\
    multiple_alignments, get_sequence_counts, check_for_removed_sequences, determine_confident_lineage,\
    evaluate_trimming_performance, parse_raxml_output, filter_placements, align_reads_to_nucs, select_query_placements,\
    summarize_placements_rpkm, write_classification_table, produce_itol_inputs, replace_contig_names, read_refpkg_tax_ids
from treesapp.jplace_utils import sub_indices_for_seq_names_jplace, jplace_parser, demultiplex_pqueries
from treesapp.clade_exclusion_evaluator import pick_taxonomic_representatives, select_rep_seqs,\
    map_seqs_to_lineages, prep_graftm_ref_files, build_graftm_package, map_headers_to_lineage, graftm_classify,\
    determine_containment, parse_distances, load_rank_depth_map, get_testable_lineages_for_rank


def info(sys_args):
    """

    """
    parser = treesapp_args.TreeSAPPArgumentParser(description="Return package, executable and refpkg information.")
    treesapp_args.add_info_arguments(parser)
    args = parser.parse_args(sys_args)
    classy.prep_logging()
    ts_info = classy.TreeSAPP("info")

    import treesapp
    import Bio
    import numpy
    import scipy
    import ete3
    import sklearn
    import joblib
    import seaborn
    import samsum
    import pyfastx
    import pygtrie
    logging.info("TreeSAPP version " + treesapp.__version__ + ".\n")

    # Write the version of all python deps
    py_deps = {"biopython": Bio.__version__,
               "ete3": ete3.__version__,
               "joblib": joblib.__version__,
               "numpy": numpy.__version__,
               "scipy": scipy.__version__,
               "scikit-learn": sklearn.__version__,
               "seaborn": seaborn.__version__,
               "samsum": samsum.__version__,
               "pyfastx": pyfastx.version(),
               "pygtrie": pygtrie.__version__}

    logging.info("Python package dependency versions:\n\t" +
                 "\n\t".join([k + ": " + v for k, v in py_deps.items()]) + "\n")

    # Write the version of executable deps
    ts_info.furnish_with_arguments(args)
    if args.refpkg_dir:
        ts_info.refpkg_dir = args.refpkg_dir
    logging.info(utilities.executable_dependency_versions(ts_info.executables))

    if args.verbose:
        refpkg_dict = file_parsers.gather_ref_packages(ts_info.refpkg_dir)
        refpkg_summary_str = "\t".join(["Name", "Code-name",
                                        "Molecule", "Tree builder", "RefPkg-type", "Leaf nodes",
                                        "Description", "Created", "Last-updated"]) + "\n"
        for refpkg_name in sorted(refpkg_dict, key=lambda x: refpkg_dict[x].prefix):
            refpkg = refpkg_dict[refpkg_name]  # type: ReferencePackage
            refpkg_summary_str += "\t".join([refpkg_name, refpkg.refpkg_code,
                                             refpkg.molecule, refpkg.tree_tool, refpkg.kind, str(refpkg.num_seqs),
                                             refpkg.description, refpkg.date, refpkg.update]
                                            ) + "\n"
        logging.info(refpkg_summary_str)

    return


def package(sys_args):
    """
    Perform specific operations on JSON-formatted single-file reference packages

    :param sys_args: Arguments from the command-line
    :return: None
    """
    pkg_usage = """
treesapp package <subcommand> <attributes> [<args>]
** Subcommands include:
view        Print reference package attributes to the console
edit        Change reference package attributes
**
Use '-h' to get subcommand-specific help, e.g.
"""
    parser = treesapp_args.TreeSAPPArgumentParser(description='Facilitate operations on reference packages')
    parser.add_argument("subcommand", nargs='?', choices=["view", "edit"],
                        help="A subcommand specifying the type of operation to perform")
    args = parser.parse_args(sys_args[0:1])
    if not args.subcommand:
        sys.stderr.write(pkg_usage)
        sys.exit(1)

    refpkg = ReferencePackage()

    treesapp_args.add_package_arguments(parser, refpkg.get_public_attributes())
    args = parser.parse_args(sys_args)

    if not args.output:
        args.output = os.path.dirname(args.pkg_path)
    if not os.path.isdir(args.output):
        os.mkdir(args.output)

    classy.prep_logging(os.path.join(args.output, 'TreeSAPP_package_log.txt'))

    refpkg.f__json = args.pkg_path
    refpkg.slurp()

    if args.subcommand == "view":
        view(refpkg, args.attributes)
    elif args.subcommand == "edit":
        edit(refpkg, args.attributes, args.output, args.overwrite)
    else:
        logging.error("Unrecognized command: '{}'.\n{}\n".format(args.subcommand, pkg_usage))
        sys.exit(1)

    return


def train(sys_args):
    parser = treesapp_args.TreeSAPPArgumentParser(description='Model evolutionary distances across taxonomic ranks.')
    treesapp_args.add_trainer_arguments(parser)
    args = parser.parse_args(sys_args)

    ts_trainer = classy.PhyTrainer()
    ts_trainer.furnish_with_arguments(args)
    ts_trainer.check_previous_output(args.overwrite)

    log_file_name = args.output + os.sep + "TreeSAPP_trainer_log.txt"
    classy.prep_logging(log_file_name, args.verbose)
    logging.info("\n##\t\t\tTrain taxonomic rank-placement distance model\t\t\t##\n")

    treesapp_args.check_parser_arguments(args, sys_args)
    treesapp_args.check_trainer_arguments(ts_trainer, args)
    ts_trainer.validate_continue(args)
    ts_trainer.ref_pkg.validate()
    ts_trainer.ref_pkg.disband(os.path.join(ts_trainer.var_output_dir, ts_trainer.ref_pkg.prefix + "_RefPkg"))

    train_seqs = fasta.FASTA(ts_trainer.input_sequences)

    if ts_trainer.stage_status("search"):
        # Read the FASTA into a dictionary - homologous sequences will be extracted from this
        train_seqs.fasta_dict = fasta.format_read_fasta(ts_trainer.input_sequences, ts_trainer.molecule_type)
        train_seqs.header_registry = fasta.register_headers(fasta.get_headers(ts_trainer.input_sequences))

        logging.info("Searching for domain sequences... ")
        hmm_domtbl_files = wrapper.run_hmmsearch(ts_trainer.executables["hmmsearch"],
                                                 ts_trainer.ref_pkg.f__search_profile,
                                                 ts_trainer.input_sequences,
                                                 ts_trainer.var_output_dir)
        logging.info("done.\n")
        hmm_matches = file_parsers.parse_domain_tables(args, hmm_domtbl_files)
        marker_gene_dict = dict()
        hmm_extract_seqs = utilities.extract_hmm_matches(hmm_matches, train_seqs.fasta_dict, train_seqs.header_registry)
        for k, v in hmm_extract_seqs.items():
            marker_gene_dict.update(v)
        logging.info(train_seqs.summarize_fasta_sequences())
        fasta.write_new_fasta(marker_gene_dict, ts_trainer.hmm_purified_seqs)
        marker_gene_dict.clear()
        utilities.hmm_pile(hmm_matches)
    else:
        train_seqs.load_fasta()
        train_seqs.change_dict_keys("formatted")
        ts_trainer.hmm_purified_seqs = ts_trainer.input_sequences

    ts_trainer.fetch_entrez_lineages(train_seqs, args.molecule, args.acc_to_taxid)
    rank_depth_map = ts_trainer.ref_pkg.taxa_trie.accepted_ranks_depths
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
                                                                                       ts_trainer.var_output_dir,
                                                                                       ts_trainer.training_ranks,
                                                                                       args.num_threads)
        # Write the tab-delimited file with metadata included for each placement
        placement_trainer.write_placement_table(pqueries, ts_trainer.placement_table, ts_trainer.ref_pkg.prefix)
    else:
        pfit_array = placement_trainer.complete_regression(taxa_evo_dists, ts_trainer.training_ranks)
        if pfit_array:
            logging.info("Placement distance model complete.\n")
        else:
            logging.info("Unable to complete phylogenetic distance and rank correlation.\n")

    # Reformat the pqueries dictionary for classifier training and testing
    refpkg_pqueries = placement_trainer.flatten_pquery_dict(pqueries, ts_trainer.ref_pkg.prefix)
    tp_names = {ts_trainer.ref_pkg.prefix:
                    [pquery.contig_name for pquery in refpkg_pqueries[ts_trainer.ref_pkg.prefix]]}
    # Train the one-class SVM model
    refpkg_classifiers = train_classification_filter(refpkg_pqueries, tp_names,
                                                     refpkg_map={ts_trainer.ref_pkg.prefix: ts_trainer.ref_pkg})

    if ts_trainer.stage_status("update"):
        ts_trainer.ref_pkg.pfit = pfit_array
        if not ts_trainer.ref_pkg.pfit:
            logging.warning("Linear regression parameters could not be estimated. " +
                            "Taxonomic ranks will not be distance-adjusted during classification for this package.\n")
            ts_trainer.ref_pkg.pfit = [0.0, 7.0]

        ts_trainer.ref_pkg.svc = refpkg_classifiers[ts_trainer.ref_pkg.prefix]
        ts_trainer.ref_pkg.f__json = os.path.join(ts_trainer.final_output_dir,
                                                  os.path.basename(ts_trainer.ref_pkg.f__json))

        ts_trainer.ref_pkg.validate()
        ts_trainer.ref_pkg.pickle_package()

    # Write the text file containing distances used in the regression analysis
    with open(ts_trainer.placement_summary, 'w') as out_handler:
        trained_string = "Regression parameters = " + re.sub(' ', '', str(pfit_array)) + "\n"
        for rank in sorted(rank_depth_map, key=lambda x: rank_depth_map[x]):
            trained_string += "# " + rank + "\n"
            if rank in taxa_evo_dists:
                trained_string += str(sorted(taxa_evo_dists[rank], key=float)) + "\n"
            trained_string += "\n"
        out_handler.write(trained_string)

    if args.delete and os.path.isdir(ts_trainer.var_output_dir):
        shutil.rmtree(ts_trainer.var_output_dir)

    return


def create(sys_args):
    parser = treesapp_args.TreeSAPPArgumentParser(description="Create a reference package for TreeSAPP.")
    treesapp_args.add_create_arguments(parser)
    args = parser.parse_args(sys_args)

    ts_create = classy.Creator()
    ts_create.furnish_with_arguments(args)
    ts_create.check_previous_output(args.overwrite)

    log_file_name = args.output + os.sep + "TreeSAPP_create_" + args.refpkg_name + "_log.txt"
    classy.prep_logging(log_file_name, args.verbose)
    logging.info("\n##\t\t\tCreating TreeSAPP reference package\t\t\t##\n")

    treesapp_args.check_parser_arguments(args, sys_args)
    treesapp_args.check_create_arguments(ts_create, args)
    ts_create.validate_continue(args)

    # Populate the final TreeSAPP reference file paths with the proper output directory
    ts_create.ref_pkg.update_file_names()
    ts_create.ref_pkg.change_file_paths(ts_create.final_output_dir)
    ts_create.ref_pkg.cmd = ' '.join(["treesapp", "create"] + sys_args)

    ref_seqs = fasta.FASTA(args.input)

    if ts_create.stage_status("search"):
        profile_match_dict = dict()
        # Read the FASTA into a dictionary - homologous sequences will be extracted from this
        ref_seqs.fasta_dict = fasta.format_read_fasta(args.input, ts_create.ref_pkg.molecule)
        ref_seqs.header_registry = fasta.register_headers(fasta.get_headers(args.input))
        logging.debug("Raw, unfiltered sequence summary:\n" + ref_seqs.summarize_fasta_sequences())

        logging.info("Searching for domain sequences... ")
        hmm_domtbl_files = wrapper.run_hmmsearch(ts_create.executables["hmmsearch"], ts_create.hmm_profile,
                                                 args.input, ts_create.var_output_dir, args.num_threads)
        logging.info("done.\n")
        hmm_matches = file_parsers.parse_domain_tables(args, hmm_domtbl_files)
        for k, v in utilities.extract_hmm_matches(hmm_matches, ref_seqs.fasta_dict, ref_seqs.header_registry).items():
            profile_match_dict.update(v)
        fasta.write_new_fasta(profile_match_dict, ts_create.hmm_purified_seqs)
        profile_match_dict.clear()

        utilities.hmm_pile(hmm_matches)
    else:
        ts_create.hmm_purified_seqs = ts_create.input_sequences

    ##
    # Synchronize records between fasta_dict and header_registry (e.g. short ones may be removed by format_read_fasta())
    ##
    ref_seqs.file = ts_create.hmm_purified_seqs
    ref_seqs.fasta_dict = fasta.format_read_fasta(fasta_input=ref_seqs.file, molecule=ts_create.ref_pkg.molecule,
                                                  min_seq_length=args.min_seq_length)
    ref_seqs.header_registry = fasta.register_headers(fasta.get_headers(ref_seqs.file))
    ref_seqs.synchronize_seqs_n_headers()
    # Get rid of ambiguity or unusual characters
    ref_seqs.replace_ambiguity_chars(ts_create.ref_pkg.molecule)
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
    fasta_records = ts_create.fetch_entrez_lineages(ref_seqs, ts_create.ref_pkg.molecule,
                                                    args.acc_to_taxid, args.seq_names_to_taxa)
    entrez_utils.fill_ref_seq_lineages(fasta_records, ts_create.seq_lineage_map)
    create_refpkg.strip_rank_prefix_from_organisms(fasta_records, ts_create.ref_pkg.taxa_trie)
    prefilter_ref_seqs = entrez_utils.entrez_record_snapshot(fasta_records)

    if ts_create.stage_status("clean"):
        # Remove the sequences failing 'filter' and/or only retain the sequences in 'screen'
        fasta_records = create_refpkg.screen_filter_taxa(fasta_records, args.screen, args.filter, ref_seqs.amendments)
        # Remove the sequence records with low resolution lineages, according to args.min_taxonomic_rank
        # TODO: Replace this function with one offered by TaxonomicHierarchy
        fasta_records = create_refpkg.remove_by_truncated_lineages(fasta_records,
                                                                   args.min_taxonomic_rank, ref_seqs.amendments)
        # Ensure there are no records with redundant headers and sequences
        fasta_records = classy.dedup_records(ref_seqs, fasta_records)

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
                                      ts_create.uclust_prefix, ts_create.ref_pkg.pid)
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
            create_refpkg.present_cluster_rep_options(cluster_dict, fasta_records, ref_seqs.header_registry,
                                                      ref_seqs.amendments)
        elif ts_create.uc and args.headless:
            create_refpkg.finalize_cluster_reps(cluster_dict, fasta_records, ref_seqs.header_registry)
        else:
            pass

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

        ts_create.ref_pkg.lineage_ids = create_refpkg.lineages_to_dict(fasta_replace_dict, args.taxa_lca)

        postfilter_ref_seqs = entrez_utils.entrez_record_snapshot(fasta_records)
        filtered_ref_seqs = utilities.dict_diff(prefilter_ref_seqs, postfilter_ref_seqs)
        logging.debug("{0} references before and {1} remaining after filtering.\n".format(len(prefilter_ref_seqs),
                                                                                          len(postfilter_ref_seqs)))
        ts_create.ref_pkg.taxa_trie.jetison_taxa_from_hierarchy(filtered_ref_seqs)

        taxonomic_summary = create_refpkg.summarize_reference_taxa(fasta_replace_dict, ts_create.ref_pkg.taxa_trie,
                                                                   args.taxa_lca)
        logging.info(taxonomic_summary)

        ##
        # Perform multiple sequence alignment
        ##
        if args.multiple_alignment:
            create_refpkg.create_new_ref_fasta(ts_create.unaln_ref_fasta, fasta_replace_dict, True)
        else:
            create_refpkg.create_new_ref_fasta(ts_create.unaln_ref_fasta, fasta_replace_dict)

        if ts_create.ref_pkg.molecule == 'rrna':
            create_refpkg.generate_cm_data(args, ts_create.unaln_ref_fasta)
            args.multiple_alignment = True
        elif args.multiple_alignment is False:
            logging.info("Aligning the sequences using MAFFT... ")
            create_refpkg.run_mafft(ts_create.executables["mafft"],
                                    ts_create.unaln_ref_fasta, ts_create.ref_pkg.f__msa, args.num_threads)
            logging.info("done.\n")
        else:
            pass
        ref_seqs.file = ts_create.ref_pkg.f__msa
        ref_seqs.load_fasta()
        ts_create.ref_pkg.num_seqs = ref_seqs.n_seqs()
        n_rows, n_cols = fasta.multiple_alignment_dimensions(mfa_file=ts_create.ref_pkg.f__msa,
                                                             seq_dict=ref_seqs.fasta_dict)
        logging.debug("Reference alignment contains " +
                      str(n_rows) + " sequences with " +
                      str(n_cols) + " character positions.\n")

        ##
        # Build the HMM profile from the aligned reference FASTA file
        ##
        if ts_create.ref_pkg.molecule == "rrna":
            # A .cm file has already been generated, no need for HMM
            pass
        else:
            wrapper.build_hmm_profile(ts_create.executables["hmmbuild"],
                                      ts_create.ref_pkg.f__msa,
                                      ts_create.ref_pkg.f__profile)
            ts_create.ref_pkg.hmm_length()

        ts_create.ref_pkg.band()
        ts_create.ref_pkg.dereplicate_hmm(dereplication_rank="genus",
                                          hmmbuild_exe=ts_create.executables["hmmbuild"],
                                          mafft_exe=ts_create.executables["mafft"],
                                          n_threads=args.num_threads, intermediates_dir=ts_create.var_output_dir)

        ##
        # Optionally trim with BMGE and create the Phylip multiple alignment file
        ##
        dict_for_phy = dict()
        if args.trim_align:
            trimmed_mfa_files = wrapper.filter_multiple_alignments(ts_create.executables,
                                                                   {ts_create.ref_pkg.refpkg_code:
                                                                    [ts_create.ref_pkg.f__msa]},
                                                                   {ts_create.ref_pkg.refpkg_code:
                                                                    ts_create.ref_pkg})
            trimmed_mfa_file = trimmed_mfa_files[ts_create.ref_pkg.refpkg_code]
            unique_ref_headers = set(ref_seqs.fasta_dict.keys())
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
            # NOTE: only a single trimmed-MSA file in the dictionary
            for trimmed_msa_file in qc_ma_dict:
                dict_for_phy = qc_ma_dict[trimmed_msa_file]
                os.remove(trimmed_msa_file)
        else:
            dict_for_phy.update(ref_seqs.fasta_dict)

        phy_dict = utilities.reformat_fasta_to_phy(dict_for_phy)
        utilities.write_phy_file(ts_create.phylip_file, phy_dict)

        ##
        # Build the tree using either RAxML-NG or FastTree
        ##
        ts_create.ref_pkg.infer_phylogeny(ts_create.phylip_file, ts_create.executables, ts_create.phy_dir,
                                          args.bootstraps, args.num_threads, args.raxml_model)

    ts_create.ref_pkg.band()
    # Build the regression model of placement distances to taxonomic ranks
    trainer_cmd = ["-i", ts_create.filtered_fasta,
                   "-r", ts_create.ref_pkg.f__json,
                   "-o", ts_create.training_dir,
                   "-m", ts_create.ref_pkg.molecule,
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
        ts_create.ref_pkg.f__json = os.path.join(ts_create.var_output_dir, "placement_trainer", "final_outputs",
                                                 ts_create.ref_pkg.prefix + ts_create.ref_pkg.refpkg_suffix)
        ts_create.ref_pkg.slurp()
        ts_create.ref_pkg.validate()
        ts_create.ref_pkg.change_file_paths(ts_create.final_output_dir)
        ts_create.ref_pkg.pickle_package()

    ts_create.remove_intermediates(args.delete)
    ts_create.print_terminal_commands()

    return


def update(sys_args):
    parser = treesapp_args.TreeSAPPArgumentParser(description='Update a TreeSAPP reference package with assigned sequences.')
    treesapp_args.add_update_arguments(parser)
    args = parser.parse_args(sys_args)

    ts_updater = classy.Updater()
    ts_updater.furnish_with_arguments(args)
    ts_updater.check_previous_output(args.overwrite)

    log_file_name = args.output + os.sep + "TreeSAPP_updater_log.txt"
    classy.prep_logging(log_file_name, args.verbose)
    logging.info("\n##\t\t\tUpdating TreeSAPP reference package\t\t\t##\n")

    treesapp_args.check_parser_arguments(args, sys_args)
    treesapp_args.check_updater_arguments(ts_updater, args)
    ts_updater.validate_continue(args)
    ts_updater.ref_pkg.validate()
    ref_seq_lineage_info = ts_updater.ref_pkg.generate_tree_leaf_references_from_refpkg()

    ##
    # Pull out sequences from TreeSAPP output
    ##
    classified_fasta = fasta.FASTA(ts_updater.query_sequences)  # These are the classified sequences
    classified_fasta.load_fasta()
    classified_lines = file_parsers.read_marker_classification_table(ts_updater.assignment_table)
    candidate_update_seqs = update_refpkg.filter_by_lwr(classified_lines, args.min_lwr)
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
    hmm_length = utilities.get_hmm_length(ts_updater.ref_pkg.f__profile)
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
    querying_classified_fasta = fasta.FASTA("")
    querying_classified_fasta.clone(classified_fasta)

    if args.skip_assign:
        name_map = update_refpkg.strip_assigment_pattern(querying_classified_fasta.get_seq_names(),
                                                         ts_updater.ref_pkg.prefix)
        querying_classified_fasta.synchronize_seqs_n_headers()
        querying_classified_fasta.swap_headers(name_map)
        fasta_records = ts_updater.fetch_entrez_lineages(ref_seqs=querying_classified_fasta, molecule=args.molecule,
                                                         seqs_to_lineage=ts_updater.seq_names_to_taxa)
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

    ref_fasta = fasta.FASTA(ts_updater.ref_pkg.f__msa)
    ref_fasta.load_fasta()
    # Update the original reference headers using info from the tax_ids file
    ref_fasta.swap_headers(ref_header_map)
    ref_fasta.custom_lineage_headers(ref_seq_lineage_map)

    classified_fasta.update(ref_fasta.fasta_dict, False)
    classified_fasta.unalign()
    
    if args.resolve:
        ##
        # The purpose of this block is to remove any former candidate reference sequences from the ref_fasta object
        # that have a more truncated lineage that the new candidate reference sequences in classified_fasta
        ##
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
            logging.warning("{} original reference sequences removed while resolving:\n\t"
                            "{}\n".format(len(collapsed), "\n\t".join(collapsed)))

        # Ensure the EntrezRecord with the most resolved lineage is the representative
        update_refpkg.resolve_cluster_lineages(cluster_dict, entrez_records, ts_updater.ref_pkg.taxa_trie)

        if args.headless:
            create_refpkg.finalize_cluster_reps(cluster_dict, entrez_records, classified_fasta.header_registry)
        else:
            # Allow user to select the representative sequence based on organism name, sequence length and similarity
            create_refpkg.present_cluster_rep_options(cluster_dict, entrez_records, classified_fasta.header_registry,
                                                      classified_fasta.amendments, True)

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
                  "--num_procs", str(args.num_threads),
                  "--bootstraps", str(args.bootstraps)]
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
    ts_updater.updated_refpkg.f__json = ts_updater.updated_refpkg_path
    ts_updater.updated_refpkg.slurp()
    ts_updater.update_refpkg_fields()

    return


def layer(sys_args):
    # STAGE 1: Prompt the user and prepare files and lists for the pipeline
    parser = treesapp_args.TreeSAPPArgumentParser(description="This script is generally used for layering extra "
                                                              "annotations such as Subgroup or Metabolism on TreeSAPP "
                                                              "outputs. This is accomplished by adding an extra column "
                                                              "(to all rows) of an existing marker_contig_map.tsv and "
                                                              "annotating the relevant sequences.")
    treesapp_args.add_layer_arguments(parser)
    args = parser.parse_args(sys_args)

    ts_layer = classy.Layerer()

    log_file_name = args.output + os.sep + "TreeSAPP_layer_log.txt"
    classy.prep_logging(log_file_name, args.verbose)
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
    master_dat, field_order = annotate_extra.parse_marker_classification_table(ts_layer.final_output_dir +
                                                                               "marker_contig_map.tsv")
    refpkg_dict = file_parsers.gather_ref_packages(ts_layer.refpkg_dir)
    tree_numbers_translation = read_refpkg_tax_ids(refpkg_dict)

    # structure of master dat:
    # {"Sequence_1": {"Field1": x, "Field2": y, "Extra": n},
    #  "Sequence_2": {"Field1": i, "Field2": j, "Extra": n}}
    for annot_f in ts_layer.annot_files:
        # Determine the marker being annotated
        data_type, refpkg = "", ""
        for refpkg_name in refpkg_dict:
            annot_marker_re = re.compile(r"^{0}_(\w+).txt$".format(refpkg_name))
            if annot_marker_re.match(os.path.basename(annot_f)):
                data_type = annot_marker_re.match(os.path.basename(annot_f)).group(1)
                refpkg = refpkg_name
                break
            else:
                data_type, refpkg = "", ""
        if refpkg not in master_dat.keys():
            continue
        if refpkg and data_type:
            unique_markers_annotated.add(refpkg)
            if data_type not in marker_subgroups:
                marker_subgroups[data_type] = dict()
                internal_nodes[data_type] = dict()
            marker_subgroups[data_type][refpkg], internal_nodes[data_type][refpkg] = file_parsers.read_colours_file(annot_f,
                                                                                                                    refpkg)
        else:
            logging.warning("Unable to parse the reference package name and/or annotation type from {}.\n"
                            "Is it possible this reference package is not in {}?".format(annot_f, ts_layer.refpkg_dir))
    # Instantiate every query sequence in marker_contig_map with an empty string for each data_type
    for data_type in marker_subgroups:
        for refpkg_name in master_dat:
            for assignment in master_dat[refpkg_name]:  # type: annotate_extra.ClassifiedSequence
                assignment.layers[data_type] = "NA"
    # Update the field_order dictionary with new fields
    field_acc = len(field_order)
    for new_datum in sorted(marker_subgroups.keys()):
        field_order[field_acc] = new_datum
        field_acc += 1

    # Load the query sequence annotations
    for data_type in marker_subgroups:
        logging.info("Annotating '%s' classifications for the following reference package(s):\n" % data_type)
        if data_type not in marker_tree_info:
            marker_tree_info[data_type] = dict()
        for refpkg_name in unique_markers_annotated:
            jplace = os.path.join(ts_layer.treesapp_output, "iTOL_output", refpkg_name,
                                  refpkg_name + "_complete_profile.jplace")

            if refpkg_name in marker_subgroups[data_type]:
                logging.info("\t" + refpkg_name + "\n")
                # Create the dictionary mapping an internal node to all leaves
                internal_node_map = entish.map_internal_nodes_leaves(jplace_parser(jplace).tree)
                # Routine for exchanging any organism designations for their respective node number
                taxa_map = tree_numbers_translation[refpkg_name]

                if internal_nodes[data_type][refpkg_name]:
                    clusters = annotate_extra.names_for_nodes(marker_subgroups[data_type][refpkg_name],
                                                              internal_node_map,
                                                              taxa_map)
                else:
                    # Convert the leaf node ranges to internal nodes for consistency
                    clusters = utilities.convert_outer_to_inner_nodes(marker_subgroups[data_type][refpkg_name],
                                                                      internal_node_map)

                marker_tree_info[data_type][refpkg_name], leaves_in_clusters = annotate_extra.annotate_internal_nodes(internal_node_map,
                                                                                                                 clusters)
                diff = len(taxa_map) - len(leaves_in_clusters)
                if diff != 0:
                    unannotated = set()
                    for inode in internal_node_map:
                        for leaf in internal_node_map[inode]:
                            if leaf not in leaves_in_clusters:
                                unannotated.add(str(leaf))
                    logging.warning("The following leaf nodes were not mapped to annotation groups:\n" +
                                    "\t" + ', '.join(sorted(unannotated, key=lambda x: int(x.split('_')[0]))) + "\n")
            else:
                pass
    marker_subgroups.clear()
    master_dat = annotate_extra.map_queries_to_annotations(marker_tree_info, master_dat)
    annotate_extra.write_classification_table(ts_layer.final_output_dir, field_order, master_dat)

    return


def assign(sys_args):
    # STAGE 1: Prompt the user and prepare files and lists for the pipeline
    parser = treesapp_args.TreeSAPPArgumentParser(description='Classify sequences through evolutionary placement.')
    treesapp_args.add_classify_arguments(parser)
    args = parser.parse_args(sys_args)

    ts_assign = classy.Assigner()
    ts_assign.furnish_with_arguments(args)
    ts_assign.check_previous_output(args.overwrite)

    log_file_name = args.output + os.sep + "TreeSAPP_classify_log.txt"
    classy.prep_logging(log_file_name, args.verbose)
    logging.info("\n##\t\t\t\tAssigning sequences with TreeSAPP\t\t\t\t##\n\n")

    treesapp_args.check_parser_arguments(args, sys_args)
    treesapp_args.check_classify_arguments(ts_assign, args)
    ts_assign.validate_continue(args)

    refpkg_dict = file_parsers.gather_ref_packages(ts_assign.refpkg_dir, ts_assign.target_refpkgs)
    prep_reference_packages_for_assign(refpkg_dict, ts_assign.var_output_dir)
    ref_alignment_dimensions = get_alignment_dims(refpkg_dict)
    tree_numbers_translation = read_refpkg_tax_ids(refpkg_dict)

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
        hmm_domtbl_files = wrapper.hmmsearch_orfs(ts_assign.executables["hmmsearch"],
                                                  refpkg_dict, ts_assign.formatted_input,
                                                  ts_assign.var_output_dir, args.num_threads, args.max_e)
        hmm_matches = file_parsers.parse_domain_tables(args, hmm_domtbl_files)
        extracted_seq_dict, numeric_contig_index = bin_hmm_matches(hmm_matches, query_seqs.fasta_dict)
        numeric_contig_index = replace_contig_names(numeric_contig_index, query_seqs)
        homolog_seq_files = write_grouped_fastas(extracted_seq_dict, numeric_contig_index,
                                                 refpkg_dict, ts_assign.var_output_dir)
        # TODO: Replace this merge_fasta_dicts_by_index with FASTA - only necessary for writing the classified sequences
        extracted_seq_dict = fasta.merge_fasta_dicts_by_index(extracted_seq_dict, numeric_contig_index)
        delete_files(args.delete, ts_assign.var_output_dir, 1)
    ##
    # STAGE 4: Run hmmalign or PaPaRa, and optionally BMGE, to produce the MSAs required to for the ML estimations
    ##
    combined_msa_files = dict()
    split_msa_files = dict()
    if ts_assign.stage_status("align"):
        create_ref_phy_files(refpkg_dict, ts_assign.var_output_dir, homolog_seq_files, ref_alignment_dimensions)
        concatenated_msa_files = multiple_alignments(ts_assign.executables, ts_assign.refpkg_dir,
                                                     ts_assign.var_output_dir, homolog_seq_files, refpkg_dict,
                                                     "hmmalign", args.num_threads)
        file_type = utilities.find_msa_type(concatenated_msa_files)
        alignment_length_dict = get_sequence_counts(concatenated_msa_files, ref_alignment_dimensions,
                                                    args.verbose, file_type)

        if args.trim_align:
            tool = "BMGE"
            trimmed_mfa_files = wrapper.filter_multiple_alignments(ts_assign.executables, concatenated_msa_files,
                                                                   refpkg_dict, args.num_threads, tool)
            qc_ma_dict = check_for_removed_sequences(ts_assign.aln_dir, trimmed_mfa_files, concatenated_msa_files,
                                                     refpkg_dict, args.min_seq_length)
            evaluate_trimming_performance(qc_ma_dict, alignment_length_dict, concatenated_msa_files, tool)
            combined_msa_files.update(qc_ma_dict)
        else:
            combined_msa_files.update(concatenated_msa_files)

        # Subset the multiple alignment of reference sequences and queries to just contain query sequences
        MSAs = namedtuple("MSAs", "ref query")
        for denominator in combined_msa_files:
            split_msa_files[denominator] = []
            for combined_msa in combined_msa_files[denominator]:
                split_msa = MSAs(os.path.dirname(combined_msa) + os.sep +
                                 os.path.basename('.'.join(combined_msa.split('.')[:-1])) + "_references.mfa",
                                 os.path.dirname(combined_msa) + os.sep +
                                 os.path.basename('.'.join(combined_msa.split('.')[:-1])) + "_queries.mfa")
                fasta.split_combined_ref_query_fasta(combined_msa, split_msa.query, split_msa.ref)
                split_msa_files[denominator].append(split_msa)
        combined_msa_files.clear()
        delete_files(args.delete, ts_assign.var_output_dir, 3)

    ##
    # STAGE 5: Run EPA-ng to compute the ML estimations
    ##
    if ts_assign.stage_status("place"):
        wrapper.launch_evolutionary_placement_queries(ts_assign.executables, split_msa_files, refpkg_dict,
                                                      ts_assign.var_output_dir, args.num_threads)
        sub_indices_for_seq_names_jplace(ts_assign.var_output_dir, numeric_contig_index, refpkg_dict)

    if ts_assign.stage_status("classify"):
        tree_saps, itol_data = parse_raxml_output(ts_assign.var_output_dir, refpkg_dict)
        select_query_placements(tree_saps)
        filter_placements(tree_saps, refpkg_dict, ts_assign.clf, ts_assign.tree_dir, args.min_likelihood)
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
                summarize_placements_rpkm(tree_saps, abundance_dict, refpkg_dict, ts_assign.final_output_dir)

        abundify_tree_saps(tree_saps, abundance_dict)
        assign_out = ts_assign.final_output_dir + os.sep + "marker_contig_map.tsv"
        determine_confident_lineage(tree_saps, tree_numbers_translation, refpkg_dict)
        write_classification_table(tree_saps, ts_assign.sample_prefix, assign_out)
        produce_itol_inputs(tree_saps, refpkg_dict, itol_data, ts_assign.output_dir, ts_assign.refpkg_dir)
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
    parser = treesapp_args.TreeSAPPArgumentParser(description="Calculate classified sequence abundances from read coverage.")
    treesapp_args.add_abundance_arguments(parser)
    args = parser.parse_args(sys_args)

    ts_abund = classy.Abundance()
    ts_abund.furnish_with_arguments(args)
    abundance_dict = {}

    log_file_name = args.output + os.sep + "TreeSAPP_abundance_log.txt"
    classy.prep_logging(log_file_name, args.verbose)
    logging.info("\n##\t\t\tCalculating abundance of classified sequences\t\t\t##\n")

    treesapp_args.check_parser_arguments(args, sys_args)
    ts_abund.check_arguments(args)
    refpkg_dict = file_parsers.gather_ref_packages(ts_abund.refpkg_dir)
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

    delete_files(args.delete, ts_abund.var_output_dir, 4)

    # TODO: Index each TreeProtein's abundance by the dataset name, write a new row for each dataset's abundance
    if args.report != "nothing" and os.path.isfile(ts_abund.classifications):
        assignments = file_parsers.read_marker_classification_table(ts_abund.classifications)
        # Convert assignments to TreeProtein instances
        tree_saps = assignments_to_treesaps(assignments, refpkg_dict)
        summarize_placements_rpkm(tree_saps, abundance_dict, refpkg_dict, ts_abund.final_output_dir)
        write_classification_table(tree_saps, ts_abund.sample_prefix, ts_abund.classifications)

    return abundance_dict


def purity(sys_args):
    """

    :param sys_args:
    :return:
    """
    parser = treesapp_args.TreeSAPPArgumentParser(description="Validate the functional purity of a reference package.")
    treesapp_args.add_purity_arguments(parser)
    args = parser.parse_args(sys_args)

    ts_purity = classy.Purity()
    ts_purity.furnish_with_arguments(args)
    ts_purity.check_previous_output(args.overwrite)

    log_file_name = args.output + os.sep + "TreeSAPP_purity_log.txt"
    classy.prep_logging(log_file_name, args.verbose)
    logging.info("\n##\t\t\tBeginning purity analysis\t\t\t##\n")

    treesapp_args.check_parser_arguments(args, sys_args)
    treesapp_args.check_purity_arguments(ts_purity, args)
    ts_purity.validate_continue(args)

    # Load FASTA data
    ref_seqs = fasta.FASTA(ts_purity.input_sequences)
    ref_seqs.load_fasta()
    ref_seqs.change_dict_keys()
    ref_seqs.unalign()
    fasta.write_new_fasta(ref_seqs.fasta_dict, ts_purity.formatted_input)

    if ts_purity.stage_status("assign"):
        assign_args = ["-i", ts_purity.formatted_input, "-o", ts_purity.assign_dir,
                       "-m", ts_purity.molecule_type, "-n", str(args.num_threads),
                       "-t", ts_purity.ref_pkg.prefix, "--refpkg_dir", ts_purity.refpkg_dir,
                       "--overwrite", "--delete", "--no_svm"]
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

        logging.info("\nSummarizing assignments for reference package " + ts_purity.ref_pkg.prefix + "\n")
        # If an information table was provided, map the metadata to classified markers
        if ts_purity.metadata_file:
            metadat_dict.update(ts_purity.load_metadata())
        # Identify the number of sequences that are descendents of each orthologous group
        jplace_file = os.sep.join([ts_purity.assign_dir, "iTOL_output", ts_purity.ref_pkg.prefix,
                                   ts_purity.ref_pkg.prefix + "_complete_profile.jplace"])
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
    Provide it a FASTA file for which it will determine the taxonomic lineage for each sequence
    and run all taxonomic representative sequences with TreeSAPP then analyze via clade exclusion

    :return:
    """
    parser = treesapp_args.TreeSAPPArgumentParser(description='Evaluate classification performance using clade-exclusion analysis.')
    treesapp_args.add_evaluate_arguments(parser)
    args = parser.parse_args(sys_args)

    ts_evaluate = classy.Evaluator()
    ts_evaluate.furnish_with_arguments(args)
    ts_evaluate.check_previous_output(args.overwrite)

    log_file_name = args.output + os.sep + "TreeSAPP_evaluation_log.txt"
    classy.prep_logging(log_file_name, args.verbose)
    logging.info("\n##\t\t\tBeginning clade exclusion analysis\t\t\t##\n")

    treesapp_args.check_parser_arguments(args, sys_args)
    treesapp_args.check_evaluate_arguments(ts_evaluate, args)
    ts_evaluate.validate_continue(args)
    load_rank_depth_map(ts_evaluate)

    ref_leaves = ts_evaluate.ref_pkg.generate_tree_leaf_references_from_refpkg()
    ref_lineages = {leaf.number: leaf.lineage for leaf in ref_leaves}

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

        # validate_ref_package_files(refpkg, ts_evaluate.output_dir)

        ranks = {"domain": 0, "phylum": 1, "class": 2, "order": 3, "family": 4, "genus": 5, "species": 6}
        for rank in args.taxon_rank:
            depth = ranks[rank]
            for lineage in get_testable_lineages_for_rank(ref_lineages, rep_accession_lineage_map, rank):
                # Select representative sequences belonging to the taxon being tested
                taxon_rep_seqs = select_rep_seqs(representative_seqs, fasta_records, lineage)
                # Decide whether to continue analyzing taxon based on number of query sequences
                if len(taxon_rep_seqs.keys()) == 0:
                    logging.debug("No query sequences for {}.\n".format(lineage))
                    continue

                # Continuing with classification
                # Refpkg input files in ts_evaluate.var_output_dir/refpkg_name/rank_tax/
                # Refpkg built in ts_evaluate.var_output_dir/refpkg_name/rank_tax/{refpkg_name}_{rank_tax}.gpkg/
                taxon = re.sub(r"([ /])", '_', lineage.split("; ")[-1])
                intermediates_path = os.path.join(ts_evaluate.var_output_dir, taxon) + os.sep
                if not os.path.isdir(intermediates_path):
                    os.makedirs(intermediates_path)

                logging.info("Classifications for the {} '{}' put {}\n".format(rank, taxon, intermediates_path))
                test_obj = ts_evaluate.new_taxa_test(lineage)
                test_obj.queries = taxon_rep_seqs.keys()
                test_rep_taxa_fasta = intermediates_path + taxon + ".fa"
                classifier_output = intermediates_path + args.tool + "_output" + os.sep

                if args.tool in ["graftm", "diamond"]:
                    test_refpkg_prefix = ts_evaluate.ref_pkg.prefix + '_' + taxon
                    tax_ids_file = intermediates_path + os.sep + ts_evaluate.ref_pkg.prefix + "_taxonomy.csv"
                    classification_table = classifier_output + taxon + os.sep + taxon + "_read_tax.tsv"
                    gpkg_path = intermediates_path + test_refpkg_prefix + ".gpkg"

                    if not os.path.isfile(classification_table):
                        # GraftM refpkg input paths:
                        filtered_gpkg_tax_ids = intermediates_path + "tax_ids_" + ts_evaluate.ref_pkg.prefix + ".txt"
                        filtered_mfa = intermediates_path + ts_evaluate.ref_pkg.prefix + ".mfa"
                        filtered_fasta = intermediates_path + ts_evaluate.ref_pkg.prefix + ".fa"
                        # GraftM refpkg output files:
                        gpkg_refpkg_path = gpkg_path + os.sep + test_refpkg_prefix + ".gpkg.refpkg" + os.sep
                        gpkg_tax_ids_file = gpkg_refpkg_path + ts_evaluate.ref_pkg.prefix + "_taxonomy.csv"

                        # Copy reference files, then exclude all clades belonging to the taxon being tested
                        prep_graftm_ref_files(intermediate_dir=intermediates_path,
                                              target_taxon=lineage,
                                              refpkg=ts_evaluate.ref_pkg,
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
                    test_obj.assignments = {ts_evaluate.ref_pkg.prefix: graftm_assignments}
                    test_obj.filter_assignments(ts_evaluate.ref_pkg.prefix)
                else:
                    clade_exclusion_json = intermediates_path + ts_evaluate.ref_pkg.prefix + ts_evaluate.ref_pkg.refpkg_suffix
                    ce_refpkg = ts_evaluate.ref_pkg.clone(clade_exclusion_json)
                    classification_table = classifier_output + "final_outputs" + os.sep + "marker_contig_map.tsv"

                    if not os.path.isfile(classification_table):
                        # Copy reference files, then exclude all clades belonging to the taxon being tested

                        ce_refpkg.exclude_clade_from_ref_files(intermediates_path, lineage,
                                                               ts_evaluate.executables, args.fresh)
                        # Write the query sequences
                        fasta.write_new_fasta(taxon_rep_seqs, test_rep_taxa_fasta)
                        assign_args = ["-i", test_rep_taxa_fasta, "-o", classifier_output,
                                       "--refpkg_dir", os.path.dirname(ce_refpkg.f__json),
                                       "-m", ts_evaluate.molecule_type, "-n", str(args.num_threads),
                                       "--min_seq_length", str(min_seq_length),
                                       "--overwrite", "--delete", "--no_svm"]
                        if args.trim_align:
                            assign_args.append("--trim_align")
                        try:
                            assign(assign_args)
                        except:  # Just in case treesapp assign fails, just continue
                            pass

                        if not os.path.isfile(classification_table):
                            # The TaxonTest object is maintained for record-keeping (to track # queries & classifieds)
                            logging.warning("TreeSAPP did not generate output for '{}'. Skipping.\n".format(lineage))
                            shutil.rmtree(classifier_output)
                            continue
                    else:
                        # Valid number of queries and these sequences have already been classified
                        pass

                    test_obj.taxonomic_tree = ce_refpkg.all_possible_assignments()
                    if os.path.isfile(classification_table):
                        assigned_lines = file_parsers.read_marker_classification_table(classification_table)
                        test_obj.assignments = file_parsers.parse_assignments(assigned_lines)
                        test_obj.filter_assignments(ts_evaluate.ref_pkg.prefix)
                        test_obj.distances = parse_distances(assigned_lines)
                    else:
                        logging.error("marker_contig_map.tsv is missing from output directory '" +
                                      os.path.dirname(classification_table) + "'\n" +
                                      "Please remove this directory and re-run.\n")
                        sys.exit(21)

    if ts_evaluate.stage_status("calculate"):
        # everything has been prepared, only need to parse the classifications and map lineages
        logging.info("Finishing up the mapping of classified, filtered taxonomic sequences.\n")
        for rank in sorted(ts_evaluate.taxa_tests):
            for test_obj in ts_evaluate.taxa_tests[rank]:  # type: classy.TaxonTest
                if test_obj.assignments:
                    marker_assignments = map_headers_to_lineage(test_obj.assignments, fasta_records)
                    # Return the number of correct, classified, and total sequences of that taxon at the current rank
                    # Identify the excluded rank for each query sequence
                    if len(marker_assignments) == 0:
                        logging.debug("No '{}' sequences were classified.\n".format(test_obj.taxon))
                        continue

                    for marker in marker_assignments:
                        ts_evaluate.markers.add(marker)
                    rank_assignments = lca_calculations.identify_excluded_clade(marker_assignments,
                                                                                test_obj.taxonomic_tree,
                                                                                ts_evaluate.ref_pkg.prefix)
                    for a_rank in rank_assignments:
                        if a_rank != rank and len(rank_assignments[a_rank]) > 0:
                            logging.warning("{}-level clade excluded but optimal classification found to be {}-level.\n"
                                            "Assignments were: {}\n".format(rank, a_rank, rank_assignments[a_rank]))
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

        if ts_evaluate.ref_pkg.prefix not in ts_evaluate.markers:
            logging.error("No sequences were classified as {}.\n".format(ts_evaluate.ref_pkg.prefix))
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

    if args.delete and os.path.isdir(ts_evaluate.var_output_dir):
        shutil.rmtree(ts_evaluate.var_output_dir)

    return
