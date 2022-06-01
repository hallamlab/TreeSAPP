import logging
import sys
import re
import os
import shutil

import tqdm
from joblib import dump as jdump
from joblib import load as jload

from treesapp import seq_clustering
from treesapp import entrez_utils
from treesapp import file_parsers
from treesapp import fasta
from treesapp import utilities
from treesapp import wrapper
from treesapp import entish
from treesapp import lca_calculations
from treesapp import placement_trainer
from treesapp import annotate_extra
from treesapp import treesapp_args
from treesapp import classy
from treesapp import logger
from treesapp import jplace_utils
from treesapp import training_utils
from treesapp import phylogeny_painting as paint
from treesapp import refpkg as ts_ref_pkg
from treesapp import phylo_seq as ts_phylo_seq
from treesapp import clade_exclusion_evaluator as ts_clade_ex
from treesapp import assign as ts_assign_mod
from treesapp import create_refpkg as ts_create_mod
from treesapp import update_refpkg as ts_update_mod
from treesapp import hmmer_tbl_parser

LOGGER = logging.getLogger(logger.logger_name())


def info(sys_args):
    """

    """
    parser = treesapp_args.TreeSAPPArgumentParser(description="Return package, executable and refpkg information.")
    treesapp_args.add_info_arguments(parser)
    args = parser.parse_args(sys_args)
    logger.prep_logging()
    ts_info = classy.TreeSAPP("info")

    import treesapp
    import Bio
    import numpy
    import packaging
    import pygtrie
    import scipy
    import ete3
    import sklearn
    import joblib
    import seaborn
    import samsum
    import pyfastx
    import tqdm
    LOGGER.info("TreeSAPP version " + treesapp.__version__ + ".\n")

    # Write the version of all python deps
    py_deps = {"biopython": Bio.__version__,
               "ete3": ete3.__version__,
               "joblib": joblib.__version__,
               "numpy": numpy.__version__,
               "packaging": packaging.__version__,
               "pyfastx": pyfastx.version(),
               "pygtrie": pygtrie.__version__,
               "samsum": samsum.__version__,
               "scikit-learn": sklearn.__version__,
               "scipy": scipy.__version__,
               "seaborn": seaborn.__version__,
               "tqdm": tqdm.__version__}

    LOGGER.info("Python package dependency versions:\n\t" +
                "\n\t".join([k + ": " + v for k, v in py_deps.items()]) + "\n")

    # Write the version of executable deps
    ts_info.furnish_with_arguments(args)
    if args.refpkg_dir:
        ts_info.refpkg_dir = args.refpkg_dir
    LOGGER.info(utilities.executable_dependency_versions(ts_info.executables))

    if args.verbose:
        refpkg_dict = ts_ref_pkg.gather_ref_packages(ts_info.refpkg_dir)
        refpkg_summary_str = "\t".join(["Name",
                                        "Molecule", "Tree builder", "RefPkg-type", "Leaf nodes",
                                        "Description", "Created", "Last-updated"]) + "\n"
        for refpkg_name in sorted(refpkg_dict, key=lambda x: refpkg_dict[x].prefix):
            refpkg = refpkg_dict[refpkg_name]  # type: ts_ref_pkg.ReferencePackage
            refpkg_summary_str += "\t".join([refpkg_name,
                                             refpkg.molecule, refpkg.tree_tool, refpkg.kind, str(refpkg.num_seqs),
                                             refpkg.description, refpkg.date, refpkg.update]
                                            ) + "\n"
        LOGGER.info(refpkg_summary_str)

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
rename      Change the name (prefix attribute) of a reference package
**
Use '-h' to get subcommand-specific help, e.g. 'treesapp package view -h'
"""
    parser = treesapp_args.TreeSAPPArgumentParser(description='Facilitate operations on reference packages')
    parser.add_argument("subcommand", nargs='?', choices=["view", "edit", "rename"],
                        help="A subcommand specifying the type of operation to perform. "
                             "`treesapp package rename` must be followed only by 'prefix' and the new value. "
                             "Example: `treesapp package rename prefix Xyz -r path/to/Xyz_build.pkl`")
    args = parser.parse_args(sys_args[0:1])
    if not args.subcommand:
        sys.stderr.write(pkg_usage)
        sys.exit(1)

    ref_pkg = ts_ref_pkg.ReferencePackage()

    treesapp_args.add_package_arguments(parser, ref_pkg.get_public_attributes())
    args = parser.parse_args(sys_args)

    if args.output:
        log_dir = args.output
        if not os.path.isdir(args.output):
            os.mkdir(args.output)
    else:
        log_dir = "./"

    logger.prep_logging(log_file=os.path.join(log_dir, 'TreeSAPP_package_log.txt'))

    for refpkg_pkl in args.pkg_path:
        ref_pkg.f__pkl = refpkg_pkl
        ref_pkg.slurp()

        if args.output:
            output_dir = args.output
        else:
            output_dir = os.path.dirname(ref_pkg.f__pkl)

        if args.subcommand == "view":
            ts_ref_pkg.view(ref_pkg, args.attributes)
        elif args.subcommand == "edit":
            ts_ref_pkg.edit(ref_pkg, args.attributes, output_dir,
                            overwrite=args.overwrite, phenotypes=args.phenotypes, reset=args.reset, join=args.join)
        elif args.subcommand == "rename":
            ts_ref_pkg.rename(ref_pkg, args.attributes, output_dir, args.overwrite)
        else:
            LOGGER.error("Unrecognized command: '{}'.\n{}\n".format(args.subcommand, pkg_usage))
            sys.exit(1)

    return


def train(sys_args):
    parser = treesapp_args.TreeSAPPArgumentParser(description='Model evolutionary distances across taxonomic ranks.')
    treesapp_args.add_trainer_arguments(parser)
    args = parser.parse_args(sys_args)

    ts_trainer = training_utils.PhyTrainer()
    ts_trainer.furnish_with_arguments(args)
    ts_trainer.check_previous_output(args.overwrite)

    log_file_name = args.output + os.sep + "TreeSAPP_trainer_log.txt"
    logger.prep_logging(log_file_name, args.verbose)
    LOGGER.info("\n##\t\t\tTrain taxonomic rank-placement distance model\t\t\t##\n")

    treesapp_args.check_parser_arguments(args, sys_args)
    ts_trainer.check_trainer_arguments(args)
    ts_trainer.set_file_paths()
    ts_trainer.decide_stage(args)
    ts_trainer.ref_pkg.disband(os.path.join(ts_trainer.var_output_dir, ts_trainer.ref_pkg.prefix + "_RefPkg"))

    if args.annot_map and args.classifier == "bin":
        ts_trainer.pkg_dbname_dict = file_parsers.read_annotation_mapping_file(args.annot_map)

    train_seqs = fasta.FASTA(ts_trainer.input_sequences)
    ##
    # STAGE 1: Optionally validate and reformat the input FASTA
    ##
    if ts_trainer.stage_status("clean"):
        LOGGER.info("Reading and formatting {}... ".format(ts_trainer.input_sequences))
        train_seqs.header_registry = fasta.register_headers(fasta.format_fasta(fasta_input=ts_trainer.input_sequences,
                                                                               molecule="prot",
                                                                               output_fasta=ts_trainer.formatted_input))
        LOGGER.info("done.\n")
        ts_trainer.increment_stage_dir()
    else:
        ts_trainer.formatted_input = ts_trainer.input_sequences
        train_seqs.header_registry = fasta.register_headers(fasta.get_headers(ts_trainer.formatted_input), True)
    train_seqs.file = ts_trainer.formatted_input

    ##
    # STAGE 2: Run hmmsearch on the query sequences to search for reference package homologs
    ##
    if ts_trainer.stage_status("search"):
        LOGGER.info("Searching for homologous sequences with hmmsearch... ")
        hmm_domtbl_file = wrapper.run_hmmsearch(ts_trainer.executables["hmmsearch"],
                                                ts_trainer.ref_pkg.f__search_profile,
                                                ts_trainer.formatted_input,
                                                ts_trainer.stage_output_dir)
        LOGGER.info("done.\n")
        thresholds = hmmer_tbl_parser.prep_args_for_parsing(args)
        hmm_matches = file_parsers.parse_domain_tables(thresholds, {(ts_trainer.ref_pkg.prefix, "ORFs"):
                                                                        hmm_domtbl_file})
        ts_assign_mod.load_homologs(hmm_matches, ts_trainer.formatted_input, train_seqs)

        LOGGER.info(train_seqs.summarize_fasta_sequences())
        fasta.write_new_fasta(train_seqs.fasta_dict, ts_trainer.hmm_purified_seqs)
        train_seqs.file = ts_trainer.hmm_purified_seqs
        utilities.hmm_pile(hmm_matches)
        ts_trainer.increment_stage_dir()
    else:
        if not os.path.isfile(ts_trainer.hmm_purified_seqs):
            ts_trainer.hmm_purified_seqs = ts_trainer.input_sequences
            train_seqs.file = ts_trainer.hmm_purified_seqs
            train_seqs.fasta_dict = fasta.read_fasta_to_dict(train_seqs.file)
            if not train_seqs.fasta_dict:
                LOGGER.error("No sequences were detected in the HMM-purified FASTA '{}'.\n".format(train_seqs.file))
                sys.exit(13)
        else:
            train_seqs.file = ts_trainer.hmm_purified_seqs
            train_seqs.load_fasta()
        LOGGER.info("Profile HMM homology search skipped. Using all sequences in {}.\n".format(train_seqs.file))
        train_seqs.change_dict_keys("original")

    training_utils.summarize_query_classes(set(ts_trainer.pkg_dbname_dict.keys()), set(train_seqs.fasta_dict.keys()))

    ##
    # STAGE 3: Download the taxonomic lineages for each query sequence
    ##
    entrez_record_dict = ts_trainer.fetch_entrez_lineages(train_seqs,
                                                          molecule=ts_trainer.ref_pkg.molecule,
                                                          acc_to_taxid=args.acc_to_taxid,
                                                          seqs_to_lineage=args.seq_names_to_taxa)

    entrez_utils.fill_ref_seq_lineages(entrez_record_dict, ts_trainer.seq_lineage_map, complete=True)
    ref_leaf_nodes = ts_phylo_seq.convert_entrez_to_tree_leaf_references(entrez_record_dict)
    ts_trainer.ref_pkg.taxa_trie.feed_leaf_nodes(ref_leaf_nodes)
    entrez_utils.sync_record_and_hierarchy_lineages(ref_leaf_nodes, entrez_record_dict)
    ts_trainer.ref_pkg.taxa_trie.validate_rank_prefixes()
    ts_trainer.ref_pkg.taxa_trie.build_multifurcating_trie()
    rank_depth_map = ts_trainer.ref_pkg.taxa_trie.accepted_ranks_depths
    taxa_evo_dists = dict()

    query_seq_records = {}
    for index, e_record in entrez_record_dict.items():  # type: entrez_utils.EntrezRecord
        query_seq_records[e_record.description] = e_record

    # Goal is to use the distances already calculated but re-print
    if ts_trainer.stage_status("place"):
        clade_ex_pqueries = dict()
        if os.path.isfile(ts_trainer.placement_summary):
            # Read the summary file and pull the phylogenetic distances for each rank
            taxa_evo_dists = placement_trainer.read_placement_summary(ts_trainer.placement_summary)
            # Remove any ranks that are not to be used in this estimation
            estimated_ranks = set(taxa_evo_dists.keys())
            for rank_key in estimated_ranks.difference(set(ts_trainer.training_ranks.keys())):
                taxa_evo_dists.pop(rank_key)

        if len(set(ts_trainer.training_ranks).difference(set(taxa_evo_dists))) > 0 or len(clade_ex_pqueries) == 0:
            clade_ex_pqueries = placement_trainer.gen_cladex_data(train_seqs.file, ts_trainer.executables,
                                                                  ts_trainer.ref_pkg, ts_trainer.seq_lineage_map,
                                                                  ts_trainer.stage_output_dir,
                                                                  max_examples=args.max_examples,
                                                                  num_threads=args.num_threads)
            if not clade_ex_pqueries:
                return
            jdump(value=clade_ex_pqueries, filename=os.path.join(ts_trainer.clade_ex_pquery_pkl))

        # Write the tab-delimited file with metadata included for each placement
        placement_trainer.write_placement_table(clade_ex_pqueries,
                                                ts_trainer.placement_table, ts_trainer.ref_pkg.prefix)

        LOGGER.info("Generating placement data without clade exclusion for SVM... ")
        LOGGER.disabled = True

        assign_prefix = os.path.join(ts_trainer.stage_output_dir, ts_trainer.ref_pkg.prefix + "_assign")
        assign_params = ["-i", train_seqs.file,
                         "-o", assign_prefix,
                         "--num_procs", str(args.num_threads),
                         "--refpkg_dir", os.path.dirname(ts_trainer.ref_pkg.f__pkl),
                         "--targets", ts_trainer.ref_pkg.prefix,
                         "--molecule", ts_trainer.ref_pkg.molecule,
                         "--delete"]
        if args.trim_align:
            assign_params.append("--trim_align")

        try:
            ts_assign_mod.assign(assign_params)
            plain_pqueries = ts_phylo_seq.assignments_to_pqueries(file_parsers.read_classification_table(
                os.path.join(assign_prefix, "final_outputs", "classifications.tsv")))
        except (SystemExit, IOError):
            LOGGER.info("failed.\n")
            LOGGER.error("treesapp assign did not complete successfully.\n")
            sys.exit(11)
        jdump(value=plain_pqueries, filename=os.path.join(ts_trainer.plain_pquery_pkl))

        # Re-enable LOGGER at the previous level
        LOGGER.disabled = False
        LOGGER.info("done.\n")
        ts_trainer.increment_stage_dir()
    else:
        LOGGER.info("Phylogenetic placement stage is being skipped. Reading saved pickles... ")
        # Load saved pquery instances from a previous run
        clade_ex_pqueries = jload(filename=ts_trainer.clade_ex_pquery_pkl)
        plain_pqueries = jload(filename=ts_trainer.plain_pquery_pkl)
        LOGGER.info("done.\n")

    taxa_evo_dists = placement_trainer.evo_dists_from_pqueries(clade_ex_pqueries, ts_trainer.training_ranks)
    ts_trainer.pqueries.update(clade_ex_pqueries)
    ts_trainer.pqueries.update(plain_pqueries)

    if ts_trainer.stage_status("train"):
        # Finish up the linear regression model
        ts_trainer.ref_pkg.pfit = placement_trainer.complete_regression(taxa_evo_dists, ts_trainer.training_ranks)
        if ts_trainer.ref_pkg.pfit:
            LOGGER.info("Placement distance model complete.\n")
        else:
            LOGGER.info("Unable to complete phylogenetic distance and rank correlation.\n")

        # Reformat the dictionary containing PQuery instances for classifier training and testing
        refpkg_pqueries = placement_trainer.flatten_pquery_dict(ts_trainer.pqueries, ts_trainer.ref_pkg.prefix)
        tp_names = {}
        fp_names = {}
        if args.classifier == "occ":
            # The individual instances are of PQuery type
            tp_names.update({ts_trainer.ref_pkg.prefix:
                                 [pquery.seq_name for pquery in refpkg_pqueries[ts_trainer.ref_pkg.prefix]]})

            for pquery in refpkg_pqueries[ts_trainer.ref_pkg.prefix]:
                if not pquery.rank:
                    pquery.rank = "species"
        else:
            tp, fp, fn = training_utils.bin_headers(refpkg_pqueries, ts_trainer.pkg_dbname_dict, query_seq_records,
                                                    {ts_trainer.ref_pkg.prefix: ts_trainer.ref_pkg})
            # The individual instances are of QuerySequence type
            for refpkg_name, tp_query_seqs in tp.items():
                tp_names[refpkg_name] = [qseq.seq_name for qseq in tp_query_seqs]
            for refpkg_name, fp_query_seqs in fp.items():
                fp_names[refpkg_name] = [qseq.seq_name for qseq in fp_query_seqs]

        LOGGER.info("Extracting features from TreeSAPP classifications... ")
        training_df = training_utils.load_training_data_frame(pqueries=refpkg_pqueries,
                                                              refpkg_map={ts_trainer.ref_pkg.prefix:
                                                                              ts_trainer.ref_pkg},
                                                              refpkg_positive_annots=tp_names)
        LOGGER.info("done.\n")
        training_utils.train_classifier_from_dataframe(training_df, args.kernel,
                                                       grid_search=args.grid_search,
                                                       num_threads=args.num_threads,
                                                       tsne=ts_trainer.tsne_plot)
        ts_trainer.increment_stage_dir()
    else:
        training_df = training_utils.load_training_data_frame({}, {}, {})

    if ts_trainer.stage_status("update"):
        if not ts_trainer.ref_pkg.pfit:
            LOGGER.warning("Linear regression parameters could not be estimated. " +
                           "Taxonomic ranks will not be distance-adjusted during classification for this package.\n")
            ts_trainer.ref_pkg.pfit = [0.0, 7.0]

        ts_trainer.ref_pkg.training_df = training_df
        ts_trainer.ref_pkg.f__pkl = os.path.join(ts_trainer.final_output_dir,
                                                 os.path.basename(ts_trainer.ref_pkg.f__pkl))

        ts_trainer.ref_pkg.validate()
        ts_trainer.ref_pkg.pickle_package()

    # Write the text file containing distances used in the regression analysis
    with open(ts_trainer.placement_summary, 'w') as out_handler:
        trained_string = "Regression parameters = " + re.sub(' ', '', str(ts_trainer.ref_pkg.pfit)) + "\n"
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
    logger.prep_logging(log_file_name, args.verbose)
    LOGGER.info("\n##\t\t\tCreating TreeSAPP reference package\t\t\t##\n")

    treesapp_args.check_parser_arguments(args, sys_args)
    ts_create_mod.check_create_arguments(ts_create, args)
    ts_create.decide_stage(args)

    # Populate the final TreeSAPP reference file paths with the proper output directory
    ts_create.ref_pkg.update_file_names()
    ts_create.ref_pkg.change_file_paths(ts_create.final_output_dir)
    ts_create.ref_pkg.cmd = ' '.join(["treesapp", "create"] + sys_args)

    ref_seqs = fasta.FASTA(ts_create.input_sequences)
    # Read the FASTA into a dictionary - homologous sequences will be extracted from this
    ref_seqs.load_fasta(format_it=True, molecule=ts_create.ref_pkg.molecule)
    if ts_create.stage_status("deduplicate"):
        seq_clustering.dereplicate_by_clustering(fasta_inst=ref_seqs,
                                                 prop_similarity=0.999,
                                                 mmseqs_exe=ts_create.executables["mmseqs"],
                                                 tmp_dir=ts_create.stage_output_dir,
                                                 num_threads=args.num_threads)
        ts_create.input_sequences = ts_create.stage_output_dir + "deduplicated.fasta"
        fasta.write_new_fasta(fasta_dict=ref_seqs.fasta_dict, fasta_name=ts_create.input_sequences)

    if ts_create.stage_status("search"):
        profile_match_dict = dict()
        LOGGER.debug("Raw, unfiltered sequence summary:\n" + ref_seqs.summarize_fasta_sequences())

        LOGGER.info("Searching for domain sequences... ")
        hmm_domtbl_file = wrapper.run_hmmsearch(ts_create.executables["hmmsearch"], ts_create.hmm_profile,
                                                ts_create.input_sequences, ts_create.var_output_dir, args.num_threads)
        LOGGER.info("done.\n")
        thresholds = hmmer_tbl_parser.prep_args_for_parsing(args)
        hmm_matches = file_parsers.parse_domain_tables(thresholds, {(ts_create.ref_pkg.prefix, "ORFs"):
                                                                        hmm_domtbl_file})
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
    ref_seqs.add_accession_to_headers()
    # Get rid of ambiguity or unusual characters
    ref_seqs.replace_ambiguity_chars(ts_create.ref_pkg.molecule)
    LOGGER.info("Sequence summary:\n" + ref_seqs.summarize_fasta_sequences())

    ##
    # If there are sequences that needs to be guaranteed to be included,
    #  add them now as its easier to work with more sequences than repeat everything
    ##
    if args.guarantee:
        ref_seqs.fasta_update(args.guarantee)
        ref_seqs.change_dict_keys("original")

    ##
    # Save all sequence names in the header registry as EntrezRecord instances
    # Using the accession-lineage map (if available), map the sequence names to their respective lineages
    # Proceed with creating the Entrez-queries for sequences lacking lineage information
    ##
    fasta_records = ts_create.fetch_entrez_lineages(ref_seqs, ts_create.ref_pkg.molecule,
                                                    args.acc_to_taxid, args.seq_names_to_taxa)
    entrez_utils.fill_ref_seq_lineages(fasta_records, ts_create.seq_lineage_map)
    ref_leaf_nodes = ts_phylo_seq.convert_entrez_to_tree_leaf_references(fasta_records)
    ts_create.ref_pkg.taxa_trie.feed_leaf_nodes(ref_leaf_nodes)
    entrez_utils.sync_record_and_hierarchy_lineages(ref_leaf_nodes, fasta_records)
    entrez_utils.repair_lineages(fasta_records, ts_create.ref_pkg.taxa_trie)
    entrez_utils.fill_entrez_record_taxon_rank(fasta_records, ts_create.ref_pkg.taxa_trie)
    ts_create.ref_pkg.taxa_trie.validate_rank_prefixes()
    ts_create.ref_pkg.taxa_trie.build_multifurcating_trie()
    ts_create_mod.strip_rank_prefix_from_organisms(fasta_records, ts_create.ref_pkg.taxa_trie)

    prefilter_ref_seqs = entrez_utils.entrez_record_snapshot(fasta_records)

    if ts_create.stage_status("clean"):
        # Remove the sequences failing 'filter' and/or only retain the sequences in 'screen'
        fasta_records = ts_create_mod.screen_filter_taxa(fasta_records, args.screen, args.filter, ref_seqs.amendments)
        # Remove the sequence records with low resolution lineages, according to args.min_taxonomic_rank
        fasta_records = ts_create_mod.remove_by_truncated_lineages(fasta_records, args.min_taxonomic_rank,
                                                                   ts_create.ref_pkg.taxa_trie, ref_seqs.amendments)

        if len(fasta_records.keys()) < 2:
            LOGGER.error("{} sequences post-homology + taxonomy filtering\n".format(len(fasta_records)))
            sys.exit(11)
        # Write a new FASTA file containing the sequences that passed the homology and taxonomy filters

    ref_seqs.file = ts_create.filtered_fasta
    # NOTE: original header must be used as this is being passed to train
    ref_seqs.unalign()
    ref_seqs.change_dict_keys("original")
    filtered_headers = [ref_seqs.header_registry[num_id].original for num_id in fasta_records]
    # ref_seqs.keep_only(filtered_headers)  # Currently avoiding this as it causes a KeyError for guaranteed seqs
    fasta.write_new_fasta(fasta_dict=ref_seqs.fasta_dict, fasta_name=ref_seqs.file, headers=filtered_headers)

    ##
    # Optionally cluster the input sequences using MMSeqs' linclust at the specified identity
    ##
    if ts_create.stage_status("cluster"):
        # TODO: replace with fasta.dereplicate_by_clustering()
        pre_cluster = ref_seqs.n_seqs()
        ref_seqs.change_dict_keys("num_id")
        # Write a FASTA for clustering containing the formatted headers since
        # not all clustering tools + versions keep whole header - spaces are replaced with underscores
        fasta.write_new_fasta(fasta_dict=ref_seqs.fasta_dict,
                              fasta_name=ts_create.cluster_input,
                              headers=list(fasta_records.keys()))
        if args.cluster:
            wrapper.cluster_sequences(software_path=ts_create.executables["mmseqs"],
                                      fasta_input=ts_create.cluster_input, output_prefix=ts_create.clusters_prefix,
                                      similarity=ts_create.ref_pkg.pid, num_threads=args.num_threads)
            ts_create.clusters_table = ts_create.clusters_prefix + "_cluster.tsv"
            cluster_alignments = ts_create.clusters_prefix + "_cluster_aln.tsv"

            cluster_dict = seq_clustering.create_mmseqs_clusters(ts_create.clusters_table, cluster_alignments)

            # Revert headers in cluster_dict from 'formatted' back to 'original'
            fasta.rename_cluster_headers(cluster_dict, ref_seqs.header_registry)
            LOGGER.debug("\t{} sequence clusters\n".format(len(cluster_dict.keys())))
            ##
            # Calculate LCA of each cluster to represent the taxonomy of the representative sequence
            ##
            ts_create_mod.find_cluster_lca(cluster_dict, fasta_records, ref_seqs.header_registry)
        else:
            cluster_dict = None

        ##
        # Swap sequences in 'guarantee' for the representatives, creating new clusters
        ##
        if args.guarantee and ts_create.clusters_table:
            # We don't want to make the tree redundant so instead of simply adding the sequences in guarantee,
            #  we will swap them for their respective representative sequences.
            # All important sequences become representative, even if multiple are in the same cluster
            very_important_seqs = set([ref_seqs.header_registry[num].original for num in ref_seqs.amendments])
            cluster_dict = ts_create_mod.guarantee_ref_seqs(cluster_dict, very_important_seqs)

        ##
        # Set the cluster-specific values for ReferenceSequence objects
        ##
        if ts_create.clusters_table and not args.headless:
            # Allow user to select the representative sequence based on organism name, sequence length and similarity
            ts_create_mod.present_cluster_rep_options(cluster_dict, fasta_records, ref_seqs.header_registry,
                                                      important_seqs=ref_seqs.amendments)
        elif ts_create.clusters_table and args.headless:
            ts_create_mod.finalize_cluster_reps(cluster_dict, fasta_records, ref_seqs.header_registry)
        else:
            pass
        ref_seqs.keep_only(header_subset=[x for x in fasta_records if fasta_records[x].cluster_rep])
        post_cluster = ref_seqs.n_seqs()
        ts_create.overcluster_warning(pre_cluster, post_cluster)

    if ts_create.stage_status("build"):
        if args.od_seq:
            ts_create_mod.remove_outlier_sequences(fasta_records,
                                                   ts_create.executables["OD-seq"], ts_create.executables["mafft"],
                                                   ts_create.var_output_dir, args.num_threads)

        # This precautionary measure is for `create` called from `update` and reference seqs have the assign signature
        accession_ids = [fasta_records[num_id].accession for num_id in fasta_records]
        name_map = ts_update_mod.strip_assigment_pattern(accession_ids, ts_create.ref_pkg.prefix)
        for num_id in fasta_records:
            record = fasta_records[num_id]
            record.accession = name_map[record.accession]
        ##
        # Re-order the fasta_records by their lineages (not phylogenetic, just alphabetical sort)
        # Remove the cluster members since they will no longer be used
        ##
        fasta_replace_dict = ts_create_mod.order_dict_by_lineage(fasta_records)

        ts_create.ref_pkg.lineage_ids = ts_create_mod.lineages_to_dict(fasta_replace_dict, args.taxa_lca)

        postfilter_ref_seqs = entrez_utils.entrez_record_snapshot(fasta_records)
        filtered_ref_seqs = utilities.dict_diff(prefilter_ref_seqs, postfilter_ref_seqs)
        LOGGER.debug("{0} references before and {1} remaining after filtering.\n".format(len(prefilter_ref_seqs),
                                                                                         len(postfilter_ref_seqs)))
        ts_create.ref_pkg.taxa_trie.jetison_taxa_from_hierarchy(filtered_ref_seqs)

        taxonomic_summary = ts_create_mod.summarize_reference_taxa(fasta_replace_dict, ts_create.ref_pkg.taxa_trie,
                                                                   args.taxa_lca)
        LOGGER.info(taxonomic_summary)

        ##
        # Perform multiple sequence alignment
        ##
        if args.multiple_alignment:
            ts_create_mod.create_new_ref_fasta(ts_create.unaln_ref_fasta, fasta_replace_dict, True)
        else:
            ts_create_mod.create_new_ref_fasta(ts_create.unaln_ref_fasta, fasta_replace_dict)

        if args.multiple_alignment is False:
            LOGGER.info("Aligning the sequences using MAFFT... ")
            wrapper.run_mafft(ts_create.executables["mafft"], ts_create.unaln_ref_fasta,
                              ts_create.ref_pkg.f__msa, args.num_threads)
            LOGGER.info("done.\n")
        else:
            pass
        ref_seqs.file = ts_create.ref_pkg.f__msa
        ref_seqs.load_fasta()
        ts_create.ref_pkg.num_seqs = ref_seqs.n_seqs()
        n_rows, n_cols = fasta.multiple_alignment_dimensions(mfa_file=ts_create.ref_pkg.f__msa,
                                                             seq_dict=ref_seqs.fasta_dict)
        LOGGER.debug("Reference alignment contains {} sequences with {} character positions.\n".format(n_rows, n_cols))

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

        ts_create.ref_pkg.band()
        ts_create.ref_pkg.hmm_length()
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
            LOGGER.debug("Number of sequences discarded: " + summary_str + "\n")
            if len(qc_ma_dict.keys()) == 0:
                # At least one of the reference sequences were discarded and therefore this package is invalid.
                LOGGER.error("Trimming removed reference sequences. This could indicate non-homologous sequences.\n" +
                             "Please improve sequence quality-control and/or rerun without the '--trim_align' flag.\n")
                sys.exit(13)
            elif len(qc_ma_dict.keys()) > 1:
                LOGGER.error("Multiple trimmed alignment files are found when only one is expected:\n" +
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
        ts_create.determine_model(ts_create.ref_pkg)
        best_tree = ts_create.ref_pkg.infer_phylogeny(ts_create.phylip_file, ts_create.executables, ts_create.phy_dir,
                                                      args.num_threads)
    else:
        best_tree = "{}{}.{}.nwk".format(ts_create.phy_dir, ts_create.ref_pkg.prefix, ts_create.ref_pkg.tree_tool)

    if ts_create.stage_status("evaluate"):
        # Evaluate the model parameters with RAxML-NG. Output is required by EPA-NG
        wrapper.model_parameters(ts_create.executables["raxml-ng"],
                                 ts_create.phylip_file, best_tree, ts_create.phy_dir + ts_create.ref_pkg.prefix,
                                 ts_create.ref_pkg.sub_model, args.num_threads)
    ts_create.ref_pkg.recover_raxmlng_model_outputs(ts_create.phy_dir)
    ts_create.ref_pkg.recover_raxmlng_tree_outputs(ts_create.phy_dir)

    if ts_create.stage_status("support"):
        # Perform non-parametric bootstrapping with RAxML-NG and calculate branch support values from bootstraps
        wrapper.support_tree_raxml(raxml_exe=ts_create.executables["raxml-ng"], model=ts_create.ref_pkg.sub_model,
                                   ref_tree=best_tree, ref_msa=ts_create.phylip_file,
                                   tree_prefix=os.path.join(ts_create.phy_dir, ts_create.ref_pkg.prefix),
                                   mre=True, n_bootstraps=args.bootstraps, num_threads=args.num_threads)
        ts_create.ref_pkg.recover_raxmlng_supports(ts_create.phy_dir)

    ts_create.ref_pkg.band()
    # Build the regression model of placement distances to taxonomic ranks
    trainer_cmd = ts_create_mod.formulate_train_command(input_seqs=ts_create.filtered_fasta,
                                                        ref_pkg=ts_create.ref_pkg,
                                                        output_dir=ts_create.training_dir,
                                                        acc_to_lin=ts_create.acc_to_lin,
                                                        args=args)

    if ts_create.stage_status("train"):
        train(trainer_cmd)
        trained_refpkg = os.path.join(ts_create.var_output_dir, "placement_trainer", "final_outputs",
                                      ts_create.ref_pkg.prefix + ts_create.ref_pkg.refpkg_suffix)
    else:
        LOGGER.info("Skipping training:\n$ treesapp train {}.\n".format(' '.join(trainer_cmd)))
        trained_refpkg = ''

    ##
    # Finish validating the file and append the reference package build parameters to the master table
    ##
    if ts_create.stage_status("update"):
        if os.path.isfile(trained_refpkg):
            ts_create.ref_pkg.f__pkl = trained_refpkg
        ts_create.ref_pkg.slurp()
        ts_create.ref_pkg.validate()
        ts_create.ref_pkg.change_file_paths(ts_create.final_output_dir)
        ts_create.ref_pkg.pickle_package()

    ts_create.remove_intermediates(args.delete)
    ts_create.print_terminal_commands()

    return


def update(sys_args):
    parser = treesapp_args.TreeSAPPArgumentParser(description='Update a reference package with assigned sequences.')
    treesapp_args.add_update_arguments(parser)
    args = parser.parse_args(sys_args)

    ts_updater = ts_update_mod.Updater()
    ts_updater.furnish_with_arguments(args)
    ts_updater.check_previous_output(args.overwrite)

    log_file_name = args.output + os.sep + "TreeSAPP_update_log.txt"
    logger.prep_logging(log_file_name, args.verbose)
    LOGGER.info("\n##\t\t\tUpdating TreeSAPP reference package\t\t\t##\n")

    treesapp_args.check_parser_arguments(args, sys_args)
    treesapp_args.check_updater_arguments(ts_updater, args)
    ts_updater.decide_stage(args)
    ts_updater.ref_pkg.validate()
    ref_seq_lineage_info = ts_updater.ref_pkg.generate_tree_leaf_references_from_refpkg()

    ##
    # Pull out sequences from TreeSAPP output
    ##
    classified_fasta = fasta.FASTA(ts_updater.query_sequences)  # These are the classified sequences
    classified_fasta.load_fasta()
    classified_fasta.add_accession_to_headers(ts_updater.ref_pkg.prefix)
    pqueries = file_parsers.load_classified_sequences_from_assign_output(assign_output_dir=ts_updater.treesapp_output,
                                                                         assigner_cls=ts_assign_mod.Assigner(),
                                                                         refpkg_name=ts_updater.ref_pkg.prefix)
    candidate_update_pqueries = ts_update_mod.filter_by_placement_thresholds(pqueries,
                                                                             args.min_lwr, args.max_pd, args.max_evo)
    classified_targets = [pq.place_name for pq in candidate_update_pqueries]
    # Remove classified sequences that are already in the reference package
    ts_update_mod.drop_queries_by_accession(classified_targets, ref_seq_lineage_info)

    if len(classified_targets) == 0:
        LOGGER.error("No new candidate reference sequences. Skipping update.\n")
        return
    classified_fasta.change_dict_keys("original")

    ##
    # Filter out sequences that shouldn't be used in the update: different refpkg, too short, low LWR, etc.
    ##
    classified_fasta.keep_only(classified_targets)
    LOGGER.info(classified_fasta.summarize_fasta_sequences())
    ts_updater.min_length = ts_update_mod.decide_length_filter(ts_updater.ref_pkg, args.min_seq_length)
    classified_fasta.remove_shorter_than(ts_updater.min_length)
    if classified_fasta.n_seqs() == 0:
        LOGGER.error("No classified sequences exceed minimum length threshold of {}.\n".format(ts_updater.min_length))
        return

    ##
    # Add lineages - use taxa if provided with a table mapping contigs to taxa, TreeSAPP-assigned taxonomy otherwise
    ##
    classified_seq_lineage_map = dict()
    querying_classified_fasta = fasta.FASTA("")
    querying_classified_fasta.clone(classified_fasta)

    if args.skip_assign:
        name_map = ts_update_mod.strip_assigment_pattern(querying_classified_fasta.get_seq_names(),
                                                         ts_updater.ref_pkg.prefix)
        querying_classified_fasta.synchronize_seqs_n_headers()
        querying_classified_fasta.swap_headers(name_map)
        fasta_records = ts_updater.fetch_entrez_lineages(ref_seqs=querying_classified_fasta,
                                                         molecule=ts_updater.ref_pkg.molecule,
                                                         seqs_to_lineage=ts_updater.seq_names_to_taxa)
        entrez_utils.fill_ref_seq_lineages(fasta_records, classified_seq_lineage_map)
        ref_leaf_nodes = ts_phylo_seq.convert_entrez_to_tree_leaf_references(fasta_records)
        ts_updater.ref_pkg.taxa_trie.feed_leaf_nodes(ref_leaf_nodes)
        entrez_utils.sync_record_and_hierarchy_lineages(ref_leaf_nodes, fasta_records)
        ts_updater.ref_pkg.taxa_trie.validate_rank_prefixes()
        ts_updater.ref_pkg.taxa_trie.build_multifurcating_trie()
        classified_seq_lineage_map.update(
            ts_update_mod.guided_header_lineage_map(header_registry=querying_classified_fasta.header_registry,
                                                    entrez_records=fasta_records))
        querying_classified_fasta.synchronize_seqs_n_headers()
    else:
        # Map candidate reference sequence names to their TreeSAPP-assigned taxonomies
        classified_seq_lineage_map.update(
            {pq.place_name.split(' ')[0]: pq.recommended_lineage for pq in candidate_update_pqueries}
        )

    ref_header_map = {leaf.number + '_' + ts_updater.ref_pkg.prefix: leaf.description for leaf in ref_seq_lineage_info}
    ref_header_map = ts_update_mod.reformat_ref_seq_descriptions(ref_header_map)
    ref_name_lineage_map = {ref_header_map[leaf.number + '_' + ts_updater.ref_pkg.prefix]:
                                leaf.lineage for leaf in ref_seq_lineage_info}
    ref_accession_lineage_map = {ref_header_map[leaf.number + '_' + ts_updater.ref_pkg.prefix].split(' ')[0]:
                                     leaf.lineage for leaf in ref_seq_lineage_info}
    num_assigned_candidates = len(classified_seq_lineage_map)
    num_ref_seqs = len(ref_name_lineage_map)

    novel = set(classified_seq_lineage_map).difference(set(ref_accession_lineage_map.keys()))

    classified_seq_lineage_map.update(ref_accession_lineage_map)
    diff = num_ref_seqs + num_assigned_candidates - len(classified_seq_lineage_map)
    if diff > 0:
        LOGGER.warning("{} candidate sequences are already in the reference package."
                       " These will be excluded from any further analysis.\n".format(diff))
        # Remove the classified sequences that are redundant with the reference set
        classified_fasta.change_dict_keys("accession")
        classified_fasta.keep_only(list(novel))
        classified_fasta.change_dict_keys("original")
    elif diff < 0:
        LOGGER.error("Something's not adding up between the reference (%d), candidate (%d) and complete (%d) "
                     "sequence collections. Reference and candidate should sum to equal complete.\n" %
                     (num_ref_seqs, num_assigned_candidates, len(classified_seq_lineage_map)))
        sys.exit(13)

    ts_update_mod.validate_mixed_lineages(classified_seq_lineage_map)
    utilities.prepend_deep_rank(classified_seq_lineage_map)

    utilities.write_dict_to_table(classified_seq_lineage_map, ts_updater.lineage_map_file)

    ref_fasta = fasta.FASTA(ts_updater.ref_pkg.f__msa)
    ref_fasta.load_fasta()
    # Update the original reference headers using info from the tax_ids file
    ref_fasta.swap_headers(ref_header_map)
    ref_fasta.custom_lineage_headers(ref_name_lineage_map)

    combined_fasta = fasta.FASTA("")
    combined_fasta.clone(classified_fasta)
    combined_fasta.fasta_join(ref_fasta)
    combined_fasta.unalign()

    if args.resolve:
        ##
        # The purpose of this block is to remove any former candidate reference sequences from the ref_fasta object
        # that have a more truncated lineage that the new candidate reference sequences in combined_fasta
        ##
        combined_fasta.change_dict_keys("num_id")
        # Write a FASTA for clustering containing the formatted headers since
        # not all clustering tools + versions keep whole header - spaces are replaced with underscores
        fasta.write_new_fasta(fasta_dict=combined_fasta.fasta_dict,
                              fasta_name=ts_updater.cluster_input)
        wrapper.cluster_sequences(ts_updater.executables["mmseqs"], ts_updater.cluster_input,
                                  ts_updater.clusters_prefix, ts_updater.prop_sim)
        clusters_table = ts_updater.clusters_prefix + "_cluster.tsv"
        cluster_alignments = ts_updater.clusters_prefix + "_cluster_aln.tsv"

        cluster_dict = seq_clustering.create_mmseqs_clusters(clusters_tbl=clusters_table, aln_tbl=cluster_alignments)

        # Revert headers in cluster_dict from 'formatted' back to 'original'
        fasta.rename_cluster_headers(cluster_dict, combined_fasta.header_registry)
        LOGGER.debug("\t" + str(len(cluster_dict.keys())) + " sequence clusters\n")

        # Calculate LCA of each cluster to represent the taxonomy of the representative sequence
        entrez_records = ts_update_mod.simulate_entrez_records(combined_fasta, classified_seq_lineage_map)
        ts_create_mod.find_cluster_lca(cluster_dict, entrez_records, combined_fasta.header_registry)

        # Ensure centroids are the original reference sequences and skip clusters with identical lineages
        ts_update_mod.prefilter_clusters(cluster_dict,
                                         lineage_lookup={er.rebuild_header(): er.lineage for (num_id, er) in
                                                         entrez_records.items()},
                                         priority=list(ref_fasta.original_header_map().keys()))

        # Ensure the EntrezRecord with the most resolved lineage is the representative
        ts_update_mod.resolve_cluster_lineages(cluster_dict, entrez_records, ts_updater.ref_pkg.taxa_trie)

        if args.headless:
            ts_create_mod.finalize_cluster_reps(cluster_dict, entrez_records, combined_fasta.header_registry)
        else:
            # Allow user to select the representative sequence based on organism name, sequence length and similarity
            ts_create_mod.present_cluster_rep_options(cluster_dict, entrez_records, combined_fasta.header_registry,
                                                      important_seqs=combined_fasta.amendments, each_lineage=True)

        # Remove sequences that were replaced by resolve from ts_updater.old_ref_fasta
        still_repping = []
        refs_resolved = []
        for num_id, ref_seq in entrez_records.items():  # type: (str, entrez_utils.EntrezRecord)
            if ref_seq.cluster_rep:
                still_repping.append(ref_seq.rebuild_header())
            elif ref_seq.rebuild_header() in set(ref_fasta.fasta_dict.keys()):
                refs_resolved.append(ref_seq)
            else:
                pass

        try:
            classified_fasta.keep_only(still_repping, superset=True)
        except AssertionError:
            LOGGER.warning("No assigned sequences were retained for updating the reference package. Stopping now.\n")
            return

        ref_fasta.keep_only(still_repping, superset=True)
        combined_fasta.change_dict_keys("original")
        combined_fasta.keep_only(still_repping)

        if refs_resolved:
            LOGGER.info("{} reference sequences were resolved by updating sequences:\n\t"
                        "{}\n".format(len(refs_resolved),
                                      "\n\t".join([ref_seq.accession + ' ' + ref_seq.description
                                                   for ref_seq in refs_resolved])))

    if classified_fasta.n_seqs() > 0:
        LOGGER.info("{} assigned sequence(s) will be used in the update.\n".format(classified_fasta.n_seqs()))

    # Write only the sequences that have been properly classified
    combined_fasta.change_dict_keys("original")
    fasta.write_new_fasta(combined_fasta.fasta_dict, ts_updater.combined_fasta)
    fasta.write_new_fasta(ref_fasta.fasta_dict, ts_updater.old_ref_fasta)

    ##
    # Call create to create a new, updated reference package where the new sequences are guaranteed
    ##
    create_cmd = ts_update_mod.formulate_create_command(ts_updater, args, final_stage="support")
    create(create_cmd)
    ts_updater.updated_refpkg.f__pkl = ts_updater.updated_refpkg_path
    ts_updater.updated_refpkg.slurp()

    if ts_updater.stage_status("train"):
        train(ts_create_mod.formulate_train_command(input_seqs=ts_updater.input_sequences,
                                                    ref_pkg=ts_updater.updated_refpkg,
                                                    output_dir=ts_updater.training_dir,
                                                    args=args))
        ts_updater.updated_refpkg.f__pkl = os.path.join(ts_updater.training_dir, "final_outputs",
                                                        ts_updater.ref_pkg.prefix + ts_updater.ref_pkg.refpkg_suffix)
        ts_updater.updated_refpkg.slurp()
    else:
        ts_updater.updated_refpkg.training_df = ts_updater.ref_pkg.training_df
        ts_updater.updated_refpkg.pfit = ts_updater.ref_pkg.pfit
        ts_updater.updated_refpkg.svc = ts_updater.ref_pkg.svc

    ts_update_mod.update_features(old_refpkg=ts_updater.ref_pkg, new_refpkg=ts_updater.updated_refpkg)

    ##
    # Summarize some key parts of the new reference package, compared to the old one
    ##
    ts_updater.updated_refpkg.validate()
    ts_updater.update_refpkg_fields(output_dir=os.path.dirname(ts_updater.updated_refpkg_path))

    return


def colour(sys_args):
    parser = treesapp_args.TreeSAPPArgumentParser(description="Generates colour style and strip files for visualizing "
                                                              "a reference package's phylogeny in iTOL based on "
                                                              "taxonomic or phenotypic data.")
    treesapp_args.add_colour_arguments(parser)
    args = parser.parse_args(sys_args)

    log_file_name = os.path.join(args.output, "TreeSAPP_colour_log.txt")
    logger.prep_logging(log_file_name, args.verbose)
    LOGGER.info("\n##\t\t\tPainting a phylogeny\t\t\t##\n")

    treesapp_args.check_parser_arguments(args, sys_args)

    ts_painter = paint.PhyPainter()
    ts_painter.primer(args)

    # Find the taxa that should be coloured for each reference package
    for refpkg_name, ref_pkg in ts_painter.refpkg_dict.items():  # type: (str, ts_ref_pkg.ReferencePackage)
        if args.attribute != 'taxonomy':
            taxon_leaf_map = {}
            try:
                clade_annots = ref_pkg.feature_annotations[ts_painter.feature_name]
            except KeyError:
                LOGGER.warning("Reference package '{}' doesn't have the '{}' feature annotated. "
                               "It is being skipped\n".format(refpkg_name, ts_painter.feature_name))
                continue
            for ca in clade_annots:
                taxon_leaf_map.update({ca.name: list(ca.members)})
        else:
            ts_painter.find_rank_depth(ref_pkg, ref_pkg.taxa_trie.accepted_ranks_depths[ts_painter.rank])

            taxon_leaf_map, unique_taxa = ref_pkg.map_rank_representatives_to_leaves(rank_name=ts_painter.rank)
            LOGGER.info("{}: {} unique taxa.\n".format(ref_pkg.prefix, len(taxon_leaf_map)))
            internal_node_map = ref_pkg.get_internal_node_leaf_map()
            ts_painter.num_taxa = len(taxon_leaf_map)
            ts_painter.num_seqs = len(unique_taxa)

            # Begin filtering leaf nodes
            if args.taxa_filter:
                taxa = ts_painter.filter_unwanted_taxa(taxon_leaf_map, unique_taxa, args.taxa_filter)
                ts_painter.remove_taxa_from_colours(taxon_leaf_map, unique_taxa, taxa)
            if args.no_poly:
                taxa = ts_painter.filter_polyphyletic_groups(taxon_leaf_map=taxon_leaf_map,
                                                             internal_node_map=internal_node_map)
                ts_painter.remove_taxa_from_colours(taxon_leaf_map, unique_taxa, taxa)
            if args.min_prop:
                taxa = ts_painter.filter_rare_groups(taxon_leaf_map, ref_pkg.num_seqs, args.min_prop)
                ts_painter.remove_taxa_from_colours(taxon_leaf_map, unique_taxa, taxa)

        ts_painter.refpkg_leaf_nodes_to_colour[refpkg_name] = taxon_leaf_map
        # Find the intersection or union between reference packages analyzed so far
        ts_painter.harmonize_taxa_colours(taxon_leaf_map, args.set_op)

    if len(ts_painter.refpkg_leaf_nodes_to_colour.keys()) == 0:
        LOGGER.error("Unable to colour phylogenies by '{}' - attributes were not found in reference packages.\n"
                     "".format(args.attribute))
        raise AssertionError

    # Sort the nodes by their internal node order
    taxa_order = paint.order_taxa(taxa_to_colour=ts_painter.taxa_to_colour,
                                  taxon_leaf_map=ts_painter.refpkg_leaf_nodes_to_colour[ref_pkg.prefix],
                                  leaf_order=ref_pkg.leaf_node_order(),
                                  method=ts_painter.order_method)

    # Determine the palette to use for taxa across all reference packages
    palette_taxa_map = ts_painter.map_colours_to_taxa(taxa_order)

    # Create the iTOL colour files
    for refpkg_name, ref_pkg in ts_painter.refpkg_dict.items():  # type: (str, ts_ref_pkg.ReferencePackage)
        taxon_leaf_map = ts_painter.refpkg_leaf_nodes_to_colour[refpkg_name]
        ts_painter.add_unknowns_to_feature_leaf_map(taxon_leaf_map, ref_pkg)
        style_file = os.path.join(ts_painter.output_dir,
                                  "{}_{}_colours_style.txt".format(ref_pkg.prefix, ts_painter.feature_name))
        strip_file = os.path.join(ts_painter.output_dir,
                                  "{}_{}_colour_strip.txt".format(ref_pkg.prefix, ts_painter.feature_name))

        # Find the minimum set of monophyletic internal nodes for each taxon
        taxa_clades = ts_painter.find_mono_clades(taxon_leaf_map, ref_pkg)
        taxa_ranges = paint.convert_clades_to_ranges(taxa_clades, ref_pkg.leaf_node_order())

        paint.write_colours_styles(taxon_leaf_map, palette_taxa_map, style_output=style_file)

        paint.write_colour_strip(taxa_ranges, palette_taxa_map,
                                 colour_strip_file=strip_file,
                                 data_label=ts_painter.feature_name)

    return


def layer(sys_args):
    # STAGE 1: Prompt the user and prepare files and lists for the pipeline
    parser = treesapp_args.TreeSAPPArgumentParser(description="This script adds extra feature annotations, such as "
                                                              "Subgroup and Metabolic Pathway, to an existing "
                                                              "classification table made by treesapp assign. "
                                                              "A new column is bound to the table for each feature.")
    treesapp_args.add_layer_arguments(parser)
    args = parser.parse_args(sys_args)

    ts_layer = annotate_extra.Layerer()

    log_file_name = args.output + os.sep + "TreeSAPP_layer_log.txt"
    logger.prep_logging(log_file_name, args.verbose)
    LOGGER.info("\n##\t\t\t\tLayering extra annotations on TreeSAPP classifications\t\t\t\t##\n\n")
    LOGGER.info("Arguments used:\n" + ' '.join(sys_args) + "\n")

    ts_layer.check_arguments(args)

    ##
    # Worklow:
    #   1. Slurp up reference packages into a dictionary
    #   2. Read the classifications.tsv file from the output directory to create the master data structure
    #   3. For each of the colours_styles files provided (potentially multiple for the same marker):
    #       3.1) Add the annotation variable to master_dat for every sequence (instantiate with "NA")
    #       3.2) Read the .jplace file for every sequence classified as marker
    #       3.3) Add the annotation information to every sequence classified as marker in master_dat
    #   4. Write the new classification file called "layered_classifications.tsv"
    ##
    marker_subgroups = set()
    unique_markers_annotated = set()
    marker_tree_info = dict()
    master_dat, field_order = annotate_extra.parse_marker_classification_table(ts_layer.final_output_dir +
                                                                               ts_layer.classification_tbl_name)
    refpkg_dict = ts_ref_pkg.gather_ref_packages(ts_layer.refpkg_dir)

    # structure of master dat:
    # {"Sequence_1": {"Field1": x, "Field2": y, "Extra": n},
    #  "Sequence_2": {"Field1": i, "Field2": j, "Extra": n}}
    for refpkg_prefix, ref_pkg in refpkg_dict.items():  # type: (str, ts_ref_pkg.ReferencePackage)
        if refpkg_prefix not in master_dat.keys():
            continue
        if len(ref_pkg.feature_annotations) > 0:
            unique_markers_annotated.add(refpkg_prefix)
            for data_type in ref_pkg.feature_annotations:
                marker_subgroups.add(data_type)

    # Instantiate every query sequence in classifications with an empty string for each data_type
    for data_type in marker_subgroups:
        for refpkg_name in master_dat:
            for assignment in master_dat[refpkg_name]:  # type: annotate_extra.ClassifiedSequence
                assignment.layers[data_type] = "NA"

    # Update the field_order dictionary with new fields
    field_acc = len(field_order)
    for new_datum in sorted(marker_subgroups):
        field_order[field_acc] = new_datum
        field_acc += 1

    # Load the query sequence annotations
    for data_type in marker_subgroups:
        LOGGER.info("Annotating '%s' classifications for the following reference package(s):\n" % data_type)
        if data_type not in marker_tree_info:
            marker_tree_info[data_type] = dict()
        for refpkg_name in unique_markers_annotated:  # type: str
            jplace = os.path.join(ts_layer.treesapp_output, "iTOL_output", refpkg_name,
                                  refpkg_name + "_complete_profile.jplace")
            ref_pkg = refpkg_dict[refpkg_name]  # type: ts_ref_pkg.ReferencePackage
            if data_type in ref_pkg.feature_annotations:
                LOGGER.info("\t" + refpkg_name + "\n")
                # Create the dictionary mapping an internal node to all leaves
                internal_node_map = entish.map_internal_nodes_leaves(jplace_utils.jplace_parser(jplace).tree)
                # Routine for exchanging any organism designations for their respective node number
                taxa_map = ref_pkg.generate_tree_leaf_references_from_refpkg()

                annotated_edges, leaves_in_clusters = annotate_extra.annotate_internal_nodes(internal_node_map,
                                                                                             ref_pkg.feature_annotations[
                                                                                                 data_type])
                marker_tree_info[data_type][refpkg_name] = annotated_edges

                diff = len(taxa_map) - len(leaves_in_clusters)
                if diff != 0:
                    unannotated = set()
                    for inode in internal_node_map:
                        for leaf in internal_node_map[inode]:
                            if leaf not in leaves_in_clusters:
                                unannotated.add(str(leaf))
                    LOGGER.warning("{} leaf nodes were not mapped to annotation groups. "
                                   "More information can be found in the log.\n".format(len(unannotated)))
                    LOGGER.debug("The following leaf nodes were not mapped to annotation groups:\n" +
                                 "\t" + ', '.join(sorted(unannotated, key=lambda x: int(x.split('_')[0]))) + "\n")
            else:
                pass
    marker_subgroups.clear()
    annotate_extra.map_queries_to_annotations(marker_tree_info, master_dat, join=True)
    annotate_extra.write_classification_table(ts_layer.layered_table, field_order, master_dat)

    return


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
    logger.prep_logging(log_file_name, args.verbose)
    LOGGER.info("\n##\t\t\tBeginning purity analysis\t\t\t##\n")

    treesapp_args.check_parser_arguments(args, sys_args)
    ts_purity.check_purity_arguments(args)
    ts_purity.decide_stage(args)

    # Load FASTA data
    if ts_purity.stage_status("clean"):
        ref_seqs = fasta.FASTA(ts_purity.input_sequences)
        ref_seqs.load_fasta()
        ref_seqs.change_dict_keys()
        ref_seqs.unalign()
        fasta.write_new_fasta(ref_seqs.fasta_dict, ts_purity.formatted_input)
        ts_purity.increment_stage_dir()

    if ts_purity.stage_status("assign"):
        assign_args = ["-i", ts_purity.formatted_input, "-o", ts_purity.stage_output_dir,
                       "-m", ts_purity.molecule_type, "-n", str(args.num_threads),
                       "-t", ts_purity.ref_pkg.prefix, "--refpkg_dir", ts_purity.refpkg_dir,
                       "--overwrite", "--delete"]
        if args.trim_align:
            assign_args.append("--trim_align")

        try:
            ts_assign_mod.assign(assign_args)
        except:  # Just in case treesapp assign fails, just continue
            LOGGER.error("TreeSAPP failed.\n")
        ts_purity.increment_stage_dir()

    if ts_purity.stage_status("summarize"):
        metadat_dict = dict()
        # Parse classification table and identify the groups that were assigned
        if os.path.isfile(ts_purity.classifications):
            assigned_lines = file_parsers.read_classification_table(ts_purity.classifications)
            ts_purity.assignments = file_parsers.parse_assignments(assigned_lines)
        else:
            LOGGER.error("{} is missing from output directory '{}'\n"
                         "Please remove this directory and re-run.\n"
                         "".format(ts_purity.classification_tbl_name,
                                   os.path.dirname(ts_purity.classifications)))
            sys.exit(5)

        LOGGER.info("\nSummarizing assignments for reference package " + ts_purity.ref_pkg.prefix + "\n")
        # If an information table was provided, map the metadata to classified markers
        if ts_purity.metadata_file:
            metadat_dict.update(ts_purity.load_metadata())
        # Identify the number of sequences that are descendents of each orthologous group
        jplace_data = jplace_utils.jplace_parser(ts_purity.assign_jplace_file)
        jplace_data.pqueries = jplace_utils.demultiplex_pqueries(jplace_data)
        node_map = entish.map_internal_nodes_leaves(jplace_data.tree)
        labelled_ref_tree = ts_purity.ref_pkg.taxonomically_label_tree()
        for pquery in jplace_data.pqueries:  # type: ts_phylo_seq.PQuery
            pquery.process_max_weight_placement(labelled_ref_tree)
        ortholog_map = ts_purity.assign_leaves_to_orthologs(jplace_data.pqueries, node_map)
        ts_purity.summarize_groups_assigned(ortholog_map, metadat_dict)

        # Write each sequence name that can be assigned to an ortholog to the log
        summary_str = ""
        for ortholog_name in sorted(ortholog_map, key=lambda x: len(ortholog_map[x])):
            summary_str += ortholog_name + ":\n\t"
            summary_str += "\n\t".join(ortholog_map[ortholog_name]) + "\n"
        LOGGER.debug(summary_str)

    return


def evaluate(sys_args):
    """
    Provide it a FASTA file for which it will determine the taxonomic lineage for each sequence
    and run all taxonomic representative sequences with TreeSAPP then analyze via clade exclusion

    :return:
    """
    parser = treesapp_args.TreeSAPPArgumentParser(
        description='Evaluate classification performance using clade-exclusion analysis.')
    treesapp_args.add_evaluate_arguments(parser)
    args = parser.parse_args(sys_args)

    ts_evaluate = classy.Evaluator()
    ts_evaluate.furnish_with_arguments(args)
    ts_evaluate.check_previous_output(args.overwrite)

    log_file_name = args.output + os.sep + "TreeSAPP_evaluation_log.txt"
    logger.prep_logging(log_file_name, args.verbose)
    LOGGER.info("\n##\t\t\tBeginning clade exclusion analysis\t\t\t##\n")

    treesapp_args.check_parser_arguments(args, sys_args)
    treesapp_args.check_evaluate_arguments(ts_evaluate, args)
    ts_evaluate.decide_stage(args)
    ts_evaluate.rank_depth_map = {ts_evaluate.ref_pkg.taxa_trie.accepted_ranks_depths[rank_name]: rank_name for
                                  rank_name in ts_evaluate.ref_pkg.taxa_trie.accepted_ranks_depths}

    ref_leaves = ts_evaluate.ref_pkg.generate_tree_leaf_references_from_refpkg()
    ref_lineages = {leaf.number: leaf.lineage for leaf in ref_leaves}

    # Load FASTA data
    query_fasta = fasta.FASTA(ts_evaluate.input_sequences)
    query_fasta.load_fasta(format_it=True, molecule=ts_evaluate.molecule_type)
    if args.length:
        query_fasta.trim_to_length(args.length)

    fasta_records = ts_evaluate.fetch_entrez_lineages(query_fasta, ts_evaluate.molecule_type, args.acc_to_taxid)
    entrez_utils.fill_ref_seq_lineages(fasta_records, ts_evaluate.seq_lineage_map)
    query_leaf_nodes = ts_phylo_seq.convert_entrez_to_tree_leaf_references(fasta_records)
    ts_evaluate.ref_pkg.taxa_trie.feed_leaf_nodes(query_leaf_nodes)
    entrez_utils.sync_record_and_hierarchy_lineages(query_leaf_nodes, fasta_records)
    ts_evaluate.ref_pkg.taxa_trie.validate_rank_prefixes()
    ts_evaluate.ref_pkg.taxa_trie.build_multifurcating_trie()

    LOGGER.info("Selecting representative sequences for each taxon... ")
    # Filter the sequences from redundant taxonomic lineages, picking up to 5 representative sequences
    representative_seqs, ts_evaluate.taxa_filter = ts_clade_ex.pick_taxonomic_representatives(fasta_records,
                                                                                              ts_evaluate.taxa_filter)
    LOGGER.info("done.\n")
    LOGGER.info("\t{} representative sequences may be used by TreeSAPP evaluate.\n".format(len(representative_seqs)))

    deduplicated_fasta_dict = ts_clade_ex.select_rep_seqs(representative_seqs,
                                                          test_sequences=fasta_records,
                                                          taxon_hierarchy=ts_evaluate.ref_pkg.taxa_trie)
    fasta.write_new_fasta(deduplicated_fasta_dict, ts_evaluate.test_rep_taxa_fasta)
    rep_accession_lineage_map = ts_clade_ex.map_seqs_to_lineages(ts_evaluate.seq_lineage_map, deduplicated_fasta_dict)

    # Ensure that both lineage sets are rooted
    ts_clade_ex.check_lineage_compatibility(ref_lineages, ts_evaluate.ref_pkg.taxa_trie)
    ts_clade_ex.check_lineage_compatibility(rep_accession_lineage_map, ts_evaluate.ref_pkg.taxa_trie)
    # Checkpoint three: We have accessions linked to taxa, and sequences to analyze with TreeSAPP, but not classified
    if ts_evaluate.stage_status("classify"):
        # Run TreeSAPP against the provided tax_ids file and the unique taxa FASTA file
        for rank in args.taxon_rank:
            candidate_lineages = ts_clade_ex.get_testable_lineages_for_rank(ref_lineages, rep_accession_lineage_map,
                                                                            rank, ts_evaluate.ref_pkg.taxa_trie)
            p_bar = tqdm.tqdm(ncols=100, desc="Evaluating {}".format(rank), total=len(candidate_lineages))
            for lineage in candidate_lineages:
                # Select representative sequences belonging to the taxon being tested
                taxon_rep_seqs = ts_clade_ex.select_rep_seqs(representative_seqs, fasta_records,
                                                             taxon_hierarchy=ts_evaluate.ref_pkg.taxa_trie,
                                                             target_lineage=lineage)
                # Decide whether to continue analyzing taxon based on number of query sequences
                if len(taxon_rep_seqs.keys()) == 0:
                    LOGGER.debug("No query sequences for {}.\n".format(lineage))
                    p_bar.update()
                    continue

                # Continuing with classification
                tt_obj = ts_evaluate.new_taxa_test(lineage, args.tool)
                tt_obj.queries = taxon_rep_seqs.keys()

                LOGGER.debug("Classifications for '{}' put in {}\n".format(tt_obj.taxon, tt_obj.intermediates_dir))
                if args.tool in ["graftm", "diamond"]:
                    ts_clade_ex.run_clade_exclusion_graftm(tt_obj, taxon_rep_seqs, ts_evaluate.ref_pkg,
                                                           graftm_classifier=args.tool,
                                                           num_threads=args.num_threads,
                                                           executables=ts_evaluate.executables)
                else:
                    ts_clade_ex.run_clade_exclusion_treesapp(tt_obj, taxon_rep_seqs, ts_evaluate.ref_pkg,
                                                             executables=ts_evaluate.executables,
                                                             fresh=args.fresh,
                                                             cl_assign_args=args,
                                                             min_seq_length=ts_evaluate.min_seq_length)
                p_bar.update()
            p_bar.close()

    if ts_evaluate.stage_status("calculate"):
        # everything has been prepared, only need to parse the classifications and map lineages
        LOGGER.info("Finishing up the mapping of classified, filtered taxonomic sequences.\n")
        for rank in sorted(ts_evaluate.taxa_tests):
            for tt_obj in ts_evaluate.taxa_tests[rank]:  # type: classy.TaxonTest
                if tt_obj.assignments:
                    marker_assignments = ts_clade_ex.map_headers_to_lineage(tt_obj.assignments, fasta_records)
                    # Return the number of correct, classified, and total sequences of that taxon at the current rank
                    # Identify the excluded rank for each query sequence
                    if len(marker_assignments) == 0:
                        LOGGER.debug("No '{}' sequences were classified.\n".format(tt_obj.taxon))
                        continue

                    for marker in marker_assignments:
                        ts_evaluate.markers.add(marker)

                    if ts_evaluate.ref_pkg.prefix not in marker_assignments:
                        LOGGER.debug("No sequences assigned as " + ts_evaluate.ref_pkg.prefix + "\n")
                        rank_assignments = {}
                    else:
                        rank_assignments = lca_calculations.identify_excluded_clade(
                            marker_assignments[ts_evaluate.ref_pkg.prefix],
                            tt_obj.taxonomic_tree)

                    for a_rank in rank_assignments:
                        if a_rank != rank and len(rank_assignments[a_rank]) > 0:
                            LOGGER.warning("{}-level clade excluded but optimal classification found to be {}-level.\n"
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
            LOGGER.debug("\n\t" + str(ts_evaluate.taxa_filter["Classified"] -
                                      ts_evaluate.taxa_filter["Unique_taxa"]) +
                         " duplicate query taxonomies removed.\n")

        if ts_evaluate.taxa_filter["Unclassified"] > 0:
            LOGGER.debug("\t" + str(ts_evaluate.taxa_filter["Unclassified"]) +
                         " query sequences with unclassified taxonomies were removed.\n" +
                         "This is not a problem, its just they have 'unclassified' somewhere in their lineages\n" +
                         "(e.g. Unclassified Bacteria) and this is not good for assessing placement accuracy.\n\n")

        if ts_evaluate.ref_pkg.prefix not in ts_evaluate.markers:
            LOGGER.error("No sequences were classified as {}.\n".format(ts_evaluate.ref_pkg.prefix))
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
        containment_strings = ts_clade_ex.determine_containment(ts_evaluate)
        ts_evaluate.write_containment_table(containment_strings, args.tool)

    if args.delete and os.path.isdir(ts_evaluate.var_output_dir):
        shutil.rmtree(ts_evaluate.var_output_dir)

    return
