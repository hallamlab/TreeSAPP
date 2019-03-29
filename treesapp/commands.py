
import logging
import sys
import re
import os
from .file_parsers import parse_ref_build_params, read_species_translation_files
from .fasta import format_read_fasta
from .treesapp_args import TreeSAPPArgumentParser, check_parser_arguments, add_classify_arguments
from .utilities import executable_dependency_versions, find_executables
import treesapp.classify as ts_classify


def info(args):
    """

    """
    parser = TreeSAPPArgumentParser(description="Return package, executable and refpkg information.")
    args = parser.parse_args(args)

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
                 "\n\t".join([k + ": " + v for k, v in py_deps]) + "\n")

    # Write the version of executable deps
    exe_dict = find_executables(args)
    logging.info(executable_dependency_versions(exe_dict))

    # TODO: Write relevant reference package information (e.g. codes, gene names, descriptions)

    return


def create():
    return


def evaluate():
    return


def classify():
    sys.stdout.write("\n##\t\t\t\tTreeSAPP classify\t\t\t\t##\n\n")
    sys.stdout.flush()
    # STAGE 1: Prompt the user and prepare files and lists for the pipeline
    parser = TreeSAPPArgumentParser(description='Taxonomically classify sequences through evolutionary placement.')
    add_classify_arguments(parser)
    args = parser.parse_args()

    args = check_parser_arguments(args)
    args = check_previous_output(args)

    marker_build_dict = parse_ref_build_params(args)
    tree_numbers_translation = read_species_translation_files(args, marker_build_dict)
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
        args.formatted_input_file = args.output_dir_var + input_multi_fasta + "_formatted.fasta"
        logging.debug("Writing formatted FASTA file to " + args.formatted_input_file + "... ")
        formatted_fasta_files = write_new_fasta(formatted_fasta_dict, args.formatted_input_file)
        logging.debug("done.\n")
        ref_alignment_dimensions = get_alignment_dims(args, marker_build_dict)

        # STAGE 3: Run hmmsearch on the query sequences to search for marker homologs
        hmm_domtbl_files = hmmsearch_orfs(args, marker_build_dict)
        hmm_matches = parse_domain_tables(args, hmm_domtbl_files)
        extracted_seq_dict, numeric_contig_index = extract_hmm_matches(hmm_matches, formatted_fasta_dict)
        homolog_seq_files = write_grouped_fastas(extracted_seq_dict, numeric_contig_index,
                                                 marker_build_dict, args.output_dir_var)

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
        orf_nuc_fasta = args.output_dir_final + sample_name + "_classified_seqs.fna"
        if not os.path.isfile(orf_nuc_fasta):
            logging.info("Creating nucleotide FASTA file of classified sequences '" + orf_nuc_fasta + "'... ")
            genome_nuc_genes_file = args.output_dir_final + sample_name + "_ORFs.fna"
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
            abundance_dict = read_rpkm(rpkm_output_file)
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
    logging.info("TreeSAPP has finished successfully.\n")

    return


def update():
    return


def train():
    return