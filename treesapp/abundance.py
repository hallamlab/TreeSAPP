import os
import logging

from samsum import commands as samsum_cmd

from treesapp import classy
from treesapp import wrapper
from treesapp import treesapp_args
from treesapp import phylo_seq
from treesapp import file_parsers


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
    ts_abund.check_previous_output(args.overwrite)

    log_file_name = args.output + os.sep + "TreeSAPP_abundance_log.txt"
    classy.prep_logging(log_file_name, args.verbose)
    logging.info("\n##\t\t\tCalculating abundance of classified sequences\t\t\t##\n")

    treesapp_args.check_parser_arguments(args, sys_args)
    ts_abund.check_arguments(args)
    ts_abund.decide_stage(args)

    if ts_abund.stage_status("align_map"):
        wrapper.align_reads_to_nucs(ts_abund.executables["bwa"], ts_abund.classified_nuc_seqs,
                                    ts_abund.stage_output_dir, args.reads, args.pairing, args.reverse,
                                    args.num_threads)
        ts_abund.increment_stage_dir()

    if ts_abund.stage_status("sam_sum"):
        if os.path.isfile(ts_abund.aln_file):
            ref_seq_abunds = samsum_cmd.ref_sequence_abundances(aln_file=ts_abund.aln_file,
                                                                seq_file=ts_abund.classified_nuc_seqs,
                                                                min_aln=10, p_cov=50, map_qual=1, multireads=False)
            abundance_dict = {ref_seq.name: ref_seq.fpkm for ref_seq in ref_seq_abunds.values()}
            ref_seq_abunds.clear()
        else:
            logging.warning("SAM file '%s' was not generated.\n" % ts_abund.aln_file)
            return {}
        ts_abund.increment_stage_dir()
    else:
        abundance_dict = {}
        logging.warning("Skipping samsum normalized relative abundance calculation in treesapp abundance.\n")

    ts_abund.delete_intermediates(args.delete)

    # TODO: Index each PQuery's abundance by the dataset name, write a new row for each dataset's abundance
    if args.report != "nothing" and os.path.isfile(ts_abund.classifications):
        pqueries = file_parsers.load_classified_sequences_from_assign_output(ts_abund.output_dir)
        phylo_seq.abundify_tree_saps(pqueries, abundance_dict)
        file_parsers.write_classification_table(pqueries, ts_abund.sample_prefix, ts_abund.classifications)

    return abundance_dict
