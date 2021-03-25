import os
import sys
import logging

from samsum import commands as samsum_cmd

from treesapp import classy
from treesapp import wrapper
from treesapp import treesapp_args
from treesapp import phylo_seq
from treesapp import file_parsers


def abundance(sys_args) -> dict:
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

    log_file_name = args.output + os.sep + "TreeSAPP_abundance_log.txt"
    classy.prep_logging(log_file_name, args.verbose)
    logging.info("\n##\t\t\tCalculating abundance of classified sequences\t\t\t##\n")

    treesapp_args.check_parser_arguments(args, sys_args)
    ts_abund.check_arguments(args)
    ts_abund.decide_stage(args)

    abundance_dict = {}
    rev_reads = None

    while args.reads:
        fwd_reads = args.reads.pop(0)
        ts_abund.strip_file_to_sample_name(fwd_reads)
        ts_abund.aln_file = ts_abund.stage_lookup("align_map").dir_path + ts_abund.sample_prefix + ".sam"
        logging.info("Working on sample '{}'.\n".format(ts_abund.sample_prefix))
        if args.pairing == 'pe' and args.reverse:
            rev_reads = args.reverse.pop(0)

        if ts_abund.stage_status("align_map"):
            wrapper.align_reads_to_nucs(ts_abund.executables["bwa"], ts_abund.ref_nuc_seqs,
                                        ts_abund.stage_lookup("align_map").dir_path, fwd_reads, args.pairing,
                                        reverse=rev_reads, sam_file=ts_abund.aln_file, num_threads=args.num_threads)
            ts_abund.increment_stage_dir()

        if ts_abund.stage_status("sam_sum"):
            if os.path.isfile(ts_abund.aln_file):
                ref_seq_abunds = samsum_cmd.ref_sequence_abundances(aln_file=ts_abund.aln_file,
                                                                    seq_file=ts_abund.ref_nuc_seqs,
                                                                    min_aln=10, p_cov=50, map_qual=1, multireads=False)
                if args.metric == "fpkm":
                    abundance_dict[ts_abund.sample_prefix] = {ref_seq.name: ref_seq.fpkm for ref_seq
                                                              in ref_seq_abunds.values()}
                elif args.metric == "tpm":
                    abundance_dict[ts_abund.sample_prefix] = {ref_seq.name: ref_seq.tpm for ref_seq
                                                              in ref_seq_abunds.values()}
                else:
                    logging.error("Unrecognized normalization metric '{}'.\n".format(args.metric))
                    sys.exit(7)

                ref_seq_abunds.clear()
            else:
                logging.warning("SAM file '%s' was not generated.\n" % ts_abund.aln_file)
                return {}
            ts_abund.increment_stage_dir()
        else:
            abundance_dict = {}
            logging.warning("Skipping samsum normalized relative abundance calculation in treesapp abundance.\n")

    ts_abund.delete_intermediates(args.delete)

    if args.report != "nothing" and os.path.isfile(ts_abund.classifications):
        refpkg_pqueries = file_parsers.load_classified_sequences_from_assign_output(ts_abund.output_dir)

        # Collect the names of all classified PQueries
        classified_seqs = set()
        for rp_pqs in refpkg_pqueries.values():
            classified_seqs.update({pq.seq_name for pq in rp_pqs})

        # Select a unique set of PQuery instances to use in the updated classification table
        for refpkg_name, pqueries in refpkg_pqueries.items():
            unique_pqueries = []
            seq_name_list = list({pq.place_name for pq in pqueries})
            for pquery in pqueries:
                i = 0
                while i < len(seq_name_list):
                    if pquery.place_name == seq_name_list[i]:
                        unique_pqueries.append(pquery)
                        seq_name_list.pop(i)
                        i = len(seq_name_list)
                    i += 1
            refpkg_pqueries[refpkg_name] = unique_pqueries

        # Fill PQuery.abundance attribute
        for sample_name, abundance_map in abundance_dict.items():
            # Filter the abundance_map dictionary to just the classified sequences
            absentee = set(abundance_map.keys()).difference(classified_seqs)
            while absentee:
                abundance_map.pop(absentee.pop())
            phylo_seq.abundify_tree_saps(refpkg_pqueries, abundance_map)
            file_parsers.write_classification_table(refpkg_pqueries, sample_name, ts_abund.classifications,
                                                    append=ts_abund.append_abundance)
            ts_abund.append_abundance = True

    return abundance_dict
