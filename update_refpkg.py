import os
import shutil
import sys
import logging
import re
import random
import subprocess
from external_command_interface import launch_write_command
from entrez_utils import get_multiple_lineages, entrez_records_to_accession_lineage_map
from utilities import reformat_string, return_sequence_info_groups
from fasta import format_read_fasta, get_headers, get_header_format, write_new_fasta
from classy import MarkerBuild, ReferenceSequence

# TODO: Refactor this entire class - loads of redundant code


class CreateFuncTreeUtility:
    """
    Output is the directory to write the outputs for the updated tree
    InputData is the path to the TreeSAPP output folder containing, various_outputs/ and final_RAxML_outputs/
    RefTree is the second column in cog_list.tsv for the gene to update
    Cluster is a flag indicating whether the protein sequences for the RefTree in InputData is to be clustered at 97%
    """
    def __init__(self, input_data, ref_marker: MarkerBuild):
        if os.path.isabs(input_data):
            self.InputData = input_data
        else:
            self.InputData = os.getcwd() + os.sep + input_data

        if self.InputData[-1] == '/':
            self.InputData = self.InputData[:-1]

        self.Output = self.InputData + os.sep + "updated_" + ref_marker.denominator + "_tree" + os.sep
        self.Denominator = ref_marker.denominator
        self.marker_molecule = ref_marker.molecule  # prot, dna, or rrna
        self.COG = ref_marker.cog
        self.raxml_model = ref_marker.model  # Use the original model rather than determining the best model again
        self.ContigDict = dict()  # Used for storing the final candidate reference sequences
        self.names = list()
        self.header_id_map = dict()  # Used for mapping the original header of the sequence to internal numeric header
        self.master_reference_index = dict()  # Used for storing ReferenceSequence objects, indexed by numeric IDs
        self.cluster_id = ref_marker.pid  # The percent similarity the original reference sequences were clustered at
        # Automatically remove the last attempt at updating the reference tree
        if os.path.isdir(self.Output):
            shutil.rmtree(self.Output)
        try:
            os.makedirs(self.Output)
        except IOError:
            logging.error("Unable to make the directory " + str(self.Output) + "\n")
            sys.exit(17)

    def get_raxml_files_for_ref(self):
        """
        Uses self.InputData to find all the RAxML_outputs for each protein sequence for self.RefTree
        :return: list of file names with a protein sequence of self.RefTree
        """
        raxml_outputs = os.listdir(self.InputData + os.sep + "final_RAxML_outputs" + os.sep)
        for placement in raxml_outputs:
            ref_tree = os.path.basename(placement).split('_')[0]
            if ref_tree == self.Denominator:
                suffix = re.sub("%s_" % ref_tree, '', placement)
                predicted_orf = re.sub("_RAxML_parsed.txt", '', suffix)
                self.names.append(predicted_orf)
        return

    def load_new_refs_fasta(self, args, centroids_fasta, ref_organism_lineage_info):
        # Determine whether the numerical IDs are a series and sorted
        acc = 0
        ref_numeric_ids = list()
        for leaf in ref_organism_lineage_info[self.Denominator]:
            ref_numeric_ids.append(int(leaf.number))
            ref_seq = ReferenceSequence()
            ref_seq.lineage = leaf.lineage
            ref_seq.organism = leaf.description.split(' | ')[0]
            ref_seq.accession = leaf.description.split(' | ')[1]
            ref_seq.short_id = '>' + leaf.number + '_' + self.COG
            self.master_reference_index[leaf.number] = ref_seq
            acc += 1

        # Read the FASTA to get headers and sequences
        centroids_fasta_dict = format_read_fasta(centroids_fasta, self.marker_molecule, args.output)

        # Create the final contig dictionary of new internal TreeSAPP headers (keys) and sequences (values)
        additional = acc
        for new_ref_seq in centroids_fasta_dict:
            if acc == sorted(ref_numeric_ids, key=int)[-1]:
                additional += 1
                # These sequences form a series, therefore continue the series for new reference sequences
                internal_id = '>' + str(additional) + '_' + self.COG
            else:
                rfive = ''.join(str(x) for x in random.sample(range(10), 5))
                while rfive in ref_numeric_ids:
                    rfive = ''.join(str(x) for x in random.sample(range(10), 5))
                internal_id = '>' + rfive + '_' + self.COG
            ref_seq = ReferenceSequence()
            ref_seq.short_id = internal_id
            # Map the new reference headers to their numerical IDs
            self.header_id_map[internal_id] = new_ref_seq
            self.ContigDict[internal_id] = centroids_fasta_dict[new_ref_seq]
            ref_seq.sequence = centroids_fasta_dict[new_ref_seq]
            self.master_reference_index[str(additional)] = ref_seq

        if additional == acc:
            # The reference sequence identifiers are random
            logging.warning("Numerical TreeSAPP identifiers for " + self.Denominator +
                            "are not in the format of a sequential series!\n" +
                            "Generating random numerical unique identifiers for the new sequence(s).\n")

        logging.debug("\t" + str(len(self.header_id_map)) + " new " + self.COG + " reference sequences.\n")

        return

    def update_tax_ids(self, args, ref_organism_lineage_info, assignments):
        """
        Write the number, organism and accession ID, if possible
        :param args:
        :param ref_organism_lineage_info: A dictionary mapping the
        :param assignments: A dictionary containing marker name as keys, and header: [lineages] mappings as values
        :return:
        """
        logging.info("Writing updated tax_ids file... ")

        tree_taxa_string = ""
        original_to_formatted_header_map = dict()
        unclassified_seqs = list()

        if self.Denominator not in ref_organism_lineage_info.keys():
            raise ValueError(self.Denominator + " not included in data from tax_ids files!\n")

        # Load the original reference sequences first since this shouldn't change
        for leaf in ref_organism_lineage_info[self.Denominator]:
            if leaf.lineage:
                tree_taxa_string += '\t'.join([str(leaf.number), leaf.description, leaf.lineage]) + "\n"
            else:
                logging.debug("Unable to retrieve lineage information for sequence " + str(leaf.number) + "\n")

        # Build a map of original headers to formatted ones
        original_candidate_headers = get_headers(args.fasta_input)
        for original in original_candidate_headers:
            original_to_formatted_header_map[reformat_string(original)] = original

        # Now figure out what to do with the new reference sequences
        for reference in sorted(self.ContigDict.keys()):
            num_id = reference[1:].split('_')[0]
            reformatted_header = self.header_id_map[reference]
            header = original_to_formatted_header_map[reformatted_header]
            header_format_re, header_db, header_molecule = '', '', ''
            description, lineage, accession, organism = '', '', '', ''
            # Check to see if this header contains an accession value or if it is a contig name with no useful info
            try:
                header_format_re, header_db, header_molecule = get_header_format(header, self.Denominator)
            except AssertionError:
                # Just a contig name, so it will need to be formatted with explicit lineage
                for candidate_header_name in assignments[self.COG]:
                    if header == candidate_header_name:
                        assigned_lineage = assignments[self.COG][candidate_header_name]
                        # Determine if this is a cellular organism or not
                        if assigned_lineage.split(';')[0] in ["Bacteria", "Archaea", "Eukaryota"]:
                            assigned_lineage = "cellular organisms; " + assigned_lineage
                        annotated_header = ' '.join([header,
                                                     "lineage=" + assigned_lineage,
                                                     "[" + re.sub('>', '', header),
                                                     self.InputData.split(os.sep)[-1],
                                                     assigned_lineage.split(';')[-1] + "]"])
                        header_format_re, header_db, header_molecule = get_header_format(annotated_header, self.Denominator)
                        header = annotated_header
                        break
                if not header_format_re:
                    unclassified_seqs.append(header)
            if header_format_re and header_db and header_molecule:
                sequence_info = header_format_re.match(header)
                if sequence_info:
                    accession, organism, locus, description, lineage = return_sequence_info_groups(sequence_info, header_db, header)
                    if lineage:
                        pass
                    elif accession:
                        entrez_records = get_multiple_lineages([accession], header_molecule)
                        accession_lineage_map = entrez_records_to_accession_lineage_map(entrez_records)
                        # Should only be one...
                        for tuple_key in accession_lineage_map:
                            lineage = accession_lineage_map[tuple_key]["lineage"]
                        try:
                            taxonomic_lineage = lineage.split('; ')
                            # Try to get the species name
                            if len(taxonomic_lineage) >= 8:
                                description = ' '.join(taxonomic_lineage[6:8])
                            # We'll settle for the Phylum
                            elif len(taxonomic_lineage) > 3:
                                description = taxonomic_lineage[3]
                            else:
                                description = taxonomic_lineage[-1]
                        except ValueError:
                            logging.warning("Attempt to parse species from lineage failed for:\n" + lineage + "\n")
                            description = sequence_info.group(2)
                    else:
                        logging.error("TreeSAPP is unsure what to do with header '" + header + "'\n")
                        sys.exit(17)
            else:
                # This sequence is probably an uninformative contig name
                lineage = "unclassified sequences; "
                description = self.InputData.split('/')[-1]

            if not description or not lineage:
                logging.warning("Description is unavailable for sequence '" + header + "'\n")
            tree_taxa_string += num_id + "\t" + description + " | " + accession + "\t" + lineage + "\n"
            self.master_reference_index[num_id].organism = organism
            self.master_reference_index[num_id].description = description
            self.master_reference_index[num_id].accession = accession
            self.master_reference_index[num_id].lineage = lineage

        # Write the new TreeSAPP numerical IDs, descriptions and lineages
        tree_taxa_list = self.Output + "tax_ids_" + self.COG + ".txt"
        tree_tax_list_handle = open(tree_taxa_list, "w")
        tree_tax_list_handle.write(tree_taxa_string)
        tree_tax_list_handle.close()

        logging.info("done.\n")

        return

    def swap_tree_names(self, tree, tree_swap_name):
        """
        Function used for replacing unique identifiers in a NEWICK tree file
        :param tree: The tree with leaf names that need to be replaced
        :param tree_swap_name: Name of the output tree file
        :return:
        """
        try:
            old_tree = open(tree, 'r')
        except IOError:
            raise IOError("Unable to open " + tree + " for reading!")
        try:
            new_tree = open(tree_swap_name, 'w')
        except IOError:
            raise IOError("Unable to open " + tree_swap_name + " for writing!")

        newick_tree = old_tree.readlines()
        old_tree.close()

        if len(newick_tree) > 1:
            raise AssertionError("ERROR: " + tree + " should only contain one line of text to be a NEWICK tree!")
        else:
            newick_tree = str(newick_tree[0])

        new_tree.write(re.sub('_' + self.COG, '', newick_tree) + "\n")
        new_tree.close()

        return

    def align_multiple_sequences(self, unaligned_ref_seqs, args):
        """
        Call MAFFT to perform a multiple sequence alignment of the reference sequences and the
        gene sequences identified by TreeSAPP

        :param args: Command-line argument object from get_options and check_parser_arguments
        :param unaligned_ref_seqs:
        :return: Name of the FASTA file containing the MSA
        """

        if args.verbose:
            logging.info("Aligning the reference and identified " + self.COG + " sequences using MAFFT... ")

        # Combine the reference and candidate sequence dictionaries
        unaligned_ref_seqs.update(self.ContigDict)
        ref_unaligned = self.Output + self.COG + "_gap_removed.fa"
        write_new_fasta(unaligned_ref_seqs, ref_unaligned)

        aligned_fasta = self.Output + self.COG + "_d_aligned.fasta"

        mafft_align_command = [args.executables["mafft"]]
        mafft_align_command += ["--maxiterate", str(1000)]
        mafft_align_command += ["--thread", str(args.num_threads)]
        if len(self.ContigDict) > 700:
            mafft_align_command.append("--auto")
        else:
            mafft_align_command.append("--localpair")
        mafft_align_command += [ref_unaligned, '1>' + aligned_fasta]
        mafft_align_command += ["2>", "/dev/null"]

        stdout, mafft_proc_returncode = launch_write_command(mafft_align_command, False)

        if mafft_proc_returncode != 0:
            logging.error("Multiple sequence alignment using " + args.executables["mafft"] +
                          " did not complete successfully! Command used:\n" + ' '.join(mafft_align_command) + "\n")
            sys.exit(17)

        logging.info("done.\n")

        return aligned_fasta

    def execute_raxml(self, phylip_file, raxml_destination_folder, args):
        os.makedirs(raxml_destination_folder)
        # No difference between this command and that in create_treesapp_ref_data
        raxml_command = [args.executables["raxmlHPC"], '-m', self.raxml_model]
        # Run RAxML using multiple threads, if CPUs available
        raxml_command += ['-T', str(int(args.num_threads))]
        if args.bootstraps == 0:
            nboot = "autoMR"
        else:
            nboot = str(args.bootstraps)
        raxml_command += ['-s', phylip_file,
                          '-f', 'a',
                          '-x', str(12345),
                          '-p', str(12345),
                          '-#', nboot,
                          '-n', self.COG,
                          '-w', raxml_destination_folder]

        logging.debug("RAxML command:\n\t" + ' '.join(raxml_command) + "\n")
        logging.info("Building Maximum-Likelihood tree with RAxML")

        raxml_pro = subprocess.Popen(' '.join(raxml_command), shell=True, preexec_fn=os.setsid)
        raxml_pro.wait()

        # logging.info("done.\n")

        return

