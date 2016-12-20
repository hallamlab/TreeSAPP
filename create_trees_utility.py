#!/usr/bin/python

import string
import os
import re
import sys
import random


class CreateTrees:

    def __init__(self, inputData, refData_chosen, RAxML_out, var_out):
        # inputData = second column of cog_list.txt  (e.g. mcrA)
        # refData_chosen = first column of cog_list.txt (e.g. mcrAref)
        self.inputData = inputData
        self.refData_chosen = refData_chosen
        self.RAxML_out = RAxML_out
        self.var_out = var_out
        
    def write_query_fasta(self, list_various_files, aa_dictionary, query_align_file, input_molecular_type):
        
        print ">>>>> Writing to file ", query_align_file
        query_align_file_handle = open(query_align_file, "w")
        
        cog_id = ""
        
        for varFile in list_various_files:
            COG_ref_match = re.match("\S*(COG\d+)", self.refData_chosen)
            
            if input_molecular_type == "P":
                
                if COG_ref_match:
                    cog_id = COG_ref_match.group(1)
                    
                # predicted_GEBA_ORF = re.match("g_(\S+)_%s(_\d+)*.fa" % (cog_id), varFile)
                predicted_ORF = re.match("%s_(\S+)_%s(_\d+)*.fa" % (self.inputData, self.refData_chosen), varFile)
                predicted_COG_ORF = re.match("g_(\S+)_%s(_\d+)*.fa" % cog_id, varFile)
                
                if predicted_ORF:
                    match_contig = predicted_ORF.group(1)
                    
                    if aa_dictionary.has_key(match_contig):
                        query_align_file_handle.write(">")
                        query_align_file_handle.write(match_contig)
                        query_align_file_handle.write("\n")
                
                        orf_file_to_open = self.var_out + "/" + varFile
                        orf_file_handle = open(orf_file_to_open, "rb")
                
                        orf_file_lines = orf_file_handle.readlines()
                
                        for each_orf_line in orf_file_lines:
                            each_orf_line = string.strip(each_orf_line)
                            # fields = string.split(each_orf_line,"\t")
                    
                            # sequence = fields[4]
                            match_fasta = re.match("([A-Za-z]+)",each_orf_line)
                    
                            if match_fasta:
                                sequence = match_fasta.group(1)
                                query_align_file_handle.write(sequence)
                                query_align_file_handle.write("\n")
                    
                        orf_file_handle.close()
                
                elif predicted_COG_ORF:
                    match_contig = predicted_COG_ORF.group(1)
                    
                    if aa_dictionary.has_key(match_contig):
                        query_align_file_handle.write(">")
                        query_align_file_handle.write(match_contig)
                        query_align_file_handle.write("\n")
                
                        orf_file_to_open = self.var_out + "/" + varFile
                        orf_file_handle = open(orf_file_to_open, "rb")
                
                        orf_file_lines = orf_file_handle.readlines()
                
                        for each_orf_line in orf_file_lines:
                            each_orf_line = string.strip(each_orf_line)
                            # fields = string.split(each_orf_line,"\t")
                    
                            # sequence = fields[4]
                            match_fasta = re.match("([A-Za-z]+)",each_orf_line)
                    
                            if match_fasta:
                                sequence = match_fasta.group(1)
                                query_align_file_handle.write(sequence)
                                query_align_file_handle.write("\n")
                    
                        orf_file_handle.close()
            elif input_molecular_type == "N":
                predicted_ORF = re.match("%s_(\S+)_%s(_\d+)*.fa" % (self.inputData, self.refData_chosen), varFile)
                predicted_COG_ORF = re.match("p_(\S+)_%s(_\d+)*.fa" % self.refData_chosen, varFile)
                predicted_GEBA_ORF = re.match("g_(\S+)_%s(_\d+)*.fa" % self.refData_chosen, varFile)
            
                if predicted_ORF:
                    match_contig = predicted_ORF.group(1)
                    if aa_dictionary.has_key(match_contig):
                    
                        query_align_file_handle.write(">")
                        query_align_file_handle.write(match_contig)
                        query_align_file_handle.write("\n")
                    
                        orf_file_to_open = self.var_out + "/" + varFile
                        orf_file_handle = open(orf_file_to_open, "rb")
                    
                        orf_file_handle.readline()
                    
                        sequence = orf_file_handle.readline()
                        sequence = string.strip(sequence)
                    
                        sequence = re.sub("[-\.*]","",sequence)

                        query_align_file_handle.write(sequence)
                        query_align_file_handle.write("\n")
                
                        orf_file_handle.close()
        
                elif predicted_COG_ORF or predicted_GEBA_ORF:
                    match_contig = predicted_COG_ORF.group(1)
                    if aa_dictionary.has_key(match_contig):
                    
                        query_align_file_handle.write(">")
                        query_align_file_handle.write(match_contig)
                        query_align_file_handle.write("\n")
                    
                        orf_file_to_open = self.var_out + "/" + varFile
                        orf_file_handle = open(orf_file_to_open, "rb")
                    
                        orf_file_handle.readline()
                    
                        sequence = orf_file_handle.readline()
                        sequence = string.strip(sequence)
                    
                        sequence = re.sub("[-\.*]", "", sequence)

                        query_align_file_handle.write(sequence)
                        query_align_file_handle.write("\n")
                
                        orf_file_handle.close()       
        query_align_file_handle.close()

    # Runs UCLUST #
    def run_uclust(self, query_align_file, output_folder, uclustID):
        
        print "******************* Executing UCLUST with identity value of %s per cent *******************" % uclustID
        
        usort_input = query_align_file
        usort_output = "usort_" + output_folder
        usort_log = "usort_" + output_folder + ".log"
        
        uclust_sort_command = "/Applications/usearch --sort %s --output %s  --log %s" % \
                              (usort_input, usort_output, usort_log)
        
        print uclust_sort_command
        
        os.system(uclust_sort_command)
        
        uclust_input = usort_output
        uclust_output = "uclust_" + output_folder + ".fasta"
        uclust_uc = "uclust_" + output_folder + ".uc"
        uclust_log = "uclust_" + output_folder + ".log"
        
        uclust_cluster_command = "/Applications/usearch --cluster %s --id %s --seedsout %s --uc %s --log %s" % \
                                 (uclust_input, uclustID, uclust_output, uclust_uc, uclust_log)
        
        print uclust_cluster_command
        
        os.system(uclust_cluster_command)
    
    # Reading the FASTA file of representative sequences of each cluster and generating the list of sequences #
    def cluster_reps(self, uclust_output):    
        
        uclust_fasta_handle = open(uclust_output, "rb")
        
        uclust_fasta_lines = uclust_fasta_handle.readlines()
        
        uc_list_map = {}
        
        for each_uc_fasta_ln in uclust_fasta_lines:
            each_uc_fasta_ln = string.strip(each_uc_fasta_ln)
            
            fasta_tag_match = re.match(">(\S+)", each_uc_fasta_ln)
            
            if fasta_tag_match:
                each_na_acc = fasta_tag_match.group(1)
                uc_list_map[each_na_acc] = 1
        
        return uc_list_map

    # subroutine write_reference_names: Generate the NAMEs file for all COG reference groups #
    def write_reference_names(self, output_folder):
        tax_id_map = {}
    
        ref_names = output_folder + "_" + "%s_ref.names" % (self.refData_chosen)
        ref_names_handle = open(ref_names, "w")
    
        # check to see whether we need to use the COG alignment files from GEBA
        
        GEBA_ref_match = re.match("g_COG(\d+)",self.refData_chosen)
        
        if GEBA_ref_match:
            cog_number = GEBA_ref_match.group(1)
            cog_id = "COG" + cog_number
            ref_alignment_handle = open("geba_alignment_data/%s.fa" % cog_id, "rb")
            
            ref_align_lines = ref_alignment_handle.readlines()
    
            for each_ref_align_line in ref_align_lines:
                each_ref_align_line = string.strip(each_ref_align_line)
                header_match = re.match("^>(\d+)_%s" % (cog_id), each_ref_align_line)
        
                if header_match:
                    tax_id = header_match.group(1)
                    header_trimmed = re.sub("^>", "", each_ref_align_line)
                    tax_id_map[tax_id] = header_trimmed
    
            ref_alignment_handle.close()
        else:
            ref_alignment_handle = open("alignment_data/%s.fa" % self.refData_chosen, "rb")
    
            ref_align_lines = ref_alignment_handle.readlines()
    
            for each_ref_align_line in ref_align_lines:
                each_ref_align_line = string.strip(each_ref_align_line)
                header_match = re.match("^>(\d+)_%s" % self.refData_chosen, each_ref_align_line)
        
                if header_match:
                    tax_id = header_match.group(1)
                    header_trimmed = re.sub("^>", "", each_ref_align_line)
                    tax_id_map[tax_id] = header_trimmed
    
            ref_alignment_handle.close()
    
        # Handle tax ids for COG here #
        cog_input_match = re.match("COG\d+", self.refData_chosen)
        geba_cog_match  = re.match("g_COG\d+", self.refData_chosen)
        
        if cog_input_match:
            ref_tax_ids_handle = open("other_data/tax_ids_nr.txt","rb")
        elif geba_cog_match:
            ref_tax_ids_handle = open("other_data/tax_ids_geba_tree.txt","rb")
        else:
            ref_tax_ids_handle = open("other_data/tax_ids_%s.txt" % self.refData_chosen, "rb")
    
        ref_tax_ids_lines = ref_tax_ids_handle.readlines()
    
        for each_ref_tax_ids_line in ref_tax_ids_lines:
            each_ref_tax_ids_line = string.strip(each_ref_tax_ids_line)
        
            ids_desc_match = re.match("(\S+)\t(\S+(\s+\S+)*)",each_ref_tax_ids_line)
        
            if ids_desc_match:
                ids = ids_desc_match.group(1)
                description = ids_desc_match.group(2)
            
                if tax_id_map.has_key(ids):
                    write_names = tax_id_map[ids] + "\t" + description
                    ref_names_handle.write(write_names)
                    ref_names_handle.write("\n")
    
        ref_tax_ids_handle.close()
        ref_names_handle.close()

    def write_unaligned_ref_fasta(self, output_folder, ref_align):
        ref_align_gap_removed = output_folder + "_" + self.refData_chosen + "_gap_removed.fa"
    
        ref_align_handle = open(ref_align, "rb")
        ref_align_gap_rm_handle = open(ref_align_gap_removed, "w")
        
        first_fas_line = ref_align_handle.readline()
        first_fas_line = string.strip(first_fas_line)
        
        first_header_match = re.match("^>", first_fas_line)
        
        if first_header_match:
            ref_align_gap_rm_handle.write(first_fas_line + "\n")
        
        fasta_in_lines = ref_align_handle.readlines()
        
        alignment_gap_removed = ""
        
        for each_fas_line in fasta_in_lines:
            each_fas_line = string.strip(each_fas_line)
            
            fasta_header_match = re.match("^>", each_fas_line)
            
            if fasta_header_match:
                ref_align_gap_rm_handle.write(alignment_gap_removed + "\n")
                ref_align_gap_rm_handle.write(each_fas_line + "\n")
                
                alignment_gap_removed = ""
            else:
                    
                alignment_gap_removed += each_fas_line
        
                if re.search("[\-]+", alignment_gap_removed):
                    alignment_gap_removed = re.sub("-", "", alignment_gap_removed)
            
        ref_align_gap_rm_handle.write(alignment_gap_removed + "\n")

        ref_align_gap_rm_handle.close()
        ref_align_handle.close()
    
    # TODO: scan the unaligned reference sequences for lines of X's (like in case of COGs)
    
    def scan_unaligned_ref_fasta(self,  output_folder, ref_align_gap_removed):    
        ref_align_handle = open(ref_align_gap_removed, "rb")
        ref_align_gap_rm_scan = output_folder + "_" + self.refData_chosen + "_gap_rm_scan.fa"
        
        ref_align_scan_handle = open(ref_align_gap_rm_scan, "w")
        
        fasta_map = {}
        
        fasta_in_lines = ref_align_handle.readlines()
        
        sequence_id = ""
        sequence = ""
        
        for each_fas_line in fasta_in_lines:
            each_fas_line = string.strip(each_fas_line)
            
            fasta_header_match = re.match("^>(\S+)", each_fas_line)
            
            if fasta_header_match:
                sequence_id = each_fas_line
            else:
                sequence = each_fas_line
                
            fasta_map[sequence_id] = sequence
        
        for each_sequence_id in fasta_map.keys():
            each_sequence = fasta_map[each_sequence_id]
            
            line_of_X_s = re.match("^X((X)+)*$", each_sequence)
            
            if not line_of_X_s:
                ref_align_scan_handle.write(each_sequence_id + "\n")
                ref_align_scan_handle.write(each_sequence+"\n")
        
        ref_align_handle.close()
        ref_align_scan_handle.close()
    
    def create_random_names(self, random_map, output_folder):
        orig_id_desc_dict = {}
        concat_names = output_folder + "_" + self.refData_chosen + "_concat.names"
        
        concat_names_handle = open(concat_names,"rb")
        
        concat_names_lines = concat_names_handle.readlines()
        
        for each_cnl in concat_names_lines:
            each_cnl = string.strip(each_cnl)
            
            tab_separated_match = re.match("(\S+)\t(.*)$", each_cnl)
            
            if tab_separated_match:
                orig_id = tab_separated_match.group(1)
                desc = tab_separated_match.group(2)
                
                orig_id_desc_dict[orig_id] = desc
        
        concat_names_handle.close()
        
        concat_random_names = output_folder + "_" + self.refData_chosen + "_concat_rand.names"
        concat_rand_names_handle = open(concat_random_names, "w")
        
        for original_id in random_map.keys():
            names_line = random_map[original_id] + "\t" + orig_id_desc_dict[original_id]
            concat_rand_names_handle.write(names_line)
            concat_rand_names_handle.write("\n")
        
        concat_rand_names_handle.close()
    
    def randomize_fasta_id(self, fasta, output_folder):
        original_random_dict = {}
        
        fasta_handle = open(fasta, "rb")
        fasta_random = output_folder + "_" + self.refData_chosen + "_concat_rfive.fasta"
        
        fasta_random_handle = open(fasta_random, "w")
        
        fasta_lines = fasta_handle.readlines()
        
        for each_fasta_line in fasta_lines:
            each_fasta_line = string.strip(each_fasta_line)
            header_match = re.match("^>(\S+)", each_fasta_line)
            
            if header_match:
                random_five = "ID"
                random_five += ''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(5))
                
                original_id = header_match.group(1)
                replaced_header = re.sub(original_id, random_five, each_fasta_line)
                original_random_dict[original_id] = random_five
                
                fasta_random_handle.write(replaced_header)
                fasta_random_handle.write("\n")
            else:
                fasta_random_handle.write(each_fasta_line)
                fasta_random_handle.write("\n")
        
        fasta_random_handle.close()
        fasta_handle.close()
        
        self.create_random_names(original_random_dict,output_folder)

    def align_sequences(self, align_option, runUclust, output_folder, ref_align):
        
        query_fasta = ""
        
        ### Default alignment ###
        if align_option == "d":
            self.write_unaligned_ref_fasta(output_folder, ref_align)
            ref_align_gap_removed = output_folder + "_" + self.refData_chosen + "_gap_removed.fa"
            
            self.scan_unaligned_ref_fasta(output_folder, ref_align_gap_removed)
            ref_align_gap_rm_scan = output_folder + "_" + self.refData_chosen + "_gap_rm_scan.fa"
            
            concat_fasta = output_folder + "_" + self.refData_chosen + "_concat.fasta"
        
            if runUclust == "y":
                query_fasta = "uclust_" + output_folder + ".fasta"
            else:
                query_fasta = output_folder + "_" + self.refData_chosen + "_unaligned.fasta"
            
            os.system('cat %s %s > %s' % (query_fasta, ref_align_gap_rm_scan , concat_fasta))
            
            aligned_fasta = output_folder + "_" + self.refData_chosen + "_d_aligned.fasta"
            
            fasta_random = output_folder + "_" + self.refData_chosen + "_concat_rfive.fasta"
            
            muscle_align_command = "muscle -in %s -out %s" % (concat_fasta, aligned_fasta)
 
            print muscle_align_command,"\n"
            os.system(muscle_align_command)
            
            self.randomize_fasta_id(aligned_fasta,output_folder)
        
        # Profile-Profile alignment #
        elif align_option == "p":
            
            if runUclust == "y":
                query_fasta = "uclust_" + output_folder + ".fasta"
            else:
                query_fasta = output_folder + "_" + self.refData_chosen + "_unaligned.fasta"
            
            query_align = output_folder + "_" + self.refData_chosen + "_query_aligned.fasta"
            
            muscle_align_command = "muscle -in %s -out %s" % (query_fasta, query_align)
            
            print muscle_align_command,"\n"
            os.system(muscle_align_command)
            
            aligned_fasta = output_folder + "_" + self.refData_chosen + "_p_aligned.fasta"
            muscle_align_command = "muscle -profile -in1 %s -in2 %s -out %s" % (query_align, ref_align, aligned_fasta)
            
            print muscle_align_command,"\n"
            os.system(muscle_align_command)
            
            self.randomize_fasta_id(aligned_fasta,output_folder)
        else:
            sys.exit("You need to specify the alignment method")

    # subroutine retrieve_data_size: return number of sequences in FASTA file #
    def retrieve_data_size(self, aligned_fasta):
        num_seqs = 0
        
        fasta_file_handle = open(aligned_fasta, "rb")
        
        fasta_lines = fasta_file_handle.readlines()
        
        for each_fa_line in fasta_lines:
            if re.search(">", each_fa_line):
                num_seqs += 1
        
        return num_seqs

    def execute_gblocks(self, aligned_fasta):
        data_size = self.retrieve_data_size(aligned_fasta) 
        min_flank_pos = str(0.55 * data_size)
        gblock_command = "sub_binaries/Gblocks %s -t=p -s=y -u=n -p=t -b3=15 -b4=3 -b5=h -b2=%s" \
                         % (aligned_fasta, min_flank_pos)
        
        print gblock_command, "\n"
        
        os.system(gblock_command)

    def execute_RAxML(self, phylip_file, output_folder):
        raxml_destination_folder = output_folder + "phy_files_%s_" % self.refData_chosen
        os.system('mkdir %s' % raxml_destination_folder)
        
        if self.inputData == "a":
            raxml_command = "sub_binaries/raxmlHPC-PTHREADS" \
                            " -f a -x 12345 -# 100 -m GTRGAMMA" \
                            " -s %s -n %s -w %s -T 8  -p 8" \
                            % (phylip_file, self.refData_chosen, raxml_destination_folder)
        else:
            raxml_command = "sub_binaries/raxmlHPC-PTHREADS " \
                            "-f a -x 12345 -# 100 -m PROTGAMMAWAG" \
                            " -s %s -n %s -w %s -T 8  -p 8" \
                            % (phylip_file, self.refData_chosen, raxml_destination_folder)
    
        print raxml_command, "\n"
        
        os.system(raxml_command)
        
        bestTree = raxml_destination_folder + "/RAxML_bestTree." + self.refData_chosen
        bootstrapTree = raxml_destination_folder + "/RAxML_bipartitions." + self.refData_chosen
    
        bestTree_nameswap = output_folder + "_" + self.refData_chosen + "_best.tree"
        bootstrap_nameswap = output_folder + "_" + self.refData_chosen + "_bootstrap.tree"
        
        concat_random_names = output_folder + "_" + self.refData_chosen + "_concat_rand.names"
        
        os.system('swapTreeNames.pl -t %s -l %s -o %s' % (bestTree, concat_random_names, bestTree_nameswap))
        os.system('swapTreeNames.pl -t %s -l %s -o %s' % (bootstrapTree, concat_random_names, bootstrap_nameswap))
