#!/usr/bin/python

from optparse import OptionParser
from time import gmtime, strftime

import create_trees_utility
import string
import os
import re
import sys


class FuncTreeUtility(create_trees_utility.CreateTrees):

    def __init__(self, inputData, ref_data_chosen, RAxML_out, var_out):
        create_trees_utility.CreateTrees.__init__(self, inputData, ref_data_chosen, RAxML_out, var_out)

if __name__ == '__main__':

    parser = OptionParser(usage='%prog ...')
    parser.add_option("-m", "--mltreemap_output",
                      help="MLTreeMap output")
    parser.add_option("-o", "--output",
                      help="Output directory")
    parser.add_option("-n", "--names",
                      help="Your input NAMES file for sequences to be added (NOT the reference sequences)")
    parser.add_option("-d", "--data_set",
                      help="Target reference data set (refer to the second column of cog_list.txt)")
    parser.add_option("-v", "--version",
                      help="Version of MLTreeMap used (N for Nucleotide or P for Peptide)")
    parser.add_option("-u", "--uclust",
                      help="Run uclust? (y or n)")
    parser.add_option("-g", "--gap_removal",
                      help="Remove minor gaps using Gblocks? (y or n)")
    parser.add_option("-i", "--identity",
                      help="Sequence identity value to be used in uclust "
                           "(in decimal value - e.g. 0.97 would be 97% identify")
    parser.add_option("-a", "--alignment_mode",
                      help="Alignment mode: \"d\" for default and \"p\" for profile-profile alignment")

    # Note:
    # For Names files, here is the required format

    # Column_1 Column_2 (two columns are separated by tab)
    # Column_1: Name of the sequence used in the MLTreeMap analysis and detected in the specific functional marker data
    # Column_2: This can be pretty much anything that describes the sequence and even the same name as in Column_1

    # parser.add_option("-r","--retrieval",help="Is program running on retrieval mode? (Y or N)")
    # parser.add_option("-l","--list",help="If running on retrieval mode, provide a list of sequences to be retrieved")

    (options, args) = parser.parse_args()

    inputFolder = options.mltreemap_output
    outputFolder = options.output
    inputNames = options.names
    inputData = options.data_set
    version = options.version
    runUclust = options.uclust

    RAxML_outputs = inputFolder + "/final_RAxML_outputs"
    various_outputs = inputFolder + "/various_outputs"
    other_data = "other_data"

    contig_aa_dict = {}

    listRAxML_Files = os.listdir(RAxML_outputs)

    cog_list = other_data + "/" + "cog_list.txt"

    cog_list_handle = open(cog_list, "rb")

    cog_list_lines = cog_list_handle.readlines()

    ref_data_chosen = ""

    for each_cog_list_line in cog_list_lines:
        each_cog_list_line = string.strip(each_cog_list_line)

        func_phylo_match = re.match("(\S+)\t(\S+)\t(\S+(\s+\S+)*)", each_cog_list_line)

        if func_phylo_match:
            refName = func_phylo_match.group(1)
            refCode = func_phylo_match.group(2)
            refDescrpition = func_phylo_match.group(3)

            if refCode == inputData:
                ref_data_chosen = refName

    cog_list_handle.close()
    # inputData is the second column of cog_list.txt
    # ref_data_chosen is the first column of cog_list.txt
    tree_util = FuncTreeUtility(inputData, ref_data_chosen, RAxML_outputs, various_outputs)

    print "**************** Reading RAxML files to generate list of contigs to be aligned *****************\n"

    acc_list_map = {}

    names_handle = open(inputNames, "rb")
    names_lines = names_handle.readlines()

    # Pulls the contig names from names file
    for each_names_ln in names_lines:
        each_names_ln = string.strip(each_names_ln)
        two_cols_match = re.match("(\S+)\t(\S+)", each_names_ln)

        if two_cols_match:
            each_na_acc = two_cols_match.group(1)
            each_na_desc = two_cols_match.group(2)
            acc_list_map[each_na_acc] = each_na_desc

    # Gets the
    for raxml_fileName in listRAxML_Files:
        # check to see whether the input reference data selected is one of 40 COG groups:
        cog_input_data = re.match("COG\d+", inputData)
        # check to see whether the input reference data selected is one of 40 COG groups:
        geba_input_data = re.match("g_COG\d+", inputData)

        if cog_input_data:
            fileName_match = re.match("p_(\S+)_RAxML_parsed.txt", raxml_fileName)
        elif geba_input_data:
            fileName_match = re.match("g_(\S+)_RAxML_parsed.txt", raxml_fileName)
        else:
            fileName_match = re.match("%s_(\S+)_RAxML_parsed.txt" % inputData, raxml_fileName)
        
        if fileName_match:
            contig = fileName_match.group(1)
                    
            if contig in acc_list_map.keys():
                # print "Processing following sequences: ", contig
                contig_aa_dict[contig] = acc_list_map[contig]

    list_various_files = os.listdir(various_outputs)
    
    # Writing input FASTA file #
    print "**************** Generating FASTA files of sequences to be aligned *****************\n"
            
    new_align_file = outputFolder + "_" + ref_data_chosen + "_unaligned.fasta"
    
    tree_util.write_query_fasta(list_various_files, contig_aa_dict, new_align_file, version)
    
    # Write the reference NAMES file #
    
    tree_util.write_reference_names(outputFolder)
 
    # Generate the finalized NAMEs file. If UCLUST was executed, only include the names of the cluster reprsentatives
    
    if runUclust == "y":
        
        uclustID = options.identity
        tree_util.run_uclust(new_align_file, outputFolder, uclustID)
        
        uclustID_fasta = "uclust_" + outputFolder + ".fasta"
        
        cluster_rep_dict = tree_util.cluster_reps(uclustID_fasta)
        
        inputNames_handle = open(inputNames, "rb")
        inputNames_lines = inputNames_handle.readlines()
        
        select_names = "selected_" + inputNames
        select_names_handle = open(select_names, "w")
        
        for each_names_line in inputNames_lines:
            each_names_line = string.strip(each_names_line)
            
            acc_and_desc_match = re.match("(\S+)\t(\S+)", each_names_line)
            
            if acc_and_desc_match:
                acc = acc_and_desc_match.group(1)
                desc = acc_and_desc_match.group(2)
                
                if acc in cluster_rep_dict.keys():
                    select_write = acc + "\t" + desc
                    select_names_handle.write(select_write)
                    select_names_handle.write("\n")
       
        select_names_handle.close()
        inputNames_handle.close()
        
        cat_command = 'cat %s_%s_ref.names %s > %s_%s_concat.names' % \
                      (outputFolder, ref_data_chosen, select_names, outputFolder, ref_data_chosen)
        print cat_command
        os.system(cat_command)
    else:
        select_names = "selected_" + inputNames
        select_names_handle = open(select_names, "w")
        
        for acc_key in contig_aa_dict.keys():
            desc = contig_aa_dict[acc_key]
            select_write = acc_key + "\t" + desc
            select_names_handle.write(select_write)
            select_names_handle.write("\n")
        
        select_names_handle.close()
        
        cat_command = 'cat %s_%s_ref.names %s > %s_%s_concat.names' % \
                      (outputFolder, ref_data_chosen, select_names, outputFolder, ref_data_chosen)        
        print cat_command
        os.system(cat_command)

    print "Final names file generated: %s_%s_concat.names" % (outputFolder, ref_data_chosen)

    # Align the sequences together #
    
    align_option = options.alignment_mode

    # check and see if we need to use COG alignments from GEBA folder:
    GEBA_ref_match = re.match("g_COG(\d+)", ref_data_chosen)
    
    if GEBA_ref_match:
        COG_number = GEBA_ref_match.group(1)
        COG_alignment = "COG" + COG_number

        ref_align = "geba_alignment_data/" + COG_alignment + ".fa"
    else:
        ref_align = "alignment_data/" + ref_data_chosen + ".fa"
        
    aligned_fasta = ""
    
    tree_util.align_sequences(align_option, runUclust, outputFolder, ref_align)
    
    if align_option == "d":
        aligned_fasta = outputFolder + "_" + ref_data_chosen + "_d_aligned.fasta"
    elif align_option == "p":
        aligned_fasta = outputFolder + "_" + ref_data_chosen + "_p_aligned.fasta"
    
    remove_gaps = options.gap_removal
    
    phylip_file = outputFolder + "_" + "%s.phy" % ref_data_chosen
    
    fasta_random = outputFolder + "_" + ref_data_chosen + "_concat_rfive.fasta"

    # Execute RAxML #
    if remove_gaps == "y":
        tree_util.execute_gblocks(fasta_random)
        os.system('cp %s-gb %s' % (fasta_random, fasta_random))
        
        os.system('java -cp readseq.jar run -a -f=12 %s' % fasta_random)
        
        os.system('mv %s.phylip %s' % (fasta_random, phylip_file))
    else:
        os.system('java -cp readseq.jar run -a -f=12 %s' % fasta_random)
        
        os.system('mv %s.phylip %s' % (fasta_random, phylip_file))
    
    tree_util.execute_RAxML(phylip_file, outputFolder)
    
    # Organize Output Files #
    
    project_folder = outputFolder + "_"
    time = strftime("%d_%b_%Y_%H_%M", gmtime())
        
    project_folder += str(time)
        
    os.system('mkdir %s' % project_folder)
    
    os.system('mkdir alignment_files')
    os.system('mv %s %s alignment_files' % (aligned_fasta, phylip_file))
    os.system('mv alignment_files %s' % project_folder)
    
    RAxML_destination_folder = outputFolder + "phy_files_%s_" % ref_data_chosen
    
    bestTree_nameswap = outputFolder + "_" + ref_data_chosen + "_best.tree"
    bootstrap_nameswap = outputFolder + "_" + ref_data_chosen + "_bootstrap.tree"
        
    os.system('mkdir final_tree_files')
    os.system('mv %s %s final_tree_files' % (bestTree_nameswap, bootstrap_nameswap))
    # Move the final bootstrap, non-bootstrap best trees into destination folder:
    os.system('mv final_tree_files %s' % project_folder)
    # Move the RAxML output folder into destination folder:
    os.system('mv %s %s' % (RAxML_destination_folder, project_folder))
    
    prefix = outputFolder + "_" + ref_data_chosen
    
    os.system('mv %s* %s' % (prefix, project_folder))
    
    if runUclust == "y":
        os.system('mkdir %s_uclust' % outputFolder)
    
        os.system('mv uclust_%s %s_uclust' % (outputFolder, outputFolder))
        os.system('mv usort_%s %s_uclust' % (outputFolder, outputFolder))
    
        os.system('mv %s_uclust %s' % (outputFolder, project_folder))
