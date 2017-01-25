#!/usr/bin/python

from optparse import OptionParser
from time import gmtime, strftime
from decimal import *

#import create_trees_utility
import string
import os
import shutil
import re
import sys
import random

#class func_tree_utility(create_trees_utility.create_trees):

#    def __init__(self, inputData, refData_chosen, RAxML_out, var_out):
#        create_trees_utility.create_trees.__init__(self, inputData, refData_chosen, RAxML_out, var_out)
        
if __name__ == '__main__':

    parser = OptionParser(usage='%prog ...')
    parser.add_option("-f","--fasta_file",help="MLTreeMap output")
    parser.add_option("-t","--table",help="Tab-separated columns listing sequences to be used in reference data")
    parser.add_option("-r","--ref_name",help="Name of reference data (refer to the first column of 'cog_list.txt' under '#functional cogs' section)")
    
    (options,args) = parser.parse_args()

    inputFasta = options.fasta_file
    inputTable = options.table
    inputRefName = options.ref_name
    

    def create_new_fasta(in_fasta, out_fasta, dictionary, header_type):
        missing_counter = 0
        
        out_fasta_handle = open(out_fasta,"w")
    
        in_fasta_handle = open(in_fasta,"rb")
    
        first_line = in_fasta_handle.readline()
        first_line = string.strip(first_line)
    
        if header_type == "ncbi":
            header_match  = re.match(">gi\|(\d+)\|(\w+)\|(\S+(\.\d+)*)\|(.*)$",first_line)
    
            if header_match:
                accession = header_match.group(3)
                if dictionary.has_key(accession):
                   out_fasta_handle.write(">%s\n"%(dictionary[accession]))
                   missing_counter = 0
                else:
                   missing_counter = 1
        elif header_type == "mltree":
            header_match  = re.match(">(r_\d+)$",first_line)
            
            if header_match:
                id = header_match.group(1)
                
                if dictionary.has_key(id):
                    out_fasta_handle.write(">%s\n"%(dictionary[id]))
                    missing_counter = 0
                else:
                    missing_counter = 1
        fasta_lines = in_fasta_handle.readlines()
    
        sequences_write = ""
        
        for each_fa_line in fasta_lines:
            each_fa_line = string.strip(each_fa_line)
        
            #header_match = re.match(">gi\|(\d+)\|(\w+)\|(\S+\.\d+)\|(.*)$",each_fa_line)

            if header_type == "ncbi":
                header_match  = re.match(">gi\|(\d+)\|(\w+)\|(\S+(\.\d+)*)\|(.*)$",each_fa_line)
    
                if header_match:
                    accession = header_match.group(3)
                    if dictionary.has_key(accession):
                        if missing_counter == 0:
                            out_fasta_handle.write(sequences_write + "\n")
                        
                        out_fasta_handle.write(">%s\n"%(dictionary[accession]))
                        
                        missing_counter = 0
                        
                        sequences_write = ""
                    else:
                        if missing_counter == 0:
                            out_fasta_handle.write(sequences_write + "\n")
                        missing_counter = 1
                else:
                    sequences_write += each_fa_line
            elif header_type == "mltree":
                header_match  = re.match(">(r_\d+)$",each_fa_line)
            
                if header_match:
                    id = header_match.group(1)
        
                    if dictionary.has_key(id):
                        if missing_counter == 0:
                            out_fasta_handle.write(sequences_write + "\n")
                            
                        out_fasta_handle.write(">%s\n"%(dictionary[id]))
                        
                        missing_counter == 0
                        
                        sequences_write = ""
                    else:
                        if missing_counter == 0:
                            out_fasta_handle.write(sequences_write + "\n")
                        missing_counter = 1
                else:
                    sequences_write += each_fa_line
        
        if missing_counter == 0:    
            out_fasta_handle.write(sequences_write + "\n")
        
        out_fasta_handle.close()
        in_fasta_handle.close()
    
    ### End of subroutine create_new_fasta
    
    table_handle = open(inputTable,"rb")
    
    fasta_repl_dict = {}
    fasta_mltree_repl_dict = {}
    tree_repl_dict = {}
    mltree_dict = {}
    tree_name_dict = {}
    
    table_lines = table_handle.readlines()
    
    print "******************** Generating files for formatting purpose: ********************"
    
    for tbl_line in table_lines:
        tbl_line = string.strip(tbl_line)
        
        column_match = re.match("(\S+(\.\d+)*)\t(r_(\d+))\t(\d+_%s)\t(.+)$"%(inputRefName),tbl_line)
        
        if column_match:
            accession = column_match.group(1)
            short_id = column_match.group(3)
            tree_id = column_match.group(4)
            mltree_id = column_match.group(5)
            description = column_match.group(6)
    
            fasta_repl_dict[accession] = short_id
            fasta_mltree_repl_dict[short_id] = mltree_id
            tree_repl_dict[short_id] = description
            mltree_dict[short_id] =tree_id
            tree_name_dict[tree_id] = description
    
    table_handle.close()
    
    fasta_replace_names = inputRefName + "_fasta_replace.names"
    fasta_mltree_names = inputRefName + "_fasta_mltree.names"
    tree_replace_names = inputRefName + "_tree_replace.names"
    
    tree_taxa_list = "tax_ids_%s.txt" % (inputRefName)
    
    tree_tax_list_handle = open(tree_taxa_list,"w")
    
    for mltree_id_key in sorted(tree_name_dict.keys()):
        tree_tax_list_handle.write("%s\t%s\n"%(mltree_id_key, tree_name_dict[mltree_id_key]))
    tree_tax_list_handle.close()
    
    print "******************** %s generated ********************\n"%(tree_taxa_list)
    
    tree_names_list = "%s_tree_replace.names" % (inputRefName)
    
    tree_names_list_handle = open(tree_names_list,"w")
    
    for short_id_key in sorted(mltree_dict.keys()):
        tree_names_list_handle.write("%s:\t%s:\n"%(short_id_key, mltree_dict[short_id_key]))
    
    tree_names_list_handle.close()
    
    fasta_replaced = inputRefName + ".fc.repl.fasta"
    
    create_new_fasta(inputFasta, fasta_replaced, fasta_repl_dict, "ncbi")
    
    print "******************** %s generated ********************\n"%(tree_names_list)
    
    
    print "******************** FASTA file, %s generated ********************\n"%(fasta_replaced)
    
    print "******************** Aligning the sequences using MUSCLE ********************\n"
    
    fasta_replaced_align = inputRefName+".fc.repl.aligned.fasta"
    
    muscle_align_command = "muscle -in %s -out %s" % (fasta_replaced, fasta_replaced_align)
    
    print muscle_align_command,"\n"
    os.system(muscle_align_command)
    
    fasta_mltree = inputRefName + ".fa"
    
    create_new_fasta(fasta_replaced_align, fasta_mltree, fasta_mltree_repl_dict, "mltree")
    
    print "******************** FASTA file, %s generated ********************\n"%(fasta_mltree)
    
    makeblastdb_command = "makeblastdb -in %s -dbtype prot -input_type fasta -out %s" % (fasta_mltree, fasta_mltree)
    os.system(makeblastdb_command)
    
    print "******************** BLAST DB for %s generated ********************\n"% (inputRefName)
    
    hmm_build_command = "hmmbuild -s %s.hmm %s" % (inputRefName, fasta_mltree)
    os.system(hmm_build_command)
    
    print "******************** HMM file for %s generated ********************\n"% (inputRefName)    
    
    phylip_command = "java -cp readseq.jar run -a -f=12 %s" % (fasta_replaced_align)
    os.system(phylip_command)
    
    phylip_file = inputRefName + ".phy"
    os.system('mv %s.phylip %s' % (fasta_replaced_align, phylip_file))
    
    raxml_out = "%s_phy_files" % (inputRefName)
    os.system("mkdir %s" % (raxml_out))
    
    raxml_command = "raxmlHPC-PTHREADS -f a -x 12345 -# 100 -m PROTGAMMAWAG -s %s -n %s -w %s -T 4" % (phylip_file, inputRefName, raxml_out)
    os.system(raxml_command)
    
    
    tree_to_swap = "%s/RAxML_bestTree.%s" % (raxml_out, inputRefName)
    final_mltree = "%s_tree.txt" % (inputRefName)
    
    swapTree_command = "swapTreeNames.pl -t %s -l %s -o %s" % (tree_to_swap, tree_names_list, final_mltree)
    os.system(swapTree_command)
    
    
    final_output_folder = "MLTreeMap_files_%s" % (inputRefName)
    os.system("mkdir %s" % (final_output_folder))
    
    os.system("mv %s.fa %s.fa.p* %s" % (inputRefName, inputRefName, final_output_folder))
    os.system("mv %s.hmm %s %s %s" % (inputRefName, tree_taxa_list, final_mltree, final_output_folder))
    