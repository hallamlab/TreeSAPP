#!/usr/bin/python


##########                   mltreemap.py by Young C. Song                    ##########
##########                Originally developed by Manuel Stark                ##########
##########            Re-written on Young C. Song on August 3, 2013           ##########

##########     Program written to detect phylogenetic/functional              ##########
##########     marker genes from FASTA input files and phylogenetically       ##########
##########     place them in optimal loci of their respective reference tree. ##########



from optparse import OptionParser
from Bio import Entrez

import mltreemap_utility
import string
import os
import shutil
import re
import sys
import random

class mltree_utility(mltreemap_utility.mltree_workflow):

    def __init__(self, alignData, hmmData, treeData):
        mltreemap_utility.mltree_workflow.__init__(self, alignData, hmmData, treeData) 

if __name__ == '__main__':

    parser = OptionParser(usage='%prog ...')
    parser.add_option("-f", "--fasta", help="FASTA file")
    parser.add_option("-t", "--seq_type", help="Sequence type: Nucleotide (n) or Protein (p): (default - n)")
    parser.add_option("-r", "--ref_search", help="Search all reference data: all (a) or selected set (s): (default - a)")
    parser.add_option("-p", "--phylo_tree", help="Phylogenetic tree: MLTreeMap reference tree (m); GEBA tree (g) or Fungi tree (f): (default - m)")
    parser.add_option("-b", "--blast_bit", help="Minimum bitscore for BLAST: (default 60)")
    #parser.add_option("-c", "--raxml_cutoff", help="Minimum cutoff of amino acid lengths for RAxML: (default 50)")
    #parser.add_option("-a", "--algorithm", help="RAxML algorith: Maximum Likelihood (m) or Maximum Parsimony (v): (default m) ")
    parser.add_option("-o", "--output",help="Output directory")
    
    (options,args) = parser.parse_args()
    
    inputFASTA = options.fasta
    seqTypeChosen = options.seq_type
    
    outDir = options.output
    
    logFile = "log.txt"
    logHandle = open(logFile, 'w')
    
    newDir = ''
    outDirExists = 0
    
    if os.path.exists(outDir):
        print "WARNING: Directory %s already exists:"
        overWrite = raw_input("Please select the following ('o' for over-write; 'n' for new directory; and 'q' for exit): ") # In python 3, raw_input is replaced by input
        
        if overWrite is "o":
            os.system("rm -r %s" % (outDir))
            os.system("mkdir %s" % (outDir))
        elif overWrite is "n":
            newDir = raw_input("Please insert the name of new output directory: ")
            outDir = newDir
            os.system("mkdir %s" % (outDir))
            outDirExists = 1
        elif overWrite is "q":
            print "Exising program ..."
            sys.exit()
        else:
            print "You have entered the wrong option."
            print "Exiting program ..."
            sys.exit()
    else:
        os.system("mkdir %s" % outDir)
    
    #################### Part A: Check the input file and make sure it is in correct format                                                  ####################
    #################### Part A1 (Lines 64 to 91): Check the input file and make sure it is in correct format                                ####################
    #################### Part A2 (Lines 93 to 117): Insert the sequences from input FASTA into a dictionary for further processing           ####################
    #################### Part A3 (Lines 119 to 404): Sequence format check (illegal and repleceables charactes in ID, mol. type consistency) ####################
    
    print "************************************ Checking Input File for Format ************************************"
    
    #################### Part A1: Check and see if the input file is indeed a FASTA file ####################
    seqCount = 0 # Counting the number of sequences for splitting purpose
    
    logHandle.write("** Checking the input file, %s for format:\n" % (inputFASTA))
    
    fasSeqDict = {} # Contain the sequence ID and sequence for format checking
 
    fastaCheckHandle = open(inputFASTA, "rb")
    
    isHeader = 0 # Checking to see whether the line contains sequence ID
    
    # Let's check the first line and see if this starts with '>' (we are checking to see whether the input is a FASTA file or not)
    headerString = '' # If the first line starts with '>', we assume that input is a FASTA file and start adding our seq. ID and seq. String into a dictionary
    seqString = ''
    
    firstLine = fastaCheckHandle.readline()
    firstLine = string.strip(firstLine)
    
    headerMatch = re.match("(^>\s*\S+(\s+\S+)*$)", firstLine)
    
    if headerMatch:  # Sequences is probably a FASTA file, although we still need to do some more checking later
        headerString = headerMatch.group(1)
        seqCount += 1
    else:
        print "The first line of %s does not start with \'>\'." % (inputFASTA)
        print "Are you sure the input is a FASTA file?"
        sys.exit()
    #################### End of Part A1 ####################
    
    #################### Part A2: Insert the sequences from input fasta into a dictionary for further processing ####################
    # Now, reading the rest of the FASTA files and processing them accordingly.
    fastaLines = fastaCheckHandle.readlines()
    
    for i, eachFasLine in enumerate(fastaLines):
        eachFasLine = string.strip(eachFasLine)
        
        headerMatch = re.match("^(>\s*\S+(\s+\S+)*$)", eachFasLine)
        
        if headerMatch:
            if headerString in fasSeqDict.keys(): # Need to check if there are sequences with duplicate ID
                logHandle.write("\tERROR: A sequence with ID, %s already exists.\n" % (headerMatch.group(1)))
                logHandle.write("\tMake sure you don't have sequences with duplicate ID.\n")
                
                print "An ERROR was reported during the input file formatting procedure."
                print "Please review 'log.txt' for details."
                sys.exit()
            else:
                fasSeqDict[headerString] = seqString
                headerString = headerMatch.group(1)
                seqString = ''
                seqCount += 1
        else:
            seqString += eachFasLine
        
        if i == len(fastaLines) - 1:  # If we hit the last line, then we insert the last sequence into the dictionary
            fasSeqDict[headerString] = seqString
    #################### End of Part A2 ####################
   
    #################### Part A3: Sequence format check ####################
    # Now that we have the sequences in the dictionary, we can check the formatting quality for each sequence.
    # We also make some adjustments to the sequence ID and generate a new dictionary containing the adjusted sequence IDs
    
    finalSeqDict = {} # This is going to be a final list of sequences that will be analyzed subsequently after filtering.
    
    for eachId in sorted(fasSeqDict.keys()):
        logHandle.write("\t* Checking %s for correct ID format:\n" % (eachId))
        
        identicalSeqs = 0 # This is a switch that indicates whether sequence identical to that we are currently processing is present in the final list of sequences.
                          # Refer to lines 166 to 178 for detailed usage of this switch. 
        
        # Check the length of the string.  If it is greater than 50, we report and exit the program
        # If it is shorter than 50 characters, we move on and check for presence of illegal characters and "replaceable" characters
        
        strLength = len(eachId)
        
        if strLength > 50:
            logHandle.write("\t\tERROR: Sequence ID %s has string length greater than 50\n" % (eachId))
            logHandle.write("\t\tPlease reduce the length to 50 or less.\n")
            logHandle.write("\t\tSequence description should be written on the .names files\n")
            
            print "An ERROR was reported during the input file formatting procedure."
            print "Please review 'log.txt' for details."
            
            sys.exit()
        else:
            # If the length of sequence ID is less than equal to 50 than we check the ID strings for presence of illegal characters.
            
            checkIllegal = re.search("[\$\(\)\'\"&]+", eachId) # Illegal characters are following:
                                                               # '$' (dollar sign); round brackets; single and double quotation marks; and '&' (ampersand)
            
            if checkIllegal:
                # If illegal characters found, we report this to the users.
                
                logHandle.write("\t\tERROR: Sequence ID %s has one of the illegal characters\n" % (eachId))
                logHandle.write("\t\tSequence ID format test FAILED ...\n")
                logHandle.write("\t\tPlease make sure that sequence ID does not contain any of the following characters: $ ( ) \' \"\n")
                logHandle.write("\t\t\t$ or & (dollar sign or ampersand)\n")
                logHandle.write("\t\t\t( ) (round brackets)\n")
                logHandle.write("\t\t\t \'\' or \"\" (single or double-quotation marks)\n")
                
                print "An ERROR was reported during the input file formatting procedure."
                print "Please review 'log.txt' for details."
            
                sys.exit()
            else:
                # If no illegal characters are found, we then check for characters that could be replaced by '_' (underscore).
                #seqTypeChosen = options.seq_type
                
                eachSeq = fasSeqDict[eachId]
                
                checkReplace = re.search("[\s\.,;:]+" ,eachId)
                    
                # We also want to check whether the users are inserting correct type of sequences (e.g. if nucleotide selected, then nucleotide sequences should be inserted).
                
                nucMatch = re.match("^[acgtnxACGTNX]+$", eachSeq)
                protMatch = re.match("[rndeqghilkmfpswyvRNDEQGHILKMFPSWYV]*", eachSeq) # amino acid regex...this will change soon to something more efficient
                
                if checkReplace:
                    # If replaceable characters found, than let's first report this, and then replace theses characters with '_'.
                    logHandle.write("\t\tWARNING: Sequence ID %s contains one of following: . (period), , (comma), ; (semi-colon), : (colon), or space.\n" % (eachId))
                    logHandle.write("\t\tThese characters will be replaced by \"_\"\n")
                    
                    print "WARNING: Sequence ID %s contains one of following: . (period), , (comma), ; (semi-colon), : (colon), or space." % (eachId)
                    print "These characters will be replaced by \"_\""
                    print ""
                    
                    ######## TO DO: Replace characters and put the formatted ID string into the dictionary
                    
                    # Also, if nucleotide is selected as molecular type...
                    if seqTypeChosen is "n":
                        
                        # and the user did indeed insert the nucleotide sequences in FASTA file,
                        if nucMatch:
                            
                            logHandle.write("\t\t\t* Sequence content format test PASSED ...\n" % (eachTempId))
                            
                            # then check and see if identical sequences are already present in our list of sequences.
                            if len(finalSeqDict) == 0:
                                
                                # If the final list of sequences empty, then insert our first sequence into the list.
                                finalSeqDict[eachId] = eachSeq
                                
                            else:
                                
                                # If the final list of sequence is not empty, then we need to compare our current sequence with those already in the final list.
                                for eachTempId in sorted(finalSeqDict.keys()):
                                    if finalSeqDict[eachTempId] == eachSeq:
                                        
                                        # If there already is sequence that is identical to the one we are processing now, turn on the switch that indicates
                                        # the presence of identital sequence, and report this as warning.
                                        identicalSeqs += 1
                                        
                                        logHandle.write("\t\t\tWARNING: Sequences %s and %s are identical\n" % (eachTempId, eachId))
                                        logHandle.write("\t\t\tAs such, only %s will be used in the analysis\n" % (eachTempId))
                                        
                                        print "WARNING: Sequences %s and %s are identical" % (eachTempId, eachId)
                                        print "As such, only %s will be used in the analysis" % (eachTempId)
                                        print ""
                            
                            # If the switch is not on, then we check the sequence for presence of illegal characters.
                            # If no illegal characters found in sequence, can insert the sequence into the final list.            
                            if identicalSeqs == 0:
                                checkIllegal = re.search("[\$\(\)\'\"&]+", eachSeq)
                                
                                if checkIllegal:
                                    logHandle.write("\t\tERROR: Sequence ID %s has one of the illegal characters\n" % (eachId))
                                    logHandle.write("\t\tSequence ID format test FAILED ...\n")
                                    logHandle.write("\t\tPlease make sure that sequence ID does not contain any of the following characters: $ ( ) \' \"\n")
                                    logHandle.write("\t\t\t$ or & (dollar sign or ampersand)\n")
                                    logHandle.write("\t\t\t( ) (round brackets)\n")
                                    logHandle.write("\t\t\t \'\' or \"\" (single or double-quotation marks)\n")
                
                                    print "An ERROR was reported during the input file formatting procedure."
                                    print "Please review 'log.txt' for details."
                                    
                                    sys.exit()
                                else:
                                    finalSeqDict[eachId] = eachSeq
                        else:
                            # However, if the user selected nucleotide as molecular types, and have inserted FASTA sequences consisting of
                            # something other than nucleotide, then we report this as error and exit the program.
                            
                            logHandle.write("\t\t\tERROR: You chose to analyse nucleotide, and there's something other than ACGT in sequence, %s\n" % (eachId))
                            logHandle.write("\t\t\tMake sure the input FASTA file only has ACGT's for each sequence.\n")
                            logHandle.write("\t\t\tSequence content format test FAILED ...\n")
                            
                            print "An ERROR was reported during the input file formatting procedure."
                            print "Please review 'log.txt' for details."
                            
                            sys.exit()
                    
                    # If protein is selected as molecular type...
                    elif seqTypeChosen is "p":
                        
                        # and if the user did indeed insert the protein sequences,
                        if protMatch:
                            
                            logHandle.write("\t\t* Sequence content format test PASSED ...\n")
                            
                            # then check and see if identical sequences are already present in our list of sequences.
                            if len(finalSeqDict) == 0:
                                
                                # If the final list of sequences empty, then insert our first sequence into the list.
                                finalSeqDict[eachId] = eachSeq
                                
                            else:
                                
                                # If the final list of sequence is not empty, then we need to compare our current sequence with those already in the final list.
                                for eachTempId in sorted(finalSeqDict.keys()):
                                    if finalSeqDict[eachTempId] == eachSeq:
                                        
                                        # If there already is sequence that is identical to the one we are processing now, turn on the switch that indicates
                                        # the presence of identital sequence, and report this as warning.
                                        identicalSeqs += 1
                                        
                                        logHandle.write("\t\t\tWARNING: Sequences %s and %s are identical\n" % (eachTempId, eachId))
                                        logHandle.write("\t\t\tAs such, only %s will be used in the analysis\n"% (eachTempId))
                                        
                                        print "WARNING: Sequences %s and %s are identical" % (eachTempId, eachId)
                                        print "As such, only %s will be used in the analysis" % (eachTempId)
                                        print ""
                            
                            # If the switch is not on, then we check the sequence for presence of illegal characters.
                            # If no illegal characters found in sequence, can insert the sequence into the final list.            
                            if identicalSeqs == 0:
                                checkIllegal = re.search("[\$\(\)\'\"&]+", eachSeq)
                                
                                if checkIllegal:
                                    logHandle.write("\t\tERROR: Sequence ID %s has one of the illegal characters\n" % (eachId))
                                    logHandle.write("\t\tSequence ID format test FAILED ...\n")
                                    logHandle.write("\t\tPlease make sure that sequence ID does not contain any of the following characters: $ ( ) \' \"\n")
                                    logHandle.write("\t\t\t$ or & (dollar sign or ampersand)\n")
                                    logHandle.write("\t\t\t( ) (round brackets)\n")
                                    logHandle.write("\t\t\t \'\' or \"\" (single or double-quotation marks)\n")
                
                                    print "An ERROR was reported during the input file formatting procedure."
                                    print "Please review 'log.txt' for details."
                                    
                                    sys.exit()
                                else:
                                    finalSeqDict[eachId] = eachSeq
                                
                        elif nucMatch:
                            # However, if the user selected protein as molecular types, and have inserted FASTA sequences consisting of
                            # nucleotides, then we report this as error and exit the program.
                            
                            logHandle.write("\t\t\tERROR: You chose to analyse protein sequences, and there are nucleotides in sequence, %s\n" % (eachId))
                            logHandle.write("\t\t\tMake sure the input FASTA file only has amino acids for each sequence.\n")
                            logHandle.write("\t\t\tSequence content format test FAILED ...\n")
                            
                            print "An ERROR was reported during the input file formatting procedure."
                            print "Please review 'log.txt' for details."
                            
                            sys.exit()
                        
                        else:
                            # However, if the user selected protein as molecular types, and have inserted FASTA sequences consisting of
                            # something other than nucleotides and amino acids (e.g. illegal characters), then we report this as error and exit the program.
                            
                            logHandle.write("\t\t\tERROR: You chose to analyse protein sequences, and there's something other than amino acid letters in sequence, %s\n" % (eachId))
                            logHandle.write("\t\t\tMake sure the input FASTA file only has amino acids for each sequence.\n")
                            logHandle.write("\t\t\tSequence content format test FAILED ...\n")
                            
                            print "An ERROR was reported during the input file formatting procedure."
                            print "Please review 'log.txt' for details."
                            
                            sys.exit()
                        
                else:
                    # If there are no illega/replaceable characters in the ID string ,then we check for duplicates and molecular type consistency as done previously
                    # in lines 149 to 235.
                    logHandle.write("\t\tSequence ID format test PASSED ...\n")
                    
                    # Also, if nucleotide is selected as molecular type...
                    if seqTypeChosen is "n":
                        
                        # and the user did indeed insert the nucleotide sequences in FASTA file,
                        if nucMatch:
                            logHandle.write("\t\t\t* Sequence content format test PASSED ...\n")
                            
                            # then check and see if identical sequences are already present in our list of sequences.
                            if len(finalSeqDict) == 0:
                                
                                # If the final list of sequences empty, then insert our first sequence into the list.
                                finalSeqDict[eachId] = eachSeq
                                
                            else:
                                
                                # If the final list of sequence is not empty, then we need to compare our current sequence with those already in the final list.
                                for eachTempId in sorted(finalSeqDict.keys()):
                                    if finalSeqDict[eachTempId] == eachSeq:
                                        
                                        # If there already is sequence that is identical to the one we are processing now, turn on the switch that indicates
                                        # the presence of identital sequence, and report this as warning.
                                        identicalSeqs += 1
                                        
                                        logHandle.write("\t\t\tWARNING: Sequences %s and %s are identical\n" % (eachTempId, eachId))
                                        logHandle.write("\t\t\tAs such, only %s will be used in the analysis\n" % (eachTempId))
                                        
                                        print "WARNING: Sequences %s and %s are identical" % (eachTempId, eachId)
                                        print "As such, only %s will be used in the analysis" % (eachTempId)
                                        print ""
                            
                            # If the switch is not on, then we check the sequence for presence of illegal characters.
                            # If no illegal characters found in sequence, can insert the sequence into the final list.            
                            if identicalSeqs == 0:
                                checkIllegal = re.search("[\$\(\)\'\"&]+", eachSeq)
                                
                                if checkIllegal:
                                    logHandle.write("\t\tERROR: Sequence ID %s has one of the illegal characters\n" % (eachId))
                                    logHandle.write("\t\tSequence ID format test FAILED ...\n")
                                    logHandle.write("\t\tPlease make sure that sequence ID does not contain any of the following characters: $ ( ) \' \"\n")
                                    logHandle.write("\t\t\t$ or & (dollar sign or ampersand)\n")
                                    logHandle.write("\t\t\t( ) (round brackets)\n")
                                    logHandle.write("\t\t\t \'\' or \"\" (single or double-quotation marks)\n")
                
                                    print "An ERROR was reported during the input file formatting procedure."
                                    print "Please review 'log.txt' for details."
                                    
                                    sys.exit()
                                else:
                                    finalSeqDict[eachId] = eachSeq
                        else:
                            logHandle.write("\t\t\tERROR: You chose to analyse nucleotide, and there's something other than ACGT in sequence, %s\n" % (eachId))
                            logHandle.write("\t\t\tMake sure the input FASTA file only has ACGT's for each sequence.\n")
                            logHandle.write("\t\t\tSequence content format test FAILED ...\n")
                            
                            print "An ERROR was reported during the input file formatting procedure."
                            print "Please review 'log.txt' for details."
                            
                            sys.exit()
                            
                    # If protein is selected as molecular type...
                    elif seqTypeChosen is "p":
                        
                        # and if the user did indeed insert the protein sequences,
                        if protMatch:
                            
                            logHandle.write("\t\t* Sequence content format test PASSED ...\n")
                            
                            # then check and see if identical sequences are already present in our list of sequences.
                            if len(finalSeqDict) == 0:
                                
                                # If the final list of sequences empty, then insert our first sequence into the list.
                                finalSeqDict[eachId] = eachSeq
                                
                            else:
                                
                                # If the final list of sequence is not empty, then we need to compare our current sequence with those already in the final list.
                                for eachTempId in sorted(finalSeqDict.keys()):
                                    if finalSeqDict[eachTempId] == eachSeq:
                                        
                                        # If there already is sequence that is identical to the one we are processing now, turn on the switch that indicates
                                        # the presence of identital sequence, and report this as warning.
                                        identicalSeqs += 1
                                        
                                        logHandle.write("\t\t\tWARNING: Sequences %s and %s are identical\n" % (eachTempId, eachId))
                                        logHandle.write("\t\t\tAs such, only %s will be used in the analysis\n" % (eachTempId))
                                        
                                        print "WARNING: Sequences %s and %s are identical" % (eachTempId, eachId)
                                        print "As such, only %s will be used in the analysis" % (eachTempId)
                                        print ""
                            
                            # If the switch is not on, then we check the sequence for presence of illegal characters.
                            # If no illegal characters found in sequence, can insert the sequence into the final list.            
                            if identicalSeqs == 0:
                                checkIllegal = re.search("[\$\(\)\'\"&]+", eachSeq)
                                
                                if checkIllegal:
                                    logHandle.write("\t\tERROR: Sequence ID %s has one of the illegal characters\n" % (eachId))
                                    logHandle.write("\t\tSequence ID format test FAILED ...\n")
                                    logHandle.write("\t\tPlease make sure that sequence ID does not contain any of the following characters: $ ( ) \' \"\n")
                                    logHandle.write("\t\t\t$ or & (dollar sign or ampersand)\n")
                                    logHandle.write("\t\t\t( ) (round brackets)\n")
                                    logHandle.write("\t\t\t \'\' or \"\" (single or double-quotation marks)\n")
                
                                    print "An ERROR was reported during the input file formatting procedure."
                                    print "Please review 'log.txt' for details."
                                    
                                    sys.exit()
                                else:
                                    finalSeqDict[eachId] = eachSeq
                                
                        elif nucMatch:
                            # However, if the user selected protein as molecular types, and have inserted FASTA sequences consisting of
                            # nucleotides, then we report this as error and exit the program.
                            
                            logHandle.write("\t\t\tERROR: You chose to analyse protein sequences, and there are nucleotides in sequence, %s\n" % (eachId))
                            logHandle.write("\t\t\tMake sure the input FASTA file only has amino acids for each sequence.\n")
                            logHandle.write("\t\t\tSequence content format test FAILED ...\n")
                            
                            print "An ERROR was reported during the input file formatting procedure."
                            print "Please review 'log.txt' for details."
                                    
                            sys.exit()
                        
                        else:
                            # However, if the user selected protein as molecular types, and have inserted FASTA sequences consisting of
                            # something other than nucleotides and amino acids (e.g. illegal characters), then we report this as error and exit the program.
                            
                            logHandle.write("\t\t\tERROR: You chose to analyse protein sequences, and there's something other than amino acid letters in sequence, %s\n" % (eachId))
                            logHandle.write("\t\t\tMake sure the input FASTA file only has amino acids for each sequence.\n")
                            logHandle.write("\t\t\tSequence content format test FAILED ...\n")
                            
                            print "An ERROR was reported during the input file formatting procedure."
                            print "Please review 'log.txt' for details."
                            
                            sys.exit()
    #################### End of Part A3 ####################
    logHandle.close()                
    fastaCheckHandle.close()

    #################### End of Part A ####################
   
    #################### Part B: Split the FASTA files into multiple FASTA files ####################
    #print len(finalSeqDict)
    
    print "************************************        Executing BLAST        ************************************"
    
    outBLAST = outDir + "/blast_out"
     
    os.system("mkdir %s" % outBLAST)
    
    if len(finalSeqDict) < 50:
        fastaSplit = outBLAST + "/" + "_0" + ".fasta"
        
        fasSplitHandle = open(fastaSplit, 'w')
        
        for finalSeqId in sorted(finalSeqDict.keys()):
            fasSplitHandle.write(finalSeqId + "\n")
            fasSplitHandle.write(finalSeqDict[finalSeqId] + "\n")
        
        fasSplitHandle.close()
    else:
        #splitCounter = 0
        for i, finalSeqId in enumerate(finalSeqDict.keys()):
            #print i, finalSeqId
            
            if i % 50 == 0: # Change this number from 50 to 500
                splitCounter = i/50
                fastaSplit = outBLAST + "/" + "_" + str(splitCounter) + ".fasta"
                fasSplitHandle = open(fastaSplit, 'w')
                   
            if ((i + 1) % 50 != 0):
                #fasSplitHandle.write(str(i) + " " + finalSeqId + "\n")
                fasSplitHandle.write(finalSeqId + "\n")
                fasSplitHandle.write(finalSeqDict[finalSeqId] + "\n")
            else:
                #fasSplitHandle.write(str(i) + " " + finalSeqId + "\n")
                fasSplitHandle.write(finalSeqId + "\n")
                fasSplitHandle.write(finalSeqDict[finalSeqId] + "\n")
                fasSplitHandle.close()
            
            if i == len(fastaLines) - 1:
                #fasSplitHandle.write(str(i) + " " + finalSeqId + "\n")
                fasSplitHandle.write(finalSeqId + "\n")
                fasSplitHandle.write(finalSeqDict[finalSeqId] + "\n")
                fasSplitHandle.close()
    #################### End Part B  ####################
    
    #################### Part C: Execute BLAST  ####################
    
    treeType = options.phylo_tree
    
    alignData = ''
    hmmData   = ''
    treeData  = "data/tree_data"
        
    if treeType is 'm':
        alignData = "data/alignment_data"
        hmmData   = "data/hmm_data"
    elif treeType is 'g':
        alignData = "data/geba_alignment_data"
        hmmData   = "data/geba_hmm_data"
    elif treeType is 'f':
        alignData = "data/fungi_alignment_data"
        hmmData   = "data/fungi_hmm_data"
    
    mltreeUtility = mltree_utility(alignData, hmmData, treeData)
    
    refSelect = options.ref_search
    
    #mltreeUtility.searchRange(seqTypeChosen, refSelect)
    mltreeUtility.runBLAST(seqTypeChosen, refSelect, outDir, alignData)
    
    print "************************************     Reading/Processing BLAST Output    ************************************"
    
    inBitScore = options.blast_bit
    
    if inBitScore == None:
        print "BLAST bit score was not defined in the input parameter"
        print "As a result, BLAST bit score is set to 60"
        
        inBitScore = 60
        blastResult = mltreeUtility.readBLAST(inBitScore, refSelect, outBLAST)
    else:
        print "BLAST bit score was set to %s by user" % (inBitScore)
        
        blastResult = mltreeUtility.readBLAST(inBitScore, refSelect, outBLAST)
    
    outGeneWise = outDir + "/genewise_out"
    
    os.system("mkdir %s" % outGeneWise)
    
    print "************************************     Preparing GeneWise Input Files    ************************************"
    
    coordShortSeq = mltreeUtility.shortSeqForGenwise(inputFASTA, blastResult, fasSeqDict, outGeneWise)
        
    print len(coordShortSeq)
        
    contigCrdDict = coordShortSeq[0]
    shortSeqDict = coordShortSeq[1]
    
    print "************************************     Executing GeneWise    ************************************"
    # We need to execute GeneWise for each filtered BLAST query and collect the information regarding the output file for each run
    geneWiseOutFiles = mltreeUtility.runGeneWise(blastResult, outGeneWise, shortSeqDict, hmmData)
    
    # The output files of GeneWise are processed here
    mltreeUtility.parseGeneWise(geneWiseOutFiles, contigCrdDict)
    
    
    
                
        
        
        #print fastaSplit
        
    #def splitFASTA():
    #    print "hello"
        