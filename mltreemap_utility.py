#!/usr/bin/python

##########       mltreemap_utility.py by Young C. Song        ##########
##########              Written on August 3, 2013             ##########
##########       functional utility class for mltreemap.py    ##########
##########       required to run dependcy tools for mltreemap ##########

from optparse import OptionParser
from time import gmtime, strftime
from os import listdir

from decimal import *
import string
import os
import shutil
import re
import sys
import random
import multiprocessing

class mltree_workflow:

    def __init__(self, alignData, hmmData, treeData):
        self.alignData = alignData
        self.hmmData   = hmmData
        self.treeData  = treeData
    
    
    ########## Function dataRange: Read cog_list.txt and return set of all available reference data ##########
    def dataSet (self, refSelect):
        refDict = {}
        
        ref_list = ''
        
        if refSelect is "a":
            ref_list = self.treeData + "/cog_list.txt"
        elif refSelect is "s":
            ref_list = self.treeData + "/reference_selection.txt"
        
        listHandle = open(ref_list,'rb')
        
        for refLine in listHandle.readlines():
            refLine = string.strip(refLine)
            title = re.match("^#(\S+\s+)*\S+", refLine)
            
            if not title:
                phyloCog = re.match("^COG(\d){4}$", refLine)                  # phylogenetic marker (COGs): line only has COG IDs
                nonCog   = re.match("^(\S+\t){2}((\S+\s+)*\S+)", refLine)     # non-cog IDs: contains alignment file prefix, reference code and description separated by tab
                
                if phyloCog:
                    refDict[refLine] = "p"
                elif nonCog:
                    nonCogFields = refLine.split("\t")
                    refCodeFile = nonCogFields[0] # e.g. mcrAref, 16srRNA
                    refCode     = nonCogFields[1] # e.g. mcrA, a
                    
                    refDict[refCodeFile] = refCode
        
        listHandle.close()
        
        return refDict
    ########## End of Function dataRange ##########
    
    ########## Function searchRange: Determines which set of reference data that will be searched against in the current analysis ##########
    def refDB (self, seqTypeChosen, refSelect, alignData):
        availData = self.dataSet(refSelect)
        
        refSearch = ""
        blastDB = {}
            
        # If using nucleotide sequences as input, we search against all reference data sets.
        if seqTypeChosen is "n":
            
            rnaSearch = alignData + "/16srRNA.fa " + alignData + "/18srRNA.fa"
            blastDB["rRNA"] = rnaSearch
            
            protSearch = ""
            
            for i, refData in enumerate(sorted(availData.keys())):
                if not (refData == '16srRNA') and not (refData == '18srRNA'):
                    if i == len(availData) - 1:
                        protSearch += alignData + "/" +  refData + ".fa"
                        blastDB["prot"] = protSearch
                    else:
                        protSearch += alignData + "/" + refData + ".fa "
            
                        
            
        # If using protein sequences as input, we don't want to be searching against 16S and 18S rRNA data
        elif seqTypeChosen is "p":
            
            protSearch = ""
            
            for i, refData in enumerate(sorted(availData.keys())):
                if not (refData == '16srRNA') and not (refData == '18srRNA'):
                    if i == len(availData) - 1:
                        protSearch += alignData + "/" +  refData + ".fa"
                        blastDB["prot"] = protSearch
                    else:
                        protSearch  += alignData + "/" +  refData + ".fa "
            
            
        
        return blastDB
    ########## End of function searchRange ##########
    
    ########## Function calcOverlap: calculates overlap between two reads being compared ##########
    def calcOverlap(self, infoDict):   
        
        baseStart  = infoDict["baseStart"]
        baseEnd    = infoDict["baseEnd"]
        checkStart = infoDict["checkStart"]
        checkEnd   = infoDict["checkEnd"]
        
        overlap = 0;
    
        #now check, if they overlap:
        if ((baseStart <= checkStart) and (checkStart <= baseEnd) and (baseEnd <= checkEnd)):
            # Base     ----------       
            # Check        ----------   
            overlap = int(baseEnd) - int(checkStart)
            #print "1. overlap is %s" % (str(overlap))
        elif ((baseStart <= checkStart) and (checkEnd <= baseEnd)):
            # Base     ----------       
            # Check        ---          
            overlap = int(checkEnd) - int(checkStart)
            #print "2. overlap is %s" % (str(overlap))
        elif ((checkStart <= baseStart) and (baseStart <= checkEnd) and (checkEnd <= baseEnd)):
            # Base           ------     
            # Check   ----------        
            overlap = int(checkEnd) - int(baseStart)
            #print "3. overlap is %s" % (str(overlap))
        elif ((checkStart <= baseStart) and (baseEnd <= checkEnd)):
            # Base           ------     
            # Check   ----------------   
            overlap = int(baseEnd) - int(baseStart)
            #print "4. overlap is %s" % (str(overlap))
            
        return overlap
    ########## End of function calcOverlap ##########
    
    ########## Function runBLAST: Executes blast command, using the reference data information collected using previous functions ##########
    def runBLAST (self, seqTypeChosen, refSelect, outDir, alignData):
        #print "Hi this is detectBLASTDB"
        
        blastType = ''
        blastDB = self.refDB(seqTypeChosen, refSelect, alignData)
        
        protDB = ''
        rnaDB = ''
        
        #numProcessors = multiprocessing.cpu_count()
        #print numProcessors
        
        if seqTypeChosen is "n":
            blastType = "blastx"
            blastNuc  = "blastn"
        elif seqTypeChosen is "p":
            blastType = "blastp"
        
        blastFolder = "%s/blast_out" % (outDir)
        
        for eachFile in listdir(blastFolder):
            fastaFile = re.match("^(_\d+)\.fasta$", eachFile)
            if fastaFile:
                blastCommand = "sub_binaries/blastall"
                
                blastIn = blastFolder + "/" + eachFile
                
                
                if seqTypeChosen is "n":
                    blastProt = "blastx"
                    blastNuc = "blastn"
                    
                    blastProtOut = blastFolder + "/" + fastaFile.group(1) + "_prot.blastout"
                    blastNucOut = blastFolder + "/" + fastaFile.group(1) + "_rRNA.blastout"
                    
                    #blastProt += " -p %s -i %s -M BLOSUM62 -d \"%s\" -e 0.01 -v 20000 -b 20000 -z 1000000 -m 8 > %s" % (blastType, blastIn, blastDB, blastOut)
                    blastProt = blastCommand + " -p %s -i %s -M BLOSUM62 -d \"%s\" -e 0.01 -v 10 -b 10 -z 1000000 -m 8 > %s" % (blastProt, blastIn, blastDB["prot"], blastProtOut)  # This is for the test purpose
                    print blastProt
                    print ""
                    os.system(blastProt)
                    
                    blastNuc = blastCommand + " -p %s -i %s -M BLOSUM62 -d \"%s\" -e 0.01 -v 10 -b 10 -z 1000000 -m 8 > %s" % (blastNuc, blastIn, blastDB["rRNA"], blastNucOut)
                    print blastNuc
                    print ""
                    os.system(blastNuc)
                
                elif seqTypeChosen is "p":
                    blastProt = "blastp"
                    
                    blastProtOut = blastFolder + "/" + fastaFile.group(1) + "_prot.blastout"
                    
                    blastCommand += " -p %s -i %s -M BLOSUM62 -d \"%s\" -e 0.01 -v 10 -b 10 -z 1000000 -m 8 > %s" % (blastType, blastIn, blastDB["prot"], blastProtOut)  # This is for the test purpose
                    print blastCommand
                    print ""
                    
                    os.system(blastProt)
    ########## End of function runBLAST ##########
                    
    ########## Function readBLAST: Read and parse BLAST output ##########
    def readBLAST(self, inBitScore, refSelect, outBLAST):
        #print "**",inBitScore
        #blastFolder = "%s/blast_out" % (outDir)
        
        blastHitDict = {}
        blastFilterDict = {}
        
        ##### We are first going to read the blast hits, will collect the following information:
        ##### Query ID (combination of query ID and integer)
        ##### Matching reference sequence
        ##### Query start and end
        ##### Sequence direction
        ##### ... and a integer (0 or 1) indicating the validity of BLAST hit.
        
        for eachFile in listdir(outBLAST):
            
            blastOutMatch = re.match("^(_\d+)_\S+\.blastout$", eachFile)
            
            if blastOutMatch:
                blastHitId = 0
                
                blastOut = outBLAST + "/" + eachFile
                #print "Opening %s" % blastOut
                
                blastHandle = open(blastOut, "rb")
        
                blastRead = blastHandle.readlines()
        
                for blastLine in blastRead:
                    blastLine = blastLine.strip()
            
                    (query, subject, percID, alignLength, mismatches, gapOpen, qStart, qEnd, rStart, rEnd, eVal, bitScore) = blastLine.split("\t")
            
                    if float(bitScore) <= float(inBitScore):
                        continue
                    
                    direction = "forward"
                    
                    if int(rStart) > int(rEnd):
                        #print rStart,"\trEnd: ",rEnd
                        tempVar = rStart
                        rStart = rEnd
                        rEnd = tempVar
                        direction = "reverse"
                        
                    if int(qStart) > int(qEnd):
                        tempVar = qStart
                        qStart = qEnd
                        qEnd = tempVar
                        
                        if (direction is "reverse"):
                            #print query, "\t", subject, "\trefStart: ", rStart, "\trefEnd: ", rEnd, "\tqStart: ", qStart, "\tqEnd: ",qEnd
                            print "ERROR: Parsing error with the BLAST results. Please notify the authors\n"
                            sys.exit()
                            
                        direction = "reverse"
                    
                    refSeqs = ""
                    
                    refMatch = re.match("\d+_(\S{7})$", subject)
                    
                    if refMatch:
                        refSeqs = refMatch.group(1)
                    else:
                        print "ERROR: Cannot find the appropriate reference sequence for sequence %s" % (query)
                        sys.exit()
                    
                    #validity = 1
                    
                    #queryBlastHit = query + ":" + str(blastHitId)
                    queryBlastHit = query
                    
                    #blastList[blastHitId] = (str(bitScore), refSeqs, str(qStart), str(qEnd), direction, str(validity))
                    #print queryBlastHit + "\t" + str(blastHitId)
                    
                    if not queryBlastHit in blastHitDict:
                        blastHitId = 0
                        blastHitDict[queryBlastHit] = {}
                    
                    if not blastHitId in blastHitDict[queryBlastHit]:
                        blastHitDict[queryBlastHit][blastHitId] = {}
                    
                    blastHitDict[queryBlastHit][blastHitId]["bitScore"]  = str(bitScore)
                    blastHitDict[queryBlastHit][blastHitId]["refSeqs"]   = refSeqs
                    blastHitDict[queryBlastHit][blastHitId]["qStart"]    = str(qStart)
                    blastHitDict[queryBlastHit][blastHitId]["qEnd"]      = str(qEnd)
                    blastHitDict[queryBlastHit][blastHitId]["direction"] = direction
                    blastHitDict[queryBlastHit][blastHitId]["validity"]  = 1
                    
                    blastHitId += 1
        
                blastHandle.close()
                
                ##### Let's purify the BLAST hits
        
                refList = self.dataSet(refSelect)
                
                for baseBlastHit in sorted(blastHitDict.keys()):
                    #print baseBlastHit
                    filterBlastID = 0
            
                    for blastHitId in sorted(blastHitDict[baseBlastHit].keys()):
                        #print blastHitId
                
                        baseBit   = blastHitDict[baseBlastHit][blastHitId]["bitScore"]
                        baseRef   = blastHitDict[baseBlastHit][blastHitId]["refSeqs"]
                        baseStart = blastHitDict[baseBlastHit][blastHitId]["qStart"]
                        baseEnd   = blastHitDict[baseBlastHit][blastHitId]["qEnd"]
                        baseDir   = blastHitDict[baseBlastHit][blastHitId]["direction"]
                        baseValid = blastHitDict[baseBlastHit][blastHitId]["validity"]
                
                        if float(baseBit) <= float(inBitScore):
                            continue
         
                        if not refList.has_key(baseRef):
                            print "%s does not exist" % (baseRef)
                            blastHitDict[baseBlastHit][blastHitId]["validity"] = 0
                
                        for checkHitId in sorted(blastHitDict[baseBlastHit].keys()):
                    
                            if blastHitId == checkHitId:
                                continue
                    
                            checkBit   = blastHitDict[baseBlastHit][checkHitId]["bitScore"]
                            checkRef   = blastHitDict[baseBlastHit][checkHitId]["refSeqs"]
                            checkStart = blastHitDict[baseBlastHit][checkHitId]["qStart"]
                            checkEnd   = blastHitDict[baseBlastHit][checkHitId]["qEnd"]
                            checkDir   = blastHitDict[baseBlastHit][checkHitId]["direction"]
                            checkValid = blastHitDict[baseBlastHit][checkHitId]["validity"]
                    
                            baseLength = int(baseEnd) - int(baseStart)
                            #print "** baseLength is %s" % str(baseLength)
                    
                            checkLength = int(checkEnd) - int(checkStart)
                            #print "&& checkLength is %s" % str(checkLength)
        
                            infoDict = {}
              
                            infoDict["baseStart"] = baseStart
                            infoDict["baseEnd"] = baseEnd
                            infoDict["checkStart"] = checkStart
                            infoDict["checkEnd"] = checkEnd
                
                            overlap = self.calcOverlap(infoDict)
                            #print "...and overlap is %s" % (overlap)
                
                            if (overlap):
                                if (((overlap / baseLength) > 0.5) and (float(baseBit) < float(checkBit))):
                                    #print "Ia. ratio is %s" % (str(overlap/baseLength))
                                    #print "Ib. checkBit is greater than baseBit"
                            
                                    blastHitDict[baseBlastHit][blastHitId]["validity"] = 0
                        
                                elif (((overlap / checkLength) > 0.5) and (float(checkBit) < float(baseBit))):
                                    #print "IIa. ratio is %s" % (str(overlap/checkLength))
                                    #print "IIb. baseBit is greater than check"
                            
                                    blastHitDict[baseBlastHit][checkHitId]["validity"] = 0
                        
                                elif ((baseStart == checkStart) and (baseEnd == checkEnd)): #if both are the same keep only the one with the smaller identifier.
                                    if (int(checkHitId) > int(blastHitId)):
                                        #print "IIIa. baseStart is checkStart and baseEnd is checkEnd"
                                        #print "IIIb. checkID is greater than queryID"
                                
                                        blastHitDict[baseBlastHit][checkHitId]["validity"] = 0
        
                                    else:
                                        #print "IV"
                                        blastHitDict[baseBlastHit][blastHitId]["validity"] = 0
                 
                        #print "*** after ----> %s\t%s\t%s\t%s\t%s\t%s\t%s" % (baseBlastHit, str(baseBit), baseRef, str(baseStart), str(baseEnd), baseDir, str(baseValid))
                
                        if(blastHitDict[baseBlastHit][blastHitId]["validity"] == 1):
                            
                            if not baseBlastHit in blastFilterDict:
                                blastFilterDict[baseBlastHit] = {}
                                
                            if not filterBlastID in blastFilterDict[baseBlastHit]:
                                blastFilterDict[baseBlastHit][filterBlastID] = {}
                                
                            blastFilterDict[baseBlastHit][filterBlastID]["bitScore"]  = baseBit
                            blastFilterDict[baseBlastHit][filterBlastID]["refSeqs"]   = baseRef
                            blastFilterDict[baseBlastHit][filterBlastID]["qStart"]    = baseStart
                            blastFilterDict[baseBlastHit][filterBlastID]["qEnd"]      = baseEnd
                            blastFilterDict[baseBlastHit][filterBlastID]["direction"] = baseDir
                            blastFilterDict[baseBlastHit][filterBlastID]["already_placed"]  = 0
                            
                            filterBlastID += 1
        
        # Writing the filter results to output files. Need to sort the filter dictionary by query ID then by the bitscore.
        for purifiedID in sorted(blastFilterDict.keys()):
            blastPurifiedFile   = outBLAST + "/" + purifiedID + "_blast_result_purified.txt"
            
            blastPurifiedHandle = open(blastPurifiedFile, "w")
            
            blastSortDict = {}
            
            for filterBlastID in sorted(blastFilterDict[purifiedID].keys()):
                
                filterBlastScore = blastFilterDict[purifiedID][filterBlastID]["bitScore"]
                
                if not filterBlastScore in blastSortDict:
                    blastSortDict[str(filterBlastScore)] = {}
                    
                if not filterBlastID in blastSortDict[str(filterBlastScore)]:
                    blastSortDict[str(filterBlastScore)][filterBlastID] = 1
            
            for sortBitScore in sorted(blastSortDict.keys()):
                for sortID in sorted(blastSortDict[sortBitScore].keys()):
                    
                    refSeq    = blastFilterDict[purifiedID][sortID]["refSeqs"]
                    start     = blastFilterDict[purifiedID][sortID]["qStart"]
                    end       = blastFilterDict[purifiedID][sortID]["qEnd"]
                    direction = blastFilterDict[purifiedID][sortID]["direction"]
                    
                    blastPurifiedHandle.write(purifiedID + "\t" + refSeq + "\t" + str(start) + "\t" + str(end) + "\t" + direction + "\t" + str(sortBitScore) + "\n")
                    
            blastPurifiedHandle.close()
        
        return blastFilterDict
    ########## End of function readBLAST ##########
        
    ########## Function shortSeqForGenwise: Read and parse BLAST output ##########
    def shortSeqForGenwise(self, inputFASTA, filterBLAST, fasSeqDict, outGeneWise):
        
        flankingLength = 0
        
        intermContigCoordDict = {}
        contigCoordDict = {}
        shortenedSeqDict = {}
    
        for eachQuery in sorted(filterBLAST.keys()):
            for eachIndex in sorted(filterBLAST[eachQuery].keys()):
                start     = filterBLAST[eachQuery][eachIndex]["qStart"] 
                end       = filterBLAST[eachQuery][eachIndex]["qEnd"]
                placed    = filterBLAST[eachQuery][eachIndex]["already_placed"]
                
                if placed == 1:
                    continue
                
                filterBLAST[eachQuery][eachIndex]["already_placed"] = 1
                
                baseStart = int(start) - flankingLength
                baseEnd   =int(end) + flankingLength
                
                numBLASTHitsPerQuery = len(sorted(filterBLAST[eachQuery]))
                
                for checkIndex in range(numBLASTHitsPerQuery):
                    
                    if int(filterBLAST[eachQuery][checkIndex]["already_placed"]) == 1:
                        continue
                    
                    checkStart = int(filterBLAST[eachQuery][checkIndex]["qStart"]) - flankingLength
                    checkEnd  = int(filterBLAST[eachQuery][checkIndex]["qEnd"]) - flankingLength
                    
                    #now check, if they overlap. Note: the small subroutine "calculate_overlap" can't be used here because here we merge_stuff.
                    
                    if((baseStart <= checkStart) and (checkStart <= baseEnd) and (baseEnd <= checkEnd)):
                        # Base     ----------       
                        # Check        ----------   
                        baseEnd = checkEnd
                        filterBLAST[eachQuery][checkIndex]["already_placed"] = 1
                        checkIndex -= 1
                        continue
                    elif((baseStart <= checkStart) and (checkEnd <= baseEnd)):
                        # Base     ----------       
                        # Check        ---
                        filterBLAST[eachQuery][checkIndex]["already_placed"] = 1
                        checkIndex -= 1
                        continue
                    elif((checkStart <= baseStart) and (baseStart <= checkEnd) and (checkEnd <= baseEnd)):
                        # Base           ------     
                        # Check   ----------
                        baseStart = checkStart
                        filterBLAST[eachQuery][checkIndex]["already_placed"] = 1
                        checkIndex -= 1
                        continue
                    elif((checkStart <= baseStart) and (baseEnd <= checkEnd)):
                        # Base           ------     
                        # Check   ----------------
                        baseStart = checkStart
                        baseEnd = checkEnd
                        filterBLAST[eachQuery][checkIndex]["already_placed"] = 1
                        checkIndex -= 1
                        continue
                    
                refseq  = filterBLAST[eachQuery][eachIndex]["refSeqs"]
                
                if not eachQuery in intermContigCoordDict:
                    intermContigCoordDict[eachQuery] = {}
                
                if not baseStart in intermContigCoordDict[eachQuery]:
                    intermContigCoordDict[eachQuery][baseStart] = {}
                    
                if not baseEnd in intermContigCoordDict[eachQuery][baseStart]:
                    intermContigCoordDict[eachQuery][baseStart][baseEnd] = refseq
                
        for fastaID in fasSeqDict.keys():
            fasIdNoHeader = fastaID.replace(">","")
            
            if fasIdNoHeader in intermContigCoordDict.keys():
                
                seqLength    = len(fasSeqDict[fastaID])
                nucPosList   = list(fasSeqDict[fastaID])
                
                shortenedSeq = ""
                refSeq       = ""
                
                ##start searching for the information to shorten the file.
                for startB in sorted(intermContigCoordDict[fasIdNoHeader].keys()):
                    for endB in sorted(intermContigCoordDict[fasIdNoHeader][startB].keys()):
                        # Correct start and end positions as required
                        refSeq = intermContigCoordDict[fasIdNoHeader][startB][endB]
                        
                        if startB < 0:
                            startB = 0
                        
                        if endB >= seqLength:
                            endB = seqLength - 1
                            
                        #Note: Genewise (GW) positions start with 1, Blast (B) positions with 0 -> thus differenciate between start_B and start_GW
                        shortenedStartGW = len(shortenedSeq) + 1
                        count = -1
                        
                        for eachNucPos in nucPosList:
                            count += 1
                            
                            if not ((count >= startB) and (count <= endB)):
                                continue
                            
                            shortenedSeq += eachNucPos
                        
                        shortenedEndGW = len(shortenedSeq)
                        addFactor = startB + 1 - shortenedEndGW
                        
                        if not fasIdNoHeader in contigCoordDict:
                            contigCoordDict[fasIdNoHeader] = {}
                        
                        if not shortenedStartGW in contigCoordDict[fasIdNoHeader]:
                            contigCoordDict[fasIdNoHeader][shortenedStartGW] = {}
                        
                        contigCoordDict[fasIdNoHeader][shortenedStartGW][shortenedEndGW] = addFactor
                        
                        # Generate a file containing shortened sequence
                        sequenceFile = outGeneWise + "/" + fasIdNoHeader + "_sequence.txt"
                        
                        seqFileHandle = open(sequenceFile, "w")
                        seqFileHandle.write(">%s\n" % (fasIdNoHeader))
                        seqFileHandle.write(shortenedSeq + "\n")
                        
                        seqFileHandle.close()
                        
                        # Generate slight different version of file containing shortened sequence
                        #$prefix_for_qsub = "a" if ($$user_options{-c} eq "s");  # changes the name based on option to use clustered computing or not
                        
                        seqShortenedFile = outGeneWise + "/" + fasIdNoHeader + "_sequence_shortened.txt"
                        
                        seqShortenedHandle = open(seqShortenedFile, "w")
                        seqShortenedHandle.write(">%s\t%s\n" % (fasIdNoHeader, refSeq))
                        seqShortenedHandle.write(shortenedSeq + "\n")
                        seqShortenedHandle.close()
                        
                        shortenedSeqDict[seqShortenedFile] = fasIdNoHeader
    
            
        return (contigCoordDict, shortenedSeqDict)                
    ########## End of function shortSeqForGenwise ##########
        
    ########## Function runGenwise: Execute GeneWise ##########                
    def runGeneWise (self, filterBLAST, outGeneWise, shortSeqDict, hmmData):                    
                        
        gwOutFileDict = {}
        
        # First we create list of genewise output files
        for shortSeqFile in sorted(shortSeqDict.keys()):
            
            querySeq = shortSeqDict[shortSeqFile]
            
            for seqID in sorted(filterBLAST[querySeq].keys()):
                
                refSeq = filterBLAST[querySeq][seqID]["refSeqs"]
                
                geneWiseOutFile = outGeneWise + "/" + querySeq + "_" + refSeq + "_genewise.txt"
                
                if not querySeq in gwOutFileDict:
                    gwOutFileDict[querySeq] = {}
                        
                if not geneWiseOutFile in gwOutFileDict[querySeq]:
                    gwOutFileDict[querySeq][geneWiseOutFile] = (refSeq, shortSeqFile)
                
                
        # Now we create genewise command and execute them
        for querySeq in sorted(gwOutFileDict.keys()):
            for geneWiseOutFile in sorted(gwOutFileDict[querySeq].keys()):
                refShortSeq = gwOutFileDict[querySeq][geneWiseOutFile]
                
                refSeq       = refShortSeq[0]
                shortSeqFile = refShortSeq[1]
                
                geneWiseCommand = "sub_binaries/genewise %s/%s.hmm %s" % (hmmData, refSeq, shortSeqFile)
                geneWiseCommand += " -init local -quiet -gene data/genewise_support_files/human.gf -matrix data/genewise_support_files/blosum62.bla"
                geneWiseCommand += " -codon data/genewise_support_files/codon.table -hmmer -subs 0.01 -indel 0.01 -gap 11 -ext 1 -both -pep -sum"
                geneWiseCommand += " > %s" % (geneWiseOutFile)
                
                print geneWiseCommand
                os.system(geneWiseCommand)
        
        return gwOutFileDict
    ########## End of function runGenwise ##########
    
    ########## Function parseGenwise: Parse and process GeneWise output files to generate input files for HMM align ##########
    def parseGeneWise (self, geneWiseOutFiles, contigCrdDict):
        
        geneWiseSummDict =  {}
        
        for querySeq in sorted(geneWiseOutFiles.keys()):
            print querySeq
            
            geneWiseResultsRaw = {}
            geneWiseResults = {}
            minOneHit = 0
            count = 0
            
            for geneWiseOutFile in sorted(geneWiseOutFiles[querySeq].keys()):
                print ">>",geneWiseOutFile
                
                gwOutFileHandle = open(geneWiseOutFile, "rb")
                
                headerCount = 0
                seqCount = -1
            
                gwFileLines = gwOutFileHandle.readlines()
                
                for fileLine in gwFileLines:
                    fileLine = string.strip(fileLine)
                    
                    bitScore    = ""
                    refAsQuery  = ""
                    start       = ""
                    end         = ""
                        
                    if re.match("\A\d", fileLine):
                        #print fileLine
                        wiseStats = fileLine.split()
                        #print wiseStats
                        bitScore    = wiseStats[0]
                        refAsQuery  = wiseStats[1]
                        start       = wiseStats[5]
                        end         = wiseStats[6]
                        
                        if refAsQuery:
                            print "yeee hawwww"
                            minOneHit = 1
                        else:
                            print "seems like something is missing here"
                        
                        direction = "forward"
                        
                        if (start > end):
                            tempVar = start
                            start = end
                            end = tempVar
                            direction = "reverse"
                            
                        #print bitScore, refAsQuery, start, end
                        
                        #correct the positions (Genewise has been run on a shortened sequence, thus calculate the true positions)
                        
                        for startCoord in sorted(contigCrdDict[querySeq].keys()):
                            if (start >= startCoord):
                                #print "%s is greater than or equal to %s" % (str(start), str(startCoord))
                                for endCoord in sorted(contigCrdDict[querySeq][startCoord].keys()):
                                    #print "end is %s and end coord is %s" % (str(end), str(endCoord))
                                    
                                    if (end <= endCoord):
                                        #print "%s is greater than or equal to %s" % (str(end), str(endCoord))
                                        addFactor = contigCrdDict[querySeq][startCoord][endCoord]
                                        #print "addFactor: ", addFactor
                                        start += addFactor
                                        end += addFactor
                                        break
                        
                        if not querySeq in geneWiseResultsRaw:
                            geneWiseResultsRaw[querySeq] = {}
                            
                        if not geneWiseOutFile in geneWiseResultsRaw[querySeq]:
                            geneWiseResultsRaw[querySeq][geneWiseOutFile] = {}
                        
                        if not seqCount in geneWiseResultsRaw[querySeq][geneWiseOutFile]:
                            geneWiseResultsRaw[querySeq][geneWiseOutFile][headerCount] = {}
                            
                        geneWiseResultsRaw[querySeq][geneWiseOutFile][headerCount]["start"]     = start
                        geneWiseResultsRaw[querySeq][geneWiseOutFile][headerCount]["end"]       = end
                        geneWiseResultsRaw[querySeq][geneWiseOutFile][headerCount]["ref"]       = refAsQuery
                        geneWiseResultsRaw[querySeq][geneWiseOutFile][headerCount]["bitscore"]  = bitScore
                        geneWiseResultsRaw[querySeq][geneWiseOutFile][headerCount]["direction"] = direction
                        headerCount += 1
                    
                    elif re.match("\A>", fileLine):
                        seqCount += 1
                        
                        #print "---*", fileLine, "\t", str(seqCount)
                        
                        if not querySeq in geneWiseResultsRaw:
                            geneWiseResultsRaw[querySeq] = {}
                        
                        if not geneWiseOutFile in geneWiseResultsRaw[querySeq]:
                            geneWiseResultsRaw[querySeq][geneWiseOutFile] = {}
                        
                        if not seqCount in geneWiseResultsRaw[querySeq][geneWiseOutFile]:
                            geneWiseResultsRaw[querySeq][geneWiseOutFile][seqCount] = {}
                        
                        geneWiseResultsRaw[querySeq][geneWiseOutFile][seqCount]["sequence"] = ""
                    
                    elif (re.match("\A\w", fileLine) and not (re.match("\ABits", fileLine)) and not (re.match("\AMaking", fileLine))):
                        #print "--->", fileLine
                        
                        if not querySeq in geneWiseResultsRaw:
                            geneWiseResultsRaw[querySeq] = {}
                        
                        if not geneWiseOutFile in geneWiseResultsRaw[querySeq]:
                            geneWiseResultsRaw[querySeq][geneWiseOutFile] = {}
                        
                        if not seqCount in geneWiseResultsRaw[querySeq][geneWiseOutFile]:
                            geneWiseResultsRaw[querySeq][geneWiseOutFile][seqCount] = {}
                        
                        geneWiseResultsRaw[querySeq][geneWiseOutFile][seqCount]["sequence"] += fileLine
                    
                        
                gwOutFileHandle.close()        
            
            
            #done
            #do the purifying step.
            
            if not minOneHit:
                continue
            
            if querySeq is "AHHF1215_g1":
                testRef   = geneWiseResultsRaw[querySeq][geneWiseOutFileBase][baseCount]["ref"]
                print "testing testing",type(testRef)
            
            for geneWiseOutFileBase in sorted(geneWiseOutFiles[querySeq].keys()):
                for baseCount in sorted(geneWiseResultsRaw[querySeq][geneWiseOutFileBase]):
                    #print "---->>>", str(baseCount)
                    
                    baseStart = geneWiseResultsRaw[querySeq][geneWiseOutFileBase][baseCount]["start"]
                    baseEnd   = geneWiseResultsRaw[querySeq][geneWiseOutFileBase][baseCount]["end"]
                    baseRef   = geneWiseResultsRaw[querySeq][geneWiseOutFileBase][baseCount]["ref"]
                    baseBit   = geneWiseResultsRaw[querySeq][geneWiseOutFileBase][baseCount]["bitscore"]
                    baseDir   = geneWiseResultsRaw[querySeq][geneWiseOutFileBase][baseCount]["direction"]
                    baseSeq   = geneWiseResultsRaw[querySeq][geneWiseOutFileBase][baseCount]["sequence"]
                    
                    print baseStart,"\t",baseEnd,"\t",baseRef,"\t",baseBit,"\t",baseDir,"\t",baseSeq
                    
                    if baseRef is None:
                        print "baseRef is empty"
                        
                    if ((baseRef == "") and (baseStart == "") and (baseEnd == "")):
                        errorString = "ERROR: the file %s cannot be parsed!" % (geneWiseOutFileBase)
                        errorString += "Please contact the authors about it.  As a quick solution to the problem, try to remove the sequence %s, which produced this hit, from your input file." % (querySeq)
                        sys.exit(errorString)
                   
            #        unless ((defined $base_cog) && (defined $base_start) && (defined $base_end)) {
            #            my $errorstring = "ERROR: the file \"$base_genewise_outputfile\" cannot be parsed!\n";
            #            $errorstring .= "Please contact the outhors about it. As a quick solution to the problem, try to remove ";
            #            $errorstring .= "the sequence, which produced this hit, from your input file.\n";
            #            die "$errorstring";
            #    
            #        }
            #        my $base_length = $base_end - $base_start;
            #        my $is_valid = 1;
            #        
            #        foreach my $check_genewise_outputfile (sort {$a cmp $b} keys %{$genewise_results_raw{$contig}}) {
            #            foreach my $check_count (sort {$a cmp $b} keys %{$genewise_results_raw{$contig}{$check_genewise_outputfile}}) {
            #                next if ($base_count == $check_count);
            #                my $check_start = $genewise_results_raw{$contig}{$check_genewise_outputfile}{$check_count}{"start"};
            #                my $check_end = $genewise_results_raw{$contig}{$check_genewise_outputfile}{$check_count}{"end"};
            #                my $check_cog = $genewise_results_raw{$contig}{$check_genewise_outputfile}{$check_count}{"cog"};
            #                unless ((defined $check_cog) && (defined $check_start) && (defined $check_end)) {
            #                    my $errorstring = "ERROR: the file \"$check_genewise_outputfile\" cannot be parsed!\n";
            #                    $errorstring .= "Please contact the outhors about it. As a quick solution to the problem, try to remove ";
            #                    $errorstring .= "the sequence, which produced this hit, from your input file.\n";
            #                    die "$errorstring"; 
            #                }
            #                my $check_length = $check_end - $check_start;
            #                my %info = ();    
            #                $info{"base"}{"start"} = $base_start;
            #                $info{"base"}{"end"} = $base_end;
            #                $info{"check"}{"start"} = $check_start;
            #                $info{"check"}{"end"} = $check_end;
            #                                
            #                my $overlap = $small_subroutines->calculate_overlap(\%info);
            #                
            #                #ok, now we have all needed information. purify.
            #                
            #                if (($overlap / $base_length) > 0.5) {
            #                    if ($base_cog eq $check_cog) {
            #                        if ($base_length < $check_length) {
            #                            # the maior difference between the hits is the length. Keep the longer.
            #                            $is_valid = 0;
            #                        }
            #                    } elsif ($base_length < ($check_length / 2)) {
            #                        #it's not the same cog. Thus only skip this one if it is <1/2 the lenght of the other...
            #                        $is_valid = 0;
            #                    }
            #                }
            #                
            #                if ($is_valid && ($base_cog eq $check_cog)) {
            #                    #ok, there was no overlap. But i want also remove sidehits of the same COG 
            #                    if ($base_length < ($check_length * 0.7)) {
            #                        $is_valid = 0;
            #                    }                       
            #                }        
            #                #done                       
            #            }
            #        }
            #        
            #        if ($is_valid) {
            #            $genewise_results{$contig}{$count}{"start"} = $base_start;
            #            $genewise_results{$contig}{$count}{"end"} = $base_end;
            #            $genewise_results{$contig}{$count}{"cog"} = $base_cog;
            #            $genewise_results{$contig}{$count}{"direction"} = $base_direction;
            #            $genewise_results{$contig}{$count}{"sequence"} = $base_sequence;
            #            $count++;
            #        }
            #    }
            #}
            
            
       
    
    
    
    
    
    ########## End of function parseGenwise ##########
        