#!/sw/bin/python
##############################
##
## MLTreeMap TSN v. 0.0
##
## To Do List
##        BLAST
##            Auto detect AA vs Nt
##        GeneWise
##            Extracts marker genes
##            1000 bp flanking
##        hmmalign
##            Align marker genes vs reference
##        Gblocks
##            Remove minor gaps
##        RAxML
##        
##        How to deal with AA vs DNA vs rRNA
##        Visualization
##
##############################

try:
    import argparse
    import sys
    import os
    from os import path
    import shutil
    import re
    import glob
    import subprocess
    import time
except:
    print """ Could not load some user defined  module functions"""
    print """ """
    print traceback.print_exc(10)
    sys.exit(3)


def os_type():
    x = sys.platform
    if x:

        hits = re.search(r'darwin', x, re.I)
        if hits :
          return 'mac'
     
        hits = re.search(r'win', x, re.I)
        if hits :
          return 'win'

        hits = re.search(r'linux', x, re.I)
        if hits:
          return 'linux'
#
# Copy Perl's ability to autovivify
#
def pathDelim():
     ostype = os_type()
     if ostype == 'win':
         return "\\"
 
     if ostype in ['linux', 'mac']:
         return "/"

PATHDELIM =  str(pathDelim())

class Autovivify(dict):
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value


#class BookKeeping:
    #we will implement a singleton class for bookkeeping

     

def getParser(): 

    #
    # Collect options from user
    #

    parser = argparse.ArgumentParser(description='TK Takes a sequence(s) and, using a Maximum Likelihood algorithm, attempts to place it in an appropriate reference tree.')
    parser.add_argument('-i', '--input', required=True, help='your sequence input file')
    parser.add_argument('-b', '--bootstraps', default=0, type=int, help='the number of Bootstrap replicates')
    parser.add_argument('-c', '--cluster', default=0, choices=[0,'s'], help='use a computer cluster? (0 = no cluster; s = sun grid)')
    parser.add_argument('-f', '--phylogeny', default='v', choices=['v','p'], help='RAxML algorithm (v = Maximum Likelihood; p = Maximum Parsimony)')
    parser.add_argument('-g', '--gblocks', default=50, type=int, help='minimal sequence length after Gblocks')
    parser.add_argument('-l', '--filelength', default=2000, type=int, help='long input files will be splint into files containing L sequences each')
    parser.add_argument('-m', '--memory', default=0, type=int, help='minimum memory on a sungrid_cluster in GB')
    parser.add_argument('-o', '--output', default='output/', help='output directory')
    parser.add_argument('-s', '--bitscore', default=60, type=int, help='minimum bitscore for the blast hits')
    parser.add_argument('-t', '--reftree', default='p', choices=['p','g','i'], help='phylogenetic reference tree (p = MLTreeMap reference tree; g = GEBA reference tree; i = fungi tree)')
    return parser

def checkParserArguments(parser):

    #
    # Ensure files contain more than 0 sequences
    #

    args = parser.parse_args()
    if args.filelength <= 0:
        sys.exit('Input files require a positive number of sequences!')
    
    #
    # Set the reference data file prefix and the reference tree name
    #
    
    if args.reftree == 'g':
        args.reference_data_prefix = 'geba_'
        args.reference_tree = 'geba.tree'
    elif args.reftree == 'i':
        args.reference_data_prefix = 'fungi_'
        args.reference_tree = 'fungitr_tree.txt'
    else:
        args.reference_data_prefix = ''
        args.reference_tree = 'MLTreeMap_reference.tree'

    return args

#
# Prompt the user how to deal with the output directory if it already exists
#

def removePreviousOutput(args):
    while os.path.isdir(args.output):
        print('WARNING: Your output directory "' + args.output + '" already exists!')
        print('Overwrite [1], quit [2], or change directory [3]?')
        answer = raw_input()
        answer = int(answer)
        while not answer == 1 and not answer == 2 and not answer == 3:
            answer = raw_input('Invalid input. Please choose 1, 2, or 3.\n')
            answer = int(answer)
        if answer == 1:
            print('Do you really want to overwrite the old output directory?')
            print('All data in it will be lost!')
            answer2 = raw_input('Yes [y] or no [n]?\n')
            while not answer2 == 'y' and not answer2 == 'n':
                answer2 = raw_input('Invalid input. Please choose y or n.\n')
            if answer2 == 'y':
                shutil.rmtree(args.output)
            else:
                sys.exit('Exit MLTreeMap\n')
        elif answer == 2:
            sys.exit('Exit MLTreeMap\n')
        else:
            args.output = raw_input('Please enter the path to the new directory.\n')
    
    #
    # Create the output directories
    #
    
    if (not re.search('/$', args.output)):
        args.output = args.output + '/'
    os.makedirs(args.output)
    os.mkdir(args.output + 'various_outputs/')
    os.mkdir(args.output + 'final_RAxML_outputs/')
    os.mkdir(args.output + 'final_outputs/')

    return args

def createCogList(args):
    
    #
    # Create list of MLTreeMap COGs
    #
    
    cog_list = Autovivify()
    text_of_analysis_type = Autovivify()
    alignment_set = args.reftree
    kind_of_cog = ''
    
    #
    # For each line in the COG list file...
    #
    
    cogInputList = open('data/tree_data/cog_list.txt', 'r')
     
    cogList = [ x.strip() for x in cogInputList.readlines() ] 
    for cogInput in cogList:
        
        #
        # Get the kind of COG if cogInput is a header line
        #
        
        if (re.match(r'\A#(.*)', cogInput)):
            kind_of_cog = re.match(r'\A#(.*)', cogInput).group(1)
            continue
        
        #
        # Add data to COG list based on the kind of COG it is
        #
        
        if (kind_of_cog == 'phylogenetic_cogs'):
            cog_list[kind_of_cog][cogInput] = alignment_set
            cog_list['all_cogs'][cogInput] = alignment_set
            text_inset = ''
            if (alignment_set == 'g'):
                text_inset = ' based on the GEBA reference'
            if (alignment_set == 'i'):
                text_inset = ' focusing only on fungi'
            text_of_analysis_type[alignment_set] = 'Phylogenetic analysis' + text_inset + ':'
        elif (kind_of_cog == 'phylogenetic_rRNA_cogs'):
            cog, denominator, text = cogInput.split('\t')
            cog_list[kind_of_cog][cog] = denominator
            cog_list['all_cogs'][cog] = denominator
            text_of_analysis_type[denominator] = 'Phylogenetic analysis, ' + text + ':'
        elif (kind_of_cog == 'functional_cogs'):
            cog, denominator, text = cogInput.split('\t')
            cog_list[kind_of_cog][cog] = denominator
            cog_list['all_cogs'][cog] = denominator
            text_of_analysis_type[denominator] = 'Functional analysis, ' + text + ':'
            
    #
    # Close the COG list file
    #
    
    cogInputList.close()

    return (cog_list, text_of_analysis_type)

def splitFastaInput(args):

    #
    # Confirm input file is a fasta file
    #

    input = open(args.input, 'r')
    
    if (not input.read(1) == '>'):
        sys.exit('ERROR: Your file does not appear to be a proper FASTA file!\n')
    
    #
    # Unread the '>' to prevent problems later
    #
    
    input.seek(-1,1)
    
    #
    # Determine the output file names
    # Open the output files
    #
    
    if (re.match(r'\A.*\/(.*)', args.input)):
        inputFileName = re.match(r'\A.*\/(.*)', args.input).group(1)
    else:
        inputFileName = args.input
    
    outputSplit = open(args.output + PATHDELIM + 'various_outputs' + PATHDELIM + inputFileName + '_0.txt', 'w')
    outputFormatted = open(args.output + PATHDELIM +  'various_outputs' + PATHDELIM + inputFileName + '_formatted.txt', 'w')
    args.formatted_input_file = args.output + PATHDELIM +  'various_outputs' + PATHDELIM + inputFileName + '_formatted.txt'
    args.output_directory_var = args.output + PATHDELIM + 'various_outputs' 
    countFiles = 0
    countSequences = 0
    
    #
    # Iterate through the input file...
    #
    
    countTotal = 0
    countNucleotides = 0
    countXN = 0
    countUndef = 0
    splitFiles = []
    
    for line in input:
    
        if (re.search('\A>', line)):
        
            countSequences += 1
        
            #
            # Replace all non a-z, A-Z, 0-9, or . characters with a _
            # Then replace the initial _ with a >
            #
        
            line = re.sub(r'[^a-zA-Z0-9.\r\n]', '_', line)
            line = re.sub(r'\A_', '>', line)
            
            #
            # RAxML can only work with file names having length <= 125
            # Check that the sequence name length is <= 100
            # If sequence name length is > 100, limit the file name length to 100
            #
            
            if (line.__len__() > 100):
                line = line[0:100]
        
            #
            # Split the file if countSequences > the number of sequences per file specified by the user
            #
            
            if (countSequences >= args.filelength):
               countSequences = 0
               splitFiles.append(args.output + PATHDELIM + 'various_outputs' + PATHDELIM + inputFileName + '_%d.txt' %(countFiles))
               countFiles += 1
               outputSplit.close()
               outputSplit = open(args.output + PATHDELIM +  'various_outputs' + PATHDELIM + inputFileName + '_%d.txt' %(countFiles), 'w')
        else:
        
            #
            # Remove all non-characters from the sequence
            #
            
            re.sub(r'[^a-zA-Z]','', line)
            
            #
            # Count the number of {atcg} and {xn} in all the sequences
            #
    
            characters = []
            characters = list(line)
            
            for character in characters:
                countTotal += 1
                if (re.match(r'[acgtACGT]', character)):
                    countNucleotides += 1
                elif (re.match(r'[xnXN]', character)):
                    countXN += 1
                else:
                    countUndef += 1
        
        #
        # Write the lines to the appropriate files
        #
        
        outputSplit.write(line)
        outputFormatted.write(line)
    
    #
    # Close the files
    #
    
    input.close()
    outputSplit.close()
    outputFormatted.close()
        
    #
    # If splitFiles is empty, add the only file to splitFiles
    #
    
    if not splitFiles:
        splitFiles.append(args.output + 'various_outputs/' + inputFileName + '_%d.txt' %(countFiles))
    
    #
    # Exit the program if character count is 0
    #
    
    if (countTotal == 0):
        sys.exit('ERROR: Your input file appears to be corrupted. No sequences were found!\n')
    
    #
    # Exit the program if all sequences are composed only of X or N
    #
    
    elif (countXN == countTotal):
        sys.exit('ERROR: Your sequence(s) contain only X or N!\n')
    
    #
    # Exit the program if less than half of the characters are nucleotides
    # The 1.0 is to cast it as a float
    #
    
    elif (float(countNucleotides / (countTotal * 1.0)) < 0.5):
        sys.exit('ERROR: Your sequence(s) most likely contain no DNA!\n')

    return splitFiles

def createBlastDBList(args):

    #
    # Create list of databases for each blastx and blastn
    #

    blastxDB = []
    blastnDB = []
    
    for file in glob.glob('data/' + args.reference_data_prefix + 'alignment_data/*.fa'):
        file.rstrip('\r\n')
        if (not re.match(r'\Adata/' + args.reference_data_prefix + 'alignment_data/\._', file)):
            if (re.match(r'.*rRNA\.fa\Z', file)):
                blastnDB.append(file)
            elif (re.match(r'.*\.fa\Z', file) and not re.match(r'rRNA', file)):
                blastxDB.append(file)

    return (blastxDB, blastnDB)

def runBlast(args, splitFiles, blastxDB, blastnDB):
    
    #
    # For each file containing a maximum of the specified number of sequences...
    #

    for splitFile in splitFiles:
        
        #
        # Ensure splitFile is a .txt file; save file name if so, die otherwise
        #
        
        blastInputFileName = ''
        
        if (not re.match(r'\A.+/(.+)\.txt\Z', splitFile)):
            sys.exit('ERROR: Something is wrong with the directory of the BLAST input file!\n')
        else:
            blastInputFileName = re.match(r'\A.+/(.+)\.txt\Z', splitFile).group(1)
    
        #
        # BLAST splitFile against each blastx DB
        #

        command = 'sub_binaries/blastall -p blastx -i ' + splitFile + ' -M BLOSUM62 -d "'
        for db in blastxDB:
             command += db + ' '
        command += '" -e 0.01 -v 20000 -b 20000 -z 1000000 -m 8 '
        command += '> ' + args.output + 'various_outputs/' + blastInputFileName + '.BLAST_results_raw.txt'
        os.system(command)
        
        #
        # BLAST splitFile against each blastn DB
        #
        
        command = 'sub_binaries/blastall -p blastn -i ' + splitFile + ' -M BLOSUM62 -d "'
        for db in blastnDB:
            command += db + ' '
        command += '" -e 0.01 -v 20000 -b 20000 -z 1000000 -m 8 '
        command += '> ' + args.output + 'various_outputs/' + blastInputFileName + '.rRNA_BLAST_results_raw.txt'
        os.system(command)
        
        #
        # Remove the BLAST input file
        #
        
        if path.exists(splitFile):
           os.remove(splitFile)
        
    #
    # Remove empty BLAST result raw files; store non-empty files in a list
    #
    
def readBlastResults(output):
    rawBlastResultFiles = []
    
    for file in glob.glob(output + '/various_outputs/*BLAST_results_raw.txt'):
        file.rstrip('\r\n')
        if path.getsize(file) <= 0:
            shutil.rmtree(file)
        else:
            rawBlastResultFiles.append(file)
    
    return rawBlastResultFiles

def parseBlastResults(args, rawBlastResultFiles, cog_list):

    counter=0
    purifiedBlastHits =Autovivify()
    for file in sorted(rawBlastResultFiles):
        try:     
           blastResults = open(file, 'r')
        except IOError:
           print "ERROR: Cannot open BLAST outputfile " + file
           continue

        contigs = Autovivify()
        identifier = 0
        for line in blastResults:
            # Clear the variables referencing the contig, COG, query start, query end, reference start, reference end, and bitscore
            # Interpret the BLAST hit, and assign the details accordingly
            tempContig, tempDetailedCOG, _, _, _, _, tempQStart, tempQEnd, tempRStart, tempREnd, _, tempBitScore = line.split('\t')
            tempREnd = int(tempREnd)
            tempRStart = int(tempRStart)
            tempQEnd = int(tempQEnd)
            tempQStart = int(tempQStart)
            tempBitScore = float(tempBitScore)

            # Skip to next BLAST hit if bit score is less than user-defined minimum
            if (tempBitScore <= args.bitscore):
                continue
            # Determine the direction of the hit relative to the reference
            #
            direction = 'forward'
            if tempRStart > tempREnd:
                temp = tempRStart
                tempRStart = tempREnd
                tempREnd = temp
                direction = 'reverse'
            if tempQStart > tempQEnd:
                temp = tempQStart
                tempQStart = tempQEnd
                tempQEnd = temp
                if (direction == 'reverse'):
                    sys.exit('ERROR: Parsing error with the BLAST results. Please notify the authors about ' + tempContig + ' at ' + tempDetailedCOG + 'q('+tempQEnd+'..'+tempQStart+'),r('+tempREnd+'..'+tempRStart+')')
                direction = 'reverse'
            
            # Trim COG name to last 7 characters of detailed COG name
            # TK - This will be important to note in the user's manual, especially if we enable people to add their own COGs later
            #
            if re.match(r'.*(.{7})\Z', tempDetailedCOG):
                tempCOG = re.match(r'.*(.{7})\Z', tempDetailedCOG).group(1)
            else:
                sys.exit('ERROR: Could not detect the COG of sequence ' + tempDetailedCOG)
            

            #
            # Save contig details to the list
            #
            contigs[tempContig][identifier]['bitscore'] = tempBitScore
            contigs[tempContig][identifier]['cog'] = tempCOG
            contigs[tempContig][identifier]['seq_start'] = tempQStart
            contigs[tempContig][identifier]['seq_end'] = tempQEnd
            contigs[tempContig][identifier]['direction'] = direction
            contigs[tempContig][identifier]['validity'] = True
            identifier += 1
        #
        # Close the file
        #
        blastResults.close()
        
        #
        # Purify the BLAST hits
        #
 
        #
        # For each contig sorted by their stringwise comparison...
        #
        for contig in sorted(contigs.keys()):
            identifier = 0
            #
            # For each blast result for that contig...
            print contig + '------------------------>'
            for base_blast_result_raw_identifier in sorted(contigs[contig].keys()):
                base_bitscore = contigs[contig][base_blast_result_raw_identifier]['bitscore']
                base_cog = contigs[contig][base_blast_result_raw_identifier]['cog']
                base_start = contigs[contig][base_blast_result_raw_identifier]['seq_start']
                base_end = contigs[contig][base_blast_result_raw_identifier]['seq_end']
                direction = contigs[contig][base_blast_result_raw_identifier]['direction']
                base_length = base_end - base_start # TK Why not +1?
                
                # Skip if base_bitscore is less than user specified minimum bitscore
                if (base_bitscore < args.bitscore):
                    continue
                #
                # Set validity to 0 if COG is not in list of MLTreeMap COGs
                if not base_cog in cog_list['all_cogs']:
                    contigs[contig][base_blast_result_raw_identifier]['validity'] = False
                
                #
                # Compare the BLAST hit (base) against all others
                # There may be several opinions about how to do this. This way is based on the original MLTreeMap
                # ----A----  --C--
                #        ---B---
                # A kills B, B kills C. (Another approach would be to let C live, but the original MLTreeMap authors don't expect C to be useful)
                #
                
                for check_blast_result_raw_identifier in sorted(contigs[contig]):
                    check_bitscore = contigs[contig][check_blast_result_raw_identifier]['bitscore']
                    check_cog = contigs[contig][check_blast_result_raw_identifier]['cog']
                    check_start = contigs[contig][check_blast_result_raw_identifier]['seq_start']
                    check_end = contigs[contig][check_blast_result_raw_identifier]['seq_end']
                    check_length = check_end - check_start # TK Why not +1?
                    
                    # Don't compare base hit against itself; skip to next iteration
                    if base_blast_result_raw_identifier == check_blast_result_raw_identifier:
                        continue
                    
                    #
                    # Compare the base and check BLAST hits
                    #
                    overlap = 0
                    
                    if base_start <= check_start and check_start <= base_end and base_end <= check_end:
                        overlap = base_end - check_start # TK Why not +1?
                    elif base_start <= check_start and check_end <= base_end:
                        # Base     --------
                        # Check        --
                        overlap = check_end - check_start # TK Why not +1?
                    elif check_start <= base_start and base_start <= check_end and check_end <= base_end:
                        # Base         -----
                        # Check    -----
                        overlap = check_end - base_start # TK Why not +1?
                    elif check_start <= base_start and base_end <= check_end:
                        # Base       --
                        # Check    --------
                        overlap = base_end - base_start #TK Why not +1?
                    counter +=1
                    #
                    # Check for validity for hits with overlap
                    #
                    if overlap > 0:
                        if overlap  > 0.5*base_length and base_bitscore < check_bitscore:
                            contigs[contig][base_blast_result_raw_identifier]['validity'] = False
                        elif overlap > 0.5*check_length and check_bitscore < base_bitscore:
                            contigs[contig][check_blast_result_raw_identifier]['validity'] = False
                        elif base_start == check_start and base_end == check_end:
                            # If both are the same, keep only the one with the smaller identifier
                           if check_blast_result_raw_identifier > base_blast_result_raw_identifier:
                                contigs[contig][check_blast_result_raw_identifier]['validity'] = False
                           else:
                                contigs[contig][base_blast_result_raw_identifier]['validity'] = False
                    
                    #
                    # Save purified hits for valid base hits
                    
                if contigs[contig][base_blast_result_raw_identifier]['validity']:
                     print contig + ' ' + str(base_blast_result_raw_identifier) + '  ' + str(base_bitscore)

                     purifiedBlastHits[contig][identifier]['bitscore'] = base_bitscore
                     purifiedBlastHits[contig][identifier]['cog'] = base_cog
                     purifiedBlastHits[contig][identifier]['start'] = base_start
                     purifiedBlastHits[contig][identifier]['end'] = base_end
                     purifiedBlastHits[contig][identifier]['direction'] = direction
                     purifiedBlastHits[contig][identifier]['is_already_placed'] = False
                     identifier += 1
    
    #
    # Print the BLAST results for each contig
    #
    print "purified list " + str(len(purifiedBlastHits.keys()))
    print '\n'.join(sorted(purifiedBlastHits.keys()))
    for contig in sorted(purifiedBlastHits.keys()):
        outfile = args.output + PATHDELIM + 'various_outputs' +  PATHDELIM + contig + '_blast_result_purified.txt'
        out = open(outfile, 'w')
        sorting_hash = {}

        for identifier in sorted(purifiedBlastHits[contig].keys()):
            if not purifiedBlastHits[contig][identifier]['bitscore'] in sorting_hash:
               sorting_hash[purifiedBlastHits[contig][identifier]['bitscore']] = {}
            sorting_hash[purifiedBlastHits[contig][identifier]['bitscore']][identifier] = 1

        for bitscore in sorted(sorting_hash.keys(), reverse=True):
            for identifier in sorted(sorting_hash[bitscore]):
                out.write(contig + '\t' + str(purifiedBlastHits[contig][identifier]['start']) + '\t' +\
                str(purifiedBlastHits[contig][identifier]['end']) + '\t' +\
                str(purifiedBlastHits[contig][identifier]['direction']) + '\t' +\
                purifiedBlastHits[contig][identifier]['cog'] + '\t' + str( bitscore) + '\n')
        out.close()
    return purifiedBlastHits

def produceGenewiseFiles(args, blast_hits_purified):
    flanking_length = 1000; # Recommended: 1000
    prae_contig_coordinates = Autovivify()
    contig_coordinates = Autovivify()
    shortened_sequence_files = {}

    for contig in sorted(blast_hits_purified.keys()):
        for base_identifier in sorted(blast_hits_purified[contig].keys()):
            #
            # Skip rRNA hits for now (we work with them later)
            #
            if re.search("rRNA", blast_hits_purified[contig][base_identifier]['cog']):
                continue
            #
            # Skip hits which have already been placed; otherwise, mark them as placed
            #
            if blast_hits_purified[contig][base_identifier]['is_already_placed']:
                continue

            blast_hits_purified[contig][base_identifier]['is_already_placed'] = True
            base_start = blast_hits_purified[contig][base_identifier]['start'] - flanking_length
            base_end = blast_hits_purified[contig][base_identifier]['end'] + flanking_length
            nr_of_blast_hits = len(blast_hits_purified[contig].keys())
            check_identifier =0
            while check_identifier < nr_of_blast_hits:
                #
                # Skip rRNA hits for now (we work with them later)
                #
                if re.search(r'rRNA', blast_hits_purified[contig][check_identifier]['cog']):
                    check_identifier +=1
                    continue
                #
                # Skip hits which have already been placed; otherwise, mark them as placed
                #
                if blast_hits_purified[contig][check_identifier]['is_already_placed']:
                    check_identifier +=1
                    continue

                check_start = blast_hits_purified[contig][check_identifier]['start'] - flanking_length
                check_end = blast_hits_purified[contig][check_identifier]['end'] + flanking_length
                #
                # Check for overlap
                #
                if base_start <= check_start and check_start <= base_end and base_end <= check_end:
                    # Base  --------
                    # Check     --------
                    base_end = check_end
                    blast_hits_purified[contig][check_identifier]['is_already_placed'] = True
                    check_identifier = 0
                    continue
                elif base_start <= check_start and check_end <= base_end:
                    # Base  --------
                    # Check   ----
                    blast_hits_purified[contig][check_identifier]['is_already_placed'] = True
                    check_identifier = 0
                    continue
                elif check_start <= base_start and base_start <= check_end and check_end <= base_end:
                    # Base      --------
                    # Check --------
                    base_start = check_start
                    blast_hits_purified[contig][check_identifier]['is_already_placed'] = True
                    check_identifier = 0
                    continue
                elif check_start <= base_start and base_end <= check_end:
                    # Base    ----
                    # Check --------
                    base_start = check_start
                    base_end = check_end
                    blast_hits_purified[contig][check_identifier]['is_already_placed'] = True
                    check_identifier = 0
                    continue
                check_identifier += 1

            prae_contig_coordinates[contig][base_start][base_end] = 1
    #
    # Produce the input files for Genewise
    #
    input = open(args.formatted_input_file, 'r')
    contig_name = ''
    sequence = ''

    line = 'x'
    while line:
        line= input.readline()
        line =  line.strip()
        line = re.sub(r'\s', '_', line)
        searchmatch =re.search(r'\A>(.+)', line)

        if searchmatch or not line:
            if not line:
               sequence += line
            if contig_name in prae_contig_coordinates:
                sequence_length = len(sequence)
                shortened_sequence="" 
                #start searching for the information to shorten the file.
                for start_B in sorted(prae_contig_coordinates[contig_name].keys()) :
                    for end_B in sorted(prae_contig_coordinates[contig_name][start_B].keys()) :
                         #ok, now we have all information about the hit. Correct start and end if needed: 
                        if start_B < 0:
                           start_B = 0 

                        if end_B >= sequence_length:
                           end_B = sequence_length -1 
      
                        #Note: Genewise (GW) positions start with 1, Blast (B) positions with 0 -> thus differenciate between start_B and start_GW
                        shortened_start_GW = len(shortened_sequence) + 1; 
                        count = -1;
                        for nucleotide in sequence: 
                            count += 1     
                            if not (count >= start_B and count <= end_B):
                               continue
                            shortened_sequence += nucleotide;
      
                        shortened_end_GW = len(shortened_sequence)
                        addition_factor = (start_B + 1) - shortened_start_GW #$start_B + 1 == $start_GW
                        contig_coordinates[contig_name][shortened_start_GW][shortened_end_GW] = addition_factor
        
        
                try:
                    with open(args.output_directory_var + PATHDELIM + contig_name + "_sequence.txt", 'w') as f:
                       fprintf(f, "%s\n", ">"+ contig_name + "\n" + sequence)
                    f.close()
                except:
                    print  "ERROR: Can't create " + args.output_directory_var + PATHDELIM + contig_name + "_sequence.txt!"; 


                try:   
                   with open(args.output_directory_var + PATHDELIM + contig_name + "_sequence_shortened.txt", 'w') as f:
                      fprintf(f, "%s\n",">" + contig_name + "\n" + shortened_sequence)
                   f.close()
                   shortened_sequence_files[args.output_directory_var + PATHDELIM +  contig_name + "_sequence_shortened.txt"]=contig_name
                except:
                   print "ERROR: Can't create " + args.output_directory_var + PATHDELIM +  contig_name +"_sequence_shortened.txt!"; 

            if searchmatch:
               contig_name = searchmatch.group(1)
            sequence = ""
        else:
            sequence += line
    input.close()
    return contig_coordinates, shortened_sequence_files



def fprintf(file, fmt, *args): 
   file.write(fmt % args)


def main(argv):
    parser = getParser()
    args = checkParserArguments(parser)

    #removePreviousOutput(args)
    
    #this creates the list of the marker COGs 
    (cogList, textOfAnalysisType) = createCogList(args)

    #splits the input files to smaller files for blasting
    splitFiles = splitFastaInput(args)

    # get the appropriate type of blast DBS
    (blastxDB, blastnDB) = createBlastDBList(args)

    print('Run BLAST')
    #runBlast(args, splitFiles, blastxDB, blastnDB)

    blastResults =  readBlastResults(args.output)

    print blastResults
    print 'types of cogs'
    print cogList.keys()

    print 'all_cogs ' + str(len(cogList['all_cogs'] ))
    print 'functional_cogs ' + str(len(cogList['functional_cogs'] ))
    print 'functional_cogs ', cogList['functional_cogs'] 

    print 'phylogenetic_rRNA_cogs ' + str( len(cogList['phylogenetic_rRNA_cogs']))
    print 'phylogenetic_rRNA_cogs ',cogList['phylogenetic_rRNA_cogs']

    purified_blast_results = parseBlastResults(args, blastResults, cogList)

    contig_coordinates, shortened_sequence_files = produceGenewiseFiles(args, purified_blast_results)
    print contig_coordinates
if __name__ == "__main__":
   main(sys.argv[1:])

