/*  RAxML-VI-HPC (version 2.2) a program for sequential and parallel estimation of phylogenetic trees 
 *  Copyright August 2006 by Alexandros Stamatakis
 *
 *  Partially derived from
 *  fastDNAml, a program for estimation of phylogenetic trees from sequences by Gary J. Olsen
 *  
 *  and 
 *
 *  Programs of the PHYLIP package by Joe Felsenstein.
 *
 *  This program is free software; you may redistribute it and/or modify its
 *  under the terms of the GNU General Public License as published by the Free
 *  Software Foundation; either version 2 of the License, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *  for more details.
 * 
 *
 *  For any other enquiries send an Email to Alexandros Stamatakis
 *  Alexandros.Stamatakis@epfl.ch
 *
 *  When publishing work that is based on the results from RAxML-VI-HPC please cite:
 *
 *  Alexandros Stamatakis:"RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses with thousands of taxa and mixed models". 
 *  Bioinformatics 2006; doi: 10.1093/bioinformatics/btl446
 */


#ifdef PARALLEL
int numOfWorkers;
#endif

int processID;
infoList iList;
FILE   *INFILE;

int Thorough = 0;

char run_id[128] = "", 
  workdir[1024] = "", 
  seq_file[1024] = "", 
  tree_file[1024]="", 
  weightFileName[1024] = "", 
  modelFileName[1024] = "", 
  excludeFileName[1024] = "",
  bootStrapFile[1024] = "", 
  permFileName[1024] = "", 
  resultFileName[1024] = "", 
  logFileName[1024] = "", 
  checkpointFileName[1024] = "", 
  infoFileName[1024] = "", 
  randomFileName[1024] = "",   
  bootstrapFileName[1024] = "", 
  bipartitionsFileName[1024] = "",
  ratesFileName[1024] = "", 
  perSiteLLsFileName[1024] = "", 
  lengthFileName[1024] = "", 
  lengthFileNameModel[1024] = "",
  proteinModelFileName[1024] = "",
  secondaryStructureFileName[1024] = "";

char *protModels[12] = {"DAYHOFF", "DCMUT", "JTT", "MTREV", "WAG", "RTREV", "CPREV", "VT", "BLOSUM62", "MTMAM", "LG", "GTR"};

char inverseMeaningBINARY[4] = {'_', '0', '1', '-'};
char inverseMeaningDNA[16]   = {'_', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', '-'};
char inverseMeaningPROT[23]  = {'A','R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 
			       'T', 'W', 'Y', 'V', 'B', 'Z', '-'};

const unsigned int bitVectorSecondary[256] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 
					      10, 11, 12, 13, 14, 15, 0, 16, 32, 48, 64, 80, 96, 112, 128, 144, 160, 176, 192, 
					      208, 224, 240, 0, 17, 34, 51, 68, 85, 102, 119, 136, 153, 170, 187, 204, 221, 238, 
					      255, 0, 256, 512, 768, 1024, 1280, 1536, 1792, 2048, 2304, 2560, 2816, 3072, 3328, 
					      3584, 3840, 0, 257, 514, 771, 1028, 1285, 1542, 1799, 2056, 2313, 2570, 2827, 3084, 
					      3341, 3598, 3855, 0, 272, 544, 816, 1088, 1360, 1632, 1904, 2176, 2448, 2720, 2992, 
					      3264, 3536, 3808, 4080, 0, 273, 546, 819, 1092, 1365, 1638, 1911, 2184, 2457, 2730, 
					      3003, 3276, 3549, 3822, 4095, 0, 4096, 8192, 12288, 16384, 20480, 24576, 28672, 32768, 
					      36864, 40960, 45056, 49152, 53248, 57344, 61440, 0, 4097, 8194, 12291, 16388, 20485, 24582, 
					      28679, 32776, 36873, 40970, 45067, 49164, 53261, 57358, 61455, 0, 4112, 8224, 12336, 16448, 
					      20560, 24672, 28784, 32896, 37008, 41120, 45232, 49344, 53456, 57568, 61680, 0, 4113, 8226, 
					      12339, 16452, 20565, 24678, 28791, 32904, 37017, 41130, 45243, 49356, 53469, 57582, 61695, 
					      0, 4352, 8704, 13056, 17408, 21760, 26112, 30464, 34816, 39168, 43520, 47872, 52224, 56576, 
					      60928, 65280, 0, 4353, 8706, 13059, 17412, 21765, 26118, 30471, 34824, 39177, 43530, 47883, 
					      52236, 56589, 60942, 65295, 0, 4368, 8736, 13104, 17472, 21840, 26208, 30576, 34944, 39312, 
					      43680, 48048, 52416, 56784, 61152, 65520, 0, 4369, 8738, 13107, 17476, 21845, 26214, 30583, 
					      34952, 39321, 43690, 48059, 52428, 56797, 61166, 65535};

const unsigned int mask32[32] = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 
					262144, 524288, 1048576, 2097152, 4194304, 8388608, 16777216, 33554432, 67108864, 134217728, 
					268435456, 536870912, 1073741824, 2147483648};

double masterTime;
int partCount = 0;
int optimizeRateCategoryInvocations = 1;


#ifdef _USE_OMP
volatile int             NumberOfThreads;
#endif


#ifdef _USE_PTHREADS
volatile int             jobCycle;
volatile int             threadJob;
volatile int             NumberOfThreads;
volatile double          *reductionBuffer;
volatile double          *reductionBufferTwo;
volatile double          *reductionBufferThree;
volatile int             *reductionBufferParsimony;
volatile int             *barrierBuffer;
#endif
