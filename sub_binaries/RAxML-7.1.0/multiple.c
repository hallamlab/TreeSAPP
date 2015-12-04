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


#ifndef WIN32
#include <unistd.h>
#endif

#include <math.h>
#include <time.h> 
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include "axml.h"


extern int  optimizeRatesInvocations;
extern int  optimizeRateCategoryInvocations;
extern int  optimizeAlphaInvocations;
extern int  optimizeInvarInvocations;
extern int  checkPointCounter;
extern int  Thorough;
extern int  partCount;
extern char tree_file[1024];
extern const int mask32[32];
extern double masterTime;

extern FILE   *INFILE, *permutationFile, *logFile, *infoFile;

extern char seq_file[1024];
extern char permFileName[1024], resultFileName[1024], 
  logFileName[1024], checkpointFileName[1024], infoFileName[1024], run_id[128], workdir[1024], bootStrapFile[1024], bootstrapFileName[1024], bipartitionsFileName[1024]; 




void catToGamma(tree *tr, analdef *adef)
{
  assert(tr->rateHetModel == CAT);  
  
  if(adef->useInvariant)
    tr->rateHetModel = GAMMA_I;
  else
    tr->rateHetModel = GAMMA;

  
#ifdef _USE_PTHREADS
  masterBarrier(THREAD_CAT_TO_GAMMA, tr); 
#endif
}

void gammaToCat(tree *tr)
{

  assert(tr->rateHetModel == GAMMA || tr->rateHetModel == GAMMA_I);
   
  tr->rateHetModel = CAT;

 
#ifdef _USE_PTHREADS
  masterBarrier(THREAD_GAMMA_TO_CAT, tr); 
#endif  
}


static void singleBootstrap(tree *tr, int i, analdef *adef, rawdata *rdta, cruncheddata *cdta)
{
  tr->treeID = i;
  tr->checkPointCounter = 0;
     
  computeNextReplicate(tr, &adef->boot, (int*)NULL, (int*)NULL, FALSE);
  
  initModel(tr, rdta, cdta, adef);                                   
         
  getStartingTree(tr, adef);

  computeBIGRAPID(tr, adef, TRUE);

  if(adef->bootstrapBranchLengths)
    {
      switch(tr->rateHetModel)
	{
	case GAMMA:
	case GAMMA_I:      
	  modOpt(tr, adef, TRUE, adef->likelihoodEpsilon);		      	    	    
	  break;
	case CAT:    	      	 		  
	  tr->likelihood = unlikely;	       
	  catToGamma(tr, adef);	  
	  initModel(tr, rdta, cdta, adef);	  	  	  
	  modOpt(tr, adef, TRUE, adef->likelihoodEpsilon);	 
	  gammaToCat(tr); 		  	      	        	
	  break;
	default:
	  assert(0);
	}
    }  

#ifndef PARALLEL
   printBootstrapResult(tr, adef, TRUE);
#endif     
}





/***************************** EXPERIMENTAL FUNCTIONS ********************************************************/




static int compareTopolRell(const void *p1, const void *p2)
{
  topolRELL **rc1 = (topolRELL **)p1;
  topolRELL **rc2 = (topolRELL **)p2;
 
  double i = (*rc1)->likelihood;
  double j = (*rc2)->likelihood;
  
  if (i > j)
    return (-1);
  if (i < j)
    return (1);
  return (0);
}


void fixModelIndices(tree *tr, int endsite)
{
  int model, i;

  assert(tr->NumberOfModels > 0);   

  tr->partitionData[0].lower = 0;
     
  model = tr->model[0];
  i = 1;

  while(i < endsite)
    {
      if(tr->model[i] != model)
	{	      
	  tr->partitionData[model].upper = i;
	  tr->partitionData[model + 1].lower = i;
	  model = tr->model[i];
	}
      i++;
    }       
  
  tr->partitionData[tr->NumberOfModels - 1].upper = endsite;
  
  for(model = 0; model < tr->NumberOfModels; model++)    
    tr->partitionData[model].width = tr->partitionData[model].upper -  tr->partitionData[model].lower;
 
  
#ifndef _USE_PTHREADS
  for(model = 0; model < tr->NumberOfModels; model++)
    {
      int 
	j,
	lower =  tr->partitionData[model].lower;
      
      /* SOS what about sumBuffer? */
      /* tr->partitionData[model].sumBuffer    = &tr->sumBuffer[offset]; */
      tr->partitionData[model].perSiteLL    = &tr->perSiteLL[lower];
      tr->partitionData[model].wr           = &tr->cdta->wr[lower];
      tr->partitionData[model].wr2          = &tr->cdta->wr2[lower];
      tr->partitionData[model].wgt          = &tr->cdta->aliaswgt[lower];
      tr->partitionData[model].invariant    = &tr->invariant[lower];
      tr->partitionData[model].rateCategory = &tr->cdta->rateCategory[lower];
      

      for(j = 1; j <= tr->mxtips; j++)
	tr->partitionData[model].yVector[j] = &(tr->yVector[j][tr->partitionData[model].lower]);

    }
#else
  masterBarrier(THREAD_FIX_MODEL_INDICES, tr);
#endif
}

void reductionCleanup(tree *tr, int *originalRateCategories, int *originalInvariant)
{
  int j;

  tr->cdta->endsite = tr->originalCrunchedLength;

  memcpy(tr->cdta->aliaswgt, tr->originalWeights, sizeof(int) * tr->cdta->endsite);
  memcpy(tr->model, tr->originalModel, sizeof(int) * tr->cdta->endsite);
  memcpy(tr->dataVector, tr->originalDataVector,  sizeof(int) * tr->cdta->endsite);

  memcpy(tr->cdta->rateCategory, originalRateCategories, sizeof(int) * tr->cdta->endsite);
  memcpy(tr->invariant,          originalInvariant,      sizeof(int) * tr->cdta->endsite); 
  
  for (j = 0; j < tr->originalCrunchedLength; j++) 
    {	
      double temp, wtemp;     
      temp = tr->cdta->patrat[originalRateCategories[j]];
      tr->cdta->wr[j]  = wtemp = temp * tr->cdta->aliaswgt[j];
      tr->cdta->wr2[j] = temp * wtemp;
    }                           
      
  memcpy(tr->rdta->y0, tr->rdta->yBUF, tr->rdta->numsp * tr->cdta->endsite * sizeof(char));  
      
  tr->cdta->endsite = tr->originalCrunchedLength;
  fixModelIndices(tr, tr->originalCrunchedLength);      
}






void computeNextReplicate(tree *tr, long *randomSeed, int *originalRateCategories, int *originalInvariant, boolean isRapid)
{ 
  int pos, nonzero, j, model, w;   
  int *weightBuffer, endsite;                
  int *weights, i, l;  
#ifdef PARALLEL
  long seed;
#endif 

  for(j = 0; j < tr->originalCrunchedLength; j++)
    tr->cdta->aliaswgt[j] = 0;

      	  
  for(model = 0; model < tr->NumberOfModels; model++)
    {
      nonzero = 0;        	 
      
      for (j = 0; j < tr->originalCrunchedLength; j++)  
	{
	  if(tr->originalModel[j] == model)
	    nonzero += tr->originalWeights[j];
	}				          
      
      weightBuffer = (int *)calloc(nonzero, sizeof(int));	 
#ifdef PARALLEL 
      seed = (long) gettimeSrand();
      for (j = 0; j < nonzero; j++)
	weightBuffer[(int) (nonzero*randum(& seed))]++;
#else       
      for (j = 0; j < nonzero; j++)
	weightBuffer[(int) (nonzero*randum(randomSeed))]++;                  
#endif
      
      pos = 0;	      
      
      for(j = 0; j < tr->originalCrunchedLength; j++) 
	{
	  if(model == tr->originalModel[j])
	    {
	      for(w = 0; w < tr->originalWeights[j]; w++)	  	  	 
		{
		  tr->cdta->aliaswgt[j] += weightBuffer[pos];
		  pos++;		      
		}				   
	    }
	}  

      free(weightBuffer);	  
      
    }       

  endsite = 0;
  
  for (j = 0; j < tr->originalCrunchedLength; j++) 
    {	      
      if(tr->cdta->aliaswgt[j] > 0)
	endsite++;
      if(isRapid)
	{
	  double temp, wtemp;
	  temp = tr->cdta->patrat[originalRateCategories[j]];
	  tr->cdta->wr[j]  = wtemp = temp * tr->cdta->aliaswgt[j];
	  tr->cdta->wr2[j] = temp * wtemp;
	}
    }          
  
  weights = tr->cdta->aliaswgt;

  for(i = 0; i < tr->rdta->numsp; i++)
    {     
      unsigned char *yPos    = &(tr->rdta->y0[tr->originalCrunchedLength * i]);
      unsigned char *origSeq = &(tr->rdta->yBUF[tr->originalCrunchedLength * i]);
      int l, j;
      
      for(j = 0, l = 0; j < tr->originalCrunchedLength; j++)      
	if(tr->cdta->aliaswgt[j] > 0)	  	    
	  yPos[l++] = origSeq[j];	                   
    }

  for(j = 0, l = 0; j < tr->originalCrunchedLength; j++)
    {     
      if(weights[j])	
	{
	  tr->cdta->aliaswgt[l]     = tr->cdta->aliaswgt[j];
	  tr->dataVector[l]         = tr->originalDataVector[j];
	  tr->model[l]              = tr->originalModel[j];

	  if(isRapid)
	    {
	      tr->cdta->wr[l]           = tr->cdta->wr[j];
	      tr->cdta->wr2[l]          = tr->cdta->wr2[j];	 
	      tr->cdta->rateCategory[l] = originalRateCategories[j];
	      tr->invariant[l]          = originalInvariant[j];
	    }
	  l++;
	}
    }

  tr->cdta->endsite = endsite;
  fixModelIndices(tr, endsite);
}


  




static pInfo *allocParams(tree *tr)
{
  int i;
  pInfo *partBuffer = (pInfo*)malloc(sizeof(pInfo) * tr->NumberOfModels);

  for(i = 0; i < tr->NumberOfModels; i++)
    {
      switch(tr->partitionData[i].dataType)
	{
	case DNA_DATA:
	  partBuffer[i].EIGN = (double*)malloc(3 * sizeof(double));
	  partBuffer[i].EV   = (double*)malloc(16 * sizeof(double));
	  partBuffer[i].EI   = (double*)malloc(12 * sizeof(double));	  
	  partBuffer[i].substRates = (double *)malloc(5 * sizeof(double));	  
	  partBuffer[i].frequencies =  (double*)malloc(4 * sizeof(double));	  
	  partBuffer[i].tipVector   = (double *)malloc(64 * sizeof(double));
	  break;
	case AA_DATA:
	  partBuffer[i].EIGN = (double*)malloc(19 * sizeof(double));
	  partBuffer[i].EV   = (double*)malloc(400 * sizeof(double));
	  partBuffer[i].EI   = (double*)malloc(380 * sizeof(double));	  
	  partBuffer[i].tipVector   = (double *)malloc(460 * sizeof(double)); 	  
	  partBuffer[i].substRates = (double *)malloc(190 * sizeof(double));	  
	  partBuffer[i].frequencies = (double*)malloc(20 * sizeof(double));
	  break;
	case BINARY_DATA:
	  partBuffer[i].EIGN = (double*)malloc(1 * sizeof(double));
	  partBuffer[i].EV   = (double*)malloc(4 * sizeof(double));
	  partBuffer[i].EI   = (double*)malloc(2 * sizeof(double));	  
	  partBuffer[i].tipVector   = (double *)malloc(8 * sizeof(double)); 	  
	  partBuffer[i].substRates = (double *)malloc(1 * sizeof(double));	  
	  partBuffer[i].frequencies = (double*)malloc(2 * sizeof(double));
	  break;
	case SECONDARY_DATA:
	  partBuffer[i].EIGN = (double*)malloc(15 * sizeof(double));
	  partBuffer[i].EV   = (double*)malloc(256 * sizeof(double));
	  partBuffer[i].EI   = (double*)malloc(240 * sizeof(double));	  
	  partBuffer[i].tipVector   = (double *)malloc(4096 * sizeof(double)); 	  
	  partBuffer[i].substRates = (double *)malloc(120 * sizeof(double));	  
	  partBuffer[i].frequencies = (double*)malloc(16 * sizeof(double));
	  break;
	case SECONDARY_DATA_6:
	  partBuffer[i].EIGN = (double*)malloc(5 * sizeof(double));
	  partBuffer[i].EV   = (double*)malloc(36 * sizeof(double));
	  partBuffer[i].EI   = (double*)malloc(30 * sizeof(double));	  
	  partBuffer[i].tipVector   = (double *)malloc(384 * sizeof(double)); 	  
	  partBuffer[i].substRates = (double *)malloc(15 * sizeof(double));	  
	  partBuffer[i].frequencies = (double*)malloc(6 * sizeof(double));
	  break; 
	case SECONDARY_DATA_7:
	  partBuffer[i].EIGN = (double*)malloc(6 * sizeof(double));
	  partBuffer[i].EV   = (double*)malloc(49 * sizeof(double));
	  partBuffer[i].EI   = (double*)malloc(42 * sizeof(double));	  
	  partBuffer[i].tipVector   = (double *)malloc(896 * sizeof(double)); 	  
	  partBuffer[i].substRates = (double *)malloc(21 * sizeof(double));	  
	  partBuffer[i].frequencies = (double*)malloc(7 * sizeof(double));
	  break;
	default:
	  assert(0);
	}
    }

  return partBuffer;
      
}

static void freeParams(int numberOfModels, pInfo *partBuffer)
{
  int i;

  for(i = 0; i < numberOfModels; i++)
    {
      free(partBuffer[i].EIGN); 
      free(partBuffer[i].EV);   
      free(partBuffer[i].EI);   
      free(partBuffer[i].substRates);
      free(partBuffer[i].frequencies); 
      free(partBuffer[i].tipVector);  
    }
      
}

static void copyParams(int numberOfModels, pInfo *dst, pInfo *src, tree *tr)
{
  int i;

  assert(src != dst);

  for(i = 0; i < numberOfModels; i++)
    {
      dst[i].dataType = src[i].dataType;

      switch(src[i].dataType)
	{
	case DNA_DATA:
	  memcpy(dst[i].EIGN,        src[i].EIGN,        3 * sizeof(double));
	  memcpy(dst[i].EV,          src[i].EV,          16 * sizeof(double));
	  memcpy(dst[i].EI,          src[i].EI,          12 * sizeof(double));	  
	  memcpy(dst[i].substRates,  src[i].substRates,  5 * sizeof(double));	  
	  memcpy(dst[i].frequencies, src[i].frequencies, 4 * sizeof(double));	  
	  memcpy(dst[i].tipVector,   src[i].tipVector,   64 * sizeof(double));
	  break;
	case AA_DATA:
	  memcpy(dst[i].EIGN,        src[i].EIGN,    19 * sizeof(double));
	  memcpy(dst[i].EV,          src[i].EV,     400 * sizeof(double));
	  memcpy(dst[i].EI,          src[i].EI,     380 * sizeof(double));	  
	  memcpy(dst[i].tipVector,   src[i].tipVector, 460 * sizeof(double)); 	  
	  memcpy(dst[i].substRates,  src[i].substRates, 190 * sizeof(double));	  
	  memcpy(dst[i].frequencies, src[i].frequencies, 20 * sizeof(double));
	  break;
	case BINARY_DATA:
	  memcpy(dst[i].EIGN,        src[i].EIGN, 1 * sizeof(double));
	  memcpy(dst[i].EV,          src[i].EV, 4 * sizeof(double));
	  memcpy(dst[i].EI,          src[i].EI, 2 * sizeof(double));	  
	  memcpy(dst[i].tipVector,   src[i].tipVector, 8 * sizeof(double)); 	  
	  memcpy(dst[i].substRates,  src[i].substRates, 1 * sizeof(double));	  
	  memcpy(dst[i].frequencies, src[i].frequencies, 2 * sizeof(double));
	  break;
	case SECONDARY_DATA:
	  memcpy(dst[i].EIGN,        src[i].EIGN, 15 * sizeof(double));
	  memcpy(dst[i].EV,          src[i].EV, 256 * sizeof(double));
	  memcpy(dst[i].EI,          src[i].EI, 240 * sizeof(double));	  
	  memcpy(dst[i].tipVector,   src[i].tipVector, 4096 * sizeof(double)); 	  
	  memcpy(dst[i].substRates,  src[i].substRates, 120 * sizeof(double));	  
	  memcpy(dst[i].frequencies, src[i].frequencies, 16 * sizeof(double));
	  break;
	case SECONDARY_DATA_6:
	  memcpy(dst[i].EIGN,        src[i].EIGN, 5 * sizeof(double));
	  memcpy(dst[i].EV,          src[i].EV, 36 * sizeof(double));
	  memcpy(dst[i].EI,          src[i].EI, 30 * sizeof(double));	  
	  memcpy(dst[i].tipVector,   src[i].tipVector, 384 * sizeof(double)); 	  
	  memcpy(dst[i].substRates,  src[i].substRates, 15 * sizeof(double));	  
	  memcpy(dst[i].frequencies, src[i].frequencies, 6 * sizeof(double));
	  break;
	case SECONDARY_DATA_7:
	  memcpy(dst[i].EIGN,        src[i].EIGN, 6 * sizeof(double));
	  memcpy(dst[i].EV,          src[i].EV, 49 * sizeof(double));
	  memcpy(dst[i].EI,          src[i].EI, 42 * sizeof(double));	  
	  memcpy(dst[i].tipVector,   src[i].tipVector, 896 * sizeof(double)); 	  
	  memcpy(dst[i].substRates,  src[i].substRates, 21 * sizeof(double));	  
	  memcpy(dst[i].frequencies, src[i].frequencies, 7 * sizeof(double));
	  break;
	default:
	  assert(0);
	}
    }
  
#ifdef _USE_PTHREADS
  masterBarrier(THREAD_COPY_PARAMS, tr);
#endif    

}




#ifndef PARALLEL




void doAllInOne(tree *tr, analdef *adef)
{
  int i, n, sites, bestIndex, model, bootstrapsPerformed;
  double loopTime; 
  int      *originalRateCategories;
  int      *originalInvariant;
  int      slowSearches, fastEvery = 5;
  int treeVectorLength = -1;
  topolRELL_LIST *rl;  
  double bestLH, mlTime, overallTime;  
  long radiusSeed = adef->rapidBoot;
  FILE *infoFile, *f;
  char bestTreeFileName[1024];  
  hashtable *h = (hashtable*)NULL;
  unsigned int **bitVectors = (unsigned int**)NULL;
  boolean bootStopIt = FALSE;
  double pearsonAverage = 0.0;
  pInfo *catParams         = allocParams(tr);
  pInfo *gammaParams = allocParams(tr);
  int vLength;


  n = adef->multipleRuns; 

  if(adef->bootStopping)
    {    
      h = initHashTable(tr->mxtips * 100, tr->mxtips);

      treeVectorLength = adef->multipleRuns;
      
      bitVectors = initBitVector(tr, &vLength);     
    }

  rl = (topolRELL_LIST *)malloc(sizeof(topolRELL_LIST));
  initTL(rl, tr, n);
     
  originalRateCategories = (int*)malloc(tr->cdta->endsite * sizeof(int));      
  originalInvariant      = (int*)malloc(tr->cdta->endsite * sizeof(int));

  sites = tr->cdta->endsite;             

  initModel(tr, tr->rdta, tr->cdta, adef);
 
  for(i = 0; i < n && !bootStopIt; i++)
    {                
      tr->treeID = i;
      tr->checkPointCounter = 0;
         
      loopTime = gettime();      
     
      if(i % 10 == 0)
	{
	  if(i > 0)
	    {	    
	      reductionCleanup(tr, originalRateCategories, originalInvariant);
	    }

	  

	  if(adef->grouping || adef->constraint)
	    {
	      FILE *f = myfopen(tree_file, "r");	

	      assert(adef->restart);
	      partCount = 0;
	      if (! treeReadLenMULT(f, tr, adef))
		exit(-1);
	      /*printf("Constraint read: %d\n", partCount);*/
	      fclose(f);
	    }
	  else
	    makeParsimonyTree(tr, adef);
	  
	  tr->likelihood = unlikely;
	  if(i == 0)
	    {
	      double t;
	          
	      onlyInitrav(tr, tr->start);
	      treeEvaluate(tr, 1);	     	
	     	      
	      t = gettime();    

	      modOpt(tr, adef, FALSE, 5.0);	    
 
	      printf("\nTime for BS model parameter optimization %f\n", gettime() - t);
	      infoFile = myfopen(infoFileName, "a");
	      fprintf(infoFile, "\nTime for BS model parameter optimization %f\n", gettime() - t);
	      fclose(infoFile);
	      
	      memcpy(originalRateCategories, tr->cdta->rateCategory, sizeof(int) * tr->cdta->endsite);
	      memcpy(originalInvariant,      tr->invariant,          sizeof(int) * tr->cdta->endsite);

	      if(adef->bootstrapBranchLengths)
		{
		  if(tr->rateHetModel == CAT)
		    {
		      copyParams(tr->NumberOfModels, catParams, tr->partitionData, tr);		      
		      assert(tr->cdta->endsite == tr->originalCrunchedLength);		 
		      catToGamma(tr, adef);
		      modOpt(tr, adef, TRUE, adef->likelihoodEpsilon);
		      copyParams(tr->NumberOfModels, gammaParams, tr->partitionData, tr);		      
		      gammaToCat(tr);
		      copyParams(tr->NumberOfModels, tr->partitionData, catParams, tr);		      
		    }
		  else
		    {		  
		      assert(tr->cdta->endsite == tr->originalCrunchedLength);		 		     		     		      		     
		    }
		}
	    }	  	  
	}

      computeNextReplicate(tr, &adef->rapidBoot, originalRateCategories, originalInvariant, TRUE); 
      resetBranches(tr);

      evaluateGenericInitrav(tr, tr->start);
    
      treeEvaluate(tr, 1);    	             
     

      computeBOOTRAPID(tr, adef, &radiusSeed);                        	  
      saveTL(rl, tr, i);

      if(adef->bootstrapBranchLengths)
	{
	  double lh = tr->likelihood;
	  int    endsite;
	 
	  if(tr->rateHetModel == CAT)
	    {
	      copyParams(tr->NumberOfModels, tr->partitionData, gammaParams, tr);	      
	      endsite = tr->cdta->endsite;
	      tr->cdta->endsite = tr->originalCrunchedLength;
	      catToGamma(tr, adef);
	      tr->cdta->endsite = endsite;
	      
	      resetBranches(tr);
	      treeEvaluate(tr, 2.0);
	  
	      endsite = tr->cdta->endsite;
	      tr->cdta->endsite = tr->originalCrunchedLength;
	      gammaToCat(tr);
	      tr->cdta->endsite = endsite;	 	    
	
	      copyParams(tr->NumberOfModels, tr->partitionData, catParams, tr);	      
	      tr->likelihood = lh;
	    }
	  else
	    {	     
	      treeEvaluate(tr, 2.0);
	      tr->likelihood = lh;
	    }
	}
      
      printBootstrapResult(tr, adef, TRUE); 

      loopTime = gettime() - loopTime; 
      writeInfoFile(adef, tr, loopTime); 
     
      if(adef->bootStopping)	
	bootStopIt = bootStop(tr, h, i, &pearsonAverage, bitVectors, treeVectorLength, vLength);
    }   

  bootstrapsPerformed = i;

  freeParams(tr->NumberOfModels, catParams);
  free(catParams);

  freeParams(tr->NumberOfModels, gammaParams);
  free(gammaParams);

  if(adef->bootStopping)
    {
      freeBitVectors(bitVectors, 2 * tr->mxtips);
      free(bitVectors);
      freeHashTable(h);
      free(h);
    }

  if(!adef->allInOne)
    {      
      double t = gettime() - masterTime;
      infoFile = myfopen(infoFileName, "a");

      printf("\n\n");
      fprintf(infoFile, "\n\n");

      if(adef->bootStopping)
	{
	  switch(tr->bootStopCriterion)
	    {
	    case FREQUENCY_STOP:
	      printf("Stopped Rapid BS search after %d replicates with FC Bootstopping criterion\n", bootstrapsPerformed);
	      printf("Pearson Average of %d random splits: %f\n",BOOTSTOP_PERMUTATIONS , pearsonAverage);
	      fprintf(infoFile, "Stopped Rapid BS search after %d replicates with Bootstopping criterion\n", 
		      bootstrapsPerformed);	  
	      fprintf(infoFile, "Pearson Average of %d random splits: %f\n", BOOTSTOP_PERMUTATIONS, pearsonAverage);
	      break;
	    case WC_STOP:
	      printf("Stopped Rapid BS search after %d replicates with WC Bootstopping criterion\n", bootstrapsPerformed);
	      printf("WRF Average of %d random splits: %f\n", BOOTSTOP_PERMUTATIONS, pearsonAverage);
	      fprintf(infoFile, "Stopped Rapid BS search after %d replicates with Bootstopping criterion\n", 
		      bootstrapsPerformed);	  
	      fprintf(infoFile, "WRF Average of %d random splits: %f\n", BOOTSTOP_PERMUTATIONS, pearsonAverage);
	      break;
	    default:
	      assert(0);
	    }
	}

      printf("Overall Time for %d Rapid Bootstraps %f\n", bootstrapsPerformed, t);
      fprintf(infoFile, "Overall Time for %d Rapid Bootstraps %f\n", bootstrapsPerformed, t);

      printf("Average Time per Rapid Bootstrap %f\n", (double)(t/((double)bootstrapsPerformed)));
      fprintf(infoFile, "Average Time per Rapid Bootstrap %f\n", (double)(t/((double)bootstrapsPerformed)));

      printf("All %d bootstrapped trees written to: %s\n", bootstrapsPerformed, bootstrapFileName);
      fprintf(infoFile, "All %d bootstrapped trees written to: %s\n", bootstrapsPerformed, bootstrapFileName);           
     
      fclose(infoFile);
      exit(0);
    }
  else
    {
      double t = gettime() - masterTime;
      infoFile = myfopen(infoFileName, "a");

      printf("\n\n");
      fprintf(infoFile, "\n\n");

      if(adef->bootStopping)
	{
	  switch(tr->bootStopCriterion)
	    {
	    case FREQUENCY_STOP:
	      printf("Stopped Rapid BS search after %d replicates with FC Bootstopping criterion\n", bootstrapsPerformed);
	      printf("Pearson Average of %d random splits: %f\n",BOOTSTOP_PERMUTATIONS , pearsonAverage);
	      fprintf(infoFile, "Stopped Rapid BS search after %d replicates with Bootstopping criterion\n", 
		      bootstrapsPerformed);	  
	      fprintf(infoFile, "Pearson Average of %d random splits: %f\n", BOOTSTOP_PERMUTATIONS, pearsonAverage);
	      break;
	    case WC_STOP:
	      printf("Stopped Rapid BS search after %d replicates with WC Bootstopping criterion\n", bootstrapsPerformed);
	      printf("WRF Average of %d random splits: %f\n", BOOTSTOP_PERMUTATIONS, pearsonAverage);
	      fprintf(infoFile, "Stopped Rapid BS search after %d replicates with Bootstopping criterion\n", 
		      bootstrapsPerformed);	  
	      fprintf(infoFile, "WRF Average of %d random splits: %f\n", BOOTSTOP_PERMUTATIONS, pearsonAverage);
	      break;
	    default:
	      assert(0);
	    }
	}

      printf("Overall Time for %d Rapid Bootstraps %f\n", bootstrapsPerformed, t);
      fprintf(infoFile, "Overall Time for %d Rapid Bootstraps %f\n", bootstrapsPerformed, t);

      printf("Average Time per Rapid Bootstrap %f\n", (double)(t/((double)bootstrapsPerformed)));                
      fprintf(infoFile, "Average Time per Rapid Bootstrap %f\n", (double)(t/((double)bootstrapsPerformed)));	           
     
      fclose(infoFile);
    }
  
  /* ML-search */ 

  mlTime = gettime();

  infoFile = myfopen(infoFileName, "a");
  printf("\nStarting ML Search ...\n\n");
  fprintf(infoFile, "\nStarting ML Search ...\n\n");
  fclose(infoFile);

  /***CLEAN UP reduction stuff */  

  reductionCleanup(tr, originalRateCategories, originalInvariant);  

  /****/     	   
  
  restoreTL(rl, tr, 0);

  resetBranches(tr);

  evaluateGenericInitrav(tr, tr->start);   
  
  modOpt(tr, adef, TRUE, adef->likelihoodEpsilon);    

  if(adef->bootStopping)
    {
      if(bootstrapsPerformed <= 100)
	fastEvery = 5;
      else
	fastEvery = bootstrapsPerformed / 20;     
    }
  else    
    fastEvery = 5;    

  for(i = 0; i < bootstrapsPerformed; i++)
    {            
      rl->t[i]->likelihood = unlikely;
    
      if(i % fastEvery == 0)
	{	  	 
	  restoreTL(rl, tr, i); 	 	    	   	
	  
	  resetBranches(tr);	 

	  evaluateGenericInitrav(tr, tr->start);
	  	  
	  treeEvaluate(tr, 1); 		 
	  
	  /*printf("F%d %f -> ", i, tr->likelihood);*/

	  optimizeRAPID(tr, adef);	  			         

	  /*printf(" %f\n", tr->likelihood);*/

	  saveTL(rl, tr, i);   	 
	}    
    }     

  infoFile = myfopen(infoFileName, "a");
  printf("Fast ML optimization finished\n\n");
  fprintf(infoFile, "Fast ML optimization finished\n\n");
  fclose(infoFile);  
  qsort(&(rl->t[0]), bootstrapsPerformed, sizeof(topolRELL*), compareTopolRell);
       
  restoreTL(rl, tr, 0);

  resetBranches(tr);

  evaluateGenericInitrav(tr, tr->start);

  modOpt(tr, adef, TRUE, adef->likelihoodEpsilon);     
  
  slowSearches = bootstrapsPerformed / 5;
  if(bootstrapsPerformed % 5 != 0)
    slowSearches++;

  slowSearches  = MIN(slowSearches, 10); 

  for(i = 0; i < slowSearches; i++)
    {           
      restoreTL(rl, tr, i);     
      rl->t[i]->likelihood = unlikely;  
      
      evaluateGenericInitrav(tr, tr->start);

      treeEvaluate(tr, 1.0);   
      /*printf("S%d %f -> ", i, tr->likelihood);*/
      thoroughOptimization(tr, adef, rl, i); 
      /*printf(" %f\n", tr->likelihood);*/
   }

  infoFile = myfopen(infoFileName, "a");
  printf("Slow ML optimization finished\n\n");
  fprintf(infoFile, "Slow ML optimization finished\n\n");
  fclose(infoFile);

  /*************************************************************************************************************/  
  
 

  if(tr->rateHetModel == CAT) 
    {      
      catToGamma(tr, adef);    
      modOpt(tr, adef, TRUE, adef->likelihoodEpsilon); 
    }

  bestIndex = -1;
  bestLH = unlikely;
    
  for(i = 0; i < slowSearches; i++)
    { 
      restoreTL(rl, tr, i);
      resetBranches(tr);

      evaluateGenericInitrav(tr, tr->start);

      treeEvaluate(tr, 2);

      infoFile = myfopen(infoFileName, "a"); 

      printf("Slow ML Search %d Likelihood: %f\n", i, tr->likelihood);
      fprintf(infoFile, "Slow ML Search %d Likelihood: %f\n", i, tr->likelihood);     
      fclose(infoFile);
      if(tr->likelihood > bestLH)
	{
	  bestLH = tr->likelihood;
	  bestIndex = i;
	}
    }
  
  restoreTL(rl, tr, bestIndex);
  resetBranches(tr);

  evaluateGenericInitrav(tr, tr->start);
 
  treeEvaluate(tr, 2); 
         
  Thorough = 1;
  tr->doCutoff = FALSE;  
  
  treeOptimizeThorough(tr, 1, 10);
  evaluateGenericInitrav(tr, tr->start);
  /*printf("%f \n", tr->likelihood);*/
  modOpt(tr, adef, TRUE, adef->likelihoodEpsilon);
 

  infoFile = myfopen(infoFileName, "a"); 
  printf("\nFinal ML Optimization Likelihood: %f\n", tr->likelihood);
  fprintf(infoFile, "\nFinal ML Optimization Likelihood: %f\n", tr->likelihood);
  
  printf("\nModel Information:\n\n");
  fprintf(infoFile, "\nModel Information:\n\n");

  for(model = 0; model < tr->NumberOfModels; model++)		    		    
    {
      double tl;
      char typeOfData[1024];
      
      switch(tr->partitionData[model].dataType)
	{
	case AA_DATA:
	  strcpy(typeOfData,"AA");
	  break;
	case DNA_DATA:
	  strcpy(typeOfData,"DNA");
	  break;
	case BINARY_DATA:
	  strcpy(typeOfData,"BINARY/MORPHOLOGICAL");
	  break;
	case SECONDARY_DATA:
	  strcpy(typeOfData,"SECONDARY 16 STATE MODEL");
	  break;
	case SECONDARY_DATA_6:
	  strcpy(typeOfData,"SECONDARY 6 STATE MODEL");
	  break;  
	case SECONDARY_DATA_7:
	  strcpy(typeOfData,"SECONDARY 7 STATE MODEL");
	  break;
	default:
	  assert(0);
	}

      fprintf(infoFile, "Model Parameters of Partition %d, Name: %s, Type of Data: %s\n", 
	      model, tr->partitionData[model].partitionName, typeOfData);
      fprintf(infoFile, "alpha: %f\n", tr->partitionData[model].alpha);

      printf("Model Parameters of Partition %d, Name: %s, Type of Data: %s\n", 
	     model, tr->partitionData[model].partitionName, typeOfData);
      printf("alpha: %f\n", tr->partitionData[model].alpha);
      
      if(adef->useInvariant)
	{
	  fprintf(infoFile, "invar: %f\n", tr->partitionData[model].propInvariant);    
	  printf("invar: %f\n", tr->partitionData[model].propInvariant);    
	}
                 
      if(adef->perGeneBranchLengths)
	tl = treeLength(tr, model);
      else
	tl = treeLength(tr, 0);

      fprintf(infoFile, "Tree-Length: %f\n", tl);    
      printf("Tree-Length: %f\n", tl);       

      

      switch(tr->partitionData[model].dataType)
	{
	case AA_DATA:
	  break;
	case DNA_DATA:
	  {
	    int k;
	    char *names[6] = {"a<->c", "a<->g", "a<->t", "c<->g", "c<->t", "g<->t"};	 
	    for(k = 0; k < DNA_RATES; k++)			    
	      {
		fprintf(infoFile, "rate %s: %f\n", names[k], tr->partitionData[model].substRates[k]);			    
		printf(           "rate %s: %f\n", names[k], tr->partitionData[model].substRates[k]);
	      }
	
	    fprintf(infoFile, "rate %s: %f\n", names[5], 1.0);
	    printf("rate %s: %f\n", names[5], 1.0);
	  }      
	  break;
	case SECONDARY_DATA_6:	  
	  /* TODO add printout here */	  
	case SECONDARY_DATA_7:
	case SECONDARY_DATA:
	case BINARY_DATA:	 
	  break;
	default:
	  assert(0);
	}

      fprintf(infoFile, "\n");
      printf("\n");
    }		    		  
  

  fclose(infoFile); 
  
  
  strcpy(bestTreeFileName, workdir); 
  strcat(bestTreeFileName, "RAxML_bestTree.");
  strcat(bestTreeFileName,         run_id);
   
  Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE, TRUE, adef, SUMMARIZE_LH);
  f = myfopen(bestTreeFileName, "w");
  fprintf(f, "%s", tr->tree_string);
  fclose(f);

  if(adef->perGeneBranchLengths)
    printTreePerGene(tr, adef, bestTreeFileName, "w");

  infoFile = myfopen(infoFileName, "a"); 
  overallTime = gettime() - masterTime;
  mlTime    = gettime() - mlTime;
  printf("\nML search took %f secs or %f hours\n", mlTime, mlTime / 3600.0);
  fprintf(infoFile, "\nML search took %f secs or %f hours\n", mlTime, mlTime / 3600.0);
  printf("\nCombined Bootstrap and ML search took %f secs or %f hours\n", overallTime, overallTime / 3600.0);
  fprintf(infoFile, "\nCombined Bootstrap and ML search took %f secs or %f hours\n", overallTime, overallTime / 3600.0);  

  printf("\nDrawing Bootstrap Support Values on best-scoring ML tree ...\n\n");
  fprintf(infoFile, "Drawing Bootstrap Support Values on best-scoring ML tree ...\n\n");  

  fclose(infoFile);
    
  freeTL(rl, tr);   
  free(rl);       
  
  calcBipartitions(tr, adef, bestTreeFileName, bootstrapFileName);  

  overallTime = gettime() - masterTime;
  printf("Program execution info written to %s\n", infoFileName);
  printf("All %d bootstrapped trees written to: %s\n\n", adef->multipleRuns, bootstrapFileName);
  printf("Best-scoring ML tree written to: %s\n\n", bestTreeFileName);
  if(adef->perGeneBranchLengths && tr->NumberOfModels > 1)    
    printf("Per-Partition branch lengths of best-scoring ML tree written to %s.PARTITION.0 to  %s.PARTITION.%d\n\n", bestTreeFileName,  bestTreeFileName, 
	   tr->NumberOfModels - 1);    
  printf("Best-scoring ML tree with support values written to: %s\n\n", bipartitionsFileName);
  printf("Overall execution time for full ML analysis: %f secs or %f hours or %f days\n\n", overallTime, overallTime/3600.0, overallTime/86400.0);
  
  infoFile = myfopen(infoFileName, "a");
  fprintf(infoFile, "All %d bootstrapped trees written to: %s\n\n", adef->multipleRuns, bootstrapFileName);
  fprintf(infoFile, "Best-scoring ML tree written to: %s\n\n", bestTreeFileName);
  if(adef->perGeneBranchLengths && tr->NumberOfModels > 1)    
    fprintf(infoFile, "Per-Partition branch lengths of best-scoring ML tree written to %s.PARTITION.0 to  %s.PARTITION.%d\n\n", bestTreeFileName,  bestTreeFileName, 
	   tr->NumberOfModels - 1);    
  fprintf(infoFile, "Best-scoring ML tree with support values written to: %s\n\n", bipartitionsFileName);
  fprintf(infoFile, "Overall execution time for full ML analysis: %f secs or %f hours or %f days\n\n", overallTime, overallTime/3600.0, overallTime/86400.0);
  fclose(infoFile);
   
  exit(0); 
}


/*******************************************EXPERIMENTAL FUNCTIONS END *****************************************************/





void doBootstrap(tree *tr, analdef *adef, rawdata *rdta, cruncheddata *cdta)
{
  int 
    i, 
    n, 
    treeVectorLength = -1, 
    vLength = -1;
  double loopTime, pearsonAverage;
  hashtable *h = (hashtable*)NULL;
  unsigned int **bitVectors = (unsigned int **)NULL;
  boolean bootStopIt = FALSE; 

  n = adef->multipleRuns; 

  if(adef->bootStopping)
    {    
      h = initHashTable(tr->mxtips * 100, tr->mxtips);
      bitVectors = initBitVector(tr, &vLength);    

      treeVectorLength = adef->multipleRuns;             
    }           

  for(i = 0; i < n && !bootStopIt; i++)
    {    
      loopTime = gettime();
                    
      singleBootstrap(tr, i, adef, rdta, cdta);              
      loopTime = gettime() - loopTime;     
      writeInfoFile(adef, tr, loopTime);  

      if(adef->bootStopping)	
	bootStopIt = bootStop(tr, h, i, &pearsonAverage, bitVectors, treeVectorLength, vLength);
    }      

  if(adef->bootStopping)
    {
      freeBitVectors(bitVectors, 2 * tr->mxtips);
      free(bitVectors);
      freeHashTable(h);
      free(h);
    }
  
   if(adef->bootStopping)
     { 
       FILE *infoFile = myfopen(infoFileName, "a");
       
       if(bootStopIt)
	 {
	   switch(tr->bootStopCriterion)
	     {
	     case FREQUENCY_STOP:
	       printf("Stopped Rapid BS search after %d replicates with FC Bootstopping criterion\n", i);
	       printf("Pearson Average of %d random splits: %f\n",BOOTSTOP_PERMUTATIONS , pearsonAverage);
	       fprintf(infoFile, "Stopped Rapid BS search after %d replicates with Bootstopping criterion\n", i);	  
	       fprintf(infoFile, "Pearson Average of %d random splits: %f\n", BOOTSTOP_PERMUTATIONS, pearsonAverage);
	       break;
	     case WC_STOP:
	       printf("Stopped Rapid BS search after %d replicates with WC Bootstopping criterion\n", i);
	       printf("WRF Average of %d random splits: %f\n", BOOTSTOP_PERMUTATIONS, pearsonAverage);
	       fprintf(infoFile, "Stopped Rapid BS search after %d replicates with Bootstopping criterion\n", i);	  
	       fprintf(infoFile, "WRF Average of %d random splits: %f\n", BOOTSTOP_PERMUTATIONS, pearsonAverage);
	       break;
	     default:
	       assert(0);
	     }
	 }
       else
	 {
	   switch(tr->bootStopCriterion)
	     {
	     case FREQUENCY_STOP:
	       printf("Rapid BS search did not converge after %d replicates with FC Bootstopping criterion\n", i);
	       printf("Pearson Average of %d random splits: %f\n",BOOTSTOP_PERMUTATIONS , pearsonAverage);
	       fprintf(infoFile, "Rapid BS search did not converge after %d replicates with Bootstopping criterion\n", i);	  
	       fprintf(infoFile, "Pearson Average of %d random splits: %f\n", BOOTSTOP_PERMUTATIONS, pearsonAverage);
	       break;
	     case WC_STOP:
	       printf("Rapid BS search did not converge after %d replicates with WC Bootstopping criterion\n", i);
	       printf("WRF Average of %d random splits: %f\n", BOOTSTOP_PERMUTATIONS, pearsonAverage);
	       fprintf(infoFile, "Stopped Rapid BS search after %d replicates with Bootstopping criterion\n", i);	  
	       fprintf(infoFile, "WRF Average of %d random splits: %f\n", BOOTSTOP_PERMUTATIONS, pearsonAverage);
	       break;
	     default:
	       assert(0);
	     }

	 }
       fclose(infoFile);
     }
}

void doInference(tree *tr, analdef *adef, rawdata *rdta, cruncheddata *cdta)
{
  int i, n;
  double loopTime;
  topolRELL_LIST *rl; 
  int 
    best = -1,
    newBest = -1;
  double 
    bestLH = unlikely; 
  FILE *f;
  char bestTreeFileName[1024]; 
  double overallTime;

  n = adef->multipleRuns;
       
  rl = (topolRELL_LIST *)malloc(sizeof(topolRELL_LIST));
  initTL(rl, tr, n);

  for(i = 0; i < n; i++)
    {     
      tr->treeID = i;
      tr->checkPointCounter = 0;
         
      loopTime = gettime();
                                       
      initModel(tr, rdta, cdta, adef);      

      getStartingTree(tr, adef);                      
            
      computeBIGRAPID(tr, adef, TRUE);  

      if(tr->likelihood > bestLH)
	{
	  best = i;
	  bestLH = tr->likelihood;
	}
     
      saveTL(rl, tr, i);

      loopTime = gettime() - loopTime; 
      writeInfoFile(adef, tr, loopTime);
     
    }
    
  f = myfopen(infoFileName, "a");
  printf("\n\nConducting final model optimizations on all %d trees under GAMMA-based models ....\n\n", n);
  fprintf(f, "\n\nConducting final model optimizations on all %d trees under GAMMA-based models ....\n\n", n);
  fclose(f);

  assert(best >= 0);

  if(tr->rateHetModel == GAMMA ||  tr->rateHetModel == GAMMA_I)
    {
      restoreTL(rl, tr, best);
      onlyInitrav(tr, tr->start);
      modOpt(tr, adef, FALSE, adef->likelihoodEpsilon);  
      bestLH = tr->likelihood;
      tr->likelihoods[best] = tr->likelihood;
      saveTL(rl, tr, best);
      tr->treeID = best; 
      printResult(tr, adef, TRUE);
      newBest = best;      
      
      for(i = 0; i < n; i++)
	{	  
	  if(i != best)
	    {
	      restoreTL(rl, tr, i);
	      onlyInitrav(tr, tr->start);
	      treeEvaluate(tr, 1);
	      tr->likelihoods[i] = tr->likelihood;

	      if(tr->likelihood > bestLH)
		{
		  newBest = i;
		  bestLH = tr->likelihood;		  
		  saveTL(rl, tr, i);
		}
	      tr->treeID = i;
	      printResult(tr, adef, TRUE);
	    }
	  f = myfopen(infoFileName, "a");
	  if(n == 1)
	    {
	      printf("Inference[%d] final GAMMA-based Likelihood: %f tree written to file %s\n", i, tr->likelihoods[i], resultFileName);
	      fprintf(f, "Inference[%d] final GAMMA-based Likelihood: %f tree written to file %s\n", i, tr->likelihoods[i], resultFileName);
	    }
	  else
	    {
	      printf("Inference[%d] final GAMMA-based Likelihood: %f tree written to file %s.RUN.%d\n", i, tr->likelihoods[i], resultFileName, i);
	      fprintf(f, "Inference[%d] final GAMMA-based Likelihood: %f tree written to file %s.RUN.%d\n", i, tr->likelihoods[i], resultFileName, i);
	    }
	  fclose(f);	 
	}    
    }
  else
    {     
      catToGamma(tr, adef);

      for(i = 0; i < n; i++)
	rl->t[i]->likelihood = unlikely;
      
      initModel(tr, rdta, cdta, adef);
      
      restoreTL(rl, tr, best);      

      resetBranches(tr);
      onlyInitrav(tr, tr->start);
      modOpt(tr, adef, TRUE, adef->likelihoodEpsilon);      
      tr->likelihoods[best] = tr->likelihood;
      bestLH = tr->likelihood;     
      saveTL(rl, tr, best);
      tr->treeID = best;
      printResult(tr, adef, TRUE);
      newBest = best;

      for(i = 0; i < n; i++)
	{
	  if(i != best)
	    {
	      restoreTL(rl, tr, i);	    
	      resetBranches(tr);
	      onlyInitrav(tr, tr->start);
	      treeEvaluate(tr, 2);
	      tr->likelihoods[i] = tr->likelihood;

	      if(tr->likelihood > bestLH)
		{ 
		  newBest = i;
		  bestLH = tr->likelihood;		
		  saveTL(rl, tr, i);	  
		}
	      tr->treeID = i;
	      printResult(tr, adef, TRUE);
	    } 
	  f = myfopen(infoFileName, "a");
	  if(n == 1)
	    {
	      printf("Inference[%d] final GAMMA-based Likelihood: %f tree written to file %s\n", i, tr->likelihoods[i], resultFileName);
	      fprintf(f, "Inference[%d] final GAMMA-based Likelihood: %f tree written to file %s\n", i, tr->likelihoods[i], resultFileName);
	    }
	  else
	    {
	      printf("Inference[%d] final GAMMA-based Likelihood: %f tree written to file %s.RUN.%d\n", i, tr->likelihoods[i], resultFileName, i);
	      fprintf(f, "Inference[%d] final GAMMA-based Likelihood: %f tree written to file %s.RUN.%d\n", i, tr->likelihoods[i], resultFileName, i);
	    }
	  fclose(f);	  
	}     
    }     
  assert(newBest >= 0);

  restoreTL(rl, tr, newBest);
  onlyInitrav(tr, tr->start);
  f = myfopen(infoFileName, "a");
  printf("\n\nStarting final GAMMA-based thorough Optimization on tree %d likelihood %f .... \n\n", newBest, tr->likelihoods[newBest]);
  fprintf(f, "\n\nStarting final GAMMA-based thorough Optimization on tree %d likelihood %f .... \n\n", newBest, tr->likelihoods[newBest]);
  fclose(f);

  Thorough = 1;
  tr->doCutoff = FALSE; 
  treeOptimizeThorough(tr, 1, 10); 
  evaluateGenericInitrav(tr, tr->start);

  f = myfopen(infoFileName, "a");
  printf("Final GAMMA-based Score of best tree %f\n\n", tr->likelihood); 
  fprintf(f, "Final GAMMA-based Score of best tree %f\n\n", tr->likelihood);  
  fclose(f);

  strcpy(bestTreeFileName, workdir); 
  strcat(bestTreeFileName, "RAxML_bestTree.");
  strcat(bestTreeFileName,         run_id);
   
  Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE, TRUE, adef, SUMMARIZE_LH);
  f = myfopen(bestTreeFileName, "w");
  fprintf(f, "%s", tr->tree_string);
  fclose(f);

  if(adef->perGeneBranchLengths)
    printTreePerGene(tr, adef, bestTreeFileName, "w");
  
  overallTime = gettime() - masterTime;
  printf("Program execution info written to %s\n", infoFileName); 
  printf("Best-scoring ML tree written to: %s\n\n", bestTreeFileName);

  if(adef->perGeneBranchLengths && tr->NumberOfModels > 1)    
    printf("Per-Partition branch lengths of best-scoring ML tree written to %s.PARTITION.0 to  %s.PARTITION.%d\n\n", bestTreeFileName,  bestTreeFileName, 
	   tr->NumberOfModels - 1);     
  printf("Overall execution time: %f secs or %f hours or %f days\n\n", overallTime, overallTime/3600.0, overallTime/86400.0);
  
  f = myfopen(infoFileName, "a");
  fprintf(f, "Program execution info written to %s\n", infoFileName); 
  fprintf(f, "Best-scoring ML tree written to: %s\n\n", bestTreeFileName);
  
  if(adef->perGeneBranchLengths && tr->NumberOfModels > 1)    
    fprintf(f, "Per-Partition branch lengths of best-scoring ML tree written to %s.PARTITION.0 to  %s.PARTITION.%d\n\n", bestTreeFileName,  bestTreeFileName, 
	   tr->NumberOfModels - 1);     
  fprintf(f, "Overall execution time: %f secs or %f hours or %f days\n\n", overallTime, overallTime/3600.0, overallTime/86400.0);
  fclose(f);

  freeTL(rl, tr);   
  free(rl); 
  exit(0);
}

#else



#include <mpi.h>

extern int processID;
extern int numOfWorkers;

static void sendTree(tree *tr, analdef *adef, double t, boolean finalPrint, int tag)
{
  int bufferSize, i, bufCount;
  double *buffer;
  char *tree_ptr;

  bufferSize = tr->treeStringLength + 4 + tr->NumberOfModels + tr->NumberOfModels;

  buffer = (double *)malloc(sizeof(double) * bufferSize);
  
  bufCount = 0;
  
  buffer[bufCount++] = (double) adef->bestTrav;
  buffer[bufCount++] = (double) tr->treeID;
  buffer[bufCount++] = tr->likelihood;
  buffer[bufCount++] = t;

  for(i = 0; i < tr->NumberOfModels; i++)        
    buffer[bufCount++] = tr->partitionData[i].alpha;

  for(i = 0; i < tr->NumberOfModels; i++)        
    buffer[bufCount++] = tr->partitionData[i].propInvariant;
    
    
  if(adef->boot || adef->rapidBoot)
    {
     if(adef->bootstrapBranchLengths)
       Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE, finalPrint, adef, SUMMARIZE_LH);
     else
       Tree2String(tr->tree_string, tr, tr->start->back, FALSE, TRUE, FALSE, FALSE, finalPrint, adef, NO_BRANCHES);
    }
  else
    {
      /*if((adef->model == M_GTRCAT || adef->model == M_PROTCAT) && (adef->useMixedModel == 0))
	Tree2String(tr->tree_string, tr, tr->start->back, FALSE, TRUE, FALSE, FALSE, finalPrint, adef, NO_BRANCHES);
	else*/
      Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE, finalPrint, adef, SUMMARIZE_LH);
    }

  tree_ptr = tr->tree_string;

  while(*tree_ptr != ';')    
    buffer[bufCount++] = (double)*tree_ptr++;        
 
  buffer[bufCount++] = (double)(';');
  buffer[bufCount++] = (double)('\n');

  MPI_Send(buffer, bufferSize, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
  
  free(buffer);
}

static void receiveTree(tree *tr, analdef *adef, int workerID, double *t, int tag)
{
  int bufferSize, i, bufCount;
  double *buffer, *buf_ptr;
  char *tree_ptr, content;
  MPI_Status msgStatus; 

  bufferSize = tr->treeStringLength + 4 + tr->NumberOfModels + tr->NumberOfModels;

  buffer = (double *)malloc(sizeof(double) * bufferSize);

  MPI_Recv(buffer, bufferSize, MPI_DOUBLE, workerID, tag, MPI_COMM_WORLD, &msgStatus);
  
  bufCount = 0;
  
  adef->bestTrav = (int)buffer[bufCount++]; 
  tr->treeID     = (int) buffer[bufCount++];
  tr->likelihood = buffer[bufCount++];
  *t = buffer[bufCount++];

  tr->likelihoods[tr->treeID] = tr->likelihood;

  for(i = 0; i < tr->NumberOfModels; i++)    
    tr->partitionData[i].alpha = buffer[bufCount++];

  for(i = 0; i < tr->NumberOfModels; i++)    
    tr->partitionData[i].propInvariant = buffer[bufCount++];

  buf_ptr = &buffer[bufCount];
  tree_ptr = tr->tree_string;

  while((content = (char)(buffer[bufCount++])) != ';')
    {      
      *tree_ptr++ = content;
    }
  
  *tree_ptr++ = ';';
  *tree_ptr++ = '\n';
#ifdef DEBUG
  printf("Received tree %s\n", tr->tree_string);
#endif 
  free(buffer);
}




void doBootstrap(tree *tr, analdef *adef, rawdata *rdta, cruncheddata *cdta)
{
  int i, n, dummy;
  double loopTime;
  MPI_Status msgStatus; 

  n = adef->multipleRuns;
          
  if(processID == 0)
    {
      int jobsSent = 0;
      int jobsReceived = n;

      while(jobsReceived > 0)
	{
	  MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &msgStatus);
	  switch(msgStatus.MPI_TAG)
	    {
	    case JOB_REQUEST:
#ifdef DEBUG
	      printf("Master receiving work request from worker %d\n",  msgStatus.MPI_SOURCE);
#endif	      
	      MPI_Recv(&dummy, 1, MPI_INT, msgStatus.MPI_SOURCE, JOB_REQUEST, MPI_COMM_WORLD, &msgStatus);
	       if(jobsSent < n)
		 {
		   MPI_Send(&jobsSent, 1, MPI_INT, msgStatus.MPI_SOURCE, COMPUTE_TREE, MPI_COMM_WORLD);
#ifdef DEBUG
		   printf("Master sending job %d to worker %d\n",  jobsSent, msgStatus.MPI_SOURCE);
#endif
		   jobsSent++;
		 }
	       break;
	    case TREE:
#ifdef DEBUG
	      printf("--------> Master receiving tree from worker %d\n",  msgStatus.MPI_SOURCE);	
#endif
	      receiveTree(tr, adef, msgStatus.MPI_SOURCE, &loopTime, TREE);	     	   	      
	      printBootstrapResult(tr, adef, TRUE);
	      printf("Bootstrap[%d] completed\n", tr->treeID);	 
	      writeInfoFile(adef, tr, loopTime);
	      jobsReceived--;
	      if(jobsSent < n)
		{
		  MPI_Send(&jobsSent, 1, MPI_INT, msgStatus.MPI_SOURCE, COMPUTE_TREE, MPI_COMM_WORLD);
#ifdef DEBUG
		  printf("Master sending job %d to worker %d\n",  jobsSent, msgStatus.MPI_SOURCE);
#endif
		  jobsSent++;
		}
	      break;
	    }
	}
      
       for(i = 1; i < numOfWorkers; i++)
	{
	  MPI_Send(&dummy, 1, MPI_INT, i, FINALIZE, MPI_COMM_WORLD);
#ifdef DEBUG
	  printf("Master sending FINALIZE to worker %d\n",  i);
#endif
	}
       return;
    }
  else
    {
      int treeCounter = 0;

      MPI_Send(&dummy, 1, MPI_INT, 0, JOB_REQUEST, MPI_COMM_WORLD);
#ifdef DEBUG
      printf("Worker %d sending job request to master\n",  processID);
#endif      
       while(1)
	{	
	  MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &msgStatus); 
	  	 
	  switch(msgStatus.MPI_TAG)
	    {
	    case COMPUTE_TREE: 
	      MPI_Recv(&dummy, 1, MPI_INT, 0, COMPUTE_TREE, MPI_COMM_WORLD, &msgStatus);	      
#ifdef DEBUG
	      printf("Worker %d receiving job %d from master\n",  processID, dummy);
#endif	
	      loopTime = masterTime = gettime();
	      
	      if(adef->multiBoot < 2)
		singleBootstrap(tr, dummy, adef, rdta, cdta);     
	      else
		multipleBootstrap(tr, dummy, adef, rdta, cdta);

	      treeCounter++;
	      loopTime = gettime() - loopTime;
	      sendTree(tr, adef, loopTime, TRUE, TREE);
	      break;
	    case FINALIZE:
	      MPI_Recv(&dummy, 1, MPI_INT, 0, FINALIZE, MPI_COMM_WORLD, &msgStatus);
#ifdef DEBUG
	      printf("Worker %d receiving FINALIZE %d\n",  processID);
#endif
	      return;
	    }
	}
    }  
}

void doInference(tree *tr, analdef *adef, rawdata *rdta, cruncheddata *cdta)
{
  int i, n, dummy;
  double loopTime;
  MPI_Status msgStatus; 

  n = adef->multipleRuns;
          
  if(processID == 0)
    {
      int jobsSent = 0;
      int jobsReceived = n;

      while(jobsReceived > 0)
	{
	  MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &msgStatus);
	  switch(msgStatus.MPI_TAG)
	    {
	    case JOB_REQUEST:
#ifdef DEBUG
	      printf("Master receiving work request from worker %d\n",  msgStatus.MPI_SOURCE);
#endif	      
	      MPI_Recv(&dummy, 1, MPI_INT, msgStatus.MPI_SOURCE, JOB_REQUEST, MPI_COMM_WORLD, &msgStatus);
	      if(jobsSent < n)
		{
		  MPI_Send(&jobsSent, 1, MPI_INT, msgStatus.MPI_SOURCE, COMPUTE_TREE, MPI_COMM_WORLD);
#ifdef DEBUG
		  printf("Master snding job %d to worker %d\n",  jobsSent, msgStatus.MPI_SOURCE);
#endif
		  jobsSent++;
		}
	       break;
	    case TREE:
#ifdef DEBUG
	      printf("--------> Master receiving tree from worker %d\n",  msgStatus.MPI_SOURCE);	
#endif
	      receiveTree(tr, adef, msgStatus.MPI_SOURCE, &loopTime, TREE);	     	   	      	      
	      printf("Inference[%d] completed\n", tr->treeID);	 
	      writeInfoFile(adef, tr, loopTime);
	      jobsReceived--;
	      if(jobsSent < n)
		{
		  MPI_Send(&jobsSent, 1, MPI_INT, msgStatus.MPI_SOURCE, COMPUTE_TREE, MPI_COMM_WORLD);
#ifdef DEBUG
		  printf("Master sending job %d to worker %d\n",  jobsSent, msgStatus.MPI_SOURCE);
#endif
		  jobsSent++;
		}
	      break;
	    }
	}
      
       for(i = 1; i < numOfWorkers; i++)
	{
	  MPI_Send(&dummy, 1, MPI_INT, i, FINALIZE, MPI_COMM_WORLD);
#ifdef DEBUG
	  printf("Master sending FINALIZE to worker %d\n",  i);
#endif
	}
       return;
    }
  else
    {
      int treeCounter = 0;

      MPI_Send(&dummy, 1, MPI_INT, 0, JOB_REQUEST, MPI_COMM_WORLD);
#ifdef DEBUG
      printf("Worker %d sending job request to master\n",  processID);
#endif      
       while(1)
	{	
	  MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &msgStatus); 
	  	 
	  switch(msgStatus.MPI_TAG)
	    {
	    case COMPUTE_TREE: 
	      MPI_Recv(&dummy, 1, MPI_INT, 0, COMPUTE_TREE, MPI_COMM_WORLD, &msgStatus);	      
#ifdef DEBUG
	      printf("Worker %d receiving job %d from master\n",  processID, dummy);
#endif
	      loopTime =  masterTime = gettime();

	      tr->treeID = dummy;
	      tr->checkPointCounter = 0;
	                    
	      
	      initModel(tr, rdta, cdta, adef); 

	      treeCounter++;

	      getStartingTree(tr, adef);  
     
	      computeBIGRAPID(tr, adef, TRUE);                     

	      if(tr->rateHetModel == GAMMA || tr->rateHetModel == GAMMA_I)
		{
		  modOpt(tr, adef, TRUE, adef->likelihoodEpsilon);			
		  printLog(tr, adef, TRUE);
		  printResult(tr, adef, TRUE);
		  loopTime = gettime() - loopTime;	 		 		
		}
	      else
		{	  
		  if(adef->useMixedModel)
		    {
		      tr->likelihood = unlikely;

		      catToGamma(tr, adef);

		      initModel(tr, rdta, cdta, adef);	  	  	  
		      modOpt(tr, adef, TRUE, adef->likelihoodEpsilon);	
		      printLog(tr, adef, TRUE);
		      printResult(tr, adef, TRUE);
		      loopTime = gettime() - loopTime;			    

		      gammaToCat(tr);		      
		    }
		  else
		    {
		      loopTime = gettime() - loopTime;       		      		      
		    }	
		}
	     
	      sendTree(tr, adef, loopTime, TRUE, TREE);
	      break;
	    case FINALIZE:
	      MPI_Recv(&dummy, 1, MPI_INT, 0, FINALIZE, MPI_COMM_WORLD, &msgStatus);
#ifdef DEBUG
	      printf("Worker %d receiving FINALIZE %d\n",  processID);
#endif
	      return;
	    }
	}
    }  
}


static void allInOneMaster(tree *tr, analdef *adef)
{  
  MPI_Status msgStatus;
  double loopTime; 
  int workers = numOfWorkers - 1;
  int i, width, dummy, whoHasBestTree = -1;
  int bsCount, bsTotal;
  int mlCount, mlTotal;
  int finished = FALSE;  
  FILE *infoFile;
  double bestLikelihood = unlikely;  

  if(adef->multipleRuns % workers == 0)
    width = adef->multipleRuns / workers;
  else
    width = (adef->multipleRuns / workers) + 1;    		      

  bsCount = 0;
  mlCount = 0;
  bsTotal = width * workers;
  mlTotal = workers;

  free(tr->likelihoods);
  tr->likelihoods = (double *)malloc((width * workers + workers) * sizeof(double));

  if(tr->rateHetModel == GAMMA || tr->rateHetModel == GAMMA_I)
    {     
      printf("\nSwitching from GAMMA to CAT for rapid Bootstrap, final ML search will be conducted under the %s model you specified\n", 
	     (adef->useInvariant)?"GAMMA+P-Invar":"GAMMA");
      infoFile = myfopen(infoFileName, "a");
      fprintf(infoFile, 
	      "\nSwitching from GAMMA to CAT for rapid Bootstrap, final ML search will be conducted under the %s model you specified\n", 
	      (adef->useInvariant)?"GAMMA+P-Invar":"GAMMA");
      fclose(infoFile);           
    }

  while(! finished)
    {
      MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &msgStatus);
      switch(msgStatus.MPI_TAG)
	{
	case BS_TREE:
	  receiveTree(tr, adef, msgStatus.MPI_SOURCE, &loopTime, BS_TREE);
#ifdef DEBUG	 
	  printf("Received BS TREE %d %f\n", bsCount, tr->likelihood);
#endif
	  printBootstrapResult(tr, adef, TRUE);
	  writeInfoFile(adef, tr, loopTime);
	  bsCount++;
	  break;
	case ML_TREE: 
	  receiveTree(tr, adef, msgStatus.MPI_SOURCE, &loopTime, ML_TREE);
	  if(tr->likelihood > bestLikelihood)
	    {
	      bestLikelihood = tr->likelihood;
	      whoHasBestTree = msgStatus.MPI_SOURCE;
	    }
#ifdef DEBUG
	  printf("Received ML TREE %d %f ID %d\n", mlCount, tr->likelihood, tr->treeID);
#endif
	  mlCount++;
	  break;
	default:
	  assert(0);
	}
      if(adef->allInOne)	
	finished = (bsCount == bsTotal && mlCount == mlTotal);
      else
	finished = (bsCount == bsTotal);

      if(adef->allInOne && bsCount == bsTotal && mlCount == 0)
	{
	  double t = gettime() - masterTime;
	  infoFile = myfopen(infoFileName, "a");

	  printf("\n\n");
	  fprintf(infoFile, "\n\n");

	  if(adef->bootStopping)
	    {
	      assert(0);	     
	    }

	  printf("Overall Time for %d Rapid Bootstraps %f\n", bsCount, t);
	  fprintf(infoFile, "Overall Time for %d Rapid Bootstraps %f\n", bsCount, t);

	  printf("Average Time per Rapid Bootstrap %f\n", (double)(t/((double)bsCount)));                
	  fprintf(infoFile, "Average Time per Rapid Bootstrap %f\n", (double)(t/((double)bsCount)));	           
     
	  fclose(infoFile);    
	}

      if(finished)
	{
	  if(adef->allInOne)
	    {
	      double overallTime;
	      char bestTreeFileName[1024];

	      strcpy(bestTreeFileName, workdir);
	      strcat(bestTreeFileName, "RAxML_bestTree.");
	      strcat(bestTreeFileName,         run_id);

	      assert(whoHasBestTree > 0 && whoHasBestTree < numOfWorkers);
#ifdef DEBUG
	      printf("worker %d xas mpest tri with %f \n", whoHasBestTree, bestLikelihood);
#endif	      
	      MPI_Send(&dummy, 1, MPI_INT, whoHasBestTree, PRINT_TREE, MPI_COMM_WORLD);
	      MPI_Recv(&dummy, 1, MPI_INT, whoHasBestTree, I_PRINTED_IT,  MPI_COMM_WORLD, &msgStatus);

	      overallTime = gettime() - masterTime;
	      printf("Program execution info written to %s\n", infoFileName);
	      printf("All %d bootstrapped trees written to: %s\n\n", bsCount, bootstrapFileName);
	      printf("Best-scoring ML tree written to: %s\n\n", bestTreeFileName);
	      if(adef->perGeneBranchLengths && tr->NumberOfModels > 1)    
		printf("Per-Partition branch lengths of best-scoring ML tree written to %s.PARTITION.0 to  %s.PARTITION.%d\n\n", 
		       bestTreeFileName,  bestTreeFileName, tr->NumberOfModels - 1);    
	      printf("Best-scoring ML tree with support values written to: %s\n\n", bipartitionsFileName);
	      printf("Overall execution time for full ML analysis: %f secs or %f hours or %f days\n\n", 
		     overallTime, overallTime/3600.0, overallTime/86400.0);
  
	      infoFile = myfopen(infoFileName, "a");
	      fprintf(infoFile, "All %d bootstrapped trees written to: %s\n\n", bsCount, bootstrapFileName);
	      fprintf(infoFile, "Best-scoring ML tree written to: %s\n\n", bestTreeFileName);
	      if(adef->perGeneBranchLengths && tr->NumberOfModels > 1)    
		fprintf(infoFile, "Per-Partition branch lengths of best-scoring ML tree written to %s.PARTITION.0 to  %s.PARTITION.%d\n\n"
			, bestTreeFileName,  bestTreeFileName, tr->NumberOfModels - 1);    
	      fprintf(infoFile, "Best-scoring ML tree with support values written to: %s\n\n", 
		      bipartitionsFileName);
	      fprintf(infoFile, "Overall execution time for full ML analysis: %f secs or %f hours or %f days\n\n", 
		      overallTime, overallTime/3600.0, overallTime/86400.0);
	      fclose(infoFile);   
	    }
	  else
	    {
	      double t = gettime() - masterTime;
	      infoFile = myfopen(infoFileName, "a");

	      printf("\n\n");
	      fprintf(infoFile, "\n\n");

	      if(adef->bootStopping)
		{
		  assert(0);
		}

	      printf("Overall Time for %d Rapid Bootstraps %f\n", bsCount, t);
	      fprintf(infoFile, "Overall Time for %d Rapid Bootstraps %f\n", bsCount, t);

	      printf("Average Time per Rapid Bootstrap %f\n", (double)(t/((double)bsCount)));
	      fprintf(infoFile, "Average Time per Rapid Bootstrap %f\n", (double)(t/((double)bsCount)));

	      printf("All %d bootstrapped trees written to: %s\n", bsCount, bootstrapFileName);
	      fprintf(infoFile, "All %d bootstrapped trees written to: %s\n", bsCount, bootstrapFileName);           
     
	      fclose(infoFile);	      
	    }

	  for(i = 1; i < numOfWorkers; i++)	    
	    MPI_Send(&dummy, 1, MPI_INT, i, FINALIZE, MPI_COMM_WORLD);	 
	}
    }
  
  MPI_Finalize();
  exit(0);
}

static void allInOneWorker(tree *tr, analdef *adef)
{
  int dummy, NumberOfLocalTrees;
  MPI_Status msgStatus; 

  if(adef->multipleRuns % (numOfWorkers - 1)  == 0)
    NumberOfLocalTrees = adef->multipleRuns / (numOfWorkers - 1);
  else
    NumberOfLocalTrees = (adef->multipleRuns / (numOfWorkers - 1)) + 1;

#ifdef DEBUG
  printf("Worker %d %d\n", processID, NumberOfLocalTrees);
#endif
  /* re-initialize adfe->rapidBoot, otherwise the workers will be doing the exact same replicates */

  adef->rapidBoot = (long)gettimeSrand();

  /* the one below is kind of an ugly fix, but who cares */

  tr->treeID = NumberOfLocalTrees * (processID - 1);

  	  
  {
    int i, n, sites, bestIndex, model, bootstrapsPerformed;
    double loopTime = 0.0; 
    int      *originalRateCategories;
    int      *originalInvariant;
    int      slowSearches, fastEvery = 5;
    topolRELL_LIST *rl;  
    double bestLH, mlTime;  
    long radiusSeed = adef->rapidBoot;
    FILE *infoFile, *f;
    char bestTreeFileName[1024];  
    modelParams *catParams   = (modelParams *)malloc(sizeof(modelParams));
    modelParams *gammaParams = (modelParams *)malloc(sizeof(modelParams));

    allocParams(catParams,   tr);
    allocParams(gammaParams, tr);
   
    n = NumberOfLocalTrees;   

    rl = (topolRELL_LIST *)malloc(sizeof(topolRELL_LIST));
    initTL(rl, tr, n);
     
    originalRateCategories = (int*)malloc(tr->cdta->endsite * sizeof(int));
    originalInvariant      = (int*)malloc(tr->cdta->endsite * sizeof(int));

    sites = tr->cdta->endsite;             

    if(tr->rateHetModel == GAMMA || tr->rateHetModel == GAMMA_I)
      {       
	tr->rateHetModel = CAT;		
	
	initModel(tr, tr->rdta, tr->cdta, adef);                
      }
 
    for(i = 0; i < n; i++)
      {                
	tr->checkPointCounter = 0;
	
	loopTime = gettime();      
	
	if(i % 10 == 0)
	  {
	    if(i > 0)
	      {	       
		reductionCleanup(tr, originalRateCategories, originalInvariant);
	      }
	    
	    makeParsimonyTree(tr, adef);
	    
	    tr->likelihood = unlikely;
	    if(i == 0)
	      {
		double t;
		
		onlyInitrav(tr, tr->start);
		
		treeEvaluate(tr, 1);	     	     
		
		t = gettime();    
		
		modOpt(tr, adef, FALSE, 5.0);	     		
		memcpy(originalRateCategories, tr->cdta->rateCategory, sizeof(int) * tr->cdta->endsite);
		memcpy(originalInvariant,      tr->invariant,          sizeof(int) * tr->cdta->endsite);

		if(adef->bootstrapBranchLengths)
		  {
		    storeParams(catParams, tr);
		    assert(tr->cdta->endsite == tr->originalCrunchedLength);
		    catToGamma(tr, adef);
		    modOpt(tr, adef, TRUE, adef->likelihoodEpsilon);
		    storeParams(gammaParams, tr);
		    gammaToCat(tr);
		    loadParams(catParams, tr);		  
		  }
	      }	  	  
	  }

	computeNextReplicate(tr, &adef->rapidBoot, originalRateCategories, originalInvariant, TRUE);       
	resetBranches(tr);
	
	evaluateGenericInitrav(tr, tr->start);    
	
	treeEvaluate(tr, 1);    	             
	
	computeBOOTRAPID(tr, adef, &radiusSeed);                        	  
	saveTL(rl, tr, i);
      
	if(adef->bootstrapBranchLengths)
	  {
	    double lh = tr->likelihood;
	    int    endsite;
	    
	    loadParams(gammaParams, tr);
	    
	    
	    endsite = tr->cdta->endsite;
	    tr->cdta->endsite = tr->originalCrunchedLength;
	    catToGamma(tr, adef);
	    tr->cdta->endsite = endsite;
	    
	    resetBranches(tr);
	    treeEvaluate(tr, 2.0);
	    
	    endsite = tr->cdta->endsite;
	    tr->cdta->endsite = tr->originalCrunchedLength;
	    gammaToCat(tr);
	    tr->cdta->endsite = endsite;	 	    
	    
	    loadParams(catParams, tr);
	    	    
	    tr->likelihood = lh;
	  }
      

	loopTime = gettime() - loopTime;
	sendTree(tr, adef, loopTime, TRUE, BS_TREE);	
     
	if(adef->bootStopping)
	  {
	    assert(0);
	    /*bootStopIt = bootStop(tr, b, i, &pearsonAverage);*/
	  }
	tr->treeID = tr->treeID + 1;
      }  
 
    bootstrapsPerformed = i;
        
    if(!adef->allInOne)
      {       	
	MPI_Recv(&dummy, 1, MPI_INT, 0, FINALIZE, MPI_COMM_WORLD, &msgStatus);
	MPI_Finalize();	 
	exit(0);
      }
    
    freeParams(catParams);
    free(catParams);

    freeParams(gammaParams);
    free(gammaParams);

    tr->treeID = NumberOfLocalTrees * (numOfWorkers - 1) + (processID - 1);
   
    mlTime = gettime();
  
#ifdef DEBUG
    printf("\nWorker %d Starting ML Search ...\n\n", processID);   
#endif

    reductionCleanup(tr, originalRateCategories, originalInvariant);    
 
    catToGamma(tr, adef);     	   
    
    restoreTL(rl, tr, 0);

    resetBranches(tr);
    
    evaluateGenericInitrav(tr, tr->start);   
    
    modOpt(tr, adef, TRUE, adef->likelihoodEpsilon);
    
    evaluateGenericInitrav(tr, tr->start);
    
    categorizeGeneric(tr, tr->start);        
    
    gammaToCat(tr);      
         
    fastEvery = 5;    

    for(i = 0; i < bootstrapsPerformed; i++)
      {            
	rl->t[i]->likelihood = unlikely;
	
	if(i % fastEvery == 0)
	  {	 
	    restoreTL(rl, tr, i); 	 	    	   
	    
	    resetBranches(tr);

	    evaluateGenericInitrav(tr, tr->start);
	    
	    treeEvaluate(tr, 1); 		 
	    
	    optimizeRAPID(tr, adef);	  		
	    saveTL(rl, tr, i);      
	  }    
      }     

#ifdef DEBUG  
    printf("Worker %d Fast ML optimization finished\n\n", processID);
#endif
   
    qsort(&(rl->t[0]), bootstrapsPerformed, sizeof(topolRELL*), compareTopolRell);
     
    catToGamma(tr, adef);  
    
    restoreTL(rl, tr, 0);

    resetBranches(tr);
    
    evaluateGenericInitrav(tr, tr->start); 
    
    modOpt(tr, adef, TRUE, adef->likelihoodEpsilon);
    
    evaluateGenericInitrav(tr, tr->start);
    
    categorizeGeneric(tr, tr->start);   
    
    gammaToCat(tr);    
  
    slowSearches = bootstrapsPerformed / 5;
    if(bootstrapsPerformed % 5 != 0)
      slowSearches++;

    slowSearches  = MIN(slowSearches, 10); 

    for(i = 0; i < slowSearches; i++)
      {           
	restoreTL(rl, tr, i);     
	rl->t[i]->likelihood = unlikely;
	
	evaluateGenericInitrav(tr, tr->start);
	
	treeEvaluate(tr, 1.0);   
	thoroughOptimization(tr, adef, rl, i);             
      }
 
#ifdef DEBUG
    printf("Worker %d Slow ML optimization finished\n\n", processID);
#endif     
    catToGamma(tr, adef);    
    
    bestIndex = -1;
    bestLH = unlikely;
    
    for(i = 0; i < slowSearches; i++)
      {      
	restoreTL(rl, tr, i);
	resetBranches(tr);

	evaluateGenericInitrav(tr, tr->start);
	
	treeEvaluate(tr, 2);
#ifdef DEBUG	
	printf("Worker %d Slow ML Search %d Likelihood: %f\n", processID, i, tr->likelihood);
#endif	
	if(tr->likelihood > bestLH)
	  {
	    bestLH = tr->likelihood;
	    bestIndex = i;
	  }
      }
    
    restoreTL(rl, tr, bestIndex);
    resetBranches(tr);

    evaluateGenericInitrav(tr, tr->start);
    
    treeEvaluate(tr, 2); 
    
    Thorough = 1;
    tr->doCutoff = FALSE;  
    
    treeOptimizeThorough(tr, 1, 10);
    modOpt(tr, adef, TRUE, adef->likelihoodEpsilon);

    sendTree(tr, adef, loopTime, TRUE, ML_TREE);
   
    while(1)
      {
	MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &msgStatus);
	switch(msgStatus.MPI_TAG)
	  {
	  case FINALIZE:	
	    MPI_Recv(&dummy, 1, MPI_INT, 0, FINALIZE, MPI_COMM_WORLD, &msgStatus);
	    MPI_Finalize();  
	    exit(0);
	    break;
	  case PRINT_TREE:
	    MPI_Recv(&dummy, 1, MPI_INT, 0, PRINT_TREE, MPI_COMM_WORLD, &msgStatus);
#ifdef DEBUG
	    printf("Ox malaka eimai o %d kai prepi na printaro to dentro %f re pousti \n", processID, tr->likelihood);
#endif
	    infoFile = myfopen(infoFileName, "a"); 
	    
	    printf("\nFinal ML Optimization Likelihood: %f\n", tr->likelihood);
	    fprintf(infoFile, "\nFinal ML Optimization Likelihood: %f\n", tr->likelihood);
	    
	    printf("\nModel Information:\n\n");
	    fprintf(infoFile, "\nModel Information:\n\n");
	    
	    for(model = 0; model < tr->NumberOfModels; model++)		    		    
	      {
		double tl;
		char typeOfData[1024];
		
		switch(tr->partitionData[model].dataType)
		  {
		  case AA_DATA:
		    strcpy(typeOfData,"AA");
		    break;
		  case DNA_DATA:
		    strcpy(typeOfData,"DNA");
		    break;
		  default:
		    assert(0);
		  }
		
		fprintf(infoFile, "Model Parameters of Partition %d, Name: %s, Type of Data: %s\n", 
			model, tr->partitionData[model].partitionName, typeOfData);
		fprintf(infoFile, "alpha: %f\n", tr->partitionData[model].alpha);
		
		printf("Model Parameters of Partition %d, Name: %s, Type of Data: %s\n", 
		       model, tr->partitionData[model].partitionName, typeOfData);
		printf("alpha: %f\n", tr->partitionData[model].alpha);
		
		if(adef->useInvariant)
		  {
		    fprintf(infoFile, "invar: %f\n", tr->partitionData[model].propInvariant);    
		    printf("invar: %f\n", tr->partitionData[model].propInvariant);    
		  }
		
		if(adef->perGeneBranchLengths)
		  tl = treeLength(tr, model);
		else
		  tl = treeLength(tr, 0);
		
		fprintf(infoFile, "Tree-Length: %f\n", tl);    
		printf("Tree-Length: %f\n", tl);       
		
		switch(tr->partitionData[model].dataType)
		  {
		  case AA_DATA:
		    break;
		  case DNA_DATA:
		    {
		      int k;
		      char *names[6] = {"a<->c", "a<->g", "a<->t", "c<->g", "c<->t", "g<->t"};	 
		      for(k = 0; k < DNA_RATES; k++)			    
			{
			  fprintf(infoFile, "rate %s: %f\n", names[k], tr->initialRates_DNA[model * DNA_RATES + k]);			    
			  printf("rate %s: %f\n", names[k], tr->initialRates_DNA[model * DNA_RATES + k]);
			}
		      
		      fprintf(infoFile, "rate %s: %f\n", names[5], 1.0);
		      printf("rate %s: %f\n", names[5], 1.0);
		    }      
		    break;
		  default:
		    assert(0);
		  }
		
		fprintf(infoFile, "\n");
		printf("\n");
	      }		    		  
	    
	    fclose(infoFile); 
	    
	    strcpy(bestTreeFileName, workdir); 
	    strcat(bestTreeFileName, "RAxML_bestTree.");
	    strcat(bestTreeFileName,         run_id);
	    
	    Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE, TRUE, adef, SUMMARIZE_LH);
	    f = myfopen(bestTreeFileName, "w");
	    fprintf(f, "%s", tr->tree_string);
	    fclose(f);

	    if(adef->perGeneBranchLengths)
	      printTreePerGene(tr, adef, bestTreeFileName, "w");
	    
	    infoFile = myfopen(infoFileName, "a"); 
	    
	    
	    printf("\nDrawing Bootstrap Support Values on best-scoring ML tree ...\n\n");
	    fprintf(infoFile, "Drawing Bootstrap Support Values on best-scoring ML tree ...\n\n");  
	    
	    fclose(infoFile);
	    
	    freeTL(rl, tr);   
	    free(rl);       
	    
	    calcBipartitions(tr, adef, bestTreeFileName, bootstrapFileName);  
	     
	    

	    MPI_Send(&dummy, 1, MPI_INT, 0, I_PRINTED_IT, MPI_COMM_WORLD);
	    break;
	  default:
	    assert(0);
	  }
      }
  }
}


void doAllInOne(tree *tr, analdef *adef)
{
  if(processID == 0)
    allInOneMaster(tr, adef);
  else
    allInOneWorker(tr, adef);
}


#endif
