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
#include <sys/times.h>
#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h> 
#endif

#include <math.h>
#include <time.h> 
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>

#include "axml.h"

extern int Thorough;
extern infoList iList;
extern char inverseMeaningBINARY[4];
extern char inverseMeaningDNA[16];
extern char seq_file[1024];
extern char resultFileName[1024];
extern char tree_file[1024];
extern FILE *INFILE;
extern char inverseMeaningPROT[23];
extern char  workdir[1024];
extern char run_id[128];
double masterTime;

typedef struct 
{
  double *likelihoods;
  branchInfo **branch;
  int count;
  int minLH;
} 
  bestBranches;

typedef struct
{
  int length;
  unsigned char *seq;
} 
  unalignedSequences;


static double getBranch(tree *tr, double *b, double *bb)
{
  double z = 0.0;

  if(!tr->multiBranch)
    {
      assert(tr->fracchange != -1.0);
      assert(b[0] == bb[0]);
      z = b[0];
      if (z < zmin) 
	z = zmin;      	 
      if(z > zmax)
	z = zmax;
      z = -log(z) * tr->fracchange;
      return z;	
    }
  else
    {       
      int i;
      double x = 0;
      
      for(i = 0; i < tr->numBranches; i++)
	{
	  assert(b[i] == bb[i]);
	  assert(tr->partitionContributions[i] != -1.0);
	  assert(tr->fracchanges[i] != -1.0);
	  x = b[i];
	  if (x < zmin) 
	    x = zmin;      	 
	  if(x > zmax)
	    x = zmax;
	  x = -log(x) * tr->fracchanges[i];
	  z += x * tr->partitionContributions[i];
	}	
      
      return z;
    } 

}


static char *Tree2StringClassifyRec(char *treestr, tree *tr, nodeptr p, int numberOfTipsForInsertion, int *countBranches, 
				    int *inserts, boolean subtreeSummary, boolean *foundTheSubtree, int CUTOFF, boolean originalTree)
{        
  branchInfo *bInf = p->bInf;
  int        i, countQuery = 0;   

  *countBranches = *countBranches + 1;

  if(!subtreeSummary && !originalTree)
    {
      for(i = 0; i < numberOfTipsForInsertion; i++)
	if(bInf->countThem[i] > 0)
	  countQuery++;  
      
      if(countQuery > 0)
	{
	  *treestr++ = '(';   
	  for(i = 0; i < numberOfTipsForInsertion; i++)
	    {
	      if(bInf->countThem[i] > 0)
		{	      
		  sprintf(treestr,"QUERY___%s___%d", tr->nameList[inserts[i]], bInf->countThem[i]);
		  while (*treestr) treestr++;	
		  *treestr++ = ',';
		}
	    }
	}
    }

  if(isTip(p->number, tr->rdta->numsp)) 
    {
      char *nameptr = tr->nameList[p->number];  

      if(subtreeSummary)
	{
	  for(i = 0; i < numberOfTipsForInsertion; i++)
	    {
	      if(p->bInf->executeThem[i] >= CUTOFF)
		foundTheSubtree[i] = TRUE;
	      else
		foundTheSubtree[i] = FALSE;
	    }	  
	}

        
      sprintf(treestr, "%s", nameptr);    
      while (*treestr) treestr++;
    }
  else 
    {                    
      *treestr++ = '(';     
      treestr = Tree2StringClassifyRec(treestr, tr, p->next->back, numberOfTipsForInsertion, 
				       countBranches, inserts,subtreeSummary,foundTheSubtree, CUTOFF, originalTree);     
      *treestr++ = ',';
      treestr = Tree2StringClassifyRec(treestr, tr, p->next->next->back, numberOfTipsForInsertion, 
				       countBranches, inserts, subtreeSummary,foundTheSubtree, CUTOFF, originalTree);          
      *treestr++ = ')';                    

      if(subtreeSummary)
	{  
	  for(i = 0; i < numberOfTipsForInsertion; i++)
	    {
	      if(p->bInf->executeThem[i] >= CUTOFF && 
		 p->next->back->bInf->executeThem[i] < CUTOFF && 
		 p->next->next->back->bInf->executeThem[i] < CUTOFF)
		foundTheSubtree[i] = TRUE;
	      else
		foundTheSubtree[i] = FALSE;
	    }
	}
    }
   
  if(!subtreeSummary && countQuery > 0)
    {
      sprintf(treestr, ":1.0[%s]", p->bInf->branchLabel);
      while (*treestr) treestr++;
      *treestr++ = ')'; 
    }
    
  if(originalTree)
    sprintf(treestr, ":%8.20f[%s", p->bInf->originalBranchLength, p->bInf->branchLabel);
  else
    sprintf(treestr, ":1.0[%s", p->bInf->branchLabel);
  while (*treestr) treestr++;

  assert(!(subtreeSummary == TRUE && originalTree == TRUE));
  assert(!(countQuery > 0 &&  originalTree == TRUE));

  if(subtreeSummary)
    {
      for(i = 0; i < numberOfTipsForInsertion; i++)
	{
	  if(foundTheSubtree[i])
	    {
	      sprintf(treestr," %s", tr->nameList[inserts[i]]);
	      while (*treestr) treestr++;
	    }
	}
    }

  sprintf(treestr, "]");            	 	        
  while (*treestr) treestr++;

  return  treestr;
}


static int findRoot(nodeptr p,  int numberOfTipsForInsertion, int ntips)
{
  if(isTip(p->number, ntips))
    {
      int i;
      for(i = 0; i < numberOfTipsForInsertion; i++)
	if(p->bInf->countThem[i] > 0)
	  return 0;

      return p->number;
    }
  else
    {
      int left;
      assert(p == p->next->next->next);

      left = findRoot(p->next->back, numberOfTipsForInsertion, ntips);

      if(left > 0)
	return left;
      else
	return findRoot(p->next->next->back, numberOfTipsForInsertion, ntips);
    }
}

static void calcSubtree(nodeptr p, int nTips, int numberOfTipsForInsertion)
{
  int i;

  if(isTip(p->number, nTips))    
    {
      for(i = 0; i < numberOfTipsForInsertion; i++)
	p->bInf->executeThem[i] = p->bInf->countThem[i];
    
      return;
    }
  else
    {
      nodeptr q;
      assert(p == p->next->next->next);
      
      q = p->next;

      while(q != p)
	{
	  calcSubtree(q->back, nTips, numberOfTipsForInsertion);	
	  q = q->next;
	}

     
      for(i = 0; i < numberOfTipsForInsertion; i++)
	p->bInf->executeThem[i] = p->next->back->bInf->executeThem[i] + 
	  p->next->next->back->bInf->executeThem[i] + p->bInf->countThem[i];
          
      return;
    }
}


static char *Tree2StringClassify(char *treestr, tree *tr, int numberOfTipsForInsertion, int *inserts, 
				 boolean subtreeSummary, int CUTOFF, int numberOfBranches, branchInfo *bInf, int root, 
				 boolean  originalTree)
{
  nodeptr p;
  int countBranches = 0; 

  if(!subtreeSummary)
    {           
      p = tr->start->back;
      
      assert(!isTip(p->number, tr->mxtips));
      
      *treestr++ = '(';
      treestr = Tree2StringClassifyRec(treestr, tr, p->back, numberOfTipsForInsertion, &countBranches, 
				       inserts, FALSE, (boolean *)NULL, -1, originalTree);
      *treestr++ = ',';
      treestr = Tree2StringClassifyRec(treestr, tr, p->next->back, numberOfTipsForInsertion, &countBranches, 
				       inserts, FALSE, (boolean *)NULL, -1, originalTree);
      *treestr++ = ',';
      treestr = Tree2StringClassifyRec(treestr, tr, p->next->next->back, numberOfTipsForInsertion, &countBranches, 
				       inserts, FALSE, (boolean *)NULL, -1, originalTree);
      *treestr++ = ')';
      *treestr++ = ';';
      
      assert(countBranches == 2 * tr->ntips - 3);
      
      *treestr++ = '\0';
      while (*treestr) treestr++;     
      return  treestr;
    }
  else
    {
      int i, j;
      boolean *foundTheSubtree = (boolean*)malloc(sizeof(boolean) * numberOfTipsForInsertion);

      assert(root > 0);
      assert(originalTree == FALSE);
      
     
      for(i = 0; i < numberOfBranches; i++)
	for(j = 0; j < numberOfTipsForInsertion; j++)
	  bInf[i].executeThem[j] = 0;
           
      p = tr->nodep[root];
	  
      calcSubtree(p, tr->mxtips, numberOfTipsForInsertion);
      calcSubtree(p->back, tr->mxtips, numberOfTipsForInsertion);     	      

      assert(isTip(p->number, tr->mxtips));           
      assert(!isTip(p->back->number, tr->mxtips));
	  
      *treestr++ = '(';
      *treestr++ = '(';	          
      sprintf(treestr, "%s", tr->nameList[p->number]);
      countBranches++;
      while (*treestr) treestr++;
      *treestr++ = ',';
      treestr = Tree2StringClassifyRec(treestr, tr, p->back, numberOfTipsForInsertion, &countBranches, 
				       inserts, TRUE, foundTheSubtree, CUTOFF, originalTree);      
      *treestr++ = ')';
      sprintf(treestr,"[ROOT]");
      while (*treestr) treestr++;
      *treestr++ = ')';
      *treestr++ = ';';
	  
      assert(countBranches == 2 * tr->ntips - 2);
      *treestr++ = '\0';
      while (*treestr) treestr++;
      free(foundTheSubtree);
      return  treestr;   
    }
}




static void markTips(nodeptr p, int *perm, int maxTips)
{
  if(isTip(p->number, maxTips))
    {
      perm[p->number] = 1;
      return;
    }
  else
    {
      nodeptr q = p->next;
      while(q != p)
	{
	  markTips(q->back, perm, maxTips);
	  q = q->next;
	}
      
    }
}


static inline double* makeDiagPtable(double z, int dataType, double *rates, double *eign)
{
  double *diagptable = (double*)NULL;

  switch(dataType)
    {
    case DNA_DATA:
      diagptable = (double *)malloc(sizeof(double) * 16);          
      break;
    case AA_DATA:
      diagptable = (double *)malloc(sizeof(double) * 80);     
      break;
    case BINARY_DATA:
      diagptable = (double *)malloc(sizeof(double) * 8);          
      break;
    case SECONDARY_DATA:
      diagptable = (double *)malloc(sizeof(double) * 64);      	      	     
      break;
    case SECONDARY_DATA_6:
      diagptable = (double *)malloc(sizeof(double) * 24);      	      	     
      break;  
    case SECONDARY_DATA_7:
      diagptable = (double *)malloc(sizeof(double) * 28);      	      	     
      break;   
    default:
      assert(0);
    }

  calcDiagptable(z, dataType, 4, rates, eign, diagptable);
  
  return diagptable;
}

static inline double scoreGeneric(double *x2, double *x1, double *diagptable, int scale, int dataType)
{
  double  term;    
  int     j, k; 

  term = 0.0;

  switch(dataType)
    {
    case AA_DATA:
      for(j = 0; j < 4; j++)
	for(k = 0; k < 20; k++)
	  term += x1[k] * x2[j * 20 + k] * diagptable[j * 20 + k];
      break;
    case DNA_DATA:
      for(j = 0; j < 4; j++)
	for(k = 0; k < 4; k++)
	  term += x1[k] * x2[j * 4 + k] * diagptable[j * 4 + k];
      break;
    case SECONDARY_DATA:
      for(j = 0; j < 4; j++)
	for(k = 0; k < 16; k++)
	  term += x1[k] * x2[j * 16 + k] * diagptable[j * 16 + k];
      break;
    case SECONDARY_DATA_6:
      for(j = 0; j < 4; j++)
	for(k = 0; k < 6; k++)
	  term += x1[k] * x2[j * 6 + k] * diagptable[j * 6 + k];
      break;  
    case SECONDARY_DATA_7:
      for(j = 0; j < 4; j++)
	for(k = 0; k < 7; k++)
	  term += x1[k] * x2[j * 7 + k] * diagptable[j * 7 + k];
      break;  
    case BINARY_DATA:
      for(j = 0; j < 4; j++)
	for(k = 0; k < 2; k++)
	  term += x1[k] * x2[j * 2 + k] * diagptable[j * 2 + k];
      break;
    default:
      assert(0);
    }
 
  term = log(0.25 * term) + scale * log(minlikelihood);

  return term;
}





#define DIAG 1
#define LEFT 2 




static double mlAlignGeneric(tree *tr, nodeptr p, unsigned char *raw, int unaligned, double *diagptable, boolean writeAlignment, 
			     double currentBest, int nodeNumber, boolean returnAlignment, unsigned char *tempAlignment)
{
  double **scoreMat, path1, path2, result;
  int i, j, nGaps, *refExp;  
  double *al  = &tr->partitionData[0].tipVector[0];
  double *gap;
  double *ref = tr->partitionData[0].xVector[p->number - tr->mxtips -1];
  double *gapVector;  
  char **backtrack = (char **)NULL; 
  int 
    tipOffset,
    vectorOffset;

  assert(tr->rateHetModel == GAMMA);

  switch(tr->partitionData[0].dataType)
    {
    case AA_DATA:     
      gap = &tr->partitionData[0].tipVector[440];
      tipOffset    = 20;
      vectorOffset = 80;
      break;
    case BINARY_DATA:
      gap = &tr->partitionData[0].tipVector[6];
      tipOffset    = 2;
      vectorOffset = 8;
      break;
    case SECONDARY_DATA:
      gap = &tr->partitionData[0].tipVector[4080];
      tipOffset    = 16;
      vectorOffset = 64;
      break;
    case SECONDARY_DATA_6:
      gap = &tr->partitionData[0].tipVector[4080];
      tipOffset    = 6;
      vectorOffset = 24;
      break;
    case SECONDARY_DATA_7:
      gap = &tr->partitionData[0].tipVector[4080];
      tipOffset    = 7;
      vectorOffset = 28;
      break;
    case DNA_DATA:      
      gap = &tr->partitionData[0].tipVector[60];
      tipOffset    = 4;
      vectorOffset = 16;
      break;
    default:
      assert(0);
    }

  assert((writeAlignment == FALSE && returnAlignment == FALSE) || (writeAlignment != returnAlignment));




  nGaps = tr->cdta->endsite - unaligned;

  refExp = tr->partitionData[0].expVector[p->number - tr->mxtips - 1];
 
  assert(nGaps > 0);

  scoreMat = (double **)malloc(sizeof(double*) * unaligned);

  for(i = 0; i < unaligned; i++)    
    scoreMat[i] = (double *)calloc(tr->cdta->endsite, sizeof(double)); 

  if(writeAlignment || returnAlignment)    
    {
      backtrack = (char **)malloc(unaligned * sizeof(char*));   
      for(i = 0; i < unaligned; i++)    
	backtrack[i] = (char *)calloc(tr->cdta->endsite, sizeof(char));
    } 

  gapVector = (double *)calloc(nGaps, sizeof(double));

  gapVector[0] = scoreGeneric(&ref[vectorOffset * 0], gap, diagptable, refExp[0], tr->partitionData[0].dataType);
  for(j = 1; j < nGaps; j++)
    gapVector[j] = gapVector[j-1] + scoreGeneric(&ref[vectorOffset * j], gap, diagptable, refExp[j], tr->partitionData[0].dataType);
   
  scoreMat[0][0] =  scoreGeneric(&ref[vectorOffset * 0], &al[tipOffset * raw[0]], diagptable, refExp[j], tr->partitionData[0].dataType);  
  if(writeAlignment || returnAlignment)
    backtrack[0][0] = DIAG;
 
  for(j = 1; j <= nGaps; j++)
    {          
      path1 = gapVector[j-1]   + scoreGeneric(&ref[vectorOffset * j], &al[tipOffset * raw[0]], diagptable, refExp[j], tr->partitionData[0].dataType);
      path2 = scoreMat[0][j-1] + scoreGeneric(&ref[vectorOffset * j], gap, diagptable, refExp[j], tr->partitionData[0].dataType);
      scoreMat[0][j] =  MAX(path1,path2);

      if(writeAlignment || returnAlignment)
	{
	  if(path1 > path2)
	    backtrack[0][j] = DIAG;
	  else
	    backtrack[0][j] = LEFT;
	}
    }          

  for(i = 1; i < unaligned; i++)
    for(j = i; j <= i + nGaps; j++)
      {
	assert(scoreMat[i-1][j-1] != 0);
	if(scoreMat[i][j-1] == 0)
	  {
	    scoreMat[i][j] = scoreMat[i-1][j-1] + scoreGeneric(&ref[vectorOffset * j], &al[tipOffset * raw[i]], diagptable, refExp[j], tr->partitionData[0].dataType);
	    if(writeAlignment || returnAlignment)
	      backtrack[i][j] = DIAG;
	  }
	else
	  {
	    path1 = scoreMat[i-1][j-1] + scoreGeneric(&ref[vectorOffset * j], &al[tipOffset * raw[i]], diagptable, refExp[j], tr->partitionData[0].dataType);
	    path2 = scoreMat[i][j-1]   + scoreGeneric(&ref[vectorOffset * j], gap,                     diagptable, refExp[j], tr->partitionData[0].dataType);
	    scoreMat[i][j] = MAX(path1, path2);

	    if(writeAlignment || returnAlignment)
	      {
		if(path1 > path2)
		  backtrack[i][j] = DIAG;
		else
		  backtrack[i][j] = LEFT;
	      }	    
	  }	    	
      }
     
  result = scoreMat[unaligned - 1][ tr->cdta->endsite - 1];

  if((writeAlignment && result > currentBest) || returnAlignment)
    {

      unsigned char *tip;      
     

      int counterTip =  tr->cdta->endsite - 1;

      if(writeAlignment)
	tip = tr->yVector[nodeNumber];
      else
	tip = tempAlignment;

      assert(isTip(nodeNumber, tr->mxtips));

      i = unaligned - 1;
      j = tr->cdta->endsite - 1;



      while(i >= 0 && j >= 0)
	{
	  assert(backtrack[i][j] != 0);
	  if(backtrack[i][j] == DIAG)
	    {	     

	      tip[counterTip--]   = raw[i];
	      i--;
	      j--;
	    }
	  else
	    {	      

	      tip[counterTip--]   = 15;
	      j--;
	    }
	}

      while(j >= 0)
	{	 
	  tip[counterTip--]   = 15;	 
	  j--;
	}

      assert(counterTip == -1);      
    }
  
  if(writeAlignment || returnAlignment)
    {
      for(i = 0; i < unaligned; i++)
	free(backtrack[i]);

      free(backtrack);
    }
  
  for(i = 0; i < unaligned; i++)
    free(scoreMat[i]);
  free(scoreMat);
   
  free(gapVector);

  return result;
}




#define _DEBUG_INS
 



static void addInsertion(nodeptr q, double result, int heuristicInsertions, bestBranches *bBr, int j)
{ 
  int i;
  double min;

  /*printf("Adding %f\n", result);*/

  assert(heuristicInsertions > 0);

  if(bBr[j].count < heuristicInsertions)
    {
      bBr[j].likelihoods[bBr[j].count] = result;
      bBr[j].branch[bBr[j].count]      = q->bInf;
      bBr[j].count          = bBr[j].count + 1;
      
      if(bBr[j].count == heuristicInsertions)
	{      
	  min = bBr[j].likelihoods[0];
	  bBr[j].minLH = 0;

	  for(i = 1; i < heuristicInsertions; i++)
	    {
	      if(bBr[j].likelihoods[i] < min)
		{
		  min = bBr[j].likelihoods[i];
		  bBr[j].minLH = i;
		}
	    }
	  
#ifdef _DEBUG_INS
	  for(i = 0; i < heuristicInsertions; i++)
	    assert(min <=  bBr[j].likelihoods[i]);
#endif
	}

      return;
    }

  if(result < bBr[j].likelihoods[bBr[j].minLH])
    return;
  
  bBr[j].likelihoods[bBr[j].minLH] = result;
  bBr[j].branch[bBr[j].minLH]      = q->bInf;
          
	 
  min = bBr[j].likelihoods[0];
  bBr[j].minLH = 0;
  for(i = 1; i < heuristicInsertions; i++)
    {
      if(bBr[j].likelihoods[i] < min)
	{
	  min = bBr[j].likelihoods[i];
	  bBr[j].minLH = i;
	}
    }

#ifdef _DEBUG_INS
  /*{
    boolean FAIL = FALSE;

    for(i = 0; i < heuristicInsertions; i++)
      if(!((min <  bBr[j].likelihoods[i]) || ((min == bBr[j].likelihoods[i]) && (bBr[j].minLH == i))))
	FAIL = TRUE;
    
    if(FAIL)
      for(i = 0; i < heuristicInsertions; i++)
	printf("%s %d %d %1.80f %1.80f\n", bBr[j].branch[i]->branchLabel, i, bBr[j].minLH, bBr[j].likelihoods[i], bBr[j].likelihoods[bBr[j].minLH]);
  */
    
  for(i = 0; i < heuristicInsertions; i++)
    assert(min <=  bBr[j].likelihoods[i]);

    /*}*/
#endif


  return; 
}


static void testInsertThorough(tree *tr, nodeptr r, nodeptr q, int *inserts,  int numberOfTipsForInsertion,
			       double *insertLikelihoods, double *branches, nodeptr *insertNodes, boolean align, double *diagptable, 
			       boolean writeAlignment, unalignedSequences *rawSeqs, int dataType, bestBranches *bBr,
			       int heuristicInsertions, boolean bootstrap)
{
  double result;           
  double  qz[NUM_BRANCHES];
  double  zqs[NUM_BRANCHES], zrs[NUM_BRANCHES], lzqr, lzqs, lzrs, lzsum, lzq, lzr, lzs, lzmax, zz;      
  double defaultArray[NUM_BRANCHES];	
  double e1[NUM_BRANCHES], e2[NUM_BRANCHES], e3[NUM_BRANCHES];
  nodeptr  x;         
  int i, j;

  assert(!tr->grouped);

  x = q->back; 

  /* we will insert node r between nodes x and q */
  
  /* first save branch x<->q */

  for(i = 0; i < tr->numBranches; i++)    
    {
      qz[i] = q->z[i];
      defaultArray[i] = defaultz;
    }
  
  for(j = 0; j < numberOfTipsForInsertion; j++)
    {          
      if(heuristicInsertions >= 0 || (heuristicInsertions == -1 && q->bInf->executeThem[j] == 1))
	{
	  if(align)
	    {
	      double 
		tempz[NUM_BRANCHES],
		*diagptableTwo = (double *)NULL;
	      unsigned char *tempAlignment  = (unsigned char *)malloc(sizeof(unsigned char) * tr->cdta->endsite);
	      unsigned char *sequenceBuffer = (unsigned char *)malloc(sizeof(unsigned char) * tr->cdta->endsite);
	      
	      nodeptr s = tr->nodep[inserts[j]];
	      
	      assert(!bootstrap);

	      for(i = 0; i < tr->numBranches; i++)
		{		      
		  tempz[i] = sqrt(qz[i]);      
		  
		  if(tempz[i] < zmin) 
		    tempz[i] = zmin;
		  if(tempz[i] > zmax)
		    tempz[i] = zmax;
		}   
	      
	      hookup(r->next,       q, tempz, tr->numBranches);
	      hookup(r->next->next, x, tempz, tr->numBranches);	
	      newviewGeneric(tr, r);
	      
	      mlAlignGeneric(tr, r, rawSeqs[j].seq, rawSeqs[j].length, diagptable, FALSE, insertLikelihoods[j],
			     inserts[j], TRUE, tempAlignment);
	     

	      /*printf("INITIAL ALIGN\t %d %f\n", *count, result);*/
	      
	      assert(s->number == inserts[j]);
	      
	      memcpy(sequenceBuffer,          tr->yVector[inserts[j]], sizeof(char) * tr->cdta->endsite);
	      memcpy(tr->yVector[inserts[j]], tempAlignment,           sizeof(char) * tr->cdta->endsite);	 
	      
	      makenewzGeneric(tr, q, s, defaultArray, iterations, zqs, FALSE);                  
	      makenewzGeneric(tr, x, s, defaultArray, iterations, zrs, FALSE);
	      
	      for(i = 0; i < tr->numBranches; i++)
		{
		  lzqr = (qz[i] > zmin) ? log(qz[i]) : log(zmin); 
		  lzqs = (zqs[i] > zmin) ? log(zqs[i]) : log(zmin);
		  lzrs = (zrs[i] > zmin) ? log(zrs[i]) : log(zmin);
		  lzsum = 0.5 * (lzqr + lzqs + lzrs);
		  
		  lzq = lzsum - lzrs;
		  lzr = lzsum - lzqs;
		  lzs = lzsum - lzqr;
		  lzmax = log(zmax);
		  
		  if      (lzq > lzmax) {lzq = lzmax; lzr = lzqr; lzs = lzqs;} 
		  else if (lzr > lzmax) {lzr = lzmax; lzq = lzqr; lzs = lzrs;}
		  else if (lzs > lzmax) {lzs = lzmax; lzq = lzqs; lzr = lzrs;}          
		  
		  e1[i] = exp(lzq);
		  e2[i] = exp(lzr);
		  e3[i] = exp(lzs);
		}	   	 
	      
	      hookup(r->next,       q, e1, tr->numBranches);
	      hookup(r->next->next, x, e2, tr->numBranches);
	      hookup(r,             s, e3, tr->numBranches); 
	      
	      newviewGeneric(tr, r);	   	  
	      localSmooth(tr, r, smoothings);
	      result = evaluateGeneric(tr, r);
	      
	      /*printf("BR-OPT \t\t %d %f\n", *count, result);*/
	      
	      memcpy(tr->yVector[inserts[j]], sequenceBuffer, sizeof(char) * tr->cdta->endsite);
	      
	      newviewGeneric(tr, r);
	      	     
	      zz = r->z[0];
	     
	      diagptableTwo =  makeDiagPtable(zz, dataType, tr->partitionData[0].gammaRates,tr->partitionData[0].EIGN);	 
	      result = mlAlignGeneric(tr, r, rawSeqs[j].seq, rawSeqs[j].length, diagptable, writeAlignment, insertLikelihoods[j],
				      inserts[j], FALSE, (unsigned char*)NULL);
	     	     
	       
	      /*printf("RE-ALIGN \t %d %f\n", *count, result);*/

	      if(result > insertLikelihoods[j])
		{	     
		  insertLikelihoods[j] = result;
		  insertNodes[j]       = q;	     		 
		}      	   	    	    
	      	      
	      if(heuristicInsertions > 0)
		addInsertion(q, result, heuristicInsertions, bBr, j);

	      free(diagptableTwo);	      
	      free(tempAlignment);
	      free(sequenceBuffer);
	    }
	  else
	    {    
	      nodeptr s = tr->nodep[inserts[j]];	      	     
	     
	      makenewzGeneric(tr, q, s, defaultArray, iterations, zqs, FALSE);                  
	      makenewzGeneric(tr, x, s, defaultArray, iterations, zrs, FALSE);
	      
	      for(i = 0; i < tr->numBranches; i++)
		{
		  lzqr = (qz[i] > zmin) ? log(qz[i]) : log(zmin);		  
		  lzqs = (zqs[i] > zmin) ? log(zqs[i]) : log(zmin);
		  lzrs = (zrs[i] > zmin) ? log(zrs[i]) : log(zmin);
		  lzsum = 0.5 * (lzqr + lzqs + lzrs);
		  
		  lzq = lzsum - lzrs;
		  lzr = lzsum - lzqs;
		  lzs = lzsum - lzqr;
		  lzmax = log(zmax);
		  
		  if      (lzq > lzmax) {lzq = lzmax; lzr = lzqr; lzs = lzqs;} 
		  else if (lzr > lzmax) {lzr = lzmax; lzq = lzqr; lzs = lzrs;}
		  else if (lzs > lzmax) {lzs = lzmax; lzq = lzqs; lzr = lzrs;}          
		  
		  e1[i] = exp(lzq);
		  e2[i] = exp(lzr);
		  e3[i] = exp(lzs);
		}
	      
	      hookup(r->next,       q, e1, tr->numBranches);
	      hookup(r->next->next, x, e2, tr->numBranches);
	      hookup(r,             s, e3, tr->numBranches);      		     
	      
	      newviewGeneric(tr, r);
	      
	      localSmooth(tr, r, smoothings);
	      result = evaluateGeneric(tr, r);
	      
#ifdef _DEBUG_CLASSIFY	      
	      printf("%f\n", result);
#endif	      
	      
	      if(result > insertLikelihoods[j])
		{	     
		  insertLikelihoods[j] = result;
		  insertNodes[j]       = q;
		  if(bootstrap)
		    {
		      assert(branches);
		      branches[j] = getBranch(tr, r->z, r->back->z);/* + getBranch(tr, q->z, q->back->z) + getBranch(tr, x->z, x->back->z);*/
		    }
		} 
	      
	      if(heuristicInsertions > 0)
		{
		  /* printf("%f \n", result);*/
		  addInsertion(q, result, heuristicInsertions, bBr, j);
		}
	    }
	  /*hookup(q, x, qz, tr->numBranches);*/
	}
    }
  /* repair branch q <-> x */

  hookup(q, x, qz, tr->numBranches);
  
  
  r->next->next->back = r->next->back = (nodeptr) NULL; 
}

static void testInsertFast(tree *tr, nodeptr r, nodeptr q, int *inserts,  int numberOfTipsForInsertion, 
			   double *insertLikelihoods, nodeptr *insertNodes, boolean align, double *diagptable, 
			   boolean writeAlignment, unalignedSequences *rawSeqs, bestBranches *bBr,
			   int heuristicInsertions)
{
  double result;           
  double  qz[NUM_BRANCHES], z[NUM_BRANCHES];
  nodeptr  x;       
  int i;
  
  x = q->back;     	       
  
  assert(!tr->grouped);                    
  
  for(i = 0; i < tr->numBranches; i++)
    {	
      qz[i] = q->z[i];
      z[i] = sqrt(q->z[i]);      
      
      if(z[i] < zmin) 
	z[i] = zmin;
      if(z[i] > zmax)
	z[i] = zmax;
    }        
  
  hookup(r->next,       q, z, tr->numBranches);
  hookup(r->next->next, x, z, tr->numBranches);	                         
  
  newviewGeneric(tr, r);   
  
  for(i = 0; i < numberOfTipsForInsertion; i++)
    {
      if(heuristicInsertions >= 0 || (heuristicInsertions == -1 && q->bInf->executeThem[i] == 1))
	{
	  if(align)	
	    result = mlAlignGeneric(tr, r, rawSeqs[i].seq, rawSeqs[i].length, diagptable, writeAlignment, insertLikelihoods[i],
				    inserts[i], FALSE, (unsigned char*)NULL);	   	
	  else
	    {	 
	      hookupDefault(r, tr->nodep[inserts[i]], tr->numBranches);
	      result = evaluateGeneric (tr, r);
	      
	      r->back = (nodeptr) NULL;
	      tr->nodep[inserts[i]]->back = (nodeptr) NULL;
	    }
	  
#ifdef _DEBUG_CLASSIFY
	  printf("%f\n", result);
#endif
	  
	  if(result > insertLikelihoods[i])
	    {	     
	      insertLikelihoods[i] = result;
	      insertNodes[i]       = q;	     
	    } 
	  
	  if(heuristicInsertions > 0)
	    addInsertion(q, result, heuristicInsertions, bBr, i);
    
	}    
    }
  
  hookup(q, x, qz, tr->numBranches);
  
  r->next->next->back = r->next->back = (nodeptr) NULL;
}

static void addTraverseRob(tree *tr, nodeptr r, nodeptr q, int *inserts,  int numberOfTipsForInsertion,
			   double *insertLikelihoods, double *branches, nodeptr *insertNodes, boolean align, double *diagptable, 
			   boolean writeAlignment, unalignedSequences *rawSeqs, boolean thorough, int dataType, bestBranches *bBr,
			   int heuristicInsertions, boolean bootstrap)
{       
  /*{
    int i;

    for(i = 0; i < numberOfBranches; i++)
      {
	if(thorough)
	  testInsertThorough(tr, r, q, inserts,  numberOfTipsForInsertion, insertLikelihoods, insertNodes, align,
			     diagptable, writeAlignment, rawSeqs, dataType, bBr, heuristicInsertions);
	else
	  testInsertFast(tr, r, q, inserts,  numberOfTipsForInsertion, insertLikelihoods, insertNodes, align,
			 diagptable, writeAlignment, rawSeqs, bBr, heuristicInsertions);
      }
      }*/


  if(thorough)
    testInsertThorough(tr, r, q, inserts,  numberOfTipsForInsertion, insertLikelihoods, branches, insertNodes, align,
		       diagptable, writeAlignment, rawSeqs, dataType, bBr, heuristicInsertions, bootstrap);
  else
    testInsertFast(tr, r, q, inserts,  numberOfTipsForInsertion, insertLikelihoods, insertNodes, align,
		   diagptable, writeAlignment, rawSeqs, bBr, heuristicInsertions);

  if(!isTip(q->number, tr->rdta->numsp))
    {   
      nodeptr a = q->next;

      while(a != q)
	{
	  addTraverseRob(tr, r, a->back,       inserts,  numberOfTipsForInsertion, insertLikelihoods, branches,
			 insertNodes, align, diagptable, writeAlignment, rawSeqs, thorough, dataType, bBr, 
			 heuristicInsertions, bootstrap);
	  a = a->next;
	}      
    }
} 

#ifdef _USE_PTHREADS_MULTIGRAIN

static size_t getContiguousVectorLength(tree *tr)
{
  size_t length = 0;
  int model;

  for(model = 0; model < tr->NumberOfModels; model++)
    {     
      size_t realWidth = tr->partitionData[model].upper - tr->partitionData[model].lower;

      switch(tr->partitionData[model].dataType) 
	{
	case AA_DATA:
	  length += realWidth * (size_t)tr->aaIncrement;
	  break;
	case DNA_DATA:
	  length += realWidth * (size_t)tr->dnaIncrement ;
	  break;
	case BINARY_DATA:
	  length += realWidth * (size_t)tr->binaryIncrement ;
	  break;
	case SECONDARY_DATA:
	  length += realWidth * (size_t)tr->secondaryIncrement ;
	  break;
	case SECONDARY_DATA_6:
	  length += realWidth * (size_t)tr->secondaryIncrement6;
	  break; 
	case SECONDARY_DATA_7:
	  length += realWidth * (size_t)tr->secondaryIncrement7;
	  break;
	default:
	  assert(0);
	} 	
    }

  return length;
}
static size_t getContiguousScalingLength(tree *tr)
{
  size_t length = 0;
  int model;

  for(model = 0; model < tr->NumberOfModels; model++)    
    length += tr->partitionData[model].upper - tr->partitionData[model].lower;

  return length;
}


#endif


/* a little likelihood vector counter to check for correctness afterwards */

int vectCount = 0;

static void setupBranchMetaInfo(tree *tr, int *countBranches, nodeptr p, int nTips, branchInfo *bInf, size_t contiguousVectorLength, size_t contiguousScalingLength)
{

  /* our current node is a tip in the tree */

  if(isTip(p->number, nTips))    
    {

      /* link the two node pointers that define this branch to the 
	 branchInfo data structure array */

      p->bInf       = &bInf[*countBranches];
      p->back->bInf = &bInf[*countBranches];               
	

      /* store the branch length, mainly used for printing out results afterwards */

      bInf[*countBranches].originalBranchLength = getBranch(tr, p->z, p->back->z);
      

#ifdef _USE_PTHREADS_MULTIGRAIN     
      {
	/* now in here we actually fo the assignment of the memory space 
	   and THREAD_GATHER operation */

	int i;
	
	/* set a pointer to the current branchInfo array entry for easier access */
	branchInfo *b = &bInf[*countBranches];

	/* store which branch number we are working on in a respective field of the master thread */
	/* we will need to then copy this to the worker thread */

	tr->branchNumber = *countBranches;

	/* 
	   store the node number of the tree nodes to the left and the right of the current branch 
	   tip nodes are numbered from 1, 2, ...., tr->mxtips, inner nodes from tr->mxtips + 1, tr->mxtips + 2,
	   ...., 2 * tr->mxtips - 3
	   
	 */

	b->leftNodeNumber = p->number;
	b->rightNodeNumber = p->back->number;

	/* 
	   here we have a branch leading to a tip, hence we only need to assign 
	   the likelihood vector at that end of the branch that has an inner node 
	*/

	b->left  = (double*)NULL; 
	b->leftScaling = (int*)NULL;


	/* assign the contiguous likelihood vector */
	b->right = (double*)malloc(sizeof(double) * contiguousVectorLength);

	/* assign the contiguous scaling vector */
	b->rightScaling = (int*)malloc(sizeof(int) * contiguousScalingLength);

	/* increment the vector count */
	vectCount++;

	/* check if the likelihood vector of the inner node already points to the present branch,
	   i.e., if the virtual root was previsouly located in the current branch.
	   If it is not we need to re-compute it by calling newviewGeneric() which computes 
	   a new view, i.e., a new orientation of likelihood vectors which always have a rooted 
	   view of the tree towards our current virtual root that is located in the current branch we are 
	   visiting */


	if(!p->back->x)
	  newviewGeneric(tr, p->back);

	/* assert that the likelihood vector is now available at the correct node 
	   of th cyclic list of three nodeptr structures that reresent an inner 
	   node. In this cyclic list with p->back , p->back->next, and p->back->next->next, note that 
	   p->back->next->next->next == p->back only one of the ->x will be set to 1 while the remaining two 
	   will be set to zero */


	assert(p->back->x);
       
	
	/* store the per-partition branch lengths in the branch Info data structure */
	for(i = 0; i < tr->numBranches; i++)
	  b->branchLengths[i] = p->z[i];

	/* and here we now do our gather operation to copy the likelihood vector and scaling vector entries 
	   that are scattered accross the therads iunto one contiguous vector, uff */

	masterBarrier(THREAD_GATHER_LIKELIHOOD, tr);
      }     
#endif

      *countBranches =  *countBranches + 1;
      return;
    }
  else
    {
      nodeptr q;
      assert(p == p->next->next->next);

      p->bInf       = &bInf[*countBranches];
      p->back->bInf = &bInf[*countBranches];

      bInf[*countBranches].originalBranchLength = getBranch(tr, p->z, p->back->z);           
      
#ifdef _USE_PTHREADS_MULTIGRAIN     
      {
	int i;
	branchInfo *b = &bInf[*countBranches];
	tr->branchNumber = *countBranches;
	b->leftNodeNumber = p->number;
	b->rightNodeNumber = p->back->number;	
	
	
	/* p can't be a tip, because otherwise we must be in the other case of this function */

	assert(!isTip(p->number, tr->mxtips));
	{
	  /* get likelihood vector */
	  if(!p->x)
	    newviewGeneric(tr, p);
	  assert(p->x);

	  /* assign space */

	  b->left  = (double*)malloc(sizeof(double) * contiguousVectorLength);
	  b->leftScaling = (int*)malloc(sizeof(int) * contiguousScalingLength);
	  vectCount++;
	}
	 

	/* now here we have to check if the node at the other end of the branch defined by 
	   p and p->back note that p = p->back->back isn't a tip.
	   if it is a tip node we will not do anything, otherwise we will just assign 
	   memory as well */

	if(!isTip(p->back->number, tr->mxtips))
	  {
	    if(!p->back->x)
	      newviewGeneric(tr, p->back);
	    assert(p->back->x);
	    b->right = (double*)malloc(sizeof(double) * contiguousVectorLength);
	    b->rightScaling = (int*)malloc(sizeof(int) * contiguousScalingLength);
	    vectCount++;
	  }
	else
	  {
	    b->right  = (double*)NULL; 
	    b->rightScaling = (int*)NULL;	    
	  }
	
	/* store branches */

	for(i = 0; i < tr->numBranches; i++)
	  b->branchLengths[i] = p->z[i];

	/* do the thread gather operation */

	masterBarrier(THREAD_GATHER_LIKELIHOOD, tr);
	
      }
#endif 

      /* increment the branch counter */
      *countBranches =  *countBranches + 1;
      

      /* recursively descend into the left and right child 
	 of this node */

      q = p->next;

      while(q != p)
	{
	  setupBranchMetaInfo(tr, countBranches, q->back, nTips, bInf, contiguousVectorLength, contiguousScalingLength);	
	  q = q->next;
	}
     
      return;
    }
}
 
static void analyzeQuerySeqs(tree *tr, int *inserts, int numberOfTipsForInsertion)
{
  /* only for dynamic alignment version, check if there are gaps at all in the 
     query seqs */

  int 
    i, 
    model, 
    j,
    totalCount = 0;  

  for(i = 0; i < numberOfTipsForInsertion; i++)
    {
      int countGaps = 0;

      for(model = 0; model < tr->NumberOfModels; model++)
	{	  
	  unsigned char *tip = tr->partitionData[model].yVector[inserts[i]];
	  int lower = tr->partitionData[model].lower;
	  int upper = tr->partitionData[model].upper;

	  switch(tr->partitionData[model].dataType)
	    {
	    case AA_DATA:
	      for(j = lower; j < upper; j++)
		if(tip[j] == UNDETERMINED_AA)
		  countGaps++;
	      break;
	    case DNA_DATA:
	      for(j = lower; j < upper; j++)
		if(tip[j] == UNDETERMINED_DNA)
		  countGaps++;		 
	      break;
	    case BINARY_DATA:
	      for(j = lower; j < upper; j++)
		if(tip[j] == UNDETERMINED_BINARY)
		  countGaps++;		 
	      break;
	    case SECONDARY_DATA:
	      for(j = lower; j < upper; j++)
		if(tip[j] == UNDETERMINED_SECONDARY)
		  countGaps++;		 
	      break;
	    case SECONDARY_DATA_6:
	      for(j = lower; j < upper; j++)
		if(tip[j] == UNDETERMINED_SECONDARY_6)
		  countGaps++;		 
	      break; 
	    case SECONDARY_DATA_7:
	      for(j = lower; j < upper; j++)
		if(tip[j] == UNDETERMINED_SECONDARY_7)
		  countGaps++;		 
	      break;
	    default:
	      assert(0);
	    }	 
	}            
     
      if(countGaps == 0)
	{
	  printf("Query sequences %s you want to dynamically align to the reference alignment\n", tr->nameList[inserts[i]]);
	  printf("does not contain any gaps, i.e., there is nothing to align here\n");
	  totalCount++;
	}	
    }

  if(totalCount > 0)
    {
      printf("found %d query seqs for dynamic alignment that do not contain any gaps, exiting ... \n", totalCount);
      exit(-1);
    }

}

static void branchLabels(tree *tr)
{  
  char originalLabelledTreeFileName[1024];

  FILE *treeFile;    
 
  strcpy(originalLabelledTreeFileName, workdir); 
  strcat(originalLabelledTreeFileName, "RAxML_originalLabelledTree.");  
  strcat(originalLabelledTreeFileName, run_id);
  
  /* TODO need to actually think about this, a seg fault waiting to happen here */
  
  free(tr->tree_string);
  tr->treeStringLength *= 16;
  tr->tree_string  = (char*)malloc(tr->treeStringLength * sizeof(char));

  treeFile = myfopen(originalLabelledTreeFileName, "w");
  Tree2StringClassify(tr->tree_string, tr, 0, (int*)NULL, FALSE, -1, -1, (branchInfo*)NULL, -1, TRUE);
  fprintf(treeFile, "%s\n", tr->tree_string);    
  fclose(treeFile);

  exit(1);
}





void classifyML(tree *tr, analdef *adef)
{
  int
    root = -1,
    i, 
    j,  
    *perm, 
    *inserts, 
    numberOfTipsForInsertion = 0, 
    countBranches = 0,
    heuristicInsertions = 0, 
    numberOfBranches = 2 * tr->ntips - 3,
    dataType;

  bestBranches *bBr = (bestBranches *)NULL;
  
  double  
    t,
    startLH, 
    *insertLikelihoods,
    *branches;
  
  nodeptr 
    *insertNodes, 
    r, 
    q;
  
  unalignedSequences *rawSeqs = (unalignedSequences*)NULL; 
  
  branchInfo *bInf;
  
  boolean 
    align          = adef->dynamicAlignment,
    writeAlignment,
    thorough       = adef->thoroughInsertion,
    useHeuristics  = TRUE;
  
  char 
    labelledTreeFileName[1024],
    originalLabelledTreeFileName[1024],
    classificationFileName[1024],
    alignmentFileName[1024];

  FILE 
    *treeFile, 
    *classificationFile, 
    *alignmentFile;

#ifdef _USE_PTHREADS
  assert(!adef->dynamicAlignment);
#endif

  strcpy(labelledTreeFileName,         workdir);
  strcpy(originalLabelledTreeFileName, workdir);
  strcpy(classificationFileName,       workdir);
  strcpy(alignmentFileName,            workdir);
  
  strcat(labelledTreeFileName,         "RAxML_labelledTree.");
  strcat(originalLabelledTreeFileName, "RAxML_originalLabelledTree.");
  strcat(classificationFileName,       "RAxML_classification.");
  strcat(alignmentFileName,            "RAxML_alignment.");
  
  strcat(labelledTreeFileName,         run_id);
  strcat(originalLabelledTreeFileName, run_id);
  strcat(classificationFileName,       run_id);
  strcat(alignmentFileName,            run_id);

  /* dynamic alignment for GAMMA and unpartitioned datasets only */
 
  if(adef->dynamicAlignment)
    assert(tr->NumberOfModels == 1 && (tr->rateHetModel == GAMMA || tr->rateHetModel == GAMMA_I));

  if(adef->dynamicAlignment)
    writeAlignment = TRUE;
  else
    writeAlignment = FALSE;

  assert(adef->restart);
 
  /*printf("Current Tips: %d length %d\n", tr->ntips, tr->cdta->endsite);*/
  printBothOpen("\nRAxML Classification Algorithm search parameters:\n");
  printBothOpen("Insertion using dynamic Alignment: %s\n", align?"YES":"NO");
  printBothOpen("Thorough  Insertion Method with branch length optimization: %s\n", thorough?"YES":"NO");
  printBothOpen("Bootstrap insertions with heuristics: %s\n\n", (adef->rapidBoot && useHeuristics)?"YES":"NO"); 
  
  adef->outgroup = FALSE;
  tr->doCutoff   = FALSE;

  evaluateGenericInitrav(tr, tr->start);

#ifdef _DEBUG_CLASSIFY
  printf("Init: %f\n", tr->likelihood);  
#endif

  if(align)
    {
      assert(tr->NumberOfModels == 1);
      switch(tr->partitionData[0].dataType)
	{
	case AA_DATA:
	  dataType = AA_DATA;
	  break;
	case DNA_DATA:
	  dataType = DNA_DATA;
	  break;
	case BINARY_DATA:
	  dataType = BINARY_DATA;
	  break;
	case SECONDARY_DATA:
	  dataType = SECONDARY_DATA;
	  break;
	case SECONDARY_DATA_6:
	  dataType = SECONDARY_DATA_6;
	  break;
	case SECONDARY_DATA_7:
	  dataType = SECONDARY_DATA_7;
	  break;
	default:
	  assert(0);
	}
    }
  else
    dataType = -1;
     
  modOpt(tr, adef, TRUE, 1.0);    
  
  

  switch(tr->rateHetModel)
    {
    case GAMMA:
      printBothOpen("GAMMA based Likelihood score of reference tree: %f\n\n", tr->likelihood);
      break;
    case GAMMA_I:
      printBothOpen("GAMMA+P-Invar based Likelihood score of reference tree: %f\n\n", tr->likelihood);
      break; 
    case CAT:
      printBothOpen("CAT based Likelihood score of reference tree: %f\n\n", tr->likelihood);
      break;
    default:
      assert(0);
    }

  startLH = tr->likelihood; 

  perm    = (int *)calloc(tr->mxtips + 1, sizeof(int));
  inserts = (int *)calloc(tr->mxtips, sizeof(int));

  markTips(tr->start,       perm, tr->mxtips);
  markTips(tr->start->back, perm ,tr->mxtips);

  numberOfTipsForInsertion = 0;

  for(i = 1; i <= tr->mxtips; i++)
    {
      if(perm[i] == 0)
	{
	  inserts[numberOfTipsForInsertion] = i;
	  numberOfTipsForInsertion++;
	}
    }      
 
  if(align)
    { 
      /* check if there are gaps in query seqs */     

      analyzeQuerySeqs(tr, inserts, numberOfTipsForInsertion);

      rawSeqs = (unalignedSequences*)malloc(sizeof(unalignedSequences) * numberOfTipsForInsertion);

      for(i = 0; i < numberOfTipsForInsertion; i++)
	{
	  int unaligned = 0;
	  unsigned char *tip, *raw;

	  rawSeqs[i].seq = raw = (unsigned char *)malloc(sizeof(unsigned char) * tr->cdta->endsite);
	  tip = tr->yVector[inserts[i]];     

	  /* TODO this is a bit stupid, it implies that AA could be 
	     aligned with DNA data */

	  for(j = 0; j < tr->cdta->endsite; j++)
	    {
	      if(tr->dataVector[j] == DNA_DATA && tip[j] != UNDETERMINED_DNA)	    	    
		raw[unaligned++] = tip[j];	     	   
	    
	      if(tr->dataVector[j] == AA_DATA && tip[j] != UNDETERMINED_AA)	    	      
		raw[unaligned++] = tip[j];	

	      if(tr->dataVector[j] == BINARY_DATA && tip[j] != UNDETERMINED_BINARY)	    	      
		raw[unaligned++] = tip[j];

	      if(tr->dataVector[j] == SECONDARY_DATA && tip[j] != UNDETERMINED_SECONDARY)	    	      
		raw[unaligned++] = tip[j];

	      if(tr->dataVector[j] == SECONDARY_DATA_6 && tip[j] != UNDETERMINED_SECONDARY_6)	    	      
		raw[unaligned++] = tip[j];

	      if(tr->dataVector[j] == SECONDARY_DATA_7 && tip[j] != UNDETERMINED_SECONDARY_7)	    	      
		raw[unaligned++] = tip[j];	    
	    }   
	 
	  rawSeqs[i].length = unaligned;
	}
	  
    }

  printBothOpen("RAxML will classify %d Query Sequences into the %d branches of the reference tree with %d taxa\n\n",  numberOfTipsForInsertion, (2 * tr->ntips - 3), tr->ntips);

  if(adef->rapidBoot)
    printBothOpen("RAxML will execute %d Bootstrap replicates to infer classification support\n\n",  adef->multipleRuns);

  assert(numberOfTipsForInsertion == (tr->mxtips - tr->ntips));

  branches          = (double*)malloc(sizeof(double) *  numberOfTipsForInsertion);
  insertLikelihoods = (double*)malloc(sizeof(double) *  numberOfTipsForInsertion);
  insertNodes       = (nodeptr*)malloc(sizeof(nodeptr) *  numberOfTipsForInsertion);      

  bInf              = (branchInfo*)malloc(numberOfBranches * sizeof(branchInfo));

  for(i = 0; i < numberOfBranches; i++)
    {      
      bInf[i].countThem   = (int*)calloc(numberOfTipsForInsertion, sizeof(int));      
      bInf[i].executeThem = (int*)calloc(numberOfTipsForInsertion, sizeof(int));
      bInf[i].branches    = (double*)calloc(numberOfTipsForInsertion, sizeof(double));
      sprintf(bInf[i].branchLabel, "I%d", i);
    }
   


  /* setup nodes for insertion and traversals */
  /* TODO can free a large chunk of memory not required
     for query sequences */

  r = tr->nodep[(tr->nextnode)++];     
  q = findAnyTip(tr->start, tr->rdta->numsp);
  assert(isTip(q->number, tr->rdta->numsp));
  assert(!isTip(q->back->number, tr->rdta->numsp));
  /*printf("Virtual Root %s\n", tr->nameList[q->number]);*/
	 
  q = q->back;               

 

#ifdef _USE_PTHREADS_MULTIGRAIN
  {
    /* here we just compute the length, in terms of columns 
       of a full length contiguous vector AND the length 
       of the contigoous likelihood array in terms of double 
       values we will need contiguousScalingLength to assign 
       memory space for the scaling vectors and 
       contiguousVectorLength to assign space for the contiguous 
       likelihood vectors that are linked to the branch info structure
    */
       

    size_t contiguousVectorLength = getContiguousVectorLength(tr);
    size_t contiguousScalingLength = getContiguousScalingLength(tr);


    /* 
       initially we just need to store some data in the tree data structure:
       1. the pointer to the branch info array
       2. the number of branches in the reference tree

       That stuff, along with some other pointers that are already stored 
       in the tree data structure of the master thread, will then need to be written
       to the respective fields in the worker thread tree data structures when 
       function  masterBarrier(THREAD_COPY_BRANCHINFO_POINTER, tr);
       is called.       
    */
       

    tr->bInfo = bInf;
    tr->numberOfBranches = numberOfBranches;
    
   

    /* 
       here we just copy some of the values we will need from the master 
       to the worker thread tree data structures, this is kind of a broadcast 
       to the workers (just kind of) 
    */

    masterBarrier(THREAD_COPY_BRANCHINFO_POINTER, tr);
    
    /* now the follwoing function does a couple of things:
       1. It traverses the tree and links the branchInfo array 
       entries to the individual branches
       2. It assigns memory space for the likelihood and scaling vectors 
       associated to the branch info data structure 
       3. It computes the STRIDED likelihood vectors to the left and right of 
       the branches and stores them contiguously in the vectors associated 
       to the branch info data structure
    */

    setupBranchMetaInfo(tr, &countBranches, q, tr->ntips, bInf, contiguousVectorLength,  contiguousScalingLength);

    /* just print out how many contiguuous vectors were assigned */
    printf("VectCount %d\n", vectCount);

    
    /* 
       check that the correct number of vectors was assigned.
       our reference tree has tr->ntips taxa 
       hence there exists a total of tr->ntips branches leading to leaves 
       and tr->ntips - 3 branches connecting inner nodes.
       For branches that lead to tips we will just need to allocate and store 
       one likelihood vector, whereas for inner branches we need two
    */
       


    assert(vectCount == (tr->ntips + ((tr->ntips - 3) * 2)));

    /* 
       just a little routine to test that what we have implemented is correct.
       We will traverse all branches of the tree and just compute the likelihood score.
       If we have two threads, thread 0 (the master thread) will compute even branches 
       and thread 1 will compute odd branches, evidently the evaluation of the likelihood
       at every branch should yield the same score 
    */
    
    masterBarrier(THREAD_TEST_CLASSIFY, tr);   
  }

#else
  setupBranchMetaInfo(tr, &countBranches, q, tr->ntips, bInf, 0, 0);
#endif

  

  assert(countBranches == numberOfBranches); 

  if(adef->printLabelledTree)    
    branchLabels(tr);    

  t = gettime();

  if(adef->rapidBoot)
    {                 
      double  
	*diagptable = (double *)NULL, 
	at;
      int 
	replicates, 
	*originalRateCategories = (int*)malloc(tr->cdta->endsite * sizeof(int)),
	*originalInvariant      = (int*)malloc(tr->cdta->endsite * sizeof(int)); 
      boolean oldThorough = thorough;

      if(align)	
	{
	  assert(tr->NumberOfModels == 1);
	  diagptable = makeDiagPtable(defaultz, dataType, tr->partitionData[0].gammaRates,tr->partitionData[0].EIGN);	 
	}	  	
      
      assert(adef->multipleRuns > 0);

      memcpy(originalRateCategories, tr->cdta->rateCategory, sizeof(int) * tr->cdta->endsite);
      memcpy(originalInvariant,      tr->invariant,          sizeof(int) * tr->cdta->endsite);

      /* need to be very careful about using compueNextReplicate for dynamic alignment */   

      at = gettime();       

      if(useHeuristics)
	heuristicInsertions = MAX(5, (int)(numberOfBranches/10));
      else
	heuristicInsertions = 0;

      if(heuristicInsertions > 0)
	{
	  bBr = (bestBranches *)malloc(sizeof(bestBranches) * numberOfTipsForInsertion);
	  for(i = 0; i < numberOfTipsForInsertion; i++)
	    {
	      bBr[i].likelihoods = (double *)    malloc(sizeof(double)      * heuristicInsertions);
	      bBr[i].branch      = (branchInfo**)malloc(sizeof(branchInfo*) * heuristicInsertions);
	      bBr[i].count       = 0;
	      bBr[i].minLH       = -1;
	      for(j = 0; j < heuristicInsertions; j++)
		{
		  bBr[i].likelihoods[j] = unlikely;
		  bBr[i].branch[j]      = (branchInfo*)NULL;
		}
	    }
	}

      /* compute the best insertion-based alignment on the original dataset */
      /* before doing BS replicates */
	
      evaluateGenericInitrav(tr, q->back);
      
      for(j = 0; j < numberOfTipsForInsertion; j++) 
	{
	  branches[j]          = -1.0;
	  insertLikelihoods[j] = unlikely;
	  insertNodes[j]     = (nodeptr)NULL;
	}  
      
      if(!thorough && useHeuristics)
	thorough = TRUE;

     

      if(align)
	{
	  addTraverseRob(tr, r, q, inserts,  numberOfTipsForInsertion, insertLikelihoods, (double*)NULL,
			 insertNodes, align, diagptable, writeAlignment, rawSeqs, thorough, dataType, bBr, 
			 heuristicInsertions, FALSE);
	  
	  /* now copy the aligned sequences into the buf from which we draw replicates, since 
	     tr->rdta->y0 will be over-written by computeNextReplicate */
	  
	  memcpy(tr->rdta->yBUF, tr->rdta->y0, tr->originalCrunchedLength * tr->rdta->numsp * sizeof(char));
	  
	  printBothOpen("Simultaneous query sequence alignment and insertion on original (non-bootstrapped) alignment and tree finished after %f seconds\n\n", gettime() - masterTime);
	}
      else
	{
	  addTraverseRob(tr, r, q, inserts,  numberOfTipsForInsertion, insertLikelihoods, (double*)NULL,
			 insertNodes, FALSE, (double*)NULL, writeAlignment, (unalignedSequences*)NULL, 
			 thorough, dataType, bBr, heuristicInsertions, FALSE);	    
	  
	  printBothOpen("Static query sequence insertion into original (non-bootstrapped) alignment and tree finished after %f seconds\n\n", gettime() - masterTime);
	}
      thorough = oldThorough;

      if(heuristicInsertions > 0)
	{
	  for(i = 0; i < numberOfTipsForInsertion; i++)	
	    {
	      assert(bBr[i].count == heuristicInsertions);
	      for(j = 0; j < bBr[i].count; j++)				 
		bBr[i].branch[j]->executeThem[i] = 1;		  		 		
	    }
	  
	  for(i = 0; i < numberOfTipsForInsertion; i++)
	    {
	      free(bBr[i].likelihoods);
	      free(bBr[i].branch);
	    }
	  free(bBr);		 
	}	  	      

      if(heuristicInsertions > 0)
	heuristicInsertions = -1;
      else
	heuristicInsertions = 0;     

      for(replicates = 0; replicates < adef->multipleRuns; replicates++)
	{
	  double rt = gettime();

	  computeNextReplicate(tr, &adef->rapidBoot, originalRateCategories, originalInvariant, TRUE);
	  /*computeNextReplicate(tr, adef, originalRateCategories, originalInvariant);*/
	  
	  resetBranches(tr);
	  evaluateGenericInitrav(tr, q->back);	  
	  treeEvaluate(tr, 1);
#ifdef _DEBUG_CLASSIFY	  
	  printf("Replicate %d %f\n", replicates, tr->likelihood);
#endif

	  for(j = 0; j < numberOfTipsForInsertion; j++) 
	    {
	      if(thorough)
		branches[j]          = -1.0;
	      insertLikelihoods[j] = unlikely;
	      insertNodes[j]     = (nodeptr)NULL;
	    }        

	  addTraverseRob(tr, r, q, inserts,  numberOfTipsForInsertion, insertLikelihoods, branches,
			 insertNodes, FALSE /*align we set align to FALSE here*/, 
			 (double*)NULL, FALSE, (unalignedSequences*)NULL /*rawSeqs*/, 
			 thorough, dataType, bBr, heuristicInsertions, TRUE);   
	  
	  
	  for(j = 0; j < numberOfTipsForInsertion; j++)
	    {
	      assert(insertNodes[j]);
	      insertNodes[j]->bInf->countThem[j] = insertNodes[j]->bInf->countThem[j] + 1;
	      if(thorough)
		{
		  insertNodes[j]->bInf->branches[j]  = insertNodes[j]->bInf->branches[j] + branches[j];
		  assert(branches[j] != -1.0);
		}
	    }

	  printBothOpen("Classification Replicate [%d] completed within %f seconds\n", replicates, gettime() - rt);	
	}
      printBothOpen("\n\n");
    }
  else    
    { 
      double  
	*diagptable = (double *)NULL;

      /* can't compute heuristics here, only for bootstrapping */

      assert(heuristicInsertions == 0);

      evaluateGenericInitrav(tr, q->back);                   
  
      for(j = 0; j < numberOfTipsForInsertion; j++) 
	{
	  if(thorough)
	    branches[j]          = -1.0;
	  insertLikelihoods[j] = unlikely;
	  insertNodes[j]     = (nodeptr)NULL;
	}                 

     

      /* 
	 TODO need to implement mixed models etc for dyn alignment at some point
	 we first need to understand dyn alignment though!	 
       */
      
      if(align)
	{
	  assert(tr->NumberOfModels == 1);
	  diagptable = makeDiagPtable(defaultz, dataType, tr->partitionData[0].gammaRates,tr->partitionData[0].EIGN);
	}
	  
      if(align)
	addTraverseRob(tr, r, q, inserts,  numberOfTipsForInsertion, insertLikelihoods, branches,
		       insertNodes, align, diagptable, writeAlignment, rawSeqs, thorough, dataType, bBr, heuristicInsertions, TRUE);       
      else
	addTraverseRob(tr, r, q, inserts,  numberOfTipsForInsertion, insertLikelihoods, branches,
		       insertNodes, align, (double*)NULL, writeAlignment, rawSeqs, thorough, dataType, bBr, heuristicInsertions, TRUE); 

      for(j = 0; j < numberOfTipsForInsertion; j++)
	{
	  assert(insertNodes[j]);
	  insertNodes[j]->bInf->countThem[j] = insertNodes[j]->bInf->countThem[j] + 1;
	  if(thorough)
	    {
	      insertNodes[j]->bInf->branches[j]  = insertNodes[j]->bInf->branches[j] + branches[j];
	      assert(branches[j] != -1);
	    }
	}

      free(diagptable); 
    }

  printBothOpen("Overall Classification time: %f\n\n", gettime() - masterTime);

#ifdef _DEBUG_CLASSIFY  
  if(!adef->rapidBoot)
    {
      evaluateGenericInitrav(tr, q->back);
      assert(startLH == tr->likelihood);
    }
#endif
  

  /* TODO need to actually think about this, a seg fault waiting to happen here */

  free(tr->tree_string);
  tr->treeStringLength *= 16;
  tr->tree_string  = (char*)malloc(tr->treeStringLength * sizeof(char));

  
 
  treeFile = myfopen(labelledTreeFileName, "w");
  Tree2StringClassify(tr->tree_string, tr, numberOfTipsForInsertion, inserts, FALSE, -1, -1, (branchInfo*)NULL, -1, FALSE);
  fprintf(treeFile, "%s\n", tr->tree_string);    
  fclose(treeFile);
    
  treeFile = myfopen(originalLabelledTreeFileName, "w");
  Tree2StringClassify(tr->tree_string, tr, numberOfTipsForInsertion, inserts, FALSE, -1, -1, (branchInfo*)NULL, -1, TRUE);
  fprintf(treeFile, "%s\n", tr->tree_string);    
  fclose(treeFile);


  if(adef->rapidBoot)
    {
      root = findRoot(tr->start,  numberOfTipsForInsertion, tr->mxtips);           
      
      if(root == 0)
	root = findRoot(tr->start->back,  numberOfTipsForInsertion, tr->mxtips);
      
      if(root == 0)    
	printf("We are in deep shit, the tree can't be rooted\n");	      
      else
	{
	  char extendedTreeFileName[1024];
	  int cutoff, k;
	  
#ifdef _DEBUG_CLASSIFY
	  printf("Can be rooted at %s\n", tr->nameList[root]);
#endif      
	  
	  for(k = 20; k < 100; k += 5)
	    {
	      cutoff = (int)(((double)adef->multipleRuns * (double)k / 100.0) + 0.5);
	      sprintf(extendedTreeFileName, "%s.%d", labelledTreeFileName, k);
	      
#ifdef _DEBUG_CLASSIFY  
	      printf("CUTOFF %d %d\n", cutoff, k);
#endif
	      
	      treeFile = myfopen(extendedTreeFileName, "w");	  
	      Tree2StringClassify(tr->tree_string, tr, numberOfTipsForInsertion, inserts, TRUE, cutoff, numberOfBranches, bInf, root, FALSE);	     
	      fprintf(treeFile, "%s\n", tr->tree_string);
	      fclose(treeFile);

	      printBothOpen("Least common ancestor file for cutoff %d (%d) written to file %s\n", cutoff, k, extendedTreeFileName);	      
	    }
	}
    }

 

  /* TODO also store and print out BS insertion br-lens
     at some point, only makes sense for thorough */

  classificationFile = myfopen(classificationFileName, "w");
  
  for(i = 0; i < numberOfTipsForInsertion; i++)    
    for(j = 0; j < numberOfBranches; j++) 
      {       
	if(bInf[j].countThem[i] > 0)	    
	  {
	    if(thorough)
	      fprintf(classificationFile, "%s I%d %d %8.20f\n", tr->nameList[inserts[i]], j, bInf[j].countThem[i], bInf[j].branches[i] / (double)bInf[j].countThem[i]);
	    else
	      fprintf(classificationFile, "%s I%d %d\n", tr->nameList[inserts[i]], j, bInf[j].countThem[i]);
	  }
      }

  fclose(classificationFile);

  if(writeAlignment)
    {
      assert(adef->dynamicAlignment);

      alignmentFile = myfopen(alignmentFileName, "w");

      /* TODO check what happpens with all these length indicators and data structures
	 after Bootstrapping */
      /*assert(tr->cdta->endsite == tr->originalCrunchedLength);*/
            
      /* restore alignment stored in y0, that was destroyed by 
	 the generation of BS replicates */
      
      if(adef->rapidBoot)
	memcpy(tr->rdta->y0, tr->rdta->yBUF, tr->originalCrunchedLength * tr->rdta->numsp * sizeof(char));

      fprintf(alignmentFile, "%d %d\n", tr->mxtips, tr->originalCrunchedLength);

      for(i = 1; i <= tr->mxtips; i++)
	{
	  unsigned char *tipI = tr->yVector[i];

	  fprintf(alignmentFile, "%s ", tr->nameList[i]); 

	  for(j = 0; j < tr->originalCrunchedLength; j++)
	    {
	       switch(dataType)
		 {
		 case AA_DATA:
		   fprintf(alignmentFile, "%c", inverseMeaningPROT[tipI[j]]);
		   break;
		 case DNA_DATA:
		   fprintf(alignmentFile, "%c", inverseMeaningDNA[tipI[j]]);
		   break;
		 case BINARY_DATA:
		   fprintf(alignmentFile, "%c", inverseMeaningBINARY[tipI[j]]);
		   break;
		 case SECONDARY_DATA:
		 case SECONDARY_DATA_6:
		 case SECONDARY_DATA_7:
		   assert(0);
		 default:
		   assert(0);
		 }
	    }	    
	  fprintf(alignmentFile, "\n");
	}

      fclose(alignmentFile);
    }   

  printBothOpen("\n\nLabelled reference tree including branch labels and query sequences written to file: %s\n\n", labelledTreeFileName); 
  printBothOpen("Labelled reference tree with branch labels (without query sequences) written to file: %s\n\n", originalLabelledTreeFileName); 
  printBothOpen("Classification result file written to file: %s\n\n", classificationFileName);

  if(writeAlignment)    
    printBothOpen("Extended ML alignment written to file: %s\n\n", alignmentFileName);
    

  exit(0);
}
