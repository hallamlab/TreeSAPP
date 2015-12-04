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

extern int Thorough;

static int determineRearrangementSettingRapid(tree *tr,  analdef *adef, bestlist *bestT, bestlist *bt)
{
  int i, 
    maxtrav, 
    bestTrav, 
    impr, 
    index, 
    MaxFast,
    it = 0;
  double startLH; 
  boolean cutoff;  

  MaxFast = 26;

  startLH = tr->likelihood;
  
  tr->doCutoff  = FALSE;
  tr->bigCutoff = FALSE;
  tr->itCount   = 0;
     
  maxtrav = 5;

  bestTrav = maxtrav = 5;

  impr = 1;

  resetBestTree(bt); 

  while(impr && maxtrav < MaxFast)
    {	
      recallBestTree(bestT, 1, tr);     
      nodeRectifier(tr);
           
      if (maxtrav > tr->ntips - 3)  
	maxtrav = tr->ntips - 3;    
 
      tr->startLH = tr->endLH = tr->likelihood;
          
      printf("Cutoff: %d\n", tr->doCutoff);

      for(i = 1; i <= tr->mxtips + tr->mxtips - 2; i++)
	{                	 
	  index = i;	 	 

	  tr->bestOfNode = unlikely;
	  if(rearrangeBIG(tr, tr->nodep[index], 1, maxtrav))
	    {	     
	      if(tr->endLH > tr->startLH)                 	
		{		 	 	      
		  restoreTreeFast(tr);	        	  	 	  	      
		  tr->startLH = tr->endLH = tr->likelihood;	  	 	  	  	  	  	  	  	      
		}	         	       	
	    }
	}
      
      treeEvaluate(tr, 0.125);
      printf("maxtrav %d %f\n", maxtrav, tr->likelihood);
      saveBestTree(bt, tr);                                    
      
      /*if(it == 0)
	{
	  tr->doCutoff  = TRUE;
	  tr->bigCutoff = FALSE;	  
	  }*/
      it++;

      if(tr->likelihood > startLH)
	{	 
	  startLH = tr->likelihood; 	  	  	  	 
	  bestTrav = maxtrav;	 
	  impr = 1;
	}
      else
	{
	  impr = 0;
	}
      maxtrav += 5;            
    }

  recallBestTree(bt, 1, tr);   
  
  return bestTrav;     
}

void superFast(tree *tr, analdef *adef)
{
  double epsilon;

  bestlist *bestT, *bt;  
  
  bestT = (bestlist *) malloc(sizeof(bestlist));
  bestT->ninit = 0;
  initBestTree(bestT, 1, tr->mxtips);
      
  bt = (bestlist *) malloc(sizeof(bestlist));      
  bt->ninit = 0;
  initBestTree(bt, 20, tr->mxtips); 

  Thorough = 0;     

  initInfoList(50);

  tr->rateHetModel = CAT;
  adef->categories = 10;
  initModel(tr, tr->rdta, tr->cdta, adef); 

  printf("Rate Cats: %d\n",  tr->NumberOfCategories);

  makeParsimonyTreeRapid(tr, adef);

  onlyInitrav(tr, tr->start);
  treeEvaluate(tr, 1);

  epsilon = -0.001 * tr->likelihood;
  printf("Epsilon: %f\n", epsilon);

  modOpt(tr, adef, FALSE, epsilon);
  printf("Rate Cats: %d\n",  tr->NumberOfCategories);
  printf("%f\n", tr->likelihood);
  
  adef->bestTrav = determineRearrangementSettingRapid(tr, adef, bestT, bt);
  treeEvaluate(tr, 1);

  exit(0);
}
