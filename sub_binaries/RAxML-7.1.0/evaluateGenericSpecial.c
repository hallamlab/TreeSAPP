/*  RAxML-VI-HPC (version 2.2) a program for sequential and parallel estimation of phylogenetic trees 
 *  Copyright August 2006 by Alexandros Stamatakis
 *
 *  Partially derived from
 *  fastDNAml, a program for estimation of phylogenetic trees from sequences by Gary J. Olsen
 *  
 *  and 
 *
 *  Programs of the PHYLIP package by Joe Felsenstein. 
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


#ifdef _USE_PTHREADS
extern volatile double *reductionBuffer;
extern volatile int NumberOfThreads;
#endif

void calcDiagptable(double z, int data, int numberOfCategories, double *rptr, double *EIGN, double *diagptable)
{
  int i, l;
  double lz;

  if (z < zmin) 
    lz = log(zmin);
  else
    lz = log(z);

  switch(data)
    {
    case BINARY_DATA:
       {
	double lz1;
	lz1 = EIGN[0] * lz;
	for(i = 0; i <  numberOfCategories; i++)
	  {		 
	    diagptable[2 * i] = 1.0;
	    diagptable[2 * i + 1] = exp(rptr[i] * lz1);	   	    
	  }
      }
      break;
    case DNA_DATA:
      {
	double lz1, lz2, lz3;
	lz1 = EIGN[0] * lz;
	lz2 = EIGN[1] * lz;
	lz3 = EIGN[2] * lz;

	for(i = 0; i <  numberOfCategories; i++)
	  {		 
	    diagptable[4 * i] = 1.0;
	    diagptable[4 * i + 1] = exp(rptr[i] * lz1);
	    diagptable[4 * i + 2] = exp(rptr[i] * lz2);
	    diagptable[4 * i + 3] = exp(rptr[i] * lz3);
	    /*printf("%f %f %f %f\n",diagptable[4 * i],  diagptable[4 * i + 1], diagptable[4 * i + 2], diagptable[4 * i + 3]);*/
	  }
      }
      break;
    case AA_DATA:
      {
	double lza[19];

	for(l = 0; l < 19; l++)      
	  lza[l] = EIGN[l] * lz; 

	for(i = 0; i <  numberOfCategories; i++)
	  {	      	       
	    diagptable[i * 20] = 1.0;

	    for(l = 1; l < 20; l++)
	      diagptable[i * 20 + l] = exp(rptr[i] * lza[l - 1]);     	          
	  }
      }
      break;
    case SECONDARY_DATA:
      {
	double lza[15];

	for(l = 0; l < 15; l++)      
	  lza[l] = EIGN[l] * lz; 

	for(i = 0; i <  numberOfCategories; i++)
	  {	      	       
	    diagptable[i * 16] = 1.0;

	    for(l = 1; l < 16; l++)
	      diagptable[i * 16 + l] = exp(rptr[i] * lza[l - 1]);     	          
	  }
      }
      break;
    case SECONDARY_DATA_6:
      {
	double lza[5];

	for(l = 0; l < 5; l++)      
	  lza[l] = EIGN[l] * lz; 

	for(i = 0; i <  numberOfCategories; i++)
	  {	      	       
	    diagptable[i * 6] = 1.0;

	    for(l = 1; l < 6; l++)
	      diagptable[i * 6 + l] = exp(rptr[i] * lza[l - 1]);     	          
	  }
      }
      break;
    case SECONDARY_DATA_7:
      {
	double lza[6];

	for(l = 0; l < 6; l++)      
	  lza[l] = EIGN[l] * lz; 

	for(i = 0; i <  numberOfCategories; i++)
	  {	      	       
	    diagptable[i * 7] = 1.0;

	    for(l = 1; l < 7; l++)
	      diagptable[i * 7 + l] = exp(rptr[i] * lza[l - 1]);     	          
	  }
      }
      break;
    default:
      assert(0);
    }
}





static double evaluateGTRCATPROT (int *ex1, int *ex2, int *cptr, int *wptr,
				  double *x1, double *x2, double *tipVector,
				  unsigned char *tipX1, int n, double *diagptable_start, double *vector, boolean writeVector)
{
  double   sum = 0.0, term;
  double  *diagptable,  *left, *right;
  int     i, l;                           
  
  if(tipX1)
    {            
      if(writeVector)
	for (i = 0; i < n; i++) 
	  {	       	
	    left = &(tipVector[20 * tipX1[i]]);
	    right = &(x2[20 * i]);
	    
	    diagptable = &diagptable_start[20 * cptr[i]];	           	 
	    
	    for(l = 0, term = 0.0; l < 20; l++)
	      term += left[l] * right[l] * diagptable[l];	 	  	  
	    
	    term = log(term) + (ex2[i] * log(minlikelihood));
	   
	    vector[i] = term;
	    
	    sum += wptr[i] * term;
	  }         
      else
	for (i = 0; i < n; i++) 
	  {	       	
	    left = &(tipVector[20 * tipX1[i]]);
	    right = &(x2[20 * i]);
	    
	    diagptable = &diagptable_start[20 * cptr[i]];	           	 
	    
	    for(l = 0, term = 0.0; l < 20; l++)
	      term += left[l] * right[l] * diagptable[l];	 	  	  
	    
	    term = log(term) + (ex2[i] * log(minlikelihood));
	   
	    sum += wptr[i] * term;
	  }      
    }    
  else
    {
    
      for (i = 0; i < n; i++) 
	{		       	      	      
	  left  = &x1[20 * i];
	  right = &x2[20 * i];
	  
	  diagptable = &diagptable_start[20 * cptr[i]];	  	
	  
	  for(l = 0, term = 0.0; l < 20; l++)
	    term += left[l] * right[l] * diagptable[l];	
	  
	  term = log(term) + ((ex1[i] + ex2[i]) * log(minlikelihood));
	  
	  sum += wptr[i] * term;      
	}
    }
             
  return  sum;         
} 

static double evaluateGTRCATSECONDARY (int *ex1, int *ex2, int *cptr, int *wptr,
				       double *x1, double *x2, double *tipVector,
				       unsigned char *tipX1, int n, double *diagptable_start, double *vector, boolean writeVector)
{
  double   sum = 0.0, term;
  double  *diagptable,  *left, *right;
  int     i, l;                           
  
  if(tipX1)
    {            
      if(writeVector)
	for (i = 0; i < n; i++) 
	  {	       	
	    left = &(tipVector[16 * tipX1[i]]);
	    right = &(x2[16 * i]);
	    
	    diagptable = &diagptable_start[16 * cptr[i]];	           	 
	    
	    for(l = 0, term = 0.0; l < 16; l++)
	      term += left[l] * right[l] * diagptable[l];	 	  	  
	    
	    term = log(term) + (ex2[i] * log(minlikelihood));
	   
	    vector[i] = term;
	    
	    sum += wptr[i] * term;
	  }         
      else
	for (i = 0; i < n; i++) 
	  {	       	
	    left = &(tipVector[16 * tipX1[i]]);
	    right = &(x2[16 * i]);
	    
	    diagptable = &diagptable_start[16 * cptr[i]];	           	 
	    
	    for(l = 0, term = 0.0; l < 16; l++)
	      term += left[l] * right[l] * diagptable[l];	 	  	  
	    
	    term = log(term) + (ex2[i] * log(minlikelihood));
	   
	    sum += wptr[i] * term;
	  }      
    }    
  else
    {
    
      for (i = 0; i < n; i++) 
	{		       	      	      
	  left  = &x1[16 * i];
	  right = &x2[16 * i];
	  
	  diagptable = &diagptable_start[16 * cptr[i]];	  	
	  
	  for(l = 0, term = 0.0; l < 16; l++)
	    term += left[l] * right[l] * diagptable[l];	
	  
	  term = log(term) + ((ex1[i] + ex2[i]) * log(minlikelihood));
	  
	  sum += wptr[i] * term;      
	}
    }
             
  return  sum;         
} 

static double evaluateGTRCATSECONDARY_6 (int *ex1, int *ex2, int *cptr, int *wptr,
				       double *x1, double *x2, double *tipVector,
				       unsigned char *tipX1, int n, double *diagptable_start, double *vector, boolean writeVector)
{
  double   sum = 0.0, term;
  double  *diagptable,  *left, *right;
  int     i, l;                           
  
  if(tipX1)
    {            
      if(writeVector)
	for (i = 0; i < n; i++) 
	  {	       	
	    left = &(tipVector[6 * tipX1[i]]);
	    right = &(x2[6 * i]);
	    
	    diagptable = &diagptable_start[6 * cptr[i]];	           	 
	    
	    for(l = 0, term = 0.0; l < 6; l++)
	      term += left[l] * right[l] * diagptable[l];	 	  	  
	    
	    term = log(term) + (ex2[i] * log(minlikelihood));
	   
	    vector[i] = term;
	    
	    sum += wptr[i] * term;
	  }         
      else
	for (i = 0; i < n; i++) 
	  {	       	
	    left = &(tipVector[6 * tipX1[i]]);
	    right = &(x2[6 * i]);
	    
	    diagptable = &diagptable_start[6 * cptr[i]];	           	 
	    
	    for(l = 0, term = 0.0; l < 6; l++)
	      term += left[l] * right[l] * diagptable[l];	 	  	  
	    
	    term = log(term) + (ex2[i] * log(minlikelihood));
	   
	    sum += wptr[i] * term;
	  }      
    }    
  else
    {
    
      for (i = 0; i < n; i++) 
	{		       	      	      
	  left  = &x1[6 * i];
	  right = &x2[6 * i];
	  
	  diagptable = &diagptable_start[6 * cptr[i]];	  	
	  
	  for(l = 0, term = 0.0; l < 6; l++)
	    term += left[l] * right[l] * diagptable[l];	
	  
	  term = log(term) + ((ex1[i] + ex2[i]) * log(minlikelihood));
	  
	  sum += wptr[i] * term;      
	}
    }
             
  return  sum;         
} 

static double evaluateGTRCATSECONDARY_7(int *ex1, int *ex2, int *cptr, int *wptr,
					double *x1, double *x2, double *tipVector,
					unsigned char *tipX1, int n, double *diagptable_start, double *vector, boolean writeVector)
{
  double   sum = 0.0, term;
  double  *diagptable,  *left, *right;
  int     i, l;                           
  
  if(tipX1)
    {            
      if(writeVector)
	for (i = 0; i < n; i++) 
	  {	       	
	    left = &(tipVector[7 * tipX1[i]]);
	    right = &(x2[7 * i]);
	    
	    diagptable = &diagptable_start[7 * cptr[i]];	           	 
	    
	    for(l = 0, term = 0.0; l < 7; l++)
	      term += left[l] * right[l] * diagptable[l];	 	  	  
	    
	    term = log(term) + (ex2[i] * log(minlikelihood));
	   
	    vector[i] = term;
	    
	    sum += wptr[i] * term;
	  }         
      else
	for (i = 0; i < n; i++) 
	  {	       	
	    left = &(tipVector[7 * tipX1[i]]);
	    right = &(x2[7 * i]);
	    
	    diagptable = &diagptable_start[7 * cptr[i]];	           	 
	    
	    for(l = 0, term = 0.0; l < 7; l++)
	      term += left[l] * right[l] * diagptable[l];	 	  	  
	    
	    term = log(term) + (ex2[i] * log(minlikelihood));
	   
	    sum += wptr[i] * term;
	  }      
    }    
  else
    {
    
      for (i = 0; i < n; i++) 
	{		       	      	      
	  left  = &x1[7 * i];
	  right = &x2[7 * i];
	  
	  diagptable = &diagptable_start[7 * cptr[i]];	  	
	  
	  for(l = 0, term = 0.0; l < 7; l++)
	    term += left[l] * right[l] * diagptable[l];	
	  
	  term = log(term) + ((ex1[i] + ex2[i]) * log(minlikelihood));
	  
	  sum += wptr[i] * term;      
	}
    }
             
  return  sum;         
} 

static double evaluateGTRCAT_BINARY (int *ex1, int *ex2, int *cptr, int *wptr,
				     double *x1_start, double *x2_start, double *tipVector, 		      
				     unsigned char *tipX1, int n, double *diagptable_start, double *vector, boolean writeVector)
{
  double  sum = 0.0, term;       
  int     i, j;  
  double  *diagptable, *x1, *x2;                      	    
 
  if(tipX1)
    {     
      if(writeVector)
	for (i = 0; i < n; i++) 
	  {	    		   	  
	    x1 = &(tipVector[2 * tipX1[i]]);
	    x2 = &x2_start[2 * i];
	    
	    diagptable = &diagptable_start[2 * cptr[i]];	    	    	  
	    
	    for(j = 0, term = 0.0; j < 2; j++)
	      term += x1[j] * x2[j] * diagptable[j];
	     
	    term = log(term) + (ex2[i] * log(minlikelihood));	   
	    
	    vector[i] = term;	   
	    
	    sum += wptr[i] * term;
	  }	
      else
	for (i = 0; i < n; i++) 
	  {	    		   	  
	    x1 = &(tipVector[2 * tipX1[i]]);
	    x2 = &x2_start[2 * i];
	    
	    diagptable = &diagptable_start[2 * cptr[i]];	    	    	  
	    
	    for(j = 0, term = 0.0; j < 2; j++)
	      {
		/*printf("%1.40f %1.40f\n", x1[j], x2[j]);*/
		term += x1[j] * x2[j] * diagptable[j];
	      }
	    
	    /*printf("GAG %d %f\n", i, term);*/
	    term = log(term) + (ex2[i] * log(minlikelihood));	   	    	   	 	  	  	 
	    
	    sum += wptr[i] * term;
	  }	
    }               
  else
    {
      for (i = 0; i < n; i++) 
	{		          	
	  x1 = &x1_start[2 * i];
	  x2 = &x2_start[2 * i];
	  
	  diagptable = &diagptable_start[2 * cptr[i]];		  
	  
	  for(j = 0, term = 0.0; j < 2; j++)
	    term += x1[j] * x2[j] * diagptable[j];     
	  /*printf("%d %f\n", i, term);*/
	  term = log(term) + ((ex1[i] + ex2[i]) * log(minlikelihood));
	  
	  sum += wptr[i] * term;
	}	   
    }
       
  return  sum;         
} 


static double evaluateGTRGAMMA_BINARY(int *ex1, int *ex2, int *wptr,
				      double *x1_start, double *x2_start, 
				      double *tipVector, 
				      unsigned char *tipX1, const int n, double *diagptable, double *vector, boolean writeVector)
{
  double   sum = 0.0, term;    
  int     i, j, k;
  double  *x1, *x2;             

  /* C-OPT once again switch to figure out if the node to the left or right of the branch at which we compute
     the final likelihood is a tip or not */

  if(tipX1)
    {    
      if(writeVector)      
	for (i = 0; i < n; i++)
	  {
	    x1 = &(tipVector[2 * tipX1[i]]);	 
	    x2 = &x2_start[8 * i];	          	  	
	    
	    for(j = 0, term = 0.0; j < 4; j++)
	      for(k = 0; k < 2; k++)
		term += x1[k] * x2[j * 2 + k] * diagptable[j * 2 + k];	          	  	  
	    
	    term = log(0.25 * term) + ex2[i] * log(minlikelihood);	 
	    	    
	    vector[i] = term;
	    
	    sum += wptr[i] * term;
	  }	  
      else	
	/* C-OPT loop over the input data length, this is where the interesting things happen */
	for (i = 0; i < n; i++)
	  {
	    x1 = &(tipVector[2 * tipX1[i]]);	 
	    x2 = &x2_start[8 * i];	          	  	
	    
	    for(j = 0, term = 0.0; j < 4; j++)
	      for(k = 0; k < 2; k++)
		term += x1[k] * x2[j * 2 + k] * diagptable[j * 2 + k];	          	  	  	    	    

	    term = log(0.25 * term) + ex2[i] * log(minlikelihood);	 

	    sum += wptr[i] * term;
	  }	  
    }
  else
    {   
      /* C-OPT loop over the input data length, this is where the interesting things happen */
      for (i = 0; i < n; i++) 
	{	  	 	  	  
	  x1 = &x1_start[8 * i];
	  x2 = &x2_start[8 * i];	  	  
	  
	  for(j = 0, term = 0.0; j < 4; j++)
	    for(k = 0; k < 2; k++)
	      term += x1[j * 2 + k] * x2[j * 2 + k] * diagptable[j * 2 + k];	          	  	  	      
	  
	  term = log(0.25 * term) + (ex1[i] + ex2[i]) * log(minlikelihood);

	  sum += wptr[i] * term;
	}                      	
    }

  return sum;
} 

static double evaluateGTRGAMMAINVAR_BINARY (int *ex1, int *ex2, int *wptr, int *iptr,
					    double *x1_start, double *x2_start,
					    double *tipVector, double *tFreqs, double invariants,
					    unsigned char *tipX1, int n, double *diagptable, double *vector, boolean writeVector)
{ 
  int     i, j, k;
  double  *x1, *x2; 
  double 
    freqs[2], 
    scaler = 0.25 * (1.0 - invariants),
    sum = 0.0, 
    term; 

  freqs[0] = tFreqs[0] * invariants; 
  freqs[1] = tFreqs[1] * invariants; 

  if(tipX1)
    {   
      if(writeVector)
	for (i = 0; i < n; i++) 
	  {
	    x1 = &(tipVector[2 * tipX1[i]]);
	    x2 = &x2_start[8 * i];	  
	    
	    for(j = 0, term = 0.0; j < 4; j++)
	      for(k = 0; k < 2; k++)
		term += x1[k] * x2[j * 2 + k] * diagptable[j * 2 + k];
	    
	    if(iptr[i] < 2)	   
	      term = log(((scaler * term) + freqs[iptr[i]]) * pow(minlikelihood, ex2[i]));
	    else
	      term = log(scaler * term) + (ex2[i] * log(minlikelihood));	 
	    	    
	    vector[i] = term;
	    
	    sum += wptr[i] * term;
	  }	  
      else
	for (i = 0; i < n; i++) 
	  {
	    x1 = &(tipVector[2 * tipX1[i]]);
	    x2 = &x2_start[8 * i];	  
	    
	    for(j = 0, term = 0.0; j < 4; j++)
	      for(k = 0; k < 2; k++)
		term += x1[k] * x2[j * 2 + k] * diagptable[j * 2 + k];
	    
	    if(iptr[i] < 2)	   
	      term = log(((scaler * term) + freqs[iptr[i]]) * pow(minlikelihood, ex2[i]));
	    else
	      term = log(scaler * term) + (ex2[i] * log(minlikelihood));	 
	    	   
	    sum += wptr[i] * term;
	  }	  
    }
  else
    {           		

      for (i = 0; i < n; i++) 
	{	  	 	  	
	  x1 = &x1_start[8 * i];
	  x2 = &x2_start[8 * i];	  
	  
	  for(j = 0, term = 0.0; j < 4; j++)
	    for(k = 0; k < 2; k++)
	      term += x1[j * 2 + k] * x2[j * 2 + k] * diagptable[j * 2 + k];	  	 	      
	  
	  if(iptr[i] < 2)	   
	    term = log(((scaler * term) + freqs[iptr[i]]) * pow(minlikelihood, ex2[i] + ex1[i]));
	  else
	    term = log(scaler * term) + ((ex1[i] + ex2[i]) * log(minlikelihood));

	  sum += wptr[i] * term;
	}	  	                        
    }

  return  sum;
} 


static double evaluateGTRCAT (int *ex1, int *ex2, int *cptr, int *wptr,
			      double *x1_start, double *x2_start, double *tipVector, 		      
			      unsigned char *tipX1, int n, double *diagptable_start, double *vector, boolean writeVector)
{
  double  sum = 0.0, term;       
  int     i, j;  
  double  *diagptable, *x1, *x2;                      	    
 
  if(tipX1)
    {     
      if(writeVector)
	for (i = 0; i < n; i++) 
	  {	    		   	  
	    x1 = &(tipVector[4 * tipX1[i]]);
	    x2 = &x2_start[4 * i];
	    
	    diagptable = &diagptable_start[4 * cptr[i]];	    	    	  
	    
	    for(j = 0, term = 0.0; j < 4; j++)
	      term += x1[j] * x2[j] * diagptable[j];
	     
	    term = log(term) + (ex2[i] * log(minlikelihood));	   
	    
	    vector[i] = term;	   
	    
	    sum += wptr[i] * term;
	  }	
      else
	for (i = 0; i < n; i++) 
	  {		   		   
	    x1 = &(tipVector[4 * tipX1[i]]);
	    x2 = &x2_start[4 * i];
	    
	    diagptable = &diagptable_start[4 * cptr[i]];	    	    	  
	    
	    for(j = 0, term = 0.0; j < 4; j++)
	      term += x1[j] * x2[j] * diagptable[j];
	    	    
	    term = log(term) + (ex2[i] * log(minlikelihood));	   	    	   	 	  	  	 
	    
	    sum += wptr[i] * term;
	  }	
    }               
  else
    {
      for (i = 0; i < n; i++) 
	{	 	           	
	  x1 = &x1_start[4 * i];
	  x2 = &x2_start[4 * i];
	  
	  diagptable = &diagptable_start[4 * cptr[i]];		  
	  
	  for(j = 0, term = 0.0; j < 4; j++)
	    term += x1[j] * x2[j] * diagptable[j];     
	  
	  term = log(term) + ((ex1[i] + ex2[i]) * log(minlikelihood));	  
	  sum += wptr[i] * term;
	}    
    }
       
  return  sum;         
} 


/* C-OPT the function below computes the actual likelihood score and will 
   take 5% of total execution time, writeVector will always be false 
   in our case */


static double evaluateGTRGAMMA(int *ex1, int *ex2, int *wptr,
			       double *x1_start, double *x2_start, 
			       double *tipVector, 
			       unsigned char *tipX1, const int n, double *diagptable, double *vector, boolean writeVector)
{
  double   sum = 0.0, term;    
  int     i, j, k;
  double  *x1, *x2;             

  /* C-OPT once again switch to figure out if the node to the left or right of the branch at which we compute
     the final likelihood is a tip or not */

  if(tipX1)
    {    
      if(writeVector)      
	for (i = 0; i < n; i++)
	  {
	    x1 = &(tipVector[4 * tipX1[i]]);	 
	    x2 = &x2_start[16 * i];	          	  	
	    
	    for(j = 0, term = 0.0; j < 4; j++)
	      for(k = 0; k < 4; k++)
		term += x1[k] * x2[j * 4 + k] * diagptable[j * 4 + k];	          	  	  
	    
	    term = log(0.25 * term) + ex2[i] * log(minlikelihood);	 
	    	    
	    vector[i] = term;
	    
	    sum += wptr[i] * term;
	  }	  
      else	
	/* C-OPT loop over the input data length, this is where the interesting things happen */
	for (i = 0; i < n; i++)
	  {
	    x1 = &(tipVector[4 * tipX1[i]]);	 
	    x2 = &x2_start[16 * i];	          	  	
	    
	    for(j = 0, term = 0.0; j < 4; j++)
	      for(k = 0; k < 4; k++)
		term += x1[k] * x2[j * 4 + k] * diagptable[j * 4 + k];	          	  	  	    	    

	    term = log(0.25 * term) + ex2[i] * log(minlikelihood);	 

	    sum += wptr[i] * term;
	  }	  
    }
  else
    {   
      /* C-OPT loop over the input data length, this is where the interesting things happen */
      for (i = 0; i < n; i++) 
	{	  	 	  	  
	  x1 = &x1_start[16 * i];
	  x2 = &x2_start[16 * i];	  	  
	  
	  for(j = 0, term = 0.0; j < 4; j++)
	    for(k = 0; k < 4; k++)
	      term += x1[j * 4 + k] * x2[j * 4 + k] * diagptable[j * 4 + k];	          	  	  	      
	  
	  term = log(0.25 * term) + (ex1[i] + ex2[i]) * log(minlikelihood);

	  sum += wptr[i] * term;
	}                      	
    }

  return sum;
} 






static double evaluateGTRGAMMAINVAR (int *ex1, int *ex2, int *wptr, int *iptr,
				     double *x1_start, double *x2_start,
				     double *tipVector, double *tFreqs, double invariants,
				     unsigned char *tipX1, int n, double *diagptable, double *vector, boolean writeVector)
{ 
  int     i, j, k;
  double  *x1, *x2; 
  double 
    freqs[4], 
    scaler = 0.25 * (1.0 - invariants),
    sum = 0.0, 
    term; 

  freqs[0] = tFreqs[0] * invariants; 
  freqs[1] = tFreqs[1] * invariants;
  freqs[2] = tFreqs[2] * invariants;
  freqs[3] = tFreqs[3] * invariants;   

  if(tipX1)
    {   
      if(writeVector)
	for (i = 0; i < n; i++) 
	  {
	    x1 = &(tipVector[4 * tipX1[i]]);
	    x2 = &x2_start[16 * i];	  
	    
	    for(j = 0, term = 0.0; j < 4; j++)
	      for(k = 0; k < 4; k++)
		term += x1[k] * x2[j * 4 + k] * diagptable[j * 4 + k];
	    
	    if(iptr[i] < 4)	   
	      term = log(((scaler * term) + freqs[iptr[i]]) * pow(minlikelihood, ex2[i]));
	    else
	      term = log(scaler * term) + (ex2[i] * log(minlikelihood));	 
	    	    
	    vector[i] = term;
	    
	    sum += wptr[i] * term;
	  }	  
      else
	for (i = 0; i < n; i++) 
	  {
	    x1 = &(tipVector[4 * tipX1[i]]);
	    x2 = &x2_start[16 * i];	  
	    
	    for(j = 0, term = 0.0; j < 4; j++)
	      for(k = 0; k < 4; k++)
		term += x1[k] * x2[j * 4 + k] * diagptable[j * 4 + k];
	    
	    if(iptr[i] < 4)	   
	      term = log(((scaler * term) + freqs[iptr[i]]) * pow(minlikelihood, ex2[i]));
	    else
	      term = log(scaler * term) + (ex2[i] * log(minlikelihood));	 
	    	   
	    sum += wptr[i] * term;
	  }	  
    }
  else
    {           		

      for (i = 0; i < n; i++) 
	{	  	 	  	
	  x1 = &x1_start[16 * i];
	  x2 = &x2_start[16 * i];	  
	  
	  for(j = 0, term = 0.0; j < 4; j++)
	    for(k = 0; k < 4; k++)
	      term += x1[j * 4 + k] * x2[j * 4 + k] * diagptable[j * 4 + k];	  	 	      
	  
	  if(iptr[i] < 4)	   
	    term = log(((scaler * term) + freqs[iptr[i]]) * pow(minlikelihood, ex2[i] + ex1[i]));
	  else
	    term = log(scaler * term) + ((ex1[i] + ex2[i]) * log(minlikelihood));

	  sum += wptr[i] * term;
	}	  	                        
    }

  return  sum;
} 




static double evaluateGTRGAMMAPROT (int *ex1, int *ex2, int *wptr,
				    double *x1, double *x2,  
				    double *tipVector, 
				    unsigned char *tipX1, int n, double *diagptable, double *vector, boolean writeVector)
{
  double   sum = 0.0, term;        
  int     i, j, l;   
  double  *left, *right;              
  
  if(tipX1)
    {          
      if(writeVector)
	for (i = 0; i < n; i++) 
	  {
	    left = &(tipVector[20 * tipX1[i]]);	  	  
	    
	    for(j = 0, term = 0.0; j < 4; j++)
	      {
		right = &(x2[80 * i + 20 * j]);
		
	      for(l = 0; l < 20; l++)
		term += left[l] * right[l] * diagptable[j * 20 + l];	      
	      }	  
	    
	    term = log(0.25 * term) + (ex2[i] * log(minlikelihood));	   
	    	    
	    vector[i] = term;
	    
	    sum += wptr[i] * term;
	  }         
      else
	{       
	  for (i = 0; i < n; i++) 
	    {	     
	      left = &(tipVector[20 * tipX1[i]]);	  	  
	      
	      for(j = 0, term = 0.0; j < 4; j++)
		{
		  right = &(x2[80 * i + 20 * j]);
		  
		  for(l = 0; l < 20; l++)
		    term += left[l] * right[l] * diagptable[j * 20 + l];	      
		}	  
	      
	      term = log(0.25 * term) + (ex2[i] * log(minlikelihood));	   
	      
	      sum += wptr[i] * term;
	    }     	 
	}
    }              
  else
    {
      for (i = 0; i < n; i++) 
	{	  	 	             
      
	  for(j = 0, term = 0.0; j < 4; j++)
	    {
	      left  = &(x1[80 * i + 20 * j]);
	      right = &(x2[80 * i + 20 * j]);	    
	      
	      for(l = 0; l < 20; l++)
		term += left[l] * right[l] * diagptable[j * 20 + l];	
	    }
	  
	  term = log(0.25 * term) + ((ex1[i] + ex2[i])*log(minlikelihood));
	  
	  sum += wptr[i] * term;
	}         
    }
       
  return  sum;
}

static double evaluateGTRGAMMASECONDARY (int *ex1, int *ex2, int *wptr,
					 double *x1, double *x2,  
					 double *tipVector, 
					 unsigned char *tipX1, int n, double *diagptable, double *vector, boolean writeVector)
{
  double   sum = 0.0, term;        
  int     i, j, l;   
  double  *left, *right;              
  
  if(tipX1)
    {          
      if(writeVector)
	for (i = 0; i < n; i++) 
	  {
	    left = &(tipVector[16 * tipX1[i]]);	  	  
	    
	    for(j = 0, term = 0.0; j < 4; j++)
	      {
		right = &(x2[64 * i + 16 * j]);
		
		for(l = 0; l < 16; l++)
		  term += left[l] * right[l] * diagptable[j * 16 + l];	      
	      }	  
	    
	    term = log(0.25 * term) + (ex2[i] * log(minlikelihood));	   
	    	    
	    vector[i] = term;
	    
	    sum += wptr[i] * term;
	  }         
      else
	{       
	  for (i = 0; i < n; i++) 
	    {	     
	      left = &(tipVector[16 * tipX1[i]]);	  	  
	      
	      for(j = 0, term = 0.0; j < 4; j++)
		{
		  right = &(x2[64 * i + 16 * j]);
		  
		  for(l = 0; l < 16; l++)
		    term += left[l] * right[l] * diagptable[j * 16 + l];	      
		}	  
	      
	      term = log(0.25 * term) + (ex2[i] * log(minlikelihood));	   
	      
	      sum += wptr[i] * term;
	    }     	 
	}
    }              
  else
    {
      for (i = 0; i < n; i++) 
	{	  	 	             
      
	  for(j = 0, term = 0.0; j < 4; j++)
	    {
	      left  = &(x1[64 * i + 16 * j]);
	      right = &(x2[64 * i + 16 * j]);	    
	      
	      for(l = 0; l < 16; l++)
		term += left[l] * right[l] * diagptable[j * 16 + l];	
	    }
	  
	  term = log(0.25 * term) + ((ex1[i] + ex2[i])*log(minlikelihood));
	  
	  sum += wptr[i] * term;
	}         
    }
       
  return  sum;
}

static double evaluateGTRGAMMASECONDARY_6 (int *ex1, int *ex2, int *wptr,
					   double *x1, double *x2,  
					   double *tipVector, 
					   unsigned char *tipX1, int n, double *diagptable, double *vector, boolean writeVector)
{
  double   sum = 0.0, term;        
  int     i, j, l;   
  double  *left, *right;              
  
  if(tipX1)
    {          
      if(writeVector)
	for (i = 0; i < n; i++) 
	  {
	    left = &(tipVector[6 * tipX1[i]]);	  	  
	    
	    for(j = 0, term = 0.0; j < 4; j++)
	      {
		right = &(x2[24 * i + 6 * j]);
		
		for(l = 0; l < 6; l++)
		  term += left[l] * right[l] * diagptable[j * 6 + l];	      
	      }	  
	    
	    term = log(0.25 * term) + (ex2[i] * log(minlikelihood));	   
	    	    
	    vector[i] = term;
	    
	    sum += wptr[i] * term;
	  }         
      else
	{       
	  for (i = 0; i < n; i++) 
	    {	     
	      left = &(tipVector[6 * tipX1[i]]);	  	  
	      
	      for(j = 0, term = 0.0; j < 4; j++)
		{
		  right = &(x2[24 * i + 6 * j]);
		  
		  for(l = 0; l < 6; l++)
		    term += left[l] * right[l] * diagptable[j * 6 + l];	      
		}	  
	      
	      term = log(0.25 * term) + (ex2[i] * log(minlikelihood));	   
	      
	      sum += wptr[i] * term;
	    }     	 
	}
    }              
  else
    {
      for (i = 0; i < n; i++) 
	{	  	 	             
      
	  for(j = 0, term = 0.0; j < 4; j++)
	    {
	      left  = &(x1[24 * i + 6 * j]);
	      right = &(x2[24 * i + 6 * j]);	    
	      
	      for(l = 0; l < 6; l++)
		term += left[l] * right[l] * diagptable[j * 6 + l];	
	    }
	  
	  term = log(0.25 * term) + ((ex1[i] + ex2[i])*log(minlikelihood));
	  
	  sum += wptr[i] * term;
	}         
    }
       
  return  sum;
}

static double evaluateGTRGAMMASECONDARY_7 (int *ex1, int *ex2, int *wptr,
					   double *x1, double *x2,  
					   double *tipVector, 
					   unsigned char *tipX1, int n, double *diagptable, double *vector, boolean writeVector)
{
  double   sum = 0.0, term;        
  int     i, j, l;   
  double  *left, *right;              
  
  if(tipX1)
    {          
      if(writeVector)
	for (i = 0; i < n; i++) 
	  {
	    left = &(tipVector[7 * tipX1[i]]);	  	  
	    
	    for(j = 0, term = 0.0; j < 4; j++)
	      {
		right = &(x2[28 * i + 7 * j]);
		
		for(l = 0; l < 7; l++)
		  term += left[l] * right[l] * diagptable[j * 7 + l];	      
	      }	  
	    
	    term = log(0.25 * term) + (ex2[i] * log(minlikelihood));	   
	    	    
	    vector[i] = term;
	    
	    sum += wptr[i] * term;
	  }         
      else
	{       
	  for (i = 0; i < n; i++) 
	    {	     
	      left = &(tipVector[7 * tipX1[i]]);	  	  
	      
	      for(j = 0, term = 0.0; j < 4; j++)
		{
		  right = &(x2[28 * i + 7 * j]);
		  
		  for(l = 0; l < 7; l++)
		    term += left[l] * right[l] * diagptable[j * 7 + l];	      
		}	  
	      
	      term = log(0.25 * term) + (ex2[i] * log(minlikelihood));	   
	      
	      sum += wptr[i] * term;
	    }     	 
	}
    }              
  else
    {
      for (i = 0; i < n; i++) 
	{	  	 	             
      
	  for(j = 0, term = 0.0; j < 4; j++)
	    {
	      left  = &(x1[28 * i + 7 * j]);
	      right = &(x2[28 * i + 7 * j]);	    
	      
	      for(l = 0; l < 7; l++)
		term += left[l] * right[l] * diagptable[j * 7 + l];	
	    }
	  
	  term = log(0.25 * term) + ((ex1[i] + ex2[i])*log(minlikelihood));
	  
	  sum += wptr[i] * term;
	}         
    }
       
  return  sum;
}

static double evaluateGTRGAMMAPROTINVAR (int *ex1, int *ex2, int *wptr, int *iptr,
					 double *x1, double *x2, 
					 double *tipVector,double *tFreqs, double invariants,
					 unsigned char *tipX1, int n, double *diagptable, double *vector, boolean writeVector)
{
  double   
    sum = 0.0, term, freqs[20],
    scaler = 0.25 * (1.0 - invariants);        
  int     i, j, l;     
  double *left, *right;   
    
  for(i = 0; i < 20; i++)
    freqs[i] = tFreqs[i] * invariants;            	  
  
  if(tipX1)
    {    
      if(writeVector)
	for (i = 0; i < n; i++) 
	  {
	    left = &(tipVector[20 * tipX1[i]]);
	    
	    for(j = 0, term = 0.0; j < 4; j++)
	      {
		right = &(x2[80 * i + 20 * j]);
		
		for(l = 0; l < 20; l++)
		  term += left[l] * right[l] * diagptable[j * 20 + l];	      
	      }	  
	    
	    if(iptr[i] < 20)	   
	      term = log(((scaler * term) + freqs[iptr[i]]) * pow(minlikelihood, ex2[i]));
	    else
	      term = log(scaler * term) + (ex2[i] * log(minlikelihood));
	    	    
	    vector[i] = term;
	   
	    sum += wptr[i] * term;
	  }         
      else
	for (i = 0; i < n; i++) 
	  {
	    left = &(tipVector[20 * tipX1[i]]);
	    
	    for(j = 0, term = 0.0; j < 4; j++)
	      {
		right = &(x2[80 * i + 20 * j]);
		
		for(l = 0; l < 20; l++)
		  term += left[l] * right[l] * diagptable[j * 20 + l];	      
	      }	  
	    
	    if(iptr[i] < 20)	   
	      term = log(((scaler * term) + freqs[iptr[i]]) * pow(minlikelihood, ex2[i]));
	    else
	      term = log(scaler * term) + (ex2[i] * log(minlikelihood));
	    	    
	    sum += wptr[i] * term;
	  }    	
    }                
  else
    {    
      for (i = 0; i < n; i++) 
	{	  	 	       	  
	  for(j = 0, term = 0.0; j < 4; j++)
	    {
	      left  = &(x1[80 * i + 20 * j]);
	      right = &(x2[80 * i + 20 * j]);	    
	      
	      for(l = 0; l < 20; l++)
		term += left[l] * right[l] * diagptable[j * 20 + l];	
	    }
	  
	  if(iptr[i] < 20)	   
	    term = log(((scaler * term) + freqs[iptr[i]]) * pow(minlikelihood, ex2[i] + ex1[i]));
	  else
	    term = log(scaler * term) + ((ex1[i] + ex2[i]) * log(minlikelihood));
	  sum += wptr[i] * term;
	}              
    }
       
  return  sum;
}

static double evaluateGTRGAMMASECONDARYINVAR (int *ex1, int *ex2, int *wptr, int *iptr,
					      double *x1, double *x2, 
					      double *tipVector,double *tFreqs, double invariants,
					      unsigned char *tipX1, int n, double *diagptable, double *vector, boolean writeVector)
{
  double   
    sum = 0.0, term, freqs[16],
    scaler = 0.25 * (1.0 - invariants);        
  int     i, j, l;     
  double *left, *right;   
    
  for(i = 0; i < 16; i++)
    freqs[i] = tFreqs[i] * invariants;            	  
  
  if(tipX1)
    {    
      if(writeVector)
	for (i = 0; i < n; i++) 
	  {
	    left = &(tipVector[16 * tipX1[i]]);
	    
	    for(j = 0, term = 0.0; j < 4; j++)
	      {
		right = &(x2[64 * i + 16 * j]);
		
		for(l = 0; l < 16; l++)
		  term += left[l] * right[l] * diagptable[j * 16 + l];	      
	      }	  
	    
	    if(iptr[i] < 16)	   
	      term = log(((scaler * term) + freqs[iptr[i]]) * pow(minlikelihood, ex2[i]));
	    else
	      term = log(scaler * term) + (ex2[i] * log(minlikelihood));
	    	    
	    vector[i] = term;
	   
	    sum += wptr[i] * term;
	  }         
      else
	for (i = 0; i < n; i++) 
	  {
	    left = &(tipVector[16 * tipX1[i]]);
	    
	    for(j = 0, term = 0.0; j < 4; j++)
	      {
		right = &(x2[64 * i + 16 * j]);
		
		for(l = 0; l < 16; l++)
		  term += left[l] * right[l] * diagptable[j * 16 + l];	      
	      }	  
	    
	    if(iptr[i] < 16)	   
	      term = log(((scaler * term) + freqs[iptr[i]]) * pow(minlikelihood, ex2[i]));
	    else
	      term = log(scaler * term) + (ex2[i] * log(minlikelihood));
	    	    
	    sum += wptr[i] * term;
	  }    	
    }                
  else
    {    
      for (i = 0; i < n; i++) 
	{	  	 	       	  
	  for(j = 0, term = 0.0; j < 4; j++)
	    {
	      left  = &(x1[64 * i + 16 * j]);
	      right = &(x2[64 * i + 16 * j]);	    
	      
	      for(l = 0; l < 16; l++)
		term += left[l] * right[l] * diagptable[j * 16 + l];	
	    }
	  
	  if(iptr[i] < 16)	   
	    term = log(((scaler * term) + freqs[iptr[i]]) * pow(minlikelihood, ex2[i] + ex1[i]));
	  else
	    term = log(scaler * term) + ((ex1[i] + ex2[i]) * log(minlikelihood));
	  sum += wptr[i] * term;
	}              
    }
       
  return  sum;
}

static double evaluateGTRGAMMASECONDARYINVAR_6 (int *ex1, int *ex2, int *wptr, int *iptr,
						double *x1, double *x2, 
						double *tipVector,double *tFreqs, double invariants,
						unsigned char *tipX1, int n, double *diagptable, double *vector, boolean writeVector)
{
  double   
    sum = 0.0, term, freqs[6],
    scaler = 0.25 * (1.0 - invariants);        
  int     i, j, l;     
  double *left, *right;   
    
  for(i = 0; i < 6; i++)
    freqs[i] = tFreqs[i] * invariants;            	  
  
  if(tipX1)
    {    
      if(writeVector)
	for (i = 0; i < n; i++) 
	  {
	    left = &(tipVector[6 * tipX1[i]]);
	    
	    for(j = 0, term = 0.0; j < 4; j++)
	      {
		right = &(x2[24 * i + 6 * j]);
		
		for(l = 0; l < 6; l++)
		  term += left[l] * right[l] * diagptable[j * 6 + l];	      
	      }	  
	    
	    if(iptr[i] < 6)	   
	      term = log(((scaler * term) + freqs[iptr[i]]) * pow(minlikelihood, ex2[i]));
	    else
	      term = log(scaler * term) + (ex2[i] * log(minlikelihood));
	    	    
	    vector[i] = term;
	   
	    sum += wptr[i] * term;
	  }         
      else
	for (i = 0; i < n; i++) 
	  {
	    left = &(tipVector[6 * tipX1[i]]);
	    
	    for(j = 0, term = 0.0; j < 4; j++)
	      {
		right = &(x2[24 * i + 6 * j]);
		
		for(l = 0; l < 6; l++)
		  term += left[l] * right[l] * diagptable[j * 6 + l];	      
	      }	  
	    
	    if(iptr[i] < 6)	   
	      term = log(((scaler * term) + freqs[iptr[i]]) * pow(minlikelihood, ex2[i]));
	    else
	      term = log(scaler * term) + (ex2[i] * log(minlikelihood));
	    	    
	    sum += wptr[i] * term;
	  }    	
    }                
  else
    {    
      for (i = 0; i < n; i++) 
	{	  	 	       	  
	  for(j = 0, term = 0.0; j < 4; j++)
	    {
	      left  = &(x1[24 * i + 6 * j]);
	      right = &(x2[24 * i + 6 * j]);	    
	      
	      for(l = 0; l < 6; l++)
		term += left[l] * right[l] * diagptable[j * 6 + l];	
	    }
	  
	  if(iptr[i] < 6)	   
	    term = log(((scaler * term) + freqs[iptr[i]]) * pow(minlikelihood, ex2[i] + ex1[i]));
	  else
	    term = log(scaler * term) + ((ex1[i] + ex2[i]) * log(minlikelihood));
	  sum += wptr[i] * term;
	}              
    }
       
  return  sum;
}

static double evaluateGTRGAMMASECONDARYINVAR_7 (int *ex1, int *ex2, int *wptr, int *iptr,
						double *x1, double *x2, 
						double *tipVector,double *tFreqs, double invariants,
						unsigned char *tipX1, int n, double *diagptable, double *vector, boolean writeVector)
{
  double   
    sum = 0.0, term, freqs[7],
    scaler = 0.25 * (1.0 - invariants);        
  int     i, j, l;     
  double *left, *right;   
    
  for(i = 0; i < 7; i++)
    freqs[i] = tFreqs[i] * invariants;            	  
  
  if(tipX1)
    {    
      if(writeVector)
	for (i = 0; i < n; i++) 
	  {
	    left = &(tipVector[7 * tipX1[i]]);
	    
	    for(j = 0, term = 0.0; j < 4; j++)
	      {
		right = &(x2[28 * i + 7 * j]);
		
		for(l = 0; l < 7; l++)
		  term += left[l] * right[l] * diagptable[j * 7 + l];	      
	      }	  
	    
	    if(iptr[i] < 7)	   
	      term = log(((scaler * term) + freqs[iptr[i]]) * pow(minlikelihood, ex2[i]));
	    else
	      term = log(scaler * term) + (ex2[i] * log(minlikelihood));
	    	    
	    vector[i] = term;
	   
	    sum += wptr[i] * term;
	  }         
      else
	for (i = 0; i < n; i++) 
	  {
	    left = &(tipVector[7 * tipX1[i]]);
	    
	    for(j = 0, term = 0.0; j < 4; j++)
	      {
		right = &(x2[28 * i + 7 * j]);
		
		for(l = 0; l < 7; l++)
		  term += left[l] * right[l] * diagptable[j * 7 + l];	      
	      }	  
	    
	    if(iptr[i] < 7)	   
	      term = log(((scaler * term) + freqs[iptr[i]]) * pow(minlikelihood, ex2[i]));
	    else
	      term = log(scaler * term) + (ex2[i] * log(minlikelihood));
	    	    
	    sum += wptr[i] * term;
	  }    	
    }                
  else
    {    
      for (i = 0; i < n; i++) 
	{	  	 	       	  
	  for(j = 0, term = 0.0; j < 4; j++)
	    {
	      left  = &(x1[28 * i + 7 * j]);
	      right = &(x2[28 * i + 7 * j]);	    
	      
	      for(l = 0; l < 7; l++)
		term += left[l] * right[l] * diagptable[j * 7 + l];	
	    }
	  
	  if(iptr[i] < 7)	   
	    term = log(((scaler * term) + freqs[iptr[i]]) * pow(minlikelihood, ex2[i] + ex1[i]));
	  else
	    term = log(scaler * term) + ((ex1[i] + ex2[i]) * log(minlikelihood));
	  sum += wptr[i] * term;
	}              
    }
       
  return  sum;
}


double evaluateIterative(tree *tr,  boolean writeVector)
{
  double 
    result = 0.0;  
  int pNumber, qNumber, model;
  double *pz;

  pNumber = tr->td[0].ti[0].pNumber;
  qNumber = tr->td[0].ti[0].qNumber;
  pz      = tr->td[0].ti[0].qz;

  newviewIterative(tr); 

  for(model = 0; model < tr->NumberOfModels; model++)
    {            
      if(tr->executeModel[model])
	{	
	  int width = tr->partitionData[model].width;
	  double z, partitionLikelihood, *_vector;
	  int    
	    *ex1 = (int*)NULL, 
	    *ex2 = (int*)NULL;
	  double 
	    *x1_start = (double*)NULL, 
	    *x2_start = (double*)NULL;
	  unsigned char *tip = (unsigned char*)NULL;
	  if(writeVector)
	    _vector = tr->partitionData[model].perSiteLL;
	  else
	    _vector = (double*)NULL;

	  if(isTip(pNumber, tr->mxtips) || isTip(qNumber, tr->mxtips))
	    {	        	    
	      if(isTip(qNumber, tr->mxtips))
		{	
		  
		  x2_start = tr->partitionData[model].xVector[pNumber - tr->mxtips -1];
		  ex2      = tr->partitionData[model].expVector[pNumber - tr->mxtips - 1];
		  
		  tip = tr->partitionData[model].yVector[qNumber];	 	      
		}           
	      else
		{
		  x2_start = tr->partitionData[model].xVector[qNumber - tr->mxtips - 1];
		  ex2      = tr->partitionData[model].expVector[qNumber - tr->mxtips - 1];	 
		  
		  tip = tr->partitionData[model].yVector[pNumber];
		}
	    }
	  else
	    {                 
	      x1_start = tr->partitionData[model].xVector[pNumber - tr->mxtips - 1];
	      x2_start = tr->partitionData[model].xVector[qNumber - tr->mxtips - 1];
	      ex1      = tr->partitionData[model].expVector[pNumber - tr->mxtips - 1];
	      ex2      = tr->partitionData[model].expVector[qNumber - tr->mxtips - 1];     
	    }


	  if(tr->multiBranch)
	    z = pz[model];
	  else
	    z = pz[0];

	  switch(tr->partitionData[model].dataType)
	    { 
	    case BINARY_DATA:
	       switch(tr->rateHetModel)
		{
		case CAT:	    
		  {
		    double *diagptable = (double*)malloc(sizeof(double) * 2 * tr->NumberOfCategories);
		    
		    calcDiagptable(z, BINARY_DATA, tr->NumberOfCategories, tr->cdta->patrat, tr->partitionData[model].EIGN, diagptable);
		    
		    partitionLikelihood =  evaluateGTRCAT_BINARY(ex1, ex2, tr->partitionData[model].rateCategory, tr->partitionData[model].wgt,
								 x1_start, x2_start, tr->partitionData[model].tipVector, 
								 tip, width, diagptable, _vector, writeVector);

		    free(diagptable);
		  }
		  break;	  	   
		case GAMMA:	   
		  {
		    double *diagptable = (double*)malloc(sizeof(double) * 8);
		    
		    calcDiagptable(z, BINARY_DATA, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);		    		    

		    partitionLikelihood = evaluateGTRGAMMA_BINARY(ex1, ex2, tr->partitionData[model].wgt,
								  x1_start, x2_start, tr->partitionData[model].tipVector,
								  tip, width, diagptable, _vector, writeVector); 
		    free(diagptable);
		  }
		  break; 
		case GAMMA_I:
		  {
		    double *diagptable = (double*)malloc(sizeof(double) * 8);
		    
		    calcDiagptable(z, BINARY_DATA, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);

		    partitionLikelihood = evaluateGTRGAMMAINVAR_BINARY(ex1, ex2, tr->partitionData[model].wgt, tr->partitionData[model].invariant,
								       x1_start, x2_start,
								       tr->partitionData[model].tipVector, tr->partitionData[model].frequencies, 
								       tr->partitionData[model].propInvariant,
								       tip, width, diagptable, _vector, writeVector);

		    free(diagptable);
		  }
		  break;
		default:
		  assert(0);
		}
	      break;
	    case DNA_DATA:
	      switch(tr->rateHetModel)
		{
		case CAT:	    
		  {
		    double *diagptable = (double*)malloc(sizeof(double) * 4 * tr->NumberOfCategories);
		    
		    calcDiagptable(z, DNA_DATA, tr->NumberOfCategories, tr->cdta->patrat, tr->partitionData[model].EIGN, diagptable);
		    
		    partitionLikelihood =  evaluateGTRCAT(ex1, ex2, tr->partitionData[model].rateCategory, tr->partitionData[model].wgt,
							  x1_start, x2_start, tr->partitionData[model].tipVector, 
							  tip, width, diagptable, _vector, writeVector);

		    free(diagptable);
		  }
		  break;	  	   
		case GAMMA:	   
		  {
		    double *diagptable = (double*)malloc(sizeof(double) * 16);

		    calcDiagptable(z, DNA_DATA, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);		    		    

		    partitionLikelihood = evaluateGTRGAMMA(ex1, ex2, tr->partitionData[model].wgt,
							   x1_start, x2_start, tr->partitionData[model].tipVector,
							   tip, width, diagptable, _vector, writeVector); 
		    free(diagptable);
		  }
		  break; 
		case GAMMA_I:
		  {
		    double *diagptable = (double*)malloc(sizeof(double) * 16);
		    
		    calcDiagptable(z, DNA_DATA, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);

		    partitionLikelihood = evaluateGTRGAMMAINVAR(ex1, ex2, tr->partitionData[model].wgt, tr->partitionData[model].invariant,
								x1_start, x2_start,
								tr->partitionData[model].tipVector, tr->partitionData[model].frequencies, 
								tr->partitionData[model].propInvariant,
								tip, width, diagptable, _vector, writeVector);

		    free(diagptable);
		  }
		  break;
		default:
		  assert(0);
		}
	      break;
	    case AA_DATA:
	      switch(tr->rateHetModel)
		{
		case CAT:	    
		  {
		    double *diagptable = (double*)malloc(sizeof(double) * 20 * tr->NumberOfCategories);

		    calcDiagptable(z, AA_DATA, tr->NumberOfCategories, tr->cdta->patrat, tr->partitionData[model].EIGN, diagptable);

		    partitionLikelihood = evaluateGTRCATPROT(ex1, ex2, tr->partitionData[model].rateCategory, tr->partitionData[model].wgt,
							     x1_start, x2_start, tr->partitionData[model].tipVector,
							     tip, width, diagptable, _vector, writeVector);
		    free(diagptable);
		  }	     	      
		  break;	      
		case GAMMA:
		  {
		    double *diagptable = (double*)malloc(sizeof(double) * 80);

		    calcDiagptable(z, AA_DATA, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);

		    partitionLikelihood = evaluateGTRGAMMAPROT(ex1, ex2, tr->partitionData[model].wgt,
							       x1_start, x2_start, tr->partitionData[model].tipVector,
							       tip, width, diagptable, _vector, writeVector);
		    free(diagptable);		    
		  }
		  break;
		case GAMMA_I:		  	    
		  {
		    double *diagptable = (double*)malloc(sizeof(double) * 80);

		    calcDiagptable(z, AA_DATA, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);
		    
		    partitionLikelihood = evaluateGTRGAMMAPROTINVAR(ex1, ex2, tr->partitionData[model].wgt, tr->partitionData[model].invariant,
								    x1_start, x2_start, 
								    tr->partitionData[model].tipVector, tr->partitionData[model].frequencies, 
								    tr->partitionData[model].propInvariant, 
								    tip, width, diagptable, _vector, writeVector); 

		    free(diagptable);
		  }	  
		  break;
		default:
		  assert(0);
		}
	      break;
	    case SECONDARY_DATA:
	      switch(tr->rateHetModel)
		{
		case CAT:	    
		  {
		    double *diagptable = (double*)malloc(sizeof(double) * 16 * tr->NumberOfCategories);

		    calcDiagptable(z, SECONDARY_DATA, tr->NumberOfCategories, tr->cdta->patrat, tr->partitionData[model].EIGN, diagptable);

		    partitionLikelihood = evaluateGTRCATSECONDARY(ex1, ex2, tr->partitionData[model].rateCategory, tr->partitionData[model].wgt,
								  x1_start, x2_start, tr->partitionData[model].tipVector,
								  tip, width, diagptable, _vector, writeVector);
		    free(diagptable);
		  }	     	      
		  break;	      
		case GAMMA:
		  {
		    double *diagptable = (double*)malloc(sizeof(double) * 64);

		    calcDiagptable(z, SECONDARY_DATA, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);

		    partitionLikelihood = evaluateGTRGAMMASECONDARY(ex1, ex2, tr->partitionData[model].wgt,
								    x1_start, x2_start, tr->partitionData[model].tipVector,
								    tip, width, diagptable, _vector, writeVector);
		    free(diagptable);		    
		  }
		  break;
		case GAMMA_I:		  	    
		  {
		    double *diagptable = (double*)malloc(sizeof(double) * 64);

		    calcDiagptable(z, SECONDARY_DATA, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);
		    
		    partitionLikelihood = evaluateGTRGAMMASECONDARYINVAR(ex1, ex2, tr->partitionData[model].wgt, tr->partitionData[model].invariant,
									 x1_start, x2_start, 
									 tr->partitionData[model].tipVector, tr->partitionData[model].frequencies, 
									 tr->partitionData[model].propInvariant, 
									 tip, width, diagptable, _vector, writeVector);

		    free(diagptable);
		  }	  
		  break;
		default:
		  assert(0);
		}
	      break;
	    case SECONDARY_DATA_6:
	      switch(tr->rateHetModel)
		{
		case CAT:	    
		  {
		    double *diagptable = (double*)malloc(sizeof(double) * 6 * tr->NumberOfCategories);

		    calcDiagptable(z, SECONDARY_DATA_6, tr->NumberOfCategories, tr->cdta->patrat, tr->partitionData[model].EIGN, diagptable);

		    partitionLikelihood = evaluateGTRCATSECONDARY_6(ex1, ex2, tr->partitionData[model].rateCategory, tr->partitionData[model].wgt,
								    x1_start, x2_start, tr->partitionData[model].tipVector,
								    tip, width, diagptable, _vector, writeVector);
		    free(diagptable);
		  }	     	      
		  break;	      
		case GAMMA:
		  {
		    double *diagptable = (double*)malloc(sizeof(double) * 24);

		    calcDiagptable(z, SECONDARY_DATA_6, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);

		    partitionLikelihood = evaluateGTRGAMMASECONDARY_6(ex1, ex2, tr->partitionData[model].wgt,
								    x1_start, x2_start, tr->partitionData[model].tipVector,
								    tip, width, diagptable, _vector, writeVector);
		    free(diagptable);		    
		  }
		  break;
		case GAMMA_I:		  	    
		  {
		    double *diagptable = (double*)malloc(sizeof(double) * 24);

		    calcDiagptable(z, SECONDARY_DATA_6, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);
		    
		    partitionLikelihood = evaluateGTRGAMMASECONDARYINVAR_6(ex1, ex2, tr->partitionData[model].wgt, tr->partitionData[model].invariant,
									   x1_start, x2_start, 
									   tr->partitionData[model].tipVector, tr->partitionData[model].frequencies, 
									   tr->partitionData[model].propInvariant, 
									   tip, width, diagptable, _vector, writeVector);

		    free(diagptable);
		  }	  
		  break;
		default:
		  assert(0);
		}
	      break;
	    case SECONDARY_DATA_7:
	      switch(tr->rateHetModel)
		{
		case CAT:	    
		  {
		    double *diagptable = (double*)malloc(sizeof(double) * 7 * tr->NumberOfCategories);

		    calcDiagptable(z, SECONDARY_DATA_7, tr->NumberOfCategories, tr->cdta->patrat, tr->partitionData[model].EIGN, diagptable);

		    partitionLikelihood = evaluateGTRCATSECONDARY_7(ex1, ex2, tr->partitionData[model].rateCategory, tr->partitionData[model].wgt,
								    x1_start, x2_start, tr->partitionData[model].tipVector,
								    tip, width, diagptable, _vector, writeVector);
		    free(diagptable);
		  }	     	      
		  break;	      
		case GAMMA:
		  {
		    double *diagptable = (double*)malloc(sizeof(double) * 28);

		    calcDiagptable(z, SECONDARY_DATA_7, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);

		    partitionLikelihood = evaluateGTRGAMMASECONDARY_7(ex1, ex2, tr->partitionData[model].wgt,
								      x1_start, x2_start, tr->partitionData[model].tipVector,
								      tip, width, diagptable, _vector, writeVector);
		    free(diagptable);		    
		  }
		  break;
		case GAMMA_I:		  	    
		  {
		    double *diagptable = (double*)malloc(sizeof(double) * 28);

		    calcDiagptable(z, SECONDARY_DATA_7, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);
		    
		    partitionLikelihood = evaluateGTRGAMMASECONDARYINVAR_7(ex1, ex2, tr->partitionData[model].wgt, tr->partitionData[model].invariant,
									   x1_start, x2_start, 
									   tr->partitionData[model].tipVector, tr->partitionData[model].frequencies, 
									   tr->partitionData[model].propInvariant, 
									   tip, width, diagptable, _vector, writeVector);

		    free(diagptable);
		  }	  
		  break;
		default:
		  assert(0);
		}
	      break;
	    default:
	      assert(0);
	    }
	  result += partitionLikelihood;	  
	  tr->perPartitionLH[model] = partitionLikelihood;
	  /*printf("TRUE %d %f\n", model, result);*/
	}
      /*else*/
	/*printf("FALSE %d %f %f\n", model, result,  tr->perPartitionLH[model]);*/
      
    }
      
  return result;
}

#ifdef _USE_PTHREADS

double evaluateClassify(tree *tr,  branchInfo *b)
{
  double 
    result = 0.0;  
  int pNumber, qNumber, model, columnCounter, offsetCounter;
  double *pz; 
  double *_vector = (double*)NULL;
  boolean writeVector = FALSE;

  /* get the node numbers for this branch and the branch lengths */

  pNumber = b->leftNodeNumber;
  qNumber = b->rightNodeNumber;
  pz      = b->branchLengths;


  /* once again loop over all partitions */

  for(model = 0, columnCounter = 0, offsetCounter = 0; model < tr->NumberOfModels; model++)
    {                
      /* get the real contiguous width of this partition */
      int width = tr->partitionData[model].upper - tr->partitionData[model].lower;
      double z, partitionLikelihood;      
      int    
	increment, 
	*ex1 = (int*)NULL, 
	*ex2 = (int*)NULL;
      double 
	*x1_start = (double*)NULL, 
	*x2_start = (double*)NULL;
      unsigned char *tip = (unsigned char*)NULL;

      /* get the additional contiguous pointers we need to compute the likelihood */

      int *rateCategory = &tr->contiguousRateCategory[columnCounter];
      int *wgt          = &tr->contiguousWgt[columnCounter];
      int *invariant    = &tr->contiguousInvariant[columnCounter];           

      /* 
	 do the switch again to figure out if one of the nodes at the end 
	 of the current branch is a tip. we also need to make sure that 
	 we get the offsets right for this partition .....
      */

      if(isTip(pNumber, tr->mxtips) || isTip(qNumber, tr->mxtips))
	{	  	        	    
	  if(isTip(qNumber, tr->mxtips))
	    {	
	      
	      x2_start = &b->left[offsetCounter]; 
	      ex2      = &b->leftScaling[columnCounter]; 	      

	      /* actually need to set the pointer to the DNA or AA data at the tip here 
		 this is just an matrix of unsigned chars that corresponds to the input 
		 alignment 
	      */

	      tip = &tr->contiguousTips[qNumber][columnCounter];	      
	    }           
	  else
	    {
	      x2_start = &b->right[offsetCounter];        
	      ex2      = &b->rightScaling[columnCounter]; 	      

	      tip = &tr->contiguousTips[pNumber][columnCounter];
	    }
	}
      else
	{                 
	  x1_start = &b->left[offsetCounter];
	  x2_start = &b->right[offsetCounter];
	  ex1      = &b->leftScaling[columnCounter];
	  ex2      = &b->rightScaling[columnCounter];
	}
     
      /* no need to worry about this now */

      if(tr->multiBranch)
	z = pz[model];
      else
	z = pz[0];

      switch(tr->partitionData[model].dataType)
	{ 
	  /* big switch over the data type and rate heterogeneity model used */
	case BINARY_DATA:
	  switch(tr->rateHetModel)
	    {
	    case CAT:	    
	      {				
		double *diagptable = (double*)malloc(sizeof(double) * 2 * tr->NumberOfCategories);

		/* compute the transition probability matrix P(z) given z */

		calcDiagptable(z, BINARY_DATA, tr->NumberOfCategories, tr->cdta->patrat, tr->partitionData[model].EIGN, diagptable);
		
		/* now compute the likelihood score for this partition, I haven't changed anything in the function that is called  */
		
		partitionLikelihood =  evaluateGTRCAT_BINARY(ex1, ex2, rateCategory, wgt,
							     x1_start, x2_start, tr->partitionData[model].tipVector, 
							     tip, width, diagptable, _vector, writeVector);

		/* and set the likelihood vector increment appropriately */

		increment = 2;
		free(diagptable);
	      }
	      break;	  	   
	    case GAMMA:	   
	      {
		/* same stuff different model */

		double *diagptable = (double*)malloc(sizeof(double) * 8);
		
		calcDiagptable(z, BINARY_DATA, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);		    		    
		
		partitionLikelihood = evaluateGTRGAMMA_BINARY(ex1, ex2, wgt,
							      x1_start, x2_start, tr->partitionData[model].tipVector,
							      tip, width, diagptable, _vector, writeVector); 
		increment = 8;
		free(diagptable);
	      }
	      break; 
	    case GAMMA_I:
	      {
		double *diagptable = (double*)malloc(sizeof(double) * 8);
		
		calcDiagptable(z, BINARY_DATA, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);
		
		partitionLikelihood = evaluateGTRGAMMAINVAR_BINARY(ex1, ex2, wgt, invariant,
								   x1_start, x2_start,
								   tr->partitionData[model].tipVector, tr->partitionData[model].frequencies, 
								   tr->partitionData[model].propInvariant,
								   tip, width, diagptable, _vector, writeVector);
		increment = 8;
		free(diagptable);
	      }
	      break;
	    default:
	      assert(0);
	    }
	  break;
	case DNA_DATA:
	  switch(tr->rateHetModel)
	    {
	    case CAT:	    
	      {
		double *diagptable = (double*)malloc(sizeof(double) * 4 * tr->NumberOfCategories);
		
		calcDiagptable(z, DNA_DATA, tr->NumberOfCategories, tr->cdta->patrat, tr->partitionData[model].EIGN, diagptable);
		
		partitionLikelihood =  evaluateGTRCAT(ex1, ex2, rateCategory, wgt,
						      x1_start, x2_start, tr->partitionData[model].tipVector, 
						      tip, width, diagptable, _vector, writeVector);
		
		increment = 4;
		free(diagptable);
	      }
	      break;	  	   
	    case GAMMA:	   
	      {
		double *diagptable = (double*)malloc(sizeof(double) * 16);
		
		calcDiagptable(z, DNA_DATA, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);		    		    
		
		partitionLikelihood = evaluateGTRGAMMA(ex1, ex2, wgt,
						       x1_start, x2_start, tr->partitionData[model].tipVector,
						       tip, width, diagptable, _vector, writeVector); 
		increment = 16;
		free(diagptable);
	      }
	      break; 
	    case GAMMA_I:
	      {
		double *diagptable = (double*)malloc(sizeof(double) * 16);
		
		calcDiagptable(z, DNA_DATA, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);
		
		partitionLikelihood = evaluateGTRGAMMAINVAR(ex1, ex2, wgt, invariant,
							    x1_start, x2_start,
							    tr->partitionData[model].tipVector, tr->partitionData[model].frequencies, 
							    tr->partitionData[model].propInvariant,
							    tip, width, diagptable, _vector, writeVector);
		increment = 16;
		free(diagptable);
	      }
	      break;
	    default:
	      assert(0);
	    }
	  break;
	case AA_DATA:
	  switch(tr->rateHetModel)
	    {
	    case CAT:	    
	      {
		double *diagptable = (double*)malloc(sizeof(double) * 20 * tr->NumberOfCategories);
		
		calcDiagptable(z, AA_DATA, tr->NumberOfCategories, tr->cdta->patrat, tr->partitionData[model].EIGN, diagptable);
		
		partitionLikelihood = evaluateGTRCATPROT(ex1, ex2, rateCategory, wgt,
							 x1_start, x2_start, tr->partitionData[model].tipVector,
							 tip, width, diagptable, _vector, writeVector);
		increment = 20;
		free(diagptable);
	      }	     	      
	      break;	      
	    case GAMMA:
	      {
		double *diagptable = (double*)malloc(sizeof(double) * 80);
		
		calcDiagptable(z, AA_DATA, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);
		
		partitionLikelihood = evaluateGTRGAMMAPROT(ex1, ex2, wgt,
							   x1_start, x2_start, tr->partitionData[model].tipVector,
							   tip, width, diagptable, _vector, writeVector);
		increment = 80;
		free(diagptable);		    
	      }
	      break;
	    case GAMMA_I:		  	    
	      {
		double *diagptable = (double*)malloc(sizeof(double) * 80);
		
		calcDiagptable(z, AA_DATA, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);
		
		partitionLikelihood = evaluateGTRGAMMAPROTINVAR(ex1, ex2, wgt, invariant,
								x1_start, x2_start, 
								tr->partitionData[model].tipVector, tr->partitionData[model].frequencies, 
								tr->partitionData[model].propInvariant, 
								tip, width, diagptable, _vector, writeVector); 
		increment = 80;
		free(diagptable);
	      }	  
	      break;
	    default:
	      assert(0);
	    }
	  break;
	case SECONDARY_DATA:
	  switch(tr->rateHetModel)
	    {
	    case CAT:	    
	      {
		double *diagptable = (double*)malloc(sizeof(double) * 16 * tr->NumberOfCategories);
		
		calcDiagptable(z, SECONDARY_DATA, tr->NumberOfCategories, tr->cdta->patrat, tr->partitionData[model].EIGN, diagptable);
		
		partitionLikelihood = evaluateGTRCATSECONDARY(ex1, ex2, rateCategory, wgt,
							      x1_start, x2_start, tr->partitionData[model].tipVector,
							      tip, width, diagptable, _vector, writeVector);
		increment = 16;
		free(diagptable);
	      }	     	      
	      break;	      
	    case GAMMA:
	      {
		double *diagptable = (double*)malloc(sizeof(double) * 64);
		
		calcDiagptable(z, SECONDARY_DATA, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);
		
		partitionLikelihood = evaluateGTRGAMMASECONDARY(ex1, ex2, wgt,
								x1_start, x2_start, tr->partitionData[model].tipVector,
								tip, width, diagptable, _vector, writeVector);
		increment = 64;
		free(diagptable);		    
	      }
	      break;
	    case GAMMA_I:		  	    
	      {
		double *diagptable = (double*)malloc(sizeof(double) * 64);
		
		calcDiagptable(z, SECONDARY_DATA, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);
		
		partitionLikelihood = evaluateGTRGAMMASECONDARYINVAR(ex1, ex2, wgt, invariant,
								     x1_start, x2_start, 
								     tr->partitionData[model].tipVector, tr->partitionData[model].frequencies, 
								     tr->partitionData[model].propInvariant, 
								     tip, width, diagptable, _vector, writeVector);
		increment = 64;
		free(diagptable);
	      }	  
	      break;
	    default:
	      assert(0);
	    }
	  break;
	default:
	  assert(0);
	}

      /* just accumulate over partition likelihoods and also store per-partition likelihood separately */
      result += partitionLikelihood;	  
      tr->perPartitionLH[model] = partitionLikelihood;
      /*printf("TRUE %d %f\n", model, result);*/

      /* increment the offset counters after we have finished with this partition */

      columnCounter += width;
      offsetCounter += width * increment;
    }
  
  /* return the likelihood :-) */

  return result;
}

#endif




double evaluateGeneric (tree *tr, nodeptr p)
{
  volatile double result;
  nodeptr q = p->back; 
  int i;
    
  tr->td[0].ti[0].pNumber = p->number;
  tr->td[0].ti[0].qNumber = q->number;          
  
  for(i = 0; i < tr->numBranches; i++)    
    tr->td[0].ti[0].qz[i] =  q->z[i];
  
  tr->td[0].count = 1;
  if(!p->x)
    computeTraversalInfo(p, &(tr->td[0].ti[0]), &(tr->td[0].count), tr->mxtips, tr->numBranches);
  if(!q->x)
    computeTraversalInfo(q, &(tr->td[0].ti[0]), &(tr->td[0].count), tr->mxtips, tr->numBranches);  
  
#ifdef _USE_PTHREADS 
 {
    int j;
    
    masterBarrier(THREAD_EVALUATE, tr); 
    if(tr->NumberOfModels == 1)
      {
	for(i = 0, result = 0.0; i < NumberOfThreads; i++)          
	  result += reductionBuffer[i];  	  	     
      
	tr->perPartitionLH[0] = result;
      }
    else
      {
	volatile double partitionResult;
	
	result = 0.0;

	for(j = 0; j < tr->NumberOfModels; j++)
	  {
	    for(i = 0, partitionResult = 0.0; i < NumberOfThreads; i++)          	      
	      partitionResult += reductionBuffer[i * tr->NumberOfModels + j];
	    result += partitionResult;
	    tr->perPartitionLH[j] = partitionResult;
	  }
      }
  }  
#else
  result = evaluateIterative(tr, FALSE);
#endif
  
  tr->likelihood = result;    

  return result;
}

double evaluateGenericInitrav (tree *tr, nodeptr p)
{
  volatile double result;   
    
  determineFullTraversal(p, tr);
      
#ifdef _USE_PTHREADS 
  {
    int i, j;
    
    masterBarrier(THREAD_EVALUATE, tr); 

    if(tr->NumberOfModels == 1)
      {
	for(i = 0, result = 0.0; i < NumberOfThreads; i++)          
	  result += reductionBuffer[i];  	  	     
      
	tr->perPartitionLH[0] = result;
      }
    else
      {
	volatile double partitionResult;
	
	result = 0.0;

	for(j = 0; j < tr->NumberOfModels; j++)
	  {
	    for(i = 0, partitionResult = 0.0; i < NumberOfThreads; i++)          	      
	      partitionResult += reductionBuffer[i * tr->NumberOfModels + j];
	    result +=  partitionResult;
	    tr->perPartitionLH[j] = partitionResult;
	  }
      }
  }
#else
  result = evaluateIterative(tr, FALSE);
#endif
     
  tr->likelihood = result;    
    
  return result;
}


void onlyInitrav(tree *tr, nodeptr p)
{   
  determineFullTraversal(p, tr);  

#ifdef _USE_PTHREADS  
  masterBarrier(THREAD_NEWVIEW, tr);  	 
#else
  newviewIterative(tr);   
#endif     
}



void evaluateGenericVector (tree *tr, nodeptr p)
{
  nodeptr q = p->back;
  int i;  
    
  assert(isTip(p->number, tr->mxtips)  || isTip(q->number, tr->mxtips));
      
  if(isTip(q->number, tr->mxtips))
    {
      nodeptr tmp = p;
      p = q;
      q = tmp;	      
    }   
  
  tr->td[0].ti[0].pNumber = p->number;
  tr->td[0].ti[0].qNumber = q->number;

  for(i = 0; i < tr->numBranches; i++)    
    tr->td[0].ti[0].qz[i] =  q->z[i];

  tr->td[0].count = 1;
  if(!q->x)
    computeTraversalInfo(q, &(tr->td[0].ti[0]), &(tr->td[0].count), tr->mxtips, tr->numBranches); 
 
#ifdef _USE_PTHREADS
  masterBarrier(THREAD_EVALUATE_VECTOR, tr);
#else
  evaluateIterative(tr, TRUE); 
#endif
}




