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
#include <pthread.h>
extern volatile int NumberOfThreads;
#endif


/* C-OPT this function below will also be called very frequently, but should not take much time 
   to compute, nonetheless it will be interesting cache-wise */

static void makeP(double z1, double z2, double *rptr, double *EI,  double *EIGN, int numberOfCategories, double *left, double *right, int data)
{
  int i, j, k; 

  switch(data)
    {
    case BINARY_DATA:
      {
	double d1, d2;
	   
	for(i = 0; i < numberOfCategories; i++)
	  {		                      	   
	    d1 = exp(rptr[i] * EIGN[0] * z1);
	    d2 = exp(rptr[i] * EIGN[0] * z2);	                       
	    
	    /*printf("EIGN %f %f %f %f %f\n", EIGN[0], z1, z2, EI[0], EI[1]);*/

	    for(j = 0; j < 2; j++)
	      {
		left[i * 4 + j * 2] = 1.0;
		right[i * 4 + j * 2] = 1.0;
				
		left[i * 4 + j * 2 + 1]  = d1 * EI[j];
		right[i * 4 + j * 2 + 1] = d2 * EI[j];
		
		/*printf("RR %f %f %f %f\n", 	left[i * 4 + j * 2], 	right[i * 4 + j * 2],	left[i * 4 + j * 2 + 1],	right[i * 4 + j * 2 + 1]);*/
	      } 
	  }
      }
      break;
    case DNA_DATA:
      {
	double d1[3], d2[3];
	   
	for(i = 0; i < numberOfCategories; i++)
	  {		                      
	    for(j = 0; j < 3; j++)
	      {
		d1[j] = exp(rptr[i] * EIGN[j] * z1);
		d2[j] = exp(rptr[i] * EIGN[j] * z2);
	      }                   
	    
	    for(j = 0; j < 4; j++)
	      {
		left[i * 16 + j * 4] = 1.0;
		right[i * 16 + j * 4] = 1.0;
		
		for(k = 0; k < 3; k++)
		  {
		    left[i * 16 + j * 4 + k + 1]  = d1[k] * EI[3 * j + k];
		    right[i * 16 + j * 4 + k + 1] = d2[k] * EI[3 * j + k];
		  }
	      } 
	  }
      }
      break;
    case SECONDARY_DATA:
      {
	double lz1[15], lz2[15], d1[15], d2[15];
	
	for(i = 0; i < 15; i++)
	  {
	    lz1[i] = EIGN[i] * z1;
	    lz2[i] = EIGN[i] * z2;  
	  }                
       
	for(i = 0; i < numberOfCategories; i++)
	  {	     
	    for(j = 0; j < 15; j++)
	      {
		d1[j] = exp (rptr[i] * lz1[j]);
		d2[j] = exp (rptr[i] * lz2[j]);
	      }        	      	 
	   	   	 
	    for(j = 0; j < 16; j++)
	      {
		left[256 * i  + 16 * j] = 1.0;	  	      
		right[256 * i + 16 * j] = 1.0;
	 
		for(k = 1; k < 16; k++)
		  {
		    left[256 * i + 16 * j + k]  = d1[k-1] * EI[15 * j + (k-1)];	  	      
		    right[256 * i + 16 * j + k] = d2[k-1] * EI[15 * j + (k-1)];
		  }	  
	      }            	              
	  }
      }
      break;
    case SECONDARY_DATA_6:
      {
	double lz1[5], lz2[5], d1[5], d2[5];
	
	for(i = 0; i < 5; i++)
	  {
	    lz1[i] = EIGN[i] * z1;
	    lz2[i] = EIGN[i] * z2;  
	  }                
       
	for(i = 0; i < numberOfCategories; i++)
	  {	     
	    for(j = 0; j < 5; j++)
	      {
		d1[j] = exp (rptr[i] * lz1[j]);
		d2[j] = exp (rptr[i] * lz2[j]);
	      }        	      	 
	   	   	 
	    for(j = 0; j < 6; j++)
	      {
		left[36 * i  + 6 * j] = 1.0;	  	      
		right[36 * i + 6 * j] = 1.0;
	 
		for(k = 1; k < 6; k++)
		  {
		    left[36 * i + 6 * j + k]  = d1[k-1] * EI[5 * j + (k-1)];	  	      
		    right[36 * i + 6 * j + k] = d2[k-1] * EI[5 * j + (k-1)];
		  }	  
	      }            	              
	  }
      }
      break;
    case SECONDARY_DATA_7:
      {
	double lz1[6], lz2[6], d1[6], d2[6];
	
	for(i = 0; i < 6; i++)
	  {
	    lz1[i] = EIGN[i] * z1;
	    lz2[i] = EIGN[i] * z2;  
	  }                
       
	for(i = 0; i < numberOfCategories; i++)
	  {	     
	    for(j = 0; j < 6; j++)
	      {
		d1[j] = exp (rptr[i] * lz1[j]);
		d2[j] = exp (rptr[i] * lz2[j]);
	      }        	      	 
	   	   	 
	    for(j = 0; j < 7; j++)
	      {
		left[49 * i  + 7 * j] = 1.0;	  	      
		right[49 * i + 7 * j] = 1.0;
	 
		for(k = 1; k < 7; k++)
		  {
		    left[49 * i + 7 * j + k]  = d1[k-1] * EI[6 * j + (k-1)];	  	      
		    right[49 * i + 7 * j + k] = d2[k-1] * EI[6 * j + (k-1)];
		  }	  
	      }            	              
	  }
      }
      break;
    case AA_DATA:
      {
	double lz1[19], lz2[19], d1[19], d2[19];

	for(i = 0; i < 19; i++)
	  {
	    lz1[i] = EIGN[i] * z1;
	    lz2[i] = EIGN[i] * z2;  
	  }                
       
	for(i = 0; i < numberOfCategories; i++)
	  {	     
	    for(j = 0; j < 19; j++)
	      {
		d1[j] = exp (rptr[i] * lz1[j]);
		d2[j] = exp (rptr[i] * lz2[j]);
	      }        	      	 
	   	   	 
	    for(j = 0; j < 20; j++)
	      {
		left[400 * i  + 20 * j] = 1.0;	  	      
		right[400 * i + 20 * j] = 1.0;
	 
		for(k = 1; k < 20; k++)
		  {
		    left[400 * i + 20 * j + k]  = d1[k-1] * EI[19 * j + (k-1)];	  	      
		    right[400 * i + 20 * j + k] = d2[k-1] * EI[19 * j + (k-1)];
		  }	  
	      }            	              
	  }
      }
      break;
    default:
      assert(0);
    }
  
}

static void newviewGTRCAT_BINARY( int tipCase,  double *EV,  int *cptr, 
				  double *x1_start,  double *x2_start,  double *x3_start,  double *tipVector,
				  int    *ex1,  int *ex2,  int *ex3, unsigned char *tipX1, unsigned char *tipX2,
				  int n,  double *left, double *right)
{         
  double  
    *le,
    *ri,
    *x1, *x2, *x3;
  double  
    ump_x1, ump_x2, x1px2[2];  
  int i, j, k, scale;   
  
  switch(tipCase)
    {
    case TIP_TIP:
      {	  	      	         	      	
	for (i = 0; i < n; i++) 
	  {		   	    	   	   	     	   
	    x1 = &(tipVector[2 * tipX1[i]]);
	    x2 = &(tipVector[2 * tipX2[i]]);	  
	    x3 = &x3_start[2 * i];
	    /*printf("%f %f %f %f\n", x1[0], x1[1], x2[0], x2[1]);*/
	    
	    le =  &left[cptr[i] * 4];
	    ri =  &right[cptr[i] * 4];	 	  
	    
	    for(j = 0; j < 2; j++)
	      {
		ump_x1 = 0.0;
		ump_x2 = 0.0;
		for(k = 0; k < 2; k++)
		  {
		    ump_x1 += x1[k] * le[j * 2 + k];
		    ump_x2 += x2[k] * ri[j * 2 + k];
		  }
		x1px2[j] = ump_x1 * ump_x2;	     
	      }
	    
	    for(j = 0; j < 2; j++)
	      x3[j] = 0.0;
	    
	    for(j = 0; j < 2; j++)                 
	      for(k = 0; k < 2; k++)	
		x3[k] += x1px2[j] * EV[j * 2 + k];	    
	    
	    

	    ex3[i] = 0;		  	    	   	    	      
	  }          
      }
      break;
    case TIP_INNER:
      {		                                
	for (i = 0; i < n; i++) 
	  {		     		      	      
	    x1 = &(tipVector[2 * tipX1[i]]);  
	    x2 = &x2_start[2 * i];
	    x3 = &x3_start[2 * i];    
	    /*printf("%f %f %f %f\n", x1[0], x1[1], x2[0], x2[1]);*/
	    le =  &left[cptr[i] * 4];
	    ri =  &right[cptr[i] * 4];	 
	    
	    for(j = 0; j < 2; j++)
	      {
		ump_x1 = 0.0;
		ump_x2 = 0.0;
		for(k = 0; k < 2; k++)
		  {
		    ump_x1 += x1[k] * le[j * 2 + k];
		    ump_x2 += x2[k] * ri[j * 2 + k];
		  }
		x1px2[j] = ump_x1 * ump_x2;
	      }
	    
	    for(j = 0; j < 2; j++)
	      x3[j] = 0.0;
	    
	    for(j = 0; j < 2; j++)          
	      for(k = 0; k < 2; k++)	
		x3[k] +=  x1px2[j] *  EV[2 * j + k];
	    
	    ex3[i] = ex2[i];

	    scale = 1;
	    for(j = 0; j < 2 && scale; j++)
	      scale = (x3[j] < minlikelihood && x3[j] > minusminlikelihood);
	    	    
	    if(scale)
	      {	 
		for(j = 0; j < 2; j++)
		  x3[j] *= twotothe256;      	      		
		ex3[i]  += 1;			   
	      }	      	      	
	  }     	            
      }
      break;
    case INNER_INNER:     
      for (i = 0; i < n; i++) 
	{		     	    
	  x1 = &x1_start[2 * i];
	  x2 = &x2_start[2 * i];
	  x3 = &x3_start[2 * i];
	  
	  le = &left[cptr[i] * 4];
	  ri = &right[cptr[i] * 4];	 

	  for(j = 0; j < 2; j++)
	    {
	      ump_x1 = 0.0;
	      ump_x2 = 0.0;
	      for(k = 0; k < 2; k++)
		{
		  ump_x1 += x1[k] * le[j * 2 + k];
		  ump_x2 += x2[k] * ri[j * 2 + k];
		}
	      x1px2[j] = ump_x1 * ump_x2;
	    }

	  for(j = 0; j < 2; j++)
	    x3[j] = 0.0;

	  for(j = 0; j < 2; j++)          
	    for(k = 0; k < 2; k++)	
	      x3[k] +=  x1px2[j] *  EV[2 * j + k];
	  
	  ex3[i] = ex1[i] + ex2[i];
	  
	  scale = 1;
	  for(j = 0; j < 2 && scale; j++)
	    scale = (x3[j] < minlikelihood && x3[j] > minusminlikelihood);
	  
	  if(scale)
	    {	 
	      for(j = 0; j < 2; j++)
		x3[j] *= twotothe256;	    
	      ex3[i] += 1;
	    }	      	          	    
	}     	  
      break;
    default:
      assert(0);
    }  
}
  
static void newviewGTRGAMMA_BINARY(int tipCase,
				   double *x1_start, double *x2_start, double *x3_start,
				   double *EV, double *tipVector,
				   int    *ex1, int *ex2, int *ex3, unsigned char *tipX1, unsigned char *tipX2,
				   const int n, double *left, double *right
				   )
{      
  double  
    *x1, *x2, *x3;      
  double
    ump_x1,
    ump_x2, 
    x1px2[4];
  int i, j, k, l, scale;  
        

  /* C-OPT figure out if we are at an inner node who has two tips/leaves 
     as descendants TIP_TIP, a tip and another inner node as descendant 
     TIP_INNER, or two inner nodes as descendants INNER_INNER */

  switch(tipCase)
    {
    case TIP_TIP:
      {      	
	for (i = 0; i < n; i++) 
	  { 
	    x1 = &(tipVector[2 * tipX1[i]]);		     	    	    
	    x2 = &(tipVector[2 * tipX2[i]]);
	    x3 = &x3_start[i * 8];

	    for(j = 0; j < 8; j++)
	      x3[j] = 0.0;

	    for (j = 0; j < 4; j++)
	      {	    	     
		for (k = 0; k < 2; k++)
		  {
		    ump_x1 = 0.0;
		    ump_x2 = 0.0;		 		
		    
		    for (l=0; l < 2; l++)
		      {
			ump_x1 += x1[l] * left[ j*4 + k*2 + l];
			ump_x2 += x2[l] * right[j*4 + k*2 + l];
		      }
		    
		    x1px2[k] = ump_x1 * ump_x2;				 
		  }
		
		for(k = 0; k < 2; k++)
		  for (l = 0; l < 2; l++)
		    x3[j * 2 + l] +=  x1px2[k] * EV[2 * k + l];
		
	      }  	    	       	    
	    
	    ex3[i] = 0;	
	  }                
      }
      break;      
    case TIP_INNER:
      {	 	 	
	 for (i = 0; i < n; i++) 
	   {	    	    	 
	     x1 = &(tipVector[2 * tipX1[i]]);
	     x2 = &x2_start[i * 8];
	     x3 = &x3_start[i * 8];

	     for(j = 0; j < 8; j++)
	       x3[j] = 0.0;

	     for (j = 0; j < 4; j++)
	       {	    	     
		 for (k = 0; k < 2; k++)
		   {
		     ump_x1 = 0.0;
		     ump_x2 = 0.0;		 		
		     
		     for (l=0; l < 2; l++)
		       {
			 ump_x1 += x1[l] * left[ j*4 + k*2 + l];
			 ump_x2 += x2[j*2 + l] * right[j*4 + k*2 + l];
		       }
		     
		     x1px2[k] = ump_x1 * ump_x2;				 
		   }
		 
		 for(k = 0; k < 2; k++)
		   for (l = 0; l < 2; l++)
		     x3[j * 2 + l] +=  x1px2[k] * EV[2 * k + l];
		 
	       }  	    	 	     

	     ex3[i] = ex2[i];

	     scale = 1;
	     for(l = 0; scale && (l < 8); l++)
	       scale = (ABS(x3[l]) <  minlikelihood);
	    
	     if(scale)
	       {
		 for (l=0; l < 8; l++)
		   x3[l] *= twotothe256;
		 ex3[i] += 1;
	       }
	         	   
	   }
      }
      break;
    case INNER_INNER:
      
      /* C-OPT here we don't do any pre-computations 
	 This should be the most compute intensive loop of the three
	 cases here. If we have one or two tips as descendants 
	 we can take a couple of shortcuts */
      
      
     for (i = 0; i < n; i++) 
       {		
	 x1 = &x1_start[i * 8];
	 x2 = &x2_start[i * 8];
	 x3 = &x3_start[i * 8];
	 
	 for(j = 0; j < 8; j++)
	   x3[j] = 0.0;

	 for (j = 0; j < 4; j++)
	   {	    	     
	     for (k = 0; k < 2; k++)
	       {
		 ump_x1 = 0.0;
		 ump_x2 = 0.0;		 		
		
		 for (l=0; l < 2; l++)
		   {
		     ump_x1 += x1[j*2 + l] * left[ j*4 + k*2 + l];
		     ump_x2 += x2[j*2 + l] * right[j*4 + k*2 + l];
		   }
		 
		 x1px2[k] = ump_x1 * ump_x2;				 
	       }
	     
	     for(k = 0; k < 2; k++)
	       for (l = 0; l < 2; l++)
		 x3[j * 2 + l] +=  x1px2[k] * EV[2 * k + l];

	   }

	 ex3[i] = ex1[i] + ex2[i];
	 scale = 1;
	 for(l = 0; scale && (l < 8); l++)
	   scale = (ABS(x3[l]) <  minlikelihood);
	
	 if(scale)
	   {
	     for (l=0; l<8; l++)
	       x3[l] *= twotothe256;
	     ex3[i] += 1;
	   }
       }
     break;
    default:
      assert(0);
    }          
}




static void newviewGTRCAT( int tipCase,  double *EV,  int *cptr, 
			   double *x1_start,  double *x2_start,  double *x3_start,  double *tipVector,
			   int    *ex1,  int *ex2,  int *ex3, unsigned char *tipX1, unsigned char *tipX2,
			   int n,  double *left, double *right)
{         
  double  
    *le,
    *ri,
    *x1, *x2, *x3;
  double  
    ump_x1, ump_x2, x1px2[4];  
  int i, j, k, scale;   
  
  switch(tipCase)
    {
    case TIP_TIP:
      {	  	      	         	      	
	for (i = 0; i < n; i++) 
	  {		   	    	   	   	     	   
	    x1 = &(tipVector[4 * tipX1[i]]);
	    x2 = &(tipVector[4 * tipX2[i]]);	  
	    x3 = &x3_start[4 * i];	 
	    
	    le =  &left[cptr[i] * 16];
	    ri =  &right[cptr[i] * 16];	 	  
	    
	    for(j = 0; j < 4; j++)
	      {
		ump_x1 = 0.0;
		ump_x2 = 0.0;
		for(k = 0; k < 4; k++)
		  {
		    ump_x1 += x1[k] * le[j * 4 + k];
		    ump_x2 += x2[k] * ri[j * 4 + k];
		  }
		x1px2[j] = ump_x1 * ump_x2;	     
	      }
	    
	    for(j = 0; j < 4; j++)
	      x3[j] = 0.0;
	    
	    for(j = 0; j < 4; j++)                 
	      for(k = 0; k < 4; k++)	
		x3[k] += x1px2[j] * EV[j * 4 + k];	    
	    
	    ex3[i] = 0;		  	    	   	    	      
	  }          
      }
      break;
    case TIP_INNER:
      {		                                
	for (i = 0; i < n; i++) 
	  {		     		      	      
	    x1 = &(tipVector[4 * tipX1[i]]);  
	    x2 = &x2_start[4 * i];
	    x3 = &x3_start[4 * i];    
	    
	    le =  &left[cptr[i] * 16];
	    ri =  &right[cptr[i] * 16];	 
	    
	    for(j = 0; j < 4; j++)
	      {
		ump_x1 = 0.0;
		ump_x2 = 0.0;
		for(k = 0; k < 4; k++)
		  {
		    ump_x1 += x1[k] * le[j * 4 + k];
		    ump_x2 += x2[k] * ri[j * 4 + k];
		  }
		x1px2[j] = ump_x1 * ump_x2;
	      }
	    
	    for(j = 0; j < 4; j++)
	      x3[j] = 0.0;
	    
	    for(j = 0; j < 4; j++)          
	      for(k = 0; k < 4; k++)	
		x3[k] +=  x1px2[j] *  EV[4 * j + k];
	    
	    ex3[i] = ex2[i];

	    scale = 1;
	    for(j = 0; j < 4 && scale; j++)
	      scale = (x3[j] < minlikelihood && x3[j] > minusminlikelihood);
	    	    
	    if(scale)
	      {	 	       
		for(j = 0; j < 4; j++)
		  x3[j] *= twotothe256;      	      		
		ex3[i]  += 1;			   
	      }	      	      	
	  }     	            
      }
      break;
    case INNER_INNER:     
      for (i = 0; i < n; i++) 
	{		     	    
	  x1 = &x1_start[4 * i];
	  x2 = &x2_start[4 * i];
	  x3 = &x3_start[4 * i];
	  
	  le = &left[cptr[i] * 16];
	  ri = &right[cptr[i] * 16];	 

	  for(j = 0; j < 4; j++)
	    {
	      ump_x1 = 0.0;
	      ump_x2 = 0.0;
	      for(k = 0; k < 4; k++)
		{
		  ump_x1 += x1[k] * le[j * 4 + k];
		  ump_x2 += x2[k] * ri[j * 4 + k];
		}
	      x1px2[j] = ump_x1 * ump_x2;
	    }

	  for(j = 0; j < 4; j++)
	    x3[j] = 0.0;

	  for(j = 0; j < 4; j++)          
	    for(k = 0; k < 4; k++)	
	      x3[k] +=  x1px2[j] *  EV[4 * j + k];
	  
	  ex3[i] = ex1[i] + ex2[i];
	  
	  scale = 1;
	  for(j = 0; j < 4 && scale; j++)
	    scale = (x3[j] < minlikelihood && x3[j] > minusminlikelihood);
	  
	  if(scale)
	    {	     
	      for(j = 0; j < 4; j++)
		x3[j] *= twotothe256;	    
	      ex3[i] += 1;
	    }	      	          	    
	}     	  
      break;
    default:
      assert(0);
    }  
}
  


/* C-OPT: now this is a function where the program will spend most of its time, between 
   50% and 60% of total run time, hence cache optimization should start here */


static void newviewGTRGAMMA(int tipCase,
			    double *x1_start, double *x2_start, double *x3_start,
			    double *EV, double *tipVector,
			    int    *ex1, int *ex2, int *ex3, unsigned char *tipX1, unsigned char *tipX2,
			    const int n, double *left, double *right
			    )
{      
  double  
    *x1, *x2, *x3;      
  double
    buf,
    ump_x1,
    ump_x2, 
    x1px2[4];
  int i, j, k, l, scale;  
        

  /* C-OPT figure out if we are at an inner node who has two tips/leaves 
     as descendants TIP_TIP, a tip and another inner node as descendant 
     TIP_INNER, or two inner nodes as descendants INNER_INNER */

  switch(tipCase)
    {
    case TIP_TIP:
      {      
	double *uX1, umpX1[256], *uX2, umpX2[256];	 


	/* C-OPT preparatory stuff to avoid extra computations */

	for(i = 1; i < 16; i++)
	  {		   
	    x1 = &(tipVector[i * 4]);	    
	    	    
	    for(j=0; j<4; j++)
	      for(k=0; k<4; k++)
		{
		  umpX1[i*16 + j*4 + k] = 0.0;
		  umpX2[i*16 + j*4 + k] = 0.0;

		  for (l=0; l < 4; l++)
		    {
		      umpX1[i*16 + j*4 + k] += x1[l] * left[j*16 + k*4 + l];
		      umpX2[i*16 + j*4 + k] += x1[l] * right[j*16 + k*4 + l];
		    }	 
		}
	  }	

	/* C-OPT now this is the loop over the input alignment length 
	   I would expect most of the cache misses to happen in the loop 
	   below */
	

	for (i = 0; i < n; i++) 
	  {		     	    	    
	    x3 = &x3_start[i * 16];

	    uX1 = &umpX1[16 * tipX1[i]];
	    uX2 = &umpX2[16 * tipX2[i]];	    	        	    

	    for(j = 0; j < 16; j++)
	      x3[j] = 0.0;

	    for (j = 0; j < 4; j++)
	      for (k = 0; k < 4; k++)
		{		 
		  buf = uX1[j*4 + k] * uX2[j*4 + k];

		  for (l=0; l<4; l++)
		    x3[j * 4 + l] +=  buf * EV[4 * k + l];
		}	       	    
	    
	    ex3[i] = 0;	
	  }                
      }
      break;      
    case TIP_INNER:
      {	 	 
	double *uX1, umpX1[256];		 	
	 
	/* C-OPT once again, pre-computing some values */
	
	for (i = 1; i < 16; i++)
	  {	   	  
	    x1 = &(tipVector[i*4]);
	   
	    for (j = 0; j < 4; j++)
	      for (k = 0; k < 4; k++)
		{
		  umpX1[i*16 + j*4 + k] = 0.0;
		  for (l=0; l < 4; l++)
		    umpX1[i*16 + j*4 + k] += x1[l] * left[j*16 + k*4 + l];
		}
	  }

	 	    	       	
	/* C-OPT main long and compute-intensive loop again */

	 for (i = 0; i < n; i++) 
	   {	    	    	     
	     x2 = &x2_start[i * 16];
	     x3 = &x3_start[i * 16];

	     uX1 = &umpX1[16 * tipX1[i]];

	     for(j = 0; j < 16; j++)
	       x3[j] = 0.0;

	     for (j = 0; j < 4; j++)
	       {
		 for (k = 0; k < 4; k++)
		   {
		     ump_x2 = 0.0;
		     
		     for (l=0; l<4; l++)
		       ump_x2 += x2[j*4 + l] * right[j* 16 + k*4 + l];
		     x1px2[k] = uX1[j * 4 + k] * ump_x2;
		   }
		  	
		 for(k = 0; k < 4; k++)
		   for (l=0; l<4; l++)
		     x3[j * 4 + l] +=  x1px2[k] * EV[4 * k + l];		      		   
	       }	    	    	 


	     ex3[i] = ex2[i];

	     scale = 1;
	     for(l = 0; scale && (l < 16); l++)
	       scale = (ABS(x3[l]) <  minlikelihood);
	    
	     if(scale)
	       {
		 for (l=0; l<16; l++)
		   x3[l] *= twotothe256;
		 ex3[i] += 1;
	       }
	         	   
	   }
      }
      break;
    case INNER_INNER:
      
      /* C-OPT here we don't do any pre-computations 
	 This should be the most compute intensive loop of the three
	 cases here. If we have one or two tips as descendants 
	 we can take a couple of shortcuts */
      
      
     for (i = 0; i < n; i++) 
       {		
	 x1 = &x1_start[i * 16];
	 x2 = &x2_start[i * 16];
	 x3 = &x3_start[i * 16];
	 
	 for(j = 0; j < 16; j++)
	   x3[j] = 0.0;

	 for (j = 0; j < 4; j++)
	   {	    	     
	     for (k = 0; k < 4; k++)
	       {
		 ump_x1 = 0.0;
		 ump_x2 = 0.0;		 		
		
		 for (l=0; l<4; l++)
		   {
		     ump_x1 += x1[j*4 + l] * left[j*16 + k*4 +l];
		     ump_x2 += x2[j*4 + l] * right[j*16 + k*4 +l];
		   }
		 
		 x1px2[k] = ump_x1 * ump_x2;				 
	       }
	     
	     for(k = 0; k < 4; k++)
	       for (l=0; l<4; l++)
		 x3[j * 4 + l] +=  x1px2[k] * EV[4 * k + l];

	   }

	 ex3[i] = ex1[i] + ex2[i];
	 scale = 1;
	 for(l = 0; scale && (l < 16); l++)
	   scale = (ABS(x3[l]) <  minlikelihood);
	
	 if(scale)
	   {
	     for (l=0; l<16; l++)
	       x3[l] *= twotothe256;
	     ex3[i] += 1;
	   }
       }
     break;
    default:
      assert(0);
    }          
}





static void newviewGTRCATPROT(int tipCase, double *extEV,
			      int *cptr,
			      double *x1, double *x2, double *x3, double *tipVector,				
			      int    *ex1, int *ex2, int *ex3, unsigned char *tipX1, unsigned char *tipX2,
			      int n, double *left, double *right)
{
  double  
    *le, *ri, *v, *vl, *vr;
  double      
    ump_x1, ump_x2, x1px2;
  int i, l, j, scale;	                                                    
  
  switch(tipCase)	 
    {
    case TIP_TIP:       
      {	  	 	 	   	   	        
	for (i = 0; i < n; i++) 
	  {		   	    	   	   	     	  
	    le = &left[cptr[i] * 400];	      
	    ri = &right[cptr[i] * 400];

	    vl = &(tipVector[20 * tipX1[i]]);	     	       	       
	    vr = &(tipVector[20 * tipX2[i]]);
	    v  = &x3[20 * i];

	    for(l = 0; l < 20; l++)
	      v[l] = 0.0;

	    for(l = 0; l < 20; l++)
	      {
		ump_x1 = 0.0;
		ump_x2 = 0.0;	      
		
		for(j = 0; j < 20; j++)
		  {
		    ump_x1 += vl[j] * le[l * 20 + j];			       		     
		    ump_x2 += vr[j] * ri[l * 20 + j];
		  }		
		
		x1px2 = ump_x1 * ump_x2;

		for(j = 0; j < 20; j++)
		  v[j] += x1px2 * extEV[l * 20 + j];
	      }	    	    	    	  	    	    	    	   	    	    	    	    	   	    
	    	    
	    ex3[i] = 0;		  	    	   	    	 
	  }       
      }
      break;
    case TIP_INNER:    
      {		   
	for (i = 0; i < n; i++) 
	  {	
	    le = &left[cptr[i] * 400];	    	       	   
	    ri = &right[cptr[i] * 400]; 
	    
	    vl = &(tipVector[20 * tipX1[i]]);     
	    vr = &x2[20 * i];
	    v  = &x3[20 * i];
	    
	    for(l = 0; l < 20; l++)
	      v[l] = 0.0;
	    
	    for(l = 0; l < 20; l++)
	      {		   
		ump_x1 = 0.0;
		ump_x2 = 0.0;	  
		
		for(j = 0; j < 20; j++)
		  {
		    ump_x1 += vl[j] * le[l * 20 + j];			       		     
		    ump_x2 += vr[j] * ri[l * 20 + j];
		  }
		
		x1px2 = ump_x1 * ump_x2;
		
		for(j = 0; j < 20; j++)
		  v[j] += x1px2 * extEV[l * 20 + j];		   
	      }	    	      	       	       	   
	    
	    scale = 1;
	    for(l = 0; scale && (l < 20); l++)
	      scale = ((v[l] < minlikelihood) && (v[l] > minusminlikelihood));
	    
	    ex3[i] = ex2[i];		  	      	       	       
	    
	    if(scale)
	      {			  
		for(l = 0; l < 20; l++)
		  v[l] *= twotothe256;		   
		ex3[i] += 1;
	      }	            	       
	  }     	 	  
      }
      break;
    case INNER_INNER:
      for(i = 0; i < n; i++) 
	{		 	  	 	  	 
	  le = &left[cptr[i] * 400];	    	       	   
	  ri = &right[cptr[i] * 400]; 

	  vl = &x1[20 * i];     
	  vr = &x2[20 * i];
	  v = &x3[20 * i];

	  for(l = 0; l < 20; l++)
	    v[l] = 0.0;

	  for(l = 0; l < 20; l++)
	    {	     
	      ump_x1 = 0.0;
	      ump_x2 = 0.0;

	      for(j = 0; j < 20; j++)
		{
		  ump_x1 += vl[j] * le[l * 20 + j];			       		     
		  ump_x2 += vr[j] * ri[l * 20 + j];
		}
	      
	      x1px2 =  ump_x1 * ump_x2;

	      for(j = 0; j < 20; j++)
		v[j] += x1px2 * extEV[l * 20 + j];	     
	    }	    	 	  	  
	        
	   scale = 1;
	   for(l = 0; scale && (l < 20); l++)
	     scale = ((v[l] < minlikelihood) && (v[l] > minusminlikelihood));

	   ex3[i] = ex1[i] + ex2[i];		  
	   
	   if(scale)	
	     {	 	       
	       for(l = 0; l < 20; l++)
		 v[l] *= twotothe256;
	       ex3[i] += 1;
	     }	            
	} 
      break;
    default:
      assert(0);
    }   
}


static void newviewGTRCATSECONDARY(int tipCase, double *extEV,
				   int *cptr,
				   double *x1, double *x2, double *x3, double *tipVector,				
				   int    *ex1, int *ex2, int *ex3, unsigned char *tipX1, unsigned char *tipX2,
				   int n, double *left, double *right)
{
  double  
    *le, *ri, *v, *vl, *vr;
  double      
    ump_x1, ump_x2, x1px2;
  int i, l, j, scale;	                                                    
  
  switch(tipCase)	 
    {
    case TIP_TIP:       
      {	  	 	 	   	   	        
	for (i = 0; i < n; i++) 
	  {		   	    	   	   	     	  
	    le = &left[cptr[i] * 256];	      
	    ri = &right[cptr[i] * 256];

	    vl = &(tipVector[16 * tipX1[i]]);	     	       	       
	    vr = &(tipVector[16 * tipX2[i]]);
	    v  = &x3[16 * i];

	    for(l = 0; l < 16; l++)
	      v[l] = 0.0;

	    for(l = 0; l < 16; l++)
	      {
		ump_x1 = 0.0;
		ump_x2 = 0.0;	      
		
		for(j = 0; j < 16; j++)
		  {
		    ump_x1 += vl[j] * le[l * 16 + j];			       		     
		    ump_x2 += vr[j] * ri[l * 16 + j];
		  }		
		
		x1px2 = ump_x1 * ump_x2;

		for(j = 0; j < 16; j++)
		  v[j] += x1px2 * extEV[l * 16 + j];
	      }	    	    	    	  	    	    	    	   	    	    	    	    	   	    
	    	    
	    ex3[i] = 0;		  	    	   	    	 
	  }       
      }
      break;
    case TIP_INNER:    
      {		   
	for (i = 0; i < n; i++) 
	  {	
	    le = &left[cptr[i] * 256];	    	       	   
	    ri = &right[cptr[i] * 256]; 
	    
	    vl = &(tipVector[16 * tipX1[i]]);     
	    vr = &x2[16 * i];
	    v  = &x3[16 * i];
	    
	    for(l = 0; l < 16; l++)
	      v[l] = 0.0;
	    
	    for(l = 0; l < 16; l++)
	      {		   
		ump_x1 = 0.0;
		ump_x2 = 0.0;	  
		
		for(j = 0; j < 16; j++)
		  {
		    ump_x1 += vl[j] * le[l * 16 + j];			       		     
		    ump_x2 += vr[j] * ri[l * 16 + j];
		  }
		
		x1px2 = ump_x1 * ump_x2;
		
		for(j = 0; j < 16; j++)
		  v[j] += x1px2 * extEV[l * 16 + j];		   
	      }	    	      	       	       	   
	    
	    scale = 1;
	    for(l = 0; scale && (l < 16); l++)
	      scale = ((v[l] < minlikelihood) && (v[l] > minusminlikelihood));
	    
	    ex3[i] = ex2[i];		  	      	       	       
	    
	    if(scale)
	      {			  
		for(l = 0; l < 16; l++)
		  v[l] *= twotothe256;		   
		ex3[i] += 1;
	      }	            	       
	  }     	 	  
      }
      break;
    case INNER_INNER:
      for(i = 0; i < n; i++) 
	{		 	  	 	  	 
	  le = &left[cptr[i] * 256];	    	       	   
	  ri = &right[cptr[i] * 256]; 

	  vl = &x1[16 * i];     
	  vr = &x2[16 * i];
	  v = &x3[16 * i];

	  for(l = 0; l < 16; l++)
	    v[l] = 0.0;

	  for(l = 0; l < 16; l++)
	    {	     
	      ump_x1 = 0.0;
	      ump_x2 = 0.0;

	      for(j = 0; j < 16; j++)
		{
		  ump_x1 += vl[j] * le[l * 16 + j];			       		     
		  ump_x2 += vr[j] * ri[l * 16 + j];
		}
	      
	      x1px2 =  ump_x1 * ump_x2;

	      for(j = 0; j < 16; j++)
		v[j] += x1px2 * extEV[l * 16 + j];	     
	    }	    	 	  	  
	        
	   scale = 1;
	   for(l = 0; scale && (l < 16); l++)
	     scale = ((v[l] < minlikelihood) && (v[l] > minusminlikelihood));

	   ex3[i] = ex1[i] + ex2[i];		  
	   
	   if(scale)	
	     {	 	       
	       for(l = 0; l < 16; l++)
		 v[l] *= twotothe256;
	       ex3[i] += 1;
	     }	            
	} 
      break;
    default:
      assert(0);
    }   
}

static void newviewGTRCATSECONDARY_6(int tipCase, double *extEV,
				   int *cptr,
				   double *x1, double *x2, double *x3, double *tipVector,				
				   int    *ex1, int *ex2, int *ex3, unsigned char *tipX1, unsigned char *tipX2,
				   int n, double *left, double *right)
{
  double  
    *le, *ri, *v, *vl, *vr;
  double      
    ump_x1, ump_x2, x1px2;
  int i, l, j, scale;	                                                    
  
  switch(tipCase)	 
    {
    case TIP_TIP:       
      {	  	 	 	   	   	        
	for (i = 0; i < n; i++) 
	  {		   	    	   	   	     	  
	    le = &left[cptr[i] * 36];	      
	    ri = &right[cptr[i] * 36];

	    vl = &(tipVector[6 * tipX1[i]]);	     	       	       
	    vr = &(tipVector[6 * tipX2[i]]);
	    v  = &x3[6 * i];

	    for(l = 0; l < 6; l++)
	      v[l] = 0.0;

	    for(l = 0; l < 6; l++)
	      {
		ump_x1 = 0.0;
		ump_x2 = 0.0;	      
		
		for(j = 0; j < 6; j++)
		  {
		    ump_x1 += vl[j] * le[l * 6 + j];			       		     
		    ump_x2 += vr[j] * ri[l * 6 + j];
		  }		
		
		x1px2 = ump_x1 * ump_x2;

		for(j = 0; j < 6; j++)
		  v[j] += x1px2 * extEV[l * 6 + j];
	      }	    	    	    	  	    	    	    	   	    	    	    	    	   	    
	    	    
	    ex3[i] = 0;		  	    	   	    	 
	  }       
      }
      break;
    case TIP_INNER:    
      {		   
	for (i = 0; i < n; i++) 
	  {	
	    le = &left[cptr[i] * 36];	    	       	   
	    ri = &right[cptr[i] * 36]; 
	    
	    vl = &(tipVector[6 * tipX1[i]]);     
	    vr = &x2[6 * i];
	    v  = &x3[6 * i];
	    
	    for(l = 0; l < 6; l++)
	      v[l] = 0.0;
	    
	    for(l = 0; l < 6; l++)
	      {		   
		ump_x1 = 0.0;
		ump_x2 = 0.0;	  
		
		for(j = 0; j < 6; j++)
		  {
		    ump_x1 += vl[j] * le[l * 6 + j];			       		     
		    ump_x2 += vr[j] * ri[l * 6 + j];
		  }
		
		x1px2 = ump_x1 * ump_x2;
		
		for(j = 0; j < 6; j++)
		  v[j] += x1px2 * extEV[l * 6 + j];		   
	      }	    	      	       	       	   
	    
	    scale = 1;
	    for(l = 0; scale && (l < 6); l++)
	      scale = ((v[l] < minlikelihood) && (v[l] > minusminlikelihood));
	    
	    ex3[i] = ex2[i];		  	      	       	       
	    
	    if(scale)
	      {			  
		for(l = 0; l < 6; l++)
		  v[l] *= twotothe256;		   
		ex3[i] += 1;
	      }	            	       
	  }     	 	  
      }
      break;
    case INNER_INNER:
      for(i = 0; i < n; i++) 
	{		 	  	 	  	 
	  le = &left[cptr[i] * 36];	    	       	   
	  ri = &right[cptr[i] * 36]; 

	  vl = &x1[6 * i];     
	  vr = &x2[6 * i];
	  v = &x3[6 * i];

	  for(l = 0; l < 6; l++)
	    v[l] = 0.0;

	  for(l = 0; l < 6; l++)
	    {	     
	      ump_x1 = 0.0;
	      ump_x2 = 0.0;

	      for(j = 0; j < 6; j++)
		{
		  ump_x1 += vl[j] * le[l * 6 + j];			       		     
		  ump_x2 += vr[j] * ri[l * 6 + j];
		}
	      
	      x1px2 =  ump_x1 * ump_x2;

	      for(j = 0; j < 6; j++)
		v[j] += x1px2 * extEV[l * 6 + j];	     
	    }	    	 	  	  
	        
	   scale = 1;
	   for(l = 0; scale && (l < 6); l++)
	     scale = ((v[l] < minlikelihood) && (v[l] > minusminlikelihood));

	   ex3[i] = ex1[i] + ex2[i];		  
	   
	   if(scale)	
	     {	 	       
	       for(l = 0; l < 6; l++)
		 v[l] *= twotothe256;
	       ex3[i] += 1;
	     }	            
	} 
      break;
    default:
      assert(0);
    }   
}

static void newviewGTRCATSECONDARY_7(int tipCase, double *extEV,
				     int *cptr,
				     double *x1, double *x2, double *x3, double *tipVector,				
				     int    *ex1, int *ex2, int *ex3, unsigned char *tipX1, unsigned char *tipX2,
				     int n, double *left, double *right)
{
  double  
    *le, *ri, *v, *vl, *vr;
  double      
    ump_x1, ump_x2, x1px2;
  int i, l, j, scale;	                                                    
  
  switch(tipCase)	 
    {
    case TIP_TIP:       
      {	  	 	 	   	   	        
	for (i = 0; i < n; i++) 
	  {		   	    	   	   	     	  
	    le = &left[cptr[i] * 49];	      
	    ri = &right[cptr[i] * 49];

	    vl = &(tipVector[7 * tipX1[i]]);	     	       	       
	    vr = &(tipVector[7 * tipX2[i]]);
	    v  = &x3[7 * i];

	    for(l = 0; l < 7; l++)
	      v[l] = 0.0;

	    for(l = 0; l < 7; l++)
	      {
		ump_x1 = 0.0;
		ump_x2 = 0.0;	      
		
		for(j = 0; j < 7; j++)
		  {
		    ump_x1 += vl[j] * le[l * 7 + j];			       		     
		    ump_x2 += vr[j] * ri[l * 7 + j];
		  }		
		
		x1px2 = ump_x1 * ump_x2;

		for(j = 0; j < 7; j++)
		  v[j] += x1px2 * extEV[l * 7 + j];
	      }	    	    	    	  	    	    	    	   	    	    	    	    	   	    
	    	    
	    ex3[i] = 0;		  	    	   	    	 
	  }       
      }
      break;
    case TIP_INNER:    
      {		   
	for (i = 0; i < n; i++) 
	  {	
	    le = &left[cptr[i] * 49];	    	       	   
	    ri = &right[cptr[i] * 49]; 
	    
	    vl = &(tipVector[7 * tipX1[i]]);     
	    vr = &x2[7 * i];
	    v  = &x3[7 * i];
	    
	    for(l = 0; l < 7; l++)
	      v[l] = 0.0;
	    
	    for(l = 0; l < 7; l++)
	      {		   
		ump_x1 = 0.0;
		ump_x2 = 0.0;	  
		
		for(j = 0; j < 7; j++)
		  {
		    ump_x1 += vl[j] * le[l * 7 + j];			       		     
		    ump_x2 += vr[j] * ri[l * 7 + j];
		  }
		
		x1px2 = ump_x1 * ump_x2;
		
		for(j = 0; j < 7; j++)
		  v[j] += x1px2 * extEV[l * 7 + j];		   
	      }	    	      	       	       	   
	    
	    scale = 1;
	    for(l = 0; scale && (l < 7); l++)
	      scale = ((v[l] < minlikelihood) && (v[l] > minusminlikelihood));
	    
	    ex3[i] = ex2[i];		  	      	       	       
	    
	    if(scale)
	      {			  
		for(l = 0; l < 7; l++)
		  v[l] *= twotothe256;		   
		ex3[i] += 1;
	      }	            	       
	  }     	 	  
      }
      break;
    case INNER_INNER:
      for(i = 0; i < n; i++) 
	{		 	  	 	  	 
	  le = &left[cptr[i] * 49];	    	       	   
	  ri = &right[cptr[i] * 49]; 

	  vl = &x1[7 * i];     
	  vr = &x2[7 * i];
	  v = &x3[7 * i];

	  for(l = 0; l < 7; l++)
	    v[l] = 0.0;

	  for(l = 0; l < 7; l++)
	    {	     
	      ump_x1 = 0.0;
	      ump_x2 = 0.0;

	      for(j = 0; j < 7; j++)
		{
		  ump_x1 += vl[j] * le[l * 7 + j];			       		     
		  ump_x2 += vr[j] * ri[l * 7 + j];
		}
	      
	      x1px2 =  ump_x1 * ump_x2;

	      for(j = 0; j < 7; j++)
		v[j] += x1px2 * extEV[l * 7 + j];	     
	    }	    	 	  	  
	        
	   scale = 1;
	   for(l = 0; scale && (l < 7); l++)
	     scale = ((v[l] < minlikelihood) && (v[l] > minusminlikelihood));

	   ex3[i] = ex1[i] + ex2[i];		  
	   
	   if(scale)	
	     {	 	       
	       for(l = 0; l < 7; l++)
		 v[l] *= twotothe256;
	       ex3[i] += 1;
	     }	            
	} 
      break;
    default:
      assert(0);
    }   
}


static void newviewGTRGAMMAPROT(int tipCase,
				double *x1, double *x2, double *x3, double *extEV, double *tipVector,
				int    *ex1, int *ex2, int *ex3, unsigned char *tipX1, unsigned char *tipX2,
				int n, double *left, double *right)
{  
  double  *uX1, *uX2, *v;
  double x1px2;    
  int  i, j, l, k, scale;   
  double *vl, *vr, al, ar;                

 
  
  switch(tipCase)
    {      
    case TIP_TIP:
      {      
	double umpX1[1840], umpX2[1840];
	
	for(i = 0; i < 23; i++)
	  {	    	     
	    v = &(tipVector[20 * i]);	     	    
	    
	    for(k = 0; k < 80; k++)
	      {	
		umpX1[80 * i + k] = 0.0;
		umpX2[80 * i + k] = 0.0;

		for(l = 0; l < 20; l++)
		  {
		    umpX1[80 * i + k] +=  v[l] *  left[k * 20 + l];
		    umpX2[80 * i + k] +=  v[l] * right[k * 20 + l];
		  }
	      }	    	   
	  }	
	
	for(i = 0; i < n; i++) 
	  {		     
	    uX1 = &umpX1[80 * tipX1[i]];
	    uX2 = &umpX2[80 * tipX2[i]];	    

	    for(j = 0; j < 4; j++)	      
	      {
		v = &x3[i * 80 + j * 20];

		for(k = 0; k < 20; k++)
		  v[k] = 0.0;
		
		for(k = 0; k < 20; k++)
		  {		   		
		    x1px2 = uX1[j * 20 + k] * uX2[j * 20 + k];
		    for(l = 0; l < 20; l++)
		      v[l] += x1px2 * extEV[20 * k + l];
		  }				
	      }

	    ex3[i] = 0;		  		    	     	            
	  }         	
      }
      break;
    case TIP_INNER:
      {	    
	double umpX1[1840], ump_x2[20];	

	 	
	for(i = 0; i < 23; i++)
	  {	    	     
	    v = &(tipVector[20 * i]);	     	    
	    
	    for(k = 0; k < 80; k++)
	      {	
		umpX1[80 * i + k] = 0.0;
		
		for(l = 0; l < 20; l++)		  
		  umpX1[80 * i + k] +=  v[l] * left[k * 20 + l];

	      }	    	   
	  }	  	       	

	for (i = 0; i < n; i++) 
	  {		   	   
	    uX1 = &umpX1[80 * tipX1[i]];			    	    
	    
	    for(k = 0; k < 4; k++)
	      {
		v = &(x2[80 * i + k * 20]);	    
		for(l = 0; l < 20; l++)
		  {			   
		    ump_x2[l] = 0.0;

		    for(j = 0; j < 20; j++)
		      ump_x2[l] += v[j] * right[k * 400 + l * 20 + j];
		  }
		 	     	       		
		v = &(x3[80 * i + 20 * k]);			

		for(l = 0; l < 20; l++)
		  v[l] = 0;

		for(l = 0; l < 20; l++)
		  {
		    x1px2 = uX1[k * 20 + l]  * ump_x2[l];
		    for(j = 0; j < 20; j++)
		      v[j] += x1px2 * extEV[l * 20  + j];		    
		  }
	      }	
	    
	    ex3[i] = ex2[i];		
	    v = &x3[80 * i];
	    scale = 1;
	    for(l = 0; scale && (l < 80); l++)
	      scale = (ABS(v[l]) <  minlikelihood);
	    
	    if (scale)	         
	      {	     	    	 
		for(l = 0; l < 80; l++)
		  v[l] *= twotothe256;		   
		ex3[i] += 1;  
	      }	 	           	 
	  }
      }
      break;
    case INNER_INNER:	      
      for (i = 0; i < n; i++) 
       {		 	
	 for(k = 0; k < 4; k++)
	   {	   
	     vl = &(x1[80 * i + 20 * k]);
	     vr = &(x2[80 * i + 20 * k]);
	     v =  &(x3[80 * i + 20 * k]);
	     	    
	     for(l = 0; l < 20; l++)	       
	       v[l] = 0;		 
	       
	     for(l = 0; l < 20; l++)
	       {		
		 al = 0.0;
		 ar = 0.0;
		 for(j = 0; j < 20; j++)
		   {
		     al += vl[j] * left[k * 400 + l * 20 + j];
		     ar += vr[j] * right[k * 400 + l * 20 + j];
		   }
		 
		 x1px2 = al * ar;
		 for(j = 0; j < 20; j++)
		   v[j] += x1px2 * extEV[20 * l + j];		
	       }
	   }
	 
	 ex3[i] = ex1[i] + ex2[i];
	 v = &(x3[80 * i]);
	 scale = 1;
	 for(l = 0; scale && (l < 80); l++)
	   scale = ((ABS(v[l]) <  minlikelihood));
	 
	 if (scale)	         
	   {	     	    	 
	     for(l = 0; l < 80; l++)
	       v[l] *= twotothe256;		   
	     ex3[i] += 1;  
	   }		 	 	 	   
       }
      break;
    default:
      assert(0);
    }  
}


static void newviewGTRGAMMASECONDARY(int tipCase,
				     double *x1, double *x2, double *x3, double *extEV, double *tipVector,
				     int    *ex1, int *ex2, int *ex3, unsigned char *tipX1, unsigned char *tipX2,
				     int n, double *left, double *right)
{  
  double  *v;
  double x1px2;    
  int  i, j, l, k, scale;   
  double *vl, *vr, al, ar;                
   
  switch(tipCase)
    {      
    case TIP_TIP:
      {      	
	for(i = 0; i < n; i++) 
	  {		     	    	    
	    for(k = 0; k < 4; k++)
	      {	   
		vl = &(tipVector[16 * tipX1[i]]);
		vr = &(tipVector[16 * tipX2[i]]);
		v =  &(x3[64 * i + 16 * k]);
		
		for(l = 0; l < 16; l++)	       
		  v[l] = 0;		 
		
		for(l = 0; l < 16; l++)
		  {		
		    al = 0.0;
		    ar = 0.0;
		    for(j = 0; j < 16; j++)
		      {
			al += vl[j] * left[k * 256 + l * 16 + j];
			ar += vr[j] * right[k * 256 + l * 16 + j];
		      }
		    
		    x1px2 = al * ar;
		    for(j = 0; j < 16; j++)
		      v[j] += x1px2 * extEV[16 * l + j];		
		  }
	      }
	   
	    ex3[i] = 0;		  		    	     	            
	  }         	
      }
      break;
    case TIP_INNER:
      {	    	       	
	for (i = 0; i < n; i++) 
	  {		   	   	  		    	    
	    for(k = 0; k < 4; k++)
	      {	   
		vl = &(tipVector[16 * tipX1[i]]);
		vr = &(x2[64 * i + 16 * k]);
		v =  &(x3[64 * i + 16 * k]);
		
		for(l = 0; l < 16; l++)	       
		  v[l] = 0;		 
		
		for(l = 0; l < 16; l++)
		  {		
		    al = 0.0;
		    ar = 0.0;
		    for(j = 0; j < 16; j++)
		      {
			al += vl[j] * left[k * 256 + l * 16 + j];
			ar += vr[j] * right[k * 256 + l * 16 + j];
		      }
		    
		    x1px2 = al * ar;
		    for(j = 0; j < 16; j++)
		      v[j] += x1px2 * extEV[16 * l + j];		
		  }
	      }
	    
	    ex3[i] = ex2[i];
	    v = &x3[64 * i];
	    scale = 1;
	    for(l = 0; scale && (l < 64); l++)
	      scale = (ABS(v[l]) <  minlikelihood);
	    
	    if (scale)	         
	      {	     	    	 
		for(l = 0; l < 64; l++)
		  v[l] *= twotothe256;		   
		ex3[i] += 1;  
	      }	 	           	 
	  }
      }
      break;
    case INNER_INNER:	      
      for (i = 0; i < n; i++) 
       {		 	
	 for(k = 0; k < 4; k++)
	   {	   
	     vl = &(x1[64 * i + 16 * k]);
	     vr = &(x2[64 * i + 16 * k]);
	     v =  &(x3[64 * i + 16 * k]);
	     	    
	     for(l = 0; l < 16; l++)	       
	       v[l] = 0;		 
	       
	     for(l = 0; l < 16; l++)
	       {		
		 al = 0.0;
		 ar = 0.0;
		 for(j = 0; j < 16; j++)
		   {
		     al += vl[j] * left[k * 256 + l * 16 + j];
		     ar += vr[j] * right[k * 256 + l * 16 + j];
		   }
		 
		 x1px2 = al * ar;
		 for(j = 0; j < 16; j++)
		   v[j] += x1px2 * extEV[16 * l + j];		
	       }
	   }
	 
	 ex3[i] = ex1[i] + ex2[i];
	 v = &(x3[64 * i]);
	 scale = 1;
	 for(l = 0; scale && (l < 64); l++)
	   scale = ((ABS(v[l]) <  minlikelihood));
	 
	 if (scale)	         
	   {	     	    	 
	     for(l = 0; l < 64; l++)
	       v[l] *= twotothe256;		   
	     ex3[i] += 1;  
	   }		 	 	 	   
       }
      break;
    default:
      assert(0);
    }  
}


static void newviewGTRGAMMASECONDARY_6(int tipCase,
				       double *x1, double *x2, double *x3, double *extEV, double *tipVector,
				       int    *ex1, int *ex2, int *ex3, unsigned char *tipX1, unsigned char *tipX2,
				       int n, double *left, double *right)
{  
  double  *v;
  double x1px2;    
  int  i, j, l, k, scale;   
  double *vl, *vr, al, ar;                
   
  switch(tipCase)
    {      
    case TIP_TIP:
      {      	
	for(i = 0; i < n; i++) 
	  {		     	    	    
	    for(k = 0; k < 4; k++)
	      {	   
		vl = &(tipVector[6 * tipX1[i]]);
		vr = &(tipVector[6 * tipX2[i]]);
		v =  &(x3[24 * i + 6 * k]);
		
		for(l = 0; l < 6; l++)	       
		  v[l] = 0;		 
		
		for(l = 0; l < 6; l++)
		  {		
		    al = 0.0;
		    ar = 0.0;
		    for(j = 0; j < 6; j++)
		      {
			al += vl[j] * left[k * 36 + l * 6 + j];
			ar += vr[j] * right[k * 36 + l * 6 + j];
		      }
		    
		    x1px2 = al * ar;
		    for(j = 0; j < 6; j++)
		      v[j] += x1px2 * extEV[6 * l + j];		
		  }
	      }
	   
	    ex3[i] = 0;		  		    	     	            
	  }         	
      }
      break;
    case TIP_INNER:
      {	    	       	
	for (i = 0; i < n; i++) 
	  {		   	   	  		    	    
	    for(k = 0; k < 4; k++)
	      {	   
		vl = &(tipVector[6 * tipX1[i]]);
		vr = &(x2[24 * i + 6 * k]);
		v =  &(x3[24 * i + 6 * k]);
		
		for(l = 0; l < 6; l++)	       
		  v[l] = 0;		 
		
		for(l = 0; l < 6; l++)
		  {		
		    al = 0.0;
		    ar = 0.0;
		    for(j = 0; j < 6; j++)
		      {
			al += vl[j] * left[k * 36 + l * 6 + j];
			ar += vr[j] * right[k * 36 + l * 6 + j];
		      }
		    
		    x1px2 = al * ar;
		    for(j = 0; j < 6; j++)
		      v[j] += x1px2 * extEV[6 * l + j];		
		  }
	      }
	    
	    ex3[i] = ex2[i];
	    v = &x3[24 * i];
	    scale = 1;
	    for(l = 0; scale && (l < 24); l++)
	      scale = (ABS(v[l]) <  minlikelihood);
	    
	    if(scale)	         
	      {	     	    	 
		for(l = 0; l < 24; l++)
		  v[l] *= twotothe256;		   
		ex3[i] += 1;  
	      }	 	           	 
	  }
      }
      break;
    case INNER_INNER:	      
      for (i = 0; i < n; i++) 
       {		 	
	 for(k = 0; k < 4; k++)
	   {	   
	     vl = &(x1[24 * i + 6 * k]);
	     vr = &(x2[24 * i + 6 * k]);
	     v =  &(x3[24 * i + 6 * k]);
	     	    
	     for(l = 0; l < 6; l++)	       
	       v[l] = 0;		 
	       
	     for(l = 0; l < 6; l++)
	       {		
		 al = 0.0;
		 ar = 0.0;
		 for(j = 0; j < 6; j++)
		   {
		     al += vl[j] * left[k * 36 + l * 6 + j];
		     ar += vr[j] * right[k * 36 + l * 6 + j];
		   }
		 
		 x1px2 = al * ar;
		 for(j = 0; j < 6; j++)
		   v[j] += x1px2 * extEV[6 * l + j];		
	       }
	   }
	 
	 ex3[i] = ex1[i] + ex2[i];
	 v = &(x3[24 * i]);
	 scale = 1;
	 for(l = 0; scale && (l < 24); l++)
	   scale = ((ABS(v[l]) <  minlikelihood));
	 
	 if (scale)	         
	   {	     	    	 
	     for(l = 0; l < 24; l++)
	       v[l] *= twotothe256;		   
	     ex3[i] += 1;  
	   }		 	 	 	   
       }
      break;
    default:
      assert(0);
    }  
}


static void newviewGTRGAMMASECONDARY_7(int tipCase,
				       double *x1, double *x2, double *x3, double *extEV, double *tipVector,
				       int    *ex1, int *ex2, int *ex3, unsigned char *tipX1, unsigned char *tipX2,
				       int n, double *left, double *right)
{  
  double  *v;
  double x1px2;    
  int  i, j, l, k, scale;   
  double *vl, *vr, al, ar;                
   
  switch(tipCase)
    {      
    case TIP_TIP:
      {      	
	for(i = 0; i < n; i++) 
	  {		     	    	    
	    for(k = 0; k < 4; k++)
	      {	   
		vl = &(tipVector[7 * tipX1[i]]);
		vr = &(tipVector[7 * tipX2[i]]);
		v =  &(x3[28 * i + 7 * k]);
		
		for(l = 0; l < 7; l++)	       
		  v[l] = 0;		 
		
		for(l = 0; l < 7; l++)
		  {		
		    al = 0.0;
		    ar = 0.0;
		    for(j = 0; j < 7; j++)
		      {
			al += vl[j] * left[k * 49 + l * 7 + j];
			ar += vr[j] * right[k * 49 + l * 7 + j];
		      }
		    
		    x1px2 = al * ar;
		    for(j = 0; j < 7; j++)
		      v[j] += x1px2 * extEV[7 * l + j];		
		  }
	      }
	   
	    ex3[i] = 0;		  		    	     	            
	  }         	
      }
      break;
    case TIP_INNER:
      {	    	       	
	for (i = 0; i < n; i++) 
	  {		   	   	  		    	    
	    for(k = 0; k < 4; k++)
	      {	   
		vl = &(tipVector[7 * tipX1[i]]);
		vr = &(x2[28 * i + 7 * k]);
		v =  &(x3[28 * i + 7 * k]);
		
		for(l = 0; l < 7; l++)	       
		  v[l] = 0;		 
		
		for(l = 0; l < 7; l++)
		  {		
		    al = 0.0;
		    ar = 0.0;
		    for(j = 0; j < 7; j++)
		      {
			al += vl[j] * left[k * 49 + l * 7 + j];
			ar += vr[j] * right[k * 49 + l * 7 + j];
		      }
		    
		    x1px2 = al * ar;
		    for(j = 0; j < 7; j++)
		      v[j] += x1px2 * extEV[7 * l + j];		
		  }
	      }
	    
	    ex3[i] = ex2[i];
	    v = &x3[28 * i];
	    scale = 1;
	    for(l = 0; scale && (l < 28); l++)
	      scale = (ABS(v[l]) <  minlikelihood);
	    
	    if (scale)	         
	      {	     	    	 
		for(l = 0; l < 28; l++)
		  v[l] *= twotothe256;		   
		ex3[i] += 1;  
	      }	 	           	 
	  }
      }
      break;
    case INNER_INNER:	      
      for (i = 0; i < n; i++) 
       {		 	
	 for(k = 0; k < 4; k++)
	   {	   
	     vl = &(x1[28 * i + 7 * k]);
	     vr = &(x2[28 * i + 7 * k]);
	     v =  &(x3[28 * i + 7 * k]);
	     	    
	     for(l = 0; l < 7; l++)	       
	       v[l] = 0;		 
	       
	     for(l = 0; l < 7; l++)
	       {		
		 al = 0.0;
		 ar = 0.0;
		 for(j = 0; j < 7; j++)
		   {
		     al += vl[j] * left[k * 49 + l * 7 + j];
		     ar += vr[j] * right[k * 49 + l * 7 + j];
		   }
		 
		 x1px2 = al * ar;
		 for(j = 0; j < 7; j++)
		   v[j] += x1px2 * extEV[7 * l + j];		
	       }
	   }
	 
	 ex3[i] = ex1[i] + ex2[i];
	 v = &(x3[28 * i]);
	 scale = 1;
	 for(l = 0; scale && (l < 28); l++)
	   scale = ((ABS(v[l]) <  minlikelihood));
	 
	 if (scale)	         
	   {	     	    	 
	     for(l = 0; l < 28; l++)
	       v[l] *= twotothe256;		   
	     ex3[i] += 1;  
	   }		 	 	 	   
       }
      break;
    default:
      assert(0);
    }  
}







void computeTraversalInfo(nodeptr p, traversalInfo *ti, int *counter, int maxTips, int numBranches)
{
  if(isTip(p->number, maxTips))
    return;

  {     
    int i;
    nodeptr q = p->next->back;
    nodeptr r = p->next->next->back;
    
    if(isTip(r->number, maxTips) && isTip(q->number, maxTips))
      {	  
	while (! p->x)
	 {
	   if (! p->x)
	     getxnode(p); 	   
	 }

	ti[*counter].tipCase = TIP_TIP; 
	ti[*counter].pNumber = p->number;
	ti[*counter].qNumber = q->number;
	ti[*counter].rNumber = r->number;
	for(i = 0; i < numBranches; i++)
	  {
	    double z;
	    z = q->z[i];
	    z = (z > zmin) ? log(z) : log(zmin);
	    ti[*counter].qz[i] = z;

	    z = r->z[i];
	    z = (z > zmin) ? log(z) : log(zmin);
	    ti[*counter].rz[i] = z;	    
	  }     
	*counter = *counter + 1;
      }  
    else
      {
	if(isTip(r->number, maxTips) || isTip(q->number, maxTips))
	  {		
	    nodeptr tmp;

	    if(isTip(r->number, maxTips))
	      {
		tmp = r;
		r = q;
		q = tmp;
	      }

	    while ((! p->x) || (! r->x)) 
	      {	 	    
		if (! r->x) 
		  computeTraversalInfo(r, ti, counter, maxTips, numBranches);
		if (! p->x) 
		  getxnode(p);	
	      }
	    	   
	    ti[*counter].tipCase = TIP_INNER; 
	    ti[*counter].pNumber = p->number;
	    ti[*counter].qNumber = q->number;
	    ti[*counter].rNumber = r->number;
	    for(i = 0; i < numBranches; i++)
	      {
		double z;
		z = q->z[i];
		z = (z > zmin) ? log(z) : log(zmin);
		ti[*counter].qz[i] = z;
		
		z = r->z[i];
		z = (z > zmin) ? log(z) : log(zmin);
		ti[*counter].rz[i] = z;		
	      }   
	    
	    *counter = *counter + 1;
	  }
	else
	  {	 

	    while ((! p->x) || (! q->x) || (! r->x)) 
	      {
		if (! q->x) 
		  computeTraversalInfo(q, ti, counter, maxTips, numBranches);
		if (! r->x) 
		  computeTraversalInfo(r, ti, counter, maxTips, numBranches);
		if (! p->x) 
		  getxnode(p);	
	      }
   
	    ti[*counter].tipCase = INNER_INNER; 
	    ti[*counter].pNumber = p->number;
	    ti[*counter].qNumber = q->number;
	    ti[*counter].rNumber = r->number;
	    for(i = 0; i < numBranches; i++)
	      {
		double z;
		z = q->z[i];
		z = (z > zmin) ? log(z) : log(zmin);
		ti[*counter].qz[i] = z;
		
		z = r->z[i];
		z = (z > zmin) ? log(z) : log(zmin);
		ti[*counter].rz[i] = z;		
	      }   
	    
	    *counter = *counter + 1;
	  }
      }    
  }

}






void newviewIterative (tree *tr)
{	   
  traversalInfo *ti   = tr->td[0].ti;
  int i, model;   

  for(i = 1; i < tr->td[0].count; i++)
    {    
      traversalInfo *tInfo = &ti[i];                 
     
      for(model = 0; model < tr->NumberOfModels; model++)
	{    	
	  if(tr->executeModel[model])
	    { 
	      double 
		*x1_start = (double*)NULL,
		*x2_start = (double*)NULL, 
		*x3_start = (double*)NULL;
	      int 
		*ex1 = (int*)NULL, 
		*ex2 = (int*)NULL, 
		*ex3 = (int*)NULL;  
	      unsigned char 
		*tipX1 = (unsigned char *)NULL,
		*tipX2 = (unsigned char *)NULL;	     
	      double qz, rz;	      
	      int width =  tr->partitionData[model].width;
	      

	      switch(tInfo->tipCase)
		{
		case TIP_TIP:
		  tipX1    = tr->partitionData[model].yVector[tInfo->qNumber];
		  tipX2    = tr->partitionData[model].yVector[tInfo->rNumber];
		  
		  x3_start = tr->partitionData[model].xVector[tInfo->pNumber - tr->mxtips - 1];
		  ex3      = tr->partitionData[model].expVector[tInfo->pNumber - tr->mxtips - 1];
		  
		  break;
		case TIP_INNER:
		  tipX1    =  tr->partitionData[model].yVector[tInfo->qNumber];
		  
		  x2_start = tr->partitionData[model].xVector[tInfo->rNumber - tr->mxtips - 1];
		  ex2      = tr->partitionData[model].expVector[tInfo->rNumber - tr->mxtips - 1];
		  
		  x3_start = tr->partitionData[model].xVector[tInfo->pNumber - tr->mxtips - 1];
		  ex3      = tr->partitionData[model].expVector[tInfo->pNumber - tr->mxtips - 1];
		  
		  break;
		case INNER_INNER:	 
		  x1_start = tr->partitionData[model].xVector[tInfo->qNumber - tr->mxtips - 1];
		  ex1      = tr->partitionData[model].expVector[tInfo->qNumber - tr->mxtips - 1];	 
		  
		  x2_start = tr->partitionData[model].xVector[tInfo->rNumber - tr->mxtips - 1];
		  ex2      = tr->partitionData[model].expVector[tInfo->rNumber - tr->mxtips - 1];	  
		  
		  x3_start = tr->partitionData[model].xVector[tInfo->pNumber - tr->mxtips - 1];
		  ex3      = tr->partitionData[model].expVector[tInfo->pNumber - tr->mxtips - 1];
		  
		  break;
		default:
		  assert(0);
		}
	  
	      

	      if(tr->multiBranch)
		{
		  qz = tInfo->qz[model];
		  rz = tInfo->rz[model];
		}
	      else
		{
		  qz = tInfo->qz[0];
		  rz = tInfo->rz[0];
		}
	      

	      switch(tr->partitionData[model].dataType)
		{
		case BINARY_DATA:
		  switch(tr->rateHetModel)
		    {
		    case CAT:	    
		      {
			double *left  = (double *)malloc(4 * tr->NumberOfCategories * sizeof(double));                           
			double *right = (double *)malloc(4 * tr->NumberOfCategories * sizeof(double));
			
			makeP(qz, rz, tr->cdta->patrat,   tr->partitionData[model].EI, 
			      tr->partitionData[model].EIGN, tr->NumberOfCategories, 
			      left, right, BINARY_DATA);
			
			newviewGTRCAT_BINARY(tInfo->tipCase,  tr->partitionData[model].EV, tr->partitionData[model].rateCategory, 
					     x1_start, x2_start, x3_start, tr->partitionData[model].tipVector,
					     ex1, ex2, ex3, tipX1, tipX2,
					     width, left, right
					     );   
			
			free(left);
			free(right);	    
		      }
		      break;	  	   
		    case GAMMA:	 
		    case GAMMA_I:
		      {
			double *left  = (double *)malloc(16 * sizeof(double));                           
			double *right = (double *)malloc(16 * sizeof(double));	  
			
			makeP(qz, rz, tr->partitionData[model].gammaRates,  
			      tr->partitionData[model].EI, tr->partitionData[model].EIGN, 
			      4, left, right, BINARY_DATA);
			
			newviewGTRGAMMA_BINARY(tInfo->tipCase,
					       x1_start, x2_start, x3_start, tr->partitionData[model].EV, tr->partitionData[model].tipVector,
					       ex1, ex2, ex3, tipX1, tipX2,
					       width, left, right);
			
			free(left);
			free(right);
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
			double *left  = (double *)malloc(16 * tr->NumberOfCategories * sizeof(double));                           
			double *right = (double *)malloc(16 * tr->NumberOfCategories * sizeof(double));
			
			makeP(qz, rz, tr->cdta->patrat,   tr->partitionData[model].EI, 
			      tr->partitionData[model].EIGN, tr->NumberOfCategories, 
			      left, right, DNA_DATA);
			
			newviewGTRCAT(tInfo->tipCase,  tr->partitionData[model].EV, tr->partitionData[model].rateCategory, 
				      x1_start, x2_start, x3_start, tr->partitionData[model].tipVector,
				      ex1, ex2, ex3, tipX1, tipX2,
				      width, left, right
				      );   
			
			free(left);
			free(right);	    
		      }
		      break;	  	   
		    case GAMMA:	 
		    case GAMMA_I:
		      {
			double *left  = (double *)malloc(64 * sizeof(double));                           
			double *right = (double *)malloc(64 * sizeof(double));	  
			
			makeP(qz, rz, tr->partitionData[model].gammaRates,  
			      tr->partitionData[model].EI, tr->partitionData[model].EIGN, 
			      4, left, right, DNA_DATA);
			
			newviewGTRGAMMA(tInfo->tipCase,
					x1_start, x2_start, x3_start, tr->partitionData[model].EV, tr->partitionData[model].tipVector,
					ex1, ex2, ex3, tipX1, tipX2,
					width, left, right);
			
			free(left);
			free(right);
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
			double *left  = (double *)malloc(400 * tr-> NumberOfCategories * sizeof(double));
			double *right = (double *)malloc(400 * tr-> NumberOfCategories * sizeof(double));	    		
			
			makeP(qz, rz, tr->cdta->patrat, 
			      tr->partitionData[model].EI, 
			      tr->partitionData[model].EIGN, 
			      tr->NumberOfCategories, left, right, AA_DATA);
			
			newviewGTRCATPROT(tInfo->tipCase,  tr->partitionData[model].EV, tr->partitionData[model].rateCategory, 
					  x1_start, x2_start, x3_start, tr->partitionData[model].tipVector,
					  ex1, ex2, ex3, tipX1, tipX2, width, left, right);	      
			
			free(left);
			free(right);
		      }	     	      
		      break;	      
		    case GAMMA:
		    case GAMMA_I:		  	    
		      {
			double *left  = (double *)malloc(1600 * sizeof(double));
			double *right = (double *)malloc(1600 * sizeof(double));	   
			
			makeP(qz, rz, tr->partitionData[model].gammaRates,  
			      tr->partitionData[model].EI, 
			      tr->partitionData[model].EIGN, 
			      4, left, right, AA_DATA);
			
			newviewGTRGAMMAPROT(tInfo->tipCase,
					    x1_start, x2_start, x3_start,
					    tr->partitionData[model].EV, 
					    tr->partitionData[model].tipVector,
					    ex1, ex2, ex3, tipX1, tipX2,
					    width, left, right);     	            	      
			
			free(left);
			free(right);
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
			double *left  = (double *)malloc(36 * tr-> NumberOfCategories * sizeof(double));
			double *right = (double *)malloc(36 * tr-> NumberOfCategories * sizeof(double));	    		
			
			makeP(qz, rz, tr->cdta->patrat, 
			      tr->partitionData[model].EI, 
			      tr->partitionData[model].EIGN, 
			      tr->NumberOfCategories, left, right, SECONDARY_DATA_6);
			
			newviewGTRCATSECONDARY_6(tInfo->tipCase,  tr->partitionData[model].EV, tr->partitionData[model].rateCategory, 
						 x1_start, x2_start, x3_start, tr->partitionData[model].tipVector,
						 ex1, ex2, ex3, tipX1, tipX2, width, left, right);	      
			
			free(left);
			free(right);
		      }	     	      
		      break;	      
		    case GAMMA:
		    case GAMMA_I:		  	    
		      {
			double *left  = (double *)malloc(144 * sizeof(double));
			double *right = (double *)malloc(144 * sizeof(double));	   
			
			makeP(qz, rz, tr->partitionData[model].gammaRates,  
			      tr->partitionData[model].EI, 
			      tr->partitionData[model].EIGN, 
			      4, left, right, SECONDARY_DATA_6);
			
			newviewGTRGAMMASECONDARY_6(tInfo->tipCase,
						 x1_start, x2_start, x3_start,
						 tr->partitionData[model].EV, 
						 tr->partitionData[model].tipVector,
						 ex1, ex2, ex3, tipX1, tipX2,
						 width, left, right);     	            	      
			
			free(left);
			free(right);
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
			double *left  = (double *)malloc(49 * tr-> NumberOfCategories * sizeof(double));
			double *right = (double *)malloc(49 * tr-> NumberOfCategories * sizeof(double));	    		
			
			makeP(qz, rz, tr->cdta->patrat, 
			      tr->partitionData[model].EI, 
			      tr->partitionData[model].EIGN, 
			      tr->NumberOfCategories, left, right, SECONDARY_DATA_7);
			
			newviewGTRCATSECONDARY_7(tInfo->tipCase,  tr->partitionData[model].EV, tr->partitionData[model].rateCategory, 
					       x1_start, x2_start, x3_start, tr->partitionData[model].tipVector,
					       ex1, ex2, ex3, tipX1, tipX2, width, left, right);	      
			
			free(left);
			free(right);
		      }	     	      
		      break;	      
		    case GAMMA:
		    case GAMMA_I:		  	    
		      {
			double *left  = (double *)malloc(196 * sizeof(double));
			double *right = (double *)malloc(196 * sizeof(double));	   
			
			makeP(qz, rz, tr->partitionData[model].gammaRates,  
			      tr->partitionData[model].EI, 
			      tr->partitionData[model].EIGN, 
			      4, left, right, SECONDARY_DATA_7);
			
			newviewGTRGAMMASECONDARY_7(tInfo->tipCase,
						   x1_start, x2_start, x3_start,
						   tr->partitionData[model].EV, 
						   tr->partitionData[model].tipVector,
						   ex1, ex2, ex3, tipX1, tipX2,
						   width, left, right);     	            	      
			
			free(left);
			free(right);
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
			double *left  = (double *)malloc(256 * tr-> NumberOfCategories * sizeof(double));
			double *right = (double *)malloc(256 * tr-> NumberOfCategories * sizeof(double));	    		
			
			makeP(qz, rz, tr->cdta->patrat, 
			      tr->partitionData[model].EI, 
			      tr->partitionData[model].EIGN, 
			      tr->NumberOfCategories, left, right, SECONDARY_DATA);
			
			newviewGTRCATSECONDARY(tInfo->tipCase,  tr->partitionData[model].EV, tr->partitionData[model].rateCategory, 
					       x1_start, x2_start, x3_start, tr->partitionData[model].tipVector,
					       ex1, ex2, ex3, tipX1, tipX2, width, left, right);	      
			
			free(left);
			free(right);
		      }	     	      
		      break;	      
		    case GAMMA:
		    case GAMMA_I:		  	    
		      {
			double *left  = (double *)malloc(1024 * sizeof(double));
			double *right = (double *)malloc(1024 * sizeof(double));	   
			
			makeP(qz, rz, tr->partitionData[model].gammaRates,  
			      tr->partitionData[model].EI, 
			      tr->partitionData[model].EIGN, 
			      4, left, right, SECONDARY_DATA);
			
			newviewGTRGAMMASECONDARY(tInfo->tipCase,
						 x1_start, x2_start, x3_start,
						 tr->partitionData[model].EV, 
						 tr->partitionData[model].tipVector,
						 ex1, ex2, ex3, tipX1, tipX2,
						 width, left, right);     	            	      
			
			free(left);
			free(right);
		      }	  
		      break;
		    default:
		      assert(0);
		    }
		  break;
		default:
		  assert(0);
		}
	    }
	}
    }
}

void newviewGeneric (tree *tr, nodeptr p)
{
  if(isTip(p->number, tr->mxtips)) 
    return;


  { 	           
    tr->td[0].count = 1;
    computeTraversalInfo(p, &(tr->td[0].ti[0]), &(tr->td[0].count), tr->mxtips, tr->numBranches);   

    if(tr->td[0].count > 1)
      {
#ifdef _USE_PTHREADS  
	masterBarrier(THREAD_NEWVIEW, tr);         
#else
	newviewIterative(tr);   
#endif
      }   
  }
}


void newviewGenericMasked(tree *tr, nodeptr p)
{
  if(isTip(p->number, tr->mxtips)) 
    return;
  
  {
    int i;

    for(i = 0; i < tr->NumberOfModels; i++)
      {      
	if(tr->partitionConverged[i])
	  tr->executeModel[i] = FALSE;
	else
	  tr->executeModel[i] = TRUE;	
      }
 	           
    tr->td[0].count = 1;
    computeTraversalInfo(p, &(tr->td[0].ti[0]), &(tr->td[0].count), tr->mxtips, tr->numBranches);   

    if(tr->td[0].count > 1)
      {
#ifdef _USE_PTHREADS  
	masterBarrier(THREAD_NEWVIEW_MASKED, tr);         
#else
	newviewIterative(tr);   
#endif
      }   

    for(i = 0; i < tr->NumberOfModels; i++)          
      tr->executeModel[i] = TRUE;    
  }
}



