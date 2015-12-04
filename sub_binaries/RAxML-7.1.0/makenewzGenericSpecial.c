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
 *  Alexandros Stamatakis:"RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses with 
 *  thousands of taxa and mixed models". 
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
extern volatile double *reductionBufferTwo;
extern volatile int NumberOfThreads;
#endif




/*******************/

static void sumCAT_BINARY(int tipCase, double *sum, double *x1_start, double *x2_start, double *tipVector,
			  unsigned char *tipX1, unsigned char *tipX2, int n)
{
  int i, j;
  double *x1, *x2;  
  
  switch(tipCase)
    {
    case TIP_TIP:
      for (i = 0; i < n; i++) 
	{    
	  x1 = &(tipVector[2 * tipX1[i]]);
	  x2 = &(tipVector[2 * tipX2[i]]);
	  
	  for(j = 0; j < 2; j++)
	    sum[i * 2 + j]     = x1[j] * x2[j];	  		    	    	  
	}    
      break;
    case TIP_INNER:         
      for (i = 0; i < n; i++) 
	{    
	  x1 = &(tipVector[2 * tipX1[i]]);
	  x2 = &x2_start[2 * i];
	  
	  for(j = 0; j < 2; j++)
	    sum[i * 2 + j]     = x1[j] * x2[j];	 				 
	}   
      break;
    case INNER_INNER:      
      for (i = 0; i < n; i++) 
	{     
	  x1 = &x1_start[2 * i];
	  x2 = &x2_start[2 * i];
	  
	  for(j = 0; j < 2; j++)
	    sum[i * 2 + j]     = x1[j] * x2[j];	 				   
	}
      break;  
    default:
      assert(0);
    }
}

static void sumCAT(int tipCase, double *sum, double *x1_start, double *x2_start, double *tipVector,
		   unsigned char *tipX1, unsigned char *tipX2, int n)
{
  int i, j;
  double *x1, *x2;  
  
  switch(tipCase)
    {
    case TIP_TIP:
      for (i = 0; i < n; i++) 
	{    
	  x1 = &(tipVector[4 * tipX1[i]]);
	  x2 = &(tipVector[4 * tipX2[i]]);
	  
	  for(j = 0; j < 4; j++)
	    sum[i * 4 + j]     = x1[j] * x2[j];	  		    	    	  
	}    
      break;
    case TIP_INNER:         
      for (i = 0; i < n; i++) 
	{    
	  x1 = &(tipVector[4 * tipX1[i]]);
	  x2 = &x2_start[4 * i];
	  
	  for(j = 0; j < 4; j++)
	    sum[i * 4 + j]     = x1[j] * x2[j];	 				 
	}   
      break;
    case INNER_INNER:      
      for (i = 0; i < n; i++) 
	{     
	  x1 = &x1_start[4 * i];
	  x2 = &x2_start[4 * i];
	  
	  for(j = 0; j < 4; j++)
	    sum[i * 4 + j]     = x1[j] * x2[j];	 				   
	}
      break;  
    default:
      assert(0);
    }
}

static void coreGTRCAT_BINARY(int upper, int numberOfCategories, double *sum, 
			      volatile double *d1, volatile double *d2, double *wrptr, double *wr2ptr, 
			      double *rptr, double *EIGN, int *cptr, double lz)
{
  int i;
  double 
    *d, *d_start,
    tmp_0, inv_Li, dlnLidlz, d2lnLidlz2,
    dlnLdlz = 0.0,
    d2lnLdlz2 = 0.0;
  double e[2];
  double dd1;
 
  e[0] = EIGN[0];
  e[1] = EIGN[0] * EIGN[0];
  
          
  d = d_start = (double *)malloc(numberOfCategories * sizeof(double));   

  dd1 = e[0] * lz;	  
   
  for(i = 0; i < numberOfCategories; i++)            
    d[i] = exp(dd1 * rptr[i]);     

  for (i = 0; i < upper; i++) 	
    {	    	 	   	    
      d = &d_start[cptr[i]];
		
      inv_Li = sum[2 * i];	   
      inv_Li += (tmp_0 = d[0] * sum[2 * i + 1]);     
		
      inv_Li = 1.0/inv_Li;
		
      dlnLidlz   = tmp_0 * e[0];
      d2lnLidlz2 = tmp_0 * e[1];		     
      
      dlnLidlz   *= inv_Li;
      d2lnLidlz2 *= inv_Li;
            
      dlnLdlz   += wrptr[i] * dlnLidlz;
      d2lnLdlz2 += wr2ptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    }

  *d1 = dlnLdlz;
  *d2 = d2lnLdlz2;

  free(d_start);
}
static void coreGTRCAT(int upper, int numberOfCategories, double *sum, 
		       volatile double *d1, volatile double *d2, double *wrptr, double *wr2ptr, 
		       double *rptr, double *EIGN, int *cptr, double lz)
{
  int i;
  double 
    *d, *d_start,
    tmp_0, tmp_1, tmp_2, inv_Li, dlnLidlz, d2lnLidlz2,
    dlnLdlz = 0.0,
    d2lnLdlz2 = 0.0;
  double e[6];
  double dd1, dd2, dd3;
 
  e[0] = EIGN[0];
  e[1] = EIGN[0] * EIGN[0];
  e[2] = EIGN[1];
  e[3] = EIGN[1] * EIGN[1];
  e[4] = EIGN[2];
  e[5] = EIGN[2] * EIGN[2];
          
  d = d_start = (double *)malloc(numberOfCategories * 4 * sizeof(double));   

  dd1 = e[0] * lz;
  dd2 = e[2] * lz;
  dd3 = e[4] * lz;	  
   
  for(i = 0; i < numberOfCategories; i++)
    {        
      d[i * 4] = exp(dd1 * rptr[i]); 
      d[i * 4 + 1] = exp(dd2 * rptr[i]);
      d[i * 4 + 2] = exp(dd3 * rptr[i]);	    	    	  	   
    }       

  for (i = 0; i < upper; i++) 	
    {	    	 	   	    
      d = &d_start[4 * cptr[i]];
		
      inv_Li = sum[4 * i];	   
      inv_Li += (tmp_0 = d[0] * sum[4 * i + 1]);
      inv_Li += (tmp_1 = d[1] * sum[4 * i + 2]);
      inv_Li += (tmp_2 = d[2] * sum[4 * i + 3]);
		
      inv_Li = 1.0/inv_Li;
		
      dlnLidlz   = tmp_0 * e[0];
      d2lnLidlz2 = tmp_0 * e[1];
		
      dlnLidlz   += tmp_1 * e[2];
      d2lnLidlz2 += tmp_1 * e[3];
      
      dlnLidlz   += tmp_2 * e[4];	    	   	   	   
      d2lnLidlz2 += tmp_2 * e[5];
      
      dlnLidlz   *= inv_Li;
      d2lnLidlz2 *= inv_Li;
      
      
      dlnLdlz   += wrptr[i] * dlnLidlz;
      d2lnLdlz2 += wr2ptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    }

  *d1 = dlnLdlz;
  *d2 = d2lnLdlz2;

  free(d_start);
}









static void sumGTRCATPROT(int tipCase, double *sumtable, double *x1, double *x2, double *tipVector,
			  unsigned char *tipX1, unsigned char *tipX2, int n)
{
  int i, l;
  double *sum, *left, *right;

  switch(tipCase)
    {
    case TIP_TIP:          
      for (i = 0; i < n; i++) 
	{    
	  left  = &(tipVector[20 * tipX1[i]]);
	  right = &(tipVector[20 * tipX2[i]]);	
	  sum = &sumtable[20 * i];
	  
	  for(l = 0; l < 20; l++)
	    sum[l] = left[l] * right[l];				 
	}    
      break;
    case TIP_INNER:     
      for (i = 0; i < n; i++) 
	{    
	  left = &(tipVector[20 * tipX1[i]]);			       
	  right = &x2[20 * i];
	  sum = &sumtable[20 * i];

	  for(l = 0; l < 20; l++)
	    sum[l] = left[l] * right[l];	  				      
	}   
      break;
    case INNER_INNER:	     
      for (i = 0; i < n; i++) 
	{     
	  left  = &x1[20 * i];
	  right = &x2[20 * i];
	  sum = &sumtable[20 * i];
	  
	  for(l = 0; l < 20; l++)
	    sum[l] = left[l] * right[l];	  			
	}   
      break;
    default:
      assert(0);
    }
}

static void sumGTRCATSECONDARY(int tipCase, double *sumtable, double *x1, double *x2, double *tipVector,
			       unsigned char *tipX1, unsigned char *tipX2, int n)
{
  int i, l;
  double *sum, *left, *right;

  switch(tipCase)
    {
    case TIP_TIP:          
      for (i = 0; i < n; i++) 
	{    
	  left  = &(tipVector[16 * tipX1[i]]);
	  right = &(tipVector[16 * tipX2[i]]);	
	  sum = &sumtable[16 * i];
	  
	  for(l = 0; l < 16; l++)
	    sum[l] = left[l] * right[l];				 
	}    
      break;
    case TIP_INNER:     
      for (i = 0; i < n; i++) 
	{    
	  left = &(tipVector[16 * tipX1[i]]);			       
	  right = &x2[16 * i];
	  sum = &sumtable[16 * i];

	  for(l = 0; l < 16; l++)
	    sum[l] = left[l] * right[l];	  				      
	}   
      break;
    case INNER_INNER:	     
      for (i = 0; i < n; i++) 
	{     
	  left  = &x1[16 * i];
	  right = &x2[16 * i];
	  sum = &sumtable[16 * i];
	  
	  for(l = 0; l < 16; l++)
	    sum[l] = left[l] * right[l];	  			
	}   
      break;
    default:
      assert(0);
    }
}

static void sumGTRCATSECONDARY_6(int tipCase, double *sumtable, double *x1, double *x2, double *tipVector,
				 unsigned char *tipX1, unsigned char *tipX2, int n)
{
  int i, l;
  double *sum, *left, *right;

  switch(tipCase)
    {
    case TIP_TIP:          
      for (i = 0; i < n; i++) 
	{    
	  left  = &(tipVector[6 * tipX1[i]]);
	  right = &(tipVector[6 * tipX2[i]]);	
	  sum = &sumtable[6 * i];
	  
	  for(l = 0; l < 6; l++)
	    sum[l] = left[l] * right[l];				 
	}    
      break;
    case TIP_INNER:     
      for (i = 0; i < n; i++) 
	{    
	  left = &(tipVector[6 * tipX1[i]]);			       
	  right = &x2[6 * i];
	  sum = &sumtable[6 * i];

	  for(l = 0; l < 6; l++)
	    sum[l] = left[l] * right[l];	  				      
	}   
      break;
    case INNER_INNER:	     
      for (i = 0; i < n; i++) 
	{     
	  left  = &x1[6 * i];
	  right = &x2[6 * i];
	  sum = &sumtable[6 * i];
	  
	  for(l = 0; l < 6; l++)
	    sum[l] = left[l] * right[l];	  			
	}   
      break;
    default:
      assert(0);
    }
}

static void sumGTRCATSECONDARY_7(int tipCase, double *sumtable, double *x1, double *x2, double *tipVector,
				 unsigned char *tipX1, unsigned char *tipX2, int n)
{
  int i, l;
  double *sum, *left, *right;

  switch(tipCase)
    {
    case TIP_TIP:          
      for (i = 0; i < n; i++) 
	{    
	  left  = &(tipVector[7 * tipX1[i]]);
	  right = &(tipVector[7 * tipX2[i]]);	
	  sum = &sumtable[7 * i];
	  
	  for(l = 0; l < 7; l++)
	    sum[l] = left[l] * right[l];				 
	}    
      break;
    case TIP_INNER:     
      for (i = 0; i < n; i++) 
	{    
	  left = &(tipVector[7 * tipX1[i]]);			       
	  right = &x2[7 * i];
	  sum = &sumtable[7 * i];

	  for(l = 0; l < 7; l++)
	    sum[l] = left[l] * right[l];	  				      
	}   
      break;
    case INNER_INNER:	     
      for (i = 0; i < n; i++) 
	{     
	  left  = &x1[7 * i];
	  right = &x2[7 * i];
	  sum = &sumtable[7 * i];
	  
	  for(l = 0; l < 7; l++)
	    sum[l] = left[l] * right[l];	  			
	}   
      break;
    default:
      assert(0);
    }
}

static void coreGTRCATPROT(double *EIGN, double lz, int numberOfCategories, double *rptr, int *cptr, int upper,
			   double *wrptr, double *wr2ptr,  volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double *sumtable)
{
  int i, l;
  double *d1, *d_start, *sum; 
  double e[20], s[20], dd[20]; 
  double tmp;
  double inv_Li, dlnLidlz, d2lnLidlz2; 
  double  dlnLdlz = 0.0;                 
  double  d2lnLdlz2 = 0.0; 
        
  d1 = d_start = (double *)malloc(numberOfCategories * 20 * sizeof(double));   
	
  for(l = 1; l < 20; l++)
    {
      e[l]  = EIGN[l-1] * EIGN[l-1];
      s[l]  = EIGN[l-1];
      dd[l] = s[l] * lz; 
    }    
    
  for(i = 0; i < numberOfCategories; i++)
    {     
      for(l = 1; l < 20; l++)
	d1[20 * i + l] = exp(dd[l] * rptr[i]); 	      	    	  	   
    }       

  for (i = 0; i < upper; i++) 	
    {	   	 	    
      d1 = &d_start[20 * cptr[i]];
      sum = &sumtable[20 * i];

      inv_Li     = sum[0];
      dlnLidlz   = 0.0;
      d2lnLidlz2 = 0.0;
      	         	      
      for(l = 1; l < 20; l++)	      
	{
	  inv_Li     += (tmp = d1[l] * sum[l]);
	  dlnLidlz   += tmp *  s[l];
	  d2lnLidlz2 += tmp *  e[l];
	}	   	    
	      
      inv_Li = 1.0/inv_Li;  	    	    	    	    	    
	      
      dlnLidlz   *= inv_Li;
      d2lnLidlz2 *= inv_Li;

      dlnLdlz  += wrptr[i] * dlnLidlz;
      d2lnLdlz2 += wr2ptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    }

  *ext_dlnLdlz   = dlnLdlz;
  *ext_d2lnLdlz2 = d2lnLdlz2;

  free(d_start);
}

static void coreGTRCATSECONDARY(double *EIGN, double lz, int numberOfCategories, double *rptr, int *cptr, int upper,
				double *wrptr, double *wr2ptr,  volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double *sumtable)
{
  int i, l;
  double *d1, *d_start, *sum; 
  double e[16], s[16], dd[16]; 
  double tmp;
  double inv_Li, dlnLidlz, d2lnLidlz2; 
  double  dlnLdlz = 0.0;                 
  double  d2lnLdlz2 = 0.0; 
        
  d1 = d_start = (double *)malloc(numberOfCategories * 16 * sizeof(double));   
	
  for(l = 1; l < 16; l++)
    {
      e[l]  = EIGN[l-1] * EIGN[l-1];
      s[l]  = EIGN[l-1];
      dd[l] = s[l] * lz; 
    }    
    
  for(i = 0; i < numberOfCategories; i++)
    {     
      for(l = 1; l < 16; l++)
	d1[16 * i + l] = exp(dd[l] * rptr[i]); 	      	    	  	   
    }       

  for (i = 0; i < upper; i++) 	
    {	   	 	    
      d1 = &d_start[16 * cptr[i]];
      sum = &sumtable[16 * i];

      inv_Li     = sum[0];
      dlnLidlz   = 0.0;
      d2lnLidlz2 = 0.0;
      	         	      
      for(l = 1; l < 16; l++)	      
	{
	  inv_Li     += (tmp = d1[l] * sum[l]);
	  dlnLidlz   += tmp *  s[l];
	  d2lnLidlz2 += tmp *  e[l];
	}	   	    
	      
      inv_Li = 1.0/inv_Li;  	    	    	    	    	    
	      
      dlnLidlz   *= inv_Li;
      d2lnLidlz2 *= inv_Li;

      dlnLdlz  += wrptr[i] * dlnLidlz;
      d2lnLdlz2 += wr2ptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    }

  *ext_dlnLdlz   = dlnLdlz;
  *ext_d2lnLdlz2 = d2lnLdlz2;

  free(d_start);
}

static void coreGTRCATSECONDARY_6(double *EIGN, double lz, int numberOfCategories, double *rptr, int *cptr, int upper,
				  double *wrptr, double *wr2ptr,  volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double *sumtable)
{
  int i, l;
  double *d1, *d_start, *sum; 
  double e[6], s[6], dd[6]; 
  double tmp;
  double inv_Li, dlnLidlz, d2lnLidlz2; 
  double  dlnLdlz = 0.0;                 
  double  d2lnLdlz2 = 0.0; 
        
  d1 = d_start = (double *)malloc(numberOfCategories * 6 * sizeof(double));   
	
  for(l = 1; l < 6; l++)
    {
      e[l]  = EIGN[l-1] * EIGN[l-1];
      s[l]  = EIGN[l-1];
      dd[l] = s[l] * lz; 
    }    
    
  for(i = 0; i < numberOfCategories; i++)
    {     
      for(l = 1; l < 6; l++)
	d1[6 * i + l] = exp(dd[l] * rptr[i]); 	      	    	  	   
    }       

  for (i = 0; i < upper; i++) 	
    {	   	 	    
      d1 = &d_start[6 * cptr[i]];
      sum = &sumtable[6 * i];

      inv_Li     = sum[0];
      dlnLidlz   = 0.0;
      d2lnLidlz2 = 0.0;
      	         	      
      for(l = 1; l < 6; l++)	      
	{
	  inv_Li     += (tmp = d1[l] * sum[l]);
	  dlnLidlz   += tmp *  s[l];
	  d2lnLidlz2 += tmp *  e[l];
	}	   	    
	      
      inv_Li = 1.0/inv_Li;  	    	    	    	    	    
	      
      dlnLidlz   *= inv_Li;
      d2lnLidlz2 *= inv_Li;

      dlnLdlz  += wrptr[i] * dlnLidlz;
      d2lnLdlz2 += wr2ptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    }

  *ext_dlnLdlz   = dlnLdlz;
  *ext_d2lnLdlz2 = d2lnLdlz2;

  free(d_start);
}

static void coreGTRCATSECONDARY_7(double *EIGN, double lz, int numberOfCategories, double *rptr, int *cptr, int upper,
				  double *wrptr, double *wr2ptr,  volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double *sumtable)
{
  int i, l;
  double *d1, *d_start, *sum; 
  double e[7], s[7], dd[7]; 
  double tmp;
  double inv_Li, dlnLidlz, d2lnLidlz2; 
  double  dlnLdlz = 0.0;                 
  double  d2lnLdlz2 = 0.0; 
        
  d1 = d_start = (double *)malloc(numberOfCategories * 7 * sizeof(double));   
	
  for(l = 1; l < 7; l++)
    {
      e[l]  = EIGN[l-1] * EIGN[l-1];
      s[l]  = EIGN[l-1];
      dd[l] = s[l] * lz; 
    }    
    
  for(i = 0; i < numberOfCategories; i++)
    {     
      for(l = 1; l < 7; l++)
	d1[7 * i + l] = exp(dd[l] * rptr[i]); 	      	    	  	   
    }       

  for (i = 0; i < upper; i++) 	
    {	   	 	    
      d1 = &d_start[7 * cptr[i]];
      sum = &sumtable[7 * i];

      inv_Li     = sum[0];
      dlnLidlz   = 0.0;
      d2lnLidlz2 = 0.0;
      	         	      
      for(l = 1; l < 7; l++)	      
	{
	  inv_Li     += (tmp = d1[l] * sum[l]);
	  dlnLidlz   += tmp *  s[l];
	  d2lnLidlz2 += tmp *  e[l];
	}	   	    
	      
      inv_Li = 1.0/inv_Li;  	    	    	    	    	    
	      
      dlnLidlz   *= inv_Li;
      d2lnLidlz2 *= inv_Li;

      dlnLdlz  += wrptr[i] * dlnLidlz;
      d2lnLdlz2 += wr2ptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    }

  *ext_dlnLdlz   = dlnLdlz;
  *ext_d2lnLdlz2 = d2lnLdlz2;

  free(d_start);
}



static void coreGTRGAMMAINVAR_BINARY(double propInvar, double *frequencies, double *gammaRates, double *EIGN,
				     double *sumtable,  volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2,
				     int *iptr, int *wrptr, int upper, double lz)
{   
  double  *sum, diagptable[12];
  int     i, j;
    double  dlnLdlz = 0;
  double d2lnLdlz2 = 0;
  double ki, kisqr;    
  double freqs[4]; 
  double scaler =  0.25 * (1.0 - propInvar);   
  double tmp_1;
  double inv_Li, dlnLidlz, d2lnLidlz2;
  
  freqs[0] = frequencies[0] * propInvar;
  freqs[1] = frequencies[1] * propInvar; 
  
  for(i = 0; i < 4; i++)
    {
      ki = gammaRates[i];
      kisqr = ki * ki;
	      
      diagptable[i * 3]     = exp (EIGN[0] * ki * lz);           
      diagptable[i * 3 + 1] = EIGN[0] * ki;
      diagptable[i * 3 + 2] = EIGN[0] * EIGN[0] * kisqr;           
    }

  for (i = 0; i < upper; i++) 
    {	    	    	        
      sum = &(sumtable[i * 8]);

      inv_Li      = 0.0;
      dlnLidlz    = 0.0;
      d2lnLidlz2  = 0.0;
       
      for(j = 0; j < 4; j++)
	{
	  inv_Li += sum[2 * j];
	   
	  
	  tmp_1      =  diagptable[3 * j] * sum[2 * j + 1];
	  inv_Li     += tmp_1;
	  dlnLidlz   += tmp_1 * diagptable[3 * j + 1];
	  d2lnLidlz2 += tmp_1 * diagptable[3 * j + 2];	    	       	    	    
	 }	    	
      
      inv_Li *= scaler;	   
      
      if(iptr[i] < 2)	      
	inv_Li += freqs[iptr[i]];
      
      inv_Li = 1.0 / inv_Li;
      
      dlnLidlz   *= inv_Li;
      d2lnLidlz2 *= inv_Li;
      
      dlnLidlz *= scaler;
      d2lnLidlz2 *= scaler;
      
      dlnLdlz  += wrptr[i] * dlnLidlz;
      d2lnLdlz2 += wrptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    }

  *ext_dlnLdlz   = dlnLdlz;
  *ext_d2lnLdlz2 = d2lnLdlz2;	  
}



static void coreGTRGAMMAINVAR(double propInvar, double *frequencies, double *gammaRates, double *EIGN,
			      double *sumtable,  volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2,
			      int *iptr, int *wrptr, int upper, double lz)
{   
  double  *sum, diagptable[64];
  int     i, j, k;
    double  dlnLdlz = 0;
  double d2lnLdlz2 = 0;
  double ki, kisqr;    
  double freqs[4]; 
  double scaler =  0.25 * (1.0 - propInvar);   
  double tmp_1;
  double inv_Li, dlnLidlz, d2lnLidlz2;
  
  freqs[0] = frequencies[0] * propInvar;
  freqs[1] = frequencies[1] * propInvar;
  freqs[2] = frequencies[2] * propInvar;
  freqs[3] = frequencies[3] * propInvar;           
  
  for(i = 0; i < 4; i++)
    {
      ki = gammaRates[i];
      kisqr = ki * ki;
	      
      diagptable[i * 16]     = exp (EIGN[0] * ki * lz);
      diagptable[i * 16 + 1] = exp (EIGN[1] * ki * lz);
      diagptable[i * 16 + 2] = exp (EIGN[2] * ki * lz);
      
      diagptable[i * 16 + 3] = EIGN[0] * ki;
      diagptable[i * 16 + 4] = EIGN[0] * EIGN[0] * kisqr;
      
      diagptable[i * 16 + 5] = EIGN[1] * ki;
      diagptable[i * 16 + 6] = EIGN[1] * EIGN[1] * kisqr;
      
      diagptable[i * 16 + 7] = EIGN[2] * ki;
      diagptable[i * 16 + 8] = EIGN[2] * EIGN[2] * kisqr;
    }

  for (i = 0; i < upper; i++) 
    {	    	    	        
      sum = &(sumtable[i * 16]);

       inv_Li      = 0.0;
       dlnLidlz    = 0.0;
       d2lnLidlz2  = 0.0;
       
       for(j = 0; j < 4; j++)
	 {
	   inv_Li += sum[4 * j];
	   
	   for(k = 0; k < 3; k++)
	     {
	       tmp_1      =  diagptable[16 * j + k] * sum[4 * j + k + 1];
	       inv_Li     += tmp_1;
	       dlnLidlz   += tmp_1 * diagptable[16 * j + k * 2 + 3];
	       d2lnLidlz2 += tmp_1 * diagptable[16 * j + k * 2 + 4];
	     }	       	       	    	    
	 }	    	
      
      inv_Li *= scaler;	   
      
      if(iptr[i] < 4)	      
	inv_Li += freqs[iptr[i]];
      
      inv_Li = 1.0 / inv_Li;
      
      dlnLidlz   *= inv_Li;
      d2lnLidlz2 *= inv_Li;
      
      dlnLidlz *= scaler;
      d2lnLidlz2 *= scaler;
      
      dlnLdlz  += wrptr[i] * dlnLidlz;
      d2lnLdlz2 += wrptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    }

  *ext_dlnLdlz   = dlnLdlz;
  *ext_d2lnLdlz2 = d2lnLdlz2;	  
}


/* C-OPT the two functions below are used to optimize the branch lengths of the tree.
   together, they account for approx 25-30% of overall run time sumGAMMA requires less 
   overall time than coreGAMMA though. */

static void sumGAMMA_BINARY(int tipCase, double *sumtable, double *x1_start, double *x2_start, double *tipVector,
			    unsigned char *tipX1, unsigned char *tipX2, int n)
{
  double *x1, *x2, *sum;
  int i, j, k;

  /* C-OPT once again switch over possible configurations at inner node */

  switch(tipCase)
    {
    case TIP_TIP:   
      /* C-OPT main for loop overt alignment length */
      for (i = 0; i < n; i++) 
	{     
	  x1 = &(tipVector[2 * tipX1[i]]);
	  x2 = &(tipVector[2 * tipX2[i]]);
	  sum = &sumtable[i * 8];
	  
	  for(j = 0; j < 4; j++)
	    for(k = 0; k < 2; k++)
	      sum[j * 2 + k] = x1[k] * x2[k];	           
	}    
      break;
    case TIP_INNER: 
      for (i = 0; i < n; i++) 
	{     
	  x1  = &(tipVector[2 * tipX1[i]]);
	  x2  = &x2_start[8 * i];
	  sum = &sumtable[8 * i];
	  
	  for(j = 0; j < 4; j++)
	    for(k = 0; k < 2; k++)
	      sum[j * 2 + k] = x1[k] * x2[j * 2 + k];	 	  	
	}   
      break;
    case INNER_INNER:
      for (i = 0; i < n; i++) 
	{     	      
	  x1  = &x1_start[8 * i];
	  x2  = &x2_start[8 * i];
	  sum = &sumtable[8 * i];
	  
	  for(j = 0; j < 4; j++)
	    for(k = 0; k < 2; k++)
	      sum[j * 2 + k] = x1[j * 2 + k] * x2[j * 2 + k];	 	  	  
	}
      break;
    default:
      assert(0);
    }
}

static void sumGAMMA(int tipCase, double *sumtable, double *x1_start, double *x2_start, double *tipVector,
		     unsigned char *tipX1, unsigned char *tipX2, int n)
{
  double *x1, *x2, *sum;
  int i, j, k;

  /* C-OPT once again switch over possible configurations at inner node */

  switch(tipCase)
    {
    case TIP_TIP:   
      /* C-OPT main for loop overt alignment length */
      for (i = 0; i < n; i++) 
	{     
	  x1 = &(tipVector[4 * tipX1[i]]);
	  x2 = &(tipVector[4 * tipX2[i]]);
	  sum = &sumtable[i * 16];
	  
	  for(j = 0; j < 4; j++)
	    for(k = 0; k < 4; k++)
	      sum[j * 4 + k] = x1[k] * x2[k];	           
	}    
      break;
    case TIP_INNER: 
      for (i = 0; i < n; i++) 
	{     
	  x1  = &(tipVector[4 * tipX1[i]]);
	  x2  = &x2_start[16 * i];
	  sum = &sumtable[16 * i];
	  
	  for(j = 0; j < 4; j++)
	    for(k = 0; k < 4; k++)
	      sum[j * 4 + k] = x1[k] * x2[j * 4 + k];	 	  	
	}   
      break;
    case INNER_INNER:
      for (i = 0; i < n; i++) 
	{     	      
	  x1  = &x1_start[16 * i];
	  x2  = &x2_start[16 * i];
	  sum = &sumtable[16 * i];
	  
	  for(j = 0; j < 4; j++)
	    for(k = 0; k < 4; k++)
	      sum[j * 4 + k] = x1[j * 4 + k] * x2[j * 4 + k];	 	  	  
	}
      break;
    default:
      assert(0);
    }
}


/* C-OPT the function below actually also computes the 1st and 2nd derivative of the 
   ML function because we need them for the Newton-Raphson procedure that we use to optimize
   branch lengths */

static void coreGTRGAMMA_BINARY(const int upper, double *sumtable, 
				volatile double *d1,   volatile double *d2, double *EIGN, double *gammaRates, double lz, int *wrptr)
{
  int i, j;
  double 
    *diagptable, *diagptable_start, *sum, 
    tmp_1, inv_Li, dlnLidlz, d2lnLidlz2, ki, kisqr,
    dlnLdlz = 0.0,
    d2lnLdlz2 = 0.0;

  diagptable = diagptable_start = (double *)malloc(sizeof(double) * 12);    
		 
  /* C-OPT some pre-computations */
  

  for(i = 0; i < 4; i++)
    {
      ki = gammaRates[i];	 
      kisqr = ki * ki;
	      
      diagptable[i * 3]     = exp (EIGN[0] * ki * lz);            
      diagptable[i * 3 + 1] = EIGN[0] * ki;
      diagptable[i * 3 + 2] = EIGN[0] * EIGN[0] * kisqr;            
    }


  /* C-OPT main fopr loop over enbtire alignment length, that where the bulk of the 
     computations takes place in this source file */

  for (i = 0; i < upper; i++) 
    {	   	    	   	    
      diagptable = diagptable_start;
      sum = &(sumtable[i * 8]);
      
      inv_Li      = 0.0;
      dlnLidlz    = 0.0;
      d2lnLidlz2  = 0.0;
      
      for(j = 0; j < 4; j++)
	{
	  inv_Li += sum[2 * j];
	  
	  tmp_1      =  diagptable[3 * j] * sum[2 * j + 1];
	  inv_Li     += tmp_1;
	  dlnLidlz   += tmp_1 * diagptable[3 * j + 1];
	  d2lnLidlz2 += tmp_1 * diagptable[3 * j + 2];	 	       	       	    	    
	}	   	   	  
      
      inv_Li = 1.0 / inv_Li;
      
      dlnLidlz   *= inv_Li;
      d2lnLidlz2 *= inv_Li;
      
    
      dlnLdlz  += wrptr[i] * dlnLidlz;
      d2lnLdlz2 += wrptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    }
    
  *d1 = dlnLdlz;
  *d2 = d2lnLdlz2;
  
  free(diagptable_start);
}



static void coreGTRGAMMA(const int upper, double *sumtable, 
			 volatile double *d1,   volatile double *d2, double *EIGN, double *gammaRates, double lz, int *wrptr)
{
  int i, j, k;
  double 
    *diagptable, *diagptable_start, *sum, 
    tmp_1, inv_Li, dlnLidlz, d2lnLidlz2, ki, kisqr,
    dlnLdlz = 0.0,
    d2lnLdlz2 = 0.0;

  diagptable = diagptable_start = (double *)malloc(sizeof(double) * 64);    
		 
  /* C-OPT some pre-computations */
  

  for(i = 0; i < 4; i++)
    {
      ki = gammaRates[i];	 
      kisqr = ki * ki;
	      
      diagptable[i * 16]     = exp (EIGN[0] * ki * lz);
      diagptable[i * 16 + 1] = exp (EIGN[1] * ki * lz);
      diagptable[i * 16 + 2] = exp (EIGN[2] * ki * lz);
      
      diagptable[i * 16 + 3] = EIGN[0] * ki;
      diagptable[i * 16 + 4] = EIGN[0] * EIGN[0] * kisqr;
      
      diagptable[i * 16 + 5] = EIGN[1] * ki;
      diagptable[i * 16 + 6] = EIGN[1] * EIGN[1] * kisqr;
      
      diagptable[i * 16 + 7] = EIGN[2] * ki;
      diagptable[i * 16 + 8] = EIGN[2] * EIGN[2] * kisqr;
    }


  /* C-OPT main fopr loop over enbtire alignment length, that where the bulk of the 
     computations takes place in this source file */

  for (i = 0; i < upper; i++) 
    {	   	    	   	    
      diagptable = diagptable_start;
      sum = &(sumtable[i * 16]);
      
      inv_Li      = 0.0;
      dlnLidlz    = 0.0;
      d2lnLidlz2  = 0.0;
      
      for(j = 0; j < 4; j++)
	{
	  inv_Li += sum[4 * j];
	  
	  for(k = 0; k < 3; k++)
	    {
	      tmp_1      =  diagptable[16 * j + k] * sum[4 * j + k + 1];
	      inv_Li     += tmp_1;
	      dlnLidlz   += tmp_1 * diagptable[16 * j + k * 2 + 3];
	      d2lnLidlz2 += tmp_1 * diagptable[16 * j + k * 2 + 4];
	    }	       	       	    	    
	}	   	   	  
      
      inv_Li = 1.0 / inv_Li;
      
      dlnLidlz   *= inv_Li;
      d2lnLidlz2 *= inv_Li;
      
    

      dlnLdlz  += wrptr[i] * dlnLidlz;
      d2lnLdlz2 += wrptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    }
    
  *d1 = dlnLdlz;
  *d2 = d2lnLdlz2;
  
  free(diagptable_start);
}









static void sumGAMMAPROT(int tipCase, double *sumtable, double *x1, double *x2, double *tipVector,
			 unsigned char *tipX1, unsigned char *tipX2, int n)
{
  int i, l, k;
  double *left, *right, *sum; 

  switch(tipCase)
    {
    case TIP_TIP:    
      for(i = 0; i < n; i++) 
	{     
	  left  = &(tipVector[20 * tipX1[i]]);
	  right = &(tipVector[20 * tipX2[i]]);      	

	  for(l = 0; l < 4; l++)
	    {
	      sum = &sumtable[i * 80 + l * 20];
	      for(k = 0; k < 20; k++)
		sum[k] = left[k] * right[k];	      	     	    	 
	    }
	}    
      break;
    case TIP_INNER:  
      for(i = 0; i < n; i++) 
	{     
	  left = &(tipVector[20 * tipX1[i]]);		    	  	 

	  for(l = 0; l < 4; l++)
	    {
	      right = &(x2[80 * i + l * 20]);
	      sum = &sumtable[i * 80 + l * 20];

	      for(k = 0; k < 20; k++)
		sum[k] = left[k] * right[k];	      	     
	    }	    	      	
	}    
      break;
    case INNER_INNER:
      for(i = 0; i < n; i++) 
	{     		    	  		  
	  for(l = 0; l < 4; l++)
	    {
	      left  = &(x1[80 * i + l * 20]);
	      right = &(x2[80 * i + l * 20]);
	      sum   = &(sumtable[i * 80 + l * 20]);

	      for(k = 0; k < 20; k++)
		sum[k] = left[k] * right[k];		      	     
	    }	  	
	}
      break;
    default:
      assert(0);
    }
}

static void sumGAMMASECONDARY(int tipCase, double *sumtable, double *x1, double *x2, double *tipVector,
			      unsigned char *tipX1, unsigned char *tipX2, int n)
{
  int i, l, k;
  double *left, *right, *sum; 

  switch(tipCase)
    {
    case TIP_TIP:    
      for(i = 0; i < n; i++) 
	{     
	  left  = &(tipVector[16 * tipX1[i]]);
	  right = &(tipVector[16 * tipX2[i]]);      	

	  for(l = 0; l < 4; l++)
	    {
	      sum = &sumtable[i * 64 + l * 16];
	      for(k = 0; k < 16; k++)
		sum[k] = left[k] * right[k];	      	     	    	 
	    }
	}    
      break;
    case TIP_INNER:  
      for(i = 0; i < n; i++) 
	{     
	  left = &(tipVector[16 * tipX1[i]]);		    	  	 

	  for(l = 0; l < 4; l++)
	    {
	      right = &(x2[64 * i + l * 16]);
	      sum = &sumtable[i * 64 + l * 16];

	      for(k = 0; k < 16; k++)
		sum[k] = left[k] * right[k];	      	     
	    }	    	      	
	}    
      break;
    case INNER_INNER:
      for(i = 0; i < n; i++) 
	{     		    	  		  
	  for(l = 0; l < 4; l++)
	    {
	      left  = &(x1[64 * i + l * 16]);
	      right = &(x2[64 * i + l * 16]);
	      sum   = &(sumtable[i * 64 + l * 16]);

	      for(k = 0; k < 16; k++)
		sum[k] = left[k] * right[k];		      	     
	    }	  	
	}
      break;
    default:
      assert(0);
    }
}

static void sumGAMMASECONDARY_6(int tipCase, double *sumtable, double *x1, double *x2, double *tipVector,
				unsigned char *tipX1, unsigned char *tipX2, int n)
{
  int i, l, k;
  double *left, *right, *sum; 

  switch(tipCase)
    {
    case TIP_TIP:    
      for(i = 0; i < n; i++) 
	{     
	  left  = &(tipVector[6 * tipX1[i]]);
	  right = &(tipVector[6 * tipX2[i]]);      	

	  for(l = 0; l < 4; l++)
	    {
	      sum = &sumtable[i * 24 + l * 6];
	      for(k = 0; k < 6; k++)
		sum[k] = left[k] * right[k];	      	     	    	 
	    }
	}    
      break;
    case TIP_INNER:  
      for(i = 0; i < n; i++) 
	{     
	  left = &(tipVector[6 * tipX1[i]]);		    	  	 

	  for(l = 0; l < 4; l++)
	    {
	      right = &(x2[24 * i + l * 6]);
	      sum = &sumtable[i * 24 + l * 6];

	      for(k = 0; k < 6; k++)
		sum[k] = left[k] * right[k];	      	     
	    }	    	      	
	}    
      break;
    case INNER_INNER:
      for(i = 0; i < n; i++) 
	{     		    	  		  
	  for(l = 0; l < 4; l++)
	    {
	      left  = &(x1[24 * i + l * 6]);
	      right = &(x2[24 * i + l * 6]);
	      sum   = &(sumtable[i * 24 + l * 6]);

	      for(k = 0; k < 6; k++)
		sum[k] = left[k] * right[k];		      	     
	    }	  	
	}
      break;
    default:
      assert(0);
    }
}

static void sumGAMMASECONDARY_7(int tipCase, double *sumtable, double *x1, double *x2, double *tipVector,
				unsigned char *tipX1, unsigned char *tipX2, int n)
{
  int i, l, k;
  double *left, *right, *sum; 

  switch(tipCase)
    {
    case TIP_TIP:    
      for(i = 0; i < n; i++) 
	{     
	  left  = &(tipVector[7 * tipX1[i]]);
	  right = &(tipVector[7 * tipX2[i]]);      	

	  for(l = 0; l < 4; l++)
	    {
	      sum = &sumtable[i * 28 + l * 7];
	      for(k = 0; k < 7; k++)
		sum[k] = left[k] * right[k];	      	     	    	 
	    }
	}    
      break;
    case TIP_INNER:  
      for(i = 0; i < n; i++) 
	{     
	  left = &(tipVector[7 * tipX1[i]]);		    	  	 

	  for(l = 0; l < 4; l++)
	    {
	      right = &(x2[28 * i + l * 7]);
	      sum = &sumtable[i * 28 + l * 7];

	      for(k = 0; k < 7; k++)
		sum[k] = left[k] * right[k];	      	     
	    }	    	      	
	}    
      break;
    case INNER_INNER:
      for(i = 0; i < n; i++) 
	{     		    	  		  
	  for(l = 0; l < 4; l++)
	    {
	      left  = &(x1[28 * i + l * 7]);
	      right = &(x2[28 * i + l * 7]);
	      sum   = &(sumtable[i * 28 + l * 7]);

	      for(k = 0; k < 7; k++)
		sum[k] = left[k] * right[k];		      	     
	    }	  	
	}
      break;
    default:
      assert(0);
    }
}


static void coreGTRGAMMAPROT(double *gammaRates, double *EIGN, double *sumtable, int upper, int *wrptr,
			      volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double lz)
{
  double  *sum, diagptable[512];
  int     i, j, l;
  double  dlnLdlz = 0;
  double d2lnLdlz2 = 0;
  double ki, kisqr; 
  double tmp;
  double inv_Li, dlnLidlz, d2lnLidlz2; 

  for(i = 0; i < 4; i++)
    {	   
      ki = gammaRates[i];	 
      kisqr = ki * ki;           

      for(l = 1; l < 20; l++)
	{	      
	  diagptable[i * 128 + l * 4]     = exp(EIGN[l-1] * ki * lz);
	  diagptable[i * 128 + l * 4 + 1] = EIGN[l-1] * ki;
	  diagptable[i * 128 + l * 4 + 2] = EIGN[l-1] * EIGN[l-1] * kisqr;
	}			    	      
    }

  for (i = 0; i < upper; i++) 
    {	  	    	        
      sum = &sumtable[i * 80];    
      inv_Li   = 0.0;
      dlnLidlz = 0.0;
      d2lnLidlz2 = 0.0;
        
      for(j = 0; j < 4; j++)
	{	  
	  inv_Li += sum[j * 20];     
      
	  for(l = 1; l < 20; l++)
	    {		
	      inv_Li     += (tmp = diagptable[j * 128 + l * 4] * sum[j * 20 + l]);
	      dlnLidlz   +=  tmp * diagptable[j * 128 + l * 4 + 1];
	      d2lnLidlz2 +=  tmp * diagptable[j * 128 + l * 4 + 2];		
	    }
	}     
      
      inv_Li = 1.0 / inv_Li;
      
      dlnLidlz   *= inv_Li;
      d2lnLidlz2 *= inv_Li;
      
      dlnLdlz   += wrptr[i] * dlnLidlz;
      d2lnLdlz2 += wrptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    } 

  *ext_dlnLdlz   = dlnLdlz;
  *ext_d2lnLdlz2 = d2lnLdlz2; 
}

static void coreGTRGAMMASECONDARY(double *gammaRates, double *EIGN, double *sumtable, int upper, int *wrptr,
				  volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double lz)
{
  double  *sum, diagptable[256];
  int     i, j, l;
  double  dlnLdlz = 0;
  double d2lnLdlz2 = 0;
  double ki, kisqr; 
  double tmp;
  double inv_Li, dlnLidlz, d2lnLidlz2; 

  for(i = 0; i < 4; i++)
    {	   
      ki = gammaRates[i];	 
      kisqr = ki * ki;           

      for(l = 1; l < 16; l++)
	{	      
	  diagptable[i * 64 + l * 4]     = exp(EIGN[l-1] * ki * lz);
	  diagptable[i * 64 + l * 4 + 1] = EIGN[l-1] * ki;
	  diagptable[i * 64 + l * 4 + 2] = EIGN[l-1] * EIGN[l-1] * kisqr;
	}			    	      
    }

  for (i = 0; i < upper; i++) 
    {	  	    	        
      sum = &sumtable[i * 64];    
      inv_Li   = 0.0;
      dlnLidlz = 0.0;
      d2lnLidlz2 = 0.0;
        
      for(j = 0; j < 4; j++)
	{	  
	  inv_Li += sum[j * 16];     
      
	  for(l = 1; l < 16; l++)
	    {		
	      inv_Li     += (tmp = diagptable[j * 64 + l * 4] * sum[j * 16 + l]);
	      dlnLidlz   +=  tmp * diagptable[j * 64 + l * 4 + 1];
	      d2lnLidlz2 +=  tmp * diagptable[j * 64 + l * 4 + 2];		
	    }
	}     
      
      inv_Li = 1.0 / inv_Li;
      
      dlnLidlz   *= inv_Li;
      d2lnLidlz2 *= inv_Li;
      
      dlnLdlz   += wrptr[i] * dlnLidlz;
      d2lnLdlz2 += wrptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    } 

  *ext_dlnLdlz   = dlnLdlz;
  *ext_d2lnLdlz2 = d2lnLdlz2; 
}

static void coreGTRGAMMASECONDARY_6(double *gammaRates, double *EIGN, double *sumtable, int upper, int *wrptr,
				    volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double lz)
{
  double  *sum, diagptable[96];
  int     i, j, l;
  double  dlnLdlz = 0;
  double d2lnLdlz2 = 0;
  double ki, kisqr; 
  double tmp;
  double inv_Li, dlnLidlz, d2lnLidlz2; 

  for(i = 0; i < 4; i++)
    {	   
      ki = gammaRates[i];	 
      kisqr = ki * ki;           

      for(l = 1; l < 6; l++)
	{	      
	  diagptable[i * 24 + l * 4]     = exp(EIGN[l-1] * ki * lz);
	  diagptable[i * 24 + l * 4 + 1] = EIGN[l-1] * ki;
	  diagptable[i * 24 + l * 4 + 2] = EIGN[l-1] * EIGN[l-1] * kisqr;
	}			    	      
    }

  for (i = 0; i < upper; i++) 
    {	  	    	        
      sum = &sumtable[i * 24];    
      inv_Li   = 0.0;
      dlnLidlz = 0.0;
      d2lnLidlz2 = 0.0;
        
      for(j = 0; j < 4; j++)
	{	  
	  inv_Li += sum[j * 6];     
      
	  for(l = 1; l < 6; l++)
	    {		
	      inv_Li     += (tmp = diagptable[j * 24 + l * 4] * sum[j * 6 + l]);
	      dlnLidlz   +=  tmp * diagptable[j * 24 + l * 4 + 1];
	      d2lnLidlz2 +=  tmp * diagptable[j * 24 + l * 4 + 2];		
	    }
	}     
      
      inv_Li = 1.0 / inv_Li;
      
      dlnLidlz   *= inv_Li;
      d2lnLidlz2 *= inv_Li;
      
      dlnLdlz   += wrptr[i] * dlnLidlz;
      d2lnLdlz2 += wrptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    } 

  *ext_dlnLdlz   = dlnLdlz;
  *ext_d2lnLdlz2 = d2lnLdlz2; 
}

static void coreGTRGAMMASECONDARY_7(double *gammaRates, double *EIGN, double *sumtable, int upper, int *wrptr,
				    volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double lz)
{
  double  *sum, diagptable[112];
  int     i, j, l;
  double  dlnLdlz = 0;
  double d2lnLdlz2 = 0;
  double ki, kisqr; 
  double tmp;
  double inv_Li, dlnLidlz, d2lnLidlz2; 

  for(i = 0; i < 4; i++)
    {	   
      ki = gammaRates[i];	 
      kisqr = ki * ki;           

      for(l = 1; l < 7; l++)
	{	      
	  diagptable[i * 28 + l * 4]     = exp(EIGN[l-1] * ki * lz);
	  diagptable[i * 28 + l * 4 + 1] = EIGN[l-1] * ki;
	  diagptable[i * 28 + l * 4 + 2] = EIGN[l-1] * EIGN[l-1] * kisqr;
	}			    	      
    }

  for (i = 0; i < upper; i++) 
    {	  	    	        
      sum = &sumtable[i * 28];    
      inv_Li   = 0.0;
      dlnLidlz = 0.0;
      d2lnLidlz2 = 0.0;
        
      for(j = 0; j < 4; j++)
	{	  
	  inv_Li += sum[j * 7];     
      
	  for(l = 1; l < 7; l++)
	    {		
	      inv_Li     += (tmp = diagptable[j * 28 + l * 4] * sum[j * 7 + l]);
	      dlnLidlz   +=  tmp * diagptable[j * 28 + l * 4 + 1];
	      d2lnLidlz2 +=  tmp * diagptable[j * 28 + l * 4 + 2];		
	    }
	}     
      
      inv_Li = 1.0 / inv_Li;
      
      dlnLidlz   *= inv_Li;
      d2lnLidlz2 *= inv_Li;
      
      dlnLdlz   += wrptr[i] * dlnLidlz;
      d2lnLdlz2 += wrptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    } 

  *ext_dlnLdlz   = dlnLdlz;
  *ext_d2lnLdlz2 = d2lnLdlz2; 
}

static void coreGTRGAMMAPROTINVAR(double *gammaRates, double *EIGN, double *sumtable, int upper, int *wrptr,
				   volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double lz, double *frequencies,
				  double propInvar, int *iptr)
{
  double  *sum, diagptable[512];
  int     i, l, j;
  double  dlnLdlz = 0;
  double d2lnLdlz2 = 0;
  double ki, kisqr;    
  double freqs[20]; 
  double scaler =  0.25 * (1.0 - propInvar);
  double tmp;
  double inv_Li, dlnLidlz, d2lnLidlz2;

  for(i = 0; i < 20; i++)    
    freqs[i] = frequencies[i] * propInvar;
  
  for(i = 0; i < 4; i++)
    {	   
      ki = gammaRates[i];	 
      kisqr = ki * ki;           

      for(l = 1; l < 20; l++)
	{	      
	  diagptable[i * 128 + l * 4]     = exp(EIGN[l-1] * ki * lz);
	  diagptable[i * 128 + l * 4 + 1] = EIGN[l-1] * ki;
	  diagptable[i * 128 + l * 4 + 2] = EIGN[l-1] * EIGN[l-1] * kisqr;
	}			    	      
    }

   
   for(i = 0; i < upper; i++) 
     {	  	    	          
       sum = &sumtable[i * 80];    
       inv_Li   = 0.0;
       dlnLidlz = 0.0;
       d2lnLidlz2 = 0.0;
       
       for(j = 0; j < 4; j++)
	{	  
	  inv_Li += sum[j * 20];     
      
	  for(l = 1; l < 20; l++)
	    {		
	      inv_Li     += (tmp = diagptable[j * 128 + l * 4] * sum[j * 20 + l]);
	      dlnLidlz   +=  tmp * diagptable[j * 128 + l * 4 + 1];
	      d2lnLidlz2 +=  tmp * diagptable[j * 128 + l * 4 + 2];		
	    }
	}     
       
       inv_Li *= scaler;
       
       if(iptr[i] < 20)	      
	 inv_Li += freqs[iptr[i]];
       
       inv_Li = 1.0 / inv_Li;
       
       dlnLidlz   *= inv_Li;
       d2lnLidlz2 *= inv_Li;
       
       dlnLidlz *= scaler;
       d2lnLidlz2 *= scaler;
       
       dlnLdlz  += wrptr[i] * dlnLidlz;
       d2lnLdlz2 += wrptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
     }

   *ext_dlnLdlz   = dlnLdlz;
   *ext_d2lnLdlz2 = d2lnLdlz2;     	  
}


static void coreGTRGAMMASECONDARYINVAR(double *gammaRates, double *EIGN, double *sumtable, int upper, int *wrptr,
				       volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double lz, double *frequencies,
				       double propInvar, int *iptr)
{
  double  *sum, diagptable[256];
  int     i, l, j;
  double  dlnLdlz = 0;
  double d2lnLdlz2 = 0;
  double ki, kisqr;    
  double freqs[16]; 
  double scaler =  0.25 * (1.0 - propInvar);
  double tmp;
  double inv_Li, dlnLidlz, d2lnLidlz2;

  for(i = 0; i < 16; i++)    
    freqs[i] = frequencies[i] * propInvar;
  
  for(i = 0; i < 4; i++)
    {	   
      ki = gammaRates[i];	 
      kisqr = ki * ki;           

      for(l = 1; l < 16; l++)
	{	      
	  diagptable[i * 64 + l * 4]     = exp(EIGN[l-1] * ki * lz);
	  diagptable[i * 64 + l * 4 + 1] = EIGN[l-1] * ki;
	  diagptable[i * 64 + l * 4 + 2] = EIGN[l-1] * EIGN[l-1] * kisqr;
	}			    	      
    }

   
   for(i = 0; i < upper; i++) 
     {	  	    	          
       sum = &sumtable[i * 64];    
       inv_Li   = 0.0;
       dlnLidlz = 0.0;
       d2lnLidlz2 = 0.0;
       
       for(j = 0; j < 4; j++)
	{	  
	  inv_Li += sum[j * 64];     
      
	  for(l = 1; l < 16; l++)
	    {		
	      inv_Li     += (tmp = diagptable[j * 64 + l * 4] * sum[j * 16 + l]);
	      dlnLidlz   +=  tmp * diagptable[j * 64 + l * 4 + 1];
	      d2lnLidlz2 +=  tmp * diagptable[j * 64 + l * 4 + 2];		
	    }
	}     
       
       inv_Li *= scaler;
       
       if(iptr[i] < 16)	      
	 inv_Li += freqs[iptr[i]];
       
       inv_Li = 1.0 / inv_Li;
       
       dlnLidlz   *= inv_Li;
       d2lnLidlz2 *= inv_Li;
       
       dlnLidlz *= scaler;
       d2lnLidlz2 *= scaler;
       
       dlnLdlz  += wrptr[i] * dlnLidlz;
       d2lnLdlz2 += wrptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
     }

   *ext_dlnLdlz   = dlnLdlz;
   *ext_d2lnLdlz2 = d2lnLdlz2;     	  
}

static void coreGTRGAMMASECONDARYINVAR_6(double *gammaRates, double *EIGN, double *sumtable, int upper, int *wrptr,
					 volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double lz, double *frequencies,
					 double propInvar, int *iptr)
{
  double  *sum, diagptable[96];
  int     i, l, j;
  double  dlnLdlz = 0;
  double d2lnLdlz2 = 0;
  double ki, kisqr;    
  double freqs[6]; 
  double scaler =  0.25 * (1.0 - propInvar);
  double tmp;
  double inv_Li, dlnLidlz, d2lnLidlz2;

  for(i = 0; i < 6; i++)    
    freqs[i] = frequencies[i] * propInvar;
  
  for(i = 0; i < 4; i++)
    {	   
      ki = gammaRates[i];	 
      kisqr = ki * ki;           

      for(l = 1; l < 6; l++)
	{	      
	  diagptable[i * 24 + l * 4]     = exp(EIGN[l-1] * ki * lz);
	  diagptable[i * 24 + l * 4 + 1] = EIGN[l-1] * ki;
	  diagptable[i * 24 + l * 4 + 2] = EIGN[l-1] * EIGN[l-1] * kisqr;
	}			    	      
    }

   
   for(i = 0; i < upper; i++) 
     {	  	    	          
       sum = &sumtable[i * 24];    
       inv_Li   = 0.0;
       dlnLidlz = 0.0;
       d2lnLidlz2 = 0.0;
       
       for(j = 0; j < 4; j++)
	{	  
	  inv_Li += sum[j * 24];     
      
	  for(l = 1; l < 6; l++)
	    {		
	      inv_Li     += (tmp = diagptable[j * 24 + l * 4] * sum[j * 6 + l]);
	      dlnLidlz   +=  tmp * diagptable[j * 24 + l * 4 + 1];
	      d2lnLidlz2 +=  tmp * diagptable[j * 24 + l * 4 + 2];		
	    }
	}     
       
       inv_Li *= scaler;
       
       if(iptr[i] < 6)	      
	 inv_Li += freqs[iptr[i]];
       
       inv_Li = 1.0 / inv_Li;
       
       dlnLidlz   *= inv_Li;
       d2lnLidlz2 *= inv_Li;
       
       dlnLidlz *= scaler;
       d2lnLidlz2 *= scaler;
       
       dlnLdlz  += wrptr[i] * dlnLidlz;
       d2lnLdlz2 += wrptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
     }

   *ext_dlnLdlz   = dlnLdlz;
   *ext_d2lnLdlz2 = d2lnLdlz2;     	  
}

static void coreGTRGAMMASECONDARYINVAR_7(double *gammaRates, double *EIGN, double *sumtable, int upper, int *wrptr,
					 volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double lz, double *frequencies,
					 double propInvar, int *iptr)
{
  double  *sum, diagptable[112];
  int     i, l, j;
  double  dlnLdlz = 0;
  double d2lnLdlz2 = 0;
  double ki, kisqr;    
  double freqs[7]; 
  double scaler =  0.25 * (1.0 - propInvar);
  double tmp;
  double inv_Li, dlnLidlz, d2lnLidlz2;

  for(i = 0; i < 7; i++)    
    freqs[i] = frequencies[i] * propInvar;
  
  for(i = 0; i < 4; i++)
    {	   
      ki = gammaRates[i];	 
      kisqr = ki * ki;           

      for(l = 1; l < 7; l++)
	{	      
	  diagptable[i * 28 + l * 4]     = exp(EIGN[l-1] * ki * lz);
	  diagptable[i * 28 + l * 4 + 1] = EIGN[l-1] * ki;
	  diagptable[i * 28 + l * 4 + 2] = EIGN[l-1] * EIGN[l-1] * kisqr;
	}			    	      
    }

   
   for(i = 0; i < upper; i++) 
     {	  	    	          
       sum = &sumtable[i * 28];    
       inv_Li   = 0.0;
       dlnLidlz = 0.0;
       d2lnLidlz2 = 0.0;
       
       for(j = 0; j < 4; j++)
	{	  
	  inv_Li += sum[j * 28];     
      
	  for(l = 1; l < 7; l++)
	    {		
	      inv_Li     += (tmp = diagptable[j * 28 + l * 4] * sum[j * 7 + l]);
	      dlnLidlz   +=  tmp * diagptable[j * 28 + l * 4 + 1];
	      d2lnLidlz2 +=  tmp * diagptable[j * 28 + l * 4 + 2];		
	    }
	}     
       
       inv_Li *= scaler;
       
       if(iptr[i] < 7)	      
	 inv_Li += freqs[iptr[i]];
       
       inv_Li = 1.0 / inv_Li;
       
       dlnLidlz   *= inv_Li;
       d2lnLidlz2 *= inv_Li;
       
       dlnLidlz *= scaler;
       d2lnLidlz2 *= scaler;
       
       dlnLdlz  += wrptr[i] * dlnLidlz;
       d2lnLdlz2 += wrptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
     }

   *ext_dlnLdlz   = dlnLdlz;
   *ext_d2lnLdlz2 = d2lnLdlz2;     	  
}


static void getVects(tree *tr, unsigned char **tipX1, unsigned char **tipX2, double **x1_start, double **x2_start, int *tipCase, int model)
{
  int pNumber = tr->td[0].ti[0].pNumber;
  int qNumber = tr->td[0].ti[0].qNumber;

  *x1_start = (double*)NULL, 
  *x2_start = (double*)NULL;
  *tipX1 = (unsigned char*)NULL,
  *tipX2 = (unsigned char*)NULL;

  if(isTip(pNumber, tr->mxtips) || isTip(qNumber, tr->mxtips))
    {	    	    	   	               
      if(!( isTip(pNumber, tr->mxtips) && isTip(qNumber, tr->mxtips)) )
	{	 
	  *tipCase = TIP_INNER;
	  if(isTip(qNumber, tr->mxtips))
	    {	
	      *tipX1 = tr->partitionData[model].yVector[qNumber];                    
	      *x2_start = tr->partitionData[model].xVector[pNumber - tr->mxtips - 1];
	      
	    }	    
	  else
	    {
	      *tipX1 = tr->partitionData[model].yVector[pNumber];                    
	      *x2_start = tr->partitionData[model].xVector[qNumber - tr->mxtips - 1];
	    }
	}
      else
	{
	  *tipCase = TIP_TIP;
	  *tipX1 = tr->partitionData[model].yVector[pNumber];
	  *tipX2 = tr->partitionData[model].yVector[qNumber];
	}
    }
  else
    {
      *tipCase = INNER_INNER;
      
      *x1_start = tr->partitionData[model].xVector[pNumber - tr->mxtips - 1];
      *x2_start = tr->partitionData[model].xVector[qNumber - tr->mxtips - 1];
    }

}
  



void makenewzIterative(tree *tr)
{  
  int model, tipCase;  
  double   
    *x1_start, 
    *x2_start;
  unsigned char 
    *tipX1,
    *tipX2;       

  newviewIterative(tr);  
   
  for(model = 0; model < tr->NumberOfModels; model++)
    {		
      if(tr->executeModel[model])
	{
	  int width = tr->partitionData[model].width;
	  
	  getVects(tr, &tipX1, &tipX2, &x1_start, &x2_start, &tipCase, model);

	  switch(tr->partitionData[model].dataType)
	    { 
	    case BINARY_DATA:
	      switch(tr->rateHetModel)
		{
		case CAT:
		  sumCAT_BINARY(tipCase, tr->partitionData[model].sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector, tipX1, tipX2,
				width);
		  break;
		case GAMMA:		  
		case GAMMA_I:
		  sumGAMMA_BINARY(tipCase, tr->partitionData[model].sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector, tipX1, tipX2, 
				  width);
		  break;
		default:
		  assert(0);
		}
	      break;
	    case DNA_DATA:
	      switch(tr->rateHetModel)
		{
		case CAT:
		  sumCAT(tipCase, tr->partitionData[model].sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector, tipX1, tipX2,
			 width);
		  break;
		case GAMMA:		  
		case GAMMA_I:
		  sumGAMMA(tipCase, tr->partitionData[model].sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector, tipX1, tipX2, 
			   width);
		  break;
		default:
		  assert(0);
		}
	      break;
	    case AA_DATA:
	      switch(tr->rateHetModel)
		{
		case CAT:
		  sumGTRCATPROT(tipCase, tr->partitionData[model].sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector, 
				tipX1, tipX2, width);
		  break;
		case GAMMA:
		case GAMMA_I:
		  sumGAMMAPROT(tipCase,  tr->partitionData[model].sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector, 
			       tipX1, tipX2, width);
		  break;
		default:
		  assert(0);
		}
	      break;
	    case SECONDARY_DATA:
	      switch(tr->rateHetModel)
	      	{
		case CAT:
		  sumGTRCATSECONDARY(tipCase, tr->partitionData[model].sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector, 
				     tipX1, tipX2, width);
		  break;
		case GAMMA:
		case GAMMA_I:
		  sumGAMMASECONDARY(tipCase,  tr->partitionData[model].sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector, 
				    tipX1, tipX2, width);
		  break;
		default:
		  assert(0);
		}
	      break;
	    case SECONDARY_DATA_6:
	      switch(tr->rateHetModel)
	      	{
		case CAT:
		  sumGTRCATSECONDARY_6(tipCase, tr->partitionData[model].sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector, 
				       tipX1, tipX2, width);
		  break;
		case GAMMA:
		case GAMMA_I:
		  sumGAMMASECONDARY_6(tipCase,  tr->partitionData[model].sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector, 
				      tipX1, tipX2, width);
		  break;
		default:
		  assert(0);
		}
	      break;  
	    case SECONDARY_DATA_7:
	      switch(tr->rateHetModel)
	      	{
		case CAT:
		  sumGTRCATSECONDARY_7(tipCase, tr->partitionData[model].sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector, 
				       tipX1, tipX2, width);
		  break;
		case GAMMA:
		case GAMMA_I:
		  sumGAMMASECONDARY_7(tipCase,  tr->partitionData[model].sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector, 
				      tipX1, tipX2, width);
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





void execCore(tree *tr, volatile double *_dlnLdlz, volatile double *_d2lnLdlz2)
{ 
  int model, branchIndex;
  double lz;

  _dlnLdlz[0]   = 0.0;
  _d2lnLdlz2[0] = 0.0;  
 
  for(model = 0; model < tr->NumberOfModels; model++)
    {
      if(tr->executeModel[model])
	{
	  int width = tr->partitionData[model].width;
	  double *sumBuffer = tr->partitionData[model].sumBuffer;
	  volatile double
	    dlnLdlz   = 0.0,
	    d2lnLdlz2 = 0.0;

	  if(tr->multiBranch)
	    {
	      branchIndex = model;
	      lz = tr->coreLZ[model];
	      _dlnLdlz[model]   = 0.0;
	      _d2lnLdlz2[model] = 0.0;	    
	    }
	  else
	    {
	      branchIndex = 0;
	      lz = tr->coreLZ[0];	 	     
	    }	 
	     	 	      	 
	  switch(tr->partitionData[model].dataType)
	    { 
	    case BINARY_DATA:
	      switch(tr->rateHetModel)
		{
		case CAT:
		  coreGTRCAT_BINARY(width, tr->NumberOfCategories, sumBuffer, 
				    &dlnLdlz, &d2lnLdlz2, tr->partitionData[model].wr, tr->partitionData[model].wr2,
				    tr->cdta->patrat, tr->partitionData[model].EIGN,  tr->partitionData[model].rateCategory, lz);
		  break;
		case GAMMA:	
		  coreGTRGAMMA_BINARY(width, sumBuffer, 
				      &dlnLdlz, &d2lnLdlz2, tr->partitionData[model].EIGN, tr->partitionData[model].gammaRates, lz, 
				      tr->partitionData[model].wgt);		   
		  break;
		case GAMMA_I:
		  coreGTRGAMMAINVAR_BINARY(tr->partitionData[model].propInvariant, tr->partitionData[model].frequencies, 
					   tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN,
					   sumBuffer, &dlnLdlz, &d2lnLdlz2,
					   tr->partitionData[model].invariant, tr->partitionData[model].wgt, 
					   width, lz);
		  break;
		default:
		  assert(0);
		}
	      break;
	    case DNA_DATA:
	      switch(tr->rateHetModel)
		{
		case CAT:
		  coreGTRCAT(width, tr->NumberOfCategories, sumBuffer, 
			     &dlnLdlz, &d2lnLdlz2, tr->partitionData[model].wr, tr->partitionData[model].wr2,
			     tr->cdta->patrat, tr->partitionData[model].EIGN,  tr->partitionData[model].rateCategory, lz);
		  break;
		case GAMMA:	
		  coreGTRGAMMA(width, sumBuffer, 
			       &dlnLdlz, &d2lnLdlz2, tr->partitionData[model].EIGN, tr->partitionData[model].gammaRates, lz, 
			       tr->partitionData[model].wgt);		      
		  break;
		case GAMMA_I:
		  coreGTRGAMMAINVAR(tr->partitionData[model].propInvariant, tr->partitionData[model].frequencies, 
				    tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN,
				    sumBuffer, &dlnLdlz, &d2lnLdlz2,
				    tr->partitionData[model].invariant, tr->partitionData[model].wgt, 
				    width, lz);
		  break;
		default:
		  assert(0);
		}
	      break;
	    case AA_DATA:
	      switch(tr->rateHetModel)
		{
		case CAT:
		  coreGTRCATPROT(tr->partitionData[model].EIGN, lz, tr->NumberOfCategories,  &(tr->cdta->patrat[0]), 
				 tr->partitionData[model].rateCategory, width,
				 tr->partitionData[model].wr, tr->partitionData[model].wr2, 
				 &dlnLdlz, &d2lnLdlz2, 
				 sumBuffer);
		  break;
		case GAMMA:
		  coreGTRGAMMAPROT(tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, 
				   sumBuffer, width, tr->partitionData[model].wgt,
				   &dlnLdlz, &d2lnLdlz2, lz);
		  break;
		case GAMMA_I:
		  coreGTRGAMMAPROTINVAR(tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, 
					sumBuffer, width, tr->partitionData[model].wgt,
					&dlnLdlz, &d2lnLdlz2, lz, tr->partitionData[model].frequencies,
					tr->partitionData[model].propInvariant, tr->partitionData[model].invariant);
		  break;
		default:
		  assert(0);
		}
	      break;
	    case SECONDARY_DATA:
	      switch(tr->rateHetModel)
		{
		case CAT:
		  coreGTRCATSECONDARY(tr->partitionData[model].EIGN, lz, tr->NumberOfCategories,  &(tr->cdta->patrat[0]), 
				      tr->partitionData[model].rateCategory, width,
				      tr->partitionData[model].wr, tr->partitionData[model].wr2, 
				      &dlnLdlz, &d2lnLdlz2, 
				      sumBuffer);
		  break;
		case GAMMA:
		  coreGTRGAMMASECONDARY(tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, 
					sumBuffer, width, tr->partitionData[model].wgt,
					&dlnLdlz, &d2lnLdlz2, lz);
		  break;
		case GAMMA_I:
		  coreGTRGAMMASECONDARYINVAR(tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, 
					     sumBuffer, width, tr->partitionData[model].wgt,
					     &dlnLdlz, &d2lnLdlz2, lz, tr->partitionData[model].frequencies,
					     tr->partitionData[model].propInvariant, tr->partitionData[model].invariant);
		  break;
		default:
		  assert(0);
		}
	      break;
	    case SECONDARY_DATA_6:
	      switch(tr->rateHetModel)
		{
		case CAT:
		  coreGTRCATSECONDARY_6(tr->partitionData[model].EIGN, lz, tr->NumberOfCategories,  &(tr->cdta->patrat[0]), 
					tr->partitionData[model].rateCategory, width,
					tr->partitionData[model].wr, tr->partitionData[model].wr2, 
					&dlnLdlz, &d2lnLdlz2, 
					sumBuffer);
		  break;
		case GAMMA:
		  coreGTRGAMMASECONDARY_6(tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, 
					  sumBuffer, width, tr->partitionData[model].wgt,
					  &dlnLdlz, &d2lnLdlz2, lz);
		  break;
		case GAMMA_I:
		  coreGTRGAMMASECONDARYINVAR_6(tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, 
					       sumBuffer, width, tr->partitionData[model].wgt,
					       &dlnLdlz, &d2lnLdlz2, lz, tr->partitionData[model].frequencies,
					       tr->partitionData[model].propInvariant, tr->partitionData[model].invariant);
		  break;
		default:
		  assert(0);
		}
	      break;
	    case SECONDARY_DATA_7:
	      switch(tr->rateHetModel)
		{
		case CAT:
		  coreGTRCATSECONDARY_7(tr->partitionData[model].EIGN, lz, tr->NumberOfCategories,  &(tr->cdta->patrat[0]), 
					tr->partitionData[model].rateCategory, width,
					tr->partitionData[model].wr, tr->partitionData[model].wr2, 
					&dlnLdlz, &d2lnLdlz2, 
					sumBuffer);
		  break;
		case GAMMA:
		  coreGTRGAMMASECONDARY_7(tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, 
					  sumBuffer, width, tr->partitionData[model].wgt,
					  &dlnLdlz, &d2lnLdlz2, lz);
		  break;
		case GAMMA_I:
		  coreGTRGAMMASECONDARYINVAR_7(tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, 
					       sumBuffer, width, tr->partitionData[model].wgt,
					       &dlnLdlz, &d2lnLdlz2, lz, tr->partitionData[model].frequencies,
					       tr->partitionData[model].propInvariant, tr->partitionData[model].invariant);
		  break;
		default:
		  assert(0);
		}
	      break;
	    default:
	      assert(0);
	    }
	  
	  _dlnLdlz[branchIndex]   = _dlnLdlz[branchIndex]   + dlnLdlz;
	  _d2lnLdlz2[branchIndex] = _d2lnLdlz2[branchIndex] + d2lnLdlz2;	
	}
    }
	            
}


static void topLevelMakenewz(tree *tr, double *z0, int _maxiter, double *result)
{
  double   z[NUM_BRANCHES], zprev[NUM_BRANCHES], zstep[NUM_BRANCHES];    
  volatile double  dlnLdlz[NUM_BRANCHES], d2lnLdlz2[NUM_BRANCHES];
  int i, maxiter[NUM_BRANCHES], numBranches, model;
  boolean firstIteration = TRUE;
  boolean outerConverged[NUM_BRANCHES]; 
  boolean loopConverged;



  if(tr->multiBranch)
    numBranches = tr->NumberOfModels;
  else
    numBranches = 1;
 


  for(i = 0; i < numBranches; i++)
    {
      z[i] = z0[i];
      maxiter[i] = _maxiter;
      outerConverged[i] = FALSE;
      tr->curvatOK[i]       = TRUE;       
    }
      
  do 
    {
      for(i = 0; i < numBranches; i++)
	{
	  if(outerConverged[i] == FALSE && tr->curvatOK[i] == TRUE)
	    {
	      tr->curvatOK[i] = FALSE;
	            
	      zprev[i] = z[i];
	      
	      zstep[i] = (1.0 - zmax) * z[i] + zmin;
	    }
	}
	 
      for(i = 0; i < numBranches; i++)
	{
	  if(outerConverged[i] == FALSE && tr->curvatOK[i] == FALSE)
	    {
	      double lz;	  	  	    	               
	    	    
	      if (z[i] < zmin) z[i] = zmin;
	      else if (z[i] > zmax) z[i] = zmax;
	      lz    = log(z[i]);
	    	 	 	    	  
	      tr->coreLZ[i] = lz;
	    }
	}

      if(tr->multiBranch)
	{
	  for(model = 0; model < tr->NumberOfModels; model++)
	    {	     
	      if(tr->executeModel[model])
		tr->executeModel[model] = !tr->curvatOK[model];	    
	      
	    }
	}
      else
	{
	  for(model = 0; model < tr->NumberOfModels; model++)	  
	    tr->executeModel[model] = !tr->curvatOK[0];	     
	}


#ifdef _USE_PTHREADS 
      if(firstIteration)
	{
	  masterBarrier(THREAD_MAKENEWZ_FIRST, tr);
	  firstIteration = FALSE;
	}
      else
	masterBarrier(THREAD_MAKENEWZ, tr);

      if(!tr->multiBranch)
	{
	  dlnLdlz[0] = 0.0;
	  d2lnLdlz2[0] = 0.0;
	  for(i = 0; i < NumberOfThreads; i++)
	    {
	      dlnLdlz[0]   += reductionBuffer[i];
	      d2lnLdlz2[0] += reductionBufferTwo[i];
	    }
	}
      else
	{
	  int j;
	  for(j = 0; j < tr->NumberOfModels; j++)
	    {
	      dlnLdlz[j] = 0.0;
	      d2lnLdlz2[j] = 0.0;
	      for(i = 0; i < NumberOfThreads; i++)
		{
		  dlnLdlz[j]   += reductionBuffer[i * tr->NumberOfModels + j];
		  d2lnLdlz2[j] += reductionBufferTwo[i * tr->NumberOfModels + j];
		}
	    }
	}
#else     	  
      if(firstIteration)
	{
	  makenewzIterative(tr);
	  firstIteration = FALSE;
	}
      execCore(tr, dlnLdlz, d2lnLdlz2);
#endif	    

      for(i = 0; i < numBranches; i++)
	{
	  if(outerConverged[i] == FALSE && tr->curvatOK[i] == FALSE)
	    {
	      if ((d2lnLdlz2[i] >= 0.0) && (z[i] < zmax))
		zprev[i] = z[i] = 0.37 * z[i] + 0.63;  /*  Bad curvature, shorten branch */
	      else
		tr->curvatOK[i] = TRUE;
	    }
	}
	    
      for(i = 0; i < numBranches; i++)
	{
	  if(tr->curvatOK[i] == TRUE && outerConverged[i] == FALSE)
	    {
	      if (d2lnLdlz2[i] < 0.0) 
		{
		  double tantmp = -dlnLdlz[i] / d2lnLdlz2[i]; 
		  if (tantmp < 100) 
		    {
		      z[i] *= exp(tantmp);
		      if (z[i] < zmin) 			
			z[i] = zmin;	    
		
		      if (z[i] > 0.25 * zprev[i] + 0.75) 
			z[i] = 0.25 * zprev[i] + 0.75;
		    } 
		  else 		  
		    z[i] = 0.25 * zprev[i] + 0.75;		   
		}
	      if (z[i] > zmax) z[i] = zmax;
	      
	      maxiter[i] = maxiter[i] - 1;
	      if(maxiter[i] > 0 && (ABS(z[i] - zprev[i]) > zstep[i]))
		outerConverged[i] = FALSE;
	      else
		outerConverged[i] = TRUE;
	    }
	}
	      
      loopConverged = TRUE;
      for(i = 0; i < numBranches; i++)
	loopConverged = loopConverged && outerConverged[i]; 
    } 
  while (!loopConverged);
  

  for(model = 0; model < tr->NumberOfModels; model++)
    tr->executeModel[model] = TRUE;	

  for(i = 0; i < numBranches; i++)  
    result[i] = z[i];
}
			     



void makenewzGeneric(tree *tr, nodeptr p, nodeptr q, double *z0, int maxiter, double *result, boolean mask)
{      
  int i;
  boolean originalExecute[NUM_BRANCHES];
  
  tr->td[0].ti[0].pNumber = p->number;
  tr->td[0].ti[0].qNumber = q->number;        
  for(i = 0; i < tr->numBranches; i++)      
    {
      originalExecute[i] =  tr->executeModel[i];
      tr->td[0].ti[0].qz[i] =  z0[i]; 
      if(mask)
	{
	  if(tr->partitionConverged[i])
	    tr->executeModel[i] = FALSE;
	  else
	    tr->executeModel[i] = TRUE;		
	}
    }

  tr->td[0].count = 1;

  if(!p->x)
    computeTraversalInfo(p, &(tr->td[0].ti[0]), &(tr->td[0].count), tr->mxtips, tr->numBranches);
  if(!q->x)
    computeTraversalInfo(q, &(tr->td[0].ti[0]), &(tr->td[0].count), tr->mxtips, tr->numBranches);        
 
  topLevelMakenewz(tr, z0, maxiter, result);
  
  for(i = 0; i < tr->numBranches; i++)          
      tr->executeModel[i] = TRUE;
}



void makenewzGenericDistance(tree *tr, int maxiter, double *z0, double *result, int taxon1, int taxon2)
{      
  int i;
  
  assert(taxon1 != taxon2);
  assert(0 < taxon1 && taxon1 <= tr->mxtips);
  assert(0 < taxon2 && taxon2 <= tr->mxtips);

  tr->td[0].ti[0].pNumber = taxon1;
  tr->td[0].ti[0].qNumber = taxon2;   
  tr->td[0].ti[0].tipCase = TIP_TIP;
   
  for(i = 0; i < tr->numBranches; i++)          
    tr->td[0].ti[0].qz[i] =  defaultz;      /* TODO why defaultz ? */

  tr->td[0].count = 1;        
 
  topLevelMakenewz(tr, z0, maxiter, result);
}



