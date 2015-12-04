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

#ifdef WIN32
#include <direct.h>
#endif

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
#include <stdarg.h>


#ifdef PARALLEL
#include <mpi.h>
#endif



#ifdef _USE_PTHREADS
#include <pthread.h>
#endif


#include "axml.h"
#include "globalVariables.h"


/***************** UTILITY FUNCTIONS **************************/


static void printBoth(FILE *f, const char* format, ... ) 
{
  va_list args;
  va_start(args, format);
  vfprintf(f, format, args );
  va_end(args);

  va_start(args, format);
  vprintf(format, args );
  va_end(args);
}

void printBothOpen(const char* format, ... ) 
{
  FILE *f = myfopen(infoFileName, "a");

  va_list args;
  va_start(args, format);
  vfprintf(f, format, args );
  va_end(args);

  va_start(args, format);
  vprintf(format, args );
  va_end(args);

  fclose(f);
}




static boolean isCat(analdef *adef)
{
  if(adef->model == M_PROTCAT || adef->model == M_GTRCAT || adef->model == M_BINCAT)
    return TRUE;
  else
    return FALSE;
}

static boolean isGamma(analdef *adef)
{
  if(adef->model == M_PROTGAMMA || adef->model == M_GTRGAMMA || adef->model == M_BINGAMMA)
    return TRUE;
  else
    return FALSE;

}

static void setRateHetAndDataIncrement(tree *tr, analdef *adef)
{
  if(isCat(adef))
    tr->rateHetModel = CAT;
  else
    {
      if(adef->useInvariant)
	tr->rateHetModel = GAMMA_I;
      else
	tr->rateHetModel = GAMMA;
    }

  switch(tr->rateHetModel)
    {
    case GAMMA:
    case GAMMA_I:
      tr->dnaIncrement        = 16;
      tr->aaIncrement         = 80;
      tr->binaryIncrement     = 8;
      tr->secondaryIncrement  = 64;
      tr->secondaryIncrement6 = 24;
      tr->secondaryIncrement7 = 28;
      break;
    case CAT:
      if((adef->boot && !adef->bootstrapBranchLengths) || (adef->classifyML))
	{
	  tr->binaryIncrement = 2;
	  tr->dnaIncrement    = 4;
	  tr->aaIncrement     = 20;
	  tr->secondaryIncrement = 16;
	  tr->secondaryIncrement7 = 7;
	  tr->secondaryIncrement6 = 6;
	}
      else
	{
	  /* in all other cases a switch to GAMMA will take place at some point
	     of time */
	  tr->binaryIncrement = 8;
	  tr->dnaIncrement    = 16;
	  tr->aaIncrement     = 80;
	  tr->secondaryIncrement = 64;
	  tr->secondaryIncrement6 = 24;
	  tr->secondaryIncrement7 = 28;
	}
      break;
    default:
      assert(0);
    }  
}


double gettime(void)
{
#ifdef WIN32
  time_t tp;
  struct tm localtm;
  tp = time(NULL);
  localtm = *localtime(&tp);
  return 60.0*localtm.tm_min + localtm.tm_sec;
#else
  struct timeval ttime;
  gettimeofday(&ttime , NULL);
  return ttime.tv_sec + ttime.tv_usec * 0.000001;
#endif
}

int gettimeSrand(void)
{
#ifdef WIN32
  time_t tp;
  struct tm localtm;
  tp = time(NULL);
  localtm = *localtime(&tp);
  return 24*60*60*localtm.tm_yday + 60*60*localtm.tm_hour + 60*localtm.tm_min  + localtm.tm_sec;
#else
  struct timeval ttime;
  gettimeofday(&ttime , NULL);
  return ttime.tv_sec + ttime.tv_usec;
#endif
}

double randum (long  *seed)    
{
  long  sum, mult0, mult1, seed0, seed1, seed2, newseed0, newseed1, newseed2;
  double res;
  
  mult0 = 1549;
  seed0 = *seed & 4095;
  sum  = mult0 * seed0;
  newseed0 = sum & 4095;
  sum >>= 12;
  seed1 = (*seed >> 12) & 4095;
  mult1 =  406;
  sum += mult0 * seed1 + mult1 * seed0;
  newseed1 = sum & 4095;
  sum >>= 12;
  seed2 = (*seed >> 24) & 255;
  sum += mult0 * seed2 + mult1 * seed1;
  newseed2 = sum & 255;
  
  *seed = newseed2 << 24 | newseed1 << 12 | newseed0;
  res = 0.00390625 * (newseed2 + 0.000244140625 * (newseed1 + 0.000244140625 * newseed0));     
  
  return res;   
}

static int filexists(char *filename)
{
  FILE *fp;
  int res;
  fp = fopen(filename,"r");
  
  if(fp) 
    {
      res = 1;
      fclose(fp);
    }
  else 
    res = 0;
       
  return res;
} 


FILE *myfopen(const char *path, const char *mode)
{
  FILE *fp = fopen(path, mode);

  if(strcmp(mode,"r") == 0 || strcmp(mode,"rb") == 0)
    {      
      if(fp)
	return fp;
      else
	{
	  if(processID == 0)
	    printf("The file %s you want to open for reading does not exist, exiting ...\n", path);
	  errorExit(-1);
	  return (FILE *)NULL;
	}
    }
  else
    {
      if(fp)
	return fp;
      else
	{
	  if(processID == 0)
	    printf("The file %s RAxML wants to open for writing or appending can not be opened [mode: %s], exiting ...\n",
		   path, mode);
	  errorExit(-1);
	  return (FILE *)NULL;
	}
    }

  
}


int countTrees(FILE *f)
{
  int numberOfTrees = 0, ch;

  while((ch = getc(f)) != EOF)
    {
      if(ch == ';')
	numberOfTrees++;
    }	
 
  rewind(f);
  
  return numberOfTrees;
}


/********************* END UTILITY FUNCTIONS ********************/


/******************************some functions for the likelihood computation ****************************/

#ifdef WIN32
boolean isTip(int number, int maxTips)
#else
inline boolean isTip(int number, int maxTips)
#endif
{  
  assert(number > 0);

  if(number <= maxTips)
    return TRUE;
  else
    return FALSE;
}






#ifdef _MULTI_GENE
void getxsnode (nodeptr p, int model)  
{  
  assert(p->xs[model] || p->next->xs[model] || p->next->next->xs[model]);
  assert(p->xs[model] + p->next->xs[model] + p->next->next->xs[model] == 1);
  
  assert(p == p->next->next->next);

  p->xs[model] = 1;
  
  if(p->next->xs[model])
    {      
      p->next->xs[model] = 0;
      return;
    }
  else
    {
      p->next->next->xs[model] = 0;
      return;
    }  

  assert(0);
}

#endif

void getxnode (nodeptr p)  
{ 
  nodeptr  s;
 
  if ((s = p->next)->x || (s = s->next)->x) 
    {
      p->x = s->x;
      s->x = 0;
    }     
  
  assert(p->x);
}





void hookup (nodeptr p, nodeptr q, double *z, int numBranches)
{
  int i;

  p->back = q;
  q->back = p;
    
  for(i = 0; i < numBranches; i++)
    p->z[i] = q->z[i] = z[i];
}

void hookupDefault (nodeptr p, nodeptr q, int numBranches)
{
  int i;
  
  p->back = q;
  q->back = p;
    
  for(i = 0; i < numBranches; i++)
    p->z[i] = q->z[i] = defaultz;  
}


/******************************some functions for the likelihood computation ****************************/




/***********************reading and initializing input ******************/

static void getnums (rawdata *rdta)
{    
  if (fscanf(INFILE, "%d %d", & rdta->numsp, & rdta->sites) != 2) 
    {
      if(processID == 0)
	printf("ERROR: Problem reading number of species and sites\n");
      errorExit(-1);
    }
    
  if (rdta->numsp < 4) 
    {
      if(processID == 0)
	printf("TOO FEW SPECIES\n");
      errorExit(-1);
    }

  if (rdta->sites < 1) 
    {
      if(processID == 0)
	printf("TOO FEW SITES\n");
      errorExit(-1);
    }
  
  return;
}





boolean whitechar (int ch)
{ 
  return (ch == ' ' || ch == '\n' || ch == '\t' || ch == '\r'); 
}


static void uppercase (int *chptr)
{
  int  ch;
  
  ch = *chptr;
  if ((ch >= 'a' && ch <= 'i') || (ch >= 'j' && ch <= 'r')
      || (ch >= 's' && ch <= 'z'))
    *chptr = ch + 'A' - 'a';
} 




static void getyspace (rawdata *rdta)
{
  long   size;
  int    i;
  unsigned char *y0;

  if (! (rdta->y = (unsigned char **) malloc((rdta->numsp + 1) * sizeof(unsigned char *)))) 
    {
      printf("ERROR: Unable to obtain space for data array pointers\n");
      exit(-1);
    }

  size = 4 * (rdta->sites / 4 + 1);

  if (! (y0 = (unsigned char *) malloc((rdta->numsp + 1) * size * sizeof(unsigned char)))) 
    {
      printf("ERROR: Unable to obtain space for data array\n");
      exit(-1);
    }

  rdta->y0 = y0;

  for (i = 0; i <= rdta->numsp; i++) 
    {
      rdta->y[i] = y0;
      y0 += size;
    }

  return;
}

unsigned int x=123456789,y=362436069,z=21288629,w=14921776,c=0;
static unsigned int KISS32(void)
{
  unsigned int t;
  x += 545925293;
  y ^= (y<<13); y ^= (y>>17); y ^= (y<<5);
  t = z+w+c; z = w; c = (t>>31); w = t&2147483647;
  return(x+y+w);
}

static boolean setupTree (tree *tr, analdef *adef)
{
  nodeptr  p0, p, q;
  int      
    i, 
    j, 
    tips, 
    inter;

  tr->bigCutoff = FALSE;

 
  
  tr->partitionContributions = (double *)malloc(sizeof(double) * tr->NumberOfModels);

  for(i = 0; i < tr->NumberOfModels; i++)
    tr->partitionContributions[i] = -1.0;    

  tr->perPartitionLH = (double *)malloc(sizeof(double) * tr->NumberOfModels);

  for(i = 0; i < tr->NumberOfModels; i++)
    tr->perPartitionLH[i] = 0.0;
  
  if(adef->grouping)
    tr->grouped = TRUE;
  else
    tr->grouped = FALSE;
  
  if(adef->constraint)
    tr->constrained = TRUE;
  else
    tr->constrained = FALSE;
  
  tr->treeID = 0;
  
  tips  = tr->mxtips;
  inter = tr->mxtips - 1;
  
 
  
 
  tr->yVector      = (unsigned char **)  malloc((tr->mxtips + 1) * sizeof(unsigned char *));

  tr->fracchanges  = (double *)malloc(tr->NumberOfModels * sizeof(double));
  tr->likelihoods  = (double *)malloc(adef->multipleRuns * sizeof(double)); 
 
 

  tr->treeStringLength = tr->mxtips * (nmlngth+128) + 256 + tr->mxtips * 2;
  
  tr->tree_string  = (char*)malloc(tr->treeStringLength * sizeof(char));

  /*TODO, must that be so long ?*/  
 

  tr->td[0].count = 0;
  tr->td[0].ti    = (traversalInfo *)malloc(sizeof(traversalInfo) * tr->mxtips);

  for(i = 0; i < tr->NumberOfModels; i++)
    tr->fracchanges[i] = -1.0;
  tr->fracchange = -1.0;





  tr->constraintVector = (int *)malloc((2 * tr->mxtips) * sizeof(int));

  tr->nameList = (char **)malloc(sizeof(char *) * (tips + 1));    
             
  if (!(p0 = (nodeptr) malloc((tips + 3*inter) * sizeof(node)))) 
    {
      printf("ERROR: Unable to obtain sufficient tree memory\n");
      return  FALSE;
    }

  if (!(tr->nodep = (nodeptr *) malloc((2*tr->mxtips) * sizeof(nodeptr)))) 
    {
      printf("ERROR: Unable to obtain sufficient tree memory, too\n");
      return  FALSE;
    }
    
  tr->nodep[0] = (node *) NULL;    /* Use as 1-based array */

  for (i = 1; i <= tips; i++) 
    {
      p = p0++;
     
      p->hash   =  KISS32(); /* hast table stuff */
      p->x      =  0;    
      p->number =  i;
      p->next   =  p;
      p->back   = (node *)NULL; 
      p->bInf   = (branchInfo *)NULL;
#ifdef  _MULTI_GENE
      {
	int k;
	
	for(k = 0; k < tr->numBranches; k++)
	  {
	    p->xs[k]    = 0;
	    p->backs[k] = (nodeptr)NULL;
	  }
      }
#endif
      tr->nodep[i] = p;
    }

  for (i = tips + 1; i <= tips + inter; i++) 
    {
      q = (node *) NULL;
      for (j = 1; j <= 3; j++) 
	{
	  p = p0++;
	  p->x      =  0;	  
	  p->number = i;
	  p->next   = q;
	  p->bInf   = (branchInfo *)NULL;
	  p->back   = (node *) NULL;
#ifdef  _MULTI_GENE
	  {
	    int k;
		  
	    for(k = 0; k < tr->numBranches; k++)
	      {
		p->xs[k]    = 0;
		p->backs[k] = (nodeptr)NULL;
	      }
	  }
#endif
	  q = p;
	}
      p->next->next->next = p;
      tr->nodep[i] = p;
    }

  tr->likelihood  = unlikely;
  tr->start       = (node *) NULL;
  tr->ntips       = 0;
  tr->nextnode    = 0;   
  
  for(i = 0; i < tr->numBranches; i++)
    tr->partitionSmoothed[i] = FALSE;
  
  return TRUE;
} 


static void checkTaxonName(char *buffer, int len)
{
  int i;

  for(i = 0; i < len - 1; i++)
    {
      boolean valid;

      switch(buffer[i])
	{
	case '\0':  
	case '\t':  
	case '\n':  
	case '\r': 
	case ' ':
	case ':':  
	case ',':   
	case '(':   
	case ')':  
	case ';':
	case '[':
	case ']':
	  valid = FALSE;
	  break;
	default:
	  valid = TRUE;
	}
     
      if(!valid)
	{
	  printf("ERROR: Taxon Name \"%s\" is invalid at position %d, it contains illegal character %c\n", buffer, i, buffer[i]);
	  printf("Illegal characters in taxon-names are: tabulators, carriage returns, spaces, \":\", \",\", \")\", \"(\", \";\", \"]\", \"[\"\n");
	  printf("Exiting\n");
	  exit(-1);
	}

    }
  assert(buffer[len - 1] == '\0');
}

static boolean getdata(analdef *adef, rawdata *rdta, tree *tr)
{
  int   i, j, basesread, basesnew, ch, my_i, meaning;
  int   meaningAA[256], meaningDNA[256], meaningBINARY[256];
  boolean  allread, firstpass;
  char buffer[nmlngth + 2]; 
  int len;
  unsigned long total = 0;
  unsigned long gaps  = 0;  
     
  for (i = 0; i < 256; i++) 
    {
      meaningAA[i]          = -1;
      meaningDNA[i]         = -1;
      meaningBINARY[i]      = -1;      
    }

  /* AA data */

  meaningAA['A'] =  0;  /* alanine */
  meaningAA['R'] =  1;  /* arginine */
  meaningAA['N'] =  2;  /*  asparagine*/
  meaningAA['D'] =  3;  /* aspartic */
  meaningAA['C'] =  4;  /* cysteine */
  meaningAA['Q'] =  5;  /* glutamine */
  meaningAA['E'] =  6;  /* glutamic */
  meaningAA['G'] =  7;  /* glycine */
  meaningAA['H'] =  8;  /* histidine */
  meaningAA['I'] =  9;  /* isoleucine */
  meaningAA['L'] =  10; /* leucine */
  meaningAA['K'] =  11; /* lysine */
  meaningAA['M'] =  12; /* methionine */
  meaningAA['F'] =  13; /* phenylalanine */
  meaningAA['P'] =  14; /* proline */
  meaningAA['S'] =  15; /* serine */
  meaningAA['T'] =  16; /* threonine */
  meaningAA['W'] =  17; /* tryptophan */
  meaningAA['Y'] =  18; /* tyrosine */
  meaningAA['V'] =  19; /* valine */
  meaningAA['B'] =  20; /* asparagine, aspartic 2 and 3*/
  meaningAA['Z'] =  21; /*21 glutamine glutamic 5 and 6*/
  meaningAA['X'] =  22;
  meaningAA['?'] =  22;
  meaningAA['*'] =  22;
  meaningAA['-'] =  22; /* all = 1.0 */  

  /* DNA data */

  meaningDNA['A'] =  1;
  meaningDNA['B'] = 14;
  meaningDNA['C'] =  2;
  meaningDNA['D'] = 13;
  meaningDNA['G'] =  4;
  meaningDNA['H'] = 11;
  meaningDNA['K'] = 12;
  meaningDNA['M'] =  3;
  meaningDNA['N'] = 15;
  meaningDNA['O'] = 15;
  meaningDNA['R'] =  5;
  meaningDNA['S'] =  6;
  meaningDNA['T'] =  8;
  meaningDNA['U'] =  8;
  meaningDNA['V'] =  7;
  meaningDNA['W'] =  9;
  meaningDNA['X'] = 15;
  meaningDNA['Y'] = 10;     
  meaningDNA['-'] = 15;	
  meaningDNA['?'] = 15; 

  /* BINARY DATA */

  meaningBINARY['0'] = 1;
  meaningBINARY['1'] = 2;
  meaningBINARY['-'] = 3;
  meaningBINARY['?'] = 3; 


  /*******************************************************************/
  
  basesread = basesnew = 0;

  allread = FALSE;
  firstpass = TRUE;
  ch = ' ';

  while (! allread) 
    {
      for (i = 1; i <= tr->mxtips; i++) 
	{   	  
	  if (firstpass) 
	    {                      	       
	      ch = getc(INFILE);
	      while(ch == ' ' || ch == '\n' || ch == '\t' || ch == '\r') 		
		ch = getc(INFILE);		  
			      
	      my_i = 0;	      

	      do 
		{
		  buffer[my_i] = ch;		  
		  ch = getc(INFILE);		   
		  my_i++;
		  if(my_i >= nmlngth)
		    {
		      if(processID == 0)
			{
			  printf("Taxon Name to long at taxon %d, adapt constant nmlngth in\n", i);
			  printf("axml.h, current setting %d\n", nmlngth);
			}
		      errorExit(-1);
		    }		 
		}
	      while(ch !=  ' ' && ch != '\n' && ch != '\t' && ch != '\r');
			      
	      buffer[my_i] = '\0';	      
	      len = strlen(buffer) + 1;	
	      checkTaxonName(buffer, len);
	      tr->nameList[i] = (char *)malloc(sizeof(char) * len);	      
	      strcpy(tr->nameList[i], buffer);				
	    }

	  j = basesread;

	  while ((j < rdta->sites) && ((ch = getc(INFILE)) != EOF) && (ch != '\n') && (ch != '\r'))
	    {	    
	      uppercase(& ch);
	     
	      assert(tr->dataVector[j + 1] != -1);

	      switch(tr->dataVector[j + 1])
		{
		case BINARY_DATA:		 
		  meaning = meaningBINARY[ch];		  
		  break;
		case DNA_DATA:
		  meaning = meaningDNA[ch];
		  break;
		case AA_DATA:
		  meaning = meaningAA[ch];
		  break;
		case SECONDARY_DATA:
		case SECONDARY_DATA_6:
		case SECONDARY_DATA_7:
		  meaning = meaningDNA[ch];
		  break;		 
		default:
		  assert(0);
		}

	      if (meaning != -1) 
		{
		  j++;		 
		  rdta->y[i][j] = ch;
		}
	      else 
		{		
		  if(!whitechar(ch))		  
		    {
		      printf("ERROR: Bad base (%c) at site %d of sequence %d\n",
			     ch, j + 1, i);
		      return FALSE;
		    }
		}
	    }

	  if (ch == EOF) 
	    {
	      printf("ERROR: End-of-file at site %d of sequence %d\n", j + 1, i);
	      return  FALSE;
	    }

	  if (! firstpass && (j == basesread)) 
	    i--; 
	  else 
	    {
	      if (i == 1) 
		basesnew = j;
	      else 
		if (j != basesnew) 
		  {
		    printf("ERROR: Sequences out of alignment\n");		    
		    printf("%d (instead of %d) residues read in sequence %d %s\n",
			   j - basesread, basesnew - basesread, i, tr->nameList[i]);
		    return  FALSE;
		  }
	    }
	  while (ch != '\n' && ch != EOF && ch != '\r') ch = getc(INFILE);  /* flush line *//* PC-LINEBREAK*/
	}

      firstpass = FALSE;
      basesread = basesnew;
      allread = (basesread >= rdta->sites);
    }      

  for(j = 1; j <= tr->mxtips; j++)    
    for(i = 1; i <= rdta->sites; i++) 	
      {
	assert(tr->dataVector[i] != -1);
	
	switch(tr->dataVector[i])
	  {
	  case DNA_DATA:
	    meaning = meaningDNA[rdta->y[j][i]];
	    if(meaning == UNDETERMINED_DNA)
	      gaps++;
	    break;
	  case AA_DATA:
	    meaning = meaningAA[rdta->y[j][i]];
	    if(meaning == UNDETERMINED_AA)
	      gaps++;
	    break;
	  case BINARY_DATA:
	    meaning = meaningBINARY[rdta->y[j][i]];
	    if(meaning == UNDETERMINED_BINARY)
	      gaps++;
	    break;
	  case SECONDARY_DATA:
	  case SECONDARY_DATA_6:
	  case SECONDARY_DATA_7:
	    meaning = meaningDNA[rdta->y[j][i]];
	    assert(tr->secondaryStructurePairs[i - 1] != -1);
	    assert(i - 1 == tr->secondaryStructurePairs[tr->secondaryStructurePairs[i - 1]]);
	    if(meaningDNA[rdta->y[j][tr->secondaryStructurePairs[i - 1] + 1]] == UNDETERMINED_DNA &&
	       meaningDNA[rdta->y[j][tr->secondaryStructurePairs[tr->secondaryStructurePairs[i - 1]] + 1]] ==  UNDETERMINED_DNA &&
	       tr->secondaryStructurePairs[i - 1] > tr->secondaryStructurePairs[tr->secondaryStructurePairs[i - 1]])
	      gaps++;
	       
	    /*TODO 
	      if(meaning == gapValueSECONDARY)
	      gaps++;
	    */
	    break;	    

	  default:
	    assert(0);
	  }     
	
	total++;
	rdta->y[j][i] = meaning;	    
      }
  
  adef->gapyness = (double)gaps / (double)total;

  return  TRUE;
}



static void inputweights (rawdata *rdta)    
{
  int i, w, fres;
  FILE *weightFile;
  int *wv = (int *)malloc(sizeof(int) *  rdta->sites + 1);
    
  weightFile = myfopen(weightFileName, "r");
  
  i = 1;
  
  while((fres = fscanf(weightFile,"%d", &w)) != EOF)
    {
      if(!fres)
	{
	  if(processID == 0)
	    printf("error reading weight file probably encountered a non-integer weight value\n");
	  errorExit(-1);
	}
      wv[i] = w;
      i++;	
    }
    
  if(i != (rdta->sites + 1))
    {
      if(processID == 0)
	printf("number %d of weights not equal to number %d of alignment columns\n", i, rdta->sites);
      errorExit(-1);
    }
 
  for(i = 1; i <= rdta->sites; i++)     
    rdta->wgt[i] = wv[i];         
  
  fclose(weightFile);
  free(wv);
}



static void getinput(analdef *adef, rawdata *rdta, cruncheddata *cdta, tree *tr)
{ 
  int i;  

  getnums(rdta);

  tr->mxtips         = rdta->numsp;
  rdta->wgt          = (int *)    malloc((rdta->sites + 1) * sizeof(int));
  rdta->wgt2         = (int *)    malloc((rdta->sites + 1) * sizeof(int));
  cdta->alias        = (int *)    malloc((rdta->sites + 1) * sizeof(int)); 
  cdta->aliaswgt     = (int *)    malloc((rdta->sites + 1) * sizeof(int));      
  cdta->rateCategory = (int *)    malloc((rdta->sites + 1) * sizeof(int)); 
  tr->model          = (int *)    calloc((rdta->sites + 1), sizeof(int));
  tr->dataVector     = (int *)    malloc((rdta->sites + 1) * sizeof(int));
  cdta->wr           = (double *) malloc((rdta->sites + 1) * sizeof(double));  
  cdta->wr2          = (double *) malloc((rdta->sites + 1) * sizeof(double));  	
  cdta->patrat       = (double *) malloc((rdta->sites + 1) * sizeof(double));
  cdta->patratStored = (double *) malloc((rdta->sites + 1) * sizeof(double));             

      
  if(!adef->useWeightFile)
    {
      for (i = 1; i <= rdta->sites; i++) 
	rdta->wgt[i] = 1;         
    }
  else   
    {
      assert(!adef->useSecondaryStructure);
      inputweights(rdta);   
    }
  
  tr->multiBranch = 0;
  tr->numBranches = 1;

  if(adef->useMultipleModel)  
    {
      int ref;
     
      parsePartitions(adef, rdta, tr);     
     
      for(i = 1; i <= rdta->sites; i++)
	{
	  ref = tr->model[i];
	  tr->dataVector[i] = tr->partitionData[ref].dataType;
	}                               
    }
  else
    {     
      int dataType = -1;

     
      tr->partitionData  = (pInfo*)malloc(sizeof(pInfo));
      tr->partitionData[0].partitionName = (char*)malloc(128 * sizeof(char));    
      strcpy(tr->partitionData[0].partitionName, "No Name Provided");

      tr->partitionData[0].protModels = adef->proteinMatrix;
      tr->partitionData[0].protFreqs  = adef->protEmpiricalFreqs;


      tr->NumberOfModels = 1;     

      if(adef->model == M_PROTCAT || adef->model == M_PROTGAMMA)
	dataType = AA_DATA;             
      if(adef->model == M_GTRCAT || adef->model == M_GTRGAMMA)
	dataType = DNA_DATA;
      if(adef->model == M_BINCAT || adef->model == M_BINGAMMA)
	dataType = BINARY_DATA;      
      assert(dataType >= 0);
      
      /* INIT data-type, model, dataVector for good */
      /* those values will be constant throughout the */
      /* inference process */
      
      tr->partitionData[0].dataType = dataType;     
      
      for(i = 0; i <= rdta->sites; i++)
	{
	  tr->dataVector[i] = dataType;      
	  tr->model[i]      = 0;
	}
    }  

  /*printPartitions(tr);*/

  parseSecondaryStructure(tr, adef, rdta->sites);

  /*printPartitions(tr);*/

  tr->executeModel   = (boolean *)malloc(sizeof(boolean) * tr->NumberOfModels);

  for(i = 0; i < tr->NumberOfModels; i++)
    tr->executeModel[i] = TRUE;

#ifdef _MULTI_GENE
  {
    int i;

    tr->startVector = (nodeptr *)malloc(sizeof(nodeptr) * tr->NumberOfModels);
    tr->tipMissing = (char **)malloc(sizeof(char *) * (tr->mxtips + 1));
    
    for(i = 0; i <= tr->mxtips; i++)
      tr->tipMissing[i] = (char *)malloc(sizeof(char) * (tr->NumberOfModels));
  }				     
#endif
 

  getyspace(rdta);
  setupTree(tr, adef);
    
  if(!getdata(adef, rdta, tr))
    {
      printf("Problem reading alignment file \n");
      errorExit(1);
    }
        
  return;
} 



static unsigned char buildStates(int secModel, unsigned char v1, unsigned char v2)
{ 
  unsigned char new = 0;

  switch(secModel)
    {
    case SECONDARY_DATA:
      new = v1;
      new = new << 4;
      new = new | v2;
      break;
    case SECONDARY_DATA_6:
      { 
	int 
	  meaningDNA[256], 
	  i;

	const unsigned char
	  allowedStates[6][2] = {{'A','T'}, {'C', 'G'}, {'G', 'C'}, {'G','T'}, {'T', 'A'}, {'T', 'G'}};
	  /* use same ordering as in phase */
	  /* allowedStates[6][2] = {{'A','T'}, {'G', 'T'}, {'G', 'C'}, {'T','A'}, {'T', 'G'}, {'C', 'G'}}; */
	
	const unsigned char
	  finalBinaryStates[6] = {1, 2, 4, 8, 16, 32};

	unsigned char	
	  intermediateBinaryStates[6];          
	  
	int length = 6;

	for(i = 0; i < 256; i++) 
	  meaningDNA[i] = -1;

	meaningDNA['A'] =  1;
	meaningDNA['B'] = 14;
	meaningDNA['C'] =  2;
	meaningDNA['D'] = 13;
	meaningDNA['G'] =  4;
	meaningDNA['H'] = 11;
	meaningDNA['K'] = 12;
	meaningDNA['M'] =  3;
	meaningDNA['N'] = 15;
	meaningDNA['O'] = 15;
	meaningDNA['R'] =  5;
	meaningDNA['S'] =  6;
	meaningDNA['T'] =  8;
	meaningDNA['U'] =  8;
	meaningDNA['V'] =  7;
	meaningDNA['W'] =  9;
	meaningDNA['X'] = 15;
	meaningDNA['Y'] = 10;     
	meaningDNA['-'] = 15;	
	meaningDNA['?'] = 15; 
      
	/*printf("Character: %c%c\n", inverseMeaningDNA[v1],  inverseMeaningDNA[v2]);*/

	for(i = 0; i < length; i++)
	  {
	    unsigned char n1 = meaningDNA[allowedStates[i][0]];
	    unsigned char n2 = meaningDNA[allowedStates[i][1]];
	    
	    new = n1;
	    new = new << 4;
	    new = new | n2;

	    intermediateBinaryStates[i] = new;
	  }	

	new = v1;
	new = new << 4;
	new = new | v2;
      
	for(i = 0; i < length; i++)
	  {
	    if(new == intermediateBinaryStates[i])
	      break;
	  }
	if(i < length)
	  new = finalBinaryStates[i];
	else
	  {	    
	    new = 0;	
	    for(i = 0; i < length; i++)
	      {
		if(v1 & meaningDNA[allowedStates[i][0]])
		  {
		    /*printf("Adding %c%c\n", allowedStates[i][0], allowedStates[i][1]);*/
		    new |= finalBinaryStates[i];		 
		  }
		if(v2 & meaningDNA[allowedStates[i][1]])
		  {
		    /*printf("Adding %c%c\n", allowedStates[i][0], allowedStates[i][1]);*/
		    new |= finalBinaryStates[i];		  
		  }
	      }
	  }
	/*printf("Final %c%c %d\n\n", inverseMeaningDNA[v1], inverseMeaningDNA[v2], new); 	 	  */
      }
      break;
    case SECONDARY_DATA_7:
      { 
	int 
	  meaningDNA[256], 
	  i;

	const unsigned char
	  allowedStates[6][2] = {{'A','T'}, {'C', 'G'}, {'G', 'C'}, {'G','T'}, {'T', 'A'}, {'T', 'G'}};
	  /*                           1           2          4           8          16          32      */
	  /* use same ordering as in phase */
	  /* allowedStates[6][2] = {{'A','T'}, {'G', 'T'}, {'G', 'C'}, {'T','A'}, {'T', 'G'}, {'C', 'G'}}; */
	const unsigned char
	  finalBinaryStates[7] = {1, 2, 4, 8, 16, 32, 64};
	/*finalBinaryStates[7] = {1, 8, 4, 16, 32, 2, 64};*/

	unsigned char	
	  intermediateBinaryStates[7];          	        

	for(i = 0; i < 256; i++) 
	  meaningDNA[i] = -1;

	meaningDNA['A'] =  1;
	meaningDNA['B'] = 14;
	meaningDNA['C'] =  2;
	meaningDNA['D'] = 13;
	meaningDNA['G'] =  4;
	meaningDNA['H'] = 11;
	meaningDNA['K'] = 12;
	meaningDNA['M'] =  3;
	meaningDNA['N'] = 15;
	meaningDNA['O'] = 15;
	meaningDNA['R'] =  5;
	meaningDNA['S'] =  6;
	meaningDNA['T'] =  8;
	meaningDNA['U'] =  8;
	meaningDNA['V'] =  7;
	meaningDNA['W'] =  9;
	meaningDNA['X'] = 15;
	meaningDNA['Y'] = 10;     
	meaningDNA['-'] = 15;	
	meaningDNA['?'] = 15; 
      
	/*printf("Character: %c%c\n", inverseMeaningDNA[v1],  inverseMeaningDNA[v2]);*/

	for(i = 0; i < 6; i++)
	  {
	    unsigned char n1 = meaningDNA[allowedStates[i][0]];
	    unsigned char n2 = meaningDNA[allowedStates[i][1]];
	    
	    new = n1;
	    new = new << 4;
	    new = new | n2;

	    intermediateBinaryStates[i] = new;
	  }	

	new = v1;
	new = new << 4;
	new = new | v2;
      
	for(i = 0; i < 6; i++)
	  {
	    /* exact match */
	    if(new == intermediateBinaryStates[i])
	      break;
	  }
	if(i < 6)
	  new = finalBinaryStates[i];
	else
	  {	
	    /* distinguish between exact mismatches and partial mismatches */	  

	    for(i = 0; i < 6; i++)
	      if((v1 & meaningDNA[allowedStates[i][0]]) && (v2 & meaningDNA[allowedStates[i][1]]))
		break;
	    if(i < 6)
	      {
		/* printf("partial mismatch\n"); */
		
		new = 0;
		for(i = 0; i < 6; i++)
		  {		    
		    if((v1 & meaningDNA[allowedStates[i][0]]) && (v2 & meaningDNA[allowedStates[i][1]]))
		      {
			/*printf("Adding %c%c\n", allowedStates[i][0], allowedStates[i][1]);*/
			new |= finalBinaryStates[i];		 
		      }
		    else
		      new |=  finalBinaryStates[6];		   
		  }	       
	      }
	    else
	      new = finalBinaryStates[6];	   
	  }
	/* printf("Final %c%c %d\n\n", inverseMeaningDNA[v1], inverseMeaningDNA[v2], new); 	 	   */
      }
      break;
    default:
      assert(0);
    }

  return new;  

}



static void adaptRdataToSecondary(tree *tr, rawdata *rdta)
{
  int *alias = (int*)calloc(rdta->sites, sizeof(int));
  int i, j, realPosition;
  
  /*tr->secondaryStructureModel = SEC_RNA_6;*/

  /*
    unsigned char v1, v2;

    unsigned char states[16] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
    unsigned int allStates[256];
    
    for(i = 0; i <256; i++)
    allStates[i] = 0;
    
    for(i = 0; i < 16; i++)
    for(j = 0; j < 16; j++)
    {  
    if(i > 0 && j > 0)
    {
    unsigned char new;
    unsigned int reassign;
    v1 = states[i];
    v2 = states[j];
    
    new = v1;
    new = new << 4;
    new = new | v2;
    
    {
    int l;
    reassign = 0;
    
    for(l = 0; l < 4; l++)
    if((v1 >> l) & 1)
    {
    unsigned int value = (unsigned int)v2;
    value = value << (l * 4);
    reassign |= value;
    }
    }
    allStates[new] = reassign;
    }
    
    
    
    printf("%c%c %d %d\n", inverseMeaningDNA[states[i]], inverseMeaningDNA[states[j]], new, reassign);
    } 
    
    for(i = 0; i < 256; i++)
    printf("%d, ", allStates[i]);
    printf("\n");
  */

  for(i = 0; i < rdta->sites; i++)
    alias[i] = -1;

  for(i = 0, realPosition = 0; i < rdta->sites; i++)
    {
      int partner = tr->secondaryStructurePairs[i];     
      if(partner != -1)
	{
	  assert(tr->dataVector[i+1] == SECONDARY_DATA || tr->dataVector[i+1] == SECONDARY_DATA_6 || tr->dataVector[i+1] == SECONDARY_DATA_7);
	  
	  if(i < partner)
	    {	 
	      for(j = 1; j <= rdta->numsp; j++)
		{
		  unsigned char v1 = rdta->y[j][i+1];
		  unsigned char v2 = rdta->y[j][partner+1];
		  /*unsigned char new = v1;*/
		  assert(i+1 < partner+1);
		  
		  /*
		    new = new << 4;
		    new = new | v2;
		    printf("Old states: %c %c, new state %d\n", inverseMeaningDNA[v1], inverseMeaningDNA[v2], new);
		    rdta->y[j][i+1] = new;
		  */
		  
		  rdta->y[j][i+1] = buildStates(tr->dataVector[i+1], v1, v2);
		}
	      alias[realPosition] = i;
	      realPosition++;
	    }       	  
	}
      else
	{
	  alias[realPosition] = i;	
	  realPosition++;
	}
    }  
    
  assert(rdta->sites - realPosition == tr->numberOfSecondaryColumns / 2);

  rdta->sites = realPosition;

  for(i = 0; i < rdta->sites; i++)
    {
      assert(alias[i] != -1);
      tr->model[i+1]    = tr->model[alias[i]+1];
      tr->dataVector[i+1] = tr->dataVector[alias[i]+1];     
      rdta->wgt2[i+1] =  rdta->wgt2[alias[i]+1];

      for(j = 1; j <= rdta->numsp; j++)
	rdta->y[j][i+1] = rdta->y[j][alias[i]+1];

      /*printf("%d: %d %d %d\n",  i+1, tr->model[i+1], tr->dataVector[i+1], rdta->wgt2[i+1]);*/
    }
  
  /*printf("SIERRA\n");

  for(j = 1; j <= rdta->numsp; j++)
    printf("%d ", 	rdta->y[j][1]);
  printf("\n");

  for(j = 1; j <= rdta->numsp; j++)
    printf("%d ", 	rdta->y[j][2]);
    printf("\n");*/

  free(alias);
}

static void sitesort(rawdata *rdta, cruncheddata *cdta, tree *tr, analdef *adef)   
{ 
  int  gap, i, j, jj, jg, k, n, nsp;
  int  *index, *category, *superCategory;
  boolean  flip, tied;
  unsigned char  **data;
    
  /*for(i = 0; i <= rdta->sites; i++)
    {
      printf("%d %d %d\n", i, tr->dataVector[i], tr->model[i]);
      }  */
  
  if(adef->useSecondaryStructure)
    {
      assert(tr->NumberOfModels > 1 && adef->useMultipleModel);
      
      adaptRdataToSecondary(tr, rdta);
    }

  if(adef->useMultipleModel)
    {
      superCategory = tr->dataVector;
      category      = tr->model;
    }
  else
    {
      category      = (int*)NULL;
      superCategory = (int*)NULL;
    }

  index    = cdta->alias;
  data     = rdta->y;
  n        = rdta->sites;
  nsp      = rdta->numsp;
  index[0] = -1;

  
  if(adef->compressPatterns)
    {
      for (gap = n / 2; gap > 0; gap /= 2) 
	{
	  for (i = gap + 1; i <= n; i++) 
	    {
	      j = i - gap;
	      
	      do 
		{
		  jj = index[j];
		  jg = index[j+gap];
		  if(adef->useMultipleModel)
		    {	
		      assert(superCategory[jj] != -1 &&
			     superCategory[jg] != -1 &&
			     category[jj] != -1      &&
			     category[jg] != -1);
		      
		      if(superCategory[jj] > superCategory[jg])
			{
			  flip = TRUE;
			  tied = FALSE;
			}
		      else
			{
			  flip = (category[jj] >  category[jg]);
			  tied = (category[jj] == category[jg]);		    
			}
		    }
		  else
		    {		    
		      flip = 0;
		      tied = 1;
		    }
		  
		  for (k = 1; (k <= nsp) && tied; k++) 
		    {
		      flip = (data[k][jj] >  data[k][jg]);
		      tied = (data[k][jj] == data[k][jg]);
		    }
		  
		  if (flip) 
		    {
		      index[j]     = jg;
		      index[j+gap] = jj;
		      j -= gap;
		    }
		} 
	      while (flip && (j > 0));	      
	    }  
	}
    }  
}


static void sitecombcrunch (rawdata *rdta, cruncheddata *cdta, tree *tr, analdef *adef)    
{ 
  int  i, sitei, j, sitej, k;   
  boolean  tied;
  int *aliasModel; 
  int *aliasSuperModel;
    


  if(adef->useMultipleModel)
    {
      aliasSuperModel = (int*)malloc(sizeof(int) * (rdta->sites + 1));
      aliasModel      = (int*)malloc(sizeof(int) * (rdta->sites + 1));
    }
  else
    {
      aliasModel      = (int*)NULL;
      aliasSuperModel = (int*)NULL;
    }
  
  i = 0;    
  cdta->alias[0]    = cdta->alias[1];
  cdta->aliaswgt[0] = 0;   
    
  for (j = 1; j <= rdta->sites; j++) 
    {
      sitei = cdta->alias[i];
      sitej = cdta->alias[j];
      if(!adef->compressPatterns)
	tied = 0;
      else
	{	 
	  if(adef->useMultipleModel)
	    {
	      tied = (tr->model[sitei] == tr->model[sitej]);	    
	      if(tied)
		assert(tr->dataVector[sitei] == tr->dataVector[sitej]);  
	    }
	  else	      
	    tied = 1;	      
	}
      
      for (k = 1; tied && (k <= rdta->numsp); k++)
	tied = (rdta->y[k][sitei] == rdta->y[k][sitej]);
      
      if (tied) 
	{
	  cdta->aliaswgt[i] += rdta->wgt2[sitej];	   
	  if(adef->useMultipleModel)
	    {
	      aliasModel[i]      = tr->model[sitej];	  	     	  
	      aliasSuperModel[i] = tr->dataVector[sitej];
	    }
	}
      else 
	{
	  if (cdta->aliaswgt[i] > 0) i++;
	  cdta->aliaswgt[i] = rdta->wgt2[sitej];
	  cdta->alias[i] = sitej;
	  if(adef->useMultipleModel)
	    {
	      aliasModel[i]      = tr->model[sitej];		   
	      aliasSuperModel[i] = tr->dataVector[sitej];
	    }
	}
    }

  cdta->endsite = i;
  if (cdta->aliaswgt[i] > 0) cdta->endsite++;       
    
  if(adef->useMultipleModel)
    {       
      for(i = 0; i <= rdta->sites; i++)
	{
	  tr->model[i]      = aliasModel[i];	  
	  tr->dataVector[i] = aliasSuperModel[i];
	}
    }
  
  if(adef->useMultipleModel)
    {
      free(aliasModel);	
      free(aliasSuperModel);
    }

  /*
    for(i = 0; i < cdta->endsite; i++)
    {
      printf("%d %d %d\n",  i, tr->model[i], tr->dataVector[i]);
    }
  */
} 


static boolean makeweights (analdef *adef, rawdata *rdta, cruncheddata *cdta, tree *tr)    
{
  int  i;

 
  for (i = 1; i <= rdta->sites; i++)  
    rdta->wgt2[i] = rdta->wgt[i];     

  for (i = 1; i <= rdta->sites; i++)  
    cdta->alias[i] = i;
        
  sitesort(rdta, cdta, tr, adef);
  sitecombcrunch(rdta, cdta, tr, adef);
      
  return TRUE;
} 




static boolean makevalues(rawdata *rdta, cruncheddata *cdta, tree *tr, analdef *adef)   
{   
  int  i, j, model, fullSites = 0, modelCounter;   

  unsigned char *y    = (unsigned char *)malloc(rdta->numsp * cdta->endsite * sizeof(unsigned char));
  unsigned char *yBUF = (unsigned char *)malloc(rdta->numsp * cdta->endsite * sizeof(unsigned char));

  for (i = 1; i <= rdta->numsp; i++)       
    for (j = 0; j < cdta->endsite; j++) 	  	          
      y[((i - 1) * cdta->endsite) + j] = rdta->y[i][cdta->alias[j]];	                

  free(rdta->y0);
  free(rdta->y);
    
  rdta->y0 = y;
  memcpy(yBUF, y, rdta->numsp * cdta->endsite * sizeof(unsigned char));
  rdta->yBUF = yBUF;
                      
  if(!adef->useMultipleModel)
    tr->NumberOfModels = 1;

  if(adef->useMultipleModel)
    {      
      int *perm          = (int*)malloc(sizeof(int) * tr->NumberOfModels);     
      pInfo *partitionData = (pInfo*)malloc(sizeof(pInfo) * tr->NumberOfModels);

      tr->partitionData[0].lower = 0;
           
      model        = tr->model[0];
      modelCounter = 0;
      perm[modelCounter] = model;
      i            = 1;                 

      while(i <  cdta->endsite)
	{
	  if(tr->model[i] != model)
	    {	     
	      tr->partitionData[modelCounter].upper     = i;
	      tr->partitionData[modelCounter + 1].lower = i;

	      model = tr->model[i];
	      perm[modelCounter + 1] = model;
	      modelCounter++;
	    }
	  i++;
	}
      
     
      tr->partitionData[tr->NumberOfModels - 1].upper = cdta->endsite;
     
      memcpy(partitionData, tr->partitionData, tr->NumberOfModels * sizeof(pInfo));
      
      for(i = 0; i < tr->NumberOfModels; i++)
	{	 	  
	  tr->partitionData[i].dataType   = partitionData[perm[i]].dataType;
	  tr->partitionData[i].protModels = partitionData[perm[i]].protModels;
	  tr->partitionData[i].protFreqs  = partitionData[perm[i]].protFreqs;
	  tr->partitionData[i].width      =  tr->partitionData[i].upper -  tr->partitionData[i].lower;
	}            

      model        = tr->model[0];
      modelCounter = 0;
      tr->model[0] = modelCounter;    
      i            = 1;                 

      while(i < cdta->endsite)
	{
	  if(tr->model[i] != model)
	    {	     
	      model = tr->model[i];	      
	      modelCounter++;
	      tr->model[i] = modelCounter;
	    }
	  else
	    tr->model[i] = modelCounter;
	  i++;
	}            

      free(perm);       
      free(partitionData);
    }
  else
    {                
      tr->partitionData[0].lower = 0;
      tr->partitionData[0].upper = cdta->endsite; 
      tr->partitionData[0].width =  tr->partitionData[0].upper -  tr->partitionData[0].lower;
    }




  tr->rdta       = rdta;
  tr->cdta       = cdta;   
  
  tr->invariant          = (int *)malloc(cdta->endsite * sizeof(int));
  tr->originalDataVector = (int *)malloc(cdta->endsite * sizeof(int));
  tr->originalModel      = (int *)malloc(cdta->endsite * sizeof(int));
  tr->originalWeights    = (int *)malloc(cdta->endsite * sizeof(int));

  memcpy(tr->originalModel, tr->model,            cdta->endsite * sizeof(int)); 
  memcpy(tr->originalDataVector, tr->dataVector,  cdta->endsite * sizeof(int));
  memcpy(tr->originalWeights, tr->cdta->aliaswgt, cdta->endsite * sizeof(int));
  
  tr->originalCrunchedLength = tr->cdta->endsite;
  for(i = 0; i < tr->cdta->endsite; i++)
    fullSites += tr->cdta->aliaswgt[i];

  tr->fullSites = fullSites;
  
  for(i = 0; i < rdta->numsp; i++)
    tr->yVector[i + 1] = &(rdta->y0[tr->originalCrunchedLength * i]);

  return TRUE;
} 







static int sequenceSimilarity(unsigned char *tipJ, unsigned char *tipK, int n)
{
  int i;
  
  for(i = 0; i < n; i++)    
    if(*tipJ++ != *tipK++)	
      return 0;	
      
  return 1;
}

static void checkSequences(tree *tr, rawdata *rdta, analdef *adef)
{
  int n = tr->mxtips + 1; 
  int i, j;
  int *omissionList     = (int *)malloc(n * sizeof(int));
  int *undeterminedList = (int *)malloc((rdta->sites + 1)* sizeof(int));
  int *modelList        = (int *)malloc((rdta->sites + 1)* sizeof(int)); 
  int count = 0;
  int countNameDuplicates = 0;
  int countUndeterminedColumns = 0;
  int countOnlyGaps = 0;
  int modelCounter = 1;
  unsigned char *tipI, *tipJ;
  FILE *f;  


  if(processID == 0)	      
    f = myfopen(infoFileName, "a");
  else
    f = (FILE *)NULL; 

  for(i = 1; i < n; i++)       
    omissionList[i] = 0;              

  for(i = 0; i < rdta->sites + 1; i++)
    undeterminedList[i] = 0;
      
  for(i = 1; i < n; i++)
    {
      for(j = i + 1; j < n; j++)
	if(strcmp(tr->nameList[i], tr->nameList[j]) == 0)
	  {
	    countNameDuplicates++;
	    if(processID == 0)
	      {
		printf("Sequence names of taxon %d and %d are identical, they are both called %s\n", i, j, tr->nameList[i]);
		fprintf(f, "Sequence names of taxon %d and %d are identical, they are both called %s\n", i, j, tr->nameList[i]);
	      }
	  }
    }
	  
  if(countNameDuplicates > 0)
    {
      if(processID == 0)
	{
	  printf("ERROR: Found %d taxa that had equal names in the alignment, exiting...\n", countNameDuplicates);
	  fprintf(f, "ERROR: Found %d taxa that had equal names in the alignment, exiting...\n", countNameDuplicates);
	  fclose(f);
	}
      errorExit(-1);
    }

  for(i = 1; i < n; i++)
    {
      j = 1;
      
      while(j <= rdta->sites)
	{	 	      
	  if(tr->dataVector[j] == DNA_DATA && rdta->y[i][j] !=  UNDETERMINED_DNA)
	    break;
	  if(tr->dataVector[j] == AA_DATA && rdta->y[i][j] !=  UNDETERMINED_AA)
	    break;
	  if(tr->dataVector[j] == BINARY_DATA && rdta->y[i][j] !=  UNDETERMINED_BINARY)
	    break; 
	  if(tr->dataVector[j] == SECONDARY_DATA && rdta->y[i][j] !=  UNDETERMINED_SECONDARY)
	    break;
	  /* TODO SEC */
	  if(tr->dataVector[j] == SECONDARY_DATA_6 && rdta->y[i][j] !=  UNDETERMINED_SECONDARY)
	    break;
	  if(tr->dataVector[j] == SECONDARY_DATA_7 && rdta->y[i][j] !=  UNDETERMINED_SECONDARY)
	     break;
	  j++;
	}

      if(j == (rdta->sites + 1))
	{       
	  if(processID == 0)
	    {
	      printf("ERROR: Sequence %s consists entirely of undetermined values which will be treated as missing data\n",      
		     tr->nameList[i]);
	      fprintf(f, "ERROR: Sequence %s consists entirely of undetermined values which will be treated as missing data\n",      
		      tr->nameList[i]);	      
	    }
	  countOnlyGaps++;
	}
      
    }
  
  if(countOnlyGaps > 0)
    {
      if(processID == 0)
	{
	  printf("ERROR: Found %d sequences that consist entirely of undetermined values, exiting...\n", countOnlyGaps);
	  fprintf(f, "ERROR: Found %d sequences that consist entirely of undetermined values, exiting...\n", countOnlyGaps);
	  fclose(f);
	}
      errorExit(-1);
    }

  for(i = 0; i <= rdta->sites; i++)
    modelList[i] = -1;

  for(i = 1; i <= rdta->sites; i++)
    {    
      j = 1;
     
      while(j < n)
	{
	  if(tr->dataVector[i] == DNA_DATA && rdta->y[j][i] !=  UNDETERMINED_DNA)
	    break;
	  if(tr->dataVector[i] == AA_DATA && rdta->y[j][i] != UNDETERMINED_AA)
	    break;	
	  if(tr->dataVector[i] == BINARY_DATA && rdta->y[j][i] != UNDETERMINED_BINARY)
	    break;
	  if(tr->dataVector[i] == SECONDARY_DATA) /* skip this */
	    break;
	  if(tr->dataVector[i] == SECONDARY_DATA_6) /* skip this */
	    break;
	  if(tr->dataVector[i] == SECONDARY_DATA_7) /* skip this */
	    break;
	  j++;
	}
      
      if(j == n)
	{
	  undeterminedList[i] = 1;
	  if(processID == 0)
	    {
	      printf("IMPORTANT WARNING: Alignment column %d contains only undetermined values which will be treated as missing data\n", 
		     i);
	      fprintf(f, "IMPORTANT WARNING: Alignment column %d contains only undetermined values which will be treated as missing data\n", i);
	    }
	  countUndeterminedColumns++;	  
	}
      else
	{
	  if(adef->useMultipleModel)
	    {
	      modelList[modelCounter] = tr->model[i];
	      modelCounter++;
	    }
	}
    }
  

  for(i = 1; i < n; i++)
    {
      if(omissionList[i] == 0)
	{
	  tipI = &(rdta->y[i][1]);

	  for(j = i + 1; j < n; j++)
	    {
	      if(omissionList[j] == 0)
		{
		  tipJ = &(rdta->y[j][1]);
		  if(sequenceSimilarity(tipI, tipJ, rdta->sites))
		    {
		      if(processID == 0)
			{
			  printf("\n\nIMPORTANT WARNING: Sequences %s and %s are exactly identical\n", tr->nameList[i], tr->nameList[j]);
			  fprintf(f, "\n\nIMPORTANT WARNING: Sequences %s and %s are exactly identical\n", tr->nameList[i], tr->nameList[j]);
			}
		      omissionList[j] = 1;
		      count++;
		    }
		}
	    }
	}
    }

  if(count > 0 || countUndeterminedColumns > 0)
    {
      char noDupFile[2048];
      char noDupModels[2048];
              
      if(count > 0)
	{
	  if(processID == 0)
	    {
	      printf("\n");
	      
	      printf("IMPORTANT WARNING\n");
	      
	      printf("Found %d %s that %s exactly identical to other sequences in the alignment.\n", count, (count == 1)?"sequence":"sequences", (count == 1)?"is":"are");
	      printf("Normally they should be excluded from the analysis.\n\n");
	      
	      fprintf(f, "\n");
	      
	      fprintf(f, "IMPORTANT WARNING\n");
	      
	      fprintf(f, "Found %d %s that %s exactly identical to other sequences in the alignment.\n", count, (count == 1)?"sequence":"sequences", (count == 1)?"is":"are");
	      fprintf(f, "Normally they should be excluded from the analysis.\n\n");
	    }
	}
      
      if(countUndeterminedColumns > 0)
	{
	  if(processID == 0)
	    {
	      printf("\n");
	      
	      printf("IMPORTANT WARNING\n");
	      
	      printf("Found %d %s that %s only undetermined values which will be treated as missing data.\n", 
		     countUndeterminedColumns, (countUndeterminedColumns == 1)?"column":"columns", (countUndeterminedColumns == 1)?"contains":"contain");
	      printf("Normally these columns should be excluded from the analysis.\n\n");
	      
	      fprintf(f, "\n");
	      
	      fprintf(f, "IMPORTANT WARNING\n");
	      
	      fprintf(f, "Found %d %s that %s only undetermined values which will be treated as missing data.\n", 
		      countUndeterminedColumns, (countUndeterminedColumns == 1)?"column":"columns", (countUndeterminedColumns == 1)?"contains":"contain");
	      fprintf(f, "Normally these columns should be excluded from the analysis.\n\n");      	  
	    }
	}

      strcpy(noDupFile, seq_file);
      strcat(noDupFile, ".reduced");

      strcpy(noDupModels, modelFileName);
      strcat(noDupModels, ".reduced");

      if(processID == 0)
	{

	  if(adef->useMultipleModel && !filexists(noDupModels) && countUndeterminedColumns)
	    {      
	      FILE *newFile = myfopen(noDupModels, "w");

	      printf("\nJust in case you might need it, a mixed model file with \n");
	      printf("model assignments for undetermined columns removed is printed to file %s\n",noDupModels);

	      fprintf(f, "\nJust in case you might need it, a mixed model file with \n");
	      fprintf(f, "model assignments for undetermined columns removed is printed to file %s\n",noDupModels);
	      
 
	      for(i = 0; i < tr->NumberOfModels; i++)
		{
		  boolean modelStillExists = FALSE;
		  
		  for(j = 1; (j <= rdta->sites) && (!modelStillExists); j++)
		    {
		      if(modelList[j] == i)
			modelStillExists = TRUE;
		    }

		  if(modelStillExists)
		    {	  
		      int k = 1;
		      int lower, upper;
		      int parts = 0;


		      switch(tr->partitionData[i].dataType)
			{
			case AA_DATA:		      		     
			  {
			    char AAmodel[1024];
			    
			    strcpy(AAmodel, protModels[tr->partitionData[i].protModels]);
			    if(tr->partitionData[i].protFreqs)
			      strcat(AAmodel, "F");		  
			  
			    fprintf(newFile, "%s, ", AAmodel);
			  }
			  break;
			case DNA_DATA:
			  fprintf(newFile, "DNA, ");
			  break;
			case BINARY_DATA:
			  fprintf(newFile, "BIN, ");
			  break;
			default:
			  assert(0);
			}

		      fprintf(newFile, "%s = ", tr->partitionData[i].partitionName);
		      
		      while(k <= rdta->sites)
			{
			  if(modelList[k] == i)
			    {
			      lower = k;
			      while((modelList[k + 1] == i) && (k <= rdta->sites))		      			
				k++;
			      upper = k;
			      
			      if(lower == upper)		  
				{
				  if(parts == 0)
				    fprintf(newFile, "%d", lower);
				  else
				    fprintf(newFile, ",%d", lower);
				}
			      else
				{
				  if(parts == 0)
				    fprintf(newFile, "%d-%d", lower, upper);
				  else
				    fprintf(newFile, ",%d-%d", lower, upper);
				}		  
			      parts++;
			    }
			  k++;
			}
		      fprintf(newFile, "\n");
		    }		  
		}	
	      fclose(newFile);
	    }
	  else
	    {
	      if(adef->useMultipleModel)
		{
		  printf("\nA mixed model file with model assignments for undetermined\n");
		  printf("columns removed has already been printed to  file %s\n",noDupModels);

		  fprintf(f, "\nA mixed model file with model assignments for undetermined\n");
		  fprintf(f, "columns removed has already been printed to  file %s\n",noDupModels);
		}	      
	    }
	     

	  if(!filexists(noDupFile))
	    {
	      FILE *newFile;
	      
	      printf("Just in case you might need it, an alignment file with \n");
	      if(count && !countUndeterminedColumns)
		printf("sequence duplicates removed is printed to file %s\n", noDupFile);
	      if(!count && countUndeterminedColumns)
		printf("undetermined columns removed is printed to file %s\n", noDupFile);
	      if(count && countUndeterminedColumns)
		printf("sequence duplicates and undetermined columns removed is printed to file %s\n", noDupFile);
	      
	      fprintf(f, "Just in case you might need it, an alignment file with \n");
	      if(count && !countUndeterminedColumns)
		fprintf(f, "sequence duplicates removed is printed to file %s\n", noDupFile);
	      if(!count && countUndeterminedColumns)
		fprintf(f, "undetermined columns removed is printed to file %s\n", noDupFile);
	      if(count && countUndeterminedColumns)
		fprintf(f, "sequence duplicates and undetermined columns removed is printed to file %s\n", noDupFile);
	      
	      newFile = myfopen(noDupFile, "w");
	      
	      fprintf(newFile, "%d %d\n", tr->mxtips - count, rdta->sites - countUndeterminedColumns);
	      
	      for(i = 1; i < n; i++)
		{
		  if(!omissionList[i])
		    {
		      fprintf(newFile, "%s ", tr->nameList[i]);
		      tipI =  &(rdta->y[i][1]);

		      for(j = 0; j < rdta->sites; j++)
			{
			  if(undeterminedList[j + 1] == 0)
			    {
			      switch(tr->dataVector[j + 1])
				{
				case AA_DATA:
				  fprintf(newFile, "%c", inverseMeaningPROT[tipI[j]]);
				  break;
				case DNA_DATA:
				  fprintf(newFile, "%c", inverseMeaningDNA[tipI[j]]);
				  break;
				case BINARY_DATA:
				  fprintf(newFile, "%c", inverseMeaningBINARY[tipI[j]]);
				  break;
				default:
				  assert(0);
				}
			    }
			}
		      		      
		      fprintf(newFile, "\n");
		    }
		}
	      
	      fclose(newFile);	    
	    }
	  else
	    {
	      if(count && !countUndeterminedColumns)
		printf("An alignment file with sequence duplicates removed has already\n");
	      if(!count && countUndeterminedColumns)
		printf("An alignment file with undetermined columns removed has already\n");
	      if(count && countUndeterminedColumns)
		printf("An alignment file with undetermined columns and sequence duplicates removed has already\n");
	      
	      printf("been printed to file %s\n",  noDupFile);
	      
	      if(count && !countUndeterminedColumns)
		fprintf(f, "An alignment file with sequence duplicates removed has already\n");
	      if(!count && countUndeterminedColumns)
		fprintf(f, "An alignment file with undetermined columns removed has already\n");
	      if(count && countUndeterminedColumns)
		fprintf(f, "An alignment file with undetermined columns and sequence duplicates removed has already\n");
	      
	      fprintf(f, "been printed to file %s\n",  noDupFile);
	    }
	}
    }


  free(undeterminedList);
  free(omissionList);
  free(modelList);
  if(processID == 0)	      
    fclose(f);
}



static float dist(int i, int j, const int sites, const float nDouble, unsigned char **y)
{
  int k, count;  
  unsigned char *tipI = &(y[i + 1][1]);
  unsigned char *tipJ = &(y[j + 1][1]); 

  for(k = 0, count = 0; k < sites; k ++)
    if(tipI[k] == tipJ[k])
      count++;
  
  return (((float)count) * nDouble);
}

static void distArray(int i, const int sites, const float nDouble, int n, float *ref, int *omitted, unsigned  char **y)
{
  int k, l;  
  unsigned char *tipI = &(y[i + 1][1]); 
  
  for(l = 0; l < n; l++)
    {
      if((!omitted[l]) && (l != i))
	{
	  unsigned char *tipJ = &(y[l + 1][1]); 
	  int count = 0;
	  for(k = 0, count = 0; k < sites; k ++)
	    if(tipI[k] == tipJ[k])
	      count++;
	  ref[l] = (((float)count) * nDouble);
	}
    }
}



static int qtCompare(const void *p1, const void *p2)
{
 qtData *rc1 = (qtData *)p1;
 qtData *rc2 = (qtData *)p2;

  float i = rc1->val;
  float j = rc2->val;
  
  if (i > j)
    return (-1);
  if (i < j)
    return (1);
  return (0);
}


static qtList * clusterQT_LARGE(int n, float thres, int *ccc, rawdata *rdta)
{  
  int clusterCount;
  int i, j;
  int *omitted, *current, *best;  
  qtList *clusters = (qtList*)NULL;
  const float nDouble = 1.0 / (float)(rdta->sites);
  double t = gettime();  
  const int sites = rdta->sites; 
  unsigned char **y = rdta->y; 
  float *ref;
  qtData *candidates;
  
  candidates = (qtData *)malloc(sizeof(qtData) * n);
  clusters = (qtList *)malloc(sizeof(qtList) * n);
  omitted = (int*)calloc(n, sizeof(int));
  current = (int*)malloc(sizeof(int) * n);
  best    = (int*)malloc(sizeof(int) * n);            
  ref     = (float*)malloc(sizeof(float) * n);
  clusterCount = 0;


  for(i = 0; i < n; i++)
    {
      if(!omitted[i])
	{
	  int entCount = 0;	  	     
	  int countCandidates = 0;	  

	  current[entCount++] = i;	      
	  omitted[i] = 1;

	  distArray(i, sites, nDouble, n, ref, omitted, y);       
	 		 
	  for(j = 0; j < n; j++)		
	    {
	      if(!omitted[j] && i != j)
		{
		  float temp;			
		  
		  if((temp = ref[j]) >= thres)
		    {
		      candidates[countCandidates].val    = temp;
		      candidates[countCandidates].number = j;			
		      countCandidates++;
		    }
		}
	    }
		  
	  if(countCandidates > 0)
	    {
	      qsort(candidates, countCandidates, sizeof(qtData), qtCompare);	 	      
	      
	      for(j = 0; j < countCandidates; j++)
		{
		  int k;	       
		  
		  for(k = 0; k < entCount; k++)					     
		    if(dist(current[k], candidates[j].number, sites, nDouble, y) < thres)
		      break;
		  
		  if(k == entCount)
		    {
		      current[entCount++] = candidates[j].number;		     		    		     
		      omitted[candidates[j].number] = 1;
		    }	      	      	
		}
	    }
	  
	  clusters[clusterCount].entries = (int *)malloc(sizeof(int) * entCount);  	      
	  memcpy(clusters[clusterCount].entries, current, entCount * sizeof(int));
	  clusters[clusterCount++].count = entCount;
	}
    }
	    
  printf("Time %f\n", gettime() - t);
  printf("FOUND %d Clusters\n", clusterCount);  
    

   if(1)
    {
      int ver = 0;
      int check = 0;
      int total = 0;
      for(i = 0; i < n; i++)
	ver += i;
      
      for(i = 0; i < clusterCount; i++)
	{	 	  
	  {
	    int k;
	    for(j = 0; j < clusters[i].count; j++)
	      for(k = 0; k < clusters[i].count; k++)	   	    
		assert(dist(clusters[i].entries[j], clusters[i].entries[k],sites, nDouble, y)  >=  thres);		
	  }
	  
	  for(j = 0; j < clusters[i].count; j++)
	    {
	      check += clusters[i].entries[j];	     
	      total++;
	    }	
	}
      assert(ver == check);          
      printf("Total: %d\n", total);
    }



  for(i = 0; i < clusterCount; i++)
    {
      float max  = 0.0;
      int length = clusters[i].count;
      int pos    = -1;      
      int *c     = clusters[i].entries;
      int buf;

      if(length > 2)
	{
	  for(j = 0; j < length; j++)
	    {
	      int k;
	      float avg = 0.0;

	      for(k = 0; k < length; k++)
		{
		  if(j != k)		    		     
		    avg += dist(c[j], c[k], sites, nDouble, y);	
		}	      	     
	       
	      if(avg > max)
		{
		  max = avg;
		  pos = j;
		}
	    }	  
	  if(pos > 0)
	    {	      
	      buf    = c[0];
	      c[0]   = c[pos];
	      c[pos] = buf;	 
	    }
	}
      for(j = 0; j < length; j++)
	c[j] = c[j] + 1;	
    }

  free(candidates);
  free(omitted);
  free(current);
  free(best);
  free(ref);
  *ccc = clusterCount;
  return clusters;
}



static qtList * clusterQT(float **d, int n, float thres, int *ccc)
{  
  int clusterCount;
  int i, j;
  int *omitted, *current, *best;  
  int total = 0;  
  qtList *clusters = (qtList*)NULL;
  double t = gettime();
  
  clusters = (qtList *)malloc(sizeof(qtList) * n);
  omitted = (int*)calloc(n, sizeof(int));
  current = (int*)malloc(sizeof(int) * n);
  best    = (int*)malloc(sizeof(int) * n);            

  clusterCount = 0;

  while(1)
    {
      int max = -1;
      int maxPos = -1;

      for(i = 0; i < n; i++)
	{
	  if(!omitted[i])
	    {
	      int entCount = 0;	  	     
	      int *inSet = (int *)calloc(n, sizeof(int));	      
	      boolean aboveThres = TRUE;

	      current[entCount++] = i;
	      inSet[i] = 1;

	      while(aboveThres)
		{	
		  float dm = -1.0;
		  int dmPos = -1;
	 		 
		  for(j = 0; j < n; j++)		
		    if(i != j && (!omitted[j]) && (!inSet[j]) && d[i][j] > dm)
		      {
			dm    = d[i][j];
			dmPos = j;
		      }
		  
		  if(dmPos == -1)
		    aboveThres = FALSE;
		  else
		    {
		      for(j = 0; j < entCount && aboveThres; j++)			
			if(d[current[j]][dmPos] < thres)
			  aboveThres = FALSE;
		      
		      if(aboveThres)
			{
			  current[entCount++] = dmPos;
			  inSet[dmPos] = 1;
			}
		    }
		}	      	      	      
	      
	      if(entCount > max)
		{
		  max    = entCount;
		  maxPos = i;
		  memcpy(best, current, entCount * sizeof(int));
		}
	      free(inSet);
	    }
	}

      if(maxPos == -1)
	break;

      clusters[clusterCount].entries = (int *)malloc(sizeof(int) * max);
      memcpy(clusters[clusterCount].entries, best, max * sizeof(int));

      for(i = 0; i < max; i++)			 
	omitted[best[i]]    = 1;
	
      clusters[clusterCount++].count = max;               
    }
  
  printf("Time %f\n", gettime() - t);
  printf("FOUND %d Clusters\n", clusterCount);
  
  if(1)
    {
      int ver = 0;
      int check = 0;
      for(i = 0; i < n; i++)
	ver += i;
      
      for(i = 0; i < clusterCount; i++)
	{
	  /*printf("Cluster %d:", i);*/
	  
	  {
	    int k;
	    for(j = 0; j < clusters[i].count; j++)
	      for(k = 0; k < clusters[i].count; k++)	   	    
		assert(d[clusters[i].entries[j]][clusters[i].entries[k]] >=  thres);		
	  }
	  
	  for(j = 0; j < clusters[i].count; j++)
	    {
	      check += clusters[i].entries[j];
	      /*printf("%d ", clusters[i].entries[j]);*/
	      total++;
	    }
	  /*printf("\n");*/
	}
      assert(ver == check);
      
      /*printf("TOTAL: %d\n", total);*/
    }

  for(i = 0; i < clusterCount; i++)
    {
      float max  = 0.0;
      int length = clusters[i].count;
      int pos    = -1;    
      int *c     = clusters[i].entries;
      int buf;

      if(length > 2)
	{
	  for(j = 0; j < length; j++)
	    {
	      int k;
	      float avg = 0.0;
	      for(k = 0; k < length; k++)
		{		  
		  if(j != k)		    		      
		    avg += d[c[j]][c[k]];		   
		}	      	      	      

	      if(avg > max)
		{
		  max = avg;
		  pos = j;
		}
	    }	  
	  /*printf("Cluster %d length %d avg %f\n", i, length, max);*/
	  
	  if(pos > 0)
	    {
	      /*printf("Cluster %d siwtching %d <-> %d\n", i, 0, pos);*/
	      buf    = c[0];
	      c[0]   = c[pos];
	      c[pos] = buf;
	    }
	}
      for(j = 0; j < length; j++)
	c[j] = c[j] + 1;	
    }

  free(omitted);
  free(current);
  free(best);
  *ccc = clusterCount;
  return clusters;
}

static void reduceBySequenceSimilarity(tree *tr, rawdata *rdta, analdef *adef)
{
  int n = tr->mxtips + 1; 
  int i, j;
  int *omissionList     = (int *)malloc(n * sizeof(int));
  int *undeterminedList = (int *)malloc((rdta->sites + 1)* sizeof(int));
  int *modelList        = (int *)malloc((rdta->sites + 1)* sizeof(int));  
  int countNameDuplicates = 0;
  int countUndeterminedColumns = 0;
  int countOnlyGaps = 0;
  int modelCounter = 1;
  char buf[16], outName[1024];
  unsigned char *tipI;
  qtList *clusters = (qtList*)NULL;
  FILE *f, *assoc;  
  int numberOfClusters = 0;
  int nonTrivial = 0;
  
  strcpy(outName,         workdir);  
  strcat(outName,         "RAxML_reducedList.");
  strcat(outName,         run_id);

  if(processID == 0)	      
    f = myfopen(infoFileName, "a");
  else
    f = (FILE *)NULL;

 

  for(i = 1; i < n; i++)       
    omissionList[i] = 0;              

  for(i = 0; i < rdta->sites + 1; i++)
    undeterminedList[i] = 0;
      
  for(i = 1; i < n; i++)
    {
      for(j = i + 1; j < n; j++)
	if(strcmp(tr->nameList[i], tr->nameList[j]) == 0)
	  {
	    countNameDuplicates++;
	    if(processID == 0)
	      {
		printf("Sequence names of taxon %d and %d are identical, they are both called %s\n", i, j, tr->nameList[i]);
		fprintf(f, "Sequence names of taxon %d and %d are identical, they are both called %s\n", i, j, tr->nameList[i]);
	      }
	  }
    }
	  
  if(countNameDuplicates > 0)
    {
      if(processID == 0)
	{
	  printf("ERROR: Found %d taxa that had equal names in the alignment, exiting...\n", countNameDuplicates);
	  fprintf(f, "ERROR: Found %d taxa that had equal names in the alignment, exiting...\n", countNameDuplicates);
	  fclose(f);
	}
      errorExit(-1);
    }

  for(i = 1; i < n; i++)
    {
      j = 1;

      while(j <= rdta->sites)
	{
	  if(tr->dataVector[j] == DNA_DATA && rdta->y[i][j] !=  UNDETERMINED_DNA)
	    break;
	  if(tr->dataVector[j] == AA_DATA && rdta->y[i][j] != UNDETERMINED_AA)
	    break;
	  if(tr->dataVector[j] == BINARY_DATA && rdta->y[i][j] != UNDETERMINED_BINARY)
	    break;
	  assert(tr->dataVector[j] != SECONDARY_DATA && tr->dataVector[j] != SECONDARY_DATA_6 && tr->dataVector[j] != SECONDARY_DATA_7);
	  j++;
	}     

      if(j == (rdta->sites + 1))
	{       
	  if(processID == 0)
	    {
	      printf("ERROR: Sequence %s consists entirely of undetermined values which will be treated as missing data\n",      tr->nameList[i]);
	      fprintf(f, "ERROR: Sequence %s consists entirely of undetermined values which will be treated as missing data\n",      tr->nameList[i]);	      
	    }
	  countOnlyGaps++;
	}
      
    }
  
  if(countOnlyGaps > 0)
    {
      if(processID == 0)
	{
	  printf("ERROR: Found %d sequences that consist entirely of undetermined values, exiting...\n", countOnlyGaps);
	  fprintf(f, "ERROR: Found %d sequences that consist entirely of undetermined values, exiting...\n", countOnlyGaps);
	  fclose(f);
	}
      errorExit(-1);
    }

  for(i = 0; i <= rdta->sites; i++)
    modelList[i] = -1;

  for(i = 1; i <= rdta->sites; i++)
    {    
      j = 1;
     
      while(j < n)
	{
	  if(tr->dataVector[i] == DNA_DATA && rdta->y[j][i] != UNDETERMINED_DNA)
	    break;
	  if(tr->dataVector[i] == AA_DATA && rdta->y[j][i] != UNDETERMINED_AA)
	    break;	
	  if(tr->dataVector[i] == BINARY_DATA && rdta->y[j][i] != UNDETERMINED_BINARY)
	    break;
	  assert(tr->dataVector[i] != SECONDARY_DATA && tr->dataVector[j] != SECONDARY_DATA_6 && tr->dataVector[j] != SECONDARY_DATA_7);
	  j++;
	}
      
      if(j == n)
	{
	  undeterminedList[i] = 1;
	  if(processID == 0)
	    {
	      printf("IMPORTANT WARNING: Alignment column %d contains only undetermined values which will be treated as missing data\n", i);
	      fprintf(f, "IMPORTANT WARNING: Alignment column %d contains only undetermined values which will be treated as missing data\n", i);
	    }
	  countUndeterminedColumns++;	  
	}
      else
	{
	  if(adef->useMultipleModel)
	    {
	      modelList[modelCounter] = tr->model[i];
	      modelCounter++;
	    }
	}
    }  

  switch(adef->similarityFilterMode)
    {
    case SMALL_DATA:
      {
	float **d;   
	int n = tr->mxtips;
	int i, j;
	double t = gettime();
	float nDouble = 1.0 / (float)(rdta->sites);    
	int sites = rdta->sites;
	unsigned char *tipI, *tipJ;
	
	
	d = (float **)malloc(sizeof(float *) * n);
	for(i = 0; i < n; i++)
	  d[i] = (float *)malloc(sizeof(float) * n);
	
	for(i = 0; i < n; i++)
	  {
	    d[i][i] = 1.0;	
	    tipI = &(rdta->y[i + 1][1]);
	    for(j = i + 1; j < n; j++)
	      {
		int k;
		int count = 0;
		tipJ = &(rdta->y[j + 1][1]);	    
		for(k = 0; k < sites; k++)
		  if(tipJ[k] == tipI[k])
		    count++;	   	   	   
		
		d[i][j] = ((float)count * nDouble);
		d[j][i] = d[i][j];
	      }
	  }
	
	printf("DistMat %f\n", gettime() - t);
		   
	t = gettime();
	clusters = clusterQT(d, n, (float)(adef->sequenceSimilarity), &numberOfClusters);
	printf("QT %f %d\n", gettime() - t, numberOfClusters);       
      }
      break;
    case LARGE_DATA:
      {	  
	double t;

	t = gettime();
	clusters = clusterQT_LARGE(tr->mxtips, (float)(adef->sequenceSimilarity), &numberOfClusters, rdta);
	printf("QT %f %d\n", gettime() - t, numberOfClusters);             
      }
      break;
    default:
      assert(0);
    }

  assoc = myfopen(outName, "w");

  for(i = 0; i < numberOfClusters; i++)
    {
      int length = clusters[i].count;
      int *c     = clusters[i].entries;
      int j;
      
      if(length > 1)
	{	 
	  fprintf(assoc, "%s:%s", tr->nameList[c[0]], tr->nameList[c[1]]);
	  for(j = 2; j < length; j++)
	    fprintf(assoc, ",%s", tr->nameList[c[j]]);
	  fprintf(assoc, "\n");
	 
	  nonTrivial++;
	}		      		     
    }

  fclose(assoc);
  

  if(nonTrivial > 0 || countUndeterminedColumns > 0)
    {
      char noDupFile[2048];
      char noDupModels[2048];
              
      if(nonTrivial > 0)
	{
	  if(processID == 0)
	    {
	      printf("\n");	      	   
	      
	      printf("Found %d non-trival clusters, reduction to %d sequences\n", nonTrivial, numberOfClusters);
	      
	      fprintf(f, "\n");
	      
	      fprintf(f, "Found %d non-trival clusters, reduction to %d sequences\n", nonTrivial, numberOfClusters);
	    }
	}
      
      if(countUndeterminedColumns > 0)
	{
	  if(processID == 0)
	    {
	      printf("\n");
	      
	      printf("IMPORTANT WARNING\n");
	      
	      printf("Found %d %s that %s only undetermined values which will be treated as missing data.\n", 
		     countUndeterminedColumns, (countUndeterminedColumns == 1)?"column":"columns", (countUndeterminedColumns == 1)?"contains":"contain");
	      printf("Normally these columns should be excluded from the analysis.\n\n");
	      
	      fprintf(f, "\n");
	      
	      fprintf(f, "IMPORTANT WARNING\n");
	      
	      fprintf(f, "Found %d %s that %s only undetermined values which will be treated as missing data.\n", 
		      countUndeterminedColumns, (countUndeterminedColumns == 1)?"column":"columns", (countUndeterminedColumns == 1)?"contains":"contain");
	      fprintf(f, "Normally these columns should be excluded from the analysis.\n\n");      	  
	    }
	}

      sprintf(buf, "%f", adef->sequenceSimilarity);

      strcpy(noDupFile, seq_file);
      strcat(noDupFile, ".reducedBy.");
      strcat(noDupFile, buf);
      

      strcpy(noDupModels, modelFileName);
      strcat(noDupModels, ".reducedBy.");
      strcat(noDupModels, buf);
	      

      if(processID == 0)
	{
	  if(adef->useMultipleModel && !filexists(noDupModels) && countUndeterminedColumns)
	    {      
	      FILE *newFile = myfopen(noDupModels, "w");

	      printf("\nJust in case you might need it, a mixed model file with \n");
	      printf("model assignments for undetermined columns removed is printed to file %s\n",noDupModels);

	      fprintf(f, "\nJust in case you might need it, a mixed model file with \n");
	      fprintf(f, "model assignments for undetermined columns removed is printed to file %s\n",noDupModels);
	      
 
	      for(i = 0; i < tr->NumberOfModels; i++)
		{
		  boolean modelStillExists = FALSE;
		  
		  for(j = 1; (j <= rdta->sites) && (!modelStillExists); j++)
		    {
		      if(modelList[j] == i)
			modelStillExists = TRUE;
		    }

		  if(modelStillExists)
		    {	  		      
		      int k = 1;
		      int lower, upper;
		      int parts = 0;		      		     		      

		      switch(tr->partitionData[i].dataType)
			{
			case AA_DATA:		      		     
			  {
			    char AAmodel[1024];
			    
			    strcpy(AAmodel, protModels[tr->partitionData[i].protModels]);
			    if(tr->partitionData[i].protFreqs)
			      strcat(AAmodel, "F");		  
			  
			    fprintf(newFile, "%s, ", AAmodel);
			  }
			  break;
			case DNA_DATA:
			  fprintf(newFile, "DNA, ");
			  break;
			case BINARY_DATA:
			  fprintf(newFile, "BIN, ");
			  break;
			default:
			  assert(0);
			}

		      fprintf(newFile, "%s = ", tr->partitionData[i].partitionName);


		      while(k <= rdta->sites)
			{
			  if(modelList[k] == i)
			    {
			      lower = k;
			      while((modelList[k + 1] == i) && (k <= rdta->sites))		      			
				k++;
			      upper = k;
			      
			      if(lower == upper)		  
				{
				  if(parts == 0)
				    fprintf(newFile, "%d", lower);
				  else
				    fprintf(newFile, ",%d", lower);
				}
			      else
				{
				  if(parts == 0)
				    fprintf(newFile, "%d-%d", lower, upper);
				  else
				    fprintf(newFile, ",%d-%d", lower, upper);
				}		  
			      parts++;
			    }
			  k++;
			}
		      fprintf(newFile, "\n");
		    }		  
		}	
	      fclose(newFile);
	    }
	  else
	    {
	      if(adef->useMultipleModel)
		{
		  printf("\n A mixed model file with model assignments for undetermined\n");
		  printf("columns removed has already been printed to  file %s\n",noDupModels);

		  fprintf(f, "\n A mixed model file with model assignments for undetermined\n");
		  fprintf(f, "columns removed has already been printed to  file %s\n",noDupModels);
		}	      
	    }
	     

	  if(!filexists(noDupFile))
	    {
	      FILE *newFile;
	      
	      printf("Just in case you might need it, an alignment file with \n");
	      if(nonTrivial && !countUndeterminedColumns)
		printf("similar sequences removed is printed to file %s\n", noDupFile);
	      if(!nonTrivial && countUndeterminedColumns)
		printf("undetermined columns removed is printed to file %s\n", noDupFile);
	      if(nonTrivial && countUndeterminedColumns)
		printf("similar sequences and undetermined columns removed is printed to file %s\n", noDupFile);
	      
	      fprintf(f, "Just in case you might need it, an alignment file with \n");
	      if(nonTrivial && !countUndeterminedColumns)
		fprintf(f, "similar sequences removed is printed to file %s\n", noDupFile);
	      if(!nonTrivial && countUndeterminedColumns)
		fprintf(f, "undetermined columns removed is printed to file %s\n", noDupFile);
	      if(nonTrivial && countUndeterminedColumns)
		fprintf(f, "similar sequences and undetermined columns removed is printed to file %s\n", noDupFile);
	      
	      newFile = myfopen(noDupFile, "w");
	      
	      fprintf(newFile, "%d %d\n", numberOfClusters, rdta->sites - countUndeterminedColumns);
	      
	      for(i = 0; i < numberOfClusters; i++)
		{
		  
		  fprintf(newFile, "%s ", tr->nameList[clusters[i].entries[0]]);
		  tipI =  &(rdta->y[clusters[i].entries[0]][1]);

		   for(j = 0; j < rdta->sites; j++)
		     {
		       if(undeterminedList[j + 1] == 0)
			 {
			   switch(tr->dataVector[j + 1])
			     {
			     case AA_DATA:
			       fprintf(newFile, "%c", inverseMeaningPROT[tipI[j]]);
			       break;
			     case DNA_DATA:
			       fprintf(newFile, "%c", inverseMeaningDNA[tipI[j]]);
			       break;
			     case BINARY_DATA:
			       fprintf(newFile, "%c", inverseMeaningBINARY[tipI[j]]);
			       break;
			     default:
			       assert(0);
			     }
			 }
		     }
		 
		  fprintf(newFile, "\n");		
		}	      
	      fclose(newFile);	    
	    }
	  else
	    {
	      if(nonTrivial && !countUndeterminedColumns)
		printf("An alignment file with similar sequences removed has already\n");
	      if(!nonTrivial && countUndeterminedColumns)
		printf("An alignment file with undetermined columns removed has already\n");
	      if(nonTrivial && countUndeterminedColumns)
		printf("An alignment file with undetermined columns and similar sequences removed has already\n");
	      
	      printf("been printed to file %s\n",  noDupFile);
	      
	      if(nonTrivial && !countUndeterminedColumns)
		fprintf(f, "An alignment file with similar sequences removed has already\n");
	      if(!nonTrivial && countUndeterminedColumns)
		fprintf(f, "An alignment file with undetermined columns removed has already\n");
	      if(nonTrivial && countUndeterminedColumns)
		fprintf(f, "An alignment file with undetermined columns and similar sequences removed has already\n");
	      
	      fprintf(f, "been printed to file %s\n",  noDupFile);
	    }
	}
    }


  free(undeterminedList);
  free(omissionList);
  free(modelList);
  if(processID == 0)	      
    fclose(f);
}




static void generateBS(tree *tr, analdef *adef)
{
  int i, j, k, w;
  int count;
  char outName[1024], buf[16];
  FILE *of;

  assert(adef->boot != 0);

  for(i = 0; i < adef->multipleRuns; i++)
    {     
      computeNextReplicate(tr, &adef->boot, (int*)NULL, (int*)NULL, FALSE);         

      count = 0;
      for(j = 0; j < tr->cdta->endsite; j++)	
	count += tr->cdta->aliaswgt[j];    	         
      
      assert(count == tr->rdta->sites);
     
      strcpy(outName, workdir);
      strcat(outName, seq_file);
      strcat(outName, ".BS");
      sprintf(buf, "%d", i);
      strcat(outName, buf);
      printf("Printing replicate %d to %s\n", i, outName);

      of = myfopen(outName, "w");

      fprintf(of, "%d %d\n", tr->mxtips, count);  
      
      for(j = 1; j <= tr->mxtips; j++)
	{	
	  unsigned char *tip   =  tr->yVector[tr->nodep[j]->number];	   
	  fprintf(of, "%s ", tr->nameList[j]);	
	  
	  for(k = 0; k < tr->cdta->endsite; k++)
	    {
	      switch(tr->dataVector[k])
		{
		case DNA_DATA:
		  for(w = 0; w < tr->cdta->aliaswgt[k]; w++)		  
		    fprintf(of, "%c", inverseMeaningDNA[tip[k]]);	
		  break;
		case AA_DATA:
		  for(w = 0; w < tr->cdta->aliaswgt[k]; w++)
		    fprintf(of, "%c", inverseMeaningPROT[tip[k]]);
		  break;
		case BINARY_DATA:
		  for(w = 0; w < tr->cdta->aliaswgt[k]; w++)
		    fprintf(of, "%c", inverseMeaningBINARY[tip[k]]);
		  break;
		default:
		  assert(0);
		}
	    }
	  
	  fprintf(of, "\n");	
	}
      fclose(of);
    }  
}



     

static void splitMultiGene(tree *tr, rawdata *rdta)
{
  int i, l;  
  int n = rdta->sites + 1;
  int *modelFilter = (int *)malloc(sizeof(int) * n);
  int length, k;
  unsigned char *tip;
  FILE *outf;
  char outFileName[2048];
  char buf[16];
  
  for(i = 0; i < tr->NumberOfModels; i++)
    {      
      strcpy(outFileName, seq_file);
      sprintf(buf, "%d", i);
      strcat(outFileName, ".GENE.");
      strcat(outFileName, buf);
      outf = myfopen(outFileName, "w");
      length = 0;
      for(k = 1; k < n; k++)
	{
	  if(tr->model[k] == i)
	    {
	      modelFilter[k] = 1;
	      length++;
	    }
	  else
	    modelFilter[k] = -1;
	}

      fprintf(outf, "%d %d\n", rdta->numsp, length);
      
      for(l = 1; l <= rdta->numsp; l++)
	{
	  fprintf(outf, "%s ", tr->nameList[l]);

	  tip = &(rdta->y[l][0]);	    

	  for(k = 1; k < n; k++)
	    {
	      if(modelFilter[k] == 1)
		{
		  switch(tr->dataVector[k])
		    {
		    case AA_DATA:
		      fprintf(outf, "%c", inverseMeaningPROT[tip[k]]);
		      break;
		    case DNA_DATA:
		      fprintf(outf, "%c", inverseMeaningDNA[tip[k]]);
		      break;
		    case BINARY_DATA:
		      fprintf(outf, "%c", inverseMeaningBINARY[tip[k]]);
		      break;
		    default:
		      assert(0);
		    }		 
		}
	    }
	  fprintf(outf, "\n");

	}
      
      fclose(outf);

      printf("Wrote individual gene/partition alignment to file %s\n", outFileName);
    }  
            
  free(modelFilter);
  printf("Wrote all %d individual gene/partition alignments\n", tr->NumberOfModels);
  printf("Exiting normally\n");
}



/* C-OPT this is the function were the memory space for the partial likelihood 
   arrays is assigned. Those array will be accessed all the time and cause the largest 
   amount of cache misses */

#ifndef _USE_PTHREADS

static void allocNodex (tree *tr)
{
  nodeptr  p;
  size_t  
    i,
    j,
    model, 
    offset, 
    memoryRequirements = 0;
  int
    *expArray;
  
  double *likelihoodArray;

  /* TODO include if partition valid etc for Pthreads etc */
  
  /* C-OPT the code below is not too interesting, since this is just the setup 
     for some data structures with low memory footprint that help to keep track
     of the biological model used */

  for(i = 0; i < (size_t)tr->NumberOfModels; i++)
    {
      switch(tr->partitionData[i].dataType)
	{
	case DNA_DATA:
	  tr->partitionData[i].EIGN = (double*)malloc(3 * sizeof(double));
	  tr->partitionData[i].EV   = (double*)malloc(16 * sizeof(double));
	  tr->partitionData[i].EI   = (double*)malloc(12 * sizeof(double));	  
	  tr->partitionData[i].substRates = (double *)malloc(5 * sizeof(double));	  
	  tr->partitionData[i].frequencies =  (double*)malloc(4 * sizeof(double));	  
	  tr->partitionData[i].tipVector   = (double *)malloc(64 * sizeof(double));
	  tr->partitionData[i].symmetryVector = (int *)malloc(6  * sizeof(int));
	  tr->partitionData[i].frequencyGrouping = (int *)malloc(4  * sizeof(int));
	  tr->partitionData[i].nonGTR = FALSE;
	  break;
	case AA_DATA:
	  tr->partitionData[i].EIGN = (double*)malloc(19 * sizeof(double));
	  tr->partitionData[i].EV   = (double*)malloc(400 * sizeof(double));
	  tr->partitionData[i].EI   = (double*)malloc(380 * sizeof(double));	  
	  tr->partitionData[i].tipVector   = (double *)malloc(460 * sizeof(double)); 	  
	  tr->partitionData[i].substRates = (double *)malloc(190 * sizeof(double));	  
	  tr->partitionData[i].frequencies = (double*)malloc(20 * sizeof(double));
	  tr->partitionData[i].symmetryVector = (int *)malloc(190  * sizeof(int));
	  tr->partitionData[i].frequencyGrouping = (int *)malloc(20  * sizeof(int));
	  tr->partitionData[i].nonGTR = FALSE;
	  break;
	case BINARY_DATA:
	  tr->partitionData[i].EIGN = (double*)malloc(1 * sizeof(double));
	  tr->partitionData[i].EV   = (double*)malloc(4 * sizeof(double));
	  tr->partitionData[i].EI   = (double*)malloc(2 * sizeof(double));	  
	  tr->partitionData[i].tipVector   = (double *)malloc(8 * sizeof(double)); 	  
	  tr->partitionData[i].substRates = (double *)malloc(1 * sizeof(double));	  
	  tr->partitionData[i].frequencies = (double*)malloc(2 * sizeof(double));
	  tr->partitionData[i].symmetryVector = (int *)malloc(2  * sizeof(int));
	  tr->partitionData[i].frequencyGrouping = (int *)malloc(2 * sizeof(int));
	  tr->partitionData[i].nonGTR = FALSE;
	  break;
	case SECONDARY_DATA:
	  tr->partitionData[i].EIGN = (double*)malloc(15 * sizeof(double));
	  tr->partitionData[i].EV   = (double*)malloc(256 * sizeof(double));
	  tr->partitionData[i].EI   = (double*)malloc(240 * sizeof(double));	  
	  tr->partitionData[i].tipVector   = (double *)malloc(4096 * sizeof(double)); 	  
	  tr->partitionData[i].substRates = (double *)malloc(120 * sizeof(double));	  
	  tr->partitionData[i].frequencies = (double*)malloc(16 * sizeof(double));
	  tr->partitionData[i].symmetryVector = (int *)malloc(120  * sizeof(int));
	  tr->partitionData[i].frequencyGrouping = (int *)malloc(16  * sizeof(int));
	  tr->partitionData[i].nonGTR = FALSE;
	  break;
	case SECONDARY_DATA_6:
	  tr->partitionData[i].EIGN = (double*)malloc(5 * sizeof(double));
	  tr->partitionData[i].EV   = (double*)malloc(36 * sizeof(double));
	  tr->partitionData[i].EI   = (double*)malloc(30 * sizeof(double));	  
	  tr->partitionData[i].tipVector   = (double *)malloc(384 * sizeof(double)); 	  
	  tr->partitionData[i].substRates = (double *)malloc(15 * sizeof(double));	  
	  tr->partitionData[i].frequencies = (double*)malloc(6 * sizeof(double));
	  tr->partitionData[i].symmetryVector = (int *)malloc(15  * sizeof(int));
	  tr->partitionData[i].frequencyGrouping = (int *)malloc(6  * sizeof(int));
	  tr->partitionData[i].nonGTR = FALSE;
	  break;
	case SECONDARY_DATA_7:
	  tr->partitionData[i].EIGN = (double*)malloc(6 * sizeof(double));
	  tr->partitionData[i].EV   = (double*)malloc(49 * sizeof(double));
	  tr->partitionData[i].EI   = (double*)malloc(42 * sizeof(double));	  
	  tr->partitionData[i].tipVector   = (double *)malloc(896 * sizeof(double)); 	  
	  tr->partitionData[i].substRates = (double *)malloc(21 * sizeof(double));	  
	  tr->partitionData[i].frequencies = (double*)malloc(7 * sizeof(double));
	  tr->partitionData[i].symmetryVector = (int *)malloc(21  * sizeof(int));
	  tr->partitionData[i].frequencyGrouping = (int *)malloc(7  * sizeof(int));
	  tr->partitionData[i].nonGTR = FALSE;
	  break;
	default:
	  assert(0);
	}
      
      tr->partitionData[i].gammaRates = (double*)malloc(sizeof(double) * 4);
      tr->partitionData[i].yVector = (unsigned char **)malloc(sizeof(unsigned char*) * (tr->mxtips + 1));
      tr->partitionData[i].xVector = (double **)malloc(sizeof(double*) * tr->mxtips);
      tr->partitionData[i].pVector = (parsimonyVector **)malloc(sizeof(parsimonyVector*) * tr->mxtips);
      tr->partitionData[i].expVector = (int **)malloc(sizeof(int*) * tr->mxtips);
      tr->partitionData[i].mxtips  = tr->mxtips;

      for(j = 1; j <= (size_t)tr->mxtips; j++)
	tr->partitionData[i].yVector[j] = &(tr->yVector[j][tr->partitionData[i].lower]); 

    }  

  /* C-OPT the stuff below just computes the memory reqs of the particular analysis */

 
  
  for(model = 0; model < (size_t)tr->NumberOfModels; model++)
    {
      size_t width = tr->partitionData[model].upper - tr->partitionData[model].lower;

      switch(tr->partitionData[model].dataType)
	{
	case AA_DATA:
	  memoryRequirements += ((size_t)tr->aaIncrement * width);
	  break;
	case DNA_DATA:
	  memoryRequirements += ((size_t)tr->dnaIncrement * width);
	  break;
	case BINARY_DATA:
	  memoryRequirements += ((size_t)tr->binaryIncrement * width);
	  break;
	case SECONDARY_DATA:
	  memoryRequirements += ((size_t)tr->secondaryIncrement * width);
	  break;
	case SECONDARY_DATA_6:
	  memoryRequirements += ((size_t)tr->secondaryIncrement6 * width);
	  break;
	case SECONDARY_DATA_7:
	  memoryRequirements += ((size_t)tr->secondaryIncrement7 * width);
	  break;
	default:
	  assert(0);
	}
    }

  /* C-OPT this is were the main bulk of memory is allocated, the important 
     large vectores are likelihoodArray and a smaller expArray. This exparray
     is used to scale very very small values such as to avoid numerical underflow 
  */

  /* printf("MemReqs: %ld Overall Mem Reqs %ld\n", memoryRequirements, (size_t)tr->mxtips * memoryRequirements * sizeof(double)); */

  tr->perSiteLL       = (double *)malloc((size_t)tr->cdta->endsite * sizeof(double));
  assert(tr->perSiteLL != NULL);
  likelihoodArray = (double *)malloc((size_t)tr->mxtips * memoryRequirements * sizeof(double));	
  assert(likelihoodArray != NULL);
  expArray = (int *)malloc((size_t)tr->cdta->endsite * (size_t)tr->mxtips * sizeof(int));
  assert(expArray != NULL);
  tr->sumBuffer  = (double *)malloc(memoryRequirements * sizeof(double));
  assert(tr->sumBuffer != NULL);
  

  assert(4 * sizeof(double) > sizeof(parsimonyVector));

  offset = 0;

  /* C-OPT for initial testing tr->NumberOfModels will be 1 */

  for(model = 0; model < (size_t)tr->NumberOfModels; model++)
    {
      size_t lower = tr->partitionData[model].lower;
      size_t width = tr->partitionData[model].upper - lower;

      /* TODO all of this must be reset/adapted when fixModelIndices is called ! */

      tr->partitionData[model].sumBuffer    = &tr->sumBuffer[offset];
      tr->partitionData[model].perSiteLL    = &tr->perSiteLL[lower];
      tr->partitionData[model].wr           = &tr->cdta->wr[lower];
      tr->partitionData[model].wr2          = &tr->cdta->wr2[lower];
      tr->partitionData[model].wgt          = &tr->cdta->aliaswgt[lower];
      tr->partitionData[model].invariant    = &tr->invariant[lower];
      tr->partitionData[model].rateCategory = &tr->cdta->rateCategory[lower];

      switch(tr->partitionData[model].dataType)
	{
	case AA_DATA:
	  offset += ((size_t)tr->aaIncrement * width);
	  break;
	case DNA_DATA:
	  offset += ((size_t)tr->dnaIncrement * width);
	  break;
	case BINARY_DATA:
	  offset += ((size_t)tr->binaryIncrement * width);
	  break;
	case SECONDARY_DATA:
	  offset += ((size_t)tr->secondaryIncrement * width);
	  break;
	case SECONDARY_DATA_6:
	  offset += ((size_t)tr->secondaryIncrement6 * width);
	  break;
	case SECONDARY_DATA_7:
	  offset += ((size_t)tr->secondaryIncrement7 * width);
	  break;
	default:
	  assert(0);
	}
    }
  
  /* C-OPT here we just set a pointer array to point to specific, non-overlapping parts of the 
     large expArray and likelihoodArray vectors. For the time being we don't care about this 
     pVector stuff, the important pointer arrays are: tr->partitionData[model].expVector
     and tr->partitionData[model].xVector[i]. Essentially what happens here is that we assign 
     those vectors to the inner nodes of our phylogenetic tree */

  for(i = 0; i < (size_t)tr->mxtips; i++)
    {
      offset = 0;
      
      for(model = 0; model < (size_t)tr->NumberOfModels; model++)
	{
	  size_t width = tr->partitionData[model].upper - tr->partitionData[model].lower;
	  
	  /*printf("%ld %ld\n", i, i * memoryRequirements + offset);*/

	  tr->partitionData[model].expVector[i] = &expArray[i * tr->cdta->endsite + tr->partitionData[model].lower];	 
	  tr->partitionData[model].xVector[i]   = &likelihoodArray[i * memoryRequirements + offset];
	  tr->partitionData[model].pVector[i]   = (parsimonyVector *)tr->partitionData[model].xVector[i];

	  switch(tr->partitionData[model].dataType)
	    {
	    case AA_DATA:
	      offset += ((size_t)tr->aaIncrement * width);
	      break;
	    case DNA_DATA:
	      offset += ((size_t)tr->dnaIncrement * width);
	      break;
	    case BINARY_DATA:
	      offset += ((size_t)tr->binaryIncrement * width);
	      break;
	    case SECONDARY_DATA:
	      offset += ((size_t)tr->secondaryIncrement * width);
	      break;
	    case SECONDARY_DATA_6:
	      offset += ((size_t)tr->secondaryIncrement6 * width);
	      break;
	    case SECONDARY_DATA_7:
	      offset += ((size_t)tr->secondaryIncrement7 * width);
	      break;
	    default:
	      assert(0);
	    }
	}
    }


  for (i = (size_t)tr->mxtips + 1; (i <= 2*((size_t)tr->mxtips) - 2); i++) 
    {    
      p = tr->nodep[i];                
      p->x = 1;	   
      p->next->x       = 0;
      p->next->next->x = 0;        
    }      
}

#endif


static void initAdef(analdef *adef)
{
  adef->mrpEncoder             = FALSE;
  adef->useSecondaryStructure  = FALSE;
  adef->bootstrapBranchLengths = FALSE;
  adef->model                  = M_GTRCAT;  
  adef->max_rearrange          = 21;  
  adef->stepwidth              = 5;
  adef->initial                = adef->bestTrav = 10;
  adef->initialSet             = FALSE;
  adef->restart                = FALSE;
  adef->mode                   = BIG_RAPID_MODE;
  adef->categories             = 25; 
  adef->boot                   = 0;
  adef->rapidBoot              = 0;
  adef->useWeightFile          = FALSE;
  adef->checkpoints            = 0;
  adef->startingTreeOnly       = 0;
  adef->multipleRuns           = 1;
  adef->useMultipleModel       = FALSE;
  adef->likelihoodEpsilon      = 0.1;
  adef->constraint             = FALSE;
  adef->grouping               = FALSE; 
  adef->randomStartingTree     = FALSE;
  adef->parsimonySeed          = 0;
  adef->proteinMatrix          = JTT;
  adef->protEmpiricalFreqs     = 0;    
  adef->outgroup               = FALSE;
  adef->useInvariant           = FALSE;
  adef->sequenceSimilarity     = 1.0;
  adef->permuteTreeoptimize    = FALSE;
  adef->useInvariant           = FALSE; 
  adef->allInOne               = FALSE; 
  adef->likelihoodTest         = FALSE;
  adef->perGeneBranchLengths   = FALSE;    
  adef->generateBS             = FALSE; 
  adef->bootStopping           = FALSE;
  adef->gapyness               = 0.0;
  adef->similarityFilterMode   = 0;
  adef->useExcludeFile         = FALSE;
  adef->userProteinModel       = FALSE;
  adef->externalAAMatrix       = (double*)NULL;  
  adef->computeELW             = FALSE;
  adef->computeDistance        = FALSE;
  adef->classifyML             = FALSE;
  adef->thoroughInsertion      = FALSE;
  adef->dynamicAlignment       = FALSE;
  adef->compressPatterns       = TRUE;
  adef->printLabelledTree      = FALSE;

}




static int modelExists(char *model, analdef *adef)
{
  int i; 
  char thisModel[1024];
  
  /********** BINARY ********************/

   if(strcmp(model, "BINGAMMAI\0") == 0)
    {
      adef->model = M_BINGAMMA;
      adef->useInvariant = TRUE;
      return 1;
    }  

  if(strcmp(model, "BINGAMMA\0") == 0)
    {
      adef->model = M_BINGAMMA;
      adef->useInvariant = FALSE;
      return 1;
    }
  
  if(strcmp(model, "BINCAT\0") == 0)
    {
      adef->model = M_BINCAT;   
      adef->useInvariant = FALSE;
      return 1;
    }     
  
  if(strcmp(model, "BINCATI\0") == 0)
    {
      adef->model = M_BINCAT;
      adef->useInvariant = TRUE;
      return 1;
    }


  /*********** DNA **********************/

  if(strcmp(model, "GTRGAMMAI\0") == 0)
    {
      adef->model = M_GTRGAMMA;
      adef->useInvariant = TRUE;
      return 1;
    }  

  if(strcmp(model, "GTRGAMMA\0") == 0)
    {
      adef->model = M_GTRGAMMA;
      adef->useInvariant = FALSE;
      return 1;
    }
  
  if(strcmp(model, "GTRCAT\0") == 0)
    {
      adef->model = M_GTRCAT;   
      adef->useInvariant = FALSE;
      return 1;
    }     
  
  if(strcmp(model, "GTRCATI\0") == 0)
    {
      adef->model = M_GTRCAT;
      adef->useInvariant = TRUE;
      return 1;
    }
 

  /*************** AA GTR ********************/
	
  /* TODO empirical FREQS */

  if(strcmp(model, "PROTCATGTR\0") == 0)
    {
      adef->model = M_PROTCAT;
      adef->proteinMatrix = GTR;
      adef->useInvariant = FALSE;
      return 1;
    }

  if(strcmp(model, "PROTCATIGTR\0") == 0)
    {
      adef->model = M_PROTCAT;
      adef->proteinMatrix = GTR;
      adef->useInvariant = TRUE;
      return 1;
    }

  if(strcmp(model, "PROTGAMMAGTR\0") == 0)
    {
      adef->model = M_PROTGAMMA;
      adef->proteinMatrix = GTR;
      adef->useInvariant = FALSE;
      return 1;
    }

  if(strcmp(model, "PROTGAMMAIGTR\0") == 0)
    {
      adef->model = M_PROTGAMMA;
      adef->proteinMatrix = GTR;
      adef->useInvariant = TRUE;
      return 1;
    }  
  
  /****************** AA ************************/
  
  for(i = 0; i < NUM_PROT_MODELS - 1; i++)
    {
      /* check CAT */

      strcpy(thisModel, "PROTCAT");
      strcat(thisModel, protModels[i]);

      if(strcmp(model, thisModel) == 0)
	{
	  adef->model = M_PROTCAT;
	  adef->proteinMatrix = i;
	  return 1;
	}

      /* check CATF */

      strcpy(thisModel, "PROTCAT");
      strcat(thisModel, protModels[i]);
      strcat(thisModel, "F");

      if(strcmp(model, thisModel) == 0)
	{
	  adef->model = M_PROTCAT;
	  adef->proteinMatrix = i;
	  adef->protEmpiricalFreqs = 1;
	  return 1;
	}
      
      /* check CATI */

      strcpy(thisModel, "PROTCATI");
      strcat(thisModel, protModels[i]);

      if(strcmp(model, thisModel) == 0)
	{
	  adef->model = M_PROTCAT;
	  adef->proteinMatrix = i;
	  adef->useInvariant = TRUE;
	  return 1;
	}

      /* check CATIF */

      strcpy(thisModel, "PROTCATI");
      strcat(thisModel, protModels[i]);
      strcat(thisModel, "F");

      if(strcmp(model, thisModel) == 0)
	{
	  adef->model = M_PROTCAT;
	  adef->proteinMatrix = i;
	  adef->protEmpiricalFreqs = 1;
	  adef->useInvariant = TRUE;
	  return 1;
	}
     

      /****************check GAMMA ************************/

      strcpy(thisModel, "PROTGAMMA");
      strcat(thisModel, protModels[i]);

      if(strcmp(model, thisModel) == 0)
	{
	  adef->model = M_PROTGAMMA;
	  adef->proteinMatrix = i;
	  adef->useInvariant = FALSE;
	  return 1;
	}	

      /*check GAMMAI*/

      strcpy(thisModel, "PROTGAMMAI");
      strcat(thisModel, protModels[i]);

      if(strcmp(model, thisModel) == 0)
	{
	  adef->model = M_PROTGAMMA;
	  adef->proteinMatrix = i;
	  adef->useInvariant = TRUE;
	  return 1;
	}


      /* check GAMMAmodelF */

      strcpy(thisModel, "PROTGAMMA");
      strcat(thisModel, protModels[i]);
      strcat(thisModel, "F");
     
      if(strcmp(model, thisModel) == 0)
	{
	  adef->model = M_PROTGAMMA;
	  adef->proteinMatrix = i;
	  adef->protEmpiricalFreqs = 1;
	  adef->useInvariant = FALSE;
	  return 1;
	}	

      /* check GAMMAImodelF */

      strcpy(thisModel, "PROTGAMMAI");
      strcat(thisModel, protModels[i]);
      strcat(thisModel, "F");
      
      if(strcmp(model, thisModel) == 0)
	{
	  adef->model = M_PROTGAMMA;
	  adef->proteinMatrix = i;
	  adef->protEmpiricalFreqs = 1;
	  adef->useInvariant = TRUE;
	  return 1;
	}	
      
    }

  /*********************************************************************************/

  
  
  return 0;
}



static int mygetopt(int argc, char **argv, char *opts, int *optind, char **optarg)
{
  static int sp = 1;
  register int c;
  register char *cp;

  if(sp == 1)
    {
      if(*optind >= argc || argv[*optind][0] != '-' || argv[*optind][1] == '\0')
	return -1;
    }
  else 
    {
      if(strcmp(argv[*optind], "--") == 0) 
	{
	  *optind =  *optind + 1;
	  return -1;
	}
    }

  c = argv[*optind][sp];
  if(c == ':' || (cp=strchr(opts, c)) == 0) 
    {
      printf(": illegal option -- %c \n", c);
      if(argv[*optind][++sp] == '\0') 
	{
	  *optind =  *optind + 1;
	  sp = 1;
	}
      return('?');
    }
  if(*++cp == ':') 
    {
      if(argv[*optind][sp+1] != '\0')
	{
	  *optarg = &argv[*optind][sp+1];
	  *optind =  *optind + 1;
	}
      else 
	{
	  *optind =  *optind + 1;
	  if(*optind >= argc) 
	    {
	      printf(": option requires an argument -- %c\n", c);
	      sp = 1;
	      return('?');
	    } 
	  else
	    {
	      *optarg = argv[*optind];
	      *optind =  *optind + 1;
	    }
	}
      sp = 1;
    } 
  else 
    {
      if(argv[*optind][++sp] == '\0') 
	{
	  sp = 1;
	  *optind =  *optind + 1;
	}
      *optarg = 0;
    }
  return(c);
  }

static void checkOutgroups(tree *tr, analdef *adef)
{
  if(adef->outgroup)
    {
      boolean found;
      int i, j;            

      for(j = 0; j < tr->numberOfOutgroups; j++)
	{
	  found = FALSE;
	  for(i = 1; (i <= tr->mxtips) && !found; i++)
	    {
	      if(strcmp(tr->nameList[i], tr->outgroups[j]) == 0)
		{
		  tr->outgroupNums[j] = i;
		  found = TRUE;
		}
	    }
	  if(!found)
	    {
	      printf("Error, the outgroup name \"%s\" you specified can not be found in the alignment, exiting ....\n", tr->outgroups[j]);
	      errorExit(-1);
	    }
	}
    }
  
}

static void parseOutgroups(char *outgr, tree *tr)
{
  int count = 1, i, k;
  char name[nmlngth];

  i = 0;
  while(outgr[i] != '\0')
    {
      if(outgr[i] == ',')
	count++;
      i++;
    }

  tr->numberOfOutgroups = count;

  tr->outgroups = (char **)malloc(sizeof(char *) * count);

  for(i = 0; i < tr->numberOfOutgroups; i++)   
    tr->outgroups[i] = (char *)malloc(sizeof(char) * nmlngth);    

  tr->outgroupNums = (int *)malloc(sizeof(int) * count);
    
  i = 0;
  k = 0;
  count = 0;
  while(outgr[i] != '\0')
    {
      if(outgr[i] == ',')
	{	
	  name[k] = '\0';
	  strcpy(tr->outgroups[count], name);
	  count++;
	  k = 0;	 
	}
      else
	{
	  name[k] = outgr[i];
	  k++;
	}
      i++;
    }

  name[k] = '\0';
  strcpy(tr->outgroups[count], name);

  /*for(i = 0; i < tr->numberOfOutgroups; i++)
    printf("%d %s \n", i, tr->outgroups[i]);*/


  /*printf("%s \n", name);*/
}


/*********************************** OUTGROUP STUFF END *********************************************************/
static void printVersionInfo(void)
{
  printf("\nThis is %s version %s released by Alexandros Stamatakis in %s\n\n",  programName, programVersion, programDate);
}

static void printMinusFUsage(void)
{
  printf("\n");
  printf("              \"-f a\": rapid Bootstrap analysis and search for best-scoring ML tree in one program run\n");

  printf("              \"-f b\": draw bipartition information on a tree provided with \"-t\" based on multiple trees\n");
  printf("                      (e.g. form a bootstrap) in a file specifed by \"-z\"\n");

  printf("              \"-f c\": check if the alignment can be properly read by RAxML\n");

  printf("              \"-f d\": new rapid hill-climbing \n");
  printf("                      DEFAULT: ON\n");

  printf("              \"-f e\": optimize model+branch lengths for given input tree under GAMMA/GAMMAI only\n"); 

  /*printf("              \"-f f\": olaf option\n");*/

  printf("              \"-f g\": compute per site log Likelihoods for one ore more trees passed via\n");
  printf("                      \"-z\" and write them to a file that can be read by CONSEL\n");
  printf("                        WARNING: does not print likelihoods in the original column order\n");
  printf("              \"-f h\": compute log likelihood test (SH-test) between best tree passed via \"-t\"\n");
  printf("                      and a bunch of other trees passed via \"-z\" \n"); 

  /* printf("              \"-f i\": very fast quick and dirty ML search\n"); */

  printf("              \"-f j\": generate a bunch of bootstrapped alignment files from an original alignemnt file don't\n");
  printf("                      to specify a seed with \"-b\" and the number of replicates with \"-#\" \n");

  printf("              \"-f k\": a posteriori bootstopping analysis using the FC criterion for a tree file containg several bootstrap replicates passed via \"-z\" \n");

  printf("              \"-f l\": a posteriori bootstopping analysis using the WC criterion for a tree file containg several bootstrap replicates passed via \"-z\" \n");

  printf("              \"-f m\": compare bipartitions between two bunches of trees passed via \"-t\" and \"-z\" \n");
  printf("                      respectively. This will return the Pearson correlation between all bipartitions found\n");
  printf("                      in the two tree files. A file called RAxML_bipartitionFrequencies.outpuFileName\n");
  printf("                      will be printed that contains the pair-wise bipartition frequencies of the two sets\n");

  printf("              \"-f n\": compute the log likelihood score of all trees contained in a tree file provided by\n");
  printf("                      \"-z\" under GAMMA or GAMMA+P-Invar\n");

  printf("              \"-f o\": old and slower rapid hill-climbing without heuristic cutoff\n");

  printf("              \"-f p\": perform pure stepwise MP addition of new sequences to an incomplete starting tree and exit\n");
   
  printf("              \"-f q\": (usage not recommended) classify a bunch of environmental sequences into a reference tree using the fast heuristics with dynamic alignment\n");
  printf("                      you will need to start RAxML with a non-comprehensive reference tree and an alignment containing all sequences (reference + query)\n");

  printf("              \"-f r\": compute pairwise Robinson-Foulds (RF) distances between all pairs of trees in a tree file passed via \"-z\" \n");
  printf("                      if the trees have node labales represented as integer support values the program will also compute two flavors of\n");
  printf("                      the weighted Robinson-Foulds (WRF) distance\n");

  printf("              \"-f s\": split up a multi-gene partitioned alignment into the respective subalignments \n");

  printf("              \"-f t\": do randomized tree searches on one fixed starting tree\n");

  printf("              \"-f v\": classify a bunch of environmental sequences into a reference tree using the slow heuristics without dynamic alignment\n");
  printf("                      you will need to start RAxML with a non-comprehensive reference tree and an alignment containing all sequences (reference + query)\n");

  printf("              \"-f w\": compute ELW test on a bunch of trees passed via \"-z\" \n");

  printf("              \"-f x\": compute pair-wise ML distances, ML model parameters will be estimated on an MP \n");
  printf("                      starting tree or a user-defined tree passed via \"-t\", only allowed for GAMMA-based\n");
  printf("                      models of rate heterogeneity\n");

  printf("              \"-f y\": classify a bunch of environmental sequences into a reference tree using the fast heuristics without dynamic alignment\n");
  printf("                      you will need to start RAxML with a non-comprehensive reference tree and an alignment containing all sequences (reference + query)\n");

  /* printf("              \"-f R\": MRP encoder, don't use \n") */

  /* printf("              \"-f V\": will probably disappear again\n"); */

  printf("              \"-f X\": (usage not recommended) classify a bunch of environmental sequences into a reference tree using the slow heuristics with dynamic alignment\n");
  printf("                      you will need to start RAxML with a non-comprehensive reference tree and an alignment containing all sequences (reference + query)\n");

  printf("\n"); 
  printf("              DEFAULT for \"-f\": new rapid hill climbing\n");  
 
  
  printf("\n");
}


static void printREADME(void)
{
  printVersionInfo();  
  printf("\n");
  printf("Please also consult the RAxML-manual\n");
  printf("\nTo report bugs send an email to stamatak@cs.tum.edu\n");
  printf("Please send me all input files, the exact invocation, details of the HW and operating system,\n");
  printf("as well as all error messages printed to screen.\n\n\n");
  
  printf("raxmlHPC[-MPI|-PTHREADS] -s sequenceFileName -n outputFileName -m substitutionModel\n");
  printf("                         [-a weightFileName] [-A secondaryStructureSubstModel]\n");
  printf("                         [-b bootstrapRandomNumberSeed] [-B wcCriterionThreshold]\n");
  printf("                         [-c numberOfCategories] [-d]\n");
  printf("                         [-e likelihoodEpsilon] [-E excludeFileName]\n");
  printf("                         [-f a|b|c|d|e|g|h|j|k|l|m|n|o|p|q|r|s|t|v|w|x|y|X]\n");
  printf("                         [-g groupingFileName] [-h] [-i initialRearrangementSetting] [-j] [-k] \n");
  printf("                         [-l sequenceSimilarityThreshold] [-L sequenceSimilarityThreshold] [-M]\n"); 
  printf("                         [-o outGroupName1[,outGroupName2[,...]]] \n");
  printf("                         [-p parsimonyRandomSeed] [-P proteinModel]\n");
  printf("                         [-q multipleModelFileName] [-r binaryConstraintTree]\n");
  printf("                         [-S secondaryStructureFile] [-t userStartingTree]\n"); 
  printf("                         [-T numberOfThreads] [-v] [-w workingDirectory]\n");
  printf("                         [-x rapidBootstrapRandomNumberSeed] [-y] [-Y]\n");
  printf("                         [-z multipleTreesFile] [-#|-N numberOfRuns|autoWC|autoFC]\n");
  printf("\n");
  printf("      -a      Specify a column weight file name to assign individual weights to each column of \n");
  printf("              the alignment. Those weights must be integers separated by any type and number \n");
  printf("              of whitespaces whithin a separate file, see file \"example_weights\" for an example.\n");
  printf("\n");
  printf("      -A      Specify one of the secondary structure substitution models implemented in RAxML.\n");
  printf("              The same nomenclature as in the PHASE manual is used, available models: \n"); 
  printf("              S6A, S6B, S6C, S6D, S6E, S7A, S7B, S7C, S7D, S7E, S7F, S16, S16A, S16B\n");
  printf("\n");
  printf("              DEFAULT: 16-state GTR model (S16)\n");
  printf("\n");
  printf("      -b      Specify an integer number (random seed) and turn on bootstrapping\n");
  printf("\n");
  printf("              DEFAULT: OFF\n");
  printf("\n");
  printf("      -B      specify a floating point number between 0.0 and 1.0 that will be used as cutoff threshold \n");
  printf("              for the WC bootstopping criterion. The recommended setting is 0.03.\n");
  printf("\n");
  printf("              DEFAULT: 0.03 (recommended empirically determined setting)\n");
  printf("\n");
  printf("      -c      Specify number of distinct rate catgories for RAxML when modelOfEvolution\n");
  printf("              is set to GTRCAT or GTRMIX\n");
  printf("              Individual per-site rates are categorized into numberOfCategories rate \n");
  printf("              categories to accelerate computations. \n");
  printf("\n");
  printf("              DEFAULT: 25\n");
  printf("\n");
  printf("      -d      start ML optimization from random starting tree \n");
  printf("\n");
  printf("              DEFAULT: OFF\n");
  printf("\n");
  printf("      -e      set model optimization precision in log likelihood units for final\n"); 
  printf("              optimization of tree topology under MIX/MIXI or GAMMA/GAMMAI\n");
  printf("\n");
  printf("              DEFAULT: 0.1   for models not using proportion of invariant sites estimate\n");
  printf("                       0.001 for models using proportion of invariant sites estimate\n");
  printf("\n");
  printf("      -E      specify an exclude file name, that contains a specification of alignment positions you wish to exclude.\n");
  printf("              Format is similar to Nexus, the file shall contain entries like \"100-200 300-400\", to exclude a\n");
  printf("              single column write, e.g., \"100-100\", if you use a mixed model, an appropriatly adapted model file\n");
  printf("              will be written.\n"); 
  printf("\n");
  printf("      -f      select algorithm:\n");

  printMinusFUsage();
 
  printf("\n");
  printf("      -g      specify the file name of a multifurcating constraint tree\n");
  printf("              this tree does not need to be comprehensive, i.e. must not contain all taxa\n");
  printf("\n");
  printf("      -h      Display this help message.\n");  
  printf("\n");
  printf("      -i      Initial rearrangement setting for the subsequent application of topological \n");
  printf("              changes phase\n");
  printf("\n");
  printf("              DEFAULT: determined by program\n");
  printf("\n");
  printf("      -j      Specifies if checkpoints will be written by the program. If checkpoints \n");
  printf("              (intermediate tree topologies) shall be written by the program specify \"-j\"\n");
  printf("\n");
  printf("              DEFAULT: OFF\n");
  printf("\n");
  printf("      -k      Specifies that bootstrapped trees should be printed with branch lengths.\n");
  printf("              The bootstraps will run a bit longer, because model parameters will be optimized\n");
  printf("              at the end of each run. Use with CATMIX/PROTMIX or GAMMA/GAMMAI.\n");
  printf("\n");
  printf("              DEFAULT: OFF\n");   
  printf("\n");
  printf("      -l      Specify a threshold for sequence similarity clustering. RAxML will then print out an alignment\n");
  printf("              to a file called sequenceFileName.reducedBy.threshold that only contains sequences <= the\n");
  printf("              specified thresold that must be between  0.0 and 1.0. RAxML uses the QT-clustering algorithm \n");
  printf("              to perform this task. In addition, a file called RAxML_reducedList.outputFileName will be written\n");
  printf("              that contains clustering information.\n");
  printf("\n");
  printf("              DEFAULT: OFF\n");   
  printf("\n");
  printf("      -L      Same functionality as \"-l\" above, but uses a less exhasutive and thus faster clustering algorithm\n");
  printf("              This is intended for very large datasets with more than 20,000-30,000 sequences\n");
  printf("\n");
  printf("              DEFAULT: OFF\n"); 
  printf("\n");
  printf("      -m      Model of Binary (Morphological), Nucleotide or Amino Acid Substitution: \n");  
  printf("\n");
  printf("              BINARY:\n\n");
  printf("                \"-m BINCAT\"        : Optimization of site-specific\n");
  printf("                                     evolutionary rates which are categorized into numberOfCategories distinct \n");
  printf("                                     rate categories for greater computational efficiency. Final tree might be evaluated\n");
  printf("                                     automatically under BINGAMMA, depending on the tree search option\n"); 
  printf("                \"-m BINCATI\"       : Optimization of site-specific\n");
  printf("                                     evolutionary rates which are categorized into numberOfCategories distinct \n");
  printf("                                     rate categories for greater computational efficiency. Final tree might be evaluated\n");
  printf("                                     automatically under BINGAMMAI, depending on the tree search option \n"); 
  printf("                \"-m BINGAMMA\"      : GAMMA model of rate \n");
  printf("                                     heterogeneity (alpha parameter will be estimated)\n");
  printf("                \"-m BINGAMMAI\"     : Same as BINGAMMA, but with estimate of proportion of invariable sites\n");
  printf("\n");
  printf("              NUCLEOTIDES:\n\n");
  printf("                \"-m GTRCAT\"        : GTR + Optimization of substitution rates + Optimization of site-specific\n");
  printf("                                     evolutionary rates which are categorized into numberOfCategories distinct \n");
  printf("                                     rate categories for greater computational efficiency.  Final tree might be evaluated\n");
  printf("                                     under GTRGAMMA, depending on the tree search option\n");
  printf("                \"-m GTRCATI\"       : GTR + Optimization of substitution rates + Optimization of site-specific\n");
  printf("                                     evolutionary rates which are categorized into numberOfCategories distinct \n");
  printf("                                     rate categories for greater computational efficiency.  Final tree might be evaluated\n");
  printf("                                     under GTRGAMMAI, depending on the tree search option\n");
  printf("                \"-m GTRGAMMA\"      : GTR + Optimization of substitution rates + GAMMA model of rate \n");
  printf("                                     heterogeneity (alpha parameter will be estimated)\n");
  printf("                \"-m GTRGAMMAI\"     : Same as GTRGAMMA, but with estimate of proportion of invariable sites \n");  
  printf("\n");
  printf("              AMINO ACIDS:\n\n");          
  printf("                \"-m PROTCATmatrixName[F]\"        : specified AA matrix + Optimization of substitution rates + Optimization of site-specific\n");
  printf("                                                   evolutionary rates which are categorized into numberOfCategories distinct \n");
  printf("                                                   rate categories for greater computational efficiency.   Final tree might be evaluated\n");
  printf("                                                   automatically under PROTGAMMAmatrixName[f], depending on the tree search option\n");
  printf("                \"-m PROTCATImatrixName[F]\"        : specified AA matrix + Optimization of substitution rates + Optimization of site-specific\n");
  printf("                                                   evolutionary rates which are categorized into numberOfCategories distinct \n");
  printf("                                                   rate categories for greater computational efficiency.   Final tree might be evaluated\n");
  printf("                                                   automatically under PROTGAMMAImatrixName[f], depending on the tree search option\n");  
  printf("                \"-m PROTGAMMAmatrixName[F]\"      : specified AA matrix + Optimization of substitution rates + GAMMA model of rate \n");
  printf("                                                   heterogeneity (alpha parameter will be estimated)\n");
  printf("                \"-m PROTGAMMAImatrixName[F]\"     : Same as PROTGAMMAmatrixName[F], but with estimate of proportion of invariable sites \n");  
  printf("\n");
  printf("                Available AA substitution models: DAYHOFF, DCMUT, JTT, MTREV, WAG, RTREV, CPREV, VT, BLOSUM62, MTMAM, LG, GTR\n");
  printf("                With the optional \"F\" appendix you can specify if you want to use empirical base frequencies\n");
  printf("                Please note that for mixed models you can in addition specify the per-gene AA model in\n");
  printf("                the mixed model file (see manual for details). Also note that if you estimate AA GTR parameters on a partitioned\n");
  printf("                dataset, they will be linked (estimated jointly) across all partitions to avoid over-parametrization\n");
  printf("\n");
  printf("      -M      Switch on estimation of individual per-partition branch lengths. Only has effect when used in combination with \"-q\"\n");  
  printf("              Branch lengths for individual partitions will be printed to separate files\n");
  printf("              A weighted average of the branch lengths is computed by using the respective partition lengths\n");
  printf("\n"),
  printf("              DEFAULT: OFF\n");  
  printf("\n");
  printf("      -n      Specifies the name of the output file.\n"); 
  printf("\n"); 
  printf("      -o      Specify the name of a single outgrpoup or a comma-separated list of outgroups, eg \"-o Rat\" \n");
  printf("              or \"-o Rat,Mouse\", in case that multiple outgroups are not monophyletic the first name \n");
  printf("              in the list will be selected as outgroup, don't leave spaces between taxon names!\n");
  printf("\n");   
  printf("      -p      Specify a random number seed for the parsimony inferences. This allows you to reproduce your results\n");
  printf("              and will help me debug the program. This option HAS NO EFFECT in the parallel MPI version\n");
  printf("\n");
  printf("      -P      Specify the file name of a user-defined AA (Protein) substitution model. This file must contain\n");
  printf("              420 entries, the first 400 being the AA substitution rates (this must be a symmetric matrix) and the\n");
  printf("              last 20 are the empirical base frequencies\n");
  printf("\n"); 
  printf("      -q      Specify the file name which contains the assignment of models to alignment\n"); 
  printf("              partitions for multiple models of substitution. For the syntax of this file\n");
  printf("              please consult the manual.\n");    
  printf("\n");
  printf("      -r      Specify the file name of a binary constraint tree.\n");
  printf("              this tree does not need to be comprehensive, i.e. must not contain all taxa\n");
  printf("\n");  
  printf("      -s      Specify the name of the alignment data file in PHYLIP format\n");
  printf("\n");
  printf("      -S      Specify the name of a secondary structure file. The file can contain \".\" for \n");
  printf("              alignment columns that do not form part of a stem and characters \"()<>[]{}\" to define \n");
  printf("              stem regions and pseudoknots\n");
  printf("\n");
  printf("      -t      Specify a user starting tree file name in Newick format\n");
  printf("\n");
  printf("      -T      PTHREADS VERSION ONLY! Specify the number of threads you want to run.\n");
  printf("              Make sure to set \"-T\" to at most the number of CPUs you have on your machine,\n");
  printf("              otherwise, there will be a huge performance decrease!\n");
  printf("\n");
  printf("      -v      Display version information\n");  
  printf("\n");
  printf("      -w      Name of the working directory where RAxML will write its output files\n");
  printf("\n");
  printf("              DEFAULT: current directory\n"); 
  printf("\n");
  printf("      -x      Specify an integer number (random seed) and turn on rapid bootstrapping\n");
  printf("              CAUTION: unlike in version 7.0.4 RAxML will conduct rapid BS replicates under \n");
  printf("              the model of rate heterogeneity you specified via \"-m\" and not by default under CAT\n");
  printf("\n");
  printf("      -y      If you want to only compute a parsimony starting tree with RAxML specify \"-y\",\n");
  printf("              the program will exit after computation of the starting tree\n");
  printf("\n");
  printf("              DEFAULT: OFF\n");
  printf("\n");  
  printf("      -Y      Do a more thorough parsimony tree search using a parsimony ratchet and exit. \n");
  printf("              specify the number of ratchet searches via \"-#\" or \"-N\"\n");
  printf("              This has just been implemented for completeness, if you want a fast MP implementation use TNT\n");
  printf("\n");
  printf("              DEFAULT: OFF\n");
  printf("\n");
  printf("      -z      Specify the file name of a file containing multiple trees e.g. from a bootstrap\n");
  printf("              that shall be used to draw bipartition values onto a tree provided with \"-t\",\n");
  printf("              It can also be used to compute per site log likelihoods in combination with \"-f g\"\n");
  printf("              and to read a bunch of trees for a couple of other options (\"-f h\", \"-f m\", \"-f n\").\n");
  printf("\n"); 
  printf("      -#|-N   Specify the number of alternative runs on distinct starting trees\n");
  printf("              In combination with the \"-b\" option, this will invoke a multiple boostrap analysis\n");
  printf("              Note that \"-N\" has been added as an alternative since \"-#\" sometimes caused problems\n");
  printf("              with certain MPI job submission systems, since \"-#\" is often used to start comments\n");
  printf("              If you want to use the bootstopping criteria specify \"-# autoWC\" for the weighted criterion\n");
  printf("              or \"-# autoFC\" for the frequency-based criterion. This will only work in combination with \"-x\"\n");
  printf("              or \"-b\"\n");
  printf("\n");
  printf("              DEFAULT: 1 single analysis\n");
  printf("\n\n\n\n");

}




  
  


static void get_args(int argc, char *argv[], analdef *adef, tree *tr)
{
  int	optind = 1;
  int        c;
  boolean    
    bad_opt    =FALSE,
    workDirSet = FALSE;
  char       aut[256];
  char       buf[2048];
  char       *optarg;
  char       model[2048] = ""; 
  char       secondaryModel[2048] = "";
  char       modelChar;
  double likelihoodEpsilon, sequenceSimilarity, wcThreshold;
  int nameSet = 0, 
    alignmentSet = 0,    
    multipleRuns = 0,    
    constraintSet = 0,
    treeSet = 0,
    groupSet = 0,
    modelSet = 0,
    treesSet  = 0;
  long parsimonySeed = 0;
  run_id[0] = 0;
  workdir[0] = 0;
  seq_file[0] = 0;
  tree_file[0] = 0;
  model[0] = 0;
  weightFileName[0] = 0;
  modelFileName[0] = 0;

  /*********** tr inits **************/

#ifdef _USE_PTHREADS
  NumberOfThreads = 0;
#endif
  tr->bootStopCriterion = -1;
  tr->wcThreshold = 0.03;
  tr->doCutoff = TRUE;
  tr->secondaryStructureModel = SEC_16; /* default setting */
 
  /********* tr inits end*************/


  while(!bad_opt && 
	((c = mygetopt(argc,argv,"T:E:N:B:l:x:z:g:r:e:a:b:c:f:i:m:t:w:s:n:o:L:P:S:A:q:#:p:vdyjhkMY", &optind, &optarg))!=-1))
    {
    switch(c) 
      { 
      case 'A':
	{
	  const char *modelList[21] = { "S6A", "S6B", "S6C", "S6D", "S6E", "S7A", "S7B", "S7C", "S7D", "S7E", "S7F", "S16", "S16A", "S16B", "S16C", 
				      "S16D", "S16E", "S16F", "S16I", "S16J", "S16K"};
	  int i;

	  sscanf(optarg, "%s", secondaryModel);

	  for(i = 0; i < 21; i++)
	    if(strcmp(secondaryModel, modelList[i]) == 0)
	      break;

	  if(i < 21)
	    tr->secondaryStructureModel = i;	 
	  else
	    {
	      printf("The secondary structure model %s you want to use does not exist, exiting .... \n", secondaryModel);
	      errorExit(0);
	    }	       	   
	}
	break;
      case 'B':
	sscanf(optarg,"%lf", &wcThreshold);
	tr->wcThreshold = wcThreshold;
	if(wcThreshold <= 0.0 || wcThreshold >= 1.0)
	  {
	    printf("\nBootstrap threshold must be set to values between 0.0 and 1.0, you just set it to %f\n", wcThreshold);
	    printf("Are you one of those lazy guys who don't read the manual and think that stamatak@cs.tum.edu is actually the phylogenetics help desk?\n\n");
	    exit(-1);
	  }
	if(wcThreshold < 0.01 || wcThreshold > 0.05)
	  {
	    printf("\n\nWARNING, reasonable settings for Bootstopping threshold with WC criterion range between 0.01 and 0.05.\n");
	    printf("You are just setting it to %f, are you one of those lazy guys who don't read the manual?\n\n", wcThreshold);
	  }		
	break;
      case 'S':
	adef->useSecondaryStructure = TRUE;
	strcpy(secondaryStructureFileName, optarg);	
	break;
      case 'T':		
#ifdef _USE_PTHREADS
	sscanf(optarg,"%d", &NumberOfThreads);	
#else
	if(processID == 0)
	  {
	    printf("Option -T does not have any effect with the sequential or parallel MPI version.\n");
	    printf("It is used to specify the number of threads for the Pthreads-based parallelization\n");
	  }	 
#endif        
	break;
      case 'P':   
	strcpy(proteinModelFileName, optarg);
	adef->userProteinModel = TRUE;
	parseProteinModel(adef);
	break;
      case 'E':
	strcpy(excludeFileName, optarg);
	adef->useExcludeFile = TRUE;
	break;
      case 'M':
	adef->perGeneBranchLengths = TRUE;
	break;     
      case 'o':
	{
	  char *outgroups;
	  outgroups = (char*)malloc(sizeof(char) * (strlen(optarg) + 1));
	  strcpy(outgroups, optarg);
	  parseOutgroups(outgroups, tr);
	  free(outgroups);
	  adef->outgroup = TRUE;
	}
	break;	
      case 'k':
	adef->bootstrapBranchLengths = TRUE;
	break;
      case 'z':	
	strcpy(bootStrapFile, optarg);
	treesSet = 1;
	break;               
      case 'd':
	adef->randomStartingTree = TRUE;
	break;
      case 'g':
	strcpy(tree_file, optarg);
	adef->grouping = TRUE;
	adef->restart  = TRUE;      
	groupSet = 1;
	break;
      case 'r':	
	strcpy(tree_file, optarg);
	adef->restart = TRUE;
	adef->constraint = TRUE;
	constraintSet = 1;
	break;
      case 'e':      
	sscanf(optarg,"%lf", &likelihoodEpsilon);
	adef->likelihoodEpsilon = likelihoodEpsilon;      
	break;
      case 'q':
	strcpy(modelFileName,optarg);
	adef->useMultipleModel = TRUE;      
        break;       
      case 'p':
	sscanf(optarg,"%ld", &parsimonySeed);
	adef->parsimonySeed = parsimonySeed;
	break;
      case 'N':
      case '#':
	/* TODO include auto in readme */
	if(sscanf(optarg,"%d", &multipleRuns) > 0)
	  {
	    adef->multipleRuns = multipleRuns;
	  }
	else
	  {
	    if((sscanf(optarg,"%s", aut) > 0) && ((strcmp(aut, "autoFC") == 0) || (strcmp(aut, "autoWC") == 0)))
	      {
		if((strcmp(aut, "autoFC") == 0))
		  {
		    adef->bootStopping = TRUE;
		    adef->multipleRuns = 1000;
		    tr->bootStopCriterion = FREQUENCY_STOP;
		  }	    
		if((strcmp(aut, "autoWC") == 0))
		  {
		    adef->bootStopping = TRUE;
		    adef->multipleRuns = 1000;
		    tr->bootStopCriterion = WC_STOP;
		  }			       		
	      }
	    else
	      {
		if(processID == 0)
		  {
		    printf("Use -# or -N option either with an integer, e.g., -# 100 or with -# autoFC or -# autoWC\n");
		    printf("or -N 100 or  -N autoFC or -N autoWC respectively, note that auto will not work for the\n");
		    printf("MPI-based parallel version\n");
		  }		
		errorExit(0);
	      }
	  }
	break;
      case 'v':
	printVersionInfo();
	errorExit(0);
      case 'y':
	adef->startingTreeOnly = 1;
	break;
      case 'Y':		
	adef->mode = THOROUGH_PARSIMONY;
	break;
      case 'h':
	printREADME();
	errorExit(0);
      case 'j':	
	adef->checkpoints = 1;
	break;
      case 'a':
	strcpy(weightFileName,optarg);
	adef->useWeightFile = TRUE;
        break;            
      case 'b':
	sscanf(optarg,"%ld", &adef->boot);             
	break;     
      case 'x':       
	sscanf(optarg,"%ld", &adef->rapidBoot);   
	break;
      case 'c':
	sscanf(optarg, "%d", &adef->categories);       
	break;	  
      case 'l':
	sscanf(optarg,"%lf", &sequenceSimilarity);
	adef->sequenceSimilarity = sequenceSimilarity; 
	adef->mode = SEQUENCE_SIMILARITY_FILTER;
	adef->similarityFilterMode = SMALL_DATA;	
	break;
      case 'L':
	sscanf(optarg,"%lf", &sequenceSimilarity);
	adef->sequenceSimilarity = sequenceSimilarity; 
	adef->mode = SEQUENCE_SIMILARITY_FILTER;
	adef->similarityFilterMode = LARGE_DATA;		
	break;     
      case 'f': 
	sscanf(optarg, "%c", &modelChar);
	switch(modelChar)
	  {
	  case 'a':
	    adef->allInOne = TRUE;
	    adef->mode = BIG_RAPID_MODE;
	    tr->doCutoff = TRUE;
	    break;
	  case 'b': 
	    adef->mode = CALC_BIPARTITIONS; 
	    break;
	  case 'c':
	    adef->mode = CHECK_ALIGNMENT;
	    break;
	  case 'd': 
	    adef->mode = BIG_RAPID_MODE;
	    tr->doCutoff = TRUE;
	    break;
	  case 'e': 
	    adef->mode = TREE_EVALUATION; 
	    break;
	  case 'f':
	    adef->mode = OLAF_OPTION;
	    break;
	  case 'g':
	    adef->mode              = PER_SITE_LL;	  
	    break;
	  case 'h':
	    adef->mode = TREE_EVALUATION;
	    adef->likelihoodTest = TRUE;
	    break;
	  case 'i':
	    adef->mode = SUPER_FAST;
	    break;
	  case 'j':
	    adef->mode = GENERATE_BS;
	    adef->generateBS = TRUE;
	    break;
	  case 'k':	    
	    adef->mode = BOOTSTOP_ONLY;	  
	    tr->bootStopCriterion = FREQUENCY_STOP;
	    break;
	  case 'l':
	    adef->mode = BOOTSTOP_ONLY;	  
	    tr->bootStopCriterion = WC_STOP;
	    break;
	  case 'm':
	    adef->mode = COMPUTE_BIPARTITION_CORRELATION;	   
	    break;
	  case 'n':
	    adef->mode = COMPUTE_LHS;	    
	    break;
	  case 'o': 
	    adef->mode = BIG_RAPID_MODE;
	    tr->doCutoff = FALSE;
	    break;
	  case 'p':
	    adef->mode =  PARSIMONY_ADDITION;
	    break;
	  case 'r':
	    adef->mode = COMPUTE_RF_DISTANCE;
	    break;
	  case 's':
	    adef->mode = SPLIT_MULTI_GENE;
	    break;       	  	 
	  case 't':
	    adef->mode = BIG_RAPID_MODE;
	    tr->doCutoff = TRUE;
	    adef->permuteTreeoptimize = TRUE;	    
	    break;	  	  
	  case 'w':
	    workDirSet = TRUE;
	    adef->mode = COMPUTE_ELW;
	    adef->computeELW = TRUE;
	    break;
	  case 'x':
	    adef->mode = DISTANCE_MODE;
	    adef->computeDistance = TRUE;	    
	    break;
	  case 'q':
	    /* TODO include in README */	  
	    adef->mode = CLASSIFY_ML;
	    adef->classifyML = TRUE;
	    adef->thoroughInsertion = FALSE;
	    adef->dynamicAlignment  = TRUE;
	    adef->compressPatterns  = FALSE;
	    break;
	  case 'V':
	    adef->mode = CLASSIFY_ML;
	    adef->classifyML = TRUE;
	    adef->thoroughInsertion = TRUE;
	    adef->dynamicAlignment  = FALSE;
	    adef->printLabelledTree = TRUE;
	    break;
	  case 'v':
	    /* TODO README */
	    adef->mode = CLASSIFY_ML;
	    adef->classifyML = TRUE;
	    adef->thoroughInsertion = TRUE;
	    adef->dynamicAlignment  = FALSE;
	    break;
	  case 'y':
	    adef->mode = CLASSIFY_ML;
	    adef->classifyML = TRUE;	   
	    adef->thoroughInsertion = FALSE;
	    adef->dynamicAlignment  = FALSE;	   
	    break;
	  case 'X':
	    adef->mode = CLASSIFY_ML;
	    adef->classifyML = TRUE;
	    adef->thoroughInsertion = TRUE;
	    adef->dynamicAlignment  = TRUE;
	    adef->compressPatterns  = FALSE;
	    break; 	 
	  case 'R':
	    adef->mrpEncoder = TRUE;
	    break; 
	  default: 
	    {
	      if(processID == 0)
		{
		  printf("Error select one of the following algorithms via -f :\n");
		  printMinusFUsage();		 		 
		}	     
	      errorExit(-1);
	    }
	  }
	break;      
      case 'i':
	sscanf(optarg, "%d", &adef->initial);
	adef->initialSet = TRUE;
	break;     
      case 'n':
        strcpy(run_id,optarg);
	nameSet = 1;
        break;
      case 'w':
        strcpy(workdir,optarg);
        break;                 
      case 't':
	strcpy(tree_file, optarg);
	adef->restart = TRUE;
	treeSet = 1;
	break;
      case 's':
	strcpy(seq_file, optarg);
	alignmentSet = 1;
	break;
      case 'm':
	strcpy(model,optarg);
	if(modelExists(model, adef) == 0)
	  {
	    if(processID == 0)
	      {
		printf("Model %s does not exist\n\n", model);
		printf("For DNA data use: GTRCAT        or GTRGAMMA      or\n");
		printf("                  GTRMIX        or GTRMIXI       or\n");
		printf("                  GTRGAMMAI     or GTRCAT_GAMMAI or\n");
		printf("                  GTRCAT_GAMMA\n\n");
		printf("For AA data use:  PROTCATmatrixName[F]        or PROTGAMMAmatrixName[F]      or\n");
		printf("                  PROTMIXmatrixName[F]        or PROTMIXImatrixName[F]       or\n");
		printf("                  PROTGAMMAImatrixName[F]     or PROTCAT_GAMMAImatrixName[F] or\n");
		printf("                  PROTCAT_GAMMAImatrixName[F]\n\n");                
		printf("The AA substitution matrix can be one of the following: \n");
		printf("DAYHOFF, DCMUT, JTT, MTREV, WAG, RTREV, CPREV, VT, BLOSUM62, MTMAM, GTR\n\n");
		printf("With the optional \"F\" appendix you can specify if you want to use empirical base frequencies\n");
		printf("Please note that for mixed models you can in addition specify the per-gene model in\n");
		printf("the mixed model file (see manual for details)\n");	    
	      }
	    errorExit(-1);
	  }      
	else	  
	  modelSet = 1;
	break;      
      default:	   
	errorExit(-1);              
    }    
  }     

#ifdef _USE_PTHREADS  
  if(NumberOfThreads < 2)
    {
      printf("\nThe number of threads is currently set to %d\n", NumberOfThreads);
      printf("Specify the number of threads to run via -T numberOfThreads\n");
      printf("NumberOfThreads must be set to an integer value greater than 1\n\n");
      errorExit(-1);
    }
#endif 




  if(adef->computeELW)
    {
      if(processID == 0)
	{
	  if(adef->boot == 0)
	    {
	      printf("Error, you must specify a bootstrap seed via \"-b\" to compute ELW statistics\n");
	      errorExit(-1);
	    }

	  if(adef->multipleRuns < 2)
	    {
	      printf("Error, you must specify the number of BS replicates via \"-#\" or \"-N\" to compute ELW statistics\n");
	      printf("it should be larger than 1, recommended setting is 100\n");
	      errorExit(-1);
	    }
	  
	  if(!treesSet)
	    {
	      printf("Error, you must specify an input file containing several candidate trees\n");
	      printf("via \"-z\" to compute ELW statistics.\n");
	      errorExit(-1);
	    }

	  if(!isGamma(adef))
	    {
	      printf("Error ELW test can only be conducted undetr GAMMA or GAMMA+P-Invar models\n");
	      errorExit(-1);
	    }
	}
    }

 
      


  if(((!adef->boot) && (!adef->rapidBoot)) && adef->bootStopping)
    {
      if(processID == 0)
	{
	  printf("Can't use automatic bootstopping without actually doing a Bootstrap\n");
	  printf("Specify either -x randomNumberSeed (rapid) or -b randomNumberSeed (standard)\n");
	  errorExit(-1);
	}      
    }

  if(adef->boot && adef->rapidBoot)
    {
      if(processID == 0)
	{
	  printf("Can't use standard and rapid BOOTSTRAP simultaneously\n");
	  errorExit(-1);
	}
    }

  if(adef->rapidBoot && !(adef->classifyML))
    {
      if(processID == 0 && (adef->restart || treesSet) && !(groupSet || constraintSet))
	{
	  printf("Error, starting tree(s) will be ignored by rapid Bootstrapping\n");
	  errorExit(-1);
	}     
    }

  if(adef->allInOne && (adef->rapidBoot == 0))
    {
      if(processID == 0)
	{
	  printf("Error, to carry out an ML search after a rapid BS inference you must specify a random number seed with -x\n");
	  errorExit(-1);
	}
    }

  if(adef->mode == SEQUENCE_SIMILARITY_FILTER)
    {
      if(processID == 0)
	{
	  if(adef->sequenceSimilarity <= 0.0 || adef->sequenceSimilarity >= 1.0)
	    {
	      printf("\n ERROR: sequence similarity must be > 0.0 and < 1.0, exiting ...\n");
	      errorExit(-1);
	    }
	}
    }

  if(adef->mode == PER_SITE_LL)
    {                      
      if(!isGamma(adef))
	{
	  if(processID == 0)		
	    printf("\n ERROR: Computation of per-site log LHs is only allowed under GAMMA model of rate heterogeneity!\n");
	  errorExit(-1);	    
	}
	  
      if(!treesSet)
	{
	  if(processID == 0)		
	    printf("\n ERROR: For Computation of per-site log LHs you need to specify several input trees with \"-z\"\n");
	  errorExit(-1);	  
	}	
    }

  
    
  if(adef->mode == SPLIT_MULTI_GENE && (!adef->useMultipleModel))
    {
      if(processID == 0)
	{
	  printf("\n  Error, you are trying to split a multi-gene alignment into individual genes with the \"-f s\" option\n");	
	  printf("Without specifying a multiple model file with \"-q modelFileName\" \n");
	}
      errorExit(-1);
    }

  if(adef->mode == CALC_BIPARTITIONS && !treesSet)
    {
      if(processID == 0)
	printf("\n  Error, in bipartition computation mode you must specify a file containing multiple trees with the \"-z\" option\n");
      errorExit(-1);
    }

  if(adef->mode == CALC_BIPARTITIONS && !adef->restart)
    {
      if(processID == 0)
	printf("\n  Error, in bipartition computation mode you must specify a tree on which bipartition information will be drawn with the \"-t\" option\n");
      errorExit(-1);
    }

  if(!modelSet)
    { 
      if(processID == 0)
	printf("\n Error, you must specify a model of substitution with the \"-m\" option\n");
      errorExit(-1);
    }
      
  if(adef->computeDistance)
    {
      if(isCat(adef))
	{
	  if(processID == 0)
	    printf("\n Error pairwise distance computation only allowed for GAMMA-based models of rate heterogeneity\n");
	  errorExit(-1);
	}

      if(adef->restart)
	{
	  if(adef->randomStartingTree)
	    {
	      if(processID == 0)
		printf("\n Error pairwise distance computation not allowed for random starting trees\n");
	      errorExit(-1);
	    }
	  	       
	  if(adef->constraint)
	    {
	      if(processID == 0)
		printf("\n Error pairwise distance computation not allowed for binary backbone  constraint tree\n");
	      errorExit(-1);
	    }
	 
	  if(adef->grouping)
	    {
	      if(processID == 0)
		printf("\n Error pairwise distance computation not allowed for constraint tree\n");
	      errorExit(-1);
	    }
	      
	}

      if(adef->boot || adef->rapidBoot)
	{
	  if(processID == 0)
	    printf("\n Bootstrapping not implemented for pairwise distance computation\n");
	  errorExit(-1);
	}
    }

  


 

 

  if(!adef->restart && adef->mode == PARSIMONY_ADDITION)
    {
       if(processID == 0)
	 {
	   printf("\n You need to specify an incomplete binary input tree with \"-t\" to execute \n");
	   printf(" RAxML MP stepwise addition with \"-f p\"\n");
	 }
      errorExit(-1);
    }

  

  if(adef->restart && adef->randomStartingTree)
    {
      if(processID == 0)
	{
	  if(adef->constraint)
	    {
	      printf("\n Error you specified a binary constraint tree with -r AND the computation\n");
	      printf("of a random starting tree with -d for the same run\n");
	    }
	  else
	    {
	      if(adef->grouping)
		{
		  printf("\n Error you specified a multifurcating constraint tree with -g AND the computation\n");
		  printf("of a random starting tree with -d for the same run\n");
		}
	      else
		{
		  printf("\n Error you specified a starting tree with -t AND the computation\n");
		  printf("of a random starting tree with -d for the same run\n");
		}
	    }
	}
      errorExit(-1);
    }

  if(treeSet && constraintSet)
    {
      if(processID == 0)
	printf("\n Error you specified a binary constraint tree AND a starting tree for the same run\n");
      errorExit(-1);
    }


  if(treeSet && groupSet)
    {
      if(processID == 0)
	printf("\n Error you specified a multifurcating constraint tree AND a starting tree for the same run\n");
      errorExit(-1);
    }


  if(groupSet && constraintSet)
    {
      if(processID == 0)
	printf("\n Error you specified a bifurcating constraint tree AND a multifurcating constraint tree for the same run\n");
      errorExit(-1);
    }

  if(adef->restart && adef->startingTreeOnly)
    {
      if(processID == 0)
	{
	  printf("\n Error conflicting options: you want to compute only a parsimony starting tree with -y\n");
	  printf(" while you actually specified a starting tree with -t %s\n", tree_file);
	}
      errorExit(-1);
    }
  
  if((adef->mode == TREE_EVALUATION) && (!adef->restart))
    {
      if(processID == 0)
	printf("\n Error: please specify a treefile for the tree you want to evaluate with -t\n");
      errorExit(-1);
    }

#ifdef PARALLEL

  if(adef->mode == SPLIT_MULTI_GENE)
    {
      if(processID == 0)
	printf("Multi gene alignment splitting (-f s) not implemented for the MPI-Version\n");
      errorExit(-1);
    }

  if(adef->mode == TREE_EVALUATION)
    {
      if(processID == 0)
	printf("Tree Evaluation mode (-f e) not implemented for the MPI-Version\n");
      errorExit(-1);
    } 

 


   if(adef->mode == CALC_BIPARTITIONS)
     {
       if(processID == 0)
	 printf("Computation of bipartitions (-f b) not implemented for the MPI-Version\n");
       errorExit(-1);
     }
   
   if(adef->multipleRuns == 1)
     {
       if(processID == 0)
	 {
	   printf("Error: you are running the parallel MPI program but only want to compute one tree\n");
	   printf("For the MPI version you must specify a number of trees greater than 1 with the -# or -N option\n");
	 }
       errorExit(-1);
     }
#endif

   
   
   if((adef->mode == TREE_EVALUATION) && (isCat(adef)))
     {
       if(processID == 0)
	 {
	   printf("\n Error: No tree evaluation with GTRCAT/PROTCAT possible\n");
	   printf("the GTRCAT likelihood values are instable at present and should not\n");
	   printf("be used to compare trees based on ML values\n");
	 }
       errorExit(-1);
     }
  
  if(!nameSet)
    {
      if(processID == 0)
	printf("\n Error: please specify a name for this run with -n\n");
      errorExit(-1);
    }
    
  if(! alignmentSet && !adef->mrpEncoder)
    {
      if(processID == 0)
	printf("\n Error: please specify an alignment for this run with -s\n");
      errorExit(-1);
    }


  /*if(!workDirSet)*/
    {
#ifdef WIN32
      if(workdir[0]==0 || workdir[0] != '\\') 
	{
	  getcwd(buf,sizeof(buf));
	  if( buf[strlen(buf)-1] != '\\') strcat(buf,"\\");
	  strcat(buf,workdir);
	  if( buf[strlen(buf)-1] != '\\') strcat(buf,"\\");
	  strcpy(workdir,buf);
	}
#else
      if(workdir[0]==0 || workdir[0] != '/') 
	{	 
	  getcwd(buf,sizeof(buf));	  	 	  
	  if( buf[strlen(buf)-1] != '/') strcat(buf,"/");
	  strcat(buf,workdir);
	  if( buf[strlen(buf)-1] != '/') strcat(buf,"/");
	  strcpy(workdir,buf);	  
	}
#endif
    }
 
 
  return;
}




void errorExit(int e)
{
#ifdef PARALLEL
  MPI_Status msgStatus; 
  int i, dummy;

  if(processID == 0)
    {      
      for(i = 1; i < numOfWorkers; i++)	
	MPI_Send(&dummy, 1, MPI_INT, i, FINALIZE, MPI_COMM_WORLD);
  
      MPI_Finalize();
      exit(e);
    }     
  else
    {	
      MPI_Recv(&dummy, 1, MPI_INT, 0, FINALIZE, MPI_COMM_WORLD, &msgStatus);     
      MPI_Finalize();
      exit(e);
    }
#else 
  exit(e);
#endif
}



static void makeFileNames(void)
{
  int infoFileExists = 0;
#ifdef PARALLEL
  MPI_Status msgStatus; 
#endif

  strcpy(permFileName,         workdir);    
  strcpy(resultFileName,       workdir);
  strcpy(logFileName,          workdir);
  strcpy(checkpointFileName,   workdir);
  strcpy(infoFileName,         workdir);
  strcpy(randomFileName,       workdir);  
  strcpy(bootstrapFileName,    workdir);
  strcpy(bipartitionsFileName, workdir);
  strcpy(ratesFileName,        workdir);
  strcpy(lengthFileName,       workdir);
  strcpy(lengthFileNameModel,  workdir);
  strcpy( perSiteLLsFileName,  workdir);

  strcat(permFileName,         "RAxML_parsimonyTree.");
  strcat(resultFileName,       "RAxML_result.");
  strcat(logFileName,          "RAxML_log.");
  strcat(checkpointFileName,   "RAxML_checkpoint.");
  strcat(infoFileName,         "RAxML_info.");
  strcat(randomFileName,       "RAxML_randomTree."); 
  strcat(bootstrapFileName,    "RAxML_bootstrap.");  
  strcat(bipartitionsFileName, "RAxML_bipartitions.");
  strcat(ratesFileName,        "RAxML_perSiteRates.");
  strcat(lengthFileName,       "RAxML_treeLength.");
  strcat(lengthFileNameModel,  "RAxML_treeLengthModel.");
  strcat(perSiteLLsFileName,   "RAxML_perSiteLLs.");

  strcat(permFileName,         run_id);
  strcat(resultFileName,       run_id);
  strcat(logFileName,          run_id);
  strcat(checkpointFileName,   run_id);
  strcat(infoFileName,         run_id);    
  strcat(randomFileName,       run_id);   
  strcat(bootstrapFileName,    run_id); 
  strcat(bipartitionsFileName, run_id);
  strcat(ratesFileName,        run_id);
  strcat(lengthFileName,       run_id);
  strcat(lengthFileNameModel,  run_id);
  strcat(perSiteLLsFileName,   run_id);

  if(processID == 0)
    {
      infoFileExists = filexists(infoFileName);

#ifdef PARALLEL
      {
	int i;

	for(i = 1; i < numOfWorkers; i++)	
	  MPI_Send(&infoFileExists, 1, MPI_INT, i, FINALIZE, MPI_COMM_WORLD);
      }
#endif

      if(infoFileExists)
	{
	  printf("RAxML output files with the run ID <%s> already exist \n", run_id);
	  printf("in directory %s ...... exiting\n", workdir);
#ifdef PARALLEL 	  
	  MPI_Finalize();
	  exit(-1);
#else
	  exit(-1);
#endif
	}     
    }
#ifdef PARALLEL
  else
    {	
      MPI_Recv(&infoFileExists, 1, MPI_INT, 0, FINALIZE, MPI_COMM_WORLD, &msgStatus);
      if(infoFileExists)
	{	 	  
	  MPI_Finalize();
	  exit(-1);
	}    
    }
#endif
}


static void readData(analdef *adef, rawdata *rdta, cruncheddata *cdta, tree *tr)
{
  INFILE = myfopen(seq_file, "r");
   
  getinput(adef, rdta, cdta, tr); 
  
  fclose(INFILE);   
}



/***********************reading and initializing input ******************/


/********************PRINTING various INFO **************************************/


static void printModelAndProgramInfo(tree *tr, analdef *adef, int argc, char *argv[])
{ 
  if(processID == 0)
    {      
      int i, model;     
      FILE *infoFile = myfopen(infoFileName, "a");  
      char modelType[128];
      const char *secondaryModelList[21] = { "S6A (GTR)", "S6B", "S6C", "S6D", "S6E", "S7A (GTR)", "S7B", "S7C", "S7D", "S7E", "S7F", "S16 (GTR)", "S16A", "S16B", "S16C", 
					     "S16D", "S16E", "S16F", "S16I", "S16J", "S16K"};

      if(adef->useInvariant)
	strcpy(modelType, "GAMMA+P-Invar");
      else
	strcpy(modelType, "GAMMA");     
      
      printBoth(infoFile, "\n\nYou are using %s version %s released by Alexandros Stamatakis in %s\n",  programName, programVersion, programDate);     

      if(!adef->compressPatterns)
	printBoth(infoFile, "\nAlignment has %d columns\n\n",  tr->cdta->endsite);
      else
	printBoth(infoFile, "\nAlignment has %d distinct alignment patterns\n\n",  tr->cdta->endsite);
      

      if(adef->useInvariant)
	printBoth(infoFile, "Found %d invariant alignment patterns that correspond to %d columns \n", tr->numberOfInvariableColumns, tr->weightOfInvariableColumns);

      printBoth(infoFile, "Proportion of gaps and completely undetermined characters in this alignment: %f\n", adef->gapyness);
     
      switch(adef->mode)
	{       
	case THOROUGH_PARSIMONY:
	  printBoth(infoFile, "\nRAxML more exhaustive parsimony search with a ratchet.\n");
	  printBoth(infoFile, "For a faster and better implementation of MP searches please use TNT by Pablo Goloboff.\n\n");
	  break;
	case DISTANCE_MODE:
	  printBoth(infoFile, "\nRAxML Computation of pairwise distances\n\n");	 
	  break;
	case TREE_EVALUATION :   
	  printBoth(infoFile, "\nRAxML Model Optimization up to an accuracy of %f log likelihood units\n\n", adef->likelihoodEpsilon);	  
	  break;    
	case OLAF_OPTION:
	  printBoth(infoFile, "\nRAxML automatic dataset partitioning based on discrete evolutionary rates\n\n", adef->likelihoodEpsilon);	  
	  break;
	case  BIG_RAPID_MODE:
	  if(adef->rapidBoot)
	    {
	      if(adef->allInOne)
		printBoth(infoFile, "\nRAxML rapid bootstrapping and subsequent ML search\n\n");
	      else
		printBoth(infoFile,  "\nRAxML rapid bootstrapping algorithm\n\n");	
	    }
	  else	    
	    printBoth(infoFile, "\nRAxML rapid hill-climbing mode\n\n");	      	    
	  break;	 
	case CALC_BIPARTITIONS:
	  printBoth(infoFile, "\nRAxML Bipartition Computation: Drawing support values from trees in file %s onto tree in file %s\n\n", 
		    bootStrapFile, tree_file);	   
	  break;
	case PER_SITE_LL:  
	  printBoth(infoFile, "\nRAxML computation of per-site log likelihoods\n");
	  printBoth(infoFile, "\nIMPORTANT WARNING: ORDER OF LIKELIHOODS DOES NOT CORRESPOND TO COLUMN ORDER IN INPUT ALIGNEMNT!\n");
	  break;
	case PARSIMONY_ADDITION:
	  printBoth(infoFile, "\nRAxML stepwise MP addition to incomplete starting tree\n\n");	  	  	    
	  break;
	case CLASSIFY_ML:
	  printBoth(infoFile, "\nRAxML classification algorithm\n\n");	 	  
	  break;
	case GENERATE_BS:
	  printBoth(infoFile, "\nRAxML BS replicate generation\n\n");	  	 
	  break;
	case COMPUTE_ELW:
	  printBoth(infoFile, "\nRAxML ELW test\n\n");	 	  
	  break;
	case BOOTSTOP_ONLY:
	  printBoth(infoFile, "\nRAxML a posteriori Bootstrap convergence assessment\n\n");	 	  
	  break;
	case SUPER_FAST:
	  printBoth(infoFile, "\nRAxML experimental super fast\n\n");	 	 
	  break;
	case COMPUTE_LHS:
	  printBoth(infoFile, "\nRAxML computation of likelihoods for a set of trees\n\n");	  	 
	  break;
	case COMPUTE_BIPARTITION_CORRELATION:
	  printBoth(infoFile, "\nRAxML computation of bipartition support correlation on two sets of trees\n\n");	 
	  break;
	case COMPUTE_RF_DISTANCE:
	  printBoth(infoFile, "\nRAxML computation of RF distances for all pairs of trees in a set of trees\n\n");	 
	  break;
	default:
	  assert(0);
	}
                      
      if(adef->mode != THOROUGH_PARSIMONY)
	{
	  if(adef->perGeneBranchLengths)
	    printBoth(infoFile, "Using %d distinct models/data partitions with individual per partition branch length optimization\n\n\n", tr->NumberOfModels);
	  else
	    printBoth(infoFile, "Using %d distinct models/data partitions with joint branch length optimization\n\n\n", tr->NumberOfModels);	             
	}

      if(adef->mode == BIG_RAPID_MODE)
	{
	  if(adef->rapidBoot)	    
	    {
	      if(adef->allInOne)
		printBoth(infoFile, "\nExecuting %d rapid bootstrap inferences and thereafter a thorough ML search \n\n", adef->multipleRuns);	   
	      else
		printBoth(infoFile, "\nExecuting %d rapid bootstrap inferences\n\n", adef->multipleRuns);	   	    
	    }
	  else
	    {
	      if(adef->boot)			 				 
		printBoth(infoFile, "Executing %d non-parametric bootstrap inferences\n\n", adef->multipleRuns);
	      else		   
		{
		  char treeType[1024];
	      
		  if(adef->restart)
		    strcpy(treeType, "user-specifed");
		  else
		    {
		      if(adef->randomStartingTree)		
			strcpy(treeType, "distinct complete random");
		      else
			strcpy(treeType, "distinct randomized MP");
		    }
		  
		  printBoth(infoFile, "Executing %d inferences on the original alignment using %d %s trees\n\n", 
			    adef->multipleRuns, adef->multipleRuns, treeType);		  		
		}
	    }
	}
      
      if(adef->mode != THOROUGH_PARSIMONY)
	printBoth(infoFile, "All free model parameters will be estimated by RAxML\n"); 
		
      if(adef->mode != THOROUGH_PARSIMONY)
	{
	  if(tr->rateHetModel == GAMMA || tr->rateHetModel == GAMMA_I)	
	    printBoth(infoFile, "%s model of rate heteorgeneity, ML estimate of alpha-parameter\n\n", modelType);	  	   	
	  else
	    {
	      printBoth(infoFile, "ML estimate of %d per site rate categories\n\n", adef->categories);
	      if(adef->mode != CLASSIFY_ML)
		printBoth(infoFile, "Likelihood of final tree will be evaluated and optimized under %s\n\n", modelType);	 	 	  	  	 	 	   
	    }
	  
	  if(adef->mode != CLASSIFY_ML)
	    printBoth(infoFile, "%s Model parameters will be estimated up to an accuracy of %2.10f Log Likelihood units\n\n", 
		      modelType, adef->likelihoodEpsilon);
	}



      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  printBoth(infoFile, "Partition: %d\n", model);
	  printBoth(infoFile, "Alignment Patterns: %d\n", tr->partitionData[model].upper - tr->partitionData[model].lower);
	  printBoth(infoFile, "Name: %s\n", tr->partitionData[model].partitionName);	  	 

	  switch(tr->partitionData[model].dataType)
	    {
	    case DNA_DATA:
	      printBoth(infoFile, "DataType: DNA\n"); 
	      if(adef->mode != THOROUGH_PARSIMONY)
		printBoth(infoFile, "Substitution Matrix: GTR\n");	         	     	      	      	      	     
	      break;
	    case AA_DATA:	     
	      assert(tr->partitionData[model].protModels >= 0 && tr->partitionData[model].protModels < NUM_PROT_MODELS);		
	      printBoth(infoFile, "DataType: AA\n");
	      if(adef->mode != THOROUGH_PARSIMONY)
		{
		  printBoth(infoFile, "Substitution Matrix: %s\n", (adef->userProteinModel)?"External user-specified model":protModels[tr->partitionData[model].protModels]);	      
		  printBoth(infoFile, "%s Base Frequencies:\n", (tr->partitionData[model].protFreqs == 1)?"Empirical":"Fixed");
		}
	      break;
	    case BINARY_DATA:
	      printBoth(infoFile, "DataType: BINARY/MORPHOLOGICAL\n");
	      if(adef->mode != THOROUGH_PARSIMONY)
		printBoth(infoFile, "Substitution Matrix: Uncorrected\n");	     	     	     	      	      	      	     
	      break;
	    case SECONDARY_DATA:
	      printBoth(infoFile, "DataType: SECONDARY STRUCTURE\n");
	      if(adef->mode != THOROUGH_PARSIMONY)
		printBoth(infoFile, "Substitution Matrix: %s\n", secondaryModelList[tr->secondaryStructureModel]);	     	     	     	      	     	      	     
	      break;
	    case SECONDARY_DATA_6:
	      printBoth(infoFile, "DataType: SECONDARY STRUCTURE 6 STATE\n");
	      if(adef->mode != THOROUGH_PARSIMONY)
		printBoth(infoFile, "Substitution Matrix: %s\n", secondaryModelList[tr->secondaryStructureModel]);	     	     	     	      	     	      	     
	      break;
	    case SECONDARY_DATA_7:
	      printBoth(infoFile, "DataType: SECONDARY STRUCTURE 7 STATE\n");
	      if(adef->mode != THOROUGH_PARSIMONY)
		printBoth(infoFile, "Substitution Matrix: %s\n", secondaryModelList[tr->secondaryStructureModel]);	     	     	     	      	     	      	     
	      break;
	    default:
	      assert(0);
	    }
	  printBoth(infoFile, "\n\n\n");	  
	}            
      printBoth(infoFile, "\n");
      
      
      printBoth(infoFile, "RAxML was called as follows:\n\n");
      for(i = 0; i < argc; i++)
	printBoth(infoFile,"%s ", argv[i]);
      printBoth(infoFile,"\n\n\n");  
      
      fclose(infoFile);
    }
}

void printResult(tree *tr, analdef *adef, boolean finalPrint)
{
  FILE *logFile;
  char temporaryFileName[1024] = "", treeID[64] = "";

  strcpy(temporaryFileName, resultFileName);

  switch(adef->mode)
    {         
    case TREE_EVALUATION:
      

      Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE, finalPrint, adef, SUMMARIZE_LH);
            
      logFile = myfopen(temporaryFileName, "w");
      fprintf(logFile, "%s", tr->tree_string);
      fclose(logFile);          

      if(adef->perGeneBranchLengths)
	printTreePerGene(tr, adef, temporaryFileName, "w");


      break;     
    case BIG_RAPID_MODE:
      if(!adef->boot)
	{
	  if(adef->multipleRuns > 1)
	    {	  	 	  
	      sprintf(treeID, "%d", tr->treeID);	  	  
	      strcat(temporaryFileName, ".RUN.");
	      strcat(temporaryFileName, treeID);	  	 	      	
	    }
	  
	  
	  if(finalPrint)
	    {
	      switch(tr->rateHetModel)
		{
		case GAMMA:
		case GAMMA_I:
		  Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE, finalPrint, adef, 
			      SUMMARIZE_LH);
	      
		  logFile = myfopen(temporaryFileName, "w");
		  fprintf(logFile, "%s", tr->tree_string);
		  fclose(logFile);
	      
		  if(adef->perGeneBranchLengths)
		    printTreePerGene(tr, adef, temporaryFileName, "w");
		  break;
		case CAT:
		  Tree2String(tr->tree_string, tr, tr->start->back, FALSE, TRUE, FALSE, FALSE, finalPrint, adef, 
			      NO_BRANCHES);
	      
		  logFile = myfopen(temporaryFileName, "w");
		  fprintf(logFile, "%s", tr->tree_string);
		  fclose(logFile);
	      		 
		  break;
		default:
		  assert(0);
		}
	    }
	  else
	    {
	      Tree2String(tr->tree_string, tr, tr->start->back, FALSE, TRUE, FALSE, FALSE, finalPrint, adef, 
			  NO_BRANCHES);
	      logFile = myfopen(temporaryFileName, "w");
	      fprintf(logFile, "%s", tr->tree_string);
	      fclose(logFile);
	    }		 
	}
      break;       
    default:
      printf("FATAL ERROR call to printResult from undefined STATE %d\n", adef->mode);
      exit(-1);
      break;
    }
}

void printBootstrapResult(tree *tr, analdef *adef, boolean finalPrint)
{
  if(processID == 0)
    {
      FILE *logFile;
      
      if(adef->mode == BIG_RAPID_MODE && (adef->boot || adef->rapidBoot))
	{
#ifndef PARALLEL
	  if(adef->bootstrapBranchLengths)	    
	    {	      
	      Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE, finalPrint, adef, SUMMARIZE_LH);	    
	      logFile = myfopen(bootstrapFileName, "a");
	      fprintf(logFile, "%s", tr->tree_string);
	      fclose(logFile); 
	      if(adef->perGeneBranchLengths)
		printTreePerGene(tr, adef, bootstrapFileName, "a");
	    }
	  else
	    {	      	     
	      Tree2String(tr->tree_string, tr, tr->start->back, FALSE, TRUE, FALSE, FALSE, finalPrint, adef, NO_BRANCHES);
	      logFile = myfopen(bootstrapFileName, "a");
	      fprintf(logFile, "%s", tr->tree_string);
	      fclose(logFile); 
	    }
#else	  
	  logFile = myfopen(bootstrapFileName, "a");
	  fprintf(logFile, "%s", tr->tree_string);
	  fclose(logFile);     
#endif	  
	}
      else
	{
	  printf("FATAL ERROR in  printBootstrapResult\n");
	  exit(-1);	 
	}
    }
}



void printBipartitionResult(tree *tr, analdef *adef, boolean finalPrint)
{
  if(processID == 0 || adef->allInOne)
    {
      FILE *logFile;
      
      Tree2String(tr->tree_string, tr, tr->start->back, FALSE, TRUE, FALSE, TRUE, finalPrint, adef, NO_BRANCHES);
      logFile = myfopen(bipartitionsFileName, "a");
      fprintf(logFile, "%s", tr->tree_string);
      fclose(logFile);   
    }
}



void printLog(tree *tr, analdef *adef, boolean finalPrint)
{
  FILE *logFile;
  char temporaryFileName[1024] = "", checkPoints[1024] = "", treeID[64] = "";
  double lh, t;
  
  lh = tr->likelihood;
  t = gettime() - masterTime;

  strcpy(temporaryFileName, logFileName);
  strcpy(checkPoints,       checkpointFileName);

  switch(adef->mode)
    {    
    case TREE_EVALUATION:	 
      logFile = myfopen(temporaryFileName, "a");

      printf("%f %f\n", t, lh);
      fprintf(logFile, "%f %f\n", t, lh);

      fclose(logFile);     
      break;      
    case BIG_RAPID_MODE:
      if(adef->boot || adef->rapidBoot)
	{
	  /* testing only printf("%f %f\n", t, lh);*/
	  /* NOTHING PRINTED so far */
	}
      else
	{
	  if(adef->multipleRuns > 1)
	    {	  	 	  
	      sprintf(treeID, "%d", tr->treeID);	  	  
	      strcat(temporaryFileName, ".RUN.");
	      strcat(temporaryFileName, treeID);	  	 
	      
	      strcat(checkPoints, ".RUN.");
	      strcat(checkPoints, treeID);	      	      
	    }


	  if(!adef->checkpoints)
	    {
	      logFile = myfopen(temporaryFileName, "a");
#ifndef PARALLEL	      
	      /*printf("%f %1.20f\n", t, lh);*/
#endif
	      fprintf(logFile, "%f %f\n", t, lh);
	      
	      fclose(logFile);
	    }
	  else
	    {
	      logFile = myfopen(temporaryFileName, "a");
#ifndef PARALLEL	      
	      /*printf("%f %f %d\n", t, lh, tr->checkPointCounter);*/
#endif
	      fprintf(logFile, "%f %f %d\n", t, lh, tr->checkPointCounter);
	      
	      fclose(logFile);
	      
	      strcat(checkPoints, ".");

	      sprintf(treeID, "%d", tr->checkPointCounter);
	      strcat(checkPoints, treeID);
	     
	      Tree2String(tr->tree_string, tr, tr->start->back, FALSE, TRUE, FALSE, FALSE, finalPrint, adef, NO_BRANCHES);

	      logFile = myfopen(checkPoints, "a");
	      fprintf(logFile, "%s", tr->tree_string);
	      fclose(logFile);

	      tr->checkPointCounter++;
	    }
	}
      break;       
    default:
      assert(0);
    }
}



void printStartingTree(tree *tr, analdef *adef, boolean finalPrint)
{  
  if(adef->boot)
    {          
      /* not printing starting trees for bootstrap */
    }
  else
    {
      FILE *treeFile;
      char temporaryFileName[1024] = "", treeID[64] = "";
   
      Tree2String(tr->tree_string, tr, tr->start->back, FALSE, TRUE, FALSE, FALSE, finalPrint, adef, NO_BRANCHES);
          
      if(adef->randomStartingTree)	    
	strcpy(temporaryFileName, randomFileName);	    
      else
	strcpy(temporaryFileName, permFileName);

      if(adef->multipleRuns > 1)
	{	  	 	  
	  sprintf(treeID, "%d", tr->treeID);	  	  
	  strcat(temporaryFileName, ".RUN.");
	  strcat(temporaryFileName, treeID);	  	 
	}
     	  	 
      treeFile = myfopen(temporaryFileName, "a");	
      fprintf(treeFile, "%s", tr->tree_string);
      fclose(treeFile);	    	
    }
}

void writeInfoFile(analdef *adef, tree *tr, double t)
{  
  if(processID == 0)
    {
      FILE *infoFile = myfopen(infoFileName, "a");           

      switch(adef->mode)
	{
	case TREE_EVALUATION:	 	        
	  break;      
	case BIG_RAPID_MODE:
	  if(adef->boot || adef->rapidBoot)
	    {
	      if(!adef->initialSet)	       
		{
		  fprintf(infoFile, "Bootstrap[%d]: Time %f bootstrap likelihood %f, best rearrangement setting %d\n", tr->treeID, t, tr->likelihood,  adef->bestTrav);		
		  printf("Bootstrap[%d]: Time %f bootstrap likelihood %f, best rearrangement setting %d\n", tr->treeID, t, tr->likelihood,  adef->bestTrav);
		}
	      else		
		{
		  fprintf(infoFile, "Bootstrap[%d]: Time %f bootstrap likelihood %f\n", tr->treeID, t, tr->likelihood);	
		  printf("Bootstrap[%d]: Time %f bootstrap likelihood %f\n", tr->treeID, t, tr->likelihood);	
		}
	    }
	  else
	    {
	      	
	      int model;
	      char modelType[128];
	      
	      switch(tr->rateHetModel)
		{
		case GAMMA_I:
		  strcpy(modelType, "GAMMA+P-Invar");
		  break;
		case GAMMA:	   
		  strcpy(modelType, "GAMMA");		  	       
		  break;
		case CAT:
		  strcpy(modelType, "CAT");
		  break;
		default:
		  assert(0);
		}
	      
	      if(!adef->initialSet)
		{
		  printf("Inference[%d]: Time %f %s-based likelihood %f, best rearrangement setting %d\n", 
			 tr->treeID, t, modelType, tr->likelihood,  adef->bestTrav);
		  fprintf(infoFile, "Inference[%d]: Time %f %s-based likelihood %f, best rearrangement setting %d, ", 
			  tr->treeID, t, modelType, tr->likelihood,  adef->bestTrav);
		}
	      else
		{
		  printf("Inference[%d]: Time %f %s-based likelihood %f\n", 
			 tr->treeID, t, modelType, tr->likelihood);
		  fprintf(infoFile, "Inference[%d]: Time %f %s-based likelihood %f, ", 
			  tr->treeID, t, modelType, tr->likelihood);		    		  
		}
	      
	      for(model = 0; model < tr->NumberOfModels; model++)		    		    
		{
		  fprintf(infoFile, "alpha[%d]: %f ", model, tr->partitionData[model].alpha);
		  if(adef->useInvariant)
		    fprintf(infoFile, "invar[%d]: %f ", model, tr->partitionData[model].propInvariant);
#ifndef PARALLEL
		  if(tr->partitionData[model].dataType == DNA_DATA)
		    {
		      int k;
		      
		      fprintf(infoFile, "rates[%d] ac ag at cg ct gt: ",model);
		      for(k = 0; k < DNA_RATES; k++)			    
			fprintf(infoFile, "%f ", tr->partitionData[model].substRates[k]);			    
		    }
		  fprintf(infoFile, "1.0 ");
#endif
		}
	      
	      fprintf(infoFile, "\n");
	    }	
	  break;
	default:
	  assert(0);
	}

      fclose(infoFile);
    }
}

static void finalizeInfoFile(tree *tr, analdef *adef)
{
  if(processID == 0)
    {
      FILE *infoFile = myfopen(infoFileName, "a");
      double t;
      int model;

      t = gettime() - masterTime;

      switch(adef->mode)
	{       
	case TREE_EVALUATION :
	  printf("\n\nOverall Time for Tree Evaluation %f\n", t);	 
	  printf("Final GAMMA  likelihood: %f\n", tr->likelihood);

	  fprintf(infoFile, "\n\nOverall Time for Tree Evaluation %f\n", t);       
	  fprintf(infoFile, "Final GAMMA  likelihood: %f\n", tr->likelihood);

	  {
	    int 
	      params,
	      paramsBrLen;

	    if(tr->NumberOfModels == 1)
	      {
		if(adef->useInvariant)
		  {
		    params      = 1 /* INVAR */ + 5 /* RATES */ + 3 /* freqs */ + 1 /* alpha */;
		    paramsBrLen = 1 /* INVAR */ + 5 /* RATES */ + 3 /* freqs */ + 1 /* alpha */ + 
		      (2 * tr->mxtips - 3);
		  }
		else
		  {
		    params      = 5 /* RATES */ + 3 /* freqs */ + 1 /* alpha */;
		    paramsBrLen = 5 /* RATES */ + 3 /* freqs */ + 1 /* alpha */ + 
		      (2 * tr->mxtips - 3);
		  }
	      }
	    else
	      {
		if(tr->multiBranch)
		  {
		    if(adef->useInvariant)
		      {
			params      = tr->NumberOfModels * (1 /* INVAR */ + 5 /* RATES */ + 3 /* freqs */ + 1 /* alpha */);
			paramsBrLen = tr->NumberOfModels * (1 /* INVAR */ + 5 /* RATES */ + 3 /* freqs */ + 1 /* alpha */ + 
							    (2 * tr->mxtips - 3));
		      }
		    else
		      {
			params      = tr->NumberOfModels * (5 /* RATES */ + 3 /* freqs */ + 1 /* alpha */);
			paramsBrLen = tr->NumberOfModels * (5 /* RATES */ + 3 /* freqs */ + 1 /* alpha */ + 
							    (2 * tr->mxtips - 3));
		      }
		  }
		else
		  {
		    if(adef->useInvariant)
		      {
			params      = tr->NumberOfModels * (1 /* INVAR */ + 5 /* RATES */ + 3 /* freqs */ + 1 /* alpha */);
			paramsBrLen = tr->NumberOfModels * (1 /* INVAR */ + 5 /* RATES */ + 3 /* freqs */ + 1 /* alpha */) 
			  + (2 * tr->mxtips - 3);
		      }
		    else
		      {
			params      = tr->NumberOfModels * (5 /* RATES */ + 3 /* freqs */ + 1 /* alpha */);
			paramsBrLen = tr->NumberOfModels * (5 /* RATES */ + 3 /* freqs */ + 1 /* alpha */) 
			  + (2 * tr->mxtips - 3);
		      }

		  }
	      }
		
	    if(tr->partitionData[0].dataType == DNA_DATA)
	      {
		printf("Number of free parameters for AIC-TEST(BR-LEN): %d\n",    paramsBrLen);
		printf("Number of free parameters for AIC-TEST(NO-BR-LEN): %d\n", params);
		fprintf(infoFile, "Number of free parameters for AIC-TEST(BR-LEN): %d\n",    paramsBrLen);
		fprintf(infoFile, "Number of free parameters for AIC-TEST(NO-BR-LEN): %d\n", params);
	      }
	    
	  }
	
	  printf("\n\n");
	  fprintf(infoFile, "\n\n");

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
		  strcpy(typeOfData,"SECONDARY STRUCTURE 16 STATE");
		  break;
		case SECONDARY_DATA_6:
		  strcpy(typeOfData,"SECONDARY STRUCTURE 6 STATE");
		  break; 
		case SECONDARY_DATA_7:
		  strcpy(typeOfData,"SECONDARY STRUCTURE 7 STATE");
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

	      if(tr->multiBranch)
		tl = treeLength(tr, model);
	      else
		tl = treeLength(tr, 0);
	      
	      fprintf(infoFile, "Tree-Length: %f\n", tl);    
	      printf("Tree-Length: %f\n", tl);       
	      
      

	      switch(tr->partitionData[model].dataType)
		{
		case SECONDARY_DATA:
		case SECONDARY_DATA_6:
		case SECONDARY_DATA_7: 
		case BINARY_DATA:	     
		case AA_DATA:
		  break;
		case DNA_DATA:
		  {
		    int k;
		    char *names[6] = {"a<->c", "a<->g", "a<->t", "c<->g", "c<->t", "g<->t"};	 
		    for(k = 0; k < DNA_RATES; k++)			    
		      {
			fprintf(infoFile, "rate %s: %f\n", names[k], tr->partitionData[model].substRates[k]);			    
			printf("rate %s: %f\n", names[k], tr->partitionData[model].substRates[k]);
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
	  	
	  printf("Final tree written to:                 %s\n", resultFileName);  
	  printf("Execution Log File written to:         %s\n", logFileName);
	   
	  fprintf(infoFile, "Final tree written to:                 %s\n", resultFileName);  
	  fprintf(infoFile, "Execution Log File written to:         %s\n", logFileName);	  
	
	  break;  
	case  BIG_RAPID_MODE:
	  if(adef->boot)
	    {
	      printf("\n\nOverall Time for %d Bootstraps %f\n", adef->multipleRuns, t);
	      printf("\n\nAverage Time per Bootstrap %f\n", (double)(t/((double)adef->multipleRuns)));
	      printf("All %d bootstrapped trees written to: %s\n", adef->multipleRuns, bootstrapFileName);
	      
	      fprintf(infoFile, "\n\nOverall Time for %d Bootstraps %f\n", adef->multipleRuns, t);
	      fprintf(infoFile, "Average Time per Bootstrap %f\n", (double)(t/((double)adef->multipleRuns)));	     
	      fprintf(infoFile, "\n\nAll %d bootstrapped trees written to: %s\n", adef->multipleRuns, bootstrapFileName);	     
	    }
	  else
	    {
	      if(adef->multipleRuns > 1)
		{
		  double avgLH = 0;
		  double bestLH = unlikely;
		  int i, bestI  = 0;
		  
		  for(i = 0; i < adef->multipleRuns; i++)
		    {     
		      avgLH   += tr->likelihoods[i];
		      if(tr->likelihoods[i] > bestLH)
			{
			  bestLH = tr->likelihoods[i];
			  bestI  = i;
			}
		    }
		  avgLH /= ((double)adef->multipleRuns);

		  printf("\n\nOverall Time for %d Inferences %f\n", adef->multipleRuns, t);
		  printf("Average Time per Inference %f\n", (double)(t/((double)adef->multipleRuns)));		 
		  printf("Average Likelihood   : %f\n", avgLH);
		  printf("\n");
		  printf("Best Likelihood in run number %d: likelihood %f\n\n", bestI, bestLH);

		  if(adef->checkpoints)   
		    printf("Checkpoints written to:                 %s.RUN.%d.* to %d.*\n", checkpointFileName, 0, adef->multipleRuns - 1);  
		  if(!adef->restart)
		    {
		      if(adef->randomStartingTree)
			printf("Random starting trees written to:       %s.RUN.%d to %d\n", randomFileName, 0, adef->multipleRuns - 1);
		      else
			printf("Parsimony starting trees written to:    %s.RUN.%d to %d\n", permFileName, 0, adef->multipleRuns - 1);   	            
		    }					  
		  printf("Final trees written to:                 %s.RUN.%d to %d\n", resultFileName,  0, adef->multipleRuns - 1);  
		  printf("Execution Log Files written to:         %s.RUN.%d to %d\n", logFileName, 0, adef->multipleRuns - 1);   
		  printf("Execution information file written to:  %s\n", infoFileName);
		  

		  fprintf(infoFile, "\n\nOverall Time for %d Inferences %f\n", adef->multipleRuns, t);
		  fprintf(infoFile, "Average Time per Inference %f\n", (double)(t/((double)adef->multipleRuns)));		 
		  fprintf(infoFile, "Average Likelihood   : %f\n", avgLH);
		  fprintf(infoFile, "\n");
		  fprintf(infoFile, "Best Likelihood in run number %d: likelihood %f\n\n", bestI, bestLH); 
		  if(adef->checkpoints)   
		    fprintf(infoFile, "Checkpoints written to:                %s.RUN.%d.* to %d.*\n", checkpointFileName, 0, adef->multipleRuns - 1);  
		  if(!adef->restart)
		    {
		      if(adef->randomStartingTree)
			fprintf(infoFile, "Random starting trees written to:      %s.RUN.%d to %d\n", randomFileName, 0, adef->multipleRuns - 1);
		      else
			fprintf(infoFile, "Parsimony starting trees written to:   %s.RUN.%d to %d\n", permFileName, 0, adef->multipleRuns - 1);   	            
		    }					  
		  fprintf(infoFile, "Final trees written to:                %s.RUN.%d to %d\n", resultFileName,  0, adef->multipleRuns - 1);  
		  fprintf(infoFile, "Execution Log Files written to:        %s.RUN.%d to %d\n", logFileName, 0, adef->multipleRuns - 1);   
		  fprintf(infoFile, "Execution information file written to: %s\n", infoFileName);
		  		  		 		   
		}
	      else
		{
		  printf("\n\nOverall Time for 1 Inference %f\n", t);		  
		  printf("Likelihood   : %f\n", tr->likelihood);
		  printf("\n\n");	     

		  if(adef->checkpoints)   
		  printf("Checkpoints written to:                %s.*\n", checkpointFileName);  
		  if(!adef->restart)
		    {
		      if(adef->randomStartingTree)
			printf("Random starting tree written to:       %s\n", randomFileName);
		      else
			printf("Parsimony starting tree written to:    %s\n", permFileName);   	            
		    }					  
		  printf("Final tree written to:                 %s\n", resultFileName);  
		  printf("Execution Log File written to:         %s\n", logFileName);   
		  printf("Execution information file written to: %s\n",infoFileName);

		  
		  
		  fprintf(infoFile, "\n\nOverall Time for 1 Inference %f\n", t);		  
		  fprintf(infoFile, "Likelihood   : %f\n", tr->likelihood);
		  fprintf(infoFile, "\n\n");

		  if(adef->checkpoints)   
		    fprintf(infoFile, "Checkpoints written to:                %s.*\n", checkpointFileName);  
		  if(!adef->restart)
		    {
		      if(adef->randomStartingTree)
			fprintf(infoFile, "Random starting tree written to:       %s\n", randomFileName);
		      else
			fprintf(infoFile, "Parsimony starting tree written to:    %s\n", permFileName);   	            
		    }					  
		  fprintf(infoFile, "Final tree written to:                 %s\n", resultFileName);  
		  fprintf(infoFile, "Execution Log File written to:         %s\n", logFileName);   
		  fprintf(infoFile, "Execution information file written to: %s\n",infoFileName);
		  		 		  
		}
	    }
	    
	  break;	 
	case CALC_BIPARTITIONS:
	  printf("\n\nTime for Computation of Bipartitions %f\n", t);
	  printf("Tree with bipartitions written to file:  %s\n", bipartitionsFileName);
	  printf("Execution information file written to :  %s\n",infoFileName);
	

	  fprintf(infoFile, "\n\nTime for Computation of Bipartitions %f\n", t);
	  fprintf(infoFile, "Tree with bipartitions written to file:  %s\n", bipartitionsFileName);

	 
	  break;
	case PER_SITE_LL:	     	   
	  printf("\n\nTime for Optimization of per-site log likelihoods %f\n", t);
	  printf("Per-site Log Likelihoods written to File %s in Tree-Puzzle format\n",  perSiteLLsFileName);
	  printf("Execution information file written to :  %s\n",infoFileName);
	  	       
	  fprintf(infoFile, "\n\nTime for Optimization of per-site log likelihoods %f\n", t);
	  fprintf(infoFile, "Per-site Log Likelihoods written to File %s in Tree-Puzzle format\n",  perSiteLLsFileName);	
	 
	  break;
	case PARSIMONY_ADDITION:
	  printf("\n\nTime for MP stepwise addition %f\n", t);	 
	  printf("Execution information file written to :  %s\n",infoFileName);
	  printf("Complete parsimony tree written to:      %s\n", permFileName); 

	 

	  fprintf(infoFile, "\n\nTime for MP stepwise addition %f\n", t);	 
	  fprintf(infoFile, "Complete parsimony tree written to:      %s\n", permFileName); 

	 
	  break;
	default:
	  assert(0);
	}
      fclose(infoFile);
    }

}



/********************PRINTING various INFO **************************************/

/************************************************************************************/

static void computeLHTest(tree *tr, analdef *adef, char *bootStrapFileName)
{
  int 
    numberOfTrees = 0, 
    i; 
  double bestLH, currentLH, weightSum = 0.0;
  double *bestVector = (double*)malloc(sizeof(double) * tr->cdta->endsite);

  for(i = 0; i < tr->cdta->endsite; i++)
    weightSum += (double)(tr->cdta->aliaswgt[i]);

  modOpt(tr, adef, TRUE, adef->likelihoodEpsilon);
  printf("Model optimization, best Tree: %f\n", tr->likelihood);
  bestLH = tr->likelihood;


  evaluateGenericInitrav(tr, tr->start); 

  evaluateGenericVector(tr, tr->start);
  memcpy(bestVector, tr->perSiteLL, tr->cdta->endsite * sizeof(double));
      
  INFILE = myfopen(bootStrapFileName, "r");       
  numberOfTrees = countTrees(INFILE);
 
 
  printf("Found %d trees in File %s\n", numberOfTrees, bootStrapFileName);
 
  for(i = 0; i < numberOfTrees; i++)
    {              
      treeReadLen(INFILE, tr, adef);      
      treeEvaluate(tr, 2);
      tr->start = tr->nodep[1];

      evaluateGenericInitrav(tr, tr->start);         

      currentLH = tr->likelihood;
      if(currentLH > bestLH)
	{
	  printf("Better tree found %d at %f\n", i, currentLH);
	  /*exit(1);*/
	}
      /*printf("Tree %d %f\n",i, tr->likelihood);*/
     
      evaluateGenericVector(tr, tr->start);
     
      {
	 int j;
	 double temp, wtemp, sum, sum2, sd;

	 sum = 0.0;
	 sum2 = 0.0;

	 for (j = 0; j < tr->cdta->endsite; j++) 
	   {
	     temp  = bestVector[j] - tr->perSiteLL[j];
	     wtemp = tr->cdta->aliaswgt[j] * temp;
	     sum  += wtemp;
	     sum2 += wtemp * temp;
	   }
	 
	 sd = sqrt( weightSum * (sum2 - sum*sum / weightSum)
		    / (weightSum - 1) );	       
	 /* this is for a 5% p level */
	 printf("Tree: %d Likelihood: %f D(LH): %f SD: %f Significantly Worse: %s\n", i, currentLH, currentLH - bestLH, sd, (sum > 1.95996 * sd) ? "Yes" : " No");	 
      }
    }
      
  fclose(INFILE); 
  
  free(bestVector);
  exit(0);
}

static void computePerSiteLLs(tree *tr, analdef *adef, char *bootStrapFileName)
{
  int 
    numberOfTrees = 0, 
    i, 
    j;
  FILE *tlf;	   
	  
  tlf = myfopen( perSiteLLsFileName, "w");
      
  INFILE = myfopen(bootStrapFileName, "r"); 
  numberOfTrees = countTrees(INFILE);
 
 
  printf("Found %d trees in File %s\n", numberOfTrees, bootStrapFileName);
 
  fprintf(tlf, "  %d  %d\n", numberOfTrees, tr->cdta->endsite);

  for(i = 0; i < numberOfTrees; i++)
    {        
      double check = 0.0;
      treeReadLen(INFILE, tr, adef);      
      if(i == 0)	
	modOpt(tr, adef, TRUE, adef->likelihoodEpsilon);	  	
      else
	treeEvaluate(tr, 2);

      printf("Tree %d: %f\n", i, tr->likelihood);

      tr->start = tr->nodep[1];

      evaluateGenericInitrav(tr, tr->start);      
     
      evaluateGenericVector(tr, tr->start);                   

      fprintf(tlf, "tr%d\t", i + 1);
      for(j = 0; j < tr->cdta->endsite; j++)	
	{
	  int w;
	  for(w = 0; w < tr->cdta->aliaswgt[j]; w++)
	    {
	      fprintf(tlf, "%f ", tr->perSiteLL[j]);	
	      check += tr->perSiteLL[j];
	    }
	}
      printf("check: %f\n", check);
      fprintf(tlf, "\n");     
    }
      
  fclose(INFILE); 
  fclose(tlf);  
}

#ifdef _USE_PTHREADS


typedef struct 
{
  tree *tr;
  int threadNumber;
} 
  threadData;



static void computeFraction(tree *localTree, int tid, int n)
{
  int 
    i,
    model;

  for(model = 0; model < localTree->NumberOfModels; model++)
    {    
      int width = 0;

      for(i = localTree->partitionData[model].lower; i < localTree->partitionData[model].upper; i++)		 
	if(i % n == tid)
	      width++;   
      
      localTree->partitionData[model].width = width;      
    }  
}



static void threadFixModelIndices(tree *tr, tree *localTree, int tid, int n)
{
  size_t 
    model, 
    j, 
    i,
    globalCounter = 0,
    localCounter  = 0,
    offset,
    countOffset,
    myLength = 0,
    memoryRequirements = 0;
      
  for(model = 0; model < (size_t)localTree->NumberOfModels; model++)
    {
      localTree->partitionData[model].lower      = tr->partitionData[model].lower;
      localTree->partitionData[model].upper      = tr->partitionData[model].upper; 
    }

  computeFraction(localTree, tid, n);
  
  for(model = 0, offset = 0, countOffset = 0; model < (size_t)localTree->NumberOfModels; model++)
    {      
      localTree->partitionData[model].sumBuffer    = &localTree->sumBuffer[offset];
      localTree->partitionData[model].perSiteLL    = &localTree->perSiteLLPtr[countOffset];
      localTree->partitionData[model].wr           = &localTree->wrPtr[countOffset];
      localTree->partitionData[model].wr2          = &localTree->wr2Ptr[countOffset];
      localTree->partitionData[model].wgt          = &localTree->wgtPtr[countOffset];
      localTree->partitionData[model].invariant    = &localTree->invariantPtr[countOffset];
      localTree->partitionData[model].rateCategory = &localTree->rateCategoryPtr[countOffset];
      
      /*assert(localTree->partitionData[model].width > 0);*/
      
      countOffset += localTree->partitionData[model].width;	 	 	 
      
      switch(localTree->partitionData[model].dataType)
	{
	case AA_DATA:
	  offset += ((size_t)localTree->aaIncrement * localTree->partitionData[model].width);
	  break;
	case DNA_DATA:
	  offset += ((size_t)localTree->dnaIncrement * localTree->partitionData[model].width);
	  break;
	case BINARY_DATA:
	  offset += ((size_t)localTree->binaryIncrement * localTree->partitionData[model].width);
	  break;       
	case SECONDARY_DATA:
	  offset += ((size_t)localTree->secondaryIncrement * localTree->partitionData[model].width);
	  break;
	case SECONDARY_DATA_6:
	  offset += ((size_t)localTree->secondaryIncrement6 * localTree->partitionData[model].width);
	  break; 
	case SECONDARY_DATA_7:
	  offset += ((size_t)localTree->secondaryIncrement7 * localTree->partitionData[model].width);
	  break;
	default:
	  assert(0);
	}	 	 	      
    }

  myLength           = countOffset;
  memoryRequirements = offset;

  
  /* figure in data */

  
  for(i = 0; i < (size_t)localTree->mxtips; i++)
    {    
      for(model = 0, offset = 0, countOffset = 0; model < (size_t)localTree->NumberOfModels; model++)
	{	
	  size_t width = localTree->partitionData[model].width;	       
	  
	  localTree->partitionData[model].expVector[i] = &localTree->expArray[i * myLength + countOffset];	 
	  localTree->partitionData[model].xVector[i]   = &localTree->likelihoodArray[i * memoryRequirements + offset];
	  localTree->partitionData[model].pVector[i]   = (parsimonyVector *)localTree->partitionData[model].xVector[i];
	  localTree->partitionData[model].yVector[i+1]   = &localTree->y_ptr[i * myLength + countOffset];
	  
	  countOffset += width;
	  
	  switch(localTree->partitionData[model].dataType)
	    {
	    case AA_DATA:
	      offset += ((size_t)localTree->aaIncrement * width);
	      break;
	    case DNA_DATA:
	      offset += ((size_t)localTree->dnaIncrement * width);
	      break;
	    case BINARY_DATA:
	      offset += ((size_t)localTree->binaryIncrement * width);
	      break;
	    case SECONDARY_DATA:
	      offset += ((size_t)localTree->secondaryIncrement * width);
	      break;
	    case SECONDARY_DATA_6:
	      offset += ((size_t)localTree->secondaryIncrement6 * width);
	      break; 
	    case SECONDARY_DATA_7:
	      offset += ((size_t)localTree->secondaryIncrement7 * width);
	      break;
	    default:
	      assert(0);
	    }
	}         
      assert(countOffset == myLength);     
    }

  for(model = 0, globalCounter = 0; model < (size_t)localTree->NumberOfModels; model++)
    {
      for(localCounter = 0, i = (size_t)localTree->partitionData[model].lower;  i < (size_t)localTree->partitionData[model].upper; i++)
	{
	  if(i % (size_t)n == (size_t)tid)
	    {	    
	      localTree->partitionData[model].wgt[localCounter]          = tr->cdta->aliaswgt[globalCounter];
	      localTree->partitionData[model].wr[localCounter]           = tr->cdta->wr[globalCounter];
	      localTree->partitionData[model].wr2[localCounter]          = tr->cdta->wr2[globalCounter];
	      localTree->partitionData[model].invariant[localCounter]    = tr->invariant[globalCounter];
	      localTree->partitionData[model].rateCategory[localCounter] = tr->cdta->rateCategory[globalCounter];
      	      
	      for(j = 1; j <= (size_t)localTree->mxtips; j++)		 
	       localTree->partitionData[model].yVector[j][localCounter] = tr->yVector[j][globalCounter];	

	      localCounter++;
	    }
	  globalCounter++;
	}
    }
  
}


static void initPartition(tree *tr, tree *localTree, int tid)
{
  int model;

  localTree->threadID = tid;

  if(tid > 0)
    {
      int totalLength = 0;      

      localTree->originalCrunchedLength  = tr->originalCrunchedLength;
      localTree->NumberOfModels          = tr->NumberOfModels;
      localTree->mxtips                  = tr->mxtips;    
      localTree->multiBranch             = tr->multiBranch;
      localTree->numBranches             = tr->numBranches; 
      localTree->lhs                     = (double*)malloc(sizeof(double)   * localTree->originalCrunchedLength);
      localTree->executeModel            = (boolean*)malloc(sizeof(boolean) * localTree->NumberOfModels);
      localTree->perPartitionLH          = (double*)malloc(sizeof(double)   * localTree->NumberOfModels);
      
      
      localTree->partitionData = (pInfo*)malloc(sizeof(pInfo) * localTree->NumberOfModels);  
      
      /* extend for multi-branch */
      localTree->td[0].count = 0;
      localTree->td[0].ti    = (traversalInfo *)malloc(sizeof(traversalInfo) * localTree->mxtips);

      localTree->cdta               = (cruncheddata*)malloc(sizeof(cruncheddata));
      localTree->cdta->patrat       = (double*)malloc(sizeof(double) * localTree->originalCrunchedLength);
      localTree->cdta->patratStored = (double*)malloc(sizeof(double) * localTree->originalCrunchedLength);

      localTree->NumberOfCategories = tr->NumberOfCategories; 
  
      localTree->secondaryIncrement6 =  tr->secondaryIncrement6;
      localTree->secondaryIncrement7 =  tr->secondaryIncrement7;
      localTree->secondaryIncrement =  tr->secondaryIncrement;
      localTree->binaryIncrement    = tr->binaryIncrement;
      localTree->dnaIncrement       = tr->dnaIncrement;
      localTree->aaIncrement        = tr->aaIncrement;

      for(model = 0; model < localTree->NumberOfModels; model++)
	{
	  localTree->partitionData[model].dataType   = tr->partitionData[model].dataType;
	  localTree->partitionData[model].protModels = tr->partitionData[model].protModels;
	  localTree->partitionData[model].protFreqs  = tr->partitionData[model].protFreqs;      
	  localTree->partitionData[model].mxtips     = tr->partitionData[model].mxtips;
	  localTree->partitionData[model].lower      = tr->partitionData[model].lower;
	  localTree->partitionData[model].upper      = tr->partitionData[model].upper;
	  localTree->executeModel[model]             = TRUE;
	  localTree->perPartitionLH[model]          = 0.0;
	  totalLength += (localTree->partitionData[model].upper -  localTree->partitionData[model].lower);
	}     
          
      assert(totalLength == localTree->originalCrunchedLength);      
    }  

  for(model = 0; model < localTree->NumberOfModels; model++)       
    localTree->partitionData[model].width        = 0;         
}




static void allocNodex(tree *tr, int tid, int n)
{
  nodeptr  p;
  size_t  
    i,   
    model,     
    memoryRequirements = 0,    
    myLength = 0;
 
  computeFraction(tr, tid, n);

  for(i = 0; i < (size_t)tr->NumberOfModels; i++)
    {    
      switch(tr->partitionData[i].dataType)
	{
	case DNA_DATA:
	  tr->partitionData[i].EIGN = (double*)malloc(3 * sizeof(double));
	  tr->partitionData[i].EV   = (double*)malloc(16 * sizeof(double));
	  tr->partitionData[i].EI   = (double*)malloc(12 * sizeof(double));	      
	  tr->partitionData[i].substRates = (double *)malloc(5 * sizeof(double));	      
	  tr->partitionData[i].frequencies =  (double*)malloc(4 * sizeof(double));	      
	  tr->partitionData[i].tipVector   = (double *)malloc(64 * sizeof(double));
	  tr->partitionData[i].symmetryVector = (int *)malloc(6  * sizeof(int));
	  tr->partitionData[i].frequencyGrouping = (int *)malloc(4  * sizeof(int));
	  tr->partitionData[i].nonGTR = FALSE;
	  break;
	case AA_DATA:
	  tr->partitionData[i].EIGN = (double*)malloc(19 * sizeof(double));
	  tr->partitionData[i].EV   = (double*)malloc(400 * sizeof(double));
	  tr->partitionData[i].EI   = (double*)malloc(380 * sizeof(double));	      
	  tr->partitionData[i].tipVector   = (double *)malloc(460 * sizeof(double)); 	      
	  tr->partitionData[i].substRates = (double *)malloc(190 * sizeof(double));	      
	  tr->partitionData[i].frequencies = (double*)malloc(20 * sizeof(double));
	  tr->partitionData[i].symmetryVector = (int *)malloc(190  * sizeof(int));
	  tr->partitionData[i].frequencyGrouping = (int *)malloc(20  * sizeof(int));
	  tr->partitionData[i].nonGTR = FALSE;
	  break;
	case BINARY_DATA:
	  tr->partitionData[i].EIGN = (double*)malloc(1 * sizeof(double));
	  tr->partitionData[i].EV   = (double*)malloc(4 * sizeof(double));
	  tr->partitionData[i].EI   = (double*)malloc(2 * sizeof(double));	  
	  tr->partitionData[i].tipVector   = (double *)malloc(8 * sizeof(double)); 	  
	  tr->partitionData[i].substRates = (double *)malloc(1 * sizeof(double));	  
	  tr->partitionData[i].frequencies = (double*)malloc(2 * sizeof(double));
	  tr->partitionData[i].symmetryVector = (int *)malloc(2  * sizeof(int));
	  tr->partitionData[i].frequencyGrouping = (int *)malloc(2  * sizeof(int));
	  tr->partitionData[i].nonGTR = FALSE;
	  break;
	case SECONDARY_DATA:
	  tr->partitionData[i].EIGN = (double*)malloc(15 * sizeof(double));
	  tr->partitionData[i].EV   = (double*)malloc(256 * sizeof(double));
	  tr->partitionData[i].EI   = (double*)malloc(240 * sizeof(double));	  
	  tr->partitionData[i].tipVector   = (double *)malloc(4096 * sizeof(double)); 	  
	  tr->partitionData[i].substRates = (double *)malloc(120 * sizeof(double));	  
	  tr->partitionData[i].frequencies = (double*)malloc(16 * sizeof(double));
	  tr->partitionData[i].symmetryVector = (int *)malloc(120  * sizeof(int));
	  tr->partitionData[i].frequencyGrouping = (int *)malloc(16  * sizeof(int));
	  tr->partitionData[i].nonGTR = FALSE;
	  break; 
	case SECONDARY_DATA_6:
	  tr->partitionData[i].EIGN = (double*)malloc(5 * sizeof(double));
	  tr->partitionData[i].EV   = (double*)malloc(36 * sizeof(double));
	  tr->partitionData[i].EI   = (double*)malloc(30 * sizeof(double));	  
	  tr->partitionData[i].tipVector   = (double *)malloc(384 * sizeof(double)); 	  
	  tr->partitionData[i].substRates = (double *)malloc(15 * sizeof(double));	  
	  tr->partitionData[i].frequencies = (double*)malloc(6 * sizeof(double));
	  tr->partitionData[i].symmetryVector = (int *)malloc(15  * sizeof(int));
	  tr->partitionData[i].frequencyGrouping = (int *)malloc(6  * sizeof(int));
	  tr->partitionData[i].nonGTR = FALSE;
	  break;
	case SECONDARY_DATA_7:
	  tr->partitionData[i].EIGN = (double*)malloc(6 * sizeof(double));
	  tr->partitionData[i].EV   = (double*)malloc(49 * sizeof(double));
	  tr->partitionData[i].EI   = (double*)malloc(42 * sizeof(double));	  
	  tr->partitionData[i].tipVector   = (double *)malloc(896 * sizeof(double)); 	  
	  tr->partitionData[i].substRates = (double *)malloc(21 * sizeof(double));	  
	  tr->partitionData[i].frequencies = (double*)malloc(7 * sizeof(double));
	  tr->partitionData[i].symmetryVector = (int *)malloc(21  * sizeof(int));
	  tr->partitionData[i].frequencyGrouping = (int *)malloc(7  * sizeof(int));
	  tr->partitionData[i].nonGTR = FALSE;
	  break; 
	default:
	  assert(0);
	}
      
      tr->partitionData[i].gammaRates = (double*)malloc(sizeof(double) * 4);
      tr->partitionData[i].yVector = (unsigned char **)malloc(sizeof(unsigned char*) * (tr->mxtips + 1));
      tr->partitionData[i].xVector = (double **)malloc(sizeof(double*) * tr->mxtips);
      tr->partitionData[i].pVector = (parsimonyVector **)malloc(sizeof(parsimonyVector*) * tr->mxtips);
      tr->partitionData[i].expVector = (int **)malloc(sizeof(int*) * tr->mxtips);
      tr->partitionData[i].mxtips  = tr->mxtips;	       
    }

 
 
  for(model = 0; model < (size_t)tr->NumberOfModels; model++)
    {     
      size_t width = tr->partitionData[model].width;
      
      myLength += width;
      
      switch(tr->partitionData[model].dataType)
	{
	case AA_DATA:
	  memoryRequirements += ((size_t)tr->aaIncrement * width);
	  break;
	case DNA_DATA:
	  memoryRequirements += ((size_t)tr->dnaIncrement * width);
	  break;
	case BINARY_DATA:
	  memoryRequirements += ((size_t)tr->binaryIncrement * width);
	  break;
	case SECONDARY_DATA:
	  memoryRequirements += ((size_t)tr->secondaryIncrement * width);
	  break;
	case SECONDARY_DATA_6:
	  memoryRequirements += ((size_t)tr->secondaryIncrement6 * width);
	  break; 
	case SECONDARY_DATA_7:
	  memoryRequirements += ((size_t)tr->secondaryIncrement7 * width);
	  break;
	default:
	  assert(0);
	}    
    }
 
  if(tid == 0)
    {     
      tr->perSiteLL       = (double *)malloc((size_t)tr->cdta->endsite * sizeof(double));
      assert(tr->perSiteLL != NULL);
    } 

  tr->likelihoodArray = (double *)malloc((size_t)tr->mxtips * memoryRequirements * sizeof(double));	
  assert(tr->likelihoodArray != NULL);

  tr->expArray = (int *)malloc(myLength * (size_t)tr->mxtips * sizeof(int));
  assert(tr->expArray != NULL);

  tr->sumBuffer  = (double *)malloc(memoryRequirements * sizeof(double));
  assert(tr->sumBuffer != NULL);

  tr->y_ptr = (unsigned char *)malloc(myLength * (size_t)(tr->mxtips) * sizeof(unsigned char));
  assert(tr->y_ptr != NULL);
  
  assert(4 * sizeof(double) > sizeof(parsimonyVector));
 
  tr->perSiteLLPtr     = (double*) malloc(myLength * sizeof(double));

  tr->wrPtr            = (double*) malloc(myLength * sizeof(double));
  assert(tr->wrPtr != NULL);

  tr->wr2Ptr           = (double*) malloc(myLength * sizeof(double));
  assert(tr->wr2Ptr != NULL);
   
  tr->wgtPtr           = (int*)    malloc(myLength * sizeof(int));
  assert(tr->wgtPtr != NULL);
   
  tr->invariantPtr     = (int*)    malloc(myLength * sizeof(int));
  assert(tr->invariantPtr != NULL);
   
  tr->rateCategoryPtr  = (int*)    malloc(myLength * sizeof(int));
  assert(tr->rateCategoryPtr != NULL);

  if(tid == 0)
    {
      for (i = (size_t)tr->mxtips + 1; (i <= 2*((size_t)tr->mxtips) - 2); i++) 
	{    
	  p = tr->nodep[i];                
	  p->x = 1;	   
	  p->next->x       = 0;
	  p->next->next->x = 0;        
	}      
    }
}








inline static void sendTraversalInfo(tree *localTree, tree *tr)
{
  /* the one below is a hack we are re-assigning the local pointer to the global one
     the memcpy version below is just for testing and preparing the 
     fine-grained MPI BlueGene version */

  if(1)
    {
      localTree->td[0] = tr->td[0];
    }
  else
    {     
      localTree->td[0].count = tr->td[0].count;     
      memcpy(localTree->td[0].ti, tr->td[0].ti, localTree->td[0].count * sizeof(traversalInfo));
    }
}


static void collectDouble(double *dst, double *src, tree *tr, int n, int tid)
{  
  int model, i;

  for(model = 0; model < tr->NumberOfModels; model++)
    {
      for(i = tr->partitionData[model].lower; i < tr->partitionData[model].upper; i++)
	{
	  if(i % n == tid)	    
	    dst[i] = src[i];	    
	}
    }     
}

/* Zsolt : important function for Pthreads sync */

static void execFunction(tree *tr, tree *localTree, int tid, int n)
{
  double volatile result;
  int
    i,
    currentJob,
    parsimonyResult,
    model,
    localCounter,
    globalCounter;
  
  currentJob = threadJob >> 16;

  switch(currentJob)      
    {
      /* initialization only */
    case THREAD_INIT_PARTITION:     
      initPartition(tr, localTree, tid);
      break;
    case THREAD_ALLOC_LIKELIHOOD:   
      allocNodex(localTree, tid, n);    
      threadFixModelIndices(tr, localTree, tid, n);
      break;	   
    case THREAD_FIX_MODEL_INDICES:    
      threadFixModelIndices(tr, localTree, tid, n);
      break;
    case THREAD_EVALUATE_PARSIMONY:
      sendTraversalInfo(localTree, tr);
      parsimonyResult = evaluateParsimonyIterative(localTree);
      reductionBufferParsimony[tid] = parsimonyResult;
      break;
    case THREAD_NEWVIEW_PARSIMONY:
      sendTraversalInfo(localTree, tr);
      newviewParsimonyIterative(localTree);  
      break;
    case THREAD_EVALUATE:
      sendTraversalInfo(localTree, tr); 
      result = evaluateIterative(localTree, FALSE);                  

      if(localTree->NumberOfModels > 1)
	{
	  for(model = 0; model < localTree->NumberOfModels; model++)
	    reductionBuffer[tid * localTree->NumberOfModels + model] = localTree->perPartitionLH[model];
	}
      else
	reductionBuffer[tid] = result;

      if(tid > 0)
	{
	  for(model = 0; model < localTree->NumberOfModels; model++)
	    localTree->executeModel[model] = TRUE;
	}
      break;
    case THREAD_NEWVIEW_MASKED:
      sendTraversalInfo(localTree, tr);
      memcpy(localTree->executeModel, tr->executeModel, sizeof(boolean) * localTree->NumberOfModels);
      newviewIterative(localTree);
      if(tid > 0)
	{
	  for(model = 0; model < localTree->NumberOfModels; model++)
	    localTree->executeModel[model] = TRUE;
	}
      break;
    case THREAD_NEWVIEW:
      sendTraversalInfo(localTree, tr);
      newviewIterative(localTree);
      break;
    case THREAD_MAKENEWZ_FIRST:
      {
	volatile double 
	  dlnLdlz[NUM_BRANCHES], 
	  d2lnLdlz2[NUM_BRANCHES];

	sendTraversalInfo(localTree, tr);
	memcpy(localTree->coreLZ,   tr->coreLZ,   sizeof(double) *  localTree->numBranches);
	memcpy(localTree->executeModel, tr->executeModel, sizeof(boolean) * localTree->NumberOfModels);

	makenewzIterative(localTree);     
	execCore(localTree, dlnLdlz, d2lnLdlz2);
	if(!tr->multiBranch)
	  {
	    reductionBuffer[tid]    = dlnLdlz[0];
	    reductionBufferTwo[tid] = d2lnLdlz2[0];
	  }
	else
	  {
	    for(i = 0; i < localTree->NumberOfModels; i++)
	      {
		reductionBuffer[tid * localTree->NumberOfModels + i]    = dlnLdlz[i];
		reductionBufferTwo[tid * localTree->NumberOfModels + i] = d2lnLdlz2[i];
	      }
	  }

	if(tid > 0)
	  {
	    for(model = 0; model < localTree->NumberOfModels; model++)
	      localTree->executeModel[model] = TRUE;
	  }
      }
      break;
    case THREAD_MAKENEWZ:
      {
	volatile double 
	  dlnLdlz[NUM_BRANCHES], 
	  d2lnLdlz2[NUM_BRANCHES];

	memcpy(localTree->coreLZ,   tr->coreLZ,   sizeof(double) *  localTree->numBranches);
	memcpy(localTree->executeModel, tr->executeModel, sizeof(boolean) * localTree->NumberOfModels);
	execCore(localTree, dlnLdlz, d2lnLdlz2); 

	if(!tr->multiBranch)
	  {
	    reductionBuffer[tid]    = dlnLdlz[0];
	    reductionBufferTwo[tid] = d2lnLdlz2[0];
	  }
	else
	  {
	    for(i = 0; i < localTree->NumberOfModels; i++)
	      {
		reductionBuffer[tid * localTree->NumberOfModels + i]    = dlnLdlz[i];
		reductionBufferTwo[tid * localTree->NumberOfModels + i] = d2lnLdlz2[i];	      
	      }
	  }
	if(tid > 0)
	  {
	    for(model = 0; model < localTree->NumberOfModels; model++)
	      localTree->executeModel[model] = TRUE;
	  }
      }	
      break;   
    case THREAD_COPY_RATES:
      if(tid > 0)
	{
	  for(model = 0; model < localTree->NumberOfModels; model++)
	    {
	      switch(localTree->partitionData[model].dataType)
		{ 
		case DNA_DATA:
		  memcpy(localTree->partitionData[model].EIGN,        tr->partitionData[model].EIGN,        3 * sizeof(double));
		  memcpy(localTree->partitionData[model].EV,          tr->partitionData[model].EV,          16 * sizeof(double));
		  memcpy(localTree->partitionData[model].EI,          tr->partitionData[model].EI,          12 * sizeof(double));	      	     	      
		  memcpy(localTree->partitionData[model].tipVector,   tr->partitionData[model].tipVector,   64 * sizeof(double));
		  break;
		case AA_DATA:
		  memcpy(localTree->partitionData[model].EIGN,        tr->partitionData[model].EIGN,       19 * sizeof(double));
		  memcpy(localTree->partitionData[model].EV,          tr->partitionData[model].EV,          400 * sizeof(double));
		  memcpy(localTree->partitionData[model].EI,          tr->partitionData[model].EI,          380 * sizeof(double));	      	     	      
		  memcpy(localTree->partitionData[model].tipVector,   tr->partitionData[model].tipVector,   460 * sizeof(double));		  
		  break;
		case BINARY_DATA:
		  memcpy(localTree->partitionData[model].EIGN,        tr->partitionData[model].EIGN,      1 * sizeof(double));
		  memcpy(localTree->partitionData[model].EV,          tr->partitionData[model].EV,        4 * sizeof(double));
		  memcpy(localTree->partitionData[model].EI,          tr->partitionData[model].EI,        2 * sizeof(double));	      	     	      
		  memcpy(localTree->partitionData[model].tipVector,   tr->partitionData[model].tipVector, 8 * sizeof(double));		  
		  break;
		case SECONDARY_DATA:
		  memcpy(localTree->partitionData[model].EIGN,        tr->partitionData[model].EIGN,      15 * sizeof(double));
		  memcpy(localTree->partitionData[model].EV,          tr->partitionData[model].EV,        256 * sizeof(double));
		  memcpy(localTree->partitionData[model].EI,          tr->partitionData[model].EI,        240 * sizeof(double));	      	     	      
		  memcpy(localTree->partitionData[model].tipVector,   tr->partitionData[model].tipVector, 4096 * sizeof(double));		  
		  break;
		case SECONDARY_DATA_6:
		  memcpy(localTree->partitionData[model].EIGN,        tr->partitionData[model].EIGN,      5 * sizeof(double));
		  memcpy(localTree->partitionData[model].EV,          tr->partitionData[model].EV,        36 * sizeof(double));
		  memcpy(localTree->partitionData[model].EI,          tr->partitionData[model].EI,        30 * sizeof(double));	      	     	      
		  memcpy(localTree->partitionData[model].tipVector,   tr->partitionData[model].tipVector, 384 * sizeof(double));		  
		  break;
		case SECONDARY_DATA_7:
		  memcpy(localTree->partitionData[model].EIGN,        tr->partitionData[model].EIGN,      6 * sizeof(double));
		  memcpy(localTree->partitionData[model].EV,          tr->partitionData[model].EV,        49 * sizeof(double));
		  memcpy(localTree->partitionData[model].EI,          tr->partitionData[model].EI,        42 * sizeof(double));	      	     	      
		  memcpy(localTree->partitionData[model].tipVector,   tr->partitionData[model].tipVector, 896 * sizeof(double));		  
		  break;
		default:
		  assert(0);
		}
	    }
	}
      break;
    case THREAD_OPT_RATE:
      if(tid > 0)
	{       
	  memcpy(localTree->executeModel, tr->executeModel, localTree->NumberOfModels * sizeof(boolean));      

	  for(model = 0; model < localTree->NumberOfModels; model++)
	    {
	      switch(localTree->partitionData[model].dataType)
		{ 
		case DNA_DATA:
		  memcpy(localTree->partitionData[model].EIGN,        tr->partitionData[model].EIGN,        3 * sizeof(double));
		  memcpy(localTree->partitionData[model].EV,          tr->partitionData[model].EV,          16 * sizeof(double));
		  memcpy(localTree->partitionData[model].EI,          tr->partitionData[model].EI,          12 * sizeof(double));	      	     	      
		  memcpy(localTree->partitionData[model].tipVector,   tr->partitionData[model].tipVector,   64 * sizeof(double));
		  break;
		case AA_DATA:
		  memcpy(localTree->partitionData[model].EIGN,        tr->partitionData[model].EIGN,       19 * sizeof(double));
		  memcpy(localTree->partitionData[model].EV,          tr->partitionData[model].EV,          400 * sizeof(double));
		  memcpy(localTree->partitionData[model].EI,          tr->partitionData[model].EI,          380 * sizeof(double));	      	     	      
		  memcpy(localTree->partitionData[model].tipVector,   tr->partitionData[model].tipVector,   460 * sizeof(double));		  
		  break;
		case BINARY_DATA:
		  memcpy(localTree->partitionData[model].EIGN,        tr->partitionData[model].EIGN,        1 * sizeof(double));
		  memcpy(localTree->partitionData[model].EV,          tr->partitionData[model].EV,          4 * sizeof(double));
		  memcpy(localTree->partitionData[model].EI,          tr->partitionData[model].EI,          2 * sizeof(double));	      	     	      
		  memcpy(localTree->partitionData[model].tipVector,   tr->partitionData[model].tipVector,   8 * sizeof(double));		  
		  break;
		case SECONDARY_DATA:
		  memcpy(localTree->partitionData[model].EIGN,        tr->partitionData[model].EIGN,        15 * sizeof(double));
		  memcpy(localTree->partitionData[model].EV,          tr->partitionData[model].EV,          256 * sizeof(double));
		  memcpy(localTree->partitionData[model].EI,          tr->partitionData[model].EI,          240 * sizeof(double));	      	     	      
		  memcpy(localTree->partitionData[model].tipVector,   tr->partitionData[model].tipVector,   4096 * sizeof(double));		  
		  break;
		case SECONDARY_DATA_6:
		  memcpy(localTree->partitionData[model].EIGN,        tr->partitionData[model].EIGN,        5 * sizeof(double));
		  memcpy(localTree->partitionData[model].EV,          tr->partitionData[model].EV,          36 * sizeof(double));
		  memcpy(localTree->partitionData[model].EI,          tr->partitionData[model].EI,          30 * sizeof(double));	      	     	      
		  memcpy(localTree->partitionData[model].tipVector,   tr->partitionData[model].tipVector,   384 * sizeof(double));		  
		  break;
		case SECONDARY_DATA_7:
		  memcpy(localTree->partitionData[model].EIGN,        tr->partitionData[model].EIGN,        6 * sizeof(double));
		  memcpy(localTree->partitionData[model].EV,          tr->partitionData[model].EV,          49 * sizeof(double));
		  memcpy(localTree->partitionData[model].EI,          tr->partitionData[model].EI,          42 * sizeof(double));	      	     	      
		  memcpy(localTree->partitionData[model].tipVector,   tr->partitionData[model].tipVector,   896 * sizeof(double));		  
		  break;
		default:
		  assert(0);
		}
	    }
	}

      result = evaluateIterative(localTree, FALSE);      
     

      if(localTree->NumberOfModels > 1)	
	{
	  for(model = 0; model < localTree->NumberOfModels; model++)	   
	    reductionBuffer[tid * localTree->NumberOfModels + model] = localTree->perPartitionLH[model];  
	}
      else
	reductionBuffer[tid] = result;


      if(tid > 0)
	{
	  for(model = 0; model < localTree->NumberOfModels; model++)
	    localTree->executeModel[model] = TRUE;
	}
      break;
    case THREAD_COPY_INVAR:
      if(tid > 0)
	{
	  for(model = 0; model < localTree->NumberOfModels; model++)
	    localTree->partitionData[model].propInvariant = tr->partitionData[model].propInvariant;
	}
      break;
    case THREAD_OPT_INVAR:
      if(tid > 0)
	{	 
	  memcpy(localTree->executeModel, tr->executeModel, localTree->NumberOfModels * sizeof(boolean));
	  for(model = 0; model < localTree->NumberOfModels; model++)
	    localTree->partitionData[model].propInvariant = tr->partitionData[model].propInvariant;
	}
      
      result = evaluateIterative(localTree, FALSE);
     
      if(localTree->NumberOfModels > 1)
	{
	  for(model = 0; model < localTree->NumberOfModels; model++)
	    reductionBuffer[tid * localTree->NumberOfModels + model] = localTree->perPartitionLH[model];
	}
      else
	reductionBuffer[tid] = result;
            
      if(tid > 0)
	{
	  for(model = 0; model < localTree->NumberOfModels; model++)
	    localTree->executeModel[model] = TRUE;
	}
      break;
    case THREAD_COPY_ALPHA:
      if(tid > 0)
	{
	  for(model = 0; model < localTree->NumberOfModels; model++)
	    {
	      memcpy(localTree->partitionData[model].gammaRates, tr->partitionData[model].gammaRates, sizeof(double) * 4);
	      localTree->partitionData[model].alpha = tr->partitionData[model].alpha;
	    }
	}
      break;
    case THREAD_OPT_ALPHA: 
      if(tid > 0)
	{	 
	  memcpy(localTree->executeModel, tr->executeModel, localTree->NumberOfModels * sizeof(boolean));
	  for(model = 0; model < localTree->NumberOfModels; model++)	
	    memcpy(localTree->partitionData[model].gammaRates, tr->partitionData[model].gammaRates, sizeof(double) * 4);
	}
       
      result = evaluateIterative(localTree, FALSE);
     

      if(localTree->NumberOfModels > 1)
	{
	  for(model = 0; model < localTree->NumberOfModels; model++)
	    reductionBuffer[tid *  localTree->NumberOfModels + model] = localTree->perPartitionLH[model];
	}
      else
	reductionBuffer[tid] = result;

      if(tid > 0)
	{
	  for(model = 0; model < localTree->NumberOfModels; model++)
	    localTree->executeModel[model] = TRUE;
	}
      break;  
    case THREAD_RESET_MODEL:     
      if(tid > 0)
	{		  	 	  
	  for(model = 0; model < localTree->NumberOfModels; model++)
	    {
	       switch(localTree->partitionData[model].dataType)
		 { 
		 case DNA_DATA:
		   memcpy(localTree->partitionData[model].EIGN,        tr->partitionData[model].EIGN,        3 * sizeof(double));
		   memcpy(localTree->partitionData[model].EV,          tr->partitionData[model].EV,          16 * sizeof(double));
		   memcpy(localTree->partitionData[model].EI,          tr->partitionData[model].EI,          12 * sizeof(double));	      
		   memcpy(localTree->partitionData[model].substRates,  tr->partitionData[model].substRates,  5 * sizeof(double));	      
		   memcpy(localTree->partitionData[model].frequencies, tr->partitionData[model].frequencies, 4 * sizeof(double));	      
		   memcpy(localTree->partitionData[model].tipVector,   tr->partitionData[model].tipVector,   64 * sizeof(double));
		   break;
		 case AA_DATA:
		   memcpy(localTree->partitionData[model].EIGN,        tr->partitionData[model].EIGN,       19 * sizeof(double));
		   memcpy(localTree->partitionData[model].EV,          tr->partitionData[model].EV,          400 * sizeof(double));
		   memcpy(localTree->partitionData[model].EI,          tr->partitionData[model].EI,          380 * sizeof(double));	      
		   memcpy(localTree->partitionData[model].substRates,  tr->partitionData[model].substRates,  190 * sizeof(double));	      
		   memcpy(localTree->partitionData[model].frequencies, tr->partitionData[model].frequencies, 20 * sizeof(double));	      
		   memcpy(localTree->partitionData[model].tipVector,   tr->partitionData[model].tipVector,   460 * sizeof(double));		  
		   break;
		 case BINARY_DATA:
		   memcpy(localTree->partitionData[model].EIGN,        tr->partitionData[model].EIGN,        1 * sizeof(double));
		   memcpy(localTree->partitionData[model].EV,          tr->partitionData[model].EV,          4 * sizeof(double));
		   memcpy(localTree->partitionData[model].EI,          tr->partitionData[model].EI,          2 * sizeof(double));	      
		   memcpy(localTree->partitionData[model].substRates,  tr->partitionData[model].substRates,  1 * sizeof(double));	      
		   memcpy(localTree->partitionData[model].frequencies, tr->partitionData[model].frequencies, 2 * sizeof(double));	      
		   memcpy(localTree->partitionData[model].tipVector,   tr->partitionData[model].tipVector,   8 * sizeof(double));		  
		   break;
		 case SECONDARY_DATA:
		   memcpy(localTree->partitionData[model].EIGN,        tr->partitionData[model].EIGN,        15 * sizeof(double));
		   memcpy(localTree->partitionData[model].EV,          tr->partitionData[model].EV,          256 * sizeof(double));
		   memcpy(localTree->partitionData[model].EI,          tr->partitionData[model].EI,          240 * sizeof(double));	      
		   memcpy(localTree->partitionData[model].substRates,  tr->partitionData[model].substRates,  120 * sizeof(double));	      
		   memcpy(localTree->partitionData[model].frequencies, tr->partitionData[model].frequencies, 16 * sizeof(double));	      
		   memcpy(localTree->partitionData[model].tipVector,   tr->partitionData[model].tipVector,   4096 * sizeof(double));		  
		   break;
		 case SECONDARY_DATA_6:
		   memcpy(localTree->partitionData[model].EIGN,        tr->partitionData[model].EIGN,        5 * sizeof(double));
		   memcpy(localTree->partitionData[model].EV,          tr->partitionData[model].EV,          36 * sizeof(double));
		   memcpy(localTree->partitionData[model].EI,          tr->partitionData[model].EI,          30 * sizeof(double));	      
		   memcpy(localTree->partitionData[model].substRates,  tr->partitionData[model].substRates,  15 * sizeof(double));	      
		   memcpy(localTree->partitionData[model].frequencies, tr->partitionData[model].frequencies, 6 * sizeof(double));	      
		   memcpy(localTree->partitionData[model].tipVector,   tr->partitionData[model].tipVector,   384 * sizeof(double));		  
		   break; 
		 case SECONDARY_DATA_7:
		   memcpy(localTree->partitionData[model].EIGN,        tr->partitionData[model].EIGN,        6 * sizeof(double));
		   memcpy(localTree->partitionData[model].EV,          tr->partitionData[model].EV,          49 * sizeof(double));
		   memcpy(localTree->partitionData[model].EI,          tr->partitionData[model].EI,          42 * sizeof(double));	      
		   memcpy(localTree->partitionData[model].substRates,  tr->partitionData[model].substRates,  21 * sizeof(double));	      
		   memcpy(localTree->partitionData[model].frequencies, tr->partitionData[model].frequencies, 7 * sizeof(double));	      
		   memcpy(localTree->partitionData[model].tipVector,   tr->partitionData[model].tipVector,   896 * sizeof(double));		  
		   break;
		 default:
		   assert(0);
		 }
	    

	       memcpy(localTree->partitionData[model].gammaRates, tr->partitionData[model].gammaRates, sizeof(double) * 4);
	       localTree->partitionData[model].alpha = tr->partitionData[model].alpha;
	       localTree->partitionData[model].propInvariant = tr->partitionData[model].propInvariant;	     
	    }
	}
      break;
      /* Zsolt add a case here */
      /* case THREAD_COPY_BRANCHINFO_POINTER */
 
    case THREAD_COPY_INIT_MODEL:      
      if(tid > 0)
	{		  
	  localTree->NumberOfCategories = tr->NumberOfCategories;
	  localTree->rateHetModel       = tr->rateHetModel;
	  
	  for(model = 0; model < localTree->NumberOfModels; model++)
	    {
	       switch(localTree->partitionData[model].dataType)
		 { 
		 case DNA_DATA:
		   memcpy(localTree->partitionData[model].EIGN,        tr->partitionData[model].EIGN,        3 * sizeof(double));
		   memcpy(localTree->partitionData[model].EV,          tr->partitionData[model].EV,          16 * sizeof(double));
		   memcpy(localTree->partitionData[model].EI,          tr->partitionData[model].EI,          12 * sizeof(double));	      
		   memcpy(localTree->partitionData[model].substRates,  tr->partitionData[model].substRates,  5 * sizeof(double));	      
		   memcpy(localTree->partitionData[model].frequencies, tr->partitionData[model].frequencies, 4 * sizeof(double));	      
		   memcpy(localTree->partitionData[model].tipVector,   tr->partitionData[model].tipVector,   64 * sizeof(double));
		   break;
		 case AA_DATA:
		   memcpy(localTree->partitionData[model].EIGN,        tr->partitionData[model].EIGN,       19 * sizeof(double));
		   memcpy(localTree->partitionData[model].EV,          tr->partitionData[model].EV,          400 * sizeof(double));
		   memcpy(localTree->partitionData[model].EI,          tr->partitionData[model].EI,          380 * sizeof(double));	      
		   memcpy(localTree->partitionData[model].substRates,  tr->partitionData[model].substRates,  190 * sizeof(double));	      
		   memcpy(localTree->partitionData[model].frequencies, tr->partitionData[model].frequencies, 20 * sizeof(double));	      
		   memcpy(localTree->partitionData[model].tipVector,   tr->partitionData[model].tipVector,   460 * sizeof(double));		  
		   break;
		 case BINARY_DATA:
		   memcpy(localTree->partitionData[model].EIGN,        tr->partitionData[model].EIGN,        1 * sizeof(double));
		   memcpy(localTree->partitionData[model].EV,          tr->partitionData[model].EV,          4 * sizeof(double));
		   memcpy(localTree->partitionData[model].EI,          tr->partitionData[model].EI,          2 * sizeof(double));	      
		   memcpy(localTree->partitionData[model].substRates,  tr->partitionData[model].substRates,  1 * sizeof(double));	      
		   memcpy(localTree->partitionData[model].frequencies, tr->partitionData[model].frequencies, 2 * sizeof(double));	      
		   memcpy(localTree->partitionData[model].tipVector,   tr->partitionData[model].tipVector,   8 * sizeof(double));		  
		   break;
		 case SECONDARY_DATA:
		   memcpy(localTree->partitionData[model].EIGN,        tr->partitionData[model].EIGN,        15 * sizeof(double));
		   memcpy(localTree->partitionData[model].EV,          tr->partitionData[model].EV,          256 * sizeof(double));
		   memcpy(localTree->partitionData[model].EI,          tr->partitionData[model].EI,          240 * sizeof(double));	      
		   memcpy(localTree->partitionData[model].substRates,  tr->partitionData[model].substRates,  120 * sizeof(double));	      
		   memcpy(localTree->partitionData[model].frequencies, tr->partitionData[model].frequencies, 16 * sizeof(double));	      
		   memcpy(localTree->partitionData[model].tipVector,   tr->partitionData[model].tipVector,   4096 * sizeof(double));		  
		   break;
		 case SECONDARY_DATA_6:
		   memcpy(localTree->partitionData[model].EIGN,        tr->partitionData[model].EIGN,        5 * sizeof(double));
		   memcpy(localTree->partitionData[model].EV,          tr->partitionData[model].EV,          36 * sizeof(double));
		   memcpy(localTree->partitionData[model].EI,          tr->partitionData[model].EI,          30 * sizeof(double));	      
		   memcpy(localTree->partitionData[model].substRates,  tr->partitionData[model].substRates,  15 * sizeof(double));	      
		   memcpy(localTree->partitionData[model].frequencies, tr->partitionData[model].frequencies, 6 * sizeof(double));	      
		   memcpy(localTree->partitionData[model].tipVector,   tr->partitionData[model].tipVector,   384 * sizeof(double));		  
		   break; 
		 case SECONDARY_DATA_7:
		   memcpy(localTree->partitionData[model].EIGN,        tr->partitionData[model].EIGN,        6 * sizeof(double));
		   memcpy(localTree->partitionData[model].EV,          tr->partitionData[model].EV,          49 * sizeof(double));
		   memcpy(localTree->partitionData[model].EI,          tr->partitionData[model].EI,          42 * sizeof(double));	      
		   memcpy(localTree->partitionData[model].substRates,  tr->partitionData[model].substRates,  21 * sizeof(double));	      
		   memcpy(localTree->partitionData[model].frequencies, tr->partitionData[model].frequencies, 7 * sizeof(double));	      
		   memcpy(localTree->partitionData[model].tipVector,   tr->partitionData[model].tipVector,   896 * sizeof(double));		  
		   break;
		 default:
		   assert(0);
		 }
	    

	       memcpy(localTree->partitionData[model].gammaRates, tr->partitionData[model].gammaRates, sizeof(double) * 4);
	       localTree->partitionData[model].alpha = tr->partitionData[model].alpha;
	       localTree->partitionData[model].propInvariant = tr->partitionData[model].propInvariant;
	       localTree->partitionData[model].lower      = tr->partitionData[model].lower;
	       localTree->partitionData[model].upper      = tr->partitionData[model].upper;
	    }
	  memcpy(localTree->cdta->patrat,        tr->cdta->patrat,      localTree->originalCrunchedLength * sizeof(double));    
	  memcpy(localTree->cdta->patratStored, tr->cdta->patratStored, localTree->originalCrunchedLength * sizeof(double)); 
	}

      /* TODO: really need this here ? */
      
      /*computeFraction(localTree, tid, n);*/
       
       for(model = 0; model < localTree->NumberOfModels; model++)
	 {
	   int localIndex;
	   for(i = localTree->partitionData[model].lower, localIndex = 0; i <  localTree->partitionData[model].upper; i++)
	     if(i % n == tid)
	       {
		 localTree->partitionData[model].wgt[localIndex]          = tr->cdta->aliaswgt[i];
		 localTree->partitionData[model].wr[localIndex]           = tr->cdta->wr[i];
		 localTree->partitionData[model].wr2[localIndex]          = tr->cdta->wr2[i];
		 localTree->partitionData[model].invariant[localIndex]    = tr->invariant[i];
		 localTree->partitionData[model].rateCategory[localIndex] = tr->cdta->rateCategory[i];
		 localIndex++;
	       }
	   /* fix other things as well? */
	   
	   /*for(j = 1; j <= localTree->mxtips; j++)		 
	     memcpy(localTree->partitionData[model].yVector[j], &(tr->yVector[j][currentLower]), sizeof(char) * width);	*/
	 }                          
      
      break; 

    case THREAD_PARSIMONY_RATCHET:
      for(model = 0; model < localTree->NumberOfModels; model++)
	{
	  int localIndex;
	  for(i = localTree->partitionData[model].lower, localIndex = 0; i <  localTree->partitionData[model].upper; i++)
	    if(i % n == tid)
	      {
		localTree->partitionData[model].wgt[localIndex]          = tr->cdta->aliaswgt[i];	
		localIndex++;
	      }
	}
      break;
    case THREAD_RATE_CATS:
      sendTraversalInfo(localTree, tr);
      if(tid > 0)
	{
	  localTree->lower_spacing = tr->lower_spacing;
	  localTree->upper_spacing = tr->upper_spacing;
	}

      optRateCatPthreads(localTree, localTree->lower_spacing, localTree->upper_spacing, localTree->lhs, n, tid);      
      
      if(tid > 0)
	{
	  collectDouble(tr->cdta->patrat,       localTree->cdta->patrat,         localTree, n, tid);    
	  collectDouble(tr->cdta->patratStored, localTree->cdta->patratStored,   localTree, n, tid);      
	  collectDouble(tr->lhs,                localTree->lhs,                  localTree, n, tid); 	 
	}      
      break;
    case THREAD_COPY_RATE_CATS:
      if(tid > 0)
	{	 	  
	  localTree->NumberOfCategories = tr->NumberOfCategories;
	  memcpy(localTree->cdta->patrat,       tr->cdta->patrat,         localTree->originalCrunchedLength * sizeof(double));    
	  memcpy(localTree->cdta->patratStored, tr->cdta->patratStored,   localTree->originalCrunchedLength * sizeof(double)); 
	}

     
      for(model = 0; model < localTree->NumberOfModels; model++)
	{
	  for(localCounter = 0, i = localTree->partitionData[model].lower;  i < localTree->partitionData[model].upper; i++)
	    {
	      if(i % n == tid)
		{			   		 
		  localTree->partitionData[model].wr[localCounter]           = tr->cdta->wr[i];
		  localTree->partitionData[model].wr2[localCounter]          = tr->cdta->wr2[i];	   
		  localTree->partitionData[model].rateCategory[localCounter] = tr->cdta->rateCategory[i];
		  localCounter++;
		}	     
	    }
	}     
      break;
    case THREAD_CAT_TO_GAMMA:
      if(tid > 0)	
	localTree->rateHetModel = tr->rateHetModel;   
      break;
    case THREAD_GAMMA_TO_CAT:
      if(tid > 0)	
	localTree->rateHetModel = tr->rateHetModel;   
      break;
    case THREAD_EVALUATE_VECTOR:
      evaluateIterative(localTree, TRUE); 
      for(model = 0, globalCounter = 0; model < localTree->NumberOfModels; model++)
	{
	  for(localCounter = 0, i = localTree->partitionData[model].lower;  i < localTree->partitionData[model].upper; i++)
	    {
	      if(i % n == tid)
		{	
		  tr->perSiteLL[globalCounter] =  localTree->partitionData[model].perSiteLL[localCounter];		 
		  localCounter++;
		}	     
	      globalCounter++;
	    }
	}     
      break;
      /*case THREAD_COMPUTE_AVERAGE:
      if(tid > 0)
	{
	  localTree->numberOfBipartitions = tr->numberOfBipartitions;
	  localTree->v1 = tr->v1;
	  localTree->v2 = tr->v2;
	}
      threadComputeAverage(localTree, tid);
      break;
    case THREAD_COMPUTE_PEARSON:
       if(tid > 0)
	{	 
	  localTree->avg1 = tr->avg1;
	  localTree->avg2 = tr->avg2;
	}
      threadComputePearson(localTree, tid);   
      break;
    case THREAD_MAKE_VECTORS:
      if(tid > 0)
	{	 
	  localTree->v1 = tr->v1;
	  localTree->v2 = tr->v2;
	  localTree->b  = tr->b;
	  localTree->perm = tr->perm;
	  localTree->numberOfBootstopTrees = tr->numberOfBootstopTrees;	  
	}
      threadMakeVector(localTree, tid);   
      break;*/
    case THREAD_COPY_PARAMS:
      if(tid > 0)
	{
	  for(model = 0; model < localTree->NumberOfModels; model++)
	    {
	      switch(localTree->partitionData[model].dataType)
		{ 
		case DNA_DATA:
		  memcpy(localTree->partitionData[model].EIGN,        tr->partitionData[model].EIGN,        3 * sizeof(double));
		  memcpy(localTree->partitionData[model].EV,          tr->partitionData[model].EV,          16 * sizeof(double));
		  memcpy(localTree->partitionData[model].EI,          tr->partitionData[model].EI,          12 * sizeof(double));	      
		  memcpy(localTree->partitionData[model].substRates,  tr->partitionData[model].substRates,  5 * sizeof(double));	      
		  memcpy(localTree->partitionData[model].frequencies, tr->partitionData[model].frequencies, 4 * sizeof(double));	      
		  memcpy(localTree->partitionData[model].tipVector,   tr->partitionData[model].tipVector,   64 * sizeof(double));
		  break;
		case AA_DATA:
		  memcpy(localTree->partitionData[model].EIGN,        tr->partitionData[model].EIGN,       19 * sizeof(double));
		  memcpy(localTree->partitionData[model].EV,          tr->partitionData[model].EV,          400 * sizeof(double));
		  memcpy(localTree->partitionData[model].EI,          tr->partitionData[model].EI,          380 * sizeof(double));	      
		  memcpy(localTree->partitionData[model].substRates,  tr->partitionData[model].substRates,  190 * sizeof(double));	      
		  memcpy(localTree->partitionData[model].frequencies, tr->partitionData[model].frequencies, 20 * sizeof(double));	      
		  memcpy(localTree->partitionData[model].tipVector,   tr->partitionData[model].tipVector,   460 * sizeof(double));		  
		  break;
		case BINARY_DATA:
		  memcpy(localTree->partitionData[model].EIGN,        tr->partitionData[model].EIGN,        1 * sizeof(double));
		  memcpy(localTree->partitionData[model].EV,          tr->partitionData[model].EV,          4 * sizeof(double));
		  memcpy(localTree->partitionData[model].EI,          tr->partitionData[model].EI,          2 * sizeof(double));	      
		  memcpy(localTree->partitionData[model].substRates,  tr->partitionData[model].substRates,  1 * sizeof(double));	      
		  memcpy(localTree->partitionData[model].frequencies, tr->partitionData[model].frequencies, 2 * sizeof(double));	      
		  memcpy(localTree->partitionData[model].tipVector,   tr->partitionData[model].tipVector,   8 * sizeof(double));		  
		  break;
		case SECONDARY_DATA:
		  memcpy(localTree->partitionData[model].EIGN,        tr->partitionData[model].EIGN,        15 * sizeof(double));
		  memcpy(localTree->partitionData[model].EV,          tr->partitionData[model].EV,          256 * sizeof(double));
		  memcpy(localTree->partitionData[model].EI,          tr->partitionData[model].EI,          240 * sizeof(double));	      
		  memcpy(localTree->partitionData[model].substRates,  tr->partitionData[model].substRates,  120 * sizeof(double));	      
		  memcpy(localTree->partitionData[model].frequencies, tr->partitionData[model].frequencies, 16 * sizeof(double));	      
		  memcpy(localTree->partitionData[model].tipVector,   tr->partitionData[model].tipVector,   4096 * sizeof(double));		  
		  break;
		case SECONDARY_DATA_6:
		  memcpy(localTree->partitionData[model].EIGN,        tr->partitionData[model].EIGN,        5 * sizeof(double));
		  memcpy(localTree->partitionData[model].EV,          tr->partitionData[model].EV,          36 * sizeof(double));
		  memcpy(localTree->partitionData[model].EI,          tr->partitionData[model].EI,          30 * sizeof(double));	      
		  memcpy(localTree->partitionData[model].substRates,  tr->partitionData[model].substRates,  15 * sizeof(double));	      
		  memcpy(localTree->partitionData[model].frequencies, tr->partitionData[model].frequencies, 6 * sizeof(double));	      
		  memcpy(localTree->partitionData[model].tipVector,   tr->partitionData[model].tipVector,   384 * sizeof(double));		  
		  break;
		case SECONDARY_DATA_7:
		  memcpy(localTree->partitionData[model].EIGN,        tr->partitionData[model].EIGN,        6 * sizeof(double));
		  memcpy(localTree->partitionData[model].EV,          tr->partitionData[model].EV,          49 * sizeof(double));
		  memcpy(localTree->partitionData[model].EI,          tr->partitionData[model].EI,          42 * sizeof(double));	      
		  memcpy(localTree->partitionData[model].substRates,  tr->partitionData[model].substRates,  21 * sizeof(double));	      
		  memcpy(localTree->partitionData[model].frequencies, tr->partitionData[model].frequencies, 7 * sizeof(double));	      
		  memcpy(localTree->partitionData[model].tipVector,   tr->partitionData[model].tipVector,   896 * sizeof(double));		  
		  break;
		default:
		  assert(0);
		}
	    }
	}
      break;
#ifdef _USE_PTHREADS_MULTIGRAIN
    case THREAD_COPY_BRANCHINFO_POINTER:
      /* 
	 so here we just copy the adresses of tr->bInfo
	 and the integer tr->numberOfBranches to the thread-local 
	 variables of the tree data structure. only need to do this 
	 for the worker threads, i.e., threads with thread id (tid) > 0 .....
      */

      if(tid > 0)
	{
	  localTree->bInfo = tr->bInfo;
	  localTree->numberOfBranches = tr->numberOfBranches;
	}

      /* 
	 now we also set some pointers to the contiguous data structures we will need in addition
	 to the scaling and likelihood vectors in order to compute the likelihood.
	 Here we also need to copy those adresses to the respective fields of the master thread tree structure, since 
	 they are also uninitialzed at the master up to this point here ......
      */

      localTree->contiguousRateCategory = tr->cdta->rateCategory;
      localTree->contiguousWgt = tr->cdta->aliaswgt;
      localTree->contiguousInvariant = tr->invariant;
      localTree->contiguousTips = tr->yVector;
      break;
    case THREAD_GATHER_LIKELIHOOD:
      {	
	/* pointers to contiguous likelihood vectors 
	   in the branchInfo data structure, i.e., the destination of 
	   our gather operation */

	double 
	  *leftContigousVector = (double*)NULL, 
	  *rightContigousVector = (double*)NULL; 

	/* pointers to contiguous scaling vectors */ 
	int 
	  *leftContigousScalingVector = (int*)NULL,
	  *rightContigousScalingVector = (int*)NULL,
	  /* and a couple of integers variables and counters we will need */
	  globalColumnCount = 0,
	  globalCount       = 0,
	  rightNumber, 
	  leftNumber;

	/* first get the current branch Number */

	if(tid > 0)	  
	  localTree->branchNumber = tr->branchNumber;
	   	
	/* now get the node numbers that correspond to the tree nodes
	   that define the branch */
	
	leftNumber =  localTree->bInfo[localTree->branchNumber].leftNodeNumber;
	rightNumber =  localTree->bInfo[localTree->branchNumber].rightNodeNumber;

	/* set the adresses of the contiguous vectors associated to the current branch */

	leftContigousVector = localTree->bInfo[localTree->branchNumber].left;
	leftContigousScalingVector = localTree->bInfo[localTree->branchNumber].leftScaling;
	
	rightContigousVector = localTree->bInfo[localTree->branchNumber].right;
	rightContigousScalingVector = localTree->bInfo[localTree->branchNumber].rightScaling;


	/* and now loop over all partitions */

	for(model = 0; model < localTree->NumberOfModels; model++)
	  {     
	    size_t 
	      blockRequirements;
	      	
	    double
	      *leftStridedVector  =  (double *)NULL,
	      *rightStridedVector =  (double *)NULL;  

	    int 
	      *leftStridedScalingVector  =  (int *)NULL,
	      *rightStridedScalingVector =  (int *)NULL,
	      /* with the two localCounters we keep track of the respective offsets within the partition */
	      localColumnCount = 0,
	      localCount = 0;

	    /* 
	       if the left node of the branch is not a tip there must be a local strided likelihood vector 
	       available for this partition which we access via 
	       localTree->partitionData[model].xVector[]
	       and 
	       localTree->partitionData[model].expVector[]
	       and a bit of fiddlinga round with the node numbers 
	    */

	    if(!isTip(leftNumber, localTree->mxtips))
	      {
		leftStridedVector        = localTree->partitionData[model].xVector[leftNumber - localTree->mxtips - 1];
		leftStridedScalingVector = localTree->partitionData[model].expVector[leftNumber - localTree->mxtips - 1];
	      }

	    /* same thing for the right node of the branch */

	    if(!isTip(rightNumber, localTree->mxtips))
	      {
		rightStridedVector        = localTree->partitionData[model].xVector[rightNumber - localTree->mxtips - 1];
		rightStridedScalingVector = localTree->partitionData[model].expVector[rightNumber - localTree->mxtips - 1];
	      }

	    /* just make sure that not both end sof the branch are tips, if they are, something must have went badly wrong */

	    assert(!(isTip(leftNumber, localTree->mxtips) && isTip(rightNumber, localTree->mxtips)));

	    /* get the likelihood vector increment */

	    switch(localTree->partitionData[model].dataType) 
	      	{
		case AA_DATA:
		  blockRequirements = (size_t)localTree->aaIncrement;
		 break;
		case DNA_DATA:
		  blockRequirements = (size_t)localTree->dnaIncrement ;
		  break;
		case BINARY_DATA:
		  blockRequirements = (size_t)localTree->binaryIncrement ;
		  break;
		case SECONDARY_DATA:
		  blockRequirements = (size_t)localTree->secondaryIncrement ;
		  break;
		case SECONDARY_DATA_6:
		  blockRequirements = (size_t)localTree->secondaryIncrement6;
		  break;
		case SECONDARY_DATA_7:
		  blockRequirements = (size_t)localTree->secondaryIncrement7;
		  break;		 
		default:
		  assert(0);
		} 		    	  

	    /* 
	       And now do a nasty loop over the columns of the present partition 
	       As I said  localTree->partitionData[model].lower and localTree->partitionData[model].upper
	       always denote the full length and real starting and ending positions of a partition ....
	       That's why we can use this to initialize the globalColumnCount here ....
	     */

	    for(globalColumnCount = localTree->partitionData[model].lower; globalColumnCount < localTree->partitionData[model].upper; globalColumnCount++)
	      {
		/* if the current thread actually owns the present column we need to copy some stuff */

		if(globalColumnCount % n == tid)
		  {
		    /* if the left node is not a tip, i.e., if leftStridedVector != NULL do the strided copy into the contiguous vector */
		    if(leftStridedVector)
		      {
			memcpy(&leftContigousVector[globalCount], &leftStridedVector[localCount], sizeof(double) * blockRequirements);
			leftContigousScalingVector[globalColumnCount] = leftStridedScalingVector[localColumnCount];
		      }

		    /* same for right node */
		    if(rightStridedVector)
		      {
			memcpy(&rightContigousVector[globalCount], &rightStridedVector[localCount], sizeof(double) * blockRequirements);
			rightContigousScalingVector[globalColumnCount] = rightStridedScalingVector[localColumnCount];
		      }
		    /* increment the local offsets */
		    localColumnCount++;
		    localCount += blockRequirements;
		  }		
		
		/* 
		   always increment the global offset such that we know to which part of the contiguous array 
		   we need to copy next 
		*/
		   
		globalCount += blockRequirements;
	      }

	    /* do some nasty assertions to make sure that we did not screw up the indexing stuff */

	    assert(localColumnCount == localTree->partitionData[model].width);
	    assert(localCount == (localTree->partitionData[model].width * (int)blockRequirements));

	  }
      }
      break;
    case THREAD_TEST_CLASSIFY:
      {
	int i;
	/* 
	   and now this is just for testing, we just loop over the branches 
	   in bInfo and use the information gathered above to compute the likelihood 
	   over contiguous array concurrently. Here we just use a cyclic distribution again 
	   so thread 0 computes the likelihood on branches 0, 2, 4 etc and thread one 
	   on 1, 3, 5 etc. This is actually very similar to what we will then need to 
	   do for computing the actual insertions ... 
	   evaluateClassify is a function that will just use the info in the branchInfo structure 
	   to compute the likelihood. We will need a slightly modified version of this later on for 
	   computing the insertions, and unfortunately analogous functions for branch length optimization 
	   and likelihood vector re-computation .....
	*/

	for(i = 0; i < localTree->numberOfBranches; i++)
	  {
	    if(i % n == tid)	    
	      {
		branchInfo *b = &(localTree->bInfo[i]);      
		double r = evaluateClassify(localTree,  b);
		printf("%d %d %f\n", tid, i, r);
	      }
	  }
      }
      break;
#endif
    default:
      printf("Job %d\n", currentJob);
      assert(0);
    }
}



void masterBarrier(int jobType, tree *tr) 
{
  const int n = NumberOfThreads;
  int i, sum;
  
  jobCycle = !jobCycle;   
  threadJob = (jobType << 16) + jobCycle; 
  
  execFunction(tr, tr, 0, n);      

  do
    {     
      for(i = 1, sum = 1; i < n; i++)
	sum += barrierBuffer[i];
    }
  while(sum < n);

  for(i = 1; i < n; i++)
    barrierBuffer[i] = 0;
}




static void *likelihoodThread(void *tData)
{ 
  threadData *td = (threadData*)tData; 
  tree 
    *tr = td->tr, 
    *localTree = (tree *)malloc(sizeof(tree));
  int    
    myCycle = 0;

  const int n = NumberOfThreads;
  const int tid             = td->threadNumber; 

  printf("\nThis is RAxML Worker Pthread Number: %d\n", tid);
   
  while(1)
    {                
      while (myCycle == threadJob);
      myCycle = threadJob;

      execFunction(tr, localTree, tid, n);
     
      barrierBuffer[tid] = 1;
    }

  return (void*)NULL;
}

static void startPthreads(tree *tr)
{  
  pthread_t *threads;
  pthread_attr_t attr;
  int rc, t;
  threadData *tData; 
  
  jobCycle        = 0; 
  threadJob       = 0;
  
  printf("\nThis is the RAxML Master Pthread\n");    

  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);
  
  threads    = (pthread_t *)malloc(NumberOfThreads * sizeof(pthread_t));
  tData      = (threadData *)malloc(NumberOfThreads * sizeof(threadData));
  reductionBuffer          = (double *)malloc(sizeof(double) *  NumberOfThreads * tr->NumberOfModels);
  reductionBufferTwo       = (double *)malloc(sizeof(double) *  NumberOfThreads * tr->NumberOfModels);
  reductionBufferThree     = (double *)malloc(sizeof(double) *  NumberOfThreads * tr->NumberOfModels);
  reductionBufferParsimony = (int *)malloc(sizeof(int) *  NumberOfThreads);
  barrierBuffer            = (int *)malloc(sizeof(int) *  NumberOfThreads);

 


  for(t = 0; t < NumberOfThreads; t++)
    barrierBuffer[t] = 0;      

  for(t = 1; t < NumberOfThreads; t++)
    {
      tData[t].tr  = tr;
      tData[t].threadNumber = t;
      rc = pthread_create(&threads[t], &attr, likelihoodThread, (void *)(&tData[t]));
      if(rc)
	{
	  printf("ERROR; return code from pthread_create() is %d\n", rc);
	  exit(-1);
	}
    }
}



#endif


/*************************************************************************************************************************************************************/





typedef struct {
  double lh;
  int tree;
  double weight;
} elw;

static int elwCompare(const void *p1, const void *p2)
{
  elw *rc1 = (elw *)p1;
  elw *rc2 = (elw *)p2;
  
  double i = rc1->weight;
  double j = rc2->weight;
  
  if (i > j)
    return (-1);
  if (i < j)
    return (1);
  return (0);
}


static void computeAllLHs(tree *tr, analdef *adef, char *bootStrapFileName)
{
  int 
    numberOfTrees = 0,   
    i; 
  double 
    bestLH = unlikely;    
  bestlist *bestT;
  FILE *infoFile, *result;
  
  infoFile = myfopen(infoFileName, "a");
  result   = myfopen(resultFileName, "w");

  bestT = (bestlist *) malloc(sizeof(bestlist));
  bestT->ninit = 0;
  initBestTree(bestT, 1, tr->mxtips);

  INFILE = myfopen(bootStrapFileName, "r"); 
  numberOfTrees = countTrees(INFILE);
 
 
  printf("\n\nFound %d trees in File %s\n\n", numberOfTrees, bootStrapFileName);
  fprintf(infoFile, "\n\nFound %d trees in File %s\n\n", numberOfTrees, bootStrapFileName);
 
  for(i = 0; i < numberOfTrees; i++)
    {              
      treeReadLen(INFILE, tr, adef);      
      
      if(i == 0)
	{
	  modOpt(tr, adef, TRUE, adef->likelihoodEpsilon);
	  printf("Model optimization, first Tree: %f\n", tr->likelihood);
	  fprintf(infoFile, "Model optimization, first Tree: %f\n", tr->likelihood);
	  bestLH = tr->likelihood;
	  resetBranches(tr);
	}
      
      treeEvaluate(tr, 2);
      Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE, 
		  TRUE, adef, SUMMARIZE_LH);
                 
      fprintf(result, "%s", tr->tree_string);
      
      saveBestTree(bestT, tr);

      if(tr->likelihood > bestLH)		
	bestLH   = tr->likelihood;	
      printf("Tree %d Likelihood %f\n", i, tr->likelihood);
      fprintf(infoFile, "Tree %d Likelihood %f\n", i, tr->likelihood);
    }        
    
  recallBestTree(bestT, 1, tr);  
  evaluateGeneric(tr, tr->start);
  printf("Model optimization, %f <-> %f\n", bestLH, tr->likelihood); 
  fprintf(infoFile, "Model optimization, %f <-> %f\n", bestLH, tr->likelihood); 
  modOpt(tr, adef, TRUE, adef->likelihoodEpsilon);
  treeEvaluate(tr, 2);
  printf("Model optimization, %f <-> %f\n", bestLH, tr->likelihood);
  fprintf(infoFile, "Model optimization, %f <-> %f\n", bestLH, tr->likelihood); 

  printf("\nAll evaluated trees with branch lengths written to File: %s\n", resultFileName);
  fprintf(infoFile, "\nAll evaluated trees with branch lengths written to File: %s\n", resultFileName);

  fclose(INFILE); 
  fclose(infoFile);
  fclose(result);
  exit(0);
}




static void computeELW(tree *tr, analdef *adef, char *bootStrapFileName)
{
  int 
    numberOfTrees = 0,   
    i, k; 

  /* 
     double 
     bestLH = unlikely,
     elwSum = 0.0;     
  */

  FILE *infoFile;
  int *originalRateCategories = (int*)malloc(tr->cdta->endsite * sizeof(int));      
  int *originalInvariant      = (int*)malloc(tr->cdta->endsite * sizeof(int));
  long startSeed;           
  double **lhs;
  double **lhweights;
  elw *bootweights;

  infoFile = myfopen(infoFileName, "a");   

  initModel(tr, tr->rdta, tr->cdta, adef); 

  INFILE = myfopen(bootStrapFileName, "r");
  numberOfTrees = countTrees(INFILE); 

  if(numberOfTrees < 2)
    {
      printf("Error, there is only one tree in file %s which you want to use to conduct an ELW test\n", bootStrapFileName);

      exit(-1);
    }
  
  printf("\n\nFound %d trees in File %s\n\n", numberOfTrees, bootStrapFileName);
  fprintf(infoFile, "\n\nFound %d trees in File %s\n\n", numberOfTrees, bootStrapFileName);

  bootweights = (elw *)malloc(sizeof(elw) * numberOfTrees);

  lhs = (double **)malloc(sizeof(double *) * numberOfTrees);

  for(k = 0; k < numberOfTrees; k++)
    lhs[k] = (double *)malloc(sizeof(double) * adef->multipleRuns);

  lhweights = (double **)malloc(sizeof(double *) * numberOfTrees);

  for(k = 0; k < numberOfTrees; k++)
    lhweights[k] = (double *)malloc(sizeof(double) * adef->multipleRuns);
 

  treeReadLen(INFILE, tr, adef);      
  modOpt(tr, adef, TRUE, adef->likelihoodEpsilon);
  rewind(INFILE);

  /*
    This is for testing only !
    for(i = 0; i < numberOfTrees; i++)
    {
      treeReadLen(INFILE, tr, adef);
      treeEvaluate(tr, 2.0);
      bootweights[i].lh = tr->likelihood;
    }
    rewind(INFILE);
  */

  printf("Model optimization, first Tree: %f\n", tr->likelihood);
  fprintf(infoFile, "Model optimization, first Tree: %f\n", tr->likelihood);
  

  memcpy(originalRateCategories, tr->cdta->rateCategory, sizeof(int) * tr->cdta->endsite);
  memcpy(originalInvariant,      tr->invariant,          sizeof(int) * tr->cdta->endsite);

  assert(adef->boot > 0);
  /* TODO this is ugly, should be passed as param to computenextreplicate() */
  startSeed = adef->boot;
  

  for(i = 0; i < numberOfTrees; i++)
    {              
      treeReadLen(INFILE, tr, adef);      
      resetBranches(tr);
      adef->rapidBoot = startSeed;     
     
      for(k = 0; k < adef->multipleRuns; k++)
	{
	  computeNextReplicate(tr, &adef->rapidBoot, originalRateCategories, originalInvariant, TRUE);
	  /*computeNextReplicate(tr, adef, originalRateCategories, originalInvariant);*/

	  if(k == 0)
	    treeEvaluate(tr, 2.0);
	  else
	    treeEvaluate(tr, 0.5);
	  /*printf("%d %d %f\n", i, k, tr->likelihood);*/
	  lhs[i][k] = tr->likelihood;	 
	}          

      reductionCleanup(tr, originalRateCategories, originalInvariant);
    }        

  

  for(k = 0; k < adef->multipleRuns; k++)
    {
      double best = unlikely;
      double sum = 0.0;

      for(i = 0; i < numberOfTrees; i++)
	if(lhs[i][k] > best)
	  best = lhs[i][k];

      for(i = 0; i < numberOfTrees; i++)
	lhweights[i][k] = exp(lhs[i][k] - best);

      for(i = 0; i < numberOfTrees; i++)
	sum += lhweights[i][k];

      for(i = 0; i < numberOfTrees; i++)
	lhweights[i][k] = lhweights[i][k] / sum;

    }
  
  
  
  for(i = 0; i < numberOfTrees; i++)
    {
      double sum = 0.0;
      
      for(k = 0; k < adef->multipleRuns; k++)
	sum += lhweights[i][k];

      bootweights[i].weight = sum / ((double)adef->multipleRuns);     
      bootweights[i].tree   = i;
    }

  qsort(bootweights, numberOfTrees, sizeof(elw), elwCompare);


  {
    double sum = 0.0;
    
    /*printf("Tree\t Posterior Probability \t Cumulative posterior probability \t Original Likelihood\n");*/
    printf("Tree\t Posterior Probability \t Cumulative posterior probability\n");
    fprintf(infoFile, "Tree\t Posterior Probability \t Cumulative posterior probability\n");
    for(i = 0; i < numberOfTrees; i++)
      {
	 sum += bootweights[i].weight;
	 /*printf("%d\t\t %f \t\t %f \t\t\t %f\n", bootweights[i].tree, bootweights[i].weight, sum,  bootweights[i].lh);*/
	 printf("%d\t\t %f \t\t %f\n", bootweights[i].tree, bootweights[i].weight, sum); 
	 fprintf(infoFile, "%d\t\t %f \t\t %f\n", bootweights[i].tree, bootweights[i].weight, sum); 
      }
  }

  free(originalRateCategories);
  free(originalInvariant);

  fclose(INFILE); 
  fclose(infoFile); 
  exit(0);
}



static void computeDistances(tree *tr, analdef *adef)
{
  int i, j, modelCounter;
  double z0[NUM_BRANCHES];
  double result[NUM_BRANCHES];
  double t;  
  char distanceFileName[1024];

  FILE 
    *out;

  strcpy(distanceFileName,         workdir);   
  strcat(distanceFileName,         "RAxML_distances.");
  strcat(distanceFileName,         run_id);

  out = myfopen(distanceFileName, "w");

  modOpt(tr, adef, TRUE, adef->likelihoodEpsilon);
    
  printBothOpen("\nLog Likelihood Score after parameter optimization: %f\n\n", tr->likelihood); 
  printBothOpen("\nComputing pairwise ML-distances ...\n");  
  
  for(modelCounter = 0; modelCounter < tr->NumberOfModels; modelCounter++)
    z0[modelCounter] = defaultz;

  t = gettime();

  for(i = 1; i <= tr->mxtips; i++)
    for(j = i + 1; j <= tr->mxtips; j++)
      {        
	double z, x;
	
	makenewzGenericDistance(tr, 10, z0, result, i, j);

	if(tr->multiBranch) 
	  {
	    int k;
	   	    
	    for(k = 0, x = 0.0; k < tr->numBranches; k++)
	      {
		assert(tr->partitionContributions[k] != -1.0);
		assert(tr->fracchanges[k] != -1.0);
		z = result[k];
		if (z < zmin) 
		  z = zmin;      	 		
		x += (-log(z) * tr->fracchanges[k]) * tr->partitionContributions[k];
	      }	
	  }
	else
	  {
	    z = result[0];
	    if (z < zmin) 
	      z = zmin;      	 
	    x = -log(z) * tr->fracchange;
	  }
	      
	/*printf("%s-%s \t %f\n", tr->nameList[i], tr->nameList[j], x);*/
	fprintf(out, "%s %s \t %f\n", tr->nameList[i], tr->nameList[j], x);
      }

  fclose(out);

  t = gettime() - t;

  printBothOpen("\nTime for pair-wise ML distance computation of %d distances: %f seconds\n", 
		 (tr->mxtips * tr->mxtips - tr->mxtips) / 2, t);
  printBothOpen("\nDistances written to file: %s\n", distanceFileName);

 

  exit(0);
}


static void olafIsANastyGuy(tree *tr, analdef *adef)
{
  assert(tr->rateHetModel == GAMMA);

  modOpt(tr, adef, TRUE, adef->likelihoodEpsilon);
  
  printf("Likelihood %f\n", tr->likelihood);
  
  gammaToCat(tr);

  


  exit(0);
}


int main (int argc, char *argv[]) 
{   
  rawdata      *rdta;
  cruncheddata *cdta;
  tree         *tr;         
  analdef      *adef; 
  int 
    i,
    countGTR = 0,
    countOtherModel = 0;
   
#ifdef PARALLEL
  MPI_Init(&argc, &argv); 
  MPI_Comm_rank(MPI_COMM_WORLD, &processID);  
  MPI_Comm_size(MPI_COMM_WORLD, &numOfWorkers);       
  if(processID == 0)
    printf("\nThis is the RAxML MPI Master process\n");
  else
    printf("\nThis is the RAxML MPI Worker Process Number: %d\n", processID);
#else
  processID = 0;
#endif
 
  masterTime = gettime();            

  adef = (analdef *)malloc(sizeof(analdef));
  rdta = (rawdata *)malloc(sizeof(rawdata));
  cdta = (cruncheddata *)malloc(sizeof(cruncheddata));     
  tr   = (tree *)malloc(sizeof(tree));

  initAdef(adef); 
  get_args(argc,argv, adef, tr);      
  setRateHetAndDataIncrement(tr, adef);
  
  
  /* 
     This is a very ugly numerical bug fix, that intends to avoid the unaesthetic phenomena
     that can occur during model param optimization due to the dependency between parameters 
     alpha and invar which are NOT independent from each other. 
     When using P-Invar set likelihood epsilon to a lower value!  

     TODO-MIX this is very ugly !

  */
           
  if(adef->mrpEncoder)
    encodeMRP(tr, rdta);
 
  
  readData(adef, rdta, cdta, tr);   
  
  

  checkOutgroups(tr, adef);
  makeFileNames();    

  if(adef->useInvariant && adef->likelihoodEpsilon > 0.001)
    {
      FILE *info = myfopen(infoFileName, "a");

      printf("\nYou are using a proportion of Invariable sites estimate, although I don't\n");
      printf("like it. The likelihood epsilon \"-f e\" will be automatically lowered to 0.001\n");
      printf("to avoid unfavorable effects caused by simultaneous optimization of alpha and P-Invar\n");

      
      fprintf(info, "\nYou are using a proportion of Invariable sites estimate, although I don't\n");
      fprintf(info, "like it. The likelihood epsilon \"-f e\" will be automatically lowered to 0.001\n");
      fprintf(info, "to avoid unfavorable effects caused by simultaneous optimization of alpha and P-Invar\n");

      adef->likelihoodEpsilon = 0.001;      


      fclose(info);
    }

  if(adef->useExcludeFile)
    {
      handleExcludeFile(tr, adef, rdta);
      exit(0);
    }
  
  if(adef->mode != SEQUENCE_SIMILARITY_FILTER)      
    checkSequences(tr, rdta, adef);    
  else
    {     
      reduceBySequenceSimilarity(tr, rdta, adef);
      exit(0);
    }

  if(adef->mode == SPLIT_MULTI_GENE)
    {     
      splitMultiGene(tr, rdta);
      exit(0);
    }
  
  if(adef->mode == CHECK_ALIGNMENT)
    {
      printf("Alignment format can be read by RAxML \n");
      exit(0);
    }  
   			       
  makeweights(adef, rdta, cdta, tr);  
  makevalues(rdta, cdta, tr, adef);     

  for(i = 0; i < tr->NumberOfModels; i++)
    {
      if(tr->partitionData[i].dataType == AA_DATA)
	{
	  if(tr->partitionData[i].protModels == GTR)
	    countGTR++;
	  else
	    countOtherModel++;
	}
    }

  if(countGTR > 0 && countOtherModel > 0)
    {
      printf("Error, it is only allowed to conduct partitioned AA analyses\n");
      printf("with a GTR model of AA substiution, if all AA partitions are assigned\n");
      printf("the GTR model.\n\n");

      printf("The following partitions do not use GTR:\n");
      
      for(i = 0; i < tr->NumberOfModels; i++)
	{
	  if(tr->partitionData[i].dataType == AA_DATA && tr->partitionData[i].protModels != GTR)
	    printf("Partition %s\n", tr->partitionData[i].partitionName);
	}
      printf("exiting ...\n");
      errorExit(-1);
    }

  if(countGTR > 0 && tr->NumberOfModels > 1)
    {
      FILE *info = fopen(infoFileName, "a");     
	
      printBoth(info, "You are using the GTR model of AA substitution!\n");
      printBoth(info, "GTR parameters for AA substiution will automatically be estimated\n");
      printBoth(info, "jointly (GTR params will be linked) across all partitions to avoid over-parametrization!\n\n\n");
	
      fclose(info);
    }

#ifdef _USE_PTHREADS
  startPthreads(tr);
  masterBarrier(THREAD_INIT_PARTITION, tr);
  masterBarrier(THREAD_ALLOC_LIKELIHOOD, tr);
#else    
  allocNodex(tr);
#endif


 
  printModelAndProgramInfo(tr, adef, argc, argv);    

  switch(adef->mode)
    {  
    case THOROUGH_PARSIMONY:
      makeParsimonyTreeThorough(tr, adef);
      break;
    case SUPER_FAST:
      superFast(tr, adef);
      break;
    case CLASSIFY_ML:     
      initModel(tr, rdta, cdta, adef); 
      getStartingTree(tr, adef);
      exit(0);
      break;
    case GENERATE_BS:
      generateBS(tr, adef);
      exit(0);
      break;
    case COMPUTE_ELW:
      computeELW(tr, adef, bootStrapFile); 
      exit(0);
      break;
    case COMPUTE_LHS: 
      initModel(tr, rdta, cdta, adef); 
      computeAllLHs(tr, adef, bootStrapFile);
      exit(0);
      break;
    case COMPUTE_BIPARTITION_CORRELATION:
      compareBips(tr, bootStrapFile);
      exit(0);
      break;
    case COMPUTE_RF_DISTANCE:
      computeRF(tr, bootStrapFile);
      exit(0);
      break;
    case BOOTSTOP_ONLY:
      computeBootStopOnly(tr, bootStrapFile);
      exit(0);
      break;
    case DISTANCE_MODE:
      initModel(tr, rdta, cdta, adef);       
      getStartingTree(tr, adef);
      computeDistances(tr, adef);
      break;      
    case  PARSIMONY_ADDITION:
      initModel(tr, rdta, cdta, adef);  
      getStartingTree(tr, adef);
      printStartingTree(tr, adef, TRUE); 
      break;
    case PER_SITE_LL:    
      initModel(tr, rdta, cdta, adef);  
      computePerSiteLLs(tr, adef, bootStrapFile);           
      break;    
    case TREE_EVALUATION: 
      initModel(tr, rdta, cdta, adef);  
      getStartingTree(tr, adef);      
      if(adef->likelihoodTest)	
	computeLHTest(tr, adef, bootStrapFile);
      else
	{	
	  modOpt(tr, adef, TRUE, adef->likelihoodEpsilon);	       	  
	  printLog(tr, adef, TRUE);         
	  printResult(tr, adef, TRUE);     
	}
      break;
    case CALC_BIPARTITIONS: 
      initModel(tr, rdta, cdta, adef);  
      calcBipartitions(tr, adef, tree_file, bootStrapFile);   
      break;
    case BIG_RAPID_MODE:
      if(adef->boot)	
	doBootstrap(tr, adef, rdta, cdta);
      else
	{   
	  
	  if(adef->rapidBoot)
	    {
	      initModel(tr, rdta, cdta, adef); 

	      doAllInOne(tr, adef);
	    }
	  else     	    	     
	    doInference(tr, adef, rdta, cdta);          			     	   
	}
      break;
    case OLAF_OPTION:
      initModel(tr, rdta, cdta, adef);
      getStartingTree(tr, adef);
      olafIsANastyGuy(tr, adef);
      break;
    default:
      assert(0);
    }

  finalizeInfoFile(tr, adef);

#ifdef PARALLEL
  MPI_Finalize();
#endif

  return 0;
}
