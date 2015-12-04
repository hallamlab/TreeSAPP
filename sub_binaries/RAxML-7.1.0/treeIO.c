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

extern FILE *INFILE;
extern char infoFileName[1024];
extern char tree_file[1024];
extern char *likelihood_key;
extern char *ntaxa_key;
extern char *smoothed_key;
extern int partCount;
extern double masterTime;

static int countTips(nodeptr p, int numsp)
{
  if(isTip(p->number, numsp))  
    return 1;    
  {
    nodeptr q;
    int tips = 0;

    q = p->next;
    while(q != p)
      { 
	tips += countTips(q->back, numsp);
	q = q->next;
      } 
    
    return tips;
  }
}



static char *Tree2StringREC(char *treestr, tree *tr, nodeptr p, boolean printBranchLengths, boolean printNames, 
			    boolean printLikelihood, boolean rellTree, boolean finalPrint, int perGene)
{
  double  x, z;
  char  *nameptr;            
      
  if(isTip(p->number, tr->rdta->numsp)) 
    {	       	  
      if(printNames)
	{
	  nameptr = tr->nameList[p->number];     
	  sprintf(treestr, "%s", nameptr);
	}
      else
	sprintf(treestr, "%d", p->number);    
	
      while (*treestr) treestr++;
    }
  else 
    {                 	 
      *treestr++ = '(';
      treestr = Tree2StringREC(treestr, tr, p->next->back, printBranchLengths, printNames, printLikelihood, rellTree, 
			       finalPrint, perGene);
      *treestr++ = ',';
      treestr = Tree2StringREC(treestr, tr, p->next->next->back, printBranchLengths, printNames, printLikelihood, rellTree, 
			       finalPrint, perGene);
      if(p == tr->start->back) 
	{
	  *treestr++ = ',';
	  treestr = Tree2StringREC(treestr, tr, p->back, printBranchLengths, printNames, printLikelihood, rellTree, 
				   finalPrint, perGene);
	}
      *treestr++ = ')';                    
    }

  if(p == tr->start->back) 
    {	      	 
      if(printBranchLengths && !rellTree)
	sprintf(treestr, ":0.0;\n");
      else
	sprintf(treestr, ";\n");	 	  	
    }
  else 
    {                   
      if(rellTree)
	{
	  if(( !isTip(p->number, tr->rdta->numsp)) && 
	     ( !isTip(p->back->number, tr->rdta->numsp)))
	    {			      
	      assert(p->bInf != (branchInfo *)NULL);
	      sprintf(treestr, "%d:%8.20f", p->bInf->support, p->z[0]);
	    }
	  else
	    {	
	      sprintf(treestr, ":%8.20f", p->z[0]);
	    }
	}
      else
	{
	  if(printBranchLengths)
	    {
	      assert(perGene != NO_BRANCHES);
	      
	      if(!tr->multiBranch)
		{
		  assert(tr->fracchange != -1.0);
		  z = p->z[0];
		  if (z < zmin) z = zmin;      	 
		  x = -log(z) * tr->fracchange;
		  sprintf(treestr, ":%8.20f", x);	
		}
	      else
		{
		  if(perGene == SUMMARIZE_LH)
		    {
		      int i;
		      double avgX = 0;
		      
		      for(i = 0; i < tr->numBranches; i++)
			{
			  assert(tr->partitionContributions[i] != -1.0);
			  assert(tr->fracchanges[i] != -1.0);
			  z = p->z[i];
			  if (z < zmin) z = zmin;      	 
			  x = -log(z) * tr->fracchanges[i];
			  avgX += x * tr->partitionContributions[i];
			}
		      sprintf(treestr, ":%8.20f", avgX);
		    }
		  else
		    {	
		      assert(tr->fracchanges[perGene] != -1.0);
		      assert(perGene >= 0 && perGene < tr->numBranches);
		      z = p->z[perGene];
		      if (z < zmin) z = zmin;      	 
		      x = -log(z) * tr->fracchanges[perGene];
		      sprintf(treestr, ":%8.20f", x);
		    }
		}
	    }
	  else
	    sprintf(treestr, "%s", "\0");	    
	}      
    }
  
  while (*treestr) treestr++;
  return  treestr;
}


static void collectSubtrees(tree *tr, nodeptr *subtrees, int *count, int ogn)
{
  int i;
  for(i = tr->mxtips + 1; i <= tr->mxtips + tr->mxtips - 2; i++)
    {
      nodeptr p, q;
      p = tr->nodep[i];
      if(countTips(p, tr->rdta->numsp) == ogn)
	{
	  subtrees[*count] = p;
	  *count = *count + 1;
	}
      q = p->next;
      while(q != p)
	{
	  if(countTips(q, tr->rdta->numsp) == ogn)
	    {
	      subtrees[*count] = q;
	      *count = *count + 1;
	    }
	  q = q->next;
	}
    }
}

static void checkOM(nodeptr p, int *n, int *c, tree *tr)
{
  if(isTip(p->number, tr->rdta->numsp))
    {
      n[*c] = p->number;
      *c = *c + 1;     
    }
  else
    {
      nodeptr q = p->next;

      while(q != p)
	{
	  checkOM(q->back, n, c, tr);
	  q = q->next;
	}
    }
}
    
static char *rootedTreeREC(char *treestr, tree *tr, nodeptr p, boolean printBranchLengths, boolean printNames, 
			   boolean printLikelihood, boolean rellTree, 
			   boolean finalPrint, analdef *adef, int perGene)
{
  double  x, z;
  char  *nameptr;

  if(isTip(p->number, tr->rdta->numsp)) 
    {	     
      if(printNames)
	{
	  nameptr = tr->nameList[p->number];     
	  sprintf(treestr, "%s", nameptr);
	}
      else
	sprintf(treestr, "%d", p->number);
      
      while (*treestr) treestr++;
    }
  else 
    {
      *treestr++ = '(';
      treestr = rootedTreeREC(treestr, tr, p->next->back, printBranchLengths, printNames, printLikelihood, 
			      rellTree, finalPrint, adef, perGene);
      *treestr++ = ',';
      treestr = rootedTreeREC(treestr, tr, p->next->next->back, printBranchLengths, printNames, printLikelihood, 
			      rellTree, finalPrint, adef, perGene);      
      *treestr++ = ')';            
      }

  if(rellTree)
    {     
      if(!isTip(p->number, tr->rdta->numsp) && !isTip(p->back->number, tr->rdta->numsp))
	{	  
	  assert(p->bInf != (branchInfo *)NULL);
	  sprintf(treestr, "%d:%8.20f", p->bInf->support, p->z[0]);
	}
      else
	{	  
	  sprintf(treestr, ":%8.20f", p->z[0]);
	}
    }
  else
    {
      if(printBranchLengths)
	{	  
	  assert(perGene != NO_BRANCHES);

	  if(!tr->multiBranch)
	    {
	      assert(tr->fracchange != -1.0);
	      z = p->z[0];
	      if (z < zmin) z = zmin;      	 
	      x = -log(z) * tr->fracchange;
	      sprintf(treestr, ":%8.20f", x);	
	    }
	  else
	    {
	      if(perGene == SUMMARIZE_LH)
		{
		  int i;
		  double avgX = 0;
		  
		  for(i = 0; i < tr->numBranches; i++)
		    {
		      assert(tr->fracchanges[i] != -1.0);
		      assert(tr->partitionContributions[i] != -1.0);
		      z = p->z[i];
		      if (z < zmin) z = zmin;      	 
		      x = -log(z) * tr->fracchanges[i];
		      avgX += x * tr->partitionContributions[i];
		    }
		  sprintf(treestr, ":%8.20f", avgX);
		}
	      else
		{
		  assert(tr->fracchanges[perGene] != -1.0);
		  assert(perGene >= 0 && perGene < tr->numBranches);
		  z = p->z[perGene];
		  if (z < zmin) z = zmin;      	 
		  x = -log(z) * tr->fracchanges[perGene];
		  sprintf(treestr, ":%8.20f", x);
		}
	    }	  
	}
      else
	sprintf(treestr, "%s", "\0");
    }

    while (*treestr) treestr++;
    return  treestr;

}

static char *rootedTree(char *treestr, tree *tr, nodeptr p, boolean printBranchLengths, boolean printNames, 
			boolean printLikelihood, boolean rellTree, 
			boolean finalPrint, analdef *adef, int perGene)
{
  double oldz[NUM_BRANCHES];
  int i;
  
  for(i = 0; i < tr->numBranches; i++)
    oldz[i] = p->z[i];

  if(rellTree)    
    {    
      p->z[0] = p->back->z[0] = oldz[0] * 0.5;
      /*printf("%f\n",  p->z[0]);*/
    }
  else
    {
      if(printBranchLengths)
	{
	  double rz, z;
	  assert(perGene != NO_BRANCHES);

	  if(!tr->multiBranch)
	    {
	      assert(tr->fracchange != -1.0);
	      z = -log(p->z[0]) * tr->fracchange;
	      rz = exp(-(z * 0.5)/ tr->fracchange);
	      p->z[0] = p->back->z[0] = rz;
	    }
	  else
	    {
	      if(perGene == SUMMARIZE_LH)
		{				  		
		  int i;	      
		  
		  for(i = 0; i < tr->numBranches; i++)
		    {	
		      assert(tr->fracchanges[i] != -1.0);
		      z    = -log(p->z[i]) * tr->fracchanges[i];	    	      
		      rz   = exp(-(z * 0.5)/ tr->fracchanges[i]);
		      p->z[i] = p->back->z[i] = rz;		    
		    }		 
		}	     	     
	      else
		{		
		  assert(tr->fracchanges[perGene] != -1.0);
		  assert(perGene >= 0 && perGene < tr->numBranches);
		  z = -log(p->z[perGene]) * tr->fracchanges[perGene];
		  rz = exp(-(z * 0.5)/ tr->fracchanges[perGene]);
		  p->z[perGene] = p->back->z[perGene] = rz;	       	      	      
		}
	    }
	}
    }

  *treestr = '(';
  treestr++;
  treestr = rootedTreeREC(treestr, tr, p,  printBranchLengths, printNames, printLikelihood, rellTree, finalPrint, adef, perGene);
  *treestr = ',';
  treestr++;
  treestr = rootedTreeREC(treestr, tr, p->back,  printBranchLengths, printNames, printLikelihood, rellTree, finalPrint, adef, perGene);
  sprintf(treestr, ");\n");
  while (*treestr) treestr++;


  for(i = 0; i < tr->numBranches; i++)
    p->z[i] = p->back->z[i] = oldz[i];  
    
  return  treestr;
}



char *Tree2String(char *treestr, tree *tr, nodeptr p, boolean printBranchLengths, boolean printNames, boolean printLikelihood, 
		  boolean rellTree, 
		  boolean finalPrint, analdef *adef, int perGene)
{ 
  if(finalPrint && adef->outgroup)
    {
      nodeptr startNode = tr->start;

      if(tr->numberOfOutgroups > 1)
	{
	  nodeptr root;
	  nodeptr *subtrees = (nodeptr *)malloc(sizeof(nodeptr) * tr->mxtips);
	  int i, k, count = 0;
	  int *nodeNumbers = (int*)malloc(sizeof(int) * tr->numberOfOutgroups);
	  int *foundVector = (int*)malloc(sizeof(int) * tr->numberOfOutgroups);
	  boolean monophyletic = FALSE;

	  collectSubtrees(tr, subtrees, &count, tr->numberOfOutgroups);

	  /*printf("Found %d subtrees of size  %d\n", count, tr->numberOfOutgroups);*/
	  
	  for(i = 0; (i < count) && (!monophyletic); i++)
	    {
	      int l, sum, nc = 0;
	      for(k = 0; k <  tr->numberOfOutgroups; k++)
		{
		  nodeNumbers[k] = -1;
		  foundVector[k] = 0;
		}

	      checkOM(subtrees[i], nodeNumbers, &nc, tr);	      
	      
	      for(l = 0; l < tr->numberOfOutgroups; l++)
		for(k = 0; k < tr->numberOfOutgroups; k++)
		  {
		    if(nodeNumbers[l] == tr->outgroupNums[k])
		      foundVector[l] = 1;
		  }
	      
	      sum = 0;
	      for(l = 0; l < tr->numberOfOutgroups; l++)
		sum += foundVector[l];
	      
	      if(sum == tr->numberOfOutgroups)
		{	       		  
		  root = subtrees[i];
		  tr->start = root;		
		  /*printf("outgroups are monphyletic!\n");*/
		  monophyletic = TRUE;		  
		}
	      else
		{
		  if(sum > 0)
		    {
		      /*printf("outgroups are NOT monophyletic!\n");*/
		      monophyletic = FALSE;
		    }	     
		}	
	    }
	  
	  if(!monophyletic)
	    {
	      printf("WARNING, outgroups are not monophyletic, using first outgroup \"%s\"\n", tr->nameList[tr->outgroupNums[0]]);
	      printf("from the list to root the tree!\n");
	     
#ifndef PARALLEL
	      {
		FILE *infoFile = myfopen(infoFileName, "a");

		fprintf(infoFile, "\nWARNING, outgroups are not monophyletic, using first outgroup \"%s\"\n", tr->nameList[tr->outgroupNums[0]]);
		fprintf(infoFile, "from the list to root the tree!\n");
		
		fclose(infoFile);
	      }
#endif 

	      tr->start = tr->nodep[tr->outgroupNums[0]];
	      /*Tree2StringREC(treestr, tr, tr->start->back, printBranchLengths, printNames, printLikelihood, rellTree, finalPrint);*/
	      rootedTree(treestr, tr, tr->start->back, printBranchLengths, printNames, printLikelihood, rellTree, finalPrint, adef, perGene);
	    }
	  else
	    {	     
	      if(/*tr->start->tip*/ isTip(tr->start->number, tr->rdta->numsp))
		{
		  printf("Outgroup-Monophyly ERROR; tr->start is a tip \n");
		  errorExit(-1);
		}
	      if(/*tr->start->back->tip*/ isTip(tr->start->back->number, tr->rdta->numsp))
	      	{
		  printf("Outgroup-Monophyly ERROR; tr->start is a tip \n");
		  errorExit(-1);
		}
	      
	      /*	      Tree2StringREC(treestr, tr, tr->start->back, printBranchLengths, printNames, printLikelihood, rellTree, finalPrint);*/
	      rootedTree(treestr, tr, tr->start->back, printBranchLengths, printNames, printLikelihood, rellTree, finalPrint, adef, perGene);
	    }
	  
	  free(foundVector);
	  free(nodeNumbers);
	  free(subtrees);
	}
      else
	{
	  /*printf("Skipping Monophyly Check, only one outgroup\n");*/
	  tr->start = tr->nodep[tr->outgroupNums[0]];
	  /*printf("%d\n", tr->outgroupNums[0]);*/
	  /*Tree2StringREC(treestr, tr, tr->start->back, printBranchLengths, printNames, printLikelihood, rellTree, finalPrint); */
	  rootedTree(treestr, tr, tr->start->back, printBranchLengths, printNames, printLikelihood, rellTree, finalPrint, adef, perGene);
	}      

      tr->start = startNode;
    }
  else
    {
      Tree2StringREC(treestr, tr, p, printBranchLengths, printNames, printLikelihood, rellTree, finalPrint, perGene);  
    }
  while (*treestr) treestr++;
  return treestr;
}


void printTreePerGene(tree *tr, analdef *adef, char *fileName, char *permission)
{  
  FILE *treeFile;
  char extendedTreeFileName[1024];
  char buf[16];
  int i;

  assert(adef->perGeneBranchLengths);
     
  for(i = 0; i < tr->numBranches; i++)	
    {
      strcpy(extendedTreeFileName, fileName);
      sprintf(buf,"%d", i);
      strcat(extendedTreeFileName, ".PARTITION.");
      strcat(extendedTreeFileName, buf);
      /*printf("Partitiuon %d file %s\n", i, extendedTreeFileName);*/
      Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE, TRUE, adef, i);
      treeFile = myfopen(extendedTreeFileName, permission);
      fprintf(treeFile, "%s", tr->tree_string);
      fclose(treeFile);
    }  
    
}



/*=======================================================================*/
/*                         Read a tree from a file                       */
/*=======================================================================*/


/*  1.0.A  Processing of quotation marks in comment removed
 */

static int treeFinishCom (FILE *fp, char **strp)
{
  int  ch;
  
  while ((ch = getc(fp)) != EOF && ch != ']') {
    if (strp != NULL) *(*strp)++ = ch;    /* save character  */
    if (ch == '[') {                      /* nested comment; find its end */
      if ((ch = treeFinishCom(fp, strp)) == EOF)  break;
      if (strp != NULL) *(*strp)++ = ch;  /* save closing ]  */
    }
  }
  
  if (strp != NULL) **strp = '\0';        /* terminate string  */
  return  ch;
} /* treeFinishCom */


static int treeGetCh (FILE *fp)         /* get next nonblank, noncomment character */
{ /* treeGetCh */
  int  ch;

  while ((ch = getc(fp)) != EOF) {
    if (whitechar(ch)) ;
    else if (ch == '[') {                   /* comment; find its end */
      if ((ch = treeFinishCom(fp, (char **) NULL)) == EOF)  break;
    }
    else  break;
  }
  
  return  ch;
} /* treeGetCh */


static boolean treeLabelEnd (int ch)
{
  switch (ch) 
    {
    case EOF:  
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
      return TRUE;
    default:
      break;
    }
  return FALSE;
} 


static boolean  treeGetLabel (FILE *fp, char *lblPtr, int maxlen)
{
  int      ch;
  boolean  done, quoted, lblfound;

  if (--maxlen < 0) 
    lblPtr = (char *) NULL; 
  else 
    if (lblPtr == NULL) 
      maxlen = 0;

  ch = getc(fp);
  done = treeLabelEnd(ch);

  lblfound = ! done;
  quoted = (ch == '\'');
  if (quoted && ! done) 
    {
      ch = getc(fp); 
      done = (ch == EOF);
    }

  while (! done) 
    {
      if (quoted) 
	{
	  if (ch == '\'') 
	    {
	      ch = getc(fp); 
	      if (ch != '\'') 
		break;
	    }
        }
      else 
	if (treeLabelEnd(ch)) break;     

      if (--maxlen >= 0) *lblPtr++ = ch;
      ch = getc(fp);
      if (ch == EOF) break;
    }

  if (ch != EOF)  (void) ungetc(ch, fp);

  if (lblPtr != NULL) *lblPtr = '\0';

  return lblfound;
}


static boolean  treeFlushLabel (FILE *fp)
{ 
  return  treeGetLabel(fp, (char *) NULL, (int) 0);
} 


static int treeFindTipByLabel (char  *str, tree *tr)                    
{
  nodeptr  q;  
  int      i; 

  for (i = 1; i <= tr->mxtips; i++) 
    {
      q = tr->nodep[i];     

      if(! q->back) 
	{        	 	  	        
	  if(strcmp(str, tr->nameList[q->number]) == 0)
	    return i;
	
	  /*
	    Uncomment this to read NEXUS-style trees with taxon number instead of names
	    i = 0;
	    numptr = num;
	    sprintf(numptr,"%d", q->number);
	    while((found = (str[i++] == (ch = *numptr++))) && ch) ;
	    if (found) return n; 
	  */
	}
      
    }

  printf("ERROR: Cannot find tree species: %s\n", str);

  return  0;
}


static int treeFindTipName (FILE *fp, tree *tr)
{
  char    str[nmlngth+2];
  int      n;

  if(treeGetLabel(fp, str, nmlngth+2))
    n = treeFindTipByLabel(str, tr);
  else
    n = 0;
   

  return  n;
} 


static void  treeEchoContext (FILE *fp1, FILE *fp2, int n)
{ /* treeEchoContext */
  int      ch;
  boolean  waswhite;
  
  waswhite = TRUE;
  
  while (n > 0 && ((ch = getc(fp1)) != EOF)) {
    if (whitechar(ch)) {
      ch = waswhite ? '\0' : ' ';
      waswhite = TRUE;
    }
    else {
      waswhite = FALSE;
    }
    
    if (ch > '\0') {putc(ch, fp2); n--;}
  }
} /* treeEchoContext */


static boolean treeProcessLength (FILE *fp, double *dptr)
{
  int  ch;
  
  if ((ch = treeGetCh(fp)) == EOF)  return FALSE;    /*  Skip comments */
  (void) ungetc(ch, fp);
  
  if (fscanf(fp, "%lf", dptr) != 1) {
    printf("ERROR: treeProcessLength: Problem reading branch length\n");
    treeEchoContext(fp, stdout, 40);
    printf("\n");
    return  FALSE;
  }
  
  return  TRUE;
}


static int treeFlushLen (FILE  *fp)
{
  double  dummy;  
  int     ch;
  
  ch = treeGetCh(fp);
  
  if (ch == ':') 
    {
      ch = treeGetCh(fp);
      
      ungetc(ch, fp);
      if(!treeProcessLength(fp, & dummy)) return 0;
      return 1;	  
    }
  
  
  
  if (ch != EOF) (void) ungetc(ch, fp);
  return 1;
} 





static boolean treeNeedCh (FILE *fp, int c1, char *where)
{
  int  c2;
  
  if ((c2 = treeGetCh(fp)) == c1)  return TRUE;
  
  printf("ERROR: Expecting '%c' %s tree; found:", c1, where);
  if (c2 == EOF) {
    printf("End-of-File");
  }
  else {
    ungetc(c2, fp);
    treeEchoContext(fp, stdout, 40);
  }
  putchar('\n');
  return FALSE;
} 



static boolean addElementLen (FILE *fp, tree *tr, nodeptr p, boolean readBranchLengths, boolean readNodeLabels, int *lcount)
{   
  nodeptr  q;
  int      n, ch, fres;
  
  if ((ch = treeGetCh(fp)) == '(') 
    { 
      n = (tr->nextnode)++;
      if (n > 2*(tr->mxtips) - 2) 
	{
	  if (tr->rooted || n > 2*(tr->mxtips) - 1) 
	    {
	      printf("ERROR: Too many internal nodes.  Is tree rooted?\n");
	      printf("       Deepest splitting should be a trifurcation.\n");
	      return FALSE;
	    }
	  else 
	    {
	      assert(!readNodeLabels);
	      tr->rooted = TRUE;
	    }
	}
      
      q = tr->nodep[n];

      if (! addElementLen(fp, tr, q->next, readBranchLengths, readNodeLabels, lcount))        return FALSE;
      if (! treeNeedCh(fp, ',', "in"))             return FALSE;
      if (! addElementLen(fp, tr, q->next->next, readBranchLengths, readNodeLabels, lcount))  return FALSE;
      if (! treeNeedCh(fp, ')', "in"))             return FALSE;
      
      if(readNodeLabels)
	{
	  char label[64];
	  int support;

	  if(treeGetLabel (fp, label, 10))
	    {	      
	      assert(sscanf(label, "%d\n", &support) == 1);
	      /*printf("LABEL %s Number %d\n", label, support);*/
	      p->support = q->support = support;
	      *lcount = *lcount + 1;
	    }
	}
      else	
	(void) treeFlushLabel(fp);
    }
  else 
    {   
      ungetc(ch, fp);
      if ((n = treeFindTipName(fp, tr)) <= 0)          return FALSE;
      q = tr->nodep[n];
      if (tr->start->number > n)  tr->start = q;
      (tr->ntips)++;
    }
  
  if(readBranchLengths)
    {
      double branch;
      if (! treeNeedCh(fp, ':', "in"))                 return FALSE;
      if (! treeProcessLength(fp, &branch))            return FALSE;
      
      /*printf("Branch %8.20f\n", branch);*/
      hookup(p, q, &branch, tr->numBranches);
    }
  else
    {
      fres = treeFlushLen(fp);
      if(!fres) return FALSE;
      
      hookupDefault(p, q, tr->numBranches);
    }
  return TRUE;          
} 











static nodeptr uprootTree (tree *tr, nodeptr p, boolean readBranchLengths)
{
  nodeptr  q, r, s, start;
  int      n, i;              

  for(i = tr->mxtips + 1; i < 2 * tr->mxtips - 1; i++)
    assert(i == tr->nodep[i]->number);


  if (isTip(p->number, tr->rdta->numsp) || p->back) 
    {
      printf("ERROR: Unable to uproot tree.\n");
      printf("       Inappropriate node marked for removal.\n");
      assert(0);
    }

  assert(p->back == (nodeptr)NULL);

  

  
  tr->nextnode = tr->nextnode - 1;
  assert(tr->nextnode < 2 * tr->mxtips);
  
  n = tr->nextnode;               /* last internal node added */
  assert(tr->nodep[tr->nextnode]);

  if (n != tr->mxtips + tr->ntips - 1) 
    {
      printf("ERROR: Unable to uproot tree.  Inconsistent\n");
      printf("       number of tips and nodes for rooted tree.\n");
      assert(0);
    }

  q = p->next->back;                  /* remove p from tree */
  r = p->next->next->back;
  assert(p->back == (nodeptr)NULL);
    
  if(readBranchLengths)
    {
      double b[NUM_BRANCHES];
      int i;
      for(i = 0; i < tr->numBranches; i++)
	b[i] = (r->z[i] + q->z[i]);
      hookup (q, r, b, tr->numBranches);
    }
  else
    {
      hookupDefault(q, r, tr->numBranches);
    }

  if(tr->grouped)
    {
      /*printf("P-NUMBER %d Grouping %d\n", p->number, tr->constraintVector[p->number]);*/
      if(tr->constraintVector[p->number] != 0)
	{
	  printf("Root node to remove shood have top-level grouping of 0\n");
	  assert(0);
	}
    }

  /*start = (r->tip || (! q->tip)) ? r : r->next->next->back;*/

  assert(tr->rdta->numsp == tr->mxtips);
  assert(!(isTip(r->number, tr->rdta->numsp) && isTip(q->number, tr->rdta->numsp)));

  /*if(isTip(r->number, tr->rdta->numsp) || (! isTip(q->number, tr->rdta->numsp)))
    start = r;
  else
  start = r->next->next->back;*/

  /* assert(isTip(start->number, tr->rdta->numsp));*/

  assert(p->number > tr->mxtips);

  if(tr->ntips > 2 && p->number != n) 
    {    	
      q = tr->nodep[n];            /* transfer last node's conections to p */
      r = q->next;
      s = q->next->next;
      

      if(tr->grouped)
	{
	  tr->constraintVector[p->number] = tr->constraintVector[q->number];
	}
      
      hookup(p,             q->back, q->z, tr->numBranches);   /* move connections to p */
      hookup(p->next,       r->back, r->z, tr->numBranches);
      hookup(p->next->next, s->back, s->z, tr->numBranches);
      
      /*if (start->number == q->number) 
	start = start->back->back;*/
      /*q->back = r->back = s->back = (nodeptr) NULL;   */
      q->back = q->next->back = q->next->next->back = (nodeptr) NULL;
    }
  else 
    {    
      p->back = p->next->back = p->next->next->back = (nodeptr) NULL;
    }
  
  assert(tr->ntips > 2);
  start = findAnyTip(tr->nodep[tr->mxtips + 1], tr->rdta->numsp);
  
  assert(isTip(start->number, tr->rdta->numsp));
  tr->rooted = FALSE;
  return  start;
}


boolean treeReadLen (FILE *fp, tree *tr, analdef *adef)
{
  nodeptr  p;
  int      i, ch, dummy; 

  for (i = 1; i <= tr->mxtips; i++) 
    tr->nodep[i]->back = (node *) NULL;

  for(i = tr->mxtips + 1; i < 2 * tr->mxtips; i++)
    {
      tr->nodep[i]->back = (nodeptr)NULL;
      tr->nodep[i]->next->back = (nodeptr)NULL;
      tr->nodep[i]->next->next->back = (nodeptr)NULL;
      tr->nodep[i]->number = i;
      tr->nodep[i]->next->number = i;
      tr->nodep[i]->next->next->number = i;
    }
    

    
  tr->start       = tr->nodep[1];
  tr->ntips       = 0;
  tr->nextnode    = tr->mxtips + 1;      
 
  for(i = 0; i < tr->numBranches; i++)
    tr->partitionSmoothed[i] = FALSE;
  
  tr->rooted      = FALSE;     

  p = tr->nodep[(tr->nextnode)++]; 
  
  while((ch = treeGetCh(fp)) != '(');
             
  if (! addElementLen(fp, tr, p, FALSE, FALSE, &dummy))                 return FALSE;
  if (! treeNeedCh(fp, ',', "in"))                return FALSE;
  if (! addElementLen(fp, tr, p->next, FALSE, FALSE, &dummy))           return FALSE;
  if (! tr->rooted) 
    {
      if ((ch = treeGetCh(fp)) == ',') 
	{ 
	  if (! addElementLen(fp, tr, p->next->next, FALSE, FALSE, &dummy)) return FALSE;	    
	}
      else 
	{                                    /*  A rooted format */
	  tr->rooted = TRUE;
	  if (ch != EOF)  (void) ungetc(ch, fp);
	}	
    }
  else 
    {
      p->next->next->back = (nodeptr) NULL;
    }
  if (! treeNeedCh(fp, ')', "in"))                return FALSE;
  (void) treeFlushLabel(fp);
  if (! treeFlushLen(fp))                         return FALSE;
 
  if (! treeNeedCh(fp, ';', "at end of"))       return FALSE;   
  
  if (tr->rooted) 
    {     
      p->next->next->back = (nodeptr) NULL;      
      tr->start = uprootTree(tr, p->next->next, FALSE);      
      if (! tr->start)                              
	{
	  printf("FATAL ERROR UPROOTING TREE\n");
	  exit(-1);	  
	}    
    }
  else 
    {
      /*
	tr->start = p->next->next->back; 
	This is start used by treeString 
      */   
      tr->start = findAnyTip(p, tr->rdta->numsp);
    }
  
  if(tr->ntips < tr->mxtips)
    {
      if(adef->computeDistance)
	{
	  printf("Error: pairwise distance computation only allows for complete, i.e., containing all taxa\n");
	  printf("bifurcating starting trees\n");
	  exit(-1);
	}     
      if(adef->classifyML)
	{	 
	  printf("RAxML classifier Algo: You provided a reference tree with %d taxa; alignmnet has %d taxa\n", tr->ntips, tr->mxtips);
	  printf("%d query taxa will be classifed under ML\n", tr->mxtips - tr->ntips);
	  classifyML(tr, adef);	  
	}
      else
	{
	  printf("You provided an incomplete starting tree %d alignmnet has %d taxa\n", tr->ntips, tr->mxtips);	  
	  makeParsimonyTreeIncomplete(tr, adef);	 		 
	}    
    }
  else
    {
      if(adef->mode == PARSIMONY_ADDITION)
	{
	  printf("Error you want to add sequences to a trees via MP stepwise addition, but \n");
	  printf("you have provided an input tree that already contains all taxa\n");
	  exit(-1);
	}
      if(adef->classifyML)
	{
	  printf("Error you want to classify query sequences into a tree via ML, but \n");
	  printf("you have provided an input tree that already contains all taxa\n");
	  exit(-1);
	}
    }
  
  
  onlyInitrav(tr, tr->start);
  
  return TRUE;
}








int treeReadTopologyOnly (FILE *fp, tree *tr, boolean readBranches, boolean forMRP, boolean readNodeLabels)
{ 
  nodeptr  p;
  int      i, ch;   
  
  int lcount = 0;

  for (i = 1; i <= tr->mxtips; i++) 
    tr->nodep[i]->back = (node *) NULL;
  
  for(i = tr->mxtips + 1; i < 2 * tr->mxtips; i++)
    {
      tr->nodep[i]->back = (nodeptr)NULL;
      tr->nodep[i]->next->back = (nodeptr)NULL;
      tr->nodep[i]->next->next->back = (nodeptr)NULL;
      tr->nodep[i]->number = i;
      tr->nodep[i]->next->number = i;
      tr->nodep[i]->next->next->number = i;
    }
    
  tr->start       = tr->nodep[tr->mxtips];
  tr->ntips       = 0;
  tr->nextnode    = tr->mxtips + 1;
  
  for(i = 0; i < tr->numBranches; i++)
    tr->partitionSmoothed[i] = FALSE;
  
  tr->rooted      = FALSE;      
  
  p = tr->nodep[(tr->nextnode)++]; 
  
  while((ch = treeGetCh(fp)) != '(');
  
  
  if (! addElementLen(fp, tr, p, readBranches, readNodeLabels, &lcount))                 exit(-1);
  if (! treeNeedCh(fp, ',', "in"))                exit(-1);
  
  
  if (! addElementLen(fp, tr, p->next, readBranches, readNodeLabels, &lcount))           exit(-1);
  if (! tr->rooted) 
    {
      if ((ch = treeGetCh(fp)) == ',') 
	{ 
	  if (! addElementLen(fp, tr, p->next->next, readBranches, readNodeLabels, &lcount)) exit(-1);	    
	}
      else 
	{                                    
	  tr->rooted = TRUE;
	  if (ch != EOF)  (void) ungetc(ch, fp);
	}	
    }
  else 
    {
      p->next->next->back = (nodeptr) NULL;
    }
  
  
  
  if (! treeNeedCh(fp, ')', "in"))                exit(-1);
  if(readNodeLabels) 
    {
      char label[64];
      int support;
      
      if(treeGetLabel (fp, label, 10))
	{
	  assert(sscanf(label, "%d\n", &support) == 1);
	  /*printf("LABEL %s\n", label);*/
	  lcount++;
	}
      assert(p->back);
      p->support = p->back->support = support;
    }
  else
    (void) treeFlushLabel(fp);
  
  if (! treeFlushLen(fp))                         exit(-1);
  
  if (! treeNeedCh(fp, ';', "at end of"))       exit(-1);
    
  if (tr->rooted) 
    {            
      assert(!readNodeLabels);
      
      p->next->next->back = (nodeptr) NULL;      
      tr->start = uprootTree(tr, p->next->next, readBranches);
      if (! tr->start)                              
	{
	  printf("FATAL ERROR UPROOTING TREE\n");
	  exit(-1);	  
	}   		
    }
  else 
    {
      /* 
	 tr->start = p->next->next->back; 
	 This is start used by treeString 
      */
      
      tr->start = findAnyTip(p, tr->rdta->numsp);
    }  
  assert(tr->mxtips == tr->rdta->numsp);
  if(tr->ntips < tr->mxtips && !forMRP)
    assert(0);
  /* 
     should be complete trees 
     this function is only called by the 
     routines in bipartitionList.c 
  */
  
  /*printf("LABELS %d\n", lcount);*/
  
  return lcount;
} 








/********************************MULTIFURCATIONS************************************************/


static boolean  addElementLenMULT (FILE *fp, tree *tr, nodeptr p, int partitionCounter)
{ 
  nodeptr  q, r, s;
  int      n, ch, fres, rn;
  double randomResolution;
  int old;
    
  tr->constraintVector[p->number] = partitionCounter; 

  if ((ch = treeGetCh(fp)) == '(') 
    {
      partCount++;
      old = partCount;       
      
      n = (tr->nextnode)++;
      if (n > 2*(tr->mxtips) - 2) 
	{
	  if (tr->rooted || n > 2*(tr->mxtips) - 1) 
	    {
	      printf("ERROR: Too many internal nodes.  Is tree rooted?\n");
	      printf("       Deepest splitting should be a trifurcation.\n");
	      return FALSE;
	    }
	  else 
	    {
	      tr->rooted = TRUE;	    
	    }
	}
      q = tr->nodep[n];
      tr->constraintVector[q->number] = partCount;
      if (! addElementLenMULT(fp, tr, q->next, old))        return FALSE;
      if (! treeNeedCh(fp, ',', "in"))             return FALSE;
      if (! addElementLenMULT(fp, tr, q->next->next, old))  return FALSE;
                 
      hookupDefault(p, q, tr->numBranches);

      while((ch = treeGetCh(fp)) == ',')
	{ 
	  n = (tr->nextnode)++;
	  if (n > 2*(tr->mxtips) - 2) 
	    {
	      if (tr->rooted || n > 2*(tr->mxtips) - 1) 
		{
		  printf("ERROR: Too many internal nodes.  Is tree rooted?\n");
		  printf("       Deepest splitting should be a trifurcation.\n");
		  return FALSE;
		}
	      else 
		{
		  tr->rooted = TRUE;
		}
	    }
	  r = tr->nodep[n];
	  tr->constraintVector[r->number] = partCount;	  

	  rn = randomInt(10000);
	  if(rn == 0) 
	    randomResolution = 0;
	  else 
	    randomResolution = ((double)rn)/10000.0;
	   	  
#ifdef DEBUG_CONSTRAINTS
	  if(1)
#endif
#ifndef DEBUG_CONSTRAINTS
	   if(randomResolution < 0.5)
#endif
	    {	    
	      s = q->next->back;	      
	      r->back = q->next;
	      q->next->back = r;	      
	      r->next->back = s;
	      s->back = r->next;	      
	      addElementLenMULT(fp, tr, r->next->next, old);	     
	    }
	  else
	    {	  
	      s = q->next->next->back;	      
	      r->back = q->next->next;
	      q->next->next->back = r;	      
	      r->next->back = s;
	      s->back = r->next;	      
	      addElementLenMULT(fp, tr, r->next->next, old);	     
	    }	    	  	  
	}       

      if(ch != ')')
	{
	  printf("Missing /) in treeReadLenMULT\n");
	  exit(-1);	        
	}
	


      (void) treeFlushLabel(fp);
    }
  else 
    {                             
      ungetc(ch, fp);
      if ((n = treeFindTipName(fp, tr)) <= 0)          return FALSE;
      q = tr->nodep[n];      
      tr->constraintVector[q->number] = partitionCounter;
#ifdef DEBUG_CONSTRAINTS
      printf("%s\n", tr->nameList[q->number]);
#endif
      if (tr->start->number > n)  tr->start = q;
      (tr->ntips)++;
      hookupDefault(p, q, tr->numBranches);
    }
  
  fres = treeFlushLen(fp);
  if(!fres) return FALSE;
    
  return TRUE;          
} 





boolean treeReadLenMULT (FILE *fp, tree *tr, analdef *adef)
{
  nodeptr  p, r, s;
  int      i, ch, n, rn;
  int partitionCounter = 0;
  double randomResolution;

  srand((unsigned int) time(NULL));
  
  for(i = 0; i < 2 * tr->mxtips; i++)
    tr->constraintVector[i] = -1;

  for (i = 1; i <= tr->mxtips; i++) 
    tr->nodep[i]->back = (node *) NULL;

  for(i = tr->mxtips + 1; i < 2 * tr->mxtips; i++)
    {
      tr->nodep[i]->back = (nodeptr)NULL;
      tr->nodep[i]->next->back = (nodeptr)NULL;
      tr->nodep[i]->next->next->back = (nodeptr)NULL;
      tr->nodep[i]->number = i;
      tr->nodep[i]->next->number = i;
      tr->nodep[i]->next->next->number = i;
    }


  tr->start       = tr->nodep[tr->mxtips];
  tr->ntips       = 0;
  tr->nextnode    = tr->mxtips + 1;
 
  for(i = 0; i < tr->numBranches; i++)
    tr->partitionSmoothed[i] = FALSE;

  tr->rooted      = FALSE;
 
  p = tr->nodep[(tr->nextnode)++]; 
  while((ch = treeGetCh(fp)) != '(');
      
  if (! addElementLenMULT(fp, tr, p, partitionCounter))                 return FALSE;
  if (! treeNeedCh(fp, ',', "in"))                return FALSE;
  if (! addElementLenMULT(fp, tr, p->next, partitionCounter))           return FALSE;
  if (! tr->rooted) 
    {
      if ((ch = treeGetCh(fp)) == ',') 
	{       
	  if (! addElementLenMULT(fp, tr, p->next->next, partitionCounter)) return FALSE;

	  while((ch = treeGetCh(fp)) == ',')
	    { 
	      n = (tr->nextnode)++;
	      assert(n <= 2*(tr->mxtips) - 2);
	
	      r = tr->nodep[n];	
	      tr->constraintVector[r->number] = partitionCounter;	   
	      
	      rn = randomInt(10000);
	      if(rn == 0) 
		randomResolution = 0;
	      else 
		randomResolution = ((double)rn)/10000.0;

#ifdef 	DEBUG_CONSTRAINTS  
	      if(1)
#endif
#ifndef DEBUG_CONSTRAINTS
	      if(randomResolution < 0.5)
#endif
		{	
		  s = p->next->next->back;		  
		  r->back = p->next->next;
		  p->next->next->back = r;		  
		  r->next->back = s;
		  s->back = r->next;		  
		  addElementLenMULT(fp, tr, r->next->next, partitionCounter);	
		}
	      else
		{
		  s = p->next->back;		  
		  r->back = p->next;
		  p->next->back = r;		  
		  r->next->back = s;
		  s->back = r->next;		  
		  addElementLenMULT(fp, tr, r->next->next, partitionCounter);
		}
	    }	  	  	      	  

	  if(ch != ')')
	    {
	      printf("Missing /) in treeReadLenMULT\n");
	      exit(-1);	        	      	      
	    }
	  else
	    ungetc(ch, fp);
	}
      else 
	{ 
	  tr->rooted = TRUE;
	  if (ch != EOF)  (void) ungetc(ch, fp);
	}       
    }
  else 
    {
      p->next->next->back = (nodeptr) NULL;
    }
    
  if (! treeNeedCh(fp, ')', "in"))                return FALSE;
  (void) treeFlushLabel(fp);
  if (! treeFlushLen(fp))                         return FALSE;
   
  if (! treeNeedCh(fp, ';', "at end of"))       return FALSE;
  

  if (tr->rooted) 
    {   
      /* 
	 printf("ROOTED\n"); 
      */

      p->next->next->back = (nodeptr) NULL;
      tr->start = uprootTree(tr, p->next->next, FALSE);
      if (! tr->start)                              return FALSE;
    }
  else 
    {
      /* 
	 printf("UNROOTED\n");
	 tr->start = p->next->next->back;
	 This is start used by treeString 
      */

      tr->start = findAnyTip(p, tr->rdta->numsp);
    }

  
  

  if(tr->ntips < tr->mxtips)
    {
#ifdef DEBUG_CONSTRAINTS
      printf("You provided an incomplete multifurcating constraint tree %d alignmnet has %d taxa\n", tr->ntips, tr->mxtips);
#endif
     
      makeParsimonyTreeIncomplete(tr, adef);
          
   }

  if(!adef->rapidBoot)
    onlyInitrav(tr, tr->start);
  return TRUE; 
}






void getStartingTree(tree *tr, analdef *adef)
{
  tr->likelihood = unlikely;
  
  if(adef->restart) 
    {	 	     	    
      INFILE = myfopen(tree_file, "r");	
                 		
      if(!adef->grouping)
	{
	  if (! treeReadLen(INFILE, tr, adef))
	    exit(-1);
	}
      else
	{
	  partCount = 0;
	  if (! treeReadLenMULT(INFILE, tr, adef))
	    exit(-1);
	}                                                               
     
      if(adef->mode == PARSIMONY_ADDITION)
	return; 

#ifdef _MULTI_GENE
      evaluateGenericInitrav(tr, tr->start); 
#endif

#ifndef _MULTI_GENE       
      evaluateGenericInitrav(tr, tr->start); 
      /*printf("%1.40f \n", tr->likelihood);*/

      treeEvaluate(tr, 1);        
      /*printf("%1.40f \n", tr->likelihood);*/
#endif      
      fclose(INFILE);
    }
  else
    { 
      assert(adef->mode != PARSIMONY_ADDITION);

      if(adef->randomStartingTree)	  
	makeRandomTree(tr, adef);       	   	 	   	  
      else
	makeParsimonyTree(tr, adef);	   	    	      		      	

      if(adef->startingTreeOnly)
	{
	  printStartingTree(tr, adef, TRUE);
	  exit(0);
	}
      else   	         
	printStartingTree(tr, adef, FALSE);     	         
             
      evaluateGenericInitrav(tr, tr->start);                                       

      /*printf("%1.40f \n", tr->likelihood);*/
     
#ifndef _MULTI_GENE     
      treeEvaluate(tr, 1); 
#endif
      
      /* printf("%1.40f \n", tr->likelihood);     */
    }         

  tr->start = tr->nodep[1];
}
