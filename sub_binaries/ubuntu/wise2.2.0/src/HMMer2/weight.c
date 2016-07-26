/* SQUID - A C function library for biological sequence analysis
 * Copyright (C) 1992-1996 Sean R. Eddy	
 *
 *    This source code is distributed under terms of the
 *    GNU General Public License. See the files COPYING 
 *    and GNULICENSE for further details.
 *
 */

/* weight.c
 * SRE, Thu Mar  3 07:56:01 1994
 * 
 * Calculate weights for sequences in an alignment.
 */

#include <ctype.h>
#include <string.h>
#include "squid.h"

static void upweight(struct phylo_s *tree, int nseq, float *lwt, float *rwt, int node);
static void downweight(struct phylo_s *tree, int nseq, float *lwt, float *rwt, 
		       float *fwt, int node);
static float simple_distance(char *s1, char *s2);
static int    simple_diffmx(char **aseqs,int num, float ***ret_dmx);

/* Function: GSCWeights()
 * 
 * Purpose:  Use Erik's tree-based algorithm to set weights for
 *           sequences in an alignment. upweight() and downweight()
 *           are derived from Graeme Mitchison's code.
 *           
 * Args:     aseq        - array of (0..nseq-1) aligned sequences
 *           ainfo       - alignment info (including weights) 
 *                           
 * Return:   1 on success, 0 on failure.
 *           weights in ainfo are set. 
 */
void
GSCWeights(char **aseq, AINFO *ainfo)
{
  float **dmx;                 /* distance (difference) matrix */
  struct phylo_s *tree;
  float  *lwt, *rwt;           /* weight on left, right of this tree node */
  float  *fwt;                 /* final weight assigned to this node */
  int      i;
  
  /* Sanity check first
   */
  if (ainfo->nseq == 1) { ainfo->wgt[0] = 1.0; return; }

  /* I use a simple fractional difference matrix derived by
   * pairwise identity. Perhaps I should include a Poisson
   * distance correction.
   */
  MakeDiffMx(aseq, ainfo->nseq, &dmx);
  if (! Cluster(dmx, ainfo->nseq, CLUSTER_MIN, &tree))  Die("Cluster() failed");
  
  /* Allocations
   */
  if ((lwt = (float *) malloc (sizeof(float) * (2 * ainfo->nseq - 1))) == NULL ||
      (rwt = (float *) malloc (sizeof(float) * (2 * ainfo->nseq - 1))) == NULL ||
      (fwt = (float *) malloc (sizeof(float) * (2 * ainfo->nseq - 1))) == NULL)
    Die("malloc failed");
  
  /* lwt and rwt are the total branch weight to the left and
   * right of a node or sequence. They are 0..2N-2.  0..N-1 are 
   * the sequences; these have weight 0. N..2N-2 are the actual
   * tree nodes.
   */
  for (i = 0; i < ainfo->nseq; i++)
    lwt[i] = rwt[i] = 0.0;
				/* recursively calculate rwt, lwt, starting
				   at node nseq (the root) */
  upweight(tree, ainfo->nseq, lwt, rwt, ainfo->nseq);

				/* recursively distribute weight across the
				   tree */
  fwt[ainfo->nseq] = ainfo->nseq;
  downweight(tree, ainfo->nseq, lwt, rwt, fwt, ainfo->nseq);
				/* collect the weights */
  for (i = 0; i < ainfo->nseq; i++)
    ainfo->wgt[i]  = fwt[i];

  FMX2Free(dmx);
  FreePhylo(tree, ainfo->nseq);
  free(lwt); free(rwt); free(fwt);
}

static void 
upweight(struct phylo_s *tree, int nseq, float *lwt, float *rwt, int node)
{
  int ld,rd;

  ld = tree[node-nseq].left;
  if (ld >= nseq) upweight(tree, nseq, lwt, rwt, ld);
  rd = tree[node-nseq].right;
  if (rd >= nseq) upweight(tree, nseq, lwt, rwt, rd);
  lwt[node] = lwt[ld] + rwt[ld] + tree[node-nseq].lblen;
  rwt[node] = lwt[rd] + rwt[rd] + tree[node-nseq].rblen;
}


static void 
downweight(struct phylo_s *tree, int nseq, float *lwt, float *rwt, float *fwt, int node)
{
  int ld,rd;
  float lnum, rnum;

  ld = tree[node-nseq].left;
  rd = tree[node-nseq].right;
  if (lwt[node] + rwt[node] > 0.0)
    {
      fwt[ld] = fwt[node] * (lwt[node] / (lwt[node] + rwt[node]));
      fwt[rd] = fwt[node] * (rwt[node] / (lwt[node] + rwt[node]));
    }
  else
    {
      lnum = (ld >= nseq) ? tree[ld-nseq].incnum : 1.0;
      rnum = (rd >= nseq) ? tree[rd-nseq].incnum : 1.0;
      fwt[ld] = fwt[node] * lnum / (lnum + rnum);
      fwt[rd] = fwt[node] * rnum / (lnum + rnum);
    }

  if (ld >= nseq) downweight(tree, nseq, lwt, rwt, fwt, ld);
  if (rd >= nseq) downweight(tree, nseq, lwt, rwt, fwt, rd);
}




/* Function: VoronoiWeights()
 * 
 * Purpose:  Calculate weights using the scheme of Sibbald &
 *           Argos (JMB 216:813-818 1990). The scheme is
 *           slightly modified because the original algorithm
 *           actually doesn't work on gapped alignments.
 *           The sequences are assumed to be protein.
 *           
 * Args:     aseq  - array of (0..nseq-1) aligned sequences
 *           ainfo - info on the alignments, including weights  
 *
 * Return:   1 on success, 0 on failure.
 *           weights are set in ainfo
 */
void
VoronoiWeights(char **aseq, AINFO *ainfo)
{
  float **dmx;                 /* distance (difference) matrix    */
  float  *halfmin;             /* 1/2 minimum distance to other seqs */
  char   **psym;                /* symbols seen in each column     */
  int     *nsym;                /* # syms seen in each column      */
  int      symseen[27];         /* flags for observed syms         */
  char    *randseq;             /* randomly generated sequence     */
  int      acol;		/* pos in aligned columns          */
  int      idx;                 /* index in sequences              */
  int      symidx;              /* 0..25 index for symbol          */
  int      i;			/* generic counter                 */
  float   min;			/* minimum distance                */
  float   dist;		/* distance between random and real */
  float   challenge, champion; /* for resolving ties              */
  int     itscale;		/* how many iterations per seq     */
  int     iteration;           
  int     best;			/* index of nearest real sequence  */

  /* Sanity check first
   */
  if (ainfo->nseq == 1) { ainfo->wgt[0] = 1.0; return; }

  itscale = 50;

  /* Precalculate 1/2 minimum distance to other
   * sequences for each sequence
   */
  if (! simple_diffmx(aseq, ainfo->nseq, &dmx)) 
    Die("simple_diffmx() failed");
  if ((halfmin = (float *) malloc (sizeof(float) * ainfo->nseq)) == NULL)
    Die("malloc failed");
  for (idx = 0; idx < ainfo->nseq; idx++)
    {
      for (min = 1.0, i = 0; i < ainfo->nseq; i++)
	{
	  if (i == idx) continue;
	  if (dmx[idx][i] < min) min = dmx[idx][i];
	}
      halfmin[idx] = min / 2.0;
    }
  Free2DArray(dmx, ainfo->nseq);

  /* Set up the random sequence generating model.
   */
  if ((psym = (char **) malloc (ainfo->alen * sizeof(char *))) == NULL ||
      (nsym = (int *)   malloc (ainfo->alen * sizeof(int)))    == NULL)
    Die("malloc failed");
  for (acol = 0; acol < ainfo->alen; acol++)
    if ((psym[acol] = (char *) malloc (27 * sizeof(char))) == NULL)
      Die("malloc failed");

/* #ifdef ORIGINAL_SIBBALD_ALGORITHM_IS_BROKEN */
  for (acol = 0; acol < ainfo->alen; acol++)
    {
      memset(symseen, 0, sizeof(int) * 27);
      for (idx = 0; idx < ainfo->nseq; idx++)
	if (! isgap(aseq[idx][acol]))
	  {
	    if (isupper(aseq[idx][acol])) 
	      symidx = aseq[idx][acol] - 'A';
	    else
	      symidx = aseq[idx][acol] - 'a';
	    if (symidx >= 0 && symidx < 26)
	      symseen[symidx] = 1;
	  }
	else
	  symseen[26] = 1;	/* a gap */

      for (nsym[acol] = 0, i = 0; i < 26; i++)
	if (symseen[i]) 
	  {
	    psym[acol][nsym[acol]] = 'A'+i;
	    nsym[acol]++;
	  }
      if (symseen[26]) { psym[acol][nsym[acol]] = ' '; nsym[acol]++; }
    }
/* #endif ORIGINAL_SIBBALD_ALGORITHM_IS_BROKEN */

  /* Note: the original Sibbald&Argos algorithm calls for
   * bounding the sampled space using a template-like random
   * sequence generator. However, this leads to one minor
   * and one major problem. The minor problem is that
   * exceptional amino acids in a column can have a
   * significant effect by altering the amount of sampled
   * sequence space; the larger the data set, the worse
   * this problem becomes. The major problem is that 
   * there is no reasonable way to deal with gaps.
   * Gapped sequences simply inhabit a different dimensionality
   * and it's pretty painful to imagine calculating Voronoi
   * volumes when the N in your N-space is varying.
   * Note that all the examples shown by Sibbald and Argos
   * are *ungapped* examples.
   * 
   * The best way I've found to circumvent this problem is
   * just not to bound the sampled space; count gaps as
   * symbols and generate completely random sequences.
   */
#ifdef ALL_SEQUENCE_SPACE
  for (acol = 0; acol < ainfo->alen; acol++)
    {
      strcpy(psym[acol], "ACDEFGHIKLMNPQRSTVWY ");
      nsym[acol] = 21;
    }
#endif
  
  /* Sibbald and Argos algorithm:
   *   1) assign all seqs weight 0.
   *   2) generate a "random" sequence
   *   3) calculate distance to every other sequence
   *      (if we get a distance < 1/2 minimum distance
   *       to other real seqs, we can stop)
   *   4) if unique closest sequence, increment its weight 1.
   *      if multiple closest seq, choose one randomly    
   *   5) repeat 2-4 for lots of iterations
   *   6) normalize all weights to sum to nseq.
   */
  if ((randseq = (char *) malloc ((ainfo->alen+1) * sizeof(char))) == NULL)
    Die("malloc failed");

  FSet(ainfo->wgt, ainfo->nseq, 0.0);
  for (iteration = 0; iteration < itscale * ainfo->nseq; iteration++)
    {
      for (acol = 0; acol < ainfo->alen; acol++)
	randseq[acol] = (nsym[acol] == 0) ? ' ' : psym[acol][CHOOSE(nsym[acol])];
      randseq[acol] = '\0';

      champion = sre_random();
      for (min = 1.0, idx = 0; idx < ainfo->nseq; idx++)
	{
	  dist = simple_distance(aseq[idx], randseq);
	  if (dist < halfmin[idx]) 
	    { 
	      best = idx; 
	      break;      
	    } 
	  if (dist < min)          
	    { champion = sre_random(); best = idx; min = dist; } 
	  else if (dist == min)    
	    { 
	      challenge = sre_random(); 
	      if (challenge > champion)
		{ champion = challenge; best = idx; min = dist; }
	    }
	}
      ainfo->wgt[best] += 1.0;
    }

  for (idx = 0; idx < ainfo->nseq; idx++)
    ainfo->wgt[idx] = ainfo->wgt[idx] / (float) itscale;

  free(randseq);
  free(nsym);
  free(halfmin);
  Free2DArray(psym, ainfo->alen);
}


/* Function: simple_distance()
 * 
 * Purpose:  For two identical-length null-terminated strings, return
 *           the fractional difference between them. (0..1)
 *           (Gaps don't count toward anything.)
 */
static float
simple_distance(char *s1, char *s2)
{
  int diff  = 0;
  int valid = 0;

  for (; *s1 != '\0'; s1++, s2++)
    {
      if (isgap(*s1) || isgap(*s2)) continue;
      if (*s1 != *s2) diff++;
      valid++;
    }
  return (float) diff / (float) valid;
}
    
/* Function: simple_diffmx()
 * 
 * Purpose:  Given a set of flushed, aligned sequences, construct
 *           an NxN fractional difference matrix using the
 *           simple_distance rule.
 *           
 * Args:     aseqs        - flushed, aligned sequences
 *           num          - number of aseqs
 *           ret_dmx      - RETURN: difference matrix (caller must free)
 *           
 * Return:   1 on success, 0 on failure.
 */
static int
simple_diffmx(char    **aseqs,
	      int       num,
	      float ***ret_dmx)
{
  float **dmx;                 /* RETURN: distance matrix           */
  int      i,j;			/* counters over sequences           */

  /* Allocate
   */
  if ((dmx = (float **) malloc (sizeof(float *) * num)) == NULL)
    Die("malloc failed");
  for (i = 0; i < num; i++)
    if ((dmx[i] = (float *) malloc (sizeof(float) * num)) == NULL)
      Die("malloc failed");

  /* Calculate distances, symmetric matrix
   */
  for (i = 0; i < num; i++)
    for (j = i; j < num; j++)
      dmx[i][j] = dmx[j][i] = simple_distance(aseqs[i], aseqs[j]);

  /* Return
   */
  *ret_dmx = dmx;
  return 1;
}



/* Function: BlosumWeights()
 * 
 * Purpose:  Assign weights to a set of aligned sequences
 *           using the BLOSUM rule:
 *             - do single linkage clustering at some pairwise identity
 *             - in each cluster, give each sequence 1/clustsize
 *               total weight.
 *
 *           The clusters have no pairwise link >= maxid. 
 *               
 * Args:     aseqs - alignment
 *           ainfo - info on alignment, including weights
 *           maxid - fractional identity (e.g. 0.62 for BLOSUM62)
 */               
void
BlosumWeights(char **aseqs, AINFO *ainfo, float maxid)
{
  float         **dmx;          /* difference matrix */
  struct phylo_s *tree;         /* UPGMA tree        */
  float           mindiff;      /* minimum distance between clusters */
  struct intstack_s *stack;
  int             node;
  int             i;

  /* Sanity check first
   */
  if (ainfo->nseq == 1) { ainfo->wgt[0] = 1.0; return; }

  mindiff = 1.0 - maxid;
                                /* first we do a difference matrix */
  MakeDiffMx(aseqs, ainfo->nseq, &dmx);
                                /* then we build a tree */
  Cluster(dmx, ainfo->nseq, CLUSTER_MIN, &tree);

  /* Find clusters below mindiff.
   * The rule is: 
   *     -traverse the tree
   *     -if the node has dropped below mindiff, then
   *      make current node a cluster.
   */
  FSet(ainfo->wgt, ainfo->nseq, 1.0);
  stack = InitIntStack();
  PushIntStack(stack, 0);       /* push root on stack to start */
  while (PopIntStack(stack, &node))
    {
      if (tree[node].diff < mindiff)
        {                       /* we're at a cluster */
          for (i = 0; i < ainfo->nseq;  i++)
            if (tree[node].is_in[i]) 
	      ainfo->wgt[i] = 1.0 / (float) tree[node].incnum;
        }
      else                      /* we're not a cluster, keep traversing */
        {
          if (tree[node].right >= ainfo->nseq)
            PushIntStack(stack, tree[node].right - ainfo->nseq);
          if (tree[node].left >= ainfo->nseq)
            PushIntStack(stack, tree[node].left - ainfo->nseq);
        }
    }
  FreeIntStack(stack);
  FreePhylo(tree, ainfo->nseq);
  FMX2Free(dmx);
  return;
}





/* Function: FilterAlignment()
 * 
 * Purpose:  Constructs a new alignment by removing near-identical 
 *           sequences from a given alignment (where identity is 
 *           calculated *based on the alignment*).
 *           Does not affect the given alignment.
 *           Keeps earlier sequence, discards later one.
 *           
 * Args:     aseq     -- given alignment
 *           nseq     -- number of seqs in alignment
 *           ainfo    -- associated info for aseq
 *           ident    -- [0.0-1.0] above this fractional identity, remove the seq
 *           ret_anew -- RETURN: new alignment
 *           ret_nnew -- RETURN: # seqs in new alignment
 *           ret_newinfo - RETURN: info assoc with new alignment
 *                         
 * Return:   (void)
 *           ret_anew must be free'd by caller: FreeAlignment().
 */
void
FilterAlignment(char **aseq, int nseq, AINFO *ainfo, float cutoff,
		char ***ret_anew, int *ret_nnew, AINFO **ret_newinfo)
{
  char  **anew;
  int     nnew;
  AINFO  *newinfo;
  int    *list;
  float   ident;
  int     i,j, idx;
  int     remove;

				/* find which seqs to keep (list) */
				/* diff matrix; allow ragged ends */
  if ((list = (int *) malloc (sizeof(int) * nseq)) == NULL)
    Die("malloc failed");
  nnew = 0;
  for (i = 0; i < nseq; i++)
    {
      remove = FALSE;
      for (j = 0; j < nnew; j++)
	{
	  ident = PairwiseIdentity(aseq[i], aseq[list[j]]);
	  if (ident > cutoff)
	    { 
	      remove = TRUE; 
	      printf("removing %12s -- fractional identity %.2f to %s\n", 
		     ainfo->sqinfo[i].name, ident,
		     ainfo->sqinfo[list[j]].name); 
	      break; 
	    }
	}
      if (remove == FALSE) list[nnew++] = i;
    }

				/* make new alignment */
  newinfo = MallocOrDie(sizeof(AINFO));
  AllocAlignment(nnew, ainfo->alen, &anew, newinfo);
  for (idx = 0; idx < nnew; idx++)
    strcpy(anew[idx], aseq[list[idx]]);
  for (idx = 0; idx < nnew; idx++)
    SeqinfoCopy(&(newinfo->sqinfo[idx]), &(ainfo->sqinfo[list[idx]]));
				/* copy optional info */
  if (ainfo->flags & AINFO_RF) {
    newinfo->flags |= AINFO_RF;
    newinfo->rf = Strdup(ainfo->rf);
  }
  if (ainfo->flags & AINFO_CS) {
    newinfo->flags |= AINFO_CS;
    newinfo->cs = Strdup(ainfo->cs);
  }

  MingapAlignment(anew, newinfo);

  free(list);
  *ret_anew    = anew;
  *ret_nnew    = nnew;
  *ret_newinfo = newinfo;
}


/* Function: SampleAlignment()
 * 
 * Purpose:  Constructs a new, smaller alignment by sampling a given
 *           number of sequences at random. Does not change the
 *           alignment nor the order of the sequences.
 *           
 *           Not a weighting method, but this is as good a place
 *           as any, since the code is borrowed from FilterAlignment().
 *           
 *           If you ask for a sample that is larger than nseqs,
 *           it silently returns the original alignment.
 *           
 * Args:     aseq     -- given alignment
 *           nseq     -- number of seqs in alignment
 *           ainfo    -- associated info for aseq
 *           sample   -- number of seqs to sample
 *           ret_anew -- RETURN: new alignment
 *           ret_nnew -- RETURN: # seqs in new alignment
 *           ret_newinfo - RETURN: info assoc with new alignment
 *                         
 * Return:   (void)
 *           ret_anew must be free'd by caller: FreeAlignment().
 */
void
SampleAlignment(char **aseq, int nseq, AINFO *ainfo, int sample,
		char ***ret_anew, int *ret_nnew, AINFO **ret_newinfo)
{
  char  **anew;
  int     nnew;
  AINFO  *newinfo;
  int    *list;                 /* array for random selection w/o replace */
  int    *useme;                /* array of flags 0..nseq-1: TRUE to use */
  int     i, idx;
  int     len;

  /* Allocations
   */
  list  = (int *) MallocOrDie (sizeof(int) * nseq);
  useme = (int *) MallocOrDie (sizeof(int) * nseq);
  for (i = 0; i < nseq; i++)
    {
      list[i]  = i;
      useme[i] = FALSE;
    }

  /* Sanity check.
   */
  if (sample >= nseq)
    {				/* use everything and don't complain */
      for (i = 0; i < nseq; i++)
	useme[i] = TRUE;
      nnew = nseq;
    }
  else
    {				/* random selection w/o replacement */
      for (len = nseq, i = 0; i < sample; i++)
	{
	  idx = CHOOSE(len);
	  printf("chose %d: %s\n", list[idx], ainfo->sqinfo[list[idx]].name);
	  useme[list[idx]] = TRUE;
	  list[idx] = list[--len];
	}
      nnew = sample;
    }

				/* make new alignment */
  if ((anew = (char **) malloc (sizeof(char *) * nnew)) == NULL ||
      (newinfo = (AINFO *) malloc (sizeof(AINFO))) == NULL)
    Die("malloc failed");
  for (idx = 0; idx < nnew; idx++)
    if ((anew[idx] = (char *) malloc (ainfo->alen + 1)) == NULL)
      Die("malloc failed");
  if ((newinfo->sqinfo = (SQINFO *) malloc (sizeof(SQINFO) * nnew)) == NULL)
    Die("malloc failed");

  for (idx = i = 0; idx < nseq; idx++)
    if (useme[idx])
      {
	anew[i] = Strdup(aseq[idx]);
	SeqinfoCopy(&(newinfo->sqinfo[i]), &(ainfo->sqinfo[idx]));
	i++;
      }

  newinfo->nseq  = nnew;
  newinfo->alen  = ainfo->alen;
				/* copy optional info */
  if (ainfo->flags & AINFO_RF) {
    newinfo->flags |= AINFO_RF;
    newinfo->rf = Strdup(ainfo->rf);
  }
  if (ainfo->flags & AINFO_CS) {
    newinfo->flags |= AINFO_CS;
    newinfo->cs = Strdup(ainfo->cs);
  }

  MingapAlignment(anew, newinfo);

  free(list);
  free(useme);
  *ret_anew    = anew;
  *ret_nnew    = nnew;
  *ret_newinfo = newinfo;
}
