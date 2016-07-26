/* SQUID - A C function library for biological sequence analysis
 * Copyright (C) 1992-1996 Sean R. Eddy	
 *
 *    This source code is distributed under terms of the
 *    GNU General Public License. See the files COPYING 
 *    and GNULICENSE for further details.
 *
 */

/* interleaved.c 
 * I/O of interleaved format multiple alignments. 
 * Modified from selex.c
 * 
 * Fri Dec  4 17:43:24 1992, SRE:
 * Reading and writing aligned sequences to/from disk files.
 * Implements a new, broader specification of SELEX format
 * and supercedes alignio.c.
 *
 * SRE, Tue Nov  9 17:40:50 1993: 
 * major revision. #= special comments and aliinfo_s optional
 * alignment info support added. Support for #=CS (consensus
 * secondary structure), #=SS (individual secondary structure),
 * #=RF (reference coordinate system), #=SQ (per-sequence header info),
 * and #=AU ("author") added.
 *
 * SRE, Mon Jan 30 14:41:49 1995:
 * #=SA side chain % surface accessibility annotation supported
 * 
 * SRE, Mon Sep 11 09:20:08 1995
 * selex.c generalized and simplified to make interleaved.c
 * 
 * SELEX format is documented in Docs/formats.tex.
 ****************************************************************************
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <memory.h>
#include <unistd.h>		/* SunOS 4.x isn't fully ANSI-compliant. */
#include "squid.h"

#ifndef SEEK_SET
#include <unistd.h>	        /* indeed, may SunOS rot in hell */
#endif

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif


static void homogenize_gapsym(char *s, char gapsym);
static int  copy_alignment_line(char *aseq, int apos, int name_rcol, 
				char *buffer, int lcol, int rcol);
static char commentsyms[] = "#%";


int
IsInterleavedFormat(int format)
{
  if (format == kSelex || format == kClustal || format == kMSF)
    return 1;
  return 0;
}


/* Function: homogenize_gapsym()
 * 
 * Purpose:  Make gap symbols homogeneous.
 */
static void 
homogenize_gapsym(char *s, char gapsym)
{
  for (; *s != '\0'; s++)
    if (isgap(*s)) *s = gapsym; 
}
      

/* Function: copy_alignment_line()
 * 
 * Purpose:  Given a line from an alignment file, and bounds lcol,rcol
 *           on what part of it may be sequence, save the alignment into
 *           aseq starting at position apos.
 *           
 *           name_rcol is set to the rightmost column this aseqs's name
 *           occupies; if name_rcol >= lcol, we have a special case in
 *           which the name intrudes into the sequence zone.
 */
static int
copy_alignment_line(char *aseq, int apos, int name_rcol, 
		    char *buffer, int lcol, int rcol)
{
  char *s1, *s2;
  int   i;
  
  s1 = aseq + apos;
  s2 = buffer;			/* be careful that buffer doesn't end before lcol! */
  for (i = 0; i < lcol; i++)
    if (*s2) s2++;

  for (i = lcol; i <= rcol; i++)
    {
      if (*s2 == '\t') {
	Warn("TAB characters will corrupt an alignment! Please remove them first.");
	return 0;
      }

      if (name_rcol >= i)                  *s1 = '.'; /* name intrusion: pad left w/ gaps */
      else if (*s2 == '\0' || *s2 == '\n') *s1 = '.'; /* short buffer: pad right w/ gaps  */
      else			           *s1 = *s2; /* normal: copy buffer into aseq */

      s1++;
      if (*s2) s2++;
    }
  return 1;
}



/* Function: is_blankline()
 * 
 * Return TRUE if line is all whitespace.
 */
static int
is_blankline(char *buffer)
{
  for (; *buffer != '\0'; buffer++) if (! isspace(*buffer)) return 0;
  return 1;
}


/* CLUSTALV and CLUSTALW support
 * 
 * skip_header(), parse_header(), is_dataline() functions for ClustalV/W
 * interleaved multiple alignment format.
 */
static int
skip_clustal(FILE *fp)
{
  char buffer[LINEBUFLEN];
  do {
    if (fgets(buffer, LINEBUFLEN, fp) == NULL) { squid_errno = SQERR_FORMAT; return 0; }
  } while (strncmp(buffer, "CLUSTAL ", 8) != 0 ||
	   strstr(buffer, "multiple sequence alignment") == NULL);
  return 1;
} 
/*ARGSUSED*/
static int
parse_clustal(FILE *fp, AINFO *ainfo)
{
  char buffer[LINEBUFLEN];
  do {
    if (fgets(buffer, LINEBUFLEN, fp) == NULL) { squid_errno = SQERR_FORMAT; return 0; }
  } while (strncmp(buffer, "CLUSTAL ", 8) != 0 ||
	   strstr(buffer, "multiple sequence alignment") == NULL);
  return 1;
}
static int 
dataline_clustal(char *buf, char *expected_name) 
{
  while (*buf && isspace(*buf)) buf++;
  if (*buf == '\0' || strchr(commentsyms, *buf) != NULL) 
    return 0;			/* blank or comment */
  if (expected_name != NULL && strncmp(buf, expected_name, strlen(expected_name) == 0))
    return 1;			/* matches expected seq name, definitely data */
  for (; *buf != '\0'; buf++)
    {				/* Clustal has no coord lines to worry about */
      if (*buf == '*' || *buf == '.' || *buf == ':') continue;  /* possible consensus line */
      if (isalnum(*buf))              return 1;   /* name or seq character   */
      if (*buf != ' ' && isgap(*buf)) return 1;   /* possible all-gap line   */
    }
  return 0;
}


/* GCG MSF support
 * 
 * skip_header(), parse_header(), is_dataline() routines for GCG MSF alignments
 */
static int
skip_MSF(FILE *fp)
{
  char buffer[LINEBUFLEN];
  do {
    if (fgets(buffer, LINEBUFLEN, fp) == NULL) { squid_errno = SQERR_FORMAT; return 0; }
  } while (strncmp(buffer, "//", 2) != 0);
  return 1;
}
static int 
parse_MSF(FILE *fp, AINFO *ainfo)
{
  char buffer[LINEBUFLEN];
  char *sptr;
  int  nseq;
  
  /* Get first dividing line. MSF format specifies ints after MSF: and Check:
   * but we don't make sure of this 
   */
  do {
    if (fgets(buffer, LINEBUFLEN, fp) == NULL) { squid_errno = SQERR_FORMAT; return 0; }
  } while (strstr(buffer, " MSF: ") == NULL ||
	   strstr(buffer, " Check: ") == NULL ||
	   strstr(buffer, " ..") == NULL);

  /* Get names, weights from header
   */
  nseq = 0;
  /*CONSTCOND*/
  while (1) 
    {
      if (fgets(buffer, LINEBUFLEN, fp) == NULL) { squid_errno = SQERR_FORMAT; return 0; }
      if (is_blankline(buffer)) continue; 
      if (strncmp(buffer, "//", 2) == 0) break;

      sptr = strtok(buffer, WHITESPACE);
      if (sptr == NULL || strcmp(sptr, "Name:") != 0 || strstr(sptr+5, "Weight:") != 0)
	{ squid_errno = SQERR_FORMAT; return 0; }
      
      if ((sptr = strtok(NULL, WHITESPACE)) == NULL) {squid_errno=SQERR_FORMAT; return 0; }
      SetSeqinfoString(&(ainfo->sqinfo[nseq]), sptr, SQINFO_NAME);

      while (sptr != NULL && strcmp(sptr, "Weight:") != 0) sptr = strtok(NULL, WHITESPACE);
      if ((sptr = strtok(NULL, WHITESPACE)) == NULL) {squid_errno=SQERR_FORMAT; return 0; }
      if (! IsReal(sptr)) { squid_errno = SQERR_FORMAT; return 0; }
      ainfo->wgt[nseq] = atof(sptr);

      nseq++;
    }

  if (nseq != ainfo->nseq) { squid_errno = SQERR_FORMAT; return 0; }
  return 1;
}
static int
dataline_MSF(char *buf, char *expected_name)
{
  while (*buf && isspace(*buf)) buf++;
  if (*buf == '\0' || strchr(commentsyms, *buf) != NULL) 
    return 0;			/* blank or comment */
  if (expected_name != NULL && strncmp(buf, expected_name, strlen(expected_name) == 0))
    return 1;			/* matches expected seq name, definitely data */
  for (; *buf != '\0'; buf++)
    {				/* MSF has coordinate lines to worry about */
      if (isspace(*buf))              continue;   /* no info from spaces     */
      if (isalpha(*buf)||isgap(*buf)) return 1;   /* has data on it          */
    }
  return 0;
}



/* Function: ReadInterleaved()
 * 
 * Purpose:  Read multiple aligned sequences from the file seqfile.
 *           Store the alignment in aseq, and the associated info
 *           in ainfo. 
 * 
 * Args:     seqfile:        name of alignment file to read.
 *           skip_header():  routine to skip the header of the file 
 *           parse_header(): routine to parse the header of the file
 *           is_dataline():  routine to determine if a line contains data
 *           ret_aseqs:      RETURN: 2D array of aligned sequences
 *           ainfo:          RETURN: optional alignment information
 *
 * Return:   Returns 1 on success. Returns 0 on failure and sets
 *           squid_errno to indicate the cause of the failure.
 */
int
ReadInterleaved(char *seqfile, 
		int (*skip_header)(FILE *),
		int (*parse_header)(FILE *, AINFO *),
		int (*is_dataline)(char *, char *),
		char ***ret_aseqs, AINFO *ainfo)
{
  FILE    *fp;                  /* ptr to opened seqfile        */
  char     buffer[LINEBUFLEN];	/* input buffer for lines       */
  char   **aseqs;               /* aligned seqs                 */
  int      nseq;		/* number of seqs read          */
  int      alen;		/* width of alignment           */
  struct block_struc {          /** alignment data for a block: */
    int lcol;			/* furthest left aligned sym    */
    int rcol;			/* furthest right aligned sym   */
  } *blocks;
  int      blocknum;		/* number of blocks             */
  char    *sptr;                /* ptr into line during parsing */
  int      currblock;		/* index for blocks             */
  int      idx;			/* counter for seqs             */
  int      currlen;
  int      inblock;		/* TRUE if in a block of data   */
  int      pos;


			/* open the file for reading; skip header*/
  fp = fopen(seqfile, "r");
  if (fp == NULL) { squid_errno = SQERR_NOFILE; return 0; }
  if (! (*skip_header) (fp)) return 0;

  /***************************************************
   * First pass across file. 
   * Determine # of seqs and width of alignment so we can alloc.
   ***************************************************/

  blocknum = 0;
  nseq     = 0; 
  alen     = 0;
  inblock  = FALSE;
  while (!feof(fp))
    {
      				/* allocate for info about this block. */
      if (blocknum == 0)
	blocks = (struct block_struc *) MallocOrDie (sizeof(struct block_struc));
      else 
	blocks = (struct block_struc *) ReallocOrDie (blocks, (blocknum+1) * sizeof(struct block_struc));
      blocks[blocknum].lcol = LINEBUFLEN+1;
      blocks[blocknum].rcol = -1;

      idx = 0;
      /*CONSTCOND*/
      while (1)			/* breaks out when blank line or EOF is hit, see below */
      {
				/* get a data line */
	do {
	  if (fgets(buffer, LINEBUFLEN, fp) == NULL) goto BREAKOUT; /* end of file  */
	  if (inblock && is_blankline(buffer))       goto BREAKOUT; /* end of block */
	} while (! (*is_dataline)(buffer, NULL));
	
	inblock = TRUE;
	if (blocknum == 0) nseq++; /* count nseq in first block */
	idx++;			   /* count # of seqs in subsequent blocks */

				/* get rcol for this block */
	for (sptr = buffer + strlen(buffer) - 1; isspace(*sptr); sptr --)
	  ;
	if (sptr - buffer > blocks[blocknum].rcol)
	  blocks[blocknum].rcol = sptr - buffer;

				/* get lcol for this block */
	if ((sptr = strtok(buffer, WHITESPACE)) == NULL) /* name */
	  { squid_errno = SQERR_FORMAT; return 0; }
	if ((sptr = strtok(NULL, WHITESPACE)) == NULL)   /* sequence */
	  { squid_errno = SQERR_FORMAT; return 0; }
	if (sptr - buffer < blocks[blocknum].lcol) 
	  blocks[blocknum].lcol = sptr - buffer;
      }
      
    BREAKOUT:			/* end of a block */
      if (inblock) 
	{
	  if (idx != nseq) { squid_errno = SQERR_FORMAT; return 0; }
	  alen += blocks[blocknum].rcol - blocks[blocknum].lcol + 1;
	  blocknum++;
	  inblock = FALSE;
	}
    }

  /***************************************************
   * Allocations; rewind file for second pass
   ***************************************************/

  AllocAlignment(nseq, alen, &aseqs, ainfo);
  rewind(fp);
  
  /***************************************************
   * Parse file header, if any. 
   * Note that we needed to know the number of seqs
   * before attempting to parse the header, because we
   * needed to allocate the alignment and assoc. info.
   ***************************************************/

  if (! (*parse_header)(fp, ainfo)) return 0;
  
 /***************************************************
  * Second pass across file: parse in the names, aseqs.
  ***************************************************/

  currlen = 0;
  for (currblock = 0 ; currblock < blocknum; currblock++)
    {
      for (idx = 0; idx < nseq; idx++)
	{
			/* get next data line */
	  do {
	    if (fgets(buffer, LINEBUFLEN, fp) == NULL) 
	      { squid_errno = SQERR_FORMAT; return 0; }
	  } while (! (*is_dataline)(buffer, ainfo->sqinfo[idx].flags & SQINFO_NAME ?
				    ainfo->sqinfo[idx].name : NULL));
	  
	  			/* find right boundary of name */
	  sptr = buffer;
	  while (*sptr && isspace(*sptr))  sptr++;
	  if (ainfo->sqinfo[idx].flags & SQINFO_NAME) 
	    while (*sptr && !isspace(*sptr)) sptr++; 
	  else			
	    {			/* first time we've seen name */
	      pos = 0;
	      while (*sptr && !isspace(*sptr) && pos < SQINFO_NAMELEN-1) {
		ainfo->sqinfo[idx].name[pos++] = *sptr;
		sptr++;
	      }
	      ainfo->sqinfo[idx].name[pos] = '\0';
	      ainfo->sqinfo[idx].flags    |= SQINFO_NAME;
	    }

				/* parse alignment line */
	  if (! copy_alignment_line(aseqs[idx], currlen, sptr - buffer,
				    buffer, blocks[currblock].lcol, 
				    blocks[currblock].rcol))
	    { squid_errno = SQERR_FORMAT; return 0; }
	}
      currlen += blocks[currblock].rcol - blocks[currblock].lcol + 1;
    }

				/* Tidy up. */
  for (idx = 0; idx < nseq; idx++)
    {
      aseqs[idx][alen] = '\0';
      homogenize_gapsym(aseqs[idx], (char) '.');
      ainfo->sqinfo[idx].len    = DealignedLength(aseqs[idx]);
      ainfo->sqinfo[idx].flags |= SQINFO_LEN;
    }
  MingapAlignment(aseqs, ainfo);

 /***************************************************
   * Garbage collection and return
   ***************************************************/
  fclose(fp);
  free(blocks);

  *ret_aseqs = aseqs;
  return 1;
}


/* Function: ReadAlignment()
 * 
 * Purpose:  Given a seqfile name and format, hand it off to appropriate
 *           parser.
 *           
 *           Currently, squid can parse alignments from the following
 *           interleaved multiple sequence alignment formats:
 *               MSF     (U. of Wisconsin GCG package MSF format)
 *               SELEX   (NeXagen/CU Boulder SELEX format)
 *               CLUSTAL (Des Higgins' CLUSTALV and CLUSTALW programs)
 *           and can sequentially read:
 *               FASTA   (aka UCSC's "a2m")
 *               
 * Return:   1 on success; 0 on failure.
 *           Returned data should be freed by caller with FreeAlignment()
 */
int
ReadAlignment(char             *seqfile, 
	      int               format,
	      char           ***ret_aseqs,
	      struct aliinfo_s *ret_ainfo)
{
  switch (format) {
  case kMSF: 
    if (! ReadInterleaved(seqfile, skip_MSF, parse_MSF, dataline_MSF, ret_aseqs, ret_ainfo)) 
      return 0;
    break;

  case kSelex:
    if (! ReadSELEX(seqfile, ret_aseqs, ret_ainfo)) 
      return 0;
    break;

  case kClustal: 
    if (! ReadInterleaved(seqfile, skip_clustal, parse_clustal, dataline_clustal, 
			  ret_aseqs, ret_ainfo)) 
      return 0;
    break;

  case kPearson:
    ReadAlignedFASTA(seqfile, NULL, ret_aseqs, ret_ainfo);
    break;

  default: squid_errno = SQERR_FORMAT; return 0;
  }
  return 1;
}
