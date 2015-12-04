/* SQUID - A C function library for biological sequence analysis
 * Copyright (C) 1992-1996 Sean R. Eddy	
 *
 *    This source code is distributed under terms of the
 *    GNU General Public License. See the files COPYING 
 *    and GNULICENSE for further details.
 *
 */

#ifndef SQUIDH_INCLUDED
#define SQUIDH_INCLUDED

/* squid.h
 * last modified Sun Aug 15 12:05:58 1993
 * 
 * Header file for my library of sequence functions.
 * 
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/types.h>		/* only for ntohl() and friends. */
#include <netinet/in.h>		/* only for ntohl() and friends. */

/* Library version info is made available as a global to
 * any interested program. These are defined in iupac.c
 * with the other globals.
 */
extern char squid_version[];	/* version number  */
extern char squid_date[];	/* date of release */
extern int  squid_errno;	/* error codes     */

/****************************************************
 * Error codes returned by squid library functions (squid_errno)
 ****************************************************/

#define SQERR_OK        0	/* no error                     */
#define SQERR_UNKNOWN   1       /* generic error, unidentified  */
#define SQERR_NODATA    2	/* unexpectedly NULL stream     */
#define SQERR_MEM       3	/* malloc or realloc failed     */
#define SQERR_NOFILE    4	/* file not found               */
#define SQERR_FORMAT    5	/* file format not recognized   */
#define SQERR_PARAMETER 6	/* bad parameter passed to func */
#define SQERR_DIVZERO   7	/* error in sre_math.c          */

/****************************************************
 * Single sequence information
 ****************************************************/ 
#define SQINFO_NAMELEN 64
#define SQINFO_DESCLEN 128

struct seqinfo_s {
  int      flags;               /* what extra data are available         */
  char     name[SQINFO_NAMELEN];/* up to 63 characters of name           */
  char     id[SQINFO_NAMELEN];	/* up to 63 char of database identifier  */
  char     acc[SQINFO_NAMELEN]; /* up to 63 char of database accession # */
  char     desc[SQINFO_DESCLEN];/* up to 127 char of description         */
  int      len;                 /* length of this seq                    */
  int      start;		/* (1..len) start position on source seq */
  int      stop;                /* (1..len) end position on source seq   */
  int      olen;                /* original length of source seq         */
  int      type;                /* kRNA, kDNA, kAmino, or kOther         */
  char    *ss;                  /* 0..len-1 secondary structure string   */
  char    *sa;			/* 0..len-1 % side chain surface access. */
};
typedef struct seqinfo_s SQINFO;

#define SQINFO_NAME  (1 << 0)
#define SQINFO_ID    (1 << 1)
#define SQINFO_ACC   (1 << 2)
#define SQINFO_DESC  (1 << 3)
#define SQINFO_START (1 << 4)
#define SQINFO_STOP  (1 << 5)
#define SQINFO_LEN   (1 << 6)
#define SQINFO_TYPE  (1 << 7)
#define SQINFO_OLEN  (1 << 8)
#define SQINFO_SS    (1 << 9)
#define SQINFO_SA    (1 << 10)

/****************************************************
 * Database indexing (GSI index file format)
 ****************************************************/
/* A GSI (generic sequence index) file is composed of
 * recnum + nfiles + 1 records. Each record contains
 * three fields; key, file number, and disk offset.
 * Record 0 contains:
 *   [ "GSI" ]  [ nfiles ]  [ recnum ]
 * Records 1..nfiles map file names to file numbers, and contain:
 *   [ filename ] [ file number, 1..nfiles ] [ 0 (unused) ]
 * Records nfiles+1 to recnum+nfiles+1 provide disk offset
 * and file number indices for every key:
 *   [ key ] [ file number ] [ offset]
 */
typedef unsigned int   sqd_uint32; /* 32 bit integer. */
typedef unsigned short sqd_uint16; /* 16 bit integer. */ 
struct gsi_s {
  FILE        *gsifp;		/* open GSI index file            */
  sqd_uint16   nfiles;		/* number of files = 16 bit int   */
  sqd_uint32   recnum;		/* number of records = 32 bit int */
};
typedef struct gsi_s GSIFILE;

#define GSI_KEYSIZE    32         /* keys are 32 bytes long */  
#define GSI_RECSIZE    38	  /* 32 + 2 + 4 bytes       */
/* if ntohl() and friends are not available, you
 * can slip replacements in by providing sre_ntohl()
 * functions. (i.e., there is a possible portability problem here.)
 */
#define sre_ntohl(x)  ntohl(x); 
#define sre_ntohs(x)  ntohs(x);

/****************************************************
 * Sequence alphabet: see also iupac.c
 ****************************************************/
				/* IUPAC symbols defined globally in iupac.c */
struct iupactype {
  char       sym;		/* character representation */
  char       symcomp;           /* complement (regular char */
  char       code;		/* my binary rep */
  char       comp;              /* binary encoded complement */
};
extern struct iupactype iupac[];
#define IUPACSYMNUM 17

extern char    *stdcode1[];	/* 1-letter amino acid translation code */
extern char    *stdcode3[];	/* 3-letter amino acid translation code */
extern float   aafq[];		/* amino acid occurrence frequencies    */
extern char     aa_alphabet[];  /* amino acid alphabet                  */
extern int      aa_index[];     /* convert 0..19 indices to 0..26       */

				/* valid symbols in IUPAC code */
#define NUCLEOTIDES    "ACGTUNRYMKSWHBVDacgtunrymkswhbvd"
#define AMINO_ALPHABET "ACDEFGHIKLMNPQRSTVWY"
#define DNA_ALPHABET   "ACGT"
#define RNA_ALPHABET   "ACGU"
#define WHITESPACE     " \t\n"

#define isgap(c) ((c) == ' ' || (c) == '.' || (c) == '_' || (c) == '-')

/****************************************************
 * Alignment information
 ****************************************************/

/* Structure: aliinfo_s
 * 
 * Purpose:   Optional information returned from an alignment file.
 * 
 *            flags: always used. Flags for which info is valid/alloced.
 *       
 *            alen: mandatory. Alignments are always flushed right
 *                  with gaps so that all aseqs are the same length, alen.
 *                  Available for all alignment formats.
 *
 *            nseq: mandatory. Aligned seqs are indexed 0..nseq-1. 
 *                  
 *            wgt:  0..nseq-1 vector of sequence weights. Mandatory.
 *                  If not explicitly set, weights are initialized to 1.0.
 *
 *            cs:   0..alen-1, just like the alignment. Contains single-letter
 *                  secondary structure codes for consensus structure; "<>^+"
 *                  for RNA, "EHL." for protein. May be NULL if unavailable
 *                  from seqfile. Only available for SELEX format files.
 *                  
 *            rf:   0..alen-1, just like the alignment. rf is an arbitrary string
 *                  of characters, used for annotating columns. Blanks are
 *                  interpreted as non-canonical columns and anything else is
 *                  considered canonical. Only available from SELEX files.
 *                  
 *            sqinfo: mandatory. Array of 0..nseq-1 
 *                  per-sequence information structures, carrying
 *                  name, id, accession, coords.
 *                  
 */
struct aliinfo_s {		
  int               flags;      /* flags for what info is valid             */
  int               alen;	/* length of alignment (columns)            */
  int               nseq;       /* number of seqs in alignment              */
  char              au[64];	/* "author" information                     */
  float            *wgt;	/* sequence weights [0..nseq-1]             */
  char             *cs;         /* consensus secondary structure string     */
  char             *rf;         /* reference coordinate system              */
  struct seqinfo_s *sqinfo;     /* name, id, coord info for each sequence   */
};
typedef struct aliinfo_s AINFO;

#define AINFO_AUTH    (1 << 0)
#define AINFO_CS      (1 << 1)
#define AINFO_RF      (1 << 2)


/****************************************************
 * Sequence i/o: originally from Don Gilbert's readseq 
 ****************************************************/
	/* buffer size for reading in lines from sequence files*/
#define LINEBUFLEN  4096

/* sequence types parsed by Seqtype()                          */
/* note that these must match hmmAMINO and hmmNUCLEIC in HMMER */
#define kOtherSeq   0
#define kDNA        1
#define kRNA        2
#define kAmino      3

/* Sequence file formats recognized */
#define kUnknown        0   /* format not determinable */
#define kIG             1
#define kGenBank        2
#define kA2M            3	/* aligned FASTA (A2M is UCSC's terminology) */
#define kEMBL           4
#define kGCG            5
#define kStrider        6
#define kPearson        7
#define kZuker          8
#define kIdraw          9	/* idraw-style PostScript (write only)       */
#define kSelex          10	/* my flat text alignment format             */
#define kMSF		11	/* GCG MSF multiple alignment format         */
#define kPIR            12      /* PIR-CODATA format                         */
#define kRaw            13      /* unformatted, raw sequence (output only)   */
#define kSquid          14	/* my sequence database format               */
				/* 15 was kXPearson, extended Pearson        */
#define kGCGdata        16	/* GCG data library format                   */
#define kClustal        17	/* Clustal V or W multiple alignment format  */

#define kMinFormat      1	/* SRE: kUnknown doesn't count */
#define kMaxFormat      17
#define kNumFormats     (kMaxFormat + 1)
#define kNoformat       -1	/* format not tested */

struct ReadSeqVars {
  FILE   *f;
  char    sbuffer[LINEBUFLEN];	/* current line we're working on */

  char   *seq;
  SQINFO *sqinfo;	/* name, id, etc. */
  char   *sp;
  int     seqlen;
  int     maxseq;

  int     longline;		/* TRUE if last fgets() didn't have a \n */
  int     dash_equals_n;        /* a hack - affects EMBL reading, to deal with Staden
				   experiment files */
  int     do_gzip;		/* TRUE if f is a pipe from gzip -dc */
  int     do_stdin;		/* TRUE if f is stdin */
  int     format;		/* format of seqfile we're reading */


  /* The ali_ section of the structure is a hack for sequential access
   * of multiple alignment files: we read the whole alignment in,
   * and then copy it one sequence at a time into seq and sqinfo.
   * It is active if ali_num > 0. Because we keep it in the SQFILE structure,
   * ReadSeq() and friends are always reentrant for multiple seqfiles.
   */
  char  **ali_aseqs;
  AINFO   ali_ainfo;
  int     ali_curridx;		/* next index to copy and return [0..ali_num-1]*/
};
typedef struct ReadSeqVars SQFILE;


/****************************************************
 * Cluster analysis and phylogenetic tree support
 ****************************************************/ 

/* struct phylo_s - a phylogenetic tree
 *                     
 * For N sequences, there will generally be an array of 0..N-2
 * phylo_s structures representing the nodes of a tree.
 * [0] is the root. The indexes of left and
 * right children are somewhat confusing so be careful. The
 * indexes can have values of 0..2N-2. If they are 0..N-1, they 
 * represent pointers to individual sequences. If they are
 * >= N, they represent pointers to a phylo_s structure
 * at (index - N).
 */
struct phylo_s {
  int    parent;                /* index of parent, N..2N-2, or -1 for root */
  int    left;			/* index of one of the branches, 0..2N-2 */
  int    right;			/* index of other branch, 0..2N-2        */
  float  diff;			/* difference score between seqs         */
  float  lblen;      		/* left branch length                    */
  float  rblen;                 /* right branch length                   */
  char  *is_in;                 /* 0..N-1 flag array, 1 if seq included  */
  int    incnum;                /* number of seqs included at this node  */
};


/* Strategies for cluster analysis; cluster by mean distance,
 * minimum distance, or maximum distance.
 */
enum clust_strategy { CLUSTER_MEAN, CLUSTER_MAX, CLUSTER_MIN };

/****************************************************
 * Generic data structure support
 ****************************************************/

/* a struct intstack_s implements a pushdown stack for storing
 * single integers.
 */
struct intstack_s {
  int                data;
  struct intstack_s *nxt;
};

/****************************************************
 * Binary nucleotide alphabet support
 ****************************************************/

/* Binary encoding of the IUPAC code for nucleotides
 * 
 *    four-bit "word", permitting rapid degenerate matching
 *         A  C  G  T/U
 *         0  0  1  0
 */
#define NTA 8
#define NTC 4
#define NTG 2
#define NTT 1
#define NTU 1
#define NTN 15			/* A|C|G|T */
#define NTR 10			/* A|G */
#define NTY 5			/* C|T */
#define NTM 12			/* A|C */
#define NTK 3			/* G|T */
#define NTS 6			/* C|G */
#define NTW 9			/* A|T */
#define NTH 13			/* A|C|T */
#define NTB 7			/* C|G|T */
#define NTV 14			/* A|C|G */
#define NTD 11			/* A|G|T */
#define NTGAP 16		/* GAP */
#define NTEND 0			/* null string terminator */

/* ntmatch(): bitwise comparison of two nuc's 
 * note that it's sensitive to the order;
 * probe may be degenerate but target should not be 
 */
#define ntmatch(probe, target)  ((probe & target) == target)

/****************************************************
 * Support for a portable, flexible Getopt()
 ****************************************************/

/* Structure: opt_s
 * 
 * Structure for declaring options to a main().
 */
struct opt_s {
  char *name;			/* name of option, e.g. "--option1" or "-o" */
  int   single;			/* TRUE if a single letter option           */
  int   argtype;		/* for typechecking, e.g. sqdARG_INT        */
};
				/* acceptable argtype's...           */
#define sqdARG_NONE   0		/* no argument                       */
#define sqdARG_INT    1		/* something that atoi() can grok    */
#define sqdARG_FLOAT  2		/* something that atof() can grok    */
#define sqdARG_CHAR   3		/* require single character or digit */
#define sqdARG_STRING 4		/* anything goes                     */

/****************************************************
 * Support for convenient Perl-y regexp matching
 *
 * Strparse() defines and manages these. 
 * sqd_parse[0] contains the substring that matched the pattern.
 * sqd_parse[1-9] contain substrings matched with ()'s.
 ****************************************************/

extern char *sqd_parse[10];

/****************************************************
 * Miscellaneous macros and defines
 ****************************************************/

#define CHOOSE(a)     ((int) (sre_random() * (a)))
				/* must declare swapfoo to use SWAP() */
#define SWAP(a,b) {swapfoo = b; b = a; a = swapfoo;}

#define Free2DArray(ptr, n) \
{ int fooidx;\
  if (ptr != NULL) { \
    for (fooidx = 0; fooidx < (n); fooidx++) if (ptr[fooidx] != NULL) free(ptr[fooidx]);\
    free(ptr);\
  } }    

#define ScalarsEqual(a,b) (fabs((a)-(b)) < 1e-7)

#ifndef MIN
#define MIN(a,b)         ((a<b)?a:b)
#endif
#ifndef MAX
#define MAX(a,b)         ((a>b)?a:b)
#endif

/* For convenience and (one hopes) clarity in boolean tests:
 */
#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE 
#define FALSE 0
#endif

/* Somewhere, there is a universe in which Unix vendors comply
 * with the ANSI C standard. Unfortunately, it is not ours:
 */
#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS 0
#endif
#ifndef EXIT_FAILURE
#define EXIT_FAILURE 1
#endif

#include "sqfuncs.h"		/* squid function declarations */
#endif /* SQUIDH_INCLUDED */
