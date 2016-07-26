/* SQUID - A C function library for biological sequence analysis
 * Copyright (C) 1992-1996 Sean R. Eddy	
 *
 *    This source code is distributed under terms of the
 *    GNU General Public License. See the files COPYING 
 *    and GNULICENSE for further details.
 *
 */

#ifndef SQFUNCSH_INCLUDED
#define SQFUNCSH_INCLUDED
/* sqfuncs.h
 * 
 * Prototypes for squid library functions;
 * also makes a good reference list for what the package contains.
 *
 * Warning: squid is a slowly evolving beast. Some functions are
 * obsolete. Some functions are probably just wrong, dating to
 * a primordial era before I knew anything about what I was doing.
 * Some functions are both obsolete and wrong but still necessary
 * to get legacy code to compile.
 */

/* 
 * from aligneval.c
 */
extern float ComparePairAlignments(char *known1, char *known2, char *calc1, char *calc2);
extern float CompareRefPairAlignments(int *ref, char *known1, char *known2, char *calc1, char *calc2);
extern float CompareMultAlignments(char **kseqs, char **tseqs, int    N);
extern float CompareRefMultAlignments(int *ref, char **kseqs, char **tseqs, int    N);
extern float PairwiseIdentity(char *s1, char *s2);

/* 
 * from alignio.c
 */
extern void AllocAlignment(int nseq, int alen, char ***ret_aseq, AINFO *ainfo);
extern void FreeAlignment(char **aseqs, AINFO *ainfo);
extern void ReadAlignedFASTA(char *filename, char *env, 
			     char ***ret_aseq, AINFO *ainfo);
extern void WriteAlignedFASTA(FILE *fp, char **aseqs, AINFO *ainfo);
extern int  MakeAlignedString(char *aseq, int alen, char *ss, char **ret_s);
extern int  MakeDealignedString(char *aseq, int alen, char *ss, char **ret_s);
extern int  DealignedLength(char *aseq);
extern int  WritePairwiseAlignment(FILE *ofp, char *aseq1, char *name1, int spos1,
				   char *aseq2, char *name2, int spos2,
				   int **pam, int indent);
extern int  MingapAlignment(char **aseqs, AINFO *ainfo);
extern int  RandomAlignment(char **rseqs, SQINFO *sqinfo, int nseq, float pop, float pex,
			    char ***ret_aseqs, AINFO *ainfo);

/* from cluster.c
 */
extern int Cluster(float **mx, int N, enum clust_strategy mode, struct phylo_s **ret_tree);
extern struct phylo_s *AllocPhylo(int N);
extern void FreePhylo(struct phylo_s *tree, int N);
extern void MakeDiffMx(char **aseqs, int num, float ***ret_dmx);
extern void MakeIdentityMx(char **aseqs, int num, float ***ret_imx);
extern void PrintNewHampshireTree(FILE *fp, AINFO *ainfo, struct phylo_s *tree, int N);
extern void PrintPhylo(FILE *fp, AINFO *ainfo, struct phylo_s *tree, int N);

/* 
 * from dayhoff.c
 */
extern int  ParsePAMFile(FILE *fp, int ***ret_pam, float *ret_scale);
extern void ScalePAM(int **pam, int scale);


/* from file.c
 */
extern char *FileDirname(char *filename);
extern char *FileTail(char *file, int noextension);
extern char *FileConcat(char *dir, char *file);
extern FILE *EnvFileOpen(char *fname, char *env);
extern int   FileExists(char *filename);


/* from getopt.c
 */
extern int Getopt(int argc, char **argv, 
		  struct opt_s *opt, int nopts, char *usage,
		  int *ret_optind, char **ret_optname, char **ret_optarg);


/* from interleaved.c
 */
extern int IsInterleavedFormat(int format);
extern int ReadInterleaved(char *seqfile, 
			   int (*skip_header)(FILE *),
			   int (*parse_header)(FILE *, AINFO *),
			   int (*is_dataline)(char *, char *), 
			   char ***ret_aseqs, AINFO *ainfo);
extern int ReadAlignment(char *seqfile, int format, char ***ret_aseqs, AINFO *ainfo);

/* 
 * from msf.c
 */
extern void  WriteMSF(FILE *fp, char **aseqs, AINFO *ainfo);

/* from revcomp.c
 */
extern char *revcomp(char *comp, char *seq);

/* 
 * from selex.c
 */
extern int  ReadSELEX(char *seqfile, char ***ret_aseqs, AINFO *ret_aliinfo);
extern void WriteSELEX(FILE *fp, char **aseqs, AINFO *ainfo, int cpl);
extern int  DealignAseqs(char **aseqs, int num, char ***ret_rseqs);
extern int  IsSELEXFormat(char *filename);
extern int  TruncateNames(char **names, int N); /* OBSOLETE? */

/* 
 * from seqencode.c
 */
extern int seqcmp(char *s1, char *s2, int allow);
extern int seqncmp(char *s1, char *s2, int n, int allow);
extern int seqencode(char *codeseq,char *str);
extern int coded_revcomp(char *comp, char *seq);
extern int seqdecode(char *str, char *codeseq);
extern int seqndecode(char *str, char *codeseq, int n);

/* 
 * from sqerror.c
 */
extern void Die(char *format, ...);
extern void Warn(char *format, ...);

/* 
 * from sqio.c
 */
extern void  FreeSequence(char *seq, SQINFO *sqinfo);
extern int   SetSeqinfoString(SQINFO *sqinfo, char *sptr, int flag);
extern void  SeqinfoCopy(SQINFO *sq1, SQINFO *sq2);
extern void  ToDNA(char *seq);
extern void  ToRNA(char *seq);
extern int   ReadMultipleRseqs(char *seqfile, int fformat, char ***ret_rseqs, 
			       SQINFO **ret_sqinfo, int *ret_num);
extern SQFILE *SeqfileOpen(char *filename, int format, char *env);
extern void    SeqfilePosition(SQFILE *sqfp, long offset);
extern void    SeqfileRewind(SQFILE *sqfp);
extern void    SeqfileClose(SQFILE *sqfp);
extern int   ReadSeq(SQFILE *fp, int format, char **ret_seq, SQINFO *sqinfo);
extern int   GCGBinaryToSequence(char *seq, int len);
extern int   GCGchecksum(char *seq, int seqlen);
extern int   GCGMultchecksum(char **seqs, int nseq);
extern int   SeqfileFormat(char *filename, int  *ret_format, char *env);
extern int   WriteSeq(FILE *outf, int outfmt, char *seq, SQINFO *sqinfo);
extern int   Seqtype(char *seq);
extern char *SeqFormatString(int code);
extern GSIFILE *GSIOpen(char *gsifile);
extern int   GSIGetOffset(GSIFILE *gsi, char *key, char *sqfile, 
			  int *fmt, long *ret_offset);
extern void  GSIClose(GSIFILE *gsi);


/* from sre_ctype.c
 */
extern int sre_tolower(int c);
extern int sre_toupper(int c);

/* from sre_math.c
 */
extern float Gaussrandom(float mean, float stddev);
extern int   Linefit(float *x, float *y, int N, 
		     float *ret_a, float *ret_b, float *ret_r);
extern void  WeightedLinefit(float *x, float *y, float *var, int N,
			     float *ret_m, float *ret_b);
extern float  Gammln(float xx);
extern int    DNorm(double *vec, int n);
extern int    FNorm(float *vec, int n);
extern void   DScale(double *vec, int n, double scale);
extern void   FScale(float *vec, int n, float scale);
extern void   DSet(double *vec, int n, double value);
extern void   FSet(float *vec, int n, float value);
extern double DSum(double *vec, int n);
extern float  FSum(float *vec, int n);
extern void   DAdd(double *vec1, double *vec2, int n);
extern void   FAdd(float *vec1, float *vec2, int n);
extern void   DCopy(double *vec1, double *vec2, int n);
extern void   FCopy(float *vec1, float *vec2, int n);
extern int    DMax(double *vec, int n);
extern int    FMax(float  *vec, int n);
extern double DDot(double *vec1, double *vec2, int n);
extern float  FDot(float *vec1, float *vec2, int n);
extern float  **FMX2Alloc(int rows, int cols);
extern void     FMX2Free(float **mx);
extern double **DMX2Alloc(int rows, int cols);
extern void     DMX2Free(double **mx);
extern void     FMX2Multiply(float **A, float **B, float **C, int m, int p, int n);
extern float  sre_random(void);
extern void   sre_srandom(int seed);
extern int    DChoose(double *p, int n);
extern int    FChoose(float *p, int n);
extern double DLogSum(double *logp, int n);
extern float  FLogSum(float *logp, int n);
extern double IncompleteGamma(double a, double x);

/* from sre_string.c
 */
#ifdef NOSTR
extern char *strstr(char *s, char *subs);
#endif
extern char *Strdup(char *s);
extern void  StringChop(char *s);
extern int   Strinsert(char *s1, char c, int pos);
extern int   Strdelete(char *s1, int pos);
extern void  s2lower(char *s);
extern void  s2upper(char *s);
extern void *MallocOrDie(size_t size);
extern void *ReallocOrDie(void *p, size_t size);
extern int   Strparse(char *rexp, char *s, int ntok);
extern void  SqdClean(void);
extern void  StrShuffle(char *s1, char *s2);
extern void  StrReverse(char *s1, char *s2);
extern void  StrRegionalShuffle(char *s1, char *s2, int w);
extern char *RandomSequence(char *alphabet, float *p, int n, int len);


/* from stack.c
 */
extern struct intstack_s *InitIntStack(void);
extern void PushIntStack(struct intstack_s *stack, int data);
extern int  PopIntStack(struct intstack_s  *stack, int *ret_data);
extern void ReverseIntStack(struct intstack_s *stack);
extern int  FreeIntStack( struct intstack_s *stack );

/* 
 * from translate.c
 */
extern char *Translate(char *seq, char **code);

/* 
 * from types.c
 */
extern int  IsInt(char *s);
extern int  IsReal(char *s);
extern void Byteswap(char *swap, int nbytes);

/* 
 * from weight.c
 */
extern void GSCWeights(char **aseq, AINFO *ainfo);
extern void VoronoiWeights(char **aseq, AINFO *ainfo);
extern void BlosumWeights(char **aseq, AINFO *ainfo, float blosumlevel);
extern void FilterAlignment(char **aseq, int nseq, AINFO *ainfo, float cutoff,
			    char ***ret_anew, int *ret_nnew, 
			    AINFO **ret_newinfo);
extern void SampleAlignment(char **aseq, int nseq, AINFO *ainfo, int sample,
			    char ***ret_anew, int *ret_nnew, 
			    AINFO **ret_newinfo);

#endif /* SQFUNCSH_INCLUDED */
