/*  
* $Revision: 1.1.1.1 $  
* $Id: dbac.c,v 1.1.1.1 2010/03/18 22:55:58 soo Exp $  
*/  
  
typedef int BOOL ;  
#define TRUE 1  
#define FALSE 0  
/*WK  #define MAXD 1.0E300 */  
  
/*#define DOPROB*/  
  
#include "dba.h"  
#include "dbadisplay.h"  
#include "version.h"  
  
/*  
 * Dbb alignment. Idea by R. Durbin  
 */  
static Probability MATCH65 = 0.65 ;  
static Probability MATCH75 = 0.75 ;  
static Probability MATCH85 = 0.85 ;  
static Probability MATCH95 = 0.95 ;  
static Probability MATCHS = 0.8 ;  
static Probability GAP = 0.05 ;  
static Probability BLOCKOPEN = 0.01 ;  
static Probability UNMATCHED_PEN = 0.99;  
  
int kbyte = 25000;  
  
  
  
void prettyAlnBlock (AlnBlock *alb, Sequence *one, Sequence *two)  
{  
  AlnColumn *alc ;  
  AlnUnit *alu ;  
  
  for (alc = alb->start ; alc ; alc = alc->next)   
    { alu = alc->alu[0] ;  
      if (!alu->text_label) alu->text_label = "empty_label" ;  
      printf ("%2s ", alu->text_label) ;  
      if (alu->end > alu->start)   
 printf ("%5d %c", alu->end, one->seq[alu->end]) ;  
      else   
        printf ("       ") ;  
      alu = alc->alu[1] ;  
      if (alu->end > alu->start)   
 printf (" %c %5d", two->seq[alu->end], alu->end) ;  
      else   
        printf ("        ") ;  
      if (!alu->text_label) alu->text_label = "empty_label" ;  
      printf (" %2s", alu->text_label) ;  
      printf ("\n") ;  
   }  
}  
  
void alnBlockSummary (AlnBlock *alb, Sequence *x, Sequence *y)  
{  
  AlnColumn *alc ;  
  AlnUnit *xa, *ya ;  
  BOOL inBlock = FALSE ;  
 /* int type = 0; */  
/*   char * type;  
 */  /* type = malloc(8); */  
        char t[8];  
  char * type = &t[0];  
  int nMatch, nMismatch, nGap, x1, y1 ;  
  
  for (alc = alb->start ; alc ; alc = alc->next)   
    { xa = alc->alu[0] ;  
      ya = alc->alu[1] ;  
      if (*xa->text_label == 'M')  
 { if (!inBlock)  
     { inBlock = TRUE ;  
    if (xa->text_label[2] == '6') {  
     /* type = 65; */  
     /* type = MATCH65 * 100; */  
     type = "MATCH_A";  
    } else if (xa->text_label[2] == '7') {  
     /* type = 75; */  
     /* type = MATCH75 * 100; */  
     type = "MATCH_B";  
    } else if (xa->text_label[2] == '8') {  
     /* type = 85; */  
     /* type = MATCH85 * 100; */  
     type = "MATCH_C";  
    } else if (xa->text_label[2] == '9') {  
     /* type = 95; */  
     /* type = MATCH95 * 100; */  
     type = "MATCH_D";  
    }  
       nMatch = nMismatch = nGap = 0 ;  
       x1 = xa->end ;   
       y1 = ya->end ;  
     }  
   if (xa->text_label[1] == 'I' || ya->text_label[1] == 'I')  
     ++nGap ;  
   else if (x->seq[xa->end] == y->seq[ya->end])  
     ++nMatch ;  
   else  
     ++nMismatch ;  
 }  
      else if (inBlock)  
 {   
/*   printf ("%6d %6d   %6d %6d   %2d%%  %d gaps, %d\n",  
    x1, xa->start, y1, ya->start,  
    (int)((100.0*nMatch)/(nMatch + nMismatch) + 0.5),  
    nGap,  
    nMatch*MATCH + nMismatch*MISMATCH + nGap*GAP +BLOCKOPEN) ; */  
     
/*    printf ("%6d %6d   %6d %6d   %2d%%  %d gaps, %s\n",  
    x1, xa->start, y1, ya->start,  
    (int)((100.0*nMatch)/(nMatch + nMismatch) + 0.5),  
    nGap,  
   type) ; */     
            
          printf ("%6d %6d   %6d %6d   %2d%%  %d gaps, %s\n",  
    x1 + 1, xa->start + 1, y1 +1 , ya->start + 1,  
    (int)((100.0*nMatch)/(nMatch + nMismatch) + 0.5),  
    nGap,  
   type) ;  
      
  
   inBlock = FALSE ;  
 }  
    }  
    free(type);  
}  
  
  
/*  
 * Produces phylogenetic footprinting format (gapped fasta) output, with a gff line as the comment line like this:  
 * ><seqname>  DBA  Block  <start> <stop>  <bitscore>  +  .  block <block#>;id <id val>;gaps <total gaps>;type <match type>  
 */  
  
void showPFF (AlnBlock *alb, Sequence *x, Sequence *y)  
{  
  AlnColumn *alc ;  
  AlnUnit *xa, *ya ;  
  BOOL inBlock = FALSE ;  
 /* int type = 0; */  
/*  char * type;  
        type = malloc(8); */  
        char t[8];  
  char * type = &t[0];  
   
  int nMatch, nMismatch, nGap, x1, y1, xcount, ycount ;  
  int blockCount = 0;  /* used for the block numbering in the comment gff line */  
/*WK   double score; */  
  double tmpbitscore;  
  double currMatch;  
  double s = 1-2*GAP-BLOCKOPEN;  /* used for calculating score */  
    
  char xtmp[2] = " ";       /*   The Jareborg way of C programming:  */  
  char ytmp[2] = " ";       /*   these are used for converting the   */  
  char * xstr = &xtmp[0];   /*   nucleotide chars to strings when    */  
  char * ystr = &ytmp[0];   /*   building the gapped sequence        */  
    
  char * xSeq;         /* holder for the gapped x sequence */  
  char * ySeq;         /* holder for the gapped y sequence */  
  xSeq = malloc(2 * x->len); /* allocate enough memory to hold the gapped x sequence */   
  ySeq = malloc(2 * y->len); /* allocate enough memory to hold the gapped y sequence */   
    
  if (xSeq == NULL) {  
      fprintf(stderr, "Can't allocate memory for sequence x\n");  
      exit(1);  
  }  
  if (ySeq == NULL) {  
      fprintf(stderr, "Can't allocate memory for sequence y\n");  
      exit(1);  
  }  
  /* make sure the sequence strings are empty before we start */  
  strcpy (xSeq, "");   
  strcpy (ySeq, "");  
    
  for (alc = alb->start ; alc ; alc = alc->next) {   
      xa = alc->alu[0] ;  
      ya = alc->alu[1] ;  
      if (*xa->text_label == 'M') {   
        if (!inBlock) {   
            inBlock = TRUE ;  
            blockCount++;  
/*WK             score = (((1-UNMATCHED_PEN)/(UNMATCHED_PEN*UNMATCHED_PEN)) * (BLOCKOPEN/UNMATCHED_PEN)) / 4; */  
/*  
 * 010613 Nic: log2 function not available on all platforms  
 */  
/*             tmpbitscore = log2((((1-UNMATCHED_PEN)/(UNMATCHED_PEN*UNMATCHED_PEN)) * (BLOCKOPEN/UNMATCHED_PEN)) / 4); */  
            tmpbitscore = (log((((1-UNMATCHED_PEN)/(UNMATCHED_PEN*UNMATCHED_PEN)) * (BLOCKOPEN/UNMATCHED_PEN)) / 4))/log(2);  
  
     if (xa->text_label[2] == '6') {  
      /* type = 65; */  
      /* type = MATCH65 * 100; */  
      type = "MATCH_A";  
                    currMatch = MATCH65;  
     } else if (xa->text_label[2] == '7') {  
      /* type = 75; */  
      /* type = MATCH75 * 100; */  
      type = "MATCH_B";  
                    currMatch = MATCH75;  
     } else if (xa->text_label[2] == '8') {  
      /* type = 85; */  
      /* type = MATCH85 * 100; */  
      type = "MATCH_C";  
                    currMatch = MATCH85;  
     } else if (xa->text_label[2] == '9') {  
      /* type = 95; */  
      /* type = MATCH95 * 100; */  
      type = "MATCH_D";  
                    currMatch = MATCH95;  
     }  
     nMatch = nMismatch = nGap = 0 ;  
     x1 = xa->end ;   
     y1 = ya->end ;  
 }  
          
 if (xa->text_label[1] == 'I') {  
   ++nGap ;  
          xtmp[0] = x->seq[xa->end];  
          strcat(xSeq, xstr);  
          strcat(ySeq, ".");  
/*WK           score *= GAP/UNMATCHED_PEN; */  
  
/*  
 * 010613 Nic: log2 function not available on all platforms  
 */  
/*           tmpbitscore += log2(GAP/UNMATCHED_PEN);  
 */  
          tmpbitscore += (log(GAP/UNMATCHED_PEN))/log(2);  
            
        } else if (ya->text_label[1] == 'I') {  
          ++nGap;  
          strcat(xSeq, ".");  
          ytmp[0] = y->seq[ya->end];  
          strcat(ySeq, ystr);  
/*WK           score *= GAP/UNMATCHED_PEN; */  
/*  
 * 010613 Nic: log2 function not available on all platforms  
 */  
/*           tmpbitscore += log2(GAP/UNMATCHED_PEN);  
 */  
          tmpbitscore += (log(GAP/UNMATCHED_PEN))/log(2);  
        } else if (x->seq[xa->end] == y->seq[ya->end]) {  
   ++nMatch ;  
          /* printf("nucl %c\n", x->seq[xa->end]); */  
          xtmp[0] = x->seq[xa->end];  
          strcat(xSeq, xstr);  
          ytmp[0] = y->seq[ya->end];  
          strcat(ySeq, ystr);  
   /*WK   if(score<=MAXD)  
     score *= (s/(UNMATCHED_PEN*UNMATCHED_PEN)) * (currMatch/0.25);  
   */  
  
/*  
 * 010613 Nic: log2 function not available on all platforms  
 */  
/*           tmpbitscore += log2( (s/(UNMATCHED_PEN*UNMATCHED_PEN))   
                               *(currMatch/0.25));  
 */            
          tmpbitscore += (log( (s/(UNMATCHED_PEN*UNMATCHED_PEN))   
                               *(currMatch/0.25)))/log(2);  
  
        } else {  
   ++nMismatch ;  
          xtmp[0] = x->seq[xa->end];  
          strcat(xSeq, xstr);  
          ytmp[0] = y->seq[ya->end];  
          strcat(ySeq, ystr);  
/*WK           score *= (s/(UNMATCHED_PEN*UNMATCHED_PEN)) * ( ((1-currMatch)/3)/0.25 ); */  
/*  
 * 010613 Nic: log2 function not available on all platforms  
 */  
/*           tmpbitscore += log2( (s/(UNMATCHED_PEN*UNMATCHED_PEN))   
                               *( ((1-currMatch)/3)/0.25 ));  
 */            
          tmpbitscore += (log( (s/(UNMATCHED_PEN*UNMATCHED_PEN))   
                               *( ((1-currMatch)/3)/0.25 )))/log(2);  
        }  
      } else if (inBlock) {   
/*WK  double bitscore = log2(score); */  
  double bitscore = tmpbitscore;  
          
 /* print comment gff line for x */  
/*         printf (">%s\tDBA\tBlock\t%d\t%d\t%.2f\t+\t.\tblock %2d; id %d; gaps %d; type %s\n",  
                x->name,   
  x1, xa->start,  
                bitscore,   
                blockCount,   
  (int)((100.0*nMatch)/(nMatch + nMismatch) + 0.5),  
  nGap, type) ;  
  
 */  
         printf (">%s\tDBA\tBlock\t%d\t%d\t%.2f\t+\t.\tblock %2d; id %d; gaps %d; type %s\n",  
                x->name,   
  x1 + 1, xa->start + 1,  
                bitscore,   
                blockCount,   
  (int)((100.0*nMatch)/(nMatch + nMismatch) + 0.5),  
  nGap, type) ;  
  
        /*print sequence for x */  
        xcount = 0;  
        while(*xSeq != '\0') {  
          putchar(*xSeq);  
          ++xSeq;  
          if (++xcount%60 == 0) putchar ('\n');  
        }  
        putchar ('\n');  
          
 /* print comment gff line for y */  
/*         printf (">%s\tDBA\tBlock\t%d\t%d\t%.2f\t+\t.\tblock %2d; id %d; gaps %d; type %s\n",  
                y->name,   
  y1, ya->start,  
                bitscore,   
                blockCount,  
  (int)((100.0*nMatch)/(nMatch + nMismatch) + 0.5),  
  nGap, type) ; */  
  
  
  
         printf (">%s\tDBA\tBlock\t%d\t%d\t%.2f\t+\t.\tblock %2d; id %d; gaps %d; type %s\n",  
                y->name,   
  y1 +1 , ya->start + 1,  
                bitscore,  
                blockCount,  
  (int)((100.0*nMatch)/(nMatch + nMismatch) + 0.5),  
  nGap, type) ;  
  
        /*print sequence for x */  
        ycount = 0;  
        while(*ySeq != '\0') {  
          putchar(*ySeq);  
          ++ySeq;  
          if (++ycount%60 == 0) putchar ('\n');  
        }  
        putchar ('\n');  
          
        /* reset sequence holders */  
        strcpy (xSeq, "");  
        strcpy (ySeq, "");  
  
 inBlock = FALSE ;  
      }  
  }  
}  
  
  
void usage (FILE * ofp)  
{  
  fprintf (ofp, "dba version: %s\n",VERSION_NUMBER) ;  
  fprintf (ofp, "Usage: dba [options] seq1 seq2\n") ;  
  fprintf (ofp, " -matchA [%g]     match level A\n",MATCH65) ;  
  fprintf (ofp, " -matchB [%g]     match level B\n",MATCH75) ;  
  fprintf (ofp, " -matchC [%g]     match level C\n",MATCH85) ;  
  fprintf (ofp, " -matchD [%g]     match level D\n",MATCH95) ;  
  fprintf (ofp, " -gap [%g]       gap probability\n",GAP) ;  
  fprintf (ofp, " -blockopen [%g] block open probability\n",BLOCKOPEN) ;  
  fprintf (ofp, " -umatch [%g]    unmatched gap probability\n",UNMATCHED_PEN) ;  
  fprintf (ofp, " -single           use only one match level, set with -matchA [%g]\n", MATCHS) ;  
  fprintf (ofp, " -nomatchn         do not match N to any base\n") ;  
  fprintf (ofp, " -align           show alignment for computer parsing\n") ;  
  fprintf (ofp, " -pretty           show alignment for ASCII viewing\n") ;  
  fprintf (ofp, " -pff           show phylogenetic footprinting format output\n") ;  
  fprintf (ofp, "             (gapped fasta)\n") ;  
  fprintf (ofp, " -label           show label alignment\n") ;  
  fprintf (ofp, " -params           print parameters\n") ;  
  fprintf (ofp, "       -kbyte [%d]    Kilobytes allowed for main matrix memory\n",kbyte);  
  show_help_DPRunImpl(ofp);  
  show_standard_options(ofp);  
    
  exit (-1) ;  
}  
  
void show_version(FILE * ofp)  
{  
  fprintf(ofp,"dba version %s\n",VERSION_NUMBER);  
  fprintf(ofp,"  Released %s\n",RELEASE_DAY);  
  fprintf(ofp,"  Compiled %s\n",COMPILE_DATE);  
  fprintf(ofp,"dba was written by Niclas Jareborg, Ewan Birney and Richard Durbin\n");  
  fprintf(ofp,"Copyright (c) 1998,1999,2000,2001 GRL ltd. It is distributed under a Gnu Public License\n");  
  fprintf(ofp,"See GNULICENSE in source directory for more information\n");  
}  
  
  
int main (int argc, char **argv)  
{  
  Sequence *one, *two ;  
  ComplexSequence *cone, *ctwo ;  
  ComplexSequenceEvalSet *cses ;  
  DnaMatrix *dcm65 ;  
  DnaMatrix *dcm75 ;  
  DnaMatrix *dcm85 ;  
  DnaMatrix *dcm95 ;  
  DnaProbMatrix * dpm65;  
  DnaProbMatrix * dpm75;  
  DnaProbMatrix * dpm85;  
  DnaProbMatrix * dpm95;  
  int temp;  
  PackAln *pal ;  
  AlnBlock *alb ;  
  int isShowAlign = 0 ;  
  int isShowParams = 0 ;  
  int show_label_align = 0;  
  int show_pretty_align = 0;  
  int isNoMatchN = 0; /* changed Jan 28 (nic) */  
  int isSingle = 0;  
  int aSet = 0;  
  int show_pff = 0;  
  
  DnaMatchBlock * dmat;  
  DPRunImpl * dpri;  
  
  int score;  
    
  strip_out_standard_options(&argc,argv,usage,show_version);  
  
  /* run time */  
    
  dpri = new_DPRunImpl_from_argv(&argc,argv);  
  
  for (++argv, --argc ; argc > 2 ; ++argv, --argc)  
    if (!strcmp (*argv, "-align"))  
      isShowAlign = 1 ;  
    else if (!strcmp (*argv, "-params"))  
      isShowParams = 1 ;  
    else if (!strcmp (*argv, "-label"))  
      show_label_align = 1;  
    else if (!strcmp (*argv, "-pretty"))  
      show_pretty_align = 1;  
    else if (!strcmp (*argv, "-matchA"))  
       { MATCH65 = atof (*++argv) ; --argc ; aSet = 1; }     
    else if (!strcmp (*argv, "-matchB"))  
       { MATCH75 = atof (*++argv) ; --argc ; }     
    else if (!strcmp (*argv, "-matchC"))  
       { MATCH85 = atof (*++argv) ; --argc ; }     
    else if (!strcmp (*argv, "-matchD"))  
       { MATCH95 = atof (*++argv) ; --argc ; }     
    else if (!strcmp (*argv, "-umatch"))  
      { UNMATCHED_PEN = atof (*++argv) ; --argc ; }  
    else if (!strcmp (*argv, "-gap"))  
      { GAP = atof (*++argv) ; --argc ; }  
    else if (!strcmp (*argv, "-kbyte"))  
      { kbyte = atoi (*++argv) ; --argc ; }  
    else if (!strcmp (*argv, "-blockopen"))  
      { BLOCKOPEN = atof (*++argv) ; --argc ; }  
    else if (!strcmp (*argv, "-single"))  
       { isSingle = 1 ; }     
    else if (!strcmp (*argv, "-nomatchn")) /* changed Jan 28 (nic) */  
      isNoMatchN = 1 ;   
    else if (!strcmp (*argv, "-pff"))   
      show_pff = 1 ;   
    else  
      { fprintf (stderr, "option %s not recognized\n", *argv) ;  
 usage(stdout) ;  
      }  
  if (argc < 2)  
    usage(stdout) ;  
    
  if(isSingle) {  
    MATCH75 = 0;  
    MATCH85 = 0;  
    MATCH95 = 0;  
    if (!aSet) {  
      MATCH65 = MATCHS;  
    }  
  }  
  
  if (isShowParams)  
    printf ("MATCH_A = %g, MATCH_B = %g, MATCH_C = %g, MATCH_D = %g, UMATCH = %g, GAP = %g, BLOCKOPEN = %g\n",  
     MATCH65, MATCH75, MATCH85, MATCH95, UNMATCHED_PEN, GAP, BLOCKOPEN) ;  
  
  /*   
   * we should use type-safe DNA *   
   * here, but what the hell...  
   *  
   */  
  
  one = read_fasta_file_Sequence(*argv++);  
  two = read_fasta_file_Sequence(*argv++);  
  
  uppercase_Sequence(one);  
  uppercase_Sequence(two);  
  /*  
   * make sure sequences are DNA. Complain if they over 50% no ATGCN  
   */  
  
  force_to_dna_Sequence(one,0.0,&temp);  
  if( temp > (one->len/2) ) {  
    warn("Sequence %s has more than 50% of its residues no ATGCN. Are you sure its DNA?",one->name);  
  }  
  
  force_to_dna_Sequence(two,0.0,&temp);  
  if( temp > (two->len/2) ) {  
    warn("Sequence %s has more than 50% of its residues no ATGCN. Are you sure its DNA?",two->name);  
  }  
  
 one->type = SEQUENCE_DNA;  
 two->type = SEQUENCE_DNA;  
  
  /** make them into DNA complex sequences **/  
  
  cses = default_DNA_ComplexSequenceEvalSet();  
  cone = new_ComplexSequence(one,cses);  
  ctwo = new_ComplexSequence(two,cses);  
  
  if( cone == NULL || ctwo == NULL ) {  
    fatal("For some reason, unable to make complexsequences");  
  }  
  
  
  /** make DnaProbMat's **/  
  
  dpm65 = DnaProbMatrix_from_match(MATCH65,isNoMatchN == 1 ? NMaskType_BANNED : NMaskType_VARIABLE);  
  if( dpm65 == NULL ) {  
    fatal("Could not build DnaProbMatrix for MATCH65");  
  }  
  flat_null_DnaProbMatrix(dpm65);  
  
  dpm75 = DnaProbMatrix_from_match(MATCH75,isNoMatchN == 1 ? NMaskType_BANNED : NMaskType_VARIABLE);  
  if( dpm75 == NULL ) {  
    fatal("Could not build DnaProbMatrix for MATCH75");  
  }  
  flat_null_DnaProbMatrix(dpm75);  
  
  dpm85 = DnaProbMatrix_from_match(MATCH85,isNoMatchN == 1 ? NMaskType_BANNED : NMaskType_VARIABLE);  
  if( dpm85 == NULL ) {  
    fatal("Could not build DnaProbMatrix for MATCH85");  
  }  
  flat_null_DnaProbMatrix(dpm85);  
  
  dpm95 = DnaProbMatrix_from_match(MATCH95,isNoMatchN == 1 ? NMaskType_BANNED : NMaskType_VARIABLE);  
  if( dpm95 == NULL ) {  
    fatal("Could not build DnaProbMatrix for MATCH95");  
  }  
  flat_null_DnaProbMatrix(dpm95);  
  
  
  dcm65 = DnaMatrix_from_DnaProbMatrix(dpm65);  
  if( dcm65 == NULL ) {  
    fatal("Could not build DnaCompMat 65");  
  }  
  dcm75 = DnaMatrix_from_DnaProbMatrix(dpm75);  
  if( dcm75 == NULL ) {  
    fatal("Could not build DnaCompMat 75");  
  }  
  dcm85 = DnaMatrix_from_DnaProbMatrix(dpm85);  
  if( dcm85 == NULL ) {  
    fatal("Could not build DnaCompMat 85");  
  }  
  dcm95 = DnaMatrix_from_DnaProbMatrix(dpm95);  
  if( dcm95 == NULL ) {  
    fatal("Could not build DnaCompMat 95");  
  }  
  
  
  /** make the alignment **/  
  
  change_max_BaseMatrix_kbytes (kbyte) ; /* 25Mb */  
  pal = PackAln_bestmemory_DnaMatchBlock (cone, ctwo,   
       dcm65,  
       dcm75,  
       dcm85,  
       dcm95,  
    Probability2Score(GAP/UNMATCHED_PEN),   
    Probability2Score(1),  
       /* horrible hack here so that   
          single searches are not / 4 */  
  
       Probability2Score(((1-UNMATCHED_PEN)/(UNMATCHED_PEN))/ (isSingle == TRUE ? 1 : 4)),   
    Probability2Score((1.0 - (GAP + GAP + BLOCKOPEN))/(UNMATCHED_PEN*UNMATCHED_PEN)),  
    Probability2Score(BLOCKOPEN/UNMATCHED_PEN),  
       NULL,dpri    
      );  
  
#ifdef DOPROB  
  dmat = matrix_logsum_DnaMatchBlock(cone, ctwo,   
       dcm65,  
       dcm75,  
       dcm85,  
       dcm95,  
    Probability2Score(GAP/UNMATCHED_PEN),   
    Probability2Score(1),  
       /* horrible hack here so that   
          single searches are not / 4 */  
  
       Probability2Score(((1-UNMATCHED_PEN)/(UNMATCHED_PEN))/ (isSingle == TRUE ? 1 : 4)),   
    Probability2Score((1.0 - (GAP + GAP + BLOCKOPEN))/(UNMATCHED_PEN*UNMATCHED_PEN)),  
    Probability2Score(BLOCKOPEN/UNMATCHED_PEN)   );  
  
  fprintf(stdout,"And the score is %d %.2f\n",dmat->basematrix->matrix[dmat->basematrix->leni-1][dmat->basematrix->lenj-1],Score2Bits(dmat->basematrix->matrix[dmat->basematrix->leni-1][dmat->basematrix->lenj-1]));  
  
  score = score_only_logsum_DnaMatchBlock(cone, ctwo,   
       dcm65,  
       dcm75,  
       dcm85,  
       dcm95,  
    Probability2Score(GAP/UNMATCHED_PEN),   
    Probability2Score(1),  
       /* horrible hack here so that   
          single searches are not / 4 */  
  
       Probability2Score(((1-UNMATCHED_PEN)/(UNMATCHED_PEN))/ (isSingle == TRUE ? 1 : 4)),   
    Probability2Score((1.0 - (GAP + GAP + BLOCKOPEN))/(UNMATCHED_PEN*UNMATCHED_PEN)),  
    Probability2Score(BLOCKOPEN/UNMATCHED_PEN)   );  
  
  
  fprintf(stdout,"And the score is %d %.2f\n",score,Score2Bits(score));  
  
#endif  
  
  if( pal == NULL )   
    fatal("Unable to build alignment!");  
  
  
  alb = convert_PackAln_to_AlnBlock_DnaMatchBlock(pal);  
  
  if( show_label_align ) {  
    printf ("score = %.2f\n", Score2Bits(pal->score)) ;  
    mapped_ascii_AlnBlock(alb,Score2Probability,0,stdout);  
  } else if ( show_pretty_align ) {  
    printf ("score = %.2f\n", Score2Bits(pal->score)) ;  
    show_pretty_dba_align(alb,one,two,stdout);  
  } else if (isShowAlign) {  
    printf ("score = %.2f\n", Score2Bits(pal->score)) ;  
    prettyAlnBlock (alb, one, two) ;  
  } else if (show_pff) {  
    showPFF (alb, one, two) ;  
  } else {  
    printf ("score = %.2f\n", Score2Bits(pal->score)) ;  
    alnBlockSummary (alb, one, two) ;  
  }  
  return 0;  
}  
  
  
/******************* end of file ***************/  
