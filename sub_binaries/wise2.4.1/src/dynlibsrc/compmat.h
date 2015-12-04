#ifndef DYNAMITEcompmatHEADERFILE
#define DYNAMITEcompmatHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "wisebase.h"
#include "probability.h"
#include "randommodel.h"


#define CompMat_AAMATCH(comp_mat,aa1,aa2) (comp_mat->comp[aa1][aa2])
/* Object CompProb
 *
 * Descrip: The probabilistic form of CompMat
 *
 *
 */
struct Wise2_CompProb {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Probability comp[26][26];    
    char * name;     
    } ;  
/* CompProb defined */ 
#ifndef DYNAMITE_DEFINED_CompProb
typedef struct Wise2_CompProb Wise2_CompProb;
#define CompProb Wise2_CompProb
#define DYNAMITE_DEFINED_CompProb
#endif


/* Object CompMat
 *
 * Descrip: This object stores BLOSUM and PAM 
 *        comparison matrices. It stores them as
 *        scores: NB - this means probabilistically
 *        we are talking about some arbitary base of
 *        log which is really annoying.
 *
 *
 */
struct Wise2_CompMat {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Score comp[26][26];  
    char * name;    /*  if any, could be NULL */ 
    } ;  
/* CompMat defined */ 
#ifndef DYNAMITE_DEFINED_CompMat
typedef struct Wise2_CompMat Wise2_CompMat;
#define CompMat Wise2_CompMat
#define DYNAMITE_DEFINED_CompMat
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  simple_CompProb(match,rnd)
 *
 * Descrip:    Makes a simple CompProb matrix 
 *
 *
 * Arg:        match [UNKN ] Undocumented argument [Probability]
 * Arg:          rnd [UNKN ] Undocumented argument [Probability]
 *
 * Return [UNKN ]  Undocumented return value [CompProb *]
 *
 */
CompProb * Wise2_simple_CompProb(Probability match,Probability rnd);
#define simple_CompProb Wise2_simple_CompProb


/* Function:  fold_column_RandomModel_CompProb(cp,rm)
 *
 * Descrip:    Folds a random model in over the columns
 *
 *
 * Arg:        cp [UNKN ] Undocumented argument [CompProb *]
 * Arg:        rm [UNKN ] Undocumented argument [RandomModel *]
 *
 */
void Wise2_fold_column_RandomModel_CompProb(CompProb * cp,RandomModel * rm);
#define fold_column_RandomModel_CompProb Wise2_fold_column_RandomModel_CompProb


/* Function:  simple_aa_CompProb(match,set,rnd)
 *
 * Descrip:    Makes a simple CompProb with simple aa rules
 *
 *
 * Arg:        match [UNKN ] Undocumented argument [Probability]
 * Arg:          set [UNKN ] Undocumented argument [Probability]
 * Arg:          rnd [UNKN ] Undocumented argument [Probability]
 *
 * Return [UNKN ]  Undocumented return value [CompProb *]
 *
 */
CompProb * Wise2_simple_aa_CompProb(Probability match,Probability set,Probability rnd);
#define simple_aa_CompProb Wise2_simple_aa_CompProb


/* Function:  CompMat_from_CompProb(cp)
 *
 * Descrip:    Maps a CompProb to a CompMat going through
 *             Probability2Score
 *
 *
 * Arg:        cp [UNKN ] Undocumented argument [CompProb *]
 *
 * Return [UNKN ]  Undocumented return value [CompMat *]
 *
 */
CompMat * Wise2_CompMat_from_CompProb(CompProb * cp);
#define CompMat_from_CompProb Wise2_CompMat_from_CompProb


/* Function:  CompProb_from_halfbit(cm)
 *
 * Descrip:    Maps a halfbit matrix to a prob matrix by rebasing
 *             etc. 
 *
 *             *Really* not sensible!
 *
 *
 * Arg:        cm [UNKN ] Undocumented argument [CompMat *]
 *
 * Return [UNKN ]  Undocumented return value [CompProb *]
 *
 */
CompProb * Wise2_CompProb_from_halfbit(CompMat * cm);
#define CompProb_from_halfbit Wise2_CompProb_from_halfbit


/* Function:  CompMat_from_halfbit(cm)
 *
 * Descrip:    flips a halfbit based matrix (eg, blosum62) into a score
 *             based matrix just by rebasing the log etc. 
 *
 *             Not a sensible function ...
 *
 *
 *
 * Arg:        cm [UNKN ] Undocumented argument [CompMat *]
 *
 * Return [UNKN ]  Undocumented return value [CompMat *]
 *
 */
CompMat * Wise2_CompMat_from_halfbit(CompMat * cm);
#define CompMat_from_halfbit Wise2_CompMat_from_halfbit


/* Function:  factor_CompMat(cm,factor)
 *
 * Descrip:    multiples all the scores by the amount
 *
 *
 * Arg:            cm [UNKN ] compmat object [CompMat *]
 * Arg:        factor [UNKN ] amount to multiple by [int]
 *
 */
void Wise2_factor_CompMat(CompMat * cm,int factor);
#define factor_CompMat Wise2_factor_CompMat


/* Function:  fail_safe_CompMat_access(cm,aa1,aa2)
 *
 * Descrip:    gives the fail form of the macro CompMat_AAMATCH which 
 *             checks that aa1 and a2 are sensible and that cm is not NULL.
 *
 *
 * Arg:         cm [UNKN ] compmat object [CompMat *]
 * Arg:        aa1 [UNKN ] first amino acid [int]
 * Arg:        aa2 [UNKN ] second amino acid [int]
 *
 * Return [UNKN ]  Undocumented return value [Score]
 *
 */
Score Wise2_fail_safe_CompMat_access(CompMat * cm,int aa1,int aa2);
#define fail_safe_CompMat_access Wise2_fail_safe_CompMat_access


/* Function:  write_Blast_CompMat(cm,ofp)
 *
 * Descrip:    writes a protien CompMat with a standard
 *             alphabet.
 *
 *
 * Arg:         cm [UNKN ] CompMat object [CompMat *]
 * Arg:        ofp [UNKN ] file to output [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_write_Blast_CompMat(CompMat * cm,FILE * ofp);
#define write_Blast_CompMat Wise2_write_Blast_CompMat


/* Function:  write_Blast_CompMat_alphabet(cm,alphabet,ofp)
 *
 * Descrip:    actualy writes out the Blast CFormat. The alphabet is 
 *             what order you want the amino acids. If you want the 
 *             standard format use /write_Blast_CompMat
 *
 *
 * Arg:              cm [UNKN ] comp mat object [CompMat *]
 * Arg:        alphabet [UNKN ] string for alphabet to be used [char *]
 * Arg:             ofp [UNKN ] fileoutput [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_write_Blast_CompMat_alphabet(CompMat * cm,char * alphabet,FILE * ofp);
#define write_Blast_CompMat_alphabet Wise2_write_Blast_CompMat_alphabet


/* Function:  read_Blast_file_CompMat(filename)
 *
 * Descrip:    Opens file, reads matrix, closes file.
 *             calls /read_Blast_CompMat for the actual format
 *             reading. Uses /openfile to open the file,
 *             so will open from config files.
 *
 *
 * Arg:        filename [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [CompMat *]
 *
 */
CompMat * Wise2_read_Blast_file_CompMat(char * filename);
#define read_Blast_file_CompMat Wise2_read_Blast_file_CompMat


/* Function:  read_Blast_CompMat(ifp)
 *
 * Descrip:    reads a BLAST format matrix and
 *             allocates a new ComMat structure.
 *
 *
 * Arg:        ifp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [CompMat *]
 *
 */
CompMat * Wise2_read_Blast_CompMat(FILE * ifp);
#define read_Blast_CompMat Wise2_read_Blast_CompMat


/* Function:  read_Blast_file_CompProb(filename)
 *
 * Descrip:    Reads a BLAST format comp prob from file
 *
 *
 * Arg:        filename [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [CompProb *]
 *
 */
CompProb * Wise2_read_Blast_file_CompProb(char * filename);
#define read_Blast_file_CompProb Wise2_read_Blast_file_CompProb


/* Function:  read_Blast_CompProb(ifp)
 *
 * Descrip:    reads a BLAST format matrix and
 *             allocates a new CompProb structure.
 *
 *
 * Arg:        ifp [READ ] file input [FILE *]
 *
 * Return [UNKN ]  newly allocated CompProb [CompProb *]
 *
 */
CompProb * Wise2_read_Blast_CompProb(FILE * ifp);
#define read_Blast_CompProb Wise2_read_Blast_CompProb


/* Function:  blank_CompMat(void)
 *
 * Descrip:    makes a 0,0 matrix
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [CompMat *]
 *
 */
CompMat * Wise2_blank_CompMat(void);
#define blank_CompMat Wise2_blank_CompMat


/* Function:  blank_CompProb(void)
 *
 * Descrip:    makes a 1.0 prob matrix
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [CompProb *]
 *
 */
CompProb * Wise2_blank_CompProb(void);
#define blank_CompProb Wise2_blank_CompProb


/* Function:  hard_link_CompProb(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [CompProb *]
 *
 * Return [UNKN ]  Undocumented return value [CompProb *]
 *
 */
CompProb * Wise2_hard_link_CompProb(CompProb * obj);
#define hard_link_CompProb Wise2_hard_link_CompProb


/* Function:  CompProb_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [CompProb *]
 *
 */
CompProb * Wise2_CompProb_alloc(void);
#define CompProb_alloc Wise2_CompProb_alloc


/* Function:  free_CompProb(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [CompProb *]
 *
 * Return [UNKN ]  Undocumented return value [CompProb *]
 *
 */
CompProb * Wise2_free_CompProb(CompProb * obj);
#define free_CompProb Wise2_free_CompProb


/* Function:  hard_link_CompMat(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [CompMat *]
 *
 * Return [UNKN ]  Undocumented return value [CompMat *]
 *
 */
CompMat * Wise2_hard_link_CompMat(CompMat * obj);
#define hard_link_CompMat Wise2_hard_link_CompMat


/* Function:  CompMat_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [CompMat *]
 *
 */
CompMat * Wise2_CompMat_alloc(void);
#define CompMat_alloc Wise2_CompMat_alloc


/* Function:  free_CompMat(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [CompMat *]
 *
 * Return [UNKN ]  Undocumented return value [CompMat *]
 *
 */
CompMat * Wise2_free_CompMat(CompMat * obj);
#define free_CompMat Wise2_free_CompMat


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
boolean Wise2_replace_name_CompMat(CompMat * obj,char * name);
#define replace_name_CompMat Wise2_replace_name_CompMat
char * Wise2_access_name_CompMat(CompMat * obj);
#define access_name_CompMat Wise2_access_name_CompMat
boolean Wise2_bad_CompMat_alphabet(char * al);
#define bad_CompMat_alphabet Wise2_bad_CompMat_alphabet

#ifdef _cplusplus
}
#endif

#endif
