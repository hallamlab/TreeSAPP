#ifndef DYNAMITEseqerrorHEADERFILE
#define DYNAMITEseqerrorHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "wisebase.h"
#include "probability.h"
#include "aln.h"
#include "sequence.h"

enum SeqErrorType {
  SeqErrorInsertion,
  SeqErrorDeletion,
  SeqErrorSubstitution
};

#define SequenceErrorSetLISTLENGTH 128

/* Object SequenceError
 *
 * Descrip: This holds information about
 *        a particular insertion or
 *        deletion for one sequence
 *
 *
 */
struct Wise2_SequenceError {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int type;    
    int start;   
    int end;     
    char * replaced_bases;   
    Probability probability;     
    char * inserted_bases;   
    char * suspected_deletion;   
    } ;  
/* SequenceError defined */ 
#ifndef DYNAMITE_DEFINED_SequenceError
typedef struct Wise2_SequenceError Wise2_SequenceError;
#define SequenceError Wise2_SequenceError
#define DYNAMITE_DEFINED_SequenceError
#endif


/* Object SequenceErrorSet
 *
 * Descrip: This holds a set of insertion
 *        or deletions for a sequence
 *
 *
 */
struct Wise2_SequenceErrorSet {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    SequenceError ** error;  
    int len;/* len for above error  */ 
    int maxlen; /* maxlen for above error */ 
    } ;  
/* SequenceErrorSet defined */ 
#ifndef DYNAMITE_DEFINED_SequenceErrorSet
typedef struct Wise2_SequenceErrorSet Wise2_SequenceErrorSet;
#define SequenceErrorSet Wise2_SequenceErrorSet
#define DYNAMITE_DEFINED_SequenceErrorSet
#endif


/* Object ErrorSequence
 *
 * Descrip: This holds a sequence and what
 *        errors have occured in it.
 *
 *
 */
struct Wise2_ErrorSequence {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Sequence * seq;  
    SequenceErrorSet * ses;  
    } ;  
/* ErrorSequence defined */ 
#ifndef DYNAMITE_DEFINED_ErrorSequence
typedef struct Wise2_ErrorSequence Wise2_ErrorSequence;
#define ErrorSequence Wise2_ErrorSequence
#define DYNAMITE_DEFINED_ErrorSequence
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  genewise_SequenceErrorSet(als)
 *
 * Descrip:    Makes a sequence error set from standard genewise labels
 *
 *
 * Arg:        als [UNKN ] Undocumented argument [AlnSequence *]
 *
 * Return [UNKN ]  Undocumented return value [SequenceErrorSet *]
 *
 */
SequenceErrorSet * Wise2_genewise_SequenceErrorSet(AlnSequence * als );
#define genewise_SequenceErrorSet Wise2_genewise_SequenceErrorSet


/* Function:  show_SequenceErrorSet(ses,ofp)
 *
 * Descrip:    Displays a set of sequence errors in space deliminted format
 *
 *
 * Arg:        ses [UNKN ] Undocumented argument [SequenceErrorSet *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_show_SequenceErrorSet(SequenceErrorSet * ses,FILE * ofp);
#define show_SequenceErrorSet Wise2_show_SequenceErrorSet


/* Function:  show_SequenceError(se,ofp)
 *
 * Descrip:    Displays sequence error in space deliminted format
 *
 *
 * Arg:         se [UNKN ] Undocumented argument [SequenceError *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_show_SequenceError(SequenceError * se,FILE * ofp);
#define show_SequenceError Wise2_show_SequenceError


/* Function:  make_ErrorSequence(seq,subs,insertion,deletion)
 *
 * Descrip:    Makes an error sequence (DNA) with set substitution
 *             and insertion/deletion rates.
 *
 *
 * Arg:              seq [UNKN ] Undocumented argument [Sequence *]
 * Arg:             subs [UNKN ] Undocumented argument [Probability]
 * Arg:        insertion [UNKN ] Undocumented argument [Probability]
 * Arg:         deletion [UNKN ] Undocumented argument [Probability]
 *
 * Return [UNKN ]  Undocumented return value [ErrorSequence *]
 *
 */
ErrorSequence * Wise2_make_ErrorSequence(Sequence * seq,Probability subs,Probability insertion,Probability deletion);
#define make_ErrorSequence Wise2_make_ErrorSequence


/* Function:  hard_link_SequenceError(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SequenceError *]
 *
 * Return [UNKN ]  Undocumented return value [SequenceError *]
 *
 */
SequenceError * Wise2_hard_link_SequenceError(SequenceError * obj);
#define hard_link_SequenceError Wise2_hard_link_SequenceError


/* Function:  SequenceError_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SequenceError *]
 *
 */
SequenceError * Wise2_SequenceError_alloc(void);
#define SequenceError_alloc Wise2_SequenceError_alloc


/* Function:  free_SequenceError(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SequenceError *]
 *
 * Return [UNKN ]  Undocumented return value [SequenceError *]
 *
 */
SequenceError * Wise2_free_SequenceError(SequenceError * obj);
#define free_SequenceError Wise2_free_SequenceError


/* Function:  add_SequenceErrorSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [SequenceErrorSet *]
 * Arg:        add [OWNER] Object to add to the list [SequenceError *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_SequenceErrorSet(SequenceErrorSet * obj,SequenceError * add);
#define add_SequenceErrorSet Wise2_add_SequenceErrorSet


/* Function:  flush_SequenceErrorSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [SequenceErrorSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_SequenceErrorSet(SequenceErrorSet * obj);
#define flush_SequenceErrorSet Wise2_flush_SequenceErrorSet


/* Function:  SequenceErrorSet_alloc_std(void)
 *
 * Descrip:    Equivalent to SequenceErrorSet_alloc_len(SequenceErrorSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SequenceErrorSet *]
 *
 */
SequenceErrorSet * Wise2_SequenceErrorSet_alloc_std(void);
#define SequenceErrorSet_alloc_std Wise2_SequenceErrorSet_alloc_std


/* Function:  SequenceErrorSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [SequenceErrorSet *]
 *
 */
SequenceErrorSet * Wise2_SequenceErrorSet_alloc_len(int len);
#define SequenceErrorSet_alloc_len Wise2_SequenceErrorSet_alloc_len


/* Function:  hard_link_SequenceErrorSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SequenceErrorSet *]
 *
 * Return [UNKN ]  Undocumented return value [SequenceErrorSet *]
 *
 */
SequenceErrorSet * Wise2_hard_link_SequenceErrorSet(SequenceErrorSet * obj);
#define hard_link_SequenceErrorSet Wise2_hard_link_SequenceErrorSet


/* Function:  SequenceErrorSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SequenceErrorSet *]
 *
 */
SequenceErrorSet * Wise2_SequenceErrorSet_alloc(void);
#define SequenceErrorSet_alloc Wise2_SequenceErrorSet_alloc


/* Function:  free_SequenceErrorSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SequenceErrorSet *]
 *
 * Return [UNKN ]  Undocumented return value [SequenceErrorSet *]
 *
 */
SequenceErrorSet * Wise2_free_SequenceErrorSet(SequenceErrorSet * obj);
#define free_SequenceErrorSet Wise2_free_SequenceErrorSet


/* Function:  hard_link_ErrorSequence(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ErrorSequence *]
 *
 * Return [UNKN ]  Undocumented return value [ErrorSequence *]
 *
 */
ErrorSequence * Wise2_hard_link_ErrorSequence(ErrorSequence * obj);
#define hard_link_ErrorSequence Wise2_hard_link_ErrorSequence


/* Function:  ErrorSequence_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ErrorSequence *]
 *
 */
ErrorSequence * Wise2_ErrorSequence_alloc(void);
#define ErrorSequence_alloc Wise2_ErrorSequence_alloc


/* Function:  free_ErrorSequence(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ErrorSequence *]
 *
 * Return [UNKN ]  Undocumented return value [ErrorSequence *]
 *
 */
ErrorSequence * Wise2_free_ErrorSequence(ErrorSequence * obj);
#define free_ErrorSequence Wise2_free_ErrorSequence


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_SequenceErrorSet(SequenceError ** list,int i,int j) ;
#define swap_SequenceErrorSet Wise2_swap_SequenceErrorSet
void Wise2_qsort_SequenceErrorSet(SequenceError ** list,int left,int right,int (*comp)(SequenceError * ,SequenceError * ));
#define qsort_SequenceErrorSet Wise2_qsort_SequenceErrorSet
void Wise2_sort_SequenceErrorSet(SequenceErrorSet * obj,int (*comp)(SequenceError *, SequenceError *));
#define sort_SequenceErrorSet Wise2_sort_SequenceErrorSet
boolean Wise2_expand_SequenceErrorSet(SequenceErrorSet * obj,int len);
#define expand_SequenceErrorSet Wise2_expand_SequenceErrorSet

#ifdef _cplusplus
}
#endif

#endif
