#ifndef DYNAMITEshadowseqHEADERFILE
#define DYNAMITEshadowseqHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "sequence.h"

#define SHADOW_SEQUENCE_REGION_BASIC 4

#define ShadowSequenceLISTLENGTH 8

struct Wise2_ShadowSeqRegion {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int start_shadow;    
    int start_seq;   
    int len;     
    Sequence * seq;  
    } ;  
/* ShadowSeqRegion defined */ 
#ifndef DYNAMITE_DEFINED_ShadowSeqRegion
typedef struct Wise2_ShadowSeqRegion Wise2_ShadowSeqRegion;
#define ShadowSeqRegion Wise2_ShadowSeqRegion
#define DYNAMITE_DEFINED_ShadowSeqRegion
#endif


struct Wise2_ShadowSequence {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Sequence * seq;  
    char dirty;  
    ShadowSeqRegion ** region;   
    int len;/* len for above region  */ 
    int maxlen; /* maxlen for above region */ 
    } ;  
/* ShadowSequence defined */ 
#ifndef DYNAMITE_DEFINED_ShadowSequence
typedef struct Wise2_ShadowSequence Wise2_ShadowSequence;
#define ShadowSequence Wise2_ShadowSequence
#define DYNAMITE_DEFINED_ShadowSequence
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  dump_ShadowSequence(shadow,ofp)
 *
 * Descrip:    Dumps shadow sequence out to file
 *
 *
 * Arg:        shadow [UNKN ] Undocumented argument [ShadowSequence *]
 * Arg:           ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_dump_ShadowSequence(ShadowSequence * shadow,FILE * ofp);
#define dump_ShadowSequence Wise2_dump_ShadowSequence


/* Function:  add_if_possible_ShadowSequence(shadow,seq,min_extension,start_shadow,start_seq,shadow_rate)
 *
 * Descrip:    Looks as to whether there is a good extension to be made, if
 *             so, does it and adds it to the ShadowSeq. Returns 0 if extension
 *             fails, if succeeds returns end point
 *
 *
 * Arg:               shadow [UNKN ] Undocumented argument [ShadowSequence *]
 * Arg:                  seq [UNKN ] Undocumented argument [Sequence *]
 * Arg:        min_extension [UNKN ] Undocumented argument [int]
 * Arg:         start_shadow [UNKN ] Undocumented argument [int]
 * Arg:            start_seq [UNKN ] Undocumented argument [int]
 * Arg:          shadow_rate [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_add_if_possible_ShadowSequence(ShadowSequence * shadow,Sequence * seq,int min_extension,int start_shadow,int start_seq,int shadow_rate);
#define add_if_possible_ShadowSequence Wise2_add_if_possible_ShadowSequence


/* Function:  add_region_ShadowSequence(shadow,seq,start_shadow,start_seq,len)
 *
 * Descrip:    Adds a region to a ShadowSequence, extending
 *             the array if necessary
 *
 *
 * Arg:              shadow [UNKN ] Undocumented argument [ShadowSequence *]
 * Arg:                 seq [UNKN ] Undocumented argument [Sequence *]
 * Arg:        start_shadow [UNKN ] Undocumented argument [int]
 * Arg:           start_seq [UNKN ] Undocumented argument [int]
 * Arg:                 len [UNKN ] Undocumented argument [int]
 *
 */
void Wise2_add_region_ShadowSequence(ShadowSequence * shadow,Sequence * seq,int start_shadow,int start_seq,int len);
#define add_region_ShadowSequence Wise2_add_region_ShadowSequence


/* Function:  new_ShadowSequence(seq)
 *
 * Descrip:    Builds a new ShadowSequence from a Sequence
 *             doesnt hard link as memory should be handled
 *             outside of the shadowsystem
 *
 *
 * Arg:        seq [UNKN ] Undocumented argument [Sequence *]
 *
 * Return [UNKN ]  Undocumented return value [ShadowSequence *]
 *
 */
ShadowSequence * Wise2_new_ShadowSequence(Sequence * seq);
#define new_ShadowSequence Wise2_new_ShadowSequence


/* Function:  hard_link_ShadowSeqRegion(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ShadowSeqRegion *]
 *
 * Return [UNKN ]  Undocumented return value [ShadowSeqRegion *]
 *
 */
ShadowSeqRegion * Wise2_hard_link_ShadowSeqRegion(ShadowSeqRegion * obj);
#define hard_link_ShadowSeqRegion Wise2_hard_link_ShadowSeqRegion


/* Function:  ShadowSeqRegion_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ShadowSeqRegion *]
 *
 */
ShadowSeqRegion * Wise2_ShadowSeqRegion_alloc(void);
#define ShadowSeqRegion_alloc Wise2_ShadowSeqRegion_alloc


/* Function:  free_ShadowSeqRegion(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ShadowSeqRegion *]
 *
 * Return [UNKN ]  Undocumented return value [ShadowSeqRegion *]
 *
 */
ShadowSeqRegion * Wise2_free_ShadowSeqRegion(ShadowSeqRegion * obj);
#define free_ShadowSeqRegion Wise2_free_ShadowSeqRegion


/* Function:  add_ShadowSequence(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [ShadowSequence *]
 * Arg:        add [OWNER] Object to add to the list [ShadowSeqRegion *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_ShadowSequence(ShadowSequence * obj,ShadowSeqRegion * add);
#define add_ShadowSequence Wise2_add_ShadowSequence


/* Function:  flush_ShadowSequence(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [ShadowSequence *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_ShadowSequence(ShadowSequence * obj);
#define flush_ShadowSequence Wise2_flush_ShadowSequence


/* Function:  ShadowSequence_alloc_std(void)
 *
 * Descrip:    Equivalent to ShadowSequence_alloc_len(ShadowSequenceLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ShadowSequence *]
 *
 */
ShadowSequence * Wise2_ShadowSequence_alloc_std(void);
#define ShadowSequence_alloc_std Wise2_ShadowSequence_alloc_std


/* Function:  ShadowSequence_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [ShadowSequence *]
 *
 */
ShadowSequence * Wise2_ShadowSequence_alloc_len(int len);
#define ShadowSequence_alloc_len Wise2_ShadowSequence_alloc_len


/* Function:  hard_link_ShadowSequence(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ShadowSequence *]
 *
 * Return [UNKN ]  Undocumented return value [ShadowSequence *]
 *
 */
ShadowSequence * Wise2_hard_link_ShadowSequence(ShadowSequence * obj);
#define hard_link_ShadowSequence Wise2_hard_link_ShadowSequence


/* Function:  ShadowSequence_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ShadowSequence *]
 *
 */
ShadowSequence * Wise2_ShadowSequence_alloc(void);
#define ShadowSequence_alloc Wise2_ShadowSequence_alloc


/* Function:  free_ShadowSequence(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ShadowSequence *]
 *
 * Return [UNKN ]  Undocumented return value [ShadowSequence *]
 *
 */
ShadowSequence * Wise2_free_ShadowSequence(ShadowSequence * obj);
#define free_ShadowSequence Wise2_free_ShadowSequence


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_ShadowSequence(ShadowSeqRegion ** list,int i,int j) ;
#define swap_ShadowSequence Wise2_swap_ShadowSequence
void Wise2_qsort_ShadowSequence(ShadowSeqRegion ** list,int left,int right,int (*comp)(ShadowSeqRegion * ,ShadowSeqRegion * ));
#define qsort_ShadowSequence Wise2_qsort_ShadowSequence
void Wise2_sort_ShadowSequence(ShadowSequence * obj,int (*comp)(ShadowSeqRegion *, ShadowSeqRegion *));
#define sort_ShadowSequence Wise2_sort_ShadowSequence
boolean Wise2_expand_ShadowSequence(ShadowSequence * obj,int len);
#define expand_ShadowSequence Wise2_expand_ShadowSequence

#ifdef _cplusplus
}
#endif

#endif
