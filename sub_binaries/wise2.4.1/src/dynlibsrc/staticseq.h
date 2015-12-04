#ifndef DYNAMITEstaticseqHEADERFILE
#define DYNAMITEstaticseqHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "sequence.h"
#include <glib.h>

/* Object StaticSeqHolder
 *
 * Descrip: Makes sequences with optimised
 *        static memory, but the sequences
 *        cannot be freed individually
 *
 *
 */
struct Wise2_StaticSeqHolder {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    GStringChunk * gstring_chunk;    
    } ;  
/* StaticSeqHolder defined */ 
#ifndef DYNAMITE_DEFINED_StaticSeqHolder
typedef struct Wise2_StaticSeqHolder Wise2_StaticSeqHolder;
#define StaticSeqHolder Wise2_StaticSeqHolder
#define DYNAMITE_DEFINED_StaticSeqHolder
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  new_StaticSeqHolder(void)
 *
 * Descrip:    makes a new StaticSeqHolder
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [StaticSeqHolder *]
 *
 */
StaticSeqHolder * Wise2_new_StaticSeqHolder(void);
#define new_StaticSeqHolder Wise2_new_StaticSeqHolder


/* Function:  free_GStringChunk(gs)
 *
 * Descrip:    for registering glib thingy for freeing
 *
 *
 * Arg:        gs [UNKN ] Undocumented argument [GStringChunk *]
 *
 * Return [UNKN ]  Undocumented return value [GStringChunk *]
 *
 */
GStringChunk * Wise2_free_GStringChunk(GStringChunk * gs);
#define free_GStringChunk Wise2_free_GStringChunk


/* Function:  new_Sequence_StaticSeqHolder(ssh,seq)
 *
 * Descrip:    Making a new sequence from a staticseq holder - consumes
 *             the sequence object (actually recycling the shell of it)
 *
 *
 * Arg:        ssh [UNKN ] Undocumented argument [StaticSeqHolder *]
 * Arg:        seq [UNKN ] Undocumented argument [Sequence *]
 *
 * Return [UNKN ]  Undocumented return value [Sequence *]
 *
 */
Sequence * Wise2_new_Sequence_StaticSeqHolder(StaticSeqHolder * ssh,Sequence * seq);
#define new_Sequence_StaticSeqHolder Wise2_new_Sequence_StaticSeqHolder


/* Function:  hard_link_StaticSeqHolder(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [StaticSeqHolder *]
 *
 * Return [UNKN ]  Undocumented return value [StaticSeqHolder *]
 *
 */
StaticSeqHolder * Wise2_hard_link_StaticSeqHolder(StaticSeqHolder * obj);
#define hard_link_StaticSeqHolder Wise2_hard_link_StaticSeqHolder


/* Function:  StaticSeqHolder_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [StaticSeqHolder *]
 *
 */
StaticSeqHolder * Wise2_StaticSeqHolder_alloc(void);
#define StaticSeqHolder_alloc Wise2_StaticSeqHolder_alloc


/* Function:  free_StaticSeqHolder(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [StaticSeqHolder *]
 *
 * Return [UNKN ]  Undocumented return value [StaticSeqHolder *]
 *
 */
StaticSeqHolder * Wise2_free_StaticSeqHolder(StaticSeqHolder * obj);
#define free_StaticSeqHolder Wise2_free_StaticSeqHolder


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/

#ifdef _cplusplus
}
#endif

#endif
