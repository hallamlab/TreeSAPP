#ifndef DYNAMITEpairwiseshortdnaHEADERFILE
#define DYNAMITEpairwiseshortdnaHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"
#include "hsp.h"
#include "subseqhash.h"


struct Wise2_PairwiseShortDna {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    HSPset * forward;    
    HSPset * reverse;    
    } ;  
/* PairwiseShortDna defined */ 
#ifndef DYNAMITE_DEFINED_PairwiseShortDna
typedef struct Wise2_PairwiseShortDna Wise2_PairwiseShortDna;
#define PairwiseShortDna Wise2_PairwiseShortDna
#define DYNAMITE_DEFINED_PairwiseShortDna
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_PairwiseShortDna(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [PairwiseShortDna *]
 *
 * Return [UNKN ]  Undocumented return value [PairwiseShortDna *]
 *
 */
PairwiseShortDna * Wise2_hard_link_PairwiseShortDna(PairwiseShortDna * obj);
#define hard_link_PairwiseShortDna Wise2_hard_link_PairwiseShortDna


/* Function:  PairwiseShortDna_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PairwiseShortDna *]
 *
 */
PairwiseShortDna * Wise2_PairwiseShortDna_alloc(void);
#define PairwiseShortDna_alloc Wise2_PairwiseShortDna_alloc


/* Function:  free_PairwiseShortDna(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [PairwiseShortDna *]
 *
 * Return [UNKN ]  Undocumented return value [PairwiseShortDna *]
 *
 */
PairwiseShortDna * Wise2_free_PairwiseShortDna(PairwiseShortDna * obj);
#define free_PairwiseShortDna Wise2_free_PairwiseShortDna


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
boolean Wise2_process_HSP(HSPset * set,Sequence * query,int query_pos,Sequence * tseq,SeqLookupResultStruct * res_struct,CompMat * mat);
#define process_HSP Wise2_process_HSP
PairwiseShortDna * Wise2_query_to_reverse_target(Sequence * query,Sequence * target,DnaMatrix * dm,int qstart,int qend,int tstart,int tend);
#define query_to_reverse_target Wise2_query_to_reverse_target


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/

#ifdef _cplusplus
}
#endif

#endif
