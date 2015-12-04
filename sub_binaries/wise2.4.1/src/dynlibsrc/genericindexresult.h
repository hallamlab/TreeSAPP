#ifndef DYNAMITEgenericindexresultHEADERFILE
#define DYNAMITEgenericindexresultHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "seqlookup.h"


struct Wise2_GenericIndexResult {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    SeqLookupResultStruct * result;  
    int len;     
    int max_len;     
    int current_pos;     
    } ;  
/* GenericIndexResult defined */ 
#ifndef DYNAMITE_DEFINED_GenericIndexResult
typedef struct Wise2_GenericIndexResult Wise2_GenericIndexResult;
#define GenericIndexResult Wise2_GenericIndexResult
#define DYNAMITE_DEFINED_GenericIndexResult
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  next_interface_GenericIndexResult(data,prev)
 *
 * Descrip:    For interface, returns next position in SeqLookup Results
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 * Arg:        prev [UNKN ] Undocumented argument [SeqLookupResultStruct *]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupResultStruct *]
 *
 */
SeqLookupResultStruct * Wise2_next_interface_GenericIndexResult(void * data,SeqLookupResultStruct * prev) ;
#define next_interface_GenericIndexResult Wise2_next_interface_GenericIndexResult


/* Function:  is_more_interface_GenericIndexResult(data)
 *
 * Descrip:    For interface, indicates whether there is more stuff to find or not
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_is_more_interface_GenericIndexResult(void * data);
#define is_more_interface_GenericIndexResult Wise2_is_more_interface_GenericIndexResult


/* Function:  free_noop_GenericIndexResult(data)
 *
 * Descrip:    Frees the data
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 *
 */
void Wise2_free_noop_GenericIndexResult(void * data);
#define free_noop_GenericIndexResult Wise2_free_noop_GenericIndexResult


/* Function:  add_GenericIndexResult(gir,seq,pos)
 *
 * Descrip:    Adds another result to a IndexResult
 *
 *
 * Arg:        gir [UNKN ] Undocumented argument [GenericIndexResult *]
 * Arg:        seq [UNKN ] Undocumented argument [Sequence *]
 * Arg:        pos [UNKN ] Undocumented argument [int]
 *
 */
void Wise2_add_GenericIndexResult(GenericIndexResult * gir,Sequence * seq,int pos);
#define add_GenericIndexResult Wise2_add_GenericIndexResult


/* Function:  free_GenericIndexResult(p)
 *
 * Descrip:    Frees GenericIndexResults - overrides dynamite default
 *
 *
 * Arg:        p [UNKN ] Undocumented argument [GenericIndexResult *]
 *
 * Return [UNKN ]  Undocumented return value [GenericIndexResult *]
 *
 */
GenericIndexResult * Wise2_free_GenericIndexResult(GenericIndexResult * p);
#define free_GenericIndexResult Wise2_free_GenericIndexResult


/* Function:  hard_link_GenericIndexResult(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GenericIndexResult *]
 *
 * Return [UNKN ]  Undocumented return value [GenericIndexResult *]
 *
 */
GenericIndexResult * Wise2_hard_link_GenericIndexResult(GenericIndexResult * obj);
#define hard_link_GenericIndexResult Wise2_hard_link_GenericIndexResult


/* Function:  GenericIndexResult_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GenericIndexResult *]
 *
 */
GenericIndexResult * Wise2_GenericIndexResult_alloc(void);
#define GenericIndexResult_alloc Wise2_GenericIndexResult_alloc


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
