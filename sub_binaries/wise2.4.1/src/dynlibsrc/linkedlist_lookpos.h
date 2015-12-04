#ifndef DYNAMITElinkedlist_lookposHEADERFILE
#define DYNAMITElinkedlist_lookposHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "seqlookup.h"


struct Wise2_SeqLookupPos {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Sequence * seq;  
    int pos;     
    struct Wise2_SeqLookupPos * next;    
    } ;  
/* SeqLookupPos defined */ 
#ifndef DYNAMITE_DEFINED_SeqLookupPos
typedef struct Wise2_SeqLookupPos Wise2_SeqLookupPos;
#define SeqLookupPos Wise2_SeqLookupPos
#define DYNAMITE_DEFINED_SeqLookupPos
#endif


struct Wise2_SeqLookupPosResult {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    SeqLookupPos * current;  
    SeqLookupResultStruct result;    
    } ;  
/* SeqLookupPosResult defined */ 
#ifndef DYNAMITE_DEFINED_SeqLookupPosResult
typedef struct Wise2_SeqLookupPosResult Wise2_SeqLookupPosResult;
#define SeqLookupPosResult Wise2_SeqLookupPosResult
#define DYNAMITE_DEFINED_SeqLookupPosResult
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  new_linkedl_SeqLookupResultInterface(head)
 *
 * Descrip:    Makes a hash based SeqLookupResultsInterface thing
 *
 *
 * Arg:        head [UNKN ] Undocumented argument [SeqLookupPos *]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupResultInterface *]
 *
 */
SeqLookupResultInterface * Wise2_new_linkedl_SeqLookupResultInterface(SeqLookupPos * head);
#define new_linkedl_SeqLookupResultInterface Wise2_new_linkedl_SeqLookupResultInterface


/* Function:  free_linkedl_SeqLook(data)
 *
 * Descrip:    Internal function returns data...
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 *
 */
void Wise2_free_linkedl_SeqLook(void * data);
#define free_linkedl_SeqLook Wise2_free_linkedl_SeqLook


/* Function:  is_more_linkedl_SeqLook(*data)
 *
 * Descrip:    Internal function for returning whether there is more data
 *
 *
 * Arg:        *data [UNKN ] Undocumented argument [void]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_is_more_linkedl_SeqLook(void *data);
#define is_more_linkedl_SeqLook Wise2_is_more_linkedl_SeqLook


/* Function:  next_linkedl_SeqLook(data,prev)
 *
 * Descrip:    Internal function for returning the next position in the hash
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 * Arg:        prev [UNKN ] Undocumented argument [SeqLookupResultStruct *]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupResultStruct *]
 *
 */
SeqLookupResultStruct * Wise2_next_linkedl_SeqLook(void * data,SeqLookupResultStruct * prev);
#define next_linkedl_SeqLook Wise2_next_linkedl_SeqLook


/* Function:  hard_link_SeqLookupPos(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SeqLookupPos *]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupPos *]
 *
 */
SeqLookupPos * Wise2_hard_link_SeqLookupPos(SeqLookupPos * obj);
#define hard_link_SeqLookupPos Wise2_hard_link_SeqLookupPos


/* Function:  SeqLookupPos_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupPos *]
 *
 */
SeqLookupPos * Wise2_SeqLookupPos_alloc(void);
#define SeqLookupPos_alloc Wise2_SeqLookupPos_alloc


/* Function:  free_SeqLookupPos(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SeqLookupPos *]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupPos *]
 *
 */
SeqLookupPos * Wise2_free_SeqLookupPos(SeqLookupPos * obj);
#define free_SeqLookupPos Wise2_free_SeqLookupPos


/* Function:  hard_link_SeqLookupPosResult(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SeqLookupPosResult *]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupPosResult *]
 *
 */
SeqLookupPosResult * Wise2_hard_link_SeqLookupPosResult(SeqLookupPosResult * obj);
#define hard_link_SeqLookupPosResult Wise2_hard_link_SeqLookupPosResult


/* Function:  SeqLookupPosResult_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupPosResult *]
 *
 */
SeqLookupPosResult * Wise2_SeqLookupPosResult_alloc(void);
#define SeqLookupPosResult_alloc Wise2_SeqLookupPosResult_alloc


/* Function:  free_SeqLookupPosResult(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SeqLookupPosResult *]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupPosResult *]
 *
 */
SeqLookupPosResult * Wise2_free_SeqLookupPosResult(SeqLookupPosResult * obj);
#define free_SeqLookupPosResult Wise2_free_SeqLookupPosResult


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
