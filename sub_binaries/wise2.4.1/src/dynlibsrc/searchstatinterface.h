#ifndef DYNAMITEsearchstatinterfaceHEADERFILE
#define DYNAMITEsearchstatinterfaceHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "wisebase.h"
#include "sequence.h"


/* Object SearchStatInterface
 *
 * Descrip: SearchStatInterface converts raw scores into both
 *        bit score and evalues. Both must be supplied. The 
 *        function signatures are query_length,target_length,raw_score
 *
 *
 */
struct Wise2_SearchStatInterface {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    double (*calc_evalue)(void *,Sequence*,Sequence*,int,int);   
    double (*calc_bits)(void *,int,int,int); 
    char* (*attribution)(void *);    
    void (*free_data)(void *);   
    void * data;     
    } ;  
/* SearchStatInterface defined */ 
#ifndef DYNAMITE_DEFINED_SearchStatInterface
typedef struct Wise2_SearchStatInterface Wise2_SearchStatInterface;
#define SearchStatInterface Wise2_SearchStatInterface
#define DYNAMITE_DEFINED_SearchStatInterface
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_SearchStatInterface(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SearchStatInterface *]
 *
 * Return [UNKN ]  Undocumented return value [SearchStatInterface *]
 *
 */
SearchStatInterface * Wise2_hard_link_SearchStatInterface(SearchStatInterface * obj);
#define hard_link_SearchStatInterface Wise2_hard_link_SearchStatInterface


/* Function:  SearchStatInterface_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SearchStatInterface *]
 *
 */
SearchStatInterface * Wise2_SearchStatInterface_alloc(void);
#define SearchStatInterface_alloc Wise2_SearchStatInterface_alloc


/* Function:  free_SearchStatInterface(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SearchStatInterface *]
 *
 * Return [UNKN ]  Undocumented return value [SearchStatInterface *]
 *
 */
SearchStatInterface * Wise2_free_SearchStatInterface(SearchStatInterface * obj);
#define free_SearchStatInterface Wise2_free_SearchStatInterface


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
