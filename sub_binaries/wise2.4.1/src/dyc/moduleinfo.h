#ifndef DYNAMITEmoduleinfoHEADERFILE
#define DYNAMITEmoduleinfoHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "ftext.h"



struct ModuleInfo {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Ftext * ft; /*  description of module */ 
    } ;  
/* ModuleInfo defined */ 
#ifndef DYNAMITE_DEFINED_ModuleInfo
typedef struct ModuleInfo ModuleInfo;
#define DYNAMITE_DEFINED_ModuleInfo
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_ModuleInfo(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ModuleInfo *]
 *
 * Return [UNKN ]  Undocumented return value [ModuleInfo *]
 *
 */
ModuleInfo * hard_link_ModuleInfo(ModuleInfo * obj);


/* Function:  ModuleInfo_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ModuleInfo *]
 *
 */
ModuleInfo * ModuleInfo_alloc(void);


/* Function:  free_ModuleInfo(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ModuleInfo *]
 *
 * Return [UNKN ]  Undocumented return value [ModuleInfo *]
 *
 */
ModuleInfo * free_ModuleInfo(ModuleInfo * obj);


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
