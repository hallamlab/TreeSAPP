#ifndef DYNAMITEdprunimplHEADERFILE
#define DYNAMITEdprunimplHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "wisebase.h"
#include "basematrix.h"

typedef enum DPRunImplMemory {
  DPIM_Default = 543,
  DPIM_Explicit,
  DPIM_Linear 
} DPRunImplMemory; 


struct Wise2_DPRunImpl {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    DPRunImplMemory memory;  
    int kbyte_size;  
    boolean debug;   
    boolean paldebug;    
    boolean should_cache;    
    BaseMatrix * cache;  
    } ;  
/* DPRunImpl defined */ 
#ifndef DYNAMITE_DEFINED_DPRunImpl
typedef struct Wise2_DPRunImpl Wise2_DPRunImpl;
#define DPRunImpl Wise2_DPRunImpl
#define DYNAMITE_DEFINED_DPRunImpl
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  clone_DPRunImpl(dpri)
 *
 * Descrip:    Clones a DPRunImpl - particularly sensible
 *             for cached cases
 *
 *
 * Arg:        dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [DPRunImpl *]
 *
 */
DPRunImpl * Wise2_clone_DPRunImpl(DPRunImpl * dpri);
#define clone_DPRunImpl Wise2_clone_DPRunImpl


/* Function:  show_help_DPRunImpl(ofp)
 *
 * Descrip:    Shows help functions for DPRunImpl
 *
 *
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_show_help_DPRunImpl(FILE * ofp);
#define show_help_DPRunImpl Wise2_show_help_DPRunImpl


/* Function:  new_DPRunImpl_from_argv(argc,argv)
 *
 * Descrip:    Makes a DPRunImpl object from stripping from
 *             a command line
 *
 *
 * Arg:        argc [UNKN ] Undocumented argument [int *]
 * Arg:        argv [UNKN ] Undocumented argument [char **]
 *
 * Return [UNKN ]  Undocumented return value [DPRunImpl *]
 *
 */
DPRunImpl * Wise2_new_DPRunImpl_from_argv(int * argc,char ** argv);
#define new_DPRunImpl_from_argv Wise2_new_DPRunImpl_from_argv


/* Function:  hard_link_DPRunImpl(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [DPRunImpl *]
 *
 */
DPRunImpl * Wise2_hard_link_DPRunImpl(DPRunImpl * obj);
#define hard_link_DPRunImpl Wise2_hard_link_DPRunImpl


/* Function:  DPRunImpl_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DPRunImpl *]
 *
 */
DPRunImpl * Wise2_DPRunImpl_alloc(void);
#define DPRunImpl_alloc Wise2_DPRunImpl_alloc


/* Function:  free_DPRunImpl(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [DPRunImpl *]
 *
 */
DPRunImpl * Wise2_free_DPRunImpl(DPRunImpl * obj);
#define free_DPRunImpl Wise2_free_DPRunImpl


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
