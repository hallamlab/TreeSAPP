#ifndef DYNAMITEdprunimplHEADERFILE
#define DYNAMITEdprunimplHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "wisebase.h"


typedef enum DPRunImplMemory {
  DPIM_Default = 543,
  DPIM_Explicit,
  DPIM_Linear 
} DPRunImplMemory; 


struct Wise2_DPRunImpl {  
    int dynamite_hard_link;  
    DPRunImplMemory memory;  
    int kbyte_size;  
    boolean debug;   
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
