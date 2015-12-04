#ifndef DYNAMITEstandardoutHEADERFILE
#define DYNAMITEstandardoutHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"


struct Wise2_StandardOutputOptions {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    boolean show_alb;    
    boolean show_pal;    
    boolean show_cumlative_alb;  
    boolean show_cumlative_pal;  
    } ;  
/* StandardOutputOptions defined */ 
#ifndef DYNAMITE_DEFINED_StandardOutputOptions
typedef struct Wise2_StandardOutputOptions Wise2_StandardOutputOptions;
#define StandardOutputOptions Wise2_StandardOutputOptions
#define DYNAMITE_DEFINED_StandardOutputOptions
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_StandardOutputOptions(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [StandardOutputOptions *]
 *
 * Return [UNKN ]  Undocumented return value [StandardOutputOptions *]
 *
 */
StandardOutputOptions * Wise2_hard_link_StandardOutputOptions(StandardOutputOptions * obj);
#define hard_link_StandardOutputOptions Wise2_hard_link_StandardOutputOptions


/* Function:  StandardOutputOptions_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [StandardOutputOptions *]
 *
 */
StandardOutputOptions * Wise2_StandardOutputOptions_alloc(void);
#define StandardOutputOptions_alloc Wise2_StandardOutputOptions_alloc


/* Function:  free_StandardOutputOptions(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [StandardOutputOptions *]
 *
 * Return [UNKN ]  Undocumented return value [StandardOutputOptions *]
 *
 */
StandardOutputOptions * Wise2_free_StandardOutputOptions(StandardOutputOptions * obj);
#define free_StandardOutputOptions Wise2_free_StandardOutputOptions


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
void Wise2_show_help_StandardOutputOptions(FILE * ofp);
#define show_help_StandardOutputOptions Wise2_show_help_StandardOutputOptions
StandardOutputOptions * Wise2_new_StandardOutputOptions_from_argv(int * argc,char ** argv);
#define new_StandardOutputOptions_from_argv Wise2_new_StandardOutputOptions_from_argv
void Wise2_show_StandardOutputOptions(StandardOutputOptions * out,AlnBlock * alb,PackAln * pal,char * divide_str,FILE * ofp);
#define show_StandardOutputOptions Wise2_show_StandardOutputOptions


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/

#ifdef _cplusplus
}
#endif

#endif
