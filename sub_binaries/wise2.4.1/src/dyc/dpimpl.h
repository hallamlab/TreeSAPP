#ifndef DYNAMITEdpimplHEADERFILE
#define DYNAMITEdpimplHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "wisebase.h"


struct DycWarning {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    boolean warn_extern;     
    boolean warn_extern_method;  
    boolean warn_c_type;     
    } ;  
/* DycWarning defined */ 
#ifndef DYNAMITE_DEFINED_DycWarning
typedef struct DycWarning DycWarning;
#define DYNAMITE_DEFINED_DycWarning
#endif


/* Object DPImplementation
 *
 * Descrip: No Description
 *
 */
struct DPImplementation {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    boolean do_threads;  
    int protect_level;   
    int db_trace_level;  
    boolean doprob;  
    boolean doone;   
    char * calcfunc;     
    boolean largemem;    
    DycWarning * dycw;   
    boolean dydebug;     
    } ;  
/* DPImplementation defined */ 
#ifndef DYNAMITE_DEFINED_DPImplementation
typedef struct DPImplementation DPImplementation;
#define DYNAMITE_DEFINED_DPImplementation
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  show_help_DycWarning(ofp)
 *
 * Descrip:    Shows to stdout the list of options
 *             used by DycWarning
 *
 *
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void show_help_DycWarning(FILE * ofp);


/* Function:  show_help_DPImplementation(ofp)
 *
 * Descrip:    Shows to stdout the list options used
 *             by DPImplementation
 *
 *
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void show_help_DPImplementation(FILE * ofp);


/* Function:  new_DycWarning_from_argstr(argc,argv)
 *
 * Descrip:    Processes the argstring into the DycWarning stuff
 *
 *
 * Arg:        argc [UNKN ] Undocumented argument [int *]
 * Arg:        argv [UNKN ] Undocumented argument [char **]
 *
 * Return [UNKN ]  Undocumented return value [DycWarning *]
 *
 */
DycWarning * new_DycWarning_from_argstr(int * argc,char ** argv);


/* Function:  new_DPImplementation_from_argstr(argc,argv)
 *
 * Descrip:    Processes the argstring into the DPImplementation
 *             datastructure.
 *
 *
 * Arg:        argc [UNKN ] Undocumented argument [int *]
 * Arg:        argv [UNKN ] Undocumented argument [char **]
 *
 * Return [UNKN ]  Undocumented return value [DPImplementation *]
 *
 */
DPImplementation * new_DPImplementation_from_argstr(int * argc,char ** argv);


/* Function:  hard_link_DycWarning(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DycWarning *]
 *
 * Return [UNKN ]  Undocumented return value [DycWarning *]
 *
 */
DycWarning * hard_link_DycWarning(DycWarning * obj);


/* Function:  DycWarning_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DycWarning *]
 *
 */
DycWarning * DycWarning_alloc(void);


/* Function:  free_DycWarning(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DycWarning *]
 *
 * Return [UNKN ]  Undocumented return value [DycWarning *]
 *
 */
DycWarning * free_DycWarning(DycWarning * obj);


/* Function:  hard_link_DPImplementation(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DPImplementation *]
 *
 * Return [UNKN ]  Undocumented return value [DPImplementation *]
 *
 */
DPImplementation * hard_link_DPImplementation(DPImplementation * obj);


/* Function:  DPImplementation_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DPImplementation *]
 *
 */
DPImplementation * DPImplementation_alloc(void);


/* Function:  free_DPImplementation(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DPImplementation *]
 *
 * Return [UNKN ]  Undocumented return value [DPImplementation *]
 *
 */
DPImplementation * free_DPImplementation(DPImplementation * obj);


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
