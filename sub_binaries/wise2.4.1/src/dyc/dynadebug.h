#ifndef DYNAMITEdynadebugHEADERFILE
#define DYNAMITEdynadebugHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif

#include "dyna2.h"
#include "dynafunc.h"
#include "dpimpl.h"

#define ExprSetLISTLENGTH 50

struct ExprSet {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    ExprTree ** expr;    
    int len;/* len for above expr  */ 
    int maxlen; /* maxlen for above expr */ 
    } ;  
/* ExprSet defined */ 
#ifndef DYNAMITE_DEFINED_ExprSet
typedef struct ExprSet ExprSet;
#define DYNAMITE_DEFINED_ExprSet
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  add_ExprSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [ExprSet *]
 * Arg:        add [OWNER] Object to add to the list [ExprTree *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean add_ExprSet(ExprSet * obj,ExprTree * add);


/* Function:  flush_ExprSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [ExprSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_ExprSet(ExprSet * obj);


/* Function:  ExprSet_alloc_std(void)
 *
 * Descrip:    Equivalent to ExprSet_alloc_len(ExprSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ExprSet *]
 *
 */
ExprSet * ExprSet_alloc_std(void);


/* Function:  ExprSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [ExprSet *]
 *
 */
ExprSet * ExprSet_alloc_len(int len);


/* Function:  hard_link_ExprSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ExprSet *]
 *
 * Return [UNKN ]  Undocumented return value [ExprSet *]
 *
 */
ExprSet * hard_link_ExprSet(ExprSet * obj);


/* Function:  ExprSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ExprSet *]
 *
 */
ExprSet * ExprSet_alloc(void);


/* Function:  free_ExprSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ExprSet *]
 *
 * Return [UNKN ]  Undocumented return value [ExprSet *]
 *
 */
ExprSet * free_ExprSet(ExprSet * obj);


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
void write_debug_funcs(DYNFILE * dfp,GenericMatrix * gm);
ExprSet * build_ExprSet(ExprTree * root);
boolean descend_ExprTree_for_ExprSet(ExprSet * set,ExprTree * t);
void make_debug_struct_func(DYNFILE * dfp,GenericMatrix * gm);
void state_debug_func(DYNFILE * dfp,GenericMatrix * gm);
void transition_debug_func(DYNFILE * dfp,GenericMatrix * gm);
void explicit_debug_func(DYNFILE * dfp,GenericMatrix * gm);


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void swap_ExprSet(ExprTree ** list,int i,int j) ;
void qsort_ExprSet(ExprTree ** list,int left,int right,int (*comp)(ExprTree * ,ExprTree * ));
void sort_ExprSet(ExprSet * obj,int (*comp)(ExprTree *, ExprTree *));
boolean expand_ExprSet(ExprSet * obj,int len);

#ifdef _cplusplus
}
#endif

#endif
