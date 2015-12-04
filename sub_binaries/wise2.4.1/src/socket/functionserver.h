#ifndef DYNAMITEfunctionserverHEADERFILE
#define DYNAMITEfunctionserverHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "wisebase.h"
#include "simplebufferedserver.h"
#include "transferinterface.h"
#include "directsocketwrite.h"
#include "anonobj.h"

#define FunctionServerLISTLENGTH 16



struct Wise2_FunctionImplementation {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    TransferedFunctionCall * transfer;   
    AnonymousObject* (*implementation)(void * h,AnonymousObjectList*);   
    void * handle;   
    } ;  
/* FunctionImplementation defined */ 
#ifndef DYNAMITE_DEFINED_FunctionImplementation
typedef struct Wise2_FunctionImplementation Wise2_FunctionImplementation;
#define FunctionImplementation Wise2_FunctionImplementation
#define DYNAMITE_DEFINED_FunctionImplementation
#endif


struct Wise2_FunctionServer {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    FunctionImplementation ** fi;    
    int len;/* len for above fi  */ 
    int maxlen; /* maxlen for above fi */ 
    int socket;  
    } ;  
/* FunctionServer defined */ 
#ifndef DYNAMITE_DEFINED_FunctionServer
typedef struct Wise2_FunctionServer Wise2_FunctionServer;
#define FunctionServer Wise2_FunctionServer
#define DYNAMITE_DEFINED_FunctionServer
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_FunctionImplementation(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [FunctionImplementation *]
 *
 * Return [UNKN ]  Undocumented return value [FunctionImplementation *]
 *
 */
FunctionImplementation * Wise2_hard_link_FunctionImplementation(FunctionImplementation * obj);
#define hard_link_FunctionImplementation Wise2_hard_link_FunctionImplementation


/* Function:  FunctionImplementation_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [FunctionImplementation *]
 *
 */
FunctionImplementation * Wise2_FunctionImplementation_alloc(void);
#define FunctionImplementation_alloc Wise2_FunctionImplementation_alloc


/* Function:  free_FunctionImplementation(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [FunctionImplementation *]
 *
 * Return [UNKN ]  Undocumented return value [FunctionImplementation *]
 *
 */
FunctionImplementation * Wise2_free_FunctionImplementation(FunctionImplementation * obj);
#define free_FunctionImplementation Wise2_free_FunctionImplementation


/* Function:  add_FunctionServer(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [FunctionServer *]
 * Arg:        add [OWNER] Object to add to the list [FunctionImplementation *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_FunctionServer(FunctionServer * obj,FunctionImplementation * add);
#define add_FunctionServer Wise2_add_FunctionServer


/* Function:  flush_FunctionServer(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [FunctionServer *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_FunctionServer(FunctionServer * obj);
#define flush_FunctionServer Wise2_flush_FunctionServer


/* Function:  FunctionServer_alloc_std(void)
 *
 * Descrip:    Equivalent to FunctionServer_alloc_len(FunctionServerLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [FunctionServer *]
 *
 */
FunctionServer * Wise2_FunctionServer_alloc_std(void);
#define FunctionServer_alloc_std Wise2_FunctionServer_alloc_std


/* Function:  FunctionServer_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [FunctionServer *]
 *
 */
FunctionServer * Wise2_FunctionServer_alloc_len(int len);
#define FunctionServer_alloc_len Wise2_FunctionServer_alloc_len


/* Function:  hard_link_FunctionServer(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [FunctionServer *]
 *
 * Return [UNKN ]  Undocumented return value [FunctionServer *]
 *
 */
FunctionServer * Wise2_hard_link_FunctionServer(FunctionServer * obj);
#define hard_link_FunctionServer Wise2_hard_link_FunctionServer


/* Function:  FunctionServer_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [FunctionServer *]
 *
 */
FunctionServer * Wise2_FunctionServer_alloc(void);
#define FunctionServer_alloc Wise2_FunctionServer_alloc


/* Function:  free_FunctionServer(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [FunctionServer *]
 *
 * Return [UNKN ]  Undocumented return value [FunctionServer *]
 *
 */
FunctionServer * Wise2_free_FunctionServer(FunctionServer * obj);
#define free_FunctionServer Wise2_free_FunctionServer


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
void Wise2_main_loop_forking_FunctionServer(FunctionServer * fs,int verbose);
#define main_loop_forking_FunctionServer Wise2_main_loop_forking_FunctionServer
FunctionServer * Wise2_new_FunctionServer(int port);
#define new_FunctionServer Wise2_new_FunctionServer
FunctionImplementation * Wise2_test_strcat_FunctionImplementation(void);
#define test_strcat_FunctionImplementation Wise2_test_strcat_FunctionImplementation
AnonymousObject * Wise2_test_strcat_implementation(void * h,AnonymousObjectList * aol);
#define test_strcat_implementation Wise2_test_strcat_implementation


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_FunctionServer(FunctionImplementation ** list,int i,int j) ;
#define swap_FunctionServer Wise2_swap_FunctionServer
void Wise2_qsort_FunctionServer(FunctionImplementation ** list,int left,int right,int (*comp)(FunctionImplementation * ,FunctionImplementation * ));
#define qsort_FunctionServer Wise2_qsort_FunctionServer
void Wise2_sort_FunctionServer(FunctionServer * obj,int (*comp)(FunctionImplementation *, FunctionImplementation *));
#define sort_FunctionServer Wise2_sort_FunctionServer
boolean Wise2_expand_FunctionServer(FunctionServer * obj,int len);
#define expand_FunctionServer Wise2_expand_FunctionServer

#ifdef _cplusplus
}
#endif

#endif
