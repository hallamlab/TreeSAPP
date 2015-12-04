#ifndef DYNAMITEfunctionclientHEADERFILE
#define DYNAMITEfunctionclientHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "transferinterface.h"
#include "anonobj.h"
#include "directsocketwrite.h"

#define FunctionProxyCoordinatorLISTLENGTH 16

struct Wise2_FunctionProxy {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    TransferedFunctionCall * transfer;   
    } ;  
/* FunctionProxy defined */ 
#ifndef DYNAMITE_DEFINED_FunctionProxy
typedef struct Wise2_FunctionProxy Wise2_FunctionProxy;
#define FunctionProxy Wise2_FunctionProxy
#define DYNAMITE_DEFINED_FunctionProxy
#endif


struct Wise2_FunctionProxyCoordinator {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    FunctionProxy ** proxy;  
    int len;/* len for above proxy  */ 
    int maxlen; /* maxlen for above proxy */ 
    int socket;  
    char * host;     
    int port;    
    } ;  
/* FunctionProxyCoordinator defined */ 
#ifndef DYNAMITE_DEFINED_FunctionProxyCoordinator
typedef struct Wise2_FunctionProxyCoordinator Wise2_FunctionProxyCoordinator;
#define FunctionProxyCoordinator Wise2_FunctionProxyCoordinator
#define DYNAMITE_DEFINED_FunctionProxyCoordinator
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  free_FunctionProxyCoordinator(obj)
 *
 * Descrip:    provides specific destructor for FunctionProxyCoordinators
 *
 *
 * Arg:        obj [UNKN ] Undocumented argument [FunctionProxyCoordinator *]
 *
 * Return [UNKN ]  Undocumented return value [FunctionProxyCoordinator *]
 *
 */
FunctionProxyCoordinator * Wise2_free_FunctionProxyCoordinator(FunctionProxyCoordinator * obj);
#define free_FunctionProxyCoordinator Wise2_free_FunctionProxyCoordinator


/* Function:  hard_link_FunctionProxy(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [FunctionProxy *]
 *
 * Return [UNKN ]  Undocumented return value [FunctionProxy *]
 *
 */
FunctionProxy * Wise2_hard_link_FunctionProxy(FunctionProxy * obj);
#define hard_link_FunctionProxy Wise2_hard_link_FunctionProxy


/* Function:  FunctionProxy_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [FunctionProxy *]
 *
 */
FunctionProxy * Wise2_FunctionProxy_alloc(void);
#define FunctionProxy_alloc Wise2_FunctionProxy_alloc


/* Function:  free_FunctionProxy(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [FunctionProxy *]
 *
 * Return [UNKN ]  Undocumented return value [FunctionProxy *]
 *
 */
FunctionProxy * Wise2_free_FunctionProxy(FunctionProxy * obj);
#define free_FunctionProxy Wise2_free_FunctionProxy


/* Function:  add_FunctionProxyCoordinator(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [FunctionProxyCoordinator *]
 * Arg:        add [OWNER] Object to add to the list [FunctionProxy *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_FunctionProxyCoordinator(FunctionProxyCoordinator * obj,FunctionProxy * add);
#define add_FunctionProxyCoordinator Wise2_add_FunctionProxyCoordinator


/* Function:  flush_FunctionProxyCoordinator(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [FunctionProxyCoordinator *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_FunctionProxyCoordinator(FunctionProxyCoordinator * obj);
#define flush_FunctionProxyCoordinator Wise2_flush_FunctionProxyCoordinator


/* Function:  FunctionProxyCoordinator_alloc_std(void)
 *
 * Descrip:    Equivalent to FunctionProxyCoordinator_alloc_len(FunctionProxyCoordinatorLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [FunctionProxyCoordinator *]
 *
 */
FunctionProxyCoordinator * Wise2_FunctionProxyCoordinator_alloc_std(void);
#define FunctionProxyCoordinator_alloc_std Wise2_FunctionProxyCoordinator_alloc_std


/* Function:  FunctionProxyCoordinator_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [FunctionProxyCoordinator *]
 *
 */
FunctionProxyCoordinator * Wise2_FunctionProxyCoordinator_alloc_len(int len);
#define FunctionProxyCoordinator_alloc_len Wise2_FunctionProxyCoordinator_alloc_len


/* Function:  hard_link_FunctionProxyCoordinator(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [FunctionProxyCoordinator *]
 *
 * Return [UNKN ]  Undocumented return value [FunctionProxyCoordinator *]
 *
 */
FunctionProxyCoordinator * Wise2_hard_link_FunctionProxyCoordinator(FunctionProxyCoordinator * obj);
#define hard_link_FunctionProxyCoordinator Wise2_hard_link_FunctionProxyCoordinator


/* Function:  FunctionProxyCoordinator_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [FunctionProxyCoordinator *]
 *
 */
FunctionProxyCoordinator * Wise2_FunctionProxyCoordinator_alloc(void);
#define FunctionProxyCoordinator_alloc Wise2_FunctionProxyCoordinator_alloc


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
AnonymousObject * Wise2_dispatch_FunctionProxy(FunctionProxyCoordinator * fpc,char * function_call,AnonymousObjectList * args);
#define dispatch_FunctionProxy Wise2_dispatch_FunctionProxy
FunctionProxy * Wise2_new_FunctionProxy(TransferedFunctionCall * t);
#define new_FunctionProxy Wise2_new_FunctionProxy
FunctionProxyCoordinator * Wise2_new_FunctionProxyCoordinator(char * host,int port);
#define new_FunctionProxyCoordinator Wise2_new_FunctionProxyCoordinator


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_FunctionProxyCoordinator(FunctionProxy ** list,int i,int j) ;
#define swap_FunctionProxyCoordinator Wise2_swap_FunctionProxyCoordinator
void Wise2_qsort_FunctionProxyCoordinator(FunctionProxy ** list,int left,int right,int (*comp)(FunctionProxy * ,FunctionProxy * ));
#define qsort_FunctionProxyCoordinator Wise2_qsort_FunctionProxyCoordinator
void Wise2_sort_FunctionProxyCoordinator(FunctionProxyCoordinator * obj,int (*comp)(FunctionProxy *, FunctionProxy *));
#define sort_FunctionProxyCoordinator Wise2_sort_FunctionProxyCoordinator
boolean Wise2_expand_FunctionProxyCoordinator(FunctionProxyCoordinator * obj,int len);
#define expand_FunctionProxyCoordinator Wise2_expand_FunctionProxyCoordinator

#ifdef _cplusplus
}
#endif

#endif
