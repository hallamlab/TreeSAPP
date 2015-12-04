#ifndef DYNAMITEpthreadpoolHEADERFILE
#define DYNAMITEpthreadpoolHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "wisebase.h"
#include <pthread.h>


#define MAX_THREAD_NUMBER 100

#ifndef DYNAMITE_DEFINED_PThreadPool
typedef struct Wise2_PThreadPool Wise2_PThreadPool;
#define PThreadPool Wise2_PThreadPool
#define DYNAMITE_DEFINED_PThreadPool
#endif

#ifndef DYNAMITE_DEFINED_PTP_Work
typedef struct Wise2_PTP_Work Wise2_PTP_Work;
#define PTP_Work Wise2_PTP_Work
#define DYNAMITE_DEFINED_PTP_Work
#endif

/* Object PTP_Work
 *
 * Descrip: This data structure holds an individual job to be
 *        executed by the PThreadPool.
 *
 *
 */
struct Wise2_PTP_Work {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    void (*work_routine)(void * data);   
    void * data;     
    PTP_Work * next;     
    } ;  
/* PTP_Work defined */ 
#ifndef DYNAMITE_DEFINED_PTP_Work
typedef struct Wise2_PTP_Work Wise2_PTP_Work;
#define PTP_Work Wise2_PTP_Work
#define DYNAMITE_DEFINED_PTP_Work
#endif


/* Object PThreadPool
 *
 * Descrip: This datastructure is to hold a thread pool.
 *        Work can be added via 
 *
 *        This is work in progress. Dont use!
 *
 *
 */
struct Wise2_PThreadPool {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int number_of_threads;   
    pthread_t threads[MAX_THREAD_NUMBER];    
    int max_work_size;   
    int current_work_size;   
    PTP_Work * head;     
    PTP_Work * tail;     
    pthread_mutex_t * lock;  
    pthread_cond_t  * work_to_do;    
    pthread_cond_t  * queue_not_full;    
    pthread_cond_t  * queue_empty;   
    int queue_closed;    
    int shutdown;    
    } ;  
/* PThreadPool defined */ 
#ifndef DYNAMITE_DEFINED_PThreadPool
typedef struct Wise2_PThreadPool Wise2_PThreadPool;
#define PThreadPool Wise2_PThreadPool
#define DYNAMITE_DEFINED_PThreadPool
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  new_PThreadPool(no_threads)
 *
 * Descrip:    Makes a new Thread pool 
 *
 *
 * Arg:        no_threads [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [PThreadPool *]
 *
 */
PThreadPool * Wise2_new_PThreadPool(int no_threads);
#define new_PThreadPool Wise2_new_PThreadPool


/* Function:  hard_link_PTP_Work(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [PTP_Work *]
 *
 * Return [UNKN ]  Undocumented return value [PTP_Work *]
 *
 */
PTP_Work * Wise2_hard_link_PTP_Work(PTP_Work * obj);
#define hard_link_PTP_Work Wise2_hard_link_PTP_Work


/* Function:  PTP_Work_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PTP_Work *]
 *
 */
PTP_Work * Wise2_PTP_Work_alloc(void);
#define PTP_Work_alloc Wise2_PTP_Work_alloc


/* Function:  free_PTP_Work(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [PTP_Work *]
 *
 * Return [UNKN ]  Undocumented return value [PTP_Work *]
 *
 */
PTP_Work * Wise2_free_PTP_Work(PTP_Work * obj);
#define free_PTP_Work Wise2_free_PTP_Work


/* Function:  hard_link_PThreadPool(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [PThreadPool *]
 *
 * Return [UNKN ]  Undocumented return value [PThreadPool *]
 *
 */
PThreadPool * Wise2_hard_link_PThreadPool(PThreadPool * obj);
#define hard_link_PThreadPool Wise2_hard_link_PThreadPool


/* Function:  PThreadPool_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PThreadPool *]
 *
 */
PThreadPool * Wise2_PThreadPool_alloc(void);
#define PThreadPool_alloc Wise2_PThreadPool_alloc


/* Function:  free_PThreadPool(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [PThreadPool *]
 *
 * Return [UNKN ]  Undocumented return value [PThreadPool *]
 *
 */
PThreadPool * Wise2_free_PThreadPool(PThreadPool * obj);
#define free_PThreadPool Wise2_free_PThreadPool


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
