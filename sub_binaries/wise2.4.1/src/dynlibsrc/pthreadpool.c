#ifdef _cplusplus
extern "C" {
#endif
#include "pthreadpool.h"


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
# line 51 "pthreadpool.dy"
PThreadPool * new_PThreadPool(int no_threads)
{
  PThreadPool * out;

  out = PThreadPool_alloc();

  out->number_of_threads = no_threads;
  out->max_work_size = max_work;
  out->cur_queue_size = 0;
  out->head = NULL;
  out->tail = NULL;
  out->queue_closed = 0;
  out->shutdown = 0;

  out->lock = (pthread_mutex_t *) ckalloc (sizeof(pthread_mutex_t));
  out->work_to_do = (pthread_cond_t *) ckalloc (sizeof(pthread_cond_t));
  out->queue_not_full = (pthread_cond_t *) ckalloc (sizeof(pthread_cond_t));
  out->queue_empty = (pthread_cond_t *) ckalloc (sizeof(pthread_cond_t));

  if( pthread_mutex_init(out->lock,NULL) != 0 ) {
    warn("Unable to initialise lock mutex");
    return NULL;
  }

  return out;

}


# line 44 "pthreadpool.c"
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
PTP_Work * hard_link_PTP_Work(PTP_Work * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a PTP_Work object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  PTP_Work_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PTP_Work *]
 *
 */
PTP_Work * PTP_Work_alloc(void) 
{
    PTP_Work * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(PTP_Work *) ckalloc (sizeof(PTP_Work))) == NULL)    {  
      warn("PTP_Work_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->work_routine = NULL;    
    out->data = NULL;    
    out->next = NULL;    


    return out;  
}    


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
PTP_Work * free_PTP_Work(PTP_Work * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a PTP_Work obj. Should be trappable");  
      return NULL;   
      }  


#ifdef PTHREAD   
    assert(pthread_mutex_lock(&(obj->dynamite_mutex)) == 0); 
#endif   
    if( obj->dynamite_hard_link > 1)     {  
      return_early = 1;  
      obj->dynamite_hard_link--; 
      }  
#ifdef PTHREAD   
    assert(pthread_mutex_unlock(&(obj->dynamite_mutex)) == 0);   
#endif   
    if( return_early == 1)   
      return NULL;   
    /* obj->work_routine is a function pointer */ 
    if( obj->data != NULL)   
      free_void(obj->data);  
    if( obj->next != NULL)   
      free_PTP_Work(obj->next);  


    ckfree(obj); 
    return NULL; 
}    


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
PThreadPool * hard_link_PThreadPool(PThreadPool * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a PThreadPool object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  PThreadPool_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PThreadPool *]
 *
 */
PThreadPool * PThreadPool_alloc(void) 
{
    PThreadPool * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(PThreadPool *) ckalloc (sizeof(PThreadPool))) == NULL)  {  
      warn("PThreadPool_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->number_of_threads = 0;  
    /* threads[MAX_THREAD_NUMBER] is an array: no default possible */ 
    out->max_work_size = 0;  
    out->current_work_size = 0;  
    out->head = NULL;    
    out->tail = NULL;    
    out->lock = NULL;    
    out->work_to_do = NULL;  
    out->queue_not_full = NULL;  
    out->queue_empty = NULL; 
    out->queue_closed = 0;   
    out->shutdown = 0;   


    return out;  
}    


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
PThreadPool * free_PThreadPool(PThreadPool * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a PThreadPool obj. Should be trappable");   
      return NULL;   
      }  


#ifdef PTHREAD   
    assert(pthread_mutex_lock(&(obj->dynamite_mutex)) == 0); 
#endif   
    if( obj->dynamite_hard_link > 1)     {  
      return_early = 1;  
      obj->dynamite_hard_link--; 
      }  
#ifdef PTHREAD   
    assert(pthread_mutex_unlock(&(obj->dynamite_mutex)) == 0);   
#endif   
    if( return_early == 1)   
      return NULL;   
    if( obj->head != NULL)   
      free_PTP_Work(obj->head);  
    if( obj->tail != NULL)   
      free_PTP_Work(obj->tail);  
    if( obj->lock != NULL)   
      free_pthread_mutex_t(obj->lock);   
    if( obj->work_to_do != NULL) 
      free_pthread_cond_t(obj->work_to_do);  
    if( obj->queue_not_full != NULL) 
      free_pthread_cond_t(obj->queue_not_full);  
    if( obj->queue_empty != NULL)    
      free_pthread_cond_t(obj->queue_empty);     


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
