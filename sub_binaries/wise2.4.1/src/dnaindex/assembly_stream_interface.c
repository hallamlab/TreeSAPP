#ifdef _cplusplus
extern "C" {
#endif
#include "assembly_stream_interface.h"


/* Function:  free_AssemblySequenceStream(obj)
 *
 * Descrip:    provides specific deconstructor 
 *
 *
 * Arg:        obj [UNKN ] Undocumented argument [AssemblySequenceStream *]
 *
 * Return [UNKN ]  Undocumented return value [AssemblySequenceStream *]
 *
 */
# line 23 "assembly_stream_interface.dy"
AssemblySequenceStream * free_AssemblySequenceStream(AssemblySequenceStream * obj)
{
  int return_early = 0;    
  
  
  if( obj == NULL) {  
    warn("Attempting to free a NULL pointer to a AssemblySequenceStream obj. Should be trappable");    
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
    

  if( obj->handle != NULL ) {
    if( obj->free_handle == NULL ) {
      warn("In assembly stream constructor, no free function for handle. Probably will leak memory");
    } else {
      (*obj->free_handle)(obj->handle);
    }
  }

  free(obj);

  return NULL;
}

# line 52 "assembly_stream_interface.c"
/* Function:  hard_link_AssemblySequenceStream(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AssemblySequenceStream *]
 *
 * Return [UNKN ]  Undocumented return value [AssemblySequenceStream *]
 *
 */
AssemblySequenceStream * hard_link_AssemblySequenceStream(AssemblySequenceStream * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a AssemblySequenceStream object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  AssemblySequenceStream_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AssemblySequenceStream *]
 *
 */
AssemblySequenceStream * AssemblySequenceStream_alloc(void) 
{
    AssemblySequenceStream * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(AssemblySequenceStream *) ckalloc (sizeof(AssemblySequenceStream))) == NULL)    {  
      warn("AssemblySequenceStream_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->next_AssemblySequence = NULL;   
    out->free_handle = NULL; 


    return out;  
}    



#ifdef _cplusplus
}
#endif
