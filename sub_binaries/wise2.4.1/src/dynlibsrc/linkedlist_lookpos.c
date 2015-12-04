#ifdef _cplusplus
extern "C" {
#endif
#include "linkedlist_lookpos.h"




/* Function:  new_linkedl_SeqLookupResultInterface(head)
 *
 * Descrip:    Makes a hash based SeqLookupResultsInterface thing
 *
 *
 * Arg:        head [UNKN ] Undocumented argument [SeqLookupPos *]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupResultInterface *]
 *
 */
# line 27 "linkedlist_lookpos.dy"
SeqLookupResultInterface * new_linkedl_SeqLookupResultInterface(SeqLookupPos * head)
{
  SeqLookupResultInterface * out;
  SeqLookupPosResult * data;

  data = SeqLookupPosResult_alloc();
  data->current = head;

  out = SeqLookupResultInterface_alloc();

  out->next = next_linkedl_SeqLook;
  out->is_more = is_more_linkedl_SeqLook;
  out->free_data = free_linkedl_SeqLook;
  out->data = (void*)data;

  return out;
  
}


/* Function:  free_linkedl_SeqLook(data)
 *
 * Descrip:    Internal function returns data...
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 *
 */
# line 50 "linkedlist_lookpos.dy"
void free_linkedl_SeqLook(void * data)
{
  SeqLookupPosResult * posres = (SeqLookupPosResult *) data;


  free_SeqLookupPosResult(posres);

}  

/* Function:  is_more_linkedl_SeqLook(*data)
 *
 * Descrip:    Internal function for returning whether there is more data
 *
 *
 * Arg:        *data [UNKN ] Undocumented argument [void]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 62 "linkedlist_lookpos.dy"
boolean is_more_linkedl_SeqLook(void *data)
{
  SeqLookupPosResult * posres = (SeqLookupPosResult *) data;

  if( posres->current == NULL ) {
    return FALSE;
  } else {
    return TRUE;
  }

}


/* Function:  next_linkedl_SeqLook(data,prev)
 *
 * Descrip:    Internal function for returning the next position in the hash
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 * Arg:        prev [UNKN ] Undocumented argument [SeqLookupResultStruct *]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupResultStruct *]
 *
 */
# line 78 "linkedlist_lookpos.dy"
SeqLookupResultStruct * next_linkedl_SeqLook(void * data,SeqLookupResultStruct * prev)
{
  SeqLookupPosResult * posres = (SeqLookupPosResult *) data;

  if( posres->current == NULL ) {
    fatal("Overrun virtual buffer");
  }

  posres->result.seq = posres->current->seq;
  posres->result.pos = posres->current->pos;
  
  posres->current = posres->current->next;

  return &(posres->result);
}
# line 100 "linkedlist_lookpos.c"
/* Function:  hard_link_SeqLookupPos(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SeqLookupPos *]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupPos *]
 *
 */
SeqLookupPos * hard_link_SeqLookupPos(SeqLookupPos * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a SeqLookupPos object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  SeqLookupPos_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupPos *]
 *
 */
SeqLookupPos * SeqLookupPos_alloc(void) 
{
    SeqLookupPos * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(SeqLookupPos *) ckalloc (sizeof(SeqLookupPos))) == NULL)    {  
      warn("SeqLookupPos_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->seq = NULL; 
    out->pos = 0;    
    out->next = NULL;    


    return out;  
}    


/* Function:  free_SeqLookupPos(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SeqLookupPos *]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupPos *]
 *
 */
SeqLookupPos * free_SeqLookupPos(SeqLookupPos * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a SeqLookupPos obj. Should be trappable");  
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
    if( obj->seq != NULL)    
      free_Sequence(obj->seq);   
    /* Unable to make free function for obj->next */ 


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_SeqLookupPosResult(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SeqLookupPosResult *]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupPosResult *]
 *
 */
SeqLookupPosResult * hard_link_SeqLookupPosResult(SeqLookupPosResult * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a SeqLookupPosResult object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  SeqLookupPosResult_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupPosResult *]
 *
 */
SeqLookupPosResult * SeqLookupPosResult_alloc(void) 
{
    SeqLookupPosResult * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(SeqLookupPosResult *) ckalloc (sizeof(SeqLookupPosResult))) == NULL)    {  
      warn("SeqLookupPosResult_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->current = NULL; 


    return out;  
}    


/* Function:  free_SeqLookupPosResult(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SeqLookupPosResult *]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupPosResult *]
 *
 */
SeqLookupPosResult * free_SeqLookupPosResult(SeqLookupPosResult * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a SeqLookupPosResult obj. Should be trappable");    
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
    if( obj->current != NULL)    
      free_SeqLookupPos(obj->current);   
    /* obj->result is linked in */ 


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
