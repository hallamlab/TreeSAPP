#ifdef _cplusplus
extern "C" {
#endif
#include "searchstatinterface.h"


# line 6 "searchstatinterface.c"
/* Function:  hard_link_SearchStatInterface(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SearchStatInterface *]
 *
 * Return [UNKN ]  Undocumented return value [SearchStatInterface *]
 *
 */
SearchStatInterface * hard_link_SearchStatInterface(SearchStatInterface * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a SearchStatInterface object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  SearchStatInterface_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SearchStatInterface *]
 *
 */
SearchStatInterface * SearchStatInterface_alloc(void) 
{
    SearchStatInterface * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(SearchStatInterface *) ckalloc (sizeof(SearchStatInterface))) == NULL)  {  
      warn("SearchStatInterface_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->calc_evalue = NULL; 
    out->calc_bits = NULL;   
    out->attribution = NULL; 
    out->free_data = NULL;   


    return out;  
}    


/* Function:  free_SearchStatInterface(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SearchStatInterface *]
 *
 * Return [UNKN ]  Undocumented return value [SearchStatInterface *]
 *
 */
SearchStatInterface * free_SearchStatInterface(SearchStatInterface * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a SearchStatInterface obj. Should be trappable");   
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
    /* obj->calc_evalue is a function pointer */ 
    /* obj->calc_bits is a function pointer */ 
    /* obj->attribution is a function pointer */ 
    /* obj->free_data is a function pointer */ 
    /* obj->data is linked in */ 


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
