#ifdef _cplusplus
extern "C" {
#endif
#include "variable.h"



# line 7 "variable.c"
/* Function:  hard_link_MatrixVariable(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [MatrixVariable *]
 *
 * Return [UNKN ]  Undocumented return value [MatrixVariable *]
 *
 */
MatrixVariable * hard_link_MatrixVariable(MatrixVariable * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a MatrixVariable object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  MatrixVariable_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [MatrixVariable *]
 *
 */
MatrixVariable * MatrixVariable_alloc(void) 
{
    MatrixVariable * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(MatrixVariable *) ckalloc (sizeof(MatrixVariable))) == NULL)    {  
      warn("MatrixVariable_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->name = NULL;    
    out->type = 0;   
    out->dim1 = 0;   
    out->dim2 = 0;   
    out->source = NULL;  


    return out;  
}    


/* Function:  free_MatrixVariable(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [MatrixVariable *]
 *
 * Return [UNKN ]  Undocumented return value [MatrixVariable *]
 *
 */
MatrixVariable * free_MatrixVariable(MatrixVariable * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a MatrixVariable obj. Should be trappable");    
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
    if( obj->name != NULL)   
      ckfree(obj->name);     
    if( obj->source != NULL) 
      ckfree(obj->source);   
    /* obj->data is linked in */ 


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
