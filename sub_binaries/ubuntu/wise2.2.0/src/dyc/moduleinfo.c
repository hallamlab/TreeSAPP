#ifdef _cplusplus
extern "C" {
#endif
#include "moduleinfo.h"


# line 6 "moduleinfo.c"
/* Function:  hard_link_ModuleInfo(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ModuleInfo *]
 *
 * Return [UNKN ]  Undocumented return value [ModuleInfo *]
 *
 */
ModuleInfo * hard_link_ModuleInfo(ModuleInfo * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a ModuleInfo object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  ModuleInfo_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ModuleInfo *]
 *
 */
ModuleInfo * ModuleInfo_alloc(void) 
{
    ModuleInfo * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(ModuleInfo *) ckalloc (sizeof(ModuleInfo))) == NULL)    {  
      warn("ModuleInfo_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
    out->ft = NULL;  


    return out;  
}    


/* Function:  free_ModuleInfo(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ModuleInfo *]
 *
 * Return [UNKN ]  Undocumented return value [ModuleInfo *]
 *
 */
ModuleInfo * free_ModuleInfo(ModuleInfo * obj) 
{


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a ModuleInfo obj. Should be trappable");    
      return NULL;   
      }  


    if( obj->dynamite_hard_link > 1)     {  
      obj->dynamite_hard_link--; 
      return NULL;   
      }  
    if( obj->ft != NULL) 
      free_Ftext(obj->ft);   


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
