#ifdef _cplusplus
extern "C" {
#endif
#include "objectinfo.h"

# line 16 "objectinfo.dy"
int write_C_ObjectInfo(ObjectInfo * oi,FILE * ofp)
{
  return show_eddystyle_Ftext(oi->ft,"Descrip:",10,ofp,"No Description");
}
			     
# line 21 "objectinfo.dy"
ObjectInfo * read_ObjectInfo_line_func(char * line,int maxline,FILE * ifp,char * (fgets_func)(char *,int,FILE*))
{
  ObjectInfo * out;

  out = ObjectInfo_alloc();

  out->ft = read_Ftext(line,maxline,ifp,"%%",fgets_func);

  return out;
}


# line 24 "objectinfo.c"
/* Function:  hard_link_ObjectInfo(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ObjectInfo *]
 *
 * Return [UNKN ]  Undocumented return value [ObjectInfo *]
 *
 */
ObjectInfo * hard_link_ObjectInfo(ObjectInfo * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a ObjectInfo object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  ObjectInfo_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ObjectInfo *]
 *
 */
ObjectInfo * ObjectInfo_alloc(void) 
{
    ObjectInfo * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(ObjectInfo *) ckalloc (sizeof(ObjectInfo))) == NULL)    {  
      warn("ObjectInfo_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->ft = NULL;  


    return out;  
}    


/* Function:  free_ObjectInfo(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ObjectInfo *]
 *
 * Return [UNKN ]  Undocumented return value [ObjectInfo *]
 *
 */
ObjectInfo * free_ObjectInfo(ObjectInfo * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a ObjectInfo obj. Should be trappable");    
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
    if( obj->ft != NULL) 
      free_Ftext(obj->ft);   


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
