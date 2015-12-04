#ifdef _cplusplus
extern "C" {
#endif
#include "anonobj.h"


/* Function:  free_AnonymousObject(ao)
 *
 * Descrip:    frees an anonymous object
 *
 *
 * Arg:        ao [UNKN ] Undocumented argument [AnonymousObject *]
 *
 * Return [UNKN ]  Undocumented return value [AnonymousObject *]
 *
 */
# line 26 "anonobj.dy"
AnonymousObject * free_AnonymousObject(AnonymousObject * ao)
{
  int return_early = 0;

#ifdef PTHREAD   
    assert(pthread_mutex_lock(&(ao->dynamite_mutex)) == 0); 
#endif   
    if( ao->dynamite_hard_link > 1)     {  
      return_early = 1;  
      ao->dynamite_hard_link--; 
      }  
#ifdef PTHREAD   
    assert(pthread_mutex_unlock(&(ao->dynamite_mutex)) == 0);   
#endif   
    if( return_early == 1)   
      return NULL;   

  if( ao->obj != NULL ) {
    if( ao->free_object == NULL ) {
      warn("With object of %s, no registered free function, but with a present object. Likely to be leaking memory now",ao->type);
    } else {
      (*ao->free_object)(ao->obj);
    }
  }

  free(ao);

  return NULL;

}


# line 47 "anonobj.c"
/* Function:  hard_link_AnonymousObject(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AnonymousObject *]
 *
 * Return [UNKN ]  Undocumented return value [AnonymousObject *]
 *
 */
AnonymousObject * hard_link_AnonymousObject(AnonymousObject * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a AnonymousObject object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  AnonymousObject_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AnonymousObject *]
 *
 */
AnonymousObject * AnonymousObject_alloc(void) 
{
    AnonymousObject * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(AnonymousObject *) ckalloc (sizeof(AnonymousObject))) == NULL)  {  
      warn("AnonymousObject_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->type = NULL;    
    out->free_object = NULL; 


    return out;  
}    


/* Function:  swap_AnonymousObjectList(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_AnonymousObjectList
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [AnonymousObject **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_AnonymousObjectList(AnonymousObject ** list,int i,int j)  
{
    AnonymousObject * temp;  
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_AnonymousObjectList(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_AnonymousObjectList which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [AnonymousObject **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_AnonymousObjectList(AnonymousObject ** list,int left,int right,int (*comp)(AnonymousObject * ,AnonymousObject * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_AnonymousObjectList(list,left,(left+right)/2);  
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_AnonymousObjectList (list,++last,i);    
      }  
    swap_AnonymousObjectList (list,left,last);   
    qsort_AnonymousObjectList(list,left,last-1,comp);    
    qsort_AnonymousObjectList(list,last+1,right,comp);   
}    


/* Function:  sort_AnonymousObjectList(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_AnonymousObjectList
 *
 *
 * Arg:         obj [UNKN ] Object containing list [AnonymousObjectList *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_AnonymousObjectList(AnonymousObjectList * obj,int (*comp)(AnonymousObject *, AnonymousObject *)) 
{
    qsort_AnonymousObjectList(obj->anon,0,obj->len-1,comp);  
    return;  
}    


/* Function:  expand_AnonymousObjectList(obj,len)
 *
 * Descrip:    Really an internal function for add_AnonymousObjectList
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [AnonymousObjectList *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_AnonymousObjectList(AnonymousObjectList * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_AnonymousObjectList called with no need");    
      return TRUE;   
      }  


    if( (obj->anon = (AnonymousObject ** ) ckrealloc (obj->anon,sizeof(AnonymousObject *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_AnonymousObjectList, returning FALSE");  
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_AnonymousObjectList(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [AnonymousObjectList *]
 * Arg:        add [OWNER] Object to add to the list [AnonymousObject *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_AnonymousObjectList(AnonymousObjectList * obj,AnonymousObject * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_AnonymousObjectList(obj,obj->len + AnonymousObjectListLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->anon[obj->len++]=add;   
    return TRUE; 
}    


/* Function:  flush_AnonymousObjectList(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [AnonymousObjectList *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_AnonymousObjectList(AnonymousObjectList * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->anon[i] != NULL)  {  
        free_AnonymousObject(obj->anon[i]);  
        obj->anon[i] = NULL; 
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  AnonymousObjectList_alloc_std(void)
 *
 * Descrip:    Equivalent to AnonymousObjectList_alloc_len(AnonymousObjectListLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AnonymousObjectList *]
 *
 */
AnonymousObjectList * AnonymousObjectList_alloc_std(void) 
{
    return AnonymousObjectList_alloc_len(AnonymousObjectListLISTLENGTH); 
}    


/* Function:  AnonymousObjectList_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [AnonymousObjectList *]
 *
 */
AnonymousObjectList * AnonymousObjectList_alloc_len(int len) 
{
    AnonymousObjectList * out;  /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = AnonymousObjectList_alloc()) == NULL)  
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->anon = (AnonymousObject ** ) ckcalloc (len,sizeof(AnonymousObject *))) == NULL) {  
      warn("Warning, ckcalloc failed in AnonymousObjectList_alloc_len"); 
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_AnonymousObjectList(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AnonymousObjectList *]
 *
 * Return [UNKN ]  Undocumented return value [AnonymousObjectList *]
 *
 */
AnonymousObjectList * hard_link_AnonymousObjectList(AnonymousObjectList * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a AnonymousObjectList object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  AnonymousObjectList_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AnonymousObjectList *]
 *
 */
AnonymousObjectList * AnonymousObjectList_alloc(void) 
{
    AnonymousObjectList * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(AnonymousObjectList *) ckalloc (sizeof(AnonymousObjectList))) == NULL)  {  
      warn("AnonymousObjectList_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->anon = NULL;    
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_AnonymousObjectList(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AnonymousObjectList *]
 *
 * Return [UNKN ]  Undocumented return value [AnonymousObjectList *]
 *
 */
AnonymousObjectList * free_AnonymousObjectList(AnonymousObjectList * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a AnonymousObjectList obj. Should be trappable");   
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
    if( obj->anon != NULL)   {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->anon[i] != NULL)    
          free_AnonymousObject(obj->anon[i]);    
        }  
      ckfree(obj->anon); 
      }  


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
