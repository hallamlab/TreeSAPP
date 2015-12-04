#ifdef _cplusplus
extern "C" {
#endif
#include "linkmap.h"

# line 18 "linkmap.dy"
char * update_links_LinkMap(char * line,char * startid,char * stopid,char * nolinkstr,LinkSet * ls)
{
  char buffer[MAXLINE];
  char * run;
  char c;
  char * end;
  LinkMap * lm;

  for(end=line,run=line;(run = strstr(run,startid)) != NULL;) {

    /*** strcat the previous string ***/

    c = *run;
    *run =  '\0';
    strcat(buffer,end);

    /*** now put in the link system ***/

    end = strstr(run,stopid);
    c = *end;
    *end = '\0';
    if( (lm = LinkMap_from_word(run,ls)) == NULL ) {
      
    }

  }

}

# line 47 "linkmap.dy"
LinkMap * new_std_LinkMap(char * word,char * urlstub,char * textstub)
{
  char buffer[128];
  char * run;
  LinkMap * out;

  out= LinkMap_alloc();

  push_scan_and_replace_pair("$LINK",word);

  strcpy(buffer,urlstub);
  run = scan_and_replace_line(buffer);
  out->urltail = stringalloc(run);
  strcpy(buffer,textstub);
  run = scan_and_replace_line(buffer);
  out->text = stringalloc(textstub);

  pop_scan_and_replace_pair();

  return out;
}
  

# line 70 "linkmap.dy"
LinkMap * LinkMap_from_word(char * word,LinkSet * ls)
{
  int i;

  for(i=0;i<ls->len;i++) {
    if( strcmp(ls->lm[i]->word,word) == 0 )
      return lm[i];
  }

  return NULL;
}





# line 76 "linkmap.c"
/* Function:  hard_link_LinkMap(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [LinkMap *]
 *
 * Return [UNKN ]  Undocumented return value [LinkMap *]
 *
 */
LinkMap * hard_link_LinkMap(LinkMap * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a LinkMap object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  LinkMap_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LinkMap *]
 *
 */
LinkMap * LinkMap_alloc(void) 
{
    LinkMap * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(LinkMap *) ckalloc (sizeof(LinkMap))) == NULL)  {  
      warn("LinkMap_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->word = NULL;    
    out->urltail = NULL; 
    out->text = NULL;    


    return out;  
}    


/* Function:  free_LinkMap(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [LinkMap *]
 *
 * Return [UNKN ]  Undocumented return value [LinkMap *]
 *
 */
LinkMap * free_LinkMap(LinkMap * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a LinkMap obj. Should be trappable");   
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
    if( obj->word != NULL)   
      ckfree(obj->word);     
    if( obj->urltail != NULL)    
      ckfree(obj->urltail);  
    if( obj->text != NULL)   
      ckfree(obj->text);     


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_LinkSet(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_LinkSet
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [LinkMap **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_LinkSet(LinkMap ** list,int i,int j)  
{
    LinkMap * temp;  
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_LinkSet(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_LinkSet which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [LinkMap **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_LinkSet(LinkMap ** list,int left,int right,int (*comp)(LinkMap * ,LinkMap * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_LinkSet(list,left,(left+right)/2);  
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_LinkSet (list,++last,i);    
      }  
    swap_LinkSet (list,left,last);   
    qsort_LinkSet(list,left,last-1,comp);    
    qsort_LinkSet(list,last+1,right,comp);   
}    


/* Function:  sort_LinkSet(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_LinkSet
 *
 *
 * Arg:         obj [UNKN ] Object containing list [LinkSet *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_LinkSet(LinkSet * obj,int (*comp)(LinkMap *, LinkMap *)) 
{
    qsort_LinkSet(obj->lm,0,obj->len-1,comp);    
    return;  
}    


/* Function:  expand_LinkSet(obj,len)
 *
 * Descrip:    Really an internal function for add_LinkSet
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [LinkSet *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_LinkSet(LinkSet * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_LinkSet called with no need");    
      return TRUE;   
      }  


    if( (obj->lm = (LinkMap ** ) ckrealloc (obj->lm,sizeof(LinkMap *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_LinkSet, returning FALSE");  
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_LinkSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [LinkSet *]
 * Arg:        add [OWNER] Object to add to the list [LinkMap *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_LinkSet(LinkSet * obj,LinkMap * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_LinkSet(obj,obj->len + LinkSetLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->lm[obj->len++]=add; 
    return TRUE; 
}    


/* Function:  flush_LinkSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [LinkSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_LinkSet(LinkSet * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->lm[i] != NULL)    {  
        free_LinkMap(obj->lm[i]);    
        obj->lm[i] = NULL;   
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  LinkSet_alloc_std(void)
 *
 * Descrip:    Equivalent to LinkSet_alloc_len(LinkSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LinkSet *]
 *
 */
LinkSet * LinkSet_alloc_std(void) 
{
    return LinkSet_alloc_len(LinkSetLISTLENGTH); 
}    


/* Function:  LinkSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [LinkSet *]
 *
 */
LinkSet * LinkSet_alloc_len(int len) 
{
    LinkSet * out;  /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = LinkSet_alloc()) == NULL)  
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->lm = (LinkMap ** ) ckcalloc (len,sizeof(LinkMap *))) == NULL)   {  
      warn("Warning, ckcalloc failed in LinkSet_alloc_len"); 
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_LinkSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [LinkSet *]
 *
 * Return [UNKN ]  Undocumented return value [LinkSet *]
 *
 */
LinkSet * hard_link_LinkSet(LinkSet * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a LinkSet object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  LinkSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LinkSet *]
 *
 */
LinkSet * LinkSet_alloc(void) 
{
    LinkSet * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(LinkSet *) ckalloc (sizeof(LinkSet))) == NULL)  {  
      warn("LinkSet_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->lm = NULL;  
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_LinkSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [LinkSet *]
 *
 * Return [UNKN ]  Undocumented return value [LinkSet *]
 *
 */
LinkSet * free_LinkSet(LinkSet * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a LinkSet obj. Should be trappable");   
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
    if( obj->lm != NULL) {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->lm[i] != NULL)  
          free_LinkMap(obj->lm[i]);  
        }  
      ckfree(obj->lm);   
      }  


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
