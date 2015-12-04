#ifdef _cplusplus
extern "C" {
#endif
#include "kmer_count.h"

# line 24 "kmer_count.dy"
KmerCount * new_KmerCount_KmerCountAllocator(KmerCountAllocator * kca)
{

  if( kca->current_count >= KMER_COUNT_BLOCKSIZE ) {
    add_KmerCountAllocator(kca,KmerCountBlock_alloc());
    kca->current_count = 0;
  }

  return &(kca->block[kca->len-1]->count[kca->current_count++]);
}


# line 36 "kmer_count.dy"
KmerCountAllocator * new_KmerCountAllocator(void)
{
  KmerCountAllocator * out;

  out = KmerCountAllocator_alloc_std();
  add_KmerCountAllocator(out,KmerCountBlock_alloc());
  out->current_count = 0;

  return out;
}
  


# line 32 "kmer_count.c"
/* Function:  hard_link_KmerCount(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [KmerCount *]
 *
 * Return [UNKN ]  Undocumented return value [KmerCount *]
 *
 */
KmerCount * hard_link_KmerCount(KmerCount * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a KmerCount object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  KmerCount_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [KmerCount *]
 *
 */
KmerCount * KmerCount_alloc(void) 
{
    KmerCount * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(KmerCount *) ckalloc (sizeof(KmerCount))) == NULL)  {  
      warn("KmerCount_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->count = 0;  


    return out;  
}    


/* Function:  free_KmerCount(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [KmerCount *]
 *
 * Return [UNKN ]  Undocumented return value [KmerCount *]
 *
 */
KmerCount * free_KmerCount(KmerCount * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a KmerCount obj. Should be trappable"); 
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


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_KmerCountBlock(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [KmerCountBlock *]
 *
 * Return [UNKN ]  Undocumented return value [KmerCountBlock *]
 *
 */
KmerCountBlock * hard_link_KmerCountBlock(KmerCountBlock * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a KmerCountBlock object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  KmerCountBlock_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [KmerCountBlock *]
 *
 */
KmerCountBlock * KmerCountBlock_alloc(void) 
{
    KmerCountBlock * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(KmerCountBlock *) ckalloc (sizeof(KmerCountBlock))) == NULL)    {  
      warn("KmerCountBlock_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    /* count[KMER_COUNT_BLOCKSIZE] is an array: no default possible */ 


    return out;  
}    


/* Function:  free_KmerCountBlock(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [KmerCountBlock *]
 *
 * Return [UNKN ]  Undocumented return value [KmerCountBlock *]
 *
 */
KmerCountBlock * free_KmerCountBlock(KmerCountBlock * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a KmerCountBlock obj. Should be trappable");    
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


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_KmerCountAllocator(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_KmerCountAllocator
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [KmerCountBlock **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_KmerCountAllocator(KmerCountBlock ** list,int i,int j)  
{
    KmerCountBlock * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_KmerCountAllocator(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_KmerCountAllocator which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [KmerCountBlock **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_KmerCountAllocator(KmerCountBlock ** list,int left,int right,int (*comp)(KmerCountBlock * ,KmerCountBlock * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_KmerCountAllocator(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_KmerCountAllocator (list,++last,i); 
      }  
    swap_KmerCountAllocator (list,left,last);    
    qsort_KmerCountAllocator(list,left,last-1,comp); 
    qsort_KmerCountAllocator(list,last+1,right,comp);    
}    


/* Function:  sort_KmerCountAllocator(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_KmerCountAllocator
 *
 *
 * Arg:         obj [UNKN ] Object containing list [KmerCountAllocator *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_KmerCountAllocator(KmerCountAllocator * obj,int (*comp)(KmerCountBlock *, KmerCountBlock *)) 
{
    qsort_KmerCountAllocator(obj->block,0,obj->len-1,comp);  
    return;  
}    


/* Function:  expand_KmerCountAllocator(obj,len)
 *
 * Descrip:    Really an internal function for add_KmerCountAllocator
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [KmerCountAllocator *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_KmerCountAllocator(KmerCountAllocator * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_KmerCountAllocator called with no need"); 
      return TRUE;   
      }  


    if( (obj->block = (KmerCountBlock ** ) ckrealloc (obj->block,sizeof(KmerCountBlock *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_KmerCountAllocator, returning FALSE");   
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_KmerCountAllocator(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [KmerCountAllocator *]
 * Arg:        add [OWNER] Object to add to the list [KmerCountBlock *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_KmerCountAllocator(KmerCountAllocator * obj,KmerCountBlock * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_KmerCountAllocator(obj,obj->len + KmerCountAllocatorLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->block[obj->len++]=add;  
    return TRUE; 
}    


/* Function:  flush_KmerCountAllocator(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [KmerCountAllocator *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_KmerCountAllocator(KmerCountAllocator * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->block[i] != NULL) {  
        free_KmerCountBlock(obj->block[i]);  
        obj->block[i] = NULL;    
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  KmerCountAllocator_alloc_std(void)
 *
 * Descrip:    Equivalent to KmerCountAllocator_alloc_len(KmerCountAllocatorLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [KmerCountAllocator *]
 *
 */
KmerCountAllocator * KmerCountAllocator_alloc_std(void) 
{
    return KmerCountAllocator_alloc_len(KmerCountAllocatorLISTLENGTH);   
}    


/* Function:  KmerCountAllocator_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [KmerCountAllocator *]
 *
 */
KmerCountAllocator * KmerCountAllocator_alloc_len(int len) 
{
    KmerCountAllocator * out;   /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = KmerCountAllocator_alloc()) == NULL)   
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->block = (KmerCountBlock ** ) ckcalloc (len,sizeof(KmerCountBlock *))) == NULL)  {  
      warn("Warning, ckcalloc failed in KmerCountAllocator_alloc_len");  
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_KmerCountAllocator(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [KmerCountAllocator *]
 *
 * Return [UNKN ]  Undocumented return value [KmerCountAllocator *]
 *
 */
KmerCountAllocator * hard_link_KmerCountAllocator(KmerCountAllocator * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a KmerCountAllocator object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  KmerCountAllocator_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [KmerCountAllocator *]
 *
 */
KmerCountAllocator * KmerCountAllocator_alloc(void) 
{
    KmerCountAllocator * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(KmerCountAllocator *) ckalloc (sizeof(KmerCountAllocator))) == NULL)    {  
      warn("KmerCountAllocator_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->block = NULL;   
    out->len = out->maxlen = 0;  
    out->current_count = 0;  


    return out;  
}    


/* Function:  free_KmerCountAllocator(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [KmerCountAllocator *]
 *
 * Return [UNKN ]  Undocumented return value [KmerCountAllocator *]
 *
 */
KmerCountAllocator * free_KmerCountAllocator(KmerCountAllocator * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a KmerCountAllocator obj. Should be trappable");    
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
    if( obj->block != NULL)  {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->block[i] != NULL)   
          free_KmerCountBlock(obj->block[i]);    
        }  
      ckfree(obj->block);    
      }  


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
