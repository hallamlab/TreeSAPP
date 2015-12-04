#ifdef _cplusplus
extern "C" {
#endif
#include "subseqlookup.h"


/* Function:  new_array_SeqLookupInterface(array_max_size)
 *
 * Descrip:    Makes a new array Lookup system
 *
 *
 * Arg:        array_max_size [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupInterface *]
 *
 */
# line 32 "subseqlookup.dy"
SeqLookupInterface * new_array_SeqLookupInterface(int array_max_size)
{
  SeqLookupInterface * out;
  SubSeqLookup * a;
  int i;

  out = SeqLookupInterface_alloc_std();
  a = SubSeqLookup_alloc();
  a->array = calloc(array_max_size,sizeof(SeqLookupPos *));
  for(i=0;i<array_max_size;i++) {
    a->array[i] = NULL;
  }

  a->array_max = array_max_size;
  a->block = new_SeqLookupPosBlockAllocator();

  out->data = (void*) a;
  out->lookup = lookup_subseqlookup;
  out->add = add_subseqlookup;
  out->is_populated = is_populated_subseqlookup;
  out->free_data = free_subseqlookup;
  
  return out;
}


/* Function:  free_subseqlookup(data)
 *
 * Descrip:    free function for the hash
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 *
 */
# line 61 "subseqlookup.dy"
void free_subseqlookup(void * data)
{
  SubSeqLookup * look = (SubSeqLookup *)data;
 
  free_SubSeqLookup(look);
}

/* Function:  is_populated_subseqlookup(data,seq_number)
 *
 * Descrip:    tells whether this is populated or not
 *
 *
 * Arg:              data [UNKN ] Undocumented argument [void *]
 * Arg:        seq_number [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 71 "subseqlookup.dy"
boolean is_populated_subseqlookup(void * data, int seq_number)
{
  SubSeqLookup * look = (SubSeqLookup *)data;

  if( look->array[seq_number] == NULL ) {
    return FALSE;
  } else {
    return TRUE;
  }

}

/* Function:  lookup_subseqlookup(data,seq_number)
 *
 * Descrip:    Retrieves a SeqLookup position 
 *
 *
 * Arg:              data [UNKN ] Undocumented argument [void *]
 * Arg:        seq_number [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupResultInterface *]
 *
 */
# line 86 "subseqlookup.dy"
SeqLookupResultInterface * lookup_subseqlookup(void * data, int seq_number)
{
  SubSeqLookup * look = (SubSeqLookup *)data;

  return new_linkedl_SeqLookupResultInterface(look->array[seq_number]);
}

/* Function:  add_subseqlookup(data,seq_number,seq,pos)
 *
 * Descrip:    Adds a sequence/pos pair to the hash
 *
 *
 * Arg:              data [UNKN ] Undocumented argument [void *]
 * Arg:        seq_number [UNKN ] Undocumented argument [int]
 * Arg:               seq [UNKN ] Undocumented argument [Sequence *]
 * Arg:               pos [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 96 "subseqlookup.dy"
boolean add_subseqlookup(void * data,int seq_number,Sequence * seq,int pos)
{
  SubSeqLookup * look = (SubSeqLookup *)data;

  SeqLookupPos * p;

  p = new_SeqLookupPos_BlockAllocator(look->block);
  p->seq = seq;
  p->pos = pos;
  p->next = NULL;

  if( look->array[seq_number] == NULL ) {
    look->array[seq_number] = p;
  } else {
    p->next = look->array[seq_number];
    look->array[seq_number] = p;
  }

  return TRUE;
}

/* Function:  new_SeqLookupPosBlockAllocator(void)
 *
 * Descrip:    Makes a new SeqPosLookup Block Allocator
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupPosBlockAllocator *]
 *
 */
# line 120 "subseqlookup.dy"
SeqLookupPosBlockAllocator * new_SeqLookupPosBlockAllocator(void)
{
  SeqLookupPosBlockAllocator * out;
  
  out = SeqLookupPosBlockAllocator_alloc_std();

  add_SeqLookupPosBlockAllocator(out,SeqLookupPosBlock_alloc());
  out->pos = 0;
 
  return out;
}
/* Function:  new_SeqLookupPos_BlockAllocator(*bla)
 *
 * Descrip:    Returns a new SeqPosLookup
 *
 *
 * Arg:        *bla [UNKN ] Undocumented argument [SeqLookupPosBlockAllocator]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupPos *]
 *
 */
# line 134 "subseqlookup.dy"
SeqLookupPos * new_SeqLookupPos_BlockAllocator(SeqLookupPosBlockAllocator *bla)
{
  if( bla->pos+1 > LOOKUP_BLOCK_SIZE ) {
    add_SeqLookupPosBlockAllocator(bla,SeqLookupPosBlock_alloc());
    bla->pos = 0;
    fprintf(stderr,"Increasing block...\n");
  }

  bla->block[bla->len-1]->block[bla->pos].dynamite_hard_link = 1;

  bla->pos++;
  return &(bla->block[bla->len-1]->block[bla->pos-1]);
}
# line 168 "subseqlookup.c"
/* Function:  hard_link_SeqLookupPosBlock(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SeqLookupPosBlock *]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupPosBlock *]
 *
 */
SeqLookupPosBlock * hard_link_SeqLookupPosBlock(SeqLookupPosBlock * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a SeqLookupPosBlock object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  SeqLookupPosBlock_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupPosBlock *]
 *
 */
SeqLookupPosBlock * SeqLookupPosBlock_alloc(void) 
{
    SeqLookupPosBlock * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(SeqLookupPosBlock *) ckalloc (sizeof(SeqLookupPosBlock))) == NULL)  {  
      warn("SeqLookupPosBlock_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    /* block[LOOKUP_BLOCK_SIZE] is an array: no default possible */ 


    return out;  
}    


/* Function:  free_SeqLookupPosBlock(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SeqLookupPosBlock *]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupPosBlock *]
 *
 */
SeqLookupPosBlock * free_SeqLookupPosBlock(SeqLookupPosBlock * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a SeqLookupPosBlock obj. Should be trappable"); 
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


/* Function:  swap_SeqLookupPosBlockAllocator(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_SeqLookupPosBlockAllocator
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [SeqLookupPosBlock **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_SeqLookupPosBlockAllocator(SeqLookupPosBlock ** list,int i,int j)  
{
    SeqLookupPosBlock * temp;    
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_SeqLookupPosBlockAllocator(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_SeqLookupPosBlockAllocator which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [SeqLookupPosBlock **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_SeqLookupPosBlockAllocator(SeqLookupPosBlock ** list,int left,int right,int (*comp)(SeqLookupPosBlock * ,SeqLookupPosBlock * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_SeqLookupPosBlockAllocator(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_SeqLookupPosBlockAllocator (list,++last,i); 
      }  
    swap_SeqLookupPosBlockAllocator (list,left,last);    
    qsort_SeqLookupPosBlockAllocator(list,left,last-1,comp); 
    qsort_SeqLookupPosBlockAllocator(list,last+1,right,comp);    
}    


/* Function:  sort_SeqLookupPosBlockAllocator(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_SeqLookupPosBlockAllocator
 *
 *
 * Arg:         obj [UNKN ] Object containing list [SeqLookupPosBlockAllocator *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_SeqLookupPosBlockAllocator(SeqLookupPosBlockAllocator * obj,int (*comp)(SeqLookupPosBlock *, SeqLookupPosBlock *)) 
{
    qsort_SeqLookupPosBlockAllocator(obj->block,0,obj->len-1,comp);  
    return;  
}    


/* Function:  expand_SeqLookupPosBlockAllocator(obj,len)
 *
 * Descrip:    Really an internal function for add_SeqLookupPosBlockAllocator
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [SeqLookupPosBlockAllocator *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_SeqLookupPosBlockAllocator(SeqLookupPosBlockAllocator * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_SeqLookupPosBlockAllocator called with no need"); 
      return TRUE;   
      }  


    if( (obj->block = (SeqLookupPosBlock ** ) ckrealloc (obj->block,sizeof(SeqLookupPosBlock *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_SeqLookupPosBlockAllocator, returning FALSE");   
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_SeqLookupPosBlockAllocator(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [SeqLookupPosBlockAllocator *]
 * Arg:        add [OWNER] Object to add to the list [SeqLookupPosBlock *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_SeqLookupPosBlockAllocator(SeqLookupPosBlockAllocator * obj,SeqLookupPosBlock * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_SeqLookupPosBlockAllocator(obj,obj->len + SeqLookupPosBlockAllocatorLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->block[obj->len++]=add;  
    return TRUE; 
}    


/* Function:  flush_SeqLookupPosBlockAllocator(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [SeqLookupPosBlockAllocator *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_SeqLookupPosBlockAllocator(SeqLookupPosBlockAllocator * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->block[i] != NULL) {  
        free_SeqLookupPosBlock(obj->block[i]);   
        obj->block[i] = NULL;    
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  SeqLookupPosBlockAllocator_alloc_std(void)
 *
 * Descrip:    Equivalent to SeqLookupPosBlockAllocator_alloc_len(SeqLookupPosBlockAllocatorLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupPosBlockAllocator *]
 *
 */
SeqLookupPosBlockAllocator * SeqLookupPosBlockAllocator_alloc_std(void) 
{
    return SeqLookupPosBlockAllocator_alloc_len(SeqLookupPosBlockAllocatorLISTLENGTH);   
}    


/* Function:  SeqLookupPosBlockAllocator_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupPosBlockAllocator *]
 *
 */
SeqLookupPosBlockAllocator * SeqLookupPosBlockAllocator_alloc_len(int len) 
{
    SeqLookupPosBlockAllocator * out;   /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = SeqLookupPosBlockAllocator_alloc()) == NULL)   
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->block = (SeqLookupPosBlock ** ) ckcalloc (len,sizeof(SeqLookupPosBlock *))) == NULL)    {  
      warn("Warning, ckcalloc failed in SeqLookupPosBlockAllocator_alloc_len");  
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_SeqLookupPosBlockAllocator(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SeqLookupPosBlockAllocator *]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupPosBlockAllocator *]
 *
 */
SeqLookupPosBlockAllocator * hard_link_SeqLookupPosBlockAllocator(SeqLookupPosBlockAllocator * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a SeqLookupPosBlockAllocator object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  SeqLookupPosBlockAllocator_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupPosBlockAllocator *]
 *
 */
SeqLookupPosBlockAllocator * SeqLookupPosBlockAllocator_alloc(void) 
{
    SeqLookupPosBlockAllocator * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(SeqLookupPosBlockAllocator *) ckalloc (sizeof(SeqLookupPosBlockAllocator))) == NULL)    {  
      warn("SeqLookupPosBlockAllocator_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->block = NULL;   
    out->len = out->maxlen = 0;  
    out->pos = 0;    


    return out;  
}    


/* Function:  free_SeqLookupPosBlockAllocator(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SeqLookupPosBlockAllocator *]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupPosBlockAllocator *]
 *
 */
SeqLookupPosBlockAllocator * free_SeqLookupPosBlockAllocator(SeqLookupPosBlockAllocator * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a SeqLookupPosBlockAllocator obj. Should be trappable");    
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
          free_SeqLookupPosBlock(obj->block[i]); 
        }  
      ckfree(obj->block);    
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_SubSeqLookup(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SubSeqLookup *]
 *
 * Return [UNKN ]  Undocumented return value [SubSeqLookup *]
 *
 */
SubSeqLookup * hard_link_SubSeqLookup(SubSeqLookup * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a SubSeqLookup object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  SubSeqLookup_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SubSeqLookup *]
 *
 */
SubSeqLookup * SubSeqLookup_alloc(void) 
{
    SubSeqLookup * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(SubSeqLookup *) ckalloc (sizeof(SubSeqLookup))) == NULL)    {  
      warn("SubSeqLookup_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->array = NULL;   
    out->array_max = 0;  
    out->block = NULL;   


    return out;  
}    


/* Function:  free_SubSeqLookup(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SubSeqLookup *]
 *
 * Return [UNKN ]  Undocumented return value [SubSeqLookup *]
 *
 */
SubSeqLookup * free_SubSeqLookup(SubSeqLookup * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a SubSeqLookup obj. Should be trappable");  
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
    if( obj->array != NULL)  
      ckfree(obj->array);    
    if( obj->block != NULL)  
      free_SeqLookupPosBlockAllocator(obj->block);   


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
