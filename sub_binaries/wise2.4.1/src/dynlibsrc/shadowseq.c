#ifdef _cplusplus
extern "C" {
#endif
#include "shadowseq.h"



/* Function:  dump_ShadowSequence(shadow,ofp)
 *
 * Descrip:    Dumps shadow sequence out to file
 *
 *
 * Arg:        shadow [UNKN ] Undocumented argument [ShadowSequence *]
 * Arg:           ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 31 "shadowseq.dy"
void dump_ShadowSequence(ShadowSequence * shadow,FILE * ofp)
{
  int i;
  assert(shadow);
  assert(ofp);

  fprintf(ofp,"Sequence %s %d\n",shadow->seq->name,shadow->seq->len);
  for(i=0;i<shadow->len;i++) {
    fprintf(ofp,"  Shadow %s %d,%d %d\n",
	    shadow->region[i]->seq->name,
	    shadow->region[i]->start_shadow,
	    shadow->region[i]->start_seq,
	    shadow->region[i]->len);
  }
  
}


/* Function:  add_if_possible_ShadowSequence(shadow,seq,min_extension,start_shadow,start_seq,shadow_rate)
 *
 * Descrip:    Looks as to whether there is a good extension to be made, if
 *             so, does it and adds it to the ShadowSeq. Returns 0 if extension
 *             fails, if succeeds returns end point
 *
 *
 * Arg:               shadow [UNKN ] Undocumented argument [ShadowSequence *]
 * Arg:                  seq [UNKN ] Undocumented argument [Sequence *]
 * Arg:        min_extension [UNKN ] Undocumented argument [int]
 * Arg:         start_shadow [UNKN ] Undocumented argument [int]
 * Arg:            start_seq [UNKN ] Undocumented argument [int]
 * Arg:          shadow_rate [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 54 "shadowseq.dy"
int add_if_possible_ShadowSequence(ShadowSequence * shadow,Sequence * seq,int min_extension,int start_shadow,int start_seq,int shadow_rate)
{
  int i;
  int j;
  int error = 0;

  assert(shadow);
  assert(seq);

  for(i=start_shadow,j=start_seq;i<shadow->seq->len && j<seq->len;i++,j++) {
    if( shadow->seq->seq[i] != seq->seq[j] ) {
      error++;
      if( i > start_shadow && error > (int)((shadow_rate*100)/(start_shadow-i)) ) 
	break;
    }
  }

  if( i - start_shadow < min_extension ) {
    return 0;
  }

  /* we have a region */

  add_region_ShadowSequence(shadow,seq,start_shadow,start_seq,i-start_shadow);

  return j;
}


/* Function:  add_region_ShadowSequence(shadow,seq,start_shadow,start_seq,len)
 *
 * Descrip:    Adds a region to a ShadowSequence, extending
 *             the array if necessary
 *
 *
 * Arg:              shadow [UNKN ] Undocumented argument [ShadowSequence *]
 * Arg:                 seq [UNKN ] Undocumented argument [Sequence *]
 * Arg:        start_shadow [UNKN ] Undocumented argument [int]
 * Arg:           start_seq [UNKN ] Undocumented argument [int]
 * Arg:                 len [UNKN ] Undocumented argument [int]
 *
 */
# line 87 "shadowseq.dy"
void add_region_ShadowSequence(ShadowSequence * shadow,Sequence * seq,int start_shadow,int start_seq,int len)
{
  ShadowSeqRegion * ssr;


  assert(shadow);
  assert(seq);

  ssr = ShadowSeqRegion_alloc();
  ssr->seq = seq;
  ssr->start_shadow = start_shadow;
  ssr->start_seq = start_seq;
  ssr->len = len;


  if( shadow->region == NULL ) {
    if((shadow->region = (ShadowSeqRegion ** ) ckcalloc (ShadowSequenceLISTLENGTH,sizeof(ShadowSeqRegion *))) == NULL)   {  
      warn("Warning, ckcalloc failed in ShadowSequence_alloc_len");  
      return;   
    }  
    shadow->len = 0;    
    shadow->maxlen = ShadowSequenceLISTLENGTH;   
  }

  add_ShadowSequence(shadow,ssr);

  return;
}



/* Function:  new_ShadowSequence(seq)
 *
 * Descrip:    Builds a new ShadowSequence from a Sequence
 *             doesnt hard link as memory should be handled
 *             outside of the shadowsystem
 *
 *
 * Arg:        seq [UNKN ] Undocumented argument [Sequence *]
 *
 * Return [UNKN ]  Undocumented return value [ShadowSequence *]
 *
 */
# line 123 "shadowseq.dy"
ShadowSequence * new_ShadowSequence(Sequence * seq)
{
  ShadowSequence * out;
  
  assert(seq);

  out = ShadowSequence_alloc();
  assert(out);
  out->seq = seq;
  out->region = NULL;
  out->dirty = 0;
  return out;
}



# line 148 "shadowseq.c"
/* Function:  hard_link_ShadowSeqRegion(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ShadowSeqRegion *]
 *
 * Return [UNKN ]  Undocumented return value [ShadowSeqRegion *]
 *
 */
ShadowSeqRegion * hard_link_ShadowSeqRegion(ShadowSeqRegion * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a ShadowSeqRegion object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  ShadowSeqRegion_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ShadowSeqRegion *]
 *
 */
ShadowSeqRegion * ShadowSeqRegion_alloc(void) 
{
    ShadowSeqRegion * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(ShadowSeqRegion *) ckalloc (sizeof(ShadowSeqRegion))) == NULL)  {  
      warn("ShadowSeqRegion_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->start_shadow = 0;   
    out->start_seq = 0;  
    out->len = 0;    
    out->seq = NULL; 


    return out;  
}    


/* Function:  free_ShadowSeqRegion(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ShadowSeqRegion *]
 *
 * Return [UNKN ]  Undocumented return value [ShadowSeqRegion *]
 *
 */
ShadowSeqRegion * free_ShadowSeqRegion(ShadowSeqRegion * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a ShadowSeqRegion obj. Should be trappable");   
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


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_ShadowSequence(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_ShadowSequence
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [ShadowSeqRegion **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_ShadowSequence(ShadowSeqRegion ** list,int i,int j)  
{
    ShadowSeqRegion * temp;  
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_ShadowSequence(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_ShadowSequence which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [ShadowSeqRegion **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_ShadowSequence(ShadowSeqRegion ** list,int left,int right,int (*comp)(ShadowSeqRegion * ,ShadowSeqRegion * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_ShadowSequence(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_ShadowSequence (list,++last,i); 
      }  
    swap_ShadowSequence (list,left,last);    
    qsort_ShadowSequence(list,left,last-1,comp); 
    qsort_ShadowSequence(list,last+1,right,comp);    
}    


/* Function:  sort_ShadowSequence(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_ShadowSequence
 *
 *
 * Arg:         obj [UNKN ] Object containing list [ShadowSequence *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_ShadowSequence(ShadowSequence * obj,int (*comp)(ShadowSeqRegion *, ShadowSeqRegion *)) 
{
    qsort_ShadowSequence(obj->region,0,obj->len-1,comp); 
    return;  
}    


/* Function:  expand_ShadowSequence(obj,len)
 *
 * Descrip:    Really an internal function for add_ShadowSequence
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [ShadowSequence *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_ShadowSequence(ShadowSequence * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_ShadowSequence called with no need"); 
      return TRUE;   
      }  


    if( (obj->region = (ShadowSeqRegion ** ) ckrealloc (obj->region,sizeof(ShadowSeqRegion *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_ShadowSequence, returning FALSE");   
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_ShadowSequence(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [ShadowSequence *]
 * Arg:        add [OWNER] Object to add to the list [ShadowSeqRegion *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_ShadowSequence(ShadowSequence * obj,ShadowSeqRegion * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_ShadowSequence(obj,obj->len + ShadowSequenceLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->region[obj->len++]=add; 
    return TRUE; 
}    


/* Function:  flush_ShadowSequence(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [ShadowSequence *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_ShadowSequence(ShadowSequence * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->region[i] != NULL)    {  
        free_ShadowSeqRegion(obj->region[i]);    
        obj->region[i] = NULL;   
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  ShadowSequence_alloc_std(void)
 *
 * Descrip:    Equivalent to ShadowSequence_alloc_len(ShadowSequenceLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ShadowSequence *]
 *
 */
ShadowSequence * ShadowSequence_alloc_std(void) 
{
    return ShadowSequence_alloc_len(ShadowSequenceLISTLENGTH);   
}    


/* Function:  ShadowSequence_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [ShadowSequence *]
 *
 */
ShadowSequence * ShadowSequence_alloc_len(int len) 
{
    ShadowSequence * out;   /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = ShadowSequence_alloc()) == NULL)   
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->region = (ShadowSeqRegion ** ) ckcalloc (len,sizeof(ShadowSeqRegion *))) == NULL)   {  
      warn("Warning, ckcalloc failed in ShadowSequence_alloc_len");  
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_ShadowSequence(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ShadowSequence *]
 *
 * Return [UNKN ]  Undocumented return value [ShadowSequence *]
 *
 */
ShadowSequence * hard_link_ShadowSequence(ShadowSequence * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a ShadowSequence object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  ShadowSequence_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ShadowSequence *]
 *
 */
ShadowSequence * ShadowSequence_alloc(void) 
{
    ShadowSequence * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(ShadowSequence *) ckalloc (sizeof(ShadowSequence))) == NULL)    {  
      warn("ShadowSequence_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->dirty = 'u';    
    out->region = NULL;  
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_ShadowSequence(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ShadowSequence *]
 *
 * Return [UNKN ]  Undocumented return value [ShadowSequence *]
 *
 */
ShadowSequence * free_ShadowSequence(ShadowSequence * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a ShadowSequence obj. Should be trappable");    
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
    /* obj->seq is linked in */ 
    if( obj->region != NULL) {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->region[i] != NULL)  
          free_ShadowSeqRegion(obj->region[i]);  
        }  
      ckfree(obj->region);   
      }  


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
