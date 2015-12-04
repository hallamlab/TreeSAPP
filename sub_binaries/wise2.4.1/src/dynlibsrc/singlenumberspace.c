#ifdef _cplusplus
extern "C" {
#endif
#include "singlenumberspace.h"


/* Function:  lookup_ShadowSequence_SingleNumberSpace(space,pos)
 *
 * Descrip:    New return using binary choping
 *
 *
 * Arg:        space [UNKN ] Undocumented argument [SingleNumberSpace *]
 * Arg:          pos [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [SingleNumberSequence *]
 *
 */
# line 34 "singlenumberspace.dy"
SingleNumberSequence * lookup_ShadowSequence_SingleNumberSpace(SingleNumberSpace * space,int pos)
{
  int index_pos;

  assert(pos < space->current_end);

  if( space->is_static != 0 ) {
    /* can make explicit calculation for position */

    index_pos = (int) (pos / space->average_len);
    if( space->sns[index_pos]->start == -1 ) {
      /* was a padding position */
      for(index_pos--;index_pos >= 0 && space->sns[index_pos]->start == -1;index_pos--) {
	fprintf(stderr,"Regressed %d with %d [pos %d]\n",index_pos,space->sns[index_pos]->start,pos);
      }
      assert(index_pos >= 0 );
    }
	
    if( space->sns[index_pos] == NULL ) {
        fatal("Going to return %d for %d with %d - max %d\n",index_pos,pos,space->average_len,space->current_end);
    }

    return space->sns[index_pos];
  }

  
  assert(pos < space->current_end);
  
  index_pos = find_position_SingleNumberSpace(space,0,space->len-1,pos);

  
  return space->sns[index_pos];
}


/* Function:  find_position_SingleNumberSpace(space,lower,higher,position)
 *
 * Descrip:    Recursive function for finding position
 *
 *
 * Arg:           space [UNKN ] Undocumented argument [SingleNumberSpace *]
 * Arg:           lower [UNKN ] Undocumented argument [int]
 * Arg:          higher [UNKN ] Undocumented argument [int]
 * Arg:        position [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 72 "singlenumberspace.dy"
int find_position_SingleNumberSpace(SingleNumberSpace * space,int lower,int higher,int position)
{
  int mid;
  /*
  fprintf(stderr,"To find position %d Entering with.... %d to %d\n",position,lower,higher);
  */

  if( lower+2 >= higher ) {
    if( position >= space->sns[lower]->start && position < space->sns[lower]->end ) 
      return lower;
    if( position >= space->sns[higher]->start && position < space->sns[higher]->end ) 
      return higher;
  }

  mid = (lower + (int)((higher-lower)/2));

  if( position >= space->sns[mid]->start && position < space->sns[mid]->end ) 
    return mid;
  

  if( position >= space->sns[mid]->end ) {
    /* top half */
    return find_position_SingleNumberSpace(space,mid,higher,position);
  } else {
    return find_position_SingleNumberSpace(space,lower,mid,position);
  }

}


/* Function:  add_ShadowSequence_SingleNumberSpace(space,seq)
 *
 * Descrip:    Adds a sequence to a single number space, giving out the start
 *             position for this sequence
 *
 *
 * Arg:        space [UNKN ] Undocumented argument [SingleNumberSpace *]
 * Arg:          seq [UNKN ] Undocumented argument [ShadowSequence *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 106 "singlenumberspace.dy"
int add_ShadowSequence_SingleNumberSpace(SingleNumberSpace * space,ShadowSequence * seq)
{
  int seq_size = 1;
  int i;
  SingleNumberSequence * a;
  SingleNumberSequence * dummy;

 
  assert(space);
  assert(seq);

  a = SingleNumberSequence_alloc();
  a->start = space->current_end;
  



  if( space->is_static != 0 ) {
    if( seq->seq->len > space->max_length ) {
      seq_size = (seq->seq->len / space->max_length) +1;
    } else {
      seq_size = 1;
    }

    a->end   = space->current_end + (space->max_length*seq_size);
    space->average_len = space->max_length;
  } else {
    a->end   = space->current_end + seq->seq->len;
    space->average_len = a->end / space->len;
  }

  a->seq = seq;

  add_SingleNumberSpace(space,a);

  /* pad for long sequences */
  if( seq_size > 1 ) {
    for(i=1;i<seq_size;i++) {
      dummy = SingleNumberSequence_alloc();
      dummy->start = -1;
      dummy->end = -1;
      dummy->seq = NULL;
      add_SingleNumberSpace(space,dummy);
    }
  }

  space->current_end = a->end;


  return a->start;
}

/* Function:  new_SingleNumberSpace(has_maxlen,max_length)
 *
 * Descrip:    Provides a new single number space
 *
 *
 * Arg:        has_maxlen [UNKN ] Undocumented argument [int]
 * Arg:        max_length [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [SingleNumberSpace *]
 *
 */
# line 161 "singlenumberspace.dy"
SingleNumberSpace * new_SingleNumberSpace(int has_maxlen,int max_length)
{
  SingleNumberSpace * out;

  out = SingleNumberSpace_alloc_std();

  out->is_static = has_maxlen;
  out->max_length = max_length;

  return out;
  
}





# line 183 "singlenumberspace.c"
/* Function:  hard_link_SingleNumberSequence(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SingleNumberSequence *]
 *
 * Return [UNKN ]  Undocumented return value [SingleNumberSequence *]
 *
 */
SingleNumberSequence * hard_link_SingleNumberSequence(SingleNumberSequence * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a SingleNumberSequence object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  SingleNumberSequence_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SingleNumberSequence *]
 *
 */
SingleNumberSequence * SingleNumberSequence_alloc(void) 
{
    SingleNumberSequence * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(SingleNumberSequence *) ckalloc (sizeof(SingleNumberSequence))) == NULL)    {  
      warn("SingleNumberSequence_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->start = 0;  
    out->end = 0;    


    return out;  
}    


/* Function:  free_SingleNumberSequence(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SingleNumberSequence *]
 *
 * Return [UNKN ]  Undocumented return value [SingleNumberSequence *]
 *
 */
SingleNumberSequence * free_SingleNumberSequence(SingleNumberSequence * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a SingleNumberSequence obj. Should be trappable");  
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


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_SingleNumberSpace(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_SingleNumberSpace
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [SingleNumberSequence **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_SingleNumberSpace(SingleNumberSequence ** list,int i,int j)  
{
    SingleNumberSequence * temp; 
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_SingleNumberSpace(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_SingleNumberSpace which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [SingleNumberSequence **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_SingleNumberSpace(SingleNumberSequence ** list,int left,int right,int (*comp)(SingleNumberSequence * ,SingleNumberSequence * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_SingleNumberSpace(list,left,(left+right)/2);    
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_SingleNumberSpace (list,++last,i);  
      }  
    swap_SingleNumberSpace (list,left,last); 
    qsort_SingleNumberSpace(list,left,last-1,comp);  
    qsort_SingleNumberSpace(list,last+1,right,comp); 
}    


/* Function:  sort_SingleNumberSpace(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_SingleNumberSpace
 *
 *
 * Arg:         obj [UNKN ] Object containing list [SingleNumberSpace *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_SingleNumberSpace(SingleNumberSpace * obj,int (*comp)(SingleNumberSequence *, SingleNumberSequence *)) 
{
    qsort_SingleNumberSpace(obj->sns,0,obj->len-1,comp); 
    return;  
}    


/* Function:  expand_SingleNumberSpace(obj,len)
 *
 * Descrip:    Really an internal function for add_SingleNumberSpace
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [SingleNumberSpace *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_SingleNumberSpace(SingleNumberSpace * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_SingleNumberSpace called with no need");  
      return TRUE;   
      }  


    if( (obj->sns = (SingleNumberSequence ** ) ckrealloc (obj->sns,sizeof(SingleNumberSequence *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_SingleNumberSpace, returning FALSE");    
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_SingleNumberSpace(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [SingleNumberSpace *]
 * Arg:        add [OWNER] Object to add to the list [SingleNumberSequence *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_SingleNumberSpace(SingleNumberSpace * obj,SingleNumberSequence * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_SingleNumberSpace(obj,obj->len + SingleNumberSpaceLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->sns[obj->len++]=add;    
    return TRUE; 
}    


/* Function:  flush_SingleNumberSpace(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [SingleNumberSpace *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_SingleNumberSpace(SingleNumberSpace * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->sns[i] != NULL)   {  
        free_SingleNumberSequence(obj->sns[i]);  
        obj->sns[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  SingleNumberSpace_alloc_std(void)
 *
 * Descrip:    Equivalent to SingleNumberSpace_alloc_len(SingleNumberSpaceLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SingleNumberSpace *]
 *
 */
SingleNumberSpace * SingleNumberSpace_alloc_std(void) 
{
    return SingleNumberSpace_alloc_len(SingleNumberSpaceLISTLENGTH); 
}    


/* Function:  SingleNumberSpace_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [SingleNumberSpace *]
 *
 */
SingleNumberSpace * SingleNumberSpace_alloc_len(int len) 
{
    SingleNumberSpace * out;/* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = SingleNumberSpace_alloc()) == NULL)    
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->sns = (SingleNumberSequence ** ) ckcalloc (len,sizeof(SingleNumberSequence *))) == NULL)    {  
      warn("Warning, ckcalloc failed in SingleNumberSpace_alloc_len");   
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_SingleNumberSpace(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SingleNumberSpace *]
 *
 * Return [UNKN ]  Undocumented return value [SingleNumberSpace *]
 *
 */
SingleNumberSpace * hard_link_SingleNumberSpace(SingleNumberSpace * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a SingleNumberSpace object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  SingleNumberSpace_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SingleNumberSpace *]
 *
 */
SingleNumberSpace * SingleNumberSpace_alloc(void) 
{
    SingleNumberSpace * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(SingleNumberSpace *) ckalloc (sizeof(SingleNumberSpace))) == NULL)  {  
      warn("SingleNumberSpace_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->current_end = 0;    
    out->sns = NULL; 
    out->len = out->maxlen = 0;  
    out->average_len = 0;    
    out->is_static = 0;  
    out->max_length = 1000;  


    return out;  
}    


/* Function:  free_SingleNumberSpace(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SingleNumberSpace *]
 *
 * Return [UNKN ]  Undocumented return value [SingleNumberSpace *]
 *
 */
SingleNumberSpace * free_SingleNumberSpace(SingleNumberSpace * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a SingleNumberSpace obj. Should be trappable"); 
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
    if( obj->sns != NULL)    {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->sns[i] != NULL) 
          free_SingleNumberSequence(obj->sns[i]);    
        }  
      ckfree(obj->sns);  
      }  
    /* obj->last_accessed is linked in */ 


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
