#ifdef _cplusplus
extern "C" {
#endif
#include "singleseqspace.h"



/* Function:  lookup_Sequence_SinglePosSpace(space,pos)
 *
 * Descrip:    New return using binary choping
 *
 *
 * Arg:        space [UNKN ] Undocumented argument [SinglePosSpace *]
 * Arg:          pos [UNKN ] Undocumented argument [long]
 *
 * Return [UNKN ]  Undocumented return value [SinglePosSequence *]
 *
 */
# line 36 "singleseqspace.dy"
SinglePosSequence * lookup_Sequence_SinglePosSpace(SinglePosSpace * space,long pos)
{
  int index_pos;

  assert(pos < space->current_end);
  assert(pos >= 0);

  if( space->is_static != 0 ) {
    /* can make explicit calculation for position */

    index_pos = (int) (pos / space->average_len);
    if( space->sns[index_pos]->start == -1 ) {
      /* was a padding position */
      for(index_pos--;index_pos >= 0 && space->sns[index_pos]->start == -1;index_pos--) {
	/*fprintf(stderr,"Regressed %d with %d [pos %d]\n",index_pos,space->sns[index_pos]->start,pos);*/
      }
      assert(index_pos >= 0 );
    }
	
    if( space->sns[index_pos] == NULL ) {
        fatal("Going to return %d for %d with %d - max %d\n",index_pos,pos,space->average_len,space->current_end);
    }

    return space->sns[index_pos];
  }

  
  assert(pos < space->current_end);
  
  index_pos = find_position_SinglePosSpace(space,0,space->len-1,pos);

  
  return space->sns[index_pos];
}


/* Function:  find_position_SinglePosSpace(space,lower,higher,position)
 *
 * Descrip:    Recursive function for finding position
 *
 *
 * Arg:           space [UNKN ] Undocumented argument [SinglePosSpace *]
 * Arg:           lower [UNKN ] Undocumented argument [int]
 * Arg:          higher [UNKN ] Undocumented argument [int]
 * Arg:        position [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 75 "singleseqspace.dy"
int find_position_SinglePosSpace(SinglePosSpace * space,int lower,int higher,int position)
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
    return find_position_SinglePosSpace(space,mid,higher,position);
  } else {
    return find_position_SinglePosSpace(space,lower,mid,position);
  }

}


/* Function:  add_Sequence_SinglePosSpace(space,length,data)
 *
 * Descrip:    Adds a sequence to a single number space, giving out the start
 *             position for this sequence
 *
 *
 * Arg:         space [UNKN ] Undocumented argument [SinglePosSpace *]
 * Arg:        length [UNKN ] Undocumented argument [long int]
 * Arg:          data [UNKN ] Undocumented argument [void *]
 *
 * Return [UNKN ]  Undocumented return value [long int]
 *
 */
# line 109 "singleseqspace.dy"
long int add_Sequence_SinglePosSpace(SinglePosSpace * space,long int length,void * data)
{
  int seq_size = 1;
  int i;
  SinglePosSequence * a;
  SinglePosSequence * dummy;


 
  assert(space);
  assert(data);
  assert(length >= 0);

  a = SinglePosSequence_alloc();
  a->start = space->current_end;
  
  if( space->is_static != 0 ) {
    if( length > space->max_length ) {
      seq_size = (length / space->max_length) +1;
    } else {
      seq_size = 1;
    }

    a->end   = space->current_end + (space->max_length*seq_size);
    space->average_len = space->max_length;
  } else {
    a->end   = space->current_end + length;
    space->average_len = a->end / length;
  }

  a->data = data;

  add_SinglePosSpace(space,a);

  /* pad for long sequences */
  if( seq_size > 1 ) {
    for(i=1;i<seq_size;i++) {
      dummy = SinglePosSequence_alloc();
      dummy->start = -1;
      dummy->end = -1;
      dummy->data = NULL;
      add_SinglePosSpace(space,dummy);
    }
  }

  space->current_end = a->end;


  return a->start;
}

/* Function:  new_SinglePosSpace(has_maxlen,max_length)
 *
 * Descrip:    Provides a new single number space
 *
 *
 * Arg:        has_maxlen [UNKN ] Undocumented argument [int]
 * Arg:        max_length [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [SinglePosSpace *]
 *
 */
# line 163 "singleseqspace.dy"
SinglePosSpace * new_SinglePosSpace(int has_maxlen,int max_length)
{
  SinglePosSpace * out;

  out = SinglePosSpace_alloc_std();

  out->is_static = has_maxlen;
  out->max_length = max_length;

  return out;
  
}

# line 181 "singleseqspace.c"
/* Function:  hard_link_SinglePosSequence(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SinglePosSequence *]
 *
 * Return [UNKN ]  Undocumented return value [SinglePosSequence *]
 *
 */
SinglePosSequence * hard_link_SinglePosSequence(SinglePosSequence * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a SinglePosSequence object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  SinglePosSequence_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SinglePosSequence *]
 *
 */
SinglePosSequence * SinglePosSequence_alloc(void) 
{
    SinglePosSequence * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(SinglePosSequence *) ckalloc (sizeof(SinglePosSequence))) == NULL)  {  
      warn("SinglePosSequence_alloc failed ");   
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


/* Function:  free_SinglePosSequence(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SinglePosSequence *]
 *
 * Return [UNKN ]  Undocumented return value [SinglePosSequence *]
 *
 */
SinglePosSequence * free_SinglePosSequence(SinglePosSequence * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a SinglePosSequence obj. Should be trappable"); 
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
    /* obj->data is linked in */ 


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_SinglePosSpace(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_SinglePosSpace
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [SinglePosSequence **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_SinglePosSpace(SinglePosSequence ** list,int i,int j)  
{
    SinglePosSequence * temp;    
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_SinglePosSpace(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_SinglePosSpace which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [SinglePosSequence **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_SinglePosSpace(SinglePosSequence ** list,int left,int right,int (*comp)(SinglePosSequence * ,SinglePosSequence * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_SinglePosSpace(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_SinglePosSpace (list,++last,i); 
      }  
    swap_SinglePosSpace (list,left,last);    
    qsort_SinglePosSpace(list,left,last-1,comp); 
    qsort_SinglePosSpace(list,last+1,right,comp);    
}    


/* Function:  sort_SinglePosSpace(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_SinglePosSpace
 *
 *
 * Arg:         obj [UNKN ] Object containing list [SinglePosSpace *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_SinglePosSpace(SinglePosSpace * obj,int (*comp)(SinglePosSequence *, SinglePosSequence *)) 
{
    qsort_SinglePosSpace(obj->sns,0,obj->len-1,comp);    
    return;  
}    


/* Function:  expand_SinglePosSpace(obj,len)
 *
 * Descrip:    Really an internal function for add_SinglePosSpace
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [SinglePosSpace *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_SinglePosSpace(SinglePosSpace * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_SinglePosSpace called with no need"); 
      return TRUE;   
      }  


    if( (obj->sns = (SinglePosSequence ** ) ckrealloc (obj->sns,sizeof(SinglePosSequence *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_SinglePosSpace, returning FALSE");   
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_SinglePosSpace(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [SinglePosSpace *]
 * Arg:        add [OWNER] Object to add to the list [SinglePosSequence *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_SinglePosSpace(SinglePosSpace * obj,SinglePosSequence * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_SinglePosSpace(obj,obj->len + SinglePosSpaceLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->sns[obj->len++]=add;    
    return TRUE; 
}    


/* Function:  flush_SinglePosSpace(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [SinglePosSpace *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_SinglePosSpace(SinglePosSpace * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->sns[i] != NULL)   {  
        free_SinglePosSequence(obj->sns[i]); 
        obj->sns[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  SinglePosSpace_alloc_std(void)
 *
 * Descrip:    Equivalent to SinglePosSpace_alloc_len(SinglePosSpaceLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SinglePosSpace *]
 *
 */
SinglePosSpace * SinglePosSpace_alloc_std(void) 
{
    return SinglePosSpace_alloc_len(SinglePosSpaceLISTLENGTH);   
}    


/* Function:  SinglePosSpace_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [SinglePosSpace *]
 *
 */
SinglePosSpace * SinglePosSpace_alloc_len(int len) 
{
    SinglePosSpace * out;   /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = SinglePosSpace_alloc()) == NULL)   
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->sns = (SinglePosSequence ** ) ckcalloc (len,sizeof(SinglePosSequence *))) == NULL)  {  
      warn("Warning, ckcalloc failed in SinglePosSpace_alloc_len");  
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_SinglePosSpace(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SinglePosSpace *]
 *
 * Return [UNKN ]  Undocumented return value [SinglePosSpace *]
 *
 */
SinglePosSpace * hard_link_SinglePosSpace(SinglePosSpace * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a SinglePosSpace object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  SinglePosSpace_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SinglePosSpace *]
 *
 */
SinglePosSpace * SinglePosSpace_alloc(void) 
{
    SinglePosSpace * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(SinglePosSpace *) ckalloc (sizeof(SinglePosSpace))) == NULL)    {  
      warn("SinglePosSpace_alloc failed ");  
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


/* Function:  free_SinglePosSpace(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SinglePosSpace *]
 *
 * Return [UNKN ]  Undocumented return value [SinglePosSpace *]
 *
 */
SinglePosSpace * free_SinglePosSpace(SinglePosSpace * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a SinglePosSpace obj. Should be trappable");    
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
          free_SinglePosSequence(obj->sns[i]);   
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
