#ifdef _cplusplus
extern "C" {
#endif
#include "vectorindex.h"



# line 28 "vectorindex.dy"
double vector_KmerVectorPosition(KmerVectorPosition * a,KmerVectorPosition * b)
{
  int i;
  double out = 0.0;

  for(i=0;i<KMER_VECTOR_LENGTH;i++) {
    out += a->kmer_vector[i] * b->kmer_vector[i] / 256;
  }

  return out;
}

# line 40 "vectorindex.dy"
KmerVectorPosition * new_KmerVectorPosition(void)
{
  KmerVectorPosition * out;

  out = malloc(sizeof(KmerVectorPosition));

  bzero(out->kmer_vector,KMER_VECTOR_LENGTH);

  return out;
}


# line 52 "vectorindex.dy"
KmerVectorPositionSet * build_KmerVectorPositionSet(Sequence * input,long int start_pos,int step_size)
{
  int i;
  int j;
  int k;
  int base[KMER_VECTOR_SIZE];
  int temp;

  int is_forward;
  int curr;

  int forward;
  int backward;

  int vec;

  KmerVectorPosition * pos;
  KmerVectorPositionSet * out;

  assert(input != NULL);

  temp = 1;
  for(i=0;i<KMER_VECTOR_SIZE;i++) {
    base[i] = temp;
    temp = temp * 4;
  }

  out= KmerVectorPositionSet_alloc_len((input->len / step_size)+2);
  
  for(i=0;i+step_size<input->len;) {
    pos = new_KmerVectorPosition();
    pos->position = start_pos + i;
    for(j=0;j<256 - KMER_VECTOR_SIZE;j++) {
      curr = i+j;
      /* first find lexical forward/backwardness */
      is_forward = 1; /* palindromes read forward */
      for(k=0;k<KMER_VECTOR_SIZE;k++) {
	forward  = base_from_char(input->seq[curr+k]);
	backward = complement_base(base_from_char(input->seq[curr+KMER_VECTOR_SIZE-k]));
	if( forward == backward ) {
	  continue;
	}
	if( forward > backward ) {
	  is_forward = 1;
	} else {
	  is_forward = 0;
	}
	break;
      }
      
      vec = 0;
      for(k=0;k<KMER_VECTOR_SIZE;k++) {
	if( is_forward ) {
	  vec += base[k] * base_from_char(input->seq[curr+k]);
	} else {
	  vec += base[k] * complement_base(base_from_char(input->seq[curr+KMER_VECTOR_SIZE-k]));
	}
      }

      assert(vec < KMER_VECTOR_LENGTH);
      pos->kmer_vector[vec]++;

      
      fatal("Ewan has not finished this yet! Idiot!");

    }
  }


  return out;
}


# line 125 "vectorindex.dy"
KmerVectorPosition * free_KmerVectorPosition(KmerVectorPosition * vec)
{
  free(vec);
  return NULL;
}

# line 114 "vectorindex.c"
/* Function:  swap_KmerVectorPositionSet(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_KmerVectorPositionSet
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [KmerVectorPosition **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_KmerVectorPositionSet(KmerVectorPosition ** list,int i,int j)  
{
    KmerVectorPosition * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_KmerVectorPositionSet(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_KmerVectorPositionSet which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [KmerVectorPosition **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_KmerVectorPositionSet(KmerVectorPosition ** list,int left,int right,int (*comp)(KmerVectorPosition * ,KmerVectorPosition * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_KmerVectorPositionSet(list,left,(left+right)/2);    
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_KmerVectorPositionSet (list,++last,i);  
      }  
    swap_KmerVectorPositionSet (list,left,last); 
    qsort_KmerVectorPositionSet(list,left,last-1,comp);  
    qsort_KmerVectorPositionSet(list,last+1,right,comp); 
}    


/* Function:  sort_KmerVectorPositionSet(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_KmerVectorPositionSet
 *
 *
 * Arg:         obj [UNKN ] Object containing list [KmerVectorPositionSet *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_KmerVectorPositionSet(KmerVectorPositionSet * obj,int (*comp)(KmerVectorPosition *, KmerVectorPosition *)) 
{
    qsort_KmerVectorPositionSet(obj->vec,0,obj->len-1,comp); 
    return;  
}    


/* Function:  expand_KmerVectorPositionSet(obj,len)
 *
 * Descrip:    Really an internal function for add_KmerVectorPositionSet
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [KmerVectorPositionSet *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_KmerVectorPositionSet(KmerVectorPositionSet * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_KmerVectorPositionSet called with no need");  
      return TRUE;   
      }  


    if( (obj->vec = (KmerVectorPosition ** ) ckrealloc (obj->vec,sizeof(KmerVectorPosition *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_KmerVectorPositionSet, returning FALSE");    
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_KmerVectorPositionSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [KmerVectorPositionSet *]
 * Arg:        add [OWNER] Object to add to the list [KmerVectorPosition *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_KmerVectorPositionSet(KmerVectorPositionSet * obj,KmerVectorPosition * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_KmerVectorPositionSet(obj,obj->len + KmerVectorPositionSetLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->vec[obj->len++]=add;    
    return TRUE; 
}    


/* Function:  flush_KmerVectorPositionSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [KmerVectorPositionSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_KmerVectorPositionSet(KmerVectorPositionSet * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->vec[i] != NULL)   {  
        free_KmerVectorPosition(obj->vec[i]);    
        obj->vec[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  KmerVectorPositionSet_alloc_std(void)
 *
 * Descrip:    Equivalent to KmerVectorPositionSet_alloc_len(KmerVectorPositionSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [KmerVectorPositionSet *]
 *
 */
KmerVectorPositionSet * KmerVectorPositionSet_alloc_std(void) 
{
    return KmerVectorPositionSet_alloc_len(KmerVectorPositionSetLISTLENGTH); 
}    


/* Function:  KmerVectorPositionSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [KmerVectorPositionSet *]
 *
 */
KmerVectorPositionSet * KmerVectorPositionSet_alloc_len(int len) 
{
    KmerVectorPositionSet * out;/* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = KmerVectorPositionSet_alloc()) == NULL)    
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->vec = (KmerVectorPosition ** ) ckcalloc (len,sizeof(KmerVectorPosition *))) == NULL)    {  
      warn("Warning, ckcalloc failed in KmerVectorPositionSet_alloc_len");   
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_KmerVectorPositionSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [KmerVectorPositionSet *]
 *
 * Return [UNKN ]  Undocumented return value [KmerVectorPositionSet *]
 *
 */
KmerVectorPositionSet * hard_link_KmerVectorPositionSet(KmerVectorPositionSet * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a KmerVectorPositionSet object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  KmerVectorPositionSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [KmerVectorPositionSet *]
 *
 */
KmerVectorPositionSet * KmerVectorPositionSet_alloc(void) 
{
    KmerVectorPositionSet * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(KmerVectorPositionSet *) ckalloc (sizeof(KmerVectorPositionSet))) == NULL)  {  
      warn("KmerVectorPositionSet_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->vec = NULL; 
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_KmerVectorPositionSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [KmerVectorPositionSet *]
 *
 * Return [UNKN ]  Undocumented return value [KmerVectorPositionSet *]
 *
 */
KmerVectorPositionSet * free_KmerVectorPositionSet(KmerVectorPositionSet * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a KmerVectorPositionSet obj. Should be trappable"); 
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
    if( obj->vec != NULL)    {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->vec[i] != NULL) 
          free_KmerVectorPosition(obj->vec[i]);  
        }  
      ckfree(obj->vec);  
      }  


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
