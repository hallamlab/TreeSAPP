#ifdef _cplusplus
extern "C" {
#endif
#include "hsp.h"




/* Function:  on_HSP(test,query_pos,target_pos)
 *
 * Descrip:    tests whether this point is on this test HSP
 *
 *
 * Arg:              test [UNKN ] Undocumented argument [HSP *]
 * Arg:         query_pos [UNKN ] Undocumented argument [int]
 * Arg:        target_pos [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 64 "hsp.dy"
boolean on_HSP(HSP * test,int query_pos,int target_pos)
{
  if( (test->query_start - test->target_start) != query_pos - target_pos || query_pos < test->query_start || target_pos < test->target_start) {
    return FALSE;
  }

  if( query_pos - test->query_start > test->length ) {
    return FALSE;
  }
  return TRUE;
}


/* Function:  compare_HSPset_score_qsort(a,b)
 *
 * Descrip:    sorting linear HSPsets via qsort function
 *
 *
 * Arg:        a [UNKN ] Undocumented argument [const void *]
 * Arg:        b [UNKN ] Undocumented argument [const void *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 80 "hsp.dy"
int compare_HSPset_score_qsort(const void * a,const void * b)
{
  const HSPset * one; 
  const HSPset * two;

  one = *((HSPset**)a);
  two = *((HSPset**)b);


  return two->score - one->score;
}

/* Function:  compare_HSPset_score(one,two)
 *
 * Descrip:    internal function for sort linear HSPsets
 *
 *
 * Arg:        one [UNKN ] Undocumented argument [HSPset *]
 * Arg:        two [UNKN ] Undocumented argument [HSPset *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 95 "hsp.dy"
int compare_HSPset_score(HSPset * one,HSPset * two)
{
  return two->score - one->score;
}


/* Function:  sort_HSPset_by_score(*set)
 *
 * Descrip:    Sorts by score
 *
 *
 * Arg:        *set [UNKN ] Undocumented argument [HSPset]
 *
 */
# line 104 "hsp.dy"
void sort_HSPset_by_score(HSPset *set)
{
  sort_HSPset(set,compare_HSP_score);
}

/* Function:  compare_HSP_score(one,two)
 *
 * Descrip:    internal for sort by score
 *
 *
 * Arg:        one [UNKN ] Undocumented argument [HSP *]
 * Arg:        two [UNKN ] Undocumented argument [HSP *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 113 "hsp.dy"
int compare_HSP_score(HSP * one,HSP * two)
{
  return two->score - one->score;
}

/* Function:  new_dna_identical_HSP(query,target,query_pos,target_pos,target_reverse)
 *
 * Descrip:    builds a new HSP for these sequences breaking at first mismatch
 *
 *
 * Arg:                 query [UNKN ] Undocumented argument [Sequence *]
 * Arg:                target [UNKN ] Undocumented argument [Sequence *]
 * Arg:             query_pos [UNKN ] Undocumented argument [int]
 * Arg:            target_pos [UNKN ] Undocumented argument [int]
 * Arg:        target_reverse [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [HSP *]
 *
 */
# line 121 "hsp.dy"
HSP * new_dna_identical_HSP(Sequence * query,Sequence * target,int query_pos,int target_pos,int target_reverse)
{
  int i,j;
  HSP * out;

  assert(query);
  assert(target);
  assert(query_pos >= 0 && query_pos < query->len);
  assert(target_pos >= 0 && target_pos < target->len);

  out = HSP_alloc();

  out->query = hard_link_Sequence(query);
  out->target= hard_link_Sequence(target);

  i = query_pos;
  j = target_pos;

  /* upstream first */


  while( i >= 0 && j >= 0 && j <= target->len ) {
/*  fprintf(stderr," in upstream loop, looking at %d (is revd %d)\n",i,target_reverse); */
    if( target_reverse == 1 ) {
      if( toupper(query->seq[i]) != char_complement_base(toupper(target->seq[j])) ) {
	/*	fprintf(stderr,"Breaking, reversed  at %d,%d %c,%c\n",i,j,query->seq[i],char_complement_base(toupper(target->seq[j]))); */

	break;
      } else {
	i--;
	j++;
      }
    } else {
      if( toupper(query->seq[i]) != toupper(target->seq[j]) ) {
	/*	fprintf(stderr,"Breaking at %d,%d %c,%c\n",i,j,query->seq[i],target->seq[j]); */
	break;
      } else{
	i--;
	j--;
      }
    }
  }

  if( target_reverse == 1 ) {
    i++;
    j--;
  } else {
    i++;
    j++;
  }
  
  out->query_start = i;
  out->target_start = j;

  i = query_pos;
  j = target_pos;

  /* downstream next */

  while( i < query->len && j >= 0 && j < target->len ) {
    if( target_reverse == 1 ) {
      if( toupper(query->seq[i]) != char_complement_base(toupper(target->seq[j])) ) {
	i--;
	j++;
	break;
      } else {
	i++;
	j--;
      }
    } else {
      if( toupper(query->seq[i]) != toupper(target->seq[j]) ) {
	i--;
	j--;
	break;
      } else{
	i++;
	j++;
      }
    }
  }

  if( target_reverse == 1 ) {
    i++;
    j--;
  } else {
    i--;
    j--;
  }

  out->length = i - out->query_start;
  out->score = out->length;

  out->target_reverse = target_reverse;

  return out;
}

/* Function:  new_HSP(cache,query,target,query_pos,target_pos,mat,drop_off)
 *
 * Descrip:    builds a new HSP for these sequences
 *
 *
 * Arg:             cache [UNKN ] Undocumented argument [HSPCache *]
 * Arg:             query [UNKN ] Undocumented argument [Sequence *]
 * Arg:            target [UNKN ] Undocumented argument [Sequence *]
 * Arg:         query_pos [UNKN ] Undocumented argument [int]
 * Arg:        target_pos [UNKN ] Undocumented argument [int]
 * Arg:               mat [UNKN ] Undocumented argument [CompMat *]
 * Arg:          drop_off [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [HSP *]
 *
 */
# line 221 "hsp.dy"
HSP * new_HSP(HSPCache * cache,Sequence * query,Sequence * target,int query_pos,int target_pos,CompMat * mat,int drop_off)
{
  int i,j;
  int pause_i,pause_j;
  int score = 0;
  int max_score= 0;
  HSP * out;



  if( cache != NULL ) {
    out = HSP_alloc_cache(cache);
  } else {
    out = HSP_alloc();
  }

  assert(out != NULL);
  assert(query != NULL);
  assert(query->seq != NULL);
  assert(target != NULL);
  assert(target->seq != NULL);
  assert(mat != NULL);
  assert(drop_off > 0);
  assert(query->seq[query_pos] >= 'A' && query->seq[query_pos] <= 'Z');
  assert(target->seq[target_pos] >= 'A' && target->seq[target_pos] <= 'Z');


  /*  fprintf(stdout,"got new hsp with drop_off of %d\n",drop_off);*/

  out->query = hard_link_Sequence(query);
  out->target= hard_link_Sequence(target);
  out->score = -(mat->comp[query->seq[query_pos]-'A'][target->seq[target_pos]-'A']);

  pause_i = i = query_pos;
  pause_j = j = target_pos;

  /* go upstream first */
  
  for(score = 0,max_score = 0;(max_score - score) < drop_off && i >= 0 && j >= 0;i--,j--) {

    score += mat->comp[query->seq[i]-'A'][target->seq[j]-'A'];
       /*fprintf(stderr,"Examining %c,%c - score %d max_score %d\n",query->seq[i],target->seq[j],score,max_score); */

    if( score > max_score ) {
      pause_i = i;
      pause_j = j;
      max_score = score;
    }
  }
  

  out->query_start = pause_i;
  out->target_start = pause_j;
  out->score += max_score;

  /* now downstream */

  for(score=0,max_score=0,pause_i=i=query_pos,pause_j=j=target_pos;(max_score - score) < drop_off && i < query->len && j < target->len;i++,j++) {
    score += mat->comp[query->seq[i]-'A'][target->seq[j]-'A'];
    if( score > max_score ) {
      pause_i = i;
      pause_j = j;
      max_score = score;
    }
  }

  out->length = pause_i - out->query_start+1;

  out->score += max_score;

  return out;
}

/* Function:  HSP_alloc_cache(hspc)
 *
 * Descrip:    Returns new HSP, using cache if needed
 *
 *
 * Arg:        hspc [UNKN ] Undocumented argument [HSPCache *]
 *
 * Return [UNKN ]  Undocumented return value [HSP *]
 *
 */
# line 297 "hsp.dy"
HSP * HSP_alloc_cache(HSPCache * hspc)
{
  assert(hspc);

  if( hspc->len > 0 ) {
    hspc->len--;
    return hspc->cache[hspc->len];
  } 

  return HSP_alloc();
}

/* Function:  free_HSP_cache(cache,hsp)
 *
 * Descrip:    Places HSP back into cache, freeing if necessary
 *
 *
 * Arg:        cache [UNKN ] Undocumented argument [HSPCache *]
 * Arg:          hsp [UNKN ] Undocumented argument [HSP *]
 *
 * Return [UNKN ]  Undocumented return value [HSP *]
 *
 */
# line 312 "hsp.dy"
HSP * free_HSP_cache(HSPCache * cache,HSP * hsp)
{
  int return_early = 0;

  assert(cache);
  assert(hsp);

#ifdef PTHREAD   
  assert(pthread_mutex_lock(&(hsp->dynamite_mutex)) == 0); 
#endif   
  if( hsp->dynamite_hard_link > 1)     {  
    return_early = 1;  
    hsp->dynamite_hard_link--; 
  }  
#ifdef PTHREAD   
  assert(pthread_mutex_unlock(&(hsp->dynamite_mutex)) == 0);   
#endif   
  
  if( return_early == 1)   
    return NULL;   
  
  if( cache->len > cache->max_cache ) {
    return free_HSP(hsp);
  } else {
    add_HSPCache(cache,hsp);
    return NULL;
  }
}

/* Function:  new_HSPCache(maxsize)
 *
 * Descrip:    Makes a new cache
 *
 *
 * Arg:        maxsize [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [HSPCache *]
 *
 */
# line 344 "hsp.dy"
HSPCache * new_HSPCache(int maxsize)
{
  HSPCache * out;

  out = HSPCache_alloc_std();
  out->max_cache = maxsize;

  return out;
}

/* Function:  show_HSPset(s,ofp)
 *
 * Descrip:    Shows a HSP set
 *
 *
 * Arg:          s [UNKN ] Undocumented argument [HSPset *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 357 "hsp.dy"
void show_HSPset(HSPset * s,FILE * ofp)
{
  int i;

  for(i=0;i<s->len;i++) {
    if( s->hsp[i]->target_reverse == 0 ) {
      fprintf(ofp,"%5d %s\t%d\t%d\t%s\t%d\t%d\t+\n",s->hsp[i]->length,s->hsp[i]->query->name,s->hsp[i]->query_start,s->hsp[i]->query_start+s->hsp[i]->length,s->hsp[i]->target->name,s->hsp[i]->target_start,s->hsp[i]->target_start + s->hsp[i]->length);
    } else {
      fprintf(ofp,"%5d %s\t%d\t%d\t%s\t%d\t%d\t-\n",s->hsp[i]->length,s->hsp[i]->query->name,s->hsp[i]->query_start,s->hsp[i]->query_start+s->hsp[i]->length,s->hsp[i]->target->name,s->hsp[i]->target_start,s->hsp[i]->target_start - s->hsp[i]->length);
    }
  }


}

/* Function:  show_HSP(hsp,linelength,out)
 *
 * Descrip:    Shows an HSP
 *
 *
 * Arg:               hsp [UNKN ] Undocumented argument [HSP *]
 * Arg:        linelength [UNKN ] Undocumented argument [int]
 * Arg:               out [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 375 "hsp.dy"
void show_HSP(HSP * hsp,int linelength,FILE * out)
{
  int i,j;

  for(i=hsp->query_start,j=hsp->target_start;(i-hsp->query_start) < hsp->length;i++,j++) {
    fprintf(out,"%c %c\n",hsp->query->seq[i],hsp->target->seq[j]);
  }

}



# line 418 "hsp.c"
/* Function:  hard_link_HSP(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [HSP *]
 *
 * Return [UNKN ]  Undocumented return value [HSP *]
 *
 */
HSP * hard_link_HSP(HSP * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a HSP object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  HSP_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [HSP *]
 *
 */
HSP * HSP_alloc(void) 
{
    HSP * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(HSP *) ckalloc (sizeof(HSP))) == NULL)  {  
      warn("HSP_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->query = NULL;   
    out->target = NULL;  
    out->query_start = 0;    
    out->target_start = 0;   
    out->length = 0; 
    out->score = 0;  
    out->is_in_block = 0;    
    out->target_reverse = 0; 


    return out;  
}    


/* Function:  free_HSP(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [HSP *]
 *
 * Return [UNKN ]  Undocumented return value [HSP *]
 *
 */
HSP * free_HSP(HSP * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a HSP obj. Should be trappable");   
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
    if( obj->query != NULL)  
      free_Sequence(obj->query);     
    if( obj->target != NULL) 
      free_Sequence(obj->target);    


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_HSPCache(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_HSPCache
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [HSP **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_HSPCache(HSP ** list,int i,int j)  
{
    HSP * temp;  
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_HSPCache(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_HSPCache which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [HSP **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_HSPCache(HSP ** list,int left,int right,int (*comp)(HSP * ,HSP * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_HSPCache(list,left,(left+right)/2); 
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_HSPCache (list,++last,i);   
      }  
    swap_HSPCache (list,left,last);  
    qsort_HSPCache(list,left,last-1,comp);   
    qsort_HSPCache(list,last+1,right,comp);  
}    


/* Function:  sort_HSPCache(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_HSPCache
 *
 *
 * Arg:         obj [UNKN ] Object containing list [HSPCache *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_HSPCache(HSPCache * obj,int (*comp)(HSP *, HSP *)) 
{
    qsort_HSPCache(obj->cache,0,obj->len-1,comp);    
    return;  
}    


/* Function:  expand_HSPCache(obj,len)
 *
 * Descrip:    Really an internal function for add_HSPCache
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [HSPCache *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_HSPCache(HSPCache * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_HSPCache called with no need");   
      return TRUE;   
      }  


    if( (obj->cache = (HSP ** ) ckrealloc (obj->cache,sizeof(HSP *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_HSPCache, returning FALSE"); 
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_HSPCache(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [HSPCache *]
 * Arg:        add [OWNER] Object to add to the list [HSP *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_HSPCache(HSPCache * obj,HSP * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_HSPCache(obj,obj->len + HSPCacheLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->cache[obj->len++]=add;  
    return TRUE; 
}    


/* Function:  flush_HSPCache(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [HSPCache *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_HSPCache(HSPCache * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->cache[i] != NULL) {  
        free_HSP(obj->cache[i]); 
        obj->cache[i] = NULL;    
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  HSPCache_alloc_std(void)
 *
 * Descrip:    Equivalent to HSPCache_alloc_len(HSPCacheLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [HSPCache *]
 *
 */
HSPCache * HSPCache_alloc_std(void) 
{
    return HSPCache_alloc_len(HSPCacheLISTLENGTH);   
}    


/* Function:  HSPCache_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [HSPCache *]
 *
 */
HSPCache * HSPCache_alloc_len(int len) 
{
    HSPCache * out; /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = HSPCache_alloc()) == NULL) 
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->cache = (HSP ** ) ckcalloc (len,sizeof(HSP *))) == NULL)    {  
      warn("Warning, ckcalloc failed in HSPCache_alloc_len");    
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_HSPCache(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [HSPCache *]
 *
 * Return [UNKN ]  Undocumented return value [HSPCache *]
 *
 */
HSPCache * hard_link_HSPCache(HSPCache * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a HSPCache object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  HSPCache_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [HSPCache *]
 *
 */
HSPCache * HSPCache_alloc(void) 
{
    HSPCache * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(HSPCache *) ckalloc (sizeof(HSPCache))) == NULL)    {  
      warn("HSPCache_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->cache = NULL;   
    out->len = out->maxlen = 0;  
    out->max_cache = 0;  


    return out;  
}    


/* Function:  free_HSPCache(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [HSPCache *]
 *
 * Return [UNKN ]  Undocumented return value [HSPCache *]
 *
 */
HSPCache * free_HSPCache(HSPCache * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a HSPCache obj. Should be trappable");  
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
    if( obj->cache != NULL)  {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->cache[i] != NULL)   
          free_HSP(obj->cache[i]);   
        }  
      ckfree(obj->cache);    
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_HSPset(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_HSPset
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [HSP **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_HSPset(HSP ** list,int i,int j)  
{
    HSP * temp;  
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_HSPset(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_HSPset which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [HSP **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_HSPset(HSP ** list,int left,int right,int (*comp)(HSP * ,HSP * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_HSPset(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_HSPset (list,++last,i); 
      }  
    swap_HSPset (list,left,last);    
    qsort_HSPset(list,left,last-1,comp); 
    qsort_HSPset(list,last+1,right,comp);    
}    


/* Function:  sort_HSPset(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_HSPset
 *
 *
 * Arg:         obj [UNKN ] Object containing list [HSPset *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_HSPset(HSPset * obj,int (*comp)(HSP *, HSP *)) 
{
    qsort_HSPset(obj->hsp,0,obj->len-1,comp);    
    return;  
}    


/* Function:  expand_HSPset(obj,len)
 *
 * Descrip:    Really an internal function for add_HSPset
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [HSPset *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_HSPset(HSPset * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_HSPset called with no need"); 
      return TRUE;   
      }  


    if( (obj->hsp = (HSP ** ) ckrealloc (obj->hsp,sizeof(HSP *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_HSPset, returning FALSE");   
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_HSPset(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [HSPset *]
 * Arg:        add [OWNER] Object to add to the list [HSP *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_HSPset(HSPset * obj,HSP * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_HSPset(obj,obj->len + HSPsetLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->hsp[obj->len++]=add;    
    return TRUE; 
}    


/* Function:  flush_HSPset(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [HSPset *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_HSPset(HSPset * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->hsp[i] != NULL)   {  
        free_HSP(obj->hsp[i]);   
        obj->hsp[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  HSPset_alloc_std(void)
 *
 * Descrip:    Equivalent to HSPset_alloc_len(HSPsetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [HSPset *]
 *
 */
HSPset * HSPset_alloc_std(void) 
{
    return HSPset_alloc_len(HSPsetLISTLENGTH);   
}    


/* Function:  HSPset_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [HSPset *]
 *
 */
HSPset * HSPset_alloc_len(int len) 
{
    HSPset * out;   /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = HSPset_alloc()) == NULL)   
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->hsp = (HSP ** ) ckcalloc (len,sizeof(HSP *))) == NULL)  {  
      warn("Warning, ckcalloc failed in HSPset_alloc_len");  
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_HSPset(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [HSPset *]
 *
 * Return [UNKN ]  Undocumented return value [HSPset *]
 *
 */
HSPset * hard_link_HSPset(HSPset * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a HSPset object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  HSPset_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [HSPset *]
 *
 */
HSPset * HSPset_alloc(void) 
{
    HSPset * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(HSPset *) ckalloc (sizeof(HSPset))) == NULL)    {  
      warn("HSPset_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->hsp = NULL; 
    out->len = out->maxlen = 0;  
    out->score = 0;  
    out->best_score = 0; 
    out->last_accessed = -1; 


    return out;  
}    


/* Function:  free_HSPset(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [HSPset *]
 *
 * Return [UNKN ]  Undocumented return value [HSPset *]
 *
 */
HSPset * free_HSPset(HSPset * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a HSPset obj. Should be trappable");    
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
    if( obj->hsp != NULL)    {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->hsp[i] != NULL) 
          free_HSP(obj->hsp[i]); 
        }  
      ckfree(obj->hsp);  
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_LinearHSPmanager(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_LinearHSPmanager
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [HSPset **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_LinearHSPmanager(HSPset ** list,int i,int j)  
{
    HSPset * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_LinearHSPmanager(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_LinearHSPmanager which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [HSPset **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_LinearHSPmanager(HSPset ** list,int left,int right,int (*comp)(HSPset * ,HSPset * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_LinearHSPmanager(list,left,(left+right)/2); 
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_LinearHSPmanager (list,++last,i);   
      }  
    swap_LinearHSPmanager (list,left,last);  
    qsort_LinearHSPmanager(list,left,last-1,comp);   
    qsort_LinearHSPmanager(list,last+1,right,comp);  
}    


/* Function:  sort_LinearHSPmanager(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_LinearHSPmanager
 *
 *
 * Arg:         obj [UNKN ] Object containing list [LinearHSPmanager *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_LinearHSPmanager(LinearHSPmanager * obj,int (*comp)(HSPset *, HSPset *)) 
{
    qsort_LinearHSPmanager(obj->set,0,obj->len-1,comp);  
    return;  
}    


/* Function:  expand_LinearHSPmanager(obj,len)
 *
 * Descrip:    Really an internal function for add_LinearHSPmanager
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [LinearHSPmanager *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_LinearHSPmanager(LinearHSPmanager * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_LinearHSPmanager called with no need");   
      return TRUE;   
      }  


    if( (obj->set = (HSPset ** ) ckrealloc (obj->set,sizeof(HSPset *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_LinearHSPmanager, returning FALSE"); 
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_LinearHSPmanager(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [LinearHSPmanager *]
 * Arg:        add [OWNER] Object to add to the list [HSPset *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_LinearHSPmanager(LinearHSPmanager * obj,HSPset * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_LinearHSPmanager(obj,obj->len + LinearHSPmanagerLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->set[obj->len++]=add;    
    return TRUE; 
}    


/* Function:  flush_LinearHSPmanager(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [LinearHSPmanager *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_LinearHSPmanager(LinearHSPmanager * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->set[i] != NULL)   {  
        free_HSPset(obj->set[i]);    
        obj->set[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  LinearHSPmanager_alloc_std(void)
 *
 * Descrip:    Equivalent to LinearHSPmanager_alloc_len(LinearHSPmanagerLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LinearHSPmanager *]
 *
 */
LinearHSPmanager * LinearHSPmanager_alloc_std(void) 
{
    return LinearHSPmanager_alloc_len(LinearHSPmanagerLISTLENGTH);   
}    


/* Function:  LinearHSPmanager_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [LinearHSPmanager *]
 *
 */
LinearHSPmanager * LinearHSPmanager_alloc_len(int len) 
{
    LinearHSPmanager * out; /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = LinearHSPmanager_alloc()) == NULL) 
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->set = (HSPset ** ) ckcalloc (len,sizeof(HSPset *))) == NULL)    {  
      warn("Warning, ckcalloc failed in LinearHSPmanager_alloc_len");    
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_LinearHSPmanager(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [LinearHSPmanager *]
 *
 * Return [UNKN ]  Undocumented return value [LinearHSPmanager *]
 *
 */
LinearHSPmanager * hard_link_LinearHSPmanager(LinearHSPmanager * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a LinearHSPmanager object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  LinearHSPmanager_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LinearHSPmanager *]
 *
 */
LinearHSPmanager * LinearHSPmanager_alloc(void) 
{
    LinearHSPmanager * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(LinearHSPmanager *) ckalloc (sizeof(LinearHSPmanager))) == NULL)    {  
      warn("LinearHSPmanager_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->set = NULL; 
    out->len = out->maxlen = 0;  
    out->query = NULL;   
    out->mat = NULL; 
    out->min_score = 0;  
    out->width = 0;  
    out->tail = 0;   
    out->worst_hsp_score = 0;    


    return out;  
}    


/* Function:  free_LinearHSPmanager(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [LinearHSPmanager *]
 *
 * Return [UNKN ]  Undocumented return value [LinearHSPmanager *]
 *
 */
LinearHSPmanager * free_LinearHSPmanager(LinearHSPmanager * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a LinearHSPmanager obj. Should be trappable");  
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
    if( obj->set != NULL)    {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->set[i] != NULL) 
          free_HSPset(obj->set[i]);  
        }  
      ckfree(obj->set);  
      }  
    if( obj->query != NULL)  
      free_Sequence(obj->query);     
    if( obj->mat != NULL)    
      free_CompMat(obj->mat);    


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
