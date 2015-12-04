#ifdef _cplusplus
extern "C" {
#endif
#include "seqhit.h"



# line 53 "seqhit.dy"
DnaSequenceHitList * new_DnaSequenceHitList(void)
{
  DnaSequenceHitList * out;

  out = DnaSequenceHitList_alloc();
  out->forward = SegmentHitList_alloc_std();
  out->backward = SegmentHitList_alloc_std();

  return out;
}

/* Function:  show_DnaSequenceHitList(dsl,ofp)
 *
 * Descrip:    shows a DnaSequenceHitsList -
 *
 *             only really useful for debugging
 *
 *
 * Arg:        dsl [UNKN ] Undocumented argument [DnaSequenceHitList *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 69 "seqhit.dy"
void show_DnaSequenceHitList(DnaSequenceHitList * dsl,FILE * ofp)
{
  show_SegmentHitList(dsl->forward,ofp);
  show_SegmentHitList(dsl->backward,ofp);
}

# line 75 "seqhit.dy"
void show_SegmentHitList(SegmentHitList * shl,FILE * ofp)
{
  int i;
  for(i=0;i<shl->len;i++)
    show_SegmentHit(shl->seghit[i],ofp);
  
}

# line 83 "seqhit.dy"
void show_SegmentHit(SegmentHit * sh,FILE * ofp)
{
  fprintf(ofp,"%4d %4d %4d %4d %4.2f %s\n",sh->qstart,sh->qend,sh->tstart,sh->tend,sh->score,sh->name);
}

# line 88 "seqhit.dy"
void sort_SegmentHitList_by_qstart(SegmentHitList * sgl)
{
  sort_SegmentHitList(sgl,compare_SegmentHit_by_qstart);
}

# line 93 "seqhit.dy"
int compare_SegmentHit_by_qstart(SegmentHit * one,SegmentHit * two)
{
  return one->qstart - two->qstart;
}

# line 98 "seqhit.dy"
DnaSequenceHitList * read_msptmp_DnaSequenceHitList(FILE * ifp)
{
  DnaSequenceHitList * out;
  SegmentHit * sh;
  char buffer[MAXLINE];
  char copy[MAXLINE];
  int strand;

  out = new_DnaSequenceHitList();

  while( fgets(buffer,MAXLINE,ifp) != NULL ) {
    strcpy(copy,buffer);
    if( (sh = SegmentHit_from_msptmp_line(buffer,&strand)) == NULL ) {
      warn("Cannot parse [%s] as a MSP crunch line",copy);
      continue;
    }

    if( strand < 0 ) {
      add_SegmentHitList(out->backward,sh);
    } else {
      add_SegmentHitList(out->forward,sh);
    }
  }

  sort_SegmentHitList_by_qstart(out->forward);
  sort_SegmentHitList_by_qstart(out->backward);

  return out;
}

# line 128 "seqhit.dy"
boolean verify_SegmentHit(SegmentHit * sh)
{
  boolean ret = TRUE;

  if ( sh->qstart <= 0 ) {
    ret = FALSE;
    warn("Got a less than zero index in query start %d",sh->qstart);
  }
  if ( sh->qend <= 0 ) {
    ret = FALSE;
    warn("Got a less than zero index in query end %d",sh->qend);
  }
  if ( sh->tstart <= 0 ) {
    ret = FALSE;
    warn("Got a less than zero index in target start %d",sh->tstart);
  }
  if ( sh->tend <= 0 ) {
    ret = FALSE;
    warn("Got a less than zero index in target end %d",sh->tend);
  }

  return ret;
}

/* Function:  read_MSPcrunch_DnaSequenceHitList(ifp)
 *
 * Descrip:    Reads a MSPcrunch -x output file 
 *
 *
 * Arg:        ifp [UNKN ] input file to read [FILE *]
 *
 * Return [UNKN ]  newly allocated structure [DnaSequenceHitList *]
 *
 */
# line 158 "seqhit.dy"
DnaSequenceHitList * read_MSPcrunch_DnaSequenceHitList(FILE * ifp)
{
  DnaSequenceHitList * out;
  SegmentHit * sh;
  char buffer[MAXLINE];
  char copy[MAXLINE];
  int strand;

  out = new_DnaSequenceHitList();

  while( fgets(buffer,MAXLINE,ifp) != NULL ) {
    strcpy(copy,buffer);
    if( (sh = SegmentHit_from_MSPcrunch_line(buffer,&strand)) == NULL ) {
      warn("Cannot parse [%s] as a MSP crunch line",copy);
      continue;
    }
    if( verify_SegmentHit(sh) == FALSE ) {
      warn("Segment hit read in from [%s] was considered bad. Discarding",copy);
      sh = free_SegmentHit(sh);
      continue;
    }

    if( strand < 0 ) {
      add_SegmentHitList(out->backward,sh);
    } else {
      add_SegmentHitList(out->forward,sh);
    }
  }

  sort_SegmentHitList_by_qstart(out->forward);
  sort_SegmentHitList_by_qstart(out->backward);

  /*  fprintf(stderr,"We have %d forward and %d backward guys\n",out->forward->len,out->backward->len); */
  return out;
}



# line 196 "seqhit.dy"
SegmentHit * SegmentHit_from_msptmp_line(char * line,int * strand)
{
  SegmentHit * out;
  char * runner;

  out = SegmentHit_alloc();

  runner = strtok(line,spacestr);

  if( runner == NULL || is_double_string(runner,&out->score) == FALSE ) 
    goto error;

  if( (runner=strtok(NULL,spacestr)) == NULL ) 
    goto error;

  if( *runner != '(' || runner[strlen(runner)-1] != ')' ) 
    goto error;

  runner[strlen(runner)-1] = '\0';
  runner++;
  if( is_integer_string(runner,strand) == FALSE)
    goto error;

  if( (runner=strtok(NULL,spacestr)) == NULL || is_integer_string(runner,&out->qstart) == FALSE ) 
    goto error;

  if( (runner=strtok(NULL,spacestr)) == NULL || is_integer_string(runner,&out->qend) == FALSE ) 
    goto error;

  if( (runner=strtok(NULL,spacestr)) == NULL || is_integer_string(runner,&out->tstart) == FALSE ) 
    goto error;

  if( (runner=strtok(NULL,spacestr)) == NULL || is_integer_string(runner,&out->tend) == FALSE ) 
    goto error;

  if( (runner=strtok(NULL,spacestr)) == NULL ) 
    goto error;

  out->name = stringalloc(runner);

  return out;

  error :

    free_SegmentHit(out);
  return NULL;

}


# line 246 "seqhit.dy"
SegmentHit * SegmentHit_from_MSPcrunch_line(char * line,int * strand)
{
  SegmentHit * out;
  char * runner;

  out = SegmentHit_alloc();

  runner = strtok(line,spacestr);

  if( runner == NULL || is_double_string(runner,&out->score) == FALSE ) 
    goto error;

  if( (runner=strtok(NULL,spacestr)) == NULL ) 
    goto error;

  if( *runner != '(' || runner[strlen(runner)-1] != ')' ) 
    goto error;

  runner[strlen(runner)-1] = '\0';
  runner++;
  if( is_integer_string(runner,strand) == FALSE)
    goto error;

  if( (runner=strtok(NULL,spacestr)) == NULL || is_integer_string(runner,&out->qstart) == FALSE ) 
    goto error;

  if( (runner=strtok(NULL,spacestr)) == NULL || is_integer_string(runner,&out->qend) == FALSE ) 
    goto error;

  if( (runner=strtok(NULL,spacestr)) == NULL || is_integer_string(runner,&out->tstart) == FALSE ) 
    goto error;

  if( (runner=strtok(NULL,spacestr)) == NULL || is_integer_string(runner,&out->tend) == FALSE ) 
    goto error;

  if( (runner=strtok(NULL,spacestr)) == NULL ) 
    goto error;

  out->name = stringalloc(runner);

  return out;

  error :

    free_SegmentHit(out);
  return NULL;

}

  

  


# line 271 "seqhit.c"
/* Function:  hard_link_SegmentHit(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SegmentHit *]
 *
 * Return [UNKN ]  Undocumented return value [SegmentHit *]
 *
 */
SegmentHit * hard_link_SegmentHit(SegmentHit * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a SegmentHit object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  SegmentHit_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SegmentHit *]
 *
 */
SegmentHit * SegmentHit_alloc(void) 
{
    SegmentHit * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(SegmentHit *) ckalloc (sizeof(SegmentHit))) == NULL)    {  
      warn("SegmentHit_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->name = NULL;    
    out->qstart = 0; 
    out->qend = 0;   
    out->tstart = 0; 
    out->tend = 0;   
    out->score = 0;  
    out->next_hit = NULL;    


    return out;  
}    


/* Function:  free_SegmentHit(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SegmentHit *]
 *
 * Return [UNKN ]  Undocumented return value [SegmentHit *]
 *
 */
SegmentHit * free_SegmentHit(SegmentHit * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a SegmentHit obj. Should be trappable");    
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
    if( obj->name != NULL)   
      ckfree(obj->name);     
    if( obj->next_hit != NULL)   
      free_SegmentHit(obj->next_hit);    


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_SequenceHit(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SequenceHit *]
 *
 * Return [UNKN ]  Undocumented return value [SequenceHit *]
 *
 */
SequenceHit * hard_link_SequenceHit(SequenceHit * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a SequenceHit object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  SequenceHit_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SequenceHit *]
 *
 */
SequenceHit * SequenceHit_alloc(void) 
{
    SequenceHit * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(SequenceHit *) ckalloc (sizeof(SequenceHit))) == NULL)  {  
      warn("SequenceHit_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->name = NULL;    
    out->start = NULL;   


    return out;  
}    


/* Function:  free_SequenceHit(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SequenceHit *]
 *
 * Return [UNKN ]  Undocumented return value [SequenceHit *]
 *
 */
SequenceHit * free_SequenceHit(SequenceHit * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a SequenceHit obj. Should be trappable");   
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
    if( obj->name != NULL)   
      ckfree(obj->name);     
    if( obj->start != NULL)  
      free_SegmentHit(obj->start);   


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_SegmentHitList(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_SegmentHitList
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [SegmentHit **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_SegmentHitList(SegmentHit ** list,int i,int j)  
{
    SegmentHit * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_SegmentHitList(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_SegmentHitList which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [SegmentHit **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_SegmentHitList(SegmentHit ** list,int left,int right,int (*comp)(SegmentHit * ,SegmentHit * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_SegmentHitList(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_SegmentHitList (list,++last,i); 
      }  
    swap_SegmentHitList (list,left,last);    
    qsort_SegmentHitList(list,left,last-1,comp); 
    qsort_SegmentHitList(list,last+1,right,comp);    
}    


/* Function:  sort_SegmentHitList(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_SegmentHitList
 *
 *
 * Arg:         obj [UNKN ] Object containing list [SegmentHitList *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_SegmentHitList(SegmentHitList * obj,int (*comp)(SegmentHit *, SegmentHit *)) 
{
    qsort_SegmentHitList(obj->seghit,0,obj->len-1,comp); 
    return;  
}    


/* Function:  expand_SegmentHitList(obj,len)
 *
 * Descrip:    Really an internal function for add_SegmentHitList
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [SegmentHitList *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_SegmentHitList(SegmentHitList * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_SegmentHitList called with no need"); 
      return TRUE;   
      }  


    if( (obj->seghit = (SegmentHit ** ) ckrealloc (obj->seghit,sizeof(SegmentHit *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_SegmentHitList, returning FALSE");   
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_SegmentHitList(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [SegmentHitList *]
 * Arg:        add [OWNER] Object to add to the list [SegmentHit *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_SegmentHitList(SegmentHitList * obj,SegmentHit * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_SegmentHitList(obj,obj->len + SegmentHitListLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->seghit[obj->len++]=add; 
    return TRUE; 
}    


/* Function:  flush_SegmentHitList(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [SegmentHitList *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_SegmentHitList(SegmentHitList * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->seghit[i] != NULL)    {  
        free_SegmentHit(obj->seghit[i]); 
        obj->seghit[i] = NULL;   
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  SegmentHitList_alloc_std(void)
 *
 * Descrip:    Equivalent to SegmentHitList_alloc_len(SegmentHitListLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SegmentHitList *]
 *
 */
SegmentHitList * SegmentHitList_alloc_std(void) 
{
    return SegmentHitList_alloc_len(SegmentHitListLISTLENGTH);   
}    


/* Function:  SegmentHitList_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [SegmentHitList *]
 *
 */
SegmentHitList * SegmentHitList_alloc_len(int len) 
{
    SegmentHitList * out;   /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = SegmentHitList_alloc()) == NULL)   
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->seghit = (SegmentHit ** ) ckcalloc (len,sizeof(SegmentHit *))) == NULL) {  
      warn("Warning, ckcalloc failed in SegmentHitList_alloc_len");  
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_SegmentHitList(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SegmentHitList *]
 *
 * Return [UNKN ]  Undocumented return value [SegmentHitList *]
 *
 */
SegmentHitList * hard_link_SegmentHitList(SegmentHitList * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a SegmentHitList object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  SegmentHitList_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SegmentHitList *]
 *
 */
SegmentHitList * SegmentHitList_alloc(void) 
{
    SegmentHitList * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(SegmentHitList *) ckalloc (sizeof(SegmentHitList))) == NULL)    {  
      warn("SegmentHitList_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->seghit = NULL;  
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_SegmentHitList(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SegmentHitList *]
 *
 * Return [UNKN ]  Undocumented return value [SegmentHitList *]
 *
 */
SegmentHitList * free_SegmentHitList(SegmentHitList * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a SegmentHitList obj. Should be trappable");    
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
    if( obj->seghit != NULL) {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->seghit[i] != NULL)  
          free_SegmentHit(obj->seghit[i]);   
        }  
      ckfree(obj->seghit);   
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_DnaSequenceHitList(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DnaSequenceHitList *]
 *
 * Return [UNKN ]  Undocumented return value [DnaSequenceHitList *]
 *
 */
DnaSequenceHitList * hard_link_DnaSequenceHitList(DnaSequenceHitList * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a DnaSequenceHitList object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  DnaSequenceHitList_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaSequenceHitList *]
 *
 */
DnaSequenceHitList * DnaSequenceHitList_alloc(void) 
{
    DnaSequenceHitList * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(DnaSequenceHitList *) ckalloc (sizeof(DnaSequenceHitList))) == NULL)    {  
      warn("DnaSequenceHitList_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->forward = NULL; 
    out->backward = NULL;    


    return out;  
}    


/* Function:  free_DnaSequenceHitList(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DnaSequenceHitList *]
 *
 * Return [UNKN ]  Undocumented return value [DnaSequenceHitList *]
 *
 */
DnaSequenceHitList * free_DnaSequenceHitList(DnaSequenceHitList * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a DnaSequenceHitList obj. Should be trappable");    
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
    if( obj->forward != NULL)    
      free_SegmentHitList(obj->forward);     
    if( obj->backward != NULL)   
      free_SegmentHitList(obj->backward);    


    ckfree(obj); 
    return NULL; 
}    


/* Function:  replace_forward_DnaSequenceHitList(obj,forward)
 *
 * Descrip:    Replace member variable forward
 *             For use principly by API functions
 *
 *
 * Arg:            obj [UNKN ] Object holding the variable [DnaSequenceHitList *]
 * Arg:        forward [OWNER] New value of the variable [SegmentHitList *]
 *
 * Return [SOFT ]  member variable forward [boolean]
 *
 */
boolean replace_forward_DnaSequenceHitList(DnaSequenceHitList * obj,SegmentHitList * forward) 
{
    if( obj == NULL)     {  
      warn("In replacement function forward for object DnaSequenceHitList, got a NULL object");  
      return FALSE;  
      }  
    obj->forward = forward;  
    return TRUE; 
}    


/* Function:  access_forward_DnaSequenceHitList(obj)
 *
 * Descrip:    Access member variable forward
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [DnaSequenceHitList *]
 *
 * Return [SOFT ]  member variable forward [SegmentHitList *]
 *
 */
SegmentHitList * access_forward_DnaSequenceHitList(DnaSequenceHitList * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function forward for object DnaSequenceHitList, got a NULL object"); 
      return NULL;   
      }  
    return obj->forward;     
}    


/* Function:  replace_backward_DnaSequenceHitList(obj,backward)
 *
 * Descrip:    Replace member variable backward
 *             For use principly by API functions
 *
 *
 * Arg:             obj [UNKN ] Object holding the variable [DnaSequenceHitList *]
 * Arg:        backward [OWNER] New value of the variable [SegmentHitList *]
 *
 * Return [SOFT ]  member variable backward [boolean]
 *
 */
boolean replace_backward_DnaSequenceHitList(DnaSequenceHitList * obj,SegmentHitList * backward) 
{
    if( obj == NULL)     {  
      warn("In replacement function backward for object DnaSequenceHitList, got a NULL object"); 
      return FALSE;  
      }  
    obj->backward = backward;    
    return TRUE; 
}    


/* Function:  access_backward_DnaSequenceHitList(obj)
 *
 * Descrip:    Access member variable backward
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [DnaSequenceHitList *]
 *
 * Return [SOFT ]  member variable backward [SegmentHitList *]
 *
 */
SegmentHitList * access_backward_DnaSequenceHitList(DnaSequenceHitList * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function backward for object DnaSequenceHitList, got a NULL object");    
      return NULL;   
      }  
    return obj->backward;    
}    


/* Function:  access_seghit_SegmentHitList(obj,i)
 *
 * Descrip:    Access members stored in the seghit list
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the list [SegmentHitList *]
 * Arg:          i [UNKN ] Position in the list [int]
 *
 * Return [SOFT ]  Element of the list [SegmentHit *]
 *
 */
SegmentHit * access_seghit_SegmentHitList(SegmentHitList * obj,int i) 
{
    if( obj == NULL)     {  
      warn("In accessor function seghit for object SegmentHitList, got a NULL object");  
      return NULL;   
      }  
    if( obj->len <= i )  {  
      warn("In accessor function seghit for object SegmentHitList, index %%d is greater than list length %%d",i,obj->len);   
      return NULL;   
      }  
    return obj->seghit[i];   
}    


/* Function:  length_seghit_SegmentHitList(obj)
 *
 * Descrip:    discover the length of the list
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the list [SegmentHitList *]
 *
 * Return [UNKN ]  length of the list [int]
 *
 */
int length_seghit_SegmentHitList(SegmentHitList * obj) 
{
    if( obj == NULL)     {  
      warn("In length function seghit for object SegmentHitList, got a NULL object");    
      return -1;     
      }  
    return obj->len;     
}    


/* Function:  replace_name_SegmentHit(obj,name)
 *
 * Descrip:    Replace member variable name
 *             For use principly by API functions
 *
 *
 * Arg:         obj [UNKN ] Object holding the variable [SegmentHit *]
 * Arg:        name [OWNER] New value of the variable [char *]
 *
 * Return [SOFT ]  member variable name [boolean]
 *
 */
boolean replace_name_SegmentHit(SegmentHit * obj,char * name) 
{
    if( obj == NULL)     {  
      warn("In replacement function name for object SegmentHit, got a NULL object"); 
      return FALSE;  
      }  
    obj->name = name;    
    return TRUE; 
}    


/* Function:  access_name_SegmentHit(obj)
 *
 * Descrip:    Access member variable name
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [SegmentHit *]
 *
 * Return [SOFT ]  member variable name [char *]
 *
 */
char * access_name_SegmentHit(SegmentHit * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function name for object SegmentHit, got a NULL object");    
      return NULL;   
      }  
    return obj->name;    
}    


/* Function:  replace_qstart_SegmentHit(obj,qstart)
 *
 * Descrip:    Replace member variable qstart
 *             For use principly by API functions
 *
 *
 * Arg:           obj [UNKN ] Object holding the variable [SegmentHit *]
 * Arg:        qstart [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable qstart [boolean]
 *
 */
boolean replace_qstart_SegmentHit(SegmentHit * obj,int qstart) 
{
    if( obj == NULL)     {  
      warn("In replacement function qstart for object SegmentHit, got a NULL object");   
      return FALSE;  
      }  
    obj->qstart = qstart;    
    return TRUE; 
}    


/* Function:  access_qstart_SegmentHit(obj)
 *
 * Descrip:    Access member variable qstart
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [SegmentHit *]
 *
 * Return [SOFT ]  member variable qstart [int]
 *
 */
int access_qstart_SegmentHit(SegmentHit * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function qstart for object SegmentHit, got a NULL object");  
      return 0;  
      }  
    return obj->qstart;  
}    


/* Function:  replace_qend_SegmentHit(obj,qend)
 *
 * Descrip:    Replace member variable qend
 *             For use principly by API functions
 *
 *
 * Arg:         obj [UNKN ] Object holding the variable [SegmentHit *]
 * Arg:        qend [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable qend [boolean]
 *
 */
boolean replace_qend_SegmentHit(SegmentHit * obj,int qend) 
{
    if( obj == NULL)     {  
      warn("In replacement function qend for object SegmentHit, got a NULL object"); 
      return FALSE;  
      }  
    obj->qend = qend;    
    return TRUE; 
}    


/* Function:  access_qend_SegmentHit(obj)
 *
 * Descrip:    Access member variable qend
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [SegmentHit *]
 *
 * Return [SOFT ]  member variable qend [int]
 *
 */
int access_qend_SegmentHit(SegmentHit * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function qend for object SegmentHit, got a NULL object");    
      return 0;  
      }  
    return obj->qend;    
}    


/* Function:  replace_tstart_SegmentHit(obj,tstart)
 *
 * Descrip:    Replace member variable tstart
 *             For use principly by API functions
 *
 *
 * Arg:           obj [UNKN ] Object holding the variable [SegmentHit *]
 * Arg:        tstart [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable tstart [boolean]
 *
 */
boolean replace_tstart_SegmentHit(SegmentHit * obj,int tstart) 
{
    if( obj == NULL)     {  
      warn("In replacement function tstart for object SegmentHit, got a NULL object");   
      return FALSE;  
      }  
    obj->tstart = tstart;    
    return TRUE; 
}    


/* Function:  access_tstart_SegmentHit(obj)
 *
 * Descrip:    Access member variable tstart
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [SegmentHit *]
 *
 * Return [SOFT ]  member variable tstart [int]
 *
 */
int access_tstart_SegmentHit(SegmentHit * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function tstart for object SegmentHit, got a NULL object");  
      return 0;  
      }  
    return obj->tstart;  
}    


/* Function:  replace_tend_SegmentHit(obj,tend)
 *
 * Descrip:    Replace member variable tend
 *             For use principly by API functions
 *
 *
 * Arg:         obj [UNKN ] Object holding the variable [SegmentHit *]
 * Arg:        tend [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable tend [boolean]
 *
 */
boolean replace_tend_SegmentHit(SegmentHit * obj,int tend) 
{
    if( obj == NULL)     {  
      warn("In replacement function tend for object SegmentHit, got a NULL object"); 
      return FALSE;  
      }  
    obj->tend = tend;    
    return TRUE; 
}    


/* Function:  access_tend_SegmentHit(obj)
 *
 * Descrip:    Access member variable tend
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [SegmentHit *]
 *
 * Return [SOFT ]  member variable tend [int]
 *
 */
int access_tend_SegmentHit(SegmentHit * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function tend for object SegmentHit, got a NULL object");    
      return 0;  
      }  
    return obj->tend;    
}    


/* Function:  replace_score_SegmentHit(obj,score)
 *
 * Descrip:    Replace member variable score
 *             For use principly by API functions
 *
 *
 * Arg:          obj [UNKN ] Object holding the variable [SegmentHit *]
 * Arg:        score [OWNER] New value of the variable [double]
 *
 * Return [SOFT ]  member variable score [boolean]
 *
 */
boolean replace_score_SegmentHit(SegmentHit * obj,double score) 
{
    if( obj == NULL)     {  
      warn("In replacement function score for object SegmentHit, got a NULL object");    
      return FALSE;  
      }  
    obj->score = score;  
    return TRUE; 
}    


/* Function:  access_score_SegmentHit(obj)
 *
 * Descrip:    Access member variable score
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [SegmentHit *]
 *
 * Return [SOFT ]  member variable score [double]
 *
 */
double access_score_SegmentHit(SegmentHit * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function score for object SegmentHit, got a NULL object");   
      return 0;  
      }  
    return obj->score;   
}    


/* Function:  replace_next_hit_SegmentHit(obj,next_hit)
 *
 * Descrip:    Replace member variable next_hit
 *             For use principly by API functions
 *
 *
 * Arg:             obj [UNKN ] Object holding the variable [SegmentHit *]
 * Arg:        next_hit [OWNER] New value of the variable [SegmentHit *]
 *
 * Return [SOFT ]  member variable next_hit [boolean]
 *
 */
boolean replace_next_hit_SegmentHit(SegmentHit * obj,SegmentHit * next_hit) 
{
    if( obj == NULL)     {  
      warn("In replacement function next_hit for object SegmentHit, got a NULL object"); 
      return FALSE;  
      }  
    obj->next_hit = next_hit;    
    return TRUE; 
}    


/* Function:  access_next_hit_SegmentHit(obj)
 *
 * Descrip:    Access member variable next_hit
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [SegmentHit *]
 *
 * Return [SOFT ]  member variable next_hit [SegmentHit *]
 *
 */
SegmentHit * access_next_hit_SegmentHit(SegmentHit * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function next_hit for object SegmentHit, got a NULL object");    
      return NULL;   
      }  
    return obj->next_hit;    
}    



#ifdef _cplusplus
}
#endif
