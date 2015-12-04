#ifdef _cplusplus
extern "C" {
#endif
#include "seqerror.h"


/* Function:  genewise_SequenceErrorSet(als)
 *
 * Descrip:    Makes a sequence error set from standard genewise labels
 *
 *
 * Arg:        als [UNKN ] Undocumented argument [AlnSequence *]
 *
 * Return [UNKN ]  Undocumented return value [SequenceErrorSet *]
 *
 */
# line 56 "seqerror.dy"
SequenceErrorSet * genewise_SequenceErrorSet(AlnSequence * als )
{
  SequenceErrorSet * out;
  SequenceError * se;
  AlnUnit * alu;

  assert(als);
  out = SequenceErrorSet_alloc_std();
  
  for(alu=als->start;alu;alu = alu->next ) {
    if( strcmp(alu->text_label,"SEQUENCE_INSERTION") == 0 ) {
      se = SequenceError_alloc();
      add_SequenceErrorSet(out,se);
      se->start = alu->start+1;
      se->end   = alu->end;
      
      se->type = SeqErrorInsertion;
      se->inserted_bases = stringalloc("?");
      se->replaced_bases = stringalloc("NNN");
    } else if ( strcmp(alu->text_label,"SEQUENCE_DELETION") == 0 ) {
      se = SequenceError_alloc();
      se->start = alu->start+1;
      se->end   = alu->end;
      se->type = SeqErrorDeletion;
      se->replaced_bases = stringalloc("NNN");
      add_SequenceErrorSet(out,se);
    }
  }

  return out;
}



/* Function:  show_SequenceErrorSet(ses,ofp)
 *
 * Descrip:    Displays a set of sequence errors in space deliminted format
 *
 *
 * Arg:        ses [UNKN ] Undocumented argument [SequenceErrorSet *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 93 "seqerror.dy"
void show_SequenceErrorSet(SequenceErrorSet * ses,FILE * ofp)
{
  int i;

  for(i=0;i<ses->len;i++) {
    show_SequenceError(ses->error[i],ofp);
  }
}

    

/* Function:  show_SequenceError(se,ofp)
 *
 * Descrip:    Displays sequence error in space deliminted format
 *
 *
 * Arg:         se [UNKN ] Undocumented argument [SequenceError *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 107 "seqerror.dy"
void show_SequenceError(SequenceError * se,FILE * ofp)
{
  fprintf(ofp,"%c %5d %5d %10s (%10s)\n",
	  (se->type == SeqErrorInsertion ? 'I' : 
	   (se->type == SeqErrorDeletion ? 'D' : 'S')),
	  se->start,se->end,
	  se->replaced_bases == NULL ? "-none-" : se->replaced_bases,
	  se->type == SeqErrorInsertion ? se->inserted_bases : "-");
}

/* Function:  make_ErrorSequence(seq,subs,insertion,deletion)
 *
 * Descrip:    Makes an error sequence (DNA) with set substitution
 *             and insertion/deletion rates.
 *
 *
 * Arg:              seq [UNKN ] Undocumented argument [Sequence *]
 * Arg:             subs [UNKN ] Undocumented argument [Probability]
 * Arg:        insertion [UNKN ] Undocumented argument [Probability]
 * Arg:         deletion [UNKN ] Undocumented argument [Probability]
 *
 * Return [UNKN ]  Undocumented return value [ErrorSequence *]
 *
 */
# line 121 "seqerror.dy"
ErrorSequence * make_ErrorSequence(Sequence * seq,Probability subs,Probability insertion,Probability deletion)
{
  ErrorSequence * out;
  
  char * seq_buffer;
  int buf_len;
  int i,j;
  double rnd;
  char c[2];
  SequenceError * se;
  
  i=j=0;

  out = ErrorSequence_alloc();
  out->ses = SequenceErrorSet_alloc_std();
  c[1] = '\0';

  buf_len = seq->len *3;
  seq_buffer = (char*)calloc(buf_len,sizeof(char));
  
  for(i=0;i<seq->len;i++) {
    if( j >= buf_len ) {
      warn("run out of temporary buffer in error sequence");
      break;
    }

    rnd = random_0_to_1();

    if( rnd > deletion ) {
      /* yes this base is here */
      rnd = random_0_to_1();
      if( rnd < subs ) {
	rnd = random_0_to_1();
	if( rnd < 0.25 ) {
	  c[0] = 'A';
	} else if ( rnd < 0.5 ) {
	  c[0] = 'T';
	} else if ( rnd < 0.75 ) {
	  c[0] = 'G';
	} else {
	  c[0] = 'C';
	}

	se = SequenceError_alloc();
	se->type = SeqErrorSubstitution;
	se->start = j+1;
	se->end = j+1;
	se->replaced_bases = stringalloc(c);
	add_SequenceErrorSet(out->ses,se);
	seq_buffer[j++] = c[0];
      } else {
	seq_buffer[j++] = seq->seq[i];
      }
    } else {
      /* there has been a deletion */
      se = SequenceError_alloc();
      se->type = SeqErrorDeletion;
      se->start = j+1;
      se->end = j+1;
      c[0] = seq_buffer[i];
      se->suspected_deletion = stringalloc(c);
      add_SequenceErrorSet(out->ses,se);
    }

    rnd = random_0_to_1();
    
    if( rnd < insertion ) {
   
      rnd = random_0_to_1();
      if( rnd < 0.25 ) {
	c[0] = 'A';
      } else if ( rnd < 0.5 ) {
	c[0] = 'T';
      } else if ( rnd < 0.75 ) {
	c[0] = 'G';
      } else {
	c[0] = 'C';
      }

      se = SequenceError_alloc();
      se->type = SeqErrorInsertion;
      se->start = j+1;
      se->end = j+1;
      se->inserted_bases = stringalloc(c);
      add_SequenceErrorSet(out->ses,se);
      seq_buffer[j++] = c[0];
    }
  }

  seq_buffer[j] = '\0';

  out->seq = new_Sequence_from_strings(seq->name,seq_buffer);

  return out;
}
    







# line 202 "seqerror.c"
/* Function:  hard_link_SequenceError(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SequenceError *]
 *
 * Return [UNKN ]  Undocumented return value [SequenceError *]
 *
 */
SequenceError * hard_link_SequenceError(SequenceError * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a SequenceError object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  SequenceError_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SequenceError *]
 *
 */
SequenceError * SequenceError_alloc(void) 
{
    SequenceError * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(SequenceError *) ckalloc (sizeof(SequenceError))) == NULL)  {  
      warn("SequenceError_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->type = 0;   
    out->start = 0;  
    out->end = 0;    
    out->replaced_bases = NULL;  
    out->probability = 0.0;  
    out->inserted_bases = NULL;  
    out->suspected_deletion = NULL;  


    return out;  
}    


/* Function:  free_SequenceError(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SequenceError *]
 *
 * Return [UNKN ]  Undocumented return value [SequenceError *]
 *
 */
SequenceError * free_SequenceError(SequenceError * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a SequenceError obj. Should be trappable"); 
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
    if( obj->replaced_bases != NULL) 
      ckfree(obj->replaced_bases);   
    if( obj->inserted_bases != NULL) 
      ckfree(obj->inserted_bases);   
    if( obj->suspected_deletion != NULL) 
      ckfree(obj->suspected_deletion);   


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_SequenceErrorSet(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_SequenceErrorSet
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [SequenceError **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_SequenceErrorSet(SequenceError ** list,int i,int j)  
{
    SequenceError * temp;    
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_SequenceErrorSet(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_SequenceErrorSet which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [SequenceError **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_SequenceErrorSet(SequenceError ** list,int left,int right,int (*comp)(SequenceError * ,SequenceError * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_SequenceErrorSet(list,left,(left+right)/2); 
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_SequenceErrorSet (list,++last,i);   
      }  
    swap_SequenceErrorSet (list,left,last);  
    qsort_SequenceErrorSet(list,left,last-1,comp);   
    qsort_SequenceErrorSet(list,last+1,right,comp);  
}    


/* Function:  sort_SequenceErrorSet(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_SequenceErrorSet
 *
 *
 * Arg:         obj [UNKN ] Object containing list [SequenceErrorSet *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_SequenceErrorSet(SequenceErrorSet * obj,int (*comp)(SequenceError *, SequenceError *)) 
{
    qsort_SequenceErrorSet(obj->error,0,obj->len-1,comp);    
    return;  
}    


/* Function:  expand_SequenceErrorSet(obj,len)
 *
 * Descrip:    Really an internal function for add_SequenceErrorSet
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [SequenceErrorSet *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_SequenceErrorSet(SequenceErrorSet * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_SequenceErrorSet called with no need");   
      return TRUE;   
      }  


    if( (obj->error = (SequenceError ** ) ckrealloc (obj->error,sizeof(SequenceError *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_SequenceErrorSet, returning FALSE"); 
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_SequenceErrorSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [SequenceErrorSet *]
 * Arg:        add [OWNER] Object to add to the list [SequenceError *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_SequenceErrorSet(SequenceErrorSet * obj,SequenceError * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_SequenceErrorSet(obj,obj->len + SequenceErrorSetLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->error[obj->len++]=add;  
    return TRUE; 
}    


/* Function:  flush_SequenceErrorSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [SequenceErrorSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_SequenceErrorSet(SequenceErrorSet * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->error[i] != NULL) {  
        free_SequenceError(obj->error[i]);   
        obj->error[i] = NULL;    
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  SequenceErrorSet_alloc_std(void)
 *
 * Descrip:    Equivalent to SequenceErrorSet_alloc_len(SequenceErrorSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SequenceErrorSet *]
 *
 */
SequenceErrorSet * SequenceErrorSet_alloc_std(void) 
{
    return SequenceErrorSet_alloc_len(SequenceErrorSetLISTLENGTH);   
}    


/* Function:  SequenceErrorSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [SequenceErrorSet *]
 *
 */
SequenceErrorSet * SequenceErrorSet_alloc_len(int len) 
{
    SequenceErrorSet * out; /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = SequenceErrorSet_alloc()) == NULL) 
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->error = (SequenceError ** ) ckcalloc (len,sizeof(SequenceError *))) == NULL)    {  
      warn("Warning, ckcalloc failed in SequenceErrorSet_alloc_len");    
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_SequenceErrorSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SequenceErrorSet *]
 *
 * Return [UNKN ]  Undocumented return value [SequenceErrorSet *]
 *
 */
SequenceErrorSet * hard_link_SequenceErrorSet(SequenceErrorSet * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a SequenceErrorSet object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  SequenceErrorSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SequenceErrorSet *]
 *
 */
SequenceErrorSet * SequenceErrorSet_alloc(void) 
{
    SequenceErrorSet * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(SequenceErrorSet *) ckalloc (sizeof(SequenceErrorSet))) == NULL)    {  
      warn("SequenceErrorSet_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->error = NULL;   
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_SequenceErrorSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SequenceErrorSet *]
 *
 * Return [UNKN ]  Undocumented return value [SequenceErrorSet *]
 *
 */
SequenceErrorSet * free_SequenceErrorSet(SequenceErrorSet * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a SequenceErrorSet obj. Should be trappable");  
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
    if( obj->error != NULL)  {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->error[i] != NULL)   
          free_SequenceError(obj->error[i]); 
        }  
      ckfree(obj->error);    
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_ErrorSequence(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ErrorSequence *]
 *
 * Return [UNKN ]  Undocumented return value [ErrorSequence *]
 *
 */
ErrorSequence * hard_link_ErrorSequence(ErrorSequence * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a ErrorSequence object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  ErrorSequence_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ErrorSequence *]
 *
 */
ErrorSequence * ErrorSequence_alloc(void) 
{
    ErrorSequence * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(ErrorSequence *) ckalloc (sizeof(ErrorSequence))) == NULL)  {  
      warn("ErrorSequence_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->seq = NULL; 
    out->ses = NULL; 


    return out;  
}    


/* Function:  free_ErrorSequence(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ErrorSequence *]
 *
 * Return [UNKN ]  Undocumented return value [ErrorSequence *]
 *
 */
ErrorSequence * free_ErrorSequence(ErrorSequence * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a ErrorSequence obj. Should be trappable"); 
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
    if( obj->ses != NULL)    
      free_SequenceErrorSet(obj->ses);   


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
