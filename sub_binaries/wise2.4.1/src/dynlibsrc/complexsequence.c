#ifdef _cplusplus
extern "C" {
#endif
#include "complexsequence.h"


/* Function:  reversed_ComplexSequence(forward,cses)
 *
 * Descrip: No Description
 *
 * Arg:        forward [UNKN ] Undocumented argument [Sequence *]
 * Arg:           cses [UNKN ] Undocumented argument [ComplexSequenceEvalSet *]
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequence *]
 *
 */
# line 102 "complexsequence.dy"
ComplexSequence * reversed_ComplexSequence(Sequence * forward,ComplexSequenceEvalSet * cses) 
{
  ComplexSequence * out;
  register int i;
  register int j;
  register int * poi;
  Sequence * rev;
  int seqlen;

  if( forward == NULL ) {
    warn("Trying to make a complex sequence from a NULL sequence - impossible!");
    return NULL;
  }

  if( can_evaluate_this_Sequence(cses,forward) == FALSE ) {
    warn("Could not evaluate these sequences Sequence type [%d][%s] Evaluation type [%d][%s]",forward->type,Sequence_type_to_string(forward->type),cses->type,Sequence_type_to_string(cses->type));
    return NULL;
  }

  if( cses->has_been_prepared == FALSE) {
    warn("Trappable error: you have not prepared this ComplexSequenceEvalSet before using. Please do so in the future");
    prepare_ComplexSequenceEvalSet(cses);
  }

  rev = reverse_complement_Sequence(forward);

  out = ComplexSequence_alloc();
  if( out == NULL )
    return NULL;

  out->creator = hard_link_ComplexSequenceEvalSet(cses);

  out->datastore = (int *) ckcalloc((forward->len+cses->left_lookback)*cses->len,sizeof(int));

  if( out->datastore == NULL ) {
    warn("Could not allocate data pointer of length %d for ComplexSequence",forward->len*cses->len);
    free_ComplexSequence(out);
    return NULL;
  }

  out->data  = out->datastore + (cses->left_lookback * cses->len);

  for(i=0;i<cses->left_lookback;i++) {
    for(j=0;j<cses->len;j++) 
      out->datastore[(i*cses->len)+j] = cses->cse[j]->outside_score;
      }

  out->depth = cses->len;
  seqlen = forward->len;

  for(i=0,poi = out->data;i<forward->len;i++,poi = next_ComplexSequence_data(out,poi)) {
    for(j=0;j<cses->len;j++) {
      /*      fprintf(stderr,"Calling with i at %d vs %d\n",i,cses->cse[j]->left_window); */
      if( i < cses->cse[j]->right_window || (i + cses->cse[j]->left_window) >= seqlen)
	poi[j] = cses->cse[j]->outside_score;
      else poi[j] = (*cses->cse[j]->eval_func)(cses->cse[j]->data_type,cses->cse[j]->data,rev->seq+seqlen-i);
    }
  }

  
  return out;

}

/* Function:  new_ComplexSequence(seq,cses)
 *
 * Descrip:    The basic way to make a ComplexSequence. Requires that
 *             you have already built a ComplexSequenceEvalSet (such as
 *             /default_aminoacid_ComplexSequenceEvalSet).
 *
 *
 *
 * Arg:         seq [UNKN ] Sequence that the ComplexSequence is based on [Sequence *]
 * Arg:        cses [UNKN ] EvalSet that defines the functions used on the sequence [ComplexSequenceEvalSet *]
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequence *]
 *
 */
# line 175 "complexsequence.dy"
ComplexSequence * new_ComplexSequence(Sequence * seq,ComplexSequenceEvalSet * cses)
{
  ComplexSequence * out;
  register int i;
  register int j;
  register int * poi;
  
  if( seq == NULL ) {
    warn("Trying to make a complex sequence from a NULL sequence - impossible!");
    return NULL;
  }

  if( can_evaluate_this_Sequence(cses,seq) == FALSE ) {
    warn("Could not evaluate these sequences Sequence type [%d][%s] Evaluation type [%d][%s]",seq->type,Sequence_type_to_string(seq->type),cses->type,Sequence_type_to_string(cses->type));
    return NULL;
  }

  if( cses->has_been_prepared == FALSE) {
    warn("Trappable error: you have not prepared this ComplexSequenceEvalSet before using. Please do so in the future");
    prepare_ComplexSequenceEvalSet(cses);
  }


  out = ComplexSequence_alloc();
  if( out == NULL )
    return NULL;

  out->creator = hard_link_ComplexSequenceEvalSet(cses);

  out->datastore = (int *) ckcalloc((seq->len+cses->left_lookback)*cses->len,sizeof(int));

  if( out->datastore == NULL ) {
    warn("Could not allocate data pointer of length %d for ComplexSequence",seq->len*cses->len);
    free_ComplexSequence(out);
    return NULL;
  }

  out->data  = out->datastore + (cses->left_lookback * cses->len);

  for(i=0;i<cses->left_lookback;i++) {
    for(j=0;j<cses->len;j++) 
      out->datastore[(i*cses->len)+j] = cses->cse[j]->outside_score;
  }

  out->depth = cses->len;

  for(i=0,poi = out->data;i<seq->len;i++,poi = next_ComplexSequence_data(out,poi)) {
    for(j=0;j<cses->len;j++) {
      /*      fprintf(stderr,"Calling with i at %d vs %d\n",i,cses->cse[j]->left_window); */
      if( i < cses->cse[j]->left_window || (i + cses->cse[j]->right_window) >= seq->len)
	poi[j] = cses->cse[j]->outside_score;
      else poi[j] = (*cses->cse[j]->eval_func)(cses->cse[j]->data_type,cses->cse[j]->data,seq->seq+i);
    }
  }
  
  out->seq    = hard_link_Sequence(seq);

  out->length = seq->len;

  return out;

}

/* Function:  show_ComplexSequence(cs,ofp)
 *
 * Descrip:    shows complex sequence in a vaguely
 *             human form
 *
 *
 * Arg:         cs [UNKN ] Undocumented argument [ComplexSequence *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 243 "complexsequence.dy"
void show_ComplexSequence(ComplexSequence * cs,FILE * ofp)
{
  register int i;
	
  assert(cs);
  assert(ofp);
  assert(cs->seq);
  fprintf(ofp,"ComplexSequence %s\n",cs->seq->name);
  for(i=0;i<cs->length;i++) {
    show_one_position_ComplexSequence(cs,i,ofp);
  }

}

/* Function:  show_one_position_ComplexSequence(cs,pos,ofp)
 *
 * Descrip:    shows one position of a complex sequence
 *
 *
 * Arg:         cs [UNKN ] Undocumented argument [ComplexSequence *]
 * Arg:        pos [UNKN ] Undocumented argument [int]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 262 "complexsequence.dy"
void show_one_position_ComplexSequence(ComplexSequence * cs,int pos,FILE * ofp)
{
  int i;

  assert(cs);
  assert(ofp);

  fprintf(ofp,"%4d  %c [",pos,cs->seq->seq[pos]);

  for(i=0;i<cs->depth;i++) {
    if( cs->creator != NULL ) {
      if( cs->creator->cse[i]->score_type == CseScoreType_Index ) {
	fprintf(ofp,"%4d%c",ComplexSequence_data(cs,pos,i),i == cs->depth-1 ? ']' : ',');
      } else {
	fprintf(ofp,"% 6d (%0+2.2f)%c",ComplexSequence_data(cs,pos,i),Score2Bits(ComplexSequence_data(cs,pos,i)), i == cs->depth-1 ? ']' : ',');
      }
    } else {
      fprintf(ofp,"%4d%c",ComplexSequence_data(cs,pos,i),i == cs->depth-1 ? ']' : ',');
    }
  }
  fprintf(ofp,"\n");
}


  /********************************/
  /* Making complex sequences     */
  /********************************/

/* Function:  prepare_ComplexSequenceEvalSet(cses)
 *
 * Descrip:    Calculates all the necessary offset for an EvalSet.
 *             This is necessary before using it in a /new_ComplexSequence
 *             place
 *
 *
 * Arg:        cses [UNKN ] Undocumented argument [ComplexSequenceEvalSet *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 296 "complexsequence.dy"
boolean prepare_ComplexSequenceEvalSet(ComplexSequenceEvalSet * cses)
{
  register int i;
  int left_window = 0;
  int right_window = 0;
  int left_lookback = 0;

  for(i=0;i<cses->len;i++) {
    if( cses->cse[i]->right_window > right_window )
      right_window = cses->cse[i]->right_window;
    if( cses->cse[i]->left_window > left_window ) 
      left_window = cses->cse[i]->left_window;
    if( cses->cse[i]->left_lookback > left_lookback ) 
      left_lookback = cses->cse[i]->left_lookback;
  }

  cses->right_window  = right_window;
  cses->left_window   = left_window;
  cses->left_lookback = left_lookback;
  cses->has_been_prepared = TRUE;

  return TRUE;
}


/* Function:  can_evaluate_this_Sequence(cses,s)
 *
 * Descrip:    Checks that this ComplexSequenceEvalSet can be used with
 *             this Sequence. This is probably going to go defunct.
 *
 *
 *
 * Arg:        cses [UNKN ] Undocumented argument [ComplexSequenceEvalSet *]
 * Arg:           s [UNKN ] Undocumented argument [Sequence *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 327 "complexsequence.dy"
boolean can_evaluate_this_Sequence(ComplexSequenceEvalSet * cses,Sequence * s)
{
  return can_evaluate_this_type(cses,s->type);
}


/* Function:  can_evaluate_this_type(cses,type)
 *
 * Descrip:    Pretty much an internal for /can_evaluate_this_Sequence
 *
 *
 * Arg:        cses [UNKN ] Undocumented argument [ComplexSequenceEvalSet *]
 * Arg:        type [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 338 "complexsequence.dy"
boolean can_evaluate_this_type(ComplexSequenceEvalSet * cses,int type)
{
  switch (type) {
  case SEQUENCE_UNKNOWN : return FALSE;
  case SEQUENCE_DNA :
    if( cses->type == SEQUENCE_DNA || cses->type == SEQUENCE_CDNA || cses->type == SEQUENCE_EST || cses->type == SEQUENCE_GENOMIC)
      return TRUE;
    else return FALSE;
  default :
    if( cses->type == type) 
      return TRUE;
    else return FALSE;
  }
}


    

# line 297 "complexsequence.c"
/* Function:  hard_link_ComplexSequenceEval(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ComplexSequenceEval *]
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequenceEval *]
 *
 */
ComplexSequenceEval * hard_link_ComplexSequenceEval(ComplexSequenceEval * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a ComplexSequenceEval object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  ComplexSequenceEval_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequenceEval *]
 *
 */
ComplexSequenceEval * ComplexSequenceEval_alloc(void) 
{
    ComplexSequenceEval * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(ComplexSequenceEval *) ckalloc (sizeof(ComplexSequenceEval))) == NULL)  {  
      warn("ComplexSequenceEval_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->type = 0;   
    out->sequence_type = 0;  
    out->left_window = 0;    
    out->right_window = 0;   
    out->left_lookback = 0;  
    out->outside_score = 0;  
    out->data_type = 0;  
    out->eval_func = NULL;   
    out->score_type = CseScoreType_Index;    


    return out;  
}    


/* Function:  free_ComplexSequenceEval(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ComplexSequenceEval *]
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequenceEval *]
 *
 */
ComplexSequenceEval * free_ComplexSequenceEval(ComplexSequenceEval * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a ComplexSequenceEval obj. Should be trappable");   
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
    /* obj->eval_func is a function pointer */ 


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_ComplexSequenceEvalSet(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_ComplexSequenceEvalSet
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [ComplexSequenceEval **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_ComplexSequenceEvalSet(ComplexSequenceEval ** list,int i,int j)  
{
    ComplexSequenceEval * temp;  
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_ComplexSequenceEvalSet(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_ComplexSequenceEvalSet which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [ComplexSequenceEval **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_ComplexSequenceEvalSet(ComplexSequenceEval ** list,int left,int right,int (*comp)(ComplexSequenceEval * ,ComplexSequenceEval * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_ComplexSequenceEvalSet(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_ComplexSequenceEvalSet (list,++last,i); 
      }  
    swap_ComplexSequenceEvalSet (list,left,last);    
    qsort_ComplexSequenceEvalSet(list,left,last-1,comp); 
    qsort_ComplexSequenceEvalSet(list,last+1,right,comp);    
}    


/* Function:  sort_ComplexSequenceEvalSet(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_ComplexSequenceEvalSet
 *
 *
 * Arg:         obj [UNKN ] Object containing list [ComplexSequenceEvalSet *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_ComplexSequenceEvalSet(ComplexSequenceEvalSet * obj,int (*comp)(ComplexSequenceEval *, ComplexSequenceEval *)) 
{
    qsort_ComplexSequenceEvalSet(obj->cse,0,obj->len-1,comp);    
    return;  
}    


/* Function:  expand_ComplexSequenceEvalSet(obj,len)
 *
 * Descrip:    Really an internal function for add_ComplexSequenceEvalSet
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [ComplexSequenceEvalSet *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_ComplexSequenceEvalSet(ComplexSequenceEvalSet * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_ComplexSequenceEvalSet called with no need"); 
      return TRUE;   
      }  


    if( (obj->cse = (ComplexSequenceEval ** ) ckrealloc (obj->cse,sizeof(ComplexSequenceEval *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_ComplexSequenceEvalSet, returning FALSE");   
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_ComplexSequenceEvalSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [ComplexSequenceEvalSet *]
 * Arg:        add [OWNER] Object to add to the list [ComplexSequenceEval *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_ComplexSequenceEvalSet(ComplexSequenceEvalSet * obj,ComplexSequenceEval * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_ComplexSequenceEvalSet(obj,obj->len + ComplexSequenceEvalSetLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->cse[obj->len++]=add;    
    return TRUE; 
}    


/* Function:  flush_ComplexSequenceEvalSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [ComplexSequenceEvalSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_ComplexSequenceEvalSet(ComplexSequenceEvalSet * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->cse[i] != NULL)   {  
        free_ComplexSequenceEval(obj->cse[i]);   
        obj->cse[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  ComplexSequenceEvalSet_alloc_std(void)
 *
 * Descrip:    Equivalent to ComplexSequenceEvalSet_alloc_len(ComplexSequenceEvalSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequenceEvalSet *]
 *
 */
ComplexSequenceEvalSet * ComplexSequenceEvalSet_alloc_std(void) 
{
    return ComplexSequenceEvalSet_alloc_len(ComplexSequenceEvalSetLISTLENGTH);   
}    


/* Function:  ComplexSequenceEvalSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequenceEvalSet *]
 *
 */
ComplexSequenceEvalSet * ComplexSequenceEvalSet_alloc_len(int len) 
{
    ComplexSequenceEvalSet * out;   /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = ComplexSequenceEvalSet_alloc()) == NULL)   
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->cse = (ComplexSequenceEval ** ) ckcalloc (len,sizeof(ComplexSequenceEval *))) == NULL)  {  
      warn("Warning, ckcalloc failed in ComplexSequenceEvalSet_alloc_len");  
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_ComplexSequenceEvalSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ComplexSequenceEvalSet *]
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequenceEvalSet *]
 *
 */
ComplexSequenceEvalSet * hard_link_ComplexSequenceEvalSet(ComplexSequenceEvalSet * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a ComplexSequenceEvalSet object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  ComplexSequenceEvalSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequenceEvalSet *]
 *
 */
ComplexSequenceEvalSet * ComplexSequenceEvalSet_alloc(void) 
{
    ComplexSequenceEvalSet * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(ComplexSequenceEvalSet *) ckalloc (sizeof(ComplexSequenceEvalSet))) == NULL)    {  
      warn("ComplexSequenceEvalSet_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->type = 0;   
    out->has_been_prepared = FALSE;  
    out->left_window = 0;    
    out->right_window = 0;   
    out->left_lookback = 0;  
    out->cse = NULL; 
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_ComplexSequenceEvalSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ComplexSequenceEvalSet *]
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequenceEvalSet *]
 *
 */
ComplexSequenceEvalSet * free_ComplexSequenceEvalSet(ComplexSequenceEvalSet * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a ComplexSequenceEvalSet obj. Should be trappable");    
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
    if( obj->cse != NULL)    {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->cse[i] != NULL) 
          free_ComplexSequenceEval(obj->cse[i]); 
        }  
      ckfree(obj->cse);  
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_ComplexSequence(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ComplexSequence *]
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequence *]
 *
 */
ComplexSequence * hard_link_ComplexSequence(ComplexSequence * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a ComplexSequence object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  ComplexSequence_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequence *]
 *
 */
ComplexSequence * ComplexSequence_alloc(void) 
{
    ComplexSequence * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(ComplexSequence *) ckalloc (sizeof(ComplexSequence))) == NULL)  {  
      warn("ComplexSequence_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->type = 0;   
    out->seq = NULL; 
    out->datastore = NULL;   
    out->depth = 0;  
    out->length = 0; 
    out->creator = NULL; 


    return out;  
}    


/* Function:  free_ComplexSequence(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ComplexSequence *]
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequence *]
 *
 */
ComplexSequence * free_ComplexSequence(ComplexSequence * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a ComplexSequence obj. Should be trappable");   
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
    /* obj->data is linked in */ 
    if( obj->datastore != NULL)  
      ckfree(obj->datastore);    
    if( obj->creator != NULL)    
      free_ComplexSequenceEvalSet(obj->creator);     


    ckfree(obj); 
    return NULL; 
}    


/* Function:  replace_type_ComplexSequence(obj,type)
 *
 * Descrip:    Replace member variable type
 *             For use principly by API functions
 *
 *
 * Arg:         obj [UNKN ] Object holding the variable [ComplexSequence *]
 * Arg:        type [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable type [boolean]
 *
 */
boolean replace_type_ComplexSequence(ComplexSequence * obj,int type) 
{
    if( obj == NULL)     {  
      warn("In replacement function type for object ComplexSequence, got a NULL object");    
      return FALSE;  
      }  
    obj->type = type;    
    return TRUE; 
}    


/* Function:  access_type_ComplexSequence(obj)
 *
 * Descrip:    Access member variable type
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [ComplexSequence *]
 *
 * Return [SOFT ]  member variable type [int]
 *
 */
int access_type_ComplexSequence(ComplexSequence * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function type for object ComplexSequence, got a NULL object");   
      return 0;  
      }  
    return obj->type;    
}    


/* Function:  replace_seq_ComplexSequence(obj,seq)
 *
 * Descrip:    Replace member variable seq
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [ComplexSequence *]
 * Arg:        seq [OWNER] New value of the variable [Sequence *]
 *
 * Return [SOFT ]  member variable seq [boolean]
 *
 */
boolean replace_seq_ComplexSequence(ComplexSequence * obj,Sequence * seq) 
{
    if( obj == NULL)     {  
      warn("In replacement function seq for object ComplexSequence, got a NULL object"); 
      return FALSE;  
      }  
    obj->seq = seq;  
    return TRUE; 
}    


/* Function:  access_seq_ComplexSequence(obj)
 *
 * Descrip:    Access member variable seq
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [ComplexSequence *]
 *
 * Return [SOFT ]  member variable seq [Sequence *]
 *
 */
Sequence * access_seq_ComplexSequence(ComplexSequence * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function seq for object ComplexSequence, got a NULL object");    
      return NULL;   
      }  
    return obj->seq;     
}    


/* Function:  replace_type_ComplexSequenceEvalSet(obj,type)
 *
 * Descrip:    Replace member variable type
 *             For use principly by API functions
 *
 *
 * Arg:         obj [UNKN ] Object holding the variable [ComplexSequenceEvalSet *]
 * Arg:        type [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable type [boolean]
 *
 */
boolean replace_type_ComplexSequenceEvalSet(ComplexSequenceEvalSet * obj,int type) 
{
    if( obj == NULL)     {  
      warn("In replacement function type for object ComplexSequenceEvalSet, got a NULL object"); 
      return FALSE;  
      }  
    obj->type = type;    
    return TRUE; 
}    


/* Function:  access_type_ComplexSequenceEvalSet(obj)
 *
 * Descrip:    Access member variable type
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [ComplexSequenceEvalSet *]
 *
 * Return [SOFT ]  member variable type [int]
 *
 */
int access_type_ComplexSequenceEvalSet(ComplexSequenceEvalSet * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function type for object ComplexSequenceEvalSet, got a NULL object");    
      return 0;  
      }  
    return obj->type;    
}    


/* Function:  replace_has_been_prepared_ComplexSequenceEvalSet(obj,has_been_prepared)
 *
 * Descrip:    Replace member variable has_been_prepared
 *             For use principly by API functions
 *
 *
 * Arg:                      obj [UNKN ] Object holding the variable [ComplexSequenceEvalSet *]
 * Arg:        has_been_prepared [OWNER] New value of the variable [boolean]
 *
 * Return [SOFT ]  member variable has_been_prepared [boolean]
 *
 */
boolean replace_has_been_prepared_ComplexSequenceEvalSet(ComplexSequenceEvalSet * obj,boolean has_been_prepared) 
{
    if( obj == NULL)     {  
      warn("In replacement function has_been_prepared for object ComplexSequenceEvalSet, got a NULL object");    
      return FALSE;  
      }  
    obj->has_been_prepared = has_been_prepared;  
    return TRUE; 
}    


/* Function:  access_has_been_prepared_ComplexSequenceEvalSet(obj)
 *
 * Descrip:    Access member variable has_been_prepared
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [ComplexSequenceEvalSet *]
 *
 * Return [SOFT ]  member variable has_been_prepared [boolean]
 *
 */
boolean access_has_been_prepared_ComplexSequenceEvalSet(ComplexSequenceEvalSet * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function has_been_prepared for object ComplexSequenceEvalSet, got a NULL object");   
      return FALSE;  
      }  
    return obj->has_been_prepared;   
}    


/* Function:  replace_left_window_ComplexSequenceEvalSet(obj,left_window)
 *
 * Descrip:    Replace member variable left_window
 *             For use principly by API functions
 *
 *
 * Arg:                obj [UNKN ] Object holding the variable [ComplexSequenceEvalSet *]
 * Arg:        left_window [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable left_window [boolean]
 *
 */
boolean replace_left_window_ComplexSequenceEvalSet(ComplexSequenceEvalSet * obj,int left_window) 
{
    if( obj == NULL)     {  
      warn("In replacement function left_window for object ComplexSequenceEvalSet, got a NULL object");  
      return FALSE;  
      }  
    obj->left_window = left_window;  
    return TRUE; 
}    


/* Function:  access_left_window_ComplexSequenceEvalSet(obj)
 *
 * Descrip:    Access member variable left_window
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [ComplexSequenceEvalSet *]
 *
 * Return [SOFT ]  member variable left_window [int]
 *
 */
int access_left_window_ComplexSequenceEvalSet(ComplexSequenceEvalSet * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function left_window for object ComplexSequenceEvalSet, got a NULL object"); 
      return 0;  
      }  
    return obj->left_window;     
}    


/* Function:  replace_right_window_ComplexSequenceEvalSet(obj,right_window)
 *
 * Descrip:    Replace member variable right_window
 *             For use principly by API functions
 *
 *
 * Arg:                 obj [UNKN ] Object holding the variable [ComplexSequenceEvalSet *]
 * Arg:        right_window [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable right_window [boolean]
 *
 */
boolean replace_right_window_ComplexSequenceEvalSet(ComplexSequenceEvalSet * obj,int right_window) 
{
    if( obj == NULL)     {  
      warn("In replacement function right_window for object ComplexSequenceEvalSet, got a NULL object"); 
      return FALSE;  
      }  
    obj->right_window = right_window;    
    return TRUE; 
}    


/* Function:  access_right_window_ComplexSequenceEvalSet(obj)
 *
 * Descrip:    Access member variable right_window
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [ComplexSequenceEvalSet *]
 *
 * Return [SOFT ]  member variable right_window [int]
 *
 */
int access_right_window_ComplexSequenceEvalSet(ComplexSequenceEvalSet * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function right_window for object ComplexSequenceEvalSet, got a NULL object");    
      return 0;  
      }  
    return obj->right_window;    
}    


/* Function:  replace_left_lookback_ComplexSequenceEvalSet(obj,left_lookback)
 *
 * Descrip:    Replace member variable left_lookback
 *             For use principly by API functions
 *
 *
 * Arg:                  obj [UNKN ] Object holding the variable [ComplexSequenceEvalSet *]
 * Arg:        left_lookback [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable left_lookback [boolean]
 *
 */
boolean replace_left_lookback_ComplexSequenceEvalSet(ComplexSequenceEvalSet * obj,int left_lookback) 
{
    if( obj == NULL)     {  
      warn("In replacement function left_lookback for object ComplexSequenceEvalSet, got a NULL object");    
      return FALSE;  
      }  
    obj->left_lookback = left_lookback;  
    return TRUE; 
}    


/* Function:  access_left_lookback_ComplexSequenceEvalSet(obj)
 *
 * Descrip:    Access member variable left_lookback
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [ComplexSequenceEvalSet *]
 *
 * Return [SOFT ]  member variable left_lookback [int]
 *
 */
int access_left_lookback_ComplexSequenceEvalSet(ComplexSequenceEvalSet * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function left_lookback for object ComplexSequenceEvalSet, got a NULL object");   
      return 0;  
      }  
    return obj->left_lookback;   
}    



#ifdef _cplusplus
}
#endif
