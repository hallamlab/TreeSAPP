#ifdef _cplusplus
extern "C" {
#endif
#include "dnahmm.h"


# line 55 "dnahmm.dy"
DnaHmmScore * DnaHmmScore_from_DnaHmmProb(DnaHmmProb * dhp)
{
  int i;
  DnaHmmScore * out;

  out = DnaHmmScore_alloc_len(dhp->len);
  
  for(i=0;i<dhp->len;i++) {
    add_DnaHmmScore(out,DnaHmmScoreUnit_from_DnaHmmProbUnit(dhp->unit[i]));
  }

  return out;
}

# line 69 "dnahmm.dy"
void make_consensus_DnaHmmProb(DnaHmmProb * dhp)
{
  int i;
  int j;
  double prob_max;
  int base;

  dhp->consensus = calloc(dhp->len+1,sizeof(char));

  for(i=0;i<dhp->len;i++) {
    prob_max = -1.0;
    for(j=0;j<4;j++) {
      if( dhp->unit[i]->match[j] > prob_max ) {
	prob_max =  dhp->unit[i]->match[j];
	base = j;
      }
    }
    if( prob_max > 0.5 ) {
      dhp->consensus[i] = char_from_base(base);
    } else if ( prob_max > 0.3 ) {
      dhp->consensus[i] = tolower(char_from_base(base));
    } else {
      dhp->consensus[i] = '.';
    }
  }
  dhp->consensus[i] = '\0';

}

# line 98 "dnahmm.dy"
DnaHmmScoreUnit * DnaHmmScoreUnit_from_DnaHmmProbUnit(DnaHmmProbUnit * dpu)
{
  DnaHmmScoreUnit * out;

  out = DnaHmmScoreUnit_alloc();

  Probability2Score_move(dpu->match,out->match,DNA_HMM_ALPHABET);
  Probability2Score_move(dpu->insert,out->insert,DNA_HMM_ALPHABET);
  Probability2Score_move(dpu->transition,out->transition,DHMM_TRANSITION_LEN);

 
  return out;
}

# line 112 "dnahmm.dy"
DnaHmmProb * new_DnaHmmProb_from_SeqAlign_ungapped(SeqAlign * sa,double simple_pseudocount)
{
  DnaHmmProb * out;
  ColumnCount * cc;
  int i;

  assert(sa);

  out = DnaHmmProb_alloc_std();
  out->ref_align = hard_link_SeqAlign(sa);

  for(i=0;i<sa->seq[0]->len;i++) {
    cc = ColumnCount_from_SeqAlign(sa,i);
    add_DnaHmmProb(out,new_DnaHmmProbUnit_from_ColumnCount_ungapped(cc,simple_pseudocount));

    free_ColumnCount(cc);
  }
  
  out->unit[0]->transition[DHMM_START2MATCH] = 1.0;
  out->unit[out->len-1]->transition[DHMM_MATCH2END] = 1.0;


  return out;
}


# line 138 "dnahmm.dy"
DnaHmmProbUnit * new_DnaHmmProbUnit_from_ColumnCount_ungapped(ColumnCount * cc,double simple_pseudocount)
{
  DnaHmmProbUnit * out;
  double total;
  char base[] = "ATGC";
  int i;

  out = DnaHmmProbUnit_alloc();
  
  for(i=0,total = 0.0;i<4;i++)
    total += (cc->count[base[i]-'A'] + simple_pseudocount);

  
  for(i=0;i<4;i++)
    out->match[base_from_char(base[i])] = (cc->count[base[i]-'A'] + simple_pseudocount) / total;
  out->match[4] = 1.0;

  for(i=0;i<DHMM_TRANSITION_LEN;i++) {
    out->transition[i] = 0.0;
  }
  out->transition[DHMM_MATCH2MATCH] = 1.0;

  return out;
}


# line 164 "dnahmm.dy"
void set_N_DnaHmmProb(DnaHmmProb * dhp,Probability basen)
{
  int i;

  for(i=0;i<dhp->len;i++) {
    dhp->unit[i]->match[BASE_N] = basen;
    dhp->unit[i]->insert[BASE_N] = basen;
  }

}

# line 175 "dnahmm.dy"
void fold_RandomModelDNA_DnaHmmProb(DnaHmmProb * dhp,RandomModelDNA * d,Probability rnd_advance)
{
  int i;
  int j;

  for(i=0;i<dhp->len;i++) {
    for(j=0;j<DNA_HMM_ALPHABET;j++) {
      dhp->unit[i]->match[j] = dhp->unit[i]->match[j] / d->base[j];
      dhp->unit[i]->insert[j] = dhp->unit[i]->insert[j] / d->base[j];
    }

    for(j=0;j<DHMM_TRANSITION_LEN;j++) {
      dhp->unit[i]->transition[j] = dhp->unit[i]->transition[j] / rnd_advance;
    }
  }
}



# line 152 "dnahmm.c"
/* Function:  hard_link_DnaHmmProbUnit(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DnaHmmProbUnit *]
 *
 * Return [UNKN ]  Undocumented return value [DnaHmmProbUnit *]
 *
 */
DnaHmmProbUnit * hard_link_DnaHmmProbUnit(DnaHmmProbUnit * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a DnaHmmProbUnit object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  DnaHmmProbUnit_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaHmmProbUnit *]
 *
 */
DnaHmmProbUnit * DnaHmmProbUnit_alloc(void) 
{
    DnaHmmProbUnit * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(DnaHmmProbUnit *) ckalloc (sizeof(DnaHmmProbUnit))) == NULL)    {  
      warn("DnaHmmProbUnit_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    /* match[DNA_HMM_ALPHABET] is an array: no default possible */ 
    /* insert[DNA_HMM_ALPHABET] is an array: no default possible */ 
    /* transition[DHMM_TRANSITION_LEN] is an array: no default possible */ 


    return out;  
}    


/* Function:  free_DnaHmmProbUnit(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DnaHmmProbUnit *]
 *
 * Return [UNKN ]  Undocumented return value [DnaHmmProbUnit *]
 *
 */
DnaHmmProbUnit * free_DnaHmmProbUnit(DnaHmmProbUnit * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a DnaHmmProbUnit obj. Should be trappable");    
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


/* Function:  swap_DnaHmmProb(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_DnaHmmProb
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [DnaHmmProbUnit **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_DnaHmmProb(DnaHmmProbUnit ** list,int i,int j)  
{
    DnaHmmProbUnit * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_DnaHmmProb(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_DnaHmmProb which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [DnaHmmProbUnit **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_DnaHmmProb(DnaHmmProbUnit ** list,int left,int right,int (*comp)(DnaHmmProbUnit * ,DnaHmmProbUnit * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_DnaHmmProb(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_DnaHmmProb (list,++last,i); 
      }  
    swap_DnaHmmProb (list,left,last);    
    qsort_DnaHmmProb(list,left,last-1,comp); 
    qsort_DnaHmmProb(list,last+1,right,comp);    
}    


/* Function:  sort_DnaHmmProb(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_DnaHmmProb
 *
 *
 * Arg:         obj [UNKN ] Object containing list [DnaHmmProb *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_DnaHmmProb(DnaHmmProb * obj,int (*comp)(DnaHmmProbUnit *, DnaHmmProbUnit *)) 
{
    qsort_DnaHmmProb(obj->unit,0,obj->len-1,comp);   
    return;  
}    


/* Function:  expand_DnaHmmProb(obj,len)
 *
 * Descrip:    Really an internal function for add_DnaHmmProb
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [DnaHmmProb *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_DnaHmmProb(DnaHmmProb * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_DnaHmmProb called with no need"); 
      return TRUE;   
      }  


    if( (obj->unit = (DnaHmmProbUnit ** ) ckrealloc (obj->unit,sizeof(DnaHmmProbUnit *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_DnaHmmProb, returning FALSE");   
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_DnaHmmProb(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [DnaHmmProb *]
 * Arg:        add [OWNER] Object to add to the list [DnaHmmProbUnit *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_DnaHmmProb(DnaHmmProb * obj,DnaHmmProbUnit * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_DnaHmmProb(obj,obj->len + DnaHmmProbLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->unit[obj->len++]=add;   
    return TRUE; 
}    


/* Function:  flush_DnaHmmProb(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [DnaHmmProb *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_DnaHmmProb(DnaHmmProb * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->unit[i] != NULL)  {  
        free_DnaHmmProbUnit(obj->unit[i]);   
        obj->unit[i] = NULL; 
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  DnaHmmProb_alloc_std(void)
 *
 * Descrip:    Equivalent to DnaHmmProb_alloc_len(DnaHmmProbLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaHmmProb *]
 *
 */
DnaHmmProb * DnaHmmProb_alloc_std(void) 
{
    return DnaHmmProb_alloc_len(DnaHmmProbLISTLENGTH);   
}    


/* Function:  DnaHmmProb_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [DnaHmmProb *]
 *
 */
DnaHmmProb * DnaHmmProb_alloc_len(int len) 
{
    DnaHmmProb * out;   /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = DnaHmmProb_alloc()) == NULL)   
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->unit = (DnaHmmProbUnit ** ) ckcalloc (len,sizeof(DnaHmmProbUnit *))) == NULL)   {  
      warn("Warning, ckcalloc failed in DnaHmmProb_alloc_len");  
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_DnaHmmProb(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DnaHmmProb *]
 *
 * Return [UNKN ]  Undocumented return value [DnaHmmProb *]
 *
 */
DnaHmmProb * hard_link_DnaHmmProb(DnaHmmProb * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a DnaHmmProb object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  DnaHmmProb_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaHmmProb *]
 *
 */
DnaHmmProb * DnaHmmProb_alloc(void) 
{
    DnaHmmProb * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(DnaHmmProb *) ckalloc (sizeof(DnaHmmProb))) == NULL)    {  
      warn("DnaHmmProb_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->unit = NULL;    
    out->len = out->maxlen = 0;  
    out->ref_align = NULL;   
    out->consensus = NULL;   


    return out;  
}    


/* Function:  free_DnaHmmProb(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DnaHmmProb *]
 *
 * Return [UNKN ]  Undocumented return value [DnaHmmProb *]
 *
 */
DnaHmmProb * free_DnaHmmProb(DnaHmmProb * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a DnaHmmProb obj. Should be trappable");    
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
    if( obj->unit != NULL)   {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->unit[i] != NULL)    
          free_DnaHmmProbUnit(obj->unit[i]); 
        }  
      ckfree(obj->unit); 
      }  
    if( obj->ref_align != NULL)  
      free_SeqAlign(obj->ref_align);     
    if( obj->consensus != NULL)  
      ckfree(obj->consensus);    


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_DnaHmmScoreUnit(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DnaHmmScoreUnit *]
 *
 * Return [UNKN ]  Undocumented return value [DnaHmmScoreUnit *]
 *
 */
DnaHmmScoreUnit * hard_link_DnaHmmScoreUnit(DnaHmmScoreUnit * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a DnaHmmScoreUnit object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  DnaHmmScoreUnit_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaHmmScoreUnit *]
 *
 */
DnaHmmScoreUnit * DnaHmmScoreUnit_alloc(void) 
{
    DnaHmmScoreUnit * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(DnaHmmScoreUnit *) ckalloc (sizeof(DnaHmmScoreUnit))) == NULL)  {  
      warn("DnaHmmScoreUnit_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    /* match[DNA_HMM_ALPHABET] is an array: no default possible */ 
    /* insert[DNA_HMM_ALPHABET] is an array: no default possible */ 
    /* transition[DHMM_TRANSITION_LEN] is an array: no default possible */ 


    return out;  
}    


/* Function:  free_DnaHmmScoreUnit(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DnaHmmScoreUnit *]
 *
 * Return [UNKN ]  Undocumented return value [DnaHmmScoreUnit *]
 *
 */
DnaHmmScoreUnit * free_DnaHmmScoreUnit(DnaHmmScoreUnit * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a DnaHmmScoreUnit obj. Should be trappable");   
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


/* Function:  swap_DnaHmmScore(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_DnaHmmScore
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [DnaHmmScoreUnit **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_DnaHmmScore(DnaHmmScoreUnit ** list,int i,int j)  
{
    DnaHmmScoreUnit * temp;  
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_DnaHmmScore(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_DnaHmmScore which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [DnaHmmScoreUnit **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_DnaHmmScore(DnaHmmScoreUnit ** list,int left,int right,int (*comp)(DnaHmmScoreUnit * ,DnaHmmScoreUnit * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_DnaHmmScore(list,left,(left+right)/2);  
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_DnaHmmScore (list,++last,i);    
      }  
    swap_DnaHmmScore (list,left,last);   
    qsort_DnaHmmScore(list,left,last-1,comp);    
    qsort_DnaHmmScore(list,last+1,right,comp);   
}    


/* Function:  sort_DnaHmmScore(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_DnaHmmScore
 *
 *
 * Arg:         obj [UNKN ] Object containing list [DnaHmmScore *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_DnaHmmScore(DnaHmmScore * obj,int (*comp)(DnaHmmScoreUnit *, DnaHmmScoreUnit *)) 
{
    qsort_DnaHmmScore(obj->unit,0,obj->len-1,comp);  
    return;  
}    


/* Function:  expand_DnaHmmScore(obj,len)
 *
 * Descrip:    Really an internal function for add_DnaHmmScore
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [DnaHmmScore *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_DnaHmmScore(DnaHmmScore * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_DnaHmmScore called with no need");    
      return TRUE;   
      }  


    if( (obj->unit = (DnaHmmScoreUnit ** ) ckrealloc (obj->unit,sizeof(DnaHmmScoreUnit *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_DnaHmmScore, returning FALSE");  
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_DnaHmmScore(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [DnaHmmScore *]
 * Arg:        add [OWNER] Object to add to the list [DnaHmmScoreUnit *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_DnaHmmScore(DnaHmmScore * obj,DnaHmmScoreUnit * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_DnaHmmScore(obj,obj->len + DnaHmmScoreLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->unit[obj->len++]=add;   
    return TRUE; 
}    


/* Function:  flush_DnaHmmScore(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [DnaHmmScore *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_DnaHmmScore(DnaHmmScore * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->unit[i] != NULL)  {  
        free_DnaHmmScoreUnit(obj->unit[i]);  
        obj->unit[i] = NULL; 
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  DnaHmmScore_alloc_std(void)
 *
 * Descrip:    Equivalent to DnaHmmScore_alloc_len(DnaHmmScoreLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaHmmScore *]
 *
 */
DnaHmmScore * DnaHmmScore_alloc_std(void) 
{
    return DnaHmmScore_alloc_len(DnaHmmScoreLISTLENGTH); 
}    


/* Function:  DnaHmmScore_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [DnaHmmScore *]
 *
 */
DnaHmmScore * DnaHmmScore_alloc_len(int len) 
{
    DnaHmmScore * out;  /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = DnaHmmScore_alloc()) == NULL)  
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->unit = (DnaHmmScoreUnit ** ) ckcalloc (len,sizeof(DnaHmmScoreUnit *))) == NULL) {  
      warn("Warning, ckcalloc failed in DnaHmmScore_alloc_len"); 
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_DnaHmmScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DnaHmmScore *]
 *
 * Return [UNKN ]  Undocumented return value [DnaHmmScore *]
 *
 */
DnaHmmScore * hard_link_DnaHmmScore(DnaHmmScore * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a DnaHmmScore object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  DnaHmmScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaHmmScore *]
 *
 */
DnaHmmScore * DnaHmmScore_alloc(void) 
{
    DnaHmmScore * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(DnaHmmScore *) ckalloc (sizeof(DnaHmmScore))) == NULL)  {  
      warn("DnaHmmScore_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->unit = NULL;    
    out->len = out->maxlen = 0;  
    out->ref_align = NULL;   


    return out;  
}    


/* Function:  free_DnaHmmScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DnaHmmScore *]
 *
 * Return [UNKN ]  Undocumented return value [DnaHmmScore *]
 *
 */
DnaHmmScore * free_DnaHmmScore(DnaHmmScore * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a DnaHmmScore obj. Should be trappable");   
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
    if( obj->unit != NULL)   {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->unit[i] != NULL)    
          free_DnaHmmScoreUnit(obj->unit[i]);    
        }  
      ckfree(obj->unit); 
      }  
    if( obj->ref_align != NULL)  
      free_SeqAlign(obj->ref_align);     


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
