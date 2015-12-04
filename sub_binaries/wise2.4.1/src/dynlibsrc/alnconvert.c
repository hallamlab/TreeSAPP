#ifdef _cplusplus
extern "C" {
#endif
#include "alnconvert.h"

 static char * unknown_label = "UNKNOWN_LABEL";


/* Function:  add_collapse_label_AlnConvertSet(acs,label1,label2)
 *
 * Descrip:    Not a sensible function. Makes the convert with label1 and label2
 *             a collapsable label
 *
 *
 * Arg:           acs [UNKN ] Undocumented argument [AlnConvertSet *]
 * Arg:        label1 [UNKN ] Undocumented argument [char *]
 * Arg:        label2 [UNKN ] Undocumented argument [char *]
 *
 */
# line 36 "alnconvert.dy"
void add_collapse_label_AlnConvertSet(AlnConvertSet * acs,char * label1,char * label2)
{
  int i;

  for(i=0;i<acs->len;i++) 
    if( strcmp(acs->acu[i]->label1,label1) == 0 && strcmp(acs->acu[i]->label2,label2) == 0 )
       acs->acu[i]->can_collapse = TRUE;

}

/* Function:  AlnBlock_from_PackAln(acs,*pal)
 *
 * Descrip:    Takes a AlnConvertSet (acs) and a PackAln (pal) 
 *             and blindly converts it to AlnBlock. This is really
 *             an internal for a dynamite produced dy function
 *
 *
 * Arg:         acs [UNKN ] Undocumented argument [AlnConvertSet *]
 * Arg:        *pal [UNKN ] Undocumented argument [PackAln]
 *
 * Return [UNKN ]  Undocumented return value [AlnBlock  *]
 *
 */
# line 51 "alnconvert.dy"
AlnBlock  * AlnBlock_from_PackAln(AlnConvertSet * acs,PackAln *pal)
{
  AlnBlock * alb;
  AlnColumn * prev;
  AlnColumn * new;
  boolean coll;
  int i;

  alb = AlnBlock_alloc_len(2);

  add_AlnBlock(alb,AlnSequence_alloc());
  add_AlnBlock(alb,AlnSequence_alloc());

  prev = NULL;
  alb->score = pal->score;
  for(i=1;i<pal->len;i++) {
    coll = FALSE;
    new=AlnColumn_from_Pal_Convert(acs,pal->pau[i-1],pal->pau[i],prev,&coll);
    if( new == NULL ) {
      if( coll == FALSE ) {
	warn("Unrecoverable error in converting PackAln to AlnBlock... bugging out with partial alignment!");
	return alb;
      }
      else {
	/** ok, was collapsed, just loop back **/
	continue;
      }
    }
    else {
      /** got a new AlnColumn **/
      if( prev == NULL ) {
	/** then first col **/

	alb->start = new;
	alb->seq[0]->start = new->alu[0];
	alb->seq[1]->start = new->alu[1];
	prev = new;
	
      }
      else {
	prev->next = new;
	prev->alu[0]->next = new->alu[0];
	prev->alu[1]->next = new->alu[1];
	prev = new;
      }
    }

  }

  return alb;
}
    

/* Function:  AlnColumn_from_Pal_Convert(acs,before,after,prev,was_collapsed)
 *
 * Descrip:    the core of the conversion.
 *
 *
 * Arg:                  acs [UNKN ] Undocumented argument [AlnConvertSet *]
 * Arg:               before [UNKN ] Undocumented argument [PackAlnUnit *]
 * Arg:                after [UNKN ] Undocumented argument [PackAlnUnit *]
 * Arg:                 prev [UNKN ] Undocumented argument [AlnColumn *]
 * Arg:        was_collapsed [UNKN ] Undocumented argument [boolean *]
 *
 * Return [UNKN ]  Undocumented return value [AlnColumn *]
 *
 */
# line 108 "alnconvert.dy"
AlnColumn * AlnColumn_from_Pal_Convert(AlnConvertSet * acs,PackAlnUnit * before,PackAlnUnit * after,AlnColumn * prev,boolean * was_collapsed)
{
  AlnConvertUnit * acu;
  AlnColumn * alc;


  acu = AlnConvertUnit_from_state_and_offset(acs,before->state,after->state,after->i - before->i,after->j - before->j);

  if( acu == NULL) {
    warn("Between state [%d,%d,%d] and [%d,%d,%d] got no labels... labelling as UNKNOWN",before->i,before->j,before->state,after->i,after->j,after->state);
    alc = new_pairwise_AlnColumn();
    
    alc->alu[0]->start = before->i;
    alc->alu[0]->end   = after->i;

    alc->alu[1]->start = before->j;
    alc->alu[1]->end   = after->j;
    alc->alu[0]->score[0] = alc->alu[1]->score[0] = after->score;

    alc->alu[1]->text_label = alc->alu[0]->text_label = unknown_label;

    return alc;
  }

  if( acu->can_collapse == TRUE && prev != NULL && strcmp(prev->alu[0]->text_label,acu->label1) == 0 && strcmp(prev->alu[1]->text_label,acu->label2) == 0 ) {
    /*** don't return something, just add into the next one ***/
    prev->alu[0]->end = after->i;
    prev->alu[1]->end = after->j;
    prev->alu[0]->score[0] += after->score;
    prev->alu[1]->score[0] += after->score;
    if( was_collapsed != NULL ) {
      *was_collapsed = TRUE;
    }
    return NULL;
  }

  /*** else, put away this unit ***/

  alc = new_pairwise_AlnColumn();

  if( acu->is_from_special == TRUE ) {
    alc->alu[0]->start = after->i -1;
    alc->alu[0]->end   = after->i;
  } else {
    alc->alu[0]->start = before->i;
    alc->alu[0]->end   = after->i;
  }

  alc->alu[1]->start = before->j;
  alc->alu[1]->end   = after->j;
  alc->alu[0]->score[0] = alc->alu[1]->score[0] = after->score;

  alc->alu[0]->text_label = acu->label1;
  alc->alu[1]->text_label = acu->label2;

  return alc;  
}


/* Function:  AlnConvertUnit_from_state_and_offset(acs,state1,state2,offi,offj)
 *
 * Descrip:    Finds the correct AlnConvertUnit for this state,state,offi,offj 
 *             quad
 *
 *
 * Arg:           acs [UNKN ] Undocumented argument [AlnConvertSet *]
 * Arg:        state1 [UNKN ] Undocumented argument [int]
 * Arg:        state2 [UNKN ] Undocumented argument [int]
 * Arg:          offi [UNKN ] Undocumented argument [int]
 * Arg:          offj [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [AlnConvertUnit *]
 *
 */
# line 172 "alnconvert.dy"
AlnConvertUnit * AlnConvertUnit_from_state_and_offset(AlnConvertSet * acs,int state1,int state2,int offi,int offj)
{
  register int i;

  for(i=0;i<acs->len;i++) {
    if( acs->acu[i]->state1 == state1 && acs->acu[i]->state2 == state2 && (acs->acu[i]->offi == -1 || acs->acu[i]->offi == offi) && offj == acs->acu[i]->offj)
      return acs->acu[i];
  }
  return NULL;
}


# line 191 "alnconvert.c"
/* Function:  hard_link_AlnConvertUnit(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AlnConvertUnit *]
 *
 * Return [UNKN ]  Undocumented return value [AlnConvertUnit *]
 *
 */
AlnConvertUnit * hard_link_AlnConvertUnit(AlnConvertUnit * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a AlnConvertUnit object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  AlnConvertUnit_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AlnConvertUnit *]
 *
 */
AlnConvertUnit * AlnConvertUnit_alloc(void) 
{
    AlnConvertUnit * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(AlnConvertUnit *) ckalloc (sizeof(AlnConvertUnit))) == NULL)    {  
      warn("AlnConvertUnit_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->state1 = 0; 
    out->state2 = 0; 
    out->offi = 0;   
    out->offj = 0;   
    out->can_collapse = FALSE;   
    out->is_from_special = FALSE;    


    return out;  
}    


/* Function:  free_AlnConvertUnit(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AlnConvertUnit *]
 *
 * Return [UNKN ]  Undocumented return value [AlnConvertUnit *]
 *
 */
AlnConvertUnit * free_AlnConvertUnit(AlnConvertUnit * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a AlnConvertUnit obj. Should be trappable");    
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
    /* obj->label1 is linked in */ 
    /* obj->label2 is linked in */ 


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_AlnConvertSet(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_AlnConvertSet
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [AlnConvertUnit **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_AlnConvertSet(AlnConvertUnit ** list,int i,int j)  
{
    AlnConvertUnit * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_AlnConvertSet(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_AlnConvertSet which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [AlnConvertUnit **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_AlnConvertSet(AlnConvertUnit ** list,int left,int right,int (*comp)(AlnConvertUnit * ,AlnConvertUnit * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_AlnConvertSet(list,left,(left+right)/2);    
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_AlnConvertSet (list,++last,i);  
      }  
    swap_AlnConvertSet (list,left,last); 
    qsort_AlnConvertSet(list,left,last-1,comp);  
    qsort_AlnConvertSet(list,last+1,right,comp); 
}    


/* Function:  sort_AlnConvertSet(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_AlnConvertSet
 *
 *
 * Arg:         obj [UNKN ] Object containing list [AlnConvertSet *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_AlnConvertSet(AlnConvertSet * obj,int (*comp)(AlnConvertUnit *, AlnConvertUnit *)) 
{
    qsort_AlnConvertSet(obj->acu,0,obj->len-1,comp); 
    return;  
}    


/* Function:  expand_AlnConvertSet(obj,len)
 *
 * Descrip:    Really an internal function for add_AlnConvertSet
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [AlnConvertSet *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_AlnConvertSet(AlnConvertSet * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_AlnConvertSet called with no need");  
      return TRUE;   
      }  


    if( (obj->acu = (AlnConvertUnit ** ) ckrealloc (obj->acu,sizeof(AlnConvertUnit *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_AlnConvertSet, returning FALSE");    
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_AlnConvertSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [AlnConvertSet *]
 * Arg:        add [OWNER] Object to add to the list [AlnConvertUnit *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_AlnConvertSet(AlnConvertSet * obj,AlnConvertUnit * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_AlnConvertSet(obj,obj->len + AlnConvertSetLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->acu[obj->len++]=add;    
    return TRUE; 
}    


/* Function:  flush_AlnConvertSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [AlnConvertSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_AlnConvertSet(AlnConvertSet * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->acu[i] != NULL)   {  
        free_AlnConvertUnit(obj->acu[i]);    
        obj->acu[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  AlnConvertSet_alloc_std(void)
 *
 * Descrip:    Equivalent to AlnConvertSet_alloc_len(AlnConvertSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AlnConvertSet *]
 *
 */
AlnConvertSet * AlnConvertSet_alloc_std(void) 
{
    return AlnConvertSet_alloc_len(AlnConvertSetLISTLENGTH); 
}    


/* Function:  AlnConvertSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [AlnConvertSet *]
 *
 */
AlnConvertSet * AlnConvertSet_alloc_len(int len) 
{
    AlnConvertSet * out;/* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = AlnConvertSet_alloc()) == NULL)    
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->acu = (AlnConvertUnit ** ) ckcalloc (len,sizeof(AlnConvertUnit *))) == NULL)    {  
      warn("Warning, ckcalloc failed in AlnConvertSet_alloc_len");   
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_AlnConvertSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AlnConvertSet *]
 *
 * Return [UNKN ]  Undocumented return value [AlnConvertSet *]
 *
 */
AlnConvertSet * hard_link_AlnConvertSet(AlnConvertSet * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a AlnConvertSet object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  AlnConvertSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AlnConvertSet *]
 *
 */
AlnConvertSet * AlnConvertSet_alloc(void) 
{
    AlnConvertSet * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(AlnConvertSet *) ckalloc (sizeof(AlnConvertSet))) == NULL)  {  
      warn("AlnConvertSet_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->acu = NULL; 
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_AlnConvertSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AlnConvertSet *]
 *
 * Return [UNKN ]  Undocumented return value [AlnConvertSet *]
 *
 */
AlnConvertSet * free_AlnConvertSet(AlnConvertSet * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a AlnConvertSet obj. Should be trappable"); 
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
    if( obj->acu != NULL)    {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->acu[i] != NULL) 
          free_AlnConvertUnit(obj->acu[i]);  
        }  
      ckfree(obj->acu);  
      }  


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
