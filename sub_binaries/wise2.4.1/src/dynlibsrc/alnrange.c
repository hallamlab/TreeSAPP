#ifdef _cplusplus
extern "C" {
#endif

#include "alnrange.h"


/* Function:  show_AlnRangeSet(ars,ofp)
 *
 * Descrip:    shows AlnRangeSet in vaguely human form
 *
 *
 * Arg:        ars [UNKN ] Undocumented argument [AlnRangeSet *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 43 "alnrange.dy"
void show_AlnRangeSet(AlnRangeSet * ars,FILE * ofp)
{
  int i;

  for(i=0;i<ars->len;i++) {
    show_AlnRange(ars->alr[i],ofp);
  }

}

/* Function:  show_AlnRange(alr,ofp)
 *
 * Descrip:    shows AlnRange in vaguely human form
 *
 *
 * Arg:        alr [UNKN ] Undocumented argument [AlnRange *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 56 "alnrange.dy"
void show_AlnRange(AlnRange * alr,FILE * ofp)
{
  fprintf(ofp,"(%d,%d,%d)-(%d,%d,%d) [%d,%d]\n",alr->starti,alr->startj,alr->startstate,
alr->stopi,alr->stopj,alr->stopstate,alr->startscore,alr->stopscore);
}

/* Function:  sort_AlnRangeSet_by_start(ars)
 *
 * Descrip:    Sorts an AlnRangeSet by start of each AlnRange
 *
 *
 * Arg:        ars [UNKN ] Undocumented argument [AlnRangeSet *]
 *
 */
# line 65 "alnrange.dy"
void sort_AlnRangeSet_by_start(AlnRangeSet * ars)
{
  sort_AlnRangeSet(ars,compare_AlnRange_start);
}

/* Function:  compare_AlnRange_start(one,two)
 *
 * Descrip:    compares to AlnRange for /sort_AlnRangeSet_by_start
 *
 *
 * Arg:        one [UNKN ] Undocumented argument [AlnRange *]
 * Arg:        two [UNKN ] Undocumented argument [AlnRange *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 74 "alnrange.dy"
int compare_AlnRange_start(AlnRange * one,AlnRange * two)
{
  if( one->startj > two->startj) 
    return 1;
  else return -1;
}


# line 69 "alnrange.c"
/* Function:  hard_link_AlnRange(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AlnRange *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRange *]
 *
 */
AlnRange * hard_link_AlnRange(AlnRange * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a AlnRange object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  AlnRange_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AlnRange *]
 *
 */
AlnRange * AlnRange_alloc(void) 
{
    AlnRange * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(AlnRange *) ckalloc (sizeof(AlnRange))) == NULL)    {  
      warn("AlnRange_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->starti = 0; 
    out->startj = 0; 
    out->startstate = 0; 
    out->stopi = 0;  
    out->stopj = 0;  
    out->stopstate = 0;  
    out->startscore = 0; 
    out->stopscore = 0;  


    return out;  
}    


/* Function:  free_AlnRange(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AlnRange *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRange *]
 *
 */
AlnRange * free_AlnRange(AlnRange * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a AlnRange obj. Should be trappable");  
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


/* Function:  swap_AlnRangeSet(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_AlnRangeSet
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [AlnRange **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_AlnRangeSet(AlnRange ** list,int i,int j)  
{
    AlnRange * temp; 
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_AlnRangeSet(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_AlnRangeSet which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [AlnRange **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_AlnRangeSet(AlnRange ** list,int left,int right,int (*comp)(AlnRange * ,AlnRange * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_AlnRangeSet(list,left,(left+right)/2);  
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_AlnRangeSet (list,++last,i);    
      }  
    swap_AlnRangeSet (list,left,last);   
    qsort_AlnRangeSet(list,left,last-1,comp);    
    qsort_AlnRangeSet(list,last+1,right,comp);   
}    


/* Function:  sort_AlnRangeSet(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_AlnRangeSet
 *
 *
 * Arg:         obj [UNKN ] Object containing list [AlnRangeSet *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_AlnRangeSet(AlnRangeSet * obj,int (*comp)(AlnRange *, AlnRange *)) 
{
    qsort_AlnRangeSet(obj->alr,0,obj->len-1,comp);   
    return;  
}    


/* Function:  expand_AlnRangeSet(obj,len)
 *
 * Descrip:    Really an internal function for add_AlnRangeSet
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [AlnRangeSet *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_AlnRangeSet(AlnRangeSet * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_AlnRangeSet called with no need");    
      return TRUE;   
      }  


    if( (obj->alr = (AlnRange ** ) ckrealloc (obj->alr,sizeof(AlnRange *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_AlnRangeSet, returning FALSE");  
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_AlnRangeSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [AlnRangeSet *]
 * Arg:        add [OWNER] Object to add to the list [AlnRange *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_AlnRangeSet(AlnRangeSet * obj,AlnRange * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_AlnRangeSet(obj,obj->len + AlnRangeSetLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->alr[obj->len++]=add;    
    return TRUE; 
}    


/* Function:  flush_AlnRangeSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [AlnRangeSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_AlnRangeSet(AlnRangeSet * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->alr[i] != NULL)   {  
        free_AlnRange(obj->alr[i]);  
        obj->alr[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  AlnRangeSet_alloc_std(void)
 *
 * Descrip:    Equivalent to AlnRangeSet_alloc_len(AlnRangeSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * AlnRangeSet_alloc_std(void) 
{
    return AlnRangeSet_alloc_len(AlnRangeSetLISTLENGTH); 
}    


/* Function:  AlnRangeSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * AlnRangeSet_alloc_len(int len) 
{
    AlnRangeSet * out;  /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = AlnRangeSet_alloc()) == NULL)  
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->alr = (AlnRange ** ) ckcalloc (len,sizeof(AlnRange *))) == NULL)    {  
      warn("Warning, ckcalloc failed in AlnRangeSet_alloc_len"); 
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_AlnRangeSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AlnRangeSet *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * hard_link_AlnRangeSet(AlnRangeSet * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a AlnRangeSet object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  AlnRangeSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * AlnRangeSet_alloc(void) 
{
    AlnRangeSet * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(AlnRangeSet *) ckalloc (sizeof(AlnRangeSet))) == NULL)  {  
      warn("AlnRangeSet_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->score = 0;  
    out->alr = NULL; 
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_AlnRangeSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AlnRangeSet *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * free_AlnRangeSet(AlnRangeSet * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a AlnRangeSet obj. Should be trappable");   
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
    if( obj->alr != NULL)    {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->alr[i] != NULL) 
          free_AlnRange(obj->alr[i]);    
        }  
      ckfree(obj->alr);  
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  replace_starti_AlnRange(obj,starti)
 *
 * Descrip:    Replace member variable starti
 *             For use principly by API functions
 *
 *
 * Arg:           obj [UNKN ] Object holding the variable [AlnRange *]
 * Arg:        starti [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable starti [boolean]
 *
 */
boolean replace_starti_AlnRange(AlnRange * obj,int starti) 
{
    if( obj == NULL)     {  
      warn("In replacement function starti for object AlnRange, got a NULL object"); 
      return FALSE;  
      }  
    obj->starti = starti;    
    return TRUE; 
}    


/* Function:  access_starti_AlnRange(obj)
 *
 * Descrip:    Access member variable starti
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [AlnRange *]
 *
 * Return [SOFT ]  member variable starti [int]
 *
 */
int access_starti_AlnRange(AlnRange * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function starti for object AlnRange, got a NULL object");    
      return 0;  
      }  
    return obj->starti;  
}    


/* Function:  replace_startj_AlnRange(obj,startj)
 *
 * Descrip:    Replace member variable startj
 *             For use principly by API functions
 *
 *
 * Arg:           obj [UNKN ] Object holding the variable [AlnRange *]
 * Arg:        startj [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable startj [boolean]
 *
 */
boolean replace_startj_AlnRange(AlnRange * obj,int startj) 
{
    if( obj == NULL)     {  
      warn("In replacement function startj for object AlnRange, got a NULL object"); 
      return FALSE;  
      }  
    obj->startj = startj;    
    return TRUE; 
}    


/* Function:  access_startj_AlnRange(obj)
 *
 * Descrip:    Access member variable startj
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [AlnRange *]
 *
 * Return [SOFT ]  member variable startj [int]
 *
 */
int access_startj_AlnRange(AlnRange * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function startj for object AlnRange, got a NULL object");    
      return 0;  
      }  
    return obj->startj;  
}    


/* Function:  replace_startstate_AlnRange(obj,startstate)
 *
 * Descrip:    Replace member variable startstate
 *             For use principly by API functions
 *
 *
 * Arg:               obj [UNKN ] Object holding the variable [AlnRange *]
 * Arg:        startstate [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable startstate [boolean]
 *
 */
boolean replace_startstate_AlnRange(AlnRange * obj,int startstate) 
{
    if( obj == NULL)     {  
      warn("In replacement function startstate for object AlnRange, got a NULL object"); 
      return FALSE;  
      }  
    obj->startstate = startstate;    
    return TRUE; 
}    


/* Function:  access_startstate_AlnRange(obj)
 *
 * Descrip:    Access member variable startstate
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [AlnRange *]
 *
 * Return [SOFT ]  member variable startstate [int]
 *
 */
int access_startstate_AlnRange(AlnRange * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function startstate for object AlnRange, got a NULL object");    
      return 0;  
      }  
    return obj->startstate;  
}    


/* Function:  replace_stopi_AlnRange(obj,stopi)
 *
 * Descrip:    Replace member variable stopi
 *             For use principly by API functions
 *
 *
 * Arg:          obj [UNKN ] Object holding the variable [AlnRange *]
 * Arg:        stopi [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable stopi [boolean]
 *
 */
boolean replace_stopi_AlnRange(AlnRange * obj,int stopi) 
{
    if( obj == NULL)     {  
      warn("In replacement function stopi for object AlnRange, got a NULL object");  
      return FALSE;  
      }  
    obj->stopi = stopi;  
    return TRUE; 
}    


/* Function:  access_stopi_AlnRange(obj)
 *
 * Descrip:    Access member variable stopi
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [AlnRange *]
 *
 * Return [SOFT ]  member variable stopi [int]
 *
 */
int access_stopi_AlnRange(AlnRange * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function stopi for object AlnRange, got a NULL object"); 
      return 0;  
      }  
    return obj->stopi;   
}    


/* Function:  replace_stopj_AlnRange(obj,stopj)
 *
 * Descrip:    Replace member variable stopj
 *             For use principly by API functions
 *
 *
 * Arg:          obj [UNKN ] Object holding the variable [AlnRange *]
 * Arg:        stopj [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable stopj [boolean]
 *
 */
boolean replace_stopj_AlnRange(AlnRange * obj,int stopj) 
{
    if( obj == NULL)     {  
      warn("In replacement function stopj for object AlnRange, got a NULL object");  
      return FALSE;  
      }  
    obj->stopj = stopj;  
    return TRUE; 
}    


/* Function:  access_stopj_AlnRange(obj)
 *
 * Descrip:    Access member variable stopj
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [AlnRange *]
 *
 * Return [SOFT ]  member variable stopj [int]
 *
 */
int access_stopj_AlnRange(AlnRange * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function stopj for object AlnRange, got a NULL object"); 
      return 0;  
      }  
    return obj->stopj;   
}    


/* Function:  replace_stopstate_AlnRange(obj,stopstate)
 *
 * Descrip:    Replace member variable stopstate
 *             For use principly by API functions
 *
 *
 * Arg:              obj [UNKN ] Object holding the variable [AlnRange *]
 * Arg:        stopstate [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable stopstate [boolean]
 *
 */
boolean replace_stopstate_AlnRange(AlnRange * obj,int stopstate) 
{
    if( obj == NULL)     {  
      warn("In replacement function stopstate for object AlnRange, got a NULL object");  
      return FALSE;  
      }  
    obj->stopstate = stopstate;  
    return TRUE; 
}    


/* Function:  access_stopstate_AlnRange(obj)
 *
 * Descrip:    Access member variable stopstate
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [AlnRange *]
 *
 * Return [SOFT ]  member variable stopstate [int]
 *
 */
int access_stopstate_AlnRange(AlnRange * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function stopstate for object AlnRange, got a NULL object"); 
      return 0;  
      }  
    return obj->stopstate;   
}    


/* Function:  replace_startscore_AlnRange(obj,startscore)
 *
 * Descrip:    Replace member variable startscore
 *             For use principly by API functions
 *
 *
 * Arg:               obj [UNKN ] Object holding the variable [AlnRange *]
 * Arg:        startscore [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable startscore [boolean]
 *
 */
boolean replace_startscore_AlnRange(AlnRange * obj,int startscore) 
{
    if( obj == NULL)     {  
      warn("In replacement function startscore for object AlnRange, got a NULL object"); 
      return FALSE;  
      }  
    obj->startscore = startscore;    
    return TRUE; 
}    


/* Function:  access_startscore_AlnRange(obj)
 *
 * Descrip:    Access member variable startscore
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [AlnRange *]
 *
 * Return [SOFT ]  member variable startscore [int]
 *
 */
int access_startscore_AlnRange(AlnRange * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function startscore for object AlnRange, got a NULL object");    
      return 0;  
      }  
    return obj->startscore;  
}    


/* Function:  replace_stopscore_AlnRange(obj,stopscore)
 *
 * Descrip:    Replace member variable stopscore
 *             For use principly by API functions
 *
 *
 * Arg:              obj [UNKN ] Object holding the variable [AlnRange *]
 * Arg:        stopscore [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable stopscore [boolean]
 *
 */
boolean replace_stopscore_AlnRange(AlnRange * obj,int stopscore) 
{
    if( obj == NULL)     {  
      warn("In replacement function stopscore for object AlnRange, got a NULL object");  
      return FALSE;  
      }  
    obj->stopscore = stopscore;  
    return TRUE; 
}    


/* Function:  access_stopscore_AlnRange(obj)
 *
 * Descrip:    Access member variable stopscore
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [AlnRange *]
 *
 * Return [SOFT ]  member variable stopscore [int]
 *
 */
int access_stopscore_AlnRange(AlnRange * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function stopscore for object AlnRange, got a NULL object"); 
      return 0;  
      }  
    return obj->stopscore;   
}    


/* Function:  replace_score_AlnRangeSet(obj,score)
 *
 * Descrip:    Replace member variable score
 *             For use principly by API functions
 *
 *
 * Arg:          obj [UNKN ] Object holding the variable [AlnRangeSet *]
 * Arg:        score [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable score [boolean]
 *
 */
boolean replace_score_AlnRangeSet(AlnRangeSet * obj,int score) 
{
    if( obj == NULL)     {  
      warn("In replacement function score for object AlnRangeSet, got a NULL object");   
      return FALSE;  
      }  
    obj->score = score;  
    return TRUE; 
}    


/* Function:  access_score_AlnRangeSet(obj)
 *
 * Descrip:    Access member variable score
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [AlnRangeSet *]
 *
 * Return [SOFT ]  member variable score [int]
 *
 */
int access_score_AlnRangeSet(AlnRangeSet * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function score for object AlnRangeSet, got a NULL object");  
      return 0;  
      }  
    return obj->score;   
}    


/* Function:  access_alr_AlnRangeSet(obj,i)
 *
 * Descrip:    Access members stored in the alr list
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the list [AlnRangeSet *]
 * Arg:          i [UNKN ] Position in the list [int]
 *
 * Return [SOFT ]  Element of the list [AlnRange *]
 *
 */
AlnRange * access_alr_AlnRangeSet(AlnRangeSet * obj,int i) 
{
    if( obj == NULL)     {  
      warn("In accessor function alr for object AlnRangeSet, got a NULL object");    
      return NULL;   
      }  
    if( obj->len <= i )  {  
      warn("In accessor function alr for object AlnRangeSet, index %%d is greater than list length %%d",i,obj->len); 
      return NULL;   
      }  
    return obj->alr[i];  
}    


/* Function:  length_alr_AlnRangeSet(obj)
 *
 * Descrip:    discover the length of the list
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the list [AlnRangeSet *]
 *
 * Return [UNKN ]  length of the list [int]
 *
 */
int length_alr_AlnRangeSet(AlnRangeSet * obj) 
{
    if( obj == NULL)     {  
      warn("In length function alr for object AlnRangeSet, got a NULL object");  
      return -1;     
      }  
    return obj->len;     
}    



#ifdef _cplusplus
}
#endif
