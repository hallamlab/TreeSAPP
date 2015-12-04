#ifdef _cplusplus
extern "C" {
#endif

#include "syexonmodel.h"

# line 31 "syexonmodel.dy"
SyExonScore * SyExonScore_flat_model(int start,int end,Probability exit,Probability final_stay)
{
  SyExonScore * out;
  SyExonModel * m;

  m = SyExonModel_flat_model(start,end,exit,final_stay);

  out = SyExonScore_from_SyExonModel(m);

  free_SyExonModel(m);

  return out;
}


# line 46 "syexonmodel.dy"
SyExonModel * SyExonModel_flat_model(int start,int end,Probability exit,Probability final_stay)
{
  int i;
  SyExon * e;
  SyExonModel * out;
  Probability prev;

  out = SyExonModel_alloc_len(end);
  prev = 0.0;

  for(i=0;i<end;i++) {
    e = SyExon_alloc();
    add_SyExonModel(out,e);
    if( i < start ) {
      e->exit_prob = 0.0;
      e->stay_prob = 0.0;
      e->prev_prob = prev;
      prev = 1.0;
    } else {
      e->exit_prob = exit;
      e->stay_prob = 0.0;
      e->prev_prob = prev;
      prev = 1.0-exit;
    }

    if( i+3 > end ) {
      e->stay_prob = final_stay;
    }
  }

  return out;
}
    
 
# line 80 "syexonmodel.dy"
SyExonScore * SyExonScore_from_SyExonModel(SyExonModel * sym)
{
  int i;

  SyExonScore * out;

  out = SyExonScore_alloc_len(sym->len);

  for(i=0;i<sym->len;i++) {
    add_SyExonScore(out,SyExonScoreUnit_from_SyExon(sym->exon[i]));
  }

  return out;
}


# line 96 "syexonmodel.dy"
SyExonScoreUnit * SyExonScoreUnit_from_SyExon(SyExon * sye)
{
  SyExonScoreUnit * out;

  out = SyExonScoreUnit_alloc();

  out->exit_score = Probability2Score(sye->exit_prob);
  out->stay_score = Probability2Score(sye->stay_prob);
  
  return out;
}


# line 109 "syexonmodel.dy"
void dump_SyExonScore(SyExonScore * sc,FILE *ofp)
{
  int i;

  for(i=0;i<sc->len;i++) {
    fprintf(ofp,"%d exit %d stay %d\n",i,sc->exon[i]->exit_score,sc->exon[i]->stay_score);
  }

}



# line 101 "syexonmodel.c"
/* Function:  hard_link_SyExon(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SyExon *]
 *
 * Return [UNKN ]  Undocumented return value [SyExon *]
 *
 */
SyExon * hard_link_SyExon(SyExon * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a SyExon object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  SyExon_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SyExon *]
 *
 */
SyExon * SyExon_alloc(void) 
{
    SyExon * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(SyExon *) ckalloc (sizeof(SyExon))) == NULL)    {  
      warn("SyExon_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->exit_prob = 0.0;    
    out->stay_prob = 0.0;    
    out->prev_prob = 0.0;    


    return out;  
}    


/* Function:  free_SyExon(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SyExon *]
 *
 * Return [UNKN ]  Undocumented return value [SyExon *]
 *
 */
SyExon * free_SyExon(SyExon * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a SyExon obj. Should be trappable");    
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


/* Function:  hard_link_SyExonScoreUnit(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SyExonScoreUnit *]
 *
 * Return [UNKN ]  Undocumented return value [SyExonScoreUnit *]
 *
 */
SyExonScoreUnit * hard_link_SyExonScoreUnit(SyExonScoreUnit * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a SyExonScoreUnit object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  SyExonScoreUnit_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SyExonScoreUnit *]
 *
 */
SyExonScoreUnit * SyExonScoreUnit_alloc(void) 
{
    SyExonScoreUnit * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(SyExonScoreUnit *) ckalloc (sizeof(SyExonScoreUnit))) == NULL)  {  
      warn("SyExonScoreUnit_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->exit_score = 0; 
    out->stay_score = 0; 


    return out;  
}    


/* Function:  free_SyExonScoreUnit(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SyExonScoreUnit *]
 *
 * Return [UNKN ]  Undocumented return value [SyExonScoreUnit *]
 *
 */
SyExonScoreUnit * free_SyExonScoreUnit(SyExonScoreUnit * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a SyExonScoreUnit obj. Should be trappable");   
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


/* Function:  swap_SyExonModel(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_SyExonModel
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [SyExon **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_SyExonModel(SyExon ** list,int i,int j)  
{
    SyExon * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_SyExonModel(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_SyExonModel which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [SyExon **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_SyExonModel(SyExon ** list,int left,int right,int (*comp)(SyExon * ,SyExon * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_SyExonModel(list,left,(left+right)/2);  
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_SyExonModel (list,++last,i);    
      }  
    swap_SyExonModel (list,left,last);   
    qsort_SyExonModel(list,left,last-1,comp);    
    qsort_SyExonModel(list,last+1,right,comp);   
}    


/* Function:  sort_SyExonModel(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_SyExonModel
 *
 *
 * Arg:         obj [UNKN ] Object containing list [SyExonModel *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_SyExonModel(SyExonModel * obj,int (*comp)(SyExon *, SyExon *)) 
{
    qsort_SyExonModel(obj->exon,0,obj->len-1,comp);  
    return;  
}    


/* Function:  expand_SyExonModel(obj,len)
 *
 * Descrip:    Really an internal function for add_SyExonModel
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [SyExonModel *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_SyExonModel(SyExonModel * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_SyExonModel called with no need");    
      return TRUE;   
      }  


    if( (obj->exon = (SyExon ** ) ckrealloc (obj->exon,sizeof(SyExon *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_SyExonModel, returning FALSE");  
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_SyExonModel(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [SyExonModel *]
 * Arg:        add [OWNER] Object to add to the list [SyExon *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_SyExonModel(SyExonModel * obj,SyExon * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_SyExonModel(obj,obj->len + SyExonModelLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->exon[obj->len++]=add;   
    return TRUE; 
}    


/* Function:  flush_SyExonModel(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [SyExonModel *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_SyExonModel(SyExonModel * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->exon[i] != NULL)  {  
        free_SyExon(obj->exon[i]);   
        obj->exon[i] = NULL; 
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  SyExonModel_alloc_std(void)
 *
 * Descrip:    Equivalent to SyExonModel_alloc_len(SyExonModelLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SyExonModel *]
 *
 */
SyExonModel * SyExonModel_alloc_std(void) 
{
    return SyExonModel_alloc_len(SyExonModelLISTLENGTH); 
}    


/* Function:  SyExonModel_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [SyExonModel *]
 *
 */
SyExonModel * SyExonModel_alloc_len(int len) 
{
    SyExonModel * out;  /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = SyExonModel_alloc()) == NULL)  
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->exon = (SyExon ** ) ckcalloc (len,sizeof(SyExon *))) == NULL)   {  
      warn("Warning, ckcalloc failed in SyExonModel_alloc_len"); 
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_SyExonModel(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SyExonModel *]
 *
 * Return [UNKN ]  Undocumented return value [SyExonModel *]
 *
 */
SyExonModel * hard_link_SyExonModel(SyExonModel * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a SyExonModel object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  SyExonModel_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SyExonModel *]
 *
 */
SyExonModel * SyExonModel_alloc(void) 
{
    SyExonModel * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(SyExonModel *) ckalloc (sizeof(SyExonModel))) == NULL)  {  
      warn("SyExonModel_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->exon = NULL;    
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_SyExonModel(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SyExonModel *]
 *
 * Return [UNKN ]  Undocumented return value [SyExonModel *]
 *
 */
SyExonModel * free_SyExonModel(SyExonModel * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a SyExonModel obj. Should be trappable");   
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
    if( obj->exon != NULL)   {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->exon[i] != NULL)    
          free_SyExon(obj->exon[i]); 
        }  
      ckfree(obj->exon); 
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_SyExonScore(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_SyExonScore
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [SyExonScoreUnit **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_SyExonScore(SyExonScoreUnit ** list,int i,int j)  
{
    SyExonScoreUnit * temp;  
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_SyExonScore(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_SyExonScore which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [SyExonScoreUnit **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_SyExonScore(SyExonScoreUnit ** list,int left,int right,int (*comp)(SyExonScoreUnit * ,SyExonScoreUnit * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_SyExonScore(list,left,(left+right)/2);  
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_SyExonScore (list,++last,i);    
      }  
    swap_SyExonScore (list,left,last);   
    qsort_SyExonScore(list,left,last-1,comp);    
    qsort_SyExonScore(list,last+1,right,comp);   
}    


/* Function:  sort_SyExonScore(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_SyExonScore
 *
 *
 * Arg:         obj [UNKN ] Object containing list [SyExonScore *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_SyExonScore(SyExonScore * obj,int (*comp)(SyExonScoreUnit *, SyExonScoreUnit *)) 
{
    qsort_SyExonScore(obj->exon,0,obj->len-1,comp);  
    return;  
}    


/* Function:  expand_SyExonScore(obj,len)
 *
 * Descrip:    Really an internal function for add_SyExonScore
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [SyExonScore *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_SyExonScore(SyExonScore * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_SyExonScore called with no need");    
      return TRUE;   
      }  


    if( (obj->exon = (SyExonScoreUnit ** ) ckrealloc (obj->exon,sizeof(SyExonScoreUnit *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_SyExonScore, returning FALSE");  
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_SyExonScore(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [SyExonScore *]
 * Arg:        add [OWNER] Object to add to the list [SyExonScoreUnit *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_SyExonScore(SyExonScore * obj,SyExonScoreUnit * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_SyExonScore(obj,obj->len + SyExonScoreLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->exon[obj->len++]=add;   
    return TRUE; 
}    


/* Function:  flush_SyExonScore(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [SyExonScore *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_SyExonScore(SyExonScore * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->exon[i] != NULL)  {  
        free_SyExonScoreUnit(obj->exon[i]);  
        obj->exon[i] = NULL; 
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  SyExonScore_alloc_std(void)
 *
 * Descrip:    Equivalent to SyExonScore_alloc_len(SyExonScoreLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SyExonScore *]
 *
 */
SyExonScore * SyExonScore_alloc_std(void) 
{
    return SyExonScore_alloc_len(SyExonScoreLISTLENGTH); 
}    


/* Function:  SyExonScore_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [SyExonScore *]
 *
 */
SyExonScore * SyExonScore_alloc_len(int len) 
{
    SyExonScore * out;  /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = SyExonScore_alloc()) == NULL)  
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->exon = (SyExonScoreUnit ** ) ckcalloc (len,sizeof(SyExonScoreUnit *))) == NULL) {  
      warn("Warning, ckcalloc failed in SyExonScore_alloc_len"); 
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_SyExonScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SyExonScore *]
 *
 * Return [UNKN ]  Undocumented return value [SyExonScore *]
 *
 */
SyExonScore * hard_link_SyExonScore(SyExonScore * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a SyExonScore object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  SyExonScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SyExonScore *]
 *
 */
SyExonScore * SyExonScore_alloc(void) 
{
    SyExonScore * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(SyExonScore *) ckalloc (sizeof(SyExonScore))) == NULL)  {  
      warn("SyExonScore_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->exon = NULL;    
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_SyExonScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SyExonScore *]
 *
 * Return [UNKN ]  Undocumented return value [SyExonScore *]
 *
 */
SyExonScore * free_SyExonScore(SyExonScore * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a SyExonScore obj. Should be trappable");   
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
    if( obj->exon != NULL)   {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->exon[i] != NULL)    
          free_SyExonScoreUnit(obj->exon[i]);    
        }  
      ckfree(obj->exon); 
      }  


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
