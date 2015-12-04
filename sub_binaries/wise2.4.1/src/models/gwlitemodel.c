#ifdef _cplusplus
extern "C" {
#endif
#include "gwlitemodel.h"


/* Function:  GwLite_AlnBlock_surgery(alb)
 *
 * Descrip:    A pretty weird function. Takes an AlnBlock made by GwLite and
 *             performs the necessary surgery at the 3SS to make it look like
 *             the AlnBlocks produced by the other genewise models. This means
 *             it has to eat into the intron by 3 residues 
 *
 *
 * Arg:        alb [UNKN ] Undocumented argument [AlnBlock *]
 *
 */
# line 82 "gwlitemodel.dy"
void GwLite_AlnBlock_surgery(AlnBlock * alb)
{
  AlnColumn * alc;
  
  for(alc = alb->start;alc != NULL;alc = alc->next ) {
    if( strstartcmp(alc->alu[1]->text_label,"3SS") == 0 ) {
      alc->alu[1]->start -= 3;
    }
  }
}


/* Function:  GwLite_from_GeneWise(gwm)
 *
 * Descrip:    Builds a GwLite model from the GeneWise model
 *
 *
 * Arg:        gwm [UNKN ] Undocumented argument [GeneWise *]
 *
 * Return [UNKN ]  Undocumented return value [GwLite *]
 *
 */
# line 97 "gwlitemodel.dy"
GwLite * GwLite_from_GeneWise(GeneWise * gwm)
{
  GwLite * out;
  GwLiteSegment * seg;
  int i,c,codon125;
  
  out = GwLite_alloc_len(gwm->len);
  for(i=0;i<gwm->len;i++) {
    seg = GwLiteSegment_alloc();
    for(c=0;c<64;c++) {
      codon125 = codon_from_base4_codon(c);
      seg->match[c]= gwm->seg[i]->match[codon125];
      seg->insert[c]= gwm->seg[i]->insert[codon125];
    }
    
    seg->transition[GWL_MATCH2MATCH] = gwm->seg[i]->transition[GW_MATCH2MATCH];
    seg->transition[GWL_MATCH2INSERT] = gwm->seg[i]->transition[GW_MATCH2INSERT];
    seg->transition[GWL_MATCH2DELETE] = gwm->seg[i]->transition[GW_MATCH2DELETE];
    seg->transition[GWL_MATCH2END] = gwm->seg[i]->transition[GW_MATCH2END];

    seg->transition[GWL_INSERT2MATCH] = gwm->seg[i]->transition[GW_INSERT2MATCH];
    seg->transition[GWL_INSERT2INSERT] = gwm->seg[i]->transition[GW_INSERT2INSERT];
    seg->transition[GWL_INSERT2DELETE] = gwm->seg[i]->transition[GW_INSERT2DELETE];
    seg->transition[GWL_INSERT2END] = gwm->seg[i]->transition[GW_INSERT2END];

    seg->transition[GWL_DELETE2MATCH] = gwm->seg[i]->transition[GW_DELETE2MATCH];
    seg->transition[GWL_DELETE2INSERT] = gwm->seg[i]->transition[GW_DELETE2INSERT];
    seg->transition[GWL_DELETE2DELETE] = gwm->seg[i]->transition[GW_DELETE2DELETE];
    seg->transition[GWL_DELETE2END] = gwm->seg[i]->transition[GW_DELETE2END];

    seg->transition[GWL_START2MATCH] = gwm->seg[i]->transition[GW_START2MATCH];
    seg->transition[GWL_START2INSERT] = gwm->seg[i]->transition[GW_START2INSERT];
    seg->transition[GWL_START2DELETE] = gwm->seg[i]->transition[GW_START2DELETE];

    add_GwLite(out,seg);
  }

  out->name = stringalloc(gwm->name);
  return out;
}


    
/* Function:  GwLiteScore_from_GwLite(gwl)
 *
 * Descrip:    Makes a lite score from a lite probability basis
 *
 *
 * Arg:        gwl [UNKN ] Undocumented argument [GwLite *]
 *
 * Return [UNKN ]  Undocumented return value [GwLiteScore *]
 *
 */
# line 143 "gwlitemodel.dy"
GwLiteScore * GwLiteScore_from_GwLite(GwLite * gwl)
{
  GwLiteScore * out;
  GwLiteSegment * prev;
  int i;

  out = GwLiteScore_alloc_len(gwl->len);
  for(i=0,prev=NULL;i<gwl->len;i++) {
    add_GwLiteScore(out,GwLiteSegmentScore_from_GwLiteSegment(prev,gwl->seg[i]));
    prev = gwl->seg[i];
  }

  out->name = stringalloc(gwl->name);
  return out;
}
      

# line 160 "gwlitemodel.dy"
GwLiteSegmentScore * GwLiteSegmentScore_from_GwLiteSegment(GwLiteSegment *prev,GwLiteSegment * seg)
{
  GwLiteSegmentScore * out;
  
  out = GwLiteSegmentScore_alloc();

  Probability2Score_move(seg->match,out->match,GWL_EMISSION_LEN);
  Probability2Score_move(seg->insert,out->insert,GWL_EMISSION_LEN);

  if( prev != NULL ) {
    out->transition[GWL_MATCH2MATCH] = Probability2Score(prev->transition[GWL_MATCH2MATCH]);
    out->transition[GWL_INSERT2MATCH] = Probability2Score(prev->transition[GWL_INSERT2MATCH]);
    out->transition[GWL_DELETE2MATCH] = Probability2Score(prev->transition[GWL_DELETE2MATCH]);

    out->transition[GWL_MATCH2DELETE] = Probability2Score(prev->transition[GWL_MATCH2DELETE]);
    out->transition[GWL_INSERT2DELETE] = Probability2Score(prev->transition[GWL_INSERT2DELETE]);
    out->transition[GWL_DELETE2DELETE] = Probability2Score(prev->transition[GWL_DELETE2DELETE]);
  } else {

    out->transition[GWL_MATCH2MATCH] = NEGI;
    out->transition[GWL_INSERT2MATCH] = NEGI;
    out->transition[GWL_DELETE2MATCH] = NEGI;

    out->transition[GWL_MATCH2DELETE] = NEGI;
    out->transition[GWL_INSERT2DELETE] = NEGI;
    out->transition[GWL_DELETE2DELETE] = NEGI;
  }

  out->transition[GWL_MATCH2INSERT] = Probability2Score(seg->transition[GWL_MATCH2INSERT]);
  out->transition[GWL_INSERT2INSERT] = Probability2Score(seg->transition[GWL_INSERT2INSERT]);
  out->transition[GWL_DELETE2INSERT] = Probability2Score(seg->transition[GWL_DELETE2INSERT]);
  
  out->transition[GWL_START2MATCH] = Probability2Score(seg->transition[GWL_START2MATCH]);
  out->transition[GWL_START2INSERT] = Probability2Score(seg->transition[GWL_START2INSERT]);
  out->transition[GWL_START2DELETE] = Probability2Score(seg->transition[GWL_START2DELETE]);
  
  out->transition[GWL_MATCH2END] = Probability2Score(seg->transition[GWL_MATCH2END]);
  out->transition[GWL_INSERT2END] = Probability2Score(seg->transition[GWL_INSERT2END]);
  out->transition[GWL_DELETE2END] = Probability2Score(seg->transition[GWL_DELETE2END]);

  return out;
}


# line 151 "gwlitemodel.c"
/* Function:  hard_link_GwLiteSegment(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GwLiteSegment *]
 *
 * Return [UNKN ]  Undocumented return value [GwLiteSegment *]
 *
 */
GwLiteSegment * hard_link_GwLiteSegment(GwLiteSegment * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a GwLiteSegment object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  GwLiteSegment_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GwLiteSegment *]
 *
 */
GwLiteSegment * GwLiteSegment_alloc(void) 
{
    GwLiteSegment * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GwLiteSegment *) ckalloc (sizeof(GwLiteSegment))) == NULL)  {  
      warn("GwLiteSegment_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    /* match[GWL_EMISSION_LEN] is an array: no default possible */ 
    /* insert[GWL_EMISSION_LEN] is an array: no default possible */ 
    /* transition[GWL_TRANSITION_LEN] is an array: no default possible */ 


    return out;  
}    


/* Function:  free_GwLiteSegment(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GwLiteSegment *]
 *
 * Return [UNKN ]  Undocumented return value [GwLiteSegment *]
 *
 */
GwLiteSegment * free_GwLiteSegment(GwLiteSegment * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a GwLiteSegment obj. Should be trappable"); 
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


/* Function:  swap_GwLite(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_GwLite
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [GwLiteSegment **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_GwLite(GwLiteSegment ** list,int i,int j)  
{
    GwLiteSegment * temp;    
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_GwLite(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_GwLite which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [GwLiteSegment **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_GwLite(GwLiteSegment ** list,int left,int right,int (*comp)(GwLiteSegment * ,GwLiteSegment * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_GwLite(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_GwLite (list,++last,i); 
      }  
    swap_GwLite (list,left,last);    
    qsort_GwLite(list,left,last-1,comp); 
    qsort_GwLite(list,last+1,right,comp);    
}    


/* Function:  sort_GwLite(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_GwLite
 *
 *
 * Arg:         obj [UNKN ] Object containing list [GwLite *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_GwLite(GwLite * obj,int (*comp)(GwLiteSegment *, GwLiteSegment *)) 
{
    qsort_GwLite(obj->seg,0,obj->len-1,comp);    
    return;  
}    


/* Function:  expand_GwLite(obj,len)
 *
 * Descrip:    Really an internal function for add_GwLite
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GwLite *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_GwLite(GwLite * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_GwLite called with no need"); 
      return TRUE;   
      }  


    if( (obj->seg = (GwLiteSegment ** ) ckrealloc (obj->seg,sizeof(GwLiteSegment *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_GwLite, returning FALSE");   
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_GwLite(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GwLite *]
 * Arg:        add [OWNER] Object to add to the list [GwLiteSegment *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_GwLite(GwLite * obj,GwLiteSegment * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_GwLite(obj,obj->len + GwLiteLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->seg[obj->len++]=add;    
    return TRUE; 
}    


/* Function:  flush_GwLite(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [GwLite *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_GwLite(GwLite * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->seg[i] != NULL)   {  
        free_GwLiteSegment(obj->seg[i]); 
        obj->seg[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  GwLite_alloc_std(void)
 *
 * Descrip:    Equivalent to GwLite_alloc_len(GwLiteLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GwLite *]
 *
 */
GwLite * GwLite_alloc_std(void) 
{
    return GwLite_alloc_len(GwLiteLISTLENGTH);   
}    


/* Function:  GwLite_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [GwLite *]
 *
 */
GwLite * GwLite_alloc_len(int len) 
{
    GwLite * out;   /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = GwLite_alloc()) == NULL)   
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->seg = (GwLiteSegment ** ) ckcalloc (len,sizeof(GwLiteSegment *))) == NULL)  {  
      warn("Warning, ckcalloc failed in GwLite_alloc_len");  
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_GwLite(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GwLite *]
 *
 * Return [UNKN ]  Undocumented return value [GwLite *]
 *
 */
GwLite * hard_link_GwLite(GwLite * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a GwLite object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  GwLite_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GwLite *]
 *
 */
GwLite * GwLite_alloc(void) 
{
    GwLite * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GwLite *) ckalloc (sizeof(GwLite))) == NULL)    {  
      warn("GwLite_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->seg = NULL; 
    out->len = out->maxlen = 0;  
    out->name = NULL;    


    return out;  
}    


/* Function:  free_GwLite(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GwLite *]
 *
 * Return [UNKN ]  Undocumented return value [GwLite *]
 *
 */
GwLite * free_GwLite(GwLite * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a GwLite obj. Should be trappable");    
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
    if( obj->seg != NULL)    {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->seg[i] != NULL) 
          free_GwLiteSegment(obj->seg[i]);   
        }  
      ckfree(obj->seg);  
      }  
    if( obj->name != NULL)   
      ckfree(obj->name);     


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_GwLiteSegmentScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GwLiteSegmentScore *]
 *
 * Return [UNKN ]  Undocumented return value [GwLiteSegmentScore *]
 *
 */
GwLiteSegmentScore * hard_link_GwLiteSegmentScore(GwLiteSegmentScore * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a GwLiteSegmentScore object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  GwLiteSegmentScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GwLiteSegmentScore *]
 *
 */
GwLiteSegmentScore * GwLiteSegmentScore_alloc(void) 
{
    GwLiteSegmentScore * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GwLiteSegmentScore *) ckalloc (sizeof(GwLiteSegmentScore))) == NULL)    {  
      warn("GwLiteSegmentScore_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    /* match[GWL_EMISSION_LEN] is an array: no default possible */ 
    /* insert[GWL_EMISSION_LEN] is an array: no default possible */ 
    /* transition[GWL_TRANSITION_LEN] is an array: no default possible */ 


    return out;  
}    


/* Function:  free_GwLiteSegmentScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GwLiteSegmentScore *]
 *
 * Return [UNKN ]  Undocumented return value [GwLiteSegmentScore *]
 *
 */
GwLiteSegmentScore * free_GwLiteSegmentScore(GwLiteSegmentScore * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a GwLiteSegmentScore obj. Should be trappable");    
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


/* Function:  swap_GwLiteScore(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_GwLiteScore
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [GwLiteSegmentScore **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_GwLiteScore(GwLiteSegmentScore ** list,int i,int j)  
{
    GwLiteSegmentScore * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_GwLiteScore(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_GwLiteScore which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [GwLiteSegmentScore **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_GwLiteScore(GwLiteSegmentScore ** list,int left,int right,int (*comp)(GwLiteSegmentScore * ,GwLiteSegmentScore * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_GwLiteScore(list,left,(left+right)/2);  
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_GwLiteScore (list,++last,i);    
      }  
    swap_GwLiteScore (list,left,last);   
    qsort_GwLiteScore(list,left,last-1,comp);    
    qsort_GwLiteScore(list,last+1,right,comp);   
}    


/* Function:  sort_GwLiteScore(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_GwLiteScore
 *
 *
 * Arg:         obj [UNKN ] Object containing list [GwLiteScore *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_GwLiteScore(GwLiteScore * obj,int (*comp)(GwLiteSegmentScore *, GwLiteSegmentScore *)) 
{
    qsort_GwLiteScore(obj->seg,0,obj->len-1,comp);   
    return;  
}    


/* Function:  expand_GwLiteScore(obj,len)
 *
 * Descrip:    Really an internal function for add_GwLiteScore
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GwLiteScore *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_GwLiteScore(GwLiteScore * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_GwLiteScore called with no need");    
      return TRUE;   
      }  


    if( (obj->seg = (GwLiteSegmentScore ** ) ckrealloc (obj->seg,sizeof(GwLiteSegmentScore *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_GwLiteScore, returning FALSE");  
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_GwLiteScore(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GwLiteScore *]
 * Arg:        add [OWNER] Object to add to the list [GwLiteSegmentScore *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_GwLiteScore(GwLiteScore * obj,GwLiteSegmentScore * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_GwLiteScore(obj,obj->len + GwLiteScoreLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->seg[obj->len++]=add;    
    return TRUE; 
}    


/* Function:  flush_GwLiteScore(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [GwLiteScore *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_GwLiteScore(GwLiteScore * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->seg[i] != NULL)   {  
        free_GwLiteSegmentScore(obj->seg[i]);    
        obj->seg[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  GwLiteScore_alloc_std(void)
 *
 * Descrip:    Equivalent to GwLiteScore_alloc_len(GwLiteScoreLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GwLiteScore *]
 *
 */
GwLiteScore * GwLiteScore_alloc_std(void) 
{
    return GwLiteScore_alloc_len(GwLiteScoreLISTLENGTH); 
}    


/* Function:  GwLiteScore_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [GwLiteScore *]
 *
 */
GwLiteScore * GwLiteScore_alloc_len(int len) 
{
    GwLiteScore * out;  /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = GwLiteScore_alloc()) == NULL)  
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->seg = (GwLiteSegmentScore ** ) ckcalloc (len,sizeof(GwLiteSegmentScore *))) == NULL)    {  
      warn("Warning, ckcalloc failed in GwLiteScore_alloc_len"); 
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_GwLiteScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GwLiteScore *]
 *
 * Return [UNKN ]  Undocumented return value [GwLiteScore *]
 *
 */
GwLiteScore * hard_link_GwLiteScore(GwLiteScore * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a GwLiteScore object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  GwLiteScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GwLiteScore *]
 *
 */
GwLiteScore * GwLiteScore_alloc(void) 
{
    GwLiteScore * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GwLiteScore *) ckalloc (sizeof(GwLiteScore))) == NULL)  {  
      warn("GwLiteScore_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->seg = NULL; 
    out->len = out->maxlen = 0;  
    out->name = NULL;    


    return out;  
}    


/* Function:  free_GwLiteScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GwLiteScore *]
 *
 * Return [UNKN ]  Undocumented return value [GwLiteScore *]
 *
 */
GwLiteScore * free_GwLiteScore(GwLiteScore * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a GwLiteScore obj. Should be trappable");   
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
    if( obj->seg != NULL)    {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->seg[i] != NULL) 
          free_GwLiteSegmentScore(obj->seg[i]);  
        }  
      ckfree(obj->seg);  
      }  
    if( obj->name != NULL)   
      ckfree(obj->name);     


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
