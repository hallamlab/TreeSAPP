#ifdef _cplusplus
extern "C" {
#endif
#include "dpenvelope.h"

/* Function:  show_DPEnvelope(dpe,ofp)
 *
 * Descrip:    shows structure. useful for debugging
 *
 *
 * Arg:        dpe [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 31 "dpenvelope.dy"
void show_DPEnvelope(DPEnvelope * dpe,FILE * ofp)
{
  int i;

  for(i=0;i<dpe->len;i++) {
    fprintf(ofp,"Unit [%d] %d-%d %d-%d\n",i,dpe->dpu[i]->starti,dpe->dpu[i]->startj,dpe->dpu[i]->starti+dpe->dpu[i]->height,dpe->dpu[i]->startj+dpe->dpu[i]->length);
  }

}

/* Function:  is_in_DPEnvelope(dpe,i,j)
 *
 * Descrip:    Tests whether this i,j position is allowed in the
 *             DPEnvelope
 *
 *
 * Arg:        dpe [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:          i [UNKN ] Undocumented argument [int]
 * Arg:          j [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 45 "dpenvelope.dy"
boolean is_in_DPEnvelope(DPEnvelope * dpe,int i,int j)
{
  int k;

  for(k=0;k<dpe->len;k++) {
    auto DPUnit * u;
    u = dpe->dpu[k];
    /* sorted by startj position */
    if( j < u->startj ) 
      return FALSE;
  

    switch (u->type) {

    case DPENV_RECT :
      if( i >= u->starti && j >= u->startj && i <= (u->starti+u->height) && j <= (u->startj+u->length) ) 
	return TRUE;
      else 
	break;
    case DPENV_DIAG :
      warn("Can't do diagonals yet.");
      return FALSE;
    default :
      warn("Bad DPUnit type put in. Yuk. Bad error... %d",u->type);
      return FALSE;
    }
  }

  /* bad error to get here! */
  return FALSE;
}

/* Function:  prepare_DPEnvelope(dpe)
 *
 * Descrip:    Should run this before using the DPEnvelope
 *
 *
 * Arg:        dpe [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 80 "dpenvelope.dy"
boolean prepare_DPEnvelope(DPEnvelope * dpe)
{
  int i;

  for(i=0;i<dpe->len;i++)
    if( dpe->dpu[i]->type != DPENV_RECT ) {
      warn("Bad envelope type %d",dpe->dpu[i]->type);
      return FALSE;
    }
  
  sort_DPEnvelope_by_startj(dpe);
  return TRUE;
}

/* Function:  sort_DPEnvelope_by_startj(dpe)
 *
 * Descrip:    Sorts by startj
 *
 *
 * Arg:        dpe [UNKN ] Undocumented argument [DPEnvelope *]
 *
 */
# line 97 "dpenvelope.dy"
void sort_DPEnvelope_by_startj(DPEnvelope * dpe)
{
  sort_DPEnvelope(dpe,compare_DPUnit_startj);
}

/* Function:  compare_DPUnit_startj(one,two)
 *
 * Descrip:    internal for sort by startj
 *
 *
 * Arg:        one [UNKN ] Undocumented argument [DPUnit *]
 * Arg:        two [UNKN ] Undocumented argument [DPUnit *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 106 "dpenvelope.dy"
int compare_DPUnit_startj(DPUnit * one,DPUnit * two)
{
  return one->startj  - two->startj;
}







# line 123 "dpenvelope.c"
/* Function:  hard_link_DPUnit(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DPUnit *]
 *
 * Return [UNKN ]  Undocumented return value [DPUnit *]
 *
 */
DPUnit * hard_link_DPUnit(DPUnit * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a DPUnit object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  DPUnit_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DPUnit *]
 *
 */
DPUnit * DPUnit_alloc(void) 
{
    DPUnit * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(DPUnit *) ckalloc (sizeof(DPUnit))) == NULL)    {  
      warn("DPUnit_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
    out->type = DPENV_RECT;  
    out->starti = 0; 
    out->startj = 0; 
    out->height = 0; 
    out->length = 0; 


    return out;  
}    


/* Function:  free_DPUnit(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DPUnit *]
 *
 * Return [UNKN ]  Undocumented return value [DPUnit *]
 *
 */
DPUnit * free_DPUnit(DPUnit * obj) 
{


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a DPUnit obj. Should be trappable");    
      return NULL;   
      }  


    if( obj->dynamite_hard_link > 1)     {  
      obj->dynamite_hard_link--; 
      return NULL;   
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_DPEnvelope(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_DPEnvelope
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [DPUnit **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_DPEnvelope(DPUnit ** list,int i,int j)  
{
    DPUnit * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_DPEnvelope(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_DPEnvelope which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [DPUnit **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_DPEnvelope(DPUnit ** list,int left,int right,int (*comp)(DPUnit * ,DPUnit * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_DPEnvelope(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_DPEnvelope (list,++last,i); 
      }  
    swap_DPEnvelope (list,left,last);    
    qsort_DPEnvelope(list,left,last-1,comp); 
    qsort_DPEnvelope(list,last+1,right,comp);    
}    


/* Function:  sort_DPEnvelope(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_DPEnvelope
 *
 *
 * Arg:         obj [UNKN ] Object containing list [DPEnvelope *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_DPEnvelope(DPEnvelope * obj,int (*comp)(DPUnit *, DPUnit *)) 
{
    qsort_DPEnvelope(obj->dpu,0,obj->len-1,comp);    
    return;  
}    


/* Function:  expand_DPEnvelope(obj,len)
 *
 * Descrip:    Really an internal function for add_DPEnvelope
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [DPEnvelope *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_DPEnvelope(DPEnvelope * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_DPEnvelope called with no need"); 
      return TRUE;   
      }  


    if( (obj->dpu = (DPUnit ** ) ckrealloc (obj->dpu,sizeof(DPUnit *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_DPEnvelope, returning FALSE");   
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_DPEnvelope(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [DPEnvelope *]
 * Arg:        add [OWNER] Object to add to the list [DPUnit *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_DPEnvelope(DPEnvelope * obj,DPUnit * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_DPEnvelope(obj,obj->len + DPEnvelopeLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->dpu[obj->len++]=add;    
    return TRUE; 
}    


/* Function:  flush_DPEnvelope(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_DPEnvelope(DPEnvelope * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->dpu[i] != NULL)   {  
        free_DPUnit(obj->dpu[i]);    
        obj->dpu[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  DPEnvelope_alloc_std(void)
 *
 * Descrip:    Equivalent to DPEnvelope_alloc_len(DPEnvelopeLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DPEnvelope *]
 *
 */
DPEnvelope * DPEnvelope_alloc_std(void) 
{
    return DPEnvelope_alloc_len(DPEnvelopeLISTLENGTH);   
}    


/* Function:  DPEnvelope_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [DPEnvelope *]
 *
 */
DPEnvelope * DPEnvelope_alloc_len(int len) 
{
    DPEnvelope * out;   /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = DPEnvelope_alloc()) == NULL)   
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->dpu = (DPUnit ** ) ckcalloc (len,sizeof(DPUnit *))) == NULL)    {  
      warn("Warning, ckcalloc failed in DPEnvelope_alloc_len");  
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_DPEnvelope(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [DPEnvelope *]
 *
 */
DPEnvelope * hard_link_DPEnvelope(DPEnvelope * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a DPEnvelope object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  DPEnvelope_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DPEnvelope *]
 *
 */
DPEnvelope * DPEnvelope_alloc(void) 
{
    DPEnvelope * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(DPEnvelope *) ckalloc (sizeof(DPEnvelope))) == NULL)    {  
      warn("DPEnvelope_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
    out->dpu = NULL; 
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_DPEnvelope(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [DPEnvelope *]
 *
 */
DPEnvelope * free_DPEnvelope(DPEnvelope * obj) 
{
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a DPEnvelope obj. Should be trappable");    
      return NULL;   
      }  


    if( obj->dynamite_hard_link > 1)     {  
      obj->dynamite_hard_link--; 
      return NULL;   
      }  
    if( obj->dpu != NULL)    {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->dpu[i] != NULL) 
          free_DPUnit(obj->dpu[i]);  
        }  
      ckfree(obj->dpu);  
      }  


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
