#ifdef _cplusplus
extern "C" {
#endif
#include "genome_evidence.h"


/* Function:  free_GenomeEvidenceUnit(obj)
 *
 * Descrip:    Specialised deconstructor. Ensures the 
 *             data structures are freed
 *
 *
 * Arg:        obj [UNKN ] Undocumented argument [GenomeEvidenceUnit *]
 *
 * Return [UNKN ]  Undocumented return value [GenomeEvidenceUnit *]
 *
 */
# line 57 "genome_evidence.dy"
GenomeEvidenceUnit * free_GenomeEvidenceUnit(GenomeEvidenceUnit * obj)
{
  (*(obj->geu_free))(obj->data);
  free(obj);
}


# line 23 "genome_evidence.c"
/* Function:  hard_link_GenomeEvidenceUnit(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GenomeEvidenceUnit *]
 *
 * Return [UNKN ]  Undocumented return value [GenomeEvidenceUnit *]
 *
 */
GenomeEvidenceUnit * hard_link_GenomeEvidenceUnit(GenomeEvidenceUnit * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a GenomeEvidenceUnit object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  GenomeEvidenceUnit_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GenomeEvidenceUnit *]
 *
 */
GenomeEvidenceUnit * GenomeEvidenceUnit_alloc(void) 
{
    GenomeEvidenceUnit * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GenomeEvidenceUnit *) ckalloc (sizeof(GenomeEvidenceUnit))) == NULL)    {  
      warn("GenomeEvidenceUnit_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->type = 0;   
    out->data = NULL;    
    out->cds_3SS = NULL; 
    out->cds_5SS = NULL; 
    out->utr_3SS = NULL; 
    out->utr_5SS = NULL; 
    out->cds_pot = NULL; 
    out->utr_pot = NULL; 
    out->cds_intron_pot = NULL;  
    out->utr_intron_pot = NULL;  
    out->frameshift_cds = NULL;  
    out->start_pot = NULL;   
    out->stop_pot = NULL;    
    out->utr3_end = NULL;    
    out->utr5_start = NULL;  
    out->geu_free = NULL;    


    return out;  
}    


/* Function:  swap_GenomeEvidenceSet(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_GenomeEvidenceSet
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [GenomeEvidenceUnit **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_GenomeEvidenceSet(GenomeEvidenceUnit ** list,int i,int j)  
{
    GenomeEvidenceUnit * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_GenomeEvidenceSet(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_GenomeEvidenceSet which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [GenomeEvidenceUnit **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_GenomeEvidenceSet(GenomeEvidenceUnit ** list,int left,int right,int (*comp)(GenomeEvidenceUnit * ,GenomeEvidenceUnit * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_GenomeEvidenceSet(list,left,(left+right)/2);    
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_GenomeEvidenceSet (list,++last,i);  
      }  
    swap_GenomeEvidenceSet (list,left,last); 
    qsort_GenomeEvidenceSet(list,left,last-1,comp);  
    qsort_GenomeEvidenceSet(list,last+1,right,comp); 
}    


/* Function:  sort_GenomeEvidenceSet(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_GenomeEvidenceSet
 *
 *
 * Arg:         obj [UNKN ] Object containing list [GenomeEvidenceSet *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_GenomeEvidenceSet(GenomeEvidenceSet * obj,int (*comp)(GenomeEvidenceUnit *, GenomeEvidenceUnit *)) 
{
    qsort_GenomeEvidenceSet(obj->geu,0,obj->len-1,comp); 
    return;  
}    


/* Function:  expand_GenomeEvidenceSet(obj,len)
 *
 * Descrip:    Really an internal function for add_GenomeEvidenceSet
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GenomeEvidenceSet *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_GenomeEvidenceSet(GenomeEvidenceSet * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_GenomeEvidenceSet called with no need");  
      return TRUE;   
      }  


    if( (obj->geu = (GenomeEvidenceUnit ** ) ckrealloc (obj->geu,sizeof(GenomeEvidenceUnit *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_GenomeEvidenceSet, returning FALSE");    
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_GenomeEvidenceSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GenomeEvidenceSet *]
 * Arg:        add [OWNER] Object to add to the list [GenomeEvidenceUnit *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_GenomeEvidenceSet(GenomeEvidenceSet * obj,GenomeEvidenceUnit * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_GenomeEvidenceSet(obj,obj->len + GenomeEvidenceSetLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->geu[obj->len++]=add;    
    return TRUE; 
}    


/* Function:  flush_GenomeEvidenceSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [GenomeEvidenceSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_GenomeEvidenceSet(GenomeEvidenceSet * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->geu[i] != NULL)   {  
        free_GenomeEvidenceUnit(obj->geu[i]);    
        obj->geu[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  GenomeEvidenceSet_alloc_std(void)
 *
 * Descrip:    Equivalent to GenomeEvidenceSet_alloc_len(GenomeEvidenceSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GenomeEvidenceSet *]
 *
 */
GenomeEvidenceSet * GenomeEvidenceSet_alloc_std(void) 
{
    return GenomeEvidenceSet_alloc_len(GenomeEvidenceSetLISTLENGTH); 
}    


/* Function:  GenomeEvidenceSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [GenomeEvidenceSet *]
 *
 */
GenomeEvidenceSet * GenomeEvidenceSet_alloc_len(int len) 
{
    GenomeEvidenceSet * out;/* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = GenomeEvidenceSet_alloc()) == NULL)    
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->geu = (GenomeEvidenceUnit ** ) ckcalloc (len,sizeof(GenomeEvidenceUnit *))) == NULL)    {  
      warn("Warning, ckcalloc failed in GenomeEvidenceSet_alloc_len");   
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_GenomeEvidenceSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GenomeEvidenceSet *]
 *
 * Return [UNKN ]  Undocumented return value [GenomeEvidenceSet *]
 *
 */
GenomeEvidenceSet * hard_link_GenomeEvidenceSet(GenomeEvidenceSet * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a GenomeEvidenceSet object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  GenomeEvidenceSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GenomeEvidenceSet *]
 *
 */
GenomeEvidenceSet * GenomeEvidenceSet_alloc(void) 
{
    GenomeEvidenceSet * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GenomeEvidenceSet *) ckalloc (sizeof(GenomeEvidenceSet))) == NULL)  {  
      warn("GenomeEvidenceSet_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->geu = NULL; 
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_GenomeEvidenceSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GenomeEvidenceSet *]
 *
 * Return [UNKN ]  Undocumented return value [GenomeEvidenceSet *]
 *
 */
GenomeEvidenceSet * free_GenomeEvidenceSet(GenomeEvidenceSet * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a GenomeEvidenceSet obj. Should be trappable"); 
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
    if( obj->geu != NULL)    {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->geu[i] != NULL) 
          free_GenomeEvidenceUnit(obj->geu[i]);  
        }  
      ckfree(obj->geu);  
      }  


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
