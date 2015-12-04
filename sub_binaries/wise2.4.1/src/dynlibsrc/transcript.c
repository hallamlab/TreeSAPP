#ifdef _cplusplus
extern "C" {
#endif
#include "transcript.h"

/* Function:  copy_Transcript(t)
 *
 * Descrip:    Makes a completely new copy 
 *             of the transcript
 *
 *
 * Arg:        t [UNKN ] Undocumented argument [Transcript *]
 *
 * Return [UNKN ]  Undocumented return value [Transcript *]
 *
 */
# line 75 "transcript.dy"
Transcript * copy_Transcript(Transcript * t)
{
  Transcript * out;
  Exon * temp;
  int i;

  out = Transcript_alloc();

  for(i=0;i<t->ex_len;i++) {
    temp = Exon_alloc();
    temp->start = t->exon[i]->start;
    temp->end   = t->exon[i]->end;
    add_ex_Transcript(out,temp);
  }

  for(i=0;i<t->len;i++) {
    add_Transcript(out,copy_Translation(t->translation[i]));
  }

  return out;
}

/* Function:  get_cDNA_from_Transcript(trs)
 *
 * Descrip:    gets the cDNA associated with this transcript,
 *             if necessary, building it from the exon information
 *             provided.
 *
 *             returns a soft-linked object. If you want to ensure
 *             that this cDNA object remains in memory use
 *             /hard_link_cDNA on the object.
 *
 *
 * Arg:        trs [READ ] transcript to get cDNA from [Transcript *]
 *
 * Return [SOFT ]  cDNA of the transcript [cDNA *]
 *
 */
# line 109 "transcript.dy"
cDNA * get_cDNA_from_Transcript(Transcript * trs)
{
  Genomic * gn;
  Sequence * base;
  int i;
  char buffer[64];


  if( trs->cDNA != NULL) 
    return trs->cDNA;

  if( trs->parent == NULL ) {
    warn("Cannot get cDNA, as no parent Gene!");
    return NULL;
  }

  if ( (gn = get_Genomic_from_Gene(trs->parent)) == NULL  ) {
    warn("Cannot get cDNA, as cannot get Genomic sequence from Gene");
    return NULL;
  }

  base = Sequence_alloc();
  sprintf(buffer,"%s.sp",Genomic_name(gn));
  base->name = stringalloc(buffer);
  base->seq = ckcalloc(length_Transcript(trs)+1,sizeof(char));
  base->seq[0]='\0';
  

  for(i=0;i<trs->ex_len;i++) {
    strncat(base->seq,gn->baseseq->seq+trs->exon[i]->start,trs->exon[i]->end-trs->exon[i]->start);
  }
  make_len_type_Sequence(base);
  base->type = SEQUENCE_CDNA;
  trs->cDNA = cDNA_from_Sequence(base);

  return trs->cDNA;
}

/* Function:  length_Transcript(tr)
 *
 * Descrip:    returns the length by looking at the
 *             exon lengths
 *
 *
 * Arg:        tr [UNKN ] Undocumented argument [Transcript *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 151 "transcript.dy"
int length_Transcript(Transcript * tr)
{
  int len = 0;
  int i;

  for(i=0;i<tr->ex_len;i++) {
    len += tr->exon[i]->end - tr->exon[i]->start;
  }

  return len;
}
  
/* Function:  show_Transcript(tr,ofp)
 *
 * Descrip:    shows a transcript in vaguely human form
 *
 *
 * Arg:         tr [UNKN ] Undocumented argument [Transcript *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 166 "transcript.dy"
void show_Transcript(Transcript * tr,FILE * ofp)
{
  int i;	
  int k;
  SupportingFeature * sf;

  for(i=0;i<tr->ex_len;i++) {
    fprintf(ofp,"Exon %d-%d\n",tr->exon[i]->start,tr->exon[i]->end);
   
    /*
    for(k=0;k<tr->exon[i]->len;k++) {
      sf = tr->exon[i]->sf[k];
      fprintf(ofp," SF %d %d %d %d\n",sf->start,sf->end,sf->hstart,sf->hend);
    }
    */

  }

  for(i=0;i<tr->len;i++) {
    show_Translation(tr->translation[i],ofp);
  }

}



# line 146 "transcript.c"
/* Function:  hard_link_SupportingFeature(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SupportingFeature *]
 *
 * Return [UNKN ]  Undocumented return value [SupportingFeature *]
 *
 */
SupportingFeature * hard_link_SupportingFeature(SupportingFeature * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a SupportingFeature object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  SupportingFeature_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SupportingFeature *]
 *
 */
SupportingFeature * SupportingFeature_alloc(void) 
{
    SupportingFeature * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(SupportingFeature *) ckalloc (sizeof(SupportingFeature))) == NULL)  {  
      warn("SupportingFeature_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->start = 0;  
    out->end = 0;    
    out->hstart = 0; 
    out->hend = 0;   
    out->hstrand = 0;    
    out->hid = NULL; 


    return out;  
}    


/* Function:  free_SupportingFeature(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SupportingFeature *]
 *
 * Return [UNKN ]  Undocumented return value [SupportingFeature *]
 *
 */
SupportingFeature * free_SupportingFeature(SupportingFeature * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a SupportingFeature obj. Should be trappable"); 
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
    if( obj->hid != NULL)    
      ckfree(obj->hid);  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_Exon(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_Exon
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [SupportingFeature **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_Exon(SupportingFeature ** list,int i,int j)  
{
    SupportingFeature * temp;    
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_Exon(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_Exon which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [SupportingFeature **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_Exon(SupportingFeature ** list,int left,int right,int (*comp)(SupportingFeature * ,SupportingFeature * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_Exon(list,left,(left+right)/2); 
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_Exon (list,++last,i);   
      }  
    swap_Exon (list,left,last);  
    qsort_Exon(list,left,last-1,comp);   
    qsort_Exon(list,last+1,right,comp);  
}    


/* Function:  sort_Exon(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_Exon
 *
 *
 * Arg:         obj [UNKN ] Object containing list [Exon *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_Exon(Exon * obj,int (*comp)(SupportingFeature *, SupportingFeature *)) 
{
    qsort_Exon(obj->sf,0,obj->len-1,comp);   
    return;  
}    


/* Function:  expand_Exon(obj,len)
 *
 * Descrip:    Really an internal function for add_Exon
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Exon *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_Exon(Exon * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_Exon called with no need");   
      return TRUE;   
      }  


    if( (obj->sf = (SupportingFeature ** ) ckrealloc (obj->sf,sizeof(SupportingFeature *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_Exon, returning FALSE"); 
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_Exon(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Exon *]
 * Arg:        add [OWNER] Object to add to the list [SupportingFeature *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_Exon(Exon * obj,SupportingFeature * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_Exon(obj,obj->len + ExonLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->sf[obj->len++]=add; 
    return TRUE; 
}    


/* Function:  flush_Exon(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [Exon *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_Exon(Exon * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->sf[i] != NULL)    {  
        free_SupportingFeature(obj->sf[i]);  
        obj->sf[i] = NULL;   
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  Exon_alloc_std(void)
 *
 * Descrip:    Equivalent to Exon_alloc_len(ExonLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Exon *]
 *
 */
Exon * Exon_alloc_std(void) 
{
    return Exon_alloc_len(ExonLISTLENGTH);   
}    


/* Function:  Exon_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [Exon *]
 *
 */
Exon * Exon_alloc_len(int len) 
{
    Exon * out; /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = Exon_alloc()) == NULL) 
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->sf = (SupportingFeature ** ) ckcalloc (len,sizeof(SupportingFeature *))) == NULL)   {  
      warn("Warning, ckcalloc failed in Exon_alloc_len");    
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_Exon(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Exon *]
 *
 * Return [UNKN ]  Undocumented return value [Exon *]
 *
 */
Exon * hard_link_Exon(Exon * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a Exon object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  Exon_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Exon *]
 *
 */
Exon * Exon_alloc(void) 
{
    Exon * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(Exon *) ckalloc (sizeof(Exon))) == NULL)    {  
      warn("Exon_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->start = 0;  
    out->end = 0;    
    out->used = FALSE;   
    out->score = 0;  
    out->sf = NULL;  
    out->len = out->maxlen = 0;  
    out->phase = -1; 


    return out;  
}    


/* Function:  free_Exon(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Exon *]
 *
 * Return [UNKN ]  Undocumented return value [Exon *]
 *
 */
Exon * free_Exon(Exon * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a Exon obj. Should be trappable");  
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
    if( obj->sf != NULL) {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->sf[i] != NULL)  
          free_SupportingFeature(obj->sf[i]);    
        }  
      ckfree(obj->sf);   
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_ex_Transcript(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_ex_Transcript
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [Exon **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_ex_Transcript(Exon ** list,int i,int j)  
{
    Exon * temp; 
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_ex_Transcript(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_ex_Transcript which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [Exon **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_ex_Transcript(Exon ** list,int left,int right,int (*comp)(Exon * ,Exon * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_ex_Transcript(list,left,(left+right)/2);    
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_ex_Transcript (list,++last,i);  
      }  
    swap_ex_Transcript (list,left,last); 
    qsort_ex_Transcript(list,left,last-1,comp);  
    qsort_ex_Transcript(list,last+1,right,comp); 
}    


/* Function:  sort_ex_Transcript(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_ex_Transcript
 *
 *
 * Arg:         obj [UNKN ] Object containing list [Transcript *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_ex_Transcript(Transcript * obj,int (*comp)(Exon *, Exon *)) 
{
    qsort_ex_Transcript(obj->exon,0,obj->ex_len-1,comp); 
    return;  
}    


/* Function:  expand_ex_Transcript(obj,len)
 *
 * Descrip:    Really an internal function for add_ex_Transcript
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Transcript *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_ex_Transcript(Transcript * obj,int len) 
{


    if( obj->ex_maxlen > obj->ex_len )   {  
      warn("expand_Transcriptex_ called with no need");  
      return TRUE;   
      }  


    if( (obj->exon = (Exon ** ) ckrealloc (obj->exon,sizeof(Exon *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_Transcript, returning FALSE");   
      return FALSE;  
      }  
    obj->ex_maxlen = len;    
    return TRUE; 
}    


/* Function:  add_ex_Transcript(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Transcript *]
 * Arg:        add [OWNER] Object to add to the list [Exon *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_ex_Transcript(Transcript * obj,Exon * add) 
{
    if( obj->ex_len >= obj->ex_maxlen)   {  
      if( expand_ex_Transcript(obj,obj->ex_len + TranscriptLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->exon[obj->ex_len++]=add;    
    return TRUE; 
}    


/* Function:  flush_ex_Transcript(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [Transcript *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_ex_Transcript(Transcript * obj) 
{
    int i;   


    for(i=0;i<obj->ex_len;i++)   { /*for i over list length*/ 
      if( obj->exon[i] != NULL)  {  
        free_Exon(obj->exon[i]); 
        obj->exon[i] = NULL; 
        }  
      } /* end of for i over list length */ 


    obj->ex_len = 0; 
    return i;    
}    


/* Function:  swap_Transcript(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_Transcript
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [Translation **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_Transcript(Translation ** list,int i,int j)  
{
    Translation * temp;  
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_Transcript(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_Transcript which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [Translation **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_Transcript(Translation ** list,int left,int right,int (*comp)(Translation * ,Translation * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_Transcript(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_Transcript (list,++last,i); 
      }  
    swap_Transcript (list,left,last);    
    qsort_Transcript(list,left,last-1,comp); 
    qsort_Transcript(list,last+1,right,comp);    
}    


/* Function:  sort_Transcript(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_Transcript
 *
 *
 * Arg:         obj [UNKN ] Object containing list [Transcript *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_Transcript(Transcript * obj,int (*comp)(Translation *, Translation *)) 
{
    qsort_Transcript(obj->translation,0,obj->len-1,comp);    
    return;  
}    


/* Function:  expand_Transcript(obj,len)
 *
 * Descrip:    Really an internal function for add_Transcript
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Transcript *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_Transcript(Transcript * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_Transcript called with no need"); 
      return TRUE;   
      }  


    if( (obj->translation = (Translation ** ) ckrealloc (obj->translation,sizeof(Translation *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_Transcript, returning FALSE");   
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_Transcript(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Transcript *]
 * Arg:        add [OWNER] Object to add to the list [Translation *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_Transcript(Transcript * obj,Translation * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_Transcript(obj,obj->len + TranscriptLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->translation[obj->len++]=add;    
    return TRUE; 
}    


/* Function:  flush_Transcript(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [Transcript *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_Transcript(Transcript * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->translation[i] != NULL)   {  
        free_Translation(obj->translation[i]);   
        obj->translation[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  Transcript_alloc_std(void)
 *
 * Descrip:    Equivalent to Transcript_alloc_len(TranscriptLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Transcript *]
 *
 */
Transcript * Transcript_alloc_std(void) 
{
    return Transcript_alloc_len(TranscriptLISTLENGTH);   
}    


/* Function:  Transcript_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [Transcript *]
 *
 */
Transcript * Transcript_alloc_len(int len) 
{
    Transcript * out;   /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = Transcript_alloc()) == NULL)   
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->exon = (Exon ** ) ckcalloc (len,sizeof(Exon *))) == NULL)   {  
      warn("Warning, ckcalloc failed in Transcript_alloc_len");  
      return NULL;   
      }  
    out->ex_len = 0; 
    out->ex_maxlen = len;    


    if((out->translation = (Translation ** ) ckcalloc (len,sizeof(Translation *))) == NULL)  {  
      warn("Warning, ckcalloc failed in Transcript_alloc_len");  
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_Transcript(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Transcript *]
 *
 * Return [UNKN ]  Undocumented return value [Transcript *]
 *
 */
Transcript * hard_link_Transcript(Transcript * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a Transcript object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  Transcript_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Transcript *]
 *
 */
Transcript * Transcript_alloc(void) 
{
    Transcript * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(Transcript *) ckalloc (sizeof(Transcript))) == NULL)    {  
      warn("Transcript_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->exon = NULL;    
    out->ex_len = out->ex_maxlen = 0;    
    out->translation = NULL; 
    out->len = out->maxlen = 0;  
    out->cDNA = NULL;    


    return out;  
}    


/* Function:  free_Transcript(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Transcript *]
 *
 * Return [UNKN ]  Undocumented return value [Transcript *]
 *
 */
Transcript * free_Transcript(Transcript * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a Transcript obj. Should be trappable");    
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
      for(i=0;i<obj->ex_len;i++) {  
        if( obj->exon[i] != NULL)    
          free_Exon(obj->exon[i]);   
        }  
      ckfree(obj->exon); 
      }  
    /* obj->parent is linked in */ 
    if( obj->translation != NULL)    {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->translation[i] != NULL) 
          free_Translation(obj->translation[i]); 
        }  
      ckfree(obj->translation);  
      }  
    if( obj->cDNA != NULL)   
      free_cDNA(obj->cDNA);  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  replace_start_Exon(obj,start)
 *
 * Descrip:    Replace member variable start
 *             For use principly by API functions
 *
 *
 * Arg:          obj [UNKN ] Object holding the variable [Exon *]
 * Arg:        start [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable start [boolean]
 *
 */
boolean replace_start_Exon(Exon * obj,int start) 
{
    if( obj == NULL)     {  
      warn("In replacement function start for object Exon, got a NULL object");  
      return FALSE;  
      }  
    obj->start = start;  
    return TRUE; 
}    


/* Function:  access_start_Exon(obj)
 *
 * Descrip:    Access member variable start
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [Exon *]
 *
 * Return [SOFT ]  member variable start [int]
 *
 */
int access_start_Exon(Exon * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function start for object Exon, got a NULL object"); 
      return 0;  
      }  
    return obj->start;   
}    


/* Function:  replace_end_Exon(obj,end)
 *
 * Descrip:    Replace member variable end
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [Exon *]
 * Arg:        end [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable end [boolean]
 *
 */
boolean replace_end_Exon(Exon * obj,int end) 
{
    if( obj == NULL)     {  
      warn("In replacement function end for object Exon, got a NULL object");    
      return FALSE;  
      }  
    obj->end = end;  
    return TRUE; 
}    


/* Function:  access_end_Exon(obj)
 *
 * Descrip:    Access member variable end
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [Exon *]
 *
 * Return [SOFT ]  member variable end [int]
 *
 */
int access_end_Exon(Exon * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function end for object Exon, got a NULL object");   
      return 0;  
      }  
    return obj->end;     
}    


/* Function:  replace_used_Exon(obj,used)
 *
 * Descrip:    Replace member variable used
 *             For use principly by API functions
 *
 *
 * Arg:         obj [UNKN ] Object holding the variable [Exon *]
 * Arg:        used [OWNER] New value of the variable [boolean]
 *
 * Return [SOFT ]  member variable used [boolean]
 *
 */
boolean replace_used_Exon(Exon * obj,boolean used) 
{
    if( obj == NULL)     {  
      warn("In replacement function used for object Exon, got a NULL object");   
      return FALSE;  
      }  
    obj->used = used;    
    return TRUE; 
}    


/* Function:  access_used_Exon(obj)
 *
 * Descrip:    Access member variable used
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [Exon *]
 *
 * Return [SOFT ]  member variable used [boolean]
 *
 */
boolean access_used_Exon(Exon * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function used for object Exon, got a NULL object");  
      return FALSE;  
      }  
    return obj->used;    
}    


/* Function:  replace_score_Exon(obj,score)
 *
 * Descrip:    Replace member variable score
 *             For use principly by API functions
 *
 *
 * Arg:          obj [UNKN ] Object holding the variable [Exon *]
 * Arg:        score [OWNER] New value of the variable [double]
 *
 * Return [SOFT ]  member variable score [boolean]
 *
 */
boolean replace_score_Exon(Exon * obj,double score) 
{
    if( obj == NULL)     {  
      warn("In replacement function score for object Exon, got a NULL object");  
      return FALSE;  
      }  
    obj->score = score;  
    return TRUE; 
}    


/* Function:  access_score_Exon(obj)
 *
 * Descrip:    Access member variable score
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [Exon *]
 *
 * Return [SOFT ]  member variable score [double]
 *
 */
double access_score_Exon(Exon * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function score for object Exon, got a NULL object"); 
      return 0;  
      }  
    return obj->score;   
}    


/* Function:  access_sf_Exon(obj,i)
 *
 * Descrip:    Access members stored in the sf list
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the list [Exon *]
 * Arg:          i [UNKN ] Position in the list [int]
 *
 * Return [SOFT ]  Element of the list [SupportingFeature *]
 *
 */
SupportingFeature * access_sf_Exon(Exon * obj,int i) 
{
    if( obj == NULL)     {  
      warn("In accessor function sf for object Exon, got a NULL object");    
      return NULL;   
      }  
    if( obj->len <= i )  {  
      warn("In accessor function sf for object Exon, index %%d is greater than list length %%d",i,obj->len); 
      return NULL;   
      }  
    return obj->sf[i];   
}    


/* Function:  length_sf_Exon(obj)
 *
 * Descrip:    discover the length of the list
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the list [Exon *]
 *
 * Return [UNKN ]  length of the list [int]
 *
 */
int length_sf_Exon(Exon * obj) 
{
    if( obj == NULL)     {  
      warn("In length function sf for object Exon, got a NULL object");  
      return -1;     
      }  
    return obj->len;     
}    


/* Function:  replace_phase_Exon(obj,phase)
 *
 * Descrip:    Replace member variable phase
 *             For use principly by API functions
 *
 *
 * Arg:          obj [UNKN ] Object holding the variable [Exon *]
 * Arg:        phase [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable phase [boolean]
 *
 */
boolean replace_phase_Exon(Exon * obj,int phase) 
{
    if( obj == NULL)     {  
      warn("In replacement function phase for object Exon, got a NULL object");  
      return FALSE;  
      }  
    obj->phase = phase;  
    return TRUE; 
}    


/* Function:  access_phase_Exon(obj)
 *
 * Descrip:    Access member variable phase
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [Exon *]
 *
 * Return [SOFT ]  member variable phase [int]
 *
 */
int access_phase_Exon(Exon * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function phase for object Exon, got a NULL object"); 
      return 0;  
      }  
    return obj->phase;   
}    


/* Function:  access_exon_Transcript(obj,i)
 *
 * Descrip:    Access members stored in the exon list
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the list [Transcript *]
 * Arg:          i [UNKN ] Position in the list [int]
 *
 * Return [SOFT ]  Element of the list [Exon *]
 *
 */
Exon * access_exon_Transcript(Transcript * obj,int i) 
{
    if( obj == NULL)     {  
      warn("In accessor function exon for object Transcript, got a NULL object");    
      return NULL;   
      }  
    if( obj->ex_len <= i )   {  
      warn("In accessor function exon for object Transcript, index %%d is greater than list length %%d",i,obj->len); 
      return NULL;   
      }  
    return obj->exon[i];     
}    


/* Function:  length_exon_Transcript(obj)
 *
 * Descrip:    discover the length of the list
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the list [Transcript *]
 *
 * Return [UNKN ]  length of the list [int]
 *
 */
int length_exon_Transcript(Transcript * obj) 
{
    if( obj == NULL)     {  
      warn("In length function exon for object Transcript, got a NULL object");  
      return -1;     
      }  
    return obj->ex_len;  
}    


/* Function:  replace_parent_Transcript(obj,parent)
 *
 * Descrip:    Replace member variable parent
 *             For use principly by API functions
 *
 *
 * Arg:           obj [UNKN ] Object holding the variable [Transcript *]
 * Arg:        parent [OWNER] New value of the variable [Gene *]
 *
 * Return [SOFT ]  member variable parent [boolean]
 *
 */
boolean replace_parent_Transcript(Transcript * obj,Gene * parent) 
{
    if( obj == NULL)     {  
      warn("In replacement function parent for object Transcript, got a NULL object");   
      return FALSE;  
      }  
    obj->parent = parent;    
    return TRUE; 
}    


/* Function:  access_parent_Transcript(obj)
 *
 * Descrip:    Access member variable parent
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [Transcript *]
 *
 * Return [SOFT ]  member variable parent [Gene *]
 *
 */
Gene * access_parent_Transcript(Transcript * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function parent for object Transcript, got a NULL object");  
      return NULL;   
      }  
    return obj->parent;  
}    


/* Function:  access_translation_Transcript(obj,i)
 *
 * Descrip:    Access members stored in the translation list
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the list [Transcript *]
 * Arg:          i [UNKN ] Position in the list [int]
 *
 * Return [SOFT ]  Element of the list [Translation *]
 *
 */
Translation * access_translation_Transcript(Transcript * obj,int i) 
{
    if( obj == NULL)     {  
      warn("In accessor function translation for object Transcript, got a NULL object"); 
      return NULL;   
      }  
    if( obj->len <= i )  {  
      warn("In accessor function translation for object Transcript, index %%d is greater than list length %%d",i,obj->len);  
      return NULL;   
      }  
    return obj->translation[i];  
}    


/* Function:  length_translation_Transcript(obj)
 *
 * Descrip:    discover the length of the list
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the list [Transcript *]
 *
 * Return [UNKN ]  length of the list [int]
 *
 */
int length_translation_Transcript(Transcript * obj) 
{
    if( obj == NULL)     {  
      warn("In length function translation for object Transcript, got a NULL object");   
      return -1;     
      }  
    return obj->len;     
}    


/* Function:  replace_cDNA_Transcript(obj,cDNA)
 *
 * Descrip:    Replace member variable cDNA
 *             For use principly by API functions
 *
 *
 * Arg:         obj [UNKN ] Object holding the variable [Transcript *]
 * Arg:        cDNA [OWNER] New value of the variable [cDNA *]
 *
 * Return [SOFT ]  member variable cDNA [boolean]
 *
 */
boolean replace_cDNA_Transcript(Transcript * obj,cDNA * cDNA) 
{
    if( obj == NULL)     {  
      warn("In replacement function cDNA for object Transcript, got a NULL object"); 
      return FALSE;  
      }  
    obj->cDNA = cDNA;    
    return TRUE; 
}    


/* Function:  access_cDNA_Transcript(obj)
 *
 * Descrip:    Access member variable cDNA
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [Transcript *]
 *
 * Return [SOFT ]  member variable cDNA [cDNA *]
 *
 */
cDNA * access_cDNA_Transcript(Transcript * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function cDNA for object Transcript, got a NULL object");    
      return NULL;   
      }  
    return obj->cDNA;    
}    



#ifdef _cplusplus
}
#endif
