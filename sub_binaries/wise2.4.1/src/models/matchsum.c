#ifdef _cplusplus
extern "C" {
#endif
#include "matchsum.h"

/* Function:  MatchSummarySet_from_AlnBlock_estwise(alb,qname,offset,target)
 *
 * Descrip:    Builds a MatchSummarySet from a
 *             EstWise alignment. this makes
 *             alot of assumptions about the labels
 *             setc in alb, so make sure it was a 
 *             estwise alignment  - however as you
 *             can notice this is exactly the same 
 *             labels as found in genewise set
 *
 *
 * Arg:           alb [UNKN ] Undocumented argument [AlnBlock *]
 * Arg:         qname [UNKN ] Undocumented argument [char *]
 * Arg:        offset [UNKN ] Undocumented argument [int]
 * Arg:        target [UNKN ] Undocumented argument [Sequence *]
 *
 * Return [UNKN ]  Undocumented return value [MatchSummarySet *]
 *
 */
# line 60 "matchsum.dy"
MatchSummarySet * MatchSummarySet_from_AlnBlock_estwise(AlnBlock * alb,char * qname,int offset,Sequence * target)
{
  return MatchSummarySet_from_AlnBlock_genewise(alb,qname,offset,target);
}

/* Function:  MatchSummarySet_from_AlnBlock_genewise(alb,qname,protoff,target)
 *
 * Descrip:    Builds a MatchSummarySet from a
 *             GeneWise alignment. this makes
 *             alot of assumptions about the labels
 *             setc in alb, so make sure it was a 
 *             genewise alignment 
 *
 *
 * Arg:            alb [UNKN ] Undocumented argument [AlnBlock *]
 * Arg:          qname [UNKN ] Undocumented argument [char *]
 * Arg:        protoff [UNKN ] Undocumented argument [int]
 * Arg:         target [UNKN ] Undocumented argument [Sequence *]
 *
 * Return [UNKN ]  Undocumented return value [MatchSummarySet *]
 *
 */
# line 72 "matchsum.dy"
MatchSummarySet * MatchSummarySet_from_AlnBlock_genewise(AlnBlock * alb,char * qname,int protoff,Sequence * target)
{
  MatchSummarySet * out;
  MatchSummary * ms;
  AlnColumn * alc;
  int len;

  out = MatchSummarySet_alloc_std();

  alc = alb->start;

  while( (ms = MatchSummary_from_AlnColumn_genewise(alc,&alc)) != NULL ) {
    ms->qstart+= protoff-1; /*** offset in old-style offset mode ***/
    ms->qend += protoff-1;
    if( target->offset > target->end ) {
      len = ms->tend - ms->tstart;
      /*    fprintf(stderr,"Look - length %d\n",len); */
      ms->tstart = target->offset-1 - ms->tstart;
      ms->tend   = ms->tstart - len;
    } else {
      ms->tstart += target->offset-1;
      ms->tend += target->offset-1;
    }

    ms->qname = stringalloc(qname);
    ms->tname = stringalloc(target->name);
    add_MatchSummarySet(out,ms);
  }

  return out;
}
  
/* Function:  MatchSummary_from_AlnColumn_genewise(alc,end)
 *
 * Descrip:    For making each alignment
 *
 *
 * Arg:        alc [UNKN ] Undocumented argument [AlnColumn *]
 * Arg:        end [UNKN ] Undocumented argument [AlnColumn **]
 *
 * Return [UNKN ]  Undocumented return value [MatchSummary *]
 *
 */
# line 108 "matchsum.dy"
MatchSummary * MatchSummary_from_AlnColumn_genewise(AlnColumn * alc,AlnColumn ** end)
{
  MatchSummary * out;
  int score;
  

  for(;alc != NULL && is_random_AlnColumn_genewise(alc) == TRUE;alc = alc->next) 
    ;

  if( alc == NULL ) {
    *end = NULL;
    return NULL;
  }

  out = MatchSummary_alloc();

  /* NB alignment start is one before the C ordinate start system */

  out->qstart = alc->alu[0]->start+1; /* in C coords */
  out->tstart = alc->alu[1]->start+1; /* in C coords */
  score  = alc->alu[0]->score[0]; 

  for(;alc->next != NULL && is_random_AlnColumn_genewise(alc->next) == FALSE;alc = alc->next) {
    score += alc->next->alu[0]->score[0];
    if( strcmp(alc->alu[1]->text_label,"SEQUENCE_DELETION") == 0 || strcmp(alc->alu[1]->text_label,"SEQUENCE_INSERTION") == 0) 
      out->tframeshift++;
    else if ( strcmp(alc->alu[1]->text_label,"CENTRAL_INTRON") == 0 ) 
      out->tintron++; /* NB, assuming we have collapsed central intron labels! */
  }

  out->qend = alc->alu[0]->end+1; /* in C coords */
  out->tend = alc->alu[1]->end+1;
  /*score += alc->alu[0]->score[0];*/
  out->bits = Score2Bits(score); /*** hmm... this could be bad news ***/

  *end = alc->next;

  return out;
}

/* Function:  show_MatchSummary_genewise_header(ofp)
 *
 * Descrip:    Shows a header (bits Query start etc) which matches
 *             the order of the show_MatchSummarySet_genewise
 *
 *
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 152 "matchsum.dy"
void show_MatchSummary_genewise_header(FILE * ofp)
{
  fprintf(ofp,"Bits   Query         start end Target      start end   idels introns\n");
}

/* Function:  show_MatchSummarySet_genewise(mss,ofp)
 *
 * Descrip:    shows Matchsummary for genewise
 *             results, ie with target indels and introns
 *
 *
 * Arg:        mss [UNKN ] Undocumented argument [MatchSummarySet *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 161 "matchsum.dy"
void show_MatchSummarySet_genewise(MatchSummarySet * mss,FILE * ofp)
{
  int i;

  for(i=0;i<mss->len;i++)
    show_MatchSummary_genewise(mss->ms[i],ofp);
}

/* Function:  show_MatchSummary_genewise(ms,ofp)
 *
 * Descrip:    shows a single match summary
 *
 *
 * Arg:         ms [UNKN ] Undocumented argument [MatchSummary *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 173 "matchsum.dy"
void show_MatchSummary_genewise(MatchSummary * ms,FILE * ofp)
{
  if( ms->tstart < ms->tend ) 
    fprintf(ofp,"%4.2f %-12s %4d %4d %-12s %4d %4d %4d %4d\n",
	    ms->bits,
	    ms->qname,
	    ms->qstart+1, /* map to bio coords */
	    ms->qend,
	    ms->tname,
	    ms->tstart+1, /* map to 'bio coords */
	    ms->tend,
	    ms->tframeshift,
	    ms->tintron);
  else 
    fprintf(ofp,"%4.2f %-12s %4d %4d %-12s %4d %4d %4d %4d\n",
	    ms->bits,
	    ms->qname,
	    ms->qstart+1, /* map to bio coords */
	    ms->qend,
	    ms->tname,
	    ms->tstart+1, /* map to 'bio coords */
	    ms->tend+2, /* map to bio coords */
	    ms->tframeshift,
	    ms->tintron);
}


/* Function:  show_MatchSummary_estwise_header(ofp)
 *
 * Descrip:    Shows a header (bits Query start etc) which matches
 *             the order of the show_MatchSummarySet_estwise
 *
 *
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 204 "matchsum.dy"
void show_MatchSummary_estwise_header(FILE * ofp)
{
  fprintf(ofp,"Bits Query         start end Target      start end   idels\n");
}

/* Function:  show_MatchSummarySet_estwise(mss,ofp)
 *
 * Descrip:    Shows an estwise match, ie with only the target with indels
 *
 *
 * Arg:        mss [UNKN ] Undocumented argument [MatchSummarySet *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 212 "matchsum.dy"
void show_MatchSummarySet_estwise(MatchSummarySet * mss,FILE * ofp)
{
  int i;

  for(i=0;i<mss->len;i++)
    show_MatchSummary_estwise(mss->ms[i],ofp);
}

 
/* Function:  show_MatchSummary_estwise(ms,ofp)
 *
 * Descrip:    shows a single estwise match
 *
 *
 * Arg:         ms [UNKN ] Undocumented argument [MatchSummary *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 225 "matchsum.dy"
void show_MatchSummary_estwise(MatchSummary * ms,FILE * ofp)
{
  fprintf(ofp,"%4.2f %-12s %4d %4d %-12s %4d %4d %4d\n",
	  ms->bits,
	  ms->qname,
	  ms->qstart+1, /* map to bio coords */
	  ms->qend,
	  ms->tname,
	  ms->tstart+1, /* map to 'bio coords */
	  ms->tend,
	  ms->tframeshift);
}

 
# line 243 "matchsum.c"
/* Function:  hard_link_MatchSummary(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [MatchSummary *]
 *
 * Return [UNKN ]  Undocumented return value [MatchSummary *]
 *
 */
MatchSummary * hard_link_MatchSummary(MatchSummary * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a MatchSummary object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  MatchSummary_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [MatchSummary *]
 *
 */
MatchSummary * MatchSummary_alloc(void) 
{
    MatchSummary * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(MatchSummary *) ckalloc (sizeof(MatchSummary))) == NULL)    {  
      warn("MatchSummary_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->bits = 0;   
    out->qname = NULL;   
    out->tname = NULL;   
    out->qstart = 0; 
    out->qend = 0;   
    out->tstart = 0; 
    out->tend = 0;   
    out->qintron = 0;    
    out->qframeshift = 0;    
    out->tintron = 0;    
    out->tframeshift = 0;    


    return out;  
}    


/* Function:  free_MatchSummary(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [MatchSummary *]
 *
 * Return [UNKN ]  Undocumented return value [MatchSummary *]
 *
 */
MatchSummary * free_MatchSummary(MatchSummary * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a MatchSummary obj. Should be trappable");  
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
    if( obj->qname != NULL)  
      ckfree(obj->qname);    
    if( obj->tname != NULL)  
      ckfree(obj->tname);    


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_MatchSummarySet(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_MatchSummarySet
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [MatchSummary **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_MatchSummarySet(MatchSummary ** list,int i,int j)  
{
    MatchSummary * temp; 
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_MatchSummarySet(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_MatchSummarySet which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [MatchSummary **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_MatchSummarySet(MatchSummary ** list,int left,int right,int (*comp)(MatchSummary * ,MatchSummary * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_MatchSummarySet(list,left,(left+right)/2);  
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_MatchSummarySet (list,++last,i);    
      }  
    swap_MatchSummarySet (list,left,last);   
    qsort_MatchSummarySet(list,left,last-1,comp);    
    qsort_MatchSummarySet(list,last+1,right,comp);   
}    


/* Function:  sort_MatchSummarySet(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_MatchSummarySet
 *
 *
 * Arg:         obj [UNKN ] Object containing list [MatchSummarySet *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_MatchSummarySet(MatchSummarySet * obj,int (*comp)(MatchSummary *, MatchSummary *)) 
{
    qsort_MatchSummarySet(obj->ms,0,obj->len-1,comp);    
    return;  
}    


/* Function:  expand_MatchSummarySet(obj,len)
 *
 * Descrip:    Really an internal function for add_MatchSummarySet
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [MatchSummarySet *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_MatchSummarySet(MatchSummarySet * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_MatchSummarySet called with no need");    
      return TRUE;   
      }  


    if( (obj->ms = (MatchSummary ** ) ckrealloc (obj->ms,sizeof(MatchSummary *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_MatchSummarySet, returning FALSE");  
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_MatchSummarySet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [MatchSummarySet *]
 * Arg:        add [OWNER] Object to add to the list [MatchSummary *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_MatchSummarySet(MatchSummarySet * obj,MatchSummary * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_MatchSummarySet(obj,obj->len + MatchSummarySetLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->ms[obj->len++]=add; 
    return TRUE; 
}    


/* Function:  flush_MatchSummarySet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [MatchSummarySet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_MatchSummarySet(MatchSummarySet * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->ms[i] != NULL)    {  
        free_MatchSummary(obj->ms[i]);   
        obj->ms[i] = NULL;   
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  MatchSummarySet_alloc_std(void)
 *
 * Descrip:    Equivalent to MatchSummarySet_alloc_len(MatchSummarySetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [MatchSummarySet *]
 *
 */
MatchSummarySet * MatchSummarySet_alloc_std(void) 
{
    return MatchSummarySet_alloc_len(MatchSummarySetLISTLENGTH); 
}    


/* Function:  MatchSummarySet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [MatchSummarySet *]
 *
 */
MatchSummarySet * MatchSummarySet_alloc_len(int len) 
{
    MatchSummarySet * out;  /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = MatchSummarySet_alloc()) == NULL)  
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->ms = (MatchSummary ** ) ckcalloc (len,sizeof(MatchSummary *))) == NULL) {  
      warn("Warning, ckcalloc failed in MatchSummarySet_alloc_len"); 
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_MatchSummarySet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [MatchSummarySet *]
 *
 * Return [UNKN ]  Undocumented return value [MatchSummarySet *]
 *
 */
MatchSummarySet * hard_link_MatchSummarySet(MatchSummarySet * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a MatchSummarySet object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  MatchSummarySet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [MatchSummarySet *]
 *
 */
MatchSummarySet * MatchSummarySet_alloc(void) 
{
    MatchSummarySet * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(MatchSummarySet *) ckalloc (sizeof(MatchSummarySet))) == NULL)  {  
      warn("MatchSummarySet_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->ms = NULL;  
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_MatchSummarySet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [MatchSummarySet *]
 *
 * Return [UNKN ]  Undocumented return value [MatchSummarySet *]
 *
 */
MatchSummarySet * free_MatchSummarySet(MatchSummarySet * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a MatchSummarySet obj. Should be trappable");   
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
    if( obj->ms != NULL) {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->ms[i] != NULL)  
          free_MatchSummary(obj->ms[i]); 
        }  
      ckfree(obj->ms);   
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  access_ms_MatchSummarySet(obj,i)
 *
 * Descrip:    Access members stored in the ms list
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the list [MatchSummarySet *]
 * Arg:          i [UNKN ] Position in the list [int]
 *
 * Return [SOFT ]  Element of the list [MatchSummary *]
 *
 */
MatchSummary * access_ms_MatchSummarySet(MatchSummarySet * obj,int i) 
{
    if( obj == NULL)     {  
      warn("In accessor function ms for object MatchSummarySet, got a NULL object"); 
      return NULL;   
      }  
    if( obj->len <= i )  {  
      warn("In accessor function ms for object MatchSummarySet, index %%d is greater than list length %%d",i,obj->len);  
      return NULL;   
      }  
    return obj->ms[i];   
}    


/* Function:  length_ms_MatchSummarySet(obj)
 *
 * Descrip:    discover the length of the list
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the list [MatchSummarySet *]
 *
 * Return [UNKN ]  length of the list [int]
 *
 */
int length_ms_MatchSummarySet(MatchSummarySet * obj) 
{
    if( obj == NULL)     {  
      warn("In length function ms for object MatchSummarySet, got a NULL object");   
      return -1;     
      }  
    return obj->len;     
}    


/* Function:  replace_bits_MatchSummary(obj,bits)
 *
 * Descrip:    Replace member variable bits
 *             For use principly by API functions
 *
 *
 * Arg:         obj [UNKN ] Object holding the variable [MatchSummary *]
 * Arg:        bits [OWNER] New value of the variable [double]
 *
 * Return [SOFT ]  member variable bits [boolean]
 *
 */
boolean replace_bits_MatchSummary(MatchSummary * obj,double bits) 
{
    if( obj == NULL)     {  
      warn("In replacement function bits for object MatchSummary, got a NULL object");   
      return FALSE;  
      }  
    obj->bits = bits;    
    return TRUE; 
}    


/* Function:  access_bits_MatchSummary(obj)
 *
 * Descrip:    Access member variable bits
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [MatchSummary *]
 *
 * Return [SOFT ]  member variable bits [double]
 *
 */
double access_bits_MatchSummary(MatchSummary * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function bits for object MatchSummary, got a NULL object");  
      return 0;  
      }  
    return obj->bits;    
}    


/* Function:  replace_qname_MatchSummary(obj,qname)
 *
 * Descrip:    Replace member variable qname
 *             For use principly by API functions
 *
 *
 * Arg:          obj [UNKN ] Object holding the variable [MatchSummary *]
 * Arg:        qname [OWNER] New value of the variable [char *]
 *
 * Return [SOFT ]  member variable qname [boolean]
 *
 */
boolean replace_qname_MatchSummary(MatchSummary * obj,char * qname) 
{
    if( obj == NULL)     {  
      warn("In replacement function qname for object MatchSummary, got a NULL object");  
      return FALSE;  
      }  
    obj->qname = qname;  
    return TRUE; 
}    


/* Function:  access_qname_MatchSummary(obj)
 *
 * Descrip:    Access member variable qname
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [MatchSummary *]
 *
 * Return [SOFT ]  member variable qname [char *]
 *
 */
char * access_qname_MatchSummary(MatchSummary * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function qname for object MatchSummary, got a NULL object"); 
      return NULL;   
      }  
    return obj->qname;   
}    


/* Function:  replace_tname_MatchSummary(obj,tname)
 *
 * Descrip:    Replace member variable tname
 *             For use principly by API functions
 *
 *
 * Arg:          obj [UNKN ] Object holding the variable [MatchSummary *]
 * Arg:        tname [OWNER] New value of the variable [char *]
 *
 * Return [SOFT ]  member variable tname [boolean]
 *
 */
boolean replace_tname_MatchSummary(MatchSummary * obj,char * tname) 
{
    if( obj == NULL)     {  
      warn("In replacement function tname for object MatchSummary, got a NULL object");  
      return FALSE;  
      }  
    obj->tname = tname;  
    return TRUE; 
}    


/* Function:  access_tname_MatchSummary(obj)
 *
 * Descrip:    Access member variable tname
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [MatchSummary *]
 *
 * Return [SOFT ]  member variable tname [char *]
 *
 */
char * access_tname_MatchSummary(MatchSummary * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function tname for object MatchSummary, got a NULL object"); 
      return NULL;   
      }  
    return obj->tname;   
}    


/* Function:  replace_qstart_MatchSummary(obj,qstart)
 *
 * Descrip:    Replace member variable qstart
 *             For use principly by API functions
 *
 *
 * Arg:           obj [UNKN ] Object holding the variable [MatchSummary *]
 * Arg:        qstart [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable qstart [boolean]
 *
 */
boolean replace_qstart_MatchSummary(MatchSummary * obj,int qstart) 
{
    if( obj == NULL)     {  
      warn("In replacement function qstart for object MatchSummary, got a NULL object"); 
      return FALSE;  
      }  
    obj->qstart = qstart;    
    return TRUE; 
}    


/* Function:  access_qstart_MatchSummary(obj)
 *
 * Descrip:    Access member variable qstart
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [MatchSummary *]
 *
 * Return [SOFT ]  member variable qstart [int]
 *
 */
int access_qstart_MatchSummary(MatchSummary * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function qstart for object MatchSummary, got a NULL object");    
      return 0;  
      }  
    return obj->qstart;  
}    


/* Function:  replace_qend_MatchSummary(obj,qend)
 *
 * Descrip:    Replace member variable qend
 *             For use principly by API functions
 *
 *
 * Arg:         obj [UNKN ] Object holding the variable [MatchSummary *]
 * Arg:        qend [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable qend [boolean]
 *
 */
boolean replace_qend_MatchSummary(MatchSummary * obj,int qend) 
{
    if( obj == NULL)     {  
      warn("In replacement function qend for object MatchSummary, got a NULL object");   
      return FALSE;  
      }  
    obj->qend = qend;    
    return TRUE; 
}    


/* Function:  access_qend_MatchSummary(obj)
 *
 * Descrip:    Access member variable qend
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [MatchSummary *]
 *
 * Return [SOFT ]  member variable qend [int]
 *
 */
int access_qend_MatchSummary(MatchSummary * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function qend for object MatchSummary, got a NULL object");  
      return 0;  
      }  
    return obj->qend;    
}    


/* Function:  replace_tstart_MatchSummary(obj,tstart)
 *
 * Descrip:    Replace member variable tstart
 *             For use principly by API functions
 *
 *
 * Arg:           obj [UNKN ] Object holding the variable [MatchSummary *]
 * Arg:        tstart [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable tstart [boolean]
 *
 */
boolean replace_tstart_MatchSummary(MatchSummary * obj,int tstart) 
{
    if( obj == NULL)     {  
      warn("In replacement function tstart for object MatchSummary, got a NULL object"); 
      return FALSE;  
      }  
    obj->tstart = tstart;    
    return TRUE; 
}    


/* Function:  access_tstart_MatchSummary(obj)
 *
 * Descrip:    Access member variable tstart
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [MatchSummary *]
 *
 * Return [SOFT ]  member variable tstart [int]
 *
 */
int access_tstart_MatchSummary(MatchSummary * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function tstart for object MatchSummary, got a NULL object");    
      return 0;  
      }  
    return obj->tstart;  
}    


/* Function:  replace_tend_MatchSummary(obj,tend)
 *
 * Descrip:    Replace member variable tend
 *             For use principly by API functions
 *
 *
 * Arg:         obj [UNKN ] Object holding the variable [MatchSummary *]
 * Arg:        tend [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable tend [boolean]
 *
 */
boolean replace_tend_MatchSummary(MatchSummary * obj,int tend) 
{
    if( obj == NULL)     {  
      warn("In replacement function tend for object MatchSummary, got a NULL object");   
      return FALSE;  
      }  
    obj->tend = tend;    
    return TRUE; 
}    


/* Function:  access_tend_MatchSummary(obj)
 *
 * Descrip:    Access member variable tend
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [MatchSummary *]
 *
 * Return [SOFT ]  member variable tend [int]
 *
 */
int access_tend_MatchSummary(MatchSummary * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function tend for object MatchSummary, got a NULL object");  
      return 0;  
      }  
    return obj->tend;    
}    


/* Function:  replace_qintron_MatchSummary(obj,qintron)
 *
 * Descrip:    Replace member variable qintron
 *             For use principly by API functions
 *
 *
 * Arg:            obj [UNKN ] Object holding the variable [MatchSummary *]
 * Arg:        qintron [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable qintron [boolean]
 *
 */
boolean replace_qintron_MatchSummary(MatchSummary * obj,int qintron) 
{
    if( obj == NULL)     {  
      warn("In replacement function qintron for object MatchSummary, got a NULL object");    
      return FALSE;  
      }  
    obj->qintron = qintron;  
    return TRUE; 
}    


/* Function:  access_qintron_MatchSummary(obj)
 *
 * Descrip:    Access member variable qintron
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [MatchSummary *]
 *
 * Return [SOFT ]  member variable qintron [int]
 *
 */
int access_qintron_MatchSummary(MatchSummary * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function qintron for object MatchSummary, got a NULL object");   
      return 0;  
      }  
    return obj->qintron;     
}    


/* Function:  replace_qframeshift_MatchSummary(obj,qframeshift)
 *
 * Descrip:    Replace member variable qframeshift
 *             For use principly by API functions
 *
 *
 * Arg:                obj [UNKN ] Object holding the variable [MatchSummary *]
 * Arg:        qframeshift [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable qframeshift [boolean]
 *
 */
boolean replace_qframeshift_MatchSummary(MatchSummary * obj,int qframeshift) 
{
    if( obj == NULL)     {  
      warn("In replacement function qframeshift for object MatchSummary, got a NULL object");    
      return FALSE;  
      }  
    obj->qframeshift = qframeshift;  
    return TRUE; 
}    


/* Function:  access_qframeshift_MatchSummary(obj)
 *
 * Descrip:    Access member variable qframeshift
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [MatchSummary *]
 *
 * Return [SOFT ]  member variable qframeshift [int]
 *
 */
int access_qframeshift_MatchSummary(MatchSummary * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function qframeshift for object MatchSummary, got a NULL object");   
      return 0;  
      }  
    return obj->qframeshift;     
}    


/* Function:  replace_tintron_MatchSummary(obj,tintron)
 *
 * Descrip:    Replace member variable tintron
 *             For use principly by API functions
 *
 *
 * Arg:            obj [UNKN ] Object holding the variable [MatchSummary *]
 * Arg:        tintron [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable tintron [boolean]
 *
 */
boolean replace_tintron_MatchSummary(MatchSummary * obj,int tintron) 
{
    if( obj == NULL)     {  
      warn("In replacement function tintron for object MatchSummary, got a NULL object");    
      return FALSE;  
      }  
    obj->tintron = tintron;  
    return TRUE; 
}    


/* Function:  access_tintron_MatchSummary(obj)
 *
 * Descrip:    Access member variable tintron
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [MatchSummary *]
 *
 * Return [SOFT ]  member variable tintron [int]
 *
 */
int access_tintron_MatchSummary(MatchSummary * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function tintron for object MatchSummary, got a NULL object");   
      return 0;  
      }  
    return obj->tintron;     
}    


/* Function:  replace_tframeshift_MatchSummary(obj,tframeshift)
 *
 * Descrip:    Replace member variable tframeshift
 *             For use principly by API functions
 *
 *
 * Arg:                obj [UNKN ] Object holding the variable [MatchSummary *]
 * Arg:        tframeshift [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable tframeshift [boolean]
 *
 */
boolean replace_tframeshift_MatchSummary(MatchSummary * obj,int tframeshift) 
{
    if( obj == NULL)     {  
      warn("In replacement function tframeshift for object MatchSummary, got a NULL object");    
      return FALSE;  
      }  
    obj->tframeshift = tframeshift;  
    return TRUE; 
}    


/* Function:  access_tframeshift_MatchSummary(obj)
 *
 * Descrip:    Access member variable tframeshift
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [MatchSummary *]
 *
 * Return [SOFT ]  member variable tframeshift [int]
 *
 */
int access_tframeshift_MatchSummary(MatchSummary * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function tframeshift for object MatchSummary, got a NULL object");   
      return 0;  
      }  
    return obj->tframeshift;     
}    



#ifdef _cplusplus
}
#endif
