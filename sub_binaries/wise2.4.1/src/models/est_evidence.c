#ifdef _cplusplus
extern "C" {
#endif
#include "est_evidence.h"

# line 33 "est_evidence.dy"
int indicate_intron_used(GenomeEvidenceSet * set,AlnBlock * alb)
{
  AlnColumn * alc;
  int i;
  int j;
  int has_used = 0;

  for(alc=alb->start;alc!=NULL;alc=alc->next ) {
    if( strstartcmp(alc->alu[1]->text_label,"5SS") == 0 ) {
      for(i=0;i<set->len;i++) {
	auto EstEvidence * evi;
	evi = (EstEvidence*) set->geu[i]->data;
	for(j=0;j<evi->len;j++) {
	  if( abs(evi->exon[i]->end - alc->alu[1]->start) < 4 ) {
	    /* this is used */
	    if( evi->exon[i]->used != 1 ) {
	      evi->exon[i]->used = 1;
	      has_used = 1;
	    }
	  }
	}
      }
    }
  }

  return has_used;
}
      
  
# line 62 "est_evidence.dy"
GenomeEvidenceSet * read_est_evidence(FILE * ifp,CodonTable * ct)
{
  char buffer[MAXLINE];
  GenomeEvidenceUnit * geu;
  EstEvidence * evi;
  EstExon * exon;
  GenomeEvidenceSet * ges;
  EstIndel * indel;

  assert(ct);
  assert(ifp);
  ges = GenomeEvidenceSet_alloc_std();
  evi = EstEvidence_alloc_std();
  evi->ct = hard_link_CodonTable(ct);
  
  while( fgets(buffer,MAXLINE,ifp) != NULL ) {
    if( buffer[0] == '#' ) {
      continue;
    }
    if( strstartcmp(buffer,"//") == 0 ) {
      geu = new_est_GenomeEvidenceUnit(evi);
      add_GenomeEvidenceSet(ges,geu);
      evi = EstEvidence_alloc_std();
      evi->ct = hard_link_CodonTable(ct);
      continue;
    } 
    if( strstartcmp(buffer,"exon") == 0 ) {
      exon = EstExon_alloc();
      exon->intron_3_score = 0;
      if( sscanf(buffer,"exon %d %d %d",&exon->start,&exon->end,&exon->intron_3_score) < 2 ) {
	fatal("Unable to read exon line as evidence [%s]");
      }

      exon->start--;
      exon->end--;
      add_EstEvidence(evi,exon);
    } else if( strstartcmp(buffer,"cds") == 0 ) {
      exon = EstExon_alloc();
      sscanf(buffer,"cds %d %d %d",&exon->start,&exon->end,&exon->phase);
      exon->start--;
      exon->end--;
      if( exon->phase > 2 || exon->phase < 0 ) {
	fprintf(stderr,"Exon has a non clear phase - %d\n",exon->phase);
	return NULL;
      }
      exon->is_coding = TRUE;
      add_EstEvidence(evi,exon);
    } else if ( strstartcmp(buffer,"indel") == 0 ) {
      indel = EstIndel_alloc();
      sscanf(buffer,"indel %d %d",&indel->start,&indel->end);
      indel->start--;
      indel->end--;
      add_indel_EstEvidence(evi,indel);
    } else {
      fprintf(stderr,"Unable to read as est evidence - %s",buffer);
    }

  }
  if( evi->len > 0 ) {
    geu = new_est_GenomeEvidenceUnit(evi);
    add_GenomeEvidenceSet(ges,geu);
  }
  return ges;
}


# line 128 "est_evidence.dy"
GenomeEvidenceUnit * new_est_GenomeEvidenceUnit(EstEvidence * evi)
{
  GenomeEvidenceUnit * in;

  in = GenomeEvidenceUnit_alloc();

  in->cds_3SS = est_cds_3SS;
  in->cds_5SS = est_cds_5SS;
  in->utr_3SS = est_3ss;
  in->utr_5SS = est_5ss;
  in->cds_pot = est_cds_pot;
  in->utr_pot = est_utr_pot;
  in->cds_intron_pot = est_intron_pot;
  in->utr_intron_pot = est_intron_pot;
  in->geu_free = free_EstEvidence;
  in->frameshift_cds = est_cds_frameshift;
  in->stop_pot = est_stop_pot;
  in->start_pot = est_start_pot;
  in->utr3_end = est_utr3_end;
  in->utr5_start = est_utr5_start;
  in->data = (void*) evi;

  return in;
}


# line 154 "est_evidence.dy"
int est_utr5_start(void * data,ComplexSequence *seq,int jposition)
{
  EstEvidence * est;

  est = (EstEvidence *)data;
  if( est->exon[0]->start == jposition ) {
    return 0;
  } else {
    return -80000;
  }

}

# line 167 "est_evidence.dy"
int est_utr3_end(void * data,ComplexSequence *seq,int jposition)
{
  EstEvidence * est;

  est = (EstEvidence *)data;
  if( est->exon[est->len-1]->end == jposition ) {
    return 0;
  } else {
    return -80000;
  }

}

# line 180 "est_evidence.dy"
int est_start_pot(void * data,ComplexSequence *seq,int jposition)
{
  EstEvidence * est;
  int i;
  int codon;
  int atg = (BASE_A*25+BASE_T*5+BASE_G);

  est = (EstEvidence *)data;

  codon = CSEQ_GENOMIC_CODON(seq,jposition);
  if( is_stop_codon(codon,est->ct) ) {
    return -10000;
  } else if( codon == atg ) {
    return 1200;
  } else {
    return 0;
  }
}


# line 200 "est_evidence.dy"
int est_stop_pot(void * data,ComplexSequence *seq,int jposition)
{
  EstEvidence * est;
  int i;
  est = (EstEvidence *)data;

  if( is_stop_codon(CSEQ_GENOMIC_CODON(seq,jposition),est->ct) ) {
    return 100;
  } else {
    return -1000;
  }
}


# line 214 "est_evidence.dy"
int est_cds_frameshift(void * data,ComplexSequence * seq,int jposition,int jump)
{
  EstEvidence * est;
  int i;
  est = (EstEvidence *)data;


  for(i=1;i<est->indel_len;i++) {
    if( jposition >= est->indel[i]->start && jposition <= est->indel[i]->end ) {
      return 0;
    }
  }

  return -10000;
}


# line 231 "est_evidence.dy"
int est_cds_3SS(void * data,ComplexSequence *seq,int jposition,int phase)
{
  switch(phase) {
  case 0 : return est_3ss(data,seq,jposition);
  case 1 : return est_3ss(data,seq,jposition);
  case 2 : return est_3ss(data,seq,jposition);
  default : return -100000;
  }
}

# line 241 "est_evidence.dy"
int est_cds_5SS(void * data,ComplexSequence *seq,int jposition,int phase)
{
  switch(phase) {
  case 0 : return est_5ss(data,seq,jposition);
  case 1 : return est_5ss(data,seq,jposition);
  case 2 : return est_5ss(data,seq,jposition);
  default : return -100000;
  }
}


# line 252 "est_evidence.dy"
int est_intron_pot(void * data,ComplexSequence *seq,int jposition)
{
  EstEvidence * est;
  int i;
  est = (EstEvidence *)data;


  for(i=1;i<est->len;i++) {
    if( (est->exon[i-1]->end <= jposition) && (jposition <= est->exon[i]->start) ) {
      return 0;
    }
    if( (est->exon[i]->start < jposition) && (jposition < est->exon[i]->end) ) {
      return -1000;
    }
  }

  return -10;
}


# line 272 "est_evidence.dy"
int est_cds_pot(void * data,ComplexSequence *seq,int jposition)
{
  int i;
  EstEvidence * est;
  int relative_frame;

  est = (EstEvidence *)data;


  for(i=0;i<est->len;i++) {
    if( est->exon[i]->start <= jposition && jposition <= est->exon[i]->end ) {
      if( is_stop_codon(CSEQ_GENOMIC_CODON(seq,jposition),est->ct) ) {
	return -1000000;
      } else {
	if( est->exon[i]->is_coding == TRUE ) {
	  /* phase calculation. difference between start and position */

	  /* more complex than it looks due to convention of where a 
	     codon lies and the phase convention */

	  relative_frame = (jposition-est->exon[i]->start)%3;
	  if( relative_frame == 2 && est->exon[i]->phase == 0 ) {
	    return 125;
	  } else if ( relative_frame == 1 && est->exon[i]->phase == 1 ) {
	    return 125;
	  } else if ( relative_frame == 0 && est->exon[i]->phase == 2) {
	    return 125;
	  }
	} else {
	  /* not coding exon - return 45 */
	  return 80;
	}
      }
    }
  }

  /* we have to return same as stop codon penalty, otherwise
     we can just dodge stop codons using evidence lines */

  return -1000000;
}

# line 314 "est_evidence.dy"
int est_3ss(void * data,ComplexSequence *seq,int jposition)
{
  int i;
  EstEvidence * est;
  est = (EstEvidence *)data;

  if( jposition == 0 ) {
    return -10000;
  }
  
	
  for(i=0;i<est->len;i++) {
    if( jposition+1 == est->exon[i]->start ) {
      return 50;
    }

    if( abs(jposition+1 - est->exon[i]->start) < est->in_smell && seq->seq->seq[jposition] == 'G' && seq->seq->seq[jposition-1] == 'A' ) {
      return -400;
    }
    

  }
  
  return -10000;
}


# line 341 "est_evidence.dy"
int est_5ss(void * data,ComplexSequence * seq,int jposition)
{
  int i;
  EstEvidence * est;
  est = (EstEvidence *)data;

  if( jposition == 0 || jposition >= seq->seq->len+2 ) {
    return -10000;
  }
  
  for(i=0;i<est->len;i++) {
    if( jposition-1 == est->exon[i]->end ) {
      if( est->exon[i]->used == 0 ) {
	return 500+est->exon[i]->intron_3_score;
      } else {
	return 10+est->exon[i]->intron_3_score;
      }
    }
    
    if( abs(jposition-1 - est->exon[i]->end) < est->in_smell && seq->seq->seq[jposition] == 'G' && seq->seq->seq[jposition+1] == 'T' ) {
      return -400;
    }
    
  }
  
  return -10000;
}
 

# line 370 "est_evidence.dy"
int est_utr_pot(void * data,ComplexSequence *seq,int jposition)
{
  int i;
  EstEvidence * est;
  est = (EstEvidence *)data;

  for(i=0;i<est->len;i++) {
    if( est->exon[i]->start <= jposition && jposition <= est->exon[i]->end ) {
      return +10;
    }
  }

  return -10;
}






# line 377 "est_evidence.c"
/* Function:  hard_link_EstExon(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [EstExon *]
 *
 * Return [UNKN ]  Undocumented return value [EstExon *]
 *
 */
EstExon * hard_link_EstExon(EstExon * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a EstExon object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  EstExon_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [EstExon *]
 *
 */
EstExon * EstExon_alloc(void) 
{
    EstExon * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(EstExon *) ckalloc (sizeof(EstExon))) == NULL)  {  
      warn("EstExon_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->start = 0;  
    out->end = 0;    
    out->is_coding = FALSE;  
    out->phase = 0;  
    out->used = 0;   
    out->intron_3_score = 0; 


    return out;  
}    


/* Function:  free_EstExon(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [EstExon *]
 *
 * Return [UNKN ]  Undocumented return value [EstExon *]
 *
 */
EstExon * free_EstExon(EstExon * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a EstExon obj. Should be trappable");   
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


/* Function:  hard_link_EstIndel(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [EstIndel *]
 *
 * Return [UNKN ]  Undocumented return value [EstIndel *]
 *
 */
EstIndel * hard_link_EstIndel(EstIndel * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a EstIndel object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  EstIndel_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [EstIndel *]
 *
 */
EstIndel * EstIndel_alloc(void) 
{
    EstIndel * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(EstIndel *) ckalloc (sizeof(EstIndel))) == NULL)    {  
      warn("EstIndel_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->start = 0;  
    out->end = 0;    


    return out;  
}    


/* Function:  free_EstIndel(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [EstIndel *]
 *
 * Return [UNKN ]  Undocumented return value [EstIndel *]
 *
 */
EstIndel * free_EstIndel(EstIndel * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a EstIndel obj. Should be trappable");  
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


/* Function:  swap_EstEvidence(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_EstEvidence
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [EstExon    **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_EstEvidence(EstExon    ** list,int i,int j)  
{
    EstExon    * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_EstEvidence(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_EstEvidence which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [EstExon    **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_EstEvidence(EstExon    ** list,int left,int right,int (*comp)(EstExon    * ,EstExon    * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_EstEvidence(list,left,(left+right)/2);  
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_EstEvidence (list,++last,i);    
      }  
    swap_EstEvidence (list,left,last);   
    qsort_EstEvidence(list,left,last-1,comp);    
    qsort_EstEvidence(list,last+1,right,comp);   
}    


/* Function:  sort_EstEvidence(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_EstEvidence
 *
 *
 * Arg:         obj [UNKN ] Object containing list [EstEvidence *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_EstEvidence(EstEvidence * obj,int (*comp)(EstExon    *, EstExon    *)) 
{
    qsort_EstEvidence(obj->exon,0,obj->len-1,comp);  
    return;  
}    


/* Function:  expand_EstEvidence(obj,len)
 *
 * Descrip:    Really an internal function for add_EstEvidence
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [EstEvidence *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_EstEvidence(EstEvidence * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_EstEvidence called with no need");    
      return TRUE;   
      }  


    if( (obj->exon = (EstExon    ** ) ckrealloc (obj->exon,sizeof(EstExon    *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_EstEvidence, returning FALSE");  
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_EstEvidence(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [EstEvidence *]
 * Arg:        add [OWNER] Object to add to the list [EstExon    *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_EstEvidence(EstEvidence * obj,EstExon    * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_EstEvidence(obj,obj->len + EstEvidenceLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->exon[obj->len++]=add;   
    return TRUE; 
}    


/* Function:  flush_EstEvidence(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [EstEvidence *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_EstEvidence(EstEvidence * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->exon[i] != NULL)  {  
        free_EstExon(obj->exon[i]);  
        obj->exon[i] = NULL; 
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  swap_indel_EstEvidence(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_indel_EstEvidence
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [EstIndel   **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_indel_EstEvidence(EstIndel   ** list,int i,int j)  
{
    EstIndel   * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_indel_EstEvidence(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_indel_EstEvidence which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [EstIndel   **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_indel_EstEvidence(EstIndel   ** list,int left,int right,int (*comp)(EstIndel   * ,EstIndel   * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_indel_EstEvidence(list,left,(left+right)/2);    
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_indel_EstEvidence (list,++last,i);  
      }  
    swap_indel_EstEvidence (list,left,last); 
    qsort_indel_EstEvidence(list,left,last-1,comp);  
    qsort_indel_EstEvidence(list,last+1,right,comp); 
}    


/* Function:  sort_indel_EstEvidence(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_indel_EstEvidence
 *
 *
 * Arg:         obj [UNKN ] Object containing list [EstEvidence *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_indel_EstEvidence(EstEvidence * obj,int (*comp)(EstIndel   *, EstIndel   *)) 
{
    qsort_indel_EstEvidence(obj->indel,0,obj->indel_len-1,comp); 
    return;  
}    


/* Function:  expand_indel_EstEvidence(obj,len)
 *
 * Descrip:    Really an internal function for add_indel_EstEvidence
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [EstEvidence *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_indel_EstEvidence(EstEvidence * obj,int len) 
{


    if( obj->indel_maxlen > obj->indel_len )     {  
      warn("expand_EstEvidenceindel_ called with no need");  
      return TRUE;   
      }  


    if( (obj->indel = (EstIndel   ** ) ckrealloc (obj->indel,sizeof(EstIndel   *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_EstEvidence, returning FALSE");  
      return FALSE;  
      }  
    obj->indel_maxlen = len; 
    return TRUE; 
}    


/* Function:  add_indel_EstEvidence(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [EstEvidence *]
 * Arg:        add [OWNER] Object to add to the list [EstIndel   *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_indel_EstEvidence(EstEvidence * obj,EstIndel   * add) 
{
    if( obj->indel_len >= obj->indel_maxlen) {  
      if( expand_indel_EstEvidence(obj,obj->indel_len + EstEvidenceLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->indel[obj->indel_len++]=add;    
    return TRUE; 
}    


/* Function:  flush_indel_EstEvidence(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [EstEvidence *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_indel_EstEvidence(EstEvidence * obj) 
{
    int i;   


    for(i=0;i<obj->indel_len;i++)    { /*for i over list length*/ 
      if( obj->indel[i] != NULL) {  
        free_EstIndel(obj->indel[i]);    
        obj->indel[i] = NULL;    
        }  
      } /* end of for i over list length */ 


    obj->indel_len = 0;  
    return i;    
}    


/* Function:  EstEvidence_alloc_std(void)
 *
 * Descrip:    Equivalent to EstEvidence_alloc_len(EstEvidenceLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [EstEvidence *]
 *
 */
EstEvidence * EstEvidence_alloc_std(void) 
{
    return EstEvidence_alloc_len(EstEvidenceLISTLENGTH); 
}    


/* Function:  EstEvidence_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [EstEvidence *]
 *
 */
EstEvidence * EstEvidence_alloc_len(int len) 
{
    EstEvidence * out;  /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = EstEvidence_alloc()) == NULL)  
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->exon = (EstExon    ** ) ckcalloc (len,sizeof(EstExon    *))) == NULL)   {  
      warn("Warning, ckcalloc failed in EstEvidence_alloc_len"); 
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    if((out->indel = (EstIndel   ** ) ckcalloc (len,sizeof(EstIndel   *))) == NULL)  {  
      warn("Warning, ckcalloc failed in EstEvidence_alloc_len"); 
      return NULL;   
      }  
    out->indel_len = 0;  
    out->indel_maxlen = len; 


    return out;  
}    


/* Function:  hard_link_EstEvidence(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [EstEvidence *]
 *
 * Return [UNKN ]  Undocumented return value [EstEvidence *]
 *
 */
EstEvidence * hard_link_EstEvidence(EstEvidence * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a EstEvidence object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  EstEvidence_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [EstEvidence *]
 *
 */
EstEvidence * EstEvidence_alloc(void) 
{
    EstEvidence * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(EstEvidence *) ckalloc (sizeof(EstEvidence))) == NULL)  {  
      warn("EstEvidence_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->exon = NULL;    
    out->len = out->maxlen = 0;  
    out->indel = NULL;   
    out->indel_len = out->indel_maxlen = 0;  
    out->ct = NULL;  
    out->in_smell = 8;   


    return out;  
}    


/* Function:  free_EstEvidence(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [EstEvidence *]
 *
 * Return [UNKN ]  Undocumented return value [EstEvidence *]
 *
 */
EstEvidence * free_EstEvidence(EstEvidence * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a EstEvidence obj. Should be trappable");   
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
          free_EstExon(obj->exon[i]);    
        }  
      ckfree(obj->exon); 
      }  
    if( obj->indel != NULL)  {  
      for(i=0;i<obj->indel_len;i++)  {  
        if( obj->indel[i] != NULL)   
          free_EstIndel(obj->indel[i]);  
        }  
      ckfree(obj->indel);    
      }  
    if( obj->ct != NULL) 
      free_CodonTable(obj->ct);  


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
