#ifdef _cplusplus
extern "C" {
#endif
#include "transregion.h"

# line 42 "transregion.dy"
void show_help_TransFactorRegionPara(FILE * ofp)
{
  fprintf(ofp,"TransFactor Region parameters\n");
  fprintf(ofp,"  -tf_convert [hmm/window] use HMM or window for conversion, default HMM\n");
  fprintf(ofp,"  -tf_density [0.2] minimum density of bases in a cluster for window\n");
  fprintf(ofp,"  -tf_window  [60]  minimum window size of region\n");
  fprintf(ofp,"  -tf_inprob  [0.35] probability of a covered base in a region for HMM\n");
  fprintf(ofp,"  -tf_outprob [0.2]  probability of a covered base outside a region for HMM\n");
  fprintf(ofp,"  -tf_entry   [0.000001] probablity of entering a region\n");
  fprintf(ofp,"  -tf_hmm_window [60] minimum window for HMM\n");
  fprintf(ofp,"  -tf_gc_region  [1.5] Expected ration of GC dinucleotides in regions vs out\n");

  return;
}

# line 57 "transregion.dy"
TransFactorRegionPara * new_TransFactorRegionPara_from_argv(int * argc,char ** argv)
{
  TransFactorRegionPara * out;
  char * temp;

  out = TransFactorRegionPara_alloc();

  temp = strip_out_assigned_argument(argc,argv,"tf_convert");
  if( temp != NULL ) {
    if( strcmp(temp,"window") == 0 ) {
      out->type = TRANSREGION_PARA_WINDOW;
    } else if ( strcmp(temp,"hmm") == 0 ) {
      out->type = TRANSREGION_PARA_DP;
    } else {
      fatal("Could not interpret %s for tf_convert string",temp);
    }
  }

  strip_out_float_argument(argc,argv,"tf_density",&out->min_density);
  strip_out_integer_argument(argc,argv,"tf_window",&out->min_window);
  strip_out_integer_argument(argc,argv,"tf_hmm_window",&out->hmm_min_window);
  strip_out_float_argument(argc,argv,"tf_inprob",&out->in_region_prob);
  strip_out_float_argument(argc,argv,"tf_outprob",&out->out_region_prob);
  strip_out_float_argument(argc,argv,"tf_entry",&out->in_cost);
  strip_out_float_argument(argc,argv,"tf_gc_region",&out->gc_region_ratio);

  return out;
}

# line 86 "transregion.dy"
void show_TransFactorRegionSet(TransFactorRegionSet * tfrs,FILE * ofp)
{
  int i;
  int j;
  Sequence * temp;

  for(i=0;i<tfrs->len;i++) {
    fprintf(ofp,"Region\t%s\t%d\t%d\t%.2f\t%.2f\n",tfrs->target->name,tfrs->region[i]->start+1,tfrs->region[i]->end,tfrs->region[i]->density,tfrs->region[i]->bits_score);
    fprintf(ofp,"motif\n");
    for(j=0;j<tfrs->region[i]->len;j++) {
      auto TransFactorMatch * tfm = tfrs->region[i]->match[j];
      fprintf(ofp,"Motif\t%s\t%d\t%d\t%d\t%s\t%.2f\t%.*s\n",tfrs->target->name,tfm->start+1,tfm->end,tfm->strand,tfm->factor->name,tfm->bit_score,tfm->end-tfm->start,tfrs->target->seq+tfm->start);
    }
    fprintf(ofp,"end motif\n");
    temp = trunc_Sequence(tfrs->target,tfrs->region[i]->start,tfrs->region[i]->end);
    write_fasta_Sequence(temp,ofp);
    free_Sequence(temp);
    fprintf(ofp,"end region\n");
  }


}

# line 109 "transregion.dy"
TransFactorRegionSet * new_TransFactorRegionSet(TransFactorMatchSet * tfms,TransFactorRegionPara * tfrp,DPRunImpl * dpri)
{
  switch(tfrp->type) {
  case TRANSREGION_PARA_DP :
    return new_dp_TransFactorRegionSet(tfms,tfrp,dpri);
    break;
  case TRANSREGION_PARA_WINDOW :
    return new_window_TransFactorRegionSet(tfms,tfrp);
    break;
  default :
    fatal("Very weird. Bad type for TransFactorRegionPara");
  }
 
  /* can't get here... but ... */
  return NULL;
}


# line 127 "transregion.dy"
TransFactorRegionSet * new_dp_TransFactorRegionSet(TransFactorMatchSet * tfms,TransFactorRegionPara * tfrp,DPRunImpl * dpri)
{
  AlnBlock * alb;
  PackAln * pal;
  AlnColumn * alc;

  SequenceBaseCoverage * sbc;
  TransRegionModel * model;

  TransFactorRegionSet * out;
  TransFactorRegion * region;

  int i;
  int covered;
  int uncovered;

  assert(tfms);
  assert(tfrp);
  assert(dpri);

  sbc = new_SequenceBaseCoverage(tfms);

  fprintf(stderr,"Making model with %.2f vs %.2f\n",tfrp->in_region_prob,tfrp->out_region_prob);

  model = new_logodds_TransRegionModel(tfrp->in_region_prob,tfrp->out_region_prob,tfrp->in_cost,tfrp->gc_region_ratio);

  fprintf(stderr,"GC score is %d (%.2f)\n",model->gc_point,Score2Bits(model->gc_point));

  pal = PackAln_bestmemory_TransRegionMatrix(model,sbc,NULL,dpri);

  alb = convert_PackAln_to_AlnBlock_TransRegionMatrix(pal);

  free_PackAln(pal);



  out = TransFactorRegionSet_alloc_std();
  out->target = hard_link_Sequence(tfms->target);

  for(alc = alb->start;alc != NULL;alc = alc->next ) {
    if( strstr(alc->alu[1]->text_label,"REGION") != NULL ) {
      if( alc->alu[1]->end - alc->alu[1]->start < tfrp->hmm_min_window ) {
	continue;
      }

      region = TransFactorRegion_alloc_std();
      add_TransFactorRegionSet(out,region);

      region->start = alc->alu[1]->start+1;
      region->end   = alc->alu[1]->end;
      region->bits_score = Score2Bits(alc->alu[0]->score[0]);

      covered = 0;
      uncovered = 0;
      for(i=alc->alu[1]->start+1;i<alc->alu[1]->end;i++) {
	if(sbc->coverage[i] == 0 ) {
	  uncovered++;
	} else {
	  covered++;
	}
      }

      region->density = (double) covered / (double)(covered+uncovered);

      for(i=0;i<tfms->len;i++) {
	if( tfms->match[i]->end < region->start ) {
	  continue;
	}
	if( tfms->match[i]->start > region->end ) {
	  break;
	}

	add_TransFactorRegion(region,hard_link_TransFactorMatch(tfms->match[i]));
      }

    }
  }


  free_SequenceBaseCoverage(sbc);
  free_TransRegionModel(model);
  free_AlnBlock(alb);

  return out;
    
}


# line 215 "transregion.dy"
TransFactorRegionSet * new_window_TransFactorRegionSet(TransFactorMatchSet * tfms,TransFactorRegionPara * tfrp)
{
  int seqpos;
  int motifpos;
  int seq_trial;
  int motif_trial;

  int end;

  int covered_bases;
  int i;
  int temp_start;
  int last_covered_base;

  TransFactorRegionSet * out;
  TransFactorRegion * region;

  assert(tfms);
  assert(tfrp);

  out = TransFactorRegionSet_alloc_std();
  out->target = hard_link_Sequence(tfms->target);

  if( tfms->len == 0 ) {
    /* not motifs... no regions! */
    return out;
  }


  sort_by_start_TransFactorMatchSet(tfms);


  end = tfms->match[tfms->len-1]->start;

  for(seqpos = tfms->match[0]->start,motifpos = 0;seqpos < end && motifpos < tfms->len;) {
    /* see whether there is a potential region here */
    covered_bases =0;
    last_covered_base = seqpos;
    for(seq_trial = seqpos, motif_trial = motifpos; seq_trial < end && seq_trial - seqpos < tfrp->min_window;seq_trial++) {

      for(;motif_trial < tfms->len;motif_trial++ ) {
	if( seq_trial >= tfms->match[motif_trial]->start && seq_trial < tfms->match[motif_trial]->end) {
	  covered_bases++;
	  last_covered_base = seq_trial;
	  break;
	}

	if( seq_trial < tfms->match[motif_trial]->start ) {
	  break;
	}
      }


    }

    if( seq_trial - seqpos < tfrp->min_window ) {
      break;
    }


    /* seq_trial is beyond window size... */

    if( (double) covered_bases / (double) (seq_trial - seqpos) < tfrp->min_density ) {
      motifpos++;
      if( motifpos >= tfms->len ) {
	break;
      } else {
	seqpos = tfms->match[motifpos]->start;
	continue; /* next motif start point */
      }
    }

    seq_trial = last_covered_base;

    /* has ok density now. Extend until we have bad density */

    /* this is not a great extension - each new motif must be providing
       at least min_density bases */
    for(motif_trial;motif_trial < tfms->len;motif_trial++) {
      if( seq_trial >= tfms->match[motif_trial]->end ) {
	continue;
      }

	
      if( (tfms->match[motif_trial]->end - tfms->match[motif_trial]->start) / (double) (tfms->match[motif_trial]->end - seq_trial) < tfrp->min_density ) {
	break;
      } else {
	temp_start = seq_trial > tfms->match[motif_trial]->start ? seq_trial : tfms->match[motif_trial]->start;

	seq_trial = tfms->match[motif_trial]->end;
	covered_bases += (tfms->match[motif_trial]->end - temp_start);
      }
    }
  
    /* we have a match! from seqpos to seq_trial and from motifpos to motiftrial */

    region = TransFactorRegion_alloc_std();
    region->start = seqpos;
    region->end   = seq_trial;
    region->density = (double) covered_bases / (double) (seq_trial - seqpos);
    
    for(i=0;motifpos+i < motif_trial;i++) {
      add_TransFactorRegion(region,hard_link_TransFactorMatch(tfms->match[motifpos+i]));
    }

    add_TransFactorRegionSet(out,region);
    motifpos = motif_trial++;
    if( motifpos >= tfms->len ) {
      break;
    }
    seqpos   = tfms->match[motifpos]->start;

  }

  return out;

}
# line 301 "transregion.c"
/* Function:  swap_TransFactorRegion(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_TransFactorRegion
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [TransFactorMatch **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_TransFactorRegion(TransFactorMatch ** list,int i,int j)  
{
    TransFactorMatch * temp; 
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_TransFactorRegion(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_TransFactorRegion which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [TransFactorMatch **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_TransFactorRegion(TransFactorMatch ** list,int left,int right,int (*comp)(TransFactorMatch * ,TransFactorMatch * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_TransFactorRegion(list,left,(left+right)/2);    
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_TransFactorRegion (list,++last,i);  
      }  
    swap_TransFactorRegion (list,left,last); 
    qsort_TransFactorRegion(list,left,last-1,comp);  
    qsort_TransFactorRegion(list,last+1,right,comp); 
}    


/* Function:  sort_TransFactorRegion(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_TransFactorRegion
 *
 *
 * Arg:         obj [UNKN ] Object containing list [TransFactorRegion *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_TransFactorRegion(TransFactorRegion * obj,int (*comp)(TransFactorMatch *, TransFactorMatch *)) 
{
    qsort_TransFactorRegion(obj->match,0,obj->len-1,comp);   
    return;  
}    


/* Function:  expand_TransFactorRegion(obj,len)
 *
 * Descrip:    Really an internal function for add_TransFactorRegion
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [TransFactorRegion *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_TransFactorRegion(TransFactorRegion * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_TransFactorRegion called with no need");  
      return TRUE;   
      }  


    if( (obj->match = (TransFactorMatch ** ) ckrealloc (obj->match,sizeof(TransFactorMatch *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_TransFactorRegion, returning FALSE");    
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_TransFactorRegion(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [TransFactorRegion *]
 * Arg:        add [OWNER] Object to add to the list [TransFactorMatch *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_TransFactorRegion(TransFactorRegion * obj,TransFactorMatch * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_TransFactorRegion(obj,obj->len + TransFactorRegionLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->match[obj->len++]=add;  
    return TRUE; 
}    


/* Function:  flush_TransFactorRegion(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [TransFactorRegion *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_TransFactorRegion(TransFactorRegion * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->match[i] != NULL) {  
        free_TransFactorMatch(obj->match[i]);    
        obj->match[i] = NULL;    
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  TransFactorRegion_alloc_std(void)
 *
 * Descrip:    Equivalent to TransFactorRegion_alloc_len(TransFactorRegionLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TransFactorRegion *]
 *
 */
TransFactorRegion * TransFactorRegion_alloc_std(void) 
{
    return TransFactorRegion_alloc_len(TransFactorRegionLISTLENGTH); 
}    


/* Function:  TransFactorRegion_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [TransFactorRegion *]
 *
 */
TransFactorRegion * TransFactorRegion_alloc_len(int len) 
{
    TransFactorRegion * out;/* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = TransFactorRegion_alloc()) == NULL)    
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->match = (TransFactorMatch ** ) ckcalloc (len,sizeof(TransFactorMatch *))) == NULL)  {  
      warn("Warning, ckcalloc failed in TransFactorRegion_alloc_len");   
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_TransFactorRegion(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [TransFactorRegion *]
 *
 * Return [UNKN ]  Undocumented return value [TransFactorRegion *]
 *
 */
TransFactorRegion * hard_link_TransFactorRegion(TransFactorRegion * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a TransFactorRegion object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  TransFactorRegion_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TransFactorRegion *]
 *
 */
TransFactorRegion * TransFactorRegion_alloc(void) 
{
    TransFactorRegion * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(TransFactorRegion *) ckalloc (sizeof(TransFactorRegion))) == NULL)  {  
      warn("TransFactorRegion_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->start = 0;  
    out->end = 0;    
    out->match = NULL;   
    out->len = out->maxlen = 0;  
    out->density = 0;    
    out->bits_score = 0; 


    return out;  
}    


/* Function:  free_TransFactorRegion(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [TransFactorRegion *]
 *
 * Return [UNKN ]  Undocumented return value [TransFactorRegion *]
 *
 */
TransFactorRegion * free_TransFactorRegion(TransFactorRegion * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a TransFactorRegion obj. Should be trappable"); 
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
    if( obj->match != NULL)  {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->match[i] != NULL)   
          free_TransFactorMatch(obj->match[i]);  
        }  
      ckfree(obj->match);    
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_TransFactorRegionSet(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_TransFactorRegionSet
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [TransFactorRegion **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_TransFactorRegionSet(TransFactorRegion ** list,int i,int j)  
{
    TransFactorRegion * temp;    
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_TransFactorRegionSet(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_TransFactorRegionSet which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [TransFactorRegion **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_TransFactorRegionSet(TransFactorRegion ** list,int left,int right,int (*comp)(TransFactorRegion * ,TransFactorRegion * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_TransFactorRegionSet(list,left,(left+right)/2); 
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_TransFactorRegionSet (list,++last,i);   
      }  
    swap_TransFactorRegionSet (list,left,last);  
    qsort_TransFactorRegionSet(list,left,last-1,comp);   
    qsort_TransFactorRegionSet(list,last+1,right,comp);  
}    


/* Function:  sort_TransFactorRegionSet(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_TransFactorRegionSet
 *
 *
 * Arg:         obj [UNKN ] Object containing list [TransFactorRegionSet *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_TransFactorRegionSet(TransFactorRegionSet * obj,int (*comp)(TransFactorRegion *, TransFactorRegion *)) 
{
    qsort_TransFactorRegionSet(obj->region,0,obj->len-1,comp);   
    return;  
}    


/* Function:  expand_TransFactorRegionSet(obj,len)
 *
 * Descrip:    Really an internal function for add_TransFactorRegionSet
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [TransFactorRegionSet *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_TransFactorRegionSet(TransFactorRegionSet * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_TransFactorRegionSet called with no need");   
      return TRUE;   
      }  


    if( (obj->region = (TransFactorRegion ** ) ckrealloc (obj->region,sizeof(TransFactorRegion *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_TransFactorRegionSet, returning FALSE"); 
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_TransFactorRegionSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [TransFactorRegionSet *]
 * Arg:        add [OWNER] Object to add to the list [TransFactorRegion *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_TransFactorRegionSet(TransFactorRegionSet * obj,TransFactorRegion * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_TransFactorRegionSet(obj,obj->len + TransFactorRegionSetLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->region[obj->len++]=add; 
    return TRUE; 
}    


/* Function:  flush_TransFactorRegionSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [TransFactorRegionSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_TransFactorRegionSet(TransFactorRegionSet * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->region[i] != NULL)    {  
        free_TransFactorRegion(obj->region[i]);  
        obj->region[i] = NULL;   
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  TransFactorRegionSet_alloc_std(void)
 *
 * Descrip:    Equivalent to TransFactorRegionSet_alloc_len(TransFactorRegionSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TransFactorRegionSet *]
 *
 */
TransFactorRegionSet * TransFactorRegionSet_alloc_std(void) 
{
    return TransFactorRegionSet_alloc_len(TransFactorRegionSetLISTLENGTH);   
}    


/* Function:  TransFactorRegionSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [TransFactorRegionSet *]
 *
 */
TransFactorRegionSet * TransFactorRegionSet_alloc_len(int len) 
{
    TransFactorRegionSet * out; /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = TransFactorRegionSet_alloc()) == NULL) 
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->region = (TransFactorRegion ** ) ckcalloc (len,sizeof(TransFactorRegion *))) == NULL)   {  
      warn("Warning, ckcalloc failed in TransFactorRegionSet_alloc_len");    
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_TransFactorRegionSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [TransFactorRegionSet *]
 *
 * Return [UNKN ]  Undocumented return value [TransFactorRegionSet *]
 *
 */
TransFactorRegionSet * hard_link_TransFactorRegionSet(TransFactorRegionSet * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a TransFactorRegionSet object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  TransFactorRegionSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TransFactorRegionSet *]
 *
 */
TransFactorRegionSet * TransFactorRegionSet_alloc(void) 
{
    TransFactorRegionSet * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(TransFactorRegionSet *) ckalloc (sizeof(TransFactorRegionSet))) == NULL)    {  
      warn("TransFactorRegionSet_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->region = NULL;  
    out->len = out->maxlen = 0;  
    out->target = NULL;  


    return out;  
}    


/* Function:  free_TransFactorRegionSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [TransFactorRegionSet *]
 *
 * Return [UNKN ]  Undocumented return value [TransFactorRegionSet *]
 *
 */
TransFactorRegionSet * free_TransFactorRegionSet(TransFactorRegionSet * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a TransFactorRegionSet obj. Should be trappable");  
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
    if( obj->region != NULL) {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->region[i] != NULL)  
          free_TransFactorRegion(obj->region[i]);    
        }  
      ckfree(obj->region);   
      }  
    if( obj->target != NULL) 
      free_Sequence(obj->target);    


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_TransFactorRegionPara(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [TransFactorRegionPara *]
 *
 * Return [UNKN ]  Undocumented return value [TransFactorRegionPara *]
 *
 */
TransFactorRegionPara * hard_link_TransFactorRegionPara(TransFactorRegionPara * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a TransFactorRegionPara object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  TransFactorRegionPara_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TransFactorRegionPara *]
 *
 */
TransFactorRegionPara * TransFactorRegionPara_alloc(void) 
{
    TransFactorRegionPara * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(TransFactorRegionPara *) ckalloc (sizeof(TransFactorRegionPara))) == NULL)  {  
      warn("TransFactorRegionPara_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->type = TRANSREGION_PARA_DP; 
    out->min_density = 0.3;  
    out->min_window = 60;    
    out->in_region_prob = 0.35;  
    out->out_region_prob = 0.2;  
    out->in_cost = 0.000001; 
    out->hmm_min_window = 60;    
    out->gc_region_ratio = 1.5;  


    return out;  
}    


/* Function:  free_TransFactorRegionPara(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [TransFactorRegionPara *]
 *
 * Return [UNKN ]  Undocumented return value [TransFactorRegionPara *]
 *
 */
TransFactorRegionPara * free_TransFactorRegionPara(TransFactorRegionPara * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a TransFactorRegionPara obj. Should be trappable"); 
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



#ifdef _cplusplus
}
#endif
