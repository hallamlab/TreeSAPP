#ifdef _cplusplus
extern "C" {
#endif
#include "localcishit.h"


# line 49 "localcishit.dy"
HitList * HitList_from_LocalCisHitSet(LocalCisHitSet * in)
{
  HitList * out;
  HitPair * p;
  HitAln * a;
  int i;

  out = HitList_alloc_std();

  for(i=0;i<in->len;i++) {
    p = HitPair_alloc_std();
    add_HitList(out,p);
    p->query = hard_link_Sequence(in->lch[i]->query);
    p->target = hard_link_Sequence(in->lch[i]->target);

    p->raw_score = in->lch[i]->score;
    p->bit_score = Score2Bits(p->raw_score);
    p->target_reversed = in->lch[i]->target_rev;
    a = HitAln_alloc();
    add_HitPair(p,a);
    a->alb = hard_link_AlnBlock(in->lch[i]->alb);
    a->raw_score = p->raw_score;
    a->bit_score = p->bit_score;

  }

  out->write_btc_func = show_pretty_Seq_dba_align_btcanvas;

  return out;
}

# line 80 "localcishit.dy"
LocalCisHitSet * expand_to_subhits_LocalCisHitSet(LocalCisHitSet * in)
{
  int i;
  AlnBlock * new;
  AlnColumn * alc;
  LocalCisHitSet * out;
  int qstart;
  int tstart;
  int qend;
  int tend;
  int score_start;
  int score;

  out = LocalCisHitSet_alloc_std();

  for(i=0;in->len;i++) {
    for(alc=in->lch[i]->alb->start;alc != NULL;alc = alc->next ) {
      if( strcmp(alc->alu[0]->text_label,"UNMATCHED_QUERY") == 0 ) {
	break;
      }
    }
    if( alc == NULL ) {
      /* just one hit */
      add_LocalCisHitSet(out,hard_link_LocalCisHit(in->lch[i]));
    } else {

      for(alc=in->lch[i]->alb->start;alc != NULL;alc = alc->next) {
	if( strstr(alc->alu[0]->text_label,"MM") != NULL ) {
	  /* start here */
	  qstart = alc->alu[0]->start;
	  tstart = alc->alu[1]->start;
	  score = 0;
	  new = AlnBlock_alloc();
	  new->start = alc; 
	  for(;alc != NULL && strstr(alc->alu[0]->text_label,"UNMATCHED") == NULL; alc = alc->next) {
	  }
	}
      }

      fatal("Ewan has not implemented subhit expansion!");
    }
  }


  return out;
}

# line 127 "localcishit.dy"
void show_help_LocalCisHitSetPara(FILE * ofp)
{
  fprintf(ofp,"Local Hit expansion parameters\n");
  fprintf(ofp,"  -lhwindow    - sequence window given to alignment [50]\n");
  fprintf(ofp,"  -lhseed      - seed score cutoff [10.0 bits]\n");
  fprintf(ofp,"  -lhaln       - aln  score cutoff [8.0 bits]\n");
  fprintf(ofp,"  -lhscore     - sort final list by score (default by position)\n");
  fprintf(ofp,"  -lhreject [none/query/both] - overlap rejection criteria in greedy assembly [query]\n");
  fprintf(ofp,"  -lhmax    [20000] - maximum number of processed hits\n");

}

# line 139 "localcishit.dy"
LocalCisHitSetPara * new_LocalCisHitSetPara_from_argv(int * argc,char ** argv)
{
  LocalCisHitSetPara * setpara;
  char * temp;

  setpara = LocalCisHitSetPara_alloc();
  setpara->expansion_size = 50;
  setpara->seed_bit_trigger = 10.0;
  setpara->aln_cutoff = 8.0;
  setpara->sort_by_score = 0;
  setpara->max = 20000;
  setpara->type = LocalCisGreedy_Query;

  strip_out_integer_argument(argc,argv,"lhwindow",&setpara->expansion_size);
  strip_out_integer_argument(argc,argv,"lhmax",&setpara->max);
  strip_out_float_argument(argc,argv,"lhseed",&setpara->seed_bit_trigger);
  strip_out_float_argument(argc,argv,"lhaln",&setpara->aln_cutoff);
  strip_out_boolean_def_argument(argc,argv,"lhscore",&setpara->sort_by_score);

  temp = strip_out_assigned_argument(argc,argv,"lhreject");
  if( temp != NULL ) {
    if( strcmp(temp,"none") == 0) {
      setpara->type = LocalCisGreedy_None;
    } else if ( strcmp(temp,"query") == 0 ) {
      setpara->type = LocalCisGreedy_Query;
    } else if ( strcmp(temp,"both") == 0){
      setpara->type = LocalCisGreedy_Both;
    } else {
      fatal("Bad parameter for lhreject %s",temp);
    }
  }


  return setpara;
}

# line 175 "localcishit.dy"
void show_pretty_LocalCisHitSet(LocalCisHitSet * lchs,FILE * ofp)
{
  int i;


  for(i=0;i<lchs->len;i++) {
    fprintf(ofp,">%s %5d,%5d  %s %5d,%5d [%c]  Bits %.2f\n",
	    lchs->lch[i]->query->name,
	    lchs->lch[i]->start_q+1,
	    lchs->lch[i]->end_q,
	    lchs->lch[i]->target->name,
	    lchs->lch[i]->start_t+1,
	    lchs->lch[i]->end_t,
	    lchs->lch[i]->target_rev == 1 ? '-' : '+',	    
	    Score2Bits(lchs->lch[i]->score)
	    );
    fprintf(ofp,"\n");
    show_pretty_dba_align(lchs->lch[i]->alb,lchs->lch[i]->query,lchs->lch[i]->target,ofp);
  }


}


# line 199 "localcishit.dy"
void show_summary_LocalCisHitSet(LocalCisHitSet * lchs,FILE * ofp)
{
  int i;

  assert(lchs);
  assert(ofp);


  for(i=0;i<lchs->len;i++) {
    fprintf(ofp,"Query %5d,%5d  Target %5d,%5d [%c]  Bits %.2f\n",
	    lchs->lch[i]->start_q+1,
	    lchs->lch[i]->end_q,
	    lchs->lch[i]->start_t+1,
	    lchs->lch[i]->end_t,
	    lchs->lch[i]->target_rev == 1 ? '-' : '+',	    
	    Score2Bits(lchs->lch[i]->score)
	    );
  }


}

# line 221 "localcishit.dy"
LocalCisHitSet * greedy_weed_LocalCisHitSet(LocalCisHitSet * set,LocalCisHitSetPara *p)
{
  LocalCisHitSet * out;
  int i;
  int j;
  int is_valid;

  sort_LocalCisHitSet_by_score(set);

  out = LocalCisHitSet_alloc_std();

  for(i=0;i<set->len;i++) {

    if( Score2Bits(set->lch[i]->score) < p->aln_cutoff) {
      break;
    }

    is_valid = 1;

    if( p->type != LocalCisGreedy_None ) {
      for(j=0;j<out->len;j++) {
	if( is_query_overlap_LocalCisHit(set->lch[i],out->lch[j]) == 1 ) {
	  is_valid = 0;
	  break;
	}
      }
      if( is_valid == 1 && p->type == LocalCisGreedy_Both ) {
	for(j=0;j<out->len;j++) {
	  if( is_target_overlap_LocalCisHit(set->lch[i],out->lch[j]) == 1 ) {
	    is_valid = 0;
	    break;
	  }
	}
      }
    }

    if( is_valid == 1 ) {
      add_LocalCisHitSet(out,hard_link_LocalCisHit(set->lch[i]));
    }

  }

  if( p->sort_by_score == 0 ) {
    sort_LocalCisHitSet_by_start(out);
  }

  return out;
}

# line 270 "localcishit.dy"
void sort_LocalCisHitSet_by_score(LocalCisHitSet * set)
{
  sort_LocalCisHitSet(set,compare_LocalCisHit_score);
}

# line 275 "localcishit.dy"
void sort_LocalCisHitSet_by_start(LocalCisHitSet * set)
{
  sort_LocalCisHitSet(set,compare_LocalCisHit_start);
}

# line 280 "localcishit.dy"
int compare_LocalCisHit_score(LocalCisHit * a,LocalCisHit * b)
{
  return b->score - a->score;
}

# line 285 "localcishit.dy"
int compare_LocalCisHit_start(LocalCisHit * a,LocalCisHit * b)
{
  return a->start_q - b->start_q;
}

# line 290 "localcishit.dy"
int is_query_overlap_LocalCisHit(LocalCisHit * a,LocalCisHit * b)
{
  if( a->start_q > b->end_q || a->end_q < b->start_q ) {
    return 0;
  } else {       
    return 1;
  }
}

# line 299 "localcishit.dy"
int is_target_overlap_LocalCisHit(LocalCisHit * a,LocalCisHit * b)
{
  if( a->start_t > b->end_t || a->end_t < b->start_t ) {
    return 0;
  } else {       
    return 1;
  }
}


# line 309 "localcishit.dy"
int is_consumed_HSP(HSP * a,int q_start,int q_end,int t_start,int t_end)
{
  int qcentre;
  int tcentre;
  /* find centre of hit, and see if it is in the system */
  
  qcentre = a->query_start + (a->length/2);
  tcentre = a->target_start + (a->length/2);

  if( qcentre >= q_start && qcentre <= q_end &&
      tcentre >= t_start && tcentre <= t_end ) {
    return 1;
  } else {
    return 0;
  }
}

# line 326 "localcishit.dy"
LocalCisHitSet * make_LocalCisHitSet(Sequence * query,Sequence * target,Sequence * target_rev,HSPset * forward,HSPset * reverse,LocalCisHitSetPara * p,LocalCisHitScore * scorepara,TransFactorMatchSet * tfms_query,TransFactorMatchSet * tfms_target,TransFactorMatchSet * tfms_target_rev,MotifMatrixScore * mms,boolean use_motif,DPRunImpl * dpri)
{
  LocalCisHitSet * set;
  LocalCisHit * lch;
  int i,j;
  int * forward_used;
  int * reverse_used;

  int q_start;
  int q_end;
  int t_start;
  int t_end;

  int temp_q_start;
  int temp_q_end;
  int temp_t_start;
  int temp_t_end;

  assert(query);
  assert(target);
  assert(target_rev);
  assert(forward);
  assert(reverse);

  sort_HSPset_by_score(forward);
  sort_HSPset_by_score(reverse);

  set = LocalCisHitSet_alloc_std();
  forward_used = calloc(forward->len,sizeof(int));
  for(i=0;i<forward->len;i++) {
    forward_used[i] = 0;
  }


  for(i=0;i<forward->len && i < p->max;i++) {
    /*    fprintf(stderr,"Looking at hsp set %d %d %f score vs trigger %f\n",i,forward->hsp[i]->score,Score2Bits(forward->hsp[i]->score),p->seed_bit_trigger); */
    if( forward_used[i] == 0 && Score2Bits(forward->hsp[i]->score) > p->seed_bit_trigger ) {
      forward_used[i] = 1;

      q_start = forward->hsp[i]->query_start-p->expansion_size;
      q_end   = forward->hsp[i]->query_start+forward->hsp[i]->length+p->expansion_size;
      t_start = forward->hsp[i]->target_start-p->expansion_size;
      t_end   = forward->hsp[i]->target_start+forward->hsp[i]->length+p->expansion_size;
	

      for(j=i;j<forward->len && j < p->max;j++) {
	if( forward_used[j] == 1 ) {
	  continue;
	}
	if( is_consumed_HSP(forward->hsp[j],q_start,q_end,t_start,t_end) ) {
	  forward_used[j] = 1;
	  /*info("consuming forward %d by %d\n",j,i);*/

	  temp_q_start = forward->hsp[j]->query_start-p->expansion_size;
	  temp_q_end   = forward->hsp[j]->query_start+forward->hsp[i]->length+p->expansion_size;
	  temp_t_start = forward->hsp[j]->target_start-p->expansion_size;
	  temp_t_end   = forward->hsp[j]->target_start+forward->hsp[i]->length+p->expansion_size;

	  q_start = temp_q_start < q_start ? temp_q_start : q_start;
	  t_start = temp_t_start < t_start ? temp_t_start : t_start;

	  t_end = temp_t_end > t_end ? temp_t_end : t_end;
	  q_end = temp_q_end > q_end ? temp_q_end : q_end;

	}
      }

      /* now build local cis hit */


      if( use_motif == 0 ) {
	lch = make_LocalCisHit(query,target,0,
			       q_start,
			       q_end,
			       t_start,
			       t_end,
			       scorepara,dpri
			       );
      } else {
	lch = make_motif_LocalCisHit(query,target,0,
				     q_start,
				     q_end,
				     t_start,
				     t_end,
				     tfms_query,
				     tfms_target,
				     mms,dpri
				     );

      }
      add_LocalCisHitSet(set,lch);
    } else {
      break;
    }
  }


  reverse_used = calloc(reverse->len,sizeof(int));
  for(i=0;i<reverse->len;i++) {
    reverse_used[i] = 0;
  }


  for(i=0;i<reverse->len && i < p->max ;i++) {
    if( reverse_used[i] == 0 && Score2Bits(reverse->hsp[i]->score) > p->seed_bit_trigger ) {

      reverse_used[i] = 1;

      q_start = reverse->hsp[i]->query_start-p->expansion_size;
      q_end   = reverse->hsp[i]->query_start+reverse->hsp[i]->length+p->expansion_size;
      t_start = reverse->hsp[i]->target_start-p->expansion_size;
      t_end   = reverse->hsp[i]->target_start+reverse->hsp[i]->length+p->expansion_size;
	
      for(j=i;j<reverse->len && j < p->max ;j++) {
	if( reverse_used[j] == 1 ) {
	  continue;
	}
	if( is_consumed_HSP(reverse->hsp[j],q_start,q_end,t_start,t_end) ) {
	  reverse_used[j] = 1;

	  /*info("consuming reverse %d by %d\n",j,i);*/

	  temp_q_start = reverse->hsp[j]->query_start-p->expansion_size;
	  temp_q_end   = reverse->hsp[j]->query_start+reverse->hsp[i]->length+p->expansion_size;
	  temp_t_start = reverse->hsp[j]->target_start-p->expansion_size;
	  temp_t_end   = reverse->hsp[j]->target_start+reverse->hsp[i]->length+p->expansion_size;

	  q_start = temp_q_start < q_start ? temp_q_start : q_start;
	  t_start = temp_t_start < t_start ? temp_t_start : t_start;

	  t_end = temp_t_end > t_end ? temp_t_end : t_end;
	  q_end = temp_q_end > q_end ? temp_q_end : q_end;

	}
      }


      if( use_motif == 0 ) {
	lch = make_LocalCisHit(query,target_rev,1,
			       q_start,
			       q_end,
			       t_start,
			       t_end,
			       scorepara,dpri
			       );
      } else {
	lch = make_motif_LocalCisHit(query,target_rev,1,
				     q_start,
				     q_end,
				     t_start,
				     t_end,
				     tfms_query,
				     tfms_target_rev,
				     mms,dpri
				     );
      }

      add_LocalCisHitSet(set,lch);
    } else {
      break;
    }
  }

  return set;

}



# line 495 "localcishit.dy"
LocalCisHit * make_LocalCisHit(Sequence * query,Sequence * target,int is_reversed,int guess_q_start,int guess_q_end,int guess_t_start,int guess_t_end,LocalCisHitScore * score,DPRunImpl * dpri)
{
  LocalCisHit * out;	
  Sequence * temp_q;
  Sequence * temp_t;

  ComplexSequence * q_cseq;
  ComplexSequence * t_cseq;
  ComplexSequenceEvalSet * cses;

  PackAln * pal;
  AlnBlock * alb;

  int i;

  assert(query);
  assert(target);
  assert(score);

  if( guess_q_start < 0 ) {
    guess_q_start = 0;
  }
  if( guess_q_end >= query->len ) {
    guess_q_end = query->len-1;
  }

  if( guess_t_start < 0 ) {
    guess_t_start = 0;
  }
  if( guess_t_end >= target->len ) {
    guess_t_end = target->len-1;
  }

				   

  
  temp_q = trunc_Sequence(query,guess_q_start,guess_q_end);
  temp_t = trunc_Sequence(target,guess_t_start,guess_t_end);

  cses = default_DNA_ComplexSequenceEvalSet();  
  q_cseq = new_ComplexSequence(temp_q,cses);  
  t_cseq = new_ComplexSequence(temp_t,cses);  

  assert(q_cseq);
  assert(t_cseq);

  pal = PackAln_bestmemory_LocalDnaMatchBlock(q_cseq,t_cseq,score,NULL,dpri);

  assert(pal);


  for(i=0;i<pal->len;i++) {
    pal->pau[i]->i = pal->pau[i]->i + guess_q_start;
    pal->pau[i]->j = pal->pau[i]->j + guess_t_start;
  }

  alb = convert_PackAln_to_AlnBlock_LocalDnaMatchBlock(pal);

  out = LocalCisHit_alloc();

  out->start_q = pal->pau[0]->i;
  out->end_q   = pal->pau[pal->len-2]->i;

  out->start_t = pal->pau[0]->j;
  out->end_t   = pal->pau[pal->len-1]->j;
  out->score   = pal->score;
  out->target_rev = is_reversed;
  out->alb = alb;
  out->query = hard_link_Sequence(query);
  out->target = hard_link_Sequence(target);

  free_PackAln(pal);

  return out;
}

# line 571 "localcishit.dy"
LocalCisHit * make_motif_LocalCisHit(Sequence * query,Sequence * target,int is_reversed,int guess_q_start,int guess_q_end,int guess_t_start,int guess_t_end,TransFactorMatchSet * tfms_query,TransFactorMatchSet * tfms_target,MotifMatrixScore * mms,DPRunImpl * dpri)
{
  LocalCisHit * out;	
  Sequence * temp_q;
  Sequence * temp_t;

  ComplexSequence * q_cseq;
  ComplexSequence * t_cseq;
  ComplexSequenceEvalSet * cses;

  MotifConsMatrix * mcm;

  PackAln * pal;
  AlnBlock * alb;

  int i;

  assert(query != NULL);
  assert(target != NULL);
  assert(mms != NULL);
  assert(tfms_query != NULL);
  assert(tfms_target != NULL);


  if( guess_q_start < 0 ) {
    guess_q_start = 0;
  }
  if( guess_q_end >= query->len ) {
    guess_q_end = query->len-1;
  }

  if( guess_t_start < 0 ) {
    guess_t_start = 0;
  }
  if( guess_t_end >= target->len ) {
    guess_t_end = target->len-1;
  }

				   

  
  temp_q = trunc_Sequence(query,guess_q_start,guess_q_end);
  temp_t = trunc_Sequence(target,guess_t_start,guess_t_end);

  cses = default_DNA_ComplexSequenceEvalSet();  
  q_cseq = new_ComplexSequence(temp_q,cses);  
  t_cseq = new_ComplexSequence(temp_t,cses);  

  assert(q_cseq);
  assert(t_cseq);

  mcm = new_MotifConsMatrix(tfms_query,guess_q_start,guess_q_end,tfms_target,guess_t_start,guess_t_end);


  pal = PackAln_bestmemory_LocalMotifMatrix(q_cseq,t_cseq,mcm,mms,NULL,dpri);

  assert(pal);


  for(i=0;i<pal->len;i++) {
    pal->pau[i]->i = pal->pau[i]->i + guess_q_start;
    pal->pau[i]->j = pal->pau[i]->j + guess_t_start;
  }

  alb = convert_PackAln_to_AlnBlock_LocalMotifMatrix(pal);

  out = LocalCisHit_alloc();

  out->start_q = pal->pau[0]->i;
  out->end_q   = pal->pau[pal->len-2]->i;

  out->start_t = pal->pau[0]->j;
  out->end_t   = pal->pau[pal->len-1]->j;
  out->score   = pal->score;
  out->target_rev = is_reversed;
  out->alb = alb;
  out->query = hard_link_Sequence(query);
  out->target = hard_link_Sequence(target);

  free_PackAln(pal);

  return out;
}





# line 633 "localcishit.c"
/* Function:  hard_link_LocalCisHit(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [LocalCisHit *]
 *
 * Return [UNKN ]  Undocumented return value [LocalCisHit *]
 *
 */
LocalCisHit * hard_link_LocalCisHit(LocalCisHit * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a LocalCisHit object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  LocalCisHit_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LocalCisHit *]
 *
 */
LocalCisHit * LocalCisHit_alloc(void) 
{
    LocalCisHit * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(LocalCisHit *) ckalloc (sizeof(LocalCisHit))) == NULL)  {  
      warn("LocalCisHit_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->start_q = 0;    
    out->end_q = 0;  
    out->start_t = 0;    
    out->end_t = 0;  
    out->target_rev = 0; 
    out->alb = NULL; 
    out->score = 0;  
    out->query = NULL;   
    out->target = NULL;  


    return out;  
}    


/* Function:  free_LocalCisHit(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [LocalCisHit *]
 *
 * Return [UNKN ]  Undocumented return value [LocalCisHit *]
 *
 */
LocalCisHit * free_LocalCisHit(LocalCisHit * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a LocalCisHit obj. Should be trappable");   
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
    if( obj->alb != NULL)    
      free_AlnBlock(obj->alb);   
    if( obj->query != NULL)  
      free_Sequence(obj->query);     
    if( obj->target != NULL) 
      free_Sequence(obj->target);    


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_LocalCisHitSet(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_LocalCisHitSet
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [LocalCisHit **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_LocalCisHitSet(LocalCisHit ** list,int i,int j)  
{
    LocalCisHit * temp;  
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_LocalCisHitSet(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_LocalCisHitSet which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [LocalCisHit **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_LocalCisHitSet(LocalCisHit ** list,int left,int right,int (*comp)(LocalCisHit * ,LocalCisHit * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_LocalCisHitSet(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_LocalCisHitSet (list,++last,i); 
      }  
    swap_LocalCisHitSet (list,left,last);    
    qsort_LocalCisHitSet(list,left,last-1,comp); 
    qsort_LocalCisHitSet(list,last+1,right,comp);    
}    


/* Function:  sort_LocalCisHitSet(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_LocalCisHitSet
 *
 *
 * Arg:         obj [UNKN ] Object containing list [LocalCisHitSet *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_LocalCisHitSet(LocalCisHitSet * obj,int (*comp)(LocalCisHit *, LocalCisHit *)) 
{
    qsort_LocalCisHitSet(obj->lch,0,obj->len-1,comp);    
    return;  
}    


/* Function:  expand_LocalCisHitSet(obj,len)
 *
 * Descrip:    Really an internal function for add_LocalCisHitSet
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [LocalCisHitSet *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_LocalCisHitSet(LocalCisHitSet * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_LocalCisHitSet called with no need"); 
      return TRUE;   
      }  


    if( (obj->lch = (LocalCisHit ** ) ckrealloc (obj->lch,sizeof(LocalCisHit *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_LocalCisHitSet, returning FALSE");   
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_LocalCisHitSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [LocalCisHitSet *]
 * Arg:        add [OWNER] Object to add to the list [LocalCisHit *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_LocalCisHitSet(LocalCisHitSet * obj,LocalCisHit * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_LocalCisHitSet(obj,obj->len + LocalCisHitSetLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->lch[obj->len++]=add;    
    return TRUE; 
}    


/* Function:  flush_LocalCisHitSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [LocalCisHitSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_LocalCisHitSet(LocalCisHitSet * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->lch[i] != NULL)   {  
        free_LocalCisHit(obj->lch[i]);   
        obj->lch[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  LocalCisHitSet_alloc_std(void)
 *
 * Descrip:    Equivalent to LocalCisHitSet_alloc_len(LocalCisHitSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LocalCisHitSet *]
 *
 */
LocalCisHitSet * LocalCisHitSet_alloc_std(void) 
{
    return LocalCisHitSet_alloc_len(LocalCisHitSetLISTLENGTH);   
}    


/* Function:  LocalCisHitSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [LocalCisHitSet *]
 *
 */
LocalCisHitSet * LocalCisHitSet_alloc_len(int len) 
{
    LocalCisHitSet * out;   /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = LocalCisHitSet_alloc()) == NULL)   
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->lch = (LocalCisHit ** ) ckcalloc (len,sizeof(LocalCisHit *))) == NULL)  {  
      warn("Warning, ckcalloc failed in LocalCisHitSet_alloc_len");  
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_LocalCisHitSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [LocalCisHitSet *]
 *
 * Return [UNKN ]  Undocumented return value [LocalCisHitSet *]
 *
 */
LocalCisHitSet * hard_link_LocalCisHitSet(LocalCisHitSet * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a LocalCisHitSet object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  LocalCisHitSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LocalCisHitSet *]
 *
 */
LocalCisHitSet * LocalCisHitSet_alloc(void) 
{
    LocalCisHitSet * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(LocalCisHitSet *) ckalloc (sizeof(LocalCisHitSet))) == NULL)    {  
      warn("LocalCisHitSet_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->lch = NULL; 
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_LocalCisHitSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [LocalCisHitSet *]
 *
 * Return [UNKN ]  Undocumented return value [LocalCisHitSet *]
 *
 */
LocalCisHitSet * free_LocalCisHitSet(LocalCisHitSet * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a LocalCisHitSet obj. Should be trappable");    
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
    if( obj->lch != NULL)    {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->lch[i] != NULL) 
          free_LocalCisHit(obj->lch[i]); 
        }  
      ckfree(obj->lch);  
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_LocalCisHitSetPara(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [LocalCisHitSetPara *]
 *
 * Return [UNKN ]  Undocumented return value [LocalCisHitSetPara *]
 *
 */
LocalCisHitSetPara * hard_link_LocalCisHitSetPara(LocalCisHitSetPara * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a LocalCisHitSetPara object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  LocalCisHitSetPara_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LocalCisHitSetPara *]
 *
 */
LocalCisHitSetPara * LocalCisHitSetPara_alloc(void) 
{
    LocalCisHitSetPara * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(LocalCisHitSetPara *) ckalloc (sizeof(LocalCisHitSetPara))) == NULL)    {  
      warn("LocalCisHitSetPara_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->seed_bit_trigger = 0;   
    out->expansion_size = 0; 
    out->aln_cutoff = 0; 
    out->sort_by_score = FALSE;  
    out->max = 0;    
    out->type = LocalCisGreedy_Query;    


    return out;  
}    


/* Function:  free_LocalCisHitSetPara(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [LocalCisHitSetPara *]
 *
 * Return [UNKN ]  Undocumented return value [LocalCisHitSetPara *]
 *
 */
LocalCisHitSetPara * free_LocalCisHitSetPara(LocalCisHitSetPara * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a LocalCisHitSetPara obj. Should be trappable");    
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
