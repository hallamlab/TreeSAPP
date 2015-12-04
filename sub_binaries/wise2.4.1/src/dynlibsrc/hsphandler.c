#ifdef _cplusplus
extern "C" {
#endif
#include "hsphandler.h"


/* Function:  new_TopScoreManager(length)
 *
 * Descrip:    Makes a new topscore manager of a specific length
 *
 *
 * Arg:        length [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [TopScoreManager *]
 *
 */
# line 57 "hsphandler.dy"
TopScoreManager * new_TopScoreManager(int length)
{


  TopScoreManager * tsm;

  tsm = TopScoreManager_alloc();
  tsm->score = calloc(length,sizeof(int));
  tsm->current_pos = 0;
  tsm->length = length;

  tsm->worst_score = 0;
  tsm->worst_position = -1;

  return tsm;
}


/* Function:  add_score_TopScoreManager(tsm,score)
 *
 * Descrip:    Adds a top score
 *
 *
 * Arg:          tsm [UNKN ] Undocumented argument [TopScoreManager *]
 * Arg:        score [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 78 "hsphandler.dy"
boolean add_score_TopScoreManager(TopScoreManager * tsm,int score)
{
  int i;
  int temp_score;
  int temp_pos;

/*  fprintf(stderr,"Looking at %d\n",score);*/

  if( tsm->current_pos < tsm->length ) {
    tsm->score[tsm->current_pos++] = score;
    return TRUE;
  }


  if( tsm->worst_position == -1 ) {
    /* fprintf(stderr,"Recalculating...\n"); */
    /* need to recalculate top score */
    for(i=1,temp_score = tsm->score[0],temp_pos = 0;i<tsm->length;i++) {
      if( temp_score > tsm->score[i] ) {
	temp_score = tsm->score[i];
	temp_pos = i;
      }
    }
    tsm->worst_score = temp_score;
    tsm->worst_position = temp_pos;
  }

  if( score < tsm->worst_score ) {
/*    fprintf(stderr,"Ignoring... %d vs %d\n",score,tsm->worst_score); */
    return FALSE;
  }

/*  fprintf(stderr,"Resetting...\n"); */
  tsm->score[tsm->worst_position] = score;
  tsm->worst_position = -1;

  return TRUE;
}

/* Function:  truncated_simple_LinearHSPmanager(lm,para)
 *
 * Descrip:    A simpler truncation method, using diagonals
 *
 *
 * Arg:          lm [UNKN ] Undocumented argument [LinearHSPmanager *]
 * Arg:        para [UNKN ] Undocumented argument [LineariseHSPPara *]
 *
 * Return [UNKN ]  Undocumented return value [LinearHSPmanager *]
 *
 */
# line 120 "hsphandler.dy"
LinearHSPmanager * truncated_simple_LinearHSPmanager(LinearHSPmanager * lm,LineariseHSPPara * para)
{
  LinearHSPmanager * out;
  int i;
  int j;
  int new_score;

  int diag_one;
  int diag_two;

  int query_position[500];
  int k;
  int aa_count;
  int query_offset;
/*
  struct timeval t1;
  struct timeval t2;
  struct timeval t3;
*/

/*  gettimeofday(&t1,NULL);*/

  if( VERBOSITY_CHECK(2,para->verbosity) ) {
    info("Input linear management of %d entries\n",lm->len);
  }


  for(i=0;i<lm->len;i++) {
    /* rescore this set wrt to top diagonal */
    sort_HSPset_by_score(lm->set[i]);
    new_score = lm->set[i]->hsp[0]->score;
    diag_one = lm->set[i]->hsp[0]->query_start - lm->set[i]->hsp[0]->target_start;

    if( VERBOSITY_CHECK(5,para->verbosity) ) {
      info("Looking at %s with starting hsp score of %d\n",lm->set[i]->hsp[0]->target->name,new_score);
    }
    
    query_offset = lm->set[i]->hsp[0]->query_start - 200;

    for(k=0;k<500;k++) {
      query_position[k] = 0;
    }
    for(k=0;k<lm->set[i]->hsp[0]->length && lm->set[i]->hsp[0]->query_start - query_offset+k < 500;k++) {
      query_position[lm->set[i]->hsp[0]->query_start - query_offset+k] = 1;
    }
    

    for(j=1;j<lm->set[i]->len;j++) {

      if( lm->set[i]->hsp[j]->score < para->min_score ) {
	if( VERBOSITY_CHECK(8,para->verbosity) ) {
	  info("  ...not accepting hsp %d,%d due to min score",lm->set[i]->hsp[j]->query_start,lm->set[i]->hsp[j]->target_start);
	}

	break;
      }

      diag_two = lm->set[i]->hsp[j]->query_start - lm->set[i]->hsp[j]->target_start;
      if( abs(diag_one - diag_two) < para->width && abs(lm->set[i]->hsp[0]->query_start - lm->set[i]->hsp[j]->query_start) < para->tail) {

	/* this is now in the right area                     */
	/* but we need to test for overlap on existing cases */


	/*	fprintf(stderr,"  .... looking at hsp overlap\n");*/

	for(aa_count = 0,k=0;k<lm->set[i]->hsp[j]->length && lm->set[i]->hsp[j]->query_start - query_offset+k < 500;k++) {
	  if( query_position[lm->set[i]->hsp[j]->query_start - query_offset+k] == 1 ) {
	    aa_count++;
	  }
	}


	/* if we have more than 15% of this hsp "accounted for" move on */

	if(aa_count/(double)lm->set[i]->hsp[j]->length > 0.15 ) {
	  continue;
	}
	
	/* fprintf(stderr,"   .... updating HSP overlap\n");*/

	/* now set these positions as used */

	for(aa_count = 0,k=0;k<lm->set[i]->hsp[j]->length && lm->set[i]->hsp[j]->query_start - query_offset+k < 500;k++) {
	  query_position[lm->set[i]->hsp[j]->query_start - query_offset+k] = 1;
	}

	new_score += lm->set[i]->hsp[j]->score;
	if( VERBOSITY_CHECK(5,para->verbosity) ) {
	  info("  ..accepting hsp on %s %d,%d new score %d",lm->set[i]->hsp[j]->target->name,lm->set[i]->hsp[j]->query_start,lm->set[i]->hsp[j]->target_start,new_score); 
	}
      } else {
	if( VERBOSITY_CHECK(5,para->verbosity) ) {
	  info("  ...not accepting hsp %d,%d due to diagonal width",lm->set[i]->hsp[j]->query_start,lm->set[i]->hsp[j]->target_start); 

	}

      }
    }	

    if( VERBOSITY_CHECK(3,para->verbosity) ) {
      info("Looking at %s with final score of %d",lm->set[i]->hsp[0]->target->name,new_score);
    }

    lm->set[i]->score = new_score;
  }

  /*  gettimeofday(&t2,NULL); */

  if( VERBOSITY_CHECK(2,para->verbosity) ) {
    info("sorting %d items",lm->len);
  }

  /*qsort(lm->set,lm->len,sizeof(HSPset*),compare_HSPset_score_qsort);*/

  sort_LinearHSPmanager(lm,compare_HSPset_score);


  /*  gettimeofday(&t3,NULL);

  info("truncation clock point: rescoring %f : sorting %f",
	  t2.tv_sec - t1.tv_sec + ((t2.tv_usec - t1.tv_usec) * 1e-6),
	  t3.tv_sec - t2.tv_sec + ((t2.tv_usec - t2.tv_usec) * 1e-6)
	  );

  */

  out = LinearHSPmanager_alloc_std();
  out->mat = hard_link_CompMat(lm->mat);


    for(i=0;i<lm->len;i++) {
      if( VERBOSITY_CHECK(5,para->verbosity) ) {
	info("Accepting hit %s position %d with score %d",lm->set[i]->hsp[0]->target->name,i,lm->set[i]->score);
      }

      add_LinearHSPmanager(out,hard_link_HSPset(lm->set[i]));
      if( i > para->max_size ) {
	break;
      }
    }


  return out;
}



/* Function:  truncated_LinearHSPmanager(lm,max_size,min_score,width,tail)
 *
 * Descrip:    Makes a truncated linear set 
 *
 *
 * Arg:               lm [UNKN ] Undocumented argument [LinearHSPmanager *]
 * Arg:         max_size [UNKN ] Undocumented argument [int]
 * Arg:        min_score [UNKN ] Undocumented argument [int]
 * Arg:            width [UNKN ] Undocumented argument [int]
 * Arg:             tail [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [LinearHSPmanager *]
 *
 */
# line 271 "hsphandler.dy"
LinearHSPmanager * truncated_LinearHSPmanager(LinearHSPmanager * lm,int max_size,int min_score,int width,int tail)
{
  
  LinearHSPmanager * out;
  int i,j;
  int worst_score = UNFEASIBLY_LARGE_SCORE;
  int worst_j = -1;
  HSPset ** worst;
  HSPset * trial;

  assert(max_size > 0);

  out = LinearHSPmanager_alloc_len(max_size+1);
  out->mat = hard_link_CompMat(lm->mat);

  sort_LinearHSPmanager(lm,compare_HSPset_score);

  for(i=0;i<lm->len;i++) {
    /* we chew up into max_size positions */
    if( out->len < max_size ) {
      trial = new_consistent_HSPset(lm->set[i],min_score,width,tail);
      if( trial->len == 0 ) {
	continue;
      }

      add_LinearHSPmanager(out,trial);
    } else {
          /* if have not found worst, scan to find worst */
      if( worst_score == UNFEASIBLY_LARGE_SCORE ) {
	for(j=0;j<out->len;j++) {
	  if( out->set[j]->score < worst_score ) {
	    worst_score = out->set[j]->score;
	    worst = out->set+j;
	    worst_j = j;
	  }
	}
      }
      /* worst_score is now the score */
      if( lm->set[i]->score < worst_score ) {
	break; /* otta here - we don't need to look at any more! */
      }
      trial = new_consistent_HSPset(lm->set[i],min_score,width,tail);
      if( trial->len == 0 ) {
	continue;
      }


      if( trial->score > worst_score ) {
	free_HSPset(out->set[worst_j]);
	out->set[worst_j] = trial;
	worst_score = UNFEASIBLY_LARGE_SCORE;
      } else {
	free_HSPset(trial);
      }
    }
  }

  return out;

}

/* Function:  new_consistent_HSPset(set,min_score,width,tail)
 *
 * Descrip:    Makes an HSP set via heuristics
 *             to deal with low complexity regions
 *
 *
 * Arg:              set [UNKN ] Undocumented argument [HSPset *]
 * Arg:        min_score [UNKN ] Undocumented argument [int]
 * Arg:            width [UNKN ] Undocumented argument [int]
 * Arg:             tail [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [HSPset *]
 *
 */
# line 336 "hsphandler.dy"
HSPset * new_consistent_HSPset(HSPset * set,int min_score,int width,int tail)
{
  HSPset * out;
  int i,j;
  int diag_a;
  int diag_b;
  int eaten;

  /*  fprintf(stderr,"Entering consistency for %s with %d\n",set->hsp[0]->target->name,set->len); */
  out = HSPset_alloc_std();

  sort_HSPset_by_score(set);


  add_HSPset(out,hard_link_HSP(set->hsp[0]));
  out->score = set->hsp[0]->score;
  
  for(i=1;i<set->len;i++) {
    /* check against exisiting HSPs. If fits, add into big set */
    eaten = 0;

    if( set->hsp[i]->score < min_score ) {
      continue;
    }

    for(j=0;j<out->len;j++) {
      diag_a = set->hsp[i]->query_start - set->hsp[i]->target_start;
      diag_b = out->hsp[j]->query_start - out->hsp[j]->target_start;
      if( abs(diag_a - diag_b) > 2 * width ) {
	continue; /* does not match */
      }
      
      eaten = 1;
      add_HSPset(out,hard_link_HSP(set->hsp[i]));
      out->score += set->hsp[i]->score;
      break;
    }
  }

  

  return out;

}

/* Function:  new_LinearHSPmanager_simple_heuristic(hspm,para)
 *
 * Descrip:    New, simpler LinearHSPmanager with diagonal heuristics
 *
 *
 * Arg:        hspm [UNKN ] Undocumented argument [HSPmanager *]
 * Arg:        para [UNKN ] Undocumented argument [LineariseHSPPara *]
 *
 * Return [UNKN ]  Undocumented return value [LinearHSPmanager *]
 *
 */
# line 384 "hsphandler.dy"
LinearHSPmanager * new_LinearHSPmanager_simple_heuristic(HSPmanager * hspm,LineariseHSPPara * para)
{
  int i;
  int j;
  LinearHSPmanager * temp;
  LinearHSPmanager * out;
  
#ifdef LINUX_TIMER
  struct timeval t1;	
  struct timeval t2;
  struct timeval t3;

  if( VERBOSITY_CHECK(1,para->verbosity) ) {
    gettimeofday(&t1,NULL);
  }

#endif

  temp = new_LinearHSPmanager_truncate_on_score(hspm);

  if( VERBOSITY_CHECK(2,para->verbosity) ) {
    info("Have got linear hsp manager with %d entries",temp->len);
  }

  
#ifdef LINUX_TIMER
  if( VERBOSITY_CHECK(1,para->verbosity) ) {
    gettimeofday(&t2,NULL);
  }
#endif
      

  out = truncated_simple_LinearHSPmanager(temp,para);

  if( VERBOSITY_CHECK(2,para->verbosity) ) {
    info("Now got linear hsp manager with %d entries",out->len);
  }

#ifdef LINUX_TIMER


  if( VERBOSITY_CHECK(1,para->verbosity) ) {
    gettimeofday(&t3,NULL);
    info("Sort breakdown: Conversion %f : Truncation %f",
	 t2.tv_sec - t1.tv_sec + ((t2.tv_usec - t1.tv_usec) * 1e-6),
	 t3.tv_sec - t2.tv_sec + ((t2.tv_usec - t2.tv_usec) * 1e-6)
	 );
  }

#endif

  /*
  for(i=0;i<out->len;i++) {
        fprintf(stderr,"%d, got score %d, top score %d %s\n",i,out->set[i]->score,out->set[i]->hsp[0]->score,out->set[i]->hsp[0]->target->name); 
    for(j=0;j<out->set[i]->len;j++) {
       fprintf(stderr,"   HSP %d,%d score %d\n",out->set[i]->hsp[j]->query_start,
	      out->set[i]->hsp[j]->target_start,
	      out->set[i]->hsp[j]->score);
      

    }
  }

  */

  info("Currently not free'ing temporary list");

  return out;
}


/* Function:  new_LinearHSPmanager_heuristic_max(hspm,max_size)
 *
 * Descrip:    Builds a new LinearHSPmanager from hash based with heuristics 
 *
 *
 * Arg:            hspm [UNKN ] Undocumented argument [HSPmanager *]
 * Arg:        max_size [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [LinearHSPmanager *]
 *
 */
# line 458 "hsphandler.dy"
LinearHSPmanager * new_LinearHSPmanager_heuristic_max(HSPmanager * hspm,int max_size)
{
  LinearHSPmanager * temp;
  LinearHSPmanager * out;

  assert(hspm);
  assert(max_size > 0);

  temp = new_LinearHSPmanager_truncate_on_score(hspm);

  out = truncated_LinearHSPmanager(temp,max_size,30,30,40);

  /*  free_LinearHSPmanager(temp); */

  return out;
}

/* Function:  new_LinearHSPmanager_truncate_on_score(hspm)
 *
 * Descrip:    Builds a LinearHSPmanager from a hash based HSP manager, using worst score truncation
 *
 *
 * Arg:        hspm [UNKN ] Undocumented argument [HSPmanager *]
 *
 * Return [UNKN ]  Undocumented return value [LinearHSPmanager *]
 *
 */
# line 478 "hsphandler.dy"
LinearHSPmanager * new_LinearHSPmanager_truncate_on_score(HSPmanager * hspm)
{
  LinearHSPmanager * out;


  out = LinearHSPmanager_alloc_std();
  out->mat = hard_link_CompMat(hspm->mat);

  out->worst_hsp_score = hspm->tsm->worst_score;
  
  /*  for(i=0;i<hspm->tsm->length;i++) {
    fprintf(stderr,"At position %d got score %d\n",i,hspm->tsm->score[i]);
  }
  */
  

  
  /*  fprintf(stderr,"Before management, we have %d, worst score is %d\n",g_hash_table_size(hspm->target_hash),out->worst_hsp_score); */

  


  g_hash_table_foreach(hspm->target_hash,linearise_HSPset_truncate_on_score,out);

/*  sort_LinearHSPmanager(out,compare_HSPset_score); */
    
  return out;

}


/* Function:  new_LinearHSPmanager_flat(hspm)
 *
 * Descrip:    Builds a LinearHSPmanager from a hash based HSP manager
 *
 *
 * Arg:        hspm [UNKN ] Undocumented argument [HSPmanager *]
 *
 * Return [UNKN ]  Undocumented return value [LinearHSPmanager *]
 *
 */
# line 512 "hsphandler.dy"
LinearHSPmanager * new_LinearHSPmanager_flat(HSPmanager * hspm)
{
  LinearHSPmanager * out;

  out = LinearHSPmanager_alloc_std();
  out->mat = hard_link_CompMat(hspm->mat);


  g_hash_table_foreach(hspm->target_hash,linearise_HSPset_flat,out);

  sort_LinearHSPmanager(out,compare_HSPset_score);
    
  return out;

}

/* Function:  linearise_HSPset_truncate_on_score(key,value,user_data)
 *
 * Descrip:    internal function for remapping HSPs
 *
 *
 * Arg:              key [UNKN ] Undocumented argument [gpointer]
 * Arg:            value [UNKN ] Undocumented argument [gpointer]
 * Arg:        user_data [UNKN ] Undocumented argument [gpointer]
 *
 */
# line 531 "hsphandler.dy"
void linearise_HSPset_truncate_on_score(gpointer key,gpointer value,gpointer user_data)
{
  LinearHSPmanager * l = (LinearHSPmanager *) user_data;
  HSPset * s = (HSPset *) value;

  if( s->len == 1 && s->score < l->worst_hsp_score ) {
    return;
  } else {
    add_LinearHSPmanager(l,hard_link_HSPset(s));
  }
}


/* Function:  linearise_HSPset_flat(key,value,user_data)
 *
 * Descrip:    internal function for remapping HSPs with score cutoff
 *
 *
 * Arg:              key [UNKN ] Undocumented argument [gpointer]
 * Arg:            value [UNKN ] Undocumented argument [gpointer]
 * Arg:        user_data [UNKN ] Undocumented argument [gpointer]
 *
 */
# line 547 "hsphandler.dy"
void linearise_HSPset_flat(gpointer key,gpointer value,gpointer user_data)
{
  LinearHSPmanager * l = (LinearHSPmanager *) user_data;
  HSPset * s = (HSPset *) value;
  int i;
  
  s->score = 0;
  s->best_score =0;


  for(i=0;i<s->len;i++) {
    s->score += s->hsp[i]->score;
    if( s->hsp[i]->score > s->best_score ) {
      s->best_score = s->hsp[i]->score;
    }
  }

  add_LinearHSPmanager(l,hard_link_HSPset(s));

}

/* Function:  linearise_HSPset_consistent(key,value,user_data)
 *
 * Descrip:    internal function for remapping HSPs
 *
 *
 * Arg:              key [UNKN ] Undocumented argument [gpointer]
 * Arg:            value [UNKN ] Undocumented argument [gpointer]
 * Arg:        user_data [UNKN ] Undocumented argument [gpointer]
 *
 */
# line 571 "hsphandler.dy"
void linearise_HSPset_consistent(gpointer key,gpointer value,gpointer user_data)
{
  LinearHSPmanager * l = (LinearHSPmanager *) user_data;
  HSPset * s = (HSPset *) value;
  HSPset * add;

  add = new_consistent_HSPset(s,l->min_score,l->width,l->tail);
  if( add->len > 0 ) {
    add_LinearHSPmanager(l,add);
  } else {
    free_HSPset(add);
  }
}



/* Function:  new_HSPmanager(query,mat,score_drop_off)
 *
 * Descrip:    Builds a new HSPmanager for a target system
 *
 *
 * Arg:                 query [UNKN ] Undocumented argument [Sequence *]
 * Arg:                   mat [UNKN ] Undocumented argument [CompMat *]
 * Arg:        score_drop_off [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [HSPmanager *]
 *
 */
# line 590 "hsphandler.dy"
HSPmanager * new_HSPmanager(Sequence * query,CompMat * mat,int score_drop_off)
{
  HSPmanager * out;

  out = HSPmanager_alloc();

  out->query = hard_link_Sequence(query);
  if( mat == NULL ) {
    out->mat = NULL;
  } else {
    out->mat   = hard_link_CompMat(mat);
  }


  out->qs = new_QuerySeqHSP(query,mat); 
  /* out->qs = NULL; */

  out->drop_off = score_drop_off;
  out->target_hash =  g_hash_table_new(g_direct_hash,g_direct_equal);
  out->tsm = new_TopScoreManager(250);
  /*  out->cache = new_HSPCache(2000);*/
  out->cache = NULL;
  out->min_score = 25;
  out->hsp_count = 0;

  return out;
}

/* Function:  add_pair_HSPmanager(hspm,target,query_pos,target_pos)
 *
 * Descrip:    adds a new target pair, irregardless of score
 *
 *
 * Arg:              hspm [UNKN ] Undocumented argument [HSPmanager *]
 * Arg:            target [UNKN ] Undocumented argument [Sequence *]
 * Arg:         query_pos [UNKN ] Undocumented argument [int]
 * Arg:        target_pos [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 621 "hsphandler.dy"
int add_pair_HSPmanager(HSPmanager * hspm,Sequence * target,int query_pos,int target_pos)
{
  return add_pair_HSPmanager_score(hspm,target,query_pos,target_pos,-1);
}

/* Function:  add_pair_HSPmanager_score(hspm,target,query_pos,target_pos,min_score)
 *
 * Descrip:    Adds a new target pair to this HSPmanager for indexing, with a min score
 *
 *
 * Arg:              hspm [UNKN ] Undocumented argument [HSPmanager *]
 * Arg:            target [UNKN ] Undocumented argument [Sequence *]
 * Arg:         query_pos [UNKN ] Undocumented argument [int]
 * Arg:        target_pos [UNKN ] Undocumented argument [int]
 * Arg:         min_score [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 629 "hsphandler.dy"
int add_pair_HSPmanager_score(HSPmanager * hspm,Sequence * target,int query_pos,int target_pos,int min_score)
{
  HSPset * set;
  int i;
  HSP * hsp;
  boolean has_set = 0;

  /* see if this target is loaded into the manager */

  if( (set = g_hash_table_lookup(hspm->target_hash,(gpointer)target)) != NULL ) {
    has_set = 1;
  }

  /* set is now the HSPset. Ensure this position is not already accounted for
     in the set */

  if( has_set == 1 ) {
    if( set->last_accessed != -1 && ON_HSP_MACRO(set->hsp[set->last_accessed],query_pos,target_pos) ) {
      return set->hsp[set->last_accessed]->score;
    }

    for(i=0;i<set->len;i++) {
      if( ON_HSP_MACRO(set->hsp[i],query_pos,target_pos) == TRUE ) {
	set->last_accessed = i;
	return set->hsp[set->last_accessed]->score;
      } 
    }
  } 

  if( hspm->qs != NULL ) {
    hsp = new_HSP_QuerySeqHSP(hspm->cache,hspm->qs,target,query_pos,target_pos,hspm->mat,hspm->drop_off,min_score);
  } else {
    hsp = new_HSP(hspm->cache,hspm->query,target,query_pos,target_pos,hspm->mat,hspm->drop_off);
  }

  if( hsp == NULL ) {
    /* internal min score optimisation */
    return 0;
  }


  if( hsp->score < hspm->min_score ) {
    /* fprintf(stderr,"hsp being lost due to min score %d\n",hsp->score); */
    free_HSP(hsp);
    return 0;
  }

  hspm->hsp_count++;

  if( has_set == 0 )  {
    set = HSPset_alloc_std();
    g_hash_table_insert(hspm->target_hash,(gpointer)target,set);
  }


  set->score += hsp->score;
  if( hsp->score > set->best_score ) {
    set->best_score = hsp->score;
    add_score_TopScoreManager(hspm->tsm,set->best_score);
  }

 
  add_HSPset(set,hsp);
  set->last_accessed = set->len-1;

  return hsp->score;
}

/* Function:  add_new_HSP_HSPmanager(hspm,target,query_start,target_start,length,score)
 *
 * Descrip:    Adds a new HSP when all info is known
 *
 *
 * Arg:                hspm [UNKN ] Undocumented argument [HSPmanager *]
 * Arg:              target [UNKN ] Undocumented argument [Sequence *]
 * Arg:         query_start [UNKN ] Undocumented argument [int]
 * Arg:        target_start [UNKN ] Undocumented argument [int]
 * Arg:              length [UNKN ] Undocumented argument [int]
 * Arg:               score [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 700 "hsphandler.dy"
boolean add_new_HSP_HSPmanager(HSPmanager * hspm,Sequence * target,int query_start,int target_start,int length,int score)
{
  HSPset * set;
  HSP * hsp;

  /* see if this target is loaded into the manager */

  if( (set = g_hash_table_lookup(hspm->target_hash,(gpointer)target)) == NULL ) {
    set = HSPset_alloc_std();
    g_hash_table_insert(hspm->target_hash,(gpointer)target,set);
  }


  if( hspm->cache != NULL ) {
    hsp = HSP_alloc_cache(hspm->cache);
  } else {
    hsp = HSP_alloc();
  }

  hsp->query = hard_link_Sequence(hspm->query);
  hsp->target= hard_link_Sequence(target);
  hsp->score = score;
  hsp->target_start = target_start;
  hsp->query_start = query_start;
  hsp->length = length;


  add_HSPset(set,hsp);

  return TRUE;
}

/* Function:  free_ghash_HSPsets(key,value,user_data)
 *
 * Descrip:    Frees the HSPsets
 *
 *
 * Arg:              key [UNKN ] Undocumented argument [gpointer]
 * Arg:            value [UNKN ] Undocumented argument [gpointer]
 * Arg:        user_data [UNKN ] Undocumented argument [gpointer]
 *
 */
# line 735 "hsphandler.dy"
void free_ghash_HSPsets(gpointer key,gpointer value,gpointer user_data)
{
  int i;
  HSPset * val = (HSPset *) value;
  HSPCache * cache = (HSPCache*) user_data;

  if( val->dynamite_hard_link > 1 ) {
    val->dynamite_hard_link--;
    return;
  }

  if( cache != NULL ) {
    for(i=0;i<val->len;i++) {
      free_HSP_cache(cache,val->hsp[i]);
    }
    val->len = 0;
  }

  free_HSPset(val);
}

/* Function:  free_HSPmanager(h)
 *
 * Descrip:    Frees the HSPmanager
 *
 *
 * Arg:        h [UNKN ] Undocumented argument [HSPmanager *]
 *
 * Return [UNKN ]  Undocumented return value [HSPmanager *]
 *
 */
# line 760 "hsphandler.dy"
HSPmanager * free_HSPmanager(HSPmanager * h)
{
  g_hash_table_foreach(h->target_hash,free_ghash_HSPsets,h->cache);
  g_hash_table_destroy(h->target_hash);
  free_CompMat(h->mat);
  if( h->tsm != NULL ) {
    free_TopScoreManager(h->tsm);
  }
  if( h->qs != NULL ) {
    free_QuerySeqHSP(h->qs);
  }

  ckfree(h);
  return NULL;
}


/* Function:  new_HSP_QuerySeqHSP(cache,query,target,query_pos,target_pos,mat,score_drop_off,min_score)
 *
 * Descrip:    builds a new HSP for these sequences
 *
 *
 * Arg:                 cache [UNKN ] Undocumented argument [HSPCache *]
 * Arg:                 query [UNKN ] Undocumented argument [QuerySeqHSP *]
 * Arg:                target [UNKN ] Undocumented argument [Sequence *]
 * Arg:             query_pos [UNKN ] Undocumented argument [int]
 * Arg:            target_pos [UNKN ] Undocumented argument [int]
 * Arg:                   mat [UNKN ] Undocumented argument [CompMat *]
 * Arg:        score_drop_off [UNKN ] Undocumented argument [int]
 * Arg:             min_score [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [HSP *]
 *
 */
# line 780 "hsphandler.dy"
HSP * new_HSP_QuerySeqHSP(HSPCache * cache,QuerySeqHSP * query,Sequence * target,int query_pos,int target_pos,CompMat * mat,int score_drop_off,int min_score)
{
  int i,j;
  int ii,jj;
  int temp_score;
  int t;
  int pause_i,pause_j;
  int score = 0;
  HSP * out = NULL;
  const char * q_seq;
  const char * t_seq;
  const int * self_score;

  int overall_score;
  int query_start;
  int target_start;

  
  /* we count the start position twice */
  overall_score = -(query->score[query_pos][target->seq[target_pos]-'A']);

  q_seq = query->query->seq;
  t_seq = target->seq;
  self_score = query->self_score;

  pause_i = i = query_pos;
  pause_j = j = target_pos;

  /* go upstream first */
  
  for(score=0;i >= 0 && j >= 0;i--,j--) {

    /* optimise identical matches */
    if( toupper(q_seq[i]) == toupper(t_seq[j]) ) {
      score += self_score[i];
      continue;
    }

    
    t= query->score[i][target->seq[j]-'A'];

    if( t >= 0 ) {
      /* this is a positive score, we are on a high scoring run, so add and move on*/
      score += t;
    } else {
      /* negative score. i+1,j+1 was our last best score */
      for(temp_score = t,ii=i-1,jj=j-1;temp_score >= -score_drop_off && ii >= 0 && jj >= 0 && temp_score < 0;ii--,jj--) {

	temp_score += query->score[ii][target->seq[jj]-'A'];
      }
      
      if( temp_score >= 0 ) {
	/* new maximum reached */
	i = ii;
	j = jj;
	score += temp_score;
	continue; /* back to main loop */
      } else {
	/* either temp_score < -drop_off or something else */
	break;
      }
    }
  }
  
  /* set start position */

  query_start = i+1;
  target_start = j+1;
  overall_score += score;


  /* downstream */


  for(score=0,i=query_pos,j=target_pos;i < query->len && j < target->len;i++,j++) {

    /* optimise identical matches */
    if( toupper(q_seq[i]) == toupper(t_seq[j]) ) {
      score += self_score[i];
      continue;
    }


    t= query->score[i][target->seq[j]-'A'];
    /*    fprintf(stderr,"Doing %d %d, [%c,%c] off %d score %d\n",i,j,query->query->seq[i],target->seq[j],t,score); */
    if( t >= 0 ) {
      /* this is a positive score, we are on a high scoring run, so add and move on*/
      score += t;
    } else {
      /* negative score. i-1,j-1 was our last best score */
      for(temp_score = t,ii=i+1,jj=j+1;temp_score >= -score_drop_off && ii < query->len && jj < target->len && temp_score < 0;ii++,jj++) {
	temp_score += query->score[ii][target->seq[jj]-'A'];
	/*	fprintf(stderr,"Negative looping %d,%d temp %d [%c,%c]\n",ii,jj,temp_score,query->query->seq[ii],target->seq[jj]); */
	
      }
      
      if( temp_score >= 0 ) {
	/* new good score reached */
	i = ii;
	j = jj;
	score += temp_score;
	continue; /* back to main loop */
      } else {
	/* either temp_score < -drop_off or something else */
	break;
      }
    }
  }

  if( overall_score + score < min_score ) {
    return NULL;
  }

  /* we have a valid HSP. Ask for memory and return */


  if( cache != NULL ) {
    out = HSP_alloc_cache(cache);
  } else {
    out = HSP_alloc();
  }

  out->query_start = query_start;
  out->target_start = target_start;
  out->score = overall_score + score;

  out->query = hard_link_Sequence(query->query);
  out->target= hard_link_Sequence(target);
      
  /* set length and score */
  out->length = i-1 - out->query_start +1; 


  return out;
}




/* Function:  new_QuerySeqHSP(seq,mat)
 *
 * Descrip:    Builds a new QuerySeqHSP from Sequence and Matrix
 *
 *
 * Arg:        seq [UNKN ] Undocumented argument [Sequence *]
 * Arg:        mat [UNKN ] Undocumented argument [CompMat *]
 *
 * Return [UNKN ]  Undocumented return value [QuerySeqHSP *]
 *
 */
# line 922 "hsphandler.dy"
QuerySeqHSP * new_QuerySeqHSP(Sequence * seq,CompMat * mat)
{
  QuerySeqHSP * out;
  int i;

  assert(seq);
  assert(mat);

  out = QuerySeqHSP_alloc();
  out->query = hard_link_Sequence(seq);
  out->score = (int **)calloc(seq->len,sizeof(int*));
  out->self_score = (int*) calloc(seq->len,sizeof(int));
  out->len = seq->len;
  
  for(i=0;i<seq->len;i++) {
    out->score[i] = mat->comp[toupper(seq->seq[i])-'A'];
    out->self_score[i] = mat->comp[toupper(seq->seq[i])-'A'][toupper(seq->seq[i])-'A'];
  }

  return out;
}
# line 1047 "hsphandler.c"
/* Function:  hard_link_TopScoreManager(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [TopScoreManager *]
 *
 * Return [UNKN ]  Undocumented return value [TopScoreManager *]
 *
 */
TopScoreManager * hard_link_TopScoreManager(TopScoreManager * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a TopScoreManager object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  TopScoreManager_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TopScoreManager *]
 *
 */
TopScoreManager * TopScoreManager_alloc(void) 
{
    TopScoreManager * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(TopScoreManager *) ckalloc (sizeof(TopScoreManager))) == NULL)  {  
      warn("TopScoreManager_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->score = NULL;   
    out->length = 0; 
    out->current_pos = 0;    
    out->worst_score = 0;    
    out->worst_position = 0; 


    return out;  
}    


/* Function:  free_TopScoreManager(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [TopScoreManager *]
 *
 * Return [UNKN ]  Undocumented return value [TopScoreManager *]
 *
 */
TopScoreManager * free_TopScoreManager(TopScoreManager * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a TopScoreManager obj. Should be trappable");   
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
    if( obj->score != NULL)  
      ckfree(obj->score);    


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_QuerySeqHSP(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [QuerySeqHSP *]
 *
 * Return [UNKN ]  Undocumented return value [QuerySeqHSP *]
 *
 */
QuerySeqHSP * hard_link_QuerySeqHSP(QuerySeqHSP * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a QuerySeqHSP object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  QuerySeqHSP_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [QuerySeqHSP *]
 *
 */
QuerySeqHSP * QuerySeqHSP_alloc(void) 
{
    QuerySeqHSP * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(QuerySeqHSP *) ckalloc (sizeof(QuerySeqHSP))) == NULL)  {  
      warn("QuerySeqHSP_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->query = NULL;   
    out->self_score = NULL;  
    out->len = 0;    


    return out;  
}    


/* Function:  free_QuerySeqHSP(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [QuerySeqHSP *]
 *
 * Return [UNKN ]  Undocumented return value [QuerySeqHSP *]
 *
 */
QuerySeqHSP * free_QuerySeqHSP(QuerySeqHSP * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a QuerySeqHSP obj. Should be trappable");   
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
    if( obj->query != NULL)  
      free_Sequence(obj->query);     
    /* obj->score is linked in */ 
    if( obj->self_score != NULL) 
      ckfree(obj->self_score);   


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_HSPmanager(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [HSPmanager *]
 *
 * Return [UNKN ]  Undocumented return value [HSPmanager *]
 *
 */
HSPmanager * hard_link_HSPmanager(HSPmanager * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a HSPmanager object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  HSPmanager_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [HSPmanager *]
 *
 */
HSPmanager * HSPmanager_alloc(void) 
{
    HSPmanager * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(HSPmanager *) ckalloc (sizeof(HSPmanager))) == NULL)    {  
      warn("HSPmanager_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->query = NULL;   
    out->mat = NULL; 
    out->drop_off = 0;   
    out->min_score = 0;  
    out->tsm = NULL; 
    out->qs = NULL;  
    out->cache = NULL;   
    out->hsp_count = 0;  


    return out;  
}    


/* Function:  hard_link_LineariseHSPPara(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [LineariseHSPPara *]
 *
 * Return [UNKN ]  Undocumented return value [LineariseHSPPara *]
 *
 */
LineariseHSPPara * hard_link_LineariseHSPPara(LineariseHSPPara * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a LineariseHSPPara object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  LineariseHSPPara_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LineariseHSPPara *]
 *
 */
LineariseHSPPara * LineariseHSPPara_alloc(void) 
{
    LineariseHSPPara * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(LineariseHSPPara *) ckalloc (sizeof(LineariseHSPPara))) == NULL)    {  
      warn("LineariseHSPPara_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->max_size = 0;   
    out->min_score = 0;  
    out->width = 0;  
    out->tail = 0;   
    out->verbosity = 0;  


    return out;  
}    


/* Function:  free_LineariseHSPPara(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [LineariseHSPPara *]
 *
 * Return [UNKN ]  Undocumented return value [LineariseHSPPara *]
 *
 */
LineariseHSPPara * free_LineariseHSPPara(LineariseHSPPara * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a LineariseHSPPara obj. Should be trappable");  
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
