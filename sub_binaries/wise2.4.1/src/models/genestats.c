#ifdef _cplusplus
extern "C" {
#endif
#include "genestats.h"


/* Function:  show_help_GeneModelParam(ofp)
 *
 * Descrip:    Shows help
 *
 *
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 97 "genestats.dy"
void show_help_GeneModelParam(FILE * ofp)
{

  fprintf(ofp,"New gene model statistics\n");
  fprintf(ofp,"  -splice_max_collar      [5.0]  maximum Bits value for a splice site \n");
  fprintf(ofp,"  -splice_min_collar      [-5.0] minimum Bits value for a splice site \n");
  fprintf(ofp,"  -splice_score_offset    [%.1f]  score offset for splice sites\n",DEFAULT_SPLICE_OFFSET_SCORE);
  fprintf(ofp,"  -[no]splice_gtag        make just gtag splice sites (default is gtag, ie no model)\n");  
  fprintf(ofp,"  -splice_gtag_prob       [0.001] probability for gt/ag \n");
  fprintf(ofp,"  -genestats              [gene.stat]\n");

}

/* Function:  show_info_GeneModelParam(p,ofp)
 *
 * Descrip:    Shows genemodel param for info
 *
 *
 * Arg:          p [UNKN ] Undocumented argument [GeneModelParam *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 113 "genestats.dy"
void show_info_GeneModelParam(GeneModelParam * p,FILE * ofp)
{
  fprintf(ofp,"Gene Parameter file: %s\n",p->gene_stats_file);
  fprintf(ofp,"Splice site model:   %s\n",p->use_gtag_splice == TRUE ? "GT/AG only" : "Position Weight Matrix");
  if( p->use_gtag_splice == FALSE ) {
    fprintf(ofp,"Pseudo counts 5'SS   %1.2f\n",p->splice5_pseudo);
    fprintf(ofp,"Pseudo counts 3'SS   %1.2f\n",p->splice3_pseudo);
    fprintf(ofp,"Minimum SS collar    %2.2f\n",p->min_collar);
    fprintf(ofp,"Maximum SS collar    %2.2f\n",p->max_collar);
    fprintf(ofp,"Splice site offset   %2.2f\n",p->score_offset);
  } else {

    fprintf(ofp,"GT/AG bits penalty   %2.2f\n",Probability2Bits(p->prob_for_gtag));
  }
  
}

/* Function:  new_GeneModelParam_from_argv(argc,argv)
 *
 * Descrip:    Makes a GeneModelParam from argv
 *
 *
 * Arg:        argc [UNKN ] Undocumented argument [int *]
 * Arg:        argv [UNKN ] Undocumented argument [char **]
 *
 * Return [UNKN ]  Undocumented return value [GeneModelParam *]
 *
 */
# line 133 "genestats.dy"
GeneModelParam * new_GeneModelParam_from_argv(int * argc,char ** argv)
{
  GeneModelParam * out;
  char * temp;

  out = std_GeneModelParam();
  
  if( (temp=strip_out_assigned_argument(argc,argv,"splice_min_collar")) != NULL ) {
    if( is_double_string(temp,&out->min_collar) == FALSE ) {
      warn("%s is not a floating point number. Can't be a splice_min_collar",temp);
      free_GeneModelParam(out);
      return NULL;
    } 
  }

  strip_out_boolean_def_argument(argc,argv,"splice_gtag",&out->use_gtag_splice);

  if( (temp=strip_out_assigned_argument(argc,argv,"splice_max_collar")) != NULL ) {
    if( is_double_string(temp,&out->max_collar) == FALSE ) {
      warn("%s is not a floating point number. Can't be a splice_max_collar",temp);
      free_GeneModelParam(out);
      return NULL;
    } 
  }

  if( (temp=strip_out_assigned_argument(argc,argv,"splice_score_offset")) != NULL ) {
    if( is_double_string(temp,&out->score_offset) == FALSE ) {
      warn("%s is not a floating point number. Can't be a splice_score_offset",temp);
      free_GeneModelParam(out);
      return NULL;
    } 
  }

  if( (temp=strip_out_assigned_argument(argc,argv,"genestats")) != NULL ) {
    if( out->gene_stats_file != NULL ) {
      ckfree(out->gene_stats_file);
    }

    out->gene_stats_file = stringalloc(temp);
  }

  if( (temp=strip_out_assigned_argument(argc,argv,"splice_gtag_prob")) != NULL ) {
    if( is_double_string(temp,&out->prob_for_gtag) == FALSE ) {
      warn("%s is not a floating pointer number. Can't be a probability for gtag",temp);
      free_GeneModelParam(out);
      return NULL;
    }
  }


  return out;

}

/* Function:  vanilla_GeneralGeneModelScore(ct,start_odds,general_odds,stop_odds)
 *
 * Descrip:    Makes a vanilla general gene model score
 *
 *
 * Arg:                  ct [UNKN ] Undocumented argument [CodonTable *]
 * Arg:          start_odds [UNKN ] Undocumented argument [Probability]
 * Arg:        general_odds [UNKN ] Undocumented argument [Probability]
 * Arg:           stop_odds [UNKN ] Undocumented argument [Probability]
 *
 * Return [UNKN ]  Undocumented return value [GeneralGeneModelScore *]
 *
 */
# line 190 "genestats.dy"
GeneralGeneModelScore * vanilla_GeneralGeneModelScore(CodonTable * ct,Probability start_odds,Probability general_odds,Probability stop_odds)
{
  GeneralGeneModelScore * out;
  RandomCodon * temp;
  out = GeneralGeneModelScore_alloc();

  temp = vanilla_start_RandomCodon(ct,start_odds);
  out->start = RandomCodonScore_from_RandomCodon(temp);
  free_RandomCodon(temp);

  temp  = vanilla_stop_RandomCodon(ct,stop_odds);
  out->stop = RandomCodonScore_from_RandomCodon(temp);
  free_RandomCodon(temp);

  temp  = vanilla_general_RandomCodon(ct,general_odds);
  out->general = RandomCodonScore_from_RandomCodon(temp);
  free_RandomCodon(temp);

  return out;
}

/* Function:  vanilla_start_RandomCodon(ct,start_codon_odd_prob)
 *
 * Descrip:    Makes a vanilla (ie, all things equal) Met based
 *             start codon rndscore
 *
 *
 * Arg:                          ct [UNKN ] Undocumented argument [CodonTable *]
 * Arg:        start_codon_odd_prob [UNKN ] Undocumented argument [Probability]
 *
 * Return [UNKN ]  Undocumented return value [RandomCodon *]
 *
 */
# line 215 "genestats.dy"
RandomCodon * vanilla_start_RandomCodon(CodonTable * ct,Probability start_codon_odd_prob)
{
  int i;
  RandomCodon * out;

  out = RandomCodon_alloc();
  
  for(i=0;i<125;i++) {
    if( aminoacid_from_codon(ct,i) == 'M' ) {
      out->codon[i] = start_codon_odd_prob;
    } else {
      out->codon[i] = 0.0;
    }
  }

  return out;
}

/* Function:  vanilla_stop_RandomCodon(ct,stop_codon_odd_prob)
 *
 * Descrip:    Makes a vanilla (ie, all things equal) * based
 *             stop codon rndscore
 *
 *
 * Arg:                         ct [UNKN ] Undocumented argument [CodonTable *]
 * Arg:        stop_codon_odd_prob [UNKN ] Undocumented argument [Probability]
 *
 * Return [UNKN ]  Undocumented return value [RandomCodon *]
 *
 */
# line 237 "genestats.dy"
RandomCodon * vanilla_stop_RandomCodon(CodonTable * ct,Probability stop_codon_odd_prob)
{
  int i;
  RandomCodon * out;

  out = RandomCodon_alloc();
  
  for(i=0;i<125;i++) {
    if( aminoacid_from_codon(ct,i) == 'X' ) {
      out->codon[i] = stop_codon_odd_prob;
    } else {
      out->codon[i] = 0.0;
    }
  }
  return out;
}

/* Function:  vanilla_general_RandomCodon(ct,coding_odd_prob)
 *
 * Descrip:    Makes a vanilla (ie, all things equal) general non-stop
 *
 *
 * Arg:                     ct [UNKN ] Undocumented argument [CodonTable *]
 * Arg:        coding_odd_prob [UNKN ] Undocumented argument [Probability]
 *
 * Return [UNKN ]  Undocumented return value [RandomCodon *]
 *
 */
# line 257 "genestats.dy"
RandomCodon * vanilla_general_RandomCodon(CodonTable * ct,Probability coding_odd_prob)
{
  int i;
  RandomCodon * out;

  out = RandomCodon_alloc();
  
  for(i=0;i<125;i++) {
    if( aminoacid_from_codon(ct,i) == 'X' ) {
      out->codon[i] = 0.0;
    } else {
      out->codon[i] = coding_odd_prob;
    }
  }

  return out;
}


/* Function:  std_GeneModelParam(void)
 *
 * Descrip:    Makes a standard GeneModelParam
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneModelParam *]
 *
 */
# line 279 "genestats.dy"
GeneModelParam * std_GeneModelParam(void)
{
  GeneModelParam * out;

  out = GeneModelParam_alloc();

  out->splice5_pseudo = 1.0;
  out->splice3_pseudo = 1.0;
  out->intron_emission_pseudo = 1.0;
  out->polyp_emission_pseudo  = 1.0;

  out->min_collar   = -5.0;
  out->max_collar   = +5.0;
  out->score_offset = DEFAULT_SPLICE_OFFSET_SCORE;
  out->gene_stats_file = stringalloc("gene.stat");
  out->use_gtag_splice = TRUE;

  out->prob_for_gtag = 0.001;

  return out;
}

/* Function:  GeneModel_from_GeneModelParam(p)
 *
 * Descrip:    Combines GeneStats_from_GeneModelParam and GeneModel_from_GeneStats
 *
 *
 * Arg:        p [UNKN ] Undocumented argument [GeneModelParam *]
 *
 * Return [UNKN ]  Undocumented return value [GeneModel *]
 *
 */
# line 304 "genestats.dy"
GeneModel * GeneModel_from_GeneModelParam(GeneModelParam * p)
{
  GeneStats * st;
  GeneModel * out;

  assert(p);

  st = GeneStats_from_GeneModelParam(p);

  assert(st);

  out = GeneModel_from_GeneStats(st,p);

  assert(out);

  free_GeneStats(st);

  return out;
}


/* Function:  GeneStats_from_GeneModelParam(p)
 *
 * Descrip:    Makes a GeneStats from GeneModelParam - basically just opening the file
 *
 *
 * Arg:        p [UNKN ] Undocumented argument [GeneModelParam *]
 *
 * Return [UNKN ]  Undocumented return value [GeneStats *]
 *
 */
# line 328 "genestats.dy"
GeneStats * GeneStats_from_GeneModelParam(GeneModelParam * p)
{
  GeneStats * gs;
  FILE * ifp;

  assert(p);
  assert(p->gene_stats_file);

  ifp = openfile(p->gene_stats_file,"r");
  if( ifp == NULL ) {
    warn("Unable to open %s  as gene stats file",p->gene_stats_file);
    return NULL;
  }

  gs = read_GeneStats(ifp);

  return gs;
}
  

/* Function:  GeneModel_from_GeneStats(gs,p)
 *
 * Descrip:    Makes a model from the stats file
 *
 *
 * Arg:        gs [UNKN ] Undocumented argument [GeneStats *]
 * Arg:         p [UNKN ] Undocumented argument [GeneModelParam *]
 *
 * Return [UNKN ]  Undocumented return value [GeneModel *]
 *
 */
# line 351 "genestats.dy"
GeneModel * GeneModel_from_GeneStats(GeneStats * gs,GeneModelParam * p)
{
  GeneModel * out;
  int i;
  double total;

  out = GeneModel_alloc();

  assert(gs);
  assert(gs->splice5);
  assert(gs->splice3);
  assert(gs->intron);
  assert(gs->rnd);
  

  for(i=0;i<64;i++) {
     out->codon[i] = gs->codon[i];
  }
  
  out->splice5 = pwmDNA_from_SeqAlign(gs->splice5,p->splice5_pseudo);
  
  /*  fprintf(stdout,"GS splice5 %d splice3 %d\n",gs->splice5,gs->splice3);*/

  fold_randommodel_pwmDNA(out->splice5,gs->rnd);

  out->splice5score = SpliceSiteScore_alloc();
  out->splice5score->score = pwmDNAScore_from_pwmDNA(out->splice5);
  out->splice5score->offset = gs->splice5_offset;
  out->splice5score->min_collar   = Probability2Score(Bits2Probability(p->min_collar));
  out->splice5score->max_collar   = Probability2Score(Bits2Probability(p->max_collar));
  out->splice5score->score_offset = Probability2Score(Bits2Probability(p->score_offset));
 

  out->splice3 = pwmDNA_from_SeqAlign(gs->splice3,p->splice3_pseudo);
  fold_randommodel_pwmDNA(out->splice3,gs->rnd);

  out->splice3score = SpliceSiteScore_alloc();
  out->splice3score->score = pwmDNAScore_from_pwmDNA(out->splice3);
  out->splice3score->offset = gs->splice3_offset;
  out->splice3score->min_collar   = Probability2Score(Bits2Probability(p->min_collar));
  out->splice3score->max_collar   = Probability2Score(Bits2Probability(p->max_collar));
  out->splice3score->score_offset = Probability2Score(Bits2Probability(p->score_offset));

  out->use_gtag_splice = p->use_gtag_splice;
  out->score_for_gtag  = Probability2Score(p->prob_for_gtag);

  out->intron  = RandomModelDNA_alloc();
  for(total = 0.0,i=0;i<4;i++)
    total += gs->intron->base[i] + p->intron_emission_pseudo;

  for(i=0;i<4;i++)
    out->intron->base[i] = (gs->intron->base[i] + p->intron_emission_pseudo)/total;

  out->intron->base[4] = 1.0;

  if( gs->polyp != NULL ) {
    out->polyp  = RandomModelDNA_alloc();
    for(total = 0.0,i=0;i<4;i++)
     total += gs->polyp->base[i] + p->polyp_emission_pseudo;

    for(i=0;i<4;i++)
       out->polyp->base[i] = (gs->polyp->base[i] + p->polyp_emission_pseudo)/total;
  }

  out->rnd = hard_link_RandomModelDNA(gs->rnd);
  return out;
}


/* Function:  show_GeneModel(gm,ofp)
 *
 * Descrip:    shows a genemodel
 *
 *
 * Arg:         gm [UNKN ] Undocumented argument [GeneModel *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 423 "genestats.dy"
void show_GeneModel(GeneModel * gm,FILE * ofp)
{
  fprintf(ofp,"Splice5\n");
  show_pwmDNA_col(gm->splice5,ofp);
  fprintf(ofp,"Splice3\n");
  show_pwmDNA_col(gm->splice3,ofp);
}

/* Function:  new_ComplexSequenceEvalSet_from_GeneModel(gm)
 *
 * Descrip:    Makes an entire ComplexSequenceEvalSet for genomic work
 *
 *
 * Arg:        gm [UNKN ] Undocumented argument [GeneModel *]
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequenceEvalSet *]
 *
 */
# line 434 "genestats.dy"
ComplexSequenceEvalSet * new_ComplexSequenceEvalSet_from_GeneModel(GeneModel * gm)
{
  ComplexSequenceEvalSet * out;

  assert(gm);
  assert(gm->splice5score);
  assert(gm->splice3score);

  out = ComplexSequenceEvalSet_alloc_len(11);

  add_ComplexSequenceEvalSet(out,base_number_ComplexSequenceEval());
  add_ComplexSequenceEvalSet(out,codon_number_ComplexSequenceEval());
  if( gm->use_gtag_splice == FALSE ) {
    add_ComplexSequenceEvalSet(out,ComplexSequenceEval_from_pwmDNAScore_splice(gm->splice5score));
    add_ComplexSequenceEvalSet(out,ComplexSequenceEval_from_pwmDNAScore_splice(gm->splice3score));
  } else {
    add_ComplexSequenceEvalSet(out,ComplexSequenceEval_for_scored_gt(&gm->score_for_gtag));
    add_ComplexSequenceEvalSet(out,ComplexSequenceEval_for_scored_ag(&gm->score_for_gtag));
  }

  add_ComplexSequenceEvalSet(out,flat_zero());
  add_ComplexSequenceEvalSet(out,flat_zero());


  out->type = SEQUENCE_GENOMIC;

  prepare_ComplexSequenceEvalSet(out);

  return out;
}


/* Function:  ComplexSequenceEval_for_scored_ag(score_for_ag)
 *
 * Descrip:    Makes a ComplexSequenceEval for a GT rule
 *
 *
 * Arg:        score_for_ag [UNKN ] Undocumented argument [Score *]
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequenceEval *]
 *
 */
# line 470 "genestats.dy"
ComplexSequenceEval * ComplexSequenceEval_for_scored_ag(Score * score_for_ag)
{
  ComplexSequenceEval * out;

  out = ComplexSequenceEval_alloc();

  out->left_window   = 3;
  out->right_window  = 3;
  out->left_lookback = 5;
  out->outside_score = NEGI;

  out->data = (void*) score_for_ag;
  out->type = SEQUENCE_GENOMIC;
  out->eval_func = scored_ag_eval_func;
  out->score_type = CseScoreType_Bits;

  return out;
}

/* Function:  scored_ag_eval_func(type,*data,seq)
 *
 * Descrip:    Function which actually does the evaluation for scored doners
 *
 *
 * Arg:         type [UNKN ] Undocumented argument [int]
 * Arg:        *data [UNKN ] Undocumented argument [void]
 * Arg:          seq [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 494 "genestats.dy"
int scored_ag_eval_func(int type,void *data,char * seq)
{
  if( *(seq-1) == 'A' && *(seq) == 'G' ) 
    return *(Score *)data;
  else return NEGI;
}

/* Function:  ComplexSequenceEval_for_scored_gt(score_for_gt)
 *
 * Descrip:    Makes a ComplexSequenceEval for a GT rule
 *
 *
 * Arg:        score_for_gt [UNKN ] Undocumented argument [Score *]
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequenceEval *]
 *
 */
# line 505 "genestats.dy"
ComplexSequenceEval * ComplexSequenceEval_for_scored_gt(Score * score_for_gt)
{
  ComplexSequenceEval * out;

  out = ComplexSequenceEval_alloc();

  out->left_window   = 3;
  out->right_window  = 7;
  out->left_lookback = 8;
  out->outside_score = NEGI;

  out->data = (void*) score_for_gt;
  out->type = SEQUENCE_GENOMIC;
  out->eval_func = scored_gt_eval_func;
  out->score_type = CseScoreType_Bits;

  return out;
}

/* Function:  scored_gt_eval_func(type,*data,seq)
 *
 * Descrip:    Function which actually does the evaluation for scored doners
 *
 *
 * Arg:         type [UNKN ] Undocumented argument [int]
 * Arg:        *data [UNKN ] Undocumented argument [void]
 * Arg:          seq [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 529 "genestats.dy"
int scored_gt_eval_func(int type,void *data,char * seq)
{
  if( *(seq) == 'G' && *(seq+1) == 'T' ) 
    return *(Score *)data;
  else return NEGI;
}


/* Function:  ComplexSequenceEval_from_pwmDNAScore_splice(score)
 *
 * Descrip:    Makes a ComplexSequenceEval for a splice site
 *             pwmdna
 *
 *
 * Arg:        score [UNKN ] Undocumented argument [SpliceSiteScore *]
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequenceEval *]
 *
 */
# line 542 "genestats.dy"
ComplexSequenceEval * ComplexSequenceEval_from_pwmDNAScore_splice(SpliceSiteScore * score)
{
  ComplexSequenceEval * out;


  /*  printf("Making CSE from %d\n",score);*/

  out = ComplexSequenceEval_alloc();

  /* shouldn't really add ones, but this is ok anyway. 
     Yukky hack due to not understanding a bug in the window 
     determination
     */

  /**
   *STILL don't know precisely what is going on down here! ***/

  /*  out->left_window  = ssm->offset + ssm->pre_splice_site +1; */
  out->left_window   =10;
  /*  out->right_window = ssm->offset + ssm->post_splice_site +1; */
  out->right_window  =10;
  out->left_lookback =10;

  out->outside_score= NEGI;
  out->data_type    = 245; /* any old key */
  out->data         = (void *) score;
  out->type         = SEQUENCE_GENOMIC;
  out->eval_func    = pwmDNA_splice_ComplexSequence_eval_func;
  out->score_type   = CseScoreType_Bits;

  return out;
}

/* Function:  pwmDNA_splice_ComplexSequence_eval_func(type,data,seq)
 *
 * Descrip:    This function is used as a pointer to function in the eval func
 *
 *             You should never be using this function yourself!
 *
 *
 * Arg:        type [UNKN ] Undocumented argument [int]
 * Arg:        data [UNKN ] Undocumented argument [void *]
 * Arg:         seq [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 581 "genestats.dy"
int pwmDNA_splice_ComplexSequence_eval_func(int type,void * data,char * seq)
{
  SpliceSiteScore * sc;
  pwmDNAScore * pds;
  int score;

  sc  = (SpliceSiteScore* ) data;
  pds = sc->score;

  /* offset is written in biological coordinates. Need to get c style coordiates */

  /* no idea what is happening here, but it works ;) */


  if( seq[0] == 'N' && seq[1] == 'N' && seq[2] == 'N' ) {
    return NEGI;
  }

  score = score_pwmDNAScore_string(pds,seq-sc->offset+1);

/*  printf("Offset is %d %c%c%c\n",sc->offset,seq[0],seq[1],seq[2]);*/

/*  fprintf(stdout,"Before collaring, %d\n",score);*/

  if( score < sc->min_collar ) {
    score = sc->min_collar;
  }
  if( score > sc->max_collar ) {
    score = sc->max_collar;
  }

/*  fprintf(stdout,"Score %d before offset\n",score);*/

  score -= sc->score_offset;

/*  fprintf(stdout,"Score %d after offset\n",score);*/

  /*  fprintf(stderr,"Scoring %d at some position\n",score);*/

  return score;
}

/* Function:  read_GeneStats(ifp)
 *
 * Descrip:    Reads a GeneStats file
 *
 *
 * Arg:        ifp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [GeneStats *]
 *
 */
# line 626 "genestats.dy"
GeneStats * read_GeneStats(FILE * ifp)
{
  char buffer[MAXLINE];
  GeneStats * out;
  SeqAlign * temp;
  char **base;
  char **brk;

  out = GeneStats_alloc();
  out->rnd = NULL;

  while( fgets(buffer,MAXLINE,ifp) != NULL ) {
    /*    fprintf(stderr,"Reading (main loop) %s",buffer); */
    if( buffer[0] == '#' )
      continue;
    
    if( buffer[0] == '%' && buffer[1] == '%' )
      break;
    
    if( strstartcmp(buffer,"splice5") == 0 ) {
      base = brk = breakstring(buffer,spacestr);
      if( *brk == NULL || *(brk+1) == NULL || is_integer_string(*(brk+1),&out->splice5_offset) == 0) {
	warn("Cannot read splice5 offset - must be splice5 <number>");
	return NULL;
      }
      ckfree(base);
      temp = read_selex_SeqAlign(ifp);
      if( temp == NULL ) {
	warn("Could not read in selex alignment for splice5");
	continue;
      }

      out->splice5 = temp;
      continue;
    }

    if( strstartcmp(buffer,"splice3") == 0 ) {

      base = brk = breakstring(buffer,spacestr);
      if( *brk == NULL || *(brk+1) == NULL || is_integer_string(*(brk+1),&out->splice3_offset) == 0) {
	warn("Cannot read splice3 offset - must be splice3 <number>");
	return NULL;
      }
      ckfree(base);

      temp = read_selex_SeqAlign(ifp);
      if( temp == NULL ) {
	warn("Could not read in selex alignment for splice5");
	continue;
      }

      out->splice3 = temp;

      continue;
    }
    
    
    if( strstartcmp(buffer,"intron_emission") == 0 ) {
      if( fgets(buffer,MAXLINE,ifp) == NULL ) {
	warn("Could not read in intron emission line");
	break;
      }
      out->intron = get_genestat_emission(buffer);
      if( fgets(buffer,MAXLINE,ifp) != NULL ) {
	continue;
      } else {
	break;
      }
    }

    if( strstartcmp(buffer,"polyp_emission") == 0 ) {
      if( fgets(buffer,MAXLINE,ifp) == NULL ) {
	warn("Could not read in polyp emission line");
	break;
      }
      out->polyp = get_genestat_emission(buffer);
      if( fgets(buffer,MAXLINE,ifp) != NULL ) {
	continue;
      } else {
	break;
      }
    }

    if( strstartcmp(buffer,"rnd_emission") == 0 ) {
      if( fgets(buffer,MAXLINE,ifp) == NULL ) {
	warn("Could not read in rnd emission line");
	break;
      }
      out->rnd = get_genestat_emission(buffer);
      if( fgets(buffer,MAXLINE,ifp) != NULL ) {
	continue;
      } else {
	break;
      }
    }

    if( strstartcmp(buffer,"rndcodon") == 0 ) {
      if( read_codon_GeneStats(out->codon,buffer,ifp) == FALSE ) {
	warn("Problem in reading codon line!");
      }
      continue;
    }

    if( isalpha(buffer[0]) ) { 
      warn("Could not read line %s in genestats reading\n",buffer);
    }
  
  }

  
  assert(out);
  assert(out->splice5);
  assert(out->splice3);
	
  return out;
}

/* Function:  get_genestat_emission(buffer)
 *
 * Descrip:    reads in the emission stuff in a genestats line
 *
 *
 * Arg:        buffer [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [RandomModelDNA *]
 *
 */
# line 747 "genestats.dy"
RandomModelDNA * get_genestat_emission(char * buffer)
{
  RandomModelDNA * out;
  int i;
  char ** base;
  char ** brk;
  double d;

  out = RandomModelDNA_alloc();

  base = brk = breakstring(buffer,spacestr);

  for(i=0;*brk != NULL && i < 5;i++, brk++){
    if( is_double_string(*brk,out->base+i) == FALSE) {
      warn("For genestat word %s, not a double in emission!",*brk);
      return FALSE;
    }
  }

  ckfree(base);

  if( i < 4 ) {
    warn("Did not read in 5 numbers for emission scores in genestats");
  }

  return out;
}

/* Function:  dump_GeneStats(st,ofp)
 *
 * Descrip:    testing function
 *
 *
 * Arg:         st [UNKN ] Undocumented argument [GeneStats *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 778 "genestats.dy"
void dump_GeneStats(GeneStats * st,FILE * ofp)
{
  int i;

  assert(st);
  assert(ofp);

  fprintf(ofp,"#\n# Dumping gene stats, wise2.2 style\n#\n");
  fprintf(ofp,"splice5\n");
  write_selex_SeqAlign(st->splice5,10,70,ofp);
  fprintf(ofp,"//\nsplice3\n");
  write_selex_SeqAlign(st->splice3,10,70,ofp);
  fprintf(ofp,"//\n");

  fprintf(ofp,"intron_emission\n");
  for(i=0;i<4;i++) {
    fprintf(ofp,"%f ",st->intron->base[i]);
  }

  fprintf(ofp,"\n");
  fprintf(ofp,"//\n");
  if( st->polyp != NULL ) {
    fprintf(ofp,"polyp_emission\n");
    for(i=0;i<4;i++) {
      fprintf(ofp,"%f ",st->polyp->base[i]);
    }
  }

  fprintf(ofp,"\n");
  fprintf(ofp,"//\n");

}


/* Function:  read_codon_GeneStats(codon_array,line,ifp)
 *
 * Descrip:    assummes codon_array is 64 positions long
 *               
 *             line should have begin consensus on it and be of MAXLINE length as it will be used as the buffer.
 *
 *             This does **not** check that you have filled up all 64 positions.
 *
 *
 * Arg:        codon_array [UNKN ] Undocumented argument [double *]
 * Arg:               line [UNKN ] Undocumented argument [char*]
 * Arg:                ifp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 819 "genestats.dy"
boolean read_codon_GeneStats(double * codon_array,char* line,FILE * ifp)
{
  boolean ret = TRUE;
  char * codon;
  char * number;


  if( strwhitestartcmp(line,"rndcodon",spacestr) != 0  ) {
    warn("In reading codon line, got no 'rndcoodon' tag [%s]",line);
    return FALSE;
  }


  while( fgets(line,MAXLINE,ifp) != NULL ) {
    if( line[0] == '#' )
      continue;

    if( strwhitestartcmp(line,"//",spacestr) == 0 )
      break;

    codon = strtok(line,spacestr);
    number = strtok(NULL,spacestr);

    if( codon == NULL ) {
      warn("Found an uncommented line in codon consensus with no leading codon word");
      continue;
    }

    if( number == NULL ) {
      warn("For codon %s, no number found",codon);
      ret = FALSE;
      continue;
    }

    if( strchr(codon,'N') != NULL ) 
      continue;

    if( is_non_ambiguous_codon_seq(codon) == FALSE ) {
      warn("Codon %s is not really a codon... problem!");
      ret = FALSE;
      continue;
    }



    codon_array[base4_codon_from_seq(codon)]= atof(number);

  }

  return ret;
}


# line 923 "genestats.c"
/* Function:  hard_link_SpliceSiteScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SpliceSiteScore *]
 *
 * Return [UNKN ]  Undocumented return value [SpliceSiteScore *]
 *
 */
SpliceSiteScore * hard_link_SpliceSiteScore(SpliceSiteScore * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a SpliceSiteScore object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  SpliceSiteScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SpliceSiteScore *]
 *
 */
SpliceSiteScore * SpliceSiteScore_alloc(void) 
{
    SpliceSiteScore * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(SpliceSiteScore *) ckalloc (sizeof(SpliceSiteScore))) == NULL)  {  
      warn("SpliceSiteScore_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->score = NULL;   
    out->offset = 0; 
    out->min_collar = 0; 
    out->max_collar = 0; 
    out->score_offset = 0;   


    return out;  
}    


/* Function:  free_SpliceSiteScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SpliceSiteScore *]
 *
 * Return [UNKN ]  Undocumented return value [SpliceSiteScore *]
 *
 */
SpliceSiteScore * free_SpliceSiteScore(SpliceSiteScore * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a SpliceSiteScore obj. Should be trappable");   
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
      free_pwmDNAScore(obj->score);  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_GeneStats(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GeneStats *]
 *
 * Return [UNKN ]  Undocumented return value [GeneStats *]
 *
 */
GeneStats * hard_link_GeneStats(GeneStats * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a GeneStats object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  GeneStats_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneStats *]
 *
 */
GeneStats * GeneStats_alloc(void) 
{
    GeneStats * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GeneStats *) ckalloc (sizeof(GeneStats))) == NULL)  {  
      warn("GeneStats_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->splice5 = NULL; 
    out->splice5_offset = 0; 
    out->splice3 = NULL; 
    out->splice3_offset = 0; 
    out->intron = NULL;  
    out->average_intron = 0; 
    out->polyp = NULL;   
    out->average_polyp = 0;  
    out->rnd = NULL; 
    /* codon[64] is an array: no default possible */ 


    return out;  
}    


/* Function:  free_GeneStats(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneStats *]
 *
 * Return [UNKN ]  Undocumented return value [GeneStats *]
 *
 */
GeneStats * free_GeneStats(GeneStats * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a GeneStats obj. Should be trappable"); 
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
    if( obj->splice5 != NULL)    
      free_SeqAlign(obj->splice5);   
    if( obj->splice3 != NULL)    
      free_SeqAlign(obj->splice3);   
    if( obj->intron != NULL) 
      free_RandomModelDNA(obj->intron);  
    if( obj->polyp != NULL)  
      free_RandomModelDNA(obj->polyp);   
    if( obj->rnd != NULL)    
      free_RandomModelDNA(obj->rnd);     


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_GeneModelParam(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GeneModelParam *]
 *
 * Return [UNKN ]  Undocumented return value [GeneModelParam *]
 *
 */
GeneModelParam * hard_link_GeneModelParam(GeneModelParam * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a GeneModelParam object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  GeneModelParam_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneModelParam *]
 *
 */
GeneModelParam * GeneModelParam_alloc(void) 
{
    GeneModelParam * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GeneModelParam *) ckalloc (sizeof(GeneModelParam))) == NULL)    {  
      warn("GeneModelParam_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->splice5_pseudo = 0; 
    out->splice3_pseudo = 0; 
    out->intron_emission_pseudo = 0; 
    out->polyp_emission_pseudo = 0;  
    out->min_collar = 0; 
    out->max_collar = 0; 
    out->score_offset = 0;   
    out->gene_stats_file = NULL; 
    out->use_gtag_splice = FALSE;    
    out->prob_for_gtag = 0;  


    return out;  
}    


/* Function:  free_GeneModelParam(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneModelParam *]
 *
 * Return [UNKN ]  Undocumented return value [GeneModelParam *]
 *
 */
GeneModelParam * free_GeneModelParam(GeneModelParam * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a GeneModelParam obj. Should be trappable");    
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
    if( obj->gene_stats_file != NULL)    
      ckfree(obj->gene_stats_file);  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_GeneModel(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GeneModel *]
 *
 * Return [UNKN ]  Undocumented return value [GeneModel *]
 *
 */
GeneModel * hard_link_GeneModel(GeneModel * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a GeneModel object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  GeneModel_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneModel *]
 *
 */
GeneModel * GeneModel_alloc(void) 
{
    GeneModel * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GeneModel *) ckalloc (sizeof(GeneModel))) == NULL)  {  
      warn("GeneModel_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->splice5 = NULL; 
    out->splice5_offset = 0; 
    out->splice3 = NULL; 
    out->splice3_offset = 0; 
    out->intron = NULL;  
    out->intron_stay_prob = 0;   
    out->polyp = NULL;   
    out->polyp_stay_prob = 0;    
    out->rnd = NULL; 
    out->splice5score = NULL;    
    out->splice3score = NULL;    
    /* codon[64] is an array: no default possible */ 
    out->use_gtag_splice = FALSE;    
    out->score_for_gtag = 0; 


    return out;  
}    


/* Function:  free_GeneModel(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneModel *]
 *
 * Return [UNKN ]  Undocumented return value [GeneModel *]
 *
 */
GeneModel * free_GeneModel(GeneModel * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a GeneModel obj. Should be trappable"); 
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
    if( obj->splice5 != NULL)    
      free_pwmDNA(obj->splice5);     
    if( obj->splice3 != NULL)    
      free_pwmDNA(obj->splice3);     
    if( obj->intron != NULL) 
      free_RandomModelDNA(obj->intron);  
    if( obj->polyp != NULL)  
      free_RandomModelDNA(obj->polyp);   
    if( obj->rnd != NULL)    
      free_RandomModelDNA(obj->rnd);     
    if( obj->splice5score != NULL)   
      free_SpliceSiteScore(obj->splice5score);   
    if( obj->splice3score != NULL)   
      free_SpliceSiteScore(obj->splice3score);   


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_GeneralGeneModelScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GeneralGeneModelScore *]
 *
 * Return [UNKN ]  Undocumented return value [GeneralGeneModelScore *]
 *
 */
GeneralGeneModelScore * hard_link_GeneralGeneModelScore(GeneralGeneModelScore * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a GeneralGeneModelScore object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  GeneralGeneModelScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneralGeneModelScore *]
 *
 */
GeneralGeneModelScore * GeneralGeneModelScore_alloc(void) 
{
    GeneralGeneModelScore * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GeneralGeneModelScore *) ckalloc (sizeof(GeneralGeneModelScore))) == NULL)  {  
      warn("GeneralGeneModelScore_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->start = NULL;   
    out->stop = NULL;    
    out->general = NULL; 


    return out;  
}    


/* Function:  free_GeneralGeneModelScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneralGeneModelScore *]
 *
 * Return [UNKN ]  Undocumented return value [GeneralGeneModelScore *]
 *
 */
GeneralGeneModelScore * free_GeneralGeneModelScore(GeneralGeneModelScore * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a GeneralGeneModelScore obj. Should be trappable"); 
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
    if( obj->start != NULL)  
      free_RandomCodonScore(obj->start);     
    if( obj->stop != NULL)   
      free_RandomCodonScore(obj->stop);  
    if( obj->general != NULL)    
      free_RandomCodonScore(obj->general);   


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
