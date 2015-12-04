#ifdef _cplusplus
extern "C" {
#endif
#include "genewisemodel.h"


/* Function:  pack_GeneWiseScore(gws)
 *
 * Descrip:    Packing up the GeneWise model into a byte structure
 *
 *
 * Arg:        gws [UNKN ] Undocumented argument [GeneWiseScore *]
 *
 * Return [UNKN ]  Undocumented return value [char *]
 *
 */
# line 110 "genewisemodel.dy"
char * pack_GeneWiseScore(GeneWiseScore * gws)
{
  char * out;

  return out;
}



/* Function:  GeneWiseScoreFlat_from_GeneWiseScore(gws)
 *
 * Descrip:    This produces a flattened GeneWiseSegment structure
 *             for use in quick implementations (memory lookup
 *             is much better due to everything being a single
 *             piece of memory).
 *
 *
 * Arg:        gws [UNKN ] Undocumented argument [GeneWiseScore *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseScoreFlat *]
 *
 */
# line 125 "genewisemodel.dy"
GeneWiseScoreFlat * GeneWiseScoreFlat_from_GeneWiseScore(GeneWiseScore * gws)
{
  int i;
  int j;
  GeneWiseScoreFlat * gwsf;

  gwsf = GeneWiseScoreFlat_alloc();
  gwsf->seg = (GeneWiseScoreSegment *) malloc (sizeof(GeneWiseScoreSegment) * gws->len);
  for(i=0;i<gws->len;i++) {
    for(j=0;j<126;j++) {
      gwsf->seg[i].match[j] = gws->seg[i]->match[j];
      gwsf->seg[i].insert[j] = gws->seg[i]->insert[j];
    }

    for(j=0;j<GW_TRANSITION_LEN;j++)
      gwsf->seg[i].transition[j] = gws->seg[i]->transition[j];
      
  }
  gwsf->len = gws->len;
  return gwsf;
}

/* Function:  free_GeneWiseScoreFlat(obj)
 *
 * Descrip:    Frees the GeneWiseScoreFlat datastructure
 *
 *             overrides the usual deconstructor
 *
 *
 * Arg:        obj [UNKN ] Undocumented argument [GeneWiseScoreFlat *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseScoreFlat *]
 *
 */
# line 153 "genewisemodel.dy"
GeneWiseScoreFlat * free_GeneWiseScoreFlat(GeneWiseScoreFlat * obj)
{
  if( obj == NULL ) {
    warn("Attempting to free a NULL GeneWiseScoreFlat object");
    return NULL;
  }

  if( obj->seg != NULL )
    ckfree(obj->seg);
  ckfree(obj);
  
  return NULL;
}


/* Function:  map_phase0_codons_AlnBlock_GeneWise(alb,gws,cseq)
 *
 * Descrip:    This function does something very
 *             sinister.
 *
 *             It maps the phase 0 introns which 
 *             have three base pairs added on the
 *             end. 
 *
 *             It actually changes the AlnBlock structure
 *
 *
 * Arg:         alb [UNKN ] Undocumented argument [AlnBlock *]
 * Arg:         gws [UNKN ] Undocumented argument [GeneWiseScore *]
 * Arg:        cseq [UNKN ] Undocumented argument [ComplexSequence *]
 *
 */
# line 178 "genewisemodel.dy"
void map_phase0_codons_AlnBlock_GeneWise(AlnBlock * alb,GeneWiseScore * gws,ComplexSequence * cseq)
{
  AlnColumn * alc;
  AlnColumn * new;

  for(alc=alb->start;alc != NULL;alc=alc->next) {
    /* don't map random columns */
    if( is_random_AlnColumn_genewise(alc) == TRUE )
      continue;
    if( strcmp(alc->alu[1]->text_label,"3SS_PHASE_0") == 0 ) {
      new = new_pairwise_AlnColumn();
      new->alu[0]->start = alc->alu[0]->start;
      new->alu[0]->end = alc->alu[0]->end;
      
      /* set old alc end to == start */

      alc->alu[0]->end = new->alu[0]->start;

      /* link label in as the same */

      new->alu[0]->text_label = alc->alu[0]->text_label;


      new->alu[1]->end = alc->alu[1]->end;

      /* new starts at -3 on the codon, which is where old ends */

      new->alu[1]->start = alc->alu[1]->end = alc->alu[1]->end -3;

      new->alu[1]->text_label = "CODON";
      
      new->alu[0]->score[0] = new->alu[1]->score[0] = gws->seg[new->alu[0]->end]->transition[GW_MATCH2MATCH] + gws->seg[new->alu[0]->end]->match[CSEQ_GENOMIC_CODON(cseq,new->alu[1]->end)];

      alc->alu[0]->score[0] -= new->alu[0]->score[0];
      alc->alu[1]->score[0] -= new->alu[0]->score[0];


      /*** add new into the AlnBlock ***/

      new->next = alc->next;
      alc->next = new;
    }
  }
}


/* Function:  flatten_balance_scores_GeneWise(gw)
 *
 * Descrip:    This function is make all the balance scores 0 (hence prob-ratio to 1).
 *
 *             Used when you are using naive models
 *
 *
 * Arg:        gw [UNKN ] genewise model to flatten [GeneWise *]
 *
 */
# line 231 "genewisemodel.dy"
void flatten_balance_scores_GeneWise(GeneWise * gw)
{
  int i;

  for(i=0;i<gw->len;i++) {
    gw->seg[i]->transition[GW_MATCH_BALANCE_5SS] = 1;
    gw->seg[i]->transition[GW_INSERT_BALANCE_5SS] = 1;
    gw->seg[i]->transition[GW_MATCH_BALANCE_3SS] = 1;
    gw->seg[i]->transition[GW_INSERT_BALANCE_3SS] = 1;
  }
}


/* Function:  GeneWise_from_ThreeStateModel_cdna(tsm,cp,cm,allN)
 *
 * Descrip:    This function makes a 
 *             GeneWise model for the estwise
 *             type algorithms
 *
 *
 * Arg:         tsm [UNKN ] Undocumented argument [ThreeStateModel *]
 * Arg:          cp [UNKN ] Undocumented argument [cDNAParser *]
 * Arg:          cm [UNKN ] Undocumented argument [CodonMapper *]
 * Arg:        allN [UNKN ] Undocumented argument [Probability]
 *
 * Return [UNKN ]  Undocumented return value [GeneWise *]
 *
 */
# line 249 "genewisemodel.dy"
GeneWise * GeneWise_from_ThreeStateModel_cdna(ThreeStateModel * tsm,cDNAParser * cp,CodonMapper * cm,Probability allN)
{
  register int i;
  GeneWise * out;
  Probability factor;

  assert(tsm);
  assert(cp);
  assert(cm);

  if( (out = GeneWise_alloc_len(tsm->len)) == NULL )
    return NULL;


  factor = (1.0 - removed_probability_from_cds_cdna(cp));

  for(i=0;i<tsm->len;i++) 
    add_GeneWise(out,GeneWiseSegment_from_ThreeStateUnit(tsm->unit[i],factor,cm,NULL,allN));


  return out;
}

/* Function:  GeneWise_from_ThreeStateModel_setfactor(tsm,factor,cm,allN)
 *
 * Descrip:    This function makes a 
 *             GeneWise model for the estwise
 *             type algorithms
 *
 *
 * Arg:           tsm [UNKN ] Undocumented argument [ThreeStateModel *]
 * Arg:        factor [UNKN ] Undocumented argument [Probability]
 * Arg:            cm [UNKN ] Undocumented argument [CodonMapper *]
 * Arg:          allN [UNKN ] Undocumented argument [Probability]
 *
 * Return [UNKN ]  Undocumented return value [GeneWise *]
 *
 */
# line 277 "genewisemodel.dy"
GeneWise * GeneWise_from_ThreeStateModel_setfactor(ThreeStateModel * tsm,Probability factor,CodonMapper * cm,Probability allN)
{
  register int i;
  GeneWise * out;

  assert(tsm  != NULL);
  assert(cm != NULL);

  if( (out = GeneWise_alloc_len(tsm->len)) == NULL )
    return NULL;

  for(i=0;i<tsm->len;i++) 
    add_GeneWise(out,GeneWiseSegment_from_ThreeStateUnit(tsm->unit[i],factor,cm,NULL,allN));


  return out;
}
  
  
/* Function:  GeneWise_from_ThreeStateModel(tsm,gp,cm,allN,gwcm)
 *
 * Descrip:    This makes a genewise model from a 
 *             threestatemodel for the genewise type
 *             algorithms.
 *
 *             Notice you have to provide the gene parameters
 *             being used
 *
 *             Stop is now not used
 *
 *
 * Arg:         tsm [UNKN ] Undocumented argument [ThreeStateModel *]
 * Arg:          gp [UNKN ] Undocumented argument [GeneParser21 *]
 * Arg:          cm [UNKN ] Undocumented argument [CodonMapper *]
 * Arg:        allN [UNKN ] Undocumented argument [Probability]
 * Arg:        gwcm [UNKN ] Undocumented argument [GeneWiseCodonModel *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWise *]
 *
 */
# line 306 "genewisemodel.dy"
GeneWise * GeneWise_from_ThreeStateModel(ThreeStateModel * tsm,GeneParser21 * gp,CodonMapper * cm,Probability allN,GeneWiseCodonModel * gwcm)
{
  register int i;
  GeneWise * out;
  Probability factor;

  assert(tsm);
  /*  assert(gp); */
  assert(cm);
  /*  assert(gwcm); */ /* can cope with null gwcm's */

  if( (out = GeneWise_alloc_len(tsm->len)) == NULL )
    return NULL;

  if( tsm->name != NULL )
    out->name = stringalloc(tsm->name);


  if( gp == NULL ) {
    factor = 1.0;
  } else {
   factor = (1.0 - removed_probability_from_cds(gp));
  }

  for(i=0;i<tsm->len;i++) {
    add_GeneWise(out,GeneWiseSegment_from_ThreeStateUnit(tsm->unit[i],factor,cm,gwcm,allN));
  }
  return out;
}

/* Function:  GeneWise_fold_in_synchronised_RandomModel(gw,rm,cm,*ct,stop_codon_background)
 *
 * Descrip:    This function places 'log-odd' scores of the
 *             genewise model assumming that the random model
 *             is a protein model with the codon mapper system
 *             added in, *and* that the path of the random model
 *             is synchronous with the query model.
 *
 *             It fudges stop codons with the stop score given
 *             as a probability.
 *
 *             In other words, this should give bits scores as
 *             if it was a protein, even though it is DNA
 *
 *
 * Arg:                           gw [UNKN ] Undocumented argument [GeneWise *]
 * Arg:                           rm [UNKN ] Undocumented argument [RandomModel *]
 * Arg:                           cm [UNKN ] Undocumented argument [CodonMapper *]
 * Arg:                          *ct [UNKN ] Undocumented argument [CodonTable]
 * Arg:        stop_codon_background [UNKN ] Undocumented argument [Probability]
 *
 */
# line 349 "genewisemodel.dy"
void GeneWise_fold_in_synchronised_RandomModel(GeneWise * gw,RandomModel * rm,CodonMapper * cm,CodonTable *ct,Probability stop_codon_background)
{
  int i;
  int j;
  double p;

  if( gw == NULL || rm == NULL || cm == NULL || ct == NULL ) {
    fatal("Null objects passed to GeneWise_fold_in_synchronised_RandomModel. Ugh!");
  }


  for(i=0;i<125;i++) {
    p = map_codon_CodonMapper(i,rm->aminoacid,cm);
    
    if( is_stop_codon(i,ct) == FALSE && p < 0.0000001 ) {
      warn("Got a close to zero probability for %d\n",i);
      p = 0.0000001;
    }

    for(j=0;j<gw->len;j++) {
      if( is_stop_codon(i,ct) == TRUE ){
	gw->seg[j]->match[i]  /= stop_codon_background;
	gw->seg[j]->insert[i] /= stop_codon_background;
      } else {
	gw->seg[j]->match[i]   = gw->seg[j]->match[i] / p;
	gw->seg[j]->insert[i]  = gw->seg[j]->insert[i] / p;
      }
    }
      
    
  }
  
  return;
}

/* Function:  check_flat_insert(gw,should_force,should_warn,ct)
 *
 * Descrip:    This function checks that the insert model is bang on
 *             zero, forcing it to zero
 *
 *             Potentially it warns for non zeros as well
 *
 *
 * Arg:                  gw [UNKN ] Undocumented argument [GeneWise *]
 * Arg:        should_force [UNKN ] Undocumented argument [boolean]
 * Arg:         should_warn [UNKN ] Undocumented argument [boolean]
 * Arg:                  ct [UNKN ] Undocumented argument [CodonTable *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 390 "genewisemodel.dy"
boolean check_flat_insert(GeneWise * gw,boolean should_force,boolean should_warn,CodonTable * ct)
{
  boolean ret = TRUE;
  int i,j;
  
  for(i=0;i<gw->len;i++) {
    for(j=0;j<125;j++) {
      if ( is_stop_codon(j,ct) ) {
	continue;
      }
      if( gw->seg[i]->insert[j] > 1.02 || gw->seg[i]->insert[j] < 0.99 ) {
	ret = FALSE;
	if( should_warn ) {
	  warn("In checking that we have a flat zero over the insert states, at %d, codon %d got %.2f",i,j,gw->seg[i]->insert[j]);
	}
	if( should_force ) {
	  gw->seg[i]->insert[j] = 1.0;
	}
      }
    }
  }

  return ret;
}

/* Function:  GeneWise_fold_in_RandomModelDNA(gw,rmd)
 *
 * Descrip:    This function folds in a simple random model
 *             (single base position) into a genewise model
 *
 *
 * Arg:         gw [UNKN ] Undocumented argument [GeneWise *]
 * Arg:        rmd [UNKN ] Undocumented argument [RandomModelDNA *]
 *
 */
# line 419 "genewisemodel.dy"
void GeneWise_fold_in_RandomModelDNA(GeneWise * gw,RandomModelDNA * rmd)
{
  register int i;
  register int j;
  Probability p;

  for(i=0;i<125;i++) {
    p = probability_of_this_codon(i,rmd);

    
    for(j=0;j<gw->len;j++) {
      gw->seg[j]->match[i]   = gw->seg[j]->match[i] / p;
      gw->seg[j]->insert[i]  = gw->seg[j]->insert[i] / p;
      
    }
    
  }
  
  return;
}

/* Function:  Protein_from_GeneWise_AlnColumn(dna,is_random_AlnColumn,col,position_in_aln,last_column,ct)
 *
 * Descrip:    Produces a protein object from a genewise/estwise
 *             style label set, setting the last retrieved column
 *
 *
 * Arg:                        dna [UNKN ] Undocumented argument [Sequence *]
 * Arg:        is_random_AlnColumn [UNKN ] Undocumented argument [NullString]
 * Arg:                        col [UNKN ] Undocumented argument [AlnColumn *]
 * Arg:            position_in_aln [UNKN ] Undocumented argument [int]
 * Arg:                last_column [UNKN ] Undocumented argument [AlnColumn **]
 * Arg:                         ct [UNKN ] Undocumented argument [CodonTable *]
 *
 * Return [UNKN ]  Undocumented return value [Protein *]
 *
 */
# line 444 "genewisemodel.dy"
Protein * Protein_from_GeneWise_AlnColumn(Sequence * dna,AlnColumn * col,int position_in_aln,AlnColumn ** last_column,CodonTable * ct,boolean (*is_random_AlnColumn)(const AlnColumn *))
{
  Sequence * out;
  char buffer[MAX_PROTEIN_GENEWISE]; /* max protein length? */
  char tempseq[4];
  int i =0;
  tempseq[3] = '\0';
  for(;col != NULL && (*is_random_AlnColumn)(col) == TRUE;col = col->next)
    ;
  if( col == NULL ) 
    return NULL;

  sprintf(buffer,"%s.pep",dna->name);
  out = empty_Sequence_from_dynamic_memory(stringalloc(buffer));
  for(;col != NULL && (*is_random_AlnColumn)(col) == FALSE;col = col->next) {
    if( i+1 >= MAX_PROTEIN_GENEWISE ) {
      buffer[i] = '\0';
      add_string_to_Sequence(out,buffer);
      i = 0;
    }

    if( strstr(col->alu[position_in_aln]->text_label,"CODON") != NULL ) {
      buffer[i++] = aminoacid_from_seq(ct,dna->seq+col->alu[position_in_aln]->start+1);
    } else if ( strstr(col->alu[position_in_aln]->text_label,"5SS_PHASE_1") != NULL ) {
      tempseq[0] = dna->seq[col->alu[position_in_aln]->start+1];
      for(col=col->next;col != NULL && strstr(col->alu[position_in_aln]->text_label,"3SS") == NULL;col=col->next)
	;
      if( col == NULL ) {
	warn("In middle of intron - got no 3'SS in making peptide translation");
	return NULL;
      }
      tempseq[1] = dna->seq[col->alu[position_in_aln]->start+4];
      tempseq[2] = dna->seq[col->alu[position_in_aln]->start+5];
      /*fprintf(stderr,"In phase 1 intron, calling %c%c%c as split codon and %c%c%c as last codon\n",tempseq[0],tempseq[1],tempseq[2],
	      dna->seq[col->alu[position_in_aln]->start+6],
	      dna->seq[col->alu[position_in_aln]->start+7],
	      dna->seq[col->alu[position_in_aln]->start+8]
	      ); */
      buffer[i++] = aminoacid_from_seq(ct,tempseq);
      /* buffer[i++] = aminoacid_from_seq(ct,dna->seq+col->alu[position_in_aln]->start+6);*/
      printf("Making a %c in phase 1 intron\n",buffer[i-1]);
    } else if ( strstr(col->alu[position_in_aln]->text_label,"5SS_PHASE_2") != NULL ) {
      tempseq[0] = dna->seq[col->alu[position_in_aln]->start+1];
      tempseq[1] = dna->seq[col->alu[position_in_aln]->start+2];
      for(col=col->next;col != NULL && strstr(col->alu[position_in_aln]->text_label,"3SS") == NULL;col=col->next)
	;
      if( col == NULL ) {
	warn("In middle of intron - got no 3'SS in making peptide translation");
	return NULL;
      }
      tempseq[2] = dna->seq[col->alu[position_in_aln]->start+4];
      buffer[i++] = aminoacid_from_seq(ct,tempseq);
      /*buffer[i++] = aminoacid_from_seq(ct,dna->seq+col->alu[position_in_aln]->start+5);*/
      printf("Making a %c in phase 2 intron\n",buffer[i-1]);
    } else if ( strstr(col->alu[position_in_aln]->text_label,"5SS_PHASE_0") != NULL ) {
      /* codon already delt with! */
      for(col=col->next;col != NULL && strstr(col->alu[position_in_aln]->text_label,"3SS") == NULL;col=col->next)
	;
      if( col == NULL ) {
	warn("In middle of intron - got no 3'SS in making peptide translation");
	return NULL;
      }
      /* buffer already sorted out. No need to provide compute */
      continue;
      buffer[i++] = aminoacid_from_seq(ct,dna->seq+col->alu[position_in_aln]->start+3);
      printf("Making a %c in phase 0 intron\n",buffer[i-1]);

      col = col->next;
    } else if ( strstr(col->alu[position_in_aln]->text_label,"SEQUENCE_DELETION") != NULL ) {
      buffer[i++] = 'X';
    } else if ( strstr(col->alu[position_in_aln]->text_label,"SEQUENCE_INSERTION") != NULL ) {
      buffer[i++] = 'X';
    } else if ( strstr(col->alu[position_in_aln]->text_label,"INSERT") != NULL ) {
      continue;
    } else {
      warn("In processing alignment to peptide, got label %s which cannot handle. Assumming X in protein translation",col->alu[position_in_aln]->text_label);
      buffer[i++] = 'X';
    }
  }

  if( last_column != NULL ) 
    *last_column = col;
  buffer[i] = '\0';

  add_string_to_Sequence(out,buffer);
  out->type = SEQUENCE_PROTEIN; /* force to protein */
  return Protein_from_Sequence(out);
}

  
/* Function:  probability_of_this_codon(codon,rmd)
 *
 * Descrip:    Helper function for getting probability of
 *             codon for a random model
 *
 *
 * Arg:        codon [UNKN ] Undocumented argument [int]
 * Arg:          rmd [UNKN ] Undocumented argument [RandomModelDNA *]
 *
 * Return [UNKN ]  Undocumented return value [Probability]
 *
 */
# line 539 "genewisemodel.dy"
Probability probability_of_this_codon(int codon,RandomModelDNA * rmd)
{
  base one;
  base two;
  base three;

  all_bases_from_codon(codon,&one,&two,&three);

  return rmd->base[one] * rmd->base[two] * rmd->base[three];
}
  
/* Function:  show_GeneWise(gw,ofp)
 *
 * Descrip:    debugging function
 *
 *
 * Arg:         gw [UNKN ] Undocumented argument [GeneWise *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 554 "genewisemodel.dy"
void show_GeneWise(GeneWise * gw,FILE * ofp)
{
  register int i;


  fprintf(stderr,"Got here at least [%d]\n",gw->len);

  for(i=0;i<gw->len;i++)
    show_GeneWiseSegment(gw->seg[i],ofp);

}

/* Function:  show_GeneWiseSegment(seg,ofp)
 *
 * Descrip:    debugging
 *
 *
 * Arg:        seg [UNKN ] Undocumented argument [GeneWiseSegment *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 570 "genewisemodel.dy"
void show_GeneWiseSegment(GeneWiseSegment * seg,FILE * ofp)
{
  fprintf(ofp,"match=\" ");
  show_Probability_array(seg->match,125,ofp); 
  fprintf(ofp,"\n");
  show_Probability_array(seg->insert,125,ofp); 
  fprintf(ofp,"\n");
  
}

/* Function:  GeneWiseSegment_from_ThreeStateUnit(tsu,factor,cm,gwcm,allN)
 *
 * Descrip:    Function which actually does the mapping from
 *             threestate model unit to genewise
 *
 *
 * Arg:           tsu [UNKN ] Undocumented argument [ThreeStateUnit *]
 * Arg:        factor [UNKN ] Undocumented argument [Probability]
 * Arg:            cm [UNKN ] Undocumented argument [CodonMapper *]
 * Arg:          gwcm [UNKN ] Undocumented argument [GeneWiseCodonModel *]
 * Arg:          allN [UNKN ] Undocumented argument [Probability]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseSegment *]
 *
 */
# line 585 "genewisemodel.dy"
GeneWiseSegment * GeneWiseSegment_from_ThreeStateUnit(ThreeStateUnit * tsu,Probability factor,CodonMapper * cm,GeneWiseCodonModel * gwcm,Probability allN)
{
  register int i;
  GeneWiseSegment * out;
  double total;
  int codon;

  
  out = GeneWiseSegment_alloc();
  if( out == NULL )
    return NULL;

  (void)set_Probability_array(out->match     , 0.0, GW_EMISSION_LEN);
  (void)set_Probability_array(out->insert    , 0.0, GW_EMISSION_LEN);
  (void)set_Probability_array(out->transition, 0.0, GW_TRANSITION_LEN);
  
  true_map_codon_array_CodonMapper(out->match,tsu->match_emission,cm);
  true_map_codon_array_CodonMapper(out->insert,tsu->insert_emission,cm);

  codon = codon_from_seq("NNN");
  out->match[codon] = allN;
  out->insert[codon] = allN;
  


  out->transition[GW_MATCH2MATCH]   = tsu->transition[TSM_MATCH2MATCH] * factor;
  out->transition[GW_MATCH2INSERT]  = tsu->transition[TSM_MATCH2INSERT] * factor;
  out->transition[GW_MATCH2DELETE]  = tsu->transition[TSM_MATCH2DELETE] * factor;
  out->transition[GW_MATCH2END]     = tsu->transition[TSM_MATCH2END] * factor;

  out->transition[GW_INSERT2MATCH]  = tsu->transition[TSM_INSERT2MATCH] * factor;
  out->transition[GW_INSERT2INSERT] = tsu->transition[TSM_INSERT2INSERT] * factor;
  out->transition[GW_INSERT2DELETE] = tsu->transition[TSM_INSERT2DELETE] * factor;
  out->transition[GW_INSERT2END]    = tsu->transition[TSM_INSERT2END] * factor;

  out->transition[GW_DELETE2MATCH]  = tsu->transition[TSM_DELETE2MATCH];
  out->transition[GW_DELETE2INSERT] = tsu->transition[TSM_DELETE2INSERT];
  out->transition[GW_DELETE2DELETE] = tsu->transition[TSM_DELETE2DELETE];
  out->transition[GW_DELETE2END]    = tsu->transition[TSM_DELETE2END];
  
  out->transition[GW_START2MATCH]   = tsu->transition[TSM_START2MATCH];
  out->transition[GW_START2INSERT]  = tsu->transition[TSM_START2INSERT];
  out->transition[GW_START2DELETE]  = tsu->transition[TSM_START2DELETE];

  /** we need 1/(sum_over_codons (match(codon)) * no(codon in 5'SS)/no(codon in cds) ) **/

  if( gwcm != NULL ) {

    total = 0.0;
    for(i=0;i<64;i++) {
      codon = codon_from_base4_codon(i);
      total += out->match[codon] * gwcm->in_donor[i]/gwcm->in_cds[i];
    }
    out->transition[GW_MATCH_BALANCE_5SS] = 1.0 / total;

    total = 0.0;
    for(i=0;i<64;i++) {
      codon = codon_from_base4_codon(i);
      total += out->insert[codon] * gwcm->in_donor[i]/gwcm->in_cds[i];
    }

    out->transition[GW_INSERT_BALANCE_5SS] = 1.0 / total;
   
    total = 0.0;
    for(i=0;i<64;i++) {
      codon = codon_from_base4_codon(i);
      total += out->match[codon] * gwcm->in_acceptor[i]/gwcm->in_cds[i];
    }
    out->transition[GW_MATCH_BALANCE_3SS] = 1.0 / total;

    total = 0.0;
    for(i=0;i<64;i++) {
      codon = codon_from_base4_codon(i);
      total += out->insert[codon] * gwcm->in_acceptor[i]/gwcm->in_cds[i];
    }
    out->transition[GW_INSERT_BALANCE_3SS] = 1.0 / total;
  } else {
    out->transition[GW_MATCH_BALANCE_5SS]  = 1.0;
    out->transition[GW_MATCH_BALANCE_3SS]  = 1.0;
    out->transition[GW_INSERT_BALANCE_5SS] = 1.0;
    out->transition[GW_INSERT_BALANCE_3SS] = 1.0;
  }

  return out;
}

/* Function:  Probability_of_codon(codon,ct,aminoacid_26_array,stop)
 *
 * Descrip:    Helper function
 *
 *
 * Arg:                     codon [UNKN ] Undocumented argument [int]
 * Arg:                        ct [UNKN ] Undocumented argument [CodonTable *]
 * Arg:        aminoacid_26_array [UNKN ] Undocumented argument [Probability *]
 * Arg:                      stop [UNKN ] Undocumented argument [Probability]
 *
 * Return [UNKN ]  Undocumented return value [Probability]
 *
 */
# line 675 "genewisemodel.dy"
Probability Probability_of_codon(int codon,CodonTable * ct,Probability * aminoacid_26_array,Probability stop)
{
  register int i,j,k;
  base base1,base2,base3;
  Probability ret = 0.0;


  if( has_random_bases(codon) == FALSE) {
    if( is_stop_codon(codon,ct) == TRUE ) {
      return stop;
    }
    else {
      fprintf(stderr,"Setting a number to %f\n",aminoacid_26_array[aminoacid_no_from_codon(ct,codon)]);
      return aminoacid_26_array[aminoacid_no_from_codon(ct,codon)];
    }
  }

  all_bases_from_codon(codon,&base1,&base2,&base3);


  /*** the first base is N ***/

  if( base1 == BASE_N && base2 != BASE_N && base3 != BASE_N) {
    for(i=0;i<4;i++ ) 
      ret += aminoacid_26_array[aminoacid_no_from_codon(ct,(i*25 + base2*5 + base3))];
    return ret;
  }

  if( base1 == BASE_N && base2 == BASE_N && base3 != BASE_N) {
    for(i=0;i<4;i++) 
      for(j=0;j<4;j++)
	ret += aminoacid_26_array[aminoacid_no_from_codon(ct,(i*25 + j*5 + base3))];
    return ret;
  }

  if( base1 == BASE_N && base2 != BASE_N && base3 == BASE_N) {
    for(i=0;i<4;i++) 
      for(k=0;k<4;k++) 
	ret += aminoacid_26_array[aminoacid_no_from_codon(ct,(i*25 + base2*5 + k))];
    return ret;
  }

  if( base1 == BASE_N && base2 == BASE_N && base3 == BASE_N) {
    for(i=0;i<4;i++) 
      for(j=0;j<4;j++) 
	for(k=0;k<4;k++)
	  ret += aminoacid_26_array[aminoacid_no_from_codon(ct,(i*25 + j*5 + k))];
    return ret;
  }


  /*** the second base is N ***/

  if( base1 != BASE_N && base2 == BASE_N && base3 != BASE_N) {
    for(j=0;j<4;j++) 
      ret += aminoacid_26_array[aminoacid_no_from_codon(ct,(base1*25 + j*5 + base3))];
    return ret;
  }


  if( base1 != BASE_N && base2 == BASE_N && base3 == BASE_N) {
    for(j=0;j<4;j++) 
      for(k=0;k<4;k++)
	ret += aminoacid_26_array[aminoacid_no_from_codon(ct,(base1*25 + j*5 + k))];
    return ret;
  }

  /*** the third base is N ***/


  if( base1 != BASE_N && base2 != BASE_N && base3 == BASE_N) {
    for(k=0;k<4;k++)
      ret += aminoacid_26_array[aminoacid_no_from_codon(ct,(base1*25 + base2*5 + k))];
    return ret;
  }

  /*** should never reach here ***/

  warn("Got to the end of Probability_of_codon without a BASE_N combination being triggered. This looks like a major problem. Codon inputted was %d",codon);

  return 0.0;
  
}

/* Function:  GeneWiseScore_from_GeneWise(gw)
 *
 * Descrip:    Makes a Score (log based) object from
 *             a probability based object
 *
 *
 * Arg:        gw [UNKN ] Undocumented argument [GeneWise *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseScore *]
 *
 */
# line 763 "genewisemodel.dy"
GeneWiseScore * GeneWiseScore_from_GeneWise(GeneWise * gw)
{
  GeneWiseScore * out;
  register int i;
  

  out = GeneWiseScore_alloc_len(gw->len);

  if(gw->name != NULL )
    out->name = stringalloc(gw->name);

  add_GeneWiseScore(out,GeneWiseScoreSegment_from_GeneWiseSegment(NULL,gw->seg[0]));

  for(i=1;i<gw->len;i++) 
    add_GeneWiseScore(out,GeneWiseScoreSegment_from_GeneWiseSegment(gw->seg[i-1],gw->seg[i]));

  return out;
}

/* Function:  GeneWiseScoreSegment_from_GeneWiseSegment(prev,seg)
 *
 * Descrip:    helper function for prob to score mapping
 *
 *
 * Arg:        prev [UNKN ] Undocumented argument [GeneWiseSegment *]
 * Arg:         seg [UNKN ] Undocumented argument [GeneWiseSegment *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseScoreSegment *]
 *
 */
# line 786 "genewisemodel.dy"
GeneWiseScoreSegment * GeneWiseScoreSegment_from_GeneWiseSegment(GeneWiseSegment * prev,GeneWiseSegment * seg)
{
  GeneWiseScoreSegment * out;

  out = GeneWiseScoreSegment_alloc();

  Probability2Score_move(seg->match,out->match,GW_EMISSION_LEN);
  Probability2Score_move(seg->insert,out->insert,GW_EMISSION_LEN);

  if( prev != NULL ) {
    out->transition[GW_MATCH2MATCH] = Probability2Score(prev->transition[GW_MATCH2MATCH]);
    out->transition[GW_INSERT2MATCH] = Probability2Score(prev->transition[GW_INSERT2MATCH]);
    out->transition[GW_DELETE2MATCH] = Probability2Score(prev->transition[GW_DELETE2MATCH]);

    out->transition[GW_MATCH2DELETE] = Probability2Score(prev->transition[GW_MATCH2DELETE]);
    out->transition[GW_INSERT2DELETE] = Probability2Score(prev->transition[GW_INSERT2DELETE]);
    out->transition[GW_DELETE2DELETE] = Probability2Score(prev->transition[GW_DELETE2DELETE]);
  } else {

    out->transition[GW_MATCH2MATCH] = NEGI;
    out->transition[GW_INSERT2MATCH] = NEGI;
    out->transition[GW_DELETE2MATCH] = NEGI;

    out->transition[GW_MATCH2DELETE] = NEGI;
    out->transition[GW_INSERT2DELETE] = NEGI;
    out->transition[GW_DELETE2DELETE] = NEGI;
  }

  out->transition[GW_MATCH2INSERT] = Probability2Score(seg->transition[GW_MATCH2INSERT]);
  out->transition[GW_INSERT2INSERT] = Probability2Score(seg->transition[GW_INSERT2INSERT]);
  out->transition[GW_DELETE2INSERT] = Probability2Score(seg->transition[GW_DELETE2INSERT]);
  
  out->transition[GW_START2MATCH] = Probability2Score(seg->transition[GW_START2MATCH]);
  out->transition[GW_START2INSERT] = Probability2Score(seg->transition[GW_START2INSERT]);
  out->transition[GW_START2DELETE] = Probability2Score(seg->transition[GW_START2DELETE]);
  
  out->transition[GW_MATCH2END] = Probability2Score(seg->transition[GW_MATCH2END]);
  out->transition[GW_INSERT2END] = Probability2Score(seg->transition[GW_INSERT2END]);
  out->transition[GW_DELETE2END] = Probability2Score(seg->transition[GW_DELETE2END]);

  out->transition[GW_MATCH_BALANCE_5SS] = Probability2Score(seg->transition[GW_MATCH_BALANCE_5SS]);
  out->transition[GW_INSERT_BALANCE_5SS] = Probability2Score(seg->transition[GW_INSERT_BALANCE_5SS]);

  out->transition[GW_MATCH_BALANCE_3SS] = Probability2Score(seg->transition[GW_MATCH_BALANCE_3SS]);
  out->transition[GW_INSERT_BALANCE_3SS] = Probability2Score(seg->transition[GW_INSERT_BALANCE_3SS]);
  
  return out;
}





# line 867 "genewisemodel.c"
/* Function:  hard_link_GeneWiseSegment(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GeneWiseSegment *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseSegment *]
 *
 */
GeneWiseSegment * hard_link_GeneWiseSegment(GeneWiseSegment * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a GeneWiseSegment object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  GeneWiseSegment_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseSegment *]
 *
 */
GeneWiseSegment * GeneWiseSegment_alloc(void) 
{
    GeneWiseSegment * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GeneWiseSegment *) ckalloc (sizeof(GeneWiseSegment))) == NULL)  {  
      warn("GeneWiseSegment_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    /* match[GW_EMISSION_LEN] is an array: no default possible */ 
    /* insert[GW_EMISSION_LEN] is an array: no default possible */ 
    /* transition[GW_TRANSITION_LEN] is an array: no default possible */ 


    return out;  
}    


/* Function:  free_GeneWiseSegment(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneWiseSegment *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseSegment *]
 *
 */
GeneWiseSegment * free_GeneWiseSegment(GeneWiseSegment * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a GeneWiseSegment obj. Should be trappable");   
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


/* Function:  swap_GeneWise(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_GeneWise
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [GeneWiseSegment **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_GeneWise(GeneWiseSegment ** list,int i,int j)  
{
    GeneWiseSegment * temp;  
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_GeneWise(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_GeneWise which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [GeneWiseSegment **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_GeneWise(GeneWiseSegment ** list,int left,int right,int (*comp)(GeneWiseSegment * ,GeneWiseSegment * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_GeneWise(list,left,(left+right)/2); 
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_GeneWise (list,++last,i);   
      }  
    swap_GeneWise (list,left,last);  
    qsort_GeneWise(list,left,last-1,comp);   
    qsort_GeneWise(list,last+1,right,comp);  
}    


/* Function:  sort_GeneWise(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_GeneWise
 *
 *
 * Arg:         obj [UNKN ] Object containing list [GeneWise *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_GeneWise(GeneWise * obj,int (*comp)(GeneWiseSegment *, GeneWiseSegment *)) 
{
    qsort_GeneWise(obj->seg,0,obj->len-1,comp);  
    return;  
}    


/* Function:  expand_GeneWise(obj,len)
 *
 * Descrip:    Really an internal function for add_GeneWise
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GeneWise *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_GeneWise(GeneWise * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_GeneWise called with no need");   
      return TRUE;   
      }  


    if( (obj->seg = (GeneWiseSegment ** ) ckrealloc (obj->seg,sizeof(GeneWiseSegment *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_GeneWise, returning FALSE"); 
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_GeneWise(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GeneWise *]
 * Arg:        add [OWNER] Object to add to the list [GeneWiseSegment *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_GeneWise(GeneWise * obj,GeneWiseSegment * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_GeneWise(obj,obj->len + GeneWiseLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->seg[obj->len++]=add;    
    return TRUE; 
}    


/* Function:  flush_GeneWise(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [GeneWise *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_GeneWise(GeneWise * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->seg[i] != NULL)   {  
        free_GeneWiseSegment(obj->seg[i]);   
        obj->seg[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  GeneWise_alloc_std(void)
 *
 * Descrip:    Equivalent to GeneWise_alloc_len(GeneWiseLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneWise *]
 *
 */
GeneWise * GeneWise_alloc_std(void) 
{
    return GeneWise_alloc_len(GeneWiseLISTLENGTH);   
}    


/* Function:  GeneWise_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [GeneWise *]
 *
 */
GeneWise * GeneWise_alloc_len(int len) 
{
    GeneWise * out; /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = GeneWise_alloc()) == NULL) 
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->seg = (GeneWiseSegment ** ) ckcalloc (len,sizeof(GeneWiseSegment *))) == NULL)  {  
      warn("Warning, ckcalloc failed in GeneWise_alloc_len");    
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_GeneWise(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GeneWise *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWise *]
 *
 */
GeneWise * hard_link_GeneWise(GeneWise * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a GeneWise object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  GeneWise_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneWise *]
 *
 */
GeneWise * GeneWise_alloc(void) 
{
    GeneWise * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GeneWise *) ckalloc (sizeof(GeneWise))) == NULL)    {  
      warn("GeneWise_alloc failed ");    
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


/* Function:  free_GeneWise(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneWise *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWise *]
 *
 */
GeneWise * free_GeneWise(GeneWise * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a GeneWise obj. Should be trappable");  
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
          free_GeneWiseSegment(obj->seg[i]); 
        }  
      ckfree(obj->seg);  
      }  
    if( obj->name != NULL)   
      ckfree(obj->name);     


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_GeneWiseScoreSegment(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GeneWiseScoreSegment *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseScoreSegment *]
 *
 */
GeneWiseScoreSegment * hard_link_GeneWiseScoreSegment(GeneWiseScoreSegment * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a GeneWiseScoreSegment object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  GeneWiseScoreSegment_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseScoreSegment *]
 *
 */
GeneWiseScoreSegment * GeneWiseScoreSegment_alloc(void) 
{
    GeneWiseScoreSegment * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GeneWiseScoreSegment *) ckalloc (sizeof(GeneWiseScoreSegment))) == NULL)    {  
      warn("GeneWiseScoreSegment_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    /* match[GW_EMISSION_LEN] is an array: no default possible */ 
    /* insert[GW_EMISSION_LEN] is an array: no default possible */ 
    /* transition[GW_TRANSITION_LEN] is an array: no default possible */ 


    return out;  
}    


/* Function:  free_GeneWiseScoreSegment(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneWiseScoreSegment *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseScoreSegment *]
 *
 */
GeneWiseScoreSegment * free_GeneWiseScoreSegment(GeneWiseScoreSegment * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a GeneWiseScoreSegment obj. Should be trappable");  
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


/* Function:  swap_GeneWiseScore(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_GeneWiseScore
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [GeneWiseScoreSegment **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_GeneWiseScore(GeneWiseScoreSegment ** list,int i,int j)  
{
    GeneWiseScoreSegment * temp; 
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_GeneWiseScore(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_GeneWiseScore which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [GeneWiseScoreSegment **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_GeneWiseScore(GeneWiseScoreSegment ** list,int left,int right,int (*comp)(GeneWiseScoreSegment * ,GeneWiseScoreSegment * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_GeneWiseScore(list,left,(left+right)/2);    
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_GeneWiseScore (list,++last,i);  
      }  
    swap_GeneWiseScore (list,left,last); 
    qsort_GeneWiseScore(list,left,last-1,comp);  
    qsort_GeneWiseScore(list,last+1,right,comp); 
}    


/* Function:  sort_GeneWiseScore(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_GeneWiseScore
 *
 *
 * Arg:         obj [UNKN ] Object containing list [GeneWiseScore *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_GeneWiseScore(GeneWiseScore * obj,int (*comp)(GeneWiseScoreSegment *, GeneWiseScoreSegment *)) 
{
    qsort_GeneWiseScore(obj->seg,0,obj->len-1,comp); 
    return;  
}    


/* Function:  expand_GeneWiseScore(obj,len)
 *
 * Descrip:    Really an internal function for add_GeneWiseScore
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GeneWiseScore *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_GeneWiseScore(GeneWiseScore * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_GeneWiseScore called with no need");  
      return TRUE;   
      }  


    if( (obj->seg = (GeneWiseScoreSegment ** ) ckrealloc (obj->seg,sizeof(GeneWiseScoreSegment *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_GeneWiseScore, returning FALSE");    
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_GeneWiseScore(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GeneWiseScore *]
 * Arg:        add [OWNER] Object to add to the list [GeneWiseScoreSegment *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_GeneWiseScore(GeneWiseScore * obj,GeneWiseScoreSegment * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_GeneWiseScore(obj,obj->len + GeneWiseScoreLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->seg[obj->len++]=add;    
    return TRUE; 
}    


/* Function:  flush_GeneWiseScore(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [GeneWiseScore *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_GeneWiseScore(GeneWiseScore * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->seg[i] != NULL)   {  
        free_GeneWiseScoreSegment(obj->seg[i]);  
        obj->seg[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  GeneWiseScore_alloc_std(void)
 *
 * Descrip:    Equivalent to GeneWiseScore_alloc_len(GeneWiseScoreLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseScore *]
 *
 */
GeneWiseScore * GeneWiseScore_alloc_std(void) 
{
    return GeneWiseScore_alloc_len(GeneWiseScoreLISTLENGTH); 
}    


/* Function:  GeneWiseScore_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseScore *]
 *
 */
GeneWiseScore * GeneWiseScore_alloc_len(int len) 
{
    GeneWiseScore * out;/* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = GeneWiseScore_alloc()) == NULL)    
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->seg = (GeneWiseScoreSegment ** ) ckcalloc (len,sizeof(GeneWiseScoreSegment *))) == NULL)    {  
      warn("Warning, ckcalloc failed in GeneWiseScore_alloc_len");   
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_GeneWiseScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GeneWiseScore *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseScore *]
 *
 */
GeneWiseScore * hard_link_GeneWiseScore(GeneWiseScore * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a GeneWiseScore object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  GeneWiseScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseScore *]
 *
 */
GeneWiseScore * GeneWiseScore_alloc(void) 
{
    GeneWiseScore * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GeneWiseScore *) ckalloc (sizeof(GeneWiseScore))) == NULL)  {  
      warn("GeneWiseScore_alloc failed ");   
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


/* Function:  free_GeneWiseScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneWiseScore *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseScore *]
 *
 */
GeneWiseScore * free_GeneWiseScore(GeneWiseScore * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a GeneWiseScore obj. Should be trappable"); 
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
          free_GeneWiseScoreSegment(obj->seg[i]);    
        }  
      ckfree(obj->seg);  
      }  
    if( obj->name != NULL)   
      ckfree(obj->name);     


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_GeneWiseScoreFlat(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GeneWiseScoreFlat *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseScoreFlat *]
 *
 */
GeneWiseScoreFlat * hard_link_GeneWiseScoreFlat(GeneWiseScoreFlat * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a GeneWiseScoreFlat object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  GeneWiseScoreFlat_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseScoreFlat *]
 *
 */
GeneWiseScoreFlat * GeneWiseScoreFlat_alloc(void) 
{
    GeneWiseScoreFlat * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GeneWiseScoreFlat *) ckalloc (sizeof(GeneWiseScoreFlat))) == NULL)  {  
      warn("GeneWiseScoreFlat_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->seg = NULL; 
    out->len = 0;    


    return out;  
}    



#ifdef _cplusplus
}
#endif
