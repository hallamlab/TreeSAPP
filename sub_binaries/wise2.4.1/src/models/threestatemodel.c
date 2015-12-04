#ifdef _cplusplus
extern "C" {
#endif
#include "threestatemodel.h"

 char * std_alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";


/* Function:  information_from_ThreeStateUnit(tsu,rm)
 *
 * Descrip:     Gets the information content (K-L divergence) vs a background for a position
 *
 *
 * Arg:        tsu [UNKN ] Undocumented argument [ThreeStateUnit *]
 * Arg:         rm [UNKN ] Undocumented argument [RandomModel *]
 *
 * Return [UNKN ]  Undocumented return value [double]
 *
 */
# line 117 "threestatemodel.dy"
double information_from_ThreeStateUnit(ThreeStateUnit * tsu,RandomModel * rm)
{
  int i;
  double info = 0.0;
  double p;

  for(i=0;i<ALPHABET_SIZE;i++) {
    if( rm->aminoacid[i] < 1.0 ) {
      p = rm->aminoacid[i]*Probability2Bits(rm->aminoacid[i]/tsu->match_emission[i]);
      info += p;
      /*  fprintf(stderr,"For amino acid %d, adding %f %f\n",i,p,rm->aminoacid[i]); */
    }
  }
  
  return info;
}

/* Function:  threestatemodel_mode_from_string(mode)
 *
 * Descrip:    gets out the mode from a string
 *
 *
 * Arg:        mode [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 137 "threestatemodel.dy"
int threestatemodel_mode_from_string(char * mode)
{
  int t;

  t = get_number_from_slashed_string(mode,"default/global/local/wing/endbias");

  switch (t) {
  case 0 : return TSM_default;
  case 1 : return TSM_global;
  case 2 : return TSM_local;
  case 3 : return TSM_wing;
  case 4 : return TSM_endbiased;
  default : warn("%s is not a valid TSM mode!",mode); return TSM_unknown;
  }

  return TSM_unknown;
}

/* Function:  set_startend_policy_ThreeStateModel(tsm,mode,wing_length,internal_bias)
 *
 * Descrip:    Sets the start/end policy on the basis of the mode
 *
 *
 * Arg:                  tsm [UNKN ] Undocumented argument [ThreeStateModel *]
 * Arg:                 mode [UNKN ] Undocumented argument [TSM_StartEndMode]
 * Arg:          wing_length [UNKN ] Undocumented argument [int]
 * Arg:        internal_bias [UNKN ] Undocumented argument [Probability]
 *
 */
# line 158 "threestatemodel.dy"
void set_startend_policy_ThreeStateModel(ThreeStateModel * tsm,TSM_StartEndMode mode,int wing_length,Probability internal_bias)
{

  switch(mode) {
  case TSM_default   : return;
  case TSM_global    : force_weighted_local_model(tsm,1.0,1.0,1.0); return;
  case TSM_local     : force_hmmfs_ThreeStateModel(tsm); return;
  case TSM_wing      : force_wing_local_model(tsm,0.75,wing_length); return;
  case TSM_endbiased : force_endbias_model(tsm,1.0,internal_bias); return; 
  case TSM_unknown   :
  default :
    warn("No valid mode passed into set_startend_policy");
  }

}

/* Function:  force_endbias_model(tsm,startend,internal)
 *
 * Descrip:    Makes start/end on probability and
 *             internal another
 *
 *             not probabilistically correct
 *
 *
 * Arg:             tsm [UNKN ] Undocumented argument [ThreeStateModel *]
 * Arg:        startend [UNKN ] Undocumented argument [double]
 * Arg:        internal [UNKN ] Undocumented argument [double]
 *
 */
# line 180 "threestatemodel.dy"
void force_endbias_model(ThreeStateModel * tsm,double startend,double internal)
{
  int i;

  if( tsm == NULL ) {
    warn("In force_endbias_model, got a NULL threestate model! Problem!");
    return;
  }

  for(i=0;i<tsm->len;i++) {
    tsm->unit[i]->transition[TSM_START2MATCH]  = internal;
    tsm->unit[i]->transition[TSM_START2DELETE] = internal;
    tsm->unit[i]->transition[TSM_START2INSERT] = internal;
    tsm->unit[i]->transition[TSM_MATCH2END]    = internal;
    tsm->unit[i]->transition[TSM_INSERT2END]   = internal;
    tsm->unit[i]->transition[TSM_DELETE2END]   = internal;
  }

  tsm->unit[0]->transition[TSM_START2MATCH] = startend;
  tsm->unit[0]->transition[TSM_START2DELETE] = startend;

  tsm->unit[tsm->len-1]->transition[TSM_MATCH2END] = startend;
  tsm->unit[tsm->len-1]->transition[TSM_DELETE2END] = startend;
} 


/* Function:  force_global_model(tsm,prob_into_model)
 *
 * Descrip:    Makes start at position 0 and end at position end,
 *             no other positions being valid
 *
 *
 *
 * Arg:                    tsm [UNKN ] ThreeStateModel to be 'forced' [ThreeStateModel *]
 * Arg:        prob_into_model [UNKN ] Probability to start the model: for true global will be 1.0 [double]
 *
 */
# line 214 "threestatemodel.dy"
void force_global_model(ThreeStateModel * tsm,double prob_into_model) 
{
  force_weighted_local_model(tsm,prob_into_model,1.0,1.0);
}

/* Function:  force_wing_local_model(tsm,terminus_prob,wing_length)
 *
 * Descrip:    Places the balanace of the start probability
 *             at the start and then the rest spread evenly
 *             descending over the wing stretch of sequences
 *
 *             Same with the end probability
 *
 *
 * Arg:                  tsm [UNKN ] Undocumented argument [ThreeStateModel *]
 * Arg:        terminus_prob [UNKN ] the amount of probability to put on the real start and end [double]
 * Arg:          wing_length [UNKN ] the rest of the probability spread over this distance into the wing [int]
 *
 */
# line 229 "threestatemodel.dy"
void force_wing_local_model(ThreeStateModel * tsm,double terminus_prob,int wing_length)
{
  int i;
  double k;
  

  /*
   * k is the amount to take off for each step along the wing
   * 
   * it comes from solving sum_(terminus_prob - nk) = 1.0 
   * with n going from 0 to wing_length
   *
   * I make this k = 2(terminus_prob)/wing_length - 2/(wing_length+1)wing_length
   *
   */

  k = (2.0 * terminus_prob / wing_length) - 2.0 / ((wing_length + 1)*wing_length);

  if( k < 0.0 ) {
    warn("Weird - got k less than zero in force wing model");
    return;
  }


  if( tsm == NULL ) {
    warn("In force_wing_local_model, got a NULL threestate model! Problem!");
    return;
  }

  
  

  for(i=0;i<tsm->len && i < wing_length ;i++) {
    tsm->unit[i]->transition[TSM_START2MATCH] = (terminus_prob) - (i * k) ;
    tsm->unit[i]->transition[TSM_START2DELETE] = (terminus_prob) - (i * k) ;
    tsm->unit[i]->transition[TSM_START2INSERT] = (terminus_prob) - ( i * k) ;

    tsm->unit[tsm->len-1-i]->transition[TSM_MATCH2END] = (terminus_prob) - (i * k) ;
    tsm->unit[tsm->len-1-i]->transition[TSM_INSERT2END] = (terminus_prob) - ( i * k) ;
    tsm->unit[tsm->len-1-i]->transition[TSM_DELETE2END] = (terminus_prob) - (i * k)  ;

  }
} 


/* Function:  force_weighted_local_model(tsm,prob_into_model,ratio_start,ratio_end)
 *
 * Descrip:    places the ratio of probability to start/end,
 *             and then distributes the rest over the start/end
 *
 *
 *
 * Arg:                    tsm [UNKN ] ThreeStateModel to be 'forced' [ThreeStateModel *]
 * Arg:        prob_into_model [UNKN ] Probability to start the model: for true global will be 1.0 [double]
 * Arg:            ratio_start [UNKN ] ratio of prob to unit 0 to the rest (1.0 means all goes to start) [double]
 * Arg:              ratio_end [UNKN ] ratio of prob to unit (last) to the rest (1.0 means all goes to the end) [double]
 *
 */
# line 284 "threestatemodel.dy"
void force_weighted_local_model(ThreeStateModel * tsm,double prob_into_model,double ratio_start,double ratio_end) 
{
  int i;

  if( tsm == NULL ) {
    warn("In force_weighted_local_model, got a NULL threestate model! Problem!");
    return;
  }

  for(i=0;i<tsm->len;i++) {
    tsm->unit[i]->transition[TSM_START2MATCH] = prob_into_model * (1.0 - ratio_start) / (tsm->len-1) ;
    tsm->unit[i]->transition[TSM_START2DELETE] = prob_into_model * (1.0 - ratio_start) / (tsm->len-1);
    tsm->unit[i]->transition[TSM_START2INSERT] = prob_into_model * (1.0 - ratio_start)/ (tsm->len-1) ;
    tsm->unit[i]->transition[TSM_MATCH2END] =  (1.0 - ratio_end) / (tsm->len-1);
    tsm->unit[i]->transition[TSM_INSERT2END] = (1.0 - ratio_end) / (tsm->len-1);
    tsm->unit[i]->transition[TSM_DELETE2END] = (1.0 - ratio_end)/ (tsm->len-1);
  }

  tsm->unit[0]->transition[TSM_START2MATCH] = prob_into_model * (ratio_start);
  tsm->unit[0]->transition[TSM_START2DELETE] = prob_into_model * (ratio_start);
  /*  tsm->unit[0]->transition[TSM_START2INSERT] = prob_into_model * (ratio_start); */
  tsm->unit[tsm->len-1]->transition[TSM_MATCH2END] = (ratio_end);
  tsm->unit[tsm->len-1]->transition[TSM_DELETE2END] = (ratio_end);
} 

/* Function:  force_hmmfs_ThreeStateModel(tsm)
 *
 * Descrip:    places the probability at start end to precisely match
 *             hmmfs code.
 *
 *
 *
 * Arg:        tsm [UNKN ] ThreeStateModel to be 'forced' [ThreeStateModel *]
 *
 */
# line 316 "threestatemodel.dy"
void force_hmmfs_ThreeStateModel(ThreeStateModel * tsm)
{
  int i;
  double prob;
  
  prob = 1.0 - (1000. / 1001. );

  for(i=0;i<tsm->len;i++) {
    tsm->unit[i]->transition[TSM_START2MATCH] =  prob*0.5/(2.0 * (tsm->len-1)) ;
    tsm->unit[i]->transition[TSM_START2DELETE] = 0.0;
    tsm->unit[i]->transition[TSM_START2INSERT] = 0.0;
    tsm->unit[i]->transition[TSM_MATCH2END] =  1.0 / (tsm->len-1); 
    tsm->unit[i]->transition[TSM_INSERT2END] = 0.0;
    tsm->unit[i]->transition[TSM_DELETE2END] = 0.0;
  }

  tsm->unit[0]->transition[TSM_START2MATCH] = prob*0.5/2.0;
  tsm->unit[tsm->len-1]->transition[TSM_MATCH2END] = 1.0;
} 


/* Function:  ThreeStateScore_from_ThreeStateModel(tsm)
 *
 * Descrip:    Converts the three probability form to the score form.
 *
 *             There is real complications to this due to the fact that the prob
 *             model is a "push" model, that is MATCH2MATCH[i] means the match from i
 *             to i+1, whereas the DP routines rely on a "pull" model, ie
 *             MATCH2MATCH[i] is i-1 to i.
 *
 *             This routines does the conversion from push to pull, essentially offsetting
 *             the Match2 and Delete2 movements.
 *
 *
 * Arg:        tsm [READ ] ThreeStateModel probability form [ThreeStateModel *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateScore     *]
 *
 */
# line 350 "threestatemodel.dy"
ThreeStateScore     * ThreeStateScore_from_ThreeStateModel(ThreeStateModel * tsm)
{
  register int i;
  ThreeStateScore * tss;

  tss = ThreeStateScore_alloc_len(tsm->len);

  add_ThreeStateScore(tss,ThreeStateScoreUnit_from_ThreeStateUnit(NULL,tsm->unit[0]));
  for(i=1;i<tsm->len;i++) 
    add_ThreeStateScore(tss,ThreeStateScoreUnit_from_ThreeStateUnit(tsm->unit[i-1],tsm->unit[i]));

  tss->name = stringalloc(tsm->name);
  if( tsm->accession != NULL )
    tss->accession = stringalloc(tsm->accession);

  return tss;
}

/* Function:  ThreeStateScoreUnit_from_ThreeStateUnit(prev,tsu)
 *
 * Descrip:    Converts a three state model unit to a score unit
 *
 *             becuase of the push to  pull conversion needed, it needs the previous
 *             unit. for the first unit this is NULL and 0s are placed in the correct
 *             places in the transitions
 *
 *
 * Arg:        prev [UNKN ] Undocumented argument [ThreeStateUnit *]
 * Arg:         tsu [UNKN ] Undocumented argument [ThreeStateUnit *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateScoreUnit *]
 *
 */
# line 376 "threestatemodel.dy"
ThreeStateScoreUnit * ThreeStateScoreUnit_from_ThreeStateUnit(ThreeStateUnit * prev,ThreeStateUnit * tsu)
{
  ThreeStateScoreUnit * out;

  out = ThreeStateScoreUnit_alloc();
  if( out == NULL )
    return NULL;

  Probability2Score_move(tsu->match_emission,out->match,ALPHABET_SIZE);
  Probability2Score_move(tsu->insert_emission,out->insert,ALPHABET_SIZE);


  if( prev != NULL ) {
    out->trans[TSM_MATCH2MATCH] = Probability2Score(prev->transition[TSM_MATCH2MATCH]);
    out->trans[TSM_INSERT2MATCH] = Probability2Score(prev->transition[TSM_INSERT2MATCH]);
    out->trans[TSM_DELETE2MATCH] = Probability2Score(prev->transition[TSM_DELETE2MATCH]);

    out->trans[TSM_MATCH2DELETE] = Probability2Score(prev->transition[TSM_MATCH2DELETE]);
    out->trans[TSM_INSERT2DELETE] = Probability2Score(prev->transition[TSM_INSERT2DELETE]);
    out->trans[TSM_DELETE2DELETE] = Probability2Score(prev->transition[TSM_DELETE2DELETE]);
  } else {

    out->trans[TSM_MATCH2MATCH] = NEGI;
    out->trans[TSM_INSERT2MATCH] = NEGI;
    out->trans[TSM_DELETE2MATCH] = NEGI;

    out->trans[TSM_MATCH2DELETE] = NEGI;
    out->trans[TSM_INSERT2DELETE] = NEGI;
    out->trans[TSM_DELETE2DELETE] = NEGI;
  }

  out->trans[TSM_MATCH2INSERT] = Probability2Score(tsu->transition[TSM_MATCH2INSERT]);
  out->trans[TSM_INSERT2INSERT] = Probability2Score(tsu->transition[TSM_INSERT2INSERT]);
  out->trans[TSM_DELETE2INSERT] = Probability2Score(tsu->transition[TSM_DELETE2INSERT]);
  
  out->trans[TSM_START2MATCH] = Probability2Score(tsu->transition[TSM_START2MATCH]);
  out->trans[TSM_START2INSERT] = Probability2Score(tsu->transition[TSM_START2INSERT]);
  out->trans[TSM_START2DELETE] = Probability2Score(tsu->transition[TSM_START2DELETE]);
  
  out->trans[TSM_MATCH2END] = Probability2Score(tsu->transition[TSM_MATCH2END]);
  out->trans[TSM_INSERT2END] = Probability2Score(tsu->transition[TSM_INSERT2END]);
  out->trans[TSM_DELETE2END] = Probability2Score(tsu->transition[TSM_DELETE2END]);

  return out;
}


/* Function:  show_ThreeStateModel(tsm,ofp)
 *
 * Descrip:    shows pretty ugly debugging type format.
 *
 *             Pretty useless
 *
 *
 * Arg:        tsm [UNKN ] Undocumented argument [ThreeStateModel *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 428 "threestatemodel.dy"
void show_ThreeStateModel(ThreeStateModel * tsm,FILE * ofp)
{
  register int i;

  for(i=0;i<tsm->len;i++) {
    fprintf(ofp,"Position %d %c\n",i,tsm->unit[i]->display_char);
    show_ThreeStateUnit(tsm->unit[i],ofp);
  }
}

/* Function:  show_ThreeStateUnit(tsu,ofp)
 *
 * Descrip:    for debugging problems
 *
 *
 * Arg:        tsu [UNKN ] Undocumented argument [ThreeStateUnit *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 442 "threestatemodel.dy"
void show_ThreeStateUnit(ThreeStateUnit * tsu,FILE * ofp)
{
  fprintf(ofp,"Match ");
  show_Probability_array(tsu->match_emission,26,ofp);
  fprintf(ofp,"\nInsert ");
  show_Probability_array(tsu->insert_emission,26,ofp);
  fprintf(ofp,"\nTransition ");
  show_Probability_array(tsu->transition,TRANSITION_LEN,ofp);
  fprintf(ofp,"\n");
  
}

/* Function:  pseudo_Protein_from_ThreeStateModel(tsm)
 *
 * Descrip:    Makes a protein sequence out of the display characters.
 *             Not very useful!
 *
 *
 *
 * Arg:        tsm [UNKN ] Undocumented argument [ThreeStateModel *]
 *
 * Return [UNKN ]  Undocumented return value [Protein *]
 *
 */
# line 459 "threestatemodel.dy"
Protein * pseudo_Protein_from_ThreeStateModel(ThreeStateModel * tsm)
{
  int i;

  Sequence * seq;

  seq = Sequence_alloc();
  seq->name = stringalloc(tsm->name);
  seq->seq = ckcalloc(tsm->len+1,sizeof(char));

  for(i=0;i<tsm->len;i++) {
    seq->seq[i] = tsm->unit[i]->display_char;
  }
  seq->seq[i]='\0';
  make_len_type_Sequence(seq);
  seq->type = SEQUENCE_PROTEIN;

  return Protein_from_Sequence(seq);
}



/* Function:  ThreeStateModel_from_half_bit_Sequence(pro,mat,rm,gap,ext)
 *
 * Descrip:    Makes a local three-state-model from a sequence.  this is scary
 *             hackery, assumming that the matrix is half-bits and normalising in a
 *             *very* wrong way to get "probabilities" out.
 *
 *             Works though
 *
 *
 * Arg:        pro [READ ] protein sequence [Protein *]
 * Arg:        mat [READ ] comparison matrix to use [CompMat *]
 * Arg:         rm [READ ] random model which you assumme the matrix was built with [RandomModel *]
 * Arg:        gap [READ ] gap open penalty [int]
 * Arg:        ext [READ ] gap ext penalty [int]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateModel *]
 *
 */
# line 494 "threestatemodel.dy"
ThreeStateModel * ThreeStateModel_from_half_bit_Sequence(Protein * pro,CompMat * mat,RandomModel * rm,int gap,int ext)
{
  Sequence * seq;
  register int i;
  ThreeStateModel * out;


  
  if( pro == NULL || mat == NULL || rm == NULL ) {
    warn("you have passed through NULL objects in trying to make TSM from sequence");
    return NULL;
  }

  if( gap > 0 || ext > 0 ) {
    warn("You have passed in gap and extension penalties of > 0 Gap %d Ext %d. Expecting them to be negated. Giving you back an error!",gap,ext);
    return NULL;
  }


  seq = pro->baseseq;

  out = ThreeStateModel_alloc_len(seq->len);

  if( seq->name != NULL )
    out->name = stringalloc(seq->name);
  else out->name = stringalloc("NoName");
  

  out->rm = hard_link_RandomModel(rm);

  for(i=0;i<seq->len;i++) { 

    add_ThreeStateModel(out,ThreeStateUnit_from_half_bit_aminoacid(seq->seq[i],mat,rm,gap,ext));
  }
		       
  return out;
}

/* Function:  ThreeStateUnit_from_half_bit_aminoacid(aa,mat,rm,gap,ext)
 *
 * Descrip:    The internal protein->hmm conversion routine.
 *
 *
 * Arg:         aa [UNKN ] Undocumented argument [char]
 * Arg:        mat [UNKN ] Undocumented argument [CompMat *]
 * Arg:         rm [UNKN ] Undocumented argument [RandomModel *]
 * Arg:        gap [UNKN ] Undocumented argument [int]
 * Arg:        ext [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateUnit *]
 *
 */
# line 536 "threestatemodel.dy"
ThreeStateUnit *  ThreeStateUnit_from_half_bit_aminoacid(char aa,CompMat * mat,RandomModel * rm,int gap,int ext)
{
  ThreeStateUnit * tsu;
  register int i;
  Probability go;
  Probability ge;
  Probability rnd_error;

  tsu = ThreeStateUnit_alloc();
  aa = toupper((int)aa);
  /* maybe do more checking on aa? */


  for(i=0;i<26;i++) {
    tsu->match_emission[i]  = halfbit2Probability(fail_safe_CompMat_access(mat,aa-'A',i)) * rm->aminoacid[i];
    tsu->insert_emission[i] = rm->aminoacid[i];
  }

  rnd_error = renormalise_Probability_array(tsu->match_emission,26);

  if( fabs(rnd_error) > 0.5 ) {
    warn("Very bad problem - got a stupid error %.2f after renormalisation from the null model with the match state on aa %c\n",rnd_error,aa);
  }

  rnd_error = renormalise_Probability_array(tsu->insert_emission,26);
  if( fabs(rnd_error) > 0.5 ) {
    warn("Very bad problem - got a stupid error %.2f after renormalisation from the null model\n",rnd_error);
  }


  ge = halfbit2Probability(ext);

  go = halfbit2Probability(gap);

  /** 1-ge is the gap closing penalty. So we have to remove that from the go **/

  go = go / (1-ge);

  tsu->transition[TSM_MATCH2INSERT] = tsu->transition[TSM_MATCH2DELETE] = go;
  tsu->transition[TSM_INSERT2MATCH] = tsu->transition[TSM_DELETE2MATCH] = 1-ge;
  tsu->transition[TSM_INSERT2INSERT] = tsu->transition[TSM_DELETE2DELETE] = ge;
  tsu->transition[TSM_INSERT2DELETE] = tsu->transition[TSM_DELETE2INSERT] = 0.0;
  tsu->transition[TSM_MATCH2MATCH]  = 1-(2*go);
  tsu->transition[TSM_START2MATCH]  = 1.0;
  tsu->transition[TSM_MATCH2END] = 1.0;

  tsu->display_char = aa;

  return tsu;
}
    

/* Function:  global_ThreeStateModel_from_half_bit_Sequence(pro,mat,rm,gap,ext)
 *
 * Descrip:    Makes a global three-state-model from a sequence.
 *             Like the local version, this is scary hackery, assumming
 *             that the matrix is half-bits and normalising in a *very*
 *             wrong way to get "probabilities" out.
 *
 *             Works though
 *
 *
 * Arg:        pro [READ ] protein sequence [Protein *]
 * Arg:        mat [READ ] comparison matrix to use [CompMat *]
 * Arg:         rm [READ ] random model which you assumme the matrix was built with [RandomModel *]
 * Arg:        gap [READ ] gap open penalty [int]
 * Arg:        ext [READ ] gap ext penalty [int]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateModel *]
 *
 */
# line 602 "threestatemodel.dy"
ThreeStateModel * global_ThreeStateModel_from_half_bit_Sequence(Protein * pro,CompMat * mat,RandomModel * rm,int gap,int ext)
{

  register int i;
  ThreeStateModel * out;
  ThreeStateUnit * temp;
  Sequence * seq;

  if( pro == NULL || mat == NULL || rm == NULL ) {
    warn("you have passed through NULL objects in trying to make TSM from sequence");
    return NULL;
  }

  seq = pro->baseseq;


  out = ThreeStateModel_alloc_len(seq->len);

  if( seq->name != NULL )
    out->name = stringalloc(seq->name);
  else out->name = stringalloc("NoName");

  for(i=0;i<seq->len;i++) { 
    temp = ThreeStateUnit_from_half_bit_aminoacid(seq->seq[i],mat,rm,gap,ext);
    temp->transition[TSM_START2MATCH] = 0.0;
    temp->transition[TSM_MATCH2END] = 0.0;
    add_ThreeStateModel(out,temp);
  }
  
  out->unit[0]->transition[TSM_START2MATCH] = 1.0;
  out->unit[out->len-1]->transition[TSM_MATCH2END] = 1.0;
  out->rm = hard_link_RandomModel(rm);

		       
  return out;
}


/* Function:  display_char_for_ThreeStateUnit(tsu)
 *
 * Descrip:    Makes a display character for
 *             this threestateunit by finding the
 *             most probable match emission
 *
 *
 * Arg:        tsu [UNKN ] Undocumented argument [ThreeStateUnit *]
 *
 * Return [UNKN ]  Undocumented return value [char]
 *
 */
# line 646 "threestatemodel.dy"
char display_char_for_ThreeStateUnit(ThreeStateUnit * tsu) 
{
  register int i;
  Probability p;
  int c;

  p = tsu->match_emission[0];
  c = 0;

  for(i=1;i<ALPHABET_SIZE;i++)
    if( tsu->match_emission[i] > p ) {
      p =  tsu->match_emission[i];
      c = i;
    }

  return c + 'A';
}

/* Function:  display_char_in_ThreeStateModel(tsm)
 *
 * Descrip:    Makes the display chars in the tsm come from
 *             the most likely amino acid in the match position
 *
 *
 * Arg:        tsm [UNKN ] Undocumented argument [ThreeStateModel *]
 *
 */
# line 668 "threestatemodel.dy"
void display_char_in_ThreeStateModel(ThreeStateModel * tsm)
{
  register int i;

  for(i=0;i<tsm->len;i++)
    tsm->unit[i]->display_char = display_char_for_ThreeStateUnit(tsm->unit[i]);
}

/* Function:  fold_RandomModel_into_ThreeStateModel(tsm,rm)
 *
 * Descrip:    divides the emission and insertion scores
 *             by the random model for use in log-odd
 *             implementations
 *
 *
 * Arg:        tsm [UNKN ] Undocumented argument [ThreeStateModel *]
 * Arg:         rm [UNKN ] Undocumented argument [RandomModel *]
 *
 */
# line 681 "threestatemodel.dy"
void fold_RandomModel_into_ThreeStateModel(ThreeStateModel * tsm,RandomModel * rm)
{
  register int i;
  register int j;

  assert(tsm);
  assert(rm);
  for(i=0;i<tsm->len;i++) {
    auto ThreeStateUnit * tsu;
    tsu = tsm->unit[i];

    for(j=0;j<26;j++) {
      if( rm->aminoacid[j] < 0.00000001 ) {
	warn("While trying to fold in random model, amino acid %d [%c] was below zero, ignoring",j,'A'+j);
	continue;
      }
      
      tsu->match_emission[j]  = tsu->match_emission[j] / rm->aminoacid[j];
      tsu->insert_emission[j] = tsu->insert_emission[j] / rm->aminoacid[j];
    }
  }
}

/* Function:  add_sensible_start_end_global_for_HMMer(tsm)
 *
 * Descrip:    added start 0.5, 0.25,0.25 for hmmls type
 *
 *             Deprecated
 *
 *
 * Arg:        tsm [UNKN ] Undocumented argument [ThreeStateModel *]
 *
 */
# line 709 "threestatemodel.dy"
void add_sensible_start_end_global_for_HMMer(ThreeStateModel * tsm)
{
  ThreeStateUnit * end;
  
  tsm->unit[0]->transition[TSM_START2MATCH] = 0.5;
  tsm->unit[0]->transition[TSM_START2INSERT] = 0.25;
  tsm->unit[0]->transition[TSM_START2DELETE] = 0.25;
  
  /*** these should not add to one due to sum in states ***/

  end = tsm->unit[tsm->len-1];

  end->transition[TSM_MATCH2END] =  (1 - end->transition[TSM_MATCH2INSERT]) ;
  end->transition[TSM_DELETE2END] = (1 - end->transition[TSM_DELETE2INSERT]) ;
  end->transition[TSM_INSERT2END] = (1 - end->transition[TSM_INSERT2INSERT]);
    

}

/* Function:  add_sensible_start_end_looping_for_HMMer(tsm,loop)
 *
 * Descrip:    adds start/end for loop probabilities coming
 *             in
 *
 *             Deprecated
 *
 *
 * Arg:         tsm [UNKN ] Undocumented argument [ThreeStateModel *]
 * Arg:        loop [UNKN ] Undocumented argument [Probability]
 *
 */
# line 734 "threestatemodel.dy"
void add_sensible_start_end_looping_for_HMMer(ThreeStateModel * tsm,Probability loop)
{
  tsm->unit[0]->transition[TSM_START2MATCH] = loop * 0.5;
  tsm->unit[0]->transition[TSM_START2INSERT] = loop * 0.25;
  tsm->unit[0]->transition[TSM_START2DELETE] = loop * 0.25;
  
  /*** these should not add to one due to sum in states ***/

  tsm->unit[tsm->len-1]->transition[TSM_MATCH2END] = 0.5;
  tsm->unit[tsm->len-1]->transition[TSM_DELETE2END] = 0.25;
  tsm->unit[tsm->len-1]->transition[TSM_INSERT2END] = 0.25;
    
}



  /**** I/O *****/



/* Function:  write_HMF_ThreeStateModel(tsm,ofp)
 *
 * Descrip:    writes the HMF format
 *
 *
 * Arg:        tsm [UNKN ] Undocumented argument [ThreeStateModel *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 757 "threestatemodel.dy"
void write_HMF_ThreeStateModel(ThreeStateModel * tsm,FILE * ofp)
{
  int i;

  fprintf(ofp,"ID %s\n",tsm->name != NULL ? tsm->name : "NoName");
  fprintf(ofp,"TR %4.2f\n",tsm->threshold);
  fprintf(ofp,"CC HMF format for HMM, each three state block starts UN after MO line\n");
  fprintf(ofp,"CC Alphabet A-Z\n");
  fprintf(ofp,"MO\n");

  for(i=0;i<tsm->len;i++) 
    write_HMF_ThreeStateUnit(tsm->unit[i],ofp);

  fprintf(ofp,"//\n");
}


/* Function:  write_HMF_ThreeStateUnit(tsu,ofp)
 *
 * Descrip:    writes each Unit line
 *
 *
 * Arg:        tsu [UNKN ] Undocumented argument [ThreeStateUnit *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 778 "threestatemodel.dy"
void write_HMF_ThreeStateUnit(ThreeStateUnit * tsu,FILE * ofp)
{
  fprintf(ofp,"UN %c\n",tsu->display_char);
  fprintf(ofp,"MA ");
  show_Probability_array_exp(tsu->match_emission,ALPHABET_SIZE,ofp);
  fprintf(ofp,"\nIN ");
  show_Probability_array_exp(tsu->insert_emission,ALPHABET_SIZE,ofp);
  fprintf(ofp,"\nTR ");
  show_Probability_array_exp(tsu->transition,TRANSITION_LEN,ofp);
  fprintf(ofp,"\n");
}

/* Function:  read_HMF_ThreeStateModel(ifp)
 *
 * Descrip:    reads the HMF format. Leaves the file ready
 *             to read the next format
 *
 *
 * Arg:        ifp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateModel *]
 *
 */
# line 794 "threestatemodel.dy"
ThreeStateModel * read_HMF_ThreeStateModel(FILE * ifp)
{
  ThreeStateUnit * temp;
  ThreeStateModel * out;
  char buffer[MAXLINE];
  char * name = NULL;
  char * runner;

  if( feof(ifp) || ferror(ifp) )
    return NULL;

  while( fgets(buffer,MAXLINE,ifp) != NULL ) {
    chop_newline(buffer);
    if( buffer[0] == '#' ) {
      continue;
    }
    if( strstartcmp(buffer,"ID") == 0 ) {
      name = stringalloc(strtok(buffer+3,spacestr));
      break;
    }
    warn("In reading HMF format, before ID line, ignoring [%s]",buffer);
  }

  if( feof(ifp) || ferror(ifp) || name == NULL)
    return NULL;

  
  
  out = ThreeStateModel_alloc_std();
  out->name = name;

  while( fgets(buffer,MAXLINE,ifp) != NULL ) {
    chop_newline(buffer);
    if( buffer[0] == '#' ) {
      continue;
    }
    if( strstartcmp(buffer,"MO") == 0 )
      break;
    else if( strstartcmp(buffer,"CC") == 0 )
      continue;
    else if( strstartcmp(buffer,"TR") == 0 ) {
      runner = strtok(buffer+2,spacestr);
      if( is_double_string(runner,&out->threshold) == FALSE) {
	warn("Could not interpret [%s] as a threshold",runner);
      }
    }
    else if( strstartcmp(buffer,"ID") == 0 ) {
      warn("Already got a name for [%s], ignoring [%s]",out->name,buffer);
    } else {
      warn("Did not understand %s",buffer);
    }
  }

  while( fgets(buffer,MAXLINE,ifp) != NULL ) {
    if( strstartcmp(buffer,"UN") == 0 )
      break;
    chop_newline(buffer);
    warn("You have a problem here: [%s] between MO and UN lines",buffer);
  }

  while( (temp = read_HMF_ThreeStateUnit(buffer,ifp)) != NULL ) {
    add_ThreeStateModel(out,temp);
    if( strstartcmp(buffer,"//") == 0 )
      break;
  }

  display_char_in_ThreeStateModel(out);

  return out;
}
  

/* Function:  read_HMF_ThreeStateUnit(line,ifp)
 *
 * Descrip:    reads a unit set of lines: 
 *
 *             Line is used as the buffer, so should be MAXLINE length. At the
 *             end, it will have a different characters (probably the next unit
 *             line) in the line buffer, so you should test that
 *
 *
 * Arg:        line [WRITE] pointer to a MAXLINE length buffer. [char *]
 * Arg:         ifp [UNKN ] input file [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateUnit *]
 *
 */
# line 877 "threestatemodel.dy"
ThreeStateUnit * read_HMF_ThreeStateUnit(char * line,FILE * ifp)
{
  ThreeStateUnit * out;
  char * runner;

  if( strstartcmp(line,"UN") != 0 ) {
    warn("Attempting to read a ThreeStateUnit with a non UN line: [%s]",line);
    return NULL;
  }

  out = blank_ThreeStateUnit();
  /* display char */

  runner = strtok(line+2,spacestr);
  if( runner != NULL ) {
    out->display_char = *runner;
  }

  while ( fgets(line,MAXLINE,ifp) != NULL ) {
    if( strstartcmp(line,"//") == 0 )
      break;
    if( strstartcmp(line,"UN") == 0 )
      break;
    if( strstartcmp(line,"MA") == 0 ) {
      if( read_Probability_array(out->match_emission,ALPHABET_SIZE,line+3) == FALSE ) {
	warn("Unable to read Match emission line");
      }
    } else if ( strstartcmp(line,"IN") == 0) {
      if( read_Probability_array(out->insert_emission,ALPHABET_SIZE,line+3) == FALSE ) {
	warn("Unable to read Insert emission line");
      }
    } else if ( strstartcmp(line,"TR") == 0 ) {
      if( read_Probability_array(out->transition,TRANSITION_LEN,line+3) == FALSE ) {
	warn("Unable to read Insert emission line");
      }
    } else {
      chop_newline(line);
      warn("Could not understand line [%s] in reading threestateunit",line);
    }

  }

  return out;
}

/* Function:  write_HMMer_1_7_ascii_ThreeStateModel(tsm,ofp)
 *
 * Descrip:    writes a HMMer version 1.7 (also ok with 1.8) file
 *
 *
 * Arg:        tsm [UNKN ] Undocumented argument [ThreeStateModel *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 925 "threestatemodel.dy"
void write_HMMer_1_7_ascii_ThreeStateModel(ThreeStateModel * tsm,FILE * ofp)
{
  register int i;

  fprintf(ofp,"# HMM v1.7\n");
  fprintf(ofp,"%d     # Length of model\n",tsm->len);
  fprintf(ofp,"%d     # length of alphabet\n",26);
  fprintf(ofp,"3      # alphabet type\n");
  fprintf(ofp,"%s     # alphabet\n",std_alphabet);
  fprintf(ofp,"no     # optional\nno      # optional \n");
  for(i=0;i<tsm->len;i++) 
    write_HMMer_1_7_ascii_ThreeStateUnit(tsm->unit[i],i,ofp);

}

/* Function:  write_HMMer_1_7_ascii_ThreeStateUnit(tsu,no,ofp)
 *
 * Descrip:    the internal for the ThreeStateModel write
 *
 *
 * Arg:        tsu [UNKN ] Undocumented argument [ThreeStateUnit *]
 * Arg:         no [UNKN ] Undocumented argument [int]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 944 "threestatemodel.dy"
void write_HMMer_1_7_ascii_ThreeStateUnit(ThreeStateUnit * tsu,int no,FILE * ofp)
{
  register int i;

  fprintf(ofp,"###MATCH_STATE %d ( ) ( )\n",no);
  fprintf(ofp,"%f\t# t_m%d\n",tsu->transition[TSM_MATCH2MATCH],no+1);
  fprintf(ofp,"%f\t# t_d%d\n",tsu->transition[TSM_MATCH2DELETE],no+1);
  fprintf(ofp,"%f\t# t_i%d\n",tsu->transition[TSM_MATCH2INSERT],no);
  for(i=0;i<26;i++) {
    fprintf(ofp,"%f\t# Symbol %c probability\n",tsu->match_emission[i],'A'+i);
  }

  fprintf(ofp,"###DELETE_STATE %d\n",no);
  fprintf(ofp,"%f\t# t_m%d\n",tsu->transition[TSM_DELETE2MATCH],no+1);
  fprintf(ofp,"%f\t# t_d%d\n",tsu->transition[TSM_DELETE2DELETE],no+1);
  fprintf(ofp,"%f\t# t_i%d\n",tsu->transition[TSM_DELETE2INSERT],no);

  fprintf(ofp,"###INSERT_STATE %d\n",no);
  fprintf(ofp,"%f\t# t_m%d\n",tsu->transition[TSM_INSERT2MATCH],no+1);
  fprintf(ofp,"%f\t# t_d%d\n",tsu->transition[TSM_INSERT2DELETE],no+1);
  fprintf(ofp,"%f\t# t_i%d\n",tsu->transition[TSM_INSERT2INSERT],no);
  for(i=0;i<26;i++) {
    fprintf(ofp,"%f\t# Symbol %c probability\n",tsu->insert_emission[i],'A'+i);
  }

}

/* Function:  read_HMMer_1_7_ascii_file(filename)
 *
 * Descrip:    reads a HMMer ascii version 1.7 (1.8) file from filename.
 *
 *
 *
 * Arg:        filename [UNKN ] the name fo the hmmer file [char *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateModel *]
 *
 */
# line 977 "threestatemodel.dy"
ThreeStateModel * read_HMMer_1_7_ascii_file(char * filename)
{
  FILE * ifp;
  ThreeStateModel * out;

  ifp = openfile(filename,"r");
  if( ifp == NULL ) {
    warn("Could not open file %s for read HMMer 17 ",filename);
    return NULL;
  }

  out = read_HMMer_1_7_ascii(ifp);

  ckfree(out->name);
  if( strcmp(filename,"-") == 0 ) {
    /* comes from stdin, we're fucked to find the name */
    out->name = stringalloc("HMMerModel");
  } else {
    if( strchr(filename,'/') != NULL) {
      filename = filename + strlen(filename) -1;
      for(;*filename != '/';filename--)
	;
      filename++;
    }
    out->name = stringalloc(filename);
  }

  fclose(ifp);

  return out;
}


/* Function:  convert_push_model_to_pull_model(tsm)
 *
 * Descrip:    Not sure if this is needed now, as score model will be pull, and
 *             probability model will be push.
 *
 *
 *
 * Arg:        tsm [UNKN ] Undocumented argument [ThreeStateModel *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateModel *]
 *
 */
# line 1016 "threestatemodel.dy"
ThreeStateModel * convert_push_model_to_pull_model(ThreeStateModel * tsm)
{
  register int i;
  ThreeStateModel * out;


  out = ThreeStateModel_alloc_len(tsm->len-1);

  for(i=1;i<tsm->len;i++) 
    add_ThreeStateModel(out,convert_push_unit_to_pull_unit(tsm->unit[i-1],tsm->unit[i]));

  if( tsm->name != NULL)
    out->name = stringalloc(tsm->name);

  return out;
}

/* Function:  convert_push_unit_to_pull_unit(prev,this)
 *
 * Descrip:    Not sure if this is needed now, as score model will be pull, and
 *             probability model will be push.
 *
 *
 *
 * Arg:        prev [UNKN ] Undocumented argument [ThreeStateUnit *]
 * Arg:        this [UNKN ] Undocumented argument [ThreeStateUnit *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateUnit  *]
 *
 */
# line 1039 "threestatemodel.dy"
ThreeStateUnit  * convert_push_unit_to_pull_unit(ThreeStateUnit * prev,ThreeStateUnit * this)
{
  ThreeStateUnit * out;


  out = ThreeStateUnit_alloc();

  
  Probability_move(this->match_emission ,out->match_emission ,ALPHABET_SIZE);
  Probability_move(this->insert_emission,out->insert_emission,ALPHABET_SIZE);

  set_Probability_array(out->transition,0.0,TRANSITION_LEN);

  out->transition[TSM_MATCH2MATCH]   = prev->transition[TSM_MATCH2MATCH];
  out->transition[TSM_DELETE2MATCH]  = prev->transition[TSM_DELETE2MATCH];
  out->transition[TSM_INSERT2MATCH]  = prev->transition[TSM_INSERT2MATCH];
  out->transition[TSM_MATCH2INSERT]  = this->transition[TSM_MATCH2INSERT];
  out->transition[TSM_INSERT2INSERT] = this->transition[TSM_INSERT2INSERT];
  out->transition[TSM_DELETE2INSERT] = this->transition[TSM_DELETE2INSERT];
  out->transition[TSM_MATCH2DELETE]  = prev->transition[TSM_MATCH2DELETE];
  out->transition[TSM_INSERT2DELETE] = prev->transition[TSM_INSERT2DELETE];
  out->transition[TSM_DELETE2DELETE] = prev->transition[TSM_DELETE2DELETE];

  /*** not dealing with start/end things here ***/

  return out;
}
 

/* Function:  read_HMMer_1_7_ascii(ifp)
 *
 * Descrip:    Basic function to read HMMer version 1.7(1.8) files.
 *
 *
 * Arg:        ifp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateModel *]
 *
 */
# line 1071 "threestatemodel.dy"
ThreeStateModel * read_HMMer_1_7_ascii(FILE * ifp) 
{
  char buffer[MAXLINE];
  char * runner;
  char * alphabet;
  int length;
  int count;
  ThreeStateModel * out;
  ThreeStateUnit  * temp;
  
  /*** read top of model   ***/
  /*** then we will loop   ***/
  /*** over all the states ***/


  if( fgets(buffer,MAXLINE,ifp) == NULL ) {
    warn("Trying to read HMMer ascii file, but no first line!");
    return NULL;
  }

  if( strstartcmp(buffer,"# HMM") != 0 ) {
    warn("In reading HMMer ascii file, expecting first line to start [# HMM], starts as [%s]... can't handle",buffer);
    return NULL;
  }

  if( strstartcmp(buffer,"# HMM v1.7") != 0)  {
    warn("In reading HMMer ascii file, expecting a version 1.7 file, got a version [%s] file... trying to read",buffer);
  }

  /*** length next ***/

  if( fgets(buffer,MAXLINE,ifp) == NULL ) {
    warn("Trying to read HMMer ascii file, length line (2nd line) not there");
    return NULL;
  }

  if( (runner=strtok(buffer,spacestr)) == NULL ) {
    warn("Trying to read HMMer ascii file, length line (2nd line) has no length!");
    return NULL;
  }

  length = atoi(runner);
  if( length <= 0 || length > 4000 ) {
    warn("Picked up length [%s] which was interpreted as the number [%d], which seems unfeasible. S'ok... can work without the length",runner,length);
    length = 100;
  }

  /*** get the alphabet length... no point in checking ***/


  if( fgets(buffer,MAXLINE,ifp) == NULL ) {
    warn("Trying to read HMMer ascii file, length of alphabet (3rd line) not there");
    return NULL;
  }

  if( fgets(buffer,MAXLINE,ifp) == NULL ) {
    warn("Trying to read HMMer ascii file, type of alphabet (4th line) not there");
    return NULL;
  }
     
  
  if( fgets(buffer,MAXLINE,ifp) == NULL ) {
    warn("Trying to read HMMer ascii file, alphabet (5th line) not there");
    return NULL;
  }

  if( (runner=strtok(buffer,spacestr)) == NULL ) {
    warn("HMMer ascii read, alphabet line (5th line) has no alphabet!");
    return NULL;
  }

  alphabet = stringalloc(runner);

  
  if( fgets(buffer,MAXLINE,ifp) == NULL ) {
    warn("Trying to read HMMer ascii file, reference (6th line) not there");
    return NULL;
  }
  
  if( fgets(buffer,MAXLINE,ifp) == NULL ) {
    warn("Trying to read HMMer ascii file, reference (7th line) not there");
    return NULL;
  }
  

  /*** ok, ready to rock, we have the length and alphabet, ***/
  /*** now loop over the states ...                        ***/

  out = ThreeStateModel_alloc_len(length);

  for(count=0;count < length+1 && !feof(ifp) && !(ferror(ifp));count++ ) {
    temp = read_HMMer_1_7_ThreeStateUnit(alphabet,ifp);
    if( temp == NULL ) {
      warn("In HMMer read, On unit number [%d], unable to get unit...",count);
      break;
    }
    if( count == 0 )
      continue; /** skip first "dummy" state **/

    else add_ThreeStateModel(out,temp);
  }

  

  out->name = stringalloc("HMMer Model");
  out->alphabet = alphabet;

  /*** by default, make it a global model ***/

  force_global_model(out,1.0);

  return out;
}

/* Function:  read_HMMer_1_7_ThreeStateUnit(alphabet,ifp)
 *
 * Descrip:    Function to read a single unit from a 1.7(1.8) file
 *
 *
 *
 * Arg:        alphabet [UNKN ] Undocumented argument [char *]
 * Arg:             ifp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateUnit *]
 *
 */
# line 1190 "threestatemodel.dy"
ThreeStateUnit * read_HMMer_1_7_ThreeStateUnit(char * alphabet,FILE * ifp) 
{
  char buffer[MAXLINE];
  char * runner;
  char * al;
  ThreeStateUnit * out;

  if( fgets(buffer,MAXLINE,ifp) == NULL ) {
    /*** could be the end of the model. Return silently for the next stage ***/
    return NULL;
  }

  if( strstartcmp(buffer,"###MATCH_STATE") != 0 ) {
    warn("Start of ThreeStateUnit HMMr read... does not match [###MATCH_STATE], instead... [%s]",buffer);
    return NULL;
  }

  out = blank_ThreeStateUnit();
  if( out == NULL )
    return NULL; /** warning already issued **/

  /*** ok read in the first three transitions ***/
  /*** transitions for MATCH                  ***/

  if( fgets(buffer,MAXLINE,ifp) == NULL ) {
    warn("In reading HMMr ThreeStateUnit, got no m2m line");
    free_ThreeStateUnit(out);
    return NULL;
  }

  if( (runner = strtok(buffer,spacestr)) == NULL ) {
    warn("In reading HMMr ThreeStateUnit, got m2m line, but no transition worked out");
    return NULL;
  }

  out->transition[TSM_MATCH2MATCH] = atof(runner);
  

  if( fgets(buffer,MAXLINE,ifp) == NULL ) {
    warn("In reading HMMr ThreeStateUnit, got no m2d line");
    free_ThreeStateUnit(out);
    return NULL;
  }

  if( (runner = strtok(buffer,spacestr)) == NULL ) {
    warn("In reading HMMr ThreeStateUnit, got m2d line, but no transition worked out");
    return NULL;
  }

  out->transition[TSM_MATCH2DELETE] = atof(runner);

  if( fgets(buffer,MAXLINE,ifp) == NULL ) {
    warn("In reading HMMr ThreeStateUnit, got no m2i line");
    free_ThreeStateUnit(out);
    return NULL;
  }

  if( (runner = strtok(buffer,spacestr)) == NULL ) {
    warn("In reading HMMr ThreeStateUnit, got m2i line, but no transition worked out");
    return NULL;
  }

  out->transition[TSM_MATCH2INSERT] = atof(runner);
  
  /*** match emission probs  ***/

  for(al=alphabet;*al;al++) {
    if( fgets(buffer,MAXLINE,ifp) == NULL ) {
      warn("In reading HMMr threestateunit, the file ended on the match emission of %c",*al);
      free_ThreeStateUnit(out);
      return NULL;
    }
    runner = strtok(buffer,spacestr);
    /**** ASSUMMING a 26 alphabet now... Aaaaah ****/
    out->match_emission[ toupper((int)*al)-'A'] = atof(runner);
  }


  /*** now to delete probabilities ***/
  

  if( fgets(buffer,MAXLINE,ifp) == NULL ) {
    warn("In reading HMMr ThreeStateUnit, got no start line!");
    return NULL;
  }

  if( strstartcmp(buffer,"###DELETE_STATE") != 0 ) {
    warn("Start of delete ThreeStateUnit HMMr read... does not match [###DELETE_STATE], instead... [%s]",buffer);
    return NULL;
  }


  /*** ok read in the first three transitions ***/
  /*** transitions for DELETE                  ***/

  if( fgets(buffer,MAXLINE,ifp) == NULL ) {
    warn("In reading HMMr ThreeStateUnit, got no d2m line");
    free_ThreeStateUnit(out);
    return NULL;
  }

  if( (runner = strtok(buffer,spacestr)) == NULL ) {
    warn("In reading HMMr ThreeStateUnit, got d2m line, but no transition worked out");
    return NULL;
  }

  out->transition[TSM_DELETE2MATCH] = atof(runner);
  

  if( fgets(buffer,MAXLINE,ifp) == NULL ) {
    warn("In reading HMMr ThreeStateUnit, got no d2d line");
    free_ThreeStateUnit(out);
    return NULL;
  }

  if( (runner = strtok(buffer,spacestr)) == NULL ) {
    warn("In reading HMMr ThreeStateUnit, got d2d line, but no transition worked out");
    return NULL;
  }

  out->transition[TSM_DELETE2DELETE] = atof(runner);

  if( fgets(buffer,MAXLINE,ifp) == NULL ) {
    warn("In reading HMMr ThreeStateUnit, got no d2i line");
    free_ThreeStateUnit(out);
    return NULL;
  }

  if( (runner = strtok(buffer,spacestr)) == NULL ) {
    warn("In reading HMMr ThreeStateUnit, got d2i line, but no transition worked out");
    return NULL;
  }

  out->transition[TSM_DELETE2INSERT] = atof(runner);



  /***** now insert probs *****/

  if( fgets(buffer,MAXLINE,ifp) == NULL ) {
    warn("In reading HMMr ThreeStateUnit, got no start line!");
    return NULL;
  }

  if( strstartcmp(buffer,"###INSERT_STATE") != 0 ) {
    warn("Start of ThreeStateUnit HMMr read... does not match [###INSERT_STATE], instead... [%s]",buffer);
    return NULL;
  }

  /*** ok read in the first three transitions ***/
  /*** transitions for INSERT                 ***/

  if( fgets(buffer,MAXLINE,ifp) == NULL ) {
    warn("In reading HMMr ThreeStateUnit, got no i2m line");
    free_ThreeStateUnit(out);
    return NULL;
  }

  if( (runner = strtok(buffer,spacestr)) == NULL ) {
    warn("In reading HMMr ThreeStateUnit, got i2m line, but no transition worked out");
    return NULL;
  }

  out->transition[TSM_INSERT2MATCH] = atof(runner);
  

  if( fgets(buffer,MAXLINE,ifp) == NULL ) {
    warn("In reading HMMr ThreeStateUnit, got no i2d line");
    free_ThreeStateUnit(out);
    return NULL;
  }

  if( (runner = strtok(buffer,spacestr)) == NULL ) {
    warn("In reading HMMr ThreeStateUnit, got i2d line, but no transition worked out");
    return NULL;
  }

  out->transition[TSM_INSERT2DELETE] = atof(runner);

  if( fgets(buffer,MAXLINE,ifp) == NULL ) {
    warn("In reading HMMr ThreeStateUnit, got no i2i line");
    free_ThreeStateUnit(out);
    return NULL;
  }

  if( (runner = strtok(buffer,spacestr)) == NULL ) {
    warn("In reading HMMr ThreeStateUnit, got i2i line, but no transition worked out");
    return NULL;
  }

  out->transition[TSM_INSERT2INSERT] = atof(runner);
  
  /*** insert emission probs  ***/

  for(al=alphabet;*al;al++) {
    if( fgets(buffer,MAXLINE,ifp) == NULL ) {
      warn("In reading HMMr threestateunit, the file ended on the insert emission of %c",*al);
      free_ThreeStateUnit(out);
      return NULL;
    }
    runner = strtok(buffer,spacestr);
    /**** ASSUMMING a 26 alphabet now... Aaaaah ****/
    out->insert_emission[ toupper((int)*al)-'A'] = atof(runner);
  }


  return out;
}


/* Function:  blank_ThreeStateUnit(void)
 *
 * Descrip:    Makes a set to 0.0 unit
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateUnit *]
 *
 */
# line 1404 "threestatemodel.dy"
ThreeStateUnit * blank_ThreeStateUnit(void)
{
  ThreeStateUnit * out;
  register int i;

  out = ThreeStateUnit_alloc();
  if(out == NULL ) 
    return NULL;

  for(i=0;i<26;i++) 
    out->match_emission[i] = out->insert_emission[i] = 0.0;
  
  for(i=0;i<TRANSITION_LEN;i++)
    out->transition[i] = 0.0;

  return out;
}


# line 1476 "threestatemodel.c"
/* Function:  hard_link_ThreeStateUnit(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ThreeStateUnit *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateUnit *]
 *
 */
ThreeStateUnit * hard_link_ThreeStateUnit(ThreeStateUnit * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a ThreeStateUnit object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  ThreeStateUnit_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateUnit *]
 *
 */
ThreeStateUnit * ThreeStateUnit_alloc(void) 
{
    ThreeStateUnit * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(ThreeStateUnit *) ckalloc (sizeof(ThreeStateUnit))) == NULL)    {  
      warn("ThreeStateUnit_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    /* match_emission[ALPHABET_SIZE] is an array: no default possible */ 
    /* insert_emission[ALPHABET_SIZE] is an array: no default possible */ 
    /* transition[TRANSITION_LEN] is an array: no default possible */ 
    out->display_char = 'u'; 


    return out;  
}    


/* Function:  free_ThreeStateUnit(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ThreeStateUnit *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateUnit *]
 *
 */
ThreeStateUnit * free_ThreeStateUnit(ThreeStateUnit * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a ThreeStateUnit obj. Should be trappable");    
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


/* Function:  swap_ThreeStateModel(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_ThreeStateModel
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [ThreeStateUnit **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_ThreeStateModel(ThreeStateUnit ** list,int i,int j)  
{
    ThreeStateUnit * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_ThreeStateModel(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_ThreeStateModel which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [ThreeStateUnit **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_ThreeStateModel(ThreeStateUnit ** list,int left,int right,int (*comp)(ThreeStateUnit * ,ThreeStateUnit * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_ThreeStateModel(list,left,(left+right)/2);  
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_ThreeStateModel (list,++last,i);    
      }  
    swap_ThreeStateModel (list,left,last);   
    qsort_ThreeStateModel(list,left,last-1,comp);    
    qsort_ThreeStateModel(list,last+1,right,comp);   
}    


/* Function:  sort_ThreeStateModel(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_ThreeStateModel
 *
 *
 * Arg:         obj [UNKN ] Object containing list [ThreeStateModel *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_ThreeStateModel(ThreeStateModel * obj,int (*comp)(ThreeStateUnit *, ThreeStateUnit *)) 
{
    qsort_ThreeStateModel(obj->unit,0,obj->len-1,comp);  
    return;  
}    


/* Function:  expand_ThreeStateModel(obj,len)
 *
 * Descrip:    Really an internal function for add_ThreeStateModel
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [ThreeStateModel *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_ThreeStateModel(ThreeStateModel * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_ThreeStateModel called with no need");    
      return TRUE;   
      }  


    if( (obj->unit = (ThreeStateUnit ** ) ckrealloc (obj->unit,sizeof(ThreeStateUnit *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_ThreeStateModel, returning FALSE");  
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_ThreeStateModel(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [ThreeStateModel *]
 * Arg:        add [OWNER] Object to add to the list [ThreeStateUnit *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_ThreeStateModel(ThreeStateModel * obj,ThreeStateUnit * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_ThreeStateModel(obj,obj->len + ThreeStateModelLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->unit[obj->len++]=add;   
    return TRUE; 
}    


/* Function:  flush_ThreeStateModel(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [ThreeStateModel *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_ThreeStateModel(ThreeStateModel * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->unit[i] != NULL)  {  
        free_ThreeStateUnit(obj->unit[i]);   
        obj->unit[i] = NULL; 
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  ThreeStateModel_alloc_std(void)
 *
 * Descrip:    Equivalent to ThreeStateModel_alloc_len(ThreeStateModelLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateModel *]
 *
 */
ThreeStateModel * ThreeStateModel_alloc_std(void) 
{
    return ThreeStateModel_alloc_len(ThreeStateModelLISTLENGTH); 
}    


/* Function:  ThreeStateModel_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateModel *]
 *
 */
ThreeStateModel * ThreeStateModel_alloc_len(int len) 
{
    ThreeStateModel * out;  /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = ThreeStateModel_alloc()) == NULL)  
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->unit = (ThreeStateUnit ** ) ckcalloc (len,sizeof(ThreeStateUnit *))) == NULL)   {  
      warn("Warning, ckcalloc failed in ThreeStateModel_alloc_len"); 
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_ThreeStateModel(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ThreeStateModel *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateModel *]
 *
 */
ThreeStateModel * hard_link_ThreeStateModel(ThreeStateModel * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a ThreeStateModel object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  ThreeStateModel_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateModel *]
 *
 */
ThreeStateModel * ThreeStateModel_alloc(void) 
{
    ThreeStateModel * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(ThreeStateModel *) ckalloc (sizeof(ThreeStateModel))) == NULL)  {  
      warn("ThreeStateModel_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->name = NULL;    
    out->unit = NULL;    
    out->len = out->maxlen = 0;  
    out->alphabet = NULL;    
    out->accession = NULL;   
    out->threshold = 0;  
    out->rm = NULL;  


    return out;  
}    


/* Function:  free_ThreeStateModel(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ThreeStateModel *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateModel *]
 *
 */
ThreeStateModel * free_ThreeStateModel(ThreeStateModel * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a ThreeStateModel obj. Should be trappable");   
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
    if( obj->name != NULL)   
      ckfree(obj->name);     
    if( obj->unit != NULL)   {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->unit[i] != NULL)    
          free_ThreeStateUnit(obj->unit[i]); 
        }  
      ckfree(obj->unit); 
      }  
    if( obj->alphabet != NULL)   
      ckfree(obj->alphabet);     
    if( obj->accession != NULL)  
      ckfree(obj->accession);    
    if( obj->rm != NULL) 
      free_RandomModel(obj->rm);     


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_ThreeStateScoreUnit(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ThreeStateScoreUnit *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateScoreUnit *]
 *
 */
ThreeStateScoreUnit * hard_link_ThreeStateScoreUnit(ThreeStateScoreUnit * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a ThreeStateScoreUnit object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  ThreeStateScoreUnit_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateScoreUnit *]
 *
 */
ThreeStateScoreUnit * ThreeStateScoreUnit_alloc(void) 
{
    ThreeStateScoreUnit * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(ThreeStateScoreUnit *) ckalloc (sizeof(ThreeStateScoreUnit))) == NULL)  {  
      warn("ThreeStateScoreUnit_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    /* match[ALPHABET_SIZE] is an array: no default possible */ 
    /* insert[ALPHABET_SIZE] is an array: no default possible */ 
    /* trans[TRANSITION_LEN] is an array: no default possible */ 


    return out;  
}    


/* Function:  free_ThreeStateScoreUnit(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ThreeStateScoreUnit *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateScoreUnit *]
 *
 */
ThreeStateScoreUnit * free_ThreeStateScoreUnit(ThreeStateScoreUnit * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a ThreeStateScoreUnit obj. Should be trappable");   
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


/* Function:  swap_ThreeStateScore(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_ThreeStateScore
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [ThreeStateScoreUnit **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_ThreeStateScore(ThreeStateScoreUnit ** list,int i,int j)  
{
    ThreeStateScoreUnit * temp;  
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_ThreeStateScore(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_ThreeStateScore which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [ThreeStateScoreUnit **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_ThreeStateScore(ThreeStateScoreUnit ** list,int left,int right,int (*comp)(ThreeStateScoreUnit * ,ThreeStateScoreUnit * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_ThreeStateScore(list,left,(left+right)/2);  
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_ThreeStateScore (list,++last,i);    
      }  
    swap_ThreeStateScore (list,left,last);   
    qsort_ThreeStateScore(list,left,last-1,comp);    
    qsort_ThreeStateScore(list,last+1,right,comp);   
}    


/* Function:  sort_ThreeStateScore(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_ThreeStateScore
 *
 *
 * Arg:         obj [UNKN ] Object containing list [ThreeStateScore *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_ThreeStateScore(ThreeStateScore * obj,int (*comp)(ThreeStateScoreUnit *, ThreeStateScoreUnit *)) 
{
    qsort_ThreeStateScore(obj->unit,0,obj->len-1,comp);  
    return;  
}    


/* Function:  expand_ThreeStateScore(obj,len)
 *
 * Descrip:    Really an internal function for add_ThreeStateScore
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [ThreeStateScore *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_ThreeStateScore(ThreeStateScore * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_ThreeStateScore called with no need");    
      return TRUE;   
      }  


    if( (obj->unit = (ThreeStateScoreUnit ** ) ckrealloc (obj->unit,sizeof(ThreeStateScoreUnit *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_ThreeStateScore, returning FALSE");  
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_ThreeStateScore(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [ThreeStateScore *]
 * Arg:        add [OWNER] Object to add to the list [ThreeStateScoreUnit *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_ThreeStateScore(ThreeStateScore * obj,ThreeStateScoreUnit * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_ThreeStateScore(obj,obj->len + ThreeStateScoreLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->unit[obj->len++]=add;   
    return TRUE; 
}    


/* Function:  flush_ThreeStateScore(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [ThreeStateScore *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_ThreeStateScore(ThreeStateScore * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->unit[i] != NULL)  {  
        free_ThreeStateScoreUnit(obj->unit[i]);  
        obj->unit[i] = NULL; 
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  ThreeStateScore_alloc_std(void)
 *
 * Descrip:    Equivalent to ThreeStateScore_alloc_len(ThreeStateScoreLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateScore *]
 *
 */
ThreeStateScore * ThreeStateScore_alloc_std(void) 
{
    return ThreeStateScore_alloc_len(ThreeStateScoreLISTLENGTH); 
}    


/* Function:  ThreeStateScore_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateScore *]
 *
 */
ThreeStateScore * ThreeStateScore_alloc_len(int len) 
{
    ThreeStateScore * out;  /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = ThreeStateScore_alloc()) == NULL)  
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->unit = (ThreeStateScoreUnit ** ) ckcalloc (len,sizeof(ThreeStateScoreUnit *))) == NULL) {  
      warn("Warning, ckcalloc failed in ThreeStateScore_alloc_len"); 
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_ThreeStateScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ThreeStateScore *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateScore *]
 *
 */
ThreeStateScore * hard_link_ThreeStateScore(ThreeStateScore * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a ThreeStateScore object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  ThreeStateScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateScore *]
 *
 */
ThreeStateScore * ThreeStateScore_alloc(void) 
{
    ThreeStateScore * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(ThreeStateScore *) ckalloc (sizeof(ThreeStateScore))) == NULL)  {  
      warn("ThreeStateScore_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->unit = NULL;    
    out->len = out->maxlen = 0;  
    out->name = NULL;    
    out->accession = NULL;   


    return out;  
}    


/* Function:  free_ThreeStateScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ThreeStateScore *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateScore *]
 *
 */
ThreeStateScore * free_ThreeStateScore(ThreeStateScore * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a ThreeStateScore obj. Should be trappable");   
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
    if( obj->unit != NULL)   {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->unit[i] != NULL)    
          free_ThreeStateScoreUnit(obj->unit[i]);    
        }  
      ckfree(obj->unit); 
      }  
    if( obj->name != NULL)   
      ckfree(obj->name);     
    if( obj->accession != NULL)  
      ckfree(obj->accession);    


    ckfree(obj); 
    return NULL; 
}    


/* Function:  replace_name_ThreeStateModel(obj,name)
 *
 * Descrip:    Replace member variable name
 *             For use principly by API functions
 *
 *
 * Arg:         obj [UNKN ] Object holding the variable [ThreeStateModel *]
 * Arg:        name [OWNER] New value of the variable [char *]
 *
 * Return [SOFT ]  member variable name [boolean]
 *
 */
boolean replace_name_ThreeStateModel(ThreeStateModel * obj,char * name) 
{
    if( obj == NULL)     {  
      warn("In replacement function name for object ThreeStateModel, got a NULL object");    
      return FALSE;  
      }  
    obj->name = name;    
    return TRUE; 
}    


/* Function:  access_name_ThreeStateModel(obj)
 *
 * Descrip:    Access member variable name
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [ThreeStateModel *]
 *
 * Return [SOFT ]  member variable name [char *]
 *
 */
char * access_name_ThreeStateModel(ThreeStateModel * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function name for object ThreeStateModel, got a NULL object");   
      return NULL;   
      }  
    return obj->name;    
}    


/* Function:  access_unit_ThreeStateModel(obj,i)
 *
 * Descrip:    Access members stored in the unit list
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the list [ThreeStateModel *]
 * Arg:          i [UNKN ] Position in the list [int]
 *
 * Return [SOFT ]  Element of the list [ThreeStateUnit *]
 *
 */
ThreeStateUnit * access_unit_ThreeStateModel(ThreeStateModel * obj,int i) 
{
    if( obj == NULL)     {  
      warn("In accessor function unit for object ThreeStateModel, got a NULL object");   
      return NULL;   
      }  
    if( obj->len <= i )  {  
      warn("In accessor function unit for object ThreeStateModel, index %%d is greater than list length %%d",i,obj->len);    
      return NULL;   
      }  
    return obj->unit[i];     
}    


/* Function:  length_unit_ThreeStateModel(obj)
 *
 * Descrip:    discover the length of the list
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the list [ThreeStateModel *]
 *
 * Return [UNKN ]  length of the list [int]
 *
 */
int length_unit_ThreeStateModel(ThreeStateModel * obj) 
{
    if( obj == NULL)     {  
      warn("In length function unit for object ThreeStateModel, got a NULL object"); 
      return -1;     
      }  
    return obj->len;     
}    


/* Function:  replace_alphabet_ThreeStateModel(obj,alphabet)
 *
 * Descrip:    Replace member variable alphabet
 *             For use principly by API functions
 *
 *
 * Arg:             obj [UNKN ] Object holding the variable [ThreeStateModel *]
 * Arg:        alphabet [OWNER] New value of the variable [char *]
 *
 * Return [SOFT ]  member variable alphabet [boolean]
 *
 */
boolean replace_alphabet_ThreeStateModel(ThreeStateModel * obj,char * alphabet) 
{
    if( obj == NULL)     {  
      warn("In replacement function alphabet for object ThreeStateModel, got a NULL object");    
      return FALSE;  
      }  
    obj->alphabet = alphabet;    
    return TRUE; 
}    


/* Function:  access_alphabet_ThreeStateModel(obj)
 *
 * Descrip:    Access member variable alphabet
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [ThreeStateModel *]
 *
 * Return [SOFT ]  member variable alphabet [char *]
 *
 */
char * access_alphabet_ThreeStateModel(ThreeStateModel * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function alphabet for object ThreeStateModel, got a NULL object");   
      return NULL;   
      }  
    return obj->alphabet;    
}    


/* Function:  replace_accession_ThreeStateModel(obj,accession)
 *
 * Descrip:    Replace member variable accession
 *             For use principly by API functions
 *
 *
 * Arg:              obj [UNKN ] Object holding the variable [ThreeStateModel *]
 * Arg:        accession [OWNER] New value of the variable [char *]
 *
 * Return [SOFT ]  member variable accession [boolean]
 *
 */
boolean replace_accession_ThreeStateModel(ThreeStateModel * obj,char * accession) 
{
    if( obj == NULL)     {  
      warn("In replacement function accession for object ThreeStateModel, got a NULL object");   
      return FALSE;  
      }  
    obj->accession = accession;  
    return TRUE; 
}    


/* Function:  access_accession_ThreeStateModel(obj)
 *
 * Descrip:    Access member variable accession
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [ThreeStateModel *]
 *
 * Return [SOFT ]  member variable accession [char *]
 *
 */
char * access_accession_ThreeStateModel(ThreeStateModel * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function accession for object ThreeStateModel, got a NULL object");  
      return NULL;   
      }  
    return obj->accession;   
}    


/* Function:  replace_threshold_ThreeStateModel(obj,threshold)
 *
 * Descrip:    Replace member variable threshold
 *             For use principly by API functions
 *
 *
 * Arg:              obj [UNKN ] Object holding the variable [ThreeStateModel *]
 * Arg:        threshold [OWNER] New value of the variable [double]
 *
 * Return [SOFT ]  member variable threshold [boolean]
 *
 */
boolean replace_threshold_ThreeStateModel(ThreeStateModel * obj,double threshold) 
{
    if( obj == NULL)     {  
      warn("In replacement function threshold for object ThreeStateModel, got a NULL object");   
      return FALSE;  
      }  
    obj->threshold = threshold;  
    return TRUE; 
}    


/* Function:  access_threshold_ThreeStateModel(obj)
 *
 * Descrip:    Access member variable threshold
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [ThreeStateModel *]
 *
 * Return [SOFT ]  member variable threshold [double]
 *
 */
double access_threshold_ThreeStateModel(ThreeStateModel * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function threshold for object ThreeStateModel, got a NULL object");  
      return 0;  
      }  
    return obj->threshold;   
}    


/* Function:  replace_rm_ThreeStateModel(obj,rm)
 *
 * Descrip:    Replace member variable rm
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [ThreeStateModel *]
 * Arg:         rm [OWNER] New value of the variable [RandomModel *]
 *
 * Return [SOFT ]  member variable rm [boolean]
 *
 */
boolean replace_rm_ThreeStateModel(ThreeStateModel * obj,RandomModel * rm) 
{
    if( obj == NULL)     {  
      warn("In replacement function rm for object ThreeStateModel, got a NULL object");  
      return FALSE;  
      }  
    obj->rm = rm;    
    return TRUE; 
}    


/* Function:  access_rm_ThreeStateModel(obj)
 *
 * Descrip:    Access member variable rm
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [ThreeStateModel *]
 *
 * Return [SOFT ]  member variable rm [RandomModel *]
 *
 */
RandomModel * access_rm_ThreeStateModel(ThreeStateModel * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function rm for object ThreeStateModel, got a NULL object"); 
      return NULL;   
      }  
    return obj->rm;  
}    


/* Function:  replace_display_char_ThreeStateUnit(obj,display_char)
 *
 * Descrip:    Replace member variable display_char
 *             For use principly by API functions
 *
 *
 * Arg:                 obj [UNKN ] Object holding the variable [ThreeStateUnit *]
 * Arg:        display_char [OWNER] New value of the variable [char]
 *
 * Return [SOFT ]  member variable display_char [boolean]
 *
 */
boolean replace_display_char_ThreeStateUnit(ThreeStateUnit * obj,char display_char) 
{
    if( obj == NULL)     {  
      warn("In replacement function display_char for object ThreeStateUnit, got a NULL object"); 
      return FALSE;  
      }  
    obj->display_char = display_char;    
    return TRUE; 
}    


/* Function:  access_display_char_ThreeStateUnit(obj)
 *
 * Descrip:    Access member variable display_char
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [ThreeStateUnit *]
 *
 * Return [SOFT ]  member variable display_char [char]
 *
 */
char access_display_char_ThreeStateUnit(ThreeStateUnit * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function display_char for object ThreeStateUnit, got a NULL object");    
      return 'u';    
      }  
    return obj->display_char;    
}    



#ifdef _cplusplus
}
#endif
