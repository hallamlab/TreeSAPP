#ifdef _cplusplus
extern "C" {
#endif
#include "estwrap.h"


/* Function:  string_from_alg_estwrap(alg_type)
 *
 * Descrip:    Returns the string form of the algorithm
 *
 *
 * Arg:        alg_type [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [char *]
 *
 */
# line 36 "estwrap.dy"
char * string_from_alg_estwrap(int alg_type)
{
  switch(alg_type) {
  case ESTWISE_3 : return "333";
  case ESTSLIM_3 : return "312";
  case ESTLOOP_3 : return "333L";
  case ESTQUICK_3 : return "312Q";
  default : return "No algorithm type specd";
  }
}

/* Function:  alg_estwrap_from_string(str)
 *
 * Descrip:    This function returns the algorithm type
 *             for an est search from the string
 *
 *
 * Arg:        str [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 51 "estwrap.dy"
int alg_estwrap_from_string(char * str)
{
  int t;

  t = get_number_from_slashed_string(str,"333/333L/333F/312/312L/312Q");

  switch (t) {
  case 0 : return ESTWISE_3;
  case 1 : return ESTLOOP_3;
  case 2 : return ESTFRAG_3;
  case 3 : return ESTSLIM_3;
  case 4 : return ESTSLIM_L;
  case 5 : return ESTQUICK_3;
  default : warn("Cannot convert string %s into a valid estwise algorithm type\n",str);
    return -1;
  }
}


/* Function:  AlnBlock_from_Protein_estwise_wrap(pro,cdna,cp,cm,ct,comp,gap,ext,is_global,rmd,alg,rm,use_syn,allN,dpri,palpoi)
 *
 * Descrip:    This function is the guts for the est single alignment
 *             mode.
 *
 *             It uses /AlnBlock_from_TSM_estwise_wrap for the
 *             heavy part of the call
 *
 *
 * Arg:              pro [READ ] protein to be used in the comparison [Protein *]
 * Arg:             cdna [READ ] cdna to be compared to [cDNA *]
 * Arg:               cp [READ ] cdna parser indicating insertion deletion probabilities [cDNAParser *]
 * Arg:               cm [READ ] codon mapper indicating substitution errors etc [CodonMapper *]
 * Arg:               ct [READ ] codon table for the codon->amino acid mappings [CodonTable *]
 * Arg:             comp [READ ] comparison matrix to use [CompMat *]
 * Arg:              gap [UNKN ] gap penalty [int]
 * Arg:              ext [UNKN ] extension penalty [int]
 * Arg:        is_global [UNKN ] if true, global start-end in protein is used [boolean]
 * Arg:              rmd [UNKN ] random model of dna to use [RandomModelDNA *]
 * Arg:              alg [UNKN ] est algorithm type to use [int]
 * Arg:               rm [UNKN ] random protein model for use with compmat [RandomModel *]
 * Arg:          use_syn [UNKN ] if true, uses a synchronous coding model [boolean]
 * Arg:             allN [UNKN ] Undocumented argument [Probability]
 * Arg:             dpri [UNKN ] Undocumented argument [DPRunImpl *]
 * Arg:           palpoi [WRITE] the raw packed alignment output if wanted [PackAln **]
 *
 * Return [UNKN ]  Undocumented return value [AlnBlock *]
 *
 */
# line 92 "estwrap.dy"
AlnBlock * AlnBlock_from_Protein_estwise_wrap(Protein * pro,cDNA * cdna,cDNAParser * cp,CodonMapper * cm,CodonTable * ct,CompMat * comp,int gap,int ext,boolean is_global,RandomModelDNA * rmd,int alg,RandomModel * rm,boolean use_syn,Probability allN,DPRunImpl * dpri,PackAln ** palpoi)
{
  ThreeStateModel * tsm;
  AlnBlock * out;

  if( pro == NULL || cdna == NULL || comp == NULL || rm == NULL || dpri == NULL){
    warn("trappable error in PackAln from protein sequence vs cDNA, passed some NULL objects, Complain!");
    return NULL;
  }

  rm = default_RandomModel();
  
  if( is_global == TRUE) 
    tsm = global_ThreeStateModel_from_half_bit_Sequence(pro,comp,rm,gap,ext);
  else
    tsm = ThreeStateModel_from_half_bit_Sequence(pro,comp,rm,gap,ext);
  
  out = AlnBlock_from_TSM_estwise_wrap(tsm,cdna,cp,cm,ct,rmd,alg,use_syn,FALSE,allN,dpri,palpoi);

  free_ThreeStateModel(tsm);
  free_RandomModel(rm);

  return out;

}


  


/* Function:  AlnBlock_from_TSM_estwise_wrap(tsm,cdna,cp,cm,ct,rmd,alg,use_syn,force_flat_insert,allN,dpri,palpoi)
 *
 * Descrip:    This function is the basic wrap for Protein models
 *             vs cDNA sequences.
 *
 *
 * Arg:                      tsm [READ ] threestatemodel to be compared to the dna [ThreeStateModel *]
 * Arg:                     cdna [READ ] cdna to be compared to [cDNA *]
 * Arg:                       cp [READ ] cdna parser indicating insertion deletion probabilities [cDNAParser *]
 * Arg:                       cm [READ ] codon mapper indicating substitution errors etc [CodonMapper *]
 * Arg:                       ct [READ ] codon table for the codon->amino acid mappings [CodonTable *]
 * Arg:                      rmd [UNKN ] random model of dna to use [RandomModelDNA *]
 * Arg:                      alg [UNKN ] est algorithm type to use [int]
 * Arg:                  use_syn [UNKN ] if true, uses a synchronous coding model [boolean]
 * Arg:        force_flat_insert [UNKN ] Undocumented argument [boolean]
 * Arg:                     allN [UNKN ] Undocumented argument [Probability]
 * Arg:                     dpri [UNKN ] Undocumented argument [DPRunImpl *]
 * Arg:                   palpoi [WRITE] the raw packed alignment output if wanted [PackAln **]
 *
 * Return [UNKN ]  Undocumented return value [AlnBlock *]
 *
 */
# line 136 "estwrap.dy"
AlnBlock * AlnBlock_from_TSM_estwise_wrap(ThreeStateModel * tsm,cDNA * cdna,cDNAParser * cp,CodonMapper * cm,CodonTable * ct,RandomModelDNA * rmd,int alg,boolean use_syn,boolean force_flat_insert,Probability allN,DPRunImpl * dpri,PackAln ** palpoi)
{
  AlnBlock * out;
  PackAln * pal;

  cDNAParserScore * cps = NULL ;
  GeneWise * gw = NULL ;
  GeneWiseScore * gws = NULL ;

  ComplexSequence * cs = NULL ;
  ComplexSequenceEvalSet * cses = NULL ;

  

  if( tsm == NULL || cdna == NULL || cp == NULL || rmd == NULL || dpri == NULL){
    warn("trappable error in AlnBlock estwise wrap, passed some NULL objects, Complain!");
    return NULL;
  }

  if( (gw=GeneWise_from_ThreeStateModel_cdna(tsm,cp,cm,allN)) == NULL) {
    warn("Unable to make GeneWise model in estwise wrap");
    goto exit;
  }

  if( use_syn == TRUE ) {
    if( tsm->rm == NULL ) {
      warn("A three state model with no random model! Ugh!");
      goto exit;
    }

    GeneWise_fold_in_synchronised_RandomModel(gw,tsm->rm,cm,ct,0.5);
    if( force_flat_insert == TRUE ) {
      check_flat_insert(gw,TRUE,FALSE,ct);
    }
  } else {
    GeneWise_fold_in_RandomModelDNA(gw,rmd);
  }

  if( (gws = GeneWiseScore_from_GeneWise(gw)) == NULL) {
    warn("Unable to make GeneWiseScore model in estwise wrap");
    goto exit;
  }


  if( (cps = cDNAParserScore_from_cDNAParser(cp)) == NULL ) {
    warn("Unable to make cDNAParserScore in estwise wrap");
    goto exit;
  }

  cses = default_cDNA_ComplexSequenceEvalSet();

  cs = new_ComplexSequence(cdna->baseseq,cses);
  
  cses = free_ComplexSequenceEvalSet(cses);

  switch(alg) {
  case ESTWISE_3 :
    pal = PackAln_bestmemory_EstWise3(gws,cs,cps,NULL,dpri);
    out = convert_PackAln_to_AlnBlock_EstWise3(pal);
    break;
  case ESTLOOP_3 :
    pal = PackAln_bestmemory_EstLoop3(gws,cs,cps,NULL,dpri);
    out = convert_PackAln_to_AlnBlock_EstLoop3(pal);
    break;
  case ESTSLIM_3 :
    pal = PackAln_bestmemory_EstSlim3(gws,cs,cps,NULL,dpri);
    out = convert_PackAln_to_AlnBlock_EstLoop3(pal);
    break;
  case ESTSLIM_L :
    pal = PackAln_bestmemory_EstSlimLoop3(gws,cs,cps,NULL,dpri);
    out = convert_PackAln_to_AlnBlock_EstSlimLoop3(pal);
    break;
  case ESTFRAG_3 :
    pal = PackAln_bestmemory_EstFrag3(gws,cs,cps,0,0,NULL,dpri);
    out = convert_PackAln_to_AlnBlock_EstFrag3(pal);
    break;
  default :
    warn("No algorithm type specified. Not good news!");
    goto exit;
  }


  if( palpoi != NULL) {
    *palpoi = pal;
  } else {
    free_PackAln(pal);
  }
  
  goto exit;


  exit :
    
  if( cps != NULL ) 
    free_cDNAParserScore(cps);
  
  if( gw != NULL )
    free_GeneWise(gw);
  
  if( gws != NULL)
    free_GeneWiseScore(gws);
  
  if( cs != NULL)
    free_ComplexSequence(cs);

  if( cses != NULL)
    free_ComplexSequenceEvalSet(cses);
  
  return out;
}


/* Function:  Hscore_from_TSM_estwise(tdb,cdb,cp,cm,rmd,use_syn,alg,bits_cutoff,allN,flat_insert,report_level,die_on_error,dbsi)
 *
 * Descrip:    Runs a database search for the estwise set
 *             of algorithms
 *
 *
 * Arg:                 tdb [READ ] a three state model database [ThreeStateDB *]
 * Arg:                 cdb [READ ] a dna sequence database [cDNADB *]
 * Arg:                  cp [READ ] the codon parser for this comparison [cDNAParser *]
 * Arg:                  cm [READ ] the codon mapper for this comparison [CodonMapper *]
 * Arg:                 rmd [READ ] random model used for the dna sequence comparison [RandomModelDNA *]
 * Arg:             use_syn [UNKN ] whether a synchronous coding model should be used or not [boolean]
 * Arg:                 alg [UNKN ] algorithm to use [int]
 * Arg:         bits_cutoff [UNKN ] Undocumented argument [double]
 * Arg:                allN [UNKN ] Undocumented argument [Probability]
 * Arg:         flat_insert [UNKN ] Undocumented argument [boolean]
 * Arg:        report_level [UNKN ] Undocumented argument [int]
 * Arg:        die_on_error [UNKN ] if true, dies if there is an error [boolean]
 * Arg:                dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 *
 * Return [OWNER]  a newly allocated Hscore structure of the search [Hscore *]
 *
 */
# line 262 "estwrap.dy"
Hscore * Hscore_from_TSM_estwise(ThreeStateDB * tdb,cDNADB * cdb,cDNAParser * cp,CodonMapper * cm,RandomModelDNA * rmd,boolean use_syn,int alg,double bits_cutoff,Probability allN,boolean flat_insert,int report_level,boolean die_on_error,DBSearchImpl * dbsi)
{
  Hscore * out = NULL;
  GeneWiseDB * gwdb;
  cDNAParserScore * cps = NULL ;
  GeneWiseQuickDB * gwq;

  Search_Return_Type ret;
  
  ret = SEARCH_ERROR;

  if( alg == ESTQUICK_3 && tdb->type != TSMDB_SINGLE ) {
    warn("Can only currently use estquick in a single mode search");
    return NULL;
  }


  gwdb = new_GeneWiseDB_cdna(tdb,cp,cm,rmd,use_syn,allN,flat_insert);
  if( gwdb == NULL ) {
    warn("Could not build a new GeneWiseDB from the objects provided. Exiting without completing the search");
    goto exit;
  }


  if( (cps = cDNAParserScore_from_cDNAParser(cp)) == NULL ) {
    warn("Unable to make cDNAParserScore in estwise wrap");
    goto exit;
  }

  /*** allocate Hscore structure ***/

  out = std_bits_Hscore(bits_cutoff,report_level);

  switch(alg) {

  case ESTWISE_3 :
    ret = search_EstWise3(dbsi,out,gwdb,cdb,cps);
    break;
  case ESTSLIM_3 :
    ret = search_EstSlim3(dbsi,out,gwdb,cdb,cps);
    break;

  case ESTQUICK_3 :
    gwq = GeneWiseQuickDB_from_GeneWiseDB(gwdb);
    ret = search_EstQuick3(dbsi,out,gwq,cdb,cps);
    free_GeneWiseQuickDB(gwq);
    break;
  default :
    warn("A major problem. No valid algorithm type passed in");
    goto exit;
  }

  goto exit;


  exit :
  
    if( gwdb != NULL ) {
      free_GeneWiseDB(gwdb);
    }
    
  if( cps != NULL ) 
    free_cDNAParserScore(cps);

   

  return out;
}

/* Function:  write_mul_estwise_AlnBlock(alb,ct,ofp)
 *
 * Descrip:    writes an mul format protein multiple
 *             alignment from an AlnBlock with the first
 *             sequence an HMM/protein and ignored, and
 *             the second and subsequent sequences cdna sequences
 *             which are then translated into proteins
 *
 *             This relies considerably on the alb being made
 *             correctly, and if it is not, then god help you.
 *
 *             the estwisedb programs makes the alb correctly
 *
 *
 * Arg:        alb [UNKN ] Undocumented argument [AlnBlock *]
 * Arg:         ct [UNKN ] Undocumented argument [CodonTable *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 343 "estwrap.dy"
void write_mul_estwise_AlnBlock(AlnBlock * alb,CodonTable * ct,FILE * ofp)
{
  char namebuffer[128];
  AlnSequence * als;
  AlnUnit *ale;
  AlnColumn * alc;
  Sequence * seq;
  int i;


  assert(alb);
  assert(ct);

  for(i=1;i<alb->len;i++) {
    als = alb->seq[i];
    if( als->data == NULL ) {
      warn("For sequence %d in the estwise alnblock, no attached sequence, and so cannot write. Skipping",i);
      continue;
    }

    seq = (Sequence *) als->data; /* scared? I am! */
    for(alc = alb->start,ale=NULL;alc->next != NULL;alc = alc->next)
      if( strstr(alc->alu[i]->text_label,"CODON") != NULL ) {
	ale = alc->alu[i];
      }

    if( ale == NULL ) {
      warn("Unable to find even a codon matching this. Exiting for sequence %s in mul output",seq->name);
      continue;
    }

    
    /* fprintf(stdout,"Ale is %d-%d %s\n",ale->start,ale->end,ale->text_label);*/

    if( is_reversed_Sequence(seq) )
      sprintf(namebuffer,"%s/%d-%d",seq->name,ale->end,als->start->start+1);
    else
      sprintf(namebuffer,"%s/%d-%d",seq->name,als->start->start+2,ale->end+1);

    fprintf(ofp,"%-30s ",namebuffer);
    for(alc = alb->start;alc != NULL;alc = alc->next ) {
      if( strstr(alc->alu[i]->text_label,"CODON") != NULL ) {
	fputc(aminoacid_from_seq(ct,seq->seq+alc->alu[i]->start+1),ofp);
      } else if( strstr(alc->alu[i]->text_label,"INSERT") != NULL ) {
	fputc('-',ofp);
      } else {
	fputc('X',ofp);
      }
    }


    fputc('\n',ofp);
  }
}

# line 410 "estwrap.c"

#ifdef _cplusplus
}
#endif
