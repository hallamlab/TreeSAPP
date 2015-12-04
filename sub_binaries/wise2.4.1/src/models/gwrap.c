#ifdef _cplusplus
extern "C" {
#endif
#include "gwrap.h"
    
  
/* Function:  GeneParameter21_wrap(gf,subs_error,indel_error,rmd,use_modelled_codon,use_modelled_splice,tie_intron_prob,ct,rnd_loop,cds_loop,rnd_to_model,link_loop,link_to_model)
 *
 * Descrip:    A general wrap over the production of parameters for the
 *             GeneWise programs. The geneparameter21 holds all the parameters,
 *             and can be approximated for the 6:23 and 4:21 algorithms
 *
 *             This function is the best way to make a GeneParameter21 object
 *             as all the different options for how to make it or modify its
 *             contents are laid out as arguments to this function
 *
 *
 * Arg:                         gf [READ ] Gene Frequency data structure, holding counts for splice sites etc [GeneFrequency21 *]
 * Arg:                 subs_error [UNKN ] substitution error on the dna sequence [double]
 * Arg:                indel_error [UNKN ] rough estimate of the insertion/deletion per base error rate [double]
 * Arg:                        rmd [UNKN ] the random model of the DNA that is used [RandomModelDNA *]
 * Arg:         use_modelled_codon [UNKN ] if TRUE, model codon frequency [boolean]
 * Arg:        use_modelled_splice [UNKN ] if TRUE, make splice models from gf parameters [boolean]
 * Arg:            tie_intron_prob [UNKN ] Undocumented argument [boolean]
 * Arg:                         ct [UNKN ] codon table which is used for codon->aa mapping [CodonTable *]
 * Arg:                   rnd_loop [UNKN ] Undocumented argument [Probability]
 * Arg:                   cds_loop [UNKN ] Undocumented argument [Probability]
 * Arg:               rnd_to_model [UNKN ] Undocumented argument [Probability]
 * Arg:                  link_loop [UNKN ] Undocumented argument [Probability]
 * Arg:              link_to_model [UNKN ] Undocumented argument [Probability]
 *
 * Return [UNKN ]  A newly allocated structure [GeneParameter21 *]
 *
 */
# line 71 "gwrap.dy"
GeneParameter21 * GeneParameter21_wrap(GeneFrequency21 * gf,double subs_error,double indel_error,RandomModelDNA * rmd,boolean use_modelled_codon,boolean use_modelled_splice,boolean tie_intron_prob,CodonTable * ct,Probability rnd_loop,Probability cds_loop,Probability rnd_to_model,Probability link_loop,Probability link_to_model)
{
  GeneParameter21 * out;
  int i;

  out = GeneParameter21_from_GeneFrequency21(gf,ct,rmd,rnd_loop,cds_loop,rnd_to_model,link_loop,link_to_model);


  if( use_modelled_codon == FALSE ) {
    out->cm = free_CodonMapper(out->cm);
    out->cm = flat_CodonMapper(ct);
  }

  if( use_modelled_splice == FALSE )  {
    out->cses = free_ComplexSequenceEvalSet(out->cses);
    out->cses = default_genomic_ComplexSequenceEvalSet();
    out->modelled_splice = FALSE;
  }

  if( tie_intron_prob == TRUE ) {
    for(i=0;i<5;i++) 
      out->gp->central[i] = rmd->base[i];
  }

  /*** errors ***/

  sprinkle_errors_over_CodonMapper(out->cm,subs_error);

  add_flat_error_probabilities_GeneParser21(out->gp,indel_error);

  GeneParser21_fold_in_RandomModelDNA(out->gp,rmd);

  fold_in_RandomModelDNA_into_RandomCodon(out->rc,rmd);

  return out;
}

/* Function:  AlnBlock_from_protein_genewise_wrap(protein,pg,is_global,dna,comp,gap,ext,gpara,rmd,intergenic,alg,use_syn,rm,allN,startendmode,dpri,pal,gwp)
 *
 * Descrip:    A function which aligns a Protein sequecne to a Genomic sequence
 *             under the Comparison matrix comp and the gene paras in gpara.
 *
 *             This is the best function for accessing GeneWise functionality
 *             for a protein to dna comparison, allowing for introns.
 *
 *             To make the protein object, you will first read in a generic
 *             sequence object using something like read_fasta_Sequence and
 *             then convert it to a protein object using new_Protein_from_Sequence
 *
 *             To make the genomic object, you will first read in a generic
 *             sequence object using something like read_fasta_Sequence and
 *             then convert it to a genomic object using new_Genomic_from_Sequence
 *
 *             To make a CompMat object you will use read_Blast_file_CompMat
 *             from the compmat module. It is likely, if the Wise2 enviroment
 *             has been set up correctly that read_Blast_file_CompMat("blosum62.bla")
 *             will be fine. You should at the moment only use halfbit matrices
 *             (blosum62 is one such matrix)
 *
 *             To make the necessary random modules use the default construtors
 *             in the randommodel module
 *
 *             To make the gene parameter object use the GeneParameter21_wrap
 *             function found in this module. It will need GeneFrequencies
 *             read in using the read_GeneFrequency21_file function in
 *             the genefrequency module.  Again if Wise2 has been set up
 *             correctly, read_GeneFrequency21_file("human.gf") should work
 *
 *             To again a valid algorithm type use gwrap_alg_type_from_string
 *             found in this module. gwrap_alg_type_from_string("623") would
 *             be a good choice
 *
 *
 *             This function basically makes a threestatemodel (standard HMM) from
 *             the protein and the comparison matrix with the *scary* assumption that
 *             the comparison matrix is in half bit form. It then calls 
 *              /AlnBlock_from_TSM_genewise_wrap to do the nasty stuff. 
 *
 *
 * Arg:             protein [UNKN ] protein sequence used in the comparison [Protein *]
 * Arg:                  pg [READ ] Potential gene - could be NULL - if rough exon positions are known  [NullString]
 * Arg:           is_global [UNKN ] has now become flag for local/global/end-biased switch [NullString]
 * Arg:                 dna [UNKN ] genomic DNA sequence used  [Genomic *]
 * Arg:                comp [UNKN ] protein comparison matrix *in half bits* [CompMat *]
 * Arg:                 gap [UNKN ] gap penalty (negative) [int]
 * Arg:                 ext [UNKN ] extension penalty (negative) [int]
 * Arg:               gpara [UNKN ] Gene parameters. [GeneParameter21 *]
 * Arg:                 rmd [UNKN ] models to be compared to [RandomModelDNA *]
 * Arg:          intergenic [UNKN ] model of random dna between genes [RandomModelDNA *]
 * Arg:                 alg [UNKN ] algorithm type [int]
 * Arg:             use_syn [UNKN ] Undocumented argument [boolean]
 * Arg:                  rm [UNKN ] Undocumented argument [RandomModel *]
 * Arg:                allN [UNKN ] Undocumented argument [Probability]
 * Arg:        startendmode [UNKN ] Undocumented argument [TSM_StartEndMode]
 * Arg:                dpri [UNKN ] Undocumented argument [DPRunImpl *]
 * Arg:                 pal [WRITE] Raw alginment to be saved if non-NULL [PackAln **]
 * Arg:                 gwp [UNKN ] Undocumented argument [GeneWiseRunPara *]
 *
 * Return [UNKN ]  Undocumented return value [AlnBlock *]
 *
 */
# line 161 "gwrap.dy"
AlnBlock * AlnBlock_from_protein_genewise_wrap(Protein * protein,Genomic * dna,CompMat * comp,int gap,int ext,GeneParameter21 * gpara,RandomModelDNA * rmd,RandomModelDNA * intergenic,int alg,boolean use_syn,RandomModel * rm,Probability allN,TSM_StartEndMode startendmode,DPRunImpl * dpri,PackAln ** pal,GeneWiseRunPara * gwp)
{
  ThreeStateModel * tsm;
  RandomModel * rm2;
  AlnBlock * out;
  DPEnvelope * dpenv = NULL;

  if( protein == NULL || dna == NULL || comp == NULL || gpara == NULL || rmd == NULL ){
    warn("trappable error in PackAln from protein sequence, passed some NULL objects, Complain!");
    return NULL;
  }

  assert(dpri);

  rm2 = default_RandomModel();
  
  if( gwp->use_hsp == TRUE ) {
    dpenv = DPEnvelope_from_protein_gen(protein->baseseq,dna->baseseq,comp,gpara->ct,gwp);
  }

  if( startendmode == TSM_default ) {
    startendmode = TSM_endbiased;
  }

  tsm = ThreeStateModel_from_half_bit_Sequence(protein,comp,rm2,gap,ext);
  set_startend_policy_ThreeStateModel(tsm,startendmode,30,halfbit2Probability(-15));
  
  out = AlnBlock_from_TSM_genewise_wrap(tsm,dna,gpara,rmd,intergenic,use_syn,alg,allN,1,dpri,pal,dpenv);

  free_ThreeStateModel(tsm);
  free_RandomModel(rm2);

  return out;

}


/* Function:  AlnBlock_from_TSM_genewise_wrap(tsm,gen,gpara,rmd,intergenic,use_syn,alg,allN,flat_insert,dpri,palpoi,dpenv)
 *
 * Descrip:    A function which aligns a protein HMM (as found
 *             in my threestatemodel structure) to a genomic DNA 
 *             sequence. 
 *
 *             At the moment you are unlikely to be reading in the
 *             HMM structure yourself, so this is not something
 *             you will be doing.
 *
 *             The core algorithms for each method are found in
 *             genewise21/geneloop21 etc files. 
 *
 *
 *
 * Arg:                tsm [UNKN ] protein TSM to be used in the comparison [ThreeStateModel *]
 * Arg:                gen [UNKN ] genomic DNA sequence used  [Genomic *]
 * Arg:              gpara [UNKN ] Gene parameters. [GeneParameter21 *]
 * Arg:                rmd [UNKN ] models to be compared to [RandomModelDNA *]
 * Arg:         intergenic [UNKN ] model of random dna between genes [RandomModelDNA *]
 * Arg:            use_syn [UNKN ] use a synchronous null model [boolean]
 * Arg:                alg [UNKN ] algorithm type [int]
 * Arg:               allN [UNKN ] Undocumented argument [Probability]
 * Arg:        flat_insert [UNKN ] Undocumented argument [boolean]
 * Arg:               dpri [UNKN ] Undocumented argument [DPRunImpl *]
 * Arg:             palpoi [WRITE] Raw alginment to be saved if non-NULL [PackAln **]
 * Arg:              dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [AlnBlock *]
 *
 */
# line 220 "gwrap.dy"
AlnBlock * AlnBlock_from_TSM_genewise_wrap(ThreeStateModel * tsm,Genomic * gen,GeneParameter21 * gpara,RandomModelDNA * rmd,RandomModelDNA * intergenic,boolean use_syn,int alg,Probability allN,boolean flat_insert,DPRunImpl * dpri,PackAln ** palpoi,DPEnvelope * dpenv)
{
  AlnBlock * out = NULL;
  PackAln * pal = NULL;
  ComplexSequence * cs = NULL;
  GeneWise * gw = NULL;
  GeneWiseScore * gws = NULL;
  RandomCodonScore * rcs = NULL ;
  GeneParser21Score  * gps = NULL;
  GeneParser4Score * gp4s = NULL;
  RandomModelDNAScore * ids = NULL;
  GeneralGeneModelScore * ggms = NULL; /* for stretch models */


  Sequence * dna;
  cDNAParserScore * cps = NULL; /* for estwise type algorithms */
  GwLite * gwl = NULL;
  GwLiteScore * gwls = NULL;
  ComplexSequenceEval * tempcse;
  ComplexSequenceEvalSet * cses;
  dna = gen->baseseq;
  

  assert(tsm);
  assert(gen);
  assert(gpara);
  assert(rmd);
  assert(gpara->rc);
  assert(dpri);


  /*show_Genomic(gen,stderr);*/

  /*show_GeneParser21(gpara->gp,stderr); */



  if( tsm == NULL || dna == NULL || gpara == NULL || rmd == NULL){
    warn("trappable error in PackAln from TSM  sequence, passed some NULL objects, Complain!");
    return NULL;
  }

  /*** prepare cses ***/

  if( prepare_ComplexSequenceEvalSet(gpara->cses) == FALSE ) {
    warn("Unable to prepare complexsequenceevalset in TMS2DNA wrap");
    goto exit;
  }

  if( alg == GWWRAP_333 || alg == GWWRAP_333L ) {
    cses = default_cDNA_ComplexSequenceEvalSet();
    cs = new_ComplexSequence(gen->baseseq,cses);
    free_ComplexSequenceEvalSet(cses);
  } else if ( alg == GWWRAP_6LITE ) {
    /* yup. This is scary. */
    tempcse = gpara->cses->cse[1];
    gpara->cses->cse[1] = codon64_number_ComplexSequenceEval();
    cs = new_ComplexSequence(gen->baseseq,gpara->cses);
    free_ComplexSequenceEval(gpara->cses->cse[1]);
    gpara->cses->cse[1] = tempcse;
  } else {
    if( (cs=evaluate_ComplexSequence_Genomic(gen,gpara->cses,0,Probability2Score(0.01))) == FALSE ) {
      warn("Unable to make ComplexSequence in TMS2DNA wrap");
      goto exit;
    }
  }

  /*show_ComplexSequence(cs,stderr);*/


  if( (gw=GeneWise_from_ThreeStateModel(tsm,gpara->gp,gpara->cm,allN,gpara->gwcm)) == NULL) {
    warn("Unable to make GeneWise model");
    goto exit;
  }

  flatten_balance_scores_GeneWise(gw);

	

  /*  show_GeneWiseSegment(gw->seg[0],stderr); */
  
  if( use_syn == TRUE ) {
    if( tsm->rm == NULL ) {
      warn("Ugh - a threestatemodel without a random model. Not in this code matey");
      goto exit;
    }


    GeneWise_fold_in_synchronised_RandomModel(gw,tsm->rm,gpara->cm,gpara->ct,0.5);
    flatten_RandomCodon(gpara->rc);
  } else {
    GeneWise_fold_in_RandomModelDNA(gw,rmd);
  }

  if( alg == GWWRAP_6LITE ) {
    gwl = GwLite_from_GeneWise(gw);
    gwls = GwLiteScore_from_GwLite(gwl);
  }


  if( flat_insert == TRUE ) {
    check_flat_insert(gw,1,0,gpara->cm->ct);
  }

  if( (gws = GeneWiseScore_from_GeneWise(gw)) == NULL) {
    warn("Unable to make GeneWiseScore model");
    goto exit;
  }


  if( (gps = GeneParser21Score_from_GeneParser21(gpara->gp)) == NULL) {
    warn("Unable to make GeneParserScore model");
    goto exit;
  }


  if( (rcs = RandomCodonScore_from_RandomCodon(gpara->rc)) == NULL) {
    warn("Unable to make RandomCodonScore model");
    goto exit;
  }

  ids = folded_RandomModelDNAScore_from_2RMD(intergenic,rmd);

  gp4s = GeneParser4Score_from_GeneParser21Score(gps);

  gp4s->transition[GP4_INTRON2CDS] = Probability2Score(halfbit2Probability((-12+(5*-2))));

  if( alg == GWWRAP_333 || alg == GWWRAP_333L ) {
    cps = cDNAParserScore_from_GeneParser21Score(gps);
  }

  if( alg == GWWRAP_623S ) {
    ggms = vanilla_GeneralGeneModelScore(gpara->ct,Bits2Probability(10),1.0,Bits2Probability(10));
  }


  switch(alg) {
  case GWWRAP_2193 :

    pal = PackAln_bestmemory_GeneWise21(gws,cs,gps,rcs,ids,dpenv,dpri);

    out = convert_PackAln_to_AlnBlock_GeneWise21(pal);
    break;

  case GWWRAP_2193I :

    warn("Algorithm currently disabled! Sorry!");
    break;

    /**
    pal = PackAln_dc_build_GeneLinker21(gws,cs,gps,rcs,ids);
    if(pal == NULL)
      goto exit;
    out = convert_PackAln_to_AlnBlock_GeneLinker21(pal,NULL);
    break;
    **/

  case GWWRAP_2193L :

    pal = PackAln_bestmemory_GeneLoop21(gws,cs,gps,rcs,ids,dpenv,dpri);
    if(pal == NULL)
      goto exit;
    out = convert_PackAln_to_AlnBlock_GeneLoop21(pal);
    break;
  case GWWRAP_623L :

    pal = PackAln_bestmemory_GeneLoop6(gws,cs,gp4s,dpenv,dpri);
    if(pal == NULL)
      goto exit;
    out = convert_PackAln_to_AlnBlock_GeneLoop6(pal);
    break;
  case GWWRAP_623 :

    /*
      if( *palpoi != NULL ) {
      warn("For genewise623, Using given packaln, not calculating!");
      pal = *palpoi;
      } else {
      pal = PackAln_bestmemory_GeneWise6(gws,cs,gp4s,dpenv,dpri);
      }
    */

    gp4s->transition[GP4_INTRON2INTRON] = 0;

    /*fprintf(stderr,"Got intron2intron score of %d\n",gp4s->transition[GP4_INTRON2INTRON]);*/

    pal = PackAln_bestmemory_GeneWise6(gws,cs,gp4s,dpenv,dpri);
      
    if(pal == NULL)
      goto exit;
    out = convert_PackAln_to_AlnBlock_GeneWise6(pal);
    break;
  case GWWRAP_623S :

    pal = PackAln_bestmemory_GeneStretch6(gws,cs,gp4s,ggms,dpenv,dpri);
    if(pal == NULL)
      goto exit;
    out = convert_PackAln_to_AlnBlock_GeneStretch6(pal);
    break;

  case GWWRAP_333 :
    cps = cDNAParserScore_from_GeneParser21Score(gps);
    pal = PackAln_bestmemory_EstWise3(gws,cs,cps,dpenv,dpri);
    if(pal == NULL)
      goto exit;
    out = convert_PackAln_to_AlnBlock_EstWise3(pal);
    break;

  case GWWRAP_333L :
    cps = cDNAParserScore_from_GeneParser21Score(gps);
    pal = PackAln_bestmemory_EstLoop3(gws,cs,cps,dpenv,dpri);
    if(pal == NULL)
      goto exit;
    out = convert_PackAln_to_AlnBlock_EstLoop3(pal);
    break;
    
  case GWWRAP_421 :
    
    pal = PackAln_bestmemory_GeneWise4(gws,cs,gp4s,dpenv,dpri);
    if(pal == NULL)
      goto exit;
    out = convert_PackAln_to_AlnBlock_GeneWise4(pal);
    break;
    
  case GWWRAP_6LITE :

    pal = PackAln_bestmemory_GeneLiteModel(gwls,cs,gp4s,dpenv,dpri);
    if(pal == NULL)
      goto exit;
    out = convert_PackAln_to_AlnBlock_GeneLiteModel(pal);
    GwLite_AlnBlock_surgery(out);
    break;

  default :
    warn("A major problem. No valid algorithm type passed in");
    goto exit;
  }

  map_phase0_codons_AlnBlock_GeneWise(out,gws,cs);

  if( palpoi != NULL ) {
    *palpoi = pal;
    pal = NULL;
  }

  goto exit;


  

  exit :

  
  if(pal != NULL)
    pal = free_PackAln(pal);
  if( cps != NULL )
    cps = free_cDNAParserScore(cps);
  if( ids != NULL )
    ids = free_RandomModelDNAScore(ids);
  if(cs != NULL )
    cs = free_ComplexSequence(cs);
  if(gw != NULL )
    gw = free_GeneWise(gw);
  if(gws != NULL )
    gws = free_GeneWiseScore(gws);
  if(gps != NULL )
    free_GeneParser21Score(gps);
  if(rcs != NULL )
    rcs = free_RandomCodonScore(rcs);
  if(gp4s != NULL)
    gp4s = free_GeneParser4Score(gp4s);
  if( ggms != NULL ) 
    ggms = free_GeneralGeneModelScore(ggms);

  return out;
}

/* Function:  Hscore_from_TSM_genewise(tdb,gdb,gpara,rmd,intergenic,use_syn,alg,bits_cutoff,allN,report_level,die_on_error,flat_insert,dbsi)
 *
 * Descrip:    Runs a database search of the genewise algorithm. 
 *
 *             This makes a high score object which you can then use 
 *             to retrieve enteries as well as print out the top score (!)
 *
 *
 * Arg:                 tdb [READ ] a database of profileHMMs  [ThreeStateDB *]
 * Arg:                 gdb [READ ] a database of genomic sequence [GenomicDB *]
 * Arg:               gpara [READ ] geneparameters [GeneParameter21 *]
 * Arg:                 rmd [READ ] random model to be compared with in non syn mode [RandomModelDNA *]
 * Arg:          intergenic [READ ] random model of intergenic DNA (usually the same as rmd) [RandomModelDNA *]
 * Arg:             use_syn [UNKN ] use synchronous random model [boolean]
 * Arg:                 alg [UNKN ] algorithm type [int]
 * Arg:         bits_cutoff [UNKN ] cutoff in bits of the scores to store [double]
 * Arg:                allN [UNKN ] Undocumented argument [Probability]
 * Arg:        report_level [UNKN ] stagger rate of reporting progress on stderr  -1 means never [int]
 * Arg:        die_on_error [UNKN ] if true, exits on error (not used at the moment) [boolean]
 * Arg:         flat_insert [UNKN ] Undocumented argument [boolean]
 * Arg:                dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 *
 * Return [UNKN ]  a new Hscore object of the entire db search [Hscore *]
 *
 */
# line 515 "gwrap.dy"
Hscore * Hscore_from_TSM_genewise(ThreeStateDB * tdb,GenomicDB * gdb,GeneParameter21 * gpara,RandomModelDNA * rmd,RandomModelDNA * intergenic,boolean use_syn,int alg,double bits_cutoff,Probability allN,int report_level,boolean die_on_error,boolean flat_insert,DBSearchImpl * dbsi)
{
  Hscore * out = NULL;
  GeneWiseDB * gwdb = NULL;
  cDNADB * cdb = NULL;
  cDNAParserScore * cps= NULL;
  GeneParser21Score * gps = NULL;
  GeneParser4Score * gp4s = NULL;
  RandomCodonScore * rcs = NULL;
  RandomModelDNAScore * ids = NULL;
  cDNA * temp;
  Search_Return_Type ret;
  ComplexSequenceEval * tempcse;
  
  ret = SEARCH_ERROR;

  gwdb = new_GeneWiseDB(tdb,gpara,rmd,use_syn,allN);
  gwdb->flat_insert = flat_insert;
  if( gwdb == NULL ) {
    warn("Could not build a new GeneWiseDB from the objects provided. Exiting without completing the search");
    goto exit;
  }


  if( (gps = GeneParser21Score_from_GeneParser21(gpara->gp)) == NULL) {
    warn("Unable to make GeneParserScore model");
    goto exit;
  }
  

  if( (rcs = RandomCodonScore_from_RandomCodon(gpara->rc)) == NULL) {
    warn("Unable to make RandomCodonScore model");
    goto exit;
  }

  ids = folded_RandomModelDNAScore_from_2RMD(intergenic,rmd);

  gp4s = GeneParser4Score_from_GeneParser21Score(gps);

  if( alg == GWWRAP_333 || alg == GWWRAP_333L ) {
    /* could be a single dna sequence */
    if( gdb->is_single_seq == TRUE ) {
      temp = cDNA_from_Sequence(hard_link_Sequence(gdb->forw->seq));
      cdb = new_cDNADB_from_single_seq(temp);
      free_cDNA(temp); /* hard linked by database */
    } else {
      cdb = new_cDNADB(gdb->sdb);
    }
  }  

  /*** allocate Hscore structure ***/

  out = std_bits_Hscore(bits_cutoff,report_level);

  switch(alg) {
  case GWWRAP_2193 :

    ret = Wise2_search_GeneWise21(dbsi,out,gwdb,gdb,gps,rcs,ids);
    break;

  case GWWRAP_2193I :

    warn("Algorithm currently disabled! Sorry!");
    break;


  case GWWRAP_2193L :

    warn("Unable to do db looping searches (sort of - pointless - )!");
    break;
  case GWWRAP_623L :
    warn("Unable to do db looping searches (sort of - pointless - )!");
    break;

  case GWWRAP_623 :

    ret = Wise2_search_GeneWise6(dbsi,out,gwdb,gdb,gp4s);
    break;

  case GWWRAP_6LITE :

    ret = Wise2_search_GeneLiteModel(dbsi,out,gwdb,gdb,gp4s);
    break;

  case GWWRAP_333 :
    cps = cDNAParserScore_from_GeneParser21Score(gps);
    ret = search_EstWise3(dbsi,out,gwdb,cdb,cps);

    break;

  case GWWRAP_333L :
    warn("Unable to do db looping searches (sort of - pointless - )!");
    break;
    

  case GWWRAP_421 :

    ret = Wise2_search_GeneWise4(dbsi,out,gwdb,gdb,gp4s);
    break;


  default :
    warn("A major problem. No valid algorithm type passed in");
    goto exit;
  }

  goto exit;




  exit :

    /* for 6LITE leaking a tiny amount of memory. Oh well... */
  if( ids != NULL )
    ids = free_RandomModelDNAScore(ids);
  if( cps != NULL ) 
    free_cDNAParserScore(cps);
  if( cdb != NULL ) 
    free_cDNADB(cdb);
  if(gps != NULL )
    free_GeneParser21Score(gps);
  if(rcs != NULL )
    rcs = free_RandomCodonScore(rcs);
  if(gp4s != NULL)
    gp4s = free_GeneParser4Score(gp4s);
  if( gwdb != NULL ) {
    free_GeneWiseDB(gwdb);
  }

  if( die_on_error == TRUE  && ret == SEARCH_ERROR) {
    if( out != NULL ) {
      free_Hscore(out);
    } 
    return NULL;
  }


  return out;
}

/* Function:  cDNAParserScore_from_GeneParser21Score(gps)
 *
 * Descrip:    Makes a cdna parser from a genewise parser. Basically
 *             copies the indel penalties across.
 *
 *
 * Arg:        gps [UNKN ] Undocumented argument [GeneParser21Score *]
 *
 * Return [UNKN ]  Undocumented return value [cDNAParserScore *]
 *
 */
# line 661 "gwrap.dy"
cDNAParserScore * cDNAParserScore_from_GeneParser21Score(GeneParser21Score * gps)
{
  cDNAParserScore * out;

  out = cDNAParserScore_alloc();

  out->trans[PCD_INSERT_2_BASE] = gps->transition[GP21_INSERT_2_BASE];
  out->trans[PCD_INSERT_1_BASE] = gps->transition[GP21_INSERT_1_BASE];
  out->trans[PCD_DELETE_2_BASE] = gps->transition[GP21_DELETE_2_BASE];
  out->trans[PCD_DELETE_1_BASE] = gps->transition[GP21_DELETE_1_BASE];


  return out;
}

/* Function:  gwrap_alg_type_from_string(str)
 *
 * Descrip:    Gives you the integer interpretation from
 *             the string, which is one of
 *             2193 2193L, 623, 623L, 421, 2193LINK
 *
 *             This integer can then be passed into routines
 *             like AlnBlock_from_protein_genewise_wrap
 *
 *
 *
 * Arg:        str [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 685 "gwrap.dy"
int gwrap_alg_type_from_string(char * str)
{
  int t;

  t = get_number_from_slashed_string(str,"2193/2193L/623/623L/421/2193LINK/333/333L/6LITE/623S/623P");

  
  switch (t) {
  case 0 : return GWWRAP_2193;
  case 1 : return GWWRAP_2193L;
  case 2 : return GWWRAP_623;
  case 3 : return GWWRAP_623L;
  case 4 : return GWWRAP_421;
  case 5 : return GWWRAP_2193I;
  case 6 : return GWWRAP_333;
  case 7 : return GWWRAP_333L;
  case 8 : return GWWRAP_6LITE;
  case 9 : return GWWRAP_623S;
  case 10 : return GWWRAP_623P;
  default : warn("Cannot convert string %s into a valid genewise algorithm type\n",str);
    return -1;
  }
}


# line 707 "gwrap.c"

#ifdef _cplusplus
}
#endif
