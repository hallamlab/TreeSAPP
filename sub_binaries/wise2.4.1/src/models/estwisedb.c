
/* has to be before the others due to nasty namespace clashes */
#define WISE2_CROSS_HMMER2
#include "wise2xhmmer2.h"

#include "dyna.h"
#include "version.h" 


char * program_name = "estwisedb";
/*
 * program specific includes
 */


#include "estwrap.h"
#include "genedisplay.h"
#include "matchsum.h"


#define WISE2_DEFINED
/*#define WISE2_CROSS_HMMER1
#include "hmmio.h" 
*/

/* cross-file Wise2<->HMMER1 */


/*
 * program specific variables
 */

enum PC_SEARCH_MODE {
  PC_SEARCH_S2DB,
  PC_SEARCH_DB2S,
  PC_SEARCH_DB2DB 
};
int search_mode;

GeneWiseDB * gwdb = NULL;


char * dna_seq_file  = NULL;
cDNADB * cdb   = NULL;
SequenceDB * sdb  = NULL;
cDNA * cdna     = NULL;
boolean use_single_dna = FALSE;
boolean do_forward_only = FALSE;

char * protein_file = NULL;

Protein * pro   = NULL;
char * hmm_file       = NULL;
ThreeStateModel * tsm = NULL;
char * hmm_name       = NULL;

ThreeStateDB * tsmdb  = NULL;

GeneWiseScore * gws = NULL;

Hscore * hs       = NULL;

CodonMapper * cm;
cDNAParser * cps;

/** different protein possibilities **/
boolean use_single_pro = FALSE;
boolean use_db_pro     = FALSE;
boolean use_tsm        = FALSE;
boolean use_pfam1      = FALSE;
boolean use_pfam2      = FALSE;


char * qstart_str    = NULL;
int qstart           = -1;

char * qend_str      = NULL;
int qend             = -1;


char * matrix_file      = "BLOSUM62.bla";
CompMat * mat           = NULL;

char * gap_str          = "12";
int gap                 = 12;

char * ext_str          = "2";
int ext                 = 2;

char * codon_file       = NULL;
CodonTable * ct         = NULL;

char * output_file      = "-";
FILE * ofp              = NULL;

char * report_str       = NULL;
int report_stagger      = -1;

RandomModelDNA * rmd    = NULL;

char * subs_string      = "0.01";
double subs_error       = 0.01;

char * indel_string      = "0.01";
double indel_error       = 0.01;


char * allN_string      = "1.0";
Probability allN        = 1.0;

boolean flat_insert     = FALSE;

char * startend_string   = "default";
int startend             = TSM_default;

char * null_string       = "syn";
boolean use_syn          = TRUE;

DBSearchImpl * dbsi      = NULL;
int alg                  = ESTSLIM_3;
char * alg_str           = NULL;

DPRunImpl * dpri         = NULL;

int aln_alg              = ESTWISE_3;
char * aln_alg_str       = NULL;

int aln_number           = 50;
char * aln_number_str    = "50";

double  search_cutoff    = 20.00;
char * search_cutoff_str   = "20.00"; 

double evalue_search_cutoff = -1.0;
char * evalue_search_str = NULL;

char * kbyte_str         = NULL;
int kbyte                = 10000; /* will be reset in build_defaults */

boolean show_histogram   = TRUE;

boolean show_PackAln     = FALSE;
boolean show_AlnBlock    = FALSE;
boolean show_pretty      = FALSE;
boolean show_pep         = FALSE;
boolean show_match_sum   = FALSE;
boolean show_para        = FALSE;

boolean do_complete_analysis = FALSE;
boolean make_anchored_aln = FALSE;

char * main_block_str      = "50";
int main_block           = 50;

char * divide_str        = "//";

Probability rnd_loop      = 0.99;
Probability cds_loop      = 0.97;
Probability rnd_to_model  = (1 - 0.99) / 3;
Probability link_loop     = 0.98;
Probability link_to_model = (1- 0.98) / 3;

AlnBlock * alb;
PackAln  * pal;

MatchSummarySet * mss;

RandomModel * rm;

boolean show_output(void)
{
  int i,k;
  ThreeStateModel * temptsm;
  AlnBlock * alb;
  PackAln * pal;
  MatchSummarySet * mss;
  Protein * ps;
  cDNA * cdna;
  double bits;
  boolean fitted_res = FALSE;
  AlnBlockList * alist;
  AlnBlock * anchored;
  SequenceSet * set;
  AlnColumn * alt;
  Protein * trans;

  /* sort by bit score first */

  sort_Hscore_by_score(hs);

  if( search_mode == PC_SEARCH_S2DB ) {
    if( hs->his == NULL || hs->his->total < 1000 ) {
	info("Cannot fit histogram to a db smaller than 1,000");
	fprintf(ofp,"[Warning: Can't fit histogram to a db smaller than 1,000]\n\n");
	show_histogram = FALSE;
    } else {
      fitted_res = TRUE;
      fit_Hscore_to_EVD(hs,20);
    }
  }

  /* deal with initialising anchored alignment.
   * Could be done for either single HMMs or single proteins,
   * but we will only do it for HMMs at the moment
   */

  if( make_anchored_aln == TRUE ) {
    if( tsm == NULL ) {
      warn("Attempting to make an achored alignment without a HMM. impossible!");
      make_anchored_aln = FALSE;
    } else {
      anchored = single_unit_AlnBlock(tsm->len,"MATCH_STATE");
      set = SequenceSet_alloc_std();
   }
  }

  /* dofus catcher */
  if( aln_alg != alg ) {
    fprintf(ofp,"\n#\n#WARNING!\n#\n# Your alignment algorithm is different from your search algorithm.\n# This is probably quite sensible but will lead to differing scores.\n# Use the search score as an indicator of the significance of the match\n# Read the docs for more information\n#\n");
  }

  fprintf(ofp,"\n\n#High Score list\n");
  fprintf(ofp,"#Protein ID                 DNA Str  ID                        Bits Evalue\n");  
  fprintf(ofp,"--------------------------------------------------------------------------\n");

  for(i=0;i<hs->len;i++) {
    bits = Score2Bits(hs->ds[i]->score);
    if( bits < search_cutoff ) {
      break;
    }

    if( fitted_res == TRUE && evalue_search_str != NULL ) {
      if( hs->ds[i]->evalue > evalue_search_cutoff ) 
	break;
    }

    if( fitted_res == TRUE) 
      fprintf(ofp,"Protein %-20sDNA [%c] %-24s %.2f %.2g\n",hs->ds[i]->query->name,hs->ds[i]->target->is_reversed == TRUE ? '-' : '+',hs->ds[i]->target->name,bits,hs->ds[i]->evalue);
    else
      fprintf(ofp,"Protein %-20sDNA [%c] %-24s %.2f\n",hs->ds[i]->query->name,hs->ds[i]->target->is_reversed == TRUE ? '-' : '+',hs->ds[i]->target->name,bits);

  }

  if( search_mode == PC_SEARCH_S2DB && show_histogram == TRUE ) {
    fprintf(ofp,"\n\n#Histogram\n");
    fprintf(ofp,"-----------------------------------------------------------------------\n");
    PrintASCIIHistogram(hs->his,ofp);
  }

  fprintf(ofp,"\n\n#Alignments\n");
  fprintf(ofp,"-----------------------------------------------------------------------\n");

  for(i=0;i<hs->len;i++) {
    bits = Score2Bits(hs->ds[i]->score);
    if( bits < search_cutoff ) {
      break;
    }
    if( i >= aln_number ) {
      break;
    }

    if( fitted_res == TRUE && evalue_search_str != NULL ) {
      if( hs->ds[i]->evalue > evalue_search_cutoff ) 
	break;
    }

    
    fprintf(ofp,"\n\n>Results for %s vs %s (%s) [%d]\n",hs->ds[i]->query->name,hs->ds[i]->target->name,hs->ds[i]->target->is_reversed == TRUE ? "reverse" : "forward",i+1 );

    cdna = get_cDNA_from_cDNADB(cdb,hs->ds[i]->target);
    temptsm = indexed_ThreeStateModel_ThreeStateDB(tsmdb,hs->ds[i]->query);


    alb = AlnBlock_from_TSM_estwise_wrap(temptsm,cdna,cps,cm,ct,rmd,aln_alg,use_syn,allN,flat_insert,dpri,&pal);

    if( alb == NULL ) {
      warn("Got a NULL alignment. Exiting now due to presumed problems");
      fprintf(ofp,"\n\n*Got a NULL alignment. Exiting now due to presumed problems*\n\n");
      return FALSE;
    }


 
    if( use_single_pro == FALSE) 
      mss = MatchSummarySet_from_AlnBlock_genewise(alb,temptsm->name,1,cdna->baseseq);
    else
      mss = MatchSummarySet_from_AlnBlock_genewise(alb,pro->baseseq->name,pro->baseseq->offset,cdna->baseseq);

    
    if( show_pretty == TRUE ) {

      fprintf(ofp,"\n%s output\nScore %4.2f bits over entire alignment.\nThis will be different from per-alignment scores. See manual for details\nFor computer parsable output, try %s -help or read the manual\n",program_name,Score2Bits(pal->score),program_name);
      
      if( use_syn == FALSE ) {
	fprintf(ofp,"Scores as bits over a flat simple random model\n\n");
      } else {
	fprintf(ofp,"Scores as bits over a synchronous coding model\n\n");
      }
      
      ps = pseudo_Protein_from_ThreeStateModel(temptsm);
      protcdna_ascii_display(alb,ps->baseseq->seq,ps->baseseq->name,ps->baseseq->offset,cdna,ct,15,main_block,TRUE,ofp);

      
      free_Protein(ps);

      fprintf(ofp,"%s\n",divide_str);
      
    }

    if( show_match_sum == TRUE ) {
      show_MatchSummary_genewise_header(ofp);
      show_MatchSummarySet_genewise(mss,ofp);
      fprintf(ofp,"%s\n",divide_str);
    }
    

    if( show_pep == TRUE ) {
      alt = alb->start;
      for(;alt != NULL;) {
	trans = Protein_from_GeneWise_AlnColumn(cdna->baseseq,alt,1,&alt,ct,is_random_AlnColumn_genewise);
	if ( trans == NULL ) 
	  break;
	write_fasta_Sequence(trans->baseseq,ofp);
	free_Protein(trans);
      }
      fprintf(ofp,"%s\n",divide_str);
    }

    if( show_AlnBlock == TRUE ) {
      mapped_ascii_AlnBlock(alb,Score2Bits,0,ofp);
      fprintf(ofp,"%s\n",divide_str);
    }
    
    if( show_PackAln == TRUE ) {
      show_simple_PackAln(pal,ofp);
      fprintf(ofp,"%s\n",divide_str);
    }

    /*
     * This goes at the end because it destroys the alb structure
     */

    if( make_anchored_aln == TRUE ) {
      /* attach sequence to als in alb, so we have it for later use */
      alb->seq[1]->data = (void *) cdna->baseseq;
      /* add to SequenceSet so we can destroy the memory */
      add_SequenceSet(set,hard_link_Sequence(cdna->baseseq));

      alist = split_AlnBlock(alb,is_random_AlnColumn_genewise);

      for(k=0;k<alist->len;k++) {
	/* actually produce the anchored alignment */
	/*mapped_ascii_AlnBlock(alist->alb[k],Score2Bits,stderr);*/
	add_to_anchored_AlnBlock(anchored,alist->alb[k]);

	/*	dump_ascii_AlnBlock(anchored,stderr);*/
      }
    }

    alb = free_AlnBlock(alb);
    pal = free_PackAln(pal);
    mss = free_MatchSummarySet(mss);
    cdna = free_cDNA(cdna);
    temptsm = free_ThreeStateModel(temptsm);

  }

  if( do_complete_analysis == TRUE ) {
    fprintf(ofp,"\n\n#Complete Analysis\n");
    fprintf(ofp,"-------------------------------------------------------------\n\n");
    
    /* ok - end of loop over relevant hits. If we have an
     * anchored alignment, print it out!
     */
    if( make_anchored_aln == TRUE ) {
      /*dump_ascii_AlnBlock(anchored,stderr);*/
      write_mul_estwise_AlnBlock(anchored,ct,ofp);
      fprintf(ofp,"%s\n",divide_str);
    }
  }


  return TRUE;
}
    
boolean search_db(void)
{

  info("Starting search...");

  hs = Hscore_from_TSM_estwise(tsmdb,cdb,cps,cm,rmd,use_syn,alg,search_cutoff,allN,flat_insert,report_stagger,FALSE,dbsi);

  if( hs == NULL ) {
    return FALSE;
  }

  return TRUE;
}

boolean build_db_objects(void)
{

  if( use_single_dna == TRUE ) {
    cdb = new_cDNADB_from_single_seq(cdna);
  } else {
    cdb = new_cDNADB(sdb);
  }
  if( do_forward_only == TRUE ) {
    cdb->forward_only = TRUE;
  }

  return TRUE;
}


boolean build_objects(void)
{
  boolean ret = TRUE;
  Protein * pro_temp;
  SequenceDB * psdb;



  startend = threestatemodel_mode_from_string(startend_string);
  if( startend == TSM_unknown ) {
    warn("String %s was unable to converted into a start/end policy\n",startend_string);
    ret = FALSE;
  }

  if( use_single_dna == TRUE ) {
    cdna = read_fasta_file_cDNA(dna_seq_file);
    if( cdna == NULL ) {
      warn("Could not open single dna sequence in %s",dna_seq_file);
      ret = FALSE;
    }
  } else {
    sdb = single_fasta_SequenceDB(dna_seq_file);
    
 
    if( sdb == NULL ) {
      warn("Could not build a sequence database on %s",dna_seq_file);
      ret = FALSE;
    }
  }

  rm = default_RandomModel();


  if( (mat = read_Blast_file_CompMat(matrix_file)) == NULL) {
    if( use_tsm == TRUE ) {
      info("I could not read the Comparison matrix file in %s; however, you are using a HMM so it is not needed. Please set the WISECONFIGDIR or WISEPERSONALDIR variable correctly to prevent this message.",matrix_file);
    } else {
      warn("Could not read Comparison matrix file in %s",matrix_file);
      ret = FALSE;
    }
  }
      
  if( is_integer_string(gap_str,&gap) == FALSE ) {
    warn("Could not get gap string number %s",gap_str);
    ret = FALSE;
  }

  if( is_integer_string(ext_str,&ext) == FALSE ) {
    warn("Could not get ext string number %s",ext_str);
    ret = FALSE;
  }

  if( qstart_str != NULL ) {
    if( is_integer_string(qstart_str,&qstart) == FALSE || qstart < 0) {
      warn("Could not make %s out as query start",qstart);
      ret = FALSE;
    }
  }

  if( qend_str != NULL ) {
    if( is_integer_string(qend_str,&qend) == FALSE || qend < 0) {
      warn("Could not make %s out as query end",qend);
      ret = FALSE;
    }
  }


  if( aln_number_str != NULL ) {
    if( is_integer_string(aln_number_str,&aln_number) == FALSE || aln_number < 0) {
      warn("Weird aln number string %s...\n",aln_number_str);
      ret = FALSE;
    }
  }

  if( report_str != NULL ) {
    if( is_integer_string(report_str,&report_stagger) == FALSE ) {
      warn("Weird report stagger asked for %s",report_str);
      ret = FALSE;
    }
  }


  if( use_pfam1 == TRUE ) {
    tsmdb = new_PfamHmmer1DB_ThreeStateDB(protein_file);
    if( set_search_type_ThreeStateDB(tsmdb,startend_string) == FALSE) {
      warn("Unable to set global/local switch on threestatedb");
      ret = FALSE;
    }

  } else if ( use_pfam2 == TRUE ) {
    tsmdb = HMMer2_ThreeStateDB(protein_file);
    if( set_search_type_ThreeStateDB(tsmdb,startend_string) == FALSE) {
      warn("Unable to set global/local switch on threestatedb");
      ret = FALSE;
    }

  } else if ( use_tsm == TRUE) {
    /** using a HMM **/

    tsm = HMMer2_read_ThreeStateModel(protein_file);

    if( tsm == NULL ) {
      warn("Could not read hmm from %s\n",protein_file);
      ret = FALSE;
    }  else {

      display_char_in_ThreeStateModel(tsm);
      if( hmm_name != NULL ) {
	if( tsm->name != NULL ) 
	  ckfree(tsm->name);
	tsm->name = stringalloc(hmm_name);
      } else {
	if( tsm->name == NULL ) {
	  tsm->name = stringalloc(protein_file);
	}
      }

      
      
      /** have to set start/end **/

      set_startend_policy_ThreeStateModel(tsm,startend,15,0.2);
      tsmdb = new_single_ThreeStateDB(tsm,rm);
      if( tsmdb == NULL ) {
	warn("Could not build a threestatemodel database from a single tsm. Weird!");
	ret = FALSE;
      }
    } /* end of else tsm != NULL */
  } /* end of else is tsm */
  else if( use_single_pro ) {


    if( startend != TSM_default && startend != TSM_global && startend != TSM_local ) {
      warn("Proteins can only have local/global startend policies set, not %s",startend_string);
      ret = FALSE;
    }

    if( (pro = read_fasta_file_Protein(protein_file)) == NULL ) {
      ret = FALSE;
      warn("Could not read Protein sequence in %s",protein_file);
    } else {
      if( qstart != -1 || qend != -1 ) {
	if( qstart == -1 )
	  qstart = 0;
	if( qend == -1 ) 
	  qend = pro->baseseq->len;

	pro_temp = truncate_Protein(pro,qstart-1,qend);
	if( pro_temp == NULL ){
	  ret = FALSE;
	} else {
	  free_Protein(pro);
	  pro = pro_temp;
	}
      }


      if( startend == TSM_global) 
	tsm = global_ThreeStateModel_from_half_bit_Sequence(pro,mat,rm,-gap,-ext);
      else
	tsm = ThreeStateModel_from_half_bit_Sequence(pro,mat,rm,-gap,-ext);

      if( tsm == NULL ) {
	warn("Could not build ThreeStateModel from a single protein sequence...");
	ret = FALSE; 
      } else {
	tsmdb = new_single_ThreeStateDB(tsm,rm);
	if( tsmdb == NULL ) {
	  warn("Could not build a threestatemodel database from a single tsm. Weird!");
	  ret = FALSE;
	}
      } /* end of could build a TSM */
    } /* else is a real protein */  

  } /* end of else is single protein */
  else if (use_db_pro == TRUE ) {
    psdb = single_fasta_SequenceDB(protein_file);
    tsmdb = new_proteindb_ThreeStateDB(psdb,mat,-gap,-ext);
    free_SequenceDB(psdb);
  }
  else {
    warn("No protein input file! Yikes!");
  }

  /***
  if( use_tsm == FALSE ) {
  } else {
  ****/


  if( main_block_str != NULL ) {
    if( is_integer_string(main_block_str,&main_block) == FALSE ) {
      warn("Could not get maximum main_block number %s",main_block_str);
      ret = FALSE;
    }
  }


  if( evalue_search_str != NULL && is_double_string(evalue_search_str,&evalue_search_cutoff) == FALSE ) {
    warn("Could not convert %s to a double",evalue_search_str);
    ret = FALSE;
  }
  
  if( is_double_string(search_cutoff_str,&search_cutoff) == FALSE ) {
    warn("Could not convert %s to a double",search_cutoff_str);
    ret = FALSE;
  }


  if( is_double_string(subs_string,&subs_error) == FALSE ) {
    warn("Could not convert %s to a double",subs_error);
    ret = FALSE;
  }

  if( is_double_string(indel_string,&indel_error) == FALSE ) {
    warn("Could not convert %s to a double",indel_error);
    ret = FALSE;
  }


  if( is_double_string(allN_string,&allN) == FALSE ) {
    warn("Could not convert %s to a double",allN_string);
    ret = FALSE;
  }
  


  if( strcmp(null_string,"syn") == 0 ) {
    use_syn = TRUE;
  } else if ( strcmp(null_string,"flat") == 0 ) {
    use_syn = FALSE;
  } else {
    warn("Cannot interpret [%s] as a null model string\n",null_string);
    ret = FALSE;
  }

   
  if( alg_str != NULL ) {
    alg = alg_estwrap_from_string(alg_str);
  } else {
    alg_str = "312";
    alg = alg_estwrap_from_string(alg_str);
  }

  if( aln_alg_str != NULL ) {
    aln_alg = alg_estwrap_from_string(aln_alg_str);
  } else {
    /* if it is a protein, don't loop */
    if( use_single_pro == TRUE || use_db_pro == TRUE ) 
      aln_alg_str = "333";
    else 
      aln_alg_str = "333L";
    aln_alg = alg_estwrap_from_string(aln_alg_str);
  }


  if( (rm = default_RandomModel()) == NULL) {
    warn("Could not make default random model\n");
    ret = FALSE;
  }

  if( (ct = read_CodonTable_file(codon_file)) == NULL) {
    ret = FALSE;
    warn("Could not read codon table file in %s",codon_file);
  }

  if( (ofp = openfile(output_file,"W")) ==  NULL) {
    warn("Could not open %s as an output file",output_file);
    ret = FALSE;
  }

  rmd = RandomModelDNA_std();


  cps = flat_cDNAParser(indel_error);
  cm = flat_CodonMapper(ct);
  sprinkle_errors_over_CodonMapper(cm,subs_error);

  return ret;

}

void free_objects(void)
{
  if( gwdb != NULL ) 
    gwdb = free_GeneWiseDB(gwdb);
  if( cdb != NULL ) 
    cdb = free_cDNADB(cdb);
  if( sdb != NULL )
    sdb = free_SequenceDB(sdb);
  if( cdna != NULL ) 
    cdna = free_cDNA(cdna);
  if( pro != NULL )
    pro = free_Protein(pro);
  if( tsm != NULL )
    tsm = free_ThreeStateModel(tsm);
  if( tsmdb != NULL )
    tsmdb = free_ThreeStateDB(tsmdb);
  if( gws != NULL )
    gws = free_GeneWiseScore(gws);
  if( hs != NULL )
    hs = free_Hscore(hs);
  if( cm != NULL )
    cm = free_CodonMapper(cm);
  if( cps != NULL )
    cps = free_cDNAParser(cps);

}

void show_short_help(void)
{
  fprintf(stdout,"%s (%s)\n",program_name,VERSION_NUMBER);
  fprintf(stdout,"This program is freely distributed under a GPL. See -version for more info\n");
  fprintf(stdout,"Copyright (c) GRL limited: portions of the code are from separate copyrights\n\n");
  fprintf(stdout,"spcwise <protein-input> <dna-input> in fasta format\n");
  fprintf(stdout," Options. In any order, '-' as filename (for any input) means stdin\nDon't use stdin for databases, as on-the-fly indexing is used\n");
  fprintf(stdout," Protein type  [-protein,-prodb,-hmmer,-pfam] [default - protein]\n");
  fprintf(stdout," Dna type      [-dnas,-dnadb] [default - dnadb]\n");
  fprintf(stdout," Dna           [-tfor]\n");
  fprintf(stdout," Protein  [-s,-t,-g,-e,-m]\n HMM      [-hmmer,-hname]\n");
  fprintf(stdout," Model    [-codon,-subs,-indel,-null]\n Alg      [-kbyte,-alg,-aalg,-aln,-noh]\n");
  fprintf(stdout," Output   [-pretty,-alb,-pal,-block,-divide]\n");
  fprintf(stdout," Standard [-help,-version,-silent,-quiet,-errorlog]\n");
  fprintf(stdout,"\nFor more help go %s -help.\n",program_name);
  fprintf(stdout,"\nSee WWW help at http://www.sanger.ac.uk/Software/Wise2/\n");
  exit(63);   
}

void show_help(FILE * ofp)
{
  fprintf(ofp,"%s (%s)\n",program_name,VERSION_NUMBER);
  fprintf(ofp,"%s <protein-input> <dna-input>\n",program_name);
  /* program specific help */
  fprintf(ofp,"Protein input type\n");
  fprintf(ofp,"  -protein  [default] single protein\n");
  fprintf(ofp,"  -prodb    protein fasta format db\n");
  fprintf(ofp,"  -pfam     pfam hmm library \n");
  fprintf(ofp,"  -pfam2    pfam style model directory (2.1) \n");
  fprintf(ofp,"  -hmmer    single hmmer 1.x HMM\n");
  fprintf(ofp,"DNA input type\n");
  fprintf(ofp,"  -dnadb    [default] dna fasta database\n");
  fprintf(ofp,"  -dnas     a single dna fasta sequence\n");
  fprintf(ofp,"DNA sequence options\n");
  fprintf(ofp,"  -tfor     search forward strands only\n");
  fprintf(ofp,"Protein comparison options\n");
  fprintf(ofp,"  -gap      [%3d]  gap penalty\n",gap);
  fprintf(ofp,"  -ext      [%3d]  extension penalty\n",ext);
  fprintf(ofp,"  -matrix   [%s]  Comparison matrix\n",matrix_file);
  fprintf(ofp,"HMM options\n");
  fprintf(ofp,"  -hname           For single hmms, use this as the name, not filename\n");
  fprintf(ofp,"Model options\n");
  fprintf(ofp,"  -init   [%s] [default/global/local/wing] start-end policy\n",startend_string);
  fprintf(ofp,"  -codon  [%s]  Codon file\n",codon_file);
  fprintf(ofp,"  -subs   [%2.2g] Substitution error rate\n",subs_error);
  fprintf(ofp,"  -indel  [%2.2g] Insertion/deletion error rate\n",indel_error);
  fprintf(ofp,"  -null   [syn/flat]   Random Model as synchronous or flat [default syn]\n");
  fprintf(ofp,"  -alln   [%s]   Probability of matching a NNN codon\n",allN_string);
  fprintf(ofp,"  -flati         Flat insert probabilities\n");
  fprintf(ofp,"Algorithm options\n");
  fprintf(ofp,"  -alg    [333/312]         Algorithm used for searching [default %s]\n",string_from_alg_estwrap(alg));
  fprintf(ofp,"  -aalg   [312/333/333L]    Algorithm used for alignment [default %s]\n",string_from_alg_estwrap(aln_alg));
  fprintf(ofp,"  -cut    [%.2f]   Bits cutoff for reporting in search algorithm\n",search_cutoff);
  fprintf(ofp,"  -ecut   [n/a]    Evalue cutoff for single protein vs DNA searches.\n");
  fprintf(ofp,"  -aln    [%d]   Max number of alignments (even if above cut)\n",aln_number);
  fprintf(ofp,"  -nohis           Don't show histogram on single protein/hmm vs DNA search\n");
  fprintf(ofp,"  -report [0]      Issue a report every x comparisons (default 0 comparisons)\n");
  fprintf(ofp,"Output options for each alignment [default -pretty -para]\n");
  fprintf(ofp,"  -pretty          show pretty ascii output\n");
  fprintf(ofp,"  -para            show parameters\n");
  fprintf(ofp,"  -pep             show protein translation, splicing frameshifts\n");
  fprintf(ofp,"  -mul             protein mul format alignments [only for one HMM vs DNA db]\n");
  fprintf(ofp,"  -sum             show summary output\n");
  fprintf(ofp,"  -alb             show logical AlnBlock alignment\n");
  fprintf(ofp,"  -pal             show raw matrix alignment\n");
  fprintf(ofp,"  -block  [%s]     Length of main block in pretty output\n",main_block_str);
  fprintf(ofp,"  -divide [%s]     divide string for multiple outputs\n",divide_str);

  show_help_DBSearchImpl(ofp);
  show_help_DPRunImpl(ofp);
  show_standard_options(ofp);

  fprintf(ofp,"\nSee WWW help at http://www.sanger.ac.uk/Software/Wise2/\n");
  exit(63);   
}


boolean show_header(FILE * ofp)
{
  fprintf(ofp,"-------------------------------------------------------------\n");
  fprintf(ofp,"Wise2 - database searching mode\n");
  fprintf(ofp,"Program: %s version: %s released: %s\n",program_name,VERSION_NUMBER,RELEASE_DAY);
  fprintf(ofp,"This program is freely distributed under a Gnu Public License.\n");
  fprintf(ofp,"   See -version for more info on copyright\n");
  fprintf(ofp,"Bugs and credits to Ewan Birney <birney@sanger.ac.uk>\n");
  fprintf(ofp,"-------------------------------------------------------------\n\n");
  fprintf(ofp,"Algorithm type:        EstWise\n");
  fprintf(ofp,"Search algorithm:      %s\n",alg_str);
  fprintf(ofp,"Implementation:        %s\n",impl_string_DBSearchImpl(dbsi));
  fprintf(ofp,"Search mode:           %s\n",search_mode == PC_SEARCH_S2DB ? "Single protein vs cdna db" : search_mode == PC_SEARCH_DB2S ? "Single cdna vs protein db" : "Protein db vs cdna db");
  fprintf(ofp,"Protein info from:     %s\n",protein_file);
  fprintf(ofp,"Dna info from:         %s\n",dna_seq_file);
  if( use_single_pro == TRUE || use_db_pro == TRUE ) {
    fprintf(ofp,"Comp Matrix:           %s\n",matrix_file);
    fprintf(ofp,"Gap open:              %d\n",gap);
    fprintf(ofp,"Gap extension:         %d\n",ext);
  }
  fprintf(ofp,"Start/End              %s\n",startend_string);
  fprintf(ofp,"Codon Table:           %s\n",codon_file);
  fprintf(ofp,"Subs error:            %2.2g\n",subs_error);
  fprintf(ofp,"Indel error:           %2.2g\n",indel_error);
  fprintf(ofp,"Null model:            %s\n",use_syn == FALSE ? "flat" : "synchronous");
  fprintf(ofp,"Protein Insertion:     %s\n",flat_insert == TRUE ? "flat" : "modelled");
  fprintf(ofp,"Alignment Alg          %s\n",aln_alg_str);
  
  return TRUE;
}

void show_version(FILE * ofp)
{
  fprintf(ofp,"%s\nVersion: %s\nReleased: %s\nCompiled: %s\n",program_name,VERSION_NUMBER,RELEASE_DAY,COMPILE_DATE);
  fprintf(ofp,"\nThis program is freely distributed under a Gnu Public License\n");
  fprintf(ofp,"The source code is copyright (c) GRL 1998 and others\n");
  fprintf(ofp,"There is no warranty, implied or otherwise on the performance of this program\n");
  fprintf(ofp,"For more information read the GNULICENSE file in the distribution\n\n");
  fprintf(ofp,"Credits: Ewan Birney <birney@sanger.ac.uk> wrote the core code.\n");
  fprintf(ofp,"         Portions of this code was from HMMer1, written by Sean Eddy\n");
  fprintf(ofp,"         Portions of this code was from HMMer2, written by Sean Eddy\n");
  exit(63);   
}


void build_defaults(void)
{
  codon_file = "codon.table";
  matrix_file = "BLOSUM62.bla";
  

}

int main(int argc,char ** argv) 
{
  int i;
  char * temp;

  build_defaults();

  bootstrap_HMMer2();
  
  strip_out_standard_options(&argc,argv,show_help,show_version);

  if( (temp = strip_out_assigned_argument(&argc,argv,"gap")) != NULL )
    gap_str = temp;

  if( (temp = strip_out_assigned_argument(&argc,argv,"g")) != NULL )
    gap_str = temp;

  if( (temp = strip_out_assigned_argument(&argc,argv,"ext")) != NULL )
    ext_str = temp;

  if( (temp = strip_out_assigned_argument(&argc,argv,"e")) != NULL )
    ext_str = temp;

  if( (temp = strip_out_assigned_argument(&argc,argv,"matrix")) != NULL )
    matrix_file = temp;

  if( (temp = strip_out_assigned_argument(&argc,argv,"m")) != NULL )
    matrix_file = temp;

  if( (temp = strip_out_assigned_argument(&argc,argv,"s")) != NULL )
    qstart_str = temp;

  if( (temp = strip_out_assigned_argument(&argc,argv,"t")) != NULL )
    qend_str = temp;

  if( (temp = strip_out_assigned_argument(&argc,argv,"aln")) != NULL )
    aln_number_str = temp;

  if( (temp = strip_out_assigned_argument(&argc,argv,"codon")) != NULL )
    codon_file = temp;


  if( (temp = strip_out_assigned_argument(&argc,argv,"alg")) != NULL )
    alg_str = temp;

  if( (temp = strip_out_assigned_argument(&argc,argv,"aalg")) != NULL )
    aln_alg_str = temp;

  if( (temp = strip_out_assigned_argument(&argc,argv,"cut")) != NULL )
    search_cutoff_str = temp;


  if( (temp = strip_out_assigned_argument(&argc,argv,"ecut")) != NULL )
    evalue_search_str = temp;

  if( (temp = strip_out_assigned_argument(&argc,argv,"subs")) != NULL )
    subs_string = temp;

  if( (temp = strip_out_assigned_argument(&argc,argv,"indel")) != NULL )
    indel_string = temp;

  if( (temp = strip_out_assigned_argument(&argc,argv,"init")) != NULL )
    startend_string = temp;

  if( (temp = strip_out_assigned_argument(&argc,argv,"alln")) != NULL )
    allN_string = temp;

  if( (temp = strip_out_assigned_argument(&argc,argv,"null")) != NULL )
    null_string = temp;

  if( (strip_out_boolean_argument(&argc,argv,"dnas")) == TRUE )
    use_single_dna = TRUE;

  if( (strip_out_boolean_argument(&argc,argv,"dnadb")) == TRUE )
    use_single_dna = FALSE;

  if( (strip_out_boolean_argument(&argc,argv,"tfor")) == TRUE )
    do_forward_only = TRUE;

  if( (strip_out_boolean_argument(&argc,argv,"flati")) == TRUE )
    flat_insert = TRUE;

  if( (strip_out_boolean_argument(&argc,argv,"hmmer")) == TRUE )
    use_tsm = TRUE;

  if( (strip_out_boolean_argument(&argc,argv,"pfam2")) == TRUE )
    use_pfam1 = TRUE;

  if( (strip_out_boolean_argument(&argc,argv,"pfam")) == TRUE )
    use_pfam2 = TRUE;

  if( (strip_out_boolean_argument(&argc,argv,"protein")) == TRUE )
    use_single_pro = TRUE;

  if( (strip_out_boolean_argument(&argc,argv,"prodb")) == TRUE )
    use_db_pro = TRUE;


  if( (temp = strip_out_assigned_argument(&argc,argv,"hname")) != NULL )
    hmm_name = temp;

  if( (strip_out_boolean_argument(&argc,argv,"nohis")) != FALSE )
    show_histogram = FALSE;

  if( (strip_out_boolean_argument(&argc,argv,"pretty")) != FALSE )
    show_pretty = TRUE;

  if( (strip_out_boolean_argument(&argc,argv,"pep")) != FALSE )
    show_pep = TRUE;

  if( (strip_out_boolean_argument(&argc,argv,"mul")) != FALSE )
    make_anchored_aln = TRUE;

  if( (strip_out_boolean_argument(&argc,argv,"para")) != FALSE )
    show_para = TRUE;

  if( (strip_out_boolean_argument(&argc,argv,"sum")) != FALSE )
    show_match_sum = TRUE;

  if( (strip_out_boolean_argument(&argc,argv,"alb")) != FALSE )
    show_AlnBlock = TRUE;

  if( (strip_out_boolean_argument(&argc,argv,"pal")) != FALSE )
    show_PackAln = TRUE;

  if( (temp = strip_out_assigned_argument(&argc,argv,"divide")) != NULL )
    divide_str = temp;

  if( (temp = strip_out_assigned_argument(&argc,argv,"block")) != NULL )
    main_block_str = temp;


  if( (temp = strip_out_assigned_argument(&argc,argv,"report")) != NULL )
    report_str = temp;

  dbsi = new_DBSearchImpl_from_argv(&argc,argv);
  
  dpri = new_DPRunImpl_from_argv(&argc,argv);


  strip_out_remaining_options_with_warning(&argc,argv);
  

  if( argc !=  3 ) {
    warn("Wrong number of arguments (expect 2)!\n");
    if( argc > 1 ){
      warn("Arg line looked like (after option processing)");
      for(i=1;i<argc;i++) {
	fprintf(stderr,"   %s\n",argv[i]);
      }
    }

    show_short_help();
  }

  if( show_pretty == FALSE && show_AlnBlock == FALSE && show_PackAln == FALSE && show_pep == FALSE ) {
    show_pretty = TRUE;
    show_para = TRUE;
  }

  if( use_db_pro == FALSE && use_single_pro == FALSE && use_tsm == FALSE && use_pfam1 == FALSE && use_pfam2 == FALSE ) {
    use_single_pro = TRUE;
  }

  if( use_single_pro == TRUE || use_tsm == TRUE ) {
    if( use_single_dna == TRUE ) 
      fatal("one on one search. Shouldn't you use pcwise?");
    search_mode = PC_SEARCH_S2DB;
  } else {
    if( use_single_dna == TRUE ) 
      search_mode = PC_SEARCH_DB2S;
    else 
      search_mode = PC_SEARCH_DB2DB;
  }

  if( evalue_search_str != NULL && search_mode != PC_SEARCH_S2DB ) {
    fatal("Trying to set a evalue cutoff on a non evalue based search. you can only use evalues in a protein HMM vs DNA database search (sorry!)");
  }

  if( make_anchored_aln == TRUE && search_mode != PC_SEARCH_S2DB ) {
    fatal("Trying to make an anchored alignment and not in single search mode");
  }

  if( make_anchored_aln == TRUE) {
    do_complete_analysis = TRUE;
  }

  /* pick up remaining args and do it */

    
  dna_seq_file = argv[2];
  protein_file = argv[1];

  if( build_objects() == FALSE) 
    fatal("Could not build objects!");

  if( build_db_objects() == FALSE) 
    fatal("Could not build database-ready objects!");


  show_header(stdout);

  if( search_db() == FALSE) 
    warn("Could not search database");


  show_output();


  free_objects();


  return 0;
}






