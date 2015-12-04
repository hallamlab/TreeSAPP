
/* has to be before the others due to nasty namespace clashes */
#define WISE2_CROSS_HMMER2
#include "wise2xhmmer2.h"

#include "dyna.h"
#include "version.h" 

#include "gwrap.h"

char * program_name = "genewisedb";

/*
 * program specific includes
 */


#include "gwrap.h"
#include "genedisplay.h"
#include "matchsum.h"



/*
 * program specific variables
 */

enum PG_SEARCH_MODE {
  PG_SEARCH_S2DB,
  PG_SEARCH_DB2S,
  PG_SEARCH_DB2DB 
};
int search_mode;

GeneWiseDB * gwdb = NULL;


char * dna_seq_file  = NULL;
GenomicDB * gdb   = NULL;
SequenceDB * sdb  = NULL;
Genomic * gen     = NULL;
boolean use_single_dna = FALSE;

int dna_seqdb_start = -1;
int dna_seqdb_end   = -1;

char * protein_file = NULL;

Protein * pro   = NULL;
char * hmm_file       = NULL;
ThreeStateModel * tsm = NULL;
char * hmm_name       = NULL;

ThreeStateDB * tsmdb  = NULL;

int pro_seqdb_start = -1;
int pro_seqdb_end   = -1;

GeneWiseScore * gws = NULL;

Hscore * hs       = NULL;

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

char * gene_file        = NULL;
GeneFrequency21 * gf    = NULL;

char * matrix_file      = "BLOSUM62.bla";
CompMat * mat           = NULL;

char * gap_str          = "12";
int gap                 = 12;

char * ext_str          = "2";
int ext                 = 2;

char * length_of_N_str  = "10";
int length_of_N         = 10;

char * prob_in_rep      = "0.001";
double rep_prob         = 0.001  ;

char * codon_file       = NULL;
CodonTable * ct         = NULL;

char * output_file      = "-";
FILE * ofp              = NULL;


char * report_str       = NULL;
int report_stagger      = -1;


RandomModelDNA * rmd    = NULL;

char * subs_string      = "0.00001";
double subs_error       = 0.00001;

char * indel_string      = "0.00001";
double indel_error       = 0.00001;

char * cfreq_string      = "flat";
boolean model_codon      = FALSE;

char * splice_string     = "model";
boolean model_splice     = TRUE;

boolean flat_insert      = TRUE;

char * startend_string   = "default";
int startend             = TSM_default;


char * allN_string      = "0.9";
Probability allN        = 0.9;

char * null_string       = "syn";
boolean use_syn          = FALSE;

char * intron_string     = "tied";
boolean use_tied_model   = FALSE;

DBSearchImpl * dbsi      = NULL;
int alg                  = GWWRAP_623;
char * alg_str           = NULL;

DPRunImpl * dpri         = NULL;

int aln_alg              = GWWRAP_623;
char * aln_alg_str       = NULL;

int aln_number           = 50;
char * aln_number_str    = "50";

double  search_cutoff    = 20.00;
char * search_cutoff_str   = "20.00"; 

double evalue_search_cutoff = -1.0;
char * evalue_search_str = NULL;


boolean show_histogram   = TRUE;

boolean show_PackAln     = FALSE;
boolean show_AlnBlock    = FALSE;
boolean show_ace         = FALSE;
boolean show_gff         = FALSE;
boolean show_trans       = FALSE;
boolean show_pep         = FALSE;
boolean show_cdna        = FALSE;
boolean show_pretty      = FALSE;
boolean show_pretty_gene = FALSE;
boolean show_gene_plain  = FALSE;
boolean show_match_sum   = FALSE;
boolean show_para        = FALSE;
boolean show_overlap     = FALSE;

boolean pseudo           = FALSE;

char * main_block_str      = "50";
int main_block           = 50;

char * divide_str        = "//";

boolean complete         = FALSE; /* complete analysis or not */
boolean complete_trans   = FALSE; 
boolean complete_gene    = FALSE; 
boolean complete_cdna    = FALSE; 
boolean complete_ace    = FALSE; 
boolean show_diana       = FALSE;
boolean show_acehalf     = FALSE;
boolean show_embl        = FALSE;

Probability rnd_loop      = 0.99;
Probability cds_loop      = 0.97;
Probability rnd_to_model  = (1 - 0.99) / 3;
Probability link_loop     = 0.98;
Probability link_to_model = (1- 0.98) / 3;

AlnBlock * alb;
PackAln  * pal;

GenomicRegion * gr;
GenomicRegion * embl;

MatchSummarySet * mss;

RandomModel * rm;

GeneParameter21   * gpara;
GeneParser21Score * gps;
GeneParser4Score  * gp4s;
RandomCodonScore  * rcs;
RandomModelDNAScore * ids;

boolean show_output(void)
{
  int i,j,k;
  Genomic * gent;
  ThreeStateModel * temptsm;
  AlnBlock * alb;
  PackAln * pal;
  GenomicRegion * gr;
  GenomicRegion * cgr;
  MatchSummarySet * mss;
  Protein * ps;
  Protein * trans;
  cDNA * cdna;
  double bits;
  double aln_cutoff;
  boolean fitted_res = FALSE;
  AlnColumn * alt;


  /* sort by bit score first */

  sort_Hscore_by_score(hs);


  aln_cutoff = search_cutoff; /* for the moment */

  if( search_mode == PG_SEARCH_S2DB ) {
    if( hs->his == NULL || hs->his->total < 1000 ) {
	info("Cannot fit histogram to a db smaller than 1,000");
	fprintf(ofp,"[Warning: Can't fit histogram to a db smaller than 1,000]\n\n");
	fitted_res = FALSE;
	show_histogram = FALSE;
    } else {
      fit_Hscore_to_EVD(hs,20);
      fitted_res = TRUE;
    }
  }


  /* if you think this is scary, you are absolutely right! */

  if( alg == GWWRAP_6LITE && aln_alg != GWWRAP_6LITE) {
    /* memory leak here. Apologies */
    gpara->cses->cse[1] = codon_number_ComplexSequenceEval();
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

/*    fprintf(ofp,"Sequence %-24s %c %.2f\n",hs->ds[i]->target->name,hs->ds[i]->target->is_reversed == TRUE ? '-' : '+',bits); */

  }


  if( search_mode == PG_SEARCH_S2DB && show_histogram == TRUE ) {
    fprintf(ofp,"\n\n#Histogram\n");
    fprintf(ofp,"-----------------------------------------------------------------------\n");
    PrintASCIIHistogram(hs->his,ofp);
  }


  fprintf(ofp,"\n\n#Alignments\n");
  fprintf(ofp,"-----------------------------------------------------------------------\n");

  if( complete == TRUE ) {
    cgr = new_GenomicRegion(gen);
  }


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

    
    fprintf(ofp,"\n>Results for %s vs %s (%s) [%d]\n",hs->ds[i]->query->name,hs->ds[i]->target->name,hs->ds[i]->target->is_reversed == TRUE ? "reverse" : "forward",i );

    gent = get_Genomic_from_GenomicDB(gdb,hs->ds[i]->target);
    if( gent == NULL ) {
      warn("Unable to retrieve %s for alignment. Skipping...",hs->ds[i]->target->name);
      continue;
    }

    temptsm = indexed_ThreeStateModel_ThreeStateDB(tsmdb,hs->ds[i]->query);

    alb = AlnBlock_from_TSM_genewise_wrap(temptsm,gent,gpara,rmd,rmd,use_syn,aln_alg,allN,flat_insert,dpri,&pal,NULL);
   

    gr = new_GenomicRegion(gent);

    add_Genes_to_GenomicRegion_GeneWise(gr,gent->baseseq->offset,gent->baseseq->end,alb,NULL,pseudo,NULL);

    /* put in their seqname */
    for(j=0;j<gr->len;j++) {
      if( temptsm->accession != NULL ) 
	gr->gene[j]->seqname = stringalloc(temptsm->accession);
      else 
	gr->gene[j]->seqname = stringalloc(temptsm->name);
    }


    if( complete == TRUE ) {
      /* copy genes over into cgr region */
      for(j=0;j<gr->len;j++) {
	if( gr->gene[j]->bits > aln_cutoff ) 
	  add_Gene_to_GenomicRegion(cgr,hard_link_Gene(gr->gene[j]));
      }
    }
  
    if( use_single_pro == FALSE || pro == NULL) 
      mss = MatchSummarySet_from_AlnBlock_genewise(alb,temptsm->name,1,gent->baseseq);
    else
      mss = MatchSummarySet_from_AlnBlock_genewise(alb,pro->baseseq->name,pro->baseseq->offset,gent->baseseq);

    
    if( show_pretty == TRUE ) {

      fprintf(ofp,"\n%s output\nScore %4.2f bits over entire alignment.\nThis will be different from per-alignment scores. See manual for details\nFor computer parsable output, try %s -help or read the manual\n",program_name,Score2Bits(pal->score),program_name);
            
      if( alg == GWWRAP_2193L || alg == GWWRAP_2193) {
	fprintf(ofp,"Entrie alignment score contains unseen 'random' score segments\nYou should only use the per-alignments score printed below\nfor the bits score of the alignment\n\n");
      }
      
      if( use_syn == FALSE ) {
	fprintf(ofp,"Scores as bits over a flat simple random model\n\n");
      } else {
	fprintf(ofp,"Scores as bits over a synchronous coding model\n\n");
      }
      
      ps = pseudo_Protein_from_ThreeStateModel(temptsm);
      protgene_ascii_display(alb,ps->baseseq->seq,ps->baseseq->name,ps->baseseq->offset,gent,ct,15,main_block,TRUE,ofp);
      free_Protein(ps);
      
      fprintf(ofp,"%s\n",divide_str);
      
    }

    if( show_match_sum == TRUE ) {
      show_MatchSummary_genewise_header(ofp);
      show_MatchSummarySet_genewise(mss,ofp);
      fprintf(ofp,"%s\n",divide_str);
    }
    
    if( show_pretty_gene == TRUE ) {
      show_pretty_GenomicRegion(gr,0,ofp);
      fprintf(ofp,"%s\n",divide_str);
    }
    
    
    
    if( show_trans == TRUE ) {
      for(k=0;k<gr->len;k++) {
	if( gr->gene[k]->ispseudo == TRUE ) {
	  fprintf(ofp,"#Gene %d is a pseudo gene - no translation possible\n",k);
	} else {
	  trans = get_Protein_from_Translation(gr->gene[k]->transcript[0]->translation[0],ct);
	  write_fasta_Sequence(trans->baseseq,ofp);
	}
      } 
      fprintf(ofp,"%s\n",divide_str);
    }


    if( show_pep == TRUE ) {
      alt = alb->start;
      for(;alt != NULL;) {
	trans = Protein_from_GeneWise_AlnColumn(gent->baseseq,alt,1,&alt,ct,is_random_AlnColumn_genewise);
	if ( trans == NULL ) 
	  break;
	write_fasta_Sequence(trans->baseseq,ofp);
	free_Protein(trans);
      }
      fprintf(ofp,"%s\n",divide_str);
    }
    
    
    if( show_cdna == TRUE ) {
      for(k=0;k<gr->len;k++) {
	cdna = get_cDNA_from_Transcript(gr->gene[k]->transcript[0]);
	write_fasta_Sequence(cdna->baseseq,ofp);
      } 
      fprintf(ofp,"%s\n",divide_str);
    }
    
    
    if( show_ace == TRUE ) {
      show_ace_GenomicRegion(gr,gent->baseseq->name,ofp);
      fprintf(ofp,"%s\n",divide_str);
    }
    
    if( show_gff == TRUE ) {
      show_GFF_GenomicRegion(gr,gent->baseseq->name,"GeneWise",ofp);
      fprintf(ofp,"%s\n",divide_str);
    }
    
    if( show_gene_plain == TRUE ) {
      show_GenomicRegion(gr,ofp);
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

    alb = free_AlnBlock(alb);
    pal = free_PackAln(pal);
    mss = free_MatchSummarySet(mss);
    gr = free_GenomicRegion(gr);
    gent = free_Genomic(gent);
    temptsm = free_ThreeStateModel(temptsm);

  }

  if( complete == TRUE ) {

    fprintf(ofp,"\n\n#Complete Analysis\n");
    fprintf(ofp,"-------------------------------------------------------------\n\n");
    

    sort_GenomicRegion_absolute(cgr);
    if( complete_ace  == TRUE ) {
      show_ace_GenomicRegion(cgr,cgr->genomic->baseseq->name,ofp);
      fprintf(ofp,"%s\n",divide_str);
    }

    if( complete_gene == TRUE ) {
      show_pretty_GenomicRegion(cgr,0,ofp);
      fprintf(ofp,"%s\n",divide_str);
    }

    if( complete_cdna == TRUE ) {
      dump_transcripts_GenomicRegion(cgr,ofp);
      fprintf(ofp,"%s\n",divide_str);
    }

    if( complete_trans == TRUE ) {
      dump_translations_GenomicRegion(cgr,ct,ofp);
      fprintf(ofp,"%s\n",divide_str);
    }

    if( show_diana == TRUE ) {
      write_Diana_FT_GenomicRegion(cgr,ofp);
    }

    if( show_embl == TRUE ) {
      write_Embl_FT_GenomicRegion(cgr,ofp);
    }

    if( show_acehalf == TRUE ) {
      show_halfwise_GenomicRegion(cgr,cgr->genomic->baseseq->name,"HALFWISE","PFAM",1,"Pfam-Sanger",ofp);
    }

  }

  return TRUE;
}

boolean show_header(FILE * ofp)
{
  fprintf(ofp,"Wise2 - database searching mode\n");
  fprintf(ofp,"Program: %s version: %s released: %s\n",program_name,VERSION_NUMBER,RELEASE_DAY);
  fprintf(ofp,"This program is freely distributed under a Gnu Public License.\n");
  fprintf(ofp,"   See -version for more info on copyright\n");
  fprintf(ofp,"Bugs and credits to Ewan Birney <birney@sanger.ac.uk>\n");
  fprintf(ofp,"-----------------------------------------------------\n\n");
  fprintf(ofp,"Algorithm type:        GeneWise\n");
  fprintf(ofp,"Search algorithm used: %s\n",alg_str);
  fprintf(ofp,"Implementation:        %s\n",impl_string_DBSearchImpl(dbsi));
  fprintf(ofp,"Search mode:           %s\n",search_mode == PG_SEARCH_S2DB ? "Single protein vs genomic db" : search_mode == PG_SEARCH_DB2S ? "Single genomic vs protein db" : "Protein db vs genomic db");
  fprintf(ofp,"Protein info from:     %s\n",protein_file);
  fprintf(ofp,"Dna info from:         %s\n",dna_seq_file);
  if( use_single_pro == TRUE || use_db_pro == TRUE )  {
    fprintf(ofp,"Comp Matrix:           %s\n",matrix_file);
    fprintf(ofp,"Gap open:              %d\n",gap);
    fprintf(ofp,"Gap extension:         %d\n",ext);
  }
  fprintf(ofp,"Start/End (protein)    %s\n",startend_string);
  fprintf(ofp,"Gene Paras:            %s\n",gene_file);
  fprintf(ofp,"Codon Table:           %s\n",codon_file);
  fprintf(ofp,"Subs error:            %2.2g\n",subs_error);
  fprintf(ofp,"Indel error:           %2.2g\n",indel_error);
  fprintf(ofp,"Model splice?          %s\n",splice_string);
  fprintf(ofp,"Model codon bias?      %s\n",cfreq_string);
  fprintf(ofp,"Model intron bias?     %s\n",intron_string);
  fprintf(ofp,"Null model             %s\n",null_string);
  fprintf(ofp,"Alignment Alg          %s\n",aln_alg_str);

  return TRUE;
}
  
boolean search_db(void)
{

  info("Starting search...");

  hs = Hscore_from_TSM_genewise(tsmdb,gdb,gpara,rmd,rmd,use_syn,alg,search_cutoff,allN,report_stagger,FALSE,flat_insert,dbsi);

  if( hs == NULL ) {
    return FALSE;
  }

  return TRUE;
}

boolean build_db_objects(void)
{

  gpara = GeneParameter21_wrap(gf,subs_error,indel_error,rmd,model_codon,model_splice,use_tied_model,ct,rnd_loop,cds_loop,rnd_to_model,link_loop,link_to_model);


  if( gpara == NULL ) {
    warn("Sorry - could not build gene parameters. Must be a bug of some sort");
    return FALSE;
  }

  if( prepare_ComplexSequenceEvalSet(gpara->cses) == FALSE ) {
    warn("Unable to prepare complexsequenceevalset");
    return FALSE;
  }

  if( alg == GWWRAP_6LITE ) {
    /* memory leak here. Apologies */
    gpara->cses->cse[1] = codon64_number_ComplexSequenceEval();
  }

  if( use_single_dna == TRUE ) {
    gdb = new_GenomicDB_from_single_seq(gen,gpara->cses,Probability2Score(rep_prob));
  } else {
    gdb = new_GenomicDB(sdb,gpara->cses,length_of_N,Probability2Score(rep_prob));
  }

  
  if( gdb == NULL ) {
    return FALSE;
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
    gen = read_fasta_file_Genomic(dna_seq_file,length_of_N);
    if( gen == NULL ) {
      warn("Could not open single dna sequence in %s",dna_seq_file);
      ret = FALSE;
    }
  } else {
    sdb = single_fasta_SequenceDB(dna_seq_file);
    
 
    if( sdb == NULL ) {
      warn("Could not build a sequence database on %s",dna_seq_file);
      ret = FALSE;
    }

    if( dna_seqdb_start != -1 && dna_seqdb_end != -1 ) {
      sdb->seq_start = dna_seqdb_start;
      sdb->seq_end   = dna_seqdb_end;
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

    /*tsm = read_HMMer_1_7_ascii_file(hmm_file);*/
    /*    tsm = Wise2_read_ThreeStateModel_from_hmmer1_file(protein_file);*/
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

      
      
      set_startend_policy_ThreeStateModel(tsm,startend,30,0.2);
      
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

    if( set_search_type_ThreeStateDB(tsmdb,startend_string) == FALSE) {
      warn("Unable to set global/local switch on threestatedb");
      ret = FALSE;
    }
    free_SequenceDB(psdb);
  }
  else {
    warn("No protein input file! Yikes!");
  }


  if( pro_seqdb_start != -1 && pro_seqdb_end != -1 ) {
    tsmdb->hmm_model_start = pro_seqdb_start;
    tsmdb->hmm_model_end   = pro_seqdb_end;
  }


  if( main_block_str != NULL ) {
    if( is_integer_string(main_block_str,&main_block) == FALSE ) {
      warn("Could not get maximum main_block number %s",main_block_str);
      ret = FALSE;
    }
  }
   
  if( alg_str != NULL ) {
    alg = gwrap_alg_type_from_string(alg_str);
  } else {
    alg_str = "623";
    alg = gwrap_alg_type_from_string(alg_str);
  }

  if( aln_alg_str != NULL ) {
    aln_alg = gwrap_alg_type_from_string(aln_alg_str);
  } else {
    if( use_single_pro == TRUE || use_db_pro == TRUE) {
      aln_alg_str = "623";
    } else {
      aln_alg_str = "623L";
    }
    aln_alg = gwrap_alg_type_from_string(aln_alg_str);
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

  
  if( strcmp(cfreq_string,"model") == 0 ) {
    model_codon = TRUE;
  } else if ( strcmp(cfreq_string,"flat") == 0 ) {
    model_codon = FALSE;
  } else {
    warn("Cannot interpret [%s] as a codon modelling parameter\n",cfreq_string);
    ret = FALSE;
  }
  

  if( strcmp(splice_string,"model") == 0 ) {
    model_splice = TRUE;
  } else if ( strcmp(splice_string,"flat") == 0 ) {
    model_splice = FALSE;
  } else {
    warn("Cannot interpret [%s] as a splice modelling parameter\n",splice_string);
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

  if( strcmp(intron_string,"model") == 0 ) {
    use_tied_model = FALSE;
  } else if ( strcmp(intron_string,"tied") == 0 ) {
    use_tied_model = TRUE;
  } else {
    warn("Cannot interpret [%s] as a intron tieing switch\n",intron_string);
    ret = FALSE;
  }



  if( (rm = default_RandomModel()) == NULL) {
    warn("Could not make default random model\n");
    ret = FALSE;
  }

  if( (gf = read_GeneFrequency21_file(gene_file)) == NULL) {
    ret = FALSE;
    warn("Could not read a GeneFrequency file in %s",gene_file);
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
  return ret;

}


void show_short_help(void)
{
  fprintf(stdout,"%s (%s)\n",program_name,VERSION_NUMBER);
  fprintf(stdout,"This program is freely distributed under a GPL. See -version for more info\n");
  fprintf(stdout,"Copyright (c) GRL limited: portions of the code are from separate copyrights\n\n");
  fprintf(stdout,"swisepg <protein-input> <dna-input> in fasta format\n");
  fprintf(stdout," Options. In any order, '-' as filename (for any input) means stdin\nDon't use stdin for databases, as on-the-fly indexing is used\n");
  fprintf(stdout," Protein type  [-protein,-prodb,-hmmer,-pfam] [default - protein]\n");
  fprintf(stdout," Dna type      [-dnas,-dnadb] [default - dnadb]\n");
  fprintf(stdout," Protein  [-s,-t,-g,-e,-m]\n HMM      [-hmmer,-hname]\n");
  fprintf(stdout," Model    [-codon,-gene,-cfreq,-splice,-subs,-indel,-intron,-null]\n Alg      [-kbyte,-alg,-aalg,-aln,-noh]\n");
  fprintf(stdout," Output   [-pretty,-genes,-para,-sum,-cdna,-trans,-ace,]\n");
  fprintf(stdout,"  ..cont  [-gff,-gener,-alb,-pal,-block,-divide]\n");
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
  fprintf(ofp,"  -pfam     pfam hmm library\n");
  fprintf(ofp,"  -pfam2    pfam old style model directory (2.1) \n");
  fprintf(ofp,"  -hmmer    single hmmer HMM (version 2 compatible)\n");
  fprintf(ofp,"  -pro_db_start start position in protein/hmm database\n");
  fprintf(ofp,"  -pro_db_end   end   position in protein/hmm database\n");
  fprintf(ofp,"DNA input type\n");
  fprintf(ofp,"  -dnadb    [default] dna fasta database\n");
  fprintf(ofp,"  -dnas     a single dna fasta sequence\n");
  fprintf(ofp,"  -dna_db_start start position in dna database\n");
  fprintf(ofp,"  -dna_db_end   end   position in dna database\n");
  fprintf(ofp,"Protein comparison options\n");
  fprintf(ofp,"  -gap      [%3d]  gap penalty\n",gap);
  fprintf(ofp,"  -ext      [%3d]  extension penalty\n",ext);
  fprintf(ofp,"  -matrix   [%s]  Comparison matrix\n",matrix_file);
  fprintf(ofp,"HMM options\n");
  fprintf(ofp,"  -hname           For single hmms, use this as the name, not filename\n");
  fprintf(ofp,"Gene Model options\n");
  fprintf(ofp,"  -init   [%s]  [default/global/local/wing/endbias] start-end policy\n",startend_string);
  fprintf(ofp,"  -codon  [%s]  Codon file\n",codon_file);
  fprintf(ofp,"  -gene   [%s]  Gene parameter file\n",gene_file);
  fprintf(ofp,"  -subs   [%2.2g] Substitution error rate\n",subs_error);
  fprintf(ofp,"  -indel  [%2.2g] Insertion/deletion error rate\n",indel_error);
  fprintf(ofp,"  -cfreq  [model/flat] Using codon bias or not?     [default flat]\n");
  fprintf(ofp,"  -splice [model/flat] Using splice model or GT/AG? [default model]\n");
  fprintf(ofp,"  -intron [model/tied] Use tied model for introns   [default tied]\n");
  fprintf(ofp,"  -null   [syn/flat]   Random Model as synchronous or flat [default syn]\n");
  fprintf(ofp,"  -insert [model/flat] Use protein insert model     [default flat]\n");

  fprintf(ofp,"Algorithm options\n");
  fprintf(ofp,"  -alg    [623/2193/]            Algorithm used for searching [default 623]\n");
  fprintf(ofp,"  -aalg   [623/623L/2193/2193L]  Algorithm used for alignment [default 623/623L]\n");
  fprintf(ofp,"  -cut    [%.2f]   Bits cutoff for reporting in search algorithm\n",search_cutoff);
  fprintf(ofp,"  -ecut   [n/a]    Evalue cutoff for single protein vs DNA searches.\n");
  fprintf(ofp,"  -aln    [%d]   Max number of alignments (even if above cut)\n",aln_number);
  fprintf(ofp,"  -alln   [%s]   Probability of matching a NNN codon (only for single HMM vs DNAdb)\n",allN_string);
  fprintf(ofp,"  -nohis           Don't show histogram on single protein/hmm vs DNA search\n");
  fprintf(ofp,"  -report [0]      Issue a report every x comparisons (default 0 comparisons)\n");
  fprintf(ofp,"Output options [default -pretty -para]\n");
  fprintf(ofp,"  -pretty          show pretty ascii output\n");
  fprintf(ofp,"  -genes           show gene structure\n");
  fprintf(ofp,"  -para            show parameters\n");
  fprintf(ofp,"  -sum             show summary output\n");
  fprintf(ofp,"  -cdna            show cDNA\n");
  fprintf(ofp,"  -trans           show protein translation\n");
  fprintf(ofp,"  -ace             ace file gene structure\n");
  fprintf(ofp,"  -gff             Gene Feature Format file\n");
  fprintf(ofp,"  -gener           raw gene structure\n");
  fprintf(ofp,"  -alb             show logical AlnBlock alignment\n");
  fprintf(ofp,"  -pal             show raw matrix alignment\n");
  fprintf(ofp,"  -block  [%s]     Length of main block in pretty output\n",main_block_str);
  fprintf(ofp,"  -divide [%s]     divide string for multiple outputs\n",divide_str);
  fprintf(ofp,"Output for complete analysis (only available for single dna seq vs proteindb)\n");
  fprintf(ofp,"  -pseudo          Mark genes with frameshifts as pseudogenes\n");
  fprintf(ofp,"  -ctrans          provide all translations\n");
  fprintf(ofp,"  -ccdna           provide all cdna\n");
  fprintf(ofp,"  -cgene           provide all gene structures\n");
  fprintf(ofp,"  -cace            provide all gene structures in ace format\n");
  fprintf(ofp,"  -cdiana          provide all gene structures in diana EMBL FT format\n");
  fprintf(ofp,"  -cembl           provide all gene structures in EMBL FT format\n");
  fprintf(ofp,"  -caceh           provide all gene structures in halfwise type ace format\n");

  show_help_DBSearchImpl(ofp);
  show_help_DPRunImpl(ofp);

  show_standard_options(ofp);
  fprintf(ofp,"\nSee WWW help at http://www.sanger.ac.uk/Software/Wise2/\n");
  exit(63);   
}

void show_version(FILE * ofp)
{
  fprintf(ofp,"%s\nVersion: %s\nReleased: %s\nCompiled: %s\n",program_name,VERSION_NUMBER,RELEASE_DAY,COMPILE_DATE);
  fprintf(ofp,"\nThis program is freely distributed under a Gnu Public License\n");
  fprintf(ofp,"The source code is copyright (c) GRL 1998 and others\n");
  fprintf(ofp,"There is no warranty, implied or otherwise on the performance of this program\n");
  fprintf(ofp,"For more information read the GNULICENSE file in the distribution\n\n");
  fprintf(ofp,"Credits: Ewan Birney <birney@sanger.ac.uk> wrote the core code.\n");
  fprintf(ofp,"         Portions of this code was from HMMer2, written by Sean Eddy\n");
  exit(63);   
}


void build_defaults(void)
{
  gene_file = "human.gf";
  codon_file = "codon.table";
  matrix_file = "BLOSUM62.bla";
  

}

int main(int argc,char ** argv) 
{
  int i;
  char * temp;

  build_defaults();

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

  if( (temp = strip_out_assigned_argument(&argc,argv,"gene")) != NULL )
    gene_file = temp;

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

  if( (temp = strip_out_assigned_argument(&argc,argv,"alln")) != NULL )
    allN_string = temp;

  if( (temp = strip_out_assigned_argument(&argc,argv,"cfreq")) != NULL )
    cfreq_string = temp;

  if( (temp = strip_out_assigned_argument(&argc,argv,"splice")) != NULL )
    splice_string = temp;

  if( (temp = strip_out_assigned_argument(&argc,argv,"init")) != NULL )
    startend_string = temp;

  if( (temp = strip_out_assigned_argument(&argc,argv,"null")) != NULL )
    null_string = temp;

  if( (temp = strip_out_assigned_argument(&argc,argv,"intron")) != NULL )
    intron_string = temp;

  if( (temp = strip_out_assigned_argument(&argc,argv,"insert")) != NULL ) {
    if( strcmp(temp,"flat") == 0 ) {
      flat_insert = TRUE;
    } else {
      flat_insert = FALSE;
    }
  }
    

  if( (temp = strip_out_assigned_argument(&argc,argv,"report")) != NULL )
    report_str = temp;

  pseudo = strip_out_boolean_argument(&argc,argv,"pseudo");

  if( (strip_out_boolean_argument(&argc,argv,"dnas")) == TRUE )
    use_single_dna = TRUE;

  if( (strip_out_boolean_argument(&argc,argv,"dnadb")) == TRUE )
    use_single_dna = FALSE;

  strip_out_integer_argument(&argc,argv,"dna_db_start",&dna_seqdb_start);
  strip_out_integer_argument(&argc,argv,"dna_db_end",&dna_seqdb_end);

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

  strip_out_integer_argument(&argc,argv,"pro_db_start",&pro_seqdb_start);
  strip_out_integer_argument(&argc,argv,"pro_db_end",&pro_seqdb_end);



  if( (strip_out_boolean_argument(&argc,argv,"intie")) == TRUE )
    use_tied_model = TRUE;

  if( (temp = strip_out_assigned_argument(&argc,argv,"hname")) != NULL )
    hmm_name = temp;

  if( (strip_out_boolean_argument(&argc,argv,"nohis")) != FALSE )
    show_histogram = FALSE;

  if( (strip_out_boolean_argument(&argc,argv,"pretty")) != FALSE )
    show_pretty = TRUE;

  if( (strip_out_boolean_argument(&argc,argv,"gff")) != FALSE )
    show_gff = TRUE;


  if( (strip_out_boolean_argument(&argc,argv,"genes")) != FALSE )
    show_pretty_gene = TRUE;

  if( (strip_out_boolean_argument(&argc,argv,"para")) != FALSE )
    show_para = TRUE;

  if( (strip_out_boolean_argument(&argc,argv,"trans")) != FALSE )
    show_trans = TRUE;

  if( (strip_out_boolean_argument(&argc,argv,"pep")) != FALSE )
    show_pep = TRUE;

  if( (strip_out_boolean_argument(&argc,argv,"cdna")) != FALSE )
    show_cdna = TRUE;

  if( (strip_out_boolean_argument(&argc,argv,"sum")) != FALSE )
    show_match_sum = TRUE;

  if( (strip_out_boolean_argument(&argc,argv,"alb")) != FALSE )
    show_AlnBlock = TRUE;

  if( (strip_out_boolean_argument(&argc,argv,"ace")) != FALSE )
    show_ace = TRUE;

  if( (strip_out_boolean_argument(&argc,argv,"pal")) != FALSE )
    show_PackAln = TRUE;

  if( (strip_out_boolean_argument(&argc,argv,"gener")) != FALSE )
    show_gene_plain = TRUE;

  if( (strip_out_boolean_argument(&argc,argv,"over")) != FALSE )
    show_overlap = TRUE;

  if( (temp = strip_out_assigned_argument(&argc,argv,"divide")) != NULL )
    divide_str = temp;

  if( (temp = strip_out_assigned_argument(&argc,argv,"block")) != NULL )
    main_block_str = temp;


  complete_trans = strip_out_boolean_argument(&argc,argv,"ctrans") ;
  complete_gene  = strip_out_boolean_argument(&argc,argv,"cgene") ;
  complete_ace   = strip_out_boolean_argument(&argc,argv,"cace") ;
  complete_cdna  = strip_out_boolean_argument(&argc,argv,"ccdna"); 

  if( (strip_out_boolean_argument(&argc,argv,"cdiana")) != FALSE )
    show_diana = TRUE;

  if( (strip_out_boolean_argument(&argc,argv,"cembl")) != FALSE )
    show_embl = TRUE;

  if( (strip_out_boolean_argument(&argc,argv,"caceh")) != FALSE )
    show_acehalf = TRUE;


  dbsi = new_DBSearchImpl_from_argv(&argc,argv);
  
  dpri = new_DPRunImpl_from_argv(&argc,argv);

  assert(dpri);
  assert(dbsi);

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

  if( show_gff == FALSE && show_overlap == FALSE && show_pretty_gene == FALSE && show_match_sum == FALSE && show_ace == FALSE && show_gene_plain == FALSE && show_pretty == FALSE && show_AlnBlock == FALSE && show_PackAln == FALSE && show_pep == FALSE) {
    show_pretty = TRUE;
    show_para = TRUE;
  }

  if( use_db_pro == FALSE && use_single_pro == FALSE && use_tsm == FALSE && use_pfam1 == FALSE && use_pfam2 == FALSE) {
    use_single_pro = TRUE;
  }

  if( use_single_pro == TRUE || use_tsm == TRUE ) {
    if( use_single_dna == TRUE ) 
      fatal("one on one search. Shouldn't you use genewise?");
    search_mode = PG_SEARCH_S2DB;
  } else {
    if( use_single_dna == TRUE ) 
      search_mode = PG_SEARCH_DB2S;
    else 
      search_mode = PG_SEARCH_DB2DB;
  }

  if( evalue_search_str != NULL && search_mode != PG_SEARCH_S2DB ) {
    fatal("Trying to set a evalue cutoff on a non evalue based search. you can only use evalues in a protein HMM vs DNA database search (sorry!)");
  }


  if( complete_trans == TRUE || complete_gene == TRUE || complete_cdna == TRUE || complete_ace == TRUE || show_diana == TRUE || show_acehalf == TRUE || show_embl == TRUE) {
    if( search_mode != PG_SEARCH_DB2S ) {
      fatal("Complete analysis is only available on single dna vs proteindb comparisons - did you use -dnas?");
    }
    complete = TRUE;
  }

    
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

  return 0;
}













