#include "est_evidence.h"
#include "genomewise9.h"
#include "geneutil.h"
#include "version.h"



char * program_name = "genomewise";


void debug_genomewise(AlnBlock * alb,GenomeEvidenceSet * ges,CodonTable * ct,Sequence * gen,FILE * ofp);

void show_utr_exon_genomewise(AlnBlock * alb,FILE * ofp);

void show_version(FILE * ofp)
{
  fprintf(ofp,"%s\nVersion: %s\nReleased: %s\nCompiled: %s\n",program_name,VERSION_NUMBER,RELEASE_DAY,COMPILE_DATE);
  fprintf(ofp,"\nThis program is freely distributed under a Gnu Public License\n");
  fprintf(ofp,"The source code is copyright (c) EMBL and others\n");
  fprintf(ofp,"There is no warranty, implied or otherwise on the performance of this program\n");
  fprintf(ofp,"For more information read the GNULICENSE file in the distribution\n\n");
  fprintf(ofp,"Credits: Ewan Birney <birney@ebi.ac.uk>\n");
  exit(63);   
}

void show_help(FILE * ofp)
{
  fprintf(ofp,"%s genomic-fasta-file evidence-file\n",program_name);
  fprintf(ofp,"   ** Genomewise is designed to work with the Ensembl EST build system\n");
  fprintf(ofp,"   ** Although you can reuse it directly, alot of the magic occurs in \n");
  fprintf(ofp,"   ** the Ensembl Runnable/RunnableDB system behind this. see www.ensembl.org \n\n");

  fprintf(ofp,"   evidence file should have exon,cds and indel lines separated by //\n");
  fprintf(ofp,"   between predictions (multiple predictions ok)\n");
  fprintf(ofp,"      exon start end -- means exon prediction, no phase restriction\n");
  fprintf(ofp,"      cds  start end phase -- means exon prediction, only in that phase\n");
  fprintf(ofp,"      indel start end -- allow frameshifting in this area\n");
  fprintf(ofp,"eg - \n");
  fprintf(ofp,"exon 120 340\n");
  fprintf(ofp,"exon 560 591\n");
  fprintf(ofp,"//\n");
  fprintf(ofp,"cds  12  56 0\n");
  fprintf(ofp,"cds  70  80 1\n");
  fprintf(ofp,"\n\nOPTIONS (can occur anywhere on the command line\n");
  fprintf(ofp,"Scoring\n");
  fprintf(ofp,"   -start  <number> no start codon penalty 30\n");
  fprintf(ofp,"   -stop   <number> no stop codon penalty 200\n");
  fprintf(ofp,"   -gene   <number> new gene cost 5000\n");
  fprintf(ofp,"   -switch <number> evidence switch cost 100\n");
  fprintf(ofp,"   -smell  <number> smell space used for out-phase splice sites 8\n");
  fprintf(ofp,"Output\n");
  fprintf(ofp,"   -[no]genes       show gene structure (default yes)\n");
  fprintf(ofp,"   -[no]geneutr     show gene structure with utrs (default yes)\n");
  fprintf(ofp,"   -[no]trans       show protein translation (default yes)\n");
  fprintf(ofp,"   -[no]gff         show gff (default yes)\n");
  fprintf(ofp,"   -[no]alb         show aln block format (default no)\n");
  fprintf(ofp,"   -[no]debug       show debug format (default no)\n");
  fprintf(ofp,"   -kbyte  <number> number of kbytes available for main memory build (10,000 default)\n");

  show_help_DPRunImpl(ofp);

  show_standard_options(ofp);
}


double id(int i)
{
  return (double)i;
}


int main(int argc,char ** argv)
{
  Sequence   * gen;
  Genomic    * genomic;
  CodonTable * ct = NULL;
  GenomeEvidenceSet * ges = NULL;
  RandomCodonScore * rcs;
  FILE * ifp = NULL;
  ComplexSequence * cs = NULL;
  ComplexSequenceEvalSet * cses = NULL;
  AlnBlock * alb;
  PackAln * pal;
  GenomicRegion * gr;
  int i;
  Protein * trans;
  cDNA    * cdna;
  int kbyte                = 10000;
  int stop_codon_pen  = 200;
  int start_codon_pen = 30;
  int new_gene        = 5000;
  int switch_cost     = 100;
  int smell           = 8;
  DPRunImpl * dpri = NULL;
    
  EstEvidence * est;

  boolean show_trans = TRUE;
  boolean show_cdna  = FALSE;
  boolean show_genes = TRUE;
  boolean show_alb   = FALSE;
  boolean show_pal   = FALSE;
  boolean show_gff   = TRUE;
  boolean show_debug = FALSE;
  boolean show_geneu = TRUE;
  char * divide_string = "//";

  strip_out_boolean_def_argument(&argc,argv,"geneutr",&show_geneu);
  strip_out_boolean_def_argument(&argc,argv,"genes",&show_genes);
  strip_out_boolean_def_argument(&argc,argv,"trans",&show_trans);
  strip_out_boolean_def_argument(&argc,argv,"gff",&show_gff);
  strip_out_boolean_def_argument(&argc,argv,"alb",&show_alb);
  strip_out_boolean_def_argument(&argc,argv,"pal",&show_pal);
  strip_out_boolean_def_argument(&argc,argv,"debug",&show_debug);
  strip_out_boolean_def_argument(&argc,argv,"cdna",&show_cdna);
  strip_out_integer_argument(&argc,argv,"stop",&stop_codon_pen);
  strip_out_integer_argument(&argc,argv,"start",&start_codon_pen);
  strip_out_integer_argument(&argc,argv,"gene",&new_gene);
  strip_out_integer_argument(&argc,argv,"switch",&switch_cost);
  strip_out_integer_argument(&argc,argv,"smell",&smell);
  
  dpri = new_DPRunImpl_from_argv(&argc,argv);
  if( dpri == NULL ) {
    fatal("Unable to build DPRun implementation. Bad arguments");
  }


  strip_out_standard_options(&argc,argv,show_help,show_version);
  if( argc != 3 ) {
    show_help(stdout);
    exit(12);
  }

    

  ct  = read_CodonTable_file("codon.table");
  gen = read_fasta_file_Sequence(argv[1]);
  ifp = openfile(argv[2],"r");
  ges = read_est_evidence(ifp,ct);

  for(i=0;i<ges->len;i++) {
    est = (EstEvidence *) ges->geu[i]->data;
    est->in_smell = smell;
  }


  rcs= RandomCodonScore_alloc();
  for(i=0;i<125;i++) {
    if( is_stop_codon(i,ct) ) {
      rcs->codon[i] = -1000000;
    } else {
      rcs->codon[i] = 0;
    }
    /*    fprintf(stderr,"Got %d for %d\n",rcs->codon[i],i); */
  }

 

  cses = default_genomic_ComplexSequenceEvalSet();
  cs   = new_ComplexSequence(gen,cses);

 
  pal  = PackAln_bestmemory_GenomeWise9(ges,cs,-switch_cost,-new_gene,-start_codon_pen,-stop_codon_pen,rcs,NULL,dpri);
  alb  = convert_PackAln_to_AlnBlock_GenomeWise9(pal);


  genomic = Genomic_from_Sequence(gen);
  gr = new_GenomicRegion(genomic);

  add_Genes_to_GenomicRegion_GeneWise(gr,1,gen->len,alb,gen->name,0,NULL);

  if( show_genes ) {
    show_pretty_GenomicRegion(gr,0,stdout);
    fprintf(stdout,"%s\n",divide_string);
  }

  if( show_gff ) {
    show_GFF_GenomicRegion(gr,gen->name,"genomwise",stdout);
    fprintf(stdout,"%s\n",divide_string);
  }

  if( show_trans ) {
    for(i=0;i<gr->len;i++) {
      if( gr->gene[i]->ispseudo == TRUE ) {
	fprintf(stdout,"#Gene %d is a pseudo gene - no translation possible\n",i);
      } else {
	trans = get_Protein_from_Translation(gr->gene[i]->transcript[0]->translation[0],ct);
	write_fasta_Sequence(trans->baseseq,stdout);
      }
    } 
    fprintf(stdout,"%s\n",divide_string);
  }

  if( show_cdna ) {
    for(i=0;i<gr->len;i++) {
      cdna = get_cDNA_from_Transcript(gr->gene[i]->transcript[0]);
      write_fasta_Sequence(cdna->baseseq,stdout);
    } 
    fprintf(stdout,"%s\n",divide_string);
  }

  if( show_geneu ) {
    show_utr_exon_genomewise(alb,stdout);
    fprintf(stdout,"%s\n",divide_string);
  }

  if( show_alb ) {
    mapped_ascii_AlnBlock(alb,id,1,stdout);
    fprintf(stdout,"%s\n",divide_string);
  }

  if( show_debug ) {
    debug_genomewise(alb,ges,ct,gen,stdout);
    fprintf(stdout,"%s\n",divide_string);
  }
    
  if( show_pal ) {
    show_simple_PackAln(pal,stdout);
    fprintf(stdout,"%s\n",divide_string);
  }

  return 0;
}


void debug_genomewise(AlnBlock * alb,GenomeEvidenceSet * ges,CodonTable * ct,Sequence * gen,FILE * ofp)
{
  AlnColumn *alc;
  int cstart;


  for(alc=alb->start;alc != NULL;alc = alc->next ) {
    fprintf(ofp,"%4d %12s %12s [%3d][%5d %5d] ",alc->alu[1]->score[0],alc->alu[0]->text_label,alc->alu[1]->text_label,alc->alu[0]->start,alc->alu[1]->start+1,alc->alu[1]->end);
    if( strstartcmp(alc->alu[1]->text_label,"CODON") == 0 ) { 
      cstart = alc->alu[1]->start+1;
      fprintf(ofp,"%c%c%c  %c\n",gen->seq[cstart],gen->seq[cstart+1],gen->seq[cstart+2],aminoacid_from_seq(ct,gen->seq+cstart));
    } else {
      fprintf(ofp,"\n");
    }
  }
    
}

#define GW_EXON_TYPE_UTR5 45
#define GW_EXON_TYPE_CDS  46
#define GW_EXON_TYPE_UTR3 47
#define GW_EXON_TYPE_NONE 48

int  exon_type_AlnColumn_genomewise(AlnColumn * alc) 
{
  if( strcmp(alc->alu[1]->text_label,"CODON") == 0 ) {
    return GW_EXON_TYPE_CDS;
  }

  if( strcmp(alc->alu[1]->text_label,"UTR5") == 0 ) {
    return GW_EXON_TYPE_UTR5;
  }

  if( strcmp(alc->alu[1]->text_label,"UTR3") == 0 || strcmp(alc->alu[1]->text_label,"STOP_CODON") == 0) {
    return GW_EXON_TYPE_UTR3;
  }

  return GW_EXON_TYPE_NONE;
}


void show_utr_exon_genomewise(AlnBlock * alb,FILE * ofp)
{
  AlnColumn * alc;
  int exon_start;
  int exon_end;
  int is_start;
  int phase;
  int endphase;
  int is_3ss;

  for(alc=alb->start;alc != NULL;) {
    /* find the first exon */
    for(;alc != NULL && exon_type_AlnColumn_genomewise(alc) != GW_EXON_TYPE_UTR5;alc = alc->next)
      ;
    if( alc == NULL ) {
      break;
    }
    fprintf(ofp,"Gene\n");

    if( alc != NULL && exon_type_AlnColumn_genomewise(alc) == GW_EXON_TYPE_UTR5 ) {
      while( alc != NULL ) { /* while loop goes over all 5UTRs */
	exon_start = alc->alu[1]->start+2;
	for(;alc != NULL && exon_type_AlnColumn_genomewise(alc) == GW_EXON_TYPE_UTR5;alc = alc->next ) {
	  ;
	} 
        /*	fprintf(stderr,"Broken out with %s\n",alc->alu[1]->text_label); */

	if( strcmp(alc->alu[1]->text_label,"UTR5_INTRON") ==0 ) {
	  /* ntron. should be +2-1 at the end of this, goes to 1*/
	  fprintf(ofp,"  utr5 %d %d\n",exon_start,alc->alu[1]->start+1);
	  /* now loop through the intron */
	  for(;alc != NULL && strcmp(alc->alu[1]->text_label,"UTR5_INTRON") == 0;alc = alc->next ) {
	    ;
	  }
	  if( alc == NULL || exon_type_AlnColumn_genomewise(alc) != GW_EXON_TYPE_UTR5 ) {
	    break; /* while loop */
	  } else{
	    continue; /* another utr5 exon */
	  }
	} else {
	  /* print this guy and break */
	  fprintf(ofp,"  utr5 %d %d\n",exon_start,alc->alu[1]->start+1);
	  break;
	}
      }
    }


    /* we now should be at a CDS column */

    if( alc != NULL && exon_type_AlnColumn_genomewise(alc) == GW_EXON_TYPE_CDS ) {
      is_start = 1;
      while( alc != NULL ) { /* while loop goes over all 5UTRs */
	/* fprintf(stderr,"Entering codoing loop with %s\n",alc->alu[1]->text_label); */

	exon_start = alc->alu[1]->start+2;
	if( strstr(alc->alu[1]->text_label,"3SS") != NULL ) {
	  is_3ss = 1;
	  if( strstr(alc->alu[1]->text_label,"1") != NULL ) {
	    phase = 1;
	  } else if ( strstr(alc->alu[1]->text_label,"2") != NULL ) {
	    phase = 2;
	  } else {
	    phase = 0;
	  }
	  alc = alc->next;
	} else {
	  is_3ss = 0;
	  phase = 0;
	}

	
	if( phase == 1 ) {
	  exon_start += 3;
	} else if ( phase == 2) {
	  exon_start += 3;
	} else if ( is_3ss ) {
	  /* phase 0 and spliced needs adjusting */
	  exon_start += 3;
	} 

	  

	for(;alc != NULL && exon_type_AlnColumn_genomewise(alc) == GW_EXON_TYPE_CDS ;alc = alc->next ) {
	  ;
	} 

	if( strstr(alc->alu[1]->text_label,"5SS") != NULL ) {

	  exon_end = alc->alu[1]->start+1;

	  if( strstr(alc->alu[1]->text_label,"1") != NULL ) {
	    endphase = 1;
	  } else if ( strstr(alc->alu[1]->text_label,"2") != NULL ) {
	    endphase = 2;
	  } else {
	    endphase = 0;
	  }

	  if( endphase == 1 ) {
	    exon_end += 1;
	  } else if( endphase == 2 ) {
	    exon_end += 2;
	  } /* no change for phase 0 */
	  
	  /* intron. should be +1-1 at the end of this, goes to 0*/
	  fprintf(ofp,"  cds %d %d phase %d\n",exon_start,exon_end,phase);

	  /* now loop through the intron */
	  for(alc= alc->next;alc != NULL && strstr(alc->alu[1]->text_label,"3SS") == NULL;alc = alc->next ) {
	    ;
	  }
	  if( alc == NULL || strstr(alc->alu[1]->text_label,"3SS") == NULL ) {
	    break; /* while loop */
	  } else{
	    continue; /* another cds exon */
	  }
	} else {
	  fprintf(ofp,"  cds %d %d phase %d\n",exon_start,alc->alu[1]->start+1,phase);
	  break;
	}
      }
    }

   
    if( alc != NULL && exon_type_AlnColumn_genomewise(alc) == GW_EXON_TYPE_UTR3 ) {
      is_start = 1;
      while( alc != NULL ) { /* while loop goes over all 5UTRs */
	exon_start = alc->alu[1]->start+2;
	for(;alc != NULL && exon_type_AlnColumn_genomewise(alc) == GW_EXON_TYPE_UTR3 ;alc = alc->next ) {
	  ;
	} 


	if( strstr(alc->alu[1]->text_label,"INTRON") != NULL ) {
	  /* intron. should be +2-1 at the end of this, goes to 1*/
	  fprintf(ofp,"  utr3 %d %d\n",exon_start,alc->alu[1]->start+1);
	  /* now loop through the intron */
	  for(;alc != NULL && strstr(alc->alu[1]->text_label,"INTRON") != NULL;alc = alc->next ) {
	    ;
	  }
	  if( alc == NULL || exon_type_AlnColumn_genomewise(alc) != GW_EXON_TYPE_UTR3 ) {
	    break; /* while loop */
	  } else{
	    continue; /* another utr5 exon */
	  }
	} else {
	  fprintf(ofp,"  utr3 %d %d\n",exon_start,alc->alu[1]->start+1);
	  break;
	}
      }
    }

    fprintf(ofp,"End\n");
    /* back to next gene */

  }


}
  






