#include "est_evidence.h"
#include "genomewise9.h"
#include "geneutil.h"
#include "geneoutput.h"
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

  show_help_GeneOutputPara(ofp);

  show_help_DPRunImpl(ofp);



  show_standard_options(ofp);
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


  int kbyte                = 10000;
  int stop_codon_pen  = 200;
  int start_codon_pen = 30;
  int new_gene        = 5000;
  int switch_cost     = 100;
  int smell           = 8;
  DPRunImpl * dpri = NULL;
  GeneOutputPara * geneout = NULL;
  GeneOutputData data;

  EstEvidence * est;

  boolean show_debug = FALSE;

  char * divide_string = "//";

  strip_out_boolean_def_argument(&argc,argv,"debug",&show_debug);

  strip_out_integer_argument(&argc,argv,"stop",&stop_codon_pen);
  strip_out_integer_argument(&argc,argv,"start",&start_codon_pen);
  strip_out_integer_argument(&argc,argv,"gene",&new_gene);
  strip_out_integer_argument(&argc,argv,"switch",&switch_cost);
  strip_out_integer_argument(&argc,argv,"smell",&smell);
  
  dpri = new_DPRunImpl_from_argv(&argc,argv);

  if( dpri == NULL ) {
    fatal("Unable to build DPRun implementation. Bad arguments");
  }

  geneout = new_GeneOutputPara_from_argv(&argc,argv);
  assert(geneout);


  strip_out_standard_options(&argc,argv,show_help,show_version);
  if( argc != 3 ) {
    show_help(stdout);
    exit(12);
  }

    

  ct  = read_CodonTable_file("codon.table");
  
  /*  fprintf(stderr,"Codon table is %d\n",ct);*/

  gen = read_fasta_file_Sequence(argv[1]);
  ifp = openfile(argv[2],"r");
  ges = read_est_evidence(ifp,ct);

  for(i=0;i<ges->len;i++) {
    est = (EstEvidence *) (ges->geu[i]->data);
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

  data.pal = pal;
  data.alb = alb;
  data.gr  = gr;
  data.gen = genomic;
  data.ct  = ct;

  show_GeneOutput(&data,geneout,stdout);

  if( show_debug ) {
    debug_genomewise(alb,ges,ct,gen,stdout);
    fprintf(stdout,"%s\n",geneout->divide_string);
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
      for (cstart = alc->alu[1]->start+1; cstart <= alc->alu[1]->end; cstart++) {
        fprintf(ofp,"%c",gen->seq[cstart]);
      }
      fprintf(ofp,"\n");
    }
  }
    
}

  






