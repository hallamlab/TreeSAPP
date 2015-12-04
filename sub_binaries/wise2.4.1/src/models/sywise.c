#include "sywise20.h"
#include "version.h"
#include "genestats.h"
#include "geneutil.h"
#include "standardout.h"


char * program_name = "sywise";

char * codon_table = "codon.table";

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
  fprintf(ofp,"%s pairwise-alignment-as-fasta\n",program_name);

  show_help_GeneModelParam(ofp);

  show_help_ShowGenomicRegionOptions(ofp);

  show_help_DPRunImpl(ofp);

  show_help_StandardOutputOptions(ofp);

  show_standard_options(ofp);
}

PairBaseModelScore * make_PairBaseModelScore(void)
{
  PairBaseModelScore * out;
  PairBaseModel * m;

  m = simple_PairBaseModel(0.6,0.4,0.9);

  out= new_PairBaseModelScore(m);

  free_PairBaseModel(m);

/*
  out = zero_PairBaseModelScore();
*/

  show_PairBaseModelScore(out,stdout);
  return out;
}

CodonTable * ct;


PairBaseCodonModelScore * make_PairBaseCodonModelScore(CompProb * cp)
{
  PairBaseCodonModelScore * out;
  PairBaseCodonModel * m;
  FILE * ifp;


  if( ct == NULL ) {
    fatal("Could not read codon table");
  }

  assert(cp);

/*
  m = very_simple_PairBaseCodonModel(0.8,1.0/20.0,0.2,0.2,ct); 
*/

  m = make_flat_PairBaseCodonModel(cp,0.1,0.01,ct);


  
  diagonal_tweak_PairBaseCodonModel(m,2.0,16.0,1.0/4.0);
  

  out = new_PairBaseCodonModelScore(m);

/*  flatten_diagonal_PairBaseCodonModelScore(out,ct); */

  free_PairBaseCodonModel(m);

  /*  show_PairBaseCodonModelScore(out,ct,stdout); */

  return out;
}


void show_score_sequence(AlnBlock * alb,PairBaseSeq * pbs,PairBaseModelScore * m,FILE * ofp)
{
  AlnColumn * alc;
  char seq1[20];
  char seq2[20];
  int match_score = 0;


  for(alc=alb->start;alc != NULL; alc=alc->next) {
    if( strcmp(alc->alu[1]->text_label,"CODON") == 0 ) {
      seq1[0] = char_for_base(anchor_base_from_pairbase(pbs->seq[alc->alu[1]->start+1]));
      seq1[1] = char_for_base(anchor_base_from_pairbase(pbs->seq[alc->alu[1]->start+2]));
      seq1[2] = char_for_base(anchor_base_from_pairbase(pbs->seq[alc->alu[1]->start+3]));

      seq2[0] = char_for_base(informant_base_from_pairbase(pbs->seq[alc->alu[1]->start+1]));
      seq2[1] = char_for_base(informant_base_from_pairbase(pbs->seq[alc->alu[1]->start+2]));
      seq2[2] = char_for_base(informant_base_from_pairbase(pbs->seq[alc->alu[1]->start+3]));

      seq1[3] = '\0';
      seq2[3] = '\0';

      match_score = m->base[pbs->seq[alc->alu[1]->start+1]];
      match_score += m->base[pbs->seq[alc->alu[1]->start+2]];
      match_score += m->base[pbs->seq[alc->alu[1]->start+3]];

      fprintf(ofp,"%s %5d %5d %s [%c] vs %s [%c] %-.2f %-.2f\n",alc->alu[1]->text_label,alc->alu[0]->start+1,alc->alu[1]->start+1,seq1,aminoacid_from_seq(ct,seq1),seq2,aminoacid_from_seq(ct,seq2),Score2Bits(alc->alu[0]->score[0]),Score2Bits(match_score));
    }
    if( strstr(alc->alu[1]->text_label,"SS") != NULL ) {

      seq1[0] = char_for_base(anchor_base_from_pairbase(pbs->seq[alc->alu[1]->start-3]));
      seq1[1] = char_for_base(anchor_base_from_pairbase(pbs->seq[alc->alu[1]->start-2]));
      seq1[2] = char_for_base(anchor_base_from_pairbase(pbs->seq[alc->alu[1]->start-1]));
      seq1[3] = char_for_base(anchor_base_from_pairbase(pbs->seq[alc->alu[1]->start]));
      seq1[4] = char_for_base(anchor_base_from_pairbase(pbs->seq[alc->alu[1]->start+1]));
      seq1[5] = char_for_base(anchor_base_from_pairbase(pbs->seq[alc->alu[1]->start+2]));
      seq1[6] = char_for_base(anchor_base_from_pairbase(pbs->seq[alc->alu[1]->start+3]));
      seq1[7] = char_for_base(anchor_base_from_pairbase(pbs->seq[alc->alu[1]->start+4]));
      seq1[8] = char_for_base(anchor_base_from_pairbase(pbs->seq[alc->alu[1]->start+5]));
      seq1[9] = char_for_base(anchor_base_from_pairbase(pbs->seq[alc->alu[1]->start+6]));
      seq1[10] = char_for_base(anchor_base_from_pairbase(pbs->seq[alc->alu[1]->start+7]));

      seq2[0] = char_for_base(informant_base_from_pairbase(pbs->seq[alc->alu[1]->start-3]));
      seq2[1] = char_for_base(informant_base_from_pairbase(pbs->seq[alc->alu[1]->start-2]));
      seq2[2] = char_for_base(informant_base_from_pairbase(pbs->seq[alc->alu[1]->start-1]));
      seq2[3] = char_for_base(informant_base_from_pairbase(pbs->seq[alc->alu[1]->start]));
      seq2[4] = char_for_base(informant_base_from_pairbase(pbs->seq[alc->alu[1]->start+1]));
      seq2[5] = char_for_base(informant_base_from_pairbase(pbs->seq[alc->alu[1]->start+2]));
      seq2[6] = char_for_base(informant_base_from_pairbase(pbs->seq[alc->alu[1]->start+3]));
      seq2[7] = char_for_base(informant_base_from_pairbase(pbs->seq[alc->alu[1]->start+4]));
      seq2[8] = char_for_base(informant_base_from_pairbase(pbs->seq[alc->alu[1]->start+5]));
      seq2[9] = char_for_base(informant_base_from_pairbase(pbs->seq[alc->alu[1]->start+6]));
      seq2[10] = char_for_base(informant_base_from_pairbase(pbs->seq[alc->alu[1]->start+7]));

      seq1[11] = '\0';
      seq2[11] = '\0';


      fprintf(ofp,"%12s %5d %5d %.2f\n",alc->alu[1]->text_label,alc->alu[0]->start+1,alc->alu[1]->start+1,Score2Bits(alc->alu[0]->score[0]));
      fprintf(ofp,"     %s\n     %s\n",seq1,seq2);
    }
  }

}


double id(int i)
{
  return (double)i;
}

int main(int argc,char ** argv)
{
  int i;

  DPRunImpl * dpri = NULL;
  GeneModelParam * gmp = NULL;
  GeneModel * gm = NULL;

  FILE * ifp;
  SeqAlign   * al;
  PairBaseSeq * pbs;

  ComplexSequenceEval * splice5;
  ComplexSequenceEval * splice3;
  ComplexSequence * cseq;


  CompMat * score_mat;
  CompProb * comp_prob;
  RandomModel * rm;

  PairBaseCodonModelScore * codon_score;
  PairBaseModelScore* nonc_score;

  PairBaseCodonModelScore * start;
  PairBaseCodonModelScore * stop;


  SyExonScore * exonscore;

  PackAln * pal;
  AlnBlock * alb;

  Genomic * genomic;
  GenomicRegion * gr;
  GenomicRegion * gr2;
  Protein * trans;

  StandardOutputOptions * std_opt;
  ShowGenomicRegionOptions * sgro;
  
  char * dump_packaln = NULL;
  char * read_packaln = NULL;
  FILE * packifp = NULL;

  boolean show_trans    = 1;
  boolean show_gene_raw = 0;



  ct = read_CodonTable_file(codon_table);
/*
  score_mat = read_Blast_file_CompMat("blosum62.bla");
  comp_prob = CompProb_from_halfbit(score_mat);
*/
  rm = default_RandomModel();

  comp_prob = read_Blast_file_CompProb("wag85");

  fold_column_RandomModel_CompProb(comp_prob,rm);

  dpri = new_DPRunImpl_from_argv(&argc,argv);
  if( dpri == NULL ) {
    fatal("Unable to build DPRun implementation. Bad arguments");
  }

  gmp = new_GeneModelParam_from_argv(&argc,argv);

  std_opt = new_StandardOutputOptions_from_argv(&argc,argv);
  sgro = new_ShowGenomicRegionOptions_from_argv(&argc,argv);

  
  dump_packaln = strip_out_assigned_argument(&argc,argv,"dump");
  read_packaln = strip_out_assigned_argument(&argc,argv,"recover");

  strip_out_standard_options(&argc,argv,show_help,show_version);
  if( argc != 2 ) {
    show_help(stdout);
    exit(12);
  }

  
  if((gm=GeneModel_from_GeneModelParam(gmp)) == NULL ) {
    fatal("Could not build gene model");
  }

  codon_score = make_PairBaseCodonModelScore(comp_prob);
  nonc_score  = make_PairBaseModelScore();

  splice5 = ComplexSequenceEval_from_pwmDNAScore_splice(gm->splice5score);
  splice3 = ComplexSequenceEval_from_pwmDNAScore_splice(gm->splice3score);

  if((ifp = openfile(argv[1],"r")) == NULL ) {
    fatal("Could not open file %s",argv[1]);
  }

  al = read_fasta_SeqAlign(ifp);
  
  assert(al);
  assert(al->len == 2);
  assert(al->seq[0]->len > 0);
  assert(al->seq[1]->len > 0);

/*  write_fasta_SeqAlign(al,stdout);*/

  
  pbs = new_PairBaseSeq_SeqAlign(al);

  if( read_packaln == NULL ) {
    cseq = ComplexSequence_from_PairBaseSeq(pbs,splice5,splice3);    
  }

  start = make_start_PairBaseCodonModelScore(ct);
  stop  = make_stop_PairBaseCodonModelScore(ct);


/*  show_PairBaseCodonModelScore(stop,ct,stdout); */
  
/*
  for(i=0;i<pbs->anchor->len;i++) {
    printf("%3d  %c For %-6d %-6d %c Rev %-6d %-6d\n",i,pbs->anchor->seq[i],
	   CSEQ_PAIR_5SS(cseq,i),CSEQ_PAIR_3SS(cseq,i),
	   char_complement_base(pbs->anchor->seq[i]),
	   CSEQ_REV_PAIR_5SS(cseq,i),CSEQ_REV_PAIR_3SS(cseq,i));
  }
*/


  /*  show_ComplexSequence(cseq,stdout);

  */


  exonscore = SyExonScore_flat_model(100,150,0.1,1.0);
  /*
  for(i=0;i<cseq->length;i++) {
    fprintf(stdout,"%d PairSeq is %d score %d\n",i,CSEQ_PAIR_PAIRBASE(cseq,i),nonc_score->base[CSEQ_PAIR_PAIRBASE(cseq,i)]);
  }
  exit(0);
  */

  if( read_packaln != NULL ) {
    packifp = openfile(read_packaln,"r");
    if( packifp == NULL ) {
      fatal("File %s is unopenable - ignoring dump command",dump_packaln);
    } else {
      pal = read_simple_PackAln(packifp);
    }
  } else {
    pal = PackAln_bestmemory_SyWise20(exonscore,cseq,codon_score,nonc_score,start,stop,Probability2Score(1.0/100.0),Probability2Score(1.0/10000.0),Probability2Score(1.0/10.0),NULL,dpri);
  }

  alb = convert_PackAln_to_AlnBlock_SyWise20(pal);


  if( dump_packaln != NULL ) {
    packifp = openfile(dump_packaln,"w");
    if( packifp == NULL ) {
      warn("File %s is unopenable - ignoring dump command",dump_packaln);
    } else {
      show_simple_PackAln(pal,packifp);
    }
  }

  show_score_sequence(alb,pbs,nonc_score,stdout);
/*
  show_StandardOutputOptions(std_opt,alb,pal,"//",stdout);
*/
  genomic = Genomic_from_Sequence(al->seq[0]);
  gr = new_GenomicRegion(genomic);
  gr2 = new_GenomicRegion(genomic);

  add_Genes_to_GenomicRegion_new(gr,alb);

  
  show_GenomicRegionOptions(sgro,gr,ct,"//",stdout);
  
  return 0;
}




