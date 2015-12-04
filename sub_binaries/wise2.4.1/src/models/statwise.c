#include "statwise10.h"
#include "version.h"
#include "genestats.h"
#include "geneutil.h"


char * program_name = "statwise";

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
  fprintf(ofp,"%s sequence-as-fasta\n",program_name);

  show_help_GeneModelParam(ofp);

  show_help_DPRunImpl(ofp);

  show_standard_options(ofp);
}

CodonTable * ct;

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

  Sequence * seq;

  RandomCodon * rc;
  RandomModelDNA * rmd;
  RandomCodonScore * rcs;


  ComplexSequenceEval * splice5;
  ComplexSequenceEval * splice3;
  ComplexSequenceEvalSet * cses;
  ComplexSequence * cseq;


  SyExonScore * exonscore;

  PackAln * pal;
  AlnBlock * alb;

  Genomic * genomic;
  GenomicRegion * gr;
  Protein * trans;

  dpri = new_DPRunImpl_from_argv(&argc,argv);
  if( dpri == NULL ) {
    fatal("Unable to build DPRun implementation. Bad arguments");
  }

  gmp = new_GeneModelParam_from_argv(&argc,argv);

  ct= read_CodonTable_file("codon.table");

  strip_out_standard_options(&argc,argv,show_help,show_version);
  if( argc != 2 ) {
    show_help(stdout);
    exit(12);
  }

  
  if((gm=GeneModel_from_GeneModelParam(gmp)) == NULL ) {
    fatal("Could not build gene model");
  }


  seq = read_fasta_file_Sequence(argv[1]);
  
  assert(seq);

  cses = new_ComplexSequenceEvalSet_from_GeneModel(gm);

  cseq = new_ComplexSequence(seq,cses);

  rc = flat_RandomCodon(ct);
  rmd = RandomModelDNA_std();

  fold_in_RandomModelDNA_into_RandomCodon(rc,rmd);
  rcs = RandomCodonScore_from_RandomCodon(rc);

  exonscore = SyExonScore_flat_model(200,250,0.1,0.1);
  /*
  for(i=0;i<cseq->length;i++) {
    fprintf(stdout,"%d PairSeq is %d score %d\n",i,CSEQ_PAIR_PAIRBASE(cseq,i),nonc_score->base[CSEQ_PAIR_PAIRBASE(cseq,i)]);
  }
  exit(0);
  */
/*
  show_RandomCodonScore(rcs,stdout);


  for(i=3;i<seq->len;i++) {
    fprintf(stdout,"seq %d is %c with score %d\n",i,aminoacid_from_seq(ct,seq->seq+i-2),rcs->codon[CSEQ_GENOMIC_CODON(cseq,i)]);
  }

  exit(0);
*/

  pal = PackAln_bestmemory_StatWise10(exonscore,cseq,rcs,Probability2Score(1.0/10.0),Probability2Score(1.0/10.0),NULL,dpri);
  alb = convert_PackAln_to_AlnBlock_StatWise10(pal);

  mapped_ascii_AlnBlock(alb,id,1,stdout);

  genomic = Genomic_from_Sequence(seq);
  gr = new_GenomicRegion(genomic);

  add_Genes_to_GenomicRegion_GeneWise(gr,1,seq->len,alb,"bollocks",0,NULL);


  for(i=0;i<gr->len;i++) {
    if( gr->gene[i]->ispseudo == TRUE ) {
      fprintf(stdout,"#Gene %d is a pseudo gene - no translation possible\n",i);
    } else {
      trans = get_Protein_from_Translation(gr->gene[i]->transcript[0]->translation[0],ct);
      write_fasta_Sequence(trans->baseseq,stdout);
    }
  } 


  
  return 0;
}




