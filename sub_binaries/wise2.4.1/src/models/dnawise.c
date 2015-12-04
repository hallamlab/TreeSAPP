#include "dnahmmdp.h"
#include "version.h"
#include "seqaligndisplay.h"

char * program_name = "dnawise";


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
  fprintf(ofp,"%s query_alignment target_sequence\n",program_name);

  show_help_DPRunImpl(ofp);

  show_standard_options(ofp);

}




int main(int argc,char ** argv)
{
  DPRunImpl * dpri = NULL;

  SeqAlign * query;
  Sequence * target;
  
  ComplexSequence * cseq;
  ComplexSequenceEvalSet * cses;

  DnaHmmProb * dhp;
  DnaHmmScore * dhs;
  RandomModelDNA * rmd;

  PackAln * pal;
  AlnBlock * alb;

  dpri      = new_DPRunImpl_from_argv(&argc,argv);


  strip_out_standard_options(&argc,argv,show_help,show_version);
  if( argc != 3 ) {
    show_help(stdout);
    exit(12);
  }

  query = read_selex_SeqAlign_file(argv[1]);
  target = read_fasta_file_Sequence(argv[2]);

  assert(query);
  assert(target);

  cses = default_DNA_ComplexSequenceEvalSet();
  dhp  = new_DnaHmmProb_from_SeqAlign_ungapped(query,1.0);
  rmd  = RandomModelDNA_std();
  cseq = new_ComplexSequence(target,cses);

  make_consensus_DnaHmmProb(dhp);

  fold_RandomModelDNA_DnaHmmProb(dhp,rmd,1.0);

  dhs = DnaHmmScore_from_DnaHmmProb(dhp);
  
  pal = PackAln_bestmemory_DnaHmmMatrix(dhs,cseq,Probability2Score(1.0),NULL,dpri);

  alb = convert_PackAln_to_AlnBlock_DnaHmmMatrix(pal);

  write_pretty_str_align(alb,"align",dhp->consensus,target->name,target->seq,12,50,stdout);


}
