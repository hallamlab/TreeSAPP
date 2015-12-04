#include "version.h"
#include "largeblockdp.h"
#include "seqaligndisplay.h"


char * program_name = "lba";

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
  fprintf(ofp,"%s query-seq target-seq\n",program_name);
  fprintf(ofp,"Parameters (all as probabilities)\n");
  fprintf(ofp,"  -matchp [0.75] probability of match in block\n");
  fprintf(ofp,"  -gap    [0.1]  probability of opening a gap in a block\n");
  fprintf(ofp,"  -block_open [0.0005] probability of a block occuring\n");
  fprintf(ofp,"  -block_ext [0.8] probability ratio of block extension vs null\n");
  fprintf(ofp,"  -undual    [1.0] probability ratio of un-matched extension vs null, both seqs\n");
  fprintf(ofp,"  -unsingle   [1.0] probability ratio of un-matched extension vs null, one seq\n");
  fprintf(ofp,"Output format\n");
  fprintf(ofp,"  -[no]alb    show alb format (default no)\n");
  fprintf(ofp,"  -[no]pretty show pretty format (default no)\n");


  show_help_DPRunImpl(ofp);


  show_standard_options(ofp);

}



int main(int argc,char ** argv)
{
  Sequence * one;
  Sequence * two;
  PackAln  * pal;
  AlnBlock * alb;
  DnaProbMatrix * dmp;
  DnaMatrix * mat;

  Probability match    = 0.75;
  Probability gap_open = 0.1;
  Probability gap_ext  = 0.3;
  Probability block_open = 0.0005;
  Probability un_dual  = 1.0;
  Probability un_single = 1.0;
  Probability real_ext  = 0.6;

  boolean show_pretty = TRUE;
  boolean show_alb    = FALSE;

  ComplexSequence *cone, *ctwo;  
  ComplexSequenceEvalSet *cses;  

  DPRunImpl * dpri;


  dpri = new_DPRunImpl_from_argv(&argc,argv);

  strip_out_float_argument(&argc,argv,"matchp",&match);
  strip_out_float_argument(&argc,argv,"gap",&gap_open);
  strip_out_float_argument(&argc,argv,"ext",&gap_ext);
  strip_out_float_argument(&argc,argv,"block_open",&block_open);
  strip_out_float_argument(&argc,argv,"block_ext",&real_ext);
  strip_out_float_argument(&argc,argv,"undual",&un_dual);
  strip_out_float_argument(&argc,argv,"unsingle",&un_single);

  strip_out_boolean_def_argument(&argc,argv,"pretty",&show_pretty);
  strip_out_boolean_def_argument(&argc,argv,"alb",&show_alb);

  strip_out_standard_options(&argc,argv,show_help,show_version);


  if( argc != 3 ) {
    show_help(stdout);
    exit(12);
  }

  one = read_fasta_file_Sequence(argv[1]);
  two = read_fasta_file_Sequence(argv[2]);

  dmp = DnaProbMatrix_from_match(match,NMaskType_BANNED);
  flat_null_DnaProbMatrix(dmp);  
  mat = DnaMatrix_from_DnaProbMatrix(dmp);  

  cses = default_DNA_ComplexSequenceEvalSet();  
  cone = new_ComplexSequence(one,cses);  
  ctwo = new_ComplexSequence(two,cses);  


  pal = PackAln_bestmemory_LargeBlockAligner(cone,ctwo,mat,
					     Probability2Score(real_ext),
					     Probability2Score(block_open),
					     Probability2Score(un_dual),
					     Probability2Score(un_single),
					     Probability2Score(gap_open),
					     Probability2Score(gap_ext),NULL,dpri);

  alb = convert_PackAln_to_AlnBlock_LargeBlockAligner(pal);

  if( show_pretty == TRUE )
    write_pretty_seq_align(alb,one,two,12,60,stdout);

  if( show_alb == TRUE )
    dump_ascii_AlnBlock(alb,stdout);

 
}
