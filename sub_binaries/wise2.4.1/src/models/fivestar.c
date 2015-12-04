#define WISE2_CROSS_HMMER2
#include "wise2xhmmer2.h"


#include "fivestate.h"
#include "version.h"
#include "seqaligndisplay.h"
#include "standardout.h"
#include "threestatedp.h"


char * program_name = "fivestar";

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
  fprintf(ofp,"%s hmmer-file protein-file\n",program_name);
  
  fprintf(ofp,"[no]five   do fivestar method\n");

  show_help_StandardOutputOptions(ofp);

  show_help_DPRunImpl(ofp);

  show_standard_options(ofp);
}





int main(int argc,char **argv)
{
  char * hmmer_file;
  char * sequence_file;

  ThreeStateModel * tsm;
  ThreeStateScore * tss;

  FiveStateFrameSet * frame;
  FiveStateModel * fsm;
  FiveStateScore * fss;

  RandomModel * rm;

  Sequence * seq;
  Protein * tseq;

  ComplexSequence * cseq;
  ComplexSequenceEvalSet * cses;

  PackAln * pal;
  AlnBlock * alb;
  
  boolean do_fivestar = FALSE;

  int i;

  DPRunImpl * dpri = NULL;
  StandardOutputOptions * std_opt = NULL;
  rm = default_RandomModel();

  dpri = new_DPRunImpl_from_argv(&argc,argv);
  std_opt = new_StandardOutputOptions_from_argv(&argc,argv);

  if( dpri == NULL ) {
    fatal("Unable to build DPRun implementation. Bad arguments");
  }

  strip_out_boolean_def_argument(&argc,argv,"five",&do_fivestar);


  strip_out_standard_options(&argc,argv,show_help,show_version);
  if( argc != 3 ) {
    show_help(stdout);
    exit(12);
  }


  if( do_fivestar ) {
    frame = read_FiveStateFrameSet_file(argv[1],"block.str");
    if( frame == NULL ) 
      fatal("Unable to make FiveStateModel from context %s, block.str file",argv[1]);
    fsm   = FiveStateModel_from_FiveStateFrameSet(frame);
    /*    dump_FiveStateModel(fsm,stdout); */
    fsm->name = stringalloc(argv[1]);
  } else {
    tsm = HMMer2_read_ThreeStateModel(argv[1]);
    fsm = FiveStateModel_from_flat_ThreeStateModel(tsm);
  }
  fold_RandomModel_into_FiveStateModel(fsm,rm);  

  /* converts probabilities to integers for calculation */
  fss = FiveStateScore_from_FiveStateModel(fsm);

  /*  dump_FiveStateScore(fss,stdout); */
  


  seq = read_fasta_file_Sequence(argv[2]);

  assert(seq);

  cses = default_aminoacid_ComplexSequenceEvalSet();
  cseq = new_ComplexSequence(seq,cses);


  /*
   * Actually do the alignment. If fivestar is switched on,
   * use the FiveStateProtein match
   */

  if( do_fivestar ) {
    /* PackAln_bestmemory is Viterbi algorithm */

    pal = PackAln_bestmemory_FiveStateProtein(fss,cseq,NULL,dpri);
    alb = convert_PackAln_to_AlnBlock_FiveStateProtein(pal);
    tseq = pseudo_Protein_from_FiveStateModel(fsm);
  } else {
    fold_RandomModel_into_ThreeStateModel(tsm,tsm->rm);
    tss = ThreeStateScore_from_ThreeStateModel(tsm);
    pal = PackAln_bestmemory_ThreeStateProtein(tss,cseq,NULL,dpri);
    alb = convert_PackAln_to_AlnBlock_ThreeStateProtein(pal);
    tseq = pseudo_Protein_from_ThreeStateModel(tsm);
  }

  show_StandardOutputOptions(std_opt,alb,pal,"//",stdout);

  printf("Score %d Bits %.3f\n",pal->score,Score2Bits(pal->score));

  write_pretty_str_align(alb,tseq->baseseq->name,tseq->baseseq->seq,seq->name,seq->seq,12,50,stdout);

  return 0;

}
