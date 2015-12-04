#define WISE2_CROSS_HMMER2
#include "wise2xhmmer2.h"


#include "fivestate.h"
#include "version.h"
#include "seqaligndisplay.h"
#include "hitlist.h"



char * program_name = "fivestarscan";

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
  fprintf(ofp,"%s fivestar-directory protein-file-fasta\n",program_name);

  fprintf(ofp,"  -ga     gathering cutoff (bits)\n");
  
  show_help_DPRunImpl(ofp);

  show_help_HitListOutputImpl(ofp);

  show_standard_options(ofp);
  
}

int main(int argc,char **argv)
{
  FiveStateFrameSet * frame;
  
  FiveStateModel * fsm;
  FiveStateScore * fss;

  RandomModel * rm;

  double gathering_cutoff = 0.0;
  double bits;
  int i;

  FILE * seqin;
  Sequence * seq;
  
  DPRunImpl * dpri;

  HitList * hl;
  HitPair * hp;
  HitAln  * ha;
  HitListOutputImpl * hloi;

  PackAln * pal;
  AlnBlock * alb;

  ComplexSequence * cseq;
  ComplexSequenceEvalSet * cses;
  

  strip_out_float_argument(&argc,argv,"ga",&gathering_cutoff);

  hloi = new_HitListOutputImpl_from_argv(&argc,argv);

  strip_out_standard_options(&argc,argv,show_help,show_version);
  if( argc != 3 ) {
    show_help(stdout);
    exit(12);
  }

  rm = default_RandomModel();

  frame = read_FiveStateFrameSet_file(argv[1],"block.str");
  if( frame == NULL ) 
    fatal("Unable to make FiveStateModel from context %s, block.str file",argv[1]);
  fsm   = FiveStateModel_from_FiveStateFrameSet(frame);
  /*    dump_FiveStateModel(fsm,stdout); */
  fsm->name = stringalloc(argv[1]);
  
  fold_RandomModel_into_FiveStateModel(fsm,rm);  

  /* converts probabilities to integers for calculation */
  fss = FiveStateScore_from_FiveStateModel(fsm);
  
  seqin = fopen(argv[2],"r");
  if( seqin == NULL ) {
    fatal("Unable to open %s",argv[2]);
  }

  cses = default_aminoacid_ComplexSequenceEvalSet();

  hl = HitList_alloc_std();

  while( (seq = read_fasta_Sequence(seqin)) ) {
    cseq = new_ComplexSequence(seq,cses);

    pal = PackAln_bestmemory_FiveStateProtein(fss,cseq,NULL,dpri);
    alb = convert_PackAln_to_AlnBlock_FiveStateProtein(pal);

    hp = HitPair_alloc_std();
    add_HitList(hl,hp);

    hp->query = seq;
    hp->bit_score = Score2Bits(pal->score);
    hp->raw_score = pal->score;
    
    
  }

  show_HitList_HitListOutputImpl(hloi,hl,stdout);

  return 0;
}
