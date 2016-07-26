#include "dnaalign.h"
#include "seqaligndisplay.h"


main(int argc,char ** argv)
{
  Sequence * one;
  Sequence * two;
  DnaMatrix * mat;
  DnaStartEnd * dse;
  AlnBlock * alb;

  one = read_fasta_file_Sequence(argv[1]);
  two = read_fasta_file_Sequence(argv[2]);

  mat = identity_DnaMatrix(4,-1);

  dse = DnaStartEnd_from_policy("local");

  alb = make_align_dnaalign(one,two,mat,dse,-5,-1,-5,-1);

  dump_ascii_AlnBlock(alb,stdout);

  write_pretty_seq_align(alb,one,two,12,60,stdout);
}

