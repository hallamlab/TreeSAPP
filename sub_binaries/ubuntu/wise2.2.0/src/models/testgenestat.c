#include "genestats.h"
#include <stdio.h>
#include <stdlib.h>

#include "sequence.h"


int main(int argc,char ** argv)
{
  Sequence * seq;

  ComplexSequence * cs;
  ComplexSequenceEvalSet * cset;


  GeneStats * st;
  GeneModel * gm;
  GeneModelParam * gp;
  FILE * ifp;

  gp = std_GeneModelParam();


  seq = read_fasta_file_Sequence("../../test_data/human.genomic");
  
  ifp = openfile("gene.stat","r");
  
  st = read_GeneStats(ifp); 
  
  /*  dump_GeneStats(st,stdout); */

  fflush(stdout);

  gm = GeneModel_from_GeneStats(st,gp);

  show_GeneModel(gm,stdout);

  fflush(stdout);

  cset = new_ComplexSequenceEvalSet_from_GeneModel(gm);

  cs = new_ComplexSequence(seq,cset);

  show_ComplexSequence(cs,stdout);

}
