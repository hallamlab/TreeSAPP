#include "genewisehsp.h"



int main(int argc,char **argv)
{
  Sequence * a;
  Sequence * b;
  DPEnvelope * dpenv;
  GeneWiseRunPara * p;
  CompMat * mat;
  CodonTable * ct;
  

  a = read_fasta_file_Sequence(argv[1]);
  b = read_fasta_file_Sequence(argv[2]);

  p = GeneWiseRunPara_alloc();

  mat = read_Blast_file_CompMat("blosum62.bla");
  ct = read_CodonTable_file("codon.table");

  dpenv = DPEnvelope_from_protein_gen(a,b,mat,ct,p);
  
  show_DPEnvelope(dpenv,stdout);

}
