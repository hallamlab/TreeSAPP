#include "aligngenemodel.h"




int main(int argc,char **argv)
{
  int i;
  AlignGeneModelParam * agmp;
  GeneStats * gs;
  GeneModelParam * gmp = NULL;
  CompProb * comp_prob;
  DnaProbMatrix * dm;
  CodonTable * ct;
  RandomModel * rm;

  Sequence * test;

  ct = read_CodonTable_file("codon.table");

  rm = default_RandomModel();

  comp_prob = read_Blast_file_CompProb("wag85");

  gmp = new_GeneModelParam_from_argv(&argc,argv);

  dm = DnaProbMatrix_from_match(0.8,NMaskType_VARIABLE);
  
  if((gs=GeneStats_from_GeneModelParam(gmp)) == NULL ) {
    fatal("Could not build gene stats");
  }
  
  agmp = std_AlignGeneModelParam(comp_prob,dm,ct,gs);

  test = read_fasta_file_Sequence(argv[1]);
  
  assert(test);
  
  for(i=0;i<test->len;i++) {
    fprintf(stdout,"%c ss5 %.6f ss3 %.6f\n",test->seq[i],prob_SpliceSiteProb(agmp->ss5,test,i),prob_SpliceSiteProb(agmp->ss3,test,i));
  }


}

