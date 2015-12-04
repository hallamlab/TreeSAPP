#include "aligngenemodel.h"
#include <stdio.h>



int main(int argc,char ** argv)
{
  AlignGeneModelParam * agmp;
  SeqAlign * sal;
  Sequence * seq;

  CompProb * cp;
  CodonTable * ct;
  DnaProbMatrix * dm;

  dm = DnaProbMatrix_from_match(0.8,NMaskType_VARIABLE);

  ct = read_CodonTable_file("codon.table");

  cp = read_Blast_file_CompProb("wag85");

  sal = SeqAlign_alloc_std();

  seq = new_Sequence_from_strings("new1","ATGGGG");
  add_SeqAlign(sal,seq);

  seq = new_Sequence_from_strings("new2","ATGGGT");
  add_SeqAlign(sal,seq);

  agmp = std_AlignGeneModelParam(cp,dm,ct);

  printf("Codon 1 %f vs %f Codon 2 %f v %f \n",
	 coding_probability_AlignGeneModel(sal,2,agmp),
	 non_coding_probability_AlignGeneModel(sal,2,agmp),
	 coding_probability_AlignGeneModel(sal,5,agmp),
	 non_coding_probability_AlignGeneModel(sal,5,agmp)
    );


}
