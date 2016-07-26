#include "randomdb.h"


main()
{
  int i;
  Sequence * out;
  RandomModel * rm;
  RandomProteinDB * rdp;
  RandomDNADB * rddp;
  RandomModelDNA * rmd;

  rm = default_RandomModel();
  rmd = RandomModelDNA_std();

  rdp = new_flat_RandomProteinDB(rm,50);


  for(i=0;i<5;i++) {
    out = Sequence_from_RandomProteinDB(rdp);
    write_fasta_Sequence(out,stdout);
    free_Sequence(out);
  }

  rddp = new_flat_RandomDNADB(rmd,50);


  for(i=0;i<5;i++) {
    out = Sequence_from_RandomDNADB(rddp);
    write_fasta_Sequence(out,stdout);
    free_Sequence(out);
  }

}
  
  
