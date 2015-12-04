#include "kmer_assembly.h"
#include "kmer_direct.h"
#include "largeseqreader.h"




int main(int argc,char ** argv)
{
  KmerIndexInterface * kii;
  KmerAssemblyIndex * kai;
  SinglePosSpace * sps;
  Sequence * seq;


  sps = new_SinglePosSpace(1,200000000);
  kii = new_KmerDirectIndex(16);

  kai = new_KmerAssemblyIndex(kii,sps);
  
  while( (seq = read_large_dna_Sequence(stdin,1000000)) != NULL ) {
    add_Sequence_KmerAssemblyIndex(kai,seq,1000000);
  }

  show_extensive_stats_KmerAssemblyIndex(kai,stdout);

  return 0;


}
