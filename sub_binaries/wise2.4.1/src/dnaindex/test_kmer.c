#include "kmer_assembly.h"
#include "kmer_direct.h"
#include "largeseqreader.h"
#include "kmer_assembly_untangler.h"
#include "kmer_assembly_error.h"
#include "kmer_assembly_contig.h"
#include "assembly.h"

int main(int argc,char ** argv)
{
  KmerIndexInterface * kii;
  KmerAssemblyIndex * kai;
  SinglePosSpace * sps;

  AssemblySequence * aseq;
  AssemblySequence * mirror;
  
  Assembly * assembly;
  KmerAssemblyContigPara * p;
  

  p = KmerAssemblyContigPara_alloc();
  p->minimum_len = 0;
  p->minimum_depth = 1;


  sps = new_SinglePosSpace(1,2000);
  kii = new_KmerDirectIndex(12);

  kai = new_KmerAssemblyIndex(kii,sps);
  
  while( (aseq = read_plain_fasta_AssemblySequence(stdin,1000000,stderr)) != NULL ) {
    add_AssemblySequence_KmerAssemblyIndex(kai,aseq,1000000);
    mirror = mirrored_AssemblySequence(aseq);
    add_AssemblySequence_KmerAssemblyIndex(kai,mirror,1000000); 
  }

  resolve_forward_errors_KmerAssembly(kai,1);

  remove_errors_KmerAssemblyIndex(kai,1);

  untangle_KmerAssembly(kai);

  assembly = Assembly_from_KmerAssemblyIndex(kai,p);

  show_extensive_stats_KmerAssemblyIndex(kai,stdout);
  
  dump_contigs_as_fasta_Assembly(assembly,stdout);

}
