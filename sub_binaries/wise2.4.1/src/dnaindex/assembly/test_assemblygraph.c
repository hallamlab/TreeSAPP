#include "assemblygraph.h"
#include "assembly_stream_fasta.h"
#include "kmer_glib_index.h"
#include "kmer_hash.h"
#include "assemblystats.h"


void test_basic_graph(char * filename,int kmer,int use_hash)
{
  AssemblyGraph * ag;
  AssemblyGraphStats * astat;
  AssemblySequenceStream * as;
  Assembly * assem;
  FILE * ifp;
  KmerIndexInterface * kii;
  long int load;


  ifp = openfile(filename,"r");
  assert(ifp);

  
  as = plain_fasta_AssemblySequenceStream(ifp);

  if( use_hash == 0 ) {
    kii = new_interface_KmerGlibIndex(kmer);
  } else {
    kii = new_KmerHashIndex(kmer);
  }


  ag = new_AssemblyGraph(kii,500);

  load = load_AssemblyGraph(ag,as);

  astat = new_AssemblyGraphStats(ag);

  show_AssemblyGraphStats(astat,10,10,stdout);

  assert(load > 0);
  /*
  assem = Assembly_from_AssemblyGraph(ag);

  dump_contigs_as_fasta_Assembly(assem,stdout);
  */
}


int main(int argc,char ** argv)
{
  test_basic_graph(argv[1],atoi(argv[2]),0);

  return 0;
}
  
  
  
