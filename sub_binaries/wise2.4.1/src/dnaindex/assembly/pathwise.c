#include "version.h" 
#include "assemblygraph.h"
#include "assembly_stream_fasta.h"
#include "assembly_sanger_project.h"
#include "kmer_glib_index.h"
#include "kmer_hash.h"
#include "assemblystats.h"
#include "graphmanager.h"

char * program_name = "pathwise";


void show_help(FILE * ofp)
{
  fprintf(ofp,"%s (%s)\n",program_name,VERSION_NUMBER);

  fprintf(ofp,"Graph statistics display\n");
  fprintf(ofp,"  -[no]stat_raw    show stats just after sequence load\n");
  fprintf(ofp,"  -[no]stat_error  show stats just after error manip\n");
  fprintf(ofp,"  -[no]stat_tangle show stats just after de-tangling\n");
  fprintf(ofp,"Graph manipulation\n");
  fprintf(ofp,"  -[no]error_free  reads are error free. (off by default)\n");

  fprintf(ofp,"Read input type\n");
  fprintf(ofp,"  -readstream [fasta/sanger] Read input type (fasta by default)\n");
  fprintf(ofp,"  -readmax [number] maximum reads to use (good for debugging)\n");
  show_help_GraphErrorPara(ofp);

  show_help_DepthPara(ofp);

  show_help_AssemblyOutputPara(ofp);

  fprintf(ofp,"Kmer options\n");
  fprintf(ofp,"  -kmer [15] kmer size (15 is max)\n");
  fprintf(ofp,"  -[no]hash  use optimised hash (true)\n");
  show_standard_options(ofp);

}

void show_version(FILE * ofp)
{
  fprintf(ofp,"%s\nVersion: %s\nReleased: %s\nCompiled: %s\n",program_name,VERSION_NUMBER,RELEASE_DAY,COMPILE_DATE);
  fprintf(ofp,"\nThis program is freely distributed under a Gnu Public License\n");
  fprintf(ofp,"The source code is copyright (c) EBI 2003\n");
  fprintf(ofp,"There is no warranty, implied or otherwise on the performance of this program\n");
  fprintf(ofp,"For more information read the GNULICENSE file in the distribution\n\n");
  fprintf(ofp,"Credits: Ewan Birney <birney@ebi.ac.uk> wrote the core code.\n");
  exit(63);   
}


int main(int argc,char ** argv)
{
  boolean show_stats_before = 0;
  boolean show_stats_error = 0;
  boolean show_stats_tangle = 0;
  boolean show_fasta_contigs = 1;
  boolean use_hash = 1;
  boolean no_error = 0;
  
  int readmax = 0;
  int kmer = 15;
  int load;
  int cycle_count;

  GraphErrorPara * gep;
  DepthPara * dp;

  AssemblyGraph * ag;
  AssemblyGraphStats * astat;
  AssemblySequenceStream * as;
  AssemblyOutputPara * aop;
  Assembly * assem;
  FILE * ifp;
  KmerIndexInterface * kii;
  char * streamtype = NULL;

  gep = new_GraphErrorPara_from_argv(&argc,argv);

  dp = new_DepthPara_from_argv(&argc,argv);

  aop = new_AssemblyOutputPara_from_argv(&argc,argv);

  strip_out_integer_argument(&argc,argv,"kmer",&kmer);
  strip_out_boolean_def_argument(&argc,argv,"hash",&use_hash);

  strip_out_boolean_def_argument(&argc,argv,"error_free",&no_error);
  strip_out_boolean_def_argument(&argc,argv,"stat_raw",&show_stats_before);
  strip_out_boolean_def_argument(&argc,argv,"stat_error",&show_stats_error);
  strip_out_boolean_def_argument(&argc,argv,"stat_tangle",&show_stats_tangle);

  strip_out_integer_argument(&argc,argv,"readmax",&readmax);
  streamtype = strip_out_assigned_argument(&argc,argv,"readstream");

  strip_out_standard_options(&argc,argv,show_help,show_version);

  if( argc == 1 ) {
    show_help(stdout);
    exit(1);
  }


  
  if( streamtype == NULL || strcmp(streamtype,"fasta") == 0 ) {
    ifp = openfile(argv[1],"r");
    as = plain_fasta_AssemblySequenceStream(ifp);
  } else if( strcmp(streamtype,"sanger") == 0 ) {
    as = new_sanger_project_AssemblySequenceStream(argv[1],"1c");
  } else {
    fatal("Bad stream type [%s]",streamtype);
  }

  
  if( use_hash == 0 ) {
    kii = new_interface_KmerGlibIndex(kmer);
  } else {
    kii = new_KmerHashIndex(kmer);
  }
  
  
  ag = new_AssemblyGraph(kii,500);

  load = load_AssemblyGraph(ag,as,readmax);
  


  if( show_stats_before  ) {
      astat = new_AssemblyGraphStats(ag);
      show_AssemblyGraphStats(astat,10,10,stdout);
      free_AssemblyGraphStats(astat);
  }
  
  
  cycle_count = find_single_cycles_AssemblyGraph(ag);

  fprintf(stderr,"Found %d cycles\n",cycle_count);

  cycle_count = find_di_cycles_AssemblyGraph(ag);

  fprintf(stderr,"Found %d dicycles\n",cycle_count);

  /* handle errors */

  if( no_error == 0 ) {
    find_simple_errors_AssemblyGraph(ag,gep);
  }


  if( show_stats_error  ) {
      astat = new_AssemblyGraphStats(ag);
      show_AssemblyGraphStats(astat,10,10,stdout);
      free_AssemblyGraphStats(astat);
  }

  /* handle tangles */

  find_and_resolve_tangles_AssemblyGraph(ag);

  /* now resolve internal tangles */
  /*  find_and_resolve_tangles_AssemblyGraph_Para(ag,dp,1);*/



  if( show_stats_tangle  ) {
      astat = new_AssemblyGraphStats(ag);
      show_AssemblyGraphStats(astat,10,10,stdout);
      free_AssemblyGraphStats(astat);
  }


  /* handle supertangles */

  /*  find_and_jump_supertangles_AssemblyGraph(ag);*/

  if( show_fasta_contigs ) {
    assem = Assembly_from_AssemblyGraph(ag,dp);
    dump_contigs_as_fasta_Assembly(assem,aop,stdout);
  }

  return 0;

}
