#include "version.h" 
#include "assemblygraph.h"
#include "assembly_stream_fasta.h"
#include "kmer_glib_index.h"
#include "kmer_hash.h"
#include "assemblystats.h"


char * program_name = "badkmer";


void show_help(FILE * ofp)
{
  fprintf(ofp,"%s (%s)\n",program_name,VERSION_NUMBER);

  fprintf(ofp,"Output options\n");
  fprintf(ofp,"  -depth [2] depth of link\n");
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
  boolean use_hash = 1;

  BaseNode * bn;

  SinglePosSequence * sps;
  AssemblySequence * aseq;

  int i,j;
  int depth = 2;
  char buffer[256];
  
  int kmer = 15;
  int load;

  AssemblyGraph * ag;
  AssemblyGraphStats * astat;
  AssemblySequenceStream * as;

  FILE * ifp;
  KmerIndexInterface * kii;

  strip_out_integer_argument(&argc,argv,"kmer",&kmer);
  strip_out_integer_argument(&argc,argv,"depth",&depth);
  strip_out_boolean_def_argument(&argc,argv,"hash",&use_hash);
  
  strip_out_standard_options(&argc,argv,show_help,show_version);

  if( argc == 1 ) {
    show_help(stdout);
    exit(1);
  }

  ifp = openfile(argv[1],"r");
  
  as = plain_fasta_AssemblySequenceStream(ifp);
  
  if( use_hash == 0 ) {
    kii = new_interface_KmerGlibIndex(kmer);
  } else {
    kii = new_KmerHashIndex(kmer);
  }
  
  
  ag = new_AssemblyGraph(kii,500);

  load = load_AssemblyGraph(ag,as);
  


  if( show_stats_before  ) {
      astat = new_AssemblyGraphStats(ag);
      show_AssemblyGraphStats(astat,10,10,stdout);
  }
  
  for(bn=ag->bg->start;bn!= NULL;bn = bn->next_node ) {
    for(i=0;i<bn->link_len;i++) {
      /* only show links attached to the left */
      if( bn->link[i]->left == bn && bn->link[i]->sequence_label_len > depth ) {
	reverse_map_dna_number(bn->link[i]->right->dna_number,kmer,buffer);
	buffer[kmer] = '\0';
	for(j=0;j<bn->link[i]->sequence_label_len;j++) {
	  sps = lookup_Sequence_SinglePosSpace(ag->sp,bn->link[i]->sequence_label[j]);
	  aseq = (AssemblySequence*) sps->data;
	  printf("%d %s %s %ld\n",bn->link[i]->sequence_label_len,buffer,aseq->seq->name,bn->link[i]->sequence_label[j] - sps->start);
	}
      }
    }
  }

  return 0;
}
