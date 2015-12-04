#include "assemblygraph.h"
#include "assembly_stream_fasta.h"
#include "kmer_glib_index.h"
#include "kmer_hash.h"
#include "assemblystats.h"
#include "graphtangle.h"


void test_basic_tangle(char * filename,int kmer,int use_hash)
{
  char buffer[120];
  AssemblyGraph * ag;
  AssemblySequenceStream * as;
  Assembly * assem;

  FILE * ifp;
  KmerIndexInterface * kii;
  long int load;
  AssemblyPath * path;
  BaseNode * bn;
  BaseLink * link;
  int no;


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

  assert(load > 0);

  bn = ag->bg->start;
  

  /*  show_BaseGraph(ag->bg,stderr);*/


  while( bn != NULL ) {
    if( bn->link_len == 1 && (bn->link[0]->state & BaseLink_TANGLE_TOUCH) != BaseLink_TANGLE_TOUCH ) {
      /* start here */


      fprintf(stderr,"Got a unique link! %d state %d\n",bn,bn->link[0]->state);

      bn->link[0]->state |= BaseLink_TANGLE_TOUCH;
      link = bn->link[0];

      if(link->left == bn ) {
	bn = link->right;
      } else {
	bn = link->left;
      }
      
      /* walk along until we hit a tangle */

      while( bn != NULL ) {
	no =next_link_if_unique_BaseLink(link,bn,&link);
	if( no > 1 ) {
	  reverse_map_dna_number(bn->dna_number,ag->bg->kii->kmer_size,buffer);
	  buffer[ag->bg->kii->kmer_size] = '\0';
	  fprintf(stderr,"Investigating tangle around %d [%s] with %d links\n",bn,buffer,bn->link_len);

	  path = find_sequence_label_tangle_path(ag,link,bn);
	  if( path == NULL ) {
	    /* can't handle this - move on! */
	    bn = bn->next_node;
	    break;
	  }


	  show_AssemblyPath(path,ag->bg->kii->kmer_size,stderr);
	  bn = lift_tangled_AssemblyPath(ag,path);
	  fprintf(stderr,"Finished lifting\n");
	  /* need to find the new node now */
	  free_AssemblyPath(path);

	  fprintf(stderr,"handled tangle, new node %d\n",bn);
	} else if( no == 0 ) {
	  bn = bn->next_node;
	  break;
	} else {
	  if(link->left == bn ) {
	    bn = link->right;
	  } else {
	    bn = link->left;
	  }
	  
	  /* fprintf(stderr,"Progressing! with link %d %ld\n",link,link->sequence_label[0]); */
	}
      }

    } else {
      bn = bn->next_node;
    }
  }


  assem = Assembly_from_AssemblyGraph(ag);

  dump_contigs_as_fasta_Assembly(assem,stdout);


}


int main(int argc,char ** argv)
{
  if( argc < 2  || atoi(argv[2]) < 0 ) {
    fprintf(stderr,"Must have kmer size\n");
    exit(0);
  }

  test_basic_tangle(argv[1],atoi(argv[2]),1);

  
  return 0;
}
  
  
  
