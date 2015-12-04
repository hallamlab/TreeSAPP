#include "eulerindex.h"


int main(int argc,char ** argv)
{
  EulerGraph * eg;
  Sequence * dna;
  int i;
  int j;

  eg = new_EulerGraph(13);

  while( (dna = read_fasta_Sequence(stdin)) != NULL ) {
    store_Sequence_EulerGraph(eg,dna);
  }


  for(i=0;i<eg->node_len;i++) {
    if( eg->node[i] != NULL && eg->node[i]->link_len > 1) {
      for(j=0;j<eg->node[i]->link_len;j++) {
	if( eg->node[i]->link[j]->depth == 1 ) {
	  fprintf(stderr,"Link has depth %d\n",eg->node[i]->link[j]->depth);
	  if( can_resolve_error_EulerGraph(eg,eg->node[i]->link[j]) ) {
	    fprintf(stderr,"....resolving...\n");
	    if( resolve_error_EulerGraph(eg,eg->node[i]->link[j]) == FALSE ) {
	      remove_EulerLink_forward_EulerNode(eg->node[i],eg->node[i]->link[j]);
	    }
	  } else {
	    remove_EulerLink_forward_EulerNode(eg->node[i],eg->node[i]->link[j]);
	  }
	}
      }
    }
  }

  /* removing remaining single links etc */

  for(i=0;i<eg->node_len;i++) {
    if( eg->node[i] != NULL && eg->node[i]->link_len > 1) {
      for(j=0;j<eg->node[i]->link_len;j++) {
	if( eg->node[i]->link[j]->depth < 10 ) {
	    remove_EulerLink_forward_EulerNode(eg->node[i],eg->node[i]->link[j]);
	    j--;
	}
      }
    }

    if( eg->node[i] != NULL && eg->node[i]->back_len > 1) {
      
      for(j=0;j<eg->node[i]->back_len;j++) {
	if( eg->node[i]->back[j]->depth < 10 ) {
	  fprintf(stderr,"...REMOVING... %d (%d)\n",i,eg->node[i]->back_len);
	  remove_EulerLink_backward_EulerNode(eg->node[i],eg->node[i]->back[j]);
	  j--;
	  fprintf(stderr,"   ...now %d (%d)\n",i,eg->node[i]->back_len);
	}
      }
    }
  }
  

  /* here goes the tangles... */

  for(i=0;i<eg->node_len;i++) {
    if( eg->node[i] != NULL && (eg->node[i]->link_len > 1 || eg->node[i]->back_len > 1)) {
      fprintf(stderr,"BEFORE UNTANGLE:%d has leaving %d and joining %d\n",eg->node[i]->number,eg->node[i]->link_len,eg->node[i]->back_len);
      for(j=0;j<eg->node[i]->link_len;j++) {
	fprintf(stderr,"...link has %d depth\n",eg->node[i]->link[j]->depth);
      }
      for(j=0;j<eg->node[i]->back_len;j++) {
	fprintf(stderr,"...back has %d depth\n",eg->node[i]->back[j]->depth);
      }

    }
  }

  for(i=0;i<eg->node_len;i++) {
    if( eg->node[i] != NULL && eg->node[i]->link_len > 1 ) {
      fprintf(stderr,"MAIN: Attempting to untangle %d\n",i);
      untangle_from_split_EulerNode(eg,eg->node[i],100);
    }
  }

  for(i=0;i<eg->node_len;i++) {
    if( eg->node[i] != NULL && (eg->node[i]->link_len > 1 || eg->node[i]->back_len > 1)) {
      fprintf(stderr,"AFTER UNTANGLE:%d has leaving %d and joining %d\n",eg->node[i]->number,eg->node[i]->link_len,eg->node[i]->back_len);
    }
  }

  /* this artificially breaks at tangles */
  /*
  for(i=0;i<eg->node_len;i++) {
    if( eg->node[i] != NULL && eg->node[i]->link_len > 1 ) {
      warn("Breaking at tangle at %d\n",i);
      for(j=0;j<eg->node[i]->link_len;j++) {
	eg->node[i]->link[j]->next->is_leftmost = 1;
      }
    }
  }
  */

  for(i=0;i<eg->node_len;i++) {
    if( eg->node[i] != NULL && eg->node[i]->is_leftmost == 1 ) {
      dna = read_Sequence_EulerNode(eg->node[i]);
      write_fasta_Sequence(dna,stdout);
      free_Sequence(dna);
    }
  }

}
