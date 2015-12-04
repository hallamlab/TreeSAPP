#ifdef _cplusplus
extern "C" {
#endif
#include "kmer_assembly.h"



# line 56 "kmer_assembly.dy"
void show_extensive_stats_KmerAssemblyIndex(KmerAssemblyIndex * kai,FILE * ofp)
{
  int i;
  int j;

  int    forward_splits  [50];
  int    backward_splits [50];
  kmer_t depth [50];

  int link_2_lengths  [50];

  kmer_t total_numbers = 0;
  kmer_t total_links = 0;
  
  kmer_t kmer = -1;

  KmerAssemblyNode * node;
  KmerAssemblyLink * run;

  fprintf(stderr,"Intialising arrays...\n");
  for(i=0;i<50;i++) {
    forward_splits[i] = backward_splits[i] = depth[i] = link_2_lengths[i] = 0;
  }

  kmer = (*kai->kii->next_filled_kmer)(kai->kii->handle,kmer);

  for(;kmer != -1;  kmer = (*kai->kii->next_filled_kmer)(kai->kii->handle,kmer)) {
    total_numbers++;
    node = (*kai->kii->retrieve_by_kmer)(kai->kii->handle,kmer);
    for(;node != NULL;node= node->node_chain ) {
      total_links += node->next_len;
      
      for(i=0;i<node->next_len;i++) {
	if( node->next[i]->sequence_label_len > 45 ) {
	  depth[45]++;
	} else {
	  depth[node->next[i]->sequence_label_len]++;
	}
      }
      
      if( node->next_len > 1 ) {
	forward_splits[node->next_len]++;
      }
      if( node->prev_len > 1 ) {
	backward_splits[node->prev_len]++;
      }
      
      if( node->prev_len > 1 ) {
	if( node->next_len > 1 ) {
	  link_2_lengths[0]++;
	} else {
	  if( node->next_len == 1 ) {
	    run = node->next[0];
	    for(j=1;run != NULL && j < 45 ;run = run->next->next[0],j++) {
	      if( run->next->next_len > 1 ) {
		break;
	      }
	      if( run->next->next_len == 0 ) {
		break;
	      }
	    }
	    link_2_lengths[j]++;
	  }
	}
      }
    }

  }

  fprintf(ofp,"Coverage stats\n");
  fprintf(ofp,"Depth  %% of links at this depth\n");
  for(i=0;i< 50;i++){
    fprintf(ofp,"%3d  %.2f%%\n",i,(100.0)*((double)depth[i]/(double)(total_links)));
  }
  fprintf(ofp,"\nSplits Forward    Backward\n");
  for(i=1; i < 50 ;i++) {
    fprintf(ofp,"%2d  %4d    %4d\n",i,forward_splits[i],backward_splits[i]);
  }
  fprintf(ofp,"\nDistribution of stream merge length\n");
  fprintf(ofp,"Length   Number\n");
  for(i=0;i<50;i++) {
    fprintf(ofp,"%3d  %4d\n",i,link_2_lengths[i]);
  }
  

}



# line 145 "kmer_assembly.dy"
void add_AssemblySequence_KmerAssemblyIndex(KmerAssemblyIndex * kai,AssemblySequence * aseq,int report)
{
  kmer_t prev_number;
  kmer_t next_number;
  long int i;
  int j;
  char c;
  long start;
  char * seq_str;

  assert(kai != NULL);
  assert(aseq != NULL);
  assert(aseq->seq != NULL);


  start = add_Sequence_SinglePosSpace(kai->sps,aseq->seq->len,(void*)aseq);

  for(i=0;i<aseq->seq->len-kai->kii->kmer_size;i++) {
    if( report != 0 && i % report == 0 && i != 0) {
      fprintf(stderr,"Loaded %ld positions in %s\n",i,aseq->seq->name);
    }

    for(j=0;j<kai->kii->kmer_size+1;j++) {
      c = toupper(aseq->seq->seq[i+j]);
      if( c != 'A' && c != 'T' && c != 'G' && c != 'C' ) {
	break;
      }
    }
    if( j < kai->kii->kmer_size+1 ) {
      continue;
    }
    seq_str = aseq->seq->seq;

    next_number = forward_dna_number_from_string(seq_str+i+1, kai->kii->kmer_size);
    prev_number = forward_dna_number_from_string(seq_str+i, kai->kii->kmer_size);

    /*    fprintf(stderr,"Adding position %d %.*s...%ld,%ld\n",i,kai->kii->kmer_size+1,aseq->seq->seq+i,next_number,prev_number); */

    store_KmerAssemblyLink_KmerAssemblyIndex(kai,prev_number,next_number,aseq->seq->seq[i+1],start+i+1);
  }

  return;
}



# line 191 "kmer_assembly.dy"
boolean store_KmerAssemblyLink_KmerAssemblyIndex(KmerAssemblyIndex * kai,kmer_t prev_number,kmer_t next_number,char base,long label)
{
  int i;
  KmerAssemblyLink * link;
  KmerAssemblyNode * prev;
  KmerAssemblyNode * next;

  assert(kai != NULL);

  prev = (*kai->kii->retrieve_by_kmer)(kai->kii->handle,prev_number);
  
  if( prev != NULL ) {
    for(i=0;i<prev->next_len;i++) {
      if( prev->next[i]->next->number == next_number ) {
	add_sequence_label_KmerAssemblyLink(prev->next[i],label);
	return TRUE;
      }
    }
  }

  /* nope - need new link */

  next = (*kai->kii->retrieve_by_kmer)(kai->kii->handle,next_number);


  if( prev == NULL ) {
    prev = new_KmerAssemblyNode(prev_number);
    assert(prev->number == prev_number);
    
    (*kai->kii->insert_by_kmer)(kai->kii->handle,prev_number,prev);
  } else {
    assert(prev->number >= 0);
  }

  if( next == NULL ) {
    next = new_KmerAssemblyNode(next_number);
    assert(next->number == next_number);

    (*kai->kii->insert_by_kmer)(kai->kii->handle,next_number,next);
  } else {
    assert(next->number >= 0);
  }


  assert(prev != NULL);
  assert(prev->number >= 0);
  assert(prev->prev != NULL);
  assert(prev->next != NULL);

  assert(next != NULL);
  assert(next->number >= 0);
  assert(next->prev != NULL);
  assert(next->next != NULL);



  link = new_KmerAssemblyLink(base);
  add_sequence_label_KmerAssemblyLink(link,label);

  link->prev = prev;
  link->next = next;

  /*  fprintf(stderr,"Adding position to %ld : %ld with %d,%d\n",prev->number,next->number,prev->prev_len,next->next_len);*/

  add_next_KmerAssemblyNode(prev,link);
  add_prev_KmerAssemblyNode(next,link);

  return TRUE;

}

# line 262 "kmer_assembly.dy"
void show_KmerAssemblyIndex(KmerAssemblyIndex * kai,FILE * ofp)
{
  kmer_t kmer;
  KmerAssemblyNode * node;
  
  assert(kai != NULL);
  assert(kai->kii != NULL);
  assert(kai->kii->next_filled_kmer != NULL);

  kmer = -1;
  kmer = (*kai->kii->next_filled_kmer)(kai->kii->handle,kmer);

  for(;kmer != -1;  kmer = (*kai->kii->next_filled_kmer)(kai->kii->handle,kmer)) {
    node = (*kai->kii->retrieve_by_kmer)(kai->kii->handle,kmer);
    show_KmerAssemblyNode(node,kai->kii->kmer_size,0,ofp);
  }
}

# line 280 "kmer_assembly.dy"
void show_KmerAssemblyNode(KmerAssemblyNode * node,int kmer_size,int level,FILE * ofp)
{
  int i;
  int k;
  char buffer[512];

  assert(node != NULL);
  assert(ofp != NULL);

  reverse_map_dna_number(node->number,kmer_size,buffer);
  buffer[kmer_size] = '\0';

  for(k=0;k<level;k++) {
    fprintf(ofp,"  ");
  }

  fprintf(ofp,"Node %ld of sequence %s \n",node->number,buffer);

  for(i=0;i<node->prev_len;i++) {
    for(k=0;k<level;k++) {
      fprintf(ofp,"  ");
    }
    fprintf(ofp,"  ... prev ... %c, %d to %ld\n",node->prev[i]->base,node->prev[i]->sequence_label_len,node->prev[i]->prev->number);
  }

  for(i=0;i<node->next_len;i++) {
    for(k=0;k<level;k++) {
      fprintf(ofp,"  ");
    }
    fprintf(ofp,"  ... next ... %c, %d to %ld\n",node->next[i]->base,node->next[i]->sequence_label_len,node->next[i]->next->number);
  }

  if( node->node_chain != NULL) {
    fprintf(ofp,"Has chained node...\n");
    show_KmerAssemblyNode(node->node_chain,kmer_size,level+1,ofp);
  }

}


# line 320 "kmer_assembly.dy"
KmerAssemblyIndex * new_KmerAssemblyIndex(KmerIndexInterface * kii,SinglePosSpace * sps)
{
  KmerAssemblyIndex * out;

  assert(kii != NULL);
  assert(sps != NULL);

  out = (KmerAssemblyIndex*) malloc(sizeof(KmerAssemblyIndex));

  out->kii = kii;
  out->sps = sps;
  
  return out;
}

# line 335 "kmer_assembly.dy"
KmerAssemblyLink * new_KmerAssemblyLink(char base)
{
  KmerAssemblyLink * out;

  out = (KmerAssemblyLink*) malloc (sizeof(KmerAssemblyLink));

  out->sequence_label = (long*) calloc(KmerAssemblyLink_LABEL_START,sizeof(long));
  out->sequence_label_maxlen = KmerAssemblyLink_LABEL_START;
  out->sequence_label_len = 0;
  
  out->base  = base;
  out->prev  = NULL;
  out->next  = NULL;
  out->state = 0;

  return out;
}

# line 353 "kmer_assembly.dy"
void remove_sequence_label_KmerAssemblyLink(KmerAssemblyLink * kal,long label)
{
  int i;
  assert(kal != NULL);

  for(i=0;i<kal->sequence_label_len;i++) {
    if( kal->sequence_label[i] == label ) {
      kal->sequence_label[i] = kal->sequence_label[--kal->sequence_label_len];
      return;
    }
  }

  fprintf(stderr,"    ...unable to remove label %ld from link %ld (%d labels)\n",label,kal,kal->sequence_label_len);
  for(i=0;i<kal->sequence_label_len;i++) {
    fprintf(stderr,"     [%ld] is %d label\n",kal->sequence_label[i]);
  }

  return;
}

# line 373 "kmer_assembly.dy"
void add_sequence_label_KmerAssemblyLink(KmerAssemblyLink * kal,long label)
{
  assert(kal != NULL);

  if( kal->sequence_label_len >= kal->sequence_label_maxlen ) {
    if( kal->sequence_label_maxlen > KmerAssemblyLink_LABEL_LINEAR ) {
      kal->sequence_label_maxlen += KmerAssemblyLink_LABEL_LINEAR;
    } else {
      kal->sequence_label_maxlen *= 2;
    }
    kal->sequence_label = (long *) realloc (kal->sequence_label,sizeof(long)*kal->sequence_label_maxlen);
  }

  kal->sequence_label[kal->sequence_label_len++] = label;


}


# line 392 "kmer_assembly.dy"
void detach_KmerAssemblyLink(KmerAssemblyIndex * kai,KmerAssemblyLink * link)
{
  assert(kai != NULL);
  assert(link != NULL);

  if( link->prev != NULL ) {
    remove_next_KmerAssemblyNode(link->prev,link);
  }

  if( link->next != NULL ) {
    remove_prev_KmerAssemblyNode(link->next,link);
  }


}

# line 408 "kmer_assembly.dy"
void remove_next_KmerAssemblyNode(KmerAssemblyNode * node,KmerAssemblyLink * next)
{
  int i;
  assert(node != NULL);
  assert(next != NULL);

  for(i=0;i<node->next_len;i++) {
    if( node->next[i] == next ) {
      node->next[i] = node->next[--node->next_len];
      return;
    }
  }
  warn("In node %d, unable to remove %d as a next node, graph collapsing....",node,next);

  return;
}

# line 425 "kmer_assembly.dy"
void add_next_KmerAssemblyNode(KmerAssemblyNode * kan,KmerAssemblyLink * kal)
{
  assert(kan != NULL);
  assert(kal != NULL);
  assert(kan->next != NULL);

  if( kan->next_len >= kan->next_maxlen ) {
    if( kan->next_maxlen > KmerAssemblyNode_LINK_LINEAR ){
      kan->next = (KmerAssemblyLink **) realloc (kan->next,sizeof(KmerAssemblyLink*) * (kan->next_maxlen + KmerAssemblyNode_LINK_LINEAR));
      kan->next_maxlen = kan->next_maxlen + KmerAssemblyNode_LINK_LINEAR;
    } else {
      kan->next = (KmerAssemblyLink **) realloc (kan->next,sizeof(KmerAssemblyLink*) * (kan->next_maxlen * 2));
      kan->next_maxlen = kan->next_maxlen * 2;
    }
  }

  kan->next[kan->next_len++] = kal;

}

# line 445 "kmer_assembly.dy"
void remove_prev_KmerAssemblyNode(KmerAssemblyNode * node,KmerAssemblyLink * prev)
{
  int i;
  assert(node != NULL);
  assert(prev != NULL);


  for(i=0;i<node->prev_len;i++) {
    if( node->prev[i] == prev ) {
      node->prev[i] = node->prev[--node->prev_len];
      return;
    }
  }

  warn("In node %d, unable to remove %d as a prev node",node,prev);
  return;
}


# line 464 "kmer_assembly.dy"
void add_prev_KmerAssemblyNode(KmerAssemblyNode * kan,KmerAssemblyLink * kal)
{
  assert(kan != NULL);
  assert(kal != NULL);
  assert(kan->prev != NULL);

  if( kan->prev_len >= kan->prev_maxlen ) {
    if( kan->prev_maxlen > KmerAssemblyNode_LINK_LINEAR ){
      kan->prev = (KmerAssemblyLink **) realloc (kan->prev,sizeof(KmerAssemblyLink*) * (kan->prev_maxlen + KmerAssemblyNode_LINK_LINEAR));
      kan->prev_maxlen = kan->prev_maxlen + KmerAssemblyNode_LINK_LINEAR;
    } else {
      kan->prev = (KmerAssemblyLink **) realloc (kan->prev,sizeof(KmerAssemblyLink*) * (kan->prev_maxlen * 2));
      kan->prev_maxlen = kan->prev_maxlen * 2;
    }
  }

  kan->prev[kan->prev_len++] = kal;

}

# line 484 "kmer_assembly.dy"
KmerAssemblyNode * new_KmerAssemblyNode(kmer_t number)
{
  KmerAssemblyNode * out;

  out = (KmerAssemblyNode*) malloc( sizeof(KmerAssemblyNode));

  assert(out != NULL);

  out->number = number;
  out->prev = (KmerAssemblyLink **) calloc (KmerAssemblyNode_LINK_START,sizeof(KmerAssemblyLink));
  out->next = (KmerAssemblyLink **) calloc (KmerAssemblyNode_LINK_START,sizeof(KmerAssemblyLink));

  out->prev_maxlen = out->next_maxlen = KmerAssemblyNode_LINK_START;
  out->prev_len = out->next_len = 0;

  out->node_chain = NULL;

  return out;
}


# line 471 "kmer_assembly.c"

#ifdef _cplusplus
}
#endif
