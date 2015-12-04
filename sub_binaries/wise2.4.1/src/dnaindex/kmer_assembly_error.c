#ifdef _cplusplus
extern "C" {
#endif
#include "kmer_assembly_error.h"


# line 14 "kmer_assembly_error.dy"
boolean lift_indel_error_KmerAssembly(KmerAssemblyIndex * kai,KmerAssemblyPath * real,KmerAssemblyPath * error,long int * start_labels,int label_len,int error_position,int error_len)
{
  int i,j,k,true;

  for(i=0,true=0;i<error->stack_len && true<real->stack_len;i++,true++) {
    if( i == error_position ) {
      if( error_len > 0 ) {
	for(j=0;j<error_len;) {
	  i++;
	  j++;
	}
      } else {
	error_len = abs(error_len);
	for(j=0;j<error_len;) {
	  true++;
	  j++;
	}
      }
    }

    for(k=0;k<label_len;k++) {
      /* remove label from error, add to real */
      if( error->stack[i] == NULL ) {
	fprintf(stderr,"issue - error path is missing a link!");
      } else {
	remove_sequence_label_KmerAssemblyLink(error->stack[i],start_labels[k]+i);
      }

      if( real->stack[true] == NULL ) {
	fprintf(stderr,"issue - true path is missing a link!");
      } else {
	add_sequence_label_KmerAssemblyLink(real->stack[true],start_labels[k]+i);
      }
      if( error->stack[i] != NULL && error->stack[i]->sequence_label_len == 0 ) {
	detach_KmerAssemblyLink(kai,error->stack[i]);
      }
    }
  }

  return TRUE;

}

# line 57 "kmer_assembly_error.dy"
boolean lift_forward_error_KmerAssembly(KmerAssemblyIndex * kai,KmerAssemblyPath * real,KmerAssemblyPath * error,long int * start_labels,int label_len)
{
  int i,k;

  for(i=0;i<error->stack_len;i++) {
    for(k=0;k<label_len;k++) {
      /* remove label from error, add to real */
      remove_sequence_label_KmerAssemblyLink(error->stack[i],start_labels[k]+i);
      add_sequence_label_KmerAssemblyLink(real->stack[i],start_labels[k]+i);
      if( error->stack[i]->sequence_label_len == 0 ) {
	detach_KmerAssemblyLink(kai,error->stack[i]);
      }
    }
  }

  return TRUE;
}

# line 75 "kmer_assembly_error.dy"
int mark_tangles_KmerAssembly(KmerAssemblyIndex * kai)
{
  int i;
  kmer_t kmer;
  KmerAssemblyNode * node;
  char buffer[256];

  int count = 0;

  kmer = -1;
  kmer = (*kai->kii->next_filled_kmer)(kai->kii->handle,kmer);

  for(;kmer != -1;  kmer = (*kai->kii->next_filled_kmer)(kai->kii->handle,kmer)) {
    node = (*kai->kii->retrieve_by_kmer)(kai->kii->handle,kmer);
    if( node->prev_len > 1 ) {
      reverse_map_dna_number(node->number,kai->kii->kmer_size,buffer);
      buffer[kai->kii->kmer_size] = '\0';

      fprintf(stderr,"Marking node (%ld) [%s] as next tangled\n",node->number,buffer);
      count++;
      for(i=0;i<node->prev_len;i++) {
	node->prev[i]->state |= KMER_ASSEMBLY_NEXT_TANGLED;
      }
    }
   
    if( node->next_len > 1 ) {
      reverse_map_dna_number(node->number,kai->kii->kmer_size,buffer);
      buffer[kai->kii->kmer_size] = '\0';


      fprintf(stderr,"Marking node (%ld) [%s] as prev tangled\n",node->number,buffer);
      count++;
      for(i=0;i<node->next_len;i++) {
	node->next[i]->state |= KMER_ASSEMBLY_PREV_TANGLED;
      }
    }
  }

  return count;
}

# line 116 "kmer_assembly_error.dy"
int resolve_forward_errors_KmerAssembly(KmerAssemblyIndex * kai,int depth,int verbose,int max_path_enum)
{
  int i;
  int total = 0;

  kmer_t kmer;
  KmerAssemblyNode * node;

  char buffer[128];
  SinglePosSequence * sps;
  
  assert(kai != NULL);
  assert(kai->kii != NULL );
  assert(kai->kii->next_filled_kmer != NULL);

  kmer = -1;
  kmer = (*kai->kii->next_filled_kmer)(kai->kii->handle,kmer);

  for(;kmer != -1;  kmer = (*kai->kii->next_filled_kmer)(kai->kii->handle,kmer)) {
    node = (*kai->kii->retrieve_by_kmer)(kai->kii->handle,kmer);
    assert(node != NULL);

    if( node->next_len < 2 ) {
      continue;
    }

    for(i=0;i<node->next_len;i++) {
      if( node->next[i]->sequence_label_len <= depth ) {
	if( verbose ) {
	  reverse_map_dna_number(node->next[i]->next->number,kai->kii->kmer_size,buffer);
	  buffer[kai->kii->kmer_size] = '\0';
	  sps = lookup_Sequence_SinglePosSpace(kai->sps,node->next[i]->sequence_label[0]);

	  fprintf(stderr,"Attempting to resolve error [%s] from sequence %s, depth %d\n",buffer,((AssemblySequence*)sps->data)->seq->name,node->next[i]->sequence_label_len);
	}

	total += attempt_resolve_forward_error_KmerAssembly(kai,node,depth,max_path_enum);
	break;
      }
    }
  }
  
  return total;
}



# line 163 "kmer_assembly_error.dy"
void remove_errors_KmerAssemblyIndex(KmerAssemblyIndex * kai,int max_error_depth)
{
  int i;
  kmer_t kmer;
  KmerAssemblyNode * node;

  assert(kai != NULL);
  assert(kai->kii != NULL);
  assert(kai->kii->next_filled_kmer != NULL);

  kmer = -1;
  kmer = (*kai->kii->next_filled_kmer)(kai->kii->handle,kmer);

  for(;kmer != -1;  kmer = (*kai->kii->next_filled_kmer)(kai->kii->handle,kmer)) {
    node = (*kai->kii->retrieve_by_kmer)(kai->kii->handle,kmer);
    assert(node);

    if( node->next_len > 1 ) {
      for(i=0;i<node->next_len;i++) {
	if( node->next[i]->sequence_label_len <= max_error_depth ) {
	  detach_KmerAssemblyLink(kai,node->next[i]);
	}
      }
    }

    if( node->prev_len > 1 ) {
      for(i=0;i<node->prev_len;i++) {
	if( node->prev[i]->sequence_label_len <= max_error_depth ) {
	  detach_KmerAssemblyLink(kai,node->prev[i]);
	}
      }
    }
  }
    
  
}

# line 200 "kmer_assembly_error.dy"
int attempt_resolve_forward_error_KmerAssembly(KmerAssemblyIndex * kai,KmerAssemblyNode * node,int depth,int max_path_enum)
{
  int success;
  int i,j,k;
  KmerAssemblyPath * final_real;
  KmerAssemblyPath * final_error;
  char * alphabet = "ATGC";
  char buffer[256];
  char next[256];
  char selfbase;

  kmer_t new_number;
  
  int no_success = 0;
  int no_failures = 0;

  assert(kai != NULL);
  assert(node != NULL);


  final_real = new_KmerAssemblyPath();
  final_error = new_KmerAssemblyPath();

  if( node->next_len == 0 ) {
    warn("No outgoing links, therefore no errors");
    return TRUE;
  }

  if( node->next_len == 0 ) {
    warn("Only one outgoing link, so no error");
    return TRUE;
  }

  /* look for links with depth outgoing */
  for(i=0;i<node->next_len;i++) {
    if( node->next[i]->sequence_label_len <= depth ) {
      /* see whether this is an error*/

      success = 0;
      reverse_map_dna_number(node->next[i]->next->number,kai->kii->kmer_size,buffer);
      
      selfbase = buffer[kai->kii->kmer_size-1];

      /* test each base */
      for(j=0;j<4;j++) {
	if( alphabet[j] == selfbase ) {
	  continue;
	}
	buffer[kai->kii->kmer_size-1] = alphabet[j];
	new_number = forward_dna_number_from_string(buffer,kai->kii->kmer_size);
	for(k=0;k<node->next_len;k++) {
	  if( k == i ) {
	    continue;
	  }
	  if( new_number == node->next[k]->next->number ) {
	    /* potential fix */
	    final_error->stack_len = 0;
	    final_real->stack_len  = 0;

	    fprintf(stderr,"Found potential fix of %c to %c in %.*s (real:%c (%d), error:%c (%d))\n",selfbase,alphabet[j],kai->kii->kmer_size,buffer,node->next[k]->base,k,node->next[i]->base,i);

	    /* consider potential insertion of one base */

	    if( find_indel_path_KmerAssembly(node->next[k],node->next[i],-1,kai->kii->kmer_size-1,final_real,final_error,kai->kii->kmer_size+1,max_path_enum) == TRUE ) {
	      fprintf(stderr,"SUCCESS in insertion modelling\n");
	      lift_indel_error_KmerAssembly(kai,final_real,final_error,node->next[i]->sequence_label,node->next[i]->sequence_label_len,kai->kii->kmer_size-1,-1);
	      success = 1;
	      break;
	    }

	    if( find_indel_path_KmerAssembly(node->next[k],node->next[i],1,kai->kii->kmer_size-1,final_real,final_error,kai->kii->kmer_size+1,max_path_enum) == TRUE ) {
	      fprintf(stderr,"SUCCESS in deletion modelling\n");
	      lift_indel_error_KmerAssembly(kai,final_real,final_error,node->next[i]->sequence_label,node->next[i]->sequence_label_len,kai->kii->kmer_size-1,1);
	      success = 1;
	      break;
	    }
	      


	    /* otherwise substitution is all we can really suggest... */

	    if( find_error_path_KmerAssembly(node->next[k],node->next[i],alphabet[j],kai->kii->kmer_size-1,final_real,final_error,kai->kii->kmer_size+1,max_path_enum) == TRUE ) {
	      /* fix error */
	      /*
		fprintf(stderr,"Found Error: Real\n");
		show_KmerAssemblyPath(final_real,stderr);
		fprintf(stderr,"   Error\n");
		show_KmerAssemblyPath(final_error,stderr);
		fprintf(stderr,"\n");
	      */		  

	      fprintf(stderr,"SUCCESS in substitution modelling\n");

	      lift_forward_error_KmerAssembly(kai,final_real,final_error,node->next[i]->sequence_label,node->next[i]->sequence_label_len);
	      success = 1;
	      break;
	    }
	  }
	}
      }
      if( success == 0 ) {
	no_failures++;
      } else {
	no_success++;
      }
    }
  }
  
  free_KmerAssemblyPath(final_real);
  free_KmerAssemblyPath(final_error);

  return no_success;
}


# line 315 "kmer_assembly_error.dy"
boolean find_indel_path_KmerAssembly(KmerAssemblyLink *real,KmerAssemblyLink * error,int delete_length,int proposed_position,KmerAssemblyPath * final_real,KmerAssemblyPath * final_error,int max_search,int max_path_enum)
{
  int current_path = 1;

  return extend_indel_path_KmerAssembly(real,error,0,0,delete_length,proposed_position,final_real,final_error,max_search,&current_path,max_path_enum);
}

# line 322 "kmer_assembly_error.dy"
boolean extend_indel_path_KmerAssembly(KmerAssemblyLink* real,KmerAssemblyLink * error,int real_pos,int error_pos,int delete_length,int proposed_position,KmerAssemblyPath * final_real,KmerAssemblyPath * final_error,int max_search,int *current_path,int max_path_enum)
{
  int real_start;
  int error_start;
  int k;
  int i,j;

  KmerAssemblyLink * temp_real[256];
  KmerAssemblyLink * temp_error[256];

  if( max_search > 256 ) {
    fatal("Madness. Max search greater than 256 - shouldn't happen!");
  }

  if( max_path_enum < *current_path ) {
    fprintf(stderr,"Terminating due too many paths\n");
    return FALSE;
  }

  real_start = real_pos;
  error_start = error_pos;


  /*  fprintf(stderr,"Entering into extension code with real %d, error %d\n",real_pos,error_pos);*/

  while( real != NULL && error != NULL ) {

    if( proposed_position != real_pos ) {
      if( real->base != error->base ) {
	fprintf(stderr,"in considering indel (%d, path %d), real (%c) and error (%c) do not agree at position %d,%d\n",delete_length,current_path,real->base,error->base,real_pos,error_pos);
	return FALSE;
      } 
    } else {
      if( delete_length > 0 ) {
	/* we move along in the real positions by delete number */
	for(k=0;k<delete_length;k++) {
	  temp_real[real_pos] = real;
	  real_pos++;
	  real = real->next->next[0];
	  if( real == NULL ) {
	    return FALSE;
	  } 
	  if( real->next->next_len > 1 ) {
	    fprintf(stderr,"Cannot handle over branched cases in indel area\n");
	    return FALSE;
	  }
	}
      } else {
	/* we move along in the error positions by the delete number */
	delete_length = abs(delete_length);
	for(k=0;k<delete_length;k++) {
	  temp_error[error_pos] = error;
	  error_pos++;
	  error = error->next->next[0];
	  if( error == NULL ) {
	    return FALSE;
	  } 
	  if( error->next->next_len > 1 ) {
	    fprintf(stderr,"Cannot handle over branched cases in indel area\n");
	    return FALSE;
	  }
	}
      }
    }
    

    temp_real[real_pos++] = real;
    temp_error[error_pos++] = error;


    if( real_pos > max_search ) {
      /* will return TRUE after loop */
      break;
    }

    if( real == error ) {
      /* same node, therefore the streams have successfully merged */
      break;
    }


    if( real->next->next_len == 1 && error->next->next_len == 1 ) {
      real = real->next->next[0];
      error = error->next->next[0];
      continue;
    }

    /* this is a tail position */
    if( error->next->next_len == 0 ) {

      /* will return TRUE after loop */
      break;
    }


    /* recursively call for each possible path: if any return
       true, return TRUE */

    /*    fprintf(stderr,"paths are branching, real with %d, error with %d\n",real->next->next_len,error->next->next_len);*/

    if( real->next->next_len > 1 || error->next->next_len > 1 ) {
      for(i=0;i<real->next->next_len;i++) {
	for(j=0;j<error->next->next_len;j++) {
	  (*current_path)++;
	  if( extend_indel_path_KmerAssembly(real->next->next[i],error->next->next[j],real_pos,error_pos,delete_length,proposed_position,final_real,final_error,max_search,current_path,max_path_enum) == TRUE ) {
	    break;
	  }
	}
      }

      /*      fprintf(stderr,"Recursively unable to find path\n");*/
      /* unable to find a good path */

      return FALSE;
    }

    fprintf(stderr,"should be impossible to reach here, all other options have been considered");
    
    return FALSE;
  }


  /* have to put the right path into the path datastructures */

  /* we might have to extend the path datastructures */

  if( final_real->max_stack < real_pos || final_error->max_stack < error_pos ) {
    fatal("Have to extend path; should be easy, but not done yet");
  }

  if( final_real->stack_len < real_pos ) {
    /* first return */
    final_real->stack_len = real_pos;
  }

  if( final_error->stack_len < error_pos ) {
    /* first return */
    final_error->stack_len = error_pos;
  }

  for(i=real_start;i<real_pos;i++) {
    final_real->stack[i]   = temp_real[i];
  }

  for(i=error_start;i<error_pos;i++) {
    final_error->stack[i]   = temp_error[i];
  }

  return TRUE;


}


# line 476 "kmer_assembly_error.dy"
boolean find_error_path_KmerAssembly(KmerAssemblyLink * start_real,KmerAssemblyLink * start_error,char proposed_fix,int proposed_fix_position,KmerAssemblyPath * final_real,KmerAssemblyPath * final_error,int max_search,int max_path_enum)
{
  int current_path = 1;
  return extend_error_path_KmerAssembly(start_real,start_error,0,proposed_fix,proposed_fix_position,final_real,final_error,max_search,&current_path,max_path_enum);
}


# line 483 "kmer_assembly_error.dy"
boolean extend_error_path_KmerAssembly(KmerAssemblyLink * real,KmerAssemblyLink * error,int current_pos,char proposed_fix,int proposed_fix_position,KmerAssemblyPath * final_real,KmerAssemblyPath * final_error,int max_search,int * current_path,int max_path_enum)
{
  int i;
  int j;

  int temp_pos;
  int startpos = current_pos;
  KmerAssemblyLink * temp_real[256];
  KmerAssemblyLink * temp_error[256];

  if( max_search > 256 ) {
    fatal("Madness. Max search greater than 256 - shouldn't happen!");
  }


  if( max_path_enum < *current_path ) {
    fprintf(stderr,"Terminating due too many paths\n");
    return FALSE;
  }

  temp_pos = current_pos;


  /*  fprintf(stderr,"Entering into extension code with %d temp_pos\n",temp_pos);*/

  while( real != NULL && error != NULL ) {

    if( proposed_fix_position != temp_pos ) {
      if( real->base != error->base ) {
	/*	fprintf(stderr,"real (%c) and error (%c) do not agree at position %d\n",real->base,error->base,temp_pos); */
	return FALSE;
      } 
    } else if ( real->base != proposed_fix ) {
      fprintf(stderr,"In fixed position real (%c) and proposed (%c) do not agree\n",real->base,proposed_fix);
      return FALSE;
    }
 

    temp_real[temp_pos] = real;
    temp_error[temp_pos] = error;
    temp_pos++;

    if( temp_pos > max_search ) {
      /* will return TRUE after loop */
      break;
    }

    if( real->next->next_len == 1 && error->next->next_len == 1 ) {
      real = real->next->next[0];
      error = error->next->next[0];
      continue;
    }

    /* this is a tail position */
    if( error->next->next_len == 0 ) {

      /* will return TRUE after loop */
      break;
    }


    /* recursively call for each possible path: if any return
       true, return TRUE */

    /*    fprintf(stderr,"paths are branching, real with %d, error with %d\n",real->next->next_len,error->next->next_len);*/

    if( real->next->next_len > 1 || error->next->next_len > 1 ) {
      for(i=0;i<real->next->next_len;i++) {
	for(j=0;j<error->next->next_len;j++) {
	  (*current_path)++;
	  if( extend_error_path_KmerAssembly(real->next->next[i],error->next->next[j],temp_pos,proposed_fix,proposed_fix_position,final_real,final_error,max_search,current_path,max_path_enum) == TRUE ) {
	    break;
	  }
	}
      }

      /*      fprintf(stderr,"Recursively unable to find path\n");*/
      /* unable to find a good path */

      return FALSE;
    }

    fprintf(stderr,"should be impossible to reach here, all other options have been considered");
    
    return FALSE;
  }


  /* have to put the right path into the path datastructures */

  /* we might have to extend the path datastructures */

  if( final_real->max_stack < temp_pos || final_error->max_stack < temp_pos ) {
    fatal("Have to extend path; should be easy, but not done yet");
  }

  if( final_real->stack_len < temp_pos ) {
    /* first return */
    final_real->stack_len = temp_pos;
    final_error->stack_len = temp_pos;
  }

  for(i=startpos;i<temp_pos;i++) {
    final_real->stack[i]   = temp_real[i];
    final_error->stack[i]  = temp_error[i];
  }

  return TRUE;

}
# line 595 "kmer_assembly_error.c"

#ifdef _cplusplus
}
#endif
