#ifdef _cplusplus
extern "C" {
#endif
#include "eulerindex.h"


# line 54 "eulerindex.dy"
void dump_EulerGraph(EulerGraph * eg,FILE * ofp)
{
  int i;
  int j;
  char dna_str[512];

  assert(eg);
  assert(ofp);

  for(i=0;i<eg->node_len;i++) {
    if( eg->node[i] != NULL ) {
      reverse_map_dna_number(i,eg->kmer,dna_str);
      dna_str[eg->kmer] = '\0';
      fprintf(ofp,"Node %d [%s]\n",i,dna_str);
      fprintf(ofp,"  Incoming Links:\n");
      for(j=0;j<eg->node[i]->back_len;j++) {
	fprintf(ofp,"   %d [%c] from %d\n",eg->node[i]->back[j]->depth,eg->node[i]->back[j]->base,eg->node[i]->back[j]->prev->number);
      }
      fprintf(ofp,"  Outgoing Links:\n");
      for(j=0;j<eg->node[i]->link_len;j++) {
	fprintf(ofp,"   %d [%c] to %d\n",eg->node[i]->link[j]->depth,eg->node[i]->link[j]->base,eg->node[i]->link[j]->next->number);
      }
    }
  }

}

# line 81 "eulerindex.dy"
Sequence * read_Sequence_EulerNode(EulerNode * leftmost)
{
  int i;
  Sequence * out;
  EulerLink * l;

  assert(leftmost);
  assert(leftmost->link_len == 1);
  
  i=0;
  l = leftmost->link[0];

  while( 1 ) {
    i++;
    if( l->next->link_len == 1 ) {
      l = l->next->link[0];
      continue;
    }

    if( l->next->link_len > 1 ) {
      warn("Untangled repeat sequence on node %d\n",l->next->number);
      break;
    }
    if( l->next->link_len == 0 ) {
      break;
    }
  }

  out = Sequence_alloc();
  out->seq = calloc(i+1,sizeof(char));

  i=0;
  l = leftmost->link[0];

  while( 1 ) {
    out->seq[i] = l->base;
    i++;
    if( l->next->link_len == 1 ) {
      l = l->next->link[0];
      continue;
    }

    if( l->next->link_len > 1 ) {
      break;
    }
    if( l->next->link_len == 0 ) {
      break;
    }
  }

  out->seq[i] = '\0';
  out->name = stringalloc("EulerGraphSequence");

  return out;
}


# line 138 "eulerindex.dy"
boolean can_resolve_error_EulerGraph(EulerGraph *eg,EulerLink * leftmost)
{
  EulerNode * node;
  EulerLink * path;

  int len = 1;
  int error_len;

  
  assert(eg);
  assert(leftmost);

  fprintf(stderr,"starting resolving\n");
  error_len = eg->kmer;
  path = leftmost;
  
  while( len < error_len ) {

    if( path->depth != 1 ) {
      return FALSE;
    }

    if( path->next->link_len > 1 ) {
      return FALSE;
    }

    path = path->next->link[0];

    if( path == NULL ) {
      return FALSE;
    }
    len++;
  }

  return TRUE;
}


# line 176 "eulerindex.dy"
boolean resolve_error_EulerGraph(EulerGraph * eg,EulerLink * leftmost)
{
  char * dna = "ATGC";
  char dna_str[128];
  int i;
  int test_number;
  EulerLink * walk;
  int w;

  assert(eg);
  assert(leftmost);

  reverse_map_dna_number(leftmost->next->number,eg->kmer,dna_str);


  for(i=0,walk = leftmost;i<eg->kmer;i++) {
    dna_str[eg->kmer+i] = walk->base;
    walk = walk->next->link[0];
  }
  dna_str[eg->kmer+i] = '\0';

  fprintf(stderr,"Have string of %s\n",dna_str);

  for(i=0;i<4;i++) {
    dna_str[eg->kmer-1] = dna[i];
    
    fprintf(stderr,"Considering %c for fixing\n",dna[i]);

    test_number = forward_dna_number_from_string(dna_str,eg->kmer);
    if( eg->node[test_number] != NULL && eg->node[test_number]->link_len ==1 ) {
      fprintf(stderr,"Test number works, now walking...\n");
      for(walk = eg->node[test_number]->back[0],w=0; walk != NULL && w<eg->kmer;w++) {
	/*	fprintf(stderr,"Looking at %c vs %c\n",dna_str[eg->kmer+w],walk->base); */
	if( dna_str[eg->kmer+w] != walk->base ) {
	  break; 
	}
	walk = walk->next->link[0];
      }
      if( w >= eg->kmer ) {
	fix_error_EulerGraph(eg,leftmost,dna_str,eg->kmer);
	return TRUE;
      }
    }
  }
  
  return FALSE;

}

# line 225 "eulerindex.dy"
void fix_error_EulerGraph(EulerGraph * eg, EulerLink * leftmost,char * dna_str,int len)
{
  int i;
  int j;
  int number;
  int prev_number;
  EulerLink * walk;
  EulerLink * temp;
  EulerNode * prev;

  fprintf(stderr,"Fixing with %s length %d\n",dna_str,len);

  prev = leftmost->prev;
  prev_number = prev->number;

  walk = leftmost;
  for(i=0;i<len;i++) {
    number = forward_dna_number_from_string(dna_str+i,eg->kmer);
    for(j=0;j<eg->node[prev_number]->link_len;j++) {
      if( eg->node[prev_number]->link[j]->next->number == number ) {
	break;
      }
    }
    if( j >= eg->node[prev_number]->link_len ) {
      fprintf(stderr,"Problem here; fixed node doesn't have active numbers in link...");
      return;
    }

    add_label_EulerLink(eg->node[prev_number]->link[j],walk->label[0]);
    prev_number = number;

    temp = walk;
    walk = walk->next->link[0];
    
    remove_EulerLink_forward_EulerNode(temp->prev,temp);
    remove_EulerLink_backward_EulerNode(temp->next,temp);
    free_EulerLink(temp);
  }
    
}

# line 266 "eulerindex.dy"
void build_new_node_path_EulerGraph(EulerGraph * eg,EulerLink * leftmost,EulerPath * path,int * starting_labels,int length)
{
  int i;
  int l;
  int path_offset;
  EulerNode * new_node;
  EulerNode * prev_node;

  EulerLink * new_link;

  assert(eg);
  assert(leftmost);
  assert(starting_labels);
  assert(length > 0);
  
  fprintf(stderr,"Going to build new node path\n");

  new_node = new_EulerNode(leftmost->next->number);
  add_dup_EulerGraph(eg,new_node);

  
  new_link = new_EulerLink();
  new_link->prev = leftmost->prev;
  new_link->next = new_node;
  new_link->base = leftmost->base;
  new_link->depth = length;

  for(l=0;l<length;l++) {
    remove_label_EulerLink(leftmost,starting_labels[l]);
    add_label_EulerLink(new_link,starting_labels[l]);
  }

  add_link_EulerNode(new_link->prev,new_link);
  add_back_EulerNode(new_link->next,new_link);

  prev_node = new_node;

  for(i=path->current-1,path_offset = 1;i >= 0 ;i--,path_offset++) {

    new_link = new_EulerLink();
    new_link->prev = prev_node;
    new_link->next = NULL;
    new_link->base = path->stack[i]->base;
    new_link->depth = length;

    for(l=0;l<length;l++) {
      remove_label_EulerLink(path->stack[i],starting_labels[i]+path_offset);
      add_label_EulerLink(new_link,starting_labels[l]+path_offset);
    }

    if( i > 0 ) {
      new_node = new_EulerNode(path->stack[i]->next->number);
      add_dup_EulerGraph(eg,new_node);
      new_link->next = new_node;
    } else {
      /* connect to the same path */
      new_link->next = path->stack[0]->next;
    }

    add_link_EulerNode(new_link->prev,new_link);
    add_back_EulerNode(new_link->next,new_link);
    
    
    prev_node = new_node;
   
  }



}


# line 338 "eulerindex.dy"
boolean attempt_untangle_EulerPath(EulerGraph *eg,EulerPath * path,EulerLink * leftmost)
{
  int i;
  int j;
  SinglePosSequence * leftmost_sps[512];
  SinglePosSequence * rightmost_sps[512];
  int can_untangle_left[512];
  int can_untangle_right[512];
  int untangle_feasible = 0;

  int starting_label[512];
  int total_labels = 0;

  int fully_left_untangled;
  int fully_right_untangled;

  assert(eg);
  assert(path);
  assert(leftmost);

  fprintf(stderr,"Considering untangle at position %d in path\n",path->current);

  for(i=0;i<leftmost->label_len && leftmost->label[i] != -1;i++) {
    leftmost_sps[i]      = lookup_Sequence_SinglePosSpace(eg->sps,leftmost->label[i]);
    if( eg->kmer+leftmost->label[i]+3 > leftmost_sps[i]->end ) {
      can_untangle_left[i] = -2;
      fprintf(stderr,"Too close to end for resolving power, %s\n",leftmost_sps[i]->seq->name);
    } else {
      can_untangle_left[i] = -1;
    }
  }

  fprintf(stderr,"Going to look at rightmost\n");

  for(i=0;i<path->stack[0]->label_len && path->stack[0]->label[i] != -1;i++) {
    rightmost_sps[i]      = lookup_Sequence_SinglePosSpace(eg->sps,path->stack[0]->label[i]);

    if( path->stack[0]->label[i] - eg->kmer - 3 < rightmost_sps[i]->start ) {
      can_untangle_right[i] = -2;
      fprintf(stderr,"Too close to start for resolving power, %s\n",rightmost_sps[i]->seq->name);
    } else {
      can_untangle_right[i] = -1;
    }

    can_untangle_right[i] = -1;
    for(j=0;j<leftmost->label_len && leftmost->label[j] != -1;j++) {
      if( leftmost_sps[j]->seq == rightmost_sps[i]->seq && (path->current == (path->stack[0]->label[i] - leftmost->label[j])) ) {

      	fprintf(stderr,"Able to untangle label %d with (%d) diff [left %d, right %d]\n",leftmost->label[j],path->stack[0]->label[i]-leftmost->label[j],j,i);
	can_untangle_left[j] = i;
	can_untangle_right[i] = j;
	untangle_feasible = 1;

	starting_label[total_labels++] = leftmost->label[j];
      }
    }
  }
  
  if( untangle_feasible == 0 ) {
    fprintf(stderr,"No untangle feasible\n");
    return 0;
  }

  fprintf(stderr,"Abotu to build new node path....\n");

  build_new_node_path_EulerGraph(eg,leftmost,path,starting_label,total_labels);

  fully_left_untangled = 1;
  for(i=0;i<leftmost->label_len && leftmost->label[i] != -1 ;i++) {
    if( can_untangle_left[i] == -1 ) {
      fprintf(stderr,"Left: Label position %d not resolved (%s)\n",i,leftmost_sps[i]->seq->name);
      fully_left_untangled = 0;
      break;
    }
  }

  fully_right_untangled = 1;
  for(i=0;i<path->stack[0]->label_len && path->stack[0]->label[i] != -1;i++) {
    if( can_untangle_right[i] == -1 ) {
      fprintf(stderr,"Right: Label position %d not resolved (%s)\n",i,rightmost_sps[i]->seq->name);
      fully_right_untangled = 0;
      break;
    }
  }

  if( fully_left_untangled == 1) {
    fprintf(stderr,"Managed to fully left untangle path\n");
    /* remove the link on leftmost */
    remove_EulerLink_forward_EulerNode(leftmost->prev,leftmost);
    remove_EulerLink_backward_EulerNode(leftmost->next,leftmost);
  }

  if( fully_right_untangled == 1 ) {
    fprintf(stderr,"Managed to fully right untangle path\n");
    /* remove the link on leftmost */
    remove_EulerLink_forward_EulerNode(path->stack[0]->prev,path->stack[0]);
    remove_EulerLink_backward_EulerNode(path->stack[0]->next,path->stack[0]);
  }

  return fully_right_untangled;
}

# line 440 "eulerindex.dy"
boolean untangle_from_split_EulerNode(EulerGraph * eg,EulerNode * split_outgoing,int max_backtrack)
{
  EulerPath * ep;
  int i;
  int resolved = 0;
  EulerLink * split[512];
  int len;

  assert(eg);
  assert(split_outgoing);
  assert(split_outgoing->link_len > 1);

  ep = new_EulerPath();

  for(i=0;i<split_outgoing->link_len;i++) {
    split[i] = split_outgoing->link[i];
  }
  len = split_outgoing->link_len;

  for(i=0;i<len;i++) {
    if( split_outgoing->link_len == 1 ) {
      fprintf(stderr,"No need to resolve this node at %d from %d\n",i,len);
      break;
    }

    ep->current = 0;
    resolved = 0;
    untangle_EulerLink_EulerPath(eg,ep,split[i],&resolved,0,max_backtrack);
  }

  free_EulerPath(ep);

}

# line 474 "eulerindex.dy"
boolean untangle_EulerLink_EulerPath(EulerGraph * eg,EulerPath * current_path,EulerLink * current,int * resolved,int backtrack_len,int max_backtrack)
{
  int i;
  int starting_path_point;

  assert(eg);
  assert(current_path);
  assert(current);


  fprintf(stderr,"Entering untangle at link between %d and %d, backtrack length of %d on path of %d\n",current->prev->number,current->next->number,backtrack_len,current_path->current);

  if( *resolved == 1 ) {
    return TRUE;
  }

  if( backtrack_len >= max_backtrack ) {
    return TRUE;
  }

  starting_path_point = current_path->current;

  while( 1 ) {
    /* push current into the path */
    fprintf(stderr,"Pushing on to path link between %d and %d with back length of %d\n",current->prev->number,current->next->number,current->prev->back_len);

    push_EulerPath(current_path,current);
    backtrack_len++;

    if( current->prev->back_len == 0 ) {
      /* end of a stream! */
      *resolved = 1;
      return TRUE;
    }


    if( current->prev->back_len > 1 ) {
      for(i=0;i<current->prev->back_len;i++) {
	if( attempt_untangle_EulerPath(eg,current_path,current->prev->back[i]) == 1 ) {
	  *resolved = 1;
	  return TRUE;
	}
      }

      for(i=0;i<current->prev->back_len;i++) {
	/* recurse into this branch */
	untangle_EulerLink_EulerPath(eg,current_path,current->prev->back[i],resolved,backtrack_len,max_backtrack);
      }
      
      current_path->current = starting_path_point;
      /* at this point, gone down all paths, so return */
      return TRUE;
    } else {
      /* continue in this routine steping back */
      current = current->prev->back[0];
    }
  }

  current_path->current = starting_path_point;
  return TRUE;

}

# line 537 "eulerindex.dy"
boolean remove_EulerLink_forward_EulerNode(EulerNode * n,EulerLink *l)
{
  int i;

  assert(n);
  assert(l);

  for(i=0;i<n->link_len;i++) {
    if( n->link[i] == l ) {
      for(++i;i<n->link_len;i++) {
	n->link[i-1] = n->link[i];
      }
      n->link_len--;
      return TRUE;
    }
  }

  return FALSE;
}


# line 558 "eulerindex.dy"
boolean remove_EulerLink_backward_EulerNode(EulerNode * n,EulerLink *l)
{
  int i;

  assert(n);
  assert(l);

  for(i=0;i<n->back_len;i++) {
    if( n->back[i] == l ) {
      for(++i;i<n->back_len;i++) {
	n->back[i-1] = n->back[i];
      }
      n->back_len--;
      return TRUE;
    }
  }

  fprintf(stderr,"...could not remove link\n");
  return FALSE;
}



# line 581 "eulerindex.dy"
boolean remove_label_EulerLink(EulerLink * el,int label)
{
  int i;

  assert(el);

  for(i=0;i<el->label_len;i++) {
    if( el->label[i] == label ) {
      break;
    }
  }

  if( i >= el->label_len ) {
    return FALSE;
  }


  for(i++;i<el->label_len && el->label[i] != -1 ;i++) {
    el->label[i-1] = el->label[i];
  }

  el->label[i] = -1;

  return TRUE;
}


# line 608 "eulerindex.dy"
boolean store_Sequence_EulerGraph(EulerGraph * eg,Sequence * seq)
{
  int i;
  int prev_number;
  int next_number;
  
  int pos_start;

  assert(eg);
  assert(seq);

  add_EulerGraph(eg,seq);

  pos_start = add_Sequence_SinglePosSpace(eg->sps,seq);

  for(i=0;i<seq->len - eg->kmer ;i++) {
    prev_number = forward_dna_number_from_string(seq->seq+i,eg->kmer);
    next_number = forward_dna_number_from_string(seq->seq+i+1,eg->kmer);

    store_EulerLink_EulerGraph(eg,prev_number,next_number,seq->seq[i+eg->kmer],pos_start+i);
  }

}

# line 632 "eulerindex.dy"
boolean store_EulerLink_EulerGraph(EulerGraph * eg,int prev_number,int next_number,char base,int label)
{
  int i;
  EulerLink * link;
  
  assert(eg);

  if( eg->node[prev_number] != NULL ) {
    for(i=0;i<eg->node[prev_number]->link_len;i++) {
      if( eg->node[prev_number]->link[i]->next->number == next_number ) {
	eg->node[prev_number]->link[i]->depth++;
	add_label_EulerLink(eg->node[prev_number]->link[i],label);
	return TRUE;
      }
    }
  }

  /* else need new link */

  if( eg->node[prev_number] == NULL ) {
    eg->node[prev_number] = new_EulerNode(prev_number);
  }

  if( eg->node[next_number] == NULL ) {
    eg->node[next_number] = new_EulerNode(next_number);
  }


  link = new_EulerLink();
  link->depth = 1;
  link->base = base;
  add_label_EulerLink(link,label);

  link->prev = eg->node[prev_number];
  link->next = eg->node[next_number];
  
  eg->node[prev_number]->is_rightmost = 0;
  eg->node[next_number]->is_leftmost = 0;

  add_link_EulerNode(eg->node[prev_number],link);
  add_back_EulerNode(eg->node[next_number],link);
  
  return TRUE;
}


# line 678 "eulerindex.dy"
EulerLink * new_EulerLink(void)
{
  int i;
  EulerLink * out;

  out = EulerLink_alloc();
  out->label = calloc(EulerLinkLabelLength,sizeof(int));
  for(i=0;i<EulerLinkLabelLength;i++) {
    out->label[i] = -1;
  }
  out->label_len = EulerLinkLabelLength;

  return out;
}

# line 693 "eulerindex.dy"
boolean add_label_EulerLink(EulerLink * el,int label)
{
  int i;
  int old_len;

  assert(el);

  for(i=0;i<el->label_len;i++) {
    if( el->label[i] == -1 ) {
      el->label[i] = label;
      return TRUE;
    }
  }
  /* run out of space */
  old_len = el->label_len;

  if( el->label_len > EulerLinkLabelLinear ) {
    el->label = realloc(el->label,el->label_len*2*sizeof(int));
    el->label_len = el->label_len*2;
  } else {
    el->label = realloc(el->label,el->label_len*el->label_len*sizeof(int));
    el->label_len = el->label_len * el->label_len;
  }

  /* put label at old_len, and then -1s */
  el->label[old_len] = label;
  for(i=old_len+1;i<el->label_len;i++) {
    el->label[i] = -1;
  }

  return TRUE;
}

# line 726 "eulerindex.dy"
EulerNode * new_EulerNode(int number)
{
  EulerNode * out;

  out = EulerNode_alloc_std();

  out->number = number;
  out->is_leftmost = 1;
  out->is_rightmost = 1;

  return out;
}

# line 739 "eulerindex.dy"
EulerGraph * new_EulerGraph(int kmer)
{
  EulerGraph * out;
  int len = 1;
  int i;

  for(i=0;i<kmer;i++) {
    len = len *4;
  }

  out = EulerGraph_alloc_std();
  out->node = calloc(len,sizeof(EulerNode*));

  out->node_len = len;
  out->kmer = kmer;

  out->sps = new_SinglePosSpace(0,3000);
  
  return out;
}

# line 760 "eulerindex.dy"
void * push_EulerPath(EulerPath * ep,EulerLink * el)
{
  assert(ep);
  assert(el);

  if( ep->current+1 >= ep->max_stack_len ) 
    extend_EulerPath_stack(ep);

  ep->stack[ep->current++] = el;

  return;
}


# line 774 "eulerindex.dy"
EulerLink * pop_EulerPath(EulerPath * ep)
{
  assert(ep);

  if( ep->current <= 0 ) {
    warn("Stack underflow on EulerPath");
    return NULL;
  }

  return ep->stack[ep->current--];
}

# line 786 "eulerindex.dy"
void extend_EulerPath_stack(EulerPath * ep)
{
  assert(ep);

  ep->stack = (EulerLink**) realloc(ep->stack,sizeof(EulerLink*) * ep->max_stack_len * 2);
  ep->max_stack_len = ep->max_stack_len * 2;

  return;
}

# line 796 "eulerindex.dy"
EulerPath * new_EulerPath(void)
{
  EulerPath * out;

  out = EulerPath_alloc();
  out->stack = (EulerLink**) calloc(64,sizeof(EulerLink*));
  out->max_stack_len = 64;
  out->current = 0;

  return out;
}

# line 782 "eulerindex.c"
/* Function:  hard_link_EulerLink(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [EulerLink *]
 *
 * Return [UNKN ]  Undocumented return value [EulerLink *]
 *
 */
EulerLink * hard_link_EulerLink(EulerLink * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a EulerLink object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  EulerLink_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [EulerLink *]
 *
 */
EulerLink * EulerLink_alloc(void) 
{
    EulerLink * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(EulerLink *) ckalloc (sizeof(EulerLink))) == NULL)  {  
      warn("EulerLink_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->label = NULL;   
    out->label_len = 0;  
    out->depth = 0;  
    out->base = 'u'; 


    return out;  
}    


/* Function:  free_EulerLink(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [EulerLink *]
 *
 * Return [UNKN ]  Undocumented return value [EulerLink *]
 *
 */
EulerLink * free_EulerLink(EulerLink * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a EulerLink obj. Should be trappable"); 
      return NULL;   
      }  


#ifdef PTHREAD   
    assert(pthread_mutex_lock(&(obj->dynamite_mutex)) == 0); 
#endif   
    if( obj->dynamite_hard_link > 1)     {  
      return_early = 1;  
      obj->dynamite_hard_link--; 
      }  
#ifdef PTHREAD   
    assert(pthread_mutex_unlock(&(obj->dynamite_mutex)) == 0);   
#endif   
    if( return_early == 1)   
      return NULL;   
    /* obj->prev is linked in */ 
    /* obj->next is linked in */ 
    if( obj->label != NULL)  
      ckfree(obj->label);    


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_EulerPath(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [EulerPath *]
 *
 * Return [UNKN ]  Undocumented return value [EulerPath *]
 *
 */
EulerPath * hard_link_EulerPath(EulerPath * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a EulerPath object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  EulerPath_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [EulerPath *]
 *
 */
EulerPath * EulerPath_alloc(void) 
{
    EulerPath * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(EulerPath *) ckalloc (sizeof(EulerPath))) == NULL)  {  
      warn("EulerPath_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->stack = NULL;   
    out->max_stack_len = 0;  
    out->current = 0;    


    return out;  
}    


/* Function:  free_EulerPath(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [EulerPath *]
 *
 * Return [UNKN ]  Undocumented return value [EulerPath *]
 *
 */
EulerPath * free_EulerPath(EulerPath * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a EulerPath obj. Should be trappable"); 
      return NULL;   
      }  


#ifdef PTHREAD   
    assert(pthread_mutex_lock(&(obj->dynamite_mutex)) == 0); 
#endif   
    if( obj->dynamite_hard_link > 1)     {  
      return_early = 1;  
      obj->dynamite_hard_link--; 
      }  
#ifdef PTHREAD   
    assert(pthread_mutex_unlock(&(obj->dynamite_mutex)) == 0);   
#endif   
    if( return_early == 1)   
      return NULL;   
    if( obj->stack != NULL)  
      ckfree(obj->stack);    


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_link_EulerNode(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_link_EulerNode
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [EulerLink **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_link_EulerNode(EulerLink ** list,int i,int j)  
{
    EulerLink * temp;    
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_link_EulerNode(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_link_EulerNode which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [EulerLink **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_link_EulerNode(EulerLink ** list,int left,int right,int (*comp)(EulerLink * ,EulerLink * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_link_EulerNode(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_link_EulerNode (list,++last,i); 
      }  
    swap_link_EulerNode (list,left,last);    
    qsort_link_EulerNode(list,left,last-1,comp); 
    qsort_link_EulerNode(list,last+1,right,comp);    
}    


/* Function:  sort_link_EulerNode(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_link_EulerNode
 *
 *
 * Arg:         obj [UNKN ] Object containing list [EulerNode *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_link_EulerNode(EulerNode * obj,int (*comp)(EulerLink *, EulerLink *)) 
{
    qsort_link_EulerNode(obj->link,0,obj->link_len-1,comp);  
    return;  
}    


/* Function:  expand_link_EulerNode(obj,len)
 *
 * Descrip:    Really an internal function for add_link_EulerNode
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [EulerNode *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_link_EulerNode(EulerNode * obj,int len) 
{


    if( obj->link_maxlen > obj->link_len )   {  
      warn("expand_EulerNodelink_ called with no need"); 
      return TRUE;   
      }  


    if( (obj->link = (EulerLink ** ) ckrealloc (obj->link,sizeof(EulerLink *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_EulerNode, returning FALSE");    
      return FALSE;  
      }  
    obj->link_maxlen = len;  
    return TRUE; 
}    


/* Function:  add_link_EulerNode(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [EulerNode *]
 * Arg:        add [OWNER] Object to add to the list [EulerLink *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_link_EulerNode(EulerNode * obj,EulerLink * add) 
{
    if( obj->link_len >= obj->link_maxlen)   {  
      if( expand_link_EulerNode(obj,obj->link_len + EulerNodeLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->link[obj->link_len++]=add;  
    return TRUE; 
}    


/* Function:  flush_link_EulerNode(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [EulerNode *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_link_EulerNode(EulerNode * obj) 
{
    int i;   


    for(i=0;i<obj->link_len;i++) { /*for i over list length*/ 
      if( obj->link[i] != NULL)  {  
        free_EulerLink(obj->link[i]);    
        obj->link[i] = NULL; 
        }  
      } /* end of for i over list length */ 


    obj->link_len = 0;   
    return i;    
}    


/* Function:  swap_back_EulerNode(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_back_EulerNode
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [EulerLink **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_back_EulerNode(EulerLink ** list,int i,int j)  
{
    EulerLink * temp;    
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_back_EulerNode(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_back_EulerNode which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [EulerLink **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_back_EulerNode(EulerLink ** list,int left,int right,int (*comp)(EulerLink * ,EulerLink * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_back_EulerNode(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_back_EulerNode (list,++last,i); 
      }  
    swap_back_EulerNode (list,left,last);    
    qsort_back_EulerNode(list,left,last-1,comp); 
    qsort_back_EulerNode(list,last+1,right,comp);    
}    


/* Function:  sort_back_EulerNode(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_back_EulerNode
 *
 *
 * Arg:         obj [UNKN ] Object containing list [EulerNode *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_back_EulerNode(EulerNode * obj,int (*comp)(EulerLink *, EulerLink *)) 
{
    qsort_back_EulerNode(obj->back,0,obj->back_len-1,comp);  
    return;  
}    


/* Function:  expand_back_EulerNode(obj,len)
 *
 * Descrip:    Really an internal function for add_back_EulerNode
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [EulerNode *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_back_EulerNode(EulerNode * obj,int len) 
{


    if( obj->back_maxlen > obj->back_len )   {  
      warn("expand_EulerNodeback_ called with no need"); 
      return TRUE;   
      }  


    if( (obj->back = (EulerLink ** ) ckrealloc (obj->back,sizeof(EulerLink *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_EulerNode, returning FALSE");    
      return FALSE;  
      }  
    obj->back_maxlen = len;  
    return TRUE; 
}    


/* Function:  add_back_EulerNode(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [EulerNode *]
 * Arg:        add [OWNER] Object to add to the list [EulerLink *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_back_EulerNode(EulerNode * obj,EulerLink * add) 
{
    if( obj->back_len >= obj->back_maxlen)   {  
      if( expand_back_EulerNode(obj,obj->back_len + EulerNodeLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->back[obj->back_len++]=add;  
    return TRUE; 
}    


/* Function:  flush_back_EulerNode(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [EulerNode *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_back_EulerNode(EulerNode * obj) 
{
    int i;   


    for(i=0;i<obj->back_len;i++) { /*for i over list length*/ 
      if( obj->back[i] != NULL)  {  
        free_EulerLink(obj->back[i]);    
        obj->back[i] = NULL; 
        }  
      } /* end of for i over list length */ 


    obj->back_len = 0;   
    return i;    
}    


/* Function:  EulerNode_alloc_std(void)
 *
 * Descrip:    Equivalent to EulerNode_alloc_len(EulerNodeLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [EulerNode *]
 *
 */
EulerNode * EulerNode_alloc_std(void) 
{
    return EulerNode_alloc_len(EulerNodeLISTLENGTH); 
}    


/* Function:  EulerNode_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [EulerNode *]
 *
 */
EulerNode * EulerNode_alloc_len(int len) 
{
    EulerNode * out;/* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = EulerNode_alloc()) == NULL)    
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->link = (EulerLink ** ) ckcalloc (len,sizeof(EulerLink *))) == NULL) {  
      warn("Warning, ckcalloc failed in EulerNode_alloc_len");   
      return NULL;   
      }  
    out->link_len = 0;   
    out->link_maxlen = len;  


    if((out->back = (EulerLink ** ) ckcalloc (len,sizeof(EulerLink *))) == NULL) {  
      warn("Warning, ckcalloc failed in EulerNode_alloc_len");   
      return NULL;   
      }  
    out->back_len = 0;   
    out->back_maxlen = len;  


    return out;  
}    


/* Function:  hard_link_EulerNode(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [EulerNode *]
 *
 * Return [UNKN ]  Undocumented return value [EulerNode *]
 *
 */
EulerNode * hard_link_EulerNode(EulerNode * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a EulerNode object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  EulerNode_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [EulerNode *]
 *
 */
EulerNode * EulerNode_alloc(void) 
{
    EulerNode * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(EulerNode *) ckalloc (sizeof(EulerNode))) == NULL)  {  
      warn("EulerNode_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->number = 0; 
    out->link = NULL;    
    out->link_len = out->link_maxlen = 0;    
    out->back = NULL;    
    out->back_len = out->back_maxlen = 0;    
    out->is_leftmost = 'u';  
    out->is_rightmost = 'u'; 


    return out;  
}    


/* Function:  free_EulerNode(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [EulerNode *]
 *
 * Return [UNKN ]  Undocumented return value [EulerNode *]
 *
 */
EulerNode * free_EulerNode(EulerNode * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a EulerNode obj. Should be trappable"); 
      return NULL;   
      }  


#ifdef PTHREAD   
    assert(pthread_mutex_lock(&(obj->dynamite_mutex)) == 0); 
#endif   
    if( obj->dynamite_hard_link > 1)     {  
      return_early = 1;  
      obj->dynamite_hard_link--; 
      }  
#ifdef PTHREAD   
    assert(pthread_mutex_unlock(&(obj->dynamite_mutex)) == 0);   
#endif   
    if( return_early == 1)   
      return NULL;   
    if( obj->link != NULL)   {  
      for(i=0;i<obj->link_len;i++)   {  
        if( obj->link[i] != NULL)    
          free_EulerLink(obj->link[i]);  
        }  
      ckfree(obj->link); 
      }  
    if( obj->back != NULL)   {  
      for(i=0;i<obj->back_len;i++)   {  
        if( obj->back[i] != NULL)    
          free_EulerLink(obj->back[i]);  
        }  
      ckfree(obj->back); 
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_EulerGraph(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_EulerGraph
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [Sequence **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_EulerGraph(Sequence ** list,int i,int j)  
{
    Sequence * temp; 
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_EulerGraph(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_EulerGraph which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [Sequence **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_EulerGraph(Sequence ** list,int left,int right,int (*comp)(Sequence * ,Sequence * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_EulerGraph(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_EulerGraph (list,++last,i); 
      }  
    swap_EulerGraph (list,left,last);    
    qsort_EulerGraph(list,left,last-1,comp); 
    qsort_EulerGraph(list,last+1,right,comp);    
}    


/* Function:  sort_EulerGraph(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_EulerGraph
 *
 *
 * Arg:         obj [UNKN ] Object containing list [EulerGraph *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_EulerGraph(EulerGraph * obj,int (*comp)(Sequence *, Sequence *)) 
{
    qsort_EulerGraph(obj->seq,0,obj->len-1,comp);    
    return;  
}    


/* Function:  expand_EulerGraph(obj,len)
 *
 * Descrip:    Really an internal function for add_EulerGraph
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [EulerGraph *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_EulerGraph(EulerGraph * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_EulerGraph called with no need"); 
      return TRUE;   
      }  


    if( (obj->seq = (Sequence ** ) ckrealloc (obj->seq,sizeof(Sequence *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_EulerGraph, returning FALSE");   
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_EulerGraph(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [EulerGraph *]
 * Arg:        add [OWNER] Object to add to the list [Sequence *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_EulerGraph(EulerGraph * obj,Sequence * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_EulerGraph(obj,obj->len + EulerGraphLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->seq[obj->len++]=add;    
    return TRUE; 
}    


/* Function:  flush_EulerGraph(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [EulerGraph *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_EulerGraph(EulerGraph * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->seq[i] != NULL)   {  
        free_Sequence(obj->seq[i]);  
        obj->seq[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  swap_dup_EulerGraph(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_dup_EulerGraph
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [EulerNode **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_dup_EulerGraph(EulerNode ** list,int i,int j)  
{
    EulerNode * temp;    
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_dup_EulerGraph(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_dup_EulerGraph which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [EulerNode **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_dup_EulerGraph(EulerNode ** list,int left,int right,int (*comp)(EulerNode * ,EulerNode * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_dup_EulerGraph(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_dup_EulerGraph (list,++last,i); 
      }  
    swap_dup_EulerGraph (list,left,last);    
    qsort_dup_EulerGraph(list,left,last-1,comp); 
    qsort_dup_EulerGraph(list,last+1,right,comp);    
}    


/* Function:  sort_dup_EulerGraph(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_dup_EulerGraph
 *
 *
 * Arg:         obj [UNKN ] Object containing list [EulerGraph *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_dup_EulerGraph(EulerGraph * obj,int (*comp)(EulerNode *, EulerNode *)) 
{
    qsort_dup_EulerGraph(obj->dup,0,obj->dup_len-1,comp);    
    return;  
}    


/* Function:  expand_dup_EulerGraph(obj,len)
 *
 * Descrip:    Really an internal function for add_dup_EulerGraph
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [EulerGraph *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_dup_EulerGraph(EulerGraph * obj,int len) 
{


    if( obj->dup_maxlen > obj->dup_len )     {  
      warn("expand_EulerGraphdup_ called with no need"); 
      return TRUE;   
      }  


    if( (obj->dup = (EulerNode ** ) ckrealloc (obj->dup,sizeof(EulerNode *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_EulerGraph, returning FALSE");   
      return FALSE;  
      }  
    obj->dup_maxlen = len;   
    return TRUE; 
}    


/* Function:  add_dup_EulerGraph(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [EulerGraph *]
 * Arg:        add [OWNER] Object to add to the list [EulerNode *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_dup_EulerGraph(EulerGraph * obj,EulerNode * add) 
{
    if( obj->dup_len >= obj->dup_maxlen) {  
      if( expand_dup_EulerGraph(obj,obj->dup_len + EulerGraphLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->dup[obj->dup_len++]=add;    
    return TRUE; 
}    


/* Function:  flush_dup_EulerGraph(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [EulerGraph *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_dup_EulerGraph(EulerGraph * obj) 
{
    int i;   


    for(i=0;i<obj->dup_len;i++)  { /*for i over list length*/ 
      if( obj->dup[i] != NULL)   {  
        free_EulerNode(obj->dup[i]); 
        obj->dup[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->dup_len = 0;    
    return i;    
}    


/* Function:  EulerGraph_alloc_std(void)
 *
 * Descrip:    Equivalent to EulerGraph_alloc_len(EulerGraphLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [EulerGraph *]
 *
 */
EulerGraph * EulerGraph_alloc_std(void) 
{
    return EulerGraph_alloc_len(EulerGraphLISTLENGTH);   
}    


/* Function:  EulerGraph_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [EulerGraph *]
 *
 */
EulerGraph * EulerGraph_alloc_len(int len) 
{
    EulerGraph * out;   /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = EulerGraph_alloc()) == NULL)   
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->seq = (Sequence ** ) ckcalloc (len,sizeof(Sequence *))) == NULL)    {  
      warn("Warning, ckcalloc failed in EulerGraph_alloc_len");  
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    if((out->dup = (EulerNode ** ) ckcalloc (len,sizeof(EulerNode *))) == NULL)  {  
      warn("Warning, ckcalloc failed in EulerGraph_alloc_len");  
      return NULL;   
      }  
    out->dup_len = 0;    
    out->dup_maxlen = len;   


    return out;  
}    


/* Function:  hard_link_EulerGraph(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [EulerGraph *]
 *
 * Return [UNKN ]  Undocumented return value [EulerGraph *]
 *
 */
EulerGraph * hard_link_EulerGraph(EulerGraph * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a EulerGraph object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  EulerGraph_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [EulerGraph *]
 *
 */
EulerGraph * EulerGraph_alloc(void) 
{
    EulerGraph * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(EulerGraph *) ckalloc (sizeof(EulerGraph))) == NULL)    {  
      warn("EulerGraph_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->node = NULL;    
    out->node_len = 0;   
    out->kmer = 0;   
    out->seq = NULL; 
    out->len = out->maxlen = 0;  
    out->sps = NULL; 
    out->dup = NULL; 
    out->dup_len = out->dup_maxlen = 0;  


    return out;  
}    


/* Function:  free_EulerGraph(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [EulerGraph *]
 *
 * Return [UNKN ]  Undocumented return value [EulerGraph *]
 *
 */
EulerGraph * free_EulerGraph(EulerGraph * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a EulerGraph obj. Should be trappable");    
      return NULL;   
      }  


#ifdef PTHREAD   
    assert(pthread_mutex_lock(&(obj->dynamite_mutex)) == 0); 
#endif   
    if( obj->dynamite_hard_link > 1)     {  
      return_early = 1;  
      obj->dynamite_hard_link--; 
      }  
#ifdef PTHREAD   
    assert(pthread_mutex_unlock(&(obj->dynamite_mutex)) == 0);   
#endif   
    if( return_early == 1)   
      return NULL;   
    if( obj->node != NULL)   
      ckfree(obj->node);     
    if( obj->seq != NULL)    {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->seq[i] != NULL) 
          free_Sequence(obj->seq[i]);    
        }  
      ckfree(obj->seq);  
      }  
    if( obj->sps != NULL)    
      free_SinglePosSpace(obj->sps);     
    if( obj->dup != NULL)    {  
      for(i=0;i<obj->dup_len;i++)    {  
        if( obj->dup[i] != NULL) 
          free_EulerNode(obj->dup[i]);   
        }  
      ckfree(obj->dup);  
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_EulerErrorPara(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [EulerErrorPara *]
 *
 * Return [UNKN ]  Undocumented return value [EulerErrorPara *]
 *
 */
EulerErrorPara * hard_link_EulerErrorPara(EulerErrorPara * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a EulerErrorPara object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  EulerErrorPara_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [EulerErrorPara *]
 *
 */
EulerErrorPara * EulerErrorPara_alloc(void) 
{
    EulerErrorPara * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(EulerErrorPara *) ckalloc (sizeof(EulerErrorPara))) == NULL)    {  
      warn("EulerErrorPara_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->error_size = 1; 


    return out;  
}    


/* Function:  free_EulerErrorPara(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [EulerErrorPara *]
 *
 * Return [UNKN ]  Undocumented return value [EulerErrorPara *]
 *
 */
EulerErrorPara * free_EulerErrorPara(EulerErrorPara * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a EulerErrorPara obj. Should be trappable");    
      return NULL;   
      }  


#ifdef PTHREAD   
    assert(pthread_mutex_lock(&(obj->dynamite_mutex)) == 0); 
#endif   
    if( obj->dynamite_hard_link > 1)     {  
      return_early = 1;  
      obj->dynamite_hard_link--; 
      }  
#ifdef PTHREAD   
    assert(pthread_mutex_unlock(&(obj->dynamite_mutex)) == 0);   
#endif   
    if( return_early == 1)   
      return NULL;   


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
