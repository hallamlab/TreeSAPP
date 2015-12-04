#ifdef _cplusplus
extern "C" {
#endif
#include "linkindex.h"

/* Function:  extract_dna_LinkStream(ln,lnad,nmer_size)
 *
 * Descrip:    Reverse builds a DNA sequence stream from LinkStream
 *
 *
 * Arg:               ln [UNKN ] Undocumented argument [LinkStream *]
 * Arg:             lnad [UNKN ] Undocumented argument [LinkNumberArrayDebug *]
 * Arg:        nmer_size [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [Sequence *]
 *
 */
# line 56 "linkindex.dy"
Sequence * extract_dna_LinkStream(LinkStream * ln,LinkNumberArrayDebug * lnad,int nmer_size)
{
  LinkStream * runner;
  LinkStream * prev;
  LinkStream * next;
  int i = 0;
  Sequence * out;
  int flipped = 0; 
  
  assert(ln);
  i++;
  
  if( lnad->extraction != 0 ) {
    fprintf(lnad->ofp,"Extracting DNA from linkstream %d\n",ln);
  }

  ln->have_seen = 1;

  if( ln->a == NULL ) {
    if( ln->b->x == ln ) {
      runner = ln->b->y;
    } else {
      runner = ln->b->x;
    }
  } else {
    if( ln->a->x == ln ) {
      runner = ln->a->y;
    } else {
      runner = ln->a->x;
    }
  }

  prev = ln;

  while( runner != NULL ) {
    /* to find the outgoing link from here, test neither a nor b
       is not NULL and figure out the right way to go from the
       fact that the link is back to ourselves */

    if( runner->a == NULL || runner->b == NULL ) {
      /* other end of stream */
      break;
    }
    runner->have_seen = 1;
    i++;
    if( lnad->extraction > 2 ) {
      fprintf(lnad->ofp,"Extracting DNA from linkstream %d, runner %d, position count %d\n",ln,runner,i);
    }
    
    if( runner->a->x == runner && runner->a->y != prev) {
      next = runner->a->y;
    } else if ( runner->a->y == runner && runner->a->x != prev ) {
      next = runner->a->x;
    } else if( runner->b->x == runner && runner->b->y != prev) {
      next = runner->b->y;
    } else if ( runner->b->y == runner && runner->b->x != prev ) {
      next = runner->b->x;
    } else {
      fatal("Unable to move off edge!");
    }

    prev = runner;
    runner = next;
  }


  out = Sequence_alloc();
  out->seq = calloc(i+1,sizeof(char));

  i = 0;
  flipped = ln->starting_flip;
  if( ln->a == NULL ) {
    out->seq[0] = first_char_from_dnanumber(ln->number,nmer_size,flipped);
  } else {
    /* b is NULL, indicating a 3' end of a sequence, so flip the flip*/
    out->seq[0] = first_char_from_dnanumber(ln->number,nmer_size,!flipped);
  }

  if( ln->a == NULL ) {
    if( ln->b->x == ln ) {
      runner = ln->b->y;
    } else {
      runner = ln->b->x;
    }

    if( ln->b->twist == 1 ) {
      flipped = !flipped;
    }

  } else {
    if( ln->a->x == ln ) {
      runner = ln->a->y;
    } else {
      runner = ln->a->x;
    }

    /* as b is always on the reverse strand, invert the flipped sense as we read down it*/
    if( ln->a->twist == 0 ) {
      flipped = !flipped;
    }

  }



  prev = ln;

  i++;

  while( runner != NULL ) {
    /* to find the outgoing link from here, test neither a nor b
       is not NULL and figure out the right way to go from the
       fact that the link is back to ourselves */

    if( runner->a == NULL || runner->b == NULL ) {
      /* other end of stream */
      break;
    }
    out->seq[i] = first_char_from_dnanumber(runner->number,nmer_size,flipped);
    i++;

    if( runner->a->x == runner && runner->a->y != prev) {
      next = runner->a->y;
      if( runner->a->twist == 1 ) {
	flipped = !flipped;
      }
    } else if ( runner->a->y == runner && runner->a->x != prev ) {
      next = runner->a->x;
      if( runner->a->twist == 1 ) {
	flipped = !flipped;
      }
    } else if( runner->b->x == runner && runner->b->y != prev) {
      next = runner->b->y;
      if( runner->b->twist == 1 ) {
	flipped = !flipped;
      }
    } else if ( runner->b->y == runner && runner->b->x != prev ) {
      next = runner->b->x;
      if( runner->b->twist == 1 ) {
	flipped = !flipped;
      }
    } else {
      fatal("Unable to move off edge!");
    }

    prev = runner;
    runner = next;
  }

  runner->have_seen = 1;
  return out;

}



/* Function:  write_existing_LinkStream(lna,lnad,dns)
 *
 * Descrip:    Writes in read to existing streams
 *
 *
 * Arg:         lna [UNKN ] Undocumented argument [LinkNumberArray *]
 * Arg:        lnad [UNKN ] Undocumented argument [LinkNumberArrayDebug *]
 * Arg:         dns [UNKN ] Undocumented argument [DnaNumberSequence *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 215 "linkindex.dy"
boolean write_existing_LinkStream(LinkNumberArray * lna,LinkNumberArrayDebug * lnad,DnaNumberSequence * dns)
{
  int i;
  int j;
  LinkStream * stream;
  LinkStream * prev = NULL;
  LinkEdge * edge = NULL;
  int is_new = 0;
  int is_old = 0;

  if( lnad->add_stream != 0 ) {
    fprintf(lnad->ofp,"Adding %s as an existing stream\n",dns->orig->name);
  }
  
  fprintf(stderr,"Adding... on existing...\n");

  for(i=0;i<dns->len;i++) {
    if( lna->array[dns->seq[i].number] == NULL ) {
      lna->array[dns->seq[i].number] = LinkNumber_alloc_std();
      /* now add */
      stream = LinkStream_alloc();
      stream->number = dns->seq[i].number;
      stream->depth = 1;
      stream->have_seen = 0;

      stream->starting_flip = dns->seq[i].flipped;
      add_LinkNumber( lna->array[dns->seq[i].number],stream );
      if( lnad->add_stream > 2 ) {
	fprintf(lnad->ofp,"Adding new position %s, on %d [%d,%c]\n",dns->orig->name,dns->seq[i].number,i,dns->orig->seq[i]);
      }
      is_new = 1;
      if( prev != NULL ) {
	edge = new_LinkEdge(lna);
	edge->x = prev;
	edge->y = stream;

	if( dns->seq[i].flipped == dns->seq[i-1].flipped ) {
	  edge->twist = 0;
	} else {
	  edge->twist = 1;
	}
	

	if( is_old != 1 ) {
	  /* set up arbitary situation */
	  prev->b = edge;
	  stream->a = edge;
	  stream->b = NULL; /* terminating */
	} else {
	  if( prev->a == NULL ) {
	    prev->a = edge;
	  } else {
	    prev->b = edge;
	  }
	  stream->a = edge;
	  stream->b = NULL; /* terminating */
	  is_old = 0;
	}
      } else {
	/* first position, assign a to NULL */
	stream->a = NULL;
	stream->b = NULL;
      }

      fprintf(stderr,"Leaving with position %d with %d on a and %d on b\n",dns->seq[i].number,
	      lna->array[dns->seq[i].number]->stream[0]->a,
	      lna->array[dns->seq[i].number]->stream[0]->b);

	      
      prev = stream;

      continue;
    } else {
      /* as this number has a stream, and we are adding, one of the streams have to be right*/

      /* if this is first, then we know this must be unambiguous */
      if( i == 0 ) {
	if( lna->array[dns->seq[0].number]->len > 1 ) {
	  fatal("Cannot deal with ambiguous first arrays");
	}
	lna->array[dns->seq[0].number]->stream[0]->depth++;
	prev = lna->array[dns->seq[0].number]->stream[0];
	continue;
      }
      
      /* otherwise, one of these streams should make sense, if not die! */

      for(j=0;j<lna->array[dns->seq[i].number]->len;j++) {
	auto LinkStream * stream = lna->array[dns->seq[i].number]->stream[j];
	if( lnad->add_stream != 0 ) {
	  fprintf(lnad->ofp,"Adding %s as an existing stream, testing %d [%d]\n",dns->orig->name,j,lna->array[dns->seq[i].number]);
	}

	if( stream->a == NULL || stream->b == NULL) {
	  /* joining to existing stream from adding things in? Or leaving stream? */
	  stream->depth++;
	  if( is_new == 1 ) {
	    /* we are entering the existing stream here */
	    edge = new_LinkEdge(lna);
	    edge->x = prev;
	    edge->y = stream;
	    if( dns->seq[i].flipped == dns->seq[i-1].flipped ) {
	      edge->twist = 0;
	    } else {
	      edge->twist = 1;
	    }
	    /* assumptions that we make as we are entering here */
	    assert(prev);
	    assert(prev->a != NULL);
	    assert(prev->b == NULL);
	      
	    prev->b = edge;
	    if( stream->a == NULL ) {
	      stream->a = edge;
	    } else {
	      stream->b = edge;
	    }
	    is_new = 0;
	    prev = stream;
	  } else {
	    /* we are leaving the existing stream here
	       have to set up the situtation so that the next link
	       works ok 
	    */
	    is_old = 1;
	    prev = stream;
	  }
	  break; /* out of over all streams */
	}

	

	if( is_linked_LinkStream(prev,lna->array[dns->seq[i].number]->stream[j],dns->seq[i-1].flipped == dns->seq[i].flipped? 0 : 1) != NULL ) {
	  lna->array[dns->seq[i].number]->stream[j]->depth++;
	  prev = lna->array[dns->seq[i].number]->stream[j];
	  break;
	}
      }
      /* this means we have no stream position that made sense */
      if( j >= lna->array[dns->seq[i].number]->len ) {
	fatal("We have a problem. No valid stream, despite adding to existing stream");
      }

      
    } /* end of else is filled */

  } /* end of final position */


}

/* Function:  write_new_LinkStream(lna,lnad,dns)
 *
 * Descrip:    Writes in read as a new stream
 *
 *
 * Arg:         lna [UNKN ] Undocumented argument [LinkNumberArray *]
 * Arg:        lnad [UNKN ] Undocumented argument [LinkNumberArrayDebug *]
 * Arg:         dns [UNKN ] Undocumented argument [DnaNumberSequence *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 369 "linkindex.dy"
boolean write_new_LinkStream(LinkNumberArray * lna,LinkNumberArrayDebug * lnad,DnaNumberSequence * dns)
{
  LinkStream * stream;
  LinkStream * prev;
  LinkEdge * edge;
  
  int i;

  /* First position is odd */

  if( lna->array[dns->seq[0].number] == NULL ) {
    lna->array[dns->seq[0].number] = LinkNumber_alloc_std();
  }
  stream = LinkStream_alloc();
  stream->number = dns->seq[0].number;
  stream->starting_flip = dns->seq[0].flipped;
  stream->depth = 1;
  add_LinkNumber( lna->array[dns->seq[0].number],stream );
  if( lnad->add_stream > 2 ) {
    fprintf(lnad->ofp,"Adding %s, on %d\n",dns->orig->name,dns->seq[0].number);
  }
  
  stream->a = NULL;
  prev = stream;


  for(i=1;i<dns->len;i++) {
    if( lna->array[dns->seq[i].number] == NULL ) {
      lna->array[dns->seq[i].number] = LinkNumber_alloc_std();
    }
    stream = LinkStream_alloc();
    stream->number = dns->seq[i].number;
    stream->starting_flip = dns->seq[i].flipped;
    stream->depth = 1;
    add_LinkNumber( lna->array[dns->seq[i].number],stream );

    if( lnad->add_stream > 2 ) {
      fprintf(lnad->ofp,"Adding %s, on %d\n",dns->orig->name,dns->seq[i].number);
    }
    
    edge = new_LinkEdge(lna);
    edge->x = prev;
    edge->y = stream;
    if( dns->seq[i].flipped == dns->seq[i-1].flipped ) {
      edge->twist = 0;
    } else {
      edge->twist = 1;
    }

    prev->b = edge;
    stream->a = edge;

    prev = stream;
  }

  /* final guy, finish off */
  prev->b = NULL;
    

}

/* Function:  is_linked_LinkStream(l_one,l_two,twist)
 *
 * Descrip:    Returns the right edge which potentially links these two link streams
 *
 *
 * Arg:        l_one [UNKN ] Undocumented argument [LinkStream *]
 * Arg:        l_two [UNKN ] Undocumented argument [LinkStream *]
 * Arg:        twist [UNKN ] Undocumented argument [char]
 *
 * Return [UNKN ]  Undocumented return value [LinkEdge *]
 *
 */
# line 433 "linkindex.dy"
LinkEdge * is_linked_LinkStream(LinkStream * l_one,LinkStream * l_two,char twist)
{
  if( l_one->a != NULL && l_one->b != NULL ) {
    fprintf(stderr,"Testing Link %d with Edges %d,%d and links %d,%d [%d] and %d,%d [%d] to %d [%d]\n",l_one,l_one->a,l_one->b,l_one->a->x,l_one->a->y,l_one->a->twist,l_one->b->x,l_one->b->y,l_one->b->twist,l_two,twist);
  }
  
  if( l_one->a != NULL && l_one->a->x == l_one && l_one->a->y == l_two && l_one->a->twist == twist ) {
    return l_one->a;
  }

  if( l_one->a != NULL && l_one->a->y == l_one && l_one->a->x == l_two && l_one->a->twist == twist ) {
    return l_one->a;
  }

  if( l_one->b != NULL && l_one->b->x == l_one && l_one->b->y == l_two && l_one->b->twist == twist ) {
    return l_one->b;
  }

  if( l_one->b != NULL && l_one->b->y == l_one && l_one->b->x == l_two && l_one->b->twist == twist ) {
    return l_one->b;
  }

  fprintf(stderr,"Going to return NULL\n");

  return NULL;
}

/* Function:  is_new_stream_LinkNumberArray(lna,lnad,dns)
 *
 * Descrip:    Determines whether this DnaNumberSequence should be newly
 *             streamed or not
 *
 *
 * Arg:         lna [UNKN ] Undocumented argument [LinkNumberArray *]
 * Arg:        lnad [UNKN ] Undocumented argument [LinkNumberArrayDebug *]
 * Arg:         dns [UNKN ] Undocumented argument [DnaNumberSequence *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 464 "linkindex.dy"
boolean is_new_stream_LinkNumberArray(LinkNumberArray * lna,LinkNumberArrayDebug * lnad,DnaNumberSequence * dns)
{
  int i,j,k;
  boolean seen_end = 0;
  LinkStream * prev;

  assert(lna);
  assert(dns);

  /* find starting position */
  for(i=0;i<dns->len;i++) {
    if( lnad->placement_stream > 2 ) {
      fprintf(lnad->ofp,"Read %s is unique testing position %d with %d [%d]\n",dns->orig->name,i,dns->seq[i].number,lna->array[dns->seq[i].number]);
    }
    if( lna->array[dns->seq[i].number] != NULL ) 
      break;
  }
   
  if( i >= dns->len ) {
    /* virgin - return TRUE */
    if( lnad->placement_stream != 0 ) {
      fprintf(lnad->ofp,"Read %s is unique, starting new stream\n",dns->orig->name);
    }
    return TRUE;
  }

  if( lna->array[dns->seq[i].number]->len > 1 ) {
    if( lnad->placement_stream != 0 ) {
      fprintf(lnad->ofp,"Read %s starts on ambiguous position, can't cope\n",dns->orig->name);
    }
    return TRUE;
  }

  prev = lna->array[dns->seq[i].number]->stream[0];

  for(i++;i<dns->len;i++) {
    if( lna->array[dns->seq[i].number] == NULL ) 
      break;

    for(j=0;j<lna->array[dns->seq[i].number]->len;j++) {
      auto LinkStream * stream = lna->array[dns->seq[i].number]->stream[j];
      auto LinkEdge * edge = is_linked_LinkStream(stream,prev,dns->seq[i-1].flipped == dns->seq[i].flipped ? 0 : 1);

      if( edge == NULL ) {
	if( lnad->placement_stream != 0 ) {
	  fprintf(lnad->ofp,"Read %s, edge failure at position %d\n",dns->orig->name,i);
	}
	return TRUE;
      }

      if( lnad->placement_stream > 2 ) {
	fprintf(lnad->ofp,"Read %s, successfully found stream for position %d\n",dns->orig->name,i);
      }
      prev = stream;
    }
  }

  if( i < dns->len ) {
    /* reached end, on a stream! return FALSE */
    if( lnad->placement_stream != 0 ) {
      fprintf(lnad->ofp,"Read %s is streamable to end\n",dns->orig->name);
    }
    return FALSE;
  }

  /* still could be valid if NULL to end of seq */
  for(;i<dns->len;i++) {
    if( lna->array[dns->seq[i].number] == NULL ) {
      if( lnad->placement_stream != 0 ) {
	fprintf(lnad->ofp,"Read %s was streamable to end, but has a re-entry\n",dns->orig->name);
      }
      return TRUE;
    }
  }

  if( lnad->placement_stream != 0 ) {
    fprintf(lnad->ofp,"Read %s is streamable to end, but does have unique tail\n",dns->orig->name);
  }

  return FALSE;
}



/* Function:  new_LinkNumberArray(nmer_size)
 *
 * Descrip:    Returns a LinkNumberArray with the appropiate
 *             array size for this size of Nmer
 *
 *
 * Arg:        nmer_size [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [LinkNumberArray *]
 *
 */
# line 552 "linkindex.dy"
LinkNumberArray * new_LinkNumberArray(int nmer_size)
{
  int size = 1;
  int i;
  LinkNumberArray * out;

  for(i=0;i<nmer_size;i++) {
    size *= 4;
  }

  out = LinkNumberArray_alloc_std();

  out->array = calloc(size,sizeof(LinkNumber*));
  out->nmer_size = nmer_size;
  out->array_len = size;

  return out;

}


/* Function:  new_LinkEdge(lna)
 *
 * Descrip:    makes a new link edge added into the memory management
 *
 *
 * Arg:        lna [UNKN ] Undocumented argument [LinkNumberArray *]
 *
 * Return [UNKN ]  Undocumented return value [LinkEdge *]
 *
 */
# line 576 "linkindex.dy"
LinkEdge * new_LinkEdge(LinkNumberArray * lna)
{
  LinkEdge * new;

  assert(lna);
  new = LinkEdge_alloc();
  add_LinkNumberArray(lna,new);

  return new;

}


# line 593 "linkindex.c"
/* Function:  hard_link_LinkEdge(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [LinkEdge *]
 *
 * Return [UNKN ]  Undocumented return value [LinkEdge *]
 *
 */
LinkEdge * hard_link_LinkEdge(LinkEdge * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a LinkEdge object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  LinkEdge_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LinkEdge *]
 *
 */
LinkEdge * LinkEdge_alloc(void) 
{
    LinkEdge * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(LinkEdge *) ckalloc (sizeof(LinkEdge))) == NULL)    {  
      warn("LinkEdge_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->twist = 'u';    


    return out;  
}    


/* Function:  free_LinkEdge(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [LinkEdge *]
 *
 * Return [UNKN ]  Undocumented return value [LinkEdge *]
 *
 */
LinkEdge * free_LinkEdge(LinkEdge * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a LinkEdge obj. Should be trappable");  
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
    /* obj->x is linked in */ 
    /* obj->y is linked in */ 


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_LinkStream(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [LinkStream *]
 *
 * Return [UNKN ]  Undocumented return value [LinkStream *]
 *
 */
LinkStream * hard_link_LinkStream(LinkStream * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a LinkStream object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  LinkStream_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LinkStream *]
 *
 */
LinkStream * LinkStream_alloc(void) 
{
    LinkStream * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(LinkStream *) ckalloc (sizeof(LinkStream))) == NULL)    {  
      warn("LinkStream_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->number = 0; 
    out->depth = 'u';    
    out->starting_flip = 'u';    
    out->have_seen = 'u';    


    return out;  
}    


/* Function:  free_LinkStream(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [LinkStream *]
 *
 * Return [UNKN ]  Undocumented return value [LinkStream *]
 *
 */
LinkStream * free_LinkStream(LinkStream * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a LinkStream obj. Should be trappable");    
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
    /* obj->a is linked in */ 
    /* obj->b is linked in */ 


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_LinkNumber(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_LinkNumber
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [LinkStream **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_LinkNumber(LinkStream ** list,int i,int j)  
{
    LinkStream * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_LinkNumber(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_LinkNumber which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [LinkStream **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_LinkNumber(LinkStream ** list,int left,int right,int (*comp)(LinkStream * ,LinkStream * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_LinkNumber(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_LinkNumber (list,++last,i); 
      }  
    swap_LinkNumber (list,left,last);    
    qsort_LinkNumber(list,left,last-1,comp); 
    qsort_LinkNumber(list,last+1,right,comp);    
}    


/* Function:  sort_LinkNumber(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_LinkNumber
 *
 *
 * Arg:         obj [UNKN ] Object containing list [LinkNumber *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_LinkNumber(LinkNumber * obj,int (*comp)(LinkStream *, LinkStream *)) 
{
    qsort_LinkNumber(obj->stream,0,obj->len-1,comp); 
    return;  
}    


/* Function:  expand_LinkNumber(obj,len)
 *
 * Descrip:    Really an internal function for add_LinkNumber
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [LinkNumber *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_LinkNumber(LinkNumber * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_LinkNumber called with no need"); 
      return TRUE;   
      }  


    if( (obj->stream = (LinkStream ** ) ckrealloc (obj->stream,sizeof(LinkStream *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_LinkNumber, returning FALSE");   
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_LinkNumber(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [LinkNumber *]
 * Arg:        add [OWNER] Object to add to the list [LinkStream *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_LinkNumber(LinkNumber * obj,LinkStream * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_LinkNumber(obj,obj->len + LinkNumberLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->stream[obj->len++]=add; 
    return TRUE; 
}    


/* Function:  flush_LinkNumber(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [LinkNumber *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_LinkNumber(LinkNumber * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->stream[i] != NULL)    {  
        free_LinkStream(obj->stream[i]); 
        obj->stream[i] = NULL;   
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  LinkNumber_alloc_std(void)
 *
 * Descrip:    Equivalent to LinkNumber_alloc_len(LinkNumberLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LinkNumber *]
 *
 */
LinkNumber * LinkNumber_alloc_std(void) 
{
    return LinkNumber_alloc_len(LinkNumberLISTLENGTH);   
}    


/* Function:  LinkNumber_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [LinkNumber *]
 *
 */
LinkNumber * LinkNumber_alloc_len(int len) 
{
    LinkNumber * out;   /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = LinkNumber_alloc()) == NULL)   
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->stream = (LinkStream ** ) ckcalloc (len,sizeof(LinkStream *))) == NULL) {  
      warn("Warning, ckcalloc failed in LinkNumber_alloc_len");  
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_LinkNumber(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [LinkNumber *]
 *
 * Return [UNKN ]  Undocumented return value [LinkNumber *]
 *
 */
LinkNumber * hard_link_LinkNumber(LinkNumber * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a LinkNumber object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  LinkNumber_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LinkNumber *]
 *
 */
LinkNumber * LinkNumber_alloc(void) 
{
    LinkNumber * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(LinkNumber *) ckalloc (sizeof(LinkNumber))) == NULL)    {  
      warn("LinkNumber_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->stream = NULL;  
    out->len = out->maxlen = 0;  
    out->number = 0; 


    return out;  
}    


/* Function:  free_LinkNumber(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [LinkNumber *]
 *
 * Return [UNKN ]  Undocumented return value [LinkNumber *]
 *
 */
LinkNumber * free_LinkNumber(LinkNumber * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a LinkNumber obj. Should be trappable");    
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
    if( obj->stream != NULL) {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->stream[i] != NULL)  
          free_LinkStream(obj->stream[i]);   
        }  
      ckfree(obj->stream);   
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_LinkNumberArray(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_LinkNumberArray
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [LinkEdge **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_LinkNumberArray(LinkEdge ** list,int i,int j)  
{
    LinkEdge * temp; 
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_LinkNumberArray(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_LinkNumberArray which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [LinkEdge **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_LinkNumberArray(LinkEdge ** list,int left,int right,int (*comp)(LinkEdge * ,LinkEdge * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_LinkNumberArray(list,left,(left+right)/2);  
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_LinkNumberArray (list,++last,i);    
      }  
    swap_LinkNumberArray (list,left,last);   
    qsort_LinkNumberArray(list,left,last-1,comp);    
    qsort_LinkNumberArray(list,last+1,right,comp);   
}    


/* Function:  sort_LinkNumberArray(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_LinkNumberArray
 *
 *
 * Arg:         obj [UNKN ] Object containing list [LinkNumberArray *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_LinkNumberArray(LinkNumberArray * obj,int (*comp)(LinkEdge *, LinkEdge *)) 
{
    qsort_LinkNumberArray(obj->edge_set,0,obj->len-1,comp);  
    return;  
}    


/* Function:  expand_LinkNumberArray(obj,len)
 *
 * Descrip:    Really an internal function for add_LinkNumberArray
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [LinkNumberArray *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_LinkNumberArray(LinkNumberArray * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_LinkNumberArray called with no need");    
      return TRUE;   
      }  


    if( (obj->edge_set = (LinkEdge ** ) ckrealloc (obj->edge_set,sizeof(LinkEdge *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_LinkNumberArray, returning FALSE");  
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_LinkNumberArray(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [LinkNumberArray *]
 * Arg:        add [OWNER] Object to add to the list [LinkEdge *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_LinkNumberArray(LinkNumberArray * obj,LinkEdge * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_LinkNumberArray(obj,obj->len + LinkNumberArrayLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->edge_set[obj->len++]=add;   
    return TRUE; 
}    


/* Function:  flush_LinkNumberArray(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [LinkNumberArray *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_LinkNumberArray(LinkNumberArray * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->edge_set[i] != NULL)  {  
        free_LinkEdge(obj->edge_set[i]); 
        obj->edge_set[i] = NULL; 
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  LinkNumberArray_alloc_std(void)
 *
 * Descrip:    Equivalent to LinkNumberArray_alloc_len(LinkNumberArrayLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LinkNumberArray *]
 *
 */
LinkNumberArray * LinkNumberArray_alloc_std(void) 
{
    return LinkNumberArray_alloc_len(LinkNumberArrayLISTLENGTH); 
}    


/* Function:  LinkNumberArray_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [LinkNumberArray *]
 *
 */
LinkNumberArray * LinkNumberArray_alloc_len(int len) 
{
    LinkNumberArray * out;  /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = LinkNumberArray_alloc()) == NULL)  
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->edge_set = (LinkEdge ** ) ckcalloc (len,sizeof(LinkEdge *))) == NULL)   {  
      warn("Warning, ckcalloc failed in LinkNumberArray_alloc_len"); 
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_LinkNumberArray(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [LinkNumberArray *]
 *
 * Return [UNKN ]  Undocumented return value [LinkNumberArray *]
 *
 */
LinkNumberArray * hard_link_LinkNumberArray(LinkNumberArray * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a LinkNumberArray object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  LinkNumberArray_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LinkNumberArray *]
 *
 */
LinkNumberArray * LinkNumberArray_alloc(void) 
{
    LinkNumberArray * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(LinkNumberArray *) ckalloc (sizeof(LinkNumberArray))) == NULL)  {  
      warn("LinkNumberArray_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->array = NULL;   
    out->array_len = 0;  
    out->nmer_size = 0;  
    out->edge_set = NULL;    
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_LinkNumberArray(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [LinkNumberArray *]
 *
 * Return [UNKN ]  Undocumented return value [LinkNumberArray *]
 *
 */
LinkNumberArray * free_LinkNumberArray(LinkNumberArray * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a LinkNumberArray obj. Should be trappable");   
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
    if( obj->array != NULL)  
      ckfree(obj->array);    
    if( obj->edge_set != NULL)   {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->edge_set[i] != NULL)    
          free_LinkEdge(obj->edge_set[i]);   
        }  
      ckfree(obj->edge_set); 
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_LinkNumberArrayDebug(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [LinkNumberArrayDebug *]
 *
 * Return [UNKN ]  Undocumented return value [LinkNumberArrayDebug *]
 *
 */
LinkNumberArrayDebug * hard_link_LinkNumberArrayDebug(LinkNumberArrayDebug * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a LinkNumberArrayDebug object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  LinkNumberArrayDebug_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LinkNumberArrayDebug *]
 *
 */
LinkNumberArrayDebug * LinkNumberArrayDebug_alloc(void) 
{
    LinkNumberArrayDebug * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(LinkNumberArrayDebug *) ckalloc (sizeof(LinkNumberArrayDebug))) == NULL)    {  
      warn("LinkNumberArrayDebug_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->placement_stream = 'u'; 
    out->add_stream = 'u';   
    out->extraction = 'u';   


    return out;  
}    


/* Function:  free_LinkNumberArrayDebug(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [LinkNumberArrayDebug *]
 *
 * Return [UNKN ]  Undocumented return value [LinkNumberArrayDebug *]
 *
 */
LinkNumberArrayDebug * free_LinkNumberArrayDebug(LinkNumberArrayDebug * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a LinkNumberArrayDebug obj. Should be trappable");  
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
    /* obj->ofp is linked in */ 


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
