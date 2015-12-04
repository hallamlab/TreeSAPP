#ifdef _cplusplus
extern "C" {
#endif
#include "proteinstreamedindex.h"

  typedef struct ProteinStreamedClient_struct {
    GenericIndexResult result;
    SeqLookupResultInterface   interface;
    ProteinStreamedIndex * index;
  } ProteinStreamedClient;

 static ProteinStreamedClient * ProteinStreamedClient_alloc(void)
 {
   ProteinStreamedClient * out;

  out = malloc(sizeof(ProteinStreamedClient));

  return out;
 }



/* Function:  new_ProteinStreamedIndex_SeqLookupInterface(waypost)
 *
 * Descrip:    Provides a SeqLookupInterface, the common runtime plug-in for indexers
 *             using a ProteinStreamedIndex 
 *
 *
 * Arg:        waypost [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupInterface *]
 *
 */
# line 59 "proteinstreamedindex.dy"
SeqLookupInterface * new_ProteinStreamedIndex_SeqLookupInterface(int waypost)
{
  SeqLookupInterface * out;
  ProteinStreamedIndex * in;

  in = new_ProteinStreamedIndex(waypost);

  out = SeqLookupInterface_alloc_std();

  out->get_client  = get_client_interface_ProteinStreamedIndex;
  out->add_seq     = add_seq_interface_ProteinStreamedIndex;
  out->free_data   = free_interface_ProteinStreamedIndex;
  out->lookup_array_head = NULL;
  out->data = (void*) in;

  return out;
}

/* Function:  get_client_interface_ProteinStreamedIndex(data)
 *
 * Descrip:    gets client interface
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void*]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupClientInterface *]
 *
 */
# line 80 "proteinstreamedindex.dy"
SeqLookupClientInterface * get_client_interface_ProteinStreamedIndex(void* data)
{
  ProteinStreamedIndex * index = (ProteinStreamedIndex*) data;
  SeqLookupClientInterface * slci;
  ProteinStreamedClient * pcl;

  pcl = ProteinStreamedClient_alloc();

  pcl->result.result = calloc(64,sizeof(SeqLookupResultStruct));
  pcl->result.max_len = 64;
  pcl->result.len = 0;


  pcl->interface.next = next_interface_GenericIndexResult;
  pcl->interface.is_more = is_more_interface_GenericIndexResult;
  pcl->interface.free_data = free_noop_GenericIndexResult;
  pcl->interface.data = (void*) &(pcl->result);

  pcl->index = index;

  slci = SeqLookupClientInterface_alloc();
  slci->lookup       =  lookup_interface_ProteinStreamedClient;
  slci->is_populated = is_populated_interface_ProteinStreamedClient;
  slci->free_data    = free_interface_ProteinStreamedClient;
  slci->data = (void*) pcl;

  return slci;
}

/* Function:  lookup_interface_ProteinStreamedClient(data,seq_number)
 *
 * Descrip:    For lookup interface, provides a result
 *
 *
 * Arg:              data [UNKN ] Undocumented argument [void *]
 * Arg:        seq_number [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupResultInterface *]
 *
 */
# line 112 "proteinstreamedindex.dy"
SeqLookupResultInterface * lookup_interface_ProteinStreamedClient(void * data,int seq_number)
{
  ProteinStreamedClient * cli= (ProteinStreamedClient*) data;
  ProteinStreamedIndex * in = cli->index;

  /* reset pointers */
  cli->result.current_pos = 0;
  cli->result.len = 0;

  if( lookup_ProteinStreamedIndex(in,seq_number,&cli->result) == FALSE ) {
    return NULL;
  }

  return &cli->interface;
}

/* Function:  is_populated_interface_ProteinStreamedClient(data,seq_number)
 *
 * Descrip:    populated function for interface
 *
 *
 * Arg:              data [UNKN ] Undocumented argument [void *]
 * Arg:        seq_number [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 131 "proteinstreamedindex.dy"
boolean is_populated_interface_ProteinStreamedClient(void * data,int seq_number)
{
  ProteinStreamedClient * client = (ProteinStreamedClient*) data;

  if( client->index->index[seq_number] != NULL ) {
    return TRUE;
  } else {
    return FALSE;
  }
}

/* Function:  add_seq_interface_ProteinStreamedIndex(data,seq,para)
 *
 * Descrip:    add function for interface
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 * Arg:         seq [UNKN ] Undocumented argument [Sequence *]
 * Arg:        para [UNKN ] Undocumented argument [SeqLookupLoadPara *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 145 "proteinstreamedindex.dy"
boolean add_seq_interface_ProteinStreamedIndex(void * data,Sequence * seq,SeqLookupLoadPara * para)
{
  ProteinStreamedIndex * in = (ProteinStreamedIndex*) data;

  return add_Sequence_ProteinStreamedIndex(in,seq,para);
}

/* Function:  free_interface_ProteinStreamedIndex(data)
 *
 * Descrip:    for interface, frees index
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 *
 */
# line 155 "proteinstreamedindex.dy"
void free_interface_ProteinStreamedIndex(void * data)
{
  ProteinStreamedIndex * in = (ProteinStreamedIndex*) data;

  free_ProteinStreamedIndex(in);
}


/* Function:  free_interface_ProteinStreamedClient(data)
 *
 * Descrip:    Frees the client data
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 *
 */
# line 166 "proteinstreamedindex.dy"
void free_interface_ProteinStreamedClient(void * data)
{
  ProteinStreamedClient * cli = (ProteinStreamedClient *)data;

  free(cli->result.result);

  free(cli);

}
 

 static int seq_number_aa_5mer_psi(char * seq)
{
  int i;
  int ret = 0;
  int base = 1;
  int no = 0;

  for(i=0;i<5;i++) {
    no = toupper(seq[i])-'A';
    if( no > 26 || no < 0 ) {
      no = 'X'-'A';
    }
    ret += base * no;
    base = base * 26;
  }

  return ret;
}


/* Function:  lookup_ProteinStreamedIndex(in,seqno,out)
 *
 * Descrip:    Traverses index for a particular number to return hits
 *
 *
 * Arg:           in [UNKN ] Undocumented argument [ProteinStreamedIndex *]
 * Arg:        seqno [UNKN ] Undocumented argument [int]
 * Arg:          out [UNKN ] Undocumented argument [GenericIndexResult *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 200 "proteinstreamedindex.dy"
boolean lookup_ProteinStreamedIndex(ProteinStreamedIndex * in,int seqno,GenericIndexResult * out)
{
  
  /* look up this position, return NULL if nothing there */
  if( in->index[seqno] == NULL ) {
    return FALSE;
  }

  /* yes, assign the number of counts in
     walk back up looking for wayposts */
  
  /* start with self */
  
  add_lookup_GenericIndexResult_from_Pos(in->index[seqno],out,0,in->waypost_depth);

  return TRUE;
}


/* Function:  add_lookup_GenericIndexResult_from_Pos(p,psir,backtrace,waylength)
 *
 * Descrip:    Adds additional sequences to index result from a particular position
 *
 *
 * Arg:                p [UNKN ] Undocumented argument [ProteinStreamedIndexPos *]
 * Arg:             psir [UNKN ] Undocumented argument [GenericIndexResult *]
 * Arg:        backtrace [UNKN ] Undocumented argument [int]
 * Arg:        waylength [UNKN ] Undocumented argument [int]
 *
 */
# line 222 "proteinstreamedindex.dy"
void add_lookup_GenericIndexResult_from_Pos(ProteinStreamedIndexPos * p,GenericIndexResult * psir,int backtrace,int waylength)
{
  int i;


  assert(p);
  assert(psir);

  /*  printf("Entering add with backtrace %d count %d\n",backtrace,p->count);*/

  while( p != NULL && backtrace < waylength ) {
    /*    printf("Backtracing.... %d with count %d\n",backtrace,p->count); */

    /* test whether this has any wayposts, add them with backtrace */
    if( p->post_len > 0 ) {
      for(i=0;i<p->post_len;i++) {
	/*	printf("....adding another index... %s, %d\n",p->post[i].seq->name,p->post[i].pos);*/
	add_GenericIndexResult(psir,p->post[i].seq,p->post[i].pos+backtrace);
      }
      if( p->post_len == p->count ) {
	/* then we have the finished...*/
	/*	printf("... finished with %d\n",p->count); */
	break;
      }
    }
    backtrace++;

    if( p->prev != NULL ) {
      /* have to recurse */
      /* first go back up the main stream */
      add_lookup_GenericIndexResult_from_Pos(p->first_prev,psir,backtrace,waylength);
      for(i=0;i<p->prev_len;i++) {
	if( p->prev[i] == NULL ) {
	  break;
	}
	add_lookup_GenericIndexResult_from_Pos(p->prev[i],psir,backtrace,waylength);
      }
      /* can assumme this is ok */

      break;
    } else {
      p = p->first_prev;
    }
  }

  return;
}


/* Function:  add_Sequence_ProteinStreamedIndex(in,seq,para)
 *
 * Descrip:    Adds in a sequence into a ProteinStreamedIndex
 *
 *
 * Arg:          in [UNKN ] Undocumented argument [ProteinStreamedIndex *]
 * Arg:         seq [UNKN ] Undocumented argument [Sequence *]
 * Arg:        para [UNKN ] Undocumented argument [SeqLookupLoadPara *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 274 "proteinstreamedindex.dy"
boolean add_Sequence_ProteinStreamedIndex(ProteinStreamedIndex * in,Sequence * seq,SeqLookupLoadPara * para)
{
  ProteinStreamedIndexPos * prev = NULL;
  int i;
  int j;
  int no;

  assert(in);
  assert(seq);

  if( para->tile_freq > 1 ) {
    fatal("Cannot have a streamed index with a tile_freq greater than 1. Misparameterisation");
  }

  for(i=0;i<seq->len-5;i++) {
    no = seq_number_aa_5mer_psi(seq->seq+i);
    if( in->index[no] == NULL ) {
      in->index[no] = new_ProteinStreamedIndexPos();
      /* because this is a newly allocated index, we can put prev into
	 first_prev now with impunity - the first one will have this as NULL
	 and that is fine */
      in->index[no]->first_prev = prev;
    }
    in->index[no]->count++;

    if( prev != NULL ) {
      /* check that reverse pointer is not present already */
      j = 0;
      if( prev != in->index[no]->first_prev ) {
	for(j=0;j<in->index[no]->prev_len;j++) {
	  if( prev == in->index[no]->prev[j] ) {
	    break;
	  }
	}
	if( j >= in->index[no]->prev_len ) {
	  /* is not in there, add it */
	  /* potentially this could be first_prev if this came from
	     the start of another sequence */
	  if( in->index[no]->first_prev == NULL ) {
	    in->index[no]->first_prev = prev;
	  } else {
	    add_pos_ProteinStreamedIndexPos(in->index[no],prev);
	  }
	}
      }
    }
    prev = in->index[no];

    /* deal with waystations */
    if( i == 0 || i % in->waypost_depth == 0 ) {
      add_waypost_ProteinStreamedIndexPos(in->index[no],seq,i);
    }

  }

  return TRUE;
}

/* Function:  dump_ProteinStreamedIndex(in,ofp)
 *
 * Descrip:    dumps in a silly format the ProteinStreamedIndex; for debugging
 *
 *
 * Arg:         in [UNKN ] Undocumented argument [ProteinStreamedIndex *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 335 "proteinstreamedindex.dy"
void dump_ProteinStreamedIndex(ProteinStreamedIndex * in,FILE * ofp)
{
  int i,j;
  int size = 26 * 26 * 26 * 26 * 26;


  for(i=0;i<size;i++) {
    if( in->index[i] != NULL ) {
      fprintf(ofp,"Position %d occupied with %d previous\n WayStations ",i,in->index[i]->prev_len);
      for(j=0;j<in->index[i]->post_len;j++) {
	if( in->index[i] == NULL ) {
	  break;
	}
	fprintf(ofp," %s [%d],",in->index[i]->post[j].seq->name,
		in->index[i]->post[j].pos);
      }
      fprintf(ofp,"\n");
    }
  }

}

/* Function:  new_ProteinStreamedIndex(waypost_depth)
 *
 * Descrip:    Builds a new ProteinStreamedIndex
 *
 *
 * Arg:        waypost_depth [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [ProteinStreamedIndex *]
 *
 */
# line 360 "proteinstreamedindex.dy"
ProteinStreamedIndex * new_ProteinStreamedIndex(int waypost_depth)
{
  ProteinStreamedIndex * out;
  int size = 26 * 26 * 26 * 26 * 26;

  out = malloc(sizeof(ProteinStreamedIndex));
  assert(out);

  out->index = calloc(size,sizeof(ProteinStreamedIndexPos*));
  out->waypost_depth = waypost_depth;
  return out;
}

/* Function:  free_ProteinStreamedIndex(in)
 *
 * Descrip:    Release protein index
 *
 *
 * Arg:        in [UNKN ] Undocumented argument [ProteinStreamedIndex *]
 *
 */
# line 376 "proteinstreamedindex.dy"
void free_ProteinStreamedIndex(ProteinStreamedIndex * in)
{
  int i;
  int size = 26 * 26 * 26 *26 * 26;
  assert(in);

  for(i=0;i<size;i++) {
    if( in->index[i] != NULL ) {
      free_ProteinStreamedIndexPos(in->index[i]);
    }
  }

  free(in->index);
  free(in);

  return;

}


/* Function:  add_pos_ProteinStreamedIndexPos(pos,prev)
 *
 * Descrip:    Adds a position, advancing array if needed
 *
 *
 * Arg:         pos [UNKN ] Undocumented argument [ProteinStreamedIndexPos *]
 * Arg:        prev [UNKN ] Undocumented argument [ProteinStreamedIndexPos *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 399 "proteinstreamedindex.dy"
boolean add_pos_ProteinStreamedIndexPos(ProteinStreamedIndexPos * pos,ProteinStreamedIndexPos * prev)
{
  int i;

  if( pos->prev == NULL ) {
    pos->prev = calloc(PROTEINSTREAMEDINDEX_PREV_START,sizeof(ProteinStreamedIndexPos*));
    pos->prev[0] = NULL;
    pos->prev_len = PROTEINSTREAMEDINDEX_PREV_START;
  }
  
  if( pos->prev[pos->prev_len-1] != NULL ) {
    /* realloc */
    pos->prev = realloc(pos->prev,(2*pos->prev_len)*sizeof(ProteinStreamedIndexPos*));
    pos->prev_len = 2*pos->prev_len;
  }
  for(i=0;i<pos->prev_len;i++) {
    if( pos->prev[i] == NULL ) {
      pos->prev[i] = prev;
      pos->prev[i+1] = NULL;
      break;
    }
  }

  return TRUE;
}

/* Function:  add_waypost_ProteinStreamedIndexPos(pos,seq,seqpos)
 *
 * Descrip:    Adds a way post, advancing array if needed
 *
 *
 * Arg:           pos [UNKN ] Undocumented argument [ProteinStreamedIndexPos *]
 * Arg:           seq [UNKN ] Undocumented argument [Sequence *]
 * Arg:        seqpos [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 428 "proteinstreamedindex.dy"
boolean add_waypost_ProteinStreamedIndexPos(ProteinStreamedIndexPos * pos,Sequence * seq,int seqpos)
{
  if( pos->post == NULL ) {
      pos->post = calloc(PROTEINSTREAMEDINDEX_WAY_START,sizeof(ProteinStreamedIndexWayPost));

      pos->post_len = 0;
      pos->post_maxlen = PROTEINSTREAMEDINDEX_WAY_START;
  }

  if( pos->post_len >= pos->post_maxlen ) {
    pos->post = realloc(pos->post,(2*pos->post_maxlen)*sizeof(ProteinStreamedIndexWayPost));
    pos->post_maxlen = 2*pos->post_maxlen;
  }

  pos->post[pos->post_len].seq = seq;
  pos->post[pos->post_len].pos = seqpos;

  pos->post_len++;
  return TRUE;
}

/* Function:  new_ProteinStreamedIndexPos(void)
 *
 * Descrip:    Creates a new ProteinStreamedIndexPos
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ProteinStreamedIndexPos *]
 *
 */
# line 452 "proteinstreamedindex.dy"
ProteinStreamedIndexPos * new_ProteinStreamedIndexPos(void)
{
  ProteinStreamedIndexPos * out;

  out = malloc(sizeof(ProteinStreamedIndexPos));
  out->count = 0;
  out->prev     = NULL;
  out->prev_len = 0;
  out->post = NULL;
  out->post_len = 0;
  out->post_maxlen = 0;
  return out;
}

/* Function:  free_ProteinStreamedIndexPos(pos)
 *
 * Descrip:    frees a ProteinStreamedIndexPos
 *
 *
 * Arg:        pos [UNKN ] Undocumented argument [ProteinStreamedIndexPos *]
 *
 */
# line 469 "proteinstreamedindex.dy"
void free_ProteinStreamedIndexPos(ProteinStreamedIndexPos * pos)
{
  if( pos == NULL ) {
    return;
  }

  free(pos->prev);
  free(pos->post);

  free(pos);

  return;

}

# line 555 "proteinstreamedindex.c"

#ifdef _cplusplus
}
#endif
