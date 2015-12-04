#ifdef _cplusplus
extern "C" {
#endif
#include "compressed_protein_index.h"


# line 39 "compressed_protein_index.dy"
SeqLookupResultStruct * get_next_CompressedProteinClient(void * data,SeqLookupResultStruct * prev)
{
  CompressedProteinClient * cl = (CompressedProteinClient *)data;

/*  fprintf(stdout," entering with point %d max %d\n",cl->point,cl->current); */

  if( prev == NULL ) {
    assert(cl->point == 0);
  } else {
    assert(&(cl->res[cl->point-1]) == prev);
  }


  if( cl->point == cl->current ) {
    return NULL;
  } else {
    return &cl->res[cl->point++];
  }

}


# line 61 "compressed_protein_index.dy"
boolean is_more_CompressedProteinClient(void * data)
{
  CompressedProteinClient * cl = (CompressedProteinClient *)data;

  if( cl->point < cl->current ) {
    return TRUE;
  } else {
    return FALSE;
  }

}

# line 73 "compressed_protein_index.dy"
void free_SeqLookupResultInterface_CompressedProteinClient(void * data)
{
  /* no op*/
  return;
}


# line 80 "compressed_protein_index.dy"
SeqLookupResultInterface * lookup_CompressedProteinClient(void * data,int seq_number)
{
  CompressedProteinClient * cl = (CompressedProteinClient *)data;
  int i;
  CompressedProteinHead * h;
  SeqLookupResultInterface * out;
  SinglePosSequence * pos;
  int len;

  out = SeqLookupResultInterface_alloc();
  h = (CompressedProteinHead*)(*cl->index->kii->retrieve_by_kmer)(cl->index->kii->handle,(kmer_t)seq_number);
  if( h == NULL ) {
    cl->current = 0;
    cl->point = 0;
  } else {
    for(i=0;i<h->maxlen;i++) {
      if( h->positions[i] == -1 ) {
	break;
      }
    }
    assert(i < h->maxlen);
    len = i;
    
    if( cl->maxlen < len ) {
      cl->maxlen = cl->maxlen * 2;
      cl->res = realloc(cl->res,sizeof(SeqLookupResultStruct)*cl->maxlen);
    }

    for(i=0;i<len;i++) {
      pos = lookup_Sequence_SinglePosSpace(cl->index->sps,h->positions[i]);
      cl->res[i].seq = (Sequence*)pos->data;
      cl->res[i].pos = (int)(h->positions[i] - pos->start);
    }

    cl->current = len;
    cl->point = 0;
  }

  out->data = cl;
  out->next      = get_next_CompressedProteinClient;
  out->is_more   = is_more_CompressedProteinClient;
  out->free_data = free_SeqLookupResultInterface_CompressedProteinClient; 
  return out;
}

# line 125 "compressed_protein_index.dy"
boolean is_populated_CompressedProteinClient(void * data,int seq_number)
{
  CompressedProteinHead * h;
  CompressedProteinClient * cl = (CompressedProteinClient *)data;


  h = (CompressedProteinHead*)(*cl->index->kii->retrieve_by_kmer)(cl->index->kii->handle,(kmer_t)seq_number);
    
  if( h == NULL ) {
    return FALSE;
  } else {
    return TRUE;
  }
}

# line 140 "compressed_protein_index.dy"
void free_Client_CompressedProteinClient(void * data)
{
  CompressedProteinClient * cl = (CompressedProteinClient *)data;
  free(cl->res);

  free_CompressedProteinClient(cl);
}


# line 149 "compressed_protein_index.dy"
SeqLookupInterface * new_direct_CompressedProteinLookup(void)
{
  KmerIndexInterface * kii;

  kii = new_KmerDirectIndex(12);

  return new_CompressedProteinLookup(kii);

}

# line 159 "compressed_protein_index.dy"
SeqLookupInterface * new_CompressedProteinLookup(KmerIndexInterface * kii)
{
  CompressedProteinIndex * cpi;
  SeqLookupInterface * out;

  out = SeqLookupInterface_alloc_std();
  cpi = CompressedProteinIndex_alloc();

  cpi->kii = kii;
  cpi->sps = new_SinglePosSpace(1,1000);

  out->data = (void*)cpi;
  out->get_client = get_client_CompressedProteinIndex;
  out->add_seq    = add_seq_CompressedProteinIndex;
  out->lookup_array_head = NULL;
  out->add_direct_number = add_direct_number_CompressedProteinIndex;
  out->free_data = free_data_CompressedProteinIndex;
  

  return out;
}



# line 183 "compressed_protein_index.dy"
SeqLookupClientInterface * get_client_CompressedProteinIndex(void * data)
{
  SeqLookupClientInterface * out;
  CompressedProteinClient * cl;

  cl = new_CompressedProteinClient();

  out = SeqLookupClientInterface_alloc();

  cl->index = (CompressedProteinIndex*)data;

  assert(cl->index != NULL);

  out->data = (void*)cl;
  out->lookup = lookup_CompressedProteinClient;
  out->is_populated = is_populated_CompressedProteinClient;
  out->free_data = free_Client_CompressedProteinClient;

  return out;
}


# line 205 "compressed_protein_index.dy"
CompressedProteinClient * new_CompressedProteinClient(void) 
{
  CompressedProteinClient * out;

  out = CompressedProteinClient_alloc();

  out->res = calloc(COMPRESSED_PROTEIN_CLIENT_INITIAL,sizeof(SeqLookupResultStruct));
  out->current = 0;
  out->maxlen = COMPRESSED_PROTEIN_CLIENT_INITIAL;

  return out;
}


# line 219 "compressed_protein_index.dy"
boolean add_direct_number_CompressedProteinIndex(void * data,int seq_number,Sequence * target, int pos) 
{
  fatal("For compressed protein indexes, impossible to add numbers directly");

  return NULL;
}

# line 226 "compressed_protein_index.dy"
ArraySeqHead * lookup_array_head_CompressedProteinIndex(void * data,int seq_number)
{
  fatal("For compressed protein indexes, arrayseqhead lookup is not feasible");

  return NULL;
}


# line 234 "compressed_protein_index.dy"
boolean add_seq_CompressedProteinIndex(void * data,Sequence * seq,SeqLookupLoadPara * para)
{
  long int offset;
  int no;
  int i;
  void * kmer_data;
  CompressedProteinHead * h;
  CompressedProteinIndex * cpi = (CompressedProteinIndex *) data;

  assert(cpi != NULL);
  assert(seq != NULL);

  offset = add_Sequence_SinglePosSpace(cpi->sps,(long int)seq->len,(void *) seq);

  for(i=0;i<seq->len-5;) {
    no = seq_number_aa_5mer(seq->seq+i);
    if( (kmer_data = (*cpi->kii->retrieve_by_kmer)(cpi->kii->handle,(kmer_t)no)) == NULL ) {
      h = new_CompressedProteinHead();
      (*cpi->kii->insert_by_kmer)(cpi->kii->handle,(kmer_t)no,(void*)h);
    } else {
      h = (CompressedProteinHead*) kmer_data;
    }

    add_pos_CompressedProteinHead(h,(int)offset+i);
    i += para->tile_freq;
  }

  return TRUE;
}

# line 264 "compressed_protein_index.dy"
void free_data_CompressedProteinIndex(void * d)
{
  CompressedProteinIndex * cpi = (CompressedProteinIndex *)d;

  assert(cpi != NULL);

  free_CompressedProteinIndex(cpi);

}


# line 275 "compressed_protein_index.dy"
boolean add_pos_CompressedProteinHead(CompressedProteinHead * c,int pos)
{
  int i;
  assert(c != NULL);

  for(i=0;i<c->maxlen;i++) {
    if( c->positions[i] == -1 ) {
      break;
    }
  }
  
  assert(i < c->maxlen);

  if( i+2 == c->maxlen ) {
    c->maxlen += COMPRESSED_EXPAND;
    c->positions = realloc(c->positions,sizeof(int)*c->maxlen);
  }

  c->positions[i] = pos;
  c->positions[i+1] = -1;

  return TRUE;
}


# line 300 "compressed_protein_index.dy"
CompressedProteinHead * new_CompressedProteinHead(void)
{
  CompressedProteinHead * out;

  out = malloc(sizeof(CompressedProteinHead));

  assert(out != NULL);
  
  out->maxlen = COMPRESSED_INITIAL_SIZE;
  out->positions = malloc(sizeof(int)*out->maxlen);
  out->positions[0] = -1;
  return out;

}

# line 298 "compressed_protein_index.c"
/* Function:  hard_link_CompressedProteinIndex(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [CompressedProteinIndex *]
 *
 * Return [UNKN ]  Undocumented return value [CompressedProteinIndex *]
 *
 */
CompressedProteinIndex * hard_link_CompressedProteinIndex(CompressedProteinIndex * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a CompressedProteinIndex object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  CompressedProteinIndex_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [CompressedProteinIndex *]
 *
 */
CompressedProteinIndex * CompressedProteinIndex_alloc(void) 
{
    CompressedProteinIndex * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(CompressedProteinIndex *) ckalloc (sizeof(CompressedProteinIndex))) == NULL)    {  
      warn("CompressedProteinIndex_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->kii = NULL; 
    out->sps = NULL; 


    return out;  
}    


/* Function:  free_CompressedProteinIndex(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [CompressedProteinIndex *]
 *
 * Return [UNKN ]  Undocumented return value [CompressedProteinIndex *]
 *
 */
CompressedProteinIndex * free_CompressedProteinIndex(CompressedProteinIndex * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a CompressedProteinIndex obj. Should be trappable");    
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
    if( obj->kii != NULL)    
      free_KmerIndexInterface(obj->kii);     
    if( obj->sps != NULL)    
      free_SinglePosSpace(obj->sps);     


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_CompressedProteinClient(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [CompressedProteinClient *]
 *
 * Return [UNKN ]  Undocumented return value [CompressedProteinClient *]
 *
 */
CompressedProteinClient * hard_link_CompressedProteinClient(CompressedProteinClient * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a CompressedProteinClient object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  CompressedProteinClient_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [CompressedProteinClient *]
 *
 */
CompressedProteinClient * CompressedProteinClient_alloc(void) 
{
    CompressedProteinClient * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(CompressedProteinClient *) ckalloc (sizeof(CompressedProteinClient))) == NULL)  {  
      warn("CompressedProteinClient_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->maxlen = 0; 
    out->current = 0;    
    out->point = 0;  


    return out;  
}    


/* Function:  free_CompressedProteinClient(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [CompressedProteinClient *]
 *
 * Return [UNKN ]  Undocumented return value [CompressedProteinClient *]
 *
 */
CompressedProteinClient * free_CompressedProteinClient(CompressedProteinClient * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a CompressedProteinClient obj. Should be trappable");   
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
    /* obj->res is linked in */ 
    /* obj->index is linked in */ 


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
