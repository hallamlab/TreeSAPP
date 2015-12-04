#ifdef _cplusplus
extern "C" {
#endif
#include "arrayseqlookup.h"


  typedef struct ArraySeqClient_struct {
    ArraySeqHeadResults      result;
    SeqLookupResultInterface interface;
    ArraySeqLookup * array;
  } ArraySeqClient;

 static ArraySeqClient * ArraySeqClient_alloc(void)
 {
  ArraySeqClient * out;

  out = malloc(sizeof(ArraySeqClient));

  return out;
 }


/* Function:  new_ArraySeq_SeqLookupInterface(len,numb_level)
 *
 * Descrip:    Exported function - makes a new seqlookupinterface in array mode
 *
 *
 * Arg:               len [UNKN ] Undocumented argument [int]
 * Arg:        numb_level [UNKN ] Undocumented argument [long]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupInterface *]
 *
 */
# line 49 "arrayseqlookup.dy"
SeqLookupInterface * new_ArraySeq_SeqLookupInterface(int len,long numb_level)
{
  SeqLookupInterface * out;
  
  out = SeqLookupInterface_alloc_std();

  out->data = (void*) new_ArraySeqLookup(len,numb_level);

  out->get_client    = get_client_arraylookup;
  out->add_seq       = add_seq_arraylookup;
  out->free_data     = free_data_arraylookup;
  out->lookup_array_head = arrayhead_direct_lookup;

  return out;
}



/* Function:  print_array_occuypancy_ArraySeq(asl,ofp)
 *
 * Descrip:    Prints out summary statistcis to a file
 *
 *
 * Arg:        asl [UNKN ] Undocumented argument [ArraySeqLookup *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 70 "arrayseqlookup.dy"
void print_array_occuypancy_ArraySeq(ArraySeqLookup * asl,FILE * ofp)
{
  int i;

  for(i=0;i<asl->len;i++) {
    if( asl->array[i] == NULL ) {
      fprintf(ofp,"%d EMPTY\n",i);
    } else {
      fprintf(ofp,"%d full %d\n",i,asl->array[i]->current_pos);
    }
  }

}

/* Function:  get_client_arraylookup(*data)
 *
 * Descrip:    Builds a new client from Array
 *
 *
 * Arg:        *data [UNKN ] Undocumented argument [void]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupClientInterface *]
 *
 */
# line 87 "arrayseqlookup.dy"
SeqLookupClientInterface * get_client_arraylookup(void *data)
{
  ArraySeqLookup * array = (ArraySeqLookup *) data;
  ArraySeqClient * cli;
  SeqLookupClientInterface * slci;
 
  cli = ArraySeqClient_alloc();

  cli->interface.next         =  next_arrayhead_search_results;
  cli->interface.is_more      =  is_more_arrayhead_search_results;
  cli->interface.free_data    =  free_arrayhead_results;
   
  cli->interface.data = (void *) &(cli->result);
  cli->result.ipos = 0;
  cli->array = array;


  slci = SeqLookupClientInterface_alloc();
  slci->lookup       =  lookup_array_client;
  slci->is_populated = is_populated_array_client;
  slci->free_data    = free_array_client;
  slci->data = (void*) cli;

  return slci;
}

/* Function:  next_arrayhead_search_results(data,prev)
 *
 * Descrip:    Internal function for results interface
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 * Arg:        prev [UNKN ] Undocumented argument [SeqLookupResultStruct *]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupResultStruct *]
 *
 */
# line 116 "arrayseqlookup.dy"
SeqLookupResultStruct * next_arrayhead_search_results(void * data,SeqLookupResultStruct * prev)
{
  ArraySeqHeadResults * a = (ArraySeqHeadResults *) data;

  if( a->ipos >= a->head->current_pos ) {
    fatal("Overran array!");
  }

  a->str.seq = a->head->units[a->ipos].seq;
  a->str.pos = a->head->units[a->ipos].pos;

  a->ipos++;

  return &(a->str);

}

/* Function:  is_more_arrayhead_search_results(data)
 *
 * Descrip:    Internal function for results interface
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 136 "arrayseqlookup.dy"
boolean is_more_arrayhead_search_results(void * data)
{
  ArraySeqHeadResults * a = (ArraySeqHeadResults *) data;

  if( a->ipos >= a->head->current_pos ) {
    return FALSE;
  }
  return TRUE;
}



/* Function:  free_arrayhead_results(data)
 *
 * Descrip:    Internal function for results interface, which is a no-op
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 *
 */
# line 151 "arrayseqlookup.dy"
void free_arrayhead_results(void * data)
{

  return;
}


/* Function:  free_array_client(data)
 *
 * Descrip:    Internal function for client interface, which frees client specific memory
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 *
 */
# line 161 "arrayseqlookup.dy"
void free_array_client(void * data)
{
  ArraySeqClient * cli = (ArraySeqClient*) data;
  free(cli);
}

  

/* Function:  new_ArraySeqLookup(len,numb_level)
 *
 * Descrip:    makes a new ArraySeqLookup taking up
 *             to len positions
 *
 *
 * Arg:               len [UNKN ] Undocumented argument [int]
 * Arg:        numb_level [UNKN ] Undocumented argument [long]
 *
 * Return [UNKN ]  Undocumented return value [ArraySeqLookup *]
 *
 */
# line 173 "arrayseqlookup.dy"
ArraySeqLookup * new_ArraySeqLookup(int len,long numb_level)
{
  int i;
  ArraySeqLookup * out;

  out = ArraySeqLookup_alloc();
  out->array = calloc(len,sizeof(ArraySeqHead*));
  assert(out->array);
  
  for(i=0;i<len;i++) {
    out->array[i] = NULL;
  }

  out->len = len;
  out->numb_level = numb_level;


  return out;
}


/* Function:  is_populated_array_client(data,seq_number)
 *
 * Descrip:    tells whether this is populated or not
 *
 *
 * Arg:              data [UNKN ] Undocumented argument [void *]
 * Arg:        seq_number [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 197 "arrayseqlookup.dy"
boolean is_populated_array_client(void * data, int seq_number)
{
  ArraySeqClient * cli = (ArraySeqClient *)data;

  if( cli->array->array[seq_number] == NULL ) {
    return FALSE;
  } else {
    return TRUE;
  }

}

/* Function:  lookup_array_client(data,seq_number)
 *
 * Descrip:    Retrieves a SeqLookup position 
 *
 *
 * Arg:              data [UNKN ] Undocumented argument [void *]
 * Arg:        seq_number [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupResultInterface *]
 *
 */
# line 212 "arrayseqlookup.dy"
SeqLookupResultInterface * lookup_array_client(void * data, int seq_number)
{
  ArraySeqClient * cli = (ArraySeqClient *)data;

  if( cli->array->array[seq_number] == NULL ) {
    return NULL;
  }

  cli->result.ipos = 0;
  cli->result.head = cli->array->array[seq_number];

  return &(cli->interface);
}

/* Function:  arrayhead_direct_lookup(data,seq_number)
 *
 * Descrip:    For array optimised lookup hash, provides direct memory access
 *
 *
 * Arg:              data [UNKN ] Undocumented argument [void *]
 * Arg:        seq_number [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [ArraySeqHead *]
 *
 */
# line 229 "arrayseqlookup.dy"
ArraySeqHead * arrayhead_direct_lookup(void * data,int seq_number)
{
  ArraySeqLookup * look = (ArraySeqLookup *)data;

/*  fprintf(stderr,"In arrayhead direct lookup with %d\n",seq_number);*/

  return look->array[seq_number];
}


/* Function:  add_seq_arraylookup(data,seq,para)
 *
 * Descrip:    Adds a sequence/pos pair to the hash
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 * Arg:         seq [UNKN ] Undocumented argument [Sequence *]
 * Arg:        para [UNKN ] Undocumented argument [SeqLookupLoadPara *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 242 "arrayseqlookup.dy"
boolean add_seq_arraylookup(void * data,Sequence * seq,SeqLookupLoadPara * para)
{
  int i;
  ArraySeqLookup * look = (ArraySeqLookup *)data;
  int seq_number;

  assert(data != NULL);
  assert(seq != NULL);
  assert(para != NULL);

  for(i=0;i<seq->len-5;i = i+para->tile_freq) {
    seq_number = seq_number_aa_5mer(seq->seq+i);

    if( look->array[seq_number] == NULL ) {
      look->array[seq_number] =  new_ArraySeqHead();
      if( para->mark_low_complexity ) {
	look->array[seq_number]->flags = flags_from_5aa_sequence(seq->seq+i);
      }
    }

    add_ArraySeqHead(look->array[seq_number],seq,i,look->numb_level);
  }

  return TRUE;
}


/* Function:  free_data_arraylookup(data)
 *
 * Descrip:    Frees data
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 *
 */
# line 272 "arrayseqlookup.dy"
void free_data_arraylookup(void * data)
{
  ArraySeqLookup * look = (ArraySeqLookup *)data;
  int i;

  for(i=0;i<look->len;i++) {
    if( look->array[i] == NULL ) {
      continue;
    }
    ckfree(look->array[i]->units);
    ckfree(look->array[i]);
  }
  
  free_ArraySeqLookup(look);
}


/* Function:  add_ArraySeqHead(h,seq,pos,numb_level)
 *
 * Descrip:    Adds a sequence/pos pair to an ArrayHead
 *
 *
 * Arg:                 h [UNKN ] Undocumented argument [ArraySeqHead *]
 * Arg:               seq [UNKN ] Undocumented argument [Sequence *]
 * Arg:               pos [UNKN ] Undocumented argument [int]
 * Arg:        numb_level [UNKN ] Undocumented argument [long]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 292 "arrayseqlookup.dy"
boolean add_ArraySeqHead(ArraySeqHead * h,Sequence * seq,int pos,long numb_level)
{
  ArraySeqLookupUnit * temp;

  if( numb_level  > 0 && h->current_pos > numb_level ) {
    return TRUE;
  }

  /*  fprintf(stderr,"adding new sequence position with %d max pos %d\n",h->current_pos,h->max);
   */
  if( h->current_pos >= h->max ) {
    temp = h->units;
    if( h->max < ARRAYSEQL_LINEAR ) {
      h->units = realloc(h->units,(h->max*2)*sizeof(ArraySeqLookupUnit));
      h->max = h->max*2;
    } else {
      h->units = realloc(h->units,(h->max + ARRAYSEQL_LINEAR)*sizeof(ArraySeqLookupUnit));
      h->max = h->max + ARRAYSEQL_LINEAR;
    }

    if( h->units == NULL ) {
      fatal("ArraySeqLookup realloc failed trying for %d positions\n",h->max);
    }

    
/*    fprintf(stderr,"... extended to %d\n",h->max); */
  }

  h->units[h->current_pos].seq = seq;
  h->units[h->current_pos].pos = pos;

  h->current_pos++;

  return TRUE;
}


/* Function:  new_ArraySeqHead(void)
 *
 * Descrip:    Builds a new ArraySeqHead structure
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ArraySeqHead *]
 *
 */
# line 332 "arrayseqlookup.dy"
ArraySeqHead * new_ArraySeqHead(void)
{
  ArraySeqHead * out;

  out = malloc(sizeof(ArraySeqHead));

  out->units = calloc(ARRAYSEQL_BASIC,sizeof(ArraySeqLookupUnit));
  out->max = ARRAYSEQL_BASIC;
  out->current_pos = 0;
  
  return out;
}
	       


# line 416 "arrayseqlookup.c"
/* Function:  hard_link_ArraySeqLookup(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ArraySeqLookup *]
 *
 * Return [UNKN ]  Undocumented return value [ArraySeqLookup *]
 *
 */
ArraySeqLookup * hard_link_ArraySeqLookup(ArraySeqLookup * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a ArraySeqLookup object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  ArraySeqLookup_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ArraySeqLookup *]
 *
 */
ArraySeqLookup * ArraySeqLookup_alloc(void) 
{
    ArraySeqLookup * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(ArraySeqLookup *) ckalloc (sizeof(ArraySeqLookup))) == NULL)    {  
      warn("ArraySeqLookup_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->array = NULL;   
    out->len = 0;    
    out->numb_level = 0; 


    return out;  
}    


/* Function:  free_ArraySeqLookup(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ArraySeqLookup *]
 *
 * Return [UNKN ]  Undocumented return value [ArraySeqLookup *]
 *
 */
ArraySeqLookup * free_ArraySeqLookup(ArraySeqLookup * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a ArraySeqLookup obj. Should be trappable");    
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


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_ArraySeqHeadResults(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ArraySeqHeadResults *]
 *
 * Return [UNKN ]  Undocumented return value [ArraySeqHeadResults *]
 *
 */
ArraySeqHeadResults * hard_link_ArraySeqHeadResults(ArraySeqHeadResults * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a ArraySeqHeadResults object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  ArraySeqHeadResults_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ArraySeqHeadResults *]
 *
 */
ArraySeqHeadResults * ArraySeqHeadResults_alloc(void) 
{
    ArraySeqHeadResults * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(ArraySeqHeadResults *) ckalloc (sizeof(ArraySeqHeadResults))) == NULL)  {  
      warn("ArraySeqHeadResults_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->ipos = 0;   


    return out;  
}    


/* Function:  free_ArraySeqHeadResults(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ArraySeqHeadResults *]
 *
 * Return [UNKN ]  Undocumented return value [ArraySeqHeadResults *]
 *
 */
ArraySeqHeadResults * free_ArraySeqHeadResults(ArraySeqHeadResults * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a ArraySeqHeadResults obj. Should be trappable");   
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
    /* obj->head is linked in */ 
    /* obj->str is linked in */ 


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
