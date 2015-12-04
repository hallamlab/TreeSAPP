#ifdef _cplusplus
extern "C" {
#endif
#include "shadowseqindex.h"


  typedef struct ShadowSeqClient_struct {
    GenericIndexResult result;
    SeqLookupResultInterface   interface;
    ShadowSequenceIndex * index;
  } ShadowSequenceClient;

 static ShadowSequenceClient * ShadowSequenceClient_alloc(void)
 {
   ShadowSequenceClient * out;

  out = malloc(sizeof(ShadowSequenceClient));

  return out;
 }



/* Function:  new_ShadowSequenceIndex_SeqLookupInterface(shadow_len,has_maxlen,maxlen,shadow_error)
 *
 * Descrip:    Provides a SeqLookupInterface, the common runtime plug-in for indexers
 *             using a ShadowSequenceIndex
 *
 *
 * Arg:          shadow_len [UNKN ] Undocumented argument [int]
 * Arg:          has_maxlen [UNKN ] Undocumented argument [int]
 * Arg:              maxlen [UNKN ] Undocumented argument [int]
 * Arg:        shadow_error [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupInterface *]
 *
 */
# line 61 "shadowseqindex.dy"
SeqLookupInterface * new_ShadowSequenceIndex_SeqLookupInterface(int shadow_len,int has_maxlen,int maxlen,int shadow_error)
{
  SeqLookupInterface * out;
  ShadowSequenceIndex * in;

  in = new_ShadowSequenceIndex(26*26*26*26*26,shadow_len,has_maxlen,maxlen,shadow_error);

  out = SeqLookupInterface_alloc_std();

  out->get_client  = get_client_interface_ShadowSequenceIndex;
  out->add_seq     = add_seq_interface_ShadowSequenceIndex;
  out->free_data   = free_interface_ShadowSequenceIndex;
  out->lookup_array_head = NULL;
  out->data = (void*) in;

  return out;
}


/* Function:  get_client_interface_ShadowSequenceIndex(data)
 *
 * Descrip:    gets client interface
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void*]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupClientInterface *]
 *
 */
# line 83 "shadowseqindex.dy"
SeqLookupClientInterface * get_client_interface_ShadowSequenceIndex(void* data)
{
  ShadowSequenceIndex * index = (ShadowSequenceIndex*) data;
  SeqLookupClientInterface * slci;
  ShadowSequenceClient * pcl;

  pcl = ShadowSequenceClient_alloc();

  pcl->result.result = calloc(64,sizeof(SeqLookupResultStruct));
  pcl->result.max_len = 64;
  pcl->result.len = 0;


  pcl->interface.next = next_interface_GenericIndexResult;
  pcl->interface.is_more = is_more_interface_GenericIndexResult;
  pcl->interface.free_data = free_noop_GenericIndexResult;
  pcl->interface.data = (void*) &(pcl->result);

  pcl->index = index;

  slci = SeqLookupClientInterface_alloc();
  slci->lookup       = lookup_interface_ShadowSequenceClient;
  slci->is_populated = is_populated_interface_ShadowSequenceClient;
  slci->free_data    = free_interface_ShadowSequenceClient;
  slci->data = (void*) pcl;

  return slci;
}


/* Function:  lookup_interface_ShadowSequenceClient(data,seq_number)
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
# line 116 "shadowseqindex.dy"
SeqLookupResultInterface * lookup_interface_ShadowSequenceClient(void * data,int seq_number)
{
  ShadowSequenceClient * cli= (ShadowSequenceClient*) data;
  ShadowSequenceIndex * in = cli->index;

  /* reset pointers */
  cli->result.current_pos = 0;
  cli->result.len = 0;

  if( lookup_result_ShadowSeq(&cli->result,in,seq_number) == FALSE ) {
    return NULL;
  }

  return &cli->interface;
}



/* Function:  is_populated_interface_ShadowSequenceClient(data,seq_number)
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
# line 137 "shadowseqindex.dy"
boolean is_populated_interface_ShadowSequenceClient(void * data,int seq_number)
{
  ShadowSequenceClient * client = (ShadowSequenceClient*) data;

  if( client->index->array[seq_number] != NULL ) {
    return TRUE;
  } else {
    return FALSE;
  }
}


/* Function:  add_seq_interface_ShadowSequenceIndex(data,seq,para)
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
# line 152 "shadowseqindex.dy"
boolean add_seq_interface_ShadowSequenceIndex(void * data,Sequence * seq,SeqLookupLoadPara * para)
{
  ShadowSequenceIndex * in = (ShadowSequenceIndex*) data;

  return add_Sequence_ShadowSequenceIndex(in,seq,in->shadow_len);
}

/* Function:  free_interface_ShadowSequenceIndex(data)
 *
 * Descrip:    for interface, frees index
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 *
 */
# line 162 "shadowseqindex.dy"
void free_interface_ShadowSequenceIndex(void * data)
{
  ShadowSequenceIndex * in = (ShadowSequenceIndex*) data;

  free_ShadowSequenceIndex(in);
}


/* Function:  free_interface_ShadowSequenceClient(data)
 *
 * Descrip:    Frees the client data
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 *
 */
# line 173 "shadowseqindex.dy"
void free_interface_ShadowSequenceClient(void * data)
{
  ShadowSequenceClient * cli = (ShadowSequenceClient *)data;

  free(cli->result.result);

  free(cli);

}



/* Function:  lookup_result_ShadowSeq(res,in,seq_no)
 *
 * Descrip:    handles the lookup and storage for a seq_no
 *             lookup
 *
 *
 * Arg:           res [UNKN ] Undocumented argument [GenericIndexResult *]
 * Arg:            in [UNKN ] Undocumented argument [ShadowSequenceIndex *]
 * Arg:        seq_no [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 189 "shadowseqindex.dy"
boolean lookup_result_ShadowSeq(GenericIndexResult * res,ShadowSequenceIndex * in,int seq_no)
{
  int i;
  SingleNumberSequence * sns;
  int pos;

  assert(res);
  assert(in);
  
  if( in->array[seq_no] == NULL ) {
    return FALSE;
  }
  
  for(i=0;i<in->array[seq_no]->current_pos;i++) {
    pos = in->array[seq_no]->seqdb_pos[i];
    sns = lookup_ShadowSequence_SingleNumberSpace(in->space,pos);
    add_result_GenericIndexResult_ShadowSeq(res,sns->seq,pos-sns->start);
  }

  return TRUE;
}

/* Function:  add_result_GenericIndexResult_ShadowSeq(res,seq,pos)
 *
 * Descrip:    adds a particular shadow sequence position, unrolling
 *             shadowed sequences into the result
 *
 *
 * Arg:        res [UNKN ] Undocumented argument [GenericIndexResult *]
 * Arg:        seq [UNKN ] Undocumented argument [ShadowSequence *]
 * Arg:        pos [UNKN ] Undocumented argument [int]
 *
 */
# line 215 "shadowseqindex.dy"
void add_result_GenericIndexResult_ShadowSeq(GenericIndexResult * res,ShadowSequence * seq,int pos)
{
  int i;

  assert(res);
  assert(seq);
  
  add_GenericIndexResult(res,seq->seq,pos);

  for(i=0;i<seq->len;i++) {
    if( pos >= seq->region[i]->start_seq && pos < seq->region[i]->start_seq + seq->region[i]->len ) {
      add_GenericIndexResult(res,seq->region[i]->seq,seq->region[i]->start_shadow + pos - seq->region[i]->start_seq);
    }
  }
  
}



/* Function:  dump_shadow_ShadowSequenceIndex(in,ofp)
 *
 * Descrip:    Dumps information about shadows
 *
 *
 * Arg:         in [UNKN ] Undocumented argument [ShadowSequenceIndex *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 237 "shadowseqindex.dy"
void dump_shadow_ShadowSequenceIndex(ShadowSequenceIndex * in,FILE * ofp)
{
  int i;

  for(i=0;i<in->len;i++) {
    dump_ShadowSequence(in->shadow[i],ofp);
  }

}

/* Function:  dump_stats_ShadowSequenceIndex(in,ofp)
 *
 * Descrip:    Dumps useful information out of shadow sequence array
 *
 *
 * Arg:         in [UNKN ] Undocumented argument [ShadowSequenceIndex *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 250 "shadowseqindex.dy"
void dump_stats_ShadowSequenceIndex(ShadowSequenceIndex * in,FILE * ofp)
{
  int i;
  int j;
  int total_index = 0;
  int total_shadow = 0;
	
  int total_array = 0;
  int total_seq = 0;
  int total_head = 0;

  for(i=0;i<in->array_len;i++) {
    if( in->array[i] != NULL ) {
      total_index += in->array[i]->current_pos;
      total_array += in->array[i]->max;
      total_head++;
    }
  }

  for(i=0;i<in->len;i++) {
    total_seq += in->shadow[i]->seq->len;
    for(j=0;j<in->shadow[i]->len;j++) {
      total_shadow += in->shadow[i]->region[j]->len;
    }
  }

  fprintf(ofp,"Arrayed %d, Shadowed %d (Compression ratio %.2f%%)\n",
	  total_index,
	  total_shadow,
	  (double)(total_shadow)*100/(double)(total_shadow+total_index));
  fprintf(ofp,"Occupied Array Memory %d [%.2f Mbytes]\n",
	  total_index,((total_index/1000000.0)*sizeof(SHADOW_TYPE)));
  fprintf(ofp,"Allocated Array memory %d [%.2fMbytes] [%d %%]\n",
	  total_array,((total_array/1000000.0)*sizeof(SHADOW_TYPE)),
	  total_index*100/total_array);
  fprintf(ofp,"Head memory %d [%.2f Mbytes]\n",total_head,(total_head*sizeof(ShadowArraySeqHead))/100000);

  fprintf(ofp,"Sequence Memory %d [%.2f Mbytes]\n",
	  total_seq,(total_seq/1000000.0)*sizeof(char));

  show_allocator_status_IntAllocatorSet(in->ias,ofp);

}

/* Function:  add_Sequence_ShadowSequenceIndex(in,seq,min_ext)
 *
 * Descrip:    Adds a Sequence to a ShadowIndex, placing shadowed regions
 *             correctly away
 *
 *
 * Arg:             in [UNKN ] Undocumented argument [ShadowSequenceIndex *]
 * Arg:            seq [UNKN ] Undocumented argument [Sequence *]
 * Arg:        min_ext [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 298 "shadowseqindex.dy"
boolean add_Sequence_ShadowSequenceIndex(ShadowSequenceIndex * in,Sequence * seq,int min_ext)
{
  int i;
  int j;
  int no;
  int temp;
  ShadowSequence * shadow = NULL;
  int is_dirty = 0;
  int start_position;

  SingleNumberSequence * sns;
  

  assert(in);
  assert(seq);

  /* reap memory into the static system */
  seq = new_Sequence_StaticSeqHolder(in->ssh,seq);

  shadow = new_ShadowSequence(seq);
  add_ShadowSequenceIndex(in,shadow);
  
  start_position = add_ShadowSequence_SingleNumberSpace(in->space,shadow);
  
  for(i=0;i<seq->len-5;) {

    no = seq_number_aa_5mer(seq->seq+i);

    if( in->array[no] != NULL ) {
      temp = 0;
      for(j=0;j<in->array[no]->current_pos;j++) {

	sns = lookup_ShadowSequence_SingleNumberSpace(in->space,in->array[no]->seqdb_pos[j]);

	/* we can make a hard assumption that sns is not NULL, but just to
	   catch errors a little more explicity... */

	assert(sns);

	if( sns->seq->dirty == 0 && sns->seq != shadow && (temp = add_if_possible_ShadowSequence(sns->seq,seq,min_ext,in->array[no]->seqdb_pos[j]-sns->start,i,in->shadow_error)) != 0 ) {
	  i = temp;
	  
	  break;
	}
      }
      if( temp != 0 ) {
	is_dirty = 1;
	continue;
      }
    }


    if( in->array[no] == NULL ) {
      in->array[no] = new_ShadowArraySeqHead(in->ias);
    }
    add_ShadowArraySeqHead(in->ias,in->array[no],start_position+i);
    i++;
  }


  shadow->dirty = is_dirty;


  return TRUE;
    
}



/* Function:  new_ShadowSequenceIndex(len,shadow_len,has_maxlen,maxlen,shadow_error)
 *
 * Descrip:    New ShadowSequenceIndex
 *
 *
 * Arg:                 len [UNKN ] Undocumented argument [int]
 * Arg:          shadow_len [UNKN ] Undocumented argument [int]
 * Arg:          has_maxlen [UNKN ] Undocumented argument [int]
 * Arg:              maxlen [UNKN ] Undocumented argument [int]
 * Arg:        shadow_error [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [ShadowSequenceIndex *]
 *
 */
# line 370 "shadowseqindex.dy"
ShadowSequenceIndex * new_ShadowSequenceIndex(int len,int shadow_len,int has_maxlen,int maxlen,int shadow_error)
{
  int i;
  ShadowSequenceIndex * out;

  out = ShadowSequenceIndex_alloc_std();
  out->array = calloc(len,sizeof(ShadowArraySeqHead*));
  assert(out->array);

  for(i=0;i<len;i++) {
    out->array[i] = NULL;
  }

  out->array_len = len;
  out->shadow_len = len;
  out->space = new_SingleNumberSpace(has_maxlen,maxlen);
  out->shadow_error = shadow_error;
  out->ias = new_IntAllocatorSet(65);
  out->ssh = new_StaticSeqHolder();
  return out;
}


/* Function:  add_ShadowArraySeqHead(ias,h,seqdb_pos)
 *
 * Descrip:    Adds a sequence/pos pair to an ArrayHead
 *
 *
 * Arg:              ias [UNKN ] Undocumented argument [IntAllocatorSet *]
 * Arg:                h [UNKN ] Undocumented argument [ShadowArraySeqHead *]
 * Arg:        seqdb_pos [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 396 "shadowseqindex.dy"
boolean add_ShadowArraySeqHead(IntAllocatorSet * ias,ShadowArraySeqHead * h,int seqdb_pos)
{

  if( h->current_pos >= h->max ) {
    if( h->max < SHADOW_ARRAYSEQL_LINEAR ) {
      h->seqdb_pos = realloc_intarray_IntAllocatorSet(ias,h->seqdb_pos,h->max,(h->max*2));
      h->max = h->max*2;
    } else {
      h->seqdb_pos = realloc_intarray_IntAllocatorSet(ias,h->seqdb_pos,h->max,h->max + SHADOW_ARRAYSEQL_LINEAR);
      h->max = h->max + SHADOW_ARRAYSEQL_LINEAR;
    }

    if( h->seqdb_pos == NULL ) {
      fatal("ArraySeqLookup realloc failed trying for %d positions\n",h->max);
    }
/*    fprintf(stderr,"... extended to %d\n",h->max); */
  }

  h->seqdb_pos[h->current_pos] = seqdb_pos;

  h->current_pos++;

  return TRUE;
}


/* Function:  new_ShadowArraySeqHead(ias)
 *
 * Descrip:    Builds a new ArraySeqHead structure
 *
 *
 * Arg:        ias [UNKN ] Undocumented argument [IntAllocatorSet *]
 *
 * Return [UNKN ]  Undocumented return value [ShadowArraySeqHead *]
 *
 */
# line 425 "shadowseqindex.dy"
ShadowArraySeqHead * new_ShadowArraySeqHead(IntAllocatorSet * ias)
{
  ShadowArraySeqHead * out;

  out = malloc(sizeof(ShadowArraySeqHead));

  out->seqdb_pos = alloc_intarray_IntAllocatorSet(ias,SHADOW_ARRAYSEQL_BASIC);
  out->max = SHADOW_ARRAYSEQL_BASIC;
  out->current_pos = 0;
  
  return out;
}



# line 507 "shadowseqindex.c"
/* Function:  swap_ShadowSequenceIndex(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_ShadowSequenceIndex
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [ShadowSequence **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_ShadowSequenceIndex(ShadowSequence ** list,int i,int j)  
{
    ShadowSequence * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_ShadowSequenceIndex(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_ShadowSequenceIndex which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [ShadowSequence **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_ShadowSequenceIndex(ShadowSequence ** list,int left,int right,int (*comp)(ShadowSequence * ,ShadowSequence * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_ShadowSequenceIndex(list,left,(left+right)/2);  
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_ShadowSequenceIndex (list,++last,i);    
      }  
    swap_ShadowSequenceIndex (list,left,last);   
    qsort_ShadowSequenceIndex(list,left,last-1,comp);    
    qsort_ShadowSequenceIndex(list,last+1,right,comp);   
}    


/* Function:  sort_ShadowSequenceIndex(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_ShadowSequenceIndex
 *
 *
 * Arg:         obj [UNKN ] Object containing list [ShadowSequenceIndex *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_ShadowSequenceIndex(ShadowSequenceIndex * obj,int (*comp)(ShadowSequence *, ShadowSequence *)) 
{
    qsort_ShadowSequenceIndex(obj->shadow,0,obj->len-1,comp);    
    return;  
}    


/* Function:  expand_ShadowSequenceIndex(obj,len)
 *
 * Descrip:    Really an internal function for add_ShadowSequenceIndex
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [ShadowSequenceIndex *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_ShadowSequenceIndex(ShadowSequenceIndex * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_ShadowSequenceIndex called with no need");    
      return TRUE;   
      }  


    if( (obj->shadow = (ShadowSequence ** ) ckrealloc (obj->shadow,sizeof(ShadowSequence *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_ShadowSequenceIndex, returning FALSE");  
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_ShadowSequenceIndex(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [ShadowSequenceIndex *]
 * Arg:        add [OWNER] Object to add to the list [ShadowSequence *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_ShadowSequenceIndex(ShadowSequenceIndex * obj,ShadowSequence * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_ShadowSequenceIndex(obj,obj->len + ShadowSequenceIndexLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->shadow[obj->len++]=add; 
    return TRUE; 
}    


/* Function:  flush_ShadowSequenceIndex(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [ShadowSequenceIndex *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_ShadowSequenceIndex(ShadowSequenceIndex * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->shadow[i] != NULL)    {  
        free_ShadowSequence(obj->shadow[i]); 
        obj->shadow[i] = NULL;   
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  ShadowSequenceIndex_alloc_std(void)
 *
 * Descrip:    Equivalent to ShadowSequenceIndex_alloc_len(ShadowSequenceIndexLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ShadowSequenceIndex *]
 *
 */
ShadowSequenceIndex * ShadowSequenceIndex_alloc_std(void) 
{
    return ShadowSequenceIndex_alloc_len(ShadowSequenceIndexLISTLENGTH); 
}    


/* Function:  ShadowSequenceIndex_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [ShadowSequenceIndex *]
 *
 */
ShadowSequenceIndex * ShadowSequenceIndex_alloc_len(int len) 
{
    ShadowSequenceIndex * out;  /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = ShadowSequenceIndex_alloc()) == NULL)  
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->shadow = (ShadowSequence ** ) ckcalloc (len,sizeof(ShadowSequence *))) == NULL) {  
      warn("Warning, ckcalloc failed in ShadowSequenceIndex_alloc_len"); 
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_ShadowSequenceIndex(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ShadowSequenceIndex *]
 *
 * Return [UNKN ]  Undocumented return value [ShadowSequenceIndex *]
 *
 */
ShadowSequenceIndex * hard_link_ShadowSequenceIndex(ShadowSequenceIndex * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a ShadowSequenceIndex object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  ShadowSequenceIndex_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ShadowSequenceIndex *]
 *
 */
ShadowSequenceIndex * ShadowSequenceIndex_alloc(void) 
{
    ShadowSequenceIndex * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(ShadowSequenceIndex *) ckalloc (sizeof(ShadowSequenceIndex))) == NULL)  {  
      warn("ShadowSequenceIndex_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->array = NULL;   
    out->array_len = 0;  
    out->shadow = NULL;  
    out->len = out->maxlen = 0;  
    out->space = NULL;   
    out->shadow_len = 0; 
    out->shadow_error = 0;   
    out->ias = NULL; 
    out->ssh = NULL; 


    return out;  
}    


/* Function:  free_ShadowSequenceIndex(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ShadowSequenceIndex *]
 *
 * Return [UNKN ]  Undocumented return value [ShadowSequenceIndex *]
 *
 */
ShadowSequenceIndex * free_ShadowSequenceIndex(ShadowSequenceIndex * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a ShadowSequenceIndex obj. Should be trappable");   
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
    if( obj->shadow != NULL) {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->shadow[i] != NULL)  
          free_ShadowSequence(obj->shadow[i]);   
        }  
      ckfree(obj->shadow);   
      }  
    if( obj->space != NULL)  
      free_SingleNumberSpace(obj->space);    
    if( obj->ias != NULL)    
      free_IntAllocatorSet(obj->ias);    
    if( obj->ssh != NULL)    
      free_StaticSeqHolder(obj->ssh);    


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
