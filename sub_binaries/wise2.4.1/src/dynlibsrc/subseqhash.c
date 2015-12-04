#ifdef _cplusplus
extern "C" {
#endif
#include "subseqhash.h"



/* Function:  new_ghash_SeqLookupInterface(void)
 *
 * Descrip:    Makes a new Glib based SeqLookup system
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupInterface *]
 *
 */
# line 24 "subseqhash.dy"
SeqLookupInterface * new_ghash_SeqLookupInterface(void)
{
  SeqLookupInterface * out;

  out = SeqLookupInterface_alloc_std();

  out->data = (void*) g_hash_table_new(g_direct_hash,g_direct_equal);
  out->get_client = get_client_subseqhash_ghash;
  out->add_seq = add_seq_subseqhash_ghash;
  out->add_direct_number = add_direct_number_subseqhash_ghash;
  out->free_data = free_subseqhash_ghash;

  return out;
}

/* Function:  get_client_subseqhash_ghash(data)
 *
 * Descrip:    provides the interface defiition for get_client
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupClientInterface *]
 *
 */
# line 42 "subseqhash.dy"
SeqLookupClientInterface * get_client_subseqhash_ghash(void * data)
{
  SeqLookupClientInterface * sci;

  GHashTable * t = (GHashTable*) data;

  sci = SeqLookupClientInterface_alloc();
  sci->is_populated = is_populated_subseqhash_ghash;
  sci->lookup = lookup_subseqhash_ghash;
  sci->free_data = free_client_subseqhash_ghash;
  sci->data = t;



  return sci;
}

/* Function:  free_client_subseqhash_ghash(data)
 *
 * Descrip:    provides the interface definition for free_data on client, which is
 *             actually a no-op in this case
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 *
 */
# line 63 "subseqhash.dy"
void free_client_subseqhash_ghash(void * data)
{
  return;
}
  

/* Function:  free_subseqhash_ghash(data)
 *
 * Descrip:    Internal function for freeing ghash
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 *
 */
# line 72 "subseqhash.dy"
void free_subseqhash_ghash(void * data)
{
  GHashTable * h;
  
  h = (GHashTable *) data;
	
  g_hash_table_foreach_remove(h,remove_subseq_SeqLookupPos,NULL);

  g_hash_table_destroy(h);

}

/* Function:  remove_subseq_SeqLookupPos(key,value,user_data)
 *
 * Descrip:    Sub call of free
 *
 *
 * Arg:              key [UNKN ] Undocumented argument [gpointer]
 * Arg:            value [UNKN ] Undocumented argument [gpointer]
 * Arg:        user_data [UNKN ] Undocumented argument [gpointer]
 *
 * Return [UNKN ]  Undocumented return value [gboolean]
 *
 */
# line 87 "subseqhash.dy"
gboolean remove_subseq_SeqLookupPos(gpointer key,gpointer value,gpointer user_data)
{
  SeqLookupPos * s;

  s = (SeqLookupPos *) user_data;

  if( s != NULL ) {
    free_SeqLookupPos(s);
  }

  return TRUE;
}

/* Function:  is_populated_subseqhash_ghash(data,seq_number)
 *
 * Descrip:    Tells whether a position is populated or not
 *
 *
 * Arg:              data [UNKN ] Undocumented argument [void *]
 * Arg:        seq_number [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 103 "subseqhash.dy"
boolean is_populated_subseqhash_ghash(void * data, int seq_number)
{
  GHashTable * h;
  
  h = (GHashTable *) data;

  /*  fprintf(stdout,"Looking up with %d\n",h);*/

  if( g_hash_table_lookup(h,(gconstpointer)seq_number) == NULL ) {
    return FALSE;
  } else {
    return TRUE;
  }
}


/* Function:  lookup_subseqhash_ghash(data,seq_number)
 *
 * Descrip:    Retrieves a SeqLookupPos from the hash
 *
 *
 * Arg:              data [UNKN ] Undocumented argument [void *]
 * Arg:        seq_number [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupResultInterface *]
 *
 */
# line 122 "subseqhash.dy"
SeqLookupResultInterface * lookup_subseqhash_ghash(void * data,int seq_number)
{
  GHashTable * h;
  
  h = (GHashTable *) data;

  return new_linkedl_SeqLookupResultInterface((SeqLookupPos *)g_hash_table_lookup(h,(gconstpointer)seq_number));

}

/* Function:  add_seq_subseqhash_ghash(data,seq,para)
 *
 * Descrip:    Adds a sequence/pos pair to the hash with this
 *             number
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 * Arg:         seq [UNKN ] Undocumented argument [Sequence *]
 * Arg:        para [UNKN ] Undocumented argument [SeqLookupLoadPara *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 136 "subseqhash.dy"
boolean add_seq_subseqhash_ghash(void * data,Sequence * seq,SeqLookupLoadPara * para)
{
  GHashTable * h;
  SeqLookupPos * p;
  SeqLookupPos * ret;
  int i;
  int seq_number;

  h = (GHashTable *) data;

  assert(para->tile_freq > 0);
  if( para->tile_freq > 100 ) {
    info("Tile frequency load greater than 100. That's a little odd");
  }


  for(i=0;i<seq->len-5;i = i+para->tile_freq) {
    seq_number = seq_number_aa_5mer(seq->seq+i);
    
    p = SeqLookupPos_alloc();
    p->seq = seq;
    p->pos = i;
    p->next = NULL;
    
    if((ret = (SeqLookupPos *) g_hash_table_lookup(h,(gconstpointer)seq_number)) == NULL ) {
      g_hash_table_insert(h,(gpointer)seq_number,p);
    } else {
      p->next = ret->next;
      ret->next = p;
    }
  }

  return TRUE;
}

/* Function:  add_direct_number_subseqhash_ghash(data,seq_number,target,pos)
 *
 * Descrip:    Adds a direct no/seq pair to the system
 *             number
 *
 *
 * Arg:              data [UNKN ] Undocumented argument [void *]
 * Arg:        seq_number [UNKN ] Undocumented argument [int]
 * Arg:            target [UNKN ] Undocumented argument [Sequence *]
 * Arg:               pos [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 175 "subseqhash.dy"
boolean add_direct_number_subseqhash_ghash(void * data,int seq_number,Sequence * target,int pos)
{
  GHashTable * h;
  SeqLookupPos * p;
  SeqLookupPos * ret;

  h = (GHashTable *) data;

  p = SeqLookupPos_alloc();
  p->seq = target;
  p->pos = pos;
  p->next = NULL;

  if((ret = (SeqLookupPos *) g_hash_table_lookup(h,(gconstpointer)seq_number)) == NULL ) {
    g_hash_table_insert(h,(gpointer)seq_number,p);
  } else {
    p->next = ret->next;
    ret->next = p;
  }

  return TRUE;
}



# line 244 "subseqhash.c"

#ifdef _cplusplus
}
#endif
