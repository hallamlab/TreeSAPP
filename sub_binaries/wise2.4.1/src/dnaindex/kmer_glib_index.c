#ifdef _cplusplus
extern "C" {
#endif
#include "kmer_glib_index.h"







# line 28 "kmer_glib_index.dy"
KmerIndexInterface * new_interface_KmerGlibIndex(int kmer_size)
{
  KmerIndexInterface * out;
  KmerGlibIndex * kgi;

  assert(kmer_size < 16);
  

  kgi = new_KmerGlibIndex();

  out = malloc(sizeof(KmerIndexInterface));

  out->retrieve_by_kmer = retrieve_by_kmer_KmerGlibIndex;
  out->insert_by_kmer   = insert_by_kmer_KmerGlibIndex;
  out->free_handle      = free_KmerGlibIndex;
  out->next_filled_kmer = next_filled_kmer_KmerGlibIndex;

  out->handle = (void*) kgi;
  out->kmer_size = kmer_size;


  return out;

}

# line 53 "kmer_glib_index.dy"
kmer_t next_filled_kmer_KmerGlibIndex(void * handle,kmer_t prev)
{
  KmerGlibIndex * kgi = (KmerGlibIndex*) handle;
  
  if( prev == -1 ) {
    if( kgi->active_index == NULL ) {
      make_active_kmer_KmerGlibIndex(kgi);
    }
    kgi->kmer_pos = 0;
  }
  
  if( kgi->kmer_pos < kgi->active_index_length ) {
    return kgi->active_index[kgi->kmer_pos++];
  } else {
    return -1;
  }

}


# line 73 "kmer_glib_index.dy"
void * retrieve_by_kmer_KmerGlibIndex(void * handle,kmer_t kmer)
{
  KmerGlibIndex * kgi = (KmerGlibIndex*) handle;

  return (void*) g_hash_table_lookup(kgi->hash,(gconstpointer)kmer);
}

# line 80 "kmer_glib_index.dy"
boolean insert_by_kmer_KmerGlibIndex(void * handle,kmer_t kmer,void * poi)
{
  KmerGlibIndex * kgi = (KmerGlibIndex*) handle;

  g_hash_table_insert(kgi->hash,(gpointer)kmer,poi);
  
  return TRUE;
}

# line 89 "kmer_glib_index.dy"
void free_KmerGlibIndex(void * handle)
{
  KmerGlibIndex * kgi = (KmerGlibIndex*) handle;

  assert(kgi);
  
  warn("Not freeing glib hash, so leaking memory here");

  
  free(kgi);
}


# line 102 "kmer_glib_index.dy"
boolean make_active_kmer_KmerGlibIndex(void * handle)
{
  KmerGlibIndex * kgi = (KmerGlibIndex*) handle;
  int size;

  assert(kgi);

  size = g_hash_table_size(kgi->hash);
  kgi->active_index = calloc(size,sizeof(kmer_t));
  kgi->active_index_length = size;
  kgi->kmer_pos = 0;
  g_hash_table_foreach(kgi->hash,retrieve_active_kmer_KmerGlibIndex,(gpointer)kgi);

  kgi->kmer_pos = 0;

  return TRUE;
}

# line 120 "kmer_glib_index.dy"
void retrieve_active_kmer_KmerGlibIndex(gpointer key,gpointer value,gpointer user_data)
{
  KmerGlibIndex * kgi = (KmerGlibIndex*) user_data;

  kgi->active_index[kgi->kmer_pos++] = (kmer_t) key;
}


# line 128 "kmer_glib_index.dy"
KmerGlibIndex * new_KmerGlibIndex(void) 
{
  KmerGlibIndex * out;

  out = malloc(sizeof(KmerGlibIndex));
  out->hash = g_hash_table_new(g_direct_hash,g_direct_equal);
  out->active_index = NULL;
  out->kmer_pos = 0;
  out->active_index_length = 0;

  return out;
}

# line 132 "kmer_glib_index.c"

#ifdef _cplusplus
}
#endif
