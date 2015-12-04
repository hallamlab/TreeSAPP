#ifdef _cplusplus
extern "C" {
#endif
#include "kmer_direct.h"

# line 18 "kmer_direct.dy"
KmerIndexInterface * new_KmerDirectIndex(int kmer_size)
{
  KmerIndexInterface * out;
  KmerDirectIndex * kdi;
  kmer_t len = 1;
  int i;
  
  for(i=0;i<kmer_size;i++) {
    len = len *4;
  }

  kdi = malloc(sizeof(KmerDirectIndex));
  kdi->index = calloc(len,sizeof(void*));
  if( kdi->index == NULL ) {
    warn("could not build index of size %lld pointers\n",len);
    return NULL;
  }
  kdi->index_length = len;

  out = malloc(sizeof(KmerIndexInterface));

  out->retrieve_by_kmer = retrieve_by_kmer_KmerDirectIndex;
  out->insert_by_kmer = insert_by_kmer_KmerDirectIndex;
  out->free_handle = free_KmerDirectIndex;
  out->next_filled_kmer = next_filled_kmer_KmerDirectIndex;

  out->handle = (void*) kdi;
  out->kmer_size = kmer_size;


  return out;

}

# line 52 "kmer_direct.dy"
kmer_t next_filled_kmer_KmerDirectIndex(void * handle,kmer_t kmer)
{

  KmerDirectIndex * ind = (KmerDirectIndex*) handle;
  assert(ind);

  for(kmer++;kmer < ind->index_length;kmer++) {
    if( ind->index[kmer] != NULL ) {
      return kmer;
    }
  }

  return -1;
}


# line 68 "kmer_direct.dy"
void * retrieve_by_kmer_KmerDirectIndex(void * handle,kmer_t kmer)
{
  KmerDirectIndex * ind = (KmerDirectIndex*) handle;
  assert(ind != NULL);
  assert(kmer < ind->index_length);

  return ind->index[kmer];
}


# line 78 "kmer_direct.dy"
boolean insert_by_kmer_KmerDirectIndex(void * handle,kmer_t kmer,void * poi)
{
  KmerDirectIndex * in = (KmerDirectIndex*) handle;
  assert(in != NULL);
  assert(kmer < in->index_length);

  in->index[kmer] = poi;
  return TRUE;
}


# line 89 "kmer_direct.dy"
void free_KmerDirectIndex(void * handle)
{
  KmerDirectIndex * in = (KmerDirectIndex*) handle;
  assert(in);
  free(in->index);
  free(in);
  return;
}
# line 89 "kmer_direct.c"

#ifdef _cplusplus
}
#endif
