#ifndef DYNAMITEkmer_hashHEADERFILE
#define DYNAMITEkmer_hashHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "kmer_index_interface.h"

#define KMER_HASH_BITS    24
#define KMER_HASH_MASK    0x0000000000ffffffL
#define KMER_HASH_BUCKETS 16777216

/*
#define KMER_HASH_BITS    28
#define KMER_HASH_MASK    0x000000000fffffffL
#define KMER_HASH_BUCKETS (16777216*16)
*/

typedef struct KmerHashBucket {
  void   **index;
  kmer_t  *kmer;
  int      length;
  int      buflen;
} KmerHashBucket;

typedef struct KmerHashIndex {
  KmerHashBucket *bucket;
  kmer_t   min_kmer;
  kmer_t   max_kmer;
} KmerHashIndex;




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
int Wise2_crc32(int crc, const char *buf, int len);
#define crc32 Wise2_crc32
KmerIndexInterface * Wise2_new_KmerHashIndex(int kmer_size);
#define new_KmerHashIndex Wise2_new_KmerHashIndex
kmer_t Wise2_hash_kmer(kmer_t kmer);
#define hash_kmer Wise2_hash_kmer
void * Wise2_retrieve_by_kmer_KmerHashIndex(void * handle,kmer_t kmer);
#define retrieve_by_kmer_KmerHashIndex Wise2_retrieve_by_kmer_KmerHashIndex
kmer_t Wise2_next_filled_kmer_KmerHashIndex(void * handle,kmer_t kmer);
#define next_filled_kmer_KmerHashIndex Wise2_next_filled_kmer_KmerHashIndex
boolean Wise2_insert_by_kmer_KmerHashIndex(void * handle,kmer_t kmer,void * poi);
#define insert_by_kmer_KmerHashIndex Wise2_insert_by_kmer_KmerHashIndex
void Wise2_free_KmerHashIndex(void * handle);
#define free_KmerHashIndex Wise2_free_KmerHashIndex


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/

#ifdef _cplusplus
}
#endif

#endif
