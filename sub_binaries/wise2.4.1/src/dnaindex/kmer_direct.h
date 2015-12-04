#ifndef DYNAMITEkmer_directHEADERFILE
#define DYNAMITEkmer_directHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "kmer_index_interface.h"


typedef struct KmerDirectIndex {
  void ** index;
  kmer_t index_length;
} KmerDirectIndex;



    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
KmerIndexInterface * Wise2_new_KmerDirectIndex(int kmer_size);
#define new_KmerDirectIndex Wise2_new_KmerDirectIndex
kmer_t Wise2_next_filled_kmer_KmerDirectIndex(void * handle,kmer_t kmer);
#define next_filled_kmer_KmerDirectIndex Wise2_next_filled_kmer_KmerDirectIndex
void * Wise2_retrieve_by_kmer_KmerDirectIndex(void * handle,kmer_t kmer);
#define retrieve_by_kmer_KmerDirectIndex Wise2_retrieve_by_kmer_KmerDirectIndex
boolean Wise2_insert_by_kmer_KmerDirectIndex(void * handle,kmer_t kmer,void * poi);
#define insert_by_kmer_KmerDirectIndex Wise2_insert_by_kmer_KmerDirectIndex
void Wise2_free_KmerDirectIndex(void * handle);
#define free_KmerDirectIndex Wise2_free_KmerDirectIndex


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/

#ifdef _cplusplus
}
#endif

#endif
