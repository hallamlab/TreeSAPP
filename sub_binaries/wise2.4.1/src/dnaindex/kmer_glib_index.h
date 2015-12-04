#ifndef DYNAMITEkmer_glib_indexHEADERFILE
#define DYNAMITEkmer_glib_indexHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "kmer_index_interface.h"
#include "glib.h"


typedef struct KmerGlibIndex {
  GHashTable * hash;
  kmer_t index_length;
  kmer_t * active_index; 
  int kmer_pos;
  int active_index_length;
}  KmerGlibIndex;



    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
KmerIndexInterface * Wise2_new_interface_KmerGlibIndex(int kmer_size);
#define new_interface_KmerGlibIndex Wise2_new_interface_KmerGlibIndex
kmer_t Wise2_next_filled_kmer_KmerGlibIndex(void * handle,kmer_t prev);
#define next_filled_kmer_KmerGlibIndex Wise2_next_filled_kmer_KmerGlibIndex
void * Wise2_retrieve_by_kmer_KmerGlibIndex(void * handle,kmer_t kmer);
#define retrieve_by_kmer_KmerGlibIndex Wise2_retrieve_by_kmer_KmerGlibIndex
boolean Wise2_insert_by_kmer_KmerGlibIndex(void * handle,kmer_t kmer,void * poi);
#define insert_by_kmer_KmerGlibIndex Wise2_insert_by_kmer_KmerGlibIndex
void Wise2_free_KmerGlibIndex(void * handle);
#define free_KmerGlibIndex Wise2_free_KmerGlibIndex
boolean Wise2_make_active_kmer_KmerGlibIndex(void * handle);
#define make_active_kmer_KmerGlibIndex Wise2_make_active_kmer_KmerGlibIndex
void Wise2_retrieve_active_kmer_KmerGlibIndex(gpointer key,gpointer value,gpointer user_data);
#define retrieve_active_kmer_KmerGlibIndex Wise2_retrieve_active_kmer_KmerGlibIndex
KmerGlibIndex * Wise2_new_KmerGlibIndex(void) ;
#define new_KmerGlibIndex Wise2_new_KmerGlibIndex


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/

#ifdef _cplusplus
}
#endif

#endif
