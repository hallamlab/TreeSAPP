#ifndef DYNAMITEkmer_index_interfaceHEADERFILE
#define DYNAMITEkmer_index_interfaceHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"

  /* this brings in kmer_t */
#include "dnamapping.h"



typedef struct KmerIndexInterface {
  void*   (*retrieve_by_kmer)(void *,kmer_t); 
  boolean (*insert_by_kmer)(void*,kmer_t,void*); 
  void    (*free_handle)(void*); 
  kmer_t  (*next_filled_kmer)(void *,kmer_t);
  void * handle;
  int kmer_size;
} KmerIndexInterface;




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
void Wise2_free_KmerIndexInterface(KmerIndexInterface * kii);
#define free_KmerIndexInterface Wise2_free_KmerIndexInterface


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/

#ifdef _cplusplus
}
#endif

#endif
