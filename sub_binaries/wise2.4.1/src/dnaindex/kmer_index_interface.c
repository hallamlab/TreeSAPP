#ifdef _cplusplus
extern "C" {
#endif
#include "kmer_index_interface.h"

# line 26 "kmer_index_interface.dy"
void free_KmerIndexInterface(KmerIndexInterface * kii)
{
  assert(kii);

  (*kii->free_handle)(kii->handle);

  free(kii);

}






# line 21 "kmer_index_interface.c"

#ifdef _cplusplus
}
#endif
