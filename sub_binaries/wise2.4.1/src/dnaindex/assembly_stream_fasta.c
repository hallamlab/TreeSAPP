#ifdef _cplusplus
extern "C" {
#endif
#include "assembly_stream_fasta.h"


# line 17 "assembly_stream_fasta.dy"
AssemblySequenceStream * plain_fasta_AssemblySequenceStream(FILE * ifp)
{
  AssemblySequenceStream * out;

  out = AssemblySequenceStream_alloc();

  out->handle = (void*) ifp;

  out->next_AssemblySequence = next_AssemblySequence_fasta_impl;
  out->free_handle = NULL; /* it is a no-op */

  return out;
}


# line 32 "assembly_stream_fasta.dy"
AssemblySequence * next_AssemblySequence_fasta_impl(void * h)
{
  FILE * ifp = (FILE *) h;

  return read_plain_fasta_AssemblySequence(ifp,0,NULL);
  
}


# line 32 "assembly_stream_fasta.c"

#ifdef _cplusplus
}
#endif
