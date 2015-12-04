#ifndef DYNAMITEassembly_stream_fastaHEADERFILE
#define DYNAMITEassembly_stream_fastaHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "assembly_stream_interface.h"





    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
AssemblySequenceStream * Wise2_plain_fasta_AssemblySequenceStream(FILE * ifp);
#define plain_fasta_AssemblySequenceStream Wise2_plain_fasta_AssemblySequenceStream
AssemblySequence * Wise2_next_AssemblySequence_fasta_impl(void * h);
#define next_AssemblySequence_fasta_impl Wise2_next_AssemblySequence_fasta_impl


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/

#ifdef _cplusplus
}
#endif

#endif
