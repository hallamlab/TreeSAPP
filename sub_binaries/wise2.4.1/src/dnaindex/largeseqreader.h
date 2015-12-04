#ifndef DYNAMITElargeseqreaderHEADERFILE
#define DYNAMITElargeseqreaderHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "sequence.h"

#define LINEAR_LARGEFASTA_READ 50000000


    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
Sequence * Wise2_read_large_dna_Sequence_file(char * filename,int report,FILE * logfp);
#define read_large_dna_Sequence_file Wise2_read_large_dna_Sequence_file
Sequence * Wise2_read_large_dna_Sequence(FILE * ifp,int report,FILE * logfp);
#define read_large_dna_Sequence Wise2_read_large_dna_Sequence


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/

#ifdef _cplusplus
}
#endif

#endif
