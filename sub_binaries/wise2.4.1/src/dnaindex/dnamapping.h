#ifndef DYNAMITEdnamappingHEADERFILE
#define DYNAMITEdnamappingHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"

/*typedef long int kmer_t;*/
typedef long long kmer_t;

struct lexical_kmer {
  kmer_t kmer;
  char lexical_forward;
};



    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
char * Wise2_map_to_basepair_numbers(char * seq,long int len);
#define map_to_basepair_numbers Wise2_map_to_basepair_numbers
kmer_t Wise2_reverse_complement_dna_number(kmer_t number,int kmer_size);
#define reverse_complement_dna_number Wise2_reverse_complement_dna_number
void Wise2_reverse_map_dna_number(kmer_t number,int kmer_size,char * buffer);
#define reverse_map_dna_number Wise2_reverse_map_dna_number
char Wise2_lexical_last_base_from_kmer(kmer_t kmer,int kmer_size,int lexical_for);
#define lexical_last_base_from_kmer Wise2_lexical_last_base_from_kmer
struct lexical_kmer Wise2_lexical_dna_number_from_string(char * str,int kmer_size);
#define lexical_dna_number_from_string Wise2_lexical_dna_number_from_string
kmer_t Wise2_forward_dna_number_from_string(char * str,int kmer_size);
#define forward_dna_number_from_string Wise2_forward_dna_number_from_string


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/

#ifdef _cplusplus
}
#endif

#endif
