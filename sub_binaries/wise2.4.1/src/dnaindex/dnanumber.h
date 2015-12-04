#ifndef DYNAMITEdnanumberHEADERFILE
#define DYNAMITEdnanumberHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "wisebase.h"
#include "codon.h"
#include "sequence.h"

typedef struct {
  int number;
  char flipped;
} DnaNumber;

struct Wise2_DnaNumberSequence {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    DnaNumber * seq;     
    int len;     
    Sequence * orig;     
    } ;  
/* DnaNumberSequence defined */ 
#ifndef DYNAMITE_DEFINED_DnaNumberSequence
typedef struct Wise2_DnaNumberSequence Wise2_DnaNumberSequence;
#define DnaNumberSequence Wise2_DnaNumberSequence
#define DYNAMITE_DEFINED_DnaNumberSequence
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  free_DnaNumber(dn)
 *
 * Descrip:    dummy free
 *
 *
 * Arg:        dn [UNKN ] Undocumented argument [DnaNumber *]
 *
 * Return [UNKN ]  Undocumented return value [DnaNumber *]
 *
 */
DnaNumber * Wise2_free_DnaNumber(DnaNumber * dn);
#define free_DnaNumber Wise2_free_DnaNumber


/* Function:  hard_link_DnaNumberSequence(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DnaNumberSequence *]
 *
 * Return [UNKN ]  Undocumented return value [DnaNumberSequence *]
 *
 */
DnaNumberSequence * Wise2_hard_link_DnaNumberSequence(DnaNumberSequence * obj);
#define hard_link_DnaNumberSequence Wise2_hard_link_DnaNumberSequence


/* Function:  DnaNumberSequence_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaNumberSequence *]
 *
 */
DnaNumberSequence * Wise2_DnaNumberSequence_alloc(void);
#define DnaNumberSequence_alloc Wise2_DnaNumberSequence_alloc


/* Function:  free_DnaNumberSequence(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DnaNumberSequence *]
 *
 * Return [UNKN ]  Undocumented return value [DnaNumberSequence *]
 *
 */
DnaNumberSequence * Wise2_free_DnaNumberSequence(DnaNumberSequence * obj);
#define free_DnaNumberSequence Wise2_free_DnaNumberSequence


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
DnaNumberSequence * Wise2_new_DnaNumberSequence(Sequence * seq,int nmer_size);
#define new_DnaNumberSequence Wise2_new_DnaNumberSequence
char Wise2_first_char_from_dnanumber(int dnanumber,int nmer_size,int flipped);
#define first_char_from_dnanumber Wise2_first_char_from_dnanumber
DnaNumber Wise2_dna_number_from_string(char * str,int nmer_size);
#define dna_number_from_string Wise2_dna_number_from_string


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/

#ifdef _cplusplus
}
#endif

#endif
