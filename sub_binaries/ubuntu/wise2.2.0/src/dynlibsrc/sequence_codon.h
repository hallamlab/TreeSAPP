#ifndef DYNAMITEsequence_codonHEADERFILE
#define DYNAMITEsequence_codonHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "sequence.h"
#include "codon.h"




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  reverse_complement_Sequence(seq)
 *
 * Descrip:    This both complements and reverses a sequence,
 *             - a common wish!
 *
 *             The start/end are correct with respect to the start/end
 *             of the sequence (ie start = end, end = start).
 *
 *
 * Arg:        seq [READ ] Sequence to that is used to reverse (makes a new Sequence) [Sequence *]
 *
 * Return [OWNER]  new Sequence which is reversed [Sequence *]
 *
 */
Sequence * Wise2_reverse_complement_Sequence(Sequence * seq);
#define reverse_complement_Sequence Wise2_reverse_complement_Sequence


/* Function:  magic_trunc_Sequence(seq,start,end)
 *
 * Descrip:    Clever function for dna sequences.
 *
 *             When start < end, truncates normally
 *
 *             when start > end, truncates end,start and then
 *             reverse complements.
 *
 *             ie. If you have a coordinate system where reverse 
 *             sequences are labelled in reverse start/end way,
 *             then this routine produces the correct sequence.
 *
 *
 * Arg:          seq [READ ] sequence that is the source to be truncated [Sequence *]
 * Arg:        start [READ ] start point [int]
 * Arg:          end [READ ] end point [int]
 *
 * Return [OWNER]  new Sequence which is truncated/reversed [Sequence *]
 *
 */
Sequence * Wise2_magic_trunc_Sequence(Sequence * seq,int start,int end);
#define magic_trunc_Sequence Wise2_magic_trunc_Sequence


/* Function:  translate_Sequence(dna,ct)
 *
 * Descrip:    This translates a DNA sequence to a protein.
 *             It assummes that it starts at first residue
 *             (use trunc_Sequence to chop a sequence up).
 *
 *
 * Arg:        dna [READ ] DNA sequence to be translated [Sequence *]
 * Arg:         ct [READ ] Codon table to do codon->aa mapping [CodonTable *]
 *
 * Return [OWNER]  new protein sequence [Sequence *]
 *
 */
Sequence * Wise2_translate_Sequence(Sequence * dna,CodonTable * ct);
#define translate_Sequence Wise2_translate_Sequence


/* Function:  force_to_dna_Sequence(seq,fraction,number_of_conver)
 *
 * Descrip:    This 
 *              a) sees how many non ATGCN characters there are in Seq
 *              b) If the level is below fraction
 *                 a) flips non ATGC chars to N
 *                 b) writes number of conversions to number_of_conver
 *                 c) returns TRUE
 *              c) else returns FALSE
 *
 *             fraction of 0.0 means completely intolerant of errors
 *             fraction of 1.0 means completely tolerant of errors
 *
 *
 *
 * Arg:                     seq [RW   ] sequence object read and converted  [Sequence *]
 * Arg:                fraction [READ ] number 0..1 for tolerance of conversion [double]
 * Arg:        number_of_conver [WRITE] number of conversions actually made [int *]
 *
 * Return [READ ]  TRUE for conversion to DNA, FALSE if not [boolean]
 *
 */
boolean Wise2_force_to_dna_Sequence(Sequence * seq,double fraction,int * number_of_conver);
#define force_to_dna_Sequence Wise2_force_to_dna_Sequence


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/

#ifdef _cplusplus
}
#endif

#endif
