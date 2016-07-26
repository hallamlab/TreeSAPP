

/* Helper functions in the module
 *
 * Wise2_reverse_complement_Sequence
 * Wise2_magic_trunc_Sequence
 * Wise2_translate_Sequence
 *



/* These functions are not associated with an object */
/* Function:  Wise2_reverse_complement_Sequence(seq)
 *
 * Descrip:    This both complements and reverses a sequence,
 *             - a common wish!
 *
 *             The start/end are correct with respect to the start/end
 *             of the sequence (ie start = end, end = start).
 *
 *
 * Arg:        seq          Sequence to that is used to reverse (makes a new Sequence) [Wise2_Sequence *]
 *
 * Returns new Sequence which is reversed [Wise2_Sequence *]
 *
 */
Wise2_Sequence * Wise2_reverse_complement_Sequence( Wise2_Sequence * seq);

/* Function:  Wise2_magic_trunc_Sequence(seq,start,end)
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
 * Arg:        seq          sequence that is the source to be truncated [Wise2_Sequence *]
 * Arg:        start        start point [int]
 * Arg:        end          end point [int]
 *
 * Returns new Sequence which is truncated/reversed [Wise2_Sequence *]
 *
 */
Wise2_Sequence * Wise2_magic_trunc_Sequence( Wise2_Sequence * seq,int start,int end);

/* Function:  Wise2_translate_Sequence(dna,ct)
 *
 * Descrip:    This translates a DNA sequence to a protein.
 *             It assummes that it starts at first residue
 *             (use trunc_Sequence to chop a sequence up).
 *
 *
 * Arg:        dna          DNA sequence to be translated [Wise2_Sequence *]
 * Arg:        ct           Codon table to do codon->aa mapping [Wise2_CodonTable *]
 *
 * Returns new protein sequence [Wise2_Sequence *]
 *
 */
Wise2_Sequence * Wise2_translate_Sequence( Wise2_Sequence * dna,Wise2_CodonTable * ct);

