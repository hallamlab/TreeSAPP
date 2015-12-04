#ifdef _cplusplus
extern "C" {
#endif
#include "sequence_codon.h"

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
# line 29 "sequence_codon.dy"
Sequence * reverse_complement_Sequence(Sequence * seq)
{
  Sequence * out;
  int i;
  int j;


  if( is_dna_Sequence(seq) != TRUE ) {
    warn("Cannot reverse complement non-DNA sequence... type is %s",Sequence_type_to_string(seq->type));
    return NULL;
  }

  out = Sequence_from_static_memory(seq->name,seq->seq);
  
  for(j=0,i=seq->len-1;i >= 0;i--,j++) {
    out->seq[j] = char_complement_base(seq->seq[i]);
    /*fprintf(stderr,"In position %d placed %c from %c\n",j,out->seq[j],seq->seq[i]);*/
  }

  out->len = strlen(seq->seq);

  out->offset = seq->end;
  out->end   = seq->offset;
  out->type  = seq->type;
  

  return out;
}
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
# line 74 "sequence_codon.dy"
Sequence * magic_trunc_Sequence(Sequence * seq,int start,int end)
{
  Sequence * temp;
  Sequence * out;

  if( is_dna_Sequence(seq) == FALSE) {
    warn("Cannot magic truncate on a non DNA sequence... type is %s",Sequence_type_to_string(seq->type));
    return NULL;
  }

  if( start < 0 || end < 0 ) {
    warn("Attempting a magic truncation on indices which are less than zero [%d:%d]. Clearly impossible",start,end);
    return NULL;
  }

  if( start < end ) {
    return trunc_Sequence(seq,start,end);
  }
  else {
    temp = trunc_Sequence(seq,end,start);
    if( temp == NULL ) {
      warn("Unable to truncate sequence");
      return NULL;
    }


    out = reverse_complement_Sequence(temp);

    free_Sequence(temp);

    return out;
  }

}

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
# line 119 "sequence_codon.dy"
Sequence * translate_Sequence(Sequence * dna,CodonTable * ct)
{
  Sequence * out;
  int i;
  int j;
  int len;
  char * seq;
  char * name;
  char buffer[512];

  if( is_dna_Sequence(dna) == FALSE) {
    warn("Trying to make a translation from a non DNA sequence... type is [%s]",Sequence_type_to_string(dna->type));
    return NULL;
  }

  len = dna->len/3 + 1;
  
  seq = ckcalloc(len,sizeof(char));

  sprintf(buffer,"%s.tr",dna->name == NULL ? "NoNameDNASeq" : dna->name);

  name = stringalloc(buffer);

  out  = Sequence_from_dynamic_memory(name,seq);

  for(i=0,j=0;i<dna->len-3;i+=3,j++) {
    out->seq[j] = aminoacid_from_seq(ct,dna->seq+i);
  }
  out->seq[j] = '\0';

  out->type  = SEQUENCE_PROTEIN;
  out->len   = strlen(out->seq);

  return out;
}

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
# line 174 "sequence_codon.dy"
boolean force_to_dna_Sequence(Sequence * seq,double fraction,int * number_of_conver)
{
  int count =0;
  int i;

  if( seq == NULL ) {
    warn("Attempting to force a sequence with no Sequence object!\n");
    return FALSE;
  } 
  if( seq->len <= 0 ) {
    warn("Trying to make a sequence with a length of %d. Bad news!",seq->len);
    return FALSE;
  }

  
  for(i=0;i<seq->len;i++) {
    /* if it is lower case, uppercase it! */
    seq->seq[i] = (char)toupper((int)seq->seq[i]);
    
    if( !is_valid_base_char(seq->seq[i]) ) {
      count++;
    }
  }


  if( ((double)count/(double)seq->len) < fraction ) {
    seq->type = SEQUENCE_DNA;
    if( count != 0 ) {
      for(i=0;i<seq->len;i++) {
	if( !is_valid_base_char(seq->seq[i]) ) {
	  seq->seq[i] = 'N';
	}
      }
    }
    if( number_of_conver != NULL ) {
      *number_of_conver = count;
    }
    return TRUE;
  }
  else {
    if( number_of_conver != NULL ) {
      *number_of_conver = count;
    }
    return FALSE;
  }
}
# line 216 "sequence_codon.c"

#ifdef _cplusplus
}
#endif
