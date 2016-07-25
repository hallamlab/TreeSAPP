#include "seqerror.h"

int main(int argc,char ** argv)
{
  Sequence * seq;
  ErrorSequence * eseq;
  int i;
  boolean just_fasta = FALSE;
  int number = 10;
  Probability error = 0.01;

  if( strip_out_boolean_argument(&argc,argv,"fasta") ) {
    just_fasta = TRUE;
  }
  
  strip_out_integer_argument(&argc,argv,"no",&number);
  strip_out_float_argument(&argc,argv,"error",&error);

  seq = read_fasta_file_Sequence(argv[1]);

  for(i=0;i<number;i++) {
    eseq= make_ErrorSequence(seq,error,error,error);
    if( just_fasta ) {
      write_fasta_Sequence(eseq->seq,stdout);
    } else {
      show_SequenceErrorSet(eseq->ses,stdout);
      printf("#\n");
      write_fasta_Sequence(eseq->seq,stdout);
      printf("//\n");
    }
    free_ErrorSequence(eseq);
  }
  return 0;
}
