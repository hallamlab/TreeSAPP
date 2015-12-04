#include "dyna_api.h"


int main(int argc,char ** argv)
{
  Wise2_Sequence * seq;
  Wise2_Sequence * rev;

  if( argc != 2 ) {
    fprintf(stderr,"have to give an argument for a file\n");
    exit(1);
  }

  seq = Wise2_read_fasta_file_Sequence(argv[1]);

  if( seq == NULL ) {
    fprintf(stderr,"Unable to read fasta file in %s\n",argv[1]);
    exit(1);
  }
  
  rev = Wise2_reverse_complement_Sequence(seq);

  printf("Original sequence\n\n");
  Wise2_write_fasta_Sequence(seq,stdout);
  printf("Revcomp sequence\n\n");
  Wise2_write_fasta_Sequence(rev,stdout);
 
  Wise2_free_Sequence(seq);
  Wise2_free_Sequence(rev);
}
