#include "dyna.h"




int main(int argc,char ** argv)
{
  Sequence * in;
  Sequence * rev;


  while( in = read_fasta_Sequence(stdin) ) {
    in->type = SEQUENCE_DNA;
    rev = reverse_complement_Sequence(in);
    write_fasta_Sequence(rev,stdout);
    free_Sequence(in);
    free_Sequence(rev);
  }


}
