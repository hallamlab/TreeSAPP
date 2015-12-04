#include "shadowseqindex.h"


int main(int argc,char ** argv)
{
  ShadowSequenceIndex * in;
  Sequence * seq;
  int i;

  in = new_ShadowSequenceIndex(26*26*26*26*26,15,0,5000,5);

  i = 0;
  while( (seq = read_fasta_Sequence(stdin)) != NULL ) {
    add_Sequence_ShadowSequenceIndex(in,seq,15);
    if( i % 500  == 0 ) {
      fprintf(stderr,"Loaded %d ...\n",i);
    }
    if( i > 10000 ) {
      break;
    }
    i++;
  }
  
  dump_stats_ShadowSequenceIndex(in,stdout);


  return 0;
}
