#include "linkindex.h"



int main(int argc,char ** argv)
{
  char buffer[MAXLINE];
  Sequence * in;
  DnaNumberSequence * dns;
  LinkNumberArray * lna;
  LinkNumberArrayDebug lnad;
  int i;
  
  lnad.placement_stream = 10;
  lnad.add_stream = 10;
  lnad.extraction = 10;
  lnad.ofp = stderr;

  lna = new_LinkNumberArray(11);

  while( (in = read_fasta_Sequence(stdin)) != NULL ) {
    dns = new_DnaNumberSequence(in,11);

    if( is_new_stream_LinkNumberArray(lna,&lnad,dns) == TRUE ) {
      write_new_LinkStream(lna,&lnad,dns);
    } else {
      write_existing_LinkStream(lna,&lnad,dns);
    }
  }

  for(i=0;i<lna->array_len;i++) {
    if( lna->array[i] != NULL && (lna->array[i]->stream[0]->a == NULL || lna->array[i]->stream[0]->b == NULL) &&  lna->array[i]->stream[0]->have_seen == 0 ) {
      in = extract_dna_LinkStream(lna->array[i]->stream[0],&lnad,lna->nmer_size);
      sprintf(buffer,"Sequence_%d",i);
      in->name = stringalloc(buffer);
      write_fasta_Sequence(in,stdout);
    }
  }


}
