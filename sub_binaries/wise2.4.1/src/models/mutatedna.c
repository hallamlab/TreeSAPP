#include "dyna.h"




int main(int argc,char ** argv)
{
  Sequence * dna = NULL;
  RandomModelDNA * rmd;
  char * seq;
  Sequence * out;

 
  double subs  = 0.01;
  double indel = 0.002;
  int i,j;

  strip_out_float_argument(&argc,argv,"subs",&subs);
  strip_out_float_argument(&argc,argv,"indel",&indel);

  dna    = read_fasta_file_Sequence(argv[1]);

  if( dna == NULL ) {
    fatal("Must have DNA sequnece as first argument");
  }

  seq = calloc(sizeof(char),dna->len*2);

  rmd = RandomModelDNA_std();

  for(i=0,j=0;i<dna->len;i++) {
    if( random_0_to_1() < subs ) {
      info("Randomising %d\n",i);
      seq[j++] = draw_random_base_RandomModelDNA(rmd);
    } else {
      seq[j++] = dna->seq[i];
    }

    if( random_0_to_1() < indel ) {
      info("inserting at %d\n",i);
      seq[j++] = draw_random_base_RandomModelDNA(rmd);
    }
  }
  seq[j]='\0';

  out = new_Sequence_from_strings(dna->name,seq);
  write_fasta_Sequence(out,stdout);
}
