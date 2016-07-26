#include "pwmdna.h"




int main (int argc,char ** argv)
{
  SeqAlign * sa;
  pwmDNA * dna;


  sa = read_selex_SeqAlign(stdin);

  if( sa == NULL ) {
    fprintf(stderr,"Unable to read alignment on stdin");
    exit(1);
  }

  dna = pwmDNA_from_SeqAlign(sa,1.0);

  show_pwmDNA_col(dna,stdout);

}


  
