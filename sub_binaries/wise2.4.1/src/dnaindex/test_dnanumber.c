#include "dnanumber.h"




int main(int argc,char ** argv)
{
  /*
  char * forward  = "GATTGTATTGT";
  char * backward = "ACAATACAATC";
  */
  
  char * forward  = "AGATTGTATTGT";
  char * backward = "GGTTGCAATACAATC";
  
  DnaNumber for_n;
  DnaNumber back_n;

  for_n = dna_number_from_string(forward,11);
  back_n = dna_number_from_string(backward,11);


  printf("Got %s with %d (%d)  %c\n",forward,for_n.number,for_n.flipped,first_char_from_dnanumber(for_n.number,11,for_n.flipped));
  printf("Got %s with %d (%d)  %c\n",backward,back_n.number,back_n.flipped,first_char_from_dnanumber(back_n.number,11,back_n.flipped));


}

  
  
