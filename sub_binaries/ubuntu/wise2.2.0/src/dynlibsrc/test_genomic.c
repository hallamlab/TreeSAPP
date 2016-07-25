#include "genomic.h"



main(int argc,char ** argv)
{
  Genomic * gen;
  Genomic * gen2;
  Genomic * gen3;
  
  gen = read_fasta_file_Genomic(argv[1],10);

  show_Genomic(gen,stdout);

  gen2 = reverse_complement_Genomic(gen);

  show_Genomic(gen2,stdout);

  gen3 = truncate_Genomic(gen,1,50);

  show_Genomic(gen3,stdout);
  

}
