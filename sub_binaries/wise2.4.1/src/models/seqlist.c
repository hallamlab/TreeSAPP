#include "dyna.h"



int main(int argc,char ** argv)
{
  int start;
  int end;
  Sequence * in;
  int i;

  in = read_fasta_file_Sequence(argv[1]);
  start = atoi(argv[2]);
  end   = atoi(argv[3]);


  for(i = start;i<= end;i++) {
    fprintf(stdout,"%d  %c   [%c]\n",i,in->seq[i-1],char_complement_base(in->seq[i-1]));
  }

}
