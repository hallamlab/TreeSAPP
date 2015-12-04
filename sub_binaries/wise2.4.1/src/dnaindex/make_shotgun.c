#include "shotgun.h"




int main(int argc,char ** argv)
{
  Sequence * input;
  ShotgunPara p;

  p.read_length = 500;
  p.insert_size = 200;
  p.number = 4000;
  p.forward_only = 1;

  input = read_fasta_Sequence(stdin);

  generate_shotgun_reads(&p,input,stdout);

}
