#include "hsp.h"





int main(int argc,char ** argv) 
{
  char * seq1s = "MVNSNQNQNGNSNGHDDDFPQDSITEPEHMRKLFIGGLDYRTTDENLKAHRRKRRR";
  char * seq2s =                  "MHKSEAPNEPEQLRKLFIGGLSFETTDESLREHFEQWGTLTDCVVMRDPNSKRSRGFGFV";
  int query_coord = 33;
  int target_coord = 16;
   
  Sequence * seq1;
  Sequence * seq2;
  HSP * hsp;
  CompMat * mat;


  seq1 = new_Sequence_from_strings("seq1",seq1s);
  seq2 = new_Sequence_from_strings("seq2",seq2s);


  mat = read_Blast_file_CompMat("blosum62.bla");


  hsp = new_HSP(seq1,seq2,query_coord,target_coord,mat,20);

  show_HSP(hsp,80,stdout);

}
