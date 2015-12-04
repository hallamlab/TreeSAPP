#include "model_api.h"


int main(int argc,char ** argv)
{
  Wise2_Sequence * temp;
  Wise2_Protein * one;
  Wise2_Protein * two;
  Wise2_CompMat * comp;
  Wise2_AlnBlock * alb;


  if( argc != 3 ) {
    fprintf(stderr,"sw <file-in-fasta> <file-in-fasta>\n");
    exit(1);
  }

  /** read into generic sequence object */
  temp = Wise2_read_fasta_file_Sequence(argv[1]);

  if( temp == NULL ) {
    fprintf(stderr,"Unable to read fasta file in %s\n",argv[1]);
    exit(1);
  }

  
  /* make into a protein - makes sure the entire sequence looks ok for protein */
  one = Wise2_Protein_from_Sequence(temp);
  
  if( one == NULL ) {
    fprintf(stderr,"Could not make a protein from sequence %s",Wise2_access_name_Sequence(temp));
    exit(1);
  }


  /* the same for the second sequence */
  temp = Wise2_read_fasta_file_Sequence(argv[2]);

  if( temp == NULL ) {
    fprintf(stderr,"Unable to read fasta file in %s\n",argv[2]);
    exit(1);
  }

  two = Wise2_Protein_from_Sequence(temp);
  
  if( two == NULL ) {
    fprintf(stderr,"Could not make a protein from sequence %s",Wise2_access_name_Sequence(temp));
    exit(1);
  }

  /* blosum matrix reading */
  comp = Wise2_read_Blast_file_CompMat("blosum62.bla");

  alb = Wise2_Align_Proteins_SmithWaterman(one,two,comp,-12,-2);

  fprintf(stdout,"Score is %d\n",Wise2_access_score_AlnBlock(alb));

  Wise2_write_pretty_Protein_align(alb,one,two,15,50,stdout);
  
}

  
