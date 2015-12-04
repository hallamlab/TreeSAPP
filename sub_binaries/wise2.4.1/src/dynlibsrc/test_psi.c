
#include "proteinstreamedindex.h"



int main(int argc,char ** argv)
{
  Sequence * seq;
  ProteinStreamedIndex * psi;
  ProteinStreamedIndexResult * psir;
  Sequence * query;
  int i,j;

  query = read_fasta_file_Sequence(argv[1]);
  fprintf(stderr,"Read in file %s\n",query->name);
  assert(query);

  psi = new_ProteinStreamedIndex(5);

  while( (seq = read_fasta_Sequence(stdin)) != NULL ) {
    fprintf(stderr,"Read in %s... still got %s\n",seq->name,query->name);
    add_Sequence_ProteinStreamedIndex(psi,seq);
  }
  

  for(i=0;i<query->len;i++) {
    psir = lookup_ProteinStreamedIndex(psi,seq_number_aa_5mer(query->seq+i));
    if( psir == NULL ) {
      continue;
    }
    for(j=0;j<psir->len;j++) {
      fprintf(stdout,"\n----------\n");
      fprintf(stdout,"Hit:%d - %d\n%.5s\n%.5s\n",i,psir->result[j].pos,query->seq+i,psir->result[j].seq->seq+psir->result[j].pos);
    }
    free_ProteinStreamedIndexResult(psir);
  }


  /*
  dump_ProteinStreamedIndex(psi,stdout);
  */

  return 0;
}
