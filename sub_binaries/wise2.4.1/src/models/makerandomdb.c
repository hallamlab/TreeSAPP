#include "randomdb.h"
#include "wisebase.h"

int main(int argc,char ** argv)
{
  int i;
  RandomDNADB * db;
  RandomModelDNA * rndmodel;
  int len = 1000;
  Sequence * seq;
  int num = 2000;

  strip_out_integer_argument(&argc,argv,"l",&len);

  rndmodel = RandomModelDNA_std();
  db = new_flat_RandomDNADB(rndmodel,len);

  for(i=0;i<num;i++) {
    seq = Sequence_from_RandomDNADB(db);
    write_fasta_Sequence(seq,stdout);
  }

}

