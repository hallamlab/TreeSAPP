#include "sequencedb.h"


int main(int argc,char ** argv)
{
  SequenceDB * db;
  Sequence * seq;
  int i;
  int status;

  for(i=0;i<4;i++) {
    db = single_fasta_SequenceDB(argv[1]);
    for(seq = init_SequenceDB(db,&status);status != DB_RETURN_END;seq = reload_SequenceDB(seq,db,&status) )
      fprintf(stdout,"Reload %d, Sequence %s\n",i,seq->name);

    close_SequenceDB(seq,db);
    fprintf(stdout,"Closed %d\n",i);
  }



  return 0;
}
    
  
