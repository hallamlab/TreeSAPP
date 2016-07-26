#include "proteindb.h"


int main(int argc,char ** argv)
{
  ProteinDB * db;
  ComplexSequence * seq;
  int i;
  int status;

  for(i=0;i<1;i++) {
    db = single_fasta_ProteinDB(argv[1]);
    for(seq = init_ProteinDB(db,&status);status != DB_RETURN_END;seq = reload_ProteinDB(seq,db,&status) )
      fprintf(stdout,"Reload %d, Sequence %s\n",i,seq->seq->name);

    close_ProteinDB(seq,db);
    
    fprintf(stdout,"Closed %d\n",i);
    free_ProteinDB(db);
    /*display_allocated_memory("-->",stdout);*/
  }



  return 0;
}
    
  
