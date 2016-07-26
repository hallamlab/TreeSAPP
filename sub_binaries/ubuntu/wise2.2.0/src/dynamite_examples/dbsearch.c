
#include "proteinsw.h"


void show_help (void)
{
  fprintf(stdout,"dbsearch [options] <protein-seq> <protein-fasta-database>\n");
  fprintf(stdout,"Valid options are\n");
  /** add more options here sometime, eg comp matrix and gap penalty*/

  /** print out dbsearch options. We don't know here what implementations are
      either possible or how they are specified. Of course, there is the problem
      that we could clash our options with the dbsearchimpl options, but that
      is not too likely, and this makes this program future proof wrt to new
      implementations
  */

  show_help_DBSearchImpl(stdout);
}

  
int main ( int argc, char ** argv) 
{
  Hscore * out;
  DBSearchImpl * dbsi;
  Protein * temp;
  Sequence * query;
  ProteinDB * querydb;
  ProteinDB * prodb;
  CompMat * mat;

  ComplexSequence * query_cs;
  ComplexSequence * target_cs;
  ComplexSequenceEvalSet  * evalfunc;

  PackAln * pal;
  AlnBlock * alb;

  int i;


  /*
   * processes the command line, removing options
   * that it wants to in order to make the new DBSearchImpl
   *
   * The great thing about this is that this programs does not
   * care about which implementation is used, and does not know either (!)
   *
   */

  dbsi = new_DBSearchImpl_from_argv(&argc,argv);

  if( argc != 3 ) {
    show_help();
    exit(1);
  }

  /*
   * first argument is a single sequence. Read it in and make it
   * into a database
   */

  query = read_fasta_file_Sequence(argv[1]);
  if( query == NULL ) 
    fatal("Cannot read sequence in %s\n",argv[1]);
  
  querydb = new_ProteinDB_from_single_seq(query);


  /*
   * Second argument is a real database. This call is
   * a nice short cut for doing this.
   */

  prodb = single_fasta_ProteinDB(argv[2]);

  if( prodb == NULL )
    fatal("Cannot read protein database in %s\n",argv[2]);


  /*
   * This is where all the results are stored. It also
   * on-the-fly stores distribution information ready
   * to be fitted by a extreme value distribution
   */
  
  /* 10 means a score cutoff of 10, -1 means don't report on stderr search progress */

  out = std_score_Hscore(10,400);

    
  
  mat = read_Blast_file_CompMat("blosum62.bla");
  

  if( search_ProteinSW(dbsi,out,querydb,prodb,mat,-12,-2) != SEARCH_OK ) 
    fatal("Some sort of error in the database search. Dieing ungracefully");

  sort_Hscore_by_score(out);

  evalfunc = default_aminoacid_ComplexSequenceEvalSet();
  
  query_cs = new_ComplexSequence(query,evalfunc);
  if( query_cs == NULL ) {
    fatal("Unable to make a protein complex sequence from %s",query->name);
  }
  

  for(i=0;i<10 && i < out->len;i++) {
      fprintf(stdout,"Comparison to %s was %d score\n",out->ds[i]->target->name,out->ds[i]->score);

      /*
       * Retrieve the protein from the database
       */

      temp = get_Protein_from_ProteinDB(prodb,out->ds[i]->target);

      /*
       * Make a complex sequence of it - see psw for info on this
       */

      target_cs = new_ComplexSequence(temp->baseseq,evalfunc);
      if( target_cs == NULL ) {
	fatal("Unable to make a protein complex sequence from %s",temp->baseseq->name);
      }

      /*
       * Actually align it
       */

      pal = PackAln_bestmemory_ProteinSW(query_cs,target_cs,mat,-12,-2,NULL);
      
      if( pal == NULL ) {
	fatal("Unable to make an alignment from %s and %s",query->name,temp->baseseq->name);
      }
      
      alb = convert_PackAln_to_AlnBlock_ProteinSW(pal);

      write_pretty_seq_align(alb,query,temp->baseseq,15,50,stdout);
      puts("//\n");

      free_Protein(temp);
      free_ComplexSequence(target_cs);
  }
      

  return 0;
}
  



  
  
  
