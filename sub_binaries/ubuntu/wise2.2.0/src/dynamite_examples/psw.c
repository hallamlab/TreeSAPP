

/*
 * include proteinsw.h - will include the dynamite
 * produced declarations provided
 */

#include "proteinsw.h" 


/*
 * seqaligndisplay - fancy display 
 *
 */
#include "seqaligndisplay.h"

void show_help(FILE * ofp)
{
  fprintf(ofp,"\npsw <options> seq1 seq2\nBoth sequences in fasta format\n"
	  "\tOPTIONS\n"
	  "\t-g gap penalty (default 12)\n"
	  "\t-e ext penatly (default 2)\n"
	  "\t-m comp matrix (default blosum62.bla)\n"
	  "\t-r show raw output\n"
	  "\t-l show label output\n"
	  "\t-f show fancy output\n"
	  "\t   (default, -f, all outputs can be shown together\n"
	  );
}


int main(int argc,char ** argv)
{
  Sequence * query;
  Sequence * target;
  ComplexSequence * query_cs;
  ComplexSequence * target_cs;
  ComplexSequenceEvalSet  * evalfunc;
  CompMat * comp;
  char * comp_file;
  int gap = (12);
  int ext = (2);

  boolean show_raw_output = FALSE;
  boolean show_label_output = FALSE;
  boolean show_fancy_output = FALSE;
  boolean has_outputted = FALSE;

  PackAln * pal;
  AlnBlock * alb;
  
  /*
   * Process command line options
   * -h or -help gives us help
   * -g for gap value (an int) - rely on commandline error processing
   * -e for ext value (an int) - rely on commandline error processing
   * -m for matrix (a char)
   * -r - raw matrix output
   * -l - label output
   * -f - fancy output
   *
   *
   * Use calls to commandline.h functions
   *
   */
  
  if( strip_out_boolean_argument(&argc,argv,"h") == TRUE || strip_out_boolean_argument(&argc,argv,"-help") == TRUE) {
    show_help(stdout);
    exit(1);
  }

  show_raw_output = strip_out_boolean_argument(&argc,argv,"r");
  show_label_output = strip_out_boolean_argument(&argc,argv,"l");
  show_fancy_output = strip_out_boolean_argument(&argc,argv,"f");


  /** if all FALSE, set fancy to TRUE **/

  if( show_raw_output == FALSE && show_label_output == FALSE ) 
    show_fancy_output = TRUE;


  (void) strip_out_integer_argument(&argc,argv,"g",&gap);
  (void) strip_out_integer_argument(&argc,argv,"e",&ext);

  comp_file = strip_out_assigned_argument(&argc,argv,"m");
  if( comp_file == NULL)
    comp_file = "blosum62.bla";

  
  
  if( argc != 3 ) {
    warn("Must have two arguments for sequence 1 and sequence 2 %d",argc);
    show_help(stdout);
    exit(1);
  }
  
  /*
   * Read in two sequences
   */
  
  if( (query=read_fasta_file_Sequence(argv[1])) == NULL ) {
    fatal("Unable to read the sequence in file %s",argv[1]);
  }
  
  if( (target=read_fasta_file_Sequence(argv[2])) == NULL ) {
    fatal("Unable to read the sequence in file %s",argv[2]);
  }
  
  
  /*
   * Open a blosum matrix. This will be opened from WISECONFIGDIR
   * or WISEPERSONALDIR if it is not present in the current directory.
   */
  
  comp = read_Blast_file_CompMat(comp_file);
  
  if( comp == NULL ) {
    fatal("unable to read file %s",comp_file);
  }
  
  /*
   * Convert sequences to ComplexSequences: 
   * To do this we need an protein ComplexSequenceEvalSet
   *
   */
  
  evalfunc = default_aminoacid_ComplexSequenceEvalSet();
  
  query_cs = new_ComplexSequence(query,evalfunc);
  if( query_cs == NULL ) {
    fatal("Unable to make a protein complex sequence from %s",query->name);
  }
  
  target_cs = new_ComplexSequence(target,evalfunc);
  if( target_cs == NULL ) {
    fatal("Unable to make a protein complex sequence from %s",target->name);
  }
  
  /*
   * Make an alignment. I don't care about the implementation:
   * If the sequences are small enough then it should use explicit memory.
   * Long sequences should use divide and conquor methods.
   *
   * Calling PackAln_bestmemory_ProteinSW is the answer
   * This function decides on the best method considering the
   * memory and changes accordingly. It frees the matrix memory 
   * at the end as well.
   *
   */		 

  pal = PackAln_bestmemory_ProteinSW(query_cs,target_cs,comp,-gap,-ext,NULL);

  if( pal == NULL ) {
    fatal("Unable to make an alignment from %s and %s",query->name,target->name);
  }

  /*
   * ok, make other alignment forms, and be ready to show
   */



  alb = convert_PackAln_to_AlnBlock_ProteinSW(pal);


  /*
   * show output. If multiple outputs, divide using //
   */

  if( show_raw_output == TRUE ) {
    show_simple_PackAln(pal,stdout);
    puts("//\n");
  }

  if( show_label_output == TRUE ) {
    show_flat_AlnBlock(alb,stdout);
  }

  if( show_fancy_output == TRUE ) {
    write_pretty_seq_align(alb,query,target,15,50,stdout);
    puts("//\n");
  }

  /*
   * Destroy the memory.
   */	

  free_Sequence(query);
  free_Sequence(target);
  free_CompMat(comp);
  free_ComplexSequence(query_cs);
  free_ComplexSequence(target_cs);
  free_PackAln(pal);
  free_AlnBlock(alb);

  return 0;
}


