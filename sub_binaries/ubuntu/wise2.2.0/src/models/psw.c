

/*
 * include proteinsw.h - will include the dynamite
 * produced declarations provided
 */

#include "sw_wrap.h" 


#include "abc.h"
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
	  "\t-abc use the abc model\n"
	  "\t-a   a penalty for above (default 120)\n"
	  "\t-b   b penalty for above (default 10)\n"
	  "\t-c   c penalty for above (default 3)\n"
	  "\t-r show raw output\n"
	  "\t-l show label output\n"
	  "\t-f show fancy output\n"
	  "\t   (default, -f, all outputs can be shown together)\n"
	  );
  show_help_DPRunImpl(ofp);
}


int main(int argc,char ** argv)
{
  Sequence * query;
  Sequence * target;
  CompMat * comp;
  char * comp_file;
  int gap = (12);
  int ext = (2);
  int a = 120;
  int b = 10;
  int c = 3;
  ComplexSequence * query_cs;
  ComplexSequence * target_cs;
  ComplexSequenceEvalSet * evalfunc;

  boolean show_label_output = FALSE;
  boolean show_fancy_output = FALSE;
  boolean use_abc = FALSE;

  PackAln * pal;
  AlnBlock * alb;

  DPRunImpl * dpri = NULL;

  /*
   * Process command line options
   * -h or -help gives us help
   * -g for gap value (an int) - rely on commandline error processing
   * -e for ext value (an int) - rely on commandline error processing
   * -m for matrix (a char)
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

  dpri = new_DPRunImpl_from_argv(&argc,argv);
  if( dpri == NULL ) {
    fatal("Unable to build DPRun implementation. Bad arguments");
  }

  show_label_output = strip_out_boolean_argument(&argc,argv,"l");
  show_fancy_output = strip_out_boolean_argument(&argc,argv,"f");


  /** if all FALSE, set fancy to TRUE **/

  if( show_label_output == FALSE ) 
    show_fancy_output = TRUE;


  (void) strip_out_integer_argument(&argc,argv,"g",&gap);
  (void) strip_out_integer_argument(&argc,argv,"e",&ext);
  (void) strip_out_integer_argument(&argc,argv,"a",&a);
  (void) strip_out_integer_argument(&argc,argv,"b",&b);
  (void) strip_out_integer_argument(&argc,argv,"c",&c);

  use_abc = strip_out_boolean_argument(&argc,argv,"abc"); 
  
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
  
  /* if abc - factor up matrix! */

  if( use_abc == TRUE ) {
    factor_CompMat(comp,10);
  }


  /*
   * Make an alignment. I don't care about the implementation:
   * hand it over to sw_wrap function to do it
   *
   */		 

  if( use_abc ) {
    evalfunc = default_aminoacid_ComplexSequenceEvalSet();
  
    query_cs = new_ComplexSequence(query,evalfunc);
    if( query_cs == NULL )
      fatal("Cannot build cs objects!");
    target_cs = new_ComplexSequence(target,evalfunc);
    if( target_cs == NULL )
      fatal("Cannot build cs objects!");

    pal = PackAln_bestmemory_abc(query_cs,target_cs,comp,-a,-b,-c,NULL,dpri);
    alb = convert_PackAln_to_AlnBlock_abc(pal);
    free_ComplexSequence(query_cs);
    free_ComplexSequence(target_cs);
  } else {
    alb = Align_Sequences_ProteinSmithWaterman(query,target,comp,-gap,-ext,dpri);
  }


  /*
   * show output. If multiple outputs, divide using //
   */


  if( show_label_output == TRUE ) {
    show_flat_AlnBlock(alb,stdout);
    puts("//\n");
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
  free_AlnBlock(alb);

  return 0;
}







