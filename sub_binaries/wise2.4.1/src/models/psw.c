

/*
 * include proteinsw.h - will include the dynamite
 * produced declarations provided
 */

#include "sw_wrap.h" 
#include "commandline.h"
#include "version.h"

#include "abc.h"
/*
 * seqaligndisplay - fancy display 
 *
 */
#include "seqaligndisplay.h"

char * program_name = "psw";

void show_help(FILE * ofp)
{
  fprintf(ofp,"\npsw <options> seq1 seq2\nBoth sequences in fasta format\n"
	  "\tOPTIONS\n"
	  "\t-g gap penalty (default 12)\n"
	  "\t-e ext penatly (default 2)\n"
	  "\t-m comp matrix (default BLOSUM62.bla)\n"
	  "\t-abc use the abc model\n"
	  "\t-a   a penalty for above (default 120)\n"
	  "\t-b   b penalty for above (default 10)\n"
	  "\t-c   c penalty for above (default 3)\n"
	  "\t-r show raw output\n"
	  "\t-l show label output\n"
	  "\t-f show fancy output\n"
          "\t-F force psw to use sequences that seem to be DNA\n"
	  "\t-dpenv DP envelope file...\n"
	  "\t   (default, -f, all outputs can be shown together)\n"
	  );
  show_help_DPRunImpl(ofp);
  show_standard_options(ofp);
  exit(63);
}

void show_version(FILE * ofp);

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

  int ident;
  int qstart;
  int qend;
  int tstart;
  int tend;
  AlnColumn * alc;

  ComplexSequence * query_cs;
  ComplexSequence * target_cs;
  ComplexSequenceEvalSet * evalfunc;

  boolean show_label_output = FALSE;
  boolean show_fancy_output = FALSE;
  boolean use_abc = FALSE;
  boolean show_perc = FALSE;
  boolean force_protein = FALSE;

  PackAln * pal;
  AlnBlock * alb;

  DPRunImpl * dpri = NULL;
  DPEnvelope * dpenv = NULL;
  char * temp;

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
  
  strip_out_standard_options(&argc,argv,show_help,show_version);

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
  show_perc = strip_out_boolean_argument(&argc,argv,"p");
  force_protein = strip_out_boolean_argument(&argc,argv,"F");

  /** if all FALSE, set fancy to TRUE **/

  if( show_label_output == FALSE && show_perc == FALSE) 
    show_fancy_output = TRUE;


  (void) strip_out_integer_argument(&argc,argv,"g",&gap);
  (void) strip_out_integer_argument(&argc,argv,"e",&ext);
  (void) strip_out_integer_argument(&argc,argv,"a",&a);
  (void) strip_out_integer_argument(&argc,argv,"b",&b);
  (void) strip_out_integer_argument(&argc,argv,"c",&c);

  use_abc = strip_out_boolean_argument(&argc,argv,"abc"); 
  
  comp_file = strip_out_assigned_argument(&argc,argv,"m");
  if( comp_file == NULL)
    comp_file = "BLOSUM62.bla";

  if( (temp = strip_out_assigned_argument(&argc,argv,"dpenv")) != NULL ) {
    dpenv = read_DPEnvelope_file(temp);
  }
  
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

  if ( force_protein ) {
    query->type = SEQUENCE_PROTEIN;
    target->type = SEQUENCE_PROTEIN;
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



/*  show_DPEnvelope(dpenv,stderr); */


  if( use_abc ) {
    evalfunc = default_aminoacid_ComplexSequenceEvalSet();
  
    query_cs = new_ComplexSequence(query,evalfunc);
    if( query_cs == NULL )
      fatal("Cannot build cs objects!");
    target_cs = new_ComplexSequence(target,evalfunc);
    if( target_cs == NULL )
      fatal("Cannot build cs objects!");

    pal = PackAln_bestmemory_abc(query_cs,target_cs,comp,-a,-b,-c,dpenv,dpri);
    alb = convert_PackAln_to_AlnBlock_abc(pal);
    free_ComplexSequence(query_cs);
    free_ComplexSequence(target_cs);
  } else {
    alb = Align_Sequences_ProteinSmithWaterman(query,target,comp,-gap,-ext,dpenv,dpri);
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

  if( show_perc == TRUE ) {
    qstart = alb->start->alu[0]->end;
    tstart = alb->start->alu[1]->end;
    ident =0;

    for(alc = alb->start;alc != NULL;alc = alc->next) {
      if( strcmp(alc->alu[0]->text_label,"SEQUENCE") == 0 &&
	  strcmp(alc->alu[1]->text_label,"SEQUENCE") == 0 ) {

	qend = alc->alu[0]->end;
	tend = alc->alu[1]->end;
	if( query->seq[alc->alu[0]->end] == target->seq[alc->alu[1]->end] ) {
	  ident++;
	}
	

      }
    }

    fprintf(stdout,"%d %s %d:%d (%.2f) %s %d:%d (%.2f)\n",ident,
	  query->name,qstart+1,qend,(ident*100.0)/(qend-qstart),
	   target->name,tstart+1,tend,(ident*100.0)/(qend-qstart));

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

void show_version(FILE * ofp)
{
  fprintf(ofp,"%s\nVersion: %s\nReleased: %s\nCompiled: %s\n",program_name,VERSION_NUMBER,RELEASE_DAY,COMPILE_DATE);
  fprintf(ofp,"\nThis program is freely distributed under a Gnu Public License\n");
  fprintf(ofp,"The source code is copyright (c) GRL 1998 and others\n");
  fprintf(ofp,"There is no warranty, implied or otherwise on the performance of this program\n");
  fprintf(ofp,"For more information read the GNULICENSE file in the distribution\n\n");
  fprintf(ofp,"Credits: Ewan Birney <birney@sanger.ac.uk> wrote the core code.\n");
  exit(63);   
}
