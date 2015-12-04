/**
  *
  This file is part of Dynamite
	
  It is copyright (c) 1998 by Ewan Birney and the Sanger Centre.

  The dynamite source code is distributed under a GNU public license.
  Please seee GNULICENSE file in the root directory for this
  package.

  None of this code has any implied use or warranty. Please
  read the license carefully for explanations about re-distribution.

  Please get in contact with Ewan at birney@sanger.ac.uk

  This file was part of the 1.0 release, Tue Jan 20 14:53:11 1998
  *
**/


/*
 * include cdna2genomic.h - will include the dynamite
 * produced declarations provided
 */

#include "cdna2genomic.h" 


/*
 * fancy display 
 */
#include "estgendisplay.h"

void show_help(FILE * ofp)
{
  fprintf(ofp,"\nest2gen <options> est-seq genomic-seq\nBoth sequences in fasta format\n"
	  "\tOPTIONS\n"
	  "\t-g gap penalty (default 2)\n"
	  "\t-e ext penatly (default 1)\n"
	  "\t-m match score (default 4)\n"
	  "\t-n mismatch score (default -3)\n"
	  "\t-r show raw output\n"
	  "\t-l show label output\n"
	  "\t-f show fancy output\n"
	  "\t   (default, -f, all outputs can be shown together\n)"
	  );
}


int main(int argc,char ** argv)
{
  Sequence * query;
  Sequence * target;
  ComplexSequence * query_cs;
  ComplexSequence * target_cs;

  int gap = (2);
  int ext = (1);
  int match = (4);
  int mismatch = (-3);

  Score cdna_open,cdna_ext,gen_open,gen_ext,intron_open;

  DnaMatrix * dm;
  ComplexSequenceEvalSet * cses_g;
  ComplexSequenceEvalSet * cses_c;

  boolean show_raw_output = FALSE;
  boolean show_label_output = FALSE;
  boolean show_fancy_output = FALSE;
  boolean has_outputted = FALSE;

  PackAln * pal;
  AlnBlock * alb;
  
  /*
   * Process command line options
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
  (void) strip_out_integer_argument(&argc,argv,"m",&match);
  (void) strip_out_integer_argument(&argc,argv,"n",&mismatch);

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
   * build dna matrix 
   */

  dm = identity_DnaMatrix(match,mismatch);
  
  
  /*
   * Convert sequences to ComplexSequences: 
   * To do this we need a cdna ComplexSequenceEvalSet and a genomic one
   *
   * Really our genomic model should be alot more complex. The 'default' one
   * has 0 at GT----AG for 5' and 3' splice sites, and NEGI (= -infinity)
   * elsewhere. 
   *
   * We could build up something much better, using complexconsensi and 
   * other machinery, but not for now <grin>. See the genewise code if you
   * want to get scared.
   */
  
  cses_g = default_genomic_ComplexSequenceEvalSet();
  cses_c = default_cDNA_ComplexSequenceEvalSet();
  
  query_cs = new_ComplexSequence(query,cses_c);
  if( query_cs == NULL ) {
    fatal("Unable to make a protein complex sequence from %s",query->name);
  }
  
  target_cs = new_ComplexSequence(target,cses_g);
  if( target_cs == NULL ) {
    fatal("Unable to make a protein complex sequence from %s",target->name);
  }

  free_ComplexSequenceEvalSet(cses_g);
  free_ComplexSequenceEvalSet(cses_c);
  
  /*
   * Make an alignment using best memory
   *
   * cDNA2Gen has alot more parameter space than the parameters to this
   * program. Firstly we are treating errors similarly on each side of the
   * sequences (? correct). 
   *
   * Secondly there is a rather complex interaction between the gap/extension
   * of what is thought to be sequencing error and the introns. Here we have
   * one more parameter, and intron open penalty, which can be set, to prevent
   * the more permissive use of introns to 'cheat' gaps.
   *
   * One good way to parameterise all this would be to have a probabilistic
   * model of the processes, derive probabilities and then map them to ints
   * (probability.h has got these mappings, such as Probability2Score).
   *
   * but this is not done here...
   */		 

  pal = PackAln_bestmemory_cDNA2Gen(query_cs,target_cs,dm,gap,-ext,-gap,-ext,0,NULL);

  if( pal == NULL ) {
    fatal("Unable to make an alignment from %s and %s",query->name,target->name);
  }

  /*
   * ok, make other alignment forms, and be ready to show
   */

  alb = convert_PackAln_to_AlnBlock_cDNA2Gen(pal);


  /*
   * show output. If multiple outputs, divide using //
   */

  if( show_raw_output == TRUE ) {
    show_simple_PackAln(pal,stdout);
    puts("//\n");
  }

  if( show_label_output == TRUE ) {
    show_flat_AlnBlock(alb,stdout);
    puts("//\n");
  }

  if( show_fancy_output == TRUE ) {
    write_pretty_estgen_seq_align(alb,query,target,15,50,stdout);
    puts("//\n");
  }

  /*
   * Destroy the memory.
   */	

  free_Sequence(query);
  free_Sequence(target);
  free_DnaMatrix(dm);
  free_ComplexSequence(query_cs);
  free_ComplexSequence(target_cs);
  free_PackAln(pal);
  free_AlnBlock(alb);

  return 0;
}






