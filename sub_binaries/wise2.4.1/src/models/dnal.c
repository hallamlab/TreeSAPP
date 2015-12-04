
/* quick dna alignment program. Short and sweet */

#include "dnaalign.h"
#include "seqaligndisplay.h"

#include "version.h"

char * program_name = "dnal";

void show_short_help(FILE * ofp);
void show_version(FILE * ofp);

int main(int argc,char ** argv)
{
  Sequence * one;
  Sequence * two;
  DnaMatrix * mat;
  DnaStartEnd * dse;
  AlnBlock * alb;
  DPRunImpl * dpri;

  int kbyte    = 100000;
  int match    = 5;
  int mismatch = -4;
  
  int gap = 5;
  int ext = 1;
  char * policy = NULL;
  
  boolean show_pretty = TRUE;
  boolean show_alb    = FALSE;
  boolean show_pal    = FALSE;
  
  strip_out_standard_options(&argc,argv,show_short_help,show_version);

  dpri = new_DPRunImpl_from_argv(&argc,argv);

  strip_out_boolean_def_argument(&argc,argv,"pretty",&show_pretty); 
  strip_out_boolean_def_argument(&argc,argv,"alb",&show_alb); 
  strip_out_boolean_def_argument(&argc,argv,"pal",&show_pal); 

  strip_out_integer_argument(&argc,argv,"match",&match);
  strip_out_integer_argument(&argc,argv,"mis",&mismatch);
  strip_out_integer_argument(&argc,argv,"gap",&gap);
  strip_out_integer_argument(&argc,argv,"ext",&ext);
  strip_out_integer_argument(&argc,argv,"kbyte",&kbyte);
  
  if( (policy = strip_out_assigned_argument(&argc,argv,"bound")) == NULL )
    policy = "local";

  if( argc != 3 ) {
    show_short_help(stdout);
  }


  one = read_fasta_file_Sequence(argv[1]);
  two = read_fasta_file_Sequence(argv[2]);
  mat = identity_DnaMatrix(match,mismatch);
  dse = DnaStartEnd_from_policy(policy);
  change_max_BaseMatrix_kbytes(kbyte);
  if( one == NULL || two == NULL || mat == NULL || dse == NULL ) {
    fatal("Some objects are null - cannot build alignment");
  }

  alb = make_align_dnaalign(one,two,mat,dse,-gap,-ext,-gap,-ext,dpri);

  printf("Score %d\n",alb->score);

  if( alb == NULL )
    fatal("Could not build alignment!");

  if( show_pretty == TRUE )
    write_pretty_seq_align(alb,one,two,12,60,stdout);

  if( show_alb == TRUE )
    dump_ascii_AlnBlock(alb,stdout);
}


void show_short_help(FILE * ofp)
{
  fprintf(ofp,"%s (%s)\n",program_name,VERSION_NUMBER);
  fprintf(ofp,"This program is freely distributed under a GPL. See -version for more info\n");
  fprintf(ofp,"Copyright (c) GRL limited: portions of the code are from separate copyright\n\n");
  fprintf(ofp,"dnal <dna-file> <dna-file> in fasta format\n");
  fprintf(ofp," Options. In any order, '-' as filename (for any input) means stdin\n");
  fprintf(ofp," -match  [4]  Match score\n");
  fprintf(ofp," -mis    [-1] MisMatch score\n");
  fprintf(ofp," -gap    [5]  Gap penalty\n");
  fprintf(ofp," -ext    [1]  Gap extension penalty\n");
  fprintf(ofp," -bound  [local/global/edge]  Start/End policy\n");
  fprintf(ofp," -kbyte  [100000] Number of kbytes allowed in main memory\n");
  fprintf(ofp," -[no]pretty - show pretty alignment (default - yes)\n");
  fprintf(ofp," -[no]alb    - show alb    alignment (default - no)\n");
  show_help_DPRunImpl(ofp);
  show_standard_options(ofp);
  fprintf(ofp,"\nSee WWW help at http://www.sanger.ac.uk/Software/Wise2/\n");
  exit(63);   
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
