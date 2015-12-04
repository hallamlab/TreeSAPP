#include "editdistdp.h"
#include "version.h"
#include "standardout.h"
#include "seqaligndisplay.h"

char * program_name = "editdist";

void show_version(FILE * ofp)
{
  fprintf(ofp,"%s\nVersion: %s\nReleased: %s\nCompiled: %s\n",program_name,VERSION_NUMBER,RELEASE_DAY,COMPILE_DATE);
  fprintf(ofp,"\nThis program is freely distributed under a Gnu Public License\n");
  fprintf(ofp,"The source code is copyright (c) EMBL and others\n");
  fprintf(ofp,"There is no warranty, implied or otherwise on the performance of this program\n");
  fprintf(ofp,"For more information read the GNULICENSE file in the distribution\n\n");
  fprintf(ofp,"Credits: Ewan Birney <birney@ebi.ac.uk>\n");
  exit(63);   
}

void show_help(FILE * ofp)
{
  fprintf(ofp,"%s string1 string2\n",program_name);
  fprintf(ofp,"-pretty   show pretty alignment\n");
  show_help_StandardOutputOptions(ofp);

  show_help_DPRunImpl(ofp);

  show_standard_options(ofp);
}



int main (int argc,char ** argv)
{
  EditString * q;
  EditString * t;

  PackAln * pal;
  AlnBlock * alb;

  int match = 1;
  int mismatch = -1;
  int gap = -1;

  boolean show_pretty = FALSE;

  DPRunImpl * dpri = NULL;
  StandardOutputOptions * std_opt = NULL;

  strip_out_boolean_def_argument(&argc,argv,"pretty",&show_pretty);


  dpri = new_DPRunImpl_from_argv(&argc,argv);
  std_opt = new_StandardOutputOptions_from_argv(&argc,argv);

  if( dpri == NULL ) {
    fatal("Unable to build DPRun implementation. Bad arguments");
  }

  if( argc != 3 ) {
    show_help(stdout);
    exit(12);
  }

  q = new_EditString(argv[1]);
  t = new_EditString(argv[2]);

  pal = PackAln_bestmemory_EditDistance(q,t,match,mismatch,gap,NULL,dpri);
  alb = convert_PackAln_to_AlnBlock_EditDistance(pal);

  fprintf(stdout,"Score %d\n",pal->score);
  if( show_pretty == TRUE ) {
    write_pretty_str_align(alb,"string1",q->string,"string2",t->string,12,60,stdout);
  }

  show_StandardOutputOptions(std_opt,alb,pal,"//",stdout);

  return 0;
}
  
  
