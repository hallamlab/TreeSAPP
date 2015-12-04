#include "dnaprofileengine.h"
#include "version.h"


char * program_name = "promotercluster";

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
  fprintf(ofp,"%s fasta-file-for-promoter-set\n",program_name);

  show_help_DnaProfileEnginePara(ofp);

  show_standard_options(ofp);

}


int main(int argc,char ** argv)
{
  DnaProfileEnginePara * dpep;
  DnaProfileNode * root;
  DnaProfileSet * set;
  FILE * ifp;

  boolean is_four = FALSE;

  dpep = new_DnaProfileEnginePara_from_argv(&argc,argv);

  is_four = strip_out_boolean_argument(&argc,argv,"bfour");

  strip_out_standard_options(&argc,argv,show_help,show_version);


  if( argc != 2 ) {
    show_help(stdout);
    exit(12);
  }

  ifp = openfile(argv[1],"r");
  if( ifp == NULL ) {
    fatal("could not open file %s",argv[1]);
  }

  if( is_four ) {
    root = balanced_4_Sequence_fasta_stream(ifp);
  } else {
    root = simple_cascade_Sequence_fasta_stream(ifp);
  }

  populate_DnaProfileNode_from_root(root,dpep);

  /*set = filter_DnaProfileSet(root->set,0,0);*/

  show_DnaProfileSet(root->set,dpep->rm,stdout);

  return 0;
}
