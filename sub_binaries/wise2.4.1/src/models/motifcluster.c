#include "transfactor.h"
#include "transregion.h"
#include "version.h"
#include "commandline.h"

char * program_name = "motifcluster";

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
  fprintf(ofp,"%s motif-library sequence\n",program_name);
  fprintf(ofp," -lr  motif library is in laurence's format\n");
  fprintf(ofp," -ben motif library is in ben's IUPAC format\n");
  fprintf(ofp," -[no]show_match  - show all matches that pass criteria, default no\n");
  fprintf(ofp," -[no]show_region - show only dense cluster regions, default no\n");
  fprintf(ofp," -[no]show_compara - show only comparative matches, default yes\n");
  fprintf(ofp," -circular [no]   - for randomisation, number of positions to permute\n");

  show_help_TransFactorBuildPara(ofp);

  show_help_TransFactorMatchPara(ofp);

  show_help_TransFactorRegionPara(ofp);

  show_help_TransFactorComparaPara(ofp);

  show_help_DPRunImpl(ofp);

  show_standard_options(ofp);
}



int main(int argc,char ** argv)
{
  TransFactorSet * tfs;
  TransFactorSet * newtfs; /* only used as a temp when permuting */
  TransFactorBuildPara * tfb;
  TransFactorMatchSetCompara * comp;
  TransFactorMatchPara * matchp;
  TransFactorComparaPara * comparap;

  boolean use_laurence = FALSE;
  boolean use_ben = FALSE;
  boolean show_matchset = FALSE;
  boolean show_region  = FALSE;
  boolean show_compara = TRUE;

  int end_on_seq = 0;
  int seq_count = 0;

  int rotate_number = 0;

  TransFactorRegionSet * tfrs;
  TransFactorRegionPara * tfrp;

  DPRunImpl * dpri;

  FILE * ifp;
  SeqAlign * sa;

  int i;

  tfb    = new_TransFactorBuildPara_from_argv(&argc,argv);
  matchp = new_TransFactorMatchPara_from_argv(&argc,argv);
  tfrp   = new_TransFactorRegionPara_from_argv(&argc,argv);
  comparap = new_TransFactorComparaPara_from_argv(&argc,argv);
  dpri   = new_DPRunImpl_from_argv(&argc,argv);

  strip_out_boolean_def_argument(&argc,argv,"show_match",&show_matchset);
  strip_out_boolean_def_argument(&argc,argv,"show_region",&show_region);
  strip_out_integer_argument(&argc,argv,"circular",&rotate_number);


  if( strip_out_boolean_argument(&argc,argv,"lr") == TRUE ) {
    use_laurence = TRUE;
  }
  if( strip_out_boolean_argument(&argc,argv,"ben") == TRUE ) {
    use_ben = TRUE;
  }


  strip_out_standard_options(&argc,argv,show_help,show_version);


  if( argc != 3 ) {
    show_help(stdout);
    exit(12);
  }

  if( use_laurence == TRUE ) {
    tfs = read_laurence_TransFactorSet_file(argv[1]);
  } else if( use_ben == TRUE ) {
    tfs = read_ben_IUPAC_TransFactorSet_file(argv[1]);
  } else {
    tfs = read_TransFactorSet_file(argv[1]);
  }

  if( tfs->len == 0 ) {
    fatal("No transcription factors in set!");
  }

  build_TransFactorSet(tfs,tfb);

  if( rotate_number > 0 ) {
    fprintf(stdout,"PERMUTED RESULTS - %d rotation\n",rotate_number);
    newtfs = circular_permuted_TransFactorSet(tfs,rotate_number);
    free_TransFactorSet(tfs);
    tfs = newtfs;
  }

  ifp = openfile(argv[2],"r");
  if( ifp == NULL ) {
    fatal("Cannot open %s as input file",argv[2]);
  }

  while( (sa = read_fasta_SeqAlign(ifp)) ) {
    for(i=0;i<sa->len;i++) {
      sa->seq[i]->type = SEQUENCE_DNA;
    }
    
    comp = calculate_TransFactorMatchSetCompara(sa,tfs,matchp,comparap);
    
    sort_by_start_TransFactorMatchSet(comp->overall);
    
    
    fprintf(stdout,"%s %d raw hits %d conserved\n",comp->overall->target->name,comp->tfms[0]->len,comp->overall->len);
    
    fflush(stdout);
    
    if( show_compara == TRUE ) {
      show_context_match_TransFactorMatchSetCompara(comp,15,stdout);
    }
    
    
    if( show_matchset == TRUE ) {
      show_TransFactorMatchSet(comp->overall,stdout);
    }
    
    
    if( show_region == TRUE ) {
      tfrs = new_TransFactorRegionSet(comp->overall,tfrp,dpri);
      show_TransFactorRegionSet(tfrs,stdout);
      free_TransFactorRegionSet(tfrs);
    }
    
    
    free_TransFactorMatchSetCompara(comp);
    free_SeqAlign(sa);
    fprintf(stdout,"//\n");
  }


}
