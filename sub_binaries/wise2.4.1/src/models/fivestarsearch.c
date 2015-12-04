#define WISE2_CROSS_HMMER2
#include "wise2xhmmer2.h"


#include "fivestate.h"
#include "version.h"
#include "seqaligndisplay.h"


char * program_name = "fivestarsearch";

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
  fprintf(ofp,"%s fivestar-directory protein-file-fasta\n",program_name);

  fprintf(ofp,"  -ga     gathering cutoff (bits)");
  show_help_DBSearchImpl(ofp);
  
  show_standard_options(ofp);
  
}





int main(int argc,char **argv)
{
  FiveStateFrameSet * frame;
  
  FiveStateModel * fsm;
  FiveStateScore * fss;

  RandomModel * rm;

  ProteinDB * proteindb;
  DBSearchImpl * dbsi;
  Hscore * hs;

  double gathering_cutoff = 0.0;
  double bits;
  int i;

  dbsi = new_DBSearchImpl_from_argv(&argc,argv);

  strip_out_float_argument(&argc,argv,"ga",&gathering_cutoff);

  strip_out_standard_options(&argc,argv,show_help,show_version);
  if( argc != 3 ) {
    show_help(stdout);
    exit(12);
  }

  rm = default_RandomModel();

  frame = read_FiveStateFrameSet_file(argv[1],"block.str");
  if( frame == NULL ) 
    fatal("Unable to make FiveStateModel from context %s, block.str file",argv[1]);
  fsm   = FiveStateModel_from_FiveStateFrameSet(frame);
  /*    dump_FiveStateModel(fsm,stdout); */
  fsm->name = stringalloc(argv[1]);
  
  fold_RandomModel_into_FiveStateModel(fsm,rm);  

  /* converts probabilities to integers for calculation */
  fss = FiveStateScore_from_FiveStateModel(fsm);
  

  proteindb = single_fasta_ProteinDB(argv[2]);

  if( proteindb== NULL )
    fatal("Unable to make proteindb from %s",argv[2]);


  hs = std_score_Hscore(Probability2Score(gathering_cutoff)-10,-1);


  search_FiveStateProtein(dbsi,hs,fss,proteindb);



  fprintf(stdout,"\n\n#High Score list\n");
  fprintf(stdout,"#Protein ID                 DNA Str  ID                        Bits Evalue\n");  
  fprintf(stdout,"--------------------------------------------------------------------------\n");

  for(i=0;i<hs->len;i++) {
    bits = Score2Bits(hs->ds[i]->score);
    if( bits < gathering_cutoff ) {
      break;
    }


    fprintf(stdout,"Protein %-20sDNA [%c] %-24s %.2f\n",hs->ds[i]->query->name,hs->ds[i]->target->is_reversed == TRUE ? '-' : '+',hs->ds[i]->target->name,bits);
  }


}
