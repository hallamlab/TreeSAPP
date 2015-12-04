#include "localdba.h"
#include "version.h"
#include "hsp.h"
#include "subseqhash.h"

#include "localcishit.h"

#include "hitlist.h"

#include "pairwiseshortdna.h"

char * program_name = "promoterwise";


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
  fprintf(ofp,"%s query_sequence target_sequence\n",program_name);

  fprintf(ofp,"Seed restriction\n");
  fprintf(ofp,"  -align [normal/motif] use normal DBA or motif alignment [normal]\n");
  fprintf(ofp,"  -s    query start position restriction\n");
  fprintf(ofp,"  -t    query end position restriction\n");
  fprintf(ofp,"  -u    target start position restriction\n");
  fprintf(ofp,"  -v    target end position restriction\n");
  show_help_LocalCisHitSetPara(ofp);
  fprintf(ofp,"Motif Matching and TransFactor matches only for motif alignment\n");
  fprintf(ofp,"  ie, when the -align motif option is used\n");
  fprintf(ofp," -lr  motif library is in Laurence's format (default is Ewan's)\n");
  fprintf(ofp," -ben motif library is in Ben's IUPAC format (default is Ewan's)\n");
  fprintf(ofp," -motiflib [filename] motif library file name\n"); 
  show_help_MotifMatrixPara(ofp);
  show_help_TransFactorBuildPara(ofp);
  show_help_TransFactorMatchPara(ofp);
  show_help_HitListOutputImpl(ofp);
  show_help_DPRunImpl(ofp);

  show_standard_options(ofp);

}

#define ALIGN_NORMAL 78
#define ALIGN_MOTIF  79

int main(int argc,char ** argv)
{
  int type = ALIGN_NORMAL;
  DPRunImpl * dpri = NULL;
  HitList * hl;
  HitListOutputImpl * hloi;

  Sequence * query;
  Sequence * target;
  Sequence * target_rev;
  PairwiseShortDna * two;


  LocalCisHitSet * set;
  LocalCisHitSet * greedy_set;

  LocalCisHitScore * lchs;
  LocalCisHitSetPara * setpara;

  MotifMatrixPara  * mmp;
  MotifMatrixScore * mms;

  TransFactorMatchSet * tfms_query = NULL;
  TransFactorMatchSet * tfms_target = NULL;
  TransFactorMatchSet * tfms_target_rev = NULL;

  int qstart = -1;
  int qend   = -1;
  
  int tstart = -1;
  int tend   = -1;
  int i;

  char * temp;

  DnaMatrix * dm;
  DnaProbMatrix * dmp;
  
  TransFactorBuildPara * tfbp;
  TransFactorMatchPara * tfmp;

  TransFactorSet * tfs;

  char * motif_library = NULL;
  int use_laurence     = FALSE;
  int use_ben          = FALSE;

  dmp = DnaProbMatrix_from_match(0.75,NMaskType_BANNED);  
  assert(dmp);
  flat_null_DnaProbMatrix(dmp);  

  dm = DnaMatrix_from_DnaProbMatrix(dmp);
  
  dpri      = new_DPRunImpl_from_argv(&argc,argv);
  hloi      = new_HitListOutputImpl_from_argv(&argc,argv);
  setpara   = new_LocalCisHitSetPara_from_argv(&argc,argv);
  mmp       = new_MotifMatrixPara_from_argv(&argc,argv);
  tfbp      = new_TransFactorBuildPara_from_argv(&argc,argv);
  tfmp      = new_TransFactorMatchPara_from_argv(&argc,argv);

  strip_out_integer_argument(&argc,argv,"s",&qstart);
  strip_out_integer_argument(&argc,argv,"t",&qend);
  strip_out_integer_argument(&argc,argv,"u",&tstart);
  strip_out_integer_argument(&argc,argv,"v",&tend);

  temp = strip_out_assigned_argument(&argc,argv,"motiflib");
  if( temp != NULL ) {
    motif_library = stringalloc(temp);
  }

  use_laurence = strip_out_boolean_argument(&argc,argv,"lr");
  use_ben      = strip_out_boolean_argument(&argc,argv,"ben");


  temp = strip_out_assigned_argument(&argc,argv,"align");
  if( temp != NULL ) {
    if( strcmp(temp,"motif") == 0 ) {
      type = ALIGN_MOTIF;
    } else if ( strcmp(temp,"normal") == 0 ) {
      type = ALIGN_NORMAL;
    } else {
      fatal("cannot recognise string %s as align type",temp);
    }
  }

  strip_out_standard_options(&argc,argv,show_help,show_version);
  if( argc != 3 ) {
    show_help(stdout);
    exit(12);
  }

    

  lchs = standard_LocalCisHitScore(NMaskType_VARIABLE);

  query = read_fasta_file_Sequence(argv[1]);
  target = read_fasta_file_Sequence(argv[2]);

  for(i=0;i<query->len;i++) {
    query->seq[i] = toupper(query->seq[i]);
  }

  assert(query != NULL);
  assert(target != NULL);

  target_rev = reverse_complement_Sequence(target);

  mms = MotifMatrixScore_from_MotifMatrixPara(mmp);

  if( type == ALIGN_MOTIF ) {
    if( motif_library == NULL ) {
      fatal("Wanted to align with motif but not motif library. Must use -motiflib");
    }


    if( use_laurence == TRUE ) {
      tfs = read_laurence_TransFactorSet_file(motif_library);
    } else if( use_ben == TRUE ) {
      tfs = read_ben_IUPAC_TransFactorSet_file(motif_library);
    } else {
      tfs = read_TransFactorSet_file(motif_library);
    }


    build_TransFactorSet(tfs,tfbp);

    tfms_query = calculate_TransFactorMatchSet(query,tfs,tfmp);
    sort_by_start_TransFactorMatchSet(tfms_query);

    tfms_target = calculate_TransFactorMatchSet(target,tfs,tfmp);
    sort_by_start_TransFactorMatchSet(tfms_target);

    tfms_target_rev = calculate_TransFactorMatchSet(target_rev,tfs,tfmp);
    sort_by_start_TransFactorMatchSet(tfms_target);

    fprintf(stdout,"Motif Set: %d in query and %d in target\n",tfms_query->len,tfms_target->len);
  }


  if( qstart == -1 ) {
    qstart = 0;
  }
  if( qend == -1 ) {
    qend = query->len;
  }
  if( tstart == -1 ) {
    tstart = 0;
  }
  if( tend == -1 ) {
    tend = target->len;
  }

  
  two = query_to_reverse_target(query,target,dm,qstart,qend,tstart,tend);

  set = make_LocalCisHitSet(query,target,target_rev,two->forward,two->reverse,setpara,lchs,tfms_query,tfms_target,tfms_target_rev,mms,type == ALIGN_MOTIF ? 1 : 0,dpri);

  greedy_set = greedy_weed_LocalCisHitSet(set,setpara);


  hl = HitList_from_LocalCisHitSet(greedy_set);

  show_HitList_HitListOutputImpl(hloi,hl,stdout);

  return 0;
}



