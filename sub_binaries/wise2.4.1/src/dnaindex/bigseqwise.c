
#include "comparapath.h"
#include "../models/version.h"
#include "kmer_direct.h"
#include "kmer_hash.h"

char * program_name = "bigseqwise";


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
  fprintf(ofp,"%s reference-assembly-as-fasta query-assembly-as-fasta\n",program_name);
  fprintf(ofp,"  -lognumber [1000000]    number of base pairs for logging [0 for none]\n");
  fprintf(ofp,"  -kmer [12]              kmer size for kmer-graph\n");
  fprintf(ofp,"  -[no]use_hash [default no] use 28bit hash for large kmer sizes\n");
  fprintf(ofp,"  -restrict               restrict reference to only active positions\n");
  show_standard_options(ofp);
}


int main(int argc,char ** argv)
{
  FILE * ref;
  FILE * q;
  FILE * logfp;

  SetofHSPset * hspset;

  ComparaLinkStartSet * reference;
  ComparaLinkStartSet * query;

  ComparaIndex * ci;
  KmerIndexInterface * kii;

  int kmer_size = 12;
  boolean use_28bithash = 0;
  boolean target_restrict = 0;

  int lognumber = 1000000;

  strip_out_integer_argument(&argc,argv,"lognumber",&lognumber);
  strip_out_integer_argument(&argc,argv,"kmer",&kmer_size);
  strip_out_boolean_def_argument(&argc,argv,"use_hash",&use_28bithash);
  strip_out_boolean_def_argument(&argc,argv,"restrict",&target_restrict);
 
  fprintf(stderr,"Hash used is %d\n",use_28bithash);

  strip_out_standard_options(&argc,argv,show_help,show_version);
  if( argc != 4 ) {
    show_help(stdout);
    exit(12);
  }

  ref = openfile(argv[1],"r");
  if( ref == NULL) {
    fatal("Unable to open file %s",argv[1]);
  }

  q   = openfile(argv[2],"r");
  if( ref == NULL) {
    fatal("Unable to open file %s",argv[2]);
  }

  if( strcmp(argv[3],"stderr") == 0 ) {
    logfp = stderr; 
  } else {
    logfp = openfile(argv[3],"w");
    if( logfp == NULL ) {
      fatal("unable to open %s as a logfile",argv[3]);
    }
  }

  fprintf(stderr,"Hash used is %d\n",use_28bithash);


  if( use_28bithash == FALSE)  {
    kii = new_KmerDirectIndex(kmer_size);
  } else {
    fprintf(stderr,"Opening hash index...\n");
    kii = new_KmerHashIndex(kmer_size);
  }

  ci = new_ComparaIndex(kii);

  query     = add_Sequence_stream_ComparaIndex(ci,q,0,lognumber,0,logfp,"query");
  reference = add_Sequence_stream_ComparaIndex(ci,ref,1,lognumber,target_restrict,logfp,"ref");


  insert_revcom_Splines_in_set(ci,query,logfp);

  hspset = SetofHSPset_from_ComparaIndex(ci,query,logfp);

  show_SetofHSPset(hspset,stdout);


  return 0;
}
