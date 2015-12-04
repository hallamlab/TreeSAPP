#include "sequencedb.h"
#include "seqlookup.h"
#include "subseqhash.h"
#include "subseqlookup.h"
#include "hsplookupscan.h"
#include "hsptwohitscan.h"
#include "hsplookupthreaded.h"
#include "hspthreadeddb.h"
#include "hspscanruntime.h"
#include "hsp2hitscan.h"

#include "arrayseqlookup.h"
#include "proteinstreamedindex.h"

#include "compressed_protein_index.h"
#include "kmer_direct.h"

#include "../models/version.h"

#include "net_hspscan.h"

#include <sys/time.h>	/* gettimeofday() */

char * program_name = "scanwise_server";

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
  fprintf(ofp,"%s sequence_file_fasta\n",program_name);
  fprintf(ofp,"\nThis program builds a in memory index of a protein database for use by \n");
  fprintf(ofp,"the scanwise family of methods. It runs a wise-style server on the specified port\n");
  fprintf(ofp,"For most databases you need around 10GB of memory\n\n");
  fprintf(ofp," OPTIONS\n");
  fprintf(ofp,"   -port         port to bind to (default 4050)\n");
  fprintf(ofp,"   -compress     use compressed index, for large indexes\n");
  fprintf(ofp,"   -hash         use glib hash not array (far slower, but easier on the memory for small DBs)\n");
  fprintf(ofp,"   -twohit       use two hit seed HSP strategy (can be faster for big databases)\n");
  fprintf(ofp,"   -streamed     use streamed index (more compact, higher run-time)\n");
  fprintf(ofp,"   -ocfile       occupancy file for array/streamed indexes (large), for debugging\n");
  fprintf(ofp,"Search Implementation (default is to provide runtime switch vanilla/threaded/twohit)\n");
  fprintf(ofp,"   -usevanilla   use a vanilla implementation only\n");
  fprintf(ofp,"   -usethreads   use a threaded query scan implementation only\n");
  fprintf(ofp,"   -threadeddb   use a threaded database scan implementation\n");
  fprintf(ofp,"   -ordered      use ordered access implementation\n");
  fprintf(ofp,"   -threadno [2] number of threads for threaded scan implementation\n");
  fprintf(ofp,"   -drop_off     [40] hsp drop off parameter\n");
  fprintf(ofp,"   -array_numb [1000000] hard array numb level in index building\n");

 
  show_help_SeqLookupLoadPara(ofp);
  
  show_standard_options(ofp);
}



int main(int argc, char* argv[])
{
  int i;

  Wise2_HSPScanInterface * hsi;

  SequenceDB * db;
  Sequence * seq;
  int ret;

  SeqLookupInterface * sli;
  CompMat * mat;
  FILE * ofp;

  KmerIndexInterface * kii;

  int count = 0;
  int is_array = 0;
  int waystation = 4;

  int drop_off       = 40;
  int score_cutoff   = 20;

  char * ior_file;
  char * oc_file = NULL;
  char * temp;
  
  SeqLookupLoadPara * slp;

  boolean use_hash = 0;
  boolean use_stream = 0;
  boolean use_comp = 0;
  boolean use_ordered = 0;

  boolean use_threads = 0;
  boolean use_vanilla = 0;
  int threadno = 2;

  int array_numb_level = 1000000;

  boolean use_two_hit = 0;

  boolean use_threadeddb = 0;

  /* threaded database implementations are a bit special */
  HSPThreadedDatabase * tdb;

  FunctionServer * fs;
  FunctionImplementation * impl;
  int port = 4050;

  struct timeval t0, t1;

  gettimeofday(&t0, NULL);

  slp = new_SeqLookupLoadPara_from_argv(&argc,argv);

  use_hash   = strip_out_boolean_argument(&argc,argv,"hash");
  use_stream = strip_out_boolean_argument(&argc,argv,"streamed");
  use_comp   = strip_out_boolean_argument(&argc,argv,"compress");

  use_threadeddb = strip_out_boolean_argument(&argc,argv,"threadeddb");
  strip_out_integer_argument(&argc,argv,"array_numb",&array_numb_level);

  use_threads = strip_out_boolean_argument(&argc,argv,"usethreads");
  use_vanilla = strip_out_boolean_argument(&argc,argv,"usevanilla");

  strip_out_integer_argument(&argc,argv,"threadno",&threadno);

  strip_out_integer_argument(&argc,argv,"drop_off",&drop_off);

  strip_out_integer_argument(&argc,argv,"hspext_cutoff",&score_cutoff);


  use_ordered = strip_out_boolean_argument(&argc,argv,"ordered");

  use_two_hit = strip_out_boolean_argument(&argc,argv,"twohit");

  strip_out_integer_argument(&argc,argv,"waystation",&waystation);


  oc_file = strip_out_assigned_argument(&argc,argv,"ocfile");

  strip_out_integer_argument(&argc,argv,"port",&port);

  strip_out_standard_options(&argc,argv,show_help,show_version);

  if( argc != 2 ) {
    show_help(stdout);
    exit(12);
  }


  fs = new_FunctionServer(port);

  db = single_fasta_SequenceDB(argv[1]);

  mat = read_Blast_file_CompMat("blosum62.bla");


  if( use_stream == 1 ) {
    sli = new_ProteinStreamedIndex_SeqLookupInterface(waystation);
  } else if( use_hash == 1 ) {
    sli = new_ghash_SeqLookupInterface();
  } else if ( use_comp == 1 ) {
    /* this is a little weird */
    /* the kmer index system was originally designed for DNA */
    /* 12mers in kmer space is about the same as 5mers in aa space */

    kii = new_KmerDirectIndex(12);
    sli = new_CompressedProteinLookup(kii);
      
  } else if ( use_threadeddb == 1 ) {
    /* do nothing here, as we have to do this different */
  } else {
    is_array = 1;
    sli = new_ArraySeq_SeqLookupInterface(26*26*26*26*26,array_numb_level); 
  }


  if( use_threadeddb != 1 ) {
    load_SequenceDB_SeqLookupLoadPara(slp,db,sli);

    if( (is_array == 1 || use_stream == 1) && oc_file != NULL) {
      ofp = fopen(oc_file,"w");
      if( ofp != NULL ) {
	info("Printing index occupancy");
	if( is_array == 1 ) {
	  print_array_occuypancy_ArraySeq((ArraySeqLookup*)sli->data,ofp);
	} else {
	  dump_ProteinStreamedIndex((ProteinStreamedIndex*)sli->data,ofp);
	}
	info("Finished index occupancy");
	fclose(ofp);
      }
    }
    
    if( use_ordered ) {
      hsi = new_ordered_HSPScanInterface(sli,mat,40,10);
    } else if( use_threads ) {
      hsi = new_threaded_HSPScanInterface(sli,mat,40,10,threadno);
    } else if( use_two_hit ) {
      hsi = new_twohit_one_off_HSPScanInterface(sli,mat,drop_off,score_cutoff);
    } else if ( use_vanilla ) {
      hsi = Wise2_new_one_off_HSPScanInterface(sli,mat,drop_off,score_cutoff);    
    } else {
      /* provide runtime implementation */
      hsi = new_runtime_HSPScanInterface(sli,mat,drop_off,score_cutoff,threadno);
    }
  } else {
    tdb = new_HSPThreadedDatabase(threadno,array_numb_level);
    load_HSPThreadedDatabase(tdb,db,slp,mat,drop_off,score_cutoff);
    hsi = new_HSPScanInterface_from_HSPThreadedDatabase(tdb);
  }
    

  impl = new_hspscan_protein_FunctionImplementation(hsi);

  add_FunctionServer(fs,impl);

  gettimeofday(&t1, NULL);
  
  info("[server stats] startup time (s): %f",
       (t1.tv_sec - t0.tv_sec) +
       (t1.tv_usec - t0.tv_usec) * 1e-6);
  

  info("Started server on port %d",port);

  main_loop_forking_FunctionServer(fs,1);

  return 0;
}
