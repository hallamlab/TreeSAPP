#include <stdio.h>
#include "sw_wrap.h"
#include "seqlookup.h"
#include "subseqhash.h"
#include "subseqlookup.h"
#include "sequencedb.h"
#include "hsp.h"
#include "hitlist.h"
#include "hsplookupscan.h"
#include "version.h"
#include "hsp2aln_sw.h"
#include "proteinindexcons.h"
#include "hspscaninterface.h"
#include "client_multihspscan.h"

#include <sys/time.h>
#include <sys/resource.h>

#include "wise2_mott_bridge.h"

#ifdef SCAN_CORBA
#include "hspscan_corba_wrapper.h"
#include "corba_singleton.h"
#endif

#ifdef SCAN_MYSQL
#include "mysql_protein_index.h"
#endif

#ifdef SCAN_WISESERVER
#include "net_hspscan.h"
#endif

#ifdef SCAN_COMPRESS
#include "compressed_protein_index.h"
#endif


#define MICROSECOND (1.0/1000000.0)

char * program_name = "scanwisep";

typedef struct {
  boolean use_corba;
  boolean use_mysql;
  boolean use_wiseserver;
  boolean use_compress;
  boolean use_multiscan;
  char * ior_file;
  char * direct_sequence;
  char * matrix_file;
  char * host;
  char * dbname;
  char * username;
  char * password;
  char * multiscan_file;
  int step;
  int port;
} ScanWiseHSPImpl;

#ifdef SCAN_CORBA
  Wise2Corba_Singleton * sorb;
#endif

HSPScanInterface * new_HSPScanInterface_from_ScanWiseHSPImpl(ScanWiseHSPImpl * i,ProteinIndexConstructor * pic,SeqLookupLoadPara * slp)
{
  HSPScanInterface * out;
  SeqLookupInterface * sli;
  SequenceDB * db;
  Sequence * seq;
  CompMat * mat;
  int ret;
  int c;
  
  mat = read_Blast_file_CompMat(i->matrix_file);
 
  if( i->use_corba == FALSE && i->use_mysql == FALSE && i->use_wiseserver == FALSE && i->use_compress == FALSE && i->use_multiscan == FALSE) {
    if( i->direct_sequence == NULL ) {
      fatal("If no server based sequence, must have direct sequence");
    } else {
      db = single_fasta_SequenceDB(i->direct_sequence);

      sli = new_SeqLookupInterface_from_ProteinIndexConstructor(pic);

      load_SequenceDB_SeqLookupLoadPara(slp,db,sli);

      free_SequenceDB(db);

      out = new_one_off_HSPScanInterface(sli,mat,15,40);
    }
  } else if( i->use_corba == TRUE ) {
#ifdef SCAN_CORBA 
    if( i->ior_file == NULL ) {
      fatal("Corba specified, but no ior file given");
    }


    out = new_corba_HSPScan(sorb,i->ior_file,mat);
#else
    fatal("Asking for CORBA, but scanwisep was not compiled with SCAN_CORBA defined.");
#endif
  } else if ( i->use_mysql == TRUE ) {
#ifdef SCAN_MYSQL
    out = new_HSPScanInterface_MysqlProteinIndex(i->host,i->dbname,i->username,i->password,mat,i->step);
#else
    fatal("Asking for mysql, but scanwisep was not compiled with SCAN_MYSQL defined");
#endif
  } else if ( i->use_wiseserver == TRUE ) {

#ifdef SCAN_WISESERVER
    out = new_wise_transfer_HSPScanInterface(i->host,i->port);
#else
    fatal("Asking for wiseserver, but scanwisep was not compiled with SCAN_WISESERVER defined");
#endif
  } else if ( i->use_compress == TRUE ) {

#ifdef SCAN_COMPRESS
    sli = new_direct_CompressedProteinLookup();
    db = single_fasta_SequenceDB(i->direct_sequence);

    load_SequenceDB_SeqLookupLoadPara(slp,db,sli);

    free_SequenceDB(db);

    out = new_one_off_HSPScanInterface(sli,mat,15,40);

#else
    fatal("Asking for compressed, but scanwisep was not compiled with SCAN_COMPRESS defined");
#endif
  } else if( i->use_multiscan == TRUE ) {
    if( i->multiscan_file == NULL ) {
      fatal("Must provide a file for a multiple server scan");
    }
    
    out = new_multiclient_HSPScanInterface(i->multiscan_file);
  }



  assert(out != NULL);


  free_CompMat(mat); /* hard linked internally */
  return out;
}


ScanWiseHSPImpl * new_ScanWiseHSPImpl_from_argv(int * argc,char ** argv)
{
  ScanWiseHSPImpl * out;
  char * temp;

  out = malloc(sizeof(ScanWiseHSPImpl));
  out->use_corba = FALSE;
  out->use_mysql = FALSE;
  out->use_compress = FALSE;
  out->use_multiscan = FALSE;
  out->ior_file = NULL;
  out->direct_sequence = NULL;
  out->matrix_file = "BLOSUM62.bla";
  out->step = 32;
  out->host = "localhost";
  out->port = 4050;
  

  strip_out_boolean_def_argument(argc,argv,"corba",&out->use_corba);

  strip_out_boolean_def_argument(argc,argv,"mysql",&out->use_mysql);

  strip_out_boolean_def_argument(argc,argv,"wiseserver",&out->use_wiseserver);

  strip_out_boolean_def_argument(argc,argv,"compress",&out->use_compress);

  strip_out_boolean_def_argument(argc,argv,"multi",&out->use_multiscan);
  
  if( (temp = strip_out_assigned_argument(argc,argv,"iorfile")) != NULL ) {
    out->ior_file = temp;
  }

  if( (temp = strip_out_assigned_argument(argc,argv,"multiserver")) != NULL ) {
    out->multiscan_file = temp;
  }

  if( (temp = strip_out_assigned_argument(argc,argv,"scan_host")) != NULL ) {
    out->host = temp;
  }
  strip_out_integer_argument(argc,argv,"scan_port",&out->port);

  if( (temp = strip_out_assigned_argument(argc,argv,"scan_dbname")) != NULL ) {
    out->dbname = temp;
  }
  if( (temp = strip_out_assigned_argument(argc,argv,"scan_username")) != NULL ) {
    out->username = temp;
  }
  if( (temp = strip_out_assigned_argument(argc,argv,"scan_password")) != NULL ) {
    out->password = temp;
  }

  strip_out_integer_argument(argc,argv,"scan_step",&out->step);

  if( (temp = strip_out_assigned_argument(argc,argv,"seqdb")) != NULL ) {
    out->direct_sequence = temp;
  }

  if( (temp = strip_out_assigned_argument(argc,argv,"seqdbmat")) != NULL ) {
    out->matrix_file = temp;
  }


  return out;

}

void show_help_ScanWiseHSPImpl(FILE * ofp)
{
  fprintf(ofp,"HSP generation options\n");

#ifdef SCAN_CORBA
  fprintf(ofp,"  -[no]corba  Use CORBA server for index (CORBA enabled)\n");
  fprintf(ofp,"  -iorfile    For CORBA cases, IOR file to use for CORBA access\n");
#else
  fprintf(ofp,"  -[no]corba  Use CORBA server for index (CORBA not enabled. go make scanwisep_corba)\n");
#endif

#ifdef SCAN_MYSQL
  fprintf(ofp,"  -[no]mysql  Use mysql server for index (Mysql enabled)\n");
  fprintf(ofp,"  -scan_host [localhost] host name for mysql server\n");
  fprintf(ofp,"  -scan_username []      user name for mysql server\n");
  fprintf(ofp,"  -scan_password         password for mysql server\n");
  fprintf(ofp,"  -scan_dbname           database name for mysql server\n");
  fprintf(ofp,"  -scan_step [32]        seq step query block\n");
#else
  fprintf(ofp,"  -[no]mysql  Use mysql server for index (Mysql not enabled. go make scanwisep_mysql)\n");
#endif

#ifdef SCAN_WISESERVER
  fprintf(ofp,"  -[no]wiseserver  Use wise socket based server for index (enabled)\n");
  fprintf(ofp,"  -scan_host [localhost] host name for wise server\n");
  fprintf(ofp,"  -scan_port [4050]      port for wise server\n");
  fprintf(ofp,"  -[no]multi       Use multiple wiseservers at once\n");
  fprintf(ofp,"  -multiserver <filename> Filename for multiple servers, <host> <port> format\n");
#else
  fprintf(ofp,"  -[no]wiseserver  Use wise server for index (Wise server. go make scanwisep_wiseserver)\n");
#endif


  fprintf(ofp,"  -seqdb      For local cases, sequence database fasta file\n");
  fprintf(ofp,"  -seqdbmat   For local cases, comparison matrix to use\n");

}



/*
 * conversion from HSPs to HitList
 */

typedef enum HSP2HitListType {
  HSP2HitList_Ungapped = 445,
  HSP2HitList_Fulldp,
  HSP2HitList_Heuristic,
  HSP2HitList_None,
  HSP2HitList_Unknown
} HSP2HitListType;

typedef struct HSP2HitListImpl {
  HSP2HitListType type;
  boolean threaded;
  int no_threads;
} HSP2HitListImpl;


HSP2HitListType HSP2HitListType_string_convert(char * string)
{
  if( strcmp(string,"ungapped") == 0 ) {
    return HSP2HitList_Ungapped;
  } 

  if( strcmp(string,"fulldp") == 0 ) {
    return HSP2HitList_Fulldp;
  } 

  if( strcmp(string,"heuristic") == 0 ) {
    return HSP2HitList_Heuristic;
  } 

  if( strcmp(string,"none") == 0 ) {
    return HSP2HitList_None;
  } 

  return HSP2HitList_Unknown;
}


struct convert_data_fulldp {
  HitList * hl;
  DPRunImpl * dpri;
  CompMat * mat;
};

HitPair * HitPair_from_HSPset_sw(HSPset * set,DPRunImpl * dpri,CompMat * mat)
{
  HitPair * out;
  HitAln * aln;
  int i;

  out = HitPair_alloc_std();
  out->query  = hard_link_Sequence(set->hsp[0]->query);
  out->target = hard_link_Sequence(set->hsp[0]->target);
  aln = HitAln_alloc();

  aln->alb = Align_Sequences_ProteinSmithWaterman(out->query,out->target,mat,-12,-2,NULL,dpri);
  aln->raw_score = aln->alb->score;
  add_HitPair(out,aln);

  out->raw_score = aln->alb->score;

  return out;
}


HitList * HitList_from_LinearHSPmanager_sw(LinearHSPmanager * lm,DPRunImpl * dpri)
{
  HitList * out;
  int i;
  struct convert_data_fulldp conv;

  
  out = HitList_alloc_std();
  out->mat = hard_link_CompMat(lm->mat);
 
  for(i=0;i<lm->len;i++) 
    add_HitList(out,HitPair_from_HSPset_sw(lm->set[i],dpri,out->mat));
		

  return out;

}


HSP2HitListImpl * new_HSP2HitListImpl_from_argv(int * argc,char ** argv)
{
  HSP2HitListImpl * out;
  char * temp;

  out = malloc(sizeof(HSP2HitListImpl));
  out->type = HSP2HitList_Heuristic;

  if( (temp = strip_out_assigned_argument(argc,argv,"hspconvert")) != NULL ) {
    out->type = HSP2HitListType_string_convert(temp);
  }

  strip_out_boolean_def_argument(argc,argv,"hspthread",&out->threaded);

  out->no_threads = 4;
  strip_out_integer_argument(argc,argv,"hspthreadno",&out->no_threads);

  return out;

}

void show_help_HSP2HitList(FILE * ofp)
{
  fprintf(ofp,"Conversion from HSP to alignments\n");
  fprintf(ofp,"  -hspconvert    [ungapped/fulldp/heuristic] Conversion type - heuristic default\n");
  fprintf(ofp,"  -hspthread     multi-thread HSP conversion\n");
  fprintf(ofp,"  -hspthreadno   number of HSP threads (4 default)\n");
}


HitList * HitList_from_HSP_HSP2HitListImpl(HSP2HitListImpl * conv,LinearHSPmanager * lm,DPRunImpl * dpri,HSPset2HitPairPara * hsp2hit)
{
  HitList * processed = NULL;
  HitList * temp;
  int i;
  int topscore;


  assert(lm);
  assert(dpri);

  if( conv->type == HSP2HitList_Ungapped ) {
    processed = HitList_from_LinearHSPmanager(lm);
  } else if ( conv->type == HSP2HitList_None ) {
    processed = HitList_alloc_std();
    return processed;
  } else if( conv->type == HSP2HitList_Fulldp ) {
    processed = HitList_from_LinearHSPmanager_sw(lm,dpri);
  } else if( conv->type == HSP2HitList_Heuristic ) {
    if( conv->threaded == FALSE ) {
      processed = HitList_from_LinearHSPmanager_heuristic(lm,dpri,hsp2hit);
    } else  {
      processed = HitList_from_LinearHSPmanager_heuristic_threaded(lm,dpri,conv->no_threads,hsp2hit);
    }
  } else {
    fatal("Could not covert type %d for hit list conversion",conv->type);
  }

  /* if there is a best-in-genome step, do a final screen of the results */

  if( hsp2hit->best_hit == TRUE && processed->len > 2) {
    temp = HitList_alloc_len(processed->len);
    topscore = processed->pair[0]->raw_score;
    for(i=0;i<processed->len;i++) {
      if( ((double)(processed->pair[i]->raw_score - topscore))*100.0/topscore > hsp2hit->perc_hit_dropoff ) {
	break;
      } else {
	add_HitList(temp,hard_link_HitPair(processed->pair[i]));
      }
    }
    free_HitList(processed);
    processed= temp;
  }


  return processed;

}


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
  fprintf(ofp,"%s query_sequence_file_fasta\n",program_name);

  fprintf(ofp,"   -dbsize [number] effective db size for Evalue calculation [300000]\n");

  fprintf(ofp,"   -[no]mott use Mott's statistics or not (default yes)\n");

  show_help_ScanWiseHSPImpl(ofp);

  show_help_HSPScanInterfacePara(ofp);

  show_help_HSP2HitList(ofp);

  show_help_HSPset2HitPairPara(ofp);

  show_help_HitListOutputImpl(ofp);

  show_help_DPRunImpl(ofp);

  show_standard_options(ofp);


  fprintf(ofp,"The following options are only applicable to the -seqdb case\n");
  show_help_SeqLookupLoadPara(ofp);
  show_help_ProteinIndexConstructor(ofp);


}


int main(int argc,char ** argv)
{
  DPRunImpl * dpri = NULL;
  ScanWiseHSPImpl * scani = NULL;
  HSP2HitListImpl * hsp2hiti = NULL;
  HitListOutputImpl * hloi = NULL;
  ProteinIndexConstructor * pic = NULL;
 


  HSPScanInterface * hsi;
  HSPScanInterfacePara * para;
  SearchStatInterface * ssi;
  SearchStatInterface * ssl;
  SeqLookupLoadPara * slp;

  HSPset2HitPairPara * hsp2hit;
  CompMat * mat;
  SequenceDB * db;
  Sequence * seq;
  int ret;
  int i;
  int effective_db_size = 300000;
  int kk;
  
  int count = 0;

  LinearHSPmanager * lm;
  HitList * hl;

  boolean use_mott = 1;

  boolean trunc_best_hsp = 0;
  boolean verbose = 0;
  static struct rusage use;

  struct timeval t0, t1;

  gettimeofday(&t0, NULL);


  dpri      = new_DPRunImpl_from_argv(&argc,argv);

  dpri->memory = DPIM_Explicit;

  scani     = new_ScanWiseHSPImpl_from_argv(&argc,argv);
  
  hsp2hiti  = new_HSP2HitListImpl_from_argv(&argc,argv);

  hloi = new_HitListOutputImpl_from_argv(&argc,argv);

  slp = new_SeqLookupLoadPara_from_argv(&argc,argv);

  pic = new_ProteinIndexConstructor_from_argv(&argc,argv);

  hsp2hit = new_HSPset2HitPairPara_from_argv(&argc,argv);

  para = new_HSPScanInterfacePara_from_argv(&argc,argv);

  verbose = strip_out_boolean_argument(&argc,argv,"verbose") ;


  strip_out_boolean_def_argument(&argc,argv,"mott",&use_mott);

  strip_out_boolean_def_argument(&argc,argv,"besthsp",&trunc_best_hsp);

  strip_out_integer_argument(&argc,argv,"dbsize",&effective_db_size);

  

#ifdef SCAN_CORBA
  sorb = get_Wise2Corba_Singleton(&argc,argv,"orbit-local-orb");
#endif

  if( dpri == NULL ) {
    fatal("Unable to build DPRun implementation. Bad arguments");
  }

  strip_out_standard_options(&argc,argv,show_help,show_version);
  if( argc != 2 ) {
    show_help(stdout);
    exit(12);
  }

  /* ugly, but we don't want to bounce matrices around the network... */

  mat = read_Blast_file_CompMat("BLOSUM62.bla");
  
  erroroff(REPORT);

  hsi = new_HSPScanInterface_from_ScanWiseHSPImpl(scani,pic,slp);

  ssi = new_Mott_SearchStatInterface();

  ssl = new_lookup_SearchStatInterface(40,2.3);


  if( verbose ) {
    info("contacted database");
  }

  db = single_fasta_SequenceDB(argv[1]);

  if( db == NULL ) {
    fatal("Could not open sequence db...\n");
  }

  for(seq = init_SequenceDB(db,&ret); seq != NULL;seq = get_next_SequenceDB(db) ) {

	count++;

    for(i=0;i<seq->len;i++) {
      if( !isalpha(seq->seq[i]) ) {
	fatal("Sequence position %d [%c] is not valid",i,seq->seq[i]);
      }
      seq->seq[i] = toupper(seq->seq[i]);
    }

    info("Processing %s",seq->name);

    getrusage(RUSAGE_SELF,&use);
    
    /*    info("Before query %s %.3fu %.3fs\n", seq->name,
	 use.ru_utime.tv_sec + use.ru_utime.tv_usec*MICROSECOND,
	 use.ru_stime.tv_sec + use.ru_stime.tv_usec*MICROSECOND
	);
    */

    lm = (*hsi->scan_query)(hsi->data,seq,para);


    fprintf(stderr,"Got linear manager is %d entries\n",lm->len);

    if( lm->mat == NULL ) {
      lm->mat = hard_link_CompMat(mat);
    }

    getrusage(RUSAGE_SELF,&use);
    /*
    info("After query %s %.3fu %.3fs\n", seq->name,
	 use.ru_utime.tv_sec + use.ru_utime.tv_usec*MICROSECOND,
	 use.ru_stime.tv_sec + use.ru_stime.tv_usec*MICROSECOND
	);
    */
    sort_LinearHSPmanager(lm,compare_HSPset_score);


    if( trunc_best_hsp == 1 ) {
      for(kk=1;kk<lm->len;kk++) {
	free_HSPset(lm->set[kk]);
	lm->set[kk] = NULL;
      }
      lm->len = 1;
    }

    getrusage(RUSAGE_SELF,&use);
    
    /*
    info("After sort %s %.3fu %.3fs\n", seq->name,
	 use.ru_utime.tv_sec + use.ru_utime.tv_usec*MICROSECOND,
	 use.ru_stime.tv_sec + use.ru_stime.tv_usec*MICROSECOND
	);
    */
    hl   = HitList_from_HSP_HSP2HitListImpl(hsp2hiti,lm,dpri,hsp2hit);


    getrusage(RUSAGE_SELF,&use);
    /*
    info("After conversion %s %.3fu %.3fs\n", seq->name,
	 use.ru_utime.tv_sec + use.ru_utime.tv_usec*MICROSECOND,
	 use.ru_stime.tv_sec + use.ru_stime.tv_usec*MICROSECOND
	);
    */
    free_LinearHSPmanager(lm);

    if( use_mott == 1 ) {
      apply_SearchStat_to_HitList(hl,ssi,effective_db_size);
    } else {
      for(kk=0;kk<hl->len;kk++) {
	hl->pair[kk]->bit_score = hl->pair[kk]->raw_score / 2.0; 
      }
    }

    sort_HitList_by_score(hl);

    show_HitList_HitListOutputImpl(hloi,hl,stdout);

    getrusage(RUSAGE_SELF,&use);
    /*
    info("After output %s %.3fu %.3fs\n", seq->name,
	 use.ru_utime.tv_sec + use.ru_utime.tv_usec*MICROSECOND,
	 use.ru_stime.tv_sec + use.ru_stime.tv_usec*MICROSECOND
	);
    */

    free_HitList(hl);
    free_Sequence(seq);
  }
    

  free_DPRunImpl(dpri);
  free_HSPScanInterface(hsi);

  gettimeofday(&t1, NULL);
  fprintf(stderr, "[client stats] queries, time (s): %d %f\n",
                count,
		(t1.tv_sec - t0.tv_sec) +
                (t1.tv_usec - t0.tv_usec) * 1e-6);

  return 0;

}
