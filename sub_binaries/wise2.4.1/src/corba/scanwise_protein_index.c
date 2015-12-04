
#include "hspscan_server_impl.h"
#include "sequencedb.h"
#include "seqlookup.h"
#include "subseqhash.h"
#include "subseqlookup.h"
#include "hsplookupscan.h"
#include "arrayseqlookup.h"
#include "proteinstreamedindex.h"

#include "../models/version.h"

#include <stdio.h>

#include "corba_singleton.h"

char * program_name = "scanwise_protein_index";

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
  fprintf(ofp,"the scanwise family of methods. It runs a CORBA server which is identified\n");
  fprintf(ofp,"by the IOR file written out. For most databases you need around 10GB of memory\n\n");
  fprintf(ofp," OPTIONS\n");
  fprintf(ofp,"   -iorfile  file to write IOR out to [hsp.ior]\n");
  fprintf(ofp,"   -hash     use glib hash not array (far slower, but easier on the memory for small DBs)\n");
  fprintf(ofp,"   -streamed use streamed index (more compact, higher run-time)\n");
  fprintf(ofp,"   -ocfile   occupancy file for array/streamed indexes (large), for debugging\n");

  show_help_SeqLookupLoadPara(ofp);
  
  show_standard_options(ofp);
}



int
main(int argc, char* argv[])
{
  int i;
  Wise2Corba_Singleton      * sorb;
  PortableServer_ObjectId*  oid;
  Wise2HSP_HSPmanagerFactory hspf;
  Wise2_HSPScanInterface * hsi;
  SequenceDB * db;
  Sequence * seq;
  int ret;
  SeqLookupInterface * sli;
  CompMat * mat;
  FILE * ofp;
  PortableServer_POAManager pm;
  CORBA_char*               objref;
  int count = 0;
  int is_array = 0;
  int waystation = 4;

  char * ior_file;
  char * oc_file = NULL;
  char * temp;
  
  SeqLookupLoadPara * slp;

  boolean use_hash = 0;
  boolean use_stream = 0;


  slp = new_SeqLookupLoadPara_from_argv(&argc,argv);

  use_hash   = strip_out_boolean_argument(&argc,argv,"hash");
  use_stream = strip_out_boolean_argument(&argc,argv,"streamed");

  strip_out_integer_argument(&argc,argv,"waystation",&waystation);

  if( (temp = strip_out_assigned_argument(&argc,argv,"iorfile")) != NULL ) {
    ior_file = temp;
  } else {
    ior_file = "hsp.ior";
  }

  oc_file = strip_out_assigned_argument(&argc,argv,"ocfile");

  strip_out_standard_options(&argc,argv,show_help,show_version);

  if( argc != 2 ) {
    show_help(stdout);
    exit(12);
  }


  ofp = fopen(ior_file,"w");
  if( ofp == NULL ) {
    fatal("unable to open %s\n",ior_file);
  }
  fclose(ofp);

  sorb = get_Wise2Corba_Singleton(&argc,argv,"orbit-local-orb");

  db = single_fasta_SequenceDB(argv[1]);

  mat = read_Blast_file_CompMat("blosum62.bla");


  if( use_stream == 1 ) {
    sli = new_ProteinStreamedIndex_SeqLookupInterface(waystation);
  } else if( use_hash == 1 ) {
    sli = new_ghash_SeqLookupInterface();
  } else {
    is_array = 1;
    sli = new_ArraySeq_SeqLookupInterface(26*26*26*26*26); 
  }


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


  hsi = Wise2_new_one_off_HSPScanInterface(sli,mat,40,10);
  
  hspf = new_Wise2HSP_HSPmanagerFactory(sorb->poa,hsi,sorb->ev);

  
  objref  = CORBA_ORB_object_to_string(sorb->orb,hspf, sorb->ev);

  ofp = fopen(ior_file,"w");
  if( ofp == NULL ) {
    fatal("unable to open %s\n",ior_file);
  }

  fprintf(ofp, "%s", objref);
  fclose(ofp);

  pm = PortableServer_POA__get_the_POAManager(sorb->poa, sorb->ev);

  PortableServer_POAManager_activate(pm, sorb->ev);

  fprintf(stderr,"IOR written out to %s, serving...\n",ior_file);

  CORBA_ORB_run(sorb->orb, sorb->ev);
  return 0;
}
