
#include "hspscan_server_impl.h"
#include "sequencedb.h"
#include "seqlookup.h"
#include "subseqhash.h"
#include "subseqlookup.h"
#include "hsplookupscan.h"
#include "arrayseqlookup.h"

#include <stdio.h>

#include "corba_singleton.h"

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

  sorb = get_Wise2Corba_Singleton(&argc,argv,"orbit-local-orb");

  db = single_fasta_SequenceDB(argv[1]);

  mat = read_Blast_file_CompMat("blosum62.bla");


  if( strcmp(argv[2],"hash") == 0 ) {
    sli = new_ghash_SeqLookupInterface();
  } else {
    is_array = 1;
    sli = new_ArraySeq_SeqLookupInterface(26*26*26*26*26); 
  }
  
  count = 0;
  for(i=0,seq = init_SequenceDB(db,&ret); seq != NULL;seq = get_next_SequenceDB(db) ) {
    load_aa_flat_Sequence_SeqLookupInterface(sli,hard_link_Sequence(seq));
    if( i % 1000 == 0 ) {
      printf("Loaded up to %d sequences...\n",i);
    }
    i++;
    /*
    if( i > 200000 ) {
      fprintf(stderr,"FOR TESTING - breaking after 200000 seqs\n");
      break;
    }
    */
  }

  if( is_array == 1 && 0) {
    ofp = fopen("array.oc","w");
    fprintf(stderr,"Array occupany...\n");
    print_array_occuypancy_ArraySeq((ArraySeqLookup*)sli->data,ofp);
    fprintf(stderr,"Finished write\n");
    fclose(ofp);
  }


  hsi = new_one_off_HSPScanInterface(sli,mat,40,10);
  

  hspf = new_Wise2HSP_HSPmanagerFactory(sorb->poa,hsi,sorb->ev);

  
  objref  = CORBA_ORB_object_to_string(sorb->orb,hspf, sorb->ev);

  ofp = fopen("hsp.ior","w");
  fprintf(ofp, "%s", objref);
  fclose(ofp);

  pm = PortableServer_POA__get_the_POAManager(sorb->poa, sorb->ev);

  PortableServer_POAManager_activate(pm, sorb->ev);

  fprintf(stderr,"Going to start running...\n");

  CORBA_ORB_run(sorb->orb, sorb->ev);
  return 0;
}
