#include "seqlookup.h"
#include "sequencedb.h"
#include "hspscan_corba_wrapper.h"
#include "corba_singleton.h"
#include "hitlist.h"
  

int main(int argc,char ** argv)
{
  Sequence * seq;
  HSPScanInterface * hsi;
  HSPScanInterfacePara p;
  LinearHSPmanager * lm;

  HitList * hl;
  CompMat * mat;
  int ret;
  Wise2Corba_Singleton * sorb;

  

  mat = read_Blast_file_CompMat("blosum62.bla");

  sorb = get_Wise2Corba_Singleton(&argc,argv,"orbit-local-orb");

  hsi = new_corba_HSPScan(sorb,"hsp.ior",mat);

  seq = read_fasta_file_Sequence(argv[1]);

  assert(seq);

  p.max_results = 30;
  lm = (*hsi->scan_query)(hsi->data,seq,&p); 


  hl = HitList_from_LinearHSPmanager(lm);

  write_pseudoblast_HitList(hl,stdout);

}
