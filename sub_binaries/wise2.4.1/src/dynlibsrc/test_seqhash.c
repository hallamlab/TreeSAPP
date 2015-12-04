
#include "seqlookup.h"
#include "subseqhash.h"
#include "sequencedb.h"
#include "hsp.h"
#include "hitlist.h"
#include "hsplookupscan.h"


void print_hsp(gpointer key,gpointer value,gpointer user_data)
{
  Sequence * query;
  HSPset * set;
  btCanvas * btc;
  btPasteArea * btp;
  int i,j;
  char q,t;

  query = (Sequence *)user_data;
  set = (HSPset *) value;

  btc = new_Ascii_btCanvas(stdout,12,50,5,3);  

  for(i=0;i<set->len;i++) {

      
      for(j=0;j<set->hsp[i]->length;) {
	
	btp = get_reserved_left_btCanvas(btc);
	
	paste_string_btPasteArea(btp,0,0,query->name,BC_RIGHT,0);
	paste_string_btPasteArea(btp,0,2,set->hsp[i]->target->name,BC_RIGHT,0);

	free_btPasteArea(btp);

	for(;j<set->hsp[i]->length && can_get_paste_area_btCanvas(btc,1) == TRUE;j++) {
	  btp = get_paste_area_btCanvas(btc,1);
	  q = query->seq[set->hsp[i]->query_start+j];
	  t = set->hsp[i]->target->seq[set->hsp[i]->target_start+j];
	  paste_char_btPasteArea(btp,0,0,q,0);
	  paste_char_btPasteArea(btp,0,2,t,0);
	  if( q == t ) {
	    paste_char_btPasteArea(btp,0,1,t,0);
	  }

	  free_btPasteArea(btp);
	}
	advance_line_btCanvas(btc);
      }

      advance_line_btCanvas(btc);
  }

  free_btCanvas(btc);
  
}

  

int main(int argc,char ** argv)
{
  SequenceDB * db;
  Sequence * seq;
  SeqLookupInterface * sli;
  SeqLookupPos * slp;
  HSPScanInterface * hsi;
  LinearHSPmanager * lm;
  HitList * hl;
  CompMat * mat;
  int ret;
  HSPScanInterfacePara p;


  p.min_score= 30;
  p.max_results = 200;

  db = single_fasta_SequenceDB(argv[1]);

  mat = read_Blast_file_CompMat("blosum62.bla");

  sli = new_ghash_SeqLookupInterface();

  for(seq = init_SequenceDB(db,&ret); seq != NULL;seq = get_next_SequenceDB(db) ) {
    load_aa_flat_Sequence_SeqLookupInterface(sli,hard_link_Sequence(seq));
  }


  seq = read_fasta_file_Sequence(argv[2]);

  assert(seq);

  hsi = Wise2_new_one_off_HSPScanInterface(sli,mat,20,10);

/*  hspm = simple_HSPScan_scan_query((void*)hsi->data,seq); */

  lm = (*hsi->scan_query)(hsi->data,seq,&p); 

  hl = Wise2_HitList_from_LinearHSPmanager(lm);

  Wise2_write_pseudoblast_HitList(hl,stdout);

}


