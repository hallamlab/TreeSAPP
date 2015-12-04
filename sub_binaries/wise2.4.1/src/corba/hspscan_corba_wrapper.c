#include "hspscan_corba_wrapper.h"
#include "hspscan_corba.h"



typedef struct hspscan_corba_data {
  Wise2Corba_Singleton * sorb;
  Wise2HSP_HSPmanagerFactory hspf;
  Wise2_CompMat * mat;
  int verbose;
} hspscan_corba_data;



void corba_HSPScan_free_data(void * data)
{
  hspscan_corba_data * corba = (hspscan_corba_data * ) data;

  CORBA_Object_release(corba->hspf,corba->sorb->ev);

  free(corba);
  return;
}



LinearHSPmanager * corba_HSPScan_scan_query_iterator(void * data,Sequence * query,HSPScanInterfacePara * para)
{
  LinearHSPmanager * lm;
  HSPset * set;
  HSP * hsp;

  Wise2HSP_HSP_ts_iterator iter;
  Wise2HSP_HSPtargetset * ts;
  Wise2_Sequence * temp_target; 
  Wise2HSP_Sequence corba_query;

  hspscan_corba_data * hsp_data = (hspscan_corba_data*) data;
  int i,j;

  if( hsp_data->verbose == 1 ) {
    fprintf(stderr,"Loading up query... %s\n",query->name);
  }

  corba_query.name = CORBA_string_dup(query->name);
  corba_query.seq  = CORBA_string_dup(query->seq);

  iter = Wise2HSP_HSPmanagerFactory_scan_query_iterator(hsp_data->hspf,&corba_query,250,hsp_data->sorb->ev);

  if( hsp_data->verbose == 1 ) {
    fprintf(stderr,"Retrieved hsp iterator...\n");
  }


  lm = LinearHSPmanager_alloc_std();
  lm->query = hard_link_Sequence(query);
  lm->mat = hard_link_CompMat(hsp_data->mat);

  while( Wise2HSP_HSP_ts_iterator_is_more(iter,hsp_data->sorb->ev) ) {
    
    ts = Wise2HSP_HSP_ts_iterator_next(iter,hsp_data->sorb->ev);

    set = HSPset_alloc_std();
    add_LinearHSPmanager(lm,set);
    set->score = 0;
    set->best_score= 0;

    temp_target = new_Sequence_from_strings(ts->target.name,ts->target.seq);
    
    for(j=0;j<ts->hsp._length;j++) {
      hsp = HSP_alloc();
      hsp->query = hard_link_Sequence(query);
      hsp->target = hard_link_Sequence(temp_target);
      hsp->score = ts->hsp._buffer[j].score;
      hsp->query_start = ts->hsp._buffer[j].query_start;
      hsp->target_start = ts->hsp._buffer[j].target_start;
      hsp->length = ts->hsp._buffer[j].length;

      set->score += hsp->score;
      if( set->best_score < hsp->score ) {
	set->best_score = hsp->score;
      }

      add_HSPset(set,hsp);

    }
    free_Sequence(temp_target); 
  }

  if( hsp_data->verbose == 1 ) {
    fprintf(stderr,"Built Wise2 hspmanager, going to return\n");
  }
  
  return lm;

}


Wise2_HSPScanInterface * new_corba_HSPScan(Wise2Corba_Singleton * sorb,char * iorfile,CompMat * client_side_mat)
{
  hspscan_corba_data * data;
  Wise2_HSPScanInterface * out;
  FILE * ifp;
  char buffer[1024];

  assert(sorb);
  assert(iorfile);


  ifp = fopen(iorfile,"r");
  if( ifp == NULL ) {
    fprintf(stderr,"Could not open file %s!",iorfile);
    return NULL;
  }

  fgets(buffer,1023,ifp);
  fclose(ifp);

  out = Wise2_HSPScanInterface_alloc();
  data = g_new0(hspscan_corba_data,1);
  out->data = (void *) data;
  data->sorb = sorb;
  data->mat = hard_link_CompMat(client_side_mat);
  
  data->hspf = CORBA_ORB_string_to_object(sorb->orb, buffer, sorb->ev);
  data->verbose = 1;

  out->scan_query = corba_HSPScan_scan_query_iterator;
  out->free_data = corba_HSPScan_free_data;
  return out;
}

