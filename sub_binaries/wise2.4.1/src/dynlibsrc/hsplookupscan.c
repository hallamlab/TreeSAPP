#ifdef _cplusplus
extern "C" {
#endif
#include "hsplookupscan.h"
#include <sys/time.h>
#include <sys/resource.h>


/* Function:  seq_number_aa_5mer_client(seq)
 *
 * Descrip:    Function for the amino acid to number on 5mers
 *
 *
 * Arg:        seq [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 36 "hsplookupscan.dy"
int seq_number_aa_5mer_client(char * seq)
{
  int i;
  int ret = 0;
  int base = 1;
  int no = 0;

  for(i=0;i<5;i++) {
    no = toupper(seq[i])-'A';
    if( no > 26 || no < 0 ) {
      no = 'X'-'A';
    }
    ret += base * no;
    base = base * 26;
  }

  return ret;
}



/* Function:  new_simple_HSPScanInterface(sli,mat,drop_off)
 *
 * Descrip:    Builds a new simple scan interface. This
 *             does not expand the query using a matrix but
 *             rather simply scans down the query sequence
 *
 *
 * Arg:             sli [UNKN ] Undocumented argument [SeqLookupInterface *]
 * Arg:             mat [UNKN ] Undocumented argument [CompMat *]
 * Arg:        drop_off [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [HSPScanInterface *]
 *
 */
# line 62 "hsplookupscan.dy"
HSPScanInterface * new_simple_HSPScanInterface(SeqLookupInterface * sli,CompMat * mat,int drop_off)
{
  HSPScanInterface * out;
  HSPScanPara * p;

  assert(sli);
  assert(mat);
  
  out = HSPScanInterface_alloc();

  p = HSPScanPara_alloc();
  p->sli = hard_link_SeqLookupInterface(sli);
  p->mat = hard_link_CompMat(mat);
  p->drop_off = drop_off;

  out->data = (void*)p;
  out->free_data = simple_HSPScan_free;
  out->scan_query = simple_HSPScan_scan_query;

  return out;
}


/* Function:  new_one_off_HSPScanInterface(sli,mat,drop_off,score_cutoff)
 *
 * Descrip:    Builds a new simple scan interface. This
 *             does expands the query using a matrix but
 *             just be considering off by one cases
 *
 *
 * Arg:                 sli [UNKN ] Undocumented argument [SeqLookupInterface *]
 * Arg:                 mat [UNKN ] Undocumented argument [CompMat *]
 * Arg:            drop_off [UNKN ] Undocumented argument [int]
 * Arg:        score_cutoff [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [HSPScanInterface *]
 *
 */
# line 90 "hsplookupscan.dy"
HSPScanInterface * new_one_off_HSPScanInterface(SeqLookupInterface * sli,CompMat * mat,int drop_off,int score_cutoff)
{

  HSPScanInterface * out;
  HSPScanPara * p;

  assert(sli);
  assert(mat);
  
  out = HSPScanInterface_alloc();

  p = HSPScanPara_alloc();
  p->sli = hard_link_SeqLookupInterface(sli);
  p->mat = hard_link_CompMat(mat);
  p->drop_off = drop_off;
  p->score_cutoff = score_cutoff;


  out->data = (void*)p;
  out->free_data = simple_HSPScan_free;
  if( sli->lookup_array_head == NULL ) {
    out->scan_query = one_off_HSPscan_scan_query;
  } else {
    out->scan_query = one_off_HSPscan_scan_query_direct;
  }

  return out;
}

/* Function:  no_op_func(data,user_data,data2)
 *
 * Descrip:    provide a no op func
 *
 *
 * Arg:             data [UNKN ] Undocumented argument [gpointer]
 * Arg:        user_data [UNKN ] Undocumented argument [gpointer]
 * Arg:            data2 [UNKN ] Undocumented argument [gpointer]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 122 "hsplookupscan.dy"
int no_op_func(gpointer  data,gpointer  user_data,gpointer  data2)
{
  return 1;
}



/* Function:  one_off_HSPscan_scan_query_direct(data,seq,para)
 *
 * Descrip:    Simple word expansion for direct access
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 * Arg:         seq [UNKN ] Undocumented argument [Sequence *]
 * Arg:        para [UNKN ] Undocumented argument [HSPScanInterfacePara *]
 *
 * Return [UNKN ]  Undocumented return value [LinearHSPmanager *]
 *
 */
# line 132 "hsplookupscan.dy"
LinearHSPmanager * one_off_HSPscan_scan_query_direct(void * data,Sequence * seq,HSPScanInterfacePara * para)
{
  char * std_aa = "ACDEFGHIKLMNPQRSTVWY";
  HSPmanager * hspm;
  LinearHSPmanager * out;
  HSPScanPara * p = (HSPScanPara *)data;
  LineariseHSPPara * hsp_para;
  int i;
  int j;
  int k;
  int aa;
  int score;
  char newseq[5];

  int jj;
  int aa2;

  int has_seen;

  int seqno[5];
  int base[5];
  int start_base;

  int self_score[5];
  int total_score;
  int current_score;

  ArraySeqHead * head;
  int no;

#ifdef LINUX_TIMER
  static struct rusage use;
  struct timeval t1;
  struct timeval t2;
  struct timeval t3;
#endif

  assert(seq);
  assert(p);
  assert(para);
  assert(para->max_results > 0);


  hsp_para = LineariseHSPPara_alloc();
  hsp_para->verbosity   = para->verbosity;
  hsp_para->max_size    = para->max_results;
  hsp_para->min_score   = para->min_hsp_score;
  hsp_para->width       = para->hsp_link_width;
  hsp_para->tail        = para->hsp_link_length;
  

  for(i=0,start_base=1;i<5;i++) {
    base[i] = start_base;
    start_base = start_base * 26;
  }

  if( VERBOSITY_CHECK(1,para->verbosity) ) {
    info("processing sequence %s with vanilla implementation",seq->name);
  }

  if( COMPILE_VERBOSITY_FLAG == 0 && para->verbosity > 0 ) {
    info("Verbosity information has not be enable in this compile. Ignoring verbosity request");
  }



  hspm = new_HSPmanager(seq,p->mat,p->drop_off);


#ifdef LINUX_TIMER 
  if( VERBOSITY_CHECK(1,para->verbosity) ) {
    gettimeofday(&t1,NULL);
  }
#endif


  for(i=0;i<seq->len-5;i++) {

    if( VERBOSITY_CHECK(2,para->verbosity) && i%50 == 0 ) {
      info("Scanning %s at position %d with %d hits",seq->name,i,hspm->hsp_count);
    }

    head = (*p->sli->lookup_array_head)(p->sli->data,seq_number_aa_5mer_client(seq->seq+i));
    if( head != NULL ) {

      if( para->numb_level < head->current_pos  || (ARRAYHEAD_IS_LOWCOMPLEXITY(head) && para->low_numb > 0 && para->low_numb <= head->current_pos) ) {
	if( VERBOSITY_CHECK(7,para->verbosity) ) {
	  info("position %d hit direct over filled position, %d vs hard %d, low %d",i,head->current_pos,para->numb_level,para->low_numb);
	}

	continue;
      }


      for(k=0;k<head->current_pos;k++) {

	if( add_pair_HSPmanager_score(hspm,head->units[k].seq,i,head->units[k].pos,para->min_hsp_score) == TRUE ) {
	  ;
	}
      }

      /*
      if( p->use_msp_crunch == 1 && head != NULL && head->current_pos > p->msp_crunch_no ) {
	continue;
      }
      */
    }

    

    total_score = 0;
    for(score=0,j=0;j<5;j++) {

      seqno[j] = base[j]*(toupper(seq->seq[i+j])-'A');
      self_score[j] = p->mat->comp[toupper(seq->seq[i+j])-'A'][toupper(seq->seq[i+j])-'A'];
      total_score += self_score[j];
    }

    has_seen = 0;

    for(j=0;j<5;j++) {
      for(aa=0;aa<20;aa++) {
	if( toupper(seq->seq[i+j]) == std_aa[aa] ) {
	  continue;
	}

	seqno[j] = base[j]*(std_aa[aa]-'A');

	no= seqno[0]+seqno[1]+seqno[2]+seqno[3]+seqno[4];

	head = (*p->sli->lookup_array_head)(p->sli->data,no);
	if( head != NULL ) {

	  if( para->numb_level < head->current_pos  || (ARRAYHEAD_IS_LOWCOMPLEXITY(head) && para->low_numb > 0 && para->low_numb <= head->current_pos) ) {
	    if( VERBOSITY_CHECK(7,para->verbosity) ) {
	      info("position %d hit one-off over filled position, %d vs hard %d, low %d",i,head->current_pos,para->numb_level,para->low_numb);
	    }

	    continue;
	  }

	  has_seen = 1;
	  for(k=0;k<head->current_pos;k++) {
	    add_pair_HSPmanager(hspm,head->units[k].seq,i,head->units[k].pos) ;
	  }
	}

	seqno[j]  = base[j]*(toupper(seq->seq[i+j])-'A');
	newseq[j] = seq->seq[i+j];
      }


      /* really at this point we want to see if we have seen any really, really
	 long HSPs, or rather enough long HSPs to remove the need for short
	 HSPs. This is hard because good with two gaps in them will not
	 have potentially three pretty short HSPs. Hmmmm.
      */


      if( 1 || has_seen == 0 ) {

	  /* latin square optimisation. The pattern of 1,1,0,1,0 of active
	     swaps will cover all 2-off positions as the tile is progressed */

	for(j=0;j<5;j++) {
	  for(jj=1;jj<5;jj+=2) {
	    for(aa=0;aa<20;aa++) {
	      for(aa2=0;aa2<20;aa2++) {


		if( seq->seq[i+j] == std_aa[aa] || seq->seq[i+jj] == std_aa[aa2] ) {
		  continue;
		}

		current_score = total_score;
		current_score -= (self_score[j] - p->mat->comp[toupper(seq->seq[i+j])-'A'][std_aa[aa]-'A']);
		current_score -= (self_score[jj] - p->mat->comp[toupper(seq->seq[i+jj])-'A'][std_aa[aa2]-'A']);

		if( current_score < para->min_word_score ) {
		  continue;
		} else {
		  /* fprintf(stderr,"Handling %d,%d amino acid %d,%d score %d\n",j,jj,aa,aa2,current_score); */
		}

		seqno[j]  = base[j]*(std_aa[aa]-'A');
		seqno[jj] = base[jj]*(std_aa[aa2]-'A');

		no= seqno[0]+seqno[1]+seqno[2]+seqno[3]+seqno[4];
		
		head = (*p->sli->lookup_array_head)(p->sli->data,no);
		if( head != NULL ) {

		  if( para->numb_level < head->current_pos  || (ARRAYHEAD_IS_LOWCOMPLEXITY(head) && para->low_numb > 0 && para->low_numb <= head->current_pos) ) {
		    if( VERBOSITY_CHECK(7,para->verbosity) ) {
		      info("position %d hit direct over filled position, %d vs hard %d, low %d",i,head->current_pos,para->numb_level,para->low_numb);
		    }
		    
		    continue;
		  }

		  has_seen = 1;
		  for(k=0;k<head->current_pos;k++) {
		    add_pair_HSPmanager_score(hspm,head->units[k].seq,i,head->units[k].pos,para->min_hsp_score) ;
		  }
		}
		
		seqno[j]   = base[j]*(toupper(seq->seq[i+j]-'A'));
		seqno[jj]  = base[jj]*(toupper(seq->seq[i+jj]-'A'));

		newseq[j] = seq->seq[i+j];
		newseq[jj] = seq->seq[i+jj];
	      }
	    }
	  }
	}
      }
      
    }  
  }


#ifdef LINUX_TIMER
  if( VERBOSITY_CHECK(1,para->verbosity) ) {
    getrusage(RUSAGE_SELF,&use);
    
    gettimeofday(&t2,NULL);
  }

#endif

  if( para->use_protein_heuristic == TRUE ) {

    /* out = new_LinearHSPmanager_flat(hspm);*/

    out = new_LinearHSPmanager_simple_heuristic(hspm,hsp_para); 
  } else {
    out = new_LinearHSPmanager_flat(hspm);
  }

  fprintf(stderr,"Made linear manager with %d elements\n",out->len);

#ifdef LINUX_TIMER
  
  if( VERBOSITY_CHECK(1,para->verbosity) ) {
    gettimeofday(&t3,NULL);	

    getrusage(RUSAGE_SELF,&use);
    
    info("Sort clock point: Scan %f : Sort %f",
	 t2.tv_sec - t1.tv_sec + ((t2.tv_usec - t1.tv_usec) * 1e-6),
	 t3.tv_sec - t2.tv_sec + ((t2.tv_usec - t2.tv_usec) * 1e-6)
	 );
  }
#endif


  /* free_HSPmanager(hspm);*/

  return out;
}



/* Function:  one_off_HSPscan_scan_query(data,seq,para)
 *
 * Descrip:    Simple word expansion - one off score drop considered
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 * Arg:         seq [UNKN ] Undocumented argument [Sequence *]
 * Arg:        para [UNKN ] Undocumented argument [HSPScanInterfacePara *]
 *
 * Return [UNKN ]  Undocumented return value [LinearHSPmanager *]
 *
 */
# line 398 "hsplookupscan.dy"
LinearHSPmanager * one_off_HSPscan_scan_query(void * data,Sequence * seq,HSPScanInterfacePara * para)
{
  HSPmanager * hspm;
  LinearHSPmanager * out;
  HSPScanPara * p = (HSPScanPara *)data;
  int i;
  int j;
  int aa;
  int score;
  int scorepos[5];
  char newseq[5];
  SeqLookupResultInterface * slri;
  SeqLookupClientInterface * slci;

  SeqLookupResultStruct * res = NULL;

  int trace = 0;

  assert(seq);
  assert(p);
  assert(para);
  assert(para->max_results > 0);

  slci = (*p->sli->get_client)(p->sli->data);

  hspm = new_HSPmanager(seq,p->mat,p->drop_off);

  for(i=0;i<seq->len-5;i++) {
    if( trace ) 
      fprintf(stderr,"Looking at %d %.5s (straight)\n",i,seq->seq+i); 

    if( (*slci->is_populated)(slci->data,seq_number_aa_5mer(seq->seq+i)) ) {

      if( trace ) 
	fprintf(stderr,".... is populated (straight) %d,%d\n",i,seq->seq[i]); 


      slri = (*slci->lookup)(slci->data,seq_number_aa_5mer(seq->seq+i));

      if( trace ) 
	fprintf(stderr,".... have results (straight) %d,%d\n",i,seq->seq[i]); 

      res = NULL;
      for(;(*slri->is_more)(slri->data);) {    
	res = (*slri->next)(slri->data,res);
	add_pair_HSPmanager(hspm,res->seq,i,res->pos);

	if( trace ) 
	  fprintf(stderr,"...adding direct %.5s\n",res->seq->seq+res->pos); 

      }
      free_SeqLookupResultInterface(slri);
    }


    if( trace ) {
      fprintf(stderr,"Expanding the sequence\n");
    }

    for(score=0,j=0;j<5;j++) {
      newseq[j] = seq->seq[i+j];
      scorepos[j] = p->mat->comp[seq->seq[i+j]-'A'][seq->seq[i+j]-'A'];
      score += p->mat->comp[seq->seq[i+j]-'A'][seq->seq[i+j]-'A'];
    }

    for(j=0;j<5;j++) {
      for(aa=0;aa<26;aa++) {
	if( score - scorepos[j] + p->mat->comp[seq->seq[i+j]-'A'][aa] > p->score_cutoff ) {
	  newseq[j] = aa+'A';

	  if( trace ) {
	    fprintf(stderr,"Seeing if is populated on expansion %d,%d\n",j,aa);
	  }
	  

	  if( (*slci->is_populated)(slci->data,seq_number_aa_5mer(newseq)) ) {

	    if( trace )
	      fprintf(stderr,"...is populated %.5s\n",res->seq->seq+res->pos); 

	    slri = (*slci->lookup)(slci->data,seq_number_aa_5mer(newseq));
	    res = NULL;
	    for(;(*slri->is_more)(slri->data);) {
	      res = (*slri->next)(slri->data,res);
	      add_pair_HSPmanager(hspm,res->seq,i,res->pos);


	      if( trace )
		fprintf(stderr,"...adding one off %.5s\n",res->seq->seq+res->pos); 

	    }
	    free_SeqLookupResultInterface(slri);
	  }

	  newseq[j] = seq->seq[i+j];
	}
      }
    }	  
  }

  free_SeqLookupClientInterface(slci);

  if( para->use_protein_heuristic == TRUE ) {
    out = new_LinearHSPmanager_heuristic_max(hspm,para->max_results);
  } else {
    out = new_LinearHSPmanager_flat(hspm);
  }


  free_HSPmanager(hspm);

  return out;
}

/* Function:  simple_HSPScan_scan_query(data,seq,para)
 *
 * Descrip:    simple Scan function, no word expansions
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 * Arg:         seq [UNKN ] Undocumented argument [Sequence *]
 * Arg:        para [UNKN ] Undocumented argument [HSPScanInterfacePara *]
 *
 * Return [UNKN ]  Undocumented return value [LinearHSPmanager *]
 *
 */
# line 515 "hsplookupscan.dy"
LinearHSPmanager * simple_HSPScan_scan_query(void * data,Sequence * seq,HSPScanInterfacePara * para)
{
  HSPmanager * hspm;
  LinearHSPmanager * out;
  HSPScanPara * p = (HSPScanPara *) data;
  int i;
  SeqLookupResultInterface * slri;
  SeqLookupClientInterface * slci;
  SeqLookupResultStruct * res = NULL;



  hspm = new_HSPmanager(seq,p->mat,p->drop_off);

  slci = (*p->sli->get_client)(p->sli->data);
  assert(slci);

  for(i=0;i<seq->len-5;i++) {

    if( (*slci->is_populated)(slci->data,seq_number_aa_5mer(seq->seq+i)) ) {
      slri = (*slci->lookup)(slci->data,seq_number_aa_5mer(seq->seq+i));
      for(;(*slri->is_more)(slri->data);) {
	res = (*slri->next)(slri->data,res);
	add_pair_HSPmanager(hspm,res->seq,i,res->pos);
      }
      free_SeqLookupResultInterface(slri);
    }
  }

  free_SeqLookupClientInterface(slci);

  out = new_LinearHSPmanager_heuristic_max(hspm,para->max_results);

  free_HSPmanager(hspm);

  return out;
}


/* Function:  simple_HSPScan_free(data)
 *
 * Descrip:    Free function for simple scans
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 *
 */
# line 557 "hsplookupscan.dy"
void simple_HSPScan_free(void * data)
{
  HSPScanPara * p = (HSPScanPara *) data;

  free_HSPScanPara(p);
}
# line 597 "hsplookupscan.c"
/* Function:  hard_link_HSPScanPara(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [HSPScanPara *]
 *
 * Return [UNKN ]  Undocumented return value [HSPScanPara *]
 *
 */
HSPScanPara * hard_link_HSPScanPara(HSPScanPara * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a HSPScanPara object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  HSPScanPara_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [HSPScanPara *]
 *
 */
HSPScanPara * HSPScanPara_alloc(void) 
{
    HSPScanPara * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(HSPScanPara *) ckalloc (sizeof(HSPScanPara))) == NULL)  {  
      warn("HSPScanPara_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->sli = NULL; 
    out->mat = NULL; 
    out->drop_off = 0;   
    out->score_cutoff = 0;   
    out->use_msp_crunch = 1; 
    out->msp_crunch_no = 10; 
    out->seed_factor = 5;    
    out->twohit_wobble = 5;  
    out->threadno = 1;   


    return out;  
}    


/* Function:  free_HSPScanPara(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [HSPScanPara *]
 *
 * Return [UNKN ]  Undocumented return value [HSPScanPara *]
 *
 */
HSPScanPara * free_HSPScanPara(HSPScanPara * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a HSPScanPara obj. Should be trappable");   
      return NULL;   
      }  


#ifdef PTHREAD   
    assert(pthread_mutex_lock(&(obj->dynamite_mutex)) == 0); 
#endif   
    if( obj->dynamite_hard_link > 1)     {  
      return_early = 1;  
      obj->dynamite_hard_link--; 
      }  
#ifdef PTHREAD   
    assert(pthread_mutex_unlock(&(obj->dynamite_mutex)) == 0);   
#endif   
    if( return_early == 1)   
      return NULL;   
    if( obj->sli != NULL)    
      free_SeqLookupInterface(obj->sli);     
    if( obj->mat != NULL)    
      free_CompMat(obj->mat);    


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
