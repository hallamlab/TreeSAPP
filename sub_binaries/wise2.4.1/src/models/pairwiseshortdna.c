#ifdef _cplusplus
extern "C" {
#endif
#include "pairwiseshortdna.h"

# line 18 "pairwiseshortdna.dy"
boolean process_HSP(HSPset * set,Sequence * query,int query_pos,Sequence * tseq,SeqLookupResultStruct * res_struct,CompMat * mat)
{
  HSP * h;
  int k;

  for(k=0;k<set->len;k++) {
    if( on_HSP(set->hsp[k],query_pos,res_struct->pos) == TRUE ) {
      return FALSE;
    }
  }

  /* new HSP - extend and add */

  h = new_HSP(NULL,query,tseq,query_pos,res_struct->pos,mat,10*Probability2Score(0.8/0.25));
  /*  fprintf(stderr,"Processing HSP with %d %d\n",h->length,h->score); */

  add_HSPset(set,h);


  return TRUE;
}


# line 41 "pairwiseshortdna.dy"
PairwiseShortDna * query_to_reverse_target(Sequence * query,Sequence * target,DnaMatrix * dm,int qstart,int qend,int tstart,int tend)
{
  PairwiseShortDna * out;
  HSPset * forward;
  HSPset * reverse;
  CompMat * cm;
  SeqLookupInterface * sli;
  SeqLookupClientInterface * slci;
  SeqLookupResultInterface * res;
  SeqLookupResultStruct * res_struct = NULL;
  Sequence * revseq;
  int i;
  int no;

  GTree * forward_tree;
  GTree * reverse_tree;

  out = PairwiseShortDna_alloc();

  cm = new_CompMat_from_DnaMatrix_flat(dm);

  /* reverse sequence */

  revseq = reverse_complement_Sequence(target);

  forward_tree = g_tree_new(g_int_equal);
  reverse_tree = g_tree_new(g_int_equal);




  /* build a hash based lookup table */

  sli = new_ghash_SeqLookupInterface();



  assert(cm);
  assert(target);
  assert(revseq);
  assert(sli);


  for(i=tstart;i<tend-7;i++) {
    if( i%1000 == 0 ) {
      info("Loaded %d positions in %s",i,target->name);
    }
    
    no = seq_number_dna_7mer_noN(target->seq+i);
    if( no != -1 ) {
      (*sli->add_direct_number)(sli->data,no,target,i);
    }

  }

  for(i=target->len-tend;i<target->len-tstart-7;i++) {
    if( i%1000 == 0 ) {
      info("Loaded %d positions in %s",i,target->name);
    }
    
    no = seq_number_dna_7mer_noN(revseq->seq+i);
    if( no != -1 ) {
      (*sli->add_direct_number)(sli->data,no,revseq,i);
    }

  }


  forward = HSPset_alloc_std();
  reverse = HSPset_alloc_std();

  /* scan the table with the query, testing whether new hits fit into HSP */

  slci = (*sli->get_client)(sli->data);

  for(i=qstart;i<qend-7;i++) {
    if( i%1000 == 0 ) {
      info("Considered %d positions in %s [forward %d HSPs, backward %d HSPs]\n",i,query->name,forward->len,reverse->len);
	}
    no = seq_number_dna_7mer_noN(query->seq+i);
    if( no == -1 ) {
      continue;
    }

    if( (*slci->is_populated)(sli->data,no) == FALSE) {
      continue;
    }

    res = (*slci->lookup)(sli->data,no);
    while( (*res->is_more)(res->data) ) {
      res_struct = (*res->next)(res->data,res_struct);
      if( res_struct->seq == target ) {
	process_HSP(forward,query,i,target,res_struct,cm);
      } else {
	process_HSP(reverse,query,i,revseq,res_struct,cm);
      }
    }
    free_SeqLookupResultInterface(res);
  } 

  free_SeqLookupClientInterface(slci);

  out->forward = forward;
  out->reverse = reverse;

  return out;
}

# line 138 "pairwiseshortdna.c"
/* Function:  hard_link_PairwiseShortDna(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [PairwiseShortDna *]
 *
 * Return [UNKN ]  Undocumented return value [PairwiseShortDna *]
 *
 */
PairwiseShortDna * hard_link_PairwiseShortDna(PairwiseShortDna * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a PairwiseShortDna object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  PairwiseShortDna_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PairwiseShortDna *]
 *
 */
PairwiseShortDna * PairwiseShortDna_alloc(void) 
{
    PairwiseShortDna * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(PairwiseShortDna *) ckalloc (sizeof(PairwiseShortDna))) == NULL)    {  
      warn("PairwiseShortDna_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->forward = NULL; 
    out->reverse = NULL; 


    return out;  
}    


/* Function:  free_PairwiseShortDna(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [PairwiseShortDna *]
 *
 * Return [UNKN ]  Undocumented return value [PairwiseShortDna *]
 *
 */
PairwiseShortDna * free_PairwiseShortDna(PairwiseShortDna * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a PairwiseShortDna obj. Should be trappable");  
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
    if( obj->forward != NULL)    
      free_HSPset(obj->forward);     
    if( obj->reverse != NULL)    
      free_HSPset(obj->reverse);     


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
