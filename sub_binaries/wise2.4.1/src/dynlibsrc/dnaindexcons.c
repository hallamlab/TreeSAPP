#ifdef _cplusplus
extern "C" {
#endif
#include "dnaindexcons.h"

/* Function:  show_help_DnaIndexConstructor(ofp)
 *
 * Descrip:    Shows help for a DnaIndexConstructor
 *
 *
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 33 "dnaindexcons.dy"
void show_help_DnaIndexConstructor(FILE * ofp)
{
  fprintf(ofp,"DnaIndex Constructor options\n");
  fprintf(ofp,"  -dic_type [hash] type of DnaIndex\n");
  fprintf(ofp,"  -dic_word [7] index word length\n");
  fprintf(ofp,"  -dic_prob [0.85] match probabilty in extension\n");
  fprintf(ofp,"  -dic_drop_off [30] drop off in HSP extension\n");

}

/* Function:  new_DnaIndexConstructor(argc,argv)
 *
 * Descrip:    Builds new DnaIndexConstructor off Command line
 *
 *
 * Arg:        argc [UNKN ] Undocumented argument [int *]
 * Arg:        argv [UNKN ] Undocumented argument [char **]
 *
 * Return [UNKN ]  Undocumented return value [DnaIndexConstructor *]
 *
 */
# line 46 "dnaindexcons.dy"
DnaIndexConstructor * new_DnaIndexConstructor(int * argc,char ** argv)
{
  DnaIndexConstructor * out;
  char * name;
  
  out = DnaIndexConstructor_alloc();

  if( (name = strip_out_assigned_argument(argc,argv,"dic_type")) != NULL ) {
    if( strcmp(name,"hash") == 0 ) {
      out->type = DnaIndexConstructor_subseq;
    } else {
      fatal("Cannot recognise %s as a potential dna index type",name);
    }
  }

  strip_out_integer_argument(argc,argv,"dic_word",&out->index_word_length);
  strip_out_float_argument(argc,argv,"dic_prob",&out->match_prob);
  strip_out_integer_argument(argc,argv,"dic_drop_off",&out->drop_off);

  
  return out;

}

/* Function:  DnaIndex_from_DnaIndexConstructor(dic)
 *
 * Descrip:    Makes a DnaIndex from a DnaIndexConstructore
 *
 *
 * Arg:        dic [UNKN ] Undocumented argument [DnaIndexConstructor *]
 *
 * Return [UNKN ]  Undocumented return value [DnaIndex *]
 *
 */
# line 73 "dnaindexcons.dy"
DnaIndex * DnaIndex_from_DnaIndexConstructor(DnaIndexConstructor * dic)
{
  DnaIndex * out;
  DnaProbMatrix * dmp;
  DnaMatrix * dm;

  assert(dic);

  out = DnaIndex_alloc();

  switch(dic->type) {
  case DnaIndexConstructor_subseq :
    out->sli = new_ghash_SeqLookupInterface();
    break;
  default :
    warn("Unable to make DnaIndex from type %d",dic->type);
    return NULL;
  }

  out->index_word_length = dic->index_word_length;


  dmp = DnaProbMatrix_from_match(dic->match_prob,NMaskType_BANNED);  
  assert(dmp);
  flat_null_DnaProbMatrix(dmp);  

  dm = DnaMatrix_from_DnaProbMatrix(dmp);

  out->cm = new_CompMat_from_DnaMatrix_flat(dm);
  out->drop_off = dic->drop_off;

  free_DnaMatrix(dm);
  free_DnaProbMatrix(dmp);
  
  return out;
}


/* Function:  LinearHSPManager_scan_DnaIndex(di,query)
 *
 * Descrip:    provides a LinearManager for a DNA sequence 
 *
 *
 * Arg:           di [UNKN ] Undocumented argument [DnaIndex *]
 * Arg:        query [UNKN ] Undocumented argument [Sequence *]
 *
 * Return [UNKN ]  Undocumented return value [LinearHSPmanager *]
 *
 */
# line 114 "dnaindexcons.dy"
LinearHSPmanager * LinearHSPManager_scan_DnaIndex(DnaIndex * di,Sequence * query)
{
  LinearHSPmanager * out;
  HSPmanager * hspm;

  hspm = HSPmanager_scan_DnaIndex(di,query);

  out = new_LinearHSPmanager_flat(hspm);

  free_HSPmanager(hspm);

  return out;
}

/* Function:  HSPmanager_scan_DnaIndex(di,seq)
 *
 * Descrip:    Provides a HSPmanager from a scan
 *
 *
 * Arg:         di [UNKN ] Undocumented argument [DnaIndex *]
 * Arg:        seq [UNKN ] Undocumented argument [Sequence *]
 *
 * Return [UNKN ]  Undocumented return value [HSPmanager *]
 *
 */
# line 131 "dnaindexcons.dy"
HSPmanager * HSPmanager_scan_DnaIndex(DnaIndex * di,Sequence * seq)
{
  HSPmanager * out;
  SeqLookupClientInterface * slci;
  SeqLookupResultInterface * slri;
  SeqLookupResultStruct * res = NULL;

  int i;
  int no;

  assert(seq);
  assert(di);

  out = new_HSPmanager(seq,di->cm,di->drop_off);
  slci = (*di->sli->get_client)(di->sli->data);

  for(i=0;i<seq->len-di->index_word_length;i++) {
    no = seq_number_dna_Nmer_noN(seq->seq+i,di->index_word_length);
    if( (*slci->is_populated)(slci->data,no) ) {
      slri = (*slci->lookup)(slci->data,no);
      for(;(*slri->is_more)(slri->data);) {
	res = (*slri->next)(slri->data,res);
	add_pair_HSPmanager(out,res->seq,i,res->pos);
      }

      free_SeqLookupResultInterface(slri);
    }
  }

  return out;
}


/* Function:  process_dna_HSP(set,query,query_pos,tseq,res_struct,mat)
 *
 * Descrip:    processes DNA based HSPs
 *
 *
 * Arg:               set [UNKN ] Undocumented argument [HSPset *]
 * Arg:             query [UNKN ] Undocumented argument [Sequence *]
 * Arg:         query_pos [UNKN ] Undocumented argument [int]
 * Arg:              tseq [UNKN ] Undocumented argument [Sequence *]
 * Arg:        res_struct [UNKN ] Undocumented argument [SeqLookupResultStruct *]
 * Arg:               mat [UNKN ] Undocumented argument [CompMat *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 167 "dnaindexcons.dy"
boolean process_dna_HSP(HSPset * set,Sequence * query,int query_pos,Sequence * tseq,SeqLookupResultStruct * res_struct,CompMat * mat)
{
  int k;

  for(k=0;k<set->len;k++) {
    if( on_HSP(set->hsp[k],query_pos,res_struct->pos) == TRUE ) {
      return FALSE;
    }
  }

  /* new HSP - extend and add */

  add_HSPset(set,new_HSP(NULL,query,tseq,query_pos,res_struct->pos,mat,5*Probability2Score(0.8/0.25)));

  return TRUE;
}


/* Function:  load_Sequence_DnaIndex(di,seq,sllp)
 *
 * Descrip:    Loads a sequence into a DnaIndex
 *
 *
 * Arg:          di [UNKN ] Undocumented argument [DnaIndex *]
 * Arg:         seq [UNKN ] Undocumented argument [Sequence *]
 * Arg:        sllp [UNKN ] Undocumented argument [SeqLookupLoadPara *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 188 "dnaindexcons.dy"
boolean load_Sequence_DnaIndex(DnaIndex * di,Sequence * seq,SeqLookupLoadPara * sllp)
{
  int i;
  int no;

  assert(di);
  assert(seq);
  assert(di->sli);
  assert(sllp);


  for(i=0;i<seq->len - di->index_word_length;i = i+sllp->tile_freq) {
    no = seq_number_dna_Nmer_noN(seq->seq+i,di->index_word_length);
    if( no != -1 ) {
      (*di->sli->add_direct_number)(di->sli->data,no,seq,i);
    }
  }
  

  return TRUE;
}


/* Function:  seq_number_dna_Nmer_noN(seq,index_length)
 *
 * Descrip:    General DNA index number generation
 *
 *
 * Arg:                 seq [UNKN ] Undocumented argument [char *]
 * Arg:        index_length [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 214 "dnaindexcons.dy"
int seq_number_dna_Nmer_noN(char * seq,int index_length)
{
  int i;
  int ret = 0;
  int base = 1;
  int no = 0;

  for(i=0;i<index_length;i++) {
    no = base_from_char(seq[i]);
    if( no == BASE_N ) {
      return -1;
    }
    
    ret += base * no;
    base = base * 4;
  }

  return ret;
    
}



# line 269 "dnaindexcons.c"
/* Function:  hard_link_DnaIndexConstructor(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DnaIndexConstructor *]
 *
 * Return [UNKN ]  Undocumented return value [DnaIndexConstructor *]
 *
 */
DnaIndexConstructor * hard_link_DnaIndexConstructor(DnaIndexConstructor * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a DnaIndexConstructor object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  DnaIndexConstructor_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaIndexConstructor *]
 *
 */
DnaIndexConstructor * DnaIndexConstructor_alloc(void) 
{
    DnaIndexConstructor * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(DnaIndexConstructor *) ckalloc (sizeof(DnaIndexConstructor))) == NULL)  {  
      warn("DnaIndexConstructor_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->type = DnaIndexConstructor_subseq;  
    out->index_word_length = 7;  
    out->match_prob = 0.8;   
    out->drop_off = 30;  


    return out;  
}    


/* Function:  free_DnaIndexConstructor(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DnaIndexConstructor *]
 *
 * Return [UNKN ]  Undocumented return value [DnaIndexConstructor *]
 *
 */
DnaIndexConstructor * free_DnaIndexConstructor(DnaIndexConstructor * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a DnaIndexConstructor obj. Should be trappable");   
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


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_DnaIndex(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DnaIndex *]
 *
 * Return [UNKN ]  Undocumented return value [DnaIndex *]
 *
 */
DnaIndex * hard_link_DnaIndex(DnaIndex * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a DnaIndex object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  DnaIndex_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaIndex *]
 *
 */
DnaIndex * DnaIndex_alloc(void) 
{
    DnaIndex * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(DnaIndex *) ckalloc (sizeof(DnaIndex))) == NULL)    {  
      warn("DnaIndex_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->index_word_length = 0;  
    out->sli = NULL; 
    out->cm = NULL;  
    out->drop_off = 30;  


    return out;  
}    


/* Function:  free_DnaIndex(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DnaIndex *]
 *
 * Return [UNKN ]  Undocumented return value [DnaIndex *]
 *
 */
DnaIndex * free_DnaIndex(DnaIndex * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a DnaIndex obj. Should be trappable");  
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
    if( obj->cm != NULL) 
      free_CompMat(obj->cm);     


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
