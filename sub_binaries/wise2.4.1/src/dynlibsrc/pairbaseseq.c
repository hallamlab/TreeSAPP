#ifdef _cplusplus
extern "C" {
#endif
#include "pairbaseseq.h"


/* Function:  ComplexSequence_from_PairBaseSeq(pbs,splice5,splice3)
 *
 * Descrip:    Makes a ComplexSeq from a PairBaseSeq
 *
 *
 * Arg:            pbs [UNKN ] Undocumented argument [PairBaseSeq *]
 * Arg:        splice5 [UNKN ] Undocumented argument [ComplexSequenceEval *]
 * Arg:        splice3 [UNKN ] Undocumented argument [ComplexSequenceEval *]
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequence *]
 *
 */
# line 54 "pairbaseseq.dy"
ComplexSequence * ComplexSequence_from_PairBaseSeq(PairBaseSeq * pbs,ComplexSequenceEval * splice5,ComplexSequenceEval * splice3)
{
  int i;
  int * poi;
  ComplexSequence * out;
  Sequence * rev;
  pairbase_type a;
  pairbase_type g;
  pairbase_type t;

  

  assert(pbs);
  assert(pbs->anchor);
  assert(splice5);
  assert(splice3);

  a = MAKE_PAIRBASE(BASE_A,BASE_A);
  g = MAKE_PAIRBASE(BASE_G,BASE_G);
  t = MAKE_PAIRBASE(BASE_T,BASE_T);


  rev = reverse_complement_Sequence(pbs->anchor);

  out = ComplexSequence_alloc();
  if( out == NULL )
    return NULL;
  
  /*  fprintf(stderr,"Got %s as anchor\n",pbs->anchor->seq); */

  out->datastore = (int *) ckcalloc((pbs->anchor->len+10)*PAIRBASE_CS_LENGTH,sizeof(int));
  out->data  = out->datastore + (10 * PAIRBASE_CS_LENGTH);
  out->depth = PAIRBASE_CS_LENGTH;

  for(i=0;i<10;i++) {
    out->datastore[(i*PAIRBASE_CS_LENGTH)+0] = 4;
  }


  for(i=0,poi = out->data;i<pbs->anchor->len;i++,poi = next_ComplexSequence_data(out,poi)) {
    /* 0 is codon */
    poi[PAIRBASE_CS_CODON] = codon_number_func(0,NULL,pbs->anchor->seq+i);
    /* 1 is paircodon */
    poi[PAIRBASE_CS_PAIRCODON] = pairbase_codon_from_seq(pbs->seq+i-2);
    /* 2 is a pairbase */
    poi[PAIRBASE_CS_PAIRBASE] = pbs->seq[i];

    /* 3 is 5' SS */
    if( i < splice5->left_window || i+splice5->right_window > pbs->anchor->len || pbs->seq[i] != g || pbs->seq[i+1] != t) {
      poi[PAIRBASE_CS_5SS] = splice5->outside_score;
    } else {
      poi[PAIRBASE_CS_5SS] = splice5->eval_func(splice5->data_type,splice5->data,pbs->anchor->seq+i);
      poi[PAIRBASE_CS_5SS] += splice5->eval_func(splice5->data_type,splice5->data,pbs->sa->seq[1]->seq+i);
    }

    /* 4 is 3' SS */ 

    if( i < splice3->left_window || i+splice3->right_window > pbs->anchor->len || pbs->seq[i] != g || pbs->seq[i-1] != a ) {
      poi[PAIRBASE_CS_3SS] = splice3->outside_score;
    } else {
      poi[PAIRBASE_CS_3SS] = splice3->eval_func(splice3->data_type,splice3->data,pbs->anchor->seq+i);
      poi[PAIRBASE_CS_3SS] += splice3->eval_func(splice3->data_type,splice3->data,pbs->sa->seq[1]->seq+i);
    }

    /* 5 is reverse codon */
    poi[PAIRBASE_RS_CODON] = reverse_codon(codon_number_func(0,NULL,pbs->anchor->seq+i));

    /* 6 is reverse paircodon */
    poi[PAIRBASE_RS_PAIRCODON] = reverse_pairbase_codon(pairbase_codon_from_seq(pbs->seq+i-2));

    /* 7 is reverse pairbase */

    poi[PAIRBASE_RS_PAIRBASE] = complement_pairbase(pbs->seq[i]);

    /* 8 is reverse 5'SS */

    if( i < splice5->left_window || i+splice5->right_window > pbs->anchor->len ) {
      poi[PAIRBASE_RS_5SS] = splice5->outside_score;
    } else {
      poi[PAIRBASE_RS_5SS] = splice5->eval_func(splice5->data_type,splice5->data,rev->seq+rev->len-1-i);
    }

    /* 9 is reverse 3'SS */

    if( i < splice3->left_window || i+splice3->right_window > pbs->anchor->len ) {
      poi[PAIRBASE_RS_3SS] = splice3->outside_score;
    } else {
      poi[PAIRBASE_RS_3SS] = splice3->eval_func(splice3->data_type,splice3->data,rev->seq+rev->len-1-i);
    }

/*    fprintf(stdout,"Position %d got %d %d \n",i,poi[PAIRBASE_RS_5SS],poi[PAIRBASE_RS_3SS]); */

  }

  out->length = pbs->anchor->len;
  out->seq = hard_link_Sequence(pbs->anchor);
  
  free_Sequence(rev);

  return out;

}


/* Function:  SeqAlign_from_PairBaseSeq(pbs)
 *
 * Descrip:    Makes a SeqAlign from a PairBaseSeq
 *
 *
 * Arg:        pbs [UNKN ] Undocumented argument [PairBaseSeq *]
 *
 * Return [UNKN ]  Undocumented return value [SeqAlign *]
 *
 */
# line 161 "pairbaseseq.dy"
SeqAlign * SeqAlign_from_PairBaseSeq(PairBaseSeq * pbs)
{
  SeqAlign * out;
  Sequence * a;
  Sequence * b;
  int i;

  assert(pbs);



  out = SeqAlign_alloc_len(2);

  
  a = Sequence_alloc_len(pbs->len);
  b = Sequence_alloc_len(pbs->len);

  add_SeqAlign(out,a);
  add_SeqAlign(out,b);

  for(i=0;i<pbs->len;i++) {
    a->seq[i] = char_from_base(anchor_base_from_pairbase(pbs->seq[i]));
    b->seq[i] = char_from_base(informant_base_from_pairbase(pbs->seq[i]));
  }

  a->seq[i] = b->seq[i] = '\0';

  return out;
}

/* Function:  new_PairBaseSeq_SeqAlign(al)
 *
 * Descrip:    Makes a pairseq from a SeqAlign
 *
 *
 * Arg:        al [UNKN ] Undocumented argument [SeqAlign *]
 *
 * Return [UNKN ]  Undocumented return value [PairBaseSeq *]
 *
 */
# line 194 "pairbaseseq.dy"
PairBaseSeq * new_PairBaseSeq_SeqAlign(SeqAlign * al)
{
  PairBaseSeq * out;

  out = new_PairBaseSeq_strings(al->seq[0]->seq,al->seq[1]->seq);
  out->sa = hard_link_SeqAlign(al);
  
  return out;
}

/* Function:  new_PairBaseSeq_strings(a,b)
 *
 * Descrip:    Makes a pairseq from two strings
 *
 *
 * Arg:        a [UNKN ] Undocumented argument [char *]
 * Arg:        b [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [PairBaseSeq *]
 *
 */
# line 207 "pairbaseseq.dy"
PairBaseSeq * new_PairBaseSeq_strings(char * a,char * b)
{
  int i;
  int l;
  PairBaseSeq * out;
  Sequence * anchor;

  assert(a);
  assert(b);
  
  if( strlen(a) != strlen(b) ) {
    warn("Cannot build PairBaseSeq with different length strings %d %d",strlen(a),strlen(b));
    return NULL;
  }

  l = strlen(a);

  out = PairBaseSeq_alloc();

  out->seq = (pairbase_type *) calloc(l,sizeof(pairbase_type));

  for(i=0;i<l;i++) {
    out->seq[i] = MAKE_PAIRBASE(base_from_char(a[i]),base_from_char(b[i]));
  }
  anchor = new_Sequence_from_strings("anchor",a);
  out->anchor = anchor;
  out->len = l;
  return out;
}
  
# line 219 "pairbaseseq.c"
/* Function:  hard_link_PairBaseSeq(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [PairBaseSeq *]
 *
 * Return [UNKN ]  Undocumented return value [PairBaseSeq *]
 *
 */
PairBaseSeq * hard_link_PairBaseSeq(PairBaseSeq * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a PairBaseSeq object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  PairBaseSeq_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PairBaseSeq *]
 *
 */
PairBaseSeq * PairBaseSeq_alloc(void) 
{
    PairBaseSeq * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(PairBaseSeq *) ckalloc (sizeof(PairBaseSeq))) == NULL)  {  
      warn("PairBaseSeq_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->seq = NULL; 
    out->len = 0;    
    out->anchor = NULL;  
    out->sa = NULL;  


    return out;  
}    


/* Function:  free_PairBaseSeq(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [PairBaseSeq *]
 *
 * Return [UNKN ]  Undocumented return value [PairBaseSeq *]
 *
 */
PairBaseSeq * free_PairBaseSeq(PairBaseSeq * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a PairBaseSeq obj. Should be trappable");   
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
    if( obj->seq != NULL)    
      ckfree(obj->seq);  
    if( obj->anchor != NULL) 
      free_Sequence(obj->anchor);    
    if( obj->sa != NULL) 
      free_SeqAlign(obj->sa);    


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
