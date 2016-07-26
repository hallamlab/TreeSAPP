#ifdef _cplusplus
extern "C" {
#endif
#include "complexevalset.h"


/* Function:  evaluate_ComplexSequence_Genomic(gen,cses,score_for_repeat,score_for_cds_in_repeat)
 *
 * Descrip:    makes a new ComplexSequence from a Genomic
 *             sequence, using the repeat information to
 *             structure the positives and negatives.
 *
 *
 * Arg:                            gen [UNKN ] Undocumented argument [Genomic *]
 * Arg:                           cses [UNKN ] Undocumented argument [ComplexSequenceEvalSet *]
 * Arg:               score_for_repeat [UNKN ] Undocumented argument [int]
 * Arg:        score_for_cds_in_repeat [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequence *]
 *
 */
# line 55 "complexevalset.dy"
ComplexSequence * evaluate_ComplexSequence_Genomic(Genomic * gen,ComplexSequenceEvalSet * cses,int score_for_repeat,int score_for_cds_in_repeat)
{
  int i,j;
  ComplexSequence * out;

  if( cses == NULL || cses->type != SEQUENCE_GENOMIC ) {
    warn("ComplexSequenceEvalSet is not valid for genomic!");
    return NULL;
  }

  assert(gen);
  assert(gen->baseseq);
  assert(gen->baseseq->seq);

  out = new_ComplexSequence(gen->baseseq,cses);


  if( out == NULL ) {
    return NULL;
  }
  
  for(i=0;i<gen->len;i++) {
    for(j=gen->repeat[i]->start;j<gen->repeat[i]->end;j++) {
      if( j >= gen->baseseq->len ) {
	warn("Overran sequence - j position %d, sequence length %d on repeat %d",j,gen->baseseq->len,i);
	break;
      }

      CSEQ_GENOMIC_REPEAT(out,j) = score_for_repeat;
      CSEQ_GENOMIC_CDSPOT(out,j) = score_for_cds_in_repeat;
    }
  }
  /*
  for(i=0;i<gen->baseseq->len;i++) 
    fprintf(stderr,"At position %6d %d\n",CSEQ_GENOMIC_CDSPOT(out,j));
  */

  return out;
}


/* Function:  default_genomic_ComplexSequenceEvalSet(void)
 *
 * Descrip:    Makes a reasonably sensible genomic sequence
 *             eval set. Has base and codon numbers (what
 *             every good genomic sequence should do) and then
 *             /stupid_5SS and /stupid_3SS. Finally the repeat/EST
 *             regions modelled with the /flat_zero
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequenceEvalSet *]
 *
 */
# line 104 "complexevalset.dy"
ComplexSequenceEvalSet * default_genomic_ComplexSequenceEvalSet(void)
{
  ComplexSequenceEvalSet * out;


  out = ComplexSequenceEvalSet_alloc_len(11);

  add_ComplexSequenceEvalSet(out,base_number_ComplexSequenceEval());
  add_ComplexSequenceEvalSet(out,codon_number_ComplexSequenceEval());
  add_ComplexSequenceEvalSet(out,stupid_5SS());
  add_ComplexSequenceEvalSet(out,stupid_3SS());
  add_ComplexSequenceEvalSet(out,flat_zero());
  add_ComplexSequenceEvalSet(out,flat_zero());


  out->type = SEQUENCE_GENOMIC;

  prepare_ComplexSequenceEvalSet(out);

  return out;
}

/* Function:  flat_negi(void)
 *
 * Descrip:    Makes a ComplexSequenceEval that puts NEGI everywhere
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequenceEval *]
 *
 */
# line 129 "complexevalset.dy"
ComplexSequenceEval * flat_negi(void)
{
  ComplexSequenceEval * out;
  out= ComplexSequenceEval_alloc();

  out->left_window   = 0;
  out->right_window  = 0;
  out->outside_score = NEGI;
  out->eval_func     = flat_negi_eval;

  return out;
}

/* Function:  flat_negi_eval(type,*data,seq)
 *
 * Descrip:    used by flat negi
 *
 *
 * Arg:         type [UNKN ] Undocumented argument [int]
 * Arg:        *data [UNKN ] Undocumented argument [void]
 * Arg:          seq [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 146 "complexevalset.dy"
int flat_negi_eval(int type,void *data,char * seq)
{
  return NEGI;
}
    

/* Function:  flat_zero(void)
 *
 * Descrip:    Makes a ComplexSequenceEval that puts 0 everywhere
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequenceEval *]
 *
 */
# line 156 "complexevalset.dy"
ComplexSequenceEval * flat_zero(void)
{
  ComplexSequenceEval * out;
  
  out= ComplexSequenceEval_alloc();

  out->left_window   = 0;
  out->right_window  = 0;
  out->outside_score = NEGI;
  out->eval_func     = flat_zero_eval;

  return out;
}

/* Function:  flat_zero_eval(type,*data,seq)
 *
 * Descrip:    Used by /flat_zero as function actually called
 *
 *
 * Arg:         type [UNKN ] Undocumented argument [int]
 * Arg:        *data [UNKN ] Undocumented argument [void]
 * Arg:          seq [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 175 "complexevalset.dy"
int flat_zero_eval(int type,void *data,char * seq)
{
  return 0;
}


/* Function:  stupid_5SS(void)
 *
 * Descrip:    Reasonably stupid 5'SS, of 0 at GT's
 *             and NEGI elsewhere. Pretends to have a longer
 *             footprint than it needs to satisify the lookbacks
 *             of more proper genomic models
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequenceEval *]
 *
 */
# line 188 "complexevalset.dy"
ComplexSequenceEval * stupid_5SS(void)
{
  ComplexSequenceEval * out;
  
  out= ComplexSequenceEval_alloc();

  out->left_window   = 3;
  out->right_window  = 7;
  out->left_lookback = 8;
  out->outside_score = NEGI;
  out->eval_func     = stupid_5SS_eval;

  return out;
}

/* Function:  stupid_5SS_eval(type,*data,seq)
 *
 * Descrip:    Function which actually does the evaluation for
 *             5'SS
 *
 *
 * Arg:         type [UNKN ] Undocumented argument [int]
 * Arg:        *data [UNKN ] Undocumented argument [void]
 * Arg:          seq [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 209 "complexevalset.dy"
int stupid_5SS_eval(int type,void *data,char * seq)
{
  if( *(seq) == 'G' && *(seq+1) == 'T' ) 
    return 0;
  else return NEGI;
}


/* Function:  stupid_3SS(void)
 *
 * Descrip:    Reasonably stupid 3'SS, of 0 at AG's
 *             and NEGI elsewhere. Pretends to have a longer
 *             footprint than it needs to satisify the lookbacks
 *             of more proper genomic models
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequenceEval *]
 *
 */
# line 224 "complexevalset.dy"
ComplexSequenceEval * stupid_3SS(void)
{
  ComplexSequenceEval * out;
  
  out= ComplexSequenceEval_alloc();

  out->left_window   = 3;
  out->right_window  = 3;
  out->left_lookback = 5;
  out->outside_score = NEGI;
  out->eval_func     = stupid_3SS_eval;

  return out;
}


/* Function:  stupid_3SS_eval(type,*data,seq)
 *
 * Descrip:    Function which actually does the evaluation for
 *             3'SS
 *
 *
 * Arg:         type [UNKN ] Undocumented argument [int]
 * Arg:        *data [UNKN ] Undocumented argument [void]
 * Arg:          seq [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 246 "complexevalset.dy"
int stupid_3SS_eval(int type,void *data,char * seq)
{
  if( *(seq-1) == 'A' && *(seq) == 'G' ) 
    return 0;
  else return NEGI;
}



   /********************************/
   /* cDNA and base stuff          */
   /********************************/

/* Function:  default_DNA_ComplexSequenceEvalSet(void)
 *
 * Descrip:    Makes a very sensible DNA sequence
 *             eval set. You shouldn't need your own
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequenceEvalSet *]
 *
 */
# line 264 "complexevalset.dy"
ComplexSequenceEvalSet * default_DNA_ComplexSequenceEvalSet(void)
{
  ComplexSequenceEvalSet * out;


  out = ComplexSequenceEvalSet_alloc_len(1);

  add_ComplexSequenceEvalSet(out,base_number_ComplexSequenceEval());
  out->type = SEQUENCE_DNA;

  prepare_ComplexSequenceEvalSet(out);

  return out;
}


/* Function:  default_cDNA_ComplexSequenceEvalSet(void)
 *
 * Descrip:    Makes a very sensible cDNA sequence
 *             eval set. You shouldn't need your own
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequenceEvalSet *]
 *
 */
# line 285 "complexevalset.dy"
ComplexSequenceEvalSet * default_cDNA_ComplexSequenceEvalSet(void)
{
  ComplexSequenceEvalSet * out;


  out = ComplexSequenceEvalSet_alloc_len(2);

  add_ComplexSequenceEvalSet(out,base_number_ComplexSequenceEval());
  add_ComplexSequenceEvalSet(out,codon_number_ComplexSequenceEval());
  out->type = SEQUENCE_CDNA;

  prepare_ComplexSequenceEvalSet(out);

  return out;
}




/* Function:  codon_number_ComplexSequenceEval(void)
 *
 * Descrip:    ComplexSequenceEval that puts a codon
 *             number in there (0-125 inc, 125 = impossible codon).
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequenceEval *]
 *
 */
# line 309 "complexevalset.dy"
ComplexSequenceEval * codon_number_ComplexSequenceEval(void)
{
  ComplexSequenceEval * out;


  out = ComplexSequenceEval_alloc();

  if (out == NULL ) {
    return NULL;
  }
  
  out->left_window = 2;
  out->right_window = 0;
  out->outside_score = 125;
  out->eval_func = codon_number_func;


  return out;
}

/* Function:  codon_number_func(type,data,seq)
 *
 * Descrip:    Function which actually does the codon evaluation
 *
 *
 * Arg:        type [UNKN ] Undocumented argument [int]
 * Arg:        data [UNKN ] Undocumented argument [void *]
 * Arg:         seq [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 334 "complexevalset.dy"
int codon_number_func(int type,void * data,char * seq)
{
  return codon_from_seq(seq-2);
} 

  
/* Function:  codon64_number_ComplexSequenceEval(void)
 *
 * Descrip:    Function which puts a 64 base codon...
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequenceEval *]
 *
 */
# line 343 "complexevalset.dy"
ComplexSequenceEval * codon64_number_ComplexSequenceEval(void)
{
  ComplexSequenceEval * out;


  out = ComplexSequenceEval_alloc();

  if (out == NULL ) {
    return NULL;
  }
  
  out->left_window = 2;
  out->right_window = 0;
  out->outside_score = 65;
  out->eval_func = codon64_number_func;


  return out;
}

/* Function:  codon64_number_func(type,data,seq)
 *
 * Descrip:    Function which actually does the codon64 evaluation
 *
 *
 * Arg:        type [UNKN ] Undocumented argument [int]
 * Arg:        data [UNKN ] Undocumented argument [void *]
 * Arg:         seq [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 368 "complexevalset.dy"
int codon64_number_func(int type,void * data,char * seq)
{
  int codon;

  if( !is_non_ambiguous_codon_seq(seq-2) ) {
      	return 64;
  }
  codon = codon_from_seq(seq-2);

  return base4_codon_from_codon(codon);
} 



/* Function:  base_number_ComplexSequenceEval(void)
 *
 * Descrip:    ComplexSequenceEval that puts a base
 *             number (0-5 inc)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequenceEval *]
 *
 */
# line 387 "complexevalset.dy"
ComplexSequenceEval * base_number_ComplexSequenceEval(void)
{
  ComplexSequenceEval * out;


  out = ComplexSequenceEval_alloc();

  if (out == NULL ) {
    return NULL;
  }
  
  out->left_window = out->right_window = 0;
  out->eval_func = base_number_func;

  return out;
}
  

/* Function:  base_number_func(type,data,seq)
 *
 * Descrip:    Function which actually does the evaluation bases
 *
 *
 * Arg:        type [UNKN ] Undocumented argument [int]
 * Arg:        data [UNKN ] Undocumented argument [void *]
 * Arg:         seq [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 410 "complexevalset.dy"
int base_number_func(int type,void * data,char * seq)
{
  return base_from_char(*seq);
} 



  /***************************/
  /* amino acid stuff        */
  /***************************/


/* Function:  default_aminoacid_ComplexSequenceEvalSet(void)
 *
 * Descrip:    Makes a very sensible protein sequence
 *             eval set. You shouldn't need your own
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequenceEvalSet *]
 *
 */
# line 427 "complexevalset.dy"
ComplexSequenceEvalSet * default_aminoacid_ComplexSequenceEvalSet(void)
{
  ComplexSequenceEvalSet * out;


  out = ComplexSequenceEvalSet_alloc_len(1);

  add_ComplexSequenceEvalSet(out,aminoacid_number_ComplexSequenceEval());

  out->type = SEQUENCE_PROTEIN;

  prepare_ComplexSequenceEvalSet(out);

  return out;
}



/* Function:  aminoacid_number_ComplexSequenceEval(void)
 *
 * Descrip:    ComplexSequenceEval that puts a amino acid
 *             number in there (0-26)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequenceEval *]
 *
 */
# line 450 "complexevalset.dy"
ComplexSequenceEval * aminoacid_number_ComplexSequenceEval(void)
{
  ComplexSequenceEval * out;


  out = ComplexSequenceEval_alloc();

  if (out == NULL ) {
    return NULL;
  }
  
  out->left_window = out->right_window = 0;
  out->eval_func = aminoacid_number_func;

  return out;
}
  

/* Function:  aminoacid_number_func(type,data,seq)
 *
 * Descrip:    Function which actually does the evaluation for aminoacids
 *
 *
 * Arg:        type [UNKN ] Undocumented argument [int]
 * Arg:        data [UNKN ] Undocumented argument [void *]
 * Arg:         seq [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 473 "complexevalset.dy"
int aminoacid_number_func(int type,void * data,char * seq)
{
  return (int)(toupper((int)*seq)-'A');
} 


  /***************************/
  /* dna stuff               */
  /***************************/


/* Function:  default_dna_ComplexSequenceEvalSet(void)
 *
 * Descrip:    Makes a very sensible dna sequence
 *             eval set. You shouldn't need your own
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequenceEvalSet *]
 *
 */
# line 489 "complexevalset.dy"
ComplexSequenceEvalSet * default_dna_ComplexSequenceEvalSet(void)
{
  ComplexSequenceEvalSet * out;

  out = ComplexSequenceEvalSet_alloc_len(1);

  add_ComplexSequenceEvalSet(out,base_number_ComplexSequenceEval());

  out->type = SEQUENCE_DNA;

  prepare_ComplexSequenceEvalSet(out);

  return out;
}

# line 572 "complexevalset.c"

#ifdef _cplusplus
}
#endif
