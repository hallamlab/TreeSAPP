#ifndef DYNAMITEcomplexevalsetHEADERFILE
#define DYNAMITEcomplexevalsetHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "complexsequence.h"
#include "genomic.h"



/*** ok, some types to works with ***/

#define AMINOACID CSEQ_PROTEIN_AMINOACID
#define CSEQ_PROTEIN_AMINOACID(cseq,index) (cseq->data[index])


/*** for genomic ***/

#define CSEQ_G_LEN 6

#define CSEQ_GENOMIC_BASE(cseq,index)      (cseq->data[(index)*CSEQ_G_LEN])
#define CSEQ_GENOMIC_CODON(cseq,index)     (cseq->data[(index)*CSEQ_G_LEN +1])
#define CSEQ_GENOMIC_5SS(cseq,index)       (cseq->data[(index)*CSEQ_G_LEN +2])
#define CSEQ_GENOMIC_3SS(cseq,index)       (cseq->data[(index)*CSEQ_G_LEN +3])
#define CSEQ_GENOMIC_REPEAT(cseq,index)    (cseq->data[(index)*CSEQ_G_LEN +4])
#define CSEQ_GENOMIC_CDSPOT(cseq,index)    (cseq->data[(index)*CSEQ_G_LEN +5])


/*** for cDNA  ***/

#define CSEQ_C_LEN 2

#define CSEQ_CDNA_BASE(cseq,index) (cseq->data[index*CSEQ_C_LEN])
#define CSEQ_CDNA_CODON(cseq,index) (cseq->data[index*CSEQ_C_LEN +1])

#define CSEQ_D_LEN 1

#define CSEQ_DNA_BASE(cseq,index) (cseq->data[index])



    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



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
ComplexSequence * Wise2_evaluate_ComplexSequence_Genomic(Genomic * gen,ComplexSequenceEvalSet * cses,int score_for_repeat,int score_for_cds_in_repeat);
#define evaluate_ComplexSequence_Genomic Wise2_evaluate_ComplexSequence_Genomic


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
ComplexSequenceEvalSet * Wise2_default_genomic_ComplexSequenceEvalSet(void);
#define default_genomic_ComplexSequenceEvalSet Wise2_default_genomic_ComplexSequenceEvalSet


/* Function:  flat_negi(void)
 *
 * Descrip:    Makes a ComplexSequenceEval that puts NEGI everywhere
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequenceEval *]
 *
 */
ComplexSequenceEval * Wise2_flat_negi(void);
#define flat_negi Wise2_flat_negi


/* Function:  flat_zero(void)
 *
 * Descrip:    Makes a ComplexSequenceEval that puts 0 everywhere
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequenceEval *]
 *
 */
ComplexSequenceEval * Wise2_flat_zero(void);
#define flat_zero Wise2_flat_zero


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
ComplexSequenceEval * Wise2_stupid_5SS(void);
#define stupid_5SS Wise2_stupid_5SS


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
ComplexSequenceEval * Wise2_stupid_3SS(void);
#define stupid_3SS Wise2_stupid_3SS


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
ComplexSequenceEvalSet * Wise2_default_DNA_ComplexSequenceEvalSet(void);
#define default_DNA_ComplexSequenceEvalSet Wise2_default_DNA_ComplexSequenceEvalSet


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
ComplexSequenceEvalSet * Wise2_default_cDNA_ComplexSequenceEvalSet(void);
#define default_cDNA_ComplexSequenceEvalSet Wise2_default_cDNA_ComplexSequenceEvalSet


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
ComplexSequenceEval * Wise2_codon_number_ComplexSequenceEval(void);
#define codon_number_ComplexSequenceEval Wise2_codon_number_ComplexSequenceEval


/* Function:  codon64_number_ComplexSequenceEval(void)
 *
 * Descrip:    Function which puts a 64 base codon...
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequenceEval *]
 *
 */
ComplexSequenceEval * Wise2_codon64_number_ComplexSequenceEval(void);
#define codon64_number_ComplexSequenceEval Wise2_codon64_number_ComplexSequenceEval


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
ComplexSequenceEval * Wise2_base_number_ComplexSequenceEval(void);
#define base_number_ComplexSequenceEval Wise2_base_number_ComplexSequenceEval


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
ComplexSequenceEvalSet * Wise2_default_aminoacid_ComplexSequenceEvalSet(void);
#define default_aminoacid_ComplexSequenceEvalSet Wise2_default_aminoacid_ComplexSequenceEvalSet


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
ComplexSequenceEval * Wise2_aminoacid_number_ComplexSequenceEval(void);
#define aminoacid_number_ComplexSequenceEval Wise2_aminoacid_number_ComplexSequenceEval


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
ComplexSequenceEvalSet * Wise2_default_dna_ComplexSequenceEvalSet(void);
#define default_dna_ComplexSequenceEvalSet Wise2_default_dna_ComplexSequenceEvalSet


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
int Wise2_flat_negi_eval(int type,void *data,char * seq);
#define flat_negi_eval Wise2_flat_negi_eval
int Wise2_flat_zero_eval(int type,void *data,char * seq);
#define flat_zero_eval Wise2_flat_zero_eval
int Wise2_stupid_5SS_eval(int type,void *data,char * seq);
#define stupid_5SS_eval Wise2_stupid_5SS_eval
int Wise2_stupid_3SS_eval(int type,void *data,char * seq);
#define stupid_3SS_eval Wise2_stupid_3SS_eval
int Wise2_codon_number_func(int type,void * data,char * seq);
#define codon_number_func Wise2_codon_number_func
int Wise2_codon64_number_func(int type,void * data,char * seq);
#define codon64_number_func Wise2_codon64_number_func
int Wise2_base_number_func(int type,void * data,char * seq);
#define base_number_func Wise2_base_number_func
int Wise2_aminoacid_number_func(int type,void * data,char * seq);
#define aminoacid_number_func Wise2_aminoacid_number_func

#ifdef _cplusplus
}
#endif

#endif
