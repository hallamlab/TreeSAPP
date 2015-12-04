#ifndef DYNAMITEpairbaseseqHEADERFILE
#define DYNAMITEpairbaseseqHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "pairbase.h"
#include "seqalign.h"
#include "complexsequence.h"
#include "complexevalset.h"


#define PAIRBASE_CS_CODON     0
#define PAIRBASE_CS_PAIRCODON 1
#define PAIRBASE_CS_PAIRBASE  2
#define PAIRBASE_CS_5SS       3
#define PAIRBASE_CS_3SS       4
#define PAIRBASE_RS_CODON     5
#define PAIRBASE_RS_PAIRCODON 6
#define PAIRBASE_RS_PAIRBASE  7
#define PAIRBASE_RS_5SS       8
#define PAIRBASE_RS_3SS       9

#define PAIRBASE_CS_LENGTH    10

#define CSEQ_PAIR_CODON(cseq,index) (cseq->data[index*PAIRBASE_CS_LENGTH+PAIRBASE_CS_CODON])
#define CSEQ_PAIR_PAIRCODON(cseq,index) (cseq->data[index*PAIRBASE_CS_LENGTH+PAIRBASE_CS_PAIRCODON])
#define CSEQ_PAIR_PAIRBASE(cseq,index) (cseq->data[index*PAIRBASE_CS_LENGTH+PAIRBASE_CS_PAIRBASE])
#define CSEQ_PAIR_5SS(cseq,index) (cseq->data[index*PAIRBASE_CS_LENGTH+PAIRBASE_CS_5SS])
#define CSEQ_PAIR_3SS(cseq,index) (cseq->data[index*PAIRBASE_CS_LENGTH+PAIRBASE_CS_3SS])
#define CSEQ_REV_PAIR_CODON(cseq,index) (cseq->data[index*PAIRBASE_CS_LENGTH+PAIRBASE_RS_CODON])
#define CSEQ_REV_PAIR_PAIRCODON(cseq,index) (cseq->data[index*PAIRBASE_CS_LENGTH+PAIRBASE_RS_PAIRCODON])
#define CSEQ_REV_PAIR_PAIRBASE(cseq,index) (cseq->data[index*PAIRBASE_CS_LENGTH+PAIRBASE_RS_PAIRBASE])
#define CSEQ_REV_PAIR_5SS(cseq,index) (cseq->data[index*PAIRBASE_CS_LENGTH+PAIRBASE_RS_5SS])
#define CSEQ_REV_PAIR_3SS(cseq,index) (cseq->data[index*PAIRBASE_CS_LENGTH+PAIRBASE_RS_3SS])



struct Wise2_PairBaseSeq {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * seq;  
    int len;     
    Sequence * anchor;  /*  may not always be here */ 
    SeqAlign * sa;  /*  may not always be here */ 
    } ;  
/* PairBaseSeq defined */ 
#ifndef DYNAMITE_DEFINED_PairBaseSeq
typedef struct Wise2_PairBaseSeq Wise2_PairBaseSeq;
#define PairBaseSeq Wise2_PairBaseSeq
#define DYNAMITE_DEFINED_PairBaseSeq
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



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
ComplexSequence * Wise2_ComplexSequence_from_PairBaseSeq(PairBaseSeq * pbs,ComplexSequenceEval * splice5,ComplexSequenceEval * splice3);
#define ComplexSequence_from_PairBaseSeq Wise2_ComplexSequence_from_PairBaseSeq


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
SeqAlign * Wise2_SeqAlign_from_PairBaseSeq(PairBaseSeq * pbs);
#define SeqAlign_from_PairBaseSeq Wise2_SeqAlign_from_PairBaseSeq


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
PairBaseSeq * Wise2_new_PairBaseSeq_SeqAlign(SeqAlign * al);
#define new_PairBaseSeq_SeqAlign Wise2_new_PairBaseSeq_SeqAlign


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
PairBaseSeq * Wise2_new_PairBaseSeq_strings(char * a,char * b);
#define new_PairBaseSeq_strings Wise2_new_PairBaseSeq_strings


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
PairBaseSeq * Wise2_hard_link_PairBaseSeq(PairBaseSeq * obj);
#define hard_link_PairBaseSeq Wise2_hard_link_PairBaseSeq


/* Function:  PairBaseSeq_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PairBaseSeq *]
 *
 */
PairBaseSeq * Wise2_PairBaseSeq_alloc(void);
#define PairBaseSeq_alloc Wise2_PairBaseSeq_alloc


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
PairBaseSeq * Wise2_free_PairBaseSeq(PairBaseSeq * obj);
#define free_PairBaseSeq Wise2_free_PairBaseSeq


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/

#ifdef _cplusplus
}
#endif

#endif
