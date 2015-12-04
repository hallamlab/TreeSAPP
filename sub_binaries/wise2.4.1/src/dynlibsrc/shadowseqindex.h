#ifndef DYNAMITEshadowseqindexHEADERFILE
#define DYNAMITEshadowseqindexHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "seqlookup.h"
#include "shadowseq.h"
#include "genericindexresult.h"
#include "singlenumberspace.h"
#include "intallocator.h"
#include "staticseq.h"

#define SHADOW_ARRAYSEQL_BASIC  4
#define SHADOW_ARRAYSEQL_LINEAR 16

#define ShadowSequenceIndexLISTLENGTH 1024

#define SHADOW_TYPE int

typedef struct Wise2_ShadowArraySeqHead {
  SHADOW_TYPE * seqdb_pos;
  int current_pos;
  int max;
} ShadowArraySeqHead;

struct Wise2_ShadowSequenceIndex {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    ShadowArraySeqHead ** array;     
    int array_len;   
    ShadowSequence ** shadow;    
    int len;/* len for above shadow  */ 
    int maxlen; /* maxlen for above shadow */ 
    SingleNumberSpace * space;   
    int shadow_len;  
    int shadow_error;    
    IntAllocatorSet * ias;   
    StaticSeqHolder * ssh;   
    } ;  
/* ShadowSequenceIndex defined */ 
#ifndef DYNAMITE_DEFINED_ShadowSequenceIndex
typedef struct Wise2_ShadowSequenceIndex Wise2_ShadowSequenceIndex;
#define ShadowSequenceIndex Wise2_ShadowSequenceIndex
#define DYNAMITE_DEFINED_ShadowSequenceIndex
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  new_ShadowSequenceIndex_SeqLookupInterface(shadow_len,has_maxlen,maxlen,shadow_error)
 *
 * Descrip:    Provides a SeqLookupInterface, the common runtime plug-in for indexers
 *             using a ShadowSequenceIndex
 *
 *
 * Arg:          shadow_len [UNKN ] Undocumented argument [int]
 * Arg:          has_maxlen [UNKN ] Undocumented argument [int]
 * Arg:              maxlen [UNKN ] Undocumented argument [int]
 * Arg:        shadow_error [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupInterface *]
 *
 */
SeqLookupInterface * Wise2_new_ShadowSequenceIndex_SeqLookupInterface(int shadow_len,int has_maxlen,int maxlen,int shadow_error);
#define new_ShadowSequenceIndex_SeqLookupInterface Wise2_new_ShadowSequenceIndex_SeqLookupInterface


/* Function:  get_client_interface_ShadowSequenceIndex(data)
 *
 * Descrip:    gets client interface
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void*]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupClientInterface *]
 *
 */
SeqLookupClientInterface * Wise2_get_client_interface_ShadowSequenceIndex(void* data);
#define get_client_interface_ShadowSequenceIndex Wise2_get_client_interface_ShadowSequenceIndex


/* Function:  lookup_interface_ShadowSequenceClient(data,seq_number)
 *
 * Descrip:    For lookup interface, provides a result
 *
 *
 * Arg:              data [UNKN ] Undocumented argument [void *]
 * Arg:        seq_number [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupResultInterface *]
 *
 */
SeqLookupResultInterface * Wise2_lookup_interface_ShadowSequenceClient(void * data,int seq_number);
#define lookup_interface_ShadowSequenceClient Wise2_lookup_interface_ShadowSequenceClient


/* Function:  is_populated_interface_ShadowSequenceClient(data,seq_number)
 *
 * Descrip:    populated function for interface
 *
 *
 * Arg:              data [UNKN ] Undocumented argument [void *]
 * Arg:        seq_number [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_is_populated_interface_ShadowSequenceClient(void * data,int seq_number);
#define is_populated_interface_ShadowSequenceClient Wise2_is_populated_interface_ShadowSequenceClient


/* Function:  add_seq_interface_ShadowSequenceIndex(data,seq,para)
 *
 * Descrip:    add function for interface
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 * Arg:         seq [UNKN ] Undocumented argument [Sequence *]
 * Arg:        para [UNKN ] Undocumented argument [SeqLookupLoadPara *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_seq_interface_ShadowSequenceIndex(void * data,Sequence * seq,SeqLookupLoadPara * para);
#define add_seq_interface_ShadowSequenceIndex Wise2_add_seq_interface_ShadowSequenceIndex


/* Function:  free_interface_ShadowSequenceIndex(data)
 *
 * Descrip:    for interface, frees index
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 *
 */
void Wise2_free_interface_ShadowSequenceIndex(void * data);
#define free_interface_ShadowSequenceIndex Wise2_free_interface_ShadowSequenceIndex


/* Function:  free_interface_ShadowSequenceClient(data)
 *
 * Descrip:    Frees the client data
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 *
 */
void Wise2_free_interface_ShadowSequenceClient(void * data);
#define free_interface_ShadowSequenceClient Wise2_free_interface_ShadowSequenceClient


/* Function:  lookup_result_ShadowSeq(res,in,seq_no)
 *
 * Descrip:    handles the lookup and storage for a seq_no
 *             lookup
 *
 *
 * Arg:           res [UNKN ] Undocumented argument [GenericIndexResult *]
 * Arg:            in [UNKN ] Undocumented argument [ShadowSequenceIndex *]
 * Arg:        seq_no [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_lookup_result_ShadowSeq(GenericIndexResult * res,ShadowSequenceIndex * in,int seq_no);
#define lookup_result_ShadowSeq Wise2_lookup_result_ShadowSeq


/* Function:  add_result_GenericIndexResult_ShadowSeq(res,seq,pos)
 *
 * Descrip:    adds a particular shadow sequence position, unrolling
 *             shadowed sequences into the result
 *
 *
 * Arg:        res [UNKN ] Undocumented argument [GenericIndexResult *]
 * Arg:        seq [UNKN ] Undocumented argument [ShadowSequence *]
 * Arg:        pos [UNKN ] Undocumented argument [int]
 *
 */
void Wise2_add_result_GenericIndexResult_ShadowSeq(GenericIndexResult * res,ShadowSequence * seq,int pos);
#define add_result_GenericIndexResult_ShadowSeq Wise2_add_result_GenericIndexResult_ShadowSeq


/* Function:  dump_shadow_ShadowSequenceIndex(in,ofp)
 *
 * Descrip:    Dumps information about shadows
 *
 *
 * Arg:         in [UNKN ] Undocumented argument [ShadowSequenceIndex *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_dump_shadow_ShadowSequenceIndex(ShadowSequenceIndex * in,FILE * ofp);
#define dump_shadow_ShadowSequenceIndex Wise2_dump_shadow_ShadowSequenceIndex


/* Function:  dump_stats_ShadowSequenceIndex(in,ofp)
 *
 * Descrip:    Dumps useful information out of shadow sequence array
 *
 *
 * Arg:         in [UNKN ] Undocumented argument [ShadowSequenceIndex *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_dump_stats_ShadowSequenceIndex(ShadowSequenceIndex * in,FILE * ofp);
#define dump_stats_ShadowSequenceIndex Wise2_dump_stats_ShadowSequenceIndex


/* Function:  add_Sequence_ShadowSequenceIndex(in,seq,min_ext)
 *
 * Descrip:    Adds a Sequence to a ShadowIndex, placing shadowed regions
 *             correctly away
 *
 *
 * Arg:             in [UNKN ] Undocumented argument [ShadowSequenceIndex *]
 * Arg:            seq [UNKN ] Undocumented argument [Sequence *]
 * Arg:        min_ext [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_Sequence_ShadowSequenceIndex(ShadowSequenceIndex * in,Sequence * seq,int min_ext);
#define add_Sequence_ShadowSequenceIndex Wise2_add_Sequence_ShadowSequenceIndex


/* Function:  new_ShadowSequenceIndex(len,shadow_len,has_maxlen,maxlen,shadow_error)
 *
 * Descrip:    New ShadowSequenceIndex
 *
 *
 * Arg:                 len [UNKN ] Undocumented argument [int]
 * Arg:          shadow_len [UNKN ] Undocumented argument [int]
 * Arg:          has_maxlen [UNKN ] Undocumented argument [int]
 * Arg:              maxlen [UNKN ] Undocumented argument [int]
 * Arg:        shadow_error [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [ShadowSequenceIndex *]
 *
 */
ShadowSequenceIndex * Wise2_new_ShadowSequenceIndex(int len,int shadow_len,int has_maxlen,int maxlen,int shadow_error);
#define new_ShadowSequenceIndex Wise2_new_ShadowSequenceIndex


/* Function:  add_ShadowArraySeqHead(ias,h,seqdb_pos)
 *
 * Descrip:    Adds a sequence/pos pair to an ArrayHead
 *
 *
 * Arg:              ias [UNKN ] Undocumented argument [IntAllocatorSet *]
 * Arg:                h [UNKN ] Undocumented argument [ShadowArraySeqHead *]
 * Arg:        seqdb_pos [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_ShadowArraySeqHead(IntAllocatorSet * ias,ShadowArraySeqHead * h,int seqdb_pos);
#define add_ShadowArraySeqHead Wise2_add_ShadowArraySeqHead


/* Function:  new_ShadowArraySeqHead(ias)
 *
 * Descrip:    Builds a new ArraySeqHead structure
 *
 *
 * Arg:        ias [UNKN ] Undocumented argument [IntAllocatorSet *]
 *
 * Return [UNKN ]  Undocumented return value [ShadowArraySeqHead *]
 *
 */
ShadowArraySeqHead * Wise2_new_ShadowArraySeqHead(IntAllocatorSet * ias);
#define new_ShadowArraySeqHead Wise2_new_ShadowArraySeqHead


/* Function:  add_ShadowSequenceIndex(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [ShadowSequenceIndex *]
 * Arg:        add [OWNER] Object to add to the list [ShadowSequence *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_ShadowSequenceIndex(ShadowSequenceIndex * obj,ShadowSequence * add);
#define add_ShadowSequenceIndex Wise2_add_ShadowSequenceIndex


/* Function:  flush_ShadowSequenceIndex(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [ShadowSequenceIndex *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_ShadowSequenceIndex(ShadowSequenceIndex * obj);
#define flush_ShadowSequenceIndex Wise2_flush_ShadowSequenceIndex


/* Function:  ShadowSequenceIndex_alloc_std(void)
 *
 * Descrip:    Equivalent to ShadowSequenceIndex_alloc_len(ShadowSequenceIndexLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ShadowSequenceIndex *]
 *
 */
ShadowSequenceIndex * Wise2_ShadowSequenceIndex_alloc_std(void);
#define ShadowSequenceIndex_alloc_std Wise2_ShadowSequenceIndex_alloc_std


/* Function:  ShadowSequenceIndex_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [ShadowSequenceIndex *]
 *
 */
ShadowSequenceIndex * Wise2_ShadowSequenceIndex_alloc_len(int len);
#define ShadowSequenceIndex_alloc_len Wise2_ShadowSequenceIndex_alloc_len


/* Function:  hard_link_ShadowSequenceIndex(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ShadowSequenceIndex *]
 *
 * Return [UNKN ]  Undocumented return value [ShadowSequenceIndex *]
 *
 */
ShadowSequenceIndex * Wise2_hard_link_ShadowSequenceIndex(ShadowSequenceIndex * obj);
#define hard_link_ShadowSequenceIndex Wise2_hard_link_ShadowSequenceIndex


/* Function:  ShadowSequenceIndex_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ShadowSequenceIndex *]
 *
 */
ShadowSequenceIndex * Wise2_ShadowSequenceIndex_alloc(void);
#define ShadowSequenceIndex_alloc Wise2_ShadowSequenceIndex_alloc


/* Function:  free_ShadowSequenceIndex(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ShadowSequenceIndex *]
 *
 * Return [UNKN ]  Undocumented return value [ShadowSequenceIndex *]
 *
 */
ShadowSequenceIndex * Wise2_free_ShadowSequenceIndex(ShadowSequenceIndex * obj);
#define free_ShadowSequenceIndex Wise2_free_ShadowSequenceIndex


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_ShadowSequenceIndex(ShadowSequence ** list,int i,int j) ;
#define swap_ShadowSequenceIndex Wise2_swap_ShadowSequenceIndex
void Wise2_qsort_ShadowSequenceIndex(ShadowSequence ** list,int left,int right,int (*comp)(ShadowSequence * ,ShadowSequence * ));
#define qsort_ShadowSequenceIndex Wise2_qsort_ShadowSequenceIndex
void Wise2_sort_ShadowSequenceIndex(ShadowSequenceIndex * obj,int (*comp)(ShadowSequence *, ShadowSequence *));
#define sort_ShadowSequenceIndex Wise2_sort_ShadowSequenceIndex
boolean Wise2_expand_ShadowSequenceIndex(ShadowSequenceIndex * obj,int len);
#define expand_ShadowSequenceIndex Wise2_expand_ShadowSequenceIndex

#ifdef _cplusplus
}
#endif

#endif
