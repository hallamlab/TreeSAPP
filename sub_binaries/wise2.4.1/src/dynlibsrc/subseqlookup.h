#ifndef DYNAMITEsubseqlookupHEADERFILE
#define DYNAMITEsubseqlookupHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "sequence.h"
#include "seqlookup.h"
#include "linkedlist_lookpos.h"

#define LOOKUP_BLOCK_SIZE 1024*1024
#define SeqLookupPosBlockAllocatorLISTLENGTH 2048

struct Wise2_SeqLookupPosBlock {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    SeqLookupPos block[LOOKUP_BLOCK_SIZE];   
    } ;  
/* SeqLookupPosBlock defined */ 
#ifndef DYNAMITE_DEFINED_SeqLookupPosBlock
typedef struct Wise2_SeqLookupPosBlock Wise2_SeqLookupPosBlock;
#define SeqLookupPosBlock Wise2_SeqLookupPosBlock
#define DYNAMITE_DEFINED_SeqLookupPosBlock
#endif


struct Wise2_SeqLookupPosBlockAllocator {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    SeqLookupPosBlock ** block;  
    int len;/* len for above block  */ 
    int maxlen; /* maxlen for above block */ 
    int pos;     
    } ;  
/* SeqLookupPosBlockAllocator defined */ 
#ifndef DYNAMITE_DEFINED_SeqLookupPosBlockAllocator
typedef struct Wise2_SeqLookupPosBlockAllocator Wise2_SeqLookupPosBlockAllocator;
#define SeqLookupPosBlockAllocator Wise2_SeqLookupPosBlockAllocator
#define DYNAMITE_DEFINED_SeqLookupPosBlockAllocator
#endif


struct Wise2_SubSeqLookup {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    SeqLookupPos ** array;   
    int array_max;   
    SeqLookupPosBlockAllocator * block;  
    } ;  
/* SubSeqLookup defined */ 
#ifndef DYNAMITE_DEFINED_SubSeqLookup
typedef struct Wise2_SubSeqLookup Wise2_SubSeqLookup;
#define SubSeqLookup Wise2_SubSeqLookup
#define DYNAMITE_DEFINED_SubSeqLookup
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  new_array_SeqLookupInterface(array_max_size)
 *
 * Descrip:    Makes a new array Lookup system
 *
 *
 * Arg:        array_max_size [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupInterface *]
 *
 */
SeqLookupInterface * Wise2_new_array_SeqLookupInterface(int array_max_size);
#define new_array_SeqLookupInterface Wise2_new_array_SeqLookupInterface


/* Function:  free_subseqlookup(data)
 *
 * Descrip:    free function for the hash
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 *
 */
void Wise2_free_subseqlookup(void * data);
#define free_subseqlookup Wise2_free_subseqlookup


/* Function:  is_populated_subseqlookup(data,seq_number)
 *
 * Descrip:    tells whether this is populated or not
 *
 *
 * Arg:              data [UNKN ] Undocumented argument [void *]
 * Arg:        seq_number [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_is_populated_subseqlookup(void * data, int seq_number);
#define is_populated_subseqlookup Wise2_is_populated_subseqlookup


/* Function:  lookup_subseqlookup(data,seq_number)
 *
 * Descrip:    Retrieves a SeqLookup position 
 *
 *
 * Arg:              data [UNKN ] Undocumented argument [void *]
 * Arg:        seq_number [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupResultInterface *]
 *
 */
SeqLookupResultInterface * Wise2_lookup_subseqlookup(void * data, int seq_number);
#define lookup_subseqlookup Wise2_lookup_subseqlookup


/* Function:  add_subseqlookup(data,seq_number,seq,pos)
 *
 * Descrip:    Adds a sequence/pos pair to the hash
 *
 *
 * Arg:              data [UNKN ] Undocumented argument [void *]
 * Arg:        seq_number [UNKN ] Undocumented argument [int]
 * Arg:               seq [UNKN ] Undocumented argument [Sequence *]
 * Arg:               pos [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_subseqlookup(void * data,int seq_number,Sequence * seq,int pos);
#define add_subseqlookup Wise2_add_subseqlookup


/* Function:  new_SeqLookupPosBlockAllocator(void)
 *
 * Descrip:    Makes a new SeqPosLookup Block Allocator
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupPosBlockAllocator *]
 *
 */
SeqLookupPosBlockAllocator * Wise2_new_SeqLookupPosBlockAllocator(void);
#define new_SeqLookupPosBlockAllocator Wise2_new_SeqLookupPosBlockAllocator


/* Function:  new_SeqLookupPos_BlockAllocator(*bla)
 *
 * Descrip:    Returns a new SeqPosLookup
 *
 *
 * Arg:        *bla [UNKN ] Undocumented argument [SeqLookupPosBlockAllocator]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupPos *]
 *
 */
SeqLookupPos * Wise2_new_SeqLookupPos_BlockAllocator(SeqLookupPosBlockAllocator *bla);
#define new_SeqLookupPos_BlockAllocator Wise2_new_SeqLookupPos_BlockAllocator


/* Function:  hard_link_SeqLookupPosBlock(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SeqLookupPosBlock *]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupPosBlock *]
 *
 */
SeqLookupPosBlock * Wise2_hard_link_SeqLookupPosBlock(SeqLookupPosBlock * obj);
#define hard_link_SeqLookupPosBlock Wise2_hard_link_SeqLookupPosBlock


/* Function:  SeqLookupPosBlock_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupPosBlock *]
 *
 */
SeqLookupPosBlock * Wise2_SeqLookupPosBlock_alloc(void);
#define SeqLookupPosBlock_alloc Wise2_SeqLookupPosBlock_alloc


/* Function:  free_SeqLookupPosBlock(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SeqLookupPosBlock *]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupPosBlock *]
 *
 */
SeqLookupPosBlock * Wise2_free_SeqLookupPosBlock(SeqLookupPosBlock * obj);
#define free_SeqLookupPosBlock Wise2_free_SeqLookupPosBlock


/* Function:  add_SeqLookupPosBlockAllocator(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [SeqLookupPosBlockAllocator *]
 * Arg:        add [OWNER] Object to add to the list [SeqLookupPosBlock *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_SeqLookupPosBlockAllocator(SeqLookupPosBlockAllocator * obj,SeqLookupPosBlock * add);
#define add_SeqLookupPosBlockAllocator Wise2_add_SeqLookupPosBlockAllocator


/* Function:  flush_SeqLookupPosBlockAllocator(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [SeqLookupPosBlockAllocator *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_SeqLookupPosBlockAllocator(SeqLookupPosBlockAllocator * obj);
#define flush_SeqLookupPosBlockAllocator Wise2_flush_SeqLookupPosBlockAllocator


/* Function:  SeqLookupPosBlockAllocator_alloc_std(void)
 *
 * Descrip:    Equivalent to SeqLookupPosBlockAllocator_alloc_len(SeqLookupPosBlockAllocatorLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupPosBlockAllocator *]
 *
 */
SeqLookupPosBlockAllocator * Wise2_SeqLookupPosBlockAllocator_alloc_std(void);
#define SeqLookupPosBlockAllocator_alloc_std Wise2_SeqLookupPosBlockAllocator_alloc_std


/* Function:  SeqLookupPosBlockAllocator_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupPosBlockAllocator *]
 *
 */
SeqLookupPosBlockAllocator * Wise2_SeqLookupPosBlockAllocator_alloc_len(int len);
#define SeqLookupPosBlockAllocator_alloc_len Wise2_SeqLookupPosBlockAllocator_alloc_len


/* Function:  hard_link_SeqLookupPosBlockAllocator(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SeqLookupPosBlockAllocator *]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupPosBlockAllocator *]
 *
 */
SeqLookupPosBlockAllocator * Wise2_hard_link_SeqLookupPosBlockAllocator(SeqLookupPosBlockAllocator * obj);
#define hard_link_SeqLookupPosBlockAllocator Wise2_hard_link_SeqLookupPosBlockAllocator


/* Function:  SeqLookupPosBlockAllocator_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupPosBlockAllocator *]
 *
 */
SeqLookupPosBlockAllocator * Wise2_SeqLookupPosBlockAllocator_alloc(void);
#define SeqLookupPosBlockAllocator_alloc Wise2_SeqLookupPosBlockAllocator_alloc


/* Function:  free_SeqLookupPosBlockAllocator(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SeqLookupPosBlockAllocator *]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupPosBlockAllocator *]
 *
 */
SeqLookupPosBlockAllocator * Wise2_free_SeqLookupPosBlockAllocator(SeqLookupPosBlockAllocator * obj);
#define free_SeqLookupPosBlockAllocator Wise2_free_SeqLookupPosBlockAllocator


/* Function:  hard_link_SubSeqLookup(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SubSeqLookup *]
 *
 * Return [UNKN ]  Undocumented return value [SubSeqLookup *]
 *
 */
SubSeqLookup * Wise2_hard_link_SubSeqLookup(SubSeqLookup * obj);
#define hard_link_SubSeqLookup Wise2_hard_link_SubSeqLookup


/* Function:  SubSeqLookup_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SubSeqLookup *]
 *
 */
SubSeqLookup * Wise2_SubSeqLookup_alloc(void);
#define SubSeqLookup_alloc Wise2_SubSeqLookup_alloc


/* Function:  free_SubSeqLookup(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SubSeqLookup *]
 *
 * Return [UNKN ]  Undocumented return value [SubSeqLookup *]
 *
 */
SubSeqLookup * Wise2_free_SubSeqLookup(SubSeqLookup * obj);
#define free_SubSeqLookup Wise2_free_SubSeqLookup


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_SeqLookupPosBlockAllocator(SeqLookupPosBlock ** list,int i,int j) ;
#define swap_SeqLookupPosBlockAllocator Wise2_swap_SeqLookupPosBlockAllocator
void Wise2_qsort_SeqLookupPosBlockAllocator(SeqLookupPosBlock ** list,int left,int right,int (*comp)(SeqLookupPosBlock * ,SeqLookupPosBlock * ));
#define qsort_SeqLookupPosBlockAllocator Wise2_qsort_SeqLookupPosBlockAllocator
void Wise2_sort_SeqLookupPosBlockAllocator(SeqLookupPosBlockAllocator * obj,int (*comp)(SeqLookupPosBlock *, SeqLookupPosBlock *));
#define sort_SeqLookupPosBlockAllocator Wise2_sort_SeqLookupPosBlockAllocator
boolean Wise2_expand_SeqLookupPosBlockAllocator(SeqLookupPosBlockAllocator * obj,int len);
#define expand_SeqLookupPosBlockAllocator Wise2_expand_SeqLookupPosBlockAllocator

#ifdef _cplusplus
}
#endif

#endif
