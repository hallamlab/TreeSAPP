#ifndef DYNAMITEkmer_countHEADERFILE
#define DYNAMITEkmer_countHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"

#define KMER_COUNT_BLOCKSIZE 400
#define KmerCountAllocatorLISTLENGTH 128

struct Wise2_KmerCount {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int count;   
    } ;  
/* KmerCount defined */ 
#ifndef DYNAMITE_DEFINED_KmerCount
typedef struct Wise2_KmerCount Wise2_KmerCount;
#define KmerCount Wise2_KmerCount
#define DYNAMITE_DEFINED_KmerCount
#endif


struct Wise2_KmerCountBlock {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    KmerCount count[KMER_COUNT_BLOCKSIZE];   
    } ;  
/* KmerCountBlock defined */ 
#ifndef DYNAMITE_DEFINED_KmerCountBlock
typedef struct Wise2_KmerCountBlock Wise2_KmerCountBlock;
#define KmerCountBlock Wise2_KmerCountBlock
#define DYNAMITE_DEFINED_KmerCountBlock
#endif


struct Wise2_KmerCountAllocator {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    KmerCountBlock ** block;     
    int len;/* len for above block  */ 
    int maxlen; /* maxlen for above block */ 
    int current_count;   
    } ;  
/* KmerCountAllocator defined */ 
#ifndef DYNAMITE_DEFINED_KmerCountAllocator
typedef struct Wise2_KmerCountAllocator Wise2_KmerCountAllocator;
#define KmerCountAllocator Wise2_KmerCountAllocator
#define DYNAMITE_DEFINED_KmerCountAllocator
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_KmerCount(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [KmerCount *]
 *
 * Return [UNKN ]  Undocumented return value [KmerCount *]
 *
 */
KmerCount * Wise2_hard_link_KmerCount(KmerCount * obj);
#define hard_link_KmerCount Wise2_hard_link_KmerCount


/* Function:  KmerCount_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [KmerCount *]
 *
 */
KmerCount * Wise2_KmerCount_alloc(void);
#define KmerCount_alloc Wise2_KmerCount_alloc


/* Function:  free_KmerCount(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [KmerCount *]
 *
 * Return [UNKN ]  Undocumented return value [KmerCount *]
 *
 */
KmerCount * Wise2_free_KmerCount(KmerCount * obj);
#define free_KmerCount Wise2_free_KmerCount


/* Function:  hard_link_KmerCountBlock(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [KmerCountBlock *]
 *
 * Return [UNKN ]  Undocumented return value [KmerCountBlock *]
 *
 */
KmerCountBlock * Wise2_hard_link_KmerCountBlock(KmerCountBlock * obj);
#define hard_link_KmerCountBlock Wise2_hard_link_KmerCountBlock


/* Function:  KmerCountBlock_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [KmerCountBlock *]
 *
 */
KmerCountBlock * Wise2_KmerCountBlock_alloc(void);
#define KmerCountBlock_alloc Wise2_KmerCountBlock_alloc


/* Function:  free_KmerCountBlock(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [KmerCountBlock *]
 *
 * Return [UNKN ]  Undocumented return value [KmerCountBlock *]
 *
 */
KmerCountBlock * Wise2_free_KmerCountBlock(KmerCountBlock * obj);
#define free_KmerCountBlock Wise2_free_KmerCountBlock


/* Function:  add_KmerCountAllocator(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [KmerCountAllocator *]
 * Arg:        add [OWNER] Object to add to the list [KmerCountBlock *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_KmerCountAllocator(KmerCountAllocator * obj,KmerCountBlock * add);
#define add_KmerCountAllocator Wise2_add_KmerCountAllocator


/* Function:  flush_KmerCountAllocator(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [KmerCountAllocator *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_KmerCountAllocator(KmerCountAllocator * obj);
#define flush_KmerCountAllocator Wise2_flush_KmerCountAllocator


/* Function:  KmerCountAllocator_alloc_std(void)
 *
 * Descrip:    Equivalent to KmerCountAllocator_alloc_len(KmerCountAllocatorLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [KmerCountAllocator *]
 *
 */
KmerCountAllocator * Wise2_KmerCountAllocator_alloc_std(void);
#define KmerCountAllocator_alloc_std Wise2_KmerCountAllocator_alloc_std


/* Function:  KmerCountAllocator_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [KmerCountAllocator *]
 *
 */
KmerCountAllocator * Wise2_KmerCountAllocator_alloc_len(int len);
#define KmerCountAllocator_alloc_len Wise2_KmerCountAllocator_alloc_len


/* Function:  hard_link_KmerCountAllocator(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [KmerCountAllocator *]
 *
 * Return [UNKN ]  Undocumented return value [KmerCountAllocator *]
 *
 */
KmerCountAllocator * Wise2_hard_link_KmerCountAllocator(KmerCountAllocator * obj);
#define hard_link_KmerCountAllocator Wise2_hard_link_KmerCountAllocator


/* Function:  KmerCountAllocator_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [KmerCountAllocator *]
 *
 */
KmerCountAllocator * Wise2_KmerCountAllocator_alloc(void);
#define KmerCountAllocator_alloc Wise2_KmerCountAllocator_alloc


/* Function:  free_KmerCountAllocator(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [KmerCountAllocator *]
 *
 * Return [UNKN ]  Undocumented return value [KmerCountAllocator *]
 *
 */
KmerCountAllocator * Wise2_free_KmerCountAllocator(KmerCountAllocator * obj);
#define free_KmerCountAllocator Wise2_free_KmerCountAllocator


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
KmerCount * Wise2_new_KmerCount_KmerCountAllocator(KmerCountAllocator * kca);
#define new_KmerCount_KmerCountAllocator Wise2_new_KmerCount_KmerCountAllocator
KmerCountAllocator * Wise2_new_KmerCountAllocator(void);
#define new_KmerCountAllocator Wise2_new_KmerCountAllocator


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_KmerCountAllocator(KmerCountBlock ** list,int i,int j) ;
#define swap_KmerCountAllocator Wise2_swap_KmerCountAllocator
void Wise2_qsort_KmerCountAllocator(KmerCountBlock ** list,int left,int right,int (*comp)(KmerCountBlock * ,KmerCountBlock * ));
#define qsort_KmerCountAllocator Wise2_qsort_KmerCountAllocator
void Wise2_sort_KmerCountAllocator(KmerCountAllocator * obj,int (*comp)(KmerCountBlock *, KmerCountBlock *));
#define sort_KmerCountAllocator Wise2_sort_KmerCountAllocator
boolean Wise2_expand_KmerCountAllocator(KmerCountAllocator * obj,int len);
#define expand_KmerCountAllocator Wise2_expand_KmerCountAllocator

#ifdef _cplusplus
}
#endif

#endif
