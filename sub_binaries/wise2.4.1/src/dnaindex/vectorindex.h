#ifndef DYNAMITEvectorindexHEADERFILE
#define DYNAMITEvectorindexHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"

#define KMER_VECTOR_LENGTH 4069
#define KMER_VECTOR_SIZE   6

typedef struct KmerVectorPosition {
  char kmer_vector[KMER_VECTOR_LENGTH];
  long int position;
} KmerVectorPosition;

#define KmerVectorPositionSetLISTLENGTH 4096  

struct Wise2_KmerVectorPositionSet {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    KmerVectorPosition ** vec;   
    int len;/* len for above vec  */ 
    int maxlen; /* maxlen for above vec */ 
    } ;  
/* KmerVectorPositionSet defined */ 
#ifndef DYNAMITE_DEFINED_KmerVectorPositionSet
typedef struct Wise2_KmerVectorPositionSet Wise2_KmerVectorPositionSet;
#define KmerVectorPositionSet Wise2_KmerVectorPositionSet
#define DYNAMITE_DEFINED_KmerVectorPositionSet
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  add_KmerVectorPositionSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [KmerVectorPositionSet *]
 * Arg:        add [OWNER] Object to add to the list [KmerVectorPosition *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_KmerVectorPositionSet(KmerVectorPositionSet * obj,KmerVectorPosition * add);
#define add_KmerVectorPositionSet Wise2_add_KmerVectorPositionSet


/* Function:  flush_KmerVectorPositionSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [KmerVectorPositionSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_KmerVectorPositionSet(KmerVectorPositionSet * obj);
#define flush_KmerVectorPositionSet Wise2_flush_KmerVectorPositionSet


/* Function:  KmerVectorPositionSet_alloc_std(void)
 *
 * Descrip:    Equivalent to KmerVectorPositionSet_alloc_len(KmerVectorPositionSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [KmerVectorPositionSet *]
 *
 */
KmerVectorPositionSet * Wise2_KmerVectorPositionSet_alloc_std(void);
#define KmerVectorPositionSet_alloc_std Wise2_KmerVectorPositionSet_alloc_std


/* Function:  KmerVectorPositionSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [KmerVectorPositionSet *]
 *
 */
KmerVectorPositionSet * Wise2_KmerVectorPositionSet_alloc_len(int len);
#define KmerVectorPositionSet_alloc_len Wise2_KmerVectorPositionSet_alloc_len


/* Function:  hard_link_KmerVectorPositionSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [KmerVectorPositionSet *]
 *
 * Return [UNKN ]  Undocumented return value [KmerVectorPositionSet *]
 *
 */
KmerVectorPositionSet * Wise2_hard_link_KmerVectorPositionSet(KmerVectorPositionSet * obj);
#define hard_link_KmerVectorPositionSet Wise2_hard_link_KmerVectorPositionSet


/* Function:  KmerVectorPositionSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [KmerVectorPositionSet *]
 *
 */
KmerVectorPositionSet * Wise2_KmerVectorPositionSet_alloc(void);
#define KmerVectorPositionSet_alloc Wise2_KmerVectorPositionSet_alloc


/* Function:  free_KmerVectorPositionSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [KmerVectorPositionSet *]
 *
 * Return [UNKN ]  Undocumented return value [KmerVectorPositionSet *]
 *
 */
KmerVectorPositionSet * Wise2_free_KmerVectorPositionSet(KmerVectorPositionSet * obj);
#define free_KmerVectorPositionSet Wise2_free_KmerVectorPositionSet


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
double Wise2_vector_KmerVectorPosition(KmerVectorPosition * a,KmerVectorPosition * b);
#define vector_KmerVectorPosition Wise2_vector_KmerVectorPosition
KmerVectorPosition * Wise2_new_KmerVectorPosition(void);
#define new_KmerVectorPosition Wise2_new_KmerVectorPosition
KmerVectorPositionSet * Wise2_build_KmerVectorPositionSet(Sequence * input,long int start_pos,int step_size);
#define build_KmerVectorPositionSet Wise2_build_KmerVectorPositionSet
KmerVectorPosition * Wise2_free_KmerVectorPosition(KmerVectorPosition * vec);
#define free_KmerVectorPosition Wise2_free_KmerVectorPosition


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_KmerVectorPositionSet(KmerVectorPosition ** list,int i,int j) ;
#define swap_KmerVectorPositionSet Wise2_swap_KmerVectorPositionSet
void Wise2_qsort_KmerVectorPositionSet(KmerVectorPosition ** list,int left,int right,int (*comp)(KmerVectorPosition * ,KmerVectorPosition * ));
#define qsort_KmerVectorPositionSet Wise2_qsort_KmerVectorPositionSet
void Wise2_sort_KmerVectorPositionSet(KmerVectorPositionSet * obj,int (*comp)(KmerVectorPosition *, KmerVectorPosition *));
#define sort_KmerVectorPositionSet Wise2_sort_KmerVectorPositionSet
boolean Wise2_expand_KmerVectorPositionSet(KmerVectorPositionSet * obj,int len);
#define expand_KmerVectorPositionSet Wise2_expand_KmerVectorPositionSet

#ifdef _cplusplus
}
#endif

#endif
