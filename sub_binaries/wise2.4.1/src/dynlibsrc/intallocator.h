#ifndef DYNAMITEintallocatorHEADERFILE
#define DYNAMITEintallocatorHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "wisebase.h"
#include <glib.h>

typedef union int_allocator_header {
  struct  {
    union int_allocator_header * next; /* when free */
  }s ;
  int dummy; /* to ensure alignment */
} IntAllocatorHeader; 

#define IntAllocator_BLOCKSIZE  512
#define IntAllocator_MEMORY_BLOCK_SIZE 512

/*
 #define IntAllocator_PARANOIA 1
*/

struct Wise2_IntAllocator {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int size;    
    IntAllocatorHeader * start_of_free;  
    void ** allocated_blocks;    
    int max_allocated_blocks;    
    int current_allocated_block;     
    } ;  
/* IntAllocator defined */ 
#ifndef DYNAMITE_DEFINED_IntAllocator
typedef struct Wise2_IntAllocator Wise2_IntAllocator;
#define IntAllocator Wise2_IntAllocator
#define DYNAMITE_DEFINED_IntAllocator
#endif


struct Wise2_IntAllocatorSet {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    IntAllocator ** allocator_set;   
    int max_size;    
    } ;  
/* IntAllocatorSet defined */ 
#ifndef DYNAMITE_DEFINED_IntAllocatorSet
typedef struct Wise2_IntAllocatorSet Wise2_IntAllocatorSet;
#define IntAllocatorSet Wise2_IntAllocatorSet
#define DYNAMITE_DEFINED_IntAllocatorSet
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  new_IntAllocatorSet(max_size)
 *
 * Descrip:    Makes a new IntAllocatorSet up to a certain size
 *
 *
 * Arg:        max_size [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [IntAllocatorSet *]
 *
 */
IntAllocatorSet * Wise2_new_IntAllocatorSet(int max_size);
#define new_IntAllocatorSet Wise2_new_IntAllocatorSet


/* Function:  realloc_intarray_IntAllocatorSet(ias,current,old_size,new_size)
 *
 * Descrip:    reallocates a piece of memory
 *
 *
 * Arg:             ias [UNKN ] Undocumented argument [IntAllocatorSet *]
 * Arg:         current [UNKN ] Undocumented argument [int *]
 * Arg:        old_size [UNKN ] Undocumented argument [int]
 * Arg:        new_size [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int *]
 *
 */
int * Wise2_realloc_intarray_IntAllocatorSet(IntAllocatorSet * ias,int * current,int old_size,int new_size);
#define realloc_intarray_IntAllocatorSet Wise2_realloc_intarray_IntAllocatorSet


/* Function:  free_intarray_IntAllocatorSet(ias,array,size)
 *
 * Descrip:    Frees a piece of memory in a IntAllocatorSet
 *
 *
 * Arg:          ias [UNKN ] Undocumented argument [IntAllocatorSet *]
 * Arg:        array [UNKN ] Undocumented argument [int *]
 * Arg:         size [UNKN ] Undocumented argument [int]
 *
 */
void Wise2_free_intarray_IntAllocatorSet(IntAllocatorSet * ias,int * array,int size);
#define free_intarray_IntAllocatorSet Wise2_free_intarray_IntAllocatorSet


/* Function:  alloc_intarray_IntAllocatorSet(ias,size)
 *
 * Descrip:    Allocates a new piece of memory in a IntAllocatorSet
 *
 *
 * Arg:         ias [UNKN ] Undocumented argument [IntAllocatorSet *]
 * Arg:        size [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int *]
 *
 */
int * Wise2_alloc_intarray_IntAllocatorSet(IntAllocatorSet * ias,int size);
#define alloc_intarray_IntAllocatorSet Wise2_alloc_intarray_IntAllocatorSet


/* Function:  new_IntAllocator(size)
 *
 * Descrip:    Makes a new int allocator
 *
 *
 * Arg:        size [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [IntAllocator *]
 *
 */
IntAllocator * Wise2_new_IntAllocator(int size);
#define new_IntAllocator Wise2_new_IntAllocator


/* Function:  is_acyclic_IntAllocator(ia)
 *
 * Descrip:    Detect cycle
 *
 *
 * Arg:        ia [UNKN ] Undocumented argument [IntAllocator *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_is_acyclic_IntAllocator(IntAllocator * ia);
#define is_acyclic_IntAllocator Wise2_is_acyclic_IntAllocator


/* Function:  show_allocator_status_IntAllocatorSet(ias,ofp)
 *
 * Descrip:    Show status of intallocator set
 *
 *
 * Arg:        ias [UNKN ] Undocumented argument [IntAllocatorSet *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_show_allocator_status_IntAllocatorSet(IntAllocatorSet * ias,FILE * ofp);
#define show_allocator_status_IntAllocatorSet Wise2_show_allocator_status_IntAllocatorSet


/* Function:  show_allocator_status_IntAllocator(ia,ofp)
 *
 * Descrip:    Shows allocator status
 *
 *
 * Arg:         ia [UNKN ] Undocumented argument [IntAllocator *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_show_allocator_status_IntAllocator(IntAllocator * ia,FILE * ofp);
#define show_allocator_status_IntAllocator Wise2_show_allocator_status_IntAllocator


/* Function:  free_intarray_IntAllocator(ia,array)
 *
 * Descrip:    returns an integer back to the pool. NOTE:
 *             This integer * must have come from the pool otherwise
 *             there is going to be a disaster...
 *
 *
 * Arg:           ia [UNKN ] Undocumented argument [IntAllocator *]
 * Arg:        array [UNKN ] Undocumented argument [int *]
 *
 */
void Wise2_free_intarray_IntAllocator(IntAllocator * ia,int * array);
#define free_intarray_IntAllocator Wise2_free_intarray_IntAllocator


/* Function:  alloc_intarray_IntAllocator(ia)
 *
 * Descrip:    returns an integer * from this allocator
 *
 *
 * Arg:        ia [UNKN ] Undocumented argument [IntAllocator *]
 *
 * Return [UNKN ]  Undocumented return value [int *]
 *
 */
int * Wise2_alloc_intarray_IntAllocator(IntAllocator * ia);
#define alloc_intarray_IntAllocator Wise2_alloc_intarray_IntAllocator


/* Function:  allocate_new_block_IntAllocator(ia)
 *
 * Descrip:    internal function to allocate and segment
 *             a block read for use, storing the memory
 *             and segmenting it correctly. Returned pointer
 *             is the first header block
 *
 *
 * Arg:        ia [UNKN ] Undocumented argument [IntAllocator *]
 *
 * Return [UNKN ]  Undocumented return value [IntAllocatorHeader *]
 *
 */
IntAllocatorHeader * Wise2_allocate_new_block_IntAllocator(IntAllocator * ia);
#define allocate_new_block_IntAllocator Wise2_allocate_new_block_IntAllocator


/* Function:  add_new_block_to_memory_handlers_IA(ia,new_block)
 *
 * Descrip:    internal function to ensure new block is added, with growth of block array
 *             if needed
 *
 *
 * Arg:               ia [UNKN ] Undocumented argument [IntAllocator *]
 * Arg:        new_block [UNKN ] Undocumented argument [void *]
 *
 */
void Wise2_add_new_block_to_memory_handlers_IA(IntAllocator * ia,void * new_block);
#define add_new_block_to_memory_handlers_IA Wise2_add_new_block_to_memory_handlers_IA


/* Function:  hard_link_IntAllocator(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [IntAllocator *]
 *
 * Return [UNKN ]  Undocumented return value [IntAllocator *]
 *
 */
IntAllocator * Wise2_hard_link_IntAllocator(IntAllocator * obj);
#define hard_link_IntAllocator Wise2_hard_link_IntAllocator


/* Function:  IntAllocator_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [IntAllocator *]
 *
 */
IntAllocator * Wise2_IntAllocator_alloc(void);
#define IntAllocator_alloc Wise2_IntAllocator_alloc


/* Function:  free_IntAllocator(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [IntAllocator *]
 *
 * Return [UNKN ]  Undocumented return value [IntAllocator *]
 *
 */
IntAllocator * Wise2_free_IntAllocator(IntAllocator * obj);
#define free_IntAllocator Wise2_free_IntAllocator


/* Function:  hard_link_IntAllocatorSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [IntAllocatorSet *]
 *
 * Return [UNKN ]  Undocumented return value [IntAllocatorSet *]
 *
 */
IntAllocatorSet * Wise2_hard_link_IntAllocatorSet(IntAllocatorSet * obj);
#define hard_link_IntAllocatorSet Wise2_hard_link_IntAllocatorSet


/* Function:  IntAllocatorSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [IntAllocatorSet *]
 *
 */
IntAllocatorSet * Wise2_IntAllocatorSet_alloc(void);
#define IntAllocatorSet_alloc Wise2_IntAllocatorSet_alloc


/* Function:  free_IntAllocatorSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [IntAllocatorSet *]
 *
 * Return [UNKN ]  Undocumented return value [IntAllocatorSet *]
 *
 */
IntAllocatorSet * Wise2_free_IntAllocatorSet(IntAllocatorSet * obj);
#define free_IntAllocatorSet Wise2_free_IntAllocatorSet


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
