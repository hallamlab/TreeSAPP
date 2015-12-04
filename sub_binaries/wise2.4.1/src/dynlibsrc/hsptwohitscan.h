#ifndef DYNAMITEhsptwohitscanHEADERFILE
#define DYNAMITEhsptwohitscanHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "hsplookupscan.h"


#define TwoHitBufferLISTLENGTH 4096
#define TwoHitSequenceLISTLENGTH 90

#define TWOHIT_FIRST_ENTRY  56
#define TWOHIT_FIRST_STORED 57
#define TWOHIT_HANDLED      58


typedef struct TwoHitStore {
   int hit_state;
   Sequence * target;
   int target_pos;
   int query_pos;
   int diagonal;
 } TwoHitStore;


#define TWOHIT_BLOCK_SIZE   40000
#define TWOHIT_BLOCK_DEPTH  1000

typedef struct TwoHitStoreBlockAllocator {
  TwoHitStore ** block;
  int current_block;
  int current_pos;
  int block_len;
} TwoHitStoreBlockAllocator;



struct Wise2_TwoHitSequence {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    TwoHitStore ** word;     
    int len;/* len for above word  */ 
    int maxlen; /* maxlen for above word */ 
    GHashTable  * diagonal_hash;     
    } ;  
/* TwoHitSequence defined */ 
#ifndef DYNAMITE_DEFINED_TwoHitSequence
typedef struct Wise2_TwoHitSequence Wise2_TwoHitSequence;
#define TwoHitSequence Wise2_TwoHitSequence
#define DYNAMITE_DEFINED_TwoHitSequence
#endif


struct Wise2_TwoHitBuffer {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    TwoHitSequence ** buffer;    
    int len;/* len for above buffer  */ 
    int maxlen; /* maxlen for above buffer */ 
    GHashTable      * target_hash;   
    TwoHitStoreBlockAllocator * thba;    
    } ;  
/* TwoHitBuffer defined */ 
#ifndef DYNAMITE_DEFINED_TwoHitBuffer
typedef struct Wise2_TwoHitBuffer Wise2_TwoHitBuffer;
#define TwoHitBuffer Wise2_TwoHitBuffer
#define DYNAMITE_DEFINED_TwoHitBuffer
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  new_twohit_one_off_HSPScanInterface(sli,mat,drop_off,score_cutoff)
 *
 * Descrip:    Builds a twohit scan interface. This
 *             does expands the query using a matrix but
 *             just be considering off by one cases
 *
 *
 * Arg:                 sli [UNKN ] Undocumented argument [SeqLookupInterface *]
 * Arg:                 mat [UNKN ] Undocumented argument [CompMat *]
 * Arg:            drop_off [UNKN ] Undocumented argument [int]
 * Arg:        score_cutoff [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [HSPScanInterface *]
 *
 */
HSPScanInterface * Wise2_new_twohit_one_off_HSPScanInterface(SeqLookupInterface * sli,CompMat * mat,int drop_off,int score_cutoff);
#define new_twohit_one_off_HSPScanInterface Wise2_new_twohit_one_off_HSPScanInterface


/* Function:  twohit_one_off_HSPscan_scan_query_direct(data,seq,para)
 *
 * Descrip:    Word expansion with two hit semantics
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 * Arg:         seq [UNKN ] Undocumented argument [Sequence *]
 * Arg:        para [UNKN ] Undocumented argument [HSPScanInterfacePara *]
 *
 * Return [UNKN ]  Undocumented return value [LinearHSPmanager *]
 *
 */
LinearHSPmanager * Wise2_twohit_one_off_HSPscan_scan_query_direct(void * data,Sequence * seq,HSPScanInterfacePara * para);
#define twohit_one_off_HSPscan_scan_query_direct Wise2_twohit_one_off_HSPscan_scan_query_direct


/* Function:  add_to_TwoHitBuffer(thb,target,query_pos,target_pos)
 *
 * Descrip:    adds a new potential hit to the TwoHitBuffer, 
 *             returning the HitStore datastructure updated
 *             if hit is 1 then this is first entry
 *             if hit is more than 1 second
 *
 *
 * Arg:               thb [UNKN ] Undocumented argument [TwoHitBuffer *]
 * Arg:            target [UNKN ] Undocumented argument [Sequence *]
 * Arg:         query_pos [UNKN ] Undocumented argument [int]
 * Arg:        target_pos [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [TwoHitStore *]
 *
 */
TwoHitStore * Wise2_add_to_TwoHitBuffer(TwoHitBuffer * thb,Sequence * target,int query_pos,int target_pos);
#define add_to_TwoHitBuffer Wise2_add_to_TwoHitBuffer


/* Function:  new_TwoHitBuffer(void)
 *
 * Descrip:    makes a new TwoHitBuffer ready for use
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TwoHitBuffer *]
 *
 */
TwoHitBuffer * Wise2_new_TwoHitBuffer(void);
#define new_TwoHitBuffer Wise2_new_TwoHitBuffer


/* Function:  new_TwoHitSequence(void)
 *
 * Descrip:    make a new TwoHitSequence ready for use
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TwoHitSequence *]
 *
 */
TwoHitSequence * Wise2_new_TwoHitSequence(void);
#define new_TwoHitSequence Wise2_new_TwoHitSequence


/* Function:  free_TwoHitSequence(t)
 *
 * Descrip:    Frees the TwoHitSequence
 *
 *
 * Arg:        t [UNKN ] Undocumented argument [TwoHitSequence *]
 *
 * Return [UNKN ]  Undocumented return value [TwoHitSequence *]
 *
 */
TwoHitSequence * Wise2_free_TwoHitSequence(TwoHitSequence * t);
#define free_TwoHitSequence Wise2_free_TwoHitSequence


/* Function:  free_TwoHitBuffer(t)
 *
 * Descrip:    Frees the TwoHitBuffer
 *
 *
 * Arg:        t [UNKN ] Undocumented argument [TwoHitBuffer *]
 *
 * Return [UNKN ]  Undocumented return value [TwoHitBuffer *]
 *
 */
TwoHitBuffer * Wise2_free_TwoHitBuffer(TwoHitBuffer * t);
#define free_TwoHitBuffer Wise2_free_TwoHitBuffer


/* Function:  TwoHitStore_alloc(void)
 *
 * Descrip:    allocator for twohitstore
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TwoHitStore *]
 *
 */
TwoHitStore * Wise2_TwoHitStore_alloc(void);
#define TwoHitStore_alloc Wise2_TwoHitStore_alloc


/* Function:  free_TwoHitStore(t)
 *
 * Descrip:    free for twohitstore
 *
 *
 * Arg:        t [UNKN ] Undocumented argument [TwoHitStore *]
 *
 * Return [UNKN ]  Undocumented return value [TwoHitStore *]
 *
 */
TwoHitStore * Wise2_free_TwoHitStore(TwoHitStore * t);
#define free_TwoHitStore Wise2_free_TwoHitStore


/* Function:  new_TwoHitStore_from_Allocator(a)
 *
 * Descrip:    gets a new store from a block allocator
 *
 *
 * Arg:        a [UNKN ] Undocumented argument [TwoHitStoreBlockAllocator *]
 *
 * Return [UNKN ]  Undocumented return value [TwoHitStore *]
 *
 */
TwoHitStore * Wise2_new_TwoHitStore_from_Allocator(TwoHitStoreBlockAllocator * a);
#define new_TwoHitStore_from_Allocator Wise2_new_TwoHitStore_from_Allocator


/* Function:  free_TwoHitStoreBlockAllocator(a)
 *
 * Descrip:    frees a block allocator
 *
 *
 * Arg:        a [UNKN ] Undocumented argument [TwoHitStoreBlockAllocator *]
 *
 * Return [UNKN ]  Undocumented return value [TwoHitStoreBlockAllocator *]
 *
 */
TwoHitStoreBlockAllocator * Wise2_free_TwoHitStoreBlockAllocator(TwoHitStoreBlockAllocator * a);
#define free_TwoHitStoreBlockAllocator Wise2_free_TwoHitStoreBlockAllocator


/* Function:  new_TwoHitStoreBlockAllocator(void)
 *
 * Descrip:    makes a new block allocator
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TwoHitStoreBlockAllocator *]
 *
 */
TwoHitStoreBlockAllocator * Wise2_new_TwoHitStoreBlockAllocator(void);
#define new_TwoHitStoreBlockAllocator Wise2_new_TwoHitStoreBlockAllocator


/* Function:  add_TwoHitSequence(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [TwoHitSequence *]
 * Arg:        add [OWNER] Object to add to the list [TwoHitStore *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_TwoHitSequence(TwoHitSequence * obj,TwoHitStore * add);
#define add_TwoHitSequence Wise2_add_TwoHitSequence


/* Function:  flush_TwoHitSequence(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [TwoHitSequence *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_TwoHitSequence(TwoHitSequence * obj);
#define flush_TwoHitSequence Wise2_flush_TwoHitSequence


/* Function:  TwoHitSequence_alloc_std(void)
 *
 * Descrip:    Equivalent to TwoHitSequence_alloc_len(TwoHitSequenceLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TwoHitSequence *]
 *
 */
TwoHitSequence * Wise2_TwoHitSequence_alloc_std(void);
#define TwoHitSequence_alloc_std Wise2_TwoHitSequence_alloc_std


/* Function:  TwoHitSequence_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [TwoHitSequence *]
 *
 */
TwoHitSequence * Wise2_TwoHitSequence_alloc_len(int len);
#define TwoHitSequence_alloc_len Wise2_TwoHitSequence_alloc_len


/* Function:  hard_link_TwoHitSequence(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [TwoHitSequence *]
 *
 * Return [UNKN ]  Undocumented return value [TwoHitSequence *]
 *
 */
TwoHitSequence * Wise2_hard_link_TwoHitSequence(TwoHitSequence * obj);
#define hard_link_TwoHitSequence Wise2_hard_link_TwoHitSequence


/* Function:  TwoHitSequence_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TwoHitSequence *]
 *
 */
TwoHitSequence * Wise2_TwoHitSequence_alloc(void);
#define TwoHitSequence_alloc Wise2_TwoHitSequence_alloc


/* Function:  add_TwoHitBuffer(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [TwoHitBuffer *]
 * Arg:        add [OWNER] Object to add to the list [TwoHitSequence *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_TwoHitBuffer(TwoHitBuffer * obj,TwoHitSequence * add);
#define add_TwoHitBuffer Wise2_add_TwoHitBuffer


/* Function:  flush_TwoHitBuffer(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [TwoHitBuffer *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_TwoHitBuffer(TwoHitBuffer * obj);
#define flush_TwoHitBuffer Wise2_flush_TwoHitBuffer


/* Function:  TwoHitBuffer_alloc_std(void)
 *
 * Descrip:    Equivalent to TwoHitBuffer_alloc_len(TwoHitBufferLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TwoHitBuffer *]
 *
 */
TwoHitBuffer * Wise2_TwoHitBuffer_alloc_std(void);
#define TwoHitBuffer_alloc_std Wise2_TwoHitBuffer_alloc_std


/* Function:  TwoHitBuffer_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [TwoHitBuffer *]
 *
 */
TwoHitBuffer * Wise2_TwoHitBuffer_alloc_len(int len);
#define TwoHitBuffer_alloc_len Wise2_TwoHitBuffer_alloc_len


/* Function:  hard_link_TwoHitBuffer(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [TwoHitBuffer *]
 *
 * Return [UNKN ]  Undocumented return value [TwoHitBuffer *]
 *
 */
TwoHitBuffer * Wise2_hard_link_TwoHitBuffer(TwoHitBuffer * obj);
#define hard_link_TwoHitBuffer Wise2_hard_link_TwoHitBuffer


/* Function:  TwoHitBuffer_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TwoHitBuffer *]
 *
 */
TwoHitBuffer * Wise2_TwoHitBuffer_alloc(void);
#define TwoHitBuffer_alloc Wise2_TwoHitBuffer_alloc


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_TwoHitSequence(TwoHitStore ** list,int i,int j) ;
#define swap_TwoHitSequence Wise2_swap_TwoHitSequence
void Wise2_qsort_TwoHitSequence(TwoHitStore ** list,int left,int right,int (*comp)(TwoHitStore * ,TwoHitStore * ));
#define qsort_TwoHitSequence Wise2_qsort_TwoHitSequence
void Wise2_sort_TwoHitSequence(TwoHitSequence * obj,int (*comp)(TwoHitStore *, TwoHitStore *));
#define sort_TwoHitSequence Wise2_sort_TwoHitSequence
boolean Wise2_expand_TwoHitSequence(TwoHitSequence * obj,int len);
#define expand_TwoHitSequence Wise2_expand_TwoHitSequence
void Wise2_swap_TwoHitBuffer(TwoHitSequence ** list,int i,int j) ;
#define swap_TwoHitBuffer Wise2_swap_TwoHitBuffer
void Wise2_qsort_TwoHitBuffer(TwoHitSequence ** list,int left,int right,int (*comp)(TwoHitSequence * ,TwoHitSequence * ));
#define qsort_TwoHitBuffer Wise2_qsort_TwoHitBuffer
void Wise2_sort_TwoHitBuffer(TwoHitBuffer * obj,int (*comp)(TwoHitSequence *, TwoHitSequence *));
#define sort_TwoHitBuffer Wise2_sort_TwoHitBuffer
boolean Wise2_expand_TwoHitBuffer(TwoHitBuffer * obj,int len);
#define expand_TwoHitBuffer Wise2_expand_TwoHitBuffer

#ifdef _cplusplus
}
#endif

#endif
