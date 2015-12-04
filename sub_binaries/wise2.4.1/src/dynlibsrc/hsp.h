#ifndef DYNAMITEhspHEADERFILE
#define DYNAMITEhspHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "sequence.h"
#include "compmat.h"


#define HSPsetLISTLENGTH 20

#define HSP_BLOCK_SIZE 1024
#define HSPCacheLISTLENGTH 1024

#define LinearHSPmanagerLISTLENGTH 128

#define ON_HSP_MACRO(test,query_pos,target_pos) ( (((test->query_start-test->target_start)!=(query_pos-target_pos))||(query_pos<test->query_start)||(target_pos<test->target_start)||(query_pos-test->query_start>test->length)) ? FALSE : TRUE)
 

struct Wise2_HSP {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Sequence * query;    
    Sequence * target;   
    int query_start;     
    int target_start;    
    int length;  
    int score;   
    char is_in_block;    
    char target_reverse;     
    } ;  
/* HSP defined */ 
#ifndef DYNAMITE_DEFINED_HSP
typedef struct Wise2_HSP Wise2_HSP;
#define HSP Wise2_HSP
#define DYNAMITE_DEFINED_HSP
#endif


struct Wise2_HSPCache {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    HSP ** cache;    
    int len;/* len for above cache  */ 
    int maxlen; /* maxlen for above cache */ 
    int max_cache;   
    } ;  
/* HSPCache defined */ 
#ifndef DYNAMITE_DEFINED_HSPCache
typedef struct Wise2_HSPCache Wise2_HSPCache;
#define HSPCache Wise2_HSPCache
#define DYNAMITE_DEFINED_HSPCache
#endif


struct Wise2_HSPset {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    HSP ** hsp;  
    int len;/* len for above hsp  */ 
    int maxlen; /* maxlen for above hsp */ 
    int score;   
    int best_score;  
    int last_accessed;   
    } ;  
/* HSPset defined */ 
#ifndef DYNAMITE_DEFINED_HSPset
typedef struct Wise2_HSPset Wise2_HSPset;
#define HSPset Wise2_HSPset
#define DYNAMITE_DEFINED_HSPset
#endif


struct Wise2_LinearHSPmanager {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    HSPset ** set;   
    int len;/* len for above set  */ 
    int maxlen; /* maxlen for above set */ 
    Sequence * query;    
    CompMat * mat;   
    int min_score;   
    int width;   
    int tail;    
    int worst_hsp_score;     
    } ;  
/* LinearHSPmanager defined */ 
#ifndef DYNAMITE_DEFINED_LinearHSPmanager
typedef struct Wise2_LinearHSPmanager Wise2_LinearHSPmanager;
#define LinearHSPmanager Wise2_LinearHSPmanager
#define DYNAMITE_DEFINED_LinearHSPmanager
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  on_HSP(test,query_pos,target_pos)
 *
 * Descrip:    tests whether this point is on this test HSP
 *
 *
 * Arg:              test [UNKN ] Undocumented argument [HSP *]
 * Arg:         query_pos [UNKN ] Undocumented argument [int]
 * Arg:        target_pos [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_on_HSP(HSP * test,int query_pos,int target_pos);
#define on_HSP Wise2_on_HSP


/* Function:  compare_HSPset_score_qsort(a,b)
 *
 * Descrip:    sorting linear HSPsets via qsort function
 *
 *
 * Arg:        a [UNKN ] Undocumented argument [const void *]
 * Arg:        b [UNKN ] Undocumented argument [const void *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_compare_HSPset_score_qsort(const void * a,const void * b);
#define compare_HSPset_score_qsort Wise2_compare_HSPset_score_qsort


/* Function:  compare_HSPset_score(one,two)
 *
 * Descrip:    internal function for sort linear HSPsets
 *
 *
 * Arg:        one [UNKN ] Undocumented argument [HSPset *]
 * Arg:        two [UNKN ] Undocumented argument [HSPset *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_compare_HSPset_score(HSPset * one,HSPset * two);
#define compare_HSPset_score Wise2_compare_HSPset_score


/* Function:  sort_HSPset_by_score(*set)
 *
 * Descrip:    Sorts by score
 *
 *
 * Arg:        *set [UNKN ] Undocumented argument [HSPset]
 *
 */
void Wise2_sort_HSPset_by_score(HSPset *set);
#define sort_HSPset_by_score Wise2_sort_HSPset_by_score


/* Function:  new_dna_identical_HSP(query,target,query_pos,target_pos,target_reverse)
 *
 * Descrip:    builds a new HSP for these sequences breaking at first mismatch
 *
 *
 * Arg:                 query [UNKN ] Undocumented argument [Sequence *]
 * Arg:                target [UNKN ] Undocumented argument [Sequence *]
 * Arg:             query_pos [UNKN ] Undocumented argument [int]
 * Arg:            target_pos [UNKN ] Undocumented argument [int]
 * Arg:        target_reverse [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [HSP *]
 *
 */
HSP * Wise2_new_dna_identical_HSP(Sequence * query,Sequence * target,int query_pos,int target_pos,int target_reverse);
#define new_dna_identical_HSP Wise2_new_dna_identical_HSP


/* Function:  new_HSP(cache,query,target,query_pos,target_pos,mat,drop_off)
 *
 * Descrip:    builds a new HSP for these sequences
 *
 *
 * Arg:             cache [UNKN ] Undocumented argument [HSPCache *]
 * Arg:             query [UNKN ] Undocumented argument [Sequence *]
 * Arg:            target [UNKN ] Undocumented argument [Sequence *]
 * Arg:         query_pos [UNKN ] Undocumented argument [int]
 * Arg:        target_pos [UNKN ] Undocumented argument [int]
 * Arg:               mat [UNKN ] Undocumented argument [CompMat *]
 * Arg:          drop_off [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [HSP *]
 *
 */
HSP * Wise2_new_HSP(HSPCache * cache,Sequence * query,Sequence * target,int query_pos,int target_pos,CompMat * mat,int drop_off);
#define new_HSP Wise2_new_HSP


/* Function:  HSP_alloc_cache(hspc)
 *
 * Descrip:    Returns new HSP, using cache if needed
 *
 *
 * Arg:        hspc [UNKN ] Undocumented argument [HSPCache *]
 *
 * Return [UNKN ]  Undocumented return value [HSP *]
 *
 */
HSP * Wise2_HSP_alloc_cache(HSPCache * hspc);
#define HSP_alloc_cache Wise2_HSP_alloc_cache


/* Function:  free_HSP_cache(cache,hsp)
 *
 * Descrip:    Places HSP back into cache, freeing if necessary
 *
 *
 * Arg:        cache [UNKN ] Undocumented argument [HSPCache *]
 * Arg:          hsp [UNKN ] Undocumented argument [HSP *]
 *
 * Return [UNKN ]  Undocumented return value [HSP *]
 *
 */
HSP * Wise2_free_HSP_cache(HSPCache * cache,HSP * hsp);
#define free_HSP_cache Wise2_free_HSP_cache


/* Function:  new_HSPCache(maxsize)
 *
 * Descrip:    Makes a new cache
 *
 *
 * Arg:        maxsize [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [HSPCache *]
 *
 */
HSPCache * Wise2_new_HSPCache(int maxsize);
#define new_HSPCache Wise2_new_HSPCache


/* Function:  show_HSPset(s,ofp)
 *
 * Descrip:    Shows a HSP set
 *
 *
 * Arg:          s [UNKN ] Undocumented argument [HSPset *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_show_HSPset(HSPset * s,FILE * ofp);
#define show_HSPset Wise2_show_HSPset


/* Function:  show_HSP(hsp,linelength,out)
 *
 * Descrip:    Shows an HSP
 *
 *
 * Arg:               hsp [UNKN ] Undocumented argument [HSP *]
 * Arg:        linelength [UNKN ] Undocumented argument [int]
 * Arg:               out [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_show_HSP(HSP * hsp,int linelength,FILE * out);
#define show_HSP Wise2_show_HSP


/* Function:  hard_link_HSP(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [HSP *]
 *
 * Return [UNKN ]  Undocumented return value [HSP *]
 *
 */
HSP * Wise2_hard_link_HSP(HSP * obj);
#define hard_link_HSP Wise2_hard_link_HSP


/* Function:  HSP_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [HSP *]
 *
 */
HSP * Wise2_HSP_alloc(void);
#define HSP_alloc Wise2_HSP_alloc


/* Function:  free_HSP(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [HSP *]
 *
 * Return [UNKN ]  Undocumented return value [HSP *]
 *
 */
HSP * Wise2_free_HSP(HSP * obj);
#define free_HSP Wise2_free_HSP


/* Function:  add_HSPCache(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [HSPCache *]
 * Arg:        add [OWNER] Object to add to the list [HSP *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_HSPCache(HSPCache * obj,HSP * add);
#define add_HSPCache Wise2_add_HSPCache


/* Function:  flush_HSPCache(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [HSPCache *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_HSPCache(HSPCache * obj);
#define flush_HSPCache Wise2_flush_HSPCache


/* Function:  HSPCache_alloc_std(void)
 *
 * Descrip:    Equivalent to HSPCache_alloc_len(HSPCacheLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [HSPCache *]
 *
 */
HSPCache * Wise2_HSPCache_alloc_std(void);
#define HSPCache_alloc_std Wise2_HSPCache_alloc_std


/* Function:  HSPCache_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [HSPCache *]
 *
 */
HSPCache * Wise2_HSPCache_alloc_len(int len);
#define HSPCache_alloc_len Wise2_HSPCache_alloc_len


/* Function:  hard_link_HSPCache(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [HSPCache *]
 *
 * Return [UNKN ]  Undocumented return value [HSPCache *]
 *
 */
HSPCache * Wise2_hard_link_HSPCache(HSPCache * obj);
#define hard_link_HSPCache Wise2_hard_link_HSPCache


/* Function:  HSPCache_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [HSPCache *]
 *
 */
HSPCache * Wise2_HSPCache_alloc(void);
#define HSPCache_alloc Wise2_HSPCache_alloc


/* Function:  free_HSPCache(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [HSPCache *]
 *
 * Return [UNKN ]  Undocumented return value [HSPCache *]
 *
 */
HSPCache * Wise2_free_HSPCache(HSPCache * obj);
#define free_HSPCache Wise2_free_HSPCache


/* Function:  add_HSPset(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [HSPset *]
 * Arg:        add [OWNER] Object to add to the list [HSP *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_HSPset(HSPset * obj,HSP * add);
#define add_HSPset Wise2_add_HSPset


/* Function:  flush_HSPset(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [HSPset *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_HSPset(HSPset * obj);
#define flush_HSPset Wise2_flush_HSPset


/* Function:  HSPset_alloc_std(void)
 *
 * Descrip:    Equivalent to HSPset_alloc_len(HSPsetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [HSPset *]
 *
 */
HSPset * Wise2_HSPset_alloc_std(void);
#define HSPset_alloc_std Wise2_HSPset_alloc_std


/* Function:  HSPset_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [HSPset *]
 *
 */
HSPset * Wise2_HSPset_alloc_len(int len);
#define HSPset_alloc_len Wise2_HSPset_alloc_len


/* Function:  hard_link_HSPset(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [HSPset *]
 *
 * Return [UNKN ]  Undocumented return value [HSPset *]
 *
 */
HSPset * Wise2_hard_link_HSPset(HSPset * obj);
#define hard_link_HSPset Wise2_hard_link_HSPset


/* Function:  HSPset_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [HSPset *]
 *
 */
HSPset * Wise2_HSPset_alloc(void);
#define HSPset_alloc Wise2_HSPset_alloc


/* Function:  free_HSPset(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [HSPset *]
 *
 * Return [UNKN ]  Undocumented return value [HSPset *]
 *
 */
HSPset * Wise2_free_HSPset(HSPset * obj);
#define free_HSPset Wise2_free_HSPset


/* Function:  add_LinearHSPmanager(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [LinearHSPmanager *]
 * Arg:        add [OWNER] Object to add to the list [HSPset *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_LinearHSPmanager(LinearHSPmanager * obj,HSPset * add);
#define add_LinearHSPmanager Wise2_add_LinearHSPmanager


/* Function:  flush_LinearHSPmanager(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [LinearHSPmanager *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_LinearHSPmanager(LinearHSPmanager * obj);
#define flush_LinearHSPmanager Wise2_flush_LinearHSPmanager


/* Function:  LinearHSPmanager_alloc_std(void)
 *
 * Descrip:    Equivalent to LinearHSPmanager_alloc_len(LinearHSPmanagerLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LinearHSPmanager *]
 *
 */
LinearHSPmanager * Wise2_LinearHSPmanager_alloc_std(void);
#define LinearHSPmanager_alloc_std Wise2_LinearHSPmanager_alloc_std


/* Function:  LinearHSPmanager_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [LinearHSPmanager *]
 *
 */
LinearHSPmanager * Wise2_LinearHSPmanager_alloc_len(int len);
#define LinearHSPmanager_alloc_len Wise2_LinearHSPmanager_alloc_len


/* Function:  hard_link_LinearHSPmanager(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [LinearHSPmanager *]
 *
 * Return [UNKN ]  Undocumented return value [LinearHSPmanager *]
 *
 */
LinearHSPmanager * Wise2_hard_link_LinearHSPmanager(LinearHSPmanager * obj);
#define hard_link_LinearHSPmanager Wise2_hard_link_LinearHSPmanager


/* Function:  LinearHSPmanager_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LinearHSPmanager *]
 *
 */
LinearHSPmanager * Wise2_LinearHSPmanager_alloc(void);
#define LinearHSPmanager_alloc Wise2_LinearHSPmanager_alloc


/* Function:  free_LinearHSPmanager(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [LinearHSPmanager *]
 *
 * Return [UNKN ]  Undocumented return value [LinearHSPmanager *]
 *
 */
LinearHSPmanager * Wise2_free_LinearHSPmanager(LinearHSPmanager * obj);
#define free_LinearHSPmanager Wise2_free_LinearHSPmanager


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
int Wise2_compare_HSP_score(HSP * one,HSP * two);
#define compare_HSP_score Wise2_compare_HSP_score
void Wise2_swap_HSPCache(HSP ** list,int i,int j) ;
#define swap_HSPCache Wise2_swap_HSPCache
void Wise2_qsort_HSPCache(HSP ** list,int left,int right,int (*comp)(HSP * ,HSP * ));
#define qsort_HSPCache Wise2_qsort_HSPCache
void Wise2_sort_HSPCache(HSPCache * obj,int (*comp)(HSP *, HSP *));
#define sort_HSPCache Wise2_sort_HSPCache
boolean Wise2_expand_HSPCache(HSPCache * obj,int len);
#define expand_HSPCache Wise2_expand_HSPCache
void Wise2_swap_HSPset(HSP ** list,int i,int j) ;
#define swap_HSPset Wise2_swap_HSPset
void Wise2_qsort_HSPset(HSP ** list,int left,int right,int (*comp)(HSP * ,HSP * ));
#define qsort_HSPset Wise2_qsort_HSPset
void Wise2_sort_HSPset(HSPset * obj,int (*comp)(HSP *, HSP *));
#define sort_HSPset Wise2_sort_HSPset
boolean Wise2_expand_HSPset(HSPset * obj,int len);
#define expand_HSPset Wise2_expand_HSPset
void Wise2_swap_LinearHSPmanager(HSPset ** list,int i,int j) ;
#define swap_LinearHSPmanager Wise2_swap_LinearHSPmanager
void Wise2_qsort_LinearHSPmanager(HSPset ** list,int left,int right,int (*comp)(HSPset * ,HSPset * ));
#define qsort_LinearHSPmanager Wise2_qsort_LinearHSPmanager
void Wise2_sort_LinearHSPmanager(LinearHSPmanager * obj,int (*comp)(HSPset *, HSPset *));
#define sort_LinearHSPmanager Wise2_sort_LinearHSPmanager
boolean Wise2_expand_LinearHSPmanager(LinearHSPmanager * obj,int len);
#define expand_LinearHSPmanager Wise2_expand_LinearHSPmanager

#ifdef _cplusplus
}
#endif

#endif
