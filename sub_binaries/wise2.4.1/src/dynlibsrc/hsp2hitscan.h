#ifndef DYNAMITEhsp2hitscanHEADERFILE
#define DYNAMITEhsp2hitscanHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "hsplookupscan.h"


#define SeedHitSetLISTLENGTH 200
#define SeedHitManagerLISTLENGTH 500

struct Wise2_SeedHit {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int diagonal;    
    int posi;    
    int posj;    
    int used;    
    } ;  
/* SeedHit defined */ 
#ifndef DYNAMITE_DEFINED_SeedHit
typedef struct Wise2_SeedHit Wise2_SeedHit;
#define SeedHit Wise2_SeedHit
#define DYNAMITE_DEFINED_SeedHit
#endif


struct Wise2_SeedHitSet {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    SeedHit ** sh;   
    int len;/* len for above sh  */ 
    int maxlen; /* maxlen for above sh */ 
    Sequence * target;   
    int score;   
    } ;  
/* SeedHitSet defined */ 
#ifndef DYNAMITE_DEFINED_SeedHitSet
typedef struct Wise2_SeedHitSet Wise2_SeedHitSet;
#define SeedHitSet Wise2_SeedHitSet
#define DYNAMITE_DEFINED_SeedHitSet
#endif


struct Wise2_SeedHitManager {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    GHashTable * hash;   
    SeedHitSet ** set;   
    int len;/* len for above set  */ 
    int maxlen; /* maxlen for above set */ 
    } ;  
/* SeedHitManager defined */ 
#ifndef DYNAMITE_DEFINED_SeedHitManager
typedef struct Wise2_SeedHitManager Wise2_SeedHitManager;
#define SeedHitManager Wise2_SeedHitManager
#define DYNAMITE_DEFINED_SeedHitManager
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  new_SeedHit(posi,posj)
 *
 * Descrip:    Makes a new SeedHit
 *
 *
 * Arg:        posi [UNKN ] Undocumented argument [int]
 * Arg:        posj [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [SeedHit *]
 *
 */
SeedHit * Wise2_new_SeedHit(int posi,int posj);
#define new_SeedHit Wise2_new_SeedHit


/* Function:  add_SeedHit_SeedHitManager(shm,t,posi,posj,score)
 *
 * Descrip:    Adds a new position to the manger, allocating any
 *             new datastructures
 *
 *
 * Arg:          shm [UNKN ] Undocumented argument [SeedHitManager *]
 * Arg:            t [UNKN ] Undocumented argument [Sequence *]
 * Arg:         posi [UNKN ] Undocumented argument [int]
 * Arg:         posj [UNKN ] Undocumented argument [int]
 * Arg:        score [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_SeedHit_SeedHitManager(SeedHitManager * shm,Sequence * t,int posi,int posj,int score);
#define add_SeedHit_SeedHitManager Wise2_add_SeedHit_SeedHitManager


/* Function:  populate_HSP_from_SeedHitManager(query,shm,para,p)
 *
 * Descrip:    Chooses which Seeds to use for HSPs. If there are more than factor x max_results
 *             cases, restricted it to factor x max_results scored by seed hits. Then askes
 *             for two diagonal seeds within a wobble factor to start the HSP process
 *
 *
 * Arg:        query [UNKN ] Undocumented argument [Sequence *]
 * Arg:          shm [UNKN ] Undocumented argument [SeedHitManager *]
 * Arg:         para [UNKN ] Undocumented argument [HSPScanInterfacePara *]
 * Arg:            p [UNKN ] Undocumented argument [HSPScanPara *]
 *
 * Return [UNKN ]  Undocumented return value [HSPmanager *]
 *
 */
HSPmanager * Wise2_populate_HSP_from_SeedHitManager(Sequence * query,SeedHitManager * shm,HSPScanInterfacePara * para,HSPScanPara * p);
#define populate_HSP_from_SeedHitManager Wise2_populate_HSP_from_SeedHitManager


/* Function:  new_twohit_HSPScanInterface(sli,mat,drop_off,score_cutoff)
 *
 * Descrip:    Builds a 2 hit search model for protein searches
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
HSPScanInterface * Wise2_new_twohit_HSPScanInterface(SeqLookupInterface * sli,CompMat * mat,int drop_off,int score_cutoff);
#define new_twohit_HSPScanInterface Wise2_new_twohit_HSPScanInterface


/* Function:  one_off_two_hit_HSPscan_query_direct(data,seq,para)
 *
 * Descrip:    two hit approach HSPscan 
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 * Arg:         seq [UNKN ] Undocumented argument [Sequence *]
 * Arg:        para [UNKN ] Undocumented argument [HSPScanInterfacePara *]
 *
 * Return [UNKN ]  Undocumented return value [LinearHSPmanager *]
 *
 */
LinearHSPmanager * Wise2_one_off_two_hit_HSPscan_query_direct(void * data,Sequence * seq,HSPScanInterfacePara * para);
#define one_off_two_hit_HSPscan_query_direct Wise2_one_off_two_hit_HSPscan_query_direct


/* Function:  compare_SeedHitSet_score(one,two)
 *
 * Descrip:    internal for score sorting
 *
 *
 * Arg:        one [UNKN ] Undocumented argument [SeedHitSet *]
 * Arg:        two [UNKN ] Undocumented argument [SeedHitSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_compare_SeedHitSet_score(SeedHitSet * one,SeedHitSet * two);
#define compare_SeedHitSet_score Wise2_compare_SeedHitSet_score


/* Function:  compare_SeedHit_diagonal(one,two)
 *
 * Descrip:    internal for diagonal sorting
 *
 *
 * Arg:        one [UNKN ] Undocumented argument [SeedHit *]
 * Arg:        two [UNKN ] Undocumented argument [SeedHit *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_compare_SeedHit_diagonal(SeedHit * one,SeedHit * two);
#define compare_SeedHit_diagonal Wise2_compare_SeedHit_diagonal


/* Function:  sort_SeedHitManager_by_score(shm)
 *
 * Descrip:    Sorts hit managers by score (highest first)
 *
 *
 * Arg:        shm [UNKN ] Undocumented argument [SeedHitManager *]
 *
 */
void Wise2_sort_SeedHitManager_by_score(SeedHitManager * shm);
#define sort_SeedHitManager_by_score Wise2_sort_SeedHitManager_by_score


/* Function:  sort_SeedHitSet_by_diagonal(sh)
 *
 * Descrip:    Sorts SeedHit by diagonal (lowest first)
 *
 *
 * Arg:        sh [UNKN ] Undocumented argument [SeedHitSet *]
 *
 */
void Wise2_sort_SeedHitSet_by_diagonal(SeedHitSet * sh);
#define sort_SeedHitSet_by_diagonal Wise2_sort_SeedHitSet_by_diagonal


/* Function:  new_SeedHitManager(void)
 *
 * Descrip:    Makes an empty SeedHit Manager
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SeedHitManager *]
 *
 */
SeedHitManager * Wise2_new_SeedHitManager(void);
#define new_SeedHitManager Wise2_new_SeedHitManager


/* Function:  delete_SeedHitManager(shm)
 *
 * Descrip:    Free SeedHitManager
 *
 *
 *
 * Arg:        shm [UNKN ] Undocumented argument [SeedHitManager *]
 *
 * Return [UNKN ]  Undocumented return value [SeedHitManager *]
 *
 */
SeedHitManager * Wise2_delete_SeedHitManager(SeedHitManager * shm);
#define delete_SeedHitManager Wise2_delete_SeedHitManager


/* Function:  hard_link_SeedHit(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SeedHit *]
 *
 * Return [UNKN ]  Undocumented return value [SeedHit *]
 *
 */
SeedHit * Wise2_hard_link_SeedHit(SeedHit * obj);
#define hard_link_SeedHit Wise2_hard_link_SeedHit


/* Function:  SeedHit_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SeedHit *]
 *
 */
SeedHit * Wise2_SeedHit_alloc(void);
#define SeedHit_alloc Wise2_SeedHit_alloc


/* Function:  free_SeedHit(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SeedHit *]
 *
 * Return [UNKN ]  Undocumented return value [SeedHit *]
 *
 */
SeedHit * Wise2_free_SeedHit(SeedHit * obj);
#define free_SeedHit Wise2_free_SeedHit


/* Function:  add_SeedHitSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [SeedHitSet *]
 * Arg:        add [OWNER] Object to add to the list [SeedHit *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_SeedHitSet(SeedHitSet * obj,SeedHit * add);
#define add_SeedHitSet Wise2_add_SeedHitSet


/* Function:  flush_SeedHitSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [SeedHitSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_SeedHitSet(SeedHitSet * obj);
#define flush_SeedHitSet Wise2_flush_SeedHitSet


/* Function:  SeedHitSet_alloc_std(void)
 *
 * Descrip:    Equivalent to SeedHitSet_alloc_len(SeedHitSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SeedHitSet *]
 *
 */
SeedHitSet * Wise2_SeedHitSet_alloc_std(void);
#define SeedHitSet_alloc_std Wise2_SeedHitSet_alloc_std


/* Function:  SeedHitSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [SeedHitSet *]
 *
 */
SeedHitSet * Wise2_SeedHitSet_alloc_len(int len);
#define SeedHitSet_alloc_len Wise2_SeedHitSet_alloc_len


/* Function:  hard_link_SeedHitSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SeedHitSet *]
 *
 * Return [UNKN ]  Undocumented return value [SeedHitSet *]
 *
 */
SeedHitSet * Wise2_hard_link_SeedHitSet(SeedHitSet * obj);
#define hard_link_SeedHitSet Wise2_hard_link_SeedHitSet


/* Function:  SeedHitSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SeedHitSet *]
 *
 */
SeedHitSet * Wise2_SeedHitSet_alloc(void);
#define SeedHitSet_alloc Wise2_SeedHitSet_alloc


/* Function:  free_SeedHitSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SeedHitSet *]
 *
 * Return [UNKN ]  Undocumented return value [SeedHitSet *]
 *
 */
SeedHitSet * Wise2_free_SeedHitSet(SeedHitSet * obj);
#define free_SeedHitSet Wise2_free_SeedHitSet


/* Function:  add_SeedHitManager(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [SeedHitManager *]
 * Arg:        add [OWNER] Object to add to the list [SeedHitSet *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_SeedHitManager(SeedHitManager * obj,SeedHitSet * add);
#define add_SeedHitManager Wise2_add_SeedHitManager


/* Function:  flush_SeedHitManager(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [SeedHitManager *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_SeedHitManager(SeedHitManager * obj);
#define flush_SeedHitManager Wise2_flush_SeedHitManager


/* Function:  SeedHitManager_alloc_std(void)
 *
 * Descrip:    Equivalent to SeedHitManager_alloc_len(SeedHitManagerLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SeedHitManager *]
 *
 */
SeedHitManager * Wise2_SeedHitManager_alloc_std(void);
#define SeedHitManager_alloc_std Wise2_SeedHitManager_alloc_std


/* Function:  SeedHitManager_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [SeedHitManager *]
 *
 */
SeedHitManager * Wise2_SeedHitManager_alloc_len(int len);
#define SeedHitManager_alloc_len Wise2_SeedHitManager_alloc_len


/* Function:  hard_link_SeedHitManager(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SeedHitManager *]
 *
 * Return [UNKN ]  Undocumented return value [SeedHitManager *]
 *
 */
SeedHitManager * Wise2_hard_link_SeedHitManager(SeedHitManager * obj);
#define hard_link_SeedHitManager Wise2_hard_link_SeedHitManager


/* Function:  SeedHitManager_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SeedHitManager *]
 *
 */
SeedHitManager * Wise2_SeedHitManager_alloc(void);
#define SeedHitManager_alloc Wise2_SeedHitManager_alloc


/* Function:  free_SeedHitManager(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SeedHitManager *]
 *
 * Return [UNKN ]  Undocumented return value [SeedHitManager *]
 *
 */
SeedHitManager * Wise2_free_SeedHitManager(SeedHitManager * obj);
#define free_SeedHitManager Wise2_free_SeedHitManager


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_SeedHitSet(SeedHit ** list,int i,int j) ;
#define swap_SeedHitSet Wise2_swap_SeedHitSet
void Wise2_qsort_SeedHitSet(SeedHit ** list,int left,int right,int (*comp)(SeedHit * ,SeedHit * ));
#define qsort_SeedHitSet Wise2_qsort_SeedHitSet
void Wise2_sort_SeedHitSet(SeedHitSet * obj,int (*comp)(SeedHit *, SeedHit *));
#define sort_SeedHitSet Wise2_sort_SeedHitSet
boolean Wise2_expand_SeedHitSet(SeedHitSet * obj,int len);
#define expand_SeedHitSet Wise2_expand_SeedHitSet
void Wise2_swap_SeedHitManager(SeedHitSet ** list,int i,int j) ;
#define swap_SeedHitManager Wise2_swap_SeedHitManager
void Wise2_qsort_SeedHitManager(SeedHitSet ** list,int left,int right,int (*comp)(SeedHitSet * ,SeedHitSet * ));
#define qsort_SeedHitManager Wise2_qsort_SeedHitManager
void Wise2_sort_SeedHitManager(SeedHitManager * obj,int (*comp)(SeedHitSet *, SeedHitSet *));
#define sort_SeedHitManager Wise2_sort_SeedHitManager
boolean Wise2_expand_SeedHitManager(SeedHitManager * obj,int len);
#define expand_SeedHitManager Wise2_expand_SeedHitManager

#ifdef _cplusplus
}
#endif

#endif
