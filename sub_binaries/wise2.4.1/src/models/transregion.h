#ifndef DYNAMITEtransregionHEADERFILE
#define DYNAMITEtransregionHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "transfactor.h"
#include "transregiondp.h"

#define TransFactorRegionLISTLENGTH 128
#define TransFactorRegionSetLISTLENGTH 128

#define TRANSREGION_PARA_DP 93
#define TRANSREGION_PARA_WINDOW 94



struct Wise2_TransFactorRegion {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int start;   
    int end;     
    TransFactorMatch ** match;   
    int len;/* len for above match  */ 
    int maxlen; /* maxlen for above match */ 
    double density;  
    double bits_score;   
    } ;  
/* TransFactorRegion defined */ 
#ifndef DYNAMITE_DEFINED_TransFactorRegion
typedef struct Wise2_TransFactorRegion Wise2_TransFactorRegion;
#define TransFactorRegion Wise2_TransFactorRegion
#define DYNAMITE_DEFINED_TransFactorRegion
#endif


struct Wise2_TransFactorRegionSet {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    TransFactorRegion ** region;     
    int len;/* len for above region  */ 
    int maxlen; /* maxlen for above region */ 
    Sequence * target;   
    } ;  
/* TransFactorRegionSet defined */ 
#ifndef DYNAMITE_DEFINED_TransFactorRegionSet
typedef struct Wise2_TransFactorRegionSet Wise2_TransFactorRegionSet;
#define TransFactorRegionSet Wise2_TransFactorRegionSet
#define DYNAMITE_DEFINED_TransFactorRegionSet
#endif


struct Wise2_TransFactorRegionPara {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int type;    
    double min_density;  
    int min_window;  
    double in_region_prob;   
    double out_region_prob;  
    double in_cost;  
    int hmm_min_window;  
    double gc_region_ratio;  
    } ;  
/* TransFactorRegionPara defined */ 
#ifndef DYNAMITE_DEFINED_TransFactorRegionPara
typedef struct Wise2_TransFactorRegionPara Wise2_TransFactorRegionPara;
#define TransFactorRegionPara Wise2_TransFactorRegionPara
#define DYNAMITE_DEFINED_TransFactorRegionPara
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  add_TransFactorRegion(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [TransFactorRegion *]
 * Arg:        add [OWNER] Object to add to the list [TransFactorMatch *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_TransFactorRegion(TransFactorRegion * obj,TransFactorMatch * add);
#define add_TransFactorRegion Wise2_add_TransFactorRegion


/* Function:  flush_TransFactorRegion(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [TransFactorRegion *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_TransFactorRegion(TransFactorRegion * obj);
#define flush_TransFactorRegion Wise2_flush_TransFactorRegion


/* Function:  TransFactorRegion_alloc_std(void)
 *
 * Descrip:    Equivalent to TransFactorRegion_alloc_len(TransFactorRegionLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TransFactorRegion *]
 *
 */
TransFactorRegion * Wise2_TransFactorRegion_alloc_std(void);
#define TransFactorRegion_alloc_std Wise2_TransFactorRegion_alloc_std


/* Function:  TransFactorRegion_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [TransFactorRegion *]
 *
 */
TransFactorRegion * Wise2_TransFactorRegion_alloc_len(int len);
#define TransFactorRegion_alloc_len Wise2_TransFactorRegion_alloc_len


/* Function:  hard_link_TransFactorRegion(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [TransFactorRegion *]
 *
 * Return [UNKN ]  Undocumented return value [TransFactorRegion *]
 *
 */
TransFactorRegion * Wise2_hard_link_TransFactorRegion(TransFactorRegion * obj);
#define hard_link_TransFactorRegion Wise2_hard_link_TransFactorRegion


/* Function:  TransFactorRegion_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TransFactorRegion *]
 *
 */
TransFactorRegion * Wise2_TransFactorRegion_alloc(void);
#define TransFactorRegion_alloc Wise2_TransFactorRegion_alloc


/* Function:  free_TransFactorRegion(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [TransFactorRegion *]
 *
 * Return [UNKN ]  Undocumented return value [TransFactorRegion *]
 *
 */
TransFactorRegion * Wise2_free_TransFactorRegion(TransFactorRegion * obj);
#define free_TransFactorRegion Wise2_free_TransFactorRegion


/* Function:  add_TransFactorRegionSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [TransFactorRegionSet *]
 * Arg:        add [OWNER] Object to add to the list [TransFactorRegion *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_TransFactorRegionSet(TransFactorRegionSet * obj,TransFactorRegion * add);
#define add_TransFactorRegionSet Wise2_add_TransFactorRegionSet


/* Function:  flush_TransFactorRegionSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [TransFactorRegionSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_TransFactorRegionSet(TransFactorRegionSet * obj);
#define flush_TransFactorRegionSet Wise2_flush_TransFactorRegionSet


/* Function:  TransFactorRegionSet_alloc_std(void)
 *
 * Descrip:    Equivalent to TransFactorRegionSet_alloc_len(TransFactorRegionSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TransFactorRegionSet *]
 *
 */
TransFactorRegionSet * Wise2_TransFactorRegionSet_alloc_std(void);
#define TransFactorRegionSet_alloc_std Wise2_TransFactorRegionSet_alloc_std


/* Function:  TransFactorRegionSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [TransFactorRegionSet *]
 *
 */
TransFactorRegionSet * Wise2_TransFactorRegionSet_alloc_len(int len);
#define TransFactorRegionSet_alloc_len Wise2_TransFactorRegionSet_alloc_len


/* Function:  hard_link_TransFactorRegionSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [TransFactorRegionSet *]
 *
 * Return [UNKN ]  Undocumented return value [TransFactorRegionSet *]
 *
 */
TransFactorRegionSet * Wise2_hard_link_TransFactorRegionSet(TransFactorRegionSet * obj);
#define hard_link_TransFactorRegionSet Wise2_hard_link_TransFactorRegionSet


/* Function:  TransFactorRegionSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TransFactorRegionSet *]
 *
 */
TransFactorRegionSet * Wise2_TransFactorRegionSet_alloc(void);
#define TransFactorRegionSet_alloc Wise2_TransFactorRegionSet_alloc


/* Function:  free_TransFactorRegionSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [TransFactorRegionSet *]
 *
 * Return [UNKN ]  Undocumented return value [TransFactorRegionSet *]
 *
 */
TransFactorRegionSet * Wise2_free_TransFactorRegionSet(TransFactorRegionSet * obj);
#define free_TransFactorRegionSet Wise2_free_TransFactorRegionSet


/* Function:  hard_link_TransFactorRegionPara(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [TransFactorRegionPara *]
 *
 * Return [UNKN ]  Undocumented return value [TransFactorRegionPara *]
 *
 */
TransFactorRegionPara * Wise2_hard_link_TransFactorRegionPara(TransFactorRegionPara * obj);
#define hard_link_TransFactorRegionPara Wise2_hard_link_TransFactorRegionPara


/* Function:  TransFactorRegionPara_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TransFactorRegionPara *]
 *
 */
TransFactorRegionPara * Wise2_TransFactorRegionPara_alloc(void);
#define TransFactorRegionPara_alloc Wise2_TransFactorRegionPara_alloc


/* Function:  free_TransFactorRegionPara(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [TransFactorRegionPara *]
 *
 * Return [UNKN ]  Undocumented return value [TransFactorRegionPara *]
 *
 */
TransFactorRegionPara * Wise2_free_TransFactorRegionPara(TransFactorRegionPara * obj);
#define free_TransFactorRegionPara Wise2_free_TransFactorRegionPara


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
void Wise2_show_help_TransFactorRegionPara(FILE * ofp);
#define show_help_TransFactorRegionPara Wise2_show_help_TransFactorRegionPara
TransFactorRegionPara * Wise2_new_TransFactorRegionPara_from_argv(int * argc,char ** argv);
#define new_TransFactorRegionPara_from_argv Wise2_new_TransFactorRegionPara_from_argv
void Wise2_show_TransFactorRegionSet(TransFactorRegionSet * tfrs,FILE * ofp);
#define show_TransFactorRegionSet Wise2_show_TransFactorRegionSet
TransFactorRegionSet * Wise2_new_TransFactorRegionSet(TransFactorMatchSet * tfms,TransFactorRegionPara * tfrp,DPRunImpl * dpri);
#define new_TransFactorRegionSet Wise2_new_TransFactorRegionSet
TransFactorRegionSet * Wise2_new_dp_TransFactorRegionSet(TransFactorMatchSet * tfms,TransFactorRegionPara * tfrp,DPRunImpl * dpri);
#define new_dp_TransFactorRegionSet Wise2_new_dp_TransFactorRegionSet
TransFactorRegionSet * Wise2_new_window_TransFactorRegionSet(TransFactorMatchSet * tfms,TransFactorRegionPara * tfrp);
#define new_window_TransFactorRegionSet Wise2_new_window_TransFactorRegionSet


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_TransFactorRegion(TransFactorMatch ** list,int i,int j) ;
#define swap_TransFactorRegion Wise2_swap_TransFactorRegion
void Wise2_qsort_TransFactorRegion(TransFactorMatch ** list,int left,int right,int (*comp)(TransFactorMatch * ,TransFactorMatch * ));
#define qsort_TransFactorRegion Wise2_qsort_TransFactorRegion
void Wise2_sort_TransFactorRegion(TransFactorRegion * obj,int (*comp)(TransFactorMatch *, TransFactorMatch *));
#define sort_TransFactorRegion Wise2_sort_TransFactorRegion
boolean Wise2_expand_TransFactorRegion(TransFactorRegion * obj,int len);
#define expand_TransFactorRegion Wise2_expand_TransFactorRegion
void Wise2_swap_TransFactorRegionSet(TransFactorRegion ** list,int i,int j) ;
#define swap_TransFactorRegionSet Wise2_swap_TransFactorRegionSet
void Wise2_qsort_TransFactorRegionSet(TransFactorRegion ** list,int left,int right,int (*comp)(TransFactorRegion * ,TransFactorRegion * ));
#define qsort_TransFactorRegionSet Wise2_qsort_TransFactorRegionSet
void Wise2_sort_TransFactorRegionSet(TransFactorRegionSet * obj,int (*comp)(TransFactorRegion *, TransFactorRegion *));
#define sort_TransFactorRegionSet Wise2_sort_TransFactorRegionSet
boolean Wise2_expand_TransFactorRegionSet(TransFactorRegionSet * obj,int len);
#define expand_TransFactorRegionSet Wise2_expand_TransFactorRegionSet

#ifdef _cplusplus
}
#endif

#endif
