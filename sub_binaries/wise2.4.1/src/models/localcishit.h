#ifndef DYNAMITElocalcishitHEADERFILE
#define DYNAMITElocalcishitHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"
#include "hsp.h"
#include "localdba.h"
#include "dbadisplay.h"
#include "hitlist.h"

#include "motifmatrix.h"
#include "motifmatrixdp.h"

#define LocalCisHitSetLISTLENGTH 128

typedef enum LocalCisGreedyType {
  LocalCisGreedy_None = 55,
  LocalCisGreedy_Query,
  LocalCisGreedy_Both
} LocalCisGreedyType;

struct Wise2_LocalCisHit {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int start_q;     
    int end_q;   
    int start_t;     
    int end_t;   
    int target_rev;  
    AlnBlock * alb; /*  relative to absolute coords */ 
    int score;   
    Sequence * query;    
    Sequence * target;   
    } ;  
/* LocalCisHit defined */ 
#ifndef DYNAMITE_DEFINED_LocalCisHit
typedef struct Wise2_LocalCisHit Wise2_LocalCisHit;
#define LocalCisHit Wise2_LocalCisHit
#define DYNAMITE_DEFINED_LocalCisHit
#endif


struct Wise2_LocalCisHitSet {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    LocalCisHit ** lch;  
    int len;/* len for above lch  */ 
    int maxlen; /* maxlen for above lch */ 
    } ;  
/* LocalCisHitSet defined */ 
#ifndef DYNAMITE_DEFINED_LocalCisHitSet
typedef struct Wise2_LocalCisHitSet Wise2_LocalCisHitSet;
#define LocalCisHitSet Wise2_LocalCisHitSet
#define DYNAMITE_DEFINED_LocalCisHitSet
#endif


struct Wise2_LocalCisHitSetPara {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    double seed_bit_trigger;     
    int expansion_size;  
    double aln_cutoff;   
    boolean sort_by_score;   
    int max;     
    LocalCisGreedyType type;     
    } ;  
/* LocalCisHitSetPara defined */ 
#ifndef DYNAMITE_DEFINED_LocalCisHitSetPara
typedef struct Wise2_LocalCisHitSetPara Wise2_LocalCisHitSetPara;
#define LocalCisHitSetPara Wise2_LocalCisHitSetPara
#define DYNAMITE_DEFINED_LocalCisHitSetPara
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_LocalCisHit(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [LocalCisHit *]
 *
 * Return [UNKN ]  Undocumented return value [LocalCisHit *]
 *
 */
LocalCisHit * Wise2_hard_link_LocalCisHit(LocalCisHit * obj);
#define hard_link_LocalCisHit Wise2_hard_link_LocalCisHit


/* Function:  LocalCisHit_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LocalCisHit *]
 *
 */
LocalCisHit * Wise2_LocalCisHit_alloc(void);
#define LocalCisHit_alloc Wise2_LocalCisHit_alloc


/* Function:  free_LocalCisHit(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [LocalCisHit *]
 *
 * Return [UNKN ]  Undocumented return value [LocalCisHit *]
 *
 */
LocalCisHit * Wise2_free_LocalCisHit(LocalCisHit * obj);
#define free_LocalCisHit Wise2_free_LocalCisHit


/* Function:  add_LocalCisHitSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [LocalCisHitSet *]
 * Arg:        add [OWNER] Object to add to the list [LocalCisHit *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_LocalCisHitSet(LocalCisHitSet * obj,LocalCisHit * add);
#define add_LocalCisHitSet Wise2_add_LocalCisHitSet


/* Function:  flush_LocalCisHitSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [LocalCisHitSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_LocalCisHitSet(LocalCisHitSet * obj);
#define flush_LocalCisHitSet Wise2_flush_LocalCisHitSet


/* Function:  LocalCisHitSet_alloc_std(void)
 *
 * Descrip:    Equivalent to LocalCisHitSet_alloc_len(LocalCisHitSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LocalCisHitSet *]
 *
 */
LocalCisHitSet * Wise2_LocalCisHitSet_alloc_std(void);
#define LocalCisHitSet_alloc_std Wise2_LocalCisHitSet_alloc_std


/* Function:  LocalCisHitSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [LocalCisHitSet *]
 *
 */
LocalCisHitSet * Wise2_LocalCisHitSet_alloc_len(int len);
#define LocalCisHitSet_alloc_len Wise2_LocalCisHitSet_alloc_len


/* Function:  hard_link_LocalCisHitSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [LocalCisHitSet *]
 *
 * Return [UNKN ]  Undocumented return value [LocalCisHitSet *]
 *
 */
LocalCisHitSet * Wise2_hard_link_LocalCisHitSet(LocalCisHitSet * obj);
#define hard_link_LocalCisHitSet Wise2_hard_link_LocalCisHitSet


/* Function:  LocalCisHitSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LocalCisHitSet *]
 *
 */
LocalCisHitSet * Wise2_LocalCisHitSet_alloc(void);
#define LocalCisHitSet_alloc Wise2_LocalCisHitSet_alloc


/* Function:  free_LocalCisHitSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [LocalCisHitSet *]
 *
 * Return [UNKN ]  Undocumented return value [LocalCisHitSet *]
 *
 */
LocalCisHitSet * Wise2_free_LocalCisHitSet(LocalCisHitSet * obj);
#define free_LocalCisHitSet Wise2_free_LocalCisHitSet


/* Function:  hard_link_LocalCisHitSetPara(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [LocalCisHitSetPara *]
 *
 * Return [UNKN ]  Undocumented return value [LocalCisHitSetPara *]
 *
 */
LocalCisHitSetPara * Wise2_hard_link_LocalCisHitSetPara(LocalCisHitSetPara * obj);
#define hard_link_LocalCisHitSetPara Wise2_hard_link_LocalCisHitSetPara


/* Function:  LocalCisHitSetPara_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LocalCisHitSetPara *]
 *
 */
LocalCisHitSetPara * Wise2_LocalCisHitSetPara_alloc(void);
#define LocalCisHitSetPara_alloc Wise2_LocalCisHitSetPara_alloc


/* Function:  free_LocalCisHitSetPara(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [LocalCisHitSetPara *]
 *
 * Return [UNKN ]  Undocumented return value [LocalCisHitSetPara *]
 *
 */
LocalCisHitSetPara * Wise2_free_LocalCisHitSetPara(LocalCisHitSetPara * obj);
#define free_LocalCisHitSetPara Wise2_free_LocalCisHitSetPara


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
HitList * Wise2_HitList_from_LocalCisHitSet(LocalCisHitSet * in);
#define HitList_from_LocalCisHitSet Wise2_HitList_from_LocalCisHitSet
LocalCisHitSet * Wise2_expand_to_subhits_LocalCisHitSet(LocalCisHitSet * in);
#define expand_to_subhits_LocalCisHitSet Wise2_expand_to_subhits_LocalCisHitSet
void Wise2_show_help_LocalCisHitSetPara(FILE * ofp);
#define show_help_LocalCisHitSetPara Wise2_show_help_LocalCisHitSetPara
LocalCisHitSetPara * Wise2_new_LocalCisHitSetPara_from_argv(int * argc,char ** argv);
#define new_LocalCisHitSetPara_from_argv Wise2_new_LocalCisHitSetPara_from_argv
void Wise2_show_pretty_LocalCisHitSet(LocalCisHitSet * lchs,FILE * ofp);
#define show_pretty_LocalCisHitSet Wise2_show_pretty_LocalCisHitSet
void Wise2_show_summary_LocalCisHitSet(LocalCisHitSet * lchs,FILE * ofp);
#define show_summary_LocalCisHitSet Wise2_show_summary_LocalCisHitSet
LocalCisHitSet * Wise2_greedy_weed_LocalCisHitSet(LocalCisHitSet * set,LocalCisHitSetPara *p);
#define greedy_weed_LocalCisHitSet Wise2_greedy_weed_LocalCisHitSet
void Wise2_sort_LocalCisHitSet_by_score(LocalCisHitSet * set);
#define sort_LocalCisHitSet_by_score Wise2_sort_LocalCisHitSet_by_score
void Wise2_sort_LocalCisHitSet_by_start(LocalCisHitSet * set);
#define sort_LocalCisHitSet_by_start Wise2_sort_LocalCisHitSet_by_start
int Wise2_compare_LocalCisHit_score(LocalCisHit * a,LocalCisHit * b);
#define compare_LocalCisHit_score Wise2_compare_LocalCisHit_score
int Wise2_compare_LocalCisHit_start(LocalCisHit * a,LocalCisHit * b);
#define compare_LocalCisHit_start Wise2_compare_LocalCisHit_start
int Wise2_is_query_overlap_LocalCisHit(LocalCisHit * a,LocalCisHit * b);
#define is_query_overlap_LocalCisHit Wise2_is_query_overlap_LocalCisHit
int Wise2_is_target_overlap_LocalCisHit(LocalCisHit * a,LocalCisHit * b);
#define is_target_overlap_LocalCisHit Wise2_is_target_overlap_LocalCisHit
int Wise2_is_consumed_HSP(HSP * a,int q_start,int q_end,int t_start,int t_end);
#define is_consumed_HSP Wise2_is_consumed_HSP
LocalCisHitSet * Wise2_make_LocalCisHitSet(Sequence * query,Sequence * target,Sequence * target_rev,HSPset * forward,HSPset * reverse,LocalCisHitSetPara * p,LocalCisHitScore * scorepara,TransFactorMatchSet * tfms_query,TransFactorMatchSet * tfms_target,TransFactorMatchSet * tfms_target_rev,MotifMatrixScore * mms,boolean use_motif,DPRunImpl * dpri);
#define make_LocalCisHitSet Wise2_make_LocalCisHitSet
LocalCisHit * Wise2_make_LocalCisHit(Sequence * query,Sequence * target,int is_reversed,int guess_q_start,int guess_q_end,int guess_t_start,int guess_t_end,LocalCisHitScore * score,DPRunImpl * dpri);
#define make_LocalCisHit Wise2_make_LocalCisHit
LocalCisHit * Wise2_make_motif_LocalCisHit(Sequence * query,Sequence * target,int is_reversed,int guess_q_start,int guess_q_end,int guess_t_start,int guess_t_end,TransFactorMatchSet * tfms_query,TransFactorMatchSet * tfms_target,MotifMatrixScore * mms,DPRunImpl * dpri);
#define make_motif_LocalCisHit Wise2_make_motif_LocalCisHit


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_LocalCisHitSet(LocalCisHit ** list,int i,int j) ;
#define swap_LocalCisHitSet Wise2_swap_LocalCisHitSet
void Wise2_qsort_LocalCisHitSet(LocalCisHit ** list,int left,int right,int (*comp)(LocalCisHit * ,LocalCisHit * ));
#define qsort_LocalCisHitSet Wise2_qsort_LocalCisHitSet
void Wise2_sort_LocalCisHitSet(LocalCisHitSet * obj,int (*comp)(LocalCisHit *, LocalCisHit *));
#define sort_LocalCisHitSet Wise2_sort_LocalCisHitSet
boolean Wise2_expand_LocalCisHitSet(LocalCisHitSet * obj,int len);
#define expand_LocalCisHitSet Wise2_expand_LocalCisHitSet

#ifdef _cplusplus
}
#endif

#endif
