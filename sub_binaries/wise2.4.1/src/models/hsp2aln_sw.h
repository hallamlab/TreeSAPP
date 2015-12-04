#ifndef DYNAMITEhsp2aln_swHEADERFILE
#define DYNAMITEhsp2aln_swHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "sw_wrap.h"
#include "hitlist.h"
#include "hsp.h"




#define Hsp2AlnHelperLISTLENGTH 20

#ifdef PTHREAD
#include "pthread.h"


struct hsp2aln_thread_manager {
  pthread_mutex_t input_lock;
  pthread_mutex_t output_lock;
  pthread_t * pool;
  int thread_size;
  int current_pos;
  LinearHSPmanager * input;
  HitList * output;
  DPRunImpl * dpri;
  struct Wise2_HSPset2HitPairPara * para;
  int topscore;
};
#endif

struct Wise2_HSPset2HitPairPara {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int hsp_width;   
    int hsp_length;  
    int no_subalns;  
    int no_hitalns;  
    boolean best_hit;    
    double perc_hit_dropoff;     
    boolean debug;   
    int poor_score_factor;   
    int poor_score;  
    } ;  
/* HSPset2HitPairPara defined */ 
#ifndef DYNAMITE_DEFINED_HSPset2HitPairPara
typedef struct Wise2_HSPset2HitPairPara Wise2_HSPset2HitPairPara;
#define HSPset2HitPairPara Wise2_HSPset2HitPairPara
#define DYNAMITE_DEFINED_HSPset2HitPairPara
#endif


struct Wise2_Hsp2AlnHelper {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    DPEnvelope ** dpenv;     
    int len;/* len for above dpenv  */ 
    int maxlen; /* maxlen for above dpenv */ 
    } ;  
/* Hsp2AlnHelper defined */ 
#ifndef DYNAMITE_DEFINED_Hsp2AlnHelper
typedef struct Wise2_Hsp2AlnHelper Wise2_Hsp2AlnHelper;
#define Hsp2AlnHelper Wise2_Hsp2AlnHelper
#define DYNAMITE_DEFINED_Hsp2AlnHelper
#endif


struct Wise2_Hsp2AlnPara {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int align_size;  
    } ;  
/* Hsp2AlnPara defined */ 
#ifndef DYNAMITE_DEFINED_Hsp2AlnPara
typedef struct Wise2_Hsp2AlnPara Wise2_Hsp2AlnPara;
#define Hsp2AlnPara Wise2_Hsp2AlnPara
#define DYNAMITE_DEFINED_Hsp2AlnPara
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_HSPset2HitPairPara(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [HSPset2HitPairPara *]
 *
 * Return [UNKN ]  Undocumented return value [HSPset2HitPairPara *]
 *
 */
HSPset2HitPairPara * Wise2_hard_link_HSPset2HitPairPara(HSPset2HitPairPara * obj);
#define hard_link_HSPset2HitPairPara Wise2_hard_link_HSPset2HitPairPara


/* Function:  HSPset2HitPairPara_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [HSPset2HitPairPara *]
 *
 */
HSPset2HitPairPara * Wise2_HSPset2HitPairPara_alloc(void);
#define HSPset2HitPairPara_alloc Wise2_HSPset2HitPairPara_alloc


/* Function:  free_HSPset2HitPairPara(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [HSPset2HitPairPara *]
 *
 * Return [UNKN ]  Undocumented return value [HSPset2HitPairPara *]
 *
 */
HSPset2HitPairPara * Wise2_free_HSPset2HitPairPara(HSPset2HitPairPara * obj);
#define free_HSPset2HitPairPara Wise2_free_HSPset2HitPairPara


/* Function:  add_Hsp2AlnHelper(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Hsp2AlnHelper *]
 * Arg:        add [OWNER] Object to add to the list [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_Hsp2AlnHelper(Hsp2AlnHelper * obj,DPEnvelope * add);
#define add_Hsp2AlnHelper Wise2_add_Hsp2AlnHelper


/* Function:  flush_Hsp2AlnHelper(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [Hsp2AlnHelper *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_Hsp2AlnHelper(Hsp2AlnHelper * obj);
#define flush_Hsp2AlnHelper Wise2_flush_Hsp2AlnHelper


/* Function:  Hsp2AlnHelper_alloc_std(void)
 *
 * Descrip:    Equivalent to Hsp2AlnHelper_alloc_len(Hsp2AlnHelperLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Hsp2AlnHelper *]
 *
 */
Hsp2AlnHelper * Wise2_Hsp2AlnHelper_alloc_std(void);
#define Hsp2AlnHelper_alloc_std Wise2_Hsp2AlnHelper_alloc_std


/* Function:  Hsp2AlnHelper_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [Hsp2AlnHelper *]
 *
 */
Hsp2AlnHelper * Wise2_Hsp2AlnHelper_alloc_len(int len);
#define Hsp2AlnHelper_alloc_len Wise2_Hsp2AlnHelper_alloc_len


/* Function:  hard_link_Hsp2AlnHelper(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Hsp2AlnHelper *]
 *
 * Return [UNKN ]  Undocumented return value [Hsp2AlnHelper *]
 *
 */
Hsp2AlnHelper * Wise2_hard_link_Hsp2AlnHelper(Hsp2AlnHelper * obj);
#define hard_link_Hsp2AlnHelper Wise2_hard_link_Hsp2AlnHelper


/* Function:  Hsp2AlnHelper_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Hsp2AlnHelper *]
 *
 */
Hsp2AlnHelper * Wise2_Hsp2AlnHelper_alloc(void);
#define Hsp2AlnHelper_alloc Wise2_Hsp2AlnHelper_alloc


/* Function:  free_Hsp2AlnHelper(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Hsp2AlnHelper *]
 *
 * Return [UNKN ]  Undocumented return value [Hsp2AlnHelper *]
 *
 */
Hsp2AlnHelper * Wise2_free_Hsp2AlnHelper(Hsp2AlnHelper * obj);
#define free_Hsp2AlnHelper Wise2_free_Hsp2AlnHelper


/* Function:  hard_link_Hsp2AlnPara(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Hsp2AlnPara *]
 *
 * Return [UNKN ]  Undocumented return value [Hsp2AlnPara *]
 *
 */
Hsp2AlnPara * Wise2_hard_link_Hsp2AlnPara(Hsp2AlnPara * obj);
#define hard_link_Hsp2AlnPara Wise2_hard_link_Hsp2AlnPara


/* Function:  Hsp2AlnPara_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Hsp2AlnPara *]
 *
 */
Hsp2AlnPara * Wise2_Hsp2AlnPara_alloc(void);
#define Hsp2AlnPara_alloc Wise2_Hsp2AlnPara_alloc


/* Function:  free_Hsp2AlnPara(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Hsp2AlnPara *]
 *
 * Return [UNKN ]  Undocumented return value [Hsp2AlnPara *]
 *
 */
Hsp2AlnPara * Wise2_free_Hsp2AlnPara(Hsp2AlnPara * obj);
#define free_Hsp2AlnPara Wise2_free_Hsp2AlnPara


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
HitList * Wise2_HitList_from_LinearHSPmanager_heuristic_threaded(LinearHSPmanager * lm,DPRunImpl * dpri,int thr_no,HSPset2HitPairPara * para);
#define HitList_from_LinearHSPmanager_heuristic_threaded Wise2_HitList_from_LinearHSPmanager_heuristic_threaded
void * Wise2_worker_thread_LM2HitList(void * p);
#define worker_thread_LM2HitList Wise2_worker_thread_LM2HitList
HitList * Wise2_HitList_from_LinearHSPmanager_heuristic_threaded(LinearHSPmanager * lm,DPRunImpl * dpri,int thr_no,HSPset2HitPairPara * para);
#define HitList_from_LinearHSPmanager_heuristic_threaded Wise2_HitList_from_LinearHSPmanager_heuristic_threaded
HitList * Wise2_HitList_from_LinearHSPmanager_heuristic(LinearHSPmanager * lm,DPRunImpl * dpri,HSPset2HitPairPara * para);
#define HitList_from_LinearHSPmanager_heuristic Wise2_HitList_from_LinearHSPmanager_heuristic
HitPair * Wise2_HitPair_from_HSPset_heuristic(HSPset * set,DPRunImpl * dpri,CompMat * mat,HSPset2HitPairPara *p);
#define HitPair_from_HSPset_heuristic Wise2_HitPair_from_HSPset_heuristic
Hsp2AlnHelper * Wise2_build_HSP2AlnHelper(HSPset * set,int width,int tail,int min_score,int small_factor);
#define build_HSP2AlnHelper Wise2_build_HSP2AlnHelper
void Wise2_show_help_HSPset2HitPairPara(FILE * ofp);
#define show_help_HSPset2HitPairPara Wise2_show_help_HSPset2HitPairPara
HSPset2HitPairPara * Wise2_new_HSPset2HitPairPara_from_argv(int * argc,char ** argv);
#define new_HSPset2HitPairPara_from_argv Wise2_new_HSPset2HitPairPara_from_argv


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_Hsp2AlnHelper(DPEnvelope ** list,int i,int j) ;
#define swap_Hsp2AlnHelper Wise2_swap_Hsp2AlnHelper
void Wise2_qsort_Hsp2AlnHelper(DPEnvelope ** list,int left,int right,int (*comp)(DPEnvelope * ,DPEnvelope * ));
#define qsort_Hsp2AlnHelper Wise2_qsort_Hsp2AlnHelper
void Wise2_sort_Hsp2AlnHelper(Hsp2AlnHelper * obj,int (*comp)(DPEnvelope *, DPEnvelope *));
#define sort_Hsp2AlnHelper Wise2_sort_Hsp2AlnHelper
boolean Wise2_expand_Hsp2AlnHelper(Hsp2AlnHelper * obj,int len);
#define expand_Hsp2AlnHelper Wise2_expand_Hsp2AlnHelper

#ifdef _cplusplus
}
#endif

#endif
