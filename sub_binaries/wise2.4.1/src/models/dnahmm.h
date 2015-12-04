#ifndef DYNAMITEdnahmmHEADERFILE
#define DYNAMITEdnahmmHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"


enum dhmm_trans_type {
  DHMM_MATCH2MATCH  = 0,
  DHMM_MATCH2INSERT,  
  DHMM_MATCH2DELETE,
  DHMM_INSERT2MATCH,
  DHMM_INSERT2INSERT,  
  DHMM_INSERT2DELETE,
  DHMM_DELETE2MATCH,
  DHMM_DELETE2INSERT,  
  DHMM_DELETE2DELETE,  
  DHMM_START2MATCH,    
  DHMM_START2INSERT,   
  DHMM_START2DELETE,   
  DHMM_MATCH2END,      
  DHMM_INSERT2END,     
  DHMM_DELETE2END,
  DHMM_TRANSITION_LEN
};

#define DNA_HMM_ALPHABET 5
#define DnaHmmProbLISTLENGTH 64
#define DnaHmmScoreLISTLENGTH 64

struct Wise2_DnaHmmProbUnit {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Probability match[DNA_HMM_ALPHABET];     
    Probability insert[DNA_HMM_ALPHABET];    
    Probability transition[DHMM_TRANSITION_LEN];     
    } ;  
/* DnaHmmProbUnit defined */ 
#ifndef DYNAMITE_DEFINED_DnaHmmProbUnit
typedef struct Wise2_DnaHmmProbUnit Wise2_DnaHmmProbUnit;
#define DnaHmmProbUnit Wise2_DnaHmmProbUnit
#define DYNAMITE_DEFINED_DnaHmmProbUnit
#endif


struct Wise2_DnaHmmProb {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    DnaHmmProbUnit ** unit;  
    int len;/* len for above unit  */ 
    int maxlen; /* maxlen for above unit */ 
    SeqAlign * ref_align;    
    char * consensus;    
    } ;  
/* DnaHmmProb defined */ 
#ifndef DYNAMITE_DEFINED_DnaHmmProb
typedef struct Wise2_DnaHmmProb Wise2_DnaHmmProb;
#define DnaHmmProb Wise2_DnaHmmProb
#define DYNAMITE_DEFINED_DnaHmmProb
#endif


struct Wise2_DnaHmmScoreUnit {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Score match[DNA_HMM_ALPHABET];   
    Score insert[DNA_HMM_ALPHABET];  
    Score transition[DHMM_TRANSITION_LEN];   
    } ;  
/* DnaHmmScoreUnit defined */ 
#ifndef DYNAMITE_DEFINED_DnaHmmScoreUnit
typedef struct Wise2_DnaHmmScoreUnit Wise2_DnaHmmScoreUnit;
#define DnaHmmScoreUnit Wise2_DnaHmmScoreUnit
#define DYNAMITE_DEFINED_DnaHmmScoreUnit
#endif


struct Wise2_DnaHmmScore {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    DnaHmmScoreUnit ** unit;     
    int len;/* len for above unit  */ 
    int maxlen; /* maxlen for above unit */ 
    SeqAlign * ref_align;    
    } ;  
/* DnaHmmScore defined */ 
#ifndef DYNAMITE_DEFINED_DnaHmmScore
typedef struct Wise2_DnaHmmScore Wise2_DnaHmmScore;
#define DnaHmmScore Wise2_DnaHmmScore
#define DYNAMITE_DEFINED_DnaHmmScore
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_DnaHmmProbUnit(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DnaHmmProbUnit *]
 *
 * Return [UNKN ]  Undocumented return value [DnaHmmProbUnit *]
 *
 */
DnaHmmProbUnit * Wise2_hard_link_DnaHmmProbUnit(DnaHmmProbUnit * obj);
#define hard_link_DnaHmmProbUnit Wise2_hard_link_DnaHmmProbUnit


/* Function:  DnaHmmProbUnit_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaHmmProbUnit *]
 *
 */
DnaHmmProbUnit * Wise2_DnaHmmProbUnit_alloc(void);
#define DnaHmmProbUnit_alloc Wise2_DnaHmmProbUnit_alloc


/* Function:  free_DnaHmmProbUnit(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DnaHmmProbUnit *]
 *
 * Return [UNKN ]  Undocumented return value [DnaHmmProbUnit *]
 *
 */
DnaHmmProbUnit * Wise2_free_DnaHmmProbUnit(DnaHmmProbUnit * obj);
#define free_DnaHmmProbUnit Wise2_free_DnaHmmProbUnit


/* Function:  add_DnaHmmProb(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [DnaHmmProb *]
 * Arg:        add [OWNER] Object to add to the list [DnaHmmProbUnit *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_DnaHmmProb(DnaHmmProb * obj,DnaHmmProbUnit * add);
#define add_DnaHmmProb Wise2_add_DnaHmmProb


/* Function:  flush_DnaHmmProb(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [DnaHmmProb *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_DnaHmmProb(DnaHmmProb * obj);
#define flush_DnaHmmProb Wise2_flush_DnaHmmProb


/* Function:  DnaHmmProb_alloc_std(void)
 *
 * Descrip:    Equivalent to DnaHmmProb_alloc_len(DnaHmmProbLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaHmmProb *]
 *
 */
DnaHmmProb * Wise2_DnaHmmProb_alloc_std(void);
#define DnaHmmProb_alloc_std Wise2_DnaHmmProb_alloc_std


/* Function:  DnaHmmProb_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [DnaHmmProb *]
 *
 */
DnaHmmProb * Wise2_DnaHmmProb_alloc_len(int len);
#define DnaHmmProb_alloc_len Wise2_DnaHmmProb_alloc_len


/* Function:  hard_link_DnaHmmProb(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DnaHmmProb *]
 *
 * Return [UNKN ]  Undocumented return value [DnaHmmProb *]
 *
 */
DnaHmmProb * Wise2_hard_link_DnaHmmProb(DnaHmmProb * obj);
#define hard_link_DnaHmmProb Wise2_hard_link_DnaHmmProb


/* Function:  DnaHmmProb_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaHmmProb *]
 *
 */
DnaHmmProb * Wise2_DnaHmmProb_alloc(void);
#define DnaHmmProb_alloc Wise2_DnaHmmProb_alloc


/* Function:  free_DnaHmmProb(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DnaHmmProb *]
 *
 * Return [UNKN ]  Undocumented return value [DnaHmmProb *]
 *
 */
DnaHmmProb * Wise2_free_DnaHmmProb(DnaHmmProb * obj);
#define free_DnaHmmProb Wise2_free_DnaHmmProb


/* Function:  hard_link_DnaHmmScoreUnit(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DnaHmmScoreUnit *]
 *
 * Return [UNKN ]  Undocumented return value [DnaHmmScoreUnit *]
 *
 */
DnaHmmScoreUnit * Wise2_hard_link_DnaHmmScoreUnit(DnaHmmScoreUnit * obj);
#define hard_link_DnaHmmScoreUnit Wise2_hard_link_DnaHmmScoreUnit


/* Function:  DnaHmmScoreUnit_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaHmmScoreUnit *]
 *
 */
DnaHmmScoreUnit * Wise2_DnaHmmScoreUnit_alloc(void);
#define DnaHmmScoreUnit_alloc Wise2_DnaHmmScoreUnit_alloc


/* Function:  free_DnaHmmScoreUnit(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DnaHmmScoreUnit *]
 *
 * Return [UNKN ]  Undocumented return value [DnaHmmScoreUnit *]
 *
 */
DnaHmmScoreUnit * Wise2_free_DnaHmmScoreUnit(DnaHmmScoreUnit * obj);
#define free_DnaHmmScoreUnit Wise2_free_DnaHmmScoreUnit


/* Function:  add_DnaHmmScore(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [DnaHmmScore *]
 * Arg:        add [OWNER] Object to add to the list [DnaHmmScoreUnit *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_DnaHmmScore(DnaHmmScore * obj,DnaHmmScoreUnit * add);
#define add_DnaHmmScore Wise2_add_DnaHmmScore


/* Function:  flush_DnaHmmScore(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [DnaHmmScore *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_DnaHmmScore(DnaHmmScore * obj);
#define flush_DnaHmmScore Wise2_flush_DnaHmmScore


/* Function:  DnaHmmScore_alloc_std(void)
 *
 * Descrip:    Equivalent to DnaHmmScore_alloc_len(DnaHmmScoreLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaHmmScore *]
 *
 */
DnaHmmScore * Wise2_DnaHmmScore_alloc_std(void);
#define DnaHmmScore_alloc_std Wise2_DnaHmmScore_alloc_std


/* Function:  DnaHmmScore_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [DnaHmmScore *]
 *
 */
DnaHmmScore * Wise2_DnaHmmScore_alloc_len(int len);
#define DnaHmmScore_alloc_len Wise2_DnaHmmScore_alloc_len


/* Function:  hard_link_DnaHmmScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DnaHmmScore *]
 *
 * Return [UNKN ]  Undocumented return value [DnaHmmScore *]
 *
 */
DnaHmmScore * Wise2_hard_link_DnaHmmScore(DnaHmmScore * obj);
#define hard_link_DnaHmmScore Wise2_hard_link_DnaHmmScore


/* Function:  DnaHmmScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaHmmScore *]
 *
 */
DnaHmmScore * Wise2_DnaHmmScore_alloc(void);
#define DnaHmmScore_alloc Wise2_DnaHmmScore_alloc


/* Function:  free_DnaHmmScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DnaHmmScore *]
 *
 * Return [UNKN ]  Undocumented return value [DnaHmmScore *]
 *
 */
DnaHmmScore * Wise2_free_DnaHmmScore(DnaHmmScore * obj);
#define free_DnaHmmScore Wise2_free_DnaHmmScore


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
DnaHmmScore * Wise2_DnaHmmScore_from_DnaHmmProb(DnaHmmProb * dhp);
#define DnaHmmScore_from_DnaHmmProb Wise2_DnaHmmScore_from_DnaHmmProb
void Wise2_make_consensus_DnaHmmProb(DnaHmmProb * dhp);
#define make_consensus_DnaHmmProb Wise2_make_consensus_DnaHmmProb
DnaHmmScoreUnit * Wise2_DnaHmmScoreUnit_from_DnaHmmProbUnit(DnaHmmProbUnit * dpu);
#define DnaHmmScoreUnit_from_DnaHmmProbUnit Wise2_DnaHmmScoreUnit_from_DnaHmmProbUnit
DnaHmmProb * Wise2_new_DnaHmmProb_from_SeqAlign_ungapped(SeqAlign * sa,double simple_pseudocount);
#define new_DnaHmmProb_from_SeqAlign_ungapped Wise2_new_DnaHmmProb_from_SeqAlign_ungapped
DnaHmmProbUnit * Wise2_new_DnaHmmProbUnit_from_ColumnCount_ungapped(ColumnCount * cc,double simple_pseudocount);
#define new_DnaHmmProbUnit_from_ColumnCount_ungapped Wise2_new_DnaHmmProbUnit_from_ColumnCount_ungapped
void Wise2_set_N_DnaHmmProb(DnaHmmProb * dhp,Probability basen);
#define set_N_DnaHmmProb Wise2_set_N_DnaHmmProb
void Wise2_fold_RandomModelDNA_DnaHmmProb(DnaHmmProb * dhp,RandomModelDNA * d,Probability rnd_advance);
#define fold_RandomModelDNA_DnaHmmProb Wise2_fold_RandomModelDNA_DnaHmmProb


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_DnaHmmProb(DnaHmmProbUnit ** list,int i,int j) ;
#define swap_DnaHmmProb Wise2_swap_DnaHmmProb
void Wise2_qsort_DnaHmmProb(DnaHmmProbUnit ** list,int left,int right,int (*comp)(DnaHmmProbUnit * ,DnaHmmProbUnit * ));
#define qsort_DnaHmmProb Wise2_qsort_DnaHmmProb
void Wise2_sort_DnaHmmProb(DnaHmmProb * obj,int (*comp)(DnaHmmProbUnit *, DnaHmmProbUnit *));
#define sort_DnaHmmProb Wise2_sort_DnaHmmProb
boolean Wise2_expand_DnaHmmProb(DnaHmmProb * obj,int len);
#define expand_DnaHmmProb Wise2_expand_DnaHmmProb
void Wise2_swap_DnaHmmScore(DnaHmmScoreUnit ** list,int i,int j) ;
#define swap_DnaHmmScore Wise2_swap_DnaHmmScore
void Wise2_qsort_DnaHmmScore(DnaHmmScoreUnit ** list,int left,int right,int (*comp)(DnaHmmScoreUnit * ,DnaHmmScoreUnit * ));
#define qsort_DnaHmmScore Wise2_qsort_DnaHmmScore
void Wise2_sort_DnaHmmScore(DnaHmmScore * obj,int (*comp)(DnaHmmScoreUnit *, DnaHmmScoreUnit *));
#define sort_DnaHmmScore Wise2_sort_DnaHmmScore
boolean Wise2_expand_DnaHmmScore(DnaHmmScore * obj,int len);
#define expand_DnaHmmScore Wise2_expand_DnaHmmScore

#ifdef _cplusplus
}
#endif

#endif
