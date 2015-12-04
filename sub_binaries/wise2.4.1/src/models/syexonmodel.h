#ifndef DYNAMITEsyexonmodelHEADERFILE
#define DYNAMITEsyexonmodelHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"

#define SyExonModelLISTLENGTH 128
#define SyExonScoreLISTLENGTH 128

struct Wise2_SyExon {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Probability exit_prob;   
    Probability stay_prob;   
    Probability prev_prob;   
    } ;  
/* SyExon defined */ 
#ifndef DYNAMITE_DEFINED_SyExon
typedef struct Wise2_SyExon Wise2_SyExon;
#define SyExon Wise2_SyExon
#define DYNAMITE_DEFINED_SyExon
#endif


struct Wise2_SyExonScoreUnit {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Score exit_score;    
    Score stay_score;    
    } ;  
/* SyExonScoreUnit defined */ 
#ifndef DYNAMITE_DEFINED_SyExonScoreUnit
typedef struct Wise2_SyExonScoreUnit Wise2_SyExonScoreUnit;
#define SyExonScoreUnit Wise2_SyExonScoreUnit
#define DYNAMITE_DEFINED_SyExonScoreUnit
#endif


struct Wise2_SyExonModel {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    SyExon ** exon;  
    int len;/* len for above exon  */ 
    int maxlen; /* maxlen for above exon */ 
    } ;  
/* SyExonModel defined */ 
#ifndef DYNAMITE_DEFINED_SyExonModel
typedef struct Wise2_SyExonModel Wise2_SyExonModel;
#define SyExonModel Wise2_SyExonModel
#define DYNAMITE_DEFINED_SyExonModel
#endif


struct Wise2_SyExonScore {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    SyExonScoreUnit ** exon;     
    int len;/* len for above exon  */ 
    int maxlen; /* maxlen for above exon */ 
    } ;  
/* SyExonScore defined */ 
#ifndef DYNAMITE_DEFINED_SyExonScore
typedef struct Wise2_SyExonScore Wise2_SyExonScore;
#define SyExonScore Wise2_SyExonScore
#define DYNAMITE_DEFINED_SyExonScore
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_SyExon(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SyExon *]
 *
 * Return [UNKN ]  Undocumented return value [SyExon *]
 *
 */
SyExon * Wise2_hard_link_SyExon(SyExon * obj);
#define hard_link_SyExon Wise2_hard_link_SyExon


/* Function:  SyExon_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SyExon *]
 *
 */
SyExon * Wise2_SyExon_alloc(void);
#define SyExon_alloc Wise2_SyExon_alloc


/* Function:  free_SyExon(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SyExon *]
 *
 * Return [UNKN ]  Undocumented return value [SyExon *]
 *
 */
SyExon * Wise2_free_SyExon(SyExon * obj);
#define free_SyExon Wise2_free_SyExon


/* Function:  hard_link_SyExonScoreUnit(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SyExonScoreUnit *]
 *
 * Return [UNKN ]  Undocumented return value [SyExonScoreUnit *]
 *
 */
SyExonScoreUnit * Wise2_hard_link_SyExonScoreUnit(SyExonScoreUnit * obj);
#define hard_link_SyExonScoreUnit Wise2_hard_link_SyExonScoreUnit


/* Function:  SyExonScoreUnit_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SyExonScoreUnit *]
 *
 */
SyExonScoreUnit * Wise2_SyExonScoreUnit_alloc(void);
#define SyExonScoreUnit_alloc Wise2_SyExonScoreUnit_alloc


/* Function:  free_SyExonScoreUnit(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SyExonScoreUnit *]
 *
 * Return [UNKN ]  Undocumented return value [SyExonScoreUnit *]
 *
 */
SyExonScoreUnit * Wise2_free_SyExonScoreUnit(SyExonScoreUnit * obj);
#define free_SyExonScoreUnit Wise2_free_SyExonScoreUnit


/* Function:  add_SyExonModel(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [SyExonModel *]
 * Arg:        add [OWNER] Object to add to the list [SyExon *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_SyExonModel(SyExonModel * obj,SyExon * add);
#define add_SyExonModel Wise2_add_SyExonModel


/* Function:  flush_SyExonModel(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [SyExonModel *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_SyExonModel(SyExonModel * obj);
#define flush_SyExonModel Wise2_flush_SyExonModel


/* Function:  SyExonModel_alloc_std(void)
 *
 * Descrip:    Equivalent to SyExonModel_alloc_len(SyExonModelLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SyExonModel *]
 *
 */
SyExonModel * Wise2_SyExonModel_alloc_std(void);
#define SyExonModel_alloc_std Wise2_SyExonModel_alloc_std


/* Function:  SyExonModel_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [SyExonModel *]
 *
 */
SyExonModel * Wise2_SyExonModel_alloc_len(int len);
#define SyExonModel_alloc_len Wise2_SyExonModel_alloc_len


/* Function:  hard_link_SyExonModel(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SyExonModel *]
 *
 * Return [UNKN ]  Undocumented return value [SyExonModel *]
 *
 */
SyExonModel * Wise2_hard_link_SyExonModel(SyExonModel * obj);
#define hard_link_SyExonModel Wise2_hard_link_SyExonModel


/* Function:  SyExonModel_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SyExonModel *]
 *
 */
SyExonModel * Wise2_SyExonModel_alloc(void);
#define SyExonModel_alloc Wise2_SyExonModel_alloc


/* Function:  free_SyExonModel(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SyExonModel *]
 *
 * Return [UNKN ]  Undocumented return value [SyExonModel *]
 *
 */
SyExonModel * Wise2_free_SyExonModel(SyExonModel * obj);
#define free_SyExonModel Wise2_free_SyExonModel


/* Function:  add_SyExonScore(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [SyExonScore *]
 * Arg:        add [OWNER] Object to add to the list [SyExonScoreUnit *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_SyExonScore(SyExonScore * obj,SyExonScoreUnit * add);
#define add_SyExonScore Wise2_add_SyExonScore


/* Function:  flush_SyExonScore(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [SyExonScore *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_SyExonScore(SyExonScore * obj);
#define flush_SyExonScore Wise2_flush_SyExonScore


/* Function:  SyExonScore_alloc_std(void)
 *
 * Descrip:    Equivalent to SyExonScore_alloc_len(SyExonScoreLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SyExonScore *]
 *
 */
SyExonScore * Wise2_SyExonScore_alloc_std(void);
#define SyExonScore_alloc_std Wise2_SyExonScore_alloc_std


/* Function:  SyExonScore_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [SyExonScore *]
 *
 */
SyExonScore * Wise2_SyExonScore_alloc_len(int len);
#define SyExonScore_alloc_len Wise2_SyExonScore_alloc_len


/* Function:  hard_link_SyExonScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SyExonScore *]
 *
 * Return [UNKN ]  Undocumented return value [SyExonScore *]
 *
 */
SyExonScore * Wise2_hard_link_SyExonScore(SyExonScore * obj);
#define hard_link_SyExonScore Wise2_hard_link_SyExonScore


/* Function:  SyExonScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SyExonScore *]
 *
 */
SyExonScore * Wise2_SyExonScore_alloc(void);
#define SyExonScore_alloc Wise2_SyExonScore_alloc


/* Function:  free_SyExonScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SyExonScore *]
 *
 * Return [UNKN ]  Undocumented return value [SyExonScore *]
 *
 */
SyExonScore * Wise2_free_SyExonScore(SyExonScore * obj);
#define free_SyExonScore Wise2_free_SyExonScore


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
SyExonScore * Wise2_SyExonScore_flat_model(int start,int end,Probability exit,Probability final_stay);
#define SyExonScore_flat_model Wise2_SyExonScore_flat_model
SyExonModel * Wise2_SyExonModel_flat_model(int start,int end,Probability exit,Probability final_stay);
#define SyExonModel_flat_model Wise2_SyExonModel_flat_model
SyExonScore * Wise2_SyExonScore_from_SyExonModel(SyExonModel * sym);
#define SyExonScore_from_SyExonModel Wise2_SyExonScore_from_SyExonModel
SyExonScoreUnit * Wise2_SyExonScoreUnit_from_SyExon(SyExon * sye);
#define SyExonScoreUnit_from_SyExon Wise2_SyExonScoreUnit_from_SyExon
void Wise2_dump_SyExonScore(SyExonScore * sc,FILE *ofp);
#define dump_SyExonScore Wise2_dump_SyExonScore


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_SyExonModel(SyExon ** list,int i,int j) ;
#define swap_SyExonModel Wise2_swap_SyExonModel
void Wise2_qsort_SyExonModel(SyExon ** list,int left,int right,int (*comp)(SyExon * ,SyExon * ));
#define qsort_SyExonModel Wise2_qsort_SyExonModel
void Wise2_sort_SyExonModel(SyExonModel * obj,int (*comp)(SyExon *, SyExon *));
#define sort_SyExonModel Wise2_sort_SyExonModel
boolean Wise2_expand_SyExonModel(SyExonModel * obj,int len);
#define expand_SyExonModel Wise2_expand_SyExonModel
void Wise2_swap_SyExonScore(SyExonScoreUnit ** list,int i,int j) ;
#define swap_SyExonScore Wise2_swap_SyExonScore
void Wise2_qsort_SyExonScore(SyExonScoreUnit ** list,int left,int right,int (*comp)(SyExonScoreUnit * ,SyExonScoreUnit * ));
#define qsort_SyExonScore Wise2_qsort_SyExonScore
void Wise2_sort_SyExonScore(SyExonScore * obj,int (*comp)(SyExonScoreUnit *, SyExonScoreUnit *));
#define sort_SyExonScore Wise2_sort_SyExonScore
boolean Wise2_expand_SyExonScore(SyExonScore * obj,int len);
#define expand_SyExonScore Wise2_expand_SyExonScore

#ifdef _cplusplus
}
#endif

#endif
