#ifndef DYNAMITEpwmdnaHEADERFILE
#define DYNAMITEpwmdnaHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "wisebase.h"
#include "seqalign.h"
#include "randommodel.h"
#include "probability.h"

#define pwmDNALISTLENGTH 64
#define pwmDNAScoreLISTLENGTH 64
/* Object pwmColScore
 *
 * Descrip: Actual Scores for a
 *        position in a Score
 *        representation of a PWM
 *
 *
 */
struct Wise2_pwmColScore {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Score emit[5];   
    } ;  
/* pwmColScore defined */ 
#ifndef DYNAMITE_DEFINED_pwmColScore
typedef struct Wise2_pwmColScore Wise2_pwmColScore;
#define pwmColScore Wise2_pwmColScore
#define DYNAMITE_DEFINED_pwmColScore
#endif


/* Object pwmColProb
 *
 * Descrip: Actual probabilities for
 *        a position in a PWM
 *
 *
 */
struct Wise2_pwmColProb {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Probability emit[5];     
    } ;  
/* pwmColProb defined */ 
#ifndef DYNAMITE_DEFINED_pwmColProb
typedef struct Wise2_pwmColProb Wise2_pwmColProb;
#define pwmColProb Wise2_pwmColProb
#define DYNAMITE_DEFINED_pwmColProb
#endif


/* Object pwmDNA
 *
 * Descrip: This structure holds a position
 *        weight matrix as probabilities.
 *
 *        Generally you want to use the Score
 *        based system for the actual scoring,
 *        but this for manipulation of the probabilities
 *
 *        You can build this data structure from a
 *        simple sequence alignment plus a pseudocount
 *
 *
 */
struct Wise2_pwmDNA {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    pwmColProb ** pos;   
    int len;/* len for above pos  */ 
    int maxlen; /* maxlen for above pos */ 
    SeqAlign * ref_align;    
    } ;  
/* pwmDNA defined */ 
#ifndef DYNAMITE_DEFINED_pwmDNA
typedef struct Wise2_pwmDNA Wise2_pwmDNA;
#define pwmDNA Wise2_pwmDNA
#define DYNAMITE_DEFINED_pwmDNA
#endif


/* Object pwmDNAScore
 *
 * Descrip: This structure holds a position
 *        weight matrix as Scores, generally
 *        log-odded to a random model
 *
 *        This is the structure used for scoring.
 *        You make it from a pwmDNA
 *
 *
 */
struct Wise2_pwmDNAScore {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    pwmColScore ** pos;  
    int len;/* len for above pos  */ 
    int maxlen; /* maxlen for above pos */ 
    SeqAlign * ref_align;    
    } ;  
/* pwmDNAScore defined */ 
#ifndef DYNAMITE_DEFINED_pwmDNAScore
typedef struct Wise2_pwmDNAScore Wise2_pwmDNAScore;
#define pwmDNAScore Wise2_pwmDNAScore
#define DYNAMITE_DEFINED_pwmDNAScore
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  max_prob_pwmDNA(in)
 *
 * Descrip:    maximum prob of pwmDNA
 *
 *
 * Arg:        in [UNKN ] Undocumented argument [pwmDNA *]
 *
 * Return [UNKN ]  Undocumented return value [double]
 *
 */
double Wise2_max_prob_pwmDNA(pwmDNA * in);
#define max_prob_pwmDNA Wise2_max_prob_pwmDNA


/* Function:  min_prob_pwmDNA(in)
 *
 * Descrip:    minimum prob of pwmDNA
 *
 *
 * Arg:        in [UNKN ] Undocumented argument [pwmDNA *]
 *
 * Return [UNKN ]  Undocumented return value [double]
 *
 */
double Wise2_min_prob_pwmDNA(pwmDNA * in);
#define min_prob_pwmDNA Wise2_min_prob_pwmDNA


/* Function:  circular_permuted_pwmDNA(in,rotate_number)
 *
 * Descrip:    Provides a rotated pwm for randomisation
 *
 *
 * Arg:                   in [UNKN ] Undocumented argument [pwmDNA *]
 * Arg:        rotate_number [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [pwmDNA *]
 *
 */
pwmDNA * Wise2_circular_permuted_pwmDNA(pwmDNA * in,int rotate_number);
#define circular_permuted_pwmDNA Wise2_circular_permuted_pwmDNA


/* Function:  score_pwmDNAScore_Sequence(pds,s,pos)
 *
 * Descrip:    This gives back a Score from a particular sequence and
 *             position
 *
 *
 * Arg:        pds [UNKN ] Undocumented argument [pwmDNAScore *]
 * Arg:          s [UNKN ] Undocumented argument [Sequence *]
 * Arg:        pos [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [Score]
 *
 */
Score Wise2_score_pwmDNAScore_Sequence(pwmDNAScore * pds,Sequence * s,int pos);
#define score_pwmDNAScore_Sequence Wise2_score_pwmDNAScore_Sequence


/* Function:  score_pwmDNAScore_string(pds,str)
 *
 * Descrip:    This gives back a Score from a particular string
 *
 *
 * Arg:        pds [UNKN ] Undocumented argument [pwmDNAScore *]
 * Arg:        str [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [Score]
 *
 */
Score Wise2_score_pwmDNAScore_string(pwmDNAScore * pds,char * str);
#define score_pwmDNAScore_string Wise2_score_pwmDNAScore_string


/* Function:  prob_pwmDNA_string(pds,str)
 *
 * Descrip:    This gives back a Probability from a particular string
 *
 *
 * Arg:        pds [UNKN ] Undocumented argument [pwmDNA *]
 * Arg:        str [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [Probability]
 *
 */
Probability Wise2_prob_pwmDNA_string(pwmDNA * pds,char * str);
#define prob_pwmDNA_string Wise2_prob_pwmDNA_string


/* Function:  fold_randommodel_pwmDNA(pd,rmd)
 *
 * Descrip:    This folds in a randommodel into pwmDNA
 *
 *
 * Arg:         pd [UNKN ] Undocumented argument [pwmDNA *]
 * Arg:        rmd [UNKN ] Undocumented argument [RandomModelDNA *]
 *
 */
void Wise2_fold_randommodel_pwmDNA(pwmDNA * pd,RandomModelDNA * rmd);
#define fold_randommodel_pwmDNA Wise2_fold_randommodel_pwmDNA


/* Function:  pwmDNA_from_SeqAlign(sa,simple_pseudocount)
 *
 * Descrip:    This function makes a single pwmDNA from
 *             a SeqAlign
 *
 *             FIXME: This DOES NOT handle ambiguity codes well
 *
 *
 * Arg:                        sa [UNKN ] Undocumented argument [SeqAlign *]
 * Arg:        simple_pseudocount [UNKN ] Undocumented argument [double]
 *
 * Return [UNKN ]  Undocumented return value [pwmDNA *]
 *
 */
pwmDNA * Wise2_pwmDNA_from_SeqAlign(SeqAlign * sa,double simple_pseudocount);
#define pwmDNA_from_SeqAlign Wise2_pwmDNA_from_SeqAlign


/* Function:  show_pwmDNA_col(pd,ofp)
 *
 * Descrip:    Shows a columns along the page
 *
 *
 * Arg:         pd [UNKN ] Undocumented argument [pwmDNA *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_show_pwmDNA_col(pwmDNA * pd,FILE * ofp);
#define show_pwmDNA_col Wise2_show_pwmDNA_col


/* Function:  pwmDNAScore_from_pwmDNA_RandomModelDNA(pwm,rmd)
 *
 * Descrip:    This function makes score represention of a
 *             position weight matrix from a probability, 
 *             with a random model folded in on-the-fly
 *
 *
 * Arg:        pwm [UNKN ] Undocumented argument [pwmDNA *]
 * Arg:        rmd [UNKN ] Undocumented argument [RandomModelDNA *]
 *
 * Return [UNKN ]  Undocumented return value [pwmDNAScore *]
 *
 */
pwmDNAScore * Wise2_pwmDNAScore_from_pwmDNA_RandomModelDNA(pwmDNA * pwm,RandomModelDNA * rmd);
#define pwmDNAScore_from_pwmDNA_RandomModelDNA Wise2_pwmDNAScore_from_pwmDNA_RandomModelDNA


/* Function:  pwmDNAScore_from_pwmDNA(pwm)
 *
 * Descrip:    This function makes score represention of a
 *             position weight matrix from a probability
 *
 *
 * Arg:        pwm [UNKN ] Undocumented argument [pwmDNA *]
 *
 * Return [UNKN ]  Undocumented return value [pwmDNAScore *]
 *
 */
pwmDNAScore * Wise2_pwmDNAScore_from_pwmDNA(pwmDNA * pwm);
#define pwmDNAScore_from_pwmDNA Wise2_pwmDNAScore_from_pwmDNA


/* Function:  hard_link_pwmColScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [pwmColScore *]
 *
 * Return [UNKN ]  Undocumented return value [pwmColScore *]
 *
 */
pwmColScore * Wise2_hard_link_pwmColScore(pwmColScore * obj);
#define hard_link_pwmColScore Wise2_hard_link_pwmColScore


/* Function:  pwmColScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [pwmColScore *]
 *
 */
pwmColScore * Wise2_pwmColScore_alloc(void);
#define pwmColScore_alloc Wise2_pwmColScore_alloc


/* Function:  free_pwmColScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [pwmColScore *]
 *
 * Return [UNKN ]  Undocumented return value [pwmColScore *]
 *
 */
pwmColScore * Wise2_free_pwmColScore(pwmColScore * obj);
#define free_pwmColScore Wise2_free_pwmColScore


/* Function:  hard_link_pwmColProb(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [pwmColProb *]
 *
 * Return [UNKN ]  Undocumented return value [pwmColProb *]
 *
 */
pwmColProb * Wise2_hard_link_pwmColProb(pwmColProb * obj);
#define hard_link_pwmColProb Wise2_hard_link_pwmColProb


/* Function:  pwmColProb_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [pwmColProb *]
 *
 */
pwmColProb * Wise2_pwmColProb_alloc(void);
#define pwmColProb_alloc Wise2_pwmColProb_alloc


/* Function:  free_pwmColProb(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [pwmColProb *]
 *
 * Return [UNKN ]  Undocumented return value [pwmColProb *]
 *
 */
pwmColProb * Wise2_free_pwmColProb(pwmColProb * obj);
#define free_pwmColProb Wise2_free_pwmColProb


/* Function:  add_pwmDNA(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [pwmDNA *]
 * Arg:        add [OWNER] Object to add to the list [pwmColProb *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_pwmDNA(pwmDNA * obj,pwmColProb * add);
#define add_pwmDNA Wise2_add_pwmDNA


/* Function:  flush_pwmDNA(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [pwmDNA *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_pwmDNA(pwmDNA * obj);
#define flush_pwmDNA Wise2_flush_pwmDNA


/* Function:  pwmDNA_alloc_std(void)
 *
 * Descrip:    Equivalent to pwmDNA_alloc_len(pwmDNALISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [pwmDNA *]
 *
 */
pwmDNA * Wise2_pwmDNA_alloc_std(void);
#define pwmDNA_alloc_std Wise2_pwmDNA_alloc_std


/* Function:  pwmDNA_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [pwmDNA *]
 *
 */
pwmDNA * Wise2_pwmDNA_alloc_len(int len);
#define pwmDNA_alloc_len Wise2_pwmDNA_alloc_len


/* Function:  hard_link_pwmDNA(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [pwmDNA *]
 *
 * Return [UNKN ]  Undocumented return value [pwmDNA *]
 *
 */
pwmDNA * Wise2_hard_link_pwmDNA(pwmDNA * obj);
#define hard_link_pwmDNA Wise2_hard_link_pwmDNA


/* Function:  pwmDNA_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [pwmDNA *]
 *
 */
pwmDNA * Wise2_pwmDNA_alloc(void);
#define pwmDNA_alloc Wise2_pwmDNA_alloc


/* Function:  free_pwmDNA(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [pwmDNA *]
 *
 * Return [UNKN ]  Undocumented return value [pwmDNA *]
 *
 */
pwmDNA * Wise2_free_pwmDNA(pwmDNA * obj);
#define free_pwmDNA Wise2_free_pwmDNA


/* Function:  add_pwmDNAScore(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [pwmDNAScore *]
 * Arg:        add [OWNER] Object to add to the list [pwmColScore *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_pwmDNAScore(pwmDNAScore * obj,pwmColScore * add);
#define add_pwmDNAScore Wise2_add_pwmDNAScore


/* Function:  flush_pwmDNAScore(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [pwmDNAScore *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_pwmDNAScore(pwmDNAScore * obj);
#define flush_pwmDNAScore Wise2_flush_pwmDNAScore


/* Function:  pwmDNAScore_alloc_std(void)
 *
 * Descrip:    Equivalent to pwmDNAScore_alloc_len(pwmDNAScoreLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [pwmDNAScore *]
 *
 */
pwmDNAScore * Wise2_pwmDNAScore_alloc_std(void);
#define pwmDNAScore_alloc_std Wise2_pwmDNAScore_alloc_std


/* Function:  pwmDNAScore_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [pwmDNAScore *]
 *
 */
pwmDNAScore * Wise2_pwmDNAScore_alloc_len(int len);
#define pwmDNAScore_alloc_len Wise2_pwmDNAScore_alloc_len


/* Function:  hard_link_pwmDNAScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [pwmDNAScore *]
 *
 * Return [UNKN ]  Undocumented return value [pwmDNAScore *]
 *
 */
pwmDNAScore * Wise2_hard_link_pwmDNAScore(pwmDNAScore * obj);
#define hard_link_pwmDNAScore Wise2_hard_link_pwmDNAScore


/* Function:  pwmDNAScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [pwmDNAScore *]
 *
 */
pwmDNAScore * Wise2_pwmDNAScore_alloc(void);
#define pwmDNAScore_alloc Wise2_pwmDNAScore_alloc


/* Function:  free_pwmDNAScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [pwmDNAScore *]
 *
 * Return [UNKN ]  Undocumented return value [pwmDNAScore *]
 *
 */
pwmDNAScore * Wise2_free_pwmDNAScore(pwmDNAScore * obj);
#define free_pwmDNAScore Wise2_free_pwmDNAScore


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
pwmColProb * Wise2_pwmColProb_from_ColumnCount(ColumnCount * cc,double simple_pseudocount);
#define pwmColProb_from_ColumnCount Wise2_pwmColProb_from_ColumnCount
pwmColScore * Wise2_pwmColScore_from_pwmColProb(pwmColProb * p);
#define pwmColScore_from_pwmColProb Wise2_pwmColScore_from_pwmColProb
pwmColScore * Wise2_pwmColScore_from_pwmColProb_rmd(pwmColProb * p,RandomModelDNA  * rmd);
#define pwmColScore_from_pwmColProb_rmd Wise2_pwmColScore_from_pwmColProb_rmd
void Wise2_swap_pwmDNA(pwmColProb ** list,int i,int j) ;
#define swap_pwmDNA Wise2_swap_pwmDNA
void Wise2_qsort_pwmDNA(pwmColProb ** list,int left,int right,int (*comp)(pwmColProb * ,pwmColProb * ));
#define qsort_pwmDNA Wise2_qsort_pwmDNA
void Wise2_sort_pwmDNA(pwmDNA * obj,int (*comp)(pwmColProb *, pwmColProb *));
#define sort_pwmDNA Wise2_sort_pwmDNA
boolean Wise2_expand_pwmDNA(pwmDNA * obj,int len);
#define expand_pwmDNA Wise2_expand_pwmDNA
void Wise2_swap_pwmDNAScore(pwmColScore ** list,int i,int j) ;
#define swap_pwmDNAScore Wise2_swap_pwmDNAScore
void Wise2_qsort_pwmDNAScore(pwmColScore ** list,int left,int right,int (*comp)(pwmColScore * ,pwmColScore * ));
#define qsort_pwmDNAScore Wise2_qsort_pwmDNAScore
void Wise2_sort_pwmDNAScore(pwmDNAScore * obj,int (*comp)(pwmColScore *, pwmColScore *));
#define sort_pwmDNAScore Wise2_sort_pwmDNAScore
boolean Wise2_expand_pwmDNAScore(pwmDNAScore * obj,int len);
#define expand_pwmDNAScore Wise2_expand_pwmDNAScore

#ifdef _cplusplus
}
#endif

#endif
