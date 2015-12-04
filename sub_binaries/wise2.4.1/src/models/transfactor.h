#ifndef DYNAMITEtransfactorHEADERFILE
#define DYNAMITEtransfactorHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "pwmdna.h"
#include "hitlist.h"


#define TransFactorSetLISTLENGTH 64

#define TransFactorMatchSetLISTLENGTH 64
#define TransFactorMatchSetComparaLISTLENGTH 8

#define TFCOMPARA_OVERLAP_PRECISE 37
#define TFCOMPARA_OVERLAP_OVERLAP 38
#define TFCOMPARA_OVERLAP_REGION  39

#define TFM_ABSOLUTE 84
#define TFM_RELATIVE 85
#define TFM_RELATIVE_MIXED 86


struct Wise2_TransFactor {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * name;     
    SeqAlign * seed;     
    pwmDNA * pwm;    
    SeqAlign * full;     
    double max_prob;     
    double min_prob;     
    } ;  
/* TransFactor defined */ 
#ifndef DYNAMITE_DEFINED_TransFactor
typedef struct Wise2_TransFactor Wise2_TransFactor;
#define TransFactor Wise2_TransFactor
#define DYNAMITE_DEFINED_TransFactor
#endif


struct Wise2_TransFactorSet {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    TransFactor ** factor;   
    int len;/* len for above factor  */ 
    int maxlen; /* maxlen for above factor */ 
    } ;  
/* TransFactorSet defined */ 
#ifndef DYNAMITE_DEFINED_TransFactorSet
typedef struct Wise2_TransFactorSet Wise2_TransFactorSet;
#define TransFactorSet Wise2_TransFactorSet
#define DYNAMITE_DEFINED_TransFactorSet
#endif


struct Wise2_TransFactorMatch {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int start;   
    int end;     
    char strand;     
    double bit_score;    
    TransFactor * factor;    
    } ;  
/* TransFactorMatch defined */ 
#ifndef DYNAMITE_DEFINED_TransFactorMatch
typedef struct Wise2_TransFactorMatch Wise2_TransFactorMatch;
#define TransFactorMatch Wise2_TransFactorMatch
#define DYNAMITE_DEFINED_TransFactorMatch
#endif


struct Wise2_TransFactorMatchSet {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    TransFactorMatch ** match;   
    int len;/* len for above match  */ 
    int maxlen; /* maxlen for above match */ 
    Sequence * target;   
    } ;  
/* TransFactorMatchSet defined */ 
#ifndef DYNAMITE_DEFINED_TransFactorMatchSet
typedef struct Wise2_TransFactorMatchSet Wise2_TransFactorMatchSet;
#define TransFactorMatchSet Wise2_TransFactorMatchSet
#define DYNAMITE_DEFINED_TransFactorMatchSet
#endif


struct Wise2_TransFactorMatchSetCompara {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    SeqAlign * sa;   
    TransFactorMatchSet ** tfms;     
    int len;/* len for above tfms  */ 
    int maxlen; /* maxlen for above tfms */ 
    TransFactorMatchSet *  overall;  
    } ;  
/* TransFactorMatchSetCompara defined */ 
#ifndef DYNAMITE_DEFINED_TransFactorMatchSetCompara
typedef struct Wise2_TransFactorMatchSetCompara Wise2_TransFactorMatchSetCompara;
#define TransFactorMatchSetCompara Wise2_TransFactorMatchSetCompara
#define DYNAMITE_DEFINED_TransFactorMatchSetCompara
#endif


struct Wise2_TransFactorComparaPara {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int overlap_type;    
    int overlap_len;     
    int missing_seq;     
    } ;  
/* TransFactorComparaPara defined */ 
#ifndef DYNAMITE_DEFINED_TransFactorComparaPara
typedef struct Wise2_TransFactorComparaPara Wise2_TransFactorComparaPara;
#define TransFactorComparaPara Wise2_TransFactorComparaPara
#define DYNAMITE_DEFINED_TransFactorComparaPara
#endif


struct Wise2_TransFactorBuildPara {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    RandomModelDNA * rnd_dna;    
    double pseudo_count;     
    boolean warn_on_small_seq;   
    } ;  
/* TransFactorBuildPara defined */ 
#ifndef DYNAMITE_DEFINED_TransFactorBuildPara
typedef struct Wise2_TransFactorBuildPara Wise2_TransFactorBuildPara;
#define TransFactorBuildPara Wise2_TransFactorBuildPara
#define DYNAMITE_DEFINED_TransFactorBuildPara
#endif


struct Wise2_TransFactorMatchPara {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    double min_bits;     
    double min_relative;     
    double relative_prob;    
    double relative_prob_bits;   
    char type;   
    int allow_N;     
    } ;  
/* TransFactorMatchPara defined */ 
#ifndef DYNAMITE_DEFINED_TransFactorMatchPara
typedef struct Wise2_TransFactorMatchPara Wise2_TransFactorMatchPara;
#define TransFactorMatchPara Wise2_TransFactorMatchPara
#define DYNAMITE_DEFINED_TransFactorMatchPara
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_TransFactor(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [TransFactor *]
 *
 * Return [UNKN ]  Undocumented return value [TransFactor *]
 *
 */
TransFactor * Wise2_hard_link_TransFactor(TransFactor * obj);
#define hard_link_TransFactor Wise2_hard_link_TransFactor


/* Function:  TransFactor_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TransFactor *]
 *
 */
TransFactor * Wise2_TransFactor_alloc(void);
#define TransFactor_alloc Wise2_TransFactor_alloc


/* Function:  free_TransFactor(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [TransFactor *]
 *
 * Return [UNKN ]  Undocumented return value [TransFactor *]
 *
 */
TransFactor * Wise2_free_TransFactor(TransFactor * obj);
#define free_TransFactor Wise2_free_TransFactor


/* Function:  add_TransFactorSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [TransFactorSet *]
 * Arg:        add [OWNER] Object to add to the list [TransFactor *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_TransFactorSet(TransFactorSet * obj,TransFactor * add);
#define add_TransFactorSet Wise2_add_TransFactorSet


/* Function:  flush_TransFactorSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [TransFactorSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_TransFactorSet(TransFactorSet * obj);
#define flush_TransFactorSet Wise2_flush_TransFactorSet


/* Function:  TransFactorSet_alloc_std(void)
 *
 * Descrip:    Equivalent to TransFactorSet_alloc_len(TransFactorSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TransFactorSet *]
 *
 */
TransFactorSet * Wise2_TransFactorSet_alloc_std(void);
#define TransFactorSet_alloc_std Wise2_TransFactorSet_alloc_std


/* Function:  TransFactorSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [TransFactorSet *]
 *
 */
TransFactorSet * Wise2_TransFactorSet_alloc_len(int len);
#define TransFactorSet_alloc_len Wise2_TransFactorSet_alloc_len


/* Function:  hard_link_TransFactorSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [TransFactorSet *]
 *
 * Return [UNKN ]  Undocumented return value [TransFactorSet *]
 *
 */
TransFactorSet * Wise2_hard_link_TransFactorSet(TransFactorSet * obj);
#define hard_link_TransFactorSet Wise2_hard_link_TransFactorSet


/* Function:  TransFactorSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TransFactorSet *]
 *
 */
TransFactorSet * Wise2_TransFactorSet_alloc(void);
#define TransFactorSet_alloc Wise2_TransFactorSet_alloc


/* Function:  free_TransFactorSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [TransFactorSet *]
 *
 * Return [UNKN ]  Undocumented return value [TransFactorSet *]
 *
 */
TransFactorSet * Wise2_free_TransFactorSet(TransFactorSet * obj);
#define free_TransFactorSet Wise2_free_TransFactorSet


/* Function:  hard_link_TransFactorMatch(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [TransFactorMatch *]
 *
 * Return [UNKN ]  Undocumented return value [TransFactorMatch *]
 *
 */
TransFactorMatch * Wise2_hard_link_TransFactorMatch(TransFactorMatch * obj);
#define hard_link_TransFactorMatch Wise2_hard_link_TransFactorMatch


/* Function:  TransFactorMatch_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TransFactorMatch *]
 *
 */
TransFactorMatch * Wise2_TransFactorMatch_alloc(void);
#define TransFactorMatch_alloc Wise2_TransFactorMatch_alloc


/* Function:  free_TransFactorMatch(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [TransFactorMatch *]
 *
 * Return [UNKN ]  Undocumented return value [TransFactorMatch *]
 *
 */
TransFactorMatch * Wise2_free_TransFactorMatch(TransFactorMatch * obj);
#define free_TransFactorMatch Wise2_free_TransFactorMatch


/* Function:  add_TransFactorMatchSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [TransFactorMatchSet *]
 * Arg:        add [OWNER] Object to add to the list [TransFactorMatch *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_TransFactorMatchSet(TransFactorMatchSet * obj,TransFactorMatch * add);
#define add_TransFactorMatchSet Wise2_add_TransFactorMatchSet


/* Function:  flush_TransFactorMatchSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [TransFactorMatchSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_TransFactorMatchSet(TransFactorMatchSet * obj);
#define flush_TransFactorMatchSet Wise2_flush_TransFactorMatchSet


/* Function:  TransFactorMatchSet_alloc_std(void)
 *
 * Descrip:    Equivalent to TransFactorMatchSet_alloc_len(TransFactorMatchSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TransFactorMatchSet *]
 *
 */
TransFactorMatchSet * Wise2_TransFactorMatchSet_alloc_std(void);
#define TransFactorMatchSet_alloc_std Wise2_TransFactorMatchSet_alloc_std


/* Function:  TransFactorMatchSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [TransFactorMatchSet *]
 *
 */
TransFactorMatchSet * Wise2_TransFactorMatchSet_alloc_len(int len);
#define TransFactorMatchSet_alloc_len Wise2_TransFactorMatchSet_alloc_len


/* Function:  hard_link_TransFactorMatchSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [TransFactorMatchSet *]
 *
 * Return [UNKN ]  Undocumented return value [TransFactorMatchSet *]
 *
 */
TransFactorMatchSet * Wise2_hard_link_TransFactorMatchSet(TransFactorMatchSet * obj);
#define hard_link_TransFactorMatchSet Wise2_hard_link_TransFactorMatchSet


/* Function:  TransFactorMatchSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TransFactorMatchSet *]
 *
 */
TransFactorMatchSet * Wise2_TransFactorMatchSet_alloc(void);
#define TransFactorMatchSet_alloc Wise2_TransFactorMatchSet_alloc


/* Function:  free_TransFactorMatchSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [TransFactorMatchSet *]
 *
 * Return [UNKN ]  Undocumented return value [TransFactorMatchSet *]
 *
 */
TransFactorMatchSet * Wise2_free_TransFactorMatchSet(TransFactorMatchSet * obj);
#define free_TransFactorMatchSet Wise2_free_TransFactorMatchSet


/* Function:  add_TransFactorMatchSetCompara(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [TransFactorMatchSetCompara *]
 * Arg:        add [OWNER] Object to add to the list [TransFactorMatchSet *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_TransFactorMatchSetCompara(TransFactorMatchSetCompara * obj,TransFactorMatchSet * add);
#define add_TransFactorMatchSetCompara Wise2_add_TransFactorMatchSetCompara


/* Function:  flush_TransFactorMatchSetCompara(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [TransFactorMatchSetCompara *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_TransFactorMatchSetCompara(TransFactorMatchSetCompara * obj);
#define flush_TransFactorMatchSetCompara Wise2_flush_TransFactorMatchSetCompara


/* Function:  TransFactorMatchSetCompara_alloc_std(void)
 *
 * Descrip:    Equivalent to TransFactorMatchSetCompara_alloc_len(TransFactorMatchSetComparaLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TransFactorMatchSetCompara *]
 *
 */
TransFactorMatchSetCompara * Wise2_TransFactorMatchSetCompara_alloc_std(void);
#define TransFactorMatchSetCompara_alloc_std Wise2_TransFactorMatchSetCompara_alloc_std


/* Function:  TransFactorMatchSetCompara_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [TransFactorMatchSetCompara *]
 *
 */
TransFactorMatchSetCompara * Wise2_TransFactorMatchSetCompara_alloc_len(int len);
#define TransFactorMatchSetCompara_alloc_len Wise2_TransFactorMatchSetCompara_alloc_len


/* Function:  hard_link_TransFactorMatchSetCompara(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [TransFactorMatchSetCompara *]
 *
 * Return [UNKN ]  Undocumented return value [TransFactorMatchSetCompara *]
 *
 */
TransFactorMatchSetCompara * Wise2_hard_link_TransFactorMatchSetCompara(TransFactorMatchSetCompara * obj);
#define hard_link_TransFactorMatchSetCompara Wise2_hard_link_TransFactorMatchSetCompara


/* Function:  TransFactorMatchSetCompara_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TransFactorMatchSetCompara *]
 *
 */
TransFactorMatchSetCompara * Wise2_TransFactorMatchSetCompara_alloc(void);
#define TransFactorMatchSetCompara_alloc Wise2_TransFactorMatchSetCompara_alloc


/* Function:  free_TransFactorMatchSetCompara(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [TransFactorMatchSetCompara *]
 *
 * Return [UNKN ]  Undocumented return value [TransFactorMatchSetCompara *]
 *
 */
TransFactorMatchSetCompara * Wise2_free_TransFactorMatchSetCompara(TransFactorMatchSetCompara * obj);
#define free_TransFactorMatchSetCompara Wise2_free_TransFactorMatchSetCompara


/* Function:  hard_link_TransFactorComparaPara(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [TransFactorComparaPara *]
 *
 * Return [UNKN ]  Undocumented return value [TransFactorComparaPara *]
 *
 */
TransFactorComparaPara * Wise2_hard_link_TransFactorComparaPara(TransFactorComparaPara * obj);
#define hard_link_TransFactorComparaPara Wise2_hard_link_TransFactorComparaPara


/* Function:  TransFactorComparaPara_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TransFactorComparaPara *]
 *
 */
TransFactorComparaPara * Wise2_TransFactorComparaPara_alloc(void);
#define TransFactorComparaPara_alloc Wise2_TransFactorComparaPara_alloc


/* Function:  free_TransFactorComparaPara(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [TransFactorComparaPara *]
 *
 * Return [UNKN ]  Undocumented return value [TransFactorComparaPara *]
 *
 */
TransFactorComparaPara * Wise2_free_TransFactorComparaPara(TransFactorComparaPara * obj);
#define free_TransFactorComparaPara Wise2_free_TransFactorComparaPara


/* Function:  hard_link_TransFactorBuildPara(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [TransFactorBuildPara *]
 *
 * Return [UNKN ]  Undocumented return value [TransFactorBuildPara *]
 *
 */
TransFactorBuildPara * Wise2_hard_link_TransFactorBuildPara(TransFactorBuildPara * obj);
#define hard_link_TransFactorBuildPara Wise2_hard_link_TransFactorBuildPara


/* Function:  TransFactorBuildPara_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TransFactorBuildPara *]
 *
 */
TransFactorBuildPara * Wise2_TransFactorBuildPara_alloc(void);
#define TransFactorBuildPara_alloc Wise2_TransFactorBuildPara_alloc


/* Function:  free_TransFactorBuildPara(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [TransFactorBuildPara *]
 *
 * Return [UNKN ]  Undocumented return value [TransFactorBuildPara *]
 *
 */
TransFactorBuildPara * Wise2_free_TransFactorBuildPara(TransFactorBuildPara * obj);
#define free_TransFactorBuildPara Wise2_free_TransFactorBuildPara


/* Function:  hard_link_TransFactorMatchPara(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [TransFactorMatchPara *]
 *
 * Return [UNKN ]  Undocumented return value [TransFactorMatchPara *]
 *
 */
TransFactorMatchPara * Wise2_hard_link_TransFactorMatchPara(TransFactorMatchPara * obj);
#define hard_link_TransFactorMatchPara Wise2_hard_link_TransFactorMatchPara


/* Function:  TransFactorMatchPara_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TransFactorMatchPara *]
 *
 */
TransFactorMatchPara * Wise2_TransFactorMatchPara_alloc(void);
#define TransFactorMatchPara_alloc Wise2_TransFactorMatchPara_alloc


/* Function:  free_TransFactorMatchPara(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [TransFactorMatchPara *]
 *
 * Return [UNKN ]  Undocumented return value [TransFactorMatchPara *]
 *
 */
TransFactorMatchPara * Wise2_free_TransFactorMatchPara(TransFactorMatchPara * obj);
#define free_TransFactorMatchPara Wise2_free_TransFactorMatchPara


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
void Wise2_show_context_match_TransFactorMatchSetCompara(TransFactorMatchSetCompara * tfmsc,int context,FILE *ofp);
#define show_context_match_TransFactorMatchSetCompara Wise2_show_context_match_TransFactorMatchSetCompara
TransFactorSet * Wise2_circular_permuted_TransFactorSet(TransFactorSet * in,int rotate_number);
#define circular_permuted_TransFactorSet Wise2_circular_permuted_TransFactorSet
void Wise2_show_TransFactorMatchSet(TransFactorMatchSet * tfms,FILE * ofp);
#define show_TransFactorMatchSet Wise2_show_TransFactorMatchSet
void Wise2_show_help_TransFactorMatchPara(FILE * ofp);
#define show_help_TransFactorMatchPara Wise2_show_help_TransFactorMatchPara
void Wise2_show_help_TransFactorComparaPara(FILE * ofp);
#define show_help_TransFactorComparaPara Wise2_show_help_TransFactorComparaPara
TransFactorComparaPara * Wise2_new_TransFactorComparaPara_from_argv(int * argc,char ** argv);
#define new_TransFactorComparaPara_from_argv Wise2_new_TransFactorComparaPara_from_argv
TransFactorMatchSetCompara * Wise2_calculate_TransFactorMatchSetCompara(SeqAlign * sa,TransFactorSet * tfs,TransFactorMatchPara * match_para,TransFactorComparaPara * para);
#define calculate_TransFactorMatchSetCompara Wise2_calculate_TransFactorMatchSetCompara
void Wise2_show_help_TransFactorBuildPara(FILE * ofp);
#define show_help_TransFactorBuildPara Wise2_show_help_TransFactorBuildPara
TransFactorBuildPara * Wise2_new_TransFactorBuildPara_from_argv(int * argc,char ** argv);
#define new_TransFactorBuildPara_from_argv Wise2_new_TransFactorBuildPara_from_argv
TransFactorMatchPara * Wise2_new_TransFactorMatchPara_from_argv(int * argc,char ** argv);
#define new_TransFactorMatchPara_from_argv Wise2_new_TransFactorMatchPara_from_argv
boolean Wise2_build_TransFactorSet(TransFactorSet * tfs,TransFactorBuildPara * p);
#define build_TransFactorSet Wise2_build_TransFactorSet
boolean Wise2_build_pwm_TransFactor(TransFactor * tf,TransFactorBuildPara * p);
#define build_pwm_TransFactor Wise2_build_pwm_TransFactor
TransFactorMatchSet * Wise2_calculate_TransFactorMatchSet(Sequence * seq,TransFactorSet * tfs,TransFactorMatchPara * p);
#define calculate_TransFactorMatchSet Wise2_calculate_TransFactorMatchSet
double Wise2_min_prob_TransFactor(TransFactor * tf);
#define min_prob_TransFactor Wise2_min_prob_TransFactor
double Wise2_max_prob_TransFactor(TransFactor * tf);
#define max_prob_TransFactor Wise2_max_prob_TransFactor
void Wise2_write_TransFactorSet(TransFactorSet * tfs,FILE * ofp);
#define write_TransFactorSet Wise2_write_TransFactorSet
void Wise2_write_TransFactor(TransFactor * tf,FILE * ofp);
#define write_TransFactor Wise2_write_TransFactor
TransFactorSet * Wise2_read_TransFactorSet_file(char * filename);
#define read_TransFactorSet_file Wise2_read_TransFactorSet_file
TransFactorSet * Wise2_read_ben_IUPAC_TransFactorSet_file(char * filename);
#define read_ben_IUPAC_TransFactorSet_file Wise2_read_ben_IUPAC_TransFactorSet_file
TransFactorSet * Wise2_read_ben_IUPAC_TransFactorSet(FILE * ifp);
#define read_ben_IUPAC_TransFactorSet Wise2_read_ben_IUPAC_TransFactorSet
TransFactorSet * Wise2_read_laurence_TransFactorSet_file(char * filename);
#define read_laurence_TransFactorSet_file Wise2_read_laurence_TransFactorSet_file
TransFactorSet * Wise2_read_laurence_TransFactorSet(FILE * ifp);
#define read_laurence_TransFactorSet Wise2_read_laurence_TransFactorSet
TransFactorSet * Wise2_read_TransFactorSet(FILE * ifp);
#define read_TransFactorSet Wise2_read_TransFactorSet
TransFactor * Wise2_read_TransFactor(char * line,FILE * ifp);
#define read_TransFactor Wise2_read_TransFactor
int Wise2_compare_start_TransFactorMatch(TransFactorMatch * a,TransFactorMatch * b);
#define compare_start_TransFactorMatch Wise2_compare_start_TransFactorMatch
void Wise2_sort_by_start_TransFactorMatchSet(TransFactorMatchSet * t);
#define sort_by_start_TransFactorMatchSet Wise2_sort_by_start_TransFactorMatchSet


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_TransFactorSet(TransFactor ** list,int i,int j) ;
#define swap_TransFactorSet Wise2_swap_TransFactorSet
void Wise2_qsort_TransFactorSet(TransFactor ** list,int left,int right,int (*comp)(TransFactor * ,TransFactor * ));
#define qsort_TransFactorSet Wise2_qsort_TransFactorSet
void Wise2_sort_TransFactorSet(TransFactorSet * obj,int (*comp)(TransFactor *, TransFactor *));
#define sort_TransFactorSet Wise2_sort_TransFactorSet
boolean Wise2_expand_TransFactorSet(TransFactorSet * obj,int len);
#define expand_TransFactorSet Wise2_expand_TransFactorSet
void Wise2_swap_TransFactorMatchSet(TransFactorMatch ** list,int i,int j) ;
#define swap_TransFactorMatchSet Wise2_swap_TransFactorMatchSet
void Wise2_qsort_TransFactorMatchSet(TransFactorMatch ** list,int left,int right,int (*comp)(TransFactorMatch * ,TransFactorMatch * ));
#define qsort_TransFactorMatchSet Wise2_qsort_TransFactorMatchSet
void Wise2_sort_TransFactorMatchSet(TransFactorMatchSet * obj,int (*comp)(TransFactorMatch *, TransFactorMatch *));
#define sort_TransFactorMatchSet Wise2_sort_TransFactorMatchSet
boolean Wise2_expand_TransFactorMatchSet(TransFactorMatchSet * obj,int len);
#define expand_TransFactorMatchSet Wise2_expand_TransFactorMatchSet
void Wise2_swap_TransFactorMatchSetCompara(TransFactorMatchSet ** list,int i,int j) ;
#define swap_TransFactorMatchSetCompara Wise2_swap_TransFactorMatchSetCompara
void Wise2_qsort_TransFactorMatchSetCompara(TransFactorMatchSet ** list,int left,int right,int (*comp)(TransFactorMatchSet * ,TransFactorMatchSet * ));
#define qsort_TransFactorMatchSetCompara Wise2_qsort_TransFactorMatchSetCompara
void Wise2_sort_TransFactorMatchSetCompara(TransFactorMatchSetCompara * obj,int (*comp)(TransFactorMatchSet *, TransFactorMatchSet *));
#define sort_TransFactorMatchSetCompara Wise2_sort_TransFactorMatchSetCompara
boolean Wise2_expand_TransFactorMatchSetCompara(TransFactorMatchSetCompara * obj,int len);
#define expand_TransFactorMatchSetCompara Wise2_expand_TransFactorMatchSetCompara

#ifdef _cplusplus
}
#endif

#endif
