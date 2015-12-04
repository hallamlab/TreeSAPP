#ifndef DYNAMITEaligngenemodelHEADERFILE
#define DYNAMITEaligngenemodelHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"
#include "pwmdna.h"
#include "genestats.h"

typedef struct align_gene_codon {
	char seq[3];
	double weight;
} AlignGeneCodon;


struct Wise2_NonCodingSimpleModel {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Probability identical;   
    Probability one_off;     
    Probability two_off;     
    } ;  
/* NonCodingSimpleModel defined */ 
#ifndef DYNAMITE_DEFINED_NonCodingSimpleModel
typedef struct Wise2_NonCodingSimpleModel Wise2_NonCodingSimpleModel;
#define NonCodingSimpleModel Wise2_NonCodingSimpleModel
#define DYNAMITE_DEFINED_NonCodingSimpleModel
#endif


struct Wise2_SpliceSiteProb {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    pwmDNA * pwm;    
    int offset;  
    } ;  
/* SpliceSiteProb defined */ 
#ifndef DYNAMITE_DEFINED_SpliceSiteProb
typedef struct Wise2_SpliceSiteProb Wise2_SpliceSiteProb;
#define SpliceSiteProb Wise2_SpliceSiteProb
#define DYNAMITE_DEFINED_SpliceSiteProb
#endif


struct Wise2_AlignGeneColumnStore {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    AlignGeneCodon * codon;  
    int len;     
    } ;  
/* AlignGeneColumnStore defined */ 
#ifndef DYNAMITE_DEFINED_AlignGeneColumnStore
typedef struct Wise2_AlignGeneColumnStore Wise2_AlignGeneColumnStore;
#define AlignGeneColumnStore Wise2_AlignGeneColumnStore
#define DYNAMITE_DEFINED_AlignGeneColumnStore
#endif


struct Wise2_AlignGeneModel {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int len;     
    Probability * forward_coding;    
    Probability * reverse_coding;    
    Probability * splice5_forward;   
    Probability * splice3_forward;   
    Probability * splice5_reverse;   
    Probability * splice3_reverse;   
    SeqAlign * align;    
    Sequence * anchor;   
    double   * change_rate;  
    } ;  
/* AlignGeneModel defined */ 
#ifndef DYNAMITE_DEFINED_AlignGeneModel
typedef struct Wise2_AlignGeneModel Wise2_AlignGeneModel;
#define AlignGeneModel Wise2_AlignGeneModel
#define DYNAMITE_DEFINED_AlignGeneModel
#endif


struct Wise2_AlignGeneModelScore {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int len;     
    Score * forward_coding;  
    Score * reverse_coding;  
    Score * splice5_forward;     
    Score * splice3_forward;     
    Score * splice5_reverse;     
    Score * splice3_reverse;     
    SeqAlign * align;    
    Sequence * anchor;   
    } ;  
/* AlignGeneModelScore defined */ 
#ifndef DYNAMITE_DEFINED_AlignGeneModelScore
typedef struct Wise2_AlignGeneModelScore Wise2_AlignGeneModelScore;
#define AlignGeneModelScore Wise2_AlignGeneModelScore
#define DYNAMITE_DEFINED_AlignGeneModelScore
#endif


struct Wise2_AlignGeneModelParam {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    CompProb      * protein;     
    RandomModel   * rm;  
    CodonMapper   * cm;  
    CodonTable    * ct;  
    DnaProbMatrix * dm;  
    SpliceSiteProb  * ss5;   
    SpliceSiteProb  * ss3;   
    GeneStats * gs;  
    Probability   total_weight;  
    boolean       tolerate_nonanchor_stops;  
    Probability nonanchor_stop;  
    NonCodingSimpleModel * ncsm;     
    double coding_window_thres;  
    double noncoding_window_thres;   
    double coding_window_bonus;  
    double noncoding_window_pen;     
    } ;  
/* AlignGeneModelParam defined */ 
#ifndef DYNAMITE_DEFINED_AlignGeneModelParam
typedef struct Wise2_AlignGeneModelParam Wise2_AlignGeneModelParam;
#define AlignGeneModelParam Wise2_AlignGeneModelParam
#define DYNAMITE_DEFINED_AlignGeneModelParam
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_NonCodingSimpleModel(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [NonCodingSimpleModel *]
 *
 * Return [UNKN ]  Undocumented return value [NonCodingSimpleModel *]
 *
 */
NonCodingSimpleModel * Wise2_hard_link_NonCodingSimpleModel(NonCodingSimpleModel * obj);
#define hard_link_NonCodingSimpleModel Wise2_hard_link_NonCodingSimpleModel


/* Function:  NonCodingSimpleModel_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [NonCodingSimpleModel *]
 *
 */
NonCodingSimpleModel * Wise2_NonCodingSimpleModel_alloc(void);
#define NonCodingSimpleModel_alloc Wise2_NonCodingSimpleModel_alloc


/* Function:  free_NonCodingSimpleModel(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [NonCodingSimpleModel *]
 *
 * Return [UNKN ]  Undocumented return value [NonCodingSimpleModel *]
 *
 */
NonCodingSimpleModel * Wise2_free_NonCodingSimpleModel(NonCodingSimpleModel * obj);
#define free_NonCodingSimpleModel Wise2_free_NonCodingSimpleModel


/* Function:  hard_link_SpliceSiteProb(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SpliceSiteProb *]
 *
 * Return [UNKN ]  Undocumented return value [SpliceSiteProb *]
 *
 */
SpliceSiteProb * Wise2_hard_link_SpliceSiteProb(SpliceSiteProb * obj);
#define hard_link_SpliceSiteProb Wise2_hard_link_SpliceSiteProb


/* Function:  SpliceSiteProb_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SpliceSiteProb *]
 *
 */
SpliceSiteProb * Wise2_SpliceSiteProb_alloc(void);
#define SpliceSiteProb_alloc Wise2_SpliceSiteProb_alloc


/* Function:  free_SpliceSiteProb(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SpliceSiteProb *]
 *
 * Return [UNKN ]  Undocumented return value [SpliceSiteProb *]
 *
 */
SpliceSiteProb * Wise2_free_SpliceSiteProb(SpliceSiteProb * obj);
#define free_SpliceSiteProb Wise2_free_SpliceSiteProb


/* Function:  hard_link_AlignGeneColumnStore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AlignGeneColumnStore *]
 *
 * Return [UNKN ]  Undocumented return value [AlignGeneColumnStore *]
 *
 */
AlignGeneColumnStore * Wise2_hard_link_AlignGeneColumnStore(AlignGeneColumnStore * obj);
#define hard_link_AlignGeneColumnStore Wise2_hard_link_AlignGeneColumnStore


/* Function:  AlignGeneColumnStore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AlignGeneColumnStore *]
 *
 */
AlignGeneColumnStore * Wise2_AlignGeneColumnStore_alloc(void);
#define AlignGeneColumnStore_alloc Wise2_AlignGeneColumnStore_alloc


/* Function:  free_AlignGeneColumnStore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AlignGeneColumnStore *]
 *
 * Return [UNKN ]  Undocumented return value [AlignGeneColumnStore *]
 *
 */
AlignGeneColumnStore * Wise2_free_AlignGeneColumnStore(AlignGeneColumnStore * obj);
#define free_AlignGeneColumnStore Wise2_free_AlignGeneColumnStore


/* Function:  hard_link_AlignGeneModel(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AlignGeneModel *]
 *
 * Return [UNKN ]  Undocumented return value [AlignGeneModel *]
 *
 */
AlignGeneModel * Wise2_hard_link_AlignGeneModel(AlignGeneModel * obj);
#define hard_link_AlignGeneModel Wise2_hard_link_AlignGeneModel


/* Function:  AlignGeneModel_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AlignGeneModel *]
 *
 */
AlignGeneModel * Wise2_AlignGeneModel_alloc(void);
#define AlignGeneModel_alloc Wise2_AlignGeneModel_alloc


/* Function:  free_AlignGeneModel(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AlignGeneModel *]
 *
 * Return [UNKN ]  Undocumented return value [AlignGeneModel *]
 *
 */
AlignGeneModel * Wise2_free_AlignGeneModel(AlignGeneModel * obj);
#define free_AlignGeneModel Wise2_free_AlignGeneModel


/* Function:  hard_link_AlignGeneModelScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AlignGeneModelScore *]
 *
 * Return [UNKN ]  Undocumented return value [AlignGeneModelScore *]
 *
 */
AlignGeneModelScore * Wise2_hard_link_AlignGeneModelScore(AlignGeneModelScore * obj);
#define hard_link_AlignGeneModelScore Wise2_hard_link_AlignGeneModelScore


/* Function:  AlignGeneModelScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AlignGeneModelScore *]
 *
 */
AlignGeneModelScore * Wise2_AlignGeneModelScore_alloc(void);
#define AlignGeneModelScore_alloc Wise2_AlignGeneModelScore_alloc


/* Function:  free_AlignGeneModelScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AlignGeneModelScore *]
 *
 * Return [UNKN ]  Undocumented return value [AlignGeneModelScore *]
 *
 */
AlignGeneModelScore * Wise2_free_AlignGeneModelScore(AlignGeneModelScore * obj);
#define free_AlignGeneModelScore Wise2_free_AlignGeneModelScore


/* Function:  hard_link_AlignGeneModelParam(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AlignGeneModelParam *]
 *
 * Return [UNKN ]  Undocumented return value [AlignGeneModelParam *]
 *
 */
AlignGeneModelParam * Wise2_hard_link_AlignGeneModelParam(AlignGeneModelParam * obj);
#define hard_link_AlignGeneModelParam Wise2_hard_link_AlignGeneModelParam


/* Function:  AlignGeneModelParam_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AlignGeneModelParam *]
 *
 */
AlignGeneModelParam * Wise2_AlignGeneModelParam_alloc(void);
#define AlignGeneModelParam_alloc Wise2_AlignGeneModelParam_alloc


/* Function:  free_AlignGeneModelParam(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AlignGeneModelParam *]
 *
 * Return [UNKN ]  Undocumented return value [AlignGeneModelParam *]
 *
 */
AlignGeneModelParam * Wise2_free_AlignGeneModelParam(AlignGeneModelParam * obj);
#define free_AlignGeneModelParam Wise2_free_AlignGeneModelParam


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
Probability Wise2_weighted_coding_AlignGeneColumn(AlignGeneColumnStore * store,AlignGeneModelParam * agmp);
#define weighted_coding_AlignGeneColumn Wise2_weighted_coding_AlignGeneColumn
void Wise2_window_AlignGeneModel(SeqAlign * sal,AlignGeneModel * agm,AlignGeneModelParam * agmp);
#define window_AlignGeneModel Wise2_window_AlignGeneModel
NonCodingSimpleModel * Wise2_create_NonCodingSimpleModel(Probability change);
#define create_NonCodingSimpleModel Wise2_create_NonCodingSimpleModel
Probability Wise2_simple_non_coding_AlignGeneColumn(AlignGeneColumnStore * store,AlignGeneModelParam * agmp);
#define simple_non_coding_AlignGeneColumn Wise2_simple_non_coding_AlignGeneColumn
Probability Wise2_weighted_non_coding_AlignGeneColumn(AlignGeneColumnStore * store,AlignGeneModelParam * agmp);
#define weighted_non_coding_AlignGeneColumn Wise2_weighted_non_coding_AlignGeneColumn
Probability Wise2_weighted_simple_non_coding_AlignGeneColumn(AlignGeneColumnStore * store,AlignGeneModelParam * agmp);
#define weighted_simple_non_coding_AlignGeneColumn Wise2_weighted_simple_non_coding_AlignGeneColumn
Probability Wise2_simple_coding_AlignGeneColumn(AlignGeneColumnStore * store,AlignGeneModelParam * agmp);
#define simple_coding_AlignGeneColumn Wise2_simple_coding_AlignGeneColumn
boolean Wise2_fill_reverse_AlignGeneColumn(AlignGeneColumnStore * store,SeqAlign * sa,int forward_codon_end);
#define fill_reverse_AlignGeneColumn Wise2_fill_reverse_AlignGeneColumn
boolean Wise2_fill_forward_AlignGeneColumn(AlignGeneColumnStore * store,SeqAlign * sa,int codon_end);
#define fill_forward_AlignGeneColumn Wise2_fill_forward_AlignGeneColumn
AlignGeneColumnStore * Wise2_new_empty_AlignGeneColumnStore(int len);
#define new_empty_AlignGeneColumnStore Wise2_new_empty_AlignGeneColumnStore
AlignGeneCodon * Wise2_free_AlignGeneCodon(AlignGeneCodon * c);
#define free_AlignGeneCodon Wise2_free_AlignGeneCodon
AlignGeneModel * Wise2_new_AlignGeneModel(int len);
#define new_AlignGeneModel Wise2_new_AlignGeneModel
void Wise2_show_AlignGeneModel(AlignGeneModel * agm,SeqAlign * sal,CodonTable * ct,GenomicRegion * gr,FILE * ofp,AlignGeneModelParam * agmp);
#define show_AlignGeneModel Wise2_show_AlignGeneModel
AlignGeneModel * Wise2_create_AlignGeneModel(SeqAlign * sal,AlignGeneModelParam * agmp);
#define create_AlignGeneModel Wise2_create_AlignGeneModel
Probability Wise2_prob_SpliceSiteProb(SpliceSiteProb * ssp,Sequence * seq,int pos);
#define prob_SpliceSiteProb Wise2_prob_SpliceSiteProb
Probability Wise2_non_coding_probability_AlignGeneModel(SeqAlign * sal,int codon_end_pos,AlignGeneModelParam * agmp);
#define non_coding_probability_AlignGeneModel Wise2_non_coding_probability_AlignGeneModel
Probability Wise2_coding_probability_AlignGeneModel(SeqAlign * sal,int codon_end_pos,AlignGeneModelParam * agmp);
#define coding_probability_AlignGeneModel Wise2_coding_probability_AlignGeneModel
AlignGeneModelParam * Wise2_new_AlignGeneModelParam_from_argv(int * argc,char ** argv);
#define new_AlignGeneModelParam_from_argv Wise2_new_AlignGeneModelParam_from_argv
void Wise2_show_help_AlignGeneModelParam(FILE * ofp);
#define show_help_AlignGeneModelParam Wise2_show_help_AlignGeneModelParam
AlignGeneModelParam * Wise2_std_AlignGeneModelParam(CompProb * cp,DnaProbMatrix * dm,CodonTable * ct,GeneStats * gs,Probability change);
#define std_AlignGeneModelParam Wise2_std_AlignGeneModelParam
AlignGeneModelScore * Wise2_AlignGeneModelScore_from_AlignGeneModel(AlignGeneModel * agm);
#define AlignGeneModelScore_from_AlignGeneModel Wise2_AlignGeneModelScore_from_AlignGeneModel
Probability * Wise2_free_Probability(Probability * p);
#define free_Probability Wise2_free_Probability
Score * Wise2_free_Score(Score * p);
#define free_Score Wise2_free_Score


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/

#ifdef _cplusplus
}
#endif

#endif
