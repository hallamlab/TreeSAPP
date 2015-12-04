#ifndef DYNAMITEgenestatsHEADERFILE
#define DYNAMITEgenestatsHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "wisebase.h"
#include "seqalign.h"
#include "pwmdna.h"
#include "randommodel.h"
#include "complexsequence.h"
#include "complexevalset.h" /* for the standard evals */


#define DEFAULT_SPLICE_OFFSET_SCORE 1.5


struct Wise2_SpliceSiteScore {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    pwmDNAScore * score;     
    int offset;  
    Score min_collar;    
    Score max_collar;    
    Score score_offset;  
    } ;  
/* SpliceSiteScore defined */ 
#ifndef DYNAMITE_DEFINED_SpliceSiteScore
typedef struct Wise2_SpliceSiteScore Wise2_SpliceSiteScore;
#define SpliceSiteScore Wise2_SpliceSiteScore
#define DYNAMITE_DEFINED_SpliceSiteScore
#endif


/* Object GeneStats
 *
 * Descrip: This structure is to hold the
 *        new generation of gene statistics
 *        for the genewise algorithms
 *
 *
 */
struct Wise2_GeneStats {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    SeqAlign * splice5;  
    int splice5_offset;  
    SeqAlign * splice3;  
    int splice3_offset;  
    RandomModelDNA * intron;    /*  actually counts */ 
    double average_intron;   
    RandomModelDNA * polyp; /*  actually counts */ 
    double average_polyp;    
    RandomModelDNA * rnd;    
    double codon[64];    
    } ;  
/* GeneStats defined */ 
#ifndef DYNAMITE_DEFINED_GeneStats
typedef struct Wise2_GeneStats Wise2_GeneStats;
#define GeneStats Wise2_GeneStats
#define DYNAMITE_DEFINED_GeneStats
#endif


/* Object GeneModelParam
 *
 * Descrip: A small helper object containing
 *        the ways of converting the actual
 *        counts/alignments to the model
 *
 *        Here really for convience so we can
 *        keep all the code associated with the
 *        creation of the GeneModel together
 *
 *
 */
struct Wise2_GeneModelParam {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    double splice5_pseudo;   
    double splice3_pseudo;   
    double intron_emission_pseudo;   
    double polyp_emission_pseudo;    
    Bits min_collar;     
    Bits max_collar;     
    Bits score_offset;   
    char * gene_stats_file;  
    boolean use_gtag_splice;     
    double prob_for_gtag;    
    } ;  
/* GeneModelParam defined */ 
#ifndef DYNAMITE_DEFINED_GeneModelParam
typedef struct Wise2_GeneModelParam Wise2_GeneModelParam;
#define GeneModelParam Wise2_GeneModelParam
#define DYNAMITE_DEFINED_GeneModelParam
#endif


/* Object GeneModel
 *
 * Descrip: This structure is to hold the
 *        new generation of models for the
 *        genewise algorithm
 *
 *
 *
 */
struct Wise2_GeneModel {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    pwmDNA * splice5;    
    int splice5_offset;  
    pwmDNA * splice3;    
    int splice3_offset;  
    RandomModelDNA * intron;     
    double intron_stay_prob;     
    RandomModelDNA * polyp;  
    double polyp_stay_prob;  
    RandomModelDNA * rnd;    
    SpliceSiteScore * splice5score;  
    SpliceSiteScore * splice3score;  
    double codon[64];    
    boolean use_gtag_splice;     
    Score score_for_gtag;    
    } ;  
/* GeneModel defined */ 
#ifndef DYNAMITE_DEFINED_GeneModel
typedef struct Wise2_GeneModel Wise2_GeneModel;
#define GeneModel Wise2_GeneModel
#define DYNAMITE_DEFINED_GeneModel
#endif


struct Wise2_GeneralGeneModelScore {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    RandomCodonScore * start;    
    RandomCodonScore * stop;     
    RandomCodonScore * general;  
    } ;  
/* GeneralGeneModelScore defined */ 
#ifndef DYNAMITE_DEFINED_GeneralGeneModelScore
typedef struct Wise2_GeneralGeneModelScore Wise2_GeneralGeneModelScore;
#define GeneralGeneModelScore Wise2_GeneralGeneModelScore
#define DYNAMITE_DEFINED_GeneralGeneModelScore
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  show_help_GeneModelParam(ofp)
 *
 * Descrip:    Shows help
 *
 *
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_show_help_GeneModelParam(FILE * ofp);
#define show_help_GeneModelParam Wise2_show_help_GeneModelParam


/* Function:  show_info_GeneModelParam(p,ofp)
 *
 * Descrip:    Shows genemodel param for info
 *
 *
 * Arg:          p [UNKN ] Undocumented argument [GeneModelParam *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_show_info_GeneModelParam(GeneModelParam * p,FILE * ofp);
#define show_info_GeneModelParam Wise2_show_info_GeneModelParam


/* Function:  new_GeneModelParam_from_argv(argc,argv)
 *
 * Descrip:    Makes a GeneModelParam from argv
 *
 *
 * Arg:        argc [UNKN ] Undocumented argument [int *]
 * Arg:        argv [UNKN ] Undocumented argument [char **]
 *
 * Return [UNKN ]  Undocumented return value [GeneModelParam *]
 *
 */
GeneModelParam * Wise2_new_GeneModelParam_from_argv(int * argc,char ** argv);
#define new_GeneModelParam_from_argv Wise2_new_GeneModelParam_from_argv


/* Function:  vanilla_GeneralGeneModelScore(ct,start_odds,general_odds,stop_odds)
 *
 * Descrip:    Makes a vanilla general gene model score
 *
 *
 * Arg:                  ct [UNKN ] Undocumented argument [CodonTable *]
 * Arg:          start_odds [UNKN ] Undocumented argument [Probability]
 * Arg:        general_odds [UNKN ] Undocumented argument [Probability]
 * Arg:           stop_odds [UNKN ] Undocumented argument [Probability]
 *
 * Return [UNKN ]  Undocumented return value [GeneralGeneModelScore *]
 *
 */
GeneralGeneModelScore * Wise2_vanilla_GeneralGeneModelScore(CodonTable * ct,Probability start_odds,Probability general_odds,Probability stop_odds);
#define vanilla_GeneralGeneModelScore Wise2_vanilla_GeneralGeneModelScore


/* Function:  vanilla_start_RandomCodon(ct,start_codon_odd_prob)
 *
 * Descrip:    Makes a vanilla (ie, all things equal) Met based
 *             start codon rndscore
 *
 *
 * Arg:                          ct [UNKN ] Undocumented argument [CodonTable *]
 * Arg:        start_codon_odd_prob [UNKN ] Undocumented argument [Probability]
 *
 * Return [UNKN ]  Undocumented return value [RandomCodon *]
 *
 */
RandomCodon * Wise2_vanilla_start_RandomCodon(CodonTable * ct,Probability start_codon_odd_prob);
#define vanilla_start_RandomCodon Wise2_vanilla_start_RandomCodon


/* Function:  vanilla_stop_RandomCodon(ct,stop_codon_odd_prob)
 *
 * Descrip:    Makes a vanilla (ie, all things equal) * based
 *             stop codon rndscore
 *
 *
 * Arg:                         ct [UNKN ] Undocumented argument [CodonTable *]
 * Arg:        stop_codon_odd_prob [UNKN ] Undocumented argument [Probability]
 *
 * Return [UNKN ]  Undocumented return value [RandomCodon *]
 *
 */
RandomCodon * Wise2_vanilla_stop_RandomCodon(CodonTable * ct,Probability stop_codon_odd_prob);
#define vanilla_stop_RandomCodon Wise2_vanilla_stop_RandomCodon


/* Function:  vanilla_general_RandomCodon(ct,coding_odd_prob)
 *
 * Descrip:    Makes a vanilla (ie, all things equal) general non-stop
 *
 *
 * Arg:                     ct [UNKN ] Undocumented argument [CodonTable *]
 * Arg:        coding_odd_prob [UNKN ] Undocumented argument [Probability]
 *
 * Return [UNKN ]  Undocumented return value [RandomCodon *]
 *
 */
RandomCodon * Wise2_vanilla_general_RandomCodon(CodonTable * ct,Probability coding_odd_prob);
#define vanilla_general_RandomCodon Wise2_vanilla_general_RandomCodon


/* Function:  std_GeneModelParam(void)
 *
 * Descrip:    Makes a standard GeneModelParam
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneModelParam *]
 *
 */
GeneModelParam * Wise2_std_GeneModelParam(void);
#define std_GeneModelParam Wise2_std_GeneModelParam


/* Function:  GeneModel_from_GeneModelParam(p)
 *
 * Descrip:    Combines GeneStats_from_GeneModelParam and GeneModel_from_GeneStats
 *
 *
 * Arg:        p [UNKN ] Undocumented argument [GeneModelParam *]
 *
 * Return [UNKN ]  Undocumented return value [GeneModel *]
 *
 */
GeneModel * Wise2_GeneModel_from_GeneModelParam(GeneModelParam * p);
#define GeneModel_from_GeneModelParam Wise2_GeneModel_from_GeneModelParam


/* Function:  GeneStats_from_GeneModelParam(p)
 *
 * Descrip:    Makes a GeneStats from GeneModelParam - basically just opening the file
 *
 *
 * Arg:        p [UNKN ] Undocumented argument [GeneModelParam *]
 *
 * Return [UNKN ]  Undocumented return value [GeneStats *]
 *
 */
GeneStats * Wise2_GeneStats_from_GeneModelParam(GeneModelParam * p);
#define GeneStats_from_GeneModelParam Wise2_GeneStats_from_GeneModelParam


/* Function:  GeneModel_from_GeneStats(gs,p)
 *
 * Descrip:    Makes a model from the stats file
 *
 *
 * Arg:        gs [UNKN ] Undocumented argument [GeneStats *]
 * Arg:         p [UNKN ] Undocumented argument [GeneModelParam *]
 *
 * Return [UNKN ]  Undocumented return value [GeneModel *]
 *
 */
GeneModel * Wise2_GeneModel_from_GeneStats(GeneStats * gs,GeneModelParam * p);
#define GeneModel_from_GeneStats Wise2_GeneModel_from_GeneStats


/* Function:  show_GeneModel(gm,ofp)
 *
 * Descrip:    shows a genemodel
 *
 *
 * Arg:         gm [UNKN ] Undocumented argument [GeneModel *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_show_GeneModel(GeneModel * gm,FILE * ofp);
#define show_GeneModel Wise2_show_GeneModel


/* Function:  new_ComplexSequenceEvalSet_from_GeneModel(gm)
 *
 * Descrip:    Makes an entire ComplexSequenceEvalSet for genomic work
 *
 *
 * Arg:        gm [UNKN ] Undocumented argument [GeneModel *]
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequenceEvalSet *]
 *
 */
ComplexSequenceEvalSet * Wise2_new_ComplexSequenceEvalSet_from_GeneModel(GeneModel * gm);
#define new_ComplexSequenceEvalSet_from_GeneModel Wise2_new_ComplexSequenceEvalSet_from_GeneModel


/* Function:  read_GeneStats(ifp)
 *
 * Descrip:    Reads a GeneStats file
 *
 *
 * Arg:        ifp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [GeneStats *]
 *
 */
GeneStats * Wise2_read_GeneStats(FILE * ifp);
#define read_GeneStats Wise2_read_GeneStats


/* Function:  dump_GeneStats(st,ofp)
 *
 * Descrip:    testing function
 *
 *
 * Arg:         st [UNKN ] Undocumented argument [GeneStats *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_dump_GeneStats(GeneStats * st,FILE * ofp);
#define dump_GeneStats Wise2_dump_GeneStats


/* Function:  read_codon_GeneStats(codon_array,line,ifp)
 *
 * Descrip:    assummes codon_array is 64 positions long
 *               
 *             line should have begin consensus on it and be of MAXLINE length as it will be used as the buffer.
 *
 *             This does **not** check that you have filled up all 64 positions.
 *
 *
 * Arg:        codon_array [UNKN ] Undocumented argument [double *]
 * Arg:               line [UNKN ] Undocumented argument [char*]
 * Arg:                ifp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_read_codon_GeneStats(double * codon_array,char* line,FILE * ifp);
#define read_codon_GeneStats Wise2_read_codon_GeneStats


/* Function:  hard_link_SpliceSiteScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SpliceSiteScore *]
 *
 * Return [UNKN ]  Undocumented return value [SpliceSiteScore *]
 *
 */
SpliceSiteScore * Wise2_hard_link_SpliceSiteScore(SpliceSiteScore * obj);
#define hard_link_SpliceSiteScore Wise2_hard_link_SpliceSiteScore


/* Function:  SpliceSiteScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SpliceSiteScore *]
 *
 */
SpliceSiteScore * Wise2_SpliceSiteScore_alloc(void);
#define SpliceSiteScore_alloc Wise2_SpliceSiteScore_alloc


/* Function:  free_SpliceSiteScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SpliceSiteScore *]
 *
 * Return [UNKN ]  Undocumented return value [SpliceSiteScore *]
 *
 */
SpliceSiteScore * Wise2_free_SpliceSiteScore(SpliceSiteScore * obj);
#define free_SpliceSiteScore Wise2_free_SpliceSiteScore


/* Function:  hard_link_GeneStats(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GeneStats *]
 *
 * Return [UNKN ]  Undocumented return value [GeneStats *]
 *
 */
GeneStats * Wise2_hard_link_GeneStats(GeneStats * obj);
#define hard_link_GeneStats Wise2_hard_link_GeneStats


/* Function:  GeneStats_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneStats *]
 *
 */
GeneStats * Wise2_GeneStats_alloc(void);
#define GeneStats_alloc Wise2_GeneStats_alloc


/* Function:  free_GeneStats(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneStats *]
 *
 * Return [UNKN ]  Undocumented return value [GeneStats *]
 *
 */
GeneStats * Wise2_free_GeneStats(GeneStats * obj);
#define free_GeneStats Wise2_free_GeneStats


/* Function:  hard_link_GeneModelParam(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GeneModelParam *]
 *
 * Return [UNKN ]  Undocumented return value [GeneModelParam *]
 *
 */
GeneModelParam * Wise2_hard_link_GeneModelParam(GeneModelParam * obj);
#define hard_link_GeneModelParam Wise2_hard_link_GeneModelParam


/* Function:  GeneModelParam_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneModelParam *]
 *
 */
GeneModelParam * Wise2_GeneModelParam_alloc(void);
#define GeneModelParam_alloc Wise2_GeneModelParam_alloc


/* Function:  free_GeneModelParam(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneModelParam *]
 *
 * Return [UNKN ]  Undocumented return value [GeneModelParam *]
 *
 */
GeneModelParam * Wise2_free_GeneModelParam(GeneModelParam * obj);
#define free_GeneModelParam Wise2_free_GeneModelParam


/* Function:  hard_link_GeneModel(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GeneModel *]
 *
 * Return [UNKN ]  Undocumented return value [GeneModel *]
 *
 */
GeneModel * Wise2_hard_link_GeneModel(GeneModel * obj);
#define hard_link_GeneModel Wise2_hard_link_GeneModel


/* Function:  GeneModel_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneModel *]
 *
 */
GeneModel * Wise2_GeneModel_alloc(void);
#define GeneModel_alloc Wise2_GeneModel_alloc


/* Function:  free_GeneModel(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneModel *]
 *
 * Return [UNKN ]  Undocumented return value [GeneModel *]
 *
 */
GeneModel * Wise2_free_GeneModel(GeneModel * obj);
#define free_GeneModel Wise2_free_GeneModel


/* Function:  hard_link_GeneralGeneModelScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GeneralGeneModelScore *]
 *
 * Return [UNKN ]  Undocumented return value [GeneralGeneModelScore *]
 *
 */
GeneralGeneModelScore * Wise2_hard_link_GeneralGeneModelScore(GeneralGeneModelScore * obj);
#define hard_link_GeneralGeneModelScore Wise2_hard_link_GeneralGeneModelScore


/* Function:  GeneralGeneModelScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneralGeneModelScore *]
 *
 */
GeneralGeneModelScore * Wise2_GeneralGeneModelScore_alloc(void);
#define GeneralGeneModelScore_alloc Wise2_GeneralGeneModelScore_alloc


/* Function:  free_GeneralGeneModelScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneralGeneModelScore *]
 *
 * Return [UNKN ]  Undocumented return value [GeneralGeneModelScore *]
 *
 */
GeneralGeneModelScore * Wise2_free_GeneralGeneModelScore(GeneralGeneModelScore * obj);
#define free_GeneralGeneModelScore Wise2_free_GeneralGeneModelScore


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
ComplexSequenceEval * Wise2_ComplexSequenceEval_for_scored_ag(Score * score_for_ag);
#define ComplexSequenceEval_for_scored_ag Wise2_ComplexSequenceEval_for_scored_ag
int Wise2_scored_ag_eval_func(int type,void *data,char * seq);
#define scored_ag_eval_func Wise2_scored_ag_eval_func
ComplexSequenceEval * Wise2_ComplexSequenceEval_for_scored_gt(Score * score_for_gt);
#define ComplexSequenceEval_for_scored_gt Wise2_ComplexSequenceEval_for_scored_gt
int Wise2_scored_gt_eval_func(int type,void *data,char * seq);
#define scored_gt_eval_func Wise2_scored_gt_eval_func
ComplexSequenceEval * Wise2_ComplexSequenceEval_from_pwmDNAScore_splice(SpliceSiteScore * score);
#define ComplexSequenceEval_from_pwmDNAScore_splice Wise2_ComplexSequenceEval_from_pwmDNAScore_splice
int Wise2_pwmDNA_splice_ComplexSequence_eval_func(int type,void * data,char * seq);
#define pwmDNA_splice_ComplexSequence_eval_func Wise2_pwmDNA_splice_ComplexSequence_eval_func
RandomModelDNA * Wise2_get_genestat_emission(char * buffer);
#define get_genestat_emission Wise2_get_genestat_emission

#ifdef _cplusplus
}
#endif

#endif
