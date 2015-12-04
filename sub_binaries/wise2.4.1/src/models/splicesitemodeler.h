#ifndef DYNAMITEsplicesitemodelerHEADERFILE
#define DYNAMITEsplicesitemodelerHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"
#include "randommodel.h"


#define SPLICESITEMODELER_TYPE 474


struct Wise2_SpliceSiteModel {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int offset;  
    int pre_splice_site;     
    int post_splice_site;    
    int start_random;    
    int stop_random;     
    ComplexConsensi     * cc;    
    RandomModelDNAScore * rmds;  
    Score error_pos;     
    } ;  
/* SpliceSiteModel defined */ 
#ifndef DYNAMITE_DEFINED_SpliceSiteModel
typedef struct Wise2_SpliceSiteModel Wise2_SpliceSiteModel;
#define SpliceSiteModel Wise2_SpliceSiteModel
#define DYNAMITE_DEFINED_SpliceSiteModel
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_SpliceSiteModel(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SpliceSiteModel *]
 *
 * Return [UNKN ]  Undocumented return value [SpliceSiteModel *]
 *
 */
SpliceSiteModel * Wise2_hard_link_SpliceSiteModel(SpliceSiteModel * obj);
#define hard_link_SpliceSiteModel Wise2_hard_link_SpliceSiteModel


/* Function:  SpliceSiteModel_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SpliceSiteModel *]
 *
 */
SpliceSiteModel * Wise2_SpliceSiteModel_alloc(void);
#define SpliceSiteModel_alloc Wise2_SpliceSiteModel_alloc


/* Function:  free_SpliceSiteModel(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SpliceSiteModel *]
 *
 * Return [UNKN ]  Undocumented return value [SpliceSiteModel *]
 *
 */
SpliceSiteModel * Wise2_free_SpliceSiteModel(SpliceSiteModel * obj);
#define free_SpliceSiteModel Wise2_free_SpliceSiteModel


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
int Wise2_SpliceSiteModeler_ComplexSequence_eval_func(int type,void * data,char * seq);
#define SpliceSiteModeler_ComplexSequence_eval_func Wise2_SpliceSiteModeler_ComplexSequence_eval_func
ComplexSequenceEval * Wise2_ComplexSequenceEval_from_SpliceSiteModel(SpliceSiteModel * ssm);
#define ComplexSequenceEval_from_SpliceSiteModel Wise2_ComplexSequenceEval_from_SpliceSiteModel
Score Wise2_SpliceSiteModel_score(SpliceSiteModel * ssm,char * seq);
#define SpliceSiteModel_score Wise2_SpliceSiteModel_score
SpliceSiteModel * Wise2_std_5SS_SpliceSiteModel(int offset,ComplexConsensi * cc,RandomModelDNAScore * rmds);
#define std_5SS_SpliceSiteModel Wise2_std_5SS_SpliceSiteModel
SpliceSiteModel * Wise2_std_3SS_SpliceSiteModel(int offset,ComplexConsensi * cc,RandomModelDNAScore * rmds);
#define std_3SS_SpliceSiteModel Wise2_std_3SS_SpliceSiteModel
SpliceSiteModel * Wise2_new_SpliceSiteModel(int offset,int pre_length,int post_length,int start,int stop,ComplexConsensi * cc,RandomModelDNAScore * rmds,Probability error);
#define new_SpliceSiteModel Wise2_new_SpliceSiteModel


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/

#ifdef _cplusplus
}
#endif

#endif
