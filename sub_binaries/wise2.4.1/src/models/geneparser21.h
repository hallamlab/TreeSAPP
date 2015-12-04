#ifndef DYNAMITEgeneparser21HEADERFILE
#define DYNAMITEgeneparser21HEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "probability.h"
#include "dyna.h"
#include "randommodel.h"
#include "genefrequency.h"

enum {
  GP21_CDS2CENTRAL = 0,
  GP21_CENTRAL2CENTRAL,
  GP21_CENTRAL2PY,
  GP21_PY2PY,
  GP21_PY2SPACER,
  GP21_PY2CDS,
  GP21_SPACER2SPACER,
  GP21_SPACER2CDS,
  GP21_INSERT_1_BASE,
  GP21_INSERT_2_BASE,
  GP21_DELETE_1_BASE,
  GP21_DELETE_2_BASE,
  GP21_CDS2CDS, /* this is for gene parsing outside of the *Wise set */
  GP21_RND2CDS,
  GP21_CDS2RND,
  GP21_RND2RND,
  GP21_RND2MODEL,
  GP21_LINK2MODEL, /* this is for link models */
  GP21_LINK2LINK,
  GP21_LINK2RND,
  GENEPARSER21_TRANSITION_LEN
};


#define GENEPARSER21_EMISSION_LEN 5

#define GeneParser21SetLISTLENGTH 32


struct Wise2_GeneParser21 {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Probability transition[GENEPARSER21_TRANSITION_LEN];     
    Probability central[GENEPARSER21_EMISSION_LEN];  
    Probability py[GENEPARSER21_EMISSION_LEN];   
    Probability spacer[GENEPARSER21_EMISSION_LEN];   
    } ;  
/* GeneParser21 defined */ 
#ifndef DYNAMITE_DEFINED_GeneParser21
typedef struct Wise2_GeneParser21 Wise2_GeneParser21;
#define GeneParser21 Wise2_GeneParser21
#define DYNAMITE_DEFINED_GeneParser21
#endif


struct Wise2_GeneParser21Score {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Score transition[GENEPARSER21_TRANSITION_LEN];   
    Score central[GENEPARSER21_EMISSION_LEN];    
    Score py[GENEPARSER21_EMISSION_LEN];     
    Score spacer[GENEPARSER21_EMISSION_LEN];     
    } ;  
/* GeneParser21Score defined */ 
#ifndef DYNAMITE_DEFINED_GeneParser21Score
typedef struct Wise2_GeneParser21Score Wise2_GeneParser21Score;
#define GeneParser21Score Wise2_GeneParser21Score
#define DYNAMITE_DEFINED_GeneParser21Score
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_GeneParser21(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GeneParser21 *]
 *
 * Return [UNKN ]  Undocumented return value [GeneParser21 *]
 *
 */
GeneParser21 * Wise2_hard_link_GeneParser21(GeneParser21 * obj);
#define hard_link_GeneParser21 Wise2_hard_link_GeneParser21


/* Function:  GeneParser21_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneParser21 *]
 *
 */
GeneParser21 * Wise2_GeneParser21_alloc(void);
#define GeneParser21_alloc Wise2_GeneParser21_alloc


/* Function:  free_GeneParser21(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneParser21 *]
 *
 * Return [UNKN ]  Undocumented return value [GeneParser21 *]
 *
 */
GeneParser21 * Wise2_free_GeneParser21(GeneParser21 * obj);
#define free_GeneParser21 Wise2_free_GeneParser21


/* Function:  hard_link_GeneParser21Score(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GeneParser21Score *]
 *
 * Return [UNKN ]  Undocumented return value [GeneParser21Score *]
 *
 */
GeneParser21Score * Wise2_hard_link_GeneParser21Score(GeneParser21Score * obj);
#define hard_link_GeneParser21Score Wise2_hard_link_GeneParser21Score


/* Function:  GeneParser21Score_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneParser21Score *]
 *
 */
GeneParser21Score * Wise2_GeneParser21Score_alloc(void);
#define GeneParser21Score_alloc Wise2_GeneParser21Score_alloc


/* Function:  free_GeneParser21Score(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneParser21Score *]
 *
 * Return [UNKN ]  Undocumented return value [GeneParser21Score *]
 *
 */
GeneParser21Score * Wise2_free_GeneParser21Score(GeneParser21Score * obj);
#define free_GeneParser21Score Wise2_free_GeneParser21Score


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
RandomModelDNA * Wise2_fudged_mixed_RandomModelDNA_from_GeneParser21(GeneParser21 * gp,RandomModelDNA * rnd);
#define fudged_mixed_RandomModelDNA_from_GeneParser21 Wise2_fudged_mixed_RandomModelDNA_from_GeneParser21
void Wise2_add_flat_error_probabilities_GeneParser21(GeneParser21 * gp21,Probability error);
#define add_flat_error_probabilities_GeneParser21 Wise2_add_flat_error_probabilities_GeneParser21
void Wise2_add_error_probabilities_GeneParser21(GeneParser21 * gp21,Probability insert_1,Probability insert_2,Probability delete_1,Probability delete_2);
#define add_error_probabilities_GeneParser21 Wise2_add_error_probabilities_GeneParser21
GeneParser21 * Wise2_GeneParser21_from_GeneFrequency21_cds(GeneFrequency21 * gf,Probability rnd_loop,Probability cds_loop,Probability rnd_to_model,Probability link_loop,Probability link_to_model);
#define GeneParser21_from_GeneFrequency21_cds Wise2_GeneParser21_from_GeneFrequency21_cds
RandomModelDNA * Wise2_RandomModelDNA_from_central_GeneParser21(GeneParser21 *gp21);
#define RandomModelDNA_from_central_GeneParser21 Wise2_RandomModelDNA_from_central_GeneParser21
void Wise2_show_GeneParser21(GeneParser21 * gp21,FILE * ofp);
#define show_GeneParser21 Wise2_show_GeneParser21
GeneParser21 * Wise2_std_GeneParser21(void);
#define std_GeneParser21 Wise2_std_GeneParser21
Probability Wise2_removed_probability_from_cds(GeneParser21 * gp21);
#define removed_probability_from_cds Wise2_removed_probability_from_cds
void Wise2_GeneParser21_fold_in_RandomModelDNA(GeneParser21 * gp21,RandomModelDNA * rmd);
#define GeneParser21_fold_in_RandomModelDNA Wise2_GeneParser21_fold_in_RandomModelDNA
GeneParser21Score * Wise2_GeneParser21Score_from_GeneParser21(GeneParser21 * gp21);
#define GeneParser21Score_from_GeneParser21 Wise2_GeneParser21Score_from_GeneParser21


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/

#ifdef _cplusplus
}
#endif

#endif
