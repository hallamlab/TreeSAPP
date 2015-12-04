#ifndef DYNAMITEgeneparameterHEADERFILE
#define DYNAMITEgeneparameterHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "geneparser21.h"
#include "genefrequency.h"
#include "splicesitemodeler.h"
#include "complexevalset.h"
#include "probability.h"

#include "genestats.h"
#define GeneParameter21LISTLENGTH 32



struct Wise2_GeneWiseCodonModel {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    double in_donor[64];     
    double in_acceptor[64];  
    double in_cds[64];   
    } ;  
/* GeneWiseCodonModel defined */ 
#ifndef DYNAMITE_DEFINED_GeneWiseCodonModel
typedef struct Wise2_GeneWiseCodonModel Wise2_GeneWiseCodonModel;
#define GeneWiseCodonModel Wise2_GeneWiseCodonModel
#define DYNAMITE_DEFINED_GeneWiseCodonModel
#endif


/* Object GeneParameter21
 *
 * Descrip: No Description
 *
 */
struct Wise2_GeneParameter21 {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    GeneParser21 * gp;   
    CodonMapper  * cm;   
    ComplexSequenceEvalSet * cses;   
    SpliceSiteModel ** ss;  /*  held only to be free'd when GeneParser21Set is free'd */ 
    int len;/* len for above ss  */ 
    int maxlen; /* maxlen for above ss */ 
    RandomCodon  * rc;  /*  needed to soak up the odd-and-sods of genes */ 
    GeneWiseCodonModel * gwcm;   
    CodonTable   * ct;   
    boolean modelled_splice;    /*  so we can alter balance scores. */ 
    GeneModel    * gms;  
    } ;  
/* GeneParameter21 defined */ 
#ifndef DYNAMITE_DEFINED_GeneParameter21
typedef struct Wise2_GeneParameter21 Wise2_GeneParameter21;
#define GeneParameter21 Wise2_GeneParameter21
#define DYNAMITE_DEFINED_GeneParameter21
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  GeneParameter21_from_GeneModel(gm,ct,rnd_loop,cds_loop,rnd_to_model,link_loop,link_to_model,subs_error,indel_error)
 *
 * Descrip:    This actually makes the GeneParameter21 stuff from the
 *             new statistics
 *
 *
 * Arg:                   gm [UNKN ] Undocumented argument [GeneModel *]
 * Arg:                   ct [UNKN ] Undocumented argument [CodonTable *]
 * Arg:             rnd_loop [UNKN ] Undocumented argument [Probability]
 * Arg:             cds_loop [UNKN ] Undocumented argument [Probability]
 * Arg:         rnd_to_model [UNKN ] Undocumented argument [Probability]
 * Arg:            link_loop [UNKN ] Undocumented argument [Probability]
 * Arg:        link_to_model [UNKN ] Undocumented argument [Probability]
 * Arg:           subs_error [UNKN ] Undocumented argument [Probability]
 * Arg:          indel_error [UNKN ] Undocumented argument [Probability]
 *
 * Return [UNKN ]  Undocumented return value [GeneParameter21 *]
 *
 */
GeneParameter21 * Wise2_GeneParameter21_from_GeneModel(GeneModel * gm,CodonTable * ct,Probability rnd_loop,Probability cds_loop,Probability rnd_to_model,Probability link_loop,Probability link_to_model,Probability subs_error,Probability indel_error);
#define GeneParameter21_from_GeneModel Wise2_GeneParameter21_from_GeneModel


/* Function:  hard_link_GeneWiseCodonModel(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GeneWiseCodonModel *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseCodonModel *]
 *
 */
GeneWiseCodonModel * Wise2_hard_link_GeneWiseCodonModel(GeneWiseCodonModel * obj);
#define hard_link_GeneWiseCodonModel Wise2_hard_link_GeneWiseCodonModel


/* Function:  GeneWiseCodonModel_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseCodonModel *]
 *
 */
GeneWiseCodonModel * Wise2_GeneWiseCodonModel_alloc(void);
#define GeneWiseCodonModel_alloc Wise2_GeneWiseCodonModel_alloc


/* Function:  free_GeneWiseCodonModel(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneWiseCodonModel *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseCodonModel *]
 *
 */
GeneWiseCodonModel * Wise2_free_GeneWiseCodonModel(GeneWiseCodonModel * obj);
#define free_GeneWiseCodonModel Wise2_free_GeneWiseCodonModel


/* Function:  add_GeneParameter21(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GeneParameter21 *]
 * Arg:        add [OWNER] Object to add to the list [SpliceSiteModel *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_GeneParameter21(GeneParameter21 * obj,SpliceSiteModel * add);
#define add_GeneParameter21 Wise2_add_GeneParameter21


/* Function:  flush_GeneParameter21(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [GeneParameter21 *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_GeneParameter21(GeneParameter21 * obj);
#define flush_GeneParameter21 Wise2_flush_GeneParameter21


/* Function:  GeneParameter21_alloc_std(void)
 *
 * Descrip:    Equivalent to GeneParameter21_alloc_len(GeneParameter21LISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneParameter21 *]
 *
 */
GeneParameter21 * Wise2_GeneParameter21_alloc_std(void);
#define GeneParameter21_alloc_std Wise2_GeneParameter21_alloc_std


/* Function:  GeneParameter21_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [GeneParameter21 *]
 *
 */
GeneParameter21 * Wise2_GeneParameter21_alloc_len(int len);
#define GeneParameter21_alloc_len Wise2_GeneParameter21_alloc_len


/* Function:  hard_link_GeneParameter21(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GeneParameter21 *]
 *
 * Return [UNKN ]  Undocumented return value [GeneParameter21 *]
 *
 */
GeneParameter21 * Wise2_hard_link_GeneParameter21(GeneParameter21 * obj);
#define hard_link_GeneParameter21 Wise2_hard_link_GeneParameter21


/* Function:  GeneParameter21_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneParameter21 *]
 *
 */
GeneParameter21 * Wise2_GeneParameter21_alloc(void);
#define GeneParameter21_alloc Wise2_GeneParameter21_alloc


/* Function:  free_GeneParameter21(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneParameter21 *]
 *
 * Return [UNKN ]  Undocumented return value [GeneParameter21 *]
 *
 */
GeneParameter21 * Wise2_free_GeneParameter21(GeneParameter21 * obj);
#define free_GeneParameter21 Wise2_free_GeneParameter21


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
GeneWiseCodonModel * Wise2_GeneWiseCodonModel_from_GeneFrequencies(double * cds,GeneConsensus * donor,GeneConsensus * acceptor);
#define GeneWiseCodonModel_from_GeneFrequencies Wise2_GeneWiseCodonModel_from_GeneFrequencies
GeneParameter21 * Wise2_GeneParameter21_from_GeneFrequency21(GeneFrequency21 * gf,CodonTable * ct,RandomModelDNA * rmd,Probability rnd_loop,Probability cds_loop,Probability rnd_to_model,Probability link_loop,Probability link_to_model);
#define GeneParameter21_from_GeneFrequency21 Wise2_GeneParameter21_from_GeneFrequency21


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_GeneParameter21(SpliceSiteModel ** list,int i,int j) ;
#define swap_GeneParameter21 Wise2_swap_GeneParameter21
void Wise2_qsort_GeneParameter21(SpliceSiteModel ** list,int left,int right,int (*comp)(SpliceSiteModel * ,SpliceSiteModel * ));
#define qsort_GeneParameter21 Wise2_qsort_GeneParameter21
void Wise2_sort_GeneParameter21(GeneParameter21 * obj,int (*comp)(SpliceSiteModel *, SpliceSiteModel *));
#define sort_GeneParameter21 Wise2_sort_GeneParameter21
boolean Wise2_expand_GeneParameter21(GeneParameter21 * obj,int len);
#define expand_GeneParameter21 Wise2_expand_GeneParameter21

#ifdef _cplusplus
}
#endif

#endif
