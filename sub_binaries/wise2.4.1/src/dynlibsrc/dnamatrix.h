#ifndef DYNAMITEdnamatrixHEADERFILE
#define DYNAMITEdnamatrixHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"

typedef enum NMaskType {
  NMaskType_BASE,
  NMaskType_VARIABLE,
  NMaskType_EXCLUDED,
  NMaskType_BANNED
} NMaskType;

#define DnaMatrix_MATCH(mat,one,two) (mat->score[one][two])

struct Wise2_DnaProbMatrix {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Probability prob[5][5];  
    } ;  
/* DnaProbMatrix defined */ 
#ifndef DYNAMITE_DEFINED_DnaProbMatrix
typedef struct Wise2_DnaProbMatrix Wise2_DnaProbMatrix;
#define DnaProbMatrix Wise2_DnaProbMatrix
#define DYNAMITE_DEFINED_DnaProbMatrix
#endif


struct Wise2_DnaMatrix {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Score score[5][5];   
    } ;  
/* DnaMatrix defined */ 
#ifndef DYNAMITE_DEFINED_DnaMatrix
typedef struct Wise2_DnaMatrix Wise2_DnaMatrix;
#define DnaMatrix Wise2_DnaMatrix
#define DYNAMITE_DEFINED_DnaMatrix
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  new_CompMat_from_DnaMatrix_flat(dm)
 *
 * Descrip:    Builds a CompMat mapping of a DnaMatrix
 *
 *
 * Arg:        dm [UNKN ] Undocumented argument [DnaMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [CompMat *]
 *
 */
CompMat * Wise2_new_CompMat_from_DnaMatrix_flat(DnaMatrix * dm);
#define new_CompMat_from_DnaMatrix_flat Wise2_new_CompMat_from_DnaMatrix_flat


/* Function:  DnaProbMatrix_from_match(match,nmask_type)
 *
 * Descrip:    Makes a probability matrix from simple match/mismatch 
 *             probabilities.
 *
 *
 *
 * Arg:             match [UNKN ] Undocumented argument [Probability]
 * Arg:        nmask_type [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [DnaProbMatrix *]
 *
 */
DnaProbMatrix * Wise2_DnaProbMatrix_from_match(Probability match,int nmask_type);
#define DnaProbMatrix_from_match Wise2_DnaProbMatrix_from_match


/* Function:  flat_null_DnaProbMatrix(dpm)
 *
 * Descrip:    makes a odds of dpm via a 0.25 factor 
 *             into each base.
 *
 *
 * Arg:        dpm [UNKN ] Undocumented argument [DnaProbMatrix *]
 *
 */
void Wise2_flat_null_DnaProbMatrix(DnaProbMatrix * dpm);
#define flat_null_DnaProbMatrix Wise2_flat_null_DnaProbMatrix


/* Function:  DnaMatrix_from_DnaProbMatrix(dpm)
 *
 * Descrip:    Maps probabilities to scores
 *
 *
 * Arg:        dpm [UNKN ] Undocumented argument [DnaProbMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [DnaMatrix *]
 *
 */
DnaMatrix * Wise2_DnaMatrix_from_DnaProbMatrix(DnaProbMatrix * dpm);
#define DnaMatrix_from_DnaProbMatrix Wise2_DnaMatrix_from_DnaProbMatrix


/* Function:  show_DnaMatrix(dcm,ofp)
 *
 * Descrip:    Simple view of DnaMatrix
 *
 *
 * Arg:        dcm [UNKN ] Undocumented argument [DnaMatrix *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_show_DnaMatrix(DnaMatrix * dcm,FILE * ofp);
#define show_DnaMatrix Wise2_show_DnaMatrix


/* Function:  show_DnaProbMatrix(dpm,ofp)
 *
 * Descrip:    Simple view of DnaProbMatrix
 *
 *
 * Arg:        dpm [UNKN ] Undocumented argument [DnaProbMatrix *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_show_DnaProbMatrix(DnaProbMatrix * dpm,FILE * ofp);
#define show_DnaProbMatrix Wise2_show_DnaProbMatrix


/* Function:  fail_safe_DnaMatrix_access(dm,one,two)
 *
 * Descrip:    Run-time checked that one and two are ok to pass
 *             into dm as bases
 *
 *
 *
 * Arg:         dm [UNKN ] DnaMatrix to get score from [DnaMatrix *]
 * Arg:        one [UNKN ] base of one sequence [base]
 * Arg:        two [UNKN ] base of the other sequence [base]
 *
 * Return [UNKN ]  Undocumented return value [Score]
 *
 */
Score Wise2_fail_safe_DnaMatrix_access(DnaMatrix * dm,base one,base two);
#define fail_safe_DnaMatrix_access Wise2_fail_safe_DnaMatrix_access


/* Function:  identity_DnaMatrix(id_score,mismatch)
 *
 * Descrip:    makes an idenity matrix wth id_score on the leading
 *             diagonal and mismatch elsewhere.
 *
 *
 *
 * Arg:        id_score [UNKN ] score of idenities [Score]
 * Arg:        mismatch [UNKN ] score of mistmatches [Score]
 *
 * Return [UNKN ]  Undocumented return value [DnaMatrix *]
 *
 */
DnaMatrix * Wise2_identity_DnaMatrix(Score id_score,Score mismatch);
#define identity_DnaMatrix Wise2_identity_DnaMatrix


/* Function:  hard_link_DnaProbMatrix(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DnaProbMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [DnaProbMatrix *]
 *
 */
DnaProbMatrix * Wise2_hard_link_DnaProbMatrix(DnaProbMatrix * obj);
#define hard_link_DnaProbMatrix Wise2_hard_link_DnaProbMatrix


/* Function:  DnaProbMatrix_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaProbMatrix *]
 *
 */
DnaProbMatrix * Wise2_DnaProbMatrix_alloc(void);
#define DnaProbMatrix_alloc Wise2_DnaProbMatrix_alloc


/* Function:  free_DnaProbMatrix(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DnaProbMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [DnaProbMatrix *]
 *
 */
DnaProbMatrix * Wise2_free_DnaProbMatrix(DnaProbMatrix * obj);
#define free_DnaProbMatrix Wise2_free_DnaProbMatrix


/* Function:  hard_link_DnaMatrix(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DnaMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [DnaMatrix *]
 *
 */
DnaMatrix * Wise2_hard_link_DnaMatrix(DnaMatrix * obj);
#define hard_link_DnaMatrix Wise2_hard_link_DnaMatrix


/* Function:  DnaMatrix_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaMatrix *]
 *
 */
DnaMatrix * Wise2_DnaMatrix_alloc(void);
#define DnaMatrix_alloc Wise2_DnaMatrix_alloc


/* Function:  free_DnaMatrix(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DnaMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [DnaMatrix *]
 *
 */
DnaMatrix * Wise2_free_DnaMatrix(DnaMatrix * obj);
#define free_DnaMatrix Wise2_free_DnaMatrix


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/

#ifdef _cplusplus
}
#endif

#endif
