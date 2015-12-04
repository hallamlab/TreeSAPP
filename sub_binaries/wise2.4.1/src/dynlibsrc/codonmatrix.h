#ifndef DYNAMITEcodonmatrixHEADERFILE
#define DYNAMITEcodonmatrixHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"


struct Wise2_CodonMatrix {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Probability prob[125][125];  
    } ;  
/* CodonMatrix defined */ 
#ifndef DYNAMITE_DEFINED_CodonMatrix
typedef struct Wise2_CodonMatrix Wise2_CodonMatrix;
#define CodonMatrix Wise2_CodonMatrix
#define DYNAMITE_DEFINED_CodonMatrix
#endif


struct Wise2_CodonMatrixScore {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Score score[125][125];   
    } ;  
/* CodonMatrixScore defined */ 
#ifndef DYNAMITE_DEFINED_CodonMatrixScore
typedef struct Wise2_CodonMatrixScore Wise2_CodonMatrixScore;
#define CodonMatrixScore Wise2_CodonMatrixScore
#define DYNAMITE_DEFINED_CodonMatrixScore
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  naive_CodonMatrixScore_from_prob(ct,cm)
 *
 * Descrip:    Combines CodonMatrixScore_from_CodonMatrix and naive_CodonMatrix
 *
 *
 * Arg:        ct [UNKN ] Undocumented argument [CodonTable *]
 * Arg:        cm [UNKN ] Undocumented argument [CompProb *]
 *
 * Return [UNKN ]  Undocumented return value [CodonMatrixScore *]
 *
 */
CodonMatrixScore * Wise2_naive_CodonMatrixScore_from_prob(CodonTable * ct,CompProb * cm);
#define naive_CodonMatrixScore_from_prob Wise2_naive_CodonMatrixScore_from_prob


/* Function:  CodonMatrixScore_from_CodonMatrix(cm)
 *
 * Descrip:    Makes a CodonMatrixScore from a CodonMatrix
 *
 *
 * Arg:        cm [UNKN ] Undocumented argument [CodonMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [CodonMatrixScore *]
 *
 */
CodonMatrixScore * Wise2_CodonMatrixScore_from_CodonMatrix(CodonMatrix * cm);
#define CodonMatrixScore_from_CodonMatrix Wise2_CodonMatrixScore_from_CodonMatrix


/* Function:  naive_CodonMatrix(ct,comp)
 *
 * Descrip:    Builds a probability matrix
 *               No codon bias
 *               No errors
 *             N codons score 1.0, stop codons probability 0.00001
 *
 *
 * Arg:          ct [UNKN ] Undocumented argument [CodonTable *]
 * Arg:        comp [UNKN ] Undocumented argument [CompProb *]
 *
 * Return [UNKN ]  Undocumented return value [CodonMatrix *]
 *
 */
CodonMatrix * Wise2_naive_CodonMatrix(CodonTable * ct,CompProb * comp);
#define naive_CodonMatrix Wise2_naive_CodonMatrix


/* Function:  naive_CodonMatrixScore(ct,comp)
 *
 * Descrip:    Builds a codon matrix from CompMat which assummes:
 *               No codon Bias
 *               No errors
 *             N codons score 0, stop codons score ??
 *
 *
 * Arg:          ct [UNKN ] CodonTable for codon->aa mapping [CodonTable *]
 * Arg:        comp [UNKN ] Comparison matrix for the score of the individual access [CompMat *]
 *
 * Return [UNKN ]  Undocumented return value [CodonMatrixScore *]
 *
 */
CodonMatrixScore * Wise2_naive_CodonMatrixScore(CodonTable * ct,CompMat * comp);
#define naive_CodonMatrixScore Wise2_naive_CodonMatrixScore


/* Function:  show_CodonMatrixScore(cms,ct,ofp)
 *
 * Descrip:    Shows a codonmatrix
 *
 *
 * Arg:        cms [UNKN ] Undocumented argument [CodonMatrixScore *]
 * Arg:         ct [UNKN ] Undocumented argument [CodonTable *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_show_CodonMatrixScore(CodonMatrixScore * cms,CodonTable * ct,FILE * ofp);
#define show_CodonMatrixScore Wise2_show_CodonMatrixScore


/* Function:  hard_link_CodonMatrix(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [CodonMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [CodonMatrix *]
 *
 */
CodonMatrix * Wise2_hard_link_CodonMatrix(CodonMatrix * obj);
#define hard_link_CodonMatrix Wise2_hard_link_CodonMatrix


/* Function:  CodonMatrix_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [CodonMatrix *]
 *
 */
CodonMatrix * Wise2_CodonMatrix_alloc(void);
#define CodonMatrix_alloc Wise2_CodonMatrix_alloc


/* Function:  free_CodonMatrix(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [CodonMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [CodonMatrix *]
 *
 */
CodonMatrix * Wise2_free_CodonMatrix(CodonMatrix * obj);
#define free_CodonMatrix Wise2_free_CodonMatrix


/* Function:  hard_link_CodonMatrixScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [CodonMatrixScore *]
 *
 * Return [UNKN ]  Undocumented return value [CodonMatrixScore *]
 *
 */
CodonMatrixScore * Wise2_hard_link_CodonMatrixScore(CodonMatrixScore * obj);
#define hard_link_CodonMatrixScore Wise2_hard_link_CodonMatrixScore


/* Function:  CodonMatrixScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [CodonMatrixScore *]
 *
 */
CodonMatrixScore * Wise2_CodonMatrixScore_alloc(void);
#define CodonMatrixScore_alloc Wise2_CodonMatrixScore_alloc


/* Function:  free_CodonMatrixScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [CodonMatrixScore *]
 *
 * Return [UNKN ]  Undocumented return value [CodonMatrixScore *]
 *
 */
CodonMatrixScore * Wise2_free_CodonMatrixScore(CodonMatrixScore * obj);
#define free_CodonMatrixScore Wise2_free_CodonMatrixScore


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
