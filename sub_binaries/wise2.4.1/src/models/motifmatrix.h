#ifndef DYNAMITEmotifmatrixHEADERFILE
#define DYNAMITEmotifmatrixHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "transregion.h"


struct Wise2_MotifMatrixPara {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Probability comp_in_match;   
    Probability comp_out_match;  
    Probability comp_spacer;     
    Probability region_in;   
    Probability motif_indel;     
    Probability cons_indel;  
    Probability spacer_indel;    
    Probability spacer_to_cons;  
    Probability spacer_to_motif;     
    Probability spacer_duration;     
    Probability motif_duration;  
    Probability cons_duration;   
    } ;  
/* MotifMatrixPara defined */ 
#ifndef DYNAMITE_DEFINED_MotifMatrixPara
typedef struct Wise2_MotifMatrixPara Wise2_MotifMatrixPara;
#define MotifMatrixPara Wise2_MotifMatrixPara
#define DYNAMITE_DEFINED_MotifMatrixPara
#endif


struct Wise2_MotifMatrixScore {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    DnaMatrix * comp_in_motif;   
    DnaMatrix * comp_out_motif;  
    DnaMatrix * comp_spacer;     
    Score region_in;     
    Score motif_indel;   
    Score cons_indel;    
    Score spacer_indel;  
    Score spacer_to_cons;    
    Score spacer_to_motif;   
    Score spacer_duration;   
    Score motif_duration;    
    Score cons_duration;     
    } ;  
/* MotifMatrixScore defined */ 
#ifndef DYNAMITE_DEFINED_MotifMatrixScore
typedef struct Wise2_MotifMatrixScore Wise2_MotifMatrixScore;
#define MotifMatrixScore Wise2_MotifMatrixScore
#define DYNAMITE_DEFINED_MotifMatrixScore
#endif


struct Wise2_MotifConsMatrix {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char ** mat;     
    int leni;   /* leni for above mat  */ 
    int maxleni;/* max length for above pointer set */ 
    int lenj;   /* lenj for above mat  */ 
    int maxlenj;/* max length for above pointer set */ 
    } ;  
/* MotifConsMatrix defined */ 
#ifndef DYNAMITE_DEFINED_MotifConsMatrix
typedef struct Wise2_MotifConsMatrix Wise2_MotifConsMatrix;
#define MotifConsMatrix Wise2_MotifConsMatrix
#define DYNAMITE_DEFINED_MotifConsMatrix
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_MotifMatrixPara(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [MotifMatrixPara *]
 *
 * Return [UNKN ]  Undocumented return value [MotifMatrixPara *]
 *
 */
MotifMatrixPara * Wise2_hard_link_MotifMatrixPara(MotifMatrixPara * obj);
#define hard_link_MotifMatrixPara Wise2_hard_link_MotifMatrixPara


/* Function:  MotifMatrixPara_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [MotifMatrixPara *]
 *
 */
MotifMatrixPara * Wise2_MotifMatrixPara_alloc(void);
#define MotifMatrixPara_alloc Wise2_MotifMatrixPara_alloc


/* Function:  free_MotifMatrixPara(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [MotifMatrixPara *]
 *
 * Return [UNKN ]  Undocumented return value [MotifMatrixPara *]
 *
 */
MotifMatrixPara * Wise2_free_MotifMatrixPara(MotifMatrixPara * obj);
#define free_MotifMatrixPara Wise2_free_MotifMatrixPara


/* Function:  hard_link_MotifMatrixScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [MotifMatrixScore *]
 *
 * Return [UNKN ]  Undocumented return value [MotifMatrixScore *]
 *
 */
MotifMatrixScore * Wise2_hard_link_MotifMatrixScore(MotifMatrixScore * obj);
#define hard_link_MotifMatrixScore Wise2_hard_link_MotifMatrixScore


/* Function:  MotifMatrixScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [MotifMatrixScore *]
 *
 */
MotifMatrixScore * Wise2_MotifMatrixScore_alloc(void);
#define MotifMatrixScore_alloc Wise2_MotifMatrixScore_alloc


/* Function:  free_MotifMatrixScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [MotifMatrixScore *]
 *
 * Return [UNKN ]  Undocumented return value [MotifMatrixScore *]
 *
 */
MotifMatrixScore * Wise2_free_MotifMatrixScore(MotifMatrixScore * obj);
#define free_MotifMatrixScore Wise2_free_MotifMatrixScore


/* Function:  MotifConsMatrix_alloc_matrix(leni,lenj)
 *
 * Descrip:    Allocates structure and matrix
 *
 *
 * Arg:        leni [UNKN ] Length of first dimension of matrix [int]
 * Arg:        lenj [UNKN ] Length of second dimension of matrix [int]
 *
 * Return [UNKN ]  Undocumented return value [MotifConsMatrix *]
 *
 */
MotifConsMatrix * Wise2_MotifConsMatrix_alloc_matrix(int leni,int lenj);
#define MotifConsMatrix_alloc_matrix Wise2_MotifConsMatrix_alloc_matrix


/* Function:  hard_link_MotifConsMatrix(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [MotifConsMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [MotifConsMatrix *]
 *
 */
MotifConsMatrix * Wise2_hard_link_MotifConsMatrix(MotifConsMatrix * obj);
#define hard_link_MotifConsMatrix Wise2_hard_link_MotifConsMatrix


/* Function:  MotifConsMatrix_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [MotifConsMatrix *]
 *
 */
MotifConsMatrix * Wise2_MotifConsMatrix_alloc(void);
#define MotifConsMatrix_alloc Wise2_MotifConsMatrix_alloc


/* Function:  free_MotifConsMatrix(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [MotifConsMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [MotifConsMatrix *]
 *
 */
MotifConsMatrix * Wise2_free_MotifConsMatrix(MotifConsMatrix * obj);
#define free_MotifConsMatrix Wise2_free_MotifConsMatrix


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
MotifMatrixScore * Wise2_MotifMatrixScore_from_MotifMatrixPara(MotifMatrixPara * mmp);
#define MotifMatrixScore_from_MotifMatrixPara Wise2_MotifMatrixScore_from_MotifMatrixPara
MotifMatrixPara * Wise2_new_MotifMatrixPara_from_argv(int * argc,char ** argv);
#define new_MotifMatrixPara_from_argv Wise2_new_MotifMatrixPara_from_argv
void Wise2_show_help_MotifMatrixPara(FILE * ofp);
#define show_help_MotifMatrixPara Wise2_show_help_MotifMatrixPara
MotifConsMatrix * Wise2_new_MotifConsMatrix(TransFactorMatchSet * one,int starti,int endi,TransFactorMatchSet * two,int startj,int endj);
#define new_MotifConsMatrix Wise2_new_MotifConsMatrix


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
boolean Wise2_expand_MotifConsMatrix(MotifConsMatrix * obj,int leni,int lenj);
#define expand_MotifConsMatrix Wise2_expand_MotifConsMatrix

#ifdef _cplusplus
}
#endif

#endif
