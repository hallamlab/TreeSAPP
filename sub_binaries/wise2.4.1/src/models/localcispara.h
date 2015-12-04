#ifndef DYNAMITElocalcisparaHEADERFILE
#define DYNAMITElocalcisparaHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"


#define LCH_GAP 0.05   
#define LCH_BLOCKOPEN  0.01   
#define LCH_UNMATCHED_PEN 0.99  


struct Wise2_LocalCisHitProb {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    DnaProbMatrix* comp65;   
    DnaProbMatrix* comp75;   
    DnaProbMatrix* comp85;   
    DnaProbMatrix* comp95;   
    Probability g;   
    Probability u;   
    Probability v;   
    Probability s;   
    Probability b;   
    } ;  
/* LocalCisHitProb defined */ 
#ifndef DYNAMITE_DEFINED_LocalCisHitProb
typedef struct Wise2_LocalCisHitProb Wise2_LocalCisHitProb;
#define LocalCisHitProb Wise2_LocalCisHitProb
#define DYNAMITE_DEFINED_LocalCisHitProb
#endif


struct Wise2_LocalCisHitScore {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    DnaMatrix* comp65;   
    DnaMatrix* comp75;   
    DnaMatrix* comp85;   
    DnaMatrix* comp95;   
    Score g;     
    Score u;     
    Score v;     
    Score s;     
    Score b;     
    } ;  
/* LocalCisHitScore defined */ 
#ifndef DYNAMITE_DEFINED_LocalCisHitScore
typedef struct Wise2_LocalCisHitScore Wise2_LocalCisHitScore;
#define LocalCisHitScore Wise2_LocalCisHitScore
#define DYNAMITE_DEFINED_LocalCisHitScore
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_LocalCisHitProb(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [LocalCisHitProb *]
 *
 * Return [UNKN ]  Undocumented return value [LocalCisHitProb *]
 *
 */
LocalCisHitProb * Wise2_hard_link_LocalCisHitProb(LocalCisHitProb * obj);
#define hard_link_LocalCisHitProb Wise2_hard_link_LocalCisHitProb


/* Function:  LocalCisHitProb_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LocalCisHitProb *]
 *
 */
LocalCisHitProb * Wise2_LocalCisHitProb_alloc(void);
#define LocalCisHitProb_alloc Wise2_LocalCisHitProb_alloc


/* Function:  free_LocalCisHitProb(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [LocalCisHitProb *]
 *
 * Return [UNKN ]  Undocumented return value [LocalCisHitProb *]
 *
 */
LocalCisHitProb * Wise2_free_LocalCisHitProb(LocalCisHitProb * obj);
#define free_LocalCisHitProb Wise2_free_LocalCisHitProb


/* Function:  hard_link_LocalCisHitScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [LocalCisHitScore *]
 *
 * Return [UNKN ]  Undocumented return value [LocalCisHitScore *]
 *
 */
LocalCisHitScore * Wise2_hard_link_LocalCisHitScore(LocalCisHitScore * obj);
#define hard_link_LocalCisHitScore Wise2_hard_link_LocalCisHitScore


/* Function:  LocalCisHitScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LocalCisHitScore *]
 *
 */
LocalCisHitScore * Wise2_LocalCisHitScore_alloc(void);
#define LocalCisHitScore_alloc Wise2_LocalCisHitScore_alloc


/* Function:  free_LocalCisHitScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [LocalCisHitScore *]
 *
 * Return [UNKN ]  Undocumented return value [LocalCisHitScore *]
 *
 */
LocalCisHitScore * Wise2_free_LocalCisHitScore(LocalCisHitScore * obj);
#define free_LocalCisHitScore Wise2_free_LocalCisHitScore


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
LocalCisHitScore * Wise2_standard_LocalCisHitScore(NMaskType nmask);
#define standard_LocalCisHitScore Wise2_standard_LocalCisHitScore
LocalCisHitProb * Wise2_standard_LocalCisHitProb(NMaskType nmask);
#define standard_LocalCisHitProb Wise2_standard_LocalCisHitProb
LocalCisHitScore * Wise2_LocalCisHitScore_from_LocalCisHitProb(LocalCisHitProb * lchp);
#define LocalCisHitScore_from_LocalCisHitProb Wise2_LocalCisHitScore_from_LocalCisHitProb


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/

#ifdef _cplusplus
}
#endif

#endif
