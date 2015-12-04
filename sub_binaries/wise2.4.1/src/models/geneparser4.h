#ifndef DYNAMITEgeneparser4HEADERFILE
#define DYNAMITEgeneparser4HEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "wisebase.h"
#include "probability.h"
#include "geneparser21.h"

enum GeneParser4Type {
  GP4_INTRON2CDS = 0,
  GP4_INTRON2INTRON,
  GP4_DELETE_1_BASE,
  GP4_DELETE_2_BASE,
  GP4_INSERT_1_BASE,
  GP4_INSERT_2_BASE,
  GP4_LOOP2LOOP,
  GP4_LOOP2MODEL,
  GP4_TRANSITION_LEN };

struct Wise2_GeneParser4 {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Probability transition[GP4_TRANSITION_LEN];  
    Probability intron[5];   
    } ;  
/* GeneParser4 defined */ 
#ifndef DYNAMITE_DEFINED_GeneParser4
typedef struct Wise2_GeneParser4 Wise2_GeneParser4;
#define GeneParser4 Wise2_GeneParser4
#define DYNAMITE_DEFINED_GeneParser4
#endif


struct Wise2_GeneParser4Score {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Score transition[GP4_TRANSITION_LEN];    
    Score intron[5];     
    } ;  
/* GeneParser4Score defined */ 
#ifndef DYNAMITE_DEFINED_GeneParser4Score
typedef struct Wise2_GeneParser4Score Wise2_GeneParser4Score;
#define GeneParser4Score Wise2_GeneParser4Score
#define DYNAMITE_DEFINED_GeneParser4Score
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_GeneParser4(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GeneParser4 *]
 *
 * Return [UNKN ]  Undocumented return value [GeneParser4 *]
 *
 */
GeneParser4 * Wise2_hard_link_GeneParser4(GeneParser4 * obj);
#define hard_link_GeneParser4 Wise2_hard_link_GeneParser4


/* Function:  GeneParser4_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneParser4 *]
 *
 */
GeneParser4 * Wise2_GeneParser4_alloc(void);
#define GeneParser4_alloc Wise2_GeneParser4_alloc


/* Function:  free_GeneParser4(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneParser4 *]
 *
 * Return [UNKN ]  Undocumented return value [GeneParser4 *]
 *
 */
GeneParser4 * Wise2_free_GeneParser4(GeneParser4 * obj);
#define free_GeneParser4 Wise2_free_GeneParser4


/* Function:  hard_link_GeneParser4Score(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GeneParser4Score *]
 *
 * Return [UNKN ]  Undocumented return value [GeneParser4Score *]
 *
 */
GeneParser4Score * Wise2_hard_link_GeneParser4Score(GeneParser4Score * obj);
#define hard_link_GeneParser4Score Wise2_hard_link_GeneParser4Score


/* Function:  GeneParser4Score_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneParser4Score *]
 *
 */
GeneParser4Score * Wise2_GeneParser4Score_alloc(void);
#define GeneParser4Score_alloc Wise2_GeneParser4Score_alloc


/* Function:  free_GeneParser4Score(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneParser4Score *]
 *
 * Return [UNKN ]  Undocumented return value [GeneParser4Score *]
 *
 */
GeneParser4Score * Wise2_free_GeneParser4Score(GeneParser4Score * obj);
#define free_GeneParser4Score Wise2_free_GeneParser4Score


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
GeneParser4 * Wise2_std_GeneParser4(double indel,double intron2cds);
#define std_GeneParser4 Wise2_std_GeneParser4
GeneParser4Score * Wise2_GeneParser4Score_from_GeneParser21Score(GeneParser21Score * gp21s);
#define GeneParser4Score_from_GeneParser21Score Wise2_GeneParser4Score_from_GeneParser21Score
GeneParser4Score * Wise2_GeneParser4Score_from_GeneParser4(GeneParser4 * gp4);
#define GeneParser4Score_from_GeneParser4 Wise2_GeneParser4Score_from_GeneParser4


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/

#ifdef _cplusplus
}
#endif

#endif
