#ifndef DYNAMITEshotgunHEADERFILE
#define DYNAMITEshotgunHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "sequence.h"


struct Wise2_ShotgunPara {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int read_length;     
    int insert_size;     
    int number;  
    int forward_only;    
    } ;  
/* ShotgunPara defined */ 
#ifndef DYNAMITE_DEFINED_ShotgunPara
typedef struct Wise2_ShotgunPara Wise2_ShotgunPara;
#define ShotgunPara Wise2_ShotgunPara
#define DYNAMITE_DEFINED_ShotgunPara
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  generate_shotgun_reads(shot,input,ofp)
 *
 * Descrip:    Generates a file of Shotgun reads from a particular
 *             sequence randomly
 *
 *
 * Arg:         shot [UNKN ] Undocumented argument [ShotgunPara *]
 * Arg:        input [UNKN ] Undocumented argument [Sequence *]
 * Arg:          ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_generate_shotgun_reads(ShotgunPara * shot,Sequence * input,FILE * ofp);
#define generate_shotgun_reads Wise2_generate_shotgun_reads


/* Function:  hard_link_ShotgunPara(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ShotgunPara *]
 *
 * Return [UNKN ]  Undocumented return value [ShotgunPara *]
 *
 */
ShotgunPara * Wise2_hard_link_ShotgunPara(ShotgunPara * obj);
#define hard_link_ShotgunPara Wise2_hard_link_ShotgunPara


/* Function:  ShotgunPara_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ShotgunPara *]
 *
 */
ShotgunPara * Wise2_ShotgunPara_alloc(void);
#define ShotgunPara_alloc Wise2_ShotgunPara_alloc


/* Function:  free_ShotgunPara(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ShotgunPara *]
 *
 * Return [UNKN ]  Undocumented return value [ShotgunPara *]
 *
 */
ShotgunPara * Wise2_free_ShotgunPara(ShotgunPara * obj);
#define free_ShotgunPara Wise2_free_ShotgunPara


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
