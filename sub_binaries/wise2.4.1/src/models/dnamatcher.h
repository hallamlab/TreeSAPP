#ifndef DYNAMITEdnamatcherHEADERFILE
#define DYNAMITEdnamatcherHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dnaalign.h"
#include "hitlist.h"

struct Wise2_DnaMatchPara {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    DPRunImpl * dpri;    
    DnaMatrix * mat;     
    DnaStartEnd * dse;   
    int gap;     
    int ext;     
    } ;  
/* DnaMatchPara defined */ 
#ifndef DYNAMITE_DEFINED_DnaMatchPara
typedef struct Wise2_DnaMatchPara Wise2_DnaMatchPara;
#define DnaMatchPara Wise2_DnaMatchPara
#define DYNAMITE_DEFINED_DnaMatchPara
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_DnaMatchPara(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DnaMatchPara *]
 *
 * Return [UNKN ]  Undocumented return value [DnaMatchPara *]
 *
 */
DnaMatchPara * Wise2_hard_link_DnaMatchPara(DnaMatchPara * obj);
#define hard_link_DnaMatchPara Wise2_hard_link_DnaMatchPara


/* Function:  DnaMatchPara_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaMatchPara *]
 *
 */
DnaMatchPara * Wise2_DnaMatchPara_alloc(void);
#define DnaMatchPara_alloc Wise2_DnaMatchPara_alloc


/* Function:  free_DnaMatchPara(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DnaMatchPara *]
 *
 * Return [UNKN ]  Undocumented return value [DnaMatchPara *]
 *
 */
DnaMatchPara * Wise2_free_DnaMatchPara(DnaMatchPara * obj);
#define free_DnaMatchPara Wise2_free_DnaMatchPara


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
void Wise2_show_help_DnaMatchPara(FILE * ofp);
#define show_help_DnaMatchPara Wise2_show_help_DnaMatchPara
DnaMatchPara * Wise2_new_DnaMatchPara_from_argv(int * argc,char ** argv);
#define new_DnaMatchPara_from_argv Wise2_new_DnaMatchPara_from_argv
HitList * Wise2_HitList_from_Sequence_SequenceSet_DNA(Sequence * query,SequenceSet * set,DnaMatchPara * p);
#define HitList_from_Sequence_SequenceSet_DNA Wise2_HitList_from_Sequence_SequenceSet_DNA


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/

#ifdef _cplusplus
}
#endif

#endif
