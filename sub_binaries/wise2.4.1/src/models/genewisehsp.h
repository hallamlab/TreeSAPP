#ifndef DYNAMITEgenewisehspHEADERFILE
#define DYNAMITEgenewisehspHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"
#include "subseqhash.h"
#include "hsplookupscan.h"


#define GeneWiseHSPmanagerLISTLENGTH 1024

struct Wise2_GeneWiseRunPara {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    boolean use_hsp;     
    int edge_query;  
    int edge_target;     
    int splice_spread;   
    boolean debug;   
    } ;  
/* GeneWiseRunPara defined */ 
#ifndef DYNAMITE_DEFINED_GeneWiseRunPara
typedef struct Wise2_GeneWiseRunPara Wise2_GeneWiseRunPara;
#define GeneWiseRunPara Wise2_GeneWiseRunPara
#define DYNAMITE_DEFINED_GeneWiseRunPara
#endif


struct Wise2_GeneWiseHSP {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int query_start;     
    int query_end;   
    int target_start;    
    int target_end;  
    int score;   
    int frame;   
    } ;  
/* GeneWiseHSP defined */ 
#ifndef DYNAMITE_DEFINED_GeneWiseHSP
typedef struct Wise2_GeneWiseHSP Wise2_GeneWiseHSP;
#define GeneWiseHSP Wise2_GeneWiseHSP
#define DYNAMITE_DEFINED_GeneWiseHSP
#endif


struct Wise2_GeneWiseHSPmanager {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    GeneWiseHSP ** hsp;  
    int len;/* len for above hsp  */ 
    int maxlen; /* maxlen for above hsp */ 
    } ;  
/* GeneWiseHSPmanager defined */ 
#ifndef DYNAMITE_DEFINED_GeneWiseHSPmanager
typedef struct Wise2_GeneWiseHSPmanager Wise2_GeneWiseHSPmanager;
#define GeneWiseHSPmanager Wise2_GeneWiseHSPmanager
#define DYNAMITE_DEFINED_GeneWiseHSPmanager
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_GeneWiseRunPara(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GeneWiseRunPara *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseRunPara *]
 *
 */
GeneWiseRunPara * Wise2_hard_link_GeneWiseRunPara(GeneWiseRunPara * obj);
#define hard_link_GeneWiseRunPara Wise2_hard_link_GeneWiseRunPara


/* Function:  GeneWiseRunPara_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseRunPara *]
 *
 */
GeneWiseRunPara * Wise2_GeneWiseRunPara_alloc(void);
#define GeneWiseRunPara_alloc Wise2_GeneWiseRunPara_alloc


/* Function:  free_GeneWiseRunPara(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneWiseRunPara *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseRunPara *]
 *
 */
GeneWiseRunPara * Wise2_free_GeneWiseRunPara(GeneWiseRunPara * obj);
#define free_GeneWiseRunPara Wise2_free_GeneWiseRunPara


/* Function:  hard_link_GeneWiseHSP(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GeneWiseHSP *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseHSP *]
 *
 */
GeneWiseHSP * Wise2_hard_link_GeneWiseHSP(GeneWiseHSP * obj);
#define hard_link_GeneWiseHSP Wise2_hard_link_GeneWiseHSP


/* Function:  GeneWiseHSP_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseHSP *]
 *
 */
GeneWiseHSP * Wise2_GeneWiseHSP_alloc(void);
#define GeneWiseHSP_alloc Wise2_GeneWiseHSP_alloc


/* Function:  free_GeneWiseHSP(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneWiseHSP *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseHSP *]
 *
 */
GeneWiseHSP * Wise2_free_GeneWiseHSP(GeneWiseHSP * obj);
#define free_GeneWiseHSP Wise2_free_GeneWiseHSP


/* Function:  add_GeneWiseHSPmanager(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GeneWiseHSPmanager *]
 * Arg:        add [OWNER] Object to add to the list [GeneWiseHSP *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_GeneWiseHSPmanager(GeneWiseHSPmanager * obj,GeneWiseHSP * add);
#define add_GeneWiseHSPmanager Wise2_add_GeneWiseHSPmanager


/* Function:  flush_GeneWiseHSPmanager(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [GeneWiseHSPmanager *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_GeneWiseHSPmanager(GeneWiseHSPmanager * obj);
#define flush_GeneWiseHSPmanager Wise2_flush_GeneWiseHSPmanager


/* Function:  GeneWiseHSPmanager_alloc_std(void)
 *
 * Descrip:    Equivalent to GeneWiseHSPmanager_alloc_len(GeneWiseHSPmanagerLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseHSPmanager *]
 *
 */
GeneWiseHSPmanager * Wise2_GeneWiseHSPmanager_alloc_std(void);
#define GeneWiseHSPmanager_alloc_std Wise2_GeneWiseHSPmanager_alloc_std


/* Function:  GeneWiseHSPmanager_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseHSPmanager *]
 *
 */
GeneWiseHSPmanager * Wise2_GeneWiseHSPmanager_alloc_len(int len);
#define GeneWiseHSPmanager_alloc_len Wise2_GeneWiseHSPmanager_alloc_len


/* Function:  hard_link_GeneWiseHSPmanager(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GeneWiseHSPmanager *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseHSPmanager *]
 *
 */
GeneWiseHSPmanager * Wise2_hard_link_GeneWiseHSPmanager(GeneWiseHSPmanager * obj);
#define hard_link_GeneWiseHSPmanager Wise2_hard_link_GeneWiseHSPmanager


/* Function:  GeneWiseHSPmanager_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseHSPmanager *]
 *
 */
GeneWiseHSPmanager * Wise2_GeneWiseHSPmanager_alloc(void);
#define GeneWiseHSPmanager_alloc Wise2_GeneWiseHSPmanager_alloc


/* Function:  free_GeneWiseHSPmanager(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneWiseHSPmanager *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseHSPmanager *]
 *
 */
GeneWiseHSPmanager * Wise2_free_GeneWiseHSPmanager(GeneWiseHSPmanager * obj);
#define free_GeneWiseHSPmanager Wise2_free_GeneWiseHSPmanager


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
GeneWiseRunPara * Wise2_new_GeneWiseRunPara_from_argv(int * argc,char ** argv);
#define new_GeneWiseRunPara_from_argv Wise2_new_GeneWiseRunPara_from_argv
void Wise2_show_help_GeneWiseRunPara(FILE * ofp);
#define show_help_GeneWiseRunPara Wise2_show_help_GeneWiseRunPara
DPEnvelope * Wise2_DPEnvelope_from_protein_gen(Sequence * prot,Sequence * dna,CompMat * mat,CodonTable * ct,GeneWiseRunPara *p);
#define DPEnvelope_from_protein_gen Wise2_DPEnvelope_from_protein_gen
int Wise2_consistent_GeneWiseHSP(GeneWiseHSP * true,GeneWiseHSP * proposed);
#define consistent_GeneWiseHSP Wise2_consistent_GeneWiseHSP
int Wise2_compare_GeneWiseHSP_start(GeneWiseHSP * one,GeneWiseHSP * two);
#define compare_GeneWiseHSP_start Wise2_compare_GeneWiseHSP_start
void Wise2_add_GeneWiseHSPmanager_HSPset(GeneWiseHSPmanager * gwh,HSPset * set,int frame);
#define add_GeneWiseHSPmanager_HSPset Wise2_add_GeneWiseHSPmanager_HSPset


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
int Wise2_compare_GeneWiseHSP_score(GeneWiseHSP * one,GeneWiseHSP * two);
#define compare_GeneWiseHSP_score Wise2_compare_GeneWiseHSP_score
void Wise2_swap_GeneWiseHSPmanager(GeneWiseHSP ** list,int i,int j) ;
#define swap_GeneWiseHSPmanager Wise2_swap_GeneWiseHSPmanager
void Wise2_qsort_GeneWiseHSPmanager(GeneWiseHSP ** list,int left,int right,int (*comp)(GeneWiseHSP * ,GeneWiseHSP * ));
#define qsort_GeneWiseHSPmanager Wise2_qsort_GeneWiseHSPmanager
void Wise2_sort_GeneWiseHSPmanager(GeneWiseHSPmanager * obj,int (*comp)(GeneWiseHSP *, GeneWiseHSP *));
#define sort_GeneWiseHSPmanager Wise2_sort_GeneWiseHSPmanager
boolean Wise2_expand_GeneWiseHSPmanager(GeneWiseHSPmanager * obj,int len);
#define expand_GeneWiseHSPmanager Wise2_expand_GeneWiseHSPmanager

#ifdef _cplusplus
}
#endif

#endif
