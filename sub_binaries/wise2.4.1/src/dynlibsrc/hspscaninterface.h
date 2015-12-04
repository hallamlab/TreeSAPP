#ifndef DYNAMITEhspscaninterfaceHEADERFILE
#define DYNAMITEhspscaninterfaceHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "hsp.h"

#define HSPSCAN_IMPLEMENTATION_VANILLA  10
#define HSPSCAN_IMPLEMENTATION_THREADED 11
#define HSPSCAN_IMPLEMENTATION_TWOHIT   12


struct Wise2_HSPScanInterfacePara {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int min_score;   
    int max_results;     
    int flags;   
    boolean use_protein_heuristic;   
    int numb_level;  
    int word_depth;  
    int min_word_score;  
    int min_hsp_score;   
    int implementation;  
    int hsp_link_width;  
    int hsp_link_length;     
    int verbosity;   
    int low_numb;    
    int hsp_avg_ext;     
    } ;  
/* HSPScanInterfacePara defined */ 
#ifndef DYNAMITE_DEFINED_HSPScanInterfacePara
typedef struct Wise2_HSPScanInterfacePara Wise2_HSPScanInterfacePara;
#define HSPScanInterfacePara Wise2_HSPScanInterfacePara
#define DYNAMITE_DEFINED_HSPScanInterfacePara
#endif


struct Wise2_HSPScanInterface {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    LinearHSPmanager * (*scan_query)(void *,Sequence *,HSPScanInterfacePara * para); 
    void (*free_data)(void*);    
    void * data;     
    } ;  
/* HSPScanInterface defined */ 
#ifndef DYNAMITE_DEFINED_HSPScanInterface
typedef struct Wise2_HSPScanInterface Wise2_HSPScanInterface;
#define HSPScanInterface Wise2_HSPScanInterface
#define DYNAMITE_DEFINED_HSPScanInterface
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  untyped_read_HSPScanInterfacePara_from_Stream(rs)
 *
 * Descrip:    untyped version of read
 *
 *
 * Arg:        rs [UNKN ] Undocumented argument [Wise2ReadStreamInterface *]
 *
 * Return [UNKN ]  Undocumented return value [void *]
 *
 */
void * Wise2_untyped_read_HSPScanInterfacePara_from_Stream(Wise2ReadStreamInterface * rs);
#define untyped_read_HSPScanInterfacePara_from_Stream Wise2_untyped_read_HSPScanInterfacePara_from_Stream


/* Function:  typed_HSPScanInterfacePara_from_Stream(rs)
 *
 * Descrip:    reads a HSPscan interface para from a stream
 *
 *
 * Arg:        rs [UNKN ] Undocumented argument [Wise2ReadStreamInterface *]
 *
 * Return [UNKN ]  Undocumented return value [HSPScanInterfacePara *]
 *
 */
HSPScanInterfacePara * Wise2_typed_HSPScanInterfacePara_from_Stream(Wise2ReadStreamInterface * rs);
#define typed_HSPScanInterfacePara_from_Stream Wise2_typed_HSPScanInterfacePara_from_Stream


/* Function:  untyped_HSPScanInterfacePara_to_Stream(p,ws)
 *
 * Descrip:    writes out a standard HSPscan interface
 *
 *
 * Arg:         p [UNKN ] Undocumented argument [void *]
 * Arg:        ws [UNKN ] Undocumented argument [Wise2WriteStreamInterface *]
 *
 */
void Wise2_untyped_HSPScanInterfacePara_to_Stream(void * p,Wise2WriteStreamInterface * ws);
#define untyped_HSPScanInterfacePara_to_Stream Wise2_untyped_HSPScanInterfacePara_to_Stream


/* Function:  untyped_HSPScanInterfacePara_free(p)
 *
 * Descrip:    untyped free function for AnonymousObjects
 *
 *
 * Arg:        p [UNKN ] Undocumented argument [void *]
 *
 */
void Wise2_untyped_HSPScanInterfacePara_free(void * p);
#define untyped_HSPScanInterfacePara_free Wise2_untyped_HSPScanInterfacePara_free


/* Function:  show_help_HSPScanInterfacePara(ofp)
 *
 * Descrip:    help function for hspscan interface para
 *
 *
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_show_help_HSPScanInterfacePara(FILE * ofp);
#define show_help_HSPScanInterfacePara Wise2_show_help_HSPScanInterfacePara


/* Function:  new_HSPScanInterfacePara_from_argv(argc,argv)
 *
 * Descrip:    makes a hspscan interface from the command line
 *
 *
 * Arg:        argc [UNKN ] Undocumented argument [int *]
 * Arg:        argv [UNKN ] Undocumented argument [char **]
 *
 * Return [UNKN ]  Undocumented return value [HSPScanInterfacePara *]
 *
 */
HSPScanInterfacePara * Wise2_new_HSPScanInterfacePara_from_argv(int * argc,char ** argv);
#define new_HSPScanInterfacePara_from_argv Wise2_new_HSPScanInterfacePara_from_argv


/* Function:  HSPScanInterfacePara_std(void)
 *
 * Descrip:    makes a standard hsp scan interface para
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [HSPScanInterfacePara *]
 *
 */
HSPScanInterfacePara * Wise2_HSPScanInterfacePara_std(void);
#define HSPScanInterfacePara_std Wise2_HSPScanInterfacePara_std


/* Function:  free_HSPScanInterface(hsi)
 *
 * Descrip:    Frees overrides dynamite default
 *
 *
 * Arg:        hsi [UNKN ] Undocumented argument [HSPScanInterface *]
 *
 * Return [UNKN ]  Undocumented return value [HSPScanInterface *]
 *
 */
HSPScanInterface * Wise2_free_HSPScanInterface(HSPScanInterface * hsi);
#define free_HSPScanInterface Wise2_free_HSPScanInterface


/* Function:  hard_link_HSPScanInterfacePara(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [HSPScanInterfacePara *]
 *
 * Return [UNKN ]  Undocumented return value [HSPScanInterfacePara *]
 *
 */
HSPScanInterfacePara * Wise2_hard_link_HSPScanInterfacePara(HSPScanInterfacePara * obj);
#define hard_link_HSPScanInterfacePara Wise2_hard_link_HSPScanInterfacePara


/* Function:  HSPScanInterfacePara_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [HSPScanInterfacePara *]
 *
 */
HSPScanInterfacePara * Wise2_HSPScanInterfacePara_alloc(void);
#define HSPScanInterfacePara_alloc Wise2_HSPScanInterfacePara_alloc


/* Function:  free_HSPScanInterfacePara(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [HSPScanInterfacePara *]
 *
 * Return [UNKN ]  Undocumented return value [HSPScanInterfacePara *]
 *
 */
HSPScanInterfacePara * Wise2_free_HSPScanInterfacePara(HSPScanInterfacePara * obj);
#define free_HSPScanInterfacePara Wise2_free_HSPScanInterfacePara


/* Function:  hard_link_HSPScanInterface(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [HSPScanInterface *]
 *
 * Return [UNKN ]  Undocumented return value [HSPScanInterface *]
 *
 */
HSPScanInterface * Wise2_hard_link_HSPScanInterface(HSPScanInterface * obj);
#define hard_link_HSPScanInterface Wise2_hard_link_HSPScanInterface


/* Function:  HSPScanInterface_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [HSPScanInterface *]
 *
 */
HSPScanInterface * Wise2_HSPScanInterface_alloc(void);
#define HSPScanInterface_alloc Wise2_HSPScanInterface_alloc


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
