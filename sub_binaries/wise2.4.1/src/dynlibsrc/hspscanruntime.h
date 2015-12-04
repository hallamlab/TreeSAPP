#ifndef DYNAMITEhspscanruntimeHEADERFILE
#define DYNAMITEhspscanruntimeHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif

#include "hsplookupscan.h"
#include "hsplookupthreaded.h"
#include "hsptwohitscan.h"


struct Wise2_HSPScanRuntimeImpl {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    HSPScanInterface * vanilla;  
    HSPScanInterface * threaded;     
    HSPScanInterface * twohit;   
    } ;  
/* HSPScanRuntimeImpl defined */ 
#ifndef DYNAMITE_DEFINED_HSPScanRuntimeImpl
typedef struct Wise2_HSPScanRuntimeImpl Wise2_HSPScanRuntimeImpl;
#define HSPScanRuntimeImpl Wise2_HSPScanRuntimeImpl
#define DYNAMITE_DEFINED_HSPScanRuntimeImpl
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  new_runtime_HSPScanInterface(sli,mat,drop_off,score_cutoff,threadno)
 *
 * Descrip:    Makes a new function which will at runtime switch between
 *             implementation; vanilla, threaded and twohit
 *
 *
 * Arg:                 sli [UNKN ] Undocumented argument [SeqLookupInterface *]
 * Arg:                 mat [UNKN ] Undocumented argument [CompMat *]
 * Arg:            drop_off [UNKN ] Undocumented argument [int]
 * Arg:        score_cutoff [UNKN ] Undocumented argument [int]
 * Arg:            threadno [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [HSPScanInterface *]
 *
 */
HSPScanInterface * Wise2_new_runtime_HSPScanInterface(SeqLookupInterface * sli,CompMat * mat,int drop_off,int score_cutoff,int threadno);
#define new_runtime_HSPScanInterface Wise2_new_runtime_HSPScanInterface


/* Function:  scan_query_runtime_hspscan(data,seq,para)
 *
 * Descrip:    Handles runtime switching between methods for the scan query
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 * Arg:         seq [UNKN ] Undocumented argument [Sequence *]
 * Arg:        para [UNKN ] Undocumented argument [HSPScanInterfacePara *]
 *
 * Return [UNKN ]  Undocumented return value [LinearHSPmanager *]
 *
 */
LinearHSPmanager * Wise2_scan_query_runtime_hspscan(void * data,Sequence * seq,HSPScanInterfacePara * para);
#define scan_query_runtime_hspscan Wise2_scan_query_runtime_hspscan


/* Function:  free_runtime_hspscan(data)
 *
 * Descrip:    free function for runtime
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 *
 */
void Wise2_free_runtime_hspscan(void * data);
#define free_runtime_hspscan Wise2_free_runtime_hspscan


/* Function:  hard_link_HSPScanRuntimeImpl(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [HSPScanRuntimeImpl *]
 *
 * Return [UNKN ]  Undocumented return value [HSPScanRuntimeImpl *]
 *
 */
HSPScanRuntimeImpl * Wise2_hard_link_HSPScanRuntimeImpl(HSPScanRuntimeImpl * obj);
#define hard_link_HSPScanRuntimeImpl Wise2_hard_link_HSPScanRuntimeImpl


/* Function:  HSPScanRuntimeImpl_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [HSPScanRuntimeImpl *]
 *
 */
HSPScanRuntimeImpl * Wise2_HSPScanRuntimeImpl_alloc(void);
#define HSPScanRuntimeImpl_alloc Wise2_HSPScanRuntimeImpl_alloc


/* Function:  free_HSPScanRuntimeImpl(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [HSPScanRuntimeImpl *]
 *
 * Return [UNKN ]  Undocumented return value [HSPScanRuntimeImpl *]
 *
 */
HSPScanRuntimeImpl * Wise2_free_HSPScanRuntimeImpl(HSPScanRuntimeImpl * obj);
#define free_HSPScanRuntimeImpl Wise2_free_HSPScanRuntimeImpl


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
