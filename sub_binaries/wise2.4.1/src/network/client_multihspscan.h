#ifndef DYNAMITEclient_multihspscanHEADERFILE
#define DYNAMITEclient_multihspscanHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "net_hspscan.h"

#define MultiHSPScanClientLISTLENGTH 32

#define MAX_CLIENT_HSP_THREADS 128

struct Wise2_HSPScanClient {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * host;     
    int port;    
    HSPScanInterface * hspi;     
    LinearHSPmanager * lm;   
    Sequence * seq;  
    HSPScanInterfacePara * p;    
    } ;  
/* HSPScanClient defined */ 
#ifndef DYNAMITE_DEFINED_HSPScanClient
typedef struct Wise2_HSPScanClient Wise2_HSPScanClient;
#define HSPScanClient Wise2_HSPScanClient
#define DYNAMITE_DEFINED_HSPScanClient
#endif


struct Wise2_MultiHSPScanClient {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    HSPScanClient ** client;     
    int len;/* len for above client  */ 
    int maxlen; /* maxlen for above client */ 
    } ;  
/* MultiHSPScanClient defined */ 
#ifndef DYNAMITE_DEFINED_MultiHSPScanClient
typedef struct Wise2_MultiHSPScanClient Wise2_MultiHSPScanClient;
#define MultiHSPScanClient Wise2_MultiHSPScanClient
#define DYNAMITE_DEFINED_MultiHSPScanClient
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_HSPScanClient(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [HSPScanClient *]
 *
 * Return [UNKN ]  Undocumented return value [HSPScanClient *]
 *
 */
HSPScanClient * Wise2_hard_link_HSPScanClient(HSPScanClient * obj);
#define hard_link_HSPScanClient Wise2_hard_link_HSPScanClient


/* Function:  HSPScanClient_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [HSPScanClient *]
 *
 */
HSPScanClient * Wise2_HSPScanClient_alloc(void);
#define HSPScanClient_alloc Wise2_HSPScanClient_alloc


/* Function:  free_HSPScanClient(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [HSPScanClient *]
 *
 * Return [UNKN ]  Undocumented return value [HSPScanClient *]
 *
 */
HSPScanClient * Wise2_free_HSPScanClient(HSPScanClient * obj);
#define free_HSPScanClient Wise2_free_HSPScanClient


/* Function:  add_MultiHSPScanClient(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [MultiHSPScanClient *]
 * Arg:        add [OWNER] Object to add to the list [HSPScanClient *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_MultiHSPScanClient(MultiHSPScanClient * obj,HSPScanClient * add);
#define add_MultiHSPScanClient Wise2_add_MultiHSPScanClient


/* Function:  flush_MultiHSPScanClient(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [MultiHSPScanClient *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_MultiHSPScanClient(MultiHSPScanClient * obj);
#define flush_MultiHSPScanClient Wise2_flush_MultiHSPScanClient


/* Function:  MultiHSPScanClient_alloc_std(void)
 *
 * Descrip:    Equivalent to MultiHSPScanClient_alloc_len(MultiHSPScanClientLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [MultiHSPScanClient *]
 *
 */
MultiHSPScanClient * Wise2_MultiHSPScanClient_alloc_std(void);
#define MultiHSPScanClient_alloc_std Wise2_MultiHSPScanClient_alloc_std


/* Function:  MultiHSPScanClient_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [MultiHSPScanClient *]
 *
 */
MultiHSPScanClient * Wise2_MultiHSPScanClient_alloc_len(int len);
#define MultiHSPScanClient_alloc_len Wise2_MultiHSPScanClient_alloc_len


/* Function:  hard_link_MultiHSPScanClient(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [MultiHSPScanClient *]
 *
 * Return [UNKN ]  Undocumented return value [MultiHSPScanClient *]
 *
 */
MultiHSPScanClient * Wise2_hard_link_MultiHSPScanClient(MultiHSPScanClient * obj);
#define hard_link_MultiHSPScanClient Wise2_hard_link_MultiHSPScanClient


/* Function:  MultiHSPScanClient_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [MultiHSPScanClient *]
 *
 */
MultiHSPScanClient * Wise2_MultiHSPScanClient_alloc(void);
#define MultiHSPScanClient_alloc Wise2_MultiHSPScanClient_alloc


/* Function:  free_MultiHSPScanClient(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [MultiHSPScanClient *]
 *
 * Return [UNKN ]  Undocumented return value [MultiHSPScanClient *]
 *
 */
MultiHSPScanClient * Wise2_free_MultiHSPScanClient(MultiHSPScanClient * obj);
#define free_MultiHSPScanClient Wise2_free_MultiHSPScanClient


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
LinearHSPmanager * Wise2_scan_MultiHSPScanClient(void * data,Sequence * seq,HSPScanInterfacePara * para);
#define scan_MultiHSPScanClient Wise2_scan_MultiHSPScanClient
void * Wise2_client_hsp_thread_worker(void * d);
#define client_hsp_thread_worker Wise2_client_hsp_thread_worker
HSPScanInterface * Wise2_new_multiclient_HSPScanInterface(char * filename);
#define new_multiclient_HSPScanInterface Wise2_new_multiclient_HSPScanInterface
void Wise2_free_multihspscanclient(void * data);
#define free_multihspscanclient Wise2_free_multihspscanclient
MultiHSPScanClient * Wise2_new_MultiHSPScanClient_from_file(FILE * ifp);
#define new_MultiHSPScanClient_from_file Wise2_new_MultiHSPScanClient_from_file


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_MultiHSPScanClient(HSPScanClient ** list,int i,int j) ;
#define swap_MultiHSPScanClient Wise2_swap_MultiHSPScanClient
void Wise2_qsort_MultiHSPScanClient(HSPScanClient ** list,int left,int right,int (*comp)(HSPScanClient * ,HSPScanClient * ));
#define qsort_MultiHSPScanClient Wise2_qsort_MultiHSPScanClient
void Wise2_sort_MultiHSPScanClient(MultiHSPScanClient * obj,int (*comp)(HSPScanClient *, HSPScanClient *));
#define sort_MultiHSPScanClient Wise2_sort_MultiHSPScanClient
boolean Wise2_expand_MultiHSPScanClient(MultiHSPScanClient * obj,int len);
#define expand_MultiHSPScanClient Wise2_expand_MultiHSPScanClient

#ifdef _cplusplus
}
#endif

#endif
