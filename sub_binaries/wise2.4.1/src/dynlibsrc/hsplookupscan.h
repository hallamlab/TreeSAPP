#ifndef DYNAMITEhsplookupscanHEADERFILE
#define DYNAMITEhsplookupscanHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "hspscaninterface.h"
#include "seqlookup.h"
#include "hsphandler.h"
#include "glib.h"




struct Wise2_HSPScanPara {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    SeqLookupInterface * sli;    
    CompMat * mat;   
    int drop_off;    
    int score_cutoff;    
    int use_msp_crunch;  
    int msp_crunch_no;   
    int seed_factor;     
    int twohit_wobble;   
    int threadno;    
    } ;  
/* HSPScanPara defined */ 
#ifndef DYNAMITE_DEFINED_HSPScanPara
typedef struct Wise2_HSPScanPara Wise2_HSPScanPara;
#define HSPScanPara Wise2_HSPScanPara
#define DYNAMITE_DEFINED_HSPScanPara
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  seq_number_aa_5mer_client(seq)
 *
 * Descrip:    Function for the amino acid to number on 5mers
 *
 *
 * Arg:        seq [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_seq_number_aa_5mer_client(char * seq);
#define seq_number_aa_5mer_client Wise2_seq_number_aa_5mer_client


/* Function:  new_simple_HSPScanInterface(sli,mat,drop_off)
 *
 * Descrip:    Builds a new simple scan interface. This
 *             does not expand the query using a matrix but
 *             rather simply scans down the query sequence
 *
 *
 * Arg:             sli [UNKN ] Undocumented argument [SeqLookupInterface *]
 * Arg:             mat [UNKN ] Undocumented argument [CompMat *]
 * Arg:        drop_off [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [HSPScanInterface *]
 *
 */
HSPScanInterface * Wise2_new_simple_HSPScanInterface(SeqLookupInterface * sli,CompMat * mat,int drop_off);
#define new_simple_HSPScanInterface Wise2_new_simple_HSPScanInterface


/* Function:  new_one_off_HSPScanInterface(sli,mat,drop_off,score_cutoff)
 *
 * Descrip:    Builds a new simple scan interface. This
 *             does expands the query using a matrix but
 *             just be considering off by one cases
 *
 *
 * Arg:                 sli [UNKN ] Undocumented argument [SeqLookupInterface *]
 * Arg:                 mat [UNKN ] Undocumented argument [CompMat *]
 * Arg:            drop_off [UNKN ] Undocumented argument [int]
 * Arg:        score_cutoff [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [HSPScanInterface *]
 *
 */
HSPScanInterface * Wise2_new_one_off_HSPScanInterface(SeqLookupInterface * sli,CompMat * mat,int drop_off,int score_cutoff);
#define new_one_off_HSPScanInterface Wise2_new_one_off_HSPScanInterface


/* Function:  no_op_func(data,user_data,data2)
 *
 * Descrip:    provide a no op func
 *
 *
 * Arg:             data [UNKN ] Undocumented argument [gpointer]
 * Arg:        user_data [UNKN ] Undocumented argument [gpointer]
 * Arg:            data2 [UNKN ] Undocumented argument [gpointer]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_no_op_func(gpointer  data,gpointer  user_data,gpointer  data2);
#define no_op_func Wise2_no_op_func


/* Function:  one_off_HSPscan_scan_query_direct(data,seq,para)
 *
 * Descrip:    Simple word expansion for direct access
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 * Arg:         seq [UNKN ] Undocumented argument [Sequence *]
 * Arg:        para [UNKN ] Undocumented argument [HSPScanInterfacePara *]
 *
 * Return [UNKN ]  Undocumented return value [LinearHSPmanager *]
 *
 */
LinearHSPmanager * Wise2_one_off_HSPscan_scan_query_direct(void * data,Sequence * seq,HSPScanInterfacePara * para);
#define one_off_HSPscan_scan_query_direct Wise2_one_off_HSPscan_scan_query_direct


/* Function:  one_off_HSPscan_scan_query(data,seq,para)
 *
 * Descrip:    Simple word expansion - one off score drop considered
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 * Arg:         seq [UNKN ] Undocumented argument [Sequence *]
 * Arg:        para [UNKN ] Undocumented argument [HSPScanInterfacePara *]
 *
 * Return [UNKN ]  Undocumented return value [LinearHSPmanager *]
 *
 */
LinearHSPmanager * Wise2_one_off_HSPscan_scan_query(void * data,Sequence * seq,HSPScanInterfacePara * para);
#define one_off_HSPscan_scan_query Wise2_one_off_HSPscan_scan_query


/* Function:  simple_HSPScan_scan_query(data,seq,para)
 *
 * Descrip:    simple Scan function, no word expansions
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 * Arg:         seq [UNKN ] Undocumented argument [Sequence *]
 * Arg:        para [UNKN ] Undocumented argument [HSPScanInterfacePara *]
 *
 * Return [UNKN ]  Undocumented return value [LinearHSPmanager *]
 *
 */
LinearHSPmanager * Wise2_simple_HSPScan_scan_query(void * data,Sequence * seq,HSPScanInterfacePara * para);
#define simple_HSPScan_scan_query Wise2_simple_HSPScan_scan_query


/* Function:  simple_HSPScan_free(data)
 *
 * Descrip:    Free function for simple scans
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 *
 */
void Wise2_simple_HSPScan_free(void * data);
#define simple_HSPScan_free Wise2_simple_HSPScan_free


/* Function:  hard_link_HSPScanPara(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [HSPScanPara *]
 *
 * Return [UNKN ]  Undocumented return value [HSPScanPara *]
 *
 */
HSPScanPara * Wise2_hard_link_HSPScanPara(HSPScanPara * obj);
#define hard_link_HSPScanPara Wise2_hard_link_HSPScanPara


/* Function:  HSPScanPara_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [HSPScanPara *]
 *
 */
HSPScanPara * Wise2_HSPScanPara_alloc(void);
#define HSPScanPara_alloc Wise2_HSPScanPara_alloc


/* Function:  free_HSPScanPara(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [HSPScanPara *]
 *
 * Return [UNKN ]  Undocumented return value [HSPScanPara *]
 *
 */
HSPScanPara * Wise2_free_HSPScanPara(HSPScanPara * obj);
#define free_HSPScanPara Wise2_free_HSPScanPara


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
