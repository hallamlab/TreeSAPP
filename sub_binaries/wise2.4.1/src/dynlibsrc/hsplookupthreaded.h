#ifndef DYNAMITEhsplookupthreadedHEADERFILE
#define DYNAMITEhsplookupthreadedHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "hsplookupscan.h"
#include "hsphandler.h"
#include "arrayseqlookup.h"


#define MAX_HSP_THREADS 64

typedef struct ordered_pos_holder {
  Sequence * seq;
  int target_pos;
  int query_pos;
  int diagonal;
} OrderedPosHolder;

#define COMP_OPH(a,b) ((a.seq - b.seq) == 0 ? (a.diagonal - b.diagonal) == 0 ? (a.query_pos - b.query_pos) : (a.diagonal - b. diagonal) : (a.seq - b.seq))

#define COMP_OPH_POINTER(a,b) ((a->seq - b->seq) == 0 ? ((a->diagonal - b->diagonal) == 0 ? (a->query_pos - b->query_pos) : (a->diagonal - b->diagonal)) : (a->seq - b->seq))


struct Wise2_HSPLookupThreadHolder {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    HSPScanPara * p;     
    HSPmanager * hspm;   
    Sequence * seq;  
    HSPScanInterfacePara * para;     
    int start;   
    int end;     
    } ;  
/* HSPLookupThreadHolder defined */ 
#ifndef DYNAMITE_DEFINED_HSPLookupThreadHolder
typedef struct Wise2_HSPLookupThreadHolder Wise2_HSPLookupThreadHolder;
#define HSPLookupThreadHolder Wise2_HSPLookupThreadHolder
#define DYNAMITE_DEFINED_HSPLookupThreadHolder
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  new_threaded_HSPScanInterface(sli,mat,drop_off,score_cutoff,threadno)
 *
 * Descrip:    Makes the wrapper structure for
 *             threaded searches
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
HSPScanInterface * Wise2_new_threaded_HSPScanInterface(SeqLookupInterface * sli,CompMat * mat,int drop_off,int score_cutoff,int threadno);
#define new_threaded_HSPScanInterface Wise2_new_threaded_HSPScanInterface


/* Function:  new_ordered_HSPScanInterface(sli,mat,drop_off,score_cutoff)
 *
 * Descrip:    Makes the wrapper structure for ordered searches
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
HSPScanInterface * Wise2_new_ordered_HSPScanInterface(SeqLookupInterface * sli,CompMat * mat,int drop_off,int score_cutoff);
#define new_ordered_HSPScanInterface Wise2_new_ordered_HSPScanInterface


/* Function:  one_off_ordered_HSPscan_scan_query_direct(data,seq,para)
 *
 * Descrip:    function which orders memory access pattern
 *             more sensibly
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 * Arg:         seq [UNKN ] Undocumented argument [Sequence *]
 * Arg:        para [UNKN ] Undocumented argument [HSPScanInterfacePara *]
 *
 * Return [UNKN ]  Undocumented return value [LinearHSPmanager *]
 *
 */
LinearHSPmanager * Wise2_one_off_ordered_HSPscan_scan_query_direct(void * data,Sequence * seq,HSPScanInterfacePara * para);
#define one_off_ordered_HSPscan_scan_query_direct Wise2_one_off_ordered_HSPscan_scan_query_direct


/* Function:  qsort_oph_pointer(pos,left,right)
 *
 * Descrip:    internal quicksort function on oph pointer arrays
 *
 *
 * Arg:          pos [UNKN ] Undocumented argument [OrderedPosHolder **]
 * Arg:         left [UNKN ] Undocumented argument [long int]
 * Arg:        right [UNKN ] Undocumented argument [long int]
 *
 */
void Wise2_qsort_oph_pointer(OrderedPosHolder ** pos,long int left,long int right);
#define qsort_oph_pointer Wise2_qsort_oph_pointer


/* Function:  qsort_oph(oph,left,right)
 *
 * Descrip:    internal quicksort function for ordered access
 *
 *
 * Arg:          oph [UNKN ] Undocumented argument [OrderedPosHolder *]
 * Arg:         left [UNKN ] Undocumented argument [int]
 * Arg:        right [UNKN ] Undocumented argument [int]
 *
 */
void Wise2_qsort_oph(OrderedPosHolder * oph,int left,int right);
#define qsort_oph Wise2_qsort_oph


/* Function:  qsort_seqnumber_array(position,seqnumber,left,right)
 *
 * Descrip:    internal quicksort function for ordered access
 *
 *
 * Arg:         position [UNKN ] Undocumented argument [int *]
 * Arg:        seqnumber [UNKN ] Undocumented argument [int *]
 * Arg:             left [UNKN ] Undocumented argument [int]
 * Arg:            right [UNKN ] Undocumented argument [int]
 *
 */
void Wise2_qsort_seqnumber_array(int * position,int * seqnumber,int left,int right);
#define qsort_seqnumber_array Wise2_qsort_seqnumber_array


/* Function:  one_off_threaded_HSPscan_scan_query_direct(data,seq,para)
 *
 * Descrip:    Main threaded function
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 * Arg:         seq [UNKN ] Undocumented argument [Sequence *]
 * Arg:        para [UNKN ] Undocumented argument [HSPScanInterfacePara *]
 *
 * Return [UNKN ]  Undocumented return value [LinearHSPmanager *]
 *
 */
LinearHSPmanager * Wise2_one_off_threaded_HSPscan_scan_query_direct(void * data,Sequence * seq,HSPScanInterfacePara * para);
#define one_off_threaded_HSPscan_scan_query_direct Wise2_one_off_threaded_HSPscan_scan_query_direct


/* Function:  merge_HSPmanager(from,to)
 *
 * Descrip:    function to merge one HSPmanager into another
 *             one
 *
 *
 * Arg:        from [UNKN ] Undocumented argument [HSPmanager *]
 * Arg:          to [UNKN ] Undocumented argument [HSPmanager *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_merge_HSPmanager(HSPmanager * from,HSPmanager * to);
#define merge_HSPmanager Wise2_merge_HSPmanager


/* Function:  merge_HSPmanager_foreach(key,value,user_data)
 *
 * Descrip:    internal function to merge HSPmanagers
 *
 *
 * Arg:              key [UNKN ] Undocumented argument [gpointer]
 * Arg:            value [UNKN ] Undocumented argument [gpointer]
 * Arg:        user_data [UNKN ] Undocumented argument [gpointer]
 *
 */
void Wise2_merge_HSPmanager_foreach(gpointer key,gpointer value,gpointer user_data);
#define merge_HSPmanager_foreach Wise2_merge_HSPmanager_foreach


/* Function:  hsp_thread_worker(ptr)
 *
 * Descrip:    thread worker
 *
 *
 * Arg:        ptr [UNKN ] Undocumented argument [void *]
 *
 * Return [UNKN ]  Undocumented return value [void *]
 *
 */
void * Wise2_hsp_thread_worker(void * ptr);
#define hsp_thread_worker Wise2_hsp_thread_worker


/* Function:  hard_link_HSPLookupThreadHolder(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [HSPLookupThreadHolder *]
 *
 * Return [UNKN ]  Undocumented return value [HSPLookupThreadHolder *]
 *
 */
HSPLookupThreadHolder * Wise2_hard_link_HSPLookupThreadHolder(HSPLookupThreadHolder * obj);
#define hard_link_HSPLookupThreadHolder Wise2_hard_link_HSPLookupThreadHolder


/* Function:  HSPLookupThreadHolder_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [HSPLookupThreadHolder *]
 *
 */
HSPLookupThreadHolder * Wise2_HSPLookupThreadHolder_alloc(void);
#define HSPLookupThreadHolder_alloc Wise2_HSPLookupThreadHolder_alloc


/* Function:  free_HSPLookupThreadHolder(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [HSPLookupThreadHolder *]
 *
 * Return [UNKN ]  Undocumented return value [HSPLookupThreadHolder *]
 *
 */
HSPLookupThreadHolder * Wise2_free_HSPLookupThreadHolder(HSPLookupThreadHolder * obj);
#define free_HSPLookupThreadHolder Wise2_free_HSPLookupThreadHolder


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
