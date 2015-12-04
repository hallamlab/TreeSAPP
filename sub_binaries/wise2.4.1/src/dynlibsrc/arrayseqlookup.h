#ifndef DYNAMITEarrayseqlookupHEADERFILE
#define DYNAMITEarrayseqlookupHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "seqlookup.h"


#define ARRAYSEQL_BASIC  8
#define ARRAYSEQL_LINEAR 500

struct Wise2_ArraySeqLookup {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    ArraySeqHead ** array;   
    int len;     
    long numb_level;     
    } ;  
/* ArraySeqLookup defined */ 
#ifndef DYNAMITE_DEFINED_ArraySeqLookup
typedef struct Wise2_ArraySeqLookup Wise2_ArraySeqLookup;
#define ArraySeqLookup Wise2_ArraySeqLookup
#define DYNAMITE_DEFINED_ArraySeqLookup
#endif


struct Wise2_ArraySeqHeadResults {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int ipos;    
    ArraySeqHead * head;     
    SeqLookupResultStruct str;   
    } ;  
/* ArraySeqHeadResults defined */ 
#ifndef DYNAMITE_DEFINED_ArraySeqHeadResults
typedef struct Wise2_ArraySeqHeadResults Wise2_ArraySeqHeadResults;
#define ArraySeqHeadResults Wise2_ArraySeqHeadResults
#define DYNAMITE_DEFINED_ArraySeqHeadResults
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  new_ArraySeq_SeqLookupInterface(len,numb_level)
 *
 * Descrip:    Exported function - makes a new seqlookupinterface in array mode
 *
 *
 * Arg:               len [UNKN ] Undocumented argument [int]
 * Arg:        numb_level [UNKN ] Undocumented argument [long]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupInterface *]
 *
 */
SeqLookupInterface * Wise2_new_ArraySeq_SeqLookupInterface(int len,long numb_level);
#define new_ArraySeq_SeqLookupInterface Wise2_new_ArraySeq_SeqLookupInterface


/* Function:  print_array_occuypancy_ArraySeq(asl,ofp)
 *
 * Descrip:    Prints out summary statistcis to a file
 *
 *
 * Arg:        asl [UNKN ] Undocumented argument [ArraySeqLookup *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_print_array_occuypancy_ArraySeq(ArraySeqLookup * asl,FILE * ofp);
#define print_array_occuypancy_ArraySeq Wise2_print_array_occuypancy_ArraySeq


/* Function:  get_client_arraylookup(*data)
 *
 * Descrip:    Builds a new client from Array
 *
 *
 * Arg:        *data [UNKN ] Undocumented argument [void]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupClientInterface *]
 *
 */
SeqLookupClientInterface * Wise2_get_client_arraylookup(void *data);
#define get_client_arraylookup Wise2_get_client_arraylookup


/* Function:  next_arrayhead_search_results(data,prev)
 *
 * Descrip:    Internal function for results interface
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 * Arg:        prev [UNKN ] Undocumented argument [SeqLookupResultStruct *]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupResultStruct *]
 *
 */
SeqLookupResultStruct * Wise2_next_arrayhead_search_results(void * data,SeqLookupResultStruct * prev);
#define next_arrayhead_search_results Wise2_next_arrayhead_search_results


/* Function:  is_more_arrayhead_search_results(data)
 *
 * Descrip:    Internal function for results interface
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_is_more_arrayhead_search_results(void * data);
#define is_more_arrayhead_search_results Wise2_is_more_arrayhead_search_results


/* Function:  free_arrayhead_results(data)
 *
 * Descrip:    Internal function for results interface, which is a no-op
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 *
 */
void Wise2_free_arrayhead_results(void * data);
#define free_arrayhead_results Wise2_free_arrayhead_results


/* Function:  free_array_client(data)
 *
 * Descrip:    Internal function for client interface, which frees client specific memory
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 *
 */
void Wise2_free_array_client(void * data);
#define free_array_client Wise2_free_array_client


/* Function:  new_ArraySeqLookup(len,numb_level)
 *
 * Descrip:    makes a new ArraySeqLookup taking up
 *             to len positions
 *
 *
 * Arg:               len [UNKN ] Undocumented argument [int]
 * Arg:        numb_level [UNKN ] Undocumented argument [long]
 *
 * Return [UNKN ]  Undocumented return value [ArraySeqLookup *]
 *
 */
ArraySeqLookup * Wise2_new_ArraySeqLookup(int len,long numb_level);
#define new_ArraySeqLookup Wise2_new_ArraySeqLookup


/* Function:  is_populated_array_client(data,seq_number)
 *
 * Descrip:    tells whether this is populated or not
 *
 *
 * Arg:              data [UNKN ] Undocumented argument [void *]
 * Arg:        seq_number [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_is_populated_array_client(void * data, int seq_number);
#define is_populated_array_client Wise2_is_populated_array_client


/* Function:  lookup_array_client(data,seq_number)
 *
 * Descrip:    Retrieves a SeqLookup position 
 *
 *
 * Arg:              data [UNKN ] Undocumented argument [void *]
 * Arg:        seq_number [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupResultInterface *]
 *
 */
SeqLookupResultInterface * Wise2_lookup_array_client(void * data, int seq_number);
#define lookup_array_client Wise2_lookup_array_client


/* Function:  arrayhead_direct_lookup(data,seq_number)
 *
 * Descrip:    For array optimised lookup hash, provides direct memory access
 *
 *
 * Arg:              data [UNKN ] Undocumented argument [void *]
 * Arg:        seq_number [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [ArraySeqHead *]
 *
 */
ArraySeqHead * Wise2_arrayhead_direct_lookup(void * data,int seq_number);
#define arrayhead_direct_lookup Wise2_arrayhead_direct_lookup


/* Function:  add_seq_arraylookup(data,seq,para)
 *
 * Descrip:    Adds a sequence/pos pair to the hash
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 * Arg:         seq [UNKN ] Undocumented argument [Sequence *]
 * Arg:        para [UNKN ] Undocumented argument [SeqLookupLoadPara *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_seq_arraylookup(void * data,Sequence * seq,SeqLookupLoadPara * para);
#define add_seq_arraylookup Wise2_add_seq_arraylookup


/* Function:  free_data_arraylookup(data)
 *
 * Descrip:    Frees data
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 *
 */
void Wise2_free_data_arraylookup(void * data);
#define free_data_arraylookup Wise2_free_data_arraylookup


/* Function:  add_ArraySeqHead(h,seq,pos,numb_level)
 *
 * Descrip:    Adds a sequence/pos pair to an ArrayHead
 *
 *
 * Arg:                 h [UNKN ] Undocumented argument [ArraySeqHead *]
 * Arg:               seq [UNKN ] Undocumented argument [Sequence *]
 * Arg:               pos [UNKN ] Undocumented argument [int]
 * Arg:        numb_level [UNKN ] Undocumented argument [long]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_ArraySeqHead(ArraySeqHead * h,Sequence * seq,int pos,long numb_level);
#define add_ArraySeqHead Wise2_add_ArraySeqHead


/* Function:  new_ArraySeqHead(void)
 *
 * Descrip:    Builds a new ArraySeqHead structure
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ArraySeqHead *]
 *
 */
ArraySeqHead * Wise2_new_ArraySeqHead(void);
#define new_ArraySeqHead Wise2_new_ArraySeqHead


/* Function:  hard_link_ArraySeqLookup(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ArraySeqLookup *]
 *
 * Return [UNKN ]  Undocumented return value [ArraySeqLookup *]
 *
 */
ArraySeqLookup * Wise2_hard_link_ArraySeqLookup(ArraySeqLookup * obj);
#define hard_link_ArraySeqLookup Wise2_hard_link_ArraySeqLookup


/* Function:  ArraySeqLookup_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ArraySeqLookup *]
 *
 */
ArraySeqLookup * Wise2_ArraySeqLookup_alloc(void);
#define ArraySeqLookup_alloc Wise2_ArraySeqLookup_alloc


/* Function:  free_ArraySeqLookup(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ArraySeqLookup *]
 *
 * Return [UNKN ]  Undocumented return value [ArraySeqLookup *]
 *
 */
ArraySeqLookup * Wise2_free_ArraySeqLookup(ArraySeqLookup * obj);
#define free_ArraySeqLookup Wise2_free_ArraySeqLookup


/* Function:  hard_link_ArraySeqHeadResults(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ArraySeqHeadResults *]
 *
 * Return [UNKN ]  Undocumented return value [ArraySeqHeadResults *]
 *
 */
ArraySeqHeadResults * Wise2_hard_link_ArraySeqHeadResults(ArraySeqHeadResults * obj);
#define hard_link_ArraySeqHeadResults Wise2_hard_link_ArraySeqHeadResults


/* Function:  ArraySeqHeadResults_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ArraySeqHeadResults *]
 *
 */
ArraySeqHeadResults * Wise2_ArraySeqHeadResults_alloc(void);
#define ArraySeqHeadResults_alloc Wise2_ArraySeqHeadResults_alloc


/* Function:  free_ArraySeqHeadResults(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ArraySeqHeadResults *]
 *
 * Return [UNKN ]  Undocumented return value [ArraySeqHeadResults *]
 *
 */
ArraySeqHeadResults * Wise2_free_ArraySeqHeadResults(ArraySeqHeadResults * obj);
#define free_ArraySeqHeadResults Wise2_free_ArraySeqHeadResults


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
