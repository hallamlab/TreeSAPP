#ifndef DYNAMITEcompressed_protein_indexHEADERFILE
#define DYNAMITEcompressed_protein_indexHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "seqlookup.h"
#include "kmer_index_interface.h"
#include "kmer_direct.h"
#include "singleseqspace.h"


typedef struct CompressedProteinHead {
  int * positions;
  int maxlen;
} CompressedProteinHead;

#define COMPRESSED_INITIAL_SIZE 50
#define COMPRESSED_EXPAND       50


#define COMPRESSED_PROTEIN_CLIENT_INITIAL 16

struct Wise2_CompressedProteinIndex {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    KmerIndexInterface * kii;    
    SinglePosSpace * sps;    
    } ;  
/* CompressedProteinIndex defined */ 
#ifndef DYNAMITE_DEFINED_CompressedProteinIndex
typedef struct Wise2_CompressedProteinIndex Wise2_CompressedProteinIndex;
#define CompressedProteinIndex Wise2_CompressedProteinIndex
#define DYNAMITE_DEFINED_CompressedProteinIndex
#endif


struct Wise2_CompressedProteinClient {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    SeqLookupResultStruct * res;    /*  actually we'll free this on demand */ 
    int maxlen;  
    int current;     
    int point;   
    CompressedProteinIndex * index;  
    } ;  
/* CompressedProteinClient defined */ 
#ifndef DYNAMITE_DEFINED_CompressedProteinClient
typedef struct Wise2_CompressedProteinClient Wise2_CompressedProteinClient;
#define CompressedProteinClient Wise2_CompressedProteinClient
#define DYNAMITE_DEFINED_CompressedProteinClient
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_CompressedProteinIndex(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [CompressedProteinIndex *]
 *
 * Return [UNKN ]  Undocumented return value [CompressedProteinIndex *]
 *
 */
CompressedProteinIndex * Wise2_hard_link_CompressedProteinIndex(CompressedProteinIndex * obj);
#define hard_link_CompressedProteinIndex Wise2_hard_link_CompressedProteinIndex


/* Function:  CompressedProteinIndex_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [CompressedProteinIndex *]
 *
 */
CompressedProteinIndex * Wise2_CompressedProteinIndex_alloc(void);
#define CompressedProteinIndex_alloc Wise2_CompressedProteinIndex_alloc


/* Function:  free_CompressedProteinIndex(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [CompressedProteinIndex *]
 *
 * Return [UNKN ]  Undocumented return value [CompressedProteinIndex *]
 *
 */
CompressedProteinIndex * Wise2_free_CompressedProteinIndex(CompressedProteinIndex * obj);
#define free_CompressedProteinIndex Wise2_free_CompressedProteinIndex


/* Function:  hard_link_CompressedProteinClient(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [CompressedProteinClient *]
 *
 * Return [UNKN ]  Undocumented return value [CompressedProteinClient *]
 *
 */
CompressedProteinClient * Wise2_hard_link_CompressedProteinClient(CompressedProteinClient * obj);
#define hard_link_CompressedProteinClient Wise2_hard_link_CompressedProteinClient


/* Function:  CompressedProteinClient_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [CompressedProteinClient *]
 *
 */
CompressedProteinClient * Wise2_CompressedProteinClient_alloc(void);
#define CompressedProteinClient_alloc Wise2_CompressedProteinClient_alloc


/* Function:  free_CompressedProteinClient(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [CompressedProteinClient *]
 *
 * Return [UNKN ]  Undocumented return value [CompressedProteinClient *]
 *
 */
CompressedProteinClient * Wise2_free_CompressedProteinClient(CompressedProteinClient * obj);
#define free_CompressedProteinClient Wise2_free_CompressedProteinClient


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
SeqLookupResultStruct * Wise2_get_next_CompressedProteinClient(void * data,SeqLookupResultStruct * prev);
#define get_next_CompressedProteinClient Wise2_get_next_CompressedProteinClient
boolean Wise2_is_more_CompressedProteinClient(void * data);
#define is_more_CompressedProteinClient Wise2_is_more_CompressedProteinClient
void Wise2_free_SeqLookupResultInterface_CompressedProteinClient(void * data);
#define free_SeqLookupResultInterface_CompressedProteinClient Wise2_free_SeqLookupResultInterface_CompressedProteinClient
SeqLookupResultInterface * Wise2_lookup_CompressedProteinClient(void * data,int seq_number);
#define lookup_CompressedProteinClient Wise2_lookup_CompressedProteinClient
boolean Wise2_is_populated_CompressedProteinClient(void * data,int seq_number);
#define is_populated_CompressedProteinClient Wise2_is_populated_CompressedProteinClient
void Wise2_free_Client_CompressedProteinClient(void * data);
#define free_Client_CompressedProteinClient Wise2_free_Client_CompressedProteinClient
SeqLookupInterface * Wise2_new_direct_CompressedProteinLookup(void);
#define new_direct_CompressedProteinLookup Wise2_new_direct_CompressedProteinLookup
SeqLookupInterface * Wise2_new_CompressedProteinLookup(KmerIndexInterface * kii);
#define new_CompressedProteinLookup Wise2_new_CompressedProteinLookup
SeqLookupClientInterface * Wise2_get_client_CompressedProteinIndex(void * data);
#define get_client_CompressedProteinIndex Wise2_get_client_CompressedProteinIndex
CompressedProteinClient * Wise2_new_CompressedProteinClient(void) ;
#define new_CompressedProteinClient Wise2_new_CompressedProteinClient
boolean Wise2_add_direct_number_CompressedProteinIndex(void * data,int seq_number,Sequence * target, int pos) ;
#define add_direct_number_CompressedProteinIndex Wise2_add_direct_number_CompressedProteinIndex
ArraySeqHead * Wise2_lookup_array_head_CompressedProteinIndex(void * data,int seq_number);
#define lookup_array_head_CompressedProteinIndex Wise2_lookup_array_head_CompressedProteinIndex
boolean Wise2_add_seq_CompressedProteinIndex(void * data,Sequence * seq,SeqLookupLoadPara * para);
#define add_seq_CompressedProteinIndex Wise2_add_seq_CompressedProteinIndex
void Wise2_free_data_CompressedProteinIndex(void * d);
#define free_data_CompressedProteinIndex Wise2_free_data_CompressedProteinIndex
boolean Wise2_add_pos_CompressedProteinHead(CompressedProteinHead * c,int pos);
#define add_pos_CompressedProteinHead Wise2_add_pos_CompressedProteinHead
CompressedProteinHead * Wise2_new_CompressedProteinHead(void);
#define new_CompressedProteinHead Wise2_new_CompressedProteinHead


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/

#ifdef _cplusplus
}
#endif

#endif
