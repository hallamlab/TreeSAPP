#ifndef DYNAMITEdnaindexconsHEADERFILE
#define DYNAMITEdnaindexconsHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "subseqhash.h"
#include "hsp.h"
#include "hsphandler.h"
#include "dnamatrix.h"

typedef enum {
  DnaIndexConstructor_subseq = 438,
} DnaIndexConstructorType;

struct Wise2_DnaIndexConstructor {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int type;    
    int index_word_length;   
    Probability match_prob;  
    int drop_off;    
    } ;  
/* DnaIndexConstructor defined */ 
#ifndef DYNAMITE_DEFINED_DnaIndexConstructor
typedef struct Wise2_DnaIndexConstructor Wise2_DnaIndexConstructor;
#define DnaIndexConstructor Wise2_DnaIndexConstructor
#define DYNAMITE_DEFINED_DnaIndexConstructor
#endif


struct Wise2_DnaIndex {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int index_word_length;   
    SeqLookupInterface * sli;    
    CompMat * cm;    
    int drop_off;    
    } ;  
/* DnaIndex defined */ 
#ifndef DYNAMITE_DEFINED_DnaIndex
typedef struct Wise2_DnaIndex Wise2_DnaIndex;
#define DnaIndex Wise2_DnaIndex
#define DYNAMITE_DEFINED_DnaIndex
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  show_help_DnaIndexConstructor(ofp)
 *
 * Descrip:    Shows help for a DnaIndexConstructor
 *
 *
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_show_help_DnaIndexConstructor(FILE * ofp);
#define show_help_DnaIndexConstructor Wise2_show_help_DnaIndexConstructor


/* Function:  new_DnaIndexConstructor(argc,argv)
 *
 * Descrip:    Builds new DnaIndexConstructor off Command line
 *
 *
 * Arg:        argc [UNKN ] Undocumented argument [int *]
 * Arg:        argv [UNKN ] Undocumented argument [char **]
 *
 * Return [UNKN ]  Undocumented return value [DnaIndexConstructor *]
 *
 */
DnaIndexConstructor * Wise2_new_DnaIndexConstructor(int * argc,char ** argv);
#define new_DnaIndexConstructor Wise2_new_DnaIndexConstructor


/* Function:  DnaIndex_from_DnaIndexConstructor(dic)
 *
 * Descrip:    Makes a DnaIndex from a DnaIndexConstructore
 *
 *
 * Arg:        dic [UNKN ] Undocumented argument [DnaIndexConstructor *]
 *
 * Return [UNKN ]  Undocumented return value [DnaIndex *]
 *
 */
DnaIndex * Wise2_DnaIndex_from_DnaIndexConstructor(DnaIndexConstructor * dic);
#define DnaIndex_from_DnaIndexConstructor Wise2_DnaIndex_from_DnaIndexConstructor


/* Function:  LinearHSPManager_scan_DnaIndex(di,query)
 *
 * Descrip:    provides a LinearManager for a DNA sequence 
 *
 *
 * Arg:           di [UNKN ] Undocumented argument [DnaIndex *]
 * Arg:        query [UNKN ] Undocumented argument [Sequence *]
 *
 * Return [UNKN ]  Undocumented return value [LinearHSPmanager *]
 *
 */
LinearHSPmanager * Wise2_LinearHSPManager_scan_DnaIndex(DnaIndex * di,Sequence * query);
#define LinearHSPManager_scan_DnaIndex Wise2_LinearHSPManager_scan_DnaIndex


/* Function:  HSPmanager_scan_DnaIndex(di,seq)
 *
 * Descrip:    Provides a HSPmanager from a scan
 *
 *
 * Arg:         di [UNKN ] Undocumented argument [DnaIndex *]
 * Arg:        seq [UNKN ] Undocumented argument [Sequence *]
 *
 * Return [UNKN ]  Undocumented return value [HSPmanager *]
 *
 */
HSPmanager * Wise2_HSPmanager_scan_DnaIndex(DnaIndex * di,Sequence * seq);
#define HSPmanager_scan_DnaIndex Wise2_HSPmanager_scan_DnaIndex


/* Function:  process_dna_HSP(set,query,query_pos,tseq,res_struct,mat)
 *
 * Descrip:    processes DNA based HSPs
 *
 *
 * Arg:               set [UNKN ] Undocumented argument [HSPset *]
 * Arg:             query [UNKN ] Undocumented argument [Sequence *]
 * Arg:         query_pos [UNKN ] Undocumented argument [int]
 * Arg:              tseq [UNKN ] Undocumented argument [Sequence *]
 * Arg:        res_struct [UNKN ] Undocumented argument [SeqLookupResultStruct *]
 * Arg:               mat [UNKN ] Undocumented argument [CompMat *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_process_dna_HSP(HSPset * set,Sequence * query,int query_pos,Sequence * tseq,SeqLookupResultStruct * res_struct,CompMat * mat);
#define process_dna_HSP Wise2_process_dna_HSP


/* Function:  load_Sequence_DnaIndex(di,seq,sllp)
 *
 * Descrip:    Loads a sequence into a DnaIndex
 *
 *
 * Arg:          di [UNKN ] Undocumented argument [DnaIndex *]
 * Arg:         seq [UNKN ] Undocumented argument [Sequence *]
 * Arg:        sllp [UNKN ] Undocumented argument [SeqLookupLoadPara *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_load_Sequence_DnaIndex(DnaIndex * di,Sequence * seq,SeqLookupLoadPara * sllp);
#define load_Sequence_DnaIndex Wise2_load_Sequence_DnaIndex


/* Function:  seq_number_dna_Nmer_noN(seq,index_length)
 *
 * Descrip:    General DNA index number generation
 *
 *
 * Arg:                 seq [UNKN ] Undocumented argument [char *]
 * Arg:        index_length [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_seq_number_dna_Nmer_noN(char * seq,int index_length);
#define seq_number_dna_Nmer_noN Wise2_seq_number_dna_Nmer_noN


/* Function:  hard_link_DnaIndexConstructor(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DnaIndexConstructor *]
 *
 * Return [UNKN ]  Undocumented return value [DnaIndexConstructor *]
 *
 */
DnaIndexConstructor * Wise2_hard_link_DnaIndexConstructor(DnaIndexConstructor * obj);
#define hard_link_DnaIndexConstructor Wise2_hard_link_DnaIndexConstructor


/* Function:  DnaIndexConstructor_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaIndexConstructor *]
 *
 */
DnaIndexConstructor * Wise2_DnaIndexConstructor_alloc(void);
#define DnaIndexConstructor_alloc Wise2_DnaIndexConstructor_alloc


/* Function:  free_DnaIndexConstructor(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DnaIndexConstructor *]
 *
 * Return [UNKN ]  Undocumented return value [DnaIndexConstructor *]
 *
 */
DnaIndexConstructor * Wise2_free_DnaIndexConstructor(DnaIndexConstructor * obj);
#define free_DnaIndexConstructor Wise2_free_DnaIndexConstructor


/* Function:  hard_link_DnaIndex(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DnaIndex *]
 *
 * Return [UNKN ]  Undocumented return value [DnaIndex *]
 *
 */
DnaIndex * Wise2_hard_link_DnaIndex(DnaIndex * obj);
#define hard_link_DnaIndex Wise2_hard_link_DnaIndex


/* Function:  DnaIndex_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaIndex *]
 *
 */
DnaIndex * Wise2_DnaIndex_alloc(void);
#define DnaIndex_alloc Wise2_DnaIndex_alloc


/* Function:  free_DnaIndex(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DnaIndex *]
 *
 * Return [UNKN ]  Undocumented return value [DnaIndex *]
 *
 */
DnaIndex * Wise2_free_DnaIndex(DnaIndex * obj);
#define free_DnaIndex Wise2_free_DnaIndex


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
