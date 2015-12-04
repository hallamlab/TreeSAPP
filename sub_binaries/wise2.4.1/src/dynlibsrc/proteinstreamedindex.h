#ifndef DYNAMITEproteinstreamedindexHEADERFILE
#define DYNAMITEproteinstreamedindexHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "sequence.h"
#include "seqlookup.h"
#include "genericindexresult.h"

typedef struct ProteinStreamedIndexWayPost {
  Sequence * seq;
  int pos;
} ProteinStreamedIndexWayPost;

typedef struct ProteinStreamedIndexPos {
  int count;
  struct ProteinStreamedIndexPos *  first_prev;
  struct ProteinStreamedIndexPos ** prev;
  int prev_len;
  struct ProteinStreamedIndexWayPost * post;
  int post_len;
  int post_maxlen;
} ProteinStreamedIndexPos;

#define PROTEINSTREAMEDINDEX_PREV_START 4
#define PROTEINSTREAMEDINDEX_WAY_START  2

typedef struct ProteinStreamedIndex {
  ProteinStreamedIndexPos ** index;
  int len;
  int waypost_depth;
} ProteinStreamedIndex;



    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  new_ProteinStreamedIndex_SeqLookupInterface(waypost)
 *
 * Descrip:    Provides a SeqLookupInterface, the common runtime plug-in for indexers
 *             using a ProteinStreamedIndex 
 *
 *
 * Arg:        waypost [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupInterface *]
 *
 */
SeqLookupInterface * Wise2_new_ProteinStreamedIndex_SeqLookupInterface(int waypost);
#define new_ProteinStreamedIndex_SeqLookupInterface Wise2_new_ProteinStreamedIndex_SeqLookupInterface


/* Function:  get_client_interface_ProteinStreamedIndex(data)
 *
 * Descrip:    gets client interface
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void*]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupClientInterface *]
 *
 */
SeqLookupClientInterface * Wise2_get_client_interface_ProteinStreamedIndex(void* data);
#define get_client_interface_ProteinStreamedIndex Wise2_get_client_interface_ProteinStreamedIndex


/* Function:  lookup_interface_ProteinStreamedClient(data,seq_number)
 *
 * Descrip:    For lookup interface, provides a result
 *
 *
 * Arg:              data [UNKN ] Undocumented argument [void *]
 * Arg:        seq_number [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupResultInterface *]
 *
 */
SeqLookupResultInterface * Wise2_lookup_interface_ProteinStreamedClient(void * data,int seq_number);
#define lookup_interface_ProteinStreamedClient Wise2_lookup_interface_ProteinStreamedClient


/* Function:  is_populated_interface_ProteinStreamedClient(data,seq_number)
 *
 * Descrip:    populated function for interface
 *
 *
 * Arg:              data [UNKN ] Undocumented argument [void *]
 * Arg:        seq_number [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_is_populated_interface_ProteinStreamedClient(void * data,int seq_number);
#define is_populated_interface_ProteinStreamedClient Wise2_is_populated_interface_ProteinStreamedClient


/* Function:  add_seq_interface_ProteinStreamedIndex(data,seq,para)
 *
 * Descrip:    add function for interface
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 * Arg:         seq [UNKN ] Undocumented argument [Sequence *]
 * Arg:        para [UNKN ] Undocumented argument [SeqLookupLoadPara *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_seq_interface_ProteinStreamedIndex(void * data,Sequence * seq,SeqLookupLoadPara * para);
#define add_seq_interface_ProteinStreamedIndex Wise2_add_seq_interface_ProteinStreamedIndex


/* Function:  free_interface_ProteinStreamedIndex(data)
 *
 * Descrip:    for interface, frees index
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 *
 */
void Wise2_free_interface_ProteinStreamedIndex(void * data);
#define free_interface_ProteinStreamedIndex Wise2_free_interface_ProteinStreamedIndex


/* Function:  free_interface_ProteinStreamedClient(data)
 *
 * Descrip:    Frees the client data
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 *
 */
void Wise2_free_interface_ProteinStreamedClient(void * data);
#define free_interface_ProteinStreamedClient Wise2_free_interface_ProteinStreamedClient


/* Function:  lookup_ProteinStreamedIndex(in,seqno,out)
 *
 * Descrip:    Traverses index for a particular number to return hits
 *
 *
 * Arg:           in [UNKN ] Undocumented argument [ProteinStreamedIndex *]
 * Arg:        seqno [UNKN ] Undocumented argument [int]
 * Arg:          out [UNKN ] Undocumented argument [GenericIndexResult *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_lookup_ProteinStreamedIndex(ProteinStreamedIndex * in,int seqno,GenericIndexResult * out);
#define lookup_ProteinStreamedIndex Wise2_lookup_ProteinStreamedIndex


/* Function:  add_lookup_GenericIndexResult_from_Pos(p,psir,backtrace,waylength)
 *
 * Descrip:    Adds additional sequences to index result from a particular position
 *
 *
 * Arg:                p [UNKN ] Undocumented argument [ProteinStreamedIndexPos *]
 * Arg:             psir [UNKN ] Undocumented argument [GenericIndexResult *]
 * Arg:        backtrace [UNKN ] Undocumented argument [int]
 * Arg:        waylength [UNKN ] Undocumented argument [int]
 *
 */
void Wise2_add_lookup_GenericIndexResult_from_Pos(ProteinStreamedIndexPos * p,GenericIndexResult * psir,int backtrace,int waylength);
#define add_lookup_GenericIndexResult_from_Pos Wise2_add_lookup_GenericIndexResult_from_Pos


/* Function:  add_Sequence_ProteinStreamedIndex(in,seq,para)
 *
 * Descrip:    Adds in a sequence into a ProteinStreamedIndex
 *
 *
 * Arg:          in [UNKN ] Undocumented argument [ProteinStreamedIndex *]
 * Arg:         seq [UNKN ] Undocumented argument [Sequence *]
 * Arg:        para [UNKN ] Undocumented argument [SeqLookupLoadPara *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_Sequence_ProteinStreamedIndex(ProteinStreamedIndex * in,Sequence * seq,SeqLookupLoadPara * para);
#define add_Sequence_ProteinStreamedIndex Wise2_add_Sequence_ProteinStreamedIndex


/* Function:  dump_ProteinStreamedIndex(in,ofp)
 *
 * Descrip:    dumps in a silly format the ProteinStreamedIndex; for debugging
 *
 *
 * Arg:         in [UNKN ] Undocumented argument [ProteinStreamedIndex *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_dump_ProteinStreamedIndex(ProteinStreamedIndex * in,FILE * ofp);
#define dump_ProteinStreamedIndex Wise2_dump_ProteinStreamedIndex


/* Function:  new_ProteinStreamedIndex(waypost_depth)
 *
 * Descrip:    Builds a new ProteinStreamedIndex
 *
 *
 * Arg:        waypost_depth [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [ProteinStreamedIndex *]
 *
 */
ProteinStreamedIndex * Wise2_new_ProteinStreamedIndex(int waypost_depth);
#define new_ProteinStreamedIndex Wise2_new_ProteinStreamedIndex


/* Function:  free_ProteinStreamedIndex(in)
 *
 * Descrip:    Release protein index
 *
 *
 * Arg:        in [UNKN ] Undocumented argument [ProteinStreamedIndex *]
 *
 */
void Wise2_free_ProteinStreamedIndex(ProteinStreamedIndex * in);
#define free_ProteinStreamedIndex Wise2_free_ProteinStreamedIndex


/* Function:  add_pos_ProteinStreamedIndexPos(pos,prev)
 *
 * Descrip:    Adds a position, advancing array if needed
 *
 *
 * Arg:         pos [UNKN ] Undocumented argument [ProteinStreamedIndexPos *]
 * Arg:        prev [UNKN ] Undocumented argument [ProteinStreamedIndexPos *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_pos_ProteinStreamedIndexPos(ProteinStreamedIndexPos * pos,ProteinStreamedIndexPos * prev);
#define add_pos_ProteinStreamedIndexPos Wise2_add_pos_ProteinStreamedIndexPos


/* Function:  add_waypost_ProteinStreamedIndexPos(pos,seq,seqpos)
 *
 * Descrip:    Adds a way post, advancing array if needed
 *
 *
 * Arg:           pos [UNKN ] Undocumented argument [ProteinStreamedIndexPos *]
 * Arg:           seq [UNKN ] Undocumented argument [Sequence *]
 * Arg:        seqpos [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_waypost_ProteinStreamedIndexPos(ProteinStreamedIndexPos * pos,Sequence * seq,int seqpos);
#define add_waypost_ProteinStreamedIndexPos Wise2_add_waypost_ProteinStreamedIndexPos


/* Function:  new_ProteinStreamedIndexPos(void)
 *
 * Descrip:    Creates a new ProteinStreamedIndexPos
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ProteinStreamedIndexPos *]
 *
 */
ProteinStreamedIndexPos * Wise2_new_ProteinStreamedIndexPos(void);
#define new_ProteinStreamedIndexPos Wise2_new_ProteinStreamedIndexPos


/* Function:  free_ProteinStreamedIndexPos(pos)
 *
 * Descrip:    frees a ProteinStreamedIndexPos
 *
 *
 * Arg:        pos [UNKN ] Undocumented argument [ProteinStreamedIndexPos *]
 *
 */
void Wise2_free_ProteinStreamedIndexPos(ProteinStreamedIndexPos * pos);
#define free_ProteinStreamedIndexPos Wise2_free_ProteinStreamedIndexPos


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
