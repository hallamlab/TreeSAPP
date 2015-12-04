#ifndef DYNAMITEsubseqhashHEADERFILE
#define DYNAMITEsubseqhashHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "sequence.h"
#include "seqlookup.h"
#include "linkedlist_lookpos.h"
#include "glib.h"

#include "wisebase.h"




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  new_ghash_SeqLookupInterface(void)
 *
 * Descrip:    Makes a new Glib based SeqLookup system
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupInterface *]
 *
 */
SeqLookupInterface * Wise2_new_ghash_SeqLookupInterface(void);
#define new_ghash_SeqLookupInterface Wise2_new_ghash_SeqLookupInterface


/* Function:  get_client_subseqhash_ghash(data)
 *
 * Descrip:    provides the interface defiition for get_client
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupClientInterface *]
 *
 */
SeqLookupClientInterface * Wise2_get_client_subseqhash_ghash(void * data);
#define get_client_subseqhash_ghash Wise2_get_client_subseqhash_ghash


/* Function:  free_client_subseqhash_ghash(data)
 *
 * Descrip:    provides the interface definition for free_data on client, which is
 *             actually a no-op in this case
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 *
 */
void Wise2_free_client_subseqhash_ghash(void * data);
#define free_client_subseqhash_ghash Wise2_free_client_subseqhash_ghash


/* Function:  free_subseqhash_ghash(data)
 *
 * Descrip:    Internal function for freeing ghash
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 *
 */
void Wise2_free_subseqhash_ghash(void * data);
#define free_subseqhash_ghash Wise2_free_subseqhash_ghash


/* Function:  remove_subseq_SeqLookupPos(key,value,user_data)
 *
 * Descrip:    Sub call of free
 *
 *
 * Arg:              key [UNKN ] Undocumented argument [gpointer]
 * Arg:            value [UNKN ] Undocumented argument [gpointer]
 * Arg:        user_data [UNKN ] Undocumented argument [gpointer]
 *
 * Return [UNKN ]  Undocumented return value [gboolean]
 *
 */
gboolean Wise2_remove_subseq_SeqLookupPos(gpointer key,gpointer value,gpointer user_data);
#define remove_subseq_SeqLookupPos Wise2_remove_subseq_SeqLookupPos


/* Function:  is_populated_subseqhash_ghash(data,seq_number)
 *
 * Descrip:    Tells whether a position is populated or not
 *
 *
 * Arg:              data [UNKN ] Undocumented argument [void *]
 * Arg:        seq_number [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_is_populated_subseqhash_ghash(void * data, int seq_number);
#define is_populated_subseqhash_ghash Wise2_is_populated_subseqhash_ghash


/* Function:  lookup_subseqhash_ghash(data,seq_number)
 *
 * Descrip:    Retrieves a SeqLookupPos from the hash
 *
 *
 * Arg:              data [UNKN ] Undocumented argument [void *]
 * Arg:        seq_number [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupResultInterface *]
 *
 */
SeqLookupResultInterface * Wise2_lookup_subseqhash_ghash(void * data,int seq_number);
#define lookup_subseqhash_ghash Wise2_lookup_subseqhash_ghash


/* Function:  add_seq_subseqhash_ghash(data,seq,para)
 *
 * Descrip:    Adds a sequence/pos pair to the hash with this
 *             number
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 * Arg:         seq [UNKN ] Undocumented argument [Sequence *]
 * Arg:        para [UNKN ] Undocumented argument [SeqLookupLoadPara *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_seq_subseqhash_ghash(void * data,Sequence * seq,SeqLookupLoadPara * para);
#define add_seq_subseqhash_ghash Wise2_add_seq_subseqhash_ghash


/* Function:  add_direct_number_subseqhash_ghash(data,seq_number,target,pos)
 *
 * Descrip:    Adds a direct no/seq pair to the system
 *             number
 *
 *
 * Arg:              data [UNKN ] Undocumented argument [void *]
 * Arg:        seq_number [UNKN ] Undocumented argument [int]
 * Arg:            target [UNKN ] Undocumented argument [Sequence *]
 * Arg:               pos [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_direct_number_subseqhash_ghash(void * data,int seq_number,Sequence * target,int pos);
#define add_direct_number_subseqhash_ghash Wise2_add_direct_number_subseqhash_ghash


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
