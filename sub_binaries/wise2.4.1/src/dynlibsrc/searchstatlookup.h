#ifndef DYNAMITEsearchstatlookupHEADERFILE
#define DYNAMITEsearchstatlookupHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "searchstatinterface.h"
#include "probability.h"
#include "histogram.h"

struct Wise2_EVDLookup {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    double mu;   
    double lambda;   
    } ;  
/* EVDLookup defined */ 
#ifndef DYNAMITE_DEFINED_EVDLookup
typedef struct Wise2_EVDLookup Wise2_EVDLookup;
#define EVDLookup Wise2_EVDLookup
#define DYNAMITE_DEFINED_EVDLookup
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  bit_halfbit_lookup_ssi(data,query_len,target_len,raw_score)
 *
 * Descrip:    Internal function for bit conversion
 *
 *
 * Arg:              data [UNKN ] Undocumented argument [void *]
 * Arg:         query_len [UNKN ] Undocumented argument [int]
 * Arg:        target_len [UNKN ] Undocumented argument [int]
 * Arg:         raw_score [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [double]
 *
 */
double Wise2_bit_halfbit_lookup_ssi(void * data,int query_len,int target_len,int raw_score);
#define bit_halfbit_lookup_ssi Wise2_bit_halfbit_lookup_ssi


/* Function:  evalue_halfbit_lookup_ssi(data,a,b,raw_score,database_size)
 *
 * Descrip:    Internal function for evalue conversion. uses externally
 *             defined parameters for evd estimation
 *
 *
 * Arg:                 data [UNKN ] Undocumented argument [void *]
 * Arg:                    a [UNKN ] Undocumented argument [Sequence *]
 * Arg:                    b [UNKN ] Undocumented argument [Sequence *]
 * Arg:            raw_score [UNKN ] Undocumented argument [int]
 * Arg:        database_size [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [double]
 *
 */
double Wise2_evalue_halfbit_lookup_ssi(void * data,Sequence * a,Sequence * b,int raw_score,int database_size);
#define evalue_halfbit_lookup_ssi Wise2_evalue_halfbit_lookup_ssi


/* Function:  free_evdlookup_void (data)
 *
 * Descrip:    Internal function for free...
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 *
 */
void Wise2_free_evdlookup_void (void * data);
#define free_evdlookup_void  Wise2_free_evdlookup_void 


/* Function:  new_lookup_SearchStatInterface(mu,lambda)
 *
 * Descrip:    Builds an external lookup statistics package
 *
 *
 * Arg:            mu [UNKN ] Undocumented argument [double]
 * Arg:        lambda [UNKN ] Undocumented argument [double]
 *
 * Return [UNKN ]  Undocumented return value [SearchStatInterface *]
 *
 */
SearchStatInterface * Wise2_new_lookup_SearchStatInterface(double mu,double lambda);
#define new_lookup_SearchStatInterface Wise2_new_lookup_SearchStatInterface


/* Function:  hard_link_EVDLookup(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [EVDLookup *]
 *
 * Return [UNKN ]  Undocumented return value [EVDLookup *]
 *
 */
EVDLookup * Wise2_hard_link_EVDLookup(EVDLookup * obj);
#define hard_link_EVDLookup Wise2_hard_link_EVDLookup


/* Function:  EVDLookup_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [EVDLookup *]
 *
 */
EVDLookup * Wise2_EVDLookup_alloc(void);
#define EVDLookup_alloc Wise2_EVDLookup_alloc


/* Function:  free_EVDLookup(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [EVDLookup *]
 *
 * Return [UNKN ]  Undocumented return value [EVDLookup *]
 *
 */
EVDLookup * Wise2_free_EVDLookup(EVDLookup * obj);
#define free_EVDLookup Wise2_free_EVDLookup


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
