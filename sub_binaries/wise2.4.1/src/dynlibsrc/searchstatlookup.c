#ifdef _cplusplus
extern "C" {
#endif
#include "searchstatlookup.h"


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
# line 24 "searchstatlookup.dy"
double bit_halfbit_lookup_ssi(void * data,int query_len,int target_len,int raw_score)
{
  return raw_score / 2.0;
}

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
# line 33 "searchstatlookup.dy"
double evalue_halfbit_lookup_ssi(void * data,Sequence * a,Sequence * b,int raw_score,int database_size)
{
  EVDLookup * l = (EVDLookup *)data;
  return database_size * ExtremeValueP((float)raw_score,l->mu,l->lambda);
}

/* Function:  free_evdlookup_void (data)
 *
 * Descrip:    Internal function for free...
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 *
 */
# line 42 "searchstatlookup.dy"
void free_evdlookup_void (void * data)
{
  EVDLookup * l = (EVDLookup *)data;
  free_EVDLookup(l);
}

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
# line 51 "searchstatlookup.dy"
SearchStatInterface * new_lookup_SearchStatInterface(double mu,double lambda)
{
  SearchStatInterface * ssi;
  EVDLookup * el;

  ssi = SearchStatInterface_alloc();
  el = EVDLookup_alloc();
  el->mu = mu;
  el->lambda = lambda;
  
  ssi->data = (void *) el;

  ssi->calc_evalue = evalue_halfbit_lookup_ssi;
  ssi->calc_bits   = bit_halfbit_lookup_ssi;
  ssi->free_data   = free_evdlookup_void;

  return ssi;
}
# line 84 "searchstatlookup.c"
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
EVDLookup * hard_link_EVDLookup(EVDLookup * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a EVDLookup object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  EVDLookup_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [EVDLookup *]
 *
 */
EVDLookup * EVDLookup_alloc(void) 
{
    EVDLookup * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(EVDLookup *) ckalloc (sizeof(EVDLookup))) == NULL)  {  
      warn("EVDLookup_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->mu = 0; 
    out->lambda = 0; 


    return out;  
}    


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
EVDLookup * free_EVDLookup(EVDLookup * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a EVDLookup obj. Should be trappable"); 
      return NULL;   
      }  


#ifdef PTHREAD   
    assert(pthread_mutex_lock(&(obj->dynamite_mutex)) == 0); 
#endif   
    if( obj->dynamite_hard_link > 1)     {  
      return_early = 1;  
      obj->dynamite_hard_link--; 
      }  
#ifdef PTHREAD   
    assert(pthread_mutex_unlock(&(obj->dynamite_mutex)) == 0);   
#endif   
    if( return_early == 1)   
      return NULL;   


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
