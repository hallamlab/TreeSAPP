#ifdef _cplusplus
extern "C" {
#endif
#include "hspscanruntime.h"


/* Function:  new_runtime_HSPScanInterface(sli,mat,drop_off,score_cutoff,threadno)
 *
 * Descrip:    Makes a new function which will at runtime switch between
 *             implementation; vanilla, threaded and twohit
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
# line 26 "hspscanruntime.dy"
HSPScanInterface * new_runtime_HSPScanInterface(SeqLookupInterface * sli,CompMat * mat,int drop_off,int score_cutoff,int threadno)
{
  HSPScanInterface * out;
  HSPScanRuntimeImpl * rt;

  out = HSPScanInterface_alloc();

  rt = HSPScanRuntimeImpl_alloc();
  rt->vanilla  = new_one_off_HSPScanInterface(sli,mat,drop_off,score_cutoff);
  rt->threaded = new_threaded_HSPScanInterface(sli,mat,drop_off,score_cutoff,threadno);
  rt->twohit   = new_twohit_one_off_HSPScanInterface(sli,mat,drop_off,score_cutoff);

  out->data = (void*) rt;
  out->free_data = free_runtime_hspscan;
  out->scan_query = scan_query_runtime_hspscan;

  return out;
}

/* Function:  scan_query_runtime_hspscan(data,seq,para)
 *
 * Descrip:    Handles runtime switching between methods for the scan query
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 * Arg:         seq [UNKN ] Undocumented argument [Sequence *]
 * Arg:        para [UNKN ] Undocumented argument [HSPScanInterfacePara *]
 *
 * Return [UNKN ]  Undocumented return value [LinearHSPmanager *]
 *
 */
# line 48 "hspscanruntime.dy"
LinearHSPmanager * scan_query_runtime_hspscan(void * data,Sequence * seq,HSPScanInterfacePara * para)
{
  HSPScanRuntimeImpl * d;

  d = (HSPScanRuntimeImpl *) data;

  switch(para->implementation) {
  case HSPSCAN_IMPLEMENTATION_VANILLA :
    return (*d->vanilla->scan_query)(d->vanilla->data,seq,para);
  case HSPSCAN_IMPLEMENTATION_THREADED :
    return (*d->threaded->scan_query)(d->threaded->data,seq,para);
  case HSPSCAN_IMPLEMENTATION_TWOHIT  :
    return (*d->twohit->scan_query)(d->twohit->data,seq,para);
  default :
    warn("No good implementation for %d as implementation",para->implementation);
  }

  return NULL;
}


/* Function:  free_runtime_hspscan(data)
 *
 * Descrip:    free function for runtime
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 *
 */
# line 72 "hspscanruntime.dy"
void free_runtime_hspscan(void * data)
{
  HSPScanRuntimeImpl * d;

  d = (HSPScanRuntimeImpl *) data;

  free_HSPScanRuntimeImpl(d);
}


  
# line 89 "hspscanruntime.c"
/* Function:  hard_link_HSPScanRuntimeImpl(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [HSPScanRuntimeImpl *]
 *
 * Return [UNKN ]  Undocumented return value [HSPScanRuntimeImpl *]
 *
 */
HSPScanRuntimeImpl * hard_link_HSPScanRuntimeImpl(HSPScanRuntimeImpl * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a HSPScanRuntimeImpl object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  HSPScanRuntimeImpl_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [HSPScanRuntimeImpl *]
 *
 */
HSPScanRuntimeImpl * HSPScanRuntimeImpl_alloc(void) 
{
    HSPScanRuntimeImpl * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(HSPScanRuntimeImpl *) ckalloc (sizeof(HSPScanRuntimeImpl))) == NULL)    {  
      warn("HSPScanRuntimeImpl_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->vanilla = NULL; 
    out->threaded = NULL;    
    out->twohit = NULL;  


    return out;  
}    


/* Function:  free_HSPScanRuntimeImpl(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [HSPScanRuntimeImpl *]
 *
 * Return [UNKN ]  Undocumented return value [HSPScanRuntimeImpl *]
 *
 */
HSPScanRuntimeImpl * free_HSPScanRuntimeImpl(HSPScanRuntimeImpl * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a HSPScanRuntimeImpl obj. Should be trappable");    
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
    if( obj->vanilla != NULL)    
      free_HSPScanInterface(obj->vanilla);   
    if( obj->threaded != NULL)   
      free_HSPScanInterface(obj->threaded);  
    if( obj->twohit != NULL) 
      free_HSPScanInterface(obj->twohit);    


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
