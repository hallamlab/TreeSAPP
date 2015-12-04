#ifdef _cplusplus
extern "C" {
#endif
#include "client_multihspscan.h"


# line 29 "client_multihspscan.dy"
LinearHSPmanager * scan_MultiHSPScanClient(void * data,Sequence * seq,HSPScanInterfacePara * para)
{
  LinearHSPmanager * out;
  MultiHSPScanClient * mhsc;
  pthread_attr_t pat;
  pthread_t thread_pool[MAX_CLIENT_HSP_THREADS];
  int i;
  int j;
  int total_len;
  int err;

  mhsc = (MultiHSPScanClient *) data;
  assert(mhsc != NULL);
  assert(seq != NULL);
  assert(mhsc->len < MAX_CLIENT_HSP_THREADS);

  pthread_attr_init(&pat);     

#ifdef  HAS_PTHREAD_SETSCOPE 
  pthread_attr_setscope(&pat, PTHREAD_SCOPE_SYSTEM);   
#endif /* set scope */   
  /* Give thread libraries a hint that there are num of threads to run */ 
#ifdef HAS_PTHREAD_SETCONCURRENCY    
  pthread_setconcurrency(mhsc->len+1);    
#endif /* set concurrency */ 

  for(i=0,total_len = 0;i<mhsc->len;i++) {
    mhsc->client[i]->seq = seq;
    mhsc->client[i]->p   = para;
    if( (err = pthread_create(&(thread_pool[i]),&pat,client_hsp_thread_worker,(void*)mhsc->client[i])) ) {
      fatal("Unable to make thread %d with error %d",i,err);
    }
  }

  for(i=0;i<mhsc->len;i++) {
    if( pthread_join(thread_pool[i],NULL) != 0 ) {
      fatal("Unable to join thread in client hsp");
    }
    total_len += mhsc->client[i]->lm->len;
  }

  out = LinearHSPmanager_alloc_len(total_len);

  for(i=0;i<mhsc->len;i++) {
    for(j=0;j<mhsc->client[i]->lm->len;j++) {
      add_LinearHSPmanager(out,hard_link_HSPset(mhsc->client[i]->lm->set[j]));
    }
  }

  qsort(out->set,out->len,sizeof(HSPset*),compare_HSPset_score_qsort);

  return out;
}

# line 83 "client_multihspscan.dy"
void * client_hsp_thread_worker(void * d)
{
  HSPScanClient * cl;

  cl = (HSPScanClient *) d;

  assert(cl != NULL);

  if( cl->p->verbosity > 4 ) {
    info("Client side thread starting for %s %d",cl->host,cl->port);
  }


  cl->lm = (*cl->hspi->scan_query)(cl->hspi->data,cl->seq,cl->p);

  if( cl->p->verbosity > 4 ) {
    info("Client side thread ending for %s %d",cl->host,cl->port);
  }

  return 0;
}



# line 107 "client_multihspscan.dy"
HSPScanInterface * new_multiclient_HSPScanInterface(char * filename)
{
  HSPScanInterface * out;
  FILE * ifp;
  MultiHSPScanClient * mhsc;

  if( (ifp = openfile(filename,"r")) == NULL ) {
    warn("Unable to open %s as file for multi scan client",filename);
    return NULL;
  }

  mhsc = new_MultiHSPScanClient_from_file(ifp);

  assert(mhsc != NULL);

  out             = HSPScanInterface_alloc();
  out->data       = (void *) mhsc;
  out->scan_query = scan_MultiHSPScanClient;
  out->free_data  = free_multihspscanclient;
  

}

# line 130 "client_multihspscan.dy"
void free_multihspscanclient(void * data)
{
  MultiHSPScanClient * d;

  d = (MultiHSPScanClient *) data;

  free_MultiHSPScanClient(d);

}

# line 140 "client_multihspscan.dy"
MultiHSPScanClient * new_MultiHSPScanClient_from_file(FILE * ifp)
{
  MultiHSPScanClient * out;
  HSPScanClient * cl;
  char buffer[MAXLINE];
  char * name;
  char * port;
  char * r;

  out = MultiHSPScanClient_alloc_std();

  while( fgets(buffer,MAXLINE,ifp) != NULL ) {
    if( strstartcmp(buffer,"//") == 0 ) {
      break;
    }
    name = buffer;
    for(r = buffer;*r && !isspace(*r);r++) {
      ;
    }
    *r = '\0';
    r++;
    for(;*r && isspace(*r);r++) {
      ;
    }
    port = r;
    for(;*r && !isspace(*r);r++) {
      ;
    }
    *r = '\0';

    cl = HSPScanClient_alloc();
    cl->host = stringalloc(name);
    cl->port = strtol(port,NULL,0);

    cl->hspi = new_wise_transfer_HSPScanInterface(cl->host,cl->port);
    if( cl->hspi == NULL ) {
      warn("Unable to make client with %s host, %d port",cl->host,cl->port);
      free_MultiHSPScanClient(out);
      return NULL;
    }

    add_MultiHSPScanClient(out,cl);

  }
  
  return out;
}


# line 171 "client_multihspscan.c"
/* Function:  hard_link_HSPScanClient(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [HSPScanClient *]
 *
 * Return [UNKN ]  Undocumented return value [HSPScanClient *]
 *
 */
HSPScanClient * hard_link_HSPScanClient(HSPScanClient * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a HSPScanClient object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  HSPScanClient_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [HSPScanClient *]
 *
 */
HSPScanClient * HSPScanClient_alloc(void) 
{
    HSPScanClient * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(HSPScanClient *) ckalloc (sizeof(HSPScanClient))) == NULL)  {  
      warn("HSPScanClient_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->host = NULL;    
    out->port = 0;   
    out->hspi = NULL;    
    out->lm = NULL;  


    return out;  
}    


/* Function:  free_HSPScanClient(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [HSPScanClient *]
 *
 * Return [UNKN ]  Undocumented return value [HSPScanClient *]
 *
 */
HSPScanClient * free_HSPScanClient(HSPScanClient * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a HSPScanClient obj. Should be trappable"); 
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
    if( obj->host != NULL)   
      ckfree(obj->host);     
    if( obj->hspi != NULL)   
      free_HSPScanInterface(obj->hspi);  
    if( obj->lm != NULL) 
      free_LinearHSPmanager(obj->lm);    
    /* obj->seq is linked in */ 
    /* obj->p is linked in */ 


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_MultiHSPScanClient(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_MultiHSPScanClient
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [HSPScanClient **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_MultiHSPScanClient(HSPScanClient ** list,int i,int j)  
{
    HSPScanClient * temp;    
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_MultiHSPScanClient(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_MultiHSPScanClient which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [HSPScanClient **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_MultiHSPScanClient(HSPScanClient ** list,int left,int right,int (*comp)(HSPScanClient * ,HSPScanClient * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_MultiHSPScanClient(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_MultiHSPScanClient (list,++last,i); 
      }  
    swap_MultiHSPScanClient (list,left,last);    
    qsort_MultiHSPScanClient(list,left,last-1,comp); 
    qsort_MultiHSPScanClient(list,last+1,right,comp);    
}    


/* Function:  sort_MultiHSPScanClient(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_MultiHSPScanClient
 *
 *
 * Arg:         obj [UNKN ] Object containing list [MultiHSPScanClient *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_MultiHSPScanClient(MultiHSPScanClient * obj,int (*comp)(HSPScanClient *, HSPScanClient *)) 
{
    qsort_MultiHSPScanClient(obj->client,0,obj->len-1,comp); 
    return;  
}    


/* Function:  expand_MultiHSPScanClient(obj,len)
 *
 * Descrip:    Really an internal function for add_MultiHSPScanClient
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [MultiHSPScanClient *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_MultiHSPScanClient(MultiHSPScanClient * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_MultiHSPScanClient called with no need"); 
      return TRUE;   
      }  


    if( (obj->client = (HSPScanClient ** ) ckrealloc (obj->client,sizeof(HSPScanClient *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_MultiHSPScanClient, returning FALSE");   
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_MultiHSPScanClient(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [MultiHSPScanClient *]
 * Arg:        add [OWNER] Object to add to the list [HSPScanClient *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_MultiHSPScanClient(MultiHSPScanClient * obj,HSPScanClient * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_MultiHSPScanClient(obj,obj->len + MultiHSPScanClientLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->client[obj->len++]=add; 
    return TRUE; 
}    


/* Function:  flush_MultiHSPScanClient(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [MultiHSPScanClient *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_MultiHSPScanClient(MultiHSPScanClient * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->client[i] != NULL)    {  
        free_HSPScanClient(obj->client[i]);  
        obj->client[i] = NULL;   
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  MultiHSPScanClient_alloc_std(void)
 *
 * Descrip:    Equivalent to MultiHSPScanClient_alloc_len(MultiHSPScanClientLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [MultiHSPScanClient *]
 *
 */
MultiHSPScanClient * MultiHSPScanClient_alloc_std(void) 
{
    return MultiHSPScanClient_alloc_len(MultiHSPScanClientLISTLENGTH);   
}    


/* Function:  MultiHSPScanClient_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [MultiHSPScanClient *]
 *
 */
MultiHSPScanClient * MultiHSPScanClient_alloc_len(int len) 
{
    MultiHSPScanClient * out;   /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = MultiHSPScanClient_alloc()) == NULL)   
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->client = (HSPScanClient ** ) ckcalloc (len,sizeof(HSPScanClient *))) == NULL)   {  
      warn("Warning, ckcalloc failed in MultiHSPScanClient_alloc_len");  
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_MultiHSPScanClient(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [MultiHSPScanClient *]
 *
 * Return [UNKN ]  Undocumented return value [MultiHSPScanClient *]
 *
 */
MultiHSPScanClient * hard_link_MultiHSPScanClient(MultiHSPScanClient * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a MultiHSPScanClient object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  MultiHSPScanClient_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [MultiHSPScanClient *]
 *
 */
MultiHSPScanClient * MultiHSPScanClient_alloc(void) 
{
    MultiHSPScanClient * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(MultiHSPScanClient *) ckalloc (sizeof(MultiHSPScanClient))) == NULL)    {  
      warn("MultiHSPScanClient_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->client = NULL;  
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_MultiHSPScanClient(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [MultiHSPScanClient *]
 *
 * Return [UNKN ]  Undocumented return value [MultiHSPScanClient *]
 *
 */
MultiHSPScanClient * free_MultiHSPScanClient(MultiHSPScanClient * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a MultiHSPScanClient obj. Should be trappable");    
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
    if( obj->client != NULL) {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->client[i] != NULL)  
          free_HSPScanClient(obj->client[i]);    
        }  
      ckfree(obj->client);   
      }  


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
