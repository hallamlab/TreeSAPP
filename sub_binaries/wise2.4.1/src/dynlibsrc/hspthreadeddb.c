#ifdef _cplusplus
extern "C" {
#endif
#include "hspthreadeddb.h"

#define HAS_PTHREAD_SETSCOPE
#define HAS_PTHREAD_SETCONCURRENCY

/* Function:  new_HSPScanInterface_from_HSPThreadedDatabase(tdb)
 *
 * Descrip:    Makes a HSPScanInterface from a loaded threaded database
 *
 *
 * Arg:        tdb [UNKN ] Undocumented argument [HSPThreadedDatabase *]
 *
 * Return [UNKN ]  Undocumented return value [HSPScanInterface *]
 *
 */
# line 39 "hspthreadeddb.dy"
HSPScanInterface * new_HSPScanInterface_from_HSPThreadedDatabase(HSPThreadedDatabase * tdb)
{
  HSPScanInterface * out;

  assert(tdb != NULL);
  assert(tdb->seg[0] != NULL);
  assert(tdb->seg[0]->sli != NULL);

  out = HSPScanInterface_alloc();

  out->data = tdb;
  out->scan_query = scan_query_hspthreadeddb;
  out->free_data  = free_data_hspthreadeddb;

  return out;
}




/* Function:  free_data_hspthreadeddb(d)
 *
 * Descrip:    frees data 
 *
 *
 * Arg:        d [UNKN ] Undocumented argument [void *]
 *
 */
# line 62 "hspthreadeddb.dy"
void free_data_hspthreadeddb(void * d)
{
  HSPThreadedDatabase * tdb = (HSPThreadedDatabase *) d;

  free_HSPThreadedDatabase(tdb);
}



/* Function:  scan_query_hspthreadeddb(d,query,para)
 *
 * Descrip:    Does the scan query for the HSPThreadedDatabase
 *
 *
 * Arg:            d [UNKN ] Undocumented argument [void *]
 * Arg:        query [UNKN ] Undocumented argument [Sequence *]
 * Arg:         para [UNKN ] Undocumented argument [HSPScanInterfacePara *]
 *
 * Return [UNKN ]  Undocumented return value [LinearHSPmanager *]
 *
 */
# line 74 "hspthreadeddb.dy"
LinearHSPmanager * scan_query_hspthreadeddb(void * d,Sequence * query,HSPScanInterfacePara * para)
{
  HSPThreadedDatabase * tdb;
  pthread_t thread_pool[DBTHREAD_MAX_SIZE];
  pthread_attr_t pat;
  int i,j;
  int err;
  LinearHSPmanager * out;
  LinearHSPmanager * real_out;


  tdb = (HSPThreadedDatabase *) d;
  pthread_attr_init(&pat);     

#ifdef  HAS_PTHREAD_SETSCOPE 
  pthread_attr_setscope(&pat, PTHREAD_SCOPE_SYSTEM);   
#endif /* set scope */   
  /* Give thread libraries a hint that there are num of threads to run */ 
#ifdef HAS_PTHREAD_SETCONCURRENCY    
  /*  pthread_setconcurrency(tdb->len+1);*/
#endif /* set concurrency */ 

  for(i=0;i<tdb->len;i++) {
    tdb->seg[i]->query = query;
    tdb->seg[i]->run_para = para;
    tdb->seg[i]->lm = NULL;
    if( (err = pthread_create(&(thread_pool[i]),&pat,threadeddb_scan_worker,(void*)tdb->seg[i])) != 0 ) {
      fatal("Unable to make thread %d with error %d",i,err);
    }
  }

  for(i=0;i<tdb->len;i++) {
    pthread_join(thread_pool[i],NULL);
  }

  fprintf(stderr,"All threads have run\n");
  if( VERBOSITY_CHECK(4,para->verbosity) ) {
    info("All threads have run");
  }

  out = LinearHSPmanager_alloc_std();

  for(i=0;i<tdb->len;i++) {
    assert(tdb->seg[i]->lm != NULL);
    for(j=0;j<tdb->seg[i]->lm->len;j++) {
      add_LinearHSPmanager(out,hard_link_HSPset(tdb->seg[i]->lm->set[j]));
    }
  }

  qsort(out->set,out->len,sizeof(HSPset*),compare_HSPset_score_qsort);

  /* this is a bit evil */

  real_out = LinearHSPmanager_alloc_len(para->max_results);
  for(i=0;i<para->max_results && i<out->len;i++) {
    add_LinearHSPmanager(real_out,out->set[i]);
  }

  /*  out->len = para->max_results;*/

  if( VERBOSITY_CHECK(4,para->verbosity) ) {
    info("Finish combine/sort with %d elements",real_out->len);
  }


  return out;
}


/* Function:  threadeddb_scan_worker(d)
 *
 * Descrip:    internal scan function
 *
 *
 * Arg:        d [UNKN ] Undocumented argument [void *]
 *
 * Return [UNKN ]  Undocumented return value [void *]
 *
 */
# line 146 "hspthreadeddb.dy"
void * threadeddb_scan_worker(void * d)
{
  HSPDatabaseSegment * seg;

  seg = (HSPDatabaseSegment *) d;

  seg->lm = (*seg->hspi->scan_query)(seg->hspi->data,seg->query,seg->run_para);

  fprintf(stderr,"For segment %d, finished query with %d (%d) linear\n",seg,(int)seg->lm,seg->lm->len);

  return NULL;
}


/* Function:  new_HSPThreadedDatabase(segments,array_numb_level)
 *
 * Descrip:    Makes a new segmented database suitable for
 *             threading
 *
 *
 * Arg:                segments [UNKN ] Undocumented argument [int]
 * Arg:        array_numb_level [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [HSPThreadedDatabase *]
 *
 */
# line 164 "hspthreadeddb.dy"
HSPThreadedDatabase * new_HSPThreadedDatabase(int segments,int array_numb_level)
{
  HSPThreadedDatabase * out;
  int i;
  HSPDatabaseSegment * seg;

  out = HSPThreadedDatabase_alloc_len(segments);
  
  for(i=0;i<segments;i++) {
    seg = HSPDatabaseSegment_alloc();
    seg->hspi = NULL;
    seg->sli  = new_ArraySeq_SeqLookupInterface(26*26*26*26*26,array_numb_level);
                
    add_HSPThreadedDatabase(out,seg);
  }

  return out;

}

/* Function:  load_HSPThreadedDatabase(db,sdb,para,mat,drop_off,score_cutoff)
 *
 * Descrip:    Loades a segmented database
 *
 *
 * Arg:                  db [UNKN ] Undocumented argument [HSPThreadedDatabase *]
 * Arg:                 sdb [UNKN ] Undocumented argument [SequenceDB *]
 * Arg:                para [UNKN ] Undocumented argument [SeqLookupLoadPara *]
 * Arg:                 mat [UNKN ] Undocumented argument [CompMat *]
 * Arg:            drop_off [UNKN ] Undocumented argument [int]
 * Arg:        score_cutoff [UNKN ] Undocumented argument [int]
 *
 */
# line 187 "hspthreadeddb.dy"
void load_HSPThreadedDatabase(HSPThreadedDatabase * db,SequenceDB * sdb,SeqLookupLoadPara * para,CompMat * mat,int drop_off,int score_cutoff)
{
  pthread_mutex_t dblock;
  pthread_t thread_pool[DBTHREAD_MAX_SIZE];
  pthread_attr_t pat;
  int i;
  int err;

  long int count = 0;
  boolean is_first = 1;

  pthread_attr_init(&pat);     

#ifdef  HAS_PTHREAD_SETSCOPE 
  pthread_attr_setscope(&pat, PTHREAD_SCOPE_SYSTEM);   
#endif /* set scope */   

#ifdef HAS_PTHREAD_SETCONCURRENCY    
  /* needed to make sure one thread doesn't dominate the IO */
  pthread_attr_setschedpolicy(&pat,SCHED_RR);    
#endif /* set concurrency */ 

  if( pthread_mutex_init(&dblock,NULL) != 0 ) {
    fatal("Unable to make mutex for db lock");
  }

  for(i=0;i<db->len;i++) {
    db->seg[i]->load_para = para;
    db->seg[i]->loaddb = sdb;
    db->seg[i]->dblock = &dblock;
    db->seg[i]->count  = &count;
    db->seg[i]->is_first = &is_first;
  }



  for(i=0;i<db->len;i++) {
    if( (err = pthread_create(&(thread_pool[i]),&pat,threaddb_load_worker,(void*)db->seg[i])) != 0 ) {
      fatal("Unable to make thread %d with error %d",i,err);
    }
  }

  for(i=0;i<db->len;i++) {
    pthread_join(thread_pool[i],NULL);
  }

  /* attach scan interfaces to the databases */

  for(i=0;i<db->len;i++) {
    db->seg[i]->hspi = new_one_off_HSPScanInterface(db->seg[i]->sli,mat,drop_off,score_cutoff);
  }

  return;
}

/* Function:  threaddb_load_worker(d)
 *
 * Descrip:    load threaddb 
 *
 *
 * Arg:        d [UNKN ] Undocumented argument [void *]
 *
 * Return [UNKN ]  Undocumented return value [void *]
 *
 */
# line 245 "hspthreadeddb.dy"
void * threaddb_load_worker(void * d)
{
  HSPDatabaseSegment * seg = (HSPDatabaseSegment *) d;
  Sequence * seq = NULL;
  long int temp_count;
  int ret;
  long load_count = 0;

  
  while( 1 ) {
    /* try to get lock */
    if( pthread_mutex_lock(seg->dblock) != 0 ) {
      fatal("Unable to get mutex lock");
    }
    
    /* try to get this sequence */
    if( *seg->is_first == 1 ) {
      seq = init_SequenceDB(seg->loaddb,&ret);
      *seg->is_first = 0;
    }  else {
      seq = get_next_SequenceDB(seg->loaddb);
    }

    /* if seq is NULL, or count > truncate end db */
    if( seq == NULL || (seg->load_para->truncate> 0 && *seg->count > seg->load_para->truncate) ) {
      info("Thread loaded %d entries",load_count);

      pthread_mutex_unlock(seg->dblock);
      return NULL;
    }

    /* if not, up count, check if we need to report */
    (*seg->count)++;

    temp_count = *seg->count;

    /* can now remove lock */

    pthread_mutex_unlock(seg->dblock);

    if( seg->load_para->report_stagger > 0 && (temp_count % seg->load_para->report_stagger == 0) ) {
      info("Threaded db load %d sequences for sli %d at %s",temp_count,seg->sli,seq->name);
    }

    /* now we need to add this sequence */

    add_SeqLookupInterface(seg->sli,seq);

    load_count++;

    (*seg->sli->add_seq)(seg->sli->data,seq,seg->load_para);

  }

    
  return NULL;
}






# line 326 "hspthreadeddb.c"
/* Function:  hard_link_HSPDatabaseSegment(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [HSPDatabaseSegment *]
 *
 * Return [UNKN ]  Undocumented return value [HSPDatabaseSegment *]
 *
 */
HSPDatabaseSegment * hard_link_HSPDatabaseSegment(HSPDatabaseSegment * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a HSPDatabaseSegment object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  HSPDatabaseSegment_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [HSPDatabaseSegment *]
 *
 */
HSPDatabaseSegment * HSPDatabaseSegment_alloc(void) 
{
    HSPDatabaseSegment * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(HSPDatabaseSegment *) ckalloc (sizeof(HSPDatabaseSegment))) == NULL)    {  
      warn("HSPDatabaseSegment_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->hspi = NULL;    
    out->sli = NULL; 


    return out;  
}    


/* Function:  free_HSPDatabaseSegment(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [HSPDatabaseSegment *]
 *
 * Return [UNKN ]  Undocumented return value [HSPDatabaseSegment *]
 *
 */
HSPDatabaseSegment * free_HSPDatabaseSegment(HSPDatabaseSegment * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a HSPDatabaseSegment obj. Should be trappable");    
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
    if( obj->hspi != NULL)   
      free_HSPScanInterface(obj->hspi);  
    if( obj->sli != NULL)    
      free_SeqLookupInterface(obj->sli);     
    /* obj->query is linked in */ 
    /* obj->run_para is linked in */ 
    /* obj->lm is linked in */ 
    /* obj->load_para is linked in */ 
    /* obj->loaddb is linked in */ 
    /* obj->dblock is linked in */ 
    /* obj->count is linked in */ 
    /* obj->is_first is linked in */ 


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_HSPThreadedDatabase(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_HSPThreadedDatabase
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [HSPDatabaseSegment **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_HSPThreadedDatabase(HSPDatabaseSegment ** list,int i,int j)  
{
    HSPDatabaseSegment * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_HSPThreadedDatabase(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_HSPThreadedDatabase which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [HSPDatabaseSegment **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_HSPThreadedDatabase(HSPDatabaseSegment ** list,int left,int right,int (*comp)(HSPDatabaseSegment * ,HSPDatabaseSegment * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_HSPThreadedDatabase(list,left,(left+right)/2);  
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_HSPThreadedDatabase (list,++last,i);    
      }  
    swap_HSPThreadedDatabase (list,left,last);   
    qsort_HSPThreadedDatabase(list,left,last-1,comp);    
    qsort_HSPThreadedDatabase(list,last+1,right,comp);   
}    


/* Function:  sort_HSPThreadedDatabase(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_HSPThreadedDatabase
 *
 *
 * Arg:         obj [UNKN ] Object containing list [HSPThreadedDatabase *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_HSPThreadedDatabase(HSPThreadedDatabase * obj,int (*comp)(HSPDatabaseSegment *, HSPDatabaseSegment *)) 
{
    qsort_HSPThreadedDatabase(obj->seg,0,obj->len-1,comp);   
    return;  
}    


/* Function:  expand_HSPThreadedDatabase(obj,len)
 *
 * Descrip:    Really an internal function for add_HSPThreadedDatabase
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [HSPThreadedDatabase *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_HSPThreadedDatabase(HSPThreadedDatabase * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_HSPThreadedDatabase called with no need");    
      return TRUE;   
      }  


    if( (obj->seg = (HSPDatabaseSegment ** ) ckrealloc (obj->seg,sizeof(HSPDatabaseSegment *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_HSPThreadedDatabase, returning FALSE");  
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_HSPThreadedDatabase(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [HSPThreadedDatabase *]
 * Arg:        add [OWNER] Object to add to the list [HSPDatabaseSegment *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_HSPThreadedDatabase(HSPThreadedDatabase * obj,HSPDatabaseSegment * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_HSPThreadedDatabase(obj,obj->len + HSPThreadedDatabaseLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->seg[obj->len++]=add;    
    return TRUE; 
}    


/* Function:  flush_HSPThreadedDatabase(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [HSPThreadedDatabase *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_HSPThreadedDatabase(HSPThreadedDatabase * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->seg[i] != NULL)   {  
        free_HSPDatabaseSegment(obj->seg[i]);    
        obj->seg[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  HSPThreadedDatabase_alloc_std(void)
 *
 * Descrip:    Equivalent to HSPThreadedDatabase_alloc_len(HSPThreadedDatabaseLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [HSPThreadedDatabase *]
 *
 */
HSPThreadedDatabase * HSPThreadedDatabase_alloc_std(void) 
{
    return HSPThreadedDatabase_alloc_len(HSPThreadedDatabaseLISTLENGTH); 
}    


/* Function:  HSPThreadedDatabase_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [HSPThreadedDatabase *]
 *
 */
HSPThreadedDatabase * HSPThreadedDatabase_alloc_len(int len) 
{
    HSPThreadedDatabase * out;  /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = HSPThreadedDatabase_alloc()) == NULL)  
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->seg = (HSPDatabaseSegment ** ) ckcalloc (len,sizeof(HSPDatabaseSegment *))) == NULL)    {  
      warn("Warning, ckcalloc failed in HSPThreadedDatabase_alloc_len"); 
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_HSPThreadedDatabase(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [HSPThreadedDatabase *]
 *
 * Return [UNKN ]  Undocumented return value [HSPThreadedDatabase *]
 *
 */
HSPThreadedDatabase * hard_link_HSPThreadedDatabase(HSPThreadedDatabase * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a HSPThreadedDatabase object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  HSPThreadedDatabase_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [HSPThreadedDatabase *]
 *
 */
HSPThreadedDatabase * HSPThreadedDatabase_alloc(void) 
{
    HSPThreadedDatabase * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(HSPThreadedDatabase *) ckalloc (sizeof(HSPThreadedDatabase))) == NULL)  {  
      warn("HSPThreadedDatabase_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->seg = NULL; 
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_HSPThreadedDatabase(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [HSPThreadedDatabase *]
 *
 * Return [UNKN ]  Undocumented return value [HSPThreadedDatabase *]
 *
 */
HSPThreadedDatabase * free_HSPThreadedDatabase(HSPThreadedDatabase * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a HSPThreadedDatabase obj. Should be trappable");   
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
    if( obj->seg != NULL)    {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->seg[i] != NULL) 
          free_HSPDatabaseSegment(obj->seg[i]);  
        }  
      ckfree(obj->seg);  
      }  


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
