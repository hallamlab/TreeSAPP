#ifndef DYNAMITEhspthreadeddbHEADERFILE
#define DYNAMITEhspthreadeddbHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "hsplookupscan.h"
#include "arrayseqlookup.h"

#define HSPThreadedDatabaseLISTLENGTH 32
#define DBTHREAD_MAX_SIZE 64

struct Wise2_HSPDatabaseSegment {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    HSPScanInterface * hspi;     
    SeqLookupInterface * sli;    
    Sequence * query;   /*  for running */ 
    HSPScanInterfacePara * run_para;    /*  for running */ 
    LinearHSPmanager * lm;  /*  for retrieving results */ 
    SeqLookupLoadPara * load_para;  /*  for loading */ 
    SequenceDB * loaddb;    /*  for loading */ 
    pthread_mutex_t * dblock;   /*  for loading */ 
    long int * count;   /*  for loading */ 
    boolean * is_first; /*  for loading */ 
    } ;  
/* HSPDatabaseSegment defined */ 
#ifndef DYNAMITE_DEFINED_HSPDatabaseSegment
typedef struct Wise2_HSPDatabaseSegment Wise2_HSPDatabaseSegment;
#define HSPDatabaseSegment Wise2_HSPDatabaseSegment
#define DYNAMITE_DEFINED_HSPDatabaseSegment
#endif


struct Wise2_HSPThreadedDatabase {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    HSPDatabaseSegment ** seg;   
    int len;/* len for above seg  */ 
    int maxlen; /* maxlen for above seg */ 
    } ;  
/* HSPThreadedDatabase defined */ 
#ifndef DYNAMITE_DEFINED_HSPThreadedDatabase
typedef struct Wise2_HSPThreadedDatabase Wise2_HSPThreadedDatabase;
#define HSPThreadedDatabase Wise2_HSPThreadedDatabase
#define DYNAMITE_DEFINED_HSPThreadedDatabase
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



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
HSPScanInterface * Wise2_new_HSPScanInterface_from_HSPThreadedDatabase(HSPThreadedDatabase * tdb);
#define new_HSPScanInterface_from_HSPThreadedDatabase Wise2_new_HSPScanInterface_from_HSPThreadedDatabase


/* Function:  free_data_hspthreadeddb(d)
 *
 * Descrip:    frees data 
 *
 *
 * Arg:        d [UNKN ] Undocumented argument [void *]
 *
 */
void Wise2_free_data_hspthreadeddb(void * d);
#define free_data_hspthreadeddb Wise2_free_data_hspthreadeddb


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
LinearHSPmanager * Wise2_scan_query_hspthreadeddb(void * d,Sequence * query,HSPScanInterfacePara * para);
#define scan_query_hspthreadeddb Wise2_scan_query_hspthreadeddb


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
void * Wise2_threadeddb_scan_worker(void * d);
#define threadeddb_scan_worker Wise2_threadeddb_scan_worker


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
HSPThreadedDatabase * Wise2_new_HSPThreadedDatabase(int segments,int array_numb_level);
#define new_HSPThreadedDatabase Wise2_new_HSPThreadedDatabase


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
void Wise2_load_HSPThreadedDatabase(HSPThreadedDatabase * db,SequenceDB * sdb,SeqLookupLoadPara * para,CompMat * mat,int drop_off,int score_cutoff);
#define load_HSPThreadedDatabase Wise2_load_HSPThreadedDatabase


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
void * Wise2_threaddb_load_worker(void * d);
#define threaddb_load_worker Wise2_threaddb_load_worker


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
HSPDatabaseSegment * Wise2_hard_link_HSPDatabaseSegment(HSPDatabaseSegment * obj);
#define hard_link_HSPDatabaseSegment Wise2_hard_link_HSPDatabaseSegment


/* Function:  HSPDatabaseSegment_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [HSPDatabaseSegment *]
 *
 */
HSPDatabaseSegment * Wise2_HSPDatabaseSegment_alloc(void);
#define HSPDatabaseSegment_alloc Wise2_HSPDatabaseSegment_alloc


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
HSPDatabaseSegment * Wise2_free_HSPDatabaseSegment(HSPDatabaseSegment * obj);
#define free_HSPDatabaseSegment Wise2_free_HSPDatabaseSegment


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
boolean Wise2_add_HSPThreadedDatabase(HSPThreadedDatabase * obj,HSPDatabaseSegment * add);
#define add_HSPThreadedDatabase Wise2_add_HSPThreadedDatabase


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
int Wise2_flush_HSPThreadedDatabase(HSPThreadedDatabase * obj);
#define flush_HSPThreadedDatabase Wise2_flush_HSPThreadedDatabase


/* Function:  HSPThreadedDatabase_alloc_std(void)
 *
 * Descrip:    Equivalent to HSPThreadedDatabase_alloc_len(HSPThreadedDatabaseLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [HSPThreadedDatabase *]
 *
 */
HSPThreadedDatabase * Wise2_HSPThreadedDatabase_alloc_std(void);
#define HSPThreadedDatabase_alloc_std Wise2_HSPThreadedDatabase_alloc_std


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
HSPThreadedDatabase * Wise2_HSPThreadedDatabase_alloc_len(int len);
#define HSPThreadedDatabase_alloc_len Wise2_HSPThreadedDatabase_alloc_len


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
HSPThreadedDatabase * Wise2_hard_link_HSPThreadedDatabase(HSPThreadedDatabase * obj);
#define hard_link_HSPThreadedDatabase Wise2_hard_link_HSPThreadedDatabase


/* Function:  HSPThreadedDatabase_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [HSPThreadedDatabase *]
 *
 */
HSPThreadedDatabase * Wise2_HSPThreadedDatabase_alloc(void);
#define HSPThreadedDatabase_alloc Wise2_HSPThreadedDatabase_alloc


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
HSPThreadedDatabase * Wise2_free_HSPThreadedDatabase(HSPThreadedDatabase * obj);
#define free_HSPThreadedDatabase Wise2_free_HSPThreadedDatabase


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_HSPThreadedDatabase(HSPDatabaseSegment ** list,int i,int j) ;
#define swap_HSPThreadedDatabase Wise2_swap_HSPThreadedDatabase
void Wise2_qsort_HSPThreadedDatabase(HSPDatabaseSegment ** list,int left,int right,int (*comp)(HSPDatabaseSegment * ,HSPDatabaseSegment * ));
#define qsort_HSPThreadedDatabase Wise2_qsort_HSPThreadedDatabase
void Wise2_sort_HSPThreadedDatabase(HSPThreadedDatabase * obj,int (*comp)(HSPDatabaseSegment *, HSPDatabaseSegment *));
#define sort_HSPThreadedDatabase Wise2_sort_HSPThreadedDatabase
boolean Wise2_expand_HSPThreadedDatabase(HSPThreadedDatabase * obj,int len);
#define expand_HSPThreadedDatabase Wise2_expand_HSPThreadedDatabase

#ifdef _cplusplus
}
#endif

#endif
