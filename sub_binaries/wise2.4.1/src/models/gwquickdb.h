#ifndef DYNAMITEgwquickdbHEADERFILE
#define DYNAMITEgwquickdbHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "genewisemodeldb.h"


/* Object GeneWiseQuickDB
 *
 * Descrip: This is a database stream of
 *        GeneWiseScoreFlat objects, layered
 *        ontop of a vanilla genewisedb.
 *
 *        The GeneWiseScoreFlat objects give
 *        around a 10% speed up due to them
 *        being allocated as a single block
 *        of memory that then gets accessed
 *
 *        This object is a very thin layer over
 *        the genewisedb object, which itself handles
 *        the actual HMM stuff via threestatedb object.
 *        So if you want to discover how HMMs are
 *        actually being streamed into the database,
 *        look in there.
 *
 *
 */
struct Wise2_GeneWiseQuickDB {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    GeneWiseDB * gwdb;   
    } ;  
/* GeneWiseQuickDB defined */ 
#ifndef DYNAMITE_DEFINED_GeneWiseQuickDB
typedef struct Wise2_GeneWiseQuickDB Wise2_GeneWiseQuickDB;
#define GeneWiseQuickDB Wise2_GeneWiseQuickDB
#define DYNAMITE_DEFINED_GeneWiseQuickDB
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  init_GeneWiseQuickDB(gdb,return_status)
 *
 * Descrip:    inits a genewisequick database. 
 *
 *
 * Arg:                  gdb [UNKN ] Undocumented argument [GeneWiseQuickDB *]
 * Arg:        return_status [UNKN ] Undocumented argument [int *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseScoreFlat *]
 *
 */
GeneWiseScoreFlat * Wise2_init_GeneWiseQuickDB(GeneWiseQuickDB * gdb,int * return_status);
#define init_GeneWiseQuickDB Wise2_init_GeneWiseQuickDB


/* Function:  reload_GeneWiseQuickDB(prev,gdb,return_status)
 *
 * Descrip:    Reloads a genewisequick database
 *
 *
 *
 * Arg:                 prev [UNKN ] Undocumented argument [GeneWiseScoreFlat *]
 * Arg:                  gdb [UNKN ] Undocumented argument [GeneWiseQuickDB *]
 * Arg:        return_status [UNKN ] Undocumented argument [int *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseScoreFlat *]
 *
 */
GeneWiseScoreFlat * Wise2_reload_GeneWiseQuickDB(GeneWiseScoreFlat * prev,GeneWiseQuickDB * gdb,int * return_status);
#define reload_GeneWiseQuickDB Wise2_reload_GeneWiseQuickDB


/* Function:  close_GeneWiseQuickDB(gws,gdb)
 *
 * Descrip:    closes a GeneWiseDB
 *
 *
 * Arg:        gws [UNKN ] Undocumented argument [GeneWiseScoreFlat *]
 * Arg:        gdb [UNKN ] Undocumented argument [GeneWiseQuickDB *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_close_GeneWiseQuickDB(GeneWiseScoreFlat * gws,GeneWiseQuickDB * gdb);
#define close_GeneWiseQuickDB Wise2_close_GeneWiseQuickDB


/* Function:  dataentry_add_GeneWiseQuickDB(de,gws,gdb)
 *
 * Descrip:    adds dataentry stuff to a query.
 *
 *
 * Arg:         de [UNKN ] Undocumented argument [DataEntry *]
 * Arg:        gws [UNKN ] Undocumented argument [GeneWiseScoreFlat *]
 * Arg:        gdb [UNKN ] Undocumented argument [GeneWiseQuickDB *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_dataentry_add_GeneWiseQuickDB(DataEntry * de,GeneWiseScoreFlat * gws,GeneWiseQuickDB * gdb);
#define dataentry_add_GeneWiseQuickDB Wise2_dataentry_add_GeneWiseQuickDB


/* Function:  GeneWiseQuickDB_from_GeneWiseDB(gwdb)
 *
 * Descrip:    Makes a new genewisequickdb from a genewisemodeldb
 *
 *
 * Arg:        gwdb [READ ] genewisedb - hard links as it enters [GeneWiseDB *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseQuickDB *]
 *
 */
GeneWiseQuickDB * Wise2_GeneWiseQuickDB_from_GeneWiseDB(GeneWiseDB * gwdb);
#define GeneWiseQuickDB_from_GeneWiseDB Wise2_GeneWiseQuickDB_from_GeneWiseDB


/* Function:  hard_link_GeneWiseQuickDB(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GeneWiseQuickDB *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseQuickDB *]
 *
 */
GeneWiseQuickDB * Wise2_hard_link_GeneWiseQuickDB(GeneWiseQuickDB * obj);
#define hard_link_GeneWiseQuickDB Wise2_hard_link_GeneWiseQuickDB


/* Function:  GeneWiseQuickDB_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseQuickDB *]
 *
 */
GeneWiseQuickDB * Wise2_GeneWiseQuickDB_alloc(void);
#define GeneWiseQuickDB_alloc Wise2_GeneWiseQuickDB_alloc


/* Function:  free_GeneWiseQuickDB(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneWiseQuickDB *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseQuickDB *]
 *
 */
GeneWiseQuickDB * Wise2_free_GeneWiseQuickDB(GeneWiseQuickDB * obj);
#define free_GeneWiseQuickDB Wise2_free_GeneWiseQuickDB


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
