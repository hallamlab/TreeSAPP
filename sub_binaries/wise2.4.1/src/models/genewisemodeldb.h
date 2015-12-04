#ifndef DYNAMITEgenewisemodeldbHEADERFILE
#define DYNAMITEgenewisemodeldbHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"
#include "genewisemodel.h"
#include "threestatedb.h"
#include "gwlitemodel.h"


/* Object GeneWiseDB
 *
 * Descrip: This is a database of genewisemodels
 *        for database searching versions of
 *        genewise and estwise type algorithms.
 *
 *        The actual HMM database streaming
 *        happens via the ThreeStateDB (tdb).
 *        This object holds the necessary conversion
 *        system for the HMM database as it
 *        is mapped into the genewisemodel
 *
 *
 */
struct Wise2_GeneWiseDB {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    ThreeStateDB * tdb;  
    double cds_factor;   
    CodonMapper * cm;   /*  NB, hard linked. */ 
    RandomModelDNA * rmd;   /*  NB hard linked */ 
    GeneWiseScore * gws;    /*  hardlinked, only if single. */ 
    GeneParameter21 * gpara;    /*  if a genomic model - est models wont use this! */ 
    boolean is_single;   
    boolean is_syn;  
    Probability allN;    
    boolean flat_insert;     
    GwLiteScore * gwls; /*  only if single */ 
    } ;  
/* GeneWiseDB defined */ 
#ifndef DYNAMITE_DEFINED_GeneWiseDB
typedef struct Wise2_GeneWiseDB Wise2_GeneWiseDB;
#define GeneWiseDB Wise2_GeneWiseDB
#define DYNAMITE_DEFINED_GeneWiseDB
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  init_GeneWiseDB(gdb,return_status)
 *
 * Descrip:    inits a genewise database. Remember this is used
 *             for both genomic and cdna searches
 *
 *
 *
 * Arg:                  gdb [UNKN ] Undocumented argument [GeneWiseDB *]
 * Arg:        return_status [UNKN ] Undocumented argument [int *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseScore *]
 *
 */
GeneWiseScore * Wise2_init_GeneWiseDB(GeneWiseDB * gdb,int * return_status);
#define init_GeneWiseDB Wise2_init_GeneWiseDB


/* Function:  reload_GeneWiseDB(prev,gdb,return_status)
 *
 * Descrip:    Reloads a genewise database
 *
 *
 *
 * Arg:                 prev [UNKN ] Undocumented argument [GeneWiseScore *]
 * Arg:                  gdb [UNKN ] Undocumented argument [GeneWiseDB *]
 * Arg:        return_status [UNKN ] Undocumented argument [int *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseScore *]
 *
 */
GeneWiseScore * Wise2_reload_GeneWiseDB(GeneWiseScore * prev,GeneWiseDB * gdb,int * return_status);
#define reload_GeneWiseDB Wise2_reload_GeneWiseDB


/* Function:  close_GeneWiseDB(gws,gdb)
 *
 * Descrip:    closes a GeneWiseDB
 *
 *
 * Arg:        gws [UNKN ] Undocumented argument [GeneWiseScore *]
 * Arg:        gdb [UNKN ] Undocumented argument [GeneWiseDB *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_close_GeneWiseDB(GeneWiseScore * gws,GeneWiseDB * gdb);
#define close_GeneWiseDB Wise2_close_GeneWiseDB


/* Function:  dataentry_add_GeneWiseDB(de,gws,gdb)
 *
 * Descrip:    adds dataentry stuff to a query. Relies completely
 *             on threestatedb
 *
 *
 * Arg:         de [UNKN ] Undocumented argument [DataEntry *]
 * Arg:        gws [UNKN ] Undocumented argument [GeneWiseScore *]
 * Arg:        gdb [UNKN ] Undocumented argument [GeneWiseDB *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_dataentry_add_GeneWiseDB(DataEntry * de,GeneWiseScore * gws,GeneWiseDB * gdb);
#define dataentry_add_GeneWiseDB Wise2_dataentry_add_GeneWiseDB


/* Function:  new_GeneWiseDB_cdna(syn,tdb,cp,cm,rmd,use_syn,flat_insert,allN)
 *
 * Descrip:    makes a new GeneWiseDB from its component parts,
 *             assumming a cDNA db.
 *
 *             All the objects are hard-linked internally, so you can, if you
 *             wish, free them once passing them into this function
 *
 *
 * Arg:                syn [UNKN ] if ture, use a synchronous coding model vs internally stored tdb rm's [NullString]
 * Arg:                tdb [UNKN ] three state model db to use [ThreeStateDB *]
 * Arg:                 cp [UNKN ] codon parser function to remove from match state [cDNAParser *]
 * Arg:                 cm [UNKN ] Undocumented argument [CodonMapper *]
 * Arg:                rmd [UNKN ] random model (dna) [RandomModelDNA *]
 * Arg:            use_syn [UNKN ] Undocumented argument [boolean]
 * Arg:        flat_insert [UNKN ] Undocumented argument [boolean]
 * Arg:               allN [UNKN ] Undocumented argument [Probability]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseDB *]
 *
 */
GeneWiseDB * Wise2_new_GeneWiseDB_cdna(ThreeStateDB * tdb,cDNAParser * cp,CodonMapper * cm,RandomModelDNA * rmd,boolean use_syn,boolean flat_insert,Probability allN) ;
#define new_GeneWiseDB_cdna Wise2_new_GeneWiseDB_cdna


/* Function:  new_GeneWiseDB(tdb,gp,rmd,use_syn,allN)
 *
 * Descrip:    makes a new GeneWiseDB from its component parts.
 *             All the objects are hard-linked internally, so you can, if you
 *             wish, free them once passing them into this function
 *
 *
 * Arg:            tdb [UNKN ] Undocumented argument [ThreeStateDB *]
 * Arg:             gp [UNKN ] Undocumented argument [GeneParameter21 *]
 * Arg:            rmd [UNKN ] Undocumented argument [RandomModelDNA *]
 * Arg:        use_syn [UNKN ] Undocumented argument [boolean]
 * Arg:           allN [UNKN ] Undocumented argument [Probability]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseDB *]
 *
 */
GeneWiseDB * Wise2_new_GeneWiseDB(ThreeStateDB * tdb,GeneParameter21 * gp,RandomModelDNA * rmd,boolean use_syn,Probability allN);
#define new_GeneWiseDB Wise2_new_GeneWiseDB


/* Function:  new_single_GeneWiseDB(gws)
 *
 * Descrip:    makes a new GeneWiseDB from a single GeneWiseScore.
 *             It hard links it, so you should free it afterwards
 *             in its own scope.
 *
 *
 * Arg:        gws [UNKN ] Undocumented argument [GeneWiseScore *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseDB *]
 *
 */
GeneWiseDB * Wise2_new_single_GeneWiseDB(GeneWiseScore * gws);
#define new_single_GeneWiseDB Wise2_new_single_GeneWiseDB


/* Function:  init_GwLite_GeneWiseDB(gdb,return_status)
 *
 * Descrip:    inits a genewise database for gwlite models
 *
 *
 *
 *
 * Arg:                  gdb [UNKN ] Undocumented argument [GeneWiseDB *]
 * Arg:        return_status [UNKN ] Undocumented argument [int *]
 *
 * Return [UNKN ]  Undocumented return value [GwLiteScore *]
 *
 */
GwLiteScore * Wise2_init_GwLite_GeneWiseDB(GeneWiseDB * gdb,int * return_status);
#define init_GwLite_GeneWiseDB Wise2_init_GwLite_GeneWiseDB


/* Function:  reload_GwLite_GeneWiseDB(prev,gdb,return_status)
 *
 * Descrip:    Reloads a genewise database for a GwLite database
 *
 *
 *
 * Arg:                 prev [UNKN ] Undocumented argument [GwLiteScore *]
 * Arg:                  gdb [UNKN ] Undocumented argument [GeneWiseDB *]
 * Arg:        return_status [UNKN ] Undocumented argument [int *]
 *
 * Return [UNKN ]  Undocumented return value [GwLiteScore *]
 *
 */
GwLiteScore * Wise2_reload_GwLite_GeneWiseDB(GwLiteScore * prev,GeneWiseDB * gdb,int * return_status);
#define reload_GwLite_GeneWiseDB Wise2_reload_GwLite_GeneWiseDB


/* Function:  close_GwLite_GeneWiseDB(gws,gdb)
 *
 * Descrip:    closes a GeneWiseDB
 *
 *
 * Arg:        gws [UNKN ] Undocumented argument [GwLiteScore *]
 * Arg:        gdb [UNKN ] Undocumented argument [GeneWiseDB *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_close_GwLite_GeneWiseDB(GwLiteScore * gws,GeneWiseDB * gdb);
#define close_GwLite_GeneWiseDB Wise2_close_GwLite_GeneWiseDB


/* Function:  dataentry_add_GwLite_GeneWiseDB(de,gws,gdb)
 *
 * Descrip:    adds dataentry stuff to a query. Relies completely
 *             on threestatedb
 *
 *
 * Arg:         de [UNKN ] Undocumented argument [DataEntry *]
 * Arg:        gws [UNKN ] Undocumented argument [GwLiteScore *]
 * Arg:        gdb [UNKN ] Undocumented argument [GeneWiseDB *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_dataentry_add_GwLite_GeneWiseDB(DataEntry * de,GwLiteScore * gws,GeneWiseDB * gdb);
#define dataentry_add_GwLite_GeneWiseDB Wise2_dataentry_add_GwLite_GeneWiseDB


/* Function:  hard_link_GeneWiseDB(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GeneWiseDB *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseDB *]
 *
 */
GeneWiseDB * Wise2_hard_link_GeneWiseDB(GeneWiseDB * obj);
#define hard_link_GeneWiseDB Wise2_hard_link_GeneWiseDB


/* Function:  GeneWiseDB_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseDB *]
 *
 */
GeneWiseDB * Wise2_GeneWiseDB_alloc(void);
#define GeneWiseDB_alloc Wise2_GeneWiseDB_alloc


/* Function:  free_GeneWiseDB(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneWiseDB *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseDB *]
 *
 */
GeneWiseDB * Wise2_free_GeneWiseDB(GeneWiseDB * obj);
#define free_GeneWiseDB Wise2_free_GeneWiseDB


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
GeneWise * Wise2_GeneWise_from_ThreeStateModel_GDB(ThreeStateModel * tsm,GeneWiseDB * gdb);
#define GeneWise_from_ThreeStateModel_GDB Wise2_GeneWise_from_ThreeStateModel_GDB

#ifdef _cplusplus
}
#endif

#endif
