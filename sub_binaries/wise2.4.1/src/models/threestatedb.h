#ifndef DYNAMITEthreestatedbHEADERFILE
#define DYNAMITEthreestatedbHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "threestatemodel.h"
#include "pfamhmmer1db.h"


typedef enum {
  TSMDB_UNKNOWN,
  TSMDB_SINGLE,
  TSMDB_HMMER1PFAM,
  TSMDB_PROTEIN,
  TSMDB_GENERIC
} TSMDB_Type;

#ifndef DYNAMITE_DEFINED_ThreeStateDB
typedef struct Wise2_ThreeStateDB Wise2_ThreeStateDB;
#define ThreeStateDB Wise2_ThreeStateDB
#define DYNAMITE_DEFINED_ThreeStateDB
#endif

#ifndef DYNAMITE_DEFINED_PfamHmmer1DB
typedef struct Wise2_PfamHmmer1DB Wise2_PfamHmmer1DB;
#define PfamHmmer1DB Wise2_PfamHmmer1DB
#define DYNAMITE_DEFINED_PfamHmmer1DB
#endif

/* Object ThreeStateDB
 *
 * Descrip: ThreeStateDB is the object that represents
 *        a database of profile-HMMs. 
 *
 *        The object hold a variety of fields on some of which are
 *        occupied depending on the type.
 *
 *        Realistically we need a more abstract class idea, which is
 *        implemented here anyway via the generic stuff, in hacky
 *        C-style pointers to function plus a void pointer. This object
 *        therefore houses a switch system around the different types
 *        including the generic system... but as the generic function
 *        stuff was bolted on later, some things are handled with
 *        explicit datastructures. It is quite messy ;). Apologies.
 *        To be cleaned up.
 *
 *        The generic stuff was principly added in to allow a decoupling of this module
 *        from the HMMer2.o interface code which is held in wise2xhmmer.dy
 *
 *        The old static datastructure code can be 
 *        made via protein sequences which are then converted or a 
 *        Pfam 2.0 style directory + HMMs file.
 *
 *
 */
struct Wise2_ThreeStateDB {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int dbtype;  
    char * filename;     
    int type;    
    FILE * current_file;     
    RandomModel * rm;   /*  NB, this is hard-linked...  */ 
    long byte_position; /*  this is the file position for the current model */ 
    ThreeStateModel * single;   /*  for single db cases */ 
    PfamHmmer1DB * phdb;    /*  for Pfam Hmmer1 style databases. */ 
    SequenceDB   * sdb; /*  for protein databases */ 
    CompMat      * comp;    /*  for protein databases */ 
    int gap;    /*  for protein databases */ 
    int ext;    /*  for protein databases */ 
    Sequence *  seq_cache;  /*  needed for a bit of inter-function communication  */ 
    ThreeStateModel * (*reload_generic)(ThreeStateDB * tdb,int * return_status); 
    boolean (*open_generic)(ThreeStateDB * tdb); 
    boolean (*close_generic)(ThreeStateDB * tdb);    
    boolean (*dataentry_add)(ThreeStateDB * tdb,DataEntry * en); 
    boolean (*open_index_generic)(ThreeStateDB *tdb);    
    ThreeStateModel * (*index_generic)(ThreeStateDB *tdb,DataEntry *de); 
    boolean (*close_index_generic)(ThreeStateDB *tdb);   
    void * data;    /*  whatever else the damn system wants to carry around with it!  */ 
    int  hmm_model_start;    
    int  hmm_model_end;  
    int  current_no;     
    } ;  
/* ThreeStateDB defined */ 
#ifndef DYNAMITE_DEFINED_ThreeStateDB
typedef struct Wise2_ThreeStateDB Wise2_ThreeStateDB;
#define ThreeStateDB Wise2_ThreeStateDB
#define DYNAMITE_DEFINED_ThreeStateDB
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  open_for_indexing_ThreeStateDB(mdb)
 *
 * Descrip:    opens database ready for calls to 
 *             /indexed_ThreeStateModel_ThreeStateDB
 *
 *
 *
 * Arg:        mdb [UNKN ] database to open for index calls [ThreeStateDB *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_open_for_indexing_ThreeStateDB(ThreeStateDB * mdb);
#define open_for_indexing_ThreeStateDB Wise2_open_for_indexing_ThreeStateDB


/* Function:  close_for_indexing_ThreeStateDB(mdb)
 *
 * Descrip:    closes an indexable database
 *
 *
 * Arg:        mdb [UNKN ] Undocumented argument [ThreeStateDB *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_close_for_indexing_ThreeStateDB(ThreeStateDB * mdb);
#define close_for_indexing_ThreeStateDB Wise2_close_for_indexing_ThreeStateDB


/* Function:  set_search_type_ThreeStateDB(tdb,type)
 *
 * Descrip:    Set the search type of this threestatedb...
 *
 *
 * Arg:         tdb [UNKN ] Undocumented argument [ThreeStateDB *]
 * Arg:        type [UNKN ] to set can be any of the modes found in threestatemodel [char *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_set_search_type_ThreeStateDB(ThreeStateDB * tdb,char * type);
#define set_search_type_ThreeStateDB Wise2_set_search_type_ThreeStateDB


/* Function:  indexed_ThreeStateModel_ThreeStateDB(mdb,en)
 *
 * Descrip:    Retrieves a model from a database which has been opened
 *             for indexing by /open_for_indexing_ThreeStateDB
 *
 *             The index information comes from the dataentry which should 
 *             have been from a search of the ThreeStateDB.
 *
 *
 * Arg:        mdb [UNKN ] database where this is indexed [ThreeStateDB *]
 * Arg:         en [UNKN ] dataentry to pull the model from [DataEntry *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateModel *]
 *
 */
ThreeStateModel * Wise2_indexed_ThreeStateModel_ThreeStateDB(ThreeStateDB * mdb,DataEntry * en);
#define indexed_ThreeStateModel_ThreeStateDB Wise2_indexed_ThreeStateModel_ThreeStateDB


/* Function:  new_proteindb_ThreeStateDB(sdb,comp,gap,ext)
 *
 * Descrip:    makes a new ThreeStateDB from a
 *             sequencedb (better be protein!)
 *
 *
 *
 * Arg:         sdb [READ ] sequence database to use [SequenceDB *]
 * Arg:        comp [READ ] comparison matrix to use [CompMat *]
 * Arg:         gap [READ ] gap open penalty [int]
 * Arg:         ext [READ ] gap extensions penalty [int]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateDB *]
 *
 */
ThreeStateDB * Wise2_new_proteindb_ThreeStateDB(SequenceDB * sdb,CompMat * comp,int gap,int ext);
#define new_proteindb_ThreeStateDB Wise2_new_proteindb_ThreeStateDB


/* Function:  new_single_ThreeStateDB(tsm,rm)
 *
 * Descrip:    Making a new ThreeStateDB from a single
 *             model
 *
 *
 *
 * Arg:        tsm [READ ] a single ThreeStateModel [ThreeStateModel *]
 * Arg:         rm [READ ] random model to be used in comparisons.. [RandomModel *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateDB *]
 *
 */
ThreeStateDB * Wise2_new_single_ThreeStateDB(ThreeStateModel * tsm,RandomModel * rm);
#define new_single_ThreeStateDB Wise2_new_single_ThreeStateDB


/* Function:  new_PfamHmmer1DB_ThreeStateDB(dirname)
 *
 * Descrip:    Makes a new PfamHmmer1DB from a filename
 *             indicating the directory
 *
 *
 * Arg:        dirname [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateDB *]
 *
 */
ThreeStateDB * Wise2_new_PfamHmmer1DB_ThreeStateDB(char * dirname);
#define new_PfamHmmer1DB_ThreeStateDB Wise2_new_PfamHmmer1DB_ThreeStateDB


/* Function:  dataentry_add_ThreeStateDB(de,tss,mdb)
 *
 * Descrip:    This function adds the internal entry information 
 *             (eg indexing point) into the dataentry
 *
 *
 * Arg:         de [UNKN ] Undocumented argument [DataEntry *]
 * Arg:        tss [UNKN ] Undocumented argument [ThreeStateScore *]
 * Arg:        mdb [UNKN ] Undocumented argument [ThreeStateDB *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_dataentry_add_ThreeStateDB(DataEntry * de,ThreeStateScore * tss,ThreeStateDB * mdb);
#define dataentry_add_ThreeStateDB Wise2_dataentry_add_ThreeStateDB


/* Function:  open_ThreeStateDB(mdb)
 *
 * Descrip:    Open function for ThreeStateDB.
 *             An internal for this file but also
 *             used by, for example, GeneWiseDB that
 *             wants to get at the underlying models, 
 *             not the log-odds.
 *
 *
 * Arg:        mdb [UNKN ] Undocumented argument [ThreeStateDB *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_open_ThreeStateDB(ThreeStateDB * mdb);
#define open_ThreeStateDB Wise2_open_ThreeStateDB


/* Function:  read_TSM_ThreeStateDB(mdb,return_status)
 *
 * Descrip:    Reads a threestatemodel out from the 
 *             database. People will probably want the
 *             ThreeStateScore *not* the model, but some
 *             systems will want the model.
 *
 *
 * Arg:                  mdb [UNKN ] Undocumented argument [ThreeStateDB *]
 * Arg:        return_status [UNKN ] Undocumented argument [int *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateModel *]
 *
 */
ThreeStateModel * Wise2_read_TSM_ThreeStateDB(ThreeStateDB * mdb,int * return_status);
#define read_TSM_ThreeStateDB Wise2_read_TSM_ThreeStateDB


/* Function:  init_ThreeStateDB(mdb,return_status)
 *
 * Descrip:    Init function for ThreeStateDB
 *
 *             Is going to open file, read first model, complain if
 *             NULL, and convert to a score system.
 *
 *
 * Arg:                  mdb [RW   ] Model database [ThreeStateDB *]
 * Arg:        return_status [WRITE] return from database.h system [int *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateScore *]
 *
 */
ThreeStateScore * Wise2_init_ThreeStateDB(ThreeStateDB * mdb,int * return_status);
#define init_ThreeStateDB Wise2_init_ThreeStateDB


/* Function:  reload_ThreeStateDB(prev,tss,mdb,return_status)
 *
 * Descrip:    reloads the ThreeStateDB.
 *
 *             Frees the previous score system (could recycle memory).
 *             Reads database. calls END if gets NULL from read_HMF_ThreeStateModel
 *
 *
 * Arg:                 prev [UNKN ] Undocumented argument [ThreeStateScore *]
 * Arg:                  tss [UNKN ] the previous score system [NullString]
 * Arg:                  mdb [UNKN ] model database system [ThreeStateDB *]
 * Arg:        return_status [WRITE] return from database.h system [int *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateScore *]
 *
 */
ThreeStateScore * Wise2_reload_ThreeStateDB(ThreeStateScore * prev,ThreeStateDB * mdb,int * return_status);
#define reload_ThreeStateDB Wise2_reload_ThreeStateDB


/* Function:  close_ThreeStateDB(prev,mdb)
 *
 * Descrip:    closes ThreeStateDB
 *
 *             At the moment, only needs to free previous
 *             and close the file
 *
 *
 * Arg:        prev [UNKN ] the last ThreeStateScore to be freed [ThreeStateScore *]
 * Arg:         mdb [UNKN ] Model database [ThreeStateDB *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_close_ThreeStateDB(ThreeStateScore * prev,ThreeStateDB * mdb);
#define close_ThreeStateDB Wise2_close_ThreeStateDB


/* Function:  hard_link_ThreeStateDB(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ThreeStateDB *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateDB *]
 *
 */
ThreeStateDB * Wise2_hard_link_ThreeStateDB(ThreeStateDB * obj);
#define hard_link_ThreeStateDB Wise2_hard_link_ThreeStateDB


/* Function:  ThreeStateDB_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateDB *]
 *
 */
ThreeStateDB * Wise2_ThreeStateDB_alloc(void);
#define ThreeStateDB_alloc Wise2_ThreeStateDB_alloc


/* Function:  free_ThreeStateDB(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ThreeStateDB *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateDB *]
 *
 */
ThreeStateDB * Wise2_free_ThreeStateDB(ThreeStateDB * obj);
#define free_ThreeStateDB Wise2_free_ThreeStateDB


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
int Wise2_access_dbtype_ThreeStateDB(ThreeStateDB * obj);
#define access_dbtype_ThreeStateDB Wise2_access_dbtype_ThreeStateDB
boolean Wise2_replace_filename_ThreeStateDB(ThreeStateDB * obj,char * filename);
#define replace_filename_ThreeStateDB Wise2_replace_filename_ThreeStateDB
boolean Wise2_replace_dbtype_ThreeStateDB(ThreeStateDB * obj,int dbtype);
#define replace_dbtype_ThreeStateDB Wise2_replace_dbtype_ThreeStateDB
char * Wise2_access_filename_ThreeStateDB(ThreeStateDB * obj);
#define access_filename_ThreeStateDB Wise2_access_filename_ThreeStateDB

#ifdef _cplusplus
}
#endif

#endif
