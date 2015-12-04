#ifndef DYNAMITEcdnadbHEADERFILE
#define DYNAMITEcdnadbHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "sequencedb.h"
#include "cdna.h"
#include "hscore.h"
#include "complexsequence.h"
#include "complexevalset.h"

typedef enum CdnaDBErrorType {
  CDNADB_READ_THROUGH = 0,
  CDNADB_FAIL_ON_ERROR = 1
} CdnaDBErrorType;

/* Object cDNADB
 *
 * Descrip: This object hold a database of
 *        cDNA sequences.
 *
 *        You will probably use it in one of
 *        two ways
 *
 *        1 A sequence formatted database, which
 *        is provided by a /SequenceDB object
 *        is used to provide the raw sequences 
 *
 *        2 A single cDNA sequence is used.
 *
 *        In each case this database provides
 *        both the forward and reverse strands
 *        into the system.
 *
 *        Notice that what is exported are
 *        /ComplexSequence objects, not cDNA dna,
 *        as this is what is generally needed. 
 *        These are things with splice sites calculated
 *        etc. This is why for initialisation this needs
 *        a /ComplexSequenceEvalSet of the correct type.
 *
 *
 */
struct Wise2_cDNADB {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    boolean is_single_seq;   
    boolean done_forward;    
    boolean forward_only;    
    ComplexSequence * forw;  
    ComplexSequence * rev;   
    SequenceDB * sdb;    
    Sequence * current;  
    ComplexSequenceEvalSet * cses;   
    CdnaDBErrorType error_handling;  
    double error_tol;    
    } ;  
/* cDNADB defined */ 
#ifndef DYNAMITE_DEFINED_cDNADB
typedef struct Wise2_cDNADB Wise2_cDNADB;
#define cDNADB Wise2_cDNADB
#define DYNAMITE_DEFINED_cDNADB
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  show_Hscore_cDNADB(hs,ofp)
 *
 * Descrip:    shows the Hscore by the cDNADB information
 *
 *
 *
 * Arg:         hs [UNKN ] High Score structure [Hscore *]
 * Arg:        ofp [UNKN ] output file [FILE *]
 *
 */
void Wise2_show_Hscore_cDNADB(Hscore * hs,FILE * ofp);
#define show_Hscore_cDNADB Wise2_show_Hscore_cDNADB


/* Function:  get_cDNA_from_cDNADB(cdnadb,de)
 *
 * Descrip:    Gets cDNA sequence out from
 *             the cDNADB using the information stored in
 *             dataentry
 *
 *
 * Arg:        cdnadb [READ ] cDNA database [cDNADB *]
 * Arg:            de [READ ] DataEntry information  [DataEntry *]
 *
 * Return [UNKN ]  Undocumented return value [cDNA *]
 *
 */
cDNA * Wise2_get_cDNA_from_cDNADB(cDNADB * cdnadb,DataEntry * de);
#define get_cDNA_from_cDNADB Wise2_get_cDNA_from_cDNADB


/* Function:  dataentry_add_cDNADB(de,cs,cdnadb)
 *
 * Descrip:    adds information to dataentry from cDNADB
 *
 *             will eventually add file offset and format information,
 *             but this is handled by the SequenceDB mainly.
 *
 *
 * Arg:            de [UNKN ] Undocumented argument [DataEntry *]
 * Arg:            cs [UNKN ] Undocumented argument [ComplexSequence *]
 * Arg:        cdnadb [UNKN ] Undocumented argument [cDNADB *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_dataentry_add_cDNADB(DataEntry * de,ComplexSequence * cs,cDNADB * cdnadb);
#define dataentry_add_cDNADB Wise2_dataentry_add_cDNADB


/* Function:  init_cDNADB(cdnadb,return_status)
 *
 * Descrip:    top level function which opens the cDNA database
 *
 *
 * Arg:               cdnadb [UNKN ] protein database [cDNADB *]
 * Arg:        return_status [WRITE] the status of the open from database.h [int *]
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequence *]
 *
 */
ComplexSequence * Wise2_init_cDNADB(cDNADB * cdnadb,int * return_status);
#define init_cDNADB Wise2_init_cDNADB


/* Function:  reload_cDNADB(last,cdnadb,return_status)
 *
 * Descrip:    function which reloads the database
 *
 *
 * Arg:                 last [UNKN ] previous complex sequence, will be freed [ComplexSequence *]
 * Arg:               cdnadb [UNKN ] Undocumented argument [cDNADB *]
 * Arg:        return_status [WRITE] return_status of the load [int *]
 *
 * Return [OWNER]  a new ComplexSequence object [ComplexSequence *]
 *
 */
ComplexSequence * Wise2_reload_cDNADB(ComplexSequence * last,cDNADB * cdnadb,int * return_status);
#define reload_cDNADB Wise2_reload_cDNADB


/* Function:  close_cDNADB(cs,cdnadb)
 *
 * Descrip:    top level function which closes the cDNA database
 *
 *
 * Arg:            cs [UNKN ] last complex sequence  [ComplexSequence *]
 * Arg:        cdnadb [UNKN ] protein database [cDNADB *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_close_cDNADB(ComplexSequence * cs,cDNADB * cdnadb) ;
#define close_cDNADB Wise2_close_cDNADB


/* Function:  new_cDNADB_from_single_seq(seq)
 *
 * Descrip:    To make a new cDNA database
 *             from a single cDNA Sequence with a eval system
 *
 *
 * Arg:        seq [UNKN ] sequence which as placed into cDNADB structure. [cDNA *]
 *
 * Return [UNKN ]  Undocumented return value [cDNADB *]
 *
 */
cDNADB * Wise2_new_cDNADB_from_single_seq(cDNA * seq);
#define new_cDNADB_from_single_seq Wise2_new_cDNADB_from_single_seq


/* Function:  new_cDNADB(seqdb)
 *
 * Descrip:    To make a new cDNA database
 *
 *
 * Arg:        seqdb [READ ] sequence database [SequenceDB *]
 *
 * Return [UNKN ]  Undocumented return value [cDNADB *]
 *
 */
cDNADB * Wise2_new_cDNADB(SequenceDB * seqdb);
#define new_cDNADB Wise2_new_cDNADB


/* Function:  hard_link_cDNADB(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [cDNADB *]
 *
 * Return [UNKN ]  Undocumented return value [cDNADB *]
 *
 */
cDNADB * Wise2_hard_link_cDNADB(cDNADB * obj);
#define hard_link_cDNADB Wise2_hard_link_cDNADB


/* Function:  cDNADB_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [cDNADB *]
 *
 */
cDNADB * Wise2_cDNADB_alloc(void);
#define cDNADB_alloc Wise2_cDNADB_alloc


/* Function:  free_cDNADB(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [cDNADB *]
 *
 * Return [UNKN ]  Undocumented return value [cDNADB *]
 *
 */
cDNADB * Wise2_free_cDNADB(cDNADB * obj);
#define free_cDNADB Wise2_free_cDNADB


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
boolean Wise2_access_is_single_seq_cDNADB(cDNADB * obj);
#define access_is_single_seq_cDNADB Wise2_access_is_single_seq_cDNADB
ComplexSequence * Wise2_access_forw_cDNADB(cDNADB * obj);
#define access_forw_cDNADB Wise2_access_forw_cDNADB
boolean Wise2_replace_rev_cDNADB(cDNADB * obj,ComplexSequence * rev);
#define replace_rev_cDNADB Wise2_replace_rev_cDNADB
boolean Wise2_replace_done_forward_cDNADB(cDNADB * obj,boolean done_forward);
#define replace_done_forward_cDNADB Wise2_replace_done_forward_cDNADB
ComplexSequenceEvalSet * Wise2_access_cses_cDNADB(cDNADB * obj);
#define access_cses_cDNADB Wise2_access_cses_cDNADB
ComplexSequence * Wise2_access_rev_cDNADB(cDNADB * obj);
#define access_rev_cDNADB Wise2_access_rev_cDNADB
boolean Wise2_replace_forward_only_cDNADB(cDNADB * obj,boolean forward_only);
#define replace_forward_only_cDNADB Wise2_replace_forward_only_cDNADB
boolean Wise2_replace_sdb_cDNADB(cDNADB * obj,SequenceDB * sdb);
#define replace_sdb_cDNADB Wise2_replace_sdb_cDNADB
boolean Wise2_replace_forw_cDNADB(cDNADB * obj,ComplexSequence * forw);
#define replace_forw_cDNADB Wise2_replace_forw_cDNADB
SequenceDB * Wise2_access_sdb_cDNADB(cDNADB * obj);
#define access_sdb_cDNADB Wise2_access_sdb_cDNADB
boolean Wise2_access_done_forward_cDNADB(cDNADB * obj);
#define access_done_forward_cDNADB Wise2_access_done_forward_cDNADB
boolean Wise2_replace_current_cDNADB(cDNADB * obj,Sequence * current);
#define replace_current_cDNADB Wise2_replace_current_cDNADB
boolean Wise2_replace_is_single_seq_cDNADB(cDNADB * obj,boolean is_single_seq);
#define replace_is_single_seq_cDNADB Wise2_replace_is_single_seq_cDNADB
Sequence * Wise2_access_current_cDNADB(cDNADB * obj);
#define access_current_cDNADB Wise2_access_current_cDNADB
boolean Wise2_access_forward_only_cDNADB(cDNADB * obj);
#define access_forward_only_cDNADB Wise2_access_forward_only_cDNADB
boolean Wise2_replace_cses_cDNADB(cDNADB * obj,ComplexSequenceEvalSet * cses);
#define replace_cses_cDNADB Wise2_replace_cses_cDNADB
cDNADB * Wise2_new_cDNADB_from_forrev_cseq(ComplexSequence * cs,ComplexSequence * cs_rev);
#define new_cDNADB_from_forrev_cseq Wise2_new_cDNADB_from_forrev_cseq

#ifdef _cplusplus
}
#endif

#endif
