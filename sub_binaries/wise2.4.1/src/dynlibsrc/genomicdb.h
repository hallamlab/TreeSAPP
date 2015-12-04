#ifndef DYNAMITEgenomicdbHEADERFILE
#define DYNAMITEgenomicdbHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "sequencedb.h"
#include "genomic.h"
#include "hscore.h"
#include "complexsequence.h"
#include "complexevalset.h"

typedef enum GenDBErrorType {
  GENDB_READ_THROUGH = 0,
  GENDB_FAIL_ON_ERROR = 1
} GenDBErrorType;

/* Object GenomicDB
 *
 * Descrip: This object hold a database of
 *        genomic sequences.
 *
 *        You will probably use it in one of
 *        two ways
 *
 *        1 A sequence formatted database, which
 *        is provided by a /SequenceDB object
 *        is used to provide the raw sequences 
 *
 *        2 A single Genomic sequence is used.
 *
 *        In each case this database provides
 *        both the forward and reverse strands
 *        into the system.
 *
 *        Notice that what is exported are
 *        /ComplexSequence objects, not genomic dna,
 *        as this is what is generally needed. 
 *        These are things with splice sites calculated
 *        etc. This is why for initialisation this needs
 *        a /ComplexSequenceEvalSet of the correct type.
 *
 *
 */
struct Wise2_GenomicDB {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    boolean is_single_seq;   
    boolean done_forward;    
    ComplexSequence * forw;  
    ComplexSequence * rev;   
    SequenceDB * sdb;    
    Genomic * current;   
    ComplexSequenceEvalSet * cses;   
    GenDBErrorType error_handling;   
    Genomic * single;   /*  for single sequence cases, so we can 'index' on it  */ 
    Genomic * revsingle;     
    int length_of_N;     
    int repeat_in_cds_score;     
    } ;  
/* GenomicDB defined */ 
#ifndef DYNAMITE_DEFINED_GenomicDB
typedef struct Wise2_GenomicDB Wise2_GenomicDB;
#define GenomicDB Wise2_GenomicDB
#define DYNAMITE_DEFINED_GenomicDB
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  show_Hscore_GenomicDB(hs,ofp)
 *
 * Descrip:    shows the Hscore by the GenomicDB information
 *
 *
 *
 * Arg:         hs [UNKN ] High Score structure [Hscore *]
 * Arg:        ofp [UNKN ] output file [FILE *]
 *
 */
void Wise2_show_Hscore_GenomicDB(Hscore * hs,FILE * ofp);
#define show_Hscore_GenomicDB Wise2_show_Hscore_GenomicDB


/* Function:  get_Genomic_from_GenomicDB(gendb,de)
 *
 * Descrip:    Gets Genomic sequence out from
 *             the GenomicDB using the information stored in
 *             dataentry
 *
 *
 * Arg:        gendb [UNKN ] Undocumented argument [GenomicDB *]
 * Arg:           de [UNKN ] Undocumented argument [DataEntry *]
 *
 * Return [UNKN ]  Undocumented return value [Genomic *]
 *
 */
Genomic * Wise2_get_Genomic_from_GenomicDB(GenomicDB * gendb,DataEntry * de);
#define get_Genomic_from_GenomicDB Wise2_get_Genomic_from_GenomicDB


/* Function:  dataentry_add_GenomicDB(de,cs,gendb)
 *
 * Descrip:    adds information to dataentry from GenomicDB
 *
 *             will eventually add file offset and format information
 *
 *
 * Arg:           de [UNKN ] Undocumented argument [DataEntry *]
 * Arg:           cs [UNKN ] Undocumented argument [ComplexSequence *]
 * Arg:        gendb [UNKN ] Undocumented argument [GenomicDB *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_dataentry_add_GenomicDB(DataEntry * de,ComplexSequence * cs,GenomicDB * gendb);
#define dataentry_add_GenomicDB Wise2_dataentry_add_GenomicDB


/* Function:  init_GenomicDB(gendb,return_status)
 *
 * Descrip:    top level function which opens the Genomic database
 *
 *
 * Arg:                gendb [UNKN ] protein database [GenomicDB *]
 * Arg:        return_status [WRITE] the status of the open from database.h [int *]
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequence *]
 *
 */
ComplexSequence * Wise2_init_GenomicDB(GenomicDB * gendb,int * return_status);
#define init_GenomicDB Wise2_init_GenomicDB


/* Function:  reload_GenomicDB(last,gendb,return_status)
 *
 * Descrip:    function which reloads the database
 *
 *
 * Arg:                 last [UNKN ] previous complex sequence, will be freed [ComplexSequence *]
 * Arg:                gendb [UNKN ] Undocumented argument [GenomicDB *]
 * Arg:        return_status [WRITE] return_status of the load [int *]
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequence *]
 *
 */
ComplexSequence * Wise2_reload_GenomicDB(ComplexSequence * last,GenomicDB * gendb,int * return_status);
#define reload_GenomicDB Wise2_reload_GenomicDB


/* Function:  close_GenomicDB(cs,gendb)
 *
 * Descrip:    top level function which closes the genomic database
 *
 *
 * Arg:           cs [UNKN ] last complex sequence  [ComplexSequence *]
 * Arg:        gendb [UNKN ] protein database [GenomicDB *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_close_GenomicDB(ComplexSequence * cs,GenomicDB * gendb) ;
#define close_GenomicDB Wise2_close_GenomicDB


/* Function:  new_GenomicDB_from_single_seq(gen,cses,score_in_repeat_coding)
 *
 * Descrip:    To make a new genomic database
 *             from a single Genomic Sequence with a eval system
 *
 *
 * Arg:                           gen [UNKN ] sequence which as placed into GenomicDB structure. [Genomic *]
 * Arg:                          cses [UNKN ] Undocumented argument [ComplexSequenceEvalSet *]
 * Arg:        score_in_repeat_coding [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [GenomicDB *]
 *
 */
GenomicDB * Wise2_new_GenomicDB_from_single_seq(Genomic * gen,ComplexSequenceEvalSet * cses,int score_in_repeat_coding);
#define new_GenomicDB_from_single_seq Wise2_new_GenomicDB_from_single_seq


/* Function:  new_GenomicDB_from_forrev_cseq(cs,cs_rev)
 *
 * Descrip:    To make a new genomic database
 *             from a single ComplexSequence
 *
 *
 * Arg:            cs [UNKN ] complex sequence which is held. [ComplexSequence *]
 * Arg:        cs_rev [UNKN ] Undocumented argument [ComplexSequence *]
 *
 * Return [UNKN ]  Undocumented return value [GenomicDB *]
 *
 */
GenomicDB * Wise2_new_GenomicDB_from_forrev_cseq(ComplexSequence * cs,ComplexSequence * cs_rev);
#define new_GenomicDB_from_forrev_cseq Wise2_new_GenomicDB_from_forrev_cseq


/* Function:  new_GenomicDB(seqdb,cses,length_of_N,repeat_in_cds_score)
 *
 * Descrip:    To make a new genomic database
 *
 *
 * Arg:                      seqdb [UNKN ] sequence database [SequenceDB *]
 * Arg:                       cses [UNKN ] protein evaluation set [ComplexSequenceEvalSet *]
 * Arg:                length_of_N [UNKN ] Undocumented argument [int]
 * Arg:        repeat_in_cds_score [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [GenomicDB *]
 *
 */
GenomicDB * Wise2_new_GenomicDB(SequenceDB * seqdb,ComplexSequenceEvalSet * cses,int length_of_N,int repeat_in_cds_score);
#define new_GenomicDB Wise2_new_GenomicDB


/* Function:  hard_link_GenomicDB(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GenomicDB *]
 *
 * Return [UNKN ]  Undocumented return value [GenomicDB *]
 *
 */
GenomicDB * Wise2_hard_link_GenomicDB(GenomicDB * obj);
#define hard_link_GenomicDB Wise2_hard_link_GenomicDB


/* Function:  GenomicDB_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GenomicDB *]
 *
 */
GenomicDB * Wise2_GenomicDB_alloc(void);
#define GenomicDB_alloc Wise2_GenomicDB_alloc


/* Function:  free_GenomicDB(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GenomicDB *]
 *
 * Return [UNKN ]  Undocumented return value [GenomicDB *]
 *
 */
GenomicDB * Wise2_free_GenomicDB(GenomicDB * obj);
#define free_GenomicDB Wise2_free_GenomicDB


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
boolean Wise2_replace_is_single_seq_GenomicDB(GenomicDB * obj,boolean is_single_seq);
#define replace_is_single_seq_GenomicDB Wise2_replace_is_single_seq_GenomicDB
ComplexSequence * Wise2_access_rev_GenomicDB(GenomicDB * obj);
#define access_rev_GenomicDB Wise2_access_rev_GenomicDB
boolean Wise2_access_is_single_seq_GenomicDB(GenomicDB * obj);
#define access_is_single_seq_GenomicDB Wise2_access_is_single_seq_GenomicDB
boolean Wise2_replace_sdb_GenomicDB(GenomicDB * obj,SequenceDB * sdb);
#define replace_sdb_GenomicDB Wise2_replace_sdb_GenomicDB
boolean Wise2_access_done_forward_GenomicDB(GenomicDB * obj);
#define access_done_forward_GenomicDB Wise2_access_done_forward_GenomicDB
SequenceDB * Wise2_access_sdb_GenomicDB(GenomicDB * obj);
#define access_sdb_GenomicDB Wise2_access_sdb_GenomicDB
ComplexSequence * Wise2_access_forw_GenomicDB(GenomicDB * obj);
#define access_forw_GenomicDB Wise2_access_forw_GenomicDB
boolean Wise2_replace_current_GenomicDB(GenomicDB * obj,Genomic * current);
#define replace_current_GenomicDB Wise2_replace_current_GenomicDB
boolean Wise2_replace_done_forward_GenomicDB(GenomicDB * obj,boolean done_forward);
#define replace_done_forward_GenomicDB Wise2_replace_done_forward_GenomicDB
Genomic * Wise2_access_current_GenomicDB(GenomicDB * obj);
#define access_current_GenomicDB Wise2_access_current_GenomicDB
boolean Wise2_replace_rev_GenomicDB(GenomicDB * obj,ComplexSequence * rev);
#define replace_rev_GenomicDB Wise2_replace_rev_GenomicDB
boolean Wise2_replace_cses_GenomicDB(GenomicDB * obj,ComplexSequenceEvalSet * cses);
#define replace_cses_GenomicDB Wise2_replace_cses_GenomicDB
boolean Wise2_replace_forw_GenomicDB(GenomicDB * obj,ComplexSequence * forw);
#define replace_forw_GenomicDB Wise2_replace_forw_GenomicDB
ComplexSequenceEvalSet * Wise2_access_cses_GenomicDB(GenomicDB * obj);
#define access_cses_GenomicDB Wise2_access_cses_GenomicDB

#ifdef _cplusplus
}
#endif

#endif
