#ifndef DYNAMITEwise2xhmmer2HEADERFILE
#define DYNAMITEwise2xhmmer2HEADERFILE
#ifdef _cplusplus
extern "C" {
#endif

  /* include stdio etc */
#include "wisebase.h"

  /* include hmmer2 files */

  /* we have to prevent some namespace clobbering between 
     different versions of histogram */

#ifndef NO_HMMER_INCLUDES
#include "structs.h"
#include "funcs.h"
#endif

  /* include dynamite files */

  /*#include "dyna.h"*/
#include "seqalign.h"
#include "threestatemodel.h"
#include "threestatedb.h"

  /* quieten down gcc about struct p7 */

struct plan7_s;

#ifndef DYNAMITE_DEFINED_ThreeStateDB
typedef struct Wise2_ThreeStateDB Wise2_ThreeStateDB;
#define ThreeStateDB Wise2_ThreeStateDB
#define DYNAMITE_DEFINED_ThreeStateDB
#endif

#ifndef DYNAMITE_DEFINED_SeqAlign
typedef struct Wise2_SeqAlign Wise2_SeqAlign;
#define SeqAlign Wise2_SeqAlign
#define DYNAMITE_DEFINED_SeqAlign
#endif



    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  bootstrap_HMMer2(void)
 *
 * Descrip:    You need to call this function before you
 *             call any other to get the HMMer2 system working ;)
 *
 *
 *
 */
void Wise2_bootstrap_HMMer2(void);
#define bootstrap_HMMer2 Wise2_bootstrap_HMMer2


/* Function:  read_SeqAlign_HMMER(seqfile)
 *
 * Descrip:    This function read a single SeqAlign from
 *             the filename using the HMMER libraries
 *
 *
 * Arg:        seqfile [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [SeqAlign *]
 *
 */
SeqAlign * Wise2_read_SeqAlign_HMMER(char * seqfile);
#define read_SeqAlign_HMMER Wise2_read_SeqAlign_HMMER


/* Function:  HMMer2_read_ThreeStateModel(filename)
 *
 * Descrip:    This function reads a single HMM from
 *             filename for a ThreeStateModel
 *
 *
 * Arg:        filename [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateModel *]
 *
 */
ThreeStateModel * Wise2_HMMer2_read_ThreeStateModel(char * filename);
#define HMMer2_read_ThreeStateModel Wise2_HMMer2_read_ThreeStateModel


/* Function:  HMMer2_ThreeStateDB(filename)
 *
 * Descrip:    This function makes a new ThreeStateDB
 *             from a filename for the HMMer2 system.
 *
 *
 * Arg:        filename [SOFT ] filename of the HMMer db [char *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateDB *]
 *
 */
ThreeStateDB * Wise2_HMMer2_ThreeStateDB(char * filename);
#define HMMer2_ThreeStateDB Wise2_HMMer2_ThreeStateDB


/* Function:  ThreeStateModel_from_HMMer2(p7)
 *
 * Descrip:    This function converts a HMMer2 HMM into
 *             a Wise2 HMM (threestatemodel). 
 *
 *
 * Arg:        p7 [UNKN ] Undocumented argument [struct plan7_s *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateModel *]
 *
 */
ThreeStateModel * Wise2_ThreeStateModel_from_HMMer2(struct plan7_s * p7);
#define ThreeStateModel_from_HMMer2 Wise2_ThreeStateModel_from_HMMer2


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
boolean Wise2_open_for_indexing_ThreeStateDB_HMMer2(ThreeStateDB * tdb);
#define open_for_indexing_ThreeStateDB_HMMer2 Wise2_open_for_indexing_ThreeStateDB_HMMer2
boolean Wise2_close_for_indexing_ThreeStateDB_HMMer2(ThreeStateDB * tdb);
#define close_for_indexing_ThreeStateDB_HMMer2 Wise2_close_for_indexing_ThreeStateDB_HMMer2
ThreeStateModel* Wise2_index_ThreeStateDB_HMMer2(ThreeStateDB * tdb,DataEntry *de);
#define index_ThreeStateDB_HMMer2 Wise2_index_ThreeStateDB_HMMer2
boolean Wise2_dataentry_ThreeStateDB_HMMer2(ThreeStateDB * tdb,DataEntry * en);
#define dataentry_ThreeStateDB_HMMer2 Wise2_dataentry_ThreeStateDB_HMMer2
boolean Wise2_close_ThreeStateDB_HMMer2(ThreeStateDB * tdb);
#define close_ThreeStateDB_HMMer2 Wise2_close_ThreeStateDB_HMMer2
ThreeStateModel * Wise2_reload_ThreeStateDB_HMMer2(ThreeStateDB * tdb,int * return_status);
#define reload_ThreeStateDB_HMMer2 Wise2_reload_ThreeStateDB_HMMer2
boolean Wise2_open_ThreeStateDB_HMMer2(ThreeStateDB * tdb) ;
#define open_ThreeStateDB_HMMer2 Wise2_open_ThreeStateDB_HMMer2

#ifdef _cplusplus
}
#endif

#endif
