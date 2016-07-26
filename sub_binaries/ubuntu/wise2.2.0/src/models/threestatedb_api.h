

/* Functions that create, manipulate or act on ThreeStateDB
 *
 * Wise2_indexed_ThreeStateModel_ThreeStateDB
 * Wise2_new_proteindb_ThreeStateDB
 * Wise2_new_PfamHmmer1DB_ThreeStateDB
 * Wise2_new_single_ThreeStateDB
 * Wise2_hard_link_ThreeStateDB
 * Wise2_ThreeStateDB_alloc
 * Wise2_replace_dbtype_ThreeStateDB
 * Wise2_access_dbtype_ThreeStateDB
 * Wise2_replace_filename_ThreeStateDB
 * Wise2_access_filename_ThreeStateDB
 * Wise2_free_ThreeStateDB [destructor]
 *
 */

/* API for object ThreeStateDB */
/* Function:  Wise2_indexed_ThreeStateModel_ThreeStateDB(mdb,en)
 *
 * Descrip:    Retrieves a model from a database which has been opened
 *             for indexing by /open_for_indexing_ThreeStateDB
 *
 *             The index information comes from the dataentry which should 
 *             have been from a search of the ThreeStateDB.
 *
 *
 * Arg:        mdb          database where this is indexed [Wise2_ThreeStateDB *]
 * Arg:        en           dataentry to pull the model from [Wise2_DataEntry *]
 *
 * Returns Undocumented return value [Wise2_ThreeStateModel *]
 *
 */
Wise2_ThreeStateModel * Wise2_indexed_ThreeStateModel_ThreeStateDB( Wise2_ThreeStateDB * mdb,Wise2_DataEntry * en);

/* Function:  Wise2_new_proteindb_ThreeStateDB(sdb,comp,gap,ext)
 *
 * Descrip:    makes a new ThreeStateDB from a
 *             sequencedb (better be protein!)
 *
 *
 *
 * Arg:        sdb          sequence database to use [Wise2_SequenceDB *]
 * Arg:        comp         comparison matrix to use [Wise2_CompMat *]
 * Arg:        gap          gap open penalty [int]
 * Arg:        ext          gap extensions penalty [int]
 *
 * Returns Undocumented return value [Wise2_ThreeStateDB *]
 *
 */
Wise2_ThreeStateDB * Wise2_new_proteindb_ThreeStateDB( Wise2_SequenceDB * sdb,Wise2_CompMat * comp,int gap,int ext);

/* Function:  Wise2_new_PfamHmmer1DB_ThreeStateDB(dirname)
 *
 * Descrip:    Makes a new PfamHmmer1DB from a filename
 *             indicating the directory
 *
 *
 * Arg:        dirname      Undocumented argument [char *]
 *
 * Returns Undocumented return value [Wise2_ThreeStateDB *]
 *
 */
Wise2_ThreeStateDB * Wise2_new_PfamHmmer1DB_ThreeStateDB( char * dirname);

/* Function:  Wise2_new_single_ThreeStateDB(tsm,rm)
 *
 * Descrip:    Making a new ThreeStateDB from a single
 *             model
 *
 *
 *
 * Arg:        tsm          a single ThreeStateModel [Wise2_ThreeStateModel *]
 * Arg:        rm           random model to be used in comparisons.. [Wise2_RandomModel *]
 *
 * Returns Undocumented return value [Wise2_ThreeStateDB *]
 *
 */
Wise2_ThreeStateDB * Wise2_new_single_ThreeStateDB( Wise2_ThreeStateModel * tsm,Wise2_RandomModel * rm);

/* Function:  Wise2_hard_link_ThreeStateDB(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_ThreeStateDB *]
 *
 * Returns Undocumented return value [Wise2_ThreeStateDB *]
 *
 */
Wise2_ThreeStateDB * Wise2_hard_link_ThreeStateDB( Wise2_ThreeStateDB * obj);

/* Function:  Wise2_ThreeStateDB_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Returns Undocumented return value [Wise2_ThreeStateDB *]
 *
 */
Wise2_ThreeStateDB * Wise2_ThreeStateDB_alloc();

/* Function:  Wise2_replace_dbtype_ThreeStateDB(obj,dbtype)
 *
 * Descrip:    Replace member variable dbtype
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_ThreeStateDB *]
 * Arg:        dbtype       New value of the variable [int]
 *
 * Returns member variable dbtype [boolean]
 *
 */
boolean Wise2_replace_dbtype_ThreeStateDB( Wise2_ThreeStateDB * obj,int dbtype);

/* Function:  Wise2_access_dbtype_ThreeStateDB(obj)
 *
 * Descrip:    Access member variable dbtype
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_ThreeStateDB *]
 *
 * Returns member variable dbtype [int]
 *
 */
int Wise2_access_dbtype_ThreeStateDB( Wise2_ThreeStateDB * obj);

/* Function:  Wise2_replace_filename_ThreeStateDB(obj,filename)
 *
 * Descrip:    Replace member variable filename
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_ThreeStateDB *]
 * Arg:        filename     New value of the variable [char *]
 *
 * Returns member variable filename [boolean]
 *
 */
boolean Wise2_replace_filename_ThreeStateDB( Wise2_ThreeStateDB * obj,char * filename);

/* Function:  Wise2_access_filename_ThreeStateDB(obj)
 *
 * Descrip:    Access member variable filename
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_ThreeStateDB *]
 *
 * Returns member variable filename [char *]
 *
 */
char * Wise2_access_filename_ThreeStateDB( Wise2_ThreeStateDB * obj);

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_ThreeStateDB(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_ThreeStateDB *]
 *
 * Returns Undocumented return value [Wise2_ThreeStateDB *]
 *
 */
Wise2_ThreeStateDB * Wise2_free_ThreeStateDB( Wise2_ThreeStateDB * obj);

