

/* Functions that create, manipulate or act on ProteinDB
 *
 * Wise2_hard_link_ProteinDB
 * Wise2_ProteinDB_alloc
 * Wise2_replace_is_single_seq_ProteinDB
 * Wise2_access_is_single_seq_ProteinDB
 * Wise2_replace_is_random_db_ProteinDB
 * Wise2_access_is_random_db_ProteinDB
 * Wise2_replace_single_ProteinDB
 * Wise2_access_single_ProteinDB
 * Wise2_replace_sdb_ProteinDB
 * Wise2_access_sdb_ProteinDB
 * Wise2_replace_cses_ProteinDB
 * Wise2_access_cses_ProteinDB
 * Wise2_replace_rnd_ProteinDB
 * Wise2_access_rnd_ProteinDB
 * Wise2_replace_test_dna_ProteinDB
 * Wise2_access_test_dna_ProteinDB
 * Wise2_free_ProteinDB [destructor]
 *
 */



/* Helper functions in the module
 *
 * Wise2_new_ProteinDB_from_single_seq
 * Wise2_single_fasta_ProteinDB
 * Wise2_new_ProteinDB
 *

/* API for object ProteinDB */
/* Function:  Wise2_hard_link_ProteinDB(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_ProteinDB *]
 *
 * Returns Undocumented return value [Wise2_ProteinDB *]
 *
 */
Wise2_ProteinDB * Wise2_hard_link_ProteinDB( Wise2_ProteinDB * obj);

/* Function:  Wise2_ProteinDB_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Returns Undocumented return value [Wise2_ProteinDB *]
 *
 */
Wise2_ProteinDB * Wise2_ProteinDB_alloc();

/* Function:  Wise2_replace_is_single_seq_ProteinDB(obj,is_single_seq)
 *
 * Descrip:    Replace member variable is_single_seq
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_ProteinDB *]
 * Arg:        is_single_seq New value of the variable [boolean]
 *
 * Returns member variable is_single_seq [boolean]
 *
 */
boolean Wise2_replace_is_single_seq_ProteinDB( Wise2_ProteinDB * obj,boolean is_single_seq);

/* Function:  Wise2_access_is_single_seq_ProteinDB(obj)
 *
 * Descrip:    Access member variable is_single_seq
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_ProteinDB *]
 *
 * Returns member variable is_single_seq [boolean]
 *
 */
boolean Wise2_access_is_single_seq_ProteinDB( Wise2_ProteinDB * obj);

/* Function:  Wise2_replace_is_random_db_ProteinDB(obj,is_random_db)
 *
 * Descrip:    Replace member variable is_random_db
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_ProteinDB *]
 * Arg:        is_random_db New value of the variable [boolean]
 *
 * Returns member variable is_random_db [boolean]
 *
 */
boolean Wise2_replace_is_random_db_ProteinDB( Wise2_ProteinDB * obj,boolean is_random_db);

/* Function:  Wise2_access_is_random_db_ProteinDB(obj)
 *
 * Descrip:    Access member variable is_random_db
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_ProteinDB *]
 *
 * Returns member variable is_random_db [boolean]
 *
 */
boolean Wise2_access_is_random_db_ProteinDB( Wise2_ProteinDB * obj);

/* Function:  Wise2_replace_single_ProteinDB(obj,single)
 *
 * Descrip:    Replace member variable single
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_ProteinDB *]
 * Arg:        single       New value of the variable [Wise2_ComplexSequence *]
 *
 * Returns member variable single [boolean]
 *
 */
boolean Wise2_replace_single_ProteinDB( Wise2_ProteinDB * obj,Wise2_ComplexSequence * single);

/* Function:  Wise2_access_single_ProteinDB(obj)
 *
 * Descrip:    Access member variable single
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_ProteinDB *]
 *
 * Returns member variable single [Wise2_ComplexSequence *]
 *
 */
Wise2_ComplexSequence * Wise2_access_single_ProteinDB( Wise2_ProteinDB * obj);

/* Function:  Wise2_replace_sdb_ProteinDB(obj,sdb)
 *
 * Descrip:    Replace member variable sdb
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_ProteinDB *]
 * Arg:        sdb          New value of the variable [Wise2_SequenceDB *]
 *
 * Returns member variable sdb [boolean]
 *
 */
boolean Wise2_replace_sdb_ProteinDB( Wise2_ProteinDB * obj,Wise2_SequenceDB * sdb);

/* Function:  Wise2_access_sdb_ProteinDB(obj)
 *
 * Descrip:    Access member variable sdb
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_ProteinDB *]
 *
 * Returns member variable sdb [Wise2_SequenceDB *]
 *
 */
Wise2_SequenceDB * Wise2_access_sdb_ProteinDB( Wise2_ProteinDB * obj);

/* Function:  Wise2_replace_cses_ProteinDB(obj,cses)
 *
 * Descrip:    Replace member variable cses
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_ProteinDB *]
 * Arg:        cses         New value of the variable [Wise2_ComplexSequenceEvalSet *]
 *
 * Returns member variable cses [boolean]
 *
 */
boolean Wise2_replace_cses_ProteinDB( Wise2_ProteinDB * obj,Wise2_ComplexSequenceEvalSet * cses);

/* Function:  Wise2_access_cses_ProteinDB(obj)
 *
 * Descrip:    Access member variable cses
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_ProteinDB *]
 *
 * Returns member variable cses [Wise2_ComplexSequenceEvalSet *]
 *
 */
Wise2_ComplexSequenceEvalSet * Wise2_access_cses_ProteinDB( Wise2_ProteinDB * obj);

/* Function:  Wise2_replace_rnd_ProteinDB(obj,rnd)
 *
 * Descrip:    Replace member variable rnd
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_ProteinDB *]
 * Arg:        rnd          New value of the variable [Wise2_RandomProteinDB *]
 *
 * Returns member variable rnd [boolean]
 *
 */
boolean Wise2_replace_rnd_ProteinDB( Wise2_ProteinDB * obj,Wise2_RandomProteinDB * rnd);

/* Function:  Wise2_access_rnd_ProteinDB(obj)
 *
 * Descrip:    Access member variable rnd
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_ProteinDB *]
 *
 * Returns member variable rnd [Wise2_RandomProteinDB *]
 *
 */
Wise2_RandomProteinDB * Wise2_access_rnd_ProteinDB( Wise2_ProteinDB * obj);

/* Function:  Wise2_replace_test_dna_ProteinDB(obj,test_dna)
 *
 * Descrip:    Replace member variable test_dna
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_ProteinDB *]
 * Arg:        test_dna     New value of the variable [boolean]
 *
 * Returns member variable test_dna [boolean]
 *
 */
boolean Wise2_replace_test_dna_ProteinDB( Wise2_ProteinDB * obj,boolean test_dna);

/* Function:  Wise2_access_test_dna_ProteinDB(obj)
 *
 * Descrip:    Access member variable test_dna
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_ProteinDB *]
 *
 * Returns member variable test_dna [boolean]
 *
 */
boolean Wise2_access_test_dna_ProteinDB( Wise2_ProteinDB * obj);

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_ProteinDB(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_ProteinDB *]
 *
 * Returns Undocumented return value [Wise2_ProteinDB *]
 *
 */
Wise2_ProteinDB * Wise2_free_ProteinDB( Wise2_ProteinDB * obj);



/* These functions are not associated with an object */
/* Function:  Wise2_new_ProteinDB_from_single_seq(seq)
 *
 * Descrip:    To make a new protein database
 *             from a single Sequence with default amino acid mapping
 *
 *
 * Arg:        seq          sequence which as placed into ProteinDB structure. [Wise2_Sequence *]
 *
 * Returns Undocumented return value [Wise2_ProteinDB *]
 *
 */
Wise2_ProteinDB * Wise2_new_ProteinDB_from_single_seq( Wise2_Sequence * seq);

/* Function:  Wise2_single_fasta_ProteinDB(filename)
 *
 * Descrip:    pre-packed single fasta protein database
 *
 *
 *
 * Arg:        filename     name of fasta file [char *]
 *
 * Returns Undocumented return value [Wise2_ProteinDB *]
 *
 */
Wise2_ProteinDB * Wise2_single_fasta_ProteinDB( char * filename);

/* Function:  Wise2_new_ProteinDB(seqdb,cses)
 *
 * Descrip:    To make a new protein database
 *
 *
 * Arg:        seqdb        sequence database [Wise2_SequenceDB *]
 * Arg:        cses         protein evaluation set [Wise2_ComplexSequenceEvalSet *]
 *
 * Returns Undocumented return value [Wise2_ProteinDB *]
 *
 */
Wise2_ProteinDB * Wise2_new_ProteinDB( Wise2_SequenceDB * seqdb,Wise2_ComplexSequenceEvalSet * cses);

