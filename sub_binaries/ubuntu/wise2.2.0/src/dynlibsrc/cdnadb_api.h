

/* Functions that create, manipulate or act on cDNADB
 *
 * Wise2_get_cDNA_from_cDNADB
 * Wise2_hard_link_cDNADB
 * Wise2_cDNADB_alloc
 * Wise2_replace_is_single_seq_cDNADB
 * Wise2_access_is_single_seq_cDNADB
 * Wise2_replace_done_forward_cDNADB
 * Wise2_access_done_forward_cDNADB
 * Wise2_replace_forward_only_cDNADB
 * Wise2_access_forward_only_cDNADB
 * Wise2_replace_forw_cDNADB
 * Wise2_access_forw_cDNADB
 * Wise2_replace_rev_cDNADB
 * Wise2_access_rev_cDNADB
 * Wise2_replace_sdb_cDNADB
 * Wise2_access_sdb_cDNADB
 * Wise2_replace_current_cDNADB
 * Wise2_access_current_cDNADB
 * Wise2_replace_cses_cDNADB
 * Wise2_access_cses_cDNADB
 * Wise2_free_cDNADB [destructor]
 *
 */



/* Helper functions in the module
 *
 * Wise2_new_cDNADB_from_single_seq
 * Wise2_new_cDNADB
 *

/* API for object cDNADB */
/* Function:  Wise2_get_cDNA_from_cDNADB(cdnadb,de)
 *
 * Descrip:    Gets cDNA sequence out from
 *             the cDNADB using the information stored in
 *             dataentry
 *
 *
 * Arg:        cdnadb       cDNA database [Wise2_cDNADB *]
 * Arg:        de           DataEntry information  [Wise2_DataEntry *]
 *
 * Returns Undocumented return value [Wise2_cDNA *]
 *
 */
Wise2_cDNA * Wise2_get_cDNA_from_cDNADB( Wise2_cDNADB * cdnadb,Wise2_DataEntry * de);

/* Function:  Wise2_hard_link_cDNADB(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_cDNADB *]
 *
 * Returns Undocumented return value [Wise2_cDNADB *]
 *
 */
Wise2_cDNADB * Wise2_hard_link_cDNADB( Wise2_cDNADB * obj);

/* Function:  Wise2_cDNADB_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Returns Undocumented return value [Wise2_cDNADB *]
 *
 */
Wise2_cDNADB * Wise2_cDNADB_alloc();

/* Function:  Wise2_replace_is_single_seq_cDNADB(obj,is_single_seq)
 *
 * Descrip:    Replace member variable is_single_seq
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_cDNADB *]
 * Arg:        is_single_seq New value of the variable [boolean]
 *
 * Returns member variable is_single_seq [boolean]
 *
 */
boolean Wise2_replace_is_single_seq_cDNADB( Wise2_cDNADB * obj,boolean is_single_seq);

/* Function:  Wise2_access_is_single_seq_cDNADB(obj)
 *
 * Descrip:    Access member variable is_single_seq
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_cDNADB *]
 *
 * Returns member variable is_single_seq [boolean]
 *
 */
boolean Wise2_access_is_single_seq_cDNADB( Wise2_cDNADB * obj);

/* Function:  Wise2_replace_done_forward_cDNADB(obj,done_forward)
 *
 * Descrip:    Replace member variable done_forward
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_cDNADB *]
 * Arg:        done_forward New value of the variable [boolean]
 *
 * Returns member variable done_forward [boolean]
 *
 */
boolean Wise2_replace_done_forward_cDNADB( Wise2_cDNADB * obj,boolean done_forward);

/* Function:  Wise2_access_done_forward_cDNADB(obj)
 *
 * Descrip:    Access member variable done_forward
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_cDNADB *]
 *
 * Returns member variable done_forward [boolean]
 *
 */
boolean Wise2_access_done_forward_cDNADB( Wise2_cDNADB * obj);

/* Function:  Wise2_replace_forward_only_cDNADB(obj,forward_only)
 *
 * Descrip:    Replace member variable forward_only
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_cDNADB *]
 * Arg:        forward_only New value of the variable [boolean]
 *
 * Returns member variable forward_only [boolean]
 *
 */
boolean Wise2_replace_forward_only_cDNADB( Wise2_cDNADB * obj,boolean forward_only);

/* Function:  Wise2_access_forward_only_cDNADB(obj)
 *
 * Descrip:    Access member variable forward_only
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_cDNADB *]
 *
 * Returns member variable forward_only [boolean]
 *
 */
boolean Wise2_access_forward_only_cDNADB( Wise2_cDNADB * obj);

/* Function:  Wise2_replace_forw_cDNADB(obj,forw)
 *
 * Descrip:    Replace member variable forw
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_cDNADB *]
 * Arg:        forw         New value of the variable [Wise2_ComplexSequence *]
 *
 * Returns member variable forw [boolean]
 *
 */
boolean Wise2_replace_forw_cDNADB( Wise2_cDNADB * obj,Wise2_ComplexSequence * forw);

/* Function:  Wise2_access_forw_cDNADB(obj)
 *
 * Descrip:    Access member variable forw
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_cDNADB *]
 *
 * Returns member variable forw [Wise2_ComplexSequence *]
 *
 */
Wise2_ComplexSequence * Wise2_access_forw_cDNADB( Wise2_cDNADB * obj);

/* Function:  Wise2_replace_rev_cDNADB(obj,rev)
 *
 * Descrip:    Replace member variable rev
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_cDNADB *]
 * Arg:        rev          New value of the variable [Wise2_ComplexSequence *]
 *
 * Returns member variable rev [boolean]
 *
 */
boolean Wise2_replace_rev_cDNADB( Wise2_cDNADB * obj,Wise2_ComplexSequence * rev);

/* Function:  Wise2_access_rev_cDNADB(obj)
 *
 * Descrip:    Access member variable rev
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_cDNADB *]
 *
 * Returns member variable rev [Wise2_ComplexSequence *]
 *
 */
Wise2_ComplexSequence * Wise2_access_rev_cDNADB( Wise2_cDNADB * obj);

/* Function:  Wise2_replace_sdb_cDNADB(obj,sdb)
 *
 * Descrip:    Replace member variable sdb
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_cDNADB *]
 * Arg:        sdb          New value of the variable [Wise2_SequenceDB *]
 *
 * Returns member variable sdb [boolean]
 *
 */
boolean Wise2_replace_sdb_cDNADB( Wise2_cDNADB * obj,Wise2_SequenceDB * sdb);

/* Function:  Wise2_access_sdb_cDNADB(obj)
 *
 * Descrip:    Access member variable sdb
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_cDNADB *]
 *
 * Returns member variable sdb [Wise2_SequenceDB *]
 *
 */
Wise2_SequenceDB * Wise2_access_sdb_cDNADB( Wise2_cDNADB * obj);

/* Function:  Wise2_replace_current_cDNADB(obj,current)
 *
 * Descrip:    Replace member variable current
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_cDNADB *]
 * Arg:        current      New value of the variable [Wise2_Sequence *]
 *
 * Returns member variable current [boolean]
 *
 */
boolean Wise2_replace_current_cDNADB( Wise2_cDNADB * obj,Wise2_Sequence * current);

/* Function:  Wise2_access_current_cDNADB(obj)
 *
 * Descrip:    Access member variable current
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_cDNADB *]
 *
 * Returns member variable current [Wise2_Sequence *]
 *
 */
Wise2_Sequence * Wise2_access_current_cDNADB( Wise2_cDNADB * obj);

/* Function:  Wise2_replace_cses_cDNADB(obj,cses)
 *
 * Descrip:    Replace member variable cses
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_cDNADB *]
 * Arg:        cses         New value of the variable [Wise2_ComplexSequenceEvalSet *]
 *
 * Returns member variable cses [boolean]
 *
 */
boolean Wise2_replace_cses_cDNADB( Wise2_cDNADB * obj,Wise2_ComplexSequenceEvalSet * cses);

/* Function:  Wise2_access_cses_cDNADB(obj)
 *
 * Descrip:    Access member variable cses
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_cDNADB *]
 *
 * Returns member variable cses [Wise2_ComplexSequenceEvalSet *]
 *
 */
Wise2_ComplexSequenceEvalSet * Wise2_access_cses_cDNADB( Wise2_cDNADB * obj);

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_cDNADB(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_cDNADB *]
 *
 * Returns Undocumented return value [Wise2_cDNADB *]
 *
 */
Wise2_cDNADB * Wise2_free_cDNADB( Wise2_cDNADB * obj);



/* These functions are not associated with an object */
/* Function:  Wise2_new_cDNADB_from_single_seq(seq)
 *
 * Descrip:    To make a new cDNA database
 *             from a single cDNA Sequence with a eval system
 *
 *
 * Arg:        seq          sequence which as placed into cDNADB structure. [Wise2_cDNA *]
 *
 * Returns Undocumented return value [Wise2_cDNADB *]
 *
 */
Wise2_cDNADB * Wise2_new_cDNADB_from_single_seq( Wise2_cDNA * seq);

/* Function:  Wise2_new_cDNADB(seqdb)
 *
 * Descrip:    To make a new cDNA database
 *
 *
 * Arg:        seqdb        sequence database [Wise2_SequenceDB *]
 *
 * Returns Undocumented return value [Wise2_cDNADB *]
 *
 */
Wise2_cDNADB * Wise2_new_cDNADB( Wise2_SequenceDB * seqdb);

