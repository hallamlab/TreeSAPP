

/* Functions that create, manipulate or act on GenomicDB
 *
 * Wise2_get_Genomic_from_GenomicDB
 * Wise2_hard_link_GenomicDB
 * Wise2_GenomicDB_alloc
 * Wise2_replace_is_single_seq_GenomicDB
 * Wise2_access_is_single_seq_GenomicDB
 * Wise2_replace_done_forward_GenomicDB
 * Wise2_access_done_forward_GenomicDB
 * Wise2_replace_forw_GenomicDB
 * Wise2_access_forw_GenomicDB
 * Wise2_replace_rev_GenomicDB
 * Wise2_access_rev_GenomicDB
 * Wise2_replace_sdb_GenomicDB
 * Wise2_access_sdb_GenomicDB
 * Wise2_replace_current_GenomicDB
 * Wise2_access_current_GenomicDB
 * Wise2_replace_cses_GenomicDB
 * Wise2_access_cses_GenomicDB
 * Wise2_free_GenomicDB [destructor]
 *
 */



/* Helper functions in the module
 *
 * Wise2_new_GenomicDB_from_single_seq
 * Wise2_new_GenomicDB
 *

/* API for object GenomicDB */
/* Function:  Wise2_get_Genomic_from_GenomicDB(gendb,de)
 *
 * Descrip:    Gets Genomic sequence out from
 *             the GenomicDB using the information stored in
 *             dataentry
 *
 *
 * Arg:        gendb        Undocumented argument [Wise2_GenomicDB *]
 * Arg:        de           Undocumented argument [Wise2_DataEntry *]
 *
 * Returns Undocumented return value [Wise2_Genomic *]
 *
 */
Wise2_Genomic * Wise2_get_Genomic_from_GenomicDB( Wise2_GenomicDB * gendb,Wise2_DataEntry * de);

/* Function:  Wise2_hard_link_GenomicDB(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_GenomicDB *]
 *
 * Returns Undocumented return value [Wise2_GenomicDB *]
 *
 */
Wise2_GenomicDB * Wise2_hard_link_GenomicDB( Wise2_GenomicDB * obj);

/* Function:  Wise2_GenomicDB_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Returns Undocumented return value [Wise2_GenomicDB *]
 *
 */
Wise2_GenomicDB * Wise2_GenomicDB_alloc();

/* Function:  Wise2_replace_is_single_seq_GenomicDB(obj,is_single_seq)
 *
 * Descrip:    Replace member variable is_single_seq
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_GenomicDB *]
 * Arg:        is_single_seq New value of the variable [boolean]
 *
 * Returns member variable is_single_seq [boolean]
 *
 */
boolean Wise2_replace_is_single_seq_GenomicDB( Wise2_GenomicDB * obj,boolean is_single_seq);

/* Function:  Wise2_access_is_single_seq_GenomicDB(obj)
 *
 * Descrip:    Access member variable is_single_seq
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_GenomicDB *]
 *
 * Returns member variable is_single_seq [boolean]
 *
 */
boolean Wise2_access_is_single_seq_GenomicDB( Wise2_GenomicDB * obj);

/* Function:  Wise2_replace_done_forward_GenomicDB(obj,done_forward)
 *
 * Descrip:    Replace member variable done_forward
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_GenomicDB *]
 * Arg:        done_forward New value of the variable [boolean]
 *
 * Returns member variable done_forward [boolean]
 *
 */
boolean Wise2_replace_done_forward_GenomicDB( Wise2_GenomicDB * obj,boolean done_forward);

/* Function:  Wise2_access_done_forward_GenomicDB(obj)
 *
 * Descrip:    Access member variable done_forward
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_GenomicDB *]
 *
 * Returns member variable done_forward [boolean]
 *
 */
boolean Wise2_access_done_forward_GenomicDB( Wise2_GenomicDB * obj);

/* Function:  Wise2_replace_forw_GenomicDB(obj,forw)
 *
 * Descrip:    Replace member variable forw
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_GenomicDB *]
 * Arg:        forw         New value of the variable [Wise2_ComplexSequence *]
 *
 * Returns member variable forw [boolean]
 *
 */
boolean Wise2_replace_forw_GenomicDB( Wise2_GenomicDB * obj,Wise2_ComplexSequence * forw);

/* Function:  Wise2_access_forw_GenomicDB(obj)
 *
 * Descrip:    Access member variable forw
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_GenomicDB *]
 *
 * Returns member variable forw [Wise2_ComplexSequence *]
 *
 */
Wise2_ComplexSequence * Wise2_access_forw_GenomicDB( Wise2_GenomicDB * obj);

/* Function:  Wise2_replace_rev_GenomicDB(obj,rev)
 *
 * Descrip:    Replace member variable rev
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_GenomicDB *]
 * Arg:        rev          New value of the variable [Wise2_ComplexSequence *]
 *
 * Returns member variable rev [boolean]
 *
 */
boolean Wise2_replace_rev_GenomicDB( Wise2_GenomicDB * obj,Wise2_ComplexSequence * rev);

/* Function:  Wise2_access_rev_GenomicDB(obj)
 *
 * Descrip:    Access member variable rev
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_GenomicDB *]
 *
 * Returns member variable rev [Wise2_ComplexSequence *]
 *
 */
Wise2_ComplexSequence * Wise2_access_rev_GenomicDB( Wise2_GenomicDB * obj);

/* Function:  Wise2_replace_sdb_GenomicDB(obj,sdb)
 *
 * Descrip:    Replace member variable sdb
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_GenomicDB *]
 * Arg:        sdb          New value of the variable [Wise2_SequenceDB *]
 *
 * Returns member variable sdb [boolean]
 *
 */
boolean Wise2_replace_sdb_GenomicDB( Wise2_GenomicDB * obj,Wise2_SequenceDB * sdb);

/* Function:  Wise2_access_sdb_GenomicDB(obj)
 *
 * Descrip:    Access member variable sdb
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_GenomicDB *]
 *
 * Returns member variable sdb [Wise2_SequenceDB *]
 *
 */
Wise2_SequenceDB * Wise2_access_sdb_GenomicDB( Wise2_GenomicDB * obj);

/* Function:  Wise2_replace_current_GenomicDB(obj,current)
 *
 * Descrip:    Replace member variable current
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_GenomicDB *]
 * Arg:        current      New value of the variable [Wise2_Genomic *]
 *
 * Returns member variable current [boolean]
 *
 */
boolean Wise2_replace_current_GenomicDB( Wise2_GenomicDB * obj,Wise2_Genomic * current);

/* Function:  Wise2_access_current_GenomicDB(obj)
 *
 * Descrip:    Access member variable current
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_GenomicDB *]
 *
 * Returns member variable current [Wise2_Genomic *]
 *
 */
Wise2_Genomic * Wise2_access_current_GenomicDB( Wise2_GenomicDB * obj);

/* Function:  Wise2_replace_cses_GenomicDB(obj,cses)
 *
 * Descrip:    Replace member variable cses
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_GenomicDB *]
 * Arg:        cses         New value of the variable [Wise2_ComplexSequenceEvalSet *]
 *
 * Returns member variable cses [boolean]
 *
 */
boolean Wise2_replace_cses_GenomicDB( Wise2_GenomicDB * obj,Wise2_ComplexSequenceEvalSet * cses);

/* Function:  Wise2_access_cses_GenomicDB(obj)
 *
 * Descrip:    Access member variable cses
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_GenomicDB *]
 *
 * Returns member variable cses [Wise2_ComplexSequenceEvalSet *]
 *
 */
Wise2_ComplexSequenceEvalSet * Wise2_access_cses_GenomicDB( Wise2_GenomicDB * obj);

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_GenomicDB(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_GenomicDB *]
 *
 * Returns Undocumented return value [Wise2_GenomicDB *]
 *
 */
Wise2_GenomicDB * Wise2_free_GenomicDB( Wise2_GenomicDB * obj);



/* These functions are not associated with an object */
/* Function:  Wise2_new_GenomicDB_from_single_seq(gen,cses,score_in_repeat_coding)
 *
 * Descrip:    To make a new genomic database
 *             from a single Genomic Sequence with a eval system
 *
 *
 * Arg:        gen          sequence which as placed into GenomicDB structure. [Wise2_Genomic *]
 * Arg:        cses         Undocumented argument [Wise2_ComplexSequenceEvalSet *]
 * Arg:        score_in_repeat_coding Undocumented argument [int]
 *
 * Returns Undocumented return value [Wise2_GenomicDB *]
 *
 */
Wise2_GenomicDB * Wise2_new_GenomicDB_from_single_seq( Wise2_Genomic * gen,Wise2_ComplexSequenceEvalSet * cses,int score_in_repeat_coding);

/* Function:  Wise2_new_GenomicDB(seqdb,cses,length_of_N,repeat_in_cds_score)
 *
 * Descrip:    To make a new genomic database
 *
 *
 * Arg:        seqdb        sequence database [Wise2_SequenceDB *]
 * Arg:        cses         protein evaluation set [Wise2_ComplexSequenceEvalSet *]
 * Arg:        length_of_N  Undocumented argument [int]
 * Arg:        repeat_in_cds_score Undocumented argument [int]
 *
 * Returns Undocumented return value [Wise2_GenomicDB *]
 *
 */
Wise2_GenomicDB * Wise2_new_GenomicDB( Wise2_SequenceDB * seqdb,Wise2_ComplexSequenceEvalSet * cses,int length_of_N,int repeat_in_cds_score);

