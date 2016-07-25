

/* Functions that create, manipulate or act on Protein
 *
 * Wise2_Protein_from_Sequence
 * Wise2_hard_link_Protein
 * Wise2_Protein_alloc
 * Wise2_replace_baseseq_Protein
 * Wise2_access_baseseq_Protein
 * Wise2_free_Protein [destructor]
 *
 */

/* API for object Protein */
/* Function:  Wise2_Protein_from_Sequence(seq)
 *
 * Descrip:    makes a new protein from a Sequence. It 
 *             owns the Sequence memory, ie will attempt a /free_Sequence
 *             on the structure when /free_Protein is called
 *
 *             If you want to give this protein this Sequence and
 *             forget about it, then just hand it this sequence and set
 *             seq to NULL (no need to free it). If you intend to use 
 *             the sequecne object elsewhere outside of the Protein datastructure
 *             then use Protein_from_Sequence(/hard_link_Sequence(seq))
 *
 *
 *
 * Arg:        seq          Sequence to make protein from [Wise2_Sequence *]
 *
 * Returns Undocumented return value [Wise2_Protein *]
 *
 */
Wise2_Protein * Wise2_Protein_from_Sequence( Wise2_Sequence * seq);

/* Function:  Wise2_hard_link_Protein(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_Protein *]
 *
 * Returns Undocumented return value [Wise2_Protein *]
 *
 */
Wise2_Protein * Wise2_hard_link_Protein( Wise2_Protein * obj);

/* Function:  Wise2_Protein_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Returns Undocumented return value [Wise2_Protein *]
 *
 */
Wise2_Protein * Wise2_Protein_alloc();

/* Function:  Wise2_replace_baseseq_Protein(obj,baseseq)
 *
 * Descrip:    Replace member variable baseseq
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Protein *]
 * Arg:        baseseq      New value of the variable [Wise2_Sequence *]
 *
 * Returns member variable baseseq [boolean]
 *
 */
boolean Wise2_replace_baseseq_Protein( Wise2_Protein * obj,Wise2_Sequence * baseseq);

/* Function:  Wise2_access_baseseq_Protein(obj)
 *
 * Descrip:    Access member variable baseseq
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Protein *]
 *
 * Returns member variable baseseq [Wise2_Sequence *]
 *
 */
Wise2_Sequence * Wise2_access_baseseq_Protein( Wise2_Protein * obj);

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_Protein(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_Protein *]
 *
 * Returns Undocumented return value [Wise2_Protein *]
 *
 */
Wise2_Protein * Wise2_free_Protein( Wise2_Protein * obj);

