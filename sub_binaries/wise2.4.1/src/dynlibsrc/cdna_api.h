

/* Functions that create, manipulate or act on cDNA
 *
 * Wise2_truncate_cDNA
 * Wise2_read_fasta_file_cDNA
 * Wise2_cDNA_name
 * Wise2_cDNA_length
 * Wise2_cDNA_seqchar
 * Wise2_cDNA_from_Sequence
 * Wise2_hard_link_cDNA
 * Wise2_cDNA_alloc
 * Wise2_replace_baseseq_cDNA
 * Wise2_access_baseseq_cDNA
 * Wise2_free_cDNA [destructor]
 *
 */

/* API for object cDNA */
/* Function:  Wise2_truncate_cDNA(cdna,start,stop)
 *
 * Descrip:    Truncates a cDNA sequence. Basically uses
 *             the /magic_trunc_Sequence function (of course!)
 *
 *             It does not alter cdna, rather it returns a new
 *             sequence with that truncation
 *
 *
 * Arg:        cdna         cDNA that is truncated [Wise2_cDNA *]
 * Arg:        start        Undocumented argument [int]
 * Arg:        stop         Undocumented argument [int]
 *
 * Returns Undocumented return value [Wise2_cDNA *]
 *
 */
Wise2_cDNA * Wise2_truncate_cDNA( Wise2_cDNA * cdna,int start,int stop);

/* Function:  Wise2_read_fasta_file_cDNA(filename)
 *
 * Descrip:    Reads a fasta file assumming that it is cDNA. 
 *             Will complain if it is not, and return NULL.
 *
 *
 * Arg:        filename     filename to be opened and read [char *]
 *
 * Returns Undocumented return value [Wise2_cDNA *]
 *
 */
Wise2_cDNA * Wise2_read_fasta_file_cDNA( char * filename);

/* Function:  Wise2_cDNA_name(cdna)
 *
 * Descrip:    Returns the name of the cDNA
 *
 *
 * Arg:        cdna         Undocumented argument [Wise2_cDNA *]
 *
 * Returns Undocumented return value [char *]
 *
 */
char * Wise2_cDNA_name( Wise2_cDNA * cdna);

/* Function:  Wise2_cDNA_length(cdna)
 *
 * Descrip:    Returns the length of the cDNA
 *
 *
 * Arg:        cdna         Undocumented argument [Wise2_cDNA *]
 *
 * Returns Undocumented return value [int]
 *
 */
int Wise2_cDNA_length( Wise2_cDNA * cdna);

/* Function:  Wise2_cDNA_seqchar(cdna,pos)
 *
 * Descrip:    Returns sequence character at this position.
 *
 *
 * Arg:        cdna         cDNA [Wise2_cDNA *]
 * Arg:        pos          position in cDNA to get char [int]
 *
 * Returns Undocumented return value [char]
 *
 */
char Wise2_cDNA_seqchar( Wise2_cDNA * cdna,int pos);

/* Function:  Wise2_cDNA_from_Sequence(seq)
 *
 * Descrip:    makes a new cDNA from a Sequence. It 
 *             owns the Sequence memory, ie will attempt a /free_Sequence
 *             on the structure when /free_cDNA is called
 *
 *             If you want to give this cDNA this Sequence and
 *             forget about it, then just hand it this sequence and set
 *             seq to NULL (no need to free it). If you intend to use 
 *             the sequence object elsewhere outside of the cDNA datastructure
 *             then use cDNA_from_Sequence(/hard_link_Sequence(seq))
 *
 *
 *
 * Arg:        seq          Sequence to make cDNA from [Wise2_Sequence *]
 *
 * Returns Undocumented return value [Wise2_cDNA *]
 *
 */
Wise2_cDNA * Wise2_cDNA_from_Sequence( Wise2_Sequence * seq);

/* Function:  Wise2_hard_link_cDNA(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_cDNA *]
 *
 * Returns Undocumented return value [Wise2_cDNA *]
 *
 */
Wise2_cDNA * Wise2_hard_link_cDNA( Wise2_cDNA * obj);

/* Function:  Wise2_cDNA_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Returns Undocumented return value [Wise2_cDNA *]
 *
 */
Wise2_cDNA * Wise2_cDNA_alloc();

/* Function:  Wise2_replace_baseseq_cDNA(obj,baseseq)
 *
 * Descrip:    Replace member variable baseseq
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_cDNA *]
 * Arg:        baseseq      New value of the variable [Wise2_Sequence *]
 *
 * Returns member variable baseseq [boolean]
 *
 */
boolean Wise2_replace_baseseq_cDNA( Wise2_cDNA * obj,Wise2_Sequence * baseseq);

/* Function:  Wise2_access_baseseq_cDNA(obj)
 *
 * Descrip:    Access member variable baseseq
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_cDNA *]
 *
 * Returns member variable baseseq [Wise2_Sequence *]
 *
 */
Wise2_Sequence * Wise2_access_baseseq_cDNA( Wise2_cDNA * obj);

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_cDNA(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_cDNA *]
 *
 * Returns Undocumented return value [Wise2_cDNA *]
 *
 */
Wise2_cDNA * Wise2_free_cDNA( Wise2_cDNA * obj);

