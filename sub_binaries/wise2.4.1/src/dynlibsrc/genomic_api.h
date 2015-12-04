

/* Functions that create, manipulate or act on Genomic
 *
 * Wise2_truncate_Genomic
 * Wise2_Genomic_name
 * Wise2_Genomic_length
 * Wise2_Genomic_seqchar
 * Wise2_show_Genomic
 * Wise2_hard_link_Genomic
 * Wise2_Genomic_alloc_std
 * Wise2_replace_baseseq_Genomic
 * Wise2_access_baseseq_Genomic
 * Wise2_access_repeat_Genomic
 * Wise2_length_repeat_Genomic
 * Wise2_flush_Genomic
 * Wise2_add_Genomic
 * Wise2_free_Genomic [destructor]
 *
 */



/* Functions that create, manipulate or act on GenomicRepeat
 *
 * Wise2_hard_link_GenomicRepeat
 * Wise2_GenomicRepeat_alloc
 * Wise2_replace_start_GenomicRepeat
 * Wise2_access_start_GenomicRepeat
 * Wise2_replace_end_GenomicRepeat
 * Wise2_access_end_GenomicRepeat
 * Wise2_replace_type_GenomicRepeat
 * Wise2_access_type_GenomicRepeat
 * Wise2_free_GenomicRepeat [destructor]
 *
 */



/* Helper functions in the module
 *
 * Wise2_read_fasta_file_Genomic
 * Wise2_read_fasta_Genomic
 * Wise2_Genomic_from_Sequence_Nheuristic
 * Wise2_Genomic_from_Sequence
 *

/* API for object Genomic */
/* Function:  Wise2_truncate_Genomic(gen,start,stop)
 *
 * Descrip:    Truncates a Genomic sequence. Basically uses
 *             the /magic_trunc_Sequence function (of course!)
 *
 *             It does not alter gen, rather it returns a new
 *             sequence with that truncation
 *
 *             Handles repeat information correctly.
 *
 *
 * Arg:        gen          Genomic that is truncated [Wise2_Genomic *]
 * Arg:        start        Undocumented argument [int]
 * Arg:        stop         Undocumented argument [int]
 *
 * Returns Undocumented return value [Wise2_Genomic *]
 *
 */
Wise2_Genomic * Wise2_truncate_Genomic( Wise2_Genomic * gen,int start,int stop);

/* Function:  Wise2_Genomic_name(gen)
 *
 * Descrip:    Returns the name of the Genomic
 *
 *
 * Arg:        gen          Undocumented argument [Wise2_Genomic *]
 *
 * Returns Undocumented return value [char *]
 *
 */
char * Wise2_Genomic_name( Wise2_Genomic * gen);

/* Function:  Wise2_Genomic_length(gen)
 *
 * Descrip:    Returns the length of the Genomic
 *
 *
 * Arg:        gen          Undocumented argument [Wise2_Genomic *]
 *
 * Returns Undocumented return value [int]
 *
 */
int Wise2_Genomic_length( Wise2_Genomic * gen);

/* Function:  Wise2_Genomic_seqchar(gen,pos)
 *
 * Descrip:    Returns sequence character at this position.
 *
 *
 * Arg:        gen          Genomic [Wise2_Genomic *]
 * Arg:        pos          position in Genomic to get char [int]
 *
 * Returns Undocumented return value [char]
 *
 */
char Wise2_Genomic_seqchar( Wise2_Genomic * gen,int pos);

/* Function:  Wise2_show_Genomic(gen,ofp)
 *
 * Descrip:    For debugging
 *
 *
 * Arg:        gen          Undocumented argument [Wise2_Genomic *]
 * Arg:        ofp          Undocumented argument [FILE *]
 *
 * Returns Undocumented return value [void]
 *
 */
void Wise2_show_Genomic( Wise2_Genomic * gen,FILE * ofp);

/* Function:  Wise2_hard_link_Genomic(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_Genomic *]
 *
 * Returns Undocumented return value [Wise2_Genomic *]
 *
 */
Wise2_Genomic * Wise2_hard_link_Genomic( Wise2_Genomic * obj);

/* Function:  Wise2_Genomic_alloc_std(void)
 *
 * Descrip:    Equivalent to Genomic_alloc_len(GenomicLISTLENGTH)
 *
 *
 *
 * Returns Undocumented return value [Wise2_Genomic *]
 *
 */
Wise2_Genomic * Wise2_Genomic_alloc_std();

/* Function:  Wise2_replace_baseseq_Genomic(obj,baseseq)
 *
 * Descrip:    Replace member variable baseseq
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Genomic *]
 * Arg:        baseseq      New value of the variable [Wise2_Sequence *]
 *
 * Returns member variable baseseq [boolean]
 *
 */
boolean Wise2_replace_baseseq_Genomic( Wise2_Genomic * obj,Wise2_Sequence * baseseq);

/* Function:  Wise2_access_baseseq_Genomic(obj)
 *
 * Descrip:    Access member variable baseseq
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Genomic *]
 *
 * Returns member variable baseseq [Wise2_Sequence *]
 *
 */
Wise2_Sequence * Wise2_access_baseseq_Genomic( Wise2_Genomic * obj);

/* Function:  Wise2_access_repeat_Genomic(obj,i)
 *
 * Descrip:    Access members stored in the repeat list
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the list [Wise2_Genomic *]
 * Arg:        i            Position in the list [int]
 *
 * Returns Element of the list [Wise2_GenomicRepeat *]
 *
 */
Wise2_GenomicRepeat * Wise2_access_repeat_Genomic( Wise2_Genomic * obj,int i);

/* Function:  Wise2_length_repeat_Genomic(obj)
 *
 * Descrip:    discover the length of the list
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the list [Wise2_Genomic *]
 *
 * Returns length of the list [int]
 *
 */
int Wise2_length_repeat_Genomic( Wise2_Genomic * obj);

/* Function:  Wise2_flush_Genomic(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj          Object which contains the list  [Wise2_Genomic *]
 *
 * Returns Undocumented return value [int]
 *
 */
int Wise2_flush_Genomic( Wise2_Genomic * obj);

/* Function:  Wise2_add_Genomic(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj          Object which contains the list [Wise2_Genomic *]
 * Arg:        add          Object to add to the list [Wise2_GenomicRepeat *]
 *
 * Returns Undocumented return value [boolean]
 *
 */
boolean Wise2_add_Genomic( Wise2_Genomic * obj,Wise2_GenomicRepeat * add);

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_Genomic(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_Genomic *]
 *
 * Returns Undocumented return value [Wise2_Genomic *]
 *
 */
Wise2_Genomic * Wise2_free_Genomic( Wise2_Genomic * obj);

/* API for object GenomicRepeat */
/* Function:  Wise2_hard_link_GenomicRepeat(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_GenomicRepeat *]
 *
 * Returns Undocumented return value [Wise2_GenomicRepeat *]
 *
 */
Wise2_GenomicRepeat * Wise2_hard_link_GenomicRepeat( Wise2_GenomicRepeat * obj);

/* Function:  Wise2_GenomicRepeat_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Returns Undocumented return value [Wise2_GenomicRepeat *]
 *
 */
Wise2_GenomicRepeat * Wise2_GenomicRepeat_alloc();

/* Function:  Wise2_replace_start_GenomicRepeat(obj,start)
 *
 * Descrip:    Replace member variable start
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_GenomicRepeat *]
 * Arg:        start        New value of the variable [int]
 *
 * Returns member variable start [boolean]
 *
 */
boolean Wise2_replace_start_GenomicRepeat( Wise2_GenomicRepeat * obj,int start);

/* Function:  Wise2_access_start_GenomicRepeat(obj)
 *
 * Descrip:    Access member variable start
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_GenomicRepeat *]
 *
 * Returns member variable start [int]
 *
 */
int Wise2_access_start_GenomicRepeat( Wise2_GenomicRepeat * obj);

/* Function:  Wise2_replace_end_GenomicRepeat(obj,end)
 *
 * Descrip:    Replace member variable end
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_GenomicRepeat *]
 * Arg:        end          New value of the variable [int]
 *
 * Returns member variable end [boolean]
 *
 */
boolean Wise2_replace_end_GenomicRepeat( Wise2_GenomicRepeat * obj,int end);

/* Function:  Wise2_access_end_GenomicRepeat(obj)
 *
 * Descrip:    Access member variable end
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_GenomicRepeat *]
 *
 * Returns member variable end [int]
 *
 */
int Wise2_access_end_GenomicRepeat( Wise2_GenomicRepeat * obj);

/* Function:  Wise2_replace_type_GenomicRepeat(obj,type)
 *
 * Descrip:    Replace member variable type
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_GenomicRepeat *]
 * Arg:        type         New value of the variable [char *]
 *
 * Returns member variable type [boolean]
 *
 */
boolean Wise2_replace_type_GenomicRepeat( Wise2_GenomicRepeat * obj,char * type);

/* Function:  Wise2_access_type_GenomicRepeat(obj)
 *
 * Descrip:    Access member variable type
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_GenomicRepeat *]
 *
 * Returns member variable type [char *]
 *
 */
char * Wise2_access_type_GenomicRepeat( Wise2_GenomicRepeat * obj);

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_GenomicRepeat(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_GenomicRepeat *]
 *
 * Returns Undocumented return value [Wise2_GenomicRepeat *]
 *
 */
Wise2_GenomicRepeat * Wise2_free_GenomicRepeat( Wise2_GenomicRepeat * obj);



/* These functions are not associated with an object */
/* Function:  Wise2_read_fasta_file_Genomic(filename,length_of_N)
 *
 * Descrip:    Reads a fasta file assumming that it is Genomic. 
 *             Will complain if it is not, and return NULL.
 *
 *
 * Arg:        filename     filename to be opened and read [char *]
 * Arg:        length_of_N  length of N to be considered repeat. -1 means none [int]
 *
 * Returns Undocumented return value [Wise2_Genomic *]
 *
 */
Wise2_Genomic * Wise2_read_fasta_file_Genomic( char * filename,int length_of_N);

/* Function:  Wise2_read_fasta_Genomic(ifp,length_of_N)
 *
 * Descrip:    Reads a fasta file assumming that it is Genomic. 
 *             Will complain if it is not, and return NULL.
 *
 *
 * Arg:        ifp          file point to be read from [FILE *]
 * Arg:        length_of_N  length of N to be considered repeat. -1 means none [int]
 *
 * Returns Undocumented return value [Wise2_Genomic *]
 *
 */
Wise2_Genomic * Wise2_read_fasta_Genomic( FILE * ifp,int length_of_N);

/* Function:  Wise2_Genomic_from_Sequence_Nheuristic(seq,length_of_N)
 *
 * Descrip:    makes a new genomic from a Sequence, but
 *             assummes that all the N runs greater than
 *             a certain level are actually repeats.
 *
 *
 * Arg:        seq          Undocumented argument [Wise2_Sequence *]
 * Arg:        length_of_N  Undocumented argument [int]
 *
 * Returns Undocumented return value [Wise2_Genomic *]
 *
 */
Wise2_Genomic * Wise2_Genomic_from_Sequence_Nheuristic( Wise2_Sequence * seq,int length_of_N);

/* Function:  Wise2_Genomic_from_Sequence(seq)
 *
 * Descrip:    makes a new genomic from a Sequence. It 
 *             owns the Sequence memory, ie will attempt a /free_Sequence
 *             on the structure when /free_Genomic is called
 *
 *             If you want to give this genomic this Sequence and
 *             forget about it, then just hand it this sequence and set
 *             seq to NULL (no need to free it). If you intend to use 
 *             the sequence object elsewhere outside of the Genomic datastructure
 *             then use Genomic_from_Sequence(/hard_link_Sequence(seq))
 *
 *             This is part of a strict typing system, and therefore
 *             is going to convert all non ATGCNs to Ns. You will lose
 *             information here.
 *
 *
 * Arg:        seq          Sequence to make genomic from [Wise2_Sequence *]
 *
 * Returns Undocumented return value [Wise2_Genomic *]
 *
 */
Wise2_Genomic * Wise2_Genomic_from_Sequence( Wise2_Sequence * seq);

