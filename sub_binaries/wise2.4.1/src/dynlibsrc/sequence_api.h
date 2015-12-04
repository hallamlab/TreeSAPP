

/* Functions that create, manipulate or act on Sequence
 *
 * Wise2_uppercase_Sequence
 * Wise2_force_to_dna_Sequence
 * Wise2_is_reversed_Sequence
 * Wise2_translate_Sequence
 * Wise2_reverse_complement_Sequence
 * Wise2_magic_trunc_Sequence
 * Wise2_trunc_Sequence
 * Wise2_read_fasta_file_Sequence
 * Wise2_read_Sequence_EMBL_seq
 * Wise2_read_fasta_Sequence
 * Wise2_show_Sequence_residue_list
 * Wise2_write_fasta_Sequence
 * Wise2_make_len_type_Sequence
 * Wise2_hard_link_Sequence
 * Wise2_Sequence_alloc
 * Wise2_replace_name_Sequence
 * Wise2_access_name_Sequence
 * Wise2_replace_seq_Sequence
 * Wise2_access_seq_Sequence
 * Wise2_replace_len_Sequence
 * Wise2_access_len_Sequence
 * Wise2_replace_maxlen_Sequence
 * Wise2_access_maxlen_Sequence
 * Wise2_replace_offset_Sequence
 * Wise2_access_offset_Sequence
 * Wise2_replace_end_Sequence
 * Wise2_access_end_Sequence
 * Wise2_replace_type_Sequence
 * Wise2_access_type_Sequence
 * Wise2_replace_tax_id_Sequence
 * Wise2_access_tax_id_Sequence
 * Wise2_replace_weight_Sequence
 * Wise2_access_weight_Sequence
 * Wise2_replace_desc_Sequence
 * Wise2_access_desc_Sequence
 * Wise2_free_Sequence [destructor]
 *
 */



/* Functions that create, manipulate or act on SequenceSet
 *
 * Wise2_hard_link_SequenceSet
 * Wise2_SequenceSet_alloc_std
 * Wise2_access_set_SequenceSet
 * Wise2_length_set_SequenceSet
 * Wise2_flush_SequenceSet
 * Wise2_add_SequenceSet
 * Wise2_free_SequenceSet [destructor]
 *
 */



/* Helper functions in the module
 *
 * Wise2_Sequence_type_to_string
 * Wise2_new_Sequence_from_strings
 *

/* API for object Sequence */
/* Function:  Wise2_uppercase_Sequence(seq)
 *
 * Descrip:    makes all the sequence uppercase
 *
 *
 * Arg:        seq          Sequence to be uppercased [Wise2_Sequence *]
 *
 * Returns Undocumented return value [void]
 *
 */
void Wise2_uppercase_Sequence( Wise2_Sequence * seq);

/* Function:  Wise2_force_to_dna_Sequence(seq,fraction,number_of_conver)
 *
 * Descrip:    This 
 *              a) sees how many non ATGCN characters there are in Seq
 *              b) If the level is below fraction
 *                 a) flips non ATGC chars to N
 *                 b) writes number of conversions to number_of_conver
 *                 c) returns TRUE
 *              c) else returns FALSE
 *
 *             fraction of 0.0 means completely intolerant of errors
 *             fraction of 1.0 means completely tolerant of errors
 *
 *
 *
 * Arg:        seq          sequence object read and converted  [Wise2_Sequence *]
 * Arg:        fraction     number 0..1 for tolerance of conversion [double]
 * Arg:        number_of_conver number of conversions actually made [int *]
 *
 * Returns TRUE for conversion to DNA, FALSE if not [boolean]
 *
 */
boolean Wise2_force_to_dna_Sequence( Wise2_Sequence * seq,double fraction,int * number_of_conver);

/* Function:  Wise2_is_reversed_Sequence(seq)
 *
 * Descrip:    Currently the sequence object stores 
 *             reversed sequences as start > end.
 *
 *             This tests that and returns true if it is
 *
 *
 * Arg:        seq          sequence to test [Wise2_Sequence *]
 *
 * Returns Undocumented return value [boolean]
 *
 */
boolean Wise2_is_reversed_Sequence( Wise2_Sequence * seq);

/* Function:  Wise2_translate_Sequence(dna,ct)
 *
 * Descrip:    This translates a DNA sequence to a protein.
 *             It assummes that it starts at first residue
 *             (use trunc_Sequence to chop a sequence up).
 *
 *
 * Arg:        dna          DNA sequence to be translated [Wise2_Sequence *]
 * Arg:        ct           Codon table to do codon->aa mapping [Wise2_CodonTable *]
 *
 * Returns new protein sequence [Wise2_Sequence *]
 *
 */
Wise2_Sequence * Wise2_translate_Sequence( Wise2_Sequence * dna,Wise2_CodonTable * ct);

/* Function:  Wise2_reverse_complement_Sequence(seq)
 *
 * Descrip:    This both complements and reverses a sequence,
 *             - a common wish!
 *
 *             The start/end are correct with respect to the start/end
 *             of the sequence (ie start = end, end = start).
 *
 *
 * Arg:        seq          Sequence to that is used to reverse (makes a new Sequence) [Wise2_Sequence *]
 *
 * Returns new Sequence which is reversed [Wise2_Sequence *]
 *
 */
Wise2_Sequence * Wise2_reverse_complement_Sequence( Wise2_Sequence * seq);

/* Function:  Wise2_magic_trunc_Sequence(seq,start,end)
 *
 * Descrip:    Clever function for dna sequences.
 *
 *             When start < end, truncates normally
 *
 *             when start > end, truncates end,start and then
 *             reverse complements.
 *
 *             ie. If you have a coordinate system where reverse 
 *             sequences are labelled in reverse start/end way,
 *             then this routine produces the correct sequence.
 *
 *
 * Arg:        seq          sequence that is the source to be truncated [Wise2_Sequence *]
 * Arg:        start        start point [int]
 * Arg:        end          end point [int]
 *
 * Returns new Sequence which is truncated/reversed [Wise2_Sequence *]
 *
 */
Wise2_Sequence * Wise2_magic_trunc_Sequence( Wise2_Sequence * seq,int start,int end);

/* Function:  Wise2_trunc_Sequence(seq,start,end)
 *
 * Descrip:    truncates a sequence. It produces a new memory structure
 *             which is filled from sequence start to end.
 *
 *             Please notice
 *               
 *               Truncation is in C coordinates. That is
 *             the first residue is 0 and end is the number of the
 *             residue after the cut-point. In otherwords to 
 *             2 - 3 would be a single residue truncation. So - if
 *             you want to work in more usual, 'inclusive' molecular
 *             biology numbers, which start at 1, then you need to say
 *
 *               trunc_Sequence(seq,start-1,end);
 *
 *             (NB, should be (end - 1 + 1) = end for the last coordinate).
 *
 *               Truncation occurs against the *absolute* coordinate
 *             system of the Sequence, not the offset/end pair inside.
 *             So, this is a very bad error
 *              
 *               ** wrong code, and also leaks memory **
 *
 *               tru = trunc_Sequence(trunc_Sequence(seq,50,80),55,75); 
 *
 *             This the most portable way of doing this
 *
 *               temp = trunc_Sequence(seq,50,80);
 *
 *               tru  = trunc_Sequence(temp,55-temp->offset,75-temp->offset);
 *
 *               free_Sequence(temp);
 *
 *
 *
 * Arg:        seq          object holding the sequence to be truncated [Wise2_Sequence *]
 * Arg:        start        start point of truncation [int]
 * Arg:        end          end point of truncation [int]
 *
 * Returns newly allocated sequence structure [Wise2_Sequence *]
 *
 */
Wise2_Sequence * Wise2_trunc_Sequence( Wise2_Sequence * seq,int start,int end);

/* Function:  Wise2_read_fasta_file_Sequence(filename)
 *
 * Descrip:    Just a call
 *               a) open filename
 *               b) read sequence with /read_fasta_Sequence
 *               c) close file.
 *
 *
 * Arg:        filename     filename to open  [char *]
 *
 * Returns Undocumented return value [Wise2_Sequence *]
 *
 */
Wise2_Sequence * Wise2_read_fasta_file_Sequence( char * filename);

/* Function:  Wise2_read_Sequence_EMBL_seq(buffer,maxlen,ifp)
 *
 * Descrip:    reads the sequence part of an EMBL file.
 *
 *             This function can either take a file which 
 *             starts
 *
 *
 *
 * Arg:        buffer       buffer containing the first line. [char *]
 * Arg:        maxlen       length of buffer [int]
 * Arg:        ifp          input file to read from [FILE *]
 *
 * Returns Undocumented return value [Wise2_Sequence *]
 *
 */
Wise2_Sequence * Wise2_read_Sequence_EMBL_seq( char * buffer,int maxlen,FILE * ifp);

/* Function:  Wise2_read_fasta_Sequence(ifp)
 *
 * Descrip:    reads a fasta file assumming pretty 
 *             standard sanity for the file layout
 *
 *
 * Arg:        ifp          Undocumented argument [FILE *]
 *
 * Returns Undocumented return value [Wise2_Sequence *]
 *
 */
Wise2_Sequence * Wise2_read_fasta_Sequence( FILE * ifp);

/* Function:  Wise2_show_Sequence_residue_list(seq,start,end,ofp)
 *
 * Descrip:    shows a region of a sequence as
 *                124  A
 *                125  T
 *
 *             etc from start to end. The numbers
 *             are in C coordinates (ie, 0 is the first
 *             letter).
 *
 *             useful for debugging
 *
 *
 * Arg:        seq          Sequence to show [Wise2_Sequence *]
 * Arg:        start        start of list [int]
 * Arg:        end          end of list [int]
 * Arg:        ofp          Undocumented argument [FILE *]
 *
 * Returns Undocumented return value [void]
 *
 */
void Wise2_show_Sequence_residue_list( Wise2_Sequence * seq,int start,int end,FILE * ofp);

/* Function:  Wise2_write_fasta_Sequence(seq,ofp)
 *
 * Descrip:    writes a fasta file of the form
 *             >name
 *             Sequence
 *
 *
 * Arg:        seq          sequence to be written [Wise2_Sequence *]
 * Arg:        ofp          file to write to [FILE *]
 *
 * Returns Undocumented return value [void]
 *
 */
void Wise2_write_fasta_Sequence( Wise2_Sequence * seq,FILE * ofp);

/* Function:  Wise2_make_len_type_Sequence(seq)
 *
 * Descrip:    makes seq->len and seq->end match the seq->seq
 *             length number. 
 *
 *             It also checks the type of the sequence with
 *             /best_guess_type
 *
 *
 * Arg:        seq          Sequence object [Wise2_Sequence *]
 *
 * Returns Undocumented return value [void]
 *
 */
void Wise2_make_len_type_Sequence( Wise2_Sequence * seq);

/* Function:  Wise2_hard_link_Sequence(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_Sequence *]
 *
 * Returns Undocumented return value [Wise2_Sequence *]
 *
 */
Wise2_Sequence * Wise2_hard_link_Sequence( Wise2_Sequence * obj);

/* Function:  Wise2_Sequence_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Returns Undocumented return value [Wise2_Sequence *]
 *
 */
Wise2_Sequence * Wise2_Sequence_alloc();

/* Function:  Wise2_replace_name_Sequence(obj,name)
 *
 * Descrip:    Replace member variable name
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Sequence *]
 * Arg:        name         New value of the variable [char *]
 *
 * Returns member variable name [boolean]
 *
 */
boolean Wise2_replace_name_Sequence( Wise2_Sequence * obj,char * name);

/* Function:  Wise2_access_name_Sequence(obj)
 *
 * Descrip:    Access member variable name
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Sequence *]
 *
 * Returns member variable name [char *]
 *
 */
char * Wise2_access_name_Sequence( Wise2_Sequence * obj);

/* Function:  Wise2_replace_seq_Sequence(obj,seq)
 *
 * Descrip:    Replace member variable seq
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Sequence *]
 * Arg:        seq          New value of the variable [char *]
 *
 * Returns member variable seq [boolean]
 *
 */
boolean Wise2_replace_seq_Sequence( Wise2_Sequence * obj,char * seq);

/* Function:  Wise2_access_seq_Sequence(obj)
 *
 * Descrip:    Access member variable seq
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Sequence *]
 *
 * Returns member variable seq [char *]
 *
 */
char * Wise2_access_seq_Sequence( Wise2_Sequence * obj);

/* Function:  Wise2_replace_len_Sequence(obj,len)
 *
 * Descrip:    Replace member variable len
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Sequence *]
 * Arg:        len          New value of the variable [int]
 *
 * Returns member variable len [boolean]
 *
 */
boolean Wise2_replace_len_Sequence( Wise2_Sequence * obj,int len);

/* Function:  Wise2_access_len_Sequence(obj)
 *
 * Descrip:    Access member variable len
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Sequence *]
 *
 * Returns member variable len [int]
 *
 */
int Wise2_access_len_Sequence( Wise2_Sequence * obj);

/* Function:  Wise2_replace_maxlen_Sequence(obj,maxlen)
 *
 * Descrip:    Replace member variable maxlen
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Sequence *]
 * Arg:        maxlen       New value of the variable [int]
 *
 * Returns member variable maxlen [boolean]
 *
 */
boolean Wise2_replace_maxlen_Sequence( Wise2_Sequence * obj,int maxlen);

/* Function:  Wise2_access_maxlen_Sequence(obj)
 *
 * Descrip:    Access member variable maxlen
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Sequence *]
 *
 * Returns member variable maxlen [int]
 *
 */
int Wise2_access_maxlen_Sequence( Wise2_Sequence * obj);

/* Function:  Wise2_replace_offset_Sequence(obj,offset)
 *
 * Descrip:    Replace member variable offset
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Sequence *]
 * Arg:        offset       New value of the variable [int]
 *
 * Returns member variable offset [boolean]
 *
 */
boolean Wise2_replace_offset_Sequence( Wise2_Sequence * obj,int offset);

/* Function:  Wise2_access_offset_Sequence(obj)
 *
 * Descrip:    Access member variable offset
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Sequence *]
 *
 * Returns member variable offset [int]
 *
 */
int Wise2_access_offset_Sequence( Wise2_Sequence * obj);

/* Function:  Wise2_replace_end_Sequence(obj,end)
 *
 * Descrip:    Replace member variable end
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Sequence *]
 * Arg:        end          New value of the variable [int]
 *
 * Returns member variable end [boolean]
 *
 */
boolean Wise2_replace_end_Sequence( Wise2_Sequence * obj,int end);

/* Function:  Wise2_access_end_Sequence(obj)
 *
 * Descrip:    Access member variable end
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Sequence *]
 *
 * Returns member variable end [int]
 *
 */
int Wise2_access_end_Sequence( Wise2_Sequence * obj);

/* Function:  Wise2_replace_type_Sequence(obj,type)
 *
 * Descrip:    Replace member variable type
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Sequence *]
 * Arg:        type         New value of the variable [char]
 *
 * Returns member variable type [boolean]
 *
 */
boolean Wise2_replace_type_Sequence( Wise2_Sequence * obj,char type);

/* Function:  Wise2_access_type_Sequence(obj)
 *
 * Descrip:    Access member variable type
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Sequence *]
 *
 * Returns member variable type [char]
 *
 */
char Wise2_access_type_Sequence( Wise2_Sequence * obj);

/* Function:  Wise2_replace_tax_id_Sequence(obj,tax_id)
 *
 * Descrip:    Replace member variable tax_id
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Sequence *]
 * Arg:        tax_id       New value of the variable [int]
 *
 * Returns member variable tax_id [boolean]
 *
 */
boolean Wise2_replace_tax_id_Sequence( Wise2_Sequence * obj,int tax_id);

/* Function:  Wise2_access_tax_id_Sequence(obj)
 *
 * Descrip:    Access member variable tax_id
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Sequence *]
 *
 * Returns member variable tax_id [int]
 *
 */
int Wise2_access_tax_id_Sequence( Wise2_Sequence * obj);

/* Function:  Wise2_replace_weight_Sequence(obj,weight)
 *
 * Descrip:    Replace member variable weight
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Sequence *]
 * Arg:        weight       New value of the variable [double]
 *
 * Returns member variable weight [boolean]
 *
 */
boolean Wise2_replace_weight_Sequence( Wise2_Sequence * obj,double weight);

/* Function:  Wise2_access_weight_Sequence(obj)
 *
 * Descrip:    Access member variable weight
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Sequence *]
 *
 * Returns member variable weight [double]
 *
 */
double Wise2_access_weight_Sequence( Wise2_Sequence * obj);

/* Function:  Wise2_replace_desc_Sequence(obj,desc)
 *
 * Descrip:    Replace member variable desc
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Sequence *]
 * Arg:        desc         New value of the variable [char *]
 *
 * Returns member variable desc [boolean]
 *
 */
boolean Wise2_replace_desc_Sequence( Wise2_Sequence * obj,char * desc);

/* Function:  Wise2_access_desc_Sequence(obj)
 *
 * Descrip:    Access member variable desc
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Sequence *]
 *
 * Returns member variable desc [char *]
 *
 */
char * Wise2_access_desc_Sequence( Wise2_Sequence * obj);

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_Sequence(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_Sequence *]
 *
 * Returns Undocumented return value [Wise2_Sequence *]
 *
 */
Wise2_Sequence * Wise2_free_Sequence( Wise2_Sequence * obj);

/* API for object SequenceSet */
/* Function:  Wise2_hard_link_SequenceSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_SequenceSet *]
 *
 * Returns Undocumented return value [Wise2_SequenceSet *]
 *
 */
Wise2_SequenceSet * Wise2_hard_link_SequenceSet( Wise2_SequenceSet * obj);

/* Function:  Wise2_SequenceSet_alloc_std(void)
 *
 * Descrip:    Equivalent to SequenceSet_alloc_len(SequenceSetLISTLENGTH)
 *
 *
 *
 * Returns Undocumented return value [Wise2_SequenceSet *]
 *
 */
Wise2_SequenceSet * Wise2_SequenceSet_alloc_std();

/* Function:  Wise2_access_set_SequenceSet(obj,i)
 *
 * Descrip:    Access members stored in the set list
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the list [Wise2_SequenceSet *]
 * Arg:        i            Position in the list [int]
 *
 * Returns Element of the list [Wise2_Sequence *]
 *
 */
Wise2_Sequence * Wise2_access_set_SequenceSet( Wise2_SequenceSet * obj,int i);

/* Function:  Wise2_length_set_SequenceSet(obj)
 *
 * Descrip:    discover the length of the list
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the list [Wise2_SequenceSet *]
 *
 * Returns length of the list [int]
 *
 */
int Wise2_length_set_SequenceSet( Wise2_SequenceSet * obj);

/* Function:  Wise2_flush_SequenceSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj          Object which contains the list  [Wise2_SequenceSet *]
 *
 * Returns Undocumented return value [int]
 *
 */
int Wise2_flush_SequenceSet( Wise2_SequenceSet * obj);

/* Function:  Wise2_add_SequenceSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj          Object which contains the list [Wise2_SequenceSet *]
 * Arg:        add          Object to add to the list [Wise2_Sequence *]
 *
 * Returns Undocumented return value [boolean]
 *
 */
boolean Wise2_add_SequenceSet( Wise2_SequenceSet * obj,Wise2_Sequence * add);

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_SequenceSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_SequenceSet *]
 *
 * Returns Undocumented return value [Wise2_SequenceSet *]
 *
 */
Wise2_SequenceSet * Wise2_free_SequenceSet( Wise2_SequenceSet * obj);



/* These functions are not associated with an object */
/* Function:  Wise2_Sequence_type_to_string(type)
 *
 * Descrip:    Converts sequence type (SEQUENCE_*) to a string
 *
 *
 * Arg:        type         type eg SEQUENCE_PROTEIN [int]
 *
 * Returns Undocumented return value [char *]
 *
 */
char * Wise2_Sequence_type_to_string( int type);

/* Function:  Wise2_new_Sequence_from_strings(name,seq)
 *
 * Descrip:    Makes a new sequence from strings given. 
 *             Separate memory will be allocated for them
 *             and them copied into it.
 *
 *             They can be NULL, in which case 
 *             o  a dummy name SequenceName will be assigned
 *             o  No sequence placed and length of zero.
 *
 *             Though this is dangerous later on. 
 *
 *             The sequence type is calculated automatically using
 *             /best_guess_type. If you want a DNA sequence but are
 *             unsure of the content of, for example, IUPAC codes,
 *             please use /force_to_dna_Sequence before using the
 *             sequence. Most of the rest of dynamite relies on a
 *             five letter A,T,G,C,N alphabet, but this function
 *             will allow any sequence type to be stored, so please
 *             check if you want to save yourself alot of grief.
 *
 *             In perl and other interfaces, this is a much safer
 *             constructor than the raw "new" type
 *
 *
 * Arg:        name         name of sequence, memory is allocated for it. [char *]
 * Arg:        seq          char * of sequence, memory is allocated for it. [char *]
 *
 * Returns Undocumented return value [Wise2_Sequence *]
 *
 */
Wise2_Sequence * Wise2_new_Sequence_from_strings( char * name,char * seq);

