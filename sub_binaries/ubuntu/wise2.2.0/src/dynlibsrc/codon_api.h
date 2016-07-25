

/* Functions that create, manipulate or act on CodonTable
 *
 * Wise2_read_CodonTable_file
 * Wise2_read_CodonTable
 * Wise2_aminoacid_from_seq
 * Wise2_aminoacid_from_codon
 * Wise2_is_stop_codon
 * Wise2_is_valid_aminoacid
 * Wise2_hard_link_CodonTable
 * Wise2_CodonTable_alloc
 * Wise2_replace_name_CodonTable
 * Wise2_access_name_CodonTable
 * Wise2_free_CodonTable [destructor]
 *
 */



/* Helper functions in the module
 *
 * Wise2_is_non_ambiguous_codon_seq
 * Wise2_codon_from_base4_codon
 * Wise2_base4_codon_from_codon
 * Wise2_has_random_bases
 * Wise2_permute_possible_random_bases
 * Wise2_base_from_codon
 * Wise2_codon_from_seq
 * Wise2_base4_codon_from_seq
 * Wise2_char_from_base
 * Wise2_base_from_char
 * Wise2_char_complement_base
 * Wise2_complement_base
 *

/* API for object CodonTable */
/* Function:  Wise2_read_CodonTable_file(file)
 *
 * Descrip:    Opens filename, reads it as if a Ewan style
 *             codon table and closes.
 *
 *
 * Arg:        file         filename to open [char *]
 *
 * Returns A codon-table, NULL if error [Wise2_CodonTable *]
 *
 */
Wise2_CodonTable * Wise2_read_CodonTable_file( char * file);

/* Function:  Wise2_read_CodonTable(ifp)
 *
 * Descrip:    reads a codon table from a filestream in Ewan
 *             format.
 *
 *             As Ewan format is really bad and has no start/stop
 *             this will effectively read to the end of the file.
 *             Ooops.
 *
 *
 * Arg:        ifp          file input [FILE *]
 *
 * Returns Undocumented return value [Wise2_CodonTable *]
 *
 */
Wise2_CodonTable * Wise2_read_CodonTable( FILE * ifp);

/* Function:  Wise2_aminoacid_from_seq(ct,seq)
 *
 * Descrip:    Returns the amino acid for this position in the DNA sequence
 *             Takes the pointer +1 and +2 points.
 *
 *             No error checks implemented. Probably a mistake ;)
 *
 *
 * Arg:        ct           codon table [Wise2_CodonTable *]
 * Arg:        seq          pointer to DNA chars [char *]
 *
 * Returns an amino acid char (A-Z) [aa]
 *
 */
aa Wise2_aminoacid_from_seq( Wise2_CodonTable * ct,char * seq);

/* Function:  Wise2_aminoacid_from_codon(ct,c)
 *
 * Descrip:    returns amino acid for this codon number (NB codon numbers 0-125)
 *
 *
 * Arg:        ct           codon table [Wise2_CodonTable *]
 * Arg:        c            codon number [codon]
 *
 * Returns aminoacid that is this codon (X for ambiguous, * for stop) [aa]
 *
 */
aa Wise2_aminoacid_from_codon( Wise2_CodonTable * ct,codon c);

/* Function:  Wise2_is_stop_codon(c,ct)
 *
 * Descrip:    tells you whether this codon number is really a stop
 *             in this translation table
 *
 *
 * Arg:        c            codon number [codon]
 * Arg:        ct           codon table [Wise2_CodonTable *]
 *
 * Returns TRUE if is stop, FALSE otherwise [boolean]
 *
 */
boolean Wise2_is_stop_codon( codon c,Wise2_CodonTable * ct);

/* Function:  Wise2_is_valid_aminoacid(ct,c)
 *
 * Descrip:    Tells you if this letter (c) is recognised as a valid amino acid
 *             in this codon table
 *
 *
 * Arg:        ct           Codon Table [Wise2_CodonTable *]
 * Arg:        c            aminoacid [char]
 *
 * Returns TRUE if valid, FALSE if not. [boolean]
 *
 */
boolean Wise2_is_valid_aminoacid( Wise2_CodonTable * ct,char c);

/* Function:  Wise2_hard_link_CodonTable(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_CodonTable *]
 *
 * Returns Undocumented return value [Wise2_CodonTable *]
 *
 */
Wise2_CodonTable * Wise2_hard_link_CodonTable( Wise2_CodonTable * obj);

/* Function:  Wise2_CodonTable_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Returns Undocumented return value [Wise2_CodonTable *]
 *
 */
Wise2_CodonTable * Wise2_CodonTable_alloc();

/* Function:  Wise2_replace_name_CodonTable(obj,name)
 *
 * Descrip:    Replace member variable name
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_CodonTable *]
 * Arg:        name         New value of the variable [char *]
 *
 * Returns member variable name [boolean]
 *
 */
boolean Wise2_replace_name_CodonTable( Wise2_CodonTable * obj,char * name);

/* Function:  Wise2_access_name_CodonTable(obj)
 *
 * Descrip:    Access member variable name
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_CodonTable *]
 *
 * Returns member variable name [char *]
 *
 */
char * Wise2_access_name_CodonTable( Wise2_CodonTable * obj);

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_CodonTable(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_CodonTable *]
 *
 * Returns Undocumented return value [Wise2_CodonTable *]
 *
 */
Wise2_CodonTable * Wise2_free_CodonTable( Wise2_CodonTable * obj);



/* These functions are not associated with an object */
/* Function:  Wise2_is_non_ambiguous_codon_seq(seq)
 *
 * Descrip:    Tells you if this codon is a real codon
 *
 *
 * Arg:        seq          pointer to DNA sequence [char *]
 *
 * Returns TRUE if real codon, FALSE if contains N's [boolean]
 *
 */
boolean Wise2_is_non_ambiguous_codon_seq( char * seq);

/* Function:  Wise2_codon_from_base4_codon(c)
 *
 * Descrip:    maps a 0-63 codon to a 0-123 codon. Suprisingly useful.
 *
 *
 * Arg:        c            Undocumented argument [int]
 *
 * Returns Undocumented return value [codon]
 *
 */
codon Wise2_codon_from_base4_codon( int c);

/* Function:  Wise2_base4_codon_from_codon(c)
 *
 * Descrip:    maps a 0-125 codon to a 0-63 codon.
 *
 *             If ambiguous then returns 64 having issued a warning.
 *
 *
 * Arg:        c            codon 0-125 [codon]
 *
 * Returns base 4 codon (0-63) [int]
 *
 */
int Wise2_base4_codon_from_codon( codon c);

/* Function:  Wise2_has_random_bases(c)
 *
 * Descrip:    Tests to see if this codon number has any N's in it
 *
 *
 * Arg:        c            codon number 0-124 [codon]
 *
 * Returns TRUE if has N's , FALSE otherwise [boolean]
 *
 */
boolean Wise2_has_random_bases( codon c);

/* Function:  Wise2_permute_possible_random_bases(c,one,two,three)
 *
 * Descrip:    Bizarely useful function for calculating ambiguity scores.
 *
 *             This takes the codon c, and for each possible base, 
 *             if it is N, replaces it with one, two or three.
 *
 *             If the base is not N, it remains the same
 *
 *
 * Arg:        c            codon number [codon]
 * Arg:        one          base to replace first position if N [base]
 * Arg:        two          base to replace second position if N [base]
 * Arg:        three        base to replace third position if N [base]
 *
 * Returns codon number  [codon]
 *
 */
codon Wise2_permute_possible_random_bases( codon c,base one,base two,base three);

/* Function:  Wise2_base_from_codon(c,pos)
 *
 * Descrip:    Probably not the best function to use for this, but 
 *             useful. Takes a codon and with pos being 1,2,3 gives
 *             you the firt,second of third base
 *
 *
 * Arg:        c            Undocumented argument [codon]
 * Arg:        pos          Undocumented argument [int]
 *
 * Returns Undocumented return value [base]
 *
 */
base Wise2_base_from_codon( codon c,int pos);

/* Function:  Wise2_codon_from_seq(seq)
 *
 * Descrip:    takes an ASCII coded pointer to a 3 base pair
 *             sequence (it could be the part of a sequence: it only
 *             assummes that the seq points with 3 chars at pos 0,1,2 
 *             in C coordinates from seq. No NULL is required). It 
 *             ives back the codon as made from standard mapping, ie,
 *             25*base_1+5*base_2 + base3 being a number from 0-124 inc.
 *
 *
 * Arg:        seq          pointer to sequence of at least 3 chrs long. [char *]
 *
 * Returns Undocumented return value [codon]
 *
 */
codon Wise2_codon_from_seq( char * seq);

/* Function:  Wise2_base4_codon_from_seq(seq)
 *
 * Descrip:    Sometimes it is more useful to work in base64, ie, 
 *             non N. this functions does the same thing as 
 *             /codon_from_seq but produces a seq being
 *             16*base1 + 4 *base2 + base3
 *
 *
 * Arg:        seq          pointer to sequence of at least 3 chrs long [char *]
 *
 * Returns Undocumented return value [int]
 *
 */
int Wise2_base4_codon_from_seq( char * seq);

/* Function:  Wise2_char_from_base(b)
 *
 * Descrip:    maps a base number (-04 inc) to A,T,G,C,N
 *
 *
 * Arg:        b            Undocumented argument [base]
 *
 * Returns Undocumented return value [char]
 *
 */
char Wise2_char_from_base( base b);

/* Function:  Wise2_base_from_char(c)
 *
 * Descrip:    maps a char (atcgn) to number, 
 *             case insensitive, returns BASE_N
 *             if not atcgn
 *
 *
 * Arg:        c            Undocumented argument [char]
 *
 * Returns Undocumented return value [base]
 *
 */
base Wise2_base_from_char( char c);

/* Function:  Wise2_char_complement_base(c)
 *
 * Descrip:    the char equivalent of /complement_base.
 *             this gives the complement in char of a base
 *             in char. Does not check for non ATGCN
 *
 *
 * Arg:        c            Undocumented argument [char]
 *
 * Returns Undocumented return value [char]
 *
 */
char Wise2_char_complement_base( char c);

/* Function:  Wise2_complement_base(b)
 *
 * Descrip:    gives back the complement as a number
 *             ofthe base (given as a number)
 *
 *
 * Arg:        b            Undocumented argument [base]
 *
 * Returns Undocumented return value [base]
 *
 */
base Wise2_complement_base( base b);

