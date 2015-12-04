#ifndef DYNAMITEcodonHEADERFILE
#define DYNAMITEcodonHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "wisebase.h"


#define BASE_A 0
#define BASE_G 1
#define BASE_C 2
#define BASE_T 3
#define BASE_N 4

typedef short int base;
typedef short int codon;
typedef char aa;

#define aminoacid_number(c) (c-'A')

#define is_base(a) (a == 'A' || a == 'T' || a == 'C' || a == 'G' || a== 'N' || a == 'a' || a == 't' || a == 'c' || a == 'g' || a == 'n' ? 1 : 0)

/* Object CodonTable
 *
 * Descrip: The codon table provides a mapping from the 64 codons to the 20 amino
 *        acids. The rest of the modules provides assorted codon<->base<->amino
 *        acid mappings.
 *
 *        Probably the trickiest thing is that there are two different types of
 *        representations of codons. One in base 5 (N being the 5th base),
 *        providing 0-124 inclusive codon numbers.  These numbers are the ones
 *        going to be principly used in most calculations.
 *
 *        However, it is often very useful to use 0-63 numbers, for example 
 *        in the precise definition of the codon table. 
 *
 *
 */
struct Wise2_CodonTable {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    aa codon_str[125];   
    char * name;     
    } ;  
/* CodonTable defined */ 
#ifndef DYNAMITE_DEFINED_CodonTable
typedef struct Wise2_CodonTable Wise2_CodonTable;
#define CodonTable Wise2_CodonTable
#define DYNAMITE_DEFINED_CodonTable
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  read_CodonTable_file(file)
 *
 * Descrip:    Opens filename, reads it as if a Ewan style
 *             codon table and closes.
 *
 *
 * Arg:        file [READ ] filename to open [char *]
 *
 * Return [OWNER]  A codon-table, NULL if error [CodonTable *]
 *
 */
CodonTable * Wise2_read_CodonTable_file(char * file);
#define read_CodonTable_file Wise2_read_CodonTable_file


/* Function:  read_CodonTable(ifp)
 *
 * Descrip:    reads a codon table from a filestream in Ewan
 *             format.
 *
 *             As Ewan format is really bad and has no start/stop
 *             this will effectively read to the end of the file.
 *             Ooops.
 *
 *
 * Arg:        ifp [READ ] file input [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [CodonTable *]
 *
 */
CodonTable * Wise2_read_CodonTable(FILE * ifp);
#define read_CodonTable Wise2_read_CodonTable


/* Function:  aminoacid_from_seq(ct,seq)
 *
 * Descrip:    Returns the amino acid for this position in the DNA sequence
 *             Takes the pointer +1 and +2 points.
 *
 *             No error checks implemented. Probably a mistake ;)
 *
 *
 * Arg:         ct [READ ] codon table [CodonTable *]
 * Arg:        seq [READ ] pointer to DNA chars [char *]
 *
 * Return [UNKN ]  an amino acid char (A-Z) [aa]
 *
 */
aa Wise2_aminoacid_from_seq(CodonTable * ct,char * seq);
#define aminoacid_from_seq Wise2_aminoacid_from_seq


/* Function:  aminoacid_from_codon(ct,c)
 *
 * Descrip:    returns amino acid for this codon number (NB codon numbers 0-125)
 *
 *
 * Arg:        ct [READ ] codon table [CodonTable *]
 * Arg:         c [READ ] codon number [codon]
 *
 * Return [READ ]  aminoacid that is this codon (X for ambiguous, * for stop) [aa]
 *
 */
aa Wise2_aminoacid_from_codon(CodonTable * ct,codon c);
#define aminoacid_from_codon Wise2_aminoacid_from_codon


/* Function:  aminoacid_no_from_codon(ct,c)
 *
 * Descrip:    a sister function to aminoacid_from_codon:
 *             returns amino acid number (0-26) for this codon number (0-125)
 *
 *
 * Arg:        ct [READ ] codon table [CodonTable *]
 * Arg:         c [READ ] codon number [codon]
 *
 * Return [READ ]  aminoacid number [0-26] for this codon [int]
 *
 */
int Wise2_aminoacid_no_from_codon(CodonTable * ct,codon c);
#define aminoacid_no_from_codon Wise2_aminoacid_no_from_codon


/* Function:  is_stop_codon(c,ct)
 *
 * Descrip:    tells you whether this codon number is really a stop
 *             in this translation table
 *
 *
 * Arg:         c [READ ] codon number [codon]
 * Arg:        ct [READ ] codon table [CodonTable *]
 *
 * Return [UNKN ]  TRUE if is stop, FALSE otherwise [boolean]
 *
 */
boolean Wise2_is_stop_codon(codon c,CodonTable * ct);
#define is_stop_codon Wise2_is_stop_codon


/* Function:  is_non_ambiguous_codon_seq(seq)
 *
 * Descrip:    Tells you if this codon is a real codon
 *
 *
 * Arg:        seq [READ ] pointer to DNA sequence [char *]
 *
 * Return [UNKN ]  TRUE if real codon, FALSE if contains N's [boolean]
 *
 */
boolean Wise2_is_non_ambiguous_codon_seq(char * seq);
#define is_non_ambiguous_codon_seq Wise2_is_non_ambiguous_codon_seq


/* Function:  is_valid_aminoacid(ct,c)
 *
 * Descrip:    Tells you if this letter (c) is recognised as a valid amino acid
 *             in this codon table
 *
 *
 * Arg:        ct [READ ] Codon Table [CodonTable *]
 * Arg:         c [UNKN ] aminoacid [char]
 *
 * Return [UNKN ]  TRUE if valid, FALSE if not. [boolean]
 *
 */
boolean Wise2_is_valid_aminoacid(CodonTable * ct,char c);
#define is_valid_aminoacid Wise2_is_valid_aminoacid


/* Function:  is_valid_base_char(c)
 *
 * Descrip:    Tells you if the letter is A,T,C,G,N (NB, N is ok).
 *
 *
 * Arg:        c [READ ] base [char]
 *
 * Return [UNKN ]  TRUE if (ATGCN) FALSE otherwise [boolean]
 *
 */
boolean Wise2_is_valid_base_char(char c);
#define is_valid_base_char Wise2_is_valid_base_char


/* Function:  codon_from_base4_codon(c)
 *
 * Descrip:    maps a 0-63 codon to a 0-123 codon. Suprisingly useful.
 *
 *
 * Arg:        c [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [codon]
 *
 */
codon Wise2_codon_from_base4_codon(int c);
#define codon_from_base4_codon Wise2_codon_from_base4_codon


/* Function:  base4_codon_from_codon(c)
 *
 * Descrip:    maps a 0-125 codon to a 0-63 codon.
 *
 *             If ambiguous then returns 64 having issued a warning.
 *
 *
 * Arg:        c [READ ] codon 0-125 [codon]
 *
 * Return [UNKN ]  base 4 codon (0-63) [int]
 *
 */
int Wise2_base4_codon_from_codon(codon c);
#define base4_codon_from_codon Wise2_base4_codon_from_codon


/* Function:  has_random_bases(c)
 *
 * Descrip:    Tests to see if this codon number has any N's in it
 *
 *
 * Arg:        c [READ ] codon number 0-124 [codon]
 *
 * Return [UNKN ]  TRUE if has N's , FALSE otherwise [boolean]
 *
 */
boolean Wise2_has_random_bases(codon c);
#define has_random_bases Wise2_has_random_bases


/* Function:  permute_possible_random_bases(c,one,two,three)
 *
 * Descrip:    Bizarely useful function for calculating ambiguity scores.
 *
 *             This takes the codon c, and for each possible base, 
 *             if it is N, replaces it with one, two or three.
 *
 *             If the base is not N, it remains the same
 *
 *
 * Arg:            c [READ ] codon number [codon]
 * Arg:          one [READ ] base to replace first position if N [base]
 * Arg:          two [READ ] base to replace second position if N [base]
 * Arg:        three [READ ] base to replace third position if N [base]
 *
 * Return [UNKN ]  codon number  [codon]
 *
 */
codon Wise2_permute_possible_random_bases(codon c,base one,base two,base three);
#define permute_possible_random_bases Wise2_permute_possible_random_bases


/* Function:  all_bases_from_codon(c,one,two,three)
 *
 * Descrip:    Really an internal function, by useful enough to
 *             encourage outside use.
 *
 *             Takes codon c and breaks it into 3 base-numbers
 *
 *
 * Arg:            c [UNKN ] Undocumented argument [codon]
 * Arg:          one [UNKN ] Undocumented argument [base *]
 * Arg:          two [UNKN ] Undocumented argument [base *]
 * Arg:        three [UNKN ] Undocumented argument [base *]
 *
 */
void Wise2_all_bases_from_codon(codon c,base * one,base * two,base * three);
#define all_bases_from_codon Wise2_all_bases_from_codon


/* Function:  reverse_codon(c)
 *
 * Descrip:    Reverses codon. Takes a forward codon number and
 *             builds the inverted codon number
 *
 *
 * Arg:        c [UNKN ] Undocumented argument [codon]
 *
 * Return [UNKN ]  Undocumented return value [codon]
 *
 */
codon Wise2_reverse_codon(codon c);
#define reverse_codon Wise2_reverse_codon


/* Function:  base_from_codon(c,pos)
 *
 * Descrip:    Probably not the best function to use for this, but 
 *             useful. Takes a codon and with pos being 1,2,3 gives
 *             you the firt,second of third base
 *
 *
 * Arg:          c [UNKN ] Undocumented argument [codon]
 * Arg:        pos [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [base]
 *
 */
base Wise2_base_from_codon(codon c,int pos);
#define base_from_codon Wise2_base_from_codon


/* Function:  codon_from_seq(seq)
 *
 * Descrip:    takes an ASCII coded pointer to a 3 base pair
 *             sequence (it could be the part of a sequence: it only
 *             assummes that the seq points with 3 chars at pos 0,1,2 
 *             in C coordinates from seq. No NULL is required). It 
 *             ives back the codon as made from standard mapping, ie,
 *             25*base_1+5*base_2 + base3 being a number from 0-124 inc.
 *
 *
 * Arg:        seq [UNKN ] pointer to sequence of at least 3 chrs long. [char *]
 *
 * Return [UNKN ]  Undocumented return value [codon]
 *
 */
codon Wise2_codon_from_seq(char * seq);
#define codon_from_seq Wise2_codon_from_seq


/* Function:  base4_codon_from_seq(seq)
 *
 * Descrip:    Sometimes it is more useful to work in base64, ie, 
 *             non N. this functions does the same thing as 
 *             /codon_from_seq but produces a seq being
 *             16*base1 + 4 *base2 + base3
 *
 *
 * Arg:        seq [UNKN ] pointer to sequence of at least 3 chrs long [char *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_base4_codon_from_seq(char * seq);
#define base4_codon_from_seq Wise2_base4_codon_from_seq


/* Function:  char_from_base(b)
 *
 * Descrip:    maps a base number (-04 inc) to A,T,G,C,N
 *
 *
 * Arg:        b [UNKN ] Undocumented argument [base]
 *
 * Return [UNKN ]  Undocumented return value [char]
 *
 */
char Wise2_char_from_base(base b);
#define char_from_base Wise2_char_from_base


/* Function:  base_from_char(c)
 *
 * Descrip:    maps a char (atcgn) to number, 
 *             case insensitive, returns BASE_N
 *             if not atcgn
 *
 *
 * Arg:        c [UNKN ] Undocumented argument [char]
 *
 * Return [UNKN ]  Undocumented return value [base]
 *
 */
base Wise2_base_from_char(char c);
#define base_from_char Wise2_base_from_char


/* Function:  char_complement_base(c)
 *
 * Descrip:    the char equivalent of /complement_base.
 *             this gives the complement in char of a base
 *             in char. Does not check for non ATGCN
 *
 *
 * Arg:        c [UNKN ] Undocumented argument [char]
 *
 * Return [UNKN ]  Undocumented return value [char]
 *
 */
char Wise2_char_complement_base(char c);
#define char_complement_base Wise2_char_complement_base


/* Function:  complement_base(b)
 *
 * Descrip:    gives back the complement as a number
 *             ofthe base (given as a number)
 *
 *
 * Arg:        b [UNKN ] Undocumented argument [base]
 *
 * Return [UNKN ]  Undocumented return value [base]
 *
 */
base Wise2_complement_base(base b);
#define complement_base Wise2_complement_base


/* Function:  four_fold_sites_CodonTable(*ct,seq)
 *
 * Descrip:    returns the number of four fold degenerate
 *             sites in this codon
 *
 *
 * Arg:        *ct [UNKN ] Undocumented argument [CodonTable]
 * Arg:        seq [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_four_fold_sites_CodonTable(CodonTable *ct,char * seq);
#define four_fold_sites_CodonTable Wise2_four_fold_sites_CodonTable


/* Function:  hard_link_CodonTable(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [CodonTable *]
 *
 * Return [UNKN ]  Undocumented return value [CodonTable *]
 *
 */
CodonTable * Wise2_hard_link_CodonTable(CodonTable * obj);
#define hard_link_CodonTable Wise2_hard_link_CodonTable


/* Function:  CodonTable_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [CodonTable *]
 *
 */
CodonTable * Wise2_CodonTable_alloc(void);
#define CodonTable_alloc Wise2_CodonTable_alloc


/* Function:  free_CodonTable(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [CodonTable *]
 *
 * Return [UNKN ]  Undocumented return value [CodonTable *]
 *
 */
CodonTable * Wise2_free_CodonTable(CodonTable * obj);
#define free_CodonTable Wise2_free_CodonTable


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
boolean Wise2_replace_name_CodonTable(CodonTable * obj,char * name);
#define replace_name_CodonTable Wise2_replace_name_CodonTable
char * Wise2_access_name_CodonTable(CodonTable * obj);
#define access_name_CodonTable Wise2_access_name_CodonTable
char * Wise2_alloc_aminoacid_from_seq(CodonTable * ct,char * seq);
#define alloc_aminoacid_from_seq Wise2_alloc_aminoacid_from_seq

#ifdef _cplusplus
}
#endif

#endif
