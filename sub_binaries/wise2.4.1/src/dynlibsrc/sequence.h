#ifndef DYNAMITEsequenceHEADERFILE
#define DYNAMITEsequenceHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif

#include "wisebase.h"
#include "codon.h"

#ifdef LINUX
#include "posix.h"
#endif

#define SEQUENCEBLOCK 8000

#define SEQUENCE_STARTSIZE 512
#define SEQUENCE_REALLOC_LINEAR (4096*4)

enum SequenceType {
SEQUENCE_UNKNOWN = 64,
SEQUENCE_PROTEIN,
SEQUENCE_DNA,
SEQUENCE_CDNA,
SEQUENCE_GENOMIC,
SEQUENCE_EST,
SEQUENCE_RNA };

#define is_dna_SequenceType(type) (type == SEQUENCE_DNA || type == SEQUENCE_CDNA || type == SEQUENCE_GENOMIC || type == SEQUENCE_EST ? TRUE : FALSE)
#define is_rna_SequenceType(type) (type == SEQUENCE_RNA ? TRUE : FALSE)
#define is_protein_SequenceType(type) (type == SEQUENCE_PROTEIN ? TRUE : FALSE )

#define is_dna_Sequence(seq) (is_dna_SequenceType(seq->type))
#define is_rna_Sequence(seq) (is_rna_SequenceType(seq->type))
#define is_protein_Sequence(seq) (is_protein_SequenceType(seq->type))

#define SequenceSetLISTLENGTH 64
/* Object Sequence
 *
 * Descrip: This object is the basic sequence object,
 *        trying to hold little more than the 
 *        name and sequence of the DNA/protein. 
 *
 *        The len/maxlen is the actual length
 *        of the sequence (strlen(obj->seq)) and
 *        amount of memory allocated in obj->seq 
 *        mainly for parsing purposes.
 *
 *
 *
 */
struct Wise2_Sequence {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * name;    /*  name of the sequence */ 
    char * seq; /*  actual sequence */ 
    int len;    /*  length of the sequence */ 
    int maxlen; /*  internal counter, indicating how much space in seq there is */ 
    int offset; /*  start (in bio-coords) of the sequence. Not called start due to weird legacy. */ 
    int end;    /*  end (in bio-coords == C coords) of the sequence */ 
    char type;  /*  guess of protein/dna type */ 
    int tax_id; /*  taxonimic id of this */ 
    double weight;   
    char * desc;    /*  description - often this will be NULL */ 
    } ;  
/* Sequence defined */ 
#ifndef DYNAMITE_DEFINED_Sequence
typedef struct Wise2_Sequence Wise2_Sequence;
#define Sequence Wise2_Sequence
#define DYNAMITE_DEFINED_Sequence
#endif


/* Object SequenceSet
 *
 * Descrip: A list of sequences. Not a database (you should
 *        use the db stuff for that!). But useful anyway
 *
 *
 */
struct Wise2_SequenceSet {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Sequence ** set;     
    int len;/* len for above set  */ 
    int maxlen; /* maxlen for above set */ 
    } ;  
/* SequenceSet defined */ 
#ifndef DYNAMITE_DEFINED_SequenceSet
typedef struct Wise2_SequenceSet Wise2_SequenceSet;
#define SequenceSet Wise2_SequenceSet
#define DYNAMITE_DEFINED_SequenceSet
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  new_Sequence_from_strings(name,seq)
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
 * Arg:        name [READ ] name of sequence, memory is allocated for it. [char *]
 * Arg:         seq [READ ] char * of sequence, memory is allocated for it. [char *]
 *
 * Return [UNKN ]  Undocumented return value [Sequence *]
 *
 */
Sequence * Wise2_new_Sequence_from_strings(char * name,char * seq);
#define new_Sequence_from_strings Wise2_new_Sequence_from_strings


/* Function:  looks_like_accession(name)
 *
 * Descrip:    Returns true if name looks like [A-Za-z]+[0-9]+
 *             This should be an accession number 
 *
 *
 * Arg:        name [READ ] name to be tested [char *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_looks_like_accession(char * name);
#define looks_like_accession Wise2_looks_like_accession


/* Function:  make_len_type_Sequence(seq)
 *
 * Descrip:    makes seq->len and seq->end match the seq->seq
 *             length number. 
 *
 *             It also checks the type of the sequence with
 *             /best_guess_type
 *
 *
 * Arg:        seq [RW   ] Sequence object [Sequence *]
 *
 */
void Wise2_make_len_type_Sequence(Sequence * seq);
#define make_len_type_Sequence Wise2_make_len_type_Sequence


/* Function:  best_guess_type(seq)
 *
 * Descrip:    Guesses DNA or protein, by adding all
 *             the A,T,G,C up and if len < 300 && > 95% or 
 *             len > 300 && > 75% then considers
 *             it to be DNA. NB - Ns now counted.
 *
 *
 * Arg:        seq [READ ] Sequence to be guessed [Sequence *]
 *
 * Return [OWNER]  SEQUENCE_DNA or SEQUENCE_PROTEIN [int]
 *
 */
int Wise2_best_guess_type(Sequence * seq);
#define best_guess_type Wise2_best_guess_type


/* Function:  Sequence_type_to_string(type)
 *
 * Descrip:    Converts sequence type (SEQUENCE_*) to a string
 *
 *
 * Arg:        type [UNKN ] type eg SEQUENCE_PROTEIN [int]
 *
 * Return [UNKN ]  Undocumented return value [char *]
 *
 */
char * Wise2_Sequence_type_to_string(int type);
#define Sequence_type_to_string Wise2_Sequence_type_to_string


/* Function:  uppercase_Sequence(seq)
 *
 * Descrip:    makes all the sequence uppercase
 *
 *
 * Arg:        seq [RW   ] Sequence to be uppercased [Sequence *]
 *
 */
void Wise2_uppercase_Sequence(Sequence * seq);
#define uppercase_Sequence Wise2_uppercase_Sequence


/* Function:  force_to_dna_Sequence(seq,fraction,number_of_conver)
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
 * Arg:                     seq [RW   ] sequence object read and converted  [Sequence *]
 * Arg:                fraction [READ ] number 0..1 for tolerance of conversion [double]
 * Arg:        number_of_conver [WRITE] number of conversions actually made [int *]
 *
 * Return [READ ]  TRUE for conversion to DNA, FALSE if not [boolean]
 *
 */
boolean Wise2_force_to_dna_Sequence(Sequence * seq,double fraction,int * number_of_conver);
#define force_to_dna_Sequence Wise2_force_to_dna_Sequence


/* Function:  is_reversed_Sequence(seq)
 *
 * Descrip:    Currently the sequence object stores 
 *             reversed sequences as start > end.
 *
 *             This tests that and returns true if it is
 *
 *
 * Arg:        seq [READ ] sequence to test [Sequence *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_is_reversed_Sequence(Sequence * seq);
#define is_reversed_Sequence Wise2_is_reversed_Sequence


/* Function:  translate_Sequence(dna,ct)
 *
 * Descrip:    This translates a DNA sequence to a protein.
 *             It assummes that it starts at first residue
 *             (use trunc_Sequence to chop a sequence up).
 *
 *
 * Arg:        dna [READ ] DNA sequence to be translated [Sequence *]
 * Arg:         ct [READ ] Codon table to do codon->aa mapping [CodonTable *]
 *
 * Return [OWNER]  new protein sequence [Sequence *]
 *
 */
Sequence * Wise2_translate_Sequence(Sequence * dna,CodonTable * ct);
#define translate_Sequence Wise2_translate_Sequence


/* Function:  reverse_complement_Sequence(seq)
 *
 * Descrip:    This both complements and reverses a sequence,
 *             - a common wish!
 *
 *             The start/end are correct with respect to the start/end
 *             of the sequence (ie start = end, end = start).
 *
 *
 * Arg:        seq [READ ] Sequence to that is used to reverse (makes a new Sequence) [Sequence *]
 *
 * Return [OWNER]  new Sequence which is reversed [Sequence *]
 *
 */
Sequence * Wise2_reverse_complement_Sequence(Sequence * seq);
#define reverse_complement_Sequence Wise2_reverse_complement_Sequence


/* Function:  magic_trunc_Sequence(seq,start,end)
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
 * Arg:          seq [READ ] sequence that is the source to be truncated [Sequence *]
 * Arg:        start [READ ] start point [int]
 * Arg:          end [READ ] end point [int]
 *
 * Return [OWNER]  new Sequence which is truncated/reversed [Sequence *]
 *
 */
Sequence * Wise2_magic_trunc_Sequence(Sequence * seq,int start,int end);
#define magic_trunc_Sequence Wise2_magic_trunc_Sequence


/* Function:  trunc_Sequence(seq,start,end)
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
 * Arg:          seq [READ ] object holding the sequence to be truncated [Sequence *]
 * Arg:        start [READ ] start point of truncation [int]
 * Arg:          end [READ ] end point of truncation [int]
 *
 * Return [OWNER]  newly allocated sequence structure [Sequence *]
 *
 */
Sequence * Wise2_trunc_Sequence(Sequence * seq,int start,int end);
#define trunc_Sequence Wise2_trunc_Sequence


/* Function:  read_SRS_db_Sequence(datastring,srsstring)
 *
 * Descrip:    A function for you to easily specify the sequence name
 *             and the database separately. Just concatonates the two
 *             strings with : betwqeen them. Therefore you should use
 *             "swisprot-id" for example as your datastring.
 *
 *             calls /read_SRS_Sequence
 *
 *
 * Arg:        datastring [READ ] string representing the database (swissprot-id) [char *]
 * Arg:         srsstring [READ ] string for the name (eg, ROA1_HUMAN) [char *]
 *
 * Return [UNKN ]  Undocumented return value [Sequence *]
 *
 */
Sequence * Wise2_read_SRS_db_Sequence(char * datastring,char * srsstring);
#define read_SRS_db_Sequence Wise2_read_SRS_db_Sequence


/* Function:  read_SRS_Sequence(srsstring)
 *
 * Descrip:    reads SRS specified sequence. calls popoen
 *             with getz -f using srs4 syntax. Will only read
 *             the first sequence if there is more than one in the 
 *             SRS spec, and does not warn you about additional 
 *             sequences
 *
 *
 * Arg:        srsstring [READ ] srs spec'd string swissprot-id:ROA1_HUMAN [char *]
 *
 * Return [UNKN ]  Undocumented return value [Sequence *]
 *
 */
Sequence * Wise2_read_SRS_Sequence(char * srsstring);
#define read_SRS_Sequence Wise2_read_SRS_Sequence


/* Function:  read_efetch_Sequence(efetchstring)
 *
 * Descrip:    reads efetch specificed sequence. calls popen to
 *             efetch. A hack around accession numbers so that if the 
 *             thing looks like WP:acc number, calls it with -a...
 *             otherwise assummes you have both database and name in the
 *             efetchstring
 *
 *
 * Arg:        efetchstring [READ ] efetch valid string [char *]
 *
 * Return [UNKN ]  Undocumented return value [Sequence *]
 *
 */
Sequence * Wise2_read_efetch_Sequence(char * efetchstring);
#define read_efetch_Sequence Wise2_read_efetch_Sequence


/* Function:  read_fasta_file_Sequence(filename)
 *
 * Descrip:    Just a call
 *               a) open filename
 *               b) read sequence with /read_fasta_Sequence
 *               c) close file.
 *
 *
 * Arg:        filename [READ ] filename to open  [char *]
 *
 * Return [UNKN ]  Undocumented return value [Sequence *]
 *
 */
Sequence * Wise2_read_fasta_file_Sequence(char * filename);
#define read_fasta_file_Sequence Wise2_read_fasta_file_Sequence


/* Function:  read_Sequence_EMBL_seq(buffer,maxlen,ifp)
 *
 * Descrip:    reads the sequence part of an EMBL file.
 *
 *             This function can either take a file which 
 *             starts
 *
 *
 *
 * Arg:        buffer [RW   ] buffer containing the first line. [char *]
 * Arg:        maxlen [READ ] length of buffer [int]
 * Arg:           ifp [READ ] input file to read from [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [Sequence *]
 *
 */
Sequence * Wise2_read_Sequence_EMBL_seq(char * buffer,int maxlen,FILE * ifp);
#define read_Sequence_EMBL_seq Wise2_read_Sequence_EMBL_seq


/* Function:  read_fasta_Sequence(ifp)
 *
 * Descrip:    reads a fasta file assumming pretty 
 *             standard sanity for the file layout
 *
 *
 * Arg:        ifp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [Sequence *]
 *
 */
Sequence * Wise2_read_fasta_Sequence(FILE * ifp);
#define read_fasta_Sequence Wise2_read_fasta_Sequence


/* Function:  read_old_fasta_Sequence(ifp)
 *
 * Descrip:    reads the fasta file: format is
 *
 *             >name
 *             sequence
 *
 *             allocates a structure and puts in the
 *             sequence. Calls /make_len_type_Sequence to
 *             check type and length.
 *
 *             It leaves the '>' on the next fasta sequence
 *             for multiple sequence reading
 *
 *
 * Arg:        ifp [READ ] input file to read from [FILE *]
 *
 * Return [OWNER]  new Sequence structure  [Sequence *]
 *
 */
Sequence * Wise2_read_old_fasta_Sequence(FILE * ifp);
#define read_old_fasta_Sequence Wise2_read_old_fasta_Sequence


/* Function:  show_Sequence_residue_list(seq,start,end,ofp)
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
 * Arg:          seq [READ ] Sequence to show [Sequence *]
 * Arg:        start [READ ] start of list [int]
 * Arg:          end [READ ] end of list [int]
 * Arg:          ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_show_Sequence_residue_list(Sequence * seq,int start,int end,FILE * ofp);
#define show_Sequence_residue_list Wise2_show_Sequence_residue_list


/* Function:  add_string_to_Sequence(seq,more)
 *
 * Descrip:    New add_string_to_Sequence which should be much better at
 *             handling memory update
 *
 *
 * Arg:         seq [UNKN ] Undocumented argument [Sequence *]
 * Arg:        more [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_string_to_Sequence(Sequence * seq,char * more);
#define add_string_to_Sequence Wise2_add_string_to_Sequence


/* Function:  empty_Sequence_from_dynamic_memory(name)
 *
 * Descrip:    Only allocates sequence structure and name
 *
 *
 * Arg:        name [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [Sequence *]
 *
 */
Sequence * Wise2_empty_Sequence_from_dynamic_memory(char * name);
#define empty_Sequence_from_dynamic_memory Wise2_empty_Sequence_from_dynamic_memory


/* Function:  Sequence_alloc_len(len)
 *
 * Descrip:    allocates sequence structure with enough
 *             length in char for len sequence.
 *
 *
 * Arg:        len [READ ] length of blank sequene space [int]
 *
 * Return [UNKN ]  Undocumented return value [Sequence *]
 *
 */
Sequence * Wise2_Sequence_alloc_len(int len);
#define Sequence_alloc_len Wise2_Sequence_alloc_len


/* Function:  write_fasta_Sequence(seq,ofp)
 *
 * Descrip:    writes a fasta file of the form
 *             >name
 *             Sequence
 *
 *
 * Arg:        seq [READ ] sequence to be written [Sequence *]
 * Arg:        ofp [UNKN ] file to write to [FILE *]
 *
 */
void Wise2_write_fasta_Sequence(Sequence * seq,FILE * ofp);
#define write_fasta_Sequence Wise2_write_fasta_Sequence


/* Function:  read_fasta_SequenceSet(ifp)
 *
 * Descrip:    reads a fasta file as a sequence set
 *
 *
 * Arg:        ifp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [SequenceSet *]
 *
 */
SequenceSet * Wise2_read_fasta_SequenceSet(FILE * ifp);
#define read_fasta_SequenceSet Wise2_read_fasta_SequenceSet


/* Function:  read_fasta_SequenceSet_file(filename)
 *
 * Descrip:    opens file and reads in sequence set
 *
 *
 * Arg:        filename [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [SequenceSet *]
 *
 */
SequenceSet * Wise2_read_fasta_SequenceSet_file(char * filename);
#define read_fasta_SequenceSet_file Wise2_read_fasta_SequenceSet_file


/* Function:  hard_link_Sequence(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Sequence *]
 *
 * Return [UNKN ]  Undocumented return value [Sequence *]
 *
 */
Sequence * Wise2_hard_link_Sequence(Sequence * obj);
#define hard_link_Sequence Wise2_hard_link_Sequence


/* Function:  Sequence_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Sequence *]
 *
 */
Sequence * Wise2_Sequence_alloc(void);
#define Sequence_alloc Wise2_Sequence_alloc


/* Function:  free_Sequence(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Sequence *]
 *
 * Return [UNKN ]  Undocumented return value [Sequence *]
 *
 */
Sequence * Wise2_free_Sequence(Sequence * obj);
#define free_Sequence Wise2_free_Sequence


/* Function:  add_SequenceSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [SequenceSet *]
 * Arg:        add [OWNER] Object to add to the list [Sequence *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_SequenceSet(SequenceSet * obj,Sequence * add);
#define add_SequenceSet Wise2_add_SequenceSet


/* Function:  flush_SequenceSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [SequenceSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_SequenceSet(SequenceSet * obj);
#define flush_SequenceSet Wise2_flush_SequenceSet


/* Function:  SequenceSet_alloc_std(void)
 *
 * Descrip:    Equivalent to SequenceSet_alloc_len(SequenceSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SequenceSet *]
 *
 */
SequenceSet * Wise2_SequenceSet_alloc_std(void);
#define SequenceSet_alloc_std Wise2_SequenceSet_alloc_std


/* Function:  SequenceSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [SequenceSet *]
 *
 */
SequenceSet * Wise2_SequenceSet_alloc_len(int len);
#define SequenceSet_alloc_len Wise2_SequenceSet_alloc_len


/* Function:  hard_link_SequenceSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SequenceSet *]
 *
 * Return [UNKN ]  Undocumented return value [SequenceSet *]
 *
 */
SequenceSet * Wise2_hard_link_SequenceSet(SequenceSet * obj);
#define hard_link_SequenceSet Wise2_hard_link_SequenceSet


/* Function:  SequenceSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SequenceSet *]
 *
 */
SequenceSet * Wise2_SequenceSet_alloc(void);
#define SequenceSet_alloc Wise2_SequenceSet_alloc


/* Function:  free_SequenceSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SequenceSet *]
 *
 * Return [UNKN ]  Undocumented return value [SequenceSet *]
 *
 */
SequenceSet * Wise2_free_SequenceSet(SequenceSet * obj);
#define free_SequenceSet Wise2_free_SequenceSet


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
boolean Wise2_replace_len_Sequence(Sequence * obj,int len);
#define replace_len_Sequence Wise2_replace_len_Sequence
boolean Wise2_replace_offset_Sequence(Sequence * obj,int offset);
#define replace_offset_Sequence Wise2_replace_offset_Sequence
int Wise2_access_offset_Sequence(Sequence * obj);
#define access_offset_Sequence Wise2_access_offset_Sequence
int Wise2_length_set_SequenceSet(SequenceSet * obj);
#define length_set_SequenceSet Wise2_length_set_SequenceSet
boolean Wise2_replace_end_Sequence(Sequence * obj,int end);
#define replace_end_Sequence Wise2_replace_end_Sequence
boolean Wise2_replace_name_Sequence(Sequence * obj,char * name);
#define replace_name_Sequence Wise2_replace_name_Sequence
int Wise2_access_end_Sequence(Sequence * obj);
#define access_end_Sequence Wise2_access_end_Sequence
boolean Wise2_replace_seq_Sequence(Sequence * obj,char * seq);
#define replace_seq_Sequence Wise2_replace_seq_Sequence
boolean Wise2_replace_type_Sequence(Sequence * obj,char type);
#define replace_type_Sequence Wise2_replace_type_Sequence
int Wise2_access_len_Sequence(Sequence * obj);
#define access_len_Sequence Wise2_access_len_Sequence
char Wise2_access_type_Sequence(Sequence * obj);
#define access_type_Sequence Wise2_access_type_Sequence
int Wise2_access_maxlen_Sequence(Sequence * obj);
#define access_maxlen_Sequence Wise2_access_maxlen_Sequence
boolean Wise2_replace_tax_id_Sequence(Sequence * obj,int tax_id);
#define replace_tax_id_Sequence Wise2_replace_tax_id_Sequence
Sequence * Wise2_access_set_SequenceSet(SequenceSet * obj,int i);
#define access_set_SequenceSet Wise2_access_set_SequenceSet
int Wise2_access_tax_id_Sequence(Sequence * obj);
#define access_tax_id_Sequence Wise2_access_tax_id_Sequence
char * Wise2_access_seq_Sequence(Sequence * obj);
#define access_seq_Sequence Wise2_access_seq_Sequence
boolean Wise2_replace_weight_Sequence(Sequence * obj,double weight);
#define replace_weight_Sequence Wise2_replace_weight_Sequence
char * Wise2_access_name_Sequence(Sequence * obj);
#define access_name_Sequence Wise2_access_name_Sequence
char * Wise2_access_desc_Sequence(Sequence * obj);
#define access_desc_Sequence Wise2_access_desc_Sequence
double Wise2_access_weight_Sequence(Sequence * obj);
#define access_weight_Sequence Wise2_access_weight_Sequence
boolean Wise2_replace_maxlen_Sequence(Sequence * obj,int maxlen);
#define replace_maxlen_Sequence Wise2_replace_maxlen_Sequence
boolean Wise2_replace_desc_Sequence(Sequence * obj,char * desc);
#define replace_desc_Sequence Wise2_replace_desc_Sequence
boolean Wise2_add_string_to_Sequence_old(Sequence * seq,char * more);
#define add_string_to_Sequence_old Wise2_add_string_to_Sequence_old
Sequence * Wise2_Sequence_from_static_memory (char * name,char * seq);
#define Sequence_from_static_memory  Wise2_Sequence_from_static_memory 
Sequence * Wise2_Sequence_from_dynamic_memory(char * name,char * seq);
#define Sequence_from_dynamic_memory Wise2_Sequence_from_dynamic_memory
void Wise2_swap_SequenceSet(Sequence ** list,int i,int j) ;
#define swap_SequenceSet Wise2_swap_SequenceSet
void Wise2_qsort_SequenceSet(Sequence ** list,int left,int right,int (*comp)(Sequence * ,Sequence * ));
#define qsort_SequenceSet Wise2_qsort_SequenceSet
void Wise2_sort_SequenceSet(SequenceSet * obj,int (*comp)(Sequence *, Sequence *));
#define sort_SequenceSet Wise2_sort_SequenceSet
boolean Wise2_expand_SequenceSet(SequenceSet * obj,int len);
#define expand_SequenceSet Wise2_expand_SequenceSet

#ifdef _cplusplus
}
#endif

#endif
