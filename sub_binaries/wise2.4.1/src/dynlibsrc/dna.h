#ifndef DYNAMITEdnaHEADERFILE
#define DYNAMITEdnaHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "sequence.h"



struct Wise2_DNA {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Sequence * baseseq;  
    } ;  
/* DNA defined */ 
#ifndef DYNAMITE_DEFINED_DNA
typedef struct Wise2_DNA Wise2_DNA;
#define DNA Wise2_DNA
#define DYNAMITE_DEFINED_DNA
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  truncate_DNA(dna,start,stop)
 *
 * Descrip:    Truncates a DNA sequence. Basically uses
 *             the /trunc_Sequence function (of course!)
 *
 *             It does not alter dna, rather it returns a new
 *             sequence with that truncation
 *
 *
 * Arg:          dna [READ ] DNA that is truncated [DNA *]
 * Arg:        start [UNKN ] Undocumented argument [int]
 * Arg:         stop [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [DNA *]
 *
 */
DNA * Wise2_truncate_DNA(DNA * dna,int start,int stop);
#define truncate_DNA Wise2_truncate_DNA


/* Function:  read_fasta_file_DNA (filename)
 *
 * Descrip:    Reads a fasta file assumming that it is DNA. 
 *             Will complain if it is not, and return NULL.
 *
 *
 * Arg:        filename [UNKN ] filename to be opened and read [char *]
 *
 * Return [UNKN ]  Undocumented return value [DNA *]
 *
 */
DNA * Wise2_read_fasta_file_DNA (char * filename);
#define read_fasta_file_DNA  Wise2_read_fasta_file_DNA 


/* Function:  read_fasta_DNA (ifp)
 *
 * Descrip:    Reads a fasta file assumming that it is DNA. 
 *             Will complain if it is not, and return NULL.
 *
 *
 * Arg:        ifp [UNKN ] file point to be read from [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [DNA *]
 *
 */
DNA * Wise2_read_fasta_DNA (FILE * ifp);
#define read_fasta_DNA  Wise2_read_fasta_DNA 


/* Function:  read_efetch_DNA(estr)
 *
 * Descrip:    Reads a efetch specified query
 *             Uses, of course /read_efetch_Sequence
 *
 *
 * Arg:        estr [READ ] efetch string which is read [char *]
 *
 * Return [UNKN ]  Undocumented return value [DNA *]
 *
 */
DNA * Wise2_read_efetch_DNA(char * estr);
#define read_efetch_DNA Wise2_read_efetch_DNA


/* Function:  read_SRS_DNA(srsquery)
 *
 * Descrip:    Reads a SRS sequence using srs4 syntax.
 *             Uses, of course, /read_SRS_Sequence
 *
 *
 *
 * Arg:        srsquery [READ ] string query representing SRS name [char *]
 *
 * Return [UNKN ]  Undocumented return value [DNA *]
 *
 */
DNA * Wise2_read_SRS_DNA(char * srsquery);
#define read_SRS_DNA Wise2_read_SRS_DNA


/* Function:  DNA_name (dna)
 *
 * Descrip:    Returns the name of the DNA
 *
 *
 * Arg:        dna [UNKN ] Undocumented argument [DNA *]
 *
 * Return [UNKN ]  Undocumented return value [char *]
 *
 */
char * Wise2_DNA_name (DNA * dna);
#define DNA_name  Wise2_DNA_name 


/* Function:  DNA_length (dna)
 *
 * Descrip:    Returns the length of the DNA
 *
 *
 * Arg:        dna [UNKN ] Undocumented argument [DNA *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_DNA_length (DNA * dna);
#define DNA_length  Wise2_DNA_length 


/* Function:  DNA_seqchar(dna,pos)
 *
 * Descrip:    Returns sequence character at this position.
 *
 *
 * Arg:        dna [UNKN ] DNA [DNA *]
 * Arg:        pos [UNKN ] position in DNA to get char [int]
 *
 * Return [UNKN ]  Undocumented return value [char]
 *
 */
char Wise2_DNA_seqchar(DNA * dna,int pos);
#define DNA_seqchar Wise2_DNA_seqchar


/* Function:  DNA_from_Sequence(seq)
 *
 * Descrip:    makes a new DNA from a Sequence. It 
 *             owns the Sequence memory, ie will attempt a /free_Sequence
 *             on the structure when /free_DNA is called
 *
 *             If you want to give this DNA this Sequence and
 *             forget about it, then just hand it this sequence and set
 *             seq to NULL (no need to free it). If you intend to use 
 *             the sequence object elsewhere outside of the DNA datastructure
 *             then use DNA_from_Sequence(/hard_link_Sequence(seq))
 *
 *
 *
 * Arg:        seq [UNKN ] Sequence to make DNA from [Sequence *]
 *
 * Return [UNKN ]  Undocumented return value [DNA *]
 *
 */
DNA * Wise2_DNA_from_Sequence(Sequence * seq);
#define DNA_from_Sequence Wise2_DNA_from_Sequence


/* Function:  hard_link_DNA(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DNA *]
 *
 * Return [UNKN ]  Undocumented return value [DNA *]
 *
 */
DNA * Wise2_hard_link_DNA(DNA * obj);
#define hard_link_DNA Wise2_hard_link_DNA


/* Function:  DNA_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DNA *]
 *
 */
DNA * Wise2_DNA_alloc(void);
#define DNA_alloc Wise2_DNA_alloc


/* Function:  free_DNA(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DNA *]
 *
 * Return [UNKN ]  Undocumented return value [DNA *]
 *
 */
DNA * Wise2_free_DNA(DNA * obj);
#define free_DNA Wise2_free_DNA


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/

#ifdef _cplusplus
}
#endif

#endif
