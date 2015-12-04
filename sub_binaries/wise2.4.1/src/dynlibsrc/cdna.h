#ifndef DYNAMITEcdnaHEADERFILE
#define DYNAMITEcdnaHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "sequence.h"



struct Wise2_cDNA {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Sequence * baseseq;  
    } ;  
/* cDNA defined */ 
#ifndef DYNAMITE_DEFINED_cDNA
typedef struct Wise2_cDNA Wise2_cDNA;
#define cDNA Wise2_cDNA
#define DYNAMITE_DEFINED_cDNA
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  truncate_cDNA(cdna,start,stop)
 *
 * Descrip:    Truncates a cDNA sequence. Basically uses
 *             the /magic_trunc_Sequence function (of course!)
 *
 *             It does not alter cdna, rather it returns a new
 *             sequence with that truncation
 *
 *
 * Arg:         cdna [READ ] cDNA that is truncated [cDNA *]
 * Arg:        start [UNKN ] Undocumented argument [int]
 * Arg:         stop [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [cDNA *]
 *
 */
cDNA * Wise2_truncate_cDNA(cDNA * cdna,int start,int stop);
#define truncate_cDNA Wise2_truncate_cDNA


/* Function:  read_fasta_file_cDNA(filename)
 *
 * Descrip:    Reads a fasta file assumming that it is cDNA. 
 *             Will complain if it is not, and return NULL.
 *
 *
 * Arg:        filename [UNKN ] filename to be opened and read [char *]
 *
 * Return [UNKN ]  Undocumented return value [cDNA *]
 *
 */
cDNA * Wise2_read_fasta_file_cDNA(char * filename);
#define read_fasta_file_cDNA Wise2_read_fasta_file_cDNA


/* Function:  read_fasta_cDNA(ifp)
 *
 * Descrip:    Reads a fasta file assumming that it is cDNA. 
 *             Will complain if it is not, and return NULL.
 *
 *
 * Arg:        ifp [UNKN ] file point to be read from [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [cDNA *]
 *
 */
cDNA * Wise2_read_fasta_cDNA(FILE * ifp);
#define read_fasta_cDNA Wise2_read_fasta_cDNA


/* Function:  read_efetch_cDNA(estr)
 *
 * Descrip:    Reads a efetch specified query
 *             Uses, of course /read_efetch_Sequence
 *
 *
 * Arg:        estr [READ ] efetch string which is read [char *]
 *
 * Return [UNKN ]  Undocumented return value [cDNA *]
 *
 */
cDNA * Wise2_read_efetch_cDNA(char * estr);
#define read_efetch_cDNA Wise2_read_efetch_cDNA


/* Function:  read_SRS_cDNA(srsquery)
 *
 * Descrip:    Reads a SRS sequence using srs4 syntax.
 *             Uses, of course, /read_SRS_Sequence
 *
 *
 *
 * Arg:        srsquery [READ ] string query representing SRS name [char *]
 *
 * Return [UNKN ]  Undocumented return value [cDNA *]
 *
 */
cDNA * Wise2_read_SRS_cDNA(char * srsquery);
#define read_SRS_cDNA Wise2_read_SRS_cDNA


/* Function:  cDNA_name(cdna)
 *
 * Descrip:    Returns the name of the cDNA
 *
 *
 * Arg:        cdna [UNKN ] Undocumented argument [cDNA *]
 *
 * Return [UNKN ]  Undocumented return value [char *]
 *
 */
char * Wise2_cDNA_name(cDNA * cdna);
#define cDNA_name Wise2_cDNA_name


/* Function:  cDNA_length(cdna)
 *
 * Descrip:    Returns the length of the cDNA
 *
 *
 * Arg:        cdna [UNKN ] Undocumented argument [cDNA *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_cDNA_length(cDNA * cdna);
#define cDNA_length Wise2_cDNA_length


/* Function:  cDNA_seqchar(cdna,pos)
 *
 * Descrip:    Returns sequence character at this position.
 *
 *
 * Arg:        cdna [UNKN ] cDNA [cDNA *]
 * Arg:         pos [UNKN ] position in cDNA to get char [int]
 *
 * Return [UNKN ]  Undocumented return value [char]
 *
 */
char Wise2_cDNA_seqchar(cDNA * cdna,int pos);
#define cDNA_seqchar Wise2_cDNA_seqchar


/* Function:  cDNA_from_Sequence(seq)
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
 * Arg:        seq [OWNER] Sequence to make cDNA from [Sequence *]
 *
 * Return [UNKN ]  Undocumented return value [cDNA *]
 *
 */
cDNA * Wise2_cDNA_from_Sequence(Sequence * seq);
#define cDNA_from_Sequence Wise2_cDNA_from_Sequence


/* Function:  hard_link_cDNA(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [cDNA *]
 *
 * Return [UNKN ]  Undocumented return value [cDNA *]
 *
 */
cDNA * Wise2_hard_link_cDNA(cDNA * obj);
#define hard_link_cDNA Wise2_hard_link_cDNA


/* Function:  cDNA_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [cDNA *]
 *
 */
cDNA * Wise2_cDNA_alloc(void);
#define cDNA_alloc Wise2_cDNA_alloc


/* Function:  free_cDNA(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [cDNA *]
 *
 * Return [UNKN ]  Undocumented return value [cDNA *]
 *
 */
cDNA * Wise2_free_cDNA(cDNA * obj);
#define free_cDNA Wise2_free_cDNA


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
boolean Wise2_replace_baseseq_cDNA(cDNA * obj,Sequence * baseseq);
#define replace_baseseq_cDNA Wise2_replace_baseseq_cDNA
Sequence * Wise2_access_baseseq_cDNA(cDNA * obj);
#define access_baseseq_cDNA Wise2_access_baseseq_cDNA

#ifdef _cplusplus
}
#endif

#endif
