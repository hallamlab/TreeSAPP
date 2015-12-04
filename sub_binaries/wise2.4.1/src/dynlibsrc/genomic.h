#ifndef DYNAMITEgenomicHEADERFILE
#define DYNAMITEgenomicHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "sequence.h"

#define GenomicLISTLENGTH 128

struct Wise2_GenomicRepeat {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int start;   
    int end;     
    char * type;     
    } ;  
/* GenomicRepeat defined */ 
#ifndef DYNAMITE_DEFINED_GenomicRepeat
typedef struct Wise2_GenomicRepeat Wise2_GenomicRepeat;
#define GenomicRepeat Wise2_GenomicRepeat
#define DYNAMITE_DEFINED_GenomicRepeat
#endif


struct Wise2_Genomic {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Sequence * baseseq;  
    GenomicRepeat ** repeat;     
    int len;/* len for above repeat  */ 
    int maxlen; /* maxlen for above repeat */ 
    } ;  
/* Genomic defined */ 
#ifndef DYNAMITE_DEFINED_Genomic
typedef struct Wise2_Genomic Wise2_Genomic;
#define Genomic Wise2_Genomic
#define DYNAMITE_DEFINED_Genomic
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  truncate_Genomic(gen,start,stop)
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
 * Arg:          gen [READ ] Genomic that is truncated [Genomic *]
 * Arg:        start [UNKN ] Undocumented argument [int]
 * Arg:         stop [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [Genomic *]
 *
 */
Genomic * Wise2_truncate_Genomic(Genomic * gen,int start,int stop);
#define truncate_Genomic Wise2_truncate_Genomic


/* Function:  reverse_complement_Genomic(gen)
 *
 * Descrip:    Reverse Complements s Genomic sequence. 
 *
 *             Handles repeat information correctly
 *
 *
 * Arg:        gen [READ ] Genomic that is revomcp [Genomic *]
 *
 * Return [UNKN ]  Undocumented return value [Genomic *]
 *
 */
Genomic * Wise2_reverse_complement_Genomic(Genomic * gen);
#define reverse_complement_Genomic Wise2_reverse_complement_Genomic


/* Function:  read_fasta_file_Genomic(filename,length_of_N)
 *
 * Descrip:    Reads a fasta file assumming that it is Genomic. 
 *             Will complain if it is not, and return NULL.
 *
 *
 * Arg:           filename [UNKN ] filename to be opened and read [char *]
 * Arg:        length_of_N [UNKN ] length of N to be considered repeat. -1 means none [int]
 *
 * Return [UNKN ]  Undocumented return value [Genomic *]
 *
 */
Genomic * Wise2_read_fasta_file_Genomic(char * filename,int length_of_N);
#define read_fasta_file_Genomic Wise2_read_fasta_file_Genomic


/* Function:  read_fasta_Genomic(ifp,length_of_N)
 *
 * Descrip:    Reads a fasta file assumming that it is Genomic. 
 *             Will complain if it is not, and return NULL.
 *
 *
 * Arg:                ifp [UNKN ] file point to be read from [FILE *]
 * Arg:        length_of_N [UNKN ] length of N to be considered repeat. -1 means none [int]
 *
 * Return [UNKN ]  Undocumented return value [Genomic *]
 *
 */
Genomic * Wise2_read_fasta_Genomic(FILE * ifp,int length_of_N);
#define read_fasta_Genomic Wise2_read_fasta_Genomic


/* Function:  Genomic_name(gen)
 *
 * Descrip:    Returns the name of the Genomic
 *
 *
 * Arg:        gen [UNKN ] Undocumented argument [Genomic *]
 *
 * Return [UNKN ]  Undocumented return value [char *]
 *
 */
char * Wise2_Genomic_name(Genomic * gen);
#define Genomic_name Wise2_Genomic_name


/* Function:  Genomic_length(gen)
 *
 * Descrip:    Returns the length of the Genomic
 *
 *
 * Arg:        gen [UNKN ] Undocumented argument [Genomic *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_Genomic_length(Genomic * gen);
#define Genomic_length Wise2_Genomic_length


/* Function:  Genomic_seqchar(gen,pos)
 *
 * Descrip:    Returns sequence character at this position.
 *
 *
 * Arg:        gen [UNKN ] Genomic [Genomic *]
 * Arg:        pos [UNKN ] position in Genomic to get char [int]
 *
 * Return [UNKN ]  Undocumented return value [char]
 *
 */
char Wise2_Genomic_seqchar(Genomic * gen,int pos);
#define Genomic_seqchar Wise2_Genomic_seqchar


/* Function:  Genomic_from_Sequence_Nheuristic(seq,length_of_N)
 *
 * Descrip:    makes a new genomic from a Sequence, but
 *             assummes that all the N runs greater than
 *             a certain level are actually repeats.
 *
 *
 * Arg:                seq [UNKN ] Undocumented argument [Sequence *]
 * Arg:        length_of_N [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [Genomic *]
 *
 */
Genomic * Wise2_Genomic_from_Sequence_Nheuristic(Sequence * seq,int length_of_N);
#define Genomic_from_Sequence_Nheuristic Wise2_Genomic_from_Sequence_Nheuristic


/* Function:  Genomic_from_Sequence(seq)
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
 * Arg:        seq [OWNER] Sequence to make genomic from [Sequence *]
 *
 * Return [UNKN ]  Undocumented return value [Genomic *]
 *
 */
Genomic * Wise2_Genomic_from_Sequence(Sequence * seq);
#define Genomic_from_Sequence Wise2_Genomic_from_Sequence


/* Function:  show_Genomic(gen,ofp)
 *
 * Descrip:    For debugging
 *
 *
 * Arg:        gen [UNKN ] Undocumented argument [Genomic *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_show_Genomic(Genomic * gen,FILE * ofp);
#define show_Genomic Wise2_show_Genomic


/* Function:  hard_link_GenomicRepeat(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GenomicRepeat *]
 *
 * Return [UNKN ]  Undocumented return value [GenomicRepeat *]
 *
 */
GenomicRepeat * Wise2_hard_link_GenomicRepeat(GenomicRepeat * obj);
#define hard_link_GenomicRepeat Wise2_hard_link_GenomicRepeat


/* Function:  GenomicRepeat_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GenomicRepeat *]
 *
 */
GenomicRepeat * Wise2_GenomicRepeat_alloc(void);
#define GenomicRepeat_alloc Wise2_GenomicRepeat_alloc


/* Function:  free_GenomicRepeat(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GenomicRepeat *]
 *
 * Return [UNKN ]  Undocumented return value [GenomicRepeat *]
 *
 */
GenomicRepeat * Wise2_free_GenomicRepeat(GenomicRepeat * obj);
#define free_GenomicRepeat Wise2_free_GenomicRepeat


/* Function:  add_Genomic(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Genomic *]
 * Arg:        add [OWNER] Object to add to the list [GenomicRepeat *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_Genomic(Genomic * obj,GenomicRepeat * add);
#define add_Genomic Wise2_add_Genomic


/* Function:  flush_Genomic(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [Genomic *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_Genomic(Genomic * obj);
#define flush_Genomic Wise2_flush_Genomic


/* Function:  Genomic_alloc_std(void)
 *
 * Descrip:    Equivalent to Genomic_alloc_len(GenomicLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Genomic *]
 *
 */
Genomic * Wise2_Genomic_alloc_std(void);
#define Genomic_alloc_std Wise2_Genomic_alloc_std


/* Function:  Genomic_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [Genomic *]
 *
 */
Genomic * Wise2_Genomic_alloc_len(int len);
#define Genomic_alloc_len Wise2_Genomic_alloc_len


/* Function:  hard_link_Genomic(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Genomic *]
 *
 * Return [UNKN ]  Undocumented return value [Genomic *]
 *
 */
Genomic * Wise2_hard_link_Genomic(Genomic * obj);
#define hard_link_Genomic Wise2_hard_link_Genomic


/* Function:  Genomic_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Genomic *]
 *
 */
Genomic * Wise2_Genomic_alloc(void);
#define Genomic_alloc Wise2_Genomic_alloc


/* Function:  free_Genomic(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Genomic *]
 *
 * Return [UNKN ]  Undocumented return value [Genomic *]
 *
 */
Genomic * Wise2_free_Genomic(Genomic * obj);
#define free_Genomic Wise2_free_Genomic


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
char * Wise2_access_type_GenomicRepeat(GenomicRepeat * obj);
#define access_type_GenomicRepeat Wise2_access_type_GenomicRepeat
GenomicRepeat * Wise2_access_repeat_Genomic(Genomic * obj,int i);
#define access_repeat_Genomic Wise2_access_repeat_Genomic
int Wise2_length_repeat_Genomic(Genomic * obj);
#define length_repeat_Genomic Wise2_length_repeat_Genomic
boolean Wise2_replace_start_GenomicRepeat(GenomicRepeat * obj,int start);
#define replace_start_GenomicRepeat Wise2_replace_start_GenomicRepeat
boolean Wise2_replace_type_GenomicRepeat(GenomicRepeat * obj,char * type);
#define replace_type_GenomicRepeat Wise2_replace_type_GenomicRepeat
int Wise2_access_start_GenomicRepeat(GenomicRepeat * obj);
#define access_start_GenomicRepeat Wise2_access_start_GenomicRepeat
Sequence * Wise2_access_baseseq_Genomic(Genomic * obj);
#define access_baseseq_Genomic Wise2_access_baseseq_Genomic
boolean Wise2_replace_end_GenomicRepeat(GenomicRepeat * obj,int end);
#define replace_end_GenomicRepeat Wise2_replace_end_GenomicRepeat
boolean Wise2_replace_baseseq_Genomic(Genomic * obj,Sequence * baseseq);
#define replace_baseseq_Genomic Wise2_replace_baseseq_Genomic
int Wise2_access_end_GenomicRepeat(GenomicRepeat * obj);
#define access_end_GenomicRepeat Wise2_access_end_GenomicRepeat
void Wise2_swap_Genomic(GenomicRepeat ** list,int i,int j) ;
#define swap_Genomic Wise2_swap_Genomic
void Wise2_qsort_Genomic(GenomicRepeat ** list,int left,int right,int (*comp)(GenomicRepeat * ,GenomicRepeat * ));
#define qsort_Genomic Wise2_qsort_Genomic
void Wise2_sort_Genomic(Genomic * obj,int (*comp)(GenomicRepeat *, GenomicRepeat *));
#define sort_Genomic Wise2_sort_Genomic
boolean Wise2_expand_Genomic(Genomic * obj,int len);
#define expand_Genomic Wise2_expand_Genomic

#ifdef _cplusplus
}
#endif

#endif
