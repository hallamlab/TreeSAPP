#ifndef DYNAMITEseqalignHEADERFILE
#define DYNAMITEseqalignHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "sequence.h"

#define SeqAlignLISTLENGTH 64

/* Object SeqAlign
 *
 * Descrip: The most stupid and bad representation
 *        of a sequence alignment, being a list
 *        of sequences with padding characters
 *        in them
 *
 *        Not very useful for anything but reformatting
 *        and column frequencies. Dont use this
 *        for complicated alignment manipulations. Use
 *        the more hard core aln stuff
 *
 *        For reading in this data structure, you can use
 *        the wise2xhmmer2x bridge to Sean Eddy's reading
 *        code. This copes inherently with automatically
 *        detecting sequence alignment.
 *
 *        Might bridge to his output stuff sometime as well...
 *
 *
 */
struct Wise2_SeqAlign {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * name;     
    Sequence ** seq;     
    int len;/* len for above seq  */ 
    int maxlen; /* maxlen for above seq */ 
    } ;  
/* SeqAlign defined */ 
#ifndef DYNAMITE_DEFINED_SeqAlign
typedef struct Wise2_SeqAlign Wise2_SeqAlign;
#define SeqAlign Wise2_SeqAlign
#define DYNAMITE_DEFINED_SeqAlign
#endif


/* Object ColumnCount
 *
 * Descrip: counts from a column of
 *        an alignment
 *
 *
 */
struct Wise2_ColumnCount {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    double count[27];    
    } ;  
/* ColumnCount defined */ 
#ifndef DYNAMITE_DEFINED_ColumnCount
typedef struct Wise2_ColumnCount Wise2_ColumnCount;
#define ColumnCount Wise2_ColumnCount
#define DYNAMITE_DEFINED_ColumnCount
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  is_gapped_SeqAlign(sa,pos)
 *
 * Descrip:    Tells you whether either . or - characters occur at
 *             a position. Uses C-style coordinates
 *
 *
 * Arg:         sa [UNKN ] Undocumented argument [SeqAlign *]
 * Arg:        pos [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_is_gapped_SeqAlign(SeqAlign * sa,int pos);
#define is_gapped_SeqAlign Wise2_is_gapped_SeqAlign


/* Function:  trim_from_N_SeqAlign(sa)
 *
 * Descrip:    trims sequence alignments to first and last
 *             non N character (for DNA obviously)
 *
 *
 * Arg:        sa [UNKN ] Undocumented argument [SeqAlign *]
 *
 * Return [UNKN ]  Undocumented return value [SeqAlign *]
 *
 */
SeqAlign * Wise2_trim_from_N_SeqAlign(SeqAlign * sa);
#define trim_from_N_SeqAlign Wise2_trim_from_N_SeqAlign


/* Function:  reweight_SeqAlign(sal)
 *
 * Descrip:    Provides reweighting for a seqalign so all the weights
 *             add up to 1.0
 *
 *
 * Arg:        sal [UNKN ] Undocumented argument [SeqAlign *]
 *
 */
void Wise2_reweight_SeqAlign(SeqAlign * sal);
#define reweight_SeqAlign Wise2_reweight_SeqAlign


/* Function:  ColumnCount_from_SeqAlign(sa,col)
 *
 * Descrip:    Gives you a column count
 *
 *             You are supposed to have enough sense to
 *             cache this on the caller if you need it again
 *
 *
 * Arg:         sa [UNKN ] Undocumented argument [SeqAlign *]
 * Arg:        col [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [ColumnCount *]
 *
 */
ColumnCount * Wise2_ColumnCount_from_SeqAlign(SeqAlign * sa,int col);
#define ColumnCount_from_SeqAlign Wise2_ColumnCount_from_SeqAlign


/* Function:  write_selex_SeqAlign(sa,name_len,block_len,ofp)
 *
 * Descrip:    Writes selex format
 *
 *
 * Arg:               sa [UNKN ] Undocumented argument [const SeqAlign *]
 * Arg:         name_len [UNKN ] Undocumented argument [int]
 * Arg:        block_len [UNKN ] Undocumented argument [int]
 * Arg:              ofp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_write_selex_SeqAlign(const SeqAlign * sa,int name_len,int block_len,FILE * ofp);
#define write_selex_SeqAlign Wise2_write_selex_SeqAlign


/* Function:  read_fasta_SeqAlign_file(filename)
 *
 * Descrip:    Reads in fasta file opening the file
 *
 *
 * Arg:        filename [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [SeqAlign *]
 *
 */
SeqAlign * Wise2_read_fasta_SeqAlign_file(char * filename);
#define read_fasta_SeqAlign_file Wise2_read_fasta_SeqAlign_file


/* Function:  read_fasta_SeqAlign(ifp)
 *
 * Descrip:    Reads in fasta file style alignments
 *
 *
 * Arg:        ifp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [SeqAlign *]
 *
 */
SeqAlign * Wise2_read_fasta_SeqAlign(FILE * ifp);
#define read_fasta_SeqAlign Wise2_read_fasta_SeqAlign


/* Function:  write_fasta_SeqAlign(sa,ofp)
 *
 * Descrip:    writes out Fasta file
 *
 *
 * Arg:         sa [UNKN ] Undocumented argument [SeqAlign *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_write_fasta_SeqAlign(SeqAlign * sa,FILE * ofp);
#define write_fasta_SeqAlign Wise2_write_fasta_SeqAlign


/* Function:  read_selex_SeqAlign_file(filename)
 *
 * Descrip:    Open files as well
 *
 *
 * Arg:        filename [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [SeqAlign *]
 *
 */
SeqAlign * Wise2_read_selex_SeqAlign_file(char * filename);
#define read_selex_SeqAlign_file Wise2_read_selex_SeqAlign_file


/* Function:  read_selex_SeqAlign(ifp)
 *
 * Descrip:    Reads in selex (Stockholm) alignment
 *
 *             At the moment ignores all #= stuff on the
 *             sequence
 *
 *             Read HMMER documentation for a definition of
 *             the format
 *
 *
 * Arg:        ifp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [SeqAlign *]
 *
 */
SeqAlign * Wise2_read_selex_SeqAlign(FILE * ifp);
#define read_selex_SeqAlign Wise2_read_selex_SeqAlign


/* Function:  add_SeqAlign(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [SeqAlign *]
 * Arg:        add [OWNER] Object to add to the list [Sequence *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_SeqAlign(SeqAlign * obj,Sequence * add);
#define add_SeqAlign Wise2_add_SeqAlign


/* Function:  flush_SeqAlign(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [SeqAlign *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_SeqAlign(SeqAlign * obj);
#define flush_SeqAlign Wise2_flush_SeqAlign


/* Function:  SeqAlign_alloc_std(void)
 *
 * Descrip:    Equivalent to SeqAlign_alloc_len(SeqAlignLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SeqAlign *]
 *
 */
SeqAlign * Wise2_SeqAlign_alloc_std(void);
#define SeqAlign_alloc_std Wise2_SeqAlign_alloc_std


/* Function:  SeqAlign_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [SeqAlign *]
 *
 */
SeqAlign * Wise2_SeqAlign_alloc_len(int len);
#define SeqAlign_alloc_len Wise2_SeqAlign_alloc_len


/* Function:  hard_link_SeqAlign(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SeqAlign *]
 *
 * Return [UNKN ]  Undocumented return value [SeqAlign *]
 *
 */
SeqAlign * Wise2_hard_link_SeqAlign(SeqAlign * obj);
#define hard_link_SeqAlign Wise2_hard_link_SeqAlign


/* Function:  SeqAlign_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SeqAlign *]
 *
 */
SeqAlign * Wise2_SeqAlign_alloc(void);
#define SeqAlign_alloc Wise2_SeqAlign_alloc


/* Function:  free_SeqAlign(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SeqAlign *]
 *
 * Return [UNKN ]  Undocumented return value [SeqAlign *]
 *
 */
SeqAlign * Wise2_free_SeqAlign(SeqAlign * obj);
#define free_SeqAlign Wise2_free_SeqAlign


/* Function:  hard_link_ColumnCount(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ColumnCount *]
 *
 * Return [UNKN ]  Undocumented return value [ColumnCount *]
 *
 */
ColumnCount * Wise2_hard_link_ColumnCount(ColumnCount * obj);
#define hard_link_ColumnCount Wise2_hard_link_ColumnCount


/* Function:  ColumnCount_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ColumnCount *]
 *
 */
ColumnCount * Wise2_ColumnCount_alloc(void);
#define ColumnCount_alloc Wise2_ColumnCount_alloc


/* Function:  free_ColumnCount(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ColumnCount *]
 *
 * Return [UNKN ]  Undocumented return value [ColumnCount *]
 *
 */
ColumnCount * Wise2_free_ColumnCount(ColumnCount * obj);
#define free_ColumnCount Wise2_free_ColumnCount


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_SeqAlign(Sequence ** list,int i,int j) ;
#define swap_SeqAlign Wise2_swap_SeqAlign
void Wise2_qsort_SeqAlign(Sequence ** list,int left,int right,int (*comp)(Sequence * ,Sequence * ));
#define qsort_SeqAlign Wise2_qsort_SeqAlign
void Wise2_sort_SeqAlign(SeqAlign * obj,int (*comp)(Sequence *, Sequence *));
#define sort_SeqAlign Wise2_sort_SeqAlign
boolean Wise2_expand_SeqAlign(SeqAlign * obj,int len);
#define expand_SeqAlign Wise2_expand_SeqAlign

#ifdef _cplusplus
}
#endif

#endif
