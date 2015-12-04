#ifndef DYNAMITEgeneHEADERFILE
#define DYNAMITEgeneHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"

#define GeneLISTLENGTH 16

#ifndef DYNAMITE_DEFINED_GenomicRegion
typedef struct Wise2_GenomicRegion Wise2_GenomicRegion;
#define GenomicRegion Wise2_GenomicRegion
#define DYNAMITE_DEFINED_GenomicRegion
#endif

#ifndef DYNAMITE_DEFINED_Transcript
typedef struct Wise2_Transcript Wise2_Transcript;
#define Transcript Wise2_Transcript
#define DYNAMITE_DEFINED_Transcript
#endif

/* Object Gene
 *
 * Descrip: Gene is the datastructure which represents a single
 *        gene. At the moment this is considered to be a series
 *        of transcripts (the first transcript being "canonical")
 *        which are made from a certain start/stop region in
 *        genomic DNA.
 *
 *        The gene datastructure does not necessarily contain
 *        any DNA sequence. When someone askes for the gene sequence,
 *        via get_Genomic_from_Gene(), it first sees if there
 *        is anything in the Genomic * 'cache'. If this is null,
 *        it then looks at parent (the genomic region), and if
 *        that is null, complains and returns null. Otherwise it
 *        truncates its parent's dna in the correct way, places
 *        the data structure into the genomic * cache, and returns
 *        it.
 *
 *        The name, bits and seqname have put into this datastructure
 *        for convience of carrying around this information into some
 *        of the gene prediction output formats. Probabaly
 *          o they should be in transcript, not gene
 *          o they shouldn't be here at all.
 *
 *        <sigh>
 *
 *
 */
struct Wise2_Gene {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int start;   
    int end;     
    GenomicRegion * parent; /*  may not be here */ 
    Genomic * genomic;  /*  may not be here! */ 
    Transcript ** transcript;    
    int len;/* len for above transcript  */ 
    int maxlen; /* maxlen for above transcript */ 
    char * name;    /*  ugly . Need a better system */ 
    double bits;    /*  ditto... */ 
    char * seqname; /*  very bad! this is for keeping track of what sequence was used to make the gene */ 
    boolean ispseudo;   /*  is a pseudogene or not */ 
    } ;  
/* Gene defined */ 
#ifndef DYNAMITE_DEFINED_Gene
typedef struct Wise2_Gene Wise2_Gene;
#define Gene Wise2_Gene
#define DYNAMITE_DEFINED_Gene
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  reversed_Gene(g)
 *
 * Descrip:    is this gene reversed?
 *
 *
 * Arg:        g [UNKN ] Undocumented argument [Gene *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_reversed_Gene(Gene * g);
#define reversed_Gene Wise2_reversed_Gene


/* Function:  copy_Gene(g)
 *
 * Descrip:    Makes a completely fresh copy of a
 *             gene
 *
 *
 * Arg:        g [UNKN ] Undocumented argument [Gene *]
 *
 * Return [UNKN ]  Undocumented return value [Gene *]
 *
 */
Gene * Wise2_copy_Gene(Gene * g);
#define copy_Gene Wise2_copy_Gene


/* Function:  is_simple_prediction_Gene(g)
 *
 * Descrip:    Does this gene have 
 *             	a single transcript
 *             	that transcript with translation start/end 
 *             	at the ends
 *
 *
 * Arg:        g [UNKN ] Undocumented argument [Gene *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_is_simple_prediction_Gene(Gene * g);
#define is_simple_prediction_Gene Wise2_is_simple_prediction_Gene


/* Function:  get_Genomic_from_Gene(gene)
 *
 * Descrip:    Gives back a Genomic sequence type
 *             from a gene.
 *
 *
 * Arg:        gene [READ ] gene to get Genomic from [Gene *]
 *
 * Return [SOFT ]  Genomic DNA data structure [Genomic *]
 *
 */
Genomic * Wise2_get_Genomic_from_Gene(Gene * gene);
#define get_Genomic_from_Gene Wise2_get_Genomic_from_Gene


/* Function:  read_EMBL_feature_Gene(buffer,maxlen,ifp)
 *
 * Descrip:    Reads in an EMBL feature table.
 *
 *             It expects to be passed a buffer with 'FT   CDS'.
 *             or 'FT   mRNA' in it. It will then 
 *             use the buffer to read successive lines of the Feature table
 *             until it comes to the next 'meta' feature line (ie, 3 place point).
 *
 *             It will use functions in /embl module for the reading.
 *
 *
 * Arg:        buffer [UNKN ] a string with FT  CDS line in it [char *]
 * Arg:        maxlen [UNKN ] length of the buffer [int]
 * Arg:           ifp [UNKN ] file stream with the rest of the feature table in it [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [Gene *]
 *
 */
Gene * Wise2_read_EMBL_feature_Gene(char * buffer,int maxlen,FILE * ifp);
#define read_EMBL_feature_Gene Wise2_read_EMBL_feature_Gene


/* Function:  write_Embl_FT_Gene(ge,key,ofp)
 *
 * Descrip:    shows a embl feature table part
 *
 *
 * Arg:         ge [UNKN ] Undocumented argument [Gene *]
 * Arg:        key [UNKN ] Undocumented argument [char *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_write_Embl_FT_Gene(Gene * ge,char * key,FILE * ofp);
#define write_Embl_FT_Gene Wise2_write_Embl_FT_Gene


/* Function:  show_pretty_Gene(ge,show_supporting,ofp)
 *
 * Descrip:    Shows a gene in the biologically accepted form
 *
 *
 * Arg:                     ge [UNKN ] Undocumented argument [Gene *]
 * Arg:        show_supporting [UNKN ] Undocumented argument [boolean]
 * Arg:                    ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_show_pretty_Gene(Gene * ge,boolean show_supporting,FILE * ofp);
#define show_pretty_Gene Wise2_show_pretty_Gene


/* Function:  show_Gene(ge,ofp)
 *
 * Descrip:    shows a gene in a vaguely human readable form
 *
 *
 * Arg:         ge [UNKN ] Undocumented argument [Gene *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_show_Gene(Gene * ge,FILE * ofp);
#define show_Gene Wise2_show_Gene


/* Function:  add_Gene(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Gene *]
 * Arg:        add [OWNER] Object to add to the list [Transcript *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_Gene(Gene * obj,Transcript * add);
#define add_Gene Wise2_add_Gene


/* Function:  flush_Gene(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [Gene *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_Gene(Gene * obj);
#define flush_Gene Wise2_flush_Gene


/* Function:  Gene_alloc_std(void)
 *
 * Descrip:    Equivalent to Gene_alloc_len(GeneLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Gene *]
 *
 */
Gene * Wise2_Gene_alloc_std(void);
#define Gene_alloc_std Wise2_Gene_alloc_std


/* Function:  Gene_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [Gene *]
 *
 */
Gene * Wise2_Gene_alloc_len(int len);
#define Gene_alloc_len Wise2_Gene_alloc_len


/* Function:  hard_link_Gene(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Gene *]
 *
 * Return [UNKN ]  Undocumented return value [Gene *]
 *
 */
Gene * Wise2_hard_link_Gene(Gene * obj);
#define hard_link_Gene Wise2_hard_link_Gene


/* Function:  Gene_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Gene *]
 *
 */
Gene * Wise2_Gene_alloc(void);
#define Gene_alloc Wise2_Gene_alloc


/* Function:  free_Gene(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Gene *]
 *
 * Return [UNKN ]  Undocumented return value [Gene *]
 *
 */
Gene * Wise2_free_Gene(Gene * obj);
#define free_Gene Wise2_free_Gene


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
boolean Wise2_replace_genomic_Gene(Gene * obj,Genomic * genomic);
#define replace_genomic_Gene Wise2_replace_genomic_Gene
Genomic * Wise2_access_genomic_Gene(Gene * obj);
#define access_genomic_Gene Wise2_access_genomic_Gene
Transcript * Wise2_access_transcript_Gene(Gene * obj,int i);
#define access_transcript_Gene Wise2_access_transcript_Gene
int Wise2_length_transcript_Gene(Gene * obj);
#define length_transcript_Gene Wise2_length_transcript_Gene
boolean Wise2_replace_start_Gene(Gene * obj,int start);
#define replace_start_Gene Wise2_replace_start_Gene
boolean Wise2_replace_name_Gene(Gene * obj,char * name);
#define replace_name_Gene Wise2_replace_name_Gene
boolean Wise2_replace_end_Gene(Gene * obj,int end);
#define replace_end_Gene Wise2_replace_end_Gene
char * Wise2_access_name_Gene(Gene * obj);
#define access_name_Gene Wise2_access_name_Gene
boolean Wise2_replace_parent_Gene(Gene * obj,GenomicRegion * parent);
#define replace_parent_Gene Wise2_replace_parent_Gene
boolean Wise2_replace_bits_Gene(Gene * obj,double bits);
#define replace_bits_Gene Wise2_replace_bits_Gene
boolean Wise2_access_ispseudo_Gene(Gene * obj);
#define access_ispseudo_Gene Wise2_access_ispseudo_Gene
double Wise2_access_bits_Gene(Gene * obj);
#define access_bits_Gene Wise2_access_bits_Gene
int Wise2_access_end_Gene(Gene * obj);
#define access_end_Gene Wise2_access_end_Gene
boolean Wise2_replace_seqname_Gene(Gene * obj,char * seqname);
#define replace_seqname_Gene Wise2_replace_seqname_Gene
int Wise2_access_start_Gene(Gene * obj);
#define access_start_Gene Wise2_access_start_Gene
char * Wise2_access_seqname_Gene(Gene * obj);
#define access_seqname_Gene Wise2_access_seqname_Gene
GenomicRegion * Wise2_access_parent_Gene(Gene * obj);
#define access_parent_Gene Wise2_access_parent_Gene
boolean Wise2_replace_ispseudo_Gene(Gene * obj,boolean ispseudo);
#define replace_ispseudo_Gene Wise2_replace_ispseudo_Gene
void Wise2_swap_Gene(Transcript ** list,int i,int j) ;
#define swap_Gene Wise2_swap_Gene
void Wise2_qsort_Gene(Transcript ** list,int left,int right,int (*comp)(Transcript * ,Transcript * ));
#define qsort_Gene Wise2_qsort_Gene
void Wise2_sort_Gene(Gene * obj,int (*comp)(Transcript *, Transcript *));
#define sort_Gene Wise2_sort_Gene
boolean Wise2_expand_Gene(Gene * obj,int len);
#define expand_Gene Wise2_expand_Gene

#ifdef _cplusplus
}
#endif

#endif
