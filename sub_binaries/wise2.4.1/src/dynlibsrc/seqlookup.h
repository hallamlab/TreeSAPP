#ifndef DYNAMITEseqlookupHEADERFILE
#define DYNAMITEseqlookupHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif

#include "sequence.h"
#include "sequencedb.h"

#define SEQLOOKUP_5AA_SIZE (26*26*26*26*26)

#define SEQLOOKUP_LOWCOMPLEXITY 1

#define ARRAYHEAD_IS_LOWCOMPLEXITY(head) (head->flags & SEQLOOKUP_LOWCOMPLEXITY)

typedef struct SeqLookupResultStruct {
  Sequence * seq;
  int pos;
} SeqLookupResultStruct;


typedef struct Wise2_ArraySeqLookup_Unit {
  Sequence * seq;
  int pos;
} ArraySeqLookupUnit;


typedef struct Wise2_ArraySeqHead {
  ArraySeqLookupUnit * units;
  int current_pos;
  int max;
  char flags;
} ArraySeqHead;

typedef enum SeqLookupLoadType {
  SeqLookupLoad_Protein = 12,
  SeqLookupLoad_DNA
} SeqLookupLoadType;

#define SeqLookupInterfaceLISTLENGTH 4028

/* Object SeqLookupResultInterface
 *
 * Descrip: This is the final interface returned on
 *        finding an occuypied index (see below about 
 *        how to get this). It basically represents
 *        an array of SeqLookupResultStruct. The
 *        interface can choose how to manage the
 *        memory - a client should call is_more function.
 *        If this function returns TRUE, then it can call
 *        next with the previous resultstruct passed in
 *        in the prev slot. For the first call it passes
 *        in a NULL pointer. The interface can decide whether
 *        to reuse the memory or not. Finally the client
 *        calls free_data to indicate that it does not want to
 *        use the information any more
 *
 *
 */
struct Wise2_SeqLookupResultInterface {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    SeqLookupResultStruct * (*next)(void * data,SeqLookupResultStruct * prev);   
    boolean (*is_more)(void *data);  
    void (*free_data)(void * data);  
    void * data;     
    } ;  
/* SeqLookupResultInterface defined */ 
#ifndef DYNAMITE_DEFINED_SeqLookupResultInterface
typedef struct Wise2_SeqLookupResultInterface Wise2_SeqLookupResultInterface;
#define SeqLookupResultInterface Wise2_SeqLookupResultInterface
#define DYNAMITE_DEFINED_SeqLookupResultInterface
#endif


/* Object SeqLookupClientInterface
 *
 * Descrip: This is a per-client interface got from the central SeqLookup
 *        interface. A client must guarentee that only a single thread
 *        will interact with a single SeqLookupClientInterface so
 *        interfaces have a chance to sensible manage their memory
 *
 *
 */
struct Wise2_SeqLookupClientInterface {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    SeqLookupResultInterface * (*lookup)(void * data,int seq_number);    
    boolean (*is_populated)(void * data,int seq_number); 
    void (*free_data)(void * data);  
    void * data;     
    } ;  
/* SeqLookupClientInterface defined */ 
#ifndef DYNAMITE_DEFINED_SeqLookupClientInterface
typedef struct Wise2_SeqLookupClientInterface Wise2_SeqLookupClientInterface;
#define SeqLookupClientInterface Wise2_SeqLookupClientInterface
#define DYNAMITE_DEFINED_SeqLookupClientInterface
#endif


struct Wise2_SeqLookupLoadPara {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    SeqLookupLoadType type;  
    int tile_freq;   
    int report_stagger;  
    int truncate;    
    int mark_low_complexity;     
    int start_seq_load;  
    int end_seq_load;    
    } ;  
/* SeqLookupLoadPara defined */ 
#ifndef DYNAMITE_DEFINED_SeqLookupLoadPara
typedef struct Wise2_SeqLookupLoadPara Wise2_SeqLookupLoadPara;
#define SeqLookupLoadPara Wise2_SeqLookupLoadPara
#define DYNAMITE_DEFINED_SeqLookupLoadPara
#endif


/* Object SeqLookupInterface
 *
 * Descrip: This interface defines the basic SeqLookup
 *        possibilities. Its main role is to give out
 *        clientInterfaces, which must be one per thread
 *        (in contrast a client can call get_client by
 *        multiple threads, though it should do each call via
 *        single threaded). 
 *
 *
 */
struct Wise2_SeqLookupInterface {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    SeqLookupClientInterface * (*get_client)(void * data);   
    boolean (*add_seq)(void * data,Sequence * seq,SeqLookupLoadPara * para); 
    ArraySeqHead * (*lookup_array_head)(void * data,int seq_number); 
    boolean (*add_direct_number)(void *data,int seq_number,Sequence * target,int pos);   
    void (*free_data)(void * data);  
    void * data;     
    Sequence ** seq_store;   
    int len;/* len for above seq_store  */ 
    int maxlen; /* maxlen for above seq_store */ 
    } ;  
/* SeqLookupInterface defined */ 
#ifndef DYNAMITE_DEFINED_SeqLookupInterface
typedef struct Wise2_SeqLookupInterface Wise2_SeqLookupInterface;
#define SeqLookupInterface Wise2_SeqLookupInterface
#define DYNAMITE_DEFINED_SeqLookupInterface
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  free_SeqLookupInterface(sli)
 *
 * Descrip:    Frees SeqLookupInterface - overrides dynamite default
 *
 *
 * Arg:        sli [UNKN ] Undocumented argument [SeqLookupInterface *]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupInterface *]
 *
 */
SeqLookupInterface * Wise2_free_SeqLookupInterface(SeqLookupInterface * sli);
#define free_SeqLookupInterface Wise2_free_SeqLookupInterface


/* Function:  load_SequenceDB_SeqLookupLoadPara(p,db,sli)
 *
 * Descrip:    Loads a SequenceDB into a hash on the basis of the SeqLookupLoadPara
 *
 *             returns the number of sequences loaded
 *
 *
 * Arg:          p [UNKN ] Undocumented argument [SeqLookupLoadPara *]
 * Arg:         db [UNKN ] Undocumented argument [SequenceDB *]
 * Arg:        sli [UNKN ] Undocumented argument [SeqLookupInterface *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_load_SequenceDB_SeqLookupLoadPara(SeqLookupLoadPara * p,SequenceDB * db,SeqLookupInterface * sli);
#define load_SequenceDB_SeqLookupLoadPara Wise2_load_SequenceDB_SeqLookupLoadPara


/* Function:  show_help_SeqLookupLoadPara(ofp)
 *
 * Descrip:    Shows help associated with Sequence loading
 *
 *
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_show_help_SeqLookupLoadPara(FILE * ofp);
#define show_help_SeqLookupLoadPara Wise2_show_help_SeqLookupLoadPara


/* Function:  new_SeqLookupLoadPara_from_argv(argc,argv)
 *
 * Descrip:    Builds new SeqLookup load from a command line
 *
 *
 * Arg:        argc [UNKN ] Undocumented argument [int *]
 * Arg:        argv [UNKN ] Undocumented argument [char **]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupLoadPara *]
 *
 */
SeqLookupLoadPara * Wise2_new_SeqLookupLoadPara_from_argv(int * argc,char ** argv);
#define new_SeqLookupLoadPara_from_argv Wise2_new_SeqLookupLoadPara_from_argv


/* Function:  seq_number_dna_15mer_noN(seq)
 *
 * Descrip:    Function for DNA sequence to number on 15mers,
 *             Ns get mapped to -1
 *
 *
 * Arg:        seq [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_seq_number_dna_15mer_noN(char * seq);
#define seq_number_dna_15mer_noN Wise2_seq_number_dna_15mer_noN


/* Function:  seq_number_dna_7mer_noN(seq)
 *
 * Descrip:    Function for DNA sequence to number on 15mers,
 *             Ns get mapped to -1
 *
 *
 * Arg:        seq [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_seq_number_dna_7mer_noN(char * seq);
#define seq_number_dna_7mer_noN Wise2_seq_number_dna_7mer_noN


/* Function:  seq_number_aa_5mer(seq)
 *
 * Descrip:    Function for the amino acid to number on 5mers
 *
 *
 * Arg:        seq [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_seq_number_aa_5mer(char * seq);
#define seq_number_aa_5mer Wise2_seq_number_aa_5mer


/* Function:  flags_from_5aa_sequence(seq)
 *
 * Descrip:    returns simple lowcomplexity flag or 
 *             not for this sequence
 *
 *
 * Arg:        seq [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [char]
 *
 */
char Wise2_flags_from_5aa_sequence(char * seq);
#define flags_from_5aa_sequence Wise2_flags_from_5aa_sequence


/* Function:  free_SeqLookupClientInterface(sli)
 *
 * Descrip:    Frees SeqLookupClientInterface - overrides dynamite default
 *
 *
 * Arg:        sli [UNKN ] Undocumented argument [SeqLookupClientInterface *]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupClientInterface *]
 *
 */
SeqLookupClientInterface * Wise2_free_SeqLookupClientInterface(SeqLookupClientInterface * sli);
#define free_SeqLookupClientInterface Wise2_free_SeqLookupClientInterface


/* Function:  free_SeqLookupResultInterface(sli)
 *
 * Descrip:    Frees SeqLookupResultInterface - overrides dynamite default
 *
 *
 * Arg:        sli [UNKN ] Undocumented argument [SeqLookupResultInterface *]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupResultInterface *]
 *
 */
SeqLookupResultInterface * Wise2_free_SeqLookupResultInterface(SeqLookupResultInterface * sli);
#define free_SeqLookupResultInterface Wise2_free_SeqLookupResultInterface


/* Function:  hard_link_SeqLookupResultInterface(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SeqLookupResultInterface *]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupResultInterface *]
 *
 */
SeqLookupResultInterface * Wise2_hard_link_SeqLookupResultInterface(SeqLookupResultInterface * obj);
#define hard_link_SeqLookupResultInterface Wise2_hard_link_SeqLookupResultInterface


/* Function:  SeqLookupResultInterface_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupResultInterface *]
 *
 */
SeqLookupResultInterface * Wise2_SeqLookupResultInterface_alloc(void);
#define SeqLookupResultInterface_alloc Wise2_SeqLookupResultInterface_alloc


/* Function:  hard_link_SeqLookupClientInterface(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SeqLookupClientInterface *]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupClientInterface *]
 *
 */
SeqLookupClientInterface * Wise2_hard_link_SeqLookupClientInterface(SeqLookupClientInterface * obj);
#define hard_link_SeqLookupClientInterface Wise2_hard_link_SeqLookupClientInterface


/* Function:  SeqLookupClientInterface_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupClientInterface *]
 *
 */
SeqLookupClientInterface * Wise2_SeqLookupClientInterface_alloc(void);
#define SeqLookupClientInterface_alloc Wise2_SeqLookupClientInterface_alloc


/* Function:  hard_link_SeqLookupLoadPara(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SeqLookupLoadPara *]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupLoadPara *]
 *
 */
SeqLookupLoadPara * Wise2_hard_link_SeqLookupLoadPara(SeqLookupLoadPara * obj);
#define hard_link_SeqLookupLoadPara Wise2_hard_link_SeqLookupLoadPara


/* Function:  SeqLookupLoadPara_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupLoadPara *]
 *
 */
SeqLookupLoadPara * Wise2_SeqLookupLoadPara_alloc(void);
#define SeqLookupLoadPara_alloc Wise2_SeqLookupLoadPara_alloc


/* Function:  free_SeqLookupLoadPara(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SeqLookupLoadPara *]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupLoadPara *]
 *
 */
SeqLookupLoadPara * Wise2_free_SeqLookupLoadPara(SeqLookupLoadPara * obj);
#define free_SeqLookupLoadPara Wise2_free_SeqLookupLoadPara


/* Function:  add_SeqLookupInterface(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [SeqLookupInterface *]
 * Arg:        add [OWNER] Object to add to the list [Sequence *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_SeqLookupInterface(SeqLookupInterface * obj,Sequence * add);
#define add_SeqLookupInterface Wise2_add_SeqLookupInterface


/* Function:  flush_SeqLookupInterface(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [SeqLookupInterface *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_SeqLookupInterface(SeqLookupInterface * obj);
#define flush_SeqLookupInterface Wise2_flush_SeqLookupInterface


/* Function:  SeqLookupInterface_alloc_std(void)
 *
 * Descrip:    Equivalent to SeqLookupInterface_alloc_len(SeqLookupInterfaceLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupInterface *]
 *
 */
SeqLookupInterface * Wise2_SeqLookupInterface_alloc_std(void);
#define SeqLookupInterface_alloc_std Wise2_SeqLookupInterface_alloc_std


/* Function:  SeqLookupInterface_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupInterface *]
 *
 */
SeqLookupInterface * Wise2_SeqLookupInterface_alloc_len(int len);
#define SeqLookupInterface_alloc_len Wise2_SeqLookupInterface_alloc_len


/* Function:  hard_link_SeqLookupInterface(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SeqLookupInterface *]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupInterface *]
 *
 */
SeqLookupInterface * Wise2_hard_link_SeqLookupInterface(SeqLookupInterface * obj);
#define hard_link_SeqLookupInterface Wise2_hard_link_SeqLookupInterface


/* Function:  SeqLookupInterface_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupInterface *]
 *
 */
SeqLookupInterface * Wise2_SeqLookupInterface_alloc(void);
#define SeqLookupInterface_alloc Wise2_SeqLookupInterface_alloc


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_SeqLookupInterface(Sequence ** list,int i,int j) ;
#define swap_SeqLookupInterface Wise2_swap_SeqLookupInterface
void Wise2_qsort_SeqLookupInterface(Sequence ** list,int left,int right,int (*comp)(Sequence * ,Sequence * ));
#define qsort_SeqLookupInterface Wise2_qsort_SeqLookupInterface
void Wise2_sort_SeqLookupInterface(SeqLookupInterface * obj,int (*comp)(Sequence *, Sequence *));
#define sort_SeqLookupInterface Wise2_sort_SeqLookupInterface
boolean Wise2_expand_SeqLookupInterface(SeqLookupInterface * obj,int len);
#define expand_SeqLookupInterface Wise2_expand_SeqLookupInterface

#ifdef _cplusplus
}
#endif

#endif
