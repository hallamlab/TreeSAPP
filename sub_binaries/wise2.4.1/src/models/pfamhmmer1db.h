#ifndef DYNAMITEpfamhmmer1dbHEADERFILE
#define DYNAMITEpfamhmmer1dbHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "threestatemodel.h"
#define NO_HMMER_INCLUDES
#include "wise2xhmmer2.h"

#define PfamHmmer1DBLISTLENGTH 512

struct Wise2_PfamHmmer1Entry {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * entryname;    
    boolean is_random;   
    boolean is_hmmls;    
    double bits_cutoff;  
    } ;  
/* PfamHmmer1Entry defined */ 
#ifndef DYNAMITE_DEFINED_PfamHmmer1Entry
typedef struct Wise2_PfamHmmer1Entry Wise2_PfamHmmer1Entry;
#define PfamHmmer1Entry Wise2_PfamHmmer1Entry
#define DYNAMITE_DEFINED_PfamHmmer1Entry
#endif


/* Object PfamHmmer1DB
 *
 * Descrip: No Description
 *
 */
struct Wise2_PfamHmmer1DB {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    PfamHmmer1Entry ** en;   
    int len;/* len for above en  */ 
    int maxlen; /* maxlen for above en */ 
    char * dirname; /*  directory name with the models */ 
    int cur;     
    RandomModel * def;  /*  default random model */ 
    } ;  
/* PfamHmmer1DB defined */ 
#ifndef DYNAMITE_DEFINED_PfamHmmer1DB
typedef struct Wise2_PfamHmmer1DB Wise2_PfamHmmer1DB;
#define PfamHmmer1DB Wise2_PfamHmmer1DB
#define DYNAMITE_DEFINED_PfamHmmer1DB
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  ThreeStateModel_from_name_PfamHmmer1DB(phd,name)
 *
 * Descrip:    reads a named model - akin to indexing
 *
 *
 * Arg:         phd [UNKN ] Undocumented argument [PfamHmmer1DB *]
 * Arg:        name [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateModel *]
 *
 */
ThreeStateModel * Wise2_ThreeStateModel_from_name_PfamHmmer1DB(PfamHmmer1DB * phd,char * name);
#define ThreeStateModel_from_name_PfamHmmer1DB Wise2_ThreeStateModel_from_name_PfamHmmer1DB


/* Function:  read_next_TSM_PfamHmmer1DB(phd,return_status)
 *
 * Descrip:    reads the next threestatemodel
 *             for PfamHmmer1DB, placing the correct
 *             status into return_status( DB_RETURN_OK, etc).
 *
 *
 * Arg:                  phd [UNKN ] Undocumented argument [PfamHmmer1DB *]
 * Arg:        return_status [UNKN ] Undocumented argument [int *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateModel *]
 *
 */
ThreeStateModel * Wise2_read_next_TSM_PfamHmmer1DB(PfamHmmer1DB * phd,int * return_status);
#define read_next_TSM_PfamHmmer1DB Wise2_read_next_TSM_PfamHmmer1DB


/* Function:  read_TSM_from_PfamHmmer1Entry(en,dir)
 *
 * Descrip:    reads an individual HMMer model from the entry
 *             specification
 *
 *
 * Arg:         en [UNKN ] Undocumented argument [PfamHmmer1Entry *]
 * Arg:        dir [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateModel *]
 *
 */
ThreeStateModel * Wise2_read_TSM_from_PfamHmmer1Entry(PfamHmmer1Entry * en,char * dir);
#define read_TSM_from_PfamHmmer1Entry Wise2_read_TSM_from_PfamHmmer1Entry


/* Function:  PfamHmmer1DB_from_dirname(dirname)
 *
 * Descrip:    Makes a new PfamHmmer1DB from the dir name.
 *
 *             The directory should have a file HMMs which has
 *             entries like
 *
 *             rrm hmmls 12
 *
 *             with -r indicating a specialised random model.
 *
 *
 * Arg:        dirname [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [PfamHmmer1DB *]
 *
 */
PfamHmmer1DB * Wise2_PfamHmmer1DB_from_dirname(char * dirname);
#define PfamHmmer1DB_from_dirname Wise2_PfamHmmer1DB_from_dirname


/* Function:  hard_link_PfamHmmer1Entry(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [PfamHmmer1Entry *]
 *
 * Return [UNKN ]  Undocumented return value [PfamHmmer1Entry *]
 *
 */
PfamHmmer1Entry * Wise2_hard_link_PfamHmmer1Entry(PfamHmmer1Entry * obj);
#define hard_link_PfamHmmer1Entry Wise2_hard_link_PfamHmmer1Entry


/* Function:  PfamHmmer1Entry_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PfamHmmer1Entry *]
 *
 */
PfamHmmer1Entry * Wise2_PfamHmmer1Entry_alloc(void);
#define PfamHmmer1Entry_alloc Wise2_PfamHmmer1Entry_alloc


/* Function:  free_PfamHmmer1Entry(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [PfamHmmer1Entry *]
 *
 * Return [UNKN ]  Undocumented return value [PfamHmmer1Entry *]
 *
 */
PfamHmmer1Entry * Wise2_free_PfamHmmer1Entry(PfamHmmer1Entry * obj);
#define free_PfamHmmer1Entry Wise2_free_PfamHmmer1Entry


/* Function:  add_PfamHmmer1DB(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [PfamHmmer1DB *]
 * Arg:        add [OWNER] Object to add to the list [PfamHmmer1Entry *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_PfamHmmer1DB(PfamHmmer1DB * obj,PfamHmmer1Entry * add);
#define add_PfamHmmer1DB Wise2_add_PfamHmmer1DB


/* Function:  flush_PfamHmmer1DB(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [PfamHmmer1DB *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_PfamHmmer1DB(PfamHmmer1DB * obj);
#define flush_PfamHmmer1DB Wise2_flush_PfamHmmer1DB


/* Function:  PfamHmmer1DB_alloc_std(void)
 *
 * Descrip:    Equivalent to PfamHmmer1DB_alloc_len(PfamHmmer1DBLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PfamHmmer1DB *]
 *
 */
PfamHmmer1DB * Wise2_PfamHmmer1DB_alloc_std(void);
#define PfamHmmer1DB_alloc_std Wise2_PfamHmmer1DB_alloc_std


/* Function:  PfamHmmer1DB_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [PfamHmmer1DB *]
 *
 */
PfamHmmer1DB * Wise2_PfamHmmer1DB_alloc_len(int len);
#define PfamHmmer1DB_alloc_len Wise2_PfamHmmer1DB_alloc_len


/* Function:  hard_link_PfamHmmer1DB(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [PfamHmmer1DB *]
 *
 * Return [UNKN ]  Undocumented return value [PfamHmmer1DB *]
 *
 */
PfamHmmer1DB * Wise2_hard_link_PfamHmmer1DB(PfamHmmer1DB * obj);
#define hard_link_PfamHmmer1DB Wise2_hard_link_PfamHmmer1DB


/* Function:  PfamHmmer1DB_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PfamHmmer1DB *]
 *
 */
PfamHmmer1DB * Wise2_PfamHmmer1DB_alloc(void);
#define PfamHmmer1DB_alloc Wise2_PfamHmmer1DB_alloc


/* Function:  free_PfamHmmer1DB(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [PfamHmmer1DB *]
 *
 * Return [UNKN ]  Undocumented return value [PfamHmmer1DB *]
 *
 */
PfamHmmer1DB * Wise2_free_PfamHmmer1DB(PfamHmmer1DB * obj);
#define free_PfamHmmer1DB Wise2_free_PfamHmmer1DB


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
int Wise2_access_cur_PfamHmmer1DB(PfamHmmer1DB * obj);
#define access_cur_PfamHmmer1DB Wise2_access_cur_PfamHmmer1DB
boolean Wise2_replace_def_PfamHmmer1DB(PfamHmmer1DB * obj,RandomModel * def);
#define replace_def_PfamHmmer1DB Wise2_replace_def_PfamHmmer1DB
RandomModel * Wise2_access_def_PfamHmmer1DB(PfamHmmer1DB * obj);
#define access_def_PfamHmmer1DB Wise2_access_def_PfamHmmer1DB
boolean Wise2_replace_entryname_PfamHmmer1Entry(PfamHmmer1Entry * obj,char * entryname);
#define replace_entryname_PfamHmmer1Entry Wise2_replace_entryname_PfamHmmer1Entry
PfamHmmer1Entry * Wise2_access_en_PfamHmmer1DB(PfamHmmer1DB * obj,int i);
#define access_en_PfamHmmer1DB Wise2_access_en_PfamHmmer1DB
char * Wise2_access_entryname_PfamHmmer1Entry(PfamHmmer1Entry * obj);
#define access_entryname_PfamHmmer1Entry Wise2_access_entryname_PfamHmmer1Entry
boolean Wise2_replace_dirname_PfamHmmer1DB(PfamHmmer1DB * obj,char * dirname);
#define replace_dirname_PfamHmmer1DB Wise2_replace_dirname_PfamHmmer1DB
boolean Wise2_replace_is_random_PfamHmmer1Entry(PfamHmmer1Entry * obj,boolean is_random);
#define replace_is_random_PfamHmmer1Entry Wise2_replace_is_random_PfamHmmer1Entry
boolean Wise2_replace_cur_PfamHmmer1DB(PfamHmmer1DB * obj,int cur);
#define replace_cur_PfamHmmer1DB Wise2_replace_cur_PfamHmmer1DB
boolean Wise2_access_is_random_PfamHmmer1Entry(PfamHmmer1Entry * obj);
#define access_is_random_PfamHmmer1Entry Wise2_access_is_random_PfamHmmer1Entry
int Wise2_length_en_PfamHmmer1DB(PfamHmmer1DB * obj);
#define length_en_PfamHmmer1DB Wise2_length_en_PfamHmmer1DB
boolean Wise2_replace_is_hmmls_PfamHmmer1Entry(PfamHmmer1Entry * obj,boolean is_hmmls);
#define replace_is_hmmls_PfamHmmer1Entry Wise2_replace_is_hmmls_PfamHmmer1Entry
double Wise2_access_bits_cutoff_PfamHmmer1Entry(PfamHmmer1Entry * obj);
#define access_bits_cutoff_PfamHmmer1Entry Wise2_access_bits_cutoff_PfamHmmer1Entry
boolean Wise2_access_is_hmmls_PfamHmmer1Entry(PfamHmmer1Entry * obj);
#define access_is_hmmls_PfamHmmer1Entry Wise2_access_is_hmmls_PfamHmmer1Entry
char * Wise2_access_dirname_PfamHmmer1DB(PfamHmmer1DB * obj);
#define access_dirname_PfamHmmer1DB Wise2_access_dirname_PfamHmmer1DB
boolean Wise2_replace_bits_cutoff_PfamHmmer1Entry(PfamHmmer1Entry * obj,double bits_cutoff);
#define replace_bits_cutoff_PfamHmmer1Entry Wise2_replace_bits_cutoff_PfamHmmer1Entry
void Wise2_swap_PfamHmmer1DB(PfamHmmer1Entry ** list,int i,int j) ;
#define swap_PfamHmmer1DB Wise2_swap_PfamHmmer1DB
void Wise2_qsort_PfamHmmer1DB(PfamHmmer1Entry ** list,int left,int right,int (*comp)(PfamHmmer1Entry * ,PfamHmmer1Entry * ));
#define qsort_PfamHmmer1DB Wise2_qsort_PfamHmmer1DB
void Wise2_sort_PfamHmmer1DB(PfamHmmer1DB * obj,int (*comp)(PfamHmmer1Entry *, PfamHmmer1Entry *));
#define sort_PfamHmmer1DB Wise2_sort_PfamHmmer1DB
boolean Wise2_expand_PfamHmmer1DB(PfamHmmer1DB * obj,int len);
#define expand_PfamHmmer1DB Wise2_expand_PfamHmmer1DB

#ifdef _cplusplus
}
#endif

#endif
