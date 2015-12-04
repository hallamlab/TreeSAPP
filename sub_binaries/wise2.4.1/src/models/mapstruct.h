#ifndef DYNAMITEmapstructHEADERFILE
#define DYNAMITEmapstructHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"

#define MappedCloneSetLISTLENGTH 1024

struct Wise2_MappedClone {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * clone_name;   
    char * accession;    
    char * contig;   
    int start;   
    int end;     
    int seen;    
    } ;  
/* MappedClone defined */ 
#ifndef DYNAMITE_DEFINED_MappedClone
typedef struct Wise2_MappedClone Wise2_MappedClone;
#define MappedClone Wise2_MappedClone
#define DYNAMITE_DEFINED_MappedClone
#endif


struct Wise2_MappedCloneSet {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    MappedClone ** clone;    
    int len;/* len for above clone  */ 
    int maxlen; /* maxlen for above clone */ 
    int length;  
    int cursor;  
    } ;  
/* MappedCloneSet defined */ 
#ifndef DYNAMITE_DEFINED_MappedCloneSet
typedef struct Wise2_MappedCloneSet Wise2_MappedCloneSet;
#define MappedCloneSet Wise2_MappedCloneSet
#define DYNAMITE_DEFINED_MappedCloneSet
#endif


struct Wise2_MappedCloneMatch {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int ** matrix;  /*  NB i,j proper */ 
    int leni;   /* leni for above matrix  */ 
    int maxleni;/* max length for above pointer set */ 
    int lenj;   /* lenj for above matrix  */ 
    int maxlenj;/* max length for above pointer set */ 
    int * skip_iset;     
    int * skip_jset;     
    } ;  
/* MappedCloneMatch defined */ 
#ifndef DYNAMITE_DEFINED_MappedCloneMatch
typedef struct Wise2_MappedCloneMatch Wise2_MappedCloneMatch;
#define MappedCloneMatch Wise2_MappedCloneMatch
#define DYNAMITE_DEFINED_MappedCloneMatch
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  synchronise_MappedCloneSets(one,two)
 *
 * Descrip:    updates the internal seen flags for the clone sets in
 *             preparation for the dp
 *
 *
 * Arg:        one [UNKN ] Undocumented argument [MappedCloneSet *]
 * Arg:        two [UNKN ] Undocumented argument [MappedCloneSet *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_synchronise_MappedCloneSets(MappedCloneSet * one,MappedCloneSet * two);
#define synchronise_MappedCloneSets Wise2_synchronise_MappedCloneSets


/* Function:  subsection_MappedCloneSet(mcs,coord_start,coord_end,only_seen)
 *
 * Descrip:    Returns a sub-section of the MappedClone 
 *
 *
 * Arg:                mcs [UNKN ] Undocumented argument [MappedCloneSet *]
 * Arg:        coord_start [UNKN ] Undocumented argument [int]
 * Arg:          coord_end [UNKN ] Undocumented argument [int]
 * Arg:          only_seen [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [MappedCloneSet *]
 *
 */
MappedCloneSet * Wise2_subsection_MappedCloneSet(MappedCloneSet * mcs,int coord_start,int coord_end,int only_seen);
#define subsection_MappedCloneSet Wise2_subsection_MappedCloneSet


/* Function:  find_named_MappedClone(mcs,clone_name)
 *
 * Descrip:    Finds a mapped clone set with this name
 *
 *
 * Arg:               mcs [UNKN ] Undocumented argument [MappedCloneSet *]
 * Arg:        clone_name [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [MappedClone *]
 *
 */
MappedClone * Wise2_find_named_MappedClone(MappedCloneSet * mcs,char * clone_name);
#define find_named_MappedClone Wise2_find_named_MappedClone


/* Function:  start_comp_MappedClone(a,b)
 *
 * Descrip:    sorting for MappedClones
 *
 *
 * Arg:        a [UNKN ] Undocumented argument [MappedClone *]
 * Arg:        b [UNKN ] Undocumented argument [MappedClone *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_start_comp_MappedClone(MappedClone * a,MappedClone * b);
#define start_comp_MappedClone Wise2_start_comp_MappedClone


/* Function:  read_MappedCloneSet(ifp)
 *
 * Descrip:    Reads in a MappedCloneSet file
 *
 *
 * Arg:        ifp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [MappedCloneSet *]
 *
 */
MappedCloneSet * Wise2_read_MappedCloneSet(FILE * ifp) ;
#define read_MappedCloneSet Wise2_read_MappedCloneSet


/* Function:  read_MappedClone_line(line)
 *
 * Descrip:    Provides a mapped clone from a name\tstart\tend format
 *
 *
 * Arg:        line [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [MappedClone *]
 *
 */
MappedClone * Wise2_read_MappedClone_line(char * line);
#define read_MappedClone_line Wise2_read_MappedClone_line


/* Function:  hard_link_MappedClone(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [MappedClone *]
 *
 * Return [UNKN ]  Undocumented return value [MappedClone *]
 *
 */
MappedClone * Wise2_hard_link_MappedClone(MappedClone * obj);
#define hard_link_MappedClone Wise2_hard_link_MappedClone


/* Function:  MappedClone_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [MappedClone *]
 *
 */
MappedClone * Wise2_MappedClone_alloc(void);
#define MappedClone_alloc Wise2_MappedClone_alloc


/* Function:  free_MappedClone(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [MappedClone *]
 *
 * Return [UNKN ]  Undocumented return value [MappedClone *]
 *
 */
MappedClone * Wise2_free_MappedClone(MappedClone * obj);
#define free_MappedClone Wise2_free_MappedClone


/* Function:  add_MappedCloneSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [MappedCloneSet *]
 * Arg:        add [OWNER] Object to add to the list [MappedClone *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_MappedCloneSet(MappedCloneSet * obj,MappedClone * add);
#define add_MappedCloneSet Wise2_add_MappedCloneSet


/* Function:  flush_MappedCloneSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [MappedCloneSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_MappedCloneSet(MappedCloneSet * obj);
#define flush_MappedCloneSet Wise2_flush_MappedCloneSet


/* Function:  MappedCloneSet_alloc_std(void)
 *
 * Descrip:    Equivalent to MappedCloneSet_alloc_len(MappedCloneSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [MappedCloneSet *]
 *
 */
MappedCloneSet * Wise2_MappedCloneSet_alloc_std(void);
#define MappedCloneSet_alloc_std Wise2_MappedCloneSet_alloc_std


/* Function:  MappedCloneSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [MappedCloneSet *]
 *
 */
MappedCloneSet * Wise2_MappedCloneSet_alloc_len(int len);
#define MappedCloneSet_alloc_len Wise2_MappedCloneSet_alloc_len


/* Function:  hard_link_MappedCloneSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [MappedCloneSet *]
 *
 * Return [UNKN ]  Undocumented return value [MappedCloneSet *]
 *
 */
MappedCloneSet * Wise2_hard_link_MappedCloneSet(MappedCloneSet * obj);
#define hard_link_MappedCloneSet Wise2_hard_link_MappedCloneSet


/* Function:  MappedCloneSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [MappedCloneSet *]
 *
 */
MappedCloneSet * Wise2_MappedCloneSet_alloc(void);
#define MappedCloneSet_alloc Wise2_MappedCloneSet_alloc


/* Function:  free_MappedCloneSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [MappedCloneSet *]
 *
 * Return [UNKN ]  Undocumented return value [MappedCloneSet *]
 *
 */
MappedCloneSet * Wise2_free_MappedCloneSet(MappedCloneSet * obj);
#define free_MappedCloneSet Wise2_free_MappedCloneSet


/* Function:  MappedCloneMatch_alloc_matrix(leni,lenj)
 *
 * Descrip:    Allocates structure and matrix
 *
 *
 * Arg:        leni [UNKN ] Length of first dimension of matrix [int]
 * Arg:        lenj [UNKN ] Length of second dimension of matrix [int]
 *
 * Return [UNKN ]  Undocumented return value [MappedCloneMatch *]
 *
 */
MappedCloneMatch * Wise2_MappedCloneMatch_alloc_matrix(int leni,int lenj);
#define MappedCloneMatch_alloc_matrix Wise2_MappedCloneMatch_alloc_matrix


/* Function:  hard_link_MappedCloneMatch(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [MappedCloneMatch *]
 *
 * Return [UNKN ]  Undocumented return value [MappedCloneMatch *]
 *
 */
MappedCloneMatch * Wise2_hard_link_MappedCloneMatch(MappedCloneMatch * obj);
#define hard_link_MappedCloneMatch Wise2_hard_link_MappedCloneMatch


/* Function:  MappedCloneMatch_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [MappedCloneMatch *]
 *
 */
MappedCloneMatch * Wise2_MappedCloneMatch_alloc(void);
#define MappedCloneMatch_alloc Wise2_MappedCloneMatch_alloc


/* Function:  free_MappedCloneMatch(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [MappedCloneMatch *]
 *
 * Return [UNKN ]  Undocumented return value [MappedCloneMatch *]
 *
 */
MappedCloneMatch * Wise2_free_MappedCloneMatch(MappedCloneMatch * obj);
#define free_MappedCloneMatch Wise2_free_MappedCloneMatch


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
MappedCloneMatch * Wise2_new_MappedCloneMatch(MappedCloneSet * iset,MappedCloneSet * jset,int match,int mismatch);
#define new_MappedCloneMatch Wise2_new_MappedCloneMatch
int Wise2_MappedCloneSet_skip(MappedCloneSet * s,int pos,int skip_cost);
#define MappedCloneSet_skip Wise2_MappedCloneSet_skip
int Wise2_MappedCloneSet_match(MappedCloneSet * weak_query,MappedCloneSet * trusted_target,int qpos,int tpos,int spread,int match,int mismatch);
#define MappedCloneSet_match Wise2_MappedCloneSet_match
int Wise2_old_MappedCloneSet_match(MappedCloneSet * one,MappedCloneSet * two,int qpos,int tpos,int spread,int match,int mismatch);
#define old_MappedCloneSet_match Wise2_old_MappedCloneSet_match


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_MappedCloneSet(MappedClone ** list,int i,int j) ;
#define swap_MappedCloneSet Wise2_swap_MappedCloneSet
void Wise2_qsort_MappedCloneSet(MappedClone ** list,int left,int right,int (*comp)(MappedClone * ,MappedClone * ));
#define qsort_MappedCloneSet Wise2_qsort_MappedCloneSet
void Wise2_sort_MappedCloneSet(MappedCloneSet * obj,int (*comp)(MappedClone *, MappedClone *));
#define sort_MappedCloneSet Wise2_sort_MappedCloneSet
boolean Wise2_expand_MappedCloneSet(MappedCloneSet * obj,int len);
#define expand_MappedCloneSet Wise2_expand_MappedCloneSet
boolean Wise2_expand_MappedCloneMatch(MappedCloneMatch * obj,int leni,int lenj);
#define expand_MappedCloneMatch Wise2_expand_MappedCloneMatch

#ifdef _cplusplus
}
#endif

#endif
