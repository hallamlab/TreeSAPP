#ifndef DYNAMITEseqhitHEADERFILE
#define DYNAMITEseqhitHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"

#define SegmentHitListLISTLENGTH 128
#ifndef DYNAMITE_DEFINED_SegmentHit
typedef struct Wise2_SegmentHit Wise2_SegmentHit;
#define SegmentHit Wise2_SegmentHit
#define DYNAMITE_DEFINED_SegmentHit
#endif

struct Wise2_SegmentHit {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * name;     
    int qstart;  
    int qend;    
    int tstart;  
    int tend;    
    double score;    
    SegmentHit * next_hit;  /*  for collecting hits together ;) */ 
    } ;  
/* SegmentHit defined */ 
#ifndef DYNAMITE_DEFINED_SegmentHit
typedef struct Wise2_SegmentHit Wise2_SegmentHit;
#define SegmentHit Wise2_SegmentHit
#define DYNAMITE_DEFINED_SegmentHit
#endif


struct Wise2_SequenceHit {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * name;     
    SegmentHit * start;  
    } ;  
/* SequenceHit defined */ 
#ifndef DYNAMITE_DEFINED_SequenceHit
typedef struct Wise2_SequenceHit Wise2_SequenceHit;
#define SequenceHit Wise2_SequenceHit
#define DYNAMITE_DEFINED_SequenceHit
#endif


struct Wise2_SegmentHitList {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    SegmentHit ** seghit;    
    int len;/* len for above seghit  */ 
    int maxlen; /* maxlen for above seghit */ 
    } ;  
/* SegmentHitList defined */ 
#ifndef DYNAMITE_DEFINED_SegmentHitList
typedef struct Wise2_SegmentHitList Wise2_SegmentHitList;
#define SegmentHitList Wise2_SegmentHitList
#define DYNAMITE_DEFINED_SegmentHitList
#endif


struct Wise2_DnaSequenceHitList {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    SegmentHitList * forward;    
    SegmentHitList * backward;   
    } ;  
/* DnaSequenceHitList defined */ 
#ifndef DYNAMITE_DEFINED_DnaSequenceHitList
typedef struct Wise2_DnaSequenceHitList Wise2_DnaSequenceHitList;
#define DnaSequenceHitList Wise2_DnaSequenceHitList
#define DYNAMITE_DEFINED_DnaSequenceHitList
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  show_DnaSequenceHitList(dsl,ofp)
 *
 * Descrip:    shows a DnaSequenceHitsList -
 *
 *             only really useful for debugging
 *
 *
 * Arg:        dsl [UNKN ] Undocumented argument [DnaSequenceHitList *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_show_DnaSequenceHitList(DnaSequenceHitList * dsl,FILE * ofp);
#define show_DnaSequenceHitList Wise2_show_DnaSequenceHitList


/* Function:  read_MSPcrunch_DnaSequenceHitList(ifp)
 *
 * Descrip:    Reads a MSPcrunch -x output file 
 *
 *
 * Arg:        ifp [UNKN ] input file to read [FILE *]
 *
 * Return [UNKN ]  newly allocated structure [DnaSequenceHitList *]
 *
 */
DnaSequenceHitList * Wise2_read_MSPcrunch_DnaSequenceHitList(FILE * ifp);
#define read_MSPcrunch_DnaSequenceHitList Wise2_read_MSPcrunch_DnaSequenceHitList


/* Function:  hard_link_SegmentHit(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SegmentHit *]
 *
 * Return [UNKN ]  Undocumented return value [SegmentHit *]
 *
 */
SegmentHit * Wise2_hard_link_SegmentHit(SegmentHit * obj);
#define hard_link_SegmentHit Wise2_hard_link_SegmentHit


/* Function:  SegmentHit_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SegmentHit *]
 *
 */
SegmentHit * Wise2_SegmentHit_alloc(void);
#define SegmentHit_alloc Wise2_SegmentHit_alloc


/* Function:  free_SegmentHit(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SegmentHit *]
 *
 * Return [UNKN ]  Undocumented return value [SegmentHit *]
 *
 */
SegmentHit * Wise2_free_SegmentHit(SegmentHit * obj);
#define free_SegmentHit Wise2_free_SegmentHit


/* Function:  hard_link_SequenceHit(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SequenceHit *]
 *
 * Return [UNKN ]  Undocumented return value [SequenceHit *]
 *
 */
SequenceHit * Wise2_hard_link_SequenceHit(SequenceHit * obj);
#define hard_link_SequenceHit Wise2_hard_link_SequenceHit


/* Function:  SequenceHit_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SequenceHit *]
 *
 */
SequenceHit * Wise2_SequenceHit_alloc(void);
#define SequenceHit_alloc Wise2_SequenceHit_alloc


/* Function:  free_SequenceHit(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SequenceHit *]
 *
 * Return [UNKN ]  Undocumented return value [SequenceHit *]
 *
 */
SequenceHit * Wise2_free_SequenceHit(SequenceHit * obj);
#define free_SequenceHit Wise2_free_SequenceHit


/* Function:  add_SegmentHitList(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [SegmentHitList *]
 * Arg:        add [OWNER] Object to add to the list [SegmentHit *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_SegmentHitList(SegmentHitList * obj,SegmentHit * add);
#define add_SegmentHitList Wise2_add_SegmentHitList


/* Function:  flush_SegmentHitList(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [SegmentHitList *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_SegmentHitList(SegmentHitList * obj);
#define flush_SegmentHitList Wise2_flush_SegmentHitList


/* Function:  SegmentHitList_alloc_std(void)
 *
 * Descrip:    Equivalent to SegmentHitList_alloc_len(SegmentHitListLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SegmentHitList *]
 *
 */
SegmentHitList * Wise2_SegmentHitList_alloc_std(void);
#define SegmentHitList_alloc_std Wise2_SegmentHitList_alloc_std


/* Function:  SegmentHitList_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [SegmentHitList *]
 *
 */
SegmentHitList * Wise2_SegmentHitList_alloc_len(int len);
#define SegmentHitList_alloc_len Wise2_SegmentHitList_alloc_len


/* Function:  hard_link_SegmentHitList(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SegmentHitList *]
 *
 * Return [UNKN ]  Undocumented return value [SegmentHitList *]
 *
 */
SegmentHitList * Wise2_hard_link_SegmentHitList(SegmentHitList * obj);
#define hard_link_SegmentHitList Wise2_hard_link_SegmentHitList


/* Function:  SegmentHitList_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SegmentHitList *]
 *
 */
SegmentHitList * Wise2_SegmentHitList_alloc(void);
#define SegmentHitList_alloc Wise2_SegmentHitList_alloc


/* Function:  free_SegmentHitList(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SegmentHitList *]
 *
 * Return [UNKN ]  Undocumented return value [SegmentHitList *]
 *
 */
SegmentHitList * Wise2_free_SegmentHitList(SegmentHitList * obj);
#define free_SegmentHitList Wise2_free_SegmentHitList


/* Function:  hard_link_DnaSequenceHitList(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DnaSequenceHitList *]
 *
 * Return [UNKN ]  Undocumented return value [DnaSequenceHitList *]
 *
 */
DnaSequenceHitList * Wise2_hard_link_DnaSequenceHitList(DnaSequenceHitList * obj);
#define hard_link_DnaSequenceHitList Wise2_hard_link_DnaSequenceHitList


/* Function:  DnaSequenceHitList_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaSequenceHitList *]
 *
 */
DnaSequenceHitList * Wise2_DnaSequenceHitList_alloc(void);
#define DnaSequenceHitList_alloc Wise2_DnaSequenceHitList_alloc


/* Function:  free_DnaSequenceHitList(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DnaSequenceHitList *]
 *
 * Return [UNKN ]  Undocumented return value [DnaSequenceHitList *]
 *
 */
DnaSequenceHitList * Wise2_free_DnaSequenceHitList(DnaSequenceHitList * obj);
#define free_DnaSequenceHitList Wise2_free_DnaSequenceHitList


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
DnaSequenceHitList * Wise2_new_DnaSequenceHitList(void);
#define new_DnaSequenceHitList Wise2_new_DnaSequenceHitList
void Wise2_show_SegmentHitList(SegmentHitList * shl,FILE * ofp);
#define show_SegmentHitList Wise2_show_SegmentHitList
void Wise2_show_SegmentHit(SegmentHit * sh,FILE * ofp);
#define show_SegmentHit Wise2_show_SegmentHit
void Wise2_sort_SegmentHitList_by_qstart(SegmentHitList * sgl);
#define sort_SegmentHitList_by_qstart Wise2_sort_SegmentHitList_by_qstart
int Wise2_compare_SegmentHit_by_qstart(SegmentHit * one,SegmentHit * two);
#define compare_SegmentHit_by_qstart Wise2_compare_SegmentHit_by_qstart
DnaSequenceHitList * Wise2_read_msptmp_DnaSequenceHitList(FILE * ifp);
#define read_msptmp_DnaSequenceHitList Wise2_read_msptmp_DnaSequenceHitList
boolean Wise2_verify_SegmentHit(SegmentHit * sh);
#define verify_SegmentHit Wise2_verify_SegmentHit
SegmentHit * Wise2_SegmentHit_from_msptmp_line(char * line,int * strand);
#define SegmentHit_from_msptmp_line Wise2_SegmentHit_from_msptmp_line
SegmentHit * Wise2_SegmentHit_from_MSPcrunch_line(char * line,int * strand);
#define SegmentHit_from_MSPcrunch_line Wise2_SegmentHit_from_MSPcrunch_line


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
boolean Wise2_replace_forward_DnaSequenceHitList(DnaSequenceHitList * obj,SegmentHitList * forward);
#define replace_forward_DnaSequenceHitList Wise2_replace_forward_DnaSequenceHitList
boolean Wise2_replace_qstart_SegmentHit(SegmentHit * obj,int qstart);
#define replace_qstart_SegmentHit Wise2_replace_qstart_SegmentHit
SegmentHit * Wise2_access_next_hit_SegmentHit(SegmentHit * obj);
#define access_next_hit_SegmentHit Wise2_access_next_hit_SegmentHit
int Wise2_access_qstart_SegmentHit(SegmentHit * obj);
#define access_qstart_SegmentHit Wise2_access_qstart_SegmentHit
SegmentHitList * Wise2_access_forward_DnaSequenceHitList(DnaSequenceHitList * obj);
#define access_forward_DnaSequenceHitList Wise2_access_forward_DnaSequenceHitList
boolean Wise2_replace_qend_SegmentHit(SegmentHit * obj,int qend);
#define replace_qend_SegmentHit Wise2_replace_qend_SegmentHit
SegmentHitList * Wise2_access_backward_DnaSequenceHitList(DnaSequenceHitList * obj);
#define access_backward_DnaSequenceHitList Wise2_access_backward_DnaSequenceHitList
int Wise2_access_qend_SegmentHit(SegmentHit * obj);
#define access_qend_SegmentHit Wise2_access_qend_SegmentHit
int Wise2_length_seghit_SegmentHitList(SegmentHitList * obj);
#define length_seghit_SegmentHitList Wise2_length_seghit_SegmentHitList
boolean Wise2_replace_tstart_SegmentHit(SegmentHit * obj,int tstart);
#define replace_tstart_SegmentHit Wise2_replace_tstart_SegmentHit
char * Wise2_access_name_SegmentHit(SegmentHit * obj);
#define access_name_SegmentHit Wise2_access_name_SegmentHit
int Wise2_access_tstart_SegmentHit(SegmentHit * obj);
#define access_tstart_SegmentHit Wise2_access_tstart_SegmentHit
boolean Wise2_replace_backward_DnaSequenceHitList(DnaSequenceHitList * obj,SegmentHitList * backward);
#define replace_backward_DnaSequenceHitList Wise2_replace_backward_DnaSequenceHitList
boolean Wise2_replace_tend_SegmentHit(SegmentHit * obj,int tend);
#define replace_tend_SegmentHit Wise2_replace_tend_SegmentHit
boolean Wise2_replace_name_SegmentHit(SegmentHit * obj,char * name);
#define replace_name_SegmentHit Wise2_replace_name_SegmentHit
int Wise2_access_tend_SegmentHit(SegmentHit * obj);
#define access_tend_SegmentHit Wise2_access_tend_SegmentHit
SegmentHit * Wise2_access_seghit_SegmentHitList(SegmentHitList * obj,int i);
#define access_seghit_SegmentHitList Wise2_access_seghit_SegmentHitList
boolean Wise2_replace_score_SegmentHit(SegmentHit * obj,double score);
#define replace_score_SegmentHit Wise2_replace_score_SegmentHit
boolean Wise2_replace_next_hit_SegmentHit(SegmentHit * obj,SegmentHit * next_hit);
#define replace_next_hit_SegmentHit Wise2_replace_next_hit_SegmentHit
double Wise2_access_score_SegmentHit(SegmentHit * obj);
#define access_score_SegmentHit Wise2_access_score_SegmentHit
void Wise2_swap_SegmentHitList(SegmentHit ** list,int i,int j) ;
#define swap_SegmentHitList Wise2_swap_SegmentHitList
void Wise2_qsort_SegmentHitList(SegmentHit ** list,int left,int right,int (*comp)(SegmentHit * ,SegmentHit * ));
#define qsort_SegmentHitList Wise2_qsort_SegmentHitList
void Wise2_sort_SegmentHitList(SegmentHitList * obj,int (*comp)(SegmentHit *, SegmentHit *));
#define sort_SegmentHitList Wise2_sort_SegmentHitList
boolean Wise2_expand_SegmentHitList(SegmentHitList * obj,int len);
#define expand_SegmentHitList Wise2_expand_SegmentHitList

#ifdef _cplusplus
}
#endif

#endif
