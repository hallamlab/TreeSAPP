#ifndef DYNAMITEdpenvelopeHEADERFILE
#define DYNAMITEdpenvelopeHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "wisebase.h"

typedef enum dpenvelope_type {
  DPENV_RECT = 0,
  DPENV_DIAG
} dpenv_type;



#define DPEnvelopeLISTLENGTH 32
struct Wise2_DPUnit {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int type;    
    int starti;  
    int startj;  
    int height; /*  for diagonal units */ 
    int length; /*  for diagonal units */ 
    } ;  
/* DPUnit defined */ 
#ifndef DYNAMITE_DEFINED_DPUnit
typedef struct Wise2_DPUnit Wise2_DPUnit;
#define DPUnit Wise2_DPUnit
#define DYNAMITE_DEFINED_DPUnit
#endif


struct Wise2_DPEnvelope {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    DPUnit ** dpu;   
    int len;/* len for above dpu  */ 
    int maxlen; /* maxlen for above dpu */ 
    DPUnit * bbox;   
    int starti;  
    int startj;  
    int endi;    
    int endj;    
    } ;  
/* DPEnvelope defined */ 
#ifndef DYNAMITE_DEFINED_DPEnvelope
typedef struct Wise2_DPEnvelope Wise2_DPEnvelope;
#define DPEnvelope Wise2_DPEnvelope
#define DYNAMITE_DEFINED_DPEnvelope
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  overlap_DPUnit(a,b)
 *
 * Descrip:    Helper function that checks whether things overlap or not
 *
 *
 * Arg:        a [UNKN ] Undocumented argument [DPUnit *]
 * Arg:        b [UNKN ] Undocumented argument [DPUnit *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_overlap_DPUnit(DPUnit * a,DPUnit * b);
#define overlap_DPUnit Wise2_overlap_DPUnit


/* Function:  read_DPEnvelope_file(filename)
 *
 * Descrip:    Helper function that also opens the filename
 *
 *
 * Arg:        filename [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [DPEnvelope *]
 *
 */
DPEnvelope * Wise2_read_DPEnvelope_file(char * filename);
#define read_DPEnvelope_file Wise2_read_DPEnvelope_file


/* Function:  read_DPEnvelope(ifp)
 *
 * Descrip:    Reads a DPEnvelope from a file
 *
 *
 * Arg:        ifp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [DPEnvelope *]
 *
 */
DPEnvelope * Wise2_read_DPEnvelope(FILE * ifp);
#define read_DPEnvelope Wise2_read_DPEnvelope


/* Function:  show_DPEnvelope(dpe,ofp)
 *
 * Descrip:    shows structure. useful for debugging
 *
 *
 * Arg:        dpe [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_show_DPEnvelope(DPEnvelope * dpe,FILE * ofp);
#define show_DPEnvelope Wise2_show_DPEnvelope


/* Function:  is_in_DPEnvelope(dpe,i,j)
 *
 * Descrip:    Tests whether this i,j position is allowed in the
 *             DPEnvelope
 *
 *
 * Arg:        dpe [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:          i [UNKN ] Undocumented argument [int]
 * Arg:          j [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_is_in_DPEnvelope(DPEnvelope * dpe,int i,int j);
#define is_in_DPEnvelope Wise2_is_in_DPEnvelope


/* Function:  prepare_DPEnvelope(dpe)
 *
 * Descrip:    Should run this before using the DPEnvelope
 *
 *
 * Arg:        dpe [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_prepare_DPEnvelope(DPEnvelope * dpe);
#define prepare_DPEnvelope Wise2_prepare_DPEnvelope


/* Function:  sort_DPEnvelope_by_startj(dpe)
 *
 * Descrip:    Sorts by startj
 *
 *
 * Arg:        dpe [UNKN ] Undocumented argument [DPEnvelope *]
 *
 */
void Wise2_sort_DPEnvelope_by_startj(DPEnvelope * dpe);
#define sort_DPEnvelope_by_startj Wise2_sort_DPEnvelope_by_startj


/* Function:  hard_link_DPUnit(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DPUnit *]
 *
 * Return [UNKN ]  Undocumented return value [DPUnit *]
 *
 */
DPUnit * Wise2_hard_link_DPUnit(DPUnit * obj);
#define hard_link_DPUnit Wise2_hard_link_DPUnit


/* Function:  DPUnit_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DPUnit *]
 *
 */
DPUnit * Wise2_DPUnit_alloc(void);
#define DPUnit_alloc Wise2_DPUnit_alloc


/* Function:  free_DPUnit(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DPUnit *]
 *
 * Return [UNKN ]  Undocumented return value [DPUnit *]
 *
 */
DPUnit * Wise2_free_DPUnit(DPUnit * obj);
#define free_DPUnit Wise2_free_DPUnit


/* Function:  add_DPEnvelope(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [DPEnvelope *]
 * Arg:        add [OWNER] Object to add to the list [DPUnit *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_DPEnvelope(DPEnvelope * obj,DPUnit * add);
#define add_DPEnvelope Wise2_add_DPEnvelope


/* Function:  flush_DPEnvelope(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_DPEnvelope(DPEnvelope * obj);
#define flush_DPEnvelope Wise2_flush_DPEnvelope


/* Function:  DPEnvelope_alloc_std(void)
 *
 * Descrip:    Equivalent to DPEnvelope_alloc_len(DPEnvelopeLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DPEnvelope *]
 *
 */
DPEnvelope * Wise2_DPEnvelope_alloc_std(void);
#define DPEnvelope_alloc_std Wise2_DPEnvelope_alloc_std


/* Function:  DPEnvelope_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [DPEnvelope *]
 *
 */
DPEnvelope * Wise2_DPEnvelope_alloc_len(int len);
#define DPEnvelope_alloc_len Wise2_DPEnvelope_alloc_len


/* Function:  hard_link_DPEnvelope(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [DPEnvelope *]
 *
 */
DPEnvelope * Wise2_hard_link_DPEnvelope(DPEnvelope * obj);
#define hard_link_DPEnvelope Wise2_hard_link_DPEnvelope


/* Function:  DPEnvelope_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DPEnvelope *]
 *
 */
DPEnvelope * Wise2_DPEnvelope_alloc(void);
#define DPEnvelope_alloc Wise2_DPEnvelope_alloc


/* Function:  free_DPEnvelope(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [DPEnvelope *]
 *
 */
DPEnvelope * Wise2_free_DPEnvelope(DPEnvelope * obj);
#define free_DPEnvelope Wise2_free_DPEnvelope


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
int Wise2_compare_DPUnit_startj(DPUnit * one,DPUnit * two);
#define compare_DPUnit_startj Wise2_compare_DPUnit_startj
void Wise2_swap_DPEnvelope(DPUnit ** list,int i,int j) ;
#define swap_DPEnvelope Wise2_swap_DPEnvelope
void Wise2_qsort_DPEnvelope(DPUnit ** list,int left,int right,int (*comp)(DPUnit * ,DPUnit * ));
#define qsort_DPEnvelope Wise2_qsort_DPEnvelope
void Wise2_sort_DPEnvelope(DPEnvelope * obj,int (*comp)(DPUnit *, DPUnit *));
#define sort_DPEnvelope Wise2_sort_DPEnvelope
boolean Wise2_expand_DPEnvelope(DPEnvelope * obj,int len);
#define expand_DPEnvelope Wise2_expand_DPEnvelope

#ifdef _cplusplus
}
#endif

#endif
