#ifndef DYNAMITEsinglenumberspaceHEADERFILE
#define DYNAMITEsinglenumberspaceHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "shadowseq.h"

#define SingleNumberSpaceLISTLENGTH 1024

struct Wise2_SingleNumberSequence {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int start;   
    int end;     
    ShadowSequence * seq;    
    } ;  
/* SingleNumberSequence defined */ 
#ifndef DYNAMITE_DEFINED_SingleNumberSequence
typedef struct Wise2_SingleNumberSequence Wise2_SingleNumberSequence;
#define SingleNumberSequence Wise2_SingleNumberSequence
#define DYNAMITE_DEFINED_SingleNumberSequence
#endif


struct Wise2_SingleNumberSpace {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int current_end;     
    SingleNumberSequence ** sns;     
    int len;/* len for above sns  */ 
    int maxlen; /* maxlen for above sns */ 
    SingleNumberSequence * last_accessed;    
    int average_len;     
    int is_static;   
    int max_length;  
    } ;  
/* SingleNumberSpace defined */ 
#ifndef DYNAMITE_DEFINED_SingleNumberSpace
typedef struct Wise2_SingleNumberSpace Wise2_SingleNumberSpace;
#define SingleNumberSpace Wise2_SingleNumberSpace
#define DYNAMITE_DEFINED_SingleNumberSpace
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  lookup_ShadowSequence_SingleNumberSpace(space,pos)
 *
 * Descrip:    New return using binary choping
 *
 *
 * Arg:        space [UNKN ] Undocumented argument [SingleNumberSpace *]
 * Arg:          pos [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [SingleNumberSequence *]
 *
 */
SingleNumberSequence * Wise2_lookup_ShadowSequence_SingleNumberSpace(SingleNumberSpace * space,int pos);
#define lookup_ShadowSequence_SingleNumberSpace Wise2_lookup_ShadowSequence_SingleNumberSpace


/* Function:  find_position_SingleNumberSpace(space,lower,higher,position)
 *
 * Descrip:    Recursive function for finding position
 *
 *
 * Arg:           space [UNKN ] Undocumented argument [SingleNumberSpace *]
 * Arg:           lower [UNKN ] Undocumented argument [int]
 * Arg:          higher [UNKN ] Undocumented argument [int]
 * Arg:        position [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_find_position_SingleNumberSpace(SingleNumberSpace * space,int lower,int higher,int position);
#define find_position_SingleNumberSpace Wise2_find_position_SingleNumberSpace


/* Function:  add_ShadowSequence_SingleNumberSpace(space,seq)
 *
 * Descrip:    Adds a sequence to a single number space, giving out the start
 *             position for this sequence
 *
 *
 * Arg:        space [UNKN ] Undocumented argument [SingleNumberSpace *]
 * Arg:          seq [UNKN ] Undocumented argument [ShadowSequence *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_add_ShadowSequence_SingleNumberSpace(SingleNumberSpace * space,ShadowSequence * seq);
#define add_ShadowSequence_SingleNumberSpace Wise2_add_ShadowSequence_SingleNumberSpace


/* Function:  new_SingleNumberSpace(has_maxlen,max_length)
 *
 * Descrip:    Provides a new single number space
 *
 *
 * Arg:        has_maxlen [UNKN ] Undocumented argument [int]
 * Arg:        max_length [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [SingleNumberSpace *]
 *
 */
SingleNumberSpace * Wise2_new_SingleNumberSpace(int has_maxlen,int max_length);
#define new_SingleNumberSpace Wise2_new_SingleNumberSpace


/* Function:  hard_link_SingleNumberSequence(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SingleNumberSequence *]
 *
 * Return [UNKN ]  Undocumented return value [SingleNumberSequence *]
 *
 */
SingleNumberSequence * Wise2_hard_link_SingleNumberSequence(SingleNumberSequence * obj);
#define hard_link_SingleNumberSequence Wise2_hard_link_SingleNumberSequence


/* Function:  SingleNumberSequence_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SingleNumberSequence *]
 *
 */
SingleNumberSequence * Wise2_SingleNumberSequence_alloc(void);
#define SingleNumberSequence_alloc Wise2_SingleNumberSequence_alloc


/* Function:  free_SingleNumberSequence(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SingleNumberSequence *]
 *
 * Return [UNKN ]  Undocumented return value [SingleNumberSequence *]
 *
 */
SingleNumberSequence * Wise2_free_SingleNumberSequence(SingleNumberSequence * obj);
#define free_SingleNumberSequence Wise2_free_SingleNumberSequence


/* Function:  add_SingleNumberSpace(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [SingleNumberSpace *]
 * Arg:        add [OWNER] Object to add to the list [SingleNumberSequence *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_SingleNumberSpace(SingleNumberSpace * obj,SingleNumberSequence * add);
#define add_SingleNumberSpace Wise2_add_SingleNumberSpace


/* Function:  flush_SingleNumberSpace(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [SingleNumberSpace *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_SingleNumberSpace(SingleNumberSpace * obj);
#define flush_SingleNumberSpace Wise2_flush_SingleNumberSpace


/* Function:  SingleNumberSpace_alloc_std(void)
 *
 * Descrip:    Equivalent to SingleNumberSpace_alloc_len(SingleNumberSpaceLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SingleNumberSpace *]
 *
 */
SingleNumberSpace * Wise2_SingleNumberSpace_alloc_std(void);
#define SingleNumberSpace_alloc_std Wise2_SingleNumberSpace_alloc_std


/* Function:  SingleNumberSpace_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [SingleNumberSpace *]
 *
 */
SingleNumberSpace * Wise2_SingleNumberSpace_alloc_len(int len);
#define SingleNumberSpace_alloc_len Wise2_SingleNumberSpace_alloc_len


/* Function:  hard_link_SingleNumberSpace(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SingleNumberSpace *]
 *
 * Return [UNKN ]  Undocumented return value [SingleNumberSpace *]
 *
 */
SingleNumberSpace * Wise2_hard_link_SingleNumberSpace(SingleNumberSpace * obj);
#define hard_link_SingleNumberSpace Wise2_hard_link_SingleNumberSpace


/* Function:  SingleNumberSpace_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SingleNumberSpace *]
 *
 */
SingleNumberSpace * Wise2_SingleNumberSpace_alloc(void);
#define SingleNumberSpace_alloc Wise2_SingleNumberSpace_alloc


/* Function:  free_SingleNumberSpace(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SingleNumberSpace *]
 *
 * Return [UNKN ]  Undocumented return value [SingleNumberSpace *]
 *
 */
SingleNumberSpace * Wise2_free_SingleNumberSpace(SingleNumberSpace * obj);
#define free_SingleNumberSpace Wise2_free_SingleNumberSpace


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_SingleNumberSpace(SingleNumberSequence ** list,int i,int j) ;
#define swap_SingleNumberSpace Wise2_swap_SingleNumberSpace
void Wise2_qsort_SingleNumberSpace(SingleNumberSequence ** list,int left,int right,int (*comp)(SingleNumberSequence * ,SingleNumberSequence * ));
#define qsort_SingleNumberSpace Wise2_qsort_SingleNumberSpace
void Wise2_sort_SingleNumberSpace(SingleNumberSpace * obj,int (*comp)(SingleNumberSequence *, SingleNumberSequence *));
#define sort_SingleNumberSpace Wise2_sort_SingleNumberSpace
boolean Wise2_expand_SingleNumberSpace(SingleNumberSpace * obj,int len);
#define expand_SingleNumberSpace Wise2_expand_SingleNumberSpace

#ifdef _cplusplus
}
#endif

#endif
