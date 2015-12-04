#ifndef DYNAMITEsingleseqspaceHEADERFILE
#define DYNAMITEsingleseqspaceHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "sequence.h"

#define SinglePosSpaceLISTLENGTH 1024

struct Wise2_SinglePosSequence {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    long start;  
    long end;    
    void * data;     
    } ;  
/* SinglePosSequence defined */ 
#ifndef DYNAMITE_DEFINED_SinglePosSequence
typedef struct Wise2_SinglePosSequence Wise2_SinglePosSequence;
#define SinglePosSequence Wise2_SinglePosSequence
#define DYNAMITE_DEFINED_SinglePosSequence
#endif


struct Wise2_SinglePosSpace {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    long current_end;    
    SinglePosSequence ** sns;    
    int len;/* len for above sns  */ 
    int maxlen; /* maxlen for above sns */ 
    SinglePosSequence * last_accessed;   
    int average_len;     
    int is_static;   
    int max_length;  
    } ;  
/* SinglePosSpace defined */ 
#ifndef DYNAMITE_DEFINED_SinglePosSpace
typedef struct Wise2_SinglePosSpace Wise2_SinglePosSpace;
#define SinglePosSpace Wise2_SinglePosSpace
#define DYNAMITE_DEFINED_SinglePosSpace
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  lookup_Sequence_SinglePosSpace(space,pos)
 *
 * Descrip:    New return using binary choping
 *
 *
 * Arg:        space [UNKN ] Undocumented argument [SinglePosSpace *]
 * Arg:          pos [UNKN ] Undocumented argument [long]
 *
 * Return [UNKN ]  Undocumented return value [SinglePosSequence *]
 *
 */
SinglePosSequence * Wise2_lookup_Sequence_SinglePosSpace(SinglePosSpace * space,long pos);
#define lookup_Sequence_SinglePosSpace Wise2_lookup_Sequence_SinglePosSpace


/* Function:  find_position_SinglePosSpace(space,lower,higher,position)
 *
 * Descrip:    Recursive function for finding position
 *
 *
 * Arg:           space [UNKN ] Undocumented argument [SinglePosSpace *]
 * Arg:           lower [UNKN ] Undocumented argument [int]
 * Arg:          higher [UNKN ] Undocumented argument [int]
 * Arg:        position [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_find_position_SinglePosSpace(SinglePosSpace * space,int lower,int higher,int position);
#define find_position_SinglePosSpace Wise2_find_position_SinglePosSpace


/* Function:  add_Sequence_SinglePosSpace(space,length,data)
 *
 * Descrip:    Adds a sequence to a single number space, giving out the start
 *             position for this sequence
 *
 *
 * Arg:         space [UNKN ] Undocumented argument [SinglePosSpace *]
 * Arg:        length [UNKN ] Undocumented argument [long int]
 * Arg:          data [UNKN ] Undocumented argument [void *]
 *
 * Return [UNKN ]  Undocumented return value [long int]
 *
 */
long int Wise2_add_Sequence_SinglePosSpace(SinglePosSpace * space,long int length,void * data);
#define add_Sequence_SinglePosSpace Wise2_add_Sequence_SinglePosSpace


/* Function:  new_SinglePosSpace(has_maxlen,max_length)
 *
 * Descrip:    Provides a new single number space
 *
 *
 * Arg:        has_maxlen [UNKN ] Undocumented argument [int]
 * Arg:        max_length [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [SinglePosSpace *]
 *
 */
SinglePosSpace * Wise2_new_SinglePosSpace(int has_maxlen,int max_length);
#define new_SinglePosSpace Wise2_new_SinglePosSpace


/* Function:  hard_link_SinglePosSequence(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SinglePosSequence *]
 *
 * Return [UNKN ]  Undocumented return value [SinglePosSequence *]
 *
 */
SinglePosSequence * Wise2_hard_link_SinglePosSequence(SinglePosSequence * obj);
#define hard_link_SinglePosSequence Wise2_hard_link_SinglePosSequence


/* Function:  SinglePosSequence_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SinglePosSequence *]
 *
 */
SinglePosSequence * Wise2_SinglePosSequence_alloc(void);
#define SinglePosSequence_alloc Wise2_SinglePosSequence_alloc


/* Function:  free_SinglePosSequence(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SinglePosSequence *]
 *
 * Return [UNKN ]  Undocumented return value [SinglePosSequence *]
 *
 */
SinglePosSequence * Wise2_free_SinglePosSequence(SinglePosSequence * obj);
#define free_SinglePosSequence Wise2_free_SinglePosSequence


/* Function:  add_SinglePosSpace(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [SinglePosSpace *]
 * Arg:        add [OWNER] Object to add to the list [SinglePosSequence *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_SinglePosSpace(SinglePosSpace * obj,SinglePosSequence * add);
#define add_SinglePosSpace Wise2_add_SinglePosSpace


/* Function:  flush_SinglePosSpace(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [SinglePosSpace *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_SinglePosSpace(SinglePosSpace * obj);
#define flush_SinglePosSpace Wise2_flush_SinglePosSpace


/* Function:  SinglePosSpace_alloc_std(void)
 *
 * Descrip:    Equivalent to SinglePosSpace_alloc_len(SinglePosSpaceLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SinglePosSpace *]
 *
 */
SinglePosSpace * Wise2_SinglePosSpace_alloc_std(void);
#define SinglePosSpace_alloc_std Wise2_SinglePosSpace_alloc_std


/* Function:  SinglePosSpace_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [SinglePosSpace *]
 *
 */
SinglePosSpace * Wise2_SinglePosSpace_alloc_len(int len);
#define SinglePosSpace_alloc_len Wise2_SinglePosSpace_alloc_len


/* Function:  hard_link_SinglePosSpace(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SinglePosSpace *]
 *
 * Return [UNKN ]  Undocumented return value [SinglePosSpace *]
 *
 */
SinglePosSpace * Wise2_hard_link_SinglePosSpace(SinglePosSpace * obj);
#define hard_link_SinglePosSpace Wise2_hard_link_SinglePosSpace


/* Function:  SinglePosSpace_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SinglePosSpace *]
 *
 */
SinglePosSpace * Wise2_SinglePosSpace_alloc(void);
#define SinglePosSpace_alloc Wise2_SinglePosSpace_alloc


/* Function:  free_SinglePosSpace(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SinglePosSpace *]
 *
 * Return [UNKN ]  Undocumented return value [SinglePosSpace *]
 *
 */
SinglePosSpace * Wise2_free_SinglePosSpace(SinglePosSpace * obj);
#define free_SinglePosSpace Wise2_free_SinglePosSpace


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_SinglePosSpace(SinglePosSequence ** list,int i,int j) ;
#define swap_SinglePosSpace Wise2_swap_SinglePosSpace
void Wise2_qsort_SinglePosSpace(SinglePosSequence ** list,int left,int right,int (*comp)(SinglePosSequence * ,SinglePosSequence * ));
#define qsort_SinglePosSpace Wise2_qsort_SinglePosSpace
void Wise2_sort_SinglePosSpace(SinglePosSpace * obj,int (*comp)(SinglePosSequence *, SinglePosSequence *));
#define sort_SinglePosSpace Wise2_sort_SinglePosSpace
boolean Wise2_expand_SinglePosSpace(SinglePosSpace * obj,int len);
#define expand_SinglePosSpace Wise2_expand_SinglePosSpace

#ifdef _cplusplus
}
#endif

#endif
