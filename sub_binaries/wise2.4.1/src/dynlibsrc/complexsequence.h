#ifndef DYNAMITEcomplexsequenceHEADERFILE
#define DYNAMITEcomplexsequenceHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "sequence.h"
#include "probability.h" /* should move this to wisestring */


#define ComplexSequenceEvalSetLISTLENGTH 16

#define next_ComplexSequence_data(cs_struct,pointer) (pointer + cs_struct->depth)
#define ComplexSequence_data(cs,position,number)     (cs->data[cs->depth*position + number])

typedef enum CseScoreType {
  CseScoreType_Index = 673,
  CseScoreType_Bits
} CseScoreType;

/* Object ComplexSequenceEval
 *
 * Descrip: This object is best left alone (!)
 *
 *        It represents a single way of mapping
 *        a sequence to some sort of number, eg
 *        amino acids to 0-25, bases to 0-4 or
 *        splice sites to log(Prob(splice))
 *
 *        This is handled by a reasonably scary
 *        pointer-to-function method
 *
 *        You'll use collections of them in 
 *        complexsequenceevalset's
 *
 *
 */
struct Wise2_ComplexSequenceEval {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int type;    
    int sequence_type;   
    int left_window;     
    int right_window;    
    int left_lookback;   
    int outside_score;   
    void * data;     
    int data_type;  /* optional, to be used by eval_func to check things if it wants to */ 
    int (*eval_func)(int,void * ,char *);    
    CseScoreType score_type;     
    } ;  
/* ComplexSequenceEval defined */ 
#ifndef DYNAMITE_DEFINED_ComplexSequenceEval
typedef struct Wise2_ComplexSequenceEval Wise2_ComplexSequenceEval;
#define ComplexSequenceEval Wise2_ComplexSequenceEval
#define DYNAMITE_DEFINED_ComplexSequenceEval
#endif


/* Object ComplexSequenceEvalSet
 *
 * Descrip: This object holds a collection of 
 *        ComplexSequenceEvals. Its role is to
 *        define the sequence specific parts of a
 *        dynamic programming algorithm as computable
 *        functions. 
 *
 *        Ideally you should use pre-made ComplexSequenceEvalSets
 *        as it will save you alot of grief
 *
 *
 */
struct Wise2_ComplexSequenceEvalSet {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int type;    
    boolean has_been_prepared;   
    int left_window;    /*  overall sequence eval  */ 
    int right_window;   /*  overall sequence eval */ 
    int left_lookback;  /*  overall sequence eval */ 
    ComplexSequenceEval ** cse;  
    int len;/* len for above cse  */ 
    int maxlen; /* maxlen for above cse */ 
    } ;  
/* ComplexSequenceEvalSet defined */ 
#ifndef DYNAMITE_DEFINED_ComplexSequenceEvalSet
typedef struct Wise2_ComplexSequenceEvalSet Wise2_ComplexSequenceEvalSet;
#define ComplexSequenceEvalSet Wise2_ComplexSequenceEvalSet
#define DYNAMITE_DEFINED_ComplexSequenceEvalSet
#endif


/* Object ComplexSequence
 *
 * Descrip: A ComplexSequence is an abstraction of a 
 *        Sequence which can be handily made using
 *        ComplexSequenceEval functions and is efficiently
 *        laid out in memory.
 *
 *
 */
struct Wise2_ComplexSequence {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int type;    
    Sequence * seq;  
    int * data;  
    int * datastore;     
    int depth;   
    int length;  
    ComplexSequenceEvalSet * creator;   /*  what made it */ 
    } ;  
/* ComplexSequence defined */ 
#ifndef DYNAMITE_DEFINED_ComplexSequence
typedef struct Wise2_ComplexSequence Wise2_ComplexSequence;
#define ComplexSequence Wise2_ComplexSequence
#define DYNAMITE_DEFINED_ComplexSequence
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  reversed_ComplexSequence(forward,cses)
 *
 * Descrip: No Description
 *
 * Arg:        forward [UNKN ] Undocumented argument [Sequence *]
 * Arg:           cses [UNKN ] Undocumented argument [ComplexSequenceEvalSet *]
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequence *]
 *
 */
ComplexSequence * Wise2_reversed_ComplexSequence(Sequence * forward,ComplexSequenceEvalSet * cses) ;
#define reversed_ComplexSequence Wise2_reversed_ComplexSequence


/* Function:  new_ComplexSequence(seq,cses)
 *
 * Descrip:    The basic way to make a ComplexSequence. Requires that
 *             you have already built a ComplexSequenceEvalSet (such as
 *             /default_aminoacid_ComplexSequenceEvalSet).
 *
 *
 *
 * Arg:         seq [UNKN ] Sequence that the ComplexSequence is based on [Sequence *]
 * Arg:        cses [UNKN ] EvalSet that defines the functions used on the sequence [ComplexSequenceEvalSet *]
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequence *]
 *
 */
ComplexSequence * Wise2_new_ComplexSequence(Sequence * seq,ComplexSequenceEvalSet * cses);
#define new_ComplexSequence Wise2_new_ComplexSequence


/* Function:  show_ComplexSequence(cs,ofp)
 *
 * Descrip:    shows complex sequence in a vaguely
 *             human form
 *
 *
 * Arg:         cs [UNKN ] Undocumented argument [ComplexSequence *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_show_ComplexSequence(ComplexSequence * cs,FILE * ofp);
#define show_ComplexSequence Wise2_show_ComplexSequence


/* Function:  prepare_ComplexSequenceEvalSet(cses)
 *
 * Descrip:    Calculates all the necessary offset for an EvalSet.
 *             This is necessary before using it in a /new_ComplexSequence
 *             place
 *
 *
 * Arg:        cses [UNKN ] Undocumented argument [ComplexSequenceEvalSet *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_prepare_ComplexSequenceEvalSet(ComplexSequenceEvalSet * cses);
#define prepare_ComplexSequenceEvalSet Wise2_prepare_ComplexSequenceEvalSet


/* Function:  can_evaluate_this_Sequence(cses,s)
 *
 * Descrip:    Checks that this ComplexSequenceEvalSet can be used with
 *             this Sequence. This is probably going to go defunct.
 *
 *
 *
 * Arg:        cses [UNKN ] Undocumented argument [ComplexSequenceEvalSet *]
 * Arg:           s [UNKN ] Undocumented argument [Sequence *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_can_evaluate_this_Sequence(ComplexSequenceEvalSet * cses,Sequence * s);
#define can_evaluate_this_Sequence Wise2_can_evaluate_this_Sequence


/* Function:  hard_link_ComplexSequenceEval(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ComplexSequenceEval *]
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequenceEval *]
 *
 */
ComplexSequenceEval * Wise2_hard_link_ComplexSequenceEval(ComplexSequenceEval * obj);
#define hard_link_ComplexSequenceEval Wise2_hard_link_ComplexSequenceEval


/* Function:  ComplexSequenceEval_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequenceEval *]
 *
 */
ComplexSequenceEval * Wise2_ComplexSequenceEval_alloc(void);
#define ComplexSequenceEval_alloc Wise2_ComplexSequenceEval_alloc


/* Function:  free_ComplexSequenceEval(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ComplexSequenceEval *]
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequenceEval *]
 *
 */
ComplexSequenceEval * Wise2_free_ComplexSequenceEval(ComplexSequenceEval * obj);
#define free_ComplexSequenceEval Wise2_free_ComplexSequenceEval


/* Function:  add_ComplexSequenceEvalSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [ComplexSequenceEvalSet *]
 * Arg:        add [OWNER] Object to add to the list [ComplexSequenceEval *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_ComplexSequenceEvalSet(ComplexSequenceEvalSet * obj,ComplexSequenceEval * add);
#define add_ComplexSequenceEvalSet Wise2_add_ComplexSequenceEvalSet


/* Function:  flush_ComplexSequenceEvalSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [ComplexSequenceEvalSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_ComplexSequenceEvalSet(ComplexSequenceEvalSet * obj);
#define flush_ComplexSequenceEvalSet Wise2_flush_ComplexSequenceEvalSet


/* Function:  ComplexSequenceEvalSet_alloc_std(void)
 *
 * Descrip:    Equivalent to ComplexSequenceEvalSet_alloc_len(ComplexSequenceEvalSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequenceEvalSet *]
 *
 */
ComplexSequenceEvalSet * Wise2_ComplexSequenceEvalSet_alloc_std(void);
#define ComplexSequenceEvalSet_alloc_std Wise2_ComplexSequenceEvalSet_alloc_std


/* Function:  ComplexSequenceEvalSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequenceEvalSet *]
 *
 */
ComplexSequenceEvalSet * Wise2_ComplexSequenceEvalSet_alloc_len(int len);
#define ComplexSequenceEvalSet_alloc_len Wise2_ComplexSequenceEvalSet_alloc_len


/* Function:  hard_link_ComplexSequenceEvalSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ComplexSequenceEvalSet *]
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequenceEvalSet *]
 *
 */
ComplexSequenceEvalSet * Wise2_hard_link_ComplexSequenceEvalSet(ComplexSequenceEvalSet * obj);
#define hard_link_ComplexSequenceEvalSet Wise2_hard_link_ComplexSequenceEvalSet


/* Function:  ComplexSequenceEvalSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequenceEvalSet *]
 *
 */
ComplexSequenceEvalSet * Wise2_ComplexSequenceEvalSet_alloc(void);
#define ComplexSequenceEvalSet_alloc Wise2_ComplexSequenceEvalSet_alloc


/* Function:  free_ComplexSequenceEvalSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ComplexSequenceEvalSet *]
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequenceEvalSet *]
 *
 */
ComplexSequenceEvalSet * Wise2_free_ComplexSequenceEvalSet(ComplexSequenceEvalSet * obj);
#define free_ComplexSequenceEvalSet Wise2_free_ComplexSequenceEvalSet


/* Function:  hard_link_ComplexSequence(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ComplexSequence *]
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequence *]
 *
 */
ComplexSequence * Wise2_hard_link_ComplexSequence(ComplexSequence * obj);
#define hard_link_ComplexSequence Wise2_hard_link_ComplexSequence


/* Function:  ComplexSequence_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequence *]
 *
 */
ComplexSequence * Wise2_ComplexSequence_alloc(void);
#define ComplexSequence_alloc Wise2_ComplexSequence_alloc


/* Function:  free_ComplexSequence(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ComplexSequence *]
 *
 * Return [UNKN ]  Undocumented return value [ComplexSequence *]
 *
 */
ComplexSequence * Wise2_free_ComplexSequence(ComplexSequence * obj);
#define free_ComplexSequence Wise2_free_ComplexSequence


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
int Wise2_access_type_ComplexSequence(ComplexSequence * obj);
#define access_type_ComplexSequence Wise2_access_type_ComplexSequence
boolean Wise2_replace_type_ComplexSequenceEvalSet(ComplexSequenceEvalSet * obj,int type);
#define replace_type_ComplexSequenceEvalSet Wise2_replace_type_ComplexSequenceEvalSet
boolean Wise2_replace_has_been_prepared_ComplexSequenceEvalSet(ComplexSequenceEvalSet * obj,boolean has_been_prepared);
#define replace_has_been_prepared_ComplexSequenceEvalSet Wise2_replace_has_been_prepared_ComplexSequenceEvalSet
int Wise2_access_left_lookback_ComplexSequenceEvalSet(ComplexSequenceEvalSet * obj);
#define access_left_lookback_ComplexSequenceEvalSet Wise2_access_left_lookback_ComplexSequenceEvalSet
boolean Wise2_access_has_been_prepared_ComplexSequenceEvalSet(ComplexSequenceEvalSet * obj);
#define access_has_been_prepared_ComplexSequenceEvalSet Wise2_access_has_been_prepared_ComplexSequenceEvalSet
boolean Wise2_replace_seq_ComplexSequence(ComplexSequence * obj,Sequence * seq);
#define replace_seq_ComplexSequence Wise2_replace_seq_ComplexSequence
boolean Wise2_replace_left_window_ComplexSequenceEvalSet(ComplexSequenceEvalSet * obj,int left_window);
#define replace_left_window_ComplexSequenceEvalSet Wise2_replace_left_window_ComplexSequenceEvalSet
int Wise2_access_type_ComplexSequenceEvalSet(ComplexSequenceEvalSet * obj);
#define access_type_ComplexSequenceEvalSet Wise2_access_type_ComplexSequenceEvalSet
boolean Wise2_replace_left_lookback_ComplexSequenceEvalSet(ComplexSequenceEvalSet * obj,int left_lookback);
#define replace_left_lookback_ComplexSequenceEvalSet Wise2_replace_left_lookback_ComplexSequenceEvalSet
int Wise2_access_left_window_ComplexSequenceEvalSet(ComplexSequenceEvalSet * obj);
#define access_left_window_ComplexSequenceEvalSet Wise2_access_left_window_ComplexSequenceEvalSet
Sequence * Wise2_access_seq_ComplexSequence(ComplexSequence * obj);
#define access_seq_ComplexSequence Wise2_access_seq_ComplexSequence
boolean Wise2_replace_right_window_ComplexSequenceEvalSet(ComplexSequenceEvalSet * obj,int right_window);
#define replace_right_window_ComplexSequenceEvalSet Wise2_replace_right_window_ComplexSequenceEvalSet
boolean Wise2_replace_type_ComplexSequence(ComplexSequence * obj,int type);
#define replace_type_ComplexSequence Wise2_replace_type_ComplexSequence
int Wise2_access_right_window_ComplexSequenceEvalSet(ComplexSequenceEvalSet * obj);
#define access_right_window_ComplexSequenceEvalSet Wise2_access_right_window_ComplexSequenceEvalSet
void Wise2_show_one_position_ComplexSequence(ComplexSequence * cs,int pos,FILE * ofp);
#define show_one_position_ComplexSequence Wise2_show_one_position_ComplexSequence
boolean Wise2_can_evaluate_this_type(ComplexSequenceEvalSet * cses,int type);
#define can_evaluate_this_type Wise2_can_evaluate_this_type
void Wise2_swap_ComplexSequenceEvalSet(ComplexSequenceEval ** list,int i,int j) ;
#define swap_ComplexSequenceEvalSet Wise2_swap_ComplexSequenceEvalSet
void Wise2_qsort_ComplexSequenceEvalSet(ComplexSequenceEval ** list,int left,int right,int (*comp)(ComplexSequenceEval * ,ComplexSequenceEval * ));
#define qsort_ComplexSequenceEvalSet Wise2_qsort_ComplexSequenceEvalSet
void Wise2_sort_ComplexSequenceEvalSet(ComplexSequenceEvalSet * obj,int (*comp)(ComplexSequenceEval *, ComplexSequenceEval *));
#define sort_ComplexSequenceEvalSet Wise2_sort_ComplexSequenceEvalSet
boolean Wise2_expand_ComplexSequenceEvalSet(ComplexSequenceEvalSet * obj,int len);
#define expand_ComplexSequenceEvalSet Wise2_expand_ComplexSequenceEvalSet

#ifdef _cplusplus
}
#endif

#endif
