#ifndef DYNAMITEkbestsearchHEADERFILE
#define DYNAMITEkbestsearchHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna2.h"
#include "dynafunc.h" 
#include "optimiser.h"

#define TransitionSetLISTLENGTH 24
#define Transition_from_special(trans) (trans->trans_type == TRANSITION_FROM_SPECIAL || trans->trans_type == TRANSITION_FROM_START ? 1 : 0)

#define TRANSITION_NORMAL 0
#define TRANSITION_FROM_SPECIAL 1
#define TRANSITION_FROM_START 2

/* Object Transition
 *
 * Descrip: This is the data structure for one
 *        transition - it will probably replace
 *        the generic matrix type data structure
 *        eventually, but as we need this for
 *        kbest type systems, we do it here
 *
 *
 */
struct Transition {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * to;   
    char * from;     
    int offi;    
    int offj;    
    char * calc;     
    int trans_type;  
    ExprTree * expr;     
    ExprTree * expr_state;   
    } ;  
/* Transition defined */ 
#ifndef DYNAMITE_DEFINED_Transition
typedef struct Transition Transition;
#define DYNAMITE_DEFINED_Transition
#endif


/* Object TransitionSet
 *
 * Descrip: Data for an entire DP matrix- just being the transitions
 *        at the moment
 *
 *        Sometime this - or something similar to this - will
 *        replace the generic matrix datastructure, but for the
 *        moment, this suffices!
 *
 *
 */
struct TransitionSet {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Transition ** trans;     
    int len;/* len for above trans  */ 
    int maxlen; /* maxlen for above trans */ 
    } ;  
/* TransitionSet defined */ 
#ifndef DYNAMITE_DEFINED_TransitionSet
typedef struct TransitionSet TransitionSet;
#define DYNAMITE_DEFINED_TransitionSet
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  write_kbest_score_GenericMatrix(dfp,gm,sc,mts,dpi)
 *
 * Descrip:    Produces a kbest search type single score function
 *
 *             kbest algorithm used here is that each cell is 
 *             reduced to a single score + state number that it
 *             came from. The kbest heuristic is to provide only k paths
 *             onto the next position in the sequence. In our case we
 *             have k = length of model / number of states. This sort
 *             of kbest heurisitc is good because it cuts down on excessive
 *             book keeping of the alignments, by being able to store
 *             information of the state only
 *
 *
 * Arg:        dfp [UNKN ] Undocumented argument [DYNFILE *]
 * Arg:         gm [UNKN ] Undocumented argument [GenericMatrix *]
 * Arg:         sc [UNKN ] Undocumented argument [Scope *]
 * Arg:        mts [UNKN ] Undocumented argument [MethodTypeSet *]
 * Arg:        dpi [UNKN ] Undocumented argument [DPImplementation *]
 *
 */
void write_kbest_score_GenericMatrix(DYNFILE * dfp,GenericMatrix * gm,Scope * sc,MethodTypeSet * mts,DPImplementation * dpi);


/* Function:  write_kbest_block(dfp,gm,matrixtag,pointertag,specialtag,use_special,cses,mts,dpi)
 *
 * Descrip:    Produces the actual kbest scoring inner loop
 *
 *
 * Arg:                dfp [UNKN ] Undocumented argument [DYNFILE *]
 * Arg:                 gm [UNKN ] Undocumented argument [GenericMatrix *]
 * Arg:          matrixtag [UNKN ] Undocumented argument [char *]
 * Arg:         pointertag [UNKN ] Undocumented argument [char *]
 * Arg:         specialtag [UNKN ] Undocumented argument [char *]
 * Arg:        use_special [UNKN ] Undocumented argument [boolean]
 * Arg:               cses [UNKN ] Undocumented argument [CommonSubExpressionSet *]
 * Arg:                mts [UNKN ] Undocumented argument [MethodTypeSet *]
 * Arg:                dpi [UNKN ] Undocumented argument [DPImplementation *]
 *
 */
void write_kbest_block(DYNFILE * dfp,GenericMatrix * gm,char * matrixtag,char * pointertag,char * specialtag,boolean use_special,CommonSubExpressionSet * cses,MethodTypeSet * mts,DPImplementation * dpi);


/* Function:  TransitionSet_from_GenericMatrix(gm)
 *
 * Descrip:    Makes a transition set from a generic matrix
 *
 *
 * Arg:        gm [UNKN ] Undocumented argument [GenericMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [TransitionSet *]
 *
 */
TransitionSet * TransitionSet_from_GenericMatrix(GenericMatrix * gm);


/* Function:  can_kbest_GenericMatrix(gm)
 *
 * Descrip:    sees whether the generic matrix is suitable for kbest optimisations
 *
 *
 * Arg:        gm [UNKN ] Undocumented argument [GenericMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean can_kbest_GenericMatrix(GenericMatrix * gm);


/* Function:  sort_TransitionSet_offset(ts)
 *
 * Descrip:    Sorts by offi then offj
 *
 *
 * Arg:        ts [UNKN ] Undocumented argument [TransitionSet *]
 *
 */
void sort_TransitionSet_offset(TransitionSet * ts);


/* Function:  comp_Transition(two,one)
 *
 * Descrip:    comparison by offi/offj
 *
 *
 * Arg:        two [UNKN ] Undocumented argument [Transition *]
 * Arg:        one [UNKN ] Undocumented argument [Transition *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int comp_Transition(Transition * two,Transition * one);


/* Function:  show_TransitionSet(tset,ofp)
 *
 * Descrip:    Shows a transition set
 *
 *
 * Arg:        tset [UNKN ] Undocumented argument [TransitionSet *]
 * Arg:         ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void show_TransitionSet(TransitionSet * tset,FILE * ofp);


/* Function:  hard_link_Transition(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Transition *]
 *
 * Return [UNKN ]  Undocumented return value [Transition *]
 *
 */
Transition * hard_link_Transition(Transition * obj);


/* Function:  Transition_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Transition *]
 *
 */
Transition * Transition_alloc(void);


/* Function:  free_Transition(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Transition *]
 *
 * Return [UNKN ]  Undocumented return value [Transition *]
 *
 */
Transition * free_Transition(Transition * obj);


/* Function:  add_TransitionSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [TransitionSet *]
 * Arg:        add [OWNER] Object to add to the list [Transition *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean add_TransitionSet(TransitionSet * obj,Transition * add);


/* Function:  flush_TransitionSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [TransitionSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_TransitionSet(TransitionSet * obj);


/* Function:  TransitionSet_alloc_std(void)
 *
 * Descrip:    Equivalent to TransitionSet_alloc_len(TransitionSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TransitionSet *]
 *
 */
TransitionSet * TransitionSet_alloc_std(void);


/* Function:  TransitionSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [TransitionSet *]
 *
 */
TransitionSet * TransitionSet_alloc_len(int len);


/* Function:  hard_link_TransitionSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [TransitionSet *]
 *
 * Return [UNKN ]  Undocumented return value [TransitionSet *]
 *
 */
TransitionSet * hard_link_TransitionSet(TransitionSet * obj);


/* Function:  TransitionSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TransitionSet *]
 *
 */
TransitionSet * TransitionSet_alloc(void);


/* Function:  free_TransitionSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [TransitionSet *]
 *
 * Return [UNKN ]  Undocumented return value [TransitionSet *]
 *
 */
TransitionSet * free_TransitionSet(TransitionSet * obj);


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void show_Transition(Transition * trans,FILE * ofp);
void swap_TransitionSet(Transition ** list,int i,int j) ;
void qsort_TransitionSet(Transition ** list,int left,int right,int (*comp)(Transition * ,Transition * ));
void sort_TransitionSet(TransitionSet * obj,int (*comp)(Transition *, Transition *));
boolean expand_TransitionSet(TransitionSet * obj,int len);

#ifdef _cplusplus
}
#endif

#endif
