#ifndef DYNAMITEdnaprofileengineHEADERFILE
#define DYNAMITEdnaprofileengineHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dnaprofile.h"
#include "dnaprofiledp.h"
#include "localcishit.h"

#include "pairwiseshortdna.h"

#define DnaProfileNode_LEAF 67
#define DnaProfileNode_SET  68
#define DnaProfileNode_UNCALC 69

#define DnaProfileSetLISTLENGTH 128

#define DnaProfileMatchPairSetLISTLENGTH 128

struct Wise2_DnaProfileEnginePara {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    DPRunImpl * dpri;    
    RandomModelDNA * rm;     
    LocalCisHitSetPara * setpara;    
    LocalCisHitScore * lchs;     
    Probability pseudo;  
    Probability open_unmatched;  
    Probability ext_unmatched;   
    Probability gap_unmatched;   
    Probability seq_id;  
    Probability m2i;     
    Probability m2d;     
    Probability i2i;     
    Probability d2d;     
    Score min_seq_prof;  
    } ;  
/* DnaProfileEnginePara defined */ 
#ifndef DYNAMITE_DEFINED_DnaProfileEnginePara
typedef struct Wise2_DnaProfileEnginePara Wise2_DnaProfileEnginePara;
#define DnaProfileEnginePara Wise2_DnaProfileEnginePara
#define DYNAMITE_DEFINED_DnaProfileEnginePara
#endif


struct Wise2_DnaProfileSet {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    DnaProfile ** dnap;  
    int len;/* len for above dnap  */ 
    int maxlen; /* maxlen for above dnap */ 
    } ;  
/* DnaProfileSet defined */ 
#ifndef DYNAMITE_DEFINED_DnaProfileSet
typedef struct Wise2_DnaProfileSet Wise2_DnaProfileSet;
#define DnaProfileSet Wise2_DnaProfileSet
#define DYNAMITE_DEFINED_DnaProfileSet
#endif


#ifndef DYNAMITE_DEFINED_DnaProfileNode
typedef struct Wise2_DnaProfileNode Wise2_DnaProfileNode;
#define DnaProfileNode Wise2_DnaProfileNode
#define DYNAMITE_DEFINED_DnaProfileNode
#endif

struct Wise2_DnaProfileNode {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int type;    
    Sequence * leaf;     
    DnaProfileSet * set;     
    DnaProfileNode * left;   
    DnaProfileNode * right;  
    DnaProfileNode * parent;     
    } ;  
/* DnaProfileNode defined */ 
#ifndef DYNAMITE_DEFINED_DnaProfileNode
typedef struct Wise2_DnaProfileNode Wise2_DnaProfileNode;
#define DnaProfileNode Wise2_DnaProfileNode
#define DYNAMITE_DEFINED_DnaProfileNode
#endif


struct Wise2_DnaProfileMatchPair {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    DnaProfile * query;  
    DnaProfile * target;     
    AlnBlock * alb;  
    int score;   
    int accepted;    
    } ;  
/* DnaProfileMatchPair defined */ 
#ifndef DYNAMITE_DEFINED_DnaProfileMatchPair
typedef struct Wise2_DnaProfileMatchPair Wise2_DnaProfileMatchPair;
#define DnaProfileMatchPair Wise2_DnaProfileMatchPair
#define DYNAMITE_DEFINED_DnaProfileMatchPair
#endif


struct Wise2_DnaProfileMatchPairSet {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    DnaProfileMatchPair ** pair;     
    int len;/* len for above pair  */ 
    int maxlen; /* maxlen for above pair */ 
    } ;  
/* DnaProfileMatchPairSet defined */ 
#ifndef DYNAMITE_DEFINED_DnaProfileMatchPairSet
typedef struct Wise2_DnaProfileMatchPairSet Wise2_DnaProfileMatchPairSet;
#define DnaProfileMatchPairSet Wise2_DnaProfileMatchPairSet
#define DYNAMITE_DEFINED_DnaProfileMatchPairSet
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_DnaProfileEnginePara(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DnaProfileEnginePara *]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileEnginePara *]
 *
 */
DnaProfileEnginePara * Wise2_hard_link_DnaProfileEnginePara(DnaProfileEnginePara * obj);
#define hard_link_DnaProfileEnginePara Wise2_hard_link_DnaProfileEnginePara


/* Function:  DnaProfileEnginePara_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileEnginePara *]
 *
 */
DnaProfileEnginePara * Wise2_DnaProfileEnginePara_alloc(void);
#define DnaProfileEnginePara_alloc Wise2_DnaProfileEnginePara_alloc


/* Function:  free_DnaProfileEnginePara(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DnaProfileEnginePara *]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileEnginePara *]
 *
 */
DnaProfileEnginePara * Wise2_free_DnaProfileEnginePara(DnaProfileEnginePara * obj);
#define free_DnaProfileEnginePara Wise2_free_DnaProfileEnginePara


/* Function:  add_DnaProfileSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [DnaProfileSet *]
 * Arg:        add [OWNER] Object to add to the list [DnaProfile *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_DnaProfileSet(DnaProfileSet * obj,DnaProfile * add);
#define add_DnaProfileSet Wise2_add_DnaProfileSet


/* Function:  flush_DnaProfileSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [DnaProfileSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_DnaProfileSet(DnaProfileSet * obj);
#define flush_DnaProfileSet Wise2_flush_DnaProfileSet


/* Function:  DnaProfileSet_alloc_std(void)
 *
 * Descrip:    Equivalent to DnaProfileSet_alloc_len(DnaProfileSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileSet *]
 *
 */
DnaProfileSet * Wise2_DnaProfileSet_alloc_std(void);
#define DnaProfileSet_alloc_std Wise2_DnaProfileSet_alloc_std


/* Function:  DnaProfileSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileSet *]
 *
 */
DnaProfileSet * Wise2_DnaProfileSet_alloc_len(int len);
#define DnaProfileSet_alloc_len Wise2_DnaProfileSet_alloc_len


/* Function:  hard_link_DnaProfileSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DnaProfileSet *]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileSet *]
 *
 */
DnaProfileSet * Wise2_hard_link_DnaProfileSet(DnaProfileSet * obj);
#define hard_link_DnaProfileSet Wise2_hard_link_DnaProfileSet


/* Function:  DnaProfileSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileSet *]
 *
 */
DnaProfileSet * Wise2_DnaProfileSet_alloc(void);
#define DnaProfileSet_alloc Wise2_DnaProfileSet_alloc


/* Function:  free_DnaProfileSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DnaProfileSet *]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileSet *]
 *
 */
DnaProfileSet * Wise2_free_DnaProfileSet(DnaProfileSet * obj);
#define free_DnaProfileSet Wise2_free_DnaProfileSet


/* Function:  hard_link_DnaProfileNode(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DnaProfileNode *]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileNode *]
 *
 */
DnaProfileNode * Wise2_hard_link_DnaProfileNode(DnaProfileNode * obj);
#define hard_link_DnaProfileNode Wise2_hard_link_DnaProfileNode


/* Function:  DnaProfileNode_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileNode *]
 *
 */
DnaProfileNode * Wise2_DnaProfileNode_alloc(void);
#define DnaProfileNode_alloc Wise2_DnaProfileNode_alloc


/* Function:  free_DnaProfileNode(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DnaProfileNode *]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileNode *]
 *
 */
DnaProfileNode * Wise2_free_DnaProfileNode(DnaProfileNode * obj);
#define free_DnaProfileNode Wise2_free_DnaProfileNode


/* Function:  hard_link_DnaProfileMatchPair(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DnaProfileMatchPair *]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileMatchPair *]
 *
 */
DnaProfileMatchPair * Wise2_hard_link_DnaProfileMatchPair(DnaProfileMatchPair * obj);
#define hard_link_DnaProfileMatchPair Wise2_hard_link_DnaProfileMatchPair


/* Function:  DnaProfileMatchPair_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileMatchPair *]
 *
 */
DnaProfileMatchPair * Wise2_DnaProfileMatchPair_alloc(void);
#define DnaProfileMatchPair_alloc Wise2_DnaProfileMatchPair_alloc


/* Function:  free_DnaProfileMatchPair(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DnaProfileMatchPair *]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileMatchPair *]
 *
 */
DnaProfileMatchPair * Wise2_free_DnaProfileMatchPair(DnaProfileMatchPair * obj);
#define free_DnaProfileMatchPair Wise2_free_DnaProfileMatchPair


/* Function:  add_DnaProfileMatchPairSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [DnaProfileMatchPairSet *]
 * Arg:        add [OWNER] Object to add to the list [DnaProfileMatchPair *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_DnaProfileMatchPairSet(DnaProfileMatchPairSet * obj,DnaProfileMatchPair * add);
#define add_DnaProfileMatchPairSet Wise2_add_DnaProfileMatchPairSet


/* Function:  flush_DnaProfileMatchPairSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [DnaProfileMatchPairSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_DnaProfileMatchPairSet(DnaProfileMatchPairSet * obj);
#define flush_DnaProfileMatchPairSet Wise2_flush_DnaProfileMatchPairSet


/* Function:  DnaProfileMatchPairSet_alloc_std(void)
 *
 * Descrip:    Equivalent to DnaProfileMatchPairSet_alloc_len(DnaProfileMatchPairSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileMatchPairSet *]
 *
 */
DnaProfileMatchPairSet * Wise2_DnaProfileMatchPairSet_alloc_std(void);
#define DnaProfileMatchPairSet_alloc_std Wise2_DnaProfileMatchPairSet_alloc_std


/* Function:  DnaProfileMatchPairSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileMatchPairSet *]
 *
 */
DnaProfileMatchPairSet * Wise2_DnaProfileMatchPairSet_alloc_len(int len);
#define DnaProfileMatchPairSet_alloc_len Wise2_DnaProfileMatchPairSet_alloc_len


/* Function:  hard_link_DnaProfileMatchPairSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DnaProfileMatchPairSet *]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileMatchPairSet *]
 *
 */
DnaProfileMatchPairSet * Wise2_hard_link_DnaProfileMatchPairSet(DnaProfileMatchPairSet * obj);
#define hard_link_DnaProfileMatchPairSet Wise2_hard_link_DnaProfileMatchPairSet


/* Function:  DnaProfileMatchPairSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileMatchPairSet *]
 *
 */
DnaProfileMatchPairSet * Wise2_DnaProfileMatchPairSet_alloc(void);
#define DnaProfileMatchPairSet_alloc Wise2_DnaProfileMatchPairSet_alloc


/* Function:  free_DnaProfileMatchPairSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DnaProfileMatchPairSet *]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileMatchPairSet *]
 *
 */
DnaProfileMatchPairSet * Wise2_free_DnaProfileMatchPairSet(DnaProfileMatchPairSet * obj);
#define free_DnaProfileMatchPairSet Wise2_free_DnaProfileMatchPairSet


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
TransFactorSet * Wise2_TransFactorSet_from_DnaProfileSet(DnaProfileSet * in);
#define TransFactorSet_from_DnaProfileSet Wise2_TransFactorSet_from_DnaProfileSet
DnaProfileNode * Wise2_balanced_4_Sequence_fasta_stream(FILE * ifp);
#define balanced_4_Sequence_fasta_stream Wise2_balanced_4_Sequence_fasta_stream
DnaProfileNode * Wise2_simple_cascade_Sequence_fasta_stream(FILE * ifp);
#define simple_cascade_Sequence_fasta_stream Wise2_simple_cascade_Sequence_fasta_stream
DnaProfileNode * Wise2_new_leaf_DnaProfileNode(Sequence * seq);
#define new_leaf_DnaProfileNode Wise2_new_leaf_DnaProfileNode
void Wise2_populate_DnaProfileNode_from_root(DnaProfileNode * root,DnaProfileEnginePara * dpep);
#define populate_DnaProfileNode_from_root Wise2_populate_DnaProfileNode_from_root
DnaProfileSet * Wise2_join_two_DnaProfileNode(DnaProfileNode * left,DnaProfileNode * right,DnaProfileEnginePara * dpep);
#define join_two_DnaProfileNode Wise2_join_two_DnaProfileNode
DnaProfileEnginePara * Wise2_new_DnaProfileEnginePara_from_argv(int * argc,char ** argv);
#define new_DnaProfileEnginePara_from_argv Wise2_new_DnaProfileEnginePara_from_argv
void Wise2_show_help_DnaProfileEnginePara(FILE * ofp);
#define show_help_DnaProfileEnginePara Wise2_show_help_DnaProfileEnginePara
DnaProfileSet * Wise2_filter_DnaProfileSet(DnaProfileSet * in,int min_length,int min_score);
#define filter_DnaProfileSet Wise2_filter_DnaProfileSet
DnaProfileSet * Wise2_DnaProfileSet_from_leaf_leaf(Sequence * one,Sequence * two,DnaProfileEnginePara * dpep);
#define DnaProfileSet_from_leaf_leaf Wise2_DnaProfileSet_from_leaf_leaf
DnaProfileSet * Wise2_DnaProfileSet_from_node_node(DnaProfileSet * one,DnaProfileSet * two,DnaProfileEnginePara * dpep);
#define DnaProfileSet_from_node_node Wise2_DnaProfileSet_from_node_node
DnaProfileSet * Wise2_DnaProfileSet_from_leaf_node(Sequence * one,DnaProfileSet * two,DnaProfileEnginePara * dpep);
#define DnaProfileSet_from_leaf_node Wise2_DnaProfileSet_from_leaf_node
DnaProfileMatchPair * Wise2_DnaProfileMatchPair_from_DnaProfile(DnaProfile * query,DnaProfile * target,DnaProfileEnginePara * dpep);
#define DnaProfileMatchPair_from_DnaProfile Wise2_DnaProfileMatchPair_from_DnaProfile
void Wise2_sort_DnaProfileMatchPairSet_by_score(DnaProfileMatchPairSet * set);
#define sort_DnaProfileMatchPairSet_by_score Wise2_sort_DnaProfileMatchPairSet_by_score
int Wise2_compare_DnaProfileMatchPair(DnaProfileMatchPair * one,DnaProfileMatchPair * two);
#define compare_DnaProfileMatchPair Wise2_compare_DnaProfileMatchPair
void Wise2_show_DnaProfileSet(DnaProfileSet * dnaps,RandomModelDNA * rm,FILE * ofp);
#define show_DnaProfileSet Wise2_show_DnaProfileSet


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_DnaProfileSet(DnaProfile ** list,int i,int j) ;
#define swap_DnaProfileSet Wise2_swap_DnaProfileSet
void Wise2_qsort_DnaProfileSet(DnaProfile ** list,int left,int right,int (*comp)(DnaProfile * ,DnaProfile * ));
#define qsort_DnaProfileSet Wise2_qsort_DnaProfileSet
void Wise2_sort_DnaProfileSet(DnaProfileSet * obj,int (*comp)(DnaProfile *, DnaProfile *));
#define sort_DnaProfileSet Wise2_sort_DnaProfileSet
boolean Wise2_expand_DnaProfileSet(DnaProfileSet * obj,int len);
#define expand_DnaProfileSet Wise2_expand_DnaProfileSet
void Wise2_swap_DnaProfileMatchPairSet(DnaProfileMatchPair ** list,int i,int j) ;
#define swap_DnaProfileMatchPairSet Wise2_swap_DnaProfileMatchPairSet
void Wise2_qsort_DnaProfileMatchPairSet(DnaProfileMatchPair ** list,int left,int right,int (*comp)(DnaProfileMatchPair * ,DnaProfileMatchPair * ));
#define qsort_DnaProfileMatchPairSet Wise2_qsort_DnaProfileMatchPairSet
void Wise2_sort_DnaProfileMatchPairSet(DnaProfileMatchPairSet * obj,int (*comp)(DnaProfileMatchPair *, DnaProfileMatchPair *));
#define sort_DnaProfileMatchPairSet Wise2_sort_DnaProfileMatchPairSet
boolean Wise2_expand_DnaProfileMatchPairSet(DnaProfileMatchPairSet * obj,int len);
#define expand_DnaProfileMatchPairSet Wise2_expand_DnaProfileMatchPairSet

#ifdef _cplusplus
}
#endif

#endif
