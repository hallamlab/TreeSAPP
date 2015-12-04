#ifndef DYNAMITEeulerindexHEADERFILE
#define DYNAMITEeulerindexHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"
#include "singleseqspace.h"
#include "dnamapping.h"

#define EulerNodeLISTLENGTH 4
#define EulerGraphLISTLENGTH 128

#define EulerLinkLabelLength 8
#define EulerLinkLabelLinear 128

#ifndef DYNAMITE_DEFINED_EulerNode
typedef struct Wise2_EulerNode Wise2_EulerNode;
#define EulerNode Wise2_EulerNode
#define DYNAMITE_DEFINED_EulerNode
#endif

#ifndef DYNAMITE_DEFINED_EulerLink
typedef struct Wise2_EulerLink Wise2_EulerLink;
#define EulerLink Wise2_EulerLink
#define DYNAMITE_DEFINED_EulerLink
#endif

struct Wise2_EulerLink {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    EulerNode * prev;    
    EulerNode * next;    
    int * label;     
    int label_len;   
    int depth;   
    char base;   
    } ;  
/* EulerLink defined */ 
#ifndef DYNAMITE_DEFINED_EulerLink
typedef struct Wise2_EulerLink Wise2_EulerLink;
#define EulerLink Wise2_EulerLink
#define DYNAMITE_DEFINED_EulerLink
#endif


struct Wise2_EulerPath {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    EulerLink ** stack;  
    int max_stack_len;   
    int current;     
    } ;  
/* EulerPath defined */ 
#ifndef DYNAMITE_DEFINED_EulerPath
typedef struct Wise2_EulerPath Wise2_EulerPath;
#define EulerPath Wise2_EulerPath
#define DYNAMITE_DEFINED_EulerPath
#endif


struct Wise2_EulerNode {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int number;  
    EulerLink ** link;   
    int link_len;   /* len for above link  */ 
    int link_maxlen;/* maxlen for above link */ 
    EulerLink ** back;   
    int back_len;   /* len for above back  */ 
    int back_maxlen;/* maxlen for above back */ 
    char is_leftmost;    
    char is_rightmost;   
    } ;  
/* EulerNode defined */ 
#ifndef DYNAMITE_DEFINED_EulerNode
typedef struct Wise2_EulerNode Wise2_EulerNode;
#define EulerNode Wise2_EulerNode
#define DYNAMITE_DEFINED_EulerNode
#endif


struct Wise2_EulerGraph {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    EulerNode ** node;   
    int node_len;    
    int kmer;    
    Sequence ** seq;     
    int len;/* len for above seq  */ 
    int maxlen; /* maxlen for above seq */ 
    SinglePosSpace * sps;    
    EulerNode ** dup;    
    int dup_len;/* len for above dup  */ 
    int dup_maxlen; /* maxlen for above dup */ 
    } ;  
/* EulerGraph defined */ 
#ifndef DYNAMITE_DEFINED_EulerGraph
typedef struct Wise2_EulerGraph Wise2_EulerGraph;
#define EulerGraph Wise2_EulerGraph
#define DYNAMITE_DEFINED_EulerGraph
#endif


struct Wise2_EulerErrorPara {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int error_size;  
    } ;  
/* EulerErrorPara defined */ 
#ifndef DYNAMITE_DEFINED_EulerErrorPara
typedef struct Wise2_EulerErrorPara Wise2_EulerErrorPara;
#define EulerErrorPara Wise2_EulerErrorPara
#define DYNAMITE_DEFINED_EulerErrorPara
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_EulerLink(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [EulerLink *]
 *
 * Return [UNKN ]  Undocumented return value [EulerLink *]
 *
 */
EulerLink * Wise2_hard_link_EulerLink(EulerLink * obj);
#define hard_link_EulerLink Wise2_hard_link_EulerLink


/* Function:  EulerLink_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [EulerLink *]
 *
 */
EulerLink * Wise2_EulerLink_alloc(void);
#define EulerLink_alloc Wise2_EulerLink_alloc


/* Function:  free_EulerLink(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [EulerLink *]
 *
 * Return [UNKN ]  Undocumented return value [EulerLink *]
 *
 */
EulerLink * Wise2_free_EulerLink(EulerLink * obj);
#define free_EulerLink Wise2_free_EulerLink


/* Function:  hard_link_EulerPath(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [EulerPath *]
 *
 * Return [UNKN ]  Undocumented return value [EulerPath *]
 *
 */
EulerPath * Wise2_hard_link_EulerPath(EulerPath * obj);
#define hard_link_EulerPath Wise2_hard_link_EulerPath


/* Function:  EulerPath_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [EulerPath *]
 *
 */
EulerPath * Wise2_EulerPath_alloc(void);
#define EulerPath_alloc Wise2_EulerPath_alloc


/* Function:  free_EulerPath(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [EulerPath *]
 *
 * Return [UNKN ]  Undocumented return value [EulerPath *]
 *
 */
EulerPath * Wise2_free_EulerPath(EulerPath * obj);
#define free_EulerPath Wise2_free_EulerPath


/* Function:  add_link_EulerNode(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [EulerNode *]
 * Arg:        add [OWNER] Object to add to the list [EulerLink *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_link_EulerNode(EulerNode * obj,EulerLink * add);
#define add_link_EulerNode Wise2_add_link_EulerNode


/* Function:  flush_link_EulerNode(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [EulerNode *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_link_EulerNode(EulerNode * obj);
#define flush_link_EulerNode Wise2_flush_link_EulerNode


/* Function:  add_back_EulerNode(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [EulerNode *]
 * Arg:        add [OWNER] Object to add to the list [EulerLink *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_back_EulerNode(EulerNode * obj,EulerLink * add);
#define add_back_EulerNode Wise2_add_back_EulerNode


/* Function:  flush_back_EulerNode(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [EulerNode *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_back_EulerNode(EulerNode * obj);
#define flush_back_EulerNode Wise2_flush_back_EulerNode


/* Function:  EulerNode_alloc_std(void)
 *
 * Descrip:    Equivalent to EulerNode_alloc_len(EulerNodeLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [EulerNode *]
 *
 */
EulerNode * Wise2_EulerNode_alloc_std(void);
#define EulerNode_alloc_std Wise2_EulerNode_alloc_std


/* Function:  EulerNode_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [EulerNode *]
 *
 */
EulerNode * Wise2_EulerNode_alloc_len(int len);
#define EulerNode_alloc_len Wise2_EulerNode_alloc_len


/* Function:  hard_link_EulerNode(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [EulerNode *]
 *
 * Return [UNKN ]  Undocumented return value [EulerNode *]
 *
 */
EulerNode * Wise2_hard_link_EulerNode(EulerNode * obj);
#define hard_link_EulerNode Wise2_hard_link_EulerNode


/* Function:  EulerNode_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [EulerNode *]
 *
 */
EulerNode * Wise2_EulerNode_alloc(void);
#define EulerNode_alloc Wise2_EulerNode_alloc


/* Function:  free_EulerNode(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [EulerNode *]
 *
 * Return [UNKN ]  Undocumented return value [EulerNode *]
 *
 */
EulerNode * Wise2_free_EulerNode(EulerNode * obj);
#define free_EulerNode Wise2_free_EulerNode


/* Function:  add_EulerGraph(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [EulerGraph *]
 * Arg:        add [OWNER] Object to add to the list [Sequence *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_EulerGraph(EulerGraph * obj,Sequence * add);
#define add_EulerGraph Wise2_add_EulerGraph


/* Function:  flush_EulerGraph(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [EulerGraph *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_EulerGraph(EulerGraph * obj);
#define flush_EulerGraph Wise2_flush_EulerGraph


/* Function:  add_dup_EulerGraph(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [EulerGraph *]
 * Arg:        add [OWNER] Object to add to the list [EulerNode *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_dup_EulerGraph(EulerGraph * obj,EulerNode * add);
#define add_dup_EulerGraph Wise2_add_dup_EulerGraph


/* Function:  flush_dup_EulerGraph(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [EulerGraph *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_dup_EulerGraph(EulerGraph * obj);
#define flush_dup_EulerGraph Wise2_flush_dup_EulerGraph


/* Function:  EulerGraph_alloc_std(void)
 *
 * Descrip:    Equivalent to EulerGraph_alloc_len(EulerGraphLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [EulerGraph *]
 *
 */
EulerGraph * Wise2_EulerGraph_alloc_std(void);
#define EulerGraph_alloc_std Wise2_EulerGraph_alloc_std


/* Function:  EulerGraph_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [EulerGraph *]
 *
 */
EulerGraph * Wise2_EulerGraph_alloc_len(int len);
#define EulerGraph_alloc_len Wise2_EulerGraph_alloc_len


/* Function:  hard_link_EulerGraph(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [EulerGraph *]
 *
 * Return [UNKN ]  Undocumented return value [EulerGraph *]
 *
 */
EulerGraph * Wise2_hard_link_EulerGraph(EulerGraph * obj);
#define hard_link_EulerGraph Wise2_hard_link_EulerGraph


/* Function:  EulerGraph_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [EulerGraph *]
 *
 */
EulerGraph * Wise2_EulerGraph_alloc(void);
#define EulerGraph_alloc Wise2_EulerGraph_alloc


/* Function:  free_EulerGraph(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [EulerGraph *]
 *
 * Return [UNKN ]  Undocumented return value [EulerGraph *]
 *
 */
EulerGraph * Wise2_free_EulerGraph(EulerGraph * obj);
#define free_EulerGraph Wise2_free_EulerGraph


/* Function:  hard_link_EulerErrorPara(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [EulerErrorPara *]
 *
 * Return [UNKN ]  Undocumented return value [EulerErrorPara *]
 *
 */
EulerErrorPara * Wise2_hard_link_EulerErrorPara(EulerErrorPara * obj);
#define hard_link_EulerErrorPara Wise2_hard_link_EulerErrorPara


/* Function:  EulerErrorPara_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [EulerErrorPara *]
 *
 */
EulerErrorPara * Wise2_EulerErrorPara_alloc(void);
#define EulerErrorPara_alloc Wise2_EulerErrorPara_alloc


/* Function:  free_EulerErrorPara(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [EulerErrorPara *]
 *
 * Return [UNKN ]  Undocumented return value [EulerErrorPara *]
 *
 */
EulerErrorPara * Wise2_free_EulerErrorPara(EulerErrorPara * obj);
#define free_EulerErrorPara Wise2_free_EulerErrorPara


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
void Wise2_dump_EulerGraph(EulerGraph * eg,FILE * ofp);
#define dump_EulerGraph Wise2_dump_EulerGraph
Sequence * Wise2_read_Sequence_EulerNode(EulerNode * leftmost);
#define read_Sequence_EulerNode Wise2_read_Sequence_EulerNode
boolean Wise2_can_resolve_error_EulerGraph(EulerGraph *eg,EulerLink * leftmost);
#define can_resolve_error_EulerGraph Wise2_can_resolve_error_EulerGraph
boolean Wise2_resolve_error_EulerGraph(EulerGraph * eg,EulerLink * leftmost);
#define resolve_error_EulerGraph Wise2_resolve_error_EulerGraph
void Wise2_fix_error_EulerGraph(EulerGraph * eg, EulerLink * leftmost,char * dna_str,int len);
#define fix_error_EulerGraph Wise2_fix_error_EulerGraph
void Wise2_build_new_node_path_EulerGraph(EulerGraph * eg,EulerLink * leftmost,EulerPath * path,int * starting_labels,int length);
#define build_new_node_path_EulerGraph Wise2_build_new_node_path_EulerGraph
boolean Wise2_attempt_untangle_EulerPath(EulerGraph *eg,EulerPath * path,EulerLink * leftmost);
#define attempt_untangle_EulerPath Wise2_attempt_untangle_EulerPath
boolean Wise2_untangle_from_split_EulerNode(EulerGraph * eg,EulerNode * split_outgoing,int max_backtrack);
#define untangle_from_split_EulerNode Wise2_untangle_from_split_EulerNode
boolean Wise2_untangle_EulerLink_EulerPath(EulerGraph * eg,EulerPath * current_path,EulerLink * current,int * resolved,int backtrack_len,int max_backtrack);
#define untangle_EulerLink_EulerPath Wise2_untangle_EulerLink_EulerPath
boolean Wise2_remove_EulerLink_forward_EulerNode(EulerNode * n,EulerLink *l);
#define remove_EulerLink_forward_EulerNode Wise2_remove_EulerLink_forward_EulerNode
boolean Wise2_remove_EulerLink_backward_EulerNode(EulerNode * n,EulerLink *l);
#define remove_EulerLink_backward_EulerNode Wise2_remove_EulerLink_backward_EulerNode
boolean Wise2_remove_label_EulerLink(EulerLink * el,int label);
#define remove_label_EulerLink Wise2_remove_label_EulerLink
boolean Wise2_store_Sequence_EulerGraph(EulerGraph * eg,Sequence * seq);
#define store_Sequence_EulerGraph Wise2_store_Sequence_EulerGraph
boolean Wise2_store_EulerLink_EulerGraph(EulerGraph * eg,int prev_number,int next_number,char base,int label);
#define store_EulerLink_EulerGraph Wise2_store_EulerLink_EulerGraph
EulerLink * Wise2_new_EulerLink(void);
#define new_EulerLink Wise2_new_EulerLink
boolean Wise2_add_label_EulerLink(EulerLink * el,int label);
#define add_label_EulerLink Wise2_add_label_EulerLink
EulerNode * Wise2_new_EulerNode(int number);
#define new_EulerNode Wise2_new_EulerNode
EulerGraph * Wise2_new_EulerGraph(int kmer);
#define new_EulerGraph Wise2_new_EulerGraph
void * Wise2_push_EulerPath(EulerPath * ep,EulerLink * el);
#define push_EulerPath Wise2_push_EulerPath
EulerLink * Wise2_pop_EulerPath(EulerPath * ep);
#define pop_EulerPath Wise2_pop_EulerPath
void Wise2_extend_EulerPath_stack(EulerPath * ep);
#define extend_EulerPath_stack Wise2_extend_EulerPath_stack
EulerPath * Wise2_new_EulerPath(void);
#define new_EulerPath Wise2_new_EulerPath


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_link_EulerNode(EulerLink ** list,int i,int j) ;
#define swap_link_EulerNode Wise2_swap_link_EulerNode
void Wise2_qsort_link_EulerNode(EulerLink ** list,int left,int right,int (*comp)(EulerLink * ,EulerLink * ));
#define qsort_link_EulerNode Wise2_qsort_link_EulerNode
void Wise2_sort_link_EulerNode(EulerNode * obj,int (*comp)(EulerLink *, EulerLink *));
#define sort_link_EulerNode Wise2_sort_link_EulerNode
boolean Wise2_expand_link_EulerNode(EulerNode * obj,int len);
#define expand_link_EulerNode Wise2_expand_link_EulerNode
void Wise2_swap_back_EulerNode(EulerLink ** list,int i,int j) ;
#define swap_back_EulerNode Wise2_swap_back_EulerNode
void Wise2_qsort_back_EulerNode(EulerLink ** list,int left,int right,int (*comp)(EulerLink * ,EulerLink * ));
#define qsort_back_EulerNode Wise2_qsort_back_EulerNode
void Wise2_sort_back_EulerNode(EulerNode * obj,int (*comp)(EulerLink *, EulerLink *));
#define sort_back_EulerNode Wise2_sort_back_EulerNode
boolean Wise2_expand_back_EulerNode(EulerNode * obj,int len);
#define expand_back_EulerNode Wise2_expand_back_EulerNode
void Wise2_swap_EulerGraph(Sequence ** list,int i,int j) ;
#define swap_EulerGraph Wise2_swap_EulerGraph
void Wise2_qsort_EulerGraph(Sequence ** list,int left,int right,int (*comp)(Sequence * ,Sequence * ));
#define qsort_EulerGraph Wise2_qsort_EulerGraph
void Wise2_sort_EulerGraph(EulerGraph * obj,int (*comp)(Sequence *, Sequence *));
#define sort_EulerGraph Wise2_sort_EulerGraph
boolean Wise2_expand_EulerGraph(EulerGraph * obj,int len);
#define expand_EulerGraph Wise2_expand_EulerGraph
void Wise2_swap_dup_EulerGraph(EulerNode ** list,int i,int j) ;
#define swap_dup_EulerGraph Wise2_swap_dup_EulerGraph
void Wise2_qsort_dup_EulerGraph(EulerNode ** list,int left,int right,int (*comp)(EulerNode * ,EulerNode * ));
#define qsort_dup_EulerGraph Wise2_qsort_dup_EulerGraph
void Wise2_sort_dup_EulerGraph(EulerGraph * obj,int (*comp)(EulerNode *, EulerNode *));
#define sort_dup_EulerGraph Wise2_sort_dup_EulerGraph
boolean Wise2_expand_dup_EulerGraph(EulerGraph * obj,int len);
#define expand_dup_EulerGraph Wise2_expand_dup_EulerGraph

#ifdef _cplusplus
}
#endif

#endif
