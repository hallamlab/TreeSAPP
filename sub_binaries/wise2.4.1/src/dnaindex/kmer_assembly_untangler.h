#ifndef DYNAMITEkmer_assembly_untanglerHEADERFILE
#define DYNAMITEkmer_assembly_untanglerHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "kmer_assembly.h"


#define KmerAssemblyPath_MAX_STACK 1024

#define MAX_TANGLE_DEPTH 10000


struct Wise2_KmerAssemblyPath {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    KmerAssemblyLink ** stack;   
    int stack_len;   
    int max_stack;   
    boolean right_tail;  
    boolean left_tail;   
    } ;  
/* KmerAssemblyPath defined */ 
#ifndef DYNAMITE_DEFINED_KmerAssemblyPath
typedef struct Wise2_KmerAssemblyPath Wise2_KmerAssemblyPath;
#define KmerAssemblyPath Wise2_KmerAssemblyPath
#define DYNAMITE_DEFINED_KmerAssemblyPath
#endif


struct Wise2_KmerTangleResolver {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    long start_label[MAX_TANGLE_DEPTH];  
    SinglePosSequence * start_pos[MAX_TANGLE_DEPTH];     
    int distance_to_end[MAX_TANGLE_DEPTH];   
    int max_end;     
    int label_len;   
    int max_path;    
    } ;  
/* KmerTangleResolver defined */ 
#ifndef DYNAMITE_DEFINED_KmerTangleResolver
typedef struct Wise2_KmerTangleResolver Wise2_KmerTangleResolver;
#define KmerTangleResolver Wise2_KmerTangleResolver
#define DYNAMITE_DEFINED_KmerTangleResolver
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_KmerAssemblyPath(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [KmerAssemblyPath *]
 *
 * Return [UNKN ]  Undocumented return value [KmerAssemblyPath *]
 *
 */
KmerAssemblyPath * Wise2_hard_link_KmerAssemblyPath(KmerAssemblyPath * obj);
#define hard_link_KmerAssemblyPath Wise2_hard_link_KmerAssemblyPath


/* Function:  KmerAssemblyPath_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [KmerAssemblyPath *]
 *
 */
KmerAssemblyPath * Wise2_KmerAssemblyPath_alloc(void);
#define KmerAssemblyPath_alloc Wise2_KmerAssemblyPath_alloc


/* Function:  free_KmerAssemblyPath(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [KmerAssemblyPath *]
 *
 * Return [UNKN ]  Undocumented return value [KmerAssemblyPath *]
 *
 */
KmerAssemblyPath * Wise2_free_KmerAssemblyPath(KmerAssemblyPath * obj);
#define free_KmerAssemblyPath Wise2_free_KmerAssemblyPath


/* Function:  hard_link_KmerTangleResolver(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [KmerTangleResolver *]
 *
 * Return [UNKN ]  Undocumented return value [KmerTangleResolver *]
 *
 */
KmerTangleResolver * Wise2_hard_link_KmerTangleResolver(KmerTangleResolver * obj);
#define hard_link_KmerTangleResolver Wise2_hard_link_KmerTangleResolver


/* Function:  KmerTangleResolver_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [KmerTangleResolver *]
 *
 */
KmerTangleResolver * Wise2_KmerTangleResolver_alloc(void);
#define KmerTangleResolver_alloc Wise2_KmerTangleResolver_alloc


/* Function:  free_KmerTangleResolver(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [KmerTangleResolver *]
 *
 * Return [UNKN ]  Undocumented return value [KmerTangleResolver *]
 *
 */
KmerTangleResolver * Wise2_free_KmerTangleResolver(KmerTangleResolver * obj);
#define free_KmerTangleResolver Wise2_free_KmerTangleResolver


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
void Wise2_show_KmerAssemblyPath(KmerAssemblyIndex * kai,KmerAssemblyPath * kap,FILE * ofp);
#define show_KmerAssemblyPath Wise2_show_KmerAssemblyPath
void Wise2_push_KmerAssemblyPath(KmerAssemblyPath * kap,KmerAssemblyLink * kal);
#define push_KmerAssemblyPath Wise2_push_KmerAssemblyPath
KmerAssemblyLink * Wise2_pop_KmerAssemblyPath(KmerAssemblyPath * kap);
#define pop_KmerAssemblyPath Wise2_pop_KmerAssemblyPath
KmerAssemblyPath * Wise2_new_KmerAssemblyPath(void);
#define new_KmerAssemblyPath Wise2_new_KmerAssemblyPath
int Wise2_untangle_KmerAssembly(KmerAssemblyIndex * kai);
#define untangle_KmerAssembly Wise2_untangle_KmerAssembly
boolean Wise2_recursive_untangle_KmerAssembly(KmerTangleResolver * res,KmerAssemblyIndex * kai,KmerAssemblyLink * current,KmerAssemblyPath * p,int current_pos);
#define recursive_untangle_KmerAssembly Wise2_recursive_untangle_KmerAssembly
KmerTangleResolver * Wise2_new_KmerTangleResolver(KmerAssemblyIndex * kai,KmerAssemblyLink * leftmost,int max_path);
#define new_KmerTangleResolver Wise2_new_KmerTangleResolver
boolean Wise2_no_longer_active_KmerTangleResolver(KmerTangleResolver * res,KmerAssemblyLink * current,int pathlen);
#define no_longer_active_KmerTangleResolver Wise2_no_longer_active_KmerTangleResolver
boolean Wise2_is_righthand_link_KmerTangleResolver(KmerTangleResolver * res,KmerAssemblyLink * rightside,int pathlen);
#define is_righthand_link_KmerTangleResolver Wise2_is_righthand_link_KmerTangleResolver
boolean Wise2_attempt_forward_untangle_KmerAssembly(KmerAssemblyIndex * kai,KmerAssemblyLink * left_input,int max_search);
#define attempt_forward_untangle_KmerAssembly Wise2_attempt_forward_untangle_KmerAssembly
boolean Wise2_old_attempt_forward_untangle_KmerAssembly(KmerAssemblyIndex * kai,KmerAssemblyLink * left_input,int max_search);
#define old_attempt_forward_untangle_KmerAssembly Wise2_old_attempt_forward_untangle_KmerAssembly
void Wise2_lift_forward_tangled_tail(KmerAssemblyIndex * kai,KmerAssemblyPath * tail,long int * start_label,int label_len);
#define lift_forward_tangled_tail Wise2_lift_forward_tangled_tail
void Wise2_lift_backward_tangled_tail(KmerAssemblyIndex * kai,KmerAssemblyLink * new,KmerAssemblyPath * tail,int * start_label,SinglePosSequence ** positions,int label_len);
#define lift_backward_tangled_tail Wise2_lift_backward_tangled_tail
KmerAssemblyPath * Wise2_lift_forward_tangled_KmerAssemblyPath(KmerAssemblyIndex * kai,KmerAssemblyPath * kap,long int * start_label,int label_len);
#define lift_forward_tangled_KmerAssemblyPath Wise2_lift_forward_tangled_KmerAssemblyPath


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/

#ifdef _cplusplus
}
#endif

#endif
