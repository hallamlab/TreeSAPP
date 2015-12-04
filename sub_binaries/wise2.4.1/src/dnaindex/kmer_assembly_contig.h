#ifndef DYNAMITEkmer_assembly_contigHEADERFILE
#define DYNAMITEkmer_assembly_contigHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "kmer_assembly.h"
#include "assembly.h"

#define KmerAssemblyContigSetLISTLENGTH 256


#ifndef DYNAMITE_DEFINED_KmerAssemblyContig
typedef struct Wise2_KmerAssemblyContig Wise2_KmerAssemblyContig;
#define KmerAssemblyContig Wise2_KmerAssemblyContig
#define DYNAMITE_DEFINED_KmerAssemblyContig
#endif

struct Wise2_KmerAssemblyContig {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    KmerAssemblyNode * start;    
    KmerAssemblyNode * end;  
    int len;     
    boolean clean_start;     
    boolean clean_end;   
    int max_depth;   
    KmerAssemblyContig * mirror;     
    int is_mirror;   
    } ;  
/* KmerAssemblyContig defined */ 
#ifndef DYNAMITE_DEFINED_KmerAssemblyContig
typedef struct Wise2_KmerAssemblyContig Wise2_KmerAssemblyContig;
#define KmerAssemblyContig Wise2_KmerAssemblyContig
#define DYNAMITE_DEFINED_KmerAssemblyContig
#endif


struct Wise2_KmerAssemblyContigSet {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    KmerAssemblyContig ** contig;    
    int len;/* len for above contig  */ 
    int maxlen; /* maxlen for above contig */ 
    KmerAssemblyIndex * kai;     
    } ;  
/* KmerAssemblyContigSet defined */ 
#ifndef DYNAMITE_DEFINED_KmerAssemblyContigSet
typedef struct Wise2_KmerAssemblyContigSet Wise2_KmerAssemblyContigSet;
#define KmerAssemblyContigSet Wise2_KmerAssemblyContigSet
#define DYNAMITE_DEFINED_KmerAssemblyContigSet
#endif


struct Wise2_KmerAssemblyContigPara {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int minimum_len;     
    int minimum_depth;   
    } ;  
/* KmerAssemblyContigPara defined */ 
#ifndef DYNAMITE_DEFINED_KmerAssemblyContigPara
typedef struct Wise2_KmerAssemblyContigPara Wise2_KmerAssemblyContigPara;
#define KmerAssemblyContigPara Wise2_KmerAssemblyContigPara
#define DYNAMITE_DEFINED_KmerAssemblyContigPara
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_KmerAssemblyContig(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [KmerAssemblyContig *]
 *
 * Return [UNKN ]  Undocumented return value [KmerAssemblyContig *]
 *
 */
KmerAssemblyContig * Wise2_hard_link_KmerAssemblyContig(KmerAssemblyContig * obj);
#define hard_link_KmerAssemblyContig Wise2_hard_link_KmerAssemblyContig


/* Function:  KmerAssemblyContig_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [KmerAssemblyContig *]
 *
 */
KmerAssemblyContig * Wise2_KmerAssemblyContig_alloc(void);
#define KmerAssemblyContig_alloc Wise2_KmerAssemblyContig_alloc


/* Function:  free_KmerAssemblyContig(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [KmerAssemblyContig *]
 *
 * Return [UNKN ]  Undocumented return value [KmerAssemblyContig *]
 *
 */
KmerAssemblyContig * Wise2_free_KmerAssemblyContig(KmerAssemblyContig * obj);
#define free_KmerAssemblyContig Wise2_free_KmerAssemblyContig


/* Function:  add_KmerAssemblyContigSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [KmerAssemblyContigSet *]
 * Arg:        add [OWNER] Object to add to the list [KmerAssemblyContig *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_KmerAssemblyContigSet(KmerAssemblyContigSet * obj,KmerAssemblyContig * add);
#define add_KmerAssemblyContigSet Wise2_add_KmerAssemblyContigSet


/* Function:  flush_KmerAssemblyContigSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [KmerAssemblyContigSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_KmerAssemblyContigSet(KmerAssemblyContigSet * obj);
#define flush_KmerAssemblyContigSet Wise2_flush_KmerAssemblyContigSet


/* Function:  KmerAssemblyContigSet_alloc_std(void)
 *
 * Descrip:    Equivalent to KmerAssemblyContigSet_alloc_len(KmerAssemblyContigSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [KmerAssemblyContigSet *]
 *
 */
KmerAssemblyContigSet * Wise2_KmerAssemblyContigSet_alloc_std(void);
#define KmerAssemblyContigSet_alloc_std Wise2_KmerAssemblyContigSet_alloc_std


/* Function:  KmerAssemblyContigSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [KmerAssemblyContigSet *]
 *
 */
KmerAssemblyContigSet * Wise2_KmerAssemblyContigSet_alloc_len(int len);
#define KmerAssemblyContigSet_alloc_len Wise2_KmerAssemblyContigSet_alloc_len


/* Function:  hard_link_KmerAssemblyContigSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [KmerAssemblyContigSet *]
 *
 * Return [UNKN ]  Undocumented return value [KmerAssemblyContigSet *]
 *
 */
KmerAssemblyContigSet * Wise2_hard_link_KmerAssemblyContigSet(KmerAssemblyContigSet * obj);
#define hard_link_KmerAssemblyContigSet Wise2_hard_link_KmerAssemblyContigSet


/* Function:  KmerAssemblyContigSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [KmerAssemblyContigSet *]
 *
 */
KmerAssemblyContigSet * Wise2_KmerAssemblyContigSet_alloc(void);
#define KmerAssemblyContigSet_alloc Wise2_KmerAssemblyContigSet_alloc


/* Function:  free_KmerAssemblyContigSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [KmerAssemblyContigSet *]
 *
 * Return [UNKN ]  Undocumented return value [KmerAssemblyContigSet *]
 *
 */
KmerAssemblyContigSet * Wise2_free_KmerAssemblyContigSet(KmerAssemblyContigSet * obj);
#define free_KmerAssemblyContigSet Wise2_free_KmerAssemblyContigSet


/* Function:  hard_link_KmerAssemblyContigPara(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [KmerAssemblyContigPara *]
 *
 * Return [UNKN ]  Undocumented return value [KmerAssemblyContigPara *]
 *
 */
KmerAssemblyContigPara * Wise2_hard_link_KmerAssemblyContigPara(KmerAssemblyContigPara * obj);
#define hard_link_KmerAssemblyContigPara Wise2_hard_link_KmerAssemblyContigPara


/* Function:  KmerAssemblyContigPara_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [KmerAssemblyContigPara *]
 *
 */
KmerAssemblyContigPara * Wise2_KmerAssemblyContigPara_alloc(void);
#define KmerAssemblyContigPara_alloc Wise2_KmerAssemblyContigPara_alloc


/* Function:  free_KmerAssemblyContigPara(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [KmerAssemblyContigPara *]
 *
 * Return [UNKN ]  Undocumented return value [KmerAssemblyContigPara *]
 *
 */
KmerAssemblyContigPara * Wise2_free_KmerAssemblyContigPara(KmerAssemblyContigPara * obj);
#define free_KmerAssemblyContigPara Wise2_free_KmerAssemblyContigPara


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
KmerAssemblyContig * Wise2_find_strict_mirrored_KmerAssemblyContig(KmerAssemblyContigSet * kcs,KmerAssemblyContig * c);
#define find_strict_mirrored_KmerAssemblyContig Wise2_find_strict_mirrored_KmerAssemblyContig
KmerAssemblyContigSet * Wise2_KmerAssemblyContigSet_from_KmerAssemblyIndex(KmerAssemblyIndex * kai);
#define KmerAssemblyContigSet_from_KmerAssemblyIndex Wise2_KmerAssemblyContigSet_from_KmerAssemblyIndex
KmerAssemblyContig * Wise2_new_KmerAssemblyContig(KmerAssemblyNode * node);
#define new_KmerAssemblyContig Wise2_new_KmerAssemblyContig
Assembly * Wise2_Assembly_from_KmerAssemblyIndex(KmerAssemblyIndex * kai,KmerAssemblyContigPara * p);
#define Assembly_from_KmerAssemblyIndex Wise2_Assembly_from_KmerAssemblyIndex
Assembly * Wise2_Assembly_from_KmerAssemblyContigSet(KmerAssemblyContigSet * kacs,KmerAssemblyContigPara * p);
#define Assembly_from_KmerAssemblyContigSet Wise2_Assembly_from_KmerAssemblyContigSet
AssemblyContig * Wise2_AssemblyContig_from_KmerAssemblyContig(KmerAssemblyContig * kac);
#define AssemblyContig_from_KmerAssemblyContig Wise2_AssemblyContig_from_KmerAssemblyContig


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_KmerAssemblyContigSet(KmerAssemblyContig ** list,int i,int j) ;
#define swap_KmerAssemblyContigSet Wise2_swap_KmerAssemblyContigSet
void Wise2_qsort_KmerAssemblyContigSet(KmerAssemblyContig ** list,int left,int right,int (*comp)(KmerAssemblyContig * ,KmerAssemblyContig * ));
#define qsort_KmerAssemblyContigSet Wise2_qsort_KmerAssemblyContigSet
void Wise2_sort_KmerAssemblyContigSet(KmerAssemblyContigSet * obj,int (*comp)(KmerAssemblyContig *, KmerAssemblyContig *));
#define sort_KmerAssemblyContigSet Wise2_sort_KmerAssemblyContigSet
boolean Wise2_expand_KmerAssemblyContigSet(KmerAssemblyContigSet * obj,int len);
#define expand_KmerAssemblyContigSet Wise2_expand_KmerAssemblyContigSet

#ifdef _cplusplus
}
#endif

#endif
