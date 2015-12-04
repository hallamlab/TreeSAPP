#ifndef DYNAMITEcomparapathHEADERFILE
#define DYNAMITEcomparapathHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "sequence.h"
#include "largeseqreader.h"
#include "singleseqspace.h"
#include "dnamapping.h"
#include "kmer_index_interface.h"
#include "kmer_direct.h"

#include "hsp.h"


#define COMPARALINK_START 8
#define COMPARALINK_LINEAR 256

#define COMPARA_NOTHING 0
#define COMPARA_QUERY_UNIQUE    1
#define COMPARA_QUERY_MULTIPLE  2
#define COMPARA_TARGET_UNIQUE   4
#define COMPARA_TARGET_MULTIPLE 8

#define COMPARA_SPLINE_FLIPPED  16
#define COMPARA_IS_REPEATED     32


#define COMPARA_IS_JOINT_FORWARD(a) is_joint_forward(a)
#define COMPARA_IS_JOINT_REVERSE(a) is_joint_reverse(a)

#define COMPARA_IS_JOINT_FORWARD_MACRO(a) (((a->state&COMPARA_QUERY_UNIQUE) && (a->state&COMPARA_TARGET_UNIQUE) && a->spline != NULL) ? 1 : 0)

#define COMPARA_IS_JOINT_REVERSE_MACRO(a) (((a->spline != NULL) && (((a->spline->state & COMPARA_QUERY_UNIQUE) && (a->state & COMPARA_TARGET_UNIQUE)) || ((a->spline->state & COMPARA_TARGET_UNIQUE) && (a->state & COMPARA_QUERY_UNIQUE)))) ? 1 : 0)

typedef struct ComparaHead {
  struct ComparaHead * next_query;
  long int position[2];
  struct ComparaHead * spline;
  long int number;
  char size;
  char state;
} ComparaHead;

typedef struct ComparaLinkStart {
  ComparaHead * start;
  Sequence * seq;
} ComparaLinkStart;

typedef struct ComparaHeadBlockAllocator {
  ComparaHead ** block;
  int current_block;
  int current_unit;
  int unit_length;
  int block_length;
} ComparaHeadBlockAllocator;

#define COMPARAHEAD_BA_BLOCK_LENGTH 10000
#define COMPARAHEAD_BA_UNIT_LENGTH  10000

typedef struct ComparaIndex {
  KmerIndexInterface * kii;
  ComparaLinkStart ** linkstart;
  int current_link;
  int link_len;
  SinglePosSpace * sps;
  ComparaHeadBlockAllocator * blockalloc;
} ComparaIndex;


#define COMPARAINDEX_LINK_START 16
#define COMPARAINDEX_LINK_LINEAR 256

#define ComparaLinkStartSetLISTLENGTH 256
#define SetofHSPsetLISTLENGTH 256

struct Wise2_ComparaLinkStartSet {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    ComparaLinkStart ** cls;     
    int len;/* len for above cls  */ 
    int maxlen; /* maxlen for above cls */ 
    } ;  
/* ComparaLinkStartSet defined */ 
#ifndef DYNAMITE_DEFINED_ComparaLinkStartSet
typedef struct Wise2_ComparaLinkStartSet Wise2_ComparaLinkStartSet;
#define ComparaLinkStartSet Wise2_ComparaLinkStartSet
#define DYNAMITE_DEFINED_ComparaLinkStartSet
#endif


struct Wise2_SetofHSPset {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    HSPset ** hspset;    
    int len;/* len for above hspset  */ 
    int maxlen; /* maxlen for above hspset */ 
    } ;  
/* SetofHSPset defined */ 
#ifndef DYNAMITE_DEFINED_SetofHSPset
typedef struct Wise2_SetofHSPset Wise2_SetofHSPset;
#define SetofHSPset Wise2_SetofHSPset
#define DYNAMITE_DEFINED_SetofHSPset
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  add_ComparaLinkStartSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [ComparaLinkStartSet *]
 * Arg:        add [OWNER] Object to add to the list [ComparaLinkStart *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_ComparaLinkStartSet(ComparaLinkStartSet * obj,ComparaLinkStart * add);
#define add_ComparaLinkStartSet Wise2_add_ComparaLinkStartSet


/* Function:  flush_ComparaLinkStartSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [ComparaLinkStartSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_ComparaLinkStartSet(ComparaLinkStartSet * obj);
#define flush_ComparaLinkStartSet Wise2_flush_ComparaLinkStartSet


/* Function:  ComparaLinkStartSet_alloc_std(void)
 *
 * Descrip:    Equivalent to ComparaLinkStartSet_alloc_len(ComparaLinkStartSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ComparaLinkStartSet *]
 *
 */
ComparaLinkStartSet * Wise2_ComparaLinkStartSet_alloc_std(void);
#define ComparaLinkStartSet_alloc_std Wise2_ComparaLinkStartSet_alloc_std


/* Function:  ComparaLinkStartSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [ComparaLinkStartSet *]
 *
 */
ComparaLinkStartSet * Wise2_ComparaLinkStartSet_alloc_len(int len);
#define ComparaLinkStartSet_alloc_len Wise2_ComparaLinkStartSet_alloc_len


/* Function:  hard_link_ComparaLinkStartSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ComparaLinkStartSet *]
 *
 * Return [UNKN ]  Undocumented return value [ComparaLinkStartSet *]
 *
 */
ComparaLinkStartSet * Wise2_hard_link_ComparaLinkStartSet(ComparaLinkStartSet * obj);
#define hard_link_ComparaLinkStartSet Wise2_hard_link_ComparaLinkStartSet


/* Function:  ComparaLinkStartSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ComparaLinkStartSet *]
 *
 */
ComparaLinkStartSet * Wise2_ComparaLinkStartSet_alloc(void);
#define ComparaLinkStartSet_alloc Wise2_ComparaLinkStartSet_alloc


/* Function:  free_ComparaLinkStartSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ComparaLinkStartSet *]
 *
 * Return [UNKN ]  Undocumented return value [ComparaLinkStartSet *]
 *
 */
ComparaLinkStartSet * Wise2_free_ComparaLinkStartSet(ComparaLinkStartSet * obj);
#define free_ComparaLinkStartSet Wise2_free_ComparaLinkStartSet


/* Function:  add_SetofHSPset(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [SetofHSPset *]
 * Arg:        add [OWNER] Object to add to the list [HSPset *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_SetofHSPset(SetofHSPset * obj,HSPset * add);
#define add_SetofHSPset Wise2_add_SetofHSPset


/* Function:  flush_SetofHSPset(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [SetofHSPset *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_SetofHSPset(SetofHSPset * obj);
#define flush_SetofHSPset Wise2_flush_SetofHSPset


/* Function:  SetofHSPset_alloc_std(void)
 *
 * Descrip:    Equivalent to SetofHSPset_alloc_len(SetofHSPsetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SetofHSPset *]
 *
 */
SetofHSPset * Wise2_SetofHSPset_alloc_std(void);
#define SetofHSPset_alloc_std Wise2_SetofHSPset_alloc_std


/* Function:  SetofHSPset_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [SetofHSPset *]
 *
 */
SetofHSPset * Wise2_SetofHSPset_alloc_len(int len);
#define SetofHSPset_alloc_len Wise2_SetofHSPset_alloc_len


/* Function:  hard_link_SetofHSPset(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SetofHSPset *]
 *
 * Return [UNKN ]  Undocumented return value [SetofHSPset *]
 *
 */
SetofHSPset * Wise2_hard_link_SetofHSPset(SetofHSPset * obj);
#define hard_link_SetofHSPset Wise2_hard_link_SetofHSPset


/* Function:  SetofHSPset_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SetofHSPset *]
 *
 */
SetofHSPset * Wise2_SetofHSPset_alloc(void);
#define SetofHSPset_alloc Wise2_SetofHSPset_alloc


/* Function:  free_SetofHSPset(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SetofHSPset *]
 *
 * Return [UNKN ]  Undocumented return value [SetofHSPset *]
 *
 */
SetofHSPset * Wise2_free_SetofHSPset(SetofHSPset * obj);
#define free_SetofHSPset Wise2_free_SetofHSPset


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
ComparaLinkStart * Wise2_free_ComparaLinkStart(ComparaLinkStart * cls);
#define free_ComparaLinkStart Wise2_free_ComparaLinkStart
void Wise2_show_SetofHSPset(SetofHSPset * set,FILE * ofp);
#define show_SetofHSPset Wise2_show_SetofHSPset
SetofHSPset * Wise2_SetofHSPset_from_ComparaIndex(ComparaIndex * ci,ComparaLinkStartSet * clss,FILE * logfp);
#define SetofHSPset_from_ComparaIndex Wise2_SetofHSPset_from_ComparaIndex
HSPset * Wise2_HSPset_from_ComparaIndex(ComparaIndex * ci,ComparaLinkStart * cls,FILE * logfp);
#define HSPset_from_ComparaIndex Wise2_HSPset_from_ComparaIndex
boolean Wise2_is_joint_forward(ComparaHead * h);
#define is_joint_forward Wise2_is_joint_forward
boolean Wise2_is_joint_reverse(ComparaHead * h);
#define is_joint_reverse Wise2_is_joint_reverse
void Wise2_show_distrib_ComparaIndex(ComparaIndex * ci,ComparaLinkStart * cls,FILE * ofp);
#define show_distrib_ComparaIndex Wise2_show_distrib_ComparaIndex
void Wise2_show_stats_ComparaIndex(ComparaIndex * ci,ComparaLinkStart * cls,FILE * ofp);
#define show_stats_ComparaIndex Wise2_show_stats_ComparaIndex
long int Wise2_insert_revcom_Splines_in_set(ComparaIndex * ci,ComparaLinkStartSet * clss,FILE * logfp);
#define insert_revcom_Splines_in_set Wise2_insert_revcom_Splines_in_set
long int Wise2_insert_revcom_Splines(ComparaIndex * ci,ComparaLinkStart * cls,FILE * logfp);
#define insert_revcom_Splines Wise2_insert_revcom_Splines
ComparaLinkStartSet * Wise2_add_Sequence_stream_ComparaIndex(ComparaIndex * ci,FILE * ifp,boolean is_target,int lognumber,int test_rev,FILE * logfp,char * tag);
#define add_Sequence_stream_ComparaIndex Wise2_add_Sequence_stream_ComparaIndex
ComparaLinkStart * Wise2_add_Sequence_ComparaIndex(ComparaIndex * ci,Sequence * seq,boolean is_target,int lognumber,long truncate,int skipsize,int * skipflag,int test_rev,FILE * logfp);
#define add_Sequence_ComparaIndex Wise2_add_Sequence_ComparaIndex
ComparaIndex * Wise2_new_ComparaIndex(KmerIndexInterface * kii);
#define new_ComparaIndex Wise2_new_ComparaIndex
ComparaLinkStart * Wise2_new_ComparaLinkStart(Sequence * seq);
#define new_ComparaLinkStart Wise2_new_ComparaLinkStart
void Wise2_add_ComparaLinkStart_to_ComparaIndex(ComparaIndex * ci,ComparaLinkStart * cls);
#define add_ComparaLinkStart_to_ComparaIndex Wise2_add_ComparaLinkStart_to_ComparaIndex
ComparaHead * Wise2_new_ComparaHead(ComparaHeadBlockAllocator * ba);
#define new_ComparaHead Wise2_new_ComparaHead
void Wise2_new_position_in_ComparaHead(ComparaHead * h,long int position);
#define new_position_in_ComparaHead Wise2_new_position_in_ComparaHead
ComparaHeadBlockAllocator * Wise2_new_ComparaHeadBlockAllocator(int block_length,int unit_length);
#define new_ComparaHeadBlockAllocator Wise2_new_ComparaHeadBlockAllocator
ComparaHead * Wise2_new_ComparaHead_from_ComparaHeadBlockAllocator(ComparaHeadBlockAllocator * chba);
#define new_ComparaHead_from_ComparaHeadBlockAllocator Wise2_new_ComparaHead_from_ComparaHeadBlockAllocator


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_ComparaLinkStartSet(ComparaLinkStart ** list,int i,int j) ;
#define swap_ComparaLinkStartSet Wise2_swap_ComparaLinkStartSet
void Wise2_qsort_ComparaLinkStartSet(ComparaLinkStart ** list,int left,int right,int (*comp)(ComparaLinkStart * ,ComparaLinkStart * ));
#define qsort_ComparaLinkStartSet Wise2_qsort_ComparaLinkStartSet
void Wise2_sort_ComparaLinkStartSet(ComparaLinkStartSet * obj,int (*comp)(ComparaLinkStart *, ComparaLinkStart *));
#define sort_ComparaLinkStartSet Wise2_sort_ComparaLinkStartSet
boolean Wise2_expand_ComparaLinkStartSet(ComparaLinkStartSet * obj,int len);
#define expand_ComparaLinkStartSet Wise2_expand_ComparaLinkStartSet
void Wise2_swap_SetofHSPset(HSPset ** list,int i,int j) ;
#define swap_SetofHSPset Wise2_swap_SetofHSPset
void Wise2_qsort_SetofHSPset(HSPset ** list,int left,int right,int (*comp)(HSPset * ,HSPset * ));
#define qsort_SetofHSPset Wise2_qsort_SetofHSPset
void Wise2_sort_SetofHSPset(SetofHSPset * obj,int (*comp)(HSPset *, HSPset *));
#define sort_SetofHSPset Wise2_sort_SetofHSPset
boolean Wise2_expand_SetofHSPset(SetofHSPset * obj,int len);
#define expand_SetofHSPset Wise2_expand_SetofHSPset

#ifdef _cplusplus
}
#endif

#endif
