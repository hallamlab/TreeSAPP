#ifndef DYNAMITEgwlitemodelHEADERFILE
#define DYNAMITEgwlitemodelHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"
#include "genewisemodel.h"

#define GwLiteScoreLISTLENGTH 128
#define GwLiteLISTLENGTH 128


enum GwLiteTransition {
  GWL_MATCH2MATCH,
  GWL_MATCH2INSERT,
  GWL_MATCH2DELETE,
  GWL_MATCH2END,
  GWL_INSERT2MATCH,
  GWL_INSERT2INSERT,
  GWL_INSERT2DELETE,
  GWL_INSERT2END,
  GWL_DELETE2MATCH,
  GWL_DELETE2INSERT,
  GWL_DELETE2DELETE,
  GWL_DELETE2END,
  GWL_START2MATCH,
  GWL_START2INSERT,
  GWL_START2DELETE,
  GWL_TRANSITION_LEN
};

#define GWL_EMISSION_LEN 65

/* Object GwLiteSegment
 *
 * Descrip: This is a particular HMM node, with
 *        match and insert emissions in the codon space
 *        and the transitions for the genewise light model
 *
 *
 */
struct Wise2_GwLiteSegment {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Probability match[GWL_EMISSION_LEN];     
    Probability insert[GWL_EMISSION_LEN];    
    Probability transition[GWL_TRANSITION_LEN];  
    } ;  
/* GwLiteSegment defined */ 
#ifndef DYNAMITE_DEFINED_GwLiteSegment
typedef struct Wise2_GwLiteSegment Wise2_GwLiteSegment;
#define GwLiteSegment Wise2_GwLiteSegment
#define DYNAMITE_DEFINED_GwLiteSegment
#endif


/* Object GwLite
 *
 * Descrip: This is the lightweight version 
 *        of GeneWise designed for faster executation
 *        and better portability to hardware environments
 *
 *
 */
struct Wise2_GwLite {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    GwLiteSegment ** seg;    
    int len;/* len for above seg  */ 
    int maxlen; /* maxlen for above seg */ 
    char * name;     
    } ;  
/* GwLite defined */ 
#ifndef DYNAMITE_DEFINED_GwLite
typedef struct Wise2_GwLite Wise2_GwLite;
#define GwLite Wise2_GwLite
#define DYNAMITE_DEFINED_GwLite
#endif


/* Object GwLiteSegmentScore
 *
 * Descrip: This is the genewise lite scoring
 *        data structure
 *
 *
 */
struct Wise2_GwLiteSegmentScore {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Score match[GWL_EMISSION_LEN];   
    Score insert[GWL_EMISSION_LEN];  
    Score transition[GWL_TRANSITION_LEN];    
    } ;  
/* GwLiteSegmentScore defined */ 
#ifndef DYNAMITE_DEFINED_GwLiteSegmentScore
typedef struct Wise2_GwLiteSegmentScore Wise2_GwLiteSegmentScore;
#define GwLiteSegmentScore Wise2_GwLiteSegmentScore
#define DYNAMITE_DEFINED_GwLiteSegmentScore
#endif


/* Object GwLiteScore
 *
 * Descrip: This is the lightweight version 
 *        of GeneWise designed for faster executation
 *        and better portability to hardware environments
 *
 *
 */
struct Wise2_GwLiteScore {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    GwLiteSegmentScore ** seg;   
    int len;/* len for above seg  */ 
    int maxlen; /* maxlen for above seg */ 
    char * name;     
    } ;  
/* GwLiteScore defined */ 
#ifndef DYNAMITE_DEFINED_GwLiteScore
typedef struct Wise2_GwLiteScore Wise2_GwLiteScore;
#define GwLiteScore Wise2_GwLiteScore
#define DYNAMITE_DEFINED_GwLiteScore
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  GwLite_AlnBlock_surgery(alb)
 *
 * Descrip:    A pretty weird function. Takes an AlnBlock made by GwLite and
 *             performs the necessary surgery at the 3SS to make it look like
 *             the AlnBlocks produced by the other genewise models. This means
 *             it has to eat into the intron by 3 residues 
 *
 *
 * Arg:        alb [UNKN ] Undocumented argument [AlnBlock *]
 *
 */
void Wise2_GwLite_AlnBlock_surgery(AlnBlock * alb);
#define GwLite_AlnBlock_surgery Wise2_GwLite_AlnBlock_surgery


/* Function:  GwLite_from_GeneWise(gwm)
 *
 * Descrip:    Builds a GwLite model from the GeneWise model
 *
 *
 * Arg:        gwm [UNKN ] Undocumented argument [GeneWise *]
 *
 * Return [UNKN ]  Undocumented return value [GwLite *]
 *
 */
GwLite * Wise2_GwLite_from_GeneWise(GeneWise * gwm);
#define GwLite_from_GeneWise Wise2_GwLite_from_GeneWise


/* Function:  GwLiteScore_from_GwLite(gwl)
 *
 * Descrip:    Makes a lite score from a lite probability basis
 *
 *
 * Arg:        gwl [UNKN ] Undocumented argument [GwLite *]
 *
 * Return [UNKN ]  Undocumented return value [GwLiteScore *]
 *
 */
GwLiteScore * Wise2_GwLiteScore_from_GwLite(GwLite * gwl);
#define GwLiteScore_from_GwLite Wise2_GwLiteScore_from_GwLite


/* Function:  hard_link_GwLiteSegment(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GwLiteSegment *]
 *
 * Return [UNKN ]  Undocumented return value [GwLiteSegment *]
 *
 */
GwLiteSegment * Wise2_hard_link_GwLiteSegment(GwLiteSegment * obj);
#define hard_link_GwLiteSegment Wise2_hard_link_GwLiteSegment


/* Function:  GwLiteSegment_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GwLiteSegment *]
 *
 */
GwLiteSegment * Wise2_GwLiteSegment_alloc(void);
#define GwLiteSegment_alloc Wise2_GwLiteSegment_alloc


/* Function:  free_GwLiteSegment(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GwLiteSegment *]
 *
 * Return [UNKN ]  Undocumented return value [GwLiteSegment *]
 *
 */
GwLiteSegment * Wise2_free_GwLiteSegment(GwLiteSegment * obj);
#define free_GwLiteSegment Wise2_free_GwLiteSegment


/* Function:  add_GwLite(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GwLite *]
 * Arg:        add [OWNER] Object to add to the list [GwLiteSegment *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_GwLite(GwLite * obj,GwLiteSegment * add);
#define add_GwLite Wise2_add_GwLite


/* Function:  flush_GwLite(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [GwLite *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_GwLite(GwLite * obj);
#define flush_GwLite Wise2_flush_GwLite


/* Function:  GwLite_alloc_std(void)
 *
 * Descrip:    Equivalent to GwLite_alloc_len(GwLiteLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GwLite *]
 *
 */
GwLite * Wise2_GwLite_alloc_std(void);
#define GwLite_alloc_std Wise2_GwLite_alloc_std


/* Function:  GwLite_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [GwLite *]
 *
 */
GwLite * Wise2_GwLite_alloc_len(int len);
#define GwLite_alloc_len Wise2_GwLite_alloc_len


/* Function:  hard_link_GwLite(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GwLite *]
 *
 * Return [UNKN ]  Undocumented return value [GwLite *]
 *
 */
GwLite * Wise2_hard_link_GwLite(GwLite * obj);
#define hard_link_GwLite Wise2_hard_link_GwLite


/* Function:  GwLite_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GwLite *]
 *
 */
GwLite * Wise2_GwLite_alloc(void);
#define GwLite_alloc Wise2_GwLite_alloc


/* Function:  free_GwLite(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GwLite *]
 *
 * Return [UNKN ]  Undocumented return value [GwLite *]
 *
 */
GwLite * Wise2_free_GwLite(GwLite * obj);
#define free_GwLite Wise2_free_GwLite


/* Function:  hard_link_GwLiteSegmentScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GwLiteSegmentScore *]
 *
 * Return [UNKN ]  Undocumented return value [GwLiteSegmentScore *]
 *
 */
GwLiteSegmentScore * Wise2_hard_link_GwLiteSegmentScore(GwLiteSegmentScore * obj);
#define hard_link_GwLiteSegmentScore Wise2_hard_link_GwLiteSegmentScore


/* Function:  GwLiteSegmentScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GwLiteSegmentScore *]
 *
 */
GwLiteSegmentScore * Wise2_GwLiteSegmentScore_alloc(void);
#define GwLiteSegmentScore_alloc Wise2_GwLiteSegmentScore_alloc


/* Function:  free_GwLiteSegmentScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GwLiteSegmentScore *]
 *
 * Return [UNKN ]  Undocumented return value [GwLiteSegmentScore *]
 *
 */
GwLiteSegmentScore * Wise2_free_GwLiteSegmentScore(GwLiteSegmentScore * obj);
#define free_GwLiteSegmentScore Wise2_free_GwLiteSegmentScore


/* Function:  add_GwLiteScore(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GwLiteScore *]
 * Arg:        add [OWNER] Object to add to the list [GwLiteSegmentScore *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_GwLiteScore(GwLiteScore * obj,GwLiteSegmentScore * add);
#define add_GwLiteScore Wise2_add_GwLiteScore


/* Function:  flush_GwLiteScore(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [GwLiteScore *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_GwLiteScore(GwLiteScore * obj);
#define flush_GwLiteScore Wise2_flush_GwLiteScore


/* Function:  GwLiteScore_alloc_std(void)
 *
 * Descrip:    Equivalent to GwLiteScore_alloc_len(GwLiteScoreLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GwLiteScore *]
 *
 */
GwLiteScore * Wise2_GwLiteScore_alloc_std(void);
#define GwLiteScore_alloc_std Wise2_GwLiteScore_alloc_std


/* Function:  GwLiteScore_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [GwLiteScore *]
 *
 */
GwLiteScore * Wise2_GwLiteScore_alloc_len(int len);
#define GwLiteScore_alloc_len Wise2_GwLiteScore_alloc_len


/* Function:  hard_link_GwLiteScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GwLiteScore *]
 *
 * Return [UNKN ]  Undocumented return value [GwLiteScore *]
 *
 */
GwLiteScore * Wise2_hard_link_GwLiteScore(GwLiteScore * obj);
#define hard_link_GwLiteScore Wise2_hard_link_GwLiteScore


/* Function:  GwLiteScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GwLiteScore *]
 *
 */
GwLiteScore * Wise2_GwLiteScore_alloc(void);
#define GwLiteScore_alloc Wise2_GwLiteScore_alloc


/* Function:  free_GwLiteScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GwLiteScore *]
 *
 * Return [UNKN ]  Undocumented return value [GwLiteScore *]
 *
 */
GwLiteScore * Wise2_free_GwLiteScore(GwLiteScore * obj);
#define free_GwLiteScore Wise2_free_GwLiteScore


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
GwLiteSegmentScore * Wise2_GwLiteSegmentScore_from_GwLiteSegment(GwLiteSegment *prev,GwLiteSegment * seg);
#define GwLiteSegmentScore_from_GwLiteSegment Wise2_GwLiteSegmentScore_from_GwLiteSegment


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_GwLite(GwLiteSegment ** list,int i,int j) ;
#define swap_GwLite Wise2_swap_GwLite
void Wise2_qsort_GwLite(GwLiteSegment ** list,int left,int right,int (*comp)(GwLiteSegment * ,GwLiteSegment * ));
#define qsort_GwLite Wise2_qsort_GwLite
void Wise2_sort_GwLite(GwLite * obj,int (*comp)(GwLiteSegment *, GwLiteSegment *));
#define sort_GwLite Wise2_sort_GwLite
boolean Wise2_expand_GwLite(GwLite * obj,int len);
#define expand_GwLite Wise2_expand_GwLite
void Wise2_swap_GwLiteScore(GwLiteSegmentScore ** list,int i,int j) ;
#define swap_GwLiteScore Wise2_swap_GwLiteScore
void Wise2_qsort_GwLiteScore(GwLiteSegmentScore ** list,int left,int right,int (*comp)(GwLiteSegmentScore * ,GwLiteSegmentScore * ));
#define qsort_GwLiteScore Wise2_qsort_GwLiteScore
void Wise2_sort_GwLiteScore(GwLiteScore * obj,int (*comp)(GwLiteSegmentScore *, GwLiteSegmentScore *));
#define sort_GwLiteScore Wise2_sort_GwLiteScore
boolean Wise2_expand_GwLiteScore(GwLiteScore * obj,int len);
#define expand_GwLiteScore Wise2_expand_GwLiteScore

#ifdef _cplusplus
}
#endif

#endif
