#ifndef DYNAMITEhsphandlerHEADERFILE
#define DYNAMITEhsphandlerHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif

#include "sequence.h"
#include "compmat.h"
#include "glib.h"
#include "hsp.h"

#define UNFEASIBLY_LARGE_SCORE 1000000



struct Wise2_TopScoreManager {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int * score;     
    int length;  
    int current_pos;     
    int worst_score;     
    int worst_position;  
    } ;  
/* TopScoreManager defined */ 
#ifndef DYNAMITE_DEFINED_TopScoreManager
typedef struct Wise2_TopScoreManager Wise2_TopScoreManager;
#define TopScoreManager Wise2_TopScoreManager
#define DYNAMITE_DEFINED_TopScoreManager
#endif


struct Wise2_QuerySeqHSP {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Sequence * query;    
    int ** score;    
    int * self_score;    
    int len;     
    } ;  
/* QuerySeqHSP defined */ 
#ifndef DYNAMITE_DEFINED_QuerySeqHSP
typedef struct Wise2_QuerySeqHSP Wise2_QuerySeqHSP;
#define QuerySeqHSP Wise2_QuerySeqHSP
#define DYNAMITE_DEFINED_QuerySeqHSP
#endif


struct Wise2_HSPmanager {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    GHashTable * target_hash;    
    Sequence * query;    
    CompMat * mat;   
    int drop_off;    
    int min_score;   
    TopScoreManager * tsm;   
    QuerySeqHSP * qs;    
    HSPCache * cache;    
    int hsp_count;   
    } ;  
/* HSPmanager defined */ 
#ifndef DYNAMITE_DEFINED_HSPmanager
typedef struct Wise2_HSPmanager Wise2_HSPmanager;
#define HSPmanager Wise2_HSPmanager
#define DYNAMITE_DEFINED_HSPmanager
#endif


struct Wise2_LineariseHSPPara {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int max_size;    
    int min_score;   
    int width;   
    int tail;    
    int verbosity;   
    } ;  
/* LineariseHSPPara defined */ 
#ifndef DYNAMITE_DEFINED_LineariseHSPPara
typedef struct Wise2_LineariseHSPPara Wise2_LineariseHSPPara;
#define LineariseHSPPara Wise2_LineariseHSPPara
#define DYNAMITE_DEFINED_LineariseHSPPara
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  new_TopScoreManager(length)
 *
 * Descrip:    Makes a new topscore manager of a specific length
 *
 *
 * Arg:        length [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [TopScoreManager *]
 *
 */
TopScoreManager * Wise2_new_TopScoreManager(int length);
#define new_TopScoreManager Wise2_new_TopScoreManager


/* Function:  add_score_TopScoreManager(tsm,score)
 *
 * Descrip:    Adds a top score
 *
 *
 * Arg:          tsm [UNKN ] Undocumented argument [TopScoreManager *]
 * Arg:        score [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_score_TopScoreManager(TopScoreManager * tsm,int score);
#define add_score_TopScoreManager Wise2_add_score_TopScoreManager


/* Function:  truncated_simple_LinearHSPmanager(lm,para)
 *
 * Descrip:    A simpler truncation method, using diagonals
 *
 *
 * Arg:          lm [UNKN ] Undocumented argument [LinearHSPmanager *]
 * Arg:        para [UNKN ] Undocumented argument [LineariseHSPPara *]
 *
 * Return [UNKN ]  Undocumented return value [LinearHSPmanager *]
 *
 */
LinearHSPmanager * Wise2_truncated_simple_LinearHSPmanager(LinearHSPmanager * lm,LineariseHSPPara * para);
#define truncated_simple_LinearHSPmanager Wise2_truncated_simple_LinearHSPmanager


/* Function:  truncated_LinearHSPmanager(lm,max_size,min_score,width,tail)
 *
 * Descrip:    Makes a truncated linear set 
 *
 *
 * Arg:               lm [UNKN ] Undocumented argument [LinearHSPmanager *]
 * Arg:         max_size [UNKN ] Undocumented argument [int]
 * Arg:        min_score [UNKN ] Undocumented argument [int]
 * Arg:            width [UNKN ] Undocumented argument [int]
 * Arg:             tail [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [LinearHSPmanager *]
 *
 */
LinearHSPmanager * Wise2_truncated_LinearHSPmanager(LinearHSPmanager * lm,int max_size,int min_score,int width,int tail);
#define truncated_LinearHSPmanager Wise2_truncated_LinearHSPmanager


/* Function:  new_consistent_HSPset(set,min_score,width,tail)
 *
 * Descrip:    Makes an HSP set via heuristics
 *             to deal with low complexity regions
 *
 *
 * Arg:              set [UNKN ] Undocumented argument [HSPset *]
 * Arg:        min_score [UNKN ] Undocumented argument [int]
 * Arg:            width [UNKN ] Undocumented argument [int]
 * Arg:             tail [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [HSPset *]
 *
 */
HSPset * Wise2_new_consistent_HSPset(HSPset * set,int min_score,int width,int tail);
#define new_consistent_HSPset Wise2_new_consistent_HSPset


/* Function:  new_LinearHSPmanager_simple_heuristic(hspm,para)
 *
 * Descrip:    New, simpler LinearHSPmanager with diagonal heuristics
 *
 *
 * Arg:        hspm [UNKN ] Undocumented argument [HSPmanager *]
 * Arg:        para [UNKN ] Undocumented argument [LineariseHSPPara *]
 *
 * Return [UNKN ]  Undocumented return value [LinearHSPmanager *]
 *
 */
LinearHSPmanager * Wise2_new_LinearHSPmanager_simple_heuristic(HSPmanager * hspm,LineariseHSPPara * para);
#define new_LinearHSPmanager_simple_heuristic Wise2_new_LinearHSPmanager_simple_heuristic


/* Function:  new_LinearHSPmanager_heuristic_max(hspm,max_size)
 *
 * Descrip:    Builds a new LinearHSPmanager from hash based with heuristics 
 *
 *
 * Arg:            hspm [UNKN ] Undocumented argument [HSPmanager *]
 * Arg:        max_size [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [LinearHSPmanager *]
 *
 */
LinearHSPmanager * Wise2_new_LinearHSPmanager_heuristic_max(HSPmanager * hspm,int max_size);
#define new_LinearHSPmanager_heuristic_max Wise2_new_LinearHSPmanager_heuristic_max


/* Function:  new_LinearHSPmanager_truncate_on_score(hspm)
 *
 * Descrip:    Builds a LinearHSPmanager from a hash based HSP manager, using worst score truncation
 *
 *
 * Arg:        hspm [UNKN ] Undocumented argument [HSPmanager *]
 *
 * Return [UNKN ]  Undocumented return value [LinearHSPmanager *]
 *
 */
LinearHSPmanager * Wise2_new_LinearHSPmanager_truncate_on_score(HSPmanager * hspm);
#define new_LinearHSPmanager_truncate_on_score Wise2_new_LinearHSPmanager_truncate_on_score


/* Function:  new_LinearHSPmanager_flat(hspm)
 *
 * Descrip:    Builds a LinearHSPmanager from a hash based HSP manager
 *
 *
 * Arg:        hspm [UNKN ] Undocumented argument [HSPmanager *]
 *
 * Return [UNKN ]  Undocumented return value [LinearHSPmanager *]
 *
 */
LinearHSPmanager * Wise2_new_LinearHSPmanager_flat(HSPmanager * hspm);
#define new_LinearHSPmanager_flat Wise2_new_LinearHSPmanager_flat


/* Function:  linearise_HSPset_truncate_on_score(key,value,user_data)
 *
 * Descrip:    internal function for remapping HSPs
 *
 *
 * Arg:              key [UNKN ] Undocumented argument [gpointer]
 * Arg:            value [UNKN ] Undocumented argument [gpointer]
 * Arg:        user_data [UNKN ] Undocumented argument [gpointer]
 *
 */
void Wise2_linearise_HSPset_truncate_on_score(gpointer key,gpointer value,gpointer user_data);
#define linearise_HSPset_truncate_on_score Wise2_linearise_HSPset_truncate_on_score


/* Function:  linearise_HSPset_flat(key,value,user_data)
 *
 * Descrip:    internal function for remapping HSPs with score cutoff
 *
 *
 * Arg:              key [UNKN ] Undocumented argument [gpointer]
 * Arg:            value [UNKN ] Undocumented argument [gpointer]
 * Arg:        user_data [UNKN ] Undocumented argument [gpointer]
 *
 */
void Wise2_linearise_HSPset_flat(gpointer key,gpointer value,gpointer user_data);
#define linearise_HSPset_flat Wise2_linearise_HSPset_flat


/* Function:  linearise_HSPset_consistent(key,value,user_data)
 *
 * Descrip:    internal function for remapping HSPs
 *
 *
 * Arg:              key [UNKN ] Undocumented argument [gpointer]
 * Arg:            value [UNKN ] Undocumented argument [gpointer]
 * Arg:        user_data [UNKN ] Undocumented argument [gpointer]
 *
 */
void Wise2_linearise_HSPset_consistent(gpointer key,gpointer value,gpointer user_data);
#define linearise_HSPset_consistent Wise2_linearise_HSPset_consistent


/* Function:  new_HSPmanager(query,mat,score_drop_off)
 *
 * Descrip:    Builds a new HSPmanager for a target system
 *
 *
 * Arg:                 query [UNKN ] Undocumented argument [Sequence *]
 * Arg:                   mat [UNKN ] Undocumented argument [CompMat *]
 * Arg:        score_drop_off [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [HSPmanager *]
 *
 */
HSPmanager * Wise2_new_HSPmanager(Sequence * query,CompMat * mat,int score_drop_off);
#define new_HSPmanager Wise2_new_HSPmanager


/* Function:  add_pair_HSPmanager(hspm,target,query_pos,target_pos)
 *
 * Descrip:    adds a new target pair, irregardless of score
 *
 *
 * Arg:              hspm [UNKN ] Undocumented argument [HSPmanager *]
 * Arg:            target [UNKN ] Undocumented argument [Sequence *]
 * Arg:         query_pos [UNKN ] Undocumented argument [int]
 * Arg:        target_pos [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_add_pair_HSPmanager(HSPmanager * hspm,Sequence * target,int query_pos,int target_pos);
#define add_pair_HSPmanager Wise2_add_pair_HSPmanager


/* Function:  add_pair_HSPmanager_score(hspm,target,query_pos,target_pos,min_score)
 *
 * Descrip:    Adds a new target pair to this HSPmanager for indexing, with a min score
 *
 *
 * Arg:              hspm [UNKN ] Undocumented argument [HSPmanager *]
 * Arg:            target [UNKN ] Undocumented argument [Sequence *]
 * Arg:         query_pos [UNKN ] Undocumented argument [int]
 * Arg:        target_pos [UNKN ] Undocumented argument [int]
 * Arg:         min_score [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_add_pair_HSPmanager_score(HSPmanager * hspm,Sequence * target,int query_pos,int target_pos,int min_score);
#define add_pair_HSPmanager_score Wise2_add_pair_HSPmanager_score


/* Function:  add_new_HSP_HSPmanager(hspm,target,query_start,target_start,length,score)
 *
 * Descrip:    Adds a new HSP when all info is known
 *
 *
 * Arg:                hspm [UNKN ] Undocumented argument [HSPmanager *]
 * Arg:              target [UNKN ] Undocumented argument [Sequence *]
 * Arg:         query_start [UNKN ] Undocumented argument [int]
 * Arg:        target_start [UNKN ] Undocumented argument [int]
 * Arg:              length [UNKN ] Undocumented argument [int]
 * Arg:               score [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_new_HSP_HSPmanager(HSPmanager * hspm,Sequence * target,int query_start,int target_start,int length,int score);
#define add_new_HSP_HSPmanager Wise2_add_new_HSP_HSPmanager


/* Function:  free_ghash_HSPsets(key,value,user_data)
 *
 * Descrip:    Frees the HSPsets
 *
 *
 * Arg:              key [UNKN ] Undocumented argument [gpointer]
 * Arg:            value [UNKN ] Undocumented argument [gpointer]
 * Arg:        user_data [UNKN ] Undocumented argument [gpointer]
 *
 */
void Wise2_free_ghash_HSPsets(gpointer key,gpointer value,gpointer user_data);
#define free_ghash_HSPsets Wise2_free_ghash_HSPsets


/* Function:  free_HSPmanager(h)
 *
 * Descrip:    Frees the HSPmanager
 *
 *
 * Arg:        h [UNKN ] Undocumented argument [HSPmanager *]
 *
 * Return [UNKN ]  Undocumented return value [HSPmanager *]
 *
 */
HSPmanager * Wise2_free_HSPmanager(HSPmanager * h);
#define free_HSPmanager Wise2_free_HSPmanager


/* Function:  new_HSP_QuerySeqHSP(cache,query,target,query_pos,target_pos,mat,score_drop_off,min_score)
 *
 * Descrip:    builds a new HSP for these sequences
 *
 *
 * Arg:                 cache [UNKN ] Undocumented argument [HSPCache *]
 * Arg:                 query [UNKN ] Undocumented argument [QuerySeqHSP *]
 * Arg:                target [UNKN ] Undocumented argument [Sequence *]
 * Arg:             query_pos [UNKN ] Undocumented argument [int]
 * Arg:            target_pos [UNKN ] Undocumented argument [int]
 * Arg:                   mat [UNKN ] Undocumented argument [CompMat *]
 * Arg:        score_drop_off [UNKN ] Undocumented argument [int]
 * Arg:             min_score [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [HSP *]
 *
 */
HSP * Wise2_new_HSP_QuerySeqHSP(HSPCache * cache,QuerySeqHSP * query,Sequence * target,int query_pos,int target_pos,CompMat * mat,int score_drop_off,int min_score);
#define new_HSP_QuerySeqHSP Wise2_new_HSP_QuerySeqHSP


/* Function:  new_QuerySeqHSP(seq,mat)
 *
 * Descrip:    Builds a new QuerySeqHSP from Sequence and Matrix
 *
 *
 * Arg:        seq [UNKN ] Undocumented argument [Sequence *]
 * Arg:        mat [UNKN ] Undocumented argument [CompMat *]
 *
 * Return [UNKN ]  Undocumented return value [QuerySeqHSP *]
 *
 */
QuerySeqHSP * Wise2_new_QuerySeqHSP(Sequence * seq,CompMat * mat);
#define new_QuerySeqHSP Wise2_new_QuerySeqHSP


/* Function:  hard_link_TopScoreManager(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [TopScoreManager *]
 *
 * Return [UNKN ]  Undocumented return value [TopScoreManager *]
 *
 */
TopScoreManager * Wise2_hard_link_TopScoreManager(TopScoreManager * obj);
#define hard_link_TopScoreManager Wise2_hard_link_TopScoreManager


/* Function:  TopScoreManager_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TopScoreManager *]
 *
 */
TopScoreManager * Wise2_TopScoreManager_alloc(void);
#define TopScoreManager_alloc Wise2_TopScoreManager_alloc


/* Function:  free_TopScoreManager(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [TopScoreManager *]
 *
 * Return [UNKN ]  Undocumented return value [TopScoreManager *]
 *
 */
TopScoreManager * Wise2_free_TopScoreManager(TopScoreManager * obj);
#define free_TopScoreManager Wise2_free_TopScoreManager


/* Function:  hard_link_QuerySeqHSP(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [QuerySeqHSP *]
 *
 * Return [UNKN ]  Undocumented return value [QuerySeqHSP *]
 *
 */
QuerySeqHSP * Wise2_hard_link_QuerySeqHSP(QuerySeqHSP * obj);
#define hard_link_QuerySeqHSP Wise2_hard_link_QuerySeqHSP


/* Function:  QuerySeqHSP_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [QuerySeqHSP *]
 *
 */
QuerySeqHSP * Wise2_QuerySeqHSP_alloc(void);
#define QuerySeqHSP_alloc Wise2_QuerySeqHSP_alloc


/* Function:  free_QuerySeqHSP(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [QuerySeqHSP *]
 *
 * Return [UNKN ]  Undocumented return value [QuerySeqHSP *]
 *
 */
QuerySeqHSP * Wise2_free_QuerySeqHSP(QuerySeqHSP * obj);
#define free_QuerySeqHSP Wise2_free_QuerySeqHSP


/* Function:  hard_link_HSPmanager(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [HSPmanager *]
 *
 * Return [UNKN ]  Undocumented return value [HSPmanager *]
 *
 */
HSPmanager * Wise2_hard_link_HSPmanager(HSPmanager * obj);
#define hard_link_HSPmanager Wise2_hard_link_HSPmanager


/* Function:  HSPmanager_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [HSPmanager *]
 *
 */
HSPmanager * Wise2_HSPmanager_alloc(void);
#define HSPmanager_alloc Wise2_HSPmanager_alloc


/* Function:  hard_link_LineariseHSPPara(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [LineariseHSPPara *]
 *
 * Return [UNKN ]  Undocumented return value [LineariseHSPPara *]
 *
 */
LineariseHSPPara * Wise2_hard_link_LineariseHSPPara(LineariseHSPPara * obj);
#define hard_link_LineariseHSPPara Wise2_hard_link_LineariseHSPPara


/* Function:  LineariseHSPPara_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LineariseHSPPara *]
 *
 */
LineariseHSPPara * Wise2_LineariseHSPPara_alloc(void);
#define LineariseHSPPara_alloc Wise2_LineariseHSPPara_alloc


/* Function:  free_LineariseHSPPara(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [LineariseHSPPara *]
 *
 * Return [UNKN ]  Undocumented return value [LineariseHSPPara *]
 *
 */
LineariseHSPPara * Wise2_free_LineariseHSPPara(LineariseHSPPara * obj);
#define free_LineariseHSPPara Wise2_free_LineariseHSPPara


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/

#ifdef _cplusplus
}
#endif

#endif
