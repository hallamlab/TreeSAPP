#ifndef DYNAMITEdnaprofileHEADERFILE
#define DYNAMITEdnaprofileHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"
#include "transfactor.h"

typedef enum DnaProfileTransition {
  DnaProfile_M2M = 0,
  DnaProfile_M2I,
  DnaProfile_I2I,
  DnaProfile_I2M,
  DnaProfile_M2D,
  DnaProfile_D2D,
  DnaProfile_D2M,
  DnaProfile_TRANS_LENGTH
} DnaProfileTransition;

#define DnaProfileScoreLISTLENGTH 128
#define DnaProfileLISTLENGTH 128

struct Wise2_DnaProfileCol {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Probability emit[5];     
    Probability trans[DnaProfile_TRANS_LENGTH];  
    int seqalign_col;    
    } ;  
/* DnaProfileCol defined */ 
#ifndef DYNAMITE_DEFINED_DnaProfileCol
typedef struct Wise2_DnaProfileCol Wise2_DnaProfileCol;
#define DnaProfileCol Wise2_DnaProfileCol
#define DYNAMITE_DEFINED_DnaProfileCol
#endif


struct Wise2_DnaProfileColScore {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Score emit[5];   
    Score trans[DnaProfile_TRANS_LENGTH];    
    } ;  
/* DnaProfileColScore defined */ 
#ifndef DYNAMITE_DEFINED_DnaProfileColScore
typedef struct Wise2_DnaProfileColScore Wise2_DnaProfileColScore;
#define DnaProfileColScore Wise2_DnaProfileColScore
#define DYNAMITE_DEFINED_DnaProfileColScore
#endif


struct Wise2_DnaProfileScore {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    DnaProfileColScore ** col;   
    int len;/* len for above col  */ 
    int maxlen; /* maxlen for above col */ 
    } ;  
/* DnaProfileScore defined */ 
#ifndef DYNAMITE_DEFINED_DnaProfileScore
typedef struct Wise2_DnaProfileScore Wise2_DnaProfileScore;
#define DnaProfileScore Wise2_DnaProfileScore
#define DYNAMITE_DEFINED_DnaProfileScore
#endif


struct Wise2_DnaProfile {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    DnaProfileCol ** col;    
    int len;/* len for above col  */ 
    int maxlen; /* maxlen for above col */ 
    SeqAlign * sa;   
    boolean folded_random;   
    } ;  
/* DnaProfile defined */ 
#ifndef DYNAMITE_DEFINED_DnaProfile
typedef struct Wise2_DnaProfile Wise2_DnaProfile;
#define DnaProfile Wise2_DnaProfile
#define DYNAMITE_DEFINED_DnaProfile
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_DnaProfileCol(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DnaProfileCol *]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileCol *]
 *
 */
DnaProfileCol * Wise2_hard_link_DnaProfileCol(DnaProfileCol * obj);
#define hard_link_DnaProfileCol Wise2_hard_link_DnaProfileCol


/* Function:  DnaProfileCol_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileCol *]
 *
 */
DnaProfileCol * Wise2_DnaProfileCol_alloc(void);
#define DnaProfileCol_alloc Wise2_DnaProfileCol_alloc


/* Function:  free_DnaProfileCol(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DnaProfileCol *]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileCol *]
 *
 */
DnaProfileCol * Wise2_free_DnaProfileCol(DnaProfileCol * obj);
#define free_DnaProfileCol Wise2_free_DnaProfileCol


/* Function:  hard_link_DnaProfileColScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DnaProfileColScore *]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileColScore *]
 *
 */
DnaProfileColScore * Wise2_hard_link_DnaProfileColScore(DnaProfileColScore * obj);
#define hard_link_DnaProfileColScore Wise2_hard_link_DnaProfileColScore


/* Function:  DnaProfileColScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileColScore *]
 *
 */
DnaProfileColScore * Wise2_DnaProfileColScore_alloc(void);
#define DnaProfileColScore_alloc Wise2_DnaProfileColScore_alloc


/* Function:  free_DnaProfileColScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DnaProfileColScore *]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileColScore *]
 *
 */
DnaProfileColScore * Wise2_free_DnaProfileColScore(DnaProfileColScore * obj);
#define free_DnaProfileColScore Wise2_free_DnaProfileColScore


/* Function:  add_DnaProfileScore(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [DnaProfileScore *]
 * Arg:        add [OWNER] Object to add to the list [DnaProfileColScore *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_DnaProfileScore(DnaProfileScore * obj,DnaProfileColScore * add);
#define add_DnaProfileScore Wise2_add_DnaProfileScore


/* Function:  flush_DnaProfileScore(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [DnaProfileScore *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_DnaProfileScore(DnaProfileScore * obj);
#define flush_DnaProfileScore Wise2_flush_DnaProfileScore


/* Function:  DnaProfileScore_alloc_std(void)
 *
 * Descrip:    Equivalent to DnaProfileScore_alloc_len(DnaProfileScoreLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileScore *]
 *
 */
DnaProfileScore * Wise2_DnaProfileScore_alloc_std(void);
#define DnaProfileScore_alloc_std Wise2_DnaProfileScore_alloc_std


/* Function:  DnaProfileScore_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileScore *]
 *
 */
DnaProfileScore * Wise2_DnaProfileScore_alloc_len(int len);
#define DnaProfileScore_alloc_len Wise2_DnaProfileScore_alloc_len


/* Function:  hard_link_DnaProfileScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DnaProfileScore *]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileScore *]
 *
 */
DnaProfileScore * Wise2_hard_link_DnaProfileScore(DnaProfileScore * obj);
#define hard_link_DnaProfileScore Wise2_hard_link_DnaProfileScore


/* Function:  DnaProfileScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileScore *]
 *
 */
DnaProfileScore * Wise2_DnaProfileScore_alloc(void);
#define DnaProfileScore_alloc Wise2_DnaProfileScore_alloc


/* Function:  free_DnaProfileScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DnaProfileScore *]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileScore *]
 *
 */
DnaProfileScore * Wise2_free_DnaProfileScore(DnaProfileScore * obj);
#define free_DnaProfileScore Wise2_free_DnaProfileScore


/* Function:  add_DnaProfile(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [DnaProfile *]
 * Arg:        add [OWNER] Object to add to the list [DnaProfileCol *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_DnaProfile(DnaProfile * obj,DnaProfileCol * add);
#define add_DnaProfile Wise2_add_DnaProfile


/* Function:  flush_DnaProfile(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [DnaProfile *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_DnaProfile(DnaProfile * obj);
#define flush_DnaProfile Wise2_flush_DnaProfile


/* Function:  DnaProfile_alloc_std(void)
 *
 * Descrip:    Equivalent to DnaProfile_alloc_len(DnaProfileLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaProfile *]
 *
 */
DnaProfile * Wise2_DnaProfile_alloc_std(void);
#define DnaProfile_alloc_std Wise2_DnaProfile_alloc_std


/* Function:  DnaProfile_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfile *]
 *
 */
DnaProfile * Wise2_DnaProfile_alloc_len(int len);
#define DnaProfile_alloc_len Wise2_DnaProfile_alloc_len


/* Function:  hard_link_DnaProfile(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DnaProfile *]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfile *]
 *
 */
DnaProfile * Wise2_hard_link_DnaProfile(DnaProfile * obj);
#define hard_link_DnaProfile Wise2_hard_link_DnaProfile


/* Function:  DnaProfile_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaProfile *]
 *
 */
DnaProfile * Wise2_DnaProfile_alloc(void);
#define DnaProfile_alloc Wise2_DnaProfile_alloc


/* Function:  free_DnaProfile(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DnaProfile *]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfile *]
 *
 */
DnaProfile * Wise2_free_DnaProfile(DnaProfile * obj);
#define free_DnaProfile Wise2_free_DnaProfile


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
TransFactor * Wise2_new_TransFactor_from_DnaProfile(DnaProfile * dp,char * name);
#define new_TransFactor_from_DnaProfile Wise2_new_TransFactor_from_DnaProfile
double Wise2_bits_entropy_DnaProfile(DnaProfile * dnap,RandomModelDNA * background);
#define bits_entropy_DnaProfile Wise2_bits_entropy_DnaProfile
SeqAlign * Wise2_merged_SeqAlign(DnaProfile * one,DnaProfile * two,AlnBlock * alb);
#define merged_SeqAlign Wise2_merged_SeqAlign
DnaProfileScore * Wise2_DnaProfileScore_from_DnaProfile(DnaProfile * dp);
#define DnaProfileScore_from_DnaProfile Wise2_DnaProfileScore_from_DnaProfile
void Wise2_fold_RandomModel_DnaProfile(DnaProfile * dp,RandomModelDNA * rm);
#define fold_RandomModel_DnaProfile Wise2_fold_RandomModel_DnaProfile
DnaProfile * Wise2_naive_DnaProfile_from_Sequence(Sequence * seq,double id,double m2i,double m2d,double i2i,double d2d);
#define naive_DnaProfile_from_Sequence Wise2_naive_DnaProfile_from_Sequence
DnaProfile * Wise2_naive_DnaProfile_from_SeqAlign(SeqAlign * sa,double simple_pseudo,double m2i,double m2d,double i2i,double d2d);
#define naive_DnaProfile_from_SeqAlign Wise2_naive_DnaProfile_from_SeqAlign
void Wise2_show_DnaProfile(DnaProfile * dnap,RandomModelDNA * rm,FILE * ofp);
#define show_DnaProfile Wise2_show_DnaProfile


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_DnaProfileScore(DnaProfileColScore ** list,int i,int j) ;
#define swap_DnaProfileScore Wise2_swap_DnaProfileScore
void Wise2_qsort_DnaProfileScore(DnaProfileColScore ** list,int left,int right,int (*comp)(DnaProfileColScore * ,DnaProfileColScore * ));
#define qsort_DnaProfileScore Wise2_qsort_DnaProfileScore
void Wise2_sort_DnaProfileScore(DnaProfileScore * obj,int (*comp)(DnaProfileColScore *, DnaProfileColScore *));
#define sort_DnaProfileScore Wise2_sort_DnaProfileScore
boolean Wise2_expand_DnaProfileScore(DnaProfileScore * obj,int len);
#define expand_DnaProfileScore Wise2_expand_DnaProfileScore
void Wise2_swap_DnaProfile(DnaProfileCol ** list,int i,int j) ;
#define swap_DnaProfile Wise2_swap_DnaProfile
void Wise2_qsort_DnaProfile(DnaProfileCol ** list,int left,int right,int (*comp)(DnaProfileCol * ,DnaProfileCol * ));
#define qsort_DnaProfile Wise2_qsort_DnaProfile
void Wise2_sort_DnaProfile(DnaProfile * obj,int (*comp)(DnaProfileCol *, DnaProfileCol *));
#define sort_DnaProfile Wise2_sort_DnaProfile
boolean Wise2_expand_DnaProfile(DnaProfile * obj,int len);
#define expand_DnaProfile Wise2_expand_DnaProfile

#ifdef _cplusplus
}
#endif

#endif
