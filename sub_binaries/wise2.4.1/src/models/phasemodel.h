#ifndef DYNAMITEphasemodelHEADERFILE
#define DYNAMITEphasemodelHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "genewisemodel.h"


#define GenePhaseModelLISTLENGTH 128
#define GenePhaseScoreLISTLENGTH 128

#define ProteinIntronListLISTLENGTH 30

struct Wise2_ProteinIntron {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int position;    
    int phase;   
    } ;  
/* ProteinIntron defined */ 
#ifndef DYNAMITE_DEFINED_ProteinIntron
typedef struct Wise2_ProteinIntron Wise2_ProteinIntron;
#define ProteinIntron Wise2_ProteinIntron
#define DYNAMITE_DEFINED_ProteinIntron
#endif


struct Wise2_ProteinIntronList {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    ProteinIntron ** intron;     
    int len;/* len for above intron  */ 
    int maxlen; /* maxlen for above intron */ 
    } ;  
/* ProteinIntronList defined */ 
#ifndef DYNAMITE_DEFINED_ProteinIntronList
typedef struct Wise2_ProteinIntronList Wise2_ProteinIntronList;
#define ProteinIntronList Wise2_ProteinIntronList
#define DYNAMITE_DEFINED_ProteinIntronList
#endif


struct Wise2_PhasedProtein {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Sequence * protein;  
    ProteinIntronList * list;    
    } ;  
/* PhasedProtein defined */ 
#ifndef DYNAMITE_DEFINED_PhasedProtein
typedef struct Wise2_PhasedProtein Wise2_PhasedProtein;
#define PhasedProtein Wise2_PhasedProtein
#define DYNAMITE_DEFINED_PhasedProtein
#endif


struct Wise2_PhasedProteinPara {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Probability marked_intron;   
    Probability unmarked_intron;     
    int gap;     
    int ext;     
    boolean use_phase;   
    char * intron_file;  
    } ;  
/* PhasedProteinPara defined */ 
#ifndef DYNAMITE_DEFINED_PhasedProteinPara
typedef struct Wise2_PhasedProteinPara Wise2_PhasedProteinPara;
#define PhasedProteinPara Wise2_PhasedProteinPara
#define DYNAMITE_DEFINED_PhasedProteinPara
#endif


struct Wise2_PhasedHMM {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    ThreeStateModel * tsm;   
    ProteinIntronList * list;    
    } ;  
/* PhasedHMM defined */ 
#ifndef DYNAMITE_DEFINED_PhasedHMM
typedef struct Wise2_PhasedHMM Wise2_PhasedHMM;
#define PhasedHMM Wise2_PhasedHMM
#define DYNAMITE_DEFINED_PhasedHMM
#endif


struct Wise2_GenePhaseSeg {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Probability intron_0;    
    Probability intron_1;    
    Probability intron_2;    
    Probability insert_intron;   
    } ;  
/* GenePhaseSeg defined */ 
#ifndef DYNAMITE_DEFINED_GenePhaseSeg
typedef struct Wise2_GenePhaseSeg Wise2_GenePhaseSeg;
#define GenePhaseSeg Wise2_GenePhaseSeg
#define DYNAMITE_DEFINED_GenePhaseSeg
#endif


struct Wise2_GenePhaseModel {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    GenePhaseSeg ** phase;   
    int len;/* len for above phase  */ 
    int maxlen; /* maxlen for above phase */ 
    GeneWise * gw;   
    } ;  
/* GenePhaseModel defined */ 
#ifndef DYNAMITE_DEFINED_GenePhaseModel
typedef struct Wise2_GenePhaseModel Wise2_GenePhaseModel;
#define GenePhaseModel Wise2_GenePhaseModel
#define DYNAMITE_DEFINED_GenePhaseModel
#endif


struct Wise2_GenePhaseSegScore {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Score intron_0;  
    Score intron_1;  
    Score intron_2;  
    Score insert_intron;     
    } ;  
/* GenePhaseSegScore defined */ 
#ifndef DYNAMITE_DEFINED_GenePhaseSegScore
typedef struct Wise2_GenePhaseSegScore Wise2_GenePhaseSegScore;
#define GenePhaseSegScore Wise2_GenePhaseSegScore
#define DYNAMITE_DEFINED_GenePhaseSegScore
#endif


struct Wise2_GenePhaseScore {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    GenePhaseSegScore ** phase;  
    int len;/* len for above phase  */ 
    int maxlen; /* maxlen for above phase */ 
    GeneWiseScore * gws;     
    } ;  
/* GenePhaseScore defined */ 
#ifndef DYNAMITE_DEFINED_GenePhaseScore
typedef struct Wise2_GenePhaseScore Wise2_GenePhaseScore;
#define GenePhaseScore Wise2_GenePhaseScore
#define DYNAMITE_DEFINED_GenePhaseScore
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_ProteinIntron(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ProteinIntron *]
 *
 * Return [UNKN ]  Undocumented return value [ProteinIntron *]
 *
 */
ProteinIntron * Wise2_hard_link_ProteinIntron(ProteinIntron * obj);
#define hard_link_ProteinIntron Wise2_hard_link_ProteinIntron


/* Function:  ProteinIntron_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ProteinIntron *]
 *
 */
ProteinIntron * Wise2_ProteinIntron_alloc(void);
#define ProteinIntron_alloc Wise2_ProteinIntron_alloc


/* Function:  free_ProteinIntron(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ProteinIntron *]
 *
 * Return [UNKN ]  Undocumented return value [ProteinIntron *]
 *
 */
ProteinIntron * Wise2_free_ProteinIntron(ProteinIntron * obj);
#define free_ProteinIntron Wise2_free_ProteinIntron


/* Function:  add_ProteinIntronList(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [ProteinIntronList *]
 * Arg:        add [OWNER] Object to add to the list [ProteinIntron *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_ProteinIntronList(ProteinIntronList * obj,ProteinIntron * add);
#define add_ProteinIntronList Wise2_add_ProteinIntronList


/* Function:  flush_ProteinIntronList(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [ProteinIntronList *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_ProteinIntronList(ProteinIntronList * obj);
#define flush_ProteinIntronList Wise2_flush_ProteinIntronList


/* Function:  ProteinIntronList_alloc_std(void)
 *
 * Descrip:    Equivalent to ProteinIntronList_alloc_len(ProteinIntronListLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ProteinIntronList *]
 *
 */
ProteinIntronList * Wise2_ProteinIntronList_alloc_std(void);
#define ProteinIntronList_alloc_std Wise2_ProteinIntronList_alloc_std


/* Function:  ProteinIntronList_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [ProteinIntronList *]
 *
 */
ProteinIntronList * Wise2_ProteinIntronList_alloc_len(int len);
#define ProteinIntronList_alloc_len Wise2_ProteinIntronList_alloc_len


/* Function:  hard_link_ProteinIntronList(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ProteinIntronList *]
 *
 * Return [UNKN ]  Undocumented return value [ProteinIntronList *]
 *
 */
ProteinIntronList * Wise2_hard_link_ProteinIntronList(ProteinIntronList * obj);
#define hard_link_ProteinIntronList Wise2_hard_link_ProteinIntronList


/* Function:  ProteinIntronList_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ProteinIntronList *]
 *
 */
ProteinIntronList * Wise2_ProteinIntronList_alloc(void);
#define ProteinIntronList_alloc Wise2_ProteinIntronList_alloc


/* Function:  free_ProteinIntronList(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ProteinIntronList *]
 *
 * Return [UNKN ]  Undocumented return value [ProteinIntronList *]
 *
 */
ProteinIntronList * Wise2_free_ProteinIntronList(ProteinIntronList * obj);
#define free_ProteinIntronList Wise2_free_ProteinIntronList


/* Function:  hard_link_PhasedProtein(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [PhasedProtein *]
 *
 * Return [UNKN ]  Undocumented return value [PhasedProtein *]
 *
 */
PhasedProtein * Wise2_hard_link_PhasedProtein(PhasedProtein * obj);
#define hard_link_PhasedProtein Wise2_hard_link_PhasedProtein


/* Function:  PhasedProtein_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PhasedProtein *]
 *
 */
PhasedProtein * Wise2_PhasedProtein_alloc(void);
#define PhasedProtein_alloc Wise2_PhasedProtein_alloc


/* Function:  free_PhasedProtein(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [PhasedProtein *]
 *
 * Return [UNKN ]  Undocumented return value [PhasedProtein *]
 *
 */
PhasedProtein * Wise2_free_PhasedProtein(PhasedProtein * obj);
#define free_PhasedProtein Wise2_free_PhasedProtein


/* Function:  hard_link_PhasedProteinPara(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [PhasedProteinPara *]
 *
 * Return [UNKN ]  Undocumented return value [PhasedProteinPara *]
 *
 */
PhasedProteinPara * Wise2_hard_link_PhasedProteinPara(PhasedProteinPara * obj);
#define hard_link_PhasedProteinPara Wise2_hard_link_PhasedProteinPara


/* Function:  PhasedProteinPara_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PhasedProteinPara *]
 *
 */
PhasedProteinPara * Wise2_PhasedProteinPara_alloc(void);
#define PhasedProteinPara_alloc Wise2_PhasedProteinPara_alloc


/* Function:  free_PhasedProteinPara(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [PhasedProteinPara *]
 *
 * Return [UNKN ]  Undocumented return value [PhasedProteinPara *]
 *
 */
PhasedProteinPara * Wise2_free_PhasedProteinPara(PhasedProteinPara * obj);
#define free_PhasedProteinPara Wise2_free_PhasedProteinPara


/* Function:  hard_link_PhasedHMM(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [PhasedHMM *]
 *
 * Return [UNKN ]  Undocumented return value [PhasedHMM *]
 *
 */
PhasedHMM * Wise2_hard_link_PhasedHMM(PhasedHMM * obj);
#define hard_link_PhasedHMM Wise2_hard_link_PhasedHMM


/* Function:  PhasedHMM_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PhasedHMM *]
 *
 */
PhasedHMM * Wise2_PhasedHMM_alloc(void);
#define PhasedHMM_alloc Wise2_PhasedHMM_alloc


/* Function:  free_PhasedHMM(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [PhasedHMM *]
 *
 * Return [UNKN ]  Undocumented return value [PhasedHMM *]
 *
 */
PhasedHMM * Wise2_free_PhasedHMM(PhasedHMM * obj);
#define free_PhasedHMM Wise2_free_PhasedHMM


/* Function:  hard_link_GenePhaseSeg(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GenePhaseSeg *]
 *
 * Return [UNKN ]  Undocumented return value [GenePhaseSeg *]
 *
 */
GenePhaseSeg * Wise2_hard_link_GenePhaseSeg(GenePhaseSeg * obj);
#define hard_link_GenePhaseSeg Wise2_hard_link_GenePhaseSeg


/* Function:  GenePhaseSeg_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GenePhaseSeg *]
 *
 */
GenePhaseSeg * Wise2_GenePhaseSeg_alloc(void);
#define GenePhaseSeg_alloc Wise2_GenePhaseSeg_alloc


/* Function:  free_GenePhaseSeg(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GenePhaseSeg *]
 *
 * Return [UNKN ]  Undocumented return value [GenePhaseSeg *]
 *
 */
GenePhaseSeg * Wise2_free_GenePhaseSeg(GenePhaseSeg * obj);
#define free_GenePhaseSeg Wise2_free_GenePhaseSeg


/* Function:  add_GenePhaseModel(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GenePhaseModel *]
 * Arg:        add [OWNER] Object to add to the list [GenePhaseSeg *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_GenePhaseModel(GenePhaseModel * obj,GenePhaseSeg * add);
#define add_GenePhaseModel Wise2_add_GenePhaseModel


/* Function:  flush_GenePhaseModel(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [GenePhaseModel *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_GenePhaseModel(GenePhaseModel * obj);
#define flush_GenePhaseModel Wise2_flush_GenePhaseModel


/* Function:  GenePhaseModel_alloc_std(void)
 *
 * Descrip:    Equivalent to GenePhaseModel_alloc_len(GenePhaseModelLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GenePhaseModel *]
 *
 */
GenePhaseModel * Wise2_GenePhaseModel_alloc_std(void);
#define GenePhaseModel_alloc_std Wise2_GenePhaseModel_alloc_std


/* Function:  GenePhaseModel_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [GenePhaseModel *]
 *
 */
GenePhaseModel * Wise2_GenePhaseModel_alloc_len(int len);
#define GenePhaseModel_alloc_len Wise2_GenePhaseModel_alloc_len


/* Function:  hard_link_GenePhaseModel(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GenePhaseModel *]
 *
 * Return [UNKN ]  Undocumented return value [GenePhaseModel *]
 *
 */
GenePhaseModel * Wise2_hard_link_GenePhaseModel(GenePhaseModel * obj);
#define hard_link_GenePhaseModel Wise2_hard_link_GenePhaseModel


/* Function:  GenePhaseModel_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GenePhaseModel *]
 *
 */
GenePhaseModel * Wise2_GenePhaseModel_alloc(void);
#define GenePhaseModel_alloc Wise2_GenePhaseModel_alloc


/* Function:  free_GenePhaseModel(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GenePhaseModel *]
 *
 * Return [UNKN ]  Undocumented return value [GenePhaseModel *]
 *
 */
GenePhaseModel * Wise2_free_GenePhaseModel(GenePhaseModel * obj);
#define free_GenePhaseModel Wise2_free_GenePhaseModel


/* Function:  hard_link_GenePhaseSegScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GenePhaseSegScore *]
 *
 * Return [UNKN ]  Undocumented return value [GenePhaseSegScore *]
 *
 */
GenePhaseSegScore * Wise2_hard_link_GenePhaseSegScore(GenePhaseSegScore * obj);
#define hard_link_GenePhaseSegScore Wise2_hard_link_GenePhaseSegScore


/* Function:  GenePhaseSegScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GenePhaseSegScore *]
 *
 */
GenePhaseSegScore * Wise2_GenePhaseSegScore_alloc(void);
#define GenePhaseSegScore_alloc Wise2_GenePhaseSegScore_alloc


/* Function:  free_GenePhaseSegScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GenePhaseSegScore *]
 *
 * Return [UNKN ]  Undocumented return value [GenePhaseSegScore *]
 *
 */
GenePhaseSegScore * Wise2_free_GenePhaseSegScore(GenePhaseSegScore * obj);
#define free_GenePhaseSegScore Wise2_free_GenePhaseSegScore


/* Function:  add_GenePhaseScore(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GenePhaseScore *]
 * Arg:        add [OWNER] Object to add to the list [GenePhaseSegScore *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_GenePhaseScore(GenePhaseScore * obj,GenePhaseSegScore * add);
#define add_GenePhaseScore Wise2_add_GenePhaseScore


/* Function:  flush_GenePhaseScore(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [GenePhaseScore *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_GenePhaseScore(GenePhaseScore * obj);
#define flush_GenePhaseScore Wise2_flush_GenePhaseScore


/* Function:  GenePhaseScore_alloc_std(void)
 *
 * Descrip:    Equivalent to GenePhaseScore_alloc_len(GenePhaseScoreLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GenePhaseScore *]
 *
 */
GenePhaseScore * Wise2_GenePhaseScore_alloc_std(void);
#define GenePhaseScore_alloc_std Wise2_GenePhaseScore_alloc_std


/* Function:  GenePhaseScore_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [GenePhaseScore *]
 *
 */
GenePhaseScore * Wise2_GenePhaseScore_alloc_len(int len);
#define GenePhaseScore_alloc_len Wise2_GenePhaseScore_alloc_len


/* Function:  hard_link_GenePhaseScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GenePhaseScore *]
 *
 * Return [UNKN ]  Undocumented return value [GenePhaseScore *]
 *
 */
GenePhaseScore * Wise2_hard_link_GenePhaseScore(GenePhaseScore * obj);
#define hard_link_GenePhaseScore Wise2_hard_link_GenePhaseScore


/* Function:  GenePhaseScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GenePhaseScore *]
 *
 */
GenePhaseScore * Wise2_GenePhaseScore_alloc(void);
#define GenePhaseScore_alloc Wise2_GenePhaseScore_alloc


/* Function:  free_GenePhaseScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GenePhaseScore *]
 *
 * Return [UNKN ]  Undocumented return value [GenePhaseScore *]
 *
 */
GenePhaseScore * Wise2_free_GenePhaseScore(GenePhaseScore * obj);
#define free_GenePhaseScore Wise2_free_GenePhaseScore


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
ProteinIntronList * Wise2_read_ProteinIntronList(FILE * ifp);
#define read_ProteinIntronList Wise2_read_ProteinIntronList
ProteinIntronList * Wise2_read_ProteinIntronList_from_filename(char * f);
#define read_ProteinIntronList_from_filename Wise2_read_ProteinIntronList_from_filename
PhasedProteinPara * Wise2_new_PhasedProteinPara_from_argv(int * argc,char ** argv);
#define new_PhasedProteinPara_from_argv Wise2_new_PhasedProteinPara_from_argv
void Wise2_show_help_PhasedProteinPara(FILE * ofp);
#define show_help_PhasedProteinPara Wise2_show_help_PhasedProteinPara
void Wise2_write_fasta_PhasedProtein(PhasedProtein * pp,FILE * ofp);
#define write_fasta_PhasedProtein Wise2_write_fasta_PhasedProtein
PhasedProtein * Wise2_read_fasta_PhasedProtein_file(char * file);
#define read_fasta_PhasedProtein_file Wise2_read_fasta_PhasedProtein_file
PhasedProtein * Wise2_read_fasta_PhasedProtein(FILE * ifp);
#define read_fasta_PhasedProtein Wise2_read_fasta_PhasedProtein
GenePhaseModel * Wise2_GenePhaseModel_from_ThreeStateModel(ThreeStateModel * tsm,CodonMapper * cm,RandomModel * rm,CompMat * mat,PhasedProteinPara * ppp);
#define GenePhaseModel_from_ThreeStateModel Wise2_GenePhaseModel_from_ThreeStateModel
GenePhaseScore * Wise2_GenePhaseScore_from_GenePhaseModel(GenePhaseModel * gpm);
#define GenePhaseScore_from_GenePhaseModel Wise2_GenePhaseScore_from_GenePhaseModel


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_ProteinIntronList(ProteinIntron ** list,int i,int j) ;
#define swap_ProteinIntronList Wise2_swap_ProteinIntronList
void Wise2_qsort_ProteinIntronList(ProteinIntron ** list,int left,int right,int (*comp)(ProteinIntron * ,ProteinIntron * ));
#define qsort_ProteinIntronList Wise2_qsort_ProteinIntronList
void Wise2_sort_ProteinIntronList(ProteinIntronList * obj,int (*comp)(ProteinIntron *, ProteinIntron *));
#define sort_ProteinIntronList Wise2_sort_ProteinIntronList
boolean Wise2_expand_ProteinIntronList(ProteinIntronList * obj,int len);
#define expand_ProteinIntronList Wise2_expand_ProteinIntronList
void Wise2_swap_GenePhaseModel(GenePhaseSeg ** list,int i,int j) ;
#define swap_GenePhaseModel Wise2_swap_GenePhaseModel
void Wise2_qsort_GenePhaseModel(GenePhaseSeg ** list,int left,int right,int (*comp)(GenePhaseSeg * ,GenePhaseSeg * ));
#define qsort_GenePhaseModel Wise2_qsort_GenePhaseModel
void Wise2_sort_GenePhaseModel(GenePhaseModel * obj,int (*comp)(GenePhaseSeg *, GenePhaseSeg *));
#define sort_GenePhaseModel Wise2_sort_GenePhaseModel
boolean Wise2_expand_GenePhaseModel(GenePhaseModel * obj,int len);
#define expand_GenePhaseModel Wise2_expand_GenePhaseModel
void Wise2_swap_GenePhaseScore(GenePhaseSegScore ** list,int i,int j) ;
#define swap_GenePhaseScore Wise2_swap_GenePhaseScore
void Wise2_qsort_GenePhaseScore(GenePhaseSegScore ** list,int left,int right,int (*comp)(GenePhaseSegScore * ,GenePhaseSegScore * ));
#define qsort_GenePhaseScore Wise2_qsort_GenePhaseScore
void Wise2_sort_GenePhaseScore(GenePhaseScore * obj,int (*comp)(GenePhaseSegScore *, GenePhaseSegScore *));
#define sort_GenePhaseScore Wise2_sort_GenePhaseScore
boolean Wise2_expand_GenePhaseScore(GenePhaseScore * obj,int len);
#define expand_GenePhaseScore Wise2_expand_GenePhaseScore

#ifdef _cplusplus
}
#endif

#endif
