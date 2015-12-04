#ifndef DYNAMITEfivestatemodelHEADERFILE
#define DYNAMITEfivestatemodelHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "threestatemodel.h"

enum fsm_trans_type {
  FSM_MATCH2MATCH  = 0,
  FSM_MATCH2INSERT,  
  FSM_MATCH2DELETE,
  FSM_MATCH2INBOUND,
  FSM_MATCH2OUTBOUND,
  FSM_INSERT2MATCH,
  FSM_INSERT2INSERT,  
  FSM_INSERT2DELETE,
  FSM_INSERT2INBOUND,
  FSM_INSERT2OUTBOUND,
  FSM_DELETE2MATCH,
  FSM_DELETE2INSERT,  
  FSM_DELETE2DELETE, 
  FSM_DELETE2INBOUND,
  FSM_DELETE2OUTBOUND,
  FSM_START2MATCH,    
  FSM_START2INSERT,   
  FSM_START2DELETE,   
  FSM_START2INBOUND,   
  FSM_MATCH2END,      
  FSM_INSERT2END,     
  FSM_DELETE2END,
  FSM_INBOUND2MATCH,
  FSM_INBOUND2INSERT,
  FSM_INBOUND2DELETE,
  FSM_INBOUND2END,
  FSM_INBOUND2INBOUND,
  FSM_OUTBOUND2MATCH,
  FSM_OUTBOUND2INSERT,
  FSM_OUTBOUND2DELETE,
  FSM_OUTBOUND2OUTBOUND,
  FSM_OUTBOUND2INBOUND,
  FSM_OUTBOUND2END,
  FSM_TRANSITION_LENGTH
};

#define FiveStateModelLISTLENGTH 128
#define FiveStateScoreLISTLENGTH 128

#define FiveStateFrameLISTLENGTH 32
#define FiveStateFrameSetLISTLENGTH 32

/*
 * Cheeky declaration to get around #include mess.
 */
ThreeStateModel * Wise2_HMMer2_read_ThreeStateModel(char * filename);

struct Wise2_FiveStateFrame {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    ThreeStateModel ** stack;    
    int len;/* len for above stack  */ 
    int maxlen; /* maxlen for above stack */ 
    } ;  
/* FiveStateFrame defined */ 
#ifndef DYNAMITE_DEFINED_FiveStateFrame
typedef struct Wise2_FiveStateFrame Wise2_FiveStateFrame;
#define FiveStateFrame Wise2_FiveStateFrame
#define DYNAMITE_DEFINED_FiveStateFrame
#endif


struct Wise2_FiveStateFrameSet {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    FiveStateFrame ** set;   
    int len;/* len for above set  */ 
    int maxlen; /* maxlen for above set */ 
    } ;  
/* FiveStateFrameSet defined */ 
#ifndef DYNAMITE_DEFINED_FiveStateFrameSet
typedef struct Wise2_FiveStateFrameSet Wise2_FiveStateFrameSet;
#define FiveStateFrameSet Wise2_FiveStateFrameSet
#define DYNAMITE_DEFINED_FiveStateFrameSet
#endif


struct Wise2_FiveStateUnit {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Probability match[ALPHABET_SIZE];    
    Probability insert[ALPHABET_SIZE];   
    Probability transition[FSM_TRANSITION_LENGTH];   
    char display_char;   
    } ;  
/* FiveStateUnit defined */ 
#ifndef DYNAMITE_DEFINED_FiveStateUnit
typedef struct Wise2_FiveStateUnit Wise2_FiveStateUnit;
#define FiveStateUnit Wise2_FiveStateUnit
#define DYNAMITE_DEFINED_FiveStateUnit
#endif


struct Wise2_FiveStateFrameIndicator {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    ThreeStateModel * tsm;   
    int start;   
    int end;     
    } ;  
/* FiveStateFrameIndicator defined */ 
#ifndef DYNAMITE_DEFINED_FiveStateFrameIndicator
typedef struct Wise2_FiveStateFrameIndicator Wise2_FiveStateFrameIndicator;
#define FiveStateFrameIndicator Wise2_FiveStateFrameIndicator
#define DYNAMITE_DEFINED_FiveStateFrameIndicator
#endif


struct Wise2_FiveStateModel {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * name;    /*  name of the model */ 
    FiveStateUnit ** unit;   
    int len;/* len for above unit  */ 
    int maxlen; /* maxlen for above unit */ 
    FiveStateFrameIndicator ** frame;    
    int tsm_len;/* len for above frame  */ 
    int tsm_maxlen; /* maxlen for above frame */ 
    } ;  
/* FiveStateModel defined */ 
#ifndef DYNAMITE_DEFINED_FiveStateModel
typedef struct Wise2_FiveStateModel Wise2_FiveStateModel;
#define FiveStateModel Wise2_FiveStateModel
#define DYNAMITE_DEFINED_FiveStateModel
#endif


struct Wise2_FiveStateScoreUnit {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Score match[ALPHABET_SIZE];  
    Score insert[ALPHABET_SIZE];     
    Score transition[FSM_TRANSITION_LENGTH];     
    char display_char;   
    } ;  
/* FiveStateScoreUnit defined */ 
#ifndef DYNAMITE_DEFINED_FiveStateScoreUnit
typedef struct Wise2_FiveStateScoreUnit Wise2_FiveStateScoreUnit;
#define FiveStateScoreUnit Wise2_FiveStateScoreUnit
#define DYNAMITE_DEFINED_FiveStateScoreUnit
#endif


struct Wise2_FiveStateScore {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * name;    /*  name of the model */ 
    FiveStateScoreUnit ** unit;  
    int len;/* len for above unit  */ 
    int maxlen; /* maxlen for above unit */ 
    } ;  
/* FiveStateScore defined */ 
#ifndef DYNAMITE_DEFINED_FiveStateScore
typedef struct Wise2_FiveStateScore Wise2_FiveStateScore;
#define FiveStateScore Wise2_FiveStateScore
#define DYNAMITE_DEFINED_FiveStateScore
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  pseudo_Protein_from_FiveStateModel(tsm)
 *
 * Descrip:    Makes a protein sequence out of the display characters.
 *             Not very useful!
 *
 *
 *
 * Arg:        tsm [UNKN ] Undocumented argument [FiveStateModel *]
 *
 * Return [UNKN ]  Undocumented return value [Protein *]
 *
 */
Protein * Wise2_pseudo_Protein_from_FiveStateModel(FiveStateModel * tsm);
#define pseudo_Protein_from_FiveStateModel Wise2_pseudo_Protein_from_FiveStateModel


/* Function:  add_FiveStateFrame(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [FiveStateFrame *]
 * Arg:        add [OWNER] Object to add to the list [ThreeStateModel *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_FiveStateFrame(FiveStateFrame * obj,ThreeStateModel * add);
#define add_FiveStateFrame Wise2_add_FiveStateFrame


/* Function:  flush_FiveStateFrame(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [FiveStateFrame *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_FiveStateFrame(FiveStateFrame * obj);
#define flush_FiveStateFrame Wise2_flush_FiveStateFrame


/* Function:  FiveStateFrame_alloc_std(void)
 *
 * Descrip:    Equivalent to FiveStateFrame_alloc_len(FiveStateFrameLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [FiveStateFrame *]
 *
 */
FiveStateFrame * Wise2_FiveStateFrame_alloc_std(void);
#define FiveStateFrame_alloc_std Wise2_FiveStateFrame_alloc_std


/* Function:  FiveStateFrame_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [FiveStateFrame *]
 *
 */
FiveStateFrame * Wise2_FiveStateFrame_alloc_len(int len);
#define FiveStateFrame_alloc_len Wise2_FiveStateFrame_alloc_len


/* Function:  hard_link_FiveStateFrame(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [FiveStateFrame *]
 *
 * Return [UNKN ]  Undocumented return value [FiveStateFrame *]
 *
 */
FiveStateFrame * Wise2_hard_link_FiveStateFrame(FiveStateFrame * obj);
#define hard_link_FiveStateFrame Wise2_hard_link_FiveStateFrame


/* Function:  FiveStateFrame_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [FiveStateFrame *]
 *
 */
FiveStateFrame * Wise2_FiveStateFrame_alloc(void);
#define FiveStateFrame_alloc Wise2_FiveStateFrame_alloc


/* Function:  free_FiveStateFrame(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [FiveStateFrame *]
 *
 * Return [UNKN ]  Undocumented return value [FiveStateFrame *]
 *
 */
FiveStateFrame * Wise2_free_FiveStateFrame(FiveStateFrame * obj);
#define free_FiveStateFrame Wise2_free_FiveStateFrame


/* Function:  add_FiveStateFrameSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [FiveStateFrameSet *]
 * Arg:        add [OWNER] Object to add to the list [FiveStateFrame *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_FiveStateFrameSet(FiveStateFrameSet * obj,FiveStateFrame * add);
#define add_FiveStateFrameSet Wise2_add_FiveStateFrameSet


/* Function:  flush_FiveStateFrameSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [FiveStateFrameSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_FiveStateFrameSet(FiveStateFrameSet * obj);
#define flush_FiveStateFrameSet Wise2_flush_FiveStateFrameSet


/* Function:  FiveStateFrameSet_alloc_std(void)
 *
 * Descrip:    Equivalent to FiveStateFrameSet_alloc_len(FiveStateFrameSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [FiveStateFrameSet *]
 *
 */
FiveStateFrameSet * Wise2_FiveStateFrameSet_alloc_std(void);
#define FiveStateFrameSet_alloc_std Wise2_FiveStateFrameSet_alloc_std


/* Function:  FiveStateFrameSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [FiveStateFrameSet *]
 *
 */
FiveStateFrameSet * Wise2_FiveStateFrameSet_alloc_len(int len);
#define FiveStateFrameSet_alloc_len Wise2_FiveStateFrameSet_alloc_len


/* Function:  hard_link_FiveStateFrameSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [FiveStateFrameSet *]
 *
 * Return [UNKN ]  Undocumented return value [FiveStateFrameSet *]
 *
 */
FiveStateFrameSet * Wise2_hard_link_FiveStateFrameSet(FiveStateFrameSet * obj);
#define hard_link_FiveStateFrameSet Wise2_hard_link_FiveStateFrameSet


/* Function:  FiveStateFrameSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [FiveStateFrameSet *]
 *
 */
FiveStateFrameSet * Wise2_FiveStateFrameSet_alloc(void);
#define FiveStateFrameSet_alloc Wise2_FiveStateFrameSet_alloc


/* Function:  free_FiveStateFrameSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [FiveStateFrameSet *]
 *
 * Return [UNKN ]  Undocumented return value [FiveStateFrameSet *]
 *
 */
FiveStateFrameSet * Wise2_free_FiveStateFrameSet(FiveStateFrameSet * obj);
#define free_FiveStateFrameSet Wise2_free_FiveStateFrameSet


/* Function:  hard_link_FiveStateUnit(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [FiveStateUnit *]
 *
 * Return [UNKN ]  Undocumented return value [FiveStateUnit *]
 *
 */
FiveStateUnit * Wise2_hard_link_FiveStateUnit(FiveStateUnit * obj);
#define hard_link_FiveStateUnit Wise2_hard_link_FiveStateUnit


/* Function:  FiveStateUnit_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [FiveStateUnit *]
 *
 */
FiveStateUnit * Wise2_FiveStateUnit_alloc(void);
#define FiveStateUnit_alloc Wise2_FiveStateUnit_alloc


/* Function:  free_FiveStateUnit(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [FiveStateUnit *]
 *
 * Return [UNKN ]  Undocumented return value [FiveStateUnit *]
 *
 */
FiveStateUnit * Wise2_free_FiveStateUnit(FiveStateUnit * obj);
#define free_FiveStateUnit Wise2_free_FiveStateUnit


/* Function:  hard_link_FiveStateFrameIndicator(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [FiveStateFrameIndicator *]
 *
 * Return [UNKN ]  Undocumented return value [FiveStateFrameIndicator *]
 *
 */
FiveStateFrameIndicator * Wise2_hard_link_FiveStateFrameIndicator(FiveStateFrameIndicator * obj);
#define hard_link_FiveStateFrameIndicator Wise2_hard_link_FiveStateFrameIndicator


/* Function:  FiveStateFrameIndicator_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [FiveStateFrameIndicator *]
 *
 */
FiveStateFrameIndicator * Wise2_FiveStateFrameIndicator_alloc(void);
#define FiveStateFrameIndicator_alloc Wise2_FiveStateFrameIndicator_alloc


/* Function:  free_FiveStateFrameIndicator(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [FiveStateFrameIndicator *]
 *
 * Return [UNKN ]  Undocumented return value [FiveStateFrameIndicator *]
 *
 */
FiveStateFrameIndicator * Wise2_free_FiveStateFrameIndicator(FiveStateFrameIndicator * obj);
#define free_FiveStateFrameIndicator Wise2_free_FiveStateFrameIndicator


/* Function:  add_FiveStateModel(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [FiveStateModel *]
 * Arg:        add [OWNER] Object to add to the list [FiveStateUnit *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_FiveStateModel(FiveStateModel * obj,FiveStateUnit * add);
#define add_FiveStateModel Wise2_add_FiveStateModel


/* Function:  flush_FiveStateModel(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [FiveStateModel *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_FiveStateModel(FiveStateModel * obj);
#define flush_FiveStateModel Wise2_flush_FiveStateModel


/* Function:  add_tsm_FiveStateModel(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [FiveStateModel *]
 * Arg:        add [OWNER] Object to add to the list [FiveStateFrameIndicator *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_tsm_FiveStateModel(FiveStateModel * obj,FiveStateFrameIndicator * add);
#define add_tsm_FiveStateModel Wise2_add_tsm_FiveStateModel


/* Function:  flush_tsm_FiveStateModel(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [FiveStateModel *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_tsm_FiveStateModel(FiveStateModel * obj);
#define flush_tsm_FiveStateModel Wise2_flush_tsm_FiveStateModel


/* Function:  FiveStateModel_alloc_std(void)
 *
 * Descrip:    Equivalent to FiveStateModel_alloc_len(FiveStateModelLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [FiveStateModel *]
 *
 */
FiveStateModel * Wise2_FiveStateModel_alloc_std(void);
#define FiveStateModel_alloc_std Wise2_FiveStateModel_alloc_std


/* Function:  FiveStateModel_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [FiveStateModel *]
 *
 */
FiveStateModel * Wise2_FiveStateModel_alloc_len(int len);
#define FiveStateModel_alloc_len Wise2_FiveStateModel_alloc_len


/* Function:  hard_link_FiveStateModel(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [FiveStateModel *]
 *
 * Return [UNKN ]  Undocumented return value [FiveStateModel *]
 *
 */
FiveStateModel * Wise2_hard_link_FiveStateModel(FiveStateModel * obj);
#define hard_link_FiveStateModel Wise2_hard_link_FiveStateModel


/* Function:  FiveStateModel_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [FiveStateModel *]
 *
 */
FiveStateModel * Wise2_FiveStateModel_alloc(void);
#define FiveStateModel_alloc Wise2_FiveStateModel_alloc


/* Function:  free_FiveStateModel(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [FiveStateModel *]
 *
 * Return [UNKN ]  Undocumented return value [FiveStateModel *]
 *
 */
FiveStateModel * Wise2_free_FiveStateModel(FiveStateModel * obj);
#define free_FiveStateModel Wise2_free_FiveStateModel


/* Function:  hard_link_FiveStateScoreUnit(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [FiveStateScoreUnit *]
 *
 * Return [UNKN ]  Undocumented return value [FiveStateScoreUnit *]
 *
 */
FiveStateScoreUnit * Wise2_hard_link_FiveStateScoreUnit(FiveStateScoreUnit * obj);
#define hard_link_FiveStateScoreUnit Wise2_hard_link_FiveStateScoreUnit


/* Function:  FiveStateScoreUnit_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [FiveStateScoreUnit *]
 *
 */
FiveStateScoreUnit * Wise2_FiveStateScoreUnit_alloc(void);
#define FiveStateScoreUnit_alloc Wise2_FiveStateScoreUnit_alloc


/* Function:  free_FiveStateScoreUnit(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [FiveStateScoreUnit *]
 *
 * Return [UNKN ]  Undocumented return value [FiveStateScoreUnit *]
 *
 */
FiveStateScoreUnit * Wise2_free_FiveStateScoreUnit(FiveStateScoreUnit * obj);
#define free_FiveStateScoreUnit Wise2_free_FiveStateScoreUnit


/* Function:  add_FiveStateScore(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [FiveStateScore *]
 * Arg:        add [OWNER] Object to add to the list [FiveStateScoreUnit *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_FiveStateScore(FiveStateScore * obj,FiveStateScoreUnit * add);
#define add_FiveStateScore Wise2_add_FiveStateScore


/* Function:  flush_FiveStateScore(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [FiveStateScore *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_FiveStateScore(FiveStateScore * obj);
#define flush_FiveStateScore Wise2_flush_FiveStateScore


/* Function:  FiveStateScore_alloc_std(void)
 *
 * Descrip:    Equivalent to FiveStateScore_alloc_len(FiveStateScoreLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [FiveStateScore *]
 *
 */
FiveStateScore * Wise2_FiveStateScore_alloc_std(void);
#define FiveStateScore_alloc_std Wise2_FiveStateScore_alloc_std


/* Function:  FiveStateScore_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [FiveStateScore *]
 *
 */
FiveStateScore * Wise2_FiveStateScore_alloc_len(int len);
#define FiveStateScore_alloc_len Wise2_FiveStateScore_alloc_len


/* Function:  hard_link_FiveStateScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [FiveStateScore *]
 *
 * Return [UNKN ]  Undocumented return value [FiveStateScore *]
 *
 */
FiveStateScore * Wise2_hard_link_FiveStateScore(FiveStateScore * obj);
#define hard_link_FiveStateScore Wise2_hard_link_FiveStateScore


/* Function:  FiveStateScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [FiveStateScore *]
 *
 */
FiveStateScore * Wise2_FiveStateScore_alloc(void);
#define FiveStateScore_alloc Wise2_FiveStateScore_alloc


/* Function:  free_FiveStateScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [FiveStateScore *]
 *
 * Return [UNKN ]  Undocumented return value [FiveStateScore *]
 *
 */
FiveStateScore * Wise2_free_FiveStateScore(FiveStateScore * obj);
#define free_FiveStateScore Wise2_free_FiveStateScore


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
void Wise2_dump_FiveStateModel(FiveStateModel * five,FILE * ofp);
#define dump_FiveStateModel Wise2_dump_FiveStateModel
void Wise2_dump_FiveStateScore(FiveStateScore * five,FILE * ofp);
#define dump_FiveStateScore Wise2_dump_FiveStateScore
FiveStateFrameSet * Wise2_read_FiveStateFrameSet_file(char * context,char * block_str);
#define read_FiveStateFrameSet_file Wise2_read_FiveStateFrameSet_file
FiveStateFrameSet * Wise2_read_FiveStateFrameSet(char * context,FILE * ifp);
#define read_FiveStateFrameSet Wise2_read_FiveStateFrameSet
void Wise2_blank_FiveStateUnit(FiveStateUnit * u);
#define blank_FiveStateUnit Wise2_blank_FiveStateUnit
FiveStateModel * Wise2_FiveStateModel_from_FiveStateFrameSet(FiveStateFrameSet * frame);
#define FiveStateModel_from_FiveStateFrameSet Wise2_FiveStateModel_from_FiveStateFrameSet
void Wise2_fold_RandomModel_into_FiveStateModel(FiveStateModel * fsm,RandomModel * rm);
#define fold_RandomModel_into_FiveStateModel Wise2_fold_RandomModel_into_FiveStateModel
FiveStateModel * Wise2_FiveStateModel_from_flat_ThreeStateModel(ThreeStateModel * tsm);
#define FiveStateModel_from_flat_ThreeStateModel Wise2_FiveStateModel_from_flat_ThreeStateModel
FiveStateScore * Wise2_FiveStateScore_from_FiveStateModel(FiveStateModel * fsm);
#define FiveStateScore_from_FiveStateModel Wise2_FiveStateScore_from_FiveStateModel
FiveStateScoreUnit * Wise2_FiveStateScoreUnit_from_FiveStateUnit(FiveStateUnit * prev,FiveStateUnit * u);
#define FiveStateScoreUnit_from_FiveStateUnit Wise2_FiveStateScoreUnit_from_FiveStateUnit


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_FiveStateFrame(ThreeStateModel ** list,int i,int j) ;
#define swap_FiveStateFrame Wise2_swap_FiveStateFrame
void Wise2_qsort_FiveStateFrame(ThreeStateModel ** list,int left,int right,int (*comp)(ThreeStateModel * ,ThreeStateModel * ));
#define qsort_FiveStateFrame Wise2_qsort_FiveStateFrame
void Wise2_sort_FiveStateFrame(FiveStateFrame * obj,int (*comp)(ThreeStateModel *, ThreeStateModel *));
#define sort_FiveStateFrame Wise2_sort_FiveStateFrame
boolean Wise2_expand_FiveStateFrame(FiveStateFrame * obj,int len);
#define expand_FiveStateFrame Wise2_expand_FiveStateFrame
void Wise2_swap_FiveStateFrameSet(FiveStateFrame ** list,int i,int j) ;
#define swap_FiveStateFrameSet Wise2_swap_FiveStateFrameSet
void Wise2_qsort_FiveStateFrameSet(FiveStateFrame ** list,int left,int right,int (*comp)(FiveStateFrame * ,FiveStateFrame * ));
#define qsort_FiveStateFrameSet Wise2_qsort_FiveStateFrameSet
void Wise2_sort_FiveStateFrameSet(FiveStateFrameSet * obj,int (*comp)(FiveStateFrame *, FiveStateFrame *));
#define sort_FiveStateFrameSet Wise2_sort_FiveStateFrameSet
boolean Wise2_expand_FiveStateFrameSet(FiveStateFrameSet * obj,int len);
#define expand_FiveStateFrameSet Wise2_expand_FiveStateFrameSet
void Wise2_swap_FiveStateModel(FiveStateUnit ** list,int i,int j) ;
#define swap_FiveStateModel Wise2_swap_FiveStateModel
void Wise2_qsort_FiveStateModel(FiveStateUnit ** list,int left,int right,int (*comp)(FiveStateUnit * ,FiveStateUnit * ));
#define qsort_FiveStateModel Wise2_qsort_FiveStateModel
void Wise2_sort_FiveStateModel(FiveStateModel * obj,int (*comp)(FiveStateUnit *, FiveStateUnit *));
#define sort_FiveStateModel Wise2_sort_FiveStateModel
boolean Wise2_expand_FiveStateModel(FiveStateModel * obj,int len);
#define expand_FiveStateModel Wise2_expand_FiveStateModel
void Wise2_swap_tsm_FiveStateModel(FiveStateFrameIndicator ** list,int i,int j) ;
#define swap_tsm_FiveStateModel Wise2_swap_tsm_FiveStateModel
void Wise2_qsort_tsm_FiveStateModel(FiveStateFrameIndicator ** list,int left,int right,int (*comp)(FiveStateFrameIndicator * ,FiveStateFrameIndicator * ));
#define qsort_tsm_FiveStateModel Wise2_qsort_tsm_FiveStateModel
void Wise2_sort_tsm_FiveStateModel(FiveStateModel * obj,int (*comp)(FiveStateFrameIndicator *, FiveStateFrameIndicator *));
#define sort_tsm_FiveStateModel Wise2_sort_tsm_FiveStateModel
boolean Wise2_expand_tsm_FiveStateModel(FiveStateModel * obj,int len);
#define expand_tsm_FiveStateModel Wise2_expand_tsm_FiveStateModel
void Wise2_swap_FiveStateScore(FiveStateScoreUnit ** list,int i,int j) ;
#define swap_FiveStateScore Wise2_swap_FiveStateScore
void Wise2_qsort_FiveStateScore(FiveStateScoreUnit ** list,int left,int right,int (*comp)(FiveStateScoreUnit * ,FiveStateScoreUnit * ));
#define qsort_FiveStateScore Wise2_qsort_FiveStateScore
void Wise2_sort_FiveStateScore(FiveStateScore * obj,int (*comp)(FiveStateScoreUnit *, FiveStateScoreUnit *));
#define sort_FiveStateScore Wise2_sort_FiveStateScore
boolean Wise2_expand_FiveStateScore(FiveStateScore * obj,int len);
#define expand_FiveStateScore Wise2_expand_FiveStateScore

#ifdef _cplusplus
}
#endif

#endif
