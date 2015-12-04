#ifndef DYNAMITEmatrixdebugHEADERFILE
#define DYNAMITEmatrixdebugHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "basematrix.h"
#include "packaln.h"

#define DebugStateLISTLENGTH    40
#define DebugMatrixLISTLENGTH  40

typedef enum MatrixDebugBreakPoint {
	MDBP_Cursor = 321,
	MDBP_Overflow,
	MDBP_Underflow,
	MDBP_NoBreakPoint
} MatrixDebugBreakPoint_type;


struct Wise2_DebugTransition {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * fromstate;    
    int from_state_num;  
    void (*show_transition)(void * matrix,int i,int j,FILE * ofp);   
    } ;  
/* DebugTransition defined */ 
#ifndef DYNAMITE_DEFINED_DebugTransition
typedef struct Wise2_DebugTransition Wise2_DebugTransition;
#define DebugTransition Wise2_DebugTransition
#define DYNAMITE_DEFINED_DebugTransition
#endif


struct Wise2_DebugState {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * statename;    
    int state_num;   
    DebugTransition ** trans;    
    int len;/* len for above trans  */ 
    int maxlen; /* maxlen for above trans */ 
    void (*show_state)(void * matrix,int i,int j,FILE *ofp); 
    } ;  
/* DebugState defined */ 
#ifndef DYNAMITE_DEFINED_DebugState
typedef struct Wise2_DebugState Wise2_DebugState;
#define DebugState Wise2_DebugState
#define DYNAMITE_DEFINED_DebugState
#endif


struct Wise2_DebugBreakPoint {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    MatrixDebugBreakPoint_type type;     
    int posi;    
    int posj;    
    int overflow;    
    int underflow;   
    } ;  
/* DebugBreakPoint defined */ 
#ifndef DYNAMITE_DEFINED_DebugBreakPoint
typedef struct Wise2_DebugBreakPoint Wise2_DebugBreakPoint;
#define DebugBreakPoint Wise2_DebugBreakPoint
#define DYNAMITE_DEFINED_DebugBreakPoint
#endif


struct Wise2_DebugMatrix {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    void * matrix;   
    int currenti;    
    int currentj;    
    int max_score;   
    int min_score;   
    int max_score_i;     
    int max_score_j;     
    int max_score_cell;  
    DebugState ** state;     
    int len;/* len for above state  */ 
    int maxlen; /* maxlen for above state */ 
    DebugBreakPoint ** point;    
    int bp_len; /* len for above point  */ 
    int bp_maxlen;  /* maxlen for above point */ 
    boolean reset;   
    FILE * in;   
    FILE * out;  
    void (*show_cell)(void * matrix,int i,int j,FILE *ofp);  
    } ;  
/* DebugMatrix defined */ 
#ifndef DYNAMITE_DEFINED_DebugMatrix
typedef struct Wise2_DebugMatrix Wise2_DebugMatrix;
#define DebugMatrix Wise2_DebugMatrix
#define DYNAMITE_DEFINED_DebugMatrix
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  show_PackAln_Debug(de,pal,ofp)
 *
 * Descrip:    Shows PackAln in debug context
 *
 *
 * Arg:         de [UNKN ] Undocumented argument [DebugMatrix *]
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_show_PackAln_Debug(DebugMatrix * de,PackAln * pal,FILE * ofp);
#define show_PackAln_Debug Wise2_show_PackAln_Debug


/* Function:  show_DebugMatrix(de,buffer)
 *
 * Descrip:    Show function cell state transition
 *
 *
 * Arg:            de [UNKN ] Undocumented argument [DebugMatrix *]
 * Arg:        buffer [UNKN ] Undocumented argument [char *]
 *
 */
void Wise2_show_DebugMatrix(DebugMatrix * de,char * buffer);
#define show_DebugMatrix Wise2_show_DebugMatrix


/* Function:  find_DebugState(de,name)
 *
 * Descrip:    Finds a DebugState 
 *
 *
 * Arg:          de [UNKN ] Undocumented argument [DebugMatrix *]
 * Arg:        name [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [DebugState *]
 *
 */
DebugState * Wise2_find_DebugState(DebugMatrix * de,char * name);
#define find_DebugState Wise2_find_DebugState


/* Function:  find_DebugTransition(state,name)
 *
 * Descrip:    Finds a DebugTransition
 *
 *
 * Arg:        state [UNKN ] Undocumented argument [DebugState *]
 * Arg:         name [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [DebugTransition *]
 *
 */
DebugTransition * Wise2_find_DebugTransition(DebugState * state,char * name);
#define find_DebugTransition Wise2_find_DebugTransition


/* Function:  user_DebugMatrix(de)
 *
 * Descrip:    Main function to talk to user
 *
 *
 * Arg:        de [UNKN ] Undocumented argument [DebugMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_user_DebugMatrix(DebugMatrix * de);
#define user_DebugMatrix Wise2_user_DebugMatrix


/* Function:  should_break_DebugMatrix(de)
 *
 * Descrip:    Indicates whether we should break at this point. Assummes
 *             the de datastructure has been updated correctly
 *
 *
 * Arg:        de [UNKN ] Undocumented argument [DebugMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [MatrixDebugBreakPoint_type]
 *
 */
MatrixDebugBreakPoint_type Wise2_should_break_DebugMatrix(DebugMatrix * de);
#define should_break_DebugMatrix Wise2_should_break_DebugMatrix


/* Function:  std_DebugMatrix(void)
 *
 * Descrip:    Builds a "standard" Debug matrix
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DebugMatrix *]
 *
 */
DebugMatrix * Wise2_std_DebugMatrix(void);
#define std_DebugMatrix Wise2_std_DebugMatrix


/* Function:  hard_link_DebugTransition(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DebugTransition *]
 *
 * Return [UNKN ]  Undocumented return value [DebugTransition *]
 *
 */
DebugTransition * Wise2_hard_link_DebugTransition(DebugTransition * obj);
#define hard_link_DebugTransition Wise2_hard_link_DebugTransition


/* Function:  DebugTransition_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DebugTransition *]
 *
 */
DebugTransition * Wise2_DebugTransition_alloc(void);
#define DebugTransition_alloc Wise2_DebugTransition_alloc


/* Function:  free_DebugTransition(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DebugTransition *]
 *
 * Return [UNKN ]  Undocumented return value [DebugTransition *]
 *
 */
DebugTransition * Wise2_free_DebugTransition(DebugTransition * obj);
#define free_DebugTransition Wise2_free_DebugTransition


/* Function:  add_DebugState(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [DebugState *]
 * Arg:        add [OWNER] Object to add to the list [DebugTransition *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_DebugState(DebugState * obj,DebugTransition * add);
#define add_DebugState Wise2_add_DebugState


/* Function:  flush_DebugState(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [DebugState *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_DebugState(DebugState * obj);
#define flush_DebugState Wise2_flush_DebugState


/* Function:  DebugState_alloc_std(void)
 *
 * Descrip:    Equivalent to DebugState_alloc_len(DebugStateLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DebugState *]
 *
 */
DebugState * Wise2_DebugState_alloc_std(void);
#define DebugState_alloc_std Wise2_DebugState_alloc_std


/* Function:  DebugState_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [DebugState *]
 *
 */
DebugState * Wise2_DebugState_alloc_len(int len);
#define DebugState_alloc_len Wise2_DebugState_alloc_len


/* Function:  hard_link_DebugState(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DebugState *]
 *
 * Return [UNKN ]  Undocumented return value [DebugState *]
 *
 */
DebugState * Wise2_hard_link_DebugState(DebugState * obj);
#define hard_link_DebugState Wise2_hard_link_DebugState


/* Function:  DebugState_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DebugState *]
 *
 */
DebugState * Wise2_DebugState_alloc(void);
#define DebugState_alloc Wise2_DebugState_alloc


/* Function:  free_DebugState(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DebugState *]
 *
 * Return [UNKN ]  Undocumented return value [DebugState *]
 *
 */
DebugState * Wise2_free_DebugState(DebugState * obj);
#define free_DebugState Wise2_free_DebugState


/* Function:  hard_link_DebugBreakPoint(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DebugBreakPoint *]
 *
 * Return [UNKN ]  Undocumented return value [DebugBreakPoint *]
 *
 */
DebugBreakPoint * Wise2_hard_link_DebugBreakPoint(DebugBreakPoint * obj);
#define hard_link_DebugBreakPoint Wise2_hard_link_DebugBreakPoint


/* Function:  DebugBreakPoint_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DebugBreakPoint *]
 *
 */
DebugBreakPoint * Wise2_DebugBreakPoint_alloc(void);
#define DebugBreakPoint_alloc Wise2_DebugBreakPoint_alloc


/* Function:  free_DebugBreakPoint(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DebugBreakPoint *]
 *
 * Return [UNKN ]  Undocumented return value [DebugBreakPoint *]
 *
 */
DebugBreakPoint * Wise2_free_DebugBreakPoint(DebugBreakPoint * obj);
#define free_DebugBreakPoint Wise2_free_DebugBreakPoint


/* Function:  add_DebugMatrix(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [DebugMatrix *]
 * Arg:        add [OWNER] Object to add to the list [DebugState *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_DebugMatrix(DebugMatrix * obj,DebugState * add);
#define add_DebugMatrix Wise2_add_DebugMatrix


/* Function:  flush_DebugMatrix(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [DebugMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_DebugMatrix(DebugMatrix * obj);
#define flush_DebugMatrix Wise2_flush_DebugMatrix


/* Function:  add_bp_DebugMatrix(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [DebugMatrix *]
 * Arg:        add [OWNER] Object to add to the list [DebugBreakPoint *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_bp_DebugMatrix(DebugMatrix * obj,DebugBreakPoint * add);
#define add_bp_DebugMatrix Wise2_add_bp_DebugMatrix


/* Function:  flush_bp_DebugMatrix(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [DebugMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_bp_DebugMatrix(DebugMatrix * obj);
#define flush_bp_DebugMatrix Wise2_flush_bp_DebugMatrix


/* Function:  DebugMatrix_alloc_std(void)
 *
 * Descrip:    Equivalent to DebugMatrix_alloc_len(DebugMatrixLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DebugMatrix *]
 *
 */
DebugMatrix * Wise2_DebugMatrix_alloc_std(void);
#define DebugMatrix_alloc_std Wise2_DebugMatrix_alloc_std


/* Function:  DebugMatrix_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [DebugMatrix *]
 *
 */
DebugMatrix * Wise2_DebugMatrix_alloc_len(int len);
#define DebugMatrix_alloc_len Wise2_DebugMatrix_alloc_len


/* Function:  hard_link_DebugMatrix(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DebugMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [DebugMatrix *]
 *
 */
DebugMatrix * Wise2_hard_link_DebugMatrix(DebugMatrix * obj);
#define hard_link_DebugMatrix Wise2_hard_link_DebugMatrix


/* Function:  DebugMatrix_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DebugMatrix *]
 *
 */
DebugMatrix * Wise2_DebugMatrix_alloc(void);
#define DebugMatrix_alloc Wise2_DebugMatrix_alloc


/* Function:  free_DebugMatrix(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DebugMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [DebugMatrix *]
 *
 */
DebugMatrix * Wise2_free_DebugMatrix(DebugMatrix * obj);
#define free_DebugMatrix Wise2_free_DebugMatrix


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_DebugState(DebugTransition ** list,int i,int j) ;
#define swap_DebugState Wise2_swap_DebugState
void Wise2_qsort_DebugState(DebugTransition ** list,int left,int right,int (*comp)(DebugTransition * ,DebugTransition * ));
#define qsort_DebugState Wise2_qsort_DebugState
void Wise2_sort_DebugState(DebugState * obj,int (*comp)(DebugTransition *, DebugTransition *));
#define sort_DebugState Wise2_sort_DebugState
boolean Wise2_expand_DebugState(DebugState * obj,int len);
#define expand_DebugState Wise2_expand_DebugState
void Wise2_swap_DebugMatrix(DebugState ** list,int i,int j) ;
#define swap_DebugMatrix Wise2_swap_DebugMatrix
void Wise2_qsort_DebugMatrix(DebugState ** list,int left,int right,int (*comp)(DebugState * ,DebugState * ));
#define qsort_DebugMatrix Wise2_qsort_DebugMatrix
void Wise2_sort_DebugMatrix(DebugMatrix * obj,int (*comp)(DebugState *, DebugState *));
#define sort_DebugMatrix Wise2_sort_DebugMatrix
boolean Wise2_expand_DebugMatrix(DebugMatrix * obj,int len);
#define expand_DebugMatrix Wise2_expand_DebugMatrix
void Wise2_swap_bp_DebugMatrix(DebugBreakPoint ** list,int i,int j) ;
#define swap_bp_DebugMatrix Wise2_swap_bp_DebugMatrix
void Wise2_qsort_bp_DebugMatrix(DebugBreakPoint ** list,int left,int right,int (*comp)(DebugBreakPoint * ,DebugBreakPoint * ));
#define qsort_bp_DebugMatrix Wise2_qsort_bp_DebugMatrix
void Wise2_sort_bp_DebugMatrix(DebugMatrix * obj,int (*comp)(DebugBreakPoint *, DebugBreakPoint *));
#define sort_bp_DebugMatrix Wise2_sort_bp_DebugMatrix
boolean Wise2_expand_bp_DebugMatrix(DebugMatrix * obj,int len);
#define expand_bp_DebugMatrix Wise2_expand_bp_DebugMatrix

#ifdef _cplusplus
}
#endif

#endif
