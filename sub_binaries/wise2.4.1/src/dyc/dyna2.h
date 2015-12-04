#ifndef DYNAMITEdyna2HEADERFILE
#define DYNAMITEdyna2HEADERFILE
#ifdef _cplusplus
extern "C" {
#endif

#include "wisebase.h"
#include "wisec.h"
#include "exprtree.h"
#include "type.h"

#define CellSourceLISTLENGTH 64
#define CellStateLISTLENGTH  64
#define GenericMatrixLISTLENGTH 64
#define CellExprLISTLENGTH 64
#define CellSignatureSetLISTLENGTH 64


#define STANDARDGENERICMATRIX   231
#define SEARCHGENERICMATRIX     232 
#define VARIABLEGENERICMATRIX   234


enum calc_unit {
  GENERICMATRIX_IDEPCALCUNIT = 214,
  GENERICMATRIX_JDEPCALCUNIT,
  GENERICMATRIX_IJDEPCALCUNIT
} ;


/*
 * Ok, gets a bit weird: if position is less
 * than 21 then it has the bit number from the
 * bit positions so it can be warned properely
 * Look at source_bit2position function.
 */


enum source_positon {
  SOURCE_POS_ALL     = 16,
  SOURCE_POS_TOPLEFT = 32,
  SOURCE_POS_TOP     = 64,
  SOURCE_POS_LEFT    = 128,
  SOURCE_POS_BOTTOM  = 256,
  SOURCE_POS_RIGHT   = 512,
  SOURCE_POS_BOTTOMRIGHT = 1024 };

#define SOURCE_TOP_BIT    1
#define SOURCE_LEFT_BIT   2
#define SOURCE_BOTTOM_BIT 4
#define SOURCE_RIGHT_BIT  8

struct CellSignature {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int offi;    
    int offj;    
    } ;  
/* CellSignature defined */ 
#ifndef DYNAMITE_DEFINED_CellSignature
typedef struct CellSignature CellSignature;
#define DYNAMITE_DEFINED_CellSignature
#endif


struct CellSignatureSet {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    CellSignature ** sig;    
    int len;/* len for above sig  */ 
    int maxlen; /* maxlen for above sig */ 
    } ;  
/* CellSignatureSet defined */ 
#ifndef DYNAMITE_DEFINED_CellSignatureSet
typedef struct CellSignatureSet CellSignatureSet;
#define DYNAMITE_DEFINED_CellSignatureSet
#endif


struct CellSource {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * state_source;     
    int offi;    
    int offj;    
    char * calc_expr;   /*  going to replace CellExpr list possibly */ 
    char * source_expr;  
    ExprTree * etr;  
    boolean isspecial;   
    char * query_label;  
    char * target_label;     
    int  position;   
    int trans_no;   /*  unique number for this transition */ 
    int from_state_no;   
    } ;  
/* CellSource defined */ 
#ifndef DYNAMITE_DEFINED_CellSource
typedef struct CellSource CellSource;
#define DYNAMITE_DEFINED_CellSource
#endif


struct CellState {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * name;     
    char * def_score;    
    char * calc_expr;    
    char * source_expr;  
    int offi;    
    int offj;    
    boolean is_special_i;    
    boolean is_special_j;    
    boolean is_end;  
    boolean is_start;    
    boolean specialtospecial;    
    CellSource ** source;    
    int len;/* len for above source  */ 
    int maxlen; /* maxlen for above source */ 
    char * query_char;   
    char * target_char;  
    int footprint_start;     
    int footprint_end;   
    char * query_label;  
    char * target_label;     
    int position;    
    ExprTree * etr;  
    int state_number;    
    } ;  
/* CellState defined */ 
#ifndef DYNAMITE_DEFINED_CellState
typedef struct CellState CellState;
#define DYNAMITE_DEFINED_CellState
#endif


struct CollapsableLabel {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * query;    
    char * target;   
    } ;  
/* CollapsableLabel defined */ 
#ifndef DYNAMITE_DEFINED_CollapsableLabel
typedef struct CollapsableLabel CollapsableLabel;
#define DYNAMITE_DEFINED_CollapsableLabel
#endif


struct ExternVariable {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * name;     
    char * type;     
    } ;  
/* ExternVariable defined */ 
#ifndef DYNAMITE_DEFINED_ExternVariable
typedef struct ExternVariable ExternVariable;
#define DYNAMITE_DEFINED_ExternVariable
#endif


struct GenericMatrix {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * name;     
    int type;    
    CellState ** state;  
    int len;/* len for above state  */ 
    int maxlen; /* maxlen for above state */ 
    CellState ** special;    
    int spec_len;   /* len for above special  */ 
    int spec_maxlen;/* maxlen for above special */ 
    StructElement * query;   
    char * query_name;   
    char * query_len;    
    Type * qtype;    
    StructElement * target;  
    char * target_name;  
    char * target_len;   
    Type * ttype;    
    StructElement  ** resource;  
    int res_len;/* len for above resource  */ 
    int res_maxlen; /* maxlen for above resource */ 
    CollapsableLabel ** cal;     
    int cal_len;/* len for above cal  */ 
    int cal_maxlen; /* maxlen for above cal */ 
    ExternVariable   ** ev;  
    int ev_len; /* len for above ev  */ 
    int ev_maxlen;  /* maxlen for above ev */ 
    char * defscore_all_states;  
    int window_i;    
    int window_j;    
    int footprint;   
    boolean cansearch;   
    boolean canlabel;    
    boolean specialtospecial;    
    StructHolder * sh;  /*  where structure for matrix is placed */ 
    char * calcfunc;     
    Scope * sc;  
    MethodTypeSet * mts;     
    } ;  
/* GenericMatrix defined */ 
#ifndef DYNAMITE_DEFINED_GenericMatrix
typedef struct GenericMatrix GenericMatrix;
#define DYNAMITE_DEFINED_GenericMatrix
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  prepare_matrix(gm,mts,dycw,failing_errors)
 *
 * Descrip:    main function to check GenericMatrix onced parsed
 *
 *             checks
 *               state defaults
 *               state/source cross references
 *               labels
 *               calc epxressions
 *               types and type migration
 *               calc parsing
 *
 *
 * Arg:                    gm [RW   ] GenericMatrix to be checked [GenericMatrix *]
 * Arg:                   mts [READ ] Type and Method Scope [MethodTypeSet *]
 * Arg:                  dycw [UNKN ] Undocumented argument [DycWarning *]
 * Arg:        failing_errors [READ ] Calc line parser on which errors fail [ParseError]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean prepare_matrix(GenericMatrix * gm,MethodTypeSet * mts,DycWarning * dycw,ParseError failing_errors);


/* Function:  assign_source_no(gm)
 *
 * Descrip:    Adds a unique transition number for CellSource
 *
 *
 * Arg:        gm [UNKN ] Undocumented argument [GenericMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean assign_source_no(GenericMatrix * gm);


/* Function:  check_start_end(gm)
 *
 * Descrip:    checks we have a start + end (and only 1 each!)
 *             and sets start's defscore to 0
 *
 *
 * Arg:        gm [UNKN ] Undocumented argument [GenericMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean check_start_end(GenericMatrix * gm);


/* Function:  check_source_positions(gm)
 *
 * Descrip:    checks the top/bottom/left/right source
 *             positions
 *
 *
 * Arg:        gm [UNKN ] Undocumented argument [GenericMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean check_source_positions(GenericMatrix * gm);


/* Function:  read_ExternVariable_line(line)
 *
 * Descrip:    reads line like extern name="xxx" type="xxx"
 *
 *
 * Arg:        line [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [ExternVariable *]
 *
 */
ExternVariable * read_ExternVariable_line(char * line);


/* Function:  hard_link_CellSignature(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [CellSignature *]
 *
 * Return [UNKN ]  Undocumented return value [CellSignature *]
 *
 */
CellSignature * hard_link_CellSignature(CellSignature * obj);


/* Function:  CellSignature_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [CellSignature *]
 *
 */
CellSignature * CellSignature_alloc(void);


/* Function:  free_CellSignature(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [CellSignature *]
 *
 * Return [UNKN ]  Undocumented return value [CellSignature *]
 *
 */
CellSignature * free_CellSignature(CellSignature * obj);


/* Function:  add_CellSignatureSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [CellSignatureSet *]
 * Arg:        add [OWNER] Object to add to the list [CellSignature *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean add_CellSignatureSet(CellSignatureSet * obj,CellSignature * add);


/* Function:  flush_CellSignatureSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [CellSignatureSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_CellSignatureSet(CellSignatureSet * obj);


/* Function:  CellSignatureSet_alloc_std(void)
 *
 * Descrip:    Equivalent to CellSignatureSet_alloc_len(CellSignatureSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [CellSignatureSet *]
 *
 */
CellSignatureSet * CellSignatureSet_alloc_std(void);


/* Function:  CellSignatureSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [CellSignatureSet *]
 *
 */
CellSignatureSet * CellSignatureSet_alloc_len(int len);


/* Function:  hard_link_CellSignatureSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [CellSignatureSet *]
 *
 * Return [UNKN ]  Undocumented return value [CellSignatureSet *]
 *
 */
CellSignatureSet * hard_link_CellSignatureSet(CellSignatureSet * obj);


/* Function:  CellSignatureSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [CellSignatureSet *]
 *
 */
CellSignatureSet * CellSignatureSet_alloc(void);


/* Function:  free_CellSignatureSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [CellSignatureSet *]
 *
 * Return [UNKN ]  Undocumented return value [CellSignatureSet *]
 *
 */
CellSignatureSet * free_CellSignatureSet(CellSignatureSet * obj);


/* Function:  hard_link_CellSource(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [CellSource *]
 *
 * Return [UNKN ]  Undocumented return value [CellSource *]
 *
 */
CellSource * hard_link_CellSource(CellSource * obj);


/* Function:  CellSource_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [CellSource *]
 *
 */
CellSource * CellSource_alloc(void);


/* Function:  free_CellSource(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [CellSource *]
 *
 * Return [UNKN ]  Undocumented return value [CellSource *]
 *
 */
CellSource * free_CellSource(CellSource * obj);


/* Function:  add_CellState(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [CellState *]
 * Arg:        add [OWNER] Object to add to the list [CellSource *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean add_CellState(CellState * obj,CellSource * add);


/* Function:  flush_CellState(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [CellState *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_CellState(CellState * obj);


/* Function:  CellState_alloc_std(void)
 *
 * Descrip:    Equivalent to CellState_alloc_len(CellStateLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [CellState *]
 *
 */
CellState * CellState_alloc_std(void);


/* Function:  CellState_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [CellState *]
 *
 */
CellState * CellState_alloc_len(int len);


/* Function:  hard_link_CellState(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [CellState *]
 *
 * Return [UNKN ]  Undocumented return value [CellState *]
 *
 */
CellState * hard_link_CellState(CellState * obj);


/* Function:  CellState_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [CellState *]
 *
 */
CellState * CellState_alloc(void);


/* Function:  free_CellState(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [CellState *]
 *
 * Return [UNKN ]  Undocumented return value [CellState *]
 *
 */
CellState * free_CellState(CellState * obj);


/* Function:  hard_link_CollapsableLabel(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [CollapsableLabel *]
 *
 * Return [UNKN ]  Undocumented return value [CollapsableLabel *]
 *
 */
CollapsableLabel * hard_link_CollapsableLabel(CollapsableLabel * obj);


/* Function:  CollapsableLabel_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [CollapsableLabel *]
 *
 */
CollapsableLabel * CollapsableLabel_alloc(void);


/* Function:  free_CollapsableLabel(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [CollapsableLabel *]
 *
 * Return [UNKN ]  Undocumented return value [CollapsableLabel *]
 *
 */
CollapsableLabel * free_CollapsableLabel(CollapsableLabel * obj);


/* Function:  hard_link_ExternVariable(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ExternVariable *]
 *
 * Return [UNKN ]  Undocumented return value [ExternVariable *]
 *
 */
ExternVariable * hard_link_ExternVariable(ExternVariable * obj);


/* Function:  ExternVariable_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ExternVariable *]
 *
 */
ExternVariable * ExternVariable_alloc(void);


/* Function:  free_ExternVariable(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ExternVariable *]
 *
 * Return [UNKN ]  Undocumented return value [ExternVariable *]
 *
 */
ExternVariable * free_ExternVariable(ExternVariable * obj);


/* Function:  add_GenericMatrix(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GenericMatrix *]
 * Arg:        add [OWNER] Object to add to the list [CellState *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean add_GenericMatrix(GenericMatrix * obj,CellState * add);


/* Function:  flush_GenericMatrix(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [GenericMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_GenericMatrix(GenericMatrix * obj);


/* Function:  add_spec_GenericMatrix(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GenericMatrix *]
 * Arg:        add [OWNER] Object to add to the list [CellState *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean add_spec_GenericMatrix(GenericMatrix * obj,CellState * add);


/* Function:  flush_spec_GenericMatrix(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [GenericMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_spec_GenericMatrix(GenericMatrix * obj);


/* Function:  add_res_GenericMatrix(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GenericMatrix *]
 * Arg:        add [OWNER] Object to add to the list [StructElement  *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean add_res_GenericMatrix(GenericMatrix * obj,StructElement  * add);


/* Function:  flush_res_GenericMatrix(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [GenericMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_res_GenericMatrix(GenericMatrix * obj);


/* Function:  add_cal_GenericMatrix(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GenericMatrix *]
 * Arg:        add [OWNER] Object to add to the list [CollapsableLabel *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean add_cal_GenericMatrix(GenericMatrix * obj,CollapsableLabel * add);


/* Function:  flush_cal_GenericMatrix(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [GenericMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_cal_GenericMatrix(GenericMatrix * obj);


/* Function:  add_ev_GenericMatrix(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GenericMatrix *]
 * Arg:        add [OWNER] Object to add to the list [ExternVariable   *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean add_ev_GenericMatrix(GenericMatrix * obj,ExternVariable   * add);


/* Function:  flush_ev_GenericMatrix(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [GenericMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_ev_GenericMatrix(GenericMatrix * obj);


/* Function:  GenericMatrix_alloc_std(void)
 *
 * Descrip:    Equivalent to GenericMatrix_alloc_len(GenericMatrixLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GenericMatrix *]
 *
 */
GenericMatrix * GenericMatrix_alloc_std(void);


/* Function:  GenericMatrix_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [GenericMatrix *]
 *
 */
GenericMatrix * GenericMatrix_alloc_len(int len);


/* Function:  hard_link_GenericMatrix(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GenericMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [GenericMatrix *]
 *
 */
GenericMatrix * hard_link_GenericMatrix(GenericMatrix * obj);


/* Function:  GenericMatrix_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GenericMatrix *]
 *
 */
GenericMatrix * GenericMatrix_alloc(void);


/* Function:  free_GenericMatrix(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GenericMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [GenericMatrix *]
 *
 */
GenericMatrix * free_GenericMatrix(GenericMatrix * obj);


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
CellSignatureSet * CellSignatureSet_from_GenericMatrix(GenericMatrix * gm);
boolean can_do_threads(GenericMatrix * gm);
CellState * CellState_from_str(GenericMatrix * gm,char * str);
boolean check_cell_refs(GenericMatrix * gm);
boolean prepare_labels(GenericMatrix * gm);
boolean handle_names(GenericMatrix * gm);
boolean handle_type_migration(GenericMatrix * gm,MethodTypeSet * mts);
boolean perculate_state_defaults(GenericMatrix * gm);
boolean make_StructHolder_for_GenericMatrix(GenericMatrix * gm,MethodTypeSet * mts);
StructElement * StructElement_for_GenericMatrix_type(char * name,char * element);
char * length_string_from_GenericMatrix_type(char * element);
boolean can_interpret_type(char * type);
char * interpret_type(char * type);
CellState * start_CellState_from_GenericMatrix(GenericMatrix * gm);
CellState * end_CellState_from_GenericMatrix(GenericMatrix * gm);
void show_GenericMatrix(GenericMatrix * gm,char padchar,FILE * ofp);
void show_CellState(CellState * cell,char padchar,int num,FILE * ofp);
void show_CellSource(CellSource * cell,char padchar,int num,FILE * ofp);
GenericMatrix * read_GenericMatrix(FILE * ifp);
GenericMatrix * read_GenericMatrix_line(char * line,FILE * ifp);
CollapsableLabel * read_CollapsableLabel_line(char * line);
CellState  * read_CellState_line(char * line,FILE * ifp);
int source_bit2pos(int bit);
CellSource * read_CellSource_line(char * line,FILE * ifp);
char * read_calc_line(char * buffer);


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void add_GenericMatrix_Scope(Scope * sc,GenericMatrix * gm);
ParseError parse_calc_line_GenericMatrix(GenericMatrix * gm,MethodTypeSet * mts,DycWarning * dycw);
boolean calc_window(GenericMatrix * gm);
boolean calc_footprint(GenericMatrix * gm);
boolean cross_reference_state_and_source(GenericMatrix * gm);
void swap_CellSignatureSet(CellSignature ** list,int i,int j) ;
void qsort_CellSignatureSet(CellSignature ** list,int left,int right,int (*comp)(CellSignature * ,CellSignature * ));
void sort_CellSignatureSet(CellSignatureSet * obj,int (*comp)(CellSignature *, CellSignature *));
boolean expand_CellSignatureSet(CellSignatureSet * obj,int len);
void swap_CellState(CellSource ** list,int i,int j) ;
void qsort_CellState(CellSource ** list,int left,int right,int (*comp)(CellSource * ,CellSource * ));
void sort_CellState(CellState * obj,int (*comp)(CellSource *, CellSource *));
boolean expand_CellState(CellState * obj,int len);
void swap_GenericMatrix(CellState ** list,int i,int j) ;
void qsort_GenericMatrix(CellState ** list,int left,int right,int (*comp)(CellState * ,CellState * ));
void sort_GenericMatrix(GenericMatrix * obj,int (*comp)(CellState *, CellState *));
boolean expand_GenericMatrix(GenericMatrix * obj,int len);
void swap_spec_GenericMatrix(CellState ** list,int i,int j) ;
void qsort_spec_GenericMatrix(CellState ** list,int left,int right,int (*comp)(CellState * ,CellState * ));
void sort_spec_GenericMatrix(GenericMatrix * obj,int (*comp)(CellState *, CellState *));
boolean expand_spec_GenericMatrix(GenericMatrix * obj,int len);
void swap_res_GenericMatrix(StructElement  ** list,int i,int j) ;
void qsort_res_GenericMatrix(StructElement  ** list,int left,int right,int (*comp)(StructElement  * ,StructElement  * ));
void sort_res_GenericMatrix(GenericMatrix * obj,int (*comp)(StructElement  *, StructElement  *));
boolean expand_res_GenericMatrix(GenericMatrix * obj,int len);
void swap_cal_GenericMatrix(CollapsableLabel ** list,int i,int j) ;
void qsort_cal_GenericMatrix(CollapsableLabel ** list,int left,int right,int (*comp)(CollapsableLabel * ,CollapsableLabel * ));
void sort_cal_GenericMatrix(GenericMatrix * obj,int (*comp)(CollapsableLabel *, CollapsableLabel *));
boolean expand_cal_GenericMatrix(GenericMatrix * obj,int len);
void swap_ev_GenericMatrix(ExternVariable   ** list,int i,int j) ;
void qsort_ev_GenericMatrix(ExternVariable   ** list,int left,int right,int (*comp)(ExternVariable   * ,ExternVariable   * ));
void sort_ev_GenericMatrix(GenericMatrix * obj,int (*comp)(ExternVariable   *, ExternVariable   *));
boolean expand_ev_GenericMatrix(GenericMatrix * obj,int len);

#ifdef _cplusplus
}
#endif

#endif
