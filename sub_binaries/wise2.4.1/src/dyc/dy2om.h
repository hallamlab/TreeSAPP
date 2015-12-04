#ifndef DYNAMITEdy2omHEADERFILE
#define DYNAMITEdy2omHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna2.h"
#include "exprtree.h"

enum OmStateType {
	OmState_TYPE_UNKNOWN = 67,
	OmState_TYPE_STATE,
	OmState_TYPE_SEMISTATE,
	OmState_TYPE_MIDSTATE };


#define OneModelLISTLENGTH 64


typedef struct descent_return {
  int is_doable; /* 1 is yes, 0 cannot be handled by Compugen */
  int is_index;  /* 1 is yes, 0 is it is a float/number */
  int is_query;  /* when is index, 1 means in query, 0 in target */
  int is_variable;
  struct OmUnit * unit; /* if non-NULL, this is the "done" unit */
} CugenYaccReturn;



struct OmState {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * name;     
    int  type;   
    } ;  
/* OmState defined */ 
#ifndef DYNAMITE_DEFINED_OmState
typedef struct OmState OmState;
#define DYNAMITE_DEFINED_OmState
#endif


struct OmUnit {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * c_string;    /*  string of the actual C to use in making profile line */ 
    int base;   /*  offset in profile line, must be the same as t_prf or w_base */ 
    int length; /*  1 for t_prf lines. Longer for w_base lines. */ 
    int is_tprf;     
    } ;  
/* OmUnit defined */ 
#ifndef DYNAMITE_DEFINED_OmUnit
typedef struct OmUnit OmUnit;
#define DYNAMITE_DEFINED_OmUnit
#endif


struct OmTrans {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    OmState * from;  
    OmState * to;    
    int dx;  
    int dy;  
    int t_prf;   
    int dprf;    
    int t_seq;   
    int dseq;    
    int w_base;  
    int w_seq;   
    int dwx;     
    int dwy;     
    OmUnit * tprf_unit;  
    OmUnit * wbase_unit;     
    int current_base;    
    } ;  
/* OmTrans defined */ 
#ifndef DYNAMITE_DEFINED_OmTrans
typedef struct OmTrans OmTrans;
#define DYNAMITE_DEFINED_OmTrans
#endif


struct OmTransFunc {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * name;     
    Method * me;     
    } ;  
/* OmTransFunc defined */ 
#ifndef DYNAMITE_DEFINED_OmTransFunc
typedef struct OmTransFunc OmTransFunc;
#define DYNAMITE_DEFINED_OmTransFunc
#endif


struct OneModel {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * name;     
    int state_num;   
    int semistate_num;   
    int midstate_num;    
    int profile_line_size;   
    OmTrans ** trans;    
    int len;/* len for above trans  */ 
    int maxlen; /* maxlen for above trans */ 
    OmState ** state;    
    int st_len; /* len for above state  */ 
    int st_maxlen;  /* maxlen for above state */ 
    OmTransFunc ** tfunc;    
    int tf_len; /* len for above tfunc  */ 
    int tf_maxlen;  /* maxlen for above tfunc */ 
    } ;  
/* OneModel defined */ 
#ifndef DYNAMITE_DEFINED_OneModel
typedef struct OneModel OneModel;
#define DYNAMITE_DEFINED_OneModel
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  descend_ExprTree_Cugen(et,mts,is_attachable,sc,ot)
 *
 * Descrip:    Recursive parser of the yacc grammar for
 *             compugen.
 *
 *             If it is impossible for compugen to do, then
 *             is_doable == 0 and the descent should abort asap.
 *
 *             If not, at the lowest level of being able to map
 *             to a compugen number, it is mapped. 
 *
 *             targetname should be entire scope.
 *
 *
 * Arg:                   et [UNKN ] Undocumented argument [ExprTree *]
 * Arg:                  mts [UNKN ] Undocumented argument [MethodTypeSet *]
 * Arg:        is_attachable [UNKN ] Undocumented argument [boolean]
 * Arg:                   sc [UNKN ] Undocumented argument [Scope *]
 * Arg:                   ot [UNKN ] Undocumented argument [OmTrans *]
 *
 * Return [UNKN ]  Undocumented return value [CugenYaccReturn]
 *
 */
CugenYaccReturn descend_ExprTree_Cugen(ExprTree * et,MethodTypeSet * mts,boolean is_attachable,Scope * sc,OmTrans * ot);


/* Function:  write_OneModel(om,ofp)
 *
 * Descrip:    writes the model definition
 *
 *
 * Arg:         om [UNKN ] Undocumented argument [OneModel *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void write_OneModel(OneModel * om,FILE * ofp);


/* Function:  write_OmTrans(number,ot,ofp)
 *
 * Descrip:    writes one transition line
 *
 *
 * Arg:        number [UNKN ] Undocumented argument [int]
 * Arg:            ot [UNKN ] Undocumented argument [OmTrans *]
 * Arg:           ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void write_OmTrans(int number,OmTrans * ot,FILE * ofp);


/* Function:  hard_link_OmState(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [OmState *]
 *
 * Return [UNKN ]  Undocumented return value [OmState *]
 *
 */
OmState * hard_link_OmState(OmState * obj);


/* Function:  OmState_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [OmState *]
 *
 */
OmState * OmState_alloc(void);


/* Function:  free_OmState(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [OmState *]
 *
 * Return [UNKN ]  Undocumented return value [OmState *]
 *
 */
OmState * free_OmState(OmState * obj);


/* Function:  hard_link_OmUnit(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [OmUnit *]
 *
 * Return [UNKN ]  Undocumented return value [OmUnit *]
 *
 */
OmUnit * hard_link_OmUnit(OmUnit * obj);


/* Function:  OmUnit_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [OmUnit *]
 *
 */
OmUnit * OmUnit_alloc(void);


/* Function:  free_OmUnit(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [OmUnit *]
 *
 * Return [UNKN ]  Undocumented return value [OmUnit *]
 *
 */
OmUnit * free_OmUnit(OmUnit * obj);


/* Function:  hard_link_OmTrans(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [OmTrans *]
 *
 * Return [UNKN ]  Undocumented return value [OmTrans *]
 *
 */
OmTrans * hard_link_OmTrans(OmTrans * obj);


/* Function:  OmTrans_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [OmTrans *]
 *
 */
OmTrans * OmTrans_alloc(void);


/* Function:  free_OmTrans(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [OmTrans *]
 *
 * Return [UNKN ]  Undocumented return value [OmTrans *]
 *
 */
OmTrans * free_OmTrans(OmTrans * obj);


/* Function:  hard_link_OmTransFunc(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [OmTransFunc *]
 *
 * Return [UNKN ]  Undocumented return value [OmTransFunc *]
 *
 */
OmTransFunc * hard_link_OmTransFunc(OmTransFunc * obj);


/* Function:  OmTransFunc_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [OmTransFunc *]
 *
 */
OmTransFunc * OmTransFunc_alloc(void);


/* Function:  free_OmTransFunc(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [OmTransFunc *]
 *
 * Return [UNKN ]  Undocumented return value [OmTransFunc *]
 *
 */
OmTransFunc * free_OmTransFunc(OmTransFunc * obj);


/* Function:  add_OneModel(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [OneModel *]
 * Arg:        add [OWNER] Object to add to the list [OmTrans *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean add_OneModel(OneModel * obj,OmTrans * add);


/* Function:  flush_OneModel(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [OneModel *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_OneModel(OneModel * obj);


/* Function:  add_st_OneModel(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [OneModel *]
 * Arg:        add [OWNER] Object to add to the list [OmState *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean add_st_OneModel(OneModel * obj,OmState * add);


/* Function:  flush_st_OneModel(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [OneModel *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_st_OneModel(OneModel * obj);


/* Function:  add_tf_OneModel(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [OneModel *]
 * Arg:        add [OWNER] Object to add to the list [OmTransFunc *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean add_tf_OneModel(OneModel * obj,OmTransFunc * add);


/* Function:  flush_tf_OneModel(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [OneModel *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_tf_OneModel(OneModel * obj);


/* Function:  OneModel_alloc_std(void)
 *
 * Descrip:    Equivalent to OneModel_alloc_len(OneModelLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [OneModel *]
 *
 */
OneModel * OneModel_alloc_std(void);


/* Function:  OneModel_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [OneModel *]
 *
 */
OneModel * OneModel_alloc_len(int len);


/* Function:  hard_link_OneModel(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [OneModel *]
 *
 * Return [UNKN ]  Undocumented return value [OneModel *]
 *
 */
OneModel * hard_link_OneModel(OneModel * obj);


/* Function:  OneModel_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [OneModel *]
 *
 */
OneModel * OneModel_alloc(void);


/* Function:  free_OneModel(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [OneModel *]
 *
 * Return [UNKN ]  Undocumented return value [OneModel *]
 *
 */
OneModel * free_OneModel(OneModel * obj);


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
OneModel * OneModel_from_GenericMatrix(GenericMatrix * gm,MethodTypeSet * mts);
OmTrans * OmTrans_from_CellSource(int base_point,CellSource * s,MethodTypeSet * mts,Scope * sc);
OmState * OmState_from_CellState(CellState * cs);
void set_CugenYaccReturn(CugenYaccReturn * ret);


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void swap_OneModel(OmTrans ** list,int i,int j) ;
void qsort_OneModel(OmTrans ** list,int left,int right,int (*comp)(OmTrans * ,OmTrans * ));
void sort_OneModel(OneModel * obj,int (*comp)(OmTrans *, OmTrans *));
boolean expand_OneModel(OneModel * obj,int len);
void swap_st_OneModel(OmState ** list,int i,int j) ;
void qsort_st_OneModel(OmState ** list,int left,int right,int (*comp)(OmState * ,OmState * ));
void sort_st_OneModel(OneModel * obj,int (*comp)(OmState *, OmState *));
boolean expand_st_OneModel(OneModel * obj,int len);
void swap_tf_OneModel(OmTransFunc ** list,int i,int j) ;
void qsort_tf_OneModel(OmTransFunc ** list,int left,int right,int (*comp)(OmTransFunc * ,OmTransFunc * ));
void sort_tf_OneModel(OneModel * obj,int (*comp)(OmTransFunc *, OmTransFunc *));
boolean expand_tf_OneModel(OneModel * obj,int len);

#ifdef _cplusplus
}
#endif

#endif
