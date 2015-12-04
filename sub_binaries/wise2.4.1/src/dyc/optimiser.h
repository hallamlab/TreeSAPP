#ifndef DYNAMITEoptimiserHEADERFILE
#define DYNAMITEoptimiserHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna2.h"
#include "type.h"

#define CommonSubExpressionSetLISTLENGTH 512

#define CSE_I_DEP 1
#define CSE_J_DEP 2

#define IS_I_DEP_CSE(cse) ((cse->type & CSE_I_DEP) == CSE_I_DEP ? TRUE : FALSE)
#define IS_J_DEP_CSE(cse) ((cse->type & CSE_J_DEP) == CSE_J_DEP ? TRUE : FALSE)
#define IS_NON_IJ_DEP_CSE(cse) (cse->type == 0 ? TRUE : FALSE)
struct CommonSubExpression {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    ExprTree * expr;     
    int id; /*  able to set things as being 1,2,3 etc */ 
    int type;   /*  is i or j or i,j or no dependent */ 
    int number; /*  number of times this sub expression is used */ 
    } ;  
/* CommonSubExpression defined */ 
#ifndef DYNAMITE_DEFINED_CommonSubExpression
typedef struct CommonSubExpression CommonSubExpression;
#define DYNAMITE_DEFINED_CommonSubExpression
#endif


struct CommonSubExpressionSet {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    CommonSubExpression ** cse; /*  list of common sub expressions */ 
    int len;/* len for above cse  */ 
    int maxlen; /* maxlen for above cse */ 
    } ;  
/* CommonSubExpressionSet defined */ 
#ifndef DYNAMITE_DEFINED_CommonSubExpressionSet
typedef struct CommonSubExpressionSet CommonSubExpressionSet;
#define DYNAMITE_DEFINED_CommonSubExpressionSet
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  expr_ij_dependence(expr,ijdep)
 *
 * Descrip:    figures out whether this cse is i or j dependent
 *
 *
 * Arg:         expr [UNKN ] Undocumented argument [ExprTree *]
 * Arg:        ijdep [UNKN ] Undocumented argument [int *]
 *
 */
void expr_ij_dependence(ExprTree * expr,int * ijdep);


/* Function:  strcat_cses_ExprTree(epr,buffer,sc,mts,dpi)
 *
 * Descrip:    Writes code with sub expressions in the correct places
 *
 *
 * Arg:           epr [UNKN ] Undocumented argument [ExprTree *]
 * Arg:        buffer [UNKN ] Undocumented argument [char *]
 * Arg:            sc [UNKN ] Undocumented argument [Scope *]
 * Arg:           mts [UNKN ] Undocumented argument [MethodTypeSet *]
 * Arg:           dpi [UNKN ] Undocumented argument [DPImplementation *]
 *
 */
void strcat_cses_ExprTree(ExprTree * epr,char * buffer,Scope * sc,MethodTypeSet * mts,DPImplementation * dpi);


/* Function:  show_CommonSubExpressionSet(cses,ofp)
 *
 * Descrip:    Shows common sub expression set
 *
 *
 * Arg:        cses [UNKN ] Undocumented argument [CommonSubExpressionSet *]
 * Arg:         ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void show_CommonSubExpressionSet(CommonSubExpressionSet * cses,FILE * ofp);


/* Function:  find_CommonSubExpressions(gm,source_ind_promote)
 *
 * Descrip:    Makes a CommonSubExpressionSet from the
 *             GenericMatrix structure. 
 *
 *             if source_ind_promote is true, then source_independent type systems
 *             are automatically promoted as a common sub expression, regardless of 
 *             number of subexpressions found
 *
 *
 * Arg:                        gm [UNKN ] Undocumented argument [GenericMatrix *]
 * Arg:        source_ind_promote [UNKN ] Undocumented argument [boolean]
 *
 * Return [UNKN ]  Undocumented return value [CommonSubExpressionSet *]
 *
 */
CommonSubExpressionSet * find_CommonSubExpressions(GenericMatrix * gm,boolean source_ind_promote);


/* Function:  attach_expressions_to_list(list,no,start)
 *
 * Descrip:    Attaches only ETR_EXPRESSIONS into the list.
 *             Really a subroutine for find_CommonSubExpressions
 *
 *
 * Arg:         list [UNKN ] Undocumented argument [ExprTree **]
 * Arg:           no [UNKN ] Undocumented argument [int *]
 * Arg:        start [UNKN ] Undocumented argument [ExprTree *]
 *
 */
void attach_expressions_to_list(ExprTree ** list,int * no,ExprTree * start);


/* Function:  hard_link_CommonSubExpression(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [CommonSubExpression *]
 *
 * Return [UNKN ]  Undocumented return value [CommonSubExpression *]
 *
 */
CommonSubExpression * hard_link_CommonSubExpression(CommonSubExpression * obj);


/* Function:  CommonSubExpression_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [CommonSubExpression *]
 *
 */
CommonSubExpression * CommonSubExpression_alloc(void);


/* Function:  free_CommonSubExpression(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [CommonSubExpression *]
 *
 * Return [UNKN ]  Undocumented return value [CommonSubExpression *]
 *
 */
CommonSubExpression * free_CommonSubExpression(CommonSubExpression * obj);


/* Function:  add_CommonSubExpressionSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [CommonSubExpressionSet *]
 * Arg:        add [OWNER] Object to add to the list [CommonSubExpression *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean add_CommonSubExpressionSet(CommonSubExpressionSet * obj,CommonSubExpression * add);


/* Function:  flush_CommonSubExpressionSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [CommonSubExpressionSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_CommonSubExpressionSet(CommonSubExpressionSet * obj);


/* Function:  CommonSubExpressionSet_alloc_std(void)
 *
 * Descrip:    Equivalent to CommonSubExpressionSet_alloc_len(CommonSubExpressionSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [CommonSubExpressionSet *]
 *
 */
CommonSubExpressionSet * CommonSubExpressionSet_alloc_std(void);


/* Function:  CommonSubExpressionSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [CommonSubExpressionSet *]
 *
 */
CommonSubExpressionSet * CommonSubExpressionSet_alloc_len(int len);


/* Function:  hard_link_CommonSubExpressionSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [CommonSubExpressionSet *]
 *
 * Return [UNKN ]  Undocumented return value [CommonSubExpressionSet *]
 *
 */
CommonSubExpressionSet * hard_link_CommonSubExpressionSet(CommonSubExpressionSet * obj);


/* Function:  CommonSubExpressionSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [CommonSubExpressionSet *]
 *
 */
CommonSubExpressionSet * CommonSubExpressionSet_alloc(void);


/* Function:  free_CommonSubExpressionSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [CommonSubExpressionSet *]
 *
 * Return [UNKN ]  Undocumented return value [CommonSubExpressionSet *]
 *
 */
CommonSubExpressionSet * free_CommonSubExpressionSet(CommonSubExpressionSet * obj);


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
boolean cses_expr_placer(ExprTree * etr,char * buffer,void * data);
void id_CommonSubExpressionSet(CommonSubExpressionSet * cses);
boolean should_store_ExprTree_cse(ExprTree * start);
boolean identical_ExprTree(ExprTree * one,ExprTree * two);
void swap_CommonSubExpressionSet(CommonSubExpression ** list,int i,int j) ;
void qsort_CommonSubExpressionSet(CommonSubExpression ** list,int left,int right,int (*comp)(CommonSubExpression * ,CommonSubExpression * ));
void sort_CommonSubExpressionSet(CommonSubExpressionSet * obj,int (*comp)(CommonSubExpression *, CommonSubExpression *));
boolean expand_CommonSubExpressionSet(CommonSubExpressionSet * obj,int len);

#ifdef _cplusplus
}
#endif

#endif
