#ifndef DYNAMITEtypeHEADERFILE
#define DYNAMITEtypeHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "wisebase.h"
#include "wisec.h"
#include "exprtree.h"
#include "y.tab.h"
#include "method.h"

#include "input.h"
#include "dpimpl.h"

#define ScopeLISTLENGTH 16

typedef unsigned int ParseError;

enum ParseErrorType {
  PERR_SYNTAX = 1,
  PERR_COMPLEXMETHOD = 2,
  PERR_ARG_NUM_MIS = 4,
  PERR_ARG_UNTYPED = 8,
  PERR_ARG_MISTYPE = 16,
  PERR_OUT_OF_SCOPE = 32,
  PERR_METHOD_SCOPE = 64
};

enum SCOPE_TYPE {
	SCOPE_RESOURCE = 82,
	SCOPE_QUERY,
	SCOPE_TARGET,
	SCOPE_EXTERN };

/*** declare yyparse ****/

void yyparse(void);

struct ScopeUnit {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * name;     
    char * app;  
    char * type;     
    boolean isglobbed;   
    int scope_type;  
    char * no_accept;   /*  for guys which don't like other arguments */ 
    } ;  
/* ScopeUnit defined */ 
#ifndef DYNAMITE_DEFINED_ScopeUnit
typedef struct ScopeUnit ScopeUnit;
#define DYNAMITE_DEFINED_ScopeUnit
#endif


struct Scope {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    ScopeUnit ** su;     
    int len;/* len for above su  */ 
    int maxlen; /* maxlen for above su */ 
    } ;  
/* Scope defined */ 
#ifndef DYNAMITE_DEFINED_Scope
typedef struct Scope Scope;
#define DYNAMITE_DEFINED_Scope
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  copy_MethodTypeSet(mts)
 *
 * Descrip:    copies MethodTypeSet by hard-linking list
 *             members
 *
 *
 *
 * Arg:        mts [UNKN ] Undocumented argument [MethodTypeSet *]
 *
 * Return [UNKN ]  Undocumented return value [MethodTypeSet *]
 *
 */
MethodTypeSet * copy_MethodTypeSet(MethodTypeSet * mts);


/* Function:  std_Dynamite_Scope(void)
 *
 * Descrip:    sets up "standard" dynamite scope.
 *
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Scope *]
 *
 */
Scope * std_Dynamite_Scope(void);


/* Function:  allocd_calc_line(calc_line,sc,mts,dycw,pe,expr)
 *
 * Descrip:    Main function used to access the parser.
 *
 *             This parses the calc_line (without damaging it). Errors will
 *             be placed normally: You should stack errors or catch them if
 *             you want to do something funky. 
 *
 *             At the end of the day, if there was parser syntax error then you
 *             will get NULL, and pe & PERR_SYNTAX will be set. Otherwise you
 *             will get a char * (stringalloc'd). Potentially any number of
 *             PERRs could be set, including unscoped variables or methods etc.
 *
 *
 * Arg:        calc_line [READ ] line to be parsed [char *]
 * Arg:               sc [READ ] scope system to use for this area [Scope *]
 * Arg:              mts [READ ] method and type information (method scope) for this area [MethodTypeSet *]
 * Arg:             dycw [UNKN ] Undocumented argument [DycWarning *]
 * Arg:               pe [WRITE] returned parser errors [ParseError *]
 * Arg:             expr [UNKN ] Undocumented argument [ExprTree **]
 *
 * Return [UNKN ]  stringalloc'd expanded string [char *]
 *
 */
char * allocd_calc_line(char * calc_line,Scope * sc,MethodTypeSet * mts,DycWarning * dycw,ParseError * pe,ExprTree ** expr);


/* Function:  hard_link_ScopeUnit(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ScopeUnit *]
 *
 * Return [UNKN ]  Undocumented return value [ScopeUnit *]
 *
 */
ScopeUnit * hard_link_ScopeUnit(ScopeUnit * obj);


/* Function:  ScopeUnit_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ScopeUnit *]
 *
 */
ScopeUnit * ScopeUnit_alloc(void);


/* Function:  free_ScopeUnit(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ScopeUnit *]
 *
 * Return [UNKN ]  Undocumented return value [ScopeUnit *]
 *
 */
ScopeUnit * free_ScopeUnit(ScopeUnit * obj);


/* Function:  add_Scope(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Scope *]
 * Arg:        add [OWNER] Object to add to the list [ScopeUnit *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean add_Scope(Scope * obj,ScopeUnit * add);


/* Function:  flush_Scope(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [Scope *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_Scope(Scope * obj);


/* Function:  Scope_alloc_std(void)
 *
 * Descrip:    Equivalent to Scope_alloc_len(ScopeLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Scope *]
 *
 */
Scope * Scope_alloc_std(void);


/* Function:  Scope_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [Scope *]
 *
 */
Scope * Scope_alloc_len(int len);


/* Function:  hard_link_Scope(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Scope *]
 *
 * Return [UNKN ]  Undocumented return value [Scope *]
 *
 */
Scope * hard_link_Scope(Scope * obj);


/* Function:  Scope_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Scope *]
 *
 */
Scope * Scope_alloc(void);


/* Function:  free_Scope(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Scope *]
 *
 * Return [UNKN ]  Undocumented return value [Scope *]
 *
 */
Scope * free_Scope(Scope * obj);


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
void complain_ParseError_to_file(ParseError pe,FILE * ofp);
char * type_from_ExprTree(ExprTree * et,Scope * sc,MethodTypeSet * mts);
ParseError typecheck_method(ExprTree * et,Method * me,Scope * sc,MethodTypeSet * mts);
ScopeUnit * ScopeUnit_from_nat(MethodTypeSet * mts,char * name,char * app,char * type);


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
ParseError strcat_ExprTree_Scoped(ExprTree * ExprTree,char * buffer,Scope * sc,MethodTypeSet * mts,DycWarning * dycw,boolean (*finish_parsing)(struct ExprTree *,char *,void *),void * data);
ScopeUnit * ScopeUnit_from_Scope(Scope * sc,char * word);
void swap_Scope(ScopeUnit ** list,int i,int j) ;
void qsort_Scope(ScopeUnit ** list,int left,int right,int (*comp)(ScopeUnit * ,ScopeUnit * ));
void sort_Scope(Scope * obj,int (*comp)(ScopeUnit *, ScopeUnit *));
boolean expand_Scope(Scope * obj,int len);

#ifdef _cplusplus
}
#endif

#endif
