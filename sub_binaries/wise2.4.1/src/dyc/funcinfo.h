#ifndef DYNAMITEfuncinfoHEADERFILE
#define DYNAMITEfuncinfoHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "wisebase.h"
#include "ftext.h"
#include "inputfile.h"

#define FuncInfoLISTLENGTH 32

enum argtype {
  ARGTYPE_UNKNOWN = 43,
  ARGTYPE_READ,
  ARGTYPE_WRITE,
  ARGTYPE_READWRITE,
  ARGTYPE_P2FUNC,
  ARGTYPE_OWNER,
  ARGTYPE_STATIC
};

enum FI_TYPE {
  FI_CALLABLE = 12,
  FI_UNKNOWN,
  FI_INTERNAL
};

#define iscword(c) (c == '_' ? 1 : isalnum(c))


struct ArgInfo {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int  argtype;    
    boolean should_NULL;     
    char * name;     
    char * type;     
    char * desc;     
    int argpos;  
    char * func_decl;   /*  for pointers to functions, so we can get declaration correct */ 
    } ;  
/* ArgInfo defined */ 
#ifndef DYNAMITE_DEFINED_ArgInfo
typedef struct ArgInfo ArgInfo;
#define DYNAMITE_DEFINED_ArgInfo
#endif


struct ErrorInfo {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * value;    
    char * desc;     
    } ;  
/* ErrorInfo defined */ 
#ifndef DYNAMITE_DEFINED_ErrorInfo
typedef struct ErrorInfo ErrorInfo;
#define DYNAMITE_DEFINED_ErrorInfo
#endif


struct FuncInfo {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * name;     
    char * type;     
    char * complete_name;    
    char * stripped_return;  
    Ftext * ft;  
    char * error;    
    ArgInfo   ** arg;    
    int len;/* len for above arg  */ 
    int maxlen; /* maxlen for above arg */ 
    ErrorInfo ** err;    
    int err_len;/* len for above err  */ 
    int err_maxlen; /* maxlen for above err */ 
    char * sdesc;    
    ArgInfo   * ret;     
    int functype;    
    int line_in_c;   
    int infopos;     
    char * simple;   
    boolean is_hand_written;     
    } ;  
/* FuncInfo defined */ 
#ifndef DYNAMITE_DEFINED_FuncInfo
typedef struct FuncInfo FuncInfo;
#define DYNAMITE_DEFINED_FuncInfo
#endif


struct ModuleInfo {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Ftext * ft;  
    } ;  
/* ModuleInfo defined */ 
#ifndef DYNAMITE_DEFINED_ModuleInfo
typedef struct ModuleInfo ModuleInfo;
#define DYNAMITE_DEFINED_ModuleInfo
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  show_eddystyle_FuncInfo(fi,ofp)
 *
 * Descrip:    shows functions in
 *
 *              *Function:
 *              *
 *              *des
 *              *
 *              * Arg
 *              *
 *
 *             Returns number of lines printed
 *
 *
 * Arg:         fi [UNKN ] Undocumented argument [FuncInfo *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int show_eddystyle_FuncInfo(FuncInfo * fi,FILE * ofp) ;


/* Function:  hard_link_ArgInfo(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ArgInfo *]
 *
 * Return [UNKN ]  Undocumented return value [ArgInfo *]
 *
 */
ArgInfo * hard_link_ArgInfo(ArgInfo * obj);


/* Function:  ArgInfo_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ArgInfo *]
 *
 */
ArgInfo * ArgInfo_alloc(void);


/* Function:  free_ArgInfo(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ArgInfo *]
 *
 * Return [UNKN ]  Undocumented return value [ArgInfo *]
 *
 */
ArgInfo * free_ArgInfo(ArgInfo * obj);


/* Function:  hard_link_ErrorInfo(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ErrorInfo *]
 *
 * Return [UNKN ]  Undocumented return value [ErrorInfo *]
 *
 */
ErrorInfo * hard_link_ErrorInfo(ErrorInfo * obj);


/* Function:  ErrorInfo_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ErrorInfo *]
 *
 */
ErrorInfo * ErrorInfo_alloc(void);


/* Function:  free_ErrorInfo(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ErrorInfo *]
 *
 * Return [UNKN ]  Undocumented return value [ErrorInfo *]
 *
 */
ErrorInfo * free_ErrorInfo(ErrorInfo * obj);


/* Function:  add_FuncInfo(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [FuncInfo *]
 * Arg:        add [OWNER] Object to add to the list [ArgInfo   *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean add_FuncInfo(FuncInfo * obj,ArgInfo   * add);


/* Function:  flush_FuncInfo(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [FuncInfo *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_FuncInfo(FuncInfo * obj);


/* Function:  add_err_FuncInfo(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [FuncInfo *]
 * Arg:        add [OWNER] Object to add to the list [ErrorInfo *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean add_err_FuncInfo(FuncInfo * obj,ErrorInfo * add);


/* Function:  flush_err_FuncInfo(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [FuncInfo *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_err_FuncInfo(FuncInfo * obj);


/* Function:  FuncInfo_alloc_std(void)
 *
 * Descrip:    Equivalent to FuncInfo_alloc_len(FuncInfoLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [FuncInfo *]
 *
 */
FuncInfo * FuncInfo_alloc_std(void);


/* Function:  FuncInfo_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [FuncInfo *]
 *
 */
FuncInfo * FuncInfo_alloc_len(int len);


/* Function:  hard_link_FuncInfo(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [FuncInfo *]
 *
 * Return [UNKN ]  Undocumented return value [FuncInfo *]
 *
 */
FuncInfo * hard_link_FuncInfo(FuncInfo * obj);


/* Function:  FuncInfo_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [FuncInfo *]
 *
 */
FuncInfo * FuncInfo_alloc(void);


/* Function:  free_FuncInfo(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [FuncInfo *]
 *
 * Return [UNKN ]  Undocumented return value [FuncInfo *]
 *
 */
FuncInfo * free_FuncInfo(FuncInfo * obj);


/* Function:  hard_link_ModuleInfo(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ModuleInfo *]
 *
 * Return [UNKN ]  Undocumented return value [ModuleInfo *]
 *
 */
ModuleInfo * hard_link_ModuleInfo(ModuleInfo * obj);


/* Function:  ModuleInfo_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ModuleInfo *]
 *
 */
ModuleInfo * ModuleInfo_alloc(void);


/* Function:  free_ModuleInfo(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ModuleInfo *]
 *
 * Return [UNKN ]  Undocumented return value [ModuleInfo *]
 *
 */
ModuleInfo * free_ModuleInfo(ModuleInfo * obj);


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
void dump_FuncInfo(FuncInfo * fi,FILE * ofp ) ;
char * ArgType_to_string(int type);
int max_argame(FuncInfo * fi );
int show_eddystyle_ArgInfo(ArgInfo * ai,int depth,int namedepth,FILE * ofp);
void sort_FuncInfo_by_position(FuncInfo * fi);
int compare_ArgInfo_pos(ArgInfo * one,ArgInfo * two);
boolean reconcile_FuncInfo_with_funcstr(FuncInfo * fi,char * pass_str);
boolean reconcile_FuncInfo_with_pfunc(FuncInfo * fi,char * str,int pos);
FuncInfo * unknown_user_FuncInfo(char * funcstr);
boolean reconcile_FuncInfo_with_argstr(FuncInfo * fi,char * str,int pos) ;
ArgInfo  * get_ArgInfo_by_name(FuncInfo * fi,char * str);
ArgInfo * ArgInfo_in_FuncInfo_from_varstr(FuncInfo * fi,char * str,...);
FuncInfo * FuncInfo_named_from_varstr(int type,char * str, ...);
FuncInfo * FuncInfo_from_str(char * str);
ModuleInfo * read_ModuleInfo_line(char * line,FILE * ifp);
FuncInfo * read_FuncInfo_line(char * line,FILE * ifp);
int get_arg_type(char * line,boolean * should_NULL);
ArgInfo * read_ArgInfo_line(char * line);


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void swap_FuncInfo(ArgInfo   ** list,int i,int j) ;
void qsort_FuncInfo(ArgInfo   ** list,int left,int right,int (*comp)(ArgInfo   * ,ArgInfo   * ));
void sort_FuncInfo(FuncInfo * obj,int (*comp)(ArgInfo   *, ArgInfo   *));
boolean expand_FuncInfo(FuncInfo * obj,int len);
void swap_err_FuncInfo(ErrorInfo ** list,int i,int j) ;
void qsort_err_FuncInfo(ErrorInfo ** list,int left,int right,int (*comp)(ErrorInfo * ,ErrorInfo * ));
void sort_err_FuncInfo(FuncInfo * obj,int (*comp)(ErrorInfo *, ErrorInfo *));
boolean expand_err_FuncInfo(FuncInfo * obj,int len);

#ifdef _cplusplus
}
#endif

#endif
