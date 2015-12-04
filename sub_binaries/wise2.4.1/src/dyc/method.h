#ifndef DYNAMITEmethodHEADERFILE
#define DYNAMITEmethodHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "wisebase.h"
#include "input.h"
#include "wisec.h"

#define MethodLISTLENGTH 16
#define MethodTypeSetLISTLENGTH 16


enum COMPUGEN_METHOD_MAP {
	CUGEN_METHOD_UNKNOWN = 78,
	CUGEN_METHOD_TSEQ,
	CUGEN_METHOD_WINDEX };

struct MethodArg {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * type;     
    } ;  
/* MethodArg defined */ 
#ifndef DYNAMITE_DEFINED_MethodArg
typedef struct MethodArg MethodArg;
#define DYNAMITE_DEFINED_MethodArg
#endif


struct Method {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * logical;  
    char * real;     
    char * retstr;   
    MethodArg ** ma;     
    int len;/* len for above ma  */ 
    int maxlen; /* maxlen for above ma */ 
    char * cugen_map;    
    int cugen_type;  
    } ;  
/* Method defined */ 
#ifndef DYNAMITE_DEFINED_Method
typedef struct Method Method;
#define DYNAMITE_DEFINED_Method
#endif


struct Type {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    boolean is_database;     
    boolean is_thread_safe;  
    char * logical;  
    char * real;     
    char * database_type;    
    char * get_id_func;  
    char * init_func;    
    char * reload_func;  
    char * close_func;   
    char * dataentry_add;    
    char * maxlen;   
    char * hard_link_func;   
    char * free_func;    
    Input * in;  
    } ;  
/* Type defined */ 
#ifndef DYNAMITE_DEFINED_Type
typedef struct Type Type;
#define DYNAMITE_DEFINED_Type
#endif


struct MethodTypeSet {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Method ** me;    
    int me_len; /* len for above me  */ 
    int me_maxlen;  /* maxlen for above me */ 
    Type   ** ty;    
    int ty_len; /* len for above ty  */ 
    int ty_maxlen;  /* maxlen for above ty */ 
    Input  ** in;    
    int in_len; /* len for above in  */ 
    int in_maxlen;  /* maxlen for above in */ 
    } ;  
/* MethodTypeSet defined */ 
#ifndef DYNAMITE_DEFINED_MethodTypeSet
typedef struct MethodTypeSet MethodTypeSet;
#define DYNAMITE_DEFINED_MethodTypeSet
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  StructElement_from_MethodTypeSet(name,type,mts)
 *
 * Descrip:    function which handles the logical->real mapping
 *
 *             At the moment "unmappable" types get assummed to be C types,
 *             trigger a warning and return the correct thing.
 *
 *
 * Arg:        name [UNKN ] Undocumented argument [char *]
 * Arg:        type [UNKN ] Undocumented argument [char *]
 * Arg:         mts [UNKN ] Undocumented argument [MethodTypeSet *]
 *
 * Return [UNKN ]  Undocumented return value [StructElement *]
 *
 */
StructElement * StructElement_from_MethodTypeSet(char * name,char * type,MethodTypeSet * mts);


/* Function:  read_Type_line(line,ifp)
 *
 * Descrip:    reads in a type structure from a line starting
 *
 *             type
 *             etc
 *
 *
 * Arg:        line [UNKN ] first line with type [char *]
 * Arg:         ifp [UNKN ] read file [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [Type *]
 *
 */
Type * read_Type_line(char * line,FILE * ifp);


/* Function:  MethodArg_from_line(line)
 *
 * Descrip:    reads in from line like
 *
 *             arg PROTEIN
 *
 *
 * Arg:        line [UNKN ] pointer to buffer [char *]
 *
 * Return [UNKN ]  Undocumented return value [MethodArg *]
 *
 */
MethodArg * MethodArg_from_line(char * line);


/* Function:  hard_link_MethodArg(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [MethodArg *]
 *
 * Return [UNKN ]  Undocumented return value [MethodArg *]
 *
 */
MethodArg * hard_link_MethodArg(MethodArg * obj);


/* Function:  MethodArg_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [MethodArg *]
 *
 */
MethodArg * MethodArg_alloc(void);


/* Function:  free_MethodArg(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [MethodArg *]
 *
 * Return [UNKN ]  Undocumented return value [MethodArg *]
 *
 */
MethodArg * free_MethodArg(MethodArg * obj);


/* Function:  add_Method(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Method *]
 * Arg:        add [OWNER] Object to add to the list [MethodArg *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean add_Method(Method * obj,MethodArg * add);


/* Function:  flush_Method(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [Method *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_Method(Method * obj);


/* Function:  Method_alloc_std(void)
 *
 * Descrip:    Equivalent to Method_alloc_len(MethodLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Method *]
 *
 */
Method * Method_alloc_std(void);


/* Function:  Method_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [Method *]
 *
 */
Method * Method_alloc_len(int len);


/* Function:  hard_link_Method(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Method *]
 *
 * Return [UNKN ]  Undocumented return value [Method *]
 *
 */
Method * hard_link_Method(Method * obj);


/* Function:  Method_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Method *]
 *
 */
Method * Method_alloc(void);


/* Function:  free_Method(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Method *]
 *
 * Return [UNKN ]  Undocumented return value [Method *]
 *
 */
Method * free_Method(Method * obj);


/* Function:  hard_link_Type(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Type *]
 *
 * Return [UNKN ]  Undocumented return value [Type *]
 *
 */
Type * hard_link_Type(Type * obj);


/* Function:  Type_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Type *]
 *
 */
Type * Type_alloc(void);


/* Function:  free_Type(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Type *]
 *
 * Return [UNKN ]  Undocumented return value [Type *]
 *
 */
Type * free_Type(Type * obj);


/* Function:  add_me_MethodTypeSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [MethodTypeSet *]
 * Arg:        add [OWNER] Object to add to the list [Method *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean add_me_MethodTypeSet(MethodTypeSet * obj,Method * add);


/* Function:  flush_me_MethodTypeSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [MethodTypeSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_me_MethodTypeSet(MethodTypeSet * obj);


/* Function:  add_ty_MethodTypeSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [MethodTypeSet *]
 * Arg:        add [OWNER] Object to add to the list [Type   *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean add_ty_MethodTypeSet(MethodTypeSet * obj,Type   * add);


/* Function:  flush_ty_MethodTypeSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [MethodTypeSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_ty_MethodTypeSet(MethodTypeSet * obj);


/* Function:  add_in_MethodTypeSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [MethodTypeSet *]
 * Arg:        add [OWNER] Object to add to the list [Input  *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean add_in_MethodTypeSet(MethodTypeSet * obj,Input  * add);


/* Function:  flush_in_MethodTypeSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [MethodTypeSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_in_MethodTypeSet(MethodTypeSet * obj);


/* Function:  MethodTypeSet_alloc_std(void)
 *
 * Descrip:    Equivalent to MethodTypeSet_alloc_len(MethodTypeSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [MethodTypeSet *]
 *
 */
MethodTypeSet * MethodTypeSet_alloc_std(void);


/* Function:  MethodTypeSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [MethodTypeSet *]
 *
 */
MethodTypeSet * MethodTypeSet_alloc_len(int len);


/* Function:  hard_link_MethodTypeSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [MethodTypeSet *]
 *
 * Return [UNKN ]  Undocumented return value [MethodTypeSet *]
 *
 */
MethodTypeSet * hard_link_MethodTypeSet(MethodTypeSet * obj);


/* Function:  MethodTypeSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [MethodTypeSet *]
 *
 */
MethodTypeSet * MethodTypeSet_alloc(void);


/* Function:  free_MethodTypeSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [MethodTypeSet *]
 *
 * Return [UNKN ]  Undocumented return value [MethodTypeSet *]
 *
 */
MethodTypeSet * free_MethodTypeSet(MethodTypeSet * obj);


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
void show_MethodTypeSet(MethodTypeSet * mts,FILE * ofp);
void show_Method(Method * m,FILE * ofp);
void show_MethodArg(MethodArg * ma,FILE * ofp);
void show_Type(Type * ty,FILE * ofp);
MethodTypeSet * read_MethodTypeSet_filename(char * filename);
boolean read_into_MethodTypeSet_filename(MethodTypeSet * mts,char * filename);
MethodTypeSet * standard_dynamite_MethodTypeSet(void);
MethodTypeSet * empty_MethodTypeSet(void);
MethodTypeSet * read_MethodTypeSet(FILE * ifp);
boolean read_into_MethodTypeSet(MethodTypeSet * mts,FILE * ifp);
boolean is_database_type(Type * ty);
Method * read_Method_line(char * line,FILE * ifp);


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
Type * Type_from_name(MethodTypeSet * mts,char * name);
boolean compare_type(char * s,char * t);
Method * Method_from_name(MethodTypeSet * mts,char * name);
StructElement * StructElement_from_nameandtype(char * name,char * type);
void swap_Method(MethodArg ** list,int i,int j) ;
void qsort_Method(MethodArg ** list,int left,int right,int (*comp)(MethodArg * ,MethodArg * ));
void sort_Method(Method * obj,int (*comp)(MethodArg *, MethodArg *));
boolean expand_Method(Method * obj,int len);
void swap_me_MethodTypeSet(Method ** list,int i,int j) ;
void qsort_me_MethodTypeSet(Method ** list,int left,int right,int (*comp)(Method * ,Method * ));
void sort_me_MethodTypeSet(MethodTypeSet * obj,int (*comp)(Method *, Method *));
boolean expand_me_MethodTypeSet(MethodTypeSet * obj,int len);
void swap_ty_MethodTypeSet(Type   ** list,int i,int j) ;
void qsort_ty_MethodTypeSet(Type   ** list,int left,int right,int (*comp)(Type   * ,Type   * ));
void sort_ty_MethodTypeSet(MethodTypeSet * obj,int (*comp)(Type   *, Type   *));
boolean expand_ty_MethodTypeSet(MethodTypeSet * obj,int len);
void swap_in_MethodTypeSet(Input  ** list,int i,int j) ;
void qsort_in_MethodTypeSet(Input  ** list,int left,int right,int (*comp)(Input  * ,Input  * ));
void sort_in_MethodTypeSet(MethodTypeSet * obj,int (*comp)(Input  *, Input  *));
boolean expand_in_MethodTypeSet(MethodTypeSet * obj,int len);

#ifdef _cplusplus
}
#endif

#endif
