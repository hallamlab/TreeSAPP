#ifndef DYNAMITEapiHEADERFILE
#define DYNAMITEapiHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dynfile.h"
#include "objectinfo.h"
#include "wisec.h"

#define APIObjectLISTLENGTH 64
#define dynAPILISTLENGTH 64


struct APIfunc {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * name;     
    FuncInfo * fi;  /*   hard linked */ 
    boolean only_C;  
    } ;  
/* APIfunc defined */ 
#ifndef DYNAMITE_DEFINED_APIfunc
typedef struct APIfunc APIfunc;
#define DYNAMITE_DEFINED_APIfunc
#endif


struct APIObject {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * name;     
    APIfunc * destructor;    
    APIfunc * hard_linker;   
    APIfunc ** member;   
    int len;/* len for above member  */ 
    int maxlen; /* maxlen for above member */ 
    ObjectInfo * info;   
    StructHolder * sh;   
    } ;  
/* APIObject defined */ 
#ifndef DYNAMITE_DEFINED_APIObject
typedef struct APIObject APIObject;
#define DYNAMITE_DEFINED_APIObject
#endif


/* Object dynAPI
 *
 * Descrip: This Object contains the api definition
 *        for the module. The api is used to 
 *        generate a name-space clear C api and
 *        a Perl interface at the moment. 
 *
 *
 */
struct dynAPI {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    APIfunc ** non_obj;  
    int non_len;/* len for above non_obj  */ 
    int non_maxlen; /* maxlen for above non_obj */ 
    APIObject ** obj;    
    int len;/* len for above obj  */ 
    int maxlen; /* maxlen for above obj */ 
    } ;  
/* dynAPI defined */ 
#ifndef DYNAMITE_DEFINED_dynAPI
typedef struct dynAPI dynAPI;
#define DYNAMITE_DEFINED_dynAPI
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_APIfunc(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [APIfunc *]
 *
 * Return [UNKN ]  Undocumented return value [APIfunc *]
 *
 */
APIfunc * hard_link_APIfunc(APIfunc * obj);


/* Function:  APIfunc_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [APIfunc *]
 *
 */
APIfunc * APIfunc_alloc(void);


/* Function:  free_APIfunc(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [APIfunc *]
 *
 * Return [UNKN ]  Undocumented return value [APIfunc *]
 *
 */
APIfunc * free_APIfunc(APIfunc * obj);


/* Function:  add_APIObject(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [APIObject *]
 * Arg:        add [OWNER] Object to add to the list [APIfunc *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean add_APIObject(APIObject * obj,APIfunc * add);


/* Function:  flush_APIObject(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [APIObject *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_APIObject(APIObject * obj);


/* Function:  APIObject_alloc_std(void)
 *
 * Descrip:    Equivalent to APIObject_alloc_len(APIObjectLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [APIObject *]
 *
 */
APIObject * APIObject_alloc_std(void);


/* Function:  APIObject_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [APIObject *]
 *
 */
APIObject * APIObject_alloc_len(int len);


/* Function:  hard_link_APIObject(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [APIObject *]
 *
 * Return [UNKN ]  Undocumented return value [APIObject *]
 *
 */
APIObject * hard_link_APIObject(APIObject * obj);


/* Function:  APIObject_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [APIObject *]
 *
 */
APIObject * APIObject_alloc(void);


/* Function:  free_APIObject(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [APIObject *]
 *
 * Return [UNKN ]  Undocumented return value [APIObject *]
 *
 */
APIObject * free_APIObject(APIObject * obj);


/* Function:  add_non_dynAPI(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [dynAPI *]
 * Arg:        add [OWNER] Object to add to the list [APIfunc *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean add_non_dynAPI(dynAPI * obj,APIfunc * add);


/* Function:  flush_non_dynAPI(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [dynAPI *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_non_dynAPI(dynAPI * obj);


/* Function:  add_dynAPI(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [dynAPI *]
 * Arg:        add [OWNER] Object to add to the list [APIObject *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean add_dynAPI(dynAPI * obj,APIObject * add);


/* Function:  flush_dynAPI(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [dynAPI *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_dynAPI(dynAPI * obj);


/* Function:  dynAPI_alloc_std(void)
 *
 * Descrip:    Equivalent to dynAPI_alloc_len(dynAPILISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [dynAPI *]
 *
 */
dynAPI * dynAPI_alloc_std(void);


/* Function:  dynAPI_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [dynAPI *]
 *
 */
dynAPI * dynAPI_alloc_len(int len);


/* Function:  hard_link_dynAPI(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [dynAPI *]
 *
 * Return [UNKN ]  Undocumented return value [dynAPI *]
 *
 */
dynAPI * hard_link_dynAPI(dynAPI * obj);


/* Function:  dynAPI_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [dynAPI *]
 *
 */
dynAPI * dynAPI_alloc(void);


/* Function:  free_dynAPI(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [dynAPI *]
 *
 * Return [UNKN ]  Undocumented return value [dynAPI *]
 *
 */
dynAPI * free_dynAPI(dynAPI * obj);


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
boolean write_pfdoc_dynAPI(dynAPI * api,char * package_name,FILE * ofp);
boolean write_pfdoc_func_def(FuncInfo * fi,FILE * ofp);
boolean write_latex_dynAPI(dynAPI * api,char * module,char * package,FILE * ofp);
boolean write_latex_APIfunc(APIfunc * f,char * package,boolean isobj,char * objectname,FILE * ofp);
boolean write_pod_dynAPI(dynAPI * api,char * module,char * package,FILE * ofp);
void write_pod_obj_APIfunc(APIfunc * f,char * package,char * obj,FILE * ofp);
void write_pod_non_APIfunc(APIfunc * f,char * package,FILE * ofp);
boolean write_type_C_dynAPI(dynAPI * api,char * package_name,FILE * ofp);
boolean write_C_dynAPI(dynAPI * api,char * package_name,FILE * ofp);
boolean write_C_APIfunc(APIfunc * api,char * package_name,FILE * ofp);
boolean is_membasic_type_API(char * type);
boolean is_basic_type_API(char * type);
dynAPI * read_dynAPI_line(char * line,FILE * ifp);
boolean write_perl_XS_accessor_functions(APIObject * obj,char * package,FILE * ofp);
boolean write_dynAPI_accessor_functions(DYNFILE * dfp,dynAPI * api);
boolean write_XS_typemap_dynAPI(dynAPI * api,char * package_name,FILE * ofp);
boolean write_XS_dynAPI(dynAPI * api,char * package_name,FILE * ofp);
boolean write_XS_header_APIObject(dynAPI * api,char * package_name,APIObject * obj,FILE * ofp);
boolean write_XS_APIfunc(APIfunc * fu,char * package_name,FILE * ofp);
boolean reconcile_dynAPI_with_FuncInfo(dynAPI * api,DYNFILE * dfp);
APIObject * APIObject_from_line(char * line,int maxline,FILE * ifp);
APIfunc * APIfunc_from_buffer(char * line);


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
char * c2dyn_type(char * c);
void swap_APIObject(APIfunc ** list,int i,int j) ;
void qsort_APIObject(APIfunc ** list,int left,int right,int (*comp)(APIfunc * ,APIfunc * ));
void sort_APIObject(APIObject * obj,int (*comp)(APIfunc *, APIfunc *));
boolean expand_APIObject(APIObject * obj,int len);
void swap_non_dynAPI(APIfunc ** list,int i,int j) ;
void qsort_non_dynAPI(APIfunc ** list,int left,int right,int (*comp)(APIfunc * ,APIfunc * ));
void sort_non_dynAPI(dynAPI * obj,int (*comp)(APIfunc *, APIfunc *));
boolean expand_non_dynAPI(dynAPI * obj,int len);
void swap_dynAPI(APIObject ** list,int i,int j) ;
void qsort_dynAPI(APIObject ** list,int left,int right,int (*comp)(APIObject * ,APIObject * ));
void sort_dynAPI(dynAPI * obj,int (*comp)(APIObject *, APIObject *));
boolean expand_dynAPI(dynAPI * obj,int len);

#ifdef _cplusplus
}
#endif

#endif
