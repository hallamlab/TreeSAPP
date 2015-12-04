#ifndef DYNAMITEwisecHEADERFILE
#define DYNAMITEwisecHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif

#include "wisebase.h"
#include "modulefunc.h"
#include "dynfile.h"
#include "objectinfo.h"

#define StructHolderLISTLENGTH 64

#define CKN(s) (s == NULL ? "" : s)

struct StructElement {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * name;     
    char * element_type;     
    char * comment;  
    char * def;  
    boolean isapointer;  
    boolean islinked;    
    boolean islist;  
    boolean ismatrix;    
    char * len_append;   
    boolean isfunc;  
    boolean islocal;     
    boolean ishidden;    
    } ;  
/* StructElement defined */ 
#ifndef DYNAMITE_DEFINED_StructElement
typedef struct StructElement StructElement;
#define DYNAMITE_DEFINED_StructElement
#endif


/* Object StructHolder
 *
 * Descrip: This structure (which is now v. old)
 *        holds the information for a single
 *        'struct' blueprint from dynamite.
 *
 *
 */
struct StructHolder {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    StructElement ** el;     
    int len;/* len for above el  */ 
    int maxlen; /* maxlen for above el */ 
    char * name;     
    char * placing_function;     
    ObjectInfo * oi;     
    } ;  
/* StructHolder defined */ 
#ifndef DYNAMITE_DEFINED_StructHolder
typedef struct StructHolder StructHolder;
#define DYNAMITE_DEFINED_StructHolder
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_StructElement(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [StructElement *]
 *
 * Return [UNKN ]  Undocumented return value [StructElement *]
 *
 */
StructElement * hard_link_StructElement(StructElement * obj);


/* Function:  StructElement_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [StructElement *]
 *
 */
StructElement * StructElement_alloc(void);


/* Function:  free_StructElement(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [StructElement *]
 *
 * Return [UNKN ]  Undocumented return value [StructElement *]
 *
 */
StructElement * free_StructElement(StructElement * obj);


/* Function:  add_StructHolder(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [StructHolder *]
 * Arg:        add [OWNER] Object to add to the list [StructElement *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean add_StructHolder(StructHolder * obj,StructElement * add);


/* Function:  flush_StructHolder(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [StructHolder *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_StructHolder(StructHolder * obj);


/* Function:  StructHolder_alloc_std(void)
 *
 * Descrip:    Equivalent to StructHolder_alloc_len(StructHolderLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [StructHolder *]
 *
 */
StructHolder * StructHolder_alloc_std(void);


/* Function:  StructHolder_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [StructHolder *]
 *
 */
StructHolder * StructHolder_alloc_len(int len);


/* Function:  hard_link_StructHolder(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [StructHolder *]
 *
 * Return [UNKN ]  Undocumented return value [StructHolder *]
 *
 */
StructHolder * hard_link_StructHolder(StructHolder * obj);


/* Function:  StructHolder_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [StructHolder *]
 *
 */
StructHolder * StructHolder_alloc(void);


/* Function:  free_StructHolder(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [StructHolder *]
 *
 * Return [UNKN ]  Undocumented return value [StructHolder *]
 *
 */
StructHolder * free_StructHolder(StructHolder * obj);


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
char * name_from_func_type(char * ty);
void write_StructHolder_header  (DYNFILE * dfp,StructHolder * sh);
void write_StructHolder_function(DYNFILE * dfp,StructHolder * sh,ModuleFunctionList * mfl);
void write_StructHolder_typedef(DYNFILE * dfp,StructHolder * sh,FILE * ofp);
void write_StructElement_typedef(DYNFILE * dfp,StructElement * se);
char * free_function(char * type);
void write_hard_link_function(DYNFILE * dfp,StructHolder * sh);
void write_free_function(DYNFILE * dfp,StructHolder * sh);
boolean is_simple_type(StructElement * se);
boolean is_simple_type_string(char * str);
boolean is_string_type(StructElement * se);
boolean is_binary_dump_type(StructElement * se);
void write_binary_dump_function(DYNFILE * dfp,StructHolder * sh);
void write_binary_read_function(DYNFILE * dfp,StructHolder * sh);
void write_expand_matrix_function(DYNFILE * dfp,StructHolder * sh);
void write_alloc_matrix_function(DYNFILE * dfp,StructHolder * sh);
void write_list_sort_function(DYNFILE * dfp,StructHolder * sh,StructElement * temp);
void write_list_flush_function(DYNFILE * dfp,StructHolder * sh,StructElement * temp);
void write_list_add_function(DYNFILE * dfp,StructHolder * sh,StructElement * temp);
void write_list_expand_function(DYNFILE * dfp,StructHolder * sh,StructElement * temp);
void write_alloc_std_function(DYNFILE * dfp,StructHolder * sh);
void write_alloc_len_function(DYNFILE * dfp,StructHolder * sh);
char * depointer_element(char * element);
boolean isarray(StructElement * se);
void write_simplealloc_function(DYNFILE * dfp,StructHolder * sh);
boolean is_listfunc(StructHolder * sh);
boolean is_matrixfunc(StructHolder * sh);
char * def_from_element(char * str);
char * get_word_from_bang(char * str,char ** word);
boolean read_StructHolder_elements(StructHolder * sh,FILE * ifp);
void show_StructElement(StructElement * se,FILE * ofp);
void show_StructHolder(StructHolder   * sh,FILE * ofp);
StructElement * basic_add_StructElement(StructHolder * sh,char * name,char * element);


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void swap_StructHolder(StructElement ** list,int i,int j) ;
void qsort_StructHolder(StructElement ** list,int left,int right,int (*comp)(StructElement * ,StructElement * ));
void sort_StructHolder(StructHolder * obj,int (*comp)(StructElement *, StructElement *));
boolean expand_StructHolder(StructHolder * obj,int len);

#ifdef _cplusplus
}
#endif

#endif
