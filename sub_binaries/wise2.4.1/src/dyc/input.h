#ifndef DYNAMITEinputHEADERFILE
#define DYNAMITEinputHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "wisebase.h"

#define CodeLISTLENGTH 128
#define InputLISTLENGTH 16
struct Code {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char ** lines;   
    int len;/* len for above lines  */ 
    int maxlen; /* maxlen for above lines */ 
    } ;  
/* Code defined */ 
#ifndef DYNAMITE_DEFINED_Code
typedef struct Code Code;
#define DYNAMITE_DEFINED_Code
#endif


struct InDeclare {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * name;     
    char * type;     
    } ;  
/* InDeclare defined */ 
#ifndef DYNAMITE_DEFINED_InDeclare
typedef struct InDeclare InDeclare;
#define DYNAMITE_DEFINED_InDeclare
#endif


struct InRequire {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * name;     
    char * type;     
    char * help;     
    char * tag;  
    } ;  
/* InRequire defined */ 
#ifndef DYNAMITE_DEFINED_InRequire
typedef struct InRequire InRequire;
#define DYNAMITE_DEFINED_InRequire
#endif


struct Input {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    InDeclare ** decl;   
    int decl_len;   /* len for above decl  */ 
    int decl_maxlen;/* maxlen for above decl */ 
    InRequire ** req;    
    int req_len;/* len for above req  */ 
    int req_maxlen; /* maxlen for above req */ 
    Code * code;     
    } ;  
/* Input defined */ 
#ifndef DYNAMITE_DEFINED_Input
typedef struct Input Input;
#define DYNAMITE_DEFINED_Input
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  add_Code(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Code *]
 * Arg:        add [OWNER] Object to add to the list [char *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean add_Code(Code * obj,char * add);


/* Function:  flush_Code(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [Code *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_Code(Code * obj);


/* Function:  Code_alloc_std(void)
 *
 * Descrip:    Equivalent to Code_alloc_len(CodeLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Code *]
 *
 */
Code * Code_alloc_std(void);


/* Function:  Code_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [Code *]
 *
 */
Code * Code_alloc_len(int len);


/* Function:  hard_link_Code(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Code *]
 *
 * Return [UNKN ]  Undocumented return value [Code *]
 *
 */
Code * hard_link_Code(Code * obj);


/* Function:  Code_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Code *]
 *
 */
Code * Code_alloc(void);


/* Function:  free_Code(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Code *]
 *
 * Return [UNKN ]  Undocumented return value [Code *]
 *
 */
Code * free_Code(Code * obj);


/* Function:  hard_link_InDeclare(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [InDeclare *]
 *
 * Return [UNKN ]  Undocumented return value [InDeclare *]
 *
 */
InDeclare * hard_link_InDeclare(InDeclare * obj);


/* Function:  InDeclare_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [InDeclare *]
 *
 */
InDeclare * InDeclare_alloc(void);


/* Function:  free_InDeclare(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [InDeclare *]
 *
 * Return [UNKN ]  Undocumented return value [InDeclare *]
 *
 */
InDeclare * free_InDeclare(InDeclare * obj);


/* Function:  hard_link_InRequire(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [InRequire *]
 *
 * Return [UNKN ]  Undocumented return value [InRequire *]
 *
 */
InRequire * hard_link_InRequire(InRequire * obj);


/* Function:  InRequire_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [InRequire *]
 *
 */
InRequire * InRequire_alloc(void);


/* Function:  free_InRequire(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [InRequire *]
 *
 * Return [UNKN ]  Undocumented return value [InRequire *]
 *
 */
InRequire * free_InRequire(InRequire * obj);


/* Function:  add_decl_Input(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Input *]
 * Arg:        add [OWNER] Object to add to the list [InDeclare *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean add_decl_Input(Input * obj,InDeclare * add);


/* Function:  flush_decl_Input(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [Input *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_decl_Input(Input * obj);


/* Function:  add_req_Input(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Input *]
 * Arg:        add [OWNER] Object to add to the list [InRequire *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean add_req_Input(Input * obj,InRequire * add);


/* Function:  flush_req_Input(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [Input *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_req_Input(Input * obj);


/* Function:  Input_alloc_std(void)
 *
 * Descrip:    Equivalent to Input_alloc_len(InputLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Input *]
 *
 */
Input * Input_alloc_std(void);


/* Function:  Input_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [Input *]
 *
 */
Input * Input_alloc_len(int len);


/* Function:  hard_link_Input(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Input *]
 *
 * Return [UNKN ]  Undocumented return value [Input *]
 *
 */
Input * hard_link_Input(Input * obj);


/* Function:  Input_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Input *]
 *
 */
Input * Input_alloc(void);


/* Function:  free_Input(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Input *]
 *
 * Return [UNKN ]  Undocumented return value [Input *]
 *
 */
Input * free_Input(Input * obj);


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
Input * read_Input_line(char * line,FILE * ifp);
Code * read_Code_line(char * line,FILE * ifp);
InRequire * InRequire_from_line(char * line);
InDeclare * InDeclare_from_line(char * line);


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void swap_Code(char ** list,int i,int j) ;
void qsort_Code(char ** list,int left,int right,int (*comp)(char * ,char * ));
void sort_Code(Code * obj,int (*comp)(char *, char *));
boolean expand_Code(Code * obj,int len);
void swap_decl_Input(InDeclare ** list,int i,int j) ;
void qsort_decl_Input(InDeclare ** list,int left,int right,int (*comp)(InDeclare * ,InDeclare * ));
void sort_decl_Input(Input * obj,int (*comp)(InDeclare *, InDeclare *));
boolean expand_decl_Input(Input * obj,int len);
void swap_req_Input(InRequire ** list,int i,int j) ;
void qsort_req_Input(InRequire ** list,int left,int right,int (*comp)(InRequire * ,InRequire * ));
void sort_req_Input(Input * obj,int (*comp)(InRequire *, InRequire *));
boolean expand_req_Input(Input * obj,int len);

#ifdef _cplusplus
}
#endif

#endif
