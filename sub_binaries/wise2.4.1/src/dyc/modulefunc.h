#ifndef DYNAMITEmodulefuncHEADERFILE
#define DYNAMITEmodulefuncHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif

#include "wisebase.h"
#define ModuleFunctionListLISTLENGTH 64




/********************************************************/
/* Wisetools version 3 file                             */
/*                                                      */
/* this file is copyright (c) Ewan Birney 1996          */
/* this file is part of the wisetools sequence analysis */
/* package                                              */
/********************************************************/

/********************************************************/
/* This stores information about which modules do       */
/* have constructors or deconstructors, such that       */
/* they can be made correctly                           */
/*                                                      */
/********************************************************/

/***** RCS Info *****************************************/
/*
   $Id

   $Log

*/
/********************************************************/

struct ModuleFunction {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * name;     
    boolean has_cons;    
    boolean has_decons;  
    boolean has_copy;    
    } ;  
/* ModuleFunction defined */ 
#ifndef DYNAMITE_DEFINED_ModuleFunction
typedef struct ModuleFunction ModuleFunction;
#define DYNAMITE_DEFINED_ModuleFunction
#endif


struct ModuleFunctionList {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    ModuleFunction ** mf;    
    int len;/* len for above mf  */ 
    int maxlen; /* maxlen for above mf */ 
    } ;  
/* ModuleFunctionList defined */ 
#ifndef DYNAMITE_DEFINED_ModuleFunctionList
typedef struct ModuleFunctionList ModuleFunctionList;
#define DYNAMITE_DEFINED_ModuleFunctionList
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_ModuleFunction(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ModuleFunction *]
 *
 * Return [UNKN ]  Undocumented return value [ModuleFunction *]
 *
 */
ModuleFunction * hard_link_ModuleFunction(ModuleFunction * obj);


/* Function:  ModuleFunction_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ModuleFunction *]
 *
 */
ModuleFunction * ModuleFunction_alloc(void);


/* Function:  free_ModuleFunction(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ModuleFunction *]
 *
 * Return [UNKN ]  Undocumented return value [ModuleFunction *]
 *
 */
ModuleFunction * free_ModuleFunction(ModuleFunction * obj);


/* Function:  add_ModuleFunctionList(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [ModuleFunctionList *]
 * Arg:        add [OWNER] Object to add to the list [ModuleFunction *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean add_ModuleFunctionList(ModuleFunctionList * obj,ModuleFunction * add);


/* Function:  flush_ModuleFunctionList(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [ModuleFunctionList *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_ModuleFunctionList(ModuleFunctionList * obj);


/* Function:  ModuleFunctionList_alloc_std(void)
 *
 * Descrip:    Equivalent to ModuleFunctionList_alloc_len(ModuleFunctionListLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ModuleFunctionList *]
 *
 */
ModuleFunctionList * ModuleFunctionList_alloc_std(void);


/* Function:  ModuleFunctionList_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [ModuleFunctionList *]
 *
 */
ModuleFunctionList * ModuleFunctionList_alloc_len(int len);


/* Function:  hard_link_ModuleFunctionList(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ModuleFunctionList *]
 *
 * Return [UNKN ]  Undocumented return value [ModuleFunctionList *]
 *
 */
ModuleFunctionList * hard_link_ModuleFunctionList(ModuleFunctionList * obj);


/* Function:  ModuleFunctionList_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ModuleFunctionList *]
 *
 */
ModuleFunctionList * ModuleFunctionList_alloc(void);


/* Function:  free_ModuleFunctionList(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ModuleFunctionList *]
 *
 * Return [UNKN ]  Undocumented return value [ModuleFunctionList *]
 *
 */
ModuleFunctionList * free_ModuleFunctionList(ModuleFunctionList * obj);


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
ModuleFunction * get_ModuleFunction_from_name(ModuleFunctionList * mfl,char * name);
ModuleFunction * new_ModuleFunction(ModuleFunctionList * mfl,char * name) ;
void show_ModuleFunctionList(ModuleFunctionList * mfl,FILE * ofp);
void show_ModuleFunction(ModuleFunction * mf,FILE * ofp);
char * parse_and_get_module_name_from_func(char * line,boolean isalloc);


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void swap_ModuleFunctionList(ModuleFunction ** list,int i,int j) ;
void qsort_ModuleFunctionList(ModuleFunction ** list,int left,int right,int (*comp)(ModuleFunction * ,ModuleFunction * ));
void sort_ModuleFunctionList(ModuleFunctionList * obj,int (*comp)(ModuleFunction *, ModuleFunction *));
boolean expand_ModuleFunctionList(ModuleFunctionList * obj,int len);

#ifdef _cplusplus
}
#endif

#endif
