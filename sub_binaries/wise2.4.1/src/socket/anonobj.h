#ifndef DYNAMITEanonobjHEADERFILE
#define DYNAMITEanonobjHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "wisebase.h"

#define AnonymousObjectListLISTLENGTH 16


struct Wise2_AnonymousObject {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    void * obj;  
    char * type;     
    void (*free_object)(void *); 
    } ;  
/* AnonymousObject defined */ 
#ifndef DYNAMITE_DEFINED_AnonymousObject
typedef struct Wise2_AnonymousObject Wise2_AnonymousObject;
#define AnonymousObject Wise2_AnonymousObject
#define DYNAMITE_DEFINED_AnonymousObject
#endif


struct Wise2_AnonymousObjectList {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    AnonymousObject ** anon;     
    int len;/* len for above anon  */ 
    int maxlen; /* maxlen for above anon */ 
    } ;  
/* AnonymousObjectList defined */ 
#ifndef DYNAMITE_DEFINED_AnonymousObjectList
typedef struct Wise2_AnonymousObjectList Wise2_AnonymousObjectList;
#define AnonymousObjectList Wise2_AnonymousObjectList
#define DYNAMITE_DEFINED_AnonymousObjectList
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  free_AnonymousObject(ao)
 *
 * Descrip:    frees an anonymous object
 *
 *
 * Arg:        ao [UNKN ] Undocumented argument [AnonymousObject *]
 *
 * Return [UNKN ]  Undocumented return value [AnonymousObject *]
 *
 */
AnonymousObject * Wise2_free_AnonymousObject(AnonymousObject * ao);
#define free_AnonymousObject Wise2_free_AnonymousObject


/* Function:  hard_link_AnonymousObject(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AnonymousObject *]
 *
 * Return [UNKN ]  Undocumented return value [AnonymousObject *]
 *
 */
AnonymousObject * Wise2_hard_link_AnonymousObject(AnonymousObject * obj);
#define hard_link_AnonymousObject Wise2_hard_link_AnonymousObject


/* Function:  AnonymousObject_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AnonymousObject *]
 *
 */
AnonymousObject * Wise2_AnonymousObject_alloc(void);
#define AnonymousObject_alloc Wise2_AnonymousObject_alloc


/* Function:  add_AnonymousObjectList(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [AnonymousObjectList *]
 * Arg:        add [OWNER] Object to add to the list [AnonymousObject *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_AnonymousObjectList(AnonymousObjectList * obj,AnonymousObject * add);
#define add_AnonymousObjectList Wise2_add_AnonymousObjectList


/* Function:  flush_AnonymousObjectList(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [AnonymousObjectList *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_AnonymousObjectList(AnonymousObjectList * obj);
#define flush_AnonymousObjectList Wise2_flush_AnonymousObjectList


/* Function:  AnonymousObjectList_alloc_std(void)
 *
 * Descrip:    Equivalent to AnonymousObjectList_alloc_len(AnonymousObjectListLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AnonymousObjectList *]
 *
 */
AnonymousObjectList * Wise2_AnonymousObjectList_alloc_std(void);
#define AnonymousObjectList_alloc_std Wise2_AnonymousObjectList_alloc_std


/* Function:  AnonymousObjectList_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [AnonymousObjectList *]
 *
 */
AnonymousObjectList * Wise2_AnonymousObjectList_alloc_len(int len);
#define AnonymousObjectList_alloc_len Wise2_AnonymousObjectList_alloc_len


/* Function:  hard_link_AnonymousObjectList(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AnonymousObjectList *]
 *
 * Return [UNKN ]  Undocumented return value [AnonymousObjectList *]
 *
 */
AnonymousObjectList * Wise2_hard_link_AnonymousObjectList(AnonymousObjectList * obj);
#define hard_link_AnonymousObjectList Wise2_hard_link_AnonymousObjectList


/* Function:  AnonymousObjectList_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AnonymousObjectList *]
 *
 */
AnonymousObjectList * Wise2_AnonymousObjectList_alloc(void);
#define AnonymousObjectList_alloc Wise2_AnonymousObjectList_alloc


/* Function:  free_AnonymousObjectList(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AnonymousObjectList *]
 *
 * Return [UNKN ]  Undocumented return value [AnonymousObjectList *]
 *
 */
AnonymousObjectList * Wise2_free_AnonymousObjectList(AnonymousObjectList * obj);
#define free_AnonymousObjectList Wise2_free_AnonymousObjectList


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_AnonymousObjectList(AnonymousObject ** list,int i,int j) ;
#define swap_AnonymousObjectList Wise2_swap_AnonymousObjectList
void Wise2_qsort_AnonymousObjectList(AnonymousObject ** list,int left,int right,int (*comp)(AnonymousObject * ,AnonymousObject * ));
#define qsort_AnonymousObjectList Wise2_qsort_AnonymousObjectList
void Wise2_sort_AnonymousObjectList(AnonymousObjectList * obj,int (*comp)(AnonymousObject *, AnonymousObject *));
#define sort_AnonymousObjectList Wise2_sort_AnonymousObjectList
boolean Wise2_expand_AnonymousObjectList(AnonymousObjectList * obj,int len);
#define expand_AnonymousObjectList Wise2_expand_AnonymousObjectList

#ifdef _cplusplus
}
#endif

#endif
