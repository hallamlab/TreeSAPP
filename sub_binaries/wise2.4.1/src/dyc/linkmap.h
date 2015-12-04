#ifndef DYNAMITElinkmapHEADERFILE
#define DYNAMITElinkmapHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "wisebase.h"

struct LinkMap {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * word;     
    char * urltail;  
    char * text;     
    } ;  
/* LinkMap defined */ 
#ifndef DYNAMITE_DEFINED_LinkMap
typedef struct LinkMap LinkMap;
#define DYNAMITE_DEFINED_LinkMap
#endif


struct LinkSet {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    LinkMap ** lm;   
    int len;/* len for above lm  */ 
    int maxlen; /* maxlen for above lm */ 
    } ;  
/* LinkSet defined */ 
#ifndef DYNAMITE_DEFINED_LinkSet
typedef struct LinkSet LinkSet;
#define DYNAMITE_DEFINED_LinkSet
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_LinkMap(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [LinkMap *]
 *
 * Return [UNKN ]  Undocumented return value [LinkMap *]
 *
 */
LinkMap * hard_link_LinkMap(LinkMap * obj);


/* Function:  LinkMap_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LinkMap *]
 *
 */
LinkMap * LinkMap_alloc(void);


/* Function:  free_LinkMap(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [LinkMap *]
 *
 * Return [UNKN ]  Undocumented return value [LinkMap *]
 *
 */
LinkMap * free_LinkMap(LinkMap * obj);


/* Function:  add_LinkSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [LinkSet *]
 * Arg:        add [OWNER] Object to add to the list [LinkMap *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean add_LinkSet(LinkSet * obj,LinkMap * add);


/* Function:  flush_LinkSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [LinkSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_LinkSet(LinkSet * obj);


/* Function:  LinkSet_alloc_std(void)
 *
 * Descrip:    Equivalent to LinkSet_alloc_len(LinkSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LinkSet *]
 *
 */
LinkSet * LinkSet_alloc_std(void);


/* Function:  LinkSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [LinkSet *]
 *
 */
LinkSet * LinkSet_alloc_len(int len);


/* Function:  hard_link_LinkSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [LinkSet *]
 *
 * Return [UNKN ]  Undocumented return value [LinkSet *]
 *
 */
LinkSet * hard_link_LinkSet(LinkSet * obj);


/* Function:  LinkSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LinkSet *]
 *
 */
LinkSet * LinkSet_alloc(void);


/* Function:  free_LinkSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [LinkSet *]
 *
 * Return [UNKN ]  Undocumented return value [LinkSet *]
 *
 */
LinkSet * free_LinkSet(LinkSet * obj);


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
char * update_links_LinkMap(char * line,char * startid,char * stopid,char * nolinkstr,LinkSet * ls);
LinkMap * new_std_LinkMap(char * word,char * urlstub,char * textstub);
LinkMap * LinkMap_from_word(char * word,LinkSet * ls);


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void swap_LinkSet(LinkMap ** list,int i,int j) ;
void qsort_LinkSet(LinkMap ** list,int left,int right,int (*comp)(LinkMap * ,LinkMap * ));
void sort_LinkSet(LinkSet * obj,int (*comp)(LinkMap *, LinkMap *));
boolean expand_LinkSet(LinkSet * obj,int len);

#ifdef _cplusplus
}
#endif

#endif
