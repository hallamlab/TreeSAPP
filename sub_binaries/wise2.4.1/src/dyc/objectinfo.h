#ifndef DYNAMITEobjectinfoHEADERFILE
#define DYNAMITEobjectinfoHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dynfile.h"
#include "ftext.h"

struct ObjectInfo {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Ftext * ft; /*  description of Object */ 
    } ;  
/* ObjectInfo defined */ 
#ifndef DYNAMITE_DEFINED_ObjectInfo
typedef struct ObjectInfo ObjectInfo;
#define DYNAMITE_DEFINED_ObjectInfo
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_ObjectInfo(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ObjectInfo *]
 *
 * Return [UNKN ]  Undocumented return value [ObjectInfo *]
 *
 */
ObjectInfo * hard_link_ObjectInfo(ObjectInfo * obj);


/* Function:  ObjectInfo_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ObjectInfo *]
 *
 */
ObjectInfo * ObjectInfo_alloc(void);


/* Function:  free_ObjectInfo(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ObjectInfo *]
 *
 * Return [UNKN ]  Undocumented return value [ObjectInfo *]
 *
 */
ObjectInfo * free_ObjectInfo(ObjectInfo * obj);


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
int write_C_ObjectInfo(ObjectInfo * oi,FILE * ofp);
ObjectInfo * read_ObjectInfo_line_func(char * line,int maxline,FILE * ifp,char * (fgets_func)(char *,int,FILE*));


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/

#ifdef _cplusplus
}
#endif

#endif
