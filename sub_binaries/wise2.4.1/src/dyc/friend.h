#ifndef DYNAMITEfriendHEADERFILE
#define DYNAMITEfriendHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dynfile.h"

struct Friend {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * name;     
    } ;  
/* Friend defined */ 
#ifndef DYNAMITE_DEFINED_Friend
typedef struct Friend Friend;
#define DYNAMITE_DEFINED_Friend
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_Friend(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Friend *]
 *
 * Return [UNKN ]  Undocumented return value [Friend *]
 *
 */
Friend * hard_link_Friend(Friend * obj);


/* Function:  Friend_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Friend *]
 *
 */
Friend * Friend_alloc(void);


/* Function:  free_Friend(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Friend *]
 *
 * Return [UNKN ]  Undocumented return value [Friend *]
 *
 */
Friend * free_Friend(Friend * obj);


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
Friend * read_Friend_line(char * line,FILE * ifp);
void write_Friend_header(DYNFILE * dfp,Friend * fr);


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/

#ifdef _cplusplus
}
#endif

#endif
