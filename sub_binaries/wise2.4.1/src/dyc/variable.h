#ifndef DYNAMITEvariableHEADERFILE
#define DYNAMITEvariableHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "wisebase.h"
/********************************************************/
/* Wisetools version 3 file                             */
/*                                                      */
/* this file is copyright (c) Ewan Birney 1996          */
/* this file is part of the wisetools sequence analysis */
/* package                                              */
/********************************************************/

/********************************************************/
/* Variable is a structure which is shared between      */
/* compile and runtime dynamite for managing variables  */
/* for fast code. The data line is only present in      */
/* the run-time version                                 */
/********************************************************/

/***** RCS Info *****************************************/
/*
   $Id: variable.dy,v 1.1.1.1 2001/06/18 13:59:57 birney Exp $

   $Log: variable.dy,v $
   Revision 1.1.1.1  2001/06/18 13:59:57  birney
   moved wise2 to ensembl cvs repository

   Revision 1.1.1.1  1998/08/28 09:30:58  birney
   Wise2

 * Revision 1.1  1996/05/02  15:29:32  birney
 * Initial revision
 *

*/
/********************************************************/


#define ZERO_DIMENSION         152
#define ONE_FIXED_DIMENSION    153
#define ONE_QUERY_DIMENSION    154
#define TWO_FIXED_DIMENSION    155
#define TWO_QUERY_DIMENSION    156


struct MatrixVariable {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * name;     
    int type;    
    int dim1;    
    int dim2;    
    char * source;   
    void * data;     
    } ;  
/* MatrixVariable defined */ 
#ifndef DYNAMITE_DEFINED_MatrixVariable
typedef struct MatrixVariable MatrixVariable;
#define DYNAMITE_DEFINED_MatrixVariable
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_MatrixVariable(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [MatrixVariable *]
 *
 * Return [UNKN ]  Undocumented return value [MatrixVariable *]
 *
 */
MatrixVariable * hard_link_MatrixVariable(MatrixVariable * obj);


/* Function:  MatrixVariable_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [MatrixVariable *]
 *
 */
MatrixVariable * MatrixVariable_alloc(void);


/* Function:  free_MatrixVariable(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [MatrixVariable *]
 *
 * Return [UNKN ]  Undocumented return value [MatrixVariable *]
 *
 */
MatrixVariable * free_MatrixVariable(MatrixVariable * obj);


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/

#ifdef _cplusplus
}
#endif

#endif
