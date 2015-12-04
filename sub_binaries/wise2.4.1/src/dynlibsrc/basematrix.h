#ifndef DYNAMITEbasematrixHEADERFILE
#define DYNAMITEbasematrixHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif


/********************************************************/
/* Wisetools version 3 file                             */
/*                                                      */
/* this file is copyright (c) Ewan Birney 1996          */
/* this file is part of the wisetools sequence analysis */
/* package                                              */
/********************************************************/

/********************************************************/
/* BaseMatrix - structure to hold the memory of         */
/* dynamite made matrices.                              */
/* A dynamite support library                           */
/*                                                      */
/********************************************************/

/***** RCS Info *****************************************/
/*
   $Id: basematrix.dy,v 1.3 2006/03/26 18:45:43 birney Exp $

   $Log: basematrix.dy,v $
   Revision 1.3  2006/03/26 18:45:43  birney
   new basematrix memory layout

   Revision 1.2  2006/02/27 09:32:23  birney
   start of folding in steves fixes into dynamite

   Revision 1.1.1.1  2001/06/18 13:59:59  birney
   moved wise2 to ensembl cvs repository

   Revision 1.2  1999/03/08 10:47:48  birney
   added in fixes for a variety of things

   Revision 1.1.1.1  1998/08/28 09:30:47  birney
   Wise2

 * Revision 1.9  1997/11/12  13:20:18  birney
 * added documentation, cleaned up
 *
 * Revision 1.8  1997/07/24  15:34:45  birney
 * changed matrix size, and made 'IMPOSSIBLY_HIGH_SCORE'
 *
 * Revision 1.7  1997/02/11  20:42:48  birney
 * added max basematrix number and query size for it
 *
 * Revision 1.6  1997/02/03  16:04:48  birney
 * forgot to set spec_len in _and_special function. Ooops!
 *
 * Revision 1.5  1997/02/03  16:01:30  birney
 * added BaseMatrix_alloc_matrix_and_specials, but not tested.
 *
 * Revision 1.4  1996/11/12  16:41:54  birney
 * deconstructor was trying to free NULL'd pointers. (bad boy!)
 *
 * Revision 1.3  1996/11/11  22:17:22  birney
 * tweaked specials to allow clean linear memory. Hacky deconstructor
 * though now...
 *
 * Revision 1.2  1996/10/15  13:00:52  birney
 * changed memory behaviour for offsets: caught a memory leak
 *
 * Revision 1.1  1996/03/03  15:38:44  birney
 * Initial revision
 *

*/
/********************************************************/



#include "wisebase.h"

enum basematrix_types {
  BASEMATRIX_TYPE_UNKNOWN  = 67,
  BASEMATRIX_TYPE_EXPLICIT,
  BASEMATRIX_TYPE_LINEAR,
  BASEMATRIX_TYPE_SHADOW,
  BASEMATRIX_TYPE_VERYSMALL
};

/* 20 MB assummed. Could break something? */
#define COMPILE_BASEMATRIX_MAX_KB 20000


#define IMPOSSIBLY_HIGH_SCORE 500000

struct Wise2_BaseMatrix {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int type;    
    int ** matrix;  /*  NB i,j proper */ 
    int leni;   /* leni for above matrix  */ 
    int maxleni;/* max length for above pointer set */ 
    int lenj;   /* lenj for above matrix  */ 
    int maxlenj;/* max length for above pointer set */ 
    int cellsize;    
    int queryoffset;     
    int targetoffset;    
    int spec_len;    
    int ** offsetmatrix;    /*  complex - this points to an internal matrix of correct offsets (!) */ 
    int ** specmatrix;  /*  no longer linked: we have this memory specific... */ 
    int ** offsetmem;    
    int ** setmem;   
    int *  optimised_shadow;     
    } ;  
/* BaseMatrix defined */ 
#ifndef DYNAMITE_DEFINED_BaseMatrix
typedef struct Wise2_BaseMatrix Wise2_BaseMatrix;
#define BaseMatrix Wise2_BaseMatrix
#define DYNAMITE_DEFINED_BaseMatrix
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  change_max_BaseMatrix_kbytes(new_kilo_number)
 *
 * Descrip:    This is to change, at run-time the maximum level of bytes basematrix *thinks*
 *             it can use. This number is *not* used for any actual calls to basematrix
 *             allocation: it is only used with /get_max_BaseMatrix_kbytes
 *
 *
 * Arg:        new_kilo_number [UNKN ] max kilobytes allowed [int]
 *
 */
void Wise2_change_max_BaseMatrix_kbytes(int new_kilo_number);
#define change_max_BaseMatrix_kbytes Wise2_change_max_BaseMatrix_kbytes


/* Function:  get_max_BaseMatrix_kbytes(void)
 *
 * Descrip:    returns the max. number of kilobytes suggested as a limited
 *             to BaseMatrix. 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_get_max_BaseMatrix_kbytes(void);
#define get_max_BaseMatrix_kbytes Wise2_get_max_BaseMatrix_kbytes


/* Function:  can_make_explicit_matrix(leni,lenj,statesize)
 *
 * Descrip:    Just checkes that leni*lenj*statesize/1024 < max_matrix_bytes.
 *             returns TRUE if so, FALSE if not
 *
 *
 * Arg:             leni [UNKN ] Undocumented argument [int]
 * Arg:             lenj [UNKN ] Undocumented argument [int]
 * Arg:        statesize [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_can_make_explicit_matrix(int leni,int lenj,int statesize);
#define can_make_explicit_matrix Wise2_can_make_explicit_matrix


/* Function:  basematrix_type_to_string(type)
 *
 * Descrip:    turns a int type to a char string of 'printable'
 *             types.
 *
 *
 * Arg:        type [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [char *]
 *
 */
char * Wise2_basematrix_type_to_string(int type);
#define basematrix_type_to_string Wise2_basematrix_type_to_string


/* Function:  BaseMatrix_alloc_matrix_and_specials(len_spec_poin,len_point,len_array,len_spec_point,len_spec_array)
 *
 * Descrip:    This function allocates the two bits of
 *             matrix memory, of course returning a decent 
 *             NULL (with memory zapped) if it can't do it
 *
 *
 * Arg:         len_spec_poin [UNKN ] length of pointers in special matrix [NullString]
 * Arg:             len_point [UNKN ] length of pointers in main matrix [int]
 * Arg:             len_array [UNKN ] length of array in main matrix [int]
 * Arg:        len_spec_point [UNKN ] Undocumented argument [int]
 * Arg:        len_spec_array [UNKN ] length of special array [int]
 *
 * Return [UNKN ]  Undocumented return value [BaseMatrix *]
 *
 */
BaseMatrix * Wise2_BaseMatrix_alloc_matrix_and_specials(int len_point,int len_array,int len_spec_point,int len_spec_array);
#define BaseMatrix_alloc_matrix_and_specials Wise2_BaseMatrix_alloc_matrix_and_specials


/* Function:  free_BaseMatrix(obj)
 *
 * Descrip:    this is the override deconstructor for basematrix. It will
 *             free both matrix and special memory
 *
 *
 * Arg:        obj [UNKN ] Undocumented argument [BaseMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [BaseMatrix *]
 *
 */
BaseMatrix * Wise2_free_BaseMatrix(BaseMatrix * obj);
#define free_BaseMatrix Wise2_free_BaseMatrix


/* Function:  BaseMatrix_alloc_matrix(leni,lenj)
 *
 * Descrip:    Allocates structure and matrix
 *
 *
 * Arg:        leni [UNKN ] Length of first dimension of matrix [int]
 * Arg:        lenj [UNKN ] Length of second dimension of matrix [int]
 *
 * Return [UNKN ]  Undocumented return value [BaseMatrix *]
 *
 */
BaseMatrix * Wise2_BaseMatrix_alloc_matrix(int leni,int lenj);
#define BaseMatrix_alloc_matrix Wise2_BaseMatrix_alloc_matrix


/* Function:  hard_link_BaseMatrix(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [BaseMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [BaseMatrix *]
 *
 */
BaseMatrix * Wise2_hard_link_BaseMatrix(BaseMatrix * obj);
#define hard_link_BaseMatrix Wise2_hard_link_BaseMatrix


/* Function:  BaseMatrix_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [BaseMatrix *]
 *
 */
BaseMatrix * Wise2_BaseMatrix_alloc(void);
#define BaseMatrix_alloc Wise2_BaseMatrix_alloc


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
boolean Wise2_expand_BaseMatrix(BaseMatrix * obj,int leni,int lenj);
#define expand_BaseMatrix Wise2_expand_BaseMatrix

#ifdef _cplusplus
}
#endif

#endif
