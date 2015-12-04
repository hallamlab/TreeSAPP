#ifndef DYNAMITEshattermatrixHEADERFILE
#define DYNAMITEshattermatrixHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif

#include "wisebase.h"
#include "dpenvelope.h"

#define ShatterMatrixLISTLENGTH 1024

struct Wise2_ShatterMatrixComponent {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int ** mat;  
    int starti;  
    int startj;  
    int endi;    
    int endj;    
    } ;  
/* ShatterMatrixComponent defined */ 
#ifndef DYNAMITE_DEFINED_ShatterMatrixComponent
typedef struct Wise2_ShatterMatrixComponent Wise2_ShatterMatrixComponent;
#define ShatterMatrixComponent Wise2_ShatterMatrixComponent
#define DYNAMITE_DEFINED_ShatterMatrixComponent
#endif


struct Wise2_ShatterMatrix {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    ShatterMatrixComponent ** smc;   
    int len;/* len for above smc  */ 
    int maxlen; /* maxlen for above smc */ 
    int ** special;  
    int cell_length;     
    int special_length;  
    int * null_cell;     
    } ;  
/* ShatterMatrix defined */ 
#ifndef DYNAMITE_DEFINED_ShatterMatrix
typedef struct Wise2_ShatterMatrix Wise2_ShatterMatrix;
#define ShatterMatrix Wise2_ShatterMatrix
#define DYNAMITE_DEFINED_ShatterMatrix
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  fetch_cell_value_ShatterMatrix(sm,ipos,jpos,state)
 *
 * Descrip:    Gets the actual value from cell
 *
 *
 * Arg:           sm [UNKN ] Undocumented argument [ShatterMatrix *]
 * Arg:         ipos [UNKN ] Undocumented argument [int]
 * Arg:         jpos [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_fetch_cell_value_ShatterMatrix(ShatterMatrix * sm,int ipos,int jpos,int state);
#define fetch_cell_value_ShatterMatrix Wise2_fetch_cell_value_ShatterMatrix


/* Function:  fetch_cell_from_ShatterMatrix(sm,ipos,jpos)
 *
 * Descrip:    gets int * pointer in the right place for this cell
 *
 *
 * Arg:          sm [UNKN ] Undocumented argument [ShatterMatrix *]
 * Arg:        ipos [UNKN ] Undocumented argument [int]
 * Arg:        jpos [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int *]
 *
 */
int * Wise2_fetch_cell_from_ShatterMatrix(ShatterMatrix * sm,int ipos,int jpos);
#define fetch_cell_from_ShatterMatrix Wise2_fetch_cell_from_ShatterMatrix


/* Function:  new_ShatterMatrix(dpenv,cell_length,jlength,special_length)
 *
 * Descrip:    Makes a new shattermatrix - needs to know the cell length 
 *             and special length
 *
 *
 * Arg:                 dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:           cell_length [UNKN ] Undocumented argument [int]
 * Arg:               jlength [UNKN ] Undocumented argument [int]
 * Arg:        special_length [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [ShatterMatrix *]
 *
 */
ShatterMatrix * Wise2_new_ShatterMatrix(DPEnvelope * dpenv,int cell_length,int jlength,int special_length);
#define new_ShatterMatrix Wise2_new_ShatterMatrix


/* Function:  new_ShatterMatrixComponent(starti,endi,startj,endj,cell_length)
 *
 * Descrip:    Makes a new ShatterMatrixComponent - needs to know the cell length and start/end
 *
 *
 * Arg:             starti [UNKN ] Undocumented argument [int]
 * Arg:               endi [UNKN ] Undocumented argument [int]
 * Arg:             startj [UNKN ] Undocumented argument [int]
 * Arg:               endj [UNKN ] Undocumented argument [int]
 * Arg:        cell_length [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [ShatterMatrixComponent *]
 *
 */
ShatterMatrixComponent * Wise2_new_ShatterMatrixComponent(int starti,int endi,int startj,int endj,int cell_length);
#define new_ShatterMatrixComponent Wise2_new_ShatterMatrixComponent


/* Function:  hard_link_ShatterMatrixComponent(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ShatterMatrixComponent *]
 *
 * Return [UNKN ]  Undocumented return value [ShatterMatrixComponent *]
 *
 */
ShatterMatrixComponent * Wise2_hard_link_ShatterMatrixComponent(ShatterMatrixComponent * obj);
#define hard_link_ShatterMatrixComponent Wise2_hard_link_ShatterMatrixComponent


/* Function:  ShatterMatrixComponent_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ShatterMatrixComponent *]
 *
 */
ShatterMatrixComponent * Wise2_ShatterMatrixComponent_alloc(void);
#define ShatterMatrixComponent_alloc Wise2_ShatterMatrixComponent_alloc


/* Function:  free_ShatterMatrixComponent(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ShatterMatrixComponent *]
 *
 * Return [UNKN ]  Undocumented return value [ShatterMatrixComponent *]
 *
 */
ShatterMatrixComponent * Wise2_free_ShatterMatrixComponent(ShatterMatrixComponent * obj);
#define free_ShatterMatrixComponent Wise2_free_ShatterMatrixComponent


/* Function:  add_ShatterMatrix(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [ShatterMatrix *]
 * Arg:        add [OWNER] Object to add to the list [ShatterMatrixComponent *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_ShatterMatrix(ShatterMatrix * obj,ShatterMatrixComponent * add);
#define add_ShatterMatrix Wise2_add_ShatterMatrix


/* Function:  flush_ShatterMatrix(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [ShatterMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_ShatterMatrix(ShatterMatrix * obj);
#define flush_ShatterMatrix Wise2_flush_ShatterMatrix


/* Function:  ShatterMatrix_alloc_std(void)
 *
 * Descrip:    Equivalent to ShatterMatrix_alloc_len(ShatterMatrixLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ShatterMatrix *]
 *
 */
ShatterMatrix * Wise2_ShatterMatrix_alloc_std(void);
#define ShatterMatrix_alloc_std Wise2_ShatterMatrix_alloc_std


/* Function:  ShatterMatrix_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [ShatterMatrix *]
 *
 */
ShatterMatrix * Wise2_ShatterMatrix_alloc_len(int len);
#define ShatterMatrix_alloc_len Wise2_ShatterMatrix_alloc_len


/* Function:  hard_link_ShatterMatrix(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ShatterMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [ShatterMatrix *]
 *
 */
ShatterMatrix * Wise2_hard_link_ShatterMatrix(ShatterMatrix * obj);
#define hard_link_ShatterMatrix Wise2_hard_link_ShatterMatrix


/* Function:  ShatterMatrix_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ShatterMatrix *]
 *
 */
ShatterMatrix * Wise2_ShatterMatrix_alloc(void);
#define ShatterMatrix_alloc Wise2_ShatterMatrix_alloc


/* Function:  free_ShatterMatrix(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ShatterMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [ShatterMatrix *]
 *
 */
ShatterMatrix * Wise2_free_ShatterMatrix(ShatterMatrix * obj);
#define free_ShatterMatrix Wise2_free_ShatterMatrix


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_ShatterMatrix(ShatterMatrixComponent ** list,int i,int j) ;
#define swap_ShatterMatrix Wise2_swap_ShatterMatrix
void Wise2_qsort_ShatterMatrix(ShatterMatrixComponent ** list,int left,int right,int (*comp)(ShatterMatrixComponent * ,ShatterMatrixComponent * ));
#define qsort_ShatterMatrix Wise2_qsort_ShatterMatrix
void Wise2_sort_ShatterMatrix(ShatterMatrix * obj,int (*comp)(ShatterMatrixComponent *, ShatterMatrixComponent *));
#define sort_ShatterMatrix Wise2_sort_ShatterMatrix
boolean Wise2_expand_ShatterMatrix(ShatterMatrix * obj,int len);
#define expand_ShatterMatrix Wise2_expand_ShatterMatrix

#ifdef _cplusplus
}
#endif

#endif
