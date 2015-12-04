#ifdef _cplusplus
extern "C" {
#endif
#include "basematrix.h"


static int max_matrix_bytes = COMPILE_BASEMATRIX_MAX_KB;

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
# line 130 "basematrix.dy"
void change_max_BaseMatrix_kbytes(int new_kilo_number)
{
  max_matrix_bytes = new_kilo_number;
}

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
# line 139 "basematrix.dy"
int get_max_BaseMatrix_kbytes(void)
{
  return max_matrix_bytes;
}

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
# line 148 "basematrix.dy"
boolean can_make_explicit_matrix(int leni,int lenj,int statesize)
{
  if( leni*lenj*statesize/1024 > max_matrix_bytes)
    return FALSE;
  return TRUE;
}


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
# line 160 "basematrix.dy"
char * basematrix_type_to_string(int type)
{
  switch (type) {
  case  BASEMATRIX_TYPE_UNKNOWN  : return "Unknown";
  case  BASEMATRIX_TYPE_EXPLICIT : return "Explicit";
  case  BASEMATRIX_TYPE_LINEAR   : return "Linear";
  case  BASEMATRIX_TYPE_SHADOW   : return "Shadow";
  default : return "Problem in converting type!";
  }

}

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
# line 182 "basematrix.dy"
BaseMatrix * BaseMatrix_alloc_matrix_and_specials(int len_point,int len_array,int len_spec_point,int len_spec_array)
{
  register int i;
  BaseMatrix *  out;

  /* use dy matrix for main stuff */

  if( (out = BaseMatrix_alloc_matrix(len_point,len_array)) == NULL ) {
    warn("Unable to allocate %d by %d [%d] int positions in basematrix main matrix",len_point,len_array,len_point*len_array);
    return NULL;
  }

  out->spec_len = 0;

  if( (out->specmatrix = (int **) ckcalloc(len_spec_point,sizeof(int *))) == NULL ) {
    warn("Unable to allocate %d special matrix pointers in basematrix",len_spec_point);
    free_BaseMatrix(out);
    return NULL;
  }

  if ((out->specmatrix[0] = (int *)ckcalloc(len_spec_point *  len_spec_array, sizeof(int))) == NULL)
    return NULL;
  for(i = 1; i < len_spec_point; i++)
    out->specmatrix[i] = out->specmatrix[0] + (i * len_spec_array);

  out->spec_len = len_spec_point;

  return out;
}

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
# line 217 "basematrix.dy"
BaseMatrix * free_BaseMatrix(BaseMatrix * obj)
{
  int i;

  if( obj == NULL ) {
    warn("Trying to free NULL basematrix object. Should be trappable");
    return NULL;
  }

  if( obj->dynamite_hard_link > 1 ) {
    obj->dynamite_hard_link--;
    return NULL;
  }
 
  if(obj->matrix != NULL ) {
    for(i=0;i<obj->leni;i++)
      if( obj->matrix[i] != NULL ) {
	ckfree(obj->matrix[i]);
      }
    free(obj->matrix);
  }
  

  if( obj->spec_len > 0 ) {
    if( obj->specmatrix == NULL ) {
      warn("Bad karma. you have a special matrix of length %d, but a NULL specmatrix pointer. I'm not going to free it!",obj->spec_len);
    } else {
      free(obj->specmatrix[0]);
      free(obj->specmatrix);
    } /* end of else */
  } /* end of if specials */

  if( obj->offsetmem != NULL )
    ckfree(obj->offsetmem);
  if( obj->setmem != NULL )
    ckfree(obj->setmem);


  ckfree(obj);

  return NULL;
}


# line 177 "basematrix.c"
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
BaseMatrix * BaseMatrix_alloc_matrix(int leni,int lenj) 
{
    BaseMatrix * out;   /* out is exported          */ 
    register int i; /* for stepping down matrix */ 
    register int j; /* for stepping across matrix */ 


    /* Call alloc function, return NULL if NULL */ 
    if((out = BaseMatrix_alloc()) == NULL)   
      return NULL;   


    /* Allocate memory for matrix  */ 
    if((out->matrix = (int **) ckcalloc (leni,sizeof(int *))) == NULL)   {  
      warn("Memory allocation problem in matrix for BaseMatrix matrix, first pointer set");  
      ckfree(out);   
      return NULL;   
      }  


    /* Add NULL to all matrix pointers so free can be called */ 
    for(i=0;i<leni;i++)  
      out->matrix[i] = NULL;     


    /* Allocate each matrix row */ 
    for(i=0;i<leni;i++)  {  
      out->matrix[i] = (int *) ckcalloc (lenj,sizeof(int ));     
      if( out->matrix[i] == NULL)    {  
        warn("Failed alloc on %d, calling free and returning NULL",i);   
        free_BaseMatrix(out);    
        return NULL;     
        }  
      }  


    for(i=0;i<leni;i++)  {  
      for(j=0;j<lenj;j++)    
        out->matrix[i][j] = 0;   
      }  


    out->leni=out->maxleni=leni;     
    out->lenj=out->maxlenj=lenj;     


    return out;  
}    


/* Function:  expand_BaseMatrix(obj,leni,lenj)
 *
 * Descrip:    Expands matrix. Rarely used
 *
 *
 * Arg:         obj [UNKN ] Undocumented argument [BaseMatrix *]
 * Arg:        leni [UNKN ] Undocumented argument [int]
 * Arg:        lenj [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_BaseMatrix(BaseMatrix * obj,int leni,int lenj) 
{
    int i;   
    int actualj;     


    if( obj == NULL)     {  
      warn("Trying to expand a BaseMatrix but is NULL!");    
      return FALSE;  
      }  


    if( leni <= obj->maxleni && lenj <= obj->maxlenj)    
      return TRUE;   


    if( obj->maxleni < leni )    {  
      if( (obj->matrix=(int **) ckrealloc (obj->matrix,sizeof(int *)*leni)) == NULL) 
        return FALSE;    
      obj->maxleni=obj->leni=leni;   
      }  
    if( lenj > obj->maxlenj )    
      actualj = lenj;    
    else actualj = obj->maxlenj; 
    for(i=0;i<obj->leni;i++) {  
      if((obj->matrix[i] = (int *) realloc (obj->matrix[i],sizeof(int ) * actualj)) == NULL) 
        return FALSE;    
      }  
    return TRUE;     
}    


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
BaseMatrix * hard_link_BaseMatrix(BaseMatrix * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a BaseMatrix object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  BaseMatrix_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [BaseMatrix *]
 *
 */
BaseMatrix * BaseMatrix_alloc(void) 
{
    BaseMatrix * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(BaseMatrix *) ckalloc (sizeof(BaseMatrix))) == NULL)    {  
      warn("BaseMatrix_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->type = BASEMATRIX_TYPE_UNKNOWN; 
    out->matrix = NULL;  
    out->leni=out->maxleni=0;    
    out->lenj=out->maxlenj=0;    
    out->cellsize = 0;   
    out->queryoffset = 0;    
    out->targetoffset = 0;   
    out->spec_len = 0;   
    out->specmatrix = NULL;  
    out->offsetmem = NULL;   
    out->setmem = NULL;  
    out->optimised_shadow = NULL;    


    return out;  
}    



#ifdef _cplusplus
}
#endif
