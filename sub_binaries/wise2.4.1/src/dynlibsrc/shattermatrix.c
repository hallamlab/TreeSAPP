#ifdef _cplusplus
extern "C" {
#endif
#include "shattermatrix.h"

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
# line 32 "shattermatrix.dy"
int fetch_cell_value_ShatterMatrix(ShatterMatrix * sm,int ipos,int jpos,int state)
{
  int * ret;

  ret = fetch_cell_from_ShatterMatrix(sm,ipos,jpos);

  return ret[state];
}

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
# line 44 "shattermatrix.dy"
int * fetch_cell_from_ShatterMatrix(ShatterMatrix * sm,int ipos,int jpos)
{
  int i;
  int offseti;
  int offsetj;

  assert(sm);

  for(i=0;i<sm->len;i++) {
    if( ipos >= sm->smc[i]->starti && ipos < sm->smc[i]->endi &&
	jpos >= sm->smc[i]->startj && jpos < sm->smc[i]->endj ) {
      /* great! */

      /* offset into this */
      offseti = ipos - sm->smc[i]->starti;
      offsetj = jpos - sm->smc[i]->startj;

      return sm->smc[i]->mat[offseti]+(offsetj*sm->cell_length);
    }

  }

  /*
  for(i=0;i<sm->cell_length;i++) {
    if( sm->null_cell[i] !=  -100000 ) {
      fatal("At position %d,%d fetched null cell with %d\n",ipos,jpos,sm->null_cell[i]);
    }
  }
	
  fprintf(stderr," ...null cell %d,%d\n",ipos,jpos);
  */

  return sm->null_cell;
}


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
# line 84 "shattermatrix.dy"
ShatterMatrix * new_ShatterMatrix(DPEnvelope * dpenv,int cell_length,int jlength,int special_length)
{
  int i;
  ShatterMatrix * out;

  out = ShatterMatrix_alloc_len(dpenv->len);
  out->cell_length = cell_length;
  out->null_cell = calloc(cell_length,sizeof(int));
  for(i=0;i<cell_length;i++) {
    out->null_cell[i] = -100000;
  }

  for(i=0;i<dpenv->len;i++) {
    auto DPUnit * dpu = dpenv->dpu[i];
    if( dpu->type == DPENV_RECT ) {

      add_ShatterMatrix(out,new_ShatterMatrixComponent(dpu->starti,dpu->starti+dpu->height,dpu->startj,dpu->startj+dpu->length,cell_length));
    } else if ( dpu->type == DPENV_DIAG ) {
      add_ShatterMatrix(out,new_ShatterMatrixComponent(dpu->starti-dpu->height,dpu->starti+dpu->height+dpu->length,dpu->startj-dpu->height,dpu->startj+dpu->height+dpu->length,cell_length));
    } else {
      fatal("Do not understand type %d",dpu->type);
    }
  }

  out->special = calloc(special_length,sizeof(int*));
  for(i=0;i<special_length;i++) {
    out->special[i] = calloc(jlength,sizeof(int));
  }

  return out;
}


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
# line 120 "shattermatrix.dy"
ShatterMatrixComponent * new_ShatterMatrixComponent(int starti,int endi,int startj,int endj,int cell_length)
{
  ShatterMatrixComponent * out;
  int i;
  int j;
  int leni;

  assert(starti < endi);
  assert(startj < endj);
  assert(cell_length > 0);

  out = ShatterMatrixComponent_alloc();
  out->starti = starti;
  out->endi   = endi;
  out->startj = startj;
  out->endj   = endj;
  leni = endi - starti;

  /*  fprintf(stderr,"Allocating area of %d,%d %d cells\n",leni,(endj-startj),leni*(endj-startj));*/

  out->mat = calloc(leni,sizeof(int *));

  for(i=0;i<leni;i++) {
    out->mat[i] = calloc((endj-startj)*cell_length,sizeof(int));
  }
  
  return out;
}
# line 160 "shattermatrix.c"
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
ShatterMatrixComponent * hard_link_ShatterMatrixComponent(ShatterMatrixComponent * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a ShatterMatrixComponent object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  ShatterMatrixComponent_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ShatterMatrixComponent *]
 *
 */
ShatterMatrixComponent * ShatterMatrixComponent_alloc(void) 
{
    ShatterMatrixComponent * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(ShatterMatrixComponent *) ckalloc (sizeof(ShatterMatrixComponent))) == NULL)    {  
      warn("ShatterMatrixComponent_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->mat = NULL; 
    out->starti = 0; 
    out->startj = 0; 
    out->endi = 0;   
    out->endj = 0;   


    return out;  
}    


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
ShatterMatrixComponent * free_ShatterMatrixComponent(ShatterMatrixComponent * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a ShatterMatrixComponent obj. Should be trappable");    
      return NULL;   
      }  


#ifdef PTHREAD   
    assert(pthread_mutex_lock(&(obj->dynamite_mutex)) == 0); 
#endif   
    if( obj->dynamite_hard_link > 1)     {  
      return_early = 1;  
      obj->dynamite_hard_link--; 
      }  
#ifdef PTHREAD   
    assert(pthread_mutex_unlock(&(obj->dynamite_mutex)) == 0);   
#endif   
    if( return_early == 1)   
      return NULL;   
    if( obj->mat != NULL)    
      ckfree(obj->mat);  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_ShatterMatrix(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_ShatterMatrix
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [ShatterMatrixComponent **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_ShatterMatrix(ShatterMatrixComponent ** list,int i,int j)  
{
    ShatterMatrixComponent * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_ShatterMatrix(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_ShatterMatrix which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [ShatterMatrixComponent **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_ShatterMatrix(ShatterMatrixComponent ** list,int left,int right,int (*comp)(ShatterMatrixComponent * ,ShatterMatrixComponent * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_ShatterMatrix(list,left,(left+right)/2);    
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_ShatterMatrix (list,++last,i);  
      }  
    swap_ShatterMatrix (list,left,last); 
    qsort_ShatterMatrix(list,left,last-1,comp);  
    qsort_ShatterMatrix(list,last+1,right,comp); 
}    


/* Function:  sort_ShatterMatrix(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_ShatterMatrix
 *
 *
 * Arg:         obj [UNKN ] Object containing list [ShatterMatrix *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_ShatterMatrix(ShatterMatrix * obj,int (*comp)(ShatterMatrixComponent *, ShatterMatrixComponent *)) 
{
    qsort_ShatterMatrix(obj->smc,0,obj->len-1,comp); 
    return;  
}    


/* Function:  expand_ShatterMatrix(obj,len)
 *
 * Descrip:    Really an internal function for add_ShatterMatrix
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [ShatterMatrix *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_ShatterMatrix(ShatterMatrix * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_ShatterMatrix called with no need");  
      return TRUE;   
      }  


    if( (obj->smc = (ShatterMatrixComponent ** ) ckrealloc (obj->smc,sizeof(ShatterMatrixComponent *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_ShatterMatrix, returning FALSE");    
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


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
/* will expand function if necessary */ 
boolean add_ShatterMatrix(ShatterMatrix * obj,ShatterMatrixComponent * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_ShatterMatrix(obj,obj->len + ShatterMatrixLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->smc[obj->len++]=add;    
    return TRUE; 
}    


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
int flush_ShatterMatrix(ShatterMatrix * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->smc[i] != NULL)   {  
        free_ShatterMatrixComponent(obj->smc[i]);    
        obj->smc[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  ShatterMatrix_alloc_std(void)
 *
 * Descrip:    Equivalent to ShatterMatrix_alloc_len(ShatterMatrixLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ShatterMatrix *]
 *
 */
ShatterMatrix * ShatterMatrix_alloc_std(void) 
{
    return ShatterMatrix_alloc_len(ShatterMatrixLISTLENGTH); 
}    


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
ShatterMatrix * ShatterMatrix_alloc_len(int len) 
{
    ShatterMatrix * out;/* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = ShatterMatrix_alloc()) == NULL)    
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->smc = (ShatterMatrixComponent ** ) ckcalloc (len,sizeof(ShatterMatrixComponent *))) == NULL)    {  
      warn("Warning, ckcalloc failed in ShatterMatrix_alloc_len");   
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


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
ShatterMatrix * hard_link_ShatterMatrix(ShatterMatrix * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a ShatterMatrix object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  ShatterMatrix_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ShatterMatrix *]
 *
 */
ShatterMatrix * ShatterMatrix_alloc(void) 
{
    ShatterMatrix * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(ShatterMatrix *) ckalloc (sizeof(ShatterMatrix))) == NULL)  {  
      warn("ShatterMatrix_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->smc = NULL; 
    out->len = out->maxlen = 0;  
    out->special = NULL; 
    out->cell_length = 0;    
    out->special_length = 0; 
    out->null_cell = NULL;   


    return out;  
}    


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
ShatterMatrix * free_ShatterMatrix(ShatterMatrix * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a ShatterMatrix obj. Should be trappable"); 
      return NULL;   
      }  


#ifdef PTHREAD   
    assert(pthread_mutex_lock(&(obj->dynamite_mutex)) == 0); 
#endif   
    if( obj->dynamite_hard_link > 1)     {  
      return_early = 1;  
      obj->dynamite_hard_link--; 
      }  
#ifdef PTHREAD   
    assert(pthread_mutex_unlock(&(obj->dynamite_mutex)) == 0);   
#endif   
    if( return_early == 1)   
      return NULL;   
    if( obj->smc != NULL)    {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->smc[i] != NULL) 
          free_ShatterMatrixComponent(obj->smc[i]);  
        }  
      ckfree(obj->smc);  
      }  
    if( obj->special != NULL)    
      ckfree(obj->special);  
    if( obj->null_cell != NULL)  
      ckfree(obj->null_cell);    


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
