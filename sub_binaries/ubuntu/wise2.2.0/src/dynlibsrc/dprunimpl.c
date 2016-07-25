#ifdef _cplusplus
extern "C" {
#endif
#include "dprunimpl.h"


/* Function:  show_help_DPRunImpl(ofp)
 *
 * Descrip:    Shows help functions for DPRunImpl
 *
 *
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 32 "dprunimpl.dy"
void show_help_DPRunImpl(FILE * ofp)
{
  fprintf(ofp,"Dynamic programming matrix implementation\n");
  fprintf(ofp,"  -dymem      memory style [default/linear/explicit]\n");
  fprintf(ofp,"  -kbyte      memory amount to use [4000]\n");
  fprintf(ofp,"  -dydebug    drop into dynamite dp matrix debugger\n");
}

/* Function:  new_DPRunImpl_from_argv(argc,argv)
 *
 * Descrip:    Makes a DPRunImpl object from stripping from
 *             a command line
 *
 *
 * Arg:        argc [UNKN ] Undocumented argument [int *]
 * Arg:        argv [UNKN ] Undocumented argument [char **]
 *
 * Return [UNKN ]  Undocumented return value [DPRunImpl *]
 *
 */
# line 44 "dprunimpl.dy"
DPRunImpl * new_DPRunImpl_from_argv(int * argc,char ** argv)
{
  DPRunImpl * out;
  char * temp;

  out = DPRunImpl_alloc();

  if( (temp = strip_out_assigned_argument(argc,argv,"dymem")) != NULL ) {
    if( strcmp(temp,"explicit") == 0) {
      out->memory = DPIM_Explicit;
    } else if( strcmp(temp,"linear") == 0 ) {
      out->memory = DPIM_Linear;
    } else if( strcmp(temp,"default") == 0 ) {
      out->memory = DPIM_Default;
    } else {
      warn("String [%s] for dynamic memory layout is not recognised",temp);
      free_DPRunImpl(out);
      return NULL;
    }
  }

  if( (temp = strip_out_assigned_argument(argc,argv,"kbyte")) != NULL ) {
    if( is_integer_string(temp,&out->kbyte_size) == FALSE ) {
      warn("String [%s] for dynamic memory size is not recognised",temp);
      free_DPRunImpl(out);
      return NULL;
    }
  }


  if(strip_out_boolean_argument(argc,argv,"dydebug") == TRUE ) {
    out->debug  = 1;
    out->memory = DPIM_Explicit;
  }

  return out;
}

  


# line 73 "dprunimpl.c"
/* Function:  hard_link_DPRunImpl(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [DPRunImpl *]
 *
 */
DPRunImpl * hard_link_DPRunImpl(DPRunImpl * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a DPRunImpl object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  DPRunImpl_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DPRunImpl *]
 *
 */
DPRunImpl * DPRunImpl_alloc(void) 
{
    DPRunImpl * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(DPRunImpl *) ckalloc (sizeof(DPRunImpl))) == NULL)  {  
      warn("DPRunImpl_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
    out->memory = DPIM_Default;  
    out->kbyte_size = 20000; 
    out->debug = FALSE;  


    return out;  
}    


/* Function:  free_DPRunImpl(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [DPRunImpl *]
 *
 */
DPRunImpl * free_DPRunImpl(DPRunImpl * obj) 
{


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a DPRunImpl obj. Should be trappable"); 
      return NULL;   
      }  


    if( obj->dynamite_hard_link > 1)     {  
      obj->dynamite_hard_link--; 
      return NULL;   
      }  


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
