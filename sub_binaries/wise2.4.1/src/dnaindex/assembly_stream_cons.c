#ifdef _cplusplus
extern "C" {
#endif
#include "assembly_stream_cons.h"

# line 24 "assembly_stream_cons.dy"
void show_help_AssemblyStreamConstructor(FILE * ofp)
{
  fprintf(ofp,"Assembly Read Sequence input options\n");
  fprintf(ofp,"  -astype [fasta/sanger] type of input\n");
  fprintf(ofp,"  -asfile <filename>     filename for fasta input\n");
  fprintf(ofp,"  -asdir  <dirname>      directory name for sanger project-style input\n");
  fprintf(ofp,"  -asdirext [1c]         extension for reads in sanger project style, defaul 1c\n");
}

# line 33 "assembly_stream_cons.dy"
AssemblySequenceStream * new_AssemblySequenceStream_from_AssemblyStreamConstructor(AssemblyStreamConstructor * asc)
{
  FILE * ifp; /* a bit evil as we don't close this. Should fix */

  switch(asc->type) {
  case AssemblyStreamTypeFasta :
    if( asc->file_name == NULL ) {
      warn("Fasta file type needs a filename with -asfile");
      return NULL;
    }
    ifp = openfile(asc->file_name,"r");
    if( ifp == NULL ) {
      warn("Unable to open file %s for fasta reading",asc->file_name);
      return NULL;
    }
    return plain_fasta_AssemblySequenceStream(ifp);
    break;

  case AssemblyStreamTypeSanger :
    if( asc->dir_name == NULL || asc->extension == NULL ) {
      warn("Sanger project type needs directory (-asdir) and extension (-asdirext) arguments");
      return NULL;
    }
    return new_sanger_project_AssemblySequenceStream(asc->dir_name,asc->extension);
    break;

  default :
    warn("Unknown assembly stream type %d",asc->type);
  }
    
  fatal("impossible to get here");
  
  return NULL;
}


# line 69 "assembly_stream_cons.dy"
AssemblyStreamConstructor * new_AssemblyStreamConstructor_from_argv(int * argc,char ** argv)
{
  char * temp;
  AssemblyStreamConstructor * out;

  out = AssemblyStreamConstructor_alloc();

  out->type = AssemblyStreamTypeFasta;

  temp = strip_out_assigned_argument(argc,argv,"astype");
  if( temp != NULL) {
    if( strcmp(temp,"fasta") == 0 ) {
      out->type = AssemblyStreamTypeFasta;
    } else if ( strcmp(temp,"sanger") == 0 ) {
      out->type = AssemblyStreamTypeSanger;
    } else {
      warn("Unable to decipher type %s as a type of assembly stream",temp);
    }
  }

  temp = strip_out_assigned_argument(argc,argv,"asdir");
  if( temp != NULL)
    out->dir_name = stringalloc(temp);

  temp = strip_out_assigned_argument(argc,argv,"asdirext");
  if( temp != NULL)
    out->extension = stringalloc(temp);
  else
    out->extension = stringalloc("1c");

  temp = strip_out_assigned_argument(argc,argv,"asfile");
  if( temp != NULL)
    out->file_name = stringalloc(temp);

  return out;
}

# line 90 "assembly_stream_cons.c"
/* Function:  hard_link_AssemblyStreamConstructor(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AssemblyStreamConstructor *]
 *
 * Return [UNKN ]  Undocumented return value [AssemblyStreamConstructor *]
 *
 */
AssemblyStreamConstructor * hard_link_AssemblyStreamConstructor(AssemblyStreamConstructor * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a AssemblyStreamConstructor object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  AssemblyStreamConstructor_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AssemblyStreamConstructor *]
 *
 */
AssemblyStreamConstructor * AssemblyStreamConstructor_alloc(void) 
{
    AssemblyStreamConstructor * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(AssemblyStreamConstructor *) ckalloc (sizeof(AssemblyStreamConstructor))) == NULL)  {  
      warn("AssemblyStreamConstructor_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->type = 0;   
    out->file_name = NULL;   
    out->dir_name = NULL;    
    out->extension = NULL;   


    return out;  
}    


/* Function:  free_AssemblyStreamConstructor(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AssemblyStreamConstructor *]
 *
 * Return [UNKN ]  Undocumented return value [AssemblyStreamConstructor *]
 *
 */
AssemblyStreamConstructor * free_AssemblyStreamConstructor(AssemblyStreamConstructor * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a AssemblyStreamConstructor obj. Should be trappable"); 
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
    if( obj->file_name != NULL)  
      ckfree(obj->file_name);    
    if( obj->dir_name != NULL)   
      ckfree(obj->dir_name);     
    if( obj->extension != NULL)  
      ckfree(obj->extension);    


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
