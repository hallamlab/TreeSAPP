#ifdef _cplusplus
extern "C" {
#endif
#include "dpimpl.h"

/* Function:  show_help_DycWarning(ofp)
 *
 * Descrip:    Shows to stdout the list of options
 *             used by DycWarning
 *
 *
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 38 "dpimpl.dy"
void show_help_DycWarning(FILE * ofp)
{
  fprintf(ofp,"Dyc Compiler Warnings\n");
  fprintf(ofp,"  -[no]cwarn   Global switch to put on/off all warnings about C constructs\n");
  fprintf(ofp,"  -[no]extern  Warning about extern scope of names\n");
  fprintf(ofp,"  -[no]methods Warning about extern methods not being scoped/typed\n");
  fprintf(ofp,"  -[no]ctype   Warning about non logical (C) types\n");
}

/* Function:  show_help_DPImplementation(ofp)
 *
 * Descrip:    Shows to stdout the list options used
 *             by DPImplementation
 *
 *
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 51 "dpimpl.dy"
void show_help_DPImplementation(FILE * ofp)
{
  fprintf(ofp,"Dynamic Programming debugging options\n");
  fprintf(ofp,"  -g             generate dynamic programming debugging\n");
  fprintf(ofp,"Dynamic Programming Optimisations\n");
  fprintf(ofp,"  -O             switch all optimisations on\n");
  fprintf(ofp,"  -[no]largemem  Assume large malloc chunks ok (default no)\n");
  fprintf(ofp,"Database search implementation options\n");
  fprintf(ofp,"  -pthreads      generate pthread code\n");
  fprintf(ofp,"  -dbtrace <no>  database trace level code production\n");
  fprintf(ofp,"  -onemodel      generate onemodel (bioxl/g) port\n");
  fprintf(ofp,"Additional generated routines\n");
  fprintf(ofp,"  -prob          all probabilistic routines\n");
  fprintf(ofp,"  -logsum <func> function to use for summing log'd scores\n");

  show_help_DycWarning(ofp);
}

/* Function:  new_DycWarning_from_argstr(argc,argv)
 *
 * Descrip:    Processes the argstring into the DycWarning stuff
 *
 *
 * Arg:        argc [UNKN ] Undocumented argument [int *]
 * Arg:        argv [UNKN ] Undocumented argument [char **]
 *
 * Return [UNKN ]  Undocumented return value [DycWarning *]
 *
 */
# line 72 "dpimpl.dy"
DycWarning * new_DycWarning_from_argstr(int * argc,char ** argv)
{
  DycWarning * out;
  boolean cwarn = TRUE;
  
  out = DycWarning_alloc();

  strip_out_boolean_def_argument(argc,argv,"cwarn",&cwarn);

  if( cwarn == TRUE ) {
    out->warn_extern = TRUE;
    out->warn_extern_method = TRUE;
    out->warn_c_type = TRUE;
  } else {
    out->warn_extern = FALSE;
    out->warn_c_type = FALSE;
    out->warn_extern_method = FALSE;
  }


  strip_out_boolean_def_argument(argc,argv,"extern",&out->warn_extern);
  strip_out_boolean_def_argument(argc,argv,"ctype",&out->warn_c_type);
  strip_out_boolean_def_argument(argc,argv,"methods",&out->warn_extern_method);
  
  return out;
}

/* Function:  new_DPImplementation_from_argstr(argc,argv)
 *
 * Descrip:    Processes the argstring into the DPImplementation
 *             datastructure.
 *
 *
 * Arg:        argc [UNKN ] Undocumented argument [int *]
 * Arg:        argv [UNKN ] Undocumented argument [char **]
 *
 * Return [UNKN ]  Undocumented return value [DPImplementation *]
 *
 */
# line 103 "dpimpl.dy"
DPImplementation * new_DPImplementation_from_argstr(int * argc,char ** argv)
{
  DPImplementation * out;
  char * temp;

  out = DPImplementation_alloc();
  
  if( (strip_out_boolean_argument(argc,argv,"pthreads")) == TRUE ) {
    out->do_threads = TRUE;
  }
  if( (temp=strip_out_assigned_argument(argc,argv,"dbtrace")) != NULL ) {
    if( is_integer_string(temp,&out->db_trace_level) == FALSE ) {
      warn("%s is not an integer argument for dbtrace",temp);
    }
  }

  if( strip_out_boolean_argument(argc,argv,"O") == TRUE ) {
    out->largemem= TRUE;
    /* other optimisations */
  }

  strip_out_boolean_def_argument(argc,argv,"largemem",&out->largemem);

  strip_out_boolean_def_argument(argc,argv,"onemodel",&out->doone);
    
  if( strip_out_boolean_argument(argc,argv,"prob") == TRUE ) {
    out->doprob = TRUE;
  }

  if( strip_out_boolean_argument(argc,argv,"g") == TRUE ) {
    out->dydebug = TRUE;
  }

  if( (temp=strip_out_assigned_argument(argc,argv,"logsum")) != NULL ) {
    out->calcfunc = stringalloc(temp);
  } else {
    out->calcfunc = stringalloc("Probability_logsum");
  }
  out->dycw = new_DycWarning_from_argstr(argc,argv);

  /*  fprintf(stderr,"And %d is extern warning",out->dycw->warn_extern);*/
  return out;
}
  

# line 141 "dpimpl.c"
/* Function:  hard_link_DycWarning(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DycWarning *]
 *
 * Return [UNKN ]  Undocumented return value [DycWarning *]
 *
 */
DycWarning * hard_link_DycWarning(DycWarning * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a DycWarning object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  DycWarning_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DycWarning *]
 *
 */
DycWarning * DycWarning_alloc(void) 
{
    DycWarning * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(DycWarning *) ckalloc (sizeof(DycWarning))) == NULL)    {  
      warn("DycWarning_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->warn_extern = TRUE; 
    out->warn_extern_method = TRUE;  
    out->warn_c_type = TRUE; 


    return out;  
}    


/* Function:  free_DycWarning(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DycWarning *]
 *
 * Return [UNKN ]  Undocumented return value [DycWarning *]
 *
 */
DycWarning * free_DycWarning(DycWarning * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a DycWarning obj. Should be trappable");    
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


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_DPImplementation(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DPImplementation *]
 *
 * Return [UNKN ]  Undocumented return value [DPImplementation *]
 *
 */
DPImplementation * hard_link_DPImplementation(DPImplementation * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a DPImplementation object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  DPImplementation_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DPImplementation *]
 *
 */
DPImplementation * DPImplementation_alloc(void) 
{
    DPImplementation * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(DPImplementation *) ckalloc (sizeof(DPImplementation))) == NULL)    {  
      warn("DPImplementation_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->do_threads = FALSE; 
    out->protect_level = 0;  
    out->db_trace_level = 0; 
    out->doprob = FALSE; 
    out->doone = FALSE;  
    out->calcfunc = NULL;    
    out->largemem = FALSE;   
    out->dycw = NULL;    
    out->dydebug = FALSE;    


    return out;  
}    


/* Function:  free_DPImplementation(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DPImplementation *]
 *
 * Return [UNKN ]  Undocumented return value [DPImplementation *]
 *
 */
DPImplementation * free_DPImplementation(DPImplementation * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a DPImplementation obj. Should be trappable");  
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
    if( obj->calcfunc != NULL)   
      ckfree(obj->calcfunc);     
    if( obj->dycw != NULL)   
      free_DycWarning(obj->dycw);    


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
