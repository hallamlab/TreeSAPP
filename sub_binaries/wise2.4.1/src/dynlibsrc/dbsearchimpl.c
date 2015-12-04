#ifdef _cplusplus
extern "C" {
#endif
#include "dbsearchimpl.h"


/* Function:  impl_string_DBSearchImpl(dbsi)
 *
 * Descrip:    Gets a static text string out of the
 *             search implementation 
 *
 *
 * Arg:        dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 *
 * Return [SOFT ]  string of the search implementation [char *]
 *
 */
# line 89 "dbsearchimpl.dy"
char * impl_string_DBSearchImpl(DBSearchImpl * dbsi)
{
  switch(dbsi->type) {
  case DBSearchImpl_Serial : 
    return "Single Threaded processor (serial)";
  case DBSearchImpl_Pthreads :
    return "Pthreaded multiple processor";
  case DBSearchImpl_MPI :
    return "MPI distributed memory multiple processor";
  case DBSearchImpl_PVM :
    return "Parallel Virtual Machine implementation";
  default :
    return "Unknown implementation";
  }
}

/* Function:  impl_string_routine_DBSearchImpl(dbsi)
 *
 * Descrip:    Gets a static text string out of the
 *             search routine implementation 
 *
 *
 * Arg:        dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 *
 * Return [SOFT ]  string of the search implementation [char *]
 *
 */
# line 112 "dbsearchimpl.dy"
char * impl_string_routine_DBSearchImpl(DBSearchImpl * dbsi)
{
  switch(dbsi->search_routine) {
  case DBSearchImplRoutine_Exact : 
    return "Exact calculation";
  case DBSearchImplRoutine_Kbest :
    return "Kbest heuristic";
  case DBSearchImplRoutine_Forward :
    return "Forward score";
  default :
    return "Unknown routine implementation";
  }
}

/* Function:  number_of_threads_DBSearchImpl(dbsi)
 *
 * Descrip:    This gets out the recommended number of
 *             threads for a dbsi
 *
 *
 * Arg:        dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 130 "dbsearchimpl.dy"
int number_of_threads_DBSearchImpl(DBSearchImpl * dbsi) 
{
  int num;

  if( dbsi->type != DBSearchImpl_Pthreads ) {
    warn("Asking a database implementation how many threads when it is not pthreads but [%s]",impl_string_DBSearchImpl(dbsi));
  }

  if( dbsi->suggest_thread_no == -1) {
    num = DBSEARCH_NPROC;
    if (num == -1) 
      num = 2;    

    return num;
  }

  return dbsi->suggest_thread_no;
}
	 

/* Function:  new_MPI_DBSearchImpl(void)
 *
 * Descrip:    Makes a new MPI DBSearchImpl
 *
 *             For use mainly for api's who don't want
 *             to initalize the object from the command
 *             line
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DBSearchImpl *]
 *
 */
# line 157 "dbsearchimpl.dy"
DBSearchImpl * new_MPI_DBSearchImpl(void)
{
  DBSearchImpl * out;

  out = DBSearchImpl_alloc();
  out->type = DBSearchImpl_MPI;
  return out;
}

/* Function:  new_pthread_DBSearchImpl(void)
 *
 * Descrip:    Makes a new pthreaded DBSearchImpl
 *
 *             For use mainly for api's who don't want
 *             to initalize the object from the command
 *             line
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DBSearchImpl *]
 *
 */
# line 173 "dbsearchimpl.dy"
DBSearchImpl * new_pthread_DBSearchImpl(void)
{
  DBSearchImpl * out;

  out = DBSearchImpl_alloc();
  out->type = DBSearchImpl_Pthreads;
  return out;
}

/* Function:  new_serial_DBSearchImpl(void)
 *
 * Descrip:    Makes a new serial DBSearchImpl
 *
 *             For use mainly for api's who don't want
 *             to initalize the object from the command
 *             line
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DBSearchImpl *]
 *
 */
# line 189 "dbsearchimpl.dy"
DBSearchImpl * new_serial_DBSearchImpl(void)
{
  DBSearchImpl * out;

  out = DBSearchImpl_alloc();
  out->type = DBSearchImpl_Serial;
  return out;
}

/* Function:  show_help_DBSearchImpl(ofp)
 *
 * Descrip:    This shows the help for the search implementation
 *             system. 
 *
 *             It prints out lines like
 *               -pthreads  use pthreaded code
 *
 *
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 205 "dbsearchimpl.dy"
void show_help_DBSearchImpl(FILE * ofp)
{
  fprintf(ofp,"Database search implementation\n");
  fprintf(ofp,"  -serial       use serial code (single processor)\n");
  fprintf(ofp,"  -pthread      use pthread code (SMP box)\n");
  fprintf(ofp,"  -pthr_no <no> Number of threads to use\n");
  fprintf(ofp,"  -mpi          use MPI code (distributed memory system)\n");
  fprintf(ofp,"  -pvm          use parallel virtual machine search system\n");
  fprintf(ofp,"  -dbtrace <no> Trace level of the database code (for debugging only)\n");
  fprintf(ofp,"  -sroutine <type> Search type routine [exact/kbest/forward]\n");
}


/* Function:  new_DBSearchImpl_from_argv(argc,argv)
 *
 * Descrip:    This process command line arguments to
 *             make a new DBSearchImpl. This way the
 *             search implementation is relatively
 *             flexible from the program which calls it.
 *
 *
 * Arg:        argc [UNKN ] Undocumented argument [int *]
 * Arg:        argv [UNKN ] Undocumented argument [char **]
 *
 * Return [UNKN ]  Undocumented return value [DBSearchImpl *]
 *
 */
# line 224 "dbsearchimpl.dy"
DBSearchImpl * new_DBSearchImpl_from_argv(int * argc,char ** argv)
{
  DBSearchImpl * out;
  char * temp;

  out = DBSearchImpl_alloc();
  out->trace_file = stderr;
  if( (strip_out_boolean_argument(argc,argv,"mpi")) == TRUE ) {
    out->type = DBSearchImpl_MPI;
  }

  if( (strip_out_boolean_argument(argc,argv,"pthread")) == TRUE ) {
    out->type = DBSearchImpl_Pthreads;
  }

  if( (strip_out_boolean_argument(argc,argv,"pvm")) == TRUE ) {
    out->type = DBSearchImpl_PVM;
  }

  if( (strip_out_boolean_argument(argc,argv,"serial")) == TRUE ) {
    out->type = DBSearchImpl_Serial;
  }
  
  if( (temp = strip_out_assigned_argument(argc,argv,"pthr_no")) != NULL ) {
    if( is_integer_string(temp,&out->suggest_thread_no) == FALSE ) {
      warn("String [%s] for pthr_no is not an integer!",temp);
      free_DBSearchImpl(out);
      return NULL;
    }
  }

  if( (temp = strip_out_assigned_argument(argc,argv,"dbtrace")) != NULL ) {
    if( is_integer_string(temp,&out->trace_level) == FALSE ) {
      warn("String [%s] for dbtrace is not an integer!",temp);
      free_DBSearchImpl(out);
      return NULL;
    }
  }

  if( (temp = strip_out_assigned_argument(argc,argv,"sroutine")) != NULL ) {
    if( strcmp(temp,"exact") == 0) {
      out->search_routine = DBSearchImplRoutine_Exact;
    } else if( strcmp(temp,"kbest") == 0 ) {
      out->search_routine = DBSearchImplRoutine_Kbest;
    } else if( strcmp(temp,"forward") == 0 ) {
      out->search_routine = DBSearchImplRoutine_Forward;
    } else {
      warn("String [%s] for search routine is not recognised",temp);
      free_DBSearchImpl(out);
      return NULL;
    }
  }


  return out;
}
  
  







# line 251 "dbsearchimpl.c"
/* Function:  hard_link_DBSearchImpl(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DBSearchImpl *]
 *
 * Return [UNKN ]  Undocumented return value [DBSearchImpl *]
 *
 */
DBSearchImpl * hard_link_DBSearchImpl(DBSearchImpl * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a DBSearchImpl object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  DBSearchImpl_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DBSearchImpl *]
 *
 */
DBSearchImpl * DBSearchImpl_alloc(void) 
{
    DBSearchImpl * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(DBSearchImpl *) ckalloc (sizeof(DBSearchImpl))) == NULL)    {  
      warn("DBSearchImpl_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->type = DBSearchImpl_Serial; 
    out->trace_level = 0;    
    out->suggest_thread_no = (-1);   
    out->search_routine = DBSearchImplRoutine_Exact; 


    return out;  
}    


/* Function:  free_DBSearchImpl(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DBSearchImpl *]
 *
 * Return [UNKN ]  Undocumented return value [DBSearchImpl *]
 *
 */
DBSearchImpl * free_DBSearchImpl(DBSearchImpl * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a DBSearchImpl obj. Should be trappable");  
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
    /* obj->trace_file is linked in */ 


    ckfree(obj); 
    return NULL; 
}    


/* Function:  replace_type_DBSearchImpl(obj,type)
 *
 * Descrip:    Replace member variable type
 *             For use principly by API functions
 *
 *
 * Arg:         obj [UNKN ] Object holding the variable [DBSearchImpl *]
 * Arg:        type [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable type [boolean]
 *
 */
boolean replace_type_DBSearchImpl(DBSearchImpl * obj,int type) 
{
    if( obj == NULL)     {  
      warn("In replacement function type for object DBSearchImpl, got a NULL object");   
      return FALSE;  
      }  
    obj->type = type;    
    return TRUE; 
}    


/* Function:  access_type_DBSearchImpl(obj)
 *
 * Descrip:    Access member variable type
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [DBSearchImpl *]
 *
 * Return [SOFT ]  member variable type [int]
 *
 */
int access_type_DBSearchImpl(DBSearchImpl * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function type for object DBSearchImpl, got a NULL object");  
      return 0;  
      }  
    return obj->type;    
}    


/* Function:  replace_trace_level_DBSearchImpl(obj,trace_level)
 *
 * Descrip:    Replace member variable trace_level
 *             For use principly by API functions
 *
 *
 * Arg:                obj [UNKN ] Object holding the variable [DBSearchImpl *]
 * Arg:        trace_level [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable trace_level [boolean]
 *
 */
boolean replace_trace_level_DBSearchImpl(DBSearchImpl * obj,int trace_level) 
{
    if( obj == NULL)     {  
      warn("In replacement function trace_level for object DBSearchImpl, got a NULL object");    
      return FALSE;  
      }  
    obj->trace_level = trace_level;  
    return TRUE; 
}    


/* Function:  access_trace_level_DBSearchImpl(obj)
 *
 * Descrip:    Access member variable trace_level
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [DBSearchImpl *]
 *
 * Return [SOFT ]  member variable trace_level [int]
 *
 */
int access_trace_level_DBSearchImpl(DBSearchImpl * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function trace_level for object DBSearchImpl, got a NULL object");   
      return 0;  
      }  
    return obj->trace_level;     
}    


/* Function:  replace_trace_file_DBSearchImpl(obj,trace_file)
 *
 * Descrip:    Replace member variable trace_file
 *             For use principly by API functions
 *
 *
 * Arg:               obj [UNKN ] Object holding the variable [DBSearchImpl *]
 * Arg:        trace_file [OWNER] New value of the variable [FILE *]
 *
 * Return [SOFT ]  member variable trace_file [boolean]
 *
 */
boolean replace_trace_file_DBSearchImpl(DBSearchImpl * obj,FILE * trace_file) 
{
    if( obj == NULL)     {  
      warn("In replacement function trace_file for object DBSearchImpl, got a NULL object"); 
      return FALSE;  
      }  
    obj->trace_file = trace_file;    
    return TRUE; 
}    


/* Function:  access_trace_file_DBSearchImpl(obj)
 *
 * Descrip:    Access member variable trace_file
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [DBSearchImpl *]
 *
 * Return [SOFT ]  member variable trace_file [FILE *]
 *
 */
FILE * access_trace_file_DBSearchImpl(DBSearchImpl * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function trace_file for object DBSearchImpl, got a NULL object");    
      return NULL;   
      }  
    return obj->trace_file;  
}    


/* Function:  replace_suggest_thread_no_DBSearchImpl(obj,suggest_thread_no)
 *
 * Descrip:    Replace member variable suggest_thread_no
 *             For use principly by API functions
 *
 *
 * Arg:                      obj [UNKN ] Object holding the variable [DBSearchImpl *]
 * Arg:        suggest_thread_no [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable suggest_thread_no [boolean]
 *
 */
boolean replace_suggest_thread_no_DBSearchImpl(DBSearchImpl * obj,int suggest_thread_no) 
{
    if( obj == NULL)     {  
      warn("In replacement function suggest_thread_no for object DBSearchImpl, got a NULL object");  
      return FALSE;  
      }  
    obj->suggest_thread_no = suggest_thread_no;  
    return TRUE; 
}    


/* Function:  access_suggest_thread_no_DBSearchImpl(obj)
 *
 * Descrip:    Access member variable suggest_thread_no
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [DBSearchImpl *]
 *
 * Return [SOFT ]  member variable suggest_thread_no [int]
 *
 */
int access_suggest_thread_no_DBSearchImpl(DBSearchImpl * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function suggest_thread_no for object DBSearchImpl, got a NULL object"); 
      return 0;  
      }  
    return obj->suggest_thread_no;   
}    


/* Function:  replace_search_routine_DBSearchImpl(obj,search_routine)
 *
 * Descrip:    Replace member variable search_routine
 *             For use principly by API functions
 *
 *
 * Arg:                   obj [UNKN ] Object holding the variable [DBSearchImpl *]
 * Arg:        search_routine [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable search_routine [boolean]
 *
 */
boolean replace_search_routine_DBSearchImpl(DBSearchImpl * obj,int search_routine) 
{
    if( obj == NULL)     {  
      warn("In replacement function search_routine for object DBSearchImpl, got a NULL object"); 
      return FALSE;  
      }  
    obj->search_routine = search_routine;    
    return TRUE; 
}    


/* Function:  access_search_routine_DBSearchImpl(obj)
 *
 * Descrip:    Access member variable search_routine
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [DBSearchImpl *]
 *
 * Return [SOFT ]  member variable search_routine [int]
 *
 */
int access_search_routine_DBSearchImpl(DBSearchImpl * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function search_routine for object DBSearchImpl, got a NULL object");    
      return 0;  
      }  
    return obj->search_routine;  
}    



#ifdef _cplusplus
}
#endif
