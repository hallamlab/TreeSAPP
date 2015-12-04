#ifndef DYNAMITEdbsearchimplHEADERFILE
#define DYNAMITEdbsearchimplHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "wisebase.h"

enum DBSearchImpl_Type {
  DBSearchImpl_Serial,
  DBSearchImpl_Pthreads,
  DBSearchImpl_MPI,
  DBSearchImpl_PVM
 };

/* lifted from HMMer2 for this sort of work. So thanks to
 * Sean Eddy and seems like warren gish as well 
 */

/* Our problem here is that POSIX apparently doesn't specify
 * a standard for how to get sysconf() to report the number of
 * processors on-line. _SC_NPROCESSORS_ONLN is specified
 * by SVR4.0MP. Thanks to W. Gish for help here.
 */
#ifdef  _SC_NPROCESSORS_ONLN    /* Sun Solaris, Digital UNIX */
#define DBSEARCH_NPROC  sysconf(_SC_NPROCESSORS_ONLN)
#else
#ifdef _SC_NPROC_ONLN		/* Silicon Graphics IRIX */
#define DBSEARCH_NPROC  sysconf(_SC_NPROC_ONLN)
#else   /* FreeBSD, Linux don't support getting ncpu via sysconf() */
#define DBSEARCH_NPROC  -1
#endif
#endif

enum DBSearchImpl_Routine {
  DBSearchImplRoutine_Exact = 0,
  DBSearchImplRoutine_Kbest,
  DBSearchImplRoutine_Forward};

/* Object DBSearchImpl
 *
 * Descrip: No Description
 *
 */
struct Wise2_DBSearchImpl {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int type;    
    int trace_level;    /*  how much debugging information to print */ 
    FILE * trace_file;  /*  for writing out trace of db stuff */ 
    int suggest_thread_no;  /*  default, -1, means the use a call to _SC_NPROC */ 
    int search_routine; /*  routine used for the calculation, exact/kbest */ 
    } ;  
/* DBSearchImpl defined */ 
#ifndef DYNAMITE_DEFINED_DBSearchImpl
typedef struct Wise2_DBSearchImpl Wise2_DBSearchImpl;
#define DBSearchImpl Wise2_DBSearchImpl
#define DYNAMITE_DEFINED_DBSearchImpl
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



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
char * Wise2_impl_string_DBSearchImpl(DBSearchImpl * dbsi);
#define impl_string_DBSearchImpl Wise2_impl_string_DBSearchImpl


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
char * Wise2_impl_string_routine_DBSearchImpl(DBSearchImpl * dbsi);
#define impl_string_routine_DBSearchImpl Wise2_impl_string_routine_DBSearchImpl


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
int Wise2_number_of_threads_DBSearchImpl(DBSearchImpl * dbsi) ;
#define number_of_threads_DBSearchImpl Wise2_number_of_threads_DBSearchImpl


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
DBSearchImpl * Wise2_new_MPI_DBSearchImpl(void);
#define new_MPI_DBSearchImpl Wise2_new_MPI_DBSearchImpl


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
DBSearchImpl * Wise2_new_pthread_DBSearchImpl(void);
#define new_pthread_DBSearchImpl Wise2_new_pthread_DBSearchImpl


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
DBSearchImpl * Wise2_new_serial_DBSearchImpl(void);
#define new_serial_DBSearchImpl Wise2_new_serial_DBSearchImpl


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
void Wise2_show_help_DBSearchImpl(FILE * ofp);
#define show_help_DBSearchImpl Wise2_show_help_DBSearchImpl


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
DBSearchImpl * Wise2_new_DBSearchImpl_from_argv(int * argc,char ** argv);
#define new_DBSearchImpl_from_argv Wise2_new_DBSearchImpl_from_argv


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
DBSearchImpl * Wise2_hard_link_DBSearchImpl(DBSearchImpl * obj);
#define hard_link_DBSearchImpl Wise2_hard_link_DBSearchImpl


/* Function:  DBSearchImpl_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DBSearchImpl *]
 *
 */
DBSearchImpl * Wise2_DBSearchImpl_alloc(void);
#define DBSearchImpl_alloc Wise2_DBSearchImpl_alloc


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
DBSearchImpl * Wise2_free_DBSearchImpl(DBSearchImpl * obj);
#define free_DBSearchImpl Wise2_free_DBSearchImpl


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
boolean Wise2_replace_trace_file_DBSearchImpl(DBSearchImpl * obj,FILE * trace_file);
#define replace_trace_file_DBSearchImpl Wise2_replace_trace_file_DBSearchImpl
FILE * Wise2_access_trace_file_DBSearchImpl(DBSearchImpl * obj);
#define access_trace_file_DBSearchImpl Wise2_access_trace_file_DBSearchImpl
int Wise2_access_type_DBSearchImpl(DBSearchImpl * obj);
#define access_type_DBSearchImpl Wise2_access_type_DBSearchImpl
boolean Wise2_replace_suggest_thread_no_DBSearchImpl(DBSearchImpl * obj,int suggest_thread_no);
#define replace_suggest_thread_no_DBSearchImpl Wise2_replace_suggest_thread_no_DBSearchImpl
int Wise2_access_trace_level_DBSearchImpl(DBSearchImpl * obj);
#define access_trace_level_DBSearchImpl Wise2_access_trace_level_DBSearchImpl
int Wise2_access_suggest_thread_no_DBSearchImpl(DBSearchImpl * obj);
#define access_suggest_thread_no_DBSearchImpl Wise2_access_suggest_thread_no_DBSearchImpl
boolean Wise2_replace_trace_level_DBSearchImpl(DBSearchImpl * obj,int trace_level);
#define replace_trace_level_DBSearchImpl Wise2_replace_trace_level_DBSearchImpl
boolean Wise2_replace_search_routine_DBSearchImpl(DBSearchImpl * obj,int search_routine);
#define replace_search_routine_DBSearchImpl Wise2_replace_search_routine_DBSearchImpl
boolean Wise2_replace_type_DBSearchImpl(DBSearchImpl * obj,int type);
#define replace_type_DBSearchImpl Wise2_replace_type_DBSearchImpl
int Wise2_access_search_routine_DBSearchImpl(DBSearchImpl * obj);
#define access_search_routine_DBSearchImpl Wise2_access_search_routine_DBSearchImpl

#ifdef _cplusplus
}
#endif

#endif
