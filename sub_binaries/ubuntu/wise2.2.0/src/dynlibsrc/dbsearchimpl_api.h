

/* Functions that create, manipulate or act on DBSearchImpl
 *
 * Wise2_impl_string_DBSearchImpl
 * Wise2_hard_link_DBSearchImpl
 * Wise2_DBSearchImpl_alloc
 * Wise2_replace_type_DBSearchImpl
 * Wise2_access_type_DBSearchImpl
 * Wise2_replace_trace_level_DBSearchImpl
 * Wise2_access_trace_level_DBSearchImpl
 * Wise2_replace_trace_file_DBSearchImpl
 * Wise2_access_trace_file_DBSearchImpl
 * Wise2_replace_suggest_thread_no_DBSearchImpl
 * Wise2_access_suggest_thread_no_DBSearchImpl
 * Wise2_replace_search_routine_DBSearchImpl
 * Wise2_access_search_routine_DBSearchImpl
 * Wise2_free_DBSearchImpl [destructor]
 *
 */



/* Helper functions in the module
 *
 * Wise2_new_pthread_DBSearchImpl
 * Wise2_new_serial_DBSearchImpl
 *

/* API for object DBSearchImpl */
/* Function:  Wise2_impl_string_DBSearchImpl(dbsi)
 *
 * Descrip:    Gets a static text string out of the
 *             search implementation 
 *
 *
 * Arg:        dbsi         Undocumented argument [Wise2_DBSearchImpl *]
 *
 * Returns string of the search implementation [char *]
 *
 */
char * Wise2_impl_string_DBSearchImpl( Wise2_DBSearchImpl * dbsi);

/* Function:  Wise2_hard_link_DBSearchImpl(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_DBSearchImpl *]
 *
 * Returns Undocumented return value [Wise2_DBSearchImpl *]
 *
 */
Wise2_DBSearchImpl * Wise2_hard_link_DBSearchImpl( Wise2_DBSearchImpl * obj);

/* Function:  Wise2_DBSearchImpl_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Returns Undocumented return value [Wise2_DBSearchImpl *]
 *
 */
Wise2_DBSearchImpl * Wise2_DBSearchImpl_alloc();

/* Function:  Wise2_replace_type_DBSearchImpl(obj,type)
 *
 * Descrip:    Replace member variable type
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_DBSearchImpl *]
 * Arg:        type         New value of the variable [int]
 *
 * Returns member variable type [boolean]
 *
 */
boolean Wise2_replace_type_DBSearchImpl( Wise2_DBSearchImpl * obj,int type);

/* Function:  Wise2_access_type_DBSearchImpl(obj)
 *
 * Descrip:    Access member variable type
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_DBSearchImpl *]
 *
 * Returns member variable type [int]
 *
 */
int Wise2_access_type_DBSearchImpl( Wise2_DBSearchImpl * obj);

/* Function:  Wise2_replace_trace_level_DBSearchImpl(obj,trace_level)
 *
 * Descrip:    Replace member variable trace_level
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_DBSearchImpl *]
 * Arg:        trace_level  New value of the variable [int]
 *
 * Returns member variable trace_level [boolean]
 *
 */
boolean Wise2_replace_trace_level_DBSearchImpl( Wise2_DBSearchImpl * obj,int trace_level);

/* Function:  Wise2_access_trace_level_DBSearchImpl(obj)
 *
 * Descrip:    Access member variable trace_level
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_DBSearchImpl *]
 *
 * Returns member variable trace_level [int]
 *
 */
int Wise2_access_trace_level_DBSearchImpl( Wise2_DBSearchImpl * obj);

/* Function:  Wise2_replace_trace_file_DBSearchImpl(obj,trace_file)
 *
 * Descrip:    Replace member variable trace_file
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_DBSearchImpl *]
 * Arg:        trace_file   New value of the variable [FILE *]
 *
 * Returns member variable trace_file [boolean]
 *
 */
boolean Wise2_replace_trace_file_DBSearchImpl( Wise2_DBSearchImpl * obj,FILE * trace_file);

/* Function:  Wise2_access_trace_file_DBSearchImpl(obj)
 *
 * Descrip:    Access member variable trace_file
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_DBSearchImpl *]
 *
 * Returns member variable trace_file [FILE *]
 *
 */
FILE * Wise2_access_trace_file_DBSearchImpl( Wise2_DBSearchImpl * obj);

/* Function:  Wise2_replace_suggest_thread_no_DBSearchImpl(obj,suggest_thread_no)
 *
 * Descrip:    Replace member variable suggest_thread_no
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_DBSearchImpl *]
 * Arg:        suggest_thread_no New value of the variable [int]
 *
 * Returns member variable suggest_thread_no [boolean]
 *
 */
boolean Wise2_replace_suggest_thread_no_DBSearchImpl( Wise2_DBSearchImpl * obj,int suggest_thread_no);

/* Function:  Wise2_access_suggest_thread_no_DBSearchImpl(obj)
 *
 * Descrip:    Access member variable suggest_thread_no
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_DBSearchImpl *]
 *
 * Returns member variable suggest_thread_no [int]
 *
 */
int Wise2_access_suggest_thread_no_DBSearchImpl( Wise2_DBSearchImpl * obj);

/* Function:  Wise2_replace_search_routine_DBSearchImpl(obj,search_routine)
 *
 * Descrip:    Replace member variable search_routine
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_DBSearchImpl *]
 * Arg:        search_routine New value of the variable [int]
 *
 * Returns member variable search_routine [boolean]
 *
 */
boolean Wise2_replace_search_routine_DBSearchImpl( Wise2_DBSearchImpl * obj,int search_routine);

/* Function:  Wise2_access_search_routine_DBSearchImpl(obj)
 *
 * Descrip:    Access member variable search_routine
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_DBSearchImpl *]
 *
 * Returns member variable search_routine [int]
 *
 */
int Wise2_access_search_routine_DBSearchImpl( Wise2_DBSearchImpl * obj);

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_DBSearchImpl(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_DBSearchImpl *]
 *
 * Returns Undocumented return value [Wise2_DBSearchImpl *]
 *
 */
Wise2_DBSearchImpl * Wise2_free_DBSearchImpl( Wise2_DBSearchImpl * obj);



/* These functions are not associated with an object */
/* Function:  Wise2_new_pthread_DBSearchImpl(void)
 *
 * Descrip:    Makes a new pthreaded DBSearchImpl
 *
 *             For use mainly for api's who don't want
 *             to initalize the object from the command
 *             line
 *
 *
 *
 * Returns Undocumented return value [Wise2_DBSearchImpl *]
 *
 */
Wise2_DBSearchImpl * Wise2_new_pthread_DBSearchImpl();

/* Function:  Wise2_new_serial_DBSearchImpl(void)
 *
 * Descrip:    Makes a new serial DBSearchImpl
 *
 *             For use mainly for api's who don't want
 *             to initalize the object from the command
 *             line
 *
 *
 *
 * Returns Undocumented return value [Wise2_DBSearchImpl *]
 *
 */
Wise2_DBSearchImpl * Wise2_new_serial_DBSearchImpl();

