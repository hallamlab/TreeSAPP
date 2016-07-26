

/* Helper functions in the module
 *
 * Wise2_strip_out_remaining_options_with_warning
 * Wise2_strip_out_boolean_argument
 * Wise2_strip_out_assigned_argument
 * Wise2_strip_out_integer_argument
 * Wise2_strip_out_float_argument
 *



/* These functions are not associated with an object */
/* Function:  Wise2_strip_out_remaining_options_with_warning(argc,argv)
 *
 * Descrip:    This removes all remaining options, issuing warnings
 *             through warn
 *
 *
 * Arg:        argc         Undocumented argument [int *]
 * Arg:        argv         Undocumented argument [char **]
 *
 * Returns Undocumented return value [boolean]
 *
 */
boolean Wise2_strip_out_remaining_options_with_warning( int * argc,char ** argv);

/* Function:  Wise2_strip_out_boolean_argument(argc,argv,tag)
 *
 * Descrip:    This removes argument in tag from the commandline if there and
 *             returns TRUE
 *
 *             otherwise returns FALSE
 *
 *
 * Arg:        argc         argc from main declaration (call as &argc) [int *]
 * Arg:        argv         argv from main declaration [char **]
 * Arg:        tag          -[string] argument to find [char *]
 *
 * Returns Undocumented return value [boolean]
 *
 */
boolean Wise2_strip_out_boolean_argument( int * argc,char ** argv,char * tag);

/* Function:  Wise2_strip_out_assigned_argument(argc,argv,tag)
 *
 * Descrip:    This removes argument in tag from the commandline if there and
 *             returns the argument to it (in -tag arg - ie with a space).
 *
 *             otherwise returns NULL
 *
 *
 * Arg:        argc         argc from main declaration (call as &argc) [int *]
 * Arg:        argv         argv from main declaration [char **]
 * Arg:        tag          -[string] argument to find [char *]
 *
 * Returns Undocumented return value [char *]
 *
 */
char * Wise2_strip_out_assigned_argument( int * argc,char ** argv,char * tag);

/* Function:  Wise2_strip_out_integer_argument(argc,argv,tag,value)
 *
 * Descrip:    This removes argument in tag from the commandline if there, and
 *             looks at the argument as whether it is an int or not.
 *
 *             No tag - returns FALSE and leaves value alone
 *             Tag but no integer value - issues warning, returns FALSE
 *             tag but integer value    - sets value, returns true
 *
 *
 * Arg:        argc         argc from main declaration (call as &argc) [int *]
 * Arg:        argv         argv from main declaration [char **]
 * Arg:        tag          -[string] argument to find [char *]
 * Arg:        value        int pointer to write value to [int *]
 *
 * Returns Undocumented return value [boolean]
 *
 */
boolean Wise2_strip_out_integer_argument( int * argc,char ** argv,char * tag,int * value);

/* Function:  Wise2_strip_out_float_argument(argc,argv,tag,value)
 *
 * Descrip:    This removes argument in tag from the commandline if there, and
 *             looks at the argument as whether it is a double or not.
 *
 *             No tag - returns FALSE and leaves value alone
 *             Tag but no integer value - issues warning, returns FALSE
 *             tag but integer value    - sets value, returns true
 *
 *
 * Arg:        argc         argc from main declaration (call as &argc) [int *]
 * Arg:        argv         argv from main declaration [char **]
 * Arg:        tag          -[string] argument to find [char *]
 * Arg:        value        double pointer to write value to [double *]
 *
 * Returns Undocumented return value [boolean]
 *
 */
boolean Wise2_strip_out_float_argument( int * argc,char ** argv,char * tag,double * value);

