#ifndef DYNAMITEcommandlineHEADERFILE
#define DYNAMITEcommandlineHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "wisebase.h"








    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  strip_out_remaining_options_with_warning(argc,argv)
 *
 * Descrip:    This removes all remaining options, issuing warnings
 *             through warn
 *
 *
 * Arg:        argc [UNKN ] Undocumented argument [int *]
 * Arg:        argv [UNKN ] Undocumented argument [char **]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_strip_out_remaining_options_with_warning(int * argc,char ** argv);
#define strip_out_remaining_options_with_warning Wise2_strip_out_remaining_options_with_warning


/* Function:  strip_out_boolean_argument(argc,argv,tag)
 *
 * Descrip:    This removes argument in tag from the commandline if there and
 *             returns TRUE
 *
 *             otherwise returns FALSE
 *
 *
 * Arg:        argc [UNKN ] argc from main declaration (call as &argc) [int *]
 * Arg:        argv [UNKN ] argv from main declaration [char **]
 * Arg:         tag [READ ] -[string] argument to find [char *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_strip_out_boolean_argument(int * argc,char ** argv,char * tag);
#define strip_out_boolean_argument Wise2_strip_out_boolean_argument


/* Function:  strip_out_boolean_def_argument(argc,argv,tag,value)
 *
 * Descrip:    This sets a boolean argument as
 *                -tag   == TRUE
 *                -notag == FALSE
 *
 *             This is probably better than strip_out_boolean_argument
 *             as 
 *               defaults can be provided
 *               user can switch either on or off
 *
 *
 * Arg:         argc [UNKN ] argc from main declaration (call as &argc) [int *]
 * Arg:         argv [UNKN ] argv from main declaration [char **]
 * Arg:          tag [READ ] -[string] argument to find [char *]
 * Arg:        value [WRITE] boolean pointer to write value to [boolean *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_strip_out_boolean_def_argument(int * argc,char ** argv,char * tag,boolean * value);
#define strip_out_boolean_def_argument Wise2_strip_out_boolean_def_argument


/* Function:  strip_out_assigned_argument(argc,argv,tag)
 *
 * Descrip:    This removes argument in tag from the commandline if there and
 *             returns the argument to it (in -tag arg - ie with a space).
 *
 *             otherwise returns NULL
 *
 *
 * Arg:        argc [UNKN ] argc from main declaration (call as &argc) [int *]
 * Arg:        argv [UNKN ] argv from main declaration [char **]
 * Arg:         tag [READ ] -[string] argument to find [char *]
 *
 * Return [UNKN ]  Undocumented return value [char *]
 *
 */
char * Wise2_strip_out_assigned_argument(int * argc,char ** argv,char * tag);
#define strip_out_assigned_argument Wise2_strip_out_assigned_argument


/* Function:  strip_out_integer_argument(argc,argv,tag,value)
 *
 * Descrip:    This removes argument in tag from the commandline if there, and
 *             looks at the argument as whether it is an int or not.
 *
 *             No tag - returns FALSE and leaves value alone
 *             Tag but no integer value - issues warning, returns FALSE
 *             tag but integer value    - sets value, returns true
 *
 *
 * Arg:         argc [UNKN ] argc from main declaration (call as &argc) [int *]
 * Arg:         argv [UNKN ] argv from main declaration [char **]
 * Arg:          tag [READ ] -[string] argument to find [char *]
 * Arg:        value [WRITE] int pointer to write value to [int *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_strip_out_integer_argument(int * argc,char ** argv,char * tag,int * value);
#define strip_out_integer_argument Wise2_strip_out_integer_argument


/* Function:  strip_out_float_argument(argc,argv,tag,value)
 *
 * Descrip:    This removes argument in tag from the commandline if there, and
 *             looks at the argument as whether it is a double or not.
 *
 *             No tag - returns FALSE and leaves value alone
 *             Tag but no integer value - issues warning, returns FALSE
 *             tag but integer value    - sets value, returns true
 *
 *
 * Arg:         argc [UNKN ] argc from main declaration (call as &argc) [int *]
 * Arg:         argv [UNKN ] argv from main declaration [char **]
 * Arg:          tag [READ ] -[string] argument to find [char *]
 * Arg:        value [WRITE] double pointer to write value to [double *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_strip_out_float_argument(int * argc,char ** argv,char * tag,double * value);
#define strip_out_float_argument Wise2_strip_out_float_argument


/* Function:  show_standard_options(ofp)
 *
 * Descrip:    Shows standard options
 *
 *
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_show_standard_options(FILE * ofp);
#define show_standard_options Wise2_show_standard_options


/* Function:  strip_out_standard_options(argc,show_version,show_help,argv)
 *
 * Descrip:    Handles default arguments correctly, setting help,version,errors
 *             handling
 *
 *
 * Arg:                argc [UNKN ] Undocumented argument [int *]
 * Arg:        show_version [UNKN ] Undocumented argument [NullString]
 * Arg:           show_help [UNKN ] Undocumented argument [NullString]
 * Arg:                argv [UNKN ] Undocumented argument [char **]
 *
 */
void Wise2_strip_out_standard_options(int * argc,char ** argv,void (*show_help)(FILE * ofp),void (*show_version)(FILE * ofp));
#define strip_out_standard_options Wise2_strip_out_standard_options


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/

#ifdef _cplusplus
}
#endif

#endif
