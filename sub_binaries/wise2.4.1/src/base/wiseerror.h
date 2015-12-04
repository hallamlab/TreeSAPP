#ifndef DYNAMITEwiseerrorHEADERFILE
#define DYNAMITEwiseerrorHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif

#include "wisebase.h"



typedef int Flag;

#define MAXERROR 256
#define MAXERRORCALL 32
#define MAXMSGSTACKERROR 64

/* flags on errors */
#define ERRORUSE      1
#define ERRORTOSTDERR 2
#define ERRORTOLOG    4
#define ERRORTOCALL   8
#define LONGERROR     16 /* not used */


/* types of error */
#define FATAL    1
#define WARNING  2
#define PEDANTIC 4 /* deprecated */
#define INFO     8
#define REPORT   16

#define erroroff(type)       error_flag_off(type,ERRORUSE)
#define erroron(type)        error_flag_on (type,ERRORUSE)
#define errorcallon(type)    error_flag_on (type,ERRORTOCALL)
#define errorcalloff(type)   error_flag_off(type,ERRORTOCALL)
#define errorstderron(type)  error_flag_on (type,ERRORTOSTDERR)
#define errorstderroff(type) error_flag_off(type,ERRORTOSTDERR)
#define errorlogon(type)     error_flag_on (type,ERRORTOLOG)
#define errorlogoff(type)    error_flag_off(type,ERRORTOLOG)

typedef enum Error_Display_Type {
  ERROR_DISPLAY_PROGRAM = 36,
  ERROR_DISPLAY_SERVER
} ErrorDisplayType;



    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  set_error_display_style(t)
 *
 * Descrip:    sets error display style
 *
 *
 * Arg:        t [UNKN ] Undocumented argument [ErrorDisplayType]
 *
 */
void Wise2_set_error_display_style(ErrorDisplayType t);
#define set_error_display_style Wise2_set_error_display_style


/* Function:  set_log_display_string(str)
 *
 * Descrip:    sets string for log display type 
 *
 *
 * Arg:        str [UNKN ] Undocumented argument [char *]
 *
 */
void Wise2_set_log_display_string(char * str);
#define set_log_display_string Wise2_set_log_display_string


/* Function:  push_errormsg_stack_call(ecall)
 *
 * Descrip:    This adds a function call for people who want
 *             to register error handling functions
 *
 *             Probably best wrapped by a separate function
 *
 *
 * Arg:        ecall [UNKN ] Undocumented argument [NullString]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_push_errormsg_stack_call( char * (*ecall)(void));
#define push_errormsg_stack_call Wise2_push_errormsg_stack_call


/* Function:  pop_errormsg_stack(void)
 *
 * Descrip:    This removes a error message from the stack
 *
 *
 *
 */
void Wise2_pop_errormsg_stack(void);
#define pop_errormsg_stack Wise2_pop_errormsg_stack


/* Function:  add_log_filename(filename)
 *
 * Descrip:    Opens filename as the log file. 
 *
 *
 * Arg:        filename [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_log_filename(char * filename);
#define add_log_filename Wise2_add_log_filename


/* Function:  error_off(type)
 *
 * Descrip:    Really for the API. Wraps some
 *             of the error concepts..
 *
 *
 * Arg:        type [UNKN ] Undocumented argument [int]
 *
 */
void Wise2_error_off(int type);
#define error_off Wise2_error_off


/* Function:  error_on(type)
 *
 * Descrip:    Really for the API. Wraps some
 *             of the error concepts..
 *
 *
 * Arg:        type [UNKN ] Undocumented argument [int]
 *
 */
void Wise2_error_on(int type);
#define error_on Wise2_error_on


/* Function:  error_flag_on(type,f)
 *
 * Descrip:    Turns on the particular type of error flag
 *             (eg, STDERR etc). 
 *
 *
 * Arg:        type [UNKN ] Undocumented argument [int]
 * Arg:           f [UNKN ] Undocumented argument [Flag]
 *
 */
void Wise2_error_flag_on(int type,Flag f);
#define error_flag_on Wise2_error_flag_on


/* Function:  error_flag_off(type,f)
 *
 * Descrip:    Turns off the particular type of error flag
 *             (eg, STDERR etc). 
 *
 *
 * Arg:        type [UNKN ] Undocumented argument [int]
 * Arg:           f [UNKN ] Undocumented argument [Flag]
 *
 */
void Wise2_error_flag_off(int type,Flag f);
#define error_flag_off Wise2_error_flag_off


/* Function:  catch_errors(catch)
 *
 * Descrip:    This is a wrapper for the error handling
 *             system. It does the following things
 *
 *             Sets function as the function to process errors
 *
 *             Switches the INFO,ERROR and FATAL flags off on STDERR
 *             and on to ERROR CALL.
 *
 *
 * Arg:        catch [FUNCP]  [void (*catch]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_catch_errors(void (*catch)(char *,int));
#define catch_errors Wise2_catch_errors


/* Function:  stop_catching_errors(void)
 *
 * Descrip:    Switches off error catching,
 *             putting flags back on to STDERR
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_stop_catching_errors(void);
#define stop_catching_errors Wise2_stop_catching_errors


/* Function:  push_error_call()
 *
 * Descrip:    Registers this function for dealing with errors
 *
 *             Try to use catch_errors instead
 *
 *
 *
 * Arg:         [UNKN ] Undocumented argument [NullString]
 *
 */
void Wise2_push_error_call(void (* func)(char *,int));
#define push_error_call Wise2_push_error_call


/* Function:  pop_error_call(void)
 *
 * Descrip:    Discards current function for dealing with errors
 *
 *
 *
 */
void Wise2_pop_error_call(void);
#define pop_error_call Wise2_pop_error_call


/* Function:  type_to_error(type)
 *
 * Descrip:    Turns int error types to Names
 *             for display purposes.
 *
 *
 * Arg:        type [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [char *]
 *
 */
char * Wise2_type_to_error(int type);
#define type_to_error Wise2_type_to_error


/* Function:  info(msg,)
 *
 * Descrip:    Produces a 'info' error message.
 *
 *
 * Arg:        msg [UNKN ] Undocumented argument [char *]
 * Arg:            [UNKN ] Undocumented argument [.]
 *
 */
void Wise2_info(char * msg, ...);
#define info Wise2_info


/* Function:  debug_report(msg,)
 *
 * Descrip:    Produces a 'debug report' error message.
 *
 *
 * Arg:        msg [UNKN ] Undocumented argument [char *]
 * Arg:            [UNKN ] Undocumented argument [.]
 *
 */
void Wise2_debug_report(char * msg, ...);
#define debug_report Wise2_debug_report


/* Function:  warn(msg,)
 *
 * Descrip:    Produces a 'warn' error message.
 *
 *
 * Arg:        msg [UNKN ] Undocumented argument [char *]
 * Arg:            [UNKN ] Undocumented argument [.]
 *
 */
void Wise2_warn(char * msg, ...);
#define warn Wise2_warn


/* Function:  fatal(msg,)
 *
 * Descrip:    Produces a 'fatal' error message.
 *
 *
 * Arg:        msg [UNKN ] Undocumented argument [char *]
 * Arg:            [UNKN ] Undocumented argument [.]
 *
 */
void Wise2_fatal(char * msg, ...);
#define fatal Wise2_fatal


/* Function:  log_full_error(type,location,msg,)
 *
 * Descrip:    Deprecated
 *
 *             produces any of the error types
 *
 *
 * Arg:            type [UNKN ] Undocumented argument [int]
 * Arg:        location [UNKN ] Undocumented argument [int]
 * Arg:             msg [UNKN ] Undocumented argument [char *]
 * Arg:                 [UNKN ] Undocumented argument [.]
 *
 */
void Wise2_log_full_error(int type,int location,char * msg, ...);
#define log_full_error Wise2_log_full_error


/* Function:  start_reporting(msg,)
 *
 * Descrip:    Starts a % reporting run. This is the header message
 *
 *
 * Arg:        msg [UNKN ] Undocumented argument [char *]
 * Arg:            [UNKN ] Undocumented argument [.]
 *
 */
void Wise2_start_reporting(char * msg,...);
#define start_reporting Wise2_start_reporting


/* Function:  stop_reporting(void)
 *
 * Descrip:    Stops a % reporting run. 
 *
 *
 *
 */
void Wise2_stop_reporting(void);
#define stop_reporting Wise2_stop_reporting


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
boolean Wise2_push_errormsg_stack(char * msg, ...);
#define push_errormsg_stack Wise2_push_errormsg_stack
void Wise2_show_message_stack(FILE * ofp);
#define show_message_stack Wise2_show_message_stack
void Wise2_add_log_file(FILE * ofp);
#define add_log_file Wise2_add_log_file
void Wise2_show_error(Flag flag,char * othermsg,int type);
#define show_error Wise2_show_error

#ifdef _cplusplus
}
#endif

#endif
