

/* Helper functions in the module
 *
 * Wise2_catch_errors
 * Wise2_stop_catching_errors
 * Wise2_warn
 * Wise2_info
 * Wise2_fatal
 *



/* These functions are not associated with an object */
/* Function:  Wise2_catch_errors(catch)
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
 * Arg:        catch         [void (*catch]
 *
 * Returns Undocumented return value [boolean]
 *
 */
boolean Wise2_catch_errors( void (*catch)(char *,int));

/* Function:  Wise2_stop_catching_errors(void)
 *
 * Descrip:    Switches off error catching,
 *             putting flags back on to STDERR
 *
 *
 *
 * Returns Undocumented return value [boolean]
 *
 */
boolean Wise2_stop_catching_errors();

/* Function:  Wise2_warn(msg,)
 *
 * Descrip:    Produces a 'warn' error message.
 *
 *
 * Arg:        msg          Undocumented argument [char *]
 * Arg:                     Undocumented argument [Wise2_.]
 *
 * Returns Undocumented return value [void]
 *
 */
void Wise2_warn( char * msg,Wise2_. );

/* Function:  Wise2_info(msg,)
 *
 * Descrip:    Produces a 'info' error message.
 *
 *
 * Arg:        msg          Undocumented argument [char *]
 * Arg:                     Undocumented argument [Wise2_.]
 *
 * Returns Undocumented return value [void]
 *
 */
void Wise2_info( char * msg,Wise2_. );

/* Function:  Wise2_fatal(msg,)
 *
 * Descrip:    Produces a 'fatal' error message.
 *
 *
 * Arg:        msg          Undocumented argument [char *]
 * Arg:                     Undocumented argument [Wise2_.]
 *
 * Returns Undocumented return value [void]
 *
 */
void Wise2_fatal( char * msg,Wise2_. );

