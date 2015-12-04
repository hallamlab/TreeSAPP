
#include <stdio.h>

typedef char aa;
typedef int base;
typedef double Probability;
typedef double Bits;
typedef int Score;
typedef int codon;
typedef int boolean;

#define WISE2_FATAL    1
#define WISE2_WARNING  2
#define WISE2_INFO     8
#define WISE2_REPORT   16

/* Function:  Wise2_error_off(type)
 *
 * Descrip:    Really for the API. Wraps some
 *             of the error concepts..
 *
 *
 * Arg:        type [UNKN ] Undocumented argument [int]
 *
 */
void Wise2_error_off(int type);


/* Function:  Wise2_error_on(type)
 *
 * Descrip:    Really for the API. Wraps some
 *             of the error concepts..
 *
 *
 * Arg:        type [UNKN ] Undocumented argument [int]
 *
 */
void Wise2_error_on(int type);

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
void Wise2_warn( char * msg, ... );

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
void Wise2_info( char * msg, ... );

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
void Wise2_fatal( char * msg, ... );


/* These functions are not associated with an object */
/* Function:  Wise2_openfile(filename,passedprot)
 *
 * Descrip:    Every file open goes through this.
 *
 *             It opens for reading in the following order 
 *                .
 *                WISEPERSONALDIR
 *                WISECONFIGDIR
 *
 *             For writing it opens in .
 *
 *             Filenames with ~'s are expanded to HOME/filename
 *
 *
 * Arg:        filename     filename to open for read/writing [const char *]
 * Arg:        passedprot   string representing standard fopen attributes [const char *]
 *
 * Returns open'd filehandle, NULL on error [FILE *]
 *
 */
FILE * Wise2_openfile( const char * filename,const char * passedprot);

