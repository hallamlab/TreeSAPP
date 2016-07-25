#ifndef DYNAMITEwisetimeHEADERFILE
#define DYNAMITEwisetimeHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "wisebase.h"
#include <time.h>


    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  now_string(void)
 *
 * Descrip:    returns a buffer
 *             of the current date and time
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [char *]
 *
 */
char * Wise2_now_string(void);
#define now_string Wise2_now_string


/* Function:  time_stamp(ofp)
 *
 * Descrip:    puts now_string into file, with
 *             no newline
 *
 *
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_time_stamp(FILE * ofp);
#define time_stamp Wise2_time_stamp


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
