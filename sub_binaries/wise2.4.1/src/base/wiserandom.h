#ifndef DYNAMITEwiserandomHEADERFILE
#define DYNAMITEwiserandomHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "wisebase.h"
#include <time.h>
#include <limits.h>




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  init_random(void)
 *
 * Descrip:    initates the random generator to
 *             time byte...
 *
 *
 *
 */
void Wise2_init_random(void);
#define init_random Wise2_init_random


/* Function:  random_integer(l)
 *
 * Descrip:    returns an integer between 0 and l
 *             though I don't think we will get 0
 *             very often. Hmmm.
 *
 *
 * Arg:        l [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_random_integer(int l);
#define random_integer Wise2_random_integer


/* Function:  random_0_to_1(void)
 *
 * Descrip:    returns a random number between
 *             0 and 1
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [double]
 *
 */
double Wise2_random_0_to_1(void);
#define random_0_to_1 Wise2_random_0_to_1


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
