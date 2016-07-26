#ifndef DYNAMITEdynadbHEADERFILE
#define DYNAMITEdynadbHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna2.h"
#include "dynafunc.h"
#include "dpimpl.h"



    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  make_search_loop_function(dfp,gm)
 *
 * Descrip:    Makes the serial search function, which loops
 *             over databases 
 *
 *
 * Arg:        dfp [UNKN ] Undocumented argument [DYNFILE *]
 * Arg:         gm [UNKN ] Undocumented argument [GenericMatrix *]
 *
 */
void make_search_loop_function(DYNFILE * dfp,GenericMatrix * gm);


/* Function:  write_one_score_GenericMatrix(dfp,gm,dpi)
 *
 * Descrip:    Makes the score only function, which gives the score
 *             for two objects. Used in the serial and the pthreads
 *             ports
 *
 *
 * Arg:        dfp [UNKN ] Undocumented argument [DYNFILE *]
 * Arg:         gm [UNKN ] Undocumented argument [GenericMatrix *]
 * Arg:        dpi [UNKN ] Undocumented argument [DPImplementation *]
 *
 */
void write_one_score_GenericMatrix(DYNFILE * dfp,GenericMatrix * gm,DPImplementation * dpi);


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
