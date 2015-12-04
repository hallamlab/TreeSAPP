#ifndef DYNAMITEdyshatterHEADERFILE
#define DYNAMITEdyshatterHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dynafunc.h"





    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  write_shatter_functions(dfp,gm,dpi)
 *
 * Descrip:    Writes shatter functions
 *
 *
 * Arg:        dfp [UNKN ] Undocumented argument [DYNFILE *]
 * Arg:         gm [UNKN ] Undocumented argument [GenericMatrix *]
 * Arg:        dpi [UNKN ] Undocumented argument [DPImplementation *]
 *
 */
void write_shatter_functions(DYNFILE * dfp,GenericMatrix * gm,DPImplementation * dpi);


/* Function:  matrix_calculate_shatter_func(dfp,gm)
 *
 * Descrip:    for shatter
 *
 *
 * Arg:        dfp [UNKN ] Undocumented argument [DYNFILE *]
 * Arg:         gm [UNKN ] Undocumented argument [GenericMatrix *]
 *
 */
void matrix_calculate_shatter_func(DYNFILE * dfp,GenericMatrix * gm);


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
void write_shatter_access_funcs(DYNFILE * dfp,GenericMatrix * gm);
void write_shatter_read_func(DYNFILE * dfp,GenericMatrix * gm);
void write_score_block_shatter(DYNFILE * dfp,GenericMatrix * gm,char * matrixtag,char * pointertag,char * specialtag,boolean use_special);


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/

#ifdef _cplusplus
}
#endif

#endif
