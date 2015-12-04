#ifndef DYNAMITEemblHEADERFILE
#define DYNAMITEemblHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "wisebase.h"



    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  next_feature_tab_line(buffer,maxlen,ifp)
 *
 * Descrip:    internal function for reading EMBL feature
 *             tables.
 *
 *             reads in line. If id 3 in, returns FALSE. 
 *             Otherwise strips FT up to central tab region. 
 *
 *
 * Arg:        buffer [UNKN ] Undocumented argument [char *]
 * Arg:        maxlen [UNKN ] Undocumented argument [int]
 * Arg:           ifp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_next_feature_tab_line(char * buffer,int maxlen,FILE * ifp);
#define next_feature_tab_line Wise2_next_feature_tab_line


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
