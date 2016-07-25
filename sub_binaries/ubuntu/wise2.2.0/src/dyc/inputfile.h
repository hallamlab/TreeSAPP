#ifndef DYNAMITEinputfileHEADERFILE
#define DYNAMITEinputfileHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "wisebase.h"




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
boolean open_watch_file(FILE * ifp);
boolean close_watch_file(void);
char * get_watched_line(char * buffer,int max,FILE * ifp);
int get_watched_linecount(void);


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/

#ifdef _cplusplus
}
#endif

#endif
