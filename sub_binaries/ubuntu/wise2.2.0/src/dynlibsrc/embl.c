#ifdef _cplusplus
extern "C" {
#endif
#include "embl.h"


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
# line 20 "embl.dy"
boolean next_feature_tab_line(char * buffer,int maxlen,FILE * ifp)
{
  char * runner;

  if( fgets(buffer,maxlen,ifp) == NULL ) {
    buffer[0]= '\0';
    return FALSE;
  }

  if( strstartcmp(buffer,"FT") != 0) {
    return FALSE;
  }

  for(runner = buffer+2;*runner && isspace((int)*runner);runner++)
    ;

  if( *runner == '\0' ) {
    warn("A supposed EMBL feature line with nothing on it");
    return FALSE;
  }

  if( runner - buffer > 3 ) {
    memmove(buffer,runner,strlen(runner));
    return TRUE;
  }

  return FALSE;
}

# line 50 "embl.c"

#ifdef _cplusplus
}
#endif
