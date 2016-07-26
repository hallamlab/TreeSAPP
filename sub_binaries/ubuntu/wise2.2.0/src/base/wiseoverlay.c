#ifdef _cplusplus
extern "C" {
#endif
#include "wiseoverlay.h"

static boolean isinuse=FALSE;
static FILE * overlay = NULL;
static int delete_len=0;


/* Function:  start_overlay(over)
 *
 * Descrip:    actually starts on overlay system
 *             on the FILE*. Should really by stderr...
 *             nothing else makes sense.
 *
 *
 * Arg:        over [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 22 "wiseoverlay.dy"
void start_overlay(FILE * over)
{
  overlay=over;
  isinuse=TRUE;
}

/* Function:  stop_overlay(void)
 *
 * Descrip:    finishes an overlay system by putting in a
 *             newline and clearing the static variables
 *
 *
 *
 */
# line 32 "wiseoverlay.dy"
void stop_overlay(void)
{
  fprintf(overlay,"\n");
  delete_len=0;
  overlay=NULL;
  isinuse=FALSE;
}

/* Function:  print_overlay(msg,)
 *
 * Descrip:    Does the business. Deletes the previous
 *             message with \b and then prints the current
 *             one ontop. 
 *
 *             don't put in \n's otherwise this is going to 
 *             look yukky ;)
 *
 *
 * Arg:        msg [UNKN ] Undocumented argument [char *]
 * Arg:            [UNKN ] Undocumented argument [.]
 *
 */
# line 48 "wiseoverlay.dy"
void print_overlay(char * msg, ... )
{
  char buffer[512];
  va_list ap;
  register int i;	
  
  if( isinuse == FALSE)
    /* error message */
    return;
  
  
  va_start(ap,msg);
	
  for(i=0;i<delete_len;i++)
    fputc('\b',overlay);
  
  vsprintf(buffer,msg,ap);
  
  delete_len=strlen(buffer);
  
  fprintf(overlay,"%s",buffer);
  fflush(overlay);
  
  va_end(ap);
  
  return;
}

# line 81 "wiseoverlay.c"

#ifdef _cplusplus
}
#endif
