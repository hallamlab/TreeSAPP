#ifdef _cplusplus
extern "C" {
#endif
#include "wisetime.h"

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
# line 14 "wisetime.dy"
char * now_string(void)
{
  char * temp;
  char * runner;
  time_t now;

  now=time(NULL);
  temp=ctime(&now);
  runner=temp+strlen(temp);
  runner--;
  *runner='\0'; /* got rid of annoying new line */
  return temp;
}

/* Function:  time_stamp(ofp)
 *
 * Descrip:    puts now_string into file, with
 *             no newline
 *
 *
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 32 "wisetime.dy"
void time_stamp(FILE * ofp)
{
  time_t now;

  now=time(NULL);
  fputs(ctime(&now),ofp);
  
  return;
}

# line 46 "wisetime.c"

#ifdef _cplusplus
}
#endif
