#ifdef _cplusplus
extern "C" {
#endif
#include "inputfile.h"


static FILE * watch_file = NULL;
static linecount = 0;


# line 18 "inputfile.dy"
boolean open_watch_file(FILE * ifp)
{
  if( watch_file != NULL ) {
    warn("Can't watch a new file, already set!");
    return FALSE;
  }

  watch_file = ifp;
  linecount = 0;
  return TRUE;
}

# line 30 "inputfile.dy"
boolean close_watch_file(void)
{
  if( watch_file == NULL ) {
    warn("No watch file to close");
    return FALSE;
  }

  watch_file = NULL;
  return TRUE;
}

# line 41 "inputfile.dy"
char * get_watched_line(char * buffer,int max,FILE * ifp)
{
  if( ifp != watch_file) {
    warn("Trying to get a watched line without watching the file!");
    return fgets(buffer,max,ifp);
  }

  linecount++;
  return fgets(buffer,max,ifp);
}

# line 52 "inputfile.dy"
int get_watched_linecount(void)
{
  return linecount;
}


# line 54 "inputfile.c"

#ifdef _cplusplus
}
#endif
