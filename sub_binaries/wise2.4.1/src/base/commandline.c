#ifdef _cplusplus
extern "C" {
#endif
#include "commandline.h"


/* Function:  strip_out_remaining_options_with_warning(argc,argv)
 *
 * Descrip:    This removes all remaining options, issuing warnings
 *             through warn
 *
 *
 * Arg:        argc [UNKN ] Undocumented argument [int *]
 * Arg:        argv [UNKN ] Undocumented argument [char **]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 33 "commandline.dy"
boolean strip_out_remaining_options_with_warning(int * argc,char ** argv)
{
  register int i;
  boolean ret = FALSE;


  for(i=0;i<*argc;i++) {
    if( argv[i][0] != '-' )
      continue;

    /*** ignore bare '-' - could be arguments ***/

    if( argv[i][1] == '\0')
      continue;

    warn("You have an unrecognised argument, %s, removing",argv[i]);
    
    memmove(argv+i,argv+i+1,sizeof(char *)*(*argc - i -1));
    *argc -= 1;
    ret = TRUE;
    }

  return ret;
}


/* Function:  strip_out_boolean_argument(argc,argv,tag)
 *
 * Descrip:    This removes argument in tag from the commandline if there and
 *             returns TRUE
 *
 *             otherwise returns FALSE
 *
 *
 * Arg:        argc [UNKN ] argc from main declaration (call as &argc) [int *]
 * Arg:        argv [UNKN ] argv from main declaration [char **]
 * Arg:         tag [READ ] -[string] argument to find [char *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 69 "commandline.dy"
boolean strip_out_boolean_argument(int * argc,char ** argv,char * tag)
{
  register int i;

  for(i=0;i<*argc;i++) {
    if( argv[i][0] != '-' )
      continue;
    
    /** ignore '-' by themselves **/

    if( argv[i][1] == '\0')
      continue;

    /** try to strcmp with tag **/

    if( strcmp(tag,argv[i]+1) == 0 ) {
      memmove(argv+i,argv+i+1,sizeof(char *)*(*argc - i -1));

      *argc -= 1;
      
      /* return */

      return TRUE;
    }
  }

  return FALSE;
}

/* Function:  strip_out_boolean_def_argument(argc,argv,tag,value)
 *
 * Descrip:    This sets a boolean argument as
 *                -tag   == TRUE
 *                -notag == FALSE
 *
 *             This is probably better than strip_out_boolean_argument
 *             as 
 *               defaults can be provided
 *               user can switch either on or off
 *
 *
 * Arg:         argc [UNKN ] argc from main declaration (call as &argc) [int *]
 * Arg:         argv [UNKN ] argv from main declaration [char **]
 * Arg:          tag [READ ] -[string] argument to find [char *]
 * Arg:        value [WRITE] boolean pointer to write value to [boolean *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 113 "commandline.dy"
boolean strip_out_boolean_def_argument(int * argc,char ** argv,char * tag,boolean * value)
{
  int arg;
  char buffer[64];

  if( (arg = strip_out_boolean_argument(argc,argv,tag)) == TRUE ) {
    *value = TRUE;
  }

  sprintf(buffer,"no%s",tag);
  if( (arg = strip_out_boolean_argument(argc,argv,buffer)) == TRUE ) {
    *value = FALSE;
  }

  return TRUE;
}

/* Function:  strip_out_assigned_argument(argc,argv,tag)
 *
 * Descrip:    This removes argument in tag from the commandline if there and
 *             returns the argument to it (in -tag arg - ie with a space).
 *
 *             otherwise returns NULL
 *
 *
 * Arg:        argc [UNKN ] argc from main declaration (call as &argc) [int *]
 * Arg:        argv [UNKN ] argv from main declaration [char **]
 * Arg:         tag [READ ] -[string] argument to find [char *]
 *
 * Return [UNKN ]  Undocumented return value [char *]
 *
 */
# line 140 "commandline.dy"
char * strip_out_assigned_argument(int * argc,char ** argv,char * tag)
{
  register int i;
  register char * ret;


  for(i=0;i<*argc;i++) {
    if( argv[i][0] != '-' )
      continue;
    
    /** ignore '-' by themselves **/

    if( argv[i][1] == '\0')
      continue;

    /** try to strcmp with tag **/

    if( strcmp(tag,argv[i]+1) == 0 ) {
      /* pick up next argument: if it has a '-' treat that as the argument. */
      
      if( i+1 >= *argc ) {
	warn("In processing command line [%s], tag [%s] expects an argument",argv[0],tag);
	return NULL; /* give 'em nothing */
      }


      /* assign return value */

      ret = argv[i+1];

      /* ok, now remove both tag and argument from the line */

      memmove(argv+i,argv+i+2,sizeof(char *)*(*argc - i -2));

      *argc -= 2;
      
      /* return */

      return ret;
    }
  }

  return NULL;
}


/* Function:  strip_out_integer_argument(argc,argv,tag,value)
 *
 * Descrip:    This removes argument in tag from the commandline if there, and
 *             looks at the argument as whether it is an int or not.
 *
 *             No tag - returns FALSE and leaves value alone
 *             Tag but no integer value - issues warning, returns FALSE
 *             tag but integer value    - sets value, returns true
 *
 *
 * Arg:         argc [UNKN ] argc from main declaration (call as &argc) [int *]
 * Arg:         argv [UNKN ] argv from main declaration [char **]
 * Arg:          tag [READ ] -[string] argument to find [char *]
 * Arg:        value [WRITE] int pointer to write value to [int *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 199 "commandline.dy"
boolean strip_out_integer_argument(int * argc,char ** argv,char * tag,int * value)
{
  char * arg;

  if( (arg = strip_out_assigned_argument(argc,argv,tag)) == NULL )
    return FALSE;

  if( is_integer_string(arg,value) == FALSE ) {
    warn("Argument [%s] to [%s] is not an integer. Not changing the value [%d]",arg,tag,value);
    return FALSE;
  }

  return TRUE;
}

/* Function:  strip_out_float_argument(argc,argv,tag,value)
 *
 * Descrip:    This removes argument in tag from the commandline if there, and
 *             looks at the argument as whether it is a double or not.
 *
 *             No tag - returns FALSE and leaves value alone
 *             Tag but no integer value - issues warning, returns FALSE
 *             tag but integer value    - sets value, returns true
 *
 *
 * Arg:         argc [UNKN ] argc from main declaration (call as &argc) [int *]
 * Arg:         argv [UNKN ] argv from main declaration [char **]
 * Arg:          tag [READ ] -[string] argument to find [char *]
 * Arg:        value [WRITE] double pointer to write value to [double *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 227 "commandline.dy"
boolean strip_out_float_argument(int * argc,char ** argv,char * tag,double * value)
{
  char * arg;

  if( (arg = strip_out_assigned_argument(argc,argv,tag)) == NULL )
    return FALSE;

  if( is_double_string(arg,value) == FALSE ) {
    warn("Argument [%s] to [%s] is not a double. Not changing the value [%f]",arg,tag,value);
    return FALSE;
  }

  return TRUE;
}
   
/* Function:  show_standard_options(ofp)
 *
 * Descrip:    Shows standard options
 *
 *
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 245 "commandline.dy"
void show_standard_options(FILE * ofp)
{
  fprintf(ofp,"Standard options\n");
  fprintf(ofp,"  -help       help\n");
  fprintf(ofp,"  -version    show version and compile info\n");
  fprintf(ofp,"  -silent     No messages    on stderr\n");
  fprintf(ofp,"  -quiet      No report/info on stderr\n");
  fprintf(ofp,"  -erroroffstd No warning messages to stderr\n");
  fprintf(ofp,"  -errorlog   [file] Log warning messages to file\n");
  fprintf(ofp,"  -errorstyle [server/program] style of error reporting (default program)\n");
}
 

/* Function:  strip_out_standard_options(argc,show_version,show_help,argv)
 *
 * Descrip:    Handles default arguments correctly, setting help,version,errors
 *             handling
 *
 *
 * Arg:                argc [UNKN ] Undocumented argument [int *]
 * Arg:        show_version [UNKN ] Undocumented argument [NullString]
 * Arg:           show_help [UNKN ] Undocumented argument [NullString]
 * Arg:                argv [UNKN ] Undocumented argument [char **]
 *
 */
# line 262 "commandline.dy"
void strip_out_standard_options(int * argc,char ** argv,void (*show_help)(FILE * ofp),void (*show_version)(FILE * ofp))
{
  char * errlog;
  char * temp;


  if( (strip_out_boolean_argument(argc,argv,"help")) == TRUE ) {
    (*show_help)(stdout);
    exit(1);
  }

  if( (strip_out_boolean_argument(argc,argv,"version")) == TRUE ) {
    (*show_version)(stdout);
    exit(1);
  }

  if( (strip_out_boolean_argument(argc,argv,"silent")) == TRUE ) {
    erroroff(REPORT);
    erroroff(INFO);
    erroroff(WARNING);
  }

  if( (strip_out_boolean_argument(argc,argv,"quiet")) == TRUE ) {
    erroroff(REPORT);
    erroroff(INFO);
  }

  if( (strip_out_boolean_argument(argc,argv,"erroroffstd")) == TRUE ) {
    errorstderroff(WARNING);
  }

  if( (temp = strip_out_assigned_argument(argc,argv,"errorstyle")) != NULL ) {
    if( strcmp(temp,"server") == 0 ) {
      set_error_display_style(ERROR_DISPLAY_SERVER);
    } else if ( strcmp(temp,"program") == 0 ) {
      set_error_display_style(ERROR_DISPLAY_PROGRAM);
    } else {
      fprintf(stderr,"Unknown errorstyle %s",temp);
    }
  }

  set_log_display_string(argv[0]);


  if( (errlog=strip_out_assigned_argument(argc,argv,"errlog")) != NULL ) {
    if( add_log_filename(errlog) == FALSE ) {
      fatal("Could not use %s as a error log file\n",errlog);
    } else {
      warn("Logging errors to %s as well as stderr",errlog);
      errorlogon(WARNING);
    }
  }
}  

# line 331 "commandline.c"

#ifdef _cplusplus
}
#endif
