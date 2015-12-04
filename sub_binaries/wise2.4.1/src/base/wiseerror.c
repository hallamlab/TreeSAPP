#ifdef _cplusplus
extern "C" {
#endif
#include "wiseerror.h"

static Flag fatal_flag    = 3;
static Flag warning_flag  = 3;
static Flag info_flag     = 3;
static Flag report_flag   = 3;

static ErrorDisplayType display_type = ERROR_DISPLAY_PROGRAM;
static char * log_display_string = NULL;

#define flag_of_type(c) (c == FATAL ? fatal_flag : c == WARNING ? warning_flag :c == INFO ? info_flag : report_flag )


static int eventc=0;
static FILE * errlog=NULL;

static void (*error_call)(char *,int)= NULL;

static char * error_msg_stack[MAXMSGSTACKERROR];

static char * (*((error_msg_call)[MAXMSGSTACKERROR]))(void);
static int  msg_stack_no=0;

/* Function:  set_error_display_style(t)
 *
 * Descrip:    sets error display style
 *
 *
 * Arg:        t [UNKN ] Undocumented argument [ErrorDisplayType]
 *
 */
# line 85 "wiseerror.dy"
void set_error_display_style(ErrorDisplayType t)
{
  display_type = t;
}

/* Function:  set_log_display_string(str)
 *
 * Descrip:    sets string for log display type 
 *
 *
 * Arg:        str [UNKN ] Undocumented argument [char *]
 *
 */
# line 93 "wiseerror.dy"
void set_log_display_string(char * str)
{
  if( log_display_string != NULL ) {
    ckfree(log_display_string);
  }
  log_display_string = stringalloc(str);
}



/* Function:  push_errormsg_stack(msg,)
 *
 * Descrip:    This adds a function onto the error message
 *             stack for stacking errors. this is for
 *             decent parsers etc for allowing 'error scope'
 *             to be propagated down.
 *
 *             It is very very bad form to push an errormsg
 *             stack and not pop it at the end. 
 *
 *
 * Arg:        msg [UNKN ] Undocumented argument [char *]
 * Arg:            [UNKN ] Undocumented argument [.]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 113 "wiseerror.dy"
boolean push_errormsg_stack(char * msg, ...)
{
  char buffer[1024];
  va_list ap;
  
  
  
  va_start(ap,msg);
  vsprintf(buffer,msg,ap);	
  
  if( msg_stack_no >= MAXMSGSTACKERROR ) {
    warn("Too many messages held on stack, [%s] discarded\n",buffer);
    /*** still should up the number ***/
    msg_stack_no++;
    return FALSE;
  }
  error_msg_call[msg_stack_no] = NULL;
  error_msg_stack[msg_stack_no++] = stringalloc(buffer);
  return TRUE;
}	

  

/* Function:  push_errormsg_stack_call(ecall)
 *
 * Descrip:    This adds a function call for people who want
 *             to register error handling functions
 *
 *             Probably best wrapped by a separate function
 *
 *
 * Arg:        ecall [UNKN ] Undocumented argument [NullString]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 142 "wiseerror.dy"
boolean push_errormsg_stack_call( char * (*ecall)(void))
{

  if( msg_stack_no >= MAXMSGSTACKERROR ) {
    warn("Too many messages held on stack, Error message call discarded\n");
    /*** still should up the number ***/
    msg_stack_no++;
    return FALSE;
  }
  error_msg_call[msg_stack_no] = ecall;
  error_msg_stack[msg_stack_no++] = NULL;
  return TRUE;
}

/* Function:  pop_errormsg_stack(void)
 *
 * Descrip:    This removes a error message from the stack
 *
 *
 *
 */
# line 159 "wiseerror.dy"
void pop_errormsg_stack(void)
{
  if( msg_stack_no < MAXMSGSTACKERROR && msg_stack_no > 0 && error_msg_stack[msg_stack_no-1] != NULL)
    ckfree(error_msg_stack[msg_stack_no-1]);
  if( msg_stack_no > 0)
    msg_stack_no--;
}

/* Function:  show_message_stack(ofp)
 *
 * Descrip:    This shows the error message set to
 *             everyone 
 *
 *
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 172 "wiseerror.dy"
void show_message_stack(FILE * ofp)
{
  register int i;
  register int j;
  register int count;
  
  for(i=0;i<msg_stack_no;i++) {
    if( error_msg_call[i] != NULL ) {
      show_text( (*(error_msg_call[i]))(),65,ofp);
    } else {
      show_text(error_msg_stack[i],65,ofp);
    }
  }
}


/* Function:  add_log_file(ofp)
 *
 * Descrip:    Makes ofp the log file for the errors. Discards
 *             any previous one
 *
 *
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 193 "wiseerror.dy"
void add_log_file(FILE * ofp)
{
  /* you must also switch on logging errors as well! */

  errlog=ofp;
}


/* Function:  add_log_filename(filename)
 *
 * Descrip:    Opens filename as the log file. 
 *
 *
 * Arg:        filename [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 204 "wiseerror.dy"
boolean add_log_filename(char * filename)
{
  register FILE * ofp;
  if( (ofp=openfile(filename,"w")) == NULL) {
    warn("Could not open %s as a log filename",filename);
    return FALSE;
  }
  
  add_log_file(ofp);
  return TRUE;
}

/* Function:  error_off(type)
 *
 * Descrip:    Really for the API. Wraps some
 *             of the error concepts..
 *
 *
 * Arg:        type [UNKN ] Undocumented argument [int]
 *
 */
# line 220 "wiseerror.dy"
void error_off(int type)
{
  error_flag_off(type,ERRORUSE);
}

/* Function:  error_on(type)
 *
 * Descrip:    Really for the API. Wraps some
 *             of the error concepts..
 *
 *
 * Arg:        type [UNKN ] Undocumented argument [int]
 *
 */
# line 229 "wiseerror.dy"
void error_on(int type)
{
  error_flag_on(type,ERRORUSE);
}


/* Function:  error_flag_on(type,f)
 *
 * Descrip:    Turns on the particular type of error flag
 *             (eg, STDERR etc). 
 *
 *
 * Arg:        type [UNKN ] Undocumented argument [int]
 * Arg:           f [UNKN ] Undocumented argument [Flag]
 *
 */
# line 239 "wiseerror.dy"
void error_flag_on(int type,Flag f)
{
  switch (type) {
    case FATAL    : fatal_flag     = fatal_flag | f; return;
    case WARNING  : warning_flag   = warning_flag | f; return;
    case INFO     : info_flag      = info_flag | f; return;
    case REPORT   : report_flag    = report_flag | f ; return;
    default : log_full_error(WARNING,0,"In error system, tried to change flag %d which doesn't exist!",type); return;
      
    }
}

/* Function:  error_flag_off(type,f)
 *
 * Descrip:    Turns off the particular type of error flag
 *             (eg, STDERR etc). 
 *
 *
 * Arg:        type [UNKN ] Undocumented argument [int]
 * Arg:           f [UNKN ] Undocumented argument [Flag]
 *
 */
# line 255 "wiseerror.dy"
void error_flag_off(int type,Flag f)
{
  switch (type) {
    case FATAL    : fatal_flag     = fatal_flag & ~f; return;
    case WARNING  : warning_flag   = warning_flag & ~f; return;
    case INFO     : info_flag      = info_flag & ~f; return;
    case REPORT   : report_flag    = report_flag & ~f ; return;
    default : log_full_error(WARNING,0,"In error system, tried to change flag %d which doesn't exist!",type); return;
    }
}

/* Function:  catch_errors(catch)
 *
 * Descrip:    This is a wrapper for the error handling
 *             system. It does the following things
 *
 *             Sets function as the function to process errors
 *
 *             Switches the INFO,ERROR and FATAL flags off on STDERR
 *             and on to ERROR CALL.
 *
 *
 * Arg:        catch [FUNCP]  [void (*catch]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 277 "wiseerror.dy"
boolean catch_errors(void (*catch)(char *,int))
{
  erroron(INFO);
  erroron(WARNING);
  erroron(FATAL);

  errorstderroff(INFO);
  errorstderroff(WARNING);
  errorstderroff(FATAL);

  errorcallon(INFO);
  errorcallon(WARNING);
  errorcallon(FATAL);

  push_error_call(catch);

  return TRUE;
}

/* Function:  stop_catching_errors(void)
 *
 * Descrip:    Switches off error catching,
 *             putting flags back on to STDERR
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 300 "wiseerror.dy"
boolean stop_catching_errors(void)
{
  errorstderron(INFO);
  errorstderron(WARNING);
  errorstderron(FATAL);

  errorcalloff(INFO);
  errorcalloff(WARNING);
  errorcalloff(FATAL);

  pop_error_call();

  return TRUE;
}

/* Function:  push_error_call()
 *
 * Descrip:    Registers this function for dealing with errors
 *
 *             Try to use catch_errors instead
 *
 *
 *
 * Arg:         [UNKN ] Undocumented argument [NullString]
 *
 */
# line 321 "wiseerror.dy"
void push_error_call(void (* func)(char *,int))
{	
  error_call=func;
}

/* Function:  pop_error_call(void)
 *
 * Descrip:    Discards current function for dealing with errors
 *
 *
 *
 */
# line 329 "wiseerror.dy"
void pop_error_call(void)
{
  error_call=NULL;
}

/* Function:  type_to_error(type)
 *
 * Descrip:    Turns int error types to Names
 *             for display purposes.
 *
 *
 * Arg:        type [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [char *]
 *
 */
# line 338 "wiseerror.dy"
char * type_to_error(int type)
{
  switch (type) {
  case FATAL :  return "Fatal Error";
  case WARNING : return "Warning Error";
  case INFO   :  return "Information";
  case REPORT :  return "Debug Report";
  default	: return "Strange error type???";
  }
}

/* Function:  info(msg,)
 *
 * Descrip:    Produces a 'info' error message.
 *
 *
 * Arg:        msg [UNKN ] Undocumented argument [char *]
 * Arg:            [UNKN ] Undocumented argument [.]
 *
 */
# line 352 "wiseerror.dy"
void info(char * msg, ...)
{
  char buffer[1024];
  va_list ap;
  int type = INFO;
  
  va_start(ap,msg);
  vsprintf(buffer,msg,ap);
  
  eventc++;	
  
  show_error(flag_of_type(type),buffer,type);
  
  return;
}		

/* Function:  debug_report(msg,)
 *
 * Descrip:    Produces a 'debug report' error message.
 *
 *
 * Arg:        msg [UNKN ] Undocumented argument [char *]
 * Arg:            [UNKN ] Undocumented argument [.]
 *
 */
# line 371 "wiseerror.dy"
void debug_report(char * msg, ...)
{
  char buffer[1024];
  va_list ap;
  int type = REPORT;
  
  va_start(ap,msg);
  vsprintf(buffer,msg,ap);
  
  eventc++;	
  
  show_error(flag_of_type(type),buffer,type);
  
  return;
}		

/* Function:  warn(msg,)
 *
 * Descrip:    Produces a 'warn' error message.
 *
 *
 * Arg:        msg [UNKN ] Undocumented argument [char *]
 * Arg:            [UNKN ] Undocumented argument [.]
 *
 */
# line 390 "wiseerror.dy"
void warn(char * msg, ...)
{
  char buffer[1024];
  va_list ap;
  int type = WARNING;
  
  va_start(ap,msg);
  vsprintf(buffer,msg,ap);
  
  eventc++;	
  
  show_error(flag_of_type(type),buffer,type);
  
  return;
}		

/* Function:  fatal(msg,)
 *
 * Descrip:    Produces a 'fatal' error message.
 *
 *
 * Arg:        msg [UNKN ] Undocumented argument [char *]
 * Arg:            [UNKN ] Undocumented argument [.]
 *
 */
# line 409 "wiseerror.dy"
void fatal(char * msg, ...)
{
  char buffer[1024];
  va_list ap;
  int type = FATAL;
  
  va_start(ap,msg);
  vsprintf(buffer,msg,ap);
  
  eventc++;	
  
  show_error(flag_of_type(type),buffer,type);
  fputc('\n',stderr);
  exit(2);
  
  return;
}		

/* Function:  log_full_error(type,location,msg,)
 *
 * Descrip:    Deprecated
 *
 *             produces any of the error types
 *
 *
 * Arg:            type [UNKN ] Undocumented argument [int]
 * Arg:        location [UNKN ] Undocumented argument [int]
 * Arg:             msg [UNKN ] Undocumented argument [char *]
 * Arg:                 [UNKN ] Undocumented argument [.]
 *
 */
# line 432 "wiseerror.dy"
void log_full_error(int type,int location,char * msg, ...)
{
  char buffer[1024];
  va_list ap;
  
  va_start(ap,msg);
  vsprintf(buffer,msg,ap);
  
  eventc++;
  
  if(type == FATAL) {
    show_error(fatal_flag,buffer,FATAL);
    fputc('\n',stderr);
    exit(2);
  }
  
  show_error(flag_of_type(type),buffer,type);
  
  return;
}		

/* Function:  start_reporting(msg,)
 *
 * Descrip:    Starts a % reporting run. This is the header message
 *
 *
 * Arg:        msg [UNKN ] Undocumented argument [char *]
 * Arg:            [UNKN ] Undocumented argument [.]
 *
 */
# line 456 "wiseerror.dy"
void start_reporting(char * msg,...)
{
  char buffer[1024];
  va_list ap;
  
  if( !(report_flag & ERRORUSE) )
    return;
  
  va_start(ap,msg);
  vsprintf(buffer,msg,ap);
  
  fputs(buffer,stderr);
  start_overlay(stderr);
}

/* Function:  stop_reporting(void)
 *
 * Descrip:    Stops a % reporting run. 
 *
 *
 *
 */
# line 474 "wiseerror.dy"
void stop_reporting(void)
{
  
  if( !(report_flag & ERRORUSE) )
    return;
  stop_overlay();
}

/* Function:  show_error(flag,othermsg,type)
 *
 * Descrip:    Actually shows the error
 *
 *
 * Arg:            flag [UNKN ] Undocumented argument [Flag]
 * Arg:        othermsg [UNKN ] Undocumented argument [char *]
 * Arg:            type [UNKN ] Undocumented argument [int]
 *
 */
# line 486 "wiseerror.dy"
void show_error(Flag flag,char * othermsg,int type)
{
  time_t curr_time;
  struct tm * loctime;
  char time_buffer[50];

  if( !(flag & ERRORUSE) )
    return;

  
  if( type == REPORT ) {
    print_overlay(othermsg);
    return;
  }
  
  if(flag&ERRORTOSTDERR) {
    if( display_type == ERROR_DISPLAY_SERVER ) {
      curr_time = time(NULL);
      loctime = localtime(&curr_time);
      
      asctime_r(loctime,time_buffer);
      time_buffer[24] = '\0';

      fprintf(stderr,"%s: [%s] [%s] %s\n",type_to_error(type),
	      log_display_string == NULL ? "[no log string]" : log_display_string,
	      time_buffer,
	      othermsg);
    } else {
      fputs(type_to_error(type),stderr);
      fputc('\n',stderr);
      if( msg_stack_no > 0 )
	show_message_stack(stderr);
      fputs("\t",stderr);
      if( type == FATAL )
	fputs(othermsg,stderr);
      else	show_text(othermsg,70,stderr);
    }
  }
  
  if( flag&ERRORTOLOG && errlog != NULL) {
    fputs(type_to_error(type),errlog);
    fputc('\n',stderr);
    if( msg_stack_no > 0 )
      show_message_stack(errlog);
    fputs("\n\t",errlog);
    show_text(othermsg,70,errlog);
  }
	
  if( flag&ERRORTOCALL && error_call != NULL)
    {
      (*(error_call))(othermsg,type);
    }
  
  
  
  return;
}


# line 604 "wiseerror.c"

#ifdef _cplusplus
}
#endif
