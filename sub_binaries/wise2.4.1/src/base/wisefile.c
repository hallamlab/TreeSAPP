#ifdef _cplusplus
extern "C" {
#endif
#include "wisefile.h"


#ifdef UNIX
static boolean isvms = FALSE;
#else
static boolean  isvms= TRUE;
#endif

static char * systemconfigdir=NULL;
static char * personaldir=NULL;
static char * homedir=NULL;
static boolean hasloaded=FALSE;


/* Function:  set_config_dir(path,*path)
 *
 * Descrip:    Programmatically set systemconfigdir to override
 *             any value set (or not) by env.var. WISECONFIGDIR.
 *
 *             Added by arve.
 *
 *
 * Arg:         path [UNKN ] path that WISECONFIGDIR is set to [NullString]
 * Arg:        *path [UNKN ] Undocumented argument [char]
 *
 */
# line 58 "wisefile.dy"
void set_config_dir(char *path) 
{
  if (systemconfigdir != NULL ) {
    ckfree(systemconfigdir);
  }
  systemconfigdir = stringalloc(path);
}


/* Function:  myfclose(ofp)
 *
 * Descrip:    reports the fclose type etc
 *
 *
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 70 "wisefile.dy"
int myfclose(FILE * ofp)
{
  fprintf(stderr,"Closing %d\n",ofp);
#undef fclose
  fclose(ofp);
  return 1;
}

    
/* Function:  try_to_load(void)
 *
 * Descrip:    Loads up 'standard' path places
 *
 *
 *
 */
# line 83 "wisefile.dy"
void try_to_load(void)
{
  char * runner;
  char buffer[256];
  
  if( hasloaded == TRUE)
    log_full_error(PEDANTIC,0,"Trying to reload configdirs");
  
  
  if( (runner=getenv("WISECONFIGDIR")) != NULL) {
    if( runner[strlen(runner)-1] == '/' ) 
      systemconfigdir=stringalloc(runner);	
    else {
      sprintf(buffer,"%s/",runner);
      systemconfigdir=stringalloc(buffer);
    }
      
  }
  
  if( (runner=getenv("WISEPERSONALDIR")) != NULL) {
    personaldir=stringalloc(runner);
  }

  if( (runner=getenv("HOME")) != NULL) {
    homedir=stringalloc(runner);
  }
  
  hasloaded=TRUE;
}

/* Function:  remove_file(filename)
 *
 * Descrip:    silly function to provide a boolean wrapper
 *             around remove. 
 *
 *
 * Arg:        filename [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 117 "wisefile.dy"
boolean remove_file(char * filename)
{
  if( remove(filename) == 0)
    return TRUE;
  else	return FALSE;
}

/* Function:  move_file(from,to)
 *
 * Descrip:    silly function to provide a boolean wrapper
 *             around rename 
 *
 *
 * Arg:        from [UNKN ] Undocumented argument [char *]
 * Arg:          to [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 128 "wisefile.dy"
boolean move_file(char * from,char * to)
{
  if( rename(from,to) == 0)
    return TRUE;
  else	return FALSE;
}

/* Function:  append_file_to_path(buffer,len,file,path)
 *
 * Descrip:    Appends file onto path in buffer
 *
 *
 * Arg:        buffer [UNKN ] Undocumented argument [char *]
 * Arg:           len [UNKN ] Undocumented argument [int]
 * Arg:          file [UNKN ] Undocumented argument [const char *]
 * Arg:          path [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 139 "wisefile.dy"
boolean append_file_to_path(char * buffer,int len,const char * file,char * path)
{
  char * runner;

  if( 1+strlen(file)+strlen(path) > len) {
    warn("Unable to expand %s with %s due to lack of buffer space",path,file);
    return FALSE;
  }

  strcpy(buffer,path);
  runner=buffer+strlen(buffer);
	
#ifdef UNIX
  if( *runner != '/') {
    *(runner++)='/';
    *runner='\0';
  }
#endif

  strcat(buffer,file);
  
  return TRUE;
}

/* Function:  touchfile(filename)
 *
 * Descrip:    sees if filename exists
 *
 *
 * Arg:        filename [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 166 "wisefile.dy"
boolean touchfile(char * filename)
{
  FILE * temp;
  boolean retval;
  
  
  if( filename == NULL ) {
    warn("Tried to touch a NULL filename");
    return FALSE;
  }
  
  if( (temp=openfile(filename,"r")) == NULL) {
    retval=FALSE;
  }
  else	{
    fclose(temp);
    retval=TRUE;
  }
  
  
  return retval;
}

/* Function:  openfile(filename,passedprot)
 *
 * Descrip:    Every file open goes through this.
 *
 *             It opens for reading in the following order 
 *                .
 *                WISEPERSONALDIR
 *                WISECONFIGDIR
 *
 *             For writing it opens in .
 *
 *             Filenames with ~'s are expanded to HOME/filename
 *
 *
 * Arg:          filename [UNKN ] filename to open for read/writing [const char *]
 * Arg:        passedprot [UNKN ] string representing standard fopen attributes [const char *]
 *
 * Return [UNKN ]  open'd filehandle, NULL on error [FILE *]
 *
 */
# line 205 "wisefile.dy"
FILE * openfile(const char * filename,const char * passedprot)
{
  FILE * ifp;
  char buffer[MAXPATHLEN];
  char prot[12]; /* protection string longer than 12 ! */
  boolean shouldreporterror=FALSE;
  
  strncpy(prot,passedprot,12);

  if( *prot == 'R' )  {
    shouldreporterror=TRUE;
    *prot = 'r';
  }
  
  if( *prot == 'W' ){
    shouldreporterror=TRUE;
    *prot = 'w';
  }
  
  if( filename == NULL) {
      if( shouldreporterror ) 
	info("Open for file failed due to NULL filename");
      return NULL;
    }	
  
  if( strcmp(filename,"-") == 0 ) {
    if( *prot == 'r' )
      return stdin;
    else	return stdout;
  }
  
  if( (ifp=fopen(filename,prot)) != NULL) {
#ifdef FILE_DEBUG
    fprintf(stderr,"Succeeded in opening %s. Filehandle is %d\n",filename,ifp);
#endif
    return ifp;
  }
  
  if( shouldreporterror ) {
    info("Direct open for [%s,%s] failed: Error message is %s",filename,prot, ERRORSTR);
  }
  
  
  if( hasloaded == FALSE)
    try_to_load();
  
  if( isvms != TRUE && homedir != NULL && filename[0]== '~' && filename[1] == '/')
    if( append_file_to_path(buffer,MAXPATHLEN,filename+2,homedir) != FALSE )
      {
	/*fprintf(stderr,"Trying to load %s\n",buffer);*/
	if( (ifp=fopen(buffer,prot)) != NULL)
	  return ifp;
	else if( shouldreporterror ) 
	  info("Expanded tilda open for[%s,%s], expanded to %s failed: Error message is %s",filename,prot,buffer, ERRORSTR);
      }


	
  if( strchr(prot,'w') != NULL)
    return NULL;
  
  
  /* next line ABSOLUTELY relies on order of evaluation */
  if( personaldir != NULL && append_file_to_path(buffer,MAXPATHLEN,filename,personaldir) != FALSE && (ifp=fopen(buffer,prot)) != NULL)
    return ifp;
  
  else if ( shouldreporterror )
    log_full_error(INFO,0,"Expanded personal direcotry open for[%s,%s], expanded to %s failed: Error message is %s",filename,prot,buffer,ERRORSTR);
  
  
  /* next line ABSOLUTELY relies on order of evaluation */
  if( systemconfigdir != NULL && append_file_to_path(buffer,MAXPATHLEN,filename,systemconfigdir) != FALSE && (ifp=fopen(buffer,prot)) != NULL)
    return ifp;
  else if ( shouldreporterror )
    log_full_error(INFO,0,"Expanded system directory open for[%s,%s], expanded to %s failed: Error message is %s",filename,prot,buffer, ERRORSTR);
  
  
#ifdef CAREFUL
  if( looks_like_vms(filename) && isvms == FALSE)
    log_full_error(WARNING,0,"Filename %s looks like a VMS "
		   "file to me, and you've compiled as UNIX",filename);
  if( looks_like_unix(filename) && isvms == TRUE)
    log_full_error(WARNING,0,"Filename %s looks like a UNIX "
		   "file to me, and you've compiled as VMS",filename);
#endif
  
  return NULL;
  
}

/* Function:  envopenfile(envname,filename,name,env)
 *
 * Descrip:    This function basically mirrors the function in file.c
 *             in HMMer2. You call it as
 *
 *               fp = Envfile(filename,envname);
 *
 *               where envname looks like "BLASTDB" etc.
 *
 *
 *
 * Arg:         envname [READ ] enviroment variable to read from [NullString]
 * Arg:        filename [UNKN ] Undocumented argument [char *]
 * Arg:            name [READ ] filename to open [NullString]
 * Arg:             env [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  a valid file pointer or NULL [FILE *]
 *
 */
# line 308 "wisefile.dy"
FILE * envopenfile(char * filename,char * env)
{
  char * envp;
  char path [512];

  if( filename == NULL || env == NULL ) {
    warn("Passed a NULL filename or enviroment name into Envfile. Should trap this elsewhere");
    return NULL;
  }

  if( (envp = getenv(env)) == NULL ) {
    /* fail gracefully - somebody might query a number of enviroment variables */
    return NULL;
  }

  if( strlen(filename) + strlen(envp) < 490 ) {
    warn("Really long filename/enviroment variables [%s] [%s] Can't cope!",filename,envp);
    return NULL;
  }
  sprintf(path,"%s/%s",filename,envp);

  return fopen(path,"r");
}
  

# line 347 "wisefile.c"

#ifdef _cplusplus
}
#endif
