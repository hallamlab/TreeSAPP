/* file.c
 * SRE, Wed Jun 19 11:19:22 1996
 * 
 * File operation utilities, dealing with pathnames, directories,
 * and environment variables.
 * 
 * The goal is to have these be platform-independent but they
 * currently are UNIX-specific: i.e. this file is currently POSIX compliant
 * but it is NOT ANSI C compliant. (The sole offender is getenv().)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "squid.h"
#include "sqfuncs.h"

/* 
 * VMS:    #define DIRSLASH ']'
 * MacOS:  #define DIRSLASH ':'
 * DOS:    #define DIRSLASH '\\'
 * 
 * The code assumes that '.' is used for file name extensions,
 * such as "foo.bar".
 */
#define DIRSLASH '/'           /* UNIX directory paths have /foo/bar */



/* Function: FileDirname()
 * 
 * Purpose:  Returns the path from a filename:
 *             "/foo/bar/baz"  -> "/foo/bar"
 *             "foo/bar"       -> "foo" 
 *             "foo"           -> "."
 *             "/"             -> "/" 
 *           i.e. the string will be non-NULL; it will
 *           contain the string up to but not including the
 *           last '/' character; returns "." if
 *           there are no '/' characters, and returns "/"
 *           if the last slash is the first character.
 *           Modeled on Tcl's "file dirname" command.
 *
 * Args:     file - name of file "/foo/bar/baz".
 *           
 * Return:   ptr to malloc'ed string "/foo/bar".          
 */
char *
FileDirname(char *file)
{
  char *dirname;
  char *lastslash;
  int   len;
  
  lastslash = strrchr(file, DIRSLASH);
  len =  (lastslash == NULL) ? 0 : (int) (lastslash - file);
  dirname = (char *) MallocOrDie (sizeof(char) * (len+2));
  if (len > 0)                strncpy(dirname, file, len);
  else if (*file != DIRSLASH) { *dirname = '.';      len = 1; }
  else                        { *dirname = DIRSLASH; len = 1; }
  dirname[len] = '\0';
  return dirname;
}


/* Function: FileTail()
 * 
 * Purpose:  Return everything after the DIRSLASH:
 *             "/foo/bar/baz.1"  -> "baz.1"
 *             "foo/bar"         -> "bar" 
 *             "foo"             -> "foo"
 *             "/"               -> "" 
 *           If noextension is TRUE, removes a trailing ".foo" extension
 *           too.
 *           
 * Args:     file        - name of file "/foo/bar/baz.1"         
 *           noextension - TRUE to also remove extensions
 *           
 * Return:   ptr to malloc'ed string "baz.1"          
 */      
char * 
FileTail(char *file, int noextension)
{
  char *tail;
  char *lastslash;
  char *lastdot;
				/* remove directory prefix */
  lastslash = strrchr(file, DIRSLASH);
  tail = (char *) MallocOrDie (sizeof(char) * (strlen(file)+1));
  if (lastslash == NULL) strcpy(tail, file);
  else                   strcpy(tail, lastslash+1);
				/* remove trailing suffix */
  if (noextension) {
    if ((lastdot = strrchr(tail, '.')) != NULL)
      *lastdot = '\0';
  }

  return tail;
}


/* Function: FileConcat()
 * 
 * Purpose:  Concatenate a directory path and a file name,
 *           returning a pointer to a malloc'ed string with the
 *           full filename.
 */
char *
FileConcat(char *dir, char *file)
{
  char *full;

  full = (char *) MallocOrDie (sizeof(char) * (strlen(dir)+strlen(file)+2));
  if (*file == DIRSLASH) strcpy(full, file); /* file = "/foo", ignore directory. */
  else                   sprintf(full, "%s%c%s", dir, DIRSLASH, file);
  return full;
}


/* Function: EnvFileOpen()
 * Date:     Sun Feb 12 10:55:29 1995            
 * 
 * Purpose:  Open a file, given a file name and an environment
 *           variable that contains a directory path. Files
 *           are opened read-only. Does not look at current directory
 *           unless "." is explicitly in the path specified by env.
 *           
 *           For instance: 
 *             fp = EnvFileOpen("BLOSUM45", "BLASTMAT");
 *           or:
 *             fp = EnvFileOpen("swiss", "BLASTDB");  
 *             
 *           Environment variables may contain a colon-delimited
 *           list of more than one path; e.g.
 *             setenv BLASTDB /nfs/databases/foo:/nfs/databases/bar
 *             
 * Args:     fname - name of file to open
 *           env   - name of environment variable containing path
 *  
 * Return:   FILE * to open file, or NULL on failure -- same as fopen()
 */
FILE *
EnvFileOpen(char *fname, char *env)
{
  FILE *fp;
  char *path;
  char *s;                      /* ptr to indiv element in env list */
  char  full[1024];             /* constructed file name */

  if (env == NULL) return NULL;
  if ((path = Strdup(getenv(env))) == NULL) return NULL;
  
  s = strtok(path, ":");
  while (s != NULL)
    {
      if (((int) strlen(fname) + (int) strlen(s) + 2) > 1024) 
	{ free(path); return NULL; }
      sprintf(full, "%s%c%s", s, DIRSLASH, fname);
      if ((fp = fopen(full, "r")) != NULL) break;
      s = strtok(NULL, ":");
    }
  free(path);
  return fp;
}


/* Function: FileExists()
 * 
 * Purpose:  Return TRUE if filename exists.
 *           Testing fopen() is the only possible platform-independent test
 *           I'm aware of.  
 */
int
FileExists(char *filename)
{
  FILE *fp;
  if ((fp = fopen(filename, "r"))) { fclose(fp); return TRUE; }
  return FALSE;
}


