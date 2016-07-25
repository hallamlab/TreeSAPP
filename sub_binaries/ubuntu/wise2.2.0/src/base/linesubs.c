#ifdef _cplusplus
extern "C" {
#endif
#include "linesubs.h"

/* this module is a little hacky...       */
/* It has buffers for loading in lines    */
/* and will edit them with the keys given */



#define MAXEDITLINE 512
#define MAXREPLACEMENT 256

static char input[MAXEDITLINE];
static char output[MAXEDITLINE];

static KeyandReplace kandr[MAXREPLACEMENT];
static int size = 0;

static int nest=0;


# line 33 "linesubs.dy"
void scan_and_replace(KeyandReplace * kr)
{
  char * runner;
  char * run2;
  char * o;
  
  if( (runner=strstr(input,kr->key)) == NULL)
    return; /* nothing to do */
  
  if( nest > 100 ) {
    warn("scan and replaced one line more than 100 times; dropping out!");
    return;
  }
  
  /*** copy upto runner in input ***/
  
  for(run2=input,o=output; run2 < runner && *run2 != '\0';)
    *(o++)=(*(run2++));
  
  /*** copy replace into output ***/
  
  for(run2=kr->replace; *run2 != '\0';)
    *(o++)=(*(run2++));
  
  /*** find end of key ***/
  
  runner += strlen(kr->key);
  
  /*** copy rest of string into output ****/
  
  
  for(run2=runner; *run2 != '\0';)
    *(o++)=(*(run2++));
  
  /*** add terminal 0 ****/
  
  *o = '\0';
  
  /*** move back to input ***/
  
  strcpy(input,output);
  
  /*** call again in case there are more to do ***/
  
  nest++;
  scan_and_replace(kr);
  
  
  return;
}

# line 84 "linesubs.dy"
char * scan_and_replace_line(char * line)
{
  register int i;


  /*** load in line ***/
  
  strcpy(input,line);
  
  
  for(i=0;i<size;i++) {
    nest = 0;
    scan_and_replace(kandr+i);
  }
  
  return input;
}




# line 105 "linesubs.dy"
boolean push_scan_and_replace_pair(char * key,char * replace)
{
  if( size > MAXREPLACEMENT ) {
    log_full_error(WARNING,0,"Scan and replace overflow for %s key\n",key);
    return FALSE;
  }
  
  kandr[size].key = stringalloc(key);
  kandr[size].replace = stringalloc(replace);
  
  size++;
  
  return TRUE;
}

# line 120 "linesubs.dy"
boolean pop_scan_and_replace_pair(void)
{
  if( size == 0 ) {
    warn("Attempting to pop an empty scan and replace stack");
    return FALSE;
  }
  ckfree(kandr[size-1].key);
  ckfree(kandr[size-1].replace);
  
  size--;

  return TRUE;
}

# line 134 "linesubs.dy"
void flush_scan_and_replace(void)
{
  register int i;
  
  
  for(i=0;i < size;i++) {
    ckfree(kandr[i].key);
    ckfree(kandr[i].replace);
  }
  size = 0;
  
  return;
}

# line 148 "linesubs.dy"
void read_plain_scan_and_replace(char * filename)
{
  FILE * ifp;
  
  ifp = openfile(filename,"r");
  
  if( ifp == NULL )
    {
      log_full_error(WARNING,0,"Cannot open filename %s \n",filename);
      return;
    }
  
  read_scan_and_replace_file(ifp,"endofscanfile");
  
  fclose(ifp);
  
  return;
}

# line 167 "linesubs.dy"
void read_scan_and_replace_file(FILE * ifp,char * endstring)
{
  char buffer[MAXLINE];
  char ** base,**brk;
  char *  key;
  char * runner;


  while( fgets(buffer,MAXLINE,ifp) != NULL)
    {
      
      /*	fprintf(stderr,"Read line %s\n",buffer); */
      
      if( endstring != NULL && endstring[0] != '\0' && strstartcmp(buffer,endstring) == 0)
	break;
      if( buffer[0] == '#' || buffer[0] == '!')
	continue;
      base = brk = breakstring_protect(buffer,spacestr,"@");
      if( base == NULL || *brk == NULL )
	{
	  log_full_error(WARNING,0,"Picked up problem in read scan and replace file");
	  continue;
	}
      
      key = *brk;
      
      brk++;
      
      if( *brk == NULL )
	{
			log_full_error(WARNING,0,"Picked up problem in read scan and replace file key %s has no replace",key);
			continue;
	}
      
      if( **brk == '@' )
	(*brk)++;
      
      runner = (*brk) + strlen((*brk));
      runner--;
      if( *runner == '@')
	*runner = '\0';
      
      
      
      push_scan_and_replace_pair(key,*brk);
      
      ckfree(base);
      
    }
  
  
  return;
}

		






# line 225 "linesubs.c"

#ifdef _cplusplus
}
#endif
