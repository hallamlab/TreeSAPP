#ifdef _cplusplus
extern "C" {
#endif
#include "wisestring.h"

/* Function:  get_number_from_slashed_string(qstr,slashstr)
 *
 * Descrip:    handy: pass a string like "xxx/yyy/zzz"    
 *             and a query like yyy, will return 2        
 *             (ie, maps positions in slashed string with 
 *             a number from start). Returns -1 if none   
 *
 *
 * Arg:            qstr [UNKN ] Undocumented argument [char *]
 * Arg:        slashstr [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 41 "wisestring.dy"
int    get_number_from_slashed_string(char * qstr,char * slashstr)
{
  char * runner;
  char * stop;
  int ret;

  /*** if the qstr is not in there, then we can get out! **/

  if( (stop=strstr(slashstr,qstr)) == NULL ) {
    return (-1);
  }

  /** ok, run through slashstr until we hit stop, counting slashes **/

  for(ret = 0,runner=slashstr;*runner && runner < stop ;runner++) 
    if( *runner == '/')
      ret++;

  return ret;
}


/* Function:  is_integer_string(string,val)
 *
 * Descrip:    checks that strings are ints or doubles    
 *             and then converts, storing value in val    
 *             if val == NULL, will not store (!)        
 *                                                         
 *             Does use sensible library functions        
 *             strtol...                          
 *
 *
 * Arg:        string [UNKN ] Undocumented argument [char *]
 * Arg:           val [UNKN ] Undocumented argument [int *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 71 "wisestring.dy"
boolean is_integer_string(char * string,int * val)
{
  int ret;
  char * end;

  ret = strtol(string,&end,10);
  
  if( val != NULL)
    *val = ret;

  if( isalpha(*end) )
    return FALSE;

  return TRUE;
}

/* Function:  is_double_string(string,val)
 *
 * Descrip:     checks that strings are doubles    
 *              and then converts, storing value in val    
 *               if val == NULL, will not store (!)        
 *                                                         
 *              Does use sensible library functions        
 *              strtod
 *
 *
 * Arg:        string [UNKN ] Undocumented argument [char *]
 * Arg:           val [UNKN ] Undocumented argument [double *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 95 "wisestring.dy"
boolean is_double_string(char * string,double * val)
{
  double ret;
  char * end;

  ret = strtod(string,&end);

  if( val != NULL)
    *val = ret;

  if( isalpha(*end) )
    return FALSE;
  
  return TRUE;
}

/* Function:  compress_space_around_punc(buffer,punc,space)
 *
 * Descrip: No Description
 *
 * Arg:        buffer [UNKN ] Undocumented argument [char *]
 * Arg:          punc [UNKN ] Undocumented argument [char *]
 * Arg:         space [UNKN ] Undocumented argument [char *]
 *
 */
# line 118 "wisestring.dy"
void compress_space_around_punc(char * buffer,char * punc,char * space)
{
  char * runner;
  char * run2;

  runner = buffer;
  for(;;) {
    /*** find next position ***/
    for(;*runner && strchr(punc,*runner) == NULL;runner++ )
      ;

    if( *runner == '\0')
      break;

    /*** this is a puncation position      ***/
    /*** go back over any number of spaces ***/
      
    for(run2 = runner - 1;run2 >= buffer && strchr(space,*run2) != NULL;run2--)
      ;

    /*** run2 now at first non space position after runner ***/

    if( run2 != buffer && run2+1 < runner ) {
      /*** we have some spaces to remove ***/

      strcpy(run2+1,runner);
      runner = run2+1; /** runner now at first punc position **/
    }

    for(;*runner && strchr(punc,*runner) != NULL;runner++)
      ;
    if( *runner == '\0')
      break;

    run2 = runner; /*** run2 at position after last punc ***/

    for(;*runner && strchr(space,*runner) != NULL;runner++)
      ;
    if( run2 != runner ) {
      strcpy(run2,runner);
    }

    runner = run2;
  }

}
      
/* Function:  striptoprint(line)
 *
 * Descrip:    useful strip functions to remove nasty chars
 *             does not allocate memory, simply uses       
 *             given memory, but returns the line pointer
 *             so you can use it in nested function calls                                
 *
 *
 * Arg:        line [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [char *]
 *
 */
# line 171 "wisestring.dy"
char * striptoprint(char * line)
{
  char * run;
  char * base;

  if( line == NULL )
    return NULL;

  base = run = line;

  for(;*line;line++)
    if( isprint(*line) )
      *(run++) = *line;

  *run = '\0';

  return base;
}

/* Function:  stringalloc_next_quoted_string(buffer)
 *
 * Descrip:    takes str's of type <garbage> "xxxxx"   
 *             and gives back xxxxxx .             
 *             stringalloc'd piece so make sure you free it
 *
 *
 * Arg:        buffer [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [char *]
 *
 */
# line 195 "wisestring.dy"
char * stringalloc_next_quoted_string(char * buffer)
{
  char c;
  char * base;
  char * ret;

  for(;*buffer && *buffer != '"';buffer++)
    ;
  if( ! *buffer )
    return NULL;

  for(base=buffer+1,buffer++;*buffer && *buffer != '"';buffer++)
    ;
  if( !*buffer )
    return NULL;

  c = *buffer;
  *buffer = '\0';
	
  ret = stringalloc (base);
  *buffer = c;

  return ret;
}


/* Function:  strwhitestartcmp(line,str,whitespace)
 *
 * Descrip:    sees if line starts with str, ignoring whitespace
 *
 *             returns 0 if they match, to look like strcmp
 *
 *
 * Arg:              line [UNKN ] Undocumented argument [char *]
 * Arg:               str [UNKN ] Undocumented argument [char *]
 * Arg:        whitespace [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 226 "wisestring.dy"
int strwhitestartcmp(char * line,char * str,char * whitespace)
{
  while( strchr(whitespace,*line) != NULL )
    line++;

  if( *line == '\0')
    return -(int)(*str);

  for(;*line == *str && *line && *str;line++,str++)
    ;
	
  if( *str == '\0')
    return 0;

  else return (int)(*line - *str);
}

/* Function:  strwordcmp(buf,str,space)
 *
 * Descrip:    sees if buf matches str\s in perl regex, ie
 *             a word match
 *
 *             space defined \s
 *
 *             returns 0 if they match to look like strcmp
 *
 *
 * Arg:          buf [UNKN ] Undocumented argument [char *]
 * Arg:          str [UNKN ] Undocumented argument [char *]
 * Arg:        space [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 251 "wisestring.dy"
int strwordcmp(char * buf,char * str,char * space)
{
  for(;*buf && *str && *buf == *str;buf++,str++)
    ;

  /** end of buf... get out, buf not space, get out, not end of str **/

  if( *str != '\0' || strchr(space,*buf) == NULL ) 
    return 1;

  return 0;
}

/* Function:  strstartcmp(buf,str)
 *
 * Descrip:    sees if buf starts with str.
 *
 *             returns 0 if so, to mimic strcmp
 *
 *
 * Arg:        buf [UNKN ] Undocumented argument [char *]
 * Arg:        str [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 269 "wisestring.dy"
int strstartcmp(char * buf,char * str)
{
  for(;*buf && *str;str++,buf++)
    {
      if( *buf > *str)
	return 1;
      else if( *buf < *str)
	return -1;
      else	continue;
    }
  if( !(*str) )
    return 0;

  return 1;
}

/* Function:  print_numbered_line(num,ofp)
 *
 * Descrip:    prints lines like _1_________
 *
 *
 * Arg:        num [UNKN ] Undocumented argument [int]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 288 "wisestring.dy"
void print_numbered_line(int num,FILE * ofp)
{
  register int i;
  fprintf(ofp,"_%d",num);

  if( num < 10000 )
    i=4;
  if( num < 1000 )
    i=3;
  if( num < 100 )
    i=2;
  if( num < 10 )
    i=1;
  for(;i<68;i++)
    fputc('_',ofp);
  fputc('\n',ofp);
  return;
}

/* Function:  print_line(ofp)
 *
 * Descrip:    prints _______________ (70 chars)
 *
 *
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 310 "wisestring.dy"
void print_line(FILE * ofp)
{
  register int i;
  for(i=0;i<70;i++)
    fputc('_',ofp);
  fputc('\n',ofp);
  return;
}

/* Function:  only_whitespace(str,space)
 *
 * Descrip:    returns true if str is only made from space.
 *
 *             Deprecated
 *
 *
 * Arg:          str [UNKN ] Undocumented argument [char *]
 * Arg:        space [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 325 "wisestring.dy"
boolean only_whitespace(char * str,char * space)
{
  for(;*str && strchr(space,*str) != NULL;str++)
    ;
  if( *str != '\0')
    return FALSE;
  return TRUE;
}

/* Function:  chop_newline(str)
 *
 * Descrip:    removes trailing newline if present
 *
 *
 * Arg:        str [UNKN ] Undocumented argument [char *]
 *
 */
# line 337 "wisestring.dy"
void chop_newline(char * str)
{
  str += strlen(str) -1;
  if( *str == '\n' ) {
    *str = '\0';
  }
}


/* Function:  good_datastring_fromend(str)
 *
 * Descrip:    Tries to find the last 'database name' 
 *             type string from a string.
 *
 *             Does not allocate memory
 *
 *
 * Arg:        str [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [char *]
 *
 */
# line 352 "wisestring.dy"
char * good_datastring_fromend(char * str)
{
  register char * base=str;
	
  str+=strlen(str);
  for(str--;str >= base;str--)
    if( !isalnum(*str) && *str != '_' && *str != '.')
      break;
  return str == base ? base : ++str;
}

/* Function:  looks_like_vms(str)
 *
 * Descrip:    not useful
 *
 *
 * Arg:        str [UNKN ] Undocumented argument [const char *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 367 "wisestring.dy"
boolean looks_like_vms(const char * str)
{
  if ( strchr(str,':') != NULL)
    return TRUE;
  else	return FALSE;
}

/* Function:  looks_like_unix(str)
 *
 * Descrip:    not useful
 *
 *
 * Arg:        str [UNKN ] Undocumented argument [const char *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 378 "wisestring.dy"
boolean looks_like_unix(const char * str)
{
  if( strchr(str,'/') != NULL)
    return TRUE;
  else	return FALSE;
}

/* Function:  estrcasecmp(one,two)
 *
 * Descrip:    returns strcmp on the captilalised
 *             one and two bufferers (doesn't touch them).
 *
 *
 * Arg:        one [UNKN ] Undocumented argument [char *]
 * Arg:        two [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 389 "wisestring.dy"
int estrcasecmp(char *  one,char *  two)
{
  char  *  tempone;
  char  *  temptwo;
  int out;

  tempone=stringalloc(one);
  temptwo=stringalloc(two);

  capitalise(tempone);
  capitalise(temptwo);

  out=strcmp(tempone,temptwo);

  free(tempone);
  free(temptwo);

  return  out;
}

/* Function:  strend(bu,se)
 *
 * Descrip:    useless
 *
 *
 * Arg:        bu [UNKN ] Undocumented argument [char *]
 * Arg:        se [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [char *]
 *
 */
# line 413 "wisestring.dy"
char * strend(char * bu,char * se)
{
  register char * runner;

  runner=strstr(bu,se);
  if(runner == NULL)
    return NULL;

  runner+=strlen(se);

  return runner;
}


/* Function:  string_before_equality(string)
 *
 * Descrip:    Useless
 *
 *
 * Arg:        string [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [char *]
 *
 */
# line 431 "wisestring.dy"
char * string_before_equality(char * string)
{
  char * base;
  char * runner;

  runner=base=stringalloc(string);
	
  for(;*runner && *runner != '=';runner++)
    ;
  if(*runner == '\0')
    {
      free(base);
      return NULL;
    }

  for(runner--;isspace(*runner);runner--)
    ;
  runner++;
  *runner='\0';

  runner=stringalloc(base);
  free(base);

  return runner;
}

/* Function:  number_from_quoted_equality(s)
 *
 * Descrip:    supposedly gets a number from a
 *             string like xxx="12"
 *
 *
 * Arg:        s [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 461 "wisestring.dy"
int number_from_quoted_equality(char * s)
{
  char * st;
  int ret;
  char c;

  for(;*s && *s != '=';s++) 
    ;
  
  if( !*s ) 
    return 0;
  
  for(;*s && isspace(*s);s++)
    ;
  /*** start of the number now... ***/ 

  st = s;
  for(;*s && isalnum(*s);s++)
    ;
  c = *s;
  *s = '\0';
  ret = atoi(st);
  *s = c;

  return ret;
}


/* Function:  number_from_equality(string)
 *
 * Descrip:    supposedly gets a number from
 *             xxxx=12
 *
 *
 * Arg:        string [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 493 "wisestring.dy"
int number_from_equality(char * string)
{
  char * st;

  for(;*string && *string != '=';string++)
    ;
	
  if(*string == '\0')
    return 0;
	
  for(;isspace(*string) && *string;string++)
    ;
  st=string;

  for(;!isspace(*string) && *string;string++)
	
    *string='\0';

  return atoi(st);
}

/* Function:  sub_string(into,key,sub)
 *
 * Descrip:    deprecated
 *
 *
 * Arg:        into [UNKN ] Undocumented argument [char *]
 * Arg:         key [UNKN ] Undocumented argument [char *]
 * Arg:         sub [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [char *]
 *
 */
# line 518 "wisestring.dy"
char * sub_string(char * into,char * key,char * sub)
{
  char buffer[MAXLINE];
  char * runner;
  char * stop;

  runner=buffer;

  while( (stop=strstr(into,key)) != NULL)
    {
      while( into != stop )
	*(runner++)=(*(into++));
      *runner='\0';
      strcat(buffer,sub);
      runner=buffer+strlen(buffer);
      into+=strlen(key);
    }
  while( *into )
    *(runner++)=(*(into++));
  *(runner)='\0';

  return stringalloc(buffer);
}

/* Function:  string_from_quoted_equality(string)
 *
 * Descrip:    gets the string from xxx="yyy". Returns
 *             yyy allocated, and messes around with string
 *
 *
 * Arg:        string [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [char *]
 *
 */
# line 546 "wisestring.dy"
char * string_from_quoted_equality(char * string)
{
  return string_from_charred_equality(string,'"');
}

/* Function:  string_from_charred_equality(string,quote)
 *
 * Descrip:    gets the string from xxx="yyy" ,where " comes
 *             from the quote argument Returns
 *             yyy allocated, and messes around with string
 *
 *
 * Arg:        string [UNKN ] Undocumented argument [char *]
 * Arg:         quote [UNKN ] Undocumented argument [char]
 *
 * Return [UNKN ]  Undocumented return value [char *]
 *
 */
# line 556 "wisestring.dy"
char * string_from_charred_equality(char * string,char quote)
{
  char * base;

  for(;*string && *string != '=';string++)
    ;
	
  if(*string == '\0')
    return NULL;
	
  for(;*string && *string != quote;string++)
    ;
	
  if(*string == '\0')
    return NULL;

  string++;
  base=string;
	
  for(;*string && *string != quote;string++)
    ;

  *string='\0';

  return stringalloc(base);
}	

/* Function:  breakstring(string,parsestr)
 *
 * Descrip:    A call to /breakstring_protect(string,parsestr,"\"")
 *
 *
 * Arg:          string [UNKN ] Undocumented argument [char *]
 * Arg:        parsestr [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [char **]
 *
 */
# line 586 "wisestring.dy"
char ** breakstring(char * string,char * parsestr)
{
  return breakstring_protect(string,parsestr,"\"");
}

/* Function:  breakstring_protect(string,parsestr,strpair)
 *
 * Descrip:    will parse out words in string using parse  
 *             as white space, like strtok, but strings    
 *             enclosed in characters from strpair are not 
 *             taken in parsed form.                       
 *             breakstring =                               
 *             breakstring_protect(string,parse,"\"");     
 *             hence will not break in double quotes       
 *                                                          
 *             unlike strtok they return char **           
 *             which is a list of char * of words          
 *             the last being NULL'd                       
 *                                                          
 *                                                          
 *             They returned an alloc'd char ** which you  
 *             are expected to free. Standard idiom is     
 *             base=brk=breakstring(buffer,spacestr)       
 *               (NB - spacestr #defin'd above )           
 *               ... do stuff using brk                     
 *               eg *brk = first word                       
 *                  *(++brk) = next word                    
 *                last word = NULL                          
 *              cleanup by ckfree(base)                     
 *
 *
 * Arg:          string [UNKN ] Undocumented argument [char *]
 * Arg:        parsestr [UNKN ] Undocumented argument [char *]
 * Arg:         strpair [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [char **]
 *
 */
# line 615 "wisestring.dy"
char ** breakstring_protect(char * string,char * parsestr,char * strpair)
{
  char ** outstr;
  int index=1;

  outstr=(char **) ckcalloc (128,sizeof(char *));

	
  while( strchr(parsestr,*string) != NULL)
    string++;

  outstr[0]=string;

  while(*string)
    {
      if( strchr(strpair,*string) != NULL )
	{
	  auto char c;
	  c = *string;
	  for(string++;strchr(strpair,*string) == NULL && *string;string++)
	    if( *string == '\\' )
	      string++;
	  if( *string == '\0')
	    {
	      log_full_error(WARNING,0,"In breakstring_protect, reached endofline in protected string [%s]",outstr[index]);
	      outstr[++index]=NULL;
	      return outstr;
	    }
	  string++;		/* shifts on one past " */
	}
      else if( strchr(parsestr,*string) != NULL)
	{
	  while( *string && strchr(parsestr,*string) != NULL)
	    *(string++)='\0';
	  if( *string == '\0')
	    break;
	  string--;
	  *string='\0';
	  outstr[index++]=string+1;
	  string++;
	}
      else string++;
    }


  outstr[index]=NULL;
	


  return outstr;
}

/* Function:  strip_quote_chars(string,quote)
 *
 * Descrip:    removes chars in quote
 *
 *
 * Arg:        string [UNKN ] Undocumented argument [char *]
 * Arg:         quote [UNKN ] Undocumented argument [char *]
 *
 */
# line 670 "wisestring.dy"
void  strip_quote_chars(char * string,char * quote)
{
  char * base =string;
  char * run;

  if( strchr(quote,*string) != NULL )
    {
      for(string++;*string;string++)
	*(string-1)= (*string);
      *(string-1) = '\0';
    }


  run = base+strlen(base)-1;
  if( strchr(quote,*run) != NULL )
    *run= '\0';

}

/* Function:  padstring(buffer,string,maxlen)
 *
 * Descrip:    copies string into buffer, and if under maxlen,
 *             adds spaces. Does *not* put in '\0'
 *
 *
 * Arg:        buffer [UNKN ] Undocumented argument [char *]
 * Arg:        string [UNKN ] Undocumented argument [char *]
 * Arg:        maxlen [UNKN ] Undocumented argument [int]
 *
 */
# line 693 "wisestring.dy"
void padstring(char * buffer,char * string,int maxlen)
{

  for(;*string && maxlen > 0;string++,buffer++,maxlen--)
    *buffer=(*string);
  if( maxlen > 0)
    for(;maxlen > 0;maxlen--,buffer++)
      *buffer=' ';
  return;
}

/* Function:  capitalise(word)
 *
 * Descrip:    toupper's each char in word
 *
 *
 * Arg:        word [UNKN ] Undocumented argument [char *]
 *
 */
# line 707 "wisestring.dy"
void capitalise(char * word)
{
  for(;*word;word++)
    *word=toupper(*word);
}



/* Function:  show_line(line,max,*ofp)
 *
 * Descrip:    This shouws line putting a new line in every max
 *             chars, not minding word boundaries
 *
 *
 * Arg:        line [UNKN ] Undocumented argument [char *]
 * Arg:         max [UNKN ] Undocumented argument [int]
 * Arg:        *ofp [UNKN ] Undocumented argument [FILE]
 *
 */
# line 719 "wisestring.dy"
void  show_line(char * line,int max,FILE *ofp)
{
  int i;

  for(i=0;*(line+i);i++)
    {
      if( (i%max) == 0 && i != 0)
	fputc('\n',ofp);
      fputc(*(line+i),ofp);
    }
  fputc('\n',ofp);

  return;
}

/* Function:  show_text(line,max,*ofp)
 *
 * Descrip:    This shouws line putting a new line in every max
 *             chars, *minding* word boundaries
 *
 *
 * Arg:        line [UNKN ] Undocumented argument [char *]
 * Arg:         max [UNKN ] Undocumented argument [int]
 * Arg:        *ofp [UNKN ] Undocumented argument [FILE]
 *
 */
# line 738 "wisestring.dy"
void show_text(char * line,int max,FILE *ofp)
{
  char push;
  char *pushpoi;
  char * runner;

  for(runner=line;;)
    {
      for(;(int)(runner-line) < max && *runner;)
	{
	  pushpoi=runner;
	  while( *(++runner) != ' ' && *runner)
	    ;			/* NB dodgy pre ++ in above line */
	}
      if(*runner == '\0')
	{
	  if((int)(runner-line) < max)
	    fprintf(ofp,"%s\n",line);
	  else	{
		  push=(*pushpoi);
		  *pushpoi='\0';
		  fprintf(ofp,"%s\n",line);
		  *pushpoi=push;
		  line=pushpoi+1;
		  fprintf(ofp,"%s\n",line);
		}
	  break;
	}
      else 	{
		  push=(*pushpoi);
		  *pushpoi='\0';
		  fprintf(ofp,"%s\n",line);
		  *pushpoi=push;
		  line=pushpoi+1;
		  runner=line;
		}
    }

  return;
}

/* Function:  second_word_alloc(str,space)
 *
 * Descrip:    returns the second word alloc'd for
 *
 *             xxx yyyy
 *
 *             returns yyyy.
 *
 *
 * Arg:          str [UNKN ] Undocumented argument [char *]
 * Arg:        space [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [char *]
 *
 */
# line 786 "wisestring.dy"
char * second_word_alloc(char * str,char * space)
{
  char * out;
  char c;

  for(;*str && strchr(space,*str) == NULL;str++)
    ;

  if( *str == '\0') {
    warn("Can only find one word in [%s] - certainly can't alloc the second",str);
    return NULL;
  }
  
  for(str++;*str && strchr(space,*str) != NULL;str++)
    ;
  
  out=str;

  for(;*str && strchr(space,*str) == NULL;str++)
    ;

  c = *str;
  *str ='\0';
  out = stringalloc(out);

  *str= c;

  return out;
}
      
/* Function:  stringallocf(str,)
 *
 * Descrip:    Don't use this
 *
 *             sprintf's then allocs
 *
 *
 * Arg:        str [UNKN ] Undocumented argument [char *]
 * Arg:            [UNKN ] Undocumented argument [.]
 *
 * Return [UNKN ]  Undocumented return value [char *]
 *
 */
# line 821 "wisestring.dy"
char * stringallocf(char * str,...)
{
  char buffer[1024];
  va_list ap;

  va_start(ap,str);

  vsprintf(buffer,str,ap);

  va_end(ap);

  return stringalloc(buffer);
}
  
  
/* Function:  stringalloc(c)
 *
 * Descrip:    returns the allocated copy of c.
 *             Usually called strdup in other packages
 *
 *
 * Arg:        c [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [char *]
 *
 */
# line 840 "wisestring.dy"
char * stringalloc(char * c)
{
  register int i;
  register char * temp;

  if(c == NULL) {
    warn("Passed stringalloc a NULL pointer");
    return NULL;
  }

  i = strlen(c);

  if((temp=(char *)ckalloc(sizeof(char)*(i+1))) == NULL) {
    warn("stringalloc unable to allocate memory");
    return NULL;
  }

  strcpy(temp,c);

  return temp;
}

# line 1036 "wisestring.c"

#ifdef _cplusplus
}
#endif
