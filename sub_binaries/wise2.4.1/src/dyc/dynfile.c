#ifdef _cplusplus
extern "C" {
#endif
#include "dynfile.h"
#include "linesubs.h"

/* Function:  make_pdoc_output(dfp,mod_name,ofp)
 *
 * Descrip:    makes a pdoc for this file.
 *
 *
 *
 * Arg:             dfp [UNKN ] DYNFILE object [DYNFILE *]
 * Arg:        mod_name [UNKN ] Undocumented argument [char *]
 * Arg:             ofp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 56 "dynfile.dy"
boolean make_pdoc_output(DYNFILE * dfp,char * mod_name,FILE * ofp)
{
  /*int i;*/

  /*** first thing is to add the module name as a link to this file ***/

  fprintf(ofp,"&%s\n",mod_name);

  /*** at the moment we don't handle module info ***/

  /*** ok - now write out a head-list for functions ***/

  return TRUE;
}
  

/* Function:  have_got_unknown_func_DYNFILE(dfp,no,should_complain)
 *
 * Descrip:    reports back whether there are any undocumented
 *             file. if should_complain is TRUE, will issue
 *             warnings through warn. Writes the number of
 *             undocumented functions into no
 *
 *
 * Arg:                    dfp [UNKN ] DYNFILE object [DYNFILE *]
 * Arg:                     no [WRITE] pointer to some memory for the number of undoc'd funcs [int *]
 * Arg:        should_complain [UNKN ] if TRUE, will issue warn statements [boolean]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 82 "dynfile.dy"
boolean have_got_unknown_func_DYNFILE(DYNFILE * dfp,int * no,boolean should_complain)
{
  int i;
  boolean ret = FALSE;
  

  *no =0;

  for(i=0;i<dfp->len;i++) 
    if( dfp->info[i]->functype == FI_UNKNOWN) {
      ret = TRUE;
      (*no)++;
      if( should_complain == TRUE) {
	warn("You have not documented %s",dfp->info[i]->name);
      }
    }

  return ret;
}

/* Function:  write_Dynamite_minimal_func(dfp)
 *
 * Descrip:    writes basically just the #include "self.h" line
 *             in the dynamite file
 *
 *
 * Arg:        dfp [UNKN ] Undocumented argument [DYNFILE *]
 *
 */
# line 106 "dynfile.dy"
void write_Dynamite_minimal_func(DYNFILE * dfp)
{
  fprintf(dfp->func,"#include \"%s.h\"\n\n",dfp->sourceroot);
  dfp->funcpos+=2;

  return;
}

/* Function:  FuncInfo_from_name_DYNFILE(dfp,name)
 *
 * Descrip:    Finds a funcinfo with name. returns NULL
 *             if it cant find them
 *
 *
 * Arg:         dfp [UNKN ] Undocumented argument [DYNFILE *]
 * Arg:        name [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [FuncInfo *]
 *
 */
# line 118 "dynfile.dy"
FuncInfo * FuncInfo_from_name_DYNFILE(DYNFILE * dfp,char * name)
{
  int i;

  for(i=0;i<dfp->len;i++) 
    if( strcmp(dfp->info[i]->name,name) == 0 ) 
      return dfp->info[i];

  return NULL;
}

/* Function:  show_html_DYNFILE(dfp,ofp,type)
 *
 * Descrip:    Deprecated - use pdoc
 *
 *
 * Arg:         dfp [UNKN ] Undocumented argument [DYNFILE *]
 * Arg:         ofp [UNKN ] Undocumented argument [FILE *]
 * Arg:        type [UNKN ] Undocumented argument [FuncOrderType]
 *
 */
# line 133 "dynfile.dy"
void show_html_DYNFILE(DYNFILE * dfp,FILE * ofp,FuncOrderType type)
{
  int i;

  sort_DYNFILE_FuncOrderType(dfp,type);
  
  if( is_separated_types(type) == TRUE) {

    fprintf(ofp,"Externally called functions functions\n<ul>");
    for(i=0;i<dfp->len;i++) {
      if( dfp->info[i]->functype != FI_CALLABLE) {
	break;
      }
      fprintf(ofp,"<li><a href=\"#%s\"> %s </a>\n",dfp->info[i]->name,dfp->info[i]->complete_name);
    }
    fprintf(ofp,"</ul><p>Unplaced (no information about calling) functions\n<ul>");
    for(;i<dfp->len;i++) {
      if( dfp->info[i]->functype != FI_UNKNOWN) {
	break;
      }
      fprintf(ofp,"<li><a href=\"#%s\"> %s </a>\n",dfp->info[i]->name,dfp->info[i]->complete_name);
    }
    fprintf(ofp,"</ul><p>Internal (unlikely to be called directly) functions\n<ul>");
    for(;i<dfp->len;i++) {
      fprintf(ofp,"<li><a href=\"#%s\"> %s </a>\n",dfp->info[i]->name,dfp->info[i]->complete_name);
    }
    
    fprintf(ofp,"</ul><hr>Externally called functions<p>\n");
    for(i=0;i<dfp->len;i++) {
      if( dfp->info[i]->functype != FI_CALLABLE) {
	break;
      }
      fprintf(ofp,"<a name=\"%s\"> %s </a><p>Description<pre>\n",dfp->info[i]->name,dfp->info[i]->complete_name);
      dump_Ftext(dfp->info[i]->ft,ofp);
      fprintf(ofp,"</pre>\n");
    }
    
    fprintf(ofp,"<hr>Unplaced functions<p>\n");
    for(;i<dfp->len;i++) {
      if( dfp->info[i]->functype != FI_UNKNOWN) {
	break;
      }
      fprintf(ofp,"<a name=\"%s\"> %s </a><p>Description<pre>\n",dfp->info[i]->name,dfp->info[i]->complete_name);
      dump_Ftext(dfp->info[i]->ft,ofp);
      fprintf(ofp,"</pre>\n");
    }
    fprintf(ofp,"<hr>Internal functions<p>\n");
    for(;i<dfp->len;i++) {
      fprintf(ofp,"<a name=\"%s\"> %s </a><p>Description<pre>\n",dfp->info[i]->name,dfp->info[i]->complete_name);
      dump_Ftext(dfp->info[i]->ft,ofp);
      fprintf(ofp,"</pre>\n");
    }
  } else {
    warn("Can't currently handle no CUI sorted Html functions, yikes");
  }

}

/* Function:  place_func_decl(dfp,fi)
 *
 * Descrip:    Places the function declaration potential with
 *             package protection
 *
 *
 * Arg:        dfp [UNKN ] Undocumented argument [DYNFILE *]
 * Arg:         fi [UNKN ] Undocumented argument [FuncInfo *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 196 "dynfile.dy"
int place_func_decl(DYNFILE * dfp,FuncInfo * fi)
{

  fprintf(dfp->head,"%s ",fi->ret->type);
  
  fprintf(dfp->head,"%s%s;\n",dfp->package_name == NULL ? "" : dfp->package_name,fi->stripped_return);
  if( dfp->package_name != NULL ) {
    fprintf(dfp->head,"#define %s %s%s\n",fi->name,dfp->package_name,fi->name);
    return 2;
  }
  else {
    return 1;
  }
}


/* Function:  place_declarations_DYNFILE(dfp,type)
 *
 * Descrip:    writes the header declarations of functions
 *             into the header file, with documentation.
 *
 *             The FuncOrderType is either the more common,
 *             usual order or alphabetical
 *
 *
 * Arg:         dfp [UNKN ] Undocumented argument [DYNFILE *]
 * Arg:        type [UNKN ] Undocumented argument [FuncOrderType]
 *
 */
# line 219 "dynfile.dy"
void place_declarations_DYNFILE(DYNFILE * dfp,FuncOrderType type)
{
  int i;

  sort_DYNFILE_FuncOrderType(dfp,type);

  if( is_separated_types(type) == TRUE) {
    fprintf(dfp->head,"\n\n    /***************************************************/\n");
    fprintf(dfp->head,    "    /* Callable functions                              */\n");
    fprintf(dfp->head,    "    /* These are the functions you are expected to use */\n");
    fprintf(dfp->head,    "    /***************************************************/\n\n");
    for(i=0;i<dfp->len;i++) {
      if( dfp->info[i]->functype != FI_CALLABLE) {
	break;
      }
      fprintf(dfp->head,"\n\n");
      show_eddystyle_FuncInfo(dfp->info[i],dfp->head);
      place_func_decl(dfp,dfp->info[i]);

    }
    
    fprintf(dfp->head,"\n\n  /* Unplaced functions */\n");
    fprintf(dfp->head,"  /* There has been no indication of the use of these functions */\n");
    for(;i<dfp->len;i++) {
      if( dfp->info[i]->functype != FI_UNKNOWN) {
	break;
      }

      place_func_decl(dfp,dfp->info[i]);

    }
    
    fprintf(dfp->head,"\n\n    /***************************************************/\n");
    fprintf(dfp->head,    "    /* Internal functions                              */\n");
    fprintf(dfp->head,    "    /* you are not expected to have to call these      */\n");
    fprintf(dfp->head,    "    /***************************************************/\n");
    for(;i<dfp->len;i++) {
      place_func_decl(dfp,dfp->info[i]);
    }
    
  } /*** not separated types ***/
  else {
    for(i=0;i<dfp->len;i++) {
      fprintf(dfp->head,"\n\n");
      show_eddystyle_FuncInfo(dfp->info[i],dfp->head);

      place_func_decl(dfp,dfp->info[i]);

    }
  }

 
  
}

/* Function:  sort_DYNFILE_FuncOrderType(dfp,type)
 *
 * Descrip:    really a sub-routine. Sorts by the FuncOrderType
 *             which can be alphabetical or by original position
 *
 *
 * Arg:         dfp [UNKN ] Undocumented argument [DYNFILE *]
 * Arg:        type [UNKN ] Undocumented argument [FuncOrderType]
 *
 */
# line 279 "dynfile.dy"
void sort_DYNFILE_FuncOrderType(DYNFILE * dfp,FuncOrderType type)
{
  switch (type) {
  case FO_CUI_ALPHABETICAL :
    sort_DYNFILE(dfp,compare_two_FuncInfo_alpha);
    break;
  case FO_CUI_POSITION :
    sort_DYNFILE(dfp,compare_two_FuncInfo_pos);
    break;
  default :
    warn("Got an unnw FuncOrderType %d",type);
  }

  return;


}

/* Function:  is_separated_types(type)
 *
 * Descrip:    Tells whether this sorting will separate
 *             things into Callable/unknown/internal. All
 *             sorts currently do
 *
 *
 * Arg:        type [UNKN ] Undocumented argument [FuncOrderType]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 303 "dynfile.dy"
boolean is_separated_types(FuncOrderType type)
{
  if( type == FO_CUI_ALPHABETICAL || type == FO_CUI_POSITION)
    return TRUE;
  return FALSE;
}

/* Function:  compare_two_FuncInfo_alpha(one,two)
 *
 * Descrip:    comparison function for sorting by functype then alphabetical
 *
 *
 * Arg:        one [UNKN ] Undocumented argument [FuncInfo *]
 * Arg:        two [UNKN ] Undocumented argument [FuncInfo *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 314 "dynfile.dy"
int compare_two_FuncInfo_alpha(FuncInfo * one,FuncInfo * two)
{
  if( one->functype != two->functype ) {
    return one->functype - two->functype;
  }
  return strcmp(one->name,two->name);
}


/* Function:  compare_two_FuncInfo_pos(one,two)
 *
 * Descrip:    comparison function for sorting by functype then by file position
 *
 *
 * Arg:        one [UNKN ] Undocumented argument [FuncInfo *]
 * Arg:        two [UNKN ] Undocumented argument [FuncInfo *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 327 "dynfile.dy"
int compare_two_FuncInfo_pos(FuncInfo * one,FuncInfo * two)
{
  if( one->functype != two->functype ) {
    return one->functype - two->functype;
  }
  return one->infopos - two->infopos;
}

/* Function:  positionify_DYNFILE(dfp)
 *
 * Descrip:    places the positions of functions as read in the file
 *             (hopefully) into the funcinfos so can be sorted on
 *             them (and rearranged potentially)
 *
 *
 * Arg:        dfp [UNKN ] Undocumented argument [DYNFILE *]
 *
 */
# line 340 "dynfile.dy"
void positionify_DYNFILE(DYNFILE * dfp)
{
  int i;

  for(i=0;i<dfp->len;i++)
    dfp->info[i]->infopos = i;
}


# line 349 "dynfile.dy"
void show_FuncInfo_DYNFILE(DYNFILE * dfp,FuncInfo * fi)
{
  switch (dfp->desc_type) {
  case DYN_DESC_EDDY : 
    dfp->funcpos += show_eddystyle_FuncInfo(fi,dfp->func);
    break;
  default :
    warn("No func descirption like %d",dfp->desc_type);
    break;
  }

}
    
 /**
 FuncInfo * find_FuncInfo_from_name_end(DYNFILE * dfp,char * name)
 {
  int i;

  for(i=dfp->len-1;i >= 0;i--) 
    if( strcmp(dfp->info[i]->name,name) == 0 )
      return dfp->info[i];

  return NULL;
 }
 **/


  /*** really for output ***/

/* Function:  strinopen(target,probe)
 *
 * Descrip:    complex function, better written with
 *             regex's. Tests whether the probe is actually
 *             a in a "open" string, ie. complete string
 *             outside of " or other things. 
 *
 *             this is mainly for use of finding if's and for's
 *             in C statements
 *
 *             it is a horrible kludge, and hopefully now not used
 *
 *
 * Arg:        target [UNKN ] Undocumented argument [char *]
 * Arg:         probe [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [char *]
 *
 */
# line 390 "dynfile.dy"
char * strinopen(char * target,char * probe)
{
  register char * p;
  register char * t;	
  register char * start;
  
  start = target;
  
  for(;*target;target++)
    {
      if( *target == '"')
	{
	  for(target++;*target && *target != '"';target++)
	    ;
	  if( *target != '"' )
	    break;
	  target++;
	}
      if( *target == *probe )
	{
	  if( target != start )
	    if( isalpha(*(target-1)) )
	      continue;
	  for(p=probe,t=target;*p && *t;p++,t++)
	    if( *p != *t )
	      break;
	  
	  if( isalpha(*(t)) )
	    continue;
	  
	  if( *p == '\0')
	    return target;
	}
    }
  return NULL;
}


# line 428 "dynfile.dy"
boolean real_is_unfinished_expression(char * str)
{
  char * runner;
  
  runner = str+strlen(str)-1;
  
  while( isspace(*runner) && runner >= str )
    runner--;
  
  if( runner == str)
    return FALSE;
  
  if( strinopen(str,"if") != NULL || strinopen(str,"for") != NULL || strinopen(str,"while") != NULL || strinopen(str,"switch") != NULL || 
     strinopen(str,"case") != NULL || strinopen(str,"default") )
    return FALSE;
  
  if( *str == '#' )
    return FALSE;
  
  if( *runner == ';')
    return FALSE;
  if( (runner=strinopen(str,"else")) != NULL)
    {
      for(;isalpha(*runner) && *runner;runner++)
	;
      for(;!isalnum(*runner) && *runner;runner++) 
	;
      if( *runner == '\0')
	return FALSE;
      else	return TRUE;
    }
  return TRUE;
}



# line 464 "dynfile.dy"
void flush_line(DYNFILE * dfp)
{
  register char * temp;

  if( *dfp->line != '\0') {
      temp = scan_and_replace_line(dfp->line);
      fprintf(dfp->func,"%s \n",temp);
      strcpy(dfp->line,"");
      dfp->funcpos++;
    }  
}

# line 476 "dynfile.dy"
void flush_line_to_header(DYNFILE * dfp)
{
  register char * temp;

  if( *dfp->line != '\0') {
      temp = scan_and_replace_line(dfp->line);
      fprintf(dfp->head,"%s \n",temp);
      strcpy(dfp->line,"");
    }
}



/*** this function is going to "pad" the line with spaces to the next
  comment block level ***/

# line 492 "dynfile.dy"
void pad_line(DYNFILE * dfp)
{
  register int i;
  int len;
  int nlen;

  if(dfp->commentstart == 0 || strlen(dfp->line)+1 > dfp->commentstart) {
    /** ok, have to find the new level **/
    len = strlen(dfp->line);
    nlen = len + (len%4 == 0 ? 0 : 4-len%4);
    for(i=len;i<nlen;i++)
      dfp->line[i] = ' ';
    dfp->line[i]='\0';
  }
  else {
    for(i=strlen(dfp->line);i<dfp->commentstart;i++)
      dfp->line[i] = ' ';
    dfp->line[i]='\0';
  }
}

/** puts in a double break **/

# line 515 "dynfile.dy"
void add_break(DYNFILE * dfp)
{
  flush_line(dfp);
	
  fprintf(dfp->func,"\n\n");
  dfp->commentstart=0;
  dfp->funcpos+=2;

}


# line 526 "dynfile.dy"
void start_function(DYNFILE * dfp,char * str, ...)
{
  FuncInfo * fi;
  char buffer[MAXLINE];
  va_list ap;

  va_start(ap,str);
  vsprintf(buffer,str,ap);
  va_end(ap);

  /** we have no FuncInfo - make a dummy! **/

  fi = FuncInfo_from_str("Undocumented function from internal dynamite code");

  start_function_FuncInfo(fi,dfp,buffer);
  
}

# line 544 "dynfile.dy"
void start_function_FuncInfo(FuncInfo * fi,DYNFILE * dfp,char * str, ... )
{
  char buffer[MAXLINE];
  va_list ap;

  va_start(ap,str);
  vsprintf(buffer,str,ap);
  va_end(ap);

  if( reconcile_FuncInfo_with_funcstr(fi,buffer) == FALSE ) {
    warn("Internal error with function [%s]: Not reconciled with documentation");
    add_break_to_Ftext(fi->ft);
    add_line_to_Ftext(fi->ft,"WARNING: This function's internal documentation was not reconciled with its argument list");
  }

  show_FuncInfo_DYNFILE(dfp,fi);
  fi->line_in_c = dfp->funcpos;
  fi->complete_name = stringalloc(buffer);

  add_DYNFILE(dfp,fi);

  true_start_function(dfp,buffer);
}

# line 568 "dynfile.dy"
void true_start_function(DYNFILE * dfp,char * str, ... )
{
  va_list ap;
  register int i;


  /** put away whatever was there before **/
  flush_line(dfp);
  
  
  for(i=0;i<128;i++)
    dfp->tag[i]=NULL;
  
  va_start(ap,str);
  vsprintf(dfp->line,str,ap);
  va_end(ap);


  if( dfp->bracelevel > 0 ) {
    warn("In dynfile... you are attempting to start a function [%s] inside a separate brace syste. Bad bad news!",dfp->line);
  }

  flush_line(dfp);


  /*** start function ***/
  
  fprintf(dfp->func,"{\n");
  dfp->funcpos++;
  dfp->bracelevel++;
  dfp->infunc = TRUE;

  return;
}

# line 603 "dynfile.dy"
void close_function(DYNFILE * dfp)
{
  if( dfp->infunc == FALSE ){
    warn("Attempting to close a function when you are not even in one... problem surely!");
  }

  flush_line(dfp);
  
  dfp->bracelevel--;
  
  strcpy(dfp->line,"} ");
  dfp->commentstart=0;
  pad_line(dfp);
}



# line 620 "dynfile.dy"
void hang_expr(DYNFILE * dfp,char * str, ...)
{
  char buffer[MAXLINE];
  va_list ap;
  
  flush_line(dfp);
  
  dfp->bracelevel++;
  current_indent(dfp);
  dfp->bracelevel--;
  
  va_start(ap,str);
  vsprintf(buffer,str,ap);
  va_end(ap);
  strcat(dfp->line,buffer);

  if( real_is_unfinished_expression(buffer) )
    strcat(dfp->line,"; ");
  
  pad_line(dfp);
}

# line 642 "dynfile.dy"
void expr(DYNFILE * dfp,char * str, ...)
{
  char buffer[MAXLINE];
  va_list ap;

  flush_line(dfp);
  current_indent(dfp);


  va_start(ap,str);
  vsprintf(buffer,str,ap);
  va_end(ap);
  strcat(dfp->line,buffer);

  if( real_is_unfinished_expression(buffer) )
    strcat(dfp->line,"; ");

  
  pad_line(dfp);
}

# line 663 "dynfile.dy"
void macro(DYNFILE * dfp,char * str, ...)
{
  char buffer[MAXLINE];
  va_list ap;
  
  flush_line(dfp);
  
  va_start(ap,str);
  vsprintf(buffer,str,ap);
  
  if(buffer[0] != '#')
    warn("asking for macro without # as first letter");
  

  /** no ident! **/
  strcat(dfp->line,buffer);
  
  pad_line(dfp);
}


# line 684 "dynfile.dy"
void warn_expr(DYNFILE * dfp,char * str, ...)
{
  char buffer[MAXLINE];
  char buf2[MAXLINE];
  va_list ap;
  
  flush_line(dfp);
  
  current_indent(dfp);
  
  
  va_start(ap,str);
  
  vsprintf(buffer,str,ap);
  
  sprintf(buf2,"warn(\"%s\");",buffer);
  
  strcat(dfp->line,buf2);
  
  pad_line(dfp);
}


# line 707 "dynfile.dy"
void start_struct(DYNFILE * dfp,char * str, ...)
{
  char buffer[MAXLINE];
  va_list ap;
  
  
  flush_line_to_header(dfp);
  
  va_start(ap,str);
  
  vsprintf(buffer,str,ap);
  
  strcat(dfp->line,buffer);
  strcat(dfp->line," { ");
  flush_line_to_header(dfp);
  dfp->bracelevel=1;
}



# line 727 "dynfile.dy"
void close_struct(DYNFILE * dfp,char * str, ...)
{
  char buffer[MAXLINE];
  va_list ap;
  
  flush_line_to_header(dfp);
  
  va_start(ap,str);
  
  vsprintf(buffer,str,ap);
  
  current_indent(dfp);
  
  strcat(dfp->line,"} ");
  strcat(dfp->line,buffer);
  
  pad_line(dfp);
  
  dfp->bracelevel--;
  flush_line_to_header(dfp);

}
	
# line 750 "dynfile.dy"
void struct_macro(DYNFILE * dfp,char * str, ...)
{
  char buffer[MAXLINE];
  va_list ap;

  flush_line_to_header(dfp);

  va_start(ap,str);
  vsprintf(buffer,str,ap);
  va_end(ap);
  strcat(dfp->line,buffer);

  pad_line(dfp);
}

# line 765 "dynfile.dy"
void struct_expr(DYNFILE * dfp,char * str, ...)
{
  char buffer[MAXLINE];
  va_list ap;

  flush_line_to_header(dfp);
  current_indent(dfp);


  va_start(ap,str);
  vsprintf(buffer,str,ap);
  va_end(ap);
  strcat(dfp->line,buffer);

  if( real_is_unfinished_expression(buffer) )
    strcat(dfp->line,"; ");

  
  pad_line(dfp);
}



 /*** Comments ***/

# line 790 "dynfile.dy"
void add_block_comment(DYNFILE * dfp,char * str, ...)
{
  char buffer[MAXLINE];
  va_list ap;
  
  flush_line(dfp);
  
  va_start(ap,str);
  vsprintf(buffer,str,ap);
  va_end(ap);

  current_indent(dfp);
  
  strcat(dfp->line,"/* ");
  strcat(dfp->line,buffer);
  strcat(dfp->line," */");
  
  dfp->commentstart=0;
  
}

# line 811 "dynfile.dy"
void add_end_comment(DYNFILE * dfp,char * str, ...)
{
  char buffer[MAXLINE];
  va_list ap;
  

  va_start(ap,str);
  
  vsprintf(buffer,str,ap);

  strcat(dfp->line,"/* ");
  strcat(dfp->line,buffer);
  strcat(dfp->line," */");

  flush_line(dfp);
}

# line 828 "dynfile.dy"
void add_end_comment_header(DYNFILE * dfp,char * str, ...)
{
  char buffer[MAXLINE];
  va_list ap;
  

  va_start(ap,str);
  
  vsprintf(buffer,str,ap);

  strcat(dfp->line,"/* ");
  strcat(dfp->line,buffer);
  strcat(dfp->line," */");

  flush_line_to_header(dfp);
}
	

# line 846 "dynfile.dy"
void set_commentstart(DYNFILE * dfp,int s)
{
  
  dfp->commentstart=s;
}
	
			
# line 853 "dynfile.dy"
void current_indent(DYNFILE * dfp)
{
  register int i;
  
  dfp->line[0] = '\0';
  
  if( dfp->bracelevel == 0)
    return;	
  
  strcat(dfp->line,"  ");
  for(i=0;i<dfp->bracelevel;i++ ){
    strcat(dfp->line,"  ");
    }
  
}
	
# line 869 "dynfile.dy"
void startcase(DYNFILE * dfp)
{
  flush_line(dfp);
  dfp->bracelevel++;
}

# line 875 "dynfile.dy"
void closecase(DYNFILE * dfp)
{
  flush_line(dfp);
  dfp->bracelevel--;
}




# line 884 "dynfile.dy"
void startbrace_tag(DYNFILE * dfp,char * str)
{
  strcat(dfp->line," { ");
  dfp->bracelevel++;
  strcat(dfp->line,"/*");
  strcat(dfp->line,str);
  strcat(dfp->line,"*/");
  flush_line(dfp);
  if( dfp->tag[dfp->bracelevel] != NULL) {
      warn("Missalgined tags in tag brace. Bad bad bug");
      dfp->tag[dfp->bracelevel] = NULL;
    }
  
  dfp->tag[dfp->bracelevel]=stringalloc(str);
  
}

# line 901 "dynfile.dy"
void startbrace(DYNFILE * dfp)
{
  strcat(dfp->line," { ");
  flush_line(dfp);
  dfp->bracelevel++;
}

# line 908 "dynfile.dy"
void closebrace(DYNFILE * dfp)
{
  flush_line(dfp);
  current_indent(dfp);
  strcat(dfp->line,"} ");
  if( dfp->tag[dfp->bracelevel] != NULL) {
      add_end_comment(dfp,"end of %s",dfp->tag[dfp->bracelevel]);
      dfp->tag[dfp->bracelevel]=ckfree(dfp->tag[dfp->bracelevel]);
    }
  
  flush_line(dfp);
  dfp->bracelevel--;
  dfp->commentstart=0;
}


# line 924 "dynfile.dy"
void fputs_func_DYNFILE(char * str,DYNFILE * dfp)
{
  fputs(str,dfp->func);
  fputc('\n',dfp->func);
  dfp->funcpos++;
}


 /** I/O stuff **/


# line 935 "dynfile.dy"
DYNFILE * open_std_no_html_dynfile(char * name)
{
  return open_std_dynfile(name,NULL);
}

# line 940 "dynfile.dy"
DYNFILE * open_std_dynfile(char * name,char * html_directory)
{
  DYNFILE * dfp;
  char func_buf[512];
  char head_buf[512];
  char html_buf[1024];
  int len;

  if( (len=strlen(name)) > 509 ) {
    warn("I can't believe this, but you are attempting to open a DYNFILE with a name longer than 509 letter. Surely not! [%s]",name);
    return NULL;
  }

  if( html_directory != NULL && strlen(html_directory) + len > 1020 ) {
    warn("I can't believe this, but a combination of the html directory path and the dynamite name is over 1020 letters. Can't handle it [%s][%s]",html_directory,name);
    return NULL;
  }

  
  sprintf(func_buf,"%s.c",name);
  sprintf(head_buf,"%s.h",name);
  if( html_directory != NULL ) 
    sprintf(html_buf,"%s/%s.html",html_directory,name);



  dfp = DYNFILE_alloc_std();

  dfp->func = openfile(func_buf,"W");
  if( dfp->func == NULL ) {
    close_dynfile(dfp);
    warn("Could not open %s as a function file",func_buf);
    return NULL;
  }


  dfp->head = openfile(head_buf,"W");
  if( dfp->head == NULL ) {
    close_dynfile(dfp);
    warn("Could not open %s as a header file",head_buf);
    return NULL;
  }

  if( html_directory != NULL ) {
    dfp->html = openfile(html_buf,"W");
    if( dfp->html == NULL ) {
      close_dynfile(dfp);
      warn("Could not open %s as a html file",html_buf);
      return NULL;
    }
  }

  /*** put in "standard" headers etc ****/

  fprintf(dfp->head,"#ifndef DYNAMITE%sHEADERFILE\n#define DYNAMITE%sHEADERFILE\n",name,name);
  fprintf(dfp->head,"#ifdef _cplusplus\nextern \"C\" {\n#endif\n");
  fprintf(dfp->func,"#ifdef _cplusplus\nextern \"C\" {\n#endif\n");
  
  dfp->sourceroot = stringalloc(name);
 
  /*** put in standard header ***/

  dfp->funcpos = 3; /* function file on its 3rd line */


  return dfp;
}

# line 1008 "dynfile.dy"
void close_dynfile(DYNFILE * dfp)
{
  
  /*** finish C function file ***/
  flush_line(dfp);

  /*** finish C header file ***/

  fprintf(dfp->func,"\n#ifdef _cplusplus\n}\n#endif\n");
  fprintf(dfp->head,"\n#ifdef _cplusplus\n}\n#endif\n");
  fprintf(dfp->head,"\n#endif\n");



  if( dfp->html != NULL)
    fclose(dfp->html);

  if( dfp->head != NULL)
    fclose(dfp->head);

  if( dfp->func != NULL)
    fclose(dfp->func);

  dfp = free_DYNFILE(dfp);
  
}

  
    


# line 1095 "dynfile.c"
/* Function:  swap_DYNFILE(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_DYNFILE
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [FuncInfo **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_DYNFILE(FuncInfo ** list,int i,int j)  
{
    FuncInfo * temp; 
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_DYNFILE(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_DYNFILE which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [FuncInfo **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_DYNFILE(FuncInfo ** list,int left,int right,int (*comp)(FuncInfo * ,FuncInfo * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_DYNFILE(list,left,(left+right)/2);  
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_DYNFILE (list,++last,i);    
      }  
    swap_DYNFILE (list,left,last);   
    qsort_DYNFILE(list,left,last-1,comp);    
    qsort_DYNFILE(list,last+1,right,comp);   
}    


/* Function:  sort_DYNFILE(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_DYNFILE
 *
 *
 * Arg:         obj [UNKN ] Object containing list [DYNFILE *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_DYNFILE(DYNFILE * obj,int (*comp)(FuncInfo *, FuncInfo *)) 
{
    qsort_DYNFILE(obj->info,0,obj->len-1,comp);  
    return;  
}    


/* Function:  expand_DYNFILE(obj,len)
 *
 * Descrip:    Really an internal function for add_DYNFILE
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [DYNFILE *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_DYNFILE(DYNFILE * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_DYNFILE called with no need");    
      return TRUE;   
      }  


    if( (obj->info = (FuncInfo ** ) ckrealloc (obj->info,sizeof(FuncInfo *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_DYNFILE, returning FALSE");  
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_DYNFILE(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [DYNFILE *]
 * Arg:        add [OWNER] Object to add to the list [FuncInfo *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_DYNFILE(DYNFILE * obj,FuncInfo * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_DYNFILE(obj,obj->len + DYNFILELISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->info[obj->len++]=add;   
    return TRUE; 
}    


/* Function:  flush_DYNFILE(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [DYNFILE *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_DYNFILE(DYNFILE * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->info[i] != NULL)  {  
        free_FuncInfo(obj->info[i]); 
        obj->info[i] = NULL; 
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  DYNFILE_alloc_std(void)
 *
 * Descrip:    Equivalent to DYNFILE_alloc_len(DYNFILELISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DYNFILE *]
 *
 */
DYNFILE * DYNFILE_alloc_std(void) 
{
    return DYNFILE_alloc_len(DYNFILELISTLENGTH); 
}    


/* Function:  DYNFILE_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [DYNFILE *]
 *
 */
DYNFILE * DYNFILE_alloc_len(int len) 
{
    DYNFILE * out;  /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = DYNFILE_alloc()) == NULL)  
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->info = (FuncInfo ** ) ckcalloc (len,sizeof(FuncInfo *))) == NULL)   {  
      warn("Warning, ckcalloc failed in DYNFILE_alloc_len"); 
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_DYNFILE(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DYNFILE *]
 *
 * Return [UNKN ]  Undocumented return value [DYNFILE *]
 *
 */
DYNFILE * hard_link_DYNFILE(DYNFILE * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a DYNFILE object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  DYNFILE_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DYNFILE *]
 *
 */
DYNFILE * DYNFILE_alloc(void) 
{
    DYNFILE * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(DYNFILE *) ckalloc (sizeof(DYNFILE))) == NULL)  {  
      warn("DYNFILE_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->desc_type = DYN_DESC_EDDY;  
    /* line[1024] is an array: no default possible */ 
    out->commentstart = 0;   
    out->bracelevel = 0; 
    out->infunc = FALSE; 
    out->sourceroot = NULL;  
    out->info = NULL;    
    out->len = out->maxlen = 0;  
    out->mi = NULL;  
    out->funcpos = 0;    
    out->code_debug_level = 0;   
    out->should_hard_link = FALSE;   
    out->package_name = NULL;    


    return out;  
}    


/* Function:  free_DYNFILE(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DYNFILE *]
 *
 * Return [UNKN ]  Undocumented return value [DYNFILE *]
 *
 */
DYNFILE * free_DYNFILE(DYNFILE * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a DYNFILE obj. Should be trappable");   
      return NULL;   
      }  


#ifdef PTHREAD   
    assert(pthread_mutex_lock(&(obj->dynamite_mutex)) == 0); 
#endif   
    if( obj->dynamite_hard_link > 1)     {  
      return_early = 1;  
      obj->dynamite_hard_link--; 
      }  
#ifdef PTHREAD   
    assert(pthread_mutex_unlock(&(obj->dynamite_mutex)) == 0);   
#endif   
    if( return_early == 1)   
      return NULL;   
    /* obj->*tag[128] is linked in */ 
    if( obj->sourceroot != NULL) 
      ckfree(obj->sourceroot);   
    /* obj->func is linked in */ 
    /* obj->head is linked in */ 
    /* obj->html is linked in */ 
    if( obj->info != NULL)   {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->info[i] != NULL)    
          free_FuncInfo(obj->info[i]);   
        }  
      ckfree(obj->info); 
      }  
    if( obj->mi != NULL) 
      free_ModuleInfo(obj->mi);  
    if( obj->package_name != NULL)   
      ckfree(obj->package_name);     


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
