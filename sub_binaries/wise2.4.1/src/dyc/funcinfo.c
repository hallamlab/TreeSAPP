#ifdef _cplusplus
extern "C" {
#endif
#include "funcinfo.h"


# line 71 "funcinfo.dy"
void dump_FuncInfo(FuncInfo * fi,FILE * ofp ) 
{
  int i;

  fprintf(ofp,"Function info of %s\n",CKS(fi->name));
  dump_Ftext(fi->ft,ofp);

  for(i=0;i<fi->len;i++) {
    fprintf(ofp,"Argument %s type [%s] Desc %s\n",fi->arg[i]->name,fi->arg[i]->type,fi->arg[i]->desc);
  }

}

/* Function:  show_eddystyle_FuncInfo(fi,ofp)
 *
 * Descrip:    shows functions in
 *
 *              *Function:
 *              *
 *              *des
 *              *
 *              * Arg
 *              *
 *
 *             Returns number of lines printed
 *
 *
 * Arg:         fi [UNKN ] Undocumented argument [FuncInfo *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 96 "funcinfo.dy"
int show_eddystyle_FuncInfo(FuncInfo * fi,FILE * ofp) 
{
  int i;
  int len=0;
  int maxn;

  maxn = max_argame(fi );

  fprintf(ofp,"/* Function:  %s(%s",CKS(fi->name),fi->len == 0 ? "void" : fi->arg[0]->name);
  for(i=1;i<fi->len;i++)
    fprintf(ofp,",%s",fi->arg[i]->name);
  fprintf(ofp,")\n *\n");
  len++;

  len += show_eddystyle_Ftext(fi->ft,"Descrip:",15,ofp,"No Description");
  fprintf(ofp," *\n");

  for(i=0;i<fi->len;i++)
    len += show_eddystyle_ArgInfo(fi->arg[i],15,maxn,ofp);

  if( fi->ret != NULL && strcmp(fi->ret->type,"void") != 0 ) {
    fprintf(ofp," *\n * Return [%s]  %s [%s]\n",ArgType_to_string(fi->ret->argtype),fi->ret->desc,CKS(fi->ret->type));
    len +=2;
  }
  len += 2;
  fprintf(ofp," *\n */\n");
  return len;
}

# line 125 "funcinfo.dy"
char * ArgType_to_string(int type)
{
  switch(type) {
  case  ARGTYPE_READ :
    return "READ ";
  case  ARGTYPE_WRITE :
    return "WRITE";
  case   ARGTYPE_READWRITE :
    return "RW   ";
  case ARGTYPE_P2FUNC :
    return "FUNCP";
  case ARGTYPE_OWNER:
    return "OWNER";
  case ARGTYPE_STATIC :
    return "SOFT ";
  default :
    return "UNKN ";
  }
}

# line 145 "funcinfo.dy"
int max_argame(FuncInfo * fi )
{
  int i;
  int max;

  if(fi->len == 0 ) {
    return 0;
  }

  for(i=1,max=strlen(fi->arg[0]->name);i < fi->len;i++) {
    if( max < strlen(fi->arg[i]->name) ) {
      max = strlen(fi->arg[i]->name);
    }
  }

  return max;
}

# line 163 "funcinfo.dy"
int show_eddystyle_ArgInfo(ArgInfo * ai,int depth,int namedepth,FILE * ofp)
{

  fprintf(ofp," * %*s%*s [%s] %s [%s]\n",-(depth-3),"Arg:",namedepth,ai->name,ArgType_to_string(ai->argtype),ai->desc,CKS(ai->type));
  return 1;
}

# line 170 "funcinfo.dy"
void sort_FuncInfo_by_position(FuncInfo * fi)
{
  sort_FuncInfo(fi,compare_ArgInfo_pos);
  return;
}

# line 176 "funcinfo.dy"
int compare_ArgInfo_pos(ArgInfo * one,ArgInfo * two)
{
  return one->argpos- two->argpos;
}

# line 181 "funcinfo.dy"
boolean reconcile_FuncInfo_with_funcstr(FuncInfo * fi,char * pass_str)
{
  /** ok, stupid parsing for the moment **/
  boolean ret = TRUE;
  char * runner;
  char * run2;
  char * arg;
  char ** base;
  char ** brk;
  char * fstr;
  int count;

  fstr = stringalloc(pass_str);

  if( fi->complete_name == NULL ) {
    fi->complete_name = stringalloc(fstr);
  }

  /** get the return type **/

  if( strstartcmp(fstr,"const") == 0 ) {
    runner = fstr+5;
    for(;isspace(*runner);runner++)
      ;
  } else {
    runner = fstr;
  }

  if( strstartcmp(runner,"signed") == 0 ) {
    runner = fstr+6;
    for(;isspace(*runner);runner++)
      ;
  }

  if( strstartcmp(runner,"long") == 0 ) {
    runner = fstr+6;
    for(;isspace(*runner);runner++)
      ;
  } 


  if( strstartcmp(runner,"struct") == 0 ) {
    runner = fstr+6;
    for(;isspace(*runner);runner++)
      ;
  } 

  for(;isalnum(*runner) || *runner == '_' ;runner++)
    ;
  for(run2=runner;isspace(*runner);runner++)
    ;

  if( *runner == '*' ) {
    for(;*runner == '*';runner++)
      ;
  } else {
    runner = run2;
  }
  *runner = '\0';

  if( fi->ret == NULL ) {
    fi->ret = ArgInfo_alloc();

    fi->ret->name = stringalloc("return");
    fi->ret->desc = stringalloc("Undocumented return value");
  } 

  fi->ret->type = stringalloc(fstr);
  

  /** got to the end of the return value **/

  for(run2=runner+1;isspace(*run2);run2++)
    ;
  fi->stripped_return = stringalloc(run2);

  /** get the name of the function ***/
  runner = strchr(run2,'(') ;


  if( runner == NULL ) {
    warn("reconciling function [%s] - no bracket!",pass_str);
    return FALSE;
  }


  arg = runner+1;
 

  *(runner) = '\0';
  
  /** run2 now at the name **/

  if( fi->name != NULL ) {
    if( strcmp(fi->name,run2) != 0 ) {
      warn("In reconciling function different names [%s] [%s]",fi->name,run2);
    }
  }else {
    fi->name = stringalloc(run2);
  }


  /*** now, process the argument list ***/
  runner = arg + strlen(arg) -1;
  
  for(;runner > arg && *runner != ')';runner--)
    ;

  *runner = '\0';


  base = brk = breakstring_protect(arg,",","()");

  if( *brk != NULL && strcmp(*brk,"void") != 0 ) { 
    for(count=0;*brk != NULL;brk++) {
      if( reconcile_FuncInfo_with_argstr(fi,*brk,count++) == FALSE ) 
	ret = FALSE;
    }
  }
  ckfree(base);
  ckfree(fstr);

  sort_FuncInfo_by_position(fi);

  return ret;

} 

# line 309 "funcinfo.dy"
boolean reconcile_FuncInfo_with_pfunc(FuncInfo * fi,char * str,int pos)
{
  char * runner;
  char * name;
  ArgInfo * temp;
  char * held;

  /**
    This is a HUGE kludge. V.v.v. embarrasing

    ***/

  /** assumme type (*name)(type,type,type) ***/

  held = stringalloc(str);

  name = runner = strchr(str,'(');
  
  name = runner = strchr(runner,'*');

  name++;

  for(runner++;!isspace(*runner) && *runner != ')' ;runner++)
    ;

  *runner = '\0';

  if( (temp=get_ArgInfo_by_name(fi,name)) == NULL ) {
    temp = ArgInfo_alloc();
    add_FuncInfo(fi,temp);

    temp->name = stringalloc(name);
    temp->desc = stringalloc("Undocumented argument");
    temp->func_decl = held;
  } else {
    temp->type = stringalloc(str);
    temp->argtype = ARGTYPE_P2FUNC;
    temp->argpos = pos;
    temp->func_decl = held;
  }

  return TRUE;
}
  
# line 353 "funcinfo.dy"
FuncInfo * unknown_user_FuncInfo(char * funcstr)
{
  FuncInfo * fi;

  fi = FuncInfo_from_str("Unknown user-defined function");

  if( reconcile_FuncInfo_with_funcstr(fi,funcstr) == FALSE ) {
    warn("Could not reconcile [%s]... bad internal error.",funcstr);
  }

  return fi;
}
  
# line 366 "funcinfo.dy"
boolean reconcile_FuncInfo_with_argstr(FuncInfo * fi,char * str,int pos) 
{
  char * runner;
  char * name;
  ArgInfo * temp;



  if( strchr(str,'(') != NULL ) 
    return reconcile_FuncInfo_with_pfunc(fi,str,pos);
  
  for(;isspace(*str);str++) 
    ;

  runner = str + strlen(str) -1;
  
  for(;runner > str && isspace(*runner);runner--)
    ;

  *(runner+1) = '\0';

  for(;runner > str && !isspace(*runner);runner--)
    ;
  
  name = runner+1;
  if( strcmp(name,"void") == 0 )
    return TRUE;

  for(;runner > str && isspace(*runner);runner--)
    ;

  *(runner+1) = '\0';
  
  
  if( (temp=get_ArgInfo_by_name(fi,name)) == NULL ) {
    temp = ArgInfo_alloc();
    add_FuncInfo(fi,temp);

    temp->name = stringalloc(name);
    temp->desc = stringalloc("Undocumented argument");
  }

  temp->type = stringalloc(str);
  temp->argpos = pos;

  return TRUE;
}
  
# line 414 "funcinfo.dy"
ArgInfo  * get_ArgInfo_by_name(FuncInfo * fi,char * str)
{
  int i;

  for(i=0;i<fi->len;i++)
    if( strcmp(fi->arg[i]->name,str) == 0 )
      return fi->arg[i];

  return NULL;
}

# line 425 "funcinfo.dy"
ArgInfo * ArgInfo_in_FuncInfo_from_varstr(FuncInfo * fi,char * str,...)
{
  char buffer[MAXLINE];
  ArgInfo * out;
  va_list ap;

  va_start(ap,str);
  vsprintf(buffer,str,ap);
  va_end(ap);

  out = ArgInfo_alloc();
  out->name = stringalloc(buffer);

  add_FuncInfo(fi,out);

  return out;
}

# line 443 "funcinfo.dy"
FuncInfo * FuncInfo_named_from_varstr(int type,char * str, ...)
{
  char buffer[MAXLINE];
  FuncInfo * out;
  va_list ap;

  va_start(ap,str);
  vsprintf(buffer,str,ap);
  va_end(ap);

  out = FuncInfo_alloc_std();
  out->functype = type;
  out->name = stringalloc(buffer);
  out->ft = Ftext_alloc_std();
  add_Ftext(out->ft,Fblock_alloc_std());

  return out;

}

# line 463 "funcinfo.dy"
FuncInfo * FuncInfo_from_str(char * str)
{
  FuncInfo * fi;

  fi = FuncInfo_alloc_std();

  
  fi->ft = single_Ftext_from_str(str);

  return fi;
}

/** I/O from user **/

# line 477 "funcinfo.dy"
ModuleInfo * read_ModuleInfo_line(char * line,FILE * ifp)
{
  ModuleInfo * out;
  char buffer[MAXLINE];

  if( strstartcmp(line,"%module") != 0 ) {
    warn("Attempting to read module help with line starting [%30s] not %module",line);
    return NULL;
  }

  out = ModuleInfo_alloc();

  out->ft = read_Ftext(buffer,MAXLINE,ifp,"%",fgets);

  return out;
}

# line 494 "funcinfo.dy"
FuncInfo * read_FuncInfo_line(char * line,FILE * ifp)
{
  FuncInfo * out;
  ArgInfo * ari;
  char buffer[MAXLINE];
  char * runner;

  
  if( strstartcmp(line,"%func") != 0 ) {
    warn("Attempting to read in-line function help with line starting [%30s] not %func",line);
    return NULL;
  }

  out = FuncInfo_alloc_std();
  out->functype = FI_CALLABLE;
  out->ft = read_Ftext(buffer,MAXLINE,ifp,"%",get_watched_line);

  if( strstartcmp(buffer,"%%") == 0 )
    return out;

  /*** could be in any order ***/

  for(;;) {
    /*    fprintf(stderr,"Looking at [%s]\n",buffer); */
    if( feof(ifp) || ferror(ifp) ) {
      warn("End of file or file read error while in FuncInfo read. Not good!");
      break;
    } else if ( strstartcmp(buffer,"%%") == 0 ) {
      break;
    } else if ( strstartcmp(buffer,"%simple") == 0 ) {
      if( (runner=strtok(buffer+7,spacestr)) == NULL ) {
	warn("Got a simple name specification, but no name!");
      } else {
	out->simple = stringalloc(runner);
      }
      get_watched_line(buffer,MAXLINE,ifp);
    } else if( strstartcmp(buffer,"%arg") == 0) {
      while( get_watched_line(buffer,MAXLINE,ifp) != NULL ) {
	if( buffer[0] == '%' )
	  break;
	ari = read_ArgInfo_line(buffer);

	if( strcmp(ari->name,"return") == 0 ) {
	  out->ret = ari;
	  
	} else 	if( ari != NULL ) {
	  add_FuncInfo(out,ari);
	}
      }
    } else if ( strstartcmp(buffer,"%short") == 0 ) {
      for(runner=buffer;*runner && !isspace(*runner);runner++) 
	;
      for(;*runner && isspace(*runner);runner++)
	;
      out->sdesc=stringalloc(runner);
      get_watched_line(buffer,MAXLINE,ifp);

    } else if ( strstartcmp(buffer,"%type") == 0 ) {
      if( strstr(buffer,"call") != NULL ) {
	out->functype = FI_CALLABLE;
      } else if ( strstr(buffer,"int") != NULL ) {
	out->functype = FI_INTERNAL;
      }
      get_watched_line(buffer,MAXLINE,ifp);
    }

    else {
      warn("Cannot understand %% tag %20s\n",buffer);
      while( get_watched_line(buffer,MAXLINE,ifp) != NULL ) {
	if( buffer[0] == '%' )
	  break;
      }
    } 
    if( buffer[0] == '%')
      continue; /*** back to for(;;) ***/
      
    /* else */
    warn("In funcinfo line, could not understand [%s], going to skip to next %% tag",buffer);

    while( get_watched_line(buffer,MAXLINE,ifp) != NULL ) {
      chop_newline(buffer);
      if( buffer[0] == '%' )
	break;
      else {
	warn("Did not interpret line [%s]\n",buffer);
      }
    }
  }

  return out;
}

# line 586 "funcinfo.dy"
int get_arg_type(char * line,boolean * should_NULL)
{
  char * runner;

  *should_NULL = FALSE;


  for(runner=line;*runner && !isspace(*runner);runner++)
    ;
  runner--;
  if( runner-line < 3 &&  *runner == 'N') {
    *should_NULL = TRUE;
    *runner=' ';
  }
  

  if ( strwordcmp(line,"rw",spacestr) == 0 )
    return ARGTYPE_READWRITE;

  if( strwordcmp(line,"r",spacestr) == 0 ) 
    return ARGTYPE_READ;
  if( strwordcmp(line,"o",spacestr) == 0 ) 
    return ARGTYPE_OWNER;
  if( strwordcmp(line,"s",spacestr) == 0 ) 
    return ARGTYPE_STATIC;
  if( strwordcmp(line,"f",spacestr) == 0 ) 
    return ARGTYPE_P2FUNC;
  else if ( strwordcmp(line,"w",spacestr) == 0 )
    return ARGTYPE_WRITE;
  else return ARGTYPE_UNKNOWN;

}


# line 620 "funcinfo.dy"
ArgInfo * read_ArgInfo_line(char * line)
{
  ArgInfo * out;
  char * runner;
  char * fix;

  out = ArgInfo_alloc();

  for(runner=line;*runner && !isalpha(*runner) ;runner++)
    ;

  fix = runner;

  for(runner=line;*runner && iscword(*runner);runner++)
    ;

  /*** got first word ***/

  *runner = '\0';
  out->name = stringalloc(fix);

  /*** next word ***/
  for(runner++;*runner && !isalpha(*runner);runner++)
    ;

  /*** if it is a valid arg type, get it and move on ***/
 
  if( (out->argtype=get_arg_type(runner,&out->should_NULL)) != ARGTYPE_UNKNOWN) {
    for(;*runner && isalnum(*runner);runner++)
      ;
    for(;*runner && isspace(*runner);runner++)
      ;
  }

  fix = runner;
  for(;*runner && *runner != '\n';runner++)
    ;
  *runner = '\0';

  out->desc = stringalloc(fix);

  return out;
}

# line 624 "funcinfo.c"
/* Function:  hard_link_ArgInfo(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ArgInfo *]
 *
 * Return [UNKN ]  Undocumented return value [ArgInfo *]
 *
 */
ArgInfo * hard_link_ArgInfo(ArgInfo * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a ArgInfo object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  ArgInfo_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ArgInfo *]
 *
 */
ArgInfo * ArgInfo_alloc(void) 
{
    ArgInfo * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(ArgInfo *) ckalloc (sizeof(ArgInfo))) == NULL)  {  
      warn("ArgInfo_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->argtype = ARGTYPE_UNKNOWN;  
    out->should_NULL = FALSE;    
    out->name = NULL;    
    out->type = NULL;    
    out->desc = NULL;    
    out->argpos = 0; 
    out->func_decl = NULL;   


    return out;  
}    


/* Function:  free_ArgInfo(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ArgInfo *]
 *
 * Return [UNKN ]  Undocumented return value [ArgInfo *]
 *
 */
ArgInfo * free_ArgInfo(ArgInfo * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a ArgInfo obj. Should be trappable");   
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
    if( obj->name != NULL)   
      ckfree(obj->name);     
    if( obj->type != NULL)   
      ckfree(obj->type);     
    if( obj->desc != NULL)   
      ckfree(obj->desc);     
    if( obj->func_decl != NULL)  
      ckfree(obj->func_decl);    


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_ErrorInfo(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ErrorInfo *]
 *
 * Return [UNKN ]  Undocumented return value [ErrorInfo *]
 *
 */
ErrorInfo * hard_link_ErrorInfo(ErrorInfo * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a ErrorInfo object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  ErrorInfo_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ErrorInfo *]
 *
 */
ErrorInfo * ErrorInfo_alloc(void) 
{
    ErrorInfo * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(ErrorInfo *) ckalloc (sizeof(ErrorInfo))) == NULL)  {  
      warn("ErrorInfo_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->value = NULL;   
    out->desc = NULL;    


    return out;  
}    


/* Function:  free_ErrorInfo(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ErrorInfo *]
 *
 * Return [UNKN ]  Undocumented return value [ErrorInfo *]
 *
 */
ErrorInfo * free_ErrorInfo(ErrorInfo * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a ErrorInfo obj. Should be trappable"); 
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
    if( obj->value != NULL)  
      ckfree(obj->value);    
    if( obj->desc != NULL)   
      ckfree(obj->desc);     


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_FuncInfo(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_FuncInfo
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [ArgInfo   **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_FuncInfo(ArgInfo   ** list,int i,int j)  
{
    ArgInfo   * temp;    
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_FuncInfo(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_FuncInfo which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [ArgInfo   **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_FuncInfo(ArgInfo   ** list,int left,int right,int (*comp)(ArgInfo   * ,ArgInfo   * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_FuncInfo(list,left,(left+right)/2); 
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_FuncInfo (list,++last,i);   
      }  
    swap_FuncInfo (list,left,last);  
    qsort_FuncInfo(list,left,last-1,comp);   
    qsort_FuncInfo(list,last+1,right,comp);  
}    


/* Function:  sort_FuncInfo(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_FuncInfo
 *
 *
 * Arg:         obj [UNKN ] Object containing list [FuncInfo *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_FuncInfo(FuncInfo * obj,int (*comp)(ArgInfo   *, ArgInfo   *)) 
{
    qsort_FuncInfo(obj->arg,0,obj->len-1,comp);  
    return;  
}    


/* Function:  expand_FuncInfo(obj,len)
 *
 * Descrip:    Really an internal function for add_FuncInfo
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [FuncInfo *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_FuncInfo(FuncInfo * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_FuncInfo called with no need");   
      return TRUE;   
      }  


    if( (obj->arg = (ArgInfo   ** ) ckrealloc (obj->arg,sizeof(ArgInfo   *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_FuncInfo, returning FALSE"); 
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_FuncInfo(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [FuncInfo *]
 * Arg:        add [OWNER] Object to add to the list [ArgInfo   *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_FuncInfo(FuncInfo * obj,ArgInfo   * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_FuncInfo(obj,obj->len + FuncInfoLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->arg[obj->len++]=add;    
    return TRUE; 
}    


/* Function:  flush_FuncInfo(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [FuncInfo *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_FuncInfo(FuncInfo * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->arg[i] != NULL)   {  
        free_ArgInfo(obj->arg[i]);   
        obj->arg[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  swap_err_FuncInfo(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_err_FuncInfo
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [ErrorInfo **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_err_FuncInfo(ErrorInfo ** list,int i,int j)  
{
    ErrorInfo * temp;    
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_err_FuncInfo(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_err_FuncInfo which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [ErrorInfo **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_err_FuncInfo(ErrorInfo ** list,int left,int right,int (*comp)(ErrorInfo * ,ErrorInfo * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_err_FuncInfo(list,left,(left+right)/2); 
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_err_FuncInfo (list,++last,i);   
      }  
    swap_err_FuncInfo (list,left,last);  
    qsort_err_FuncInfo(list,left,last-1,comp);   
    qsort_err_FuncInfo(list,last+1,right,comp);  
}    


/* Function:  sort_err_FuncInfo(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_err_FuncInfo
 *
 *
 * Arg:         obj [UNKN ] Object containing list [FuncInfo *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_err_FuncInfo(FuncInfo * obj,int (*comp)(ErrorInfo *, ErrorInfo *)) 
{
    qsort_err_FuncInfo(obj->err,0,obj->err_len-1,comp);  
    return;  
}    


/* Function:  expand_err_FuncInfo(obj,len)
 *
 * Descrip:    Really an internal function for add_err_FuncInfo
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [FuncInfo *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_err_FuncInfo(FuncInfo * obj,int len) 
{


    if( obj->err_maxlen > obj->err_len )     {  
      warn("expand_FuncInfoerr_ called with no need");   
      return TRUE;   
      }  


    if( (obj->err = (ErrorInfo ** ) ckrealloc (obj->err,sizeof(ErrorInfo *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_FuncInfo, returning FALSE"); 
      return FALSE;  
      }  
    obj->err_maxlen = len;   
    return TRUE; 
}    


/* Function:  add_err_FuncInfo(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [FuncInfo *]
 * Arg:        add [OWNER] Object to add to the list [ErrorInfo *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_err_FuncInfo(FuncInfo * obj,ErrorInfo * add) 
{
    if( obj->err_len >= obj->err_maxlen) {  
      if( expand_err_FuncInfo(obj,obj->err_len + FuncInfoLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->err[obj->err_len++]=add;    
    return TRUE; 
}    


/* Function:  flush_err_FuncInfo(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [FuncInfo *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_err_FuncInfo(FuncInfo * obj) 
{
    int i;   


    for(i=0;i<obj->err_len;i++)  { /*for i over list length*/ 
      if( obj->err[i] != NULL)   {  
        free_ErrorInfo(obj->err[i]); 
        obj->err[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->err_len = 0;    
    return i;    
}    


/* Function:  FuncInfo_alloc_std(void)
 *
 * Descrip:    Equivalent to FuncInfo_alloc_len(FuncInfoLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [FuncInfo *]
 *
 */
FuncInfo * FuncInfo_alloc_std(void) 
{
    return FuncInfo_alloc_len(FuncInfoLISTLENGTH);   
}    


/* Function:  FuncInfo_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [FuncInfo *]
 *
 */
FuncInfo * FuncInfo_alloc_len(int len) 
{
    FuncInfo * out; /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = FuncInfo_alloc()) == NULL) 
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->arg = (ArgInfo   ** ) ckcalloc (len,sizeof(ArgInfo   *))) == NULL)  {  
      warn("Warning, ckcalloc failed in FuncInfo_alloc_len");    
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    if((out->err = (ErrorInfo ** ) ckcalloc (len,sizeof(ErrorInfo *))) == NULL)  {  
      warn("Warning, ckcalloc failed in FuncInfo_alloc_len");    
      return NULL;   
      }  
    out->err_len = 0;    
    out->err_maxlen = len;   


    return out;  
}    


/* Function:  hard_link_FuncInfo(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [FuncInfo *]
 *
 * Return [UNKN ]  Undocumented return value [FuncInfo *]
 *
 */
FuncInfo * hard_link_FuncInfo(FuncInfo * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a FuncInfo object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  FuncInfo_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [FuncInfo *]
 *
 */
FuncInfo * FuncInfo_alloc(void) 
{
    FuncInfo * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(FuncInfo *) ckalloc (sizeof(FuncInfo))) == NULL)    {  
      warn("FuncInfo_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->name = NULL;    
    out->type = NULL;    
    out->complete_name = NULL;   
    out->stripped_return = NULL; 
    out->ft = NULL;  
    out->error = NULL;   
    out->arg = NULL; 
    out->len = out->maxlen = 0;  
    out->err = NULL; 
    out->err_len = out->err_maxlen = 0;  
    out->sdesc = NULL;   
    out->ret = NULL; 
    out->functype = FI_UNKNOWN;  
    out->line_in_c = 0;  
    out->infopos = 0;    
    out->simple = NULL;  
    out->is_hand_written = FALSE;    


    return out;  
}    


/* Function:  free_FuncInfo(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [FuncInfo *]
 *
 * Return [UNKN ]  Undocumented return value [FuncInfo *]
 *
 */
FuncInfo * free_FuncInfo(FuncInfo * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a FuncInfo obj. Should be trappable");  
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
    if( obj->name != NULL)   
      ckfree(obj->name);     
    if( obj->type != NULL)   
      ckfree(obj->type);     
    if( obj->complete_name != NULL)  
      ckfree(obj->complete_name);    
    if( obj->stripped_return != NULL)    
      ckfree(obj->stripped_return);  
    if( obj->ft != NULL) 
      free_Ftext(obj->ft);   
    if( obj->error != NULL)  
      ckfree(obj->error);    
    if( obj->arg != NULL)    {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->arg[i] != NULL) 
          free_ArgInfo(obj->arg[i]); 
        }  
      ckfree(obj->arg);  
      }  
    if( obj->err != NULL)    {  
      for(i=0;i<obj->err_len;i++)    {  
        if( obj->err[i] != NULL) 
          free_ErrorInfo(obj->err[i]);   
        }  
      ckfree(obj->err);  
      }  
    if( obj->sdesc != NULL)  
      ckfree(obj->sdesc);    
    if( obj->ret != NULL)    
      free_ArgInfo(obj->ret);    
    if( obj->simple != NULL) 
      ckfree(obj->simple);   


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_ModuleInfo(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ModuleInfo *]
 *
 * Return [UNKN ]  Undocumented return value [ModuleInfo *]
 *
 */
ModuleInfo * hard_link_ModuleInfo(ModuleInfo * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a ModuleInfo object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  ModuleInfo_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ModuleInfo *]
 *
 */
ModuleInfo * ModuleInfo_alloc(void) 
{
    ModuleInfo * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(ModuleInfo *) ckalloc (sizeof(ModuleInfo))) == NULL)    {  
      warn("ModuleInfo_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->ft = NULL;  


    return out;  
}    


/* Function:  free_ModuleInfo(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ModuleInfo *]
 *
 * Return [UNKN ]  Undocumented return value [ModuleInfo *]
 *
 */
ModuleInfo * free_ModuleInfo(ModuleInfo * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a ModuleInfo obj. Should be trappable");    
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
    if( obj->ft != NULL) 
      free_Ftext(obj->ft);   


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
