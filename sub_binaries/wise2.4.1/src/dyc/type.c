#ifdef _cplusplus
extern "C" {
#endif
#include "type.h"

 char * calc_lex_string;
 int stringpos =0;
 ExprTree * root = NULL;




/* Function:  copy_MethodTypeSet(mts)
 *
 * Descrip:    copies MethodTypeSet by hard-linking list
 *             members
 *
 *
 *
 * Arg:        mts [UNKN ] Undocumented argument [MethodTypeSet *]
 *
 * Return [UNKN ]  Undocumented return value [MethodTypeSet *]
 *
 */
# line 67 "type.dy"
MethodTypeSet * copy_MethodTypeSet(MethodTypeSet * mts)
{
  int i;
  MethodTypeSet * out;

  out = MethodTypeSet_alloc_std();

  for(i=0;i<mts->me_len;i++) 
    add_me_MethodTypeSet(out,hard_link_Method(mts->me[i]));

  for(i=0;i<mts->ty_len;i++) 
    add_ty_MethodTypeSet(out,hard_link_Type(mts->ty[i]));
  
  return out;
}
			 

/* Function:  std_Dynamite_Scope(void)
 *
 * Descrip:    sets up "standard" dynamite scope.
 *
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Scope *]
 *
 */
# line 88 "type.dy"
Scope * std_Dynamite_Scope(void)
{
  Scope * sc;

  sc = Scope_alloc_std();

  /*** std set up code has moved to dyna2.dy ***/
  return sc;
}


/* Function:  allocd_calc_line(calc_line,sc,mts,dycw,pe,expr)
 *
 * Descrip:    Main function used to access the parser.
 *
 *             This parses the calc_line (without damaging it). Errors will
 *             be placed normally: You should stack errors or catch them if
 *             you want to do something funky. 
 *
 *             At the end of the day, if there was parser syntax error then you
 *             will get NULL, and pe & PERR_SYNTAX will be set. Otherwise you
 *             will get a char * (stringalloc'd). Potentially any number of
 *             PERRs could be set, including unscoped variables or methods etc.
 *
 *
 * Arg:        calc_line [READ ] line to be parsed [char *]
 * Arg:               sc [READ ] scope system to use for this area [Scope *]
 * Arg:              mts [READ ] method and type information (method scope) for this area [MethodTypeSet *]
 * Arg:             dycw [UNKN ] Undocumented argument [DycWarning *]
 * Arg:               pe [WRITE] returned parser errors [ParseError *]
 * Arg:             expr [UNKN ] Undocumented argument [ExprTree **]
 *
 * Return [UNKN ]  stringalloc'd expanded string [char *]
 *
 */
# line 117 "type.dy"
char * allocd_calc_line(char * calc_line,Scope * sc,MethodTypeSet * mts,DycWarning * dycw,ParseError * pe,ExprTree ** expr)
{
  char buffer[MAXLINE]; /** ouch - watch overflows **/

  /*** set globals that we share with lex/yacc ***/

  if( calc_line == NULL ) {
    warn("Got a NULL calc line - can hardly parse it!");
    return NULL;
  }

  calc_lex_string = calc_line;
  stringpos= 0;
  root = NULL;

  /*** parse it ***/

  yyparse();

  /*** check root, if NULL... otta here ***/

  if( root == NULL ) {
    /* silent fail. already been warned */
    (*pe) |= PERR_SYNTAX;
    return NULL;
  }


  /*** make top-level name which could be ***/

  find_toplevel_name(root);

  /*** make string, checking semantics this time ***/

  buffer[0]='\0';

  (*pe) |= strcat_ExprTree_Scoped(root,buffer,sc,mts,dycw,NULL,NULL);

  if( expr != NULL ) {
    *expr = root;
  } else {
    /* should free... don't at the moment */
  }

  return stringalloc(buffer);
}

# line 164 "type.dy"
void complain_ParseError_to_file(ParseError pe,FILE * ofp)
{
  if( pe & PERR_SYNTAX ) 
    fprintf(ofp,"Parser Syntax error on calc line\n");
  if( pe & PERR_COMPLEXMETHOD)
    fprintf(ofp,"Complex methods (ie, pointer-functions) which cannot be typed\n");
  if( pe & PERR_ARG_NUM_MIS ) 
    fprintf(ofp,"Argument number mismatches\n");
  if( pe & PERR_ARG_UNTYPED )
    fprintf(ofp,"Un-typable arguments, due to array or structure deferences\n");
  if( pe & PERR_ARG_MISTYPE ) 
    fprintf(ofp,"Mistyped arguments\n");
  if( pe & PERR_OUT_OF_SCOPE) 
    fprintf(ofp,"Variable out of scope\n");
  if( pe & PERR_METHOD_SCOPE)
    fprintf(ofp,"Method out of scope\n");
}


# line 183 "type.dy"
char * type_from_ExprTree(ExprTree * et,Scope * sc,MethodTypeSet * mts)
{
  ScopeUnit * su;
  Method * me;

  switch( et->type ) {
  case ETR_EXPRESSION :
    /*** for the moment, return NULL... could return "int" quite ok.. ***/
    return NULL;
  case ETR_NAME :
    su = ScopeUnit_from_Scope(sc,et->word);
    if( su == NULL ) {
      return NULL;
    } else {
      return su->type;
    }
    break;
  case ETR_TAG :
    if( et->nochild == 1 && et->child[0]->type == ETR_NAME) {
      return type_from_ExprTree(et->child[0],sc,mts);
    }
    return NULL;
    break;
  case ETR_METHOD :
    if( et->child[0]->nochild != 1 || et->child[0]->child[0]->type != ETR_NAME ) {
      return NULL;
    }

    me = Method_from_name(mts,et->child[0]->child[0]->word);
    if( me == NULL) {
      return NULL;
    } else {
      return me->retstr;
    }
    break;
  case ETR_STRUCTREF :
  case ETR_REFERENCE :
    return NULL;
  default :
    warn("Unable to type a Expr Node [%d]. Returing a 'blank' type",et->type);
    return NULL;
  }

}

/* Function:  strcat_ExprTree_Scoped(ExprTree,finish_parsing,buffer,sc,mts,dycw,data)
 *
 * Descrip:    Main internal recursive functions that
 *             descends the Expr tree. Should at the
 *             start be given the root node: buffer
 *             is written in the final expression, mapped
 *             etc. 
 *
 *
 * Arg:              ExprTree [UNKN ] Undocumented argument [ExprTree *]
 * Arg:        finish_parsing [UNKN ] Undocumented argument [NullString]
 * Arg:                buffer [UNKN ] Undocumented argument [char *]
 * Arg:                    sc [UNKN ] Undocumented argument [Scope *]
 * Arg:                   mts [UNKN ] Undocumented argument [MethodTypeSet *]
 * Arg:                  dycw [UNKN ] Undocumented argument [DycWarning *]
 * Arg:                  data [UNKN ] Undocumented argument [void *]
 *
 * Return [UNKN ]  Undocumented return value [ParseError]
 *
 */
# line 236 "type.dy"
ParseError strcat_ExprTree_Scoped(ExprTree * ExprTree,char * buffer,Scope * sc,MethodTypeSet * mts,DycWarning * dycw,boolean (*finish_parsing)(struct ExprTree *,char *,void *),void * data)
{
  ParseError pe = 0; 
  ScopeUnit * su;
  Method * me;
  char tempbuf[512];
  int i;


  if( finish_parsing != NULL ) {
    if( (*finish_parsing)(ExprTree,buffer,data) == TRUE ) {
      return pe;
    }
  }

  switch(ExprTree->type) {
  case ETR_NUMBER : strcat(buffer,ExprTree->word); 
        break;

  case ETR_OPERATOR : strcat(buffer,ExprTree->word); 
    break;

  case ETR_EXPRESSION : strcat(buffer,"(");
    for(i=0;i< ExprTree->nochild;i++)
      pe |= strcat_ExprTree_Scoped(ExprTree->child[i],buffer,sc,mts,dycw,finish_parsing,data);
    strcat(buffer,")");
    break;
  case ETR_STATEMENT :
    for(i=0;i< ExprTree->nochild;i++)
      pe |= strcat_ExprTree_Scoped(ExprTree->child[i],buffer,sc,mts,dycw,finish_parsing,data);
    break;
  case ETR_NAME : 
    if( (ExprTree->attrib & IS_TOPLEVEL) == IS_TOPLEVEL) {
      /* fprintf(stderr,"Got a top level name %s\n",ExprTree->word);*/
      su = ScopeUnit_from_Scope(sc,ExprTree->word);
      if( su == NULL ) {
	pe |= PERR_OUT_OF_SCOPE;
	if( dycw == NULL || dycw->warn_extern == TRUE )
	  warn("Name [%s] is out of scope: assumming extern",ExprTree->word);
      } else {
	/*	fprintf(stderr,"Ok.[%d] for name %s, going to add %s\n",ExprTree->attrib,ExprTree->word,su->app);*/
	strcat(buffer,su->app);
      }
    }
    strcat(buffer,ExprTree->word);
    break;
  case ETR_ARRAY : 
    pe |= strcat_ExprTree_Scoped(ExprTree->child[0],buffer,sc,mts,dycw,finish_parsing,data);
    strcat(buffer,"[");
    pe |= strcat_ExprTree_Scoped(ExprTree->child[1],buffer,sc,mts,dycw,finish_parsing,data);
    strcat(buffer,"]");
    break;
  case ETR_TAG : 
    for(i=0;i< ExprTree->nochild;i++)
      pe |= strcat_ExprTree_Scoped(ExprTree->child[i],buffer,sc,mts,dycw,finish_parsing,data);
    break;
  case ETR_STRUCTREF :
    pe |= strcat_ExprTree_Scoped(ExprTree->child[0],buffer,sc,mts,dycw,finish_parsing,data);
    strcat(buffer,ExprTree->child[1]->word);
    pe |= strcat_ExprTree_Scoped(ExprTree->child[2],buffer,sc,mts,dycw,finish_parsing,data);
    break;
  case ETR_REFERENCE :
    strcat(buffer,ExprTree->child[0]->word);
    pe |= strcat_ExprTree_Scoped(ExprTree->child[1],buffer,sc,mts,dycw,finish_parsing,data);
    break;
  case ETR_METHOD :
    
    if( ExprTree->child[0]->nochild != 1 || ExprTree->child[0]->child[0]->type != ETR_NAME ) {
      pe |= PERR_COMPLEXMETHOD;
      
      /*** generate error ***/
      tempbuf[0]='\0';
      strcat_ExprTree(ExprTree->child[0],tempbuf);
      warn("Method [%s] is not a pure name [%d]: this tolerated in the parser but can't be checked",tempbuf,ExprTree->child[0]->type);
      pe |= strcat_ExprTree_Scoped(ExprTree->child[0],buffer,sc,mts,dycw,finish_parsing,data);

      
    } else {
      me = Method_from_name(mts,ExprTree->child[0]->child[0]->word);
      if( me == NULL ) {
	pe |= PERR_METHOD_SCOPE;
	if( dycw == NULL || dycw->warn_extern_method == TRUE )
	  warn("Implied method [%s] has not been typed, and hence can't be type-checked or logical-mapped",ExprTree->child[0]->child[0]->word);
	pe |= strcat_ExprTree_Scoped(ExprTree->child[0],buffer,sc,mts,dycw,finish_parsing,data);
      } else { 
	pe |= typecheck_method(ExprTree->child[1],me,sc,mts);
	strcat(buffer,me->real);
	/* do nothing at the moment */
      }
    }

    strcat(buffer,"(");
    pe |= strcat_ExprTree_Scoped(ExprTree->child[1],buffer,sc,mts,dycw,finish_parsing,data);
    strcat(buffer,")");

    break;
  case ETR_COMMALIST :
    for(i=0;i<ExprTree->nochild;i++) {
      pe |= strcat_ExprTree_Scoped(ExprTree->child[i],buffer,sc,mts,dycw,finish_parsing,data);
      if( i != ExprTree->nochild-1)
	strcat(buffer,",");
    }
    break;

  default :
    pe |= PERR_SYNTAX;
      warn("In trying to make Expr string, got unobtainable type!");
  }

  return pe;
}

# line 348 "type.dy"
ParseError typecheck_method(ExprTree * et,Method * me,Scope * sc,MethodTypeSet * mts)
{
  int i;
  char * t;
  int ret = 0;
  int j;
  ScopeUnit * sunit = NULL;

  if( et->type != ETR_COMMALIST ) {
    warn("Trying to typecheck with a non-comma list system %d --- internal parser error and v.bad",et->type);
  }

  if( me->len != et->nochild ) {
    ret |= PERR_ARG_NUM_MIS;
    warn("In method [%s], expect %d arguments, given %d arguments",me->logical,me->len,et->nochild);
    return ret;
  }

  for(i=0;i<me->len;i++) {
    t = type_from_ExprTree(et->child[i],sc,mts);
    if( t == NULL ) {
      ret |= PERR_ARG_UNTYPED; 
      /* warn("For argument %d, of %s, no type information",i,me->logical); */
    } else {
      if ( compare_type(t,me->ma[i]->type) == FALSE ) {
	warn("Mis-type in argument %d of %s: wanted [%s] got [%s]",i+1,me->logical,me->ma[i]->type,t);
	ret |= PERR_ARG_MISTYPE;
      }
    }
  }

  /** do forbidden pairs if appropiate **/

  for(i=0;i<me->len;i++) {
    auto char * test_argument;
    if( et->child[i]->type == ETR_NAME ) {
      test_argument = et->child[i]->word;
      if( (sunit = ScopeUnit_from_Scope(sc,et->child[i]->word)) == NULL ) 
	continue; /** ugh - unscoped **/
    } else if ( et->child[i]->type == ETR_TAG && et->child[i]->nochild == 1 && et->child[i]->child[0]->type == ETR_NAME ) {
      test_argument = et->child[i]->child[0]->word;
      if( (sunit = ScopeUnit_from_Scope(sc,et->child[i]->child[0]->word)) == NULL ) 
	continue; /** ugh - unscoped **/
    }
    if( sunit == NULL ) {
      continue;
    }

    /** we have scope. Continue if it has no forbiddens **/
    if( sunit->no_accept == NULL )
      continue;

    for(j=i+1;j<me->len;j++) {
      auto char * word;
      word = NULL;
      if( et->child[j]->type == ETR_NAME ) {
	word = et->child[j]->word;

      } else if ( et->child[i]->type == ETR_TAG && et->child[i]->nochild == 1 && et->child[i]->child[0]->type == ETR_NAME ) {
	word = et->child[j]->child[0]->word;
      }
      
      if( word != NULL && strcmp(word,sunit->no_accept) == 0 && 0) 
	warn("For function %s, you have arguments %s and %s, which do not expect to paired directly in a function. This is just a warning that you can ignore",me->logical,word,test_argument);
    }
  }

  /*** end of forbidden pairs code ***/
      


  return ret;
}


/* Function:  ScopeUnit_from_Scope(sc,word)
 *
 * Descrip:    gets a ScopeUnit from the name
 *
 *
 * Arg:          sc [UNKN ] Undocumented argument [Scope *]
 * Arg:        word [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [ScopeUnit *]
 *
 */
# line 427 "type.dy"
ScopeUnit * ScopeUnit_from_Scope(Scope * sc,char * word)
{
  int i;

  for(i=0;i<sc->len;i++) {
    if( sc->su[i]->isglobbed == TRUE) {
      /*** hacky ***/
      if( strstartcmp(word,sc->su[i]->name) == 0 ) {
	return sc->su[i];
      }
    } else {
      if( strcmp(sc->su[i]->name,word) == 0 )
	return sc->su[i];
    }
  }

  return NULL;
}

# line 446 "type.dy"
ScopeUnit * ScopeUnit_from_nat(MethodTypeSet * mts,char * name,char * app,char * type)
{
  ScopeUnit * su;

  su = ScopeUnit_alloc();

  if( name[strlen(name)-1] == '*') {
    name[strlen(name)-1] = '\0';
    su->isglobbed = TRUE;
  }

  su->name = stringalloc(name);
  su->app = stringalloc(app);
  su->type = stringalloc(type);

  return su;
}


# line 452 "type.c"
/* Function:  hard_link_ScopeUnit(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ScopeUnit *]
 *
 * Return [UNKN ]  Undocumented return value [ScopeUnit *]
 *
 */
ScopeUnit * hard_link_ScopeUnit(ScopeUnit * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a ScopeUnit object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  ScopeUnit_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ScopeUnit *]
 *
 */
ScopeUnit * ScopeUnit_alloc(void) 
{
    ScopeUnit * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(ScopeUnit *) ckalloc (sizeof(ScopeUnit))) == NULL)  {  
      warn("ScopeUnit_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->name = NULL;    
    out->app = NULL; 
    out->type = NULL;    
    out->isglobbed = FALSE;  
    out->scope_type = SCOPE_EXTERN;  
    out->no_accept = NULL;   


    return out;  
}    


/* Function:  free_ScopeUnit(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ScopeUnit *]
 *
 * Return [UNKN ]  Undocumented return value [ScopeUnit *]
 *
 */
ScopeUnit * free_ScopeUnit(ScopeUnit * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a ScopeUnit obj. Should be trappable"); 
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
    if( obj->app != NULL)    
      ckfree(obj->app);  
    if( obj->type != NULL)   
      ckfree(obj->type);     
    if( obj->no_accept != NULL)  
      ckfree(obj->no_accept);    


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_Scope(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_Scope
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [ScopeUnit **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_Scope(ScopeUnit ** list,int i,int j)  
{
    ScopeUnit * temp;    
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_Scope(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_Scope which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [ScopeUnit **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_Scope(ScopeUnit ** list,int left,int right,int (*comp)(ScopeUnit * ,ScopeUnit * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_Scope(list,left,(left+right)/2);    
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_Scope (list,++last,i);  
      }  
    swap_Scope (list,left,last); 
    qsort_Scope(list,left,last-1,comp);  
    qsort_Scope(list,last+1,right,comp); 
}    


/* Function:  sort_Scope(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_Scope
 *
 *
 * Arg:         obj [UNKN ] Object containing list [Scope *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_Scope(Scope * obj,int (*comp)(ScopeUnit *, ScopeUnit *)) 
{
    qsort_Scope(obj->su,0,obj->len-1,comp);  
    return;  
}    


/* Function:  expand_Scope(obj,len)
 *
 * Descrip:    Really an internal function for add_Scope
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Scope *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_Scope(Scope * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_Scope called with no need");  
      return TRUE;   
      }  


    if( (obj->su = (ScopeUnit ** ) ckrealloc (obj->su,sizeof(ScopeUnit *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_Scope, returning FALSE");    
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_Scope(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Scope *]
 * Arg:        add [OWNER] Object to add to the list [ScopeUnit *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_Scope(Scope * obj,ScopeUnit * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_Scope(obj,obj->len + ScopeLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->su[obj->len++]=add; 
    return TRUE; 
}    


/* Function:  flush_Scope(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [Scope *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_Scope(Scope * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->su[i] != NULL)    {  
        free_ScopeUnit(obj->su[i]);  
        obj->su[i] = NULL;   
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  Scope_alloc_std(void)
 *
 * Descrip:    Equivalent to Scope_alloc_len(ScopeLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Scope *]
 *
 */
Scope * Scope_alloc_std(void) 
{
    return Scope_alloc_len(ScopeLISTLENGTH); 
}    


/* Function:  Scope_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [Scope *]
 *
 */
Scope * Scope_alloc_len(int len) 
{
    Scope * out;/* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = Scope_alloc()) == NULL)    
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->su = (ScopeUnit ** ) ckcalloc (len,sizeof(ScopeUnit *))) == NULL)   {  
      warn("Warning, ckcalloc failed in Scope_alloc_len");   
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_Scope(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Scope *]
 *
 * Return [UNKN ]  Undocumented return value [Scope *]
 *
 */
Scope * hard_link_Scope(Scope * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a Scope object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  Scope_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Scope *]
 *
 */
Scope * Scope_alloc(void) 
{
    Scope * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(Scope *) ckalloc (sizeof(Scope))) == NULL)  {  
      warn("Scope_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->su = NULL;  
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_Scope(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Scope *]
 *
 * Return [UNKN ]  Undocumented return value [Scope *]
 *
 */
Scope * free_Scope(Scope * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a Scope obj. Should be trappable"); 
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
    if( obj->su != NULL) {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->su[i] != NULL)  
          free_ScopeUnit(obj->su[i]);    
        }  
      ckfree(obj->su);   
      }  


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
