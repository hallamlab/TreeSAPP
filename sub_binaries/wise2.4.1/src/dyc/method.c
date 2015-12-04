#ifdef _cplusplus
extern "C" {
#endif
#include "method.h"




/* Function:  Type_from_name(mts,name)
 *
 * Descrip:    gets a type by its name
 *
 *
 * Arg:         mts [UNKN ] Undocumented argument [MethodTypeSet *]
 * Arg:        name [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [Type *]
 *
 */
# line 61 "method.dy"
Type * Type_from_name(MethodTypeSet * mts,char * name)
{
  int i;
  Type * out = NULL;

  if( mts == NULL ) {
    warn("Attempting to get a type name from a null mts. Nope!");
    return NULL;
  }

  for(i=0;i<mts->ty_len;i++) {
    if( strcmp(mts->ty[i]->logical,name) == 0 ) {
      if( out == NULL ) {
	out = mts->ty[i];
      } else {
	warn("Multiple definitions for %s - taking the last one\n",name);
	out = mts->ty[i];
      }
    }
  }

  return out;
}



/* Function:  compare_type(s,t)
 *
 * Descrip:    Essentially compares two strings, disregarding white space
 *
 *             Not ideal!!!
 *
 *
 * Arg:        s [UNKN ] Undocumented argument [char *]
 * Arg:        t [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 93 "method.dy"
boolean compare_type(char * s,char * t)
{
  for(;*s && isspace(*s);s++)
    ;
  for(;*t && isspace(*t);t++) 
    ;

  for(;*s && *t && *s == *t;) {
    for(s++;*s && isspace(*s);s++)
      ;
    for(t++;*t && isspace(*t);t++)
      ;
  }

  if( *s == '\0' && *t == '\0' )
    return TRUE;

  return FALSE;
}

/* Function:  Method_from_name(mts,name)
 *
 * Descrip:    gets a Method by its name
 *
 *
 * Arg:         mts [UNKN ] Undocumented argument [MethodTypeSet *]
 * Arg:        name [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [Method *]
 *
 */
# line 117 "method.dy"
Method * Method_from_name(MethodTypeSet * mts,char * name)
{
  int i;

  if( mts == NULL) {
    warn("Attempting to get a method name from a null mts. Nope!");
    return NULL;
  }

  for(i=0;i<mts->me_len;i++) {
    if( strcmp(mts->me[i]->logical,name) == 0 )
      return mts->me[i];
  }

  return NULL;
}


/*** I/O on methods/types ***/

# line 137 "method.dy"
void show_MethodTypeSet(MethodTypeSet * mts,FILE * ofp)
{
  int i;

  for(i=0;i<mts->me_len;i++)
    show_Method(mts->me[i],ofp);

  for(i=0;i<mts->ty_len;i++)
    show_Type(mts->ty[i],ofp);

}
  
# line 149 "method.dy"
void show_Method(Method * m,FILE * ofp)
{
  int i;

  fprintf(ofp,"Method [%s] map [%s]\n",m->logical,m->real);
  for(i=0;i<m->len;i++)
    show_MethodArg(m->ma[i],ofp);
}

# line 158 "method.dy"
void show_MethodArg(MethodArg * ma,FILE * ofp)
{
  fprintf(ofp,"Argument: %s\n",ma->type);
}

# line 163 "method.dy"
void show_Type(Type * ty,FILE * ofp)
{
  fprintf(ofp,"Type: Logial %s Real %s\n",ty->logical,ty->real);
}

/* Function:  StructElement_from_MethodTypeSet(name,type,mts)
 *
 * Descrip:    function which handles the logical->real mapping
 *
 *             At the moment "unmappable" types get assummed to be C types,
 *             trigger a warning and return the correct thing.
 *
 *
 * Arg:        name [UNKN ] Undocumented argument [char *]
 * Arg:        type [UNKN ] Undocumented argument [char *]
 * Arg:         mts [UNKN ] Undocumented argument [MethodTypeSet *]
 *
 * Return [UNKN ]  Undocumented return value [StructElement *]
 *
 */
# line 174 "method.dy"
StructElement * StructElement_from_MethodTypeSet(char * name,char * type,MethodTypeSet * mts)
{
  Type * ty;

  if( (ty = Type_from_name(mts,type)) == NULL ) {
    warn("Type [%s] is not recognised as a logical Dynamite type. Assumming it is a real C type",type);
    return StructElement_from_nameandtype(name,type);
  }

  return StructElement_from_nameandtype(name,ty->real);

}

/* Function:  StructElement_from_nameandtype(name,type)
 *
 * Descrip:    an internal for StructElement_from_MethodTypeSet. don't use otherwise please!
 *
 *
 * Arg:        name [UNKN ] Undocumented argument [char *]
 * Arg:        type [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [StructElement *]
 *
 */
# line 191 "method.dy"
StructElement * StructElement_from_nameandtype(char * name,char * type)
{
  StructElement * out;

  out = StructElement_alloc();

  out->name = stringalloc(name);
  out->element_type = stringalloc(type);

  out->islinked = TRUE;

  return out;
}
  


# line 207 "method.dy"
MethodTypeSet * read_MethodTypeSet_filename(char * filename)
{
  FILE * ifp;
  MethodTypeSet * mts;

  ifp=openfile(filename,"r");

  if( ifp == NULL ) {
    warn("Could not open [%s] for MethodTypeSet reading",filename);
    return FALSE;
  }

  mts = read_MethodTypeSet(ifp);

  fclose(ifp);

  return mts;
}

# line 226 "method.dy"
boolean read_into_MethodTypeSet_filename(MethodTypeSet * mts,char * filename)
{
  FILE * ifp;
  boolean ret;

  ifp=openfile(filename,"r");

  if( ifp == NULL ) {
    warn("Could not open [%s] for MethodTypeSet reading",filename);
    return FALSE;
  }

  ret = read_into_MethodTypeSet(mts,ifp);

  fclose(ifp);

  return ret;
}

# line 245 "method.dy"
MethodTypeSet * standard_dynamite_MethodTypeSet(void)
{
  MethodTypeSet * mts;
  Type * temp;

  mts = empty_MethodTypeSet();

  temp = Type_alloc();
  temp->logical=stringalloc("int");
  temp->real=stringalloc("int");
  add_ty_MethodTypeSet(mts,temp);

  temp = Type_alloc();
  temp->logical=stringalloc("double");
  temp->real=stringalloc("double");
  add_ty_MethodTypeSet(mts,temp);

  temp = Type_alloc();
  temp->logical=stringalloc("Score");
  temp->real=stringalloc("Score");
  add_ty_MethodTypeSet(mts,temp);


  return mts;
}


# line 272 "method.dy"
MethodTypeSet * empty_MethodTypeSet(void)
{
  return MethodTypeSet_alloc_std();
}

# line 277 "method.dy"
MethodTypeSet * read_MethodTypeSet(FILE * ifp)
{
  MethodTypeSet * mts;

  mts = empty_MethodTypeSet();

  read_into_MethodTypeSet(mts,ifp);
   
  return mts;
}

# line 288 "method.dy"
boolean read_into_MethodTypeSet(MethodTypeSet * mts,FILE * ifp)
{
  char buffer[MAXLINE];
  Method * me;
  Type * ty;
  Input * in;

  while( fgets(buffer,MAXLINE,ifp) != NULL) {
    chop_newline(buffer);
    
    if( buffer[0] == '#' || strwhitestartcmp(buffer,"#",spacestr) == 0 )
      continue;

    if( only_whitespace(buffer,spacestr) == TRUE) 
      continue;

    if( strstartcmp(buffer,"method") == 0 ) {
      if( (me=read_Method_line(buffer,ifp)) == NULL ) {
	warn("Unable to read method in line [%s] ",buffer);
      } else {
	add_me_MethodTypeSet(mts,me);
      }
    } else if ( strstartcmp(buffer,"type") == 0 ) {
      if( (ty=read_Type_line(buffer,ifp)) == NULL ) {
	warn("Unable to read type in line [%s] ",buffer);
      } else {
	add_ty_MethodTypeSet(mts,ty);
      }
    } else if ( strstartcmp(buffer,"input") == 0 ) {
      if( (in = read_Input_line(buffer,ifp)) == NULL ) {
	warn("Unable to read type in line [%s]",buffer);
      } else {
	add_in_MethodTypeSet(mts,in);
      }
    } else {
      warn("In reading only method/types got an impossible line [%s]",buffer);
    }
  }

  return TRUE;
}


# line 331 "method.dy"
boolean is_database_type(Type * ty)
{
  if( ty->database_type == NULL ||  ty->reload_func == NULL || ty->init_func == NULL || ty->close_func == NULL)
    return FALSE;

  return TRUE;
}


/* Function:  read_Type_line(line,ifp)
 *
 * Descrip:    reads in a type structure from a line starting
 *
 *             type
 *             etc
 *
 *
 * Arg:        line [UNKN ] first line with type [char *]
 * Arg:         ifp [UNKN ] read file [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [Type *]
 *
 */
# line 349 "method.dy"
Type * read_Type_line(char * line,FILE * ifp)
{
  Type * out;
  char * temp;
  char buffer[MAXLINE];

  if( strstartcmp(line,"type") != 0 ) {
    warn("Attempting to read a method with no method line!");
    return NULL;
  }


  out = Type_alloc();
  out->logical = second_word_alloc(line,spacestr);

  while( fgets(buffer,MAXLINE,ifp) != NULL ) {
    chop_newline(buffer);
    if( strstartcmp(buffer,"end") == 0 ) 
      break;
    else if( strstartcmp(buffer,"real") == 0 ) {
      temp = second_word_alloc(buffer,spacestr);
      if( out->real != NULL ) { 
	warn("For type [%s], second real replacing [%s] with [%s]",out->logical,out->real,temp);
	ckfree(out->real);
      } 
      out->real = temp;
    } else if( strstartcmp(buffer,"threadsafe") == 0 ) {
      out->is_thread_safe = TRUE;
    } else if ( strstartcmp(buffer,"dbtype") == 0 ) {
      temp = second_word_alloc(buffer,spacestr);
      if( out->database_type != NULL ) { 
	warn("For type [%s], second database type replacing [%s] with [%s]",out->logical,out->database_type,temp);
	ckfree(out->database_type);
      } 
      out->database_type = temp;
    } else if ( strstartcmp(buffer,"init") == 0 ) {
      temp = second_word_alloc(buffer,spacestr);
      if( out->init_func != NULL ) { 
	warn("For type [%s], second init replacing [%s] with [%s]",out->logical,out->init_func,temp);
	ckfree(out->init_func);
      } 
      out->init_func = temp;
    } else if ( strstartcmp(buffer,"maxlen") == 0 ) {
      temp = second_word_alloc(buffer,spacestr);
      if( out->maxlen != NULL ) { 
	warn("For type [%s], second maxlen replacing [%s] with [%s]",out->logical,out->maxlen,temp);
	ckfree(out->maxlen);
      } 
      out->maxlen = temp;
    } else if ( strstartcmp(buffer,"reload") == 0 ) {
      temp = second_word_alloc(buffer,spacestr);
      if( out->reload_func != NULL ) { 
	warn("For type [%s], second reload function replacing [%s] with [%s]",out->logical,out->reload_func,temp);
	ckfree(out->reload_func);
      } 
      out->reload_func = temp;
    } else if ( strstartcmp(buffer,"addentry") == 0 ) {
      temp = second_word_alloc(buffer,spacestr);
      if( out->dataentry_add != NULL ) { 
	warn("For type [%s], second dataentry_add function replacing [%s] with [%s]",out->logical,out->dataentry_add,temp);
	ckfree(out->dataentry_add);
      } 
      out->dataentry_add = temp;
    } else if ( strstartcmp(buffer,"close") == 0 ) {
      temp = second_word_alloc(buffer,spacestr);
      if( out->close_func != NULL ) { 
	warn("For type [%s], second close func replacing [%s] with [%s]",out->logical,out->close_func,temp);
	ckfree(out->close_func);
      } 
      out->close_func = temp;
    } else if ( strstartcmp(buffer,"hardlink") == 0 ) {
      temp = second_word_alloc(buffer,spacestr);
      if( out->hard_link_func != NULL ) { 
	warn("For type [%s], second hard_link func replacing [%s] with [%s]",out->logical,out->hard_link_func,temp);
	ckfree(out->hard_link_func);
      } 
      out->hard_link_func = temp;
    } else if ( strstartcmp(buffer,"free") == 0 ) {
      temp = second_word_alloc(buffer,spacestr);
      if( out->free_func != NULL ) { 
	warn("For type [%s], second free func replacing [%s] with [%s]",out->logical,out->free_func,temp);
	ckfree(out->free_func);
      } 
      out->free_func = temp;
    } else if ( strstartcmp(buffer,"name") == 0 ) {
      temp = second_word_alloc(buffer,spacestr);
      if( out->get_id_func != NULL ) { 
	warn("For type [%s], second get name func replacing [%s] with [%s]",out->logical,out->get_id_func,temp);
	ckfree(out->get_id_func);
      } 
      out->get_id_func = temp;
    } else {
      warn("In reading type [%s] did not understand [%s]",out->logical,buffer);
    }
  }
  if( out->is_thread_safe == TRUE ) {
    if( out->hard_link_func == NULL || out->free_func == NULL ) {
      warn("Trying to make type %s threadsafe but have not supplied hardlink and free functions",out->logical);
      out->is_thread_safe = FALSE;
    }
  }

  out->is_database = is_database_type(out);
  return out;
}
  
# line 455 "method.dy"
Method * read_Method_line(char * line,FILE * ifp)
{
  Method * out;
  char buffer[MAXLINE];
  char * temp;
  MethodArg * m;

  if( strstartcmp(line,"method") != 0 ) {
    warn("Attempting to read a method with no method line!");
    return NULL;
  }
  
  out = Method_alloc_std();
  out->logical = second_word_alloc(line,spacestr);

  while( fgets(buffer,MAXLINE,ifp) != NULL ) {
    chop_newline(buffer);
    if( strstartcmp(buffer,"end") == 0 ) 
      break;
    else if( strstartcmp(buffer,"map") == 0 ) {
      temp = second_word_alloc(buffer,spacestr);
      if( out->real != NULL ) { 
	warn("For method [%s], second map replacing [%s] with [%s]",out->logical,out->real,temp);
	ckfree(out->real);
      } 
      out->real = temp;
    } else if ( strstartcmp(buffer,"arg") == 0 ) {
      m = MethodArg_from_line(buffer);
      if( m == NULL ) {
	warn("Got a NULL method arg. Yikes!");
      } else {
	add_Method(out,m);
      }
    } else if ( strstartcmp(buffer,"return") == 0 ) {
      out->retstr = second_word_alloc(buffer,spacestr);
    } else {
      warn("In method [%s] did not understand %s",out->logical,buffer);
    }
  }

  return out;
}
  

  
/* Function:  MethodArg_from_line(line)
 *
 * Descrip:    reads in from line like
 *
 *             arg PROTEIN
 *
 *
 * Arg:        line [UNKN ] pointer to buffer [char *]
 *
 * Return [UNKN ]  Undocumented return value [MethodArg *]
 *
 */
# line 507 "method.dy"
MethodArg * MethodArg_from_line(char * line)
{
  MethodArg * out;

  if( strstartcmp(line,"arg") != 0 ) {
    warn("Attempting to read a method argument with no arg line!");
    return NULL;
  }

  out = MethodArg_alloc();

  out->type = second_word_alloc(line,spacestr);

  return out;
}


# line 527 "method.c"
/* Function:  hard_link_MethodArg(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [MethodArg *]
 *
 * Return [UNKN ]  Undocumented return value [MethodArg *]
 *
 */
MethodArg * hard_link_MethodArg(MethodArg * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a MethodArg object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  MethodArg_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [MethodArg *]
 *
 */
MethodArg * MethodArg_alloc(void) 
{
    MethodArg * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(MethodArg *) ckalloc (sizeof(MethodArg))) == NULL)  {  
      warn("MethodArg_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->type = NULL;    


    return out;  
}    


/* Function:  free_MethodArg(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [MethodArg *]
 *
 * Return [UNKN ]  Undocumented return value [MethodArg *]
 *
 */
MethodArg * free_MethodArg(MethodArg * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a MethodArg obj. Should be trappable"); 
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
    if( obj->type != NULL)   
      ckfree(obj->type);     


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_Method(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_Method
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [MethodArg **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_Method(MethodArg ** list,int i,int j)  
{
    MethodArg * temp;    
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_Method(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_Method which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [MethodArg **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_Method(MethodArg ** list,int left,int right,int (*comp)(MethodArg * ,MethodArg * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_Method(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_Method (list,++last,i); 
      }  
    swap_Method (list,left,last);    
    qsort_Method(list,left,last-1,comp); 
    qsort_Method(list,last+1,right,comp);    
}    


/* Function:  sort_Method(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_Method
 *
 *
 * Arg:         obj [UNKN ] Object containing list [Method *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_Method(Method * obj,int (*comp)(MethodArg *, MethodArg *)) 
{
    qsort_Method(obj->ma,0,obj->len-1,comp); 
    return;  
}    


/* Function:  expand_Method(obj,len)
 *
 * Descrip:    Really an internal function for add_Method
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Method *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_Method(Method * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_Method called with no need"); 
      return TRUE;   
      }  


    if( (obj->ma = (MethodArg ** ) ckrealloc (obj->ma,sizeof(MethodArg *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_Method, returning FALSE");   
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_Method(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Method *]
 * Arg:        add [OWNER] Object to add to the list [MethodArg *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_Method(Method * obj,MethodArg * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_Method(obj,obj->len + MethodLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->ma[obj->len++]=add; 
    return TRUE; 
}    


/* Function:  flush_Method(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [Method *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_Method(Method * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->ma[i] != NULL)    {  
        free_MethodArg(obj->ma[i]);  
        obj->ma[i] = NULL;   
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  Method_alloc_std(void)
 *
 * Descrip:    Equivalent to Method_alloc_len(MethodLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Method *]
 *
 */
Method * Method_alloc_std(void) 
{
    return Method_alloc_len(MethodLISTLENGTH);   
}    


/* Function:  Method_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [Method *]
 *
 */
Method * Method_alloc_len(int len) 
{
    Method * out;   /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = Method_alloc()) == NULL)   
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->ma = (MethodArg ** ) ckcalloc (len,sizeof(MethodArg *))) == NULL)   {  
      warn("Warning, ckcalloc failed in Method_alloc_len");  
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_Method(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Method *]
 *
 * Return [UNKN ]  Undocumented return value [Method *]
 *
 */
Method * hard_link_Method(Method * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a Method object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  Method_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Method *]
 *
 */
Method * Method_alloc(void) 
{
    Method * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(Method *) ckalloc (sizeof(Method))) == NULL)    {  
      warn("Method_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->logical = NULL; 
    out->real = NULL;    
    out->retstr = NULL;  
    out->ma = NULL;  
    out->len = out->maxlen = 0;  
    out->cugen_map = NULL;   
    out->cugen_type = CUGEN_METHOD_UNKNOWN;  


    return out;  
}    


/* Function:  free_Method(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Method *]
 *
 * Return [UNKN ]  Undocumented return value [Method *]
 *
 */
Method * free_Method(Method * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a Method obj. Should be trappable");    
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
    if( obj->logical != NULL)    
      ckfree(obj->logical);  
    if( obj->real != NULL)   
      ckfree(obj->real);     
    if( obj->retstr != NULL) 
      ckfree(obj->retstr);   
    if( obj->ma != NULL) {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->ma[i] != NULL)  
          free_MethodArg(obj->ma[i]);    
        }  
      ckfree(obj->ma);   
      }  
    if( obj->cugen_map != NULL)  
      ckfree(obj->cugen_map);    


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_Type(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Type *]
 *
 * Return [UNKN ]  Undocumented return value [Type *]
 *
 */
Type * hard_link_Type(Type * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a Type object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  Type_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Type *]
 *
 */
Type * Type_alloc(void) 
{
    Type * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(Type *) ckalloc (sizeof(Type))) == NULL)    {  
      warn("Type_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->is_database = FALSE;    
    out->is_thread_safe = FALSE; 
    out->logical = NULL; 
    out->real = NULL;    
    out->database_type = NULL;   
    out->get_id_func = NULL; 
    out->init_func = NULL;   
    out->reload_func = NULL; 
    out->close_func = NULL;  
    out->dataentry_add = NULL;   
    out->maxlen = NULL;  
    out->hard_link_func = NULL;  
    out->free_func = NULL;   
    out->in = NULL;  


    return out;  
}    


/* Function:  free_Type(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Type *]
 *
 * Return [UNKN ]  Undocumented return value [Type *]
 *
 */
Type * free_Type(Type * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a Type obj. Should be trappable");  
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
    if( obj->logical != NULL)    
      ckfree(obj->logical);  
    if( obj->real != NULL)   
      ckfree(obj->real);     
    if( obj->database_type != NULL)  
      ckfree(obj->database_type);    
    if( obj->get_id_func != NULL)    
      ckfree(obj->get_id_func);  
    if( obj->init_func != NULL)  
      ckfree(obj->init_func);    
    if( obj->reload_func != NULL)    
      ckfree(obj->reload_func);  
    if( obj->close_func != NULL) 
      ckfree(obj->close_func);   
    if( obj->dataentry_add != NULL)  
      ckfree(obj->dataentry_add);    
    if( obj->maxlen != NULL) 
      ckfree(obj->maxlen);   
    if( obj->hard_link_func != NULL) 
      ckfree(obj->hard_link_func);   
    if( obj->free_func != NULL)  
      ckfree(obj->free_func);    
    if( obj->in != NULL) 
      free_Input(obj->in);   


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_me_MethodTypeSet(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_me_MethodTypeSet
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [Method **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_me_MethodTypeSet(Method ** list,int i,int j)  
{
    Method * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_me_MethodTypeSet(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_me_MethodTypeSet which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [Method **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_me_MethodTypeSet(Method ** list,int left,int right,int (*comp)(Method * ,Method * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_me_MethodTypeSet(list,left,(left+right)/2); 
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_me_MethodTypeSet (list,++last,i);   
      }  
    swap_me_MethodTypeSet (list,left,last);  
    qsort_me_MethodTypeSet(list,left,last-1,comp);   
    qsort_me_MethodTypeSet(list,last+1,right,comp);  
}    


/* Function:  sort_me_MethodTypeSet(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_me_MethodTypeSet
 *
 *
 * Arg:         obj [UNKN ] Object containing list [MethodTypeSet *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_me_MethodTypeSet(MethodTypeSet * obj,int (*comp)(Method *, Method *)) 
{
    qsort_me_MethodTypeSet(obj->me,0,obj->me_len-1,comp);    
    return;  
}    


/* Function:  expand_me_MethodTypeSet(obj,len)
 *
 * Descrip:    Really an internal function for add_me_MethodTypeSet
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [MethodTypeSet *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_me_MethodTypeSet(MethodTypeSet * obj,int len) 
{


    if( obj->me_maxlen > obj->me_len )   {  
      warn("expand_MethodTypeSetme_ called with no need");   
      return TRUE;   
      }  


    if( (obj->me = (Method ** ) ckrealloc (obj->me,sizeof(Method *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_MethodTypeSet, returning FALSE");    
      return FALSE;  
      }  
    obj->me_maxlen = len;    
    return TRUE; 
}    


/* Function:  add_me_MethodTypeSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [MethodTypeSet *]
 * Arg:        add [OWNER] Object to add to the list [Method *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_me_MethodTypeSet(MethodTypeSet * obj,Method * add) 
{
    if( obj->me_len >= obj->me_maxlen)   {  
      if( expand_me_MethodTypeSet(obj,obj->me_len + MethodTypeSetLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->me[obj->me_len++]=add;  
    return TRUE; 
}    


/* Function:  flush_me_MethodTypeSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [MethodTypeSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_me_MethodTypeSet(MethodTypeSet * obj) 
{
    int i;   


    for(i=0;i<obj->me_len;i++)   { /*for i over list length*/ 
      if( obj->me[i] != NULL)    {  
        free_Method(obj->me[i]); 
        obj->me[i] = NULL;   
        }  
      } /* end of for i over list length */ 


    obj->me_len = 0; 
    return i;    
}    


/* Function:  swap_ty_MethodTypeSet(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_ty_MethodTypeSet
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [Type   **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_ty_MethodTypeSet(Type   ** list,int i,int j)  
{
    Type   * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_ty_MethodTypeSet(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_ty_MethodTypeSet which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [Type   **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_ty_MethodTypeSet(Type   ** list,int left,int right,int (*comp)(Type   * ,Type   * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_ty_MethodTypeSet(list,left,(left+right)/2); 
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_ty_MethodTypeSet (list,++last,i);   
      }  
    swap_ty_MethodTypeSet (list,left,last);  
    qsort_ty_MethodTypeSet(list,left,last-1,comp);   
    qsort_ty_MethodTypeSet(list,last+1,right,comp);  
}    


/* Function:  sort_ty_MethodTypeSet(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_ty_MethodTypeSet
 *
 *
 * Arg:         obj [UNKN ] Object containing list [MethodTypeSet *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_ty_MethodTypeSet(MethodTypeSet * obj,int (*comp)(Type   *, Type   *)) 
{
    qsort_ty_MethodTypeSet(obj->ty,0,obj->ty_len-1,comp);    
    return;  
}    


/* Function:  expand_ty_MethodTypeSet(obj,len)
 *
 * Descrip:    Really an internal function for add_ty_MethodTypeSet
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [MethodTypeSet *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_ty_MethodTypeSet(MethodTypeSet * obj,int len) 
{


    if( obj->ty_maxlen > obj->ty_len )   {  
      warn("expand_MethodTypeSetty_ called with no need");   
      return TRUE;   
      }  


    if( (obj->ty = (Type   ** ) ckrealloc (obj->ty,sizeof(Type   *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_MethodTypeSet, returning FALSE");    
      return FALSE;  
      }  
    obj->ty_maxlen = len;    
    return TRUE; 
}    


/* Function:  add_ty_MethodTypeSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [MethodTypeSet *]
 * Arg:        add [OWNER] Object to add to the list [Type   *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_ty_MethodTypeSet(MethodTypeSet * obj,Type   * add) 
{
    if( obj->ty_len >= obj->ty_maxlen)   {  
      if( expand_ty_MethodTypeSet(obj,obj->ty_len + MethodTypeSetLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->ty[obj->ty_len++]=add;  
    return TRUE; 
}    


/* Function:  flush_ty_MethodTypeSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [MethodTypeSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_ty_MethodTypeSet(MethodTypeSet * obj) 
{
    int i;   


    for(i=0;i<obj->ty_len;i++)   { /*for i over list length*/ 
      if( obj->ty[i] != NULL)    {  
        free_Type(obj->ty[i]);   
        obj->ty[i] = NULL;   
        }  
      } /* end of for i over list length */ 


    obj->ty_len = 0; 
    return i;    
}    


/* Function:  swap_in_MethodTypeSet(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_in_MethodTypeSet
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [Input  **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_in_MethodTypeSet(Input  ** list,int i,int j)  
{
    Input  * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_in_MethodTypeSet(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_in_MethodTypeSet which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [Input  **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_in_MethodTypeSet(Input  ** list,int left,int right,int (*comp)(Input  * ,Input  * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_in_MethodTypeSet(list,left,(left+right)/2); 
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_in_MethodTypeSet (list,++last,i);   
      }  
    swap_in_MethodTypeSet (list,left,last);  
    qsort_in_MethodTypeSet(list,left,last-1,comp);   
    qsort_in_MethodTypeSet(list,last+1,right,comp);  
}    


/* Function:  sort_in_MethodTypeSet(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_in_MethodTypeSet
 *
 *
 * Arg:         obj [UNKN ] Object containing list [MethodTypeSet *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_in_MethodTypeSet(MethodTypeSet * obj,int (*comp)(Input  *, Input  *)) 
{
    qsort_in_MethodTypeSet(obj->in,0,obj->in_len-1,comp);    
    return;  
}    


/* Function:  expand_in_MethodTypeSet(obj,len)
 *
 * Descrip:    Really an internal function for add_in_MethodTypeSet
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [MethodTypeSet *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_in_MethodTypeSet(MethodTypeSet * obj,int len) 
{


    if( obj->in_maxlen > obj->in_len )   {  
      warn("expand_MethodTypeSetin_ called with no need");   
      return TRUE;   
      }  


    if( (obj->in = (Input  ** ) ckrealloc (obj->in,sizeof(Input  *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_MethodTypeSet, returning FALSE");    
      return FALSE;  
      }  
    obj->in_maxlen = len;    
    return TRUE; 
}    


/* Function:  add_in_MethodTypeSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [MethodTypeSet *]
 * Arg:        add [OWNER] Object to add to the list [Input  *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_in_MethodTypeSet(MethodTypeSet * obj,Input  * add) 
{
    if( obj->in_len >= obj->in_maxlen)   {  
      if( expand_in_MethodTypeSet(obj,obj->in_len + MethodTypeSetLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->in[obj->in_len++]=add;  
    return TRUE; 
}    


/* Function:  flush_in_MethodTypeSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [MethodTypeSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_in_MethodTypeSet(MethodTypeSet * obj) 
{
    int i;   


    for(i=0;i<obj->in_len;i++)   { /*for i over list length*/ 
      if( obj->in[i] != NULL)    {  
        free_Input(obj->in[i]);  
        obj->in[i] = NULL;   
        }  
      } /* end of for i over list length */ 


    obj->in_len = 0; 
    return i;    
}    


/* Function:  MethodTypeSet_alloc_std(void)
 *
 * Descrip:    Equivalent to MethodTypeSet_alloc_len(MethodTypeSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [MethodTypeSet *]
 *
 */
MethodTypeSet * MethodTypeSet_alloc_std(void) 
{
    return MethodTypeSet_alloc_len(MethodTypeSetLISTLENGTH); 
}    


/* Function:  MethodTypeSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [MethodTypeSet *]
 *
 */
MethodTypeSet * MethodTypeSet_alloc_len(int len) 
{
    MethodTypeSet * out;/* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = MethodTypeSet_alloc()) == NULL)    
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->me = (Method ** ) ckcalloc (len,sizeof(Method *))) == NULL) {  
      warn("Warning, ckcalloc failed in MethodTypeSet_alloc_len");   
      return NULL;   
      }  
    out->me_len = 0; 
    out->me_maxlen = len;    


    if((out->ty = (Type   ** ) ckcalloc (len,sizeof(Type   *))) == NULL) {  
      warn("Warning, ckcalloc failed in MethodTypeSet_alloc_len");   
      return NULL;   
      }  
    out->ty_len = 0; 
    out->ty_maxlen = len;    


    if((out->in = (Input  ** ) ckcalloc (len,sizeof(Input  *))) == NULL) {  
      warn("Warning, ckcalloc failed in MethodTypeSet_alloc_len");   
      return NULL;   
      }  
    out->in_len = 0; 
    out->in_maxlen = len;    


    return out;  
}    


/* Function:  hard_link_MethodTypeSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [MethodTypeSet *]
 *
 * Return [UNKN ]  Undocumented return value [MethodTypeSet *]
 *
 */
MethodTypeSet * hard_link_MethodTypeSet(MethodTypeSet * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a MethodTypeSet object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  MethodTypeSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [MethodTypeSet *]
 *
 */
MethodTypeSet * MethodTypeSet_alloc(void) 
{
    MethodTypeSet * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(MethodTypeSet *) ckalloc (sizeof(MethodTypeSet))) == NULL)  {  
      warn("MethodTypeSet_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->me = NULL;  
    out->me_len = out->me_maxlen = 0;    
    out->ty = NULL;  
    out->ty_len = out->ty_maxlen = 0;    
    out->in = NULL;  
    out->in_len = out->in_maxlen = 0;    


    return out;  
}    


/* Function:  free_MethodTypeSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [MethodTypeSet *]
 *
 * Return [UNKN ]  Undocumented return value [MethodTypeSet *]
 *
 */
MethodTypeSet * free_MethodTypeSet(MethodTypeSet * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a MethodTypeSet obj. Should be trappable"); 
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
    if( obj->me != NULL) {  
      for(i=0;i<obj->me_len;i++) {  
        if( obj->me[i] != NULL)  
          free_Method(obj->me[i]);   
        }  
      ckfree(obj->me);   
      }  
    if( obj->ty != NULL) {  
      for(i=0;i<obj->ty_len;i++) {  
        if( obj->ty[i] != NULL)  
          free_Type(obj->ty[i]); 
        }  
      ckfree(obj->ty);   
      }  
    if( obj->in != NULL) {  
      for(i=0;i<obj->in_len;i++) {  
        if( obj->in[i] != NULL)  
          free_Input(obj->in[i]);    
        }  
      ckfree(obj->in);   
      }  


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
