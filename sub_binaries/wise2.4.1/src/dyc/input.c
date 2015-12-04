#ifdef _cplusplus
extern "C" {
#endif

#include "input.h"


# line 34 "input.dy"
Input * read_Input_line(char * line,FILE * ifp)
{
  Input * out;
  InRequire * re;
  InDeclare * de;
  char buffer[MAXLINE];

  out = Input_alloc_std();

  while( fgets(buffer,MAXLINE,ifp) != NULL ) {
    if( strwhitestartcmp(buffer,"#",spacestr) == 0 )
      continue;

    if( strstartcmp(buffer,"declare") == 0 ) {
      de = InDeclare_from_line(buffer);
      if( de != NULL ) {
	add_decl_Input(out,de);
      } else {
	warn("Unreadable declaration line in input");
      }
    } else if ( strstartcmp(buffer,"require") == 0 ) {
      re = InRequire_from_line(buffer);
      if( re != NULL ) {
	add_req_Input(out,re);
      } else {
	warn("Unreadable require line in input");
      }
    } else if ( strstartcmp(buffer,"code") == 0 ) {
      if( out->code != NULL ) {
	warn("Two code segments in an input!");
	out->code = free_Code(out->code);
      }
      out->code = read_Code_line(buffer,ifp);
    } else {
      warn("Could not understand line %s in input specification",buffer);
    }
  }

  return out;
}

# line 75 "input.dy"
Code * read_Code_line(char * line,FILE * ifp)
{
  char buffer[MAXLINE];
  Code * out;

  out = Code_alloc_std();

  while( fgets(buffer,MAXLINE,ifp) != NULL ) {
    if( strstartcmp(buffer,"endcode") == 0 ) 
      break;
    
  }


  return out;
}

# line 92 "input.dy"
InRequire * InRequire_from_line(char * line)
{
  InRequire * out;
  int i;
  char ** base;
  char ** brk;


  base = brk = breakstring(line,spacestr);
  if( *brk == NULL || strcmp(*brk,"require") == 0 ) {
    warn("In require - no require line!");
  }
  for(i=1;i<5;i++) {
    if( brk[i] == NULL ) {
      warn("wrong number of arguments for require line [%s]",line);
      return NULL;
    }
  }

  out = InRequire_alloc();

  out->name = stringalloc(brk[1]);
  out->type = stringalloc(brk[2]);
  strip_quote_chars(brk[3],"\"");
  out->help = stringalloc(brk[3]);
  out->tag  = stringalloc(brk[4]);

  ckfree(base);
  return out;
}

# line 123 "input.dy"
InDeclare * InDeclare_from_line(char * line)
{
  InDeclare * out;
  char * runner;
  char * name;
  char * type;


  runner = strtok(line,spacestr);
  if( strcmp(runner,"declare") != 0 ) {
    warn("You don't have a declaration line in %s",line);
    return NULL;
  }
  name = strtok(NULL,spacestr);
  type = strtok(NULL,spacestr);

  if( type == NULL ) {
    warn("In reading the line %s, we have no type. format problem",line);
    return NULL;
  }

  out = InDeclare_alloc();
  out->name = stringalloc(name);
  out->type = stringalloc(type);
  
  return out;
}

      


# line 131 "input.c"
/* Function:  swap_Code(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_Code
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [char **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_Code(char ** list,int i,int j)  
{
    char * temp; 
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_Code(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_Code which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [char **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_Code(char ** list,int left,int right,int (*comp)(char * ,char * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_Code(list,left,(left+right)/2); 
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_Code (list,++last,i);   
      }  
    swap_Code (list,left,last);  
    qsort_Code(list,left,last-1,comp);   
    qsort_Code(list,last+1,right,comp);  
}    


/* Function:  sort_Code(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_Code
 *
 *
 * Arg:         obj [UNKN ] Object containing list [Code *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_Code(Code * obj,int (*comp)(char *, char *)) 
{
    qsort_Code(obj->lines,0,obj->len-1,comp);    
    return;  
}    


/* Function:  expand_Code(obj,len)
 *
 * Descrip:    Really an internal function for add_Code
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Code *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_Code(Code * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_Code called with no need");   
      return TRUE;   
      }  


    if( (obj->lines = (char ** ) ckrealloc (obj->lines,sizeof(char *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_Code, returning FALSE"); 
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_Code(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Code *]
 * Arg:        add [OWNER] Object to add to the list [char *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_Code(Code * obj,char * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_Code(obj,obj->len + CodeLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->lines[obj->len++]=add;  
    return TRUE; 
}    


/* Function:  flush_Code(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [Code *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_Code(Code * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->lines[i] != NULL) {  
        ckfree(obj->lines[i]);   
        obj->lines[i] = NULL;    
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  Code_alloc_std(void)
 *
 * Descrip:    Equivalent to Code_alloc_len(CodeLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Code *]
 *
 */
Code * Code_alloc_std(void) 
{
    return Code_alloc_len(CodeLISTLENGTH);   
}    


/* Function:  Code_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [Code *]
 *
 */
Code * Code_alloc_len(int len) 
{
    Code * out; /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = Code_alloc()) == NULL) 
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->lines = (char ** ) ckcalloc (len,sizeof(char *))) == NULL)  {  
      warn("Warning, ckcalloc failed in Code_alloc_len");    
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_Code(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Code *]
 *
 * Return [UNKN ]  Undocumented return value [Code *]
 *
 */
Code * hard_link_Code(Code * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a Code object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  Code_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Code *]
 *
 */
Code * Code_alloc(void) 
{
    Code * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(Code *) ckalloc (sizeof(Code))) == NULL)    {  
      warn("Code_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->lines = NULL;   
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_Code(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Code *]
 *
 * Return [UNKN ]  Undocumented return value [Code *]
 *
 */
Code * free_Code(Code * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a Code obj. Should be trappable");  
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
    if( obj->lines != NULL)  {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->lines[i] != NULL)   
          ckfree(obj->lines[i]); 
        }  
      ckfree(obj->lines);    
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_InDeclare(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [InDeclare *]
 *
 * Return [UNKN ]  Undocumented return value [InDeclare *]
 *
 */
InDeclare * hard_link_InDeclare(InDeclare * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a InDeclare object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  InDeclare_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [InDeclare *]
 *
 */
InDeclare * InDeclare_alloc(void) 
{
    InDeclare * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(InDeclare *) ckalloc (sizeof(InDeclare))) == NULL)  {  
      warn("InDeclare_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->name = NULL;    
    out->type = NULL;    


    return out;  
}    


/* Function:  free_InDeclare(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [InDeclare *]
 *
 * Return [UNKN ]  Undocumented return value [InDeclare *]
 *
 */
InDeclare * free_InDeclare(InDeclare * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a InDeclare obj. Should be trappable"); 
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


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_InRequire(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [InRequire *]
 *
 * Return [UNKN ]  Undocumented return value [InRequire *]
 *
 */
InRequire * hard_link_InRequire(InRequire * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a InRequire object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  InRequire_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [InRequire *]
 *
 */
InRequire * InRequire_alloc(void) 
{
    InRequire * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(InRequire *) ckalloc (sizeof(InRequire))) == NULL)  {  
      warn("InRequire_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->name = NULL;    
    out->type = NULL;    
    out->help = NULL;    
    out->tag = NULL; 


    return out;  
}    


/* Function:  free_InRequire(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [InRequire *]
 *
 * Return [UNKN ]  Undocumented return value [InRequire *]
 *
 */
InRequire * free_InRequire(InRequire * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a InRequire obj. Should be trappable"); 
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
    if( obj->help != NULL)   
      ckfree(obj->help);     
    if( obj->tag != NULL)    
      ckfree(obj->tag);  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_decl_Input(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_decl_Input
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [InDeclare **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_decl_Input(InDeclare ** list,int i,int j)  
{
    InDeclare * temp;    
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_decl_Input(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_decl_Input which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [InDeclare **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_decl_Input(InDeclare ** list,int left,int right,int (*comp)(InDeclare * ,InDeclare * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_decl_Input(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_decl_Input (list,++last,i); 
      }  
    swap_decl_Input (list,left,last);    
    qsort_decl_Input(list,left,last-1,comp); 
    qsort_decl_Input(list,last+1,right,comp);    
}    


/* Function:  sort_decl_Input(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_decl_Input
 *
 *
 * Arg:         obj [UNKN ] Object containing list [Input *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_decl_Input(Input * obj,int (*comp)(InDeclare *, InDeclare *)) 
{
    qsort_decl_Input(obj->decl,0,obj->decl_len-1,comp);  
    return;  
}    


/* Function:  expand_decl_Input(obj,len)
 *
 * Descrip:    Really an internal function for add_decl_Input
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Input *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_decl_Input(Input * obj,int len) 
{


    if( obj->decl_maxlen > obj->decl_len )   {  
      warn("expand_Inputdecl_ called with no need"); 
      return TRUE;   
      }  


    if( (obj->decl = (InDeclare ** ) ckrealloc (obj->decl,sizeof(InDeclare *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_Input, returning FALSE");    
      return FALSE;  
      }  
    obj->decl_maxlen = len;  
    return TRUE; 
}    


/* Function:  add_decl_Input(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Input *]
 * Arg:        add [OWNER] Object to add to the list [InDeclare *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_decl_Input(Input * obj,InDeclare * add) 
{
    if( obj->decl_len >= obj->decl_maxlen)   {  
      if( expand_decl_Input(obj,obj->decl_len + InputLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->decl[obj->decl_len++]=add;  
    return TRUE; 
}    


/* Function:  flush_decl_Input(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [Input *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_decl_Input(Input * obj) 
{
    int i;   


    for(i=0;i<obj->decl_len;i++) { /*for i over list length*/ 
      if( obj->decl[i] != NULL)  {  
        free_InDeclare(obj->decl[i]);    
        obj->decl[i] = NULL; 
        }  
      } /* end of for i over list length */ 


    obj->decl_len = 0;   
    return i;    
}    


/* Function:  swap_req_Input(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_req_Input
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [InRequire **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_req_Input(InRequire ** list,int i,int j)  
{
    InRequire * temp;    
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_req_Input(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_req_Input which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [InRequire **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_req_Input(InRequire ** list,int left,int right,int (*comp)(InRequire * ,InRequire * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_req_Input(list,left,(left+right)/2);    
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_req_Input (list,++last,i);  
      }  
    swap_req_Input (list,left,last); 
    qsort_req_Input(list,left,last-1,comp);  
    qsort_req_Input(list,last+1,right,comp); 
}    


/* Function:  sort_req_Input(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_req_Input
 *
 *
 * Arg:         obj [UNKN ] Object containing list [Input *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_req_Input(Input * obj,int (*comp)(InRequire *, InRequire *)) 
{
    qsort_req_Input(obj->req,0,obj->req_len-1,comp); 
    return;  
}    


/* Function:  expand_req_Input(obj,len)
 *
 * Descrip:    Really an internal function for add_req_Input
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Input *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_req_Input(Input * obj,int len) 
{


    if( obj->req_maxlen > obj->req_len )     {  
      warn("expand_Inputreq_ called with no need");  
      return TRUE;   
      }  


    if( (obj->req = (InRequire ** ) ckrealloc (obj->req,sizeof(InRequire *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_Input, returning FALSE");    
      return FALSE;  
      }  
    obj->req_maxlen = len;   
    return TRUE; 
}    


/* Function:  add_req_Input(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Input *]
 * Arg:        add [OWNER] Object to add to the list [InRequire *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_req_Input(Input * obj,InRequire * add) 
{
    if( obj->req_len >= obj->req_maxlen) {  
      if( expand_req_Input(obj,obj->req_len + InputLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->req[obj->req_len++]=add;    
    return TRUE; 
}    


/* Function:  flush_req_Input(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [Input *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_req_Input(Input * obj) 
{
    int i;   


    for(i=0;i<obj->req_len;i++)  { /*for i over list length*/ 
      if( obj->req[i] != NULL)   {  
        free_InRequire(obj->req[i]); 
        obj->req[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->req_len = 0;    
    return i;    
}    


/* Function:  Input_alloc_std(void)
 *
 * Descrip:    Equivalent to Input_alloc_len(InputLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Input *]
 *
 */
Input * Input_alloc_std(void) 
{
    return Input_alloc_len(InputLISTLENGTH); 
}    


/* Function:  Input_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [Input *]
 *
 */
Input * Input_alloc_len(int len) 
{
    Input * out;/* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = Input_alloc()) == NULL)    
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->decl = (InDeclare ** ) ckcalloc (len,sizeof(InDeclare *))) == NULL) {  
      warn("Warning, ckcalloc failed in Input_alloc_len");   
      return NULL;   
      }  
    out->decl_len = 0;   
    out->decl_maxlen = len;  


    if((out->req = (InRequire ** ) ckcalloc (len,sizeof(InRequire *))) == NULL)  {  
      warn("Warning, ckcalloc failed in Input_alloc_len");   
      return NULL;   
      }  
    out->req_len = 0;    
    out->req_maxlen = len;   


    return out;  
}    


/* Function:  hard_link_Input(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Input *]
 *
 * Return [UNKN ]  Undocumented return value [Input *]
 *
 */
Input * hard_link_Input(Input * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a Input object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  Input_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Input *]
 *
 */
Input * Input_alloc(void) 
{
    Input * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(Input *) ckalloc (sizeof(Input))) == NULL)  {  
      warn("Input_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->decl = NULL;    
    out->decl_len = out->decl_maxlen = 0;    
    out->req = NULL; 
    out->req_len = out->req_maxlen = 0;  
    out->code = NULL;    


    return out;  
}    


/* Function:  free_Input(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Input *]
 *
 * Return [UNKN ]  Undocumented return value [Input *]
 *
 */
Input * free_Input(Input * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a Input obj. Should be trappable"); 
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
    if( obj->decl != NULL)   {  
      for(i=0;i<obj->decl_len;i++)   {  
        if( obj->decl[i] != NULL)    
          free_InDeclare(obj->decl[i]);  
        }  
      ckfree(obj->decl); 
      }  
    if( obj->req != NULL)    {  
      for(i=0;i<obj->req_len;i++)    {  
        if( obj->req[i] != NULL) 
          free_InRequire(obj->req[i]);   
        }  
      ckfree(obj->req);  
      }  
    if( obj->code != NULL)   
      free_Code(obj->code);  


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
