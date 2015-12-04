#ifdef _cplusplus
extern "C" {
#endif
#include "modulefunc.h"

# line 54 "modulefunc.dy"
ModuleFunction * get_ModuleFunction_from_name(ModuleFunctionList * mfl,char * name)
{
  register int i;

  for(i=0;i<mfl->len;i++)
    if( strcmp(mfl->mf[i]->name,name) == 0 )
      return mfl->mf[i];

  return NULL;
}

# line 65 "modulefunc.dy"
ModuleFunction * new_ModuleFunction(ModuleFunctionList * mfl,char * name) 
{
  ModuleFunction * out;


  out = ModuleFunction_alloc();
  if( out == NULL )
    return out;

  out->name = stringalloc(name);
  add_ModuleFunctionList(mfl,out);


  return out;
}

# line 81 "modulefunc.dy"
void show_ModuleFunctionList(ModuleFunctionList * mfl,FILE * ofp)
{
  register int i;

  for(i=0;i<mfl->len;i++) 
    show_ModuleFunction(mfl->mf[i],ofp);

}

# line 90 "modulefunc.dy"
void show_ModuleFunction(ModuleFunction * mf,FILE * ofp)
{
  fprintf(ofp,"Module %s Constructor: %s Deconstructor %s\n",mf->name,
mf->has_cons == TRUE ? "yes" : " no", mf->has_decons == TRUE ? "yes" : " no");
}


# line 97 "modulefunc.dy"
char * parse_and_get_module_name_from_func(char * line,boolean isalloc)
{
  char * name;
  char * next;
  char * func;
  char buffer[128]; /** max function name! **/

  name = strtok(line," \t(");
  if( name == NULL ) {
    warn("Cannot even get first name from line [%s] in parse module_name_alloc",line);
    return NULL;
  }

  if( *(name + strlen(name) - 1) == '*' )
    next = name + strlen(name) -1;
  else {
    next = strtok(NULL," \t(");
    if ( next == NULL ) {
      warn("Cannot get pointer ref from line [%s] in parse module_name_alloc [name %s]",line,name);
      return NULL;
    }
  }

  func = strtok(NULL," \t(");
  if( name == NULL ) {
    warn("Cannot get function from line [%s] in parse module_name_alloc [name %s]",line,name);
    return NULL;
  }

  if( strlen(next) > 1 || *next != '*' ) {
    warn("In parse_module_name, the pointer string [%s] was invalid for name [%s]",next,name);
    return NULL;
  }

  if( isalloc == TRUE ) {
    sprintf(buffer,"%s_alloc",name);
    if( strcmp(buffer,func) != 0 ) {
      warn("In parse_module_name, the function [%s] did not match the type-proto [%s]",func,buffer);
      return NULL;
    }
  }
  else {
    sprintf(buffer,"free_%s",name);
    if( strcmp(buffer,func) != 0 ) {
      warn("In parse_module_name, the function [%s] did not match the type-proto [%s]",func,buffer);
      return NULL;
    }
  }

  return stringalloc(name);
}

  



# line 109 "modulefunc.c"
/* Function:  hard_link_ModuleFunction(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ModuleFunction *]
 *
 * Return [UNKN ]  Undocumented return value [ModuleFunction *]
 *
 */
ModuleFunction * hard_link_ModuleFunction(ModuleFunction * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a ModuleFunction object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  ModuleFunction_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ModuleFunction *]
 *
 */
ModuleFunction * ModuleFunction_alloc(void) 
{
    ModuleFunction * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(ModuleFunction *) ckalloc (sizeof(ModuleFunction))) == NULL)    {  
      warn("ModuleFunction_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->name = NULL;    
    out->has_cons = FALSE;   
    out->has_decons = FALSE; 
    out->has_copy = FALSE;   


    return out;  
}    


/* Function:  free_ModuleFunction(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ModuleFunction *]
 *
 * Return [UNKN ]  Undocumented return value [ModuleFunction *]
 *
 */
ModuleFunction * free_ModuleFunction(ModuleFunction * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a ModuleFunction obj. Should be trappable");    
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


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_ModuleFunctionList(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_ModuleFunctionList
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [ModuleFunction **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_ModuleFunctionList(ModuleFunction ** list,int i,int j)  
{
    ModuleFunction * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_ModuleFunctionList(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_ModuleFunctionList which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [ModuleFunction **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_ModuleFunctionList(ModuleFunction ** list,int left,int right,int (*comp)(ModuleFunction * ,ModuleFunction * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_ModuleFunctionList(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_ModuleFunctionList (list,++last,i); 
      }  
    swap_ModuleFunctionList (list,left,last);    
    qsort_ModuleFunctionList(list,left,last-1,comp); 
    qsort_ModuleFunctionList(list,last+1,right,comp);    
}    


/* Function:  sort_ModuleFunctionList(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_ModuleFunctionList
 *
 *
 * Arg:         obj [UNKN ] Object containing list [ModuleFunctionList *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_ModuleFunctionList(ModuleFunctionList * obj,int (*comp)(ModuleFunction *, ModuleFunction *)) 
{
    qsort_ModuleFunctionList(obj->mf,0,obj->len-1,comp); 
    return;  
}    


/* Function:  expand_ModuleFunctionList(obj,len)
 *
 * Descrip:    Really an internal function for add_ModuleFunctionList
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [ModuleFunctionList *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_ModuleFunctionList(ModuleFunctionList * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_ModuleFunctionList called with no need"); 
      return TRUE;   
      }  


    if( (obj->mf = (ModuleFunction ** ) ckrealloc (obj->mf,sizeof(ModuleFunction *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_ModuleFunctionList, returning FALSE");   
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_ModuleFunctionList(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [ModuleFunctionList *]
 * Arg:        add [OWNER] Object to add to the list [ModuleFunction *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_ModuleFunctionList(ModuleFunctionList * obj,ModuleFunction * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_ModuleFunctionList(obj,obj->len + ModuleFunctionListLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->mf[obj->len++]=add; 
    return TRUE; 
}    


/* Function:  flush_ModuleFunctionList(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [ModuleFunctionList *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_ModuleFunctionList(ModuleFunctionList * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->mf[i] != NULL)    {  
        free_ModuleFunction(obj->mf[i]); 
        obj->mf[i] = NULL;   
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  ModuleFunctionList_alloc_std(void)
 *
 * Descrip:    Equivalent to ModuleFunctionList_alloc_len(ModuleFunctionListLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ModuleFunctionList *]
 *
 */
ModuleFunctionList * ModuleFunctionList_alloc_std(void) 
{
    return ModuleFunctionList_alloc_len(ModuleFunctionListLISTLENGTH);   
}    


/* Function:  ModuleFunctionList_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [ModuleFunctionList *]
 *
 */
ModuleFunctionList * ModuleFunctionList_alloc_len(int len) 
{
    ModuleFunctionList * out;   /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = ModuleFunctionList_alloc()) == NULL)   
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->mf = (ModuleFunction ** ) ckcalloc (len,sizeof(ModuleFunction *))) == NULL) {  
      warn("Warning, ckcalloc failed in ModuleFunctionList_alloc_len");  
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_ModuleFunctionList(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ModuleFunctionList *]
 *
 * Return [UNKN ]  Undocumented return value [ModuleFunctionList *]
 *
 */
ModuleFunctionList * hard_link_ModuleFunctionList(ModuleFunctionList * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a ModuleFunctionList object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  ModuleFunctionList_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ModuleFunctionList *]
 *
 */
ModuleFunctionList * ModuleFunctionList_alloc(void) 
{
    ModuleFunctionList * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(ModuleFunctionList *) ckalloc (sizeof(ModuleFunctionList))) == NULL)    {  
      warn("ModuleFunctionList_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->mf = NULL;  
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_ModuleFunctionList(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ModuleFunctionList *]
 *
 * Return [UNKN ]  Undocumented return value [ModuleFunctionList *]
 *
 */
ModuleFunctionList * free_ModuleFunctionList(ModuleFunctionList * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a ModuleFunctionList obj. Should be trappable");    
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
    if( obj->mf != NULL) {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->mf[i] != NULL)  
          free_ModuleFunction(obj->mf[i]);   
        }  
      ckfree(obj->mf);   
      }  


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
