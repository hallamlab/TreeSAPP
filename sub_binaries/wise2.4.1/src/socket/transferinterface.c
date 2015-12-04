#ifdef _cplusplus
extern "C" {
#endif
#include "transferinterface.h"




# line 29 "transferinterface.dy"
TransferedFunctionCall * test_stringcat_TransferedFunctionCall(void)
{
  TransferedFunctionCall * out;

  out = TransferedFunctionCall_alloc_std();

  out->name = stringalloc("stringcat");
  out->returned_type = new_string_Marshaller();

  add_TransferedFunctionCall(out,new_string_Marshaller());
  add_TransferedFunctionCall(out,new_string_Marshaller());

  return out;
}

# line 44 "transferinterface.dy"
TransferedObjectMarshaller * new_string_Marshaller(void)
{
  TransferedObjectMarshaller * out;

  out = TransferedObjectMarshaller_alloc();
  out->object_to_stream = object_to_stream_string;
  out->stream_to_object = stream_to_object_string;
  return out;
}

# line 54 "transferinterface.dy"
void object_to_stream_string(void * obj,Wise2WriteStreamInterface* write)
{
  char * string = (char *) obj;

  fprintf(stderr,"Writing to buffer with string [%s]\n",string);
  (*write->write_buffer)(write->handle,string);
  (*write->write_buffer)(write->handle,"\n//\n");
}

# line 63 "transferinterface.dy"
void * stream_to_object_string(Wise2ReadStreamInterface * read)
{
  char buffer[1024];
  char * out;
  int i;

  fprintf(stderr,"reading string...\n");

  WISE2_READ_BUFFER(buffer,1024,read);

  fprintf(stderr,"read buffer %s...\n",buffer);
  for(i=0;i<1024;i++) {
    if( buffer[i] == '\n' || buffer[i] == '\r' || buffer[i] == '\0' ) {
      break;
    }
  }

  buffer[i] = '\0';
  out = stringalloc(buffer);

  WISE2_READ_BUFFER(buffer,1024,read);

  fprintf(stderr,"discarding buffer... %s, have string %s\n",buffer,out);

  return out;
}


# line 74 "transferinterface.c"
/* Function:  hard_link_TransferedObjectMarshaller(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [TransferedObjectMarshaller *]
 *
 * Return [UNKN ]  Undocumented return value [TransferedObjectMarshaller *]
 *
 */
TransferedObjectMarshaller * hard_link_TransferedObjectMarshaller(TransferedObjectMarshaller * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a TransferedObjectMarshaller object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  TransferedObjectMarshaller_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TransferedObjectMarshaller *]
 *
 */
TransferedObjectMarshaller * TransferedObjectMarshaller_alloc(void) 
{
    TransferedObjectMarshaller * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(TransferedObjectMarshaller *) ckalloc (sizeof(TransferedObjectMarshaller))) == NULL)    {  
      warn("TransferedObjectMarshaller_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->object_to_stream = NULL;    
    out->stream_to_object = NULL;    
    out->untyped_free_object = NULL; 


    return out;  
}    


/* Function:  free_TransferedObjectMarshaller(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [TransferedObjectMarshaller *]
 *
 * Return [UNKN ]  Undocumented return value [TransferedObjectMarshaller *]
 *
 */
TransferedObjectMarshaller * free_TransferedObjectMarshaller(TransferedObjectMarshaller * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a TransferedObjectMarshaller obj. Should be trappable");    
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
    /* obj->object_to_stream is a function pointer */ 
    /* obj->stream_to_object is a function pointer */ 
    /* obj->object_type is linked in */ 
    /* obj->untyped_free_object is a function pointer */ 


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_TransferedFunctionCall(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_TransferedFunctionCall
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [TransferedObjectMarshaller **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_TransferedFunctionCall(TransferedObjectMarshaller ** list,int i,int j)  
{
    TransferedObjectMarshaller * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_TransferedFunctionCall(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_TransferedFunctionCall which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [TransferedObjectMarshaller **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_TransferedFunctionCall(TransferedObjectMarshaller ** list,int left,int right,int (*comp)(TransferedObjectMarshaller * ,TransferedObjectMarshaller * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_TransferedFunctionCall(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_TransferedFunctionCall (list,++last,i); 
      }  
    swap_TransferedFunctionCall (list,left,last);    
    qsort_TransferedFunctionCall(list,left,last-1,comp); 
    qsort_TransferedFunctionCall(list,last+1,right,comp);    
}    


/* Function:  sort_TransferedFunctionCall(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_TransferedFunctionCall
 *
 *
 * Arg:         obj [UNKN ] Object containing list [TransferedFunctionCall *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_TransferedFunctionCall(TransferedFunctionCall * obj,int (*comp)(TransferedObjectMarshaller *, TransferedObjectMarshaller *)) 
{
    qsort_TransferedFunctionCall(obj->input,0,obj->len-1,comp);  
    return;  
}    


/* Function:  expand_TransferedFunctionCall(obj,len)
 *
 * Descrip:    Really an internal function for add_TransferedFunctionCall
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [TransferedFunctionCall *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_TransferedFunctionCall(TransferedFunctionCall * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_TransferedFunctionCall called with no need"); 
      return TRUE;   
      }  


    if( (obj->input = (TransferedObjectMarshaller ** ) ckrealloc (obj->input,sizeof(TransferedObjectMarshaller *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_TransferedFunctionCall, returning FALSE");   
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_TransferedFunctionCall(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [TransferedFunctionCall *]
 * Arg:        add [OWNER] Object to add to the list [TransferedObjectMarshaller *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_TransferedFunctionCall(TransferedFunctionCall * obj,TransferedObjectMarshaller * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_TransferedFunctionCall(obj,obj->len + TransferedFunctionCallLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->input[obj->len++]=add;  
    return TRUE; 
}    


/* Function:  flush_TransferedFunctionCall(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [TransferedFunctionCall *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_TransferedFunctionCall(TransferedFunctionCall * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->input[i] != NULL) {  
        free_TransferedObjectMarshaller(obj->input[i]);  
        obj->input[i] = NULL;    
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  TransferedFunctionCall_alloc_std(void)
 *
 * Descrip:    Equivalent to TransferedFunctionCall_alloc_len(TransferedFunctionCallLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TransferedFunctionCall *]
 *
 */
TransferedFunctionCall * TransferedFunctionCall_alloc_std(void) 
{
    return TransferedFunctionCall_alloc_len(TransferedFunctionCallLISTLENGTH);   
}    


/* Function:  TransferedFunctionCall_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [TransferedFunctionCall *]
 *
 */
TransferedFunctionCall * TransferedFunctionCall_alloc_len(int len) 
{
    TransferedFunctionCall * out;   /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = TransferedFunctionCall_alloc()) == NULL)   
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->input = (TransferedObjectMarshaller ** ) ckcalloc (len,sizeof(TransferedObjectMarshaller *))) == NULL)  {  
      warn("Warning, ckcalloc failed in TransferedFunctionCall_alloc_len");  
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_TransferedFunctionCall(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [TransferedFunctionCall *]
 *
 * Return [UNKN ]  Undocumented return value [TransferedFunctionCall *]
 *
 */
TransferedFunctionCall * hard_link_TransferedFunctionCall(TransferedFunctionCall * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a TransferedFunctionCall object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  TransferedFunctionCall_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TransferedFunctionCall *]
 *
 */
TransferedFunctionCall * TransferedFunctionCall_alloc(void) 
{
    TransferedFunctionCall * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(TransferedFunctionCall *) ckalloc (sizeof(TransferedFunctionCall))) == NULL)    {  
      warn("TransferedFunctionCall_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->name = NULL;    
    out->returned_type = NULL;   
    out->input = NULL;   
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_TransferedFunctionCall(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [TransferedFunctionCall *]
 *
 * Return [UNKN ]  Undocumented return value [TransferedFunctionCall *]
 *
 */
TransferedFunctionCall * free_TransferedFunctionCall(TransferedFunctionCall * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a TransferedFunctionCall obj. Should be trappable");    
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
    if( obj->returned_type != NULL)  
      free_TransferedObjectMarshaller(obj->returned_type);   
    if( obj->input != NULL)  {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->input[i] != NULL)   
          free_TransferedObjectMarshaller(obj->input[i]);    
        }  
      ckfree(obj->input);    
      }  


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
