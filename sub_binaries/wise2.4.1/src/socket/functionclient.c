#ifdef _cplusplus
extern "C" {
#endif
#include "functionclient.h"

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>

#include <errno.h>  /* errno */



/* Function:  free_FunctionProxyCoordinator(obj)
 *
 * Descrip:    provides specific destructor for FunctionProxyCoordinators
 *
 *
 * Arg:        obj [UNKN ] Undocumented argument [FunctionProxyCoordinator *]
 *
 * Return [UNKN ]  Undocumented return value [FunctionProxyCoordinator *]
 *
 */
# line 41 "functionclient.dy"
FunctionProxyCoordinator * free_FunctionProxyCoordinator(FunctionProxyCoordinator * obj)
{
  int i;
  int return_early = 0;    
  char buffer[1024];
  
  
  if( obj == NULL) {  
    warn("Attempting to free a NULL pointer to a FunctionProxy obj. Should be trappable"); 
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
  

  /* send FINISH to server */

  strcpy(buffer,"FINISH");
  if (write(obj->socket, buffer, 6) == -1) {
    fprintf(stderr, "write() failed\n");
    fatal(strerror(errno));
  }
  
  close(obj->socket);
  
  for(i=0;i<obj->len;i++) {
    if( obj->proxy[i] != NULL ) {
      free_FunctionProxy(obj->proxy[i]);
    }
  }
  
  ckfree(obj); 
  return NULL; 
}

# line 88 "functionclient.dy"
AnonymousObject * dispatch_FunctionProxy(FunctionProxyCoordinator * fpc,char * function_call,AnonymousObjectList * args)
{
  int i,j;
  AnonymousObject * ret;
  char buffer[1024];

  FILE * ifp;
  Wise2WriteStreamInterface * ws;
  Wise2ReadStreamInterface * rs;


  assert(fpc != NULL);
  assert(function_call != NULL);
  assert(args != NULL);

  for(i=0;i<fpc->len;i++) {
    if( strcmp(fpc->proxy[i]->transfer->name,function_call) == 0 ) {
      break;
    }
  }

  if( i >= fpc->len ) {
    warn("In FunctionProxyCoordinator attached to %s,%d, function call %s has not been registered",fpc->host,fpc->port,function_call);
    return NULL;
  } else {
    fprintf(stderr,"Function proxy is ok with %s\n",function_call);
  }

  /* ok, go into cycle */

    warn("Ewan has not protected this write from resend");
    if (write(fpc->socket, function_call, strlen(function_call)) == -1) {
	fprintf(stderr, "write() failed\n");
	fatal(strerror(errno));
    }

    if (read(fpc->socket, buffer, 1024) == -1) {
	fprintf(stderr, "read() failed\n");
	fatal(strerror(errno));
    }

  if( strncmp(buffer,"NOT FOUND",9) == 0 ) {
    warn("In FunctionProxyCoordinator attached to %s,%d, function call %s is not present on server",fpc->host,fpc->port,function_call);
    return NULL;
  } else if( strncmp(buffer,"SEND",4) == 0 ) {
    fprintf(stderr,"Function call %s was accepted, sending data\n",function_call);
  } else {
    warn("In FunctionProxyCoordinator attached to %s, %d, "
	 "function call %s triggered unknown response ('%s')",
	 fpc->host, fpc->port, function_call, buffer);
    return NULL;
  }

  ws = WriteStream_from_socket(fpc->socket);

  for(j=0;j<fpc->proxy[i]->transfer->len;j++) {
    (*fpc->proxy[i]->transfer->input[j]->object_to_stream)(args->anon[j]->obj,ws);
  }
  
  if ((ifp = fdopen(fpc->socket,"r")) == NULL) {
      fprintf(stderr, "fdopen() failed\n");
      fatal(strerror(errno));
  }
  rs = ReadStream_from_FILE(ifp);

  ret = AnonymousObject_alloc();
  ret->obj = (*fpc->proxy[i]->transfer->returned_type->stream_to_object)(rs);
  ret->free_object = fpc->proxy[i]->transfer->returned_type->untyped_free_object;


  return ret;

}


# line 163 "functionclient.dy"
FunctionProxy * new_FunctionProxy(TransferedFunctionCall * t)
{
  FunctionProxy * out;

  out = FunctionProxy_alloc();
  out->transfer = t;

  return out;
}


# line 174 "functionclient.dy"
FunctionProxyCoordinator * new_FunctionProxyCoordinator(char * host,int port)
{
  FunctionProxyCoordinator * out;

  struct sockaddr_in server;
  struct hostent * hp;

  out = FunctionProxyCoordinator_alloc_std();

  out->host = stringalloc(host);
  out->port = port;

  out->socket = socket(AF_INET, SOCK_STREAM, 0);

  server.sin_family = AF_INET;
  hp = gethostbyname(host);
  bcopy( hp->h_addr, &server.sin_addr, hp->h_length);
  server.sin_port = htons(port);

  connect(out->socket, &server, sizeof(server));
  

  return out;
}



# line 189 "functionclient.c"
/* Function:  hard_link_FunctionProxy(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [FunctionProxy *]
 *
 * Return [UNKN ]  Undocumented return value [FunctionProxy *]
 *
 */
FunctionProxy * hard_link_FunctionProxy(FunctionProxy * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a FunctionProxy object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  FunctionProxy_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [FunctionProxy *]
 *
 */
FunctionProxy * FunctionProxy_alloc(void) 
{
    FunctionProxy * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(FunctionProxy *) ckalloc (sizeof(FunctionProxy))) == NULL)  {  
      warn("FunctionProxy_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->transfer = NULL;    


    return out;  
}    


/* Function:  free_FunctionProxy(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [FunctionProxy *]
 *
 * Return [UNKN ]  Undocumented return value [FunctionProxy *]
 *
 */
FunctionProxy * free_FunctionProxy(FunctionProxy * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a FunctionProxy obj. Should be trappable"); 
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
    if( obj->transfer != NULL)   
      free_TransferedFunctionCall(obj->transfer);    


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_FunctionProxyCoordinator(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_FunctionProxyCoordinator
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [FunctionProxy **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_FunctionProxyCoordinator(FunctionProxy ** list,int i,int j)  
{
    FunctionProxy * temp;    
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_FunctionProxyCoordinator(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_FunctionProxyCoordinator which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [FunctionProxy **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_FunctionProxyCoordinator(FunctionProxy ** list,int left,int right,int (*comp)(FunctionProxy * ,FunctionProxy * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_FunctionProxyCoordinator(list,left,(left+right)/2); 
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_FunctionProxyCoordinator (list,++last,i);   
      }  
    swap_FunctionProxyCoordinator (list,left,last);  
    qsort_FunctionProxyCoordinator(list,left,last-1,comp);   
    qsort_FunctionProxyCoordinator(list,last+1,right,comp);  
}    


/* Function:  sort_FunctionProxyCoordinator(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_FunctionProxyCoordinator
 *
 *
 * Arg:         obj [UNKN ] Object containing list [FunctionProxyCoordinator *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_FunctionProxyCoordinator(FunctionProxyCoordinator * obj,int (*comp)(FunctionProxy *, FunctionProxy *)) 
{
    qsort_FunctionProxyCoordinator(obj->proxy,0,obj->len-1,comp);    
    return;  
}    


/* Function:  expand_FunctionProxyCoordinator(obj,len)
 *
 * Descrip:    Really an internal function for add_FunctionProxyCoordinator
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [FunctionProxyCoordinator *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_FunctionProxyCoordinator(FunctionProxyCoordinator * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_FunctionProxyCoordinator called with no need");   
      return TRUE;   
      }  


    if( (obj->proxy = (FunctionProxy ** ) ckrealloc (obj->proxy,sizeof(FunctionProxy *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_FunctionProxyCoordinator, returning FALSE"); 
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_FunctionProxyCoordinator(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [FunctionProxyCoordinator *]
 * Arg:        add [OWNER] Object to add to the list [FunctionProxy *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_FunctionProxyCoordinator(FunctionProxyCoordinator * obj,FunctionProxy * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_FunctionProxyCoordinator(obj,obj->len + FunctionProxyCoordinatorLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->proxy[obj->len++]=add;  
    return TRUE; 
}    


/* Function:  flush_FunctionProxyCoordinator(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [FunctionProxyCoordinator *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_FunctionProxyCoordinator(FunctionProxyCoordinator * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->proxy[i] != NULL) {  
        free_FunctionProxy(obj->proxy[i]);   
        obj->proxy[i] = NULL;    
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  FunctionProxyCoordinator_alloc_std(void)
 *
 * Descrip:    Equivalent to FunctionProxyCoordinator_alloc_len(FunctionProxyCoordinatorLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [FunctionProxyCoordinator *]
 *
 */
FunctionProxyCoordinator * FunctionProxyCoordinator_alloc_std(void) 
{
    return FunctionProxyCoordinator_alloc_len(FunctionProxyCoordinatorLISTLENGTH);   
}    


/* Function:  FunctionProxyCoordinator_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [FunctionProxyCoordinator *]
 *
 */
FunctionProxyCoordinator * FunctionProxyCoordinator_alloc_len(int len) 
{
    FunctionProxyCoordinator * out; /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = FunctionProxyCoordinator_alloc()) == NULL) 
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->proxy = (FunctionProxy ** ) ckcalloc (len,sizeof(FunctionProxy *))) == NULL)    {  
      warn("Warning, ckcalloc failed in FunctionProxyCoordinator_alloc_len");    
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_FunctionProxyCoordinator(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [FunctionProxyCoordinator *]
 *
 * Return [UNKN ]  Undocumented return value [FunctionProxyCoordinator *]
 *
 */
FunctionProxyCoordinator * hard_link_FunctionProxyCoordinator(FunctionProxyCoordinator * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a FunctionProxyCoordinator object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  FunctionProxyCoordinator_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [FunctionProxyCoordinator *]
 *
 */
FunctionProxyCoordinator * FunctionProxyCoordinator_alloc(void) 
{
    FunctionProxyCoordinator * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(FunctionProxyCoordinator *) ckalloc (sizeof(FunctionProxyCoordinator))) == NULL)    {  
      warn("FunctionProxyCoordinator_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->proxy = NULL;   
    out->len = out->maxlen = 0;  
    out->socket = 0; 
    out->host = NULL;    
    out->port = 0;   


    return out;  
}    



#ifdef _cplusplus
}
#endif
