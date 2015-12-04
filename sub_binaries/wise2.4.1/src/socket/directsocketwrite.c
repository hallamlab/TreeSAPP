#ifdef _cplusplus
extern "C" {
#endif
#include "directsocketwrite.h"


# line 17 "directsocketwrite.dy"
void write_buffer_direct_socket_server_impl(void * handle,char * string)
{
  DirectSocketWrite * dsw = (DirectSocketWrite*) handle;

  write_DirectSocketWrite(dsw,string);
}

# line 24 "directsocketwrite.dy"
void write_bufferf_direct_socket_server_impl(void * handle,char * format,...)
{
  DirectSocketWrite * dsw = (DirectSocketWrite*) handle;
  char buffer[1024];
  va_list ap;


  va_start(ap,format);
  vsprintf(buffer,format,ap);	

  write_DirectSocketWrite(dsw,buffer);

  return;
}

# line 39 "directsocketwrite.dy"
void close_and_free_direct_socket_server(void * handle)
{
  DirectSocketWrite * dsw = (DirectSocketWrite*) handle;

  free_DirectSocketWrite(dsw);
  return;
}


# line 48 "directsocketwrite.dy"
Wise2WriteStreamInterface * WriteStream_from_socket(int socket)
{
  Wise2WriteStreamInterface * out;
  DirectSocketWrite * dsw;

  dsw = DirectSocketWrite_alloc();
  dsw->socket = socket;

  out = malloc(sizeof(Wise2WriteStreamInterface));

  out->write_buffer  = write_buffer_direct_socket_server_impl;
  out->write_bufferf = write_bufferf_direct_socket_server_impl;
  out->close_and_free_handle = close_and_free_direct_socket_server;

  out->handle = (void *) dsw;

  return out;
}



# line 69 "directsocketwrite.dy"
boolean write_DirectSocketWrite(DirectSocketWrite * dsw,char * string)
{
  int curr;
  int length;
  int ret;

  length = strlen(string);

  curr = 0;

  while( curr < length ) {
    ret = write(dsw->socket,string+curr,length-curr);
    if( ret == -1 ) {
      /* warn("Error sending with %d\n",errno); */
      break;
    }
    if( ret+curr >= length ) {
      break;
    }
    curr += ret;
  }
  
  return TRUE;
}


# line 89 "directsocketwrite.c"
/* Function:  hard_link_DirectSocketWrite(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DirectSocketWrite *]
 *
 * Return [UNKN ]  Undocumented return value [DirectSocketWrite *]
 *
 */
DirectSocketWrite * hard_link_DirectSocketWrite(DirectSocketWrite * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a DirectSocketWrite object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  DirectSocketWrite_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DirectSocketWrite *]
 *
 */
DirectSocketWrite * DirectSocketWrite_alloc(void) 
{
    DirectSocketWrite * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(DirectSocketWrite *) ckalloc (sizeof(DirectSocketWrite))) == NULL)  {  
      warn("DirectSocketWrite_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->socket = 0; 


    return out;  
}    


/* Function:  free_DirectSocketWrite(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DirectSocketWrite *]
 *
 * Return [UNKN ]  Undocumented return value [DirectSocketWrite *]
 *
 */
DirectSocketWrite * free_DirectSocketWrite(DirectSocketWrite * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a DirectSocketWrite obj. Should be trappable"); 
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


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
