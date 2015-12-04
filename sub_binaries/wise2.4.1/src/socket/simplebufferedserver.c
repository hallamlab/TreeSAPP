#ifdef _cplusplus
extern "C" {
#endif
#include "simplebufferedserver.h"

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <stdio.h>
#include <stdlib.h>


# line 33 "simplebufferedserver.dy"
void write_buffer_simple_buffered_server_impl(void * handle,char * string)
{
  SimpleBufferedServer * sbs = (SimpleBufferedServer*) handle;

  add_string_SimpleBufferedServer(sbs,string);
}

# line 40 "simplebufferedserver.dy"
void write_bufferf_simple_buffered_server_impl(void * handle,char * format,...)
{
  SimpleBufferedServer * sbs = (SimpleBufferedServer*) handle;
  char buffer[1024];
  va_list ap;


  va_start(ap,format);
  vsprintf(buffer,format,ap);	

  add_string_SimpleBufferedServer(sbs,buffer);

  return;
}

# line 55 "simplebufferedserver.dy"
void close_and_free_simple_buffered_server(void * handle)
{
  /* no-op, as this is not the owner */
  return;
}

# line 61 "simplebufferedserver.dy"
void add_string_SimpleBufferedServer(SimpleBufferedServer * sbs,char * unallocated_string)
{
  add_SimpleBufferedServer(sbs,stringalloc(unallocated_string));
}


# line 67 "simplebufferedserver.dy"
Wise2WriteStreamInterface * WriteStream_from_SimpleBufferedServer(SimpleBufferedServer * sbs)
{
  Wise2WriteStreamInterface * out;

  out = malloc(sizeof(Wise2WriteStreamInterface));

  out->write_buffer  = write_buffer_simple_buffered_server_impl;
  out->write_bufferf = write_bufferf_simple_buffered_server_impl;
  out->close_and_free_handle = close_and_free_simple_buffered_server;

  out->handle = (void *) sbs;

  return out;
}

# line 82 "simplebufferedserver.dy"
SimpleBufferedServer * new_SimpleBufferedServer(int port)
{
  struct sockaddr_in server;
  struct sockaddr * name;
  SimpleBufferedServer * out;
  int yes = 1;

  out = SimpleBufferedServer_alloc_std();

  out->socket = socket(AF_INET, SOCK_STREAM, 0);
  
  server.sin_family = AF_INET;
  server.sin_addr.s_addr = INADDR_ANY;
  server.sin_port = htons(port);

  setsockopt(out->socket, SOL_SOCKET, SO_REUSEADDR, &yes, sizeof yes);
  /* If setsockopt() fails, we're not too unhappy, so no check */

  if (bind(out->socket, (struct sockaddr *)&server, sizeof(server))) {
      /* Probably tried to bind to a port that was already in use... */
      fatal(strerror(errno));
  }

  return out;
}

# line 108 "simplebufferedserver.dy"
char * free_BufferedString(char * str)
{
  free(str);
}

# line 113 "simplebufferedserver.dy"
void listen_and_respond_SimpleBufferedServer(SimpleBufferedServer * sbs)
{
  int i;
  char buf[1024];
  int ret;
  int len;
  int curr;
  int ms;

  int new_socket;

  listen(sbs->socket, 5);
  
  
  while( 1 ) {
    fprintf(stderr,"blocking at accept...\n");
    ms = accept(sbs->socket, 0, 0);
    fprintf(stderr,"Accepted socket...\n");
    read(ms,buf,sizeof(buf));
    fprintf(stderr,"Read buffer...\n");
    if( strstartcmp(buf,"<WISE2-START-SIMPLE-STREAM>") == 0 ) {
      for(i=0;i<sbs->len;i++) {
	fprintf(stderr,"Going to send position %d\n",i);
	len = strlen(sbs->buffer[i]);
	curr = 0;
	while( curr < len ) {
	  fprintf(stderr,"Going to send at position %d\n",curr);
	  ret = write(ms,sbs->buffer[i]+curr,len-curr);
	  if( ret == -1 ) {
	    warn("error sending with %d\n",errno);
	  }	
	  fprintf(stderr,"...sent at position %d, sent %d bytes [%.6s]\n",curr,ret,sbs->buffer[i]);
	  if( ret+curr >= len ) {
	    break;
	  }
	  curr += ret;
	  sleep(2);
	}
      }
      buf[0] = '\0';
      write(ms,buf,1);
      close(ms);
    } else {
      fatal("Bad client initialisation string...%s\n",buf);
    }
  }

}

# line 150 "simplebufferedserver.c"
/* Function:  swap_SimpleBufferedServer(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_SimpleBufferedServer
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [BufferedString **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_SimpleBufferedServer(BufferedString ** list,int i,int j)  
{
    BufferedString * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_SimpleBufferedServer(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_SimpleBufferedServer which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [BufferedString **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_SimpleBufferedServer(BufferedString ** list,int left,int right,int (*comp)(BufferedString * ,BufferedString * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_SimpleBufferedServer(list,left,(left+right)/2); 
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_SimpleBufferedServer (list,++last,i);   
      }  
    swap_SimpleBufferedServer (list,left,last);  
    qsort_SimpleBufferedServer(list,left,last-1,comp);   
    qsort_SimpleBufferedServer(list,last+1,right,comp);  
}    


/* Function:  sort_SimpleBufferedServer(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_SimpleBufferedServer
 *
 *
 * Arg:         obj [UNKN ] Object containing list [SimpleBufferedServer *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_SimpleBufferedServer(SimpleBufferedServer * obj,int (*comp)(BufferedString *, BufferedString *)) 
{
    qsort_SimpleBufferedServer(obj->buffer,0,obj->len-1,comp);   
    return;  
}    


/* Function:  expand_SimpleBufferedServer(obj,len)
 *
 * Descrip:    Really an internal function for add_SimpleBufferedServer
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [SimpleBufferedServer *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_SimpleBufferedServer(SimpleBufferedServer * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_SimpleBufferedServer called with no need");   
      return TRUE;   
      }  


    if( (obj->buffer = (BufferedString ** ) ckrealloc (obj->buffer,sizeof(BufferedString *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_SimpleBufferedServer, returning FALSE"); 
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_SimpleBufferedServer(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [SimpleBufferedServer *]
 * Arg:        add [OWNER] Object to add to the list [BufferedString *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_SimpleBufferedServer(SimpleBufferedServer * obj,BufferedString * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_SimpleBufferedServer(obj,obj->len + SimpleBufferedServerLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->buffer[obj->len++]=add; 
    return TRUE; 
}    


/* Function:  flush_SimpleBufferedServer(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [SimpleBufferedServer *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_SimpleBufferedServer(SimpleBufferedServer * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->buffer[i] != NULL)    {  
        free_BufferedString(obj->buffer[i]); 
        obj->buffer[i] = NULL;   
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  SimpleBufferedServer_alloc_std(void)
 *
 * Descrip:    Equivalent to SimpleBufferedServer_alloc_len(SimpleBufferedServerLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SimpleBufferedServer *]
 *
 */
SimpleBufferedServer * SimpleBufferedServer_alloc_std(void) 
{
    return SimpleBufferedServer_alloc_len(SimpleBufferedServerLISTLENGTH);   
}    


/* Function:  SimpleBufferedServer_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [SimpleBufferedServer *]
 *
 */
SimpleBufferedServer * SimpleBufferedServer_alloc_len(int len) 
{
    SimpleBufferedServer * out; /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = SimpleBufferedServer_alloc()) == NULL) 
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->buffer = (BufferedString ** ) ckcalloc (len,sizeof(BufferedString *))) == NULL) {  
      warn("Warning, ckcalloc failed in SimpleBufferedServer_alloc_len");    
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_SimpleBufferedServer(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SimpleBufferedServer *]
 *
 * Return [UNKN ]  Undocumented return value [SimpleBufferedServer *]
 *
 */
SimpleBufferedServer * hard_link_SimpleBufferedServer(SimpleBufferedServer * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a SimpleBufferedServer object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  SimpleBufferedServer_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SimpleBufferedServer *]
 *
 */
SimpleBufferedServer * SimpleBufferedServer_alloc(void) 
{
    SimpleBufferedServer * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(SimpleBufferedServer *) ckalloc (sizeof(SimpleBufferedServer))) == NULL)    {  
      warn("SimpleBufferedServer_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->socket = 0; 
    out->buffer = NULL;  
    out->len = out->maxlen = 0;  
    out->current_read = 0;   


    return out;  
}    


/* Function:  free_SimpleBufferedServer(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SimpleBufferedServer *]
 *
 * Return [UNKN ]  Undocumented return value [SimpleBufferedServer *]
 *
 */
SimpleBufferedServer * free_SimpleBufferedServer(SimpleBufferedServer * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a SimpleBufferedServer obj. Should be trappable");  
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
    if( obj->buffer != NULL) {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->buffer[i] != NULL)  
          free_BufferedString(obj->buffer[i]);   
        }  
      ckfree(obj->buffer);   
      }  


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
