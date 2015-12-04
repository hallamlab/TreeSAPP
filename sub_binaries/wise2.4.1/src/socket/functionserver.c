#ifdef _cplusplus
extern "C" {
#endif
#include "functionserver.h"

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <stdio.h>
#include <stdlib.h>

#include <errno.h>	/* errno */
#include <string.h>	/* strerror() */

#include <signal.h>

#include <sys/time.h>		/* gettimeofday() */
#include <sys/resource.h>


# line 45 "functionserver.dy"
void main_loop_forking_FunctionServer(FunctionServer * fs,int verbose)
{
  int new_socket;
  char buf[1024];
  int i;

  int j;
  AnonymousObject * anon;
  AnonymousObjectList * aol;
  FILE * ifp;
  FILE * ofp;
  Wise2ReadStreamInterface  * rs;
  Wise2WriteStreamInterface * ws; 
  struct sigaction sigchange;

  ssize_t bytes_in;

  assert(fs != NULL);



  /* 
   * ignore child signals
   */

  sigchange.sa_handler = SIG_IGN;
  sigemptyset(&sigchange.sa_mask);

    
  listen(fs->socket,30);

  if( (sigaction(SIGCHLD,&sigchange,NULL) == -1)) {
    fatal("Unable to set signal handling behaviour");
  }

  while( 1 ) {
    if( verbose ) {
      info("Blocking on accept in main loop");
    } 

    new_socket = accept(fs->socket,0,0);

    if( fork() == 0 ) {
	int count = 0;
      /* forked process for client */

	struct timeval t0, t1 , t2 , t3;

	gettimeofday(&t0, NULL);

      if( verbose ) {
	info("Accepted connection");
      } 
    
      /* loop around possible multiple functions */
      while( 1 ) {
	++count;

	if ((bytes_in = read(new_socket,buf,sizeof(buf))) == -1) {
	  fatal(strerror(errno));
	}
	/* Terminate the recieved string */
	buf[bytes_in] = '\0';
	
	
	if( strncmp(buf,"FINISH",6) == 0) {
	  if( verbose ) {
	    info("Finished connection");
	  }
	  break;
	}
	
	for (i = 0; i < fs->len; i++) {
	  if (strncmp(buf, fs->fi[i]->transfer->name,
		      strlen(fs->fi[i]->transfer->name)) == 0 ) {
	    break;
	  }
	}
	
	if( i >= fs->len ) {
	  warn("Client asked for %s, unable to respond\n",buf);

	  strcpy(buf,"NOT FOUND");
	  buf[9] = '\0';
	  write(new_socket,buf,9);
	  close(new_socket);
	  /* exit out of child process */
          exit(0);
	}
	
	if( verbose ) 
	  info("Found %s for client, processing",fs->fi[i]->transfer->name);
	
	/* read the inputs, ask for data */
	
	strcpy(buf,"SEND\n\0");
	write(new_socket,buf,6);
	
	/* now inputs */
	
	aol = AnonymousObjectList_alloc_len(fs->fi[i]->transfer->len);
	ifp = fdopen(new_socket,"r");
	
	rs = ReadStream_from_FILE(ifp);
	
	for(j=0;j<fs->fi[i]->transfer->len;j++) {
	  anon = AnonymousObject_alloc();
	  anon->obj = (*fs->fi[i]->transfer->input[j]->stream_to_object)(rs);
	  anon->free_object = fs->fi[i]->transfer->input[j]->untyped_free_object;
	  
	  
	  add_AnonymousObjectList(aol,anon);
	}
	
	/* call function */
	
	gettimeofday(&t1, NULL);
	
	anon = (*fs->fi[i]->implementation)(fs->fi[i]->handle,aol);
	
	gettimeofday(&t2, NULL);

	/* send data back */
	
	if( verbose ) {
	  info("implementation successful, preparing to send result");
	}
	
	ws = WriteStream_from_socket(new_socket);
	(*fs->fi[i]->transfer->returned_type->object_to_stream)(anon->obj,ws);
	
	/* handle server side memory */
	
	free_AnonymousObjectList(aol);
	free_AnonymousObject(anon);
	
	
	strcpy(buf,"DONE\n");
	write(new_socket,buf,5);


	if( verbose ) {
	  info("sent return object");
	}
	
	
	gettimeofday(&t3, NULL);

      }
      close(new_socket);

      info("Exiting child having closed socket");

      /* client exits */
      exit(0);
    } else { /* fork != 0 */
      
      /* parent needs to close socket as well */
      close(new_socket);

      /* main server, back to accept */
    }
  } /* main loop */

  close(fs->socket);
}


# line 213 "functionserver.dy"
FunctionServer * new_FunctionServer(int port)
{
  struct sockaddr_in server;
  struct sockaddr * name;
  FunctionServer * out;
  int yes = 1;

  out = FunctionServer_alloc_std();

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

# line 239 "functionserver.dy"
FunctionImplementation * test_strcat_FunctionImplementation(void)
{
  FunctionImplementation * out;

  out = FunctionImplementation_alloc();
  out->transfer = test_stringcat_TransferedFunctionCall();
  out->implementation = test_strcat_implementation;

  return out;
}

# line 250 "functionserver.dy"
AnonymousObject * test_strcat_implementation(void * h,AnonymousObjectList * aol)
{
  char buffer[1024];
  char * one;
  char * two;
  AnonymousObject * ret;

  assert(aol != NULL);
  assert(aol->len == 2);

  one = (char*)aol->anon[0]->obj;
  two = (char*)aol->anon[1]->obj;

  ret = AnonymousObject_alloc();

  strcpy(buffer,one);
  strcat(buffer,two);

  ret->obj = stringalloc(buffer);

  return ret;
}

# line 253 "functionserver.c"
/* Function:  hard_link_FunctionImplementation(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [FunctionImplementation *]
 *
 * Return [UNKN ]  Undocumented return value [FunctionImplementation *]
 *
 */
FunctionImplementation * hard_link_FunctionImplementation(FunctionImplementation * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a FunctionImplementation object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  FunctionImplementation_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [FunctionImplementation *]
 *
 */
FunctionImplementation * FunctionImplementation_alloc(void) 
{
    FunctionImplementation * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(FunctionImplementation *) ckalloc (sizeof(FunctionImplementation))) == NULL)    {  
      warn("FunctionImplementation_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->transfer = NULL;    
    out->implementation = NULL;  


    return out;  
}    


/* Function:  free_FunctionImplementation(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [FunctionImplementation *]
 *
 * Return [UNKN ]  Undocumented return value [FunctionImplementation *]
 *
 */
FunctionImplementation * free_FunctionImplementation(FunctionImplementation * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a FunctionImplementation obj. Should be trappable");    
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
    /* obj->implementation is a function pointer */ 
    /* obj->handle is linked in */ 


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_FunctionServer(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_FunctionServer
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [FunctionImplementation **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_FunctionServer(FunctionImplementation ** list,int i,int j)  
{
    FunctionImplementation * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_FunctionServer(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_FunctionServer which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [FunctionImplementation **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_FunctionServer(FunctionImplementation ** list,int left,int right,int (*comp)(FunctionImplementation * ,FunctionImplementation * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_FunctionServer(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_FunctionServer (list,++last,i); 
      }  
    swap_FunctionServer (list,left,last);    
    qsort_FunctionServer(list,left,last-1,comp); 
    qsort_FunctionServer(list,last+1,right,comp);    
}    


/* Function:  sort_FunctionServer(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_FunctionServer
 *
 *
 * Arg:         obj [UNKN ] Object containing list [FunctionServer *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_FunctionServer(FunctionServer * obj,int (*comp)(FunctionImplementation *, FunctionImplementation *)) 
{
    qsort_FunctionServer(obj->fi,0,obj->len-1,comp); 
    return;  
}    


/* Function:  expand_FunctionServer(obj,len)
 *
 * Descrip:    Really an internal function for add_FunctionServer
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [FunctionServer *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_FunctionServer(FunctionServer * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_FunctionServer called with no need"); 
      return TRUE;   
      }  


    if( (obj->fi = (FunctionImplementation ** ) ckrealloc (obj->fi,sizeof(FunctionImplementation *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_FunctionServer, returning FALSE");   
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_FunctionServer(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [FunctionServer *]
 * Arg:        add [OWNER] Object to add to the list [FunctionImplementation *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_FunctionServer(FunctionServer * obj,FunctionImplementation * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_FunctionServer(obj,obj->len + FunctionServerLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->fi[obj->len++]=add; 
    return TRUE; 
}    


/* Function:  flush_FunctionServer(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [FunctionServer *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_FunctionServer(FunctionServer * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->fi[i] != NULL)    {  
        free_FunctionImplementation(obj->fi[i]); 
        obj->fi[i] = NULL;   
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  FunctionServer_alloc_std(void)
 *
 * Descrip:    Equivalent to FunctionServer_alloc_len(FunctionServerLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [FunctionServer *]
 *
 */
FunctionServer * FunctionServer_alloc_std(void) 
{
    return FunctionServer_alloc_len(FunctionServerLISTLENGTH);   
}    


/* Function:  FunctionServer_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [FunctionServer *]
 *
 */
FunctionServer * FunctionServer_alloc_len(int len) 
{
    FunctionServer * out;   /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = FunctionServer_alloc()) == NULL)   
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->fi = (FunctionImplementation ** ) ckcalloc (len,sizeof(FunctionImplementation *))) == NULL) {  
      warn("Warning, ckcalloc failed in FunctionServer_alloc_len");  
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_FunctionServer(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [FunctionServer *]
 *
 * Return [UNKN ]  Undocumented return value [FunctionServer *]
 *
 */
FunctionServer * hard_link_FunctionServer(FunctionServer * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a FunctionServer object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  FunctionServer_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [FunctionServer *]
 *
 */
FunctionServer * FunctionServer_alloc(void) 
{
    FunctionServer * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(FunctionServer *) ckalloc (sizeof(FunctionServer))) == NULL)    {  
      warn("FunctionServer_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->fi = NULL;  
    out->len = out->maxlen = 0;  
    out->socket = 0; 


    return out;  
}    


/* Function:  free_FunctionServer(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [FunctionServer *]
 *
 * Return [UNKN ]  Undocumented return value [FunctionServer *]
 *
 */
FunctionServer * free_FunctionServer(FunctionServer * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a FunctionServer obj. Should be trappable");    
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
    if( obj->fi != NULL) {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->fi[i] != NULL)  
          free_FunctionImplementation(obj->fi[i]);   
        }  
      ckfree(obj->fi);   
      }  


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
