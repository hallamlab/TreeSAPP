#ifndef DYNAMITEsimplesocketclientHEADERFILE
#define DYNAMITEsimplesocketclientHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "wisebase.h"

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>


typedef struct Wise2_SimpleRead_Socket {
  int socket;
  FILE * socket_ifp;
  boolean has_ended;
  char buffer[1024];
  char * curr_read;
  int curr_len;
} Wise2_SimpleRead_Socket;



    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
char * Wise2_read_buffer_simple_socket_impl(void * handle,char * buffer,int maxsize);
#define read_buffer_simple_socket_impl Wise2_read_buffer_simple_socket_impl
boolean Wise2_is_end_simple_socket_impl(void * handle);
#define is_end_simple_socket_impl Wise2_is_end_simple_socket_impl
void Wise2_close_and_free_simple_socket_impl(void * handle);
#define close_and_free_simple_socket_impl Wise2_close_and_free_simple_socket_impl
Wise2ReadStreamInterface * Wise2_open_simple_socket_for_reading(char * hostname,int port_number);
#define open_simple_socket_for_reading Wise2_open_simple_socket_for_reading


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/

#ifdef _cplusplus
}
#endif

#endif
