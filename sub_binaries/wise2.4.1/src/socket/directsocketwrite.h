#ifndef DYNAMITEdirectsocketwriteHEADERFILE
#define DYNAMITEdirectsocketwriteHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "wisebase.h"


struct Wise2_DirectSocketWrite {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int socket;  
    } ;  
/* DirectSocketWrite defined */ 
#ifndef DYNAMITE_DEFINED_DirectSocketWrite
typedef struct Wise2_DirectSocketWrite Wise2_DirectSocketWrite;
#define DirectSocketWrite Wise2_DirectSocketWrite
#define DYNAMITE_DEFINED_DirectSocketWrite
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



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
DirectSocketWrite * Wise2_hard_link_DirectSocketWrite(DirectSocketWrite * obj);
#define hard_link_DirectSocketWrite Wise2_hard_link_DirectSocketWrite


/* Function:  DirectSocketWrite_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DirectSocketWrite *]
 *
 */
DirectSocketWrite * Wise2_DirectSocketWrite_alloc(void);
#define DirectSocketWrite_alloc Wise2_DirectSocketWrite_alloc


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
DirectSocketWrite * Wise2_free_DirectSocketWrite(DirectSocketWrite * obj);
#define free_DirectSocketWrite Wise2_free_DirectSocketWrite


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
void Wise2_write_buffer_direct_socket_server_impl(void * handle,char * string);
#define write_buffer_direct_socket_server_impl Wise2_write_buffer_direct_socket_server_impl
void Wise2_write_bufferf_direct_socket_server_impl(void * handle,char * format,...);
#define write_bufferf_direct_socket_server_impl Wise2_write_bufferf_direct_socket_server_impl
void Wise2_close_and_free_direct_socket_server(void * handle);
#define close_and_free_direct_socket_server Wise2_close_and_free_direct_socket_server
Wise2WriteStreamInterface * Wise2_WriteStream_from_socket(int socket);
#define WriteStream_from_socket Wise2_WriteStream_from_socket
boolean Wise2_write_DirectSocketWrite(DirectSocketWrite * dsw,char * string);
#define write_DirectSocketWrite Wise2_write_DirectSocketWrite


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/

#ifdef _cplusplus
}
#endif

#endif
