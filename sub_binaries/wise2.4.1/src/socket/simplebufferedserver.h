#ifndef DYNAMITEsimplebufferedserverHEADERFILE
#define DYNAMITEsimplebufferedserverHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "wisebase.h"
#include <errno.h>


typedef char BufferedString;

#define SimpleBufferedServerLISTLENGTH 1024

struct Wise2_SimpleBufferedServer {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int socket;  
    BufferedString ** buffer;    
    int len;/* len for above buffer  */ 
    int maxlen; /* maxlen for above buffer */ 
    int current_read;    
    } ;  
/* SimpleBufferedServer defined */ 
#ifndef DYNAMITE_DEFINED_SimpleBufferedServer
typedef struct Wise2_SimpleBufferedServer Wise2_SimpleBufferedServer;
#define SimpleBufferedServer Wise2_SimpleBufferedServer
#define DYNAMITE_DEFINED_SimpleBufferedServer
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



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
boolean Wise2_add_SimpleBufferedServer(SimpleBufferedServer * obj,BufferedString * add);
#define add_SimpleBufferedServer Wise2_add_SimpleBufferedServer


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
int Wise2_flush_SimpleBufferedServer(SimpleBufferedServer * obj);
#define flush_SimpleBufferedServer Wise2_flush_SimpleBufferedServer


/* Function:  SimpleBufferedServer_alloc_std(void)
 *
 * Descrip:    Equivalent to SimpleBufferedServer_alloc_len(SimpleBufferedServerLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SimpleBufferedServer *]
 *
 */
SimpleBufferedServer * Wise2_SimpleBufferedServer_alloc_std(void);
#define SimpleBufferedServer_alloc_std Wise2_SimpleBufferedServer_alloc_std


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
SimpleBufferedServer * Wise2_SimpleBufferedServer_alloc_len(int len);
#define SimpleBufferedServer_alloc_len Wise2_SimpleBufferedServer_alloc_len


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
SimpleBufferedServer * Wise2_hard_link_SimpleBufferedServer(SimpleBufferedServer * obj);
#define hard_link_SimpleBufferedServer Wise2_hard_link_SimpleBufferedServer


/* Function:  SimpleBufferedServer_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SimpleBufferedServer *]
 *
 */
SimpleBufferedServer * Wise2_SimpleBufferedServer_alloc(void);
#define SimpleBufferedServer_alloc Wise2_SimpleBufferedServer_alloc


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
SimpleBufferedServer * Wise2_free_SimpleBufferedServer(SimpleBufferedServer * obj);
#define free_SimpleBufferedServer Wise2_free_SimpleBufferedServer


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
void Wise2_write_buffer_simple_buffered_server_impl(void * handle,char * string);
#define write_buffer_simple_buffered_server_impl Wise2_write_buffer_simple_buffered_server_impl
void Wise2_write_bufferf_simple_buffered_server_impl(void * handle,char * format,...);
#define write_bufferf_simple_buffered_server_impl Wise2_write_bufferf_simple_buffered_server_impl
void Wise2_close_and_free_simple_buffered_server(void * handle);
#define close_and_free_simple_buffered_server Wise2_close_and_free_simple_buffered_server
void Wise2_add_string_SimpleBufferedServer(SimpleBufferedServer * sbs,char * unallocated_string);
#define add_string_SimpleBufferedServer Wise2_add_string_SimpleBufferedServer
Wise2WriteStreamInterface * Wise2_WriteStream_from_SimpleBufferedServer(SimpleBufferedServer * sbs);
#define WriteStream_from_SimpleBufferedServer Wise2_WriteStream_from_SimpleBufferedServer
SimpleBufferedServer * Wise2_new_SimpleBufferedServer(int port);
#define new_SimpleBufferedServer Wise2_new_SimpleBufferedServer
char * Wise2_free_BufferedString(char * str);
#define free_BufferedString Wise2_free_BufferedString
void Wise2_listen_and_respond_SimpleBufferedServer(SimpleBufferedServer * sbs);
#define listen_and_respond_SimpleBufferedServer Wise2_listen_and_respond_SimpleBufferedServer


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_SimpleBufferedServer(BufferedString ** list,int i,int j) ;
#define swap_SimpleBufferedServer Wise2_swap_SimpleBufferedServer
void Wise2_qsort_SimpleBufferedServer(BufferedString ** list,int left,int right,int (*comp)(BufferedString * ,BufferedString * ));
#define qsort_SimpleBufferedServer Wise2_qsort_SimpleBufferedServer
void Wise2_sort_SimpleBufferedServer(SimpleBufferedServer * obj,int (*comp)(BufferedString *, BufferedString *));
#define sort_SimpleBufferedServer Wise2_sort_SimpleBufferedServer
boolean Wise2_expand_SimpleBufferedServer(SimpleBufferedServer * obj,int len);
#define expand_SimpleBufferedServer Wise2_expand_SimpleBufferedServer

#ifdef _cplusplus
}
#endif

#endif
