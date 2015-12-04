#ifndef DYNAMITEwisestreaminterfaceHEADERFILE
#define DYNAMITEwisestreaminterfaceHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "wisebase.h"


typedef struct Wise2_Stream_Read_Interface {
  char* (*read_buffer)(void *,char*,int);
  boolean (*is_end)(void*);
  void  (*close_and_free_handle)(void *);
  void * handle;
} Wise2ReadStreamInterface;


typedef struct Wise2_Stream_Write_Interface {
  void (*write_buffer)(void *,char*);
  void (*write_bufferf)(void *,char*,...);
  void (*close_and_free_handle)(void *);
  void * handle;
} Wise2WriteStreamInterface;


#define WISE2_READ_BUFFER(buffer,length,ri) ((*ri->read_buffer)(ri->handle,buffer,length))
#define WISE2_WRITE_STRING(string,wi) ((*wi->write_buffer)(wi->handle,string))




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  read_buffer_FILE_impl(handle,input_buffer,maxsize)
 *
 * Descrip:    Implementation function for normal files for reading
 *
 *
 * Arg:              handle [UNKN ] Undocumented argument [void *]
 * Arg:        input_buffer [UNKN ] Undocumented argument [char *]
 * Arg:             maxsize [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [char *]
 *
 */
char * Wise2_read_buffer_FILE_impl(void * handle,char * input_buffer,int maxsize);
#define read_buffer_FILE_impl Wise2_read_buffer_FILE_impl


/* Function:  is_end_FILE_impl(handle)
 *
 * Descrip:    Implementation function for normal files for end flag
 *
 *
 * Arg:        handle [UNKN ] Undocumented argument [void *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_is_end_FILE_impl(void * handle);
#define is_end_FILE_impl Wise2_is_end_FILE_impl


/* Function:  close_and_free_FILE_impl(handle)
 *
 * Descrip:    Implementation function for normal files for closing. Works for both
 *             reading and writing
 *
 *
 * Arg:        handle [UNKN ] Undocumented argument [void *]
 *
 */
void Wise2_close_and_free_FILE_impl(void * handle);
#define close_and_free_FILE_impl Wise2_close_and_free_FILE_impl


/* Function:  ReadStream_from_FILE(ifp)
 *
 * Descrip:    Makes a Wise2ReadStream interface from a normal C filehandle
 *
 *
 * Arg:        ifp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [Wise2ReadStreamInterface *]
 *
 */
Wise2ReadStreamInterface * Wise2_ReadStream_from_FILE(FILE * ifp);
#define ReadStream_from_FILE Wise2_ReadStream_from_FILE


/* Function:  ReadStream_openfile(filename)
 *
 * Descrip:    opens a file from filename and gives back a ReadStream,
 *             NULL if unable to return
 *
 *
 * Arg:        filename [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [Wise2ReadStreamInterface *]
 *
 */
Wise2ReadStreamInterface * Wise2_ReadStream_openfile(char * filename);
#define ReadStream_openfile Wise2_ReadStream_openfile


/* Function:  write_buffer_FILE_impl(handle,string)
 *
 * Descrip:    implementation for normal C files for string writing
 *
 *
 * Arg:        handle [UNKN ] Undocumented argument [void *]
 * Arg:        string [UNKN ] Undocumented argument [char *]
 *
 */
void Wise2_write_buffer_FILE_impl(void * handle,char * string);
#define write_buffer_FILE_impl Wise2_write_buffer_FILE_impl


/* Function:  write_bufferf_FILE_impl(handle,format,)
 *
 * Descrip:    implementation for normal C files for formatted writing
 *
 *
 * Arg:        handle [UNKN ] Undocumented argument [void *]
 * Arg:        format [UNKN ] Undocumented argument [char *]
 * Arg:               [UNKN ] Undocumented argument [.]
 *
 */
void Wise2_write_bufferf_FILE_impl(void * handle,char * format,...);
#define write_bufferf_FILE_impl Wise2_write_bufferf_FILE_impl


/* Function:  WriteStream_from_FILE(ofp)
 *
 * Descrip:    makes a WriteStream from a normal C FILE structure
 *
 *
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [Wise2WriteStreamInterface *]
 *
 */
Wise2WriteStreamInterface * Wise2_WriteStream_from_FILE(FILE * ofp);
#define WriteStream_from_FILE Wise2_WriteStream_from_FILE


/* Function:  WriteStream_openfile(filename)
 *
 * Descrip:    opens a file from filename and gives back a WriteStream,
 *             NULL if unable to return
 *
 *
 * Arg:        filename [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [Wise2WriteStreamInterface *]
 *
 */
Wise2WriteStreamInterface * Wise2_WriteStream_openfile(char * filename);
#define WriteStream_openfile Wise2_WriteStream_openfile


/* Function:  cat_ReadStream_into_WriteStream(rs,ws)
 *
 * Descrip:    helper function - cats one stream into another
 *
 *
 * Arg:        rs [UNKN ] Undocumented argument [Wise2ReadStreamInterface *]
 * Arg:        ws [UNKN ] Undocumented argument [Wise2WriteStreamInterface *]
 *
 */
void Wise2_cat_ReadStream_into_WriteStream(Wise2ReadStreamInterface * rs,Wise2WriteStreamInterface * ws);
#define cat_ReadStream_into_WriteStream Wise2_cat_ReadStream_into_WriteStream


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/

#ifdef _cplusplus
}
#endif

#endif
