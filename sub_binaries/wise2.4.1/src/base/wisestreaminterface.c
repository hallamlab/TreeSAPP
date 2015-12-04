#ifdef _cplusplus
extern "C" {
#endif
#include "wisestreaminterface.h"



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
# line 40 "wisestreaminterface.dy"
char * read_buffer_FILE_impl(void * handle,char * input_buffer,int maxsize)
{
  char * ret;
  FILE * ifp = (FILE *)handle;

  ret = fgets(input_buffer,maxsize,ifp);

  return ret;
}

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
# line 53 "wisestreaminterface.dy"
boolean is_end_FILE_impl(void * handle)
{
  FILE * ifp = (FILE *)handle;

  if( feof(ifp) ) {
    return TRUE;
  } else {
    return FALSE;
  }
}

/* Function:  close_and_free_FILE_impl(handle)
 *
 * Descrip:    Implementation function for normal files for closing. Works for both
 *             reading and writing
 *
 *
 * Arg:        handle [UNKN ] Undocumented argument [void *]
 *
 */
# line 68 "wisestreaminterface.dy"
void close_and_free_FILE_impl(void * handle)
{
  FILE * ifp = (FILE *)handle;

  fclose(ifp);
}

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
# line 78 "wisestreaminterface.dy"
Wise2ReadStreamInterface * ReadStream_from_FILE(FILE * ifp)
{
  Wise2ReadStreamInterface * out;

  out = malloc(sizeof(Wise2ReadStreamInterface));
  out->read_buffer = read_buffer_FILE_impl;
  out->is_end      = is_end_FILE_impl;
  out->close_and_free_handle = close_and_free_FILE_impl;
  out->handle = (void*) ifp;

  return out;
}

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
# line 95 "wisestreaminterface.dy"
Wise2ReadStreamInterface * ReadStream_openfile(char * filename)
{
  FILE * ifp;
  assert(filename != NULL);

  ifp = openfile(filename,"r");
  if( ifp == NULL ) {
    warn("Unable to open file %s",filename);
    return NULL;
  }
  
  return ReadStream_from_FILE(ifp);
}


/* Function:  write_buffer_FILE_impl(handle,string)
 *
 * Descrip:    implementation for normal C files for string writing
 *
 *
 * Arg:        handle [UNKN ] Undocumented argument [void *]
 * Arg:        string [UNKN ] Undocumented argument [char *]
 *
 */
# line 113 "wisestreaminterface.dy"
void write_buffer_FILE_impl(void * handle,char * string)
{
  FILE * ifp = (FILE *)handle;

  fputs(string,ifp);
  
  return;
}

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
# line 125 "wisestreaminterface.dy"
void write_bufferf_FILE_impl(void * handle,char * format,...)
{
  FILE * ifp = (FILE *)handle;
  char buffer[1024];
  va_list ap;


  va_start(ap,format);
  vsprintf(buffer,format,ap);	

  fputs(buffer,ifp);

  return;
}

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
# line 143 "wisestreaminterface.dy"
Wise2WriteStreamInterface * WriteStream_from_FILE(FILE * ofp)
{
  Wise2WriteStreamInterface * out;

  out = malloc(sizeof(Wise2WriteStreamInterface));
  out->write_buffer = write_buffer_FILE_impl;
  out->write_bufferf = write_bufferf_FILE_impl;
  out->close_and_free_handle = close_and_free_FILE_impl;
  out->handle = (void*) ofp;

  return out;
}

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
# line 160 "wisestreaminterface.dy"
Wise2WriteStreamInterface * WriteStream_openfile(char * filename)
{
  FILE * ofp;
  assert(filename != NULL);

  ofp = openfile(filename,"w");
  if( ofp == NULL ) {
    warn("Unable to open file %s",filename);
    return NULL;
  }
  
  return WriteStream_from_FILE(ofp);
}

/* Function:  cat_ReadStream_into_WriteStream(rs,ws)
 *
 * Descrip:    helper function - cats one stream into another
 *
 *
 * Arg:        rs [UNKN ] Undocumented argument [Wise2ReadStreamInterface *]
 * Arg:        ws [UNKN ] Undocumented argument [Wise2WriteStreamInterface *]
 *
 */
# line 177 "wisestreaminterface.dy"
void cat_ReadStream_into_WriteStream(Wise2ReadStreamInterface * rs,Wise2WriteStreamInterface * ws)
{
  char str[1024];

  while( (*rs->is_end)(rs->handle) != TRUE ) {
    if( (*rs->read_buffer)(rs->handle,str,1024) == NULL ) {
	return;
    }	
    (*ws->write_buffer)(ws->handle,str);
  }

  
}


# line 220 "wisestreaminterface.c"

#ifdef _cplusplus
}
#endif
