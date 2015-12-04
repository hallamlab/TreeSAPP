#ifndef DYNAMITEwisememmanHEADERFILE
#define DYNAMITEwisememmanHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif

#define CKALLOC_GUARD
#include "wisebase.h"
#undef  CKALLOC_GUARD

#define WISE_MEMORY_WATCH_MAX 4096



    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  display_allocated_memory(tag,ofp)
 *
 * Descrip:    Displays to filehandle all non freed memory
 *
 *
 * Arg:        tag [UNKN ] Undocumented argument [char *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_display_allocated_memory(char * tag,FILE * ofp);
#define display_allocated_memory Wise2_display_allocated_memory


/* Function:  allocate_watched_memory_file(file,lineno,no_bytes)
 *
 * Descrip:    Only if compiled with WISE_MEMORY_WATCH
 *
 *             Allocates memory and watches the usage of it,
 *             knowing the file and line number
 *
 *
 * Arg:            file [UNKN ] Undocumented argument [char *]
 * Arg:          lineno [UNKN ] Undocumented argument [int]
 * Arg:        no_bytes [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [void *]
 *
 */
void * Wise2_allocate_watched_memory_file(char * file,int lineno,int no_bytes);
#define allocate_watched_memory_file Wise2_allocate_watched_memory_file


/* Function:  allocate_watched_memory(no_bytes)
 *
 * Descrip:    Only if compiled with WISE_MEMORY_WATCH
 *
 *             Allocates memory and watches the usage of it
 *
 *
 * Arg:        no_bytes [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [void *]
 *
 */
void * Wise2_allocate_watched_memory(int no_bytes);
#define allocate_watched_memory Wise2_allocate_watched_memory


/* Function:  free_watched_memory(mem)
 *
 * Descrip:    Only if compiled with WISE_MEMORY_WATCH
 *
 *             frees memory that has been watched
 *
 *
 * Arg:        mem [UNKN ] Undocumented argument [void *]
 *
 * Return [UNKN ]  Undocumented return value [void *]
 *
 */
void * Wise2_free_watched_memory(void * mem);
#define free_watched_memory Wise2_free_watched_memory


/* Function:  ckalloc(bytes)
 *
 * Descrip:    Tries to alloc bytes of memory. Posts
 *             to warn if it fails
 *
 *
 * Arg:        bytes [UNKN ] Undocumented argument [size_t]
 *
 * Return [UNKN ]  Undocumented return value [void *]
 *
 */
void * Wise2_ckalloc(size_t bytes);
#define ckalloc Wise2_ckalloc


/* Function:  ckcalloc(len,bytes)
 *
 * Descrip:    calloc equivalent
 *
 *
 * Arg:          len [UNKN ] Undocumented argument [int]
 * Arg:        bytes [UNKN ] Undocumented argument [size_t]
 *
 * Return [UNKN ]  Undocumented return value [void *]
 *
 */
void * Wise2_ckcalloc(int len,size_t bytes);
#define ckcalloc Wise2_ckcalloc


/* Function:  ckrealloc(*ptr,bytes)
 *
 * Descrip:    realloc equivalent
 *
 *
 * Arg:         *ptr [UNKN ] Undocumented argument [void]
 * Arg:        bytes [UNKN ] Undocumented argument [size_t]
 *
 * Return [UNKN ]  Undocumented return value [void *]
 *
 */
void * Wise2_ckrealloc(void *ptr, size_t bytes);
#define ckrealloc Wise2_ckrealloc


/* Function:  ckfree(*ptr)
 *
 * Descrip:    free equivalent
 *
 *
 * Arg:        *ptr [UNKN ] Undocumented argument [void]
 *
 * Return [UNKN ]  Undocumented return value [void *]
 *
 */
void * Wise2_ckfree(void *ptr);
#define ckfree Wise2_ckfree


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
