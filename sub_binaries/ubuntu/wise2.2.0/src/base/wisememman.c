#ifdef _cplusplus
extern "C" {
#endif
#include "wisememman.h"


#ifdef WISE_MEMORY_WATCH
 struct wise_memory_watcher {
   void * position;
   int amount;
   int has_freed;
   char * file;
   int lineno;
 };

 static struct wise_memory_watcher mempool[WISE_MEMORY_WATCH_MAX];
 static int nextpool = 0;

/* Function:  display_allocated_memory(tag,ofp)
 *
 * Descrip:    Displays to filehandle all non freed memory
 *
 *
 * Arg:        tag [UNKN ] Undocumented argument [char *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 31 "wisememman.dy"
void display_allocated_memory(char * tag,FILE * ofp)
{
  int i;
  for(i=0;i < nextpool && i<WISE_MEMORY_WATCH_MAX;i++)
    if( mempool[i].has_freed == FALSE ) {
      if( mempool[i].file != NULL ) 
	fprintf(ofp,"%s [%s:%d] %d of %d bytes not free\n",tag,mempool[i].file,mempool[i].lineno,mempool[i].position,mempool[i].amount);
      else
	fprintf(ofp,"%s %d of %d bytes not free\n",tag,mempool[i].position,mempool[i].amount);
      fprintf(ofp,"Bytes: %c%c\n",*((char *)mempool[i].position),*((char *)mempool[i].position +1));
    }

}


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
# line 52 "wisememman.dy"
void * allocate_watched_memory_file(char * file,int lineno,int no_bytes)
{
  if( nextpool >= WISE_MEMORY_WATCH_MAX ) {
    warn("Memory allocation over limits for watching");
    return malloc(no_bytes);
  }

  mempool[nextpool].position = malloc(no_bytes);
  mempool[nextpool].amount = no_bytes;
  mempool[nextpool].has_freed = FALSE;
  mempool[nextpool].file = file;
  mempool[nextpool].lineno = lineno;
  

  nextpool++;

  return mempool[nextpool-1].position;
}


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
# line 77 "wisememman.dy"
void * allocate_watched_memory(int no_bytes)
{
  if( nextpool >= WISE_MEMORY_WATCH_MAX ) {
    warn("Memory allocation over limits for watching");
    return malloc(no_bytes);
  }

  mempool[nextpool].position = malloc(no_bytes);
  mempool[nextpool].amount = no_bytes;
  mempool[nextpool].has_freed = FALSE;
  mempool[nextpool].file = NULL;
  mempool[nextpool].lineno = -1;

  nextpool++;

  return mempool[nextpool-1].position;
}

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
# line 100 "wisememman.dy"
void * free_watched_memory(void * mem)
{
  int i;
  for(i=0;i < nextpool && i<WISE_MEMORY_WATCH_MAX;i++)
    if( mempool[i].position == mem ) {
      break;
    }

  if( i == nextpool || i == WISE_MEMORY_WATCH_MAX ) {
    warn("Problem! memory position %d not watched",(int)mem);
    return NULL;
  }


  free(mempool[i].position); 
  mempool[i].has_freed = TRUE;

  return NULL;
}

#endif /* WISE_MEMORY_WATCH */

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
# line 126 "wisememman.dy"
void * ckalloc(size_t bytes)
{
  register void *ret;
  extern void *calloc (size_t nelem, size_t elsize);

#ifdef WISE_MEMORY_WATCH
  /* call into the watched memory pool */
  ret = allocate_watched_memory(bytes);
  if( ret == NULL ) {
    warn("Out of memory (watched) on %d bytes\n",bytes);
    return NULL;
  } else 
    return ret;
  
#endif

  if( (ret = calloc(bytes, sizeof(char))) == NULL) {
    warn("Out of memory, on asking for %d bytes\n",bytes);
    return NULL; /*** for the moment, could fail here ***/
  } else
    return ret;	
}

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
# line 152 "wisememman.dy"
void * ckcalloc(int len,size_t bytes)
{
  return ckalloc(len*bytes);
}

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
# line 160 "wisememman.dy"
void * ckrealloc(void *ptr, size_t bytes)
{
  register void *ret;
  extern void *realloc (void *ptr, size_t size);

  if (ptr == NULL) {	
    warn("Bad call to ckrealloc, NULL pointer\n");
    return NULL;
  }
  else if( (ret = realloc(ptr, bytes)) == NULL) {
    warn("Out of memory, trying to realloc %d bytes\n",bytes);
    return NULL;
  }
  else
    return ret;	
}

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
# line 180 "wisememman.dy"
void * ckfree(void *ptr)
{

#ifdef WISE_MEMORY_WATCH
  free_watched_memory(ptr);
  return NULL;
#endif

  if (ptr == NULL)
    warn("Bad call to ckfree - NULL pointer\n");
  else {
    free(ptr);
    ptr = NULL;
  }
  return ptr;
}





# line 242 "wisememman.c"

#ifdef _cplusplus
}
#endif
