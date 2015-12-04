#ifndef DYNAMITEassembly_stream_interfaceHEADERFILE
#define DYNAMITEassembly_stream_interfaceHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "assembly.h"


struct Wise2_AssemblySequenceStream {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    AssemblySequence * (*next_AssemblySequence)(void *); 
    void (*free_handle)(void *); 
    void * handle;   
    } ;  
/* AssemblySequenceStream defined */ 
#ifndef DYNAMITE_DEFINED_AssemblySequenceStream
typedef struct Wise2_AssemblySequenceStream Wise2_AssemblySequenceStream;
#define AssemblySequenceStream Wise2_AssemblySequenceStream
#define DYNAMITE_DEFINED_AssemblySequenceStream
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  free_AssemblySequenceStream(obj)
 *
 * Descrip:    provides specific deconstructor 
 *
 *
 * Arg:        obj [UNKN ] Undocumented argument [AssemblySequenceStream *]
 *
 * Return [UNKN ]  Undocumented return value [AssemblySequenceStream *]
 *
 */
AssemblySequenceStream * Wise2_free_AssemblySequenceStream(AssemblySequenceStream * obj);
#define free_AssemblySequenceStream Wise2_free_AssemblySequenceStream


/* Function:  hard_link_AssemblySequenceStream(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AssemblySequenceStream *]
 *
 * Return [UNKN ]  Undocumented return value [AssemblySequenceStream *]
 *
 */
AssemblySequenceStream * Wise2_hard_link_AssemblySequenceStream(AssemblySequenceStream * obj);
#define hard_link_AssemblySequenceStream Wise2_hard_link_AssemblySequenceStream


/* Function:  AssemblySequenceStream_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AssemblySequenceStream *]
 *
 */
AssemblySequenceStream * Wise2_AssemblySequenceStream_alloc(void);
#define AssemblySequenceStream_alloc Wise2_AssemblySequenceStream_alloc


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
