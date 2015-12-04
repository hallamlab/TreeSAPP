#ifndef DYNAMITEassembly_sanger_projectHEADERFILE
#define DYNAMITEassembly_sanger_projectHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "assembly_stream_interface.h"

#include <sys/types.h>
#include <dirent.h>

struct Wise2_SangerProjectDirectory {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    DIR  * dir;  
    char * extension;    
    char * directory;    
    } ;  
/* SangerProjectDirectory defined */ 
#ifndef DYNAMITE_DEFINED_SangerProjectDirectory
typedef struct Wise2_SangerProjectDirectory Wise2_SangerProjectDirectory;
#define SangerProjectDirectory Wise2_SangerProjectDirectory
#define DYNAMITE_DEFINED_SangerProjectDirectory
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  free_SangerProjectDirectory(obj)
 *
 * Descrip:    deconstructor
 *
 *
 * Arg:        obj [UNKN ] Undocumented argument [SangerProjectDirectory *]
 *
 * Return [UNKN ]  Undocumented return value [SangerProjectDirectory *]
 *
 */
SangerProjectDirectory * Wise2_free_SangerProjectDirectory(SangerProjectDirectory * obj);
#define free_SangerProjectDirectory Wise2_free_SangerProjectDirectory


/* Function:  hard_link_SangerProjectDirectory(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SangerProjectDirectory *]
 *
 * Return [UNKN ]  Undocumented return value [SangerProjectDirectory *]
 *
 */
SangerProjectDirectory * Wise2_hard_link_SangerProjectDirectory(SangerProjectDirectory * obj);
#define hard_link_SangerProjectDirectory Wise2_hard_link_SangerProjectDirectory


/* Function:  SangerProjectDirectory_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SangerProjectDirectory *]
 *
 */
SangerProjectDirectory * Wise2_SangerProjectDirectory_alloc(void);
#define SangerProjectDirectory_alloc Wise2_SangerProjectDirectory_alloc


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
AssemblySequence * Wise2_next_AssemblySequence_sanger_impl(void * h);
#define next_AssemblySequence_sanger_impl Wise2_next_AssemblySequence_sanger_impl
AssemblySequence * Wise2_read_clipped_sanger_AssemblySequence(FILE * ifp);
#define read_clipped_sanger_AssemblySequence Wise2_read_clipped_sanger_AssemblySequence
void Wise2_free_handle_sanger_impl(void * h);
#define free_handle_sanger_impl Wise2_free_handle_sanger_impl
AssemblySequenceStream * Wise2_new_sanger_project_AssemblySequenceStream(char * dir_name,char * extension);
#define new_sanger_project_AssemblySequenceStream Wise2_new_sanger_project_AssemblySequenceStream
SangerProjectDirectory * Wise2_new_SangerProjectDirectory(char * dir_name,char * extension);
#define new_SangerProjectDirectory Wise2_new_SangerProjectDirectory


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/

#ifdef _cplusplus
}
#endif

#endif
