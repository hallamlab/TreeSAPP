#ifndef DYNAMITEassembly_stream_consHEADERFILE
#define DYNAMITEassembly_stream_consHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "assembly_stream_interface.h"
#include "assembly_sanger_project.h"
#include "assembly_stream_fasta.h"

typedef enum assembly_stream_type {
  AssemblyStreamTypeFasta = 791,
  AssemblyStreamTypeSanger,
} AssemblyStreamType;

struct Wise2_AssemblyStreamConstructor {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int type;    
    char * file_name;    
    char * dir_name;     
    char * extension;    
    } ;  
/* AssemblyStreamConstructor defined */ 
#ifndef DYNAMITE_DEFINED_AssemblyStreamConstructor
typedef struct Wise2_AssemblyStreamConstructor Wise2_AssemblyStreamConstructor;
#define AssemblyStreamConstructor Wise2_AssemblyStreamConstructor
#define DYNAMITE_DEFINED_AssemblyStreamConstructor
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_AssemblyStreamConstructor(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AssemblyStreamConstructor *]
 *
 * Return [UNKN ]  Undocumented return value [AssemblyStreamConstructor *]
 *
 */
AssemblyStreamConstructor * Wise2_hard_link_AssemblyStreamConstructor(AssemblyStreamConstructor * obj);
#define hard_link_AssemblyStreamConstructor Wise2_hard_link_AssemblyStreamConstructor


/* Function:  AssemblyStreamConstructor_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AssemblyStreamConstructor *]
 *
 */
AssemblyStreamConstructor * Wise2_AssemblyStreamConstructor_alloc(void);
#define AssemblyStreamConstructor_alloc Wise2_AssemblyStreamConstructor_alloc


/* Function:  free_AssemblyStreamConstructor(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AssemblyStreamConstructor *]
 *
 * Return [UNKN ]  Undocumented return value [AssemblyStreamConstructor *]
 *
 */
AssemblyStreamConstructor * Wise2_free_AssemblyStreamConstructor(AssemblyStreamConstructor * obj);
#define free_AssemblyStreamConstructor Wise2_free_AssemblyStreamConstructor


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
void Wise2_show_help_AssemblyStreamConstructor(FILE * ofp);
#define show_help_AssemblyStreamConstructor Wise2_show_help_AssemblyStreamConstructor
AssemblySequenceStream * Wise2_new_AssemblySequenceStream_from_AssemblyStreamConstructor(AssemblyStreamConstructor * asc);
#define new_AssemblySequenceStream_from_AssemblyStreamConstructor Wise2_new_AssemblySequenceStream_from_AssemblyStreamConstructor
AssemblyStreamConstructor * Wise2_new_AssemblyStreamConstructor_from_argv(int * argc,char ** argv);
#define new_AssemblyStreamConstructor_from_argv Wise2_new_AssemblyStreamConstructor_from_argv


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/

#ifdef _cplusplus
}
#endif

#endif
