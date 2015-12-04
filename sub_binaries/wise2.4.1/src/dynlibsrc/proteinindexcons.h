#ifndef DYNAMITEproteinindexconsHEADERFILE
#define DYNAMITEproteinindexconsHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif

#include "arrayseqlookup.h"
#include "subseqhash.h"
#include "proteinstreamedindex.h"
#include "shadowseqindex.h"

typedef enum ProteinIndexConstructorType {
  ProteinIndexConstructor_Array = 78,
  ProteinIndexConstructor_Hash,
  ProteinIndexConstructor_Stream,
  ProteinIndexConstructor_Shadow
} ProteinIndexConstructureType;


struct Wise2_ProteinIndexConstructor {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int type;    
    int waypost;     
    int shadowlength;    
    int has_maxlen;  
    int max_seqlen;  
    int shadow_error;    
    } ;  
/* ProteinIndexConstructor defined */ 
#ifndef DYNAMITE_DEFINED_ProteinIndexConstructor
typedef struct Wise2_ProteinIndexConstructor Wise2_ProteinIndexConstructor;
#define ProteinIndexConstructor Wise2_ProteinIndexConstructor
#define DYNAMITE_DEFINED_ProteinIndexConstructor
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  new_SeqLookupInterface_from_ProteinIndexConstructor(pic)
 *
 * Descrip:    Makes a new SeqLookupInterface from ProteinIndexConstructor
 *
 *
 * Arg:        pic [UNKN ] Undocumented argument [ProteinIndexConstructor *]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupInterface *]
 *
 */
SeqLookupInterface * Wise2_new_SeqLookupInterface_from_ProteinIndexConstructor(ProteinIndexConstructor * pic);
#define new_SeqLookupInterface_from_ProteinIndexConstructor Wise2_new_SeqLookupInterface_from_ProteinIndexConstructor


/* Function:  show_help_ProteinIndexConstructor(ofp)
 *
 * Descrip:    provides help for protein index constructor
 *
 *
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_show_help_ProteinIndexConstructor(FILE * ofp);
#define show_help_ProteinIndexConstructor Wise2_show_help_ProteinIndexConstructor


/* Function:  new_ProteinIndexConstructor_from_argv(argc,argv)
 *
 * Descrip:    Provides a ProteinIndexConstructor argument from argv
 *
 *
 * Arg:        argc [UNKN ] Undocumented argument [int *]
 * Arg:        argv [UNKN ] Undocumented argument [char **]
 *
 * Return [UNKN ]  Undocumented return value [ProteinIndexConstructor *]
 *
 */
ProteinIndexConstructor * Wise2_new_ProteinIndexConstructor_from_argv(int * argc,char ** argv);
#define new_ProteinIndexConstructor_from_argv Wise2_new_ProteinIndexConstructor_from_argv


/* Function:  hard_link_ProteinIndexConstructor(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ProteinIndexConstructor *]
 *
 * Return [UNKN ]  Undocumented return value [ProteinIndexConstructor *]
 *
 */
ProteinIndexConstructor * Wise2_hard_link_ProteinIndexConstructor(ProteinIndexConstructor * obj);
#define hard_link_ProteinIndexConstructor Wise2_hard_link_ProteinIndexConstructor


/* Function:  ProteinIndexConstructor_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ProteinIndexConstructor *]
 *
 */
ProteinIndexConstructor * Wise2_ProteinIndexConstructor_alloc(void);
#define ProteinIndexConstructor_alloc Wise2_ProteinIndexConstructor_alloc


/* Function:  free_ProteinIndexConstructor(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ProteinIndexConstructor *]
 *
 * Return [UNKN ]  Undocumented return value [ProteinIndexConstructor *]
 *
 */
ProteinIndexConstructor * Wise2_free_ProteinIndexConstructor(ProteinIndexConstructor * obj);
#define free_ProteinIndexConstructor Wise2_free_ProteinIndexConstructor


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
