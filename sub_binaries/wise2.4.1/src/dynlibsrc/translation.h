#ifndef DYNAMITEtranslationHEADERFILE
#define DYNAMITEtranslationHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"

#ifndef DYNAMITE_DEFINED_Transcript
typedef struct Wise2_Transcript Wise2_Transcript;
#define Transcript Wise2_Transcript
#define DYNAMITE_DEFINED_Transcript
#endif

/* Object Translation
 *
 * Descrip: No Description
 *
 */
struct Wise2_Translation {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int start;   
    int end;     
    Transcript * parent;     
    Protein * protein;   
    } ;  
/* Translation defined */ 
#ifndef DYNAMITE_DEFINED_Translation
typedef struct Wise2_Translation Wise2_Translation;
#define Translation Wise2_Translation
#define DYNAMITE_DEFINED_Translation
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  copy_Translation(t)
 *
 * Descrip:    Makes a complete clean copy of the translation
 *
 *
 * Arg:        t [UNKN ] Undocumented argument [Translation *]
 *
 * Return [UNKN ]  Undocumented return value [Translation *]
 *
 */
Translation * Wise2_copy_Translation(Translation * t);
#define copy_Translation Wise2_copy_Translation


/* Function:  get_Protein_from_Translation(ts,ct)
 *
 * Descrip:    Gets the protein
 *
 *
 * Arg:        ts [UNKN ] translation [Translation *]
 * Arg:        ct [UNKN ] codon table to use [CodonTable *]
 *
 * Return [SOFT ]  Protein sequence [Protein *]
 *
 */
Protein * Wise2_get_Protein_from_Translation(Translation * ts,CodonTable * ct);
#define get_Protein_from_Translation Wise2_get_Protein_from_Translation


/* Function:  show_Translation(*ts,ofp)
 *
 * Descrip:    shows a translation in vaguely human form
 *
 *
 * Arg:        *ts [UNKN ] Undocumented argument [Translation]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_show_Translation(Translation *ts,FILE * ofp);
#define show_Translation Wise2_show_Translation


/* Function:  hard_link_Translation(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Translation *]
 *
 * Return [UNKN ]  Undocumented return value [Translation *]
 *
 */
Translation * Wise2_hard_link_Translation(Translation * obj);
#define hard_link_Translation Wise2_hard_link_Translation


/* Function:  Translation_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Translation *]
 *
 */
Translation * Wise2_Translation_alloc(void);
#define Translation_alloc Wise2_Translation_alloc


/* Function:  free_Translation(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Translation *]
 *
 * Return [UNKN ]  Undocumented return value [Translation *]
 *
 */
Translation * Wise2_free_Translation(Translation * obj);
#define free_Translation Wise2_free_Translation


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
boolean Wise2_replace_start_Translation(Translation * obj,int start);
#define replace_start_Translation Wise2_replace_start_Translation
boolean Wise2_replace_parent_Translation(Translation * obj,Transcript * parent);
#define replace_parent_Translation Wise2_replace_parent_Translation
int Wise2_access_start_Translation(Translation * obj);
#define access_start_Translation Wise2_access_start_Translation
Transcript * Wise2_access_parent_Translation(Translation * obj);
#define access_parent_Translation Wise2_access_parent_Translation
int Wise2_access_end_Translation(Translation * obj);
#define access_end_Translation Wise2_access_end_Translation
boolean Wise2_replace_protein_Translation(Translation * obj,Protein * protein);
#define replace_protein_Translation Wise2_replace_protein_Translation
boolean Wise2_replace_end_Translation(Translation * obj,int end);
#define replace_end_Translation Wise2_replace_end_Translation
Protein * Wise2_access_protein_Translation(Translation * obj);
#define access_protein_Translation Wise2_access_protein_Translation

#ifdef _cplusplus
}
#endif

#endif
