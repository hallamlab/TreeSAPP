#ifndef DYNAMITEtranscriptHEADERFILE
#define DYNAMITEtranscriptHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"

#define TranscriptLISTLENGTH 32
#define ExonLISTLENGTH 128

struct Wise2_SupportingFeature {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int start;  /*  in exon coordinates */ 
    int end;    /*  in exon coordinates */ 
    int hstart;  
    int hend;    
    int hstrand;     
    char* hid;   
    } ;  
/* SupportingFeature defined */ 
#ifndef DYNAMITE_DEFINED_SupportingFeature
typedef struct Wise2_SupportingFeature Wise2_SupportingFeature;
#define SupportingFeature Wise2_SupportingFeature
#define DYNAMITE_DEFINED_SupportingFeature
#endif


struct Wise2_Exon {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int start;   
    int end;     
    boolean used;   /*  used by some prediction programs etc */ 
    double score;    
    SupportingFeature ** sf;     
    int len;/* len for above sf  */ 
    int maxlen; /* maxlen for above sf */ 
    int phase;   
    } ;  
/* Exon defined */ 
#ifndef DYNAMITE_DEFINED_Exon
typedef struct Wise2_Exon Wise2_Exon;
#define Exon Wise2_Exon
#define DYNAMITE_DEFINED_Exon
#endif


#ifndef DYNAMITE_DEFINED_Translation
typedef struct Wise2_Translation Wise2_Translation;
#define Translation Wise2_Translation
#define DYNAMITE_DEFINED_Translation
#endif

#ifndef DYNAMITE_DEFINED_Gene
typedef struct Wise2_Gene Wise2_Gene;
#define Gene Wise2_Gene
#define DYNAMITE_DEFINED_Gene
#endif

/* Object Transcript
 *
 * Descrip: No Description
 *
 */
struct Wise2_Transcript {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Exon ** exon;    
    int ex_len; /* len for above exon  */ 
    int ex_maxlen;  /* maxlen for above exon */ 
    Gene * parent;   
    Translation ** translation;  
    int len;/* len for above translation  */ 
    int maxlen; /* maxlen for above translation */ 
    cDNA * cDNA;    /*  may not be here! */ 
    } ;  
/* Transcript defined */ 
#ifndef DYNAMITE_DEFINED_Transcript
typedef struct Wise2_Transcript Wise2_Transcript;
#define Transcript Wise2_Transcript
#define DYNAMITE_DEFINED_Transcript
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  copy_Transcript(t)
 *
 * Descrip:    Makes a completely new copy 
 *             of the transcript
 *
 *
 * Arg:        t [UNKN ] Undocumented argument [Transcript *]
 *
 * Return [UNKN ]  Undocumented return value [Transcript *]
 *
 */
Transcript * Wise2_copy_Transcript(Transcript * t);
#define copy_Transcript Wise2_copy_Transcript


/* Function:  get_cDNA_from_Transcript(trs)
 *
 * Descrip:    gets the cDNA associated with this transcript,
 *             if necessary, building it from the exon information
 *             provided.
 *
 *             returns a soft-linked object. If you want to ensure
 *             that this cDNA object remains in memory use
 *             /hard_link_cDNA on the object.
 *
 *
 * Arg:        trs [READ ] transcript to get cDNA from [Transcript *]
 *
 * Return [SOFT ]  cDNA of the transcript [cDNA *]
 *
 */
cDNA * Wise2_get_cDNA_from_Transcript(Transcript * trs);
#define get_cDNA_from_Transcript Wise2_get_cDNA_from_Transcript


/* Function:  length_Transcript(tr)
 *
 * Descrip:    returns the length by looking at the
 *             exon lengths
 *
 *
 * Arg:        tr [UNKN ] Undocumented argument [Transcript *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_length_Transcript(Transcript * tr);
#define length_Transcript Wise2_length_Transcript


/* Function:  show_Transcript(tr,ofp)
 *
 * Descrip:    shows a transcript in vaguely human form
 *
 *
 * Arg:         tr [UNKN ] Undocumented argument [Transcript *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_show_Transcript(Transcript * tr,FILE * ofp);
#define show_Transcript Wise2_show_Transcript


/* Function:  hard_link_SupportingFeature(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SupportingFeature *]
 *
 * Return [UNKN ]  Undocumented return value [SupportingFeature *]
 *
 */
SupportingFeature * Wise2_hard_link_SupportingFeature(SupportingFeature * obj);
#define hard_link_SupportingFeature Wise2_hard_link_SupportingFeature


/* Function:  SupportingFeature_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SupportingFeature *]
 *
 */
SupportingFeature * Wise2_SupportingFeature_alloc(void);
#define SupportingFeature_alloc Wise2_SupportingFeature_alloc


/* Function:  free_SupportingFeature(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SupportingFeature *]
 *
 * Return [UNKN ]  Undocumented return value [SupportingFeature *]
 *
 */
SupportingFeature * Wise2_free_SupportingFeature(SupportingFeature * obj);
#define free_SupportingFeature Wise2_free_SupportingFeature


/* Function:  add_Exon(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Exon *]
 * Arg:        add [OWNER] Object to add to the list [SupportingFeature *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_Exon(Exon * obj,SupportingFeature * add);
#define add_Exon Wise2_add_Exon


/* Function:  flush_Exon(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [Exon *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_Exon(Exon * obj);
#define flush_Exon Wise2_flush_Exon


/* Function:  Exon_alloc_std(void)
 *
 * Descrip:    Equivalent to Exon_alloc_len(ExonLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Exon *]
 *
 */
Exon * Wise2_Exon_alloc_std(void);
#define Exon_alloc_std Wise2_Exon_alloc_std


/* Function:  Exon_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [Exon *]
 *
 */
Exon * Wise2_Exon_alloc_len(int len);
#define Exon_alloc_len Wise2_Exon_alloc_len


/* Function:  hard_link_Exon(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Exon *]
 *
 * Return [UNKN ]  Undocumented return value [Exon *]
 *
 */
Exon * Wise2_hard_link_Exon(Exon * obj);
#define hard_link_Exon Wise2_hard_link_Exon


/* Function:  Exon_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Exon *]
 *
 */
Exon * Wise2_Exon_alloc(void);
#define Exon_alloc Wise2_Exon_alloc


/* Function:  free_Exon(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Exon *]
 *
 * Return [UNKN ]  Undocumented return value [Exon *]
 *
 */
Exon * Wise2_free_Exon(Exon * obj);
#define free_Exon Wise2_free_Exon


/* Function:  add_ex_Transcript(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Transcript *]
 * Arg:        add [OWNER] Object to add to the list [Exon *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_ex_Transcript(Transcript * obj,Exon * add);
#define add_ex_Transcript Wise2_add_ex_Transcript


/* Function:  flush_ex_Transcript(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [Transcript *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_ex_Transcript(Transcript * obj);
#define flush_ex_Transcript Wise2_flush_ex_Transcript


/* Function:  add_Transcript(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Transcript *]
 * Arg:        add [OWNER] Object to add to the list [Translation *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_Transcript(Transcript * obj,Translation * add);
#define add_Transcript Wise2_add_Transcript


/* Function:  flush_Transcript(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [Transcript *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_Transcript(Transcript * obj);
#define flush_Transcript Wise2_flush_Transcript


/* Function:  Transcript_alloc_std(void)
 *
 * Descrip:    Equivalent to Transcript_alloc_len(TranscriptLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Transcript *]
 *
 */
Transcript * Wise2_Transcript_alloc_std(void);
#define Transcript_alloc_std Wise2_Transcript_alloc_std


/* Function:  Transcript_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [Transcript *]
 *
 */
Transcript * Wise2_Transcript_alloc_len(int len);
#define Transcript_alloc_len Wise2_Transcript_alloc_len


/* Function:  hard_link_Transcript(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Transcript *]
 *
 * Return [UNKN ]  Undocumented return value [Transcript *]
 *
 */
Transcript * Wise2_hard_link_Transcript(Transcript * obj);
#define hard_link_Transcript Wise2_hard_link_Transcript


/* Function:  Transcript_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Transcript *]
 *
 */
Transcript * Wise2_Transcript_alloc(void);
#define Transcript_alloc Wise2_Transcript_alloc


/* Function:  free_Transcript(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Transcript *]
 *
 * Return [UNKN ]  Undocumented return value [Transcript *]
 *
 */
Transcript * Wise2_free_Transcript(Transcript * obj);
#define free_Transcript Wise2_free_Transcript


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
cDNA * Wise2_access_cDNA_Transcript(Transcript * obj);
#define access_cDNA_Transcript Wise2_access_cDNA_Transcript
int Wise2_access_end_Exon(Exon * obj);
#define access_end_Exon Wise2_access_end_Exon
boolean Wise2_replace_end_Exon(Exon * obj,int end);
#define replace_end_Exon Wise2_replace_end_Exon
boolean Wise2_replace_used_Exon(Exon * obj,boolean used);
#define replace_used_Exon Wise2_replace_used_Exon
boolean Wise2_replace_cDNA_Transcript(Transcript * obj,cDNA * cDNA);
#define replace_cDNA_Transcript Wise2_replace_cDNA_Transcript
boolean Wise2_access_used_Exon(Exon * obj);
#define access_used_Exon Wise2_access_used_Exon
double Wise2_access_score_Exon(Exon * obj);
#define access_score_Exon Wise2_access_score_Exon
SupportingFeature * Wise2_access_sf_Exon(Exon * obj,int i);
#define access_sf_Exon Wise2_access_sf_Exon
boolean Wise2_replace_start_Exon(Exon * obj,int start);
#define replace_start_Exon Wise2_replace_start_Exon
int Wise2_length_sf_Exon(Exon * obj);
#define length_sf_Exon Wise2_length_sf_Exon
int Wise2_access_start_Exon(Exon * obj);
#define access_start_Exon Wise2_access_start_Exon
boolean Wise2_replace_phase_Exon(Exon * obj,int phase);
#define replace_phase_Exon Wise2_replace_phase_Exon
boolean Wise2_replace_score_Exon(Exon * obj,double score);
#define replace_score_Exon Wise2_replace_score_Exon
int Wise2_access_phase_Exon(Exon * obj);
#define access_phase_Exon Wise2_access_phase_Exon
Gene * Wise2_access_parent_Transcript(Transcript * obj);
#define access_parent_Transcript Wise2_access_parent_Transcript
Exon * Wise2_access_exon_Transcript(Transcript * obj,int i);
#define access_exon_Transcript Wise2_access_exon_Transcript
Translation * Wise2_access_translation_Transcript(Transcript * obj,int i);
#define access_translation_Transcript Wise2_access_translation_Transcript
int Wise2_length_exon_Transcript(Transcript * obj);
#define length_exon_Transcript Wise2_length_exon_Transcript
int Wise2_length_translation_Transcript(Transcript * obj);
#define length_translation_Transcript Wise2_length_translation_Transcript
boolean Wise2_replace_parent_Transcript(Transcript * obj,Gene * parent);
#define replace_parent_Transcript Wise2_replace_parent_Transcript
void Wise2_swap_Exon(SupportingFeature ** list,int i,int j) ;
#define swap_Exon Wise2_swap_Exon
void Wise2_qsort_Exon(SupportingFeature ** list,int left,int right,int (*comp)(SupportingFeature * ,SupportingFeature * ));
#define qsort_Exon Wise2_qsort_Exon
void Wise2_sort_Exon(Exon * obj,int (*comp)(SupportingFeature *, SupportingFeature *));
#define sort_Exon Wise2_sort_Exon
boolean Wise2_expand_Exon(Exon * obj,int len);
#define expand_Exon Wise2_expand_Exon
void Wise2_swap_ex_Transcript(Exon ** list,int i,int j) ;
#define swap_ex_Transcript Wise2_swap_ex_Transcript
void Wise2_qsort_ex_Transcript(Exon ** list,int left,int right,int (*comp)(Exon * ,Exon * ));
#define qsort_ex_Transcript Wise2_qsort_ex_Transcript
void Wise2_sort_ex_Transcript(Transcript * obj,int (*comp)(Exon *, Exon *));
#define sort_ex_Transcript Wise2_sort_ex_Transcript
boolean Wise2_expand_ex_Transcript(Transcript * obj,int len);
#define expand_ex_Transcript Wise2_expand_ex_Transcript
void Wise2_swap_Transcript(Translation ** list,int i,int j) ;
#define swap_Transcript Wise2_swap_Transcript
void Wise2_qsort_Transcript(Translation ** list,int left,int right,int (*comp)(Translation * ,Translation * ));
#define qsort_Transcript Wise2_qsort_Transcript
void Wise2_sort_Transcript(Transcript * obj,int (*comp)(Translation *, Translation *));
#define sort_Transcript Wise2_sort_Transcript
boolean Wise2_expand_Transcript(Transcript * obj,int len);
#define expand_Transcript Wise2_expand_Transcript

#ifdef _cplusplus
}
#endif

#endif
