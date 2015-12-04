#ifndef DYNAMITEmatchsumHEADERFILE
#define DYNAMITEmatchsumHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "gwrap.h"

#define MatchSummarySetLISTLENGTH 128

/* Object MatchSummary
 *
 * Descrip: A Match Summary has summary statistics
 *        for a single alignment, with the
 *        two start/end ranges and the number of
 *        introns and frameshifts for each
 *        sequence (obviously, if one is a protein
 *        then neither are valid!)
 *
 *
 */
struct Wise2_MatchSummary {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    double bits;     
    char * qname;    
    char * tname;    
    int qstart;  
    int qend;    
    int tstart;  
    int tend;    
    int qintron;     
    int qframeshift;     
    int tintron;     
    int tframeshift;     
    } ;  
/* MatchSummary defined */ 
#ifndef DYNAMITE_DEFINED_MatchSummary
typedef struct Wise2_MatchSummary Wise2_MatchSummary;
#define MatchSummary Wise2_MatchSummary
#define DYNAMITE_DEFINED_MatchSummary
#endif


/* Object MatchSummarySet
 *
 * Descrip: This holds a set of MatchSummaries,
 *
 *
 */
struct Wise2_MatchSummarySet {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    MatchSummary ** ms;  
    int len;/* len for above ms  */ 
    int maxlen; /* maxlen for above ms */ 
    } ;  
/* MatchSummarySet defined */ 
#ifndef DYNAMITE_DEFINED_MatchSummarySet
typedef struct Wise2_MatchSummarySet Wise2_MatchSummarySet;
#define MatchSummarySet Wise2_MatchSummarySet
#define DYNAMITE_DEFINED_MatchSummarySet
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  MatchSummarySet_from_AlnBlock_estwise(alb,qname,offset,target)
 *
 * Descrip:    Builds a MatchSummarySet from a
 *             EstWise alignment. this makes
 *             alot of assumptions about the labels
 *             setc in alb, so make sure it was a 
 *             estwise alignment  - however as you
 *             can notice this is exactly the same 
 *             labels as found in genewise set
 *
 *
 * Arg:           alb [UNKN ] Undocumented argument [AlnBlock *]
 * Arg:         qname [UNKN ] Undocumented argument [char *]
 * Arg:        offset [UNKN ] Undocumented argument [int]
 * Arg:        target [UNKN ] Undocumented argument [Sequence *]
 *
 * Return [UNKN ]  Undocumented return value [MatchSummarySet *]
 *
 */
MatchSummarySet * Wise2_MatchSummarySet_from_AlnBlock_estwise(AlnBlock * alb,char * qname,int offset,Sequence * target);
#define MatchSummarySet_from_AlnBlock_estwise Wise2_MatchSummarySet_from_AlnBlock_estwise


/* Function:  MatchSummarySet_from_AlnBlock_genewise(alb,qname,protoff,target)
 *
 * Descrip:    Builds a MatchSummarySet from a
 *             GeneWise alignment. this makes
 *             alot of assumptions about the labels
 *             setc in alb, so make sure it was a 
 *             genewise alignment 
 *
 *
 * Arg:            alb [UNKN ] Undocumented argument [AlnBlock *]
 * Arg:          qname [UNKN ] Undocumented argument [char *]
 * Arg:        protoff [UNKN ] Undocumented argument [int]
 * Arg:         target [UNKN ] Undocumented argument [Sequence *]
 *
 * Return [UNKN ]  Undocumented return value [MatchSummarySet *]
 *
 */
MatchSummarySet * Wise2_MatchSummarySet_from_AlnBlock_genewise(AlnBlock * alb,char * qname,int protoff,Sequence * target);
#define MatchSummarySet_from_AlnBlock_genewise Wise2_MatchSummarySet_from_AlnBlock_genewise


/* Function:  show_MatchSummary_genewise_header(ofp)
 *
 * Descrip:    Shows a header (bits Query start etc) which matches
 *             the order of the show_MatchSummarySet_genewise
 *
 *
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_show_MatchSummary_genewise_header(FILE * ofp);
#define show_MatchSummary_genewise_header Wise2_show_MatchSummary_genewise_header


/* Function:  show_MatchSummarySet_genewise(mss,ofp)
 *
 * Descrip:    shows Matchsummary for genewise
 *             results, ie with target indels and introns
 *
 *
 * Arg:        mss [UNKN ] Undocumented argument [MatchSummarySet *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_show_MatchSummarySet_genewise(MatchSummarySet * mss,FILE * ofp);
#define show_MatchSummarySet_genewise Wise2_show_MatchSummarySet_genewise


/* Function:  show_MatchSummary_estwise_header(ofp)
 *
 * Descrip:    Shows a header (bits Query start etc) which matches
 *             the order of the show_MatchSummarySet_estwise
 *
 *
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_show_MatchSummary_estwise_header(FILE * ofp);
#define show_MatchSummary_estwise_header Wise2_show_MatchSummary_estwise_header


/* Function:  show_MatchSummarySet_estwise(mss,ofp)
 *
 * Descrip:    Shows an estwise match, ie with only the target with indels
 *
 *
 * Arg:        mss [UNKN ] Undocumented argument [MatchSummarySet *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_show_MatchSummarySet_estwise(MatchSummarySet * mss,FILE * ofp);
#define show_MatchSummarySet_estwise Wise2_show_MatchSummarySet_estwise


/* Function:  hard_link_MatchSummary(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [MatchSummary *]
 *
 * Return [UNKN ]  Undocumented return value [MatchSummary *]
 *
 */
MatchSummary * Wise2_hard_link_MatchSummary(MatchSummary * obj);
#define hard_link_MatchSummary Wise2_hard_link_MatchSummary


/* Function:  MatchSummary_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [MatchSummary *]
 *
 */
MatchSummary * Wise2_MatchSummary_alloc(void);
#define MatchSummary_alloc Wise2_MatchSummary_alloc


/* Function:  free_MatchSummary(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [MatchSummary *]
 *
 * Return [UNKN ]  Undocumented return value [MatchSummary *]
 *
 */
MatchSummary * Wise2_free_MatchSummary(MatchSummary * obj);
#define free_MatchSummary Wise2_free_MatchSummary


/* Function:  add_MatchSummarySet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [MatchSummarySet *]
 * Arg:        add [OWNER] Object to add to the list [MatchSummary *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_MatchSummarySet(MatchSummarySet * obj,MatchSummary * add);
#define add_MatchSummarySet Wise2_add_MatchSummarySet


/* Function:  flush_MatchSummarySet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [MatchSummarySet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_MatchSummarySet(MatchSummarySet * obj);
#define flush_MatchSummarySet Wise2_flush_MatchSummarySet


/* Function:  MatchSummarySet_alloc_std(void)
 *
 * Descrip:    Equivalent to MatchSummarySet_alloc_len(MatchSummarySetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [MatchSummarySet *]
 *
 */
MatchSummarySet * Wise2_MatchSummarySet_alloc_std(void);
#define MatchSummarySet_alloc_std Wise2_MatchSummarySet_alloc_std


/* Function:  MatchSummarySet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [MatchSummarySet *]
 *
 */
MatchSummarySet * Wise2_MatchSummarySet_alloc_len(int len);
#define MatchSummarySet_alloc_len Wise2_MatchSummarySet_alloc_len


/* Function:  hard_link_MatchSummarySet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [MatchSummarySet *]
 *
 * Return [UNKN ]  Undocumented return value [MatchSummarySet *]
 *
 */
MatchSummarySet * Wise2_hard_link_MatchSummarySet(MatchSummarySet * obj);
#define hard_link_MatchSummarySet Wise2_hard_link_MatchSummarySet


/* Function:  MatchSummarySet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [MatchSummarySet *]
 *
 */
MatchSummarySet * Wise2_MatchSummarySet_alloc(void);
#define MatchSummarySet_alloc Wise2_MatchSummarySet_alloc


/* Function:  free_MatchSummarySet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [MatchSummarySet *]
 *
 * Return [UNKN ]  Undocumented return value [MatchSummarySet *]
 *
 */
MatchSummarySet * Wise2_free_MatchSummarySet(MatchSummarySet * obj);
#define free_MatchSummarySet Wise2_free_MatchSummarySet


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
MatchSummary * Wise2_access_ms_MatchSummarySet(MatchSummarySet * obj,int i);
#define access_ms_MatchSummarySet Wise2_access_ms_MatchSummarySet
boolean Wise2_replace_qstart_MatchSummary(MatchSummary * obj,int qstart);
#define replace_qstart_MatchSummary Wise2_replace_qstart_MatchSummary
int Wise2_access_qstart_MatchSummary(MatchSummary * obj);
#define access_qstart_MatchSummary Wise2_access_qstart_MatchSummary
int Wise2_access_qend_MatchSummary(MatchSummary * obj);
#define access_qend_MatchSummary Wise2_access_qend_MatchSummary
int Wise2_access_tframeshift_MatchSummary(MatchSummary * obj);
#define access_tframeshift_MatchSummary Wise2_access_tframeshift_MatchSummary
boolean Wise2_replace_tstart_MatchSummary(MatchSummary * obj,int tstart);
#define replace_tstart_MatchSummary Wise2_replace_tstart_MatchSummary
int Wise2_length_ms_MatchSummarySet(MatchSummarySet * obj);
#define length_ms_MatchSummarySet Wise2_length_ms_MatchSummarySet
int Wise2_access_tstart_MatchSummary(MatchSummary * obj);
#define access_tstart_MatchSummary Wise2_access_tstart_MatchSummary
double Wise2_access_bits_MatchSummary(MatchSummary * obj);
#define access_bits_MatchSummary Wise2_access_bits_MatchSummary
boolean Wise2_replace_tend_MatchSummary(MatchSummary * obj,int tend);
#define replace_tend_MatchSummary Wise2_replace_tend_MatchSummary
char * Wise2_access_qname_MatchSummary(MatchSummary * obj);
#define access_qname_MatchSummary Wise2_access_qname_MatchSummary
int Wise2_access_tend_MatchSummary(MatchSummary * obj);
#define access_tend_MatchSummary Wise2_access_tend_MatchSummary
char * Wise2_access_tname_MatchSummary(MatchSummary * obj);
#define access_tname_MatchSummary Wise2_access_tname_MatchSummary
boolean Wise2_replace_qintron_MatchSummary(MatchSummary * obj,int qintron);
#define replace_qintron_MatchSummary Wise2_replace_qintron_MatchSummary
boolean Wise2_replace_qend_MatchSummary(MatchSummary * obj,int qend);
#define replace_qend_MatchSummary Wise2_replace_qend_MatchSummary
int Wise2_access_qintron_MatchSummary(MatchSummary * obj);
#define access_qintron_MatchSummary Wise2_access_qintron_MatchSummary
boolean Wise2_replace_bits_MatchSummary(MatchSummary * obj,double bits);
#define replace_bits_MatchSummary Wise2_replace_bits_MatchSummary
boolean Wise2_replace_qframeshift_MatchSummary(MatchSummary * obj,int qframeshift);
#define replace_qframeshift_MatchSummary Wise2_replace_qframeshift_MatchSummary
boolean Wise2_replace_tname_MatchSummary(MatchSummary * obj,char * tname);
#define replace_tname_MatchSummary Wise2_replace_tname_MatchSummary
int Wise2_access_qframeshift_MatchSummary(MatchSummary * obj);
#define access_qframeshift_MatchSummary Wise2_access_qframeshift_MatchSummary
boolean Wise2_replace_qname_MatchSummary(MatchSummary * obj,char * qname);
#define replace_qname_MatchSummary Wise2_replace_qname_MatchSummary
boolean Wise2_replace_tintron_MatchSummary(MatchSummary * obj,int tintron);
#define replace_tintron_MatchSummary Wise2_replace_tintron_MatchSummary
boolean Wise2_replace_tframeshift_MatchSummary(MatchSummary * obj,int tframeshift);
#define replace_tframeshift_MatchSummary Wise2_replace_tframeshift_MatchSummary
int Wise2_access_tintron_MatchSummary(MatchSummary * obj);
#define access_tintron_MatchSummary Wise2_access_tintron_MatchSummary
MatchSummary * Wise2_MatchSummary_from_AlnColumn_genewise(AlnColumn * alc,AlnColumn ** end);
#define MatchSummary_from_AlnColumn_genewise Wise2_MatchSummary_from_AlnColumn_genewise
void Wise2_show_MatchSummary_genewise(MatchSummary * ms,FILE * ofp);
#define show_MatchSummary_genewise Wise2_show_MatchSummary_genewise
void Wise2_show_MatchSummary_estwise(MatchSummary * ms,FILE * ofp);
#define show_MatchSummary_estwise Wise2_show_MatchSummary_estwise
void Wise2_swap_MatchSummarySet(MatchSummary ** list,int i,int j) ;
#define swap_MatchSummarySet Wise2_swap_MatchSummarySet
void Wise2_qsort_MatchSummarySet(MatchSummary ** list,int left,int right,int (*comp)(MatchSummary * ,MatchSummary * ));
#define qsort_MatchSummarySet Wise2_qsort_MatchSummarySet
void Wise2_sort_MatchSummarySet(MatchSummarySet * obj,int (*comp)(MatchSummary *, MatchSummary *));
#define sort_MatchSummarySet Wise2_sort_MatchSummarySet
boolean Wise2_expand_MatchSummarySet(MatchSummarySet * obj,int len);
#define expand_MatchSummarySet Wise2_expand_MatchSummarySet

#ifdef _cplusplus
}
#endif

#endif
