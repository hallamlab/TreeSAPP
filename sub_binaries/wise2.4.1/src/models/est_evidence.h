#ifndef DYNAMITEest_evidenceHEADERFILE
#define DYNAMITEest_evidenceHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "genome_evidence.h"
#include "aln.h"

#define EstEvidenceLISTLENGTH 24

struct Wise2_EstExon {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int start;   
    int end;     
    boolean is_coding;   
    int phase;   
    int used;    
    int intron_3_score;  
    } ;  
/* EstExon defined */ 
#ifndef DYNAMITE_DEFINED_EstExon
typedef struct Wise2_EstExon Wise2_EstExon;
#define EstExon Wise2_EstExon
#define DYNAMITE_DEFINED_EstExon
#endif


struct Wise2_EstIndel {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int start;   
    int end;     
    } ;  
/* EstIndel defined */ 
#ifndef DYNAMITE_DEFINED_EstIndel
typedef struct Wise2_EstIndel Wise2_EstIndel;
#define EstIndel Wise2_EstIndel
#define DYNAMITE_DEFINED_EstIndel
#endif


struct Wise2_EstEvidence {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    EstExon    ** exon;  
    int len;/* len for above exon  */ 
    int maxlen; /* maxlen for above exon */ 
    EstIndel   ** indel;     
    int indel_len;  /* len for above indel  */ 
    int indel_maxlen;   /* maxlen for above indel */ 
    CodonTable * ct;     
    int in_smell;    
    } ;  
/* EstEvidence defined */ 
#ifndef DYNAMITE_DEFINED_EstEvidence
typedef struct Wise2_EstEvidence Wise2_EstEvidence;
#define EstEvidence Wise2_EstEvidence
#define DYNAMITE_DEFINED_EstEvidence
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_EstExon(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [EstExon *]
 *
 * Return [UNKN ]  Undocumented return value [EstExon *]
 *
 */
EstExon * Wise2_hard_link_EstExon(EstExon * obj);
#define hard_link_EstExon Wise2_hard_link_EstExon


/* Function:  EstExon_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [EstExon *]
 *
 */
EstExon * Wise2_EstExon_alloc(void);
#define EstExon_alloc Wise2_EstExon_alloc


/* Function:  free_EstExon(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [EstExon *]
 *
 * Return [UNKN ]  Undocumented return value [EstExon *]
 *
 */
EstExon * Wise2_free_EstExon(EstExon * obj);
#define free_EstExon Wise2_free_EstExon


/* Function:  hard_link_EstIndel(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [EstIndel *]
 *
 * Return [UNKN ]  Undocumented return value [EstIndel *]
 *
 */
EstIndel * Wise2_hard_link_EstIndel(EstIndel * obj);
#define hard_link_EstIndel Wise2_hard_link_EstIndel


/* Function:  EstIndel_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [EstIndel *]
 *
 */
EstIndel * Wise2_EstIndel_alloc(void);
#define EstIndel_alloc Wise2_EstIndel_alloc


/* Function:  free_EstIndel(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [EstIndel *]
 *
 * Return [UNKN ]  Undocumented return value [EstIndel *]
 *
 */
EstIndel * Wise2_free_EstIndel(EstIndel * obj);
#define free_EstIndel Wise2_free_EstIndel


/* Function:  add_EstEvidence(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [EstEvidence *]
 * Arg:        add [OWNER] Object to add to the list [EstExon    *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_EstEvidence(EstEvidence * obj,EstExon    * add);
#define add_EstEvidence Wise2_add_EstEvidence


/* Function:  flush_EstEvidence(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [EstEvidence *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_EstEvidence(EstEvidence * obj);
#define flush_EstEvidence Wise2_flush_EstEvidence


/* Function:  add_indel_EstEvidence(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [EstEvidence *]
 * Arg:        add [OWNER] Object to add to the list [EstIndel   *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_indel_EstEvidence(EstEvidence * obj,EstIndel   * add);
#define add_indel_EstEvidence Wise2_add_indel_EstEvidence


/* Function:  flush_indel_EstEvidence(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [EstEvidence *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_indel_EstEvidence(EstEvidence * obj);
#define flush_indel_EstEvidence Wise2_flush_indel_EstEvidence


/* Function:  EstEvidence_alloc_std(void)
 *
 * Descrip:    Equivalent to EstEvidence_alloc_len(EstEvidenceLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [EstEvidence *]
 *
 */
EstEvidence * Wise2_EstEvidence_alloc_std(void);
#define EstEvidence_alloc_std Wise2_EstEvidence_alloc_std


/* Function:  EstEvidence_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [EstEvidence *]
 *
 */
EstEvidence * Wise2_EstEvidence_alloc_len(int len);
#define EstEvidence_alloc_len Wise2_EstEvidence_alloc_len


/* Function:  hard_link_EstEvidence(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [EstEvidence *]
 *
 * Return [UNKN ]  Undocumented return value [EstEvidence *]
 *
 */
EstEvidence * Wise2_hard_link_EstEvidence(EstEvidence * obj);
#define hard_link_EstEvidence Wise2_hard_link_EstEvidence


/* Function:  EstEvidence_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [EstEvidence *]
 *
 */
EstEvidence * Wise2_EstEvidence_alloc(void);
#define EstEvidence_alloc Wise2_EstEvidence_alloc


/* Function:  free_EstEvidence(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [EstEvidence *]
 *
 * Return [UNKN ]  Undocumented return value [EstEvidence *]
 *
 */
EstEvidence * Wise2_free_EstEvidence(EstEvidence * obj);
#define free_EstEvidence Wise2_free_EstEvidence


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
int Wise2_indicate_intron_used(GenomeEvidenceSet * set,AlnBlock * alb);
#define indicate_intron_used Wise2_indicate_intron_used
GenomeEvidenceSet * Wise2_read_est_evidence(FILE * ifp,CodonTable * ct);
#define read_est_evidence Wise2_read_est_evidence
GenomeEvidenceUnit * Wise2_new_est_GenomeEvidenceUnit(EstEvidence * evi);
#define new_est_GenomeEvidenceUnit Wise2_new_est_GenomeEvidenceUnit
int Wise2_est_utr5_start(void * data,ComplexSequence *seq,int jposition);
#define est_utr5_start Wise2_est_utr5_start
int Wise2_est_utr3_end(void * data,ComplexSequence *seq,int jposition);
#define est_utr3_end Wise2_est_utr3_end
int Wise2_est_start_pot(void * data,ComplexSequence *seq,int jposition);
#define est_start_pot Wise2_est_start_pot
int Wise2_est_stop_pot(void * data,ComplexSequence *seq,int jposition);
#define est_stop_pot Wise2_est_stop_pot
int Wise2_est_cds_frameshift(void * data,ComplexSequence * seq,int jposition,int jump);
#define est_cds_frameshift Wise2_est_cds_frameshift
int Wise2_est_cds_3SS(void * data,ComplexSequence *seq,int jposition,int phase);
#define est_cds_3SS Wise2_est_cds_3SS
int Wise2_est_cds_5SS(void * data,ComplexSequence *seq,int jposition,int phase);
#define est_cds_5SS Wise2_est_cds_5SS
int Wise2_est_intron_pot(void * data,ComplexSequence *seq,int jposition);
#define est_intron_pot Wise2_est_intron_pot
int Wise2_est_cds_pot(void * data,ComplexSequence *seq,int jposition);
#define est_cds_pot Wise2_est_cds_pot
int Wise2_est_3ss(void * data,ComplexSequence *seq,int jposition);
#define est_3ss Wise2_est_3ss
int Wise2_est_5ss(void * data,ComplexSequence * seq,int jposition);
#define est_5ss Wise2_est_5ss
int Wise2_est_utr_pot(void * data,ComplexSequence *seq,int jposition);
#define est_utr_pot Wise2_est_utr_pot


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_EstEvidence(EstExon    ** list,int i,int j) ;
#define swap_EstEvidence Wise2_swap_EstEvidence
void Wise2_qsort_EstEvidence(EstExon    ** list,int left,int right,int (*comp)(EstExon    * ,EstExon    * ));
#define qsort_EstEvidence Wise2_qsort_EstEvidence
void Wise2_sort_EstEvidence(EstEvidence * obj,int (*comp)(EstExon    *, EstExon    *));
#define sort_EstEvidence Wise2_sort_EstEvidence
boolean Wise2_expand_EstEvidence(EstEvidence * obj,int len);
#define expand_EstEvidence Wise2_expand_EstEvidence
void Wise2_swap_indel_EstEvidence(EstIndel   ** list,int i,int j) ;
#define swap_indel_EstEvidence Wise2_swap_indel_EstEvidence
void Wise2_qsort_indel_EstEvidence(EstIndel   ** list,int left,int right,int (*comp)(EstIndel   * ,EstIndel   * ));
#define qsort_indel_EstEvidence Wise2_qsort_indel_EstEvidence
void Wise2_sort_indel_EstEvidence(EstEvidence * obj,int (*comp)(EstIndel   *, EstIndel   *));
#define sort_indel_EstEvidence Wise2_sort_indel_EstEvidence
boolean Wise2_expand_indel_EstEvidence(EstEvidence * obj,int len);
#define expand_indel_EstEvidence Wise2_expand_indel_EstEvidence

#ifdef _cplusplus
}
#endif

#endif
