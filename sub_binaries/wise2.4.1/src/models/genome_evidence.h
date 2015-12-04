#ifndef DYNAMITEgenome_evidenceHEADERFILE
#define DYNAMITEgenome_evidenceHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"


#define GNE_CDS_3SS(evi,i,seq,j,phase) ((*((evi)->geu[i]->cds_3SS))((evi)->geu[i]->data,seq,j,phase))
#define GNE_CDS_5SS(evi,i,seq,j,phase) ((*((evi)->geu[i]->cds_5SS))((evi)->geu[i]->data,seq,j,phase))
#define GNE_CDS(evi,i,seq,j)           ((*((evi)->geu[i]->cds_pot))((evi)->geu[i]->data,seq,j))
#define GNE_CDS_INTRON(evi,i,seq,j)    ((*((evi)->geu[i]->cds_intron_pot))((evi)->geu[i]->data,seq,j)) 
#define GNE_UTR_3SS(evi,i,seq,j)       ((*((evi)->geu[i]->utr_3SS))((evi)->geu[i]->data,seq,j))
#define GNE_UTR_5SS(evi,i,seq,j)       ((*((evi)->geu[i]->utr_5SS))((evi)->geu[i]->data,seq,j))
#define GNE_UTR(evi,i,seq,j)           ((*((evi)->geu[i]->utr_pot))((evi)->geu[i]->data,seq,j))
#define GNE_UTR_INTRON(evi,i,seq,j)    ((*((evi)->geu[i]->utr_intron_pot))((evi)->geu[i]->data,seq,j))
#define GNE_CDS_FRAMESHIFT(evi,i,seq,j,jump) ((*((evi)->geu[i]->frameshift_cds))((evi)->geu[i]->data,seq,j,jump))
#define GNE_START_CODON(evi,i,seq,j)   ((*((evi)->geu[i]->start_pot))((evi)->geu[i]->data,seq,j))
#define GNE_STOP_CODON(evi,i,seq,j)    ((*((evi)->geu[i]->stop_pot))((evi)->geu[i]->data,seq,j))
#define GNE_UTR3_END(evi,i,seq,j)      ((*((evi)->geu[i]->utr3_end))((evi)->geu[i]->data,seq,j))
#define GNE_UTR5_START(evi,i,seq,j)      ((*((evi)->geu[i]->utr5_start))((evi)->geu[i]->data,seq,j))

#define GenomeEvidenceSetLISTLENGTH 512

struct Wise2_GenomeEvidenceUnit {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int type;    
    void * data;     
    int (*cds_3SS)(void*,Wise2_ComplexSequence*,int,int);    
    int (*cds_5SS)(void*,Wise2_ComplexSequence*,int,int);    
    int (*utr_3SS)(void*,Wise2_ComplexSequence*,int);    
    int (*utr_5SS)(void*,Wise2_ComplexSequence*,int);    
    int (*cds_pot)(void*,Wise2_ComplexSequence*,int);    
    int (*utr_pot)(void*,Wise2_ComplexSequence*,int);    
    int (*cds_intron_pot)(void*,Wise2_ComplexSequence*,int); 
    int (*utr_intron_pot)(void*,Wise2_ComplexSequence*,int); 
    int (*frameshift_cds)(void*,Wise2_ComplexSequence*,int,int); 
    int (*start_pot)(void*,Wise2_ComplexSequence*,int);  
    int (*stop_pot)(void*,Wise2_ComplexSequence*,int);   
    int (*utr3_end)(void*,Wise2_ComplexSequence*,int);   
    int (*utr5_start)(void*,Wise2_ComplexSequence*,int); 
    int (*geu_free)(void *); 
    } ;  
/* GenomeEvidenceUnit defined */ 
#ifndef DYNAMITE_DEFINED_GenomeEvidenceUnit
typedef struct Wise2_GenomeEvidenceUnit Wise2_GenomeEvidenceUnit;
#define GenomeEvidenceUnit Wise2_GenomeEvidenceUnit
#define DYNAMITE_DEFINED_GenomeEvidenceUnit
#endif


struct Wise2_GenomeEvidenceSet {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    GenomeEvidenceUnit ** geu;   
    int len;/* len for above geu  */ 
    int maxlen; /* maxlen for above geu */ 
    } ;  
/* GenomeEvidenceSet defined */ 
#ifndef DYNAMITE_DEFINED_GenomeEvidenceSet
typedef struct Wise2_GenomeEvidenceSet Wise2_GenomeEvidenceSet;
#define GenomeEvidenceSet Wise2_GenomeEvidenceSet
#define DYNAMITE_DEFINED_GenomeEvidenceSet
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  free_GenomeEvidenceUnit(obj)
 *
 * Descrip:    Specialised deconstructor. Ensures the 
 *             data structures are freed
 *
 *
 * Arg:        obj [UNKN ] Undocumented argument [GenomeEvidenceUnit *]
 *
 * Return [UNKN ]  Undocumented return value [GenomeEvidenceUnit *]
 *
 */
GenomeEvidenceUnit * Wise2_free_GenomeEvidenceUnit(GenomeEvidenceUnit * obj);
#define free_GenomeEvidenceUnit Wise2_free_GenomeEvidenceUnit


/* Function:  hard_link_GenomeEvidenceUnit(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GenomeEvidenceUnit *]
 *
 * Return [UNKN ]  Undocumented return value [GenomeEvidenceUnit *]
 *
 */
GenomeEvidenceUnit * Wise2_hard_link_GenomeEvidenceUnit(GenomeEvidenceUnit * obj);
#define hard_link_GenomeEvidenceUnit Wise2_hard_link_GenomeEvidenceUnit


/* Function:  GenomeEvidenceUnit_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GenomeEvidenceUnit *]
 *
 */
GenomeEvidenceUnit * Wise2_GenomeEvidenceUnit_alloc(void);
#define GenomeEvidenceUnit_alloc Wise2_GenomeEvidenceUnit_alloc


/* Function:  add_GenomeEvidenceSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GenomeEvidenceSet *]
 * Arg:        add [OWNER] Object to add to the list [GenomeEvidenceUnit *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_GenomeEvidenceSet(GenomeEvidenceSet * obj,GenomeEvidenceUnit * add);
#define add_GenomeEvidenceSet Wise2_add_GenomeEvidenceSet


/* Function:  flush_GenomeEvidenceSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [GenomeEvidenceSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_GenomeEvidenceSet(GenomeEvidenceSet * obj);
#define flush_GenomeEvidenceSet Wise2_flush_GenomeEvidenceSet


/* Function:  GenomeEvidenceSet_alloc_std(void)
 *
 * Descrip:    Equivalent to GenomeEvidenceSet_alloc_len(GenomeEvidenceSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GenomeEvidenceSet *]
 *
 */
GenomeEvidenceSet * Wise2_GenomeEvidenceSet_alloc_std(void);
#define GenomeEvidenceSet_alloc_std Wise2_GenomeEvidenceSet_alloc_std


/* Function:  GenomeEvidenceSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [GenomeEvidenceSet *]
 *
 */
GenomeEvidenceSet * Wise2_GenomeEvidenceSet_alloc_len(int len);
#define GenomeEvidenceSet_alloc_len Wise2_GenomeEvidenceSet_alloc_len


/* Function:  hard_link_GenomeEvidenceSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GenomeEvidenceSet *]
 *
 * Return [UNKN ]  Undocumented return value [GenomeEvidenceSet *]
 *
 */
GenomeEvidenceSet * Wise2_hard_link_GenomeEvidenceSet(GenomeEvidenceSet * obj);
#define hard_link_GenomeEvidenceSet Wise2_hard_link_GenomeEvidenceSet


/* Function:  GenomeEvidenceSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GenomeEvidenceSet *]
 *
 */
GenomeEvidenceSet * Wise2_GenomeEvidenceSet_alloc(void);
#define GenomeEvidenceSet_alloc Wise2_GenomeEvidenceSet_alloc


/* Function:  free_GenomeEvidenceSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GenomeEvidenceSet *]
 *
 * Return [UNKN ]  Undocumented return value [GenomeEvidenceSet *]
 *
 */
GenomeEvidenceSet * Wise2_free_GenomeEvidenceSet(GenomeEvidenceSet * obj);
#define free_GenomeEvidenceSet Wise2_free_GenomeEvidenceSet


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_GenomeEvidenceSet(GenomeEvidenceUnit ** list,int i,int j) ;
#define swap_GenomeEvidenceSet Wise2_swap_GenomeEvidenceSet
void Wise2_qsort_GenomeEvidenceSet(GenomeEvidenceUnit ** list,int left,int right,int (*comp)(GenomeEvidenceUnit * ,GenomeEvidenceUnit * ));
#define qsort_GenomeEvidenceSet Wise2_qsort_GenomeEvidenceSet
void Wise2_sort_GenomeEvidenceSet(GenomeEvidenceSet * obj,int (*comp)(GenomeEvidenceUnit *, GenomeEvidenceUnit *));
#define sort_GenomeEvidenceSet Wise2_sort_GenomeEvidenceSet
boolean Wise2_expand_GenomeEvidenceSet(GenomeEvidenceSet * obj,int len);
#define expand_GenomeEvidenceSet Wise2_expand_GenomeEvidenceSet

#ifdef _cplusplus
}
#endif

#endif
