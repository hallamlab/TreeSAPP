#ifndef DYNAMITEassemblyHEADERFILE
#define DYNAMITEassemblyHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "sequence.h"
#include "largeseqreader.h"

#define AssemblySequenceLISTLENGTH 8
#define AssemblyContigLISTLENGTH 512
#define AssemblyLISTLENGTH 1024

#define AssemblyOpaqueTypeSetLISTLENGTH 16

struct Wise2_AssemblySequenceEdit {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char type;   
    char base;   
    } ;  
/* AssemblySequenceEdit defined */ 
#ifndef DYNAMITE_DEFINED_AssemblySequenceEdit
typedef struct Wise2_AssemblySequenceEdit Wise2_AssemblySequenceEdit;
#define AssemblySequenceEdit Wise2_AssemblySequenceEdit
#define DYNAMITE_DEFINED_AssemblySequenceEdit
#endif


#ifndef DYNAMITE_DEFINED_AssemblySequence
typedef struct Wise2_AssemblySequence Wise2_AssemblySequence;
#define AssemblySequence Wise2_AssemblySequence
#define DYNAMITE_DEFINED_AssemblySequence
#endif

struct Wise2_AssemblyOpaqueType {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int int_type;    
    char base;   
    } ;  
/* AssemblyOpaqueType defined */ 
#ifndef DYNAMITE_DEFINED_AssemblyOpaqueType
typedef struct Wise2_AssemblyOpaqueType Wise2_AssemblyOpaqueType;
#define AssemblyOpaqueType Wise2_AssemblyOpaqueType
#define DYNAMITE_DEFINED_AssemblyOpaqueType
#endif


struct Wise2_AssemblyOpaqueTypeSet {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    AssemblyOpaqueType ** type;  
    int len;/* len for above type  */ 
    int maxlen; /* maxlen for above type */ 
    } ;  
/* AssemblyOpaqueTypeSet defined */ 
#ifndef DYNAMITE_DEFINED_AssemblyOpaqueTypeSet
typedef struct Wise2_AssemblyOpaqueTypeSet Wise2_AssemblyOpaqueTypeSet;
#define AssemblyOpaqueTypeSet Wise2_AssemblyOpaqueTypeSet
#define DYNAMITE_DEFINED_AssemblyOpaqueTypeSet
#endif


struct Wise2_AssemblyOpaqueFeature {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int start;   
    int length;  
    AssemblyOpaqueType * type;   
    } ;  
/* AssemblyOpaqueFeature defined */ 
#ifndef DYNAMITE_DEFINED_AssemblyOpaqueFeature
typedef struct Wise2_AssemblyOpaqueFeature Wise2_AssemblyOpaqueFeature;
#define AssemblyOpaqueFeature Wise2_AssemblyOpaqueFeature
#define DYNAMITE_DEFINED_AssemblyOpaqueFeature
#endif


struct Wise2_AssemblySequence {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Sequence * seq;  
    double * quality;    
    AssemblySequenceEdit ** edit;    
    int len;/* len for above edit  */ 
    int maxlen; /* maxlen for above edit */ 
    char state;  
    AssemblySequence * mirror;   
    AssemblySequence * pair;     
    char mirror_seq;     
    AssemblyOpaqueFeature ** opaque;     
    int opq_len;/* len for above opaque  */ 
    int opq_maxlen; /* maxlen for above opaque */ 
    Sequence * orig;     
    int * abrev_repeat;  
    } ;  
/* AssemblySequence defined */ 
#ifndef DYNAMITE_DEFINED_AssemblySequence
typedef struct Wise2_AssemblySequence Wise2_AssemblySequence;
#define AssemblySequence Wise2_AssemblySequence
#define DYNAMITE_DEFINED_AssemblySequence
#endif


struct Wise2_AssemblySequencePlacement {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    AssemblySequence * aseq;     
    long contig_start;   
    long contig_end;     
    int seq_start;   
    int seq_end;     
    boolean ungapped;    
    } ;  
/* AssemblySequencePlacement defined */ 
#ifndef DYNAMITE_DEFINED_AssemblySequencePlacement
typedef struct Wise2_AssemblySequencePlacement Wise2_AssemblySequencePlacement;
#define AssemblySequencePlacement Wise2_AssemblySequencePlacement
#define DYNAMITE_DEFINED_AssemblySequencePlacement
#endif


struct Wise2_AssemblyContig {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Sequence * consensus;    
    AssemblySequencePlacement ** reads;  
    int len;/* len for above reads  */ 
    int maxlen; /* maxlen for above reads */ 
    boolean clean_start;     
    boolean clean_end;   
    int max_depth;   
    } ;  
/* AssemblyContig defined */ 
#ifndef DYNAMITE_DEFINED_AssemblyContig
typedef struct Wise2_AssemblyContig Wise2_AssemblyContig;
#define AssemblyContig Wise2_AssemblyContig
#define DYNAMITE_DEFINED_AssemblyContig
#endif


struct Wise2_Assembly {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    AssemblyContig ** contig;    
    int len;/* len for above contig  */ 
    int maxlen; /* maxlen for above contig */ 
    AssemblySequence ** chaff;   
    int chaff_len;  /* len for above chaff  */ 
    int chaff_maxlen;   /* maxlen for above chaff */ 
    } ;  
/* Assembly defined */ 
#ifndef DYNAMITE_DEFINED_Assembly
typedef struct Wise2_Assembly Wise2_Assembly;
#define Assembly Wise2_Assembly
#define DYNAMITE_DEFINED_Assembly
#endif


struct Wise2_AssemblyOutputPara {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int min_size;    
    } ;  
/* AssemblyOutputPara defined */ 
#ifndef DYNAMITE_DEFINED_AssemblyOutputPara
typedef struct Wise2_AssemblyOutputPara Wise2_AssemblyOutputPara;
#define AssemblyOutputPara Wise2_AssemblyOutputPara
#define DYNAMITE_DEFINED_AssemblyOutputPara
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_AssemblySequenceEdit(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AssemblySequenceEdit *]
 *
 * Return [UNKN ]  Undocumented return value [AssemblySequenceEdit *]
 *
 */
AssemblySequenceEdit * Wise2_hard_link_AssemblySequenceEdit(AssemblySequenceEdit * obj);
#define hard_link_AssemblySequenceEdit Wise2_hard_link_AssemblySequenceEdit


/* Function:  AssemblySequenceEdit_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AssemblySequenceEdit *]
 *
 */
AssemblySequenceEdit * Wise2_AssemblySequenceEdit_alloc(void);
#define AssemblySequenceEdit_alloc Wise2_AssemblySequenceEdit_alloc


/* Function:  free_AssemblySequenceEdit(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AssemblySequenceEdit *]
 *
 * Return [UNKN ]  Undocumented return value [AssemblySequenceEdit *]
 *
 */
AssemblySequenceEdit * Wise2_free_AssemblySequenceEdit(AssemblySequenceEdit * obj);
#define free_AssemblySequenceEdit Wise2_free_AssemblySequenceEdit


/* Function:  hard_link_AssemblyOpaqueType(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AssemblyOpaqueType *]
 *
 * Return [UNKN ]  Undocumented return value [AssemblyOpaqueType *]
 *
 */
AssemblyOpaqueType * Wise2_hard_link_AssemblyOpaqueType(AssemblyOpaqueType * obj);
#define hard_link_AssemblyOpaqueType Wise2_hard_link_AssemblyOpaqueType


/* Function:  AssemblyOpaqueType_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AssemblyOpaqueType *]
 *
 */
AssemblyOpaqueType * Wise2_AssemblyOpaqueType_alloc(void);
#define AssemblyOpaqueType_alloc Wise2_AssemblyOpaqueType_alloc


/* Function:  free_AssemblyOpaqueType(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AssemblyOpaqueType *]
 *
 * Return [UNKN ]  Undocumented return value [AssemblyOpaqueType *]
 *
 */
AssemblyOpaqueType * Wise2_free_AssemblyOpaqueType(AssemblyOpaqueType * obj);
#define free_AssemblyOpaqueType Wise2_free_AssemblyOpaqueType


/* Function:  add_AssemblyOpaqueTypeSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [AssemblyOpaqueTypeSet *]
 * Arg:        add [OWNER] Object to add to the list [AssemblyOpaqueType *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_AssemblyOpaqueTypeSet(AssemblyOpaqueTypeSet * obj,AssemblyOpaqueType * add);
#define add_AssemblyOpaqueTypeSet Wise2_add_AssemblyOpaqueTypeSet


/* Function:  flush_AssemblyOpaqueTypeSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [AssemblyOpaqueTypeSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_AssemblyOpaqueTypeSet(AssemblyOpaqueTypeSet * obj);
#define flush_AssemblyOpaqueTypeSet Wise2_flush_AssemblyOpaqueTypeSet


/* Function:  AssemblyOpaqueTypeSet_alloc_std(void)
 *
 * Descrip:    Equivalent to AssemblyOpaqueTypeSet_alloc_len(AssemblyOpaqueTypeSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AssemblyOpaqueTypeSet *]
 *
 */
AssemblyOpaqueTypeSet * Wise2_AssemblyOpaqueTypeSet_alloc_std(void);
#define AssemblyOpaqueTypeSet_alloc_std Wise2_AssemblyOpaqueTypeSet_alloc_std


/* Function:  AssemblyOpaqueTypeSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [AssemblyOpaqueTypeSet *]
 *
 */
AssemblyOpaqueTypeSet * Wise2_AssemblyOpaqueTypeSet_alloc_len(int len);
#define AssemblyOpaqueTypeSet_alloc_len Wise2_AssemblyOpaqueTypeSet_alloc_len


/* Function:  hard_link_AssemblyOpaqueTypeSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AssemblyOpaqueTypeSet *]
 *
 * Return [UNKN ]  Undocumented return value [AssemblyOpaqueTypeSet *]
 *
 */
AssemblyOpaqueTypeSet * Wise2_hard_link_AssemblyOpaqueTypeSet(AssemblyOpaqueTypeSet * obj);
#define hard_link_AssemblyOpaqueTypeSet Wise2_hard_link_AssemblyOpaqueTypeSet


/* Function:  AssemblyOpaqueTypeSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AssemblyOpaqueTypeSet *]
 *
 */
AssemblyOpaqueTypeSet * Wise2_AssemblyOpaqueTypeSet_alloc(void);
#define AssemblyOpaqueTypeSet_alloc Wise2_AssemblyOpaqueTypeSet_alloc


/* Function:  free_AssemblyOpaqueTypeSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AssemblyOpaqueTypeSet *]
 *
 * Return [UNKN ]  Undocumented return value [AssemblyOpaqueTypeSet *]
 *
 */
AssemblyOpaqueTypeSet * Wise2_free_AssemblyOpaqueTypeSet(AssemblyOpaqueTypeSet * obj);
#define free_AssemblyOpaqueTypeSet Wise2_free_AssemblyOpaqueTypeSet


/* Function:  hard_link_AssemblyOpaqueFeature(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AssemblyOpaqueFeature *]
 *
 * Return [UNKN ]  Undocumented return value [AssemblyOpaqueFeature *]
 *
 */
AssemblyOpaqueFeature * Wise2_hard_link_AssemblyOpaqueFeature(AssemblyOpaqueFeature * obj);
#define hard_link_AssemblyOpaqueFeature Wise2_hard_link_AssemblyOpaqueFeature


/* Function:  AssemblyOpaqueFeature_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AssemblyOpaqueFeature *]
 *
 */
AssemblyOpaqueFeature * Wise2_AssemblyOpaqueFeature_alloc(void);
#define AssemblyOpaqueFeature_alloc Wise2_AssemblyOpaqueFeature_alloc


/* Function:  free_AssemblyOpaqueFeature(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AssemblyOpaqueFeature *]
 *
 * Return [UNKN ]  Undocumented return value [AssemblyOpaqueFeature *]
 *
 */
AssemblyOpaqueFeature * Wise2_free_AssemblyOpaqueFeature(AssemblyOpaqueFeature * obj);
#define free_AssemblyOpaqueFeature Wise2_free_AssemblyOpaqueFeature


/* Function:  add_AssemblySequence(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [AssemblySequence *]
 * Arg:        add [OWNER] Object to add to the list [AssemblySequenceEdit *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_AssemblySequence(AssemblySequence * obj,AssemblySequenceEdit * add);
#define add_AssemblySequence Wise2_add_AssemblySequence


/* Function:  flush_AssemblySequence(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [AssemblySequence *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_AssemblySequence(AssemblySequence * obj);
#define flush_AssemblySequence Wise2_flush_AssemblySequence


/* Function:  add_opq_AssemblySequence(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [AssemblySequence *]
 * Arg:        add [OWNER] Object to add to the list [AssemblyOpaqueFeature *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_opq_AssemblySequence(AssemblySequence * obj,AssemblyOpaqueFeature * add);
#define add_opq_AssemblySequence Wise2_add_opq_AssemblySequence


/* Function:  flush_opq_AssemblySequence(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [AssemblySequence *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_opq_AssemblySequence(AssemblySequence * obj);
#define flush_opq_AssemblySequence Wise2_flush_opq_AssemblySequence


/* Function:  AssemblySequence_alloc_std(void)
 *
 * Descrip:    Equivalent to AssemblySequence_alloc_len(AssemblySequenceLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AssemblySequence *]
 *
 */
AssemblySequence * Wise2_AssemblySequence_alloc_std(void);
#define AssemblySequence_alloc_std Wise2_AssemblySequence_alloc_std


/* Function:  AssemblySequence_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [AssemblySequence *]
 *
 */
AssemblySequence * Wise2_AssemblySequence_alloc_len(int len);
#define AssemblySequence_alloc_len Wise2_AssemblySequence_alloc_len


/* Function:  hard_link_AssemblySequence(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AssemblySequence *]
 *
 * Return [UNKN ]  Undocumented return value [AssemblySequence *]
 *
 */
AssemblySequence * Wise2_hard_link_AssemblySequence(AssemblySequence * obj);
#define hard_link_AssemblySequence Wise2_hard_link_AssemblySequence


/* Function:  AssemblySequence_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AssemblySequence *]
 *
 */
AssemblySequence * Wise2_AssemblySequence_alloc(void);
#define AssemblySequence_alloc Wise2_AssemblySequence_alloc


/* Function:  free_AssemblySequence(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AssemblySequence *]
 *
 * Return [UNKN ]  Undocumented return value [AssemblySequence *]
 *
 */
AssemblySequence * Wise2_free_AssemblySequence(AssemblySequence * obj);
#define free_AssemblySequence Wise2_free_AssemblySequence


/* Function:  hard_link_AssemblySequencePlacement(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AssemblySequencePlacement *]
 *
 * Return [UNKN ]  Undocumented return value [AssemblySequencePlacement *]
 *
 */
AssemblySequencePlacement * Wise2_hard_link_AssemblySequencePlacement(AssemblySequencePlacement * obj);
#define hard_link_AssemblySequencePlacement Wise2_hard_link_AssemblySequencePlacement


/* Function:  AssemblySequencePlacement_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AssemblySequencePlacement *]
 *
 */
AssemblySequencePlacement * Wise2_AssemblySequencePlacement_alloc(void);
#define AssemblySequencePlacement_alloc Wise2_AssemblySequencePlacement_alloc


/* Function:  free_AssemblySequencePlacement(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AssemblySequencePlacement *]
 *
 * Return [UNKN ]  Undocumented return value [AssemblySequencePlacement *]
 *
 */
AssemblySequencePlacement * Wise2_free_AssemblySequencePlacement(AssemblySequencePlacement * obj);
#define free_AssemblySequencePlacement Wise2_free_AssemblySequencePlacement


/* Function:  add_AssemblyContig(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [AssemblyContig *]
 * Arg:        add [OWNER] Object to add to the list [AssemblySequencePlacement *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_AssemblyContig(AssemblyContig * obj,AssemblySequencePlacement * add);
#define add_AssemblyContig Wise2_add_AssemblyContig


/* Function:  flush_AssemblyContig(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [AssemblyContig *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_AssemblyContig(AssemblyContig * obj);
#define flush_AssemblyContig Wise2_flush_AssemblyContig


/* Function:  AssemblyContig_alloc_std(void)
 *
 * Descrip:    Equivalent to AssemblyContig_alloc_len(AssemblyContigLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AssemblyContig *]
 *
 */
AssemblyContig * Wise2_AssemblyContig_alloc_std(void);
#define AssemblyContig_alloc_std Wise2_AssemblyContig_alloc_std


/* Function:  AssemblyContig_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [AssemblyContig *]
 *
 */
AssemblyContig * Wise2_AssemblyContig_alloc_len(int len);
#define AssemblyContig_alloc_len Wise2_AssemblyContig_alloc_len


/* Function:  hard_link_AssemblyContig(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AssemblyContig *]
 *
 * Return [UNKN ]  Undocumented return value [AssemblyContig *]
 *
 */
AssemblyContig * Wise2_hard_link_AssemblyContig(AssemblyContig * obj);
#define hard_link_AssemblyContig Wise2_hard_link_AssemblyContig


/* Function:  AssemblyContig_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AssemblyContig *]
 *
 */
AssemblyContig * Wise2_AssemblyContig_alloc(void);
#define AssemblyContig_alloc Wise2_AssemblyContig_alloc


/* Function:  free_AssemblyContig(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AssemblyContig *]
 *
 * Return [UNKN ]  Undocumented return value [AssemblyContig *]
 *
 */
AssemblyContig * Wise2_free_AssemblyContig(AssemblyContig * obj);
#define free_AssemblyContig Wise2_free_AssemblyContig


/* Function:  add_Assembly(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Assembly *]
 * Arg:        add [OWNER] Object to add to the list [AssemblyContig *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_Assembly(Assembly * obj,AssemblyContig * add);
#define add_Assembly Wise2_add_Assembly


/* Function:  flush_Assembly(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [Assembly *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_Assembly(Assembly * obj);
#define flush_Assembly Wise2_flush_Assembly


/* Function:  add_chaff_Assembly(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Assembly *]
 * Arg:        add [OWNER] Object to add to the list [AssemblySequence *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_chaff_Assembly(Assembly * obj,AssemblySequence * add);
#define add_chaff_Assembly Wise2_add_chaff_Assembly


/* Function:  flush_chaff_Assembly(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [Assembly *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_chaff_Assembly(Assembly * obj);
#define flush_chaff_Assembly Wise2_flush_chaff_Assembly


/* Function:  Assembly_alloc_std(void)
 *
 * Descrip:    Equivalent to Assembly_alloc_len(AssemblyLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Assembly *]
 *
 */
Assembly * Wise2_Assembly_alloc_std(void);
#define Assembly_alloc_std Wise2_Assembly_alloc_std


/* Function:  Assembly_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [Assembly *]
 *
 */
Assembly * Wise2_Assembly_alloc_len(int len);
#define Assembly_alloc_len Wise2_Assembly_alloc_len


/* Function:  hard_link_Assembly(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Assembly *]
 *
 * Return [UNKN ]  Undocumented return value [Assembly *]
 *
 */
Assembly * Wise2_hard_link_Assembly(Assembly * obj);
#define hard_link_Assembly Wise2_hard_link_Assembly


/* Function:  Assembly_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Assembly *]
 *
 */
Assembly * Wise2_Assembly_alloc(void);
#define Assembly_alloc Wise2_Assembly_alloc


/* Function:  free_Assembly(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Assembly *]
 *
 * Return [UNKN ]  Undocumented return value [Assembly *]
 *
 */
Assembly * Wise2_free_Assembly(Assembly * obj);
#define free_Assembly Wise2_free_Assembly


/* Function:  hard_link_AssemblyOutputPara(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AssemblyOutputPara *]
 *
 * Return [UNKN ]  Undocumented return value [AssemblyOutputPara *]
 *
 */
AssemblyOutputPara * Wise2_hard_link_AssemblyOutputPara(AssemblyOutputPara * obj);
#define hard_link_AssemblyOutputPara Wise2_hard_link_AssemblyOutputPara


/* Function:  AssemblyOutputPara_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AssemblyOutputPara *]
 *
 */
AssemblyOutputPara * Wise2_AssemblyOutputPara_alloc(void);
#define AssemblyOutputPara_alloc Wise2_AssemblyOutputPara_alloc


/* Function:  free_AssemblyOutputPara(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AssemblyOutputPara *]
 *
 * Return [UNKN ]  Undocumented return value [AssemblyOutputPara *]
 *
 */
AssemblyOutputPara * Wise2_free_AssemblyOutputPara(AssemblyOutputPara * obj);
#define free_AssemblyOutputPara Wise2_free_AssemblyOutputPara


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
void Wise2_show_help_AssemblyOutputPara(FILE * ofp);
#define show_help_AssemblyOutputPara Wise2_show_help_AssemblyOutputPara
AssemblyOutputPara * Wise2_new_AssemblyOutputPara_from_argv(int * argc,char ** argv);
#define new_AssemblyOutputPara_from_argv Wise2_new_AssemblyOutputPara_from_argv
AssemblyOpaqueTypeSet * Wise2_homopolymer_AssemblyOpaqueTypeSet(void);
#define homopolymer_AssemblyOpaqueTypeSet Wise2_homopolymer_AssemblyOpaqueTypeSet
int Wise2_annotate_AssemblyOpaqueFeatures(AssemblyOpaqueTypeSet * aots,AssemblySequence * aseq,int kmer_size);
#define annotate_AssemblyOpaqueFeatures Wise2_annotate_AssemblyOpaqueFeatures
void Wise2_show_AssemblySequence(AssemblySequence * aseq,FILE * ofp);
#define show_AssemblySequence Wise2_show_AssemblySequence
AssemblyOpaqueType * Wise2_new_homopolymer_AssemblyOpaqueType(int int_type,char base);
#define new_homopolymer_AssemblyOpaqueType Wise2_new_homopolymer_AssemblyOpaqueType
AssemblySequence * Wise2_read_plain_fasta_AssemblySequence(FILE * ifp,int report_log,FILE * report);
#define read_plain_fasta_AssemblySequence Wise2_read_plain_fasta_AssemblySequence
AssemblySequence * Wise2_mirrored_AssemblySequence(AssemblySequence * aseq);
#define mirrored_AssemblySequence Wise2_mirrored_AssemblySequence
void Wise2_dump_contigs_as_fasta_Assembly(Assembly * assembly,AssemblyOutputPara * aop,FILE * ofp);
#define dump_contigs_as_fasta_Assembly Wise2_dump_contigs_as_fasta_Assembly


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_AssemblyOpaqueTypeSet(AssemblyOpaqueType ** list,int i,int j) ;
#define swap_AssemblyOpaqueTypeSet Wise2_swap_AssemblyOpaqueTypeSet
void Wise2_qsort_AssemblyOpaqueTypeSet(AssemblyOpaqueType ** list,int left,int right,int (*comp)(AssemblyOpaqueType * ,AssemblyOpaqueType * ));
#define qsort_AssemblyOpaqueTypeSet Wise2_qsort_AssemblyOpaqueTypeSet
void Wise2_sort_AssemblyOpaqueTypeSet(AssemblyOpaqueTypeSet * obj,int (*comp)(AssemblyOpaqueType *, AssemblyOpaqueType *));
#define sort_AssemblyOpaqueTypeSet Wise2_sort_AssemblyOpaqueTypeSet
boolean Wise2_expand_AssemblyOpaqueTypeSet(AssemblyOpaqueTypeSet * obj,int len);
#define expand_AssemblyOpaqueTypeSet Wise2_expand_AssemblyOpaqueTypeSet
void Wise2_swap_AssemblySequence(AssemblySequenceEdit ** list,int i,int j) ;
#define swap_AssemblySequence Wise2_swap_AssemblySequence
void Wise2_qsort_AssemblySequence(AssemblySequenceEdit ** list,int left,int right,int (*comp)(AssemblySequenceEdit * ,AssemblySequenceEdit * ));
#define qsort_AssemblySequence Wise2_qsort_AssemblySequence
void Wise2_sort_AssemblySequence(AssemblySequence * obj,int (*comp)(AssemblySequenceEdit *, AssemblySequenceEdit *));
#define sort_AssemblySequence Wise2_sort_AssemblySequence
boolean Wise2_expand_AssemblySequence(AssemblySequence * obj,int len);
#define expand_AssemblySequence Wise2_expand_AssemblySequence
void Wise2_swap_opq_AssemblySequence(AssemblyOpaqueFeature ** list,int i,int j) ;
#define swap_opq_AssemblySequence Wise2_swap_opq_AssemblySequence
void Wise2_qsort_opq_AssemblySequence(AssemblyOpaqueFeature ** list,int left,int right,int (*comp)(AssemblyOpaqueFeature * ,AssemblyOpaqueFeature * ));
#define qsort_opq_AssemblySequence Wise2_qsort_opq_AssemblySequence
void Wise2_sort_opq_AssemblySequence(AssemblySequence * obj,int (*comp)(AssemblyOpaqueFeature *, AssemblyOpaqueFeature *));
#define sort_opq_AssemblySequence Wise2_sort_opq_AssemblySequence
boolean Wise2_expand_opq_AssemblySequence(AssemblySequence * obj,int len);
#define expand_opq_AssemblySequence Wise2_expand_opq_AssemblySequence
void Wise2_swap_AssemblyContig(AssemblySequencePlacement ** list,int i,int j) ;
#define swap_AssemblyContig Wise2_swap_AssemblyContig
void Wise2_qsort_AssemblyContig(AssemblySequencePlacement ** list,int left,int right,int (*comp)(AssemblySequencePlacement * ,AssemblySequencePlacement * ));
#define qsort_AssemblyContig Wise2_qsort_AssemblyContig
void Wise2_sort_AssemblyContig(AssemblyContig * obj,int (*comp)(AssemblySequencePlacement *, AssemblySequencePlacement *));
#define sort_AssemblyContig Wise2_sort_AssemblyContig
boolean Wise2_expand_AssemblyContig(AssemblyContig * obj,int len);
#define expand_AssemblyContig Wise2_expand_AssemblyContig
void Wise2_swap_Assembly(AssemblyContig ** list,int i,int j) ;
#define swap_Assembly Wise2_swap_Assembly
void Wise2_qsort_Assembly(AssemblyContig ** list,int left,int right,int (*comp)(AssemblyContig * ,AssemblyContig * ));
#define qsort_Assembly Wise2_qsort_Assembly
void Wise2_sort_Assembly(Assembly * obj,int (*comp)(AssemblyContig *, AssemblyContig *));
#define sort_Assembly Wise2_sort_Assembly
boolean Wise2_expand_Assembly(Assembly * obj,int len);
#define expand_Assembly Wise2_expand_Assembly
void Wise2_swap_chaff_Assembly(AssemblySequence ** list,int i,int j) ;
#define swap_chaff_Assembly Wise2_swap_chaff_Assembly
void Wise2_qsort_chaff_Assembly(AssemblySequence ** list,int left,int right,int (*comp)(AssemblySequence * ,AssemblySequence * ));
#define qsort_chaff_Assembly Wise2_qsort_chaff_Assembly
void Wise2_sort_chaff_Assembly(Assembly * obj,int (*comp)(AssemblySequence *, AssemblySequence *));
#define sort_chaff_Assembly Wise2_sort_chaff_Assembly
boolean Wise2_expand_chaff_Assembly(Assembly * obj,int len);
#define expand_chaff_Assembly Wise2_expand_chaff_Assembly

#ifdef _cplusplus
}
#endif

#endif
