#ifndef DYNAMITElinkindexHEADERFILE
#define DYNAMITElinkindexHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dnanumber.h"


#define LinkNumberLISTLENGTH 32
#define LinkStreamLISTLENGTH 32

#define LinkNumberArrayLISTLENGTH 1024

#ifndef DYNAMITE_DEFINED_LinkNumber
typedef struct Wise2_LinkNumber Wise2_LinkNumber;
#define LinkNumber Wise2_LinkNumber
#define DYNAMITE_DEFINED_LinkNumber
#endif

#ifndef DYNAMITE_DEFINED_LinkStream
typedef struct Wise2_LinkStream Wise2_LinkStream;
#define LinkStream Wise2_LinkStream
#define DYNAMITE_DEFINED_LinkStream
#endif

#ifndef DYNAMITE_DEFINED_LinkEdge
typedef struct Wise2_LinkEdge Wise2_LinkEdge;
#define LinkEdge Wise2_LinkEdge
#define DYNAMITE_DEFINED_LinkEdge
#endif

struct Wise2_LinkEdge {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    LinkStream * x;  
    LinkStream * y;  
    char twist;  
    } ;  
/* LinkEdge defined */ 
#ifndef DYNAMITE_DEFINED_LinkEdge
typedef struct Wise2_LinkEdge Wise2_LinkEdge;
#define LinkEdge Wise2_LinkEdge
#define DYNAMITE_DEFINED_LinkEdge
#endif


struct Wise2_LinkStream {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    LinkEdge * a;    
    LinkEdge * b;    
    int number;  
    char depth;  
    char starting_flip;  
    char have_seen;  
    } ;  
/* LinkStream defined */ 
#ifndef DYNAMITE_DEFINED_LinkStream
typedef struct Wise2_LinkStream Wise2_LinkStream;
#define LinkStream Wise2_LinkStream
#define DYNAMITE_DEFINED_LinkStream
#endif


struct Wise2_LinkNumber {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    LinkStream ** stream;    
    int len;/* len for above stream  */ 
    int maxlen; /* maxlen for above stream */ 
    int number;  
    } ;  
/* LinkNumber defined */ 
#ifndef DYNAMITE_DEFINED_LinkNumber
typedef struct Wise2_LinkNumber Wise2_LinkNumber;
#define LinkNumber Wise2_LinkNumber
#define DYNAMITE_DEFINED_LinkNumber
#endif


struct Wise2_LinkNumberArray {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    LinkNumber ** array;     
    int array_len;   
    int nmer_size;   
    LinkEdge ** edge_set;    
    int len;/* len for above edge_set  */ 
    int maxlen; /* maxlen for above edge_set */ 
    } ;  
/* LinkNumberArray defined */ 
#ifndef DYNAMITE_DEFINED_LinkNumberArray
typedef struct Wise2_LinkNumberArray Wise2_LinkNumberArray;
#define LinkNumberArray Wise2_LinkNumberArray
#define DYNAMITE_DEFINED_LinkNumberArray
#endif


struct Wise2_LinkNumberArrayDebug {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char placement_stream;   
    char add_stream;     
    char extraction;     
    FILE * ofp;  
    } ;  
/* LinkNumberArrayDebug defined */ 
#ifndef DYNAMITE_DEFINED_LinkNumberArrayDebug
typedef struct Wise2_LinkNumberArrayDebug Wise2_LinkNumberArrayDebug;
#define LinkNumberArrayDebug Wise2_LinkNumberArrayDebug
#define DYNAMITE_DEFINED_LinkNumberArrayDebug
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  extract_dna_LinkStream(ln,lnad,nmer_size)
 *
 * Descrip:    Reverse builds a DNA sequence stream from LinkStream
 *
 *
 * Arg:               ln [UNKN ] Undocumented argument [LinkStream *]
 * Arg:             lnad [UNKN ] Undocumented argument [LinkNumberArrayDebug *]
 * Arg:        nmer_size [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [Sequence *]
 *
 */
Sequence * Wise2_extract_dna_LinkStream(LinkStream * ln,LinkNumberArrayDebug * lnad,int nmer_size);
#define extract_dna_LinkStream Wise2_extract_dna_LinkStream


/* Function:  write_existing_LinkStream(lna,lnad,dns)
 *
 * Descrip:    Writes in read to existing streams
 *
 *
 * Arg:         lna [UNKN ] Undocumented argument [LinkNumberArray *]
 * Arg:        lnad [UNKN ] Undocumented argument [LinkNumberArrayDebug *]
 * Arg:         dns [UNKN ] Undocumented argument [DnaNumberSequence *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_write_existing_LinkStream(LinkNumberArray * lna,LinkNumberArrayDebug * lnad,DnaNumberSequence * dns);
#define write_existing_LinkStream Wise2_write_existing_LinkStream


/* Function:  write_new_LinkStream(lna,lnad,dns)
 *
 * Descrip:    Writes in read as a new stream
 *
 *
 * Arg:         lna [UNKN ] Undocumented argument [LinkNumberArray *]
 * Arg:        lnad [UNKN ] Undocumented argument [LinkNumberArrayDebug *]
 * Arg:         dns [UNKN ] Undocumented argument [DnaNumberSequence *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_write_new_LinkStream(LinkNumberArray * lna,LinkNumberArrayDebug * lnad,DnaNumberSequence * dns);
#define write_new_LinkStream Wise2_write_new_LinkStream


/* Function:  is_linked_LinkStream(l_one,l_two,twist)
 *
 * Descrip:    Returns the right edge which potentially links these two link streams
 *
 *
 * Arg:        l_one [UNKN ] Undocumented argument [LinkStream *]
 * Arg:        l_two [UNKN ] Undocumented argument [LinkStream *]
 * Arg:        twist [UNKN ] Undocumented argument [char]
 *
 * Return [UNKN ]  Undocumented return value [LinkEdge *]
 *
 */
LinkEdge * Wise2_is_linked_LinkStream(LinkStream * l_one,LinkStream * l_two,char twist);
#define is_linked_LinkStream Wise2_is_linked_LinkStream


/* Function:  is_new_stream_LinkNumberArray(lna,lnad,dns)
 *
 * Descrip:    Determines whether this DnaNumberSequence should be newly
 *             streamed or not
 *
 *
 * Arg:         lna [UNKN ] Undocumented argument [LinkNumberArray *]
 * Arg:        lnad [UNKN ] Undocumented argument [LinkNumberArrayDebug *]
 * Arg:         dns [UNKN ] Undocumented argument [DnaNumberSequence *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_is_new_stream_LinkNumberArray(LinkNumberArray * lna,LinkNumberArrayDebug * lnad,DnaNumberSequence * dns);
#define is_new_stream_LinkNumberArray Wise2_is_new_stream_LinkNumberArray


/* Function:  new_LinkNumberArray(nmer_size)
 *
 * Descrip:    Returns a LinkNumberArray with the appropiate
 *             array size for this size of Nmer
 *
 *
 * Arg:        nmer_size [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [LinkNumberArray *]
 *
 */
LinkNumberArray * Wise2_new_LinkNumberArray(int nmer_size);
#define new_LinkNumberArray Wise2_new_LinkNumberArray


/* Function:  new_LinkEdge(lna)
 *
 * Descrip:    makes a new link edge added into the memory management
 *
 *
 * Arg:        lna [UNKN ] Undocumented argument [LinkNumberArray *]
 *
 * Return [UNKN ]  Undocumented return value [LinkEdge *]
 *
 */
LinkEdge * Wise2_new_LinkEdge(LinkNumberArray * lna);
#define new_LinkEdge Wise2_new_LinkEdge


/* Function:  hard_link_LinkEdge(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [LinkEdge *]
 *
 * Return [UNKN ]  Undocumented return value [LinkEdge *]
 *
 */
LinkEdge * Wise2_hard_link_LinkEdge(LinkEdge * obj);
#define hard_link_LinkEdge Wise2_hard_link_LinkEdge


/* Function:  LinkEdge_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LinkEdge *]
 *
 */
LinkEdge * Wise2_LinkEdge_alloc(void);
#define LinkEdge_alloc Wise2_LinkEdge_alloc


/* Function:  free_LinkEdge(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [LinkEdge *]
 *
 * Return [UNKN ]  Undocumented return value [LinkEdge *]
 *
 */
LinkEdge * Wise2_free_LinkEdge(LinkEdge * obj);
#define free_LinkEdge Wise2_free_LinkEdge


/* Function:  hard_link_LinkStream(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [LinkStream *]
 *
 * Return [UNKN ]  Undocumented return value [LinkStream *]
 *
 */
LinkStream * Wise2_hard_link_LinkStream(LinkStream * obj);
#define hard_link_LinkStream Wise2_hard_link_LinkStream


/* Function:  LinkStream_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LinkStream *]
 *
 */
LinkStream * Wise2_LinkStream_alloc(void);
#define LinkStream_alloc Wise2_LinkStream_alloc


/* Function:  free_LinkStream(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [LinkStream *]
 *
 * Return [UNKN ]  Undocumented return value [LinkStream *]
 *
 */
LinkStream * Wise2_free_LinkStream(LinkStream * obj);
#define free_LinkStream Wise2_free_LinkStream


/* Function:  add_LinkNumber(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [LinkNumber *]
 * Arg:        add [OWNER] Object to add to the list [LinkStream *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_LinkNumber(LinkNumber * obj,LinkStream * add);
#define add_LinkNumber Wise2_add_LinkNumber


/* Function:  flush_LinkNumber(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [LinkNumber *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_LinkNumber(LinkNumber * obj);
#define flush_LinkNumber Wise2_flush_LinkNumber


/* Function:  LinkNumber_alloc_std(void)
 *
 * Descrip:    Equivalent to LinkNumber_alloc_len(LinkNumberLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LinkNumber *]
 *
 */
LinkNumber * Wise2_LinkNumber_alloc_std(void);
#define LinkNumber_alloc_std Wise2_LinkNumber_alloc_std


/* Function:  LinkNumber_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [LinkNumber *]
 *
 */
LinkNumber * Wise2_LinkNumber_alloc_len(int len);
#define LinkNumber_alloc_len Wise2_LinkNumber_alloc_len


/* Function:  hard_link_LinkNumber(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [LinkNumber *]
 *
 * Return [UNKN ]  Undocumented return value [LinkNumber *]
 *
 */
LinkNumber * Wise2_hard_link_LinkNumber(LinkNumber * obj);
#define hard_link_LinkNumber Wise2_hard_link_LinkNumber


/* Function:  LinkNumber_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LinkNumber *]
 *
 */
LinkNumber * Wise2_LinkNumber_alloc(void);
#define LinkNumber_alloc Wise2_LinkNumber_alloc


/* Function:  free_LinkNumber(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [LinkNumber *]
 *
 * Return [UNKN ]  Undocumented return value [LinkNumber *]
 *
 */
LinkNumber * Wise2_free_LinkNumber(LinkNumber * obj);
#define free_LinkNumber Wise2_free_LinkNumber


/* Function:  add_LinkNumberArray(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [LinkNumberArray *]
 * Arg:        add [OWNER] Object to add to the list [LinkEdge *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_LinkNumberArray(LinkNumberArray * obj,LinkEdge * add);
#define add_LinkNumberArray Wise2_add_LinkNumberArray


/* Function:  flush_LinkNumberArray(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [LinkNumberArray *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_LinkNumberArray(LinkNumberArray * obj);
#define flush_LinkNumberArray Wise2_flush_LinkNumberArray


/* Function:  LinkNumberArray_alloc_std(void)
 *
 * Descrip:    Equivalent to LinkNumberArray_alloc_len(LinkNumberArrayLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LinkNumberArray *]
 *
 */
LinkNumberArray * Wise2_LinkNumberArray_alloc_std(void);
#define LinkNumberArray_alloc_std Wise2_LinkNumberArray_alloc_std


/* Function:  LinkNumberArray_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [LinkNumberArray *]
 *
 */
LinkNumberArray * Wise2_LinkNumberArray_alloc_len(int len);
#define LinkNumberArray_alloc_len Wise2_LinkNumberArray_alloc_len


/* Function:  hard_link_LinkNumberArray(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [LinkNumberArray *]
 *
 * Return [UNKN ]  Undocumented return value [LinkNumberArray *]
 *
 */
LinkNumberArray * Wise2_hard_link_LinkNumberArray(LinkNumberArray * obj);
#define hard_link_LinkNumberArray Wise2_hard_link_LinkNumberArray


/* Function:  LinkNumberArray_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LinkNumberArray *]
 *
 */
LinkNumberArray * Wise2_LinkNumberArray_alloc(void);
#define LinkNumberArray_alloc Wise2_LinkNumberArray_alloc


/* Function:  free_LinkNumberArray(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [LinkNumberArray *]
 *
 * Return [UNKN ]  Undocumented return value [LinkNumberArray *]
 *
 */
LinkNumberArray * Wise2_free_LinkNumberArray(LinkNumberArray * obj);
#define free_LinkNumberArray Wise2_free_LinkNumberArray


/* Function:  hard_link_LinkNumberArrayDebug(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [LinkNumberArrayDebug *]
 *
 * Return [UNKN ]  Undocumented return value [LinkNumberArrayDebug *]
 *
 */
LinkNumberArrayDebug * Wise2_hard_link_LinkNumberArrayDebug(LinkNumberArrayDebug * obj);
#define hard_link_LinkNumberArrayDebug Wise2_hard_link_LinkNumberArrayDebug


/* Function:  LinkNumberArrayDebug_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LinkNumberArrayDebug *]
 *
 */
LinkNumberArrayDebug * Wise2_LinkNumberArrayDebug_alloc(void);
#define LinkNumberArrayDebug_alloc Wise2_LinkNumberArrayDebug_alloc


/* Function:  free_LinkNumberArrayDebug(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [LinkNumberArrayDebug *]
 *
 * Return [UNKN ]  Undocumented return value [LinkNumberArrayDebug *]
 *
 */
LinkNumberArrayDebug * Wise2_free_LinkNumberArrayDebug(LinkNumberArrayDebug * obj);
#define free_LinkNumberArrayDebug Wise2_free_LinkNumberArrayDebug


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_LinkNumber(LinkStream ** list,int i,int j) ;
#define swap_LinkNumber Wise2_swap_LinkNumber
void Wise2_qsort_LinkNumber(LinkStream ** list,int left,int right,int (*comp)(LinkStream * ,LinkStream * ));
#define qsort_LinkNumber Wise2_qsort_LinkNumber
void Wise2_sort_LinkNumber(LinkNumber * obj,int (*comp)(LinkStream *, LinkStream *));
#define sort_LinkNumber Wise2_sort_LinkNumber
boolean Wise2_expand_LinkNumber(LinkNumber * obj,int len);
#define expand_LinkNumber Wise2_expand_LinkNumber
void Wise2_swap_LinkNumberArray(LinkEdge ** list,int i,int j) ;
#define swap_LinkNumberArray Wise2_swap_LinkNumberArray
void Wise2_qsort_LinkNumberArray(LinkEdge ** list,int left,int right,int (*comp)(LinkEdge * ,LinkEdge * ));
#define qsort_LinkNumberArray Wise2_qsort_LinkNumberArray
void Wise2_sort_LinkNumberArray(LinkNumberArray * obj,int (*comp)(LinkEdge *, LinkEdge *));
#define sort_LinkNumberArray Wise2_sort_LinkNumberArray
boolean Wise2_expand_LinkNumberArray(LinkNumberArray * obj,int len);
#define expand_LinkNumberArray Wise2_expand_LinkNumberArray

#ifdef _cplusplus
}
#endif

#endif
