#ifndef DYNAMITEtransferinterfaceHEADERFILE
#define DYNAMITEtransferinterfaceHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "wisebase.h"

#define TransferedFunctionCallLISTLENGTH 16

struct Wise2_TransferedObjectMarshaller {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    void (*object_to_stream)(void *,Wise2WriteStreamInterface*); 
    void * (*stream_to_object)(Wise2ReadStreamInterface*);   
    char * object_type;  
    void (*untyped_free_object)(void *); 
    } ;  
/* TransferedObjectMarshaller defined */ 
#ifndef DYNAMITE_DEFINED_TransferedObjectMarshaller
typedef struct Wise2_TransferedObjectMarshaller Wise2_TransferedObjectMarshaller;
#define TransferedObjectMarshaller Wise2_TransferedObjectMarshaller
#define DYNAMITE_DEFINED_TransferedObjectMarshaller
#endif


struct Wise2_TransferedFunctionCall {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * name;     
    TransferedObjectMarshaller * returned_type;  
    TransferedObjectMarshaller ** input;     
    int len;/* len for above input  */ 
    int maxlen; /* maxlen for above input */ 
    } ;  
/* TransferedFunctionCall defined */ 
#ifndef DYNAMITE_DEFINED_TransferedFunctionCall
typedef struct Wise2_TransferedFunctionCall Wise2_TransferedFunctionCall;
#define TransferedFunctionCall Wise2_TransferedFunctionCall
#define DYNAMITE_DEFINED_TransferedFunctionCall
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_TransferedObjectMarshaller(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [TransferedObjectMarshaller *]
 *
 * Return [UNKN ]  Undocumented return value [TransferedObjectMarshaller *]
 *
 */
TransferedObjectMarshaller * Wise2_hard_link_TransferedObjectMarshaller(TransferedObjectMarshaller * obj);
#define hard_link_TransferedObjectMarshaller Wise2_hard_link_TransferedObjectMarshaller


/* Function:  TransferedObjectMarshaller_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TransferedObjectMarshaller *]
 *
 */
TransferedObjectMarshaller * Wise2_TransferedObjectMarshaller_alloc(void);
#define TransferedObjectMarshaller_alloc Wise2_TransferedObjectMarshaller_alloc


/* Function:  free_TransferedObjectMarshaller(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [TransferedObjectMarshaller *]
 *
 * Return [UNKN ]  Undocumented return value [TransferedObjectMarshaller *]
 *
 */
TransferedObjectMarshaller * Wise2_free_TransferedObjectMarshaller(TransferedObjectMarshaller * obj);
#define free_TransferedObjectMarshaller Wise2_free_TransferedObjectMarshaller


/* Function:  add_TransferedFunctionCall(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [TransferedFunctionCall *]
 * Arg:        add [OWNER] Object to add to the list [TransferedObjectMarshaller *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_TransferedFunctionCall(TransferedFunctionCall * obj,TransferedObjectMarshaller * add);
#define add_TransferedFunctionCall Wise2_add_TransferedFunctionCall


/* Function:  flush_TransferedFunctionCall(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [TransferedFunctionCall *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_TransferedFunctionCall(TransferedFunctionCall * obj);
#define flush_TransferedFunctionCall Wise2_flush_TransferedFunctionCall


/* Function:  TransferedFunctionCall_alloc_std(void)
 *
 * Descrip:    Equivalent to TransferedFunctionCall_alloc_len(TransferedFunctionCallLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TransferedFunctionCall *]
 *
 */
TransferedFunctionCall * Wise2_TransferedFunctionCall_alloc_std(void);
#define TransferedFunctionCall_alloc_std Wise2_TransferedFunctionCall_alloc_std


/* Function:  TransferedFunctionCall_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [TransferedFunctionCall *]
 *
 */
TransferedFunctionCall * Wise2_TransferedFunctionCall_alloc_len(int len);
#define TransferedFunctionCall_alloc_len Wise2_TransferedFunctionCall_alloc_len


/* Function:  hard_link_TransferedFunctionCall(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [TransferedFunctionCall *]
 *
 * Return [UNKN ]  Undocumented return value [TransferedFunctionCall *]
 *
 */
TransferedFunctionCall * Wise2_hard_link_TransferedFunctionCall(TransferedFunctionCall * obj);
#define hard_link_TransferedFunctionCall Wise2_hard_link_TransferedFunctionCall


/* Function:  TransferedFunctionCall_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TransferedFunctionCall *]
 *
 */
TransferedFunctionCall * Wise2_TransferedFunctionCall_alloc(void);
#define TransferedFunctionCall_alloc Wise2_TransferedFunctionCall_alloc


/* Function:  free_TransferedFunctionCall(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [TransferedFunctionCall *]
 *
 * Return [UNKN ]  Undocumented return value [TransferedFunctionCall *]
 *
 */
TransferedFunctionCall * Wise2_free_TransferedFunctionCall(TransferedFunctionCall * obj);
#define free_TransferedFunctionCall Wise2_free_TransferedFunctionCall


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
TransferedFunctionCall * Wise2_test_stringcat_TransferedFunctionCall(void);
#define test_stringcat_TransferedFunctionCall Wise2_test_stringcat_TransferedFunctionCall
TransferedObjectMarshaller * Wise2_new_string_Marshaller(void);
#define new_string_Marshaller Wise2_new_string_Marshaller
void Wise2_object_to_stream_string(void * obj,Wise2WriteStreamInterface* write);
#define object_to_stream_string Wise2_object_to_stream_string
void * Wise2_stream_to_object_string(Wise2ReadStreamInterface * read);
#define stream_to_object_string Wise2_stream_to_object_string


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_TransferedFunctionCall(TransferedObjectMarshaller ** list,int i,int j) ;
#define swap_TransferedFunctionCall Wise2_swap_TransferedFunctionCall
void Wise2_qsort_TransferedFunctionCall(TransferedObjectMarshaller ** list,int left,int right,int (*comp)(TransferedObjectMarshaller * ,TransferedObjectMarshaller * ));
#define qsort_TransferedFunctionCall Wise2_qsort_TransferedFunctionCall
void Wise2_sort_TransferedFunctionCall(TransferedFunctionCall * obj,int (*comp)(TransferedObjectMarshaller *, TransferedObjectMarshaller *));
#define sort_TransferedFunctionCall Wise2_sort_TransferedFunctionCall
boolean Wise2_expand_TransferedFunctionCall(TransferedFunctionCall * obj,int len);
#define expand_TransferedFunctionCall Wise2_expand_TransferedFunctionCall

#ifdef _cplusplus
}
#endif

#endif
