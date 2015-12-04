#ifndef DYNAMITElabelmasterHEADERFILE
#define DYNAMITElabelmasterHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif

#include "wisebase.h"
#include "dyna2.h"


#define LabelLISTLENGTH 64
#define LabelMasterLISTLENGTH 64

#define LMNOTYPE    213
#define LMQUERYTYPE 214
#define LMTARGETTYPE 215

struct LabelInstance {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * state;    
    char * source;   
    int offi;    
    int offj;    
    char * calc_line;    
    } ;  
/* LabelInstance defined */ 
#ifndef DYNAMITE_DEFINED_LabelInstance
typedef struct LabelInstance LabelInstance;
#define DYNAMITE_DEFINED_LabelInstance
#endif


struct Label {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * name;     
    int type;    
    LabelInstance ** li;     
    int len;/* len for above li  */ 
    int maxlen; /* maxlen for above li */ 
    } ;  
/* Label defined */ 
#ifndef DYNAMITE_DEFINED_Label
typedef struct Label Label;
#define DYNAMITE_DEFINED_Label
#endif


struct LabelMaster {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Label ** query;  
    int query_len;  /* len for above query  */ 
    int query_maxlen;   /* maxlen for above query */ 
    Label ** target;     
    int target_len; /* len for above target  */ 
    int target_maxlen;  /* maxlen for above target */ 
    } ;  
/* LabelMaster defined */ 
#ifndef DYNAMITE_DEFINED_LabelMaster
typedef struct LabelMaster LabelMaster;
#define DYNAMITE_DEFINED_LabelMaster
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_LabelInstance(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [LabelInstance *]
 *
 * Return [UNKN ]  Undocumented return value [LabelInstance *]
 *
 */
LabelInstance * hard_link_LabelInstance(LabelInstance * obj);


/* Function:  LabelInstance_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LabelInstance *]
 *
 */
LabelInstance * LabelInstance_alloc(void);


/* Function:  free_LabelInstance(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [LabelInstance *]
 *
 * Return [UNKN ]  Undocumented return value [LabelInstance *]
 *
 */
LabelInstance * free_LabelInstance(LabelInstance * obj);


/* Function:  add_Label(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Label *]
 * Arg:        add [OWNER] Object to add to the list [LabelInstance *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean add_Label(Label * obj,LabelInstance * add);


/* Function:  flush_Label(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [Label *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_Label(Label * obj);


/* Function:  Label_alloc_std(void)
 *
 * Descrip:    Equivalent to Label_alloc_len(LabelLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Label *]
 *
 */
Label * Label_alloc_std(void);


/* Function:  Label_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [Label *]
 *
 */
Label * Label_alloc_len(int len);


/* Function:  hard_link_Label(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Label *]
 *
 * Return [UNKN ]  Undocumented return value [Label *]
 *
 */
Label * hard_link_Label(Label * obj);


/* Function:  Label_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Label *]
 *
 */
Label * Label_alloc(void);


/* Function:  free_Label(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Label *]
 *
 * Return [UNKN ]  Undocumented return value [Label *]
 *
 */
Label * free_Label(Label * obj);


/* Function:  add_query_LabelMaster(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [LabelMaster *]
 * Arg:        add [OWNER] Object to add to the list [Label *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean add_query_LabelMaster(LabelMaster * obj,Label * add);


/* Function:  flush_query_LabelMaster(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [LabelMaster *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_query_LabelMaster(LabelMaster * obj);


/* Function:  add_target_LabelMaster(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [LabelMaster *]
 * Arg:        add [OWNER] Object to add to the list [Label *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean add_target_LabelMaster(LabelMaster * obj,Label * add);


/* Function:  flush_target_LabelMaster(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [LabelMaster *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_target_LabelMaster(LabelMaster * obj);


/* Function:  LabelMaster_alloc_std(void)
 *
 * Descrip:    Equivalent to LabelMaster_alloc_len(LabelMasterLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LabelMaster *]
 *
 */
LabelMaster * LabelMaster_alloc_std(void);


/* Function:  LabelMaster_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [LabelMaster *]
 *
 */
LabelMaster * LabelMaster_alloc_len(int len);


/* Function:  hard_link_LabelMaster(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [LabelMaster *]
 *
 * Return [UNKN ]  Undocumented return value [LabelMaster *]
 *
 */
LabelMaster * hard_link_LabelMaster(LabelMaster * obj);


/* Function:  LabelMaster_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LabelMaster *]
 *
 */
LabelMaster * LabelMaster_alloc(void);


/* Function:  free_LabelMaster(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [LabelMaster *]
 *
 * Return [UNKN ]  Undocumented return value [LabelMaster *]
 *
 */
LabelMaster * free_LabelMaster(LabelMaster * obj);


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
char * quoted_string_from_Label(Label ** list,int len);
char * query_quoted_string_from_LabelMaster(LabelMaster * lm);
char * target_quoted_string_from_LabelMaster(LabelMaster * lm);
int index_for_label(char * name,Label ** list,int len);
int index_for_query_label(char * name,LabelMaster * lm);
int index_for_target_label(char * name,LabelMaster * lm);
Label * target_Label_from_name(LabelMaster * lm,char * name);
Label * query_Label_from_name(LabelMaster * lm,char * name);
LabelMaster * LabelMaster_from_GenericMatrix(GenericMatrix * gm);
boolean add_CellSource_to_LabelMaster(LabelMaster * lm,CellSource * cs,char * state);
Label * new_query_Label(char * name);
Label * new_target_Label(char * name);
LabelInstance * LabelInstance_from_CellSource(char * state,CellSource * cs);


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void swap_Label(LabelInstance ** list,int i,int j) ;
void qsort_Label(LabelInstance ** list,int left,int right,int (*comp)(LabelInstance * ,LabelInstance * ));
void sort_Label(Label * obj,int (*comp)(LabelInstance *, LabelInstance *));
boolean expand_Label(Label * obj,int len);
void swap_query_LabelMaster(Label ** list,int i,int j) ;
void qsort_query_LabelMaster(Label ** list,int left,int right,int (*comp)(Label * ,Label * ));
void sort_query_LabelMaster(LabelMaster * obj,int (*comp)(Label *, Label *));
boolean expand_query_LabelMaster(LabelMaster * obj,int len);
void swap_target_LabelMaster(Label ** list,int i,int j) ;
void qsort_target_LabelMaster(Label ** list,int left,int right,int (*comp)(Label * ,Label * ));
void sort_target_LabelMaster(LabelMaster * obj,int (*comp)(Label *, Label *));
boolean expand_target_LabelMaster(LabelMaster * obj,int len);

#ifdef _cplusplus
}
#endif

#endif
