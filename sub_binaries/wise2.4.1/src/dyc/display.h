#ifndef DYNAMITEdisplayHEADERFILE
#define DYNAMITEdisplayHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif

#include "wisebase.h"
#include "dynfile.h"


/********************************************************/
/* Wisetools version 3 file                             */
/*                                                      */
/* this file is copyright (c) Ewan Birney 1996          */
/* this file is part of the wisetools sequence analysis */
/* package                                              */
/********************************************************/

/********************************************************/
/* Display is an internal module for dynamite itself    */
/* to manage the display tag found in dynamite files    */
/*                                                      */
/*                                                      */
/********************************************************/

/***** RCS Info *****************************************/
/*
   $Id: 

   $Log: 
*/
/*********************************************************/



#define Label2DisplayLISTLENGTH    32
#define Sequence2DisplayLISTLENGTH 32
#define Aln2DisplayLISTLENGTH      32

#define ALL_OTHER_NUMBER (-10)

struct Aln2DisplayField {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * string;   
    char * length;   
    char * direction;    
    char * convert;  
    boolean is_static;   
    boolean is_string;   
    boolean is_single;   
    boolean is_sub;  
    } ;  
/* Aln2DisplayField defined */ 
#ifndef DYNAMITE_DEFINED_Aln2DisplayField
typedef struct Aln2DisplayField Aln2DisplayField;
#define DYNAMITE_DEFINED_Aln2DisplayField
#endif


struct Index2Display {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * name;     
    char * eval;     
    boolean is_start;    
    int number;  
    } ;  
/* Index2Display defined */ 
#ifndef DYNAMITE_DEFINED_Index2Display
typedef struct Index2Display Index2Display;
#define DYNAMITE_DEFINED_Index2Display
#endif


struct Label2Display {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * label;    
    char * direction;    
    Aln2DisplayField ** a2df;    
    int len;/* len for above a2df  */ 
    int maxlen; /* maxlen for above a2df */ 
    Index2Display ** i2d;    
    int ind_len;/* len for above i2d  */ 
    int ind_maxlen; /* maxlen for above i2d */ 
    boolean is_lone;     
    } ;  
/* Label2Display defined */ 
#ifndef DYNAMITE_DEFINED_Label2Display
typedef struct Label2Display Label2Display;
#define DYNAMITE_DEFINED_Label2Display
#endif


struct Sequence2Display {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int number;  
    Label2Display ** l2d;    
    int len;/* len for above l2d  */ 
    int maxlen; /* maxlen for above l2d */ 
    Index2Display ** i2d;    
    int ind_len;/* len for above i2d  */ 
    int ind_maxlen; /* maxlen for above i2d */ 
    char * name_str;     
    int length; /*  of name */ 
    boolean is_static;   
    char * direction;    
    } ;  
/* Sequence2Display defined */ 
#ifndef DYNAMITE_DEFINED_Sequence2Display
typedef struct Sequence2Display Sequence2Display;
#define DYNAMITE_DEFINED_Sequence2Display
#endif


struct Aln2DisplayResource {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * type;     
    char * name;     
    char * arg;  
    } ;  
/* Aln2DisplayResource defined */ 
#ifndef DYNAMITE_DEFINED_Aln2DisplayResource
typedef struct Aln2DisplayResource Aln2DisplayResource;
#define DYNAMITE_DEFINED_Aln2DisplayResource
#endif


struct Aln2Display {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * name;     
    Aln2DisplayResource ** a2dr;     
    int res_len;/* len for above a2dr  */ 
    int res_maxlen; /* maxlen for above a2dr */ 
    Sequence2Display ** s2d;     
    int len;/* len for above s2d  */ 
    int maxlen; /* maxlen for above s2d */ 
    Sequence2Display * other;    
    } ;  
/* Aln2Display defined */ 
#ifndef DYNAMITE_DEFINED_Aln2Display
typedef struct Aln2Display Aln2Display;
#define DYNAMITE_DEFINED_Aln2Display
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_Aln2DisplayField(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Aln2DisplayField *]
 *
 * Return [UNKN ]  Undocumented return value [Aln2DisplayField *]
 *
 */
Aln2DisplayField * hard_link_Aln2DisplayField(Aln2DisplayField * obj);


/* Function:  Aln2DisplayField_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Aln2DisplayField *]
 *
 */
Aln2DisplayField * Aln2DisplayField_alloc(void);


/* Function:  free_Aln2DisplayField(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Aln2DisplayField *]
 *
 * Return [UNKN ]  Undocumented return value [Aln2DisplayField *]
 *
 */
Aln2DisplayField * free_Aln2DisplayField(Aln2DisplayField * obj);


/* Function:  hard_link_Index2Display(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Index2Display *]
 *
 * Return [UNKN ]  Undocumented return value [Index2Display *]
 *
 */
Index2Display * hard_link_Index2Display(Index2Display * obj);


/* Function:  Index2Display_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Index2Display *]
 *
 */
Index2Display * Index2Display_alloc(void);


/* Function:  free_Index2Display(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Index2Display *]
 *
 * Return [UNKN ]  Undocumented return value [Index2Display *]
 *
 */
Index2Display * free_Index2Display(Index2Display * obj);


/* Function:  add_Label2Display(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Label2Display *]
 * Arg:        add [OWNER] Object to add to the list [Aln2DisplayField *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean add_Label2Display(Label2Display * obj,Aln2DisplayField * add);


/* Function:  flush_Label2Display(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [Label2Display *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_Label2Display(Label2Display * obj);


/* Function:  add_ind_Label2Display(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Label2Display *]
 * Arg:        add [OWNER] Object to add to the list [Index2Display *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean add_ind_Label2Display(Label2Display * obj,Index2Display * add);


/* Function:  flush_ind_Label2Display(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [Label2Display *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_ind_Label2Display(Label2Display * obj);


/* Function:  Label2Display_alloc_std(void)
 *
 * Descrip:    Equivalent to Label2Display_alloc_len(Label2DisplayLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Label2Display *]
 *
 */
Label2Display * Label2Display_alloc_std(void);


/* Function:  Label2Display_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [Label2Display *]
 *
 */
Label2Display * Label2Display_alloc_len(int len);


/* Function:  hard_link_Label2Display(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Label2Display *]
 *
 * Return [UNKN ]  Undocumented return value [Label2Display *]
 *
 */
Label2Display * hard_link_Label2Display(Label2Display * obj);


/* Function:  Label2Display_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Label2Display *]
 *
 */
Label2Display * Label2Display_alloc(void);


/* Function:  free_Label2Display(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Label2Display *]
 *
 * Return [UNKN ]  Undocumented return value [Label2Display *]
 *
 */
Label2Display * free_Label2Display(Label2Display * obj);


/* Function:  add_Sequence2Display(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Sequence2Display *]
 * Arg:        add [OWNER] Object to add to the list [Label2Display *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean add_Sequence2Display(Sequence2Display * obj,Label2Display * add);


/* Function:  flush_Sequence2Display(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [Sequence2Display *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_Sequence2Display(Sequence2Display * obj);


/* Function:  add_ind_Sequence2Display(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Sequence2Display *]
 * Arg:        add [OWNER] Object to add to the list [Index2Display *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean add_ind_Sequence2Display(Sequence2Display * obj,Index2Display * add);


/* Function:  flush_ind_Sequence2Display(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [Sequence2Display *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_ind_Sequence2Display(Sequence2Display * obj);


/* Function:  Sequence2Display_alloc_std(void)
 *
 * Descrip:    Equivalent to Sequence2Display_alloc_len(Sequence2DisplayLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Sequence2Display *]
 *
 */
Sequence2Display * Sequence2Display_alloc_std(void);


/* Function:  Sequence2Display_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [Sequence2Display *]
 *
 */
Sequence2Display * Sequence2Display_alloc_len(int len);


/* Function:  hard_link_Sequence2Display(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Sequence2Display *]
 *
 * Return [UNKN ]  Undocumented return value [Sequence2Display *]
 *
 */
Sequence2Display * hard_link_Sequence2Display(Sequence2Display * obj);


/* Function:  Sequence2Display_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Sequence2Display *]
 *
 */
Sequence2Display * Sequence2Display_alloc(void);


/* Function:  free_Sequence2Display(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Sequence2Display *]
 *
 * Return [UNKN ]  Undocumented return value [Sequence2Display *]
 *
 */
Sequence2Display * free_Sequence2Display(Sequence2Display * obj);


/* Function:  hard_link_Aln2DisplayResource(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Aln2DisplayResource *]
 *
 * Return [UNKN ]  Undocumented return value [Aln2DisplayResource *]
 *
 */
Aln2DisplayResource * hard_link_Aln2DisplayResource(Aln2DisplayResource * obj);


/* Function:  Aln2DisplayResource_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Aln2DisplayResource *]
 *
 */
Aln2DisplayResource * Aln2DisplayResource_alloc(void);


/* Function:  free_Aln2DisplayResource(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Aln2DisplayResource *]
 *
 * Return [UNKN ]  Undocumented return value [Aln2DisplayResource *]
 *
 */
Aln2DisplayResource * free_Aln2DisplayResource(Aln2DisplayResource * obj);


/* Function:  add_res_Aln2Display(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Aln2Display *]
 * Arg:        add [OWNER] Object to add to the list [Aln2DisplayResource *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean add_res_Aln2Display(Aln2Display * obj,Aln2DisplayResource * add);


/* Function:  flush_res_Aln2Display(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [Aln2Display *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_res_Aln2Display(Aln2Display * obj);


/* Function:  add_Aln2Display(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Aln2Display *]
 * Arg:        add [OWNER] Object to add to the list [Sequence2Display *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean add_Aln2Display(Aln2Display * obj,Sequence2Display * add);


/* Function:  flush_Aln2Display(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [Aln2Display *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_Aln2Display(Aln2Display * obj);


/* Function:  Aln2Display_alloc_std(void)
 *
 * Descrip:    Equivalent to Aln2Display_alloc_len(Aln2DisplayLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Aln2Display *]
 *
 */
Aln2Display * Aln2Display_alloc_std(void);


/* Function:  Aln2Display_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [Aln2Display *]
 *
 */
Aln2Display * Aln2Display_alloc_len(int len);


/* Function:  hard_link_Aln2Display(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Aln2Display *]
 *
 * Return [UNKN ]  Undocumented return value [Aln2Display *]
 *
 */
Aln2Display * hard_link_Aln2Display(Aln2Display * obj);


/* Function:  Aln2Display_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Aln2Display *]
 *
 */
Aln2Display * Aln2Display_alloc(void);


/* Function:  free_Aln2Display(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Aln2Display *]
 *
 * Return [UNKN ]  Undocumented return value [Aln2Display *]
 *
 */
Aln2Display * free_Aln2Display(Aln2Display * obj);


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
void write_Aln2Display_header(DYNFILE * dfp,Aln2Display * a2d);
void write_Aln2Display_function(DYNFILE * dfp,Aln2Display * a2d);
void write_Aln2Display_convert_func(DYNFILE * dfp,Aln2Display * a2d);
boolean prepare_Aln2Display(Aln2Display * a2d);
void number_up_Sequence_Index2Display(Aln2Display * a2d);
boolean crosslink_Index(Aln2Display * a2d);
boolean crosslink_Label2Display_Index(Label2Display * l2d,Sequence2Display * s2d);
Index2Display * Index2Display_from_name(Sequence2Display * s2d,char * name);
Aln2Display         * read_Aln2Display_line(char * line,FILE * ifp);
Sequence2Display    * read_Sequence2Display_line(char * line,FILE * ifp);
Label2Display       * read_Label2Display_line(char * line,FILE * ifp);
Aln2DisplayField    * read_Aln2DisplayField_line(char * line,FILE * ifp);
boolean read_name_line(Sequence2Display * s2d,char * name_line);
boolean put_away_Aln2DisplayField_strpair(Aln2DisplayField * a2df,char * pair);
Aln2DisplayResource * read_Aln2DisplayResource_line(char * line);
Index2Display * read_Index2Display_line(char * line);
void show_Aln2Display(Aln2Display * a2d,FILE * ofp);
void show_Aln2DisplayResource(Aln2DisplayResource * a2dr,FILE * ofp);
void show_Sequence2Display(Sequence2Display * s2d,FILE * ofp);
void show_Label2Display(Label2Display * l2d,FILE * ofp);
void show_Aln2DisplayField(Aln2DisplayField * a2df,FILE * ofp);


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void swap_Label2Display(Aln2DisplayField ** list,int i,int j) ;
void qsort_Label2Display(Aln2DisplayField ** list,int left,int right,int (*comp)(Aln2DisplayField * ,Aln2DisplayField * ));
void sort_Label2Display(Label2Display * obj,int (*comp)(Aln2DisplayField *, Aln2DisplayField *));
boolean expand_Label2Display(Label2Display * obj,int len);
void swap_ind_Label2Display(Index2Display ** list,int i,int j) ;
void qsort_ind_Label2Display(Index2Display ** list,int left,int right,int (*comp)(Index2Display * ,Index2Display * ));
void sort_ind_Label2Display(Label2Display * obj,int (*comp)(Index2Display *, Index2Display *));
boolean expand_ind_Label2Display(Label2Display * obj,int len);
void swap_Sequence2Display(Label2Display ** list,int i,int j) ;
void qsort_Sequence2Display(Label2Display ** list,int left,int right,int (*comp)(Label2Display * ,Label2Display * ));
void sort_Sequence2Display(Sequence2Display * obj,int (*comp)(Label2Display *, Label2Display *));
boolean expand_Sequence2Display(Sequence2Display * obj,int len);
void swap_ind_Sequence2Display(Index2Display ** list,int i,int j) ;
void qsort_ind_Sequence2Display(Index2Display ** list,int left,int right,int (*comp)(Index2Display * ,Index2Display * ));
void sort_ind_Sequence2Display(Sequence2Display * obj,int (*comp)(Index2Display *, Index2Display *));
boolean expand_ind_Sequence2Display(Sequence2Display * obj,int len);
void swap_res_Aln2Display(Aln2DisplayResource ** list,int i,int j) ;
void qsort_res_Aln2Display(Aln2DisplayResource ** list,int left,int right,int (*comp)(Aln2DisplayResource * ,Aln2DisplayResource * ));
void sort_res_Aln2Display(Aln2Display * obj,int (*comp)(Aln2DisplayResource *, Aln2DisplayResource *));
boolean expand_res_Aln2Display(Aln2Display * obj,int len);
void swap_Aln2Display(Sequence2Display ** list,int i,int j) ;
void qsort_Aln2Display(Sequence2Display ** list,int left,int right,int (*comp)(Sequence2Display * ,Sequence2Display * ));
void sort_Aln2Display(Aln2Display * obj,int (*comp)(Sequence2Display *, Sequence2Display *));
boolean expand_Aln2Display(Aln2Display * obj,int len);

#ifdef _cplusplus
}
#endif

#endif
