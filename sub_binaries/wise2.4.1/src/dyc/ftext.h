#ifndef DYNAMITEftextHEADERFILE
#define DYNAMITEftextHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "wisebase.h"

#define FblockLISTLENGTH 64
#define FtextLISTLENGTH 64

struct Fblock {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char ** line;    
    int len;/* len for above line  */ 
    int maxlen; /* maxlen for above line */ 
    } ;  
/* Fblock defined */ 
#ifndef DYNAMITE_DEFINED_Fblock
typedef struct Fblock Fblock;
#define DYNAMITE_DEFINED_Fblock
#endif


struct Ftext {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Fblock ** fb;    
    int len;/* len for above fb  */ 
    int maxlen; /* maxlen for above fb */ 
    } ;  
/* Ftext defined */ 
#ifndef DYNAMITE_DEFINED_Ftext
typedef struct Ftext Ftext;
#define DYNAMITE_DEFINED_Ftext
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  add_line_to_Ftext(ft,str,)
 *
 * Descrip:    adds a vsprintf'd line (below MAXLINE length!) to
 *             the Ftext.
 *
 *
 * Arg:         ft [UNKN ] Undocumented argument [Ftext *]
 * Arg:        str [UNKN ] Undocumented argument [char *]
 * Arg:            [UNKN ] Undocumented argument [.]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean add_line_to_Ftext(Ftext * ft,char * str,...);


/* Function:  add_break_to_Ftext(ft)
 *
 * Descrip:    puts in a break into the Ftext
 *
 *
 * Arg:        ft [UNKN ] Undocumented argument [Ftext *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean add_break_to_Ftext(Ftext * ft);


/* Function:  single_Ftext_from_str(str)
 *
 * Descrip:    Makes a complete Ftext from just this string:
 *
 *
 * Arg:        str [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [Ftext *]
 *
 */
Ftext * single_Ftext_from_str(char * str);


/* Function:  show_eddystyle_Ftext(ft,header,depth,ofp,blank_text)
 *
 * Descrip:    shows Ftext as a C style comment with
 *               *
 *               * 
 *               *
 *             indenting: 
 *
 *             will return number of lines printed.
 *
 *
 * Arg:                ft [READ ] Ftext to be shown [Ftext *]
 * Arg:            header [READ ] Header for the first line, eg, "description:" [char *]
 * Arg:             depth [READ ] depth of from * to text [int]
 * Arg:               ofp [UNKN ] output file [FILE *]
 * Arg:        blank_text [READ ] if non NULL, what to put if ft is empty, Can be NULL. [char *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int show_eddystyle_Ftext(Ftext * ft,char * header,int depth,FILE * ofp,char * blank_text);


/* Function:  latex_Ftext(ft,ofp)
 *
 * Descrip:    Provides a latex dump of some text.
 *
 *             Lines that start flush to the left
 *             are given as paragraphs
 *
 *             Line that are indented are made verbatim
 *
 *
 * Arg:         ft [UNKN ] Undocumented argument [Ftext *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void latex_Ftext(Ftext * ft,FILE * ofp);


/* Function:  dump_Ftext(ft,*ofp)
 *
 * Descrip:    stupid function which gives flat dump of ftext
 *
 *
 * Arg:          ft [UNKN ] Undocumented argument [Ftext *]
 * Arg:        *ofp [UNKN ] Undocumented argument [FILE]
 *
 */
void dump_Ftext(Ftext * ft,FILE *ofp);


/* Function:  dump_Ftext_pre(pre,ft,*ofp)
 *
 * Descrip:    stupid function which gives flat dump of ftext,
 *             with a start of pre 
 *
 *
 * Arg:         pre [UNKN ] Undocumented argument [char *]
 * Arg:          ft [UNKN ] Undocumented argument [Ftext *]
 * Arg:        *ofp [UNKN ] Undocumented argument [FILE]
 *
 */
void dump_Ftext_pre(char * pre,Ftext * ft,FILE *ofp);


/* Function:  read_Ftext(buffer,fgets_func,maxlen,ifp,endpoint)
 *
 * Descrip:    reads in lines until it hits endpoint.
 *             A bit of an internal: is going to use buffer and maxlen
 *             as its buffer (this is so it is not fixed to one length of
 *             buffer). You probably have abuffer in your calling function.
 *
 *
 * Arg:            buffer [WRITE] pointer to a char * buffer of maxlen that can be written in [char *]
 * Arg:        fgets_func [UNKN ] Undocumented argument [NullString]
 * Arg:            maxlen [READ ] maximum size of buffer to read [int]
 * Arg:               ifp [UNKN ] input file [FILE *]
 * Arg:          endpoint [READ ] a string which is used (using /strstartcmp) as a tag for the end of text [char *]
 *
 * Return [UNKN ]  Undocumented return value [Ftext *]
 *
 */
Ftext * read_Ftext(char * buffer,int maxlen,FILE * ifp,char * endpoint,char * (*fgets_func)(char *,int,FILE *));


/* Function:  add_Fblock(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Fblock *]
 * Arg:        add [OWNER] Object to add to the list [char *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean add_Fblock(Fblock * obj,char * add);


/* Function:  flush_Fblock(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [Fblock *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_Fblock(Fblock * obj);


/* Function:  Fblock_alloc_std(void)
 *
 * Descrip:    Equivalent to Fblock_alloc_len(FblockLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Fblock *]
 *
 */
Fblock * Fblock_alloc_std(void);


/* Function:  Fblock_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [Fblock *]
 *
 */
Fblock * Fblock_alloc_len(int len);


/* Function:  hard_link_Fblock(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Fblock *]
 *
 * Return [UNKN ]  Undocumented return value [Fblock *]
 *
 */
Fblock * hard_link_Fblock(Fblock * obj);


/* Function:  Fblock_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Fblock *]
 *
 */
Fblock * Fblock_alloc(void);


/* Function:  free_Fblock(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Fblock *]
 *
 * Return [UNKN ]  Undocumented return value [Fblock *]
 *
 */
Fblock * free_Fblock(Fblock * obj);


/* Function:  add_Ftext(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Ftext *]
 * Arg:        add [OWNER] Object to add to the list [Fblock *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean add_Ftext(Ftext * obj,Fblock * add);


/* Function:  flush_Ftext(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [Ftext *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_Ftext(Ftext * obj);


/* Function:  Ftext_alloc_std(void)
 *
 * Descrip:    Equivalent to Ftext_alloc_len(FtextLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Ftext *]
 *
 */
Ftext * Ftext_alloc_std(void);


/* Function:  Ftext_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [Ftext *]
 *
 */
Ftext * Ftext_alloc_len(int len);


/* Function:  hard_link_Ftext(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Ftext *]
 *
 * Return [UNKN ]  Undocumented return value [Ftext *]
 *
 */
Ftext * hard_link_Ftext(Ftext * obj);


/* Function:  Ftext_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Ftext *]
 *
 */
Ftext * Ftext_alloc(void);


/* Function:  free_Ftext(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Ftext *]
 *
 * Return [UNKN ]  Undocumented return value [Ftext *]
 *
 */
Ftext * free_Ftext(Ftext * obj);


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void dump_Fblock_str(char * pre,Fblock * fb,FILE * ofp);
Fblock * read_Fblock(char * buffer,int maxlen,FILE * ifp,char * endpoint,char * (*fgets_func)(char *,int,FILE *));
void swap_Fblock(char ** list,int i,int j) ;
void qsort_Fblock(char ** list,int left,int right,int (*comp)(char * ,char * ));
void sort_Fblock(Fblock * obj,int (*comp)(char *, char *));
boolean expand_Fblock(Fblock * obj,int len);
void swap_Ftext(Fblock ** list,int i,int j) ;
void qsort_Ftext(Fblock ** list,int left,int right,int (*comp)(Fblock * ,Fblock * ));
void sort_Ftext(Ftext * obj,int (*comp)(Fblock *, Fblock *));
boolean expand_Ftext(Ftext * obj,int len);

#ifdef _cplusplus
}
#endif

#endif
