#ifndef DYNAMITEdynfileHEADERFILE
#define DYNAMITEdynfileHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "wisebase.h"
#include "funcinfo.h"
#include "inputfile.h"

#define DYNFILELISTLENGTH  1024

enum desctype {
  DYN_DESC_EDDY = 56
};


typedef int FuncOrderType;

enum fotype {
  FO_CUI_ALPHABETICAL = 83,
  FO_CUI_POSITION
};

struct DYNFILE {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int desc_type;   
    char line[1024];     
    char *tag[128];  
    int commentstart;    
    int  bracelevel;     
    boolean infunc;  
    char * sourceroot;   
    FILE * func;     
    FILE * head;     
    FILE * html;     
    FuncInfo ** info;    
    int len;/* len for above info  */ 
    int maxlen; /* maxlen for above info */ 
    ModuleInfo * mi;     
    int funcpos;     
    int code_debug_level;    
    boolean should_hard_link;    
    char * package_name;     
    } ;  
/* DYNFILE defined */ 
#ifndef DYNAMITE_DEFINED_DYNFILE
typedef struct DYNFILE DYNFILE;
#define DYNAMITE_DEFINED_DYNFILE
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  make_pdoc_output(dfp,mod_name,ofp)
 *
 * Descrip:    makes a pdoc for this file.
 *
 *
 *
 * Arg:             dfp [UNKN ] DYNFILE object [DYNFILE *]
 * Arg:        mod_name [UNKN ] Undocumented argument [char *]
 * Arg:             ofp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean make_pdoc_output(DYNFILE * dfp,char * mod_name,FILE * ofp);


/* Function:  have_got_unknown_func_DYNFILE(dfp,no,should_complain)
 *
 * Descrip:    reports back whether there are any undocumented
 *             file. if should_complain is TRUE, will issue
 *             warnings through warn. Writes the number of
 *             undocumented functions into no
 *
 *
 * Arg:                    dfp [UNKN ] DYNFILE object [DYNFILE *]
 * Arg:                     no [WRITE] pointer to some memory for the number of undoc'd funcs [int *]
 * Arg:        should_complain [UNKN ] if TRUE, will issue warn statements [boolean]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean have_got_unknown_func_DYNFILE(DYNFILE * dfp,int * no,boolean should_complain);


/* Function:  write_Dynamite_minimal_func(dfp)
 *
 * Descrip:    writes basically just the #include "self.h" line
 *             in the dynamite file
 *
 *
 * Arg:        dfp [UNKN ] Undocumented argument [DYNFILE *]
 *
 */
void write_Dynamite_minimal_func(DYNFILE * dfp);


/* Function:  FuncInfo_from_name_DYNFILE(dfp,name)
 *
 * Descrip:    Finds a funcinfo with name. returns NULL
 *             if it cant find them
 *
 *
 * Arg:         dfp [UNKN ] Undocumented argument [DYNFILE *]
 * Arg:        name [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [FuncInfo *]
 *
 */
FuncInfo * FuncInfo_from_name_DYNFILE(DYNFILE * dfp,char * name);


/* Function:  place_declarations_DYNFILE(dfp,type)
 *
 * Descrip:    writes the header declarations of functions
 *             into the header file, with documentation.
 *
 *             The FuncOrderType is either the more common,
 *             usual order or alphabetical
 *
 *
 * Arg:         dfp [UNKN ] Undocumented argument [DYNFILE *]
 * Arg:        type [UNKN ] Undocumented argument [FuncOrderType]
 *
 */
void place_declarations_DYNFILE(DYNFILE * dfp,FuncOrderType type);


/* Function:  positionify_DYNFILE(dfp)
 *
 * Descrip:    places the positions of functions as read in the file
 *             (hopefully) into the funcinfos so can be sorted on
 *             them (and rearranged potentially)
 *
 *
 * Arg:        dfp [UNKN ] Undocumented argument [DYNFILE *]
 *
 */
void positionify_DYNFILE(DYNFILE * dfp);


/* Function:  add_DYNFILE(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [DYNFILE *]
 * Arg:        add [OWNER] Object to add to the list [FuncInfo *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean add_DYNFILE(DYNFILE * obj,FuncInfo * add);


/* Function:  flush_DYNFILE(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [DYNFILE *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_DYNFILE(DYNFILE * obj);


/* Function:  DYNFILE_alloc_std(void)
 *
 * Descrip:    Equivalent to DYNFILE_alloc_len(DYNFILELISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DYNFILE *]
 *
 */
DYNFILE * DYNFILE_alloc_std(void);


/* Function:  DYNFILE_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [DYNFILE *]
 *
 */
DYNFILE * DYNFILE_alloc_len(int len);


/* Function:  hard_link_DYNFILE(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DYNFILE *]
 *
 * Return [UNKN ]  Undocumented return value [DYNFILE *]
 *
 */
DYNFILE * hard_link_DYNFILE(DYNFILE * obj);


/* Function:  DYNFILE_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DYNFILE *]
 *
 */
DYNFILE * DYNFILE_alloc(void);


/* Function:  free_DYNFILE(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DYNFILE *]
 *
 * Return [UNKN ]  Undocumented return value [DYNFILE *]
 *
 */
DYNFILE * free_DYNFILE(DYNFILE * obj);


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
void show_FuncInfo_DYNFILE(DYNFILE * dfp,FuncInfo * fi);
boolean real_is_unfinished_expression(char * str);
void flush_line(DYNFILE * dfp);
void flush_line_to_header(DYNFILE * dfp);
void pad_line(DYNFILE * dfp);
void add_break(DYNFILE * dfp);
void start_function(DYNFILE * dfp,char * str, ...);
void start_function_FuncInfo(FuncInfo * fi,DYNFILE * dfp,char * str, ... );
void true_start_function(DYNFILE * dfp,char * str, ... );
void close_function(DYNFILE * dfp);
void hang_expr(DYNFILE * dfp,char * str, ...);
void expr(DYNFILE * dfp,char * str, ...);
void macro(DYNFILE * dfp,char * str, ...);
void warn_expr(DYNFILE * dfp,char * str, ...);
void start_struct(DYNFILE * dfp,char * str, ...);
void close_struct(DYNFILE * dfp,char * str, ...);
void struct_macro(DYNFILE * dfp,char * str, ...);
void struct_expr(DYNFILE * dfp,char * str, ...);
void add_block_comment(DYNFILE * dfp,char * str, ...);
void add_end_comment(DYNFILE * dfp,char * str, ...);
void add_end_comment_header(DYNFILE * dfp,char * str, ...);
void set_commentstart(DYNFILE * dfp,int s);
void current_indent(DYNFILE * dfp);
void startcase(DYNFILE * dfp);
void closecase(DYNFILE * dfp);
void startbrace_tag(DYNFILE * dfp,char * str);
void startbrace(DYNFILE * dfp);
void closebrace(DYNFILE * dfp);
void fputs_func_DYNFILE(char * str,DYNFILE * dfp);
DYNFILE * open_std_no_html_dynfile(char * name);
DYNFILE * open_std_dynfile(char * name,char * html_directory);
void close_dynfile(DYNFILE * dfp);


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void show_html_DYNFILE(DYNFILE * dfp,FILE * ofp,FuncOrderType type);
int place_func_decl(DYNFILE * dfp,FuncInfo * fi);
void sort_DYNFILE_FuncOrderType(DYNFILE * dfp,FuncOrderType type);
boolean is_separated_types(FuncOrderType type);
int compare_two_FuncInfo_alpha(FuncInfo * one,FuncInfo * two);
int compare_two_FuncInfo_pos(FuncInfo * one,FuncInfo * two);
char * strinopen(char * target,char * probe);
void swap_DYNFILE(FuncInfo ** list,int i,int j) ;
void qsort_DYNFILE(FuncInfo ** list,int left,int right,int (*comp)(FuncInfo * ,FuncInfo * ));
void sort_DYNFILE(DYNFILE * obj,int (*comp)(FuncInfo *, FuncInfo *));
boolean expand_DYNFILE(DYNFILE * obj,int len);

#ifdef _cplusplus
}
#endif

#endif
