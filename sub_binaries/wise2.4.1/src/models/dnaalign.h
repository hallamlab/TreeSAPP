#ifndef DYNAMITEdnaalignHEADERFILE
#define DYNAMITEdnaalignHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"

enum DnaStartEndEnum {
  DSE_GLOBAL_START = 0,
  DSE_GLOBAL_END,
  DSE_EDGE_START,
  DSE_EDGE_END,
  DSE_LOCAL_START,
  DSE_LOCAL_END,
  DSE_NUMBER
};

struct Wise2_DnaStartEnd {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int trans[DSE_NUMBER];  /*  start/end points possibilities */ 
    } ;  
/* DnaStartEnd defined */ 
#ifndef DYNAMITE_DEFINED_DnaStartEnd
typedef struct Wise2_DnaStartEnd Wise2_DnaStartEnd;
#define DnaStartEnd Wise2_DnaStartEnd
#define DYNAMITE_DEFINED_DnaStartEnd
#endif


struct Wise2_DnaAlign {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    BaseMatrix * basematrix;     
    ShatterMatrix * shatter;     
    int leni;    
    int lenj;    
    ComplexSequence* query;  
    ComplexSequence* target;     
    DnaMatrix* comp;     
    Score qgap;  
    Score qext;  
    Score tgap;  
    Score text;  
    DnaStartEnd * startend;  
    } ;  
/* DnaAlign defined */ 
#ifndef DYNAMITE_DEFINED_DnaAlign
typedef struct Wise2_DnaAlign Wise2_DnaAlign;
#define DnaAlign Wise2_DnaAlign
#define DYNAMITE_DEFINED_DnaAlign
#endif


struct Wise2_DnaAlign_access_func_holder {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int (*access_main)(DnaAlign*,int,int,int);   
    int (*access_special)(DnaAlign*,int,int,int);    
    } ;  
/* DnaAlign_access_func_holder defined */ 
#ifndef DYNAMITE_DEFINED_DnaAlign_access_func_holder
typedef struct Wise2_DnaAlign_access_func_holder Wise2_DnaAlign_access_func_holder;
#define DnaAlign_access_func_holder Wise2_DnaAlign_access_func_holder
#define DYNAMITE_DEFINED_DnaAlign_access_func_holder
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  make_align_dnaalign(one,two,mat,se,qgap,qext,tgap,text,dpri)
 *
 * Descrip:    Makes an alignment out of two DNA sequences
 *
 *
 * Arg:         one [READ ] first sequence to align [Sequence *]
 * Arg:         two [READ ] second sequence to align [Sequence *]
 * Arg:         mat [READ ] DnaMatrix for the matching [DnaMatrix *]
 * Arg:          se [READ ] DnaStartEnd policy [DnaStartEnd *]
 * Arg:        qgap [READ ] gap open penalty in query (one) coordinate [int]
 * Arg:        qext [READ ] gap extension penalty in query (one) coordinate [int]
 * Arg:        tgap [READ ] gap open penalty in target (two) coordinate [int]
 * Arg:        text [READ ] gap extension penalty in target (two) coordinate [int]
 * Arg:        dpri [READ ] DPRunImpl structure [DPRunImpl *]
 *
 * Return [UNKN ]  an alb structure of the alignment [AlnBlock *]
 *
 */
AlnBlock * Wise2_make_align_dnaalign(Sequence * one,Sequence * two,DnaMatrix * mat,DnaStartEnd * se,int qgap,int qext,int tgap,int text,DPRunImpl * dpri);
#define make_align_dnaalign Wise2_make_align_dnaalign


/* Function:  DnaStartEnd_from_policy(policy)
 *
 * Descrip:    Makes a DnaStartEnd from a particular string.
 *             Possible strings are:
 *
 *             local - fully local
 *             global - fully global
 *             edge - aligns only to edges
 *
 *
 * Arg:        policy [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [DnaStartEnd *]
 *
 */
DnaStartEnd * Wise2_DnaStartEnd_from_policy(char * policy);
#define DnaStartEnd_from_policy Wise2_DnaStartEnd_from_policy


/* Function:  hard_link_DnaStartEnd(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DnaStartEnd *]
 *
 * Return [UNKN ]  Undocumented return value [DnaStartEnd *]
 *
 */
DnaStartEnd * Wise2_hard_link_DnaStartEnd(DnaStartEnd * obj);
#define hard_link_DnaStartEnd Wise2_hard_link_DnaStartEnd


/* Function:  DnaStartEnd_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaStartEnd *]
 *
 */
DnaStartEnd * Wise2_DnaStartEnd_alloc(void);
#define DnaStartEnd_alloc Wise2_DnaStartEnd_alloc


/* Function:  free_DnaStartEnd(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DnaStartEnd *]
 *
 * Return [UNKN ]  Undocumented return value [DnaStartEnd *]
 *
 */
DnaStartEnd * Wise2_free_DnaStartEnd(DnaStartEnd * obj);
#define free_DnaStartEnd Wise2_free_DnaStartEnd


/* Function:  PackAln_read_Shatter_DnaAlign(mat)
 *
 * Descrip:    Reads off PackAln from shatter matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [DnaAlign *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Shatter_DnaAlign(DnaAlign * mat);
#define PackAln_read_Shatter_DnaAlign Wise2_PackAln_read_Shatter_DnaAlign


/* Function:  calculate_shatter_DnaAlign(mat,dpenv)
 *
 * Descrip:    This function calculates the DnaAlign matrix when in shatter mode
 *
 *
 * Arg:          mat [UNKN ] (null) [DnaAlign *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_shatter_DnaAlign(DnaAlign * mat,DPEnvelope * dpenv);
#define calculate_shatter_DnaAlign Wise2_calculate_shatter_DnaAlign


/* Function:  search_DnaAlign(dbsi,out,query,target,comp,qgap,qext,tgap,text,startend)
 *
 * Descrip:    This function makes a database search of DnaAlign
 *             It uses the dbsi structure to choose which implementation
 *             to use of the database searching. This way at run time you
 *             can switch between single threaded/multi-threaded or hardware
 *
 *
 * Arg:            dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 * Arg:             out [UNKN ] Undocumented argument [Hscore *]
 * Arg:           query [UNKN ] Undocumented argument [ComplexSequence*]
 * Arg:          target [UNKN ] Undocumented argument [ComplexSequence*]
 * Arg:            comp [UNKN ] Undocumented argument [DnaMatrix*]
 * Arg:            qgap [UNKN ] Undocumented argument [Score]
 * Arg:            qext [UNKN ] Undocumented argument [Score]
 * Arg:            tgap [UNKN ] Undocumented argument [Score]
 * Arg:            text [UNKN ] Undocumented argument [Score]
 * Arg:        startend [UNKN ] Undocumented argument [DnaStartEnd *]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_search_DnaAlign(DBSearchImpl * dbsi,Hscore * out,ComplexSequence* query,ComplexSequence* target ,DnaMatrix* comp,Score qgap,Score qext,Score tgap,Score text,DnaStartEnd * startend);
#define search_DnaAlign Wise2_search_DnaAlign


/* Function:  serial_search_DnaAlign(out,query,target,comp,qgap,qext,tgap,text,startend)
 *
 * Descrip:    This function makes a database search of DnaAlign
 *             It is a single processor implementation
 *
 *
 * Arg:             out [UNKN ] Undocumented argument [Hscore *]
 * Arg:           query [UNKN ] Undocumented argument [ComplexSequence*]
 * Arg:          target [UNKN ] Undocumented argument [ComplexSequence*]
 * Arg:            comp [UNKN ] Undocumented argument [DnaMatrix*]
 * Arg:            qgap [UNKN ] Undocumented argument [Score]
 * Arg:            qext [UNKN ] Undocumented argument [Score]
 * Arg:            tgap [UNKN ] Undocumented argument [Score]
 * Arg:            text [UNKN ] Undocumented argument [Score]
 * Arg:        startend [UNKN ] Undocumented argument [DnaStartEnd *]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_serial_search_DnaAlign(Hscore * out,ComplexSequence* query,ComplexSequence* target ,DnaMatrix* comp,Score qgap,Score qext,Score tgap,Score text,DnaStartEnd * startend);
#define serial_search_DnaAlign Wise2_serial_search_DnaAlign


/* Function:  PackAln_bestmemory_DnaAlign(query,target,comp,qgap,qext,tgap,text,startend,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_DnaAlign
 *
 *
 * Arg:           query [UNKN ] query data structure [ComplexSequence*]
 * Arg:          target [UNKN ] target data structure [ComplexSequence*]
 * Arg:            comp [UNKN ] Resource [DnaMatrix*]
 * Arg:            qgap [UNKN ] Resource [Score]
 * Arg:            qext [UNKN ] Resource [Score]
 * Arg:            tgap [UNKN ] Resource [Score]
 * Arg:            text [UNKN ] Resource [Score]
 * Arg:        startend [UNKN ] Resource [DnaStartEnd *]
 * Arg:           dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:            dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_bestmemory_DnaAlign(ComplexSequence* query,ComplexSequence* target ,DnaMatrix* comp,Score qgap,Score qext,Score tgap,Score text,DnaStartEnd * startend,DPEnvelope * dpenv,DPRunImpl * dpri);
#define PackAln_bestmemory_DnaAlign Wise2_PackAln_bestmemory_DnaAlign


/* Function:  allocate_Expl_DnaAlign(query,target,comp,qgap,qext,tgap,text,startend,dpri)
 *
 * Descrip:    This function allocates the DnaAlign structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_DnaAlign_only
 *
 *
 * Arg:           query [UNKN ] query data structure [ComplexSequence*]
 * Arg:          target [UNKN ] target data structure [ComplexSequence*]
 * Arg:            comp [UNKN ] Resource [DnaMatrix*]
 * Arg:            qgap [UNKN ] Resource [Score]
 * Arg:            qext [UNKN ] Resource [Score]
 * Arg:            tgap [UNKN ] Resource [Score]
 * Arg:            text [UNKN ] Resource [Score]
 * Arg:        startend [UNKN ] Resource [DnaStartEnd *]
 * Arg:            dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [DnaAlign *]
 *
 */
DnaAlign * Wise2_allocate_Expl_DnaAlign(ComplexSequence* query,ComplexSequence* target ,DnaMatrix* comp,Score qgap,Score qext,Score tgap,Score text,DnaStartEnd * startend,DPRunImpl * dpri);
#define allocate_Expl_DnaAlign Wise2_allocate_Expl_DnaAlign


/* Function:  recalculate_PackAln_DnaAlign(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by DnaAlign
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [DnaAlign *]
 *
 */
void Wise2_recalculate_PackAln_DnaAlign(PackAln * pal,DnaAlign * mat);
#define recalculate_PackAln_DnaAlign Wise2_recalculate_PackAln_DnaAlign


/* Function:  allocate_Small_DnaAlign(query,target,comp,qgap,qext,tgap,text,startend)
 *
 * Descrip:    This function allocates the DnaAlign structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_DnaAlign_only
 *
 *
 * Arg:           query [UNKN ] query data structure [ComplexSequence*]
 * Arg:          target [UNKN ] target data structure [ComplexSequence*]
 * Arg:            comp [UNKN ] Resource [DnaMatrix*]
 * Arg:            qgap [UNKN ] Resource [Score]
 * Arg:            qext [UNKN ] Resource [Score]
 * Arg:            tgap [UNKN ] Resource [Score]
 * Arg:            text [UNKN ] Resource [Score]
 * Arg:        startend [UNKN ] Resource [DnaStartEnd *]
 *
 * Return [UNKN ]  Undocumented return value [DnaAlign *]
 *
 */
DnaAlign * Wise2_allocate_Small_DnaAlign(ComplexSequence* query,ComplexSequence* target ,DnaMatrix* comp,Score qgap,Score qext,Score tgap,Score text,DnaStartEnd * startend);
#define allocate_Small_DnaAlign Wise2_allocate_Small_DnaAlign


/* Function:  PackAln_calculate_Small_DnaAlign(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for DnaAlign structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_DnaAlign 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_DnaAlign 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [DnaAlign *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_calculate_Small_DnaAlign(DnaAlign * mat,DPEnvelope * dpenv);
#define PackAln_calculate_Small_DnaAlign Wise2_PackAln_calculate_Small_DnaAlign


/* Function:  AlnRangeSet_calculate_Small_DnaAlign(mat)
 *
 * Descrip:    This function calculates an alignment for DnaAlign structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_DnaAlign 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_DnaAlign
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_DnaAlign 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [DnaAlign *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_calculate_Small_DnaAlign(DnaAlign * mat);
#define AlnRangeSet_calculate_Small_DnaAlign Wise2_AlnRangeSet_calculate_Small_DnaAlign


/* Function:  AlnRangeSet_from_DnaAlign(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for DnaAlign structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_DnaAlign 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_DnaAlign
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [DnaAlign *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_from_DnaAlign(DnaAlign * mat);
#define AlnRangeSet_from_DnaAlign Wise2_AlnRangeSet_from_DnaAlign


/* Function:  convert_PackAln_to_AlnBlock_DnaAlign(pal)
 *
 * Descrip:    Converts a path alignment to a label alignment
 *             The label alignment is probably much more useful than the path
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 *
 * Return [UNKN ]  Undocumented return value [AlnBlock *]
 *
 */
AlnBlock * Wise2_convert_PackAln_to_AlnBlock_DnaAlign(PackAln * pal);
#define convert_PackAln_to_AlnBlock_DnaAlign Wise2_convert_PackAln_to_AlnBlock_DnaAlign


/* Function:  PackAln_read_Expl_DnaAlign(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [DnaAlign *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Expl_DnaAlign(DnaAlign * mat);
#define PackAln_read_Expl_DnaAlign Wise2_PackAln_read_Expl_DnaAlign


/* Function:  PackAln_read_generic_DnaAlign(mat,h)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [DnaAlign *]
 * Arg:          h [UNKN ] Undocumented argument [DnaAlign_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_generic_DnaAlign(DnaAlign * mat,DnaAlign_access_func_holder h);
#define PackAln_read_generic_DnaAlign Wise2_PackAln_read_generic_DnaAlign


/* Function:  calculate_DnaAlign(mat)
 *
 * Descrip:    This function calculates the DnaAlign matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_DnaAlign
 *
 *
 * Arg:        mat [UNKN ] DnaAlign which contains explicit basematrix memory [DnaAlign *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_DnaAlign(DnaAlign * mat);
#define calculate_DnaAlign Wise2_calculate_DnaAlign


/* Function:  calculate_dpenv_DnaAlign(mat,dpenv)
 *
 * Descrip:    This function calculates the DnaAlign matrix when in explicit mode, subject to the envelope
 *
 *
 * Arg:          mat [UNKN ] DnaAlign which contains explicit basematrix memory [DnaAlign *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_dpenv_DnaAlign(DnaAlign * mat,DPEnvelope * dpenv);
#define calculate_dpenv_DnaAlign Wise2_calculate_dpenv_DnaAlign


/* Function:  DnaAlign_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaAlign *]
 *
 */
DnaAlign * Wise2_DnaAlign_alloc(void);
#define DnaAlign_alloc Wise2_DnaAlign_alloc


/* Function:  free_DnaAlign(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DnaAlign *]
 *
 * Return [UNKN ]  Undocumented return value [DnaAlign *]
 *
 */
DnaAlign * Wise2_free_DnaAlign(DnaAlign * obj);
#define free_DnaAlign Wise2_free_DnaAlign


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
int Wise2_DnaAlign_shatter_access_main(DnaAlign * mat,int i,int j,int state);
#define DnaAlign_shatter_access_main Wise2_DnaAlign_shatter_access_main
int Wise2_DnaAlign_shatter_access_special(DnaAlign * mat,int i,int j,int state);
#define DnaAlign_shatter_access_special Wise2_DnaAlign_shatter_access_special
int Wise2_score_only_DnaAlign(ComplexSequence* query,ComplexSequence* target ,DnaMatrix* comp,Score qgap,Score qext,Score tgap,Score text,DnaStartEnd * startend);
#define score_only_DnaAlign Wise2_score_only_DnaAlign
DnaAlign * Wise2_allocate_DnaAlign_only(ComplexSequence* query,ComplexSequence* target ,DnaMatrix* comp,Score qgap,Score qext,Score tgap,Score text,DnaStartEnd * startend);
#define allocate_DnaAlign_only Wise2_allocate_DnaAlign_only
void Wise2_init_DnaAlign(DnaAlign * mat);
#define init_DnaAlign Wise2_init_DnaAlign
AlnRange * Wise2_AlnRange_build_DnaAlign(DnaAlign * mat,int stopj,int stopspecstate,int * startj,int * startspecstate);
#define AlnRange_build_DnaAlign Wise2_AlnRange_build_DnaAlign
boolean Wise2_read_hidden_DnaAlign(DnaAlign * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out);
#define read_hidden_DnaAlign Wise2_read_hidden_DnaAlign
int Wise2_max_hidden_DnaAlign(DnaAlign * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_hidden_DnaAlign Wise2_max_hidden_DnaAlign
boolean Wise2_read_special_strip_DnaAlign(DnaAlign * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out);
#define read_special_strip_DnaAlign Wise2_read_special_strip_DnaAlign
int Wise2_max_special_strip_DnaAlign(DnaAlign * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_special_strip_DnaAlign Wise2_max_special_strip_DnaAlign
int Wise2_max_matrix_to_special_DnaAlign(DnaAlign * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_matrix_to_special_DnaAlign Wise2_max_matrix_to_special_DnaAlign
void Wise2_calculate_hidden_DnaAlign(DnaAlign * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv);
#define calculate_hidden_DnaAlign Wise2_calculate_hidden_DnaAlign
void Wise2_init_hidden_DnaAlign(DnaAlign * mat,int starti,int startj,int stopi,int stopj);
#define init_hidden_DnaAlign Wise2_init_hidden_DnaAlign
boolean Wise2_full_dc_DnaAlign(DnaAlign * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv);
#define full_dc_DnaAlign Wise2_full_dc_DnaAlign
boolean Wise2_do_dc_single_pass_DnaAlign(DnaAlign * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done);
#define do_dc_single_pass_DnaAlign Wise2_do_dc_single_pass_DnaAlign
void Wise2_push_dc_at_merge_DnaAlign(DnaAlign * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv);
#define push_dc_at_merge_DnaAlign Wise2_push_dc_at_merge_DnaAlign
void Wise2_follow_on_dc_DnaAlign(DnaAlign * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define follow_on_dc_DnaAlign Wise2_follow_on_dc_DnaAlign
void Wise2_run_up_dc_DnaAlign(DnaAlign * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define run_up_dc_DnaAlign Wise2_run_up_dc_DnaAlign
void Wise2_init_dc_DnaAlign(DnaAlign * mat);
#define init_dc_DnaAlign Wise2_init_dc_DnaAlign
int Wise2_start_end_find_end_DnaAlign(DnaAlign * mat,int * endj);
#define start_end_find_end_DnaAlign Wise2_start_end_find_end_DnaAlign
boolean Wise2_dc_optimised_start_end_calc_DnaAlign(DnaAlign *mat,DPEnvelope * dpenv);
#define dc_optimised_start_end_calc_DnaAlign Wise2_dc_optimised_start_end_calc_DnaAlign
void Wise2_init_start_end_linear_DnaAlign(DnaAlign * mat);
#define init_start_end_linear_DnaAlign Wise2_init_start_end_linear_DnaAlign
AlnConvertSet * Wise2_AlnConvertSet_DnaAlign(void);
#define AlnConvertSet_DnaAlign Wise2_AlnConvertSet_DnaAlign
int Wise2_DnaAlign_explicit_access_main(DnaAlign * mat,int i,int j,int state);
#define DnaAlign_explicit_access_main Wise2_DnaAlign_explicit_access_main
int Wise2_DnaAlign_explicit_access_special(DnaAlign * mat,int i,int j,int state);
#define DnaAlign_explicit_access_special Wise2_DnaAlign_explicit_access_special
int Wise2_find_end_DnaAlign(DnaAlign * mat,int * ri,int * rj,int * state,boolean * isspecial,DnaAlign_access_func_holder h);
#define find_end_DnaAlign Wise2_find_end_DnaAlign
void Wise2_DnaAlign_debug_show_matrix(DnaAlign * mat,int starti,int stopi,int startj,int stopj,FILE * ofp);
#define DnaAlign_debug_show_matrix Wise2_DnaAlign_debug_show_matrix
int Wise2_max_calc_DnaAlign(DnaAlign * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,DnaAlign_access_func_holder h);
#define max_calc_DnaAlign Wise2_max_calc_DnaAlign
int Wise2_max_calc_special_DnaAlign(DnaAlign * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,DnaAlign_access_func_holder h);
#define max_calc_special_DnaAlign Wise2_max_calc_special_DnaAlign

#ifdef _cplusplus
}
#endif

#endif
