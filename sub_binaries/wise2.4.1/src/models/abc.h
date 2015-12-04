#ifndef DYNAMITEabcHEADERFILE
#define DYNAMITEabcHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"



struct Wise2_abc {  
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
    CompMat* comp;   
    int a;   
    int b;   
    int c;   
    } ;  
/* abc defined */ 
#ifndef DYNAMITE_DEFINED_abc
typedef struct Wise2_abc Wise2_abc;
#define abc Wise2_abc
#define DYNAMITE_DEFINED_abc
#endif


#ifdef PTHREAD
struct thread_pool_holder_abc {  
    ComplexSequence* query; /* Query object placeholder */ 
    ProteinDB* querydb; /* Query database object */ 
    boolean query_init;  
    ComplexSequence* target;/* Target object placeholder */ 
    ProteinDB* targetdb;/* Target database object */ 
    boolean target_init; 
    CompMat* comp;   
    int a;   
    int b;   
    int c;   
    pthread_mutex_t input_lock;  
    pthread_mutex_t output_lock;     
    Hscore * out;    
    pthread_t * pool;    
    int number_of_threads;   
    boolean search_has_ended;    
    DBSearchImpl * dbsi;     
    } ;  
#endif /* PTHREAD */
struct Wise2_abc_access_func_holder {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int (*access_main)(abc*,int,int,int);    
    int (*access_special)(abc*,int,int,int); 
    } ;  
/* abc_access_func_holder defined */ 
#ifndef DYNAMITE_DEFINED_abc_access_func_holder
typedef struct Wise2_abc_access_func_holder Wise2_abc_access_func_holder;
#define abc_access_func_holder Wise2_abc_access_func_holder
#define DYNAMITE_DEFINED_abc_access_func_holder
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  PackAln_read_Shatter_abc(mat)
 *
 * Descrip:    Reads off PackAln from shatter matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [abc *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Shatter_abc(abc * mat);
#define PackAln_read_Shatter_abc Wise2_PackAln_read_Shatter_abc


/* Function:  calculate_shatter_abc(mat,dpenv)
 *
 * Descrip:    This function calculates the abc matrix when in shatter mode
 *
 *
 * Arg:          mat [UNKN ] (null) [abc *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_shatter_abc(abc * mat,DPEnvelope * dpenv);
#define calculate_shatter_abc Wise2_calculate_shatter_abc


/* Function:  search_abc(dbsi,out,querydb,targetdb,comp,a,b,c)
 *
 * Descrip:    This function makes a database search of abc
 *             It uses the dbsi structure to choose which implementation
 *             to use of the database searching. This way at run time you
 *             can switch between single threaded/multi-threaded or hardware
 *
 *
 * Arg:            dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 * Arg:             out [UNKN ] Undocumented argument [Hscore *]
 * Arg:         querydb [UNKN ] Undocumented argument [ProteinDB*]
 * Arg:        targetdb [UNKN ] Undocumented argument [ProteinDB*]
 * Arg:            comp [UNKN ] Undocumented argument [CompMat*]
 * Arg:               a [UNKN ] Undocumented argument [int]
 * Arg:               b [UNKN ] Undocumented argument [int]
 * Arg:               c [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_search_abc(DBSearchImpl * dbsi,Hscore * out,ProteinDB* querydb,ProteinDB* targetdb ,CompMat* comp,int a,int b,int c);
#define search_abc Wise2_search_abc


/* Function:  serial_search_abc(out,querydb,targetdb,comp,a,b,c)
 *
 * Descrip:    This function makes a database search of abc
 *             It is a single processor implementation
 *
 *
 * Arg:             out [UNKN ] Undocumented argument [Hscore *]
 * Arg:         querydb [UNKN ] Undocumented argument [ProteinDB*]
 * Arg:        targetdb [UNKN ] Undocumented argument [ProteinDB*]
 * Arg:            comp [UNKN ] Undocumented argument [CompMat*]
 * Arg:               a [UNKN ] Undocumented argument [int]
 * Arg:               b [UNKN ] Undocumented argument [int]
 * Arg:               c [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_serial_search_abc(Hscore * out,ProteinDB* querydb,ProteinDB* targetdb ,CompMat* comp,int a,int b,int c);
#define serial_search_abc Wise2_serial_search_abc


/* Function:  PackAln_bestmemory_abc(query,target,comp,a,b,c,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_abc
 *
 *
 * Arg:         query [UNKN ] query data structure [ComplexSequence*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:          comp [UNKN ] Resource [CompMat*]
 * Arg:             a [UNKN ] Resource [int]
 * Arg:             b [UNKN ] Resource [int]
 * Arg:             c [UNKN ] Resource [int]
 * Arg:         dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:          dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_bestmemory_abc(ComplexSequence* query,ComplexSequence* target ,CompMat* comp,int a,int b,int c,DPEnvelope * dpenv,DPRunImpl * dpri);
#define PackAln_bestmemory_abc Wise2_PackAln_bestmemory_abc


/* Function:  allocate_Expl_abc(query,target,comp,a,b,c,dpri)
 *
 * Descrip:    This function allocates the abc structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_abc_only
 *
 *
 * Arg:         query [UNKN ] query data structure [ComplexSequence*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:          comp [UNKN ] Resource [CompMat*]
 * Arg:             a [UNKN ] Resource [int]
 * Arg:             b [UNKN ] Resource [int]
 * Arg:             c [UNKN ] Resource [int]
 * Arg:          dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [abc *]
 *
 */
abc * Wise2_allocate_Expl_abc(ComplexSequence* query,ComplexSequence* target ,CompMat* comp,int a,int b,int c,DPRunImpl * dpri);
#define allocate_Expl_abc Wise2_allocate_Expl_abc


/* Function:  recalculate_PackAln_abc(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by abc
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [abc *]
 *
 */
void Wise2_recalculate_PackAln_abc(PackAln * pal,abc * mat);
#define recalculate_PackAln_abc Wise2_recalculate_PackAln_abc


/* Function:  allocate_Small_abc(query,target,comp,a,b,c)
 *
 * Descrip:    This function allocates the abc structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_abc_only
 *
 *
 * Arg:         query [UNKN ] query data structure [ComplexSequence*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:          comp [UNKN ] Resource [CompMat*]
 * Arg:             a [UNKN ] Resource [int]
 * Arg:             b [UNKN ] Resource [int]
 * Arg:             c [UNKN ] Resource [int]
 *
 * Return [UNKN ]  Undocumented return value [abc *]
 *
 */
abc * Wise2_allocate_Small_abc(ComplexSequence* query,ComplexSequence* target ,CompMat* comp,int a,int b,int c);
#define allocate_Small_abc Wise2_allocate_Small_abc


/* Function:  PackAln_calculate_Small_abc(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for abc structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_abc 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_abc 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [abc *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_calculate_Small_abc(abc * mat,DPEnvelope * dpenv);
#define PackAln_calculate_Small_abc Wise2_PackAln_calculate_Small_abc


/* Function:  AlnRangeSet_calculate_Small_abc(mat)
 *
 * Descrip:    This function calculates an alignment for abc structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_abc 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_abc
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_abc 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [abc *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_calculate_Small_abc(abc * mat);
#define AlnRangeSet_calculate_Small_abc Wise2_AlnRangeSet_calculate_Small_abc


/* Function:  AlnRangeSet_from_abc(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for abc structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_abc 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_abc
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [abc *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_from_abc(abc * mat);
#define AlnRangeSet_from_abc Wise2_AlnRangeSet_from_abc


/* Function:  convert_PackAln_to_AlnBlock_abc(pal)
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
AlnBlock * Wise2_convert_PackAln_to_AlnBlock_abc(PackAln * pal);
#define convert_PackAln_to_AlnBlock_abc Wise2_convert_PackAln_to_AlnBlock_abc


/* Function:  PackAln_read_Expl_abc(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [abc *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Expl_abc(abc * mat);
#define PackAln_read_Expl_abc Wise2_PackAln_read_Expl_abc


/* Function:  PackAln_read_generic_abc(mat,h)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [abc *]
 * Arg:          h [UNKN ] Undocumented argument [abc_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_generic_abc(abc * mat,abc_access_func_holder h);
#define PackAln_read_generic_abc Wise2_PackAln_read_generic_abc


/* Function:  calculate_abc(mat)
 *
 * Descrip:    This function calculates the abc matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_abc
 *
 *
 * Arg:        mat [UNKN ] abc which contains explicit basematrix memory [abc *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_abc(abc * mat);
#define calculate_abc Wise2_calculate_abc


/* Function:  calculate_dpenv_abc(mat,dpenv)
 *
 * Descrip:    This function calculates the abc matrix when in explicit mode, subject to the envelope
 *
 *
 * Arg:          mat [UNKN ] abc which contains explicit basematrix memory [abc *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_dpenv_abc(abc * mat,DPEnvelope * dpenv);
#define calculate_dpenv_abc Wise2_calculate_dpenv_abc


/* Function:  abc_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [abc *]
 *
 */
abc * Wise2_abc_alloc(void);
#define abc_alloc Wise2_abc_alloc


/* Function:  free_abc(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [abc *]
 *
 * Return [UNKN ]  Undocumented return value [abc *]
 *
 */
abc * Wise2_free_abc(abc * obj);
#define free_abc Wise2_free_abc


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
int Wise2_abc_shatter_access_main(abc * mat,int i,int j,int state);
#define abc_shatter_access_main Wise2_abc_shatter_access_main
int Wise2_abc_shatter_access_special(abc * mat,int i,int j,int state);
#define abc_shatter_access_special Wise2_abc_shatter_access_special
void * Wise2_thread_loop_abc(void * ptr);
#define thread_loop_abc Wise2_thread_loop_abc
int Wise2_score_only_abc(ComplexSequence* query,ComplexSequence* target ,CompMat* comp,int a,int b,int c);
#define score_only_abc Wise2_score_only_abc
abc * Wise2_allocate_abc_only(ComplexSequence* query,ComplexSequence* target ,CompMat* comp,int a,int b,int c);
#define allocate_abc_only Wise2_allocate_abc_only
void Wise2_init_abc(abc * mat);
#define init_abc Wise2_init_abc
AlnRange * Wise2_AlnRange_build_abc(abc * mat,int stopj,int stopspecstate,int * startj,int * startspecstate);
#define AlnRange_build_abc Wise2_AlnRange_build_abc
boolean Wise2_read_hidden_abc(abc * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out);
#define read_hidden_abc Wise2_read_hidden_abc
int Wise2_max_hidden_abc(abc * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_hidden_abc Wise2_max_hidden_abc
boolean Wise2_read_special_strip_abc(abc * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out);
#define read_special_strip_abc Wise2_read_special_strip_abc
int Wise2_max_special_strip_abc(abc * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_special_strip_abc Wise2_max_special_strip_abc
int Wise2_max_matrix_to_special_abc(abc * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_matrix_to_special_abc Wise2_max_matrix_to_special_abc
void Wise2_calculate_hidden_abc(abc * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv);
#define calculate_hidden_abc Wise2_calculate_hidden_abc
void Wise2_init_hidden_abc(abc * mat,int starti,int startj,int stopi,int stopj);
#define init_hidden_abc Wise2_init_hidden_abc
boolean Wise2_full_dc_abc(abc * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv);
#define full_dc_abc Wise2_full_dc_abc
boolean Wise2_do_dc_single_pass_abc(abc * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done);
#define do_dc_single_pass_abc Wise2_do_dc_single_pass_abc
void Wise2_push_dc_at_merge_abc(abc * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv);
#define push_dc_at_merge_abc Wise2_push_dc_at_merge_abc
void Wise2_follow_on_dc_abc(abc * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define follow_on_dc_abc Wise2_follow_on_dc_abc
void Wise2_run_up_dc_abc(abc * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define run_up_dc_abc Wise2_run_up_dc_abc
void Wise2_init_dc_abc(abc * mat);
#define init_dc_abc Wise2_init_dc_abc
int Wise2_start_end_find_end_abc(abc * mat,int * endj);
#define start_end_find_end_abc Wise2_start_end_find_end_abc
boolean Wise2_dc_optimised_start_end_calc_abc(abc *mat,DPEnvelope * dpenv);
#define dc_optimised_start_end_calc_abc Wise2_dc_optimised_start_end_calc_abc
void Wise2_init_start_end_linear_abc(abc * mat);
#define init_start_end_linear_abc Wise2_init_start_end_linear_abc
AlnConvertSet * Wise2_AlnConvertSet_abc(void);
#define AlnConvertSet_abc Wise2_AlnConvertSet_abc
int Wise2_abc_explicit_access_main(abc * mat,int i,int j,int state);
#define abc_explicit_access_main Wise2_abc_explicit_access_main
int Wise2_abc_explicit_access_special(abc * mat,int i,int j,int state);
#define abc_explicit_access_special Wise2_abc_explicit_access_special
int Wise2_find_end_abc(abc * mat,int * ri,int * rj,int * state,boolean * isspecial,abc_access_func_holder h);
#define find_end_abc Wise2_find_end_abc
void Wise2_abc_debug_show_matrix(abc * mat,int starti,int stopi,int startj,int stopj,FILE * ofp);
#define abc_debug_show_matrix Wise2_abc_debug_show_matrix
int Wise2_max_calc_abc(abc * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,abc_access_func_holder h);
#define max_calc_abc Wise2_max_calc_abc
int Wise2_max_calc_special_abc(abc * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,abc_access_func_holder h);
#define max_calc_special_abc Wise2_max_calc_special_abc

#ifdef _cplusplus
}
#endif

#endif
