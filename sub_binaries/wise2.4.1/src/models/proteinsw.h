#ifndef DYNAMITEproteinswHEADERFILE
#define DYNAMITEproteinswHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"

struct Wise2_ProteinSW {  
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
    int gap;     
    int ext;     
    } ;  
/* ProteinSW defined */ 
#ifndef DYNAMITE_DEFINED_ProteinSW
typedef struct Wise2_ProteinSW Wise2_ProteinSW;
#define ProteinSW Wise2_ProteinSW
#define DYNAMITE_DEFINED_ProteinSW
#endif


#ifdef PTHREAD
struct thread_pool_holder_ProteinSW {  
    ComplexSequence* query; /* Query object placeholder */ 
    ProteinDB* querydb; /* Query database object */ 
    boolean query_init;  
    ComplexSequence* target;/* Target object placeholder */ 
    ProteinDB* targetdb;/* Target database object */ 
    boolean target_init; 
    CompMat* comp;   
    int gap;     
    int ext;     
    pthread_mutex_t input_lock;  
    pthread_mutex_t output_lock;     
    Hscore * out;    
    pthread_t * pool;    
    int number_of_threads;   
    boolean search_has_ended;    
    DBSearchImpl * dbsi;     
    } ;  
#endif /* PTHREAD */
struct Wise2_ProteinSW_access_func_holder {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int (*access_main)(ProteinSW*,int,int,int);  
    int (*access_special)(ProteinSW*,int,int,int);   
    } ;  
/* ProteinSW_access_func_holder defined */ 
#ifndef DYNAMITE_DEFINED_ProteinSW_access_func_holder
typedef struct Wise2_ProteinSW_access_func_holder Wise2_ProteinSW_access_func_holder;
#define ProteinSW_access_func_holder Wise2_ProteinSW_access_func_holder
#define DYNAMITE_DEFINED_ProteinSW_access_func_holder
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  PackAln_read_Shatter_ProteinSW(mat)
 *
 * Descrip:    Reads off PackAln from shatter matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [ProteinSW *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Shatter_ProteinSW(ProteinSW * mat);
#define PackAln_read_Shatter_ProteinSW Wise2_PackAln_read_Shatter_ProteinSW


/* Function:  calculate_shatter_ProteinSW(mat,dpenv)
 *
 * Descrip:    This function calculates the ProteinSW matrix when in shatter mode
 *
 *
 * Arg:          mat [UNKN ] (null) [ProteinSW *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_shatter_ProteinSW(ProteinSW * mat,DPEnvelope * dpenv);
#define calculate_shatter_ProteinSW Wise2_calculate_shatter_ProteinSW


/* Function:  search_ProteinSW(dbsi,out,querydb,targetdb,comp,gap,ext)
 *
 * Descrip:    This function makes a database search of ProteinSW
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
 * Arg:             gap [UNKN ] Undocumented argument [int]
 * Arg:             ext [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_search_ProteinSW(DBSearchImpl * dbsi,Hscore * out,ProteinDB* querydb,ProteinDB* targetdb ,CompMat* comp,int gap,int ext);
#define search_ProteinSW Wise2_search_ProteinSW


/* Function:  serial_search_ProteinSW(out,querydb,targetdb,comp,gap,ext)
 *
 * Descrip:    This function makes a database search of ProteinSW
 *             It is a single processor implementation
 *
 *
 * Arg:             out [UNKN ] Undocumented argument [Hscore *]
 * Arg:         querydb [UNKN ] Undocumented argument [ProteinDB*]
 * Arg:        targetdb [UNKN ] Undocumented argument [ProteinDB*]
 * Arg:            comp [UNKN ] Undocumented argument [CompMat*]
 * Arg:             gap [UNKN ] Undocumented argument [int]
 * Arg:             ext [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_serial_search_ProteinSW(Hscore * out,ProteinDB* querydb,ProteinDB* targetdb ,CompMat* comp,int gap,int ext);
#define serial_search_ProteinSW Wise2_serial_search_ProteinSW


/* Function:  PackAln_bestmemory_ProteinSW(query,target,comp,gap,ext,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_ProteinSW
 *
 *
 * Arg:         query [UNKN ] query data structure [ComplexSequence*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:          comp [UNKN ] Resource [CompMat*]
 * Arg:           gap [UNKN ] Resource [int]
 * Arg:           ext [UNKN ] Resource [int]
 * Arg:         dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:          dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_bestmemory_ProteinSW(ComplexSequence* query,ComplexSequence* target ,CompMat* comp,int gap,int ext,DPEnvelope * dpenv,DPRunImpl * dpri);
#define PackAln_bestmemory_ProteinSW Wise2_PackAln_bestmemory_ProteinSW


/* Function:  allocate_Expl_ProteinSW(query,target,comp,gap,ext,dpri)
 *
 * Descrip:    This function allocates the ProteinSW structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_ProteinSW_only
 *
 *
 * Arg:         query [UNKN ] query data structure [ComplexSequence*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:          comp [UNKN ] Resource [CompMat*]
 * Arg:           gap [UNKN ] Resource [int]
 * Arg:           ext [UNKN ] Resource [int]
 * Arg:          dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [ProteinSW *]
 *
 */
ProteinSW * Wise2_allocate_Expl_ProteinSW(ComplexSequence* query,ComplexSequence* target ,CompMat* comp,int gap,int ext,DPRunImpl * dpri);
#define allocate_Expl_ProteinSW Wise2_allocate_Expl_ProteinSW


/* Function:  recalculate_PackAln_ProteinSW(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by ProteinSW
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [ProteinSW *]
 *
 */
void Wise2_recalculate_PackAln_ProteinSW(PackAln * pal,ProteinSW * mat);
#define recalculate_PackAln_ProteinSW Wise2_recalculate_PackAln_ProteinSW


/* Function:  allocate_Small_ProteinSW(query,target,comp,gap,ext)
 *
 * Descrip:    This function allocates the ProteinSW structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_ProteinSW_only
 *
 *
 * Arg:         query [UNKN ] query data structure [ComplexSequence*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:          comp [UNKN ] Resource [CompMat*]
 * Arg:           gap [UNKN ] Resource [int]
 * Arg:           ext [UNKN ] Resource [int]
 *
 * Return [UNKN ]  Undocumented return value [ProteinSW *]
 *
 */
ProteinSW * Wise2_allocate_Small_ProteinSW(ComplexSequence* query,ComplexSequence* target ,CompMat* comp,int gap,int ext);
#define allocate_Small_ProteinSW Wise2_allocate_Small_ProteinSW


/* Function:  PackAln_calculate_Small_ProteinSW(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for ProteinSW structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_ProteinSW 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_ProteinSW 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [ProteinSW *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_calculate_Small_ProteinSW(ProteinSW * mat,DPEnvelope * dpenv);
#define PackAln_calculate_Small_ProteinSW Wise2_PackAln_calculate_Small_ProteinSW


/* Function:  AlnRangeSet_calculate_Small_ProteinSW(mat)
 *
 * Descrip:    This function calculates an alignment for ProteinSW structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_ProteinSW 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_ProteinSW
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_ProteinSW 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [ProteinSW *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_calculate_Small_ProteinSW(ProteinSW * mat);
#define AlnRangeSet_calculate_Small_ProteinSW Wise2_AlnRangeSet_calculate_Small_ProteinSW


/* Function:  AlnRangeSet_from_ProteinSW(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for ProteinSW structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_ProteinSW 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_ProteinSW
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [ProteinSW *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_from_ProteinSW(ProteinSW * mat);
#define AlnRangeSet_from_ProteinSW Wise2_AlnRangeSet_from_ProteinSW


/* Function:  convert_PackAln_to_AlnBlock_ProteinSW(pal)
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
AlnBlock * Wise2_convert_PackAln_to_AlnBlock_ProteinSW(PackAln * pal);
#define convert_PackAln_to_AlnBlock_ProteinSW Wise2_convert_PackAln_to_AlnBlock_ProteinSW


/* Function:  PackAln_read_Expl_ProteinSW(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [ProteinSW *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Expl_ProteinSW(ProteinSW * mat);
#define PackAln_read_Expl_ProteinSW Wise2_PackAln_read_Expl_ProteinSW


/* Function:  PackAln_read_generic_ProteinSW(mat,h)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [ProteinSW *]
 * Arg:          h [UNKN ] Undocumented argument [ProteinSW_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_generic_ProteinSW(ProteinSW * mat,ProteinSW_access_func_holder h);
#define PackAln_read_generic_ProteinSW Wise2_PackAln_read_generic_ProteinSW


/* Function:  calculate_ProteinSW(mat)
 *
 * Descrip:    This function calculates the ProteinSW matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_ProteinSW
 *
 *
 * Arg:        mat [UNKN ] ProteinSW which contains explicit basematrix memory [ProteinSW *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_ProteinSW(ProteinSW * mat);
#define calculate_ProteinSW Wise2_calculate_ProteinSW


/* Function:  calculate_dpenv_ProteinSW(mat,dpenv)
 *
 * Descrip:    This function calculates the ProteinSW matrix when in explicit mode, subject to the envelope
 *
 *
 * Arg:          mat [UNKN ] ProteinSW which contains explicit basematrix memory [ProteinSW *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_dpenv_ProteinSW(ProteinSW * mat,DPEnvelope * dpenv);
#define calculate_dpenv_ProteinSW Wise2_calculate_dpenv_ProteinSW


/* Function:  ProteinSW_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ProteinSW *]
 *
 */
ProteinSW * Wise2_ProteinSW_alloc(void);
#define ProteinSW_alloc Wise2_ProteinSW_alloc


/* Function:  free_ProteinSW(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ProteinSW *]
 *
 * Return [UNKN ]  Undocumented return value [ProteinSW *]
 *
 */
ProteinSW * Wise2_free_ProteinSW(ProteinSW * obj);
#define free_ProteinSW Wise2_free_ProteinSW


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
int Wise2_ProteinSW_shatter_access_main(ProteinSW * mat,int i,int j,int state);
#define ProteinSW_shatter_access_main Wise2_ProteinSW_shatter_access_main
int Wise2_ProteinSW_shatter_access_special(ProteinSW * mat,int i,int j,int state);
#define ProteinSW_shatter_access_special Wise2_ProteinSW_shatter_access_special
void * Wise2_thread_loop_ProteinSW(void * ptr);
#define thread_loop_ProteinSW Wise2_thread_loop_ProteinSW
int Wise2_score_only_ProteinSW(ComplexSequence* query,ComplexSequence* target ,CompMat* comp,int gap,int ext);
#define score_only_ProteinSW Wise2_score_only_ProteinSW
ProteinSW * Wise2_allocate_ProteinSW_only(ComplexSequence* query,ComplexSequence* target ,CompMat* comp,int gap,int ext);
#define allocate_ProteinSW_only Wise2_allocate_ProteinSW_only
void Wise2_init_ProteinSW(ProteinSW * mat);
#define init_ProteinSW Wise2_init_ProteinSW
AlnRange * Wise2_AlnRange_build_ProteinSW(ProteinSW * mat,int stopj,int stopspecstate,int * startj,int * startspecstate);
#define AlnRange_build_ProteinSW Wise2_AlnRange_build_ProteinSW
boolean Wise2_read_hidden_ProteinSW(ProteinSW * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out);
#define read_hidden_ProteinSW Wise2_read_hidden_ProteinSW
int Wise2_max_hidden_ProteinSW(ProteinSW * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_hidden_ProteinSW Wise2_max_hidden_ProteinSW
boolean Wise2_read_special_strip_ProteinSW(ProteinSW * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out);
#define read_special_strip_ProteinSW Wise2_read_special_strip_ProteinSW
int Wise2_max_special_strip_ProteinSW(ProteinSW * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_special_strip_ProteinSW Wise2_max_special_strip_ProteinSW
int Wise2_max_matrix_to_special_ProteinSW(ProteinSW * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_matrix_to_special_ProteinSW Wise2_max_matrix_to_special_ProteinSW
void Wise2_calculate_hidden_ProteinSW(ProteinSW * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv);
#define calculate_hidden_ProteinSW Wise2_calculate_hidden_ProteinSW
void Wise2_init_hidden_ProteinSW(ProteinSW * mat,int starti,int startj,int stopi,int stopj);
#define init_hidden_ProteinSW Wise2_init_hidden_ProteinSW
boolean Wise2_full_dc_ProteinSW(ProteinSW * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv);
#define full_dc_ProteinSW Wise2_full_dc_ProteinSW
boolean Wise2_do_dc_single_pass_ProteinSW(ProteinSW * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done);
#define do_dc_single_pass_ProteinSW Wise2_do_dc_single_pass_ProteinSW
void Wise2_push_dc_at_merge_ProteinSW(ProteinSW * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv);
#define push_dc_at_merge_ProteinSW Wise2_push_dc_at_merge_ProteinSW
void Wise2_follow_on_dc_ProteinSW(ProteinSW * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define follow_on_dc_ProteinSW Wise2_follow_on_dc_ProteinSW
void Wise2_run_up_dc_ProteinSW(ProteinSW * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define run_up_dc_ProteinSW Wise2_run_up_dc_ProteinSW
void Wise2_init_dc_ProteinSW(ProteinSW * mat);
#define init_dc_ProteinSW Wise2_init_dc_ProteinSW
int Wise2_start_end_find_end_ProteinSW(ProteinSW * mat,int * endj);
#define start_end_find_end_ProteinSW Wise2_start_end_find_end_ProteinSW
boolean Wise2_dc_optimised_start_end_calc_ProteinSW(ProteinSW *mat,DPEnvelope * dpenv);
#define dc_optimised_start_end_calc_ProteinSW Wise2_dc_optimised_start_end_calc_ProteinSW
void Wise2_init_start_end_linear_ProteinSW(ProteinSW * mat);
#define init_start_end_linear_ProteinSW Wise2_init_start_end_linear_ProteinSW
AlnConvertSet * Wise2_AlnConvertSet_ProteinSW(void);
#define AlnConvertSet_ProteinSW Wise2_AlnConvertSet_ProteinSW
int Wise2_ProteinSW_explicit_access_main(ProteinSW * mat,int i,int j,int state);
#define ProteinSW_explicit_access_main Wise2_ProteinSW_explicit_access_main
int Wise2_ProteinSW_explicit_access_special(ProteinSW * mat,int i,int j,int state);
#define ProteinSW_explicit_access_special Wise2_ProteinSW_explicit_access_special
int Wise2_find_end_ProteinSW(ProteinSW * mat,int * ri,int * rj,int * state,boolean * isspecial,ProteinSW_access_func_holder h);
#define find_end_ProteinSW Wise2_find_end_ProteinSW
void Wise2_ProteinSW_debug_show_matrix(ProteinSW * mat,int starti,int stopi,int startj,int stopj,FILE * ofp);
#define ProteinSW_debug_show_matrix Wise2_ProteinSW_debug_show_matrix
int Wise2_max_calc_ProteinSW(ProteinSW * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,ProteinSW_access_func_holder h);
#define max_calc_ProteinSW Wise2_max_calc_ProteinSW
int Wise2_max_calc_special_ProteinSW(ProteinSW * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,ProteinSW_access_func_holder h);
#define max_calc_special_ProteinSW Wise2_max_calc_special_ProteinSW

#ifdef _cplusplus
}
#endif

#endif
