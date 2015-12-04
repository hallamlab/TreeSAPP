#ifndef DYNAMITEthreestateloopHEADERFILE
#define DYNAMITEthreestateloopHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "threestatemodel.h"


struct Wise2_ThreeStateLoop {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    BaseMatrix * basematrix;     
    ShatterMatrix * shatter;     
    int leni;    
    int lenj;    
    ThreeStateScore* query;  
    ComplexSequence* target;     
    } ;  
/* ThreeStateLoop defined */ 
#ifndef DYNAMITE_DEFINED_ThreeStateLoop
typedef struct Wise2_ThreeStateLoop Wise2_ThreeStateLoop;
#define ThreeStateLoop Wise2_ThreeStateLoop
#define DYNAMITE_DEFINED_ThreeStateLoop
#endif


#ifdef PTHREAD
struct thread_pool_holder_ThreeStateLoop {  
    ThreeStateScore* query; /* Static query data: never free'd */ 
    ComplexSequence* target;/* Target object placeholder */ 
    ProteinDB* targetdb;/* Target database object */ 
    boolean target_init; 
    pthread_mutex_t input_lock;  
    pthread_mutex_t output_lock;     
    Hscore * out;    
    pthread_t * pool;    
    int number_of_threads;   
    boolean search_has_ended;    
    DBSearchImpl * dbsi;     
    } ;  
#endif /* PTHREAD */
struct Wise2_ThreeStateLoop_access_func_holder {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int (*access_main)(ThreeStateLoop*,int,int,int); 
    int (*access_special)(ThreeStateLoop*,int,int,int);  
    } ;  
/* ThreeStateLoop_access_func_holder defined */ 
#ifndef DYNAMITE_DEFINED_ThreeStateLoop_access_func_holder
typedef struct Wise2_ThreeStateLoop_access_func_holder Wise2_ThreeStateLoop_access_func_holder;
#define ThreeStateLoop_access_func_holder Wise2_ThreeStateLoop_access_func_holder
#define DYNAMITE_DEFINED_ThreeStateLoop_access_func_holder
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  PackAln_read_Shatter_ThreeStateLoop(mat)
 *
 * Descrip:    Reads off PackAln from shatter matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [ThreeStateLoop *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Shatter_ThreeStateLoop(ThreeStateLoop * mat);
#define PackAln_read_Shatter_ThreeStateLoop Wise2_PackAln_read_Shatter_ThreeStateLoop


/* Function:  calculate_shatter_ThreeStateLoop(mat,dpenv)
 *
 * Descrip:    This function calculates the ThreeStateLoop matrix when in shatter mode
 *
 *
 * Arg:          mat [UNKN ] (null) [ThreeStateLoop *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_shatter_ThreeStateLoop(ThreeStateLoop * mat,DPEnvelope * dpenv);
#define calculate_shatter_ThreeStateLoop Wise2_calculate_shatter_ThreeStateLoop


/* Function:  search_ThreeStateLoop(dbsi,out,query,targetdb)
 *
 * Descrip:    This function makes a database search of ThreeStateLoop
 *             It uses the dbsi structure to choose which implementation
 *             to use of the database searching. This way at run time you
 *             can switch between single threaded/multi-threaded or hardware
 *
 *
 * Arg:            dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 * Arg:             out [UNKN ] Undocumented argument [Hscore *]
 * Arg:           query [UNKN ] Undocumented argument [ThreeStateScore*]
 * Arg:        targetdb [UNKN ] Undocumented argument [ProteinDB*]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_search_ThreeStateLoop(DBSearchImpl * dbsi,Hscore * out,ThreeStateScore* query,ProteinDB* targetdb );
#define search_ThreeStateLoop Wise2_search_ThreeStateLoop


/* Function:  serial_search_ThreeStateLoop(out,query,targetdb)
 *
 * Descrip:    This function makes a database search of ThreeStateLoop
 *             It is a single processor implementation
 *
 *
 * Arg:             out [UNKN ] Undocumented argument [Hscore *]
 * Arg:           query [UNKN ] Undocumented argument [ThreeStateScore*]
 * Arg:        targetdb [UNKN ] Undocumented argument [ProteinDB*]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_serial_search_ThreeStateLoop(Hscore * out,ThreeStateScore* query,ProteinDB* targetdb );
#define serial_search_ThreeStateLoop Wise2_serial_search_ThreeStateLoop


/* Function:  PackAln_bestmemory_ThreeStateLoop(query,target,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_ThreeStateLoop
 *
 *
 * Arg:         query [UNKN ] query data structure [ThreeStateScore*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:         dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:          dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_bestmemory_ThreeStateLoop(ThreeStateScore* query,ComplexSequence* target ,DPEnvelope * dpenv,DPRunImpl * dpri);
#define PackAln_bestmemory_ThreeStateLoop Wise2_PackAln_bestmemory_ThreeStateLoop


/* Function:  allocate_Expl_ThreeStateLoop(query,target,dpri)
 *
 * Descrip:    This function allocates the ThreeStateLoop structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_ThreeStateLoop_only
 *
 *
 * Arg:         query [UNKN ] query data structure [ThreeStateScore*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:          dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateLoop *]
 *
 */
ThreeStateLoop * Wise2_allocate_Expl_ThreeStateLoop(ThreeStateScore* query,ComplexSequence* target ,DPRunImpl * dpri);
#define allocate_Expl_ThreeStateLoop Wise2_allocate_Expl_ThreeStateLoop


/* Function:  recalculate_PackAln_ThreeStateLoop(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by ThreeStateLoop
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [ThreeStateLoop *]
 *
 */
void Wise2_recalculate_PackAln_ThreeStateLoop(PackAln * pal,ThreeStateLoop * mat);
#define recalculate_PackAln_ThreeStateLoop Wise2_recalculate_PackAln_ThreeStateLoop


/* Function:  allocate_Small_ThreeStateLoop(query,target)
 *
 * Descrip:    This function allocates the ThreeStateLoop structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_ThreeStateLoop_only
 *
 *
 * Arg:         query [UNKN ] query data structure [ThreeStateScore*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateLoop *]
 *
 */
ThreeStateLoop * Wise2_allocate_Small_ThreeStateLoop(ThreeStateScore* query,ComplexSequence* target );
#define allocate_Small_ThreeStateLoop Wise2_allocate_Small_ThreeStateLoop


/* Function:  PackAln_calculate_Small_ThreeStateLoop(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for ThreeStateLoop structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_ThreeStateLoop 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_ThreeStateLoop 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [ThreeStateLoop *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_calculate_Small_ThreeStateLoop(ThreeStateLoop * mat,DPEnvelope * dpenv);
#define PackAln_calculate_Small_ThreeStateLoop Wise2_PackAln_calculate_Small_ThreeStateLoop


/* Function:  AlnRangeSet_calculate_Small_ThreeStateLoop(mat)
 *
 * Descrip:    This function calculates an alignment for ThreeStateLoop structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_ThreeStateLoop 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_ThreeStateLoop
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_ThreeStateLoop 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [ThreeStateLoop *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_calculate_Small_ThreeStateLoop(ThreeStateLoop * mat);
#define AlnRangeSet_calculate_Small_ThreeStateLoop Wise2_AlnRangeSet_calculate_Small_ThreeStateLoop


/* Function:  AlnRangeSet_from_ThreeStateLoop(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for ThreeStateLoop structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_ThreeStateLoop 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_ThreeStateLoop
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [ThreeStateLoop *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_from_ThreeStateLoop(ThreeStateLoop * mat);
#define AlnRangeSet_from_ThreeStateLoop Wise2_AlnRangeSet_from_ThreeStateLoop


/* Function:  convert_PackAln_to_AlnBlock_ThreeStateLoop(pal)
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
AlnBlock * Wise2_convert_PackAln_to_AlnBlock_ThreeStateLoop(PackAln * pal);
#define convert_PackAln_to_AlnBlock_ThreeStateLoop Wise2_convert_PackAln_to_AlnBlock_ThreeStateLoop


/* Function:  PackAln_read_Expl_ThreeStateLoop(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [ThreeStateLoop *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Expl_ThreeStateLoop(ThreeStateLoop * mat);
#define PackAln_read_Expl_ThreeStateLoop Wise2_PackAln_read_Expl_ThreeStateLoop


/* Function:  PackAln_read_generic_ThreeStateLoop(mat,h)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [ThreeStateLoop *]
 * Arg:          h [UNKN ] Undocumented argument [ThreeStateLoop_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_generic_ThreeStateLoop(ThreeStateLoop * mat,ThreeStateLoop_access_func_holder h);
#define PackAln_read_generic_ThreeStateLoop Wise2_PackAln_read_generic_ThreeStateLoop


/* Function:  calculate_ThreeStateLoop(mat)
 *
 * Descrip:    This function calculates the ThreeStateLoop matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_ThreeStateLoop
 *
 *
 * Arg:        mat [UNKN ] ThreeStateLoop which contains explicit basematrix memory [ThreeStateLoop *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_ThreeStateLoop(ThreeStateLoop * mat);
#define calculate_ThreeStateLoop Wise2_calculate_ThreeStateLoop


/* Function:  calculate_dpenv_ThreeStateLoop(mat,dpenv)
 *
 * Descrip:    This function calculates the ThreeStateLoop matrix when in explicit mode, subject to the envelope
 *
 *
 * Arg:          mat [UNKN ] ThreeStateLoop which contains explicit basematrix memory [ThreeStateLoop *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_dpenv_ThreeStateLoop(ThreeStateLoop * mat,DPEnvelope * dpenv);
#define calculate_dpenv_ThreeStateLoop Wise2_calculate_dpenv_ThreeStateLoop


/* Function:  ThreeStateLoop_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateLoop *]
 *
 */
ThreeStateLoop * Wise2_ThreeStateLoop_alloc(void);
#define ThreeStateLoop_alloc Wise2_ThreeStateLoop_alloc


/* Function:  free_ThreeStateLoop(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ThreeStateLoop *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateLoop *]
 *
 */
ThreeStateLoop * Wise2_free_ThreeStateLoop(ThreeStateLoop * obj);
#define free_ThreeStateLoop Wise2_free_ThreeStateLoop


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
int Wise2_ThreeStateLoop_shatter_access_main(ThreeStateLoop * mat,int i,int j,int state);
#define ThreeStateLoop_shatter_access_main Wise2_ThreeStateLoop_shatter_access_main
int Wise2_ThreeStateLoop_shatter_access_special(ThreeStateLoop * mat,int i,int j,int state);
#define ThreeStateLoop_shatter_access_special Wise2_ThreeStateLoop_shatter_access_special
void * Wise2_thread_loop_ThreeStateLoop(void * ptr);
#define thread_loop_ThreeStateLoop Wise2_thread_loop_ThreeStateLoop
int Wise2_score_only_ThreeStateLoop(ThreeStateScore* query,ComplexSequence* target );
#define score_only_ThreeStateLoop Wise2_score_only_ThreeStateLoop
ThreeStateLoop * Wise2_allocate_ThreeStateLoop_only(ThreeStateScore* query,ComplexSequence* target );
#define allocate_ThreeStateLoop_only Wise2_allocate_ThreeStateLoop_only
void Wise2_init_ThreeStateLoop(ThreeStateLoop * mat);
#define init_ThreeStateLoop Wise2_init_ThreeStateLoop
AlnRange * Wise2_AlnRange_build_ThreeStateLoop(ThreeStateLoop * mat,int stopj,int stopspecstate,int * startj,int * startspecstate);
#define AlnRange_build_ThreeStateLoop Wise2_AlnRange_build_ThreeStateLoop
boolean Wise2_read_hidden_ThreeStateLoop(ThreeStateLoop * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out);
#define read_hidden_ThreeStateLoop Wise2_read_hidden_ThreeStateLoop
int Wise2_max_hidden_ThreeStateLoop(ThreeStateLoop * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_hidden_ThreeStateLoop Wise2_max_hidden_ThreeStateLoop
boolean Wise2_read_special_strip_ThreeStateLoop(ThreeStateLoop * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out);
#define read_special_strip_ThreeStateLoop Wise2_read_special_strip_ThreeStateLoop
int Wise2_max_special_strip_ThreeStateLoop(ThreeStateLoop * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_special_strip_ThreeStateLoop Wise2_max_special_strip_ThreeStateLoop
int Wise2_max_matrix_to_special_ThreeStateLoop(ThreeStateLoop * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_matrix_to_special_ThreeStateLoop Wise2_max_matrix_to_special_ThreeStateLoop
void Wise2_calculate_hidden_ThreeStateLoop(ThreeStateLoop * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv);
#define calculate_hidden_ThreeStateLoop Wise2_calculate_hidden_ThreeStateLoop
void Wise2_init_hidden_ThreeStateLoop(ThreeStateLoop * mat,int starti,int startj,int stopi,int stopj);
#define init_hidden_ThreeStateLoop Wise2_init_hidden_ThreeStateLoop
boolean Wise2_full_dc_ThreeStateLoop(ThreeStateLoop * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv);
#define full_dc_ThreeStateLoop Wise2_full_dc_ThreeStateLoop
boolean Wise2_do_dc_single_pass_ThreeStateLoop(ThreeStateLoop * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done);
#define do_dc_single_pass_ThreeStateLoop Wise2_do_dc_single_pass_ThreeStateLoop
void Wise2_push_dc_at_merge_ThreeStateLoop(ThreeStateLoop * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv);
#define push_dc_at_merge_ThreeStateLoop Wise2_push_dc_at_merge_ThreeStateLoop
void Wise2_follow_on_dc_ThreeStateLoop(ThreeStateLoop * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define follow_on_dc_ThreeStateLoop Wise2_follow_on_dc_ThreeStateLoop
void Wise2_run_up_dc_ThreeStateLoop(ThreeStateLoop * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define run_up_dc_ThreeStateLoop Wise2_run_up_dc_ThreeStateLoop
void Wise2_init_dc_ThreeStateLoop(ThreeStateLoop * mat);
#define init_dc_ThreeStateLoop Wise2_init_dc_ThreeStateLoop
int Wise2_start_end_find_end_ThreeStateLoop(ThreeStateLoop * mat,int * endj);
#define start_end_find_end_ThreeStateLoop Wise2_start_end_find_end_ThreeStateLoop
boolean Wise2_dc_optimised_start_end_calc_ThreeStateLoop(ThreeStateLoop *mat,DPEnvelope * dpenv);
#define dc_optimised_start_end_calc_ThreeStateLoop Wise2_dc_optimised_start_end_calc_ThreeStateLoop
void Wise2_init_start_end_linear_ThreeStateLoop(ThreeStateLoop * mat);
#define init_start_end_linear_ThreeStateLoop Wise2_init_start_end_linear_ThreeStateLoop
AlnConvertSet * Wise2_AlnConvertSet_ThreeStateLoop(void);
#define AlnConvertSet_ThreeStateLoop Wise2_AlnConvertSet_ThreeStateLoop
int Wise2_ThreeStateLoop_explicit_access_main(ThreeStateLoop * mat,int i,int j,int state);
#define ThreeStateLoop_explicit_access_main Wise2_ThreeStateLoop_explicit_access_main
int Wise2_ThreeStateLoop_explicit_access_special(ThreeStateLoop * mat,int i,int j,int state);
#define ThreeStateLoop_explicit_access_special Wise2_ThreeStateLoop_explicit_access_special
int Wise2_find_end_ThreeStateLoop(ThreeStateLoop * mat,int * ri,int * rj,int * state,boolean * isspecial,ThreeStateLoop_access_func_holder h);
#define find_end_ThreeStateLoop Wise2_find_end_ThreeStateLoop
void Wise2_ThreeStateLoop_debug_show_matrix(ThreeStateLoop * mat,int starti,int stopi,int startj,int stopj,FILE * ofp);
#define ThreeStateLoop_debug_show_matrix Wise2_ThreeStateLoop_debug_show_matrix
int Wise2_max_calc_ThreeStateLoop(ThreeStateLoop * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,ThreeStateLoop_access_func_holder h);
#define max_calc_ThreeStateLoop Wise2_max_calc_ThreeStateLoop
int Wise2_max_calc_special_ThreeStateLoop(ThreeStateLoop * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,ThreeStateLoop_access_func_holder h);
#define max_calc_special_ThreeStateLoop Wise2_max_calc_special_ThreeStateLoop

#ifdef _cplusplus
}
#endif

#endif
