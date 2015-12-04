#ifndef DYNAMITEclonewisedpHEADERFILE
#define DYNAMITEclonewisedpHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "mapstruct.h"


struct Wise2_CloneWise {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    BaseMatrix * basematrix;     
    ShatterMatrix * shatter;     
    int leni;    
    int lenj;    
    MappedCloneSet * q;  
    MappedCloneSet * t;  
    MappedCloneMatch* match;     
    Score target_skip_start;     
    Score target_skip;   
    Score query_skip_start;  
    Score query_skip;    
    int spread;  
    int target_special_s;    
    } ;  
/* CloneWise defined */ 
#ifndef DYNAMITE_DEFINED_CloneWise
typedef struct Wise2_CloneWise Wise2_CloneWise;
#define CloneWise Wise2_CloneWise
#define DYNAMITE_DEFINED_CloneWise
#endif


#ifdef PTHREAD
struct thread_pool_holder_CloneWise {  
    MappedCloneSet * q; /* Static query data: never free'd */ 
    MappedCloneSet * t; /* Static target data: never free'd */ 
    MappedCloneMatch* match;     
    Score target_skip_start;     
    Score target_skip;   
    Score query_skip_start;  
    Score query_skip;    
    int spread;  
    int target_special_s;    
    pthread_mutex_t input_lock;  
    pthread_mutex_t output_lock;     
    Hscore * out;    
    pthread_t * pool;    
    int number_of_threads;   
    boolean search_has_ended;    
    DBSearchImpl * dbsi;     
    } ;  
#endif /* PTHREAD */
struct Wise2_CloneWise_access_func_holder {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int (*access_main)(CloneWise*,int,int,int);  
    int (*access_special)(CloneWise*,int,int,int);   
    } ;  
/* CloneWise_access_func_holder defined */ 
#ifndef DYNAMITE_DEFINED_CloneWise_access_func_holder
typedef struct Wise2_CloneWise_access_func_holder Wise2_CloneWise_access_func_holder;
#define CloneWise_access_func_holder Wise2_CloneWise_access_func_holder
#define DYNAMITE_DEFINED_CloneWise_access_func_holder
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  PackAln_read_Shatter_CloneWise(mat)
 *
 * Descrip:    Reads off PackAln from shatter matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [CloneWise *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Shatter_CloneWise(CloneWise * mat);
#define PackAln_read_Shatter_CloneWise Wise2_PackAln_read_Shatter_CloneWise


/* Function:  calculate_shatter_CloneWise(mat,dpenv)
 *
 * Descrip:    This function calculates the CloneWise matrix when in shatter mode
 *
 *
 * Arg:          mat [UNKN ] (null) [CloneWise *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_shatter_CloneWise(CloneWise * mat,DPEnvelope * dpenv);
#define calculate_shatter_CloneWise Wise2_calculate_shatter_CloneWise


/* Function:  search_CloneWise(dbsi,out,q,t,match,target_skip_start,target_skip,query_skip_start,query_skip,spread,target_special_s)
 *
 * Descrip:    This function makes a database search of CloneWise
 *             It uses the dbsi structure to choose which implementation
 *             to use of the database searching. This way at run time you
 *             can switch between single threaded/multi-threaded or hardware
 *
 *
 * Arg:                     dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 * Arg:                      out [UNKN ] Undocumented argument [Hscore *]
 * Arg:                        q [UNKN ] Undocumented argument [MappedCloneSet *]
 * Arg:                        t [UNKN ] Undocumented argument [MappedCloneSet *]
 * Arg:                    match [UNKN ] Undocumented argument [MappedCloneMatch*]
 * Arg:        target_skip_start [UNKN ] Undocumented argument [Score]
 * Arg:              target_skip [UNKN ] Undocumented argument [Score]
 * Arg:         query_skip_start [UNKN ] Undocumented argument [Score]
 * Arg:               query_skip [UNKN ] Undocumented argument [Score]
 * Arg:                   spread [UNKN ] Undocumented argument [int]
 * Arg:         target_special_s [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_search_CloneWise(DBSearchImpl * dbsi,Hscore * out,MappedCloneSet * q,MappedCloneSet * t ,MappedCloneMatch* match,Score target_skip_start,Score target_skip,Score query_skip_start,Score query_skip,int spread,int target_special_s);
#define search_CloneWise Wise2_search_CloneWise


/* Function:  serial_search_CloneWise(out,q,t,match,target_skip_start,target_skip,query_skip_start,query_skip,spread,target_special_s)
 *
 * Descrip:    This function makes a database search of CloneWise
 *             It is a single processor implementation
 *
 *
 * Arg:                      out [UNKN ] Undocumented argument [Hscore *]
 * Arg:                        q [UNKN ] Undocumented argument [MappedCloneSet *]
 * Arg:                        t [UNKN ] Undocumented argument [MappedCloneSet *]
 * Arg:                    match [UNKN ] Undocumented argument [MappedCloneMatch*]
 * Arg:        target_skip_start [UNKN ] Undocumented argument [Score]
 * Arg:              target_skip [UNKN ] Undocumented argument [Score]
 * Arg:         query_skip_start [UNKN ] Undocumented argument [Score]
 * Arg:               query_skip [UNKN ] Undocumented argument [Score]
 * Arg:                   spread [UNKN ] Undocumented argument [int]
 * Arg:         target_special_s [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_serial_search_CloneWise(Hscore * out,MappedCloneSet * q,MappedCloneSet * t ,MappedCloneMatch* match,Score target_skip_start,Score target_skip,Score query_skip_start,Score query_skip,int spread,int target_special_s);
#define serial_search_CloneWise Wise2_serial_search_CloneWise


/* Function:  PackAln_bestmemory_CloneWise(q,t,match,target_skip_start,target_skip,query_skip_start,query_skip,spread,target_special_s,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_CloneWise
 *
 *
 * Arg:                        q [UNKN ] query data structure [MappedCloneSet *]
 * Arg:                        t [UNKN ] target data structure [MappedCloneSet *]
 * Arg:                    match [UNKN ] Resource [MappedCloneMatch*]
 * Arg:        target_skip_start [UNKN ] Resource [Score]
 * Arg:              target_skip [UNKN ] Resource [Score]
 * Arg:         query_skip_start [UNKN ] Resource [Score]
 * Arg:               query_skip [UNKN ] Resource [Score]
 * Arg:                   spread [UNKN ] Resource [int]
 * Arg:         target_special_s [UNKN ] Resource [int]
 * Arg:                    dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:                     dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_bestmemory_CloneWise(MappedCloneSet * q,MappedCloneSet * t ,MappedCloneMatch* match,Score target_skip_start,Score target_skip,Score query_skip_start,Score query_skip,int spread,int target_special_s,DPEnvelope * dpenv,DPRunImpl * dpri);
#define PackAln_bestmemory_CloneWise Wise2_PackAln_bestmemory_CloneWise


/* Function:  allocate_Expl_CloneWise(q,t,match,target_skip_start,target_skip,query_skip_start,query_skip,spread,target_special_s,dpri)
 *
 * Descrip:    This function allocates the CloneWise structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_CloneWise_only
 *
 *
 * Arg:                        q [UNKN ] query data structure [MappedCloneSet *]
 * Arg:                        t [UNKN ] target data structure [MappedCloneSet *]
 * Arg:                    match [UNKN ] Resource [MappedCloneMatch*]
 * Arg:        target_skip_start [UNKN ] Resource [Score]
 * Arg:              target_skip [UNKN ] Resource [Score]
 * Arg:         query_skip_start [UNKN ] Resource [Score]
 * Arg:               query_skip [UNKN ] Resource [Score]
 * Arg:                   spread [UNKN ] Resource [int]
 * Arg:         target_special_s [UNKN ] Resource [int]
 * Arg:                     dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [CloneWise *]
 *
 */
CloneWise * Wise2_allocate_Expl_CloneWise(MappedCloneSet * q,MappedCloneSet * t ,MappedCloneMatch* match,Score target_skip_start,Score target_skip,Score query_skip_start,Score query_skip,int spread,int target_special_s,DPRunImpl * dpri);
#define allocate_Expl_CloneWise Wise2_allocate_Expl_CloneWise


/* Function:  recalculate_PackAln_CloneWise(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by CloneWise
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [CloneWise *]
 *
 */
void Wise2_recalculate_PackAln_CloneWise(PackAln * pal,CloneWise * mat);
#define recalculate_PackAln_CloneWise Wise2_recalculate_PackAln_CloneWise


/* Function:  allocate_Small_CloneWise(q,t,match,target_skip_start,target_skip,query_skip_start,query_skip,spread,target_special_s)
 *
 * Descrip:    This function allocates the CloneWise structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_CloneWise_only
 *
 *
 * Arg:                        q [UNKN ] query data structure [MappedCloneSet *]
 * Arg:                        t [UNKN ] target data structure [MappedCloneSet *]
 * Arg:                    match [UNKN ] Resource [MappedCloneMatch*]
 * Arg:        target_skip_start [UNKN ] Resource [Score]
 * Arg:              target_skip [UNKN ] Resource [Score]
 * Arg:         query_skip_start [UNKN ] Resource [Score]
 * Arg:               query_skip [UNKN ] Resource [Score]
 * Arg:                   spread [UNKN ] Resource [int]
 * Arg:         target_special_s [UNKN ] Resource [int]
 *
 * Return [UNKN ]  Undocumented return value [CloneWise *]
 *
 */
CloneWise * Wise2_allocate_Small_CloneWise(MappedCloneSet * q,MappedCloneSet * t ,MappedCloneMatch* match,Score target_skip_start,Score target_skip,Score query_skip_start,Score query_skip,int spread,int target_special_s);
#define allocate_Small_CloneWise Wise2_allocate_Small_CloneWise


/* Function:  PackAln_calculate_Small_CloneWise(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for CloneWise structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_CloneWise 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_CloneWise 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [CloneWise *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_calculate_Small_CloneWise(CloneWise * mat,DPEnvelope * dpenv);
#define PackAln_calculate_Small_CloneWise Wise2_PackAln_calculate_Small_CloneWise


/* Function:  AlnRangeSet_calculate_Small_CloneWise(mat)
 *
 * Descrip:    This function calculates an alignment for CloneWise structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_CloneWise 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_CloneWise
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_CloneWise 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [CloneWise *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_calculate_Small_CloneWise(CloneWise * mat);
#define AlnRangeSet_calculate_Small_CloneWise Wise2_AlnRangeSet_calculate_Small_CloneWise


/* Function:  AlnRangeSet_from_CloneWise(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for CloneWise structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_CloneWise 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_CloneWise
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [CloneWise *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_from_CloneWise(CloneWise * mat);
#define AlnRangeSet_from_CloneWise Wise2_AlnRangeSet_from_CloneWise


/* Function:  convert_PackAln_to_AlnBlock_CloneWise(pal)
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
AlnBlock * Wise2_convert_PackAln_to_AlnBlock_CloneWise(PackAln * pal);
#define convert_PackAln_to_AlnBlock_CloneWise Wise2_convert_PackAln_to_AlnBlock_CloneWise


/* Function:  PackAln_read_Expl_CloneWise(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [CloneWise *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Expl_CloneWise(CloneWise * mat);
#define PackAln_read_Expl_CloneWise Wise2_PackAln_read_Expl_CloneWise


/* Function:  PackAln_read_generic_CloneWise(mat,h)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [CloneWise *]
 * Arg:          h [UNKN ] Undocumented argument [CloneWise_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_generic_CloneWise(CloneWise * mat,CloneWise_access_func_holder h);
#define PackAln_read_generic_CloneWise Wise2_PackAln_read_generic_CloneWise


/* Function:  calculate_CloneWise(mat)
 *
 * Descrip:    This function calculates the CloneWise matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_CloneWise
 *
 *
 * Arg:        mat [UNKN ] CloneWise which contains explicit basematrix memory [CloneWise *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_CloneWise(CloneWise * mat);
#define calculate_CloneWise Wise2_calculate_CloneWise


/* Function:  calculate_dpenv_CloneWise(mat,dpenv)
 *
 * Descrip:    This function calculates the CloneWise matrix when in explicit mode, subject to the envelope
 *
 *
 * Arg:          mat [UNKN ] CloneWise which contains explicit basematrix memory [CloneWise *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_dpenv_CloneWise(CloneWise * mat,DPEnvelope * dpenv);
#define calculate_dpenv_CloneWise Wise2_calculate_dpenv_CloneWise


/* Function:  CloneWise_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [CloneWise *]
 *
 */
CloneWise * Wise2_CloneWise_alloc(void);
#define CloneWise_alloc Wise2_CloneWise_alloc


/* Function:  free_CloneWise(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [CloneWise *]
 *
 * Return [UNKN ]  Undocumented return value [CloneWise *]
 *
 */
CloneWise * Wise2_free_CloneWise(CloneWise * obj);
#define free_CloneWise Wise2_free_CloneWise


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
int Wise2_CloneWise_shatter_access_main(CloneWise * mat,int i,int j,int state);
#define CloneWise_shatter_access_main Wise2_CloneWise_shatter_access_main
int Wise2_CloneWise_shatter_access_special(CloneWise * mat,int i,int j,int state);
#define CloneWise_shatter_access_special Wise2_CloneWise_shatter_access_special
void * Wise2_thread_loop_CloneWise(void * ptr);
#define thread_loop_CloneWise Wise2_thread_loop_CloneWise
int Wise2_score_only_CloneWise(MappedCloneSet * q,MappedCloneSet * t ,MappedCloneMatch* match,Score target_skip_start,Score target_skip,Score query_skip_start,Score query_skip,int spread,int target_special_s);
#define score_only_CloneWise Wise2_score_only_CloneWise
CloneWise * Wise2_allocate_CloneWise_only(MappedCloneSet * q,MappedCloneSet * t ,MappedCloneMatch* match,Score target_skip_start,Score target_skip,Score query_skip_start,Score query_skip,int spread,int target_special_s);
#define allocate_CloneWise_only Wise2_allocate_CloneWise_only
void Wise2_init_CloneWise(CloneWise * mat);
#define init_CloneWise Wise2_init_CloneWise
AlnRange * Wise2_AlnRange_build_CloneWise(CloneWise * mat,int stopj,int stopspecstate,int * startj,int * startspecstate);
#define AlnRange_build_CloneWise Wise2_AlnRange_build_CloneWise
boolean Wise2_read_hidden_CloneWise(CloneWise * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out);
#define read_hidden_CloneWise Wise2_read_hidden_CloneWise
int Wise2_max_hidden_CloneWise(CloneWise * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_hidden_CloneWise Wise2_max_hidden_CloneWise
boolean Wise2_read_special_strip_CloneWise(CloneWise * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out);
#define read_special_strip_CloneWise Wise2_read_special_strip_CloneWise
int Wise2_max_special_strip_CloneWise(CloneWise * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_special_strip_CloneWise Wise2_max_special_strip_CloneWise
int Wise2_max_matrix_to_special_CloneWise(CloneWise * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_matrix_to_special_CloneWise Wise2_max_matrix_to_special_CloneWise
void Wise2_calculate_hidden_CloneWise(CloneWise * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv);
#define calculate_hidden_CloneWise Wise2_calculate_hidden_CloneWise
void Wise2_init_hidden_CloneWise(CloneWise * mat,int starti,int startj,int stopi,int stopj);
#define init_hidden_CloneWise Wise2_init_hidden_CloneWise
boolean Wise2_full_dc_CloneWise(CloneWise * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv);
#define full_dc_CloneWise Wise2_full_dc_CloneWise
boolean Wise2_do_dc_single_pass_CloneWise(CloneWise * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done);
#define do_dc_single_pass_CloneWise Wise2_do_dc_single_pass_CloneWise
void Wise2_push_dc_at_merge_CloneWise(CloneWise * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv);
#define push_dc_at_merge_CloneWise Wise2_push_dc_at_merge_CloneWise
void Wise2_follow_on_dc_CloneWise(CloneWise * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define follow_on_dc_CloneWise Wise2_follow_on_dc_CloneWise
void Wise2_run_up_dc_CloneWise(CloneWise * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define run_up_dc_CloneWise Wise2_run_up_dc_CloneWise
void Wise2_init_dc_CloneWise(CloneWise * mat);
#define init_dc_CloneWise Wise2_init_dc_CloneWise
int Wise2_start_end_find_end_CloneWise(CloneWise * mat,int * endj);
#define start_end_find_end_CloneWise Wise2_start_end_find_end_CloneWise
boolean Wise2_dc_optimised_start_end_calc_CloneWise(CloneWise *mat,DPEnvelope * dpenv);
#define dc_optimised_start_end_calc_CloneWise Wise2_dc_optimised_start_end_calc_CloneWise
void Wise2_init_start_end_linear_CloneWise(CloneWise * mat);
#define init_start_end_linear_CloneWise Wise2_init_start_end_linear_CloneWise
AlnConvertSet * Wise2_AlnConvertSet_CloneWise(void);
#define AlnConvertSet_CloneWise Wise2_AlnConvertSet_CloneWise
int Wise2_CloneWise_explicit_access_main(CloneWise * mat,int i,int j,int state);
#define CloneWise_explicit_access_main Wise2_CloneWise_explicit_access_main
int Wise2_CloneWise_explicit_access_special(CloneWise * mat,int i,int j,int state);
#define CloneWise_explicit_access_special Wise2_CloneWise_explicit_access_special
int Wise2_find_end_CloneWise(CloneWise * mat,int * ri,int * rj,int * state,boolean * isspecial,CloneWise_access_func_holder h);
#define find_end_CloneWise Wise2_find_end_CloneWise
void Wise2_CloneWise_debug_show_matrix(CloneWise * mat,int starti,int stopi,int startj,int stopj,FILE * ofp);
#define CloneWise_debug_show_matrix Wise2_CloneWise_debug_show_matrix
int Wise2_max_calc_CloneWise(CloneWise * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,CloneWise_access_func_holder h);
#define max_calc_CloneWise Wise2_max_calc_CloneWise
int Wise2_max_calc_special_CloneWise(CloneWise * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,CloneWise_access_func_holder h);
#define max_calc_special_CloneWise Wise2_max_calc_special_CloneWise

#ifdef _cplusplus
}
#endif

#endif
