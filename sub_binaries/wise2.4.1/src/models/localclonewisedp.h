#ifndef DYNAMITElocalclonewisedpHEADERFILE
#define DYNAMITElocalclonewisedpHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "mapstruct.h"


struct Wise2_LocalCloneWise {  
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
/* LocalCloneWise defined */ 
#ifndef DYNAMITE_DEFINED_LocalCloneWise
typedef struct Wise2_LocalCloneWise Wise2_LocalCloneWise;
#define LocalCloneWise Wise2_LocalCloneWise
#define DYNAMITE_DEFINED_LocalCloneWise
#endif


#ifdef PTHREAD
struct thread_pool_holder_LocalCloneWise {  
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
struct Wise2_LocalCloneWise_access_func_holder {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int (*access_main)(LocalCloneWise*,int,int,int); 
    int (*access_special)(LocalCloneWise*,int,int,int);  
    } ;  
/* LocalCloneWise_access_func_holder defined */ 
#ifndef DYNAMITE_DEFINED_LocalCloneWise_access_func_holder
typedef struct Wise2_LocalCloneWise_access_func_holder Wise2_LocalCloneWise_access_func_holder;
#define LocalCloneWise_access_func_holder Wise2_LocalCloneWise_access_func_holder
#define DYNAMITE_DEFINED_LocalCloneWise_access_func_holder
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  PackAln_read_Shatter_LocalCloneWise(mat)
 *
 * Descrip:    Reads off PackAln from shatter matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [LocalCloneWise *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Shatter_LocalCloneWise(LocalCloneWise * mat);
#define PackAln_read_Shatter_LocalCloneWise Wise2_PackAln_read_Shatter_LocalCloneWise


/* Function:  calculate_shatter_LocalCloneWise(mat,dpenv)
 *
 * Descrip:    This function calculates the LocalCloneWise matrix when in shatter mode
 *
 *
 * Arg:          mat [UNKN ] (null) [LocalCloneWise *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_shatter_LocalCloneWise(LocalCloneWise * mat,DPEnvelope * dpenv);
#define calculate_shatter_LocalCloneWise Wise2_calculate_shatter_LocalCloneWise


/* Function:  search_LocalCloneWise(dbsi,out,q,t,match,target_skip_start,target_skip,query_skip_start,query_skip,spread,target_special_s)
 *
 * Descrip:    This function makes a database search of LocalCloneWise
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
Search_Return_Type Wise2_search_LocalCloneWise(DBSearchImpl * dbsi,Hscore * out,MappedCloneSet * q,MappedCloneSet * t ,MappedCloneMatch* match,Score target_skip_start,Score target_skip,Score query_skip_start,Score query_skip,int spread,int target_special_s);
#define search_LocalCloneWise Wise2_search_LocalCloneWise


/* Function:  serial_search_LocalCloneWise(out,q,t,match,target_skip_start,target_skip,query_skip_start,query_skip,spread,target_special_s)
 *
 * Descrip:    This function makes a database search of LocalCloneWise
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
Search_Return_Type Wise2_serial_search_LocalCloneWise(Hscore * out,MappedCloneSet * q,MappedCloneSet * t ,MappedCloneMatch* match,Score target_skip_start,Score target_skip,Score query_skip_start,Score query_skip,int spread,int target_special_s);
#define serial_search_LocalCloneWise Wise2_serial_search_LocalCloneWise


/* Function:  PackAln_bestmemory_LocalCloneWise(q,t,match,target_skip_start,target_skip,query_skip_start,query_skip,spread,target_special_s,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_LocalCloneWise
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
PackAln * Wise2_PackAln_bestmemory_LocalCloneWise(MappedCloneSet * q,MappedCloneSet * t ,MappedCloneMatch* match,Score target_skip_start,Score target_skip,Score query_skip_start,Score query_skip,int spread,int target_special_s,DPEnvelope * dpenv,DPRunImpl * dpri);
#define PackAln_bestmemory_LocalCloneWise Wise2_PackAln_bestmemory_LocalCloneWise


/* Function:  allocate_Expl_LocalCloneWise(q,t,match,target_skip_start,target_skip,query_skip_start,query_skip,spread,target_special_s,dpri)
 *
 * Descrip:    This function allocates the LocalCloneWise structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_LocalCloneWise_only
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
 * Return [UNKN ]  Undocumented return value [LocalCloneWise *]
 *
 */
LocalCloneWise * Wise2_allocate_Expl_LocalCloneWise(MappedCloneSet * q,MappedCloneSet * t ,MappedCloneMatch* match,Score target_skip_start,Score target_skip,Score query_skip_start,Score query_skip,int spread,int target_special_s,DPRunImpl * dpri);
#define allocate_Expl_LocalCloneWise Wise2_allocate_Expl_LocalCloneWise


/* Function:  recalculate_PackAln_LocalCloneWise(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by LocalCloneWise
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [LocalCloneWise *]
 *
 */
void Wise2_recalculate_PackAln_LocalCloneWise(PackAln * pal,LocalCloneWise * mat);
#define recalculate_PackAln_LocalCloneWise Wise2_recalculate_PackAln_LocalCloneWise


/* Function:  allocate_Small_LocalCloneWise(q,t,match,target_skip_start,target_skip,query_skip_start,query_skip,spread,target_special_s)
 *
 * Descrip:    This function allocates the LocalCloneWise structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_LocalCloneWise_only
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
 * Return [UNKN ]  Undocumented return value [LocalCloneWise *]
 *
 */
LocalCloneWise * Wise2_allocate_Small_LocalCloneWise(MappedCloneSet * q,MappedCloneSet * t ,MappedCloneMatch* match,Score target_skip_start,Score target_skip,Score query_skip_start,Score query_skip,int spread,int target_special_s);
#define allocate_Small_LocalCloneWise Wise2_allocate_Small_LocalCloneWise


/* Function:  PackAln_calculate_Small_LocalCloneWise(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for LocalCloneWise structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_LocalCloneWise 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_LocalCloneWise 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [LocalCloneWise *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_calculate_Small_LocalCloneWise(LocalCloneWise * mat,DPEnvelope * dpenv);
#define PackAln_calculate_Small_LocalCloneWise Wise2_PackAln_calculate_Small_LocalCloneWise


/* Function:  AlnRangeSet_calculate_Small_LocalCloneWise(mat)
 *
 * Descrip:    This function calculates an alignment for LocalCloneWise structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_LocalCloneWise 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_LocalCloneWise
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_LocalCloneWise 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [LocalCloneWise *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_calculate_Small_LocalCloneWise(LocalCloneWise * mat);
#define AlnRangeSet_calculate_Small_LocalCloneWise Wise2_AlnRangeSet_calculate_Small_LocalCloneWise


/* Function:  AlnRangeSet_from_LocalCloneWise(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for LocalCloneWise structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_LocalCloneWise 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_LocalCloneWise
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [LocalCloneWise *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_from_LocalCloneWise(LocalCloneWise * mat);
#define AlnRangeSet_from_LocalCloneWise Wise2_AlnRangeSet_from_LocalCloneWise


/* Function:  convert_PackAln_to_AlnBlock_LocalCloneWise(pal)
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
AlnBlock * Wise2_convert_PackAln_to_AlnBlock_LocalCloneWise(PackAln * pal);
#define convert_PackAln_to_AlnBlock_LocalCloneWise Wise2_convert_PackAln_to_AlnBlock_LocalCloneWise


/* Function:  PackAln_read_Expl_LocalCloneWise(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [LocalCloneWise *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Expl_LocalCloneWise(LocalCloneWise * mat);
#define PackAln_read_Expl_LocalCloneWise Wise2_PackAln_read_Expl_LocalCloneWise


/* Function:  PackAln_read_generic_LocalCloneWise(mat,h)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [LocalCloneWise *]
 * Arg:          h [UNKN ] Undocumented argument [LocalCloneWise_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_generic_LocalCloneWise(LocalCloneWise * mat,LocalCloneWise_access_func_holder h);
#define PackAln_read_generic_LocalCloneWise Wise2_PackAln_read_generic_LocalCloneWise


/* Function:  calculate_LocalCloneWise(mat)
 *
 * Descrip:    This function calculates the LocalCloneWise matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_LocalCloneWise
 *
 *
 * Arg:        mat [UNKN ] LocalCloneWise which contains explicit basematrix memory [LocalCloneWise *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_LocalCloneWise(LocalCloneWise * mat);
#define calculate_LocalCloneWise Wise2_calculate_LocalCloneWise


/* Function:  calculate_dpenv_LocalCloneWise(mat,dpenv)
 *
 * Descrip:    This function calculates the LocalCloneWise matrix when in explicit mode, subject to the envelope
 *
 *
 * Arg:          mat [UNKN ] LocalCloneWise which contains explicit basematrix memory [LocalCloneWise *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_dpenv_LocalCloneWise(LocalCloneWise * mat,DPEnvelope * dpenv);
#define calculate_dpenv_LocalCloneWise Wise2_calculate_dpenv_LocalCloneWise


/* Function:  LocalCloneWise_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LocalCloneWise *]
 *
 */
LocalCloneWise * Wise2_LocalCloneWise_alloc(void);
#define LocalCloneWise_alloc Wise2_LocalCloneWise_alloc


/* Function:  free_LocalCloneWise(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [LocalCloneWise *]
 *
 * Return [UNKN ]  Undocumented return value [LocalCloneWise *]
 *
 */
LocalCloneWise * Wise2_free_LocalCloneWise(LocalCloneWise * obj);
#define free_LocalCloneWise Wise2_free_LocalCloneWise


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
int Wise2_LocalCloneWise_shatter_access_main(LocalCloneWise * mat,int i,int j,int state);
#define LocalCloneWise_shatter_access_main Wise2_LocalCloneWise_shatter_access_main
int Wise2_LocalCloneWise_shatter_access_special(LocalCloneWise * mat,int i,int j,int state);
#define LocalCloneWise_shatter_access_special Wise2_LocalCloneWise_shatter_access_special
void * Wise2_thread_loop_LocalCloneWise(void * ptr);
#define thread_loop_LocalCloneWise Wise2_thread_loop_LocalCloneWise
int Wise2_score_only_LocalCloneWise(MappedCloneSet * q,MappedCloneSet * t ,MappedCloneMatch* match,Score target_skip_start,Score target_skip,Score query_skip_start,Score query_skip,int spread,int target_special_s);
#define score_only_LocalCloneWise Wise2_score_only_LocalCloneWise
LocalCloneWise * Wise2_allocate_LocalCloneWise_only(MappedCloneSet * q,MappedCloneSet * t ,MappedCloneMatch* match,Score target_skip_start,Score target_skip,Score query_skip_start,Score query_skip,int spread,int target_special_s);
#define allocate_LocalCloneWise_only Wise2_allocate_LocalCloneWise_only
void Wise2_init_LocalCloneWise(LocalCloneWise * mat);
#define init_LocalCloneWise Wise2_init_LocalCloneWise
AlnRange * Wise2_AlnRange_build_LocalCloneWise(LocalCloneWise * mat,int stopj,int stopspecstate,int * startj,int * startspecstate);
#define AlnRange_build_LocalCloneWise Wise2_AlnRange_build_LocalCloneWise
boolean Wise2_read_hidden_LocalCloneWise(LocalCloneWise * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out);
#define read_hidden_LocalCloneWise Wise2_read_hidden_LocalCloneWise
int Wise2_max_hidden_LocalCloneWise(LocalCloneWise * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_hidden_LocalCloneWise Wise2_max_hidden_LocalCloneWise
boolean Wise2_read_special_strip_LocalCloneWise(LocalCloneWise * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out);
#define read_special_strip_LocalCloneWise Wise2_read_special_strip_LocalCloneWise
int Wise2_max_special_strip_LocalCloneWise(LocalCloneWise * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_special_strip_LocalCloneWise Wise2_max_special_strip_LocalCloneWise
int Wise2_max_matrix_to_special_LocalCloneWise(LocalCloneWise * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_matrix_to_special_LocalCloneWise Wise2_max_matrix_to_special_LocalCloneWise
void Wise2_calculate_hidden_LocalCloneWise(LocalCloneWise * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv);
#define calculate_hidden_LocalCloneWise Wise2_calculate_hidden_LocalCloneWise
void Wise2_init_hidden_LocalCloneWise(LocalCloneWise * mat,int starti,int startj,int stopi,int stopj);
#define init_hidden_LocalCloneWise Wise2_init_hidden_LocalCloneWise
boolean Wise2_full_dc_LocalCloneWise(LocalCloneWise * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv);
#define full_dc_LocalCloneWise Wise2_full_dc_LocalCloneWise
boolean Wise2_do_dc_single_pass_LocalCloneWise(LocalCloneWise * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done);
#define do_dc_single_pass_LocalCloneWise Wise2_do_dc_single_pass_LocalCloneWise
void Wise2_push_dc_at_merge_LocalCloneWise(LocalCloneWise * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv);
#define push_dc_at_merge_LocalCloneWise Wise2_push_dc_at_merge_LocalCloneWise
void Wise2_follow_on_dc_LocalCloneWise(LocalCloneWise * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define follow_on_dc_LocalCloneWise Wise2_follow_on_dc_LocalCloneWise
void Wise2_run_up_dc_LocalCloneWise(LocalCloneWise * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define run_up_dc_LocalCloneWise Wise2_run_up_dc_LocalCloneWise
void Wise2_init_dc_LocalCloneWise(LocalCloneWise * mat);
#define init_dc_LocalCloneWise Wise2_init_dc_LocalCloneWise
int Wise2_start_end_find_end_LocalCloneWise(LocalCloneWise * mat,int * endj);
#define start_end_find_end_LocalCloneWise Wise2_start_end_find_end_LocalCloneWise
boolean Wise2_dc_optimised_start_end_calc_LocalCloneWise(LocalCloneWise *mat,DPEnvelope * dpenv);
#define dc_optimised_start_end_calc_LocalCloneWise Wise2_dc_optimised_start_end_calc_LocalCloneWise
void Wise2_init_start_end_linear_LocalCloneWise(LocalCloneWise * mat);
#define init_start_end_linear_LocalCloneWise Wise2_init_start_end_linear_LocalCloneWise
AlnConvertSet * Wise2_AlnConvertSet_LocalCloneWise(void);
#define AlnConvertSet_LocalCloneWise Wise2_AlnConvertSet_LocalCloneWise
int Wise2_LocalCloneWise_explicit_access_main(LocalCloneWise * mat,int i,int j,int state);
#define LocalCloneWise_explicit_access_main Wise2_LocalCloneWise_explicit_access_main
int Wise2_LocalCloneWise_explicit_access_special(LocalCloneWise * mat,int i,int j,int state);
#define LocalCloneWise_explicit_access_special Wise2_LocalCloneWise_explicit_access_special
int Wise2_find_end_LocalCloneWise(LocalCloneWise * mat,int * ri,int * rj,int * state,boolean * isspecial,LocalCloneWise_access_func_holder h);
#define find_end_LocalCloneWise Wise2_find_end_LocalCloneWise
void Wise2_LocalCloneWise_debug_show_matrix(LocalCloneWise * mat,int starti,int stopi,int startj,int stopj,FILE * ofp);
#define LocalCloneWise_debug_show_matrix Wise2_LocalCloneWise_debug_show_matrix
int Wise2_max_calc_LocalCloneWise(LocalCloneWise * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,LocalCloneWise_access_func_holder h);
#define max_calc_LocalCloneWise Wise2_max_calc_LocalCloneWise
int Wise2_max_calc_special_LocalCloneWise(LocalCloneWise * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,LocalCloneWise_access_func_holder h);
#define max_calc_special_LocalCloneWise Wise2_max_calc_special_LocalCloneWise

#ifdef _cplusplus
}
#endif

#endif
