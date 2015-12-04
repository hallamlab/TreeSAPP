#ifndef DYNAMITEfivestateHEADERFILE
#define DYNAMITEfivestateHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "fivestatemodel.h"


struct Wise2_FiveStateProtein {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    BaseMatrix * basematrix;     
    ShatterMatrix * shatter;     
    int leni;    
    int lenj;    
    FiveStateScore* query;   
    ComplexSequence* target;     
    } ;  
/* FiveStateProtein defined */ 
#ifndef DYNAMITE_DEFINED_FiveStateProtein
typedef struct Wise2_FiveStateProtein Wise2_FiveStateProtein;
#define FiveStateProtein Wise2_FiveStateProtein
#define DYNAMITE_DEFINED_FiveStateProtein
#endif


#ifdef PTHREAD
struct thread_pool_holder_FiveStateProtein {  
    FiveStateScore* query;  /* Static query data: never free'd */ 
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
struct Wise2_FiveStateProtein_access_func_holder {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int (*access_main)(FiveStateProtein*,int,int,int);   
    int (*access_special)(FiveStateProtein*,int,int,int);    
    } ;  
/* FiveStateProtein_access_func_holder defined */ 
#ifndef DYNAMITE_DEFINED_FiveStateProtein_access_func_holder
typedef struct Wise2_FiveStateProtein_access_func_holder Wise2_FiveStateProtein_access_func_holder;
#define FiveStateProtein_access_func_holder Wise2_FiveStateProtein_access_func_holder
#define DYNAMITE_DEFINED_FiveStateProtein_access_func_holder
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  PackAln_read_Shatter_FiveStateProtein(mat)
 *
 * Descrip:    Reads off PackAln from shatter matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [FiveStateProtein *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Shatter_FiveStateProtein(FiveStateProtein * mat);
#define PackAln_read_Shatter_FiveStateProtein Wise2_PackAln_read_Shatter_FiveStateProtein


/* Function:  calculate_shatter_FiveStateProtein(mat,dpenv)
 *
 * Descrip:    This function calculates the FiveStateProtein matrix when in shatter mode
 *
 *
 * Arg:          mat [UNKN ] (null) [FiveStateProtein *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_shatter_FiveStateProtein(FiveStateProtein * mat,DPEnvelope * dpenv);
#define calculate_shatter_FiveStateProtein Wise2_calculate_shatter_FiveStateProtein


/* Function:  search_FiveStateProtein(dbsi,out,query,targetdb)
 *
 * Descrip:    This function makes a database search of FiveStateProtein
 *             It uses the dbsi structure to choose which implementation
 *             to use of the database searching. This way at run time you
 *             can switch between single threaded/multi-threaded or hardware
 *
 *
 * Arg:            dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 * Arg:             out [UNKN ] Undocumented argument [Hscore *]
 * Arg:           query [UNKN ] Undocumented argument [FiveStateScore*]
 * Arg:        targetdb [UNKN ] Undocumented argument [ProteinDB*]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_search_FiveStateProtein(DBSearchImpl * dbsi,Hscore * out,FiveStateScore* query,ProteinDB* targetdb );
#define search_FiveStateProtein Wise2_search_FiveStateProtein


/* Function:  serial_search_FiveStateProtein(out,query,targetdb)
 *
 * Descrip:    This function makes a database search of FiveStateProtein
 *             It is a single processor implementation
 *
 *
 * Arg:             out [UNKN ] Undocumented argument [Hscore *]
 * Arg:           query [UNKN ] Undocumented argument [FiveStateScore*]
 * Arg:        targetdb [UNKN ] Undocumented argument [ProteinDB*]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_serial_search_FiveStateProtein(Hscore * out,FiveStateScore* query,ProteinDB* targetdb );
#define serial_search_FiveStateProtein Wise2_serial_search_FiveStateProtein


/* Function:  PackAln_bestmemory_FiveStateProtein(query,target,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_FiveStateProtein
 *
 *
 * Arg:         query [UNKN ] query data structure [FiveStateScore*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:         dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:          dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_bestmemory_FiveStateProtein(FiveStateScore* query,ComplexSequence* target ,DPEnvelope * dpenv,DPRunImpl * dpri);
#define PackAln_bestmemory_FiveStateProtein Wise2_PackAln_bestmemory_FiveStateProtein


/* Function:  allocate_Expl_FiveStateProtein(query,target,dpri)
 *
 * Descrip:    This function allocates the FiveStateProtein structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_FiveStateProtein_only
 *
 *
 * Arg:         query [UNKN ] query data structure [FiveStateScore*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:          dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [FiveStateProtein *]
 *
 */
FiveStateProtein * Wise2_allocate_Expl_FiveStateProtein(FiveStateScore* query,ComplexSequence* target ,DPRunImpl * dpri);
#define allocate_Expl_FiveStateProtein Wise2_allocate_Expl_FiveStateProtein


/* Function:  recalculate_PackAln_FiveStateProtein(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by FiveStateProtein
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [FiveStateProtein *]
 *
 */
void Wise2_recalculate_PackAln_FiveStateProtein(PackAln * pal,FiveStateProtein * mat);
#define recalculate_PackAln_FiveStateProtein Wise2_recalculate_PackAln_FiveStateProtein


/* Function:  allocate_Small_FiveStateProtein(query,target)
 *
 * Descrip:    This function allocates the FiveStateProtein structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_FiveStateProtein_only
 *
 *
 * Arg:         query [UNKN ] query data structure [FiveStateScore*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 *
 * Return [UNKN ]  Undocumented return value [FiveStateProtein *]
 *
 */
FiveStateProtein * Wise2_allocate_Small_FiveStateProtein(FiveStateScore* query,ComplexSequence* target );
#define allocate_Small_FiveStateProtein Wise2_allocate_Small_FiveStateProtein


/* Function:  PackAln_calculate_Small_FiveStateProtein(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for FiveStateProtein structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_FiveStateProtein 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_FiveStateProtein 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [FiveStateProtein *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_calculate_Small_FiveStateProtein(FiveStateProtein * mat,DPEnvelope * dpenv);
#define PackAln_calculate_Small_FiveStateProtein Wise2_PackAln_calculate_Small_FiveStateProtein


/* Function:  AlnRangeSet_calculate_Small_FiveStateProtein(mat)
 *
 * Descrip:    This function calculates an alignment for FiveStateProtein structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_FiveStateProtein 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_FiveStateProtein
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_FiveStateProtein 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [FiveStateProtein *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_calculate_Small_FiveStateProtein(FiveStateProtein * mat);
#define AlnRangeSet_calculate_Small_FiveStateProtein Wise2_AlnRangeSet_calculate_Small_FiveStateProtein


/* Function:  AlnRangeSet_from_FiveStateProtein(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for FiveStateProtein structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_FiveStateProtein 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_FiveStateProtein
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [FiveStateProtein *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_from_FiveStateProtein(FiveStateProtein * mat);
#define AlnRangeSet_from_FiveStateProtein Wise2_AlnRangeSet_from_FiveStateProtein


/* Function:  convert_PackAln_to_AlnBlock_FiveStateProtein(pal)
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
AlnBlock * Wise2_convert_PackAln_to_AlnBlock_FiveStateProtein(PackAln * pal);
#define convert_PackAln_to_AlnBlock_FiveStateProtein Wise2_convert_PackAln_to_AlnBlock_FiveStateProtein


/* Function:  PackAln_read_Expl_FiveStateProtein(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [FiveStateProtein *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Expl_FiveStateProtein(FiveStateProtein * mat);
#define PackAln_read_Expl_FiveStateProtein Wise2_PackAln_read_Expl_FiveStateProtein


/* Function:  PackAln_read_generic_FiveStateProtein(mat,h)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [FiveStateProtein *]
 * Arg:          h [UNKN ] Undocumented argument [FiveStateProtein_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_generic_FiveStateProtein(FiveStateProtein * mat,FiveStateProtein_access_func_holder h);
#define PackAln_read_generic_FiveStateProtein Wise2_PackAln_read_generic_FiveStateProtein


/* Function:  calculate_FiveStateProtein(mat)
 *
 * Descrip:    This function calculates the FiveStateProtein matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_FiveStateProtein
 *
 *
 * Arg:        mat [UNKN ] FiveStateProtein which contains explicit basematrix memory [FiveStateProtein *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_FiveStateProtein(FiveStateProtein * mat);
#define calculate_FiveStateProtein Wise2_calculate_FiveStateProtein


/* Function:  calculate_dpenv_FiveStateProtein(mat,dpenv)
 *
 * Descrip:    This function calculates the FiveStateProtein matrix when in explicit mode, subject to the envelope
 *
 *
 * Arg:          mat [UNKN ] FiveStateProtein which contains explicit basematrix memory [FiveStateProtein *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_dpenv_FiveStateProtein(FiveStateProtein * mat,DPEnvelope * dpenv);
#define calculate_dpenv_FiveStateProtein Wise2_calculate_dpenv_FiveStateProtein


/* Function:  FiveStateProtein_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [FiveStateProtein *]
 *
 */
FiveStateProtein * Wise2_FiveStateProtein_alloc(void);
#define FiveStateProtein_alloc Wise2_FiveStateProtein_alloc


/* Function:  free_FiveStateProtein(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [FiveStateProtein *]
 *
 * Return [UNKN ]  Undocumented return value [FiveStateProtein *]
 *
 */
FiveStateProtein * Wise2_free_FiveStateProtein(FiveStateProtein * obj);
#define free_FiveStateProtein Wise2_free_FiveStateProtein


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
int Wise2_FiveStateProtein_shatter_access_main(FiveStateProtein * mat,int i,int j,int state);
#define FiveStateProtein_shatter_access_main Wise2_FiveStateProtein_shatter_access_main
int Wise2_FiveStateProtein_shatter_access_special(FiveStateProtein * mat,int i,int j,int state);
#define FiveStateProtein_shatter_access_special Wise2_FiveStateProtein_shatter_access_special
void * Wise2_thread_loop_FiveStateProtein(void * ptr);
#define thread_loop_FiveStateProtein Wise2_thread_loop_FiveStateProtein
int Wise2_score_only_FiveStateProtein(FiveStateScore* query,ComplexSequence* target );
#define score_only_FiveStateProtein Wise2_score_only_FiveStateProtein
FiveStateProtein * Wise2_allocate_FiveStateProtein_only(FiveStateScore* query,ComplexSequence* target );
#define allocate_FiveStateProtein_only Wise2_allocate_FiveStateProtein_only
void Wise2_init_FiveStateProtein(FiveStateProtein * mat);
#define init_FiveStateProtein Wise2_init_FiveStateProtein
AlnRange * Wise2_AlnRange_build_FiveStateProtein(FiveStateProtein * mat,int stopj,int stopspecstate,int * startj,int * startspecstate);
#define AlnRange_build_FiveStateProtein Wise2_AlnRange_build_FiveStateProtein
boolean Wise2_read_hidden_FiveStateProtein(FiveStateProtein * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out);
#define read_hidden_FiveStateProtein Wise2_read_hidden_FiveStateProtein
int Wise2_max_hidden_FiveStateProtein(FiveStateProtein * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_hidden_FiveStateProtein Wise2_max_hidden_FiveStateProtein
boolean Wise2_read_special_strip_FiveStateProtein(FiveStateProtein * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out);
#define read_special_strip_FiveStateProtein Wise2_read_special_strip_FiveStateProtein
int Wise2_max_special_strip_FiveStateProtein(FiveStateProtein * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_special_strip_FiveStateProtein Wise2_max_special_strip_FiveStateProtein
int Wise2_max_matrix_to_special_FiveStateProtein(FiveStateProtein * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_matrix_to_special_FiveStateProtein Wise2_max_matrix_to_special_FiveStateProtein
void Wise2_calculate_hidden_FiveStateProtein(FiveStateProtein * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv);
#define calculate_hidden_FiveStateProtein Wise2_calculate_hidden_FiveStateProtein
void Wise2_init_hidden_FiveStateProtein(FiveStateProtein * mat,int starti,int startj,int stopi,int stopj);
#define init_hidden_FiveStateProtein Wise2_init_hidden_FiveStateProtein
boolean Wise2_full_dc_FiveStateProtein(FiveStateProtein * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv);
#define full_dc_FiveStateProtein Wise2_full_dc_FiveStateProtein
boolean Wise2_do_dc_single_pass_FiveStateProtein(FiveStateProtein * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done);
#define do_dc_single_pass_FiveStateProtein Wise2_do_dc_single_pass_FiveStateProtein
void Wise2_push_dc_at_merge_FiveStateProtein(FiveStateProtein * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv);
#define push_dc_at_merge_FiveStateProtein Wise2_push_dc_at_merge_FiveStateProtein
void Wise2_follow_on_dc_FiveStateProtein(FiveStateProtein * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define follow_on_dc_FiveStateProtein Wise2_follow_on_dc_FiveStateProtein
void Wise2_run_up_dc_FiveStateProtein(FiveStateProtein * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define run_up_dc_FiveStateProtein Wise2_run_up_dc_FiveStateProtein
void Wise2_init_dc_FiveStateProtein(FiveStateProtein * mat);
#define init_dc_FiveStateProtein Wise2_init_dc_FiveStateProtein
int Wise2_start_end_find_end_FiveStateProtein(FiveStateProtein * mat,int * endj);
#define start_end_find_end_FiveStateProtein Wise2_start_end_find_end_FiveStateProtein
boolean Wise2_dc_optimised_start_end_calc_FiveStateProtein(FiveStateProtein *mat,DPEnvelope * dpenv);
#define dc_optimised_start_end_calc_FiveStateProtein Wise2_dc_optimised_start_end_calc_FiveStateProtein
void Wise2_init_start_end_linear_FiveStateProtein(FiveStateProtein * mat);
#define init_start_end_linear_FiveStateProtein Wise2_init_start_end_linear_FiveStateProtein
AlnConvertSet * Wise2_AlnConvertSet_FiveStateProtein(void);
#define AlnConvertSet_FiveStateProtein Wise2_AlnConvertSet_FiveStateProtein
int Wise2_FiveStateProtein_explicit_access_main(FiveStateProtein * mat,int i,int j,int state);
#define FiveStateProtein_explicit_access_main Wise2_FiveStateProtein_explicit_access_main
int Wise2_FiveStateProtein_explicit_access_special(FiveStateProtein * mat,int i,int j,int state);
#define FiveStateProtein_explicit_access_special Wise2_FiveStateProtein_explicit_access_special
int Wise2_find_end_FiveStateProtein(FiveStateProtein * mat,int * ri,int * rj,int * state,boolean * isspecial,FiveStateProtein_access_func_holder h);
#define find_end_FiveStateProtein Wise2_find_end_FiveStateProtein
void Wise2_FiveStateProtein_debug_show_matrix(FiveStateProtein * mat,int starti,int stopi,int startj,int stopj,FILE * ofp);
#define FiveStateProtein_debug_show_matrix Wise2_FiveStateProtein_debug_show_matrix
int Wise2_max_calc_FiveStateProtein(FiveStateProtein * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,FiveStateProtein_access_func_holder h);
#define max_calc_FiveStateProtein Wise2_max_calc_FiveStateProtein
int Wise2_max_calc_special_FiveStateProtein(FiveStateProtein * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,FiveStateProtein_access_func_holder h);
#define max_calc_special_FiveStateProtein Wise2_max_calc_special_FiveStateProtein

#ifdef _cplusplus
}
#endif

#endif
