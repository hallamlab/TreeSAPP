#ifndef DYNAMITElocaldbaHEADERFILE
#define DYNAMITElocaldbaHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"
#include "localcispara.h"

struct Wise2_LocalDnaMatchBlock {  
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
    LocalCisHitScore* lchs;  
    } ;  
/* LocalDnaMatchBlock defined */ 
#ifndef DYNAMITE_DEFINED_LocalDnaMatchBlock
typedef struct Wise2_LocalDnaMatchBlock Wise2_LocalDnaMatchBlock;
#define LocalDnaMatchBlock Wise2_LocalDnaMatchBlock
#define DYNAMITE_DEFINED_LocalDnaMatchBlock
#endif


struct Wise2_LocalDnaMatchBlock_access_func_holder {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int (*access_main)(LocalDnaMatchBlock*,int,int,int); 
    int (*access_special)(LocalDnaMatchBlock*,int,int,int);  
    } ;  
/* LocalDnaMatchBlock_access_func_holder defined */ 
#ifndef DYNAMITE_DEFINED_LocalDnaMatchBlock_access_func_holder
typedef struct Wise2_LocalDnaMatchBlock_access_func_holder Wise2_LocalDnaMatchBlock_access_func_holder;
#define LocalDnaMatchBlock_access_func_holder Wise2_LocalDnaMatchBlock_access_func_holder
#define DYNAMITE_DEFINED_LocalDnaMatchBlock_access_func_holder
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  PackAln_read_Shatter_LocalDnaMatchBlock(mat)
 *
 * Descrip:    Reads off PackAln from shatter matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [LocalDnaMatchBlock *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Shatter_LocalDnaMatchBlock(LocalDnaMatchBlock * mat);
#define PackAln_read_Shatter_LocalDnaMatchBlock Wise2_PackAln_read_Shatter_LocalDnaMatchBlock


/* Function:  calculate_shatter_LocalDnaMatchBlock(mat,dpenv)
 *
 * Descrip:    This function calculates the LocalDnaMatchBlock matrix when in shatter mode
 *
 *
 * Arg:          mat [UNKN ] (null) [LocalDnaMatchBlock *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_shatter_LocalDnaMatchBlock(LocalDnaMatchBlock * mat,DPEnvelope * dpenv);
#define calculate_shatter_LocalDnaMatchBlock Wise2_calculate_shatter_LocalDnaMatchBlock


/* Function:  search_LocalDnaMatchBlock(dbsi,out,query,target,lchs)
 *
 * Descrip:    This function makes a database search of LocalDnaMatchBlock
 *             It uses the dbsi structure to choose which implementation
 *             to use of the database searching. This way at run time you
 *             can switch between single threaded/multi-threaded or hardware
 *
 *
 * Arg:          dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 * Arg:           out [UNKN ] Undocumented argument [Hscore *]
 * Arg:         query [UNKN ] Undocumented argument [ComplexSequence*]
 * Arg:        target [UNKN ] Undocumented argument [ComplexSequence*]
 * Arg:          lchs [UNKN ] Undocumented argument [LocalCisHitScore*]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_search_LocalDnaMatchBlock(DBSearchImpl * dbsi,Hscore * out,ComplexSequence* query,ComplexSequence* target ,LocalCisHitScore* lchs);
#define search_LocalDnaMatchBlock Wise2_search_LocalDnaMatchBlock


/* Function:  serial_search_LocalDnaMatchBlock(out,query,target,lchs)
 *
 * Descrip:    This function makes a database search of LocalDnaMatchBlock
 *             It is a single processor implementation
 *
 *
 * Arg:           out [UNKN ] Undocumented argument [Hscore *]
 * Arg:         query [UNKN ] Undocumented argument [ComplexSequence*]
 * Arg:        target [UNKN ] Undocumented argument [ComplexSequence*]
 * Arg:          lchs [UNKN ] Undocumented argument [LocalCisHitScore*]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_serial_search_LocalDnaMatchBlock(Hscore * out,ComplexSequence* query,ComplexSequence* target ,LocalCisHitScore* lchs);
#define serial_search_LocalDnaMatchBlock Wise2_serial_search_LocalDnaMatchBlock


/* Function:  PackAln_bestmemory_LocalDnaMatchBlock(query,target,lchs,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_LocalDnaMatchBlock
 *
 *
 * Arg:         query [UNKN ] query data structure [ComplexSequence*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:          lchs [UNKN ] Resource [LocalCisHitScore*]
 * Arg:         dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:          dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_bestmemory_LocalDnaMatchBlock(ComplexSequence* query,ComplexSequence* target ,LocalCisHitScore* lchs,DPEnvelope * dpenv,DPRunImpl * dpri);
#define PackAln_bestmemory_LocalDnaMatchBlock Wise2_PackAln_bestmemory_LocalDnaMatchBlock


/* Function:  allocate_Expl_LocalDnaMatchBlock(query,target,lchs,dpri)
 *
 * Descrip:    This function allocates the LocalDnaMatchBlock structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_LocalDnaMatchBlock_only
 *
 *
 * Arg:         query [UNKN ] query data structure [ComplexSequence*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:          lchs [UNKN ] Resource [LocalCisHitScore*]
 * Arg:          dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [LocalDnaMatchBlock *]
 *
 */
LocalDnaMatchBlock * Wise2_allocate_Expl_LocalDnaMatchBlock(ComplexSequence* query,ComplexSequence* target ,LocalCisHitScore* lchs,DPRunImpl * dpri);
#define allocate_Expl_LocalDnaMatchBlock Wise2_allocate_Expl_LocalDnaMatchBlock


/* Function:  recalculate_PackAln_LocalDnaMatchBlock(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by LocalDnaMatchBlock
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [LocalDnaMatchBlock *]
 *
 */
void Wise2_recalculate_PackAln_LocalDnaMatchBlock(PackAln * pal,LocalDnaMatchBlock * mat);
#define recalculate_PackAln_LocalDnaMatchBlock Wise2_recalculate_PackAln_LocalDnaMatchBlock


/* Function:  allocate_Small_LocalDnaMatchBlock(query,target,lchs)
 *
 * Descrip:    This function allocates the LocalDnaMatchBlock structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_LocalDnaMatchBlock_only
 *
 *
 * Arg:         query [UNKN ] query data structure [ComplexSequence*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:          lchs [UNKN ] Resource [LocalCisHitScore*]
 *
 * Return [UNKN ]  Undocumented return value [LocalDnaMatchBlock *]
 *
 */
LocalDnaMatchBlock * Wise2_allocate_Small_LocalDnaMatchBlock(ComplexSequence* query,ComplexSequence* target ,LocalCisHitScore* lchs);
#define allocate_Small_LocalDnaMatchBlock Wise2_allocate_Small_LocalDnaMatchBlock


/* Function:  PackAln_calculate_Small_LocalDnaMatchBlock(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for LocalDnaMatchBlock structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_LocalDnaMatchBlock 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_LocalDnaMatchBlock 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [LocalDnaMatchBlock *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_calculate_Small_LocalDnaMatchBlock(LocalDnaMatchBlock * mat,DPEnvelope * dpenv);
#define PackAln_calculate_Small_LocalDnaMatchBlock Wise2_PackAln_calculate_Small_LocalDnaMatchBlock


/* Function:  AlnRangeSet_calculate_Small_LocalDnaMatchBlock(mat)
 *
 * Descrip:    This function calculates an alignment for LocalDnaMatchBlock structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_LocalDnaMatchBlock 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_LocalDnaMatchBlock
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_LocalDnaMatchBlock 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [LocalDnaMatchBlock *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_calculate_Small_LocalDnaMatchBlock(LocalDnaMatchBlock * mat);
#define AlnRangeSet_calculate_Small_LocalDnaMatchBlock Wise2_AlnRangeSet_calculate_Small_LocalDnaMatchBlock


/* Function:  AlnRangeSet_from_LocalDnaMatchBlock(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for LocalDnaMatchBlock structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_LocalDnaMatchBlock 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_LocalDnaMatchBlock
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [LocalDnaMatchBlock *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_from_LocalDnaMatchBlock(LocalDnaMatchBlock * mat);
#define AlnRangeSet_from_LocalDnaMatchBlock Wise2_AlnRangeSet_from_LocalDnaMatchBlock


/* Function:  convert_PackAln_to_AlnBlock_LocalDnaMatchBlock(pal)
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
AlnBlock * Wise2_convert_PackAln_to_AlnBlock_LocalDnaMatchBlock(PackAln * pal);
#define convert_PackAln_to_AlnBlock_LocalDnaMatchBlock Wise2_convert_PackAln_to_AlnBlock_LocalDnaMatchBlock


/* Function:  PackAln_read_Expl_LocalDnaMatchBlock(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [LocalDnaMatchBlock *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Expl_LocalDnaMatchBlock(LocalDnaMatchBlock * mat);
#define PackAln_read_Expl_LocalDnaMatchBlock Wise2_PackAln_read_Expl_LocalDnaMatchBlock


/* Function:  PackAln_read_generic_LocalDnaMatchBlock(mat,h)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [LocalDnaMatchBlock *]
 * Arg:          h [UNKN ] Undocumented argument [LocalDnaMatchBlock_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_generic_LocalDnaMatchBlock(LocalDnaMatchBlock * mat,LocalDnaMatchBlock_access_func_holder h);
#define PackAln_read_generic_LocalDnaMatchBlock Wise2_PackAln_read_generic_LocalDnaMatchBlock


/* Function:  calculate_LocalDnaMatchBlock(mat)
 *
 * Descrip:    This function calculates the LocalDnaMatchBlock matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_LocalDnaMatchBlock
 *
 *
 * Arg:        mat [UNKN ] LocalDnaMatchBlock which contains explicit basematrix memory [LocalDnaMatchBlock *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_LocalDnaMatchBlock(LocalDnaMatchBlock * mat);
#define calculate_LocalDnaMatchBlock Wise2_calculate_LocalDnaMatchBlock


/* Function:  calculate_dpenv_LocalDnaMatchBlock(mat,dpenv)
 *
 * Descrip:    This function calculates the LocalDnaMatchBlock matrix when in explicit mode, subject to the envelope
 *
 *
 * Arg:          mat [UNKN ] LocalDnaMatchBlock which contains explicit basematrix memory [LocalDnaMatchBlock *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_dpenv_LocalDnaMatchBlock(LocalDnaMatchBlock * mat,DPEnvelope * dpenv);
#define calculate_dpenv_LocalDnaMatchBlock Wise2_calculate_dpenv_LocalDnaMatchBlock


/* Function:  LocalDnaMatchBlock_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LocalDnaMatchBlock *]
 *
 */
LocalDnaMatchBlock * Wise2_LocalDnaMatchBlock_alloc(void);
#define LocalDnaMatchBlock_alloc Wise2_LocalDnaMatchBlock_alloc


/* Function:  free_LocalDnaMatchBlock(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [LocalDnaMatchBlock *]
 *
 * Return [UNKN ]  Undocumented return value [LocalDnaMatchBlock *]
 *
 */
LocalDnaMatchBlock * Wise2_free_LocalDnaMatchBlock(LocalDnaMatchBlock * obj);
#define free_LocalDnaMatchBlock Wise2_free_LocalDnaMatchBlock


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
int Wise2_LocalDnaMatchBlock_shatter_access_main(LocalDnaMatchBlock * mat,int i,int j,int state);
#define LocalDnaMatchBlock_shatter_access_main Wise2_LocalDnaMatchBlock_shatter_access_main
int Wise2_LocalDnaMatchBlock_shatter_access_special(LocalDnaMatchBlock * mat,int i,int j,int state);
#define LocalDnaMatchBlock_shatter_access_special Wise2_LocalDnaMatchBlock_shatter_access_special
int Wise2_score_only_LocalDnaMatchBlock(ComplexSequence* query,ComplexSequence* target ,LocalCisHitScore* lchs);
#define score_only_LocalDnaMatchBlock Wise2_score_only_LocalDnaMatchBlock
LocalDnaMatchBlock * Wise2_allocate_LocalDnaMatchBlock_only(ComplexSequence* query,ComplexSequence* target ,LocalCisHitScore* lchs);
#define allocate_LocalDnaMatchBlock_only Wise2_allocate_LocalDnaMatchBlock_only
void Wise2_init_LocalDnaMatchBlock(LocalDnaMatchBlock * mat);
#define init_LocalDnaMatchBlock Wise2_init_LocalDnaMatchBlock
AlnRange * Wise2_AlnRange_build_LocalDnaMatchBlock(LocalDnaMatchBlock * mat,int stopj,int stopspecstate,int * startj,int * startspecstate);
#define AlnRange_build_LocalDnaMatchBlock Wise2_AlnRange_build_LocalDnaMatchBlock
boolean Wise2_read_hidden_LocalDnaMatchBlock(LocalDnaMatchBlock * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out);
#define read_hidden_LocalDnaMatchBlock Wise2_read_hidden_LocalDnaMatchBlock
int Wise2_max_hidden_LocalDnaMatchBlock(LocalDnaMatchBlock * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_hidden_LocalDnaMatchBlock Wise2_max_hidden_LocalDnaMatchBlock
boolean Wise2_read_special_strip_LocalDnaMatchBlock(LocalDnaMatchBlock * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out);
#define read_special_strip_LocalDnaMatchBlock Wise2_read_special_strip_LocalDnaMatchBlock
int Wise2_max_special_strip_LocalDnaMatchBlock(LocalDnaMatchBlock * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_special_strip_LocalDnaMatchBlock Wise2_max_special_strip_LocalDnaMatchBlock
int Wise2_max_matrix_to_special_LocalDnaMatchBlock(LocalDnaMatchBlock * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_matrix_to_special_LocalDnaMatchBlock Wise2_max_matrix_to_special_LocalDnaMatchBlock
void Wise2_calculate_hidden_LocalDnaMatchBlock(LocalDnaMatchBlock * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv);
#define calculate_hidden_LocalDnaMatchBlock Wise2_calculate_hidden_LocalDnaMatchBlock
void Wise2_init_hidden_LocalDnaMatchBlock(LocalDnaMatchBlock * mat,int starti,int startj,int stopi,int stopj);
#define init_hidden_LocalDnaMatchBlock Wise2_init_hidden_LocalDnaMatchBlock
boolean Wise2_full_dc_LocalDnaMatchBlock(LocalDnaMatchBlock * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv);
#define full_dc_LocalDnaMatchBlock Wise2_full_dc_LocalDnaMatchBlock
boolean Wise2_do_dc_single_pass_LocalDnaMatchBlock(LocalDnaMatchBlock * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done);
#define do_dc_single_pass_LocalDnaMatchBlock Wise2_do_dc_single_pass_LocalDnaMatchBlock
void Wise2_push_dc_at_merge_LocalDnaMatchBlock(LocalDnaMatchBlock * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv);
#define push_dc_at_merge_LocalDnaMatchBlock Wise2_push_dc_at_merge_LocalDnaMatchBlock
void Wise2_follow_on_dc_LocalDnaMatchBlock(LocalDnaMatchBlock * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define follow_on_dc_LocalDnaMatchBlock Wise2_follow_on_dc_LocalDnaMatchBlock
void Wise2_run_up_dc_LocalDnaMatchBlock(LocalDnaMatchBlock * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define run_up_dc_LocalDnaMatchBlock Wise2_run_up_dc_LocalDnaMatchBlock
void Wise2_init_dc_LocalDnaMatchBlock(LocalDnaMatchBlock * mat);
#define init_dc_LocalDnaMatchBlock Wise2_init_dc_LocalDnaMatchBlock
int Wise2_start_end_find_end_LocalDnaMatchBlock(LocalDnaMatchBlock * mat,int * endj);
#define start_end_find_end_LocalDnaMatchBlock Wise2_start_end_find_end_LocalDnaMatchBlock
boolean Wise2_dc_optimised_start_end_calc_LocalDnaMatchBlock(LocalDnaMatchBlock *mat,DPEnvelope * dpenv);
#define dc_optimised_start_end_calc_LocalDnaMatchBlock Wise2_dc_optimised_start_end_calc_LocalDnaMatchBlock
void Wise2_init_start_end_linear_LocalDnaMatchBlock(LocalDnaMatchBlock * mat);
#define init_start_end_linear_LocalDnaMatchBlock Wise2_init_start_end_linear_LocalDnaMatchBlock
AlnConvertSet * Wise2_AlnConvertSet_LocalDnaMatchBlock(void);
#define AlnConvertSet_LocalDnaMatchBlock Wise2_AlnConvertSet_LocalDnaMatchBlock
int Wise2_LocalDnaMatchBlock_explicit_access_main(LocalDnaMatchBlock * mat,int i,int j,int state);
#define LocalDnaMatchBlock_explicit_access_main Wise2_LocalDnaMatchBlock_explicit_access_main
int Wise2_LocalDnaMatchBlock_explicit_access_special(LocalDnaMatchBlock * mat,int i,int j,int state);
#define LocalDnaMatchBlock_explicit_access_special Wise2_LocalDnaMatchBlock_explicit_access_special
int Wise2_find_end_LocalDnaMatchBlock(LocalDnaMatchBlock * mat,int * ri,int * rj,int * state,boolean * isspecial,LocalDnaMatchBlock_access_func_holder h);
#define find_end_LocalDnaMatchBlock Wise2_find_end_LocalDnaMatchBlock
void Wise2_LocalDnaMatchBlock_debug_show_matrix(LocalDnaMatchBlock * mat,int starti,int stopi,int startj,int stopj,FILE * ofp);
#define LocalDnaMatchBlock_debug_show_matrix Wise2_LocalDnaMatchBlock_debug_show_matrix
int Wise2_max_calc_LocalDnaMatchBlock(LocalDnaMatchBlock * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,LocalDnaMatchBlock_access_func_holder h);
#define max_calc_LocalDnaMatchBlock Wise2_max_calc_LocalDnaMatchBlock
int Wise2_max_calc_special_LocalDnaMatchBlock(LocalDnaMatchBlock * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,LocalDnaMatchBlock_access_func_holder h);
#define max_calc_special_LocalDnaMatchBlock Wise2_max_calc_special_LocalDnaMatchBlock

#ifdef _cplusplus
}
#endif

#endif
