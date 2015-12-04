#ifndef DYNAMITEdbaHEADERFILE
#define DYNAMITEdbaHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"

struct Wise2_DnaMatchBlock {  
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
    DnaMatrix* comp65;   
    DnaMatrix* comp75;   
    DnaMatrix* comp85;   
    DnaMatrix* comp95;   
    Score g;     
    Score u;     
    Score v;     
    Score s;     
    Score b;     
    } ;  
/* DnaMatchBlock defined */ 
#ifndef DYNAMITE_DEFINED_DnaMatchBlock
typedef struct Wise2_DnaMatchBlock Wise2_DnaMatchBlock;
#define DnaMatchBlock Wise2_DnaMatchBlock
#define DYNAMITE_DEFINED_DnaMatchBlock
#endif


struct Wise2_DnaMatchBlock_access_func_holder {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int (*access_main)(DnaMatchBlock*,int,int,int);  
    int (*access_special)(DnaMatchBlock*,int,int,int);   
    } ;  
/* DnaMatchBlock_access_func_holder defined */ 
#ifndef DYNAMITE_DEFINED_DnaMatchBlock_access_func_holder
typedef struct Wise2_DnaMatchBlock_access_func_holder Wise2_DnaMatchBlock_access_func_holder;
#define DnaMatchBlock_access_func_holder Wise2_DnaMatchBlock_access_func_holder
#define DYNAMITE_DEFINED_DnaMatchBlock_access_func_holder
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  PackAln_read_Shatter_DnaMatchBlock(mat)
 *
 * Descrip:    Reads off PackAln from shatter matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [DnaMatchBlock *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Shatter_DnaMatchBlock(DnaMatchBlock * mat);
#define PackAln_read_Shatter_DnaMatchBlock Wise2_PackAln_read_Shatter_DnaMatchBlock


/* Function:  calculate_shatter_DnaMatchBlock(mat,dpenv)
 *
 * Descrip:    This function calculates the DnaMatchBlock matrix when in shatter mode
 *
 *
 * Arg:          mat [UNKN ] (null) [DnaMatchBlock *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_shatter_DnaMatchBlock(DnaMatchBlock * mat,DPEnvelope * dpenv);
#define calculate_shatter_DnaMatchBlock Wise2_calculate_shatter_DnaMatchBlock


/* Function:  search_DnaMatchBlock(dbsi,out,query,target,comp65,comp75,comp85,comp95,g,u,v,s,b)
 *
 * Descrip:    This function makes a database search of DnaMatchBlock
 *             It uses the dbsi structure to choose which implementation
 *             to use of the database searching. This way at run time you
 *             can switch between single threaded/multi-threaded or hardware
 *
 *
 * Arg:          dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 * Arg:           out [UNKN ] Undocumented argument [Hscore *]
 * Arg:         query [UNKN ] Undocumented argument [ComplexSequence*]
 * Arg:        target [UNKN ] Undocumented argument [ComplexSequence*]
 * Arg:        comp65 [UNKN ] Undocumented argument [DnaMatrix*]
 * Arg:        comp75 [UNKN ] Undocumented argument [DnaMatrix*]
 * Arg:        comp85 [UNKN ] Undocumented argument [DnaMatrix*]
 * Arg:        comp95 [UNKN ] Undocumented argument [DnaMatrix*]
 * Arg:             g [UNKN ] Undocumented argument [Score]
 * Arg:             u [UNKN ] Undocumented argument [Score]
 * Arg:             v [UNKN ] Undocumented argument [Score]
 * Arg:             s [UNKN ] Undocumented argument [Score]
 * Arg:             b [UNKN ] Undocumented argument [Score]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_search_DnaMatchBlock(DBSearchImpl * dbsi,Hscore * out,ComplexSequence* query,ComplexSequence* target ,DnaMatrix* comp65,DnaMatrix* comp75,DnaMatrix* comp85,DnaMatrix* comp95,Score g,Score u,Score v,Score s,Score b);
#define search_DnaMatchBlock Wise2_search_DnaMatchBlock


/* Function:  serial_search_DnaMatchBlock(out,query,target,comp65,comp75,comp85,comp95,g,u,v,s,b)
 *
 * Descrip:    This function makes a database search of DnaMatchBlock
 *             It is a single processor implementation
 *
 *
 * Arg:           out [UNKN ] Undocumented argument [Hscore *]
 * Arg:         query [UNKN ] Undocumented argument [ComplexSequence*]
 * Arg:        target [UNKN ] Undocumented argument [ComplexSequence*]
 * Arg:        comp65 [UNKN ] Undocumented argument [DnaMatrix*]
 * Arg:        comp75 [UNKN ] Undocumented argument [DnaMatrix*]
 * Arg:        comp85 [UNKN ] Undocumented argument [DnaMatrix*]
 * Arg:        comp95 [UNKN ] Undocumented argument [DnaMatrix*]
 * Arg:             g [UNKN ] Undocumented argument [Score]
 * Arg:             u [UNKN ] Undocumented argument [Score]
 * Arg:             v [UNKN ] Undocumented argument [Score]
 * Arg:             s [UNKN ] Undocumented argument [Score]
 * Arg:             b [UNKN ] Undocumented argument [Score]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_serial_search_DnaMatchBlock(Hscore * out,ComplexSequence* query,ComplexSequence* target ,DnaMatrix* comp65,DnaMatrix* comp75,DnaMatrix* comp85,DnaMatrix* comp95,Score g,Score u,Score v,Score s,Score b);
#define serial_search_DnaMatchBlock Wise2_serial_search_DnaMatchBlock


/* Function:  PackAln_bestmemory_DnaMatchBlock(query,target,comp65,comp75,comp85,comp95,g,u,v,s,b,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_DnaMatchBlock
 *
 *
 * Arg:         query [UNKN ] query data structure [ComplexSequence*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:        comp65 [UNKN ] Resource [DnaMatrix*]
 * Arg:        comp75 [UNKN ] Resource [DnaMatrix*]
 * Arg:        comp85 [UNKN ] Resource [DnaMatrix*]
 * Arg:        comp95 [UNKN ] Resource [DnaMatrix*]
 * Arg:             g [UNKN ] Resource [Score]
 * Arg:             u [UNKN ] Resource [Score]
 * Arg:             v [UNKN ] Resource [Score]
 * Arg:             s [UNKN ] Resource [Score]
 * Arg:             b [UNKN ] Resource [Score]
 * Arg:         dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:          dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_bestmemory_DnaMatchBlock(ComplexSequence* query,ComplexSequence* target ,DnaMatrix* comp65,DnaMatrix* comp75,DnaMatrix* comp85,DnaMatrix* comp95,Score g,Score u,Score v,Score s,Score b,DPEnvelope * dpenv,DPRunImpl * dpri);
#define PackAln_bestmemory_DnaMatchBlock Wise2_PackAln_bestmemory_DnaMatchBlock


/* Function:  allocate_Expl_DnaMatchBlock(query,target,comp65,comp75,comp85,comp95,g,u,v,s,b,dpri)
 *
 * Descrip:    This function allocates the DnaMatchBlock structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_DnaMatchBlock_only
 *
 *
 * Arg:         query [UNKN ] query data structure [ComplexSequence*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:        comp65 [UNKN ] Resource [DnaMatrix*]
 * Arg:        comp75 [UNKN ] Resource [DnaMatrix*]
 * Arg:        comp85 [UNKN ] Resource [DnaMatrix*]
 * Arg:        comp95 [UNKN ] Resource [DnaMatrix*]
 * Arg:             g [UNKN ] Resource [Score]
 * Arg:             u [UNKN ] Resource [Score]
 * Arg:             v [UNKN ] Resource [Score]
 * Arg:             s [UNKN ] Resource [Score]
 * Arg:             b [UNKN ] Resource [Score]
 * Arg:          dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [DnaMatchBlock *]
 *
 */
DnaMatchBlock * Wise2_allocate_Expl_DnaMatchBlock(ComplexSequence* query,ComplexSequence* target ,DnaMatrix* comp65,DnaMatrix* comp75,DnaMatrix* comp85,DnaMatrix* comp95,Score g,Score u,Score v,Score s,Score b,DPRunImpl * dpri);
#define allocate_Expl_DnaMatchBlock Wise2_allocate_Expl_DnaMatchBlock


/* Function:  recalculate_PackAln_DnaMatchBlock(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by DnaMatchBlock
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [DnaMatchBlock *]
 *
 */
void Wise2_recalculate_PackAln_DnaMatchBlock(PackAln * pal,DnaMatchBlock * mat);
#define recalculate_PackAln_DnaMatchBlock Wise2_recalculate_PackAln_DnaMatchBlock


/* Function:  allocate_Small_DnaMatchBlock(query,target,comp65,comp75,comp85,comp95,g,u,v,s,b)
 *
 * Descrip:    This function allocates the DnaMatchBlock structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_DnaMatchBlock_only
 *
 *
 * Arg:         query [UNKN ] query data structure [ComplexSequence*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:        comp65 [UNKN ] Resource [DnaMatrix*]
 * Arg:        comp75 [UNKN ] Resource [DnaMatrix*]
 * Arg:        comp85 [UNKN ] Resource [DnaMatrix*]
 * Arg:        comp95 [UNKN ] Resource [DnaMatrix*]
 * Arg:             g [UNKN ] Resource [Score]
 * Arg:             u [UNKN ] Resource [Score]
 * Arg:             v [UNKN ] Resource [Score]
 * Arg:             s [UNKN ] Resource [Score]
 * Arg:             b [UNKN ] Resource [Score]
 *
 * Return [UNKN ]  Undocumented return value [DnaMatchBlock *]
 *
 */
DnaMatchBlock * Wise2_allocate_Small_DnaMatchBlock(ComplexSequence* query,ComplexSequence* target ,DnaMatrix* comp65,DnaMatrix* comp75,DnaMatrix* comp85,DnaMatrix* comp95,Score g,Score u,Score v,Score s,Score b);
#define allocate_Small_DnaMatchBlock Wise2_allocate_Small_DnaMatchBlock


/* Function:  PackAln_calculate_Small_DnaMatchBlock(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for DnaMatchBlock structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_DnaMatchBlock 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_DnaMatchBlock 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [DnaMatchBlock *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_calculate_Small_DnaMatchBlock(DnaMatchBlock * mat,DPEnvelope * dpenv);
#define PackAln_calculate_Small_DnaMatchBlock Wise2_PackAln_calculate_Small_DnaMatchBlock


/* Function:  AlnRangeSet_calculate_Small_DnaMatchBlock(mat)
 *
 * Descrip:    This function calculates an alignment for DnaMatchBlock structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_DnaMatchBlock 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_DnaMatchBlock
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_DnaMatchBlock 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [DnaMatchBlock *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_calculate_Small_DnaMatchBlock(DnaMatchBlock * mat);
#define AlnRangeSet_calculate_Small_DnaMatchBlock Wise2_AlnRangeSet_calculate_Small_DnaMatchBlock


/* Function:  AlnRangeSet_from_DnaMatchBlock(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for DnaMatchBlock structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_DnaMatchBlock 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_DnaMatchBlock
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [DnaMatchBlock *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_from_DnaMatchBlock(DnaMatchBlock * mat);
#define AlnRangeSet_from_DnaMatchBlock Wise2_AlnRangeSet_from_DnaMatchBlock


/* Function:  convert_PackAln_to_AlnBlock_DnaMatchBlock(pal)
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
AlnBlock * Wise2_convert_PackAln_to_AlnBlock_DnaMatchBlock(PackAln * pal);
#define convert_PackAln_to_AlnBlock_DnaMatchBlock Wise2_convert_PackAln_to_AlnBlock_DnaMatchBlock


/* Function:  PackAln_read_Expl_DnaMatchBlock(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [DnaMatchBlock *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Expl_DnaMatchBlock(DnaMatchBlock * mat);
#define PackAln_read_Expl_DnaMatchBlock Wise2_PackAln_read_Expl_DnaMatchBlock


/* Function:  PackAln_read_generic_DnaMatchBlock(mat,h)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [DnaMatchBlock *]
 * Arg:          h [UNKN ] Undocumented argument [DnaMatchBlock_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_generic_DnaMatchBlock(DnaMatchBlock * mat,DnaMatchBlock_access_func_holder h);
#define PackAln_read_generic_DnaMatchBlock Wise2_PackAln_read_generic_DnaMatchBlock


/* Function:  calculate_DnaMatchBlock(mat)
 *
 * Descrip:    This function calculates the DnaMatchBlock matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_DnaMatchBlock
 *
 *
 * Arg:        mat [UNKN ] DnaMatchBlock which contains explicit basematrix memory [DnaMatchBlock *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_DnaMatchBlock(DnaMatchBlock * mat);
#define calculate_DnaMatchBlock Wise2_calculate_DnaMatchBlock


/* Function:  calculate_dpenv_DnaMatchBlock(mat,dpenv)
 *
 * Descrip:    This function calculates the DnaMatchBlock matrix when in explicit mode, subject to the envelope
 *
 *
 * Arg:          mat [UNKN ] DnaMatchBlock which contains explicit basematrix memory [DnaMatchBlock *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_dpenv_DnaMatchBlock(DnaMatchBlock * mat,DPEnvelope * dpenv);
#define calculate_dpenv_DnaMatchBlock Wise2_calculate_dpenv_DnaMatchBlock


/* Function:  DnaMatchBlock_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaMatchBlock *]
 *
 */
DnaMatchBlock * Wise2_DnaMatchBlock_alloc(void);
#define DnaMatchBlock_alloc Wise2_DnaMatchBlock_alloc


/* Function:  free_DnaMatchBlock(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DnaMatchBlock *]
 *
 * Return [UNKN ]  Undocumented return value [DnaMatchBlock *]
 *
 */
DnaMatchBlock * Wise2_free_DnaMatchBlock(DnaMatchBlock * obj);
#define free_DnaMatchBlock Wise2_free_DnaMatchBlock


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
int Wise2_DnaMatchBlock_shatter_access_main(DnaMatchBlock * mat,int i,int j,int state);
#define DnaMatchBlock_shatter_access_main Wise2_DnaMatchBlock_shatter_access_main
int Wise2_DnaMatchBlock_shatter_access_special(DnaMatchBlock * mat,int i,int j,int state);
#define DnaMatchBlock_shatter_access_special Wise2_DnaMatchBlock_shatter_access_special
int Wise2_score_only_DnaMatchBlock(ComplexSequence* query,ComplexSequence* target ,DnaMatrix* comp65,DnaMatrix* comp75,DnaMatrix* comp85,DnaMatrix* comp95,Score g,Score u,Score v,Score s,Score b);
#define score_only_DnaMatchBlock Wise2_score_only_DnaMatchBlock
DnaMatchBlock * Wise2_allocate_DnaMatchBlock_only(ComplexSequence* query,ComplexSequence* target ,DnaMatrix* comp65,DnaMatrix* comp75,DnaMatrix* comp85,DnaMatrix* comp95,Score g,Score u,Score v,Score s,Score b);
#define allocate_DnaMatchBlock_only Wise2_allocate_DnaMatchBlock_only
void Wise2_init_DnaMatchBlock(DnaMatchBlock * mat);
#define init_DnaMatchBlock Wise2_init_DnaMatchBlock
AlnRange * Wise2_AlnRange_build_DnaMatchBlock(DnaMatchBlock * mat,int stopj,int stopspecstate,int * startj,int * startspecstate);
#define AlnRange_build_DnaMatchBlock Wise2_AlnRange_build_DnaMatchBlock
boolean Wise2_read_hidden_DnaMatchBlock(DnaMatchBlock * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out);
#define read_hidden_DnaMatchBlock Wise2_read_hidden_DnaMatchBlock
int Wise2_max_hidden_DnaMatchBlock(DnaMatchBlock * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_hidden_DnaMatchBlock Wise2_max_hidden_DnaMatchBlock
boolean Wise2_read_special_strip_DnaMatchBlock(DnaMatchBlock * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out);
#define read_special_strip_DnaMatchBlock Wise2_read_special_strip_DnaMatchBlock
int Wise2_max_special_strip_DnaMatchBlock(DnaMatchBlock * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_special_strip_DnaMatchBlock Wise2_max_special_strip_DnaMatchBlock
int Wise2_max_matrix_to_special_DnaMatchBlock(DnaMatchBlock * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_matrix_to_special_DnaMatchBlock Wise2_max_matrix_to_special_DnaMatchBlock
void Wise2_calculate_hidden_DnaMatchBlock(DnaMatchBlock * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv);
#define calculate_hidden_DnaMatchBlock Wise2_calculate_hidden_DnaMatchBlock
void Wise2_init_hidden_DnaMatchBlock(DnaMatchBlock * mat,int starti,int startj,int stopi,int stopj);
#define init_hidden_DnaMatchBlock Wise2_init_hidden_DnaMatchBlock
boolean Wise2_full_dc_DnaMatchBlock(DnaMatchBlock * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv);
#define full_dc_DnaMatchBlock Wise2_full_dc_DnaMatchBlock
boolean Wise2_do_dc_single_pass_DnaMatchBlock(DnaMatchBlock * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done);
#define do_dc_single_pass_DnaMatchBlock Wise2_do_dc_single_pass_DnaMatchBlock
void Wise2_push_dc_at_merge_DnaMatchBlock(DnaMatchBlock * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv);
#define push_dc_at_merge_DnaMatchBlock Wise2_push_dc_at_merge_DnaMatchBlock
void Wise2_follow_on_dc_DnaMatchBlock(DnaMatchBlock * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define follow_on_dc_DnaMatchBlock Wise2_follow_on_dc_DnaMatchBlock
void Wise2_run_up_dc_DnaMatchBlock(DnaMatchBlock * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define run_up_dc_DnaMatchBlock Wise2_run_up_dc_DnaMatchBlock
void Wise2_init_dc_DnaMatchBlock(DnaMatchBlock * mat);
#define init_dc_DnaMatchBlock Wise2_init_dc_DnaMatchBlock
int Wise2_start_end_find_end_DnaMatchBlock(DnaMatchBlock * mat,int * endj);
#define start_end_find_end_DnaMatchBlock Wise2_start_end_find_end_DnaMatchBlock
boolean Wise2_dc_optimised_start_end_calc_DnaMatchBlock(DnaMatchBlock *mat,DPEnvelope * dpenv);
#define dc_optimised_start_end_calc_DnaMatchBlock Wise2_dc_optimised_start_end_calc_DnaMatchBlock
void Wise2_init_start_end_linear_DnaMatchBlock(DnaMatchBlock * mat);
#define init_start_end_linear_DnaMatchBlock Wise2_init_start_end_linear_DnaMatchBlock
AlnConvertSet * Wise2_AlnConvertSet_DnaMatchBlock(void);
#define AlnConvertSet_DnaMatchBlock Wise2_AlnConvertSet_DnaMatchBlock
int Wise2_DnaMatchBlock_explicit_access_main(DnaMatchBlock * mat,int i,int j,int state);
#define DnaMatchBlock_explicit_access_main Wise2_DnaMatchBlock_explicit_access_main
int Wise2_DnaMatchBlock_explicit_access_special(DnaMatchBlock * mat,int i,int j,int state);
#define DnaMatchBlock_explicit_access_special Wise2_DnaMatchBlock_explicit_access_special
int Wise2_find_end_DnaMatchBlock(DnaMatchBlock * mat,int * ri,int * rj,int * state,boolean * isspecial,DnaMatchBlock_access_func_holder h);
#define find_end_DnaMatchBlock Wise2_find_end_DnaMatchBlock
void Wise2_DnaMatchBlock_debug_show_matrix(DnaMatchBlock * mat,int starti,int stopi,int startj,int stopj,FILE * ofp);
#define DnaMatchBlock_debug_show_matrix Wise2_DnaMatchBlock_debug_show_matrix
int Wise2_max_calc_DnaMatchBlock(DnaMatchBlock * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,DnaMatchBlock_access_func_holder h);
#define max_calc_DnaMatchBlock Wise2_max_calc_DnaMatchBlock
int Wise2_max_calc_special_DnaMatchBlock(DnaMatchBlock * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,DnaMatchBlock_access_func_holder h);
#define max_calc_special_DnaMatchBlock Wise2_max_calc_special_DnaMatchBlock

#ifdef _cplusplus
}
#endif

#endif
