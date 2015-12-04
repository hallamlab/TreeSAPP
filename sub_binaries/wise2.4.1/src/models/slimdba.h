#ifndef DYNAMITEslimdbaHEADERFILE
#define DYNAMITEslimdbaHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"

struct Wise2_SlimDnaMatchBlock {  
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
    Score g;     
    Score u;     
    Score v;     
    Score s;     
    Score b;     
    } ;  
/* SlimDnaMatchBlock defined */ 
#ifndef DYNAMITE_DEFINED_SlimDnaMatchBlock
typedef struct Wise2_SlimDnaMatchBlock Wise2_SlimDnaMatchBlock;
#define SlimDnaMatchBlock Wise2_SlimDnaMatchBlock
#define DYNAMITE_DEFINED_SlimDnaMatchBlock
#endif


struct Wise2_SlimDnaMatchBlock_access_func_holder {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int (*access_main)(SlimDnaMatchBlock*,int,int,int);  
    int (*access_special)(SlimDnaMatchBlock*,int,int,int);   
    } ;  
/* SlimDnaMatchBlock_access_func_holder defined */ 
#ifndef DYNAMITE_DEFINED_SlimDnaMatchBlock_access_func_holder
typedef struct Wise2_SlimDnaMatchBlock_access_func_holder Wise2_SlimDnaMatchBlock_access_func_holder;
#define SlimDnaMatchBlock_access_func_holder Wise2_SlimDnaMatchBlock_access_func_holder
#define DYNAMITE_DEFINED_SlimDnaMatchBlock_access_func_holder
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  PackAln_read_Shatter_SlimDnaMatchBlock(mat)
 *
 * Descrip:    Reads off PackAln from shatter matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [SlimDnaMatchBlock *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Shatter_SlimDnaMatchBlock(SlimDnaMatchBlock * mat);
#define PackAln_read_Shatter_SlimDnaMatchBlock Wise2_PackAln_read_Shatter_SlimDnaMatchBlock


/* Function:  calculate_shatter_SlimDnaMatchBlock(mat,dpenv)
 *
 * Descrip:    This function calculates the SlimDnaMatchBlock matrix when in shatter mode
 *
 *
 * Arg:          mat [UNKN ] (null) [SlimDnaMatchBlock *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_shatter_SlimDnaMatchBlock(SlimDnaMatchBlock * mat,DPEnvelope * dpenv);
#define calculate_shatter_SlimDnaMatchBlock Wise2_calculate_shatter_SlimDnaMatchBlock


/* Function:  search_SlimDnaMatchBlock(dbsi,out,query,target,comp65,g,u,v,s,b)
 *
 * Descrip:    This function makes a database search of SlimDnaMatchBlock
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
 * Arg:             g [UNKN ] Undocumented argument [Score]
 * Arg:             u [UNKN ] Undocumented argument [Score]
 * Arg:             v [UNKN ] Undocumented argument [Score]
 * Arg:             s [UNKN ] Undocumented argument [Score]
 * Arg:             b [UNKN ] Undocumented argument [Score]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_search_SlimDnaMatchBlock(DBSearchImpl * dbsi,Hscore * out,ComplexSequence* query,ComplexSequence* target ,DnaMatrix* comp65,Score g,Score u,Score v,Score s,Score b);
#define search_SlimDnaMatchBlock Wise2_search_SlimDnaMatchBlock


/* Function:  serial_search_SlimDnaMatchBlock(out,query,target,comp65,g,u,v,s,b)
 *
 * Descrip:    This function makes a database search of SlimDnaMatchBlock
 *             It is a single processor implementation
 *
 *
 * Arg:           out [UNKN ] Undocumented argument [Hscore *]
 * Arg:         query [UNKN ] Undocumented argument [ComplexSequence*]
 * Arg:        target [UNKN ] Undocumented argument [ComplexSequence*]
 * Arg:        comp65 [UNKN ] Undocumented argument [DnaMatrix*]
 * Arg:             g [UNKN ] Undocumented argument [Score]
 * Arg:             u [UNKN ] Undocumented argument [Score]
 * Arg:             v [UNKN ] Undocumented argument [Score]
 * Arg:             s [UNKN ] Undocumented argument [Score]
 * Arg:             b [UNKN ] Undocumented argument [Score]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_serial_search_SlimDnaMatchBlock(Hscore * out,ComplexSequence* query,ComplexSequence* target ,DnaMatrix* comp65,Score g,Score u,Score v,Score s,Score b);
#define serial_search_SlimDnaMatchBlock Wise2_serial_search_SlimDnaMatchBlock


/* Function:  PackAln_bestmemory_SlimDnaMatchBlock(query,target,comp65,g,u,v,s,b,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_SlimDnaMatchBlock
 *
 *
 * Arg:         query [UNKN ] query data structure [ComplexSequence*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:        comp65 [UNKN ] Resource [DnaMatrix*]
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
PackAln * Wise2_PackAln_bestmemory_SlimDnaMatchBlock(ComplexSequence* query,ComplexSequence* target ,DnaMatrix* comp65,Score g,Score u,Score v,Score s,Score b,DPEnvelope * dpenv,DPRunImpl * dpri);
#define PackAln_bestmemory_SlimDnaMatchBlock Wise2_PackAln_bestmemory_SlimDnaMatchBlock


/* Function:  allocate_Expl_SlimDnaMatchBlock(query,target,comp65,g,u,v,s,b,dpri)
 *
 * Descrip:    This function allocates the SlimDnaMatchBlock structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_SlimDnaMatchBlock_only
 *
 *
 * Arg:         query [UNKN ] query data structure [ComplexSequence*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:        comp65 [UNKN ] Resource [DnaMatrix*]
 * Arg:             g [UNKN ] Resource [Score]
 * Arg:             u [UNKN ] Resource [Score]
 * Arg:             v [UNKN ] Resource [Score]
 * Arg:             s [UNKN ] Resource [Score]
 * Arg:             b [UNKN ] Resource [Score]
 * Arg:          dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [SlimDnaMatchBlock *]
 *
 */
SlimDnaMatchBlock * Wise2_allocate_Expl_SlimDnaMatchBlock(ComplexSequence* query,ComplexSequence* target ,DnaMatrix* comp65,Score g,Score u,Score v,Score s,Score b,DPRunImpl * dpri);
#define allocate_Expl_SlimDnaMatchBlock Wise2_allocate_Expl_SlimDnaMatchBlock


/* Function:  recalculate_PackAln_SlimDnaMatchBlock(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by SlimDnaMatchBlock
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [SlimDnaMatchBlock *]
 *
 */
void Wise2_recalculate_PackAln_SlimDnaMatchBlock(PackAln * pal,SlimDnaMatchBlock * mat);
#define recalculate_PackAln_SlimDnaMatchBlock Wise2_recalculate_PackAln_SlimDnaMatchBlock


/* Function:  allocate_Small_SlimDnaMatchBlock(query,target,comp65,g,u,v,s,b)
 *
 * Descrip:    This function allocates the SlimDnaMatchBlock structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_SlimDnaMatchBlock_only
 *
 *
 * Arg:         query [UNKN ] query data structure [ComplexSequence*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:        comp65 [UNKN ] Resource [DnaMatrix*]
 * Arg:             g [UNKN ] Resource [Score]
 * Arg:             u [UNKN ] Resource [Score]
 * Arg:             v [UNKN ] Resource [Score]
 * Arg:             s [UNKN ] Resource [Score]
 * Arg:             b [UNKN ] Resource [Score]
 *
 * Return [UNKN ]  Undocumented return value [SlimDnaMatchBlock *]
 *
 */
SlimDnaMatchBlock * Wise2_allocate_Small_SlimDnaMatchBlock(ComplexSequence* query,ComplexSequence* target ,DnaMatrix* comp65,Score g,Score u,Score v,Score s,Score b);
#define allocate_Small_SlimDnaMatchBlock Wise2_allocate_Small_SlimDnaMatchBlock


/* Function:  PackAln_calculate_Small_SlimDnaMatchBlock(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for SlimDnaMatchBlock structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_SlimDnaMatchBlock 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_SlimDnaMatchBlock 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [SlimDnaMatchBlock *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_calculate_Small_SlimDnaMatchBlock(SlimDnaMatchBlock * mat,DPEnvelope * dpenv);
#define PackAln_calculate_Small_SlimDnaMatchBlock Wise2_PackAln_calculate_Small_SlimDnaMatchBlock


/* Function:  AlnRangeSet_calculate_Small_SlimDnaMatchBlock(mat)
 *
 * Descrip:    This function calculates an alignment for SlimDnaMatchBlock structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_SlimDnaMatchBlock 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_SlimDnaMatchBlock
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_SlimDnaMatchBlock 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [SlimDnaMatchBlock *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_calculate_Small_SlimDnaMatchBlock(SlimDnaMatchBlock * mat);
#define AlnRangeSet_calculate_Small_SlimDnaMatchBlock Wise2_AlnRangeSet_calculate_Small_SlimDnaMatchBlock


/* Function:  AlnRangeSet_from_SlimDnaMatchBlock(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for SlimDnaMatchBlock structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_SlimDnaMatchBlock 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_SlimDnaMatchBlock
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [SlimDnaMatchBlock *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_from_SlimDnaMatchBlock(SlimDnaMatchBlock * mat);
#define AlnRangeSet_from_SlimDnaMatchBlock Wise2_AlnRangeSet_from_SlimDnaMatchBlock


/* Function:  convert_PackAln_to_AlnBlock_SlimDnaMatchBlock(pal)
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
AlnBlock * Wise2_convert_PackAln_to_AlnBlock_SlimDnaMatchBlock(PackAln * pal);
#define convert_PackAln_to_AlnBlock_SlimDnaMatchBlock Wise2_convert_PackAln_to_AlnBlock_SlimDnaMatchBlock


/* Function:  PackAln_read_Expl_SlimDnaMatchBlock(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [SlimDnaMatchBlock *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Expl_SlimDnaMatchBlock(SlimDnaMatchBlock * mat);
#define PackAln_read_Expl_SlimDnaMatchBlock Wise2_PackAln_read_Expl_SlimDnaMatchBlock


/* Function:  PackAln_read_generic_SlimDnaMatchBlock(mat,h)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [SlimDnaMatchBlock *]
 * Arg:          h [UNKN ] Undocumented argument [SlimDnaMatchBlock_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_generic_SlimDnaMatchBlock(SlimDnaMatchBlock * mat,SlimDnaMatchBlock_access_func_holder h);
#define PackAln_read_generic_SlimDnaMatchBlock Wise2_PackAln_read_generic_SlimDnaMatchBlock


/* Function:  calculate_SlimDnaMatchBlock(mat)
 *
 * Descrip:    This function calculates the SlimDnaMatchBlock matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_SlimDnaMatchBlock
 *
 *
 * Arg:        mat [UNKN ] SlimDnaMatchBlock which contains explicit basematrix memory [SlimDnaMatchBlock *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_SlimDnaMatchBlock(SlimDnaMatchBlock * mat);
#define calculate_SlimDnaMatchBlock Wise2_calculate_SlimDnaMatchBlock


/* Function:  calculate_dpenv_SlimDnaMatchBlock(mat,dpenv)
 *
 * Descrip:    This function calculates the SlimDnaMatchBlock matrix when in explicit mode, subject to the envelope
 *
 *
 * Arg:          mat [UNKN ] SlimDnaMatchBlock which contains explicit basematrix memory [SlimDnaMatchBlock *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_dpenv_SlimDnaMatchBlock(SlimDnaMatchBlock * mat,DPEnvelope * dpenv);
#define calculate_dpenv_SlimDnaMatchBlock Wise2_calculate_dpenv_SlimDnaMatchBlock


/* Function:  SlimDnaMatchBlock_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SlimDnaMatchBlock *]
 *
 */
SlimDnaMatchBlock * Wise2_SlimDnaMatchBlock_alloc(void);
#define SlimDnaMatchBlock_alloc Wise2_SlimDnaMatchBlock_alloc


/* Function:  free_SlimDnaMatchBlock(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SlimDnaMatchBlock *]
 *
 * Return [UNKN ]  Undocumented return value [SlimDnaMatchBlock *]
 *
 */
SlimDnaMatchBlock * Wise2_free_SlimDnaMatchBlock(SlimDnaMatchBlock * obj);
#define free_SlimDnaMatchBlock Wise2_free_SlimDnaMatchBlock


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
int Wise2_SlimDnaMatchBlock_shatter_access_main(SlimDnaMatchBlock * mat,int i,int j,int state);
#define SlimDnaMatchBlock_shatter_access_main Wise2_SlimDnaMatchBlock_shatter_access_main
int Wise2_SlimDnaMatchBlock_shatter_access_special(SlimDnaMatchBlock * mat,int i,int j,int state);
#define SlimDnaMatchBlock_shatter_access_special Wise2_SlimDnaMatchBlock_shatter_access_special
int Wise2_score_only_SlimDnaMatchBlock(ComplexSequence* query,ComplexSequence* target ,DnaMatrix* comp65,Score g,Score u,Score v,Score s,Score b);
#define score_only_SlimDnaMatchBlock Wise2_score_only_SlimDnaMatchBlock
SlimDnaMatchBlock * Wise2_allocate_SlimDnaMatchBlock_only(ComplexSequence* query,ComplexSequence* target ,DnaMatrix* comp65,Score g,Score u,Score v,Score s,Score b);
#define allocate_SlimDnaMatchBlock_only Wise2_allocate_SlimDnaMatchBlock_only
void Wise2_init_SlimDnaMatchBlock(SlimDnaMatchBlock * mat);
#define init_SlimDnaMatchBlock Wise2_init_SlimDnaMatchBlock
AlnRange * Wise2_AlnRange_build_SlimDnaMatchBlock(SlimDnaMatchBlock * mat,int stopj,int stopspecstate,int * startj,int * startspecstate);
#define AlnRange_build_SlimDnaMatchBlock Wise2_AlnRange_build_SlimDnaMatchBlock
boolean Wise2_read_hidden_SlimDnaMatchBlock(SlimDnaMatchBlock * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out);
#define read_hidden_SlimDnaMatchBlock Wise2_read_hidden_SlimDnaMatchBlock
int Wise2_max_hidden_SlimDnaMatchBlock(SlimDnaMatchBlock * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_hidden_SlimDnaMatchBlock Wise2_max_hidden_SlimDnaMatchBlock
boolean Wise2_read_special_strip_SlimDnaMatchBlock(SlimDnaMatchBlock * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out);
#define read_special_strip_SlimDnaMatchBlock Wise2_read_special_strip_SlimDnaMatchBlock
int Wise2_max_special_strip_SlimDnaMatchBlock(SlimDnaMatchBlock * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_special_strip_SlimDnaMatchBlock Wise2_max_special_strip_SlimDnaMatchBlock
int Wise2_max_matrix_to_special_SlimDnaMatchBlock(SlimDnaMatchBlock * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_matrix_to_special_SlimDnaMatchBlock Wise2_max_matrix_to_special_SlimDnaMatchBlock
void Wise2_calculate_hidden_SlimDnaMatchBlock(SlimDnaMatchBlock * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv);
#define calculate_hidden_SlimDnaMatchBlock Wise2_calculate_hidden_SlimDnaMatchBlock
void Wise2_init_hidden_SlimDnaMatchBlock(SlimDnaMatchBlock * mat,int starti,int startj,int stopi,int stopj);
#define init_hidden_SlimDnaMatchBlock Wise2_init_hidden_SlimDnaMatchBlock
boolean Wise2_full_dc_SlimDnaMatchBlock(SlimDnaMatchBlock * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv);
#define full_dc_SlimDnaMatchBlock Wise2_full_dc_SlimDnaMatchBlock
boolean Wise2_do_dc_single_pass_SlimDnaMatchBlock(SlimDnaMatchBlock * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done);
#define do_dc_single_pass_SlimDnaMatchBlock Wise2_do_dc_single_pass_SlimDnaMatchBlock
void Wise2_push_dc_at_merge_SlimDnaMatchBlock(SlimDnaMatchBlock * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv);
#define push_dc_at_merge_SlimDnaMatchBlock Wise2_push_dc_at_merge_SlimDnaMatchBlock
void Wise2_follow_on_dc_SlimDnaMatchBlock(SlimDnaMatchBlock * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define follow_on_dc_SlimDnaMatchBlock Wise2_follow_on_dc_SlimDnaMatchBlock
void Wise2_run_up_dc_SlimDnaMatchBlock(SlimDnaMatchBlock * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define run_up_dc_SlimDnaMatchBlock Wise2_run_up_dc_SlimDnaMatchBlock
void Wise2_init_dc_SlimDnaMatchBlock(SlimDnaMatchBlock * mat);
#define init_dc_SlimDnaMatchBlock Wise2_init_dc_SlimDnaMatchBlock
int Wise2_start_end_find_end_SlimDnaMatchBlock(SlimDnaMatchBlock * mat,int * endj);
#define start_end_find_end_SlimDnaMatchBlock Wise2_start_end_find_end_SlimDnaMatchBlock
boolean Wise2_dc_optimised_start_end_calc_SlimDnaMatchBlock(SlimDnaMatchBlock *mat,DPEnvelope * dpenv);
#define dc_optimised_start_end_calc_SlimDnaMatchBlock Wise2_dc_optimised_start_end_calc_SlimDnaMatchBlock
void Wise2_init_start_end_linear_SlimDnaMatchBlock(SlimDnaMatchBlock * mat);
#define init_start_end_linear_SlimDnaMatchBlock Wise2_init_start_end_linear_SlimDnaMatchBlock
AlnConvertSet * Wise2_AlnConvertSet_SlimDnaMatchBlock(void);
#define AlnConvertSet_SlimDnaMatchBlock Wise2_AlnConvertSet_SlimDnaMatchBlock
int Wise2_SlimDnaMatchBlock_explicit_access_main(SlimDnaMatchBlock * mat,int i,int j,int state);
#define SlimDnaMatchBlock_explicit_access_main Wise2_SlimDnaMatchBlock_explicit_access_main
int Wise2_SlimDnaMatchBlock_explicit_access_special(SlimDnaMatchBlock * mat,int i,int j,int state);
#define SlimDnaMatchBlock_explicit_access_special Wise2_SlimDnaMatchBlock_explicit_access_special
int Wise2_find_end_SlimDnaMatchBlock(SlimDnaMatchBlock * mat,int * ri,int * rj,int * state,boolean * isspecial,SlimDnaMatchBlock_access_func_holder h);
#define find_end_SlimDnaMatchBlock Wise2_find_end_SlimDnaMatchBlock
void Wise2_SlimDnaMatchBlock_debug_show_matrix(SlimDnaMatchBlock * mat,int starti,int stopi,int startj,int stopj,FILE * ofp);
#define SlimDnaMatchBlock_debug_show_matrix Wise2_SlimDnaMatchBlock_debug_show_matrix
int Wise2_max_calc_SlimDnaMatchBlock(SlimDnaMatchBlock * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,SlimDnaMatchBlock_access_func_holder h);
#define max_calc_SlimDnaMatchBlock Wise2_max_calc_SlimDnaMatchBlock
int Wise2_max_calc_special_SlimDnaMatchBlock(SlimDnaMatchBlock * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,SlimDnaMatchBlock_access_func_holder h);
#define max_calc_special_SlimDnaMatchBlock Wise2_max_calc_special_SlimDnaMatchBlock

#ifdef _cplusplus
}
#endif

#endif
