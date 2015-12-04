#ifndef DYNAMITElargeblockdpHEADERFILE
#define DYNAMITElargeblockdpHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"


struct Wise2_LargeBlockAligner {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    BaseMatrix * basematrix;     
    ShatterMatrix * shatter;     
    int leni;    
    int lenj;    
    ComplexSequence* q;  
    ComplexSequence* t;  
    DnaMatrix* dm;   
    Score real_ext;  
    Score block_open;    
    Score un_dual;   
    Score un_single;     
    Score gap_open;  
    Score gap_ext;   
    } ;  
/* LargeBlockAligner defined */ 
#ifndef DYNAMITE_DEFINED_LargeBlockAligner
typedef struct Wise2_LargeBlockAligner Wise2_LargeBlockAligner;
#define LargeBlockAligner Wise2_LargeBlockAligner
#define DYNAMITE_DEFINED_LargeBlockAligner
#endif


struct Wise2_LargeBlockAligner_access_func_holder {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int (*access_main)(LargeBlockAligner*,int,int,int);  
    int (*access_special)(LargeBlockAligner*,int,int,int);   
    } ;  
/* LargeBlockAligner_access_func_holder defined */ 
#ifndef DYNAMITE_DEFINED_LargeBlockAligner_access_func_holder
typedef struct Wise2_LargeBlockAligner_access_func_holder Wise2_LargeBlockAligner_access_func_holder;
#define LargeBlockAligner_access_func_holder Wise2_LargeBlockAligner_access_func_holder
#define DYNAMITE_DEFINED_LargeBlockAligner_access_func_holder
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  PackAln_read_Shatter_LargeBlockAligner(mat)
 *
 * Descrip:    Reads off PackAln from shatter matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [LargeBlockAligner *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Shatter_LargeBlockAligner(LargeBlockAligner * mat);
#define PackAln_read_Shatter_LargeBlockAligner Wise2_PackAln_read_Shatter_LargeBlockAligner


/* Function:  calculate_shatter_LargeBlockAligner(mat,dpenv)
 *
 * Descrip:    This function calculates the LargeBlockAligner matrix when in shatter mode
 *
 *
 * Arg:          mat [UNKN ] (null) [LargeBlockAligner *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_shatter_LargeBlockAligner(LargeBlockAligner * mat,DPEnvelope * dpenv);
#define calculate_shatter_LargeBlockAligner Wise2_calculate_shatter_LargeBlockAligner


/* Function:  search_LargeBlockAligner(dbsi,out,q,t,dm,real_ext,block_open,un_dual,un_single,gap_open,gap_ext)
 *
 * Descrip:    This function makes a database search of LargeBlockAligner
 *             It uses the dbsi structure to choose which implementation
 *             to use of the database searching. This way at run time you
 *             can switch between single threaded/multi-threaded or hardware
 *
 *
 * Arg:              dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 * Arg:               out [UNKN ] Undocumented argument [Hscore *]
 * Arg:                 q [UNKN ] Undocumented argument [ComplexSequence*]
 * Arg:                 t [UNKN ] Undocumented argument [ComplexSequence*]
 * Arg:                dm [UNKN ] Undocumented argument [DnaMatrix*]
 * Arg:          real_ext [UNKN ] Undocumented argument [Score]
 * Arg:        block_open [UNKN ] Undocumented argument [Score]
 * Arg:           un_dual [UNKN ] Undocumented argument [Score]
 * Arg:         un_single [UNKN ] Undocumented argument [Score]
 * Arg:          gap_open [UNKN ] Undocumented argument [Score]
 * Arg:           gap_ext [UNKN ] Undocumented argument [Score]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_search_LargeBlockAligner(DBSearchImpl * dbsi,Hscore * out,ComplexSequence* q,ComplexSequence* t ,DnaMatrix* dm,Score real_ext,Score block_open,Score un_dual,Score un_single,Score gap_open,Score gap_ext);
#define search_LargeBlockAligner Wise2_search_LargeBlockAligner


/* Function:  serial_search_LargeBlockAligner(out,q,t,dm,real_ext,block_open,un_dual,un_single,gap_open,gap_ext)
 *
 * Descrip:    This function makes a database search of LargeBlockAligner
 *             It is a single processor implementation
 *
 *
 * Arg:               out [UNKN ] Undocumented argument [Hscore *]
 * Arg:                 q [UNKN ] Undocumented argument [ComplexSequence*]
 * Arg:                 t [UNKN ] Undocumented argument [ComplexSequence*]
 * Arg:                dm [UNKN ] Undocumented argument [DnaMatrix*]
 * Arg:          real_ext [UNKN ] Undocumented argument [Score]
 * Arg:        block_open [UNKN ] Undocumented argument [Score]
 * Arg:           un_dual [UNKN ] Undocumented argument [Score]
 * Arg:         un_single [UNKN ] Undocumented argument [Score]
 * Arg:          gap_open [UNKN ] Undocumented argument [Score]
 * Arg:           gap_ext [UNKN ] Undocumented argument [Score]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_serial_search_LargeBlockAligner(Hscore * out,ComplexSequence* q,ComplexSequence* t ,DnaMatrix* dm,Score real_ext,Score block_open,Score un_dual,Score un_single,Score gap_open,Score gap_ext);
#define serial_search_LargeBlockAligner Wise2_serial_search_LargeBlockAligner


/* Function:  PackAln_bestmemory_LargeBlockAligner(q,t,dm,real_ext,block_open,un_dual,un_single,gap_open,gap_ext,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_LargeBlockAligner
 *
 *
 * Arg:                 q [UNKN ] query data structure [ComplexSequence*]
 * Arg:                 t [UNKN ] target data structure [ComplexSequence*]
 * Arg:                dm [UNKN ] Resource [DnaMatrix*]
 * Arg:          real_ext [UNKN ] Resource [Score]
 * Arg:        block_open [UNKN ] Resource [Score]
 * Arg:           un_dual [UNKN ] Resource [Score]
 * Arg:         un_single [UNKN ] Resource [Score]
 * Arg:          gap_open [UNKN ] Resource [Score]
 * Arg:           gap_ext [UNKN ] Resource [Score]
 * Arg:             dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:              dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_bestmemory_LargeBlockAligner(ComplexSequence* q,ComplexSequence* t ,DnaMatrix* dm,Score real_ext,Score block_open,Score un_dual,Score un_single,Score gap_open,Score gap_ext,DPEnvelope * dpenv,DPRunImpl * dpri);
#define PackAln_bestmemory_LargeBlockAligner Wise2_PackAln_bestmemory_LargeBlockAligner


/* Function:  allocate_Expl_LargeBlockAligner(q,t,dm,real_ext,block_open,un_dual,un_single,gap_open,gap_ext,dpri)
 *
 * Descrip:    This function allocates the LargeBlockAligner structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_LargeBlockAligner_only
 *
 *
 * Arg:                 q [UNKN ] query data structure [ComplexSequence*]
 * Arg:                 t [UNKN ] target data structure [ComplexSequence*]
 * Arg:                dm [UNKN ] Resource [DnaMatrix*]
 * Arg:          real_ext [UNKN ] Resource [Score]
 * Arg:        block_open [UNKN ] Resource [Score]
 * Arg:           un_dual [UNKN ] Resource [Score]
 * Arg:         un_single [UNKN ] Resource [Score]
 * Arg:          gap_open [UNKN ] Resource [Score]
 * Arg:           gap_ext [UNKN ] Resource [Score]
 * Arg:              dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [LargeBlockAligner *]
 *
 */
LargeBlockAligner * Wise2_allocate_Expl_LargeBlockAligner(ComplexSequence* q,ComplexSequence* t ,DnaMatrix* dm,Score real_ext,Score block_open,Score un_dual,Score un_single,Score gap_open,Score gap_ext,DPRunImpl * dpri);
#define allocate_Expl_LargeBlockAligner Wise2_allocate_Expl_LargeBlockAligner


/* Function:  recalculate_PackAln_LargeBlockAligner(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by LargeBlockAligner
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [LargeBlockAligner *]
 *
 */
void Wise2_recalculate_PackAln_LargeBlockAligner(PackAln * pal,LargeBlockAligner * mat);
#define recalculate_PackAln_LargeBlockAligner Wise2_recalculate_PackAln_LargeBlockAligner


/* Function:  allocate_Small_LargeBlockAligner(q,t,dm,real_ext,block_open,un_dual,un_single,gap_open,gap_ext)
 *
 * Descrip:    This function allocates the LargeBlockAligner structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_LargeBlockAligner_only
 *
 *
 * Arg:                 q [UNKN ] query data structure [ComplexSequence*]
 * Arg:                 t [UNKN ] target data structure [ComplexSequence*]
 * Arg:                dm [UNKN ] Resource [DnaMatrix*]
 * Arg:          real_ext [UNKN ] Resource [Score]
 * Arg:        block_open [UNKN ] Resource [Score]
 * Arg:           un_dual [UNKN ] Resource [Score]
 * Arg:         un_single [UNKN ] Resource [Score]
 * Arg:          gap_open [UNKN ] Resource [Score]
 * Arg:           gap_ext [UNKN ] Resource [Score]
 *
 * Return [UNKN ]  Undocumented return value [LargeBlockAligner *]
 *
 */
LargeBlockAligner * Wise2_allocate_Small_LargeBlockAligner(ComplexSequence* q,ComplexSequence* t ,DnaMatrix* dm,Score real_ext,Score block_open,Score un_dual,Score un_single,Score gap_open,Score gap_ext);
#define allocate_Small_LargeBlockAligner Wise2_allocate_Small_LargeBlockAligner


/* Function:  PackAln_calculate_Small_LargeBlockAligner(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for LargeBlockAligner structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_LargeBlockAligner 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_LargeBlockAligner 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [LargeBlockAligner *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_calculate_Small_LargeBlockAligner(LargeBlockAligner * mat,DPEnvelope * dpenv);
#define PackAln_calculate_Small_LargeBlockAligner Wise2_PackAln_calculate_Small_LargeBlockAligner


/* Function:  AlnRangeSet_calculate_Small_LargeBlockAligner(mat)
 *
 * Descrip:    This function calculates an alignment for LargeBlockAligner structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_LargeBlockAligner 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_LargeBlockAligner
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_LargeBlockAligner 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [LargeBlockAligner *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_calculate_Small_LargeBlockAligner(LargeBlockAligner * mat);
#define AlnRangeSet_calculate_Small_LargeBlockAligner Wise2_AlnRangeSet_calculate_Small_LargeBlockAligner


/* Function:  AlnRangeSet_from_LargeBlockAligner(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for LargeBlockAligner structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_LargeBlockAligner 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_LargeBlockAligner
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [LargeBlockAligner *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_from_LargeBlockAligner(LargeBlockAligner * mat);
#define AlnRangeSet_from_LargeBlockAligner Wise2_AlnRangeSet_from_LargeBlockAligner


/* Function:  convert_PackAln_to_AlnBlock_LargeBlockAligner(pal)
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
AlnBlock * Wise2_convert_PackAln_to_AlnBlock_LargeBlockAligner(PackAln * pal);
#define convert_PackAln_to_AlnBlock_LargeBlockAligner Wise2_convert_PackAln_to_AlnBlock_LargeBlockAligner


/* Function:  PackAln_read_Expl_LargeBlockAligner(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [LargeBlockAligner *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Expl_LargeBlockAligner(LargeBlockAligner * mat);
#define PackAln_read_Expl_LargeBlockAligner Wise2_PackAln_read_Expl_LargeBlockAligner


/* Function:  PackAln_read_generic_LargeBlockAligner(mat,h)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [LargeBlockAligner *]
 * Arg:          h [UNKN ] Undocumented argument [LargeBlockAligner_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_generic_LargeBlockAligner(LargeBlockAligner * mat,LargeBlockAligner_access_func_holder h);
#define PackAln_read_generic_LargeBlockAligner Wise2_PackAln_read_generic_LargeBlockAligner


/* Function:  calculate_LargeBlockAligner(mat)
 *
 * Descrip:    This function calculates the LargeBlockAligner matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_LargeBlockAligner
 *
 *
 * Arg:        mat [UNKN ] LargeBlockAligner which contains explicit basematrix memory [LargeBlockAligner *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_LargeBlockAligner(LargeBlockAligner * mat);
#define calculate_LargeBlockAligner Wise2_calculate_LargeBlockAligner


/* Function:  calculate_dpenv_LargeBlockAligner(mat,dpenv)
 *
 * Descrip:    This function calculates the LargeBlockAligner matrix when in explicit mode, subject to the envelope
 *
 *
 * Arg:          mat [UNKN ] LargeBlockAligner which contains explicit basematrix memory [LargeBlockAligner *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_dpenv_LargeBlockAligner(LargeBlockAligner * mat,DPEnvelope * dpenv);
#define calculate_dpenv_LargeBlockAligner Wise2_calculate_dpenv_LargeBlockAligner


/* Function:  LargeBlockAligner_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LargeBlockAligner *]
 *
 */
LargeBlockAligner * Wise2_LargeBlockAligner_alloc(void);
#define LargeBlockAligner_alloc Wise2_LargeBlockAligner_alloc


/* Function:  free_LargeBlockAligner(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [LargeBlockAligner *]
 *
 * Return [UNKN ]  Undocumented return value [LargeBlockAligner *]
 *
 */
LargeBlockAligner * Wise2_free_LargeBlockAligner(LargeBlockAligner * obj);
#define free_LargeBlockAligner Wise2_free_LargeBlockAligner


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
int Wise2_LargeBlockAligner_shatter_access_main(LargeBlockAligner * mat,int i,int j,int state);
#define LargeBlockAligner_shatter_access_main Wise2_LargeBlockAligner_shatter_access_main
int Wise2_LargeBlockAligner_shatter_access_special(LargeBlockAligner * mat,int i,int j,int state);
#define LargeBlockAligner_shatter_access_special Wise2_LargeBlockAligner_shatter_access_special
int Wise2_score_only_LargeBlockAligner(ComplexSequence* q,ComplexSequence* t ,DnaMatrix* dm,Score real_ext,Score block_open,Score un_dual,Score un_single,Score gap_open,Score gap_ext);
#define score_only_LargeBlockAligner Wise2_score_only_LargeBlockAligner
LargeBlockAligner * Wise2_allocate_LargeBlockAligner_only(ComplexSequence* q,ComplexSequence* t ,DnaMatrix* dm,Score real_ext,Score block_open,Score un_dual,Score un_single,Score gap_open,Score gap_ext);
#define allocate_LargeBlockAligner_only Wise2_allocate_LargeBlockAligner_only
void Wise2_init_LargeBlockAligner(LargeBlockAligner * mat);
#define init_LargeBlockAligner Wise2_init_LargeBlockAligner
AlnRange * Wise2_AlnRange_build_LargeBlockAligner(LargeBlockAligner * mat,int stopj,int stopspecstate,int * startj,int * startspecstate);
#define AlnRange_build_LargeBlockAligner Wise2_AlnRange_build_LargeBlockAligner
boolean Wise2_read_hidden_LargeBlockAligner(LargeBlockAligner * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out);
#define read_hidden_LargeBlockAligner Wise2_read_hidden_LargeBlockAligner
int Wise2_max_hidden_LargeBlockAligner(LargeBlockAligner * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_hidden_LargeBlockAligner Wise2_max_hidden_LargeBlockAligner
boolean Wise2_read_special_strip_LargeBlockAligner(LargeBlockAligner * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out);
#define read_special_strip_LargeBlockAligner Wise2_read_special_strip_LargeBlockAligner
int Wise2_max_special_strip_LargeBlockAligner(LargeBlockAligner * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_special_strip_LargeBlockAligner Wise2_max_special_strip_LargeBlockAligner
int Wise2_max_matrix_to_special_LargeBlockAligner(LargeBlockAligner * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_matrix_to_special_LargeBlockAligner Wise2_max_matrix_to_special_LargeBlockAligner
void Wise2_calculate_hidden_LargeBlockAligner(LargeBlockAligner * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv);
#define calculate_hidden_LargeBlockAligner Wise2_calculate_hidden_LargeBlockAligner
void Wise2_init_hidden_LargeBlockAligner(LargeBlockAligner * mat,int starti,int startj,int stopi,int stopj);
#define init_hidden_LargeBlockAligner Wise2_init_hidden_LargeBlockAligner
boolean Wise2_full_dc_LargeBlockAligner(LargeBlockAligner * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv);
#define full_dc_LargeBlockAligner Wise2_full_dc_LargeBlockAligner
boolean Wise2_do_dc_single_pass_LargeBlockAligner(LargeBlockAligner * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done);
#define do_dc_single_pass_LargeBlockAligner Wise2_do_dc_single_pass_LargeBlockAligner
void Wise2_push_dc_at_merge_LargeBlockAligner(LargeBlockAligner * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv);
#define push_dc_at_merge_LargeBlockAligner Wise2_push_dc_at_merge_LargeBlockAligner
void Wise2_follow_on_dc_LargeBlockAligner(LargeBlockAligner * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define follow_on_dc_LargeBlockAligner Wise2_follow_on_dc_LargeBlockAligner
void Wise2_run_up_dc_LargeBlockAligner(LargeBlockAligner * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define run_up_dc_LargeBlockAligner Wise2_run_up_dc_LargeBlockAligner
void Wise2_init_dc_LargeBlockAligner(LargeBlockAligner * mat);
#define init_dc_LargeBlockAligner Wise2_init_dc_LargeBlockAligner
int Wise2_start_end_find_end_LargeBlockAligner(LargeBlockAligner * mat,int * endj);
#define start_end_find_end_LargeBlockAligner Wise2_start_end_find_end_LargeBlockAligner
boolean Wise2_dc_optimised_start_end_calc_LargeBlockAligner(LargeBlockAligner *mat,DPEnvelope * dpenv);
#define dc_optimised_start_end_calc_LargeBlockAligner Wise2_dc_optimised_start_end_calc_LargeBlockAligner
void Wise2_init_start_end_linear_LargeBlockAligner(LargeBlockAligner * mat);
#define init_start_end_linear_LargeBlockAligner Wise2_init_start_end_linear_LargeBlockAligner
AlnConvertSet * Wise2_AlnConvertSet_LargeBlockAligner(void);
#define AlnConvertSet_LargeBlockAligner Wise2_AlnConvertSet_LargeBlockAligner
int Wise2_LargeBlockAligner_explicit_access_main(LargeBlockAligner * mat,int i,int j,int state);
#define LargeBlockAligner_explicit_access_main Wise2_LargeBlockAligner_explicit_access_main
int Wise2_LargeBlockAligner_explicit_access_special(LargeBlockAligner * mat,int i,int j,int state);
#define LargeBlockAligner_explicit_access_special Wise2_LargeBlockAligner_explicit_access_special
int Wise2_find_end_LargeBlockAligner(LargeBlockAligner * mat,int * ri,int * rj,int * state,boolean * isspecial,LargeBlockAligner_access_func_holder h);
#define find_end_LargeBlockAligner Wise2_find_end_LargeBlockAligner
void Wise2_LargeBlockAligner_debug_show_matrix(LargeBlockAligner * mat,int starti,int stopi,int startj,int stopj,FILE * ofp);
#define LargeBlockAligner_debug_show_matrix Wise2_LargeBlockAligner_debug_show_matrix
int Wise2_max_calc_LargeBlockAligner(LargeBlockAligner * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,LargeBlockAligner_access_func_holder h);
#define max_calc_LargeBlockAligner Wise2_max_calc_LargeBlockAligner
int Wise2_max_calc_special_LargeBlockAligner(LargeBlockAligner * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,LargeBlockAligner_access_func_holder h);
#define max_calc_special_LargeBlockAligner Wise2_max_calc_special_LargeBlockAligner

#ifdef _cplusplus
}
#endif

#endif
