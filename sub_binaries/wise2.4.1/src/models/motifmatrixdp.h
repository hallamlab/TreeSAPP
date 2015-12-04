#ifndef DYNAMITEmotifmatrixdpHEADERFILE
#define DYNAMITEmotifmatrixdpHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "motifmatrix.h"

#define MOTIF_MATRIX_IN(motif,i,j,score) (motif->mat[i][j] == 0 ? 0 : (score))

struct Wise2_LocalMotifMatrix {  
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
    MotifConsMatrix* motif;  
    MotifMatrixScore* ms;    
    } ;  
/* LocalMotifMatrix defined */ 
#ifndef DYNAMITE_DEFINED_LocalMotifMatrix
typedef struct Wise2_LocalMotifMatrix Wise2_LocalMotifMatrix;
#define LocalMotifMatrix Wise2_LocalMotifMatrix
#define DYNAMITE_DEFINED_LocalMotifMatrix
#endif


struct Wise2_LocalMotifMatrix_access_func_holder {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int (*access_main)(LocalMotifMatrix*,int,int,int);   
    int (*access_special)(LocalMotifMatrix*,int,int,int);    
    } ;  
/* LocalMotifMatrix_access_func_holder defined */ 
#ifndef DYNAMITE_DEFINED_LocalMotifMatrix_access_func_holder
typedef struct Wise2_LocalMotifMatrix_access_func_holder Wise2_LocalMotifMatrix_access_func_holder;
#define LocalMotifMatrix_access_func_holder Wise2_LocalMotifMatrix_access_func_holder
#define DYNAMITE_DEFINED_LocalMotifMatrix_access_func_holder
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  PackAln_read_Shatter_LocalMotifMatrix(mat)
 *
 * Descrip:    Reads off PackAln from shatter matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [LocalMotifMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Shatter_LocalMotifMatrix(LocalMotifMatrix * mat);
#define PackAln_read_Shatter_LocalMotifMatrix Wise2_PackAln_read_Shatter_LocalMotifMatrix


/* Function:  calculate_shatter_LocalMotifMatrix(mat,dpenv)
 *
 * Descrip:    This function calculates the LocalMotifMatrix matrix when in shatter mode
 *
 *
 * Arg:          mat [UNKN ] (null) [LocalMotifMatrix *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_shatter_LocalMotifMatrix(LocalMotifMatrix * mat,DPEnvelope * dpenv);
#define calculate_shatter_LocalMotifMatrix Wise2_calculate_shatter_LocalMotifMatrix


/* Function:  search_LocalMotifMatrix(dbsi,out,query,target,motif,ms)
 *
 * Descrip:    This function makes a database search of LocalMotifMatrix
 *             It uses the dbsi structure to choose which implementation
 *             to use of the database searching. This way at run time you
 *             can switch between single threaded/multi-threaded or hardware
 *
 *
 * Arg:          dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 * Arg:           out [UNKN ] Undocumented argument [Hscore *]
 * Arg:         query [UNKN ] Undocumented argument [ComplexSequence*]
 * Arg:        target [UNKN ] Undocumented argument [ComplexSequence*]
 * Arg:         motif [UNKN ] Undocumented argument [MotifConsMatrix*]
 * Arg:            ms [UNKN ] Undocumented argument [MotifMatrixScore*]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_search_LocalMotifMatrix(DBSearchImpl * dbsi,Hscore * out,ComplexSequence* query,ComplexSequence* target ,MotifConsMatrix* motif,MotifMatrixScore* ms);
#define search_LocalMotifMatrix Wise2_search_LocalMotifMatrix


/* Function:  serial_search_LocalMotifMatrix(out,query,target,motif,ms)
 *
 * Descrip:    This function makes a database search of LocalMotifMatrix
 *             It is a single processor implementation
 *
 *
 * Arg:           out [UNKN ] Undocumented argument [Hscore *]
 * Arg:         query [UNKN ] Undocumented argument [ComplexSequence*]
 * Arg:        target [UNKN ] Undocumented argument [ComplexSequence*]
 * Arg:         motif [UNKN ] Undocumented argument [MotifConsMatrix*]
 * Arg:            ms [UNKN ] Undocumented argument [MotifMatrixScore*]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_serial_search_LocalMotifMatrix(Hscore * out,ComplexSequence* query,ComplexSequence* target ,MotifConsMatrix* motif,MotifMatrixScore* ms);
#define serial_search_LocalMotifMatrix Wise2_serial_search_LocalMotifMatrix


/* Function:  PackAln_bestmemory_LocalMotifMatrix(query,target,motif,ms,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_LocalMotifMatrix
 *
 *
 * Arg:         query [UNKN ] query data structure [ComplexSequence*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:         motif [UNKN ] Resource [MotifConsMatrix*]
 * Arg:            ms [UNKN ] Resource [MotifMatrixScore*]
 * Arg:         dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:          dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_bestmemory_LocalMotifMatrix(ComplexSequence* query,ComplexSequence* target ,MotifConsMatrix* motif,MotifMatrixScore* ms,DPEnvelope * dpenv,DPRunImpl * dpri);
#define PackAln_bestmemory_LocalMotifMatrix Wise2_PackAln_bestmemory_LocalMotifMatrix


/* Function:  allocate_Expl_LocalMotifMatrix(query,target,motif,ms,dpri)
 *
 * Descrip:    This function allocates the LocalMotifMatrix structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_LocalMotifMatrix_only
 *
 *
 * Arg:         query [UNKN ] query data structure [ComplexSequence*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:         motif [UNKN ] Resource [MotifConsMatrix*]
 * Arg:            ms [UNKN ] Resource [MotifMatrixScore*]
 * Arg:          dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [LocalMotifMatrix *]
 *
 */
LocalMotifMatrix * Wise2_allocate_Expl_LocalMotifMatrix(ComplexSequence* query,ComplexSequence* target ,MotifConsMatrix* motif,MotifMatrixScore* ms,DPRunImpl * dpri);
#define allocate_Expl_LocalMotifMatrix Wise2_allocate_Expl_LocalMotifMatrix


/* Function:  recalculate_PackAln_LocalMotifMatrix(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by LocalMotifMatrix
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [LocalMotifMatrix *]
 *
 */
void Wise2_recalculate_PackAln_LocalMotifMatrix(PackAln * pal,LocalMotifMatrix * mat);
#define recalculate_PackAln_LocalMotifMatrix Wise2_recalculate_PackAln_LocalMotifMatrix


/* Function:  allocate_Small_LocalMotifMatrix(query,target,motif,ms)
 *
 * Descrip:    This function allocates the LocalMotifMatrix structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_LocalMotifMatrix_only
 *
 *
 * Arg:         query [UNKN ] query data structure [ComplexSequence*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:         motif [UNKN ] Resource [MotifConsMatrix*]
 * Arg:            ms [UNKN ] Resource [MotifMatrixScore*]
 *
 * Return [UNKN ]  Undocumented return value [LocalMotifMatrix *]
 *
 */
LocalMotifMatrix * Wise2_allocate_Small_LocalMotifMatrix(ComplexSequence* query,ComplexSequence* target ,MotifConsMatrix* motif,MotifMatrixScore* ms);
#define allocate_Small_LocalMotifMatrix Wise2_allocate_Small_LocalMotifMatrix


/* Function:  PackAln_calculate_Small_LocalMotifMatrix(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for LocalMotifMatrix structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_LocalMotifMatrix 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_LocalMotifMatrix 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [LocalMotifMatrix *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_calculate_Small_LocalMotifMatrix(LocalMotifMatrix * mat,DPEnvelope * dpenv);
#define PackAln_calculate_Small_LocalMotifMatrix Wise2_PackAln_calculate_Small_LocalMotifMatrix


/* Function:  AlnRangeSet_calculate_Small_LocalMotifMatrix(mat)
 *
 * Descrip:    This function calculates an alignment for LocalMotifMatrix structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_LocalMotifMatrix 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_LocalMotifMatrix
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_LocalMotifMatrix 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [LocalMotifMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_calculate_Small_LocalMotifMatrix(LocalMotifMatrix * mat);
#define AlnRangeSet_calculate_Small_LocalMotifMatrix Wise2_AlnRangeSet_calculate_Small_LocalMotifMatrix


/* Function:  AlnRangeSet_from_LocalMotifMatrix(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for LocalMotifMatrix structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_LocalMotifMatrix 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_LocalMotifMatrix
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [LocalMotifMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_from_LocalMotifMatrix(LocalMotifMatrix * mat);
#define AlnRangeSet_from_LocalMotifMatrix Wise2_AlnRangeSet_from_LocalMotifMatrix


/* Function:  convert_PackAln_to_AlnBlock_LocalMotifMatrix(pal)
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
AlnBlock * Wise2_convert_PackAln_to_AlnBlock_LocalMotifMatrix(PackAln * pal);
#define convert_PackAln_to_AlnBlock_LocalMotifMatrix Wise2_convert_PackAln_to_AlnBlock_LocalMotifMatrix


/* Function:  PackAln_read_Expl_LocalMotifMatrix(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [LocalMotifMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Expl_LocalMotifMatrix(LocalMotifMatrix * mat);
#define PackAln_read_Expl_LocalMotifMatrix Wise2_PackAln_read_Expl_LocalMotifMatrix


/* Function:  PackAln_read_generic_LocalMotifMatrix(mat,h)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [LocalMotifMatrix *]
 * Arg:          h [UNKN ] Undocumented argument [LocalMotifMatrix_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_generic_LocalMotifMatrix(LocalMotifMatrix * mat,LocalMotifMatrix_access_func_holder h);
#define PackAln_read_generic_LocalMotifMatrix Wise2_PackAln_read_generic_LocalMotifMatrix


/* Function:  calculate_LocalMotifMatrix(mat)
 *
 * Descrip:    This function calculates the LocalMotifMatrix matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_LocalMotifMatrix
 *
 *
 * Arg:        mat [UNKN ] LocalMotifMatrix which contains explicit basematrix memory [LocalMotifMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_LocalMotifMatrix(LocalMotifMatrix * mat);
#define calculate_LocalMotifMatrix Wise2_calculate_LocalMotifMatrix


/* Function:  calculate_dpenv_LocalMotifMatrix(mat,dpenv)
 *
 * Descrip:    This function calculates the LocalMotifMatrix matrix when in explicit mode, subject to the envelope
 *
 *
 * Arg:          mat [UNKN ] LocalMotifMatrix which contains explicit basematrix memory [LocalMotifMatrix *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_dpenv_LocalMotifMatrix(LocalMotifMatrix * mat,DPEnvelope * dpenv);
#define calculate_dpenv_LocalMotifMatrix Wise2_calculate_dpenv_LocalMotifMatrix


/* Function:  LocalMotifMatrix_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LocalMotifMatrix *]
 *
 */
LocalMotifMatrix * Wise2_LocalMotifMatrix_alloc(void);
#define LocalMotifMatrix_alloc Wise2_LocalMotifMatrix_alloc


/* Function:  free_LocalMotifMatrix(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [LocalMotifMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [LocalMotifMatrix *]
 *
 */
LocalMotifMatrix * Wise2_free_LocalMotifMatrix(LocalMotifMatrix * obj);
#define free_LocalMotifMatrix Wise2_free_LocalMotifMatrix


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
int Wise2_LocalMotifMatrix_shatter_access_main(LocalMotifMatrix * mat,int i,int j,int state);
#define LocalMotifMatrix_shatter_access_main Wise2_LocalMotifMatrix_shatter_access_main
int Wise2_LocalMotifMatrix_shatter_access_special(LocalMotifMatrix * mat,int i,int j,int state);
#define LocalMotifMatrix_shatter_access_special Wise2_LocalMotifMatrix_shatter_access_special
int Wise2_score_only_LocalMotifMatrix(ComplexSequence* query,ComplexSequence* target ,MotifConsMatrix* motif,MotifMatrixScore* ms);
#define score_only_LocalMotifMatrix Wise2_score_only_LocalMotifMatrix
LocalMotifMatrix * Wise2_allocate_LocalMotifMatrix_only(ComplexSequence* query,ComplexSequence* target ,MotifConsMatrix* motif,MotifMatrixScore* ms);
#define allocate_LocalMotifMatrix_only Wise2_allocate_LocalMotifMatrix_only
void Wise2_init_LocalMotifMatrix(LocalMotifMatrix * mat);
#define init_LocalMotifMatrix Wise2_init_LocalMotifMatrix
AlnRange * Wise2_AlnRange_build_LocalMotifMatrix(LocalMotifMatrix * mat,int stopj,int stopspecstate,int * startj,int * startspecstate);
#define AlnRange_build_LocalMotifMatrix Wise2_AlnRange_build_LocalMotifMatrix
boolean Wise2_read_hidden_LocalMotifMatrix(LocalMotifMatrix * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out);
#define read_hidden_LocalMotifMatrix Wise2_read_hidden_LocalMotifMatrix
int Wise2_max_hidden_LocalMotifMatrix(LocalMotifMatrix * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_hidden_LocalMotifMatrix Wise2_max_hidden_LocalMotifMatrix
boolean Wise2_read_special_strip_LocalMotifMatrix(LocalMotifMatrix * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out);
#define read_special_strip_LocalMotifMatrix Wise2_read_special_strip_LocalMotifMatrix
int Wise2_max_special_strip_LocalMotifMatrix(LocalMotifMatrix * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_special_strip_LocalMotifMatrix Wise2_max_special_strip_LocalMotifMatrix
int Wise2_max_matrix_to_special_LocalMotifMatrix(LocalMotifMatrix * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_matrix_to_special_LocalMotifMatrix Wise2_max_matrix_to_special_LocalMotifMatrix
void Wise2_calculate_hidden_LocalMotifMatrix(LocalMotifMatrix * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv);
#define calculate_hidden_LocalMotifMatrix Wise2_calculate_hidden_LocalMotifMatrix
void Wise2_init_hidden_LocalMotifMatrix(LocalMotifMatrix * mat,int starti,int startj,int stopi,int stopj);
#define init_hidden_LocalMotifMatrix Wise2_init_hidden_LocalMotifMatrix
boolean Wise2_full_dc_LocalMotifMatrix(LocalMotifMatrix * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv);
#define full_dc_LocalMotifMatrix Wise2_full_dc_LocalMotifMatrix
boolean Wise2_do_dc_single_pass_LocalMotifMatrix(LocalMotifMatrix * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done);
#define do_dc_single_pass_LocalMotifMatrix Wise2_do_dc_single_pass_LocalMotifMatrix
void Wise2_push_dc_at_merge_LocalMotifMatrix(LocalMotifMatrix * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv);
#define push_dc_at_merge_LocalMotifMatrix Wise2_push_dc_at_merge_LocalMotifMatrix
void Wise2_follow_on_dc_LocalMotifMatrix(LocalMotifMatrix * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define follow_on_dc_LocalMotifMatrix Wise2_follow_on_dc_LocalMotifMatrix
void Wise2_run_up_dc_LocalMotifMatrix(LocalMotifMatrix * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define run_up_dc_LocalMotifMatrix Wise2_run_up_dc_LocalMotifMatrix
void Wise2_init_dc_LocalMotifMatrix(LocalMotifMatrix * mat);
#define init_dc_LocalMotifMatrix Wise2_init_dc_LocalMotifMatrix
int Wise2_start_end_find_end_LocalMotifMatrix(LocalMotifMatrix * mat,int * endj);
#define start_end_find_end_LocalMotifMatrix Wise2_start_end_find_end_LocalMotifMatrix
boolean Wise2_dc_optimised_start_end_calc_LocalMotifMatrix(LocalMotifMatrix *mat,DPEnvelope * dpenv);
#define dc_optimised_start_end_calc_LocalMotifMatrix Wise2_dc_optimised_start_end_calc_LocalMotifMatrix
void Wise2_init_start_end_linear_LocalMotifMatrix(LocalMotifMatrix * mat);
#define init_start_end_linear_LocalMotifMatrix Wise2_init_start_end_linear_LocalMotifMatrix
AlnConvertSet * Wise2_AlnConvertSet_LocalMotifMatrix(void);
#define AlnConvertSet_LocalMotifMatrix Wise2_AlnConvertSet_LocalMotifMatrix
int Wise2_LocalMotifMatrix_explicit_access_main(LocalMotifMatrix * mat,int i,int j,int state);
#define LocalMotifMatrix_explicit_access_main Wise2_LocalMotifMatrix_explicit_access_main
int Wise2_LocalMotifMatrix_explicit_access_special(LocalMotifMatrix * mat,int i,int j,int state);
#define LocalMotifMatrix_explicit_access_special Wise2_LocalMotifMatrix_explicit_access_special
int Wise2_find_end_LocalMotifMatrix(LocalMotifMatrix * mat,int * ri,int * rj,int * state,boolean * isspecial,LocalMotifMatrix_access_func_holder h);
#define find_end_LocalMotifMatrix Wise2_find_end_LocalMotifMatrix
void Wise2_LocalMotifMatrix_debug_show_matrix(LocalMotifMatrix * mat,int starti,int stopi,int startj,int stopj,FILE * ofp);
#define LocalMotifMatrix_debug_show_matrix Wise2_LocalMotifMatrix_debug_show_matrix
int Wise2_max_calc_LocalMotifMatrix(LocalMotifMatrix * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,LocalMotifMatrix_access_func_holder h);
#define max_calc_LocalMotifMatrix Wise2_max_calc_LocalMotifMatrix
int Wise2_max_calc_special_LocalMotifMatrix(LocalMotifMatrix * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,LocalMotifMatrix_access_func_holder h);
#define max_calc_special_LocalMotifMatrix Wise2_max_calc_special_LocalMotifMatrix

#ifdef _cplusplus
}
#endif

#endif
