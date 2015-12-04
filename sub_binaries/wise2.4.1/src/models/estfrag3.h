#ifndef DYNAMITEestfrag3HEADERFILE
#define DYNAMITEestfrag3HEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"
#include "cdparser.h"
#include "genewisemodel.h"
#include "genewisemodeldb.h"




struct Wise2_EstFrag3 {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    BaseMatrix * basematrix;     
    ShatterMatrix * shatter;     
    int leni;    
    int lenj;    
    GeneWiseScore* query;    
    ComplexSequence* target;     
    cDNAParserScore * cp;    
    Score start_frag;    
    Score end_frag;  
    } ;  
/* EstFrag3 defined */ 
#ifndef DYNAMITE_DEFINED_EstFrag3
typedef struct Wise2_EstFrag3 Wise2_EstFrag3;
#define EstFrag3 Wise2_EstFrag3
#define DYNAMITE_DEFINED_EstFrag3
#endif


struct Wise2_EstFrag3_access_func_holder {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int (*access_main)(EstFrag3*,int,int,int);   
    int (*access_special)(EstFrag3*,int,int,int);    
    } ;  
/* EstFrag3_access_func_holder defined */ 
#ifndef DYNAMITE_DEFINED_EstFrag3_access_func_holder
typedef struct Wise2_EstFrag3_access_func_holder Wise2_EstFrag3_access_func_holder;
#define EstFrag3_access_func_holder Wise2_EstFrag3_access_func_holder
#define DYNAMITE_DEFINED_EstFrag3_access_func_holder
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  PackAln_read_Shatter_EstFrag3(mat)
 *
 * Descrip:    Reads off PackAln from shatter matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [EstFrag3 *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Shatter_EstFrag3(EstFrag3 * mat);
#define PackAln_read_Shatter_EstFrag3 Wise2_PackAln_read_Shatter_EstFrag3


/* Function:  calculate_shatter_EstFrag3(mat,dpenv)
 *
 * Descrip:    This function calculates the EstFrag3 matrix when in shatter mode
 *
 *
 * Arg:          mat [UNKN ] (null) [EstFrag3 *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_shatter_EstFrag3(EstFrag3 * mat,DPEnvelope * dpenv);
#define calculate_shatter_EstFrag3 Wise2_calculate_shatter_EstFrag3


/* Function:  search_EstFrag3(dbsi,out,querydb,targetdb,cp,start_frag,end_frag)
 *
 * Descrip:    This function makes a database search of EstFrag3
 *             It uses the dbsi structure to choose which implementation
 *             to use of the database searching. This way at run time you
 *             can switch between single threaded/multi-threaded or hardware
 *
 *
 * Arg:              dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 * Arg:               out [UNKN ] Undocumented argument [Hscore *]
 * Arg:           querydb [UNKN ] Undocumented argument [GeneWiseDB*]
 * Arg:          targetdb [UNKN ] Undocumented argument [cDNADB*]
 * Arg:                cp [UNKN ] Undocumented argument [cDNAParserScore *]
 * Arg:        start_frag [UNKN ] Undocumented argument [Score]
 * Arg:          end_frag [UNKN ] Undocumented argument [Score]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_search_EstFrag3(DBSearchImpl * dbsi,Hscore * out,GeneWiseDB* querydb,cDNADB* targetdb ,cDNAParserScore * cp,Score start_frag,Score end_frag);
#define search_EstFrag3 Wise2_search_EstFrag3


/* Function:  serial_search_EstFrag3(out,querydb,targetdb,cp,start_frag,end_frag)
 *
 * Descrip:    This function makes a database search of EstFrag3
 *             It is a single processor implementation
 *
 *
 * Arg:               out [UNKN ] Undocumented argument [Hscore *]
 * Arg:           querydb [UNKN ] Undocumented argument [GeneWiseDB*]
 * Arg:          targetdb [UNKN ] Undocumented argument [cDNADB*]
 * Arg:                cp [UNKN ] Undocumented argument [cDNAParserScore *]
 * Arg:        start_frag [UNKN ] Undocumented argument [Score]
 * Arg:          end_frag [UNKN ] Undocumented argument [Score]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_serial_search_EstFrag3(Hscore * out,GeneWiseDB* querydb,cDNADB* targetdb ,cDNAParserScore * cp,Score start_frag,Score end_frag);
#define serial_search_EstFrag3 Wise2_serial_search_EstFrag3


/* Function:  PackAln_bestmemory_EstFrag3(query,target,cp,start_frag,end_frag,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_EstFrag3
 *
 *
 * Arg:             query [UNKN ] query data structure [GeneWiseScore*]
 * Arg:            target [UNKN ] target data structure [ComplexSequence*]
 * Arg:                cp [UNKN ] Resource [cDNAParserScore *]
 * Arg:        start_frag [UNKN ] Resource [Score]
 * Arg:          end_frag [UNKN ] Resource [Score]
 * Arg:             dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:              dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_bestmemory_EstFrag3(GeneWiseScore* query,ComplexSequence* target ,cDNAParserScore * cp,Score start_frag,Score end_frag,DPEnvelope * dpenv,DPRunImpl * dpri);
#define PackAln_bestmemory_EstFrag3 Wise2_PackAln_bestmemory_EstFrag3


/* Function:  allocate_Expl_EstFrag3(query,target,cp,start_frag,end_frag,dpri)
 *
 * Descrip:    This function allocates the EstFrag3 structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_EstFrag3_only
 *
 *
 * Arg:             query [UNKN ] query data structure [GeneWiseScore*]
 * Arg:            target [UNKN ] target data structure [ComplexSequence*]
 * Arg:                cp [UNKN ] Resource [cDNAParserScore *]
 * Arg:        start_frag [UNKN ] Resource [Score]
 * Arg:          end_frag [UNKN ] Resource [Score]
 * Arg:              dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [EstFrag3 *]
 *
 */
EstFrag3 * Wise2_allocate_Expl_EstFrag3(GeneWiseScore* query,ComplexSequence* target ,cDNAParserScore * cp,Score start_frag,Score end_frag,DPRunImpl * dpri);
#define allocate_Expl_EstFrag3 Wise2_allocate_Expl_EstFrag3


/* Function:  recalculate_PackAln_EstFrag3(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by EstFrag3
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [EstFrag3 *]
 *
 */
void Wise2_recalculate_PackAln_EstFrag3(PackAln * pal,EstFrag3 * mat);
#define recalculate_PackAln_EstFrag3 Wise2_recalculate_PackAln_EstFrag3


/* Function:  allocate_Small_EstFrag3(query,target,cp,start_frag,end_frag)
 *
 * Descrip:    This function allocates the EstFrag3 structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_EstFrag3_only
 *
 *
 * Arg:             query [UNKN ] query data structure [GeneWiseScore*]
 * Arg:            target [UNKN ] target data structure [ComplexSequence*]
 * Arg:                cp [UNKN ] Resource [cDNAParserScore *]
 * Arg:        start_frag [UNKN ] Resource [Score]
 * Arg:          end_frag [UNKN ] Resource [Score]
 *
 * Return [UNKN ]  Undocumented return value [EstFrag3 *]
 *
 */
EstFrag3 * Wise2_allocate_Small_EstFrag3(GeneWiseScore* query,ComplexSequence* target ,cDNAParserScore * cp,Score start_frag,Score end_frag);
#define allocate_Small_EstFrag3 Wise2_allocate_Small_EstFrag3


/* Function:  PackAln_calculate_Small_EstFrag3(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for EstFrag3 structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_EstFrag3 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_EstFrag3 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [EstFrag3 *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_calculate_Small_EstFrag3(EstFrag3 * mat,DPEnvelope * dpenv);
#define PackAln_calculate_Small_EstFrag3 Wise2_PackAln_calculate_Small_EstFrag3


/* Function:  AlnRangeSet_calculate_Small_EstFrag3(mat)
 *
 * Descrip:    This function calculates an alignment for EstFrag3 structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_EstFrag3 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_EstFrag3
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_EstFrag3 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [EstFrag3 *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_calculate_Small_EstFrag3(EstFrag3 * mat);
#define AlnRangeSet_calculate_Small_EstFrag3 Wise2_AlnRangeSet_calculate_Small_EstFrag3


/* Function:  AlnRangeSet_from_EstFrag3(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for EstFrag3 structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_EstFrag3 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_EstFrag3
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [EstFrag3 *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_from_EstFrag3(EstFrag3 * mat);
#define AlnRangeSet_from_EstFrag3 Wise2_AlnRangeSet_from_EstFrag3


/* Function:  convert_PackAln_to_AlnBlock_EstFrag3(pal)
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
AlnBlock * Wise2_convert_PackAln_to_AlnBlock_EstFrag3(PackAln * pal);
#define convert_PackAln_to_AlnBlock_EstFrag3 Wise2_convert_PackAln_to_AlnBlock_EstFrag3


/* Function:  PackAln_read_Expl_EstFrag3(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [EstFrag3 *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Expl_EstFrag3(EstFrag3 * mat);
#define PackAln_read_Expl_EstFrag3 Wise2_PackAln_read_Expl_EstFrag3


/* Function:  PackAln_read_generic_EstFrag3(mat,h)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [EstFrag3 *]
 * Arg:          h [UNKN ] Undocumented argument [EstFrag3_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_generic_EstFrag3(EstFrag3 * mat,EstFrag3_access_func_holder h);
#define PackAln_read_generic_EstFrag3 Wise2_PackAln_read_generic_EstFrag3


/* Function:  calculate_EstFrag3(mat)
 *
 * Descrip:    This function calculates the EstFrag3 matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_EstFrag3
 *
 *
 * Arg:        mat [UNKN ] EstFrag3 which contains explicit basematrix memory [EstFrag3 *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_EstFrag3(EstFrag3 * mat);
#define calculate_EstFrag3 Wise2_calculate_EstFrag3


/* Function:  calculate_dpenv_EstFrag3(mat,dpenv)
 *
 * Descrip:    This function calculates the EstFrag3 matrix when in explicit mode, subject to the envelope
 *
 *
 * Arg:          mat [UNKN ] EstFrag3 which contains explicit basematrix memory [EstFrag3 *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_dpenv_EstFrag3(EstFrag3 * mat,DPEnvelope * dpenv);
#define calculate_dpenv_EstFrag3 Wise2_calculate_dpenv_EstFrag3


/* Function:  EstFrag3_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [EstFrag3 *]
 *
 */
EstFrag3 * Wise2_EstFrag3_alloc(void);
#define EstFrag3_alloc Wise2_EstFrag3_alloc


/* Function:  free_EstFrag3(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [EstFrag3 *]
 *
 * Return [UNKN ]  Undocumented return value [EstFrag3 *]
 *
 */
EstFrag3 * Wise2_free_EstFrag3(EstFrag3 * obj);
#define free_EstFrag3 Wise2_free_EstFrag3


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
int Wise2_EstFrag3_shatter_access_main(EstFrag3 * mat,int i,int j,int state);
#define EstFrag3_shatter_access_main Wise2_EstFrag3_shatter_access_main
int Wise2_EstFrag3_shatter_access_special(EstFrag3 * mat,int i,int j,int state);
#define EstFrag3_shatter_access_special Wise2_EstFrag3_shatter_access_special
int Wise2_score_only_EstFrag3(GeneWiseScore* query,ComplexSequence* target ,cDNAParserScore * cp,Score start_frag,Score end_frag);
#define score_only_EstFrag3 Wise2_score_only_EstFrag3
EstFrag3 * Wise2_allocate_EstFrag3_only(GeneWiseScore* query,ComplexSequence* target ,cDNAParserScore * cp,Score start_frag,Score end_frag);
#define allocate_EstFrag3_only Wise2_allocate_EstFrag3_only
void Wise2_init_EstFrag3(EstFrag3 * mat);
#define init_EstFrag3 Wise2_init_EstFrag3
AlnRange * Wise2_AlnRange_build_EstFrag3(EstFrag3 * mat,int stopj,int stopspecstate,int * startj,int * startspecstate);
#define AlnRange_build_EstFrag3 Wise2_AlnRange_build_EstFrag3
boolean Wise2_read_hidden_EstFrag3(EstFrag3 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out);
#define read_hidden_EstFrag3 Wise2_read_hidden_EstFrag3
int Wise2_max_hidden_EstFrag3(EstFrag3 * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_hidden_EstFrag3 Wise2_max_hidden_EstFrag3
boolean Wise2_read_special_strip_EstFrag3(EstFrag3 * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out);
#define read_special_strip_EstFrag3 Wise2_read_special_strip_EstFrag3
int Wise2_max_special_strip_EstFrag3(EstFrag3 * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_special_strip_EstFrag3 Wise2_max_special_strip_EstFrag3
int Wise2_max_matrix_to_special_EstFrag3(EstFrag3 * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_matrix_to_special_EstFrag3 Wise2_max_matrix_to_special_EstFrag3
void Wise2_calculate_hidden_EstFrag3(EstFrag3 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv);
#define calculate_hidden_EstFrag3 Wise2_calculate_hidden_EstFrag3
void Wise2_init_hidden_EstFrag3(EstFrag3 * mat,int starti,int startj,int stopi,int stopj);
#define init_hidden_EstFrag3 Wise2_init_hidden_EstFrag3
boolean Wise2_full_dc_EstFrag3(EstFrag3 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv);
#define full_dc_EstFrag3 Wise2_full_dc_EstFrag3
boolean Wise2_do_dc_single_pass_EstFrag3(EstFrag3 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done);
#define do_dc_single_pass_EstFrag3 Wise2_do_dc_single_pass_EstFrag3
void Wise2_push_dc_at_merge_EstFrag3(EstFrag3 * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv);
#define push_dc_at_merge_EstFrag3 Wise2_push_dc_at_merge_EstFrag3
void Wise2_follow_on_dc_EstFrag3(EstFrag3 * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define follow_on_dc_EstFrag3 Wise2_follow_on_dc_EstFrag3
void Wise2_run_up_dc_EstFrag3(EstFrag3 * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define run_up_dc_EstFrag3 Wise2_run_up_dc_EstFrag3
void Wise2_init_dc_EstFrag3(EstFrag3 * mat);
#define init_dc_EstFrag3 Wise2_init_dc_EstFrag3
int Wise2_start_end_find_end_EstFrag3(EstFrag3 * mat,int * endj);
#define start_end_find_end_EstFrag3 Wise2_start_end_find_end_EstFrag3
boolean Wise2_dc_optimised_start_end_calc_EstFrag3(EstFrag3 *mat,DPEnvelope * dpenv);
#define dc_optimised_start_end_calc_EstFrag3 Wise2_dc_optimised_start_end_calc_EstFrag3
void Wise2_init_start_end_linear_EstFrag3(EstFrag3 * mat);
#define init_start_end_linear_EstFrag3 Wise2_init_start_end_linear_EstFrag3
AlnConvertSet * Wise2_AlnConvertSet_EstFrag3(void);
#define AlnConvertSet_EstFrag3 Wise2_AlnConvertSet_EstFrag3
int Wise2_EstFrag3_explicit_access_main(EstFrag3 * mat,int i,int j,int state);
#define EstFrag3_explicit_access_main Wise2_EstFrag3_explicit_access_main
int Wise2_EstFrag3_explicit_access_special(EstFrag3 * mat,int i,int j,int state);
#define EstFrag3_explicit_access_special Wise2_EstFrag3_explicit_access_special
int Wise2_find_end_EstFrag3(EstFrag3 * mat,int * ri,int * rj,int * state,boolean * isspecial,EstFrag3_access_func_holder h);
#define find_end_EstFrag3 Wise2_find_end_EstFrag3
void Wise2_EstFrag3_debug_show_matrix(EstFrag3 * mat,int starti,int stopi,int startj,int stopj,FILE * ofp);
#define EstFrag3_debug_show_matrix Wise2_EstFrag3_debug_show_matrix
int Wise2_max_calc_EstFrag3(EstFrag3 * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,EstFrag3_access_func_holder h);
#define max_calc_EstFrag3 Wise2_max_calc_EstFrag3
int Wise2_max_calc_special_EstFrag3(EstFrag3 * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,EstFrag3_access_func_holder h);
#define max_calc_special_EstFrag3 Wise2_max_calc_special_EstFrag3

#ifdef _cplusplus
}
#endif

#endif
