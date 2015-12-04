#ifndef DYNAMITEgenestretch6HEADERFILE
#define DYNAMITEgenestretch6HEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"
#include "geneparser4.h"
#include "genewisemodel.h"
#include "genewisemodeldb.h"
#include "genestats.h"

struct Wise2_GeneStretch6 {  
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
    GeneParser4Score * gp;   
    GeneralGeneModelScore * general_model;   
    } ;  
/* GeneStretch6 defined */ 
#ifndef DYNAMITE_DEFINED_GeneStretch6
typedef struct Wise2_GeneStretch6 Wise2_GeneStretch6;
#define GeneStretch6 Wise2_GeneStretch6
#define DYNAMITE_DEFINED_GeneStretch6
#endif


#ifdef PTHREAD
struct thread_pool_holder_GeneStretch6 {  
    GeneWiseScore* query;   /* Query object placeholder */ 
    GeneWiseDB* querydb;/* Query database object */ 
    boolean query_init;  
    ComplexSequence* target;/* Target object placeholder */ 
    GenomicDB* targetdb;/* Target database object */ 
    boolean target_init; 
    GeneParser4Score * gp;   
    GeneralGeneModelScore * general_model;   
    pthread_mutex_t input_lock;  
    pthread_mutex_t output_lock;     
    Hscore * out;    
    pthread_t * pool;    
    int number_of_threads;   
    boolean search_has_ended;    
    DBSearchImpl * dbsi;     
    } ;  
#endif /* PTHREAD */
struct Wise2_GeneStretch6_access_func_holder {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int (*access_main)(GeneStretch6*,int,int,int);   
    int (*access_special)(GeneStretch6*,int,int,int);    
    } ;  
/* GeneStretch6_access_func_holder defined */ 
#ifndef DYNAMITE_DEFINED_GeneStretch6_access_func_holder
typedef struct Wise2_GeneStretch6_access_func_holder Wise2_GeneStretch6_access_func_holder;
#define GeneStretch6_access_func_holder Wise2_GeneStretch6_access_func_holder
#define DYNAMITE_DEFINED_GeneStretch6_access_func_holder
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  PackAln_read_Shatter_GeneStretch6(mat)
 *
 * Descrip:    Reads off PackAln from shatter matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [GeneStretch6 *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Shatter_GeneStretch6(GeneStretch6 * mat);
#define PackAln_read_Shatter_GeneStretch6 Wise2_PackAln_read_Shatter_GeneStretch6


/* Function:  calculate_shatter_GeneStretch6(mat,dpenv)
 *
 * Descrip:    This function calculates the GeneStretch6 matrix when in shatter mode
 *
 *
 * Arg:          mat [UNKN ] (null) [GeneStretch6 *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_shatter_GeneStretch6(GeneStretch6 * mat,DPEnvelope * dpenv);
#define calculate_shatter_GeneStretch6 Wise2_calculate_shatter_GeneStretch6


/* Function:  search_GeneStretch6(dbsi,out,querydb,targetdb,gp,general_model)
 *
 * Descrip:    This function makes a database search of GeneStretch6
 *             It uses the dbsi structure to choose which implementation
 *             to use of the database searching. This way at run time you
 *             can switch between single threaded/multi-threaded or hardware
 *
 *
 * Arg:                 dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 * Arg:                  out [UNKN ] Undocumented argument [Hscore *]
 * Arg:              querydb [UNKN ] Undocumented argument [GeneWiseDB*]
 * Arg:             targetdb [UNKN ] Undocumented argument [GenomicDB*]
 * Arg:                   gp [UNKN ] Undocumented argument [GeneParser4Score *]
 * Arg:        general_model [UNKN ] Undocumented argument [GeneralGeneModelScore *]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_search_GeneStretch6(DBSearchImpl * dbsi,Hscore * out,GeneWiseDB* querydb,GenomicDB* targetdb ,GeneParser4Score * gp,GeneralGeneModelScore * general_model);
#define search_GeneStretch6 Wise2_search_GeneStretch6


/* Function:  serial_search_GeneStretch6(out,querydb,targetdb,gp,general_model)
 *
 * Descrip:    This function makes a database search of GeneStretch6
 *             It is a single processor implementation
 *
 *
 * Arg:                  out [UNKN ] Undocumented argument [Hscore *]
 * Arg:              querydb [UNKN ] Undocumented argument [GeneWiseDB*]
 * Arg:             targetdb [UNKN ] Undocumented argument [GenomicDB*]
 * Arg:                   gp [UNKN ] Undocumented argument [GeneParser4Score *]
 * Arg:        general_model [UNKN ] Undocumented argument [GeneralGeneModelScore *]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_serial_search_GeneStretch6(Hscore * out,GeneWiseDB* querydb,GenomicDB* targetdb ,GeneParser4Score * gp,GeneralGeneModelScore * general_model);
#define serial_search_GeneStretch6 Wise2_serial_search_GeneStretch6


/* Function:  PackAln_bestmemory_GeneStretch6(query,target,gp,general_model,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_GeneStretch6
 *
 *
 * Arg:                query [UNKN ] query data structure [GeneWiseScore*]
 * Arg:               target [UNKN ] target data structure [ComplexSequence*]
 * Arg:                   gp [UNKN ] Resource [GeneParser4Score *]
 * Arg:        general_model [UNKN ] Resource [GeneralGeneModelScore *]
 * Arg:                dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:                 dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_bestmemory_GeneStretch6(GeneWiseScore* query,ComplexSequence* target ,GeneParser4Score * gp,GeneralGeneModelScore * general_model,DPEnvelope * dpenv,DPRunImpl * dpri);
#define PackAln_bestmemory_GeneStretch6 Wise2_PackAln_bestmemory_GeneStretch6


/* Function:  allocate_Expl_GeneStretch6(query,target,gp,general_model,dpri)
 *
 * Descrip:    This function allocates the GeneStretch6 structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_GeneStretch6_only
 *
 *
 * Arg:                query [UNKN ] query data structure [GeneWiseScore*]
 * Arg:               target [UNKN ] target data structure [ComplexSequence*]
 * Arg:                   gp [UNKN ] Resource [GeneParser4Score *]
 * Arg:        general_model [UNKN ] Resource [GeneralGeneModelScore *]
 * Arg:                 dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [GeneStretch6 *]
 *
 */
GeneStretch6 * Wise2_allocate_Expl_GeneStretch6(GeneWiseScore* query,ComplexSequence* target ,GeneParser4Score * gp,GeneralGeneModelScore * general_model,DPRunImpl * dpri);
#define allocate_Expl_GeneStretch6 Wise2_allocate_Expl_GeneStretch6


/* Function:  recalculate_PackAln_GeneStretch6(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by GeneStretch6
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [GeneStretch6 *]
 *
 */
void Wise2_recalculate_PackAln_GeneStretch6(PackAln * pal,GeneStretch6 * mat);
#define recalculate_PackAln_GeneStretch6 Wise2_recalculate_PackAln_GeneStretch6


/* Function:  allocate_Small_GeneStretch6(query,target,gp,general_model)
 *
 * Descrip:    This function allocates the GeneStretch6 structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_GeneStretch6_only
 *
 *
 * Arg:                query [UNKN ] query data structure [GeneWiseScore*]
 * Arg:               target [UNKN ] target data structure [ComplexSequence*]
 * Arg:                   gp [UNKN ] Resource [GeneParser4Score *]
 * Arg:        general_model [UNKN ] Resource [GeneralGeneModelScore *]
 *
 * Return [UNKN ]  Undocumented return value [GeneStretch6 *]
 *
 */
GeneStretch6 * Wise2_allocate_Small_GeneStretch6(GeneWiseScore* query,ComplexSequence* target ,GeneParser4Score * gp,GeneralGeneModelScore * general_model);
#define allocate_Small_GeneStretch6 Wise2_allocate_Small_GeneStretch6


/* Function:  PackAln_calculate_Small_GeneStretch6(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for GeneStretch6 structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_GeneStretch6 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_GeneStretch6 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [GeneStretch6 *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_calculate_Small_GeneStretch6(GeneStretch6 * mat,DPEnvelope * dpenv);
#define PackAln_calculate_Small_GeneStretch6 Wise2_PackAln_calculate_Small_GeneStretch6


/* Function:  AlnRangeSet_calculate_Small_GeneStretch6(mat)
 *
 * Descrip:    This function calculates an alignment for GeneStretch6 structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_GeneStretch6 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_GeneStretch6
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_GeneStretch6 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [GeneStretch6 *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_calculate_Small_GeneStretch6(GeneStretch6 * mat);
#define AlnRangeSet_calculate_Small_GeneStretch6 Wise2_AlnRangeSet_calculate_Small_GeneStretch6


/* Function:  AlnRangeSet_from_GeneStretch6(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for GeneStretch6 structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_GeneStretch6 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_GeneStretch6
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [GeneStretch6 *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_from_GeneStretch6(GeneStretch6 * mat);
#define AlnRangeSet_from_GeneStretch6 Wise2_AlnRangeSet_from_GeneStretch6


/* Function:  convert_PackAln_to_AlnBlock_GeneStretch6(pal)
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
AlnBlock * Wise2_convert_PackAln_to_AlnBlock_GeneStretch6(PackAln * pal);
#define convert_PackAln_to_AlnBlock_GeneStretch6 Wise2_convert_PackAln_to_AlnBlock_GeneStretch6


/* Function:  PackAln_read_Expl_GeneStretch6(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [GeneStretch6 *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Expl_GeneStretch6(GeneStretch6 * mat);
#define PackAln_read_Expl_GeneStretch6 Wise2_PackAln_read_Expl_GeneStretch6


/* Function:  PackAln_read_generic_GeneStretch6(mat,h)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [GeneStretch6 *]
 * Arg:          h [UNKN ] Undocumented argument [GeneStretch6_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_generic_GeneStretch6(GeneStretch6 * mat,GeneStretch6_access_func_holder h);
#define PackAln_read_generic_GeneStretch6 Wise2_PackAln_read_generic_GeneStretch6


/* Function:  calculate_GeneStretch6(mat)
 *
 * Descrip:    This function calculates the GeneStretch6 matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_GeneStretch6
 *
 *
 * Arg:        mat [UNKN ] GeneStretch6 which contains explicit basematrix memory [GeneStretch6 *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_GeneStretch6(GeneStretch6 * mat);
#define calculate_GeneStretch6 Wise2_calculate_GeneStretch6


/* Function:  calculate_dpenv_GeneStretch6(mat,dpenv)
 *
 * Descrip:    This function calculates the GeneStretch6 matrix when in explicit mode, subject to the envelope
 *
 *
 * Arg:          mat [UNKN ] GeneStretch6 which contains explicit basematrix memory [GeneStretch6 *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_dpenv_GeneStretch6(GeneStretch6 * mat,DPEnvelope * dpenv);
#define calculate_dpenv_GeneStretch6 Wise2_calculate_dpenv_GeneStretch6


/* Function:  GeneStretch6_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneStretch6 *]
 *
 */
GeneStretch6 * Wise2_GeneStretch6_alloc(void);
#define GeneStretch6_alloc Wise2_GeneStretch6_alloc


/* Function:  free_GeneStretch6(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneStretch6 *]
 *
 * Return [UNKN ]  Undocumented return value [GeneStretch6 *]
 *
 */
GeneStretch6 * Wise2_free_GeneStretch6(GeneStretch6 * obj);
#define free_GeneStretch6 Wise2_free_GeneStretch6


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
int Wise2_GeneStretch6_shatter_access_main(GeneStretch6 * mat,int i,int j,int state);
#define GeneStretch6_shatter_access_main Wise2_GeneStretch6_shatter_access_main
int Wise2_GeneStretch6_shatter_access_special(GeneStretch6 * mat,int i,int j,int state);
#define GeneStretch6_shatter_access_special Wise2_GeneStretch6_shatter_access_special
void * Wise2_thread_loop_GeneStretch6(void * ptr);
#define thread_loop_GeneStretch6 Wise2_thread_loop_GeneStretch6
int Wise2_score_only_GeneStretch6(GeneWiseScore* query,ComplexSequence* target ,GeneParser4Score * gp,GeneralGeneModelScore * general_model);
#define score_only_GeneStretch6 Wise2_score_only_GeneStretch6
GeneStretch6 * Wise2_allocate_GeneStretch6_only(GeneWiseScore* query,ComplexSequence* target ,GeneParser4Score * gp,GeneralGeneModelScore * general_model);
#define allocate_GeneStretch6_only Wise2_allocate_GeneStretch6_only
void Wise2_init_GeneStretch6(GeneStretch6 * mat);
#define init_GeneStretch6 Wise2_init_GeneStretch6
AlnRange * Wise2_AlnRange_build_GeneStretch6(GeneStretch6 * mat,int stopj,int stopspecstate,int * startj,int * startspecstate);
#define AlnRange_build_GeneStretch6 Wise2_AlnRange_build_GeneStretch6
boolean Wise2_read_hidden_GeneStretch6(GeneStretch6 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out);
#define read_hidden_GeneStretch6 Wise2_read_hidden_GeneStretch6
int Wise2_max_hidden_GeneStretch6(GeneStretch6 * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_hidden_GeneStretch6 Wise2_max_hidden_GeneStretch6
boolean Wise2_read_special_strip_GeneStretch6(GeneStretch6 * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out);
#define read_special_strip_GeneStretch6 Wise2_read_special_strip_GeneStretch6
int Wise2_max_special_strip_GeneStretch6(GeneStretch6 * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_special_strip_GeneStretch6 Wise2_max_special_strip_GeneStretch6
int Wise2_max_matrix_to_special_GeneStretch6(GeneStretch6 * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_matrix_to_special_GeneStretch6 Wise2_max_matrix_to_special_GeneStretch6
void Wise2_calculate_hidden_GeneStretch6(GeneStretch6 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv);
#define calculate_hidden_GeneStretch6 Wise2_calculate_hidden_GeneStretch6
void Wise2_init_hidden_GeneStretch6(GeneStretch6 * mat,int starti,int startj,int stopi,int stopj);
#define init_hidden_GeneStretch6 Wise2_init_hidden_GeneStretch6
boolean Wise2_full_dc_GeneStretch6(GeneStretch6 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv);
#define full_dc_GeneStretch6 Wise2_full_dc_GeneStretch6
boolean Wise2_do_dc_single_pass_GeneStretch6(GeneStretch6 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done);
#define do_dc_single_pass_GeneStretch6 Wise2_do_dc_single_pass_GeneStretch6
void Wise2_push_dc_at_merge_GeneStretch6(GeneStretch6 * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv);
#define push_dc_at_merge_GeneStretch6 Wise2_push_dc_at_merge_GeneStretch6
void Wise2_follow_on_dc_GeneStretch6(GeneStretch6 * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define follow_on_dc_GeneStretch6 Wise2_follow_on_dc_GeneStretch6
void Wise2_run_up_dc_GeneStretch6(GeneStretch6 * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define run_up_dc_GeneStretch6 Wise2_run_up_dc_GeneStretch6
void Wise2_init_dc_GeneStretch6(GeneStretch6 * mat);
#define init_dc_GeneStretch6 Wise2_init_dc_GeneStretch6
int Wise2_start_end_find_end_GeneStretch6(GeneStretch6 * mat,int * endj);
#define start_end_find_end_GeneStretch6 Wise2_start_end_find_end_GeneStretch6
boolean Wise2_dc_optimised_start_end_calc_GeneStretch6(GeneStretch6 *mat,DPEnvelope * dpenv);
#define dc_optimised_start_end_calc_GeneStretch6 Wise2_dc_optimised_start_end_calc_GeneStretch6
void Wise2_init_start_end_linear_GeneStretch6(GeneStretch6 * mat);
#define init_start_end_linear_GeneStretch6 Wise2_init_start_end_linear_GeneStretch6
AlnConvertSet * Wise2_AlnConvertSet_GeneStretch6(void);
#define AlnConvertSet_GeneStretch6 Wise2_AlnConvertSet_GeneStretch6
int Wise2_GeneStretch6_explicit_access_main(GeneStretch6 * mat,int i,int j,int state);
#define GeneStretch6_explicit_access_main Wise2_GeneStretch6_explicit_access_main
int Wise2_GeneStretch6_explicit_access_special(GeneStretch6 * mat,int i,int j,int state);
#define GeneStretch6_explicit_access_special Wise2_GeneStretch6_explicit_access_special
int Wise2_find_end_GeneStretch6(GeneStretch6 * mat,int * ri,int * rj,int * state,boolean * isspecial,GeneStretch6_access_func_holder h);
#define find_end_GeneStretch6 Wise2_find_end_GeneStretch6
void Wise2_GeneStretch6_debug_show_matrix(GeneStretch6 * mat,int starti,int stopi,int startj,int stopj,FILE * ofp);
#define GeneStretch6_debug_show_matrix Wise2_GeneStretch6_debug_show_matrix
int Wise2_max_calc_GeneStretch6(GeneStretch6 * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,GeneStretch6_access_func_holder h);
#define max_calc_GeneStretch6 Wise2_max_calc_GeneStretch6
int Wise2_max_calc_special_GeneStretch6(GeneStretch6 * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,GeneStretch6_access_func_holder h);
#define max_calc_special_GeneStretch6 Wise2_max_calc_special_GeneStretch6

#ifdef _cplusplus
}
#endif

#endif
