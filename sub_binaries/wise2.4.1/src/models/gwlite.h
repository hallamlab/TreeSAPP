#ifndef DYNAMITEgwliteHEADERFILE
#define DYNAMITEgwliteHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "gwlitemodel.h"
#include "geneparser4.h"
#include "genewisemodeldb.h"

struct Wise2_GeneLiteModel {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    BaseMatrix * basematrix;     
    ShatterMatrix * shatter;     
    int leni;    
    int lenj;    
    GwLiteScore* query;  
    ComplexSequence* target;     
    GeneParser4Score * gp;   
    } ;  
/* GeneLiteModel defined */ 
#ifndef DYNAMITE_DEFINED_GeneLiteModel
typedef struct Wise2_GeneLiteModel Wise2_GeneLiteModel;
#define GeneLiteModel Wise2_GeneLiteModel
#define DYNAMITE_DEFINED_GeneLiteModel
#endif


#ifdef PTHREAD
struct thread_pool_holder_GeneLiteModel {  
    GwLiteScore* query; /* Query object placeholder */ 
    GeneWiseDB* querydb;/* Query database object */ 
    boolean query_init;  
    ComplexSequence* target;/* Target object placeholder */ 
    GenomicDB* targetdb;/* Target database object */ 
    boolean target_init; 
    GeneParser4Score * gp;   
    pthread_mutex_t input_lock;  
    pthread_mutex_t output_lock;     
    Hscore * out;    
    pthread_t * pool;    
    int number_of_threads;   
    boolean search_has_ended;    
    DBSearchImpl * dbsi;     
    } ;  
#endif /* PTHREAD */
struct Wise2_GeneLiteModel_access_func_holder {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int (*access_main)(GeneLiteModel*,int,int,int);  
    int (*access_special)(GeneLiteModel*,int,int,int);   
    } ;  
/* GeneLiteModel_access_func_holder defined */ 
#ifndef DYNAMITE_DEFINED_GeneLiteModel_access_func_holder
typedef struct Wise2_GeneLiteModel_access_func_holder Wise2_GeneLiteModel_access_func_holder;
#define GeneLiteModel_access_func_holder Wise2_GeneLiteModel_access_func_holder
#define DYNAMITE_DEFINED_GeneLiteModel_access_func_holder
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  PackAln_read_Shatter_GeneLiteModel(mat)
 *
 * Descrip:    Reads off PackAln from shatter matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [GeneLiteModel *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Shatter_GeneLiteModel(GeneLiteModel * mat);
#define PackAln_read_Shatter_GeneLiteModel Wise2_PackAln_read_Shatter_GeneLiteModel


/* Function:  calculate_shatter_GeneLiteModel(mat,dpenv)
 *
 * Descrip:    This function calculates the GeneLiteModel matrix when in shatter mode
 *
 *
 * Arg:          mat [UNKN ] (null) [GeneLiteModel *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_shatter_GeneLiteModel(GeneLiteModel * mat,DPEnvelope * dpenv);
#define calculate_shatter_GeneLiteModel Wise2_calculate_shatter_GeneLiteModel


/* Function:  search_GeneLiteModel(dbsi,out,querydb,targetdb,gp)
 *
 * Descrip:    This function makes a database search of GeneLiteModel
 *             It uses the dbsi structure to choose which implementation
 *             to use of the database searching. This way at run time you
 *             can switch between single threaded/multi-threaded or hardware
 *
 *
 * Arg:            dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 * Arg:             out [UNKN ] Undocumented argument [Hscore *]
 * Arg:         querydb [UNKN ] Undocumented argument [GeneWiseDB*]
 * Arg:        targetdb [UNKN ] Undocumented argument [GenomicDB*]
 * Arg:              gp [UNKN ] Undocumented argument [GeneParser4Score *]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_search_GeneLiteModel(DBSearchImpl * dbsi,Hscore * out,GeneWiseDB* querydb,GenomicDB* targetdb ,GeneParser4Score * gp);
#define search_GeneLiteModel Wise2_search_GeneLiteModel


/* Function:  serial_search_GeneLiteModel(out,querydb,targetdb,gp)
 *
 * Descrip:    This function makes a database search of GeneLiteModel
 *             It is a single processor implementation
 *
 *
 * Arg:             out [UNKN ] Undocumented argument [Hscore *]
 * Arg:         querydb [UNKN ] Undocumented argument [GeneWiseDB*]
 * Arg:        targetdb [UNKN ] Undocumented argument [GenomicDB*]
 * Arg:              gp [UNKN ] Undocumented argument [GeneParser4Score *]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_serial_search_GeneLiteModel(Hscore * out,GeneWiseDB* querydb,GenomicDB* targetdb ,GeneParser4Score * gp);
#define serial_search_GeneLiteModel Wise2_serial_search_GeneLiteModel


/* Function:  PackAln_bestmemory_GeneLiteModel(query,target,gp,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_GeneLiteModel
 *
 *
 * Arg:         query [UNKN ] query data structure [GwLiteScore*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:            gp [UNKN ] Resource [GeneParser4Score *]
 * Arg:         dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:          dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_bestmemory_GeneLiteModel(GwLiteScore* query,ComplexSequence* target ,GeneParser4Score * gp,DPEnvelope * dpenv,DPRunImpl * dpri);
#define PackAln_bestmemory_GeneLiteModel Wise2_PackAln_bestmemory_GeneLiteModel


/* Function:  allocate_Expl_GeneLiteModel(query,target,gp,dpri)
 *
 * Descrip:    This function allocates the GeneLiteModel structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_GeneLiteModel_only
 *
 *
 * Arg:         query [UNKN ] query data structure [GwLiteScore*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:            gp [UNKN ] Resource [GeneParser4Score *]
 * Arg:          dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [GeneLiteModel *]
 *
 */
GeneLiteModel * Wise2_allocate_Expl_GeneLiteModel(GwLiteScore* query,ComplexSequence* target ,GeneParser4Score * gp,DPRunImpl * dpri);
#define allocate_Expl_GeneLiteModel Wise2_allocate_Expl_GeneLiteModel


/* Function:  recalculate_PackAln_GeneLiteModel(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by GeneLiteModel
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [GeneLiteModel *]
 *
 */
void Wise2_recalculate_PackAln_GeneLiteModel(PackAln * pal,GeneLiteModel * mat);
#define recalculate_PackAln_GeneLiteModel Wise2_recalculate_PackAln_GeneLiteModel


/* Function:  allocate_Small_GeneLiteModel(query,target,gp)
 *
 * Descrip:    This function allocates the GeneLiteModel structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_GeneLiteModel_only
 *
 *
 * Arg:         query [UNKN ] query data structure [GwLiteScore*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:            gp [UNKN ] Resource [GeneParser4Score *]
 *
 * Return [UNKN ]  Undocumented return value [GeneLiteModel *]
 *
 */
GeneLiteModel * Wise2_allocate_Small_GeneLiteModel(GwLiteScore* query,ComplexSequence* target ,GeneParser4Score * gp);
#define allocate_Small_GeneLiteModel Wise2_allocate_Small_GeneLiteModel


/* Function:  PackAln_calculate_Small_GeneLiteModel(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for GeneLiteModel structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_GeneLiteModel 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_GeneLiteModel 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [GeneLiteModel *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_calculate_Small_GeneLiteModel(GeneLiteModel * mat,DPEnvelope * dpenv);
#define PackAln_calculate_Small_GeneLiteModel Wise2_PackAln_calculate_Small_GeneLiteModel


/* Function:  AlnRangeSet_calculate_Small_GeneLiteModel(mat)
 *
 * Descrip:    This function calculates an alignment for GeneLiteModel structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_GeneLiteModel 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_GeneLiteModel
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_GeneLiteModel 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [GeneLiteModel *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_calculate_Small_GeneLiteModel(GeneLiteModel * mat);
#define AlnRangeSet_calculate_Small_GeneLiteModel Wise2_AlnRangeSet_calculate_Small_GeneLiteModel


/* Function:  AlnRangeSet_from_GeneLiteModel(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for GeneLiteModel structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_GeneLiteModel 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_GeneLiteModel
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [GeneLiteModel *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_from_GeneLiteModel(GeneLiteModel * mat);
#define AlnRangeSet_from_GeneLiteModel Wise2_AlnRangeSet_from_GeneLiteModel


/* Function:  convert_PackAln_to_AlnBlock_GeneLiteModel(pal)
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
AlnBlock * Wise2_convert_PackAln_to_AlnBlock_GeneLiteModel(PackAln * pal);
#define convert_PackAln_to_AlnBlock_GeneLiteModel Wise2_convert_PackAln_to_AlnBlock_GeneLiteModel


/* Function:  PackAln_read_Expl_GeneLiteModel(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [GeneLiteModel *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Expl_GeneLiteModel(GeneLiteModel * mat);
#define PackAln_read_Expl_GeneLiteModel Wise2_PackAln_read_Expl_GeneLiteModel


/* Function:  PackAln_read_generic_GeneLiteModel(mat,h)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [GeneLiteModel *]
 * Arg:          h [UNKN ] Undocumented argument [GeneLiteModel_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_generic_GeneLiteModel(GeneLiteModel * mat,GeneLiteModel_access_func_holder h);
#define PackAln_read_generic_GeneLiteModel Wise2_PackAln_read_generic_GeneLiteModel


/* Function:  calculate_GeneLiteModel(mat)
 *
 * Descrip:    This function calculates the GeneLiteModel matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_GeneLiteModel
 *
 *
 * Arg:        mat [UNKN ] GeneLiteModel which contains explicit basematrix memory [GeneLiteModel *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_GeneLiteModel(GeneLiteModel * mat);
#define calculate_GeneLiteModel Wise2_calculate_GeneLiteModel


/* Function:  calculate_dpenv_GeneLiteModel(mat,dpenv)
 *
 * Descrip:    This function calculates the GeneLiteModel matrix when in explicit mode, subject to the envelope
 *
 *
 * Arg:          mat [UNKN ] GeneLiteModel which contains explicit basematrix memory [GeneLiteModel *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_dpenv_GeneLiteModel(GeneLiteModel * mat,DPEnvelope * dpenv);
#define calculate_dpenv_GeneLiteModel Wise2_calculate_dpenv_GeneLiteModel


/* Function:  GeneLiteModel_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneLiteModel *]
 *
 */
GeneLiteModel * Wise2_GeneLiteModel_alloc(void);
#define GeneLiteModel_alloc Wise2_GeneLiteModel_alloc


/* Function:  free_GeneLiteModel(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneLiteModel *]
 *
 * Return [UNKN ]  Undocumented return value [GeneLiteModel *]
 *
 */
GeneLiteModel * Wise2_free_GeneLiteModel(GeneLiteModel * obj);
#define free_GeneLiteModel Wise2_free_GeneLiteModel


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
int Wise2_GeneLiteModel_shatter_access_main(GeneLiteModel * mat,int i,int j,int state);
#define GeneLiteModel_shatter_access_main Wise2_GeneLiteModel_shatter_access_main
int Wise2_GeneLiteModel_shatter_access_special(GeneLiteModel * mat,int i,int j,int state);
#define GeneLiteModel_shatter_access_special Wise2_GeneLiteModel_shatter_access_special
void * Wise2_thread_loop_GeneLiteModel(void * ptr);
#define thread_loop_GeneLiteModel Wise2_thread_loop_GeneLiteModel
int Wise2_score_only_GeneLiteModel(GwLiteScore* query,ComplexSequence* target ,GeneParser4Score * gp);
#define score_only_GeneLiteModel Wise2_score_only_GeneLiteModel
GeneLiteModel * Wise2_allocate_GeneLiteModel_only(GwLiteScore* query,ComplexSequence* target ,GeneParser4Score * gp);
#define allocate_GeneLiteModel_only Wise2_allocate_GeneLiteModel_only
void Wise2_init_GeneLiteModel(GeneLiteModel * mat);
#define init_GeneLiteModel Wise2_init_GeneLiteModel
AlnRange * Wise2_AlnRange_build_GeneLiteModel(GeneLiteModel * mat,int stopj,int stopspecstate,int * startj,int * startspecstate);
#define AlnRange_build_GeneLiteModel Wise2_AlnRange_build_GeneLiteModel
boolean Wise2_read_hidden_GeneLiteModel(GeneLiteModel * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out);
#define read_hidden_GeneLiteModel Wise2_read_hidden_GeneLiteModel
int Wise2_max_hidden_GeneLiteModel(GeneLiteModel * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_hidden_GeneLiteModel Wise2_max_hidden_GeneLiteModel
boolean Wise2_read_special_strip_GeneLiteModel(GeneLiteModel * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out);
#define read_special_strip_GeneLiteModel Wise2_read_special_strip_GeneLiteModel
int Wise2_max_special_strip_GeneLiteModel(GeneLiteModel * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_special_strip_GeneLiteModel Wise2_max_special_strip_GeneLiteModel
int Wise2_max_matrix_to_special_GeneLiteModel(GeneLiteModel * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_matrix_to_special_GeneLiteModel Wise2_max_matrix_to_special_GeneLiteModel
void Wise2_calculate_hidden_GeneLiteModel(GeneLiteModel * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv);
#define calculate_hidden_GeneLiteModel Wise2_calculate_hidden_GeneLiteModel
void Wise2_init_hidden_GeneLiteModel(GeneLiteModel * mat,int starti,int startj,int stopi,int stopj);
#define init_hidden_GeneLiteModel Wise2_init_hidden_GeneLiteModel
boolean Wise2_full_dc_GeneLiteModel(GeneLiteModel * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv);
#define full_dc_GeneLiteModel Wise2_full_dc_GeneLiteModel
boolean Wise2_do_dc_single_pass_GeneLiteModel(GeneLiteModel * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done);
#define do_dc_single_pass_GeneLiteModel Wise2_do_dc_single_pass_GeneLiteModel
void Wise2_push_dc_at_merge_GeneLiteModel(GeneLiteModel * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv);
#define push_dc_at_merge_GeneLiteModel Wise2_push_dc_at_merge_GeneLiteModel
void Wise2_follow_on_dc_GeneLiteModel(GeneLiteModel * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define follow_on_dc_GeneLiteModel Wise2_follow_on_dc_GeneLiteModel
void Wise2_run_up_dc_GeneLiteModel(GeneLiteModel * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define run_up_dc_GeneLiteModel Wise2_run_up_dc_GeneLiteModel
void Wise2_init_dc_GeneLiteModel(GeneLiteModel * mat);
#define init_dc_GeneLiteModel Wise2_init_dc_GeneLiteModel
int Wise2_start_end_find_end_GeneLiteModel(GeneLiteModel * mat,int * endj);
#define start_end_find_end_GeneLiteModel Wise2_start_end_find_end_GeneLiteModel
boolean Wise2_dc_optimised_start_end_calc_GeneLiteModel(GeneLiteModel *mat,DPEnvelope * dpenv);
#define dc_optimised_start_end_calc_GeneLiteModel Wise2_dc_optimised_start_end_calc_GeneLiteModel
void Wise2_init_start_end_linear_GeneLiteModel(GeneLiteModel * mat);
#define init_start_end_linear_GeneLiteModel Wise2_init_start_end_linear_GeneLiteModel
AlnConvertSet * Wise2_AlnConvertSet_GeneLiteModel(void);
#define AlnConvertSet_GeneLiteModel Wise2_AlnConvertSet_GeneLiteModel
int Wise2_GeneLiteModel_explicit_access_main(GeneLiteModel * mat,int i,int j,int state);
#define GeneLiteModel_explicit_access_main Wise2_GeneLiteModel_explicit_access_main
int Wise2_GeneLiteModel_explicit_access_special(GeneLiteModel * mat,int i,int j,int state);
#define GeneLiteModel_explicit_access_special Wise2_GeneLiteModel_explicit_access_special
int Wise2_find_end_GeneLiteModel(GeneLiteModel * mat,int * ri,int * rj,int * state,boolean * isspecial,GeneLiteModel_access_func_holder h);
#define find_end_GeneLiteModel Wise2_find_end_GeneLiteModel
void Wise2_GeneLiteModel_debug_show_matrix(GeneLiteModel * mat,int starti,int stopi,int startj,int stopj,FILE * ofp);
#define GeneLiteModel_debug_show_matrix Wise2_GeneLiteModel_debug_show_matrix
int Wise2_max_calc_GeneLiteModel(GeneLiteModel * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,GeneLiteModel_access_func_holder h);
#define max_calc_GeneLiteModel Wise2_max_calc_GeneLiteModel
int Wise2_max_calc_special_GeneLiteModel(GeneLiteModel * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,GeneLiteModel_access_func_holder h);
#define max_calc_special_GeneLiteModel Wise2_max_calc_special_GeneLiteModel

#ifdef _cplusplus
}
#endif

#endif
