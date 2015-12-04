#ifndef DYNAMITEestquick3HEADERFILE
#define DYNAMITEestquick3HEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"
#include "cdparser.h"
#include "genewisemodel.h"
#include "gwquickdb.h"



struct Wise2_EstQuick3 {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    BaseMatrix * basematrix;     
    ShatterMatrix * shatter;     
    int leni;    
    int lenj;    
    GeneWiseScoreFlat* query;    
    ComplexSequence* target;     
    cDNAParserScore * cp;    
    } ;  
/* EstQuick3 defined */ 
#ifndef DYNAMITE_DEFINED_EstQuick3
typedef struct Wise2_EstQuick3 Wise2_EstQuick3;
#define EstQuick3 Wise2_EstQuick3
#define DYNAMITE_DEFINED_EstQuick3
#endif


#ifdef PTHREAD
struct thread_pool_holder_EstQuick3 {  
    GeneWiseScoreFlat* query;   /* Query object placeholder */ 
    GeneWiseQuickDB* querydb;   /* Query database object */ 
    boolean query_init;  
    ComplexSequence* target;/* Target object placeholder */ 
    cDNADB* targetdb;   /* Target database object */ 
    boolean target_init; 
    cDNAParserScore * cp;    
    pthread_mutex_t input_lock;  
    pthread_mutex_t output_lock;     
    Hscore * out;    
    pthread_t * pool;    
    int number_of_threads;   
    boolean search_has_ended;    
    DBSearchImpl * dbsi;     
    } ;  
#endif /* PTHREAD */
struct Wise2_EstQuick3_access_func_holder {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int (*access_main)(EstQuick3*,int,int,int);  
    int (*access_special)(EstQuick3*,int,int,int);   
    } ;  
/* EstQuick3_access_func_holder defined */ 
#ifndef DYNAMITE_DEFINED_EstQuick3_access_func_holder
typedef struct Wise2_EstQuick3_access_func_holder Wise2_EstQuick3_access_func_holder;
#define EstQuick3_access_func_holder Wise2_EstQuick3_access_func_holder
#define DYNAMITE_DEFINED_EstQuick3_access_func_holder
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  PackAln_read_Shatter_EstQuick3(mat)
 *
 * Descrip:    Reads off PackAln from shatter matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [EstQuick3 *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Shatter_EstQuick3(EstQuick3 * mat);
#define PackAln_read_Shatter_EstQuick3 Wise2_PackAln_read_Shatter_EstQuick3


/* Function:  calculate_shatter_EstQuick3(mat,dpenv)
 *
 * Descrip:    This function calculates the EstQuick3 matrix when in shatter mode
 *
 *
 * Arg:          mat [UNKN ] (null) [EstQuick3 *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_shatter_EstQuick3(EstQuick3 * mat,DPEnvelope * dpenv);
#define calculate_shatter_EstQuick3 Wise2_calculate_shatter_EstQuick3


/* Function:  search_EstQuick3(dbsi,out,querydb,targetdb,cp)
 *
 * Descrip:    This function makes a database search of EstQuick3
 *             It uses the dbsi structure to choose which implementation
 *             to use of the database searching. This way at run time you
 *             can switch between single threaded/multi-threaded or hardware
 *
 *
 * Arg:            dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 * Arg:             out [UNKN ] Undocumented argument [Hscore *]
 * Arg:         querydb [UNKN ] Undocumented argument [GeneWiseQuickDB*]
 * Arg:        targetdb [UNKN ] Undocumented argument [cDNADB*]
 * Arg:              cp [UNKN ] Undocumented argument [cDNAParserScore *]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_search_EstQuick3(DBSearchImpl * dbsi,Hscore * out,GeneWiseQuickDB* querydb,cDNADB* targetdb ,cDNAParserScore * cp);
#define search_EstQuick3 Wise2_search_EstQuick3


/* Function:  serial_search_EstQuick3(out,querydb,targetdb,cp)
 *
 * Descrip:    This function makes a database search of EstQuick3
 *             It is a single processor implementation
 *
 *
 * Arg:             out [UNKN ] Undocumented argument [Hscore *]
 * Arg:         querydb [UNKN ] Undocumented argument [GeneWiseQuickDB*]
 * Arg:        targetdb [UNKN ] Undocumented argument [cDNADB*]
 * Arg:              cp [UNKN ] Undocumented argument [cDNAParserScore *]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_serial_search_EstQuick3(Hscore * out,GeneWiseQuickDB* querydb,cDNADB* targetdb ,cDNAParserScore * cp);
#define serial_search_EstQuick3 Wise2_serial_search_EstQuick3


/* Function:  PackAln_bestmemory_EstQuick3(query,target,cp,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_EstQuick3
 *
 *
 * Arg:         query [UNKN ] query data structure [GeneWiseScoreFlat*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:            cp [UNKN ] Resource [cDNAParserScore *]
 * Arg:         dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:          dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_bestmemory_EstQuick3(GeneWiseScoreFlat* query,ComplexSequence* target ,cDNAParserScore * cp,DPEnvelope * dpenv,DPRunImpl * dpri);
#define PackAln_bestmemory_EstQuick3 Wise2_PackAln_bestmemory_EstQuick3


/* Function:  allocate_Expl_EstQuick3(query,target,cp,dpri)
 *
 * Descrip:    This function allocates the EstQuick3 structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_EstQuick3_only
 *
 *
 * Arg:         query [UNKN ] query data structure [GeneWiseScoreFlat*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:            cp [UNKN ] Resource [cDNAParserScore *]
 * Arg:          dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [EstQuick3 *]
 *
 */
EstQuick3 * Wise2_allocate_Expl_EstQuick3(GeneWiseScoreFlat* query,ComplexSequence* target ,cDNAParserScore * cp,DPRunImpl * dpri);
#define allocate_Expl_EstQuick3 Wise2_allocate_Expl_EstQuick3


/* Function:  recalculate_PackAln_EstQuick3(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by EstQuick3
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [EstQuick3 *]
 *
 */
void Wise2_recalculate_PackAln_EstQuick3(PackAln * pal,EstQuick3 * mat);
#define recalculate_PackAln_EstQuick3 Wise2_recalculate_PackAln_EstQuick3


/* Function:  allocate_Small_EstQuick3(query,target,cp)
 *
 * Descrip:    This function allocates the EstQuick3 structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_EstQuick3_only
 *
 *
 * Arg:         query [UNKN ] query data structure [GeneWiseScoreFlat*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:            cp [UNKN ] Resource [cDNAParserScore *]
 *
 * Return [UNKN ]  Undocumented return value [EstQuick3 *]
 *
 */
EstQuick3 * Wise2_allocate_Small_EstQuick3(GeneWiseScoreFlat* query,ComplexSequence* target ,cDNAParserScore * cp);
#define allocate_Small_EstQuick3 Wise2_allocate_Small_EstQuick3


/* Function:  PackAln_calculate_Small_EstQuick3(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for EstQuick3 structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_EstQuick3 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_EstQuick3 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [EstQuick3 *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_calculate_Small_EstQuick3(EstQuick3 * mat,DPEnvelope * dpenv);
#define PackAln_calculate_Small_EstQuick3 Wise2_PackAln_calculate_Small_EstQuick3


/* Function:  AlnRangeSet_calculate_Small_EstQuick3(mat)
 *
 * Descrip:    This function calculates an alignment for EstQuick3 structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_EstQuick3 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_EstQuick3
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_EstQuick3 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [EstQuick3 *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_calculate_Small_EstQuick3(EstQuick3 * mat);
#define AlnRangeSet_calculate_Small_EstQuick3 Wise2_AlnRangeSet_calculate_Small_EstQuick3


/* Function:  AlnRangeSet_from_EstQuick3(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for EstQuick3 structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_EstQuick3 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_EstQuick3
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [EstQuick3 *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_from_EstQuick3(EstQuick3 * mat);
#define AlnRangeSet_from_EstQuick3 Wise2_AlnRangeSet_from_EstQuick3


/* Function:  convert_PackAln_to_AlnBlock_EstQuick3(pal)
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
AlnBlock * Wise2_convert_PackAln_to_AlnBlock_EstQuick3(PackAln * pal);
#define convert_PackAln_to_AlnBlock_EstQuick3 Wise2_convert_PackAln_to_AlnBlock_EstQuick3


/* Function:  PackAln_read_Expl_EstQuick3(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [EstQuick3 *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Expl_EstQuick3(EstQuick3 * mat);
#define PackAln_read_Expl_EstQuick3 Wise2_PackAln_read_Expl_EstQuick3


/* Function:  PackAln_read_generic_EstQuick3(mat,h)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [EstQuick3 *]
 * Arg:          h [UNKN ] Undocumented argument [EstQuick3_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_generic_EstQuick3(EstQuick3 * mat,EstQuick3_access_func_holder h);
#define PackAln_read_generic_EstQuick3 Wise2_PackAln_read_generic_EstQuick3


/* Function:  calculate_EstQuick3(mat)
 *
 * Descrip:    This function calculates the EstQuick3 matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_EstQuick3
 *
 *
 * Arg:        mat [UNKN ] EstQuick3 which contains explicit basematrix memory [EstQuick3 *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_EstQuick3(EstQuick3 * mat);
#define calculate_EstQuick3 Wise2_calculate_EstQuick3


/* Function:  calculate_dpenv_EstQuick3(mat,dpenv)
 *
 * Descrip:    This function calculates the EstQuick3 matrix when in explicit mode, subject to the envelope
 *
 *
 * Arg:          mat [UNKN ] EstQuick3 which contains explicit basematrix memory [EstQuick3 *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_dpenv_EstQuick3(EstQuick3 * mat,DPEnvelope * dpenv);
#define calculate_dpenv_EstQuick3 Wise2_calculate_dpenv_EstQuick3


/* Function:  EstQuick3_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [EstQuick3 *]
 *
 */
EstQuick3 * Wise2_EstQuick3_alloc(void);
#define EstQuick3_alloc Wise2_EstQuick3_alloc


/* Function:  free_EstQuick3(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [EstQuick3 *]
 *
 * Return [UNKN ]  Undocumented return value [EstQuick3 *]
 *
 */
EstQuick3 * Wise2_free_EstQuick3(EstQuick3 * obj);
#define free_EstQuick3 Wise2_free_EstQuick3


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
int Wise2_EstQuick3_shatter_access_main(EstQuick3 * mat,int i,int j,int state);
#define EstQuick3_shatter_access_main Wise2_EstQuick3_shatter_access_main
int Wise2_EstQuick3_shatter_access_special(EstQuick3 * mat,int i,int j,int state);
#define EstQuick3_shatter_access_special Wise2_EstQuick3_shatter_access_special
void * Wise2_thread_loop_EstQuick3(void * ptr);
#define thread_loop_EstQuick3 Wise2_thread_loop_EstQuick3
int Wise2_score_only_EstQuick3(GeneWiseScoreFlat* query,ComplexSequence* target ,cDNAParserScore * cp);
#define score_only_EstQuick3 Wise2_score_only_EstQuick3
EstQuick3 * Wise2_allocate_EstQuick3_only(GeneWiseScoreFlat* query,ComplexSequence* target ,cDNAParserScore * cp);
#define allocate_EstQuick3_only Wise2_allocate_EstQuick3_only
void Wise2_init_EstQuick3(EstQuick3 * mat);
#define init_EstQuick3 Wise2_init_EstQuick3
AlnRange * Wise2_AlnRange_build_EstQuick3(EstQuick3 * mat,int stopj,int stopspecstate,int * startj,int * startspecstate);
#define AlnRange_build_EstQuick3 Wise2_AlnRange_build_EstQuick3
boolean Wise2_read_hidden_EstQuick3(EstQuick3 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out);
#define read_hidden_EstQuick3 Wise2_read_hidden_EstQuick3
int Wise2_max_hidden_EstQuick3(EstQuick3 * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_hidden_EstQuick3 Wise2_max_hidden_EstQuick3
boolean Wise2_read_special_strip_EstQuick3(EstQuick3 * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out);
#define read_special_strip_EstQuick3 Wise2_read_special_strip_EstQuick3
int Wise2_max_special_strip_EstQuick3(EstQuick3 * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_special_strip_EstQuick3 Wise2_max_special_strip_EstQuick3
int Wise2_max_matrix_to_special_EstQuick3(EstQuick3 * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_matrix_to_special_EstQuick3 Wise2_max_matrix_to_special_EstQuick3
void Wise2_calculate_hidden_EstQuick3(EstQuick3 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv);
#define calculate_hidden_EstQuick3 Wise2_calculate_hidden_EstQuick3
void Wise2_init_hidden_EstQuick3(EstQuick3 * mat,int starti,int startj,int stopi,int stopj);
#define init_hidden_EstQuick3 Wise2_init_hidden_EstQuick3
boolean Wise2_full_dc_EstQuick3(EstQuick3 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv);
#define full_dc_EstQuick3 Wise2_full_dc_EstQuick3
boolean Wise2_do_dc_single_pass_EstQuick3(EstQuick3 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done);
#define do_dc_single_pass_EstQuick3 Wise2_do_dc_single_pass_EstQuick3
void Wise2_push_dc_at_merge_EstQuick3(EstQuick3 * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv);
#define push_dc_at_merge_EstQuick3 Wise2_push_dc_at_merge_EstQuick3
void Wise2_follow_on_dc_EstQuick3(EstQuick3 * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define follow_on_dc_EstQuick3 Wise2_follow_on_dc_EstQuick3
void Wise2_run_up_dc_EstQuick3(EstQuick3 * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define run_up_dc_EstQuick3 Wise2_run_up_dc_EstQuick3
void Wise2_init_dc_EstQuick3(EstQuick3 * mat);
#define init_dc_EstQuick3 Wise2_init_dc_EstQuick3
int Wise2_start_end_find_end_EstQuick3(EstQuick3 * mat,int * endj);
#define start_end_find_end_EstQuick3 Wise2_start_end_find_end_EstQuick3
boolean Wise2_dc_optimised_start_end_calc_EstQuick3(EstQuick3 *mat,DPEnvelope * dpenv);
#define dc_optimised_start_end_calc_EstQuick3 Wise2_dc_optimised_start_end_calc_EstQuick3
void Wise2_init_start_end_linear_EstQuick3(EstQuick3 * mat);
#define init_start_end_linear_EstQuick3 Wise2_init_start_end_linear_EstQuick3
AlnConvertSet * Wise2_AlnConvertSet_EstQuick3(void);
#define AlnConvertSet_EstQuick3 Wise2_AlnConvertSet_EstQuick3
int Wise2_EstQuick3_explicit_access_main(EstQuick3 * mat,int i,int j,int state);
#define EstQuick3_explicit_access_main Wise2_EstQuick3_explicit_access_main
int Wise2_EstQuick3_explicit_access_special(EstQuick3 * mat,int i,int j,int state);
#define EstQuick3_explicit_access_special Wise2_EstQuick3_explicit_access_special
int Wise2_find_end_EstQuick3(EstQuick3 * mat,int * ri,int * rj,int * state,boolean * isspecial,EstQuick3_access_func_holder h);
#define find_end_EstQuick3 Wise2_find_end_EstQuick3
void Wise2_EstQuick3_debug_show_matrix(EstQuick3 * mat,int starti,int stopi,int startj,int stopj,FILE * ofp);
#define EstQuick3_debug_show_matrix Wise2_EstQuick3_debug_show_matrix
int Wise2_max_calc_EstQuick3(EstQuick3 * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,EstQuick3_access_func_holder h);
#define max_calc_EstQuick3 Wise2_max_calc_EstQuick3
int Wise2_max_calc_special_EstQuick3(EstQuick3 * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,EstQuick3_access_func_holder h);
#define max_calc_special_EstQuick3 Wise2_max_calc_special_EstQuick3

#ifdef _cplusplus
}
#endif

#endif
