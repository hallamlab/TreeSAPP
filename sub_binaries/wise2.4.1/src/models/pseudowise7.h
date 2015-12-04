#ifndef DYNAMITEpseudowise7HEADERFILE
#define DYNAMITEpseudowise7HEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"
#include "geneparser4.h"
#include "genewisemodel.h"
#include "genewisemodeldb.h"

struct Wise2_PseudoWise7 {  
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
    Score dna_ext;   
    } ;  
/* PseudoWise7 defined */ 
#ifndef DYNAMITE_DEFINED_PseudoWise7
typedef struct Wise2_PseudoWise7 Wise2_PseudoWise7;
#define PseudoWise7 Wise2_PseudoWise7
#define DYNAMITE_DEFINED_PseudoWise7
#endif


#ifdef PTHREAD
struct thread_pool_holder_PseudoWise7 {  
    GeneWiseScore* query;   /* Query object placeholder */ 
    GeneWiseDB* querydb;/* Query database object */ 
    boolean query_init;  
    ComplexSequence* target;/* Target object placeholder */ 
    GenomicDB* targetdb;/* Target database object */ 
    boolean target_init; 
    GeneParser4Score * gp;   
    Score dna_ext;   
    pthread_mutex_t input_lock;  
    pthread_mutex_t output_lock;     
    Hscore * out;    
    pthread_t * pool;    
    int number_of_threads;   
    boolean search_has_ended;    
    DBSearchImpl * dbsi;     
    } ;  
#endif /* PTHREAD */
struct Wise2_PseudoWise7_access_func_holder {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int (*access_main)(PseudoWise7*,int,int,int);    
    int (*access_special)(PseudoWise7*,int,int,int); 
    } ;  
/* PseudoWise7_access_func_holder defined */ 
#ifndef DYNAMITE_DEFINED_PseudoWise7_access_func_holder
typedef struct Wise2_PseudoWise7_access_func_holder Wise2_PseudoWise7_access_func_holder;
#define PseudoWise7_access_func_holder Wise2_PseudoWise7_access_func_holder
#define DYNAMITE_DEFINED_PseudoWise7_access_func_holder
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  PackAln_read_Shatter_PseudoWise7(mat)
 *
 * Descrip:    Reads off PackAln from shatter matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [PseudoWise7 *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Shatter_PseudoWise7(PseudoWise7 * mat);
#define PackAln_read_Shatter_PseudoWise7 Wise2_PackAln_read_Shatter_PseudoWise7


/* Function:  calculate_shatter_PseudoWise7(mat,dpenv)
 *
 * Descrip:    This function calculates the PseudoWise7 matrix when in shatter mode
 *
 *
 * Arg:          mat [UNKN ] (null) [PseudoWise7 *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_shatter_PseudoWise7(PseudoWise7 * mat,DPEnvelope * dpenv);
#define calculate_shatter_PseudoWise7 Wise2_calculate_shatter_PseudoWise7


/* Function:  search_PseudoWise7(dbsi,out,querydb,targetdb,gp,dna_ext)
 *
 * Descrip:    This function makes a database search of PseudoWise7
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
 * Arg:         dna_ext [UNKN ] Undocumented argument [Score]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_search_PseudoWise7(DBSearchImpl * dbsi,Hscore * out,GeneWiseDB* querydb,GenomicDB* targetdb ,GeneParser4Score * gp,Score dna_ext);
#define search_PseudoWise7 Wise2_search_PseudoWise7


/* Function:  serial_search_PseudoWise7(out,querydb,targetdb,gp,dna_ext)
 *
 * Descrip:    This function makes a database search of PseudoWise7
 *             It is a single processor implementation
 *
 *
 * Arg:             out [UNKN ] Undocumented argument [Hscore *]
 * Arg:         querydb [UNKN ] Undocumented argument [GeneWiseDB*]
 * Arg:        targetdb [UNKN ] Undocumented argument [GenomicDB*]
 * Arg:              gp [UNKN ] Undocumented argument [GeneParser4Score *]
 * Arg:         dna_ext [UNKN ] Undocumented argument [Score]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_serial_search_PseudoWise7(Hscore * out,GeneWiseDB* querydb,GenomicDB* targetdb ,GeneParser4Score * gp,Score dna_ext);
#define serial_search_PseudoWise7 Wise2_serial_search_PseudoWise7


/* Function:  PackAln_bestmemory_PseudoWise7(query,target,gp,dna_ext,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_PseudoWise7
 *
 *
 * Arg:          query [UNKN ] query data structure [GeneWiseScore*]
 * Arg:         target [UNKN ] target data structure [ComplexSequence*]
 * Arg:             gp [UNKN ] Resource [GeneParser4Score *]
 * Arg:        dna_ext [UNKN ] Resource [Score]
 * Arg:          dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:           dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_bestmemory_PseudoWise7(GeneWiseScore* query,ComplexSequence* target ,GeneParser4Score * gp,Score dna_ext,DPEnvelope * dpenv,DPRunImpl * dpri);
#define PackAln_bestmemory_PseudoWise7 Wise2_PackAln_bestmemory_PseudoWise7


/* Function:  allocate_Expl_PseudoWise7(query,target,gp,dna_ext,dpri)
 *
 * Descrip:    This function allocates the PseudoWise7 structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_PseudoWise7_only
 *
 *
 * Arg:          query [UNKN ] query data structure [GeneWiseScore*]
 * Arg:         target [UNKN ] target data structure [ComplexSequence*]
 * Arg:             gp [UNKN ] Resource [GeneParser4Score *]
 * Arg:        dna_ext [UNKN ] Resource [Score]
 * Arg:           dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [PseudoWise7 *]
 *
 */
PseudoWise7 * Wise2_allocate_Expl_PseudoWise7(GeneWiseScore* query,ComplexSequence* target ,GeneParser4Score * gp,Score dna_ext,DPRunImpl * dpri);
#define allocate_Expl_PseudoWise7 Wise2_allocate_Expl_PseudoWise7


/* Function:  recalculate_PackAln_PseudoWise7(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by PseudoWise7
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [PseudoWise7 *]
 *
 */
void Wise2_recalculate_PackAln_PseudoWise7(PackAln * pal,PseudoWise7 * mat);
#define recalculate_PackAln_PseudoWise7 Wise2_recalculate_PackAln_PseudoWise7


/* Function:  allocate_Small_PseudoWise7(query,target,gp,dna_ext)
 *
 * Descrip:    This function allocates the PseudoWise7 structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_PseudoWise7_only
 *
 *
 * Arg:          query [UNKN ] query data structure [GeneWiseScore*]
 * Arg:         target [UNKN ] target data structure [ComplexSequence*]
 * Arg:             gp [UNKN ] Resource [GeneParser4Score *]
 * Arg:        dna_ext [UNKN ] Resource [Score]
 *
 * Return [UNKN ]  Undocumented return value [PseudoWise7 *]
 *
 */
PseudoWise7 * Wise2_allocate_Small_PseudoWise7(GeneWiseScore* query,ComplexSequence* target ,GeneParser4Score * gp,Score dna_ext);
#define allocate_Small_PseudoWise7 Wise2_allocate_Small_PseudoWise7


/* Function:  PackAln_calculate_Small_PseudoWise7(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for PseudoWise7 structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_PseudoWise7 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_PseudoWise7 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [PseudoWise7 *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_calculate_Small_PseudoWise7(PseudoWise7 * mat,DPEnvelope * dpenv);
#define PackAln_calculate_Small_PseudoWise7 Wise2_PackAln_calculate_Small_PseudoWise7


/* Function:  AlnRangeSet_calculate_Small_PseudoWise7(mat)
 *
 * Descrip:    This function calculates an alignment for PseudoWise7 structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_PseudoWise7 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_PseudoWise7
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_PseudoWise7 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [PseudoWise7 *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_calculate_Small_PseudoWise7(PseudoWise7 * mat);
#define AlnRangeSet_calculate_Small_PseudoWise7 Wise2_AlnRangeSet_calculate_Small_PseudoWise7


/* Function:  AlnRangeSet_from_PseudoWise7(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for PseudoWise7 structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_PseudoWise7 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_PseudoWise7
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [PseudoWise7 *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_from_PseudoWise7(PseudoWise7 * mat);
#define AlnRangeSet_from_PseudoWise7 Wise2_AlnRangeSet_from_PseudoWise7


/* Function:  convert_PackAln_to_AlnBlock_PseudoWise7(pal)
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
AlnBlock * Wise2_convert_PackAln_to_AlnBlock_PseudoWise7(PackAln * pal);
#define convert_PackAln_to_AlnBlock_PseudoWise7 Wise2_convert_PackAln_to_AlnBlock_PseudoWise7


/* Function:  PackAln_read_Expl_PseudoWise7(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [PseudoWise7 *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Expl_PseudoWise7(PseudoWise7 * mat);
#define PackAln_read_Expl_PseudoWise7 Wise2_PackAln_read_Expl_PseudoWise7


/* Function:  PackAln_read_generic_PseudoWise7(mat,h)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [PseudoWise7 *]
 * Arg:          h [UNKN ] Undocumented argument [PseudoWise7_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_generic_PseudoWise7(PseudoWise7 * mat,PseudoWise7_access_func_holder h);
#define PackAln_read_generic_PseudoWise7 Wise2_PackAln_read_generic_PseudoWise7


/* Function:  calculate_PseudoWise7(mat)
 *
 * Descrip:    This function calculates the PseudoWise7 matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_PseudoWise7
 *
 *
 * Arg:        mat [UNKN ] PseudoWise7 which contains explicit basematrix memory [PseudoWise7 *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_PseudoWise7(PseudoWise7 * mat);
#define calculate_PseudoWise7 Wise2_calculate_PseudoWise7


/* Function:  calculate_dpenv_PseudoWise7(mat,dpenv)
 *
 * Descrip:    This function calculates the PseudoWise7 matrix when in explicit mode, subject to the envelope
 *
 *
 * Arg:          mat [UNKN ] PseudoWise7 which contains explicit basematrix memory [PseudoWise7 *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_dpenv_PseudoWise7(PseudoWise7 * mat,DPEnvelope * dpenv);
#define calculate_dpenv_PseudoWise7 Wise2_calculate_dpenv_PseudoWise7


/* Function:  PseudoWise7_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PseudoWise7 *]
 *
 */
PseudoWise7 * Wise2_PseudoWise7_alloc(void);
#define PseudoWise7_alloc Wise2_PseudoWise7_alloc


/* Function:  free_PseudoWise7(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [PseudoWise7 *]
 *
 * Return [UNKN ]  Undocumented return value [PseudoWise7 *]
 *
 */
PseudoWise7 * Wise2_free_PseudoWise7(PseudoWise7 * obj);
#define free_PseudoWise7 Wise2_free_PseudoWise7


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
int Wise2_PseudoWise7_shatter_access_main(PseudoWise7 * mat,int i,int j,int state);
#define PseudoWise7_shatter_access_main Wise2_PseudoWise7_shatter_access_main
int Wise2_PseudoWise7_shatter_access_special(PseudoWise7 * mat,int i,int j,int state);
#define PseudoWise7_shatter_access_special Wise2_PseudoWise7_shatter_access_special
void * Wise2_thread_loop_PseudoWise7(void * ptr);
#define thread_loop_PseudoWise7 Wise2_thread_loop_PseudoWise7
int Wise2_score_only_PseudoWise7(GeneWiseScore* query,ComplexSequence* target ,GeneParser4Score * gp,Score dna_ext);
#define score_only_PseudoWise7 Wise2_score_only_PseudoWise7
PseudoWise7 * Wise2_allocate_PseudoWise7_only(GeneWiseScore* query,ComplexSequence* target ,GeneParser4Score * gp,Score dna_ext);
#define allocate_PseudoWise7_only Wise2_allocate_PseudoWise7_only
void Wise2_init_PseudoWise7(PseudoWise7 * mat);
#define init_PseudoWise7 Wise2_init_PseudoWise7
AlnRange * Wise2_AlnRange_build_PseudoWise7(PseudoWise7 * mat,int stopj,int stopspecstate,int * startj,int * startspecstate);
#define AlnRange_build_PseudoWise7 Wise2_AlnRange_build_PseudoWise7
boolean Wise2_read_hidden_PseudoWise7(PseudoWise7 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out);
#define read_hidden_PseudoWise7 Wise2_read_hidden_PseudoWise7
int Wise2_max_hidden_PseudoWise7(PseudoWise7 * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_hidden_PseudoWise7 Wise2_max_hidden_PseudoWise7
boolean Wise2_read_special_strip_PseudoWise7(PseudoWise7 * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out);
#define read_special_strip_PseudoWise7 Wise2_read_special_strip_PseudoWise7
int Wise2_max_special_strip_PseudoWise7(PseudoWise7 * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_special_strip_PseudoWise7 Wise2_max_special_strip_PseudoWise7
int Wise2_max_matrix_to_special_PseudoWise7(PseudoWise7 * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_matrix_to_special_PseudoWise7 Wise2_max_matrix_to_special_PseudoWise7
void Wise2_calculate_hidden_PseudoWise7(PseudoWise7 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv);
#define calculate_hidden_PseudoWise7 Wise2_calculate_hidden_PseudoWise7
void Wise2_init_hidden_PseudoWise7(PseudoWise7 * mat,int starti,int startj,int stopi,int stopj);
#define init_hidden_PseudoWise7 Wise2_init_hidden_PseudoWise7
boolean Wise2_full_dc_PseudoWise7(PseudoWise7 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv);
#define full_dc_PseudoWise7 Wise2_full_dc_PseudoWise7
boolean Wise2_do_dc_single_pass_PseudoWise7(PseudoWise7 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done);
#define do_dc_single_pass_PseudoWise7 Wise2_do_dc_single_pass_PseudoWise7
void Wise2_push_dc_at_merge_PseudoWise7(PseudoWise7 * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv);
#define push_dc_at_merge_PseudoWise7 Wise2_push_dc_at_merge_PseudoWise7
void Wise2_follow_on_dc_PseudoWise7(PseudoWise7 * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define follow_on_dc_PseudoWise7 Wise2_follow_on_dc_PseudoWise7
void Wise2_run_up_dc_PseudoWise7(PseudoWise7 * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define run_up_dc_PseudoWise7 Wise2_run_up_dc_PseudoWise7
void Wise2_init_dc_PseudoWise7(PseudoWise7 * mat);
#define init_dc_PseudoWise7 Wise2_init_dc_PseudoWise7
int Wise2_start_end_find_end_PseudoWise7(PseudoWise7 * mat,int * endj);
#define start_end_find_end_PseudoWise7 Wise2_start_end_find_end_PseudoWise7
boolean Wise2_dc_optimised_start_end_calc_PseudoWise7(PseudoWise7 *mat,DPEnvelope * dpenv);
#define dc_optimised_start_end_calc_PseudoWise7 Wise2_dc_optimised_start_end_calc_PseudoWise7
void Wise2_init_start_end_linear_PseudoWise7(PseudoWise7 * mat);
#define init_start_end_linear_PseudoWise7 Wise2_init_start_end_linear_PseudoWise7
AlnConvertSet * Wise2_AlnConvertSet_PseudoWise7(void);
#define AlnConvertSet_PseudoWise7 Wise2_AlnConvertSet_PseudoWise7
int Wise2_PseudoWise7_explicit_access_main(PseudoWise7 * mat,int i,int j,int state);
#define PseudoWise7_explicit_access_main Wise2_PseudoWise7_explicit_access_main
int Wise2_PseudoWise7_explicit_access_special(PseudoWise7 * mat,int i,int j,int state);
#define PseudoWise7_explicit_access_special Wise2_PseudoWise7_explicit_access_special
int Wise2_find_end_PseudoWise7(PseudoWise7 * mat,int * ri,int * rj,int * state,boolean * isspecial,PseudoWise7_access_func_holder h);
#define find_end_PseudoWise7 Wise2_find_end_PseudoWise7
void Wise2_PseudoWise7_debug_show_matrix(PseudoWise7 * mat,int starti,int stopi,int startj,int stopj,FILE * ofp);
#define PseudoWise7_debug_show_matrix Wise2_PseudoWise7_debug_show_matrix
int Wise2_max_calc_PseudoWise7(PseudoWise7 * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,PseudoWise7_access_func_holder h);
#define max_calc_PseudoWise7 Wise2_max_calc_PseudoWise7
int Wise2_max_calc_special_PseudoWise7(PseudoWise7 * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,PseudoWise7_access_func_holder h);
#define max_calc_special_PseudoWise7 Wise2_max_calc_special_PseudoWise7

#ifdef _cplusplus
}
#endif

#endif
