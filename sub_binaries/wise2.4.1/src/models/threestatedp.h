#ifndef DYNAMITEthreestatedpHEADERFILE
#define DYNAMITEthreestatedpHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "threestatemodel.h"


struct Wise2_ThreeStateProtein {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    BaseMatrix * basematrix;     
    ShatterMatrix * shatter;     
    int leni;    
    int lenj;    
    ThreeStateScore* query;  
    ComplexSequence* target;     
    } ;  
/* ThreeStateProtein defined */ 
#ifndef DYNAMITE_DEFINED_ThreeStateProtein
typedef struct Wise2_ThreeStateProtein Wise2_ThreeStateProtein;
#define ThreeStateProtein Wise2_ThreeStateProtein
#define DYNAMITE_DEFINED_ThreeStateProtein
#endif


#ifdef PTHREAD
struct thread_pool_holder_ThreeStateProtein {  
    ThreeStateScore* query; /* Static query data: never free'd */ 
    ComplexSequence* target;/* Target object placeholder */ 
    ProteinDB* targetdb;/* Target database object */ 
    boolean target_init; 
    pthread_mutex_t input_lock;  
    pthread_mutex_t output_lock;     
    Hscore * out;    
    pthread_t * pool;    
    int number_of_threads;   
    boolean search_has_ended;    
    DBSearchImpl * dbsi;     
    } ;  
#endif /* PTHREAD */
struct Wise2_ThreeStateProtein_access_func_holder {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int (*access_main)(ThreeStateProtein*,int,int,int);  
    int (*access_special)(ThreeStateProtein*,int,int,int);   
    } ;  
/* ThreeStateProtein_access_func_holder defined */ 
#ifndef DYNAMITE_DEFINED_ThreeStateProtein_access_func_holder
typedef struct Wise2_ThreeStateProtein_access_func_holder Wise2_ThreeStateProtein_access_func_holder;
#define ThreeStateProtein_access_func_holder Wise2_ThreeStateProtein_access_func_holder
#define DYNAMITE_DEFINED_ThreeStateProtein_access_func_holder
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  PackAln_read_Shatter_ThreeStateProtein(mat)
 *
 * Descrip:    Reads off PackAln from shatter matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [ThreeStateProtein *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Shatter_ThreeStateProtein(ThreeStateProtein * mat);
#define PackAln_read_Shatter_ThreeStateProtein Wise2_PackAln_read_Shatter_ThreeStateProtein


/* Function:  calculate_shatter_ThreeStateProtein(mat,dpenv)
 *
 * Descrip:    This function calculates the ThreeStateProtein matrix when in shatter mode
 *
 *
 * Arg:          mat [UNKN ] (null) [ThreeStateProtein *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_shatter_ThreeStateProtein(ThreeStateProtein * mat,DPEnvelope * dpenv);
#define calculate_shatter_ThreeStateProtein Wise2_calculate_shatter_ThreeStateProtein


/* Function:  search_ThreeStateProtein(dbsi,out,query,targetdb)
 *
 * Descrip:    This function makes a database search of ThreeStateProtein
 *             It uses the dbsi structure to choose which implementation
 *             to use of the database searching. This way at run time you
 *             can switch between single threaded/multi-threaded or hardware
 *
 *
 * Arg:            dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 * Arg:             out [UNKN ] Undocumented argument [Hscore *]
 * Arg:           query [UNKN ] Undocumented argument [ThreeStateScore*]
 * Arg:        targetdb [UNKN ] Undocumented argument [ProteinDB*]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_search_ThreeStateProtein(DBSearchImpl * dbsi,Hscore * out,ThreeStateScore* query,ProteinDB* targetdb );
#define search_ThreeStateProtein Wise2_search_ThreeStateProtein


/* Function:  serial_search_ThreeStateProtein(out,query,targetdb)
 *
 * Descrip:    This function makes a database search of ThreeStateProtein
 *             It is a single processor implementation
 *
 *
 * Arg:             out [UNKN ] Undocumented argument [Hscore *]
 * Arg:           query [UNKN ] Undocumented argument [ThreeStateScore*]
 * Arg:        targetdb [UNKN ] Undocumented argument [ProteinDB*]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_serial_search_ThreeStateProtein(Hscore * out,ThreeStateScore* query,ProteinDB* targetdb );
#define serial_search_ThreeStateProtein Wise2_serial_search_ThreeStateProtein


/* Function:  PackAln_bestmemory_ThreeStateProtein(query,target,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_ThreeStateProtein
 *
 *
 * Arg:         query [UNKN ] query data structure [ThreeStateScore*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:         dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:          dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_bestmemory_ThreeStateProtein(ThreeStateScore* query,ComplexSequence* target ,DPEnvelope * dpenv,DPRunImpl * dpri);
#define PackAln_bestmemory_ThreeStateProtein Wise2_PackAln_bestmemory_ThreeStateProtein


/* Function:  allocate_Expl_ThreeStateProtein(query,target,dpri)
 *
 * Descrip:    This function allocates the ThreeStateProtein structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_ThreeStateProtein_only
 *
 *
 * Arg:         query [UNKN ] query data structure [ThreeStateScore*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:          dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateProtein *]
 *
 */
ThreeStateProtein * Wise2_allocate_Expl_ThreeStateProtein(ThreeStateScore* query,ComplexSequence* target ,DPRunImpl * dpri);
#define allocate_Expl_ThreeStateProtein Wise2_allocate_Expl_ThreeStateProtein


/* Function:  recalculate_PackAln_ThreeStateProtein(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by ThreeStateProtein
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [ThreeStateProtein *]
 *
 */
void Wise2_recalculate_PackAln_ThreeStateProtein(PackAln * pal,ThreeStateProtein * mat);
#define recalculate_PackAln_ThreeStateProtein Wise2_recalculate_PackAln_ThreeStateProtein


/* Function:  allocate_Small_ThreeStateProtein(query,target)
 *
 * Descrip:    This function allocates the ThreeStateProtein structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_ThreeStateProtein_only
 *
 *
 * Arg:         query [UNKN ] query data structure [ThreeStateScore*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateProtein *]
 *
 */
ThreeStateProtein * Wise2_allocate_Small_ThreeStateProtein(ThreeStateScore* query,ComplexSequence* target );
#define allocate_Small_ThreeStateProtein Wise2_allocate_Small_ThreeStateProtein


/* Function:  PackAln_calculate_Small_ThreeStateProtein(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for ThreeStateProtein structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_ThreeStateProtein 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_ThreeStateProtein 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [ThreeStateProtein *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_calculate_Small_ThreeStateProtein(ThreeStateProtein * mat,DPEnvelope * dpenv);
#define PackAln_calculate_Small_ThreeStateProtein Wise2_PackAln_calculate_Small_ThreeStateProtein


/* Function:  AlnRangeSet_calculate_Small_ThreeStateProtein(mat)
 *
 * Descrip:    This function calculates an alignment for ThreeStateProtein structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_ThreeStateProtein 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_ThreeStateProtein
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_ThreeStateProtein 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [ThreeStateProtein *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_calculate_Small_ThreeStateProtein(ThreeStateProtein * mat);
#define AlnRangeSet_calculate_Small_ThreeStateProtein Wise2_AlnRangeSet_calculate_Small_ThreeStateProtein


/* Function:  AlnRangeSet_from_ThreeStateProtein(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for ThreeStateProtein structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_ThreeStateProtein 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_ThreeStateProtein
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [ThreeStateProtein *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_from_ThreeStateProtein(ThreeStateProtein * mat);
#define AlnRangeSet_from_ThreeStateProtein Wise2_AlnRangeSet_from_ThreeStateProtein


/* Function:  convert_PackAln_to_AlnBlock_ThreeStateProtein(pal)
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
AlnBlock * Wise2_convert_PackAln_to_AlnBlock_ThreeStateProtein(PackAln * pal);
#define convert_PackAln_to_AlnBlock_ThreeStateProtein Wise2_convert_PackAln_to_AlnBlock_ThreeStateProtein


/* Function:  PackAln_read_Expl_ThreeStateProtein(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [ThreeStateProtein *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Expl_ThreeStateProtein(ThreeStateProtein * mat);
#define PackAln_read_Expl_ThreeStateProtein Wise2_PackAln_read_Expl_ThreeStateProtein


/* Function:  PackAln_read_generic_ThreeStateProtein(mat,h)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [ThreeStateProtein *]
 * Arg:          h [UNKN ] Undocumented argument [ThreeStateProtein_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_generic_ThreeStateProtein(ThreeStateProtein * mat,ThreeStateProtein_access_func_holder h);
#define PackAln_read_generic_ThreeStateProtein Wise2_PackAln_read_generic_ThreeStateProtein


/* Function:  calculate_ThreeStateProtein(mat)
 *
 * Descrip:    This function calculates the ThreeStateProtein matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_ThreeStateProtein
 *
 *
 * Arg:        mat [UNKN ] ThreeStateProtein which contains explicit basematrix memory [ThreeStateProtein *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_ThreeStateProtein(ThreeStateProtein * mat);
#define calculate_ThreeStateProtein Wise2_calculate_ThreeStateProtein


/* Function:  calculate_dpenv_ThreeStateProtein(mat,dpenv)
 *
 * Descrip:    This function calculates the ThreeStateProtein matrix when in explicit mode, subject to the envelope
 *
 *
 * Arg:          mat [UNKN ] ThreeStateProtein which contains explicit basematrix memory [ThreeStateProtein *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_dpenv_ThreeStateProtein(ThreeStateProtein * mat,DPEnvelope * dpenv);
#define calculate_dpenv_ThreeStateProtein Wise2_calculate_dpenv_ThreeStateProtein


/* Function:  ThreeStateProtein_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateProtein *]
 *
 */
ThreeStateProtein * Wise2_ThreeStateProtein_alloc(void);
#define ThreeStateProtein_alloc Wise2_ThreeStateProtein_alloc


/* Function:  free_ThreeStateProtein(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ThreeStateProtein *]
 *
 * Return [UNKN ]  Undocumented return value [ThreeStateProtein *]
 *
 */
ThreeStateProtein * Wise2_free_ThreeStateProtein(ThreeStateProtein * obj);
#define free_ThreeStateProtein Wise2_free_ThreeStateProtein


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
int Wise2_ThreeStateProtein_shatter_access_main(ThreeStateProtein * mat,int i,int j,int state);
#define ThreeStateProtein_shatter_access_main Wise2_ThreeStateProtein_shatter_access_main
int Wise2_ThreeStateProtein_shatter_access_special(ThreeStateProtein * mat,int i,int j,int state);
#define ThreeStateProtein_shatter_access_special Wise2_ThreeStateProtein_shatter_access_special
void * Wise2_thread_loop_ThreeStateProtein(void * ptr);
#define thread_loop_ThreeStateProtein Wise2_thread_loop_ThreeStateProtein
int Wise2_score_only_ThreeStateProtein(ThreeStateScore* query,ComplexSequence* target );
#define score_only_ThreeStateProtein Wise2_score_only_ThreeStateProtein
ThreeStateProtein * Wise2_allocate_ThreeStateProtein_only(ThreeStateScore* query,ComplexSequence* target );
#define allocate_ThreeStateProtein_only Wise2_allocate_ThreeStateProtein_only
void Wise2_init_ThreeStateProtein(ThreeStateProtein * mat);
#define init_ThreeStateProtein Wise2_init_ThreeStateProtein
AlnRange * Wise2_AlnRange_build_ThreeStateProtein(ThreeStateProtein * mat,int stopj,int stopspecstate,int * startj,int * startspecstate);
#define AlnRange_build_ThreeStateProtein Wise2_AlnRange_build_ThreeStateProtein
boolean Wise2_read_hidden_ThreeStateProtein(ThreeStateProtein * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out);
#define read_hidden_ThreeStateProtein Wise2_read_hidden_ThreeStateProtein
int Wise2_max_hidden_ThreeStateProtein(ThreeStateProtein * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_hidden_ThreeStateProtein Wise2_max_hidden_ThreeStateProtein
boolean Wise2_read_special_strip_ThreeStateProtein(ThreeStateProtein * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out);
#define read_special_strip_ThreeStateProtein Wise2_read_special_strip_ThreeStateProtein
int Wise2_max_special_strip_ThreeStateProtein(ThreeStateProtein * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_special_strip_ThreeStateProtein Wise2_max_special_strip_ThreeStateProtein
int Wise2_max_matrix_to_special_ThreeStateProtein(ThreeStateProtein * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_matrix_to_special_ThreeStateProtein Wise2_max_matrix_to_special_ThreeStateProtein
void Wise2_calculate_hidden_ThreeStateProtein(ThreeStateProtein * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv);
#define calculate_hidden_ThreeStateProtein Wise2_calculate_hidden_ThreeStateProtein
void Wise2_init_hidden_ThreeStateProtein(ThreeStateProtein * mat,int starti,int startj,int stopi,int stopj);
#define init_hidden_ThreeStateProtein Wise2_init_hidden_ThreeStateProtein
boolean Wise2_full_dc_ThreeStateProtein(ThreeStateProtein * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv);
#define full_dc_ThreeStateProtein Wise2_full_dc_ThreeStateProtein
boolean Wise2_do_dc_single_pass_ThreeStateProtein(ThreeStateProtein * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done);
#define do_dc_single_pass_ThreeStateProtein Wise2_do_dc_single_pass_ThreeStateProtein
void Wise2_push_dc_at_merge_ThreeStateProtein(ThreeStateProtein * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv);
#define push_dc_at_merge_ThreeStateProtein Wise2_push_dc_at_merge_ThreeStateProtein
void Wise2_follow_on_dc_ThreeStateProtein(ThreeStateProtein * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define follow_on_dc_ThreeStateProtein Wise2_follow_on_dc_ThreeStateProtein
void Wise2_run_up_dc_ThreeStateProtein(ThreeStateProtein * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define run_up_dc_ThreeStateProtein Wise2_run_up_dc_ThreeStateProtein
void Wise2_init_dc_ThreeStateProtein(ThreeStateProtein * mat);
#define init_dc_ThreeStateProtein Wise2_init_dc_ThreeStateProtein
int Wise2_start_end_find_end_ThreeStateProtein(ThreeStateProtein * mat,int * endj);
#define start_end_find_end_ThreeStateProtein Wise2_start_end_find_end_ThreeStateProtein
boolean Wise2_dc_optimised_start_end_calc_ThreeStateProtein(ThreeStateProtein *mat,DPEnvelope * dpenv);
#define dc_optimised_start_end_calc_ThreeStateProtein Wise2_dc_optimised_start_end_calc_ThreeStateProtein
void Wise2_init_start_end_linear_ThreeStateProtein(ThreeStateProtein * mat);
#define init_start_end_linear_ThreeStateProtein Wise2_init_start_end_linear_ThreeStateProtein
AlnConvertSet * Wise2_AlnConvertSet_ThreeStateProtein(void);
#define AlnConvertSet_ThreeStateProtein Wise2_AlnConvertSet_ThreeStateProtein
int Wise2_ThreeStateProtein_explicit_access_main(ThreeStateProtein * mat,int i,int j,int state);
#define ThreeStateProtein_explicit_access_main Wise2_ThreeStateProtein_explicit_access_main
int Wise2_ThreeStateProtein_explicit_access_special(ThreeStateProtein * mat,int i,int j,int state);
#define ThreeStateProtein_explicit_access_special Wise2_ThreeStateProtein_explicit_access_special
int Wise2_find_end_ThreeStateProtein(ThreeStateProtein * mat,int * ri,int * rj,int * state,boolean * isspecial,ThreeStateProtein_access_func_holder h);
#define find_end_ThreeStateProtein Wise2_find_end_ThreeStateProtein
void Wise2_ThreeStateProtein_debug_show_matrix(ThreeStateProtein * mat,int starti,int stopi,int startj,int stopj,FILE * ofp);
#define ThreeStateProtein_debug_show_matrix Wise2_ThreeStateProtein_debug_show_matrix
int Wise2_max_calc_ThreeStateProtein(ThreeStateProtein * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,ThreeStateProtein_access_func_holder h);
#define max_calc_ThreeStateProtein Wise2_max_calc_ThreeStateProtein
int Wise2_max_calc_special_ThreeStateProtein(ThreeStateProtein * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,ThreeStateProtein_access_func_holder h);
#define max_calc_special_ThreeStateProtein Wise2_max_calc_special_ThreeStateProtein

#ifdef _cplusplus
}
#endif

#endif
