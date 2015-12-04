#ifndef DYNAMITEalignwisedpHEADERFILE
#define DYNAMITEalignwisedpHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "aligngenemodel.h"



struct Wise2_AlignGeneModelFrame {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int len;     
    } ;  
/* AlignGeneModelFrame defined */ 
#ifndef DYNAMITE_DEFINED_AlignGeneModelFrame
typedef struct Wise2_AlignGeneModelFrame Wise2_AlignGeneModelFrame;
#define AlignGeneModelFrame Wise2_AlignGeneModelFrame
#define DYNAMITE_DEFINED_AlignGeneModelFrame
#endif


struct Wise2_AlignWise {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    BaseMatrix * basematrix;     
    ShatterMatrix * shatter;     
    int leni;    
    int lenj;    
    AlignGeneModelFrame* model;  
    AlignGeneModelScore* align;  
    Score intronopen;    
    Score geneopen;  
    } ;  
/* AlignWise defined */ 
#ifndef DYNAMITE_DEFINED_AlignWise
typedef struct Wise2_AlignWise Wise2_AlignWise;
#define AlignWise Wise2_AlignWise
#define DYNAMITE_DEFINED_AlignWise
#endif


#ifdef PTHREAD
struct thread_pool_holder_AlignWise {  
    AlignGeneModelFrame* model; /* Static query data: never free'd */ 
    AlignGeneModelScore* align; /* Static target data: never free'd */ 
    Score intronopen;    
    Score geneopen;  
    pthread_mutex_t input_lock;  
    pthread_mutex_t output_lock;     
    Hscore * out;    
    pthread_t * pool;    
    int number_of_threads;   
    boolean search_has_ended;    
    DBSearchImpl * dbsi;     
    } ;  
#endif /* PTHREAD */
struct Wise2_AlignWise_access_func_holder {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int (*access_main)(AlignWise*,int,int,int);  
    int (*access_special)(AlignWise*,int,int,int);   
    } ;  
/* AlignWise_access_func_holder defined */ 
#ifndef DYNAMITE_DEFINED_AlignWise_access_func_holder
typedef struct Wise2_AlignWise_access_func_holder Wise2_AlignWise_access_func_holder;
#define AlignWise_access_func_holder Wise2_AlignWise_access_func_holder
#define DYNAMITE_DEFINED_AlignWise_access_func_holder
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_AlignGeneModelFrame(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AlignGeneModelFrame *]
 *
 * Return [UNKN ]  Undocumented return value [AlignGeneModelFrame *]
 *
 */
AlignGeneModelFrame * Wise2_hard_link_AlignGeneModelFrame(AlignGeneModelFrame * obj);
#define hard_link_AlignGeneModelFrame Wise2_hard_link_AlignGeneModelFrame


/* Function:  AlignGeneModelFrame_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AlignGeneModelFrame *]
 *
 */
AlignGeneModelFrame * Wise2_AlignGeneModelFrame_alloc(void);
#define AlignGeneModelFrame_alloc Wise2_AlignGeneModelFrame_alloc


/* Function:  free_AlignGeneModelFrame(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AlignGeneModelFrame *]
 *
 * Return [UNKN ]  Undocumented return value [AlignGeneModelFrame *]
 *
 */
AlignGeneModelFrame * Wise2_free_AlignGeneModelFrame(AlignGeneModelFrame * obj);
#define free_AlignGeneModelFrame Wise2_free_AlignGeneModelFrame


/* Function:  PackAln_read_Shatter_AlignWise(mat)
 *
 * Descrip:    Reads off PackAln from shatter matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [AlignWise *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Shatter_AlignWise(AlignWise * mat);
#define PackAln_read_Shatter_AlignWise Wise2_PackAln_read_Shatter_AlignWise


/* Function:  calculate_shatter_AlignWise(mat,dpenv)
 *
 * Descrip:    This function calculates the AlignWise matrix when in shatter mode
 *
 *
 * Arg:          mat [UNKN ] (null) [AlignWise *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_shatter_AlignWise(AlignWise * mat,DPEnvelope * dpenv);
#define calculate_shatter_AlignWise Wise2_calculate_shatter_AlignWise


/* Function:  search_AlignWise(dbsi,out,model,align,intronopen,geneopen)
 *
 * Descrip:    This function makes a database search of AlignWise
 *             It uses the dbsi structure to choose which implementation
 *             to use of the database searching. This way at run time you
 *             can switch between single threaded/multi-threaded or hardware
 *
 *
 * Arg:              dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 * Arg:               out [UNKN ] Undocumented argument [Hscore *]
 * Arg:             model [UNKN ] Undocumented argument [AlignGeneModelFrame*]
 * Arg:             align [UNKN ] Undocumented argument [AlignGeneModelScore*]
 * Arg:        intronopen [UNKN ] Undocumented argument [Score]
 * Arg:          geneopen [UNKN ] Undocumented argument [Score]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_search_AlignWise(DBSearchImpl * dbsi,Hscore * out,AlignGeneModelFrame* model,AlignGeneModelScore* align ,Score intronopen,Score geneopen);
#define search_AlignWise Wise2_search_AlignWise


/* Function:  serial_search_AlignWise(out,model,align,intronopen,geneopen)
 *
 * Descrip:    This function makes a database search of AlignWise
 *             It is a single processor implementation
 *
 *
 * Arg:               out [UNKN ] Undocumented argument [Hscore *]
 * Arg:             model [UNKN ] Undocumented argument [AlignGeneModelFrame*]
 * Arg:             align [UNKN ] Undocumented argument [AlignGeneModelScore*]
 * Arg:        intronopen [UNKN ] Undocumented argument [Score]
 * Arg:          geneopen [UNKN ] Undocumented argument [Score]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_serial_search_AlignWise(Hscore * out,AlignGeneModelFrame* model,AlignGeneModelScore* align ,Score intronopen,Score geneopen);
#define serial_search_AlignWise Wise2_serial_search_AlignWise


/* Function:  PackAln_bestmemory_AlignWise(model,align,intronopen,geneopen,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_AlignWise
 *
 *
 * Arg:             model [UNKN ] query data structure [AlignGeneModelFrame*]
 * Arg:             align [UNKN ] target data structure [AlignGeneModelScore*]
 * Arg:        intronopen [UNKN ] Resource [Score]
 * Arg:          geneopen [UNKN ] Resource [Score]
 * Arg:             dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:              dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_bestmemory_AlignWise(AlignGeneModelFrame* model,AlignGeneModelScore* align ,Score intronopen,Score geneopen,DPEnvelope * dpenv,DPRunImpl * dpri);
#define PackAln_bestmemory_AlignWise Wise2_PackAln_bestmemory_AlignWise


/* Function:  allocate_Expl_AlignWise(model,align,intronopen,geneopen,dpri)
 *
 * Descrip:    This function allocates the AlignWise structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_AlignWise_only
 *
 *
 * Arg:             model [UNKN ] query data structure [AlignGeneModelFrame*]
 * Arg:             align [UNKN ] target data structure [AlignGeneModelScore*]
 * Arg:        intronopen [UNKN ] Resource [Score]
 * Arg:          geneopen [UNKN ] Resource [Score]
 * Arg:              dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [AlignWise *]
 *
 */
AlignWise * Wise2_allocate_Expl_AlignWise(AlignGeneModelFrame* model,AlignGeneModelScore* align ,Score intronopen,Score geneopen,DPRunImpl * dpri);
#define allocate_Expl_AlignWise Wise2_allocate_Expl_AlignWise


/* Function:  recalculate_PackAln_AlignWise(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by AlignWise
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [AlignWise *]
 *
 */
void Wise2_recalculate_PackAln_AlignWise(PackAln * pal,AlignWise * mat);
#define recalculate_PackAln_AlignWise Wise2_recalculate_PackAln_AlignWise


/* Function:  allocate_Small_AlignWise(model,align,intronopen,geneopen)
 *
 * Descrip:    This function allocates the AlignWise structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_AlignWise_only
 *
 *
 * Arg:             model [UNKN ] query data structure [AlignGeneModelFrame*]
 * Arg:             align [UNKN ] target data structure [AlignGeneModelScore*]
 * Arg:        intronopen [UNKN ] Resource [Score]
 * Arg:          geneopen [UNKN ] Resource [Score]
 *
 * Return [UNKN ]  Undocumented return value [AlignWise *]
 *
 */
AlignWise * Wise2_allocate_Small_AlignWise(AlignGeneModelFrame* model,AlignGeneModelScore* align ,Score intronopen,Score geneopen);
#define allocate_Small_AlignWise Wise2_allocate_Small_AlignWise


/* Function:  PackAln_calculate_Small_AlignWise(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for AlignWise structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_AlignWise 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_AlignWise 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [AlignWise *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_calculate_Small_AlignWise(AlignWise * mat,DPEnvelope * dpenv);
#define PackAln_calculate_Small_AlignWise Wise2_PackAln_calculate_Small_AlignWise


/* Function:  AlnRangeSet_calculate_Small_AlignWise(mat)
 *
 * Descrip:    This function calculates an alignment for AlignWise structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_AlignWise 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_AlignWise
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_AlignWise 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [AlignWise *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_calculate_Small_AlignWise(AlignWise * mat);
#define AlnRangeSet_calculate_Small_AlignWise Wise2_AlnRangeSet_calculate_Small_AlignWise


/* Function:  AlnRangeSet_from_AlignWise(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for AlignWise structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_AlignWise 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_AlignWise
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [AlignWise *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_from_AlignWise(AlignWise * mat);
#define AlnRangeSet_from_AlignWise Wise2_AlnRangeSet_from_AlignWise


/* Function:  convert_PackAln_to_AlnBlock_AlignWise(pal)
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
AlnBlock * Wise2_convert_PackAln_to_AlnBlock_AlignWise(PackAln * pal);
#define convert_PackAln_to_AlnBlock_AlignWise Wise2_convert_PackAln_to_AlnBlock_AlignWise


/* Function:  PackAln_read_Expl_AlignWise(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [AlignWise *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Expl_AlignWise(AlignWise * mat);
#define PackAln_read_Expl_AlignWise Wise2_PackAln_read_Expl_AlignWise


/* Function:  PackAln_read_generic_AlignWise(mat,h)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [AlignWise *]
 * Arg:          h [UNKN ] Undocumented argument [AlignWise_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_generic_AlignWise(AlignWise * mat,AlignWise_access_func_holder h);
#define PackAln_read_generic_AlignWise Wise2_PackAln_read_generic_AlignWise


/* Function:  calculate_AlignWise(mat)
 *
 * Descrip:    This function calculates the AlignWise matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_AlignWise
 *
 *
 * Arg:        mat [UNKN ] AlignWise which contains explicit basematrix memory [AlignWise *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_AlignWise(AlignWise * mat);
#define calculate_AlignWise Wise2_calculate_AlignWise


/* Function:  calculate_dpenv_AlignWise(mat,dpenv)
 *
 * Descrip:    This function calculates the AlignWise matrix when in explicit mode, subject to the envelope
 *
 *
 * Arg:          mat [UNKN ] AlignWise which contains explicit basematrix memory [AlignWise *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_dpenv_AlignWise(AlignWise * mat,DPEnvelope * dpenv);
#define calculate_dpenv_AlignWise Wise2_calculate_dpenv_AlignWise


/* Function:  AlignWise_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AlignWise *]
 *
 */
AlignWise * Wise2_AlignWise_alloc(void);
#define AlignWise_alloc Wise2_AlignWise_alloc


/* Function:  free_AlignWise(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AlignWise *]
 *
 * Return [UNKN ]  Undocumented return value [AlignWise *]
 *
 */
AlignWise * Wise2_free_AlignWise(AlignWise * obj);
#define free_AlignWise Wise2_free_AlignWise


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
int Wise2_AlignWise_shatter_access_main(AlignWise * mat,int i,int j,int state);
#define AlignWise_shatter_access_main Wise2_AlignWise_shatter_access_main
int Wise2_AlignWise_shatter_access_special(AlignWise * mat,int i,int j,int state);
#define AlignWise_shatter_access_special Wise2_AlignWise_shatter_access_special
void * Wise2_thread_loop_AlignWise(void * ptr);
#define thread_loop_AlignWise Wise2_thread_loop_AlignWise
int Wise2_score_only_AlignWise(AlignGeneModelFrame* model,AlignGeneModelScore* align ,Score intronopen,Score geneopen);
#define score_only_AlignWise Wise2_score_only_AlignWise
AlignWise * Wise2_allocate_AlignWise_only(AlignGeneModelFrame* model,AlignGeneModelScore* align ,Score intronopen,Score geneopen);
#define allocate_AlignWise_only Wise2_allocate_AlignWise_only
void Wise2_init_AlignWise(AlignWise * mat);
#define init_AlignWise Wise2_init_AlignWise
AlnRange * Wise2_AlnRange_build_AlignWise(AlignWise * mat,int stopj,int stopspecstate,int * startj,int * startspecstate);
#define AlnRange_build_AlignWise Wise2_AlnRange_build_AlignWise
boolean Wise2_read_hidden_AlignWise(AlignWise * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out);
#define read_hidden_AlignWise Wise2_read_hidden_AlignWise
int Wise2_max_hidden_AlignWise(AlignWise * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_hidden_AlignWise Wise2_max_hidden_AlignWise
boolean Wise2_read_special_strip_AlignWise(AlignWise * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out);
#define read_special_strip_AlignWise Wise2_read_special_strip_AlignWise
int Wise2_max_special_strip_AlignWise(AlignWise * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_special_strip_AlignWise Wise2_max_special_strip_AlignWise
int Wise2_max_matrix_to_special_AlignWise(AlignWise * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_matrix_to_special_AlignWise Wise2_max_matrix_to_special_AlignWise
void Wise2_calculate_hidden_AlignWise(AlignWise * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv);
#define calculate_hidden_AlignWise Wise2_calculate_hidden_AlignWise
void Wise2_init_hidden_AlignWise(AlignWise * mat,int starti,int startj,int stopi,int stopj);
#define init_hidden_AlignWise Wise2_init_hidden_AlignWise
boolean Wise2_full_dc_AlignWise(AlignWise * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv);
#define full_dc_AlignWise Wise2_full_dc_AlignWise
boolean Wise2_do_dc_single_pass_AlignWise(AlignWise * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done);
#define do_dc_single_pass_AlignWise Wise2_do_dc_single_pass_AlignWise
void Wise2_push_dc_at_merge_AlignWise(AlignWise * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv);
#define push_dc_at_merge_AlignWise Wise2_push_dc_at_merge_AlignWise
void Wise2_follow_on_dc_AlignWise(AlignWise * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define follow_on_dc_AlignWise Wise2_follow_on_dc_AlignWise
void Wise2_run_up_dc_AlignWise(AlignWise * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define run_up_dc_AlignWise Wise2_run_up_dc_AlignWise
void Wise2_init_dc_AlignWise(AlignWise * mat);
#define init_dc_AlignWise Wise2_init_dc_AlignWise
int Wise2_start_end_find_end_AlignWise(AlignWise * mat,int * endj);
#define start_end_find_end_AlignWise Wise2_start_end_find_end_AlignWise
boolean Wise2_dc_optimised_start_end_calc_AlignWise(AlignWise *mat,DPEnvelope * dpenv);
#define dc_optimised_start_end_calc_AlignWise Wise2_dc_optimised_start_end_calc_AlignWise
void Wise2_init_start_end_linear_AlignWise(AlignWise * mat);
#define init_start_end_linear_AlignWise Wise2_init_start_end_linear_AlignWise
AlnConvertSet * Wise2_AlnConvertSet_AlignWise(void);
#define AlnConvertSet_AlignWise Wise2_AlnConvertSet_AlignWise
int Wise2_AlignWise_explicit_access_main(AlignWise * mat,int i,int j,int state);
#define AlignWise_explicit_access_main Wise2_AlignWise_explicit_access_main
int Wise2_AlignWise_explicit_access_special(AlignWise * mat,int i,int j,int state);
#define AlignWise_explicit_access_special Wise2_AlignWise_explicit_access_special
int Wise2_find_end_AlignWise(AlignWise * mat,int * ri,int * rj,int * state,boolean * isspecial,AlignWise_access_func_holder h);
#define find_end_AlignWise Wise2_find_end_AlignWise
void Wise2_AlignWise_debug_show_matrix(AlignWise * mat,int starti,int stopi,int startj,int stopj,FILE * ofp);
#define AlignWise_debug_show_matrix Wise2_AlignWise_debug_show_matrix
int Wise2_max_calc_AlignWise(AlignWise * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,AlignWise_access_func_holder h);
#define max_calc_AlignWise Wise2_max_calc_AlignWise
int Wise2_max_calc_special_AlignWise(AlignWise * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,AlignWise_access_func_holder h);
#define max_calc_special_AlignWise Wise2_max_calc_special_AlignWise

#ifdef _cplusplus
}
#endif

#endif
