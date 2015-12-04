#ifndef DYNAMITEpbaHEADERFILE
#define DYNAMITEpbaHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"

struct Wise2_ProteinBlockAligner {  
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
    CompMat* m;  
    Score bentry;    
    Score bexit;     
    Score bfor_trans;    
    Score b_self_trans;  
    Score b3exit;    
    } ;  
/* ProteinBlockAligner defined */ 
#ifndef DYNAMITE_DEFINED_ProteinBlockAligner
typedef struct Wise2_ProteinBlockAligner Wise2_ProteinBlockAligner;
#define ProteinBlockAligner Wise2_ProteinBlockAligner
#define DYNAMITE_DEFINED_ProteinBlockAligner
#endif


#ifdef PTHREAD
struct thread_pool_holder_ProteinBlockAligner {  
    ComplexSequence* q; /* Query object placeholder */ 
    ProteinDB* querydb; /* Query database object */ 
    boolean query_init;  
    ComplexSequence* t; /* Target object placeholder */ 
    ProteinDB* targetdb;/* Target database object */ 
    boolean target_init; 
    CompMat* m;  
    Score bentry;    
    Score bexit;     
    Score bfor_trans;    
    Score b_self_trans;  
    Score b3exit;    
    pthread_mutex_t input_lock;  
    pthread_mutex_t output_lock;     
    Hscore * out;    
    pthread_t * pool;    
    int number_of_threads;   
    boolean search_has_ended;    
    DBSearchImpl * dbsi;     
    } ;  
#endif /* PTHREAD */
struct Wise2_ProteinBlockAligner_access_func_holder {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int (*access_main)(ProteinBlockAligner*,int,int,int);    
    int (*access_special)(ProteinBlockAligner*,int,int,int); 
    } ;  
/* ProteinBlockAligner_access_func_holder defined */ 
#ifndef DYNAMITE_DEFINED_ProteinBlockAligner_access_func_holder
typedef struct Wise2_ProteinBlockAligner_access_func_holder Wise2_ProteinBlockAligner_access_func_holder;
#define ProteinBlockAligner_access_func_holder Wise2_ProteinBlockAligner_access_func_holder
#define DYNAMITE_DEFINED_ProteinBlockAligner_access_func_holder
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  PackAln_read_Shatter_ProteinBlockAligner(mat)
 *
 * Descrip:    Reads off PackAln from shatter matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [ProteinBlockAligner *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Shatter_ProteinBlockAligner(ProteinBlockAligner * mat);
#define PackAln_read_Shatter_ProteinBlockAligner Wise2_PackAln_read_Shatter_ProteinBlockAligner


/* Function:  calculate_shatter_ProteinBlockAligner(mat,dpenv)
 *
 * Descrip:    This function calculates the ProteinBlockAligner matrix when in shatter mode
 *
 *
 * Arg:          mat [UNKN ] (null) [ProteinBlockAligner *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_shatter_ProteinBlockAligner(ProteinBlockAligner * mat,DPEnvelope * dpenv);
#define calculate_shatter_ProteinBlockAligner Wise2_calculate_shatter_ProteinBlockAligner


/* Function:  search_ProteinBlockAligner(dbsi,out,querydb,targetdb,m,bentry,bexit,bfor_trans,b_self_trans,b3exit)
 *
 * Descrip:    This function makes a database search of ProteinBlockAligner
 *             It uses the dbsi structure to choose which implementation
 *             to use of the database searching. This way at run time you
 *             can switch between single threaded/multi-threaded or hardware
 *
 *
 * Arg:                dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 * Arg:                 out [UNKN ] Undocumented argument [Hscore *]
 * Arg:             querydb [UNKN ] Undocumented argument [ProteinDB*]
 * Arg:            targetdb [UNKN ] Undocumented argument [ProteinDB*]
 * Arg:                   m [UNKN ] Undocumented argument [CompMat*]
 * Arg:              bentry [UNKN ] Undocumented argument [Score]
 * Arg:               bexit [UNKN ] Undocumented argument [Score]
 * Arg:          bfor_trans [UNKN ] Undocumented argument [Score]
 * Arg:        b_self_trans [UNKN ] Undocumented argument [Score]
 * Arg:              b3exit [UNKN ] Undocumented argument [Score]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_search_ProteinBlockAligner(DBSearchImpl * dbsi,Hscore * out,ProteinDB* querydb,ProteinDB* targetdb ,CompMat* m,Score bentry,Score bexit,Score bfor_trans,Score b_self_trans,Score b3exit);
#define search_ProteinBlockAligner Wise2_search_ProteinBlockAligner


/* Function:  serial_search_ProteinBlockAligner(out,querydb,targetdb,m,bentry,bexit,bfor_trans,b_self_trans,b3exit)
 *
 * Descrip:    This function makes a database search of ProteinBlockAligner
 *             It is a single processor implementation
 *
 *
 * Arg:                 out [UNKN ] Undocumented argument [Hscore *]
 * Arg:             querydb [UNKN ] Undocumented argument [ProteinDB*]
 * Arg:            targetdb [UNKN ] Undocumented argument [ProteinDB*]
 * Arg:                   m [UNKN ] Undocumented argument [CompMat*]
 * Arg:              bentry [UNKN ] Undocumented argument [Score]
 * Arg:               bexit [UNKN ] Undocumented argument [Score]
 * Arg:          bfor_trans [UNKN ] Undocumented argument [Score]
 * Arg:        b_self_trans [UNKN ] Undocumented argument [Score]
 * Arg:              b3exit [UNKN ] Undocumented argument [Score]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_serial_search_ProteinBlockAligner(Hscore * out,ProteinDB* querydb,ProteinDB* targetdb ,CompMat* m,Score bentry,Score bexit,Score bfor_trans,Score b_self_trans,Score b3exit);
#define serial_search_ProteinBlockAligner Wise2_serial_search_ProteinBlockAligner


/* Function:  PackAln_bestmemory_ProteinBlockAligner(q,t,m,bentry,bexit,bfor_trans,b_self_trans,b3exit,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_ProteinBlockAligner
 *
 *
 * Arg:                   q [UNKN ] query data structure [ComplexSequence*]
 * Arg:                   t [UNKN ] target data structure [ComplexSequence*]
 * Arg:                   m [UNKN ] Resource [CompMat*]
 * Arg:              bentry [UNKN ] Resource [Score]
 * Arg:               bexit [UNKN ] Resource [Score]
 * Arg:          bfor_trans [UNKN ] Resource [Score]
 * Arg:        b_self_trans [UNKN ] Resource [Score]
 * Arg:              b3exit [UNKN ] Resource [Score]
 * Arg:               dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:                dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_bestmemory_ProteinBlockAligner(ComplexSequence* q,ComplexSequence* t ,CompMat* m,Score bentry,Score bexit,Score bfor_trans,Score b_self_trans,Score b3exit,DPEnvelope * dpenv,DPRunImpl * dpri);
#define PackAln_bestmemory_ProteinBlockAligner Wise2_PackAln_bestmemory_ProteinBlockAligner


/* Function:  allocate_Expl_ProteinBlockAligner(q,t,m,bentry,bexit,bfor_trans,b_self_trans,b3exit,dpri)
 *
 * Descrip:    This function allocates the ProteinBlockAligner structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_ProteinBlockAligner_only
 *
 *
 * Arg:                   q [UNKN ] query data structure [ComplexSequence*]
 * Arg:                   t [UNKN ] target data structure [ComplexSequence*]
 * Arg:                   m [UNKN ] Resource [CompMat*]
 * Arg:              bentry [UNKN ] Resource [Score]
 * Arg:               bexit [UNKN ] Resource [Score]
 * Arg:          bfor_trans [UNKN ] Resource [Score]
 * Arg:        b_self_trans [UNKN ] Resource [Score]
 * Arg:              b3exit [UNKN ] Resource [Score]
 * Arg:                dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [ProteinBlockAligner *]
 *
 */
ProteinBlockAligner * Wise2_allocate_Expl_ProteinBlockAligner(ComplexSequence* q,ComplexSequence* t ,CompMat* m,Score bentry,Score bexit,Score bfor_trans,Score b_self_trans,Score b3exit,DPRunImpl * dpri);
#define allocate_Expl_ProteinBlockAligner Wise2_allocate_Expl_ProteinBlockAligner


/* Function:  recalculate_PackAln_ProteinBlockAligner(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by ProteinBlockAligner
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [ProteinBlockAligner *]
 *
 */
void Wise2_recalculate_PackAln_ProteinBlockAligner(PackAln * pal,ProteinBlockAligner * mat);
#define recalculate_PackAln_ProteinBlockAligner Wise2_recalculate_PackAln_ProteinBlockAligner


/* Function:  allocate_Small_ProteinBlockAligner(q,t,m,bentry,bexit,bfor_trans,b_self_trans,b3exit)
 *
 * Descrip:    This function allocates the ProteinBlockAligner structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_ProteinBlockAligner_only
 *
 *
 * Arg:                   q [UNKN ] query data structure [ComplexSequence*]
 * Arg:                   t [UNKN ] target data structure [ComplexSequence*]
 * Arg:                   m [UNKN ] Resource [CompMat*]
 * Arg:              bentry [UNKN ] Resource [Score]
 * Arg:               bexit [UNKN ] Resource [Score]
 * Arg:          bfor_trans [UNKN ] Resource [Score]
 * Arg:        b_self_trans [UNKN ] Resource [Score]
 * Arg:              b3exit [UNKN ] Resource [Score]
 *
 * Return [UNKN ]  Undocumented return value [ProteinBlockAligner *]
 *
 */
ProteinBlockAligner * Wise2_allocate_Small_ProteinBlockAligner(ComplexSequence* q,ComplexSequence* t ,CompMat* m,Score bentry,Score bexit,Score bfor_trans,Score b_self_trans,Score b3exit);
#define allocate_Small_ProteinBlockAligner Wise2_allocate_Small_ProteinBlockAligner


/* Function:  PackAln_calculate_Small_ProteinBlockAligner(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for ProteinBlockAligner structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_ProteinBlockAligner 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_ProteinBlockAligner 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [ProteinBlockAligner *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_calculate_Small_ProteinBlockAligner(ProteinBlockAligner * mat,DPEnvelope * dpenv);
#define PackAln_calculate_Small_ProteinBlockAligner Wise2_PackAln_calculate_Small_ProteinBlockAligner


/* Function:  AlnRangeSet_calculate_Small_ProteinBlockAligner(mat)
 *
 * Descrip:    This function calculates an alignment for ProteinBlockAligner structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_ProteinBlockAligner 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_ProteinBlockAligner
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_ProteinBlockAligner 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [ProteinBlockAligner *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_calculate_Small_ProteinBlockAligner(ProteinBlockAligner * mat);
#define AlnRangeSet_calculate_Small_ProteinBlockAligner Wise2_AlnRangeSet_calculate_Small_ProteinBlockAligner


/* Function:  AlnRangeSet_from_ProteinBlockAligner(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for ProteinBlockAligner structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_ProteinBlockAligner 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_ProteinBlockAligner
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [ProteinBlockAligner *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_from_ProteinBlockAligner(ProteinBlockAligner * mat);
#define AlnRangeSet_from_ProteinBlockAligner Wise2_AlnRangeSet_from_ProteinBlockAligner


/* Function:  convert_PackAln_to_AlnBlock_ProteinBlockAligner(pal)
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
AlnBlock * Wise2_convert_PackAln_to_AlnBlock_ProteinBlockAligner(PackAln * pal);
#define convert_PackAln_to_AlnBlock_ProteinBlockAligner Wise2_convert_PackAln_to_AlnBlock_ProteinBlockAligner


/* Function:  PackAln_read_Expl_ProteinBlockAligner(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [ProteinBlockAligner *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Expl_ProteinBlockAligner(ProteinBlockAligner * mat);
#define PackAln_read_Expl_ProteinBlockAligner Wise2_PackAln_read_Expl_ProteinBlockAligner


/* Function:  PackAln_read_generic_ProteinBlockAligner(mat,h)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [ProteinBlockAligner *]
 * Arg:          h [UNKN ] Undocumented argument [ProteinBlockAligner_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_generic_ProteinBlockAligner(ProteinBlockAligner * mat,ProteinBlockAligner_access_func_holder h);
#define PackAln_read_generic_ProteinBlockAligner Wise2_PackAln_read_generic_ProteinBlockAligner


/* Function:  calculate_ProteinBlockAligner(mat)
 *
 * Descrip:    This function calculates the ProteinBlockAligner matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_ProteinBlockAligner
 *
 *
 * Arg:        mat [UNKN ] ProteinBlockAligner which contains explicit basematrix memory [ProteinBlockAligner *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_ProteinBlockAligner(ProteinBlockAligner * mat);
#define calculate_ProteinBlockAligner Wise2_calculate_ProteinBlockAligner


/* Function:  calculate_dpenv_ProteinBlockAligner(mat,dpenv)
 *
 * Descrip:    This function calculates the ProteinBlockAligner matrix when in explicit mode, subject to the envelope
 *
 *
 * Arg:          mat [UNKN ] ProteinBlockAligner which contains explicit basematrix memory [ProteinBlockAligner *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_dpenv_ProteinBlockAligner(ProteinBlockAligner * mat,DPEnvelope * dpenv);
#define calculate_dpenv_ProteinBlockAligner Wise2_calculate_dpenv_ProteinBlockAligner


/* Function:  ProteinBlockAligner_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ProteinBlockAligner *]
 *
 */
ProteinBlockAligner * Wise2_ProteinBlockAligner_alloc(void);
#define ProteinBlockAligner_alloc Wise2_ProteinBlockAligner_alloc


/* Function:  free_ProteinBlockAligner(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ProteinBlockAligner *]
 *
 * Return [UNKN ]  Undocumented return value [ProteinBlockAligner *]
 *
 */
ProteinBlockAligner * Wise2_free_ProteinBlockAligner(ProteinBlockAligner * obj);
#define free_ProteinBlockAligner Wise2_free_ProteinBlockAligner


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
int Wise2_ProteinBlockAligner_shatter_access_main(ProteinBlockAligner * mat,int i,int j,int state);
#define ProteinBlockAligner_shatter_access_main Wise2_ProteinBlockAligner_shatter_access_main
int Wise2_ProteinBlockAligner_shatter_access_special(ProteinBlockAligner * mat,int i,int j,int state);
#define ProteinBlockAligner_shatter_access_special Wise2_ProteinBlockAligner_shatter_access_special
void * Wise2_thread_loop_ProteinBlockAligner(void * ptr);
#define thread_loop_ProteinBlockAligner Wise2_thread_loop_ProteinBlockAligner
int Wise2_score_only_ProteinBlockAligner(ComplexSequence* q,ComplexSequence* t ,CompMat* m,Score bentry,Score bexit,Score bfor_trans,Score b_self_trans,Score b3exit);
#define score_only_ProteinBlockAligner Wise2_score_only_ProteinBlockAligner
ProteinBlockAligner * Wise2_allocate_ProteinBlockAligner_only(ComplexSequence* q,ComplexSequence* t ,CompMat* m,Score bentry,Score bexit,Score bfor_trans,Score b_self_trans,Score b3exit);
#define allocate_ProteinBlockAligner_only Wise2_allocate_ProteinBlockAligner_only
void Wise2_init_ProteinBlockAligner(ProteinBlockAligner * mat);
#define init_ProteinBlockAligner Wise2_init_ProteinBlockAligner
AlnRange * Wise2_AlnRange_build_ProteinBlockAligner(ProteinBlockAligner * mat,int stopj,int stopspecstate,int * startj,int * startspecstate);
#define AlnRange_build_ProteinBlockAligner Wise2_AlnRange_build_ProteinBlockAligner
boolean Wise2_read_hidden_ProteinBlockAligner(ProteinBlockAligner * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out);
#define read_hidden_ProteinBlockAligner Wise2_read_hidden_ProteinBlockAligner
int Wise2_max_hidden_ProteinBlockAligner(ProteinBlockAligner * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_hidden_ProteinBlockAligner Wise2_max_hidden_ProteinBlockAligner
boolean Wise2_read_special_strip_ProteinBlockAligner(ProteinBlockAligner * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out);
#define read_special_strip_ProteinBlockAligner Wise2_read_special_strip_ProteinBlockAligner
int Wise2_max_special_strip_ProteinBlockAligner(ProteinBlockAligner * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_special_strip_ProteinBlockAligner Wise2_max_special_strip_ProteinBlockAligner
int Wise2_max_matrix_to_special_ProteinBlockAligner(ProteinBlockAligner * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_matrix_to_special_ProteinBlockAligner Wise2_max_matrix_to_special_ProteinBlockAligner
void Wise2_calculate_hidden_ProteinBlockAligner(ProteinBlockAligner * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv);
#define calculate_hidden_ProteinBlockAligner Wise2_calculate_hidden_ProteinBlockAligner
void Wise2_init_hidden_ProteinBlockAligner(ProteinBlockAligner * mat,int starti,int startj,int stopi,int stopj);
#define init_hidden_ProteinBlockAligner Wise2_init_hidden_ProteinBlockAligner
boolean Wise2_full_dc_ProteinBlockAligner(ProteinBlockAligner * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv);
#define full_dc_ProteinBlockAligner Wise2_full_dc_ProteinBlockAligner
boolean Wise2_do_dc_single_pass_ProteinBlockAligner(ProteinBlockAligner * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done);
#define do_dc_single_pass_ProteinBlockAligner Wise2_do_dc_single_pass_ProteinBlockAligner
void Wise2_push_dc_at_merge_ProteinBlockAligner(ProteinBlockAligner * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv);
#define push_dc_at_merge_ProteinBlockAligner Wise2_push_dc_at_merge_ProteinBlockAligner
void Wise2_follow_on_dc_ProteinBlockAligner(ProteinBlockAligner * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define follow_on_dc_ProteinBlockAligner Wise2_follow_on_dc_ProteinBlockAligner
void Wise2_run_up_dc_ProteinBlockAligner(ProteinBlockAligner * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define run_up_dc_ProteinBlockAligner Wise2_run_up_dc_ProteinBlockAligner
void Wise2_init_dc_ProteinBlockAligner(ProteinBlockAligner * mat);
#define init_dc_ProteinBlockAligner Wise2_init_dc_ProteinBlockAligner
int Wise2_start_end_find_end_ProteinBlockAligner(ProteinBlockAligner * mat,int * endj);
#define start_end_find_end_ProteinBlockAligner Wise2_start_end_find_end_ProteinBlockAligner
boolean Wise2_dc_optimised_start_end_calc_ProteinBlockAligner(ProteinBlockAligner *mat,DPEnvelope * dpenv);
#define dc_optimised_start_end_calc_ProteinBlockAligner Wise2_dc_optimised_start_end_calc_ProteinBlockAligner
void Wise2_init_start_end_linear_ProteinBlockAligner(ProteinBlockAligner * mat);
#define init_start_end_linear_ProteinBlockAligner Wise2_init_start_end_linear_ProteinBlockAligner
AlnConvertSet * Wise2_AlnConvertSet_ProteinBlockAligner(void);
#define AlnConvertSet_ProteinBlockAligner Wise2_AlnConvertSet_ProteinBlockAligner
int Wise2_ProteinBlockAligner_explicit_access_main(ProteinBlockAligner * mat,int i,int j,int state);
#define ProteinBlockAligner_explicit_access_main Wise2_ProteinBlockAligner_explicit_access_main
int Wise2_ProteinBlockAligner_explicit_access_special(ProteinBlockAligner * mat,int i,int j,int state);
#define ProteinBlockAligner_explicit_access_special Wise2_ProteinBlockAligner_explicit_access_special
int Wise2_find_end_ProteinBlockAligner(ProteinBlockAligner * mat,int * ri,int * rj,int * state,boolean * isspecial,ProteinBlockAligner_access_func_holder h);
#define find_end_ProteinBlockAligner Wise2_find_end_ProteinBlockAligner
void Wise2_ProteinBlockAligner_debug_show_matrix(ProteinBlockAligner * mat,int starti,int stopi,int startj,int stopj,FILE * ofp);
#define ProteinBlockAligner_debug_show_matrix Wise2_ProteinBlockAligner_debug_show_matrix
int Wise2_max_calc_ProteinBlockAligner(ProteinBlockAligner * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,ProteinBlockAligner_access_func_holder h);
#define max_calc_ProteinBlockAligner Wise2_max_calc_ProteinBlockAligner
int Wise2_max_calc_special_ProteinBlockAligner(ProteinBlockAligner * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,ProteinBlockAligner_access_func_holder h);
#define max_calc_special_ProteinBlockAligner Wise2_max_calc_special_ProteinBlockAligner

#ifdef _cplusplus
}
#endif

#endif
