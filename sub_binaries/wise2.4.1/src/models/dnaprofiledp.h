#ifndef DYNAMITEdnaprofiledpHEADERFILE
#define DYNAMITEdnaprofiledpHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dnaprofile.h"


struct Wise2_DnaProfileMatchScore {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Score ** score;  
    int leni;   /* leni for above score  */ 
    int maxleni;/* max length for above pointer set */ 
    int lenj;   /* lenj for above score  */ 
    int maxlenj;/* max length for above pointer set */ 
    } ;  
/* DnaProfileMatchScore defined */ 
#ifndef DYNAMITE_DEFINED_DnaProfileMatchScore
typedef struct Wise2_DnaProfileMatchScore Wise2_DnaProfileMatchScore;
#define DnaProfileMatchScore Wise2_DnaProfileMatchScore
#define DYNAMITE_DEFINED_DnaProfileMatchScore
#endif


struct Wise2_DnaProfileMat {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    BaseMatrix * basematrix;     
    ShatterMatrix * shatter;     
    int leni;    
    int lenj;    
    DnaProfileScore * q;     
    DnaProfileScore * t;     
    DnaProfileMatchScore* m;     
    Score open_unmatched;    
    Score ext_unmatched;     
    Score gap_unmatched;     
    } ;  
/* DnaProfileMat defined */ 
#ifndef DYNAMITE_DEFINED_DnaProfileMat
typedef struct Wise2_DnaProfileMat Wise2_DnaProfileMat;
#define DnaProfileMat Wise2_DnaProfileMat
#define DYNAMITE_DEFINED_DnaProfileMat
#endif


#ifdef PTHREAD
struct thread_pool_holder_DnaProfileMat {  
    DnaProfileScore * q;/* Static query data: never free'd */ 
    DnaProfileScore * t;/* Static target data: never free'd */ 
    DnaProfileMatchScore* m;     
    Score open_unmatched;    
    Score ext_unmatched;     
    Score gap_unmatched;     
    pthread_mutex_t input_lock;  
    pthread_mutex_t output_lock;     
    Hscore * out;    
    pthread_t * pool;    
    int number_of_threads;   
    boolean search_has_ended;    
    DBSearchImpl * dbsi;     
    } ;  
#endif /* PTHREAD */
struct Wise2_DnaProfileMat_access_func_holder {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int (*access_main)(DnaProfileMat*,int,int,int);  
    int (*access_special)(DnaProfileMat*,int,int,int);   
    } ;  
/* DnaProfileMat_access_func_holder defined */ 
#ifndef DYNAMITE_DEFINED_DnaProfileMat_access_func_holder
typedef struct Wise2_DnaProfileMat_access_func_holder Wise2_DnaProfileMat_access_func_holder;
#define DnaProfileMat_access_func_holder Wise2_DnaProfileMat_access_func_holder
#define DYNAMITE_DEFINED_DnaProfileMat_access_func_holder
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  DnaProfileMatchScore_alloc_matrix(leni,lenj)
 *
 * Descrip:    Allocates structure and matrix
 *
 *
 * Arg:        leni [UNKN ] Length of first dimension of matrix [int]
 * Arg:        lenj [UNKN ] Length of second dimension of matrix [int]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileMatchScore *]
 *
 */
DnaProfileMatchScore * Wise2_DnaProfileMatchScore_alloc_matrix(int leni,int lenj);
#define DnaProfileMatchScore_alloc_matrix Wise2_DnaProfileMatchScore_alloc_matrix


/* Function:  hard_link_DnaProfileMatchScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DnaProfileMatchScore *]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileMatchScore *]
 *
 */
DnaProfileMatchScore * Wise2_hard_link_DnaProfileMatchScore(DnaProfileMatchScore * obj);
#define hard_link_DnaProfileMatchScore Wise2_hard_link_DnaProfileMatchScore


/* Function:  DnaProfileMatchScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileMatchScore *]
 *
 */
DnaProfileMatchScore * Wise2_DnaProfileMatchScore_alloc(void);
#define DnaProfileMatchScore_alloc Wise2_DnaProfileMatchScore_alloc


/* Function:  free_DnaProfileMatchScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DnaProfileMatchScore *]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileMatchScore *]
 *
 */
DnaProfileMatchScore * Wise2_free_DnaProfileMatchScore(DnaProfileMatchScore * obj);
#define free_DnaProfileMatchScore Wise2_free_DnaProfileMatchScore


/* Function:  PackAln_read_Shatter_DnaProfileMat(mat)
 *
 * Descrip:    Reads off PackAln from shatter matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [DnaProfileMat *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Shatter_DnaProfileMat(DnaProfileMat * mat);
#define PackAln_read_Shatter_DnaProfileMat Wise2_PackAln_read_Shatter_DnaProfileMat


/* Function:  calculate_shatter_DnaProfileMat(mat,dpenv)
 *
 * Descrip:    This function calculates the DnaProfileMat matrix when in shatter mode
 *
 *
 * Arg:          mat [UNKN ] (null) [DnaProfileMat *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_shatter_DnaProfileMat(DnaProfileMat * mat,DPEnvelope * dpenv);
#define calculate_shatter_DnaProfileMat Wise2_calculate_shatter_DnaProfileMat


/* Function:  search_DnaProfileMat(dbsi,out,q,t,m,open_unmatched,ext_unmatched,gap_unmatched)
 *
 * Descrip:    This function makes a database search of DnaProfileMat
 *             It uses the dbsi structure to choose which implementation
 *             to use of the database searching. This way at run time you
 *             can switch between single threaded/multi-threaded or hardware
 *
 *
 * Arg:                  dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 * Arg:                   out [UNKN ] Undocumented argument [Hscore *]
 * Arg:                     q [UNKN ] Undocumented argument [DnaProfileScore *]
 * Arg:                     t [UNKN ] Undocumented argument [DnaProfileScore *]
 * Arg:                     m [UNKN ] Undocumented argument [DnaProfileMatchScore*]
 * Arg:        open_unmatched [UNKN ] Undocumented argument [Score]
 * Arg:         ext_unmatched [UNKN ] Undocumented argument [Score]
 * Arg:         gap_unmatched [UNKN ] Undocumented argument [Score]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_search_DnaProfileMat(DBSearchImpl * dbsi,Hscore * out,DnaProfileScore * q,DnaProfileScore * t ,DnaProfileMatchScore* m,Score open_unmatched,Score ext_unmatched,Score gap_unmatched);
#define search_DnaProfileMat Wise2_search_DnaProfileMat


/* Function:  serial_search_DnaProfileMat(out,q,t,m,open_unmatched,ext_unmatched,gap_unmatched)
 *
 * Descrip:    This function makes a database search of DnaProfileMat
 *             It is a single processor implementation
 *
 *
 * Arg:                   out [UNKN ] Undocumented argument [Hscore *]
 * Arg:                     q [UNKN ] Undocumented argument [DnaProfileScore *]
 * Arg:                     t [UNKN ] Undocumented argument [DnaProfileScore *]
 * Arg:                     m [UNKN ] Undocumented argument [DnaProfileMatchScore*]
 * Arg:        open_unmatched [UNKN ] Undocumented argument [Score]
 * Arg:         ext_unmatched [UNKN ] Undocumented argument [Score]
 * Arg:         gap_unmatched [UNKN ] Undocumented argument [Score]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_serial_search_DnaProfileMat(Hscore * out,DnaProfileScore * q,DnaProfileScore * t ,DnaProfileMatchScore* m,Score open_unmatched,Score ext_unmatched,Score gap_unmatched);
#define serial_search_DnaProfileMat Wise2_serial_search_DnaProfileMat


/* Function:  PackAln_bestmemory_DnaProfileMat(q,t,m,open_unmatched,ext_unmatched,gap_unmatched,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_DnaProfileMat
 *
 *
 * Arg:                     q [UNKN ] query data structure [DnaProfileScore *]
 * Arg:                     t [UNKN ] target data structure [DnaProfileScore *]
 * Arg:                     m [UNKN ] Resource [DnaProfileMatchScore*]
 * Arg:        open_unmatched [UNKN ] Resource [Score]
 * Arg:         ext_unmatched [UNKN ] Resource [Score]
 * Arg:         gap_unmatched [UNKN ] Resource [Score]
 * Arg:                 dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:                  dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_bestmemory_DnaProfileMat(DnaProfileScore * q,DnaProfileScore * t ,DnaProfileMatchScore* m,Score open_unmatched,Score ext_unmatched,Score gap_unmatched,DPEnvelope * dpenv,DPRunImpl * dpri);
#define PackAln_bestmemory_DnaProfileMat Wise2_PackAln_bestmemory_DnaProfileMat


/* Function:  allocate_Expl_DnaProfileMat(q,t,m,open_unmatched,ext_unmatched,gap_unmatched,dpri)
 *
 * Descrip:    This function allocates the DnaProfileMat structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_DnaProfileMat_only
 *
 *
 * Arg:                     q [UNKN ] query data structure [DnaProfileScore *]
 * Arg:                     t [UNKN ] target data structure [DnaProfileScore *]
 * Arg:                     m [UNKN ] Resource [DnaProfileMatchScore*]
 * Arg:        open_unmatched [UNKN ] Resource [Score]
 * Arg:         ext_unmatched [UNKN ] Resource [Score]
 * Arg:         gap_unmatched [UNKN ] Resource [Score]
 * Arg:                  dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileMat *]
 *
 */
DnaProfileMat * Wise2_allocate_Expl_DnaProfileMat(DnaProfileScore * q,DnaProfileScore * t ,DnaProfileMatchScore* m,Score open_unmatched,Score ext_unmatched,Score gap_unmatched,DPRunImpl * dpri);
#define allocate_Expl_DnaProfileMat Wise2_allocate_Expl_DnaProfileMat


/* Function:  recalculate_PackAln_DnaProfileMat(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by DnaProfileMat
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [DnaProfileMat *]
 *
 */
void Wise2_recalculate_PackAln_DnaProfileMat(PackAln * pal,DnaProfileMat * mat);
#define recalculate_PackAln_DnaProfileMat Wise2_recalculate_PackAln_DnaProfileMat


/* Function:  allocate_Small_DnaProfileMat(q,t,m,open_unmatched,ext_unmatched,gap_unmatched)
 *
 * Descrip:    This function allocates the DnaProfileMat structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_DnaProfileMat_only
 *
 *
 * Arg:                     q [UNKN ] query data structure [DnaProfileScore *]
 * Arg:                     t [UNKN ] target data structure [DnaProfileScore *]
 * Arg:                     m [UNKN ] Resource [DnaProfileMatchScore*]
 * Arg:        open_unmatched [UNKN ] Resource [Score]
 * Arg:         ext_unmatched [UNKN ] Resource [Score]
 * Arg:         gap_unmatched [UNKN ] Resource [Score]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileMat *]
 *
 */
DnaProfileMat * Wise2_allocate_Small_DnaProfileMat(DnaProfileScore * q,DnaProfileScore * t ,DnaProfileMatchScore* m,Score open_unmatched,Score ext_unmatched,Score gap_unmatched);
#define allocate_Small_DnaProfileMat Wise2_allocate_Small_DnaProfileMat


/* Function:  PackAln_calculate_Small_DnaProfileMat(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for DnaProfileMat structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_DnaProfileMat 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_DnaProfileMat 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [DnaProfileMat *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_calculate_Small_DnaProfileMat(DnaProfileMat * mat,DPEnvelope * dpenv);
#define PackAln_calculate_Small_DnaProfileMat Wise2_PackAln_calculate_Small_DnaProfileMat


/* Function:  AlnRangeSet_calculate_Small_DnaProfileMat(mat)
 *
 * Descrip:    This function calculates an alignment for DnaProfileMat structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_DnaProfileMat 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_DnaProfileMat
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_DnaProfileMat 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [DnaProfileMat *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_calculate_Small_DnaProfileMat(DnaProfileMat * mat);
#define AlnRangeSet_calculate_Small_DnaProfileMat Wise2_AlnRangeSet_calculate_Small_DnaProfileMat


/* Function:  AlnRangeSet_from_DnaProfileMat(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for DnaProfileMat structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_DnaProfileMat 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_DnaProfileMat
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [DnaProfileMat *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_from_DnaProfileMat(DnaProfileMat * mat);
#define AlnRangeSet_from_DnaProfileMat Wise2_AlnRangeSet_from_DnaProfileMat


/* Function:  convert_PackAln_to_AlnBlock_DnaProfileMat(pal)
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
AlnBlock * Wise2_convert_PackAln_to_AlnBlock_DnaProfileMat(PackAln * pal);
#define convert_PackAln_to_AlnBlock_DnaProfileMat Wise2_convert_PackAln_to_AlnBlock_DnaProfileMat


/* Function:  PackAln_read_Expl_DnaProfileMat(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [DnaProfileMat *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Expl_DnaProfileMat(DnaProfileMat * mat);
#define PackAln_read_Expl_DnaProfileMat Wise2_PackAln_read_Expl_DnaProfileMat


/* Function:  PackAln_read_generic_DnaProfileMat(mat,h)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [DnaProfileMat *]
 * Arg:          h [UNKN ] Undocumented argument [DnaProfileMat_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_generic_DnaProfileMat(DnaProfileMat * mat,DnaProfileMat_access_func_holder h);
#define PackAln_read_generic_DnaProfileMat Wise2_PackAln_read_generic_DnaProfileMat


/* Function:  calculate_DnaProfileMat(mat)
 *
 * Descrip:    This function calculates the DnaProfileMat matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_DnaProfileMat
 *
 *
 * Arg:        mat [UNKN ] DnaProfileMat which contains explicit basematrix memory [DnaProfileMat *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_DnaProfileMat(DnaProfileMat * mat);
#define calculate_DnaProfileMat Wise2_calculate_DnaProfileMat


/* Function:  calculate_dpenv_DnaProfileMat(mat,dpenv)
 *
 * Descrip:    This function calculates the DnaProfileMat matrix when in explicit mode, subject to the envelope
 *
 *
 * Arg:          mat [UNKN ] DnaProfileMat which contains explicit basematrix memory [DnaProfileMat *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_dpenv_DnaProfileMat(DnaProfileMat * mat,DPEnvelope * dpenv);
#define calculate_dpenv_DnaProfileMat Wise2_calculate_dpenv_DnaProfileMat


/* Function:  DnaProfileMat_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileMat *]
 *
 */
DnaProfileMat * Wise2_DnaProfileMat_alloc(void);
#define DnaProfileMat_alloc Wise2_DnaProfileMat_alloc


/* Function:  free_DnaProfileMat(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DnaProfileMat *]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileMat *]
 *
 */
DnaProfileMat * Wise2_free_DnaProfileMat(DnaProfileMat * obj);
#define free_DnaProfileMat Wise2_free_DnaProfileMat


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
DnaProfileMatchScore * Wise2_new_ALLR_DnaProfileMatchScore(DnaProfile * q,DnaProfile * t);
#define new_ALLR_DnaProfileMatchScore Wise2_new_ALLR_DnaProfileMatchScore
DnaProfileMatchScore * Wise2_new_DnaProfileMatchScore(DnaProfileScore * q,DnaProfileScore * t);
#define new_DnaProfileMatchScore Wise2_new_DnaProfileMatchScore


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
boolean Wise2_expand_DnaProfileMatchScore(DnaProfileMatchScore * obj,int leni,int lenj);
#define expand_DnaProfileMatchScore Wise2_expand_DnaProfileMatchScore
int Wise2_DnaProfileMat_shatter_access_main(DnaProfileMat * mat,int i,int j,int state);
#define DnaProfileMat_shatter_access_main Wise2_DnaProfileMat_shatter_access_main
int Wise2_DnaProfileMat_shatter_access_special(DnaProfileMat * mat,int i,int j,int state);
#define DnaProfileMat_shatter_access_special Wise2_DnaProfileMat_shatter_access_special
void * Wise2_thread_loop_DnaProfileMat(void * ptr);
#define thread_loop_DnaProfileMat Wise2_thread_loop_DnaProfileMat
int Wise2_score_only_DnaProfileMat(DnaProfileScore * q,DnaProfileScore * t ,DnaProfileMatchScore* m,Score open_unmatched,Score ext_unmatched,Score gap_unmatched);
#define score_only_DnaProfileMat Wise2_score_only_DnaProfileMat
DnaProfileMat * Wise2_allocate_DnaProfileMat_only(DnaProfileScore * q,DnaProfileScore * t ,DnaProfileMatchScore* m,Score open_unmatched,Score ext_unmatched,Score gap_unmatched);
#define allocate_DnaProfileMat_only Wise2_allocate_DnaProfileMat_only
void Wise2_init_DnaProfileMat(DnaProfileMat * mat);
#define init_DnaProfileMat Wise2_init_DnaProfileMat
AlnRange * Wise2_AlnRange_build_DnaProfileMat(DnaProfileMat * mat,int stopj,int stopspecstate,int * startj,int * startspecstate);
#define AlnRange_build_DnaProfileMat Wise2_AlnRange_build_DnaProfileMat
boolean Wise2_read_hidden_DnaProfileMat(DnaProfileMat * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out);
#define read_hidden_DnaProfileMat Wise2_read_hidden_DnaProfileMat
int Wise2_max_hidden_DnaProfileMat(DnaProfileMat * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_hidden_DnaProfileMat Wise2_max_hidden_DnaProfileMat
boolean Wise2_read_special_strip_DnaProfileMat(DnaProfileMat * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out);
#define read_special_strip_DnaProfileMat Wise2_read_special_strip_DnaProfileMat
int Wise2_max_special_strip_DnaProfileMat(DnaProfileMat * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_special_strip_DnaProfileMat Wise2_max_special_strip_DnaProfileMat
int Wise2_max_matrix_to_special_DnaProfileMat(DnaProfileMat * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_matrix_to_special_DnaProfileMat Wise2_max_matrix_to_special_DnaProfileMat
void Wise2_calculate_hidden_DnaProfileMat(DnaProfileMat * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv);
#define calculate_hidden_DnaProfileMat Wise2_calculate_hidden_DnaProfileMat
void Wise2_init_hidden_DnaProfileMat(DnaProfileMat * mat,int starti,int startj,int stopi,int stopj);
#define init_hidden_DnaProfileMat Wise2_init_hidden_DnaProfileMat
boolean Wise2_full_dc_DnaProfileMat(DnaProfileMat * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv);
#define full_dc_DnaProfileMat Wise2_full_dc_DnaProfileMat
boolean Wise2_do_dc_single_pass_DnaProfileMat(DnaProfileMat * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done);
#define do_dc_single_pass_DnaProfileMat Wise2_do_dc_single_pass_DnaProfileMat
void Wise2_push_dc_at_merge_DnaProfileMat(DnaProfileMat * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv);
#define push_dc_at_merge_DnaProfileMat Wise2_push_dc_at_merge_DnaProfileMat
void Wise2_follow_on_dc_DnaProfileMat(DnaProfileMat * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define follow_on_dc_DnaProfileMat Wise2_follow_on_dc_DnaProfileMat
void Wise2_run_up_dc_DnaProfileMat(DnaProfileMat * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define run_up_dc_DnaProfileMat Wise2_run_up_dc_DnaProfileMat
void Wise2_init_dc_DnaProfileMat(DnaProfileMat * mat);
#define init_dc_DnaProfileMat Wise2_init_dc_DnaProfileMat
int Wise2_start_end_find_end_DnaProfileMat(DnaProfileMat * mat,int * endj);
#define start_end_find_end_DnaProfileMat Wise2_start_end_find_end_DnaProfileMat
boolean Wise2_dc_optimised_start_end_calc_DnaProfileMat(DnaProfileMat *mat,DPEnvelope * dpenv);
#define dc_optimised_start_end_calc_DnaProfileMat Wise2_dc_optimised_start_end_calc_DnaProfileMat
void Wise2_init_start_end_linear_DnaProfileMat(DnaProfileMat * mat);
#define init_start_end_linear_DnaProfileMat Wise2_init_start_end_linear_DnaProfileMat
AlnConvertSet * Wise2_AlnConvertSet_DnaProfileMat(void);
#define AlnConvertSet_DnaProfileMat Wise2_AlnConvertSet_DnaProfileMat
int Wise2_DnaProfileMat_explicit_access_main(DnaProfileMat * mat,int i,int j,int state);
#define DnaProfileMat_explicit_access_main Wise2_DnaProfileMat_explicit_access_main
int Wise2_DnaProfileMat_explicit_access_special(DnaProfileMat * mat,int i,int j,int state);
#define DnaProfileMat_explicit_access_special Wise2_DnaProfileMat_explicit_access_special
int Wise2_find_end_DnaProfileMat(DnaProfileMat * mat,int * ri,int * rj,int * state,boolean * isspecial,DnaProfileMat_access_func_holder h);
#define find_end_DnaProfileMat Wise2_find_end_DnaProfileMat
void Wise2_DnaProfileMat_debug_show_matrix(DnaProfileMat * mat,int starti,int stopi,int startj,int stopj,FILE * ofp);
#define DnaProfileMat_debug_show_matrix Wise2_DnaProfileMat_debug_show_matrix
int Wise2_max_calc_DnaProfileMat(DnaProfileMat * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,DnaProfileMat_access_func_holder h);
#define max_calc_DnaProfileMat Wise2_max_calc_DnaProfileMat
int Wise2_max_calc_special_DnaProfileMat(DnaProfileMat * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,DnaProfileMat_access_func_holder h);
#define max_calc_special_DnaProfileMat Wise2_max_calc_special_DnaProfileMat

#ifdef _cplusplus
}
#endif

#endif
