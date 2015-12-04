#ifndef DYNAMITEtransregiondpHEADERFILE
#define DYNAMITEtransregiondpHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "transfactor.h"
#include "dyna.h"

#define TRANS_REGION_MODEL_SCORE(seq,j,model) (seq->coverage[j] == 0 ? model->region_notpresent : model->region_present)
#define TRANS_NON_REGION_MODEL_SCORE(seq,j,model) (seq->coverage[j] == 0 ? model->nonregion_notpresent : model->nonregion_present)
#define TRANS_REGION_GC_SCORE(seq,j,model) (seq->gc_points[j] == 1 ? model->gc_point : model->gc_non_point)

struct Wise2_SequenceBaseCoverage {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * coverage;     
    char * gc_points;    
    int len;     
    Sequence * seq;  
    } ;  
/* SequenceBaseCoverage defined */ 
#ifndef DYNAMITE_DEFINED_SequenceBaseCoverage
typedef struct Wise2_SequenceBaseCoverage Wise2_SequenceBaseCoverage;
#define SequenceBaseCoverage Wise2_SequenceBaseCoverage
#define DYNAMITE_DEFINED_SequenceBaseCoverage
#endif


struct Wise2_TransRegionModel {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int len;     
    Score region_present;    
    Score region_notpresent;     
    Score nonregion_present;     
    Score nonregion_notpresent;  
    Score region_start;  
    Score gc_point;  
    Score gc_non_point;  
    } ;  
/* TransRegionModel defined */ 
#ifndef DYNAMITE_DEFINED_TransRegionModel
typedef struct Wise2_TransRegionModel Wise2_TransRegionModel;
#define TransRegionModel Wise2_TransRegionModel
#define DYNAMITE_DEFINED_TransRegionModel
#endif


struct Wise2_TransRegionMatrix {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    BaseMatrix * basematrix;     
    ShatterMatrix * shatter;     
    int leni;    
    int lenj;    
    TransRegionModel* model;     
    SequenceBaseCoverage* seq;   
    } ;  
/* TransRegionMatrix defined */ 
#ifndef DYNAMITE_DEFINED_TransRegionMatrix
typedef struct Wise2_TransRegionMatrix Wise2_TransRegionMatrix;
#define TransRegionMatrix Wise2_TransRegionMatrix
#define DYNAMITE_DEFINED_TransRegionMatrix
#endif


#ifdef PTHREAD
struct thread_pool_holder_TransRegionMatrix {  
    TransRegionModel* model;/* Static query data: never free'd */ 
    SequenceBaseCoverage* seq;  /* Static target data: never free'd */ 
    pthread_mutex_t input_lock;  
    pthread_mutex_t output_lock;     
    Hscore * out;    
    pthread_t * pool;    
    int number_of_threads;   
    boolean search_has_ended;    
    DBSearchImpl * dbsi;     
    } ;  
#endif /* PTHREAD */
struct Wise2_TransRegionMatrix_access_func_holder {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int (*access_main)(TransRegionMatrix*,int,int,int);  
    int (*access_special)(TransRegionMatrix*,int,int,int);   
    } ;  
/* TransRegionMatrix_access_func_holder defined */ 
#ifndef DYNAMITE_DEFINED_TransRegionMatrix_access_func_holder
typedef struct Wise2_TransRegionMatrix_access_func_holder Wise2_TransRegionMatrix_access_func_holder;
#define TransRegionMatrix_access_func_holder Wise2_TransRegionMatrix_access_func_holder
#define DYNAMITE_DEFINED_TransRegionMatrix_access_func_holder
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_SequenceBaseCoverage(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SequenceBaseCoverage *]
 *
 * Return [UNKN ]  Undocumented return value [SequenceBaseCoverage *]
 *
 */
SequenceBaseCoverage * Wise2_hard_link_SequenceBaseCoverage(SequenceBaseCoverage * obj);
#define hard_link_SequenceBaseCoverage Wise2_hard_link_SequenceBaseCoverage


/* Function:  SequenceBaseCoverage_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SequenceBaseCoverage *]
 *
 */
SequenceBaseCoverage * Wise2_SequenceBaseCoverage_alloc(void);
#define SequenceBaseCoverage_alloc Wise2_SequenceBaseCoverage_alloc


/* Function:  free_SequenceBaseCoverage(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SequenceBaseCoverage *]
 *
 * Return [UNKN ]  Undocumented return value [SequenceBaseCoverage *]
 *
 */
SequenceBaseCoverage * Wise2_free_SequenceBaseCoverage(SequenceBaseCoverage * obj);
#define free_SequenceBaseCoverage Wise2_free_SequenceBaseCoverage


/* Function:  hard_link_TransRegionModel(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [TransRegionModel *]
 *
 * Return [UNKN ]  Undocumented return value [TransRegionModel *]
 *
 */
TransRegionModel * Wise2_hard_link_TransRegionModel(TransRegionModel * obj);
#define hard_link_TransRegionModel Wise2_hard_link_TransRegionModel


/* Function:  TransRegionModel_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TransRegionModel *]
 *
 */
TransRegionModel * Wise2_TransRegionModel_alloc(void);
#define TransRegionModel_alloc Wise2_TransRegionModel_alloc


/* Function:  free_TransRegionModel(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [TransRegionModel *]
 *
 * Return [UNKN ]  Undocumented return value [TransRegionModel *]
 *
 */
TransRegionModel * Wise2_free_TransRegionModel(TransRegionModel * obj);
#define free_TransRegionModel Wise2_free_TransRegionModel


/* Function:  PackAln_read_Shatter_TransRegionMatrix(mat)
 *
 * Descrip:    Reads off PackAln from shatter matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [TransRegionMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Shatter_TransRegionMatrix(TransRegionMatrix * mat);
#define PackAln_read_Shatter_TransRegionMatrix Wise2_PackAln_read_Shatter_TransRegionMatrix


/* Function:  calculate_shatter_TransRegionMatrix(mat,dpenv)
 *
 * Descrip:    This function calculates the TransRegionMatrix matrix when in shatter mode
 *
 *
 * Arg:          mat [UNKN ] (null) [TransRegionMatrix *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_shatter_TransRegionMatrix(TransRegionMatrix * mat,DPEnvelope * dpenv);
#define calculate_shatter_TransRegionMatrix Wise2_calculate_shatter_TransRegionMatrix


/* Function:  search_TransRegionMatrix(dbsi,out,model,seq)
 *
 * Descrip:    This function makes a database search of TransRegionMatrix
 *             It uses the dbsi structure to choose which implementation
 *             to use of the database searching. This way at run time you
 *             can switch between single threaded/multi-threaded or hardware
 *
 *
 * Arg:         dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 * Arg:          out [UNKN ] Undocumented argument [Hscore *]
 * Arg:        model [UNKN ] Undocumented argument [TransRegionModel*]
 * Arg:          seq [UNKN ] Undocumented argument [SequenceBaseCoverage*]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_search_TransRegionMatrix(DBSearchImpl * dbsi,Hscore * out,TransRegionModel* model,SequenceBaseCoverage* seq );
#define search_TransRegionMatrix Wise2_search_TransRegionMatrix


/* Function:  serial_search_TransRegionMatrix(out,model,seq)
 *
 * Descrip:    This function makes a database search of TransRegionMatrix
 *             It is a single processor implementation
 *
 *
 * Arg:          out [UNKN ] Undocumented argument [Hscore *]
 * Arg:        model [UNKN ] Undocumented argument [TransRegionModel*]
 * Arg:          seq [UNKN ] Undocumented argument [SequenceBaseCoverage*]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_serial_search_TransRegionMatrix(Hscore * out,TransRegionModel* model,SequenceBaseCoverage* seq );
#define serial_search_TransRegionMatrix Wise2_serial_search_TransRegionMatrix


/* Function:  PackAln_bestmemory_TransRegionMatrix(model,seq,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_TransRegionMatrix
 *
 *
 * Arg:        model [UNKN ] query data structure [TransRegionModel*]
 * Arg:          seq [UNKN ] target data structure [SequenceBaseCoverage*]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:         dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_bestmemory_TransRegionMatrix(TransRegionModel* model,SequenceBaseCoverage* seq ,DPEnvelope * dpenv,DPRunImpl * dpri);
#define PackAln_bestmemory_TransRegionMatrix Wise2_PackAln_bestmemory_TransRegionMatrix


/* Function:  allocate_Expl_TransRegionMatrix(model,seq,dpri)
 *
 * Descrip:    This function allocates the TransRegionMatrix structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_TransRegionMatrix_only
 *
 *
 * Arg:        model [UNKN ] query data structure [TransRegionModel*]
 * Arg:          seq [UNKN ] target data structure [SequenceBaseCoverage*]
 * Arg:         dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [TransRegionMatrix *]
 *
 */
TransRegionMatrix * Wise2_allocate_Expl_TransRegionMatrix(TransRegionModel* model,SequenceBaseCoverage* seq ,DPRunImpl * dpri);
#define allocate_Expl_TransRegionMatrix Wise2_allocate_Expl_TransRegionMatrix


/* Function:  recalculate_PackAln_TransRegionMatrix(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by TransRegionMatrix
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [TransRegionMatrix *]
 *
 */
void Wise2_recalculate_PackAln_TransRegionMatrix(PackAln * pal,TransRegionMatrix * mat);
#define recalculate_PackAln_TransRegionMatrix Wise2_recalculate_PackAln_TransRegionMatrix


/* Function:  allocate_Small_TransRegionMatrix(model,seq)
 *
 * Descrip:    This function allocates the TransRegionMatrix structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_TransRegionMatrix_only
 *
 *
 * Arg:        model [UNKN ] query data structure [TransRegionModel*]
 * Arg:          seq [UNKN ] target data structure [SequenceBaseCoverage*]
 *
 * Return [UNKN ]  Undocumented return value [TransRegionMatrix *]
 *
 */
TransRegionMatrix * Wise2_allocate_Small_TransRegionMatrix(TransRegionModel* model,SequenceBaseCoverage* seq );
#define allocate_Small_TransRegionMatrix Wise2_allocate_Small_TransRegionMatrix


/* Function:  PackAln_calculate_Small_TransRegionMatrix(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for TransRegionMatrix structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_TransRegionMatrix 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_TransRegionMatrix 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [TransRegionMatrix *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_calculate_Small_TransRegionMatrix(TransRegionMatrix * mat,DPEnvelope * dpenv);
#define PackAln_calculate_Small_TransRegionMatrix Wise2_PackAln_calculate_Small_TransRegionMatrix


/* Function:  AlnRangeSet_calculate_Small_TransRegionMatrix(mat)
 *
 * Descrip:    This function calculates an alignment for TransRegionMatrix structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_TransRegionMatrix 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_TransRegionMatrix
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_TransRegionMatrix 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [TransRegionMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_calculate_Small_TransRegionMatrix(TransRegionMatrix * mat);
#define AlnRangeSet_calculate_Small_TransRegionMatrix Wise2_AlnRangeSet_calculate_Small_TransRegionMatrix


/* Function:  AlnRangeSet_from_TransRegionMatrix(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for TransRegionMatrix structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_TransRegionMatrix 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_TransRegionMatrix
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [TransRegionMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_from_TransRegionMatrix(TransRegionMatrix * mat);
#define AlnRangeSet_from_TransRegionMatrix Wise2_AlnRangeSet_from_TransRegionMatrix


/* Function:  convert_PackAln_to_AlnBlock_TransRegionMatrix(pal)
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
AlnBlock * Wise2_convert_PackAln_to_AlnBlock_TransRegionMatrix(PackAln * pal);
#define convert_PackAln_to_AlnBlock_TransRegionMatrix Wise2_convert_PackAln_to_AlnBlock_TransRegionMatrix


/* Function:  PackAln_read_Expl_TransRegionMatrix(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [TransRegionMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Expl_TransRegionMatrix(TransRegionMatrix * mat);
#define PackAln_read_Expl_TransRegionMatrix Wise2_PackAln_read_Expl_TransRegionMatrix


/* Function:  PackAln_read_generic_TransRegionMatrix(mat,h)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [TransRegionMatrix *]
 * Arg:          h [UNKN ] Undocumented argument [TransRegionMatrix_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_generic_TransRegionMatrix(TransRegionMatrix * mat,TransRegionMatrix_access_func_holder h);
#define PackAln_read_generic_TransRegionMatrix Wise2_PackAln_read_generic_TransRegionMatrix


/* Function:  calculate_TransRegionMatrix(mat)
 *
 * Descrip:    This function calculates the TransRegionMatrix matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_TransRegionMatrix
 *
 *
 * Arg:        mat [UNKN ] TransRegionMatrix which contains explicit basematrix memory [TransRegionMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_TransRegionMatrix(TransRegionMatrix * mat);
#define calculate_TransRegionMatrix Wise2_calculate_TransRegionMatrix


/* Function:  calculate_dpenv_TransRegionMatrix(mat,dpenv)
 *
 * Descrip:    This function calculates the TransRegionMatrix matrix when in explicit mode, subject to the envelope
 *
 *
 * Arg:          mat [UNKN ] TransRegionMatrix which contains explicit basematrix memory [TransRegionMatrix *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_dpenv_TransRegionMatrix(TransRegionMatrix * mat,DPEnvelope * dpenv);
#define calculate_dpenv_TransRegionMatrix Wise2_calculate_dpenv_TransRegionMatrix


/* Function:  TransRegionMatrix_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TransRegionMatrix *]
 *
 */
TransRegionMatrix * Wise2_TransRegionMatrix_alloc(void);
#define TransRegionMatrix_alloc Wise2_TransRegionMatrix_alloc


/* Function:  free_TransRegionMatrix(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [TransRegionMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [TransRegionMatrix *]
 *
 */
TransRegionMatrix * Wise2_free_TransRegionMatrix(TransRegionMatrix * obj);
#define free_TransRegionMatrix Wise2_free_TransRegionMatrix


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
SequenceBaseCoverage * Wise2_new_SequenceBaseCoverage(TransFactorMatchSet * tfms);
#define new_SequenceBaseCoverage Wise2_new_SequenceBaseCoverage
TransRegionModel * Wise2_new_logodds_TransRegionModel(double in_region_prob,double out_region_prob,double in_cost,double gc_region_ratio);
#define new_logodds_TransRegionModel Wise2_new_logodds_TransRegionModel


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
int Wise2_TransRegionMatrix_shatter_access_main(TransRegionMatrix * mat,int i,int j,int state);
#define TransRegionMatrix_shatter_access_main Wise2_TransRegionMatrix_shatter_access_main
int Wise2_TransRegionMatrix_shatter_access_special(TransRegionMatrix * mat,int i,int j,int state);
#define TransRegionMatrix_shatter_access_special Wise2_TransRegionMatrix_shatter_access_special
void * Wise2_thread_loop_TransRegionMatrix(void * ptr);
#define thread_loop_TransRegionMatrix Wise2_thread_loop_TransRegionMatrix
int Wise2_score_only_TransRegionMatrix(TransRegionModel* model,SequenceBaseCoverage* seq );
#define score_only_TransRegionMatrix Wise2_score_only_TransRegionMatrix
TransRegionMatrix * Wise2_allocate_TransRegionMatrix_only(TransRegionModel* model,SequenceBaseCoverage* seq );
#define allocate_TransRegionMatrix_only Wise2_allocate_TransRegionMatrix_only
void Wise2_init_TransRegionMatrix(TransRegionMatrix * mat);
#define init_TransRegionMatrix Wise2_init_TransRegionMatrix
AlnRange * Wise2_AlnRange_build_TransRegionMatrix(TransRegionMatrix * mat,int stopj,int stopspecstate,int * startj,int * startspecstate);
#define AlnRange_build_TransRegionMatrix Wise2_AlnRange_build_TransRegionMatrix
boolean Wise2_read_hidden_TransRegionMatrix(TransRegionMatrix * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out);
#define read_hidden_TransRegionMatrix Wise2_read_hidden_TransRegionMatrix
int Wise2_max_hidden_TransRegionMatrix(TransRegionMatrix * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_hidden_TransRegionMatrix Wise2_max_hidden_TransRegionMatrix
boolean Wise2_read_special_strip_TransRegionMatrix(TransRegionMatrix * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out);
#define read_special_strip_TransRegionMatrix Wise2_read_special_strip_TransRegionMatrix
int Wise2_max_special_strip_TransRegionMatrix(TransRegionMatrix * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_special_strip_TransRegionMatrix Wise2_max_special_strip_TransRegionMatrix
int Wise2_max_matrix_to_special_TransRegionMatrix(TransRegionMatrix * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_matrix_to_special_TransRegionMatrix Wise2_max_matrix_to_special_TransRegionMatrix
void Wise2_calculate_hidden_TransRegionMatrix(TransRegionMatrix * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv);
#define calculate_hidden_TransRegionMatrix Wise2_calculate_hidden_TransRegionMatrix
void Wise2_init_hidden_TransRegionMatrix(TransRegionMatrix * mat,int starti,int startj,int stopi,int stopj);
#define init_hidden_TransRegionMatrix Wise2_init_hidden_TransRegionMatrix
boolean Wise2_full_dc_TransRegionMatrix(TransRegionMatrix * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv);
#define full_dc_TransRegionMatrix Wise2_full_dc_TransRegionMatrix
boolean Wise2_do_dc_single_pass_TransRegionMatrix(TransRegionMatrix * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done);
#define do_dc_single_pass_TransRegionMatrix Wise2_do_dc_single_pass_TransRegionMatrix
void Wise2_push_dc_at_merge_TransRegionMatrix(TransRegionMatrix * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv);
#define push_dc_at_merge_TransRegionMatrix Wise2_push_dc_at_merge_TransRegionMatrix
void Wise2_follow_on_dc_TransRegionMatrix(TransRegionMatrix * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define follow_on_dc_TransRegionMatrix Wise2_follow_on_dc_TransRegionMatrix
void Wise2_run_up_dc_TransRegionMatrix(TransRegionMatrix * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define run_up_dc_TransRegionMatrix Wise2_run_up_dc_TransRegionMatrix
void Wise2_init_dc_TransRegionMatrix(TransRegionMatrix * mat);
#define init_dc_TransRegionMatrix Wise2_init_dc_TransRegionMatrix
int Wise2_start_end_find_end_TransRegionMatrix(TransRegionMatrix * mat,int * endj);
#define start_end_find_end_TransRegionMatrix Wise2_start_end_find_end_TransRegionMatrix
boolean Wise2_dc_optimised_start_end_calc_TransRegionMatrix(TransRegionMatrix *mat,DPEnvelope * dpenv);
#define dc_optimised_start_end_calc_TransRegionMatrix Wise2_dc_optimised_start_end_calc_TransRegionMatrix
void Wise2_init_start_end_linear_TransRegionMatrix(TransRegionMatrix * mat);
#define init_start_end_linear_TransRegionMatrix Wise2_init_start_end_linear_TransRegionMatrix
AlnConvertSet * Wise2_AlnConvertSet_TransRegionMatrix(void);
#define AlnConvertSet_TransRegionMatrix Wise2_AlnConvertSet_TransRegionMatrix
int Wise2_TransRegionMatrix_explicit_access_main(TransRegionMatrix * mat,int i,int j,int state);
#define TransRegionMatrix_explicit_access_main Wise2_TransRegionMatrix_explicit_access_main
int Wise2_TransRegionMatrix_explicit_access_special(TransRegionMatrix * mat,int i,int j,int state);
#define TransRegionMatrix_explicit_access_special Wise2_TransRegionMatrix_explicit_access_special
int Wise2_find_end_TransRegionMatrix(TransRegionMatrix * mat,int * ri,int * rj,int * state,boolean * isspecial,TransRegionMatrix_access_func_holder h);
#define find_end_TransRegionMatrix Wise2_find_end_TransRegionMatrix
void Wise2_TransRegionMatrix_debug_show_matrix(TransRegionMatrix * mat,int starti,int stopi,int startj,int stopj,FILE * ofp);
#define TransRegionMatrix_debug_show_matrix Wise2_TransRegionMatrix_debug_show_matrix
int Wise2_max_calc_TransRegionMatrix(TransRegionMatrix * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,TransRegionMatrix_access_func_holder h);
#define max_calc_TransRegionMatrix Wise2_max_calc_TransRegionMatrix
int Wise2_max_calc_special_TransRegionMatrix(TransRegionMatrix * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,TransRegionMatrix_access_func_holder h);
#define max_calc_special_TransRegionMatrix Wise2_max_calc_special_TransRegionMatrix

#ifdef _cplusplus
}
#endif

#endif
