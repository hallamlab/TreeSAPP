#ifndef DYNAMITEsywise20HEADERFILE
#define DYNAMITEsywise20HEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"
#include "pairbaseseq.h"
#include "pairbase.h"
#include "syexonmodel.h"


struct Wise2_SyWise20 {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    BaseMatrix * basematrix;     
    ShatterMatrix * shatter;     
    int leni;    
    int lenj;    
    SyExonScore* exonmodel;  
    ComplexSequence* seq;    
    PairBaseCodonModelScore* codon;  
    PairBaseModelScore* nonc;    
    PairBaseCodonModelScore* start;  
    PairBaseCodonModelScore* stop;   
    Score intron_open;   
    Score gene_open;     
    Score nonc_cost;     
    } ;  
/* SyWise20 defined */ 
#ifndef DYNAMITE_DEFINED_SyWise20
typedef struct Wise2_SyWise20 Wise2_SyWise20;
#define SyWise20 Wise2_SyWise20
#define DYNAMITE_DEFINED_SyWise20
#endif


#ifdef PTHREAD
struct thread_pool_holder_SyWise20 {  
    SyExonScore* exonmodel; /* Static query data: never free'd */ 
    ComplexSequence* seq;   /* Static target data: never free'd */ 
    PairBaseCodonModelScore* codon;  
    PairBaseModelScore* nonc;    
    PairBaseCodonModelScore* start;  
    PairBaseCodonModelScore* stop;   
    Score intron_open;   
    Score gene_open;     
    Score nonc_cost;     
    pthread_mutex_t input_lock;  
    pthread_mutex_t output_lock;     
    Hscore * out;    
    pthread_t * pool;    
    int number_of_threads;   
    boolean search_has_ended;    
    DBSearchImpl * dbsi;     
    } ;  
#endif /* PTHREAD */
struct Wise2_SyWise20_access_func_holder {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int (*access_main)(SyWise20*,int,int,int);   
    int (*access_special)(SyWise20*,int,int,int);    
    } ;  
/* SyWise20_access_func_holder defined */ 
#ifndef DYNAMITE_DEFINED_SyWise20_access_func_holder
typedef struct Wise2_SyWise20_access_func_holder Wise2_SyWise20_access_func_holder;
#define SyWise20_access_func_holder Wise2_SyWise20_access_func_holder
#define DYNAMITE_DEFINED_SyWise20_access_func_holder
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  PackAln_read_Shatter_SyWise20(mat)
 *
 * Descrip:    Reads off PackAln from shatter matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [SyWise20 *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Shatter_SyWise20(SyWise20 * mat);
#define PackAln_read_Shatter_SyWise20 Wise2_PackAln_read_Shatter_SyWise20


/* Function:  calculate_shatter_SyWise20(mat,dpenv)
 *
 * Descrip:    This function calculates the SyWise20 matrix when in shatter mode
 *
 *
 * Arg:          mat [UNKN ] (null) [SyWise20 *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_shatter_SyWise20(SyWise20 * mat,DPEnvelope * dpenv);
#define calculate_shatter_SyWise20 Wise2_calculate_shatter_SyWise20


/* Function:  search_SyWise20(dbsi,out,exonmodel,seq,codon,nonc,start,stop,intron_open,gene_open,nonc_cost)
 *
 * Descrip:    This function makes a database search of SyWise20
 *             It uses the dbsi structure to choose which implementation
 *             to use of the database searching. This way at run time you
 *             can switch between single threaded/multi-threaded or hardware
 *
 *
 * Arg:               dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 * Arg:                out [UNKN ] Undocumented argument [Hscore *]
 * Arg:          exonmodel [UNKN ] Undocumented argument [SyExonScore*]
 * Arg:                seq [UNKN ] Undocumented argument [ComplexSequence*]
 * Arg:              codon [UNKN ] Undocumented argument [PairBaseCodonModelScore*]
 * Arg:               nonc [UNKN ] Undocumented argument [PairBaseModelScore*]
 * Arg:              start [UNKN ] Undocumented argument [PairBaseCodonModelScore*]
 * Arg:               stop [UNKN ] Undocumented argument [PairBaseCodonModelScore*]
 * Arg:        intron_open [UNKN ] Undocumented argument [Score]
 * Arg:          gene_open [UNKN ] Undocumented argument [Score]
 * Arg:          nonc_cost [UNKN ] Undocumented argument [Score]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_search_SyWise20(DBSearchImpl * dbsi,Hscore * out,SyExonScore* exonmodel,ComplexSequence* seq ,PairBaseCodonModelScore* codon,PairBaseModelScore* nonc,PairBaseCodonModelScore* start,PairBaseCodonModelScore* stop,Score intron_open,Score gene_open,Score nonc_cost);
#define search_SyWise20 Wise2_search_SyWise20


/* Function:  serial_search_SyWise20(out,exonmodel,seq,codon,nonc,start,stop,intron_open,gene_open,nonc_cost)
 *
 * Descrip:    This function makes a database search of SyWise20
 *             It is a single processor implementation
 *
 *
 * Arg:                out [UNKN ] Undocumented argument [Hscore *]
 * Arg:          exonmodel [UNKN ] Undocumented argument [SyExonScore*]
 * Arg:                seq [UNKN ] Undocumented argument [ComplexSequence*]
 * Arg:              codon [UNKN ] Undocumented argument [PairBaseCodonModelScore*]
 * Arg:               nonc [UNKN ] Undocumented argument [PairBaseModelScore*]
 * Arg:              start [UNKN ] Undocumented argument [PairBaseCodonModelScore*]
 * Arg:               stop [UNKN ] Undocumented argument [PairBaseCodonModelScore*]
 * Arg:        intron_open [UNKN ] Undocumented argument [Score]
 * Arg:          gene_open [UNKN ] Undocumented argument [Score]
 * Arg:          nonc_cost [UNKN ] Undocumented argument [Score]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_serial_search_SyWise20(Hscore * out,SyExonScore* exonmodel,ComplexSequence* seq ,PairBaseCodonModelScore* codon,PairBaseModelScore* nonc,PairBaseCodonModelScore* start,PairBaseCodonModelScore* stop,Score intron_open,Score gene_open,Score nonc_cost);
#define serial_search_SyWise20 Wise2_serial_search_SyWise20


/* Function:  PackAln_bestmemory_SyWise20(exonmodel,seq,codon,nonc,start,stop,intron_open,gene_open,nonc_cost,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_SyWise20
 *
 *
 * Arg:          exonmodel [UNKN ] query data structure [SyExonScore*]
 * Arg:                seq [UNKN ] target data structure [ComplexSequence*]
 * Arg:              codon [UNKN ] Resource [PairBaseCodonModelScore*]
 * Arg:               nonc [UNKN ] Resource [PairBaseModelScore*]
 * Arg:              start [UNKN ] Resource [PairBaseCodonModelScore*]
 * Arg:               stop [UNKN ] Resource [PairBaseCodonModelScore*]
 * Arg:        intron_open [UNKN ] Resource [Score]
 * Arg:          gene_open [UNKN ] Resource [Score]
 * Arg:          nonc_cost [UNKN ] Resource [Score]
 * Arg:              dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:               dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_bestmemory_SyWise20(SyExonScore* exonmodel,ComplexSequence* seq ,PairBaseCodonModelScore* codon,PairBaseModelScore* nonc,PairBaseCodonModelScore* start,PairBaseCodonModelScore* stop,Score intron_open,Score gene_open,Score nonc_cost,DPEnvelope * dpenv,DPRunImpl * dpri);
#define PackAln_bestmemory_SyWise20 Wise2_PackAln_bestmemory_SyWise20


/* Function:  allocate_Expl_SyWise20(exonmodel,seq,codon,nonc,start,stop,intron_open,gene_open,nonc_cost,dpri)
 *
 * Descrip:    This function allocates the SyWise20 structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_SyWise20_only
 *
 *
 * Arg:          exonmodel [UNKN ] query data structure [SyExonScore*]
 * Arg:                seq [UNKN ] target data structure [ComplexSequence*]
 * Arg:              codon [UNKN ] Resource [PairBaseCodonModelScore*]
 * Arg:               nonc [UNKN ] Resource [PairBaseModelScore*]
 * Arg:              start [UNKN ] Resource [PairBaseCodonModelScore*]
 * Arg:               stop [UNKN ] Resource [PairBaseCodonModelScore*]
 * Arg:        intron_open [UNKN ] Resource [Score]
 * Arg:          gene_open [UNKN ] Resource [Score]
 * Arg:          nonc_cost [UNKN ] Resource [Score]
 * Arg:               dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [SyWise20 *]
 *
 */
SyWise20 * Wise2_allocate_Expl_SyWise20(SyExonScore* exonmodel,ComplexSequence* seq ,PairBaseCodonModelScore* codon,PairBaseModelScore* nonc,PairBaseCodonModelScore* start,PairBaseCodonModelScore* stop,Score intron_open,Score gene_open,Score nonc_cost,DPRunImpl * dpri);
#define allocate_Expl_SyWise20 Wise2_allocate_Expl_SyWise20


/* Function:  recalculate_PackAln_SyWise20(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by SyWise20
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [SyWise20 *]
 *
 */
void Wise2_recalculate_PackAln_SyWise20(PackAln * pal,SyWise20 * mat);
#define recalculate_PackAln_SyWise20 Wise2_recalculate_PackAln_SyWise20


/* Function:  allocate_Small_SyWise20(exonmodel,seq,codon,nonc,start,stop,intron_open,gene_open,nonc_cost)
 *
 * Descrip:    This function allocates the SyWise20 structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_SyWise20_only
 *
 *
 * Arg:          exonmodel [UNKN ] query data structure [SyExonScore*]
 * Arg:                seq [UNKN ] target data structure [ComplexSequence*]
 * Arg:              codon [UNKN ] Resource [PairBaseCodonModelScore*]
 * Arg:               nonc [UNKN ] Resource [PairBaseModelScore*]
 * Arg:              start [UNKN ] Resource [PairBaseCodonModelScore*]
 * Arg:               stop [UNKN ] Resource [PairBaseCodonModelScore*]
 * Arg:        intron_open [UNKN ] Resource [Score]
 * Arg:          gene_open [UNKN ] Resource [Score]
 * Arg:          nonc_cost [UNKN ] Resource [Score]
 *
 * Return [UNKN ]  Undocumented return value [SyWise20 *]
 *
 */
SyWise20 * Wise2_allocate_Small_SyWise20(SyExonScore* exonmodel,ComplexSequence* seq ,PairBaseCodonModelScore* codon,PairBaseModelScore* nonc,PairBaseCodonModelScore* start,PairBaseCodonModelScore* stop,Score intron_open,Score gene_open,Score nonc_cost);
#define allocate_Small_SyWise20 Wise2_allocate_Small_SyWise20


/* Function:  PackAln_calculate_Small_SyWise20(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for SyWise20 structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_SyWise20 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_SyWise20 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [SyWise20 *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_calculate_Small_SyWise20(SyWise20 * mat,DPEnvelope * dpenv);
#define PackAln_calculate_Small_SyWise20 Wise2_PackAln_calculate_Small_SyWise20


/* Function:  AlnRangeSet_calculate_Small_SyWise20(mat)
 *
 * Descrip:    This function calculates an alignment for SyWise20 structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_SyWise20 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_SyWise20
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_SyWise20 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [SyWise20 *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_calculate_Small_SyWise20(SyWise20 * mat);
#define AlnRangeSet_calculate_Small_SyWise20 Wise2_AlnRangeSet_calculate_Small_SyWise20


/* Function:  AlnRangeSet_from_SyWise20(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for SyWise20 structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_SyWise20 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_SyWise20
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [SyWise20 *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_from_SyWise20(SyWise20 * mat);
#define AlnRangeSet_from_SyWise20 Wise2_AlnRangeSet_from_SyWise20


/* Function:  convert_PackAln_to_AlnBlock_SyWise20(pal)
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
AlnBlock * Wise2_convert_PackAln_to_AlnBlock_SyWise20(PackAln * pal);
#define convert_PackAln_to_AlnBlock_SyWise20 Wise2_convert_PackAln_to_AlnBlock_SyWise20


/* Function:  PackAln_read_Expl_SyWise20(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [SyWise20 *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Expl_SyWise20(SyWise20 * mat);
#define PackAln_read_Expl_SyWise20 Wise2_PackAln_read_Expl_SyWise20


/* Function:  PackAln_read_generic_SyWise20(mat,h)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [SyWise20 *]
 * Arg:          h [UNKN ] Undocumented argument [SyWise20_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_generic_SyWise20(SyWise20 * mat,SyWise20_access_func_holder h);
#define PackAln_read_generic_SyWise20 Wise2_PackAln_read_generic_SyWise20


/* Function:  calculate_SyWise20(mat)
 *
 * Descrip:    This function calculates the SyWise20 matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_SyWise20
 *
 *
 * Arg:        mat [UNKN ] SyWise20 which contains explicit basematrix memory [SyWise20 *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_SyWise20(SyWise20 * mat);
#define calculate_SyWise20 Wise2_calculate_SyWise20


/* Function:  calculate_dpenv_SyWise20(mat,dpenv)
 *
 * Descrip:    This function calculates the SyWise20 matrix when in explicit mode, subject to the envelope
 *
 *
 * Arg:          mat [UNKN ] SyWise20 which contains explicit basematrix memory [SyWise20 *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_dpenv_SyWise20(SyWise20 * mat,DPEnvelope * dpenv);
#define calculate_dpenv_SyWise20 Wise2_calculate_dpenv_SyWise20


/* Function:  SyWise20_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SyWise20 *]
 *
 */
SyWise20 * Wise2_SyWise20_alloc(void);
#define SyWise20_alloc Wise2_SyWise20_alloc


/* Function:  free_SyWise20(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SyWise20 *]
 *
 * Return [UNKN ]  Undocumented return value [SyWise20 *]
 *
 */
SyWise20 * Wise2_free_SyWise20(SyWise20 * obj);
#define free_SyWise20 Wise2_free_SyWise20


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
int Wise2_SyWise20_shatter_access_main(SyWise20 * mat,int i,int j,int state);
#define SyWise20_shatter_access_main Wise2_SyWise20_shatter_access_main
int Wise2_SyWise20_shatter_access_special(SyWise20 * mat,int i,int j,int state);
#define SyWise20_shatter_access_special Wise2_SyWise20_shatter_access_special
void * Wise2_thread_loop_SyWise20(void * ptr);
#define thread_loop_SyWise20 Wise2_thread_loop_SyWise20
int Wise2_score_only_SyWise20(SyExonScore* exonmodel,ComplexSequence* seq ,PairBaseCodonModelScore* codon,PairBaseModelScore* nonc,PairBaseCodonModelScore* start,PairBaseCodonModelScore* stop,Score intron_open,Score gene_open,Score nonc_cost);
#define score_only_SyWise20 Wise2_score_only_SyWise20
SyWise20 * Wise2_allocate_SyWise20_only(SyExonScore* exonmodel,ComplexSequence* seq ,PairBaseCodonModelScore* codon,PairBaseModelScore* nonc,PairBaseCodonModelScore* start,PairBaseCodonModelScore* stop,Score intron_open,Score gene_open,Score nonc_cost);
#define allocate_SyWise20_only Wise2_allocate_SyWise20_only
void Wise2_init_SyWise20(SyWise20 * mat);
#define init_SyWise20 Wise2_init_SyWise20
AlnRange * Wise2_AlnRange_build_SyWise20(SyWise20 * mat,int stopj,int stopspecstate,int * startj,int * startspecstate);
#define AlnRange_build_SyWise20 Wise2_AlnRange_build_SyWise20
boolean Wise2_read_hidden_SyWise20(SyWise20 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out);
#define read_hidden_SyWise20 Wise2_read_hidden_SyWise20
int Wise2_max_hidden_SyWise20(SyWise20 * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_hidden_SyWise20 Wise2_max_hidden_SyWise20
boolean Wise2_read_special_strip_SyWise20(SyWise20 * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out);
#define read_special_strip_SyWise20 Wise2_read_special_strip_SyWise20
int Wise2_max_special_strip_SyWise20(SyWise20 * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_special_strip_SyWise20 Wise2_max_special_strip_SyWise20
int Wise2_max_matrix_to_special_SyWise20(SyWise20 * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_matrix_to_special_SyWise20 Wise2_max_matrix_to_special_SyWise20
void Wise2_calculate_hidden_SyWise20(SyWise20 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv);
#define calculate_hidden_SyWise20 Wise2_calculate_hidden_SyWise20
void Wise2_init_hidden_SyWise20(SyWise20 * mat,int starti,int startj,int stopi,int stopj);
#define init_hidden_SyWise20 Wise2_init_hidden_SyWise20
boolean Wise2_full_dc_SyWise20(SyWise20 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv);
#define full_dc_SyWise20 Wise2_full_dc_SyWise20
boolean Wise2_do_dc_single_pass_SyWise20(SyWise20 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done);
#define do_dc_single_pass_SyWise20 Wise2_do_dc_single_pass_SyWise20
void Wise2_push_dc_at_merge_SyWise20(SyWise20 * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv);
#define push_dc_at_merge_SyWise20 Wise2_push_dc_at_merge_SyWise20
void Wise2_follow_on_dc_SyWise20(SyWise20 * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define follow_on_dc_SyWise20 Wise2_follow_on_dc_SyWise20
void Wise2_run_up_dc_SyWise20(SyWise20 * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define run_up_dc_SyWise20 Wise2_run_up_dc_SyWise20
void Wise2_init_dc_SyWise20(SyWise20 * mat);
#define init_dc_SyWise20 Wise2_init_dc_SyWise20
int Wise2_start_end_find_end_SyWise20(SyWise20 * mat,int * endj);
#define start_end_find_end_SyWise20 Wise2_start_end_find_end_SyWise20
boolean Wise2_dc_optimised_start_end_calc_SyWise20(SyWise20 *mat,DPEnvelope * dpenv);
#define dc_optimised_start_end_calc_SyWise20 Wise2_dc_optimised_start_end_calc_SyWise20
void Wise2_init_start_end_linear_SyWise20(SyWise20 * mat);
#define init_start_end_linear_SyWise20 Wise2_init_start_end_linear_SyWise20
AlnConvertSet * Wise2_AlnConvertSet_SyWise20(void);
#define AlnConvertSet_SyWise20 Wise2_AlnConvertSet_SyWise20
int Wise2_SyWise20_explicit_access_main(SyWise20 * mat,int i,int j,int state);
#define SyWise20_explicit_access_main Wise2_SyWise20_explicit_access_main
int Wise2_SyWise20_explicit_access_special(SyWise20 * mat,int i,int j,int state);
#define SyWise20_explicit_access_special Wise2_SyWise20_explicit_access_special
int Wise2_find_end_SyWise20(SyWise20 * mat,int * ri,int * rj,int * state,boolean * isspecial,SyWise20_access_func_holder h);
#define find_end_SyWise20 Wise2_find_end_SyWise20
void Wise2_SyWise20_debug_show_matrix(SyWise20 * mat,int starti,int stopi,int startj,int stopj,FILE * ofp);
#define SyWise20_debug_show_matrix Wise2_SyWise20_debug_show_matrix
int Wise2_max_calc_SyWise20(SyWise20 * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,SyWise20_access_func_holder h);
#define max_calc_SyWise20 Wise2_max_calc_SyWise20
int Wise2_max_calc_special_SyWise20(SyWise20 * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,SyWise20_access_func_holder h);
#define max_calc_special_SyWise20 Wise2_max_calc_special_SyWise20

#ifdef _cplusplus
}
#endif

#endif
