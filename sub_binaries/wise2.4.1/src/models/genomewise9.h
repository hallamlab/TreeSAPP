#ifndef DYNAMITEgenomewise9HEADERFILE
#define DYNAMITEgenomewise9HEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"
#include "genome_evidence.h"

struct Wise2_GenomeWise9 {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    BaseMatrix * basematrix;     
    ShatterMatrix * shatter;     
    int leni;    
    int lenj;    
    GenomeEvidenceSet* evi;  
    ComplexSequence* gen;    
    int switchcost;  
    int newgenecost;     
    int non_start_codon;     
    int non_stop_codon;  
    RandomCodonScore * rndcodon;     
    } ;  
/* GenomeWise9 defined */ 
#ifndef DYNAMITE_DEFINED_GenomeWise9
typedef struct Wise2_GenomeWise9 Wise2_GenomeWise9;
#define GenomeWise9 Wise2_GenomeWise9
#define DYNAMITE_DEFINED_GenomeWise9
#endif


#ifdef PTHREAD
struct thread_pool_holder_GenomeWise9 {  
    GenomeEvidenceSet* evi; /* Static query data: never free'd */ 
    ComplexSequence* gen;   /* Target object placeholder */ 
    GenomicDB* targetdb;/* Target database object */ 
    boolean target_init; 
    int switchcost;  
    int newgenecost;     
    int non_start_codon;     
    int non_stop_codon;  
    RandomCodonScore * rndcodon;     
    pthread_mutex_t input_lock;  
    pthread_mutex_t output_lock;     
    Hscore * out;    
    pthread_t * pool;    
    int number_of_threads;   
    boolean search_has_ended;    
    DBSearchImpl * dbsi;     
    } ;  
#endif /* PTHREAD */
struct Wise2_GenomeWise9_access_func_holder {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int (*access_main)(GenomeWise9*,int,int,int);    
    int (*access_special)(GenomeWise9*,int,int,int); 
    } ;  
/* GenomeWise9_access_func_holder defined */ 
#ifndef DYNAMITE_DEFINED_GenomeWise9_access_func_holder
typedef struct Wise2_GenomeWise9_access_func_holder Wise2_GenomeWise9_access_func_holder;
#define GenomeWise9_access_func_holder Wise2_GenomeWise9_access_func_holder
#define DYNAMITE_DEFINED_GenomeWise9_access_func_holder
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  PackAln_read_Shatter_GenomeWise9(mat)
 *
 * Descrip:    Reads off PackAln from shatter matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [GenomeWise9 *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Shatter_GenomeWise9(GenomeWise9 * mat);
#define PackAln_read_Shatter_GenomeWise9 Wise2_PackAln_read_Shatter_GenomeWise9


/* Function:  calculate_shatter_GenomeWise9(mat,dpenv)
 *
 * Descrip:    This function calculates the GenomeWise9 matrix when in shatter mode
 *
 *
 * Arg:          mat [UNKN ] (null) [GenomeWise9 *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_shatter_GenomeWise9(GenomeWise9 * mat,DPEnvelope * dpenv);
#define calculate_shatter_GenomeWise9 Wise2_calculate_shatter_GenomeWise9


/* Function:  search_GenomeWise9(dbsi,out,evi,targetdb,switchcost,newgenecost,non_start_codon,non_stop_codon,rndcodon)
 *
 * Descrip:    This function makes a database search of GenomeWise9
 *             It uses the dbsi structure to choose which implementation
 *             to use of the database searching. This way at run time you
 *             can switch between single threaded/multi-threaded or hardware
 *
 *
 * Arg:                   dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 * Arg:                    out [UNKN ] Undocumented argument [Hscore *]
 * Arg:                    evi [UNKN ] Undocumented argument [GenomeEvidenceSet*]
 * Arg:               targetdb [UNKN ] Undocumented argument [GenomicDB*]
 * Arg:             switchcost [UNKN ] Undocumented argument [int]
 * Arg:            newgenecost [UNKN ] Undocumented argument [int]
 * Arg:        non_start_codon [UNKN ] Undocumented argument [int]
 * Arg:         non_stop_codon [UNKN ] Undocumented argument [int]
 * Arg:               rndcodon [UNKN ] Undocumented argument [RandomCodonScore *]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_search_GenomeWise9(DBSearchImpl * dbsi,Hscore * out,GenomeEvidenceSet* evi,GenomicDB* targetdb ,int switchcost,int newgenecost,int non_start_codon,int non_stop_codon,RandomCodonScore * rndcodon);
#define search_GenomeWise9 Wise2_search_GenomeWise9


/* Function:  serial_search_GenomeWise9(out,evi,targetdb,switchcost,newgenecost,non_start_codon,non_stop_codon,rndcodon)
 *
 * Descrip:    This function makes a database search of GenomeWise9
 *             It is a single processor implementation
 *
 *
 * Arg:                    out [UNKN ] Undocumented argument [Hscore *]
 * Arg:                    evi [UNKN ] Undocumented argument [GenomeEvidenceSet*]
 * Arg:               targetdb [UNKN ] Undocumented argument [GenomicDB*]
 * Arg:             switchcost [UNKN ] Undocumented argument [int]
 * Arg:            newgenecost [UNKN ] Undocumented argument [int]
 * Arg:        non_start_codon [UNKN ] Undocumented argument [int]
 * Arg:         non_stop_codon [UNKN ] Undocumented argument [int]
 * Arg:               rndcodon [UNKN ] Undocumented argument [RandomCodonScore *]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_serial_search_GenomeWise9(Hscore * out,GenomeEvidenceSet* evi,GenomicDB* targetdb ,int switchcost,int newgenecost,int non_start_codon,int non_stop_codon,RandomCodonScore * rndcodon);
#define serial_search_GenomeWise9 Wise2_serial_search_GenomeWise9


/* Function:  PackAln_bestmemory_GenomeWise9(evi,gen,switchcost,newgenecost,non_start_codon,non_stop_codon,rndcodon,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_GenomeWise9
 *
 *
 * Arg:                    evi [UNKN ] query data structure [GenomeEvidenceSet*]
 * Arg:                    gen [UNKN ] target data structure [ComplexSequence*]
 * Arg:             switchcost [UNKN ] Resource [int]
 * Arg:            newgenecost [UNKN ] Resource [int]
 * Arg:        non_start_codon [UNKN ] Resource [int]
 * Arg:         non_stop_codon [UNKN ] Resource [int]
 * Arg:               rndcodon [UNKN ] Resource [RandomCodonScore *]
 * Arg:                  dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:                   dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_bestmemory_GenomeWise9(GenomeEvidenceSet* evi,ComplexSequence* gen ,int switchcost,int newgenecost,int non_start_codon,int non_stop_codon,RandomCodonScore * rndcodon,DPEnvelope * dpenv,DPRunImpl * dpri);
#define PackAln_bestmemory_GenomeWise9 Wise2_PackAln_bestmemory_GenomeWise9


/* Function:  allocate_Expl_GenomeWise9(evi,gen,switchcost,newgenecost,non_start_codon,non_stop_codon,rndcodon,dpri)
 *
 * Descrip:    This function allocates the GenomeWise9 structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_GenomeWise9_only
 *
 *
 * Arg:                    evi [UNKN ] query data structure [GenomeEvidenceSet*]
 * Arg:                    gen [UNKN ] target data structure [ComplexSequence*]
 * Arg:             switchcost [UNKN ] Resource [int]
 * Arg:            newgenecost [UNKN ] Resource [int]
 * Arg:        non_start_codon [UNKN ] Resource [int]
 * Arg:         non_stop_codon [UNKN ] Resource [int]
 * Arg:               rndcodon [UNKN ] Resource [RandomCodonScore *]
 * Arg:                   dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [GenomeWise9 *]
 *
 */
GenomeWise9 * Wise2_allocate_Expl_GenomeWise9(GenomeEvidenceSet* evi,ComplexSequence* gen ,int switchcost,int newgenecost,int non_start_codon,int non_stop_codon,RandomCodonScore * rndcodon,DPRunImpl * dpri);
#define allocate_Expl_GenomeWise9 Wise2_allocate_Expl_GenomeWise9


/* Function:  recalculate_PackAln_GenomeWise9(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by GenomeWise9
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [GenomeWise9 *]
 *
 */
void Wise2_recalculate_PackAln_GenomeWise9(PackAln * pal,GenomeWise9 * mat);
#define recalculate_PackAln_GenomeWise9 Wise2_recalculate_PackAln_GenomeWise9


/* Function:  allocate_Small_GenomeWise9(evi,gen,switchcost,newgenecost,non_start_codon,non_stop_codon,rndcodon)
 *
 * Descrip:    This function allocates the GenomeWise9 structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_GenomeWise9_only
 *
 *
 * Arg:                    evi [UNKN ] query data structure [GenomeEvidenceSet*]
 * Arg:                    gen [UNKN ] target data structure [ComplexSequence*]
 * Arg:             switchcost [UNKN ] Resource [int]
 * Arg:            newgenecost [UNKN ] Resource [int]
 * Arg:        non_start_codon [UNKN ] Resource [int]
 * Arg:         non_stop_codon [UNKN ] Resource [int]
 * Arg:               rndcodon [UNKN ] Resource [RandomCodonScore *]
 *
 * Return [UNKN ]  Undocumented return value [GenomeWise9 *]
 *
 */
GenomeWise9 * Wise2_allocate_Small_GenomeWise9(GenomeEvidenceSet* evi,ComplexSequence* gen ,int switchcost,int newgenecost,int non_start_codon,int non_stop_codon,RandomCodonScore * rndcodon);
#define allocate_Small_GenomeWise9 Wise2_allocate_Small_GenomeWise9


/* Function:  PackAln_calculate_Small_GenomeWise9(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for GenomeWise9 structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_GenomeWise9 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_GenomeWise9 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [GenomeWise9 *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_calculate_Small_GenomeWise9(GenomeWise9 * mat,DPEnvelope * dpenv);
#define PackAln_calculate_Small_GenomeWise9 Wise2_PackAln_calculate_Small_GenomeWise9


/* Function:  AlnRangeSet_calculate_Small_GenomeWise9(mat)
 *
 * Descrip:    This function calculates an alignment for GenomeWise9 structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_GenomeWise9 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_GenomeWise9
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_GenomeWise9 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [GenomeWise9 *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_calculate_Small_GenomeWise9(GenomeWise9 * mat);
#define AlnRangeSet_calculate_Small_GenomeWise9 Wise2_AlnRangeSet_calculate_Small_GenomeWise9


/* Function:  AlnRangeSet_from_GenomeWise9(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for GenomeWise9 structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_GenomeWise9 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_GenomeWise9
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [GenomeWise9 *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_from_GenomeWise9(GenomeWise9 * mat);
#define AlnRangeSet_from_GenomeWise9 Wise2_AlnRangeSet_from_GenomeWise9


/* Function:  convert_PackAln_to_AlnBlock_GenomeWise9(pal)
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
AlnBlock * Wise2_convert_PackAln_to_AlnBlock_GenomeWise9(PackAln * pal);
#define convert_PackAln_to_AlnBlock_GenomeWise9 Wise2_convert_PackAln_to_AlnBlock_GenomeWise9


/* Function:  PackAln_read_Expl_GenomeWise9(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [GenomeWise9 *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Expl_GenomeWise9(GenomeWise9 * mat);
#define PackAln_read_Expl_GenomeWise9 Wise2_PackAln_read_Expl_GenomeWise9


/* Function:  PackAln_read_generic_GenomeWise9(mat,h)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [GenomeWise9 *]
 * Arg:          h [UNKN ] Undocumented argument [GenomeWise9_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_generic_GenomeWise9(GenomeWise9 * mat,GenomeWise9_access_func_holder h);
#define PackAln_read_generic_GenomeWise9 Wise2_PackAln_read_generic_GenomeWise9


/* Function:  calculate_GenomeWise9(mat)
 *
 * Descrip:    This function calculates the GenomeWise9 matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_GenomeWise9
 *
 *
 * Arg:        mat [UNKN ] GenomeWise9 which contains explicit basematrix memory [GenomeWise9 *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_GenomeWise9(GenomeWise9 * mat);
#define calculate_GenomeWise9 Wise2_calculate_GenomeWise9


/* Function:  calculate_dpenv_GenomeWise9(mat,dpenv)
 *
 * Descrip:    This function calculates the GenomeWise9 matrix when in explicit mode, subject to the envelope
 *
 *
 * Arg:          mat [UNKN ] GenomeWise9 which contains explicit basematrix memory [GenomeWise9 *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_dpenv_GenomeWise9(GenomeWise9 * mat,DPEnvelope * dpenv);
#define calculate_dpenv_GenomeWise9 Wise2_calculate_dpenv_GenomeWise9


/* Function:  GenomeWise9_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GenomeWise9 *]
 *
 */
GenomeWise9 * Wise2_GenomeWise9_alloc(void);
#define GenomeWise9_alloc Wise2_GenomeWise9_alloc


/* Function:  free_GenomeWise9(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GenomeWise9 *]
 *
 * Return [UNKN ]  Undocumented return value [GenomeWise9 *]
 *
 */
GenomeWise9 * Wise2_free_GenomeWise9(GenomeWise9 * obj);
#define free_GenomeWise9 Wise2_free_GenomeWise9


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
int Wise2_GenomeWise9_shatter_access_main(GenomeWise9 * mat,int i,int j,int state);
#define GenomeWise9_shatter_access_main Wise2_GenomeWise9_shatter_access_main
int Wise2_GenomeWise9_shatter_access_special(GenomeWise9 * mat,int i,int j,int state);
#define GenomeWise9_shatter_access_special Wise2_GenomeWise9_shatter_access_special
void * Wise2_thread_loop_GenomeWise9(void * ptr);
#define thread_loop_GenomeWise9 Wise2_thread_loop_GenomeWise9
int Wise2_score_only_GenomeWise9(GenomeEvidenceSet* evi,ComplexSequence* gen ,int switchcost,int newgenecost,int non_start_codon,int non_stop_codon,RandomCodonScore * rndcodon);
#define score_only_GenomeWise9 Wise2_score_only_GenomeWise9
GenomeWise9 * Wise2_allocate_GenomeWise9_only(GenomeEvidenceSet* evi,ComplexSequence* gen ,int switchcost,int newgenecost,int non_start_codon,int non_stop_codon,RandomCodonScore * rndcodon);
#define allocate_GenomeWise9_only Wise2_allocate_GenomeWise9_only
void Wise2_init_GenomeWise9(GenomeWise9 * mat);
#define init_GenomeWise9 Wise2_init_GenomeWise9
AlnRange * Wise2_AlnRange_build_GenomeWise9(GenomeWise9 * mat,int stopj,int stopspecstate,int * startj,int * startspecstate);
#define AlnRange_build_GenomeWise9 Wise2_AlnRange_build_GenomeWise9
boolean Wise2_read_hidden_GenomeWise9(GenomeWise9 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out);
#define read_hidden_GenomeWise9 Wise2_read_hidden_GenomeWise9
int Wise2_max_hidden_GenomeWise9(GenomeWise9 * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_hidden_GenomeWise9 Wise2_max_hidden_GenomeWise9
boolean Wise2_read_special_strip_GenomeWise9(GenomeWise9 * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out);
#define read_special_strip_GenomeWise9 Wise2_read_special_strip_GenomeWise9
int Wise2_max_special_strip_GenomeWise9(GenomeWise9 * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_special_strip_GenomeWise9 Wise2_max_special_strip_GenomeWise9
int Wise2_max_matrix_to_special_GenomeWise9(GenomeWise9 * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_matrix_to_special_GenomeWise9 Wise2_max_matrix_to_special_GenomeWise9
void Wise2_calculate_hidden_GenomeWise9(GenomeWise9 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv);
#define calculate_hidden_GenomeWise9 Wise2_calculate_hidden_GenomeWise9
void Wise2_init_hidden_GenomeWise9(GenomeWise9 * mat,int starti,int startj,int stopi,int stopj);
#define init_hidden_GenomeWise9 Wise2_init_hidden_GenomeWise9
boolean Wise2_full_dc_GenomeWise9(GenomeWise9 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv);
#define full_dc_GenomeWise9 Wise2_full_dc_GenomeWise9
boolean Wise2_do_dc_single_pass_GenomeWise9(GenomeWise9 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done);
#define do_dc_single_pass_GenomeWise9 Wise2_do_dc_single_pass_GenomeWise9
void Wise2_push_dc_at_merge_GenomeWise9(GenomeWise9 * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv);
#define push_dc_at_merge_GenomeWise9 Wise2_push_dc_at_merge_GenomeWise9
void Wise2_follow_on_dc_GenomeWise9(GenomeWise9 * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define follow_on_dc_GenomeWise9 Wise2_follow_on_dc_GenomeWise9
void Wise2_run_up_dc_GenomeWise9(GenomeWise9 * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define run_up_dc_GenomeWise9 Wise2_run_up_dc_GenomeWise9
void Wise2_init_dc_GenomeWise9(GenomeWise9 * mat);
#define init_dc_GenomeWise9 Wise2_init_dc_GenomeWise9
int Wise2_start_end_find_end_GenomeWise9(GenomeWise9 * mat,int * endj);
#define start_end_find_end_GenomeWise9 Wise2_start_end_find_end_GenomeWise9
boolean Wise2_dc_optimised_start_end_calc_GenomeWise9(GenomeWise9 *mat,DPEnvelope * dpenv);
#define dc_optimised_start_end_calc_GenomeWise9 Wise2_dc_optimised_start_end_calc_GenomeWise9
void Wise2_init_start_end_linear_GenomeWise9(GenomeWise9 * mat);
#define init_start_end_linear_GenomeWise9 Wise2_init_start_end_linear_GenomeWise9
AlnConvertSet * Wise2_AlnConvertSet_GenomeWise9(void);
#define AlnConvertSet_GenomeWise9 Wise2_AlnConvertSet_GenomeWise9
int Wise2_GenomeWise9_explicit_access_main(GenomeWise9 * mat,int i,int j,int state);
#define GenomeWise9_explicit_access_main Wise2_GenomeWise9_explicit_access_main
int Wise2_GenomeWise9_explicit_access_special(GenomeWise9 * mat,int i,int j,int state);
#define GenomeWise9_explicit_access_special Wise2_GenomeWise9_explicit_access_special
int Wise2_find_end_GenomeWise9(GenomeWise9 * mat,int * ri,int * rj,int * state,boolean * isspecial,GenomeWise9_access_func_holder h);
#define find_end_GenomeWise9 Wise2_find_end_GenomeWise9
void Wise2_GenomeWise9_debug_show_matrix(GenomeWise9 * mat,int starti,int stopi,int startj,int stopj,FILE * ofp);
#define GenomeWise9_debug_show_matrix Wise2_GenomeWise9_debug_show_matrix
int Wise2_max_calc_GenomeWise9(GenomeWise9 * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,GenomeWise9_access_func_holder h);
#define max_calc_GenomeWise9 Wise2_max_calc_GenomeWise9
int Wise2_max_calc_special_GenomeWise9(GenomeWise9 * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,GenomeWise9_access_func_holder h);
#define max_calc_special_GenomeWise9 Wise2_max_calc_special_GenomeWise9

#ifdef _cplusplus
}
#endif

#endif
