#ifndef DYNAMITEcdnawise10HEADERFILE
#define DYNAMITEcdnawise10HEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"
#include "geneparser4.h"


struct Wise2_CdnaWise10 {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    BaseMatrix * basematrix;     
    ShatterMatrix * shatter;     
    int leni;    
    int lenj;    
    ComplexSequence* query;  
    ComplexSequence* target;     
    GeneParser4Score * gp;   
    CodonMatrixScore * sc;   
    RandomCodonScore * rndcodon;     
    DnaMatrix* utr;  
    Score gap_open;  
    Score gap_ext;   
    Score utr_gap;   
    Score intron_gap;    
    } ;  
/* CdnaWise10 defined */ 
#ifndef DYNAMITE_DEFINED_CdnaWise10
typedef struct Wise2_CdnaWise10 Wise2_CdnaWise10;
#define CdnaWise10 Wise2_CdnaWise10
#define DYNAMITE_DEFINED_CdnaWise10
#endif


#ifdef PTHREAD
struct thread_pool_holder_CdnaWise10 {  
    ComplexSequence* query; /* Query object placeholder */ 
    cDNADB* querydb;/* Query database object */ 
    boolean query_init;  
    ComplexSequence* target;/* Target object placeholder */ 
    GenomicDB* targetdb;/* Target database object */ 
    boolean target_init; 
    GeneParser4Score * gp;   
    CodonMatrixScore * sc;   
    RandomCodonScore * rndcodon;     
    DnaMatrix* utr;  
    Score gap_open;  
    Score gap_ext;   
    Score utr_gap;   
    Score intron_gap;    
    pthread_mutex_t input_lock;  
    pthread_mutex_t output_lock;     
    Hscore * out;    
    pthread_t * pool;    
    int number_of_threads;   
    boolean search_has_ended;    
    DBSearchImpl * dbsi;     
    } ;  
#endif /* PTHREAD */
struct Wise2_CdnaWise10_access_func_holder {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int (*access_main)(CdnaWise10*,int,int,int); 
    int (*access_special)(CdnaWise10*,int,int,int);  
    } ;  
/* CdnaWise10_access_func_holder defined */ 
#ifndef DYNAMITE_DEFINED_CdnaWise10_access_func_holder
typedef struct Wise2_CdnaWise10_access_func_holder Wise2_CdnaWise10_access_func_holder;
#define CdnaWise10_access_func_holder Wise2_CdnaWise10_access_func_holder
#define DYNAMITE_DEFINED_CdnaWise10_access_func_holder
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  PackAln_read_Shatter_CdnaWise10(mat)
 *
 * Descrip:    Reads off PackAln from shatter matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [CdnaWise10 *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Shatter_CdnaWise10(CdnaWise10 * mat);
#define PackAln_read_Shatter_CdnaWise10 Wise2_PackAln_read_Shatter_CdnaWise10


/* Function:  calculate_shatter_CdnaWise10(mat,dpenv)
 *
 * Descrip:    This function calculates the CdnaWise10 matrix when in shatter mode
 *
 *
 * Arg:          mat [UNKN ] (null) [CdnaWise10 *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_shatter_CdnaWise10(CdnaWise10 * mat,DPEnvelope * dpenv);
#define calculate_shatter_CdnaWise10 Wise2_calculate_shatter_CdnaWise10


/* Function:  search_CdnaWise10(dbsi,out,querydb,targetdb,gp,sc,rndcodon,utr,gap_open,gap_ext,utr_gap,intron_gap)
 *
 * Descrip:    This function makes a database search of CdnaWise10
 *             It uses the dbsi structure to choose which implementation
 *             to use of the database searching. This way at run time you
 *             can switch between single threaded/multi-threaded or hardware
 *
 *
 * Arg:              dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 * Arg:               out [UNKN ] Undocumented argument [Hscore *]
 * Arg:           querydb [UNKN ] Undocumented argument [cDNADB*]
 * Arg:          targetdb [UNKN ] Undocumented argument [GenomicDB*]
 * Arg:                gp [UNKN ] Undocumented argument [GeneParser4Score *]
 * Arg:                sc [UNKN ] Undocumented argument [CodonMatrixScore *]
 * Arg:          rndcodon [UNKN ] Undocumented argument [RandomCodonScore *]
 * Arg:               utr [UNKN ] Undocumented argument [DnaMatrix*]
 * Arg:          gap_open [UNKN ] Undocumented argument [Score]
 * Arg:           gap_ext [UNKN ] Undocumented argument [Score]
 * Arg:           utr_gap [UNKN ] Undocumented argument [Score]
 * Arg:        intron_gap [UNKN ] Undocumented argument [Score]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_search_CdnaWise10(DBSearchImpl * dbsi,Hscore * out,cDNADB* querydb,GenomicDB* targetdb ,GeneParser4Score * gp,CodonMatrixScore * sc,RandomCodonScore * rndcodon,DnaMatrix* utr,Score gap_open,Score gap_ext,Score utr_gap,Score intron_gap);
#define search_CdnaWise10 Wise2_search_CdnaWise10


/* Function:  serial_search_CdnaWise10(out,querydb,targetdb,gp,sc,rndcodon,utr,gap_open,gap_ext,utr_gap,intron_gap)
 *
 * Descrip:    This function makes a database search of CdnaWise10
 *             It is a single processor implementation
 *
 *
 * Arg:               out [UNKN ] Undocumented argument [Hscore *]
 * Arg:           querydb [UNKN ] Undocumented argument [cDNADB*]
 * Arg:          targetdb [UNKN ] Undocumented argument [GenomicDB*]
 * Arg:                gp [UNKN ] Undocumented argument [GeneParser4Score *]
 * Arg:                sc [UNKN ] Undocumented argument [CodonMatrixScore *]
 * Arg:          rndcodon [UNKN ] Undocumented argument [RandomCodonScore *]
 * Arg:               utr [UNKN ] Undocumented argument [DnaMatrix*]
 * Arg:          gap_open [UNKN ] Undocumented argument [Score]
 * Arg:           gap_ext [UNKN ] Undocumented argument [Score]
 * Arg:           utr_gap [UNKN ] Undocumented argument [Score]
 * Arg:        intron_gap [UNKN ] Undocumented argument [Score]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_serial_search_CdnaWise10(Hscore * out,cDNADB* querydb,GenomicDB* targetdb ,GeneParser4Score * gp,CodonMatrixScore * sc,RandomCodonScore * rndcodon,DnaMatrix* utr,Score gap_open,Score gap_ext,Score utr_gap,Score intron_gap);
#define serial_search_CdnaWise10 Wise2_serial_search_CdnaWise10


/* Function:  PackAln_bestmemory_CdnaWise10(query,target,gp,sc,rndcodon,utr,gap_open,gap_ext,utr_gap,intron_gap,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_CdnaWise10
 *
 *
 * Arg:             query [UNKN ] query data structure [ComplexSequence*]
 * Arg:            target [UNKN ] target data structure [ComplexSequence*]
 * Arg:                gp [UNKN ] Resource [GeneParser4Score *]
 * Arg:                sc [UNKN ] Resource [CodonMatrixScore *]
 * Arg:          rndcodon [UNKN ] Resource [RandomCodonScore *]
 * Arg:               utr [UNKN ] Resource [DnaMatrix*]
 * Arg:          gap_open [UNKN ] Resource [Score]
 * Arg:           gap_ext [UNKN ] Resource [Score]
 * Arg:           utr_gap [UNKN ] Resource [Score]
 * Arg:        intron_gap [UNKN ] Resource [Score]
 * Arg:             dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:              dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_bestmemory_CdnaWise10(ComplexSequence* query,ComplexSequence* target ,GeneParser4Score * gp,CodonMatrixScore * sc,RandomCodonScore * rndcodon,DnaMatrix* utr,Score gap_open,Score gap_ext,Score utr_gap,Score intron_gap,DPEnvelope * dpenv,DPRunImpl * dpri);
#define PackAln_bestmemory_CdnaWise10 Wise2_PackAln_bestmemory_CdnaWise10


/* Function:  allocate_Expl_CdnaWise10(query,target,gp,sc,rndcodon,utr,gap_open,gap_ext,utr_gap,intron_gap,dpri)
 *
 * Descrip:    This function allocates the CdnaWise10 structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_CdnaWise10_only
 *
 *
 * Arg:             query [UNKN ] query data structure [ComplexSequence*]
 * Arg:            target [UNKN ] target data structure [ComplexSequence*]
 * Arg:                gp [UNKN ] Resource [GeneParser4Score *]
 * Arg:                sc [UNKN ] Resource [CodonMatrixScore *]
 * Arg:          rndcodon [UNKN ] Resource [RandomCodonScore *]
 * Arg:               utr [UNKN ] Resource [DnaMatrix*]
 * Arg:          gap_open [UNKN ] Resource [Score]
 * Arg:           gap_ext [UNKN ] Resource [Score]
 * Arg:           utr_gap [UNKN ] Resource [Score]
 * Arg:        intron_gap [UNKN ] Resource [Score]
 * Arg:              dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [CdnaWise10 *]
 *
 */
CdnaWise10 * Wise2_allocate_Expl_CdnaWise10(ComplexSequence* query,ComplexSequence* target ,GeneParser4Score * gp,CodonMatrixScore * sc,RandomCodonScore * rndcodon,DnaMatrix* utr,Score gap_open,Score gap_ext,Score utr_gap,Score intron_gap,DPRunImpl * dpri);
#define allocate_Expl_CdnaWise10 Wise2_allocate_Expl_CdnaWise10


/* Function:  recalculate_PackAln_CdnaWise10(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by CdnaWise10
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [CdnaWise10 *]
 *
 */
void Wise2_recalculate_PackAln_CdnaWise10(PackAln * pal,CdnaWise10 * mat);
#define recalculate_PackAln_CdnaWise10 Wise2_recalculate_PackAln_CdnaWise10


/* Function:  allocate_Small_CdnaWise10(query,target,gp,sc,rndcodon,utr,gap_open,gap_ext,utr_gap,intron_gap)
 *
 * Descrip:    This function allocates the CdnaWise10 structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_CdnaWise10_only
 *
 *
 * Arg:             query [UNKN ] query data structure [ComplexSequence*]
 * Arg:            target [UNKN ] target data structure [ComplexSequence*]
 * Arg:                gp [UNKN ] Resource [GeneParser4Score *]
 * Arg:                sc [UNKN ] Resource [CodonMatrixScore *]
 * Arg:          rndcodon [UNKN ] Resource [RandomCodonScore *]
 * Arg:               utr [UNKN ] Resource [DnaMatrix*]
 * Arg:          gap_open [UNKN ] Resource [Score]
 * Arg:           gap_ext [UNKN ] Resource [Score]
 * Arg:           utr_gap [UNKN ] Resource [Score]
 * Arg:        intron_gap [UNKN ] Resource [Score]
 *
 * Return [UNKN ]  Undocumented return value [CdnaWise10 *]
 *
 */
CdnaWise10 * Wise2_allocate_Small_CdnaWise10(ComplexSequence* query,ComplexSequence* target ,GeneParser4Score * gp,CodonMatrixScore * sc,RandomCodonScore * rndcodon,DnaMatrix* utr,Score gap_open,Score gap_ext,Score utr_gap,Score intron_gap);
#define allocate_Small_CdnaWise10 Wise2_allocate_Small_CdnaWise10


/* Function:  PackAln_calculate_Small_CdnaWise10(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for CdnaWise10 structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_CdnaWise10 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_CdnaWise10 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [CdnaWise10 *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_calculate_Small_CdnaWise10(CdnaWise10 * mat,DPEnvelope * dpenv);
#define PackAln_calculate_Small_CdnaWise10 Wise2_PackAln_calculate_Small_CdnaWise10


/* Function:  AlnRangeSet_calculate_Small_CdnaWise10(mat)
 *
 * Descrip:    This function calculates an alignment for CdnaWise10 structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_CdnaWise10 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_CdnaWise10
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_CdnaWise10 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [CdnaWise10 *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_calculate_Small_CdnaWise10(CdnaWise10 * mat);
#define AlnRangeSet_calculate_Small_CdnaWise10 Wise2_AlnRangeSet_calculate_Small_CdnaWise10


/* Function:  AlnRangeSet_from_CdnaWise10(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for CdnaWise10 structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_CdnaWise10 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_CdnaWise10
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [CdnaWise10 *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_from_CdnaWise10(CdnaWise10 * mat);
#define AlnRangeSet_from_CdnaWise10 Wise2_AlnRangeSet_from_CdnaWise10


/* Function:  convert_PackAln_to_AlnBlock_CdnaWise10(pal)
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
AlnBlock * Wise2_convert_PackAln_to_AlnBlock_CdnaWise10(PackAln * pal);
#define convert_PackAln_to_AlnBlock_CdnaWise10 Wise2_convert_PackAln_to_AlnBlock_CdnaWise10


/* Function:  PackAln_read_Expl_CdnaWise10(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [CdnaWise10 *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Expl_CdnaWise10(CdnaWise10 * mat);
#define PackAln_read_Expl_CdnaWise10 Wise2_PackAln_read_Expl_CdnaWise10


/* Function:  PackAln_read_generic_CdnaWise10(mat,h)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [CdnaWise10 *]
 * Arg:          h [UNKN ] Undocumented argument [CdnaWise10_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_generic_CdnaWise10(CdnaWise10 * mat,CdnaWise10_access_func_holder h);
#define PackAln_read_generic_CdnaWise10 Wise2_PackAln_read_generic_CdnaWise10


/* Function:  calculate_CdnaWise10(mat)
 *
 * Descrip:    This function calculates the CdnaWise10 matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_CdnaWise10
 *
 *
 * Arg:        mat [UNKN ] CdnaWise10 which contains explicit basematrix memory [CdnaWise10 *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_CdnaWise10(CdnaWise10 * mat);
#define calculate_CdnaWise10 Wise2_calculate_CdnaWise10


/* Function:  calculate_dpenv_CdnaWise10(mat,dpenv)
 *
 * Descrip:    This function calculates the CdnaWise10 matrix when in explicit mode, subject to the envelope
 *
 *
 * Arg:          mat [UNKN ] CdnaWise10 which contains explicit basematrix memory [CdnaWise10 *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_dpenv_CdnaWise10(CdnaWise10 * mat,DPEnvelope * dpenv);
#define calculate_dpenv_CdnaWise10 Wise2_calculate_dpenv_CdnaWise10


/* Function:  CdnaWise10_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [CdnaWise10 *]
 *
 */
CdnaWise10 * Wise2_CdnaWise10_alloc(void);
#define CdnaWise10_alloc Wise2_CdnaWise10_alloc


/* Function:  free_CdnaWise10(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [CdnaWise10 *]
 *
 * Return [UNKN ]  Undocumented return value [CdnaWise10 *]
 *
 */
CdnaWise10 * Wise2_free_CdnaWise10(CdnaWise10 * obj);
#define free_CdnaWise10 Wise2_free_CdnaWise10


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
int Wise2_CdnaWise10_shatter_access_main(CdnaWise10 * mat,int i,int j,int state);
#define CdnaWise10_shatter_access_main Wise2_CdnaWise10_shatter_access_main
int Wise2_CdnaWise10_shatter_access_special(CdnaWise10 * mat,int i,int j,int state);
#define CdnaWise10_shatter_access_special Wise2_CdnaWise10_shatter_access_special
void * Wise2_thread_loop_CdnaWise10(void * ptr);
#define thread_loop_CdnaWise10 Wise2_thread_loop_CdnaWise10
int Wise2_score_only_CdnaWise10(ComplexSequence* query,ComplexSequence* target ,GeneParser4Score * gp,CodonMatrixScore * sc,RandomCodonScore * rndcodon,DnaMatrix* utr,Score gap_open,Score gap_ext,Score utr_gap,Score intron_gap);
#define score_only_CdnaWise10 Wise2_score_only_CdnaWise10
CdnaWise10 * Wise2_allocate_CdnaWise10_only(ComplexSequence* query,ComplexSequence* target ,GeneParser4Score * gp,CodonMatrixScore * sc,RandomCodonScore * rndcodon,DnaMatrix* utr,Score gap_open,Score gap_ext,Score utr_gap,Score intron_gap);
#define allocate_CdnaWise10_only Wise2_allocate_CdnaWise10_only
void Wise2_init_CdnaWise10(CdnaWise10 * mat);
#define init_CdnaWise10 Wise2_init_CdnaWise10
AlnRange * Wise2_AlnRange_build_CdnaWise10(CdnaWise10 * mat,int stopj,int stopspecstate,int * startj,int * startspecstate);
#define AlnRange_build_CdnaWise10 Wise2_AlnRange_build_CdnaWise10
boolean Wise2_read_hidden_CdnaWise10(CdnaWise10 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out);
#define read_hidden_CdnaWise10 Wise2_read_hidden_CdnaWise10
int Wise2_max_hidden_CdnaWise10(CdnaWise10 * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_hidden_CdnaWise10 Wise2_max_hidden_CdnaWise10
boolean Wise2_read_special_strip_CdnaWise10(CdnaWise10 * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out);
#define read_special_strip_CdnaWise10 Wise2_read_special_strip_CdnaWise10
int Wise2_max_special_strip_CdnaWise10(CdnaWise10 * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_special_strip_CdnaWise10 Wise2_max_special_strip_CdnaWise10
int Wise2_max_matrix_to_special_CdnaWise10(CdnaWise10 * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_matrix_to_special_CdnaWise10 Wise2_max_matrix_to_special_CdnaWise10
void Wise2_calculate_hidden_CdnaWise10(CdnaWise10 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv);
#define calculate_hidden_CdnaWise10 Wise2_calculate_hidden_CdnaWise10
void Wise2_init_hidden_CdnaWise10(CdnaWise10 * mat,int starti,int startj,int stopi,int stopj);
#define init_hidden_CdnaWise10 Wise2_init_hidden_CdnaWise10
boolean Wise2_full_dc_CdnaWise10(CdnaWise10 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv);
#define full_dc_CdnaWise10 Wise2_full_dc_CdnaWise10
boolean Wise2_do_dc_single_pass_CdnaWise10(CdnaWise10 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done);
#define do_dc_single_pass_CdnaWise10 Wise2_do_dc_single_pass_CdnaWise10
void Wise2_push_dc_at_merge_CdnaWise10(CdnaWise10 * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv);
#define push_dc_at_merge_CdnaWise10 Wise2_push_dc_at_merge_CdnaWise10
void Wise2_follow_on_dc_CdnaWise10(CdnaWise10 * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define follow_on_dc_CdnaWise10 Wise2_follow_on_dc_CdnaWise10
void Wise2_run_up_dc_CdnaWise10(CdnaWise10 * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define run_up_dc_CdnaWise10 Wise2_run_up_dc_CdnaWise10
void Wise2_init_dc_CdnaWise10(CdnaWise10 * mat);
#define init_dc_CdnaWise10 Wise2_init_dc_CdnaWise10
int Wise2_start_end_find_end_CdnaWise10(CdnaWise10 * mat,int * endj);
#define start_end_find_end_CdnaWise10 Wise2_start_end_find_end_CdnaWise10
boolean Wise2_dc_optimised_start_end_calc_CdnaWise10(CdnaWise10 *mat,DPEnvelope * dpenv);
#define dc_optimised_start_end_calc_CdnaWise10 Wise2_dc_optimised_start_end_calc_CdnaWise10
void Wise2_init_start_end_linear_CdnaWise10(CdnaWise10 * mat);
#define init_start_end_linear_CdnaWise10 Wise2_init_start_end_linear_CdnaWise10
AlnConvertSet * Wise2_AlnConvertSet_CdnaWise10(void);
#define AlnConvertSet_CdnaWise10 Wise2_AlnConvertSet_CdnaWise10
int Wise2_CdnaWise10_explicit_access_main(CdnaWise10 * mat,int i,int j,int state);
#define CdnaWise10_explicit_access_main Wise2_CdnaWise10_explicit_access_main
int Wise2_CdnaWise10_explicit_access_special(CdnaWise10 * mat,int i,int j,int state);
#define CdnaWise10_explicit_access_special Wise2_CdnaWise10_explicit_access_special
int Wise2_find_end_CdnaWise10(CdnaWise10 * mat,int * ri,int * rj,int * state,boolean * isspecial,CdnaWise10_access_func_holder h);
#define find_end_CdnaWise10 Wise2_find_end_CdnaWise10
void Wise2_CdnaWise10_debug_show_matrix(CdnaWise10 * mat,int starti,int stopi,int startj,int stopj,FILE * ofp);
#define CdnaWise10_debug_show_matrix Wise2_CdnaWise10_debug_show_matrix
int Wise2_max_calc_CdnaWise10(CdnaWise10 * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,CdnaWise10_access_func_holder h);
#define max_calc_CdnaWise10 Wise2_max_calc_CdnaWise10
int Wise2_max_calc_special_CdnaWise10(CdnaWise10 * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,CdnaWise10_access_func_holder h);
#define max_calc_special_CdnaWise10 Wise2_max_calc_special_CdnaWise10

#ifdef _cplusplus
}
#endif

#endif
