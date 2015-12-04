#ifndef DYNAMITEbigdbaHEADERFILE
#define DYNAMITEbigdbaHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"

struct Wise2_BigDnaMatchBlock {  
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
    DnaMatrix* comp55;   
    DnaMatrix* comp65;   
    DnaMatrix* comp75;   
    DnaMatrix* comp85;   
    DnaMatrix* comp95;   
    Score g;     
    Score g65;   
    Score g55;   
    Score u;     
    Score v;     
    Score s;     
    Score b;     
    } ;  
/* BigDnaMatchBlock defined */ 
#ifndef DYNAMITE_DEFINED_BigDnaMatchBlock
typedef struct Wise2_BigDnaMatchBlock Wise2_BigDnaMatchBlock;
#define BigDnaMatchBlock Wise2_BigDnaMatchBlock
#define DYNAMITE_DEFINED_BigDnaMatchBlock
#endif


struct Wise2_BigDnaMatchBlock_access_func_holder {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int (*access_main)(BigDnaMatchBlock*,int,int,int);   
    int (*access_special)(BigDnaMatchBlock*,int,int,int);    
    } ;  
/* BigDnaMatchBlock_access_func_holder defined */ 
#ifndef DYNAMITE_DEFINED_BigDnaMatchBlock_access_func_holder
typedef struct Wise2_BigDnaMatchBlock_access_func_holder Wise2_BigDnaMatchBlock_access_func_holder;
#define BigDnaMatchBlock_access_func_holder Wise2_BigDnaMatchBlock_access_func_holder
#define DYNAMITE_DEFINED_BigDnaMatchBlock_access_func_holder
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  PackAln_read_Shatter_BigDnaMatchBlock(mat)
 *
 * Descrip:    Reads off PackAln from shatter matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [BigDnaMatchBlock *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Shatter_BigDnaMatchBlock(BigDnaMatchBlock * mat);
#define PackAln_read_Shatter_BigDnaMatchBlock Wise2_PackAln_read_Shatter_BigDnaMatchBlock


/* Function:  calculate_shatter_BigDnaMatchBlock(mat,dpenv)
 *
 * Descrip:    This function calculates the BigDnaMatchBlock matrix when in shatter mode
 *
 *
 * Arg:          mat [UNKN ] (null) [BigDnaMatchBlock *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_shatter_BigDnaMatchBlock(BigDnaMatchBlock * mat,DPEnvelope * dpenv);
#define calculate_shatter_BigDnaMatchBlock Wise2_calculate_shatter_BigDnaMatchBlock


/* Function:  search_BigDnaMatchBlock(dbsi,out,query,target,comp55,comp65,comp75,comp85,comp95,g,g65,g55,u,v,s,b)
 *
 * Descrip:    This function makes a database search of BigDnaMatchBlock
 *             It uses the dbsi structure to choose which implementation
 *             to use of the database searching. This way at run time you
 *             can switch between single threaded/multi-threaded or hardware
 *
 *
 * Arg:          dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 * Arg:           out [UNKN ] Undocumented argument [Hscore *]
 * Arg:         query [UNKN ] Undocumented argument [ComplexSequence*]
 * Arg:        target [UNKN ] Undocumented argument [ComplexSequence*]
 * Arg:        comp55 [UNKN ] Undocumented argument [DnaMatrix*]
 * Arg:        comp65 [UNKN ] Undocumented argument [DnaMatrix*]
 * Arg:        comp75 [UNKN ] Undocumented argument [DnaMatrix*]
 * Arg:        comp85 [UNKN ] Undocumented argument [DnaMatrix*]
 * Arg:        comp95 [UNKN ] Undocumented argument [DnaMatrix*]
 * Arg:             g [UNKN ] Undocumented argument [Score]
 * Arg:           g65 [UNKN ] Undocumented argument [Score]
 * Arg:           g55 [UNKN ] Undocumented argument [Score]
 * Arg:             u [UNKN ] Undocumented argument [Score]
 * Arg:             v [UNKN ] Undocumented argument [Score]
 * Arg:             s [UNKN ] Undocumented argument [Score]
 * Arg:             b [UNKN ] Undocumented argument [Score]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_search_BigDnaMatchBlock(DBSearchImpl * dbsi,Hscore * out,ComplexSequence* query,ComplexSequence* target ,DnaMatrix* comp55,DnaMatrix* comp65,DnaMatrix* comp75,DnaMatrix* comp85,DnaMatrix* comp95,Score g,Score g65,Score g55,Score u,Score v,Score s,Score b);
#define search_BigDnaMatchBlock Wise2_search_BigDnaMatchBlock


/* Function:  serial_search_BigDnaMatchBlock(out,query,target,comp55,comp65,comp75,comp85,comp95,g,g65,g55,u,v,s,b)
 *
 * Descrip:    This function makes a database search of BigDnaMatchBlock
 *             It is a single processor implementation
 *
 *
 * Arg:           out [UNKN ] Undocumented argument [Hscore *]
 * Arg:         query [UNKN ] Undocumented argument [ComplexSequence*]
 * Arg:        target [UNKN ] Undocumented argument [ComplexSequence*]
 * Arg:        comp55 [UNKN ] Undocumented argument [DnaMatrix*]
 * Arg:        comp65 [UNKN ] Undocumented argument [DnaMatrix*]
 * Arg:        comp75 [UNKN ] Undocumented argument [DnaMatrix*]
 * Arg:        comp85 [UNKN ] Undocumented argument [DnaMatrix*]
 * Arg:        comp95 [UNKN ] Undocumented argument [DnaMatrix*]
 * Arg:             g [UNKN ] Undocumented argument [Score]
 * Arg:           g65 [UNKN ] Undocumented argument [Score]
 * Arg:           g55 [UNKN ] Undocumented argument [Score]
 * Arg:             u [UNKN ] Undocumented argument [Score]
 * Arg:             v [UNKN ] Undocumented argument [Score]
 * Arg:             s [UNKN ] Undocumented argument [Score]
 * Arg:             b [UNKN ] Undocumented argument [Score]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_serial_search_BigDnaMatchBlock(Hscore * out,ComplexSequence* query,ComplexSequence* target ,DnaMatrix* comp55,DnaMatrix* comp65,DnaMatrix* comp75,DnaMatrix* comp85,DnaMatrix* comp95,Score g,Score g65,Score g55,Score u,Score v,Score s,Score b);
#define serial_search_BigDnaMatchBlock Wise2_serial_search_BigDnaMatchBlock


/* Function:  PackAln_bestmemory_BigDnaMatchBlock(query,target,comp55,comp65,comp75,comp85,comp95,g,g65,g55,u,v,s,b,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_BigDnaMatchBlock
 *
 *
 * Arg:         query [UNKN ] query data structure [ComplexSequence*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:        comp55 [UNKN ] Resource [DnaMatrix*]
 * Arg:        comp65 [UNKN ] Resource [DnaMatrix*]
 * Arg:        comp75 [UNKN ] Resource [DnaMatrix*]
 * Arg:        comp85 [UNKN ] Resource [DnaMatrix*]
 * Arg:        comp95 [UNKN ] Resource [DnaMatrix*]
 * Arg:             g [UNKN ] Resource [Score]
 * Arg:           g65 [UNKN ] Resource [Score]
 * Arg:           g55 [UNKN ] Resource [Score]
 * Arg:             u [UNKN ] Resource [Score]
 * Arg:             v [UNKN ] Resource [Score]
 * Arg:             s [UNKN ] Resource [Score]
 * Arg:             b [UNKN ] Resource [Score]
 * Arg:         dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:          dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_bestmemory_BigDnaMatchBlock(ComplexSequence* query,ComplexSequence* target ,DnaMatrix* comp55,DnaMatrix* comp65,DnaMatrix* comp75,DnaMatrix* comp85,DnaMatrix* comp95,Score g,Score g65,Score g55,Score u,Score v,Score s,Score b,DPEnvelope * dpenv,DPRunImpl * dpri);
#define PackAln_bestmemory_BigDnaMatchBlock Wise2_PackAln_bestmemory_BigDnaMatchBlock


/* Function:  allocate_Expl_BigDnaMatchBlock(query,target,comp55,comp65,comp75,comp85,comp95,g,g65,g55,u,v,s,b,dpri)
 *
 * Descrip:    This function allocates the BigDnaMatchBlock structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_BigDnaMatchBlock_only
 *
 *
 * Arg:         query [UNKN ] query data structure [ComplexSequence*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:        comp55 [UNKN ] Resource [DnaMatrix*]
 * Arg:        comp65 [UNKN ] Resource [DnaMatrix*]
 * Arg:        comp75 [UNKN ] Resource [DnaMatrix*]
 * Arg:        comp85 [UNKN ] Resource [DnaMatrix*]
 * Arg:        comp95 [UNKN ] Resource [DnaMatrix*]
 * Arg:             g [UNKN ] Resource [Score]
 * Arg:           g65 [UNKN ] Resource [Score]
 * Arg:           g55 [UNKN ] Resource [Score]
 * Arg:             u [UNKN ] Resource [Score]
 * Arg:             v [UNKN ] Resource [Score]
 * Arg:             s [UNKN ] Resource [Score]
 * Arg:             b [UNKN ] Resource [Score]
 * Arg:          dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [BigDnaMatchBlock *]
 *
 */
BigDnaMatchBlock * Wise2_allocate_Expl_BigDnaMatchBlock(ComplexSequence* query,ComplexSequence* target ,DnaMatrix* comp55,DnaMatrix* comp65,DnaMatrix* comp75,DnaMatrix* comp85,DnaMatrix* comp95,Score g,Score g65,Score g55,Score u,Score v,Score s,Score b,DPRunImpl * dpri);
#define allocate_Expl_BigDnaMatchBlock Wise2_allocate_Expl_BigDnaMatchBlock


/* Function:  recalculate_PackAln_BigDnaMatchBlock(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by BigDnaMatchBlock
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [BigDnaMatchBlock *]
 *
 */
void Wise2_recalculate_PackAln_BigDnaMatchBlock(PackAln * pal,BigDnaMatchBlock * mat);
#define recalculate_PackAln_BigDnaMatchBlock Wise2_recalculate_PackAln_BigDnaMatchBlock


/* Function:  allocate_Small_BigDnaMatchBlock(query,target,comp55,comp65,comp75,comp85,comp95,g,g65,g55,u,v,s,b)
 *
 * Descrip:    This function allocates the BigDnaMatchBlock structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_BigDnaMatchBlock_only
 *
 *
 * Arg:         query [UNKN ] query data structure [ComplexSequence*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:        comp55 [UNKN ] Resource [DnaMatrix*]
 * Arg:        comp65 [UNKN ] Resource [DnaMatrix*]
 * Arg:        comp75 [UNKN ] Resource [DnaMatrix*]
 * Arg:        comp85 [UNKN ] Resource [DnaMatrix*]
 * Arg:        comp95 [UNKN ] Resource [DnaMatrix*]
 * Arg:             g [UNKN ] Resource [Score]
 * Arg:           g65 [UNKN ] Resource [Score]
 * Arg:           g55 [UNKN ] Resource [Score]
 * Arg:             u [UNKN ] Resource [Score]
 * Arg:             v [UNKN ] Resource [Score]
 * Arg:             s [UNKN ] Resource [Score]
 * Arg:             b [UNKN ] Resource [Score]
 *
 * Return [UNKN ]  Undocumented return value [BigDnaMatchBlock *]
 *
 */
BigDnaMatchBlock * Wise2_allocate_Small_BigDnaMatchBlock(ComplexSequence* query,ComplexSequence* target ,DnaMatrix* comp55,DnaMatrix* comp65,DnaMatrix* comp75,DnaMatrix* comp85,DnaMatrix* comp95,Score g,Score g65,Score g55,Score u,Score v,Score s,Score b);
#define allocate_Small_BigDnaMatchBlock Wise2_allocate_Small_BigDnaMatchBlock


/* Function:  PackAln_calculate_Small_BigDnaMatchBlock(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for BigDnaMatchBlock structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_BigDnaMatchBlock 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_BigDnaMatchBlock 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [BigDnaMatchBlock *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_calculate_Small_BigDnaMatchBlock(BigDnaMatchBlock * mat,DPEnvelope * dpenv);
#define PackAln_calculate_Small_BigDnaMatchBlock Wise2_PackAln_calculate_Small_BigDnaMatchBlock


/* Function:  AlnRangeSet_calculate_Small_BigDnaMatchBlock(mat)
 *
 * Descrip:    This function calculates an alignment for BigDnaMatchBlock structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_BigDnaMatchBlock 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_BigDnaMatchBlock
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_BigDnaMatchBlock 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [BigDnaMatchBlock *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_calculate_Small_BigDnaMatchBlock(BigDnaMatchBlock * mat);
#define AlnRangeSet_calculate_Small_BigDnaMatchBlock Wise2_AlnRangeSet_calculate_Small_BigDnaMatchBlock


/* Function:  AlnRangeSet_from_BigDnaMatchBlock(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for BigDnaMatchBlock structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_BigDnaMatchBlock 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_BigDnaMatchBlock
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [BigDnaMatchBlock *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_from_BigDnaMatchBlock(BigDnaMatchBlock * mat);
#define AlnRangeSet_from_BigDnaMatchBlock Wise2_AlnRangeSet_from_BigDnaMatchBlock


/* Function:  convert_PackAln_to_AlnBlock_BigDnaMatchBlock(pal)
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
AlnBlock * Wise2_convert_PackAln_to_AlnBlock_BigDnaMatchBlock(PackAln * pal);
#define convert_PackAln_to_AlnBlock_BigDnaMatchBlock Wise2_convert_PackAln_to_AlnBlock_BigDnaMatchBlock


/* Function:  PackAln_read_Expl_BigDnaMatchBlock(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [BigDnaMatchBlock *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Expl_BigDnaMatchBlock(BigDnaMatchBlock * mat);
#define PackAln_read_Expl_BigDnaMatchBlock Wise2_PackAln_read_Expl_BigDnaMatchBlock


/* Function:  PackAln_read_generic_BigDnaMatchBlock(mat,h)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [BigDnaMatchBlock *]
 * Arg:          h [UNKN ] Undocumented argument [BigDnaMatchBlock_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_generic_BigDnaMatchBlock(BigDnaMatchBlock * mat,BigDnaMatchBlock_access_func_holder h);
#define PackAln_read_generic_BigDnaMatchBlock Wise2_PackAln_read_generic_BigDnaMatchBlock


/* Function:  calculate_BigDnaMatchBlock(mat)
 *
 * Descrip:    This function calculates the BigDnaMatchBlock matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_BigDnaMatchBlock
 *
 *
 * Arg:        mat [UNKN ] BigDnaMatchBlock which contains explicit basematrix memory [BigDnaMatchBlock *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_BigDnaMatchBlock(BigDnaMatchBlock * mat);
#define calculate_BigDnaMatchBlock Wise2_calculate_BigDnaMatchBlock


/* Function:  calculate_dpenv_BigDnaMatchBlock(mat,dpenv)
 *
 * Descrip:    This function calculates the BigDnaMatchBlock matrix when in explicit mode, subject to the envelope
 *
 *
 * Arg:          mat [UNKN ] BigDnaMatchBlock which contains explicit basematrix memory [BigDnaMatchBlock *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_dpenv_BigDnaMatchBlock(BigDnaMatchBlock * mat,DPEnvelope * dpenv);
#define calculate_dpenv_BigDnaMatchBlock Wise2_calculate_dpenv_BigDnaMatchBlock


/* Function:  BigDnaMatchBlock_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [BigDnaMatchBlock *]
 *
 */
BigDnaMatchBlock * Wise2_BigDnaMatchBlock_alloc(void);
#define BigDnaMatchBlock_alloc Wise2_BigDnaMatchBlock_alloc


/* Function:  free_BigDnaMatchBlock(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [BigDnaMatchBlock *]
 *
 * Return [UNKN ]  Undocumented return value [BigDnaMatchBlock *]
 *
 */
BigDnaMatchBlock * Wise2_free_BigDnaMatchBlock(BigDnaMatchBlock * obj);
#define free_BigDnaMatchBlock Wise2_free_BigDnaMatchBlock


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
int Wise2_BigDnaMatchBlock_shatter_access_main(BigDnaMatchBlock * mat,int i,int j,int state);
#define BigDnaMatchBlock_shatter_access_main Wise2_BigDnaMatchBlock_shatter_access_main
int Wise2_BigDnaMatchBlock_shatter_access_special(BigDnaMatchBlock * mat,int i,int j,int state);
#define BigDnaMatchBlock_shatter_access_special Wise2_BigDnaMatchBlock_shatter_access_special
int Wise2_score_only_BigDnaMatchBlock(ComplexSequence* query,ComplexSequence* target ,DnaMatrix* comp55,DnaMatrix* comp65,DnaMatrix* comp75,DnaMatrix* comp85,DnaMatrix* comp95,Score g,Score g65,Score g55,Score u,Score v,Score s,Score b);
#define score_only_BigDnaMatchBlock Wise2_score_only_BigDnaMatchBlock
BigDnaMatchBlock * Wise2_allocate_BigDnaMatchBlock_only(ComplexSequence* query,ComplexSequence* target ,DnaMatrix* comp55,DnaMatrix* comp65,DnaMatrix* comp75,DnaMatrix* comp85,DnaMatrix* comp95,Score g,Score g65,Score g55,Score u,Score v,Score s,Score b);
#define allocate_BigDnaMatchBlock_only Wise2_allocate_BigDnaMatchBlock_only
void Wise2_init_BigDnaMatchBlock(BigDnaMatchBlock * mat);
#define init_BigDnaMatchBlock Wise2_init_BigDnaMatchBlock
AlnRange * Wise2_AlnRange_build_BigDnaMatchBlock(BigDnaMatchBlock * mat,int stopj,int stopspecstate,int * startj,int * startspecstate);
#define AlnRange_build_BigDnaMatchBlock Wise2_AlnRange_build_BigDnaMatchBlock
boolean Wise2_read_hidden_BigDnaMatchBlock(BigDnaMatchBlock * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out);
#define read_hidden_BigDnaMatchBlock Wise2_read_hidden_BigDnaMatchBlock
int Wise2_max_hidden_BigDnaMatchBlock(BigDnaMatchBlock * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_hidden_BigDnaMatchBlock Wise2_max_hidden_BigDnaMatchBlock
boolean Wise2_read_special_strip_BigDnaMatchBlock(BigDnaMatchBlock * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out);
#define read_special_strip_BigDnaMatchBlock Wise2_read_special_strip_BigDnaMatchBlock
int Wise2_max_special_strip_BigDnaMatchBlock(BigDnaMatchBlock * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_special_strip_BigDnaMatchBlock Wise2_max_special_strip_BigDnaMatchBlock
int Wise2_max_matrix_to_special_BigDnaMatchBlock(BigDnaMatchBlock * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_matrix_to_special_BigDnaMatchBlock Wise2_max_matrix_to_special_BigDnaMatchBlock
void Wise2_calculate_hidden_BigDnaMatchBlock(BigDnaMatchBlock * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv);
#define calculate_hidden_BigDnaMatchBlock Wise2_calculate_hidden_BigDnaMatchBlock
void Wise2_init_hidden_BigDnaMatchBlock(BigDnaMatchBlock * mat,int starti,int startj,int stopi,int stopj);
#define init_hidden_BigDnaMatchBlock Wise2_init_hidden_BigDnaMatchBlock
boolean Wise2_full_dc_BigDnaMatchBlock(BigDnaMatchBlock * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv);
#define full_dc_BigDnaMatchBlock Wise2_full_dc_BigDnaMatchBlock
boolean Wise2_do_dc_single_pass_BigDnaMatchBlock(BigDnaMatchBlock * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done);
#define do_dc_single_pass_BigDnaMatchBlock Wise2_do_dc_single_pass_BigDnaMatchBlock
void Wise2_push_dc_at_merge_BigDnaMatchBlock(BigDnaMatchBlock * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv);
#define push_dc_at_merge_BigDnaMatchBlock Wise2_push_dc_at_merge_BigDnaMatchBlock
void Wise2_follow_on_dc_BigDnaMatchBlock(BigDnaMatchBlock * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define follow_on_dc_BigDnaMatchBlock Wise2_follow_on_dc_BigDnaMatchBlock
void Wise2_run_up_dc_BigDnaMatchBlock(BigDnaMatchBlock * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define run_up_dc_BigDnaMatchBlock Wise2_run_up_dc_BigDnaMatchBlock
void Wise2_init_dc_BigDnaMatchBlock(BigDnaMatchBlock * mat);
#define init_dc_BigDnaMatchBlock Wise2_init_dc_BigDnaMatchBlock
int Wise2_start_end_find_end_BigDnaMatchBlock(BigDnaMatchBlock * mat,int * endj);
#define start_end_find_end_BigDnaMatchBlock Wise2_start_end_find_end_BigDnaMatchBlock
boolean Wise2_dc_optimised_start_end_calc_BigDnaMatchBlock(BigDnaMatchBlock *mat,DPEnvelope * dpenv);
#define dc_optimised_start_end_calc_BigDnaMatchBlock Wise2_dc_optimised_start_end_calc_BigDnaMatchBlock
void Wise2_init_start_end_linear_BigDnaMatchBlock(BigDnaMatchBlock * mat);
#define init_start_end_linear_BigDnaMatchBlock Wise2_init_start_end_linear_BigDnaMatchBlock
AlnConvertSet * Wise2_AlnConvertSet_BigDnaMatchBlock(void);
#define AlnConvertSet_BigDnaMatchBlock Wise2_AlnConvertSet_BigDnaMatchBlock
int Wise2_BigDnaMatchBlock_explicit_access_main(BigDnaMatchBlock * mat,int i,int j,int state);
#define BigDnaMatchBlock_explicit_access_main Wise2_BigDnaMatchBlock_explicit_access_main
int Wise2_BigDnaMatchBlock_explicit_access_special(BigDnaMatchBlock * mat,int i,int j,int state);
#define BigDnaMatchBlock_explicit_access_special Wise2_BigDnaMatchBlock_explicit_access_special
int Wise2_find_end_BigDnaMatchBlock(BigDnaMatchBlock * mat,int * ri,int * rj,int * state,boolean * isspecial,BigDnaMatchBlock_access_func_holder h);
#define find_end_BigDnaMatchBlock Wise2_find_end_BigDnaMatchBlock
void Wise2_BigDnaMatchBlock_debug_show_matrix(BigDnaMatchBlock * mat,int starti,int stopi,int startj,int stopj,FILE * ofp);
#define BigDnaMatchBlock_debug_show_matrix Wise2_BigDnaMatchBlock_debug_show_matrix
int Wise2_max_calc_BigDnaMatchBlock(BigDnaMatchBlock * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,BigDnaMatchBlock_access_func_holder h);
#define max_calc_BigDnaMatchBlock Wise2_max_calc_BigDnaMatchBlock
int Wise2_max_calc_special_BigDnaMatchBlock(BigDnaMatchBlock * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,BigDnaMatchBlock_access_func_holder h);
#define max_calc_special_BigDnaMatchBlock Wise2_max_calc_special_BigDnaMatchBlock

#ifdef _cplusplus
}
#endif

#endif
