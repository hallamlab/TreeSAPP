#ifndef DYNAMITEeditdistdpHEADERFILE
#define DYNAMITEeditdistdpHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"



#define EDITMATCH(string1c,string2c,match,mismatch) (string1c == string2c ? match : mismatch)

struct Wise2_EditString {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * string;   
    int len;     
    } ;  
/* EditString defined */ 
#ifndef DYNAMITE_DEFINED_EditString
typedef struct Wise2_EditString Wise2_EditString;
#define EditString Wise2_EditString
#define DYNAMITE_DEFINED_EditString
#endif


struct Wise2_EditDistance {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    BaseMatrix * basematrix;     
    ShatterMatrix * shatter;     
    int leni;    
    int lenj;    
    EditString* q;   
    EditString* t;   
    int match;   
    int mismatch;    
    int gap;     
    } ;  
/* EditDistance defined */ 
#ifndef DYNAMITE_DEFINED_EditDistance
typedef struct Wise2_EditDistance Wise2_EditDistance;
#define EditDistance Wise2_EditDistance
#define DYNAMITE_DEFINED_EditDistance
#endif


#ifdef PTHREAD
struct thread_pool_holder_EditDistance {  
    EditString* q;  /* Static query data: never free'd */ 
    EditString* t;  /* Static target data: never free'd */ 
    int match;   
    int mismatch;    
    int gap;     
    pthread_mutex_t input_lock;  
    pthread_mutex_t output_lock;     
    Hscore * out;    
    pthread_t * pool;    
    int number_of_threads;   
    boolean search_has_ended;    
    DBSearchImpl * dbsi;     
    } ;  
#endif /* PTHREAD */
struct Wise2_EditDistance_access_func_holder {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int (*access_main)(EditDistance*,int,int,int);   
    int (*access_special)(EditDistance*,int,int,int);    
    } ;  
/* EditDistance_access_func_holder defined */ 
#ifndef DYNAMITE_DEFINED_EditDistance_access_func_holder
typedef struct Wise2_EditDistance_access_func_holder Wise2_EditDistance_access_func_holder;
#define EditDistance_access_func_holder Wise2_EditDistance_access_func_holder
#define DYNAMITE_DEFINED_EditDistance_access_func_holder
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_EditString(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [EditString *]
 *
 * Return [UNKN ]  Undocumented return value [EditString *]
 *
 */
EditString * Wise2_hard_link_EditString(EditString * obj);
#define hard_link_EditString Wise2_hard_link_EditString


/* Function:  EditString_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [EditString *]
 *
 */
EditString * Wise2_EditString_alloc(void);
#define EditString_alloc Wise2_EditString_alloc


/* Function:  free_EditString(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [EditString *]
 *
 * Return [UNKN ]  Undocumented return value [EditString *]
 *
 */
EditString * Wise2_free_EditString(EditString * obj);
#define free_EditString Wise2_free_EditString


/* Function:  PackAln_read_Shatter_EditDistance(mat)
 *
 * Descrip:    Reads off PackAln from shatter matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [EditDistance *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Shatter_EditDistance(EditDistance * mat);
#define PackAln_read_Shatter_EditDistance Wise2_PackAln_read_Shatter_EditDistance


/* Function:  calculate_shatter_EditDistance(mat,dpenv)
 *
 * Descrip:    This function calculates the EditDistance matrix when in shatter mode
 *
 *
 * Arg:          mat [UNKN ] (null) [EditDistance *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_shatter_EditDistance(EditDistance * mat,DPEnvelope * dpenv);
#define calculate_shatter_EditDistance Wise2_calculate_shatter_EditDistance


/* Function:  search_EditDistance(dbsi,out,q,t,match,mismatch,gap)
 *
 * Descrip:    This function makes a database search of EditDistance
 *             It uses the dbsi structure to choose which implementation
 *             to use of the database searching. This way at run time you
 *             can switch between single threaded/multi-threaded or hardware
 *
 *
 * Arg:            dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 * Arg:             out [UNKN ] Undocumented argument [Hscore *]
 * Arg:               q [UNKN ] Undocumented argument [EditString*]
 * Arg:               t [UNKN ] Undocumented argument [EditString*]
 * Arg:           match [UNKN ] Undocumented argument [int]
 * Arg:        mismatch [UNKN ] Undocumented argument [int]
 * Arg:             gap [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_search_EditDistance(DBSearchImpl * dbsi,Hscore * out,EditString* q,EditString* t ,int match,int mismatch,int gap);
#define search_EditDistance Wise2_search_EditDistance


/* Function:  serial_search_EditDistance(out,q,t,match,mismatch,gap)
 *
 * Descrip:    This function makes a database search of EditDistance
 *             It is a single processor implementation
 *
 *
 * Arg:             out [UNKN ] Undocumented argument [Hscore *]
 * Arg:               q [UNKN ] Undocumented argument [EditString*]
 * Arg:               t [UNKN ] Undocumented argument [EditString*]
 * Arg:           match [UNKN ] Undocumented argument [int]
 * Arg:        mismatch [UNKN ] Undocumented argument [int]
 * Arg:             gap [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_serial_search_EditDistance(Hscore * out,EditString* q,EditString* t ,int match,int mismatch,int gap);
#define serial_search_EditDistance Wise2_serial_search_EditDistance


/* Function:  PackAln_bestmemory_EditDistance(q,t,match,mismatch,gap,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_EditDistance
 *
 *
 * Arg:               q [UNKN ] query data structure [EditString*]
 * Arg:               t [UNKN ] target data structure [EditString*]
 * Arg:           match [UNKN ] Resource [int]
 * Arg:        mismatch [UNKN ] Resource [int]
 * Arg:             gap [UNKN ] Resource [int]
 * Arg:           dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:            dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_bestmemory_EditDistance(EditString* q,EditString* t ,int match,int mismatch,int gap,DPEnvelope * dpenv,DPRunImpl * dpri);
#define PackAln_bestmemory_EditDistance Wise2_PackAln_bestmemory_EditDistance


/* Function:  allocate_Expl_EditDistance(q,t,match,mismatch,gap,dpri)
 *
 * Descrip:    This function allocates the EditDistance structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_EditDistance_only
 *
 *
 * Arg:               q [UNKN ] query data structure [EditString*]
 * Arg:               t [UNKN ] target data structure [EditString*]
 * Arg:           match [UNKN ] Resource [int]
 * Arg:        mismatch [UNKN ] Resource [int]
 * Arg:             gap [UNKN ] Resource [int]
 * Arg:            dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [EditDistance *]
 *
 */
EditDistance * Wise2_allocate_Expl_EditDistance(EditString* q,EditString* t ,int match,int mismatch,int gap,DPRunImpl * dpri);
#define allocate_Expl_EditDistance Wise2_allocate_Expl_EditDistance


/* Function:  recalculate_PackAln_EditDistance(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by EditDistance
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [EditDistance *]
 *
 */
void Wise2_recalculate_PackAln_EditDistance(PackAln * pal,EditDistance * mat);
#define recalculate_PackAln_EditDistance Wise2_recalculate_PackAln_EditDistance


/* Function:  allocate_Small_EditDistance(q,t,match,mismatch,gap)
 *
 * Descrip:    This function allocates the EditDistance structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_EditDistance_only
 *
 *
 * Arg:               q [UNKN ] query data structure [EditString*]
 * Arg:               t [UNKN ] target data structure [EditString*]
 * Arg:           match [UNKN ] Resource [int]
 * Arg:        mismatch [UNKN ] Resource [int]
 * Arg:             gap [UNKN ] Resource [int]
 *
 * Return [UNKN ]  Undocumented return value [EditDistance *]
 *
 */
EditDistance * Wise2_allocate_Small_EditDistance(EditString* q,EditString* t ,int match,int mismatch,int gap);
#define allocate_Small_EditDistance Wise2_allocate_Small_EditDistance


/* Function:  PackAln_calculate_Small_EditDistance(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for EditDistance structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_EditDistance 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_EditDistance 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [EditDistance *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_calculate_Small_EditDistance(EditDistance * mat,DPEnvelope * dpenv);
#define PackAln_calculate_Small_EditDistance Wise2_PackAln_calculate_Small_EditDistance


/* Function:  AlnRangeSet_calculate_Small_EditDistance(mat)
 *
 * Descrip:    This function calculates an alignment for EditDistance structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_EditDistance 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_EditDistance
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_EditDistance 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [EditDistance *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_calculate_Small_EditDistance(EditDistance * mat);
#define AlnRangeSet_calculate_Small_EditDistance Wise2_AlnRangeSet_calculate_Small_EditDistance


/* Function:  AlnRangeSet_from_EditDistance(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for EditDistance structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_EditDistance 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_EditDistance
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [EditDistance *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_from_EditDistance(EditDistance * mat);
#define AlnRangeSet_from_EditDistance Wise2_AlnRangeSet_from_EditDistance


/* Function:  convert_PackAln_to_AlnBlock_EditDistance(pal)
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
AlnBlock * Wise2_convert_PackAln_to_AlnBlock_EditDistance(PackAln * pal);
#define convert_PackAln_to_AlnBlock_EditDistance Wise2_convert_PackAln_to_AlnBlock_EditDistance


/* Function:  PackAln_read_Expl_EditDistance(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [EditDistance *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Expl_EditDistance(EditDistance * mat);
#define PackAln_read_Expl_EditDistance Wise2_PackAln_read_Expl_EditDistance


/* Function:  PackAln_read_generic_EditDistance(mat,h)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [EditDistance *]
 * Arg:          h [UNKN ] Undocumented argument [EditDistance_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_generic_EditDistance(EditDistance * mat,EditDistance_access_func_holder h);
#define PackAln_read_generic_EditDistance Wise2_PackAln_read_generic_EditDistance


/* Function:  calculate_EditDistance(mat)
 *
 * Descrip:    This function calculates the EditDistance matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_EditDistance
 *
 *
 * Arg:        mat [UNKN ] EditDistance which contains explicit basematrix memory [EditDistance *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_EditDistance(EditDistance * mat);
#define calculate_EditDistance Wise2_calculate_EditDistance


/* Function:  calculate_dpenv_EditDistance(mat,dpenv)
 *
 * Descrip:    This function calculates the EditDistance matrix when in explicit mode, subject to the envelope
 *
 *
 * Arg:          mat [UNKN ] EditDistance which contains explicit basematrix memory [EditDistance *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_dpenv_EditDistance(EditDistance * mat,DPEnvelope * dpenv);
#define calculate_dpenv_EditDistance Wise2_calculate_dpenv_EditDistance


/* Function:  EditDistance_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [EditDistance *]
 *
 */
EditDistance * Wise2_EditDistance_alloc(void);
#define EditDistance_alloc Wise2_EditDistance_alloc


/* Function:  free_EditDistance(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [EditDistance *]
 *
 * Return [UNKN ]  Undocumented return value [EditDistance *]
 *
 */
EditDistance * Wise2_free_EditDistance(EditDistance * obj);
#define free_EditDistance Wise2_free_EditDistance


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
EditString * Wise2_new_EditString(char * string);
#define new_EditString Wise2_new_EditString


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
int Wise2_EditDistance_shatter_access_main(EditDistance * mat,int i,int j,int state);
#define EditDistance_shatter_access_main Wise2_EditDistance_shatter_access_main
int Wise2_EditDistance_shatter_access_special(EditDistance * mat,int i,int j,int state);
#define EditDistance_shatter_access_special Wise2_EditDistance_shatter_access_special
void * Wise2_thread_loop_EditDistance(void * ptr);
#define thread_loop_EditDistance Wise2_thread_loop_EditDistance
int Wise2_score_only_EditDistance(EditString* q,EditString* t ,int match,int mismatch,int gap);
#define score_only_EditDistance Wise2_score_only_EditDistance
EditDistance * Wise2_allocate_EditDistance_only(EditString* q,EditString* t ,int match,int mismatch,int gap);
#define allocate_EditDistance_only Wise2_allocate_EditDistance_only
void Wise2_init_EditDistance(EditDistance * mat);
#define init_EditDistance Wise2_init_EditDistance
AlnRange * Wise2_AlnRange_build_EditDistance(EditDistance * mat,int stopj,int stopspecstate,int * startj,int * startspecstate);
#define AlnRange_build_EditDistance Wise2_AlnRange_build_EditDistance
boolean Wise2_read_hidden_EditDistance(EditDistance * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out);
#define read_hidden_EditDistance Wise2_read_hidden_EditDistance
int Wise2_max_hidden_EditDistance(EditDistance * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_hidden_EditDistance Wise2_max_hidden_EditDistance
boolean Wise2_read_special_strip_EditDistance(EditDistance * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out);
#define read_special_strip_EditDistance Wise2_read_special_strip_EditDistance
int Wise2_max_special_strip_EditDistance(EditDistance * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_special_strip_EditDistance Wise2_max_special_strip_EditDistance
int Wise2_max_matrix_to_special_EditDistance(EditDistance * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_matrix_to_special_EditDistance Wise2_max_matrix_to_special_EditDistance
void Wise2_calculate_hidden_EditDistance(EditDistance * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv);
#define calculate_hidden_EditDistance Wise2_calculate_hidden_EditDistance
void Wise2_init_hidden_EditDistance(EditDistance * mat,int starti,int startj,int stopi,int stopj);
#define init_hidden_EditDistance Wise2_init_hidden_EditDistance
boolean Wise2_full_dc_EditDistance(EditDistance * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv);
#define full_dc_EditDistance Wise2_full_dc_EditDistance
boolean Wise2_do_dc_single_pass_EditDistance(EditDistance * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done);
#define do_dc_single_pass_EditDistance Wise2_do_dc_single_pass_EditDistance
void Wise2_push_dc_at_merge_EditDistance(EditDistance * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv);
#define push_dc_at_merge_EditDistance Wise2_push_dc_at_merge_EditDistance
void Wise2_follow_on_dc_EditDistance(EditDistance * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define follow_on_dc_EditDistance Wise2_follow_on_dc_EditDistance
void Wise2_run_up_dc_EditDistance(EditDistance * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define run_up_dc_EditDistance Wise2_run_up_dc_EditDistance
void Wise2_init_dc_EditDistance(EditDistance * mat);
#define init_dc_EditDistance Wise2_init_dc_EditDistance
int Wise2_start_end_find_end_EditDistance(EditDistance * mat,int * endj);
#define start_end_find_end_EditDistance Wise2_start_end_find_end_EditDistance
boolean Wise2_dc_optimised_start_end_calc_EditDistance(EditDistance *mat,DPEnvelope * dpenv);
#define dc_optimised_start_end_calc_EditDistance Wise2_dc_optimised_start_end_calc_EditDistance
void Wise2_init_start_end_linear_EditDistance(EditDistance * mat);
#define init_start_end_linear_EditDistance Wise2_init_start_end_linear_EditDistance
AlnConvertSet * Wise2_AlnConvertSet_EditDistance(void);
#define AlnConvertSet_EditDistance Wise2_AlnConvertSet_EditDistance
int Wise2_EditDistance_explicit_access_main(EditDistance * mat,int i,int j,int state);
#define EditDistance_explicit_access_main Wise2_EditDistance_explicit_access_main
int Wise2_EditDistance_explicit_access_special(EditDistance * mat,int i,int j,int state);
#define EditDistance_explicit_access_special Wise2_EditDistance_explicit_access_special
int Wise2_find_end_EditDistance(EditDistance * mat,int * ri,int * rj,int * state,boolean * isspecial,EditDistance_access_func_holder h);
#define find_end_EditDistance Wise2_find_end_EditDistance
void Wise2_EditDistance_debug_show_matrix(EditDistance * mat,int starti,int stopi,int startj,int stopj,FILE * ofp);
#define EditDistance_debug_show_matrix Wise2_EditDistance_debug_show_matrix
int Wise2_max_calc_EditDistance(EditDistance * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,EditDistance_access_func_holder h);
#define max_calc_EditDistance Wise2_max_calc_EditDistance
int Wise2_max_calc_special_EditDistance(EditDistance * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,EditDistance_access_func_holder h);
#define max_calc_special_EditDistance Wise2_max_calc_special_EditDistance

#ifdef _cplusplus
}
#endif

#endif
