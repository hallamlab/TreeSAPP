#ifndef DYNAMITEdynashadowHEADERFILE
#define DYNAMITEdynashadowHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna2.h"







    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  alloc_small_func_GenericMatrix(dfp,gm)
 *
 * Descrip:    makes the allocate_Small_xxx function.
 *             This calls allocate_xxx_only function
 *             (made by /write_safe_alloc_function)
 *             and then allocates basematrix stuff as well.
 *
 *
 * Arg:        dfp [UNKN ] dynamite file pointer [DYNFILE *]
 * Arg:         gm [READ ] generic matrix structure [const GenericMatrix *]
 *
 */
void alloc_small_func_GenericMatrix(DYNFILE * dfp,const GenericMatrix * gm);


/* Function:  make_small_calculate_func(dfp,gm)
 *
 * Descrip:    make the calculate function for
 *             small PackAln system
 *
 *
 * Arg:        dfp [UNKN ] Undocumented argument [DYNFILE *]
 * Arg:         gm [UNKN ] Undocumented argument [GenericMatrix *]
 *
 */
void make_small_calculate_func(DYNFILE * dfp,GenericMatrix * gm);


/* Function:  one_shot_AlnRangeSet_func(dfp,gm)
 *
 * Descrip:    makes AlnRangeSet from small memory
 *
 *
 * Arg:        dfp [UNKN ] Undocumented argument [DYNFILE *]
 * Arg:         gm [UNKN ] Undocumented argument [GenericMatrix *]
 *
 */
void one_shot_AlnRangeSet_func(DYNFILE * dfp,GenericMatrix * gm);


/* Function:  write_dc_PackAln_build_func(dfp,gm)
 *
 * Descrip:    This functions is now defunct
 *
 *
 * Arg:        dfp [UNKN ] Undocumented argument [DYNFILE *]
 * Arg:         gm [UNKN ] Undocumented argument [GenericMatrix *]
 *
 */
void write_dc_PackAln_build_func(DYNFILE * dfp,GenericMatrix * gm);


/* Function:  write_full_dc_func(dfp,gm)
 *
 * Descrip:    writes the main divide and conquor routine
 *
 *
 * Arg:        dfp [UNKN ] Undocumented argument [DYNFILE *]
 * Arg:         gm [UNKN ] Undocumented argument [GenericMatrix *]
 *
 */
void write_full_dc_func(DYNFILE * dfp,GenericMatrix * gm);


/* Function:  write_shadow_dc_alloc(dfp,gm)
 *
 * Descrip:    Defunct (I think)
 *
 *
 * Arg:        dfp [UNKN ] Undocumented argument [DYNFILE *]
 * Arg:         gm [UNKN ] Undocumented argument [GenericMatrix *]
 *
 */
void write_shadow_dc_alloc(DYNFILE * dfp,GenericMatrix * gm);


/* Function:  heavy_optimised_shadow_GenericMatrix(dfp,gm)
 *
 * Descrip:    heavily optimised start end calc function
 *
 *
 * Arg:        dfp [UNKN ] Undocumented argument [DYNFILE *]
 * Arg:         gm [UNKN ] Undocumented argument [GenericMatrix *]
 *
 */
void heavy_optimised_shadow_GenericMatrix(DYNFILE * dfp,GenericMatrix * gm);


/* Function:  optimised_shadow_GenericMatrix(dfp,*gm)
 *
 * Descrip:    Makes the optimised shadow matrix routine,
 *             worked out by Steve Searle - memory access
 *             is put into one array so that the routine
 *             is faster
 *
 *
 * Arg:        dfp [UNKN ] Undocumented argument [DYNFILE *]
 * Arg:        *gm [UNKN ] Undocumented argument [const GenericMatrix]
 *
 */
void optimised_shadow_GenericMatrix(DYNFILE * dfp,const GenericMatrix *gm);


/* Function:  write_main_shadow_block(dfp,gm,matrixtag,pointer_tag,specialtag,shadow_main_tag,shadow_special_tag,shadow_length,shadow_on_special,use_special,debug,use_shadow_pointer)
 *
 * Descrip:    The core inner loop for shadow methods. Pretty terrifying stuff in here.
 *
 *             Shadow is considered to be in the memory shadow_main_tag and shadow_special_tag
 *             and usually has the form MATRIX_TYPE_MAIN/SPECIAL_SP. 
 *
 *             shadow has positions 0-length-1, as defined by shadow_length. These are
 *             indexed by k.
 *
 *             Either shadow_on_special is false
 *             Some other routine has to place the shadow pointers. This routine just
 *             propagates the shadow pointers around
 *
 *             Or if shadow_on_special  is true
 *             This routines pushes when it leaves special and also pushes when it 
 *             enters special. This way, special states have the necessary information
 *             to know where in the main matrix they were made from. Means there must be 7 shadow
 *             positions.
 *
 *             In addition, as a complete mess, shadow on special needs to have score independent
 *             calcs added onto each movement, otherwise the scores are incorrect. So... this
 *             is not *IDEAL* in the slightest!!!
 *
 *
 * Arg:                       dfp [UNKN ] Undocumented argument [DYNFILE *]
 * Arg:                        gm [UNKN ] Undocumented argument [GenericMatrix *]
 * Arg:                 matrixtag [UNKN ] Undocumented argument [char *]
 * Arg:               pointer_tag [UNKN ] Undocumented argument [char *]
 * Arg:                specialtag [UNKN ] Undocumented argument [char *]
 * Arg:           shadow_main_tag [UNKN ] Undocumented argument [char *]
 * Arg:        shadow_special_tag [UNKN ] Undocumented argument [char *]
 * Arg:             shadow_length [UNKN ] Undocumented argument [int]
 * Arg:         shadow_on_special [UNKN ] Undocumented argument [boolean]
 * Arg:               use_special [UNKN ] Undocumented argument [boolean]
 * Arg:                     debug [UNKN ] Undocumented argument [int]
 * Arg:        use_shadow_pointer [UNKN ] Undocumented argument [int]
 *
 */
void write_main_shadow_block(DYNFILE * dfp,GenericMatrix * gm,char * matrixtag,char * pointer_tag,char * specialtag,char * shadow_main_tag,char * shadow_special_tag,int shadow_length,boolean shadow_on_special,boolean use_special,int debug,int use_shadow_pointer);


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
void write_dc_functions(DYNFILE * dfp,GenericMatrix * gm);
void write_AlnRangeSet_build_func(DYNFILE * dfp,GenericMatrix * gm);
void write_AlnRange_build_func(DYNFILE * dfp,GenericMatrix * gm);
void write_single_dc_pass_func(DYNFILE * dfp,GenericMatrix * gm);
void write_push_dc_func(DYNFILE * dfp,GenericMatrix * gm);
void write_init_dc_func(DYNFILE * dfp,GenericMatrix * gm);
void write_up_to_dc_func(DYNFILE * dfp,GenericMatrix * gm);
void write_follow_on_dc_func(DYNFILE * dfp,GenericMatrix * gm);
void write_special_strip_read_func(DYNFILE * dfp,GenericMatrix * gm);
void write_hidden_read_func(DYNFILE * dfp,GenericMatrix * gm);
void write_hidden_calc_func(DYNFILE * dfp,GenericMatrix * gm);
void write_hidden_init_func(DYNFILE * dfp,GenericMatrix * gm);
void write_hidden_max_func(DYNFILE * dfp,GenericMatrix * gm);
void write_matrix_to_special_max_func(DYNFILE * dfp,GenericMatrix * gm);
void write_special_strip_max_func(DYNFILE * dfp,GenericMatrix * gm);
void write_shadow_dc_macros(DYNFILE * dfp,GenericMatrix * gm);
void write_start_end_init(DYNFILE * dfp,GenericMatrix * gm);
void write_start_end_find_end(DYNFILE * dfp,GenericMatrix * gm);
void write_start_end_build(DYNFILE * dfp,GenericMatrix * gm);
void write_start_end_macros(DYNFILE * dfp,GenericMatrix * gm);
void write_shadow_start_end_alloc(DYNFILE * dfp,GenericMatrix * gm);
void write_special_shadow_block(DYNFILE * dfp,GenericMatrix * gm,char * matrix,char * pointer_tag,char * special,char *  shadow_main_tag,char * shadow_special_tag,int shadow_length,boolean shadow_on_special,int debug);


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/

#ifdef _cplusplus
}
#endif

#endif
