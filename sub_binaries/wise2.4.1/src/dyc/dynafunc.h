#ifndef DYNAMITEdynafuncHEADERFILE
#define DYNAMITEdynafuncHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif

#include "dynfile.h"
#include "dyna2.h"
#include "labelmaster.h"
#include "linesubs.h"
#include "dynashadow.h"
#include "dynadb.h"
#include "dpimpl.h"
#include "dbthread.h"
#include "probal.h"
#include "dynadebug.h"



    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  recalculate_PackAln_func(dfp,gm)
 *
 * Descrip:    makes the recalculate packaln function
 *
 *
 * Arg:        dfp [UNKN ] Undocumented argument [DYNFILE *]
 * Arg:         gm [UNKN ] Undocumented argument [const GenericMatrix *]
 *
 */
void recalculate_PackAln_func(DYNFILE * dfp,const GenericMatrix * gm);


/* Function:  alloc_expl_func_GenericMatrix(dfp,gm)
 *
 * Descrip:    makes the allocate_Expl_xxx function.
 *             This calls allocate_xxx_only function
 *             (made by /write_safe_alloc_function)
 *             and then allocates basematrix stuff as well.
 *
 *
 * Arg:        dfp [UNKN ] dynamite file pointer [DYNFILE *]
 * Arg:         gm [READ ] generic matrix structure [const GenericMatrix *]
 *
 */
void alloc_expl_func_GenericMatrix(DYNFILE * dfp,const GenericMatrix * gm);


/* Function:  write_safe_alloc_function(dfp,gm)
 *
 * Descrip:    produces the allocate_xxx_only function,
 *             which allocates the matrix structure, checks
 *             resources which it can check, but does NOT 
 *             allocate basematrix area
 *
 *             This function will be called by allocate_Expl_xxxx 
 *             and allocate_Small_xxxx etc.
 *
 *
 * Arg:        dfp [UNKN ] dynmaite file pointer to func/head [DYNFILE *]
 * Arg:         gm [UNKN ] generic matrix structure [GenericMatrix *]
 *
 */
void write_safe_alloc_function(DYNFILE * dfp,GenericMatrix * gm);


/* Function:  add_args_GenericMatrix_FuncInfo(fi,gm)
 *
 * Descrip:    Information partner to /get_argstr_GenericMatrix.
 *             Loads up the arglist information for the GenericMatrix,
 *             ie query-type query-name "what information" etc.
 *
 *
 *
 * Arg:        fi [WRITE] FuncInfo structure to add ArgInfo to  [FuncInfo *]
 * Arg:        gm [READ ] generic matrix structure to read from [const GenericMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean add_args_GenericMatrix_FuncInfo(FuncInfo * fi,const GenericMatrix * gm);


/* Function:  get_chainstr_GenericMatrix(gm)
 *
 * Descrip:    makes an the argument calling string which
 *             is compatible with the arg_str from
 *             /get_argstr_GenericMatrix
 *
 *             eg "query,target,comp_mat"
 *
 *
 * Arg:        gm [READ ] structure holding generic matrix [const GenericMatrix *]
 *
 * Return [UNKN ]  allocated string (must free) with chained-args [char *]
 *
 */
char * get_chainstr_GenericMatrix(const GenericMatrix * gm);


/* Function:  get_argstr_GenericMatrix(gm)
 *
 * Descrip:    makes the argument list for a generic matrix alloc, ie
 *             query-type query-name,target-type target-name etc
 *
 *             To load up info, use /add_arg_GenericMatrix_FuncInfo
 *             To chain up along function calls, use /get_chainstr_GenericMatrix
 *
 *
 * Arg:        gm [READ ] structure holding generic matrix [const GenericMatrix *]
 *
 * Return [UNKN ]  allocated string (must free) with args [char *]
 *
 */
char * get_argstr_GenericMatrix(const GenericMatrix * gm);


/* Function:  init_matrix_func(dfp,gm)
 *
 * Descrip:    init explicit matrix
 *
 *
 * Arg:        dfp [UNKN ] Undocumented argument [DYNFILE *]
 * Arg:         gm [UNKN ] Undocumented argument [GenericMatrix *]
 *
 */
void init_matrix_func(DYNFILE * dfp,GenericMatrix * gm);


/* Function:  matrix_calculate_func(dfp,gm)
 *
 * Descrip:    makes calculate_xxx functions, which
 *             is for explicit matrix implementations
 *
 *
 * Arg:        dfp [UNKN ] Undocumented argument [DYNFILE *]
 * Arg:         gm [UNKN ] Undocumented argument [GenericMatrix *]
 *
 */
void matrix_calculate_func(DYNFILE * dfp,GenericMatrix * gm);


/* Function:  matrix_calculate_func_dpenv(dfp,gm)
 *
 * Descrip:    makes calculate_xxx functions, which
 *             is for explicit matrix implementations
 *
 *
 * Arg:        dfp [UNKN ] Undocumented argument [DYNFILE *]
 * Arg:         gm [UNKN ] Undocumented argument [GenericMatrix *]
 *
 */
void matrix_calculate_func_dpenv(DYNFILE * dfp,GenericMatrix * gm);


/* Function:  do_cell_function(dfp,gm)
 *
 * Descrip:    useful function for debugging, but not
 *             currently used
 *
 *
 * Arg:        dfp [UNKN ] Undocumented argument [DYNFILE *]
 * Arg:         gm [UNKN ] Undocumented argument [GenericMatrix *]
 *
 */
void do_cell_function(DYNFILE * dfp,GenericMatrix * gm);


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
void one_shot_aln_func(DYNFILE * dfp,GenericMatrix * gm,DPImplementation * dpi);
void write_dotter_dump(DYNFILE * dfp,GenericMatrix * gm);
void write_memory_macros(DYNFILE * dfp,GenericMatrix * gm);
void write_search_macros(DYNFILE * dfp,const GenericMatrix * gm,DPImplementation * dpi);
void write_expl_access_funcs(DYNFILE * dfp,GenericMatrix * gm);
void write_expl_read_func(DYNFILE * dfp,GenericMatrix * gm);
void write_basic_read_func(DYNFILE * dfp,GenericMatrix * gm);
void write_pal_to_ars_func(DYNFILE * dfp,GenericMatrix * gm);
void write_alncconvert_make_func(DYNFILE * dfp,GenericMatrix * gm);
void write_aln_conversion_func(DYNFILE * dfp,GenericMatrix * gm);
void find_end_func(DYNFILE * dfp,GenericMatrix * gm);
void debug_func(DYNFILE * dfp,GenericMatrix * gm);
void write_special_max_calc_func_debug(DYNFILE * dfp,GenericMatrix * gm,int debug);
void write_special_max_calc_func(DYNFILE * dfp,GenericMatrix * gm);
void write_max_calc_block(DYNFILE * dfp,GenericMatrix * gm,char * matrix_tag,char * special_tag,boolean use_special,boolean use_holder);
void write_max_calc_func(DYNFILE * dfp,GenericMatrix * gm);
char * source_allowed_statement(int position,int offi,int offj);
void write_score_block_debug(DYNFILE * dfp,GenericMatrix * gm,char * matrixtag,char * pointertag,char * specialtag,boolean use_special,int debug);
void write_score_block(DYNFILE * dfp,GenericMatrix * gm,char * matrixtag,char * pointertag,char * specialtag,boolean use_special);
void do_special_function(DYNFILE * dfp,GenericMatrix * gm);
void write_search_distributor_func(DYNFILE * dfp,GenericMatrix * gm,DPImplementation * dpi);
void write_GenericMatrix_header(DYNFILE * dfp,GenericMatrix * gm,DPImplementation * dpi);
void write_GenericMatrix_func(DYNFILE * dfp,GenericMatrix * gm,DPImplementation * dpi);


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void write_special_block(DYNFILE * dfp,GenericMatrix * gm,char * matrix,char * special,char * bestscore);

#ifdef _cplusplus
}
#endif

#endif
