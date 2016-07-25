#ifndef DYNAMITEdbadisplayHEADERFILE
#define DYNAMITEdbadisplayHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"





    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  show_pretty_dba_align(alb,one,two,ofp)
 *
 * Descrip:    Shows an alignment of from the dba algorithm in
 *             pretty formatted ascii text
 *
 *
 * Arg:        alb [UNKN ] Undocumented argument [AlnBlock *]
 * Arg:        one [UNKN ] Undocumented argument [Sequence *]
 * Arg:        two [UNKN ] Undocumented argument [Sequence *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_show_pretty_dba_align(AlnBlock * alb,Sequence * one,Sequence * two,FILE * ofp);
#define show_pretty_dba_align Wise2_show_pretty_dba_align


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
boolean Wise2_is_unmatched_block(AlnColumn * alc);
#define is_unmatched_block Wise2_is_unmatched_block
boolean Wise2_show_pretty_dba_align_btcanvas(btCanvas * btc,AlnBlock * alb,Sequence * one,Sequence * two);
#define show_pretty_dba_align_btcanvas Wise2_show_pretty_dba_align_btcanvas

#ifdef _cplusplus
}
#endif

#endif
