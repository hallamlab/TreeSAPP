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


/* Function:  show_pretty_Seq_dba_align_btcanvas(alb,one,two,btc)
 *
 * Descrip:    Different func signature for show_pretty
 *
 *
 * Arg:        alb [UNKN ] Undocumented argument [AlnBlock *]
 * Arg:        one [UNKN ] Undocumented argument [Sequence *]
 * Arg:        two [UNKN ] Undocumented argument [Sequence *]
 * Arg:        btc [UNKN ] Undocumented argument [btCanvas *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_show_pretty_Seq_dba_align_btcanvas(AlnBlock * alb,Sequence * one,Sequence * two,btCanvas * btc);
#define show_pretty_Seq_dba_align_btcanvas Wise2_show_pretty_Seq_dba_align_btcanvas


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
boolean Wise2_is_unmatched_block(AlnColumn * alc);
#define is_unmatched_block Wise2_is_unmatched_block
boolean Wise2_show_pretty_dba_align_btcanvas(btCanvas * btc,AlnBlock * alb,Sequence * one,Sequence * two,boolean blast_compatible);
#define show_pretty_dba_align_btcanvas Wise2_show_pretty_dba_align_btcanvas

#ifdef _cplusplus
}
#endif

#endif
