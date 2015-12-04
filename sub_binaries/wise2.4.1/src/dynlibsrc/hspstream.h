#ifndef DYNAMITEhspstreamHEADERFILE
#define DYNAMITEhspstreamHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "../dynlibsrc/hsp.h"
#include "sequencestream.h"




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  untyped_LinearHSPmanager_to_Stream(lm,ws)
 *
 * Descrip:    untyped linear hsp manager to stream
 *
 *
 * Arg:        lm [UNKN ] Undocumented argument [void *]
 * Arg:        ws [UNKN ] Undocumented argument [Wise2WriteStreamInterface *]
 *
 */
void Wise2_untyped_LinearHSPmanager_to_Stream(void * lm,Wise2WriteStreamInterface * ws);
#define untyped_LinearHSPmanager_to_Stream Wise2_untyped_LinearHSPmanager_to_Stream


/* Function:  typed_LinearHSPmanager_to_Stream(lm,ws)
 *
 * Descrip:    typed linear hsp manager to stream
 *
 *             writes out sequneces first, with query as first
 *             sequence and then each target in turn, and then all the
 *             hsps in turn, only writing down the target name
 *
 *
 * Arg:        lm [UNKN ] Undocumented argument [LinearHSPmanager *]
 * Arg:        ws [UNKN ] Undocumented argument [Wise2WriteStreamInterface *]
 *
 */
void Wise2_typed_LinearHSPmanager_to_Stream(LinearHSPmanager * lm,Wise2WriteStreamInterface * ws);
#define typed_LinearHSPmanager_to_Stream Wise2_typed_LinearHSPmanager_to_Stream


/* Function:  untyped_LinearHSPmanager_from_Stream(rs)
 *
 * Descrip:    untyped read
 *
 *
 * Arg:        rs [UNKN ] Undocumented argument [Wise2ReadStreamInterface *]
 *
 * Return [UNKN ]  Undocumented return value [void *]
 *
 */
void * Wise2_untyped_LinearHSPmanager_from_Stream(Wise2ReadStreamInterface * rs);
#define untyped_LinearHSPmanager_from_Stream Wise2_untyped_LinearHSPmanager_from_Stream


/* Function:  typed_LinearHSPmanager_from_Stream(rs)
 *
 * Descrip:    typed linear hsp manager from stream
 *
 *             assummes sequences coming in first, with query as first sequence
 *             and then hsps one after the other
 *
 *
 * Arg:        rs [UNKN ] Undocumented argument [Wise2ReadStreamInterface *]
 *
 * Return [UNKN ]  Undocumented return value [LinearHSPmanager *]
 *
 */
LinearHSPmanager * Wise2_typed_LinearHSPmanager_from_Stream(Wise2ReadStreamInterface * rs);
#define typed_LinearHSPmanager_from_Stream Wise2_typed_LinearHSPmanager_from_Stream


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/

#ifdef _cplusplus
}
#endif

#endif
