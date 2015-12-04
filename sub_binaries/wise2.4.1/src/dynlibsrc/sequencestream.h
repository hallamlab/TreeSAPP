#ifndef DYNAMITEsequencestreamHEADERFILE
#define DYNAMITEsequencestreamHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "sequence.h"



    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  untyped_SequenceSet_from_Stream(rs)
 *
 * Descrip:    untyped version of reading from Wise2ReadStreamInterface
 *
 *
 * Arg:        rs [UNKN ] Undocumented argument [Wise2ReadStreamInterface *]
 *
 * Return [UNKN ]  Undocumented return value [void *]
 *
 */
void * Wise2_untyped_SequenceSet_from_Stream(Wise2ReadStreamInterface * rs);
#define untyped_SequenceSet_from_Stream Wise2_untyped_SequenceSet_from_Stream


/* Function:  typed_SequenceSet_from_Stream(rs)
 *
 * Descrip:    Typed version of reading from Stream, making a SequenceSet
 *
 *
 * Arg:        rs [UNKN ] Undocumented argument [Wise2ReadStreamInterface *]
 *
 * Return [UNKN ]  Undocumented return value [SequenceSet *]
 *
 */
SequenceSet * Wise2_typed_SequenceSet_from_Stream(Wise2ReadStreamInterface * rs);
#define typed_SequenceSet_from_Stream Wise2_typed_SequenceSet_from_Stream


/* Function:  typed_write_SequenceSet_to_Stream(set,ws)
 *
 * Descrip:    typed version of writing a sequence set to a stream
 *
 *
 * Arg:        set [UNKN ] Undocumented argument [SequenceSet *]
 * Arg:         ws [UNKN ] Undocumented argument [Wise2WriteStreamInterface *]
 *
 */
void Wise2_typed_write_SequenceSet_to_Stream(SequenceSet * set,Wise2WriteStreamInterface * ws);
#define typed_write_SequenceSet_to_Stream Wise2_typed_write_SequenceSet_to_Stream


/* Function:  untyped_write_SequenceSet_to_Stream(s,ws)
 *
 * Descrip:    untyped version for writing a sequence set to a stream
 *
 *
 * Arg:         s [UNKN ] Undocumented argument [void *]
 * Arg:        ws [UNKN ] Undocumented argument [Wise2WriteStreamInterface *]
 *
 */
void Wise2_untyped_write_SequenceSet_to_Stream(void * s,Wise2WriteStreamInterface * ws);
#define untyped_write_SequenceSet_to_Stream Wise2_untyped_write_SequenceSet_to_Stream


/* Function:  typed_write_one_Sequence_to_Stream(seq,ws)
 *
 * Descrip:    internal function for writing out one sequence
 *
 *
 * Arg:        seq [UNKN ] Undocumented argument [Sequence *]
 * Arg:         ws [UNKN ] Undocumented argument [Wise2WriteStreamInterface *]
 *
 */
void Wise2_typed_write_one_Sequence_to_Stream(Sequence * seq,Wise2WriteStreamInterface * ws);
#define typed_write_one_Sequence_to_Stream Wise2_typed_write_one_Sequence_to_Stream


/* Function:  untyped_write_Sequence_to_Stream(seq,ws)
 *
 * Descrip:    untyped write one sequence stream
 *
 *
 * Arg:        seq [UNKN ] Undocumented argument [void *]
 * Arg:         ws [UNKN ] Undocumented argument [Wise2WriteStreamInterface *]
 *
 */
void Wise2_untyped_write_Sequence_to_Stream(void * seq,Wise2WriteStreamInterface * ws);
#define untyped_write_Sequence_to_Stream Wise2_untyped_write_Sequence_to_Stream


/* Function:  untyped_read_Sequence_from_Stream(rs)
 *
 * Descrip:    reading one sequence from stream
 *
 *
 * Arg:        rs [UNKN ] Undocumented argument [Wise2ReadStreamInterface *]
 *
 * Return [UNKN ]  Undocumented return value [void *]
 *
 */
void * Wise2_untyped_read_Sequence_from_Stream(Wise2ReadStreamInterface * rs);
#define untyped_read_Sequence_from_Stream Wise2_untyped_read_Sequence_from_Stream


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
