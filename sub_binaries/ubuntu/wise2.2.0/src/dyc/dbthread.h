#ifndef DYNAMITEdbthreadHEADERFILE
#define DYNAMITEdbthreadHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna2.h"
#include "dpimpl.h"




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  write_thread_struct(dfp,dpi,gm)
 *
 * Descrip:    writes out thread structure
 *             for the use in the loop code
 *
 *
 * Arg:        dfp [UNKN ] Undocumented argument [DYNFILE *]
 * Arg:        dpi [UNKN ] Undocumented argument [DPImplementation *]
 * Arg:         gm [UNKN ] Undocumented argument [GenericMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean write_thread_struct(DYNFILE * dfp,DPImplementation * dpi,GenericMatrix * gm);


/* Function:  write_thread_loop(dfp,dpi,gm)
 *
 * Descrip:    Major loop for each thread.
 *
 *             This function is really complicated by the way
 *             dynamite databases work - in particular the init
 *             function returns an object, so we have to keep
 *             track for each database if we need to init or not.
 *
 *             Secondly we have to get out the information for 
 *             each query and target regardless of whether they will
 *             score over the default or not. To do this we need 
 *             to rely on the datascore storage allocator to be
 *             able to accept frees (annoying or what!).
 *
 *             Here is the pseudo code for both a query or target db
 *
 *             forever {
 *               get input lock
 *               if( search ended ) {
 *                 unlock input lock
 *                 break;
 *               }
 *
 *               get datascore datatstructure
 *
 *               if( query not init )
 *                 init query
 *               hard link query
 *               read query info into datastructure
 *
 *               if( target not init )
 *                 init target
 *               else 
 *                 reload target
 *
 *               if ( end of the target database )
 *                 free query data structure  -- we have got it hard_linked!
 *                 reload query database into holder
 *                 if( end of the query database )
 *                    flag search finished
 *                    unlock input lock
 *                 else
 *                    close target database
 *                    set target to be initiated -- by next thread 
 *                    unlock input lock
 *               else
 *                 unlock input lock
 *
 *               do score
 *
 *               if( should store )
 *                 get output lock
 *                 add data score
 *               else
 *                 get input lock (otherwise clash with datascore storage allocator)
 *                 return data score back to storage
 *                 unlock input lock
 *             }
 *
 *
 * Arg:        dfp [UNKN ] Undocumented argument [DYNFILE *]
 * Arg:        dpi [UNKN ] Undocumented argument [DPImplementation *]
 * Arg:         gm [UNKN ] Undocumented argument [GenericMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean write_thread_loop(DYNFILE * dfp,DPImplementation * dpi,GenericMatrix * gm);


/* Function:  get_pre_chainstr_GenericMatrix(pre,gm)
 *
 * Descrip:    makes an the argument calling string which
 *             is compatible with the arg_str from
 *             get_argstr_GenericMatrix, but with a
 *             pre placement in it
 *
 *             eg "query,target,holder->comp_mat"
 *
 *
 * Arg:        pre [UNKN ] Undocumented argument [const char *]
 * Arg:         gm [READ ] structure holding generic matrix [const GenericMatrix *]
 *
 * Return [UNKN ]  allocated string (must free) with chained-args [char *]
 *
 */
char * get_pre_chainstr_GenericMatrix(const char * pre,const GenericMatrix * gm);


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
