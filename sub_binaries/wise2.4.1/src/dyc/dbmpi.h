#ifndef DYNAMITEdbmpiHEADERFILE
#define DYNAMITEdbmpiHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif

#include "dyna2.h"
#include "dpimpl.h"



    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  write_MPI_master_loop(dfp,dpi,gm)
 *
 * Descrip:    The master mpi process loop.
 *
 *             The essential idea here is that there is a queue of processes waiting 
 *             for jobs; they indicate their readiness by sending a message to the master.
 *             The master then reads in the score (if necessary), and sends off a job to 
 *             the slave (if there are more jobs to be done).
 *
 *             Here is the pseudo code:
 *
 *             check number of processes, and allocate space for (n-1) datascore objects
 *
 *             initialize the databases
 *
 *             forever {
 *
 *               Receive a message
 *               if the message contains data
 *                 process it
 *
 *               if there are no more jobs to be done
 *                 send a message to the slave, telling it to die
 *
 *               if all the slaves have been notified
 *                 Terminate MPI and return.
 *
 *               if we have not done everything
 *                 Pack data, and dispatch; do bookkeeping, reload next target
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
boolean write_MPI_master_loop(DYNFILE * dfp,DPImplementation * dpi,GenericMatrix * gm);


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
boolean write_MPI_slave_loop(DYNFILE * dfp,DPImplementation * dpi,GenericMatrix * gm);


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/

#ifdef _cplusplus
}
#endif

#endif
