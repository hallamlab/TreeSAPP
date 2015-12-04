

/* Helper functions in the module
 *
 * Wise2_change_max_BaseMatrix_kbytes
 * Wise2_get_max_BaseMatrix_kbytes
 *



/* These functions are not associated with an object */
/* Function:  Wise2_change_max_BaseMatrix_kbytes(new_kilo_number)
 *
 * Descrip:    This is to change, at run-time the maximum level of bytes basematrix *thinks*
 *             it can use. This number is *not* used for any actual calls to basematrix
 *             allocation: it is only used with /get_max_BaseMatrix_kbytes
 *
 *
 * Arg:        new_kilo_number max kilobytes allowed [int]
 *
 * Returns Undocumented return value [void]
 *
 */
void Wise2_change_max_BaseMatrix_kbytes( int new_kilo_number);

/* Function:  Wise2_get_max_BaseMatrix_kbytes(void)
 *
 * Descrip:    returns the max. number of kilobytes suggested as a limited
 *             to BaseMatrix. 
 *
 *
 *
 * Returns Undocumented return value [int]
 *
 */
int Wise2_get_max_BaseMatrix_kbytes();

