

/* Helper functions in the module
 *
 * Wise2_openfile
 *



/* These functions are not associated with an object */
/* Function:  Wise2_openfile(filename,passedprot)
 *
 * Descrip:    Every file open goes through this.
 *
 *             It opens for reading in the following order 
 *                .
 *                WISEPERSONALDIR
 *                WISECONFIGDIR
 *
 *             For writing it opens in .
 *
 *             Filenames with ~'s are expanded to HOME/filename
 *
 *
 * Arg:        filename     filename to open for read/writing [const char *]
 * Arg:        passedprot   string representing standard fopen attributes [const char *]
 *
 * Returns open'd filehandle, NULL on error [FILE *]
 *
 */
FILE * Wise2_openfile( const char * filename,const char * passedprot);

