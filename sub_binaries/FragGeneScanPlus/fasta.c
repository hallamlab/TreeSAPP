/* Simple API for FASTA file reading
 * for Bio5495/BME537 Computational Molecular Biology
 * SRE, Sun Sep  8 05:35:11 2002 [AA2721, transatlantic]
 * CVS $Id: fasta.c,v 1.1 2003/10/05 18:43:39 eddy Exp $
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "fasta.h"

/* Function: OpenFASTA(), ReadFASTA(), CloseFASTA().
 * Date:     SRE, Sun Sep  8 06:39:26 2002 [AA2721, transatlantic]
 *
 * Purpose:  A very rudimentary FASTA file reading API. Designed
 *           for simplicity and clarity, not for robustness.
 *
 *           The API is:
 *
 *           FASTAFILE *ffp;
 *           char      *seq;
 *           char      *name;
 *           int        seqlen;
 *
 *           ffp = OpenFASTA(seqfile);
 *           while (ReadFASTA(ffp, &seq, &name, &seqlen)
 *           {
 *             do stuff with sequence;
 *             free(name);
 *             free(seq);
 *           }
 *           CloseFASTA(ffp);
 *
 * Args:
 *           seqfile   - name of a FASTA file to open.
 *           seq       - RETURN: one sequence
 *           name      - RETURN: name of the sequence
 *           seqlen    - RETURN: length of the sequence in residues
 *           ffp       - ptr to a FASTAFILE object.
 *
 * Commentary:
 *           The basic problem with reading FASTA files is that there is
 *           no end-of-record indicator. When you're reading sequence n,
 *           you don't know you're done until you've read the header line
 *           for sequence n+1, which you won't parse 'til later (when
 *           you're reading in the sequence n+1). One common trick for
 *           this is to implement a one-line "lookahead" buffer that you
 *           can peek at, before parsing later.
 *
 *           This buffer is kept in a small structure (a FASTAFILE), rather
 *           than in a static char[] in the function. This allows
 *           us to have multiple FASTA files open at once. The static approach
 *           would only allow us to have one file open at a time. ANSI C
 *           predates the widespread use of parallel programming. It was
 *           not overly concerned about the drawbacks of statics. Today,
 *           though, you should keep in mind that you may someday want to
 *           turn your program into a multithreaded, parallel program, and
 *           all functions in parallelized code must be "reentrant": able to
 *           be called a second time - with different arguments,
 *           and while the code in the first function call is still executing! -
 *           without overwriting or corrupting any static storage in the
 *           function. Statics have fewer uses now (for example, to
 *           test that some initialization code for a function is run once
 *           and only once.)
 *
 * Limitations:
 *           There is no error handling, for clarity's sake. Also,
 *           the parser is brittle. Improper FASTA files (for instance,
 *           blank lines between records) will cause unexpected
 *           behavior. Real file parsers are more complex.
 *           In real life, they have to deal with absolutely anything the user might
 *           pass as a "FASTA file"; and either parse it correctly,
 *           or detect that it's an invalid format and fail cleanly.
 *
 *           Lines are read in from the file using ANSI C's fgets(). fgets()
 *           requires a maximum buffer length (here, STRINGLEN, which is
 *           defined as 512 in bio5495.h). Some FASTA files have very long
 *           description lines, however; notably the NCBI NR database. Static
 *           limitations on things like line or sequence lengths should be
 *           avoided. An example of a replacement for fgets() that dynamically
 *           allocates its buffer size and allows any line length is
 *           SQUID's sre_fgets().
*
*           We use ANSI C's strtok() to parse the sequence name out of the line.
*           strtok() is deprecated in modern programs because it is not threadsafe.
*           (See comments above.) An example of a threadsafe version is
*           SQUID's sre_strtok().
*
* Returns:
*           OpenFASTA() returns a FASTAFILE pointer, or NULL on failure (for
    *           instance, if the file doesn't exist, or isn't readable).
*
*           ReadFASTA() returns 1 on success, or a 0 if there are no
*           more sequences to read in the file.
*
*           CloseFASTA() "always succeeds" and returns void.
*/
FASTAFILE * OpenFASTA(char *seqfile) {

    FASTAFILE *ffp;
    ffp = malloc(sizeof(FASTAFILE));

    if(strcmp(seqfile, "stdin") == 0) {
        ffp->fp = stdin;
    } else {
        ffp->fp = fopen(seqfile, "r");              /* Assume seqfile exists & readable!   */
        if (ffp->fp == NULL) {
            free(ffp);
            return NULL;
        }
        if ((fgets(ffp->buffer, STRINGLEN, ffp->fp)) == NULL) {
            free(ffp);
            return NULL;
        }
    }

    return ffp;
}

int
ReadFASTA(FASTAFILE *ffp, char **ret_seq, char **ret_name, int *ret_L) {
    char *s;
    static char *name=0;
    static char *seq =0;
    int   n;
    int   nalloc;

    /* Peek at the lookahead buffer; see if it appears to be a valid FASTA descline.
    */
    if (ffp->buffer[0] != '>') return 0;

    /* Parse out the name: the first non-whitespace token after the >
    */
    s  = strtok(ffp->buffer+1, " \t\n");
    //name = malloc(sizeof(char) * (strlen(s)+1));
    if( name==0)
        name = malloc(sizeof(char) * 1024);
    strcpy(name, s);

    /* Everything else 'til the next descline is the sequence.
     * Note the idiom for dynamic reallocation of seq as we
     * read more characters, so we don't have to assume a maximum
     * sequence length.
     */
    if( seq==0)
        seq = malloc(sizeof(char) * 1024);     /* allocate seq in blocks of 128 residues */
    nalloc = 128;
    n = 0;
    while (fgets(ffp->buffer, STRINGLEN, ffp->fp)) {
        if (ffp->buffer[0] == '>') break;	/* a-ha, we've reached the next descline */

        for (s = ffp->buffer; *s != '\0'; s++) {
            if (! isalpha(*s)) continue;  /* accept any alphabetic character */

            seq[n] = *s;                  /* store the character, bump length n */
            n++;
            if (nalloc == n) {        /* are we out of room in seq? if so, expand */
                /* (remember, need space for the final '\0')*/
                nalloc += 128;
                seq = realloc(seq, sizeof(char) * nalloc);
            }
        }
    }
    seq[n] = '\0';

    *ret_name = name;
    *ret_seq  = seq;
    *ret_L    = n;
    return 1;
}

void
CloseFASTA(FASTAFILE *ffp) {
    fclose(ffp->fp);
    free(ffp);
}
