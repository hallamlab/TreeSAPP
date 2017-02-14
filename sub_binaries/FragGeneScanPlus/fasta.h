/* fasta.h
 * Declarations for simple FASTA i/o library
 * SRE, Sun Sep  8 05:37:38 2002 [AA2721, transatlantic]
 * CVS $Id: fasta.h,v 1.1 2003/10/05 18:43:39 eddy Exp $
 */

#ifndef FASTA_H
#define FASTA_H

#include <stdio.h>

#define STRINGLEN 1000000

typedef struct fastafile_s {
    FILE *fp;
    char  buffer[STRINGLEN];
} FASTAFILE;

extern FASTAFILE *OpenFASTA(char *seqfile);
extern int        ReadFASTA(FASTAFILE *fp, char **ret_seq, char **ret_name, int *ret_L);
extern void       CloseFASTA(FASTAFILE *ffp);

#endif
