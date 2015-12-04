/************************************************************
 * HMMER - Biological sequence analysis with profile HMMs
 * Copyright (C) 1992-1998 Washington University School of Medicine
 *
 *   This source code is distributed under the terms of the 
 *   GNU General Public License. See the files COPYING and 
 *   GNULICENSE for details.
 *    
 ************************************************************/

/* version.h
 * Version and copyright stamp information for the package.
 * RCS $Header: /nfs/ensembl/cvsroot/wise2/src/HMMer2/version.h,v 1.1.1.1 2001/06/18 13:59:52 birney Exp $
 */


#define RELEASE     "2.0"
#define RELEASEDATE "June 1998"

#ifndef LINTING			/* compiled-in stamp for all binaries  */
static char *copyright = "\
   HMMER -- Biological sequence analysis with profile HMMs.\n\
   Copyright (C) 1992-1998 Washington University School of Medicine.\n\
\n\
   HMMER is distributed under the terms of the GNU General Public\n\
   License. See the files COPYING and GNULICENSE for details.";
#endif /* LINTING */
