
#ifndef TOPBASEHEADER
#define TOPBASEHEADER


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdarg.h>
#include <limits.h>
#include <float.h>
#include <assert.h>

#include <unistd.h>	  


/**** OK some system wide defines now - used all over the place ****/

#define MAXLINE 512 /* generalised maximum input line */
#define MAXBINARYDUMP 1024 /*** ok could be tricky here... ****/

#ifndef BOOLEANDEFINED

typedef int boolean;
#define TRUE 1
#define FALSE 0

#define BOOLEANDEFINED
#endif

#ifdef PTHREAD
#include <pthread.h>
#endif

/**** include the rest of the base files ****/

#include "wisestring.h"
#include "wisefile.h"
#include "wiseconfig.h"
#include "wisetime.h"
#include "wiserandom.h"
#include "wisememman.h"  /* memory manager - ckalloc/ckfree etc */
#include "wiseerror.h"
#include "wiseoverlay.h"
#include "commandline.h"

#ifdef WISE_MEMORY_WATCH
#ifndef CKALLOC_GUARD
#undef  ckalloc
#define ckalloc(byte) allocate_watched_memory_file(__FILE__,__LINE__,byte)
#endif
#endif


#endif /* TOP_BASE.H loaded */
