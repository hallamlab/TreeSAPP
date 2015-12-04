
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
#include <time.h>

#include <unistd.h>	  


#ifdef COMPILE_VERBOSITY
#define COMPILE_VERBOSITY_FLAG 1
#else
#define COMPILE_VERBOSITY_FLAG 0
#endif

#define VERBOSITY_CHECK(level,verbosity) ((COMPILE_VERBOSITY_FLAG == 1)&(verbosity > level))


#define PTHREAD
/**** OK some system wide defines now - used all over the place ****/

#define MAXLINE 512 /* generalised maximum input line */
#define MAXBINARYDUMP 1024 /*** ok could be tricky here... ****/

#ifndef BOOLEANDEFINED

typedef int boolean;
#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif


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

/* stream wrappers */
#include "wisestreaminterface.h"


#ifdef WISE_MEMORY_WATCH
#ifndef CKALLOC_GUARD
#undef  ckalloc
#define ckalloc(byte) allocate_watched_memory_file(__FILE__,__LINE__,byte)
#endif
#endif

#define error warn

#endif /* TOP_BASE.H loaded */
