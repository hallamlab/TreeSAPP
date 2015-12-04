


.SUFFIXES :  .dy



NBASEOBJ = wiseconfig.o \
	wisestring.o \
	wiseerror.o \
	wisememman.o \
	wisefile.o \
	wiserandom.o \
	wisetime.o \
	wiseoverlay.o \
	wisestreaminterface.o \
	commandline.o \
	linesubs.o



#CFLAGS = -c  -g -DUNIX -DWISE_MEMORY_WATCH
CFLAGS = -c  -O2 -DUNIX -pthread
LFLAGS =  -g  -lm 
DFLAGS = -n Wise2_ -a _api.h -b _api.t -l -D
AR_OPTIONS = ru

RANLIB_NEEDED = 0

#
# Base library build - call wisebase
#

libwisebase.a : $(NBASEOBJ) 
	ar $(AR_OPTIONS) libwisebase.a $?
	-if test -x /bin/ranlib; then /bin/ranlib libwisebase.a; else exit 0; fi
	-if test -x /usr/bin/ranlib; then /usr/bin/ranlib libwisebase.a; else exit 0; fi

init : 
	dyc $(DFLAGS) *.dy

.c.o :
	$(CC) $(CFLAGS) $?

.dy.c : 
	dyc $(DFLAGS)  $<	

.dy.o :
	dyc $(DFLAGS) $<	
	$(CC) $(CFLAGS) $*.c

#
# pro2DNAmatrix.c pro2DNAmatrix.h : pro2DNAmatrix.dy
#	dy pro2DNAmatrix.dy







