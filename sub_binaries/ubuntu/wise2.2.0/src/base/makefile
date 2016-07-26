# makefile for wise tools
#
# there are various little tricky things here...
# OSF/1 machine should be compiled with -DOSFHACK
# Dynamic.c has to be compiled with either
# -DSEARCHLINK for searchwise or with nothing for pairwise.
#
# I usualy do this by going
# make pairwise
# cc -c -DUNIX -DSEARCHLINK dynamic.c - NB for OSF uses, also -DOSFHACK
# make searchwise
#
# ewan  birney@molbiol.ox.ac.uk
#


.SUFFIXES :  .dy



NBASEOBJ = wiseconfig.o \
	wisestring.o \
	wiseerror.o \
	wisememman.o \
	wisefile.o \
	wiserandom.o \
	wisetime.o \
	wiseoverlay.o \
	commandline.o \
	linesubs.o



#CFLAGS = -c  -g -DUNIX -DWISE_MEMORY_WATCH
CFLAGS = -c  -g -DUNIX 
LFLAGS =  -g  -lm 
DFLAGS = -n Wise2_ -a _api.h -b _api.t -l -D
AR_OPTIONS = ru

RANLIB_NEEDED = 0

#
# Base library build - call wisebase
#

libwisebase.a : $(NBASEOBJ) 
	ar $(AR_OPTIONS) libwisebase.a $?

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







