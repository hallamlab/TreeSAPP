

.SUFFIXES :  .dy

CC = cc
CFLAGS =   -O2 -g3 -c -pthread  -DUNIX  -I../base/ -I../dynlibsrc `glib-config --cflags`
LFLAGS =  -g -L../base/ -L../dynlibsrc  -lm `glib-config --libs`
AR_OPTIONS = ru

RANLIB_NEEDED = 0


.c.o :
	$(CC) $(CFLAGS) $(INCFLAGS) $?

DFLAGS = -l -n Wise2_ -a _api.h -b _api.t -latex -perl

.dy.c : 
	dyc  $(DFLAGS) $<	

.dy.o :
	dyc $(DFLAGS) $<	
	$(CC) $(CFLAGS) $*.c

init:
	dyc $(DFLAGS) *.dy
