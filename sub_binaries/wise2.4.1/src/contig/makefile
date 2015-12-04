
.SUFFIXES :  .dy

.c.o :
	$(CC) $(CFLAGS) $(INCFLAGS) $?

INCFLAGS = -I../base/ -I../dynlibsrc/
DFLAGS = -l -D -n Wise2_ -a _api.h -b _api.t -latex -perl -pthreads -dbtrace 5 -nocwarn
#DFLAGS = -l -O -D -n Wise2_ -a _api.h -b _api.t -latex -perl -pthreads -dbtrace 5
#DFLAGS = -l -n Wise2_ 

LFLAGS =  -g  -L../base/ -L../dynlibsrc/


OBJ = contig.o overlap.o


testoverlap : testoverlap.o $(OBJ)
	$(CC) -o testoverlap testoverlap.o contig.o overlap.o -lwisebase -ldyna $(LFLAGS)

CC = cc
CFLAGS = -c -g  -I../base/ -I../dynlibsrc/

.dy.c : 
	dyc  $(DFLAGS) $<	

