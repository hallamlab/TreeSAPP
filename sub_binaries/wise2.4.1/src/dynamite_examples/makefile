

#
# Minimal dynamite makefile
# set suffix rules .dy->.c .c->.o
#
#
#

CC = cc
CFLAGS = -O -I../dynlibsrc/ -I../base
LFLAGS = -L../dynlibsrc/ -L../base/  -ldyna -lwisebase -lm

.SUFFIXES : .dy


.dy.c :
	../dyc/dyc $<

.dy.o :
	../dyc/dyc $<
	$(CC) -c $(CFLAGS) $*.c


.c.o :
	$(CC) -c $(CFLAGS) $<


psw : proteinsw.o seqaligndisplay.o psw.o
	$(CC) -o psw psw.o proteinsw.o seqaligndisplay.o $(LFLAGS)

est2gen : cdna2genomic.o estgendisplay.o est2gen.o
	$(CC) -o est2gen est2gen.o cdna2genomic.o estgendisplay.o $(LFLAGS)

dbsearch : proteinsw.o dbsearch.o seqaligndisplay.o 
	$(CC) -o dbsearch dbsearch.o proteinsw.o seqaligndisplay.o $(LFLAGS)

dbthread : dbsearch.o seqaligndisplay.o
	../dyc/dyc -pthreads proteinsw.dy
	$(CC) -c $(CFLAGS) -DPTHREAD proteinsw.c
	$(CC) -o dbsearch dbsearch.o proteinsw.o seqaligndisplay.o $(LFLAGS) -lpthread


