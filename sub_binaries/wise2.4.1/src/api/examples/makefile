

CC = cc
CFLAGS = -c -I../include
LIBDIRS = -L../libs

all : revcomp sw

revcomp : revcomp.o
	$(CC) -o revcomp revcomp.o $(LIBDIRS) -ldyna -lwisebase  -lm

sw : sw.o
	$(CC) -o sw sw.o $(LIBDIRS) -lmodel -ldyna -lwisebase -lm

clean : 
	rm -f *.o

.c.o :
	$(CC) $(CFLAGS) $?
