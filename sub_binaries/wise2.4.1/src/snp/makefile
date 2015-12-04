
.SUFFIXES :  .dy


CC = cc
#CFLAGS = -c -O2  -pg -I../base/ -I../dynlibsrc/
#CFLAGS = -c -O2 -DPTHREAD -DHAS_PTHREAD_SETSCOPE  -DUNIX  -I../base/ -I../dynlibsrc/ -I.
CFLAGS = -c -g3 -O2  -DUNIX  -I../base/ -I../dynlibsrc/ -I. -pthread `glib-config --cflags`

LFLAGS =  -g -L../base/ -L../dynlibsrc/ -lpthread `glib-config --libs` -lpthread -lm
AR_OPTIONS = ru

RANLIB_NEEDED = 0

.c.o :
	$(CC) $(CFLAGS) $(INCFLAGS) $?

test_est : test_est.o genotype.o person.o locus_framework.o locus_model_estimators.o frequency_count.o 
	$(CC) -o test_est test_est.o genotype.o person.o locus_framework.o locus_model_estimators.o frequency_count.o ../dynlibsrc/libdyna.a ../base/libwisebase.a $(LFLAGS)

INCFLAGS = -I../base/ -I../dynlibsrc/
DFLAGS = -l -D -n Wise2_ -a _api.h -b _api.t -pthreads -dbtrace 5 -nocwarn

.dy.c : 
	dyc  $(DFLAGS) $<	

.dy.o :
	dyc $(DFLAGS) $<	
	$(CC) $(CFLAGS) $*.c

init : 
	dyc $(DFLAGS) *.dy

