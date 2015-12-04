
.SUFFIXES :  .dy

CC = cc
INCFLAGS = -I../base/  -I../socket -I../dynlibsrc -I../dnaindex
CFLAGS =    -O2  -c -pthread  -DUNIX  -I../base/ -I../socket -I../dynlibsrc -I../dnaindex `glib-config --cflags`
LFLAGS =  -g -L../base/ -L../socket -L../dynlibsrc -L../dnaindex -lm `glib-config --libs`
AR_OPTIONS = ru

RANLIB_NEEDED = 0

all : scanwise_server client_multihspscan.o


scanwise_server : wise_proteinindex_server.o net_hspscan.o ../dynlibsrc/libdyna.a ../dynlibsrc/libdyna_glib.a ../dnaindex/compressed_protein_index.o ../dnaindex/kmer_index_interface.o ../dnaindex/singleseqspace.o ../dnaindex/kmer_direct.o 
	$(CC) -g -o scanwise_server wise_proteinindex_server.o net_hspscan.o ../dnaindex/compressed_protein_index.o ../dnaindex/kmer_index_interface.o ../dnaindex/singleseqspace.o ../dnaindex/kmer_direct.o  -ldyna_glib -ldyna -lwisesocket -lwisebase $(LFLAGS) -lpthread 


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
