

.SUFFIXES :  .dy


ASSEMBLY_OBJ  = kmer_assembly.o \
	kmer_index_interface.o \
	kmer_direct.o \
	kmer_hash.o \
	kmer_count.o \
	kmer_glib_index.o \
	singleseqspace.o \
	dnamapping.o \
	largeseqreader.o \
	kmer_assembly_untangler.o \
	kmer_assembly_contig.o \
	kmer_assembly_error.o \
	assembly.o \
	assembly_stream_interface.o \
	assembly_stream_fasta.o \
	assembly_sanger_project.o \
	assembly_stream_cons.o \
        compressed_protein_index.o


all: $(ASSEMBLY_OBJ)


kwise : kwise.o $(ASSEMBLY_OBJ)
	$(CC) -o kwise kwise.o $(ASSEMBLY_OBJ)  -ldyna -lwisebase $(LFLAGS) `glib-config --libs`

arraymatcher : arraymatcher.o $(ASSEMBLY_OBJ) 
	$(CC) -o arraymatcher arraymatcher.o $(ASSEMBLY_OBJ)  -ldyna -lwisebase $(LFLAGS) `glib-config --libs`


findbad_kmer : findbad_kmer.o $(ASSEMBLY_OBJ)
	$(CC) -o findbad_kmer findbad_kmer.o $(ASSEMBLY_OBJ)  -ldyna -lwisebase $(LFLAGS) `glib-config --libs`

make_shotgun : make_shotgun.o shotgun.o
	$(CC) -o make_shotgun make_shotgun.o shotgun.o -ldyna -lwisebase $(LFLAGS)

test_dnanumber : test_dnanumber.o dnanumber.o
	$(CC) -o test_dnanumber test_dnanumber.o dnanumber.o -ldyna -lwisebase $(LFLAGS)

test_linkindex : test_linkindex.o linkindex.o dnanumber.o
	$(CC) -o test_linkindex test_linkindex.o linkindex.o dnanumber.o -ldyna -lwisebase $(LFLAGS)

test_forwardindex : test_forwardindex.o forwardindex.o 
	$(CC) -o test_forwardindex test_forwardindex.o forwardindex.o -ldyna -lwisebase $(LFLAGS)

test_euler : test_euler.o eulerindex.o singleseqspace.o
	$(CC) -o test_euler test_euler.o eulerindex.o singleseqspace.o -ldyna -lwisebase $(LFLAGS)

test_kmer : test_kmer.o kmer_assembly.o kmer_index_interface.o kmer_direct.o singleseqspace.o dnamapping.o largeseqreader.o  kmer_assembly_untangler.o kmer_assembly_contig.o kmer_assembly_error.o assembly.o 
	$(CC) -o test_kmer  test_kmer.o kmer_assembly.o kmer_index_interface.o kmer_direct.o singleseqspace.o dnamapping.o  largeseqreader.o kmer_assembly_untangler.o   kmer_assembly_error.o kmer_assembly_contig.o assembly.o  -ldyna -lwisebase $(LFLAGS)

kmer_stat : kmer_stat.o kmer_assembly.o kmer_index_interface.o kmer_direct.o singleseqspace.o dnamapping.o largeseqreader.o 
	$(CC) -o kmer_stat kmer_stat.o kmer_assembly.o kmer_index_interface.o kmer_direct.o singleseqspace.o dnamapping.o  largeseqreader.o -ldyna -lwisebase $(LFLAGS)

test_compara : test_compara.o comparapath.o dnamapping.o singleseqspace.o chhash.o largeseqreader.o kmer_direct.o kmer_hash.o
	$(CC) -o test_compara test_compara.o comparapath.o dnamapping.o singleseqspace.o chhash.o largeseqreader.o kmer_direct.o kmer_hash.o -ldyna -lwisebase $(LFLAGS)

bigseqwise : bigseqwise.o comparapath.o dnamapping.o singleseqspace.o chhash.o largeseqreader.o kmer_direct.o kmer_hash.o
	$(CC) -o bigseqwise bigseqwise.o comparapath.o dnamapping.o singleseqspace.o chhash.o largeseqreader.o kmer_direct.o kmer_hash.o -ldyna -lwisebase $(LFLAGS)




show_euler : show_euler.o eulerindex.o singleseqspace.o
	$(CC) -o show_euler show_euler.o eulerindex.o singleseqspace.o -ldyna -lwisebase $(LFLAGS)

CC = cc
CFLAGS = -pthread -c  -DHAS_PTHREAD_SETSCOPE -g -DUNIX -I../base/ -I../dynlibsrc/ `glib-config --cflags`


LFLAGS =  -g -L../base/ -L../dynlibsrc/ -lpthread `glib-config --libs` -lm


.c.o :
	$(CC) $(CFLAGS) $(INCFLAGS) $?

INCFLAGS = -I../base/ -I../dynlibsrc/ -I../models
DFLAGS = -l -n Wise2_ -pthreads -dbtrace 5 -nocwarn

.dy.c : 
	dyc  $(DFLAGS) $<	

.dy.o :
	dyc $(DFLAGS) $<	
	$(CC) $(CFLAGS) $*.c

init : 
	dyc $(DFLAGS) *.dy
