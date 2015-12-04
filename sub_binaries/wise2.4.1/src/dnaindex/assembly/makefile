

.SUFFIXES :  .dy

test_basegraph : test_basegraph.o basegraph.o
	$(CC) -o test_basegraph test_basegraph.o basegraph.o ../dnamapping.o ../kmer_index_interface.o ../kmer_glib_index.o $(LFLAGS) -ldyna -lwisebase

test_assemblygraph : test_assemblygraph.o assemblygraph.o assemblystats.o basegraph.o
	$(CC) -o test_assemblygraph test_assemblygraph.o  assemblystats.o basegraph.o assemblygraph.o ../dnamapping.o ../kmer_index_interface.o ../kmer_glib_index.o ../kmer_hash.o ../assembly.o ../assembly_stream_fasta.o ../largeseqreader.o ../singleseqspace.o ../assembly_stream_interface.o $(LFLAGS) -ldyna -lwisebase


test_graphtangle : test_graphtangle.o assemblygraph.o assemblystats.o basegraph.o graphtangle.o assemblypath.o 
	$(CC) -o test_graphtangle test_graphtangle.o  assemblystats.o basegraph.o assemblygraph.o graphtangle.o assemblypath.o ../dnamapping.o ../kmer_index_interface.o ../kmer_glib_index.o ../assembly.o ../assembly_stream_fasta.o ../largeseqreader.o ../singleseqspace.o ../assembly_stream_interface.o ../kmer_hash.o $(LFLAGS) -ldyna -lwisebase


pathwise : pathwise.o supertangle.o assemblygraph.o assemblystats.o basegraph.o graphtangle.o assemblypath.o grapherror.o graphmanager.o depthmanager.o 
	$(CC) -o pathwise pathwise.o  supertangle.o depthmanager.o assemblystats.o basegraph.o assemblygraph.o ../dnamapping.o ../kmer_index_interface.o graphtangle.o assemblypath.o grapherror.o graphmanager.o  ../kmer_glib_index.o ../kmer_hash.o ../assembly.o ../assembly_stream_fasta.o ../largeseqreader.o ../singleseqspace.o ../assembly_stream_interface.o ../assembly_sanger_project.o $(LFLAGS) -ldyna -lwisebase


badkmer : badkmer.o assemblygraph.o assemblystats.o basegraph.o
	$(CC) -o badkmer badkmer.o  assemblystats.o basegraph.o assemblygraph.o ../dnamapping.o ../kmer_index_interface.o ../kmer_glib_index.o ../kmer_hash.o ../assembly.o ../assembly_stream_fasta.o ../largeseqreader.o ../singleseqspace.o ../assembly_stream_interface.o $(LFLAGS) -ldyna -lwisebase

arraymatcher : arraymatcher.o assemblygraph.o assemblystats.o basegraph.o
	$(CC) -o arraymatcher arraymatcher.o  assemblystats.o basegraph.o assemblygraph.o ../dnamapping.o ../kmer_index_interface.o ../kmer_glib_index.o ../kmer_hash.o ../assembly.o ../assembly_stream_fasta.o ../largeseqreader.o ../singleseqspace.o ../assembly_stream_interface.o $(LFLAGS) -ldyna -lwisebase



CC = cc
CFLAGS = -Wall -pthread -c  -DHAS_PTHREAD_SETSCOPE -g -DUNIX -I../../base/ -I../../dynlibsrc/ -I../ -I../../models/ `glib-config --cflags`


LFLAGS =  -g -L../../base/ -L../../dynlibsrc/ -lpthread `glib-config --libs` -lm


.c.o :
	$(CC) $(CFLAGS) $(INCFLAGS) $?

INCFLAGS = -I../base/ -I../dynlibsrc/
DFLAGS = -l -n Wise2_ -pthreads -dbtrace 5 -nocwarn

.dy.c : 
	dyc  $(DFLAGS) $<	

.dy.o :
	dyc $(DFLAGS) $<	
	$(CC) $(CFLAGS) $*.c

init : 
	dyc $(DFLAGS) *.dy
