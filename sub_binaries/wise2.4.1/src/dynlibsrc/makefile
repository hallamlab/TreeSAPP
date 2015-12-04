

.SUFFIXES :  .dy


DYNGLIBOBJ = subseqhash.o \
	intallocator.o \
	proteinstreamedindex.o \
	shadowseq.o \
	shadowseqindex.o \
	hsphandler.o \
	hspscaninterface.o \
	hsp2hitscan.o \
	hsplookupscan.o \
	hsplookupthreaded.o \
	hspthreadeddb.o \
	hspscanruntime.o \
	hsptwohitscan.o \
	proteinindexcons.o \
	dnaindexcons.o \
	staticseq.o




DYNHELPEROBJ = packaln.o \
	aln.o \
	dnamatrix.o \
	probability.o \
	alnrange.o \
	alnconvert.o \
	basematrix.o \
	shattermatrix.o \
	matrixdebug.o \
	dpenvelope.o \
	dbsearchimpl.o \
	dprunimpl.o \
	complexsequence.o \
	complexevalset.o \
	complexconsensi.o \
	sequence.o \
	sequencestream.o \
	seqalign.o \
	hitlist.o \
	hsp.o \
	hspstream.o \
	codon.o \
	compmat.o \
	codonmatrix.o \
	codonmapper.o \
	sequencedb.o \
	hscore.o \
	seqlookup.o \
	arrayseqlookup.o \
	genericindexresult.o \
	linkedlist_lookpos.o \
	singlenumberspace.o \
	histogram.o \
	searchstatinterface.o \
	searchstatlookup.o \
	proteindb.o \
	protein.o \
	pairbase.o \
	pairbaseseq.o \
	genomicdb.o \
	randommodel.o \
	randomdb.o \
	genomic.o \
	cdna.o \
	cdnadb.o \
	dna.o \
	embl.o \
	genomicregion.o \
	gene.o \
	transcript.o \
	translation.o \
	btcanvas.o \
	asciibtcanvas.o \
	dynlibcross.o




CC = cc
INCFLAGS = -I../base/
#CFLAGS = -Wall -g -pedantic -c  -DUNIX -I../base/
CFLAGS =  -Wall -pthread -g3 -DCOMPILE_VERBOSITY -DPTHREAD -O2 -c  -DUNIX -I../base/ `glib-config --cflags`
#CFLAGS =  -Wall -g3 -pg -DPTHREAD -O2  -c  -DUNIX -I../base/ `glib-config --cflags`
#CFLAGS =   -O2 -c -pthread  -DUNIX  -I../base/ `glib-config --cflags`
#CFLAGS =   -O -c -DUNIX -DNOERROR -I../base/
LFLAGS =  -g -L../base/ -lm 
AR_OPTIONS = ru

RANLIB_NEEDED = 0

#
# Base library build - call wisebase
#

all : libdyna.a libdyna_glib.a


libdyndb.a : $(DATABASEOBJ)
	ar $(AR_OPTIONS) libdyndb.a $?



libdyna.a : $(DYNHELPEROBJ)
	ar $(AR_OPTIONS) libdyna.a $?
	-if test -x /bin/ranlib; then /bin/ranlib libdyna.a; else exit 0; fi
	-if test -x /usr/bin/ranlib; then /usr/bin/ranlib libdyna.a; else exit 0; fi


libdyna_glib.a : $(DYNGLIBOBJ)
	ar $(AR_OPTIONS) libdyna_glib.a $(DYNGLIBOBJ)
	-if test -x /bin/ranlib; then /bin/ranlib libdyna_glib.a; else exit 0; fi
	-if test -x /usr/bin/ranlib; then /usr/bin/ranlib libdyna_glib.a; else exit 0; fi


libline.a : $(INTER) 
	ar $(AR_OPTIONS) libline.a $?
#	if ( $(RANLIB_NEEDED)) then ranlib libline.a; fi


libprob.a : $(PROBOBJ)
	ar $(AR_OPTIONS) libprob.a $?

init : 
	dyc $(DFLAGS) *.dy

#
#
#

test_genomic : genomic.o test_genomic.o sequence.o codon.o
	$(CC) -o test_genomic test_genomic.o genomic.o sequence.o codon.o ../base/libwisebase.a -lm

test_tree : tree.o test_tree.o
	$(CC) -o test_tree test_tree.o tree.o libdyna.a ../base/libwisebase.a -lm

test_intallocator : intallocator.o test_intallocator.o 
	$(CC) -o test_intallocator intallocator.o test_intallocator.o  ../base/libwisebase.a -lm

test_psi : test_psi.o proteinstreamedindex.o 
	$(CC) -o test_psi test_psi.o  proteinstreamedindex.o ./libdyna.a ../base/libwisebase.a -lm

test_shadowindex : test_shadowindex.o shadowseqindex.o shadowseq.o libdyna.a
	$(CC) -g -o test_shadowindex test_shadowindex.o shadowseqindex.o shadowseq.o ./libdyna.a ../base/libwisebase.a -lm

test_seqhash : test_seqhash.o subseqhash.o hsp.o hitlist.o hspscaninterface.o hsplookupscan.o libdyna.a
	$(CC) -g -o test_seqhash `glib-config --libs` test_seqhash.o subseqhash.o hsp.o  hitlist.o hspscaninterface.o hsplookupscan.o libdyna.a ../base/libwisebase.a -lm  

test_hsp : test_hsp.o codon.o sequence.o hsp.o compmat.o probability.o 
	$(CC) -g -o test_hsp test_hsp.o codon.o sequence.o hsp.o compmat.o probability.o  ../base/libwisebase.a -lm `glib-config --libs`

test_topscore : test_topscore.o codon.o sequence.o hsp.o compmat.o probability.o 
	$(CC) -g -o test_topscore test_topscore.o codon.o sequence.o hsp.o compmat.o probability.o  ../base/libwisebase.a -lm `glib-config --libs`


testgendb : probability.o genomicdb.o genomic.o testgendb.o sequence.o codon.o sequencedb.o complexsequence.o complexevalset.o
	$(CC) -o testgendb probability.o genomicdb.o genomic.o testgendb.o sequence.o codon.o sequencedb.o complexsequence.o complexevalset.o ../base/libwisebase.a -lm

testseqalign : seqalign.o testseqalign.o sequence.o codon.o
	$(CC) -o testseqalign testseqalign.o seqalign.o sequence.o codon.o ../base/libwisebase.a -lm

testanc : libdyna.a
	$(CC) -o testanc ../models/testanc.o ../models/proteinsw.o ../models/abc.o ../models/pba.o ../models/sw_wrap.o libdyna.a ../base/libwisebase.a -lm

testdb : testdb.o sequencedb.o sequence.o codon.o probability.o
	$(CC) -o testdb testdb.o sequencedb.o sequence.o codon.o probability.o ../base/libwisebase.a -lm

testpdb : testpdb.o protein.o proteindb.o complexevalset.o complexsequence.o sequencedb.o sequence.o codon.o probability.o randomdb.o randommodel.o histogram.o 
	$(CC) -o testpdb testpdb.o protein.o proteindb.o complexsequence.o sequencedb.o sequence.o codon.o probability.o randomdb.o randommodel.o complexevalset.o histogram.o ../base/libwisebase.a -lm

testpb : testpairbase.o pairbase.o pairbaseseq.o sequence.o seqalign.o codon.o probability.o complexsequence.o complexevalset.o
	$(CC) -g -o testpb testpairbase.o pairbase.o pairbaseseq.o sequence.o codon.o seqalign.o probability.o complexsequence.o complexevalset.o ../base/libwisebase.a -lm

.c.o :
	$(CC) $(CFLAGS) $(INCFLAGS) $?

DFLAGS = -l -D -F -n Wise2_ -a _api.h -b _api.t -latex -perl

.dy.c : 
	dyc  $(DFLAGS) $<	

.dy.o :
	dyc $(DFLAGS) $<	
	$(CC) $(CFLAGS) $*.c

#
# pro2DNAmatrix.c pro2DNAmatrix.h : pro2DNAmatrix.dy
#	dy pro2DNAmatrix.dy





test_randomdb : test_randomdb.o libdyna.a
	$(CC) -o test_randomdb test_randomdb.o libdyna.a -lwisebase -L../base -lm









