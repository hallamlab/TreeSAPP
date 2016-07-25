

.SUFFIXES :  .dy




DYNHELPEROBJ = packaln.o \
	aln.o \
	dnamatrix.o \
	probability.o \
	alnrange.o \
	alnconvert.o \
	basematrix.o \
	matrixdebug.o \
	dpenvelope.o \
	dbsearchimpl.o \
	dprunimpl.o \
	complexsequence.o \
	complexevalset.o \
	complexconsensi.o \
	sequence.o \
	seqalign.o \
	codon.o \
	compmat.o \
	codonmatrix.o \
	codonmapper.o \
	sequencedb.o \
	hscore.o \
	histogram.o \
	proteindb.o \
	protein.o \
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
#CFLAGS =  -g -c  -DUNIX -DWISE_MEMORY_WATCH -I../base/
CFLAGS =  -g -c  -DUNIX  -I../base/
#CFLAGS =   -O -c -DUNIX -DNOERROR -I../base/
LFLAGS =  -g -L../base/ -lm 
AR_OPTIONS = ru

RANLIB_NEEDED = 0

#
# Base library build - call wisebase
#

all : libdyna.a 


libdyndb.a : $(DATABASEOBJ)
	ar $(AR_OPTIONS) libdyndb.a $?
#	if ( $(RANLIB_NEEDED)) then ranlib libline.a; fi

libdyna.a : $(DYNHELPEROBJ)
	ar $(AR_OPTIONS) libdyna.a $?

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

testpb : testpairbase.o pairbase.o pairbaseseq.o sequence.o seqalign.o codon.o
	$(CC) -o testpb testpairbase.o pairbase.o pairbaseseq.o sequence.o codon.o seqalign.o ../base/libwisebase.a -lm

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









