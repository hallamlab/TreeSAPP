CC = g++ -Wall
CCFLAG= -m64
CFLAGS= -O3 -m64 -w
FLAGS= -lm -lpthread  

TS_BIN_DIR=treesapp/sub_binaries
RPKM_SRC= $(TS_BIN_DIR)/rpkm_calculator/
RPKM_SOURCES= $(RPKM_SRC)utilities.cpp $(RPKM_SRC)rpkm.cpp $(RPKM_SRC)helper.cpp $(RPKM_SRC)fastareader.cpp $(RPKM_SRC)matchoutputparser.cpp
RPKM_OBJECTS= $(RPKM_SOURCES:.cpp=.o)
RPKM_HEADERS= $(RPKM_SOURCES:.cpp=.h)
HMMSEARCH_EXE := $(shell which hmmsearch)
ODSEQ_EXE := $(shell which OD-Seq)

%.o: %.cpp   $(RPKM_SOURCES) $(RPKM_SRC)types.h
	$(CC) $(CCFLAG)  $< -c -o $@  

all: rpkm hmmer odseq

rpkm: $(RPKM_OBJECTS) $(RPKM_HEADERS) $(RPKM_SRC)types.h
	$(CC) $(CCFLAG) $(RPKM_OBJECTS) -o $(TS_BIN_DIR)/rpkm

hmmer:
ifeq ($(HMMSEARCH_EXE),)
	curl -LJ0 --output hmmer-3.3.tar.gz http://eddylab.org/software/hmmer/hmmer-3.3.tar.gz
	tar -xzf hmmer-3.3.tar.gz; cd hmmer-3.3/
	./configure --prefix $(TS_BIN_DIR)/; make -f Makefile
	cd -
else
	@echo HMMER found
endif

odseq:
ifeq ($(ODSEQ_EXE),)
	curl -LJ0 --output od-seq.tar.gz http://www.bioinf.ucd.ie/download/od-seq.tar.gz
	tar -xzf od-seq.tar.gz; cd OD-Seq/
	g++ -fopenmp -o OD-seq \
	AliReader.cpp Bootstrap.cpp DistCalc.cpp DistMatReader.cpp DistWriter.cpp FastaWriter.cpp IQR.cpp ODseq.cpp \
	PairwiseAl.cpp Protein.cpp ResultWriter.cpp runtimeargs.cpp util.cpp
	cd -
else
	@echo OD-Seq found
endif

clean:
	rm -rf $(RPKM_OBJECTS) rpkm $(TS_BIN_DIR)/rpkm  _tree_parser*.so _fasta_reader*.so
	rm -r OD-Seq/

install:
	mv rpkm $(TS_BIN_DIR)
	mv OD-Seq/OD-Seq $(TS_BIN_DIR)/
	mv prodigal.linux $(TS_BIN_DIR)/

