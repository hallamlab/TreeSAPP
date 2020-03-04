CC = g++ -Wall
CCFLAG= -m64
CFLAGS= -O3 -m64 -w
FLAGS= -lm -lpthread  

TS_BIN_DIR=treesapp/sub_binaries

RPKM_SRC= $(TS_BIN_DIR)/rpkm_calculator
RPKM_SOURCES= $(RPKM_SRC)/utilities.cpp $(RPKM_SRC)/rpkm.cpp $(RPKM_SRC)/helper.cpp $(RPKM_SRC)/fastareader.cpp $(RPKM_SRC)/matchoutputparser.cpp
RPKM_OBJECTS= $(RPKM_SOURCES:.cpp=.o)
RPKM_HEADERS= $(RPKM_SOURCES:.cpp=.h)

HMMSEARCH_EXE := $(shell which hmmsearch)

ODSEQ_EXE := $(shell which OD-Seq)
ODSEQ_SOURCES= OD-Seq/AliReader.cpp OD-Seq/Bootstrap.cpp OD-Seq/DistCalc.cpp OD-Seq/DistMatReader.cpp \
	OD-Seq/DistWriter.cpp OD-Seq/FastaWriter.cpp OD-Seq/IQR.cpp OD-Seq/ODseq.cpp OD-Seq/PairwiseAl.cpp \
	OD-Seq/Protein.cpp OD-Seq/ResultWriter.cpp OD-Seq/runtimeargs.cpp OD-Seq/util.cpp

%.o: %.cpp $(RPKM_SOURCES) $(RPKM_SRC)/types.h
	$(CC) $(CCFLAG)  $< -c -o $@  

all: rpkm hmmer odseq

rpkm: $(RPKM_OBJECTS) $(RPKM_HEADERS) $(RPKM_SRC)/types.h
	$(CC) $(CCFLAG) $(RPKM_OBJECTS) -o $(TS_BIN_DIR)/rpkm

hmmer:
ifeq ($(HMMSEARCH_EXE),)
	curl -LJ0 --output hmmer-3.3.tar.gz http://eddylab.org/software/hmmer/hmmer-3.3.tar.gz
	tar -xzf hmmer-3.3.tar.gz; cd hmmer-3.3/; ./configure; make -f Makefile; cd -
	rm hmmer-3.3.tar.gz
else
	@echo HMMER found
endif

odseq:
ifeq ($(ODSEQ_EXE),)
	curl -LJ0 --output od-seq.tar.gz http://www.bioinf.ucd.ie/download/od-seq.tar.gz
	tar -xzf od-seq.tar.gz;
	g++ -fopenmp -o $(TS_BIN_DIR)/OD-seq $(ODSEQ_SOURCES)
	rm od-seq.tar.gz
else
	@echo OD-Seq found
endif

clean:
	rm -rf $(RPKM_OBJECTS) rpkm $(TS_BIN_DIR)/rpkm  _tree_parser*.so _fasta_reader*.so
	rm -r OD-Seq/
	rm -r hmmer-3.3/

install:
	cp rpkm $(TS_BIN_DIR)
	cp OD-Seq/OD-Seq $(TS_BIN_DIR)/
	cp hmmer-3.3/src/hmmsearch hmmer-3.3/src/hmmbuild hmmer-3.3/src/hmmalign hmmer-3.3/src/hmmfetch $(TS_BIN_DIR)/

