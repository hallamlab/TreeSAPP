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

ODSEQ_EXE := $(shell which OD-seq)
ODSEQ_SOURCES= OD-Seq/AliReader.cpp OD-Seq/Bootstrap.cpp OD-Seq/DistCalc.cpp OD-Seq/DistMatReader.cpp \
	OD-Seq/DistWriter.cpp OD-Seq/FastaWriter.cpp OD-Seq/IQR.cpp OD-Seq/ODseq.cpp OD-Seq/PairwiseAl.cpp \
	OD-Seq/Protein.cpp OD-Seq/ResultWriter.cpp OD-Seq/runtimeargs.cpp OD-Seq/util.cpp

%.o: %.cpp $(RPKM_SOURCES) $(RPKM_SRC)/types.h
	$(CC) $(CCFLAG)  $< -c -o $@  

all: rpkm hmmer odseq

rpkm: $(RPKM_OBJECTS) $(RPKM_HEADERS) $(RPKM_SRC)/types.h
	$(CC) $(CCFLAG) $(RPKM_OBJECTS) -o $(RPKM_SRC)/rpkm

hmmer:
ifeq ($(HMMSEARCH_EXE),)
	curl -LJ0 --output hmmer-3.3.tar.gz http://eddylab.org/software/hmmer/hmmer-3.3.tar.gz
	tar -xzf hmmer-3.3.tar.gz; cd hmmer-3.3/; ./configure; make -f Makefile; cd -
	rm hmmer-3.3.tar.gz
	HMMER_DIR = "hmmer-3.3/"
else
	@echo "HMMER found ($(HMMSEARCH_EXE))"
    HMMER_DIR := $(dir $(HMMSEARCH_EXE))
endif

odseq:
ifeq ($(ODSEQ_EXE),)
	curl -LJ0 --output od-seq.tar.gz http://www.bioinf.ucd.ie/download/od-seq.tar.gz
	tar -xzf od-seq.tar.gz; g++ -fopenmp -o OD-Seq/OD-seq $(ODSEQ_SOURCES)
	rm od-seq.tar.gz
	ODSEQ_DIR = "OD-Seq/"
else
	@echo "OD-seq found ($(ODSEQ_EXE))"
    ODSEQ_DIR := $(dir $(ODSEQ_EXE))
endif

clean:
	rm -rf $(RPKM_OBJECTS) $(RPKM_SRC)/rpkm  _tree_parser*.so _fasta_reader*.so
	rm -r OD-Seq/
	rm -r hmmer-3.3/

install:
	cp $(RPKM_SRC)/rpkm $(TS_BIN_DIR)/
	cp $(ODSEQ_DIR)/OD-seq $(TS_BIN_DIR)/
	cp $(HMMER_DIR)/hmmsearch $(HMMER_DIR)/hmmbuild $(HMMER_DIR)/hmmalign $(TS_BIN_DIR)/

