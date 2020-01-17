CC = g++ -Wall
CCFLAG= -m64
CFLAGS= -O3 -m64 -w
FLAGS= -lm -lpthread  

RPKM_SRC= treesapp/sub_binaries/rpkm_calculator/
RPKM_SOURCES= $(RPKM_SRC)utilities.cpp $(RPKM_SRC)rpkm.cpp $(RPKM_SRC)helper.cpp $(RPKM_SRC)fastareader.cpp $(RPKM_SRC)matchoutputparser.cpp
RPKM_OBJECTS= $(RPKM_SOURCES:.cpp=.o)
RPKM_HEADERS= $(RPKM_SOURCES:.cpp=.h)

%.o: %.cpp   $(RPKM_SOURCES) $(RPKM_SRC)types.h
	$(CC) $(CCFLAG)  $< -c -o $@  

all: rpkm 

rpkm: $(RPKM_OBJECTS) $(RPKM_HEADERS) $(RPKM_SRC)types.h
	$(CC) $(CCFLAG) $(RPKM_OBJECTS) -o treesapp/sub_binaries/rpkm
	
clean:
	rm -rf $(RPKM_OBJECTS) rpkm treesapp/sub_binaries/rpkm  _tree_parser*.so _fasta_reader*.so

install:
	mv rpkm treesapp/sub_binaries/
	python3 setup.py install

# hmmbuild, hmmalign, raxmlHPC, tree_parser
