CC = g++ -Wall
CCFLAG= -m64
CFLAGS= -O3 -m64 -w
FLAGS= -lm -lpthread  

RPKM_SRC= sub_binaries/rpkm_calculator/
RPKM_SOURCES= $(RPKM_SRC)utilities.c++ $(RPKM_SRC)rpkm.c++ $(RPKM_SRC)helper.c++ $(RPKM_SRC)fastareader.c++ $(RPKM_SRC)matchoutputparser.c++
RPKM_OBJECTS= $(RPKM_SOURCES:.c++=.o)
RPKM_HEADERS= $(RPKM_SOURCES:.c++=.h)

%.o: %.c++   $(RPKM_SOURCES) $(RPKM_SRC)types.h
	$(CC) $(CCFLAG)  $< -c -o $@  

all: rpkm 

rpkm: $(RPKM_OBJECTS) $(RPKM_HEADERS) $(RPKM_SRC)types.h
	$(CC) $(CCFLAG) $(RPKM_OBJECTS) -o rpkm
	
clean:
	rm -rf $(RPKM_OBJECTS) rpkm sub_binaries/rpkm  _tree_parser*.so _fasta_reader*.so

install:
	mv rpkm sub_binaries/
	python3 setup.py build_ext --inplace

# hmmbuild, hmmalign, Gblocks, genewise, raxmlHPC, tree_parser
