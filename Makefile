CC = g++ -Wall
CCFLAG= -m64
CFLAGS= -O3 -m64 -w
FLAGS= -lm -lpthread  

RPKM_SRC= sub_binaries/rpkm_calculator/
RPKM_SOURCES= $(RPKM_SRC)utilities.c++ $(RPKM_SRC)rpkm.c++ $(RPKM_SRC)helper.c++ $(RPKM_SRC)fastareader.c++ $(RPKM_SRC)matchoutputparser.c++
RPKM_OBJECTS= $(RPKM_SOURCES:.c++=.o)
RPKM_HEADERS= $(RPKM_SOURCES:.c++=.h)

FGS_SRC= sub_binaries/FragGeneScanPlus/
FGS_SOURCES=	$(FGS_SRC)util_lib.c $(FGS_SRC)hmm_lib.c $(FGS_SRC)run_hmm.c $(FGS_SRC)fasta.c
FGS_OBJECTS=	$(FGS_SRC)util_lib.o $(FGS_SRC)hmm_lib.o $(FGS_SRC)run_hmm.o $(FGS_SRC)fasta.o
FGS_HEADERS=	$(FGS_SRC)util_lib.h $(FGS_SRC)fasta.h $(FGS_SRC)run_hmm.h


%.o: %.c++   $(RPKM_SOURCES) $(RPKM_SRC)types.h
	$(CC) $(CCFLAG)  $< -c -o $@  

%.o: %.c $(FGS_SOURCES)
	gcc $(CFLAGS) -c $< -o $@ 

all: rpkm fgs 

rpkm: $(RPKM_OBJECTS) $(RPKM_HEADERS) $(RPKM_SRC)types.h
	$(CC) $(CCFLAG) $(RPKM_OBJECTS) -o rpkm
	
fgs:  $(FGS_OBJECTS) $(FGS_HEADERS)
	echo $(FGS_OBJECTS)
	gcc $(CFLAGS) $(FGS_OBJECTS) -o FGS+ $(FLAGS)

clean:
	rm -rf $(RPKM_OBJECTS) rpkm $(FGS_OBJECTS) FGS+

install:
	mv rpkm sub_binaries/
	mv FGS+ sub_binaries/

