CC = g++ -Wall # -g
#CC = g++ -pg
#CC = g++ -Wall 
CCFLAGS=  -m64


PROG = rpkm
RPKM_SRC= sub_binaries/rpkm_calculator/
SOURCES= $(RPKM_SRC)utilities.c++ $(RPKM_SRC)rpkm.c++ $(RPKM_SRC)helper.c++ $(RPKM_SRC)fastareader.c++ $(RPKM_SRC)matchoutputparser.c++
OBJECTS= $(SOURCES:.c++=.o)
HEADERS= $(SOURCES:.c++=.h)

%.o: %.c++   $(SOURCES) $(RPKM_SRC)types.h
	$(CC) $(CCFLAGS)  $< -c -o $@  

all: $(PROG)

clean:
	rm -rf $(OBJECTS) $(PROG)

$(PROG): $(OBJECTS) $(HEADERS) $(RPKM_SRC)types.h
	$(CC) $(CCFLAGS) $(OBJECTS) -o $(PROG)
	
install:
	mv rpkm sub_binaries/
