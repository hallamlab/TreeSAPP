

all : libmott.a

clean :
	rm -f *.[oa]

libmott.a : mott_api.o gaplib.o wise2_mott_bridge.o
	ar $(AR_OPTIONS) libmott.a $?
	-if test -x /bin/ranlib; then /bin/ranlib libmott.a; else exit 0; fi
	-if test -x /usr/bin/ranlib; then /usr/bin/ranlib libmott.a; else exit 0; fi


wise2_mott_bridge.o : wise2_mott_bridge.c
	$(CC) $(CFLAGS) -I../../dynlibsrc -I../../base wise2_mott_bridge.c

CC = cc

CFLAGS =  -Wall -g3 -DPTHREAD -O2 -c 
LFLAGS =  -g -L../base/ -lm 
AR_OPTIONS = ru


