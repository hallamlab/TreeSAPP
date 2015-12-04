
.SUFFIXES :  .dy

CC = cc
INCFLAGS = -I../base/
CFLAGS =   -O2 -g -c -pthread  -DUNIX  -I../base/
LFLAGS =  -g -L../base/ -lm 
AR_OPTIONS = ru

RANLIB_NEEDED = 0


SOCKET_OBJ = functionserver.o \
	functionclient.o \
	anonobj.o \
	transferinterface.o \
	directsocketwrite.o 


libwisesocket.a : $(SOCKET_OBJ)
	ar $(AR_OPTIONS) libwisesocket.a $?
	-if test -x /bin/ranlib; then /bin/ranlib libwisesocket.a; else exit 0; fi
	-if test -x /usr/bin/ranlib; then /usr/bin/ranlib libwisesocket.a; else exit 0; fi



test_server : test_cat_server.o simplebufferedserver.o 
	$(CC) -g -o test_server test_cat_server.o simplebufferedserver.o  ../base/wisestreaminterface.o ../base/libwisebase.a 

test_strcat_server : test_strcat_server.o functionserver.o transferinterface.o directsocketwrite.o anonobj.o
	$(CC) -g -o test_strcat_server  test_strcat_server.o functionserver.o transferinterface.o directsocketwrite.o anonobj.o ../base/wisestreaminterface.o ../base/libwisebase.a 

test_strcat_client : test_strcat_client.o functionclient.o anonobj.o transferinterface.o directsocketwrite.o 
	$(CC) -g -o test_strcat_client  test_strcat_client.o anonobj.o functionclient.o transferinterface.o directsocketwrite.o  ../base/wisestreaminterface.o ../base/libwisebase.a 


test_client : test_cat_client.o simplesocketclient.o 
	$(CC) -g -o test_client test_cat_client.o simplesocketclient.o  ../base/wisestreaminterface.o ../base/libwisebase.a 

.c.o :
	$(CC) $(CFLAGS) $(INCFLAGS) $?

init : 
	dyc $(DFLAGS) *.dy

DFLAGS = -l -n Wise2_ -a _api.h -b _api.t -latex -perl

.dy.c : 
	dyc  $(DFLAGS) $<	

.dy.o :
	dyc $(DFLAGS) $<	
	$(CC) $(CFLAGS) $*.c
