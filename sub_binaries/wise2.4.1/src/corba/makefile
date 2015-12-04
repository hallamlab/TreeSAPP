

CC = cc
INCFLAGS = -I../base/ -I../dynlibsrc
#gCFLAGS = -Wall -g -pedantic -c  -DUNIX -I../base/
#CFLAGS =  -O2 -c  -DUNIX  -I../base/ `glib-config --cflags` -I/usr/local/include/orbit-1.0/
CFLAGS =  -pthread -c  -DUNIX  -I../base/ -I../dynlibsrc `orbit-config --cflags client`
LFLAGS =   -L../base/ -lm -lpthread
AR_OPTIONS = ru



hspscan_corba.c : hspscan_corba.idl
	orbit-idl hspscan_corba.idl

test_server : test_server.o hspscan_server_impl.o hspscan_corba-common.o hspscan_corba-skels.o hspscan_corba-stubs.o corba_singleton.o 
	cc -g -o test_server  test_server.o hspscan_server_impl.o  hspscan_corba-common.o hspscan_corba-skels.o hspscan_corba-stubs.o  corba_singleton.o ../dynlibsrc/libdyna.a ../base/libwisebase.a ../dynlibsrc/hsp.o ../dynlibsrc/subseqhash.o ../dynlibsrc/linkedlist_lookpos.o   ../dynlibsrc/libdyna_glib.a `glib-config --libs` -lm -lORBit -lIIOP -lORBitutil -lpthread

scanwise_protein_index : scanwise_protein_index.o hspscan_server_impl.o hspscan_corba-common.o hspscan_corba-skels.o hspscan_corba-stubs.o corba_singleton.o 
	cc -o scanwise_protein_index  scanwise_protein_index.o hspscan_server_impl.o  hspscan_corba-common.o hspscan_corba-skels.o hspscan_corba-stubs.o  corba_singleton.o  ../dynlibsrc/libdyna_glib.a ../dynlibsrc/libdyna.a ../base/libwisebase.a `glib-config --libs` -lm -lORBit -lIIOP -lORBitutil -lpthread


test_client : test_client.o hspscan_corba-common.o hspscan_corba-stubs.o 
	cc -g -o test_client  test_client.o  hspscan_corba-common.o hspscan_corba-stubs.o ../dynlibsrc/libdyna.a ../base/libwisebase.a `glib-config --libs` -lm -lORBit -lIIOP -lORBitutil -lpthread


test_wrapper : test_wrapper.o hspscan_corba-common.o hspscan_corba-stubs.o hspscan_corba_wrapper.o corba_singleton.o 
	cc -g -o test_wrapper test_wrapper.o  hspscan_corba-common.o hspscan_corba-stubs.o   corba_singleton.o hspscan_corba_wrapper.o ../dynlibsrc/hsp.o ../dynlibsrc/libdyna.a ../base/libwisebase.a  `glib-config --libs` -lm -lORBit -lIIOP -lORBitutil -lpthread

.c.o :
	$(CC) $(CFLAGS) $(INCFLAGS) $?


