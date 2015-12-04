


#
# This is the top level makefile for Wise2. It drives the underlying
# makefile - for normal usage you do not need to cd into the underlying
# directories.
#

#
# all is the target for the binaries in the system
#
# clean does remove all .o and .a files and moves the current
# binaries over to oldbin
#
#
#

# You may want to change these lines if you have
# a non ansi C compiler

# Change the cc line to gcc if you get compiler errors
# Most likely you have a non ANSI C compiler
# CC =  /opt/intel_cc_80/bin/icc
CC = cc

# These are the CFLAGS to use. These days Wise2 is inherently
# pthreaded, and everything is compiled with threads. However
# the main algorithms don't really take advantage of pthreads,
# it is more about making sure the libaries are thread safe


# Intel icc flags
# CFLAGS = -c  -O3 -axWK -ipo `glib-config --cflags`

# alpha flags (-pthread)
# CFLAGS = -c -O3 -pthread `glib-config --cflags`

# normal linux/bsd/mac flags
CFLAGS = -c -O3  `glib-config --cflags`


EXTRALIBS = -lm


#
# To compile with a different version of HMMER, change the following
# lines. This is not recommended, because sean and I diverged a while
# back on indexing stratgies. At some point I will need to port the later
# HMMer into Wise2
#



HMMER_INCLUDE = ../HMMer2/
HMMER_LIBS    = ../HMMer2/
HMMER_DEFINE  = HMMER_INTERNAL

#HMMER_DEFINE = HMMER_2_1_2=1
#HMMER_INCLUDE = /nfs/disk21/birney/prog/hmmer-2.1.2/src/
#HMMER_LIBS    = /nfs/disk21/birney/prog/hmmer-2.1.2/src/


#
# You shouldn't want to change anything below this.
#

start     :
	@echo 'to make the wise2 development directory read INSTALL\n\nmake all will try to make the package\n'

boot : 
	cd bootstrap ; $(MAKE) dyc
	cd dyc ;  ../bootstrap/dyc -l *.dy 
	cd base ; ../bootstrap/dyc -n Wise2_ -l *.dy 
	@echo '\n\nDon't worry about the above errors: now make newdyc'

newdyc : 
	(cd base   ; $(MAKE) CC="$(CC)" CFLAGS="$(CFLAGS)" libwisebase.a )
	(cd dyc   ; $(MAKE) CC="$(CC)" CFLAGS="$(CFLAGS)" dyc )
	@echo '\n\nNow move dyc/dyc to somewhere on your path'
	@echo 'setenv WISECONFIGDIR to wisecfg and make init followed by make all'

all      : realall
	$(MAKE) bin

realall  :
	(cd base   ; $(MAKE) CC="$(CC)" CFLAGS="$(CFLAGS)" libwisebase.a )
	(cd HMMer2 ; $(MAKE) CC="$(CC)" CFLAGS="$(CFLAGS)" libhmmer.a )
	(cd dynlibsrc ; $(MAKE) CC="$(CC)" CFLAGS="$(CFLAGS)" libdyna.a )
	(cd dynlibsrc ; $(MAKE) CC="$(CC)" CFLAGS="$(CFLAGS)" libdyna_glib.a )
	(cd external ;  $(MAKE) CC="$(CC)" CFLAGS="$(CFLAGS)" all )
	(cd socket ; $(MAKE) CC="$(CC)" CFLAGS="$(CFLAGS)" libwisesocket.a )
	(cd dnaindex ; $(MAKE) CC="$(CC)" CFLAGS="$(CFLAGS)" all )
	(cd network ; $(MAKE)  CC="$(CC)" CFLAGS="$(CFLAGS)" all )
	(cd models  ; $(MAKE) CC="$(CC)" CFLAGS="$(CFLAGS)" EXTRALIBS="$(EXTRALIBS)" HMMER_DEFINE="$(HMMER_DEFINE)" HMMER_INCLUDE="$(HMMER_INCLUDE)" HMMER_LIBS="$(HMMER_LIBS)" all )

strict :
	(cd base   ; $(MAKE) CC="gcc" CFLAGS="-g -c -Wall -DNOERROR" libwisebase.a )
	(cd HMMer2 ; $(MAKE) CC="gcc" CFLAGS="-g -c -Wall" libhmmer.a )
	(cd dynlibsrc ; $(MAKE) CC="gcc" CFLAGS="-g -c -Wall" libdyna.a )
	(cd models  ; $(MAKE) CC="gcc" CFLAGS="-g -c -Wall" all )

debug :
	(cd base   ; $(MAKE) CC="cc" CFLAGS="-g -c -DNOERROR" libwisebase.a )
	(cd HMMer2 ; $(MAKE) CC="cc" CFLAGS="-g -c" libhmmer.a )
	(cd dynlibsrc ; $(MAKE) CC="cc" CFLAGS="-g -c" libdyna.a )
	(cd models  ; $(MAKE) CC="cc" CFLAGS="-g -c" all )

bin : 
	mkdir bin
	cp models/pswdb models/psw models/genewisedb models/estwisedb models/estwise models/genewise models/dba models/dnal models/promoterwise network/scanwise_server models/scanwise ./bin
	csh welcome.csh

libs : 
	(cd base   ; $(MAKE) CC="$(CC)" CFLAGS="$(CFLAGS)" libwisebase.a )
	(cd HMMer2 ; $(MAKE) CC="$(CC)" CFLAGS="$(CFLAGS)" libhmmer.a )
	(cd dynlibsrc ; $(MAKE) CC="$(CC)" CFLAGS="$(CFLAGS)" libdyna.a )
	(cd models  ; $(MAKE) CC="$(CC)" CFLAGS="$(CFLAGS)" libmodel.a )

init :
	cd base ; $(MAKE) init
	cd dynlibsrc; $(MAKE) init
	cd models ; $(MAKE) init
	cd dyc ; $(MAKE) init
	cd socket ; $(MAKE) init
	cd dnaindex ; $(MAKE) init
	cd network; $(MAKE) init

clean : 
	cd external ; $(MAKE) clean
	-if test -d dynlibsrc; then cd dynlibsrc ; rm -f *.[oa]; fi
	-if test -d models; then cd models ; rm -f *.[oa]; fi
	-if test -d base; then cd base ; rm -f *.[oa]; fi
	-if test -d socket; then cd socket ; rm -f *.[oa]; fi
	-if test -d dnaindex; then cd dnaindex ; rm -f *.[oa]; fi
	-if test -d network; then cd network ; rm -f *.[oa]; fi
	-if test -d dyc; then cd dyc ; rm -f *.[oa]; fi
	-if test -d HMMer2; then cd HMMer2 ; rm -f *.[oa]; fi
	-if test -d perl; then cd perl/Wise2/libs ; rm -f *.[oa]; fi
	-if test -x perl/Wise2/Makefile; then cd perl/Wise2/ ; $(MAKE) clean; fi
	-if test -d oldbin; then rm -rf oldbin; fi
	-if test -d bin; then echo 'moving binaries to oldbin'; mv -f bin oldbin; fi


api  : libs apimake

apimake :
	cd api/include ; $(MAKE) apih
	cd api/libs ; $(MAKE) apilibs
	@echo 'Now cd into api/examples and try some out!'

perl : libs perlmake

perlmake :
	cd perl ; $(MAKE) int libs
	cd perl/Wise2 ; perl Makefile.PL
	cd perl/Wise2 ; $(MAKE) 
	@echo 'Now cd into perl/Wise2 and make install, once becoming root/perl master'

test : test_dummy

test_dummy :
	cd test ; perl testman.pl -cf ../wisecfg < short.test



removedyc :
	(cd base ; rm -f *.dy)
	(cd dynlibsrc ; rm -f *.dy)
	(cd models ; rm -f *.dy)
	rm -rf dyc
	rm -rf bootstrap

docs : realdocs

realdocs :
	(cd docs ; $(MAKE) api)
	(cd docs ; $(MAKE) progdoc )


dycdist :
	(cd base ; rm -f *.dy)
	(cd dynlibsrc ; rm -f *.dy)
	(rm -rf models)
	(rm -rf api)
	(rm -rf HMMer2)
	(rm -rf perl)
	(rm -rf test)
	(rm -rf INSTALL)
	(rm welcome.csh)
	(rm -rf example)
	(cd dyc ; rm -f *.dy)
	(rm -rf bootstrap)
	(mv -f DYNAMITE_README README)
	(mv -f makefile.dynamite makefile)
	(mv dynamite_examples examples)

