# Makefile for HMMER: Source directory
# RCS $Id: Makefile,v 1.1.1.1 2001/06/18 13:59:47 birney Exp $
##########
# HMMER - Biological sequence analysis with HMMs
# Copyright (C) 1992-1998 Sean R. Eddy
#
#   This source code is distributed under the terms of the 
#   GNU General Public License. See the files COPYING and 
#   LICENSE for details.
#    
###########

## your compiler and compiler flags
#
CC     = cc
CFLAGS = -g 

## machine specific definitions
#
MDEFS = 

## archiving command, and ranlib command if you need one.
#  In general, you shouldn't need to change these, and they're
#  only used for building the Shiva testsuite anyway.
#
AR     = ar rcv 
RANLIB = ranlib



#######
## You should not need to modify below this line
#######
SHELL  = /bin/sh
LIBS   = -lm

PROGS = hmmalign hmmbuild hmmcalibrate hmmconvert hmmemit hmmpfam hmmsearch 

HDRS =  config.h funcs.h globals.h structs.h version.h

SQUIDHDRS = gnuregex.h squid.h sqfuncs.h 

SRC =   hmmalign.c    hmmbuild.c        hmmcalibrate.c    hmmconvert.c\
        hmmpfam.c     hmmsearch.c\
	alphabet.c    core_algorithms.c debug.c           emit.c\
	emulation.c   histogram.c       hmmio.c           mathsupport.c\
        masks.c       misc.c            modelmakers.c     plan7.c\
	plan9.c	      prior.c           tophits.c         trace.c

SQUIDSRC = aligneval.c  alignio.c   cluster.c   dayhoff.c\
	   file.c       getopt.c    gnuregex.c  interleaved.c\
	   iupac.c      msf.c       revcomp.c   selex.c\
	   sqerror.c    sqio.c      sre_ctype.c sre_math.c\
	   sre_string.c stack.c     translate.c types.c\
	   weight.c

OBJ =  	alphabet.o    core_algorithms.o debug.o	          emit.o\
	emulation.o   histogram.o       hmmio.o           mathsupport.o\
        masks.o       misc.o            modelmakers.o     plan7.o\
	plan9.o       prior.o           tophits.o         trace.o\
        aligneval.o   alignio.o         cluster.o         dayhoff.o\
	file.o        getopt.o          gnuregex.o        interleaved.o\
	iupac.o       msf.o             revcomp.o         selex.o\
	sqerror.o     sqio.o            sre_ctype.o       sre_math.o\
	sre_string.o  stack.o           translate.o       types.o\
	weight.o

.c.o:
	$(CC) $(CFLAGS) $(MDEFS) -c $<		

#######
## Targets defining how to make HMMER executables.
#######

all: 	$(PROGS)

hmmalign: $(OBJ) hmmalign.o
	$(CC) $(CFLAGS) $(MDEFS) -o $@ hmmalign.o $(OBJ) $(LIBS) 

hmmbuild: $(OBJ) hmmbuild.o
	$(CC) $(CFLAGS) $(MDEFS) -o $@ hmmbuild.o $(OBJ) $(LIBS) 

hmmcalibrate: $(OBJ) hmmcalibrate.o
	$(CC) $(CFLAGS) $(MDEFS) -o $@ hmmcalibrate.o $(OBJ) $(LIBS) 

hmmconvert: $(OBJ) hmmconvert.o
	$(CC) $(CFLAGS) $(MDEFS) -o $@ hmmconvert.o $(OBJ) $(LIBS) 

hmmemit: $(OBJ) hmmemit.o
	$(CC) $(CFLAGS) $(MDEFS) -o $@ hmmemit.o $(OBJ) $(LIBS) 

hmmpfam: $(OBJ) hmmpfam.o
	$(CC) $(CFLAGS) $(MDEFS) -o $@ hmmpfam.o $(OBJ) $(LIBS) 

hmmsearch: $(OBJ) hmmsearch.o
	$(CC) $(CFLAGS) $(MDEFS) -o $@ hmmsearch.o $(OBJ) $(LIBS) 

#######
## Target used in making Shiva test suite
#######

libhmmer.a: $(OBJ)
	$(AR) libhmmer.a $(OBJ)
	-if test -x /bin/ranlib; then /bin/ranlib libhmmer.a; else exit 0; fi
	-if test -x /usr/bin/ranlib; then /usr/bin/ranlib libhmmer.a; else exit 0; fi
	-if test -x $(RANLIB); then $(RANLIB) libhmmer.a; else exit 0; fi
	chmod 644 libhmmer.a


#######
## Targets defining different build styles;
#######

# Make a "production" copy
prodv:
	make clean
	make CFLAGS="-O2" MDEFS="" LIBS="-lm"

# Make a debugging copy
dbv:
	make clean	
	make CFLAGS="-g -Wall -pedantic" MDEFS="" LIBS="-lm"

# Make a memory debugging copy
memv: 
	make clean	
	make CFLAGS="-g -Wall" MDEFS="-DMEMDEBUG -L." LIBS="-lm -ldbmalloc"


# Make an Electric Fence version
efencev:
	make clean
	make CFLAGS="-g -Wall" MDEFS="" LIBS="-lm -lefence"

# Make a gprof version
gprofv:
	make clean
	make CFLAGS="-pg" MDEFS="" LIBS="-lm" 


#######
## Targets only used in development, irrelevant to distribution copy
#######

# Extract squid files from RCS.
squid:  
	for file in $(SQUIDHDRS); do\
	   co $(HOME)/src/squid/RCS/$$file,v $$file;\
	done
	for file in $(SQUIDSRC); do\
	   co $(HOME)/src/squid/RCS/$$file,v $$file;\
	done

#######
## Miscellaneous
#######

clean:
	-rm -f *.o *~ Makefile.bak core $(PROGS) TAGS gmon.out libhmmer.a

tags:
	etags -t $(SRC) $(HDRS) $(READMES) 

lint:
	lint $(MDEFS) $(SRC) 


